#include <pcut/model.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace pcut {
void LocalSpace::validate() const {
    if (charges.empty() || charges[0] != 0 || !std::isfinite(vacuum_energy))
        throw std::invalid_argument("local space needs basis-zero vacuum and finite reference energy");
    for (std::size_t i = 1; i < charges.size(); ++i)
        if (charges[i] < 0 || charges[i] > 1'000'000)
            throw std::invalid_argument("local charges must be nonnegative integers <= 1000000");
    if ((!particles.empty() && particles.size() != charges.size()) ||
        (!parity.empty() && parity.size() != charges.size()))
        throw std::invalid_argument("local metadata dimension mismatch");
    for (std::size_t i = 0; i < charges.size(); ++i) {
        if (!particles.empty() && (particles[i] < 0 || particles[i] > 1'000'000))
            throw std::invalid_argument("invalid local particle number");
        if (!parity.empty() && (parity[i] > 1 ||
            (!particles.empty() && parity[i] != static_cast<unsigned>(particles[i] % 2))))
            throw std::invalid_argument("invalid or inconsistent local parity");
    }
}
void LocalSpace::require_product_vacuum() const {
    validate();
    if (std::count(charges.begin(), charges.end(), 0) != 1)
        throw std::invalid_argument("driver requires a unique local charge-zero product vacuum");
}
namespace {
std::vector<std::pair<std::size_t,std::size_t>> gather_inversions(
    std::size_t size, const std::vector<std::size_t>& support) {
    std::vector<std::size_t> order = support;
    std::set<std::size_t> used;
    for (auto s : support) if (s >= size || !used.insert(s).second)
        throw std::invalid_argument("invalid graded embedding support");
    for (std::size_t s = 0; s < size; ++s) if (!used.contains(s)) order.push_back(s);
    std::vector<std::pair<std::size_t,std::size_t>> inversions;
    for (std::size_t i=0; i<size; ++i) for (std::size_t j=i+1; j<size; ++j)
        if (order[i] > order[j]) inversions.emplace_back(order[i],order[j]);
    return inversions;
}
}
int graded_gather_sign(const std::vector<LocalSpace>& spaces, const std::vector<unsigned>& local,
                       const std::vector<std::size_t>& support) {
    if (local.size() != spaces.size()) throw std::invalid_argument("graded state size mismatch");
    for (std::size_t s=0; s<spaces.size(); ++s) {
        spaces[s].validate();
        if (local[s] >= spaces[s].charges.size() || spaces[s].parity.empty())
            throw std::invalid_argument("graded embedding requires parity metadata");
    }
    unsigned parity=0;
    for (auto [a,b] : gather_inversions(spaces.size(),support))
        parity ^= spaces[a].parity[local[a]] & spaces[b].parity[local[b]];
    return parity ? -1 : 1;
}
ClusterModel::ClusterModel(std::vector<LocalSpace> spaces, std::vector<LocalTerm> terms,
                           double gap, double zero_tolerance)
    : spaces_(std::move(spaces)), gap_(gap) {
    if (!std::isfinite(gap) || gap <= 0 || !std::isfinite(zero_tolerance) || zero_tolerance < 0)
        throw std::invalid_argument("gap must be positive and tolerance nonnegative, both finite");
    if (spaces_.size() > 1000) throw std::length_error("too many local spaces");
    for (const auto& space : spaces_) {
        space.validate();
        if (space.particles.empty()) conserves_particles_ = false;
        strides_.push_back(dimension_);
        if (dimension_ > std::numeric_limits<State>::max() / space.charges.size())
            throw std::length_error("tensor basis exceeds 64-bit state encoding");
        dimension_ *= space.charges.size();
    }
    const bool graded=std::any_of(spaces_.begin(),spaces_.end(),[](const auto& s) { return !s.parity.empty(); });
    if (graded) for (const auto& space : spaces_) if (space.parity.empty()) operator_grading_valid_=false;
    std::set<int> changes;
    for (const auto& term : terms) {
        if (term.sites.empty()) throw std::invalid_argument("local term must have support");
        std::set<std::size_t> unique;
        CompiledTerm compiled;
        compiled.sites = term.sites;
        compiled.fermionic = term.fermionic;
        fermionic_ |= term.fermionic;
        if (term.fermionic) {
            for (const auto& space : spaces_) if (space.parity.empty())
                throw std::invalid_argument("fermionic terms require parity on all sites");
            compiled.inversions = gather_inversions(sites(), term.sites);
        }
        State dim = 1;
        for (auto s : term.sites) {
            if (s >= sites() || !unique.insert(s).second) throw std::invalid_argument("invalid or repeated operator leg");
            compiled.local_stride.push_back(dim);
            dim *= spaces_[s].charges.size();
        }
        if (dim > static_cast<State>(std::numeric_limits<Eigen::Index>::max()) ||
            term.matrix.rows() != static_cast<Eigen::Index>(dim) || term.matrix.cols() != static_cast<Eigen::Index>(dim))
            throw std::invalid_argument("local matrix dimension does not match tensor legs");
        if (!term.matrix.allFinite() || !term.matrix.isApprox(term.matrix.adjoint(), 1e-12))
            throw std::invalid_argument("local Hamiltonian terms must be finite and Hermitian");
        std::vector<int> q(static_cast<std::size_t>(dim), 0);
        std::vector<unsigned> parity(static_cast<std::size_t>(dim), 0);
        std::vector<int> number(static_cast<std::size_t>(dim), 0);
        for (State i = 0; i < dim; ++i)
            for (std::size_t leg = 0; leg < term.sites.size(); ++leg) {
                const auto& space = spaces_[term.sites[leg]];
                const auto local = (i / compiled.local_stride[leg]) % space.charges.size();
                q[i] += space.charges[local];
                if (!space.particles.empty()) number[i] += space.particles[local];
                if (term.fermionic) parity[i] ^= space.parity[local];
            }
        for (State col = 0; col < dim; ++col) for (State row = 0; row < dim; ++row) {
            const Complex v = term.matrix(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
            if (std::abs(v) <= zero_tolerance) continue;
            if (term.fermionic && parity[row] != parity[col])
                throw std::invalid_argument("fermionic Hamiltonian terms must preserve total parity");
            if (number[row] != number[col]) conserves_particles_ = false;
            if (graded && !term.fermionic) for (std::size_t leg=0;leg<term.sites.size();++leg) {
                const auto& space=spaces_[term.sites[leg]];
                if (!space.parity.empty() &&
                    space.parity[(row/compiled.local_stride[leg])%space.charges.size()] !=
                    space.parity[(col/compiled.local_stride[leg])%space.charges.size()])
                    operator_grading_valid_=false;
            }
            const int change = q[row] - q[col];
            auto& columns = compiled.by_change[change];
            if (columns.empty()) columns.resize(static_cast<std::size_t>(dim));
            columns[col].push_back({row, v});
            changes.insert(change);
        }
        terms_.push_back(std::move(compiled));
    }
    changes_.assign(changes.begin(), changes.end());
}
double ClusterModel::vacuum_energy() const noexcept {
    double e = 0;
    for (const auto& s : spaces_) e += s.vacuum_energy;
    return e;
}
State ClusterModel::encode(const std::vector<unsigned>& local) const {
    if (local.size() != sites()) throw std::invalid_argument("state has wrong site count");
    State state = 0;
    for (std::size_t i = 0; i < sites(); ++i) {
        if (local[i] >= spaces_[i].charges.size()) throw std::out_of_range("local basis state");
        state += local[i] * strides_[i];
    }
    return state;
}
std::vector<unsigned> ClusterModel::decode(State state) const {
    if (state >= dimension_) throw std::out_of_range("tensor basis state");
    std::vector<unsigned> result(sites());
    for (std::size_t i = 0; i < sites(); ++i)
        result[i] = static_cast<unsigned>((state / strides_[i]) % spaces_[i].charges.size());
    return result;
}
int ClusterModel::charge(State state) const {
    const auto local = decode(state);
    int q = 0;
    for (std::size_t i = 0; i < sites(); ++i) q += spaces_[i].charges[local[i]];
    return q;
}
int ClusterModel::particle_number(State state) const {
    const auto local = decode(state);
    int n = 0;
    for (std::size_t s=0; s<sites(); ++s) {
        if (spaces_[s].particles.empty()) throw std::invalid_argument("particle metadata missing");
        n += spaces_[s].particles[local[s]];
    }
    return n;
}
void ClusterModel::require_product_vacuum() const {
    for (const auto& space : spaces_) space.require_product_vacuum();
}
void ClusterModel::require_operator_grading() const {
    if (!operator_grading_valid_)
        throw std::invalid_argument("operator linking requires complete parity metadata and graded parity-changing terms");
}
SparseState ClusterModel::apply(int change, const SparseState& input, std::size_t max_states) const {
    SparseState result;
    for (const auto& [state, amplitude] : input) {
        if (state >= dimension_ || !std::isfinite(amplitude.real()) || !std::isfinite(amplitude.imag()))
            throw std::invalid_argument("invalid sparse input state");
        if (amplitude == Complex{}) continue;
        for (const auto& term : terms_) {
            const auto block = term.by_change.find(change);
            if (block == term.by_change.end()) continue;
            State column = 0, removed = 0;
            for (std::size_t leg = 0; leg < term.sites.size(); ++leg) {
                const auto s = term.sites[leg];
                const State value = (state / strides_[s]) % spaces_[s].charges.size();
                column += value * term.local_stride[leg];
                removed += value * strides_[s];
            }
            for (const auto& transition : block->second[column]) {
                State output = state - removed;
                for (std::size_t leg = 0; leg < term.sites.size(); ++leg) {
                    const auto s = term.sites[leg];
                    output += ((transition.output / term.local_stride[leg]) % spaces_[s].charges.size()) * strides_[s];
                }
                unsigned parity = 0;
                if (term.fermionic) for (auto [a,b] : term.inversions) {
                    const auto pa = [&](State x, std::size_t site) {
                        return spaces_[site].parity[(x / strides_[site]) % spaces_[site].charges.size()];
                    };
                    parity ^= (pa(state,a) & pa(state,b)) ^ (pa(output,a) & pa(output,b));
                }
                result[output] += (parity ? -amplitude : amplitude) * transition.value;
                if (result.size() > max_states) throw std::length_error("sparse intermediate-state budget exceeded");
            }
        }
    }
    // Remove exact cancellations only; no numerical Q-sector or amplitude truncation.
    for (auto it = result.begin(); it != result.end();)
        if (it->second == Complex{}) it = result.erase(it); else ++it;
    return result;
}
Matrix ClusterModel::dense_hamiltonian(double lambda, State max_dimension) const {
    if (!std::isfinite(lambda)) throw std::invalid_argument("lambda must be finite");
    if (dimension_ > max_dimension || dimension_ > static_cast<State>(std::numeric_limits<Eigen::Index>::max()))
        throw std::length_error("dense validation matrix exceeds dimension limit");
    Matrix h = Matrix::Zero(static_cast<Eigen::Index>(dimension_), static_cast<Eigen::Index>(dimension_));
    for (State i = 0; i < dimension_; ++i) {
        h(static_cast<Eigen::Index>(i), static_cast<Eigen::Index>(i)) = vacuum_energy() + gap_ * charge(i);
        for (int m : changes_) for (const auto& [j, v] : apply(m, {{i,1.0}}))
            h(static_cast<Eigen::Index>(j), static_cast<Eigen::Index>(i)) += lambda * v;
    }
    return h;
}
Complex evaluate(const Series& series, double lambda) {
    if (!std::isfinite(lambda)) throw std::invalid_argument("lambda must be finite");
    Complex value = 0;
    for (auto it = series.rbegin(); it != series.rend(); ++it) value = value * lambda + *it;
    return value;
}
}
