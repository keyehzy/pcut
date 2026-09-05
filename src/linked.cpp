#include <pcut/linked.hpp>
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace pcut {
namespace {
void check_order(const ClusterCatalog& catalog, unsigned order) {
    if (catalog.max_edges() < order) throw std::invalid_argument("catalog must include clusters through perturbation order");
}
}
ScalarExpansion linked_scalar(const ClusterCatalog& catalog, unsigned order,
                               const ClusterEvaluator& evaluate_cluster, Complex reference_per_cell) {
    check_order(catalog, order);
    if (!std::isfinite(reference_per_cell.real()) || !std::isfinite(reference_per_cell.imag()))
        throw std::invalid_argument("nonfinite reference energy");
    ScalarExpansion result{Series(order+1), {}};
    result.per_cell[0] = reference_per_cell;
    for (const auto& entry : catalog.entries()) {
        Series weight = evaluate_cluster(cluster_model(catalog.lattice(), entry.edges));
        if (weight.size() != order+1 || weight[0] != Complex{})
            throw std::invalid_argument("scalar callback must return order+1 corrections with zero constant");
        for (const auto& c : weight) if (!std::isfinite(c.real()) || !std::isfinite(c.imag()))
            throw std::invalid_argument("nonfinite scalar cluster value");
        for (const auto& sub : entry.subclusters)
            for (unsigned n = 1; n <= order; ++n) weight[n] -= result.weights[sub.index][n];
        for (unsigned n = 1; n <= order; ++n) result.per_cell[n] += weight[n];
        result.weights.push_back(std::move(weight));
    }
    return result;
}
std::vector<int> charge_changes(const PeriodicLattice& lattice) {
    lattice.validate();
    std::set<int> changes;
    for (std::size_t type = 0; type < lattice.interactions.size(); ++type) {
        const auto model = cluster_model(lattice, {{type, Coordinate(lattice.dimension,0)}});
        changes.insert(model.changes().begin(),model.changes().end());
    }
    return {changes.begin(),changes.end()};
}
LinkedResult linked_expand(const ClusterCatalog& catalog, const EffectiveOperator& effective, LinkedOptions options) {
    const auto order = effective.order();
    check_order(catalog,order);
    const auto& lattice = catalog.lattice();
    LinkedResult result;
    result.dimension = lattice.dimension;
    result.energy_per_cell.resize(order+1);
    for (std::size_t b = 0; b < lattice.cell.size(); ++b) {
        result.energy_per_cell[0] += lattice.cell[b].vacuum_energy;
        if (options.one_particle) for (unsigned local = 1; local < lattice.cell[b].charges.size(); ++local)
            if (lattice.cell[b].charges[local] == 1) result.flavors.push_back({b,local});
    }
    for (std::size_t a = 0; a < result.flavors.size(); ++a) {
        Series s(order+1); s[0] = lattice.gap;
        result.hopping.emplace(Hopping{a,a,Coordinate(lattice.dimension,0)},std::move(s));
    }
    struct Weight { std::vector<ParticleState> basis; std::vector<Matrix> h1; };
    std::vector<Weight> weights;
    for (const auto& entry : catalog.entries()) {
        const auto model = cluster_model(lattice, entry.edges);
        const auto raw_energy = effective.vacuum(model, options.solver);
        auto energy = raw_energy;
        energy[0] = 0;
        Weight weight;
        if (options.one_particle) {
            weight.basis = one_particle_basis(model);
            std::vector<State> states;
            for (const auto& p : weight.basis) states.push_back(p.state);
            weight.h1 = effective.block(model, states, options.solver);
            weight.h1[0].setZero(); // bare gap already embedded once per local flavor
            for (unsigned n = 1; n <= order; ++n) weight.h1[n].diagonal().array() -= raw_energy[n];
        }
        for (const auto& sub : entry.subclusters) {
            for (unsigned n = 1; n <= order; ++n) energy[n] -= result.vacuum_weights[sub.index][n];
            if (options.one_particle) {
                const auto& child = weights[sub.index];
                std::vector<Eigen::Index> map;
                for (const auto& p : child.basis) {
                    const auto parent_site = sub.vertex_map[p.site];
                    const auto it = std::find_if(weight.basis.begin(),weight.basis.end(),[&](const ParticleState& s) {
                        return s.site == parent_site && s.local == p.local;
                    });
                    if (it == weight.basis.end()) throw std::logic_error("one-particle embedding mismatch");
                    map.push_back(static_cast<Eigen::Index>(it-weight.basis.begin()));
                }
                for (unsigned n = 1; n <= order; ++n)
                    for (std::size_t i = 0; i < map.size(); ++i) for (std::size_t j = 0; j < map.size(); ++j)
                        weight.h1[n](map[i],map[j]) -= child.h1[n](static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j));
            }
        }
        for (unsigned n = 1; n <= order; ++n) result.energy_per_cell[n] += energy[n];
        result.vacuum_weights.push_back(std::move(energy));
        if (options.one_particle) {
            for (std::size_t i = 0; i < weight.basis.size(); ++i) for (std::size_t j = 0; j < weight.basis.size(); ++j) {
                const auto& out = weight.basis[i]; const auto& in = weight.basis[j];
                const auto& out_site = entry.sites[out.site]; const auto& in_site = entry.sites[in.site];
                const auto out_flavor = static_cast<std::size_t>(std::lower_bound(result.flavors.begin(),result.flavors.end(),Flavor{out_site.basis,out.local})-result.flavors.begin());
                const auto in_flavor = static_cast<std::size_t>(std::lower_bound(result.flavors.begin(),result.flavors.end(),Flavor{in_site.basis,in.local})-result.flavors.begin());
                Coordinate displacement(lattice.dimension);
                for (unsigned d = 0; d < lattice.dimension; ++d) {
                    const auto delta = static_cast<long long>(out_site.cell[d])-in_site.cell[d];
                    if (delta < std::numeric_limits<int>::min() || delta > std::numeric_limits<int>::max())
                        throw std::overflow_error("hopping displacement overflow");
                    displacement[d] = static_cast<int>(delta);
                }
                auto& series = result.hopping[{out_flavor,in_flavor,std::move(displacement)}];
                if (series.empty()) series.resize(order+1);
                for (unsigned n = 1; n <= order; ++n)
                    series[n] += weight.h1[n](static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j));
            }
        }
        weights.push_back(std::move(weight));
    }
    return result;
}
std::vector<Matrix> LinkedResult::bloch_series(const std::vector<double>& k) const {
    if (k.size() != dimension || !std::all_of(k.begin(),k.end(),[](double x){return std::isfinite(x);}))
        throw std::invalid_argument("invalid momentum dimension or value");
    const auto size = static_cast<Eigen::Index>(flavors.size());
    std::vector<Matrix> result(energy_per_cell.size(), Matrix::Zero(size,size));
    for (const auto& [h,s] : hopping) {
        double phase = 0;
        for (unsigned d = 0; d < dimension; ++d) phase -= k[d]*h.displacement[d];
        const auto factor = std::polar(1.0,phase);
        for (std::size_t n = 0; n < s.size(); ++n)
            result[n](static_cast<Eigen::Index>(h.output),static_cast<Eigen::Index>(h.input)) += s[n]*factor;
    }
    return result;
}
Matrix LinkedResult::bloch(const std::vector<double>& k, double lambda) const {
    if (!std::isfinite(lambda)) throw std::invalid_argument("lambda must be finite");
    const auto coefficients = bloch_series(k);
    Matrix result = Matrix::Zero(static_cast<Eigen::Index>(flavors.size()),static_cast<Eigen::Index>(flavors.size()));
    for (auto it = coefficients.rbegin(); it != coefficients.rend(); ++it) result = lambda*result + *it;
    return result;
}
}
