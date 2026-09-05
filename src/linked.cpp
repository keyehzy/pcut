#include <pcut/linked.hpp>
#include "detail.hpp"
#include <algorithm>
#include <cmath>
#include <limits>
#include <set>
#include <stdexcept>

namespace pcut {
namespace {
void check_order(const WhiteGraphExpansion& catalog, unsigned order) {
    if (catalog.max_edges() < order) throw std::invalid_argument("catalog must include clusters through perturbation order");
}
}
ScalarExpansion linked_scalar(const WhiteGraphExpansion& catalog, unsigned order,
                               const ScalarEvaluator& evaluator, Complex reference_per_cell) {
    check_order(catalog, order);
    if (!std::isfinite(reference_per_cell.real()) || !std::isfinite(reference_per_cell.imag()))
        throw std::invalid_argument("nonfinite reference energy");
    ScalarExpansion result{Series(order+1), {}};
    result.per_cell[0]=reference_per_cell;
    std::vector<SymbolicSeries> weights;
    for (std::size_t graph=0;graph<catalog.graphs().size();++graph) {
        const auto& entry=catalog.graphs()[graph];
        if (entry.graph.edges.size()>order) break;
        auto weight=*catalog.cache()->scalar(entry,catalog.structure().gap,order,evaluator);
        for (const auto& sub : catalog.graph_subclusters(graph)) for (unsigned n=1;n<=order;++n)
            add_polynomial(weight[n],weights[sub.graph][n],-1,&sub.map.channels);
        weights.push_back(std::move(weight));
    }
    for (std::size_t index=0;index<catalog.embeddings().size();++index) {
        const auto& embedding=catalog.embeddings()[index];
        if (embedding.edges.size()>order) break;
        auto value=substitute(weights[embedding.graph],catalog.couplings(index));
        for (unsigned n=1;n<=order;++n) result.per_cell[n]+=value[n];
        result.weights.push_back(std::move(value));
    }
    return result;
}
std::vector<int> charge_changes(const PeriodicLattice& lattice) {
    lattice.validate();
    std::set<int> changes;
    for (const auto& interaction : lattice.interactions) {
        std::vector<LocalSpace> spaces;
        std::vector<std::size_t> legs;
        for (const auto& site : interaction.legs) {
            legs.push_back(spaces.size());
            spaces.push_back(lattice.cell[site.basis]);
        }
        std::vector<LocalTerm> terms;
        for (const auto& channel : interaction.channels)
            terms.push_back({legs,channel.op.matrix,channel.op.fermionic});
        // Compile unbound channels separately, preserving letters even when
        // numerical couplings vanish or cancel.
        const ClusterModel model(std::move(spaces),std::move(terms),lattice.gap);
        changes.insert(model.changes().begin(),model.changes().end());
    }
    return {changes.begin(),changes.end()};
}
LinkedResult linked_expand(const WhiteGraphExpansion& catalog, const EffectiveOperator& effective, LinkedOptions options) {
    const auto order = effective.order();
    check_order(catalog,order);
    const auto& lattice = catalog.structure();
    for (const auto& space : lattice.cell) space.require_product_vacuum();
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
    for (std::size_t index=0;index<catalog.embeddings().size();++index) {
        const auto& entry=catalog.embeddings()[index];
        if (entry.edges.size()>order) break;
        const auto& model = catalog.structural_model(index);
        std::vector<ParticleState> basis;
        std::vector<State> states{0};
        if (options.one_particle) {
            basis=one_particle_basis(model);
            for (const auto& p : basis) states.push_back(p.state);
        }
        const auto raw=catalog.block(index,effective,states,true);
        Series energy(order+1);
        for (unsigned n=1;n<=order;++n) energy[n]=raw[n](0,0);
        for (unsigned n = 1; n <= order; ++n) result.energy_per_cell[n] += energy[n];
        result.vacuum_weights.push_back(std::move(energy));
        if (options.one_particle) {
            for (std::size_t i = 0; i < basis.size(); ++i) for (std::size_t j = 0; j < basis.size(); ++j) {
                const auto& out = basis[i]; const auto& in = basis[j];
                const auto& out_site = entry.sites[out.site]; const auto& in_site = entry.sites[in.site];
                const auto out_flavor = static_cast<std::size_t>(std::lower_bound(result.flavors.begin(),result.flavors.end(),Flavor{out_site.basis,out.local})-result.flavors.begin());
                const auto in_flavor = static_cast<std::size_t>(std::lower_bound(result.flavors.begin(),result.flavors.end(),Flavor{in_site.basis,in.local})-result.flavors.begin());
                auto displacement=detail::translate(out_site.cell,in_site.cell,-1);
                Hopping key{out_flavor,in_flavor,std::move(displacement)};
                auto& series = result.hopping[key];
                if (series.empty()) series.resize(order+1);
                for (unsigned n = 1; n <= order; ++n)
                    series[n] += raw[n](static_cast<Eigen::Index>(i+1),static_cast<Eigen::Index>(j+1))
                        - (i==j ? raw[n](0,0) : Complex{});
            }
        }
    }
    return result;
}
std::vector<Matrix> LinkedResult::bloch_series(const std::vector<double>& k) const {
    if (k.size() != dimension || !std::all_of(k.begin(),k.end(),[](double x){return std::isfinite(x);}))
        throw std::invalid_argument("invalid momentum dimension or value");
    auto result=detail::matrix_series(flavors.size(),energy_per_cell.size());
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
