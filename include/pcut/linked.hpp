#pragma once
#include <pcut/effective.hpp>
#include <pcut/lattice.hpp>
#include <functional>
#include <map>

namespace pcut {
using ClusterEvaluator = std::function<Series(const ClusterModel&)>;
struct ScalarExpansion {
    Series per_cell;
    std::vector<Series> weights;
};
// For any user-supplied cluster-additive scalar. Evaluator returns perturbative
// corrections (order zero must be zero); caller supplies the on-site reference.
[[nodiscard]] ScalarExpansion linked_scalar(const ClusterCatalog& catalog, unsigned order,
                                             const ClusterEvaluator& evaluate_cluster,
                                             Complex reference_per_cell = 0);
struct Flavor {
    std::size_t basis;
    unsigned local;
    auto operator<=>(const Flavor&) const = default;
};
struct Hopping {
    std::size_t output;
    std::size_t input;
    Coordinate displacement; // output cell - input cell
    auto operator<=>(const Hopping&) const = default;
};
struct LinkedOptions {
    bool one_particle = true;
};
struct LinkedResult {
    Series energy_per_cell;
    std::vector<Flavor> flavors;
    std::map<Hopping, Series> hopping;
    std::vector<Series> vacuum_weights;
    // k is in reciprocal unit-cell coordinates (radians). Fourier convention:
    // |k,a> = sum_R exp(i k.R)|R,a>, so H_ab(k)=sum_d t_ab(d)exp(-i k.d).
    [[nodiscard]] std::vector<Matrix> bloch_series(const std::vector<double>& k) const;
    [[nodiscard]] Matrix bloch(const std::vector<double>& k, double lambda) const;
    unsigned dimension = 0;
};
[[nodiscard]] std::vector<int> charge_changes(const PeriodicLattice& lattice);
[[nodiscard]] LinkedResult linked_expand(const ClusterCatalog& catalog,
                                         const EffectiveOperator& effective,
                                         LinkedOptions options = {});
}
