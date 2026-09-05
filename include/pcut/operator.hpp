#pragma once
#include <pcut/linked.hpp>
#include <optional>

namespace pcut {
struct OperatorOptions {
    std::size_t max_basis = 4096;
    std::size_t max_basis_visits = 1'000'000;
    // Counts complex matrix entries across all orders (and all stored weights).
    std::size_t max_matrix_elements = 32'000'000;
    // Bounds parent basis * child basis * number of orders per embedding.
    std::size_t max_embedding_work = 100'000'000;
    SolverOptions solver;
};
struct OperatorBlock {
    std::vector<LocalSpace> spaces;
    std::vector<State> basis; // original mixed-radix encodings, ascending
    std::vector<Matrix> coefficients;
    [[nodiscard]] Matrix evaluate(double lambda) const;
};
// External Q=0 states only; intermediate sectors are unrestricted. Optional
// particle selection is for finite evaluations, NEVER for linked weights.
[[nodiscard]] std::vector<State> zero_charge_basis(const ClusterModel& model,
    std::optional<int> particles = {}, std::size_t max_basis = 4096,
    std::size_t max_visits = 1'000'000);
[[nodiscard]] OperatorBlock zero_charge_operator(const ClusterModel& model,
    const EffectiveOperator& effective, std::optional<int> particles = {}, OperatorOptions options = {});
// Add scale * (child operator tensor identity spectators), with graded site
// permutation signs whenever parity metadata is present. Supports selected
// parent sectors; the child's basis must contain the entire Q=0 manifold.
void add_embedded_operator(OperatorBlock& parent, const OperatorBlock& child,
    const std::vector<std::size_t>& vertex_map, Complex scale = 1,
    std::size_t max_work = 100'000'000);
struct OperatorWeight {
    Cluster edges;
    std::vector<Site> sites;
    OperatorBlock block;
};
struct LinkedOperator {
    PeriodicLattice lattice;
    unsigned order = 0;
    // Each connected colored edge animal is summed over ALL translations once.
    // These are infinite-lattice operator couplings, with identity spectators.
    std::vector<OperatorWeight> weights;
};
// Dedicated operator-valued linked expansion in the entire tensor Q=0 manifold.
// Reference constants are assembled separately, not subtracted as vacuum energy.
[[nodiscard]] LinkedOperator linked_zero_charge(const ClusterCatalog& catalog,
    const EffectiveOperator& effective, OperatorOptions options = {});
// Restrict the translation sum to a supplied finite set of colored edges.
// sites specifies the Fock ordering and may contain isolated spectator sites.
// This is an OPEN finite-system evaluation, not an infinite-lattice spectrum.
[[nodiscard]] OperatorBlock assemble_operator(const LinkedOperator& linked,
    const std::vector<Site>& sites, const Cluster& edges,
    std::optional<int> particles = {}, OperatorOptions options = {});
}
