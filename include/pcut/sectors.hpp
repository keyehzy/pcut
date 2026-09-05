#pragma once
#include <pcut/linked.hpp>

namespace pcut {
struct SectorOptions {
    std::size_t max_basis = 1024;
    std::size_t max_contractions = 10'000'000;
    std::size_t max_kernels = 1'000'000;
    SolverOptions solver;
};
struct IrreducibleSectors {
    std::vector<State> basis; // all product states with Q <= max_charge, ordered by Q then encoding
    std::vector<Matrix> kernels;
};
// Hubbard-operator kernels: recursively subtract every unchanged occupied
// spectator subset. For local charges {0,1,...,1}, these are the n-particle
// irreducible H_n. More general local charges give charge-conserving conversion
// kernels. Off-site tensor Hubbard operators commute; this is not a graded
// fermionic normal-ordering convention.
[[nodiscard]] IrreducibleSectors irreducible_sectors(const ClusterModel& model,
                                                   const EffectiveOperator& effective,
                                                   unsigned max_charge,
                                                   SectorOptions options = {});
struct Excitation {
    Site site;
    unsigned local;
    auto operator<=>(const Excitation&) const = default;
};
struct Kernel {
    std::vector<Excitation> output;
    std::vector<Excitation> input;
    auto operator<=>(const Kernel&) const = default;
};
struct LinkedSectors {
    Series energy_per_cell;
    // Translation-normalized coefficients of products of local creation and
    // annihilation Hubbard operators; no factorials or permutation factors.
    // Sum each key over all unit-cell translations to reconstruct H_eff.
    std::map<Kernel, Series> kernels;
};
[[nodiscard]] LinkedSectors linked_expand_sectors(const ClusterCatalog& catalog,
                                                 const EffectiveOperator& effective,
                                                 unsigned max_charge,
                                                 SectorOptions options = {});
}
