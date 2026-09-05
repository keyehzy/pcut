#pragma once
#include <pcut/coefficients.hpp>
#include <pcut/polynomial.hpp>

namespace pcut {
// Immutable shared-prefix program for sparse applications of the universal series.
// Safe to share between threads; all evaluation scratch storage is local.
class EffectiveOperator {
public:
    explicit EffectiveOperator(const Coefficients& coefficients);
    [[nodiscard]] unsigned order() const noexcept { return order_; }
    [[nodiscard]] std::vector<SparseState> apply(const ClusterModel& model, State input) const;
    [[nodiscard]] std::vector<Matrix> block(const ClusterModel& model, const std::vector<State>& basis) const;
    [[nodiscard]] Series vacuum(const ClusterModel& model) const;
    [[nodiscard]] SymbolicBlock symbolic_block(const ClusterModel& model, const std::vector<State>& basis,
        const std::vector<std::size_t>& channel_edges = {}) const;
    [[nodiscard]] const std::string& identity() const noexcept { return identity_; }
private:
    std::string identity_;
    struct Node {
        std::map<int, std::size_t> children;
        Rational coefficient = 0;
        unsigned depth = 0;
    };
    unsigned order_;
    std::vector<int> changes_;
    std::vector<Node> nodes_;
    struct FactoredEdge { int change; std::size_t target; Rational scale; };
    struct FactoredNode { Rational coefficient; std::vector<FactoredEdge> edges; };
    std::vector<std::vector<FactoredNode>> layers_;
    Rational root_scale_=1;
    void factor_program();
};
struct ParticleState {
    std::size_t site;
    unsigned local;
    State state;
};
[[nodiscard]] std::vector<ParticleState> one_particle_basis(const ClusterModel& model);
}
