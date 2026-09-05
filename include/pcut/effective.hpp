#pragma once
#include <pcut/coefficients.hpp>
#include <pcut/model.hpp>

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
private:
    struct Node {
        std::map<int, std::size_t> children;
        double coefficient = 0;
        unsigned depth = 0;
    };
    unsigned order_;
    std::vector<int> changes_;
    std::vector<Node> nodes_;
};
struct ParticleState {
    std::size_t site;
    unsigned local;
    State state;
};
[[nodiscard]] std::vector<ParticleState> one_particle_basis(const ClusterModel& model);
}
