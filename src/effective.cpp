#include <pcut/effective.hpp>
#include "detail.hpp"
#include <algorithm>
#include <cmath>
#include <functional>
#include <set>
#include <stdexcept>

namespace pcut {
EffectiveOperator::EffectiveOperator(const Coefficients& coefficients)
    : order_(coefficients.order()), changes_(coefficients.changes()), nodes_(1) {
    for (const auto& [word, c] : coefficients.terms()) {
        std::size_t node = 0;
        unsigned depth = 0;
        for (auto it = word.rbegin(); it != word.rend(); ++it) {
            ++depth;
            const auto existing = nodes_[node].children.find(*it);
            if (existing != nodes_[node].children.end()) node = existing->second;
            else {
                const auto next = nodes_.size();
                nodes_[node].children.emplace(*it, next);
                nodes_.push_back({{}, 0.0, depth});
                node = next;
            }
        }
        nodes_[node].coefficient = c.convert_to<double>();
    }
}
std::vector<SparseState> EffectiveOperator::apply(const ClusterModel& model, State input,
                                                 SolverOptions options) const {
    if (input >= model.dimension()) throw std::out_of_range("effective input state");
    if (options.max_states == 0) throw std::invalid_argument("state budget must be positive");
    for (int m : model.changes()) if (!std::binary_search(changes_.begin(), changes_.end(), m))
        throw std::invalid_argument("coefficient alphabet does not cover model charge changes");
    std::vector<SparseState> output(order_ + 1);
    output[0][input] = model.vacuum_energy() + model.gap() * model.charge(input);
    std::function<void(std::size_t,const SparseState&)> visit = [&](std::size_t index, const SparseState& state) {
        const auto& node = nodes_[index];
        if (node.coefficient != 0) {
            for (const auto& [s, v] : state) output[node.depth][s] += node.coefficient * v;
            if (output[node.depth].size() > options.max_states) throw std::length_error("output state budget exceeded");
        }
        for (const auto& [change, next] : node.children) {
            auto transformed = model.apply_scaled(change, state, options.max_states, model.gap());
            if (!transformed.empty()) visit(next, transformed);
        }
    };
    visit(0, {{input, 1.0}});
    // Traverse with V/Delta, then restore units after summing each order.
    for (unsigned n=1;n<=order_;++n) for (auto& [s,v] : output[n]) {
        v *= model.gap();
        if (!std::isfinite(v.real()) || !std::isfinite(v.imag()))
            throw std::overflow_error("pCUT energy scale overflow");
    }
    return output;
}
std::vector<Matrix> EffectiveOperator::block(const ClusterModel& model, const std::vector<State>& basis,
                                            SolverOptions options) const {
    const std::set<State> unique(basis.begin(), basis.end());
    if (unique.size() != basis.size()) throw std::invalid_argument("duplicate block basis state");
    for (auto state : basis) if (state>=model.dimension()) throw std::out_of_range("block basis state");
    auto result=detail::matrix_series(basis.size(),static_cast<std::size_t>(order_)+1,options.max_matrix_elements);
    for (std::size_t col = 0; col < basis.size(); ++col) {
        const auto action = apply(model, basis[col], options);
        for (unsigned n = 0; n <= order_; ++n)
            for (std::size_t row = 0; row < basis.size(); ++row) {
                const auto it = action[n].find(basis[row]);
                if (it != action[n].end()) result[n](static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col)) = it->second;
            }
    }
    return result;
}
Series EffectiveOperator::vacuum(const ClusterModel& model, SolverOptions options) const {
    model.require_product_vacuum();
    const auto action = apply(model, 0, options);
    Series result(order_+1);
    for (unsigned n = 0; n <= order_; ++n) {
        const auto it = action[n].find(0);
        if (it != action[n].end()) result[n] = it->second;
    }
    return result;
}
std::vector<ParticleState> one_particle_basis(const ClusterModel& model) {
    model.require_product_vacuum();
    std::vector<ParticleState> basis;
    std::vector<unsigned> local(model.sites(), 0);
    for (std::size_t site = 0; site < model.sites(); ++site) {
        const auto& space = model.spaces()[site];
        for (unsigned l = 1; l < space.charges.size(); ++l) if (space.charges[l] == 1) {
            local[site] = l;
            basis.push_back({site, l, model.encode(local)});
        }
        local[site] = 0;
    }
    return basis;
}
}
