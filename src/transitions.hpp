#pragma once
#include <pcut/model.hpp>

namespace pcut::detail {
// Enumerate each channel separately in physical tensor order. Both numerical
// and symbolic accumulation use these output encodings and Fock signs.
struct ModelTransitions {
    template<class Emit>
    static void enumerate(const ClusterModel& model, int change, State state, Emit&& emit) {
        for (std::size_t channel=0;channel<model.terms_.size();++channel) {
            const auto& term=model.terms_[channel];
            const auto block = term.by_change->find(change);
            if (block == term.by_change->end()) continue;
            State column = 0, removed = 0;
            for (std::size_t leg = 0; leg < term.sites.size(); ++leg) {
                const auto s = term.sites[leg];
                const State value = (state / model.strides_[s]) % model.spaces_[s].charges.size();
                column += value * term.local_stride[leg];
                removed += value * model.strides_[s];
            }
            for (const auto& transition : block->second[column]) {
                State output = state - removed;
                for (std::size_t leg = 0; leg < term.sites.size(); ++leg) {
                    const auto s = term.sites[leg];
                    output += ((transition.output / term.local_stride[leg]) % model.spaces_[s].charges.size()) * model.strides_[s];
                }
                unsigned parity = 0;
                if (term.fermionic) for (auto [a,b] : term.inversions) {
                    const auto pa = [&](State x, std::size_t site) {
                        return model.spaces_[site].parity[(x / model.strides_[site]) % model.spaces_[site].charges.size()];
                    };
                    parity ^= (pa(state,a) & pa(state,b)) ^ (pa(output,a) & pa(output,b));
                }
                emit(output, parity ? -transition.value : transition.value, channel);
            }
        }
    }
};
}
