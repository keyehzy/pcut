#pragma once
#include <pcut/polynomial.hpp>
#include <bit>
#include <boost/unordered/unordered_flat_map.hpp>
#include <limits>
#include <stdexcept>

namespace pcut::detail {
// Exact exponent interning. Packed keys are lossless for this order/alphabet;
// wider alphabets use sparse Monomial keys. Support pruning is purely algebraic.
class MonomialRegistry {
public:
    static constexpr auto absent=std::numeric_limits<std::size_t>::max();
    MonomialRegistry(unsigned order, std::size_t channels, const std::vector<std::size_t>& channel_edges)
        : order_(order), channels_(channels), channel_edges_(channel_edges),
          bits_(static_cast<unsigned>(std::bit_width(order))), packed_(bits_==0 || channels<=64/bits_) {
        if (!channel_edges.empty()) {
            if (channel_edges.size()!=channels) throw std::invalid_argument("channel edge map size");
            for (auto e : channel_edges) {
                if (e>=64) throw std::invalid_argument("edge support exceeds 64");
                required_ |= std::uint64_t{1}<<e;
            }
        }
    }
    [[nodiscard]] std::uint64_t required_support() const { return required_; }
    [[nodiscard]] std::uint64_t support(std::size_t id) const { return monomials[id].support; }
    [[nodiscard]] const Monomial& value(std::size_t id) const { return monomials[id].value; }
    std::size_t multiply(std::size_t id, std::size_t channel) {
        if (monomials[id].products.empty()) monomials[id].products.resize(channels_,absent);
        if (monomials[id].products[channel]!=absent) return monomials[id].products[channel];
        if (monomials[id].degree>=64) throw std::length_error("monomial degree exceeds 64");
        const auto mask=monomials[id].support | (required_ ? std::uint64_t{1}<<channel_edges_[channel] : 0);
        const auto degree=monomials[id].degree+1;
        if (static_cast<unsigned>(std::popcount(required_ & ~mask))>order_-degree) return absent;
        std::size_t next;
        if (packed_) {
            const auto code=monomials[id].code+(std::uint64_t{1}<<(bits_*channel));
            auto [it,fresh]=packed_ids.try_emplace(code,monomials.size()); next=it->second;
            if (fresh) monomials.push_back({monomials[id].value.multiplied(channel),code,mask,degree,{}});
        } else {
            auto value=monomials[id].value.multiplied(channel);
            auto [it,fresh]=sparse_ids.try_emplace(value,monomials.size()); next=it->second;
            if (fresh) monomials.push_back({std::move(value),0,mask,degree,{}});
        }
        monomials[id].products[channel]=next;
        return next;
    }
private:
    unsigned order_;
    std::size_t channels_;
    const std::vector<std::size_t>& channel_edges_;
    unsigned bits_;
    bool packed_;
    std::uint64_t required_=0;
    struct MonomialInfo { Monomial value; std::uint64_t code, support; unsigned degree; std::vector<std::size_t> products; };
    std::vector<MonomialInfo> monomials{{{},0,0,0,{}}};
    boost::unordered_flat_map<std::uint64_t,std::size_t> packed_ids{{0,0}};
    std::map<Monomial,std::size_t> sparse_ids{{{},0}};
};
}
