#pragma once
#include <pcut/model.hpp>
#include <compare>
#include <bit>
#include <boost/container/small_vector.hpp>
#include <boost/container/flat_map.hpp>

namespace pcut {
// Immutable sparse exponents. Sixteen four-bit exponent lanes fit inline;
// larger indices/exponents use a shared sparse vector, with no loss of range.
class Monomial {
public:
    using Powers=boost::container::small_vector<std::pair<std::size_t,unsigned>,4>;
    Monomial() = default;
    explicit Monomial(std::initializer_list<std::pair<std::size_t,unsigned>> powers);
    bool operator<(const Monomial& other) const;
    bool operator==(const Monomial& other) const;
    [[nodiscard]] unsigned degree() const;
    [[nodiscard]] Monomial multiplied(std::size_t variable) const;
    [[nodiscard]] Monomial relabel(const std::vector<std::size_t>& variables) const;
    template<class Function> void for_each_power(Function&& function) const {
        if (sparse_) { for (const auto& [v,n] : *sparse_) function(v,n); }
        else {
            auto bits=packed_;
            while (bits) {
                const auto variable=static_cast<unsigned>(std::countr_zero(bits))/4;
                function(static_cast<std::size_t>(variable),static_cast<unsigned>((bits>>(4*variable))&15));
                bits &= ~(std::uint64_t{15}<<(4*variable));
            }
        }
    }
private:
    explicit Monomial(Powers powers);
    std::uint64_t packed_=0;
    std::shared_ptr<const Powers> sparse_;
};
using Polynomial = boost::container::flat_map<Monomial, Complex>;
using SymbolicSeries = std::vector<Polynomial>;
struct SymbolicBlock {
    std::vector<State> basis;
    // Sparse (output state, input state) matrix for each lambda order.
    std::vector<std::map<std::pair<State, State>, Polynomial>> coefficients;
};
[[nodiscard]] Complex substitute(const Polynomial& polynomial, const std::vector<double>& values);
[[nodiscard]] Series substitute(const SymbolicSeries& series, const std::vector<double>& values);
void add_polynomial(Polynomial& target, const Polynomial& source, Complex scale = 1,
                    const std::vector<std::size_t>* variable_map = nullptr);
}
