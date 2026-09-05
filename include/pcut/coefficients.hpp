#pragma once
#include <boost/multiprecision/cpp_int.hpp>
#include <cstddef>
#include <map>
#include <vector>

namespace pcut {
using Rational = boost::multiprecision::cpp_rational;
using Word = std::vector<int>;

// C(m_1,...,m_k) multiplies T_{m_1}...T_{m_k}; rightmost acts first.
class Coefficients {
public:
    explicit Coefficients(std::vector<int> changes, unsigned order);
    [[nodiscard]] unsigned order() const noexcept { return order_; }
    [[nodiscard]] const std::vector<int>& changes() const noexcept { return changes_; }
    [[nodiscard]] const std::map<Word, Rational>& terms() const noexcept { return terms_; }
    [[nodiscard]] Rational at(const Word& word) const;
private:
    unsigned order_;
    std::vector<int> changes_;
    std::map<Word, Rational> terms_;
};
}
