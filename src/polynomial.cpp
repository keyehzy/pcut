#include <pcut/polynomial.hpp>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace pcut {
Monomial::Monomial(std::initializer_list<std::pair<std::size_t,unsigned>> powers) : Monomial(Powers(powers)) {}
Monomial::Monomial(Powers powers) {
    std::sort(powers.begin(),powers.end());
    unsigned degree=0;
    bool compact=true;
    for (std::size_t i=0;i<powers.size();++i) {
        const auto [v,n]=powers[i];
        if (!n || n>64 || degree>64-n || (i && powers[i-1].first==v))
            throw std::invalid_argument("invalid monomial powers or degree above 64");
        degree+=n; compact &= v<16 && n<16;
    }
    if (compact) for (const auto& [v,n] : powers) packed_ |= static_cast<std::uint64_t>(n)<<(4*v);
    else sparse_=std::make_shared<const Powers>(std::move(powers));
}
bool Monomial::operator<(const Monomial& other) const {
    if (static_cast<bool>(sparse_)!=static_cast<bool>(other.sparse_)) return !sparse_;
    return sparse_ ? *sparse_<*other.sparse_ : packed_<other.packed_;
}
bool Monomial::operator==(const Monomial& other) const {
    if (static_cast<bool>(sparse_)!=static_cast<bool>(other.sparse_)) return false;
    return sparse_ ? *sparse_==*other.sparse_ : packed_==other.packed_;
}
unsigned Monomial::degree() const {
    unsigned result=0;
    for_each_power([&](std::size_t,unsigned n) { result+=n; });
    return result;
}
Monomial Monomial::multiplied(std::size_t variable) const {
    if (degree()==64) throw std::length_error("monomial degree exceeds 64");
    if (!sparse_ && variable<16 && ((packed_>>(4*variable))&15)<15) {
        auto result=*this; result.packed_+=std::uint64_t{1}<<(4*variable); return result;
    }
    Powers powers;
    for_each_power([&](std::size_t v,unsigned n) { powers.emplace_back(v,n); });
    auto it=std::lower_bound(powers.begin(),powers.end(),variable,[](const auto& p,auto v) { return p.first<v; });
    if (it!=powers.end() && it->first==variable) ++it->second;
    else powers.insert(it,{variable,1});
    return Monomial(std::move(powers));
}
Monomial Monomial::relabel(const std::vector<std::size_t>& variables) const {
    Monomial result;
    for_each_power([&](std::size_t v,unsigned n) { for (unsigned i=0;i<n;++i) result=result.multiplied(variables.at(v)); });
    return result;
}
void add_polynomial(Polynomial& target, const Polynomial& source, Complex scale,
                    const std::vector<std::size_t>* variables) {
    if (&target==&source) { const auto copy=source; add_polynomial(target,copy,scale,variables); return; }
    for (const auto& [m,c] : source) {
        const auto key=variables ? m.relabel(*variables) : m;
        auto& value=target[key]; value+=scale*c;
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
            throw std::overflow_error("symbolic matrix element overflow");
        if (value==Complex{}) target.erase(key);
    }
}
Complex substitute(const Polynomial& polynomial, const std::vector<double>& values) {
    for (double x : values) if (!std::isfinite(x)) throw std::invalid_argument("nonfinite coupling");
    Complex result=0;
    for (const auto& [m,c] : polynomial) {
        Complex value=c;
        m.for_each_power([&](std::size_t v,unsigned n) { value*=std::pow(values.at(v),n); });
        result+=value;
    }
    if (!std::isfinite(result.real()) || !std::isfinite(result.imag()))
        throw std::overflow_error("coupling substitution overflow");
    return result;
}
Series substitute(const SymbolicSeries& series, const std::vector<double>& values) {
    Series result;
    for (const auto& p : series) result.push_back(substitute(p,values));
    return result;
}
}
