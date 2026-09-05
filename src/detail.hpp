#pragma once
#include <pcut/lattice.hpp>
#include <algorithm>
#include <limits>
#include <set>
#include <stdexcept>

namespace pcut::detail {
inline std::size_t matrix_size(std::size_t basis, std::size_t orders, std::size_t budget) {
    const auto index_max=static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max());
    budget=std::min(budget,std::numeric_limits<std::size_t>::max()/sizeof(Complex));
    if (!orders || basis>index_max || (basis &&
        (basis>index_max/basis || basis>budget/basis || orders>budget/basis/basis)))
        throw std::length_error("matrix-element budget or Eigen indexing exceeded");
    return basis*basis*orders;
}
inline std::vector<Matrix> matrix_series(std::size_t basis, std::size_t orders, std::size_t budget) {
    matrix_size(basis,orders,budget);
    std::vector<Matrix> result;
    result.reserve(orders);
    for (std::size_t n=0;n<orders;++n)
        result.emplace_back(Matrix::Zero(static_cast<Eigen::Index>(basis),static_cast<Eigen::Index>(basis)));
    return result;
}
struct StorageBudget {
    std::size_t remaining;
    void take(std::size_t count) {
        if (count>remaining) throw std::length_error("linked coefficient storage budget exceeded");
        remaining-=count;
    }
    void matrices(std::size_t basis, std::size_t orders) {
        take(matrix_size(basis,orders,remaining));
    }
};
inline Coordinate translate(const Coordinate& a, const Coordinate& b, int sign=1) {
    if (a.size()!=b.size()) throw std::invalid_argument("coordinate dimension mismatch");
    Coordinate out(a.size());
    for (std::size_t d=0;d<a.size();++d) {
        const auto x=static_cast<long long>(a[d])+static_cast<long long>(sign)*b[d];
        if (x<std::numeric_limits<int>::min() || x>std::numeric_limits<int>::max())
            throw std::overflow_error("coordinate overflow");
        out[d]=static_cast<int>(x);
    }
    return out;
}
inline void subtract_mapped(std::vector<Matrix>& parent, const std::vector<Matrix>& child,
                            const std::vector<Eigen::Index>& map) {
    for (std::size_t n=1;n<parent.size();++n)
        for (std::size_t i=0;i<map.size();++i) for (std::size_t j=0;j<map.size();++j)
            parent[n](map[i],map[j])-=child[n](static_cast<Eigen::Index>(i),static_cast<Eigen::Index>(j));
}
inline std::vector<std::pair<std::size_t,std::size_t>> gather_inversions(
    std::size_t size, const std::vector<std::size_t>& support) {
    std::vector<std::size_t> order=support;
    std::set<std::size_t> used;
    for (auto s : support) if (s>=size || !used.insert(s).second)
        throw std::invalid_argument("invalid graded embedding support");
    for (std::size_t s=0;s<size;++s) if (!used.contains(s)) order.push_back(s);
    std::vector<std::pair<std::size_t,std::size_t>> inversions;
    for (std::size_t i=0;i<size;++i) for (std::size_t j=i+1;j<size;++j)
        if (order[i]>order[j]) inversions.emplace_back(order[i],order[j]);
    return inversions;
}
inline int gather_sign(const std::vector<LocalSpace>& spaces, const std::vector<unsigned>& local,
                       const std::vector<std::pair<std::size_t,std::size_t>>& inversions) {
    unsigned parity=0;
    for (auto [a,b] : inversions) parity^=spaces[a].parity[local[a]] & spaces[b].parity[local[b]];
    return parity ? -1 : 1;
}
}
