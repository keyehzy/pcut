#include <pcut/effective.hpp>
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace pcut {
void EffectiveOperator::factor_program() {
    // Exact row-space factorization of the suffix program. A feature is either
    // a terminal coefficient or (next letter, next-layer basis vector).
    // Gaussian elimination shares linearly dependent continuations, not merely
    // identical strings. Every operation here remains a Boost rational.
    using Row=std::map<std::size_t,Rational>;
    std::vector<Row> representations(nodes_.size());
    layers_.resize(order_+1);
    for (unsigned depth=order_+1;depth-->0;) {
        const auto width=depth<order_ ? layers_[depth+1].size() : 0;
        if (!changes_.empty() && width>(std::numeric_limits<std::size_t>::max()-1)/changes_.size())
            throw std::length_error("coefficient program feature count overflow");
        std::map<std::size_t,std::size_t> pivots;
        std::vector<Row> rows;
        for (std::size_t i=0;i<nodes_.size();++i) if (nodes_[i].depth==depth) {
            Row row;
            if (nodes_[i].coefficient!=0) row[0]=nodes_[i].coefficient;
            for (const auto& [change,child] : nodes_[i].children) {
                const auto letter=static_cast<std::size_t>(std::lower_bound(changes_.begin(),changes_.end(),change)-changes_.begin());
                for (const auto& [basis,scale] : representations[child]) row[1+letter*width+basis]+=scale;
            }
            auto& representation=representations[i];
            while (!row.empty()) {
                const auto pivot=row.begin()->first;
                const Rational scale=row.begin()->second;
                if (scale==0) { row.erase(row.begin()); continue; }
                const auto found=pivots.find(pivot);
                if (found==pivots.end()) {
                    const auto basis=rows.size(); pivots.emplace(pivot,basis);
                    for (auto& [feature,c] : row) { (void)feature; c/=scale; }
                    representation[basis]+=scale;
                    rows.push_back(std::move(row)); break;
                }
                representation[found->second]+=scale;
                for (const auto& [feature,c] : rows[found->second]) {
                    auto& value=row[feature]; value-=scale*c;
                    if (value==0) row.erase(feature);
                }
            }
        }
        for (const auto& row : rows) {
            FactoredNode node;
            for (const auto& [feature,c] : row) {
                if (feature==0) node.coefficient=c;
                else node.edges.push_back({changes_[(feature-1)/width],(feature-1)%width,c});
            }
            layers_[depth].push_back(std::move(node));
        }
    }
    if (!representations[0].empty()) root_scale_=representations[0].begin()->second;
}
}
