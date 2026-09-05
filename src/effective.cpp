#include <pcut/effective.hpp>
#include "detail.hpp"
#include <algorithm>
#include <bit>
#include <boost/unordered/unordered_flat_map.hpp>
#include <cmath>
#include <functional>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>

namespace pcut {
EffectiveOperator::EffectiveOperator(const Coefficients& coefficients)
    : order_(coefficients.order()), changes_(coefficients.changes()), nodes_(1) {
    std::ostringstream identity;
    identity << order_ << ':';
    for (int change : changes_) identity << change << ',';
    for (const auto& [word, c] : coefficients.terms()) {
        identity << ';';
        for (int change : word) identity << change << ',';
        identity << '=' << c;
        std::size_t node = 0;
        unsigned depth = 0;
        for (auto it = word.rbegin(); it != word.rend(); ++it) {
            ++depth;
            const auto existing = nodes_[node].children.find(*it);
            if (existing != nodes_[node].children.end()) node = existing->second;
            else {
                const auto next = nodes_.size();
                nodes_[node].children.emplace(*it, next);
                nodes_.push_back({{}, 0, depth});
                node = next;
            }
        }
        nodes_[node].coefficient = c;
    }
    identity_=identity.str();
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
std::vector<SparseState> EffectiveOperator::apply(const ClusterModel& model, State input) const {
    if (input >= model.dimension()) throw std::out_of_range("effective input state");
    for (int m : model.changes()) if (!std::binary_search(changes_.begin(), changes_.end(), m))
        throw std::invalid_argument("coefficient alphabet does not cover model charge changes");
    std::vector<SparseState> output(order_ + 1);
    output[0][input] = model.vacuum_energy() + model.gap() * model.charge(input);
    std::function<void(std::size_t,const SparseState&)> visit = [&](std::size_t index, const SparseState& state) {
        const auto& node = nodes_[index];
        if (node.coefficient != 0) {
            const double coefficient=node.coefficient.convert_to<double>();
            for (const auto& [s, v] : state) output[node.depth][s] += coefficient * v;
        }
        for (const auto& [change, next] : node.children) {
            auto transformed = model.apply_scaled(change, state, model.gap());
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
std::vector<Matrix> EffectiveOperator::block(const ClusterModel& model, const std::vector<State>& basis) const {
    const std::set<State> unique(basis.begin(), basis.end());
    if (unique.size() != basis.size()) throw std::invalid_argument("duplicate block basis state");
    for (auto state : basis) if (state>=model.dimension()) throw std::out_of_range("block basis state");
    auto result=detail::matrix_series(basis.size(),static_cast<std::size_t>(order_)+1);
    for (std::size_t col = 0; col < basis.size(); ++col) {
        const auto action = apply(model, basis[col]);
        for (unsigned n = 0; n <= order_; ++n)
            for (std::size_t row = 0; row < basis.size(); ++row) {
                const auto it = action[n].find(basis[row]);
                if (it != action[n].end()) result[n](static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col)) = it->second;
            }
    }
    return result;
}
Series EffectiveOperator::vacuum(const ClusterModel& model) const {
    model.require_product_vacuum();
    const auto action = apply(model, 0);
    Series result(order_+1);
    for (unsigned n = 0; n <= order_; ++n) {
        const auto it = action[n].find(0);
        if (it != action[n].end()) result[n] = it->second;
    }
    return result;
}
SymbolicBlock EffectiveOperator::symbolic_block(const ClusterModel& model, const std::vector<State>& basis,
    const std::vector<std::size_t>& channel_edges) const {
    auto evaluate=[&]<class Scalar>() -> SymbolicBlock {
        auto finite=[](Scalar value) {
            if constexpr (std::is_same_v<Scalar,double>) return std::isfinite(value);
            else return std::isfinite(value.real()) && std::isfinite(value.imag());
        };
        std::uint64_t required=0;
        if (!channel_edges.empty()) {
            if (channel_edges.size()!=model.terms_.size()) throw std::invalid_argument("channel edge map size");
            for (auto e : channel_edges) {
                if (e>=64) throw std::invalid_argument("edge support exceeds 64");
                required |= std::uint64_t{1}<<e;
            }
        }
        // Intern exponent vectors once. Most programs fit losslessly in a uint64;
        // larger alphabets use the same exact registry with compact sparse keys.
        const auto bits=static_cast<unsigned>(std::bit_width(order_));
        const bool packed=bits==0 || model.terms_.size()<=64/bits;
        struct MonomialInfo { Monomial value; std::uint64_t code, support; unsigned degree; std::vector<std::size_t> products; };
        std::vector<MonomialInfo> monomials{{{},0,0,0,{}}};
        boost::unordered_flat_map<std::uint64_t,std::size_t> packed_ids{{0,0}};
        std::map<Monomial,std::size_t> sparse_ids{{{},0}};
        const auto absent=std::numeric_limits<std::size_t>::max();
        auto multiply=[&](std::size_t id,std::size_t channel) {
            if (monomials[id].products.empty()) monomials[id].products.resize(model.terms_.size(),absent);
            if (monomials[id].products[channel]!=absent) return monomials[id].products[channel];
            if (monomials[id].degree>=64) throw std::length_error("monomial degree exceeds 64");
            const auto mask=monomials[id].support | (required ? std::uint64_t{1}<<channel_edges[channel] : 0);
            const auto degree=monomials[id].degree+1;
            if (static_cast<unsigned>(std::popcount(required & ~mask))>order_-degree) return absent;
            std::size_t next;
            if (packed) {
                const auto code=monomials[id].code+(std::uint64_t{1}<<(bits*channel));
                auto [it,fresh]=packed_ids.try_emplace(code,monomials.size()); next=it->second;
                if (fresh) monomials.push_back({monomials[id].value.multiplied(channel),code,mask,degree,{}});
            } else {
                auto value=monomials[id].value.multiplied(channel);
                auto [it,fresh]=sparse_ids.try_emplace(value,monomials.size()); next=it->second;
                if (fresh) monomials.push_back({std::move(value),0,mask,degree,{}});
            }
            monomials[id].products[channel]=next;
            return next;
        };
        struct Transition { State out; Scalar value; std::size_t channel; };
        boost::unordered_flat_map<std::pair<int,State>,std::vector<Transition>> transition_cache;
        auto transitions=[&](int change,State in) -> const std::vector<Transition>& {
            auto [it,fresh]=transition_cache.try_emplace({change,in});
            if (fresh) for (std::size_t c=0;c<model.terms_.size();++c) {
                const auto& term=model.terms_[c];
                const auto block=term.by_change->find(change);
                if (block==term.by_change->end()) continue;
                State column=0,removed=0;
                for (std::size_t leg=0;leg<term.sites.size();++leg) {
                    const auto site=term.sites[leg];
                    const auto local=(in/model.strides_[site])%model.spaces_[site].charges.size();
                    column+=local*term.local_stride[leg]; removed+=local*model.strides_[site];
                }
                for (const auto& transition : block->second[column]) {
                    State out=in-removed;
                    for (std::size_t leg=0;leg<term.sites.size();++leg) {
                        const auto site=term.sites[leg];
                        out+=((transition.output/term.local_stride[leg])%model.spaces_[site].charges.size())*model.strides_[site];
                    }
                    unsigned parity=0;
                    if (term.fermionic) for (const auto& [a,b] : term.inversions) {
                        const auto pa=[&](State state,std::size_t site) { return model.spaces_[site].parity[(state/model.strides_[site])%model.spaces_[site].charges.size()]; };
                        parity^=(pa(in,a)&pa(in,b))^(pa(out,a)&pa(out,b));
                    }
                    const auto value=(parity ? -transition.value : transition.value)/model.gap();
                    if constexpr (std::is_same_v<Scalar,double>) it->second.push_back({out,value.real(),c});
                    else it->second.push_back({out,value,c});
                }
            }
            return it->second;
        };
        using Scratch=boost::unordered_flat_map<State,boost::unordered_flat_map<std::size_t,Scalar>>;
        struct NumericalEdge { int change; std::size_t target; double scale; };
        struct NumericalNode { double coefficient; std::vector<NumericalEdge> edges; };
        std::vector<std::vector<NumericalNode>> program;
        for (const auto& layer : layers_) {
            program.emplace_back();
            for (const auto& node : layer) {
                NumericalNode out{node.coefficient.convert_to<double>(),{}};
                for (const auto& edge : node.edges) out.edges.push_back({edge.change,edge.target,edge.scale.convert_to<double>()});
                program.back().push_back(std::move(out));
            }
        }
        const std::set<State> external(basis.begin(),basis.end());
        if (external.size()!=basis.size()) throw std::invalid_argument("duplicate symbolic external state");
        for (int m : model.changes()) if (!std::binary_search(changes_.begin(),changes_.end(),m))
            throw std::invalid_argument("coefficient alphabet does not cover operator channels");
        SymbolicBlock output{basis,{}};
        output.coefficients.resize(order_+1);
        for (State input : basis) {
            transition_cache.clear();
            if (input>=model.dimension()) throw std::out_of_range("symbolic external state");
            if (!required) output.coefficients[0][{input,input}][Monomial{}]=model.vacuum_energy()+model.gap()*model.charge(input);
            std::vector<Scratch> accumulated(order_+1);
            std::vector<Scratch> current(program[0].size());
            if (!current.empty()) current[0][input][0]=root_scale_.convert_to<double>();
            for (unsigned depth=0;depth<=order_;++depth) {
                std::vector<Scratch> next(depth<order_ ? program[depth+1].size() : 0);
                for (std::size_t node=0;node<current.size();++node) {
                    const auto& state=current[node];
                    const auto& instruction=program[depth][node];
                    if (instruction.coefficient!=0) for (const auto& [out,p] : state) if (external.contains(out))
                        for (const auto& [id,c] : p) if (monomials[id].support==required)
                            accumulated[depth][out][id]+=instruction.coefficient*c;
                    for (const auto& edge : instruction.edges) {
                        auto& transformed=next[edge.target];
                        for (const auto& [in,p] : state) for (const auto& t : transitions(edge.change,in)) {
                            if (depth+1==order_ && !external.contains(t.out)) continue;
                            decltype(&transformed.begin()->second) dest=nullptr;
                            for (const auto& [m,c] : p) {
                                const auto id=multiply(m,t.channel);
                                if (id==absent) continue;
                                if (!dest) dest=&transformed[t.out];
                                auto& target=(*dest)[id]; target+=edge.scale*t.value*c;
                                if (!finite(target))
                                    throw std::overflow_error("symbolic intermediate amplitude overflow");
                            }
                        }
                    }
                    // Release processed states before advancing the layer.
                    current[node]=Scratch{};
                }
                for (auto& transformed : next) for (auto it=transformed.begin();it!=transformed.end();) {
                    auto element=it++;
                    auto& p=element->second;
                    for (auto jt=p.begin();jt!=p.end();) { auto term=jt++; if (term->second==Scalar{}) p.erase(term); }
                    if (p.empty()) transformed.erase(element);
                }
                current=std::move(next);
            }
            for (unsigned n=1;n<=order_;++n) for (const auto& [out,p] : accumulated[n]) {
                std::vector<std::pair<Monomial,Complex>> terms;
                terms.reserve(p.size());
                for (const auto& [id,c] : p) if (c!=Scalar{}) terms.emplace_back(monomials[id].value,c);
                std::sort(terms.begin(),terms.end(),[](const auto& a,const auto& b) { return a.first<b.first; });
                if (!terms.empty()) output.coefficients[n][{out,input}]=Polynomial(boost::container::ordered_unique_range,terms.begin(),terms.end());
            }

        }
        // Restore energy units only after summing dimensionless model products.
        for (unsigned n=1;n<=order_;++n) for (auto& [states,p] : output.coefficients[n]) {
            (void)states;
            for (auto& [monomial,value] : p) {
                (void)monomial; value*=model.gap();
                if (!std::isfinite(value.real()) || !std::isfinite(value.imag()))
                    throw std::overflow_error("symbolic pCUT energy scale overflow");
            }
        }
        return output;
    };
    bool real=true;
    for (const auto& term : model.terms_) for (const auto& [change,columns] : *term.by_change) {
        (void)change;
        for (const auto& column : columns) for (const auto& t : column) real &= t.value.imag()==0;
    }
    if (real) return evaluate.template operator()<double>();
    return evaluate.template operator()<Complex>();
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
