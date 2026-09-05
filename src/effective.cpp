#include <pcut/effective.hpp>
#include "detail.hpp"
#include "transitions.hpp"
#include "monomial_registry.hpp"
#include <algorithm>
#include <boost/unordered/unordered_flat_map.hpp>
#include <cmath>
#include <functional>
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
    factor_program();
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
        detail::MonomialRegistry monomials(order_,model.terms_.size(),channel_edges);
        const auto required=monomials.required_support();
        struct Transition { State out; Scalar value; std::size_t channel; };
        boost::unordered_flat_map<std::pair<int,State>,std::vector<Transition>> transition_cache;
        auto transitions=[&](int change,State in) -> const std::vector<Transition>& {
            auto [it,fresh]=transition_cache.try_emplace({change,in});
            if (fresh) detail::ModelTransitions::enumerate(model,change,in,[&](State out,Complex amplitude,std::size_t channel) {
                const auto value=amplitude/model.gap();
                if constexpr (std::is_same_v<Scalar,double>) it->second.push_back({out,value.real(),channel});
                else it->second.push_back({out,value,channel});
            });
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
                        for (const auto& [id,c] : p) if (monomials.support(id)==required)
                            accumulated[depth][out][id]+=instruction.coefficient*c;
                    for (const auto& edge : instruction.edges) {
                        auto& transformed=next[edge.target];
                        for (const auto& [in,p] : state) for (const auto& t : transitions(edge.change,in)) {
                            if (depth+1==order_ && !external.contains(t.out)) continue;
                            decltype(&transformed.begin()->second) dest=nullptr;
                            for (const auto& [m,c] : p) {
                                const auto id=monomials.multiply(m,t.channel);
                                if (id==detail::MonomialRegistry::absent) continue;
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
                for (const auto& [id,c] : p) if (c!=Scalar{}) terms.emplace_back(monomials.value(id),c);
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
