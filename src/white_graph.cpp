#include <pcut/white_graph.hpp>
#include "detail.hpp"
#include "incidence.hpp"
#include "graph_structure.hpp"
#include <boost/multiprecision/cpp_int.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace pcut {
namespace {
using detail::number;
using detail::real;
using detail::string;
using detail::sequence;
using detail::space_key;
using detail::channel_key;
using detail::channel_order;
using detail::edge_structure;
std::string evaluation_key(const WhiteGraph& graph) {
    std::string key; number(key,graph.spaces.size());
    for (const auto& space : graph.spaces) string(key,space_key(space));
    number(key,graph.edges.size());
    for (const auto& edge : graph.edges) {
        sequence(key,edge.legs); number(key,edge.channels.size());
        for (const auto& channel : edge.channels) string(key,channel_key(channel));
    }
    return key;
}
// Exact interning is local to a construction. IDs accelerate serialization and
// memoization; persistent graph/cache identities still contain full structures.
struct Signatures {
    struct EdgeInfo { GraphEdge edge; std::string key; std::vector<std::size_t> order; };
    std::vector<std::pair<LocalSpace,std::string>> spaces;
    std::vector<EdgeInfo> edges;
    std::size_t space(const LocalSpace& s) {
        for (std::size_t i=0;i<spaces.size();++i) {
            const auto& a=spaces[i].first;
            if (a.name==s.name && a.charges==s.charges && a.particles==s.particles && a.parity==s.parity && a.vacuum_energy==s.vacuum_energy) return i;
        }
        spaces.emplace_back(s,space_key(s)); return spaces.size()-1;
    }
    std::size_t edge(const GraphEdge& e) {
        for (std::size_t i=0;i<edges.size();++i) {
            const auto& a=edges[i].edge;
            if (a.legs.size()!=e.legs.size() || a.channels.size()!=e.channels.size()) continue;
            bool equal=true;
            for (std::size_t c=0;c<e.channels.size() && equal;++c) {
                const auto& x=a.channels[c]; const auto& y=e.channels[c];
                equal=x.fermionic==y.fermionic && x.matrix.rows()==y.matrix.rows() && x.matrix.cols()==y.matrix.cols() && x.matrix==y.matrix;
            }
            if (equal) return i;
        }
        edges.push_back({e,edge_structure(e),channel_order(e)}); return edges.size()-1;
    }
};
bool graph_connected(const WhiteGraph& g) {
    if (g.edges.empty()) return false;
    std::set<std::size_t> seen(g.edges.front().legs.begin(),g.edges.front().legs.end());
    for (std::size_t pass=0;pass<g.edges.size();++pass) for (const auto& e : g.edges)
        if (std::any_of(e.legs.begin(),e.legs.end(),[&](auto v) { return seen.contains(v); }))
            seen.insert(e.legs.begin(),e.legs.end());
    return seen.size()==g.spaces.size();
}
std::vector<std::size_t> offsets(const WhiteGraph& g) {
    std::vector<std::size_t> result{0};
    for (const auto& e : g.edges) {
        if (e.channels.size()>std::numeric_limits<std::size_t>::max()-result.back()) throw std::length_error("channel count overflow");
        result.push_back(result.back()+e.channels.size());
    }
    return result;
}

}
ClusterModel WhiteGraph::model(double gap) const {
    std::vector<LocalTerm> terms;
    for (const auto& edge : edges) for (const auto& channel : edge.channels)
        terms.push_back({edge.legs,channel.matrix,channel.fermionic});
    return {spaces,std::move(terms),gap};
}
std::size_t WhiteGraph::variables() const { return offsets(*this).back(); }
static CanonicalGraph canonicalize_impl(const WhiteGraph& g, bool validate, Signatures& pool) {
    if (g.edges.empty() || g.edges.size()>64 || g.spaces.empty() || g.spaces.size()>1000)
        throw std::invalid_argument("invalid abstract graph size");
    if (validate) {
    // Validate local operators on their support without constructing a full graph tensor space.
    for (const auto& space : g.spaces) space.validate();
    for (const auto& edge : g.edges) {
        std::vector<LocalSpace> local; std::vector<std::size_t> legs; std::set<std::size_t> unique;
        if (edge.channels.empty()) throw std::invalid_argument("edge needs channels");
        for (auto v : edge.legs) {
            if (v>=g.spaces.size() || !unique.insert(v).second) throw std::invalid_argument("invalid abstract leg");
            legs.push_back(local.size()); local.push_back(g.spaces[v]);
        }
        for (const auto& c : edge.channels) { const ClusterModel check(local,{{legs,c.matrix,c.fermionic}}); (void)check; }
    }
    }
    if (!graph_connected(g)) throw std::invalid_argument("abstract graph must be connected without isolated vertices");
    std::vector<std::string> labels;
    std::vector<std::size_t> space_ids;
    for (const auto& s : g.spaces) { const auto id=pool.space(s); space_ids.push_back(id); labels.push_back(pool.spaces[id].second); }
    std::vector<std::string> structures;
    std::vector<std::vector<std::size_t>> orders;
    for (const auto& e : g.edges) { const auto id=pool.edge(e); structures.push_back(pool.edges[id].key); orders.push_back(pool.edges[id].order); }
    const auto incidence=detail::encode_incidence(g,labels,structures);
    // Restriction to physical vertices is onto. Its kernel only permutes
    // identical ordered edge occurrences. Channels live in the edge colors.
    // Count the kernel exactly; only the quotient must fit size_t.
    boost::multiprecision::cpp_int kernel=1;
    std::map<std::string,std::size_t> edge_multiplicities;
    for (std::size_t e=0;e<g.edges.size();++e) {
        auto key=structures[e]; sequence(key,g.edges[e].legs);
        kernel*=++edge_multiplicities[key];
    }
    const auto labeling=detail::label_incidence(incidence);
    std::vector<std::size_t> permutation;
    for (int node : labeling.canonical_to_input)
        if (static_cast<std::size_t>(node)<g.spaces.size()) permutation.push_back(static_cast<std::size_t>(node));
    boost::multiprecision::cpp_int group=1;
    for (int index : labeling.group_indices) group*=index;
    if (group%kernel!=0) throw std::logic_error("incidence automorphism kernel mismatch");
    group/=kernel;
    if (group>std::numeric_limits<std::size_t>::max()) throw std::length_error("automorphism count overflow");
    const auto source_offsets=offsets(g);
    std::vector<std::size_t> inverse(g.spaces.size());
    CanonicalGraph candidate;
    candidate.map.vertices=permutation;
    number(candidate.key,g.spaces.size());
    for (std::size_t v=0;v<permutation.size();++v) {
        inverse[permutation[v]]=v; candidate.graph.spaces.push_back(g.spaces[permutation[v]]);
        string(candidate.key,pool.spaces[space_ids[permutation[v]]].second);
    }
    // Once vertex labels are canonical, exact structural sorting recovers
    // edge/channel maps. Equal occurrences use their input index as a tie-break;
    // it affects only the occurrence map, never the serialized identity.
    std::vector<std::pair<std::string,std::size_t>> edges;
    for (std::size_t e=0;e<g.edges.size();++e) {
        auto key=structures[e];
        for (auto v : g.edges[e].legs) number(key,inverse[v]);
        edges.emplace_back(std::move(key),e);
    }
    std::sort(edges.begin(),edges.end()); number(candidate.key,edges.size());
    for (const auto& [key,e] : edges) {
        string(candidate.key,key); candidate.map.edges.push_back(e);
        GraphEdge edge;
        for (auto v : g.edges[e].legs) edge.legs.push_back(inverse[v]);
        for (auto c : orders[e]) { edge.channels.push_back(g.edges[e].channels[c]); candidate.map.channels.push_back(source_offsets[e]+c); }
        candidate.graph.edges.push_back(std::move(edge));
    }
    candidate.vertex_automorphisms=group.convert_to<std::size_t>();
    return candidate;
}
CanonicalGraph canonicalize(const WhiteGraph& graph) { Signatures pool; return canonicalize_impl(graph,true,pool); }
struct WhiteGraphExpansion::Plan {
    PeriodicLattice lattice;
    mutable std::once_flag graph_subclusters_once, embedding_subclusters_once;
    mutable std::vector<std::vector<GraphSubcluster>> graph_subclusters;
    mutable std::vector<std::vector<EmbeddedSubcluster>> embedding_subclusters;
    std::vector<CanonicalGraph> graphs;
    std::vector<GraphEmbedding> embeddings;
    std::vector<ClusterModel> canonical_models, physical_models;
    struct Readout {
        std::shared_ptr<const SymbolicBlock> block;
        struct Entry { Eigen::Index row, column; int sign; const Polynomial* polynomial; };
        std::vector<std::vector<Entry>> entries;
    };
    mutable std::mutex mutex;
    using BasisKey=std::tuple<std::size_t,std::vector<State>,bool>;
    mutable std::map<std::string,std::map<BasisKey,Readout>> readouts;
};
const std::vector<CanonicalGraph>& WhiteGraphExpansion::graphs() const noexcept { return plan_->graphs; }
const std::vector<GraphEmbedding>& WhiteGraphExpansion::embeddings() const noexcept { return plan_->embeddings; }
const PeriodicLattice& WhiteGraphExpansion::structure() const noexcept { return plan_->lattice; }
PeriodicLattice WhiteGraphExpansion::bound_lattice() const {
    auto result=structure();
    for (std::size_t t=0;t<ratios_.size();++t)
        for (std::size_t c=0;c<ratios_[t].size();++c) result.interactions[t].channels[c].coupling=ratios_[t][c];
    return result;
}
WhiteGraphExpansion WhiteGraphExpansion::bind(const std::vector<std::vector<double>>& ratios) const {
    return {plan_,max_edges_,cache_,ratios};
}
WhiteGraphExpansion::WhiteGraphExpansion(std::shared_ptr<const Plan> plan,unsigned max_edges,
    std::shared_ptr<GraphCache> cache,const std::vector<std::vector<double>>& ratios)
    : max_edges_(max_edges), cache_(std::move(cache)), plan_(std::move(plan)), ratios_(ratios) {
    if (ratios.size()!=structure().interactions.size()) throw std::invalid_argument("interaction coupling count mismatch");
    for (std::size_t t=0;t<ratios.size();++t) {
        if (ratios[t].size()!=structure().interactions[t].channels.size()) throw std::invalid_argument("channel coupling count mismatch");
        for (double x : ratios[t]) if (!std::isfinite(x)) throw std::invalid_argument("nonfinite coupling");
    }
    couplings_.reserve(embeddings().size());
    for (const auto& e : embeddings()) {
        std::vector<double> physical, canonical;
        for (const auto& edge : e.edges) for (double x : ratios[edge.type]) physical.push_back(x);
        for (auto c : e.map.channels) canonical.push_back(physical[c]);
        couplings_.push_back(std::move(canonical));
    }
}
const std::vector<double>& WhiteGraphExpansion::couplings(std::size_t index) const {
    return couplings_.at(index);
}
const ClusterModel& WhiteGraphExpansion::structural_model(std::size_t index) const {
    return plan_->physical_models.at(index);
}
WhiteGraphExpansion::WhiteGraphExpansion(PeriodicLattice lattice,unsigned max_edges,std::shared_ptr<GraphCache> cache)
    : max_edges_(max_edges), cache_(std::move(cache)) {
    if (max_edges>64 || !cache_) throw std::invalid_argument("invalid white graph order or cache");
    lattice.validate();
    const auto ratios=lattice.couplings();
    for (auto& interaction : lattice.interactions) for (auto& c : interaction.channels) c.coupling=1;
    auto plan=std::make_shared<Plan>();
    plan->lattice=std::move(lattice);
    const auto& lattice_=plan->lattice;
    auto& graphs_=plan->graphs;
    auto& embeddings_=plan->embeddings;
    Signatures signatures;
    std::map<std::string,CanonicalGraph> canonical_memo;
    auto canonicalize_graph=[&](const WhiteGraph& graph) {
        std::string key;
        number(key,graph.spaces.size());
        for (const auto& space : graph.spaces) number(key,signatures.space(space));
        number(key,graph.edges.size());
        for (const auto& edge : graph.edges) { number(key,signatures.edge(edge)); sequence(key,edge.legs); }
        auto it=canonical_memo.find(key);
        if (it==canonical_memo.end()) it=canonical_memo.emplace(std::move(key),canonicalize_impl(graph,false,signatures)).first;
        return it->second;
    };
    // Grow translation classes of injective physical edge sets as embedding
    // witnesses. Canonicalize their uncolored structures before any evaluation.
    // Every connected set has a connected deletion, so growth is exhaustive on
    // the infinite lattice, including duplicate templates and hyperedges.
    std::map<std::string,std::size_t> indices;
    std::set<Cluster> level;
    if (max_edges) for (std::size_t t=0;t<lattice_.interactions.size();++t)
        level.insert({{t,Coordinate(lattice_.dimension,0)}});
    for (unsigned size=1;size<=max_edges && !level.empty();++size) {
        for (const auto& edges : level) {
            GraphEmbedding embedding;
            embedding.edges=edges; embedding.sites=vertices(lattice_,edges);
            WhiteGraph graph;
            for (const auto& site : embedding.sites) graph.spaces.push_back(lattice_.cell[site.basis]);
            for (const auto& edge : edges) {
                GraphEdge abstract;
                const auto& interaction=lattice_.interactions[edge.type];
                for (auto site : interaction.legs) {
                    site.cell=detail::translate(site.cell,edge.origin);
                    abstract.legs.push_back(static_cast<std::size_t>(std::lower_bound(embedding.sites.begin(),embedding.sites.end(),site)-embedding.sites.begin()));
                }
                for (const auto& c : interaction.channels) abstract.channels.push_back(c.op);
                graph.edges.push_back(std::move(abstract));
            }
            auto canonical=canonicalize_graph(graph);
            auto [it,fresh]=indices.emplace(canonical.key,graphs_.size());
            embedding.graph=it->second; embedding.map=canonical.map;
            if (fresh) graphs_.push_back(std::move(canonical));
            embeddings_.push_back(std::move(embedding));
        }
        if (size==max_edges) break;
        std::set<Cluster> next;
        for (const auto& edges : level) for (const auto& site : vertices(lattice_,edges))
            for (std::size_t t=0;t<lattice_.interactions.size();++t) for (const auto& leg : lattice_.interactions[t].legs)
                if (site.basis==leg.basis) {
                    const Edge edge{t,detail::translate(site.cell,leg.cell,-1)};
                    if (std::binary_search(edges.begin(),edges.end(),edge)) continue;
                    auto grown=edges; grown.push_back(edge); next.insert(normalize(std::move(grown)).cluster);
                }
        level=std::move(next);
    }
    std::sort(embeddings_.begin(),embeddings_.end(),[](const auto& a,const auto& b) {
        return a.edges.size()==b.edges.size() ? a.edges<b.edges : a.edges.size()<b.edges.size();
    });
    for (const auto& graph : graphs_) plan->canonical_models.push_back(graph.graph.model(lattice_.gap));
    for (const auto& e : embeddings_)
        plan->physical_models.push_back(plan->canonical_models[e.graph].reordered(e.map.vertices));
    plan_=std::move(plan);
    *this=bind(ratios);
}
const std::vector<GraphSubcluster>& WhiteGraphExpansion::graph_subclusters(std::size_t index) const {
    (void)graphs().at(index);
    std::call_once(plan_->graph_subclusters_once,[&] {
        const auto& graphs_=graphs();
        std::vector<std::vector<GraphSubcluster>> result(graphs_.size());
        std::map<std::string,std::size_t> indices;
        for (std::size_t i=0;i<graphs_.size();++i) indices.emplace(graphs_[i].key,i);
        Signatures signatures;
        auto canonicalize_graph=[&](const WhiteGraph& graph) { return canonicalize_impl(graph,false,signatures); };
        for (std::size_t i=0;i<graphs_.size();++i) {
            const auto& entry=graphs_[i];
            const auto& graph=entry.graph;
            const auto parent_offsets=offsets(graph);
            std::set<std::vector<std::size_t>> visited;
            std::vector<std::vector<std::size_t>> pending(1);
            pending[0].resize(graph.edges.size()); std::iota(pending[0].begin(),pending[0].end(),0);
            for (std::size_t pos=0;pos<pending.size();++pos) {
                const auto subset=pending[pos];
                if (subset.size()==1) continue;
                for (std::size_t erase=0;erase<subset.size();++erase) {
                    auto selected=subset; selected.erase(selected.begin()+static_cast<std::ptrdiff_t>(erase));
                    if (!visited.insert(selected).second) continue;
                    std::set<std::size_t> vertices;
                    for (auto e : selected) vertices.insert(graph.edges[e].legs.begin(),graph.edges[e].legs.end());
                    std::vector<std::size_t> vertex_map(vertices.begin(),vertices.end()), variable_map;
                    WhiteGraph child;
                    for (auto v : vertex_map) child.spaces.push_back(graph.spaces[v]);
                    for (auto e : selected) {
                        auto edge=graph.edges[e];
                        for (auto& v : edge.legs) v=static_cast<std::size_t>(std::lower_bound(vertex_map.begin(),vertex_map.end(),v)-vertex_map.begin());
                        for (std::size_t c=parent_offsets[e];c<parent_offsets[e+1];++c) variable_map.push_back(c);
                        child.edges.push_back(std::move(edge));
                    }
                    if (!graph_connected(child)) continue;
                    pending.push_back(selected);
                    const auto c=canonicalize_graph(child);
                    GraphMap map;
                    for (auto v : c.map.vertices) map.vertices.push_back(vertex_map[v]);
                    for (auto e : c.map.edges) map.edges.push_back(selected[e]);
                    for (auto v : c.map.channels) map.channels.push_back(variable_map[v]);
                    result[i].push_back({indices.at(c.key),std::move(map)});
                }
            }
        }
        plan_->graph_subclusters=std::move(result);
    });
    return plan_->graph_subclusters[index];
}
const std::vector<EmbeddedSubcluster>& WhiteGraphExpansion::embedding_subclusters(std::size_t index) const {
    (void)embeddings().at(index);
    std::call_once(plan_->embedding_subclusters_once,[&] {
        const auto& embeddings_=embeddings();
        const auto& lattice_=structure();
        std::vector<std::vector<EmbeddedSubcluster>> result(embeddings_.size());
        std::map<Cluster,std::size_t> physical_indices;
        for (std::size_t i=0;i<embeddings_.size();++i) physical_indices.emplace(embeddings_[i].edges,i);
        for (std::size_t i=0;i<embeddings_.size();++i) {
            const auto& parent=embeddings_[i];
            for (const auto& sub : graph_subclusters(parent.graph)) {
                Cluster edges;
                for (auto e : sub.map.edges) edges.push_back(parent.edges[parent.map.edges[e]]);
                const auto normalized=normalize(edges);
                const auto child_index=physical_indices.at(normalized.cluster);
                const auto& child=embeddings_[child_index];
                GraphMap map;
                for (auto site : child.sites) {
                    site.cell=detail::translate(site.cell,normalized.shift);
                    map.vertices.push_back(static_cast<std::size_t>(std::lower_bound(parent.sites.begin(),parent.sites.end(),site)-parent.sites.begin()));
                }
                std::vector<std::size_t> parent_offsets{0};
                for (const auto& e : parent.edges) parent_offsets.push_back(parent_offsets.back()+lattice_.interactions[e.type].channels.size());
                for (auto edge : child.edges) {
                    edge.origin=detail::translate(edge.origin,normalized.shift);
                    const auto e=static_cast<std::size_t>(std::lower_bound(parent.edges.begin(),parent.edges.end(),edge)-parent.edges.begin());
                    map.edges.push_back(e);
                    for (std::size_t c=parent_offsets[e];c<parent_offsets[e+1];++c) map.channels.push_back(c);
                }
                result[i].push_back({child_index,std::move(map)});
            }
        }
        plan_->embedding_subclusters=std::move(result);
    });
    return plan_->embedding_subclusters[index];
}
ScalarEvaluator::ScalarEvaluator(Function function) : function_(std::make_shared<const Function>(std::move(function))) {
    if (!*function_) throw std::invalid_argument("empty symbolic scalar evaluator");
}
std::shared_ptr<const SymbolicBlock> GraphCache::block(const CanonicalGraph& graph,double gap,
    const EffectiveOperator& effective,const std::vector<State>& basis,bool linked) {
    return compiled_block(graph,gap,effective,basis,linked,nullptr);
}
std::shared_ptr<const SymbolicBlock> GraphCache::compiled_block(const CanonicalGraph& graph,double gap,
    const EffectiveOperator& effective,const std::vector<State>& basis,bool linked,const ClusterModel* model) {
    std::string key; number(key,linked); string(key,evaluation_key(graph.graph)); real(key,gap); string(key,effective.identity()); sequence(key,basis);
    std::lock_guard lock(mutex_);
    if (const auto it=blocks_.find(key);it!=blocks_.end()) { ++hits_; return it->second; }
    std::vector<std::size_t> channel_edges;
    if (linked) for (std::size_t e=0;e<graph.graph.edges.size();++e)
        for (std::size_t c=0;c<graph.graph.edges[e].channels.size();++c) channel_edges.push_back(e);
    std::optional<ClusterModel> owned;
    if (!model) { owned.emplace(graph.graph.model(gap)); model=&*owned; }
    auto value=std::make_shared<const SymbolicBlock>(effective.symbolic_block(*model,basis,channel_edges));
    blocks_.emplace(std::move(key),value); ++evaluations_; return value;
}
std::shared_ptr<const SymbolicSeries> GraphCache::scalar(const CanonicalGraph& graph,double gap,
    unsigned order,const ScalarEvaluator& evaluator) {
    if (order>64) throw std::invalid_argument("scalar order exceeds 64");
    std::string key; string(key,evaluation_key(graph.graph)); real(key,gap); number(key,order);
    number(key,reinterpret_cast<std::uintptr_t>(evaluator.function().get()));
    std::lock_guard lock(mutex_);
    if (const auto it=scalars_.find(key);it!=scalars_.end()) { ++hits_; return it->second; }
    auto result=(*evaluator.function())(graph.graph,gap,order);
    if (result.size()!=order+1 || !result[0].empty()) throw std::invalid_argument("scalar callback needs order+1 corrections and empty constant");
    for (std::size_t n=1;n<result.size();++n) for (const auto& [m,c] : result[n]) {
        if (m.degree()!=n || !std::isfinite(c.real()) || !std::isfinite(c.imag())) throw std::invalid_argument("invalid scalar monomial degree or amplitude");
        m.for_each_power([&](std::size_t v,unsigned) { if (v>=graph.graph.variables()) throw std::invalid_argument("invalid scalar variable"); });
    }
    auto value=std::make_shared<const SymbolicSeries>(std::move(result));
    scalars_.emplace(std::move(key),value); evaluators_.push_back(evaluator.function()); ++evaluations_; return value;
}
void GraphCache::record_hit() { std::lock_guard lock(mutex_); ++hits_; }
std::size_t GraphCache::evaluations() const { std::lock_guard lock(mutex_); return evaluations_; }
std::size_t GraphCache::hits() const { std::lock_guard lock(mutex_); return hits_; }
std::vector<Matrix> WhiteGraphExpansion::block(std::size_t index,const EffectiveOperator& effective,
    const std::vector<State>& basis,bool linked) const {
    const auto& physical=structural_model(index);
    const auto& embedding=embeddings()[index];
    if (physical.fermionic()) physical.require_operator_grading();
    if (effective.order()>max_edges_) throw std::invalid_argument("graph expansion does not cover order");
    std::lock_guard lock(plan_->mutex);
    auto& program=plan_->readouts[effective.identity()];
    const Plan::BasisKey key{index,basis,linked};
    auto found=program.find(key);
    if (found==program.end()) {
        const auto& graph=graphs().at(embedding.graph);
        const auto& canonical=plan_->canonical_models[embedding.graph];
        std::map<State,std::pair<Eigen::Index,int>> indices;
        const auto inversions=detail::gather_inversions(physical.sites(),embedding.map.vertices);
        for (std::size_t i=0;i<basis.size();++i) {
            const auto local=physical.decode(basis[i]);
            std::vector<unsigned> reordered;
            for (auto v : embedding.map.vertices) reordered.push_back(local.at(v));
            const auto state=canonical.encode(reordered);
            const int sign=canonical.fermionic() ? detail::gather_sign(physical.spaces(),local,inversions) : 1;
            if (!indices.emplace(state,std::pair{static_cast<Eigen::Index>(i),sign}).second)
                throw std::invalid_argument("duplicate external state");
        }
        std::vector<State> requested;
        for (const auto& [state,position] : indices) { (void)position; requested.push_back(state); }
        Plan::Readout readout;
        readout.block=cache_->compiled_block(graph,structure().gap,effective,requested,linked,&canonical);
        readout.entries.resize(effective.order()+1);
        for (std::size_t n=0;n<readout.entries.size();++n) for (const auto& [states,p] : readout.block->coefficients[n]) {
            const auto [i,si]=indices.at(states.first); const auto [j,sj]=indices.at(states.second);
            readout.entries[n].push_back({i,j,si*sj,&p});
        }
        found=program.emplace(key,std::move(readout)).first;
    } else cache_->record_hit();
    auto result=detail::matrix_series(basis.size(),effective.order()+1);
    const auto& values=couplings_[index];
    for (std::size_t n=0;n<result.size();++n) for (const auto& e : found->second.entries[n])
        result[n](e.row,e.column)=static_cast<double>(e.sign)*substitute(*e.polynomial,values);
    return result;
}
}
