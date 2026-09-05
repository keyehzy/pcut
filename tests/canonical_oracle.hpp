// Exhaustive starting-revision algorithm, test-only; use on small graphs.
#pragma once
#include <pcut/white_graph.hpp>

#include <algorithm>
#include <bit>
#include <functional>
#include <limits>
#include <numeric>
#include <set>
#include <stdexcept>

namespace pcut {
namespace test_oracle {
namespace {
void number(std::string& key, std::uint64_t n) {
    char bytes[8];
    for (unsigned i=0;i<8;++i) bytes[i]=static_cast<char>((n>>(8*i))&255);
    key.append(bytes,8);
}
void real(std::string& key,double x) { number(key,std::bit_cast<std::uint64_t>(x==0 ? 0.0 : x)); }
void string(std::string& key,const std::string& s) { number(key,s.size()); key+=s; }
template<class T> void sequence(std::string& key,const std::vector<T>& s) {
    number(key,s.size()); for (auto x : s) number(key,static_cast<std::uint64_t>(x));
}
std::string space_key(const LocalSpace& space) {
    std::string key;
    string(key,space.name); sequence(key,space.charges); sequence(key,space.particles);
    sequence(key,space.parity); real(key,space.vacuum_energy); return key;
}
std::string channel_key(const OperatorChannel& op) {
    std::string key;
    number(key,op.fermionic); number(key,static_cast<std::uint64_t>(op.matrix.rows()));
    number(key,static_cast<std::uint64_t>(op.matrix.cols()));
    for (Eigen::Index i=0;i<op.matrix.size();++i) { real(key,op.matrix.data()[i].real()); real(key,op.matrix.data()[i].imag()); }
    return key;
}
std::vector<std::size_t> channel_order(const GraphEdge& edge) {
    std::vector<std::size_t> order(edge.channels.size()); std::iota(order.begin(),order.end(),0);
    std::stable_sort(order.begin(),order.end(),[&](auto a,auto b) { return channel_key(edge.channels[a])<channel_key(edge.channels[b]); });
    return order;
}
std::string edge_structure(const GraphEdge& edge) {
    std::string key; number(key,edge.legs.size()); number(key,edge.channels.size());
    for (auto c : channel_order(edge)) string(key,channel_key(edge.channels[c]));
    return key;
}
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
inline CanonicalGraph canonicalize(const WhiteGraph& g) {
    if (g.edges.empty() || g.edges.size()>64 || g.spaces.empty() || g.spaces.size()>1000)
        throw std::invalid_argument("invalid abstract graph size");
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
    if (!graph_connected(g)) throw std::invalid_argument("abstract graph must be connected without isolated vertices");
    std::vector<std::string> space_keys;
    for (const auto& space : g.spaces) space_keys.push_back(space_key(space));
    auto labels=space_keys;
    std::vector<std::string> structures;
    std::vector<std::vector<std::size_t>> orders;
    for (const auto& edge : g.edges) {
        structures.push_back(edge_structure(edge)); orders.push_back(channel_order(edge));
    }
    auto unique_structures=structures;
    std::sort(unique_structures.begin(),unique_structures.end());
    unique_structures.erase(std::unique(unique_structures.begin(),unique_structures.end()),unique_structures.end());
    std::vector<std::size_t> structure_ids;
    for (const auto& key : structures) structure_ids.push_back(static_cast<std::size_t>(std::lower_bound(unique_structures.begin(),unique_structures.end(),key)-unique_structures.begin()));
    // Ordered incidence color refinement. Full permutation search within the
    // resulting cells makes this exact even when refinement cannot distinguish graphs.
    for (std::size_t pass=0;pass<g.spaces.size();++pass) {
        std::map<std::string,std::size_t> ids;
        for (const auto& l : labels) ids.emplace(l,0);
        std::size_t next=0; for (auto& [l,id] : ids) { (void)l; id=next++; }
        std::vector<std::string> refined(g.spaces.size());
        for (std::size_t v=0;v<g.spaces.size();++v) {
            number(refined[v],ids.at(labels[v]));
            std::vector<std::string> incident;
            for (std::size_t ei=0;ei<g.edges.size();++ei) for (std::size_t leg=0;leg<g.edges[ei].legs.size();++leg) if (g.edges[ei].legs[leg]==v) {
                const auto& e=g.edges[ei];
                std::string signature; number(signature,structure_ids[ei]); number(signature,leg);
                for (auto w : e.legs) number(signature,ids.at(labels[w]));
                incident.push_back(std::move(signature));
            }
            std::sort(incident.begin(),incident.end());
            for (const auto& signature : incident) string(refined[v],signature);
        }
        labels=std::move(refined);
    }
    std::map<std::string,std::vector<std::size_t>> cells;
    for (std::size_t v=0;v<labels.size();++v) cells[labels[v]].push_back(v);
    std::vector<std::vector<std::size_t>> groups;
    for (auto& [key,cell] : cells) { (void)key; groups.push_back(std::move(cell)); }
    CanonicalGraph best;
    const auto source_offsets=offsets(g);
    std::vector<std::size_t> permutation;
    std::function<void(std::size_t)> search=[&](std::size_t group) {
        if (group<groups.size()) {
            auto cell=groups[group];
            do {
                permutation.insert(permutation.end(),cell.begin(),cell.end()); search(group+1);
                permutation.resize(permutation.size()-cell.size());
            } while (std::next_permutation(cell.begin(),cell.end()));
            return;
        }
        std::vector<std::size_t> inverse(g.spaces.size());
        CanonicalGraph candidate;
        candidate.map.vertices=permutation;
        number(candidate.key,g.spaces.size());
        for (std::size_t v=0;v<permutation.size();++v) {
            inverse[permutation[v]]=v; candidate.graph.spaces.push_back(g.spaces[permutation[v]]);
            string(candidate.key,space_keys[permutation[v]]);
        }
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
        if (best.key.empty() || candidate.key<best.key) { best=std::move(candidate); best.vertex_automorphisms=1; }
        else if (candidate.key==best.key) {
            if (best.vertex_automorphisms==std::numeric_limits<std::size_t>::max()) throw std::length_error("automorphism count overflow");
            ++best.vertex_automorphisms;
        }
    };
    search(0); return best;
}

} // namespace test_oracle
} // namespace pcut
