#include "incidence.hpp"
#include "graph_structure.hpp"
#include "nauty_support.hpp"
#include <algorithm>
#include <limits>
#include <numeric>
#include <stdexcept>

extern "C" void pcut_nauty_release_unowned();

namespace pcut::detail {
int Incidence::node(std::string color) {
    // nauty uses signed int counts and expressions including n+2.
    if (adjacent.size()>=static_cast<std::size_t>(std::numeric_limits<int>::max()-2))
        throw std::length_error("incidence vertex count overflow");
    const auto id=static_cast<int>(adjacent.size());
    adjacent.emplace_back(); colors.push_back(std::move(color)); return id;
}
void Incidence::join(int a,int b) {
    adjacent[static_cast<std::size_t>(a)].push_back(b);
    adjacent[static_cast<std::size_t>(b)].push_back(a);
}
Incidence encode_incidence(const WhiteGraph& graph,
                           const std::vector<std::string>& labels,
                           const std::vector<std::string>& structures) {
    Incidence result;
    for (const auto& label : labels) result.node("V"+label);
    for (std::size_t e=0;e<graph.edges.size();++e) {
        const auto& edge=graph.edges[e];
        // The full sorted channel multiset is already in this color. Channel
        // leaves add no information; occurrence maps are reconstructed later.
        const int occurrence=result.node("E"+structures[e]);
        for (std::size_t leg=0;leg<edge.legs.size();++leg) {
            std::string color="L"; number(color,leg);
            const int port=result.node(std::move(color));
            result.join(occurrence,port);
            result.join(port,static_cast<int>(edge.legs[leg]));
        }
    }
    return result;
}
Incidence encode_incidence(const WhiteGraph& graph) {
    std::vector<std::string> labels,structures;
    for (const auto& space : graph.spaces) labels.push_back(space_key(space));
    for (const auto& edge : graph.edges) structures.push_back(edge_structure(edge));
    return encode_incidence(graph,labels,structures);
}
namespace {
thread_local std::vector<int>* indices=nullptr;
void level(int*,int*,int,int*,statsblk*,int,int index,int,int cells,int,int n) {
    if (cells!=n) indices->push_back(index); // reserved n entries before entering C
}
}
NautyScratch::~NautyScratch() {
    nauty_freedyn(); nautil_freedyn(); nausparse_freedyn(); naugraph_freedyn();
    pcut_nauty_release_unowned(); indices=nullptr;
}
sparsegraph SparseStorage::view() {
    sparsegraph result{};
    result.nv=static_cast<int>(starts.size()); result.nde=neighbors.size();
    result.v=starts.data(); result.vlen=starts.size();
    result.d=degrees.data(); result.dlen=degrees.size();
    result.e=neighbors.data(); result.elen=neighbors.size();
    return result;
}
PreparedIncidence::PreparedIncidence(const Incidence& input) {
    if (input.adjacent.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()/1000)*WORDSIZE-WORDSIZE)
        throw std::length_error("nauty workspace index overflow");
    const auto n=static_cast<int>(input.adjacent.size());
    if (!n) throw std::invalid_argument("empty incidence graph");
    partition.resize(static_cast<std::size_t>(n),1);
    lab.resize(static_cast<std::size_t>(n));
    std::iota(lab.begin(),lab.end(),0);
    std::stable_sort(lab.begin(),lab.end(),[&](int a,int b) {
        return input.colors[static_cast<std::size_t>(a)]<input.colors[static_cast<std::size_t>(b)];
    });
    auto& [starts,degrees,neighbors]=graph;
    for (int i=0;i<n;++i) {
        const auto v=static_cast<std::size_t>(i);
        if (i==n-1 || input.colors[static_cast<std::size_t>(lab[v])]!=
            input.colors[static_cast<std::size_t>(lab[v+1])]) partition[v]=0;
        starts.push_back(neighbors.size());
        if (input.adjacent[v].size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("incidence degree overflow");
        degrees.push_back(static_cast<int>(input.adjacent[v].size()));
        neighbors.insert(neighbors.end(),input.adjacent[v].begin(),input.adjacent[v].end());
        std::sort(neighbors.begin()+static_cast<std::ptrdiff_t>(starts.back()),neighbors.end());
    }
}
IncidenceLabel label_incidence(const Incidence& input) {
    PreparedIncidence prepared(input);
    auto canonical_storage=prepared.graph;
    auto graph=prepared.graph.view(), canonical=canonical_storage.view();
    std::vector<int> orbits(prepared.lab.size());
    IncidenceLabel result;
    result.canonical_to_input=std::move(prepared.lab);
    result.group_indices.reserve(result.canonical_to_input.size());
    DEFAULTOPTIONS_SPARSEGRAPH(options);
    options.getcanon=TRUE; options.defaultptn=FALSE; options.userlevelproc=level;
    options.schreier=FALSE; // deterministic exact stabilizer-index path
    statsblk stats{};
    NautyScratch scratch;
    indices=&result.group_indices;
    sparsenauty(&graph,result.canonical_to_input.data(),prepared.partition.data(),orbits.data(),&options,&stats,&canonical);
    indices=nullptr;
    if (stats.errstatus) throw std::runtime_error("nauty canonical labeling failed");
    return result;
}
} // namespace pcut::detail
