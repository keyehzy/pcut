#include "incidence.hpp"
#include <nauty.h>
#include <nausparse.h>
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
namespace {
thread_local std::vector<int>* indices=nullptr;
void level(int*,int*,int,int*,statsblk*,int,int index,int,int cells,int,int n) {
    if (cells!=n) indices->push_back(index); // reserved n entries before entering C
}
struct Scratch {
    ~Scratch() { nauty_freedyn(); nautil_freedyn(); nausparse_freedyn(); naugraph_freedyn(); pcut_nauty_release_unowned(); indices=nullptr; }
};
}
IncidenceLabel label_incidence(const Incidence& input) {
    const auto n=static_cast<int>(input.adjacent.size());
    if (input.adjacent.size()>static_cast<std::size_t>(std::numeric_limits<int>::max()/1000)*WORDSIZE-WORDSIZE)
        throw std::length_error("nauty workspace index overflow");
    if (!n) throw std::invalid_argument("empty incidence graph");
    std::vector<std::size_t> starts;
    std::vector<int> degrees, neighbors, partition(static_cast<std::size_t>(n),1), orbits(static_cast<std::size_t>(n));
    IncidenceLabel result;
    result.canonical_to_input.resize(static_cast<std::size_t>(n));
    result.group_indices.reserve(static_cast<std::size_t>(n));
    std::iota(result.canonical_to_input.begin(),result.canonical_to_input.end(),0);
    std::stable_sort(result.canonical_to_input.begin(),result.canonical_to_input.end(),[&](int a,int b) {
        return input.colors[static_cast<std::size_t>(a)]<input.colors[static_cast<std::size_t>(b)];
    });
    for (int i=0;i<n;++i) {
        const auto v=static_cast<std::size_t>(i);
        if (i==n-1 || input.colors[static_cast<std::size_t>(result.canonical_to_input[v])]!=
            input.colors[static_cast<std::size_t>(result.canonical_to_input[v+1])]) partition[v]=0;
        starts.push_back(neighbors.size());
        if (input.adjacent[v].size()>static_cast<std::size_t>(std::numeric_limits<int>::max()))
            throw std::length_error("incidence degree overflow");
        degrees.push_back(static_cast<int>(input.adjacent[v].size()));
        neighbors.insert(neighbors.end(),input.adjacent[v].begin(),input.adjacent[v].end());
        std::sort(neighbors.begin()+static_cast<std::ptrdiff_t>(starts.back()),neighbors.end());
    }
    sparsegraph graph{};
    graph.nv=n; graph.nde=neighbors.size(); graph.v=starts.data(); graph.d=degrees.data(); graph.e=neighbors.data();
    // Supply the canonical graph storage ourselves; upstream never owns it.
    auto canonical_starts=starts;
    auto canonical_degrees=degrees, canonical_neighbors=neighbors;
    sparsegraph canonical{};
    canonical.nv=n; canonical.nde=neighbors.size();
    canonical.v=canonical_starts.data(); canonical.vlen=canonical_starts.size();
    canonical.d=canonical_degrees.data(); canonical.dlen=canonical_degrees.size();
    canonical.e=canonical_neighbors.data(); canonical.elen=canonical_neighbors.size();
    DEFAULTOPTIONS_SPARSEGRAPH(options);
    options.getcanon=TRUE; options.defaultptn=FALSE; options.userlevelproc=level;
    options.schreier=FALSE; // deterministic exact stabilizer-index path
    statsblk stats{};
    Scratch scratch;
    indices=&result.group_indices;
    sparsenauty(&graph,result.canonical_to_input.data(),partition.data(),orbits.data(),&options,&stats,&canonical);
    indices=nullptr;
    if (stats.errstatus) throw std::runtime_error("nauty canonical labeling failed");
    return result;
}
} // namespace pcut::detail
