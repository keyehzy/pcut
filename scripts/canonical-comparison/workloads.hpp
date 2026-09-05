#pragma once
#include <pcut/pcut.hpp>
#include <algorithm>
#include <string>
using namespace pcut;
inline PeriodicLattice lattice_for(const std::string& name) {
    if (name=="dimer") return models::dimerized_chain(0.17);
    if (name=="hubbard") return models::hubbard_square();
    auto lattice=models::ising_chain(); lattice.dimension=2;
    const auto channels=lattice.interactions[0].channels;
    lattice.interactions.clear();
    const std::vector<Coordinate> displacements=name=="four_color" ?
        std::vector<Coordinate>{{1,0},{0,1},{1,1},{1,-1}} : std::vector<Coordinate>{{1,0},{0,1}};
    for (auto d : displacements) lattice.interactions.push_back({{{{0,0},0},{d,0}},channels});
    return lattice;
}
inline std::vector<WhiteGraph> workload(const std::string& name) {
    if (name=="square" || name=="four_color" || name=="hubbard" || name=="dimer") {
        const auto lattice=lattice_for(name);
        const WhiteGraphExpansion catalog(lattice,name=="dimer" ? 6 : name=="four_color" ? 3 : 4);
        std::vector<WhiteGraph> graphs;
        // Recover the same sorted physical inputs for both engines, regardless
        // of their different canonical ordering. Geometry/couplings are omitted.
        for (const auto& embedding : catalog.embeddings()) {
            const auto sites=vertices(lattice,embedding.edges);
            WhiteGraph g;
            for (const auto& site : sites) g.spaces.push_back(lattice.cell[site.basis]);
            for (const auto& edge : embedding.edges) {
                const auto& t=lattice.interactions[edge.type]; GraphEdge e;
                for (auto leg : t.legs) {
                    for (std::size_t d=0;d<leg.cell.size();++d) leg.cell[d]+=edge.origin[d];
                    e.legs.push_back(static_cast<std::size_t>(std::lower_bound(sites.begin(),sites.end(),leg)-sites.begin()));
                }
                for (const auto& c : t.channels) e.channels.push_back(c.op);
                g.edges.push_back(std::move(e));
            }
            graphs.push_back(std::move(g));
        }
        return graphs;
    }
    WhiteGraph g;
    const std::size_t n=name=="star20" ? 21 : name=="cycle64" ? 64 : name=="star8" ? 9 : name=="complete7" ? 7 : name=="duplicates" ? 2 : 8;
    g.spaces.resize(n,models::ising_chain().cell[0]);
    const auto add=[&](std::size_t a,std::size_t b) { g.edges.push_back({{a,b},{{Matrix::Identity(4,4),false}}}); };
    if (name.starts_with("star")) for (std::size_t v=1;v<n;++v) add(0,v);
    else if (name=="complete7") { for (std::size_t a=0;a<n;++a) for (std::size_t b=0;b<n;++b) if (a!=b) add(a,b); }
    else if (name=="duplicates") { for (unsigned e=0;e<32;++e) add(0,1); for (auto& e : g.edges) e.channels.resize(8,e.channels[0]); }
    else for (std::size_t v=0;v<n;++v) add(v,(v+1)%n);
    return {g};
}
