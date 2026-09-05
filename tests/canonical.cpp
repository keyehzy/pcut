#include <pcut/pcut.hpp>
#include <catch2/catch_test_macros.hpp>
#include "canonical_oracle.hpp"
#include <algorithm>
#include <future>
#include <numeric>
#include <random>
using namespace pcut;
namespace {
WhiteGraph make_graph(std::size_t n,const std::vector<std::vector<std::size_t>>& legs) {
    WhiteGraph g; g.spaces.resize(n,models::ising_chain().cell[0]);
    for (const auto& l : legs) g.edges.push_back({l,{{Matrix::Identity(1u<<l.size(),1u<<l.size()),false}}});
    return g;
}
void check_map(const WhiteGraph& input,const CanonicalGraph& c) {
    auto check_permutation=[](auto actual,std::size_t n) {
        std::sort(actual.begin(),actual.end()); std::vector<std::size_t> expected(n);
        std::iota(expected.begin(),expected.end(),0); REQUIRE(actual==expected);
    };
    check_permutation(c.map.vertices,input.spaces.size());
    check_permutation(c.map.edges,input.edges.size());
    check_permutation(c.map.channels,input.variables());
    std::vector<std::size_t> offsets{0};
    for (const auto& e : input.edges) offsets.push_back(offsets.back()+e.channels.size());
    for (std::size_t v=0;v<c.graph.spaces.size();++v)
        REQUIRE(test_oracle::space_key(c.graph.spaces[v])==test_oracle::space_key(input.spaces[c.map.vertices[v]]));
    std::size_t flat=0;
    for (std::size_t e=0;e<c.graph.edges.size();++e) {
        const auto source=c.map.edges[e]; const auto& edge=c.graph.edges[e];
        REQUIRE(edge.legs.size()==input.edges[source].legs.size());
        for (std::size_t l=0;l<edge.legs.size();++l)
            REQUIRE(c.map.vertices[edge.legs[l]]==input.edges[source].legs[l]);
        for (const auto& channel : edge.channels) {
            const auto variable=c.map.channels[flat++];
            REQUIRE(variable>=offsets[source]); REQUIRE(variable<offsets[source+1]);
            REQUIRE(test_oracle::channel_key(channel)==test_oracle::channel_key(input.edges[source].channels[variable-offsets[source]]));
        }
    }
}
WhiteGraph shuffled(WhiteGraph g,std::mt19937& rng) {
    std::vector<std::size_t> p(g.spaces.size()); std::iota(p.begin(),p.end(),0); std::shuffle(p.begin(),p.end(),rng);
    auto old=g.spaces; for (std::size_t v=0;v<p.size();++v) g.spaces[p[v]]=old[v];
    std::shuffle(g.edges.begin(),g.edges.end(),rng);
    for (auto& e : g.edges) { for (auto& v : e.legs) v=p[v]; std::shuffle(e.channels.begin(),e.channels.end(),rng); }
    return g;
}
void compare(const WhiteGraph& g,std::mt19937& rng) {
    const auto oracle=test_oracle::canonicalize(g), current=canonicalize(g);
    REQUIRE(test_oracle::canonicalize(current.graph).key==oracle.key);
    REQUIRE(current.vertex_automorphisms==oracle.vertex_automorphisms);
    check_map(g,current);
    for (int repeat=0;repeat<6;++repeat) {
        const auto input=shuffled(g,rng); const auto candidate=canonicalize(input);
        REQUIRE(candidate.key==current.key);
        REQUIRE(candidate.vertex_automorphisms==oracle.vertex_automorphisms);
        check_map(input,candidate);
    }
}
}
TEST_CASE("Incidence labels agree with exhaustive equivalence and reconstruct every occurrence", "[white][canonical]") {
    std::mt19937 rng(72819);
    // All connected directed simple graphs on three vertices: includes pairs
    // that refinement alone cannot serve as an identity for in larger catalogs.
    std::map<std::string,std::string> forward,backward;
    for (unsigned bits=1;bits<64;++bits) {
        std::vector<std::vector<std::size_t>> legs;
        unsigned bit=0;
        for (std::size_t a=0;a<3;++a) for (std::size_t b=0;b<3;++b) if (a!=b) {
            if (bits&(1u<<bit)) legs.push_back({a,b}); ++bit;
        }
        auto g=make_graph(3,legs); if (!test_oracle::graph_connected(g)) continue;
        compare(g,rng);
        const auto old=test_oracle::canonicalize(g).key, now=canonicalize(g).key;
        if (forward.contains(old)) REQUIRE(forward.at(old)==now); else forward[old]=now;
        if (backward.contains(now)) REQUIRE(backward.at(now)==old); else backward[now]=old;
    }
    for (unsigned sample=0;sample<100;++sample) {
        const std::size_t n=2+rng()%5; std::vector<std::vector<std::size_t>> legs;
        for (std::size_t v=1;v<n;++v) legs.push_back({rng()%v,v});
        legs.push_back({rng()%n});
        if (n>=3) legs.push_back({2,0,1});
        auto g=make_graph(n,legs);
        if (sample%2) g.spaces[rng()%n].name="species B";
        if (sample%3==0) g.edges.push_back(g.edges[rng()%g.edges.size()]);
        for (auto& e : g.edges) {
            e.channels.push_back(e.channels[0]);
            auto different=e.channels[0]; different.matrix(0,0)=2; e.channels.push_back(different);
        }
        compare(g,rng);
    }
    // Connected regular graphs with identical equitable colors but different
    // structure: triangular prism and K3,3 (both directions on every bond).
    auto regular=[](std::vector<std::vector<std::size_t>> bonds) {
        auto directed=bonds; for (auto e : bonds) { std::reverse(e.begin(),e.end()); directed.push_back(e); }
        return make_graph(6,directed);
    };
    const auto prism=regular({{0,1},{1,2},{2,0},{3,4},{4,5},{5,3},{0,3},{1,4},{2,5}});
    const auto bipartite=regular({{0,3},{0,4},{0,5},{1,3},{1,4},{1,5},{2,3},{2,4},{2,5}});
    compare(prism,rng); compare(bipartite,rng);
    REQUIRE(canonicalize(prism).key!=canonicalize(bipartite).key);
    auto unequal=make_graph(3,{{0,1,2},{1}});
    unequal.spaces[1].charges={0,1,2};
    unequal.edges[0].channels[0].matrix=Matrix::Identity(12,12);
    unequal.edges[1].channels[0].matrix=Matrix::Identity(3,3);
    compare(unequal,rng);
    for (auto lattice : {models::ising_chain(),models::dimerized_chain(0.17),models::hubbard_square()}) {
        const WhiteGraphExpansion catalog(lattice,4);
        for (const auto& c : catalog.graphs()) compare(c.graph,rng);
    }
}
TEST_CASE("Incidence colors preserve exact physical attributes and ordered hyperedge roles", "[white][canonical]") {
    auto g=make_graph(3,{{0,1,2},{0},{1,2}});
    for (auto& s : g.spaces) { s.parity={0,1}; s.particles={0,1}; }
    const auto original=canonicalize(g).key;
    auto changed=g; changed.spaces[0].name+="x"; REQUIRE(canonicalize(changed).key!=original);
    changed=g; changed.spaces[0].charges[1]=2; REQUIRE(canonicalize(changed).key!=original);
    changed=g; changed.spaces[0].vacuum_energy=0.25; REQUIRE(canonicalize(changed).key!=original);
    changed=g; changed.spaces[0].particles[1]=3; REQUIRE(canonicalize(changed).key!=original);
    changed=g; for (auto& s : changed.spaces) s.particles.clear();
    const auto parity_original=canonicalize(changed).key;
    changed.spaces[0].parity[1]=0; REQUIRE(canonicalize(changed).key!=parity_original);
    changed=g; changed.edges[0].channels[0].fermionic=true; REQUIRE(canonicalize(changed).key!=original);
    changed=g; changed.edges[0].channels[0].matrix(0,0)=std::nextafter(1.0,2.0); REQUIRE(canonicalize(changed).key!=original);
    changed=g; changed.edges[0].channels[0].matrix(0,1)={0,0.2}; changed.edges[0].channels[0].matrix(1,0)={0,-0.2}; REQUIRE(canonicalize(changed).key!=original);
    changed=g; std::swap(changed.edges[0].legs[0],changed.edges[0].legs[1]); REQUIRE(canonicalize(changed).key!=original);
    changed=g; changed.edges.push_back(g.edges[0]); REQUIRE(canonicalize(changed).key!=original);
    changed=g; changed.edges[0].channels.push_back(g.edges[0].channels[0]); REQUIRE(canonicalize(changed).key!=original);
    // Equal local dimensions alone are insufficient to identify species.
    changed=g; changed.spaces[1].name="other"; std::mt19937 rng(1); compare(changed,rng);
}
TEST_CASE("Vertex group orders remove huge gadget kernels before exact overflow checking", "[white][canonical]") {
    auto parallel=make_graph(2,{{0,1}});
    parallel.edges.resize(64,parallel.edges[0]);
    for (auto& e : parallel.edges) e.channels.resize(24,e.channels[0]);
    const auto duplicate=canonicalize(parallel);
    REQUIRE(duplicate.vertex_automorphisms==1); // 64! edge kernel
    REQUIRE(duplicate.graph.variables()==1536);
    check_map(parallel,duplicate);
    for (std::size_t leaves : {4u,12u,20u,21u}) {
        std::vector<std::vector<std::size_t>> legs;
        for (std::size_t v=1;v<=leaves;++v) legs.push_back({0,v});
        auto star=make_graph(leaves+1,legs);
        if (leaves==21 || (sizeof(std::size_t)<8 && leaves>=13)) {
            REQUIRE_THROWS_AS(canonicalize(star),std::length_error);
        } else {
            std::size_t expected=1; for (std::size_t i=2;i<=leaves;++i) expected*=i;
            REQUIRE(canonicalize(star).vertex_automorphisms==expected);
        }
    }
    auto graph=make_graph(5,{{0,1},{0,2},{0,3},{0,4}});
    const auto expected=canonicalize(graph).key;
    std::vector<std::future<void>> jobs;
    for (unsigned worker=0;worker<4;++worker) jobs.push_back(std::async(std::launch::async,[&,worker] {
        std::mt19937 rng(worker); for (int i=0;i<30;++i)
            if (canonicalize(shuffled(graph,rng)).key!=expected) throw std::runtime_error("concurrent labeling mismatch");
    }));
    for (auto& job : jobs) REQUIRE_NOTHROW(job.get());
}
