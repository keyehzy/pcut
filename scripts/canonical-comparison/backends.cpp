#include "workloads.hpp"
#include "../../tests/canonical_oracle.hpp"
#include "incidence.hpp"
#include <nauty.h>
#include <traces.h>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sys/resource.h>
extern "C" void pcut_nauty_release_unowned();
namespace {
struct Prepared {
    std::vector<std::size_t> starts;
    std::vector<int> degrees,neighbors,lab,partition;
    Prepared(const WhiteGraph& g) {
        detail::Incidence graph;
        for (const auto& s : g.spaces) graph.node("V"+test_oracle::space_key(s));
        for (const auto& e : g.edges) {
            const auto occurrence=graph.node("E"+test_oracle::edge_structure(e));
            for (std::size_t l=0;l<e.legs.size();++l) {
                std::string color="L"; test_oracle::number(color,l);
                const auto port=graph.node(color); graph.join(port,occurrence); graph.join(port,static_cast<int>(e.legs[l]));
            }
            for (std::size_t c=0;c<e.channels.size();++c) {
                std::string color="C"; test_oracle::number(color,c);
                graph.join(occurrence,graph.node(std::move(color)));
            }
        }
        lab.resize(graph.colors.size()); std::iota(lab.begin(),lab.end(),0);
        std::stable_sort(lab.begin(),lab.end(),[&](int a,int b) { return graph.colors[a]<graph.colors[b]; });
        partition.resize(lab.size(),1);
        for (std::size_t i=0;i<lab.size();++i) {
            if (i+1==lab.size() || graph.colors[lab[i]]!=graph.colors[lab[i+1]]) partition[i]=0;
            starts.push_back(neighbors.size()); degrees.push_back(static_cast<int>(graph.adjacent[i].size()));
            auto adjacent=graph.adjacent[i]; std::sort(adjacent.begin(),adjacent.end());
            neighbors.insert(neighbors.end(),adjacent.begin(),adjacent.end());
        }
    }
    std::vector<int> solve(bool traces) {
        auto labels=lab,ptn=partition,orbits=lab,cd=degrees,ce=neighbors;
        auto cv=starts;
        sparsegraph g{},c{};
        g.nv=static_cast<int>(lab.size()); g.nde=neighbors.size(); g.v=starts.data(); g.d=degrees.data(); g.e=neighbors.data();
        c.nv=g.nv; c.nde=g.nde; c.v=cv.data(); c.vlen=cv.size(); c.d=cd.data(); c.dlen=cd.size(); c.e=ce.data(); c.elen=ce.size();
        if (traces) {
            DEFAULTOPTIONS_TRACES(options); options.getcanon=TRUE; options.defaultptn=FALSE;
            TracesStats stats{}; Traces(&g,labels.data(),ptn.data(),orbits.data(),&options,&stats,&c);
            if (stats.errstatus) throw std::runtime_error("Traces error");
            traces_freedyn(); schreier_freedyn();
        } else {
            DEFAULTOPTIONS_SPARSEGRAPH(options); options.getcanon=TRUE; options.defaultptn=FALSE; options.schreier=FALSE;
            statsblk stats{}; sparsenauty(&g,labels.data(),ptn.data(),orbits.data(),&options,&stats,&c);
            if (stats.errstatus) throw std::runtime_error("nauty error");
        }
        nauty_freedyn(); nautil_freedyn(); nausparse_freedyn(); naugraph_freedyn(); pcut_nauty_release_unowned();
        // Serialize adjacency rows without sorting away the canonical vertex order.
        std::vector<int> key;
        for (std::size_t i=0;i<cd.size();++i) {
            key.push_back(cd[i]); std::sort(ce.begin()+cv[i],ce.begin()+cv[i]+cd[i]);
            key.insert(key.end(),ce.begin()+cv[i],ce.begin()+cv[i]+cd[i]);
        }
        return key;
    }
};
}
int main(int argc,char** argv) {
    if (argc!=4) return 2;
    const bool traces=std::string(argv[1])=="traces";
    std::vector<Prepared> graphs;
    for (auto g : workload(argv[2])) {
        graphs.emplace_back(g);
        const auto expected=graphs.back().solve(traces);
        std::reverse(g.spaces.begin(),g.spaces.end()); std::reverse(g.edges.begin(),g.edges.end());
        for (auto& e : g.edges) { for (auto& v : e.legs) v=g.spaces.size()-1-v; std::reverse(e.channels.begin(),e.channels.end()); }
        Prepared permuted(g);
        if (permuted.solve(traces)!=expected) throw std::runtime_error("backend permutation mismatch");
    }
    const auto repeats=std::stoul(argv[3]); std::size_t checksum=0;
    const auto start=std::chrono::steady_clock::now();
    for (std::size_t r=0;r<repeats;++r) for (auto& g : graphs) checksum+=g.solve(traces).size();
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    rusage usage{}; getrusage(RUSAGE_SELF,&usage);
    std::cout<<std::setprecision(17)<<seconds<<' '<<usage.ru_maxrss<<' '<<graphs.size()*repeats<<' '<<checksum<<'\n';
}
