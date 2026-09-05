#include "workloads.hpp"
#include "nauty_support.hpp"
#include <traces.h>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sys/resource.h>
namespace {
struct BackendScratch {
    detail::NautyScratch nauty;
    bool traces;
    ~BackendScratch() { if (traces) { traces_freedyn(); schreier_freedyn(); } }
};
struct Prepared {
    detail::PreparedIncidence input;
    explicit Prepared(const WhiteGraph& g) : input(detail::encode_incidence(g)) {}
    std::vector<int> solve(bool traces) {
        auto labels=input.lab,ptn=input.partition,orbits=input.lab;
        auto canonical=input.graph;
        auto g=input.graph.view(),c=canonical.view();
        BackendScratch scratch{{},traces};
        if (traces) {
            DEFAULTOPTIONS_TRACES(options); options.getcanon=TRUE; options.defaultptn=FALSE;
            TracesStats stats{}; Traces(&g,labels.data(),ptn.data(),orbits.data(),&options,&stats,&c);
            if (stats.errstatus) throw std::runtime_error("Traces error");
        } else {
            DEFAULTOPTIONS_SPARSEGRAPH(options); options.getcanon=TRUE; options.defaultptn=FALSE; options.schreier=FALSE;
            statsblk stats{}; sparsenauty(&g,labels.data(),ptn.data(),orbits.data(),&options,&stats,&c);
            if (stats.errstatus) throw std::runtime_error("nauty error");
        }
        // Serialize adjacency rows without sorting away the canonical vertex order.
        auto& [cv,cd,ce]=canonical;
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
