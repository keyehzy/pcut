#include <pcut/pcut.hpp>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sys/resource.h>
#include "workloads.hpp"
using namespace pcut;
int main(int argc,char** argv) {
    if (argc!=3) return 2;
    auto graphs=workload(argv[1]);
    const auto start=std::chrono::steady_clock::now();
    std::size_t bytes=0,automorphisms=0;
    const auto repeats=std::stoul(argv[2]);
    for (std::size_t r=0;r<repeats;++r) for (const auto& g : graphs) {
        const auto c=canonicalize(g); bytes+=c.key.size(); automorphisms+=c.vertex_automorphisms;
    }
    const double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    rusage usage{}; getrusage(RUSAGE_SELF,&usage);
    std::cout<<std::setprecision(17)<<seconds<<' '<<usage.ru_maxrss<<' '<<graphs.size()*repeats<<' '<<bytes<<' '<<automorphisms<<'\n';
}
