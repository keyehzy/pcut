#include "workloads.hpp"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <sys/resource.h>
int main(int argc,char** argv) {
    if (argc!=3) return 2;
    const std::string name=argv[1];
    const auto lattice=lattice_for(name);
    const auto start=std::chrono::steady_clock::now();
    const WhiteGraphExpansion catalog(lattice,name=="dimer" ? 6 : name=="four_color" ? 3 : 4);
    const auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
    rusage usage{}; getrusage(RUSAGE_SELF,&usage);
    std::cout<<std::setprecision(17)<<seconds<<' '<<usage.ru_maxrss<<" 1 "<<catalog.graphs().size()<<' '<<catalog.embeddings().size()<<'\n';
}
