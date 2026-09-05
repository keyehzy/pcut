#include <pcut/pcut.hpp>
#include <chrono>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <sys/resource.h>

using namespace pcut;
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point start) { return std::chrono::duration<double>(Clock::now()-start).count(); }
PeriodicLattice model(const std::string& name,unsigned sweep) {
    if (name=="dimer") return models::dimerized_chain(0.17+0.001*sweep);
    if (name=="hubbard") {
        auto lattice=models::hubbard_square();
#ifdef PREVIOUS
        lattice.interactions[0].matrix*=0.7+0.001*sweep;
        lattice.interactions[1].matrix*=-0.2+0.001*sweep;
#else
        lattice.interactions[0].channels[0].coupling=0.7+0.001*sweep;
        lattice.interactions[1].channels[0].coupling=-0.2+0.001*sweep;
#endif
        return lattice;
    }
    auto lattice=models::ising_chain(); lattice.dimension=2;
#ifdef PREVIOUS
    const auto matrix=lattice.interactions[0].matrix;
#else
    const auto channels=lattice.interactions[0].channels;
#endif
    lattice.interactions.clear();
    const std::vector<Coordinate> displacements=name=="four_color" ?
        std::vector<Coordinate>{{1,0},{0,1},{1,1},{1,-1}} : std::vector<Coordinate>{{1,0},{0,1}};
    for (std::size_t i=0;i<displacements.size();++i) {
        const double ratio=0.7-0.3*static_cast<double>(i)+0.001*sweep;
#ifdef PREVIOUS
        lattice.interactions.push_back({{{{0,0},0},{displacements[i],0}},ratio*matrix});
#else
        auto bound=channels; bound[0].coupling=ratio;
        lattice.interactions.push_back({{{{0,0},0},{displacements[i],0}},bound});
#endif
    }
    return lattice;
}
int main(int argc,char** argv) {
    if (argc!=5) return 2;
    const std::string name=argv[1]; const auto order=static_cast<unsigned>(std::stoul(argv[2]));
    const auto count=static_cast<unsigned>(std::stoul(argv[3])); const bool particles=std::stoi(argv[4])!=0;
    const auto all=Clock::now();
    const auto lattice=model(name,0);
    const auto coefficient_start=Clock::now();
    const Coefficients coefficients(charge_changes(lattice),order);
    const EffectiveOperator effective(coefficients);
    const double coefficient_time=seconds(coefficient_start);
    const auto setup_start=Clock::now();
#ifdef PREVIOUS
    const auto topology=std::make_shared<ClusterTopology>(lattice,order);
#else
    const auto cache=std::make_shared<GraphCache>();
    const WhiteGraphExpansion topology(lattice,order,cache);
#endif
    const double setup_time=seconds(setup_start);
    std::cout << std::setprecision(17) << "setup " << coefficient_time << ' ' << setup_time << '\n';
    for (unsigned i=0;i<count;++i) {
        const auto bind_start=Clock::now();
#ifdef PREVIOUS
        const ClusterCatalog expansion(model(name,i),topology);
        const auto embeddings=expansion.entries().size(), graphs=embeddings;
        const std::size_t new_evaluations=embeddings;
#else
        const auto expansion=topology.bind(model(name,i).couplings());
        const auto embeddings=expansion.embeddings().size(), graphs=expansion.graphs().size();
        const auto before=cache->evaluations();
#endif
        const double bind_time=seconds(bind_start);
        const auto link_start=Clock::now();
        std::vector<double> output;
        if (name=="hubbard") {
            const auto result=linked_zero_charge(expansion,effective);
            output.resize(3*(order+1));
            for (const auto& weight : result.weights) for (unsigned n=0;n<=order;++n) {
                output[3*n]+=weight.block.coefficients[n].squaredNorm();
                output[3*n+1]+=weight.block.coefficients[n].trace().real();
                output[3*n+2]+=weight.block.coefficients[n].trace().imag();
            }
        } else {
            const auto result=linked_expand(expansion,effective,{particles});
            for (const auto x : result.energy_per_cell) { output.push_back(x.real()); output.push_back(x.imag()); }
            if (particles) {
                const auto series=result.bloch_series(std::vector<double>(lattice.dimension,0.37));
                for (const auto& m : series) for (Eigen::Index j=0;j<m.size();++j) {
                    output.push_back(m.data()[j].real()); output.push_back(m.data()[j].imag());
                }
            }
        }
        const double link_time=seconds(link_start);
#ifndef PREVIOUS
        const auto new_evaluations=cache->evaluations()-before;
#endif
        if (i==0) std::cout << "cold " << seconds(all) << '\n';
        std::cout << "step " << i << ' ' << bind_time << ' ' << link_time << ' '
                  << graphs << ' ' << embeddings << ' ' << new_evaluations;
        for (double x : output) std::cout << ' ' << x;
        std::cout << '\n' << std::flush;
    }
    rusage usage{}; getrusage(RUSAGE_SELF,&usage);
    std::cout << "total " << seconds(all) << ' ' << usage.ru_maxrss << '\n';
    return 0;
}
