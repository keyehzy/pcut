#include <pcut/pcut.hpp>
#include <iomanip>
#include <iostream>

int main() {
    auto lattice=pcut::models::ising_chain();
    lattice.dimension=2;
    const auto channels=lattice.interactions[0].channels;
    lattice.interactions={{{{{0,0},0},{{1,0},0}},channels},
                          {{{{0,0},0},{{0,1},0}},channels}};
    const pcut::EffectiveOperator effective(pcut::Coefficients(pcut::charge_changes(lattice),3));
    const auto cache=std::make_shared<pcut::GraphCache>();
    const pcut::WhiteGraphExpansion topology(lattice,3,cache);
    std::cout << "Jx,Jy,energy_lambda_0.1,graphs,embeddings,new_graph_evaluations\n" << std::setprecision(12);
    for (const auto& ratios : {std::pair{0.7,-0.2},std::pair{0.0,1.3},std::pair{1.0,1.0}}) {
        lattice.interactions[0].channels[0].coupling=ratios.first;
        lattice.interactions[1].channels[0].coupling=ratios.second;
        const auto expansion=topology.bind(lattice.couplings());
        const auto before=cache->evaluations();
        const auto result=pcut::linked_expand(expansion,effective);
        std::cout << ratios.first << ',' << ratios.second << ','
                  << pcut::evaluate(result.energy_per_cell,0.1).real() << ','
                  << expansion.graphs().size() << ',' << expansion.embeddings().size() << ','
                  << cache->evaluations()-before << '\n';
    }
}
