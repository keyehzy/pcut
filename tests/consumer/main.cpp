#include <pcut/pcut.hpp>
#include <cmath>
int main() {
    const auto lattice=pcut::models::ising_chain();
    const pcut::Coefficients coefficients(pcut::charge_changes(lattice),2);
    const auto result=pcut::linked_expand(pcut::WhiteGraphExpansion(lattice,2),pcut::EffectiveOperator(coefficients));
    if (std::abs(result.energy_per_cell[2].real()+0.5)>=1e-12) return 1;
    const auto hubbard=pcut::models::hubbard_chain();
    const pcut::EffectiveOperator effective(pcut::Coefficients({-1,0,1},2));
    const auto operators=pcut::linked_zero_charge(pcut::WhiteGraphExpansion(hubbard,2),effective);
    const pcut::Cluster edge{{0,{0}}};
    const auto half=pcut::assemble_operator(operators,pcut::vertices(hubbard,edge),edge,2);
    return half.basis.size()==4 && std::abs(half.coefficients[2].trace().real()+4)<1e-12 ? 0 : 1;
}
