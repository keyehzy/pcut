#include <pcut/pcut.hpp>
#include <cmath>
int main() {
    const auto lattice=pcut::models::ising_chain();
    const pcut::Coefficients coefficients(pcut::charge_changes(lattice),2);
    const auto result=pcut::linked_expand(pcut::ClusterCatalog(lattice,2),pcut::EffectiveOperator(coefficients));
    return std::abs(result.energy_per_cell[2].real()+0.5)<1e-12 ? 0 : 1;
}
