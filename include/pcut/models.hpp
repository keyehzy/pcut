#pragma once
#include <pcut/lattice.hpp>

namespace pcut::models {
// Dimer basis: singlet, t+, t0, t-. Energies -3J/4, J/4, J/4, J/4.
[[nodiscard]] LocalSpace spin_dimer(double J = 1.0);
// V on neighboring dimers = J [S_R.S_L + alpha(S_L.S_L + S_R.S_R)].
[[nodiscard]] Matrix dimer_bond(double alpha, double J = 1.0);
[[nodiscard]] PeriodicLattice dimerized_chain(double alpha = 0.0, double J = 1.0);
// H = sum_i n_i - lambda * coupling sum_<ij> X_i X_j.
[[nodiscard]] PeriodicLattice ising_chain(double coupling = 1.0);
}
