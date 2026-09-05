#pragma once
#include <pcut/lattice.hpp>

namespace pcut::models {
// Local order: empty, up, down, up-down = c_up^dagger c_down^dagger |0>.
[[nodiscard]] LocalSpace hubbard_site();
// Signed dimensionless amplitude a: a * sum_sigma(c0^dagger c1 + h.c.).
[[nodiscard]] Matrix hubbard_hopping(double amplitude = -1.0);
// H/U = Q + (t/U) V. Return dimensionless models with gap=1.
[[nodiscard]] PeriodicLattice hubbard_chain(double amplitude = -1.0);
[[nodiscard]] PeriodicLattice hubbard_square(double amplitude = -1.0);
// Dimer basis: singlet, t+, t0, t-. Energies -3J/4, J/4, J/4, J/4.
[[nodiscard]] LocalSpace spin_dimer(double J = 1.0);
// V on neighboring dimers = J [S_R.S_L + alpha(S_L.S_L + S_R.S_R)].
[[nodiscard]] Matrix dimer_bond(double alpha, double J = 1.0);
[[nodiscard]] PeriodicLattice dimerized_chain(double alpha = 0.0, double J = 1.0);
// H = sum_i n_i - lambda * coupling sum_<ij> X_i X_j.
[[nodiscard]] PeriodicLattice ising_chain(double coupling = 1.0);
}
