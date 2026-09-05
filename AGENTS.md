# Working on pcut

This is a C++20 numerical perturbative continuous unitary transformation (pCUT)
and linked-cluster expansion library. Read README.md and docs/design.md first.
The original research sources and their provenance live in references/.

## Physics contracts

- Standard pCUT assumes H = E_ref + Delta Q + lambda V, Delta > 0,
  integer nonnegative Q, and a product vacuum for the supplied LCE drivers.
  Do not claim support for unperturbed spectra outside that class.
- Universal coefficients use the quasiparticle/sign generator, not the Wegner
  generator. Preserve operator order: the rightmost T acts first.
- Keep universal coefficients exact; convert to floating point only when
  evaluating model matrix elements. Never infer coefficients by numerical fits.
- Subtract the vacuum contribution before linking one-particle matrix elements.
  Raw fixed-particle Hamiltonian blocks are not cluster additive.
- Count every embedded connected subcluster with its multiplicity. Preserve
  interaction type and ordered local-operator legs. White graphs are out of scope.
- Never use a finite-size approximation as an infinite-lattice linked result.

## Engineering

- Use Eigen for linear algebra, Boost.Multiprecision for exact rational arithmetic,
  and Catch2/CTest for tests. Prefer existing libraries over replacement utilities.
- Public interfaces are in include/pcut; implementations in src. Keep model
  examples separate from model-independent algorithms. Check invalid inputs and
  resource limits explicitly. No silent truncation of intermediate Q sectors.
- Maintain meaningful tests: exact coefficient fixtures, analytic toy models,
  linked-cluster cancellation and embedding checks, and literature benchmarks.
- Build and test with `cmake --preset release`, `cmake --build --preset release`,
  and `ctest --preset release`. Use the sanitize preset for memory/UB validation.
- Document units, basis ordering, series convention, limits, and provenance.
  Do not edit downloaded reference sources. Do not commit build products.
