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
  operator channels, statistics, site species and ordered local-operator legs.
  Use canonical white graphs with explicit vertex/edge/channel maps; numerical
  couplings enter only at embedding. Preserve duplicate interaction templates.
- Cache immutable symbolic graph evaluations by the complete physical/operator
  context and exact coefficient program. Never identify graphs by geometry or
  numerical coupling values, and never substitute a hash for an exact identity.
- Never use a finite-size approximation as an infinite-lattice linked result.

## Engineering

- Use Eigen for linear algebra, Boost.Multiprecision for exact rational arithmetic,
  and Catch2/CTest for tests. Prefer existing libraries over replacement utilities.
- Public interfaces are in include/pcut; implementations in src. Keep model
  examples separate from model-independent algorithms. Check invalid inputs and
  representation limits explicitly. Let allocation failures propagate; do not add
  configurable storage or work budgets. No silent truncation of intermediate Q sectors.
- Maintain meaningful tests: exact coefficient fixtures, analytic toy models,
  linked-cluster cancellation and embedding checks, and literature benchmarks.
- Build and test with `cmake --preset release`, `cmake --build --preset release`,
  and `ctest --preset release`. Use the sanitize preset for memory/UB validation.
- Document units, basis ordering, series convention, limits, and provenance.
  Do not edit downloaded reference sources. Do not commit build products.

## Other
- Do not keep any backwards compatibility, this is a single user library and will
  likely slow down our fast-paced development.
- Prefer cloning/download sources locally (e.g. to /tmp) instead of open indivisual
  pages online.