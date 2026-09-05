# pcut

A C++20 numerical perturbative continuous unitary transformation (pCUT) and
linked-cluster expansion library. It generates universal coefficients exactly,
applies local operators to sparse product states, and embeds connected cluster
weights directly into an infinite periodic lattice. White graphs are not used.

The implementation supports the standard pCUT class

\[
H(\lambda)=E_{\rm ref}+\Delta Q+\lambda V,\qquad \Delta>0,
\]

where `Q` has nonnegative integer local eigenvalues. The vacuum and particle
drivers require a unique product vacuum; the dedicated operator-valued driver
also supports a degenerate local charge-zero manifold.
Models supply local Hilbert spaces and complex Hermitian interaction matrices.
Interactions may have arbitrary finite support, ordered tensor legs, different
coupling strengths, and different site species. Periodic lattices can have
multiple sites per cell and arbitrary spatial dimension within configured limits.
Non-equidistant spectra require a different perturbative scheme.

## Build and run

Requires a C++20 compiler, CMake 3.25+, Ninja, Eigen 3.4+, and Boost 1.74+.
Catch2 3 is used if installed; otherwise CMake downloads the pinned, checksummed
Catch2 3.8.1 source archive. Eigen and Boost are system dependencies.

```sh
# macOS
brew install cmake ninja eigen boost
# Ubuntu/Debian alternative
sudo apt-get install cmake ninja-build g++ libeigen3-dev libboost-dev

cmake --preset release
cmake --build --preset release -j 6
ctest --preset release

build/release/pcut_chain --order 6 --alpha 0.17 --lambda 0.3 --k 0
build/release/pcut_chain --order 8 --alpha 0.17 --vacuum-only
```

The example prints CSV coefficients and an evaluated truncated series. `--help`
lists options. `alpha` is the frustration ratio in the dimer expansion, not the
original chain's `alpha_0`. Momenta use dimer spacing; energy includes the physical
reference `-3/4` per dimer (`-3/8` per spin), with strong exchange `J=1`.
The value at `k=0` is a dispersion value; the code does not assume that this is
the global gap for every frustration or coupling.

```sh
cmake --preset sanitize
cmake --build --preset sanitize -j 6
ctest --preset sanitize
```

For an offline build, install Catch2 3 or disable tests with `-DBUILD_TESTING=OFF`.
For a nonstandard dependency prefix, pass `-DCMAKE_PREFIX_PATH=/your/prefix`.

## Library use

```cpp
#include <pcut/pcut.hpp>

const auto lattice = pcut::models::dimerized_chain(0.17);
const unsigned order = 6;
const pcut::Coefficients coefficients(pcut::charge_changes(lattice), order);
const pcut::EffectiveOperator effective(coefficients);
const pcut::ClusterCatalog clusters(lattice, order);
const auto result = pcut::linked_expand(clusters, effective);

const auto energy = pcut::evaluate(result.energy_per_cell, 0.3);
const auto bands = result.bloch({0.4}, 0.3); // Hermitian matrix, one row per flavor
```

`Series[n]` is the coefficient of `lambda^n`, without factorials. `Coefficients`
uses Boost arbitrary-precision rationals; model evaluation uses complex doubles.
Tables and effective programs can be reused across geometries and coupling sweeps
that share a charge-change alphabet.

For the repulsive spinful Hubbard model, `models::hubbard_chain()` and
`models::hubbard_square()` use `H/U = Q + (t/U) V`. The dedicated
`linked_zero_charge` API returns operators throughout the no-doublon manifold;
`assemble_operator` selects electron number in a finite assembled problem.
Half filling is `Ne/L=1`, quarter filling is `Ne/L=1/2`. Projected hopping,
exchange, three-site processes and all contributions through fourth order,
including square ring exchange, are validated. See the [Hubbard guide](docs/hubbard.md)
for conventions, public APIs, resource limits, and infinite-coupling versus
finite-spectrum distinctions.

```sh
build/release/pcut_hubbard --lattice chain --filling half --order 4
build/release/pcut_hubbard --lattice square --filling quarter --order 4
```

For new models, define `LocalSpace`, `Interaction`, and `PeriodicLattice`.
For finite systems, use `ClusterModel` directly. `EffectiveOperator::apply` and
`block` work in any chosen charge sector, retaining all reachable intermediate
sectors. `linked_expand` computes vacuum energy and one-particle hopping;
`linked_expand_sectors(clusters, effective, max_charge)` additionally produces
irreducible multi-particle and charge-conversion kernels. `linked_scalar` accepts
a custom cluster-additive scalar evaluator. See [the model guide](docs/models.md)
and [the algorithm description](docs/design.md) for contracts and examples.

```sh
cmake --install build/release --prefix "$PWD/build/install"
```

Downstream CMake projects use `find_package(pcut CONFIG REQUIRED)` and link
`pcut::pcut`. An example package consumer lives in `tests/consumer/`.

## Validation and scope

The tests check 602 exact reference coefficients, analytic two-level and Ising
models, Hubbard operators through fourth order at half and quarter filling,
finite-system exact diagonalization, graph counts and embedded
multiplicities, complex hopping, multi-particle kernels, and the dimerized /
frustrated chain's vacuum and dispersion series through sixth order.
[Validation details](docs/validation.md) explain units and a factor-of-two
inconsistency in the reference's hopping notation. The reproducible benchmark
script also checks the vacuum series through eighth order.

This is a finite-order expansion, without resummation or a convergence guarantee
at large `lambda`. Computational cost grows exponentially with perturbation order,
charge alphabet, and cluster complexity. Explicit budgets reject oversized jobs.
The default universal-word budget accommodates order eight for `{-2,-1,0,1,2}`;
the library accepts larger user budgets. Sparse tensor-state IDs are 64-bit.
Finite local matrices must be supplied, including any truncation of bosonic spaces
and statistics metadata for graded fermionic terms. Couplings are numerical ratios
multiplying one formal expansion parameter; symbolic multivariate polynomials, transformed
observables, graph-isomorphism caching, and white graphs are not implemented.

## References

Original LaTeX sources, archives, and SHA-256 checksums are in
[references/](references/README.md).

- [Coester and Schmidt, Optimizing linked cluster expansions by white graphs](https://arxiv.org/abs/1505.02975).
- [Knetter and Uhrig, Perturbation Theory by Flow Equations: Dimerized and Frustrated S=1/2 Chain](https://arxiv.org/abs/cond-mat/9906243).
- [Knetter, Schmidt and Uhrig, The Structure of Operators in Effective Particle-Conserving Models](https://arxiv.org/abs/cond-mat/0306333).

- [Chernyshev et al., Higher order effective low-energy theories](https://arxiv.org/abs/cond-mat/0407255).
- [Delannoy et al., Néel order, ring exchange and charge fluctuations in the half-filled Hubbard model](https://arxiv.org/abs/cond-mat/0412033).
