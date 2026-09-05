# Validation and reproducible benchmarks

## Automated checks

Run `ctest --preset release` and `ctest --preset sanitize`. The suite includes:

- All 602 tabulated exact C(word) entries extracted from the 1999 reference's
  Appendix C, through order six, plus the adjoint symmetry of every generated
  nonzero term. The data are independent of the recurrence implementation.
- A complex two-level Hamiltonian with analytic coefficients through order six,
  including a non-unit unperturbed gap and a nonzero reference energy.
- General charge changes, tensor-leg order, multi-site operators, Hermiticity,
  and state-encoding and dense-matrix representation limits.
- A finite three-dimer spectrum compared with Eigen exact diagonalization. The
  fourth-order vacuum truncation error decreases with the expected fifth-order
  scaling as lambda is halved.
- Connected square-lattice bond-animal counts 2, 6, 22, 88 through four edges,
  including the plaquette; chain subinterval multiplicities; scalar cancellation;
  on-site interactions; and multiple cell sites with complex directed hopping.
- Independent low-order colored edge-animal enumeration checks white-graph
  embedding sets and linked vacuum coefficients with unequal/zero couplings.
  Arbitrary vertex permutations, star automorphisms, straight/bent paths,
  nonisomorphic graphs, channel reordering, ordered hyperedges, parallel templates,
  multiple species, monomial/subcluster maps, and graded Hubbard blocks are checked.
  Cache counters verify reuse and invalidation without timing assertions.
- The transverse-field Ising chain with `H=sum n-lambda sum X_i X_(i+1)`.
  Its exact dispersion is `sqrt(1+4 lambda^2-4 lambda cos k)`. Energy and
  dispersion agree through fourth order at multiple momenta.
- General-sector linking agrees with the independent specialized Q=1 driver.
  Spectator subtraction removes disconnected dressed particles. Three-body
  interactions and charge-two / two-charge-one conversion kernels are retained.
  The Ising second-neighbor hopping `-lambda^2/2` and its Jordan-Wigner correlated
  hopping `+lambda^2` are both reproduced.
- Dimer-chain ground-state energy through sixth order at alpha=0, 0.17, 0.5.
  On-site and all nonzero-distance dispersion coefficients through sixth order
  match Appendix D at alpha=0, 0.23, 0.5 for all three triplet flavors, using the
  normalization discussed below.

The local verification used Apple Clang 17 on arm64 macOS with Eigen and Boost
from Homebrew. Both Release and address/undefined-behavior sanitizer configurations
pass. A separate CMake consumer successfully builds and runs against the installed
`pcut::pcut` package. Linux Release and sanitizer checks are configured in GitHub
Actions; they have not been run remotely as part of this local implementation.

## Performance comparison

Seven alternating matched comparisons with `a22f963` show complete coupling
sweeps faster by 8.0× (square Ising), 4.4× (four-color Ising), 4.9× (Hubbard),
and 19.5× (dimer). Hubbard and dimer also have clear cold-run improvements;
Ising cold ranges overlap and peak-memory medians increased. See the
[optimization report](performance.md) for timings, variability, memory and cache
counts. The [earlier review](performance-review.md) records the initial redesign's
regressions.

## Dimer-chain conventions

The implemented Hamiltonian, with dimer i containing spins L and R, is

```
H = J sum_i S_(i,L).S_(i,R)
  + J lambda sum_i [S_(i,R).S_(i+1,L)
                   + alpha (S_(i,L).S_(i+1,L) + S_(i,R).S_(i+1,R))].
```

`Q` counts triplets, `Delta=J`, and each dimer contributes `-3J/4` to the bare
vacuum energy. The paper's ground-state expression `E_grund` is per spin and
shifts the unperturbed vacuum to zero. To compare, restore `-3J/8` per spin and
multiply by two to obtain energy per dimer. Its expansion variable is
`lambda_bar=lambda/4`.

For alpha=0, our energy per dimer starts

```
E0/J = -3/4 - 3 lambda^2/32 - 3 lambda^3/128
       - 13 lambda^4/2048 - 89 lambda^5/24576
       - 463 lambda^6/196608 + O(lambda^7).
```

At alpha=1/2 all vacuum corrections vanish: the dimer-singlet product remains an
exact eigenstate. This is checked coefficient by coefficient and does not assume
that the state stays the ground state at arbitrary large coupling.

For the original bond-alternating spin-chain parameters, the reference defines
`J=J0(1+delta)`, `lambda=(1-delta)/(1+delta)`, and
`alpha=alpha0/(1-delta)`. The executable takes **alpha**, not alpha0. A fixed-alpha
lambda series is therefore not a fixed-alpha0 dimerization sweep.

## Hopping normalization in the supplied reference

There is an internal factor-of-two inconsistency in cond-mat/9906243v1. Its
`hopp_coef` and `disp_2` define directed hopping `a_r` and
`omega(k)=a_0-E0+2 sum_(r>0) a_r cos(r k)`. But Appendix D starts with
`a_1=-2(1-2alpha) lambda_bar`, which is twice the direct matrix element

```
<t,s| V |s,t> = -(1-2alpha)/4,
```

in units J=1. Our spin-matrix construction independently gives this matrix
element. Every Appendix D nonzero-distance polynomial through order six agrees
with **2 times** our directed hopping; the on-site polynomial agrees without that
factor. We therefore interpret the appendix's nonzero-distance entries as cosine
amplitudes. This is an inference from the supplied source and independent spin
algebra, not a claim that an external erratum has been located.

The fixture `tests/data/chain_hopping.tsv` preserves the appendix polynomials.
Tests explicitly divide the r>0 entries by two; the library applies no corrective
factor. For example, at alpha=0 the computed k=0 dispersion begins
`1-lambda/2-3 lambda^2/8+lambda^3/32`, consistent with this interpretation.

## Eighth-order benchmark

```sh
python3 scripts/benchmark.py --executable build/release/pcut_chain \
  --output build/benchmark.json
python3 scripts/verify_references.py
```

The benchmark compares the vacuum series to the original expression through
order eight and records coefficients, absolute errors, and wall-clock times.
The current white-graph result is stored in [benchmark.json](benchmark.json).
All four runs passed; the largest vacuum coefficient error was below `1e-16`,
and the eighth-order run at alpha=0.17 had error below `1.4e-17`. The report
records wall-clock times from serial runs after sanitizer checks finished;
these are measurements, not portable performance guarantees. Symbolic monomial
growth can make high-order multi-channel calculations more costly even when
graph isomorphism reduces the number of evaluations. The square-lattice sweep
example demonstrates 10 evaluations for 30 colored embeddings through order
three, with zero new graph evaluations on later sweeps, independent of timing.

These checks establish the finite-order coefficients and embedding conventions.
They do not establish convergence at lambda=1, behavior across a phase transition,
or accuracy of unimplemented resummation methods.

## Hubbard validation

`tests/hubbard.cpp` checks the full no-doublon operator through fourth order.
Independent global Fock creation matrices validate nonadjacent and reversed-leg
hopping, exchange and three-site correlated hopping. The dimer matches its
analytic singlet and one-electron Hamiltonians. Operator subtraction reconstructs
connected and disconnected finite systems, including permuted square embeddings
at half and quarter filling. All low-energy eigenvalues of independently built
four-site Hubbard Hamiltonians converge to the truncated expansion as `t/U`
decreases.

Exact rational word fixtures match Chernyshev et al.'s fourth-order CT1 result
with the sign-generator rotation `gamma=1/4`. The half-filled matrix also matches
Delannoy et al.'s `Hs4`. Infinite-square weights yield `J1=4x²-24x⁴`,
`J2=J3=4x⁴`, `Jc=80x⁴` in units U, by Pauli trace projection, without fitting.
The full connected plaquette matrix, including pair corrections and its constant,
is checked. [The Hubbard guide](hubbard.md) specifies the effective-Hamiltonian
convention and the distinction between infinite couplings and finite spectra.
