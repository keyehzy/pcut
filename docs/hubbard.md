# Hubbard strong-coupling expansion

The repulsive single-band model is

```
H = U Q + t V,  U > 0,
Q = sum_i n_i,up n_i,down,
V = -sum_<ij>,sigma (c_i,sigma^dagger c_j,sigma + h.c.).
```

`models::hubbard_chain()` and `models::hubbard_square()` supply dimensionless
nearest-neighbor hopping with signed amplitude `-1`. Their optional amplitude
argument changes that dimensionless signed ratio, not `U`. Every undirected bond
occurs once, and its matrix includes both directions and both spins. They return
`gap=1`: compute `H_eff/U = sum_n A[n] (t/U)^n` and multiply the evaluated matrix
by the physical positive `U`. Thus order one has units `t`, order two `t²/U`,
and order four `t⁴/U³`. There are no factorials. Do not set `gap=U` while also
passing `t/U` as the expansion parameter with dimensionless hopping.

The exact `Coefficients({-1,0,1}, order)` engine uses the quasiparticle/sign
flow generator. Only matrix-element evaluation converts its rationals to double.
No fitting or finite-size extrapolation is used. The low-energy block retains
**every** state with no doublons, including empty sites, either spin, all spin
sectors and all electron numbers. Degenerate `Q=0` does not require a unique
reference spin state or a spin gap. The unperturbed charge ladder remains integer
and nonnegative; arbitrary non-equidistant unperturbed spectra are unsupported.

## Local basis and fermionic embedding

`LocalSpace` now allows degenerate charge zero and optional `particles` and
`parity` arrays of the same length as `charges`. Hubbard's ordered local basis is

| index | state | Q | electrons | parity |
| --- | --- | --- | --- | --- |
| 0 | empty | 0 | 0 | 0 |
| 1 | up | 0 | 1 | 1 |
| 2 | down | 0 | 1 | 1 |
| 3 | up-down | 1 | 2 | 0 |

The orbital order is `(site 0 up, site 0 down, site 1 up, site 1 down, ...)`.
A state is the product of occupied creation operators in **increasing orbital
order**, acting on the empty Fock vacuum; locally `|up-down> = c_up† c_down† |0>`.
Site zero is the least significant mixed-radix digit. For Hubbard the state ID
is also the bit mask of occupied orbitals. An annihilator/creator contributes
`(-1)^(number of occupied preceding orbitals)`. Matrix columns are inputs and
rows are outputs. In an operator word, the rightmost operator acts first.

Set `LocalTerm::fermionic` / `OperatorChannel::fermionic` for an even Hamiltonian
matrix supplied in the Fock order of its **ordered tensor legs**. The library
gathers those legs to the front, followed by spectators in original order.
For each inversion of two odd occupied local states it contributes a minus sign,
on both the input and the output. This is a graded permutation of whole site
Fock spaces; it supplies nonadjacent parity strings and handles arbitrary site
permutations, also in two dimensions. Fermionic matrices must conserve total
parity, and every cluster site must have parity metadata. Do not include an
additional global Jordan-Wigner string in a matrix marked `fermionic`.

When parity metadata is present, the operator-valued driver requires it on every
site and rejects unmarked tensor terms that change any individual leg's parity.
Parity-preserving density/spin terms can use either embedding convention.

Ordinary tensor matrices keep `fermionic=false`, preserving the existing spin
models. Explicitly supplied strings remain usable in that convention. The old
vacuum, one-particle and tensor-sector drivers require unique local charge zero
and reject Hubbard spaces. The tensor-sector driver also rejects graded terms;
its spectator contraction is not a fermionic normal-ordering rule.

## Operator-valued linking

`zero_charge_operator` evaluates `P H_eff P`, where `P` is the full tensor product
of local charge-zero projectors. It never truncates intermediate charge sectors.
`linked_zero_charge` reuses a canonical white-graph symbolic block for each
connected physical embedding `C`, then forms

```
W(C) = P_C H_eff(C) P_C - E_ref(C) I
       - sum_(proper connected embedded S subset C) embed(W(S), S -> C).
```

`embed` includes the identity on **all** spectator states, with graded permutation
signs. Each actual edge subset is subtracted with its own vertex map; equivalent
subsets at different positions retain multiplicity. Only the bare order-zero
reference constant is treated separately. There is no interacting vacuum-energy
subtraction and no particle-irreducible contraction.

Why this is additive: for disjoint systems A and B, even operators on A and B
commute. Starting with `H_A + H_B`, the sign generator is `eta_A + eta_B` because
each local matrix element's total charge change is its subsystem charge change.
The flow therefore stays additive order by order. Since charges are nonnegative,
`P_(A union B) = P_A tensor P_B`, and the projected effective Hamiltonian is
`H_eff,A tensor I_B + I_A tensor H_eff,B` under the graded identification of Fock
spaces. Connected operator subtraction applies directly to this object. A block
at a fixed electron number on each *subcluster* does not have this tensor-product
property and cannot replace it.

The output stores translation classes of operator matrices and their colored
edge supports. Sum each weight over all unit-cell translations once. The catalog
retains orientations separately, so no automorphism division is needed. Clusters
through `n` distinct edges contain every contribution of order `n`. These are
infinite-lattice effective **couplings**, independent of filling; they are not
an infinite-lattice spectrum or a finite-size estimate of one.

## Reusable API and finite evaluation

```cpp
#include <pcut/pcut.hpp>

const auto lattice = pcut::models::hubbard_square();
const pcut::Coefficients coefficients({-1,0,1}, 4);
const pcut::EffectiveOperator effective(coefficients);
const pcut::WhiteGraphExpansion catalog(lattice, 4);
const auto linked = pcut::linked_zero_charge(catalog, effective);

// One open square; template 0 is horizontal, template 1 vertical.
const pcut::Cluster edges{{0,{0,0}}, {1,{1,0}}, {0,{0,1}}, {1,{0,0}}};
const auto sites = pcut::vertices(lattice, edges);
const int electrons = 2; // quarter filling: Ne/L=1/2; half filling is 4
const auto finite = pcut::assemble_operator(linked, sites, edges, electrons);
const double U = 8.0, t = 0.4;
const pcut::Matrix H = U * finite.evaluate(t/U);
```

`assemble_operator` restricts the translation sum to embeddings whose colored
edges all occur in the supplied finite system. The `sites` vector specifies its
Fock order, may be permuted, and may include isolated spectators. Electron number
is selected only here, after retaining full local spaces in every weight.
`OperatorBlock::basis` gives the ascending original state IDs for interpreting
rows and columns. Empty particle sectors produce zero-dimensional blocks.
Particle selection requires metadata and particle-number conservation.

For a direct finite calculation, construct a `ClusterModel` with
`LocalTerm{{i,j}, models::hubbard_hopping(), true}` on each bond and call
`zero_charge_operator(model, effective, electrons)`. The reusable
`add_embedded_operator` adds an even child operator with an arbitrary injective
site map and identity spectators. Its child basis must contain the entire
ordered charge-zero manifold; a parent may be restricted to a particle sector.

```sh
cmake --preset release
cmake --build --preset release
build/release/pcut_hubbard --lattice chain --filling half --order 4 --ratio 0.05 --U 8
build/release/pcut_hubbard --lattice chain --filling quarter --order 4 --ratio 0.05 --U 8
build/release/pcut_hubbard --lattice square --filling half --order 4 --ratio 0.05 --U 8
build/release/pcut_hubbard --lattice square --filling quarter --order 4 --ratio 0.05 --U 8
```

The executable builds infinite-lattice weights, reports the finite assembled
coefficient norms, and diagonalizes a four-site open chain or plaquette. Its
lowest eigenvalue is that of a **finite truncated Hamiltonian**, not a Taylor
series for a thermodynamic ground-state energy. The library does not choose a
phase or solve the thermodynamic projected many-electron Hamiltonian.

## Fourth-order convention and validation

At second order the full operator contains projected hopping, the exchange
`4t²/U (S_i.S_j - n_i n_j/4)`, and the three-site term

```
-(t²/U) sum_(i != k neighbors of j, sigma)
 P [c_i,sigma† n_j,-sigma c_k,sigma
    - c_i,sigma† c_j,-sigma† c_j,sigma c_k,-sigma] P.
```

All third- and fourth-order terms are evaluated as operator words, including
charge-motion terms away from half filling. They are not replaced by a spin-only
ansatz. Half filling on the square lattice gives the published spin coefficients
`J1 = 4t²/U - 24t⁴/U³`, `J2 = J3 = 4t⁴/U³`, and `Jc = 80t⁴/U³` multiplying

```
(S_i.S_j)(S_k.S_l) + (S_i.S_l)(S_j.S_k) - (S_i.S_k)(S_j.S_l)
```

on each cyclically ordered square `(i,j,k,l)`. The library retains constant
shifts; the fully polarized half-filled state has zero energy. In particular,
the connected fourth-order square weight is `80 * ring - 4 * sum_(six pairs)
S_a.S_b + I`. This distinguishes the spin-product `Jc=80` convention from
coefficients multiplying a cyclic permutation operator.

Away from half filling, compare effective Hamiltonians in the same unitary
convention. In dimensionless units let `A=T_-1 T_1`. Relative to Chernyshev et al.'s
`CT1` expression (`AM` in the downloaded source), the sign-generator result is

```
H4(sign) = H4(CT1) + 1/4 * (T0² A + A T0² - 2 T0 A T0).
```

This is their charge-preserving rotation with
`S0 = (T0 A - A T0)/4`, using `H' = exp(-S0) H exp(S0)` and implicit order-three
scaling of `S0`. The difference vanishes in the half-filled external manifold.
The tests check these word coefficients **exactly** and compare independently
constructed full-Fock products against the complete matrix, including doped
sectors. They also verify the half-filled expression `Hs4` of Delannoy et al.
See [primary-source provenance](../references/README.md).

Additional tests cover nonadjacent hopping and reversed legs, permuted 2D
embeddings, the analytic dimer singlet `(U-sqrt(U²+16t²))/2`, the complete
second-order three-site formula, disconnected operator cancellation and repeated
subcluster embeddings. Independent finite Hubbard diagonalizations compare
**all** low-energy eigenvalues at `Ne=2,4` on four-site chains and squares; errors
decrease at the expected order as `t/U` is halved. Infinite-square spin couplings
are extracted by orthogonal Pauli traces of linked weights, not numerical fits.

## Computational cost and scope

The operator APIs impose no configurable basis, storage, traversal or embedding
work budgets. They retain complete requested bases and linked weights; allocation
failures propagate instead of discarding states, sectors, words or clusters.
State IDs must fit 64 bits, and dense matrix sizes must fit Eigen indexing and
byte-size representation. See [design.md](design.md) for implementation limits.

Dense external operator matrices scale as `9^L` for Hubbard, so high orders become
costly. Callers choose calculations appropriate to their memory and available run
time. The core APIs accept higher orders; Hubbard physics is validated through
fourth order and the example limits its order to 1–4. Finite assembly supports
explicit open lattice edge sets, not periodic boundary identifications or a
large-system sparse many-body eigensolver. Custom hopping geometries can be
supplied through `PeriodicLattice` or direct `ClusterModel` terms.

No resummation or convergence guarantee at large `t/U`, transformed observables,
or thermodynamic phase solver is supplied. White-graph caching supports arbitrary
numerical hopping ratios with the same fixed graded channel operators. Preserve the chosen
unitary convention when using these couplings in another many-body solver.
