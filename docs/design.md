# Algorithms and contracts

## Universal pCUT

Write `V = sum_m T_m` with `[Q,T_m]=m T_m`. We use the quasiparticle generator
`eta_ij = sign(Q_i-Q_j) H_ij`. For a word `w` with total charge `M(w)`, the
coefficient functions satisfy the 1999 reference's equation `dgl`:

```
dF_w/dl = -abs(M(w))*F_w
           + sum_(w=u.v) [sign(M(u))-sign(M(v))] F_u F_v.
```

Length-one initial functions are `exp(-abs(m)*l)`; longer ones start at zero.
`coefficients.cpp` represents each function as a finite sum of polynomials in `l`
times `exp(-r*l)`. Products, integration, and cancellation use exact
Boost.Multiprecision rationals. The effective coefficient is the constant
zero-decay term at infinity for `M(w)=0`. A nonconstant zero-decay term is an
internal error. Words outside the preserved block bandwidth are omitted.

Coefficients are generated for the actual sorted change alphabet, which need not
be consecutive or bounded by two. The implementation generates coefficients
rather than fitting spectra or numerically stopping a flow at a finite time.

A reversed-word trie shares rightmost operator applications in the numerical
finite-cluster evaluator. The symbolic evaluator further factors its continuation
program: at each depth, a sparse exact-rational row represents a terminal
coefficient and its `(next charge change, next-layer basis)` continuations.
Gaussian elimination produces a basis for that row space. Propagating into these
basis vectors combines linearly dependent continuations before further local
applications. This is an exact algebraic compilation of the original coefficient
table, not a fit; cache identity retains the complete original rational program.

Numerical and symbolic evaluation share one internal transition enumerator for
local-state extraction, physical output encoding and fermionic signs. Accumulation
remains separate; exact program factorization and monomial interning have separate
implementations. Only currently nonzero product states propagate through compiled
local transition lists. Every reachable intermediate charge sector is retained; only the final
external basis restricts output. Each local transition is divided by `Delta`
before multiplication, and the final sum restores one factor of `Delta`.
This avoids unit-dependent overflow/underflow in powers of physical energies;
dimensionless model products must still fit double precision. Rational program
weights are converted when evaluating model amplitudes. Exactly real operators
use real scratch arithmetic; operators with any nonzero imaginary entry use
complex arithmetic. Only exact zero amplitudes are removed.

Operator tables and coefficient programs are immutable. Simultaneous evaluations
use independent scratch storage. Reached transitions are memoized within each
external input column, so temporary storage is released/reused between columns.
The default drivers execute serially.

## Local Hilbert spaces

`LocalSpace::charges[0]` is zero and every charge is a nonnegative integer.
Vacuum-dependent drivers explicitly require every other charge to be positive;
the operator-valued charge-zero driver admits local degeneracy.
The physical on-site spectrum is `vacuum_energy + Delta*charges[i]`, with a common
positive `Delta`. Thus degeneracies and missing levels are allowed; arbitrary
incommensurate spacings are not. A tensor state uses mixed-radix 64-bit encoding.
The entire tensor dimension must fit that encoding even when only a few states
are reached. No full tensor-space matrix is allocated in the pCUT solver.

Each `LocalTerm` matrix is finite and Hermitian. Tensor leg zero is the least
significant local index. This convention applies both to matrix construction and
embedding, including reversed legs and unequal local dimensions. Optional
`particles` and `parity` metadata specify electron number and Z2 grading.
Even matrices marked `fermionic` use ordered-leg Fock bases and graded gathering
permutations to insert spectator signs. Unmarked terms use ordinary tensor
embedding. The library does not choose a bosonic cutoff. The complete Hubbard
Fock convention is in [hubbard.md](hubbard.md).

`ClusterModel::dense_hamiltonian` is a dense validation helper;
Eigen diagonalization is used only in tests and for the small Bloch matrices.

## Canonical abstract white graphs and infinite embeddings

An interaction template contains ordered sites and a list of
`CoupledChannel{OperatorChannel{matrix, fermionic}, coupling}` values. Its
Hamiltonian is `lambda * sum_channel coupling * matrix`. Ratios are real and
finite, including zero and negative values; phases belong in the Hermitian
channel matrices. Distinct physical templates remain distinct edge occurrences,
even when their support and channels coincide. Duplicate templates intentionally
add the interaction again. Multiple channels on one template share an edge but
have separate formal variables.

`WhiteGraphExpansion` grows translation classes of connected injective physical
edge sets on the infinite lattice and uses them as witnesses for canonical
abstract white graphs. Adding each incident interaction template grows every
connected set: every connected hyperedge set with more than one edge has a
connected single-edge deletion. Templates retain their distinct occurrence IDs,
including duplicates; on-site, parallel and higher-arity edges are supported.
Canonicalization removes coordinates and numerical couplings before evaluation.
Witness growth avoids generating abstract candidates that cannot embed, and
collects every physical embedding while sharing its canonical graph evaluation.
There is no finite periodic box or finite-size approximation.

Canonicalization uses sparse **nauty 2.9.3** on a lossless vertex-colored
incidence graph. Four disjoint color namespaces encode physical vertices,
edge occurrences, ordered-leg ports and channel occurrences:

- A physical vertex carries the exact `LocalSpace` serialization: species name,
  ordered charges, reference energy, particle numbers and parity.
- Each edge occurrence has its own node, colored by arity and the complete
  sorted multiset of channel matrices/statistics.
- Each ordered leg has a port adjacent to its edge node and physical vertex;
  its color contains the leg position. On-site and arbitrary-arity edges use
  exactly the same construction.
- Each channel occurrence has a separate leaf attached to its edge, colored by
  its slot in that edge's structurally sorted channel list. The parent color
  contains the full list of exact matrix dimensions, real/imaginary double
  values and fermionic flags, so a slot identifies the operator losslessly.
  Even equal channels have distinct slot colors. Input reordering of equal
  channels changes only their occurrence map; it cannot change vertex
  equivalence. Duplicate edge templates remain separate nodes.

Colors are sorted by their full strings to form nauty's ordered initial
partition. Coordinates and coupling ratios never enter it. Refinement and hashes
are not identities. Nauty's canonical `lab` supplies the canonical-to-input
physical vertex order. Given those vertex labels, sorting complete edge
structures and ordered legs, then exact channel structures, recovers explicit
canonical-to-input vertex, edge and flattened-channel maps. Input occurrence
indices break ties only between identical edges/channels; the structural
serialization is unchanged by those ties. The persistent key contains every
local-space and operator signature and every ordered leg, never just an
incidence hash, interned ID, or refinement signature. Canonical ordering differs
from the former exhaustive algorithm and is not a stable file format.

The incidence group's restriction to physical vertices is onto the white graph's
vertex automorphism group: every vertex automorphism extends by matching equal
edge occurrences and the corresponding sorted channel slots. Its kernel has
order `product_(identical ordered edges) multiplicity!`. Sorted channel-slot
colors remove internal channel permutations without restricting any physical
vertex automorphism. Once physical vertices and edge occurrences are fixed,
ordered ports and channel slots have no freedom. We collect
nauty's integer stabilizer indices using `userlevelproc`, multiply them with
Boost arbitrary-precision integers, and divide by that exact kernel. We verify
divisibility and only then check whether the vertex count fits `size_t`.
Floating-point group-size statistics are never used. A huge gadget group is
therefore allowed when its vertex quotient fits. Neither this count nor the
kernel is used to normalize embeddings or subclusters.

The exhaustive algorithm exists only in the small-graph test oracle. Production
has no fallback search or backend selection switch. The native GCC/Clang build,
allocation-failure propagation, TLS, pinned dependency, installed package and
offline build are documented in [dependencies.md](dependencies.md). Backend
selection and measured trade-offs are in [performance.md](performance.md).

Each completed physical edge set is normalized by its minimum edge origin.
Deduplicating these exact sets removes translation and automorphism overcounting.
Each stored `GraphEmbedding` therefore occurs once per unit cell. Its sites are
sorted physical lattice sites; canonicalization supplies explicit vertex, edge
and flattened channel maps retaining ordered legs and physical displacements.
Exact species/operator signatures are interned during construction, and compact
local IDs accelerate serialization and repeated canonicalization. Persistent graph
identity still contains full structural serialization, never just those IDs or
a hash. Operator validation occurs at the public boundary; generated candidates
reuse already validated structure.

An expansion owns a shared immutable plan containing the lattice/operator
structure, graphs, embeddings and compiled canonical/physical models.
`bind(ratios)` constructs numerical coupling storage directly and shares the plan.
Embedding indices provide constant-time coupling and model access.
`structure()` exposes the shared lattice with unit couplings; `bound_lattice()`
materializes a bound copy only when requested. Physical models share
immutable transition tables through validated site permutations, including
updated fermionic gathering inversions. Requested external bases and exact
coefficient programs select synchronized reusable sparse readout plans.

Subtraction maps are initialized lazily and synchronized across bindings: custom
scalar linking requests canonical maps, and general-sector linking additionally
requests physical maps. Support-projection drivers construct neither.
Every proper connected edge subset is generated by recursive deletion, retaining
separate subsets even when they canonicalize to the same child. Each
`GraphSubcluster` carries canonical child-to-parent vertex/edge/channel maps.
These maps also construct the physical subcluster maps used in linking. Entries
are ordered by edge count. An order-n operator product uses at most n distinct
edges, so graphs through n edges suffice, and drivers skip excess catalog orders.

## Symbolic evaluation and cache validity

`Monomial` is an immutable sparse exponent vector. Sixteen four-bit exponent
lanes fit inline in 64 bits; larger variable indices or powers use a shared sparse
vector. This is lossless: variables retain `size_t` range and total degree remains
at most 64. Cached `Polynomial` objects store sorted contiguous monomial/amplitude
pairs. The symbolic scratch evaluator interns monomials, caches variable
multiplications and propagates integer monomial IDs in Boost flat hash tables.
Its interning key is a lossless packed exponent vector when representable and an
exact sparse vector otherwise. Allocation failures propagate.

Applying a channel increments its formal exponent. Formal degree equals lambda
order. There is no interpolation, numerical pruning, fit or work/storage budget.
Universal coefficients and program factorization use Boost exact rationals;
model amplitudes use floating point with the energy scaling described above.

`GraphCache` owns immutable `SymbolicBlock` results. The exact cache identity
contains the actual ordered graph spaces/operators, gap, exact coefficient
program (alphabet, order and rational terms), and requested external basis.
Numerical couplings and lattice coordinates are excluded. Changes to local
spaces, reference energies, statistics, matrices, order or sectors select a new
entry automatically. External state encodings and variable assignments are
mapped consistently; fermionic matrix elements receive both input and output
Fock permutation signs. Before this readout, graphs containing any fermionic
channel require complete parity metadata and ordinary channels preserving parity
on every individual site. General mixtures remain available for direct finite
evaluation, which applies each term with its supplied convention. Requested bases are sorted in canonical coordinates
before cache lookup. `evaluations()` counts cache misses that completed, while
`hits()` records reuse; no timing assumptions are needed to demonstrate savings.
Caches synchronize evaluation and lookup and keep results alive by shared
ownership. They do not evict results; allocation failures propagate.

The vacuum/one-particle and charge-zero operator drivers evaluate the exact
**full-edge support projection**: retain only monomials containing at least one
channel from every graph edge. For an additive operator/kernel, mapped connected
subcluster subtraction is precisely this projection. A monomial supported on a
proper connected edge subset is removed once by its corresponding weight;
disconnected supports cancel by additivity of the Hamiltonian and sign generator.
This eliminates repeated numerical mapped subtraction in these drivers.

During propagation, a prefix is skipped only if its missing edge count exceeds
the number of remaining operators. This is an algebraic support/degree test,
independent of coupling values and amplitudes, and never a Q-sector cutoff.
The one-particle driver still subtracts the projected vacuum contribution from
the projected one-particle diagonal before embedding. Projection commutes with
that linear subtraction. The general sector driver retains raw blocks and its
explicit spectator and connected-subcluster subtraction. Custom scalar linking
retains polynomial subcluster subtraction with explicit channel maps.

Raw blocks and support projections have separate exact cache identities. Sparse
readout plans traverse stored nonzero matrix entries, with precomputed physical
row/column positions and both graded signs; absent entries require no sparse-map
search. Binding substitutes numerical ratios only at this stage. Algebraically
equivalent summation/factorization can change roundoff residues, which are kept
unless exactly zero.

A `ScalarEvaluator` owns an immutable, isomorphism-covariant callback context.
Copies share its cache identity; create a new evaluator when captured parameters
change. The callback receives `(WhiteGraph, gap, order)` and must return
`SymbolicSeries(order+1)` with an empty constant polynomial and degree n monomials
at order n. It can use `EffectiveOperator::symbolic_block(graph.model(gap), basis)`
to obtain effective matrix elements. Cluster additivity and covariance remain
the caller's semantic responsibility. Cache callbacks must not recursively call
the same cache. Coupling-dependent numerical callbacks have been removed.

## Linking and particle irreducibility

For an additive scalar `P`, compute `W(C)=P(C)-sum_(S proper connected subset C) W(S)`.
Each embedded subset occurs once in this sum. Summing weights over translation
classes yields the correction per unit cell. On-site reference energies are
embedded separately. A scalar callback must return symbolic corrections with an empty polynomial at
order zero and must actually be cluster additive; this semantic property is the
caller's responsibility.

For the specialized one-particle driver, first form
`K_C(i,j)=<i|H_eff|j>-delta_ij E0(C)`. Subtract embedded **irreducible** subcluster
matrices from `K_C`, never raw one-particle blocks. The resulting hopping series
is accumulated by output flavor, input flavor, and displacement. The Fourier
convention is `|k,a>=sum_R exp(i k.R)|R,a>` and
`H_ab(k)=sum_d t_ab(d)exp(-i k.d)`. The bare gap is included once per local flavor.

For general sectors, enumerate external product states with `Q <= max_charge`.
This restriction applies to external states only. Let A and B be their occupied
site/flavor configurations. The finite-cluster Hubbard kernel is

```
K(A,B) = <A|H_eff|B> - sum_(nonempty S subset common(A,B)) K(A\S,B\S).
```

A common spectator is an occupied site with the same local state in A and B.
Lower-charge columns are processed first. For local charges `{0,1,...,1}`, these
are exactly the hard-core n-particle irreducible terms described in the 2003
reference. For other charges, kernels can convert different numbers of occupied
sites while conserving Q. They use tensor Hubbard operators, not a graded
fermionic normal-ordering convention.

`linked_expand_sectors` subtracts embedded connected subcluster kernels and sums
translation-normalized creation/annihilation configurations. No factorials are
needed: configurations are site ordered, with at most one local flavor per site.
It reconstructs the effective Hamiltonian on sectors through the requested Q;
higher-sector interactions require a correspondingly larger external basis.
There are no configurable storage or work budgets. All requested external states,
spectator subsets, embedded subclusters and translation kernels are retained or
processed in full. Callers are responsible for choosing feasible calculations;
allocation failures propagate, and expensive computations have no built-in time
limit. Matrix allocation checks still reject multiplication overflow and sizes
that cannot fit Eigen indexing or byte-size representation. Matrices are
constructed individually without a temporary zero-block copy.
Tiny roundoff residues are retained, not silently pruned.

## Degenerate charge-zero operator linking

`linked_zero_charge` stores complete Q=0 operator matrices on each connected
physical embedding of a canonical white graph. It applies the full-edge support projection, equivalent to subtracting every
embedded connected-subcluster weight with identity spectators and graded signs,
without vacuum or particle-irreducible subtraction. Private embedding plans cache basis indices,
decoded states and the gathering permutation; internal loops reuse validated
bases and avoid public validation within matrix-element traversal. Even disjoint
subsystem Hamiltonians and their sign generators add, and their Q=0 projectors factor, so this operator is cluster
additive. Bare reference constants are assembled separately. Electron number is
selected only in finite evaluation, never separately in linked subclusters.
See [hubbard.md](hubbard.md) for the proof, units, fourth-order unitary convention,
public APIs and validation.

## Explicit limitations

- One formal expansion parameter, with arbitrary numerical coupling ratios
  substituted into sparse formal graph/channel monomials.
- Finite local Hilbert spaces and an integer equidistant H0 charge ladder.
  A unique product vacuum is required for the vacuum/particle drivers.
- Sparse tensor dimensions must fit 64-bit state encoding; dense matrices must
  fit Eigen indexing and byte-size representation. Supported coefficient recursion
  and cluster orders remain limited to 64. Existing local-charge, site-count and
  geometry bounds also remain implementation limits. No storage or work budgets
  restrict otherwise supported calculations, and states or clusters are never
  dropped to fit memory. Exact-rational arithmetic and double-valued matrix
  elements have different precision contracts.
- Large local matrices are currently supplied densely and compiled into sparse
  transitions; a matrix-free local-operator frontend is a possible extension.
- Nauty prunes canonical-labeling search using refinement and automorphisms;
  difficult incidence graphs and large duplicate gadget groups can still be
  expensive. Gadget counts/degrees must fit nauty's signed indices and its
  `1000 * SETWORDSNEEDED(n)` workspace length; these representation checks
  throw rather than truncate. Ordered-leg symmetries of particular matrices
  are not inferred. Graph/channel monomial growth can be exponential.
- No resummation, transformed observables, automatic statistics, or built-in
  distributed execution. Fermionic
  statistics are explicit metadata, not inferred from arbitrary tensor matrices.

## Provenance

The white-graph generalized monomials, canonical graph decomposition, and late
coupling substitution follow Coester and Schmidt (2015), sections III–IV in
`references/sources/1505.02975/white_graphs.tex`. Ordered hyperedges, species,
parallel channels, graded state mappings and cache context serialization are
implementation extensions. Downloaded research sources are unchanged.

Canonical labeling follows the official nauty/Traces guide; see
[dependency provenance](dependencies.md). The incidence encoding, exact kernel
quotient and physical occurrence maps are pcut integration code.
