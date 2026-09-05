# Adding a model

The minimal input is a list of local charges, a reference energy, and interaction
matrices. The unperturbed gap is shared by all sites. This example has one
hard-core excitation per site and a nearest-neighbor interaction:

```cpp
pcut::LocalSpace spin{{0,1}, 0.0, "spin"};
pcut::Matrix v = pcut::Matrix::Zero(4,4);
// Tensor index = left + 2*right. Pair creation and its Hermitian conjugate.
v(3,0) = 0.7;
v(0,3) = 0.7;
// Hopping with a complex phase and its Hermitian conjugate.
v(2,1) = pcut::Complex(0,0.2);
v(1,2) = pcut::Complex(0,-0.2);

pcut::PeriodicLattice lattice{
    1, {spin},
    {{{{{0},0}, {{1},0}}, {{{v, false}, 1.0}}}},
    1.0
};
pcut::Coefficients coefficients(pcut::charge_changes(lattice), 4);
pcut::EffectiveOperator effective(coefficients);
pcut::WhiteGraphExpansion clusters(lattice, 4);
auto one_particle = pcut::linked_expand(clusters, effective);
auto two_particle = pcut::linked_expand_sectors(clusters, effective, 2);
```

A `Site` contains `{cell_coordinates, basis_index}`. For a square lattice, use
`dimension=2` and templates with offsets `{{0,0},0}` to `{{1,0},0}` and
`{{0,0},0}` to `{{0,1},0}`. Each physical translated interaction appears once:
reverse bonds should not be added again unless they are distinct terms whose
sum is intended. Hermitian conjugation belongs inside each local matrix.
Different basis indices may have different dimensions and charges. A one-site
perturbation is a template with one leg; an n-site term uses n distinct legs.

The tensor leg order matters. For leg dimensions `d0,d1,d2`, the local matrix
index is `i0 + d0*i1 + d0*d1*i2`. With Eigen's conventional Kronecker ordering,
use `kroneckerProduct(second, first)` for this two-leg convention. The dimer-chain
implementation demonstrates the construction from spin matrices.

For a finite graph, create `ClusterModel(spaces, terms, gap)` directly:

```cpp
pcut::ClusterModel model({spin,spin}, {{{0,1},v}}, 1.0);
auto energy = effective.vacuum(model);
auto block = effective.block(model, {model.encode({1,0}), model.encode({0,1})});
auto all_kernels = pcut::irreducible_sectors(model, effective, 2);
```

`block` returns exactly the requested matrix elements; it does not assert that
the supplied basis spans a complete Q sector. Use `apply` to inspect all output
states. Neither API projects intermediate states onto the external sector.

The general linked result stores `Kernel{output,input}` keys. Each excitation
contains a lattice site and nonzero local-state index. With `e^dagger_(i,a)`
equal to `|a>_i<0|`, a key gives the coefficient of the product of its creation
operators followed by its annihilation operators, summed over cell translations.
When a site occurs in both lists, the local factor is `|a><b|`. Other occupied
sites are spectators. `energy_per_cell` stores the vacuum term separately.

Use separate channels for independently varying operator structures. For example,
`{{{pair_matrix, false}, pair_ratio}, {{hop_matrix, false}, hop_ratio}}`
represents two channels on the same ordered support. Matrices are fixed Hermitian
operators in physical energy units; ratios are dimensionless real numbers and may
be unequal, negative or zero. A matrix's complex phases remain inside the matrix.
All channels multiply the same lambda. Duplicate templates and duplicate channels
are separate occurrences whose Hamiltonian contributions add.

For coupling sweeps, reuse the exact coefficient program and the expansion plan:

```cpp
const pcut::WhiteGraphExpansion plan(lattice, 4);
for (double coupling : {0.0, 0.4, 0.6}) {
    auto bound = lattice;
    bound.interactions[0].channels[0].coupling = coupling;
    const auto graphs = plan.bind(bound.couplings());
    const auto result = pcut::linked_expand(graphs, effective);
    // plan.cache()->evaluations() stays fixed after the first evaluation.
    // Graphs, embeddings, compiled models and readout mappings are shared.
}
```

`graphs.graphs()` exposes canonical abstract graphs and all connected-subgraph
maps. `graphs.embeddings()` exposes their infinite-lattice realizations with
physical vertex/edge/channel maps. Geometry does not enter graph/cache identity;
the same cache can therefore serve a different lattice. Local species are
identified by `LocalSpace::name`; dimensions, ordered charges, reference energies,
particle/parity metadata and fixed operator matrices must also match. Changing
any of these, the gap, coefficient program or external basis causes a new cache
entry. `bind` accepts one row per physical interaction and one finite real ratio per
channel, including zeros. It cannot change operator structure or geometry;
construct a new expansion and optionally share its `GraphCache` for those changes.
Bindings own their ratios and leave earlier bindings unchanged. Coupling-only
changes reuse immutable symbolic blocks. Changing channel
list order is supported by explicit channel maps. Additional leg symmetries of
particular matrices are not automatically inferred.

The charge alphabet is the union of individual **operator channel** alphabets,
including channels whose coupling is zero. Thus a coupling sweep can reuse the
same coefficients without accidentally dropping a virtual process. To construct
a finite numerical model, use `cluster_model(lattice, edges)` or `ClusterModel`
directly. `WhiteGraphExpansion::block(index, effective, basis)` evaluates the
cached symbolic block in a requested physical basis; pass an index into
`embeddings()`. Indices are shared across bindings of the same plan. Fermionic
input/output permutations include graded signs. When any channel is fermionic,
ordinary channels must preserve parity on each individual site; incompatible
mixed conventions are rejected before readout. The optional fourth argument `true` requests
only monomials touching every edge; it is an algebraic support projection, and
one-particle callers must still subtract its vacuum contribution. Raw blocks and
support projections have separate cache contexts. `structural_model(index)`
provides a compiled unbound model in physical site order, with terms in canonical
variable order, for basis and statistics inspection; use `cluster_model(graphs.bound_lattice(), embedding.edges)` for a numerical
Hamiltonian with bound couplings. `structure()` returns the shared immutable
lattice with unit couplings; `bound_lattice()` explicitly materializes a bound lattice
copy. `couplings(index)` accesses canonical channel ratios in constant time.
Connected subtraction maps are requested through `graph_subclusters(graph)`
or `embedding_subclusters(index)` and constructed lazily, preserving every
embedded occurrence.

Custom scalars use a reusable `ScalarEvaluator` with a pure callback returning
formal corrections. This example is the additive sum of all channel strengths:

```cpp
const pcut::ScalarEvaluator sum_channels(
    [](const pcut::WhiteGraph& graph, double gap, unsigned order) {
        (void)gap;
        pcut::SymbolicSeries series(order + 1);
        if (order) for (std::size_t v = 0; v < graph.variables(); ++v)
            series[1][pcut::Monomial{}.multiplied(v)] = 1.0;
        return series;
    });
const auto total = pcut::linked_scalar(clusters, 4, sum_channels, 0.0);
```

The callback must be covariant under graph relabeling and cluster additive. At
order n every monomial has total degree n; order zero is an empty polynomial,
with the physical reference per cell passed separately. Captured context must
remain immutable; construct a new evaluator when it changes. Copies share cache
identity. To evaluate a pCUT-derived scalar, capture an `EffectiveOperator` and
extract elements of `effective.symbolic_block(graph.model(gap), basis)` inside
the callback. Raw excited-sector energies and raw block traces are generally not
cluster additive. Vacuum/spectator subtraction is still required where relevant.

`build/release/pcut_white_sweep` demonstrates two unequal square-lattice coupling
colors and reports graph evaluations, embeddings and cache reuse without timing
assertions. The dimer-chain factory similarly separates its base and frustration
operators into two fixed channels; alpha is a coupling ratio.

Calculations have no configurable storage or work budgets. Choose orders and
external sectors that fit your available resources; allocation failures propagate.
Tensor encodings and dense matrix sizes are still checked for representability.
See [design.md](design.md) for implementation limits.
Only exact zeros are removed; there is no amplitude-pruning tolerance.

## Repulsive Hubbard model

For degenerate local charge zero, fermionic signs, operator-valued linking and
half-/quarter-filling examples, use the [Hubbard guide](hubbard.md). The existing
vacuum and tensor-sector examples above retain their product-vacuum contract.
