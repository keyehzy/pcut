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
    {{{{{0},0}, {{1},0}}, v}},
    1.0
};
pcut::Coefficients coefficients(pcut::charge_changes(lattice), 4);
pcut::EffectiveOperator effective(coefficients);
pcut::ClusterCatalog clusters(lattice, 4);
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

A custom `linked_scalar` callback is useful for other additive scalar quantities.
It must return `order+1` coefficients, with zero constant term; pass the reference
per cell separately. Raw excited-sector energies and raw block traces are not
usually cluster additive and must not be used without irreducible subtraction.

For coupling sweeps, reuse the coefficient table, effective operator and geometry:

```cpp
const auto topology = std::make_shared<pcut::ClusterTopology>(lattice, 4);
for (double coupling : {0.2, 0.4, 0.6}) {
    auto bound = lattice;
    bound.interactions[0].matrix = coupling * v;
    const pcut::ClusterCatalog catalog(std::move(bound), topology);
    const auto result = pcut::linked_expand(catalog, effective);
}
```

An existing catalog exposes its shared geometry with `catalog.topology()`.
Rebinding checks spatial dimension, unit-cell site count and every ordered
interaction leg, including interaction type. Matrices, gap, reference energies
and local spaces belong to each numerical binding and are validated there.
`catalog.model(edges)` shares compiled transition tables by interaction type;
use it for repeated finite cluster evaluations too. Models retain shared tables
safely after the catalog is destroyed. The free `cluster_model(lattice, edges)`
remains a direct finite-model constructor and compiles the supplied terms.
The model's actual charge alphabet must be covered by the coefficient table.
Zero-valued couplings may reduce the actual alphabet, so using the union alphabet
across a sweep avoids rebuilding coefficients. This is still a univariate series
in lambda; no polynomial interpolation or symbolic dependence is implied.

Dense blocks use `SolverOptions::max_matrix_elements`, counting all orders;
linked options also count retained weights and output series cumulatively.
Both default to 32 million complex entries (about 512 MB of payload). When
raising limits for a large calculation, raise both the linked limit and its
nested solver limit as needed. Bloch series accept an optional matrix-entry
budget. See [design.md](design.md) for accounting and scratch-memory limits.
Only exact zeros are removed; there is no amplitude-pruning tolerance.

## Repulsive Hubbard model

For degenerate local charge zero, fermionic signs, operator-valued linking and
half-/quarter-filling examples, use the [Hubbard guide](hubbard.md). The existing
vacuum and tensor-sector examples above retain their product-vacuum contract.
