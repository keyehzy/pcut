# White-graph canonicalizer measurements

The production exhaustive vertex search has been replaced by sparse **nauty
2.9.3**. The production encoding retains every edge/channel occurrence, full exact
physical signatures, ordered legs and explicit maps. Canonical order may change;
full structural and cache identities remain exact. See [design.md](design.md)
and [dependency provenance](dependencies.md).

These measurements compare the starting revision
`a4ea23c88325141eed0e2f2cff387546214b976f` with the measured implementation. The earlier
engine redesign comparison against `a22f963` is preserved separately in
[performance-engine.md](performance-engine.md).

The retained timings below predate removal of redundant channel leaves. The
current duplicate fixture has 98 incidence nodes instead of 354 and still
retains all 256 channel occurrences. Production and backend tooling now share
encoding, partition construction, sparse storage and nauty cleanup. The tables
below retain the original seven-repeat measurements.

## Method and scope

Measured 2026-09-05 on macOS 15.7.7 arm64, Apple Clang 17, Release
`-O3 -DNDEBUG`, with matching compiler/architecture/Eigen/Boost settings. Seven
independent processes per engine/workload ran serially with alternating order,
after release and sanitizer workloads finished. Every sample from the final
run is included. The [machine-readable report](performance-canonical.json)
contains the source fingerprint, compiler/build metadata, all raw scalar samples,
medians, MADs and ranges. No timing assertions or work budgets were introduced.

The real catalogs contain 118 square-Ising inputs (order 4), 276 four-color-Ising
inputs (order 3), 118 square-Hubbard inputs (order 4), and six dimer inputs (order
6). Inputs come from identical sorted physical embeddings; geometry and ratios
are removed. Each catalog is canonicalized 100 times per isolated process.
Backend timing uses the same prepared incidence graphs and includes partition
reset, canonical output, adjacency readout and scratch cleanup; it excludes
physical validation/serialization, encoding, maps and the exact group quotient.
Complete API timing includes all of those operations.

Standalone topology processes construct one catalog and sample peak RSS before
coefficient generation or linking. End-to-end sweeps use the original four
comparison workloads: 21 points for both Ising models and dimer vacuum/particle,
six points for the full Hubbard Q=0 operator. Cold time includes model,
coefficients, topology, first binding and linking. Sweep time includes all points
and output formatting. Peak RSS is whole-process memory, not incremental nauty
scratch; backend/API RSS includes untimed catalog/input setup. It includes
allocator and dense result storage for end-to-end runs.

## Backend choice

Median microseconds per prepared incidence graph; RSS is MiB.

| Workload | nauty µs | Traces µs | nauty / Traces peak RSS |
| --- | ---: | ---: | ---: |
| Square Ising, order 4 | 1.15 | 2.31 | 2.97 / 3.12 |
| Four-color Ising, order 3 | 1.00 | 2.16 | 3.23 / 3.31 |
| Square Hubbard, order 4 | 1.15 | 2.28 | 7.89 / 8.05 |
| Dimer, order 6 | 1.05 | 2.37 | 3.36 / 3.48 |
| Star, 8 leaves | 7.36 | 2.62 | 1.64 / 1.78 |
| Directed 8-cycle | 2.06 | 2.60 | 1.62 / 1.78 |
| Complete bidirected graph, 7 vertices | 32.66 | 59.35 | 1.69 / 1.92 |
| 32 duplicate edges × 8 duplicate channels | 346.02 | 11.80 | 2.02 / 2.23 |
| Star, 20 leaves | 54.45 | 5.31 | 1.67 / 1.84 |
| Directed 64-cycle | 14.61 | 13.41 | 1.80 / 2.03 |

Sparse nauty is about 2–2.3× faster on the four real catalogs. Traces is faster
on the stars, larger cycle and duplicate stress case; it is **not** uniformly
slower. We select nauty because of the catalog results and its documented
integer `userlevelproc` stabilizer indices, which support exact group order
without another permutation-group implementation. Traces' documented group-size
statistics are approximate and are insufficient for exact overflow checks.
Both backends are tested on identical colored inputs; their canonical ordering
need not agree. Only nauty is linked into production.

## Complete canonicalization and topology

The complete API includes validation, exact signatures, encoding, labeling,
vertex automorphism count, structural serialization and every map. Medians are
microseconds per call. Speedup is baseline/current; below 1 means a regression.

| Workload | Exhaustive µs | nauty API µs | Speedup | Previous / current RSS (MiB) |
| --- | ---: | ---: | ---: | ---: |
| Square Ising, order 4 | 11.15 | 9.11 | 1.22× | 2.92 / 3.16 |
| Four-color Ising, order 3 | 8.04 | 7.28 | 1.10× | 3.23 / 3.41 |
| Square Hubbard, order 4 | 20.58 | 18.37 | 1.12× | 8.08 / 8.09 |
| Dimer, order 6 | 39.93 | 41.21 | 0.97× | 3.12 / 3.34 |
| Star, 8 leaves | 103351.71 | 88.71 | 1165.08× | 1.92 / 1.83 |
| Directed 8-cycle | 95250.38 | 79.96 | 1191.26× | 1.88 / 1.84 |
| Complete bidirected graph, 7 vertices | 57049.75 | 203.12 | 280.86× | 1.78 / 1.94 |
| 32 duplicate edges × 8 duplicate channels | 233.04 | 653.83 | 0.36× | 2.31 / 2.55 |

Standalone topology construction measures time and process peak memory without
coefficient generation or linking. Times are milliseconds.

| Workload | Previous ms | Current ms | Previous / current peak RSS (MiB) |
| --- | ---: | ---: | ---: |
| Square Ising, order 4 | 1.979 | 1.689 | 2.75 / 2.88 |
| Four-color Ising, order 3 | 1.358 | 1.282 | 2.66 / 2.80 |
| Square Hubbard, order 4 | 2.634 | 2.419 | 7.69 / 7.70 |
| Dimer, order 6 | 0.315 | 0.343 | 3.06 / 3.27 |

The star/cycle/complete cases demonstrate removal of the exhaustive factorial
search, with roughly 280–1,200× complete-API improvements at these measured
sizes. No baseline timing is claimed for the backend-only star-20 or cycle-64.
Dimer API time regresses by 3.2%, and its standalone topology median by 8.6%
(about 27 µs). The old refinement already makes these small oriented paths
trivial; incidence construction and full gadget labeling add fixed work.
Duplicate-heavy input was 2.8× slower than the baseline: two physical
vertices are trivial for the old search, while that encoding processed 354 incidence
vertices and a `32!` edge kernel. The current encoding has 98 nodes; the exact
edge kernel is unchanged.

## Regression investigation

The initial encoding used full operator-signature colors on channel leaves.
It duplicated matrix serialization and introduced a channel-permutation kernel
`(8!)^32` in the duplicate stress case. The initial complete-API medians were
14.53 ms for duplicates, 27.08 µs for Hubbard and 58.49 µs for dimer. The initial
seven-repeat summary and source fingerprint are retained in
[performance-canonical-initial.json](performance-canonical-initial.json).

The measured encoding used distinct sorted-channel slot colors. The parent edge
still contains the complete sorted operator list, so each slot identifies its
operator exactly. All channel occurrences and maps survive, but internal
channel permutations disappear. Every physical vertex automorphism still
extends through corresponding slots; only duplicate-edge permutations remain
in the exact kernel quotient. This removes repeated operator color bytes and
most of the duplicate stress cost without changing graph identity semantics,
lattice growth, symbolic evaluation or embedding multiplicities. Final oracle,
map, channel-cancellation, cache-reuse and graded-readout tests pass.

## End-to-end trade-offs

Times are milliseconds; each workload uses the same complete coupling sequence
on both engines. These topology times follow coefficient compilation, unlike
the standalone topology processes above.

| Workload | Topology previous / current | Cold previous / current | Sweep previous / current |
| --- | ---: | ---: | ---: |
| Square Ising, order 4 | 1.470 / 1.275 | 3.261 / 3.117 | 6.352 / 6.217 |
| Four-color Ising, order 3 | 1.406 / 1.301 | 2.526 / 2.398 | 7.912 / 7.796 |
| Square Hubbard, order 4 | 2.674 / 2.436 | 93.970 / 92.724 | 366.433 / 363.914 |
| Dimer, order 6 | 0.296 / 0.346 | 224.377 / 225.481 | 233.536 / 234.685 |

All four cold-time and sweep-time ranges overlap between engines. The small
median differences do **not** establish a reliable broad end-to-end gain or
regression. Canonicalization is a small part of these complete workloads;
coefficient compilation and symbolic/readout work remain substantial.

Peak process RSS, median [minimum, maximum], in MiB:

| Workload | Previous | Current |
| --- | ---: | ---: |
| Square Ising, order 4 | 2.97 [2.97, 3.23] | 3.05 [3.05, 3.58] |
| Four-color Ising, order 3 | 2.91 [2.91, 3.34] | 2.98 [2.98, 3.41] |
| Square Hubbard, order 4 | 777.12 [647.31, 852.05] | 766.41 [666.11, 829.47] |
| Dimer, order 6 | 18.31 [16.08, 21.92] | 19.89 [16.19, 22.80] |

There is no demonstrated end-to-end memory improvement. Small Ising medians
increase by 0.078 MiB, and dimer by 1.58 MiB. Hubbard's median decreases, but
its broad overlapping range does not establish a reduction. Standalone
topology and API RSS also include input storage and allocator behavior; these
measurements do not identify the contribution of each allocation. Inspection
confirms nauty scratch is freed after each call and cache contexts/evaluation
counts are unchanged; no extra persistent canonicalizer cache was added.

## Validation and reproduction

Release and sanitize configure/build/CTest workflows each pass all 53 tests.
Additional checks cover the exhaustive oracle on small and real catalog graphs,
arbitrary vertex/edge/channel reorderings, connected regular graphs with equal
refinement colors, all three maps, exact species/operator metadata, unequal
local dimensions, on-site/higher-arity edges, duplicate templates/channels,
vertex groups through `20!`, overflow at `21!`, huge gadget kernels, concurrent
calls, and recovery from every injected upstream allocation failure in the
fixture. Graded fermionic blocks, linked cancellation, subcluster multiplicity,
coefficient fixtures and literature benchmarks remain covered. The strengthened
identical-channel test additionally passes under both presets: opposite ratios
cancel and channel reordering reuses cached full blocks.

All final sweep outputs agree at `1e-10` absolute/relative tolerance; maximum
scaled discrepancy is **4.14e-16**. Graph, embedding and evaluation counts match
exactly. First/later graph evaluations remain 30/0 for square and Hubbard,
13/0 for four-color, and 6/0 for dimer. Scalar norms/traces in the timing harness
supplement, rather than replace, full-matrix correctness tests.

Default checksum-verified dependency fetching, offline source override,
relocated static-package consumption and shared-package consumption were
validated. Downloaded research reference checksums still match. Remaining
limitations are native GCC/Clang/POSIX-sh integration (validated here with Apple
Clang), untested cross/MSVC builds, nauty's integer representation limits,
expensive hard graph families/duplicate-edge kernels, and the pre-existing
exponential symbolic/cluster growth. No ordering compatibility, finite-size
approximation, graph-work budget, or truncation was added.

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
cmake --preset sanitize
cmake --build --preset sanitize
ctest --preset sanitize
# Run serially after validation finishes:
python3 scripts/canonical-comparison/run.py
```

See the [comparison tooling instructions](../scripts/canonical-comparison/README.md)
for workload details, scope and raw outputs under `build/canonical-comparison/`.
All scalar samples and variability are retained in the versioned JSON report;
build products remain unversioned.
