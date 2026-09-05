# White-graph performance review (before optimization)

This is the historical review of the initial redesign. The implemented changes
and new measurements are in [performance.md](performance.md).

The current redesign is **not an end-to-end performance improvement** over
commit `a22f963`. It adds reusable symbolic evaluations and reduces graph solves,
but all four measured workflows are slower, including warm coupling sweeps.
The reduction in evaluation counts is real; treating it as evidence of reduced
overall runtime would be incorrect.

## Matched measurements

Both versions were built with Apple Clang 17, `Release`, `-O3 -DNDEBUG`, arm64,
and identical Eigen/Boost installations. Each workload ran in three independent
processes per engine, serially with alternating engine order; no sanitizer or
other pcut benchmark ran concurrently. Both engines reuse their coefficient
program. The previous engine also reuses its existing `ClusterTopology`, while
the current engine reuses `GraphCache`, matching their respective supported
coupling-sweep APIs. Every point changes coupling ratios. Initial setup and
binding are included in cold timings. Warm timings include binding/generation
and linking for each later point. Tables show medians across processes.

| Workload | Previous cold | Current cold | Cold slowdown | Previous warm/point | Current warm/point | Warm slowdown |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Square Ising, order 4 | 0.004621 s | 0.022487 s | 4.87× | 0.003077 s | 0.016402 s | 5.33× |
| Four-color Ising, order 3 | 0.003419 s | 0.007795 s | 2.28× | 0.001993 s | 0.007138 s | 3.58× |
| Square Hubbard Q=0, order 4 | 0.326762 s | 0.864275 s | 2.64× | 0.339028 s | 0.695618 s | 2.05× |
| Dimer vacuum + one particle, order 6 | 0.250987 s | 26.689863 s | 106.34× | 0.209428 s | 0.575657 s | 2.75× |

There were 21 coupling points per Ising/dimer process and six per Hubbard
process. The complete 21-point dimer sweep increased from 4.45 s to 38.24 s;
its median process peak RSS increased from 12.1 MiB to 128.7 MiB. Peak RSS medians
were higher for all four new-engine workloads, although Hubbard memory samples
varied substantially and their ranges overlap. These are local measurements,
not a universal claim about every possible model/order or machine.

The square models use two directional templates. The four-color Ising model
uses displacements `(1,0)`, `(0,1)`, `(1,1)`, `(1,-1)` with distinct numerical
ratios. Ising and dimer runs request vacuum and one-particle results. Hubbard
runs request the entire charge-zero operator and do not include finite assembly
or diagonalization in the timed region.

Energy/Bloch series and Hubbard per-order operator norms/traces agree within
`1e-10` relative/absolute tolerance at all measured points. The largest scaled
difference in these compared outputs is `7.8e-13`. This comparison supplements
the existing correctness suite; norms/traces alone are not a complete matrix
equality test.

## What improved, and why the total regressed

- Canonical evaluation reuse works. At order four the square models reduce
  118 physical cluster evaluations to 30 canonical evaluations; the four-color
  order-three example reduces 276 to 13. Later coupling points perform zero
  new symbolic graph evaluations.
- Cached dimer **linking alone** improves from about 0.209 s to 0.019 s per point.
  However, reconstructing `WhiteGraphExpansion` costs about 0.556 s at every
  coupling point. The previous topology binding takes about 0.000016 s.
  Consequently, the usable warm end-to-end API is still slower.
- The chain has six canonical graphs and six embeddings through order six, so
  graph isomorphism provides no reduction for a single chain calculation. Its
  first symbolic solve pays the cost of independent edge/channel polynomials,
  while the previous engine applies numerically combined bond matrices.
- Source inspection shows additional repeated work: the expansion constructor
  regenerates/canonicalizes graphs and embeddings for each sweep point;
  `WhiteGraphExpansion::block` recompiles physical and canonical models before
  cache lookup; built-in drivers construct another physical model per embedding;
  and each cache lookup serializes the full structural/program identity.
  The monomial hot path copies a tree-based exponent map on every multiplication.
  These are identified optimization opportunities, not individually isolated
  profiler measurements.
- The cache stores raw graph blocks. Built-in linking still substitutes,
  subtracts subclusters, and performs matrix embedding for every physical
  realization. Zero graph-evaluation misses do not make this work free.

See [white_graph.cpp](../src/white_graph.cpp),
[effective.cpp](../src/effective.cpp), and [polynomial.cpp](../src/polynomial.cpp).

## Sampled bottlenecks

Additional macOS `sample` profiles used five-second windows with a one-millisecond
sampling interval on the Release binaries. These percentages describe sampled
main-thread stacks in those windows, not exact whole-run phase fractions.

- **Cold order-six dimer:** almost the entire window is inside symbolic block
  evaluation. About 48% of leaf samples are allocation/free machinery; another
  27% are in the symbolic traversal body. Monomial multiplication copies a
  tree-based exponent map, and nested state/polynomial maps add allocation,
  comparison, and destruction costs. Local operator application is much smaller
  in this sample. This model has no reduction in graph count to offset its
  additional formal-channel work.
- **Repeated order-six dimer graph construction:** about 89% of leaf samples
  are string `push_back` or its call stub. Graph refinement and channel ordering
  repeatedly serialize full dense operator matrices into structural signatures.
  The measured constructor bottleneck is therefore repeated serialization;
  these samples do not establish permutation search as the dominant cost.
  Construction consumes about 0.556 s of each 0.576 s warm point (97%).
- **Warm order-four Hubbard:** about 64% of leaf samples are in
  `WhiteGraphExpansion::block`, whose dense row/column loop searches the cached
  sparse map for every matrix position, including absent entries. Dense matrix
  initialization accounts for another 5%, and matrix embedding for about 5%.
  Structural string construction accounts for about 11% across the sampled
  warm points. Recompilation is present but is not the dominant cost here;
  the cache lookup subtree itself accounts for less than 1% of samples.

The immediate targets are reusable graph/embedding plans and precomputed
structural signatures, a compact monomial representation with fewer temporary
allocations, and cached-block readout that iterates existing sparse entries.
Their actual gains must be measured after implementation. Raw profiles are
`build/engine-comparison/cold-dimer-sample.txt`, `graph-setup-sample.txt`, and
`warm-hubbard-sample.txt` in the same directory.

## Consequence

Do not describe this revision as a demonstrated speed or memory optimization.
Before doing so, retain reusable graph/embedding plans across coupling bindings,
reuse compiled models and state/variable maps, and reduce symbolic allocation
and unnecessary raw-block work. Verify those changes with matched cold runs,
complete coupling sweeps, and memory measurements. Graph evaluation counts should
remain a structural regression test, alongside those end-to-end measurements.

The numerical summary and ranges are retained in
[performance-comparison.json](performance-comparison.json). The session's
comparison harness, CMake project, raw outputs and full JSON results are in
`build/engine-comparison/`. The previous implementation was built in an isolated
detached worktree at `/tmp/pcut-before-white-graphs`; the current library was not
changed during this review. The harness can be rerun with
`python3 build/engine-comparison/run.py`, followed by
`python3 build/engine-comparison/summarize.py`.
