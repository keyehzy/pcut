# Canonicalizer comparison

Run from the repository root after correctness/sanitizer jobs have finished:

```sh
python3 scripts/canonical-comparison/run.py
```

The script checks out the clean starting revision
`a4ea23c88325141eed0e2f2cff387546214b976f` in `/tmp/pcut-canonical-baseline`,
configures both releases with matched compiler/flags/Eigen/Boost settings,
then builds a common harness. Seven fresh processes per engine and workload run
serially with alternating engine order; every sample is retained. The script
checks physics outputs at `1e-10` absolute/relative tolerance and requires
identical graph, embedding and evaluation counts. Both versions use the same
modern binding API. The old engine-comparison harness is reused without its
historical `PREVIOUS` compatibility define.

Four measurements have distinct scopes:

1. `backends` compares sparse nauty and Traces on identical **prepared** colored
   incidence graphs. It includes copying/resetting partitions, canonical output,
   an adjacency readout and scratch cleanup. It excludes physical serialization,
   validation, encoding, full maps and the exact group quotient. Each backend
   also checks invariance under a reversed vertex/edge/channel input before the
   timer. These timings alone do not establish pcut-level improvements.
2. `isolated` measures the complete public `canonicalize`: input validation,
   signatures, incidence construction (current), labeling, exact vertex group
   order, full structural serialization and maps. Catalog creation is outside
   the timer. Physical inputs are reconstructed from sorted infinite-lattice
   embeddings, so both engines get identical input ordering and matrices.
3. `topology` measures one complete topology constructor in a fresh process,
   with peak RSS sampled before any coefficient generation or linking. Model
   creation precedes the timer. This isolates topology-process peak memory.
4. `sweep` measures topology construction, complete cold time, all sweep points,
   and process peak RSS using the shared end-to-end harness. Cold includes model
   and coefficient setup, topology, binding, linking and result summaries.

RSS is **process** peak resident memory, including untimed setup/input retention
for backend and isolated runs; it is not an incremental scratch allocation
measurement. macOS bytes and Linux KiB are converted to MiB. Reports retain
median, MAD, minimum, maximum and all raw scalar samples. Timing assertions are
not CTest correctness checks.

Real catalogs are square Ising order 4 (118 physical inputs), four-color Ising
order 3 (276), square Hubbard order 4 (118), and dimer order 6 (6). Isolated and
backend runs repeat each catalog 100 times. Stress graphs are outward star with
8 leaves, directed 8-cycle, complete bidirected graph on 7 vertices, and 32
identical ordered parallel edges with 8 identical channels apiece. Complete API
stress calls run once per process; backend calls repeat 100 times. Backend-only
star with 20 leaves and directed 64-cycle also repeat 100 times. We do not run
those two larger cases with the exhaustive baseline and claim no measured
baseline speedup for them. Stress sizes are explicit benchmark workloads,
not production limits or timeouts.

The end-to-end workloads/couplings match `scripts/engine-comparison`: square
Ising order 4, four-color Ising order 3, dimer vacuum/particle order 6, each with
21 ratios, and square Hubbard Q=0 order 4 with six ratios. Physics summaries
supplement the full-matrix, graded-sign and linked-cancellation correctness tests.

Raw outputs and complete per-point results go to `build/canonical-comparison/`.
`docs/performance-canonical.json` retains build metadata, source fingerprint,
all scalar samples and summary variability. Dependency downloads and binaries
remain in build directories or `/tmp` and are not versioned.
