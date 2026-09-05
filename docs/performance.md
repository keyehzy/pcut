# White-graph optimization results

The optimized engine delivers a clear end-to-end coupling-sweep improvement over
`a22f963` on all four comparison workloads. The original stricter requirement of
an established cold-run improvement everywhere **and no peak-memory regression
is not met**: the Ising cold ranges overlap, and peak RSS medians increased.
The implemented trade-off accepts that overhead in exchange for substantial
complete-sweep gains, as agreed during this optimization.

## Matched end-to-end measurements

Measured on 2026-09-05, macOS 15.7.7 arm64, Apple Clang 17, Release
`-O3 -DNDEBUG`, with identical Eigen and Boost 1.90 installations and matched
compiler/architecture settings. Seven independent processes per engine/workload
ran serially, alternating engine order, after all correctness and sanitizer
workloads finished. No samples were discarded. Both engines reuse coefficient
programs and immutable topology. Current bindings also reuse graph evaluations
and sparse readout maps. Cold time includes initial model creation, coefficient
program compilation, topology, first binding and linking. Complete sweep time
includes every point, binding, linking and output formatting. Process wall time
is also retained in the raw results.

The cases and coupling sequences match the original review: 21 points for both
Ising cases and dimer, six for Hubbard. Ising and dimer include vacuum and
one-particle results. Hubbard includes the complete charge-zero operator and its
per-order summaries, with no finite assembly or diagonalization in the timed
region. Square templates have displacements `(1,0)` and `(0,1)`; four-color adds
`(1,1)` and `(1,-1)`. The harness records every binding separately.

Times below are process medians, in seconds.

| Workload | Previous cold | Optimized cold | Previous complete sweep | Optimized complete sweep | Sweep speedup |
| --- | ---: | ---: | ---: | ---: | ---: |
| Square Ising, order 4 | 0.004654 | 0.004433 | 0.065176 | 0.008141 | 8.01× |
| Four-color Ising, order 3 | 0.003459 | 0.003304 | 0.042937 | 0.009660 | 4.44× |
| Square Hubbard Q=0, order 4 | 0.325816 | 0.097005 | 1.963909 | 0.398109 | 4.93× |
| Dimer vacuum + one particle, order 6 | 0.242629 | 0.217250 | 4.425363 | 0.226763 | 19.52× |

For **every complete sweep**, the slowest optimized sample is faster than the
fastest baseline sample, so the gains exceed the observed variability. Hubbard
and dimer also have disjoint cold-time ranges (3.36× and 1.12× median speedups).
Ising cold medians improve by about 5%, but overlapping ranges do not establish
a reliable cold advantage.

## Variability and peak memory

Each time entry is `median ± MAD [minimum, maximum]`, in seconds. MAD is median
absolute deviation, not a confidence interval. RSS is process peak resident
memory in MiB, shown as `median [minimum, maximum]`. Full precision, warm-point
statistics and RSS MADs are retained in [performance-optimized.json](performance-optimized.json).

| Workload / engine | Cold time | Complete sweep | Peak RSS (MiB) |
| --- | ---: | ---: | ---: |
| square / previous | 0.004654 ± 0.000086 [0.004538, 0.011644] | 0.065176 ± 0.000714 [0.064341, 0.095378] | 2.42 [2.02, 4.55] |
| square / optimized | 0.004433 ± 0.000037 [0.004362, 0.011348] | 0.008141 ± 0.000091 [0.008050, 0.020615] | 3.28 [3.28, 3.42] |
| four_color / previous | 0.003459 ± 0.000013 [0.003384, 0.003515] | 0.042937 ± 0.000562 [0.042019, 0.043527] | 2.36 [2.36, 3.09] |
| four_color / optimized | 0.003304 ± 0.000034 [0.003242, 0.003540] | 0.009660 ± 0.000066 [0.009563, 0.010012] | 3.14 [3.14, 3.41] |
| hubbard / previous | 0.325816 ± 0.004602 [0.321214, 0.341115] | 1.963909 ± 0.028951 [1.923373, 2.023851] | 741.98 [681.09, 977.23] |
| hubbard / optimized | 0.097005 ± 0.001163 [0.095369, 0.099746] | 0.398109 ± 0.022925 [0.364629, 0.434111] | 812.64 [526.36, 898.00] |
| dimer / previous | 0.242629 ± 0.000652 [0.240611, 0.249040] | 4.425363 ± 0.014828 [4.400389, 4.458020] | 12.05 [11.88, 12.22] |
| dimer / optimized | 0.217250 ± 0.001569 [0.213660, 0.219720] | 0.226763 ± 0.001441 [0.222730, 0.229248] | 21.17 [16.11, 22.97] |

Peak memory remains a cost of this design. Dimer grows from 12.05 to 21.17 MiB
median; the small Ising cases add roughly 0.8 MiB. Hubbard's median rises from
741.98 to 812.64 MiB, with broad overlapping ranges. RSS includes dense result
matrices and allocator behavior as well as cached symbolic results; these data
do not isolate the cache's contribution. There is no claim of memory improvement.

## Reuse counts and implementation

| Workload | Previous evaluations: first / remaining sweep | Optimized evaluations: first / remaining sweep |
| --- | ---: | ---: |
| square | 118 / 2360 | 30 / 0 |
| four_color | 276 / 5520 | 13 / 0 |
| hubbard | 118 / 590 | 30 / 0 |
| dimer | 6 / 120 | 6 / 0 |

The implementation changes are general engine operations:

- `WhiteGraphExpansion::bind(lattice.couplings())` shares immutable canonical
  graphs, infinite embeddings, connected-subcluster maps, compiled models and
  sparse readout plans. Numerical ratios remain binding data.
- Infinite edge-set growth supplies embedding witnesses for canonical white
  graphs. Exact structural signatures and canonicalization results are reused
  during construction; ordered legs, species, parallel templates and channels
  remain distinct. Persistent identities contain full exact structures.
- Monomials use inline exponent packing with an unrestricted sparse fallback
  within the existing degree/index representation. Symbolic scratch uses interned
  monomial IDs and Boost flat hash tables; cached polynomials use sorted
  contiguous storage. Reached transitions are memoized per input column.
- Exact-rational row-space factorization shares linearly dependent coefficient
  continuations. The original exact program remains the cache identity. Real
  operators use real scratch arithmetic; complex operators retain complex values.
- Vacuum/one-particle and Q=0 drivers use an exact full-edge monomial projection
  to eliminate proper-subcluster subtraction. One-particle vacuum subtraction
  remains mandatory. The raw-block, general-sector and custom scalar paths remain
  available with their full semantics.
- Cached block readout traverses stored entries using persistent physical row,
  column and graded-sign maps, without dense searches through sparse maps.

See [design.md](design.md) for the support-projection argument, exact program
compilation, basis order, units, limits and cache contracts. Boost 1.81+ is now
required for its flat hash containers. No compatibility wrappers, configurable
work/storage budgets, numerical pruning or benchmark-specific branches were added.

## Correctness and reproduction

All 47 CTest cases pass after configure/build with both `release` and `sanitize`
presets. The suite retains exact coefficient fixtures, independent analytic and
literature results, intermediate-charge checks, graded embedding and linking
checks. New tests cover shared binding ownership and invalid ratios, compact
monomial fallbacks, wide complex channel alphabets, and full Hubbard linked
matrices through fourth order against independent mapped subtraction at unequal
and zero couplings. Raw and support-projected blocks have separate exact cache
contexts; cache-count tests reflect that distinction.

Every original comparison output agrees at `1e-10` absolute/relative tolerance;
the largest scaled discrepancy is 3.88e-13. Hubbard norms/traces supplement the
full-matrix correctness tests and are not treated as proof of matrix equality.
The independent dimer benchmark through eighth order also passes, with maximum
absolute vacuum-coefficient error 9.72e-17; its output is in
[benchmark.json](benchmark.json). Timing assertions are absent from CTest.

```sh
cmake --preset release
cmake --build --preset release
ctest --preset release
cmake --preset sanitize
cmake --build --preset sanitize
ctest --preset sanitize
# Run these serially after the tests finish:
python3 scripts/benchmark.py
python3 scripts/engine-comparison/run.py
```

The comparison script builds both libraries and a common harness, verifies the
clean baseline commit, and retains raw outputs, per-point values and complete
measurement JSON in `build/engine-comparison-optimized/`. See the
[tooling instructions](../scripts/engine-comparison/README.md). The compact report
contains the measured source fingerprint and build settings. Build products are
not versioned. The [earlier review](performance-review.md) remains a historical
record of the initial redesign's regressions.
