# End-to-end white-graph comparison

Run from the repository root with Python 3.10+:

```sh
python3 scripts/engine-comparison/run.py
```

The script creates a detached `a22f963` worktree if needed, rejects a dirty or
incorrect baseline, and configures/builds both libraries with matching Release
compiler, flags, architecture, Eigen and Boost settings. It then builds the same
measurement source against each API. `--baseline`, `--output` and `--repeats`
select the worktree, artifact directory and repetition count (default seven).
These are measurement controls, not library work or storage budgets.

Run after correctness/sanitizer tests finish, with no competing benchmarks or
builds. Each workload uses serial independent processes, alternating which engine
runs first. No result is discarded. The four workloads, orders, point counts and
coupling sequences match the original performance review. Both engines reuse
coefficient programs and immutable topology. Current bindings additionally share
canonical evaluations and sparse readout maps. Cold timing includes all initial
setup and first-point work; sweep timing includes every binding/link and output
formatting. Program setup, binding and linking are also recorded separately.

`results.json` and individual text outputs retain all measurements and numerical
summaries. `summary.json` reports median, minimum, maximum and median absolute
deviation (MAD) for cold time, warm-point time, total sweep and process peak RSS.
macOS reports `ru_maxrss` in bytes; Linux reports KiB. The script checks every
compared numerical summary at `1e-10` absolute/relative tolerance. Hubbard
norms/traces are supplementary checks: the correctness suite independently
compares full linked matrices through fourth order at unequal and zero couplings.
Timing assertions are deliberately absent from CTest.

Optional macOS profiling, separate from timing runs:

```sh
sample PID 5 1 -file build/engine-comparison-optimized/profile.txt
```

Keep build products, raw process logs and large measurement files in `build/`.
The report in `docs/performance.md` links the retained compact summary and explains
any cold-time or memory overhead alongside the sweep benefit.
