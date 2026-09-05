# Dependencies and canonical-labeling provenance

pcut uses Eigen 3.4+ for linear algebra, Boost 1.81+ for exact arithmetic and
containers, and Catch2 3 for tests. Catch2 falls back to its checksummed 3.8.1
archive. Eigen and Boost remain system dependencies.

## nauty 2.9.3

The private canonical-labeling dependency is the official **nauty/Traces 2.9.3**
source release, dated January 1, 2026. Production uses sparse nauty, with
64-bit setwords, dynamically sized graphs and thread-local scratch. Traces 2.2
is built only by the comparison tooling.

- Archive: [nauty2_9_3.tar.gz](https://pallini.di.uniroma1.it/nauty2_9_3.tar.gz)
- SHA-256: `9fc4edae04f88a0f5883985be3b39cf7f898fd6cc96e96b9ee25452743cc1b5b`
- [Official introduction](https://pallini.di.uniroma1.it/Introduction.html):
  colored isomorphism, canonical labeling, individualization/refinement and
  stabilizer-chain group order.
- [Official guide](https://pallini.di.uniroma1.it/nug29.pdf), sections 3–9:
  sparse graphs, ordered initial partitions, `lab` (canonical-to-input), and
  `userlevelproc`; section 16: build configuration and TLS.
  The current guide identifies itself as version 2.9.0, July 17, 2025.
  The separately downloaded guide has SHA-256
  `6123215d9caec8c2d0427583fc4f329934983e86883a2ea4aa15b7014cafda58`.
- Brendan D. McKay and Adolfo Piperno, *Practical Graph Isomorphism, II*,
  Journal of Symbolic Computation 60 (2014), 94–112, linked by the official site.

The introduction, guide and archive were downloaded into `/tmp` for inspection
on 2026-09-05. The inspected files include `nauty.c`, `nausparse.c`, `nauty-h.in`,
`traces.h`, `traces.c`, `configure.ac` and the upstream makefile. Downloaded
research references and dependency sources are not edited.

`cmake/Nauty.cmake` uses CMake FetchContent with an exact archive checksum. It
runs upstream `configure` in the binary tree to generate platform headers,
then builds only the required C implementation files and pcut's allocation
adapter. Configuration uses the selected C compiler, TLS and `WORDSIZE=64`.
It is cached by source location/compiler/flags. The supported integration is a
native GCC/Clang C and C++ toolchain with a POSIX `sh`; MSVC and cross compilation
have not been validated. It does not require an installed nauty package,
Autoconf, or a separately invoked upstream make build.

The build leaves upstream files unchanged. A forced-include adapter selects the
generated header and overrides allocation macros. It checks allocation sizes,
throws `std::bad_alloc` on allocation failure, and allocates replacement scratch
before freeing old storage. An intrusive, thread-local allocation list recovers
partially constructed upstream nodes on exceptions. The C objects are compiled
with exception unwinding enabled. All scratch is freed at the end of each call.
Nauty's randomized Schreier option is disabled. Fault-injection tests fail every
upstream allocation reached by a symmetric graph, check for leaks, and then
canonicalize again in the same thread. Test builds instrument the allocator in
`pcut_nauty` itself, so injection also reaches a shared pcut library on macOS.
The test requires at least one injected exception and rejects an incorrect key
immediately. Address/UB sanitizers cover the C
implementation as well as the adapter and library.

The installed package includes `libpcut_nauty` and exports its private link
relationship through `pcut::pcut`. Downstream users need no nauty headers,
downloads or source directories. Upstream `COPYRIGHT` and `LICENSE-2.0.txt`
are installed under `share/pcut/nauty`; the selected files are under Apache 2.0.
A relocated installation is checked with `tests/consumer`.

For an offline build, extract the **verified pinned archive** locally and use
`-DFETCHCONTENT_SOURCE_DIR_NAUTY=/absolute/path/nauty2_9_3`. Install Catch2 3,
supply its FetchContent source directory, or disable tests. Source-directory
overrides bypass FetchContent's download checksum, so verify the archive first.
No downloaded dependency or generated build product is versioned in pcut.
