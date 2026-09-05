# Research sources

Downloaded from arXiv on 2026-09-05. Original archives are in `sources/*.tar`
(the arXiv response may actually be gzip); extracted sources are alongside them.
SHA-256 hashes are in `checksums.json`. These are third-party research documents,
retaining their authors' rights; they are not part of the library's own code.

| arXiv | Authors and title | Source URL | Main LaTeX |
| --- | --- | --- | --- |
| 1505.02975v1 | Coester and Schmidt, Optimizing linked cluster expansions by white graphs | https://arxiv.org/src/1505.02975 | sources/1505.02975/white_graphs.tex |
| cond-mat/9906243v1 | Knetter and Uhrig, Perturbation Theory by Flow Equations: Dimerized and Frustrated S=1/2 Chain | https://arxiv.org/src/cond-mat/9906243 | sources/cond-mat-9906243/kuvers5.tex |
| cond-mat/0306333v1 | Knetter, Schmidt and Uhrig, The Structure of Operators in Effective Particle-Conserving Models | https://arxiv.org/src/cond-mat/0306333 | sources/cond-mat-0306333/perstruc-gu2a.tex |

The 1999 source supplies the sign-generator recurrence (`dgl`), the effective
Hamiltonian (`effhamilfin`), exact coefficient tables, the ground-state series
(`E_grund`), and hopping coefficients (Appendix D). The 2003 source supplies
particle-irreducible subtraction (`rekurs`) and the cluster-additivity argument.
The 2015 source provides LCE context; its white-graph optimization is deferred.

## Hubbard sources

Downloaded from the following primary arXiv source endpoints on 2026-09-05.
The exact response bytes and extracted text are retained unmodified; the hashes
in `checksums.json` identify these snapshots of the unversioned endpoints.
Existing reference archives and extracted sources were not changed.

| Authors and source | Download URL | Original response / extracted LaTeX |
| --- | --- | --- |
| Chernyshev, Galanakis, Phillips, Rozhkov, Tremblay, *Higher order effective low-energy theories*; published as *Higher order corrections to effective low-energy theories for strongly correlated electron systems*, Phys. Rev. B 70, 235111 (2004), DOI 10.1103/PhysRevB.70.235111 | https://arxiv.org/src/cond-mat/0407255 | `sources/cond-mat-0407255.tar` (gzip of one LaTeX file, not a tar) / `sources/cond-mat-0407255/main.tex` |
| Delannoy, Gingras, Holdsworth, Tremblay, *Néel order, ring exchange, and charge fluctuations in the half-filled Hubbard model*, Phys. Rev. B 72, 115114 (2005), DOI 10.1103/PhysRevB.72.115114 | https://arxiv.org/src/cond-mat/0412033 | `sources/cond-mat-0412033.tar` (gzip tar) / `sources/cond-mat-0412033/dght.tex` |

Chernyshev's labels `AM`, `S0`, and `Canonical_T_explicit` supply the complete
fourth-order doped Hamiltonian and its unitary-convention freedom. The exact
sign-generator coefficients agree with CT1 plus the documented `gamma=1/4`
rotation. Delannoy's `Hs4` supplies the half-filled operator-word expression;
`Hso` supplies the square spin Hamiltonian and `J1`, `J2`, `J3`, `Jc`. Constants
omitted in the spin-only presentation are retained in our matrices and fixed by
the zero-energy fully polarized state. No downloaded text was edited to match
our conventions. Run `python3 scripts/verify_references.py` to check both the
response hashes and extracted content.
