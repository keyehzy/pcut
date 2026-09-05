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
