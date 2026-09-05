#!/usr/bin/env python3
"""Matched a4ea23c comparison; run serially after release/sanitizer validation."""
from pathlib import Path
import sys

SOURCE = Path(__file__).resolve().parent
sys.path.insert(0, str(SOURCE.parent))
import comparison_common as common
from comparison_workloads import CASES

BASE = 'a4ea23c88325141eed0e2f2cff387546214b976f'
KINDS = ['backends', 'isolated', 'topology', 'sweep']
STRESS = ['star8', 'cycle8', 'complete7', 'duplicates']


def main():
    args = common.arguments(__doc__, '/tmp/pcut-canonical-baseline', 'canonical-comparison')
    report = common.build(args, BASE, SOURCE, backends=True)
    report.update(source_sha256=common.source_digest(SOURCE), repeats=args.repeats, runs=[])
    for kind in KINDS:
        cases = CASES if kind in ['sweep', 'topology'] else [
            (name, 0, 0, 0) for name in [c[0] for c in CASES]+STRESS+(['star20', 'cycle64'] if kind == 'backends' else [])]
        for name, order, count, particles in cases:
            engines = ['nauty', 'traces'] if kind == 'backends' else ['previous', 'current']
            iterations = 1 if kind == 'topology' else (100 if name not in STRESS or kind == 'backends' else 1)
            def command(engine):
                if kind == 'backends':
                    return [args.output/'bin/backends', engine, name, iterations]
                return [args.output/'bin'/f'{engine}_{kind}', name,
                        *([order, count, particles] if kind == 'sweep' else [iterations])]
            parse = common.parse_sweep if kind == 'sweep' else common.parse_scalar
            for repeat, pair in common.samples(args, engines, command, parse,
                    dict(kind=kind, case=name), f'{kind}-{name}', report):
                if kind == 'sweep':
                    error = common.compare_sweeps(pair, ['embeddings', 'graphs', 'evaluations'])
                    report.setdefault('equivalence', []).append(dict(case=name, repeat=repeat, scaled_error=error))
                elif pair[engines[0]]['checksum'] != pair[engines[1]]['checksum']:
                    raise ValueError(f'{kind} {name} checksum mismatch')
    report['summary'] = []
    for kind in KINDS:
        for case in sorted({r['case'] for r in report['runs'] if r['kind'] == kind}):
            rows = [r for r in report['runs'] if r['kind'] == kind and r['case'] == case]
            summary = dict(kind=kind, case=case)
            for engine in sorted({r['engine'] for r in rows}):
                selected = [r for r in rows if r['engine'] == engine]
                keys = ['peak_mib']+(['topology_seconds', 'cold_seconds', 'sweep_seconds'] if kind == 'sweep' else ['seconds'])
                summary[engine] = common.summarize(selected, keys)
                if kind != 'sweep':
                    summary[engine]['calls'] = selected[0]['calls']
            report['summary'].append(summary)
    common.write_json(args.output/'results.json', report)
    # Preserve scalar samples in the versioned report; full physics stays in build/.
    compact = {k: v for k, v in report.items() if k != 'runs'}
    compact['runs'] = [{k: v for k, v in r.items() if k != 'steps'} for r in report['runs']]
    common.write_json(common.ROOT/'docs/performance-canonical.json', compact)


if __name__ == '__main__':
    main()
