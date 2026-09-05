#!/usr/bin/env python3
"""Build matched engines and measure serial alternating full coupling sweeps."""
from pathlib import Path
import sys

SOURCE = Path(__file__).resolve().parent
sys.path.insert(0, str(SOURCE.parent))
import comparison_common as common
from comparison_workloads import CASES

BASE = 'a22f963'


def main():
    args = common.arguments(__doc__, '/tmp/pcut-before-white-graphs', 'engine-comparison-optimized')
    report = common.build(args, BASE, SOURCE)
    report.update(current_source_sha256=common.source_digest(SOURCE),
                  method='Serial alternating independent processes. Both reuse programs and immutable topology. '
                         'Cold includes initial model, coefficient program, topology, binding, linking and summaries. '
                         'Sweep includes every point and output formatting. No concurrent pcut workloads.', runs=[])
    for name, order, count, particles in CASES:
        def command(engine):
            return [args.output/'bin'/engine, name, order, count, particles]
        for repeat, pair in common.samples(args, ['previous', 'current'], command, common.parse_sweep,
                dict(model=name, order=order, points=count), name, report):
            error = common.compare_sweeps(pair, ['embeddings'])
            report.setdefault('equivalence', []).append(dict(model=name, repeat=repeat, scaled_error=error))
            common.write_json(args.output/'results.json', report)
    summary = {k: report[k] for k in ['baseline_commit', 'current_source_sha256', 'platform', 'build', 'compiler', 'method']}
    summary.update(repeats=args.repeats, cases=[])
    for name, order, count, _ in CASES:
        case = dict(model=name, order=order, points=count)
        for engine in ['previous', 'current']:
            rows = [r for r in report['runs'] if r['model'] == name and r['engine'] == engine]
            case[engine] = common.summarize(rows, ['cold_seconds', 'warm_seconds', 'sweep_seconds', 'peak_mib'])
            case[engine]['first_evaluations'] = rows[0]['steps'][0]['evaluations']
            case[engine]['later_evaluations'] = sum(s['evaluations'] for s in rows[0]['steps'][1:])
        case['scaled_output_error'] = max(r['scaled_error'] for r in report['equivalence'] if r['model'] == name)
        for metric in ['sweep', 'cold']:
            case[f'{metric}_speedup'] = case['previous'][f'{metric}_seconds']['median']/case['current'][f'{metric}_seconds']['median']
        summary['cases'].append(case)
    common.write_json(args.output/'summary.json', summary)


if __name__ == '__main__':
    main()
