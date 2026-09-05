#!/usr/bin/env python3
"""Build matched engines and measure serial alternating full coupling sweeps."""
import argparse
import hashlib
import json
import math
import pathlib
import platform
import statistics
import subprocess
import time

SOURCE = pathlib.Path(__file__).resolve().parent
ROOT = SOURCE.parent.parent
CASES = [('square', 4, 21, 1), ('four_color', 3, 21, 1),
         ('hubbard', 4, 6, 0), ('dimer', 6, 21, 1)]


def command(args, cwd=ROOT):
    subprocess.run(list(map(str, args)), cwd=cwd, check=True)


def cmake_cache(path):
    result = {}
    for line in path.read_text().splitlines():
        if not line.startswith(('#', '//')) and ':' in line and '=' in line:
            key, value = line.split('=', 1)
            result[key.split(':')[0]] = value
    return result


def describe(values):
    center = statistics.median(values)
    return dict(median=center, minimum=min(values), maximum=max(values),
                mad=statistics.median(abs(x-center) for x in values))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', type=pathlib.Path, default=pathlib.Path('/tmp/pcut-before-white-graphs'))
    parser.add_argument('--output', type=pathlib.Path, default=ROOT/'build/engine-comparison-optimized')
    parser.add_argument('--repeats', type=int, default=7)
    args = parser.parse_args()
    if args.repeats < 3:
        parser.error('at least three independent repetitions are required')
    args.output = args.output.resolve()
    args.baseline = args.baseline.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    if not args.baseline.exists():
        command(['git', 'worktree', 'add', '--detach', args.baseline, 'a22f963'])
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=args.baseline, text=True).strip()
    expected = subprocess.check_output(['git', 'rev-parse', 'a22f963'], cwd=ROOT, text=True).strip()
    if head != expected or subprocess.check_output(['git', 'status', '--porcelain'], cwd=args.baseline, text=True):
        raise RuntimeError('baseline must be a clean checkout of a22f963')
    command(['cmake', '--preset', 'release'])
    cache = cmake_cache(ROOT/'build/release/CMakeCache.txt')
    keys = ['CMAKE_CXX_COMPILER', 'CMAKE_CXX_FLAGS', 'CMAKE_CXX_FLAGS_RELEASE',
            'CMAKE_OSX_ARCHITECTURES', 'CMAKE_OSX_DEPLOYMENT_TARGET', 'Eigen3_DIR', 'Boost_DIR']
    matched = {key: cache[key] for key in keys if key in cache}
    definitions = [f'-D{k}={v}' for k, v in matched.items()]
    command(['cmake', '--preset', 'release', *definitions], cwd=args.baseline)
    command(['cmake', '--build', '--preset', 'release', '-j', '6'])
    command(['cmake', '--build', '--preset', 'release', '-j', '6'], cwd=args.baseline)
    command(['cmake', '-S', SOURCE, '-B', args.output/'bin', '-G', 'Ninja',
             '-DCMAKE_BUILD_TYPE=Release', f'-Dprevious_ROOT={args.baseline}',
             f'-Dcurrent_ROOT={ROOT}', *definitions])
    command(['cmake', '--build', args.output/'bin', '-j', '6'])
    digest = hashlib.sha256()
    for path in sorted([*(ROOT/'src').glob('*.cpp'), *(ROOT/'include/pcut').glob('*.hpp'), SOURCE/'measure.cpp']):
        digest.update(str(path.relative_to(ROOT)).encode())
        digest.update(path.read_bytes())
    report = dict(baseline_commit=expected, current_source_sha256=digest.hexdigest(), platform=platform.platform(), build=matched,
                  compiler=subprocess.check_output([cache['CMAKE_CXX_COMPILER'], '--version'], text=True),
                  method='Serial alternating independent processes. Both reuse programs and immutable topology. '
                         'Cold includes initial model, coefficient program, topology, binding, linking and summaries. '
                         'Sweep includes every point and output formatting. No concurrent pcut workloads.', runs=[])
    for name, order, count, particles in CASES:
        for repeat in range(args.repeats):
            paired = {}
            for engine in (['previous', 'current'] if repeat % 2 == 0 else ['current', 'previous']):
                start = time.perf_counter()
                output = subprocess.check_output([str(args.output/'bin'/engine), name,
                                                  str(order), str(count), str(particles)], text=True)
                (args.output/f'{name}-{engine}-{repeat}.txt').write_text(output)
                row = dict(model=name, engine=engine, repeat=repeat, order=order, points=count,
                           process_wall_seconds=time.perf_counter()-start, steps=[])
                for line in output.splitlines():
                    f = line.split()
                    if f[0] == 'setup':
                        row['coefficient_seconds'], row['topology_seconds'] = map(float, f[1:])
                    elif f[0] == 'cold':
                        row['cold_seconds'] = float(f[1])
                    elif f[0] == 'step':
                        row['steps'].append(dict(index=int(f[1]), bind_seconds=float(f[2]),
                            link_seconds=float(f[3]), graphs=int(f[4]), embeddings=int(f[5]),
                            evaluations=int(f[6]), values=list(map(float, f[7:]))))
                    elif f[0] == 'total':
                        row['sweep_seconds'] = float(f[1])
                        row['peak_mib'] = int(f[2])/(2**20 if platform.system() == 'Darwin' else 1024)
                row['warm_seconds'] = statistics.median(s['bind_seconds']+s['link_seconds'] for s in row['steps'][1:])
                report['runs'].append(row)
                paired[engine] = row
                print(name, engine, repeat, 'cold', row['cold_seconds'], 'sweep', row['sweep_seconds'],
                      'MiB', row['peak_mib'], flush=True)
            error = 0.0
            for old, new in zip(paired['previous']['steps'], paired['current']['steps'], strict=True):
                assert old['embeddings'] == new['embeddings']
                for x, y in zip(old['values'], new['values'], strict=True):
                    assert math.isclose(x, y, rel_tol=1e-10, abs_tol=1e-10), (name, x, y)
                    error = max(error, abs(x-y)/max(1, abs(x), abs(y)))
            report.setdefault('equivalence', []).append(dict(model=name, repeat=repeat, scaled_error=error))
            (args.output/'results.json').write_text(json.dumps(report, indent=2)+'\n')
    summary = {k: report[k] for k in ['baseline_commit', 'current_source_sha256', 'platform', 'build', 'compiler', 'method']}
    summary['repeats'] = args.repeats
    summary['cases'] = []
    for name, order, count, _ in CASES:
        case = dict(model=name, order=order, points=count)
        for engine in ['previous', 'current']:
            rows = [r for r in report['runs'] if r['model'] == name and r['engine'] == engine]
            case[engine] = {key: describe([r[key] for r in rows]) for key in
                            ['cold_seconds', 'warm_seconds', 'sweep_seconds', 'peak_mib']}
            case[engine]['first_evaluations'] = rows[0]['steps'][0]['evaluations']
            case[engine]['later_evaluations'] = sum(s['evaluations'] for s in rows[0]['steps'][1:])
        case['scaled_output_error'] = max(r['scaled_error'] for r in report['equivalence'] if r['model'] == name)
        case['sweep_speedup'] = case['previous']['sweep_seconds']['median']/case['current']['sweep_seconds']['median']
        case['cold_speedup'] = case['previous']['cold_seconds']['median']/case['current']['cold_seconds']['median']
        summary['cases'].append(case)
    (args.output/'summary.json').write_text(json.dumps(summary, indent=2)+'\n')


if __name__ == '__main__':
    main()
