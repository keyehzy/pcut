"""Shared build, sampling, validation and reporting for historical comparisons."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import platform
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parent.parent


def command(args, cwd=ROOT):
    subprocess.run(list(map(str, args)), cwd=cwd, check=True)


def arguments(description, baseline, output):
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument('--baseline', type=Path, default=Path(baseline))
    parser.add_argument('--output', type=Path, default=ROOT/'build'/output)
    parser.add_argument('--repeats', type=int, default=7)
    args = parser.parse_args()
    if args.repeats < 3:
        parser.error('at least three independent repetitions are required')
    args.output = args.output.resolve()
    args.baseline = args.baseline.resolve()
    args.output.mkdir(parents=True, exist_ok=True)
    return args


def cmake_cache(path):
    result = {}
    for line in path.read_text().splitlines():
        if not line.startswith(('#', '//')) and ':' in line and '=' in line:
            key, value = line.split('=', 1)
            result[key.split(':')[0]] = value
    return result


def build(args, revision, source, *, backends=False):
    expected = subprocess.check_output(['git', 'rev-parse', revision], cwd=ROOT, text=True).strip()
    if not args.baseline.exists():
        command(['git', 'worktree', 'add', '--detach', args.baseline, expected])
    head = subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=args.baseline, text=True).strip()
    if head != expected or subprocess.check_output(['git', 'status', '--porcelain'], cwd=args.baseline, text=True):
        raise RuntimeError(f'baseline must be a clean checkout of {revision}')
    command(['cmake', '--preset', 'release'])
    cache = cmake_cache(ROOT/'build/release/CMakeCache.txt')
    keys = ['CMAKE_C_COMPILER', 'CMAKE_CXX_COMPILER', 'CMAKE_CXX_FLAGS',
            'CMAKE_CXX_FLAGS_RELEASE', 'CMAKE_C_FLAGS', 'CMAKE_C_FLAGS_RELEASE',
            'CMAKE_OSX_ARCHITECTURES', 'CMAKE_OSX_DEPLOYMENT_TARGET', 'Eigen3_DIR', 'Boost_DIR']
    matched = {key: cache[key] for key in keys if key in cache}
    matched['BUILD_SHARED_LIBS'] = 'OFF'
    # The harnesses link archives explicitly, even if the user's release was shared.
    definitions = [f'-D{k}={v}' for k, v in matched.items()]
    command(['cmake', '--preset', 'release', *definitions])
    command(['cmake', '--preset', 'release', *definitions], cwd=args.baseline)
    for tree in [ROOT, args.baseline]:
        command(['cmake', '--build', '--preset', 'release', '-j', '6'], cwd=tree)
    extra = []
    if backends:
        nauty_source = cache.get('FETCHCONTENT_SOURCE_DIR_NAUTY') or ROOT/'build/release/_deps/nauty-src'
        extra = [f'-Dnauty_SOURCE={nauty_source}', f'-Dnauty_BINARY={ROOT}/build/release/_deps/nauty-build']
    command(['cmake', '-S', source, '-B', args.output/'bin', '-G', 'Ninja',
             '-DCMAKE_BUILD_TYPE=Release', f'-Dprevious_ROOT={args.baseline}',
             f'-Dcurrent_ROOT={ROOT}', *extra, *definitions])
    command(['cmake', '--build', args.output/'bin', '-j', '6'])
    return dict(baseline_commit=expected, platform=platform.platform(), build=matched,
                compiler=subprocess.check_output([cache['CMAKE_CXX_COMPILER'], '--version'], text=True))


def source_digest(source):
    paths = [*(ROOT/'src').glob('*'), *(ROOT/'include/pcut').glob('*.hpp'),
             *(ROOT/'cmake').glob('*'), ROOT/'CMakeLists.txt',
             *Path(__file__).parent.glob('comparison_*.py'),
             ROOT/'scripts/engine-comparison/measure.cpp', *source.glob('*')]
    digest = hashlib.sha256()
    for path in sorted(set(paths)):
        if path.is_file():
            digest.update(str(path.relative_to(ROOT)).encode())
            digest.update(path.read_bytes())
    return digest.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def describe(values):
    center = statistics.median(values)
    return dict(median=center, minimum=min(values), maximum=max(values),
                mad=statistics.median(abs(x-center) for x in values))


def summarize(rows, keys):
    return {key: describe([row[key] for row in rows]) for key in keys}


def peak_mib(value):
    return int(value)/(2**20 if platform.system() == 'Darwin' else 1024)


def parse_scalar(output):
    fields = output.split()
    return dict(seconds=float(fields[0]), peak_mib=peak_mib(fields[1]),
                calls=int(fields[2]), checksum=fields[3:])


def parse_sweep(output):
    row = dict(steps=[])
    for line in output.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == 'setup':
            row['coefficient_seconds'], row['topology_seconds'] = map(float, fields[1:])
        elif fields[0] == 'cold':
            row['cold_seconds'] = float(fields[1])
        elif fields[0] == 'step':
            row['steps'].append(dict(index=int(fields[1]), bind_seconds=float(fields[2]),
                link_seconds=float(fields[3]), graphs=int(fields[4]), embeddings=int(fields[5]),
                evaluations=int(fields[6]), values=list(map(float, fields[7:]))))
        elif fields[0] == 'total':
            row['sweep_seconds'] = float(fields[1])
            row['peak_mib'] = peak_mib(fields[2])
    for key in ['coefficient_seconds', 'topology_seconds', 'cold_seconds', 'sweep_seconds', 'peak_mib']:
        if key not in row:
            raise ValueError(f'missing sweep field: {key}')
    row['warm_seconds'] = statistics.median(s['bind_seconds']+s['link_seconds'] for s in row['steps'][1:])
    return row


def samples(args, engines, make_command, parse, metadata, filename, report):
    """Yield complete pairs; each subprocess finishes before the next starts."""
    for repeat in range(args.repeats):
        pair = {}
        for engine in (engines if repeat % 2 == 0 else engines[::-1]):
            start = time.perf_counter()
            output = subprocess.check_output(list(map(str, make_command(engine))), text=True)
            wall = time.perf_counter()-start
            (args.output/f'{filename}-{engine}-{repeat}.txt').write_text(output)
            row = dict(metadata, engine=engine, repeat=repeat, process_wall_seconds=wall, **parse(output))
            report['runs'].append(row)
            pair[engine] = row
            write_json(args.output/'results.json', report)
            print(filename, engine, repeat, row.get('seconds', row.get('sweep_seconds')), flush=True)
        yield repeat, pair


def compare_sweeps(pair, counts):
    error = 0.0
    for old, new in zip(pair['previous']['steps'], pair['current']['steps'], strict=True):
        for key in ['index', *counts]:
            if old[key] != new[key]:
                raise ValueError(f'sweep {key} mismatch: {old[key]} != {new[key]}')
        for x, y in zip(old['values'], new['values'], strict=True):
            if not math.isclose(x, y, rel_tol=1e-10, abs_tol=1e-10):
                raise ValueError(f'physics mismatch: {x} != {y}')
            error = max(error, abs(x-y)/max(1, abs(x), abs(y)))
    return error
