#!/usr/bin/env python3
"""Reproduce numerical dimer-chain benchmarks against cond-mat/9906243 Eq. E_grund."""
import argparse
import csv
import io
import json
from pathlib import Path
import platform
import subprocess
import time


def literature(alpha, order):
    a = alpha
    b = (1 - 2*a)**2
    polynomials = {
        2: 3/4,
        3: 3/4 + 3*a/2,
        4: 13/16 + 27*a/4 - 3*a*a/4,
        5: 89/48 + 311*a/24 + 93*a*a/4 - 45*a**3/2,
        6: 463/96 + 227*a/9 + 1307*a*a/12 - 42*a**3 - 159*a**4/2,
        7: 81557/6912 + 257909*a/3456 + 215995*a*a/864
           + 173579*a**3/432 - 14865*a**4/16 + 879*a**5/8,
        8: 414359/12960 + 139801*a/648 + 8477587*a*a/12960
           + 152558*a**3/81 - 2774357*a**4/1620 - 4002*a**5 + 4527*a**6/2,
    }
    return [-0.75, 0.0] + [-2*b*polynomials[n]/4**n for n in range(2, order+1)]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', default='build/release/pcut_chain')
    parser.add_argument('--output', default='build/benchmark.json')
    args = parser.parse_args()
    results = []
    for alpha, order, vacuum_only in [(0.0, 6, False), (0.17, 6, False),
                                       (0.5, 6, False), (0.17, 8, True)]:
        command = [args.executable, '--order', str(order), '--alpha', str(alpha),
                   '--lambda', '0.3', '--k', '0']
        if vacuum_only:
            command.append('--vacuum-only')
        start = time.perf_counter()
        process = subprocess.run(command, check=True, capture_output=True, text=True)
        elapsed = time.perf_counter() - start
        lines = '\n'.join(line for line in process.stdout.splitlines() if not line.startswith('#'))
        rows = list(csv.DictReader(io.StringIO(lines)))
        computed = [float(row['energy_per_dimer']) for row in rows]
        expected = literature(alpha, order)
        if len(computed) != order+1:
            raise RuntimeError('unexpected coefficient count')
        error = max(abs(x-y) for x, y in zip(computed, expected))
        if error > 5e-11:
            raise RuntimeError(f'literature mismatch alpha={alpha}: {error}')
        results.append(dict(alpha=alpha, order=order, one_particle=not vacuum_only,
                            elapsed_seconds=elapsed, max_absolute_coefficient_error=error,
                            energy_per_dimer=computed, literature_energy_per_dimer=expected,
                            omega_k0=[float(r['omega_k']) for r in rows] if not vacuum_only else None))
        print(f'alpha={alpha:g} order={order}: max error {error:.3g}, elapsed {elapsed:.3f} s')
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(dict(platform=platform.platform(), machine=platform.machine(),
                                     executable=args.executable, strong_exchange_J=1,
                                     reference='cond-mat/9906243v1, Eq. E_grund', results=results), indent=2)+'\n')


if __name__ == '__main__':
    main()
