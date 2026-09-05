#!/usr/bin/env python3
"""Matched a4ea23c comparison; run serially after release/sanitizer validation."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import platform
import statistics
import subprocess
ROOT=Path(__file__).resolve().parents[2]
SOURCE=Path(__file__).resolve().parent
BASE='a4ea23c88325141eed0e2f2cff387546214b976f'
spec=importlib.util.spec_from_file_location('common',ROOT/'scripts/engine-comparison/run.py')
common=importlib.util.module_from_spec(spec); spec.loader.exec_module(common)
def run(args,cwd=ROOT): subprocess.run(list(map(str,args)),cwd=cwd,check=True)
def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline',type=Path,default=Path('/tmp/pcut-canonical-baseline'))
    parser.add_argument('--output',type=Path,default=ROOT/'build/canonical-comparison')
    parser.add_argument('--repeats',type=int,default=7)
    args=parser.parse_args(); args.output=args.output.resolve(); args.baseline=args.baseline.resolve()
    args.output.mkdir(parents=True,exist_ok=True)
    if not args.baseline.exists(): run(['git','worktree','add','--detach',args.baseline,BASE])
    if subprocess.check_output(['git','rev-parse','HEAD'],cwd=args.baseline,text=True).strip()!=BASE or subprocess.check_output(['git','status','--porcelain'],cwd=args.baseline,text=True):
        raise RuntimeError('baseline must be the clean starting revision')
    run(['cmake','--preset','release'])
    cache=common.cmake_cache(ROOT/'build/release/CMakeCache.txt')
    keys=['CMAKE_C_COMPILER','CMAKE_CXX_COMPILER','CMAKE_CXX_FLAGS','CMAKE_CXX_FLAGS_RELEASE','CMAKE_C_FLAGS','CMAKE_C_FLAGS_RELEASE','CMAKE_OSX_ARCHITECTURES','CMAKE_OSX_DEPLOYMENT_TARGET','Eigen3_DIR','Boost_DIR']
    matched={k:cache[k] for k in keys if k in cache}
    definitions=[f'-D{k}={v}' for k,v in matched.items()]
    run(['cmake','--preset','release',*definitions],args.baseline)
    run(['cmake','--build','--preset','release','-j','6'],args.baseline)
    run(['cmake','--build','--preset','release','-j','6'])
    nauty_source=cache.get('FETCHCONTENT_SOURCE_DIR_NAUTY') or str(ROOT/'build/release/_deps/nauty-src')
    run(['cmake','-S',SOURCE,'-B',args.output/'bin','-G','Ninja','-DCMAKE_BUILD_TYPE=Release',f'-Dprevious_ROOT={args.baseline}',f'-Dcurrent_ROOT={ROOT}',f'-Dnauty_SOURCE={nauty_source}',f'-Dnauty_BINARY={ROOT}/build/release/_deps/nauty-build',*definitions])
    run(['cmake','--build',args.output/'bin','-j','6'])
    digest=hashlib.sha256()
    paths=sorted([*(ROOT/'src').glob('*'),*(ROOT/'include/pcut').glob('*.hpp'),*(ROOT/'cmake').glob('*'),ROOT/'CMakeLists.txt',ROOT/'tests/canonical_oracle.hpp',ROOT/'scripts/engine-comparison/measure.cpp',*SOURCE.glob('*')])
    for path in paths:
        if path.is_file(): digest.update(str(path.relative_to(ROOT)).encode()); digest.update(path.read_bytes())
    report=dict(baseline_commit=BASE,source_sha256=digest.hexdigest(),platform=platform.platform(),compiler=subprocess.check_output([cache['CMAKE_CXX_COMPILER'],'--version'],text=True),build=matched,repeats=args.repeats,runs=[])
    rss_scale=2**20 if platform.system()=='Darwin' else 1024
    def record(row):
        report['runs'].append(row)
        (args.output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
        print(row['kind'],row['case'],row['engine'],row.get('seconds',row.get('sweep_seconds')),flush=True)
    for kind in ['backends','isolated','topology','sweep']:
        cases=common.CASES if kind in ['sweep','topology'] else [(c,0,0,0) for c in ['square','four_color','hubbard','dimer','star8','cycle8','complete7','duplicates']+(['star20','cycle64'] if kind=='backends' else [])]
        for name,order,count,particles in cases:
            engines=['nauty','traces'] if kind=='backends' else ['previous','current']
            for repeat in range(args.repeats):
                pair={}
                for engine in (engines if repeat%2==0 else engines[::-1]):
                    iterations=1 if kind=='topology' else (100 if name in ['square','four_color','hubbard','dimer'] else (100 if kind=='backends' else 1))
                    command=([args.output/'bin/backends',engine,name,iterations] if kind=='backends' else [args.output/'bin'/f'{engine}_{kind}',name,*([order,count,particles] if kind=='sweep' else [iterations])])
                    output=subprocess.check_output(list(map(str,command)),text=True)
                    (args.output/f'{kind}-{name}-{engine}-{repeat}.txt').write_text(output)
                    row=dict(kind=kind,case=name,engine=engine,repeat=repeat)
                    if kind!='sweep':
                        f=output.split(); row.update(seconds=float(f[0]),peak_mib=int(f[1])/rss_scale,calls=int(f[2]),checksum=f[3:])
                    else:
                        row['steps']=[]
                        for line in output.splitlines():
                            f=line.split()
                            if f[0]=='setup': row['coefficient_seconds'],row['topology_seconds']=map(float,f[1:])
                            elif f[0]=='cold': row['cold_seconds']=float(f[1])
                            elif f[0]=='total': row['sweep_seconds']=float(f[1]); row['peak_mib']=int(f[2])/rss_scale
                            elif f[0]=='step': row['steps'].append(dict(index=int(f[1]),bind_seconds=float(f[2]),link_seconds=float(f[3]),graphs=int(f[4]),embeddings=int(f[5]),evaluations=int(f[6]),values=list(map(float,f[7:]))))
                    pair[engine]=row; record(row)
                if kind=='sweep':
                    error=0
                    for old,new in zip(pair['previous']['steps'],pair['current']['steps'],strict=True):
                        assert old['embeddings']==new['embeddings'] and old['graphs']==new['graphs'] and old['evaluations']==new['evaluations']
                        for x,y in zip(old['values'],new['values'],strict=True):
                            assert math.isclose(x,y,rel_tol=1e-10,abs_tol=1e-10),(name,x,y)
                            error=max(error,abs(x-y)/max(1,abs(x),abs(y)))
                    report.setdefault('equivalence',[]).append(dict(case=name,repeat=repeat,scaled_error=error))
                else: assert pair[engines[0]]['checksum']==pair[engines[1]]['checksum']
    report['summary']=[]
    for kind in ['backends','isolated','topology','sweep']:
        cases=sorted({r['case'] for r in report['runs'] if r['kind']==kind})
        for case in cases:
            rows=[r for r in report['runs'] if r['kind']==kind and r['case']==case]
            summary=dict(kind=kind,case=case)
            for engine in sorted({r['engine'] for r in rows}):
                selected=[r for r in rows if r['engine']==engine]
                keys=['peak_mib']+(['topology_seconds','cold_seconds','sweep_seconds'] if kind=='sweep' else ['seconds'])
                summary[engine]={k:common.describe([r[k] for r in selected]) for k in keys}
                if kind!='sweep': summary[engine]['calls']=selected[0]['calls']
            report['summary'].append(summary)
    (args.output/'results.json').write_text(json.dumps(report,indent=2)+'\n')
    # Retain all raw scalar timings and variability in the versioned report;
    # per-point physics output remains available in the build directory.
    compact={k:v for k,v in report.items() if k!='runs'}
    compact['runs']=[{k:v for k,v in r.items() if k!='steps'} for r in report['runs']]
    (ROOT/'docs/performance-canonical.json').write_text(json.dumps(compact,indent=2)+'\n')
if __name__=='__main__': main()
