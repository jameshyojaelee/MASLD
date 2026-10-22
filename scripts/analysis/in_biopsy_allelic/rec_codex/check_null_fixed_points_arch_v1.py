#!/usr/bin/env python3
"""Evaluate two saved synthetic null points; no optimizer or resampling."""
import argparse
import ast
import contextlib
from dataclasses import replace
import hashlib
import importlib
import io
import json
import os
from pathlib import Path
import platform
import sys

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
PROOF=REC/'repair_p3_gamma_order_draw_0_22083807'
REPEAT=REC/'repair_p3_gamma_order_family_draw_0_22084113'
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py','repair_p3_bb_gamma_difference_v3.py')


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True);ap.add_argument('--inspect',action='store_true');a=ap.parse_args()
    require(sha(Path(__file__))==a.expected_runner_sha and sha(a.pins)==a.expected_pins_sha,'Runner/pins differ')
    pins=json.loads(a.pins.read_text())
    for p,h in pins.items():
        path=Path(p);require(path.is_absolute() and path.resolve().is_relative_to(ROOT.resolve()) and sha(path)==h,'Input differs')
        if path.suffix=='.py':ast.parse(path.read_text())
    required=[*(HERE/n for n in SOURCES),*(d/n for d in (PROOF,REPEAT) for n in ('p3_000.json','draw_counts.json')),
        *(FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'))]
    require(all(str(p) in pins for p in required),'Input/source not pinned')
    if a.inspect:print(json.dumps(dict(points=2,BB_arithmetics=2,no_fits=True)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
        'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
        'Threads not fixed')
    out=REC/('null_fixed_points_arch_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        q=archive/p.name;q.write_bytes(p.read_bytes());require(sha(q)==(a.expected_runner_sha if p==Path(__file__) else pins[str(p)]),
            'Archive differs')
    import numpy as np
    import scipy
    sys.path.insert(0,str(archive));loader=importlib.import_module('bfix_v1')
    original=importlib.import_module('likelihood_v1');gamma=importlib.import_module('repair_p3_bb_gamma_difference_v3')
    design,settings=loader.load();template=design.model()
    read=lambda p:json.loads(p.read_text());counts=read(PROOF/'draw_counts.json')
    require(counts==read(REPEAT/'draw_counts.json') and len(counts['counts'])==len(template.rows)==173,
        'Saved synthetic draw counts differ')
    require([[r.individual,settings['gene_order'][r.gene]] for r in template.rows]==counts['full_row_keys'],
        'Saved row order differs')
    model=original.Model([replace(row,a=count) for row,count in zip(template.rows,counts['counts'])],gene_count=8)
    require(model.size==28 and model.Z==7,'Synthetic model shape differs')
    bb=original.bb_terms;results=[]
    for label,d in (('proof',PROOF),('repeat',REPEAT)):
        fit=read(d/'p3_000.json')['fit'];p=np.asarray(fit['null']['parameters'],dtype=float);p.setflags(write=False)
        for method,function in (('original',bb),('Gamma3',gamma.stable_bb_terms)):
            try:
                original.bb_terms=function
                ll,g,h=model.evaluate(p);ml,mg=model.marginal(p,0.)
            finally:original.bb_terms=bb
            require(np.isfinite(ll) and np.isfinite(g).all() and np.isfinite(h).all()
                and np.isfinite(ml) and np.isfinite(mg).all(),'Fixed-point arrays invalid')
            require(ll==ml and np.array_equal(g,mg[:-1]),'Evaluate/boundary marginal coordinates differ')
            bounds=model.bounds();projected=g.copy()
            for j,(lo,hi) in enumerate(bounds):
                if (lo is not None and lo==hi) or (lo is not None and p[j]<=lo+1e-8 and g[j]<0) or (
                    hi is not None and p[j]>=hi-1e-8 and g[j]>0):projected[j]=0.
            results.append(dict(point=label,BB_arithmetic=method,parameters=p.tolist(),loglik=ll,
                gradient=g.tolist(),hessian=h.tolist(),boundary_score=float(mg[-1]),
                projected_gradient=projected.tolist(),max_projected_gradient=float(np.max(np.abs(projected))),
                original_saved_LL=fit['null']['loglik'],original_saved_gradient=fit['null']['gradient']))
    runtime=io.StringIO()
    with contextlib.redirect_stdout(runtime):np.show_config();np.show_runtime()
    (out/'numpy_runtime.txt').write_text(runtime.getvalue())
    cpu=[line for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')]
    report=dict(completed=True,node=os.environ.get('SLURMD_NODENAME'),cpu_model=cpu[0] if cpu else None,
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        source_and_input_sha256=pins,results=results,no_fit_or_resampling=True,no_human_data=True)
    for p,h in pins.items():require(sha(Path(p))==h,'Input changed during evaluation')
    with (out/'results.json').open('x') as f:json.dump(report,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    print(json.dumps(dict(completed=True,node=report['node'],points=2,BB_arithmetics=2)))


if __name__=='__main__':main()
