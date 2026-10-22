#!/usr/bin/env python3
"""Two separate synthetic null refits using the checked Gamma BB precision.

The original fitter, cold/warm starts, tie and convergence rule are inherited.
After restoring original BB arithmetic, independently evaluate the endpoint.
No existing draw is relabeled and no positive integral or CRT is run.
"""
import argparse
import ast
from dataclasses import replace
import hashlib
import importlib
import json
import os
from pathlib import Path
import sys
import time

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
WARM=REC/'repair_p3_gamma_observed_batched_GL_22080549/observed_null.json'
POINTS={0:(REC/'repair_p3_gamma_order_family_draw_0_22084113','p3_000.json'),
    11:(REC/'repair_p3_gamma_GL64_draw_11_22082995','p3_011.json')}
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py','repair_p3_bb_gamma_difference_v3.py')


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--draw',type=int,choices=(0,11),required=True)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True);ap.add_argument('--inspect',action='store_true');a=ap.parse_args()
    require(sha(Path(__file__))==a.expected_runner_sha and sha(a.pins)==a.expected_pins_sha,'Runner/pins differ')
    pins=json.loads(a.pins.read_text());read=lambda p:json.loads(p.read_text())
    for path,h in pins.items():
        p=Path(path);require(p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,'Input differs')
        if p.suffix=='.py':ast.parse(p.read_text())
    required=[WARM,*(HERE/n for n in SOURCES),*(d/n for d,f in POINTS.values() for n in (f,'draw_counts.json')),
        *(FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'))]
    require(all(str(p) in pins for p in required),'Missing source/input pin')
    require(pins[str(HERE/'likelihood_v1.py')]=='41ebc2531eb6ff208bf3291f061b5d251267a8781ca05a4be4b1af304039cbc6'
        and pins[str(HERE/'repair_p3_bb_gamma_difference_v3.py')]=='546e7ea2f1db21e6fcdd6aed72204d3bbe48c2d818996ab21de5f5ab5215bcf3',
        'Previously checked BB/fitter sources required')
    warm=read(WARM)['fit'];require(warm['converged'],'Observed original warm point incomplete')
    if a.inspect:print(json.dumps(dict(draw=a.draw,null_refits=1,starts='original warm and cold',positive_integrals=0)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
        'Single numerical thread required')
    out=REC/('Gamma_null_refit_draw'+str(a.draw)+'_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        q=archive/p.name;q.write_bytes(p.read_bytes());require(sha(q)==(a.expected_runner_sha if p==Path(__file__) else pins[str(p)]),'Archive differs')
    import numpy as np
    import scipy
    sys.path.insert(0,str(archive));original=importlib.import_module('likelihood_v1');loader=importlib.import_module('bfix_v1')
    gamma=importlib.import_module('repair_p3_bb_gamma_difference_v3');design,settings=loader.load();template=design.model()
    require(loader.FIXTURE.resolve()==FIXTURE.resolve(),'Fixture loader path differs')
    directory,filename=POINTS[a.draw];counts=read(directory/'draw_counts.json');previous=read(directory/filename)['fit']['null']
    require(counts['full_row_keys']==[[r.individual,settings['gene_order'][r.gene]] for r in template.rows]
        and len(counts['counts'])==len(template.rows)==173,'Count row order differs')
    model=original.Model([replace(r,a=c) for r,c in zip(template.rows,counts['counts'])],gene_count=8)
    require(model.G==8 and model.Z==7 and model.size==28,'Synthetic model dimension differs')
    start=np.asarray(warm['parameters']);start.setflags(write=False);before=start.tobytes();rows_before=repr(model.rows)
    bb=original.bb_terms;began=time.perf_counter()
    def safe(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,dict):return {k:safe(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [safe(v) for v in x]
        return x
    def save(name,value):
        with (out/(name+'.json')).open('x') as f:json.dump(safe(value),f,indent=2,allow_nan=False);f.flush();os.fsync(f.fileno())
    save('protocol',dict(draw=a.draw,source_and_input_sha256=pins,node=os.environ.get('SLURMD_NODENAME'),
        environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
        scientific_convergence_cutoff=1e-6,starts_tie_and_fitter_unchanged=True,
        changed_only_BB_special_function_precision=True,no_existing_result_reclassification=True,
        no_positive_integral=True,no_CRT=True,no_calibration=True,no_human_data=True))
    try:
        original.bb_terms=gamma.stable_bb_terms
        fit=model.fit(start=start)
    finally:original.bb_terms=bb
    seconds=time.perf_counter()-began;save('Gamma_fit',fit)
    p=np.asarray(fit['parameters']);ll,g,h=model.evaluate(p);bounds=model.bounds();pg=g.copy()
    require(p.shape==g.shape==(28,) and h.shape==(28,28) and np.isfinite(p).all()
        and np.isfinite(ll) and np.isfinite(g).all() and np.isfinite(h).all(),'Nonfinite original endpoint')
    for j,(lo,hi) in enumerate(bounds):
        if (lo is not None and lo==hi) or (lo is not None and p[j]<=lo+1e-8 and g[j]<0) or (hi is not None and p[j]>=hi-1e-8 and g[j]>0):pg[j]=0.
    maximum=float(np.max(np.abs(pg)));difference=float(ll-fit['loglik'])
    gradient_ok=bool(np.allclose(g,fit['gradient'],atol=1e-8,rtol=1e-6))
    expected_keep='cold' if fit['start_logliks']['cold']-fit['start_logliks']['warm']>=1e-6 else 'warm'
    require(fit['kept_start']==expected_keep and start.tobytes()==before and repr(model.rows)==rows_before
        and original.bb_terms is bb,'Original starts/tie/inputs changed')
    passed=bool(fit['converged'] and maximum<1e-6 and abs(difference)<=1e-6 and gradient_ok)
    save('original_endpoint',dict(parameters=p,loglik=ll,gradient=g,hessian=h,projected_gradient=pg,
        max_projected_gradient=maximum,Gamma_LL_difference=difference,
        Gamma_gradient_max_difference=float(np.max(np.abs(g-fit['gradient']))),Gamma_gradient_agreement=gradient_ok))
    for path,hsh in pins.items():require(sha(Path(path))==hsh,'Input changed during refit')
    save('summary',dict(completed=True,draw=a.draw,passed_original_endpoint_check=passed,fit_seconds=seconds,
        Gamma_converged=fit['converged'],original_max_projected_gradient=maximum,
        previous_failed_null_LL_difference=float(ll-previous['loglik']),
        previous_failed_null_kappa_difference=float(p[16]-previous['parameters'][16]),
        no_controlled_speed_comparison=True,no_positive_integral=True,no_existing_result_reclassification=True,
        no_fifty_draw_agreement=True,no_calibration=True,no_human_validation=True))
    print(json.dumps(dict(completed=True,draw=a.draw,passed_original_endpoint_check=passed,seconds=seconds)))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
