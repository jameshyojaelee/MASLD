#!/usr/bin/env python3
"""Replay the saved draw11 warm point through the original solver first."""
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
POINT=REC/'repair_p3_gamma_order_polish_draw_11_failed11_22084421'
PREVIOUS=REC/'repair_p3_gamma_GL64_draw_11_22082995/p3_011.json'
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
    'repair_p3_bb_gamma_difference_v3.py','repair_p3_score_integral_polish.py',
    'repair_p3_gamma_order_polish_model_v2.py')


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True);ap.add_argument('--inspect',action='store_true');a=ap.parse_args()
    require(sha(Path(__file__))==a.expected_runner_sha and sha(a.pins)==a.expected_pins_sha,'Runner/pins differ')
    pins=json.loads(a.pins.read_text());read=lambda p:json.loads(p.read_text())
    for path,h in pins.items():
        p=Path(path);require(p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,'Input differs')
        if p.suffix=='.py':ast.parse(p.read_text())
    required=[*(HERE/n for n in SOURCES),PREVIOUS,*(POINT/n for n in ('null_polish_evidence.json','draw_counts.json')),
        *(FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'))]
    require(all(str(p) in pins for p in required),'Missing source/input pin')
    old=read(POINT/'null_polish_evidence.json')['trials'];require(len(old)==1 and not old[0]['converged'],
        'Saved nonconverged warm path differs')
    start_values=old[0]['history'][0]['parameters'];require(len(start_values)==28,'Warm dimension differs')
    if a.inspect:print(json.dumps(dict(saved_warm_point=True,original_constrained_solver_first=True,no_full_fit=True)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
        'Single numerical thread required')
    out=REC/('original_first_warm_polish_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        q=archive/p.name;q.write_bytes(p.read_bytes());require(sha(q)==(a.expected_runner_sha if p==Path(__file__) else pins[str(p)]),'Archive differs')
    import numpy as np
    sys.path.insert(0,str(archive));original=importlib.import_module('likelihood_v1');loader=importlib.import_module('bfix_v1')
    gamma=importlib.import_module('repair_p3_bb_gamma_difference_v3');polish=importlib.import_module('repair_p3_score_integral_polish')
    revised=importlib.import_module('repair_p3_gamma_order_polish_model_v2');design,settings=loader.load();template=design.model()
    counts=read(POINT/'draw_counts.json');require(len(counts['counts'])==len(template.rows)==173
        and counts['full_row_keys']==[[r.individual,settings['gene_order'][r.gene]] for r in template.rows],'Count row order differs')
    Model=revised.model_class(original.Model,original,gamma,polish,
        expected_polish_sha=pins[str(HERE/'repair_p3_score_integral_polish.py')])
    model=Model([replace(r,a=c) for r,c in zip(template.rows,counts['counts'])],gene_count=8)
    start=np.asarray(start_values);start.setflags(write=False);before=start.tobytes();rows=repr(model.rows)
    bounds=[b if keep else (start[j],start[j]) for j,(b,keep) in enumerate(zip(model.bounds(),model.informative()))]
    def evaluate(p):
        ll,g,h=model.evaluate(p);pg=g.copy()
        for j,(lo,hi) in enumerate(bounds):
            if (lo is not None and lo==hi) or (lo is not None and p[j]<=lo+1e-8 and g[j]<0) or (hi is not None and p[j]>=hi-1e-8 and g[j]>0):pg[j]=0.
        require(np.isfinite(ll) and np.isfinite(g).all() and np.isfinite(h).all(),'Nonfinite original evaluation')
        return dict(parameters=p.tolist(),loglik=float(ll),gradient=g.tolist(),hessian=h.tolist(),
            projected_gradient=pg.tolist(),max_projected_gradient=float(np.max(np.abs(pg))))
    initial=evaluate(start);began=time.perf_counter();end,_=model.constrained_polish(start,bounds)
    seconds=time.perf_counter()-began;evidence=model.null_polish_evidence[-1]
    baseline=evaluate(np.asarray(evidence['original_constrained_baseline_parameters']));result=evaluate(end)
    previous=read(PREVIOUS)['fit']['null'];passed=bool(evidence['original_constrained_solver_retained']
        and evidence['converged'] and result['max_projected_gradient']<1e-6
        and abs(result['loglik']-previous['loglik'])<=1e-6)
    require(start.tobytes()==before and repr(model.rows)==rows,'Fixed inputs changed')
    for path,h in pins.items():require(sha(Path(path))==h,'Input changed during polish')
    with (out/'summary.json').open('x') as f:
        json.dump(dict(completed=True,passed=passed,node=os.environ.get('SLURMD_NODENAME'),seconds=seconds,
            initial=initial,original_solver_baseline=baseline,after_refinement=result,evidence=evidence,
            previous_failed_warm_LL_difference=result['loglik']-previous['loglik'],source_and_input_sha256=pins,
            no_full_fit=True,no_existing_result_reclassification=True,no_calibration=True,no_human_data=True),f,indent=2,allow_nan=False)
        f.flush();os.fsync(f.fileno())
    print(json.dumps(dict(completed=True,passed=passed,max_projected_gradient=result['max_projected_gradient'])))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
