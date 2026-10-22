#!/usr/bin/env python3
"""Separate numerical polish trials at saved null points; no new full fits."""
import argparse
import ast
from dataclasses import replace
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
    'repair_p3_bb_gamma_difference_v3.py','repair_p3_score_integral_polish.py')
POINTS=(('proof0',REC/'repair_p3_gamma_order_draw_0_22083807','p3_000.json'),
    ('repeat0',REC/'repair_p3_gamma_order_family_draw_0_22084113','p3_000.json'),
    ('failed11',REC/'repair_p3_gamma_GL64_draw_11_22082995','p3_011.json'))


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
    required=[*(HERE/n for n in SOURCES),*(d/n for _,d,f in POINTS for n in (f,'draw_counts.json')),
        *(FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'))]
    require(all(str(p) in pins for p in required),'Input/source not pinned')
    if a.inspect:print(json.dumps(dict(saved_points=3,new_full_fits=0)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
        'Individual compute allocation required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
        'Threads not fixed')
    out=REC/('score_integral_polish_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        q=archive/p.name;q.write_bytes(p.read_bytes());require(sha(q)==(a.expected_runner_sha if p==Path(__file__) else pins[str(p)]),
            'Archive differs')
    import numpy as np
    import scipy
    sys.path.insert(0,str(archive));loader=importlib.import_module('bfix_v1')
    original=importlib.import_module('likelihood_v1');gamma=importlib.import_module('repair_p3_bb_gamma_difference_v3')
    candidate=importlib.import_module('repair_p3_score_integral_polish');design,settings=loader.load();template=design.model()
    read=lambda p:json.loads(p.read_text());began=time.perf_counter();records=[]
    def save(name,value):
        with (out/(name+'.json')).open('x') as f:json.dump(value,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    save('protocol',dict(source_and_input_sha256=pins,node=os.environ.get('SLURMD_NODENAME'),
        environment=dict(python=platform.python_version(),numpy=np.__version__,scipy=scipy.__version__),
        scientific_convergence_cutoff=1e-6,original_polish_iterations=100,candidate_polish_iterations=20,
        no_full_fit=True,no_original_result_reclassification=True,no_calibration=True,no_human_data=True))
    def score_max(model,p):
        ll,g,h=model.evaluate(p);bounds=model.bounds();pg=g.copy()
        for j,(lo,hi) in enumerate(bounds):
            if ((lo is not None and lo==hi) or (lo is not None and p[j]<=lo+1e-8 and g[j]<0)
                or (hi is not None and p[j]>=hi-1e-8 and g[j]>0)):pg[j]=0.
        return dict(parameters=p.tolist(),loglik=ll,gradient=g.tolist(),hessian=h.tolist(),
            projected_gradient=pg.tolist(),max_projected_gradient=float(np.max(np.abs(pg))))
    def original_increment(model,p,d,order):
        nodes,weights=np.polynomial.legendre.leggauss(order);parts=[]
        for x,w in zip(nodes,weights):parts.append(float(.5*w)*float(np.dot(model.evaluate(p+(.5+.5*x)*d,need_hessian=False)[1],d)))
        return math.fsum(parts),math.fsum(abs(x) for x in parts)
    for label,d,filename in POINTS:
        counts=read(d/'draw_counts.json');require(len(counts['counts'])==len(template.rows)==173,'Count vector differs')
        require([[r.individual,settings['gene_order'][r.gene]] for r in template.rows]==counts['full_row_keys'],'Row order differs')
        model=original.Model([replace(row,a=count) for row,count in zip(template.rows,counts['counts'])],gene_count=8)
        start=np.asarray(read(d/filename)['fit']['null']['parameters']);start.setflags(write=False);before=score_max(model,start)
        bounds=[b if keep else (start[j],start[j])
            for j,(b,keep) in enumerate(zip(model.bounds(),model.informative()))]
        before_boundary=float(model.marginal(start,0.)[1][-1])
        print(json.dumps(dict(event='original_polish_started',point=label)),flush=True)
        old,old_history=model.constrained_polish(start,bounds);old_result=score_max(model,old)
        baseline_boundary=float(model.marginal(old,0.)[1][-1])
        print(json.dumps(dict(event='score_integral_polish_started',point=label)),flush=True)
        new,evidence=candidate.polish(model,start,bounds,original,gamma);new_result=score_max(model,new)
        independent=[]
        for iteration in evidence['history']:
            for trial in iteration['trials']:
                if not trial.get('accepted'):continue
                p=np.asarray(iteration['parameters']);delta=np.asarray(trial['parameters'])-p
                coarse,_=original_increment(model,p,delta,16);fine,absolute=original_increment(model,p,delta,32)
                lower=fine-8*abs(fine-coarse)-64*np.finfo(float).eps*absolute
                noise_proxy=abs(fine-trial['score_integral32'])
                independent.append(dict(original_score_integral16=coarse,original_score_integral32=fine,
                    original_empirical_gain_lower_estimate=lower,Gamma_score_integral32=trial['score_integral32'],
                    arithmetic_difference=fine-trial['score_integral32'],arithmetic_noise_proxy=noise_proxy,
                    both_gains_exceed_32_times_observed_arithmetic_difference=bool(
                        min(lower,trial['empirical_gain_lower_estimate'])>32*noise_proxy),
                    empirical_not_certified=True))
        _,boundary=model.marginal(new,0.)
        expected_repair=label!='proof0'
        passed=(evidence['converged'] and new_result['max_projected_gradient']<1e-6
            and abs(new_result['loglik']-before['loglik'])<=1e-6 and abs(new[16]-start[16])<=1e-4
            and (not expected_repair or (before['max_projected_gradient']>=1e-6 and len(independent)>0))
            and all(x['both_gains_exceed_32_times_observed_arithmetic_difference'] for x in independent))
        record=dict(point=label,before=before,original_polish=old_result,original_history=old_history,
            candidate=new_result,candidate_evidence=evidence,independent_original_score_integrals=independent,
            boundary_variance_score=float(boundary[-1]),before_boundary_variance_score=before_boundary,
            baseline_boundary_variance_score=baseline_boundary,boundary_compatible=bool(boundary[-1]<=1e-6),
            expected_repair=expected_repair,accepted_steps=len(independent),candidate_passed=bool(passed),
            LL_difference=new_result['loglik']-before['loglik'],kappa_difference=float(new[16]-start[16]))
        save(label,record);records.append(record)
        print(json.dumps(dict(point=label,candidate_passed=bool(passed),max_projected_gradient=new_result['max_projected_gradient'])),flush=True)
    for p,h in pins.items():require(sha(Path(p))==h,'Input changed during evaluation')
    passed=all(r['candidate_passed'] for r in records)
    save('summary',dict(completed=True,points=3,all_candidate_points_passed=passed,
        accepted_steps={r['point']:r['accepted_steps'] for r in records},
        all_boundary_scores_compatible=all(r['boundary_compatible'] for r in records),seconds=time.perf_counter()-began,
        no_full_fit=True,no_original_result_reclassification=True,no_calibration=True))
    print(json.dumps(dict(completed=True,all_candidate_points_passed=passed)))
    if not passed:raise SystemExit(1)


if __name__=='__main__':main()
