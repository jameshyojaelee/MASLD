#!/usr/bin/env python3
"""Two saved synthetic gene integrals, unchanged final precision budgets."""
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
POINT=REC/'repair_p3_gamma_order_polish_draw_0_amd0_22084420'


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--gene',type=int,choices=(2,5),required=True)
    ap.add_argument('--pins',type=Path,required=True);ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True);ap.add_argument('--inspect',action='store_true');a=ap.parse_args()
    require(sha(Path(__file__))==a.expected_runner_sha and sha(a.pins)==a.expected_pins_sha,'Runner/pins differ')
    pins=json.loads(a.pins.read_text())
    for p,h in pins.items():
        p=Path(p);require(p.resolve().is_relative_to(ROOT.resolve()) and sha(p)==h,'Input differs')
        if p.suffix=='.py':ast.parse(p.read_text())
    required=[HERE/'repair_p3_local_gl_initial64.py',*(POINT/n for n in ('recipe_summary.json','draw_counts.json','protocol.json')),
        *(POINT/'executed_sources'/n for n in ('likelihood_v1.py','bfix_v1.py','design_v1.py',
            'repair_p3_integral.py','repair_p3_local_moments.py','repair_p3_gamma_batched_conditional.py',
            'repair_p3_bb_gamma_difference_v3.py')),
        *(FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'))]
    require(all(str(p) in pins for p in required),'Missing source/input pin')
    failure=json.loads((POINT/'recipe_summary.json').read_text())['fitting_integration_failure']
    failed={r['gene']:r for r in failure['all_worker_records'] if not r['success']}
    require(set(failed)=={2,5} and failure['call']==14
        and all(r['message']=='GL order refinement exceeds unchanged precision allowance' for r in failed.values()),
        'Saved refusal differs')
    if a.inspect:print(json.dumps(dict(gene=a.gene,saved_point=True,no_fit=True,final_budget_unchanged=True)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),'Individual compute job required')
    require(all(os.environ.get(k)=='1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),
        'Single numerical thread required')
    out=REC/('GL_initial64_refusal_gene'+str(a.gene)+'_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for path,h in pins.items():
        path=Path(path)
        if path.suffix=='.py':
            q=archive/path.name;require(not q.exists(),'Duplicate source basename');q.write_bytes(path.read_bytes());require(sha(q)==h,'Archive differs')
    import numpy as np
    sys.path.insert(0,str(archive));original=importlib.import_module('likelihood_v1');loader=importlib.import_module('bfix_v1')
    conditional=importlib.import_module('repair_p3_gamma_batched_conditional')
    gamma=importlib.import_module('repair_p3_bb_gamma_difference_v3')
    kernel=importlib.import_module('repair_p3_integral');gl=importlib.import_module('repair_p3_local_gl_initial64')
    design,settings=loader.load();template=design.model();counts=json.loads((POINT/'draw_counts.json').read_text())
    require(counts['full_row_keys']==[[r.individual,settings['gene_order'][r.gene]] for r in template.rows]
        and len(counts['counts'])==len(template.rows)==173,'Count row order differs')
    model=original.Model([replace(r,a=c) for r,c in zip(template.rows,counts['counts'])],gene_count=8)
    rows=tuple(r for r in model.rows if r.gene==a.gene);p=np.asarray(failure['parameters']);p.setflags(write=False)
    v=failure['variance'];before=p.tobytes();Context=conditional.context_class(archive);context=Context(model,p,rows)
    started=time.perf_counter();ll,g,h,report=gl.integrate_gene(model,p,v,rows,context.callback,
        with_hessian=False,expected_kernel_sha=sha(archive/'repair_p3_integral.py'),
        kernel_path=archive/'repair_p3_integral.py',settings=kernel.Settings(),batch_size=64)
    seconds=time.perf_counter()-started;e=report['empirical_refinement_plus_tail_allowances']
    require(np.isfinite(ll) and g.shape==(29,) and np.isfinite(g).all() and h is None
        and e[0]<=1e-9 and max(e[1])<=1e-8 and e[2] is None,'Final precision budget unresolved')
    indices=report['augmented_global_indices'];old=failed[a.gene]['refused_gene_ladders'][-1]
    old_difference=max(abs(x-y) for x,y in zip(g[indices],old['gradient']))
    native=[];bb=original.bb_terms
    selected=original.Model(rows,gene_count=8)
    try:
        original.bb_terms=gamma.stable_bb_terms
        for order in (40,80):
            nl,ng=selected.marginal(p,v,nodes=order)
            native.append(dict(order=order,logmass=float(nl),gradient=ng.tolist(),
                LL_difference=float(nl-ll),gradient_max_difference=float(np.max(np.abs(ng-g)))))
    finally:original.bb_terms=bb
    require(p.tobytes()==before,'Fixed parameters changed')
    for path,hsh in pins.items():require(sha(Path(path))==hsh,'Input changed during check')
    def safe(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,dict):return {k:safe(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [safe(v) for v in x]
        return x
    result=dict(completed=True,gene=a.gene,parameters=p,variance=v,rows=len(rows),logmass=float(ll),gradient=g,
        report=report,GL_seconds=seconds,old_initial64_gradient_max_difference=old_difference,
        native40_80=native,native_comparison_passed=all(abs(x['LL_difference'])<=1e-8
            and x['gradient_max_difference']<=1e-8 for x in native),source_and_input_sha256=pins,
        node=os.environ.get('SLURMD_NODENAME'),unchanged_final_precision=True,empirical_not_certified=True,
        no_fit=True,no_calibration=True,no_human_data=True)
    with (out/'summary.json').open('x') as f:json.dump(safe(result),f,indent=2,allow_nan=False);f.flush();os.fsync(f.fileno())
    print(json.dumps(dict(completed=True,gene=a.gene,native_comparison_passed=result['native_comparison_passed'])))
    if not result['native_comparison_passed']:raise SystemExit(1)


if __name__=='__main__':main()
