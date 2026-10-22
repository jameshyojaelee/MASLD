#!/usr/bin/env python3
"""Paired serial/parallel original fit_single checks at saved synthetic points.

Endpoint agreement checks scheduling only. A shared failed convergence flag
does not become a successful fit. No CRT, positive integral or human analysis.
"""
import argparse
import ast
from dataclasses import replace
import hashlib
import importlib
import json
import math
import os
from pathlib import Path
import sys
import time

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
POINTS = {0: (REC/'Gamma_null_refit_draw0_22084585/Gamma_fit.json',
              REC/'repair_p3_gamma_order_polish_draw_0_amd0_22084420/draw_counts.json'),
          11: (REC/'Gamma_null_refit_draw11_22084590/Gamma_fit.json',
               REC/'repair_p3_gamma_order_polish_draw_11_failed11_22084421/draw_counts.json')}
SOURCES = ('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py','repair_p3_parallel_null_shifts_v1.py')


def require(ok, message):
    if not ok: raise RuntimeError(message)


def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--draw',type=int,choices=(0,11),required=True)
    ap.add_argument('--pins',type=Path,required=True)
    ap.add_argument('--expected-pins-sha',required=True)
    ap.add_argument('--expected-runner-sha',required=True)
    ap.add_argument('--inspect',action='store_true')
    a = ap.parse_args()
    require(sha(Path(__file__)) == a.expected_runner_sha and sha(a.pins) == a.expected_pins_sha,'Runner/pins differ')
    pins = json.loads(a.pins.read_text())
    for path,digest in pins.items():
        p=Path(path); require(p.resolve().is_relative_to(ROOT.resolve()) and sha(p) == digest,'Input differs')
        if p.suffix == '.py': ast.parse(p.read_text())
    required=[*(HERE/n for n in SOURCES),*POINTS[a.draw],
        *(FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'))]
    require(all(str(p) in pins for p in required),'Missing source/input pin')
    if a.inspect:
        print(json.dumps(dict(draw=a.draw,shifts=[0.,-.1,.1],workers=3,unchanged_post_selection_AST=True,no_CRT=True)))
        return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID')
        and os.environ.get('SLURM_CPUS_PER_TASK') == '4','Individual four-CPU compute required')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'Threads differ')
    out=REC/f"parallel_null_shift_check_draw{a.draw}_{os.environ['SLURM_JOB_ID']}"
    out.mkdir(exist_ok=False); archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        q=archive/p.name;q.write_bytes(p.read_bytes())
        require(sha(q) == (a.expected_runner_sha if p == Path(__file__) else pins[str(p)]),'Archive differs')
    import numpy as np
    import scipy
    sys.path.insert(0,str(archive))
    modules={n[:-3]:importlib.import_module(n[:-3]) for n in SOURCES}
    require(all(Path(m.__file__).resolve().parent == archive.resolve() for m in modules.values()),'Import escaped archive')
    original,loader=modules['likelihood_v1'],modules['bfix_v1']
    require(loader.FIXTURE.resolve() == FIXTURE.resolve(),'Fixture differs')
    design,settings=loader.load();template=design.model()
    read=lambda p:json.loads(p.read_text())
    previous=read(POINTS[a.draw][0]); counts=read(POINTS[a.draw][1]);p=np.asarray(previous['parameters'])
    require(p.shape == (28,) and np.isfinite(p).all() and len(counts['counts']) == len(template.rows) == 173
        and counts['full_row_keys'] == [[r.individual,settings['gene_order'][r.gene]] for r in template.rows],'Synthetic input order differs')
    rows=[replace(r,a=c) for r,c in zip(template.rows,counts['counts'])]
    serial=original.Model(rows,gene_count=8)
    Parallel=modules['repair_p3_parallel_null_shifts_v1'].model_class(original.Model,original)
    parallel=Parallel(rows,gene_count=8)
    def safe(x):
        if isinstance(x,np.ndarray):return x.tolist()
        if isinstance(x,np.generic):return x.item()
        if isinstance(x,dict):return {k:safe(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [safe(v) for v in x]
        return x
    def save(name,value):
        with (out/(name+'.json')).open('x') as f:
            json.dump(safe(value),f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    save('protocol',dict(draw=a.draw,parameters=p,source_and_input_sha256=pins,
        environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),node=os.environ.get('SLURMD_NODENAME'),
        generated_method=Parallel.original_AST_parallel_change,serial_first=True,no_positive_integral=True,no_CRT=True,
        no_existing_draw_reclassification=True,no_calibration=True,no_human_validation=True))
    before,row_before=p.tobytes(),repr(rows); results=[]
    for label,model in (('serial',serial),('parallel',parallel)):
        began=time.perf_counter(); fit=model.fit_single(p)
        seconds=time.perf_counter()-began
        require(np.isfinite(fit['loglik']) and np.isfinite(fit['parameters']).all() and np.isfinite(fit['gradient']).all()
            and np.isfinite(fit['hessian']).all(),'Nonfinite original null endpoint')
        save(label,dict(fit=fit,seconds=seconds));results.append((fit,seconds))
    left,right=results[0][0],results[1][0]
    differences=dict(LL=abs(left['loglik']-right['loglik']),parameters=float(np.max(np.abs(left['parameters']-right['parameters']))),
        gradient=float(np.max(np.abs(left['gradient']-right['gradient']))),hessian=float(np.max(np.abs(left['hessian']-right['hessian']))),
        kappa=abs(float(left['parameters'][16]-right['parameters'][16])))
    exact=all(np.asarray(left[k]).tobytes() == np.asarray(right[k]).tobytes() for k in ('parameters','loglik','gradient','hessian','converged'))
    passed=(differences['LL'] <= 1e-6 and differences['kappa'] <= 1e-4
        and np.allclose(left['gradient'],right['gradient'],atol=1e-8,rtol=1e-6)
        and left['converged'] == right['converged'])
    require(p.tobytes() == before and repr(rows) == row_before and original.bb_terms.__name__ == 'bb_terms','Inputs/arithmetic changed')
    for path,digest in pins.items():require(sha(Path(path)) == digest,'Input changed during check')
    save('summary',dict(completed=True,draw=a.draw,endpoint_agreement_passed=bool(passed),all_endpoint_fields_byte_exact=exact,
        differences=differences,serial_seconds=results[0][1],parallel_seconds=results[1][1],
        serial_converged=left['converged'],parallel_converged=right['converged'],
        one_paired_scheduling_check_not_full_fit_speed=True,no_CRT=True,no_calibration=True,no_human_validation=True))
    if not passed:raise SystemExit(1)


if __name__ == '__main__':main()
