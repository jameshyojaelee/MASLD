#!/usr/bin/env python3
"""Diagnose order refusal from a saved converged draw0 near-null point.

This is not a replay of the lost pilot null endpoint and not a full P3 refit.
It runs the original marginal single-start method without repeating a null fit.
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
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXTURE=ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-a-bfix-20260930T115521Z'
NULL=REC/'repair_p3_gamma_order_draw_0_22083807/p3_000.json'
COUNTS=REC/'standard_hermite_pilot_draw0_22084966/draw_counts.json'
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py','repair_p3_local_moments.py',
         'repair_p3_pathwise_moments_v1.py','repair_p3_longdouble_conditional_v1.py',
         'repair_p3_standard_hermite_pathwise_v2.py','repair_p3_standard_hermite_model_v1.py')


def require(ok,message):
    if not ok:raise RuntimeError(message)


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--pins',type=Path,required=True)
    ap.add_argument('--expected-pins-sha',required=True);ap.add_argument('--expected-runner-sha',required=True)
    ap.add_argument('--inspect',action='store_true');a=ap.parse_args()
    require(sha(Path(__file__)) == a.expected_runner_sha and sha(a.pins) == a.expected_pins_sha,'Runner/pins differ')
    pins=json.loads(a.pins.read_text())
    for path,digest in pins.items():
        p=Path(path);require(p.resolve().is_relative_to(ROOT.resolve()) and sha(p) == digest,'Input differs')
        if p.suffix == '.py':ast.parse(p.read_text())
    require(all(str(p) in pins for p in [*(HERE/n for n in SOURCES),NULL,COUNTS,
        *(FIXTURE/n for n in ('globals.json','cohorts.tsv','individuals.tsv','rows.tsv'))]),'Missing source/input pin')
    read=lambda p:json.loads(p.read_text());null=read(NULL)['fit']['null'];counts=read(COUNTS)
    require(null['converged'] and len(null['parameters']) == 28 and len(counts['counts']) == 173,'Saved scope differs')
    if a.inspect:print(json.dumps(dict(saved_null_converged=True,single_start_only=True,no_full_refit=True)));return
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID')
        and os.environ.get('SLURM_CPUS_PER_TASK') == '8','Individual eight-CPU compute required')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS')),'Threads differ')
    out=REC/('standard_hermite_refusal_check_'+os.environ['SLURM_JOB_ID']);out.mkdir(exist_ok=False)
    archive=out/'executed_sources';archive.mkdir()
    for p in (Path(__file__),*(HERE/n for n in SOURCES)):
        q=archive/p.name;q.write_bytes(p.read_bytes())
        require(sha(q) == (a.expected_runner_sha if p == Path(__file__) else pins[str(p)]),'Archive differs')
    import numpy as np
    import scipy
    sys.path.insert(0,str(archive));modules={n[:-3]:importlib.import_module(n[:-3]) for n in SOURCES}
    require(all(Path(m.__file__).resolve().parent == archive.resolve() for m in modules.values()),'Import escaped archive')
    loader,original=modules['bfix_v1'],modules['likelihood_v1']
    require(loader.FIXTURE.resolve() == FIXTURE.resolve(),'Fixture differs')
    design,settings=loader.load();template=design.model()
    require(counts['full_row_keys'] == [[r.individual,settings['gene_order'][r.gene]] for r in template.rows],'Count order differs')
    rows=[replace(r,a=c) for r,c in zip(template.rows,counts['counts'])]
    started=time.perf_counter()
    def safe(x):
        if isinstance(x,np.ndarray):return safe(x.tolist())
        if isinstance(x,np.generic):return safe(x.item())
        if isinstance(x,float) and not math.isfinite(x):return {'nonfinite_float':x.hex()}
        if isinstance(x,dict):return {k:safe(v) for k,v in x.items()}
        if isinstance(x,(list,tuple)):return [safe(v) for v in x]
        return x
    def save(name,value):
        with (out/(name+'.json')).open('x') as f:json.dump(safe(value),f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    def event(fields):
        with (out/'progress.jsonl').open('a') as f:f.write(json.dumps(safe(dict(seconds=time.perf_counter()-started,**fields)))+'\n');f.flush()
    wrapper=modules['repair_p3_standard_hermite_model_v1'];Model=wrapper.model_class(original,
        modules['repair_p3_local_moments'],modules['repair_p3_pathwise_moments_v1'],
        modules['repair_p3_longdouble_conditional_v1'],modules['repair_p3_standard_hermite_pathwise_v2'],event_sink=event)
    model=Model(rows,gene_count=8);start=np.r_[null['parameters'],0.]
    save('protocol',dict(source_and_input_sha256=pins,start=start,rows=len(rows),
        environment=dict(python=sys.version,numpy=np.__version__,scipy=scipy.__version__),
        saved_near_null_point_not_actual_failed_pilot_endpoint=True,no_full_refit=True,no_CRT=True,no_calibration=True))
    try:
        fit=model.fit_marginal_single(start,None)
    except wrapper.HermiteRefusal as exc:
        failure=model.last_integration_failure;require(failure is not None,'Refusal point missing')
        save('refusal',dict(message=str(exc),failure=failure,traceback=traceback.format_exc()))
        reproduced=True
    else:
        save('single_start_fit',fit);reproduced=False
    for path,digest in pins.items():require(sha(Path(path)) == digest,'Input changed during diagnostic')
    save('summary',dict(completed=True,refusal_reproduced=reproduced,seconds=time.perf_counter()-started,
        positive_marginal_calls=getattr(model,'marginal_call_count',0),no_existing_draw_reclassification=True,
        no_full_refit=True,no_CRT=True,no_calibration=True,no_human_validation=True))


if __name__ == '__main__':main()
