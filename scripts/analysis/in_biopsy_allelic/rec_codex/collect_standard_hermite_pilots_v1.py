#!/usr/bin/env python3
"""Read the two fixed synthetic pilots without replacing old draw results."""
import argparse
import ast
import hashlib
import json
import math
import os
from pathlib import Path
import struct
import subprocess

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
LAUNCH = REC/'standard_hermite_pilot_launch_20261007T103340Z'
ROSTER = ((0,'22084966'),(11,'22084967'))
RUNNER_SHA = 'c5810b4e264f912b114343eb02a2cfd1a3b986a30f7af767015c29681fa7a03a'
PINS_SHA = 'ef15156306adf5ecd6fc34e71f7a9bd3d230353b45cd0a00e57cd8d650b38513'


def require(ok, message):
    if not ok: raise RuntimeError(message)


def finite(x): return type(x) in (int,float) and math.isfinite(x)


def project(p, g, bounds, tolerance):
    require(len(p) == len(g) == len(bounds) and all(finite(x) for x in (*p,*g)),'Invalid point/gradient/bounds')
    return [0. if ((lo is not None and lo == hi) or (lo is not None and x<=lo+tolerance and y<0)
        or (hi is not None and x>=hi-tolerance and y>0)) else y for x,y,(lo,hi) in zip(p,g,bounds)]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    require(a.out.resolve().parent == REC.resolve(),'New direct child output required')
    a.out.mkdir(exist_ok=False);inputs={};read_bytes=0
    def body(p):
        nonlocal read_bytes
        require(p.resolve().is_relative_to(ROOT.resolve()) and not p.is_symlink(),'Input escaped repository')
        size=p.stat().st_size
        require(size <= 32*1024**2 and read_bytes+size <= 128*1024**2,'Read cap exceeded')
        value=p.read_bytes();read_bytes+=len(value);inputs[str(p)]=hashlib.sha256(value).hexdigest();return value
    def read(p):
        def bad(x):raise ValueError('Nonfinite JSON constant: '+x)
        return json.loads(body(p),parse_constant=bad)
    pins=read(LAUNCH/'pins.json');runner=body(LAUNCH/'run_standard_hermite_pilot_v1.py')
    require(inputs[str(LAUNCH/'pins.json')] == PINS_SHA and hashlib.sha256(runner).hexdigest() == RUNNER_SHA,'Fixed pilot launch differs')
    sources=next(ast.literal_eval(n.value) for n in ast.parse(runner).body if isinstance(n,ast.Assign)
        and any(isinstance(t,ast.Name) and t.id == 'SOURCES' for t in n.targets))
    for path,digest in pins.items():require(hashlib.sha256(body(Path(path))).hexdigest() == digest,'Pinned input differs')
    status=subprocess.check_output(['sacct','-j',','.join(j for _,j in ROSTER),'-n','-P',
        '--format=JobIDRaw,State,ExitCode,Elapsed,MaxRSS'],text=True)
    (a.out/'sacct.txt').write_text(status)
    states={row[0]:row[1:] for line in status.splitlines() if (row:=line.split('|')) and row[0]}
    records=[];issues=[]
    for draw,job in ROSTER:
        receipt=read(LAUNCH/f'draw{draw}_submission.json')
        require(receipt['draw'] == draw and receipt['jobid'] == job and receipt['runner_sha256'] == RUNNER_SHA
            and receipt['pins_sha256'] == PINS_SHA,'Pilot roster differs')
        state=states.get(job,['MISSING']);d=REC/f'standard_hermite_pilot_draw{draw}_{job}'
        record=dict(draw=draw,jobid=job,state=state,batch_state=states.get(job+'.batch'),completed=False)
        records.append(record)
        if state[0] not in ('RUNNING','PENDING','COMPLETED'):issues.append(dict(draw=draw,scheduler_state=state))
        if not d.exists():continue
        for n in sources:require(hashlib.sha256(body(d/'executed_sources'/n)).hexdigest() == pins[str(HERE/n)],'Executed helper differs')
        require(hashlib.sha256(body(d/'executed_sources/run_standard_hermite_pilot_v1.py')).hexdigest() == RUNNER_SHA,'Executed runner differs')
        protocol=read(d/'protocol.json')
        require(protocol['draw'] == draw and protocol['B'] == 50 and protocol['source_and_input_sha256'] == pins
            and protocol['original_draw_body_unchanged'] and protocol['original_generating_null_parameters']
            and protocol['original_starts_and_tie'] and protocol['original_convergence_cutoff'] == 1e-6
            and protocol['no_existing_draw_reclassification'] and protocol['no_backend_fallback']
            and protocol['no_calibration'] and protocol['no_human_data'],'Pilot protocol differs')
        if (d/'preflight.json').exists():
            pre=read(d/'preflight.json');require(pre['full_reference_differences'] == [0.,0.,0.]
                and pre['score_and_Hessian_modes_same_LL_gradient'] and pre['original_zero_boundary_exact'],'Preflight failed')
            record['preflight_passed']=True
        if (d/'draw_counts.json').exists():
            counts=read(d/'draw_counts.json');vector=counts['counts']
            require(len(vector) == 173 and all(type(x) is int for x in vector)
                and hashlib.sha256(struct.pack('<173q',*vector)).hexdigest() == counts['count_vector_sha256']
                and counts['exact_previous_counts'],'Counts/hash differ')
            old=REC/('repair_p3_gamma_order_polish_draw_0_amd0_22084420' if draw == 0 else
                     'repair_p3_gamma_order_polish_draw_11_failed11_22084421')/'draw_counts.json'
            previous=read(old);require(vector == previous['counts'] and counts['full_row_keys'] == previous['full_row_keys'],'Old draw counts differ')
            record['original_counts_verified']=True
        if (d/'progress.jsonl').exists():
            data=body(d/'progress.jsonl');lines=data.splitlines()
            if lines and not data.endswith(b'\n'):lines.pop()
            events=[json.loads(x) for x in lines]
            record['last_complete_progress']=events[-1] if events else None
            record['saved_positive_calls']=sum(e['event'] == 'marginal_evaluation' for e in events)
        if (d/'failure.json').exists():record['failure']=read(d/'failure.json')
        if state[:2] != ['COMPLETED','0:0'] or not (d/'summary.json').exists():
            if state[0] == 'COMPLETED':issues.append(dict(draw=draw,error='Completed job lacks successful complete output'))
            continue
        summary=read(d/'summary.json');raw=read(d/f'p3_{draw:03d}.json');fit=raw['fit']
        require(summary['completed'] and summary['draw'] == raw['b'] == draw,'Completed identity differs')
        failed=fit is None or not (fit['converged'] and fit['null']['converged'])
        require(summary['failed_draw'] == failed and summary['LR'] == raw['LR'] and summary['exceed'] == raw['exceed'] is True,'Original failure/exceedance differs')
        require(all(e['known_numerical_refusal'] for e in summary['caught_exceptions']),'Unqualified caught exception')
        require((fit is None) == bool(summary['caught_exceptions']),'Unexplained missing fit')
        record.update(completed=True,failed_draw=failed,LR=summary['LR'],seconds=summary['seconds'])
        if failed:require(summary['LR'] == {'nonfinite_float':'nan'},'Failed LR must remain NaN')
        else:
            lr=0. if fit['v'] == 0. else max(0.,2*(fit['loglik']-fit['null']['loglik']))
            require(summary['LR'] == lr,'LR arithmetic differs')
        if fit is not None:
            diag=read(d/'returned_diagnostic.json');point=fit['parameters']+[fit['v']]
            require(point == diag['parameters'] and len(point) == 29 and len(fit['null']['parameters']) == 28,'Returned coordinates differ')
            pg=project(point,diag['gradient'],diag['bounds'],1e-10);maximum=max(map(abs,pg))
            null_pg=max(map(abs,project(fit['null']['parameters'],fit['null']['gradient'],diag['bounds'][:-1],1e-8)))
            require(maximum == diag['max_projected_gradient'],'Saved marginal score projection differs')
            null=fit['null'];expected='cold' if null['start_logliks']['cold']-null['start_logliks']['warm'] >= 1e-6 else 'warm'
            require(null['kept_start'] == expected,'Null tie rule differs')
            agreement=abs(diag['loglik']-fit['loglik']) <= 1e-6 and all(abs(x-y) <= 1e-8+1e-6*abs(y) for x,y in zip(diag['gradient'],fit['gradient']))
            record.update(marginal_max_projected_gradient=maximum,null_max_projected_gradient=null_pg,
                LL_gradient_recomputation_agrees=agreement,null_kept_start=expected,null_start_converged=null['start_converged'],
                marginal_start_converged=fit['start_converged'])
            if not failed and not (maximum < 1e-6 and null_pg < 1e-6 and agreement):issues.append(dict(draw=draw,error='Successful recipe fails fresh score/LL checks'))
    completed=sum(r['completed'] for r in records);failed=sum(r.get('failed_draw',False) for r in records)
    result=dict(completed=completed,required=2,execution_complete=completed == 2,converged_pilots=completed-failed,
        prescribed_failures=failed,all_pilots_qualified=(not issues and failed == 0) if completed == 2 else None,
        records=records,issues=issues,input_sha256=inputs,no_fifty_draw_agreement=True,no_calibration=True,no_human_validation=True)
    with (a.out/'collection.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False);f.write('\n');f.flush();os.fsync(f.fileno())
    print(json.dumps({k:v for k,v in result.items() if k not in ('records','issues','input_sha256')}))
    if issues or failed:raise SystemExit(1)


if __name__ == '__main__':main()
