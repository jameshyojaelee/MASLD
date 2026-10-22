#!/usr/bin/env python3
"""Compare three bounded full-fit polish trials, never substitute fifty draws."""
import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import struct
import subprocess

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
LAUNCH=REC/'repair_p3_gamma_order_polish_launch_20261007T084132Z'
PINS_SHA='66e8d8f261591314ccda9de7384c842bc3fec0de460a6688ccc76383786266f0'
RUNNER_SHA='1a9f1e307d4fff6c2e28f3148b8d531e25d2586e63444e8abca415901a80e38c'
JOBS=(('intel0',0,'22084419'),('amd0',0,'22084420'),('failed11',11,'22084421'))


def require(ok,message):
    if not ok:raise RuntimeError(message)


def finite(x):return type(x) in (int,float) and math.isfinite(x)


def projected(p,g,bounds,tolerance):
    require(len(p)==len(g)==len(bounds) and all(finite(x) for x in (*p,*g)),'Invalid point/score')
    return [0. if (lo is not None and lo==hi) or (lo is not None and x<=lo+tolerance and s<0)
        or (hi is not None and x>=hi-tolerance and s>0) else s for x,s,(lo,hi) in zip(p,g,bounds)]


def compare(fit,lr,reference,reference_lr):
    require(len(fit['gradient'])==len(reference['gradient'])==29,'Comparison gradient dimension differs')
    g_ok=all(abs(a-b)<=1e-8+1e-6*abs(b) for a,b in zip(fit['gradient'],reference['gradient']))
    ng_ok=all(abs(a-b)<=1e-8+1e-6*abs(b) for a,b in zip(fit['null']['gradient'],reference['null']['gradient']))
    ll=fit['loglik']-reference['loglik'];nl=fit['null']['loglik']-reference['null']['loglik']
    kappa=fit['parameters'][16]-reference['parameters'][16]
    tau=math.sqrt(fit['v'])-math.sqrt(reference['v']);delta_lr=lr-reference_lr
    passed=(all(finite(x) for x in (ll,nl,kappa,tau,delta_lr)) and g_ok and ng_ok
        and abs(ll)<=1e-6 and abs(nl)<=1e-6 and abs(kappa)<=1e-4
        and (abs(tau)<=1e-3 if reference_lr>1 else abs(delta_lr)<=1e-6))
    return dict(passed=passed,LL_difference=ll,null_LL_difference=nl,kappa_difference=kappa,
        tau_difference=tau,LR_difference=delta_lr,gradient_agreement=g_ok,null_gradient_agreement=ng_ok,
        gradient_max_difference=max(abs(a-b) for a,b in zip(fit['gradient'],reference['gradient'])),
        nuisance_max_difference=max(abs(a-b) for a,b in zip(fit['parameters'],reference['parameters'])),
        gradient_atol=1e-8,gradient_rtol=1e-6,LL_atol=1e-6,kappa_atol=1e-4)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True);a=ap.parse_args()
    require(a.out.resolve().parent==REC.resolve(),'New direct child of REC required');a.out.mkdir(exist_ok=False)
    digests={};bytes_read=0
    def body(p):
        nonlocal bytes_read
        require(p.resolve().is_relative_to(ROOT.resolve()) and not p.is_symlink(),'Input escaped repository')
        size=p.stat().st_size;require(size<=32*1024**2 and bytes_read+size<=256*1024**2,'Read cap exceeded')
        b=p.read_bytes();bytes_read+=len(b);digests[str(p)]=hashlib.sha256(b).hexdigest();return b
    def read(p):
        def invalid(x):raise ValueError('Nonfinite JSON literal: '+x)
        return json.loads(body(p),parse_constant=invalid)
    pins=read(LAUNCH/'pins.json');require(digests[str(LAUNCH/'pins.json')]==PINS_SHA,'Input pins differ')
    sources=[Path(p) for p in pins if Path(p).parent==HERE and Path(p).suffix=='.py']
    require(len(sources)==20,'Source family differs')
    for p,h in pins.items():require(hashlib.sha256(body(Path(p))).hexdigest()==h,'Pinned input changed')
    rows_file=next(Path(p) for p in pins if Path(p).name=='rows.tsv')
    row_keys=[[r['individual_id'],r['gene_id']] for r in csv.DictReader(io.StringIO(body(rows_file).decode()),delimiter='\t')]
    require(len(row_keys)==173,'Fixture row order differs')
    runner=LAUNCH/'repair_p3_gamma_order_polish_draw.py';require(hashlib.sha256(body(runner)).hexdigest()==RUNNER_SHA,'Runner differs')
    state_text=subprocess.check_output(['sacct','-j',','.join(j for _,_,j in JOBS),'-n','-P',
        '--format=JobIDRaw,State,ExitCode,Elapsed,MaxRSS'],text=True);(a.out/'sacct.txt').write_text(state_text)
    steps={p[0]:p[1:] for line in state_text.splitlines() if (p:=line.split('|'))[0]}
    states={job:state for job,state in steps.items() if job.isdigit()}
    records=[];issues=[];completed={}
    proof=read(REC/'repair_p3_gamma_order_draw_0_22083807/p3_000.json')
    for trial,draw,job in JOBS:
        directory=REC/f'repair_p3_gamma_order_polish_draw_{draw}_{trial}_{job}';state=states.get(job,['MISSING'])
        r=dict(trial=trial,draw=draw,jobid=job,state=state,
            batch_MaxRSS=steps.get(job+'.batch',[None]*4)[3],completed=False);records.append(r)
        submission=read(LAUNCH/(trial+'_submission.json'));require(submission['job_id']==job and submission['draw']==draw,'Submission differs')
        if state[0] not in ('RUNNING','PENDING','CONFIGURING','COMPLETING','SUSPENDED','COMPLETED'):
            issues.append(dict(trial=trial,error='scheduler_not_successful',state=state))
        if not directory.exists():continue
        for source in sources:
            require(hashlib.sha256(body(directory/'executed_sources'/source.name)).hexdigest()==pins[str(source)],'Executed helper differs')
        require(hashlib.sha256(body(directory/'executed_sources'/runner.name)).hexdigest()==RUNNER_SHA,'Executed runner differs')
        protocol=read(directory/'protocol.json');claim=read(LAUNCH/f'draw_{draw:03d}_{trial}.claim.json')
        require(protocol['trial']==trial and protocol['draw']==draw and protocol['source_and_input_sha256']==pins
            and protocol['backend']=='adaptive_order_GL_with_null_score_integral_polish'
            and protocol['separate_bounded_trial'] and protocol['existing_fifty_draws_unchanged']
            and protocol['original_draw_body_unchanged'] and protocol['section4_from_null_unchanged']
            and protocol['no_held_results_read'] and protocol['no_calibration'],'Executed protocol differs')
        require(claim==dict(draw=draw,jobid=job,output=str(directory)),'Claim differs')
        progress=body(directory/'progress.jsonl');events=[json.loads(x) for x in progress.split(b'\n')[:-1]]
        r['last_saved_progress']=events[-1] if events else None
        r['saved_positive_calls']=sum(e['event']=='marginal_evaluation' for e in events)
        r['saved_null_polish_calls']=sum(e['event']=='null_score_integral_polish' for e in events)
        if state[:2]!=['COMPLETED','0:0']:
            if (directory/'failure.json').is_file():r['execution_failure']=read(directory/'failure.json')
            if state[0]=='COMPLETED':issues.append(dict(trial=trial,error='terminal_nonzero_exit'))
            continue
        required=('summary.json','recipe_summary.json',f'p3_{draw:03d}.json','draw_counts.json','null_polish_evidence.json')
        missing=[n for n in required if not (directory/n).is_file()]
        if missing:issues.append(dict(trial=trial,error='terminal_output_missing',missing=missing));continue
        summary=read(directory/'summary.json');recipe=read(directory/'recipe_summary.json');raw=read(directory/f'p3_{draw:03d}.json')
        require(summary['completed'] and recipe['completed_original_draw_body'] and raw['b']==draw,'Recipe identity differs')
        fit=raw['fit'];failed=fit is None or not (fit['converged'] and fit['null']['converged'])
        require(failed==summary['failed_draw']==recipe['original_recipe_failed_draw'],'Failure accounting differs')
        require(summary['caught_exceptions']==recipe['caught_exceptions'] and not summary['unknown_caught_exception_present'],
            'Unknown exception in completed draw')
        counts=read(directory/'draw_counts.json');vector=counts['counts']
        require(counts['full_row_keys']==row_keys,'Count row keys differ')
        require(len(vector)==173 and all(type(x) is int for x in vector),'Count vector invalid')
        require(hashlib.sha256(struct.pack('<173q',*vector)).hexdigest()==counts['count_vector_sha256']
            ==counts['original_baseline_count_vector_sha256'] and counts['exact_original_baseline_counts_equal'],'Counts differ')
        evidence=read(directory/'null_polish_evidence.json');accepted=0
        for call in evidence['trials']:
            for it in call['history']:
                for t in it['trials']:
                    if not t.get('accepted'):continue
                    accepted+=1
                    require(all(finite(t[k]) and t[k]<(1-.5*t['damping'])*it[k]
                        for k in ('Gamma_max_projected_gradient','original_max_projected_gradient')),'Accepted score decrease differs')
                    proxy=abs(t['original_score_integral32']-t['score_integral32'])
                    require(proxy==t['arithmetic_noise_proxy'] and min(t['empirical_gain_lower_estimate'],
                        t['original_empirical_gain_lower_estimate'])>32*proxy,'Accepted empirical gain check differs')
        failure=summary['fitting_integration_failure']
        r.update(completed=True,failed_draw=failed,accepted_null_polish_steps=accepted,
            caught_numerical_exceptions=summary['caught_exceptions'],
            fitting_integration_failure_summary={k:failure[k] for k in ('call','gene','failure_phase','message','variance')
                if k in failure} if failure is not None else None)
        if failed:issues.append(dict(trial=trial,error='prescribed_recipe_failed'));continue
        if summary['diagnostic_failure'] is not None or not summary['diagnostic_completed']:
            issues.append(dict(trial=trial,error='post_fit_diagnostic_failed'))
            continue
        if not (directory/'GL_fit_diagnostic.json').is_file():
            issues.append(dict(trial=trial,error='final_diagnostic_missing'));continue
        lr=0. if fit['v']==0 else max(0.,2*(fit['loglik']-fit['null']['loglik']))
        require(raw['LR']==summary['LR']==recipe['LR']==lr and raw['exceed'] is True,'Original LR arithmetic differs')
        d=read(directory/'GL_fit_diagnostic.json');point=fit['parameters']+[fit['v']]
        require(d['parameters']==point and len(point)==len(d['gradient'])==len(d['bounds'])==len(fit['gradient'])==29,
            'Diagnostic identity differs')
        pg=projected(point,d['gradient'],d['bounds'],1e-10)
        ng=projected(fit['null']['parameters'],fit['null']['gradient'],d['bounds'][:-1],1e-8)
        require(pg==d['projected_gradient'] and max(map(abs,pg))==d['max_projected_gradient'],'Projected diagnostic differs')
        agreement=abs(d['loglik']-fit['loglik'])<=1e-8+1e-8*abs(fit['loglik']) and all(
            abs(x-y)<=1e-8+1e-6*abs(y) for x,y in zip(d['gradient'],fit['gradient']))
        score_ok=max(map(abs,pg))<1e-6 and max(map(abs,ng))<1e-6 and agreement
        r.update(LR=lr,variance=fit['v'],max_projected_gradient=max(map(abs,pg)),
            null_max_projected_gradient=max(map(abs,ng)),null_gradient_saved_not_recomputed=True,score_checks_passed=score_ok)
        if not score_ok:issues.append(dict(trial=trial,error='final_score_check_failed'))
        if fit['v']>0:
            e=d['positive_integration_error_allowances'];ok=(e is not None and finite(e['logmass']) and 0<=e['logmass']<=1e-8
                and len(e['gradient'])==29 and all(finite(x) and 0<=x<=1e-7 for x in e['gradient']))
            r['positive_empirical_allowances_passed']=ok
            if not ok:issues.append(dict(trial=trial,error='positive_allowances_unresolved'))
        if draw==0:
            c=compare(fit,lr,proof['fit'],proof['LR']);r['comparison_to_prior_converged_draw0']=c
            if not c['passed']:issues.append(dict(trial=trial,error='prior_draw0_agreement_failed'))
        completed[trial]=(fit,lr)
    cross=None
    if all(k in completed for k in ('intel0','amd0')):
        cross=compare(*completed['amd0'],*completed['intel0'])
        if not cross['passed']:issues.append(dict(error='cross_architecture_draw0_agreement_failed'))
    exercised=any(r.get('accepted_null_polish_steps',0)>0 for r in records)
    if len(completed)==3 and not exercised:issues.append(dict(error='polish_never_exercised'))
    report=dict(completed_draws=sum(r['completed'] for r in records),required_draws=3,
        converged_full_fits=len(completed),
        all_three_full_fits_passed=(len(completed)==3 and not issues),records=records,issues=issues,
        cross_architecture_draw0=cross,polish_exercised=exercised,
        collector_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        input_sha256=digests,bytes_read=bytes_read,
        empirical_gain_sign_not_certified=True,no_fifty_draw_agreement=True,no_calibration=True,no_human_validation=True)
    (a.out/'collection.json').write_text(json.dumps(report,indent=2,allow_nan=False)+'\n')
    print(json.dumps(dict(completed=report['completed_draws'],required=3,issues=len(issues),all_three_full_fits_passed=report['all_three_full_fits_passed'])))
    if issues:raise SystemExit(1)


if __name__=='__main__':main()
