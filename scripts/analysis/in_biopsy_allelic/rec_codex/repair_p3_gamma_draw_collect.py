#!/usr/bin/env python3
"""Collect the fixed fifty native-Gamma P3 draws without held-result access.

Recompute original statistics and separate diagnostic eligibility. Never infer
completion from a partial fit or turn a numerical refusal into a returned fit.
Only bounded JSON, source bytes and live scheduler state are inspected.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
LAUNCH=REC/'repair_p3_gamma_draw_v3_launch_20261007T034152Z'
OBSERVED=REC/'repair_p3_gamma_observed_native_22080548'
RUNNER_SHA='9a3c978dccc685cf648afa35e30c419a4b4b2eb9944a6391c4fb1480b2d16b81'
PINS_SHA='56a131cdd0ed5c1f77657a94055ba082ac2d0bd1360907ec710073c107251c9c'
CAP=128*1024**2
SOURCES=('likelihood_v1.py','design_v1.py','bfix_v1.py','tests_v1.py',
    'check_bfix_tests_v13.py','repair_p3_bb_gamma_difference_v3.py',
    'repair_p3_gamma_parallel_model.py','repair_p3_gamma_native_model.py')


def require(ok,message):
    if not ok:raise ValueError(message)


def finite(x):
    return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    require(args.out.resolve().is_relative_to(REC),'Output escaped synthetic execution tree')
    args.out.mkdir(exist_ok=False);bytes_read=0;digests={}
    def body(p):
        nonlocal bytes_read
        require(p.resolve().is_relative_to(ROOT.resolve()) and not p.is_symlink(),'Input escaped repository')
        size=p.stat().st_size;require(size<=8*1024**2 and bytes_read+size<=CAP,'Input byte cap exceeded')
        value=p.read_bytes();bytes_read+=len(value);require(len(value)==size,'Input changed size during collection')
        digests[str(p)]=hashlib.sha256(value).hexdigest();return value
    def read(p):
        def invalid(x):raise ValueError('Nonfinite JSON constant: '+x)
        return json.loads(body(p),parse_constant=invalid)
    try:
        jobs=read(LAUNCH/'jobs.json');pins=read(LAUNCH/'pins.json')
        require(len(jobs)==50 and [r['draw'] for r in jobs]==list(range(50))
            and [r['jobid'] for r in jobs]==[str(n) for n in range(22081569,22081619)],'Fixed fifty-draw roster differs')
        require(digests[str(LAUNCH/'pins.json')]==PINS_SHA,'Declared input family differs')
        observed_protocol=read(OBSERVED/'protocol.json')
        require(digests[str(OBSERVED/'protocol.json')]==pins[str(OBSERVED/'protocol.json')],
                'Observed-fit numerical environment reference differs')
        require(hashlib.sha256(body(LAUNCH/'repair_p3_gamma_draw.py')).hexdigest()==RUNNER_SHA,'Executed runner differs')
        for p,h in pins.items():require(hashlib.sha256(body(Path(p))).hexdigest()==h,'Pinned source/input differs')
        status=subprocess.check_output(['sacct','-j',','.join(r['jobid'] for r in jobs),
            '--noheader','--parsable2','--format=JobIDRaw,State,ExitCode,Elapsed,MaxRSS'],text=True)
        (args.out/'sacct.txt').write_text(status)
        states={parts[0]:parts[1:] for line in status.splitlines() if (parts:=line.split('|'))[0].isdigit()}
        records=[];incomplete=[];issues=[]
        for job in jobs:
            draw=job['draw'];directory=Path(job['output_directory'])
            require(directory==REC/f"repair_p3_gamma_draw_{draw}_{job['jobid']}",'Output identity differs')
            state=states.get(job['jobid'],['MISSING','MISSING'])
            if state[:2]!=['COMPLETED','0:0']:
                record=dict(draw=draw,jobid=job['jobid'],state=state)
                for n in ('failure','recipe_summary','diagnostic_failure'):
                    if (directory/(n+'.json')).exists():record[n]=read(directory/(n+'.json'))
                incomplete.append(record);continue
            missing=[n for n in ('protocol.json','summary.json','recipe_summary.json',
                f'p3_{draw:03d}.json','draw_counts.json') if not (directory/n).is_file()]
            if missing:
                incomplete.append(dict(draw=draw,jobid=job['jobid'],state=state,missing_terminal_files=missing));continue
            protocol=read(directory/'protocol.json');summary=read(directory/'summary.json')
            recipe=read(directory/'recipe_summary.json');raw=read(directory/f'p3_{draw:03d}.json')
            counts=read(directory/'draw_counts.json');claim=read(LAUNCH/f'draw_{draw:03d}.claim.json')
            require(protocol['draw']==summary['draw']==recipe['draw']==raw['b']==claim['draw']==draw
                and claim['jobid']==job['jobid'] and claim['output']==str(directory),'Draw identity differs')
            require(protocol['B']==summary['B']==recipe['B']==50 and protocol['observed_LR']==0.
                and protocol['source_and_input_sha256']==pins,'Draw protocol differs')
            for name in SOURCES:
                h=pins[str(ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'/name)]
                require(hashlib.sha256(body(directory/'executed_sources'/name)).hexdigest()==h,
                    'Executed numerical source differs')
            require(hashlib.sha256(body(directory/'executed_sources/repair_p3_gamma_draw.py')).hexdigest()==RUNNER_SHA,
                    'Executed draw runner differs')
            require(summary['completed'] is True and recipe['completed_original_draw_body'] is True,'Completed draw lacks recipe')
            fit=raw['fit'];caught=recipe['caught_exceptions']
            require((fit is None)==bool(caught) and summary['caught_exceptions']==caught,'Refusal record differs')
            failed=fit is None or not (fit['converged'] and fit['null']['converged'])
            require(failed==summary['failed_draw']==recipe['original_recipe_failed_draw'],'Original failure count differs')
            require(recipe['unknown_caught_exception_present']==any(not r['known_numerical_refusal'] for r in caught),
                    'Unknown exception classification differs')
            if failed:
                require(raw['LR']==recipe['LR']==summary['LR']=={'nonfinite_float':'nan'},'Failed draw lost its unavailable LR')
            else:
                require(finite(fit['v']) and fit['v']>=0 and finite(fit['loglik']) and finite(fit['null']['loglik']),
                    'Converged draw has invalid parameters/LL')
                lr=0. if fit['v']==0. else max(0.,2*(fit['loglik']-fit['null']['loglik']))
                require(raw['LR']==recipe['LR']==summary['LR']==lr,'Original LR arithmetic differs')
            require(raw['exceed'] is True and summary['exceed'] is True and recipe['exceed'] is True,
                    'Original observed-LR-zero exceedance differs')
            if counts['counts'] is not None:
                vector=counts['counts'];require(len(vector)==173 and all(type(x) is int for x in vector),'Draw count vector invalid')
                digest=hashlib.sha256(struct.pack('<173q',*vector)).hexdigest()
                require(digest==counts['count_vector_sha256']==counts['original_baseline_count_vector_sha256']
                    and counts['exact_original_baseline_counts_equal'] is True,'Draw differs from original baseline counts')
            elif fit is not None:raise ValueError('Returned fit lacks exact count vector')
            diagnostics=None;diagnostic_pass=False;null_maximum=None;forty_status='not_available'
            if fit is not None and (directory/'returned_fit_diagnostic.json').exists():
                diagnostics=read(directory/'returned_fit_diagnostic.json');p=diagnostics['parameters'];g=diagnostics['gradient']
                bounds=diagnostics['bounds'];require(len(p)==len(g)==len(bounds)==29
                    and all(finite(x) for x in (*p,*g)),'Invalid saved projected-score diagnostic')
                require(p==fit['parameters']+[fit['v']],'Diagnostic parameters differ from fitted LR point')
                require(finite(diagnostics['loglik']) and finite(fit['loglik']),'Invalid recomputed fit likelihood')
                require(diagnostics['native_converged']==fit['converged'],'Native convergence flag differs')
                fit_LL_agree=abs(diagnostics['loglik']-fit['loglik'])<=1e-8+1e-8*abs(fit['loglik'])
                saved_g=fit['gradient']
                fit_score_agree=(len(saved_g)==29 and all(finite(x) for x in saved_g)
                    and all(abs(a-b)<=1e-8+1e-6*abs(b) for a,b in zip(g,saved_g)))
                projected=[0. if ((lo is not None and lo==hi) or (lo is not None and x<=lo+1e-10 and score<0)
                    or (hi is not None and x>=hi-1e-10 and score>0)) else score
                    for x,score,(lo,hi) in zip(p,g,bounds)]
                maximum=max(map(abs,projected));require(projected==diagnostics['projected_gradient']
                    and maximum==diagnostics['max_projected_gradient'],'Projected gradient differs')
                nullp=fit['null']['parameters'];nullg=fit['null']['gradient']
                if len(nullp)==len(nullg)==28 and all(finite(x) for x in (*nullp,*nullg)):
                    null_projected=[0. if ((lo is not None and lo==hi) or (lo is not None and x<=lo+1e-8 and score<0)
                        or (hi is not None and x>=hi-1e-8 and score>0)) else score
                        for x,score,(lo,hi) in zip(nullp,nullg,bounds[:-1])]
                    null_maximum=max(map(abs,null_projected))
                forty=diagnostics['native40'];order_pass=False
                if forty['available']:
                    forty_finite=finite(forty['loglik']) and len(forty['gradient'])==29 and all(finite(x) for x in forty['gradient'])
                    if forty_finite:
                        dl=forty['loglik']-diagnostics['loglik'];dg=max(abs(a-b) for a,b in zip(forty['gradient'],g))
                        require(dl==forty['LL_difference'] and dg==forty['score_max_difference'],'Native40 difference differs')
                        order_pass=(abs(dl)<=1e-8+1e-8*abs(forty['loglik'])
                            and all(abs(a-b)<=1e-8+1e-6*abs(b) for a,b in zip(forty['gradient'],g)))
                    require(order_pass==forty['LL_score_agree'] and forty['order_refinement_applicable']==(fit['v']>0),
                            'Native40 criterion or boundary applicability differs')
                    forty_status='returned_finite' if forty_finite else 'returned_nonfinite'
                elif fit['v']>0:forty_status='refused'
                if fit['v']==0:forty_status='not_applicable_zero_boundary'
                diagnostic_pass=(maximum<1e-6 and null_maximum is not None and null_maximum<1e-6
                    and fit_LL_agree and fit_score_agree and (fit['v']==0 or order_pass))
            local=[]
            if recipe['unknown_caught_exception_present']:local.append('unknown_original_caught_exception')
            if not failed and not diagnostic_pass:local.append('returned_fit_numerical_diagnostic_not_passed')
            if summary['diagnostic_failure'] is not None:local.append('post_recipe_diagnostic_failure')
            if local:issues.append(dict(draw=draw,issues=local))
            records.append(dict(draw=draw,jobid=job['jobid'],failed_draw=failed,LR=raw['LR'],exceed=True,
                recipe=recipe,summary=summary,diagnostic_pass=diagnostic_pass,
                null_max_projected_gradient_from_saved_score=null_maximum,native40_status=forty_status,issues=local))
        complete=len(records)==50
        report=dict(all50_jobs_completed_successfully=complete,completed_draws=len(records),required_draws=50,
            environment=observed_protocol['environment'],fixture_participants=60,synthetic_fixture=True,
            variance_coordinate='v=tau**2',
            scheduler_states=dict(Counter(states.get(r['jobid'],['MISSING'])[0] for r in jobs)),
            incomplete_jobs=incomplete,draw_records=records,numerical_or_integrity_issues=issues,
            no_independent_fifty_draw_agreement=True,no_calibration=True,no_human_validation=True)
        if complete:
            p=(1+sum(r['exceed'] for r in records))/51
            report['original_recipe_P3']=dict(B=50,p=p,mc_se=math.sqrt(p*(1-p)/50),observed_LR=0.,
                failed=[r['draw'] for r in records if r['failed_draw']],draw_LR=[r['LR'] for r in records],
                exceedance=[r['exceed'] for r in records],note='Synthetic agreement fixture, not calibration.')
            report['all_recipe_successful_fit_diagnostics_passed']=all(r['diagnostic_pass'] for r in records if not r['failed_draw'])
            report['no_unknown_caught_exceptions']=not any('unknown_original_caught_exception' in r['issues'] for r in records)
        report['input_sha256']=digests;report['bytes_read']=bytes_read
        with (args.out/'collection.json').open('x') as f:json.dump(report,f,allow_nan=False,indent=2);f.write('\n')
        print(json.dumps(dict(complete=complete,completed_draws=len(records),incomplete=len(incomplete),issues=len(issues))))
    except Exception:
        (args.out/'failure.txt').write_text(traceback.format_exc());raise


if __name__=='__main__':main()
