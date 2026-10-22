#!/usr/bin/env python3
"""Read the two fixed GL pilots; never pool them with the native fifty draws."""
import argparse
import ast
from collections import Counter
import csv
import hashlib
import io
import json
import math
from pathlib import Path
import struct
import subprocess
import traceback

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
LAUNCH=REC/'repair_p3_gamma_GL_draw_pilot_v3_launch_20261007T045626Z'
RUNNER_SHA='66342d0caf36d3284ab515384adf3b43a063c6140d146b91112d39379b3edc38'
PINS_SHA='98b981f5f77709eff30ea69bfb5408a38ab643ec7b2d04ecd50bad93a086d5d0'


def require(ok,message):
    if not ok:raise ValueError(message)


def finite(x):
    return type(x) in (int,float) and math.isfinite(x)


def projected(p,g,bounds,tolerance):
    require(len(p)==len(g)==len(bounds) and all(finite(x) for x in (*p,*g)),
            'Invalid point/score/bounds dimensions')
    return [0. if ((lo is not None and lo==hi) or (lo is not None and x<=lo+tolerance and s<0)
        or (hi is not None and x>=hi-tolerance and s>0)) else s for x,s,(lo,hi) in zip(p,g,bounds)]


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    require(args.out.resolve().parent==REC.resolve(),'Output must be a new direct child of REC')
    args.out.mkdir(exist_ok=False);digests={};bytes_read=0
    def body(p):
        nonlocal bytes_read
        require(p.resolve().is_relative_to(ROOT.resolve()) and not p.is_symlink(),'Input escaped repository')
        size=p.stat().st_size;require(size<=32*1024**2 and bytes_read+size<=128*1024**2,'Read cap exceeded')
        with p.open('rb') as f:data=f.read(size)
        require(len(data)==size,'Input shrank while reading');bytes_read+=len(data)
        digests[str(p)]=hashlib.sha256(data).hexdigest();return data
    def read(p):
        def invalid(x):raise ValueError('Nonfinite JSON constant: '+x)
        return json.loads(body(p),parse_constant=invalid)
    try:
        jobs=read(LAUNCH/'jobs.json');pins=read(LAUNCH/'pins.json')
        require([(j['draw'],j['jobid']) for j in jobs]==[(0,'22082078'),(18,'22082079')],'Fixed pilot roster differs')
        require(digests[str(LAUNCH/'pins.json')]==PINS_SHA,'Pinned family differs')
        runner=body(LAUNCH/'repair_p3_gamma_GL_draw.py');require(hashlib.sha256(runner).hexdigest()==RUNNER_SHA,'Runner differs')
        tree=ast.parse(runner);sources=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign)
            and any(isinstance(t,ast.Name) and t.id=='SOURCES' for t in n.targets))
        require(len(sources)==len(set(sources))==16,'Executed source family differs')
        for path,digest in pins.items():require(hashlib.sha256(body(Path(path))).hexdigest()==digest,'Pinned input differs')
        row_file=next(Path(p) for p in pins if Path(p).name=='rows.tsv')
        row_keys=[[r['individual_id'],r['gene_id']] for r in csv.DictReader(io.StringIO(body(row_file).decode()),delimiter='\t')]
        require(len(row_keys)==173,'Supplied row-key count differs')
        status=subprocess.check_output(['sacct','-j',','.join(j['jobid'] for j in jobs),'--noheader','--parsable2',
            '--format=JobIDRaw,State,ExitCode,Elapsed,MaxRSS'],text=True)
        (args.out/'sacct.txt').write_text(status)
        states={p[0]:p[1:] for line in status.splitlines() if (p:=line.split('|'))[0].isdigit()}
        records=[];issues=[]
        for job in jobs:
            draw=job['draw'];directory=Path(job['output_directory']);state=states.get(job['jobid'],['MISSING'])
            require(directory==REC/f"repair_p3_gamma_GL_draw_{draw}_{job['jobid']}",'Output identity differs')
            record=dict(draw=draw,jobid=job['jobid'],state=state,completed_record=False)
            records.append(record)
            terminal_error=state[0].split()[0].rstrip('+') in ('FAILED','TIMEOUT','OUT_OF_MEMORY',
                'CANCELLED','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED')
            if terminal_error:issues.append(dict(draw=draw,issues=['scheduler_not_successful'],state=state))
            if state[0]=='MISSING':issues.append(dict(draw=draw,issues=['scheduler_state_missing']))
            if not directory.exists():
                record['output_directory_missing']=True
                if state[0]=='COMPLETED':issues.append(dict(draw=draw,issues=['terminal_output_directory_missing']))
                continue
            for name in sources:
                require(hashlib.sha256(body(directory/'executed_sources'/name)).hexdigest()==pins[str(HERE/name)],
                        'Archived helper differs')
            require(hashlib.sha256(body(directory/'executed_sources/repair_p3_gamma_GL_draw.py')).hexdigest()==RUNNER_SHA,
                    'Archived runner differs')
            if (directory/'protocol.json').exists():
                protocol=read(directory/'protocol.json');claim=read(LAUNCH/f'draw_{draw:03d}.claim.json')
                require(protocol['draw']==draw and protocol['B']==50 and protocol['backend']=='batched_GL'
                    and protocol['pilot_draws']==[0,18] and protocol['observed_LR']==0.
                    and protocol['source_and_input_sha256']==pins and protocol['no_backend_fallback'],
                    'Pilot protocol differs')
                require(claim==dict(draw=draw,jobid=job['jobid'],output=str(directory)),'Claim differs')
            progress=directory/'progress.jsonl'
            if progress.exists():
                data=body(progress);lines=data.splitlines();events=[]
                # A live final line may still be in flight; complete newline records only.
                for line in lines[:-1]:events.append(json.loads(line))
                if data.endswith(b'\n') and lines:events.append(json.loads(lines[-1]))
                completed=[e for e in events if e['event']=='marginal_evaluation']
                record['saved_completed_positive_evaluations_lower_bound']=len(completed)
                record['last_saved_progress']=events[-1] if events else None
            if terminal_error and (directory/'failure.json').exists():record['failure']=read(directory/'failure.json')
            required=['protocol.json','summary.json','recipe_summary.json',f'p3_{draw:03d}.json','draw_counts.json']
            missing=[n for n in required if not (directory/n).is_file()]
            if state[:2]!=['COMPLETED','0:0'] or missing:
                record['missing_files']=missing
                if state[0]=='COMPLETED':issues.append(dict(draw=draw,issues=['terminal_job_incomplete_or_unsuccessful'],missing=missing))
                continue
            summary=read(directory/'summary.json');recipe=read(directory/'recipe_summary.json')
            raw=read(directory/f'p3_{draw:03d}.json');counts=read(directory/'draw_counts.json')
            require(summary['completed'] and recipe['completed_original_draw_body']
                and summary['draw']==recipe['draw']==raw['b']==draw and summary['B']==recipe['B']==50,
                'Completed recipe identity differs')
            fit=raw['fit'];caught=recipe['caught_exceptions']
            require((fit is None)==bool(caught) and summary['caught_exceptions']==caught
                and not recipe['unknown_caught_exception_present'] and not summary['unknown_caught_exception_present']
                and all(e['known_numerical_refusal'] for e in caught),'Unqualified original exception')
            failed=fit is None or not (fit['converged'] and fit['null']['converged'])
            require(failed==summary['failed_draw']==recipe['original_recipe_failed_draw'],'Failure classification differs')
            if failed:lr={'nonfinite_float':'nan'}
            else:
                require(finite(fit['v']) and fit['v']>=0 and finite(fit['loglik']) and finite(fit['null']['loglik']),
                        'Converged point invalid')
                lr=0. if fit['v']==0. else max(0.,2*(fit['loglik']-fit['null']['loglik']))
            require(raw['LR']==recipe['LR']==summary['LR']==lr and raw['exceed'] is True
                and recipe['exceed'] is True and summary['exceed'] is True,'Original LR/exceedance arithmetic differs')
            require(counts['full_row_keys']==row_keys,'Draw row keys differ')
            vector=counts['counts']
            if vector is None:
                require(fit is None and caught and all(e['where'] in ('conditional_draw','draw_predictors',
                    'draw_model_constructor') for e in caught) and counts['count_vector_sha256'] is None
                    and counts['exact_original_baseline_counts_equal'] is None,'Unexplained unavailable draw counts')
                record['counts_unavailable_before_model_constructor_completed']=True
            else:
                require(len(vector)==173 and all(type(x) is int for x in vector),'Invalid draw counts')
                digest=hashlib.sha256(struct.pack('<173q',*vector)).hexdigest()
                require(digest==counts['count_vector_sha256']==counts['original_baseline_count_vector_sha256']
                    and counts['exact_original_baseline_counts_equal'],'Original count reconstruction differs')
            require(summary['fitting_integration_failure']==recipe['fitting_integration_failure'],
                    'Fitting refusal details differ')
            record.update(completed_record=True,prescribed_recipe_failed_under_GL=failed,LR=lr,caught_exceptions=caught,
                fitting_integration_failure=summary['fitting_integration_failure'],
                positive_marginal_call_count=summary['positive_marginal_call_count'],
                positive_marginal_call_count_includes_final_GL_diagnostic=True,score_checks_passed=False,
                native40_non_gating=True)
            diagnostic_path=directory/'GL_fit_diagnostic.json'
            if fit is not None and diagnostic_path.exists():
                d=read(diagnostic_path);p=d['parameters'];g=d['gradient'];bounds=d['bounds']
                require(len(p)==29 and p==fit['parameters']+[fit['v']] and d['native_converged']==fit['converged'],
                        'GL diagnostic differs from fitted LR point')
                pg=projected(p,g,bounds,1e-10);maximum=max(map(abs,pg))
                require(pg==d['projected_gradient'] and maximum==d['max_projected_gradient'],'GL projection differs')
                null=fit['null'];require(len(null['parameters'])==len(null['gradient'])==28,'Null dimension differs')
                null_max=max(map(abs,projected(null['parameters'],null['gradient'],bounds[:-1],1e-8)))
                require(finite(d['loglik']) and finite(fit['loglik']) and len(fit['gradient'])==29
                    and all(finite(x) for x in fit['gradient']),'Invalid saved LL/score')
                fit_agreement=abs(d['loglik']-fit['loglik'])<=1e-8+1e-8*abs(fit['loglik']) and all(
                    abs(a-b)<=1e-8+1e-6*abs(b) for a,b in zip(g,fit['gradient']))
                allowances=d['positive_integration_error_allowances']
                allowances_pass=None
                if fit['v']>0:
                    allowances_pass=(allowances is not None and finite(allowances['logmass']) and 0<=allowances['logmass']<=1e-8
                        and len(allowances['gradient'])==29 and all(finite(x) and 0<=x<=1e-7 for x in allowances['gradient']))
                record.update(GL_max_projected_gradient=maximum,null_max_projected_gradient=null_max,
                    null_score_is_saved_not_recomputed=True,null_projection_bound_tolerance=1e-8,
                    fit_LL_and_score_recomputation_agree=fit_agreement,
                    score_checks_passed=maximum<1e-6 and null_max<1e-6 and fit_agreement,
                    saved_GL_error_allowances=allowances if fit['v']>0 else None,
                    GL_allowance_check_applicable=fit['v']>0,GL_allowance_check_passed=allowances_pass,
                    allowances_not_independently_certified=True)
            comparison_path=directory/'returned_fit_diagnostic.json'
            if comparison_path.exists():
                comparison=read(comparison_path);forty=comparison['native40'];record['native40']=forty
                d=read(diagnostic_path);require(comparison==dict(d,native40=forty),'GL diagnostic changed during comparison')
                if forty['available']:
                    valid=finite(forty['loglik']) and len(forty['gradient'])==29 and all(finite(x) for x in forty['gradient'])
                    record['native40_finite']=valid
                    require(forty['order_refinement_applicable']==(fit['v']>0),'Native comparison applicability differs')
                    if valid:
                        dl=forty['loglik']-d['loglik'];dg=max(abs(a-b) for a,b in zip(forty['gradient'],d['gradient']))
                        require(dl==forty['LL_difference'] and dg==forty['score_max_difference'],
                                'Native comparison arithmetic differs')
                        agree=abs(dl)<=1e-8+1e-8*abs(forty['loglik']) and all(
                            abs(a-b)<=1e-8+1e-6*abs(b) for a,b in zip(forty['gradient'],d['gradient']))
                        require(agree==forty['LL_score_agree'],'Native comparison criterion differs')
                    else:require(forty['LL_score_agree'] is False,'Nonfinite native comparison claimed agreement')
            local=[]
            if not failed and not record['score_checks_passed']:local.append('successful_recipe_GL_score_check_not_passed')
            if record.get('GL_allowance_check_passed') is False:local.append('saved_GL_allowance_check_not_passed')
            if summary['diagnostic_failure'] is not None:local.append('post_recipe_diagnostic_failure')
            if local:issues.append(dict(draw=draw,issues=local))
        collector_sha=hashlib.sha256(body(Path(__file__).resolve())).hexdigest()
        report=dict(completed_pilots=sum(r['completed_record'] for r in records),required_pilots=2,
            collector_sha256=collector_sha,
            recipe_failed_pilots=sum(r.get('prescribed_recipe_failed_under_GL',False) for r in records),
            scheduler_states=dict(Counter(states.get(j['jobid'],['MISSING'])[0] for j in jobs)),
            records=records,issues=issues,input_sha256=digests,bytes_read=bytes_read,
            synthetic_fixture_participants=60,variance_coordinate='v=tau**2',
            no_fifty_draw_agreement=True,no_native_draw_replacement=True,no_calibration=True,no_human_validation=True)
        with (args.out/'collection.json').open('x') as f:json.dump(report,f,allow_nan=False,indent=2);f.write('\n')
        print(json.dumps(dict(completed=report['completed_pilots'],required=2,
            prescribed_recipe_failures_under_GL=report['recipe_failed_pilots'],issues=len(issues))))
    except Exception:
        (args.out/'failure.txt').write_text(traceback.format_exc());raise


if __name__=='__main__':main()
