#!/usr/bin/env python3
"""Collect all fifteen fixed Gamma GL mass checks; no fit or calibration claim."""
import argparse
import ast
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import struct
import subprocess
import sys

ROOT=Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE=ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
REC=ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
LAUNCH=REC/'repair_p3_gamma_GL_five_gene_v3_launch_20261007T051812Z'
REFERENCE=REC/'repair_p3_five_gene_accuracy_launch_20261005T143525Z/five_gene_accuracy_collection.json'
RUNNER_SHA='ac62db5d5a6bae283ddd73bad09f6b6fa72ce5362282f6e6082bdcba26a6803b'
PINS_SHA='e02e4b3ec5d37cf474c3f640749015c84d3bef107f6fd831d5730e3563235245'
SOURCES=('likelihood_v1.py','repair_p3_bb_gamma_difference_v3.py','check_quadrature_v1.py',
    'repair_p3_integral.py','repair_p3_composite_gl.py')


def require(ok,message):
    if not ok:raise ValueError(message)


def finite(x):return type(x) in (int,float) and math.isfinite(x)


def main():
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--out',type=Path,required=True);args=ap.parse_args()
    require(args.out.resolve().parent==REC.resolve(),'Output must be a new direct child of REC')
    args.out.mkdir(exist_ok=False);digests={};bytes_read=0
    def body(p):
        nonlocal bytes_read
        require(p.resolve().is_relative_to(ROOT.resolve()) and not p.is_symlink(),'Input escaped repository')
        data=p.read_bytes();require(len(data)<=32*1024**2 and bytes_read+len(data)<=256*1024**2,'Read cap exceeded')
        bytes_read+=len(data);digests[str(p)]=hashlib.sha256(data).hexdigest();return data
    def read(p):
        def invalid(x):raise ValueError('Nonfinite JSON constant: '+x)
        return json.loads(body(p),parse_constant=invalid)
    jobs=read(LAUNCH/'jobs.json');pins=read(LAUNCH/'pins.json')
    require(digests[str(LAUNCH/'pins.json')]==PINS_SHA,'Pinned input family differs')
    require([j['jobid'] for j in jobs]==[str(n) for n in range(22082310,22082325)]
        and [(j['gene'],j['tau_multiple']) for j in jobs]==[(g,m) for m in (.25,1.,4.) for g in range(5)],
        'Fixed fifteen-point roster differs')
    require(hashlib.sha256(body(LAUNCH/'repair_p3_gamma_GL_five_gene_check.py')).hexdigest()==RUNNER_SHA,'Runner differs')
    for p,h in pins.items():require(hashlib.sha256(body(Path(p))).hexdigest()==h,'Pinned source/reference differs')
    settings_class=next(n for n in ast.parse(body(HERE/'repair_p3_integral.py')).body
        if isinstance(n,ast.ClassDef) and n.name=='Settings')
    settings={n.target.id:ast.literal_eval(n.value) for n in settings_class.body if isinstance(n,ast.AnnAssign)}
    reference=read(REFERENCE);require(reference['execution_complete'],'Independent scalar reference incomplete')
    refs={(x['gene'],x['tau_multiple']):x for x in reference['results']}
    require(len(refs)==15,'Independent scalar family differs')
    fixture=read(REC/('repair_p3_five_gene_accuracy_'+reference['results'][0]['job_id'])/'fixture.json')
    counts=fixture['drawn_counts'];require(len(counts)==805 and all(type(x) is int for x in counts),'Reference counts differ')
    count_sha=hashlib.sha256(struct.pack('<805q',*counts)).hexdigest()
    status=subprocess.check_output(['sacct','-j',','.join(j['jobid'] for j in jobs),'--noheader','--parsable2',
        '--format=JobIDRaw,State,ExitCode,Elapsed,MaxRSS'],text=True);(args.out/'sacct.txt').write_text(status)
    states={p[0]:p[1:] for line in status.splitlines() if (p:=line.split('|'))[0].isdigit()}
    records=[];issues=[]
    for job in jobs:
        gene=job['gene'];multiple=job['tau_multiple'];d=Path(job['output_directory']);state=states.get(job['jobid'],['MISSING'])
        require(d==REC/('repair_p3_gamma_GL_five_gene_check_'+job['jobid']),'Output identity differs')
        record=dict(gene=gene,tau_multiple=multiple,jobid=job['jobid'],state=state,completed=False);records.append(record)
        if state[:2]!=['COMPLETED','0:0']:
            if state[0] not in ('RUNNING','PENDING','CONFIGURING','COMPLETING','SUSPENDED'):
                issues.append(dict(jobid=job['jobid'],reason='scheduler_not_successful',state=state))
            continue
        rp=d/f'gene{gene}_tau{.05*multiple}.json'
        missing=[str(p) for p in (rp,d/'protocol.json',d/'summary.json',d/'native_GH.json') if not p.is_file()]
        if missing:issues.append(dict(jobid=job['jobid'],reason='terminal_files_missing',files=missing));continue
        protocol=read(d/'protocol.json');summary=read(d/'summary.json');answer=read(rp);native=read(d/'native_GH.json')
        require(protocol['source_and_input_sha256']==pins and protocol['seed']==20260923 and protocol['required_points']==15
            and protocol['requested_gene']==gene and protocol['requested_tau_multiple']==multiple
            and protocol['full_drawn_count_vector_sha256']==count_sha
            and protocol['executed_GL_settings']==settings,'Executed fixture/point/settings differ')
        require(protocol['rows_per_gene']==[5,25,75,200,500] and protocol['dense_G5_GL_not_G8_batched_scope_admission'],
                'Executed synthetic scope differs')
        for name in SOURCES:
            require(hashlib.sha256(body(d/'executed_sources'/name)).hexdigest()==pins[str(HERE/name)],'Archived helper differs')
        require(hashlib.sha256(body(d/'executed_sources/repair_p3_gamma_GL_five_gene_check.py')).hexdigest()==RUNNER_SHA,
                'Archived runner differs')
        require(summary['completed'] and summary['points']==1 and summary['required_family_points']==15
            and summary['requested_gene']==gene and summary['requested_tau_multiple']==multiple,'Summary scope differs')
        ref=refs[gene,multiple]['record'];require(answer['gene']==gene and answer['tau']==.05*multiple
            and answer['variance']==(.05*multiple)**2 and answer['scalar']==ref['scalar'],'Scalar/point differs')
        require(native=={k:v for k,v in answer.items() if k!='GL'},'Native output changed during GL')
        require([o['order'] for o in answer['orders']]==[20,40],'Native orders differ')
        scalar=ref['scalar']['expanded']['logmass'];require(finite(scalar),'Scalar mass nonfinite')
        record['scalar_relative_LL_absolute_tolerance']=1e-8*abs(scalar)
        record['native_scalar_differences']={}
        native_ok=True
        for order in answer['orders']:
            if not order['returned']:native_ok=False;continue
            dl=order['logmass']-scalar;ok=abs(dl)<=1e-8*abs(scalar)
            require(finite(order['logmass']) and len(order['gradient'])==17 and all(finite(x) for x in order['gradient'])
                and dl==order['scalar_difference'] and ok==order['scalar_relative_LL_pass'],
                    'Native scalar comparison differs')
            native_ok=native_ok and ok
            record['native_scalar_differences'][str(order['order'])]=dl
        if all(o['returned'] for o in answer['orders']):
            dl=answer['orders'][0]['logmass']-answer['orders'][1]['logmass'];ok=abs(dl)<=1e-8
            require(dl==answer['GH20_GH40_difference'] and ok==answer['GH20_GH40_absolute_pass'],'Native refinement differs')
            native_ok=native_ok and ok
        gl=answer['GL'];require(gl['returned']==summary['GL_returned'] and type(gl['Gamma_BB_calls']) is int
            and gl['Gamma_BB_calls']>0,'GL execution/return evidence differs')
        gl_ok=False
        if gl['returned']:
            g=gl['gradient'];h=gl['hessian'];require(finite(gl['logmass']) and len(g)==17 and len(h)==17
                and all(finite(x) for x in g) and all(len(row)==17 and all(finite(x) for x in row) for row in h),
                'GL arrays invalid')
            raw=read(d/'GL_returned_raw.json');require(raw=={k:gl[k] for k in ('logmass','gradient','hessian','report')},
                'Raw returned GL arrays differ')
            report=gl['report'];refined=report['refined']
            require(refined['logmass']==gl['logmass'] and refined['gradient']==g and refined['hessian']==h,
                    'Returned GL arrays differ from final integration report')
            evaluations=report['likelihood_evaluations'];rows=[5,25,75,200,500][gene]
            require(type(evaluations) is int and evaluations>0 and
                2*rows*evaluations<=gl['Gamma_BB_calls']<=3*rows*evaluations,
                'Gamma BB counter outside observed-plus-selection call bounds')
            require(all(len(x)==17 for x in (report['gradient_refinement_difference'],refined['gradient_error_allowance']))
                and all(len(x)==17 and all(len(row)==17 for row in x) for x in
                    (report['hessian_refinement_difference'],refined['hessian_error_allowance'])),
                'Empirical refinement array dimensions differ')
            components=[report['logmass_refinement_difference'],refined['logmass_error_allowance'],
                *report['gradient_refinement_difference'],*refined['gradient_error_allowance'],
                *(x for row in report['hessian_refinement_difference'] for x in row),
                *(x for row in refined['hessian_error_allowance'] for x in row)]
            require(all(finite(x) and x>=0 for x in components),'Invalid empirical refinement allowances')
            errors=dict(logmass=report['logmass_refinement_difference']+refined['logmass_error_allowance'],
                gradient=max(a+b for a,b in zip(report['gradient_refinement_difference'],refined['gradient_error_allowance'])),
                hessian=max(a+b for ra,rb in zip(report['hessian_refinement_difference'],refined['hessian_error_allowance'])
                    for a,b in zip(ra,rb)))
            require(errors['logmass']<=settings['logmass_tolerance'] and errors['gradient']<=settings['gradient_tolerance']
                and errors['hessian']<=settings['hessian_tolerance'],'Empirical GL refinement allowance exceeds unchanged settings')
            dl=gl['logmass']-scalar;gl_ok=abs(dl)<=1e-8*abs(scalar)
            require(dl==gl['scalar_difference'] and gl_ok==gl['scalar_relative_LL_pass'],'GL scalar comparison differs')
            inactive=[j for j in range(10) if j not in (gene,5+gene)]
            require(all(g[j]==0 and all(x==0 for x in h[j]) and all(row[j]==0 for row in h) for j in inactive),
                    'Cross-gene derivative layout differs')
            record.update(GL_scalar_difference=dl,GL_scalar_relative_LL_pass=gl_ok,
                finite17scores_and17x17Hessian=True,cross_gene_derivatives_exact_zero=True,
                Gamma_BB_counter_within_per_node_support_call_bounds=True,
                empirical_refinement_allowance_maxima=errors,
                Hessian_max_asymmetry=max(abs(h[i][j]-h[j][i]) for i in range(17) for j in range(17)))
        else:
            require(gl['known_numerical_refusal'] and summary['GL_known_numerical_refusal'],
                    'GL refusal not classified as numerical')
            record['GL_refusal']=gl
            record['GL_scalar_relative_LL_pass']=False
        passed=native_ok and gl_ok
        require(summary['selected_point_likelihood_criteria_passed']==passed
            and summary['GL_only_scalar_relative_LL_pass']==gl.get('scalar_relative_LL_pass')
            and summary['returned_native_integrals']==sum(o['returned'] for o in answer['orders']),
            'Selected-point criterion/count differs')
        record.update(completed=True,GL_returned=gl['returned'],Gamma_BB_calls=gl['Gamma_BB_calls'],
            native_criteria_passed=native_ok,point_likelihood_criteria_passed=passed)
    complete=all(x['completed'] for x in records)
    failed_jobs=[r['jobid'] for r in records if r['state'][0].split()[0].rstrip('+') in
        ('FAILED','TIMEOUT','OUT_OF_MEMORY','CANCELLED','NODE_FAIL','PREEMPTED','BOOT_FAIL','DEADLINE','REVOKED')]
    collector_sha=hashlib.sha256(body(Path(__file__).resolve())).hexdigest()
    report=dict(complete_family=complete,completed_points=sum(x['completed'] for x in records),required_points=15,
        all15_likelihood_criteria_passed=(all(x['point_likelihood_criteria_passed'] for x in records) if complete else None),
        all15_GL_only_scalar_relative_LL_passed=(all(x['GL_scalar_relative_LL_pass'] for x in records) if complete else None),
        terminal_failure_jobids=failed_jobs,family_assessment_status=('incomplete_due_to_job_or_record_error' if issues
            else 'complete' if complete else 'in_progress'),
        records=records,issues=issues,scheduler_states=dict(Counter(states.get(j['jobid'],['MISSING'])[0] for j in jobs)),
        collector_sha256=collector_sha,input_sha256=digests,bytes_read=bytes_read,
        scalar_and_quadrature_error_estimates_not_certified=True,no_independent_gradient_or_Hessian_accuracy=True,
        no_G8_batched_scope_generalization=True,no_fits=True,no_calibration=True,no_human_validation=True)
    with (args.out/'collection.json').open('x') as f:json.dump(report,f,allow_nan=False,indent=2);f.write('\n')
    print(json.dumps(dict(complete=complete,completed=report['completed_points'],required=15,issues=len(issues))))
    if issues:sys.exit(1)


if __name__=='__main__':main()
