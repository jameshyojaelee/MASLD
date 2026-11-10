#!/usr/bin/env python3
"""Submit individual authorized refits, at most six outstanding at a time."""
import argparse
import datetime as dt
import json
from pathlib import Path
import subprocess
import time


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--once',action='store_true')
    p.add_argument('--memory-gib',type=int,default=32)
    p.add_argument('--walltime',default='06:00:00')
    p.add_argument('--account',default='nslab')
    p.add_argument('--qos',default='nslab')
    p.add_argument('--job-name',default='test')
    args=p.parse_args()
    a=Path(__file__).resolve().parents[1]
    candidate=args.candidate.resolve()
    marker=candidate/'prediction_baseline_pilot/baseline_reproduction.json'
    assert json.loads(marker.read_text())['status']=='passed','Baseline pilot must pass first'
    spec=json.loads((candidate/'prediction_inputs_v2/specification.json').read_text())
    path=candidate/'prediction_jobs.json'
    if path.exists():
        plan=json.loads(path.read_text())
    else:
        jobs=[]
        for scenario in spec['scenarios']:
            name=scenario['name']
            if name=='baseline': continue
            for fold in range(5):
                reuse=bool(scenario.get('omit_donor') and spec['outer_map'][scenario['omit_donor']]==fold)
                jobs.append(dict(scenario=name,fold=fold,reuse_original=reuse,
                    output=str(candidate/f'prediction_refits/{name}/fold{fold}'),
                    reason='omitted_donor_is_held_out_so_training_and_evaluation_unchanged' if reuse else 'refit_required'))
        plan=dict(jobs=jobs,maximum_outstanding_jobs=6,cpus_per_job=8,memory_GiB_per_job=64,
                  walltime='90:00:00',array=False,scientific_adoption=False)
        path.write_text(json.dumps(plan,indent=2))
    previous=dict(memory_GiB_per_job=plan['memory_GiB_per_job'],walltime=plan['walltime'])
    if previous != dict(memory_GiB_per_job=args.memory_gib,walltime=args.walltime):
        plan.setdefault('resource_adjustments',[]).append(dict(previous=previous,
            memory_GiB_per_job=args.memory_gib,walltime=args.walltime,
            reason='Seven completed refits peak below 19 GiB and finish within 65 minutes',
            changed_utc=dt.datetime.now(dt.timezone.utc).isoformat()))
    plan['memory_GiB_per_job']=args.memory_gib
    plan['walltime']=args.walltime
    while True:
        ids=[j['job_id'] for j in plan['jobs'] if j.get('job_id')]
        alive=set()
        if ids:
            text=subprocess.check_output(['squeue','-h','-j',','.join(ids),'-o','%i'],text=True)
            alive=set(text.split())
        completed=0; failures=[]
        for job in plan['jobs']:
            if job['reuse_original']:
                completed+=1; continue
            marker=Path(job['output'])/'completion.json'
            if marker.exists():
                record=json.loads(marker.read_text())
                if record.get('status')=='completed' and not record.get('technical_pilot',True):
                    completed+=1; job['status']='completed'; continue
            if job.get('job_id') and job['job_id'] not in alive:
                # Accounting can lag briefly after submission, but never resubmit
                # automatically into an existing output or hide a failed fit.
                submitted=dt.datetime.fromisoformat(job['submitted_utc'])
                if (dt.datetime.now(dt.timezone.utc)-submitted).total_seconds()>120:
                    failures.append(job); job['status']='failed_or_missing_completion'
        if failures:
            plan['status']='requires_failure_review'; plan['failures']=failures
            path.write_text(json.dumps(plan,indent=2)); raise RuntimeError('Inspect failed refits before continuing')
        slots=6-len(alive)
        for job in plan['jobs']:
            if slots<=0: break
            if job['reuse_original'] or job.get('job_id'): continue
            output=Path(job['output']); output.parent.mkdir(parents=True,exist_ok=True)
            assert not output.exists(),'Never overwrite a refit'
            cmd=['sbatch','--parsable',f'--account={args.account}',f'--qos={args.qos}',
                 f'--mem={args.memory_gib}G',f'--time={args.walltime}',
                 f'--job-name={args.job_name}',
                 str(a/'run_final_analysis.sbatch'),str(a/'scripts/final_prediction_refits.py'),'fold',
                 '--root',str(candidate/'prediction_inputs_v2'),'--scenario',job['scenario'],
                 '--fold',str(job['fold']),'--output',str(output),'--workers','8']
            jid=subprocess.check_output(cmd,text=True).strip().split(';')[0]
            job.update(job_id=jid,command=cmd,status='submitted',submitted_utc=dt.datetime.now(dt.timezone.utc).isoformat())
            print('Submitted',jid,job['scenario'],job['fold'],flush=True)
            slots-=1
            path.write_text(json.dumps(plan,indent=2))
        plan['status']='completed' if completed==len(plan['jobs']) else 'running'
        plan['completed_or_exact_reuse']=completed
        plan['updated_utc']=dt.datetime.now(dt.timezone.utc).isoformat()
        path.write_text(json.dumps(plan,indent=2))
        print(plan['updated_utc'],plan['status'],completed,'/',len(plan['jobs']),flush=True)
        if plan['status']=='completed' or args.once: break
        time.sleep(45)


if __name__=='__main__': main()
