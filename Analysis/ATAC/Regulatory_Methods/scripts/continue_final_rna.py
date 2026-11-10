#!/usr/bin/env python3
"""Continue the authorized joint RNA reconstruction after pilot correctness checks."""
import argparse
import datetime as dt
import json
from pathlib import Path
import subprocess
import time

import pandas as pd
from check_final_rna_recount import check


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--account',default='nslab')
    p.add_argument('--qos',default='nslab')
    p.add_argument('--job-name',default='test')
    p.add_argument('--maximum-outstanding',type=int,default=4)
    args=p.parse_args(); c=args.candidate.resolve(); a=Path(__file__).resolve().parents[1]
    pilot=c/'rna_recovery/JB_303_attempt2'
    while not (pilot/'completion.json').exists():
        active=subprocess.check_output(['squeue','-h','-j','22407389','-o','%i'],text=True).strip()
        if not active: raise RuntimeError('RNA pilot failed: inspect job 22407389 before scaling')
        time.sleep(45)
    check(pilot)
    source=json.loads((c/'rna_source_audit/source_workbook.json').read_text())
    libraries=sorted({r[1] for rows in source.values() for r in rows[1:] if r[1]})
    path=c/'rna_jobs.json'
    if path.exists(): plan=json.loads(path.read_text())
    else:
        plan=dict(jobs=[dict(library=lib,output=str(c/f'rna_recovery/{lib}')) for lib in libraries if lib!='JB_303'],
                  maximum_outstanding=4,cpus_per_job=8,memory_GiB_per_job=64,walltime='90:00:00',array=False,
                  pilot=dict(library='JB_303',job_id='22407389',output=str(pilot),scientific_checks='passed'))
        path.write_text(json.dumps(plan,indent=2))
    assert args.maximum_outstanding > 0
    plan['maximum_outstanding']=args.maximum_outstanding
    while True:
        ids=[j['job_id'] for j in plan['jobs'] if j.get('job_id')]
        alive=set(subprocess.check_output(['squeue','-h','-j',','.join(ids),'-o','%i'],text=True).split()) if ids else set()
        completed=0
        for job in plan['jobs']:
            directory=Path(job['output'])
            if (directory/'completion.json').exists():
                if not (directory/'scientific_checks.json').exists(): check(directory)
                assert json.loads((directory/'scientific_checks.json').read_text())['status']=='passed'
                job['status']='completed'; completed+=1;continue
            if job.get('job_id') and job['job_id'] not in alive:
                submitted=dt.datetime.fromisoformat(job['submitted_utc'])
                if (dt.datetime.now(dt.timezone.utc)-submitted).total_seconds()>120:
                    job['status']='failed_requires_review'; path.write_text(json.dumps(plan,indent=2))
                    raise RuntimeError(f"RNA library {job['library']} failed; do not sum its run matrices")
        slots=args.maximum_outstanding-len(alive)
        for job in plan['jobs']:
            if slots<=0: break
            if job.get('job_id'):continue
            assert not Path(job['output']).exists()
            cmd=['sbatch','--parsable',f'--account={args.account}',f'--qos={args.qos}',
                 f'--job-name={args.job_name}',str(a/'run_final_analysis.sbatch'),
                 str(a/'scripts/final_rna_recovery.py'),'--candidate',str(c),'--library',job['library'],'--cpus','8','--memory','56']
            jid=subprocess.check_output(cmd,text=True).strip().split(';')[0]
            job.update(job_id=jid,command=cmd,status='submitted',submitted_utc=dt.datetime.now(dt.timezone.utc).isoformat())
            path.write_text(json.dumps(plan,indent=2)); slots-=1
            print('Submitted RNA',jid,job['library'],flush=True)
        plan['completed_libraries']=completed+1
        plan['status']='completed' if completed==len(plan['jobs']) else 'running'
        path.write_text(json.dumps(plan,indent=2))
        print(dt.datetime.now(dt.timezone.utc).isoformat(),plan['status'],completed+1,'/18',flush=True)
        if plan['status']=='completed':break
        time.sleep(45)


if __name__=='__main__':main()
