#!/usr/bin/env python3
"""Run remaining candidate branches after their actual inputs pass checks."""
import argparse
import datetime as dt
import json
from pathlib import Path
import subprocess
import time
from final_annotation_paths import annotation_paths


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--account',default='nslab')
    p.add_argument('--qos',default='nslab')
    p.add_argument('--job-name',default='test')
    args=p.parse_args();c=args.candidate.resolve();a=Path(__file__).resolve().parents[1]
    path=c/'downstream_jobs.json'
    jobs=json.loads(path.read_text()) if path.exists() else {}
    def ready(stage):
        q=c/stage/'completion.json'
        return q.exists() and json.loads(q.read_text()).get('status')=='completed'
    def submit(name,command):
        if name in jobs:return
        assert not (c/name).exists(),f'Existing incomplete {name}; inspect before retry'
        cmd=['sbatch','--parsable',f'--account={args.account}',f'--qos={args.qos}',
             '--cpus-per-task=4','--mem=32G',f'--job-name={args.job_name}',str(a/'run_final_analysis.sbatch')]+command
        jid=subprocess.check_output(cmd,text=True).strip().split(';')[0]
        jobs[name]=dict(job_id=jid,command=cmd,submitted_utc=dt.datetime.now(dt.timezone.utc).isoformat())
        path.write_text(json.dumps(jobs,indent=2));print('Submitted',name,jid,flush=True)
    py=lambda name:[str(a/'scripts'/name),'--candidate',str(c)]
    rscript='/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript'
    while True:
        active=set()
        ids=[v['job_id'] for v in jobs.values()]
        if ids:active=set(subprocess.check_output(['squeue','-h','-j',','.join(ids),'-o','%i'],text=True).split())
        for name,job in jobs.items():
            marker=c/name/'completion.json' if name!='finalize' else c/'numerical_completion.json'
            if marker.exists() and json.loads(marker.read_text()).get('status')=='completed':continue
            age=(dt.datetime.now(dt.timezone.utc)-dt.datetime.fromisoformat(job['submitted_utc'])).total_seconds()
            if job['job_id'] not in active and age>120:raise RuntimeError('Inspect failed downstream job '+name+' '+job['job_id'])
        rna_plan=c/'rna_jobs.json'
        if rna_plan.exists() and json.loads(rna_plan.read_text()).get('status')=='completed':
            command=py('final_rna_pseudobulk.py')+['--output',str(c/'rna_pseudobulk')]
            key=c/'verified_cross_assay_pairing.tsv'
            if key.exists():command+=['--verified-pairing',str(key)]
            submit('rna_pseudobulk',command)
        if ready('rna_pseudobulk'):
            submit('rna_models',['--command',rscript,str(a/'scripts/final_rna_models.R'),str(a),str(c/'rna_pseudobulk'),str(c/'rna_models')])
        if ready('rna_models'):
            submit('donor_associations',['--command',rscript,str(a/'scripts/final_donor_associations.R'),str(a),str(c),str(c/'donor_associations')])
            submit('disease_context',py('final_candidate_summary.py')+['--output',str(c/'disease_context')])
        pred_plan=c/'prediction_jobs.json'
        if pred_plan.exists() and json.loads(pred_plan.read_text()).get('status')=='completed':
            submit('prediction_summary',py('final_prediction_summary.py')+['--output',str(c/'prediction_summary')])
        if all(ready(stage) for stage in ['rna_models','disease_context','prediction_summary']):
            submit('panels',py('render_final_candidate.py')+['--context',str(c/'disease_context'),'--output',str(c/'panels')])
        motifs=all(ready(str(annotation_paths(c)['motif_context']/f'{cohort}_{lineage}')) for cohort in ['GSE244832','GSE281367'] for lineage in ['hepatocyte','stellate'])
        if motifs and all(ready(stage) for stage in ['rna_models','donor_associations','disease_context','prediction_summary','panels']):
            if 'finalize' not in jobs:
                cmd=['sbatch','--parsable',f'--account={args.account}',f'--qos={args.qos}',
                     '--cpus-per-task=4','--mem=32G',f'--job-name={args.job_name}',str(a/'run_final_analysis.sbatch')]+py('finalize_final_candidate.py')
                jid=subprocess.check_output(cmd,text=True).strip().split(';')[0]
                jobs['finalize']=dict(job_id=jid,command=cmd,submitted_utc=dt.datetime.now(dt.timezone.utc).isoformat())
                path.write_text(json.dumps(jobs,indent=2));print('Submitted final conclusion',jid,flush=True)
        if (c/'numerical_completion.json').exists():
            submit('independent_numerical_review',py('verify_final_candidate.py'))
            if ready('independent_numerical_review'):
                print('Numerical candidate and independent reconstruction complete; final PDF and authority review remains required',flush=True)
                break
        print(dt.datetime.now(dt.timezone.utc).isoformat(),'waiting for qualified branch inputs',flush=True)
        time.sleep(45)


if __name__=='__main__':main()
