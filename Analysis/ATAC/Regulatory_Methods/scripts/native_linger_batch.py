#!/usr/bin/env python3
"""Frozen-source, four-worker maximum LINGER DPM1 donor-training experiment.

Runs only native training/link extraction. Completed links still require the
common heldout-donor readout; prior independence remains unverified.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''): h.update(block)
    return h.hexdigest()


def write_json(path,value):
    with Path(path).open('x') as f: json.dump(value,f,indent=2)


def check_units(stage,fit,feature_genes):
    import numpy as np
    import pandas as pd
    info=json.loads((stage/'stage.json').read_text())
    meta=pd.read_csv(stage/'metacells.tsv',sep='\t',dtype=str)
    if set(meta.donor_id)&set(info['excluded_donors']): raise ValueError('Excluded donor in training')
    if set(meta.groupby('donor_id').size())!={2}: raise ValueError('Unequal donor input budget')
    rna=pd.read_csv(stage/'rna.tsv',sep='\t'); atac=pd.read_csv(stage/'atac.tsv',sep='\t')
    actual=pd.read_csv(fit/'native/Exp.txt',sep='\t',header=None).to_numpy()
    if not np.allclose(actual.ravel(),np.log2(1+rna[info['gene_id']]),rtol=1e-12,atol=1e-12):
        raise ValueError('Native response transform mismatch')
    openness=pd.read_csv(fit/'native/Openness.txt',sep='\t',header=None).to_numpy()
    if not np.allclose(openness,atac.to_numpy().T,rtol=1e-12,atol=1e-12):
        raise ValueError('Native ATAC units/order mismatch')
    tfs=pd.read_csv(fit/'native/TFName.txt',sep='\t',header=None)[0]
    tfexp=pd.read_csv(fit/'native/TFexp.txt',sep='\t',header=None).to_numpy()
    genes=pd.read_csv(feature_genes,sep='\t'); mapping=dict(zip(genes.gene_id,genes.gene_symbol))
    positions={mapping[g]:i for i,g in enumerate(rna.columns)}
    expected=np.log2(1+rna.iloc[:,[positions[tf] for tf in tfs]].to_numpy().T)
    if not np.allclose(tfexp,expected,rtol=1e-12,atol=1e-12):
        raise ValueError('Native full-TF transform/order mismatch')
    seed_audit=json.loads((fit/'seed_and_background_audit.json').read_text())
    if seed_audit['native_SEED']!=42 or seed_audit['native_GLOBAL_SEED']!=42 or seed_audit['python_hash_seed_at_process_start']!='42':
        raise ValueError('Unexpected native seed defaults or startup hash seed')
    if any(c['seed']!=42 or c['torch_initial_seed']!=42 for c in seed_audit['seed_calls']):
        raise ValueError('Native seed override differs from42')
    backgrounds=seed_audit['backgrounds']
    if len(backgrounds)!=1 or len(set(backgrounds[0]['row_indices']))!=min(50,len(meta)):
        raise ValueError('SHAP background sampling mismatch')
    if set(backgrounds[0]['donor_ids'])&set(info['excluded_donors']): raise ValueError('Heldout background donor')
    prior=seed_audit['prior_alignments']
    if len(prior)!=1 or not prior[0]['features'] or not all(r['match'] and r['expected_identity']==r['actual_identity'] for r in prior[0]['features']):
        raise ValueError('Prior full-feature identity audit failed')
    if any((r['feature_type']=='RE')!=(r['prior_column0']>=prior[0]['bulk_TFs']) for r in prior[0]['features']):
        raise ValueError('TF/RE prior block boundary mismatch')
    if info['control']=='pairing_shuffle':
        paired=pd.read_csv(stage.parent/'paired/rna.tsv',sep='\t')
        for _,idx in meta.groupby('donor_id').indices.items():
            digest=lambda x:sorted(hashlib.sha256(np.ascontiguousarray(v).tobytes()).hexdigest() for v in x)
            if digest(rna.iloc[idx].to_numpy())!=digest(paired.iloc[idx].to_numpy()):
                raise ValueError('RNA shuffle did not preserve donor RNA rows')
        if sha(stage/'atac.tsv')!=sha(stage.parent/'paired/atac.tsv'):
            raise ValueError('RNA shuffle changed ATAC')
    links=pd.read_csv(fit/'links.tsv',sep='\t')
    if links.duplicated(['gene_id','peak_id']).any() or not np.isfinite(links.score).all():
        raise ValueError('Invalid cis link scores')
    if not set(links.gene_id).issubset({info['gene_id']}) or not set(links.peak_id).issubset(atac.columns):
        raise ValueError('Result feature identities mismatch')
    return {'status':'PASS','native_TFs':len(tfs),'training_donors':meta.donor_id.nunique(),
            'training_aggregates':len(meta),'native_links':len(links),'donor_ids':sorted(set(meta.donor_id))}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['paired-root','annotation','prior-root','model-preflight','python','prepare-python','output']:
        p.add_argument('--'+name,type=Path,required=True)
    p.add_argument('--workers',type=int,default=4); p.add_argument('--timeout',type=int,default=600)
    a=p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('Existing compute allocation required')
    if not 1<=a.workers<=4 or not 1<=a.timeout<=600: raise ValueError('Worker/time cap exceeded')
    for name,value in vars(a).items():
        if isinstance(value,Path): setattr(a,name,value.resolve())
    a.output.mkdir(parents=True,exist_ok=False)
    snapshot=a.output/'source_snapshot'; snapshot.mkdir()
    sources={}
    for name in ['native_paired_methods.py','native_paired_methods.R','network_compare.py','native_linger_batch.py']:
        original=Path(__file__).parent/name
        before=sha(original); shutil.copy2(original,snapshot/name)
        if before!=sha(snapshot/name) or before!=sha(original): raise RuntimeError('Source changed during snapshot')
        sources[name]=before
    write_json(a.output/'source_sha256.json',sources)
    for name in ['paired-root','annotation','prior-root','model-preflight','python','prepare-python']:
        if not getattr(a,name.replace('-','_')).exists(): raise FileNotFoundError(name)
    manifest={str(path):sha(path) for path in [a.annotation,a.model_preflight,
      a.paired_root/'cells.tsv',a.paired_root/'genes.tsv',a.paired_root/'peaks.tsv',
      a.paired_root/'ATAC_counts.mtx.gz',a.paired_root/'rna_counts.npz',
      a.paired_root/'paired_export.json',a.paired_root/'RNA_export.json']}
    write_json(a.output/'input_sha256.json',manifest)
    package_root=a.python.parent.parent/'lib/python3.10/site-packages/LingerGRN'
    package_sources={name:sha(package_root/name) for name in ['preprocess.py','pseudo_bulk.py','LINGER_tr.py','LL_net.py']}
    write_json(a.output/'native_package_source_sha256.json',package_sources)
    env=os.environ.copy()
    env.update(PYTHONNOUSERSITE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1',CUDA_VISIBLE_DEVICES='',SLURM_CPUS_PER_TASK='1',PYTHONUNBUFFERED='1',PYTHONHASHSEED='42')
    env['PATH']=str(a.python.parent)+':'+env.get('PATH','')
    env['LD_LIBRARY_PATH']=str(a.python.parent.parent/'lib')
    driver=snapshot/'native_paired_methods.py'; inputs=a.output/'inputs'
    prepare_env=env.copy()
    prepare_env['PATH']=str(a.prepare_python.parent)+':'+os.environ.get('PATH','')
    prepare_env['LD_LIBRARY_PATH']=str(a.prepare_python.parent.parent/'lib')
    with (a.output/'prepare.log').open('x') as log:
        subprocess.run([str(a.prepare_python),str(driver),'prepare-linger-stage','--paired-root',str(a.paired_root),
          '--genes',str(a.annotation),'--target','ENSG00000000419','--stage','all','--output',str(inputs)],
          env=prepare_env,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=600)
    stage_names=[f'outer{i}'+suffix for i in range(5) for suffix in ['', '_inner0','_inner1','_inner2']]
    jobs=[(name,control) for name in stage_names for control in ['paired','pairing_shuffle']]
    def worker(job):
        name,control=job; stage=inputs/name/control; output=a.output/'fits'/name/control
        output.parent.mkdir(parents=True,exist_ok=True)
        logpath=a.output/'fits'/name/(control+'.log')
        started=time.time()
        try:
            input_hashes={f:sha(stage/f) for f in ['rna.tsv','atac.tsv','metacells.tsv','stage.json']}
            with logpath.open('x') as log:
                subprocess.run([str(a.python),str(driver),'fit-linger-stage','--stage',str(stage),
                  '--genes',str(a.paired_root/'genes.tsv'),'--peaks',str(a.paired_root/'peaks.tsv'),
                  '--prior-root',str(a.prior_root),'--model-preflight',str(a.model_preflight),
                  '--chromosome','chr20','--shap-background','training_cap','--output',str(output)],env=env,stdout=log,stderr=subprocess.STDOUT,
                  check=True,timeout=a.timeout)
            checks=check_units(stage,output,a.paired_root/'genes.tsv')
            if any(sha(stage/f)!=h for f,h in input_hashes.items()): raise RuntimeError('Stage inputs changed')
            write_json(output/'native_units_check.json',checks)
            write_json(output/'source_and_input_sha256.json',{'source':sources,'inputs':input_hashes})
            combined=a.output/'stages'/name/control; combined.mkdir(parents=True)
            shutil.copy2(output/'links.tsv',combined/'links.tsv')
            write_json(combined/'links.training_donors.json',{'donor_ids':checks['donor_ids']})
            write_json(combined/'completion.json',{'status':'completed','scope':'native_training_links_only',
              'input_units':'LINGER_per_nucleus_transform_then_mean','native_units_check':'PASS'})
            return {'stage':name,'control':control,'status':'completed','seconds':time.time()-started,
                    'links':str((combined/'links.tsv').resolve()),'n_links':checks['native_links']}
        except Exception as e:
            return {'stage':name,'control':control,'status':'failed','seconds':time.time()-started,'error':repr(e)}
    results=[]
    with ThreadPoolExecutor(max_workers=a.workers) as pool:
        for future in as_completed([pool.submit(worker,job) for job in jobs]):
            result=future.result(); results.append(result)
            with (a.output/'progress.jsonl').open('a') as f: f.write(json.dumps(result)+'\n')
            print(json.dumps(result),flush=True)
    complete=len(results)==40 and all(r['status']=='completed' for r in results)
    if any(sha(snapshot/name)!=h for name,h in sources.items()) or any(sha(package_root/name)!=h for name,h in package_sources.items()):
        complete=False
        results.append({'status':'failed','error':'frozen_or_native_package_source_changed'})
    for name in stage_names:
        audits=[]
        for control in ['paired','pairing_shuffle']:
            path=a.output/'fits'/name/control/'seed_and_background_audit.json'
            if path.is_file(): audits.append(json.loads(path.read_text()))
        if len(audits)==2 and audits[0]['backgrounds'][0]['row_indices']!=audits[1]['backgrounds'][0]['row_indices']:
            complete=False; results.append({'status':'failed','stage':name,'error':'paired_shuffle_background_mismatch'})
    config={'external_stage_links':{},'external_stage_roots':{},
            'external_method_controls':{'LINGER_native_bounded':'LINGER_native_bounded_pairing_shuffle'}}
    for control in ['paired','pairing_shuffle']:
        method='LINGER_native_bounded' if control=='paired' else 'LINGER_native_bounded_pairing_shuffle'
        mapping={}
        for r in results:
            if r['control']==control and r['status']=='completed':
                name=r['stage']; outer=name.split('_')[0].replace('outer','')
                mapping[outer+(':'+name.split('_')[1] if '_' in name else ':outer')]=r['links']
        config['external_stage_links'][method]=mapping
        config['external_stage_roots'][method]=str(a.output)
    write_json(a.output/'external_stage_config.json',config)
    write_json(a.output/'completion.json',{'status':'completed' if complete else 'incomplete',
      'scope':'DPM1_only_bounded_donor_balanced_native_training','training_fits':len(results),
      'failures':[r for r in results if r['status']!='completed'],'workers':a.workers,
      'RNA_only_control':'not_implemented','prior_only_control':'not_implemented',
      'prior_independence':'unverified','heldout_readout':'not_run_by_this_script'})
    if not complete: raise SystemExit(1)


if __name__=='__main__': main()
