"""Independent compute-only reconstruction of saved program/PCA deviance.

Only the completed99-participant development outputs and their fixed H3
fixture are read. No producer functions, model fit, RNA or receiving data.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import shutil
import time
import traceback

import numpy as np


ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
SOURCE = REC/'codex_program_relative_h3_21999301'
FIX = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture'
FIX_SHA = '612f0cb079ba611157663799c7cbb312aeedf8f32b25cba3ddaf171e1f9fbc29'
SOURCE_HASHES = {
    'COMPLETE': '81b2bd4ea98c8db66554fbc8d7637a1a69a130f331feb732b75caab4c4868fd5',
    'source_receipt.json': '2f38c1df8563d960a15ad84080c3c47bfe0fc6d792533028a0e8ece6b6e434f7',
    'summary.json': '4dd320f5fbbe5453774862d8061be508cedd650b9f7896bf01a94d4a47f34c41',
    'development_results.json': '5ed7e152018c000776da992409d6da9dcadd5cfa9ef61764a00878357cd32c6a',
    'all99_donor_results.tsv': 'c44ae31c8c5abf2d6aa3f1d2e81d227383fa79b8f0f68f3367607946919b0ac7',
    'oof_program_profiles.npy': '4d38849d67e55305762a97ffcc7df6620fa63fd9f8681156bddfd74658b8bde8',
    'oof_pca_profiles.npy': '177b5717a755ab54d2a7dd238ea77ba65dde41963d3ca94defbf1ca606058f12',
    'participant_outer_folds.tsv': '9ddfaf9c58b58a5b239eebbbc19b214210543fb4f60f3f5a7c563a70fbebb0d1',
    'participant_axis.tsv': '7e5be1a16df505e5daad451398e8b1a4b3cbbc48a32cd1381113e1ae4b7347bc',
    'h3k27ac_feature_axis.tsv': 'cbe35aeb188e531283a1bcd636ba14de4cfbf215dfc544b99e3ebb97f9831986',
    'codex_program_relative_h3.py': '04c23621c948fd3cb9002cb5520fa6cd05602b1e673bd1b2fb7d464085305727',
    'prospective_protocol.json': '5f14bd9d84bbc353c67378df8173602fd3bc7825cc72fc1a82175bfafd03342c',
    'fold_0_receipt.json': 'b9e9b7e6e19170dc50a3f0a22522a5e67b34474791c40703d166d59aa990e45c',
    'fold_1_receipt.json': '50bab6026684b951d8632012f0c07bd27d539cb9e3ea95d0c4302de1ff306dd6',
    'fold_2_receipt.json': '4da58c96b41f5da7eed6400af6849b122a6ca5084c176fff08e2796098150f93',
    'fold_3_receipt.json': 'f0444c59bf67dce98f6af3f793e668016642c91e1cf78887c8d5850a5d9a3acf',
    'fold_4_receipt.json': '6b3e6e4f6fdc17d3c185d5ce58cf5d7ee7f3544ecd8fec4238ac70cb5b1fc64a',
}
FIX_HASHES = {
    'molecular/ARTIFACTS.json': '2a60ef5d5d904f5f833833ebef75e7946e7db3b6b9aaf502cb84ddedf97eb456',
    'folds/ARTIFACTS.json': '9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3',
    'molecular/h3k27ac_counts.npy': '9a4fa736a6dd5de0cd44e51428f2bf19e8dae9adb7a279e5b6c8f1c14103cd22',
    'molecular/h3k27ac_feature_axis.tsv': SOURCE_HASHES['h3k27ac_feature_axis.tsv'],
    'molecular/participant_axis.tsv': SOURCE_HASHES['participant_axis.tsv'],
    'folds/participant_outer_folds.tsv': SOURCE_HASHES['participant_outer_folds.tsv'],
}
STATE = {'phase': 'startup', 'out': None}


def require(okay, message):
    if not okay:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1024*1024), b''):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(reader.fieldnames and len(set(reader.fieldnames)) == len(reader.fieldnames), 'Bad TSV header')
        rows = list(reader)
    require(all(None not in row and all(v is not None for v in row.values()) for row in rows), 'Bad TSV width')
    return rows


def save(name, value):
    with (STATE['out']/(name+'.json')).open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def validate_hashes():
    require(sha(FIX/'ARTIFACTS.json') == FIX_SHA, 'Development fixture manifest changed')
    for name, digest in SOURCE_HASHES.items():
        require(sha(SOURCE/name) == digest, 'Completed output changed: '+name)
    manifest = json.loads((FIX/'ARTIFACTS.json').read_text())
    indexed = {r['path']:r['sha256'] for r in manifest['artifacts']}
    for name, digest in FIX_HASHES.items():
        require(indexed.get(name) == digest and sha(FIX/name) == digest, 'Fixed H3/axis/fold artifact changed: '+name)


def profile_floor(probability):
    # Additive pseudomass on probabilities, then renormalize, not max-clipping
    # counts or adding one read. The same fixed floor applies to both baselines.
    x = probability+1e-12
    return x/x.sum()


def divergence(q, total, prediction):
    require(np.isfinite(prediction).all() and (prediction > 0).all()
            and abs(prediction.sum()-1) <= 1e-12, 'Invalid saved/baseline probability profile')
    use = q > 0
    return float(2*total*np.sum(q[use]*(np.log(q[use])-np.log(prediction[use]))))


def main():
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Compute allocation required')
    out = REC/('codex_program_deviance_verify_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    STATE.update(out=out, phase='input_hashes')
    started = time.perf_counter()
    validate_hashes()
    wrapper = Path(__file__).resolve()
    launcher = wrapper.with_name('run_codex_program_deviance_verify.sbatch')
    owned = {str(p):sha(p) for p in (wrapper,launcher)}
    for path in (wrapper,launcher): shutil.copyfile(path,out/path.name)
    save('protocol', dict(source=str(SOURCE), fixture=str(FIX), source_sha256=SOURCE_HASHES,
                          fixture_manifest_sha256=FIX_SHA, opened_fixture_sha256=FIX_HASHES,
                          owned_sha256=owned, seed=20260930, bootstrap_replicates=10000,
                          quantile_method='linear', participants=99, features=96460,
                          formula='D=2*N*sum(q>0) q*log(q/p); q=C/N',
                          denominator='Fold-training pooled count MLE; additive1e-12 probability floor then normalize',
                          deviance_atol=1e-6,deviance_rtol=1e-12,skill_ci_atol=1e-12,
                          producer_functions_imported=False, fits=False, receiving_values=False,
                          python=platform.python_version(), numpy=np.__version__))
    summary = json.loads((SOURCE/'summary.json').read_text())
    expected = json.loads((SOURCE/'development_results.json').read_text())
    require(summary['status']=='success' and summary['participants']==99
            and summary['receiving_evaluated'] is False and (SOURCE/'COMPLETE').read_text()=='success\n',
            'Producer did not complete fixed development experiment')
    require(summary['development_results_sha256']==SOURCE_HASHES['development_results.json']
            and expected['oof_program_sha256']==SOURCE_HASHES['oof_program_profiles.npy']
            and expected['oof_pca_sha256']==SOURCE_HASHES['oof_pca_profiles.npy'], 'OOF receipt chain differs')
    STATE['phase']='identity_and_fold_join'
    people=read_tsv(FIX/'molecular/participant_axis.tsv')
    folds=read_tsv(FIX/'folds/participant_outer_folds.tsv')
    axis=read_tsv(FIX/'molecular/h3k27ac_feature_axis.tsv')
    saved=read_tsv(SOURCE/'all99_donor_results.tsv')
    ids=[r['participant_id'] for r in people]
    require(len(ids)==99 and len(set(ids))==99
            and [int(r['participant_index']) for r in people]==list(range(99)), 'Participant axis not unique ordered99')
    require(ids==[r['participant_id'] for r in folds]==[r['participant_id'] for r in saved]
            and [int(r['participant_index']) for r in saved]==list(range(99)), 'OOF/target participant order differs')
    foldids=np.asarray([int(r['outer_fold']) for r in folds])
    require(np.isin(foldids,np.arange(5)).all()
            and [int((foldids==f).sum()) for f in range(5)]==[21,21,21,19,17]
            and np.array_equal(foldids,[int(r['outer_fold']) for r in saved]), 'Fixed folds differ')
    require(len(axis)==96460 and [int(r['h3k27ac_feature_index']) for r in axis]==list(range(96460))
            and len({r['opaque_source_feature_key'] for r in axis})==96460,'Target axis differs')
    for f in range(5):
        receipt=json.loads((SOURCE/('fold_'+str(f)+'_receipt.json')).read_text())
        require(receipt['status']=='success' and receipt['fold']==f
                and receipt['training_indices']==np.flatnonzero(foldids!=f).tolist()
                and receipt['query_indices']==np.flatnonzero(foldids==f).tolist(), 'Fold training/query receipt differs')
    STATE['phase']='development_arrays'
    counts=np.load(FIX/'molecular/h3k27ac_counts.npy',mmap_mode='r',allow_pickle=False)
    predictions={name:np.load(SOURCE/('oof_'+name+'_profiles.npy'),mmap_mode='r',allow_pickle=False)
                 for name in ('program','pca')}
    require(counts.shape==(99,96460) and counts.dtype==np.uint32,'Fixed native H3 count units/shape differ')
    require(all(p.shape==counts.shape and p.dtype==np.float64 for p in predictions.values()),'Saved OOF fullfloat64 shapes differ')
    total=counts.sum(axis=1,dtype=np.uint64)
    require((total>0).all(),'Zero target library; no donor may be dropped')
    q=np.asarray(counts,dtype=np.float64)/total[:,None]
    metric_rows=[]
    vectors={name:np.empty(99,dtype=np.float64) for name in
             ('program_deviance','pca_deviance','pooled_mle_deviance','macro_mean_deviance')}
    STATE['phase']='independent_deviance'
    for f in range(5):
        train=foldids!=f
        pooled=counts[train].sum(axis=0,dtype=np.float64)
        pooled=profile_floor(pooled/pooled.sum())
        macro=profile_floor(q[train].mean(axis=0))
        for i in np.flatnonzero(foldids==f):
            for name,p in predictions.items():
                vectors[name+'_deviance'][i]=divergence(q[i],total[i],p[i])
            vectors['pooled_mle_deviance'][i]=divergence(q[i],total[i],pooled)
            vectors['macro_mean_deviance'][i]=divergence(q[i],total[i],macro)
    require(all(np.isfinite(v).all() for v in vectors.values())
            and (vectors['pooled_mle_deviance']>0).all(),'Undefined native normalized deviance')
    denominator=vectors['pooled_mle_deviance']
    vectors['program_skill']=1-vectors['program_deviance']/denominator
    vectors['pca_skill']=1-vectors['pca_deviance']/denominator
    vectors['normalized_program_minus_pca_skill']=(vectors['pca_deviance']-vectors['program_deviance'])/denominator
    differences={}
    for name,actual in vectors.items():
        reference=np.asarray([float(r[name]) for r in saved])
        require(np.isfinite(reference).all(),'Saved metric nonfinite')
        differences[name]=float(np.max(np.abs(actual-reference)))
    save('metric_max_absolute_differences',differences)
    for name,actual in vectors.items():
        atol=1e-6 if name.endswith('_deviance') else 1e-12
        require(np.allclose(actual,[float(r[name]) for r in saved],rtol=1e-12,atol=atol),'Metric reconstruction differs: '+name)
    require(np.array_equal(total,[int(r['h3_library_total_evaluator_only']) for r in saved]),'Exact evaluator count totals differ')
    for i in range(99):
        metric_rows.append(dict(participant_index=i,participant_id=ids[i],outer_fold=int(foldids[i]),
                                **{name:float(v[i]) for name,v in vectors.items()}))
    with (out/'reconstructed_all99_metrics.tsv').open('x',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(metric_rows[0]),delimiter='\t')
        writer.writeheader();writer.writerows(metric_rows)
    STATE['phase']='paired_bootstrap'
    delta=vectors['normalized_program_minus_pca_skill']
    rng=np.random.default_rng(20260930)
    samples=delta[rng.integers(0,99,size=(10000,99))].mean(axis=1)
    ci=np.quantile(samples,[.025,.975],method='linear')
    result=dict(primary_mean_skill_difference=float(delta.mean()),paired_conditional_percentile95=ci.tolist(),
                candidate_mean_pooled_skill=float(vectors['program_skill'].mean()),
                pca_mean_pooled_skill=float(vectors['pca_skill'].mean()))
    result['statistical_advancement_rule_passes']=bool(ci[0]>0 and result['candidate_mean_pooled_skill']>0)
    save('reconstructed_development_results',result)
    for name in ('primary_mean_skill_difference','paired_conditional_percentile95',
                 'candidate_mean_pooled_skill','pca_mean_pooled_skill'):
        require(np.allclose(result[name],expected[name],rtol=0,atol=1e-12),'Development summary differs: '+name)
    require(expected['bootstrap_seed']==20260930 and expected['bootstrap_replicates']==10000
            and result['statistical_advancement_rule_passes']==expected['statistical_advancement_rule_passes'],
            'Fixed bootstrap or descriptive decision differs')
    STATE['phase']='postread_hashes'
    validate_hashes()
    require(all(sha(Path(p))==d for p,d in owned.items()),'Verification source changed')
    save('summary',dict(status='success',participants=99,all_metrics_reconstructed=True,
                        paired_bootstrap_reconstructed=True,elapsed_seconds=time.perf_counter()-started,
                        max_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                        candidate_fit_or_tuning=False,external_validation=False,
                        qualification='Conditional development uncertainty for saved cross-validation predictions; overlapping training folds, retraining and source/assay uncertainty are not included.'))
    print('success: saved99 development metrics independently reconstructed; no fit or receiver')


if __name__=='__main__':
    try:
        main()
    except Exception:
        if STATE['out'] is not None:
            save('failure',dict(phase=STATE['phase'],traceback=traceback.format_exc(),partials_retained=True))
        raise
