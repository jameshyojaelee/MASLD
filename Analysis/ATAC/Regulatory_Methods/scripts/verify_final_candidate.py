#!/usr/bin/env python3
"""Independent rederivation of final donor losses and fixed-fold invariants."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

A=Path(__file__).resolve().parents[1]
WIDE=A/'results/20261024/ridge_grid_sensitivity_r1'
METHODS=['promoter_only','distance','cicero']
KEYS=['outer_fold','gene_id','donor_id']


def read(path,original):
    chunks=[]
    for frame in pd.read_csv(path,sep='\t',chunksize=100000,dtype={'donor_id':str}):
        keep=frame.method.isin(METHODS)
        if original:keep &= frame.grid.eq('wide')
        chunks.append(frame.loc[keep,KEYS+['method','standardized_squared_error']])
    frame=pd.concat(chunks,ignore_index=True)
    assert not frame.duplicated(KEYS+['method']).any()
    wide=frame.set_index(KEYS+['method']).standardized_squared_error.unstack('method')
    return wide.reindex(columns=METHODS).dropna()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--candidate',type=Path,required=True)
    args=p.parse_args();c=args.candidate.resolve();assert os.environ.get('SLURM_JOB_ID')
    spec=json.loads((c/'prediction_inputs_v2/specification.json').read_text())
    plan=json.loads((c/'prediction_jobs.json').read_text())
    expected=pd.read_csv(c/'prediction_summary/comparisons.tsv',sep='\t')
    expected_changes=pd.read_csv(c/'prediction_summary/paired_scenario_change_summary.tsv',sep='\t')
    baseline=pd.concat([read(WIDE/f'hepatocyte/fold{f}/predictions.tsv',True) for f in range(5)])
    maximum=0.;reconstructed=[];reconstructed_changes=[]
    for name in ['baseline']+[s['name'] for s in spec['scenarios'] if s['name']!='baseline']:
        if name=='baseline':current=baseline
        else:
            parts=[]
            for f in range(5):
                job=next(j for j in plan['jobs'] if j['scenario']==name and j['fold']==f)
                if job.get('reuse_original'):
                    parts.append(baseline[baseline.index.get_level_values('outer_fold')==f]);continue
                directory=Path(job['output'])
                record=json.loads((directory/'completion.json').read_text())
                scenario=record['scenario']
                outer=pd.read_csv(directory/'outer_folds.tsv',sep='\t',dtype={'donor_id':str})
                inner=pd.read_csv(directory/'inner_folds.tsv',sep='\t',dtype={'donor_id':str})
                assert outer.groupby('donor_id').outer_fold.nunique().eq(1).all()
                assert all(spec['outer_map'][d]==v for d,v in zip(outer.donor_id,outer.outer_fold))
                assert all(spec['inner_map'][str(f)][d]==v for d,v in zip(inner.donor_id,inner.inner_fold))
                expected_training={d for d,v in spec['outer_map'].items() if v!=f}-{scenario.get('omit_donor','')}
                assert set(inner.donor_id)==expected_training
                parts.append(read(directory/'predictions.tsv',False))
            current=pd.concat(parts)
            current=current.loc[current.index.intersection(baseline.index)]
        donors=sorted(set(current.index.get_level_values('donor_id')))
        original_current=baseline.loc[current.index]
        assert len(donors)==39
        for comparator in ['promoter_only','distance']:
            differences=[];base_losses=[]
            for d in donors:
                local=current.xs(d,level='donor_id')
                base_losses.append(float(np.mean(local[comparator].to_numpy())))
                differences.append(float(np.mean(local[comparator].to_numpy()-local.cicero.to_numpy())))
            values=np.array(differences)
            target=expected[(expected.scenario==name)&(expected.baseline==comparator)].iloc[0]
            salt=int(hashlib.sha256((name+comparator).encode()).hexdigest()[:8],16)
            rng=np.random.default_rng((42+salt)%2**32)
            bootstrap=np.array([np.mean(values[rng.integers(0,len(values),size=len(values))]) for _ in range(2000)])
            lo,hi=np.quantile(bootstrap,[.025,.975])
            actual=np.array([values.mean(),lo,hi,100*values.mean()/np.mean(base_losses)])
            saved=target[['mean_gain','ci_low','ci_high','relative_error_reduction_percent']].to_numpy(float)
            delta=float(np.max(np.abs(actual-saved)));maximum=max(maximum,delta)
            assert np.allclose(actual,saved,rtol=1e-10,atol=1e-11),(name,comparator,actual,saved)
            assert int(target.n_positive_donors)==int((values>0).sum()) and int(target.n_gene_donor_fold_pairs)==len(current)
            reconstructed.append(dict(scenario=name,baseline=comparator,mean_gain=float(values.mean()),ci_low=float(lo),ci_high=float(hi),maximum_difference=delta))
            if name!='baseline':
                paired=[]
                for d,scenario_gain in zip(donors,values):
                    local=current.xs(d,level='donor_id')
                    original=original_current.xs(d,level='donor_id').reindex(local.index)
                    original_gain=np.mean(original[comparator].to_numpy()-original.cicero.to_numpy())
                    paired.append(original_gain-scenario_gain)
                paired=np.asarray(paired)
                salt=int(hashlib.sha256((name+comparator+'paired_change').encode()).hexdigest()[:8],16)
                rng=np.random.default_rng((42+salt)%2**32)
                draws=np.array([np.mean(paired[rng.integers(0,len(paired),size=len(paired))]) for _ in range(2000)])
                lo,hi=np.quantile(draws,[.025,.975])
                target=expected_changes[(expected_changes.scenario==name)&(expected_changes.baseline==comparator)].iloc[0]
                actual=np.array([paired.mean(),lo,hi])
                saved=target[['original_gain_minus_scenario_gain','ci_low','ci_high']].to_numpy(float)
                delta=float(np.max(np.abs(actual-saved)));maximum=max(maximum,delta)
                assert np.allclose(actual,saved,rtol=1e-10,atol=1e-11),(name,comparator,'paired_change',actual,saved)
                reconstructed_changes.append(dict(scenario=name,baseline=comparator,original_gain_minus_scenario_gain=float(paired.mean()),ci_low=float(lo),ci_high=float(hi),maximum_difference=delta))
    out=c/'independent_numerical_review';out.mkdir(exist_ok=False)
    pd.DataFrame(reconstructed).to_csv(out/'reconstructed_comparisons.tsv',sep='\t',index=False)
    pd.DataFrame(reconstructed_changes).to_csv(out/'reconstructed_paired_changes.tsv',sep='\t',index=False)
    (out/'completion.json').write_text(json.dumps(dict(status='completed',n_comparisons=len(reconstructed),
        n_donors=39,maximum_difference=maximum,n_paired_changes=len(reconstructed_changes),outer_and_inner_folds_verified=True,training_only_omission_verified=True,
        conditioning='saved_fitted_predictions_not_model_refitting_bootstrap',scientific_adoption=False),indent=2))
    print('PASS: independently reconstructed',len(reconstructed),'comparisons; maximum difference',maximum,flush=True)


if __name__=='__main__':main()
