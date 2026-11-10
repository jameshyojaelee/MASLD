#!/usr/bin/env python3
"""Reconstruct donor-weighted robustness results on common evaluable keys."""
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


def read_predictions(path, original=False):
    parts=[]
    for chunk in pd.read_csv(path,sep='\t',chunksize=100000,dtype={'donor_id':str}):
        chunk=chunk[chunk.method.isin(METHODS)]
        if original: chunk=chunk[chunk.grid.eq('wide')]
        parts.append(chunk)
    result=pd.concat(parts,ignore_index=True)
    assert not result.duplicated(KEYS+['method']).any()
    assert np.isfinite(result.standardized_squared_error).all()
    return result


def matched(frame):
    wide=frame.pivot(index=KEYS,columns='method',values='standardized_squared_error')
    coverage=dict(n_union=len(wide),n_common=int(wide[METHODS].notna().all(axis=1).sum()))
    return wide.dropna(subset=METHODS),coverage


def donor_losses(frame):
    return frame.groupby(level='donor_id')[METHODS].mean()


def write_gene_exclusions(frame, original, genes, scenario, output):
    """List unscored gene/fold combinations without inventing exclusion causes."""
    rows=[]
    universe=set(genes)
    for fold in range(5):
        local=frame[frame.outer_fold.eq(fold)]
        base=original[original.outer_fold.eq(fold)]
        available={method:set(local[local.method.eq(method)].gene_id) for method in METHODS}
        base_common=set.intersection(*(set(base[base.method.eq(method)].gene_id) for method in METHODS))
        assert set.union(*available.values()).issubset(universe)
        # Saved losses already average aggregates within a donor; each scored
        # gene/method must retain every held donor before gene-level coverage.
        assert local.groupby(['gene_id','method']).donor_id.nunique().eq(base.donor_id.nunique()).all()
        table=pd.DataFrame(dict(gene_id=genes,outer_fold=fold,scenario=scenario))
        for method in METHODS:table['scored_'+method]=table.gene_id.isin(available[method])
        table['scored_by_all_original_models']=table.gene_id.isin(base_common)
        common=table[['scored_'+method for method in METHODS]].all(axis=1)
        table['exclusion_reason']=np.where(~common,'not_scored_by_all_three_scenario_models',
            'not_scored_by_all_three_original_models')
        rows.append(table[~(common&table.scored_by_all_original_models)])
    pd.concat(rows,ignore_index=True).to_csv(output/f'{scenario}_excluded_genes.tsv.gz',sep='\t',index=False)


def estimate(donor,scenario,baseline,keys_count):
    values=(donor[baseline]-donor.cicero).to_numpy()
    salt=int(hashlib.sha256((scenario+baseline).encode()).hexdigest()[:8],16)
    rng=np.random.default_rng((42+salt)%2**32)
    means=np.array([rng.choice(values,len(values),replace=True).mean() for _ in range(2000)])
    lo,hi=np.quantile(means,[.025,.975])
    return dict(scenario=scenario,baseline=baseline,n_donors=len(values),n_gene_donor_fold_pairs=keys_count,
        mean_gain=float(values.mean()),ci_low=float(lo),ci_high=float(hi),
        relative_error_reduction_percent=float(100*values.mean()/donor[baseline].mean()),
        n_positive_donors=int((values>0).sum()),inference='descriptive_2000_donor_bootstrap_conditional_fixed_fits_no_new_significance_calls')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--candidate',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    args.output.mkdir(exist_ok=False)
    spec=json.loads((args.candidate/'prediction_inputs_v2/specification.json').read_text())
    jobs=json.loads((args.candidate/'prediction_jobs.json').read_text())
    baseline_parts=[]; tuning=[]
    for fold in range(5):
        part=read_predictions(WIDE/f'hepatocyte/fold{fold}/predictions.tsv',True)
        baseline_parts.append(part)
        t=part[['outer_fold','gene_id','method','ridge_alpha']].drop_duplicates()
        t['scenario']='baseline'; tuning.append(t)
    baseline= pd.concat(baseline_parts,ignore_index=True)
    genes=pd.read_csv(WIDE/'hepatocyte/prepared/genes.tsv',sep='\t').gene_id
    assert genes.is_unique and genes.notna().all()
    write_gene_exclusions(baseline,baseline,genes,'baseline',args.output)
    bwide,bcov=matched(baseline)
    bdonor=donor_losses(bwide)
    rows=[estimate(bdonor,'baseline',b,len(bwide)) for b in ['promoter_only','distance']]
    donor_rows=[bdonor.reset_index().assign(scenario='baseline',comparison_universe='within_scenario_common_three_models')]
    coverage=[dict(scenario='baseline',**bcov)]
    changes=[]; key_counts=[]
    for scenario in spec['scenarios']:
        name=scenario['name']
        if name=='baseline': continue
        parts=[]
        for fold in range(5):
            record=next(j for j in jobs['jobs'] if j['scenario']==name and j['fold']==fold)
            if record.get('reuse_original'):
                reused=baseline[baseline.outer_fold.eq(fold)].copy()
                parts.append(reused)
                t=reused[['outer_fold','gene_id','method','ridge_alpha']].drop_duplicates()
                t['scenario']=name
                tuning.append(t)
                continue
            directory=Path(record['output'])
            completion=json.loads((directory/'completion.json').read_text())
            assert completion['status']=='completed' and not completion['technical_pilot']
            assert completion['ridge_alphas']==spec['ridge_alphas']
            part=read_predictions(directory/'predictions.tsv')
            parts.append(part)
            t=pd.read_csv(directory/'tuning.tsv',sep='\t'); t['scenario']=name
            tuning.append(t)
        frame=pd.concat(parts,ignore_index=True)
        write_gene_exclusions(frame,baseline,genes,name,args.output)
        wide,cov=matched(frame); donor=donor_losses(wide)
        assert set(donor.index)==set(bdonor.index), 'A scenario lost held-out donors'
        coverage.append(dict(scenario=name,**cov))
        common=wide.index.intersection(bwide.index)
        # Paired scenario-vs-baseline changes use precisely the same evaluable keys.
        a=donor_losses(bwide.loc[common]); b=donor_losses(wide.loc[common])
        rows.extend(estimate(b,name,comparator,len(common)) for comparator in ['promoter_only','distance'])
        donor_rows.append(b.reset_index().assign(scenario=name,comparison_universe='baseline_scenario_common_keys_and_three_models'))
        key_counts.append(dict(scenario=name,n_common_keys=len(common),n_baseline_keys=len(bwide),n_scenario_keys=len(wide)))
        for comparator in ['promoter_only','distance']:
            differences=(a[comparator]-a.cicero)-(b[comparator]-b.cicero)
            for d,v in differences.items():
                changes.append(dict(scenario=name,baseline=comparator,donor_id=d,
                                    original_gain_minus_scenario_gain=float(v)))
        common.to_frame(index=False).to_csv(args.output/f'{name}_matched_keys.tsv.gz',sep='\t',index=False)
    result=pd.DataFrame(rows)
    result.to_csv(args.output/'comparisons.tsv',sep='\t',index=False)
    pd.concat(donor_rows).to_csv(args.output/'donor_scores.tsv',sep='\t',index=False)
    change_table=pd.DataFrame(changes)
    change_table.to_csv(args.output/'paired_scenario_changes.tsv',sep='\t',index=False)
    change_summary=[]
    for (scenario,comparator),group in change_table.groupby(['scenario','baseline'],sort=True):
        values=group.original_gain_minus_scenario_gain.to_numpy()
        salt=int(hashlib.sha256((scenario+comparator+'paired_change').encode()).hexdigest()[:8],16)
        rng=np.random.default_rng((42+salt)%2**32)
        samples=np.array([rng.choice(values,len(values),replace=True).mean() for _ in range(2000)])
        lo,hi=np.quantile(samples,[.025,.975])
        change_summary.append(dict(scenario=scenario,baseline=comparator,n_donors=len(values),
            original_gain_minus_scenario_gain=float(values.mean()),ci_low=float(lo),ci_high=float(hi),
            interpretation='paired_descriptive_donor_bootstrap_on_identical_gene_donor_fold_keys'))
    pd.DataFrame(change_summary).to_csv(args.output/'paired_scenario_change_summary.tsv',sep='\t',index=False)
    pd.DataFrame(coverage).to_csv(args.output/'coverage.tsv',sep='\t',index=False)
    pd.DataFrame(key_counts).to_csv(args.output/'scenario_key_intersections.tsv',sep='\t',index=False)
    tune=pd.concat(tuning,ignore_index=True)
    tune['upper_boundary']=tune.ridge_alpha.eq(1e8)
    tune.groupby(['scenario','method']).agg(n_fits=('ridge_alpha','size'),upper_boundary_fraction=('upper_boundary','mean')).reset_index().to_csv(args.output/'tuning_boundaries.tsv',sep='\t',index=False)
    unshuffled=result[~result.scenario.str.contains('shuffle')]
    robust=bool(unshuffled.mean_gain.gt(0).all())
    report=['Hepatocyte Cicero: fixed wider-grid internal robustness analysis',
            'Donors are inference units; intervals condition on fitted predictions and omit retraining uncertainty.',
            'Every scenario retains original outer and inner donor assignments; targeted omissions affect training only.',
            'Prediction retains the frozen wider-grid signed-nonzero Cicero selection rule. Positive-only disease-cohort annotation maps are a separate descriptive analysis.',
            'Pairing shuffles preserve donor RNA distributions; within-donor aggregate targets and training variance may change.',
            'No new confirmatory P values or multiplicity family. Existing Holm results remain separate.',
            'Robust internal predictive support: '+str(robust), '']
    for row in result.to_dict('records'):
        report.append(f"{row['scenario']} / Cicero vs {row['baseline']}: gain {row['mean_gain']:.8g} [{row['ci_low']:.8g}, {row['ci_high']:.8g}], {row['relative_error_reduction_percent']:.4f}% lower MSE; {row['n_positive_donors']}/{row['n_donors']} positive donors.")
    report += ['', 'Pairing controls: original minus shuffled Cicero gain on identical evaluable keys. Positive values favor correct pairing.']
    for row in change_summary:
        if row['scenario'].startswith('pairing_shuffle_'):
            report.append(f"{row['scenario']} / vs {row['baseline']}: paired-minus-shuffled gain {row['original_gain_minus_scenario_gain']:.8g} [{row['ci_low']:.8g}, {row['ci_high']:.8g}].")
    shuffle_distance=[row for row in change_summary if row['scenario'].startswith('pairing_shuffle_') and row['baseline']=='distance']
    assert len(shuffle_distance)==3
    if all(row['original_gain_minus_scenario_gain']<=0 for row in shuffle_distance):
        report.append('Correct pairing does not increase the point-estimate Cicero gain over distance in any of the three shuffles. These comparisons do not demonstrate added value from correct nucleus-level pairing for Cicero beyond proximity.')
    report.append('Shuffling preserves donor marginals, so persistence of a gain does not rule out donor-level chromatin-RNA association.')
    report += ['', 'This cohort has no qualified MASH contrast and is already used for development. No independent MASLD validation, causal regulation, or individual-edge validation follows.']
    (args.output/'REPORT.txt').write_text('\n'.join(report)+'\n')
    (args.output/'completion.json').write_text(json.dumps(dict(status='completed',robust_internal_predictive_support=robust,
        n_scenarios=result.scenario.nunique(),scientific_adoption=False),indent=2))


if __name__=='__main__': main()
