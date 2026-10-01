#!/usr/bin/env python3
"""Reconstruct B2 error from source counts and saved predictions independently."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--b2', type=Path, required=True)
    p.add_argument('--inference', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Run scientific checks on a compute node')
    a.out.mkdir(parents=True, exist_ok=False)
    fixture = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture'
    counts = np.load(fixture/'molecular/h3k27ac_counts.npy', mmap_mode='r')
    axis = pd.read_csv(fixture/'molecular/participant_axis.tsv', sep='\t')
    folds = pd.read_csv(fixture/'folds/participant_outer_folds.tsv', sep='\t').set_index('participant_id').loc[axis.participant_id, 'outer_fold'].to_numpy()
    regions = pd.read_csv(a.b2/'regions.tsv', sep='\t')
    observed = np.log2(1 + counts[:, regions.region_index].astype(float) * 1e6 / np.sum(counts, axis=1, dtype=np.float64)[:, None])
    metrics = pd.read_csv(a.b2/'metrics.tsv', sep='\t')
    predictions = np.load(a.b2/'internal_predictions.npz', allow_pickle=False)

    def pair_error(error):
        # Enumerate real distinct donor pairs. Pair counts are not biological n.
        total = 0.
        for f in np.unique(folds):
            ids = np.flatnonzero(folds == f)
            first, second = np.triu_indices(len(ids), k=1)
            total += len(ids)/len(folds) * np.square(error[ids[first]] - error[ids[second]]).mean()
        return total

    checked = []
    for row in metrics.itertuples():
        keep = np.ones(len(regions), dtype=bool) if row.evaluation == 'held_donors_trained_regions' else (regions.region_role == 'held').to_numpy()
        target = observed[:, keep]
        pred = predictions[row.evaluation+'|'+row.model].astype(float)
        mse = np.square(pred-target).mean()
        pair_mse = pair_error(pred-target)
        skill = 1-pair_mse/pair_error(target)
        np.testing.assert_allclose(mse, row.mse_log2cpm, rtol=1e-5, atol=2e-6)
        np.testing.assert_allclose(pair_mse, row.donor_difference_mse_log2cpm, rtol=1e-5, atol=2e-6)
        np.testing.assert_allclose(skill, row.donor_difference_skill_vs_zero, rtol=1e-5, atol=2e-6)
        checked.append(dict(evaluation=row.evaluation, model=row.model, mse_log2cpm=mse,
                            donor_difference_mse_log2cpm=pair_mse, donor_difference_skill_vs_zero=skill))
    pd.DataFrame(checked).to_csv(a.out/'B2_independently_reconstructed.tsv', sep='\t', index=False)
    inference = {}
    for name in ('variant', 'haplotype'):
        result = json.loads((a.inference/(name+'.json')).read_text())
        identities = [(r['assay'], r['track']) for r in result['effects']]
        assert identities and len(identities) == len(set(identities))
        assert result['repeat_max_abs_track_change'] <= 1e-5
        for r in result['effects']:
            expected = float(np.float32(r['alt_value'])-np.float32(r['ref_value']))
            assert r['effect'] == expected and r['kind'] == 'predicted'
            assert r['track_metadata']['name'] == r['track']
        inference[name] = dict(tracks=len(identities), variants=len(result['variants']),
                               repeat_max_abs=result['repeat_max_abs_track_change'],
                               identities_unique=True, fp32_subtraction_matches=True)
    report = dict(status='pass', B2_models_reconstructed=len(checked), biological_n=99,
                  B2_native_target='log2(1+CPM), full source-library denominator',
                  B2_error='explicit within-fold donor pairs, donor-balanced fold weights',
                  saved_prediction_precision='float32; error tolerance 2e-6 absolute plus 1e-5 relative',
                  inference=inference, protected_outcomes_read=False)
    (a.out/'checks.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
