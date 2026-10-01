#!/usr/bin/env python3
"""Reconstruct the common representation population and its reported metrics."""
import argparse
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[3]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--frozen', type=Path, required=True)
    p.add_argument('--out', type=Path, required=True)
    a = p.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        raise RuntimeError('Scientific checks require an allocated compute node')
    a.out.mkdir(parents=True, exist_ok=False)
    source = ROOT/'GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz'
    labels = pd.read_csv(source, sep='\t').set_index('lead_variant_id')
    manifest = pd.read_csv(a.frozen/'2048/manifest.tsv', sep='\t')
    other = pd.read_csv(a.frozen/'16384/manifest.tsv', sep='\t')
    assert manifest.key.tolist() == other.key.tolist()
    keep = np.ones(len(manifest), dtype=bool)
    coverage = []
    for length in (2048,16384):
        with np.load(a.frozen/str(length)/'coverage.npz', allow_pickle=False) as z:
            for pooling in ('variant','symmetric','target'):
                assert z[pooling].dtype == bool and len(z[pooling]) == len(manifest)
                keep &= z[pooling]
                coverage.append(dict(length=length,pooling=pooling,eligible_rows=int(z[pooling].sum())))
    expected = manifest.loc[keep].copy()
    observed = pd.read_csv(a.frozen/'comparisons/predictions.tsv.gz', sep='\t')
    assert observed.key.tolist() == expected.key.tolist()
    assert not observed.key.duplicated().any()
    actual_labels = labels.loc['chr'+observed.key]
    np.testing.assert_allclose(observed.beta_alt, actual_labels.beta_alt, rtol=1e-6, atol=1e-7)
    assert observed.heldout_fold.tolist() == actual_labels.heldout_fold.tolist()
    assert observed.block_1mb.tolist() == actual_labels.block_1mb.tolist()
    assert observed.groupby('block_1mb').heldout_fold.nunique().max() == 1
    reported = pd.read_csv(a.frozen/'comparisons/performance.tsv', sep='\t')
    failures = pd.read_csv(a.frozen/'comparisons/failed_fits.tsv', sep='\t')
    assert len(reported) == 12 and len(failures) == 0
    checked = []
    for r in reported.itertuples():
        pred = observed[r.recipe].to_numpy(float)
        assert np.isfinite(pred).all()
        # Recreate the producer's documented float32 label/prediction storage.
        target = observed.beta_alt.to_numpy(np.float32)
        rmse = np.sqrt(np.mean((pred.astype(np.float32)-target)**2))
        ranks = [spearmanr(pred[observed.heldout_fold == f], target[observed.heldout_fold == f]).statistic for f in range(5)]
        # The established endpoint averages correlations in Fisher-z space.
        macro = np.tanh(np.mean(np.arctanh(np.clip(ranks,-0.999999,0.999999))))
        np.testing.assert_allclose(rmse, r.RMSE_beta_units, rtol=1e-6, atol=1e-7)
        np.testing.assert_allclose(macro, r.macro_spearman, rtol=1e-6, atol=1e-7)
        checked.append(dict(recipe=r.recipe,macro_spearman=macro,RMSE_beta_units=float(rmse)))
    pd.DataFrame(checked).to_csv(a.out/'independent_performance.tsv',sep='\t',index=False)
    report = dict(status='pass',source_label_rows=len(labels),matched_rows=len(observed),
                  blocks=observed.block_1mb.nunique(),recipes_checked=len(checked),coverage=coverage,
                  biological_unit='source participant-based caQTL estimates; prediction comparison by locus',
                  orientation='source ALT-dosage convention checked separately in data evidence record',
                  macro_definition='tanh(mean(arctanh(per_fold_Spearman)))',
                  native_comparison='pending complete matched native scores',protected_outcomes_read=False)
    (a.out/'checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__ == '__main__':
    main()
