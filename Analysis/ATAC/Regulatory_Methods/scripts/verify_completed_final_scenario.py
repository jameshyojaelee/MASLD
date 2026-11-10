#!/usr/bin/env python3
"""Independent interim donor-loss reconstruction; no full robustness conclusion."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd

from verify_final_candidate import WIDE, read


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--scenario', required=True)
    parser.add_argument('--summary', type=Path, required=True)
    args = parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    saved = pd.read_csv(args.summary/'comparisons.tsv', sep='\t')
    parts = []
    for fold in range(5):
        original = read(WIDE/f'hepatocyte/fold{fold}/predictions.tsv', True)
        current = read(args.candidate/f'prediction_refits/{args.scenario}/fold{fold}/predictions.tsv', False)
        parts.append(current.loc[current.index.intersection(original.index)])
    current = pd.concat(parts)
    rows = []
    for comparator in ['promoter_only', 'distance']:
        donor_ids = sorted(set(current.index.get_level_values('donor_id')))
        assert len(donor_ids) == 39
        differences = []; baseline_losses = []
        for donor in donor_ids:
            local = current.xs(donor, level='donor_id')
            differences.append(np.mean(local[comparator].to_numpy()-local.cicero.to_numpy()))
            baseline_losses.append(np.mean(local[comparator].to_numpy()))
        values = np.array(differences)
        salt = int(hashlib.sha256((args.scenario+comparator).encode()).hexdigest()[:8], 16)
        rng = np.random.default_rng((42+salt)%2**32)
        draws = np.array([np.mean(values[rng.integers(0, len(values), size=len(values))]) for _ in range(2000)])
        lo, hi = np.quantile(draws, [.025, .975])
        result = np.array([values.mean(), lo, hi, 100*values.mean()/np.mean(baseline_losses)])
        target = saved[(saved.scenario==args.scenario)&(saved.baseline==comparator)].iloc[0]
        expected = target[['mean_gain','ci_low','ci_high','relative_error_reduction_percent']].to_numpy(float)
        assert np.allclose(result, expected, rtol=1e-10, atol=1e-11)
        assert int(target.n_positive_donors) == int((values>0).sum())
        assert int(target.n_gene_donor_fold_pairs) == len(current)
        rows.append(dict(baseline=comparator, maximum_difference=float(np.max(np.abs(result-expected)))))
    output = args.summary/f'{args.scenario}_independent_check.json'
    assert not output.exists()
    output.write_text(json.dumps(dict(status='passed', n_donors=39, checks=rows,
                                      full_robustness_conclusion_available=False), indent=2))
    print(output.read_text(), flush=True)


if __name__ == '__main__':
    main()
