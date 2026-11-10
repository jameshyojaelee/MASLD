#!/usr/bin/env python3
"""Descriptive interim review of scenarios with all five folds complete."""
import argparse
import json
import os
from pathlib import Path

import pandas as pd

from final_prediction_summary import WIDE, donor_losses, estimate, matched, read_predictions


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    args.output.mkdir(exist_ok=False)
    plan = json.loads((args.candidate/'prediction_jobs.json').read_text())
    groups = {}
    for job in plan['jobs']:
        groups.setdefault(job['scenario'], []).append(job)
    rows = []; loss_rows = []; coverage = []
    for name, jobs in groups.items():
        ready = all(j.get('reuse_original') or (Path(j['output'])/'completion.json').exists() for j in jobs)
        if len(jobs) != 5 or not ready:
            continue
        scenario_parts = []; baseline_parts = []
        for job in sorted(jobs, key=lambda j: j['fold']):
            fold = job['fold']
            base, _ = matched(read_predictions(WIDE/f'hepatocyte/fold{fold}/predictions.tsv', True))
            if job.get('reuse_original'):
                current = base
            else:
                record = json.loads((Path(job['output'])/'completion.json').read_text())
                assert record['status'] == 'completed' and not record['technical_pilot']
                current, _ = matched(read_predictions(Path(job['output'])/'predictions.tsv'))
            common = base.index.intersection(current.index)
            coverage.append(dict(scenario=name, outer_fold=fold, n_matched_keys=len(common)))
            scenario_parts.append(current.loc[common])
            baseline_parts.append(base.loc[common])
        current = pd.concat(scenario_parts)
        original = pd.concat(baseline_parts)
        losses = donor_losses(current)
        original_losses = donor_losses(original)
        assert len(losses) == 39 and losses.index.equals(original_losses.index)
        for comparator in ['promoter_only', 'distance']:
            row = estimate(losses, name, comparator, len(current))
            row['original_gain_on_same_keys'] = float((original_losses[comparator]-original_losses.cicero).mean())
            rows.append(row)
        loss_rows.append(losses.reset_index().assign(scenario=name))
    assert rows, 'No full scenario is ready'
    pd.DataFrame(rows).to_csv(args.output/'comparisons.tsv', sep='\t', index=False)
    pd.concat(loss_rows).to_csv(args.output/'donor_losses.tsv', sep='\t', index=False)
    pd.DataFrame(coverage).to_csv(args.output/'coverage.tsv', sep='\t', index=False)
    (args.output/'completion.json').write_text(json.dumps(dict(status='completed', interim=True,
        n_scenarios=len(loss_rows), full_robustness_conclusion_available=False), indent=2))
    print(pd.DataFrame(rows).to_string(index=False), flush=True)


if __name__ == '__main__':
    main()
