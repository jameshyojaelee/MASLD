#!/usr/bin/env python3
"""Preserve fitted maps and derive the prespecified positive-only annotations."""
import argparse
import json
import os
from pathlib import Path

import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--candidate', type=Path, required=True)
    args = parser.parse_args()
    assert os.environ.get('SLURM_JOB_ID')
    target = args.candidate / 'disease_links_positive'
    target.mkdir(exist_ok=False)
    for source in sorted((args.candidate / 'disease_links').iterdir()):
        if not source.is_dir():
            continue
        output = target / source.name
        output.mkdir()
        for path in source.iterdir():
            if path.name not in ['links.tsv.gz', 'coverage.tsv', 'completion.json']:
                (output / path.name).symlink_to(path.resolve())
        links = pd.read_csv(source / 'links.tsv.gz', sep='\t')
        keep = links.method.ne('cicero') | links.score.gt(0)
        selected = links.loc[keep].copy()
        assert selected.loc[selected.method.eq('cicero'), 'score'].gt(0).all()
        assert selected.groupby(['gene_id', 'method']).size().le(10).all()
        selected.to_csv(output / 'links.tsv.gz', sep='\t', index=False)
        coverage = pd.read_csv(source / 'coverage.tsv', sep='\t')
        counts = selected[selected.method.eq('cicero')].groupby('gene_id').size()
        coverage['cicero_links'] = coverage.gene_id.map(counts).fillna(0).astype(int)
        assert coverage.cicero_links.sum() == selected.method.eq('cicero').sum()
        coverage.to_csv(output / 'coverage.tsv', sep='\t', index=False)
        record = json.loads((source / 'completion.json').read_text())
        record.update(n_links=len(selected), source=str(source.resolve()),
                      n_nonpositive_cicero_excluded=int((~keep).sum()),
                      cicero_annotation_rule='top_ten_positive_only',
                      prediction_feature_selection_changed=False)
        (output / 'completion.json').write_text(json.dumps(record, indent=2))
        print(source.name, record['n_links'], record['n_nonpositive_cicero_excluded'], flush=True)


if __name__ == '__main__':
    main()
