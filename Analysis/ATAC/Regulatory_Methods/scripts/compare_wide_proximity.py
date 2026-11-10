"""Compare wide-grid link readouts on the three-method gene intersection."""
import argparse
import importlib.util
import json
import os
from pathlib import Path
import numpy as np
import pandas as pd

p = argparse.ArgumentParser()
p.add_argument('--root', type=Path, required=True)
p.add_argument('--output', type=Path, required=True)
a = p.parse_args()
assert os.environ.get('SLURM_JOB_ID'), 'Compute allocation required'
assert not a.output.exists(), 'Refusing overwrite'
spec = importlib.util.spec_from_file_location('check', a.root / 'summary_source/check_wide_ridge_results.py')
check = importlib.util.module_from_spec(spec)
spec.loader.exec_module(check)
# Reuse the independently verified bootstrap implementation with a contrast-specific seed.
original_seed = check.seed_for
check.seed_for = lambda method: original_seed(method + '_versus_distance_common_three_methods')
methods = ['distance', 'cicero', 'directnet_adapted']
results, donor_rows, coverage, membership = [], [], [], []
for lineage, expected in [('hepatocyte', 39), ('stellate', 18)]:
    pieces = []
    for fold in range(5):
        directory = a.root / lineage / f'fold{fold}'
        assert json.loads((directory / 'completion.json').read_text())['status'] == 'completed'
        fold_n = 0
        for frame in check.blocks(directory / 'predictions.tsv'):
            frame = frame[frame.grid.eq('wide') & frame.method.isin(methods)]
            assert frame.outer_fold.eq(fold).all() and frame.endpoint.eq('rna').all()
            keys = ['gene_id', 'donor_id']
            assert not frame.duplicated(keys + ['method']).any()
            losses = frame.pivot(index=keys, columns='method', values='standardized_squared_error').dropna()
            if losses.empty:
                continue
            # All three methods must score the identical target for every retained pair.
            for field in ['observed', 'n_metacells', 'target_unit']:
                target = frame.pivot(index=keys, columns='method', values=field).loc[losses.index]
                assert target.nunique(axis=1).eq(1).all(), field
            assert np.isfinite(losses.to_numpy()).all()
            assert (losses.to_numpy() >= 0).all()
            genes = sorted(set(losses.index.get_level_values('gene_id')))
            held = set(frame.donor_id)
            assert losses.groupby(level='gene_id').size().eq(len(held)).all()
            fold_n += len(genes)
            membership.extend(dict(lineage=lineage, outer_fold=fold, gene_id=g) for g in genes)
            losses = losses.reset_index()
            losses['outer_fold'] = fold
            pieces.append(losses)
        coverage.append(dict(lineage=lineage, outer_fold=fold, common_genes=fold_n))
        print(lineage, fold, fold_n, flush=True)
    pairs = pd.concat(pieces, ignore_index=True)
    assert not pairs.duplicated(['gene_id', 'donor_id']).any()
    assert pairs.groupby('donor_id').outer_fold.nunique().eq(1).all()
    scores = pairs.groupby('donor_id')[methods].mean()
    assert len(scores) == expected
    for method in methods[1:]:
        difference = scores.distance - scores[method]
        stat = check.inference(difference.to_numpy(), method)
        # Independent arithmetic route for the donor mean differences.
        direct = (pairs.distance - pairs[method]).groupby(pairs.donor_id).mean()
        assert np.allclose(difference, direct, rtol=1e-10, atol=1e-12)
        results.append(dict(lineage=lineage, method=method, baseline='distance',
                            n_donors=expected, common_genes=pairs.gene_id.nunique(),
                            common_gene_donor_pairs=len(pairs),
                            baseline_mse=scores.distance.mean(), candidate_mse=scores[method].mean(),
                            relative_mse_reduction=stat['improvement']/scores.distance.mean(),
                            positive_donors=int((difference > 0).sum()), **stat))
        for donor, value in difference.items():
            donor_rows.append(dict(lineage=lineage, method=method, donor_id=donor,
                                   gain=value, distance_mse=scores.loc[donor, 'distance'],
                                   candidate_mse=scores.loc[donor, method]))
    del pairs, pieces
for row, adjusted in zip(results, check.holm([r['pvalue'] for r in results])):
    row['holm_pvalue'] = adjusted
    row['supported'] = row['ci_low'] > 0 and adjusted < .05
a.output.mkdir(parents=True)
for name, rows in [('comparisons', results), ('donor_scores', donor_rows), ('fold_coverage', coverage), ('common_genes', membership)]:
    pd.DataFrame(rows).to_csv(a.output / (name + '.tsv'), sep='\t', index=False)
lines = ['Wider-grid RNA prediction beyond genomic proximity', '',
         'All three methods use the identical gene intersection within each original donor fold.',
         'Both comparisons use the existing wider-grid fits; no refitting or new tuning.',
         'Positive gain is distance minus candidate MSE in training RNA-variance units.',
         'Donors are the inference units. Bootstrap: 2,000; sign flips: 10,000; seed: 42 with contrast-specific deterministic seeds; Holm family: four.', '']
for r in results:
    lines.append(f"{r['lineage']} {r['method']}: gain {r['improvement']:.9g} [{r['ci_low']:.9g}, {r['ci_high']:.9g}], relative reduction {100*r['relative_mse_reduction']:.4f}%, Holm p={r['holm_pvalue']:.9g}; n={r['n_donors']}, genes={r['common_genes']}, positive donors={r['positive_donors']}; supported={r['supported']}.")
lines += ['', 'This is a post hoc internal comparison on previously evaluated donors, conditional on fitted cross-validation predictions. It omits model-refitting uncertainty; sign flips assume independent symmetric donor differences. Mesenchymal is the broad source population, not purified stellate cells.',
          'Beating distance does not establish improvement over promoter-only prediction, independent MASH replication or causal regulatory links. The separate wider-grid promoter comparisons remain unchanged.']
(a.output / 'REPORT.txt').write_text('\n'.join(lines) + '\n')
(a.output / 'completion.json').write_text(json.dumps(dict(status='completed', family_size=4, same_three_method_gene_universe=True, seed=42, bootstrap_reps=2000, signflip_reps=10000, source=str(a.root.resolve()), scientific_adoption=False), indent=2)+'\n')
