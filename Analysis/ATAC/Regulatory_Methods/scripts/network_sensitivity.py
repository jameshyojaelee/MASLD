#!/usr/bin/env python3
"""Donor-preserving network sensitivity refits and matched baseline comparisons.

Default writes a reviewable run plan. --execute runs sequentially on compute.
Pairing shuffle changes RNA cell identities within donor/lineage before aggregation;
it preserves donor expression distributions and tests nucleus-level pairing value.
"""
import argparse
import hashlib
import importlib.util
import json
import csv
import sqlite3
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse


def load_network():
    path = Path(__file__).with_name('network_compare.py')
    spec = importlib.util.spec_from_file_location('regulatory_network', path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def shuffle_index(cells, seed):
    """Permutation never crosses cohort, lineage or biological donor."""
    index = np.arange(len(cells))
    for key, rows in cells.groupby(['cohort', 'lineage', 'donor_id'], sort=True).indices.items():
        offset = int(hashlib.sha256(str(key).encode()).hexdigest()[:8], 16)
        index[rows] = np.random.default_rng((seed + offset) % 2**32).permutation(rows)
    return index


def fold_map(path):
    frame = pd.read_csv(path, sep='\t', dtype={'donor_id': str})
    if (frame.groupby('donor_id').outer_fold.nunique() != 1).any():
        raise ValueError('A baseline donor occurs in multiple outer folds')
    return frame.groupby('donor_id').outer_fold.first().to_dict()


def compare_tables(baseline, candidate, output, seed=42, bootstrap_reps=2000):
    """Disk-backed joins bound RAM; bootstrap conditions on fixed fitted predictions."""
    output.mkdir(exist_ok=False)
    rows_per_chunk = 50000
    with tempfile.TemporaryDirectory(prefix='sensitivity-', dir=output) as temporary:
        connection = sqlite3.connect(str(Path(temporary) / 'comparison.sqlite'))
        connection.execute('PRAGMA cache_size=-32768')
        connection.execute('PRAGMA temp_store=FILE')
        for kind, columns in [('links', ['outer_fold','gene_id','method','peak_id','score']), ('predictions', ['outer_fold','gene_id','endpoint','method','donor_id','squared_error'])]:
            for tag, root in [('a', baseline), ('b', candidate)]:
                for chunk in pd.read_csv(root / f'{kind}.tsv', sep='\t', usecols=columns, chunksize=rows_per_chunk, dtype={'donor_id':str, 'gene_id':str}):
                    chunk.to_sql(kind + tag, connection, if_exists='append', index=False)
                keys = ['outer_fold','gene_id','method','peak_id'] if kind == 'links' else ['outer_fold','gene_id','endpoint','method','donor_id']
                connection.execute(f"CREATE UNIQUE INDEX {kind}{tag}_keys ON {kind}{tag} ({','.join(keys)})")
            connection.commit()
        refitted_methods = {row[0] for row in connection.execute('SELECT DISTINCT method FROM predictionsb')}
        edge_fields = ['outer_fold','gene_id','method','baseline_edges','sensitivity_edges','common_edges','top_edge_jaccard','common_edge_spearman','coverage_state']
        with (output / 'edge_stability.tsv').open('w') as handle:
            writer = csv.DictWriter(handle, fieldnames=edge_fields, delimiter='\t'); writer.writeheader()
            keys = connection.execute('SELECT outer_fold,gene_id,method FROM linksa UNION SELECT outer_fold,gene_id,method FROM linksb ORDER BY 1,2,3')
            for key in keys:
                scores = [dict(connection.execute(f'SELECT peak_id,score FROM links{tag} WHERE outer_fold=? AND gene_id=? AND method=?', key)) for tag in ['a','b']]
                la, lb = map(set, scores)
                common = sorted(la & lb)
                if key[2] not in refitted_methods:
                    writer.writerow(dict(zip(edge_fields[:3], key), baseline_edges=len(la),
                        sensitivity_edges=len(lb), common_edges=np.nan, top_edge_jaccard=np.nan,
                        common_edge_spearman=np.nan, coverage_state='method_not_refit'))
                    continue
                correlation = pd.Series([scores[0][p] for p in common]).rank().corr(pd.Series([scores[1][p] for p in common]).rank()) if len(common) > 1 else np.nan
                writer.writerow(dict(zip(edge_fields[:3], key), baseline_edges=len(la), sensitivity_edges=len(lb), common_edges=len(common), top_edge_jaccard=len(common)/len(la|lb) if la|lb else np.nan, common_edge_spearman=correlation, coverage_state='both' if la and lb else 'baseline_only' if la else 'sensitivity_only'))
        join_keys = ['outer_fold','gene_id','endpoint','method','donor_id']
        on = ' AND '.join(f'a.{k}=b.{k}' for k in join_keys)
        selected = ','.join('a.'+k for k in join_keys)
        query = f"SELECT {selected},a.squared_error,b.squared_error,CASE WHEN b.donor_id IS NULL THEN 'left_only' ELSE 'both' END FROM predictionsa a LEFT JOIN predictionsb b ON {on} UNION ALL SELECT {','.join('b.'+k for k in join_keys)},NULL,b.squared_error,'right_only' FROM predictionsb b LEFT JOIN predictionsa a ON {on} WHERE a.donor_id IS NULL"
        totals, coverage = {}, {}
        fields = join_keys + ['squared_error_baseline','squared_error_sensitivity','_merge','baseline_minus_sensitivity_squared_error']
        with (output / 'prediction_matched_coverage.tsv').open('w') as handle:
            writer = csv.writer(handle, delimiter='\t'); writer.writerow(fields)
            for row in connection.execute(query):
                fold, gene, endpoint, method, donor, left, right, state = row
                difference = left-right if state == 'both' else None
                writer.writerow([*row, difference])
                key = (endpoint, method, state)
                coverage[key] = coverage.get(key, 0) + 1
                if state == 'both':
                    key = (endpoint, method, donor)
                    total, count = totals.get(key, (0., 0))
                    totals[key] = (total+difference, count+1)
        connection.close()
    donor_rows = [dict(endpoint=k[0], method=k[1], donor_id=k[2], n_common_gene_fold_pairs=n, baseline_minus_sensitivity_squared_error=total/n, estimand='mean_error_difference_over_identical_evaluable_gene_fold_pairs_within_donor') for k,(total,n) in sorted(totals.items())]
    donor = pd.DataFrame(donor_rows, columns=['endpoint','method','donor_id','n_common_gene_fold_pairs','baseline_minus_sensitivity_squared_error','estimand'])
    donor.to_csv(output / 'donor_prediction_difference.tsv', sep='\t', index=False)
    pd.DataFrame([dict(endpoint=k[0],method=k[1],_merge=k[2],gene_donor_fold_rows=n) for k,n in sorted(coverage.items())]).to_csv(output/'prediction_coverage.tsv',sep='\t',index=False)
    bootstrap = []
    for (endpoint, method), group in donor.groupby(['endpoint','method'], sort=True):
        values = group.baseline_minus_sensitivity_squared_error.to_numpy()
        offset = int(hashlib.sha256(str((endpoint,method)).encode()).hexdigest()[:8],16)
        rng = np.random.default_rng((seed+offset) % 2**32)
        means = np.array([rng.choice(values,size=len(values),replace=True).mean() for _ in range(bootstrap_reps)])
        low, high = np.quantile(means,[.025,.975]) if len(values)>1 else (np.nan,np.nan)
        bootstrap.append(dict(endpoint=endpoint,method=method,n_donors=len(values),mean_baseline_minus_sensitivity_squared_error=values.mean(),ci_low=low,ci_high=high,bootstrap_reps=bootstrap_reps,seed=seed,inference='exploratory_equal_donor_bootstrap_conditional_on_fixed_predictions_no_primary_Holm_calls'))
    pd.DataFrame(bootstrap).to_csv(output/'donor_bootstrap_difference.tsv',sep='\t',index=False)
    (output/'memory_and_estimand.json').write_text(json.dumps({'join':'disk_backed_SQLite','input_chunk_rows':rows_per_chunk,'SQLite_cache_MiB':32,'live_memory':'one input chunk, one gene-fold-method edge set, donor sum/count accumulators; independent of total prediction rows','bootstrap_unit':'donor','conditioning':'fixed trained predictions and matched gene coverage; no model refit in bootstrap','inference':'exploratory; no new primary multiplicity calls'},indent=2)+'\n')


def refit(module, counts, cells, peaks, genes, rna, rna_genes, settings, seed, baseline, output):
    aggregate, meta, groups, eligible = module.build_metacells(counts, cells, seed, int(settings.get('metacell_size', 50)), int(settings.get('max_cells_per_donor', 1000)))
    mapping = fold_map(baseline / 'fold_assignments.tsv')
    if meta.empty or not set(meta.donor_id).issubset(mapping):
        raise ValueError('Eligible sensitivity donors must have baseline outer folds')
    if len(set(mapping[d] for d in meta.donor_id)) < 2:
        raise ValueError('Fewer than two surviving outer folds')
    paired = None if rna is None else sparse.vstack([sparse.csr_matrix(rna[g].sum(axis=0)) for g in groups], format='csr')
    fit_settings = dict(settings, outer_fold_map=mapping)
    module.fit_group(aggregate, meta, peaks, genes, output, fit_settings, seed, paired, rna_genes)
    eligible.to_csv(output / 'eligibility.tsv', sep='\t', index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['config', 'input-root', 'genes', 'baseline', 'output']:
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--rna-root', type=Path)
    parser.add_argument('--stages', nargs='+', choices=['seeds', 'caps', 'omission', 'shuffle'], default=['seeds', 'caps', 'omission', 'shuffle'])
    parser.add_argument('--seeds', nargs='+', type=int, default=[43, 44])
    parser.add_argument('--caps', nargs='+', type=int, default=[500, 250])
    parser.add_argument('--omit-donors', nargs='+', help='Default: every donor; changes training and evaluation donor population')
    parser.add_argument('--execute', action='store_true')
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(args.output)
    completion = args.baseline / 'completion.json'
    if not completion.exists() or json.loads(completion.read_text()).get('status') != 'completed' or any(args.baseline.glob('INVALIDATED*.json')):
        raise ValueError('Sensitivity requires a completed baseline run')
    config = json.loads(args.config.read_text())
    seed = int(config.get('seed', 42))
    settings = dict(config.get('network', {}))
    baseline_settings = args.baseline / 'run_settings.json'
    if baseline_settings.exists():
        recorded = json.loads(baseline_settings.read_text())
        settings.update({k: v for k, v in recorded.items() if k != 'seed'})
        seed = int(recorded.get('seed', seed))
    # Native stage links were learned from a particular donor roster and seed.
    # Sensitivity refits need independently fitted native maps, unavailable here.
    excluded_external_methods = sorted(settings.pop('external_stage_links', {}))
    settings.pop('external_stage_roots', None)
    settings.pop('external_method_controls', None)
    cells = pd.read_csv(args.input_root / 'cells.tsv', sep='\t', dtype=str)
    baseline_donors = set()
    for folds in args.baseline.glob('*/fold_assignments.tsv'):
        baseline_donors.update(fold_map(folds))
    if not baseline_donors:
        raise ValueError('Baseline lacks qualified donor fold assignments')
    cells = cells[cells.donor_id.isin(baseline_donors)]
    scenarios = []
    if 'seeds' in args.stages:
        scenarios += [dict(name=f'seed_{s}', seed=s) for s in args.seeds if s != seed]
    if 'caps' in args.stages:
        if any(c < 100 or c % 50 or c >= int(settings.get('max_cells_per_donor', 1000)) for c in args.caps):
            raise ValueError('Caps must be lower whole-cell multiples of 50, at least 100')
        scenarios += [dict(name=f'cap_{c}', cap=c) for c in args.caps]
    if 'omission' in args.stages:
        donors = args.omit_donors or sorted(cells.donor_id.unique())
        if not set(donors).issubset(set(cells.donor_id)):
            raise ValueError('Unknown omitted donor')
        pairs = cells[cells.donor_id.isin(donors)][['cohort', 'donor_id']].drop_duplicates().sort_values(['cohort', 'donor_id'])
        scenarios += [dict(name='omit_' + hashlib.sha256((cohort + ':' + donor).encode()).hexdigest()[:12], omit_donor=donor, omit_cohort=cohort) for cohort, donor in pairs.itertuples(index=False, name=None)]
    if 'shuffle' in args.stages:
        if not args.rna_root:
            raise ValueError('Pairing shuffle requires genuine paired RNA')
        scenarios += [dict(name='pairing_shuffle', shuffle=True)]
    args.output.mkdir(parents=True, exist_ok=False)
    plan = dict(scenarios=scenarios, baseline=str(args.baseline.resolve()), inference_unit='donor', outer_folds='fixed_to_baseline_donor_assignments', inner_folds='refit_and_seeded_within_outer_training_donors', estimand='prediction_and_link_sensitivity_conditional_on_fixed_peak_and_gene_universe', omission_estimand='refit_on_remaining_donors_with_omitted_donor_absent_from_both_training_and_evaluation', shuffle_estimand='within_donor_nucleus_pairing_value_preserving_donor_RNA_distributions', native_LINGER_TRIPOD_RNA_only_prior_only='not_run_by_this_orchestrator; dedicated_native_adapters_and_qualified_inputs_required', executed=args.execute)
    plan['analysis_scope'] = settings.get('analysis_scope', 'candidate_sensitivity_no_adoption')
    plan['excluded_external_methods'] = excluded_external_methods
    plan['pilot'] = bool(settings.get('pilot', False) or settings.get('source_pilot', False))
    plan['source_sha256'] = {str(path.resolve()): hashlib.sha256(path.read_bytes()).hexdigest() for path in [Path(__file__), args.config, args.genes, Path(__file__).with_name('network_compare.py')]}
    (args.output / 'sensitivity_plan.json').write_text(json.dumps(plan, indent=2) + '\n')
    if not args.execute:
        return
    module = load_network()
    module.require_compute()
    counts, loaded_cells, peaks = module.load_input(args.input_root, seed, int(settings.get('max_cells_per_donor', 1000)))
    # Sparse chunk inputs may preselect cells; changing seed then measures aggregation/model stability conditional on those cells.
    plan['input_cell_universe'] = 'loaded_once_with_primary_seed; sparse_chunk_preselection_is_fixed'
    cells = loaded_cells
    rna, rna_genes = None, None
    if args.rna_root:
        if set(cells.cohort) != {'GSE296875'}:
            raise ValueError('Only GSE296875 has locally qualified genuine paired nuclei')
        rcells = module.read_table(args.rna_root / 'cells.tsv')
        if not module.cell_keys(cells).equals(module.cell_keys(rcells)):
            raise ValueError('RNA and ATAC exact cell order must agree')
        rna = module.read_matrix(args.rna_root, ['rna_counts.npz', 'RNA_counts.mtx.gz'])
        rna_genes = module.read_table(args.rna_root / 'genes.tsv')
        if rna.shape != (len(cells), len(rna_genes)) or rna_genes.gene_id.duplicated().any():
            raise ValueError('Invalid RNA dimensions or duplicate gene identities')
    genes = module.read_table(args.genes)
    failures = []
    for scenario in scenarios:
        root = args.output / scenario['name']
        root.mkdir()
        local_settings = dict(settings)
        (root / 'run_settings.json').write_text(json.dumps({'scenario': scenario, 'baseline_settings': settings, 'fixed_outer_folds': True, 'scientific_adoption': False}, indent=2) + '\n')
        local_settings['max_cells_per_donor'] = scenario.get('cap', settings.get('max_cells_per_donor', 1000))
        for (cohort, lineage), idx in cells.groupby(['cohort', 'lineage'], sort=True).indices.items():
            name = f'{cohort}_{lineage}'
            baseline = args.baseline / name
            if not (baseline / 'predictions.tsv').exists():
                continue
            qualified_donors = fold_map(baseline / 'fold_assignments.tsv')
            idx = np.array([i for i in idx if cells.iloc[i].donor_id in qualified_donors and not (cohort == scenario.get('omit_cohort') and cells.iloc[i].donor_id == scenario.get('omit_donor'))])
            local = cells.iloc[idx].reset_index(drop=True)
            local_rna = None if rna is None else rna[idx]
            if scenario.get('shuffle'):
                permutation = shuffle_index(local, seed)
                pd.DataFrame({'ATAC_cell_id': local.cell_id, 'RNA_cell_id': local.iloc[permutation].cell_id.to_numpy(), 'donor_id': local.donor_id}).to_csv(root / f'{name}_pairing.tsv', sep='\t', index=False)
                local_rna = local_rna[permutation]
            x = counts[idx].copy()
            if settings.get('representation', 'binary') == 'binary':
                x.eliminate_zeros()
                x.data[:] = 1
            target = root / name
            try:
                refit(module, x, local, peaks, genes, local_rna, rna_genes, local_settings, scenario.get('seed', seed), baseline, target)
                compare_tables(baseline, target, root / f'{name}_stability', seed=seed)
            except Exception as exc:
                failures.append({'scenario': scenario['name'], 'context': name, 'reason': str(exc)})
                (root / f'{name}_failure.json').write_text(json.dumps(failures[-1], indent=2) + '\n')
        (root / 'completion.json').write_text(json.dumps({'status': 'partial' if any(f['scenario'] == scenario['name'] for f in failures) else 'completed', 'scientific_adoption': False}) + '\n')
    (args.output / 'sensitivity_plan.json').write_text(json.dumps(plan, indent=2) + '\n')
    (args.output / 'failures.json').write_text(json.dumps(failures, indent=2) + '\n')
    (args.output / 'completion.json').write_text(json.dumps({'status': 'partial' if failures else 'completed', 'scientific_adoption': False}) + '\n')


if __name__ == '__main__':
    main()
