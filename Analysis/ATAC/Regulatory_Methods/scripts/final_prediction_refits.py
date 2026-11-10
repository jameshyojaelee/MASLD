#!/usr/bin/env python3
"""Fixed-fold, wider-grid hepatocyte refits for the final MASLD candidate.

Omission affects training only. All scenario folds inherit donor assignments
from the original seed-42 analysis, including inner folds. No outcome is used
to choose genes, cells, folds, or regularization ranges.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.util
import json
import multiprocessing as mp
import os
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd
from scipy import sparse
from threadpoolctl import threadpool_limits

from promoter_lookup import PromoterLookup
from ridge_path import tune_path

A = Path(__file__).resolve().parents[1]
WIDE = A / 'results/20261024/ridge_grid_sensitivity_r1'
PREPARED = WIDE / 'hepatocyte/prepared'
GRID = [0.1, 1., 10., 100., 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]
CORE = None
CTX = None
LIMIT = None


def dump(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def load_core():
    global CORE
    source = PREPARED / 'source_snapshot/network_compare.py'
    spec = importlib.util.spec_from_file_location('final_frozen_core', source)
    CORE = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(CORE)
    CORE.require_compute()
    receipt = json.loads((PREPARED / 'preparation.json').read_text())
    expected = [v for k, v in receipt['source_inventory'].items() if Path(k).name == source.name]
    assert expected == [hashlib.sha256(source.read_bytes()).hexdigest()]
    return receipt


def permutation(cells, seed):
    result = np.arange(len(cells))
    for key, rows in cells.groupby(['cohort', 'lineage', 'donor_id'], sort=True).indices.items():
        offset = int(hashlib.sha256(str(key).encode()).hexdigest()[:8], 16)
        result[rows] = np.random.default_rng((seed + offset) % 2**32).permutation(rows)
    return result


def retained_permutation(cells, groups, seed):
    """Shuffle only evaluated nuclei; preserve every donor's aggregate RNA sum."""
    retained = np.concatenate(groups)
    assert len(set(retained)) == len(retained)
    result = np.arange(len(cells))
    result[retained] = retained[permutation(cells.iloc[retained].reset_index(drop=True), seed)]
    return result


def prepare(root):
    receipt = load_core()
    root.mkdir(exist_ok=False)
    cells_root = Path(receipt['input_root'])
    counts, cells, peaks = CORE.load_input(cells_root, 42, 1000)
    rna_cells = CORE.read_table(cells_root / 'cells.tsv')
    assert CORE.cell_keys(cells).equals(CORE.cell_keys(rna_cells))
    rna = CORE.read_matrix(cells_root, ['rna_counts.npz', 'RNA_counts.mtx.gz'])
    assert rna.shape[0] == len(cells)
    counts = counts.copy(); counts.eliminate_zeros(); counts.data[:] = 1
    folds = CORE.read_table(PREPARED / 'fold_assignments.tsv')
    outer_map = folds.groupby('donor_id').outer_fold.first().astype(int).to_dict()
    assert folds.groupby('donor_id').outer_fold.nunique().eq(1).all()
    old_meta = CORE.read_table(PREPARED / 'metadata.tsv')
    old_donors = old_meta.donor_id.to_numpy()
    inner_map = {}
    for fold in range(5):
        idx = np.flatnonzero(folds.outer_fold.to_numpy(int) != fold)
        inner = CORE.donor_folds(old_donors[idx], 3, CORE.stable_seed(42, ('inner', fold)))
        inner_map[str(fold)] = dict(zip(old_donors[idx], map(int, inner)))
    scores = pd.read_csv(A / 'results/20261024/wide_proximity_comparison_r1/donor_scores.tsv', sep='\t', dtype={'donor_id': str})
    scores = scores[(scores.lineage == 'hepatocyte') & (scores.method == 'cicero')].copy()
    scores['absolute_influence'] = np.abs(scores.gain - scores.gain.mean()) / (len(scores) - 1)
    selected = scores.sort_values(['absolute_influence', 'donor_id'], ascending=[False, True]).head(3)
    selected.to_csv(root / 'targeted_donor_selection.tsv', sep='\t', index=False)
    scenarios = [dict(name='baseline', aggregation_seed=42, cap=1000),
                 dict(name='aggregation_43', aggregation_seed=43, cap=1000),
                 dict(name='aggregation_44', aggregation_seed=44, cap=1000),
                 dict(name='cap_500', aggregation_seed=42, cap=500)]
    scenarios += [dict(name=f'pairing_shuffle_{s}', aggregation_seed=42, cap=1000, shuffle_seed=s) for s in [43,44,45]]
    scenarios += [dict(name=f'omit_training_{d}', aggregation_seed=42, cap=1000, omit_donor=d) for d in selected.donor_id]
    for scenario in scenarios:
        target = root / scenario['name']; target.mkdir()
        aggregate, meta, groups, eligibility = CORE.build_metacells(counts, cells, scenario['aggregation_seed'], 50, scenario['cap'])
        assert set(meta.donor_id) == set(outer_map), 'A robustness scenario changed the donor roster'
        perm = retained_permutation(cells, groups, scenario['shuffle_seed']) if 'shuffle_seed' in scenario else np.arange(len(cells))
        paired = sparse.vstack([sparse.csr_matrix(rna[perm[group]].sum(axis=0)) for group in groups], format='csr')
        assert all(len(set(cells.iloc[group].donor_id)) == 1 for group in groups)
        assert all(len(set(cells.iloc[perm[group]].donor_id)) == 1 for group in groups)
        if scenario['name'] == 'baseline':
            original = sparse.load_npz(PREPARED / 'counts.npz')
            original_rna = sparse.load_npz(PREPARED / 'rna.npz')
            assert aggregate.shape == original.shape and (aggregate != original).nnz == 0
            assert paired.shape == original_rna.shape and (paired != original_rna).nnz == 0
            assert CORE.table_digest(meta) == CORE.table_digest(old_meta)
        sparse.save_npz(target / 'counts.npz', aggregate)
        sparse.save_npz(target / 'rna.npz', paired)
        meta.to_csv(target / 'metadata.tsv', sep='\t', index=False)
        eligibility.to_csv(target / 'eligibility.tsv', sep='\t', index=False)
        memberships = []
        for i, group in enumerate(groups):
            for n in group:
                memberships.append(dict(metacell_id=meta.iloc[i].metacell_id, donor_id=cells.iloc[n].donor_id,
                                        ATAC_cell_id=cells.iloc[n].cell_id, RNA_cell_id=cells.iloc[perm[n]].cell_id))
        pd.DataFrame(memberships).to_csv(target / 'membership.tsv.gz', sep='\t', index=False)
        scenario.update(n_donors=meta.donor_id.nunique(), n_metacells=len(meta),
                        counts_sha256=CORE.matrix_digest(aggregate), rna_sha256=CORE.matrix_digest(paired),
                        metadata_sha256=CORE.table_digest(meta))
        dump(target / 'scenario.json', scenario)
        print(json.dumps(scenario), flush=True)
    dump(root / 'specification.json', dict(status='completed', scenarios=scenarios, outer_map=outer_map,
        inner_map=inner_map, original_prepared=str(PREPARED), ridge_alphas=GRID,
        inference='conditional_donor_bootstrap_descriptive_no_new_significance_family',
        omission='remove_from_training_only_keep_original_test_donors', scientific_adoption=False))


def repair_shuffles(root, output):
    """New-directory correction: conserve RNA marginals in retained nuclei too."""
    load_core()
    output.mkdir(exist_ok=False)
    spec = json.loads((root / 'specification.json').read_text())
    receipt = json.loads((PREPARED / 'preparation.json').read_text())
    source = Path(receipt['input_root'])
    cells = CORE.read_table(source / 'cells.tsv')
    rna = CORE.read_matrix(source, ['rna_counts.npz', 'RNA_counts.mtx.gz'])
    positions = dict(zip(cells.cell_id, range(len(cells))))
    membership = pd.read_csv(root / 'baseline/membership.tsv.gz', sep='\t', dtype=str)
    metadata = CORE.read_table(root / 'baseline/metadata.tsv')
    grouped = {k:g.ATAC_cell_id.map(positions).to_numpy(int) for k,g in membership.groupby('metacell_id',sort=False)}
    groups = [grouped[k] for k in metadata.metacell_id]
    baseline = sparse.load_npz(root / 'baseline/rna.npz')
    for scenario in spec['scenarios']:
        old = root / scenario['name']; new = output / scenario['name']
        if 'shuffle_seed' not in scenario:
            new.symlink_to(old.resolve(), target_is_directory=True)
            continue
        new.mkdir()
        for name in ['counts.npz','metadata.tsv','eligibility.tsv']:
            (new/name).symlink_to((old/name).resolve())
        perm = retained_permutation(cells, groups, scenario['shuffle_seed'])
        paired = sparse.vstack([sparse.csr_matrix(rna[perm[g]].sum(axis=0)) for g in groups],format='csr')
        for donor in metadata.donor_id.unique():
            idx = np.flatnonzero(metadata.donor_id.eq(donor))
            assert np.array_equal(np.asarray(paired[idx].sum(axis=0)),np.asarray(baseline[idx].sum(axis=0)))
        sparse.save_npz(new/'rna.npz',paired)
        updated=[]
        for i,g in enumerate(groups):
            for n in g:
                updated.append(dict(metacell_id=metadata.iloc[i].metacell_id,donor_id=cells.iloc[n].donor_id,
                    ATAC_cell_id=cells.iloc[n].cell_id,RNA_cell_id=cells.iloc[perm[n]].cell_id))
        pd.DataFrame(updated).to_csv(new/'membership.tsv.gz',sep='\t',index=False)
        scenario.update(rna_sha256=CORE.matrix_digest(paired),retained_donor_RNA_count_conservation=True)
        dump(new/'scenario.json',scenario)
    spec['predecessor']=str(root.resolve())
    spec['shuffle_definition']='permute_within_donor_lineage_retained_nuclei_only; donor_RNA_UMIs_exactly_conserved'
    dump(output/'specification.json',spec)
    print('PASS: retained donor RNA totals conserved for all three shuffles',flush=True)


def init_worker():
    global LIMIT
    LIMIT = threadpool_limits(limits=1)


def fixed_donor_splits(donors, outer_map, inner_map, fold, omitted=''):
    outer = np.array([outer_map[d] for d in donors], int)
    train = np.flatnonzero((outer != fold) & (donors != omitted))
    test = np.flatnonzero(outer == fold)
    inner = np.array([inner_map[str(fold)][d] for d in donors[train]], int)
    assert len(set(inner)) == 3 and not set(donors[train]) & set(donors[test])
    assert omitted not in set(donors[train])
    return outer, train, test, inner


def fit_gene(job):
    number, gene = job
    c = CTX; gid = gene['gene_id']
    rows, links, tuning = [], [], []
    if gid not in c['rna_lookup']:
        return number, gid, rows, links, tuning
    prom, candidates = c['geometry'].query(gene, 250000, 2000)
    if not len(prom) or not len(candidates):
        return number, gid, rows, links, tuning
    promoter = c['promoter'].promoter(prom)
    train, test, inner, donors = [c[k] for k in ['train','test','inner','donors']]
    if np.var(promoter[train]) <= 1e-12:
        return number, gid, rows, links, tuning
    y = CORE.log_columns(c['rna'], [c['rna_lookup'][gid]]).ravel()
    weights = CORE.donor_weights(donors[train])
    variance = float(np.average((y[train] - np.average(y[train], weights=weights))**2, weights=weights))
    if variance <= 1e-12:
        return number, gid, rows, links, tuning
    for method in ['promoter_only', 'distance', 'cicero']:
        if method == 'promoter_only':
            design = promoter[:, None]
            designs = {k: design for k in sorted(set(inner))}
        else:
            chosen, scores = CORE.selected_indices(method, candidates, prom, c['peaks'], gene, c['training_x'],
                promoter[train], donors[train], CORE.stable_seed(42, (c['fold'], gid)), 10, c['cicero']['outer'])
            if not len(chosen):
                continue
            design = CORE.design_matrix(c['feature_x'], chosen, promoter, 'rna')
            links += [dict(outer_fold=c['fold'], gene_id=gid, method=method, rank=rank,
                           peak_id=c['peaks'].iloc[i].peak_id, score=float(score))
                      for rank, (i, score) in enumerate(zip(chosen, scores), 1)]
            designs = {}
            for k in sorted(set(inner)):
                it = train[inner != k]
                selected, _ = CORE.selected_indices(method, candidates, prom, c['peaks'], gene, c['inner_x'][k],
                    promoter[it], donors[it], CORE.stable_seed(42, (c['fold'], k, gid)), 10, c['cicero'][f'inner{k}'])
                designs[k] = CORE.design_matrix(c['feature_x'], selected, promoter, 'rna')
        losses = tune_path(CORE, designs, y, donors, train, inner, GRID)
        alpha = min(GRID, key=lambda a: (losses[a], a))
        prediction = CORE.weighted_ridge(design[train], y[train], donors[train], alpha, design[test])
        CORE.save_predictions(rows, c['fold'], gid, 'rna', method, test, c['metadata'], y, prediction, variance, alpha)
        tuning.append(dict(outer_fold=c['fold'], gene_id=gid, method=method, ridge_alpha=alpha,
                           inner_mse=losses[alpha], upper_boundary_mse=losses[GRID[-1]], n_features=design.shape[1]))
    return number, gid, rows, links, tuning


def run_fold(root, scenario_name, fold, output, workers, max_genes):
    global CTX
    receipt = load_core()
    started = time.monotonic()
    spec = json.loads((root / 'specification.json').read_text())
    assert spec['ridge_alphas'] == GRID
    scenario = json.loads((root / scenario_name / 'scenario.json').read_text())
    output.mkdir(exist_ok=False)
    source = root / scenario_name
    aggregate = sparse.load_npz(source / 'counts.npz')
    rna = sparse.load_npz(source / 'rna.npz')
    meta = CORE.read_table(source / 'metadata.tsv')
    assert CORE.matrix_digest(aggregate) == scenario['counts_sha256']
    assert CORE.matrix_digest(rna) == scenario['rna_sha256']
    donors = meta.donor_id.to_numpy()
    outer, train, test, inner = fixed_donor_splits(donors, spec['outer_map'], spec['inner_map'], fold,
                                                 scenario.get('omit_donor', ''))
    pd.DataFrame(dict(metacell_id=meta.metacell_id, donor_id=donors, outer_fold=outer)).to_csv(output / 'outer_folds.tsv', sep='\t', index=False)
    pd.DataFrame(dict(metacell_id=meta.iloc[train].metacell_id, donor_id=donors[train], inner_fold=inner)).to_csv(output / 'inner_folds.tsv', sep='\t', index=False)
    peaks = CORE.read_table(PREPARED / 'peaks.tsv')
    peaks[['start0','end']] = peaks[['start0','end']].astype(np.int64)
    genes = CORE.read_table(PREPARED / 'genes.tsv'); genes.tss0 = genes.tss0.astype(np.int64)
    rgenes = CORE.read_table(PREPARED / 'rna_genes.tsv')
    # All ATAC-only map reuse requires exact aggregate counts, peak and donor order.
    sys.path.insert(0, str(A / 'scripts'))
    import ridge_grid_sensitivity as wider
    wider.CORE = CORE
    cache = {}; cache_checks = []
    stages = [('outer', train)] + [(f'inner{k}', train[inner != k]) for k in sorted(set(inner))]
    base_count_sha = next(s['counts_sha256'] for s in spec['scenarios'] if s['name'] == 'baseline')
    for tag, indices in stages:
        old = WIDE / f'hepatocyte/fold{fold}/{tag}/source_verification.json'
        can_reuse = scenario['counts_sha256'] == base_count_sha and not scenario.get('omit_donor')
        if can_reuse:
            old_root = Path(json.loads(old.read_text())['source'])
            cache[tag], checked = wider.cached_cicero(old_root, output / tag, aggregate, peaks, meta, indices, 100)
            checked['reason'] = 'ATAC_only_method_exact_counts_and_training_donor_order_unchanged'
            cache_checks.append(checked)
        else:
            nonzero = np.flatnonzero(np.asarray(aggregate[indices].sum(axis=0)).ravel() > 0)
            print(f'fit Cicero {scenario_name} fold={fold} stage={tag} donors={len(set(donors[indices]))}', flush=True)
            cache[tag] = CORE.cicero_scores(aggregate[indices][:,nonzero], peaks.iloc[nonzero], output / tag,
                receipt['settings']['cicero_rscript'], CORE.stable_seed(42, (fold, tag)), 500000,
                distance_samples=100, donor_ids=donors[indices], workers=workers)
            log = (output / tag / 'run.log').read_text()
            import re
            assert re.findall(r'Models with errors:\s*(\d+)', log) == ['0']
    x = CORE.normalize_sparse(aggregate)
    CTX = dict(peaks=peaks, geometry=CORE.PeakGeometry(peaks), metadata=meta, donors=donors,
               train=train, test=test, inner=inner, fold=fold, promoter=PromoterLookup(x), feature_x=x.tocsc(),
               training_x=x[train].tocsc(), inner_x={k:x[train[inner != k]].tocsc() for k in sorted(set(inner))},
               rna=CORE.normalize_sparse(rna).tocsc(), rna_lookup={g:i for i,g in enumerate(rgenes.gene_id)}, cicero=cache)
    jobs = list(enumerate(genes.sort_values('gene_id').to_dict('records'), 1))
    if max_genes:
        # A deterministic technical subset with existing baseline predictions.
        pilot_ids = set()
        for chunk in pd.read_csv(WIDE / f'hepatocyte/fold{fold}/predictions.tsv', sep='\t', chunksize=100000):
            pilot_ids.update(chunk.gene_id.unique())
            if len(pilot_ids) >= max_genes:
                break
        jobs = [j for j in jobs if j[1]['gene_id'] in pilot_ids][:max_genes]
    columns = ['outer_fold','gene_id','endpoint','method','donor_id','n_metacells','observed','predicted','squared_error','standardized_squared_error','ridge_alpha']
    predictions = CORE.TableStream(output / 'predictions.tsv', columns)
    links = CORE.TableStream(output / 'links.tsv', ['outer_fold','gene_id','method','rank','peak_id','score'])
    tuning = CORE.TableStream(output / 'tuning.tsv', ['outer_fold','gene_id','method','ridge_alpha','inner_mse','upper_boundary_mse','n_features'])
    pool = mp.get_context('fork').Pool(workers, initializer=init_worker)
    try:
        for number, gid, rows, edges, tune in pool.imap(fit_gene, jobs, chunksize=4):
            predictions.append(rows); links.append(edges); tuning.append(tune)
            if number % 200 == 0 or max_genes:
                print(f'readout {scenario_name} fold={fold} gene={number} rows={predictions.rows}', flush=True)
        pool.close(); pool.join()
    except BaseException:
        pool.terminate(); pool.join(); raise
    finally:
        predictions.close(); links.close(); tuning.close()
    if scenario_name == 'baseline':
        keys = ['outer_fold','gene_id','endpoint','method','donor_id']
        actual = pd.read_csv(output / 'predictions.tsv', sep='\t', dtype={'donor_id':str})
        expected = []
        gids = set(actual.gene_id)
        for chunk in pd.read_csv(WIDE / f'hepatocyte/fold{fold}/predictions.tsv', sep='\t', chunksize=100000, dtype={'donor_id':str}):
            expected.append(chunk[chunk.grid.eq('wide') & chunk.method.isin(['promoter_only','distance','cicero']) & chunk.gene_id.isin(gids)])
        expected = pd.concat(expected).sort_values(keys).reset_index(drop=True)
        actual = actual.sort_values(keys).reset_index(drop=True)
        assert actual[keys].equals(expected[keys])
        numeric = ['observed','predicted','squared_error','standardized_squared_error','ridge_alpha']
        assert np.allclose(actual[numeric], expected[numeric], rtol=1e-10, atol=1e-11)
        dump(output / 'baseline_reproduction.json', dict(status='passed', rows=len(actual),
             max_absolute_difference=float(np.max(np.abs(actual[numeric].to_numpy()-expected[numeric].to_numpy())))))
    dump(output / 'completion.json', dict(status='completed', scenario=scenario, outer_fold=fold,
        n_held_donors=len(set(donors[test])), n_training_donors=len(set(donors[train])),
        elapsed_seconds=time.monotonic()-started, predictions=predictions.rows, technical_pilot=bool(max_genes),
        cache_checks=cache_checks, ridge_alphas=GRID, scientific_adoption=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage', choices=['prepare','repair-shuffles','fold'])
    parser.add_argument('--root', required=True, type=Path)
    parser.add_argument('--scenario')
    parser.add_argument('--fold', type=int, choices=range(5))
    parser.add_argument('--output', type=Path)
    parser.add_argument('--workers', type=int, default=4)
    parser.add_argument('--max-genes', type=int, default=0)
    args = parser.parse_args()
    if args.stage == 'prepare':
        prepare(args.root)
    elif args.stage == 'repair-shuffles':
        repair_shuffles(args.root,args.output)
    else:
        assert args.scenario and args.fold is not None and args.output
        run_fold(args.root, args.scenario, args.fold, args.output, args.workers, args.max_genes)


if __name__ == '__main__':
    main()
