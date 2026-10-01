#!/usr/bin/env python3
"""Frozen well-held lineage RNA to common-fragment ATAC development producer.

Preparation does not authorize execution. Main requires an independent pinned
execution decision. LineageScorer is an RNA/nucleus/lineage-only callable.
"""
import argparse
import contextlib
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time

import numpy as np

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ALIGNED = REC / 'codex_gse296875_aligned_native_pseudobulk_21999585'
COMMON = REC / 'codex_gse296875_common_fragment_counts_22000534'
GLOBAL = REC / 'codex_gse296875_donor_global_raw_22000941'
RECIPE = REC / '20261001_codex_gse296875_common_fragment_rna_atac_prespec.json'
RECIPE_SHA = '17ab84770eb3bfac46ea50d0a1400e52a70ffea1bbb7683e45b00a40e499403c'
GENE_HASH = 'afbcc02143e68d1fbcda9631917f1344adfd3ac8efb752e2ebb7f044df934e60'
REGION_HASH = '50a52d097a5232edcf136dd53d6549ca2d1aa5dce9a57cf055ff10b7b2e0f3e3'
G, R, U, K = 36601, 306706, 273, 20
ALPHA = 10.0
CAP = 64 * 1024**3
RESERVE = 1024**2
ARMS = ('lineage_mean', 'technical_only', 'donor_global_RNA', 'target_lineage_RNA')
PARTITIONS = ('full', 'halfA', 'halfB')
MODES = (('full', 'full', 'full'), ('halfA_RNA_to_halfB_ATAC', 'halfA', 'halfB'),
         ('halfB_RNA_to_halfA_ATAC', 'halfB', 'halfA'))
GUARDS = {
    RECIPE: RECIPE_SHA,
    ROOT / 'Analysis/SingleCell/continual_integration/reference/gse296875_source_lock_v44.json': '6be686b4187309f038097cf73d3f2a8dd34e72d666bc92733a53e99a1ccbb620',
    ALIGNED / 'artifacts.json': '3c4102ecde36825237b7d31e6ff4c4395570dc97ee222e274429833cc000d8a5',
    ALIGNED / 'summary.json': 'dd1986b643e1e5ae2cec54d255fb8ed7619c6db668541194a4032277d8c00226',
    ALIGNED / 'unit_identity.tsv': '5c0a1475cb5cbf71ae4e25f31652a822d4b773c4cd8f5f962bcc44b5bb80be68',
    ALIGNED / 'genes.tsv': '2a33d4d339e2c6db695525ae2655af5156ccda008300bd992e306d48f01d38f0',
    ALIGNED / 'native_peak_geometry.tsv': REGION_HASH,
    ALIGNED / 'cell_identity.tsv': '438f972c0883837fc254260a82f058caf95d634643d2a2b13156ba4089b848d4',
    COMMON / 'artifacts.json': '8735af6f215fc82be928f2349cb8182d3d76e834304cce8eb32e858f5880a960',
    COMMON / 'summary.json': 'bd64347bc066cbd182695c466b957ddd3c6efc689b53f8db2d557aed34321f8e',
    COMMON / 'unit_identity.tsv': '509dadfb442445c28eb4139d2f3a9699d6734c82a3ac6442b29e2d767417a964',
    COMMON / 'units.tsv': 'babe4cd3bc0f4420c6303526d9ebbe342ae69ca792237a2c05b2bf09aa780f2a',
    GLOBAL / 'artifacts.json': '5629741635252c1869cc68da59c365833ad8912fa7f46d62fef25335beead111',
    GLOBAL / 'summary.json': '9684d83681b898ac507680e2df39a4781b42daccbd2abfac6342faa88031fb19',
    GLOBAL / 'donors.tsv': '5263e2ba9c65ca6bb8e47d12542b55df1e2faa0c67a1d71d5012f292e9328aad',
}


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def pin(path, expected, size=None):
    path = Path(path)
    require(path.is_file(), 'required source missing')
    if size is not None:
        require(path.stat().st_size == size, 'source size differs')
    observed = sha(path)
    require(observed == expected, 'source SHA differs')
    return {'path': str(path), 'sha256': observed, 'bytes': path.stat().st_size}


def load_json(path):
    require(Path(path).stat().st_size <= 1024**2, 'metadata JSON exceeds limit')
    return json.loads(Path(path).read_text())


def table(path):
    with Path(path).open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


class Output:
    def __init__(self, directory, logs):
        self.directory, self.logs = directory, logs
        self.reserved = sum(p.stat().st_size for p in directory.rglob('*') if p.is_file())
        self.start = time.monotonic()

    def reserve(self, size):
        require(self.reserved + size <= CAP - RESERVE, '64GiB saved-output cap exceeded')
        require(time.monotonic() - self.start < 88 * 3600, '88h producer time cap exceeded')
        self.reserved += size

    def path(self, name):
        path = self.directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        require(not path.exists(), 'output already exists')
        return path

    def json(self, name, value):
        body = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
        self.reserve(len(body))
        self.path(name).write_bytes(body)

    def tsv(self, name, rows):
        require(bool(rows), 'empty table must have explicit unavailable receipt')
        fields = list(rows[0])
        stream = io.StringIO()
        writer = csv.DictWriter(stream, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
        body = stream.getvalue().encode()
        self.reserve(len(body))
        self.path(name).write_bytes(body)

    def array(self, name, values):
        self.reserve(values.nbytes + 4096)
        path = self.path(name)
        with path.open('xb') as stream:
            np.save(stream, values, allow_pickle=False)
        return path

    def parameters(self, name, values):
        self.reserve(sum(v.nbytes for v in values.values()) + 65536)
        path = self.path(name)
        with path.open('xb') as stream:
            np.savez(stream, **values)
        return path

    def empty_oof(self, name):
        self.reserve(U * R * 8 + 4096)
        path = self.path(name)
        mapped = np.lib.format.open_memmap(path, mode='w+', dtype='<f8', shape=(U, R))
        for row in range(U):
            mapped[row] = np.nan  # explicit unmeasured storage sentinel, never biological zero
        mapped.flush()
        mapped._mmap.close()
        return path

    def check(self):
        actual = sum(p.stat().st_size for p in self.directory.rglob('*') if p.is_file())
        actual += sum(p.stat().st_size for p in self.logs if p.exists())
        require(actual <= CAP - RESERVE and actual <= self.reserved + RESERVE,
                'inclusive output accounting differs or exceeds cap')
        require(time.monotonic() - self.start < 88 * 3600, '88h producer time cap exceeded')
        return actual


def write_oof(path, indices, values):
    mapped = np.load(path, mmap_mode='r+', allow_pickle=False)
    mapped[indices] = values
    mapped.flush()
    mapped._mmap.close()


def normalize_counts(counts):
    require(counts.ndim == 2 and counts.shape[1] == G and counts.dtype.kind in 'iu',
            'complete ordered integer RNA axis required')
    require(np.all(counts >= 0), 'negative RNA counts')
    totals = np.array([sum(int(v) for v in row) for row in counts], dtype=object)
    require(all(0 < total < 2**53 for total in totals), 'RNA denominator undefined or inexact float64 range')
    totals = np.asarray(totals, dtype=np.float64)
    result = np.log1p(np.asarray(counts, dtype=np.float64) * (1e6 / totals[:, None]))
    require(np.isfinite(result).all(), 'RNA normalization nonfinite')
    return result, totals


def technical(counts, nuclei):
    _, totals = normalize_counts(counts)
    require(np.all(nuclei > 0), 'nucleus input undefined')
    return np.column_stack((np.log1p(totals), np.log1p(nuclei), np.count_nonzero(counts, axis=1) / G))


def standard_fit(values):
    mean = np.mean(values, axis=0)
    sd = np.std(values, axis=0, ddof=1)
    # A repeated float can have a rounded nonzero SD; compare actual rows for
    # exact constancy, without imposing a nearzero numerical threshold.
    active = np.any(values != values[0], axis=0)
    require(np.isfinite(mean).all() and np.isfinite(sd).all() and np.all(sd[active] > 0),
            'nonconstant training column has undefined sample SD')
    safe = np.where(active, sd, 1.0)
    standardized = (values - mean) / safe
    standardized[:, ~active] = 0.0
    return standardized, mean, safe, active


def pca_fit(values):
    mean = np.mean(values, axis=0)
    centered = values - mean
    _, singular, vt = np.linalg.svd(centered, full_matrices=False)
    # Only structural centered-matrix availability; no numerical-rank mask.
    available = min(K, values.shape[0] - 1, G)
    loadings = np.zeros((K, G), dtype=np.float64)
    used = available
    loadings[:used] = vt[:used]
    for row in range(used):
        maximal_gene = int(np.argmax(np.abs(loadings[row])))
        if loadings[row, maximal_gene] < 0:
            loadings[row] *= -1
    scores, score_mean, score_sd, active = standard_fit(centered @ loadings.T)
    active[used:] = False
    scores[:, ~active] = 0.0
    return scores, {'gene_mean': mean, 'loadings': loadings, 'score_mean': score_mean,
        'score_sd': score_sd, 'score_active': active, 'singular_values': singular,
        'structurally_available': np.array(available),
        'exact_tied_adjacent_singular_pairs': np.flatnonzero(singular[:-1] == singular[1:])}


def pca_apply(values, params):
    scores = ((values - params['gene_mean']) @ params['loadings'].T - params['score_mean']) / params['score_sd']
    scores[:, ~params['score_active']] = 0.0
    require(np.isfinite(scores).all(), 'frozen PCA projection nonfinite')
    return scores


def fit_ridge(features, targets):
    x_mean, y_mean = np.mean(features, axis=0), np.mean(targets, axis=0)
    if features.shape[1] == 0:
        return {'coefficients': np.zeros((0, R)), 'intercept': y_mean, 'feature_mean': x_mean}, 0.0
    x = features - x_mean
    gram = x.T @ x + ALPHA * np.eye(features.shape[1])
    rhs = x.T @ (targets - y_mean)
    coefficients = np.linalg.solve(gram, rhs)
    residual = float(np.max(np.abs(gram @ coefficients - rhs)))
    require(residual <= 5e-10 * (1.0 + float(np.max(np.abs(rhs)))), 'ridge normal equations failed')
    intercept = y_mean - x_mean @ coefficients
    require(np.isfinite(coefficients).all() and np.isfinite(intercept).all(), 'ridge parameters nonfinite')
    return {'coefficients': coefficients, 'intercept': intercept, 'feature_mean': x_mean}, residual


def predict(features, params):
    raw = features @ params['coefficients'] + params['intercept']
    require(np.isfinite(raw).all(), 'raw prediction nonfinite')
    return raw, np.maximum(raw, 0.0)


def pack(ridge, tech_state, pca_state=None):
    result = dict(ridge)
    if tech_state is not None:
        for name, value in zip(('technical_mean', 'technical_sd', 'technical_active'), tech_state):
            result[name] = value
    if pca_state is not None:
        result.update({'PCA_' + name: value for name, value in pca_state.items()})
    return result


def target_rows(atac, denominator, indices):
    d = denominator[indices]
    require(np.all(d > 0) and np.all(d < 2**53), 'ATAC denominator undefined or inexact float64 range')
    counts = np.asarray(atac[indices], dtype=np.float64)
    require(np.all(counts >= 0), 'negative ATAC counts')
    result = np.log1p(counts * (1e6 / d[:, None]))
    require(np.isfinite(result).all(), 'ATAC normalization nonfinite')
    return result


def identities(recipe):
    units = table(ALIGNED / 'unit_identity.tsv')
    measured = table(COMMON / 'units.tsv')
    counted_identity = table(COMMON / 'unit_identity.tsv')
    donors = table(GLOBAL / 'donors.tsv')
    require(len(units) == U and len(measured) == U and len(counted_identity) == U and len(donors) == 39,
            'source identity row counts differ')
    lineages = recipe['source_identity']['ordered_author_lineages']
    roster = recipe['source_identity']['well_donor_roster']
    donor_order = [d for ids in roster.values() for d in ids]
    require([d['donor_id'] for d in donors] == donor_order, 'global donor order differs')
    for index, donor in enumerate(donors):
        require(int(donor['donor_index']) == index and donor['donor_id'] in roster[donor['well_id']],
                'global donor identity differs')
    keys = ('unit_index', 'donor_id', 'well_id', 'native_source_celltype', 'nuclei', 'halfA_nuclei', 'halfB_nuclei')
    for index, unit in enumerate(units):
        require(int(unit['unit_index']) == index and unit['donor_id'] == donor_order[index // 7] and
                unit['native_source_celltype'] == lineages[index % 7], 'unit order differs')
        require(unit['donor_id'] in roster[unit['well_id']], 'well ownership differs')
        require(all(unit[k] == measured[index][k] == counted_identity[index][k] for k in keys),
                'RNA and common-fragment unit identities differ')
        require(int(unit['nuclei']) == int(unit['halfA_nuclei']) + int(unit['halfB_nuclei']), 'half membership differs')
    full = np.array([int(u['nuclei']) > 0 for u in units])
    halves = np.array([int(u['halfA_nuclei']) > 0 and int(u['halfB_nuclei']) > 0 for u in units])
    require(int(full.sum()) == 272 and int(halves.sum()) == 270, 'fixed sampled roster differs')
    require([(u['donor_id'], u['native_source_celltype']) for i, u in enumerate(units) if not full[i]] == [('733', 'B cells')],
            'unsampled identity differs')
    require([i for i in range(U) if full[i] and not halves[i]] == [91, 224], 'fixed half exclusions differ')
    gene_rows = table(ALIGNED / 'genes.tsv')
    require(len(gene_rows) == G and [int(g['gene_index']) for g in gene_rows] == list(range(G)), 'gene axis differs')
    gene_digest = hashlib.sha256()
    for row in gene_rows:
        gene_digest.update(row['ensembl_id'].encode() + b'\0' + row['native_gene_symbol'].encode() + b'\n')
    require(gene_digest.hexdigest() == GENE_HASH, 'gene semantic order differs')
    return units, measured, donors, gene_rows, full, halves


def guarded_arrays(inventory):
    rna_artifacts = load_json(ALIGNED / 'artifacts.json')['artifacts']
    atac_artifacts = load_json(COMMON / 'artifacts.json')
    global_artifacts = load_json(GLOBAL / 'artifacts.json')
    rna, atac, global_rna = {}, {}, {}
    for partition in PARTITIONS:
        specs = ((ALIGNED, f'RNA_{partition}_counts_uint64.npy', rna_artifacts, 'sha256', (U, G), '<u8', rna),
            (COMMON, f'ATAC_fragment_overlap_{partition}_int32.npy', atac_artifacts, 'computed_sha256', (U, R), '<i4', atac),
            (GLOBAL, f'RNA_donor_{partition}_counts_uint64.npy', global_artifacts, 'sha256', (39, G), '<u8', global_rna))
        for directory, name, artifacts, hash_key, shape, dtype, destination in specs:
            path, receipt = directory / name, artifacts[name]
            inventory[str(path)] = pin(path, receipt[hash_key], receipt['bytes'])
            array = np.load(path, mmap_mode='r', allow_pickle=False)
            require(array.shape == shape and array.dtype.str == dtype and array.flags.c_contiguous,
                    'source array shape/dtype/order differs')
            destination[partition] = array
    # Every integer conservation operation is bounded and performed on compute.
    for start in range(0, U, 8):
        end = min(U, start + 8)
        a, b, full = (rna[p][start:end] for p in ('halfA', 'halfB', 'full'))
        for counts in (a, b, full):
            require(all(sum(int(v) for v in row) < 2**53 for row in counts), 'RNA row exceeds exact integer range')
        require(np.array_equal(full, a + b), 'RNA full does not equal raw half sum')
        a, b, full = (np.asarray(atac[p][start:end], dtype=np.int64) for p in ('halfA', 'halfB', 'full'))
        require(np.all(a >= 0) and np.all(b >= 0) and np.array_equal(full, a + b), 'ATAC raw conservation differs')
    for partition in PARTITIONS:
        for donor in range(39):
            native = rna[partition][donor * 7:(donor + 1) * 7]
            require(sum(sum(int(v) for v in row) for row in native) < 2**53, 'global RNA integer range exceeded')
            require(np.array_equal(np.sum(native, axis=0, dtype=np.uint64), global_rna[partition][donor]),
                    'global RNA is not the raw donor sum over all seven lineages')
    require(np.array_equal(global_rna['full'], global_rna['halfA'] + global_rna['halfB']), 'global RNA half conservation differs')
    return rna, atac, global_rna


def input_state(units, measured, rna, global_rna):
    totals = {p: np.array([sum(int(v) for v in row) for row in rna[p]], dtype=np.float64) for p in PARTITIONS}
    global_totals = {p: np.array([sum(int(v) for v in row) for row in global_rna[p]], dtype=np.float64) for p in PARTITIONS}
    nuclei = {p: np.array([int(u['nuclei' if p == 'full' else p + '_nuclei']) for u in units], dtype=np.int64) for p in PARTITIONS}
    d = {p: np.array([int(u[p + '_autosomal_fragment_records']) for u in measured], dtype=np.int64) for p in PARTITIONS}
    require(np.all(d['full'] == d['halfA'] + d['halfB']), 'autosomal denominators do not conserve across halves')
    require(all(np.all(values >= 0) and np.all(values < 2**53) for values in d.values()), 'invalid autosomal denominator')
    return totals, global_totals, nuclei, d


def aggregate_losses(mode, losses, eligible, units, donor_order, wells, bootstrap_draws, out):
    donor_losses = np.empty((39, 4), dtype=np.float64)
    donor_rows, lineage_rows, well_rows = [], [], []
    for j, donor in enumerate(donor_order):
        indices = [i for i, u in enumerate(units) if eligible[i] and u['donor_id'] == donor]
        require(indices and np.isfinite(losses[indices]).all(), 'donor loss incomplete')
        donor_losses[j] = np.mean(losses[indices], axis=0)
        row = dict(mode=mode, donor_id=donor, well_id=units[indices[0]]['well_id'], sampled_lineages=len(indices))
        row.update({arm + '_MSE': float(donor_losses[j, a]) for a, arm in enumerate(ARMS)})
        row.update({('lineage_minus_' + arm): float(donor_losses[j, 3] - donor_losses[j, a]) for a, arm in enumerate(ARMS[:3])})
        donor_rows.append(row)
    for lineage in range(7):
        indices = [i for i in range(U) if eligible[i] and i % 7 == lineage]
        row = dict(mode=mode, lineage=units[lineage]['native_source_celltype'], donor_n=len(indices))
        row.update({arm + '_MSE': float(np.mean(losses[indices, a])) for a, arm in enumerate(ARMS)})
        row.update({('lineage_minus_' + arm): row[ARMS[3] + '_MSE'] - row[arm + '_MSE'] for arm in ARMS[:3]})
        lineage_rows.append(row)
    well_sums, well_n = [], []
    for well in wells:
        js = [j for j, donor in enumerate(donor_order) if next(u['well_id'] for u in units if u['donor_id'] == donor) == well]
        well_sums.append(np.sum(donor_losses[js], axis=0))
        well_n.append(len(js))
        row = dict(mode=mode, well_id=well, donor_n=len(js))
        row.update({arm + '_MSE': float(np.mean(donor_losses[js, a])) for a, arm in enumerate(ARMS)})
        row.update({('lineage_minus_' + arm): row[ARMS[3] + '_MSE'] - row[arm + '_MSE'] for arm in ARMS[:3]})
        well_rows.append(row)
    point = np.mean(donor_losses, axis=0)
    # Resampled well multiplicities retain their unequal donor counts: no average of well means.
    draws = np.asarray(well_sums)[bootstrap_draws].sum(axis=1) / np.asarray(well_n)[bootstrap_draws].sum(axis=1)[:, None]
    intervals = {arm: np.quantile(draws[:, a], [0.025, 0.975], method='linear').tolist() for a, arm in enumerate(ARMS)}
    contrasts = {}
    for a, arm in enumerate(ARMS[:3]):
        contrasts[arm] = {'lineage_minus_comparator': float(point[3] - point[a]),
            'conditional_95_percentile_interval': np.quantile(draws[:, 3] - draws[:, a], [0.025, 0.975], method='linear').tolist()}
    out.tsv(mode + '/donor_losses_and_effects.tsv', donor_rows)
    out.tsv(mode + '/lineage_losses_and_effects.tsv', lineage_rows)
    out.tsv(mode + '/well_losses_and_effects.tsv', well_rows)
    result = {'complete_fixed_metadata_roster': True, 'sampled_units': int(eligible.sum()), 'biological_donors': 39,
        'well_clusters': 8, 'arm_losses': {arm: float(point[a]) for a, arm in enumerate(ARMS)},
        'arm_conditional_intervals': intervals, 'contrasts': contrasts,
        'all_three_necessary_point_directions_negative': bool(all(v['lineage_minus_comparator'] < 0 for v in contrasts.values())),
        'bootstrap_seed': 20261001, 'bootstrap_draws': 10000, 'bootstrap_reweights_resampled_donor_count': True,
        'conditional_fixed_loss_not_refit_uncertainty': True, 'calibrated_population_p_value': None,
        'independent_validation_or_model_admission': False}
    out.json(mode + '/aggregate_summary.json', result)
    return result


def full_diagnostics(indices, truth, predictions, target_features, target_ridge, units, fold, lineage, out):
    pairs, substitutions, recipient_effects, checks = [], [], [], []
    n = len(indices)
    for left in range(n):
        for right in range(left + 1, n):
            observed = truth[left] - truth[right]
            for arm, clipped in predictions.items():
                estimated = clipped[left] - clipped[right]
                if arm == 'lineage_mean':
                    require(np.all(estimated == 0.0), 'fixed lineage-mean difference is not zero')
                error = estimated - observed
                pairs.append({'well_id': fold, 'lineage': lineage, 'left_donor': units[indices[left]]['donor_id'],
                    'right_donor': units[indices[right]]['donor_id'], 'arm': arm,
                    'difference_MSE': float(np.mean(error**2)), 'difference_bias': float(np.mean(error)),
                    'observed_squared_difference': float(np.mean(observed**2)),
                    'predicted_squared_difference': float(np.mean(estimated**2)),
                    'observed_predicted_cross_product': float(np.mean(observed * estimated))})
    if n >= 2:
        for arm, clipped in predictions.items():
            errors = clipped - truth
            rhs = float(2.0 / (n - 1) * np.sum(np.mean((errors - np.mean(errors, axis=0))**2, axis=1)))
            lhs = float(np.mean([row['difference_MSE'] for row in pairs if row['arm'] == arm]))
            require(abs(lhs - rhs) <= 1e-10 * (1 + abs(rhs)), 'pair error cancellation identity differs')
            checks.append({'well_id': fold, 'lineage': lineage, 'arm': arm, 'pair_n': n * (n - 1) // 2,
                'mean_pair_MSE': lhs, 'centered_error_identity_MSE': rhs, 'identity_checked_after_clipping': True})
    for recipient in range(n):
        wrong_losses = []
        for other in range(n):
            if other == recipient:
                continue
            features = target_features[recipient:recipient + 1].copy()
            features[:, 3:] = target_features[other, 3:]  # unchanged recipient technical vector
            _, wrong = predict(features, target_ridge)
            loss = float(np.mean((wrong[0] - truth[recipient])**2))
            wrong_losses.append(loss)
            substitutions.append({'well_id': fold, 'lineage': lineage, 'recipient_donor': units[indices[recipient]]['donor_id'],
                'wrong_RNA_donor': units[indices[other]]['donor_id'], 'wrong_donor_MSE': loss,
                'recipient_technical_vector_unchanged': True})
        own = float(np.mean((predictions['target_lineage_RNA'][recipient] - truth[recipient])**2))
        recipient_effects.append({'well_id': fold, 'lineage': lineage, 'recipient_donor': units[indices[recipient]]['donor_id'],
            'other_metadata_donors': len(wrong_losses), 'diagnostic_available': bool(wrong_losses),
            'own_MSE': own, 'all_wrong_donor_mean_MSE': float(np.mean(wrong_losses)) if wrong_losses else None,
            'own_minus_all_wrong_mean': own - float(np.mean(wrong_losses)) if wrong_losses else None})
    prefix = 'full/diagnostics/' + fold + '/lineage' + str(indices[0] % 7)
    if pairs:
        out.tsv(prefix + '_all_pairs.tsv', pairs)
        out.tsv(prefix + '_pair_identity_checks.tsv', checks)
    if substitutions:
        out.tsv(prefix + '_all_wrong_donors.tsv', substitutions)
    out.tsv(prefix + '_recipient_substitution_effects.tsv', recipient_effects)
    return pairs, recipient_effects


def diagnostic_aggregates(pairs, recipients, units, out):
    metrics = ('difference_MSE', 'difference_bias', 'observed_squared_difference',
               'predicted_squared_difference', 'observed_predicted_cross_product')
    groups = []
    for well in dict.fromkeys(u['well_id'] for u in units):
        for lineage in dict.fromkeys(u['native_source_celltype'] for u in units):
            subset = [r for r in pairs if r['well_id'] == well and r['lineage'] == lineage]
            for arm in ARMS:
                values = [r for r in subset if r['arm'] == arm]
                if values:
                    row = {'well_id': well, 'lineage': lineage, 'arm': arm, 'pair_n': len(values)}
                    row.update({metric: float(np.mean([r[metric] for r in values])) for metric in metrics})
                    groups.append(row)
    out.tsv('full/diagnostics/lineage_well_pair_means.tsv', groups)
    well_groups = []
    for well in dict.fromkeys(u['well_id'] for u in units):
        for arm in ARMS:
            values = [r for r in groups if r['well_id'] == well and r['arm'] == arm]
            row = {'well_id': well, 'arm': arm, 'available_lineages': len(values)}
            row.update({metric: float(np.mean([r[metric] for r in values])) for metric in metrics})
            well_groups.append(row)
    out.tsv('full/diagnostics/well_pair_means.tsv', well_groups)
    out.json('full/diagnostics/pair_summary.json', {'weighting': 'pairs within lineage/well, lineages within well, then wells',
        'arm_summaries': {arm: {metric: float(np.mean([r[metric] for r in well_groups if r['arm'] == arm]))
                               for metric in metrics} for arm in ARMS},
        'zero_difference_baseline_is_lineage_mean': True, 'signed_difference_orientation': 'earlier minus later frozen donor order',
        'pairs_not_independent_donors': True, 'no_calibration_fitted': True})
    donor_effects, lineage_effects, well_effects = [], [], []
    for donor in dict.fromkeys(u['donor_id'] for u in units):
        values = [r for r in recipients if r['recipient_donor'] == donor]
        complete = bool(values) and all(r['diagnostic_available'] for r in values)
        donor_effects.append({'donor_id': donor, 'well_id': values[0]['well_id'], 'sampled_lineages': len(values),
            'diagnostic_available': complete, 'own_minus_wrong_mean': float(np.mean([r['own_minus_all_wrong_mean'] for r in values])) if complete else None})
    for lineage in dict.fromkeys(u['native_source_celltype'] for u in units):
        values = [r for r in recipients if r['lineage'] == lineage]
        complete = all(r['diagnostic_available'] for r in values)
        lineage_effects.append({'lineage': lineage, 'donor_n': len(values), 'diagnostic_available': complete,
            'own_minus_wrong_mean': float(np.mean([r['own_minus_all_wrong_mean'] for r in values])) if complete else None})
    for well in dict.fromkeys(u['well_id'] for u in units):
        values = [r for r in donor_effects if r['well_id'] == well]
        complete = all(r['diagnostic_available'] for r in values)
        well_effects.append({'well_id': well, 'donor_n': len(values), 'diagnostic_available': complete,
            'donor_equal_own_minus_wrong_mean': float(np.mean([r['own_minus_wrong_mean'] for r in values])) if complete else None})
    out.tsv('full/diagnostics/donor_substitution_effects.tsv', donor_effects)
    out.tsv('full/diagnostics/lineage_substitution_effects.tsv', lineage_effects)
    out.tsv('full/diagnostics/well_substitution_effects.tsv', well_effects)
    complete = all(r['diagnostic_available'] for r in donor_effects)
    out.json('full/diagnostics/substitution_summary.json', {'all_metadata_recipients_reported': True,
        'all_other_eligible_same_well_lineage_donors_enumerated': True, 'diagnostic_available': complete,
        'donor_equal_own_minus_wrong_mean': float(np.mean([r['own_minus_wrong_mean'] for r in donor_effects])) if complete else None,
        'fitting_or_advancement_endpoint': False, 'causal_or_randomized_exposure': False})


def train_mode(mode, rna_partition, atac_partition, eligible, units, rna, atac, global_rna, nuclei, denominator, out, final=False):
    folds = ('all_training_donors',) if final else tuple(dict.fromkeys(u['well_id'] for u in units))
    losses = np.full((U, len(ARMS)), np.nan)
    oof = {} if final else {arm: {kind: out.empty_oof(mode + '/' + arm + '_' + kind + '_OOF.npy')
                                  for kind in ('raw', 'clipped')} for arm in ARMS}
    pair_rows, recipient_rows, fitted_metadata, final_models = [], [], [], {}
    for lineage_index in range(7):
        lineage = units[lineage_index]['native_source_celltype']
        all_indices = np.array([i for i in range(U) if eligible[i] and i % 7 == lineage_index], dtype=np.int64)
        for fold in folds:
            training = all_indices if final else np.array([i for i in all_indices if units[i]['well_id'] != fold])
            held = np.array([], dtype=np.int64) if final else np.array([i for i in all_indices if units[i]['well_id'] == fold])
            require(len(training) > K and (final or len(held) > 0), 'metadata fold insufficient or incomplete')
            # Both PCA representations use exactly these lineage-eligible training donors, each once.
            training_global_rows = training // 7
            require(len(set(training_global_rows.tolist())) == len(training), 'global donor training rows repeated')
            tech_values = technical(rna[rna_partition][training], nuclei[rna_partition][training])
            train_tech, tech_mean, tech_sd, tech_active = standard_fit(tech_values)
            train_target_log, _ = normalize_counts(rna[rna_partition][training])
            train_global_log, _ = normalize_counts(global_rna[rna_partition][training_global_rows])
            train_target_pc, target_pca = pca_fit(train_target_log)
            train_global_pc, global_pca = pca_fit(train_global_log)
            y_train = target_rows(atac[atac_partition], denominator[atac_partition], training)
            feature_sets = {'lineage_mean': np.zeros((len(training), 0)), 'technical_only': train_tech,
                'donor_global_RNA': np.column_stack((train_tech, train_global_pc)),
                'target_lineage_RNA': np.column_stack((train_tech, train_target_pc))}
            held_features, predictions = {}, {}
            if not final:
                held_tech_values = technical(rna[rna_partition][held], nuclei[rna_partition][held])
                held_tech = (held_tech_values - tech_mean) / tech_sd
                held_tech[:, ~tech_active] = 0.0
                held_target_log, _ = normalize_counts(rna[rna_partition][held])
                held_global_log, _ = normalize_counts(global_rna[rna_partition][held // 7])
                held_features = {'lineage_mean': np.zeros((len(held), 0)), 'technical_only': held_tech,
                    'donor_global_RNA': np.column_stack((held_tech, pca_apply(held_global_log, global_pca))),
                    'target_lineage_RNA': np.column_stack((held_tech, pca_apply(held_target_log, target_pca)))}
                y_held = target_rows(atac[atac_partition], denominator[atac_partition], held)
            fold_prefix = mode + '/parameters/lineage' + str(lineage_index) + '/' + fold
            arm_metadata, target_ridge = {}, None
            for a, arm in enumerate(ARMS):
                ridge, residual = fit_ridge(feature_sets[arm], y_train)
                pca = global_pca if arm == 'donor_global_RNA' else target_pca if arm == 'target_lineage_RNA' else None
                parameters = pack(ridge, None if arm == 'lineage_mean' else (tech_mean, tech_sd, tech_active), pca)
                path = out.parameters(fold_prefix + '/' + arm + '.npz', parameters)
                arm_metadata[arm] = {'parameter_path': str(path.relative_to(out.directory)),
                    'parameter_sha256': sha(path), 'features': feature_sets[arm].shape[1],
                    'normal_equation_max_abs_residual': residual,
                    'PCA_structurally_available_slots': int(pca['structurally_available']) if pca is not None else None,
                    'PCA_exact_zero_variance_slots': np.flatnonzero(~pca['score_active']).tolist() if pca is not None else None,
                    'PCA_exact_adjacent_singular_ties': pca['exact_tied_adjacent_singular_pairs'].tolist() if pca is not None else None}
                if final:
                    final_models[lineage + '/' + arm] = arm_metadata[arm]
                else:
                    raw, clipped = predict(held_features[arm], ridge)
                    losses[held, a] = np.mean((clipped - y_held)**2, axis=1)
                    write_oof(oof[arm]['raw'], held, raw)
                    write_oof(oof[arm]['clipped'], held, clipped)
                    predictions[arm] = clipped
                    if arm == 'target_lineage_RNA':
                        target_ridge = ridge
            metadata = {'mode': mode, 'lineage': lineage, 'held_well': fold,
                'training_unit_indices': training.tolist(), 'held_unit_indices': held.tolist(),
                'training_donor_ids': [units[i]['donor_id'] for i in training],
                'held_donor_ids': [units[i]['donor_id'] for i in held],
                'global_and_lineage_PCA_identical_training_donors': True,
                'raw_donor_global_sum_before_normalization': True,
                'gene_variance_scaling_before_PCA': False, 'alpha_sum_SSE': ALPHA,
                'output_regions': R, 'PCA_slots': K, 'arms': arm_metadata}
            out.json(fold_prefix + '/fold_contract.json', metadata)
            fitted_metadata.append(metadata)
            if mode == 'full' and not final:
                pairs, recipients = full_diagnostics(held, y_held, predictions, held_features['target_lineage_RNA'],
                                                       target_ridge, units, fold, lineage, out)
                pair_rows.extend(pairs)
                recipient_rows.extend(recipients)
            out.check()
            print(json.dumps({'phase': 'fitted_fixed_fold', 'mode': mode, 'lineage_index': lineage_index,
                              'fold': fold, 'training_donors': len(training), 'held_donors': len(held)}), flush=True)
    if not final:
        require(np.isfinite(losses[eligible]).all() and np.isnan(losses[~eligible]).all(), 'OOF loss identity incomplete')
        unit_rows = []
        for i, unit in enumerate(units):
            row = {'unit_index': i, 'donor_id': unit['donor_id'], 'well_id': unit['well_id'],
                'lineage': unit['native_source_celltype'], 'metadata_eligible': bool(eligible[i]),
                'reason_if_ineligible': '' if eligible[i] else 'fixed_metadata_unsampled_or_empty_half'}
            row.update({arm + '_MSE': float(losses[i, a]) if eligible[i] else None for a, arm in enumerate(ARMS)})
            unit_rows.append(row)
        out.tsv(mode + '/all273_unit_losses.tsv', unit_rows)
        for arm in ARMS:
            for kind in ('raw', 'clipped'):
                array = np.load(oof[arm][kind], mmap_mode='r', allow_pickle=False)
                for i in range(U):
                    require(np.isfinite(array[i]).all() if eligible[i] else np.isnan(array[i]).all(),
                            'saved OOF prediction identity incomplete')
                    if eligible[i] and kind == 'clipped':
                        require(np.all(array[i] >= 0), 'saved prediction projection differs')
                array._mmap.close()
        if mode == 'full':
            diagnostic_aggregates(pair_rows, recipient_rows, units, out)
    return losses, final_models


class LineageScorer:
    """Apply a final frozen model using complete native RNA counts and a fixed label.

    No ATAC, H3, clinical, histology or receiving centering inputs are accepted.
    Construction/prediction are numerical compute operations, not login checks.
    """
    def __init__(self, model_directory, lineage, arm='target_lineage_RNA'):
        self.directory = Path(model_directory)
        self.contract = load_json(self.directory / 'final_scoring_contract.json')
        require(self.contract['recipe_sha256'] == RECIPE_SHA and lineage in self.contract['lineages'] and arm in ARMS,
                'unknown model recipe, lineage or arm')
        self.lineage, self.arm = lineage, arm
        self.gene_ids = tuple(row['ensembl_id'] for row in table(self.directory / 'genes.tsv'))
        require(sha(self.directory / 'genes.tsv') == self.contract['gene_axis_file_sha256'] and len(self.gene_ids) == G,
                'scorer gene identity differs')
        entry = self.contract['models'][lineage + '/' + arm]
        path = self.directory / entry['parameter_path']
        require(path.resolve().is_relative_to(self.directory.resolve()), 'parameter path escapes model directory')
        pin(path, entry['parameter_sha256'])
        with np.load(path, allow_pickle=False) as archive:
            self.parameters = {name: archive[name] for name in archive.files}

    def predict(self, lineage_counts, nucleus_count, gene_ids, global_counts=None, global_gene_ids=None):
        require(tuple(gene_ids) == self.gene_ids, 'prediction requires complete original ordered gene IDs')
        require(isinstance(nucleus_count, (int, np.integer)) and not isinstance(nucleus_count, (bool, np.bool_)) and
                0 < nucleus_count < 2**53, 'positive exact nucleus count required')
        values = np.asarray(lineage_counts)
        require(values.shape == (G,) and values.dtype.kind in 'iu' and np.all(values >= 0), 'complete integer lineage RNA required')
        # Required RNA observability is validated even for the fixed-mean comparator.
        target_log, _ = normalize_counts(values[None, :])
        if self.arm == 'lineage_mean':
            features = np.zeros((1, 0))
        else:
            tech = technical(values[None, :], np.array([nucleus_count]))
            tech = (tech - self.parameters['technical_mean']) / self.parameters['technical_sd']
            tech[:, ~self.parameters['technical_active']] = 0.0
            features = tech
            if self.arm in ('donor_global_RNA', 'target_lineage_RNA'):
                if self.arm == 'donor_global_RNA':
                    require(global_counts is not None and global_gene_ids is not None and tuple(global_gene_ids) == self.gene_ids,
                            'global arm needs complete donor-global raw RNA and identical ordered IDs')
                    global_values = np.asarray(global_counts)
                    require(global_values.shape == (G,) and global_values.dtype.kind in 'iu' and np.all(global_values >= values),
                            'global counts must be a nonnegative integer donor sum including target lineage')
                    input_log, _ = normalize_counts(global_values[None, :])
                else:
                    input_log = target_log
                pca = {name[4:]: value for name, value in self.parameters.items() if name.startswith('PCA_')}
                features = np.column_stack((tech, pca_apply(input_log, pca)))
        raw, clipped = predict(features, self.parameters)
        return {'raw': raw[0], 'clipped': clipped[0], 'lineage': self.lineage, 'arm': self.arm,
                'region_order_sha256': REGION_HASH, 'prediction_units': self.contract['prediction_units']}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-script-sha256', required=True)
    parser.add_argument('--execution-decision', type=Path, required=True)
    parser.add_argument('--execution-decision-sha256', required=True)
    args = parser.parse_args()
    job = os.environ.get('SLURM_JOB_ID', '')
    require(re.fullmatch(r'[1-9][0-9]*', job) is not None and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'individual compute SLURM allocation required')
    source = Path(__file__).resolve()
    launcher = source.with_name('run_codex_gse296875_lineage_model.sbatch')
    inventory = {str(source): pin(source, args.expected_script_sha256)}
    inventory[str(args.execution_decision)] = pin(args.execution_decision, args.execution_decision_sha256)
    decision = load_json(args.execution_decision)
    require(decision.get('source_fitting_authorized') is True and decision.get('prespec_sha256') == RECIPE_SHA and
            decision.get('prespec_path') == str(RECIPE.relative_to(ROOT)) and
            decision.get('external_prediction_or_evaluation_authorized') is False and
            decision.get('scientific_authority_changes_or_final_adoption_authorized') is False,
            'explicit matching development-only fitting decision required')
    expected_sources = {
        'native_lineage_directory': str(ALIGNED.relative_to(ROOT)),
        'native_lineage_summary_sha256': GUARDS[ALIGNED / 'summary.json'],
        'common_fragment_directory': str(COMMON.relative_to(ROOT)),
        'common_fragment_summary_sha256': GUARDS[COMMON / 'summary.json'],
        'donor_global_directory': str(GLOBAL.relative_to(ROOT)),
        'donor_global_summary_sha256': GUARDS[GLOBAL / 'summary.json'],
        'donor_global_artifacts_sha256': GUARDS[GLOBAL / 'artifacts.json'],
    }
    require(decision.get('source_inputs') == expected_sources, 'execution-decision source inputs differ')
    expected_design = {'biological_n': 39, 'well_held_folds': 8, 'stored_units': U,
        'sampled_full_units': 272, 'two_nonempty_half_units': 270, 'RNA_genes': G,
        'ATAC_windows': R, 'RNA_PCs_per_arm': K, 'ridge_sum_SSE_alpha': ALPHA,
        'bootstrap_draws': 10000, 'bootstrap_seed': 20261001, 'arms': list(ARMS)}
    require(all(decision.get('fixed_design', {}).get(key) == value for key, value in expected_design.items()),
            'execution-decision fixed design differs')
    inventory.update({str(path): pin(path, expected) for path, expected in GUARDS.items()})
    inventory[str(launcher)] = pin(launcher, sha(launcher))
    recipe = load_json(RECIPE)
    require(recipe['authorization']['fitting_authorized_by_this_file'] is False and
            recipe['splits_and_training']['PCA_components'] == K and
            recipe['splits_and_training']['ridge']['alpha'] == ALPHA, 'frozen recipe differs')
    directory = REC / ('codex_gse296875_lineage_model_' + job)
    owned = False
    phase = 'fresh_output'
    out = None
    try:
        directory.mkdir(mode=0o750)
        owned = True
        logs = [REC / ('1595_lineage_model_' + job + suffix) for suffix in ('.out', '.err')]
        out = Output(directory, logs)
        for path in (source, launcher, RECIPE, args.execution_decision):
            out.reserve(path.stat().st_size)
            shutil.copyfile(path, out.path('executed_sources/' + path.name))
        config = io.StringIO()
        with contextlib.redirect_stdout(config):
            np.show_config()
        out.json('runtime_inventory.json', {'Python': sys.version, 'NumPy': np.__version__, 'NumPy_configuration': config.getvalue(),
            'Python_binary_sha256': sha(Path(sys.executable).resolve()), 'NumPy_path': np.__file__,
            'thread_environment': {key: os.environ.get(key) for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')},
            'model_randomness': False, 'bootstrap_generator': 'NumPy default_rng PCG64', 'bootstrap_seed': 20261001})
        phase = 'metadata_identity'
        units, measured, donors, genes, full, halves = identities(recipe)
        out.tsv('fixed_unit_identity.tsv', units)
        out.reserve((ALIGNED / 'genes.tsv').stat().st_size)
        shutil.copyfile(ALIGNED / 'genes.tsv', out.path('genes.tsv'))
        out.reserve((ALIGNED / 'native_peak_geometry.tsv').stat().st_size)
        shutil.copyfile(ALIGNED / 'native_peak_geometry.tsv', out.path('native_peak_geometry.tsv'))
        phase = 'source_array_hashes_and_integer_conservation'
        rna, atac, global_rna = guarded_arrays(inventory)
        totals, global_totals, nuclei, denominator = input_state(units, measured, rna, global_rna)
        out.json('source_inventory_before.json', inventory)
        phase = 'fixed_roster_measurability'
        states = {}
        for mode, rp, ap in MODES:
            eligible = full if mode == 'full' else halves
            rows = []
            for i, unit in enumerate(units):
                reasons = []
                if eligible[i]:
                    if totals[rp][i] <= 0: reasons.append('zero_target_lineage_RNA_denominator')
                    if global_totals[rp][i // 7] <= 0: reasons.append('zero_global_RNA_denominator')
                    if denominator[ap][i] <= 0: reasons.append('zero_target_autosomal_fragment_denominator')
                rows.append({'unit_index': i, 'donor_id': unit['donor_id'], 'well_id': unit['well_id'],
                    'lineage': unit['native_source_celltype'], 'metadata_eligible': bool(eligible[i]),
                    'required_input_and_target_defined': bool(eligible[i] and not reasons), 'failure_reasons': '|'.join(reasons)})
            states[mode] = not any(r['metadata_eligible'] and not r['required_input_and_target_defined'] for r in rows)
            out.tsv(mode + '/measurability_all273_units.tsv', rows)
        out.json('measurability_summary.json', {'modes_complete': states, 'full_fixed272_primary_available': states['full'],
            'no_value_based_row_drop': True, 'null_storage_units_not_zero_observations': True})
        results, final_models = {}, {}
        if states['full']:
            draws = np.random.default_rng(20261001).integers(0, 8, size=(10000, 8), dtype=np.int16)
            out.array('bootstrap_eight_well_indices.npy', draws)
            donor_order = [d['donor_id'] for d in donors]
            wells = tuple(recipe['source_identity']['well_donor_roster'])
            for mode, rp, ap in MODES:
                if not states[mode]:
                    results[mode] = {'complete_fixed_metadata_roster': False, 'unavailable_due_to_undefined_measurement': True,
                                     'no_value_based_removal': True}
                    continue
                phase = 'fixed_fold_models/' + mode
                eligible = full if mode == 'full' else halves
                losses, _ = train_mode(mode, rp, ap, eligible, units, rna, atac, global_rna, nuclei, denominator, out)
                results[mode] = aggregate_losses(mode, losses, eligible, units, donor_order, wells, draws, out)
            phase = 'paired_half_target_precision_descriptive'
            agreement = []
            for i in np.flatnonzero(halves):
                defined = denominator['halfA'][i] > 0 and denominator['halfB'][i] > 0
                row = {'unit_index': int(i), 'donor_id': units[i]['donor_id'], 'well_id': units[i]['well_id'],
                       'lineage': units[i]['native_source_celltype'], 'agreement_available': bool(defined)}
                if defined:
                    a = target_rows(atac['halfA'], denominator['halfA'], np.array([i]))[0]
                    b = target_rows(atac['halfB'], denominator['halfB'], np.array([i]))[0]
                    row.update({'half_target_MSE': float(np.mean((a - b)**2)), 'half_target_mean_difference': float(np.mean(a - b)),
                                'halfA_mean': float(np.mean(a)), 'halfB_mean': float(np.mean(b))})
                else:
                    row.update({'half_target_MSE': None, 'half_target_mean_difference': None, 'halfA_mean': None, 'halfB_mean': None})
                agreement.append(row)
            out.tsv('paired_half_target_agreement_all270.tsv', agreement)
            phase = 'final_all_sampled_donor_lineage_scorers'
            _, final_models = train_mode('final_full', 'full', 'full', full, units, rna, atac, global_rna, nuclei, denominator, out, final=True)
            out.json('final_scoring_contract.json', {'recipe_sha256': RECIPE_SHA, 'producer_sha256': args.expected_script_sha256,
                'lineages': recipe['source_identity']['ordered_author_lineages'], 'arms': ARMS, 'models': final_models,
                'gene_axis_file_sha256': GUARDS[ALIGNED / 'genes.tsv'], 'ordered_gene_ID_symbol_sha256': GENE_HASH,
                'ordered_region_axis_sha256': REGION_HASH, 'RNA_genes': G, 'counted_windows': R,
                'prediction_units': 'ln(1+1e6*unique_fragment_body_overlap_count/autosomal_unique_fragment_count)',
                'required_inputs': 'complete ordered target-lineage integer RNA, positive nucleus count, fixed compatible lineage label; global arm additionally full raw donor-global RNA',
                'prediction_forbidden_inputs': ['ATAC', 'H3', 'clinical', 'histology', 'genotype', 'allelic'],
                'projection': 'elementwise maximum(raw,0) before evaluation', 'source_exposed_development_only': True,
                'external_evaluation_or_model_admission': False, 'source_filter_and_annotation_compatibility_required': True})
        else:
            results['full'] = {'complete_fixed_metadata_roster': False, 'primary_unmeasurable': True,
                               'no_fits_run_to_rescue_roster': True}
        phase = 'final_source_guards_and_accounting'
        after = {path: pin(path, record['sha256'], record['bytes']) for path, record in inventory.items()}
        require(after == inventory, 'sources changed during development execution')
        out.json('source_inventory_after.json', after)
        artifacts = {str(p.relative_to(directory)): {'bytes': p.stat().st_size, 'sha256': sha(p)}
                     for p in directory.rglob('*') if p.is_file()}
        out.json('artifacts.json', artifacts)
        out.json('summary.json', {'development_execution_completed': True, 'SLURM_job_id': job, 'recipe_sha256': RECIPE_SHA,
            'producer_sha256': args.expected_script_sha256, 'execution_decision_sha256': args.execution_decision_sha256,
            'complete_primary_available': states['full'], 'results': results, 'final_scorers': len(final_models),
            'full_sampled_units': 272, 'two_half_metadata_eligible_units': 270, 'biological_donors': 39, 'well_folds': 8,
            'all_regions_direct_ridge': R, 'PCA_gene_axis': G, 'PCA_components': K, 'alpha_sum_SSE': ALPHA,
            'both_PCA_arms_identical_training_donor_rows': True, 'prediction_projection_applied_before_all_metrics': True,
            'conditional_bootstrap_not_refit_or_population_pvalue': True, 'source_inventory_before_after_identical': True,
            'saved_bytes_before_summary': out.check(), 'saved_output_cap_bytes_inclusive': CAP,
            'source_exposed_development': True, 'native_ChromatinAssay_count_semantics_equivalence_claimed': False,
            'fixed117_programs_modified': False, 'independently_validated_or_usable_model_claimed': False,
            'clinical_H3_genotype_allelic_or_receiving_values_read': False, 'external_model_evaluation': False,
            'external_acquisition': False, 'frozen_recipe_unchanged': True})
        print(json.dumps({'completed': True, 'complete_primary_available': states['full'], 'saved_bytes': out.check(),
                          'summary_sha256': sha(directory / 'summary.json')}), flush=True)
        out.check()
    except Exception as exc:
        failure = {'completed': False, 'phase': phase, 'exception_type': type(exc).__name__,
                   'reason': str(exc) if isinstance(exc, RuntimeError) else type(exc).__name__,
                   'partials_preserved': True, 'fixed_roster_or_axis_changed': False, 'no_model_admission': True}
        if owned:
            try:
                with (directory / 'failure.json').open('x') as stream:
                    json.dump(failure, stream, indent=2)
                    stream.write('\n')
            except Exception:
                pass
        print(json.dumps(failure), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
