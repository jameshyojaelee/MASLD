#!/usr/bin/env python3
"""Fixed donor-global source development; numerical execution requires SLURM."""
import argparse
import ast
import contextlib
import csv
import hashlib
import io
import json
import os
from pathlib import Path
import shutil
import struct
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
HERE = ROOT / 'scripts/analysis/in_biopsy_allelic/rec_codex'
INPUT = REC / 'codex_gse296875_donor_global_raw_22000941'
PRESPEC = REC / '20261001_codex_gse296875_donor_global_rna_atac_prespec.json'
DECISION = REC / '20261001_codex_donor_global_source_execution_decision.json'
LAUNCHER = HERE / 'run_codex_gse296875_donor_global_model.sbatch'
PARTITIONS = ('full', 'halfA', 'halfB')
LABELS = ('B cells', 'Cholangiocytes', 'Hepatocytes', 'Kupffer', 'LSEC', 'Mesenchymal', 'NK-T')
ARMS = ('donor_mean', 'technical_only', 'donor_global_RNA_PC', 'technical_source_fractions')
Q = (0, 3, 23, 9)
DIRECTIONS = (('full_RNA_to_full_ATAC', 'full', 'full'),
              ('halfA_RNA_to_halfB_ATAC', 'halfA', 'halfB'),
              ('halfB_RNA_to_halfA_ATAC', 'halfB', 'halfA'))
DENOMS = ('primary_fragment_records', 'autosomal_fragment_records', 'union_overlapped_fragments',
          'readSupport_diagnostic', 'summed_peak_counts')
N, G, R, PCS, ALPHA, CHUNK = 39, 36601, 306706, 20, 10.0, 8192
CAP = 8 * 1024**3
ACTIVE_OUTPUT = None
ACTIVE_CAP = CAP
GUARDS = {
    PRESPEC: '88b9e4239b08cc1ccbe78b092e8eacd395778453a496fd8915b4acc9e2729ae1',
    DECISION: '5a8e616c7a4f4f7cbe7c88b8c9c64e8939eff009bef4afc1d53568f318bad326',
    INPUT / 'summary.json': '9684d83681b898ac507680e2df39a4781b42daccbd2abfac6342faa88031fb19',
    INPUT / 'artifacts.json': '5629741635252c1869cc68da59c365833ad8912fa7f46d62fef25335beead111',
    INPUT / 'manifest.json': 'd0f9e816e48d4784672242506eac56b074a0737c450f529c0c9836f93b8807c0',
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(path):
    require(path.is_file() and not path.is_symlink(), 'Missing regular source: ' + str(path))
    s = path.stat()
    return dict(bytes=s.st_size, inode=s.st_ino, device=s.st_dev, mtime_ns=s.st_mtime_ns, ctime_ns=s.st_ctime_ns)


def size(output):
    return sum(p.stat().st_size for p in output.rglob('*') if p.is_file())


def cap(output):
    require(size(output) <= CAP, 'Inclusive saved output exceeds 8 GiB')


def write_json(path, value):
    partial = path.with_name(path.name + '.partial')
    require(not path.exists() and not partial.exists(), 'Refusing JSON overwrite')
    with partial.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    if ACTIVE_OUTPUT is not None:
        require(size(ACTIVE_OUTPUT) <= ACTIVE_CAP, 'Inclusive output cap exceeded before JSON completion')
    partial.rename(path)


def write_table(path, rows):
    with path.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def read_table(path, columns):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(set(columns) <= set(reader.fieldnames or ()), 'Source metadata schema differs')
        return [{key: row[key] for key in columns} for row in reader]


def header(path):
    with path.open('rb') as handle:
        require(handle.read(8) == b'\x93NUMPY\x01\x00', 'NPY version differs')
        length = struct.unpack('<H', handle.read(2))[0]
        require(0 < length <= 65536, 'NPY header exceeds bound')
        parsed = ast.literal_eval(handle.read(length).decode('latin1').strip())
        require(set(parsed) == {'descr', 'shape', 'fortran_order'}, 'NPY header schema differs')
        offset = handle.tell()
        handle.seek(0)
        digest = hashlib.sha256(handle.read(offset)).hexdigest()
    return dict(dtype=parsed['descr'], shape=list(parsed['shape']), fortran_order=parsed['fortran_order'],
                offset=offset, header_sha256=digest)


def save_array(path, values):
    require(not path.exists(), 'Refusing array overwrite')
    partial = path.with_name(path.name + '.partial')
    with partial.open('xb') as handle:
        np.save(handle, values, allow_pickle=False)
    partial.rename(path)


def new_map(path, shape):
    require(not path.exists(), 'Refusing mapped output overwrite')
    return np.lib.format.open_memmap(path, mode='w+', dtype='<f8', shape=shape)


def close_map(value):
    value.flush() if value.mode != 'r' else None
    value._mmap.close()


def scale_fit(values):
    mean = values.mean(axis=0)
    sd = values.std(axis=0, ddof=1)
    constant = np.all(values == values[0], axis=0)
    sd[constant] = 0.0
    return mean, sd


def scaled(values, mean, sd):
    result = np.zeros_like(values, dtype=np.float64)
    np.divide(values - mean, sd, out=result, where=sd > 0)
    require(np.isfinite(result).all(), 'Nonfinite standardized predictors')
    return result


def fit_transform(x, tech, fractions, train):
    """Every fitted quantity below sees training rows only."""
    require(len(train) > 1 and len(set(map(int, train))) == len(train), 'Invalid training identity')
    xt = x[train]
    gene_mean = xt.mean(axis=0)
    centered = xt - gene_mean
    constant_genes = np.all(xt == xt[0], axis=0)
    centered[:, constant_genes] = 0.0
    _, singular, vt = np.linalg.svd(centered, full_matrices=False)
    available = min(PCS, len(singular), len(train) - 1)
    loadings = np.zeros((x.shape[1], PCS), dtype=np.float64)
    loadings[:, :available] = vt[:available].T
    # These coordinates are mathematically zero; remove LAPACK roundoff only.
    loadings[constant_genes] = 0.0
    sign_indices = np.argmax(np.abs(loadings), axis=0)
    for component in range(available):
        if loadings[sign_indices[component], component] < 0:
            loadings[:, component] *= -1
    scores = (xt - gene_mean) @ loadings
    pc_mean, pc_sd = scale_fit(scores)
    pc_sd[available:] = 0.0
    tech_mean, tech_sd = scale_fit(tech[train])
    fraction_mean, fraction_sd = scale_fit(fractions[train])
    return dict(gene_mean=gene_mean, loadings=loadings, singular_values=singular,
                available_PC_slots=np.array(available), sign_gene_indices=sign_indices,
                PC_mean=pc_mean, PC_sd=pc_sd, technical_mean=tech_mean, technical_sd=tech_sd,
                fraction_mean=fraction_mean, fraction_sd=fraction_sd, training_indices=np.array(train))


def apply_transform(x, tech, fractions, transform):
    t = scaled(tech, transform['technical_mean'], transform['technical_sd'])
    p = scaled((x - transform['gene_mean']) @ transform['loadings'], transform['PC_mean'], transform['PC_sd'])
    f = scaled(fractions, transform['fraction_mean'], transform['fraction_sd'])
    return (np.zeros((len(x), 0)), t, np.column_stack((t, p)), np.column_stack((t, f)))


def ridge_operator(xtrain):
    mean = xtrain.mean(axis=0)
    centered = xtrain - mean
    operator = np.linalg.solve(centered.T @ centered + ALPHA * np.eye(xtrain.shape[1]), centered.T)
    return mean, operator


def ridge_head(ytrain, feature_mean, operator):
    target_mean = ytrain.mean(axis=0)
    beta = operator @ (ytrain - target_mean)
    intercept = target_mean - feature_mean @ beta
    return intercept, beta


def bootstrap_means(losses, well_indices, draws):
    # Sum donor contributions, then divide by the draw's total number of donors.
    sums = np.stack([losses[index].sum(axis=0) for index in well_indices])
    counts = np.array([len(index) for index in well_indices])
    return sums[draws].sum(axis=1) / counts[draws].sum(axis=1)[:, None]


def selfchecks():
    # Extreme changes in held inputs/targets must leave transforms and heads invariant.
    x = np.array([[0., 2., 1.], [1., 2., 2.], [3., 2., 0.], [7., 99., 8.]])
    tech = np.array([[2., 1., 0.], [2., 2., 1.], [2., 3., 2.], [9., 8., 7.]])
    fractions = np.column_stack((tech[:, 1:], tech[:, 1:], tech[:, 1:]))
    train = np.array([0, 1, 2])
    a = fit_transform(x, tech, fractions, train)
    changed = x.copy()
    changed[3] = 1e9
    changed_tech, changed_frac = tech.copy(), fractions.copy()
    changed_tech[3], changed_frac[3] = 1e9, 1e9
    b = fit_transform(changed, changed_tech, changed_frac, train)
    require(all(np.array_equal(a[k], b[k]) for k in a), 'Held inputs entered fitted transforms')
    design = apply_transform(x, tech, fractions, a)[2]
    require(np.all(design[train, 0] == 0) and a['technical_sd'][0] == 0,
            'Training constant technical variable not mapped to zero')
    require(np.all(a['loadings'][1] == 0) and np.all(design[:, 3 + 2:] == 0),
            'Constant RNA gene or unavailable PC slots mishandled')
    y = np.array([[1., 0.], [2., 1.], [4., 2.], [99., 99.]])
    center, operator = ridge_operator(design[train])
    head = ridge_head(y[train], center, operator)
    y[3] = -1e9
    other = ridge_head(y[train], center, operator)
    require(all(np.array_equal(u, v) for u, v in zip(head, other)), 'Held targets entered fitted heads')
    # Independently hand-derived sum-SSE ridge with an unpenalized intercept.
    centered_x = np.array([[-1.], [0.], [1.]])
    y_check = np.array([[-2.], [1.], [4.]])
    ridge_mean, ridge_map = ridge_operator(centered_x)
    a0, beta0 = ridge_head(y_check, ridge_mean, ridge_map)
    a7, beta7 = ridge_head(y_check + 7, ridge_mean, ridge_map)
    require(np.allclose(beta0, 0.5, rtol=0, atol=1e-15) and np.array_equal(a0, np.array([1.])) and
            np.array_equal(beta7, beta0) and np.array_equal(a7, np.array([8.])),
            'Sum-SSE alpha10 or unpenalized intercept differs from hand-derived result')
    # Unequal wells: the mixed draw is (1+3+9)/3, not mean(2,9).
    losses = np.array([[1.], [3.], [9.]])
    draws = np.array([[0, 1], [0, 0], [1, 1]])
    actual = bootstrap_means(losses, (np.array([0, 1]), np.array([2])), draws)[:, 0]
    require(np.array_equal(actual, np.array([13 / 3, 2., 9.])), 'Cluster bootstrap donor weighting differs')
    require(np.array_equal(np.maximum(np.array([-1., 0., 2.]), 0), np.array([0., 0., 2.])), 'Fixed clipping differs')
    return dict(passed=True, held_inputs_and_targets_excluded=True, constant_features_zero=True,
                fixed_20_PC_slots_padded=True, unequal_well_bootstrap_donor_weighting=True, clipping=True,
                hand_derived_ridge_beta_0_5_intercept_1_and_shift_7=True,
                source_raw_half_conservation='Separately required on all real source arrays before fitting.')


def source_metadata(recipe, artifacts, summary):
    fields = ('donor_index', 'donor_id', 'well_id', 'stored_lineage_buckets', 'sampled_lineage_buckets') + tuple(
        p + '_' + k for p in PARTITIONS for k in ('nuclei', 'RNA_total_UMIs') + DENOMS)
    donors = read_table(INPUT / 'donors.tsv', fields)
    require(len(donors) == N and [int(d['donor_index']) for d in donors] == list(range(N)), 'Donor axis differs')
    ids = [d['donor_id'] for d in donors]
    roster = recipe['source_identity']['well_donor_roster']
    require(ids == sorted({d for group in roster.values() for d in group}, key=int), 'Exact39 donor IDs/order differ')
    wells = tuple(sorted(roster))
    folds = []
    for well in wells:
        held = [i for i, d in enumerate(donors) if d['well_id'] == well]
        train = [i for i, d in enumerate(donors) if d['well_id'] != well]
        require([ids[i] for i in held] == roster[well] and len(train) == recipe['splits_and_model']['training_donor_counts'][well],
                'Whole-well membership differs')
        require(not set(train) & set(held) and set(train) | set(held) == set(range(N)), 'Train/held donor exclusion differs')
        folds.append(dict(well_id=well, train_indices=train, held_indices=held,
                          train_donor_ids=[ids[i] for i in train], held_donor_ids=[ids[i] for i in held]))
    require(len(wells) == 8 and sorted(i for fold in folds for i in fold['held_indices']) == list(range(N)), 'Incomplete eight folds')
    for d in donors:
        require(int(d['stored_lineage_buckets']) == 7 and int(d['sampled_lineage_buckets']) == (6 if d['donor_id'] == '733' else 7),
                'Missing sampled lineage source bucket')
        for field in ('nuclei', 'RNA_total_UMIs') + DENOMS:
            require(int(d['full_' + field]) == int(d['halfA_' + field]) + int(d['halfB_' + field]), 'Raw donor metadata full=A+B differs')
        for p in PARTITIONS:
            for field in ('nuclei', 'RNA_total_UMIs', 'autosomal_fragment_records'):
                require(int(d[p + '_' + field]) > 0, 'Undefined fixed donor input/target; complete result withheld')
    for p in PARTITIONS:
        for field in ('nuclei', 'RNA_total_UMIs') + DENOMS:
            require(sum(int(d[p + '_' + field]) for d in donors) == summary['global_raw_metadata'][p][field], 'Global raw metadata differs')
    fractions = read_table(INPUT / 'lineage_nucleus_fractions.tsv',
                          ('donor_index', 'donor_id', 'well_id', 'unit_index', 'native_source_celltype', 'imputed_measurement') +
                          tuple(p + '_' + field for p in PARTITIONS for field in ('nuclei', 'nucleus_fraction')))
    require(len(fractions) == 273 and [(u['donor_id'], u['native_source_celltype']) for u in fractions] ==
            [(d, label) for d in ids for label in LABELS] and [int(u['unit_index']) for u in fractions] == list(range(273)),
            'Complete native fraction Cartesian identity differs')
    require(all(u['imputed_measurement'] == 'False' for u in fractions), 'Imputed source fraction measurement')
    f = {}
    for p in PARTITIONS:
        values = np.zeros((N, 7), dtype=np.float64)
        for i, d in enumerate(donors):
            group = fractions[7 * i:7 * (i + 1)]
            require(all(u['donor_index'] == str(i) and u['well_id'] == d['well_id'] for u in group), 'Fraction donor/well join differs')
            require(sum(int(u[p + '_nuclei']) for u in group) == int(d[p + '_nuclei']), 'Fraction roster is incomplete')
            for k, u in enumerate(group):
                value = int(u[p + '_nuclei']) / int(d[p + '_nuclei'])
                require(value == float(u[p + '_nucleus_fraction']), 'Frozen input fraction differs')
                values[i, k] = value
        require(np.allclose(values.sum(axis=1), 1, rtol=0, atol=1e-15), 'Seven source nucleus fractions do not sum to one')
        f[p] = values[:, :6].copy()
    genes = read_table(INPUT / 'genes.tsv', ('gene_index', 'ensembl_id', 'native_gene_symbol'))
    require(len(genes) == G and [int(g['gene_index']) for g in genes] == list(range(G)) and
            len({g['ensembl_id'] for g in genes}) == G, 'Complete native ordered gene axis differs')
    h = hashlib.sha256()
    for g in genes:
        h.update(g['ensembl_id'].encode() + b'\0' + g['native_gene_symbol'].encode() + b'\n')
    require(h.hexdigest() == recipe['source_identity']['rna_gene_id_symbol_order_sha256'], 'Gene identity/symbol codec differs')
    require(artifacts['native_peak_geometry.tsv']['sha256'] == recipe['source_identity']['native_peak_geometry_sha256'], 'Frozen window axis differs')
    return donors, folds, f


def check_raw_arrays(donors, artifacts):
    records, xs, technical = [], {}, {}
    for modality, columns in (('RNA', G), ('ATAC', R)):
        maps = []
        for p in PARTITIONS:
            name = f'{modality}_donor_{p}_counts_uint64.npy'
            path = INPUT / name
            expected = artifacts[name]
            before = fingerprint(path)
            h = header(path)
            require(h == expected['header'] and h['shape'] == [N, columns] and h['dtype'] == '<u8' and
                    h['fortran_order'] is False and before['bytes'] == expected['bytes'] == 128 + N * columns * 8,
                    'Aggregate NPY header/dtype/size differs')
            require(sha(path) == expected['sha256'] and fingerprint(path) == before, 'Aggregate array hash/stat differs')
            maps.append(np.load(path, mmap_mode='r', allow_pickle=False))
            records.append(dict(path=str(path), sha256=expected['sha256'], stat_before=before, header=h))
        totals = [np.zeros(N, dtype=np.uint64) for _ in PARTITIONS]
        for start in range(0, columns, CHUNK):
            stop = min(columns, start + CHUNK)
            blocks = [m[:, start:stop] for m in maps]
            require(np.all(blocks[1] <= np.iinfo(np.uint64).max - blocks[2]) and
                    np.array_equal(blocks[0], blocks[1] + blocks[2]), 'Persisted raw feature full=A+B differs')
            for i, block in enumerate(blocks):
                require(int(block.max()) < 2**53, 'Raw integer magnitude cannot be represented exactly in float64')
                addition = block.sum(axis=1, dtype=np.uint64)
                require(np.all(totals[i] <= np.iinfo(np.uint64).max - addition), 'Source row sum uint64 overflow')
                totals[i] += addition
        for i, p in enumerate(PARTITIONS):
            field = p + ('_RNA_total_UMIs' if modality == 'RNA' else '_summed_peak_counts')
            require(totals[i].tolist() == [int(d[field]) for d in donors] and int(totals[i].max()) < 2**53,
                    'Persisted source count row total differs')
            if modality == 'RNA':
                raw = np.array(maps[i], dtype=np.float64)
                xs[p] = np.log1p(raw * (1e6 / totals[i].astype(np.float64))[:, None])
                technical[p] = np.column_stack((np.log1p(totals[i].astype(np.float64)),
                    np.log1p([int(d[p + '_nuclei']) for d in donors]), (raw > 0).sum(axis=1) / G))
                require(np.isfinite(xs[p]).all() and np.isfinite(technical[p]).all(), 'Nonfinite fixed RNA transform')
        del blocks
        for m in maps:
            close_map(m)
    return records, xs, technical


def save_transform(directory, transform):
    # Numeric arrays only; predictor names and identities remain in JSON/TSV.
    path = directory / 'transform.npz'
    with path.open('xb') as handle:
        np.savez(handle, **transform)


def fit_heads(output, directory, designs, target, train, held, oof, filled):
    require(not set(map(int, train)) & set(map(int, held)), 'Held donors entered head training')
    for arm_index, (arm, design) in enumerate(zip(ARMS, designs)):
        require(design.shape == (N, Q[arm_index]), 'Frozen head design dimension differs')
        feature_mean, operator = ridge_operator(design[train]) if Q[arm_index] else (np.zeros(0), np.zeros((0, len(train))))
        coefficient_path = directory / (arm + '_coefficients.npy')
        intercept_path = directory / (arm + '_intercept.npy')
        coefficient = new_map(coefficient_path, (Q[arm_index], R)) if Q[arm_index] else None
        intercepts = new_map(intercept_path, (R,))
        save_array(directory / (arm + '_head_feature_mean.npy'), feature_mean)
        for start in range(0, R, CHUNK):
            stop = min(R, start + CHUNK)
            intercept, beta = ridge_head(target[train, start:stop], feature_mean, operator)
            require(np.isfinite(intercept).all() and np.isfinite(beta).all(), 'Nonfinite fitted ridge head')
            if coefficient is not None:
                coefficient[:, start:stop] = beta
            intercepts[start:stop] = intercept
            if oof is not None:
                raw = intercept[None, :] + designs[arm_index][held] @ beta
                require(np.isfinite(raw).all(), 'Nonfinite held prediction')
                oof[arm][0][np.ix_(held, np.arange(start, stop))] = raw
                oof[arm][1][np.ix_(held, np.arange(start, stop))] = np.maximum(raw, 0.0)
        if coefficient is not None:
            close_map(coefficient)
        else:
            save_array(coefficient_path, np.zeros((0, R), dtype=np.float64))
        close_map(intercepts)
        if oof is not None:
            require(not filled[arm][held].any(), 'Held donor predicted more than once')
            filled[arm][held] = True
            for value in oof[arm]:
                value.flush()
        cap(output)


def load_deployment(bundle):
    """Load saved full-source parameters, never targets or a new fitted transform."""
    directory = bundle / 'full_source_deployment'
    with np.load(directory / 'transform.npz', allow_pickle=False) as stored:
        transform = {key: stored[key] for key in stored.files}
    require(transform['loadings'].shape == (G, PCS) and transform['training_indices'].tolist() == list(range(N)),
            'Full-source deployment transform differs')
    genes = read_table(bundle / 'genes.tsv', ('gene_index', 'ensembl_id'))
    require([int(g['gene_index']) for g in genes] == list(range(G)), 'Deployment gene identity axis differs')
    return dict(transform=transform, gene_ids=tuple(g['ensembl_id'] for g in genes), directory=directory)


def predict_donor_global(raw_rna_UMIs, nucleus_count, qualified_gene_ids, deployment, assay_qualification):
    """One admitted raw snRNA donor; returns three operational clipped profiles.

    The caller supplies an explicit reviewed assay admission record. This routine
    checks its input contract, not the truth of claimed biological provenance.
    No lineage labels, targets, chromatin or clinical predictor fields are used.
    """
    require(os.environ.get('SLURM_JOB_ID', '').isdecimal(), 'Prediction requires compute in this environment')
    raw = np.asarray(raw_rna_UMIs)
    require(raw.shape == (G,) and raw.dtype.kind in 'iu' and (raw >= 0).all(), 'Complete nonnegative integer RNA vector required')
    require(tuple(qualified_gene_ids) == deployment['gene_ids'], 'Ordered stable gene IDs differ; no missing-gene filling')
    require(isinstance(nucleus_count, int) and not isinstance(nucleus_count, bool) and nucleus_count > 0, 'Positive pooled nucleus count required')
    require(assay_qualification.get('admitted') is True and assay_qualification.get('assay') == 'raw_snRNA_UMI' and
            isinstance(assay_qualification.get('admission_record_sha256'), str) and
            len(assay_qualification['admission_record_sha256']) == 64,
            'Reviewed assay/reference/filter admission record required; arbitrary RNA modalities unadmitted')
    total = sum(map(int, raw))
    require(0 < total < 2**53 and int(raw.max()) < 2**53, 'Undefined or unsupported RNA integer denominator')
    x = np.log1p(raw.astype(np.float64)[None, :] * (1e6 / total))
    tech = np.array([[np.log1p(total), np.log1p(nucleus_count), int((raw > 0).sum()) / G]])
    transform = deployment['transform']
    t = scaled(tech, transform['technical_mean'], transform['technical_sd'])
    p = scaled((x - transform['gene_mean']) @ transform['loadings'], transform['PC_mean'], transform['PC_sd'])
    designs = (np.zeros((1, 0)), t, np.column_stack((t, p)))
    predictions = {}
    for arm_index, arm in enumerate(ARMS[:3]):
        directory = deployment['directory']
        intercept = np.load(directory / (arm + '_intercept.npy'), mmap_mode='r', allow_pickle=False)
        beta = np.load(directory / (arm + '_coefficients.npy'), mmap_mode='r', allow_pickle=False)
        require(intercept.shape == (R,) and beta.shape == (Q[arm_index], R), 'Saved operational head dimension differs')
        values = np.empty(R, dtype=np.float64)
        for start in range(0, R, CHUNK):
            stop = min(R, start + CHUNK)
            values[start:stop] = np.maximum(intercept[start:stop] + (designs[arm_index] @ beta[:, start:stop])[0], 0)
        require(np.isfinite(values).all(), 'Nonfinite deployment predictions')
        close_map(intercept)
        # Empty mean-head arrays have no payload and are not memory mapped.
        if isinstance(beta, np.memmap):
            close_map(beta)
        predictions[arm] = values
    return predictions, dict(genes=G, windows=R, RNA_total_UMIs=total, nuclei=nucleus_count,
             log_base='natural', units='ln1p fragments per million all autosomal unique fragments',
             elementwise_nonnegative_clipping=True, no_refit=True, no_target_or_lineage_predictor_inputs=True,
             assay_qualification=assay_qualification, provenance_truth_established_by_this_callable=False)


def check_deployment_reconstruction(output, xs, tech, fractions, donors):
    deployment = load_deployment(output)
    designs = apply_transform(xs['full'], tech['full'], fractions['full'], deployment['transform'])
    raw = np.load(INPUT / 'RNA_donor_full_counts_uint64.npy', mmap_mode='r', allow_pickle=False)
    maxdiff = {arm: 0.0 for arm in ARMS[:3]}
    for i, donor in enumerate(donors):
        predicted, _ = predict_donor_global(raw[i], int(donor['full_nuclei']), deployment['gene_ids'], deployment,
                    dict(admitted=True, assay='raw_snRNA_UMI', admission_record_sha256=GUARDS[PRESPEC],
                         scope='Frozen observed source input only; no receiving admission.'))
        for j, arm in enumerate(ARMS[:3]):
            directory = deployment['directory']
            intercept = np.load(directory / (arm + '_intercept.npy'), mmap_mode='r', allow_pickle=False)
            beta = np.load(directory / (arm + '_coefficients.npy'), mmap_mode='r', allow_pickle=False)
            for start in range(0, R, CHUNK):
                stop = min(R, start + CHUNK)
                direct = np.maximum(intercept[start:stop] + (designs[j][i:i+1] @ beta[:, start:stop])[0], 0)
                difference = float(np.abs(direct - predicted[arm][start:stop]).max())
                maxdiff[arm] = max(maxdiff[arm], difference)
                require(np.allclose(direct, predicted[arm][start:stop], rtol=1e-12, atol=1e-12), 'Saved inference callable reconstruction differs')
            close_map(intercept)
            if isinstance(beta, np.memmap):
                close_map(beta)
    close_map(raw)
    return dict(passed=True, all_source_donors_checked=N, operational_heads=3, output_windows=R,
                compared_entries=N * 3 * R, max_absolute_difference=maxdiff, relative_tolerance=1e-12, absolute_tolerance=1e-12,
                saved_parameters_only=True, no_refit_or_accuracy_selection=True)


def inference_main(args):
    require(os.environ.get('SLURM_JOB_ID', '').isdecimal(), 'Inference requires compute')
    require(args.expected_helper_sha and sha(Path(__file__)) == args.expected_helper_sha, 'Inference helper source drift')
    require(args.bundle and args.bundle_summary_sha and args.bundle_artifacts_sha and args.admission_json and args.admission_sha,
            'Prediction requires immutable completed bundle and separate reviewed input admission')
    bundle = Path(args.bundle)
    require(sha(bundle / 'summary.json') == args.bundle_summary_sha and sha(bundle / 'artifacts.json') == args.bundle_artifacts_sha,
            'Completed model bundle metadata differs')
    require(json.loads((bundle / 'summary.json').read_text())['source_donor_global_model_completed'], 'Incomplete model bundle')
    artifacts = json.loads((bundle / 'artifacts.json').read_text())
    used = ['genes.tsv', 'full_source_deployment/transform.npz', 'full_source_deployment/interface.json']
    used += [f'full_source_deployment/{arm}_{kind}.npy' for arm in ARMS[:3] for kind in ('intercept','coefficients')]
    snapshots = {}
    for name in used:
        path = bundle / name
        before = fingerprint(path)
        require(before['bytes'] == artifacts[name]['bytes'] and sha(path) == artifacts[name]['sha256'] and fingerprint(path) == before,
                'Immutable inference parameter/axis differs')
        snapshots[name] = before
    admission_path = Path(args.admission_json)
    require(sha(admission_path) == args.admission_sha, 'Input admission record differs')
    admission = json.loads(admission_path.read_text())
    query_path = Path(args.predict_json)
    require(admission.get('prediction_authorized') is True and admission.get('input_sha256') == sha(query_path) and
            admission.get('model_summary_sha256') == args.bundle_summary_sha, 'Prediction input not admitted to this frozen model')
    require(query_path.stat().st_size <= 2 * 1024**2, 'Query JSON exceeds input bound')
    query_before = fingerprint(query_path)
    query = json.loads(query_path.read_text())
    require(set(query) == {'raw_rna_UMIs','nucleus_count','ensembl_ids'}, 'Unexpected prediction fields')
    global np
    import numpy as np
    predictions, metadata = predict_donor_global(query['raw_rna_UMIs'], query['nucleus_count'], query['ensembl_ids'],
            load_deployment(bundle), dict(admitted=True, assay=admission.get('assay'), admission_record_sha256=args.admission_sha))
    output = REC / ('codex_gse296875_donor_global_prediction_' + os.environ['SLURM_JOB_ID'])
    output.mkdir(mode=0o750)
    global ACTIVE_OUTPUT, ACTIVE_CAP
    ACTIVE_OUTPUT, ACTIVE_CAP = output, 16 * 1024**2
    for arm, values in predictions.items():
        save_array(output / (arm + '_clipped_log1p_FPM.npy'), values)
    for name, before in snapshots.items():
        require(fingerprint(bundle / name) == before and sha(bundle / name) == artifacts[name]['sha256'], 'Deployment parameter changed during prediction')
    require(fingerprint(query_path) == query_before and sha(query_path) == admission['input_sha256'] and
            sha(admission_path) == args.admission_sha, 'Input/admission changed during prediction')
    write_json(output / 'summary.json', dict(prediction_completed=True, model_summary_sha256=args.bundle_summary_sha,
               model_artifacts_sha256=args.bundle_artifacts_sha, input_sha256=admission['input_sha256'], **metadata))
    require(size(output) <= 16 * 1024**2, 'Prediction output exceeds16MiB')


def evaluate(output, directory, target, donors, folds, draws):
    require(draws.shape == (10000, 8) and draws.dtype.kind in 'iu' and
            np.all((draws >= 0) & (draws < 8)), 'Frozen bootstrap draw shape/index bounds differ')
    require(len(folds) == 8 and sorted(i for fold in folds for i in fold['held_indices']) == list(range(N)),
            'Bootstrap wells are not the exact held donor partition')
    losses = np.zeros((N, len(ARMS)))
    projected_counts = []
    for arm_index, arm in enumerate(ARMS):
        raw = np.load(directory / (arm + '_OOF_raw.npy'), mmap_mode='r', allow_pickle=False)
        clipped = np.load(directory / (arm + '_OOF_clipped.npy'), mmap_mode='r', allow_pickle=False)
        negative = 0
        for start in range(0, R, CHUNK):
            stop = min(R, start + CHUNK)
            require(np.isfinite(raw[:, start:stop]).all() and np.isfinite(clipped[:, start:stop]).all() and
                    np.array_equal(clipped[:, start:stop], np.maximum(raw[:, start:stop], 0)), 'Persisted prediction projection differs')
            negative += int((raw[:, start:stop] < 0).sum())
            losses[:, arm_index] += np.square(target[:, start:stop] - clipped[:, start:stop]).sum(axis=1)
        close_map(raw)
        close_map(clipped)
        projected_counts.append(negative)
    losses /= R
    require(np.isfinite(losses).all(), 'Nonfinite donor MSE')
    boot = bootstrap_means(losses, [np.array(f['held_indices']) for f in folds], draws)
    contrasts = ((2, 1, 'RNA_PC_minus_technical'), (2, 0, 'RNA_PC_minus_mean'), (2, 3, 'RNA_PC_minus_source_fractions'))
    rows = []
    for i, donor in enumerate(donors):
        row = dict(donor_index=i, donor_id=donor['donor_id'], well_id=donor['well_id'])
        row.update({arm + '_MSE': float(losses[i, j]) for j, arm in enumerate(ARMS)})
        row.update({name: float(losses[i, a] - losses[i, b]) for a, b, name in contrasts})
        rows.append(row)
    write_table(directory / 'donor_losses.tsv', rows)
    well_rows = []
    for fold in folds:
        ix = fold['held_indices']
        row = dict(well_id=fold['well_id'], donors=len(ix))
        row.update({arm + '_MSE': float(losses[ix, j].mean()) for j, arm in enumerate(ARMS)})
        row.update({name: float((losses[ix, a] - losses[ix, b]).mean()) for a, b, name in contrasts})
        well_rows.append(row)
    write_table(directory / 'well_losses.tsv', well_rows)
    metrics = dict(donors=N, regions=R, wells=8, equal_donor_MSE=True,
                   absolute_MSE={arm: float(losses[:, i].mean()) for i, arm in enumerate(ARMS)},
                   raw_negative_entries_clipped=dict(zip(ARMS, projected_counts)), contrasts={},
                   bootstrap=dict(draws=10000, seed=20261001, bit_generator='PCG64',
                                  algorithm='Resample eight whole wells; divide summed donor losses by resampled donor count.',
                                  percentile_method='linear', descriptive_conditional_only=True, calibrated_p_value=None),
                   source_only_fraction_diagnostic=True, independent_validation=False,
                   necessary_descriptive_source_accuracy_guards_passed=bool(
                       losses[:, 2].mean() < losses[:, 1].mean() and losses[:, 2].mean() < losses[:, 0].mean()),
                   accuracy_guards_establish_independent_validation_or_practical_utility=False,
                   unconditional_training_uncertainty_included=False)
    for a, b, name in contrasts:
        values = boot[:, a] - boot[:, b]
        save_array(directory / (name + '_bootstrap.npy'), values)
        metrics['contrasts'][name] = dict(donor_equal_mean=float((losses[:, a] - losses[:, b]).mean()),
                                        descriptive_conditional_95_percentile_interval=np.percentile(values, [2.5, 97.5], method='linear').tolist())
    write_json(directory / 'metrics.json', metrics)
    cap(output)
    return metrics


def environment_inventory():
    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        np.show_config()
    paths = {Path(sys.executable).resolve()}
    for name, module in tuple(sys.modules.items()):
        if name == 'numpy' or name.startswith('numpy.'):
            path = getattr(module, '__file__', None)
            if path and Path(path).is_file():
                paths.add(Path(path).resolve())
    # Include loaded BLAS/LAPACK shared libraries, not only Python wrappers.
    maps = Path('/proc/self/maps')
    if maps.exists():
        for line in maps.read_text().splitlines():
            path = line.split()[-1]
            if path.startswith('/') and any(k in path.lower() for k in ('openblas', 'lapack', 'mkl', 'numpy')) and Path(path).is_file():
                paths.add(Path(path).resolve())
    return dict(python=sys.version, executable=sys.executable, numpy_version=np.__version__, numpy_config=stream.getvalue(),
                numerical_file_hashes={str(p): sha(p) for p in sorted(paths)},
                threads={k: os.environ.get(k) for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')},
                SVD='float64 numpy.linalg.svd full_matrices=False; exact20 slots; first max-absolute loading sign >=0',
                tied_20th_singular_space='No universal bitwise portability claimed; no accuracy-based tie selection.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true')
    parser.add_argument('--expected-helper-sha')
    parser.add_argument('--predict-json')
    parser.add_argument('--bundle')
    parser.add_argument('--bundle-summary-sha')
    parser.add_argument('--bundle-artifacts-sha')
    parser.add_argument('--admission-json')
    parser.add_argument('--admission-sha')
    args = parser.parse_args()
    if args.inspect:
        print(json.dumps(dict(input=str(INPUT), donors=N, genes=G, windows=R, directions=DIRECTIONS, arms=ARMS,
                              fixed_PC_slots=PCS, alpha_sum_SSE=ALPHA, cap_bytes=CAP, numerical_execution=False)))
        return
    if args.predict_json:
        inference_main(args)
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdecimal() and os.environ.get('SLURM_CPUS_PER_TASK') == '1', 'Requires one-CPU compute job')
    require(all(os.environ.get(k) == '1' for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS')), 'Numerical thread cap differs')
    require(args.expected_helper_sha and sha(Path(__file__)) == args.expected_helper_sha, 'Reviewed helper source drift')
    output = REC / ('codex_gse296875_donor_global_model_' + os.environ['SLURM_JOB_ID'])
    output.mkdir(mode=0o750)
    global ACTIVE_OUTPUT, ACTIVE_CAP
    ACTIVE_OUTPUT, ACTIVE_CAP = output, CAP
    started = time.monotonic()
    try:
        archive = output / 'executed_sources'
        archive.mkdir()
        source_pins = {}
        for source in (Path(__file__), LAUNCHER):
            data = source.read_bytes()
            (archive / source.name).write_bytes(data)
            source_pins[str(source)] = hashlib.sha256(data).hexdigest()
        require(source_pins[str(Path(__file__))] == args.expected_helper_sha, 'Source changed before archive')
        source_stats = {}
        for path, expected in GUARDS.items():
            before = fingerprint(path)
            require(sha(path) == expected and fingerprint(path) == before, 'Immutable authority/completion metadata differs')
            source_stats[str(path)] = before
        recipe = json.loads(PRESPEC.read_text())
        decision = json.loads(DECISION.read_text())
        require(not recipe['authorization']['fitting_authorized'] and decision['source_fitting_authorized'] is True and
                decision['prespec_sha256'] == GUARDS[PRESPEC] and ROOT / decision['source_input_directory'] == INPUT and
                decision['source_input_summary_sha256'] == GUARDS[INPUT / 'summary.json'] and
                decision['source_input_artifacts_sha256'] == GUARDS[INPUT / 'artifacts.json'] and
                not decision['external_prediction_or_evaluation_authorized'], 'Separate source execution authority differs')
        summary = json.loads((INPUT / 'summary.json').read_text())
        require(summary['raw_donor_global_aggregation_completed'] and summary['donors'] == N and summary['RNA_genes'] == G and
                summary['ATAC_peaks'] == R and summary['sampled_source_buckets'] == 272 and summary['source_nuclei'] == 68398 and
                summary['persisted_per_donor_per_feature_full_equals_halves'] and summary['all_source_hashes_verified_before_after'] and
                summary['source_and_output_row_global_totals_reconciled'] and summary['no_imputation_or_source_bucket_selection'],
                'Incomplete or changed raw aggregation')
        artifacts = json.loads((INPUT / 'artifacts.json').read_text())
        for name in ('donors.tsv','lineage_nucleus_fractions.tsv','genes.tsv','native_peak_geometry.tsv'):
            path = INPUT / name
            before = fingerprint(path)
            require(before['bytes'] == artifacts[name]['bytes'] and sha(path) == artifacts[name]['sha256'] and fingerprint(path) == before,
                    'Immutable source identity axis differs')
            source_stats[str(path)] = before
        # Import numerical code only after compute/source/authority guards and own archive.
        global np
        import numpy as np
        checks = selfchecks()
        write_json(output / 'selfchecks.json', checks)
        donors, folds, fractions = source_metadata(recipe, artifacts, summary)
        arrays, xs, tech = check_raw_arrays(donors, artifacts)
        write_json(output / 'source_raw_conservation.json', dict(real_all_features_full_equals_halves=True,
                   real_all_donor_row_totals_match=True, raw_conservation_precedes_all_transforms=True,
                   all_full_and_half_inputs_and_autosomal_denominators_defined=True, nuclei=68398, donors=N))
        rng = np.random.Generator(np.random.PCG64(20261001))
        draws = rng.integers(0, 8, size=(10000, 8))
        save_array(output / 'bootstrap_well_indices.npy', draws)
        initial_environment = environment_inventory()
        write_json(output / 'manifest.json', dict(prespec_sha256=GUARDS[PRESPEC], execution_decision_sha256=GUARDS[DECISION],
                   source_pins=source_pins, aggregate_metadata_pins={str(p):v for p,v in GUARDS.items()}, source_stat_before=source_stats,
                   aggregate_arrays=arrays, ordered_donor_ids=[d['donor_id'] for d in donors], folds=folds,
                   environment=initial_environment, argv=sys.argv, resources=dict(cpus=1,memory_GiB=4,hours=90),
                   output_cap_bytes_inclusive=CAP, no_external_H3_genotype_clinical_reads=True,
                   operational_arms=list(ARMS[:3]), source_only_diagnostic_arm=ARMS[3],
                   fraction_order=list(LABELS[:6]), omitted_fraction=LABELS[6],
                   technical_order=['ln1p_RNA_UMIs','ln1p_nuclei','detected_gene_fraction'],
                   primary_complete39_donors=True, losses_after_fixed_nonnegative_clipping=True))
        for name in ('donors.tsv','genes.tsv','native_peak_geometry.tsv','lineage_nucleus_fractions.tsv'):
            shutil.copyfile(INPUT / name, output / name)
        results = {}
        for direction, input_partition, target_partition in DIRECTIONS:
            directory = output / direction
            directory.mkdir()
            fold_dir = directory / 'folds'
            fold_dir.mkdir()
            raw = np.load(INPUT / f'ATAC_donor_{target_partition}_counts_uint64.npy', mmap_mode='r', allow_pickle=False)
            target = np.array(raw, dtype=np.float64)
            close_map(raw)
            target *= (1e6 / np.array([int(d[target_partition + '_autosomal_fragment_records']) for d in donors]))[:, None]
            np.log1p(target, out=target)
            require(np.isfinite(target).all(), 'Nonfinite fixed complete target')
            oof, filled = {}, {}
            for arm in ARMS:
                oof[arm] = tuple(new_map(directory / (arm + '_OOF_' + kind + '.npy.partial'), (N, R)) for kind in ('raw','clipped'))
                for value in oof[arm]:
                    value[:] = np.nan
                filled[arm] = np.zeros(N, dtype=bool)
            for fold in folds:
                train, held = np.array(fold['train_indices']), np.array(fold['held_indices'])
                transform = fit_transform(xs[input_partition], tech[input_partition], fractions[input_partition], train)
                require(np.array_equal(transform['training_indices'], train) and not np.intersect1d(train, held).size,
                        'Held donor included in fitted transform')
                designs = apply_transform(xs[input_partition], tech[input_partition], fractions[input_partition], transform)
                destination = fold_dir / fold['well_id']
                destination.mkdir()
                save_transform(destination, transform)
                write_json(destination / 'membership.json', dict(**fold, input_partition=input_partition,
                           target_partition=target_partition, target_head_training_indices=train.tolist(),
                           transform_training_indices=train.tolist(), no_held_data_used_to_fit=True))
                fit_heads(output, destination, designs, target, train, held, oof, filled)
                write_json(destination / 'completed.json', dict(completed=True, all_four_heads=True,
                           held_donors=len(held), training_donors=len(train), direct_outputs=R, alpha_sum_SSE=ALPHA))
                cap(output)
                print(json.dumps(dict(direction=direction, well=fold['well_id'], completed_folds=len(list(fold_dir.glob('*/completed.json'))),
                                      folds_total=8, held_donors=len(held))), flush=True)
            require(all(values.all() for values in filled.values()), 'Incomplete donor OOF prediction coverage')
            for arm in ARMS:
                for kind, value in zip(('raw','clipped'), oof[arm]):
                    close_map(value)
                    (directory / (arm + '_OOF_' + kind + '.npy.partial')).rename(directory / (arm + '_OOF_' + kind + '.npy'))
            results[direction] = evaluate(output, directory, target, donors, folds, draws)
            if input_partition == 'full':
                deployment = output / 'full_source_deployment'
                deployment.mkdir()
                train = np.arange(N)
                transform = fit_transform(xs['full'], tech['full'], fractions['full'], train)
                save_transform(deployment, transform)
                designs = apply_transform(xs['full'], tech['full'], fractions['full'], transform)
                fit_heads(output, deployment, designs, target, train, np.zeros(0, dtype=int), None, None)
                write_json(deployment / 'interface.json', dict(operational_heads=list(ARMS[:3]),
                           source_only_labeled_fraction_diagnostic_head=ARMS[3], training_donors=[d['donor_id'] for d in donors],
                           RNA_input='Complete36601 native raw integer snRNA UMI genes plus pooled nucleus number and admitted assay provenance.',
                           RNA_normalization='ln1p(1e6*UMI/sum_complete36601_UMI)', target_units='ln1p fragments per million all autosomal unique fragments',
                           target_axis='All306706 frozen native windows', clip_all_heads_at_zero=True,
                           external_use_or_compatibility_authorized=False, practical_or_independently_validated_model_claimed=False))
            del target, oof
        write_json(output / 'deployment_reconstruction.json', check_deployment_reconstruction(output, xs, tech, fractions, donors))
        for record in arrays:
            path = Path(record['path'])
            require(fingerprint(path) == record['stat_before'] and sha(path) == record['sha256'] and fingerprint(path) == record['stat_before'],
                    'Input raw array changed during source fitting')
        for path, expected in GUARDS.items():
            require(fingerprint(path) == source_stats[str(path)] and sha(path) == expected, 'Frozen authority/metadata changed during fitting')
        for name in ('donors.tsv','lineage_nucleus_fractions.tsv','genes.tsv','native_peak_geometry.tsv'):
            path = INPUT / name
            require(fingerprint(path) == source_stats[str(path)] and sha(path) == artifacts[name]['sha256'], 'Input identity axis changed during fitting')
        final_environment = environment_inventory()
        require(all(final_environment['numerical_file_hashes'].get(p) == v for p,v in initial_environment['numerical_file_hashes'].items()),
                'Imported numerical source/binary drift')
        for path, expected in source_pins.items():
            require(sha(Path(path)) == expected, 'Reviewed source changed during execution')
        write_json(output / 'final_environment.json', final_environment)
        output_artifacts = {}
        for path in sorted(output.rglob('*')):
            if path.is_file():
                entry = dict(bytes=path.stat().st_size, sha256=sha(path))
                if path.suffix == '.npy':
                    entry['header'] = header(path)
                output_artifacts[str(path.relative_to(output))] = entry
        write_json(output / 'artifacts.json', output_artifacts)
        cap(output)
        write_json(output / 'summary.json', dict(source_donor_global_model_completed=True, source_development_only=True,
                   donors=N, wells=8, windows=R, genes=G, fold_directions=3, complete_folds=24, operational_arms=list(ARMS[:3]),
                   source_only_fraction_diagnostic_arm=ARMS[3], primary=results[DIRECTIONS[0][0]],
                   cross_half_diagnostics={name:results[name] for name,_,_ in DIRECTIONS[1:]},
                   persisted_raw_source_full_equals_halves_verified=True, all_held_donors_excluded_from_transforms_and_heads=True,
                   complete_full_source_deployment_parameters_saved=True, source_hashes_verified_before_after=True,
                   bootstrap_descriptive_conditional_only=True, calibrated_population_p_values=None,
                   independent_validation_or_practical_utility_claimed=False, external_evaluation_run=False,
                   interval_limits=recipe['uncertainty']['limits'], parallel_hypothesis_policy=recipe['uncertainty']['parallel_hypothesis_policy'],
                   tied_20th_singular_space_qualification=final_environment['tied_20th_singular_space'],
                   elapsed_seconds=time.monotonic()-started, output_cap_bytes_inclusive=CAP,
                   saved_bytes_before_summary=size(output), final_inclusive_cap_checked=True))
        cap(output)
        print(json.dumps(dict(source_donor_global_model_completed=True, saved_bytes=size(output))), flush=True)
    except BaseException:
        if not (output / 'failure.json').exists():
            write_json(output / 'failure.json', dict(completed=False, traceback=traceback.format_exc(), partials_preserved=True))
        raise


if __name__ == '__main__':
    main()
