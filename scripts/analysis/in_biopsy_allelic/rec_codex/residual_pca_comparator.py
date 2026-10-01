#!/usr/bin/env python3
"""Fixed ordinary residual-H3 comparator; compute only, no receiving H3."""
import ast
import csv
import hashlib
import json
import os
import platform
import shutil
from pathlib import Path

import numpy as np
import sklearn
from sklearn.utils.extmath import randomized_svd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BENCH = ROOT / 'Analysis/MASLD_Model_Benchmark'
REC = BENCH / 'executions/codex-rec-20260929T142434Z'
FIX = BENCH / 'executions/model-data-064-21079902/fixture'
TARGET = REC / 'training_h3_target_transform_21997855'
CAP = 500_000_000


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def table(path):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        rows = list(reader)
    require(rows and all(None not in r and all(v is not None for v in r.values()) for r in rows), 'TSV width/rows differ')
    return rows


def numeric(path, names, extra=None):
    nodes = [n for n in ast.parse(path.read_text()).body if isinstance(n, ast.FunctionDef) and n.name in names]
    require({n.name for n in nodes} == set(names), 'missing source function')
    env = dict(np=np, randomized_svd=randomized_svd, PairedBaselineFitError=ValueError,
               LIBRARY_SCALE=1e6, RNA_VARIANCE_FEATURES=4096, RNA_PCS=20, H3_PCS=30)
    env.update(extra or {})
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), env)
    return env


def fit_head(log_rna, residual, base, seed):
    indices = base['top_variance_indices'](log_rna, 4096)
    scaled, _, mean, scale = base['standardized_columns'](log_rna[:, indices], log_rna[:, indices])
    _, _, right = np.linalg.svd(scaled, full_matrices=False)
    rna_loadings = right[:20]
    require(len(indices) == 4096 and rna_loadings.shape == (20, 4096), 'RNA PCA dimensions differ')
    x = scaled @ rna_loadings.T
    y, loadings, target_mean = base['target_pca'](residual, seed=seed)
    require(loadings.shape == (30, 96460), 'residual PCA dimensions differ')
    xm, ym = x.mean(0), y.mean(0)
    xc = x - xm
    coef = np.linalg.solve(xc.T @ xc + 10 * np.eye(20), xc.T @ (y - ym))
    model = dict(rna_indices=indices, rna_mean=mean, rna_scale=scale,
                rna_loadings=rna_loadings, ridge_xmean=xm, ridge_ymean=ym,
                ridge_coef=coef, target_mean=target_mean, target_loadings=loadings)
    reference = target_mean + base['ridge_predict'](x, x[:3], y, alpha=10) @ loadings
    require(np.allclose(reference, predict_log(log_rna[:3], model), rtol=1e-10, atol=1e-12), 'source ridge factor arithmetic differs')
    return model


def predict_log(log_rna, model):
    x = (log_rna[:, model['rna_indices']] - model['rna_mean']) / model['rna_scale']
    x = x @ model['rna_loadings'].T
    y = (x - model['ridge_xmean']) @ model['ridge_coef'] + model['ridge_ymean']
    return model['target_mean'] + y @ model['target_loadings']


def predict_rna(rna, stable_ids, model):
    require(np.array_equal(np.asarray(stable_ids), model['stable_gene_id']), 'RNA feature identity/order differs')
    values = np.asarray(rna, dtype=np.float64)
    require(values.ndim == 2 and values.shape[1] == 42163 and np.isfinite(values).all()
            and (values >= 0).all() and (values.sum(1) > 0).all(), 'invalid complete native RNA')
    result = predict_log(np.log1p(values / values.sum(1, keepdims=True) * 1e6), model)
    require(result.shape == (len(values), 96460) and np.isfinite(result).all(), 'invalid residual predictions')
    return result


def cap(out):
    size = sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
    require(size <= CAP, 'output cap exceeded; partial outputs retained')
    return size


def run(out, protocol):
    guards = {ROOT / p: h for p, h in protocol['input_sha256'].items()}
    for path, expected in guards.items():
        require(sha(path) == expected, 'frozen input differs: ' + str(path))
    manifest = json.loads((FIX / 'ARTIFACTS.json').read_text())
    for item in manifest['artifacts']:
        guards[FIX / item['path']] = item['sha256']
    for path, expected in guards.items():
        require(sha(path) == expected, 'fixture differs: ' + str(path))
    source = BENCH / 'scripts/fit_predict_gse267145_paired_baselines.py'
    stable = BENCH / 'executions/chromatin-stable-rrr-forms-20260908T192024Z/transfer_functions.py'
    base = numeric(source, ['library_log1p', 'top_variance_indices', 'standardized_columns', 'target_pca', 'ridge_predict'])
    target = numeric(stable, ['conc_block', 'residualise'])
    people = table(FIX / 'molecular/participant_axis.tsv')
    axis = table(FIX / 'molecular/rna_feature_axis.tsv')
    regions = table(FIX / 'molecular/h3k27ac_feature_axis.tsv')
    folds = table(FIX / 'folds/participant_outer_folds.tsv')
    require(len(people) == 99 and len({r['participant_id'] for r in people}) == 99
            and [int(r['participant_index']) for r in people] == list(range(99)), 'participant axis differs')
    require(all(r['pairing'] == 'same_sample_different_aliquot' and r['rna_observation_state'] == r['h3k27ac_observation_state'] == 'observed' for r in people), 'pair admission differs')
    require([r['participant_id'] for r in people] == [r['participant_id'] for r in folds], 'fold join differs')
    ids = np.array([r['stable_gene_id'] for r in axis])
    keys = np.array([r['opaque_source_feature_key'] for r in regions])
    require(len(set(ids)) == 42163 and [int(r['rna_feature_index']) for r in axis] == list(range(42163)), 'RNA axis differs')
    require(len(set(keys)) == 96460 and [int(r['h3k27ac_feature_index']) for r in regions] == list(range(96460)), 'region axis differs')
    f = np.array([int(r['outer_fold']) for r in folds])
    require([int((f == i).sum()) for i in range(5)] == [21, 21, 21, 19, 17], 'fold allocation differs')
    for name in ['rna', 'h3k27ac']:
        mask = np.load(FIX / ('molecular/' + name + '_observed_mask.npy'), allow_pickle=False)
        require(mask.shape == (99,) and mask.dtype == np.bool_ and mask.all(), 'observability differs')
    rna = np.load(FIX / 'molecular/rna_values.npy', allow_pickle=False)
    counts = np.load(FIX / 'molecular/h3k27ac_counts.npy', allow_pickle=False)
    require(rna.shape == (99, 42163) and rna.dtype == np.float64, 'RNA matrix differs')
    require(counts.shape == (99, 96460) and counts.dtype == np.uint32 and (counts.sum(1) > 0).all(), 'H3 development counts differ')
    log_rna = base['library_log1p'](rna)
    raw = np.log2(1 + counts / counts.sum(1, keepdims=True) * 1e6)
    descriptors = target['conc_block'](counts.astype(np.float64))
    require(np.isfinite(descriptors).all(), 'invalid concentration descriptors')
    design = np.column_stack([np.ones(99), descriptors])
    b99 = np.load(TARGET / 'B99_concentration_coefficients.npy', allow_pickle=False)
    fixed_mean = np.load(TARGET / 'training_residual_mean.npy', allow_pickle=False)
    coefficient = np.linalg.lstsq(design, raw, rcond=None)[0]
    require(b99.shape == (5, 96460) and fixed_mean.shape == (96460,)
            and b99.dtype == fixed_mean.dtype == np.float64
            and np.isfinite(b99).all() and np.isfinite(fixed_mean).all(), 'fixed target dimensions/precision differ')
    require(np.allclose(coefficient, b99, rtol=0, atol=1e-9), 'all99 coefficients do not reproduce exact frozen B99')
    residual99 = raw - design @ b99
    source_residual = target['residualise'](raw, descriptors, np.arange(99))
    require(np.allclose(residual99, source_residual, rtol=0, atol=1e-9)
            and np.allclose(residual99.mean(0), fixed_mean, rtol=0, atol=1e-9), 'fixed target/mean equivalence failed')
    oof = np.lib.format.open_memmap(out / 'oof_residual_predictions.npy', mode='w+', dtype=np.float64, shape=(99, 96460))
    oof[:] = np.nan
    donor_rows = [None] * 99
    for fold in range(5):
        tr, te = np.flatnonzero(f != fold), np.flatnonzero(f == fold)
        residual = target['residualise'](raw, descriptors, tr)
        model = fit_head(log_rna[tr], residual[tr], base, 267145 + fold)
        pred = predict_log(log_rna[te], model)
        require(np.isfinite(pred).all(), 'nonfinite fold prediction')
        oof[te] = pred
        oof.flush()
        mse = np.mean((residual[te] - pred)**2, axis=1)
        mean_mse = np.mean((residual[te] - residual[tr].mean(0))**2, axis=1)
        for j, i in enumerate(te):
            donor_rows[i] = dict(participant_id=people[i]['participant_id'], outer_fold=fold,
                                 pca_mse=float(mse[j]), training_mean_mse=float(mean_mse[j]))
        write_json(out / ('fold_' + str(fold) + '.json'), dict(training_indices=tr.tolist(), query_indices=te.tolist(),
                   target_seed=267145 + fold, rna_selected_indices=model['rna_indices'].tolist(), B_fit_on_training_only=True))
        cap(out)
    require(all(r is not None for r in donor_rows) and np.isfinite(oof).all(), 'incomplete OOF results')
    write_json(out / 'development_donor_mse.json', donor_rows)
    delta = np.array([r['training_mean_mse'] - r['pca_mse'] for r in donor_rows])
    draws = np.random.default_rng(20260930).integers(0, 99, size=(10000, 99))
    write_json(out / 'development_summary.json', dict(mean_mse_reduction_vs_foldmean=float(delta.mean()),
               conditional_descriptive_percentile95=np.quantile(delta[draws].mean(1), [.025, .975]).tolist(),
               selection=False, export_unconditional=True, independently_validated=False))
    model = fit_head(log_rna, residual99, base, 267145)
    model.update(stable_gene_id=ids, region_key=keys, fixed_training_residual_mean=fixed_mean)
    require(all(np.isfinite(v).all() for v in model.values() if v.dtype.kind == 'f'), 'invalid model factors')
    np.savez(out / 'all99_residual_pca.npz', **model)
    reference = predict_log(log_rna[:3], model)
    require(np.allclose(reference, predict_rna(rna[:3], ids, model), rtol=1e-10, atol=1e-12), 'callable arithmetic differs')
    receiving = []
    for item in protocol['receiving_rna']:
        rows = table(ROOT / item['path'])
        require(list(rows[0]) == ['gene_id', item['donor']] and [r['gene_id'] for r in rows] == ids.tolist(), 'receiving RNA axis differs')
        receiving.append([float(r[item['donor']]) for r in rows])
    prediction = predict_rna(np.asarray(receiving), ids, model)
    np.save(out / 'frozen17_residual_pca_predictions.npy', prediction, allow_pickle=False)
    write_json(out / 'receiving_prediction_axis.json', dict(donors=[r['donor'] for r in protocol['receiving_rna']],
               region_axis_sha256=protocol['input_sha256'][str((FIX / 'molecular/h3k27ac_feature_axis.tsv').relative_to(ROOT))],
               receiving_H3_read=False, accuracy_evaluated=False, no_receiver_fit=True))
    for path, expected in guards.items():
        require(sha(path) == expected, 'frozen input changed during fit: ' + str(path))
    write_json(out / 'summary.json', dict(status='success', source_target_max_abs=float(np.max(np.abs(residual99 - source_residual))),
               coefficient_max_abs=float(np.max(np.abs(coefficient - b99))), output_bytes_before_summary=cap(out),
               model_sha256=sha(out / 'all99_residual_pca.npz'), prediction_sha256=sha(out / 'frozen17_residual_pca_predictions.npy'),
               environment=dict(python=platform.python_version(), numpy=np.__version__, sklearn=sklearn.__version__),
               receiving_accuracy_evaluated=False, existing_primary_changed=False))


def main():
    job = os.environ.get('SLURM_JOB_ID', '')
    require(job.isdigit(), 'SLURM compute required')
    out = REC / ('residual_pca_comparator_' + job)
    out.mkdir(exist_ok=False)
    protocol_path = HERE / 'residual_pca_comparator_protocol.json'
    protocol = json.loads(protocol_path.read_text())
    own = [Path(__file__), protocol_path, HERE / 'run_residual_pca_comparator.sbatch']
    initial = {p: sha(p) for p in own}
    for p in own:
        shutil.copyfile(p, out / p.name)
    write_json(out / 'execution_sources.json', {str(p): h for p, h in initial.items()})
    try:
        run(out, protocol)
        require(all(sha(p) == h for p, h in initial.items()), 'own source changed during run')
        cap(out)
        (out / 'COMPLETE').write_text('success\n')
    except Exception as exc:
        write_json(out / 'failure.json', dict(error_type=type(exc).__name__, message=str(exc), partial_outputs_retained=True))
        raise


if __name__ == '__main__':
    main()
