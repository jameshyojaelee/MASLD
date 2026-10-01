#!/usr/bin/env python3
"""Compute-only all99 extension of two historical relative-H3 recipes.

No receiving observations, residual targets, clinical fields or model selection.
Exports remain internal research derivatives subject to the source's terms.
"""
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

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
BENCH = ROOT / 'Analysis/MASLD_Model_Benchmark'
REC = BENCH / 'executions/codex-rec-20260929T142434Z'
FIXTURE = BENCH / 'executions/model-data-064-21079902/fixture'
SOURCES = {
    'c3': (BENCH / 'executions/codebase-review-20260905T010217Z/agents/C3/rna_conditioned_h3k27ac_review.py',
           '679bd8bea371dca1fdd423e041fe8592b19f9193fbb4f2b8f44ab9271d36ca5c'),
    'baseline': (BENCH / 'scripts/fit_predict_gse267145_paired_baselines.py',
                 '106892cd03de40281ad21ed49737750b5608935c56db3de276d65bddf5a55c4d'),
    'task': (BENCH / 'config/evaluation/paired_bulk_rna_h3k27ac_task.toml',
             '8bf05fe693c787ffdd8fe11443e019e668ecd03318c708e390785edd837b1a28'),
}
FIXTURE_SHA = '612f0cb079ba611157663799c7cbb312aeedf8f32b25cba3ddaf171e1f9fbc29'
FINAL_PCA_SEED = 267145
MODEL_CAP = 80_000_000


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def write_json(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def numeric_source(path, names, constants=None):
    """Execute only named numeric definitions; never source imports/load/main."""
    tree = ast.parse(path.read_text())
    nodes = [node for node in tree.body
             if isinstance(node, ast.FunctionDef) and node.name in names]
    require({node.name for node in nodes} == set(names), 'numeric source definitions missing')
    env = {'np': np, 'randomized_svd': randomized_svd,
           'PairedBaselineFitError': ValueError}
    env.update(constants or {})
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), env)
    return env


def table(path):
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def predict_relative(rna_values, archive, method):
    """RNA-only prediction; caller must separately verify fixed feature order.

    Input is complete native 42163 fractional-count axis. Denominator is its
    row sum, not an external exposure, TPM or total over a larger gene axis.
    """
    values = np.asarray(rna_values, dtype=np.float64)
    require(values.ndim == 2 and values.shape[1] == 42163,
            'RNA axis shape differs')
    require(np.isfinite(values).all() and (values >= 0).all(), 'invalid RNA values')
    totals = values.sum(1, keepdims=True)
    require((totals > 0).all(), 'nonpositive RNA library')
    cpm = values / totals * 1e6
    if method == 'c3_rank4':
        x = (np.log2(cpm + 1) - archive['c3_rna_mean']) / archive['c3_rna_scale']
        logv = (x @ archive['c3_rna_factor']) @ archive['c3_h3_factor'] + archive['c3_h3_mean']
        p = np.maximum(np.exp2(logv) - 1, 0)
    elif method == 'pca_ridge':
        x = np.log1p(cpm)[:, archive['pca_rna_indices']]
        x = (x - archive['pca_rna_mean']) / archive['pca_rna_scale']
        scores = x @ archive['pca_rna_loadings'].T
        scores = (scores - archive['pca_ridge_xmean']) @ archive['pca_ridge_coef'] + archive['pca_ridge_ymean']
        p = np.maximum(archive['pca_h3_mean'] + scores @ archive['pca_h3_loadings'], 0)**2
        p += 1e-12
    else:
        raise ValueError('unknown method')
    require(np.isfinite(p).all() and (p.sum(1) > 0).all(), 'invalid relative prediction')
    return p / p.sum(1, keepdims=True)


def main():
    job = os.environ.get('SLURM_JOB_ID', '')
    require(job.isdigit(), 'must run under SLURM')
    out = REC / ('codex_relative_h3_heads_' + job)
    out.mkdir(exist_ok=False)
    source_receipt = {}
    for name, (path, expected) in SOURCES.items():
        require(sha(path) == expected, 'source hash differs: ' + name)
        source_receipt[name] = {'path': str(path), 'sha256': expected}
        shutil.copyfile(path, out / path.name)
    require(sha(FIXTURE / 'ARTIFACTS.json') == FIXTURE_SHA, 'fixture manifest differs')
    manifest = json.loads((FIXTURE / 'ARTIFACTS.json').read_text())
    protocol = {
        'source': source_receipt, 'fixture_path': str(FIXTURE),
        'fixture_manifest_sha256': FIXTURE_SHA, 'fixture_artifacts': manifest['artifacts'],
        'fit_people': 99, 'rna_features': 42163, 'h3_regions': 96460,
        'candidate': 'historical C3 rank4, log2 CPM RNA and unresidualized log2 CPM H3; ddof0 RNA scale; ridge 0.001*Smax^2',
        'comparator': 'historical top4096 RNA ddof1 variance; ln1p CPM; ddof1 scale; RNA PC20; Hellinger H3 PC30 randomized SVD n_iter7; ridge alpha10',
        'final_pca_seed': FINAL_PCA_SEED,
        'seed_note': 'new fixed all99 seed; historical outer-fold seeds were 267145+fold',
        'baselines': ['native_pooled_count_mle', 'historical_training_mean_h3_profile', 'c3_inverse_mean_logcpm'],
        'normalization_axis': 'complete frozen modeled RNA row sum and fixed H3 region row sum, separately',
        'precision': 'float64 factors, no dense 42163x96460 coefficient',
        'model_byte_cap': MODEL_CAP,
        'interpretation': 'internal final all99 extensions; historical held-person development evidence does not evaluate these final artifacts',
        'receiving_target': 'not evaluated; separate prospective relative count-composition task and evaluator review required',
        'qualifications': 'no novelty, independent-person, private-PISCES equivalence or CUTRUN-to-ChIP equivalence claim; unchanged residual MSE task is separate',
        'receiving_family': 'unregistered here; candidate and baselines remain separate, no target-based chooser',
        'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'sklearn': sklearn.__version__},
        'helper_sha256': sha(Path(__file__)),
    }
    # Protocol is preserved before opening molecular arrays or fitting.
    write_json(out / 'protocol.json', protocol)
    for item in manifest['artifacts']:
        require(sha(FIXTURE / item['path']) == item['sha256'], 'fixture artifact hash differs: ' + item['path'])
    people = table(FIXTURE / 'molecular/participant_axis.tsv')
    rna_axis = table(FIXTURE / 'molecular/rna_feature_axis.tsv')
    h3_axis = table(FIXTURE / 'molecular/h3k27ac_feature_axis.tsv')
    folds = table(FIXTURE / 'folds/participant_outer_folds.tsv')
    require(len(people) == 99 and len({r['participant_id'] for r in people}) == 99, 'participant axis differs')
    require([int(r['participant_index']) for r in people] == list(range(99)), 'participant order differs')
    require(all(r['pairing'] == 'same_sample_different_aliquot' and r['rna_observation_state'] == r['h3k27ac_observation_state'] == 'observed' for r in people), 'pair observation contract differs')
    require([r['participant_id'] for r in people] == [r['participant_id'] for r in folds], 'fold order differs')
    require([sum(int(r['outer_fold']) == i for r in folds) for i in range(5)] == [21, 21, 21, 19, 17], 'fold counts differ')
    for rows, field, n in [(rna_axis, 'rna_feature_index', 42163), (h3_axis, 'h3k27ac_feature_index', 96460)]:
        require([int(r[field]) for r in rows] == list(range(n)), 'feature axis differs')
    require(len({r['stable_gene_id'] for r in rna_axis}) == 42163, 'RNA identifiers not unique')
    require(len({r['opaque_source_feature_key'] for r in h3_axis}) == 96460, 'H3 identifiers not unique')
    for name in ['rna', 'h3k27ac']:
        mask = np.load(FIXTURE / ('molecular/' + name + '_observed_mask.npy'), allow_pickle=False)
        require(mask.dtype == np.bool_ and mask.shape == (99,) and mask.all(), 'observation mask differs')
    rna = np.load(FIXTURE / 'molecular/rna_values.npy', allow_pickle=False)
    h3 = np.load(FIXTURE / 'molecular/h3k27ac_counts.npy', allow_pickle=False)
    require(rna.shape == (99, 42163) and rna.dtype == np.float64, 'RNA matrix contract differs')
    require(h3.shape == (99, 96460) and h3.dtype == np.uint32, 'H3 matrix contract differs')
    require(np.isfinite(rna).all() and (rna >= 0).all() and (rna.sum(1) > 0).all() and (h3.sum(1) > 0).all(), 'molecular libraries invalid')
    c3 = numeric_source(SOURCES['c3'][0], ['logcpm', 'to_props', 'rrr_fit_predict'])
    base = numeric_source(SOURCES['baseline'][0], ['library_log1p', 'h3_hellinger', 'normalize_profile', 'top_variance_indices', 'standardized_columns', 'rna_pca_scores', 'target_pca', 'ridge_predict'],
                          {'LIBRARY_SCALE': 1e6, 'PROFILE_FLOOR': 1e-12, 'RNA_VARIANCE_FEATURES': 4096, 'RNA_PCS': 20, 'H3_PCS': 30})
    x, y = c3['logcpm'](rna), c3['logcpm'](h3.astype(np.float64))
    xm, ym, xs = x.mean(0), y.mean(0), x.std(0)
    xs[xs == 0] = 1
    u, s, vt = np.linalg.svd((x - xm) / xs, full_matrices=False)
    lam = 1e-3 * s[0]**2
    d = s / (s**2 + lam)
    w = u.T @ (y - ym)
    _, sf, vft = np.linalg.svd((u * (s*d)) @ w, full_matrices=False)
    require(len(sf) >= 4, 'rank4 unavailable')
    vr = vft[:4].T
    result = {'c3_rna_mean': xm, 'c3_rna_scale': xs, 'c3_h3_mean': ym,
              'c3_rna_factor': (vt.T*d) @ (w@vr), 'c3_h3_factor': vr.T}
    # Exact-source comparison verifies exported arithmetic, not held accuracy.
    ref_c3 = c3['to_props'](c3['rrr_fit_predict'](x, y, x[:3], 4))
    require(np.allclose(predict_relative(rna[:3], result, 'c3_rank4'), ref_c3, rtol=1e-9, atol=1e-12), 'C3 factor reproduction differs')
    del x, y, u, s, vt, w, vr, vft
    xn = base['library_log1p'](rna)
    indices = base['top_variance_indices'](xn, 4096)
    scaled, _, mean, scale = base['standardized_columns'](xn[:, indices], xn[:, indices])
    _, _, right = np.linalg.svd(scaled, full_matrices=False)
    loadings = right[:20]
    xp = scaled @ loadings.T
    source_xp, _, source_indices = base['rna_pca_scores'](xn, xn[:3])
    require(np.array_equal(indices, source_indices) and np.allclose(xp, source_xp, rtol=1e-10, atol=1e-10), 'RNA PCA source reproduction differs')
    hellinger, totals = base['h3_hellinger'](h3)
    yp, target_loadings, target_mean = base['target_pca'](hellinger, seed=FINAL_PCA_SEED)
    xmean, ymean = xp.mean(0), yp.mean(0)
    centered = xp - xmean
    coef = np.linalg.solve(centered.T@centered + 10*np.eye(20), centered.T@(yp-ymean))
    result.update(pca_rna_indices=indices, pca_rna_mean=mean, pca_rna_scale=scale,
                  pca_rna_loadings=loadings, pca_ridge_xmean=xmean, pca_ridge_ymean=ymean,
                  pca_ridge_coef=coef, pca_h3_mean=target_mean, pca_h3_loadings=target_loadings)
    ref_pca = base['normalize_profile'](target_mean + base['ridge_predict'](xp, xp[:3], yp, alpha=10)@target_loadings)
    require(np.allclose(predict_relative(rna[:3], result, 'pca_ridge'), ref_pca, rtol=1e-10, atol=1e-12), 'PCA factor reproduction differs')
    pooled = h3.sum(0, dtype=np.float64)
    result['native_pooled_count_mle'] = pooled / pooled.sum()
    mean_profile = (hellinger**2).mean(0)
    mean_profile /= mean_profile.sum()
    result['historical_training_mean_h3_profile'] = base['normalize_profile'](np.sqrt(mean_profile)[None])[0]
    result['c3_inverse_mean_logcpm'] = c3['to_props'](ym[None])[0]
    for name in ['native_pooled_count_mle', 'historical_training_mean_h3_profile', 'c3_inverse_mean_logcpm']:
        require((result[name] >= 0).all() and np.isclose(result[name].sum(), 1, rtol=0, atol=1e-12), 'baseline profile invalid')
    require(all(v.dtype == np.float64 and np.isfinite(v).all() for k, v in result.items() if k != 'pca_rna_indices'), 'export precision or finiteness differs')
    require(sum(v.nbytes for v in result.values()) < MODEL_CAP, 'model byte budget exceeded')
    np.savez(out / 'relative_h3_heads.npz', **result)
    with np.load(out / 'relative_h3_heads.npz', allow_pickle=False) as archive:
        for method in ['c3_rank4', 'pca_ridge']:
            p = predict_relative(rna, archive, method)
            require(np.allclose(p.sum(1), 1, atol=1e-12) and (p >= 0).all(), 'invalid exported profile')
    for path in ['molecular/participant_axis.tsv', 'molecular/rna_feature_axis.tsv', 'molecular/h3k27ac_feature_axis.tsv', 'folds/participant_outer_folds.tsv']:
        shutil.copyfile(FIXTURE / path, out / Path(path).name)
    # Recheck fixture hashes after fitting, preserving partial outputs on failure.
    for item in manifest['artifacts']:
        require(sha(FIXTURE / item['path']) == item['sha256'], 'fixture changed during fit')
    write_json(out / 'summary.json', {'status': 'success', 'participants': 99,
               'rna_features': 42163, 'h3_regions': 96460, 'final_pca_seed': FINAL_PCA_SEED,
               'model_sha256': sha(out / 'relative_h3_heads.npz'),
               'model_bytes': (out / 'relative_h3_heads.npz').stat().st_size,
               'protocol_sha256': sha(out / 'protocol.json'), 'receiving_evaluated': False,
               'checks': 'hashes/axes/observability/source factor reproduction/float64/probability normalization passed'})
    (out / 'COMPLETE').write_text('success\n')
    print('success: all99 relative H3 rank4 and PCA-ridge exports; receiving evaluation not run')


if __name__ == '__main__':
    main()
