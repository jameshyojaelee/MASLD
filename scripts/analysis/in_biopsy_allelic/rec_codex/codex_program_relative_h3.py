#!/usr/bin/env python3
"""Compute-only fixed-program relative-H3 development experiment.

All numerical arrays are development GSE267145. No receiving data are opened.
Outputs and failures are retained in a new job directory; no head is selected.
"""
import ast
import csv
import hashlib
import json
import os
import platform
import shutil
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import sklearn
from sklearn.utils.extmath import randomized_svd

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
BENCH = ROOT / 'Analysis/MASLD_Model_Benchmark'
REC = BENCH / 'executions/codex-rec-20260929T142434Z'
FIX = BENCH / 'executions/model-data-064-21079902/fixture'
FIX_SHA = '612f0cb079ba611157663799c7cbb312aeedf8f32b25cba3ddaf171e1f9fbc29'
PROTOCOL_SHA = '5f14bd9d84bbc353c67378df8173602fd3bc7825cc72fc1a82175bfafd03342c'
SOURCES = {
    'membership': (ROOT / 'Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv', 'feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b'),
    'mapping': (BENCH / 'executions/gse267145-reference-crosswalk-21065154/rna_gene_crosswalk.tsv', '98ba04de29bf9220da3aa5b5c4a85d8ead4f50be6d848dfb879b1a2ab2005e57'),
    'baseline': (BENCH / 'scripts/fit_predict_gse267145_paired_baselines.py', '106892cd03de40281ad21ed49737750b5608935c56db3de276d65bddf5a55c4d'),
    'atlas_score': (ROOT / 'scripts/manuscript/program_observability_map/analysis_lib.R', 'a6439390506398be9283eb335b3866c45e9a510fa87aec07f1bf7fc21accb71b'),
    'task': (BENCH / 'config/evaluation/paired_bulk_rna_h3k27ac_task.toml', '8bf05fe693c787ffdd8fe11443e019e668ecd03318c708e390785edd837b1a28'),
}
CAP = 250_000_000
RUN_STATE = {'out': None, 'phase': 'startup'}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(reader.fieldnames and len(reader.fieldnames) == len(set(reader.fieldnames)), 'TSV header invalid')
        rows = list(reader)
    require(all(None not in row and all(v is not None for v in row.values()) for row in rows), 'TSV row width differs')
    return rows


def write_tsv(path, rows):
    require(bool(rows), 'cannot write empty result table')
    with path.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(rows)


def write_json(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
        handle.write('\n')


def check_cap(out):
    size = sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
    require(size <= CAP, 'aggregate output cap exceeded; retained outputs require review')
    return size


def source_functions(path):
    names = {'library_log1p', 'h3_hellinger', 'normalize_profile', 'top_variance_indices',
             'standardized_columns', 'rna_pca_scores', 'target_pca', 'ridge_predict'}
    nodes = [node for node in ast.parse(path.read_text()).body
             if isinstance(node, ast.FunctionDef) and node.name in names]
    require({node.name for node in nodes} == names, 'historical numeric functions missing')
    env = {'np': np, 'randomized_svd': randomized_svd, 'PairedBaselineFitError': ValueError,
           'LIBRARY_SCALE': 1e6, 'PROFILE_FLOOR': 1e-12,
           'RNA_VARIANCE_FEATURES': 4096, 'RNA_PCS': 20, 'H3_PCS': 30}
    exec(compile(ast.Module(body=nodes, type_ignores=[]), str(path), 'exec'), env)
    return env


def admission(members, mapping, axis):
    require(len(members) == 7093, 'frozen membership count differs')
    axis_ids = [r['stable_gene_id'] for r in axis]
    require(len(axis_ids) == 42163 and len(set(axis_ids)) == 42163, 'RNA stable axis differs')
    gene_index = {gene: i for i, gene in enumerate(axis_ids)}
    by_symbol = defaultdict(set)
    admitted_mapping = []
    for row in mapping:
        if row['stable_gene_id'] in gene_index:
            require(row['allowed_project_input'] == 'true' and row['mapping_state'] == 'stable_id_exact_v98_and_v49', 'modeled RNA mapping differs')
            admitted_mapping.append(row)
            by_symbol[row['gencode49_gene_name']].add(row['stable_gene_id'])
    require(len(admitted_mapping) == 42163 and len({r['stable_gene_id'] for r in admitted_mapping}) == 42163, 'mapping axis join not complete unique')
    require({r['matrix_gene_id'] for r in admitted_mapping} == {r['matrix_gene_id'] for r in axis}, 'matrix IDs differ from mapping')
    require({(r['matrix_gene_id'], r['stable_gene_id']) for r in admitted_mapping} == {(r['matrix_gene_id'], r['stable_gene_id']) for r in axis}, 'native matrix-to-stable pairs differ')
    by_program = defaultdict(list)
    flags = Counter()
    joined = []
    for member in members:
        weight = float(member['original_l1_weight'])
        require(np.isfinite(weight) and weight > 0, 'original biological weight invalid')
        symbol = member['mapped_symbol']
        state = ('mapping_unconfirmed' if member['mapped_symbol_status'] != 'gencode_v49_unique_symbol_confirmed'
                 else 'absent_fixed_axis' if not by_symbol[symbol]
                 else 'ambiguous_fixed_axis' if len(by_symbol[symbol]) != 1 else 'unique')
        flags[state] += 1
        gene = next(iter(by_symbol[symbol])) if state == 'unique' else ''
        row = dict(member, native_stable_gene_id=gene, state=state,
                   native_gene_index=gene_index[gene] if gene else '')
        joined.append(row)
        by_program[member['program_uid']].append(row)
    require(len(by_program) == 117 and flags == Counter(unique=6468, mapping_unconfirmed=611, absent_fixed_axis=14), 'fixed mapping census differs')
    coverage = []
    for uid in sorted(by_program):
        rows = by_program[uid]
        total = sum(float(r['original_l1_weight']) for r in rows)
        observed = sum(float(r['original_l1_weight']) for r in rows if r['state'] == 'unique')
        require(np.isclose(total, 1, rtol=0, atol=1e-12), 'original L1 total differs')
        coverage.append({'program_uid': uid, 'original_total_l1': total,
                         'observed_original_l1': observed, 'coverage': observed/total,
                         'observed_members': sum(r['state'] == 'unique' for r in rows),
                         'total_members': len(rows), 'admitted': int(observed/total >= .60)})
    uids = [r['program_uid'] for r in coverage if r['admitted']]
    require(len(uids) == 115, 'fixed115 program admission differs')
    pindex = {uid: i for i, uid in enumerate(uids)}
    weights = defaultdict(float)
    for row in joined:
        if row['program_uid'] in pindex and row['state'] == 'unique':
            weights[(pindex[row['program_uid']], int(row['native_gene_index']))] += float(row['original_l1_weight'])
    keys = sorted(weights)
    model_mapping = {'weight_program_index': np.array([k[0] for k in keys], dtype=np.int32),
                     'weight_gene_index': np.array([k[1] for k in keys], dtype=np.int32),
                     'original_weights': np.array([weights[k] for k in keys], dtype=np.float64),
                     'program_uid': np.array(uids), 'stable_gene_id': np.array(axis_ids)}
    return joined, coverage, model_mapping


def program_contribution(z, mapping):
    scores = np.zeros((len(z), 115), dtype=np.float64)
    for p in range(115):
        rows = mapping['weight_program_index'] == p
        scores[:, p] = z[:, mapping['weight_gene_index'][rows]] @ mapping['original_weights'][rows]
    require(np.isfinite(scores).all(), 'program contribution invalid')
    return scores


def program_inputs(training_log, query_log, mapping, numeric):
    train_z, query_z, gene_mean, gene_scale = numeric['standardized_columns'](training_log, query_log)
    train_q = program_contribution(train_z, mapping)
    query_q = program_contribution(query_z, mapping)
    train_x, query_x, q_mean, q_scale = numeric['standardized_columns'](train_q, query_q)
    require(train_x.shape[1] == 115, 'program coordinates dropped')
    return train_x, query_x, {'gene_mean': gene_mean, 'gene_scale': gene_scale,
                            'program_mean': q_mean, 'program_scale': q_scale}


def native_profile_floor(profile):
    p = np.asarray(profile, dtype=np.float64)
    require(p.shape[-1] == 96460 and np.isfinite(p).all() and (p >= 0).all(), 'baseline profile invalid')
    p = p + 1e-12
    return p / p.sum(axis=-1, keepdims=True)


def deviance(counts, profiles):
    p = np.broadcast_to(profiles, counts.shape)
    require(np.isfinite(p).all() and (p > 0).all() and np.allclose(p.sum(1), 1, rtol=0, atol=1e-12), 'deviance profile invalid')
    expected = counts.sum(1, dtype=np.float64)[:, None] * p
    terms = np.zeros(counts.shape, dtype=np.float64)
    nz = counts > 0
    terms[nz] = counts[nz] * np.log(counts[nz] / expected[nz])
    return 2 * terms.sum(1)


def predict_program(rna, stable_gene_ids, model):
    """Arithmetic callable, not an external validated inference claim."""
    require(np.array_equal(np.asarray(stable_gene_ids), model['stable_gene_id']), 'RNA identifier/order differs')
    x = np.asarray(rna, dtype=np.float64)
    require(x.ndim == 2 and x.shape[1] == 42163 and np.isfinite(x).all() and (x >= 0).all(), 'RNA estimate axis/units invalid')
    total = x.sum(1, keepdims=True)
    require((total > 0).all(), 'nonpositive RNA native-axis total')
    logx = np.log1p(x / total * 1e6)
    z = (logx - model['gene_mean']) / model['gene_scale']
    q = (program_contribution(z, model) - model['program_mean']) / model['program_scale']
    scores = (q - model['ridge_xmean']) @ model['ridge_coef'] + model['ridge_ymean']
    h = model['h3_mean'] + scores @ model['h3_loadings']
    p = np.maximum(h, 0)**2 + 1e-12
    require(np.isfinite(p).all() and (p.sum(1) > 0).all(), 'program profile invalid')
    return p / p.sum(1, keepdims=True)


def main():
    job = os.environ.get('SLURM_JOB_ID', '')
    require(job.isdigit(), 'compute-only: SLURM required')
    out = REC / ('codex_program_relative_h3_' + job)
    out.mkdir(exist_ok=False)
    RUN_STATE.update(out=out, phase='source_snapshots')
    protocol_path = Path(__file__).with_name('codex_program_relative_h3_protocol.json')
    launcher_path = Path(__file__).with_name('run_codex_program_relative_h3.sbatch')
    require(sha(protocol_path) == PROTOCOL_SHA, 'prospective protocol hash differs')
    protocol = json.loads(protocol_path.read_text())
    require(protocol['protocol_id'] == 'codex_program_relative_h3_development_v1' and protocol['aggregate_output_cap_bytes'] == CAP, 'protocol differs')
    source_receipt = {}
    for name, (path, expected) in SOURCES.items():
        require(sha(path) == expected, 'source differs: ' + name)
        shutil.copyfile(path, out / path.name)
        require(sha(out / path.name) == expected, 'source snapshot differs: ' + name)
        source_receipt[name] = {'path': str(path), 'sha256': expected}
    require(sha(FIX / 'ARTIFACTS.json') == FIX_SHA, 'fixture manifest differs')
    manifest = json.loads((FIX / 'ARTIFACTS.json').read_text())
    shutil.copyfile(protocol_path, out / 'prospective_protocol.json')
    own_receipt = {str(p): sha(p) for p in [Path(__file__), launcher_path, protocol_path]}
    shutil.copyfile(Path(__file__), out / Path(__file__).name)
    shutil.copyfile(launcher_path, out / launcher_path.name)
    for path, initial_sha in own_receipt.items():
        snapshot = out / ('prospective_protocol.json' if Path(path) == protocol_path else Path(path).name)
        require(sha(snapshot) == initial_sha, 'owned source snapshot differs: ' + path)
    write_json(out / 'source_receipt.json', {'sources': source_receipt,
               'fixture': str(FIX), 'fixture_manifest_sha256': FIX_SHA,
               'fixture_artifacts': manifest['artifacts'],
               'protocol_sha256': sha(protocol_path), 'helper_sha256': sha(Path(__file__)),
               'owned_initial_hashes': own_receipt,
               'environment': {'python': platform.python_version(), 'numpy': np.__version__, 'sklearn': sklearn.__version__},
               'opened_before_molecular_arrays': True, 'receiving_values': False})
    for item in manifest['artifacts']:
        require(sha(FIX / item['path']) == item['sha256'], 'fixture artifact differs: ' + item['path'])
    RUN_STATE['phase'] = 'metadata_admission'
    people = read_tsv(FIX / 'molecular/participant_axis.tsv')
    axis = read_tsv(FIX / 'molecular/rna_feature_axis.tsv')
    h3axis = read_tsv(FIX / 'molecular/h3k27ac_feature_axis.tsv')
    folds = read_tsv(FIX / 'folds/participant_outer_folds.tsv')
    require(len(people) == 99 and len({r['participant_id'] for r in people}) == 99, 'participant axis differs')
    require([int(r['participant_index']) for r in people] == list(range(99)), 'participant row order differs')
    require(all(r['pairing'] == 'same_sample_different_aliquot' and r['rna_observation_state'] == r['h3k27ac_observation_state'] == 'observed' for r in people), 'pair observation contract differs')
    require([r['participant_id'] for r in people] == [r['participant_id'] for r in folds], 'fold join differs')
    foldids = np.array([int(r['outer_fold']) for r in folds], dtype=np.int32)
    require([int((foldids == f).sum()) for f in range(5)] == [21, 21, 21, 19, 17], 'fold roster differs')
    require([int(r['rna_feature_index']) for r in axis] == list(range(42163)), 'RNA feature order differs')
    require([int(r['h3k27ac_feature_index']) for r in h3axis] == list(range(96460)) and len({r['opaque_source_feature_key'] for r in h3axis}) == 96460, 'H3 feature axis differs')
    joined, coverage, mapping = admission(read_tsv(SOURCES['membership'][0]), read_tsv(SOURCES['mapping'][0]), axis)
    write_tsv(out / 'program_membership_join.tsv', joined)
    write_tsv(out / 'program_coverage.tsv', coverage)
    for path in ['molecular/participant_axis.tsv', 'molecular/rna_feature_axis.tsv', 'molecular/h3k27ac_feature_axis.tsv', 'folds/participant_outer_folds.tsv']:
        shutil.copyfile(FIX / path, out / Path(path).name)
    for name in ['rna', 'h3k27ac']:
        mask = np.load(FIX / ('molecular/' + name + '_observed_mask.npy'), allow_pickle=False)
        require(mask.shape == (99,) and mask.dtype == np.bool_ and mask.all(), 'native observation mask differs')
    RUN_STATE['phase'] = 'development_array_validation'
    rna = np.load(FIX / 'molecular/rna_values.npy', allow_pickle=False)
    h3 = np.load(FIX / 'molecular/h3k27ac_counts.npy', allow_pickle=False)
    require(rna.shape == (99, 42163) and rna.dtype == np.float64 and np.isfinite(rna).all() and (rna >= 0).all() and (rna.sum(1) > 0).all(), 'RNA array invalid')
    require(h3.shape == (99, 96460) and h3.dtype == np.uint32 and (h3.sum(1) > 0).all(), 'H3 array invalid')
    numeric = source_functions(out / SOURCES['baseline'][0].name)
    log_rna = numeric['library_log1p'](rna)
    # Native row transforms use no other donor. Every learned parameter is fold-local.
    oof_program = np.lib.format.open_memmap(out / 'oof_program_profiles.npy', mode='w+', dtype=np.float64, shape=(99, 96460))
    oof_pca = np.lib.format.open_memmap(out / 'oof_pca_profiles.npy', mode='w+', dtype=np.float64, shape=(99, 96460))
    oof_program[:] = np.nan
    oof_pca[:] = np.nan
    oof_program.flush()
    oof_pca.flush()
    results = [None] * 99
    for fold in range(5):
        RUN_STATE['phase'] = 'development_fold_' + str(fold)
        train, query = np.flatnonzero(foldids != fold), np.flatnonzero(foldids == fold)
        x, qx, _ = program_inputs(log_rna[train], log_rna[query], mapping, numeric)
        hellinger, _ = numeric['h3_hellinger'](h3[train])
        target_scores, target_loadings, target_mean = numeric['target_pca'](hellinger, seed=267145+fold)
        require(target_loadings.shape == (30, 96460), 'target PCA rank differs')
        p_program = numeric['normalize_profile'](target_mean + numeric['ridge_predict'](x, qx, target_scores, alpha=10) @ target_loadings)
        px, pq, selected = numeric['rna_pca_scores'](log_rna[train], log_rna[query])
        require(px.shape[1] == 20 and len(selected) == 4096, 'historical RNA PCA recipe not executable')
        p_pca = numeric['normalize_profile'](target_mean + numeric['ridge_predict'](px, pq, target_scores, alpha=10) @ target_loadings)
        pooled = h3[train].sum(0, dtype=np.float64)
        pooled = native_profile_floor(pooled / pooled.sum())
        macro = native_profile_floor((hellinger**2).mean(0))
        d_program, d_pca = deviance(h3[query], p_program), deviance(h3[query], p_pca)
        d_pooled, d_macro = deviance(h3[query], pooled), deviance(h3[query], macro)
        require(np.isfinite(d_pooled).all() and (d_pooled > 0).all(), 'undefined pooled normalization; no donors may be dropped')
        require(np.isfinite(d_program).all() and np.isfinite(d_pca).all() and np.isfinite(d_macro).all(), 'nonfinite donor deviance')
        oof_program[query], oof_pca[query] = p_program, p_pca
        oof_program.flush()
        oof_pca.flush()
        fold_results = []
        for j, i in enumerate(query):
            row = {'participant_index': int(i), 'participant_id': people[i]['participant_id'], 'outer_fold': fold,
                   'h3_library_total_evaluator_only': int(h3[i].sum(dtype=np.uint64)),
                   'program_deviance': float(d_program[j]), 'pca_deviance': float(d_pca[j]),
                   'pooled_mle_deviance': float(d_pooled[j]), 'macro_mean_deviance': float(d_macro[j]),
                   'program_skill': float(1-d_program[j]/d_pooled[j]),
                   'pca_skill': float(1-d_pca[j]/d_pooled[j]),
                   'normalized_program_minus_pca_skill': float((d_pca[j]-d_program[j])/d_pooled[j])}
            results[i] = row
            fold_results.append(row)
        write_tsv(out / ('fold_' + str(fold) + '_donor_results.tsv'), fold_results)
        write_json(out / ('fold_' + str(fold) + '_receipt.json'), {'status': 'success', 'fold': fold,
                   'training_indices': train.tolist(), 'query_indices': query.tolist(),
                   'target_seed': 267145+fold, 'rna_pca_selected_indices': selected.tolist(),
                   'same_target_basis': True, 'no_selection': True,
                   'program_prediction_sha256': hashlib.sha256(p_program.tobytes()).hexdigest(),
                   'pca_prediction_sha256': hashlib.sha256(p_pca.tobytes()).hexdigest()})
        check_cap(out)
    require(all(r is not None for r in results) and np.isfinite(oof_program).all() and np.isfinite(oof_pca).all(), 'incomplete99 outputs; no all99 export allowed')
    require(np.allclose(oof_program.sum(1), 1, rtol=0, atol=1e-12) and np.allclose(oof_pca.sum(1), 1, rtol=0, atol=1e-12), 'OOF profiles not normalized')
    write_tsv(out / 'all99_donor_results.tsv', results)
    RUN_STATE['phase'] = 'conditional_development_summary'
    delta = np.array([r['normalized_program_minus_pca_skill'] for r in results])
    rng = np.random.default_rng(20260930)
    draws = delta[rng.integers(0, 99, size=(10000, 99))].mean(1)
    ci = np.quantile(draws, [.025, .975])
    skill = float(np.mean([r['program_skill'] for r in results]))
    development = {'participants': 99, 'primary_mean_skill_difference': float(delta.mean()),
                   'paired_conditional_percentile95': ci.tolist(), 'bootstrap_seed': 20260930,
                   'bootstrap_replicates': 10000, 'candidate_mean_pooled_skill': skill,
                   'pca_mean_pooled_skill': float(np.mean([r['pca_skill'] for r in results])),
                   'statistical_advancement_rule_passes': bool(ci[0] > 0 and skill > 0),
                   'external_or_program_specificity_claim': False,
                   'all99_export_regardless_of_result': True,
                   'predictor_dimensions': {'program': 115, 'rna_pca': 20},
                   'different_effective_complexity': True,
                   'skill_denominator': 'native fold pooled-count MLE with common additive profile floor; not historical macro-mean-baseline skill',
                   'interval_label': 'conditional percentile95 interval; descriptive development screening only',
                   'conditional_interval_limits': protocol['qualifications'],
                   'oof_program_sha256': sha(out / 'oof_program_profiles.npy'),
                   'oof_pca_sha256': sha(out / 'oof_pca_profiles.npy')}
    # Save the entire development result before any final all99 fit.
    write_json(out / 'development_results.json', development)
    RUN_STATE['phase'] = 'unconditional_all99_export'
    x, _, transforms = program_inputs(log_rna, log_rna, mapping, numeric)
    hellinger, _ = numeric['h3_hellinger'](h3)
    yp, loadings, mean = numeric['target_pca'](hellinger, seed=267145)
    xm, ym = x.mean(0), yp.mean(0)
    xc = x - xm
    coef = np.linalg.solve(xc.T@xc + 10*np.eye(115), xc.T@(yp-ym))
    model = dict(mapping, **transforms, ridge_xmean=xm, ridge_ymean=ym,
                 ridge_coef=coef, h3_loadings=loadings, h3_mean=mean)
    for key, value in model.items():
        if value.dtype.kind == 'f':
            require(value.dtype == np.float64 and np.isfinite(value).all(), 'final model precision/finiteness differs: ' + key)
    require(sum(v.nbytes for v in model.values()) < 100_000_000, 'factor storage exceeds100MB')
    np.savez(out / 'all99_program_relative_h3.npz', **model)
    exact = numeric['normalize_profile'](mean + numeric['ridge_predict'](x, x[:3], yp, alpha=10)@loadings)
    with np.load(out / 'all99_program_relative_h3.npz', allow_pickle=False) as archive:
        reproduced = predict_program(rna[:3], mapping['stable_gene_id'], archive)
    require(np.allclose(exact, reproduced, rtol=1e-10, atol=1e-12), 'all99 factor arithmetic reproduction differs; not an accuracy test')
    RUN_STATE['phase'] = 'postfit_source_and_output_verification'
    for name, (path, expected) in SOURCES.items():
        require(sha(path) == expected, 'source changed during experiment: ' + name)
    require(sha(protocol_path) == sha(out / 'prospective_protocol.json'), 'protocol changed during experiment')
    require(sha(Path(__file__)) == sha(out / Path(__file__).name), 'helper changed during experiment')
    for path, initial_sha in own_receipt.items():
        require(sha(Path(path)) == initial_sha, 'owned source changed during experiment: ' + path)
    require(sha(FIX / 'ARTIFACTS.json') == FIX_SHA, 'fixture manifest changed during experiment')
    for item in manifest['artifacts']:
        require(sha(FIX / item['path']) == item['sha256'], 'fixture changed during experiment')
    total_bytes = check_cap(out)
    write_json(out / 'summary.json', {'status': 'success', 'participants': 99, 'programs': 115,
               'final_model_sha256': sha(out / 'all99_program_relative_h3.npz'),
               'development_results_sha256': sha(out / 'development_results.json'),
               'output_bytes_before_summary': total_bytes, 'receiving_evaluated': False,
               'official_task_or_adoption_changed': False})
    (out / 'COMPLETE').write_text('success\n')
    RUN_STATE['phase'] = 'success'
    print('success:99 donor outputs and final115-program head retained; no receiving evaluation')


if __name__ == '__main__':
    try:
        main()
    except Exception as error:
        failed_out = RUN_STATE['out']
        if failed_out is not None:
            write_json(failed_out / 'failure.json', {'status': 'failure',
                       'phase': RUN_STATE['phase'], 'exception': type(error).__name__,
                       'message': str(error), 'partials_retained': True})
        raise
