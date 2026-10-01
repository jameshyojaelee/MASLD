"""Future compute-only, fixed descriptive comparison of paired native masks.

Metadata admission and a saved methods contract precede numeric beta access.
Never opens mixed-fold labels or imports p/q values. No inference/API calls.
"""
import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[4]
BASE = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
MANIFEST = BASE/'target_region_manifest_21995444/development_target_manifest_1mb.tsv.gz'
SCORING = BASE/'target_readout_development_21997286'
EXPORT = BASE/'target_readout_effect_export_21997874'
EFFECTS = EXPORT/'development_target_effects.tsv.gz'
PROTOCOL = Path(__file__).with_name('target_readout_development_protocol.json')
PROTOCOL_SHA = 'f1e1c8451ba4bca0d0b7dc0228c3cb383541b6e900439a6b90f5e528cda0505d'
MANIFEST_SHA = 'ef47cc5b60aa2aab542c686772a620575af3f90e1567eac7c17beb1ac2c63eef'
EFFECTS_SHA = 'dc05d98a67b980ba7c34a64799d1b76d39d8dc58ba64d93c7b485908d731a7ad'
SOURCE_SHA = '72bfc08cb05dc0ff8d4e69580b1d20484ea6f9f2a3bb45336c8336edea1e9916'
SCORER_SHA = '586ecee076917c52098263e3313c170da587fc835ab370141c552e57a9831ca0'
HELPER_SHA = '3b8ceaac1cfd1118a12d81795250a666083b8682d7f0bdc6a55af9cf302f70ba'
TRACKS = ['UBERON:0001114 ATAC-seq', 'UBERON:0001115 ATAC-seq', 'UBERON:0002107 ATAC-seq']
IDENTITY = ['key','chr','pos_hg38','ref','alt','peak_id','peak_start0','peak_end0','heldout_fold']
SCORE_META = ['key','peak_id','heldout_fold','geometry_stratum','peak_width_bp','sequence_sha256']
MODELS = ['zero','training_mean','center_ols','peak_ols','center_peak_ols']


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def metadata(path, columns):
    frame = pd.read_csv(path, sep='\t', usecols=columns, dtype=str, keep_default_na=False)
    require(not frame.key.duplicated().any(), 'Duplicate metadata keys: '+str(path))
    require(frame.heldout_fold.isin(['1','2','3','4']).all(), 'Forbidden development fold: '+str(path))
    return frame


def save_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False)+'\n')


def correlation(y, pred):
    if np.ptp(y) == 0 or np.ptp(pred) == 0:
        return None
    value = float(np.corrcoef(y, pred)[0,1])
    return value if np.isfinite(value) else None


def fit_ols(train_x, train_y, test_x, columns):
    coef, _, rank, singular = np.linalg.lstsq(train_x, train_y, rcond=None)
    require(rank == train_x.shape[1], 'Rank-deficient OLS: '+','.join(columns))
    require(np.isfinite(coef).all(), 'Nonfinite OLS coefficients.')
    pred = test_x @ coef
    require(np.isfinite(pred).all(), 'Nonfinite OLS predictions.')
    return pred, {'columns':columns, 'coefficients':coef.tolist(), 'rank':int(rank),
                  'singular_values':singular.tolist(), 'training_rows':len(train_y)}


def main():
    require(bool(os.environ.get('SLURM_JOB_ID')), 'Compute allocation required.')
    # Completion is necessary before opening the growing prediction table.
    summary = read_json(SCORING/'summary.json')
    receipt = read_json(SCORING/'receipt.json')
    admission = read_json(SCORING/'sequence_admission.json')
    export = read_json(EXPORT/'summary.json')
    protocol = read_json(PROTOCOL)
    require(sha(PROTOCOL) == PROTOCOL_SHA, 'Frozen methods contract changed.')
    require(summary['status'] == 'development_paired_readouts_complete', 'Scoring unfinished.')
    require([summary[k] for k in ['selected','admitted','scored','rejected']] == [25477,24430,24430,1047],
            'Fixed complete scoring population changed.')
    require(summary['outcomes_read'] is False and summary['accuracy_evaluated'] is False,
            'Prediction producer used outcomes.')
    require(receipt['run_mode'] == 'full_development' and receipt['selected_rows'] == 25477,
            'Wrong prediction campaign.')
    require(receipt['manifest_sha256'] == MANIFEST_SHA and receipt['script_sha256'] == SCORER_SHA
            and receipt['helper_sha256'] == HELPER_SHA, 'Prediction source identity changed.')
    require(receipt['tracks'] == TRACKS and receipt['aggregation'] == 'DIFF_LOG2_SUM'
            and receipt['requested_outputs'] == ['ATAC','DNASE'], 'Native readout definition changed.')
    require(receipt['invariant_checks'] == 'passed' and receipt['outcomes_read'] is False,
            'Prediction invariant checks missing.')
    require(export['status'] == 'development_only_effect_export_complete'
            and export['covered_rows'] == 25477 and export['source_rows'] == 25748
            and export['geometry_excluded_rows'] == 271, 'Wrong effect export.')
    require(export['input_sha256']['source'] == SOURCE_SHA
            and export['input_sha256']['manifest'] == MANIFEST_SHA
            and export['output_sha256'] == EFFECTS_SHA, 'Effect export source changed.')
    require(export['fold0_source_opened'] is False and export['p_or_q_values_imported'] is False
            and export['fits_performed'] is False and export['metrics_computed'] is False,
            'Effect export boundaries changed.')
    require(protocol['training_folds'] == [2,3,4] and protocol['held_development_fold'] == 1
            and protocol['bootstrap']['replicates'] == 2000 and protocol['bootstrap']['seed'] == 20260930
            and list(protocol['models']) == MODELS, 'Fixed methods contract changed.')
    require(sha(MANIFEST) == MANIFEST_SHA and sha(EFFECTS) == EFFECTS_SHA, 'Fixed input hash mismatch.')
    paths = {'manifest':MANIFEST, 'effects':EFFECTS, 'scores':SCORING/'scores.tsv',
             'scoring_summary':SCORING/'summary.json', 'scoring_receipt':SCORING/'receipt.json',
             'sequence_admission':SCORING/'sequence_admission.json', 'effect_summary':EXPORT/'summary.json',
             'protocol':PROTOCOL, 'scoring_executed_source':SCORING/'executed_source.py'}
    hashes = {name:sha(path) for name,path in paths.items()}
    require(hashes['scoring_executed_source'] == SCORER_SHA, 'Executed scorer source mismatch.')
    manifest = metadata(MANIFEST, IDENTITY+['variant_id','assembly','geometry_stratum','peak_width_bp',
                                           'window_start0','window_end0','target_readout_state','coordinate_state'])
    effects_meta = metadata(EFFECTS, IDENTITY+['lead_variant_id','block_1mb'])
    scores_meta = metadata(paths['scores'], SCORE_META)
    require(len(manifest) == 25477 and len(effects_meta) == 25477 and len(scores_meta) == 24430,
            'Unexpected metadata row counts.')
    repeats = summary['same_process_repeats']
    require(len(repeats) == 3 and [r['key'] for r in repeats] == scores_meta.key.iloc[:3].tolist(),
            'Fixed same-process repeats missing or changed.')
    tolerance = protocol['same_process_repeat_absolute_tolerance']
    require(all(np.isfinite(r[field]) and abs(r[field]) <= tolerance
                for r in repeats for field in ['center_delta','peak_delta']),
            'Same-process repeated readout differs beyond numerical tolerance.')
    require(set(manifest.key) == set(effects_meta.key), 'Effects/manifest key sets differ.')
    effect_order = effects_meta.key.to_numpy(copy=True)
    effect_folds = effects_meta.heldout_fold.to_numpy(copy=True)
    effects_meta = effects_meta.set_index('key').loc[manifest.key].reset_index()
    for name in IDENTITY:
        require(np.array_equal(manifest[name], effects_meta[name]), 'Effects/manifest mismatch: '+name)
    require(np.array_equal(manifest.variant_id, effects_meta.lead_variant_id), 'Source variant identity differs.')
    require(manifest.assembly.eq('GRCh38').all() and manifest.target_readout_state.eq('whole_peak_contained').all()
            and manifest.coordinate_state.eq('geometry_agrees').all(), 'Native peak/assembly guard failed.')
    start = pd.to_numeric(manifest.peak_start0, errors='raise').to_numpy(np.int64)
    end = pd.to_numeric(manifest.peak_end0, errors='raise').to_numpy(np.int64)
    pos = pd.to_numeric(manifest.pos_hg38, errors='raise').to_numpy(np.int64)
    window_start = pd.to_numeric(manifest.window_start0, errors='raise').to_numpy(np.int64)
    window_end = pd.to_numeric(manifest.window_end0, errors='raise').to_numpy(np.int64)
    require((end > start).all() and (start >= window_start).all() and (end <= window_end).all()
            and np.array_equal(window_start, pos-1-524288)
            and np.array_equal(window_end-window_start, np.full(len(manifest),1048576)), 'Peak/window coordinates differ.')
    require(manifest.ref.str.fullmatch('[ACGT]').all() and manifest.alt.str.fullmatch('[ACGT]').all()
            and manifest.ref.ne(manifest.alt).all(), 'Native REF/ALT invalid.')
    require(manifest.groupby('chr').heldout_fold.nunique().max() == 1, 'Chromosome crosses folds.')
    require(admission['admitted'] == 24430 and admission['replacements'] == 0
            and admission['outcomes_read'] is False, 'Sequence admission changed.')
    rejected = pd.DataFrame(admission['rejected'])
    require(len(rejected) == 1047 and not rejected.key.duplicated().any(), 'Invalid sequence exclusion keys.')
    require((~rejected.acgt_ok | ~rejected.reference_match).all(), 'Sequence exclusion has no failed guard.')
    require(set(scores_meta.key).isdisjoint(rejected.key)
            and set(scores_meta.key) | set(rejected.key) == set(manifest.key), 'Common/excluded partition differs.')
    expected = manifest.set_index('key').loc[scores_meta.key].reset_index()
    for name in ['key','peak_id','heldout_fold','geometry_stratum','peak_width_bp']:
        require(np.array_equal(expected[name], scores_meta[name]), 'Score/manifest mismatch: '+name)
    require(scores_meta.sequence_sha256.str.fullmatch('[0-9a-f]{64}').all(), 'Invalid sequence hash.')
    rejected_expected = manifest.set_index('key').loc[rejected.key]
    require(np.array_equal(rejected_expected.geometry_stratum, rejected.geometry_stratum),
            'Sequence exclusion stratum changed.')
    require(np.array_equal(pd.to_numeric(expected.peak_width_bp).to_numpy(),
                           pd.to_numeric(expected.peak_end0).to_numpy()-pd.to_numeric(expected.peak_start0).to_numpy()),
            'Peak width differs from native interval.')
    coverage = manifest[['key','heldout_fold','geometry_stratum']].copy()
    coverage['common_rows'] = coverage.key.isin(scores_meta.key).astype(int)
    coverage['sequence_excluded_rows'] = coverage.key.isin(rejected.key).astype(int)
    coverage['fixed_manifest_rows'] = 1
    coverage = coverage.groupby(['heldout_fold','geometry_stratum'], as_index=False)[
        ['fixed_manifest_rows','common_rows','sequence_excluded_rows']].sum()
    for name,path in paths.items():
        require(sha(path) == hashes[name], 'Input changed during metadata admission: '+name)
    out = Path(os.environ['CODEX_REC_OUTPUT'])/('target_readout_development_comparison_'+os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    # The complete methods and all input hashes are written BEFORE beta access.
    save_json(out/'protocol.json', {'methods':protocol, 'input_paths':{k:str(v) for k,v in paths.items()},
                                  'input_sha256':hashes, 'script_sha256':sha(__file__),
                                  'metadata_admission':'passed', 'numeric_beta_read':False})
    (out/'executed_source.py').write_bytes(Path(__file__).read_bytes())
    coverage.to_csv(out/'coverage_by_fold_stratum.tsv', sep='\t', index=False)
    rejected[['key','geometry_stratum','acgt_ok','reference_match']].to_csv(
        out/'sequence_excluded_metadata.tsv', sep='\t', index=False)
    # Import only allowed beta, preserving the source ALT orientation without flips.
    effect_strings = metadata(EFFECTS, ['key','heldout_fold','beta_alt'])
    require(np.array_equal(effect_strings.key, effect_order)
            and np.array_equal(effect_strings.heldout_fold, effect_folds), 'Effect pass identity/order changed.')
    require(set(effect_strings.key) == set(manifest.key), 'Effect pass keys changed.')
    labels = effect_strings.set_index('key').loc[expected.key]
    require(np.array_equal(labels.heldout_fold, expected.heldout_fold), 'Effect pass folds changed.')
    beta = pd.to_numeric(labels.beta_alt, errors='raise').to_numpy(float)
    score_columns = ['center_score','peak_score']+[f'{mask}_track{i}' for mask in ['center','peak'] for i in range(3)]
    values = pd.read_csv(paths['scores'], sep='\t', usecols=['key']+score_columns)
    require(np.array_equal(values.key, expected.key), 'Numeric score keys/order changed.')
    require(np.isfinite(beta).all() and np.isfinite(values[score_columns].to_numpy(float)).all(), 'Nonfinite numeric input.')
    for mask in ['center','peak']:
        require(np.allclose(values[mask+'_score'],values[[f'{mask}_track{i}' for i in range(3)]].mean(axis=1),
                            atol=1e-12,rtol=1e-12), 'Three-track mean changed: '+mask)
    for name,path in paths.items():
        require(sha(path) == hashes[name], 'Input changed during numeric read: '+name)
    fold = expected.heldout_fold.to_numpy(int)
    train, held = np.isin(fold,[2,3,4]), fold == 1
    require(train.any() and held.any() and (train | held).all(), 'Incomplete training/held partition.')
    center, peak = values.center_score.to_numpy(float), values.peak_score.to_numpy(float)
    preds = {'zero':np.zeros(held.sum()), 'training_mean':np.full(held.sum(), beta[train].mean())}
    coefficients = {'zero':{'constant':0.0}, 'training_mean':{'constant':float(beta[train].mean()),
                                                          'training_rows':int(train.sum())}}
    for name, features in [('center_ols',[center]),('peak_ols',[peak]),('center_peak_ols',[center,peak])]:
        x = np.column_stack([np.ones(len(beta)),*features])
        columns = ['intercept']+({'center_ols':['center'],'peak_ols':['peak'],
                                 'center_peak_ols':['center','peak']}[name])
        preds[name], coefficients[name] = fit_ols(x[train],beta[train],x[held],columns)
    y = beta[held]
    rows = expected.loc[held,IDENTITY+['geometry_stratum','peak_width_bp']].reset_index(drop=True)
    rows['center_score'], rows['peak_score'], rows['beta_alt'] = center[held], peak[held], y
    metrics = {}
    for name,pred in preds.items():
        rows[name+'_prediction'] = pred
        rows[name+'_squared_error'] = (y-pred)**2
        rows[name+'_absolute_error'] = np.abs(y-pred)
        metrics[name] = {'MSE':float(np.mean((y-pred)**2)), 'MAE':float(np.mean(np.abs(y-pred))),
                         'Pearson':correlation(y,pred)}
    grouped = rows.groupby('chr',sort=True)
    per_chr = grouped.size().rename('rows').to_frame()
    for name in MODELS:
        per_chr[name+'_SSE'] = grouped[name+'_squared_error'].sum()
        per_chr[name+'_SAE'] = grouped[name+'_absolute_error'].sum()
    require(len(per_chr) >= 2, 'Insufficient chromosome clusters for bootstrap.')
    rng = np.random.default_rng(20260930)
    draws = rng.integers(0,len(per_chr),size=(2000,len(per_chr)))
    draw_n = per_chr.rows.to_numpy()[draws].sum(axis=1)
    contrasts = []
    for contrast in protocol['paired_contrasts']:
        ref, candidate = contrast['reference'], contrast['candidate']
        require(ref in MODELS and candidate in MODELS, 'Unknown contrast model.')
        per_chr[ref+'_minus_'+candidate+'_SSE'] = per_chr[ref+'_SSE']-per_chr[candidate+'_SSE']
        ref_sse = per_chr[ref+'_SSE'].to_numpy()[draws].sum(axis=1)
        cand_sse = per_chr[candidate+'_SSE'].to_numpy()[draws].sum(axis=1)
        delta = (ref_sse-cand_sse)/draw_n
        ratio_valid = bool((ref_sse > 0).all())
        ratios = cand_sse/ref_sse if ratio_valid else None
        ratio = metrics[candidate]['MSE']/metrics[ref]['MSE'] if metrics[ref]['MSE'] > 0 else None
        contrasts.append({**contrast, 'MSE_difference_reference_minus_candidate':metrics[ref]['MSE']-metrics[candidate]['MSE'],
                          'MSE_difference_percentile95':np.quantile(delta,[.025,.975]).tolist(),
                          'candidate_reference_MSE_ratio':ratio,
                          'ratio_percentile95':np.quantile(ratios,[.025,.975]).tolist() if ratio_valid else None,
                          'relative_MSE_reduction':1-ratio if ratio is not None else None,
                          'relative_reduction_percentile95':np.quantile(1-ratios,[.025,.975]).tolist() if ratio_valid else None,
                          'bootstrap_zero_reference_SSE_draws':int((ref_sse == 0).sum())})
    rows.to_csv(out/'held_development_row_predictions.tsv.gz',sep='\t',index=False,
                compression={'method':'gzip','mtime':0})
    per_chr.reset_index().to_csv(out/'held_development_per_chromosome.tsv',sep='\t',index=False)
    save_json(out/'coefficients.json',coefficients)
    save_json(out/'summary.json',{'status':'descriptive_development_comparison_complete',
        'fixed_manifest_rows':25477,'common_rows':24430,'sequence_excluded_rows':1047,
        'training_rows':int(train.sum()),'held_development_rows':int(held.sum()),
        'held_chromosomes':list(per_chr.index),'held_metrics':metrics,'paired_contrasts':contrasts,
        'source_pregeometry_fold_counts':export['source_fold_counts'],
        'geometry_excluded_fold_counts':export['excluded_fold_counts'],
        'input_sha256':hashes,'beta_unit':protocol['target'], 'bootstrap':protocol['bootstrap'],
        'independent_validation':False,'p_or_q_values_imported':False,'fold0_opened':False,
        'margin_or_pass_criterion':None,'biological_target_claims':False,
        'python':sys.version,'numpy':np.__version__,'pandas':pd.__version__,
        'output_sha256':{path.name:sha(path) for path in sorted(out.iterdir()) if path.is_file()}})
    print(json.dumps({'status':'descriptive_development_comparison_complete','common_rows':24430,
                      'training_rows':int(train.sum()),'held_development_rows':int(held.sum())}),flush=True)


if __name__ == '__main__':
    main()
