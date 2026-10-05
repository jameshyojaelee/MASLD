#!/usr/bin/env python3
"""Two descriptive fixed-PCA comparisons against saved primary residual target."""
import argparse
import csv
import hashlib
import json
import os
import shutil
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
PRIMARY = None  # Exact successful replacement primary job supplied prospectively.
PRIMARY_JOB = None
PRIMARY_SOURCE_SHA = '7ae9d401789880bd8ca40cbdbe75b54d7e73f43d1008048f6f390ee6308c9895'
H3_PROCESSING = REC / 'codex_receiving_h3_reads_storage_recovery_22026732'
H3_SOURCE_SHA = 'b71d1f1bdd080bc157d24b2d121d558048bea6d5b835c212f99c22e1cf8811c3'
DECLARATION_SHA = 'bf305fdfb60c8383787b7041193b41ff950ba351332d0f7093c25efb7f458b46'
PCA = REC / 'residual_pca_comparator_21999355'
CAP = 10_000_000


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''):
            h.update(b)
    return h.hexdigest()


def read_json(path):
    return json.loads(path.read_text())


def write_json(path, value):
    with path.open('x') as f:
        json.dump(value, f, indent=2, sort_keys=True, allow_nan=False)
        f.write('\n')


def run(out, declaration):
    require(PRIMARY is not None and PRIMARY_JOB is not None and not (PRIMARY/'failure.json').exists(),
            'Require successful exact replacement primary job')
    original_source_path = 'scripts/analysis/in_biopsy_allelic/rec_codex/codex_receiving_profile_evaluation.py'
    original_source_sha = 'e12c2e617faf1c35e3093ec715bb5f1af2c018eb8cba92869ce393c330efff33'
    require(declaration['primary_job']==21998897
            and declaration['primary_path']==str((REC/'codex_receiving_profile_evaluation_21998897').relative_to(ROOT))
            and declaration['input_sha256'][original_source_path]==original_source_sha,
            'Original descriptive declaration metadata differs')
    guards = {ROOT / p: h for p, h in declaration['input_sha256'].items()}
    # Only the cancelled primary job/source binding changes; fixed declaration stays on disk.
    guards[ROOT/original_source_path] = PRIMARY_SOURCE_SHA
    write_json(out/'primary_input_binding.json', dict(original_primary_job=21998897,
               replacement_primary_job=PRIMARY_JOB, replacement_primary_path=str(PRIMARY),
               original_primary_source_sha256=original_source_sha,
               replacement_primary_source_sha256=PRIMARY_SOURCE_SHA,
               H3_processing_directory=str(H3_PROCESSING), H3_recovery_source_sha256=H3_SOURCE_SHA,
               changed_metadata_fields=['primary_job','primary_path','primary source hash',
                                        'H3 processing path','H3 processing source hash'],
               endpoints_heads_donors_regions_seed_unchanged=True))
    for p, h in guards.items():
        require(sha(p) == h, 'frozen input differs: ' + str(p))
    # The future primary receipt is not assigned an invented presubmission hash.
    # Bind its complete outputs to the specified job, reviewed source and schema.
    primary_summary_sha = sha(PRIMARY/'summary.json')
    summary = read_json(PRIMARY / 'summary.json')
    require(summary['status'] == 'complete_fixed_transport_evaluation'
            and summary['biological_n'] == 17 and summary['regions'] == 96460
            and summary['inputs_unchanged_after_evaluation'] is True
            and summary['model_fitted'] is False
            and summary['independently_validated'] is False
            and summary['source_person_disjointness_verified'] is False
            and summary['primary_family'] == 'one registered paired mean-error contrast', 'primary receipt contract differs')
    for relative, expected in summary['output_sha256'].items():
        p = PRIMARY / relative
        require(p.resolve().is_relative_to(PRIMARY.resolve()), 'primary receipt path escapes directory')
        require(sha(p) == expected, 'primary output receipt hash differs: ' + relative)
        guards[p] = expected
    require({'fixed_B99_residual_target.npy', 'protocol.json', 'paired_profile_results.json'} <= set(summary['output_sha256']), 'primary output schema missing')
    require(sha(PRIMARY / 'protocol.json') == summary['protocol_sha256'], 'primary protocol hash differs')
    primary_protocol = read_json(PRIMARY / 'protocol.json')
    donors = declaration['donors']
    require(primary_protocol['donors'] == donors and primary_protocol['biological_n'] == 17
            and primary_protocol['regions'] == 96460
            and primary_protocol['models_refitted'] is False
            and primary_protocol['receiving_target_selected'] is False
            and primary_protocol['donors_filtered_by_diagnostics'] is False
            and primary_protocol['source_person_disjointness'] == 'unverified'
            and primary_protocol['target'] == 'log2 modeled-region CPM minus fixed [1, native concentration descriptors] @ B99'
            and primary_protocol['primary'] == 'mean donor MSE(global rrr) minus mean donor MSE(local rrr_offset_cis)', 'primary target/axis changed')
    require(primary_protocol.get('H3_processing_directory')==str(H3_PROCESSING)
            and primary_protocol.get('H3_recovery_source_sha256')==H3_SOURCE_SHA
            and not (H3_PROCESSING/'failure.json').exists(),
            'Primary receiving processing origin changed')
    processor=HERE/'codex_receiving_h3_reads_storage_recovery.py'
    processing_protocol=H3_PROCESSING/'protocol.json'
    processing_protocol_sha=primary_protocol['H3_protocol_sha256']
    require(primary_protocol['input_sha256'].get(str(processor))==H3_SOURCE_SHA
            and primary_protocol['input_sha256'].get(str(processing_protocol))==processing_protocol_sha
            and sha(processor)==H3_SOURCE_SHA and sha(processing_protocol)==processing_protocol_sha,
            'New H3 processor bytes or protocol not bound by primary')
    receiving_protocol=read_json(processing_protocol)
    require(receiving_protocol.get('recovery_source_sha256')==H3_SOURCE_SHA
            and receiving_protocol['input_sha256'].get(str(processor))==H3_SOURCE_SHA,
            'New H3 processor not bound by processing protocol')
    guards[processor]=H3_SOURCE_SHA
    guards[processing_protocol]=processing_protocol_sha
    for item in declaration['heads']:
        p = ROOT / item['path']
        require(primary_protocol['input_sha256'].get(str(p)) == item['sha256'], 'head not frozen in primary input receipt')
    for path in declaration['primary_source_paths']:
        p = ROOT / path
        require(primary_protocol['input_sha256'].get(str(p)) == guards[p], 'primary reviewed source/target identity differs')
    guards[PRIMARY / 'summary.json'] = primary_summary_sha
    axis = read_json(PCA / 'receiving_prediction_axis.json')
    require(axis['donors'] == donors and axis['region_axis_sha256'] == declaration['region_axis_sha256']
            and not axis['accuracy_evaluated'] and not axis['receiving_H3_read'] and axis['no_receiver_fit'], 'PCA axis/prediction contract differs')
    pca_summary = read_json(PCA / 'summary.json')
    require(pca_summary['status'] == 'success' and not pca_summary['receiving_accuracy_evaluated']
            and not pca_summary['existing_primary_changed']
            and pca_summary['prediction_sha256'] == declaration['pca_prediction_sha256'], 'PCA completion differs')
    with (ROOT / declaration['region_axis_path']).open() as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    regions = [r['opaque_source_feature_key'] for r in rows]
    require(len(regions) == len(set(regions)) == 96460
            and [int(r['h3k27ac_feature_index']) for r in rows] == list(range(96460)), 'region order differs')
    write_json(out / 'input_receipt.json', dict(frozen_input_sha256={str(p): h for p, h in guards.items()},
               primary_job=PRIMARY_JOB, primary_target_sha256=summary['output_sha256']['fixed_B99_residual_target.npy'],
               primary_protocol_sha256=summary['protocol_sha256'], accuracy_opened_after_receipt=True))
    for p,h in guards.items():
        require(sha(p)==h,'input changed before accuracy: '+str(p))
    target = np.load(PRIMARY / 'fixed_B99_residual_target.npy', allow_pickle=False)
    pca = np.load(PCA / 'frozen17_residual_pca_predictions.npy', allow_pickle=False)
    require(target.shape == pca.shape == (17, 96460) and target.dtype == pca.dtype == np.float64
            and np.isfinite(target).all() and np.isfinite(pca).all(), 'target/PCA numeric schema differs')
    predictions = {'global': [], 'local': []}
    for item in declaration['heads']:
        with np.load(ROOT / item['path'], allow_pickle=False) as a:
            require(set(a.files) == {'profile','region_key','sample_id','form','transport','training_mean'}, 'head archive schema differs')
            require(a['sample_id'].astype(str).tolist() == [item['donor']]
                    and a['region_key'].astype(str).tolist() == regions
                    and a['form'].item() == item['form'] and a['transport'].item() == 'raw', 'head donor/region/form differs')
            values = a['profile']
            require(values.shape == (1, 96460) and np.isfinite(values).all(), 'head profile invalid')
            predictions[item['arm']].append(values[0].copy())
    mse = {'pca': np.mean((target - pca)**2, axis=1)}
    for arm in ('global', 'local'):
        require([x['donor'] for x in declaration['heads'] if x['arm'] == arm] == donors, 'head donor sequence differs')
        mse[arm] = np.mean((target - np.stack(predictions[arm]))**2, axis=1)
    require(all(np.isfinite(v).all() for v in mse.values()), 'MSE overflow; no donor omission')
    primary_result = read_json(PRIMARY / 'paired_profile_results.json')
    require(primary_result['donor_id'] == donors and primary_result['biological_n'] == 17
            and primary_result['regions'] == 96460, 'primary result donor/region schema differs')
    for arm in ('global', 'local'):
        saved = np.asarray(primary_result['donor_mse'][arm], dtype=np.float64)
        require(saved.shape == (17,) and np.isfinite(saved).all()
                and np.allclose(mse[arm], saved, rtol=1e-10, atol=1e-12), 'independent MSE does not reproduce primary: ' + arm)
    require(np.isclose(float(np.mean(mse['global'] - mse['local'])),
                       primary_result['mean_paired_error_reduction'], rtol=1e-10, atol=1e-12),
            'independent donor averaging does not reproduce primary contrast')
    common = np.random.default_rng(20260930).integers(0, 17, size=(10000, 17))
    contrasts = {}
    for arm in ('global', 'local'):
        delta = mse['pca'] - mse[arm]
        boot = delta[common].mean(1)
        contrasts[arm] = dict(mean_pca_minus_existing_mse=float(delta.mean()), donor_difference=delta.tolist(),
                             descriptive_conditional_percentile95=np.quantile(boot, [.025, .975], method='linear').tolist())
    write_json(out / 'descriptive_results.json', dict(donors=donors, regions=96460,
               donor_mse={k: v.tolist() for k, v in mse.items()}, contrasts=contrasts,
               units='squared concentration-residual log2CPM averaged over all96460regions then equally over17donors',
               positive_contrast_means='existing head has lower MSE than PCA',
               bootstrap=dict(replicates=10000, seed=20260930, unit='participant', common_indices=True, method='linear percentile95',
                              limits='conditional on models/B99/counting/labels; omits source refits, technical error and assay shift; no calibrated coverage claim'),
               family='exactly two descriptive supplemental contrasts; unadjusted descriptive intervals, not confirmatory; no p-values, selection or primary change',
               source_person_disjointness='unverified', independently_validated=False,
               qualification='approximate native RNA; CUTRUN-trained residual target transported to single-read ChIP; not absolute physical mark abundance'))
    for p, h in guards.items():
        require(sha(p) == h, 'input changed during evaluation: ' + str(p))
    size = sum(p.stat().st_size for p in out.rglob('*') if p.is_file())
    require(size <= CAP, 'supplement cap exceeded; outputs retained')
    write_json(out / 'summary.json', dict(status='complete_descriptive_supplement', donors=17, regions=96460,
               result_sha256=sha(out / 'descriptive_results.json'), primary_target_sha256=summary['output_sha256']['fixed_B99_residual_target.npy'],
               primary_output_unchanged=True, primary_per_donor_mse_and_mean_contrast_reproduced=True,
               no_receiver_fit=True, no_pvalues_or_selection=True, independently_validated=False))


def main():
    global PRIMARY, PRIMARY_JOB
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--primary-job-id',type=int,required=True)
    parser.add_argument('--expected-primary-source-sha256',required=True)
    args=parser.parse_args()
    require(args.primary_job_id==22026740
            and args.expected_primary_source_sha256==PRIMARY_SOURCE_SHA,
            'Require replacement primary job and reviewed source')
    PRIMARY_JOB=args.primary_job_id
    PRIMARY=REC/('codex_receiving_profile_evaluation_'+str(PRIMARY_JOB))
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'compute-only SLURM required')
    out = REC / ('residual_pca_supplement_' + os.environ['SLURM_JOB_ID'])
    out.mkdir(exist_ok=False)
    paths = [Path(__file__), HERE / 'residual_pca_supplement_protocol.json', HERE / 'run_residual_pca_supplement.sbatch']
    require(sha(paths[1])==DECLARATION_SHA, 'Original fixed descriptive declaration changed')
    initial = {p: sha(p) for p in paths}
    for p in paths:
        shutil.copyfile(p, out / p.name)
    write_json(out / 'execution_sources.json', {str(p): h for p, h in initial.items()})
    try:
        run(out, read_json(paths[1]))
        require(all(sha(p) == h for p, h in initial.items()), 'supplement source changed')
        (out / 'COMPLETE').write_text('success\n')
    except BaseException as exc:
        write_json(out / 'failure.json', dict(error_type=type(exc).__name__, message=str(exc), outputs_retained=True))
        raise


if __name__ == '__main__':
    main()
