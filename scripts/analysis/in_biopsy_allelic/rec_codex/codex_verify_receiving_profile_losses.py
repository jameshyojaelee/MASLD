"""Independent fixed H3 loss arithmetic; complete-primary compute execution only.

No producer functions, raw receiving counts, descriptors, fitting or calibration.
Signed target/predictions are used unchanged. This does not verify assay truth.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import sys
import traceback
import zipfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
PRIMARY = REC/'codex_receiving_profile_evaluation_22026740'
PRIMARY_SOURCE = HERE/'codex_receiving_profile_evaluation.py'
PRIMARY_SOURCE_SHA = '7ae9d401789880bd8ca40cbdbe75b54d7e73f43d1008048f6f390ee6308c9895'
RECOVERY = REC/'codex_receiving_h3_reads_storage_recovery_22026732'
RECOVERY_SOURCE = HERE/'codex_receiving_h3_reads_storage_recovery.py'
RECOVERY_SOURCE_SHA = 'b71d1f1bdd080bc157d24b2d121d558048bea6d5b835c212f99c22e1cf8811c3'
ORIGINAL_PROTOCOL = REC/'codex_receiving_h3_reads_21998882/protocol.json'
ORIGINAL_PROTOCOL_SHA = 'd6394b1f4a938a1f5681b727287d4d7676e67af2f53c9329585edc62af19a609'
AXIS = ROOT/'Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture/molecular/h3k27ac_feature_axis.tsv'
AXIS_SHA = 'cbe35aeb188e531283a1bcd636ba14de4cfbf215dfc544b99e3ebb97f9831986'
TRAINING_MEAN = REC/'training_h3_target_transform_21997855/training_residual_mean.npy'
MEAN_SHA = 'a16179d29b2a45ccb374987796ab75393528c54fc38f255b30faa6af255b7fa3'
DONORS = ('B1','B7','B8','B9','B10','B12','B15','B19','B21','B22',
          'B24','B26','B36','B38','B41','B46','B47')
ARMS = ('global','local','fixed_training_mean','zero_secondary')
REGIONS, SEED, DRAWS, CHUNK = 96460, 20260930, 10000, 8192
CAP, ATOL, RTOL = 4*1024*1024, 1e-12, 1e-12


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda:handle.read(1<<20), b''):
            h.update(block)
    return h.hexdigest()


def stat_identity(path):
    s = path.stat()
    return dict(bytes=s.st_size,mtime_ns=s.st_mtime_ns,inode=s.st_ino,device=s.st_dev)


class Consumed:
    def __init__(self):
        self.receipts = {}

    def add(self, path, expected):
        require(re.fullmatch('[0-9a-f]{64}',expected or '') is not None,'Missing valid input hash')
        before = stat_identity(path)
        require(sha(path)==expected and stat_identity(path)==before,'Input changed/hash differs: '+str(path))
        record = dict(sha256=expected,stat=before)
        if path in self.receipts:
            require(self.receipts[path]==record,'Conflicting input receipt')
        self.receipts[path] = record

    def json(self, path):
        require(path in self.receipts and path.stat().st_size<=10_000_000,'JSON not admitted/bounded')
        return json.loads(path.read_text())

    def recheck(self):
        for path, record in self.receipts.items():
            require(stat_identity(path)==record['stat'] and sha(path)==record['sha256']
                    and stat_identity(path)==record['stat'],'Consumed input changed: '+str(path))


def prediction_paths():
    dirs = {'B1':REC/'codex_native_rna_pilot_21998663'}
    for job, donors in ((21998814,('B7','B8','B9','B10')),
                       (21998815,('B12','B15','B19','B21')),
                       (21998816,('B22','B24','B26','B36')),
                       (21998817,('B38','B41','B46','B47'))):
        for donor in donors:
            dirs[donor] = REC/('codex_native_rna_batch_'+str(job))/donor
    return {(d,a):dirs[d]/(a+'_prediction.npz') for d in DONORS for a in ('global','local')}


def admit(consumed, expected_summary_sha):
    require(not (PRIMARY/'failure.json').exists() and not (RECOVERY/'failure.json').exists(),
            'Require completed primary and recovery without failure receipts')
    consumed.add(PRIMARY/'summary.json',expected_summary_sha)
    summary = consumed.json(PRIMARY/'summary.json')
    require(summary.get('status')=='complete_fixed_transport_evaluation'
            and summary.get('biological_n')==17 and summary.get('regions')==REGIONS
            and summary.get('inputs_unchanged_after_evaluation') is True
            and summary.get('model_fitted') is False and summary.get('independently_validated') is False
            and summary.get('source_person_disjointness_verified') is False
            and summary.get('primary_family')=='one registered paired mean-error contrast',
            'Primary completion does not admit fixed seventeen-person comparison')
    outputs = summary['output_sha256']
    for name in ('protocol.json','fixed_B99_residual_target.npy','paired_profile_results.json',
                 'executed_sources/codex_receiving_profile_evaluation.py',
                 'executed_sources/codex_receiving_h3_reads_storage_recovery.py'):
        require(name in outputs,'Required allowed primary output missing')
        consumed.add(PRIMARY/name,outputs[name])
    require(outputs['protocol.json']==summary['protocol_sha256']
            and outputs['executed_sources/codex_receiving_profile_evaluation.py']==PRIMARY_SOURCE_SHA
            and outputs['executed_sources/codex_receiving_h3_reads_storage_recovery.py']==RECOVERY_SOURCE_SHA,
            'Primary executed source/protocol differs')
    protocol = consumed.json(PRIMARY/'protocol.json')
    require(tuple(protocol['donors'])==DONORS and protocol['biological_n']==17
            and protocol['regions']==REGIONS and protocol['bootstrap_draws']==DRAWS
            and protocol['seed']==SEED and protocol['models_refitted'] is False
            and protocol['receiving_target_selected'] is False and protocol['donors_filtered_by_diagnostics'] is False
            and protocol['source_person_disjointness']=='unverified'
            and protocol['primary']=='mean donor MSE(global rrr) minus mean donor MSE(local rrr_offset_cis)'
            and protocol['target']=='log2 modeled-region CPM minus fixed [1, native concentration descriptors] @ B99'
            and protocol['H3_processing_directory']==str(RECOVERY)
            and protocol['H3_recovery_source_sha256']==RECOVERY_SOURCE_SHA,
            'Primary scope, target or bootstrap contract differs')
    bound = protocol['input_sha256']
    require(bound.get(str(PRIMARY_SOURCE))==PRIMARY_SOURCE_SHA,'Primary live source binding missing')
    consumed.add(PRIMARY_SOURCE,PRIMARY_SOURCE_SHA)
    require(bound.get(str(RECOVERY_SOURCE))==RECOVERY_SOURCE_SHA,'New H3 processor source binding missing')
    consumed.add(RECOVERY_SOURCE,RECOVERY_SOURCE_SHA)
    for name in ('summary.json','protocol.json'):
        consumed.add(RECOVERY/name,bound[str(RECOVERY/name)])
    recovery_summary = consumed.json(RECOVERY/'summary.json')
    recovery_protocol = consumed.json(RECOVERY/'protocol.json')
    require(recovery_summary['all17_reads_processed'] is True
            and tuple(recovery_summary['completed'])==DONORS and recovery_summary['biological_n']==17
            and recovery_summary['predictive_accuracy_evaluated'] is False
            and recovery_summary['protocol_sha256']==bound[str(RECOVERY/'protocol.json')]
            and recovery_summary['recovery_source_sha256']==RECOVERY_SOURCE_SHA,'Receiving processing incomplete')
    require(recovery_protocol['RNA_predictions_frozen_before_H3'] is True
            and tuple(recovery_protocol['fixed_donors'])==DONORS
            and recovery_protocol['model_fit'] is False and recovery_protocol['predictive_accuracy_evaluated'] is False
            and recovery_protocol.get('recovery_source_sha256')==RECOVERY_SOURCE_SHA
            and recovery_protocol['input_sha256'].get(str(RECOVERY_SOURCE))==RECOVERY_SOURCE_SHA,
            'Recovery pre-H3 prediction boundary differs')
    require(bound.get(str(ORIGINAL_PROTOCOL))==ORIGINAL_PROTOCOL_SHA
            and recovery_protocol['reused_input_sha256'].get(str(ORIGINAL_PROTOCOL))==ORIGINAL_PROTOCOL_SHA,
            'Original pre-H3 protocol binding missing')
    consumed.add(ORIGINAL_PROTOCOL,ORIGINAL_PROTOCOL_SHA)
    original = consumed.json(ORIGINAL_PROTOCOL)
    require(original['RNA_predictions_frozen_before_H3'] is True and tuple(original['fixed_donors'])==DONORS
            and original['model_fit'] is False,'Original pre-H3 prediction boundary differs')
    paths = prediction_paths()
    for path in paths.values():
        expected = bound.get(str(path))
        require(expected==recovery_protocol['input_sha256'].get(str(path))
                and expected==original['input_sha256'].get(str(path)),
                'Prediction not identically frozen before original H3 acquisition')
        consumed.add(path,expected)
    require(bound.get(str(AXIS))==AXIS_SHA and bound.get(str(TRAINING_MEAN))==MEAN_SHA,
            'Fixed region/baseline identity not admitted by primary')
    consumed.add(AXIS,AXIS_SHA);consumed.add(TRAINING_MEAN,MEAN_SHA)
    consumed.recheck()
    # All consumed allowed output/source/prediction hashes precede array access.
    return summary, paths


def read_header(handle, np):
    version = np.lib.format.read_magic(handle)
    require(version in ((1,0),(2,0)),'Unsupported NPY header version')
    reader = np.lib.format.read_array_header_1_0 if version==(1,0) else np.lib.format.read_array_header_2_0
    shape, fortran, dtype = reader(handle)
    require(not dtype.hasobject,'Object array refused')
    return dict(version=list(version),shape=list(shape),fortran=bool(fortran),dtype=dtype.str)


def file_header(path, np):
    with path.open('rb') as handle:
        return read_header(handle,np)


def archive_headers(path, np):
    info = {}
    with zipfile.ZipFile(path) as archive:
        names = archive.namelist()
        require(len(names)==len(set(names)) and set(names)=={
            'profile.npy','region_key.npy','sample_id.npy','form.npy','transport.npy','training_mean.npy'},
            'Prediction NPZ member schema differs')
        for name in names:
            member = archive.getinfo(name)
            require(member.file_size<=100_000_000,'Prediction member exceeds declared bound')
            with archive.open(name) as handle:
                info[name] = dict(npy=read_header(handle,np),bytes=member.file_size,
                                 compressed_bytes=member.compress_size,crc32=member.CRC)
    return info


def squared_error_mean(target, prediction):
    parts = []
    for start in range(0,REGIONS,CHUNK):
        stop = min(start+CHUNK,REGIONS)
        errors = target[start:stop] if prediction is None else target[start:stop]-prediction[start:stop]
        parts.append(math.fsum(float(x)*float(x) for x in errors))
    return math.fsum(parts)/REGIONS


def check_close(actual, expected, label):
    require(math.isfinite(float(actual)) and math.isfinite(float(expected))
            and abs(float(actual)-float(expected))<=ATOL+RTOL*abs(float(expected)),
            'Saved arithmetic differs: '+label)


def linear_percentile(ordered, fraction):
    position = (len(ordered)-1)*fraction
    lower, upper = int(math.floor(position)), int(math.ceil(position))
    weight = position-lower
    return (1-weight)*ordered[lower]+weight*ordered[upper]


def verify(consumed, summary, paths):
    import numpy as np
    require(np.__version__==summary['numpy'],'NumPy RNG environment differs from primary')
    with AXIS.open(newline='') as handle:
        rows = list(csv.DictReader(handle,delimiter='\t'))
    keys = tuple(row['opaque_source_feature_key'] for row in rows)
    require(len(keys)==len(set(keys))==REGIONS
            and [int(r['h3k27ac_feature_index']) for r in rows]==list(range(REGIONS)),
            'Fixed metadata region order differs')
    target_path = PRIMARY/'fixed_B99_residual_target.npy'
    headers = {str(target_path):file_header(target_path,np),str(TRAINING_MEAN):file_header(TRAINING_MEAN,np)}
    target = np.load(target_path,mmap_mode='r',allow_pickle=False)
    mean = np.load(TRAINING_MEAN,mmap_mode='r',allow_pickle=False)
    require(target.shape==(17,REGIONS) and mean.shape==(REGIONS,) and target.dtype==mean.dtype==np.dtype('float64')
            and np.isfinite(target).all() and np.isfinite(mean).all(),'Saved target/baseline schema differs')
    losses = {arm:[] for arm in ARMS}
    for index, donor in enumerate(DONORS):
        observed = target[index]
        losses['fixed_training_mean'].append(squared_error_mean(observed,mean))
        losses['zero_secondary'].append(squared_error_mean(observed,None))
        for arm, form in (('global','rrr'),('local','rrr_offset_cis')):
            path = paths[(donor,arm)]
            headers[str(path)] = archive_headers(path,np)
            with np.load(path,allow_pickle=False) as archive:
                require(tuple(archive['sample_id'].astype(str))==(donor,)
                        and tuple(archive['region_key'].astype(str))==keys
                        and archive['form'].shape==archive['transport'].shape==()
                        and archive['form'].item()==form and archive['transport'].item()=='raw',
                        'Prediction donor/region/form/transport differs')
                profile, baseline = archive['profile'], archive['training_mean']
                require(profile.shape==(1,REGIONS) and profile.dtype==np.dtype('float64')
                        and np.isfinite(profile).all() and baseline.shape==(REGIONS,)
                        and baseline.dtype==np.dtype('float64') and np.array_equal(baseline,mean),
                        'Frozen profile/baseline numeric schema differs')
                losses[arm].append(squared_error_mean(observed,profile[0]))
    cohort = {arm:math.fsum(losses[arm])/17 for arm in ARMS}
    delta = [losses['global'][i]-losses['local'][i] for i in range(17)]
    point = math.fsum(delta)/17
    rng = np.random.Generator(np.random.PCG64(SEED))
    draws = rng.integers(0,17,size=(DRAWS,17),dtype=np.int64)
    bootstrap = []
    for draw in draws:
        multiplicities = [0]*17
        for donor in draw:
            multiplicities[int(donor)] += 1
        bootstrap.append(math.fsum(delta[i]*multiplicities[i] for i in range(17))/17)
    ordered = sorted(bootstrap)
    interval = [linear_percentile(ordered,.025),linear_percentile(ordered,.975)]
    saved = consumed.json(PRIMARY/'paired_profile_results.json')
    require(tuple(saved['donor_id'])==DONORS and saved['biological_n']==17 and saved['regions']==REGIONS
            and saved['fixed_region_axis_sha256']==AXIS_SHA and set(saved['donor_mse'])==set(ARMS)
            and set(saved['cohort_mse'])==set(ARMS) and saved['bootstrap']['draws']==DRAWS
            and saved['bootstrap']['seed']==SEED and saved['bootstrap']['unit']=='biological donor'
            and saved['bootstrap']['paired_indices_shared'] is True
            and saved['bootstrap']['method']=='percentile, linear quantiles'
            and saved['target_or_prediction_recalibrated'] is False and saved['model_selected'] is False
            and saved['independently_validated'] is False,'Saved loss/bootstrap scope differs')
    residuals = []
    for arm in ARMS:
        require(len(saved['donor_mse'][arm])==17,'Saved donor loss count differs')
        for i in range(17):
            check_close(losses[arm][i],saved['donor_mse'][arm][i],arm+'/'+DONORS[i])
            residuals.append(abs(losses[arm][i]-saved['donor_mse'][arm][i]))
        check_close(cohort[arm],saved['cohort_mse'][arm],arm+'/cohort')
    require(len(saved['donor_paired_error_reduction'])==17
            and len(saved['paired_error_reduction_percentile95'])==2,'Saved contrast dimensions differ')
    for i in range(17):
        check_close(delta[i],saved['donor_paired_error_reduction'][i],'paired/'+DONORS[i])
    check_close(point,saved['mean_paired_error_reduction'],'paired/cohort')
    for i in range(2):
        check_close(interval[i],saved['paired_error_reduction_percentile95'][i],'paired/percentile'+str(i))
    baseline_loss = cohort['fixed_training_mean']
    skills = {arm:None if baseline_loss==0 else 1-cohort[arm]/baseline_loss for arm in ('global','local')}
    require(set(saved['descriptive_skill_vs_fixed_training_mean'])=={'global','local'},
            'Existing descriptive skill fields differ')
    for arm, value in skills.items():
        expected = saved['descriptive_skill_vs_fixed_training_mean'][arm]
        if value is None:
            require(expected is None,'Zero baseline acquired a defined skill')
        else:
            check_close(value,expected,'existing_skill/'+arm)
    consumed.recheck()
    for path_text, before in headers.items():
        path = Path(path_text)
        after = archive_headers(path,np) if path.suffix=='.npz' else file_header(path,np)
        require(after==before,'Consumed array header changed: '+path_text)
    require(not (PRIMARY/'failure.json').exists() and not (RECOVERY/'failure.json').exists(),
            'Failure receipt appeared during verification')
    return dict(donors=list(DONORS),regions=REGIONS,donor_mse=losses,cohort_mse=cohort,
                donor_paired_error_reduction=delta,mean_paired_error_reduction=point,
                paired_error_reduction_percentile95=interval,max_donor_mse_absolute_difference=max(residuals),
                descriptive_skill_vs_fixed_training_mean=skills,
                bootstrap=dict(draws=DRAWS,seed=SEED,unit='biological donor',paired_indices_shared=True,
                    algorithm='PCG64 int64 donor draws; explicit multiplicities; manually sorted linear percentiles',
                    draw_bytes_sha256=hashlib.sha256(draws.tobytes(order='C')).hexdigest()),
                numerical_agreement_tolerance=dict(absolute=ATOL,relative=RTOL,
                    meaning='Floating-point reconstruction tolerance, not a scientific success margin'),
                headers=headers,environment=dict(python=sys.version,numpy=np.__version__,
                    slurm_job_id=os.environ['SLURM_JOB_ID']))


def write_json(path, value):
    data = (json.dumps(value,indent=2,sort_keys=True,allow_nan=False)+'\n').encode()
    require(len(data)<=CAP,'Receipt exceeds bounded saved-output allowance')
    with path.open('xb') as handle:
        handle.write(data)


def saved_bytes(output):
    total = sum(p.stat().st_size for p in output.rglob('*') if p.is_file())
    for suffix in ('out','err'):
        log = REC/('1596_h3_verify_'+os.environ['SLURM_JOB_ID']+'.'+suffix)
        if log.exists():
            total += log.stat().st_size
    return total


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-source-sha256',required=True)
    parser.add_argument('--expected-primary-summary-sha256',required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID','').isdigit(),'Verification requires allocated compute')
    consumed = Consumed()
    consumed.add(Path(__file__),args.expected_source_sha256)
    launcher = HERE/'run_codex_verify_receiving_profile_losses.sbatch'
    consumed.add(launcher,sha(launcher))
    output = REC/('codex_verify_receiving_profile_losses_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    try:
        summary, paths = admit(consumed,args.expected_primary_summary_sha256)
        result = verify(consumed,summary,paths)
        write_json(output/'loss_reconstruction.json',result)
        write_json(output/'consumed_inputs.json',{str(p):r for p,r in consumed.receipts.items()})
        require(saved_bytes(output)<=CAP,'Saved verification bytes exceed four MiB')
        consumed.recheck()
        write_json(output/'summary.json',dict(status='saved_primary_loss_arithmetic_verified',primary_job=22026740,
            biological_n=17,regions=REGIONS,all_four_losses_verified=True,paired_point_and_interval_verified=True,
            inputs_and_headers_unchanged=True,source_arithmetic_verified_only=True,independently_validated=False,
            person_disjointness_verified=False,assay_truth_verified=False,B99_target_construction_reverified=False,
            raw_counts_or_descriptors_read=False,producer_functions_imported=False,model_fitted=False,
            calibration_or_selection_run=False,no_new_endpoints_or_pvalues=True,
            expected_primary_summary_sha256=args.expected_primary_summary_sha256,
            verifier_source_sha256=args.expected_source_sha256,
            output_sha256={p.name:sha(p) for p in output.iterdir() if p.is_file()},
            saved_bytes_before_summary=saved_bytes(output),cap_bytes=CAP,
            interpretation='Only arithmetic on saved B99 residual targets and pre-H3 frozen predictions; no verification of target construction, assay transport or participant independence'))
        require(saved_bytes(output)<=CAP,'Final saved verification bytes exceed four MiB')
        print(json.dumps(dict(status='saved_primary_loss_arithmetic_verified',biological_n=17,regions=REGIONS)),flush=True)
    except BaseException as exc:
        write_json(output/'failure.json',dict(error_type=type(exc).__name__,message=str(exc),
            traceback=traceback.format_exc(),outputs_preserved=True))
        raise


if __name__=='__main__':
    main()
