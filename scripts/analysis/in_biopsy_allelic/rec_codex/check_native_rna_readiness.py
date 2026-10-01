"""Compute-only completion check of B1 plus the four fixed RNA batches.

Inspect receipts, byte hashes/sizes and NPZ headers/identity metadata only.
Never deserialize expression counts, numerical predictions, H3 or outcomes.
"""
import argparse
import ast
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import traceback
import zipfile


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
B1 = REC / 'codex_native_rna_pilot_21998663'
INDEX = REC / 'codex_native_rna_index_21998573'
ROSTER = REC / 'paired_h3k27ac_raw_roster_21997742/raw_files.tsv'
FIX = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/model-data-064-21079902/fixture/molecular'
DONORS = ('B1', 'B7', 'B8', 'B9', 'B10', 'B12', 'B15', 'B19', 'B21',
          'B22', 'B24', 'B26', 'B36', 'B38', 'B41', 'B46', 'B47')
GROUPS = {'21998814': DONORS[1:5], '21998815': DONORS[5:9],
          '21998816': DONORS[9:13], '21998817': DONORS[13:17]}
ARMS = {'global': 'rrr', 'local': 'rrr_offset_cis'}
RECIPE = {
    'codex_native_rna_pilot.py': '413b9617160e0f472682aecc21a629002fae0f54556e1244ffc3aac044e232c5',
    'codex_build_native_rna_index.py': '8d6e64a8d270fc7133bdb5a07e56a2963bb305762b417212d01193dba8627c40',
    'export_training_h3_target_transform.py': '3499fb927373a04c7fc7121bd89598709774deb38f652b26d5506275e09ee907',
    'codex_native_tximport_gene_counts.R': '8aec882b16a9202348f0232c9f7d4365dd0c237386b2f98cea6cfb0219e86582',
}
BATCH_SHA = 'a4980388792007611dac5cae3c5c094ebdb120835cab61dea1863dd199d64247'
BATCH_LAUNCHER_SHA = 'cad120f728cbc29b1b4452e38a3585f2013d5bb8c7b383736fe46d07b0993827'
SMALL = {
    B1 / 'summary.json': '97d80729969fb1f075d9a18c4134a8d9a94cd570d6f7643dbcfac6ef25c41b46',
    B1 / 'protocol.json': '2bc6fdf97604daece80259b79f0992c9b2ac81f1af94251ffa2fc5f3abcdc0f8',
    B1 / 'geneCounts.report.json': '42b1c7a563d7ee4c2cf36a088eb80536bdf585e001aedec8abce8b7c885c73a4',
    INDEX / 'summary.json': '3461a87af753436971a1340b04eeabdd9bd1ba87169b14869c219fef4a6214f8',
    INDEX / 'protocol.json': '1388892c71c7949bf667b8e0cdf756669a053126daa0c59af4fbc4f4203fec05',
    ROSTER: 'f91c00fb37387020a624a02a388ff00d440efc7fa5bdac5ef41da2a8b1ae7e78',
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def fingerprint(path):
    info = path.stat()
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    require(path.is_file() and path.stat().st_size <= 5_000_000,
            'Missing or oversized receipt: ' + str(path))
    return json.loads(path.read_text())


def static_setup():
    for path, expected in SMALL.items():
        require(sha(path) == expected, 'Fixed small receipt/roster changed: ' + str(path))
    for name, expected in RECIPE.items():
        require(sha(HERE / name) == expected, 'Successful B1 recipe source changed: ' + name)
    baseline = read_json(B1 / 'protocol.json')
    groups = {}
    for job, donors in GROUPS.items():
        folder = REC / ('codex_native_rna_batch_' + job)
        receipt = read_json(folder / 'protocol.json')
        require(tuple(receipt['donors']) == donors and receipt['slurm_job_id'] == job,
                'Actual batch donor/job identity differs')
        require(receipt['source_sha256'] == RECIPE and receipt['copied_helper_sha256'] == RECIPE
                and receipt['wrapper_sha256'] == BATCH_SHA
                and receipt['launcher_sha256'] == BATCH_LAUNCHER_SHA,
                'Actual executed batch recipe differs')
        require(receipt['input_sha256'] == baseline['input_sha256']
                and receipt['index_summary_sha256'] == baseline['index_summary_sha256']
                and receipt['native_index_metadata_check'] == baseline['native_index_metadata_check'],
                'Actual batch reference/model metadata differs from B1')
        require(receipt['salmon_version'] == baseline['salmon_version'] == 'salmon 1.10.3'
                and receipt['tximport_version'] == '1.34.0' and receipt['jsonlite_version'] == '2.0.0'
                and receipt['Rscript_sha256'] == '750e8efcde0d04740fa754e5dc3cf22b66f5c176c7898a73d18e5b7d9635ad02',
                'Native environment metadata differs')
        groups[job] = (folder, receipt)
    require(tuple(donor for donors in GROUPS.values() for donor in donors) == DONORS[1:],
            'Fixed remaining donor source order changed')
    return baseline, groups


def terminal_jobs():
    binary = shutil.which('sacct')
    require(binary is not None, 'SLURM accounting binary unavailable; terminal completion unverified')
    jobs = ('21998663', *GROUPS)
    command = [binary, '-X', '-j', ','.join(jobs), '--format=JobIDRaw,State,ExitCode',
               '--parsable2', '--noheader']
    response = subprocess.run(command, capture_output=True, text=True, timeout=60, check=True)
    observed = {}
    for line in response.stdout.splitlines():
        fields = [field.strip() for field in line.split('|')]
        if len(fields) == 3 and fields[0] in jobs:
            require(fields[0] not in observed, 'Duplicate scheduler accounting record')
            observed[fields[0]] = dict(state=fields[1], exit_code=fields[2])
    require(set(observed) == set(jobs), 'Missing producer accounting record; completion unverified')
    for job, value in observed.items():
        require(value == {'state': 'COMPLETED', 'exit_code': '0:0'},
                'Producer is not terminal successful: ' + job + ' ' + str(value))
    return dict(command=command, jobs=observed)


def npy_header(handle, size):
    require(handle.read(6) == b'\x93NUMPY', 'NPZ member is not NPY')
    version = tuple(handle.read(2))
    if version == (1, 0):
        raw = handle.read(2)
        require(len(raw) == 2, 'Truncated NPY length')
        length = struct.unpack('<H', raw)[0]
        encoding = 'latin1'
    else:
        require(version in ((2, 0), (3, 0)), 'Unsupported NPY header version')
        raw = handle.read(4)
        require(len(raw) == 4, 'Truncated NPY length')
        length = struct.unpack('<I', raw)[0]
        encoding = 'utf8' if version == (3, 0) else 'latin1'
    require(0 < length <= 65536, 'Invalid NPY header length')
    raw = handle.read(length)
    require(len(raw) == length, 'Truncated NPY header')
    header = ast.literal_eval(raw.decode(encoding).strip())
    require(set(header) == {'descr', 'fortran_order', 'shape'}
            and isinstance(header['descr'], str) and isinstance(header['fortran_order'], bool)
            and isinstance(header['shape'], tuple)
            and all(type(value) is int and value >= 0 for value in header['shape']), 'Invalid NPY schema')
    require(not re.search(r'O|V', header['descr']), 'Object/structured payload forbidden; no pickle allowed')
    header['version'] = version
    header['payload_offset'] = handle.tell()
    header['uncompressed_size'] = size
    return header


def identity_array(archive, name, expected_shape):
    info = archive.getinfo(name)
    require(info.file_size <= 50_000_000, 'Oversized identity metadata member')
    with archive.open(info) as handle:
        header = npy_header(handle, info.file_size)
        require(header['shape'] == expected_shape, 'Identity metadata shape differs: ' + name)
        match = re.fullmatch(r'([<>])U([1-9][0-9]*)', header['descr'])
        require(match is not None, 'Identity metadata must be nonobject Unicode: ' + name)
        width = int(match[2])
        require(width <= 1024, 'Oversized Unicode identity width')
        count = math.prod(expected_shape)
        expected_bytes = count * width * 4
        require(header['payload_offset'] + expected_bytes == info.file_size,
                'Identity NPY byte count differs')
        payload = handle.read(expected_bytes)
        require(len(payload) == expected_bytes and handle.read(1) == b'', 'Truncated/extra identity payload')
        codec = 'utf-32-le' if match[1] == '<' else 'utf-32-be'
        values = tuple(payload[i:i+width*4].decode(codec).rstrip('\x00')
                       for i in range(0, len(payload), width*4))
    return values, header


def prediction_metadata(path, donor, form, regions, expected_count=96460):
    with zipfile.ZipFile(path) as archive:
        expected = {'profile.npy', 'training_mean.npy', 'region_key.npy', 'sample_id.npy',
                    'form.npy', 'transport.npy'}
        require(set(archive.namelist()) == expected and len(archive.namelist()) == len(expected),
                'Native prediction NPZ members differ')
        headers = {}
        for name, shape in (('profile.npy', (1, expected_count)), ('training_mean.npy', (expected_count,))):
            info = archive.getinfo(name)
            with archive.open(info) as handle:
                header = npy_header(handle, info.file_size)
            require(header['shape'] == shape and header['descr'] in ('<f8', '>f8')
                    and header['payload_offset'] + math.prod(shape)*8 == info.file_size,
                    'Numerical prediction header/byte length differs')
            headers[name] = header
            # No numerical payload bytes are requested or deserialized here.
        actual, headers['region_key.npy'] = identity_array(archive, 'region_key.npy', (expected_count,))
        require(actual == regions and len(set(actual)) == expected_count, 'Prediction region identity/order differs')
        actual, headers['sample_id.npy'] = identity_array(archive, 'sample_id.npy', (1,))
        require(actual == (donor,), 'Prediction sample identity differs')
        actual, headers['form.npy'] = identity_array(archive, 'form.npy', ())
        require(actual == (form,), 'Prediction model form differs')
        actual, headers['transport.npy'] = identity_array(archive, 'transport.npy', ())
        require(actual == ('raw',), 'Prediction transport differs')
    return dict(headers=headers, regions=expected_count, donor=donor, form=form, transport='raw',
                numerical_values_deserialized=False, finite_values_rechecked=False)


def check_metadata_only():
    """Check completed B1 identity/header metadata and six tiny refusal cases."""
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Metadata archive checks require compute')
    output = REC / ('codex_native_rna_metadata_check_' + os.environ['SLURM_JOB_ID'])
    require(not output.exists(), 'Refusing metadata-check output overwrite')
    output.mkdir(exist_ok=False)
    (output / 'executed_source.py').write_bytes(Path(__file__).read_bytes())
    (output / 'executed_launcher.sbatch').write_bytes((HERE / 'run_check_native_rna_readiness.sbatch').read_bytes())
    phase = 'completed_B1_receipts'
    try:
        require(sha(B1 / 'summary.json') == SMALL[B1 / 'summary.json'], 'Frozen B1 summary changed')
        summary = read_json(B1 / 'summary.json')
        axis = FIX / 'h3k27ac_feature_axis.tsv'
        protocol = read_json(B1 / 'protocol.json')
        require(sha(B1 / 'protocol.json') == SMALL[B1 / 'protocol.json'], 'Frozen B1 protocol changed')
        require(sha(axis) == protocol['input_sha256'][str(axis)], 'Frozen region metadata changed')
        with axis.open(newline='') as handle:
            regions = tuple(row['opaque_source_feature_key'] for row in csv.DictReader(handle, delimiter='\t'))
        require(len(regions) == len(set(regions)) == 96460, 'B1 full region axis differs')
        actual = {}
        phase = 'completed_B1_prediction_headers_and_identities'
        for arm, form in ARMS.items():
            path = B1 / (arm + '_prediction.npz')
            require(sha(path) == summary['predictions'][arm]['sha256'], 'Completed B1 prediction bytes changed')
            actual[arm] = prediction_metadata(path, 'B1', form, regions)
        phase = 'tiny_synthetic_archive_refusals'

        def npy_bytes(descr, shape, payload):
            header = repr(dict(descr=descr, fortran_order=False, shape=shape)).encode('latin1')
            padding = (16 - ((10 + len(header) + 1) % 16)) % 16
            header += b' '*padding + b'\n'
            return b'\x93NUMPY\x01\x00' + struct.pack('<H', len(header)) + header + payload

        def unicode_bytes(values, shape):
            width = max(len(value) for value in values)
            payload = b''.join(value.ljust(width, '\x00').encode('utf-32-le') for value in values)
            return npy_bytes('<U' + str(width), shape, payload)

        miniature = ('chr1:1-10', 'chr1:20-30', 'chr2:1-5')
        valid = {'profile.npy': npy_bytes('<f8', (1, 3), b'\x00'*24),
                 'training_mean.npy': npy_bytes('<f8', (3,), b'\x00'*24),
                 'region_key.npy': unicode_bytes(miniature, (3,)),
                 'sample_id.npy': unicode_bytes(('S0',), (1,)),
                 'form.npy': unicode_bytes(('rrr',), ()),
                 'transport.npy': unicode_bytes(('raw',), ())}

        def write_archive(name, members):
            path = output / (name + '.npz')
            with zipfile.ZipFile(path, 'x', compression=zipfile.ZIP_DEFLATED) as archive:
                for member, data in members.items():
                    archive.writestr(member, data)
            return path

        positive = prediction_metadata(write_archive('valid_tiny', valid), 'S0', 'rrr', miniature,
                                       expected_count=3)
        cases = {}
        changed = valid.copy(); changed['profile.npy'] = npy_bytes('<f8', (1, 2), b'\x00'*16)
        cases['wrong_profile_shape'] = changed
        changed = valid.copy(); del changed['training_mean.npy']
        cases['missing_required_member'] = changed
        changed = valid.copy(); changed['region_key.npy'] = unicode_bytes(miniature[::-1], (3,))
        cases['wrong_region_order'] = changed
        changed = valid.copy(); changed['sample_id.npy'] = unicode_bytes(('S1',), (1,))
        cases['wrong_donor_identity'] = changed
        changed = valid.copy(); changed['sample_id.npy'] = npy_bytes('|O', (1,), b'')
        cases['object_dtype_forbidden'] = changed
        changed = valid.copy(); changed['region_key.npy'] = changed['region_key.npy'][:-4]
        cases['truncated_identity_payload'] = changed
        refused = {}
        for name, members in cases.items():
            path = write_archive(name, members)
            try:
                prediction_metadata(path, 'S0', 'rrr', miniature, expected_count=3)
            except ValueError as error:
                refused[name] = str(error)
            else:
                raise AssertionError('Invalid archive accepted: ' + name)
        result = dict(metadata_decoder_checks_passed=True, completed_B1_arms=actual,
                      tiny_positive_archive=positive, invalid_cases_refused=refused,
                      numerical_values_deserialized=False, raw_fastqs_rehashed=False,
                      expression_tables_read=False, receiving_H3_read=False,
                      all17_readiness_evaluated=False, downstream_action_invoked=False,
                      checker_sha256=sha(output / 'executed_source.py'), python=sys.version,
                      slurm_job_id=os.environ['SLURM_JOB_ID'])
        (output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    except BaseException as error:
        (output / 'failure.json').write_text(json.dumps(dict(metadata_decoder_checks_passed=False,
            phase=phase, exception_type=type(error).__name__, message=str(error),
            traceback=traceback.format_exc(), numerical_values_deserialized=False,
            receiving_H3_read=False, all17_readiness_evaluated=False), indent=2) + '\n')
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    modes = parser.add_mutually_exclusive_group()
    modes.add_argument('--inspect', action='store_true', help='small static source/JSON setup only; no scheduler, inventories or output hashes')
    modes.add_argument('--check-metadata', action='store_true', help='compute-only B1 metadata and tiny invalid-archive checks; no all17 readiness')
    args = parser.parse_args()
    if args.check_metadata:
        check_metadata_only()
        return
    if args.inspect:
        static_setup()
        print(json.dumps(dict(expected_donor_order=DONORS, producer_jobs=['21998663', *GROUPS],
                              expected_native_transcripts=227368, expected_modeled_genes=42163,
                              expected_prediction_shape=[1, 96460], readiness_evaluated=False,
                              scheduler_called=False, numerical_values_deserialized=False), indent=2))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Byte hashes and final inventories require compute')
    output = REC / ('codex_native_rna_readiness_' + os.environ['SLURM_JOB_ID'])
    require(not output.exists(), 'Refusing readiness output overwrite')
    output.mkdir(exist_ok=False)
    (output / 'executed_source.py').write_bytes(Path(__file__).read_bytes())
    (output / 'executed_launcher.sbatch').write_bytes((HERE / 'run_check_native_rna_readiness.sbatch').read_bytes())
    fingerprints, verified = {}, {}
    phase, donor = 'fixed_small_metadata', None

    def hash_file(path, expected=None):
        path = Path(path)
        require(path.is_file(), 'Expected file absent: ' + str(path))
        before = fingerprint(path)
        if str(path) not in verified:
            actual = sha(path)
            require(fingerprint(path) == before, 'File changed during byte hashing: ' + str(path))
            verified[str(path)] = actual
            fingerprints[str(path)] = before
        require(fingerprints[str(path)] == before, 'File changed after prior check: ' + str(path))
        if expected is not None:
            require(verified[str(path)] == expected, 'Actual file hash differs from saved receipt: ' + str(path))
        return verified[str(path)]

    def inventory(folder):
        files = {}
        for path in sorted(folder.rglob('*')):
            require(not path.is_symlink(), 'Unexpected symlink in receiving output')
            if path.is_file():
                files[str(path)] = dict(bytes=path.stat().st_size)
        return dict(bytes=sum(value['bytes'] for value in files.values()), files=files)

    try:
        baseline, groups = static_setup()
        phase = 'producer_terminal_completion'
        scheduler = terminal_jobs()
        phase = 'frozen_source_reference_hashes'
        fixed = {Path(path): expected for path, expected in baseline['input_sha256'].items()}
        fixed.update(SMALL)
        fixed.update({HERE / name: expected for name, expected in RECIPE.items()})
        fixed[HERE / 'codex_native_rna_batch.py'] = BATCH_SHA
        fixed[HERE / 'run_codex_native_rna_batch.sbatch'] = BATCH_LAUNCHER_SHA
        index_summary = read_json(INDEX / 'summary.json')
        index_protocol = read_json(INDEX / 'protocol.json')
        require(index_summary['index_built'] and index_summary['index_returncode'] == 0
                and index_summary['protocol_sha256'] == SMALL[INDEX / 'protocol.json'],
                'Successful piloted index receipt differs')
        # Use the piloted index's existing reference receipts/metadata, without
        # introducing another whole-genome or search-payload hashing workflow.
        require(index_protocol['transcripts'] == 227368 and index_protocol['decoys'] == 194
                and index_protocol['k'] == 31 and index_protocol['keep_duplicates']
                and not index_protocol['poly_a_clipping'], 'Piloted index recipe differs')
        fixed.update({INDEX / 'index' / name: expected for name, expected in index_summary['index_metadata_sha256'].items()})
        for path, expected in fixed.items():
            hash_file(path, expected)
        for folder, receipt in groups.values():
            for name, expected in RECIPE.items():
                hash_file(folder / 'executed_sources' / name, expected)
            hash_file(folder / 'executed_sources/codex_native_rna_batch.py', BATCH_SHA)
            hash_file(folder / 'executed_sources/run_codex_native_rna_batch.sbatch', BATCH_LAUNCHER_SHA)
            binary = Path(baseline['quant_command'][0]).with_name('Rscript')
            hash_file(binary, receipt['Rscript_sha256'])
        with ROSTER.open(newline='') as handle:
            rows = list(csv.DictReader(handle, delimiter='\t'))
        with (FIX / 'h3k27ac_feature_axis.tsv').open(newline='') as handle:
            regions = tuple(row['opaque_source_feature_key'] for row in csv.DictReader(handle, delimiter='\t'))
        require(len(regions) == len(set(regions)) == 96460, 'Frozen full region metadata differs')
        phase = 'terminal_batch_summaries'
        locations = {'B1': (B1, baseline, None)}
        group_results = {}
        for job, (folder, receipt) in groups.items():
            final = read_json(folder / 'summary.json')
            state = read_json(folder / 'batch_state.json')
            require(final == state and final['group_complete'] is True
                    and final['cohort_complete'] is False and tuple(final['donors']) == GROUPS[job]
                    and set(final['status']) == set(GROUPS[job]) and 'fatal_exception' not in final,
                    'Batch does not have matching terminal successful summaries')
            require(not final['receiving_H3_read'] and not final['model_fitted']
                    and not final['predictive_accuracy_evaluated'], 'Batch contains an unauthorized evaluation')
            for donor in GROUPS[job]:
                item = final['status'][donor]
                require(item['state'] == 'succeeded' and Path(item['output']).resolve() == (folder / donor).resolve(),
                        'Fixed donor does not have terminal successful output')
                hash_file(folder / donor / 'summary.json', item['summary_sha256'])
                locations[donor] = (folder / donor, receipt, receipt['donor_caps'][donor])
            group_results[job] = dict(directory=str(folder), saved_cap=receipt['group_cap_bytes'])
        require(tuple(locations) == DONORS and len(locations) == 17, 'Final donor population/order differs')
        results = []
        for donor in DONORS:
            phase = 'donor_native_axis_prediction_metadata'
            folder, receipt, cap = locations[donor]
            summary = read_json(folder / 'summary.json')
            protocol = read_json(folder / 'protocol.json')
            report = read_json(folder / 'geneCounts.report.json')
            hash_file(folder / 'protocol.json', summary['protocol_sha256'])
            require(summary['donor'] == protocol['donor'] == report['sample'] == donor
                    and summary['biological_n'] == 1 and summary['acquired_files'] == 16
                    and summary['native_quantification_completed'] and summary['frozen_predictions_completed'],
                    'Native donor completion metadata differs')
            require(summary['modeled_genes'] == 42163 and summary['native_genes'] == 60623
                    and summary['normalized_native_genes'] == 60623
                    and not summary['rows_fabricated'] and not summary['biological_absence_asserted']
                    and summary['denominator'] == 'all42163 modeled native gene estimated counts, unchanged',
                    'Full native/modeled-axis metadata differs')
            require(not summary['receiving_H3_read']
                    and not summary['model_fitted'] and not summary['predictive_accuracy_evaluated']
                    and not summary['independently_validated'], 'Donor metadata makes unauthorized validation claim')
            require(report['native_transcripts'] == 227368 and report['native_genes'] == 60623
                    and report['countsFromAbundance'] == 'no' and not report['ignoreTxVersion']
                    and not report['ignoreAfterBar'] and not report['inferential_replicates_imported']
                    and not report['rows_padded_or_removed'] and report['full_model_gene_extraction_delegated']
                    and not report['biological_absence_asserted'] and not report['private_PISCES_equivalence_asserted'],
                    'Exact native transcript/gene aggregation metadata differs')
            require(report['script_sha256'] == RECIPE['codex_native_tximport_gene_counts.R']
                    and report['versions']['tximport'] == '1.34.0' and report['versions']['jsonlite'] == '2.0.0',
                    'Native aggregation source/environment differs')
            if donor == 'B1':
                cap = protocol['raw_cap_and_saved_output_cap_bytes']
            rr = sorted((row for row in rows if row['donor'] == donor and row['assay'] == 'RNA'),
                        key=lambda row: (row['run_accession'], int(row['file_number'])))
            require(len(rr) == 16 and len({row['run_accession'] for row in rr}) == 8
                    and protocol['technical_runs'] == sorted({row['run_accession'] for row in rr})
                    and protocol['expected_download_bytes'] == sum(int(row['fastq_bytes']) for row in rr),
                    'Donor frozen RNA roster differs')
            acquisitions = read_json(folder / 'acquisition_receipts.json')
            require(len(acquisitions) == 16, 'Donor acquisition receipt incomplete')
            for row, acquired in zip(rr, acquisitions):
                require(all(acquired[key] == value for key, value in row.items())
                        and acquired['verified_bytes'] == int(row['fastq_bytes'])
                        and acquired['gzip_crc_verified'] is True, 'RNA acquisition receipt differs from frozen metadata')
                expected_path = folder / 'raw' / Path(row['fastq_url']).name
                require(Path(acquired['file']).resolve() == expected_path.resolve(), 'Raw artifact path escaped donor')
                require(expected_path.stat().st_size == int(row['fastq_bytes']), 'Raw saved size differs')
                require(re.fullmatch(r'[0-9a-f]{64}', acquired['sha256']) is not None,
                        'Malformed producer raw SHA256 receipt')
                # Producer receipt/size metadata only. Do not rehash FASTQs.
            require(Path(report['inputs']['quant']).resolve() == (folder / 'quant/quant.sf').resolve(),
                    'Native quantifier input path differs')
            hash_file(folder / 'quant/quant.sf', report['inputs']['quant_sha256'])
            mapping = Path(next(path for path in baseline['input_sha256'] if path.endswith('/transcript_gene_axis.tsv')))
            require(Path(report['inputs']['tx2gene']).resolve() == mapping.resolve()
                    and report['inputs']['tx2gene_sha256'] == baseline['input_sha256'][str(mapping)],
                    'Native aggregation mapping identity differs')
            require(protocol['input_sha256'] == baseline['input_sha256'], 'Donor frozen recipe input hashes differ')
            hash_file(Path(report['inputs']['tx2gene']), report['inputs']['tx2gene_sha256'])
            hash_file(folder / 'geneCounts.tsv', report['output_sha256']['gene_counts'])
            hash_file(folder / 'geneCounts.tsv.native_zero_transcripts.tsv', report['output_sha256']['native_zero_transcripts'])
            require((folder / 'modeled_geneCounts.tsv').is_file(), 'Modeled counts artifact absent')
            hash_file(folder / 'modeled_geneCounts.tsv')
            require(set(summary['predictions']) == set(ARMS), 'Prediction arms omitted or added')
            arms = {}
            for arm, form in ARMS.items():
                saved = summary['predictions'][arm]
                path = folder / (arm + '_prediction.npz')
                require(Path(saved['path']).resolve() == path.resolve() and saved['form'] == form
                        and saved['transport'] == 'raw', 'Frozen prediction path/form/transport differs')
                hash_file(path, saved['sha256'])
                arms[arm] = prediction_metadata(path, donor, form, regions)
            final_inventory = inventory(folder)
            require(final_inventory['bytes'] <= cap, 'Final donor output including summary exceeds declared cap')
            results.append(dict(donor=donor, directory=str(folder), final_bytes=final_inventory['bytes'],
                                saved_cap=cap, native_transcripts=227368, native_genes=60623,
                                modeled_genes=42163, arms=arms, inventory=final_inventory))
            print(json.dumps(dict(donor=donor, checked_donors=len(results), total_donors=17)), flush=True)
        phase = 'final_group_inventory_and_stability'
        for job, (folder, receipt) in groups.items():
            current = inventory(folder)
            require(current['bytes'] <= receipt['group_cap_bytes'],
                    'Final batch output including summaries exceeds declared cap')
            group_results[job]['final_inventory'] = current
        for path, previous in fingerprints.items():
            require(fingerprint(Path(path)) == previous, 'Checked artifact changed before completion: ' + path)
        scheduler_after = terminal_jobs()
        require(scheduler_after['jobs'] == scheduler['jobs'], 'Producer terminal states changed during check')
        result = dict(all17_rna_frozen_predictions_ready=True, fixed_donor_order=DONORS,
                      scheduler=scheduler_after, donors=results, groups=group_results,
                      verified_file_sha256=verified, producer_finite_value_admission_used=True,
                      raw_fastqs_rehashed=False,
                      numerical_values_deserialized=False, receiving_H3_read=False,
                      H3_acquisition_or_evaluation_invoked=False, model_fitted=False,
                      predictive_accuracy_evaluated=False, independently_validated=False,
                      algorithmic_zeros='Native zero estimates retained; no per-gene positivity rule or biological-absence claim',
                      index_limit='Frozen reference/metadata verified; original producer has no full search-payload hash record',
                      identity_limit='Fixed public donor-labelled axis; biological independence and private measurement equivalence not established',
                      checker_sha256=sha(output / 'executed_source.py'), python=sys.version,
                      slurm_job_id=os.environ['SLURM_JOB_ID'])
        (output / 'summary.json').write_text(json.dumps(result, indent=2) + '\n')
    except BaseException as error:
        failure = dict(all17_rna_frozen_predictions_ready=False, phase=phase, donor=donor,
                       exception_type=type(error).__name__, message=str(error),
                       traceback=traceback.format_exc(), verified_file_sha256=verified,
                       numerical_values_deserialized=False, receiving_H3_read=False,
                       H3_acquisition_or_evaluation_invoked=False, outputs_modified=False)
        (output / 'failure.json').write_text(json.dumps(failure, indent=2) + '\n')
        raise


if __name__ == '__main__':
    main()
