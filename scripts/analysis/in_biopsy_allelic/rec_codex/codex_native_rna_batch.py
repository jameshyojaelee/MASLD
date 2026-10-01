"""Four fixed serial RNA donors per allocation using the successful B1 recipe.

No receiving H3, model fitting, comparison, padding or donor substitution.
Every requested donor retains a terminal success/failure entry and outputs.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback
import types


HERE = Path(__file__).resolve().parent
ROOT = Path(__file__).resolve().parents[4]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
GROUPS = (('B7', 'B8', 'B9', 'B10'), ('B12', 'B15', 'B19', 'B21'),
          ('B22', 'B24', 'B26', 'B36'), ('B38', 'B41', 'B46', 'B47'))
PAIRED = {'B1', *(d for group in GROUPS for d in group)}
SOURCES = {
    'codex_native_rna_pilot.py': '413b9617160e0f472682aecc21a629002fae0f54556e1244ffc3aac044e232c5',
    'codex_build_native_rna_index.py': '8d6e64a8d270fc7133bdb5a07e56a2963bb305762b417212d01193dba8627c40',
    'export_training_h3_target_transform.py': '3499fb927373a04c7fc7121bd89598709774deb38f652b26d5506275e09ee907',
    'codex_native_tximport_gene_counts.R': '8aec882b16a9202348f0232c9f7d4365dd0c237386b2f98cea6cfb0219e86582',
}
EXTRA_PER_DONOR = 2_000_000_000
RECEIPT_BUFFER = 20_000_000


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def require(condition, message):
    if not condition:
        raise ValueError(message)


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2)+'\n')


def load_pilot(snapshot):
    for name, digest in SOURCES.items():
        source = HERE / name
        require(sha(source) == digest, 'Successful recipe source changed: '+name)
        (snapshot / name).write_bytes(source.read_bytes())
    # Execute the copied bytes, keeping the original __file__ only to preserve
    # the successful pilot's fixed reference/adapter paths. main() is not called.
    module = types.ModuleType('codex_frozen_successful_b1_recipe')
    module.__file__ = str(HERE / 'codex_native_rna_pilot.py')
    sys.modules[module.__name__] = module
    exec(compile((snapshot / 'codex_native_rna_pilot.py').read_bytes(), module.__file__, 'exec'), module.__dict__)
    return module


def validate_fixed(pilot):
    for name, digest in SOURCES.items():
        require(sha(HERE / name) == digest, 'Recipe changed during group: '+name)
    for path, digest in pilot.FROZEN.items():
        require(pilot.sha(path) == digest, 'Frozen input changed: '+str(path))
    summary = json.loads((pilot.INDEX_DIR / 'summary.json').read_text())
    protocol = json.loads((pilot.INDEX_DIR / 'protocol.json').read_text())
    require(summary['index_built'] and summary['index_returncode'] == 0
            and pilot.sha(pilot.INDEX_DIR / 'protocol.json') == summary['protocol_sha256'], 'RNA index incomplete')
    require(protocol['transcripts'] == 227368 and protocol['decoys'] == 194
            and protocol['k'] == 31 and protocol['keep_duplicates']
            and not protocol['poly_a_clipping'], 'RNA reference recipe changed')
    for name, digest in summary['index_metadata_sha256'].items():
        require(pilot.sha(pilot.INDEX_DIR / 'index' / name) == digest, 'RNA index metadata changed')
    native = pilot.qualify_index_names_lengths()
    version = subprocess.run([str(pilot.ENV / 'salmon'), '--version'],
                             capture_output=True, text=True, check=True).stdout.strip()
    require(version == 'salmon 1.10.3', 'Salmon version changed')
    return native, version, pilot.sha(pilot.INDEX_DIR / 'summary.json')


def donor_rows(roster, donor):
    rows = sorted((r for r in roster if r['assay'] == 'RNA' and r['donor'] == donor),
                  key=lambda r: (r['run_accession'], int(r['file_number'])))
    runs = sorted({r['run_accession'] for r in rows})
    require(runs and len(rows) == 2*len(runs), 'Missing or duplicate RNA mates: '+donor)
    require(len({r['fastq_url'] for r in rows}) == len(rows), 'Repeated RNA file URL: '+donor)
    for run in runs:
        pair = [r for r in rows if r['run_accession'] == run]
        require([r['file_number'] for r in pair] == ['1', '2']
                and len({r['source_biosample'] for r in pair}) == 1
                and all(r['source_biosample'] for r in pair), 'RNA mate identity differs: '+donor)
    expected = sum(int(r['fastq_bytes']) for r in rows)
    require(expected > 0 and all(int(r['fastq_bytes']) > 0 for r in rows), 'Invalid RNA source size')
    return rows, runs, expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--donors', nargs='+', required=True)
    args = parser.parse_args()
    donors = tuple(args.donors)
    require(donors in GROUPS, 'Choose exactly one ordered fixed four-donor group')
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'RNA acquisition requires compute')
    require(int(os.environ['SLURM_CPUS_PER_TASK']) >= 8, 'Requires eight allocated CPUs')
    output = REC / ('codex_native_rna_batch_'+os.environ['SLURM_JOB_ID'])
    require(not output.exists(), 'Refusing output overwrite')
    output.mkdir()
    snapshot = output / 'executed_sources'
    snapshot.mkdir()
    (snapshot / Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    launcher = HERE / 'run_codex_native_rna_batch.sbatch'
    (snapshot / launcher.name).write_bytes(launcher.read_bytes())
    state = {'donors': list(donors), 'status': {d: {'state': 'not_attempted'} for d in donors},
             'group_complete': False, 'cohort_complete': False, 'receiving_H3_read': False,
             'model_fitted': False, 'predictive_accuracy_evaluated': False}
    def save_state():
        write_json(output / 'batch_state.json', state)
    save_state()
    try:
        pilot = load_pilot(snapshot)
        native, version, index_summary_sha = validate_fixed(pilot)
        with pilot.ROSTER.open(newline='') as handle:
            roster = list(csv.DictReader(handle, delimiter='\t'))
        paired = ({r['donor'] for r in roster if r['assay'] == 'RNA'}
                  & {r['donor'] for r in roster if r['assay'] == 'H3K27ac'})
        require(paired == PAIRED, 'Exact fixed paired17 identity roster differs')
        plans = {d: donor_rows(roster, d) for d in donors}
        caps = {d: plans[d][2]+EXTRA_PER_DONOR for d in donors}
        group_cap = sum(caps.values())+RECEIPT_BUFFER
        r_inventory = subprocess.run([str(pilot.ENV / 'Rscript'), '-e',
            'cat(R.version.string,"\\n");cat(as.character(packageVersion("tximport")),"\\n");'
            'cat(as.character(packageVersion("jsonlite")),"\\n")'],
            capture_output=True, text=True, check=True).stdout.splitlines()
        r_inventory = [line.strip() for line in r_inventory if line.strip()]
        require(len(r_inventory) == 3 and r_inventory[1:] == ['1.34.0', '2.0.0'], 'Native R package environment differs')
        receipt = dict(donors=donors, fixed_groups=GROUPS, fixed_paired_donors=sorted(PAIRED),
                       biological_n=len(donors), serial_within_allocation=True,
                       expected_raw_bytes={d: plans[d][2] for d in donors}, donor_caps=caps,
                       group_cap_bytes=group_cap, source_sha256=SOURCES,
                       wrapper_sha256=sha(Path(__file__)), launcher_sha256=sha(launcher),
                       copied_helper_sha256={name: sha(snapshot / name) for name in SOURCES},
                       input_sha256={str(p): digest for p, digest in pilot.FROZEN.items()},
                       native_index_metadata_check=native, index_summary_sha256=index_summary_sha,
                       salmon_version=version, python=sys.version, numpy=pilot.np.__version__,
                       R_version=r_inventory[0], tximport_version=r_inventory[1], jsonlite_version=r_inventory[2],
                       Rscript_sha256=sha(pilot.ENV / 'Rscript'),
                       slurm_job_id=os.environ['SLURM_JOB_ID'], cohort_complete=False,
                       receiving_H3_read=False, model_fitted=False, predictive_accuracy_evaluated=False,
                       private_PISCES_equivalence=False, biological_absence_asserted=False,
                       missing_rows_supplied=False, seed='Successful pilot defaults; no resampling/bias flags',
                       failure_policy='Preserve every donor failure and output; no retries/substitution/complete-case omission')
        write_json(output / 'protocol.json', receipt)
        for donor in donors:
            phase = 'metadata_precheck'
            donor_output = output / donor
            state['status'][donor] = {'state': 'running', 'output': str(donor_output)}
            save_state()
            try:
                native, version, current_index_summary_sha = validate_fixed(pilot)
                require(current_index_summary_sha == index_summary_sha, 'Index receipt changed during group')
                pilot.DONOR, pilot.CAP = donor, caps[donor]
                rows, runs, expected = plans[donor]
                donor_output.mkdir(exist_ok=False)
                raw = donor_output / 'raw'
                raw.mkdir()
                paths = {r['fastq_url']: raw / Path(r['fastq_url']).name for r in rows}
                require(len(set(paths.values())) == len(rows), 'RNA basename collision')
                command = [str(pilot.ENV / 'salmon'), 'quant', '-i', str(pilot.INDEX_DIR / 'index'), '-l', 'A',
                           '-1', *[str(paths[r['fastq_url']]) for r in rows if r['file_number'] == '1'],
                           '-2', *[str(paths[r['fastq_url']]) for r in rows if r['file_number'] == '2'],
                           '-p', str(pilot.THREADS), '-o', str(donor_output / 'quant')]
                counts = donor_output / 'geneCounts.tsv'
                adapter_command = [str(pilot.ENV / 'Rscript'), str(pilot.ADAPTER),
                                   '--quant', str(donor_output / 'quant/quant.sf'), '--tx2gene', str(pilot.MAPPING),
                                   '--sample', donor, '--out', str(counts),
                                   '--report', str(donor_output / 'geneCounts.report.json')]
                donor_receipt = receipt.copy()
                donor_receipt.update(donor=donor, biological_n=1,
                                     technical_runs=runs, expected_download_bytes=expected,
                                     quant_command=command, adapter_command=adapter_command,
                                     command_sha256=hashlib.sha256(json.dumps([command, adapter_command]).encode()).hexdigest(),
                                     input_contract='Complete227368 native transcript estimates then exact42163 modeled denominator; raw frozen global/local')
                write_json(donor_output / 'protocol.json', donor_receipt)
                phase = 'acquisition'
                acquired = []
                for row in rows:
                    acquired.append(pilot.acquire(row, paths[row['fastq_url']]))
                    require(pilot.output_bytes(donor_output) <= caps[donor], 'Donor saved-output cap exceeded')
                    require(pilot.output_bytes(output) <= group_cap, 'Group saved-output cap exceeded')
                    write_json(donor_output / 'acquisition_receipts.json', acquired)
                    print(json.dumps(dict(donor=donor, acquired_files=len(acquired), expected_files=len(rows))), flush=True)
                phase = 'native_quantification'
                pilot.run_checked(command, donor_output / 'salmon_quant.log', donor_output)
                phase = 'native_gene_aggregation'
                pilot.run_checked(adapter_command, donor_output / 'tximport.log', donor_output)
                phase = 'frozen_raw_predictions'
                prediction = pilot.predict_native(counts, donor_output)
                phase = 'final_metadata_check'
                validate_fixed(pilot)
                require(pilot.sha(pilot.INDEX_DIR / 'summary.json') == index_summary_sha, 'Index receipt changed')
                require(sha(Path(__file__)) == receipt['wrapper_sha256'] and sha(launcher) == receipt['launcher_sha256'],
                        'Batch source changed during execution')
                require(sha(pilot.ENV / 'Rscript') == receipt['Rscript_sha256'], 'Rscript changed during group')
                require(pilot.output_bytes(donor_output) <= caps[donor]
                        and pilot.output_bytes(output) <= group_cap, 'Completed saved-output cap exceeded')
                result = dict(**prediction, donor=donor, biological_n=1, acquired_files=len(acquired),
                              saved_output_bytes=pilot.output_bytes(donor_output),
                              protocol_sha256=sha(donor_output / 'protocol.json'),
                              native_quantification_completed=True, frozen_predictions_completed=True,
                              receiving_H3_read=False, model_fitted=False,
                              predictive_accuracy_evaluated=False, independently_validated=False)
                write_json(donor_output / 'summary.json', result)
                state['status'][donor] = {'state': 'succeeded', 'output': str(donor_output),
                                          'summary_sha256': sha(donor_output / 'summary.json')}
            except Exception:
                failure = dict(state='failed', phase=phase, output=str(donor_output), traceback=traceback.format_exc())
                state['status'][donor] = failure
                write_json(output / (donor+'_failure.json'), failure)
            save_state()
        state['group_complete'] = all(state['status'][d]['state'] == 'succeeded' for d in donors)
        state['saved_output_bytes'] = pilot.output_bytes(output)
        save_state()
        write_json(output / 'summary.json', state)
        require(state['group_complete'], 'One or more fixed donors failed; all failures retained')
    except BaseException:
        state['fatal_exception'] = traceback.format_exc()
        for donor in donors:
            previous = state['status'][donor]['state']
            if previous in ('running', 'not_attempted'):
                failure = dict(state='failed', phase='group_fatal',
                               attempted=previous == 'running', output=str(output / donor),
                               traceback=state['fatal_exception'])
                state['status'][donor] = failure
                write_json(output / (donor+'_failure.json'), failure)
        save_state()
        raise


if __name__ == '__main__':
    main()
