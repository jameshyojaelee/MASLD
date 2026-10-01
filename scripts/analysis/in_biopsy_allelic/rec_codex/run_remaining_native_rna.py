"""Apply the piloted native RNA recipe to the remaining sixteen fixed donors.

No receiving H3 measurements, fitting, calibration or predictive loss.
The completed B1 pilot is hash-qualified and reused without overwriting it.
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
ROOT = HERE.parents[3]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
B1 = REC / 'codex_native_rna_pilot_21998663'
INDEX = REC / 'codex_native_rna_index_21998573'
ROSTER = REC / 'paired_h3k27ac_raw_roster_21997742/raw_files.tsv'
DONORS = ('B1', 'B7', 'B8', 'B9', 'B10', 'B12', 'B15', 'B19', 'B21',
          'B22', 'B24', 'B26', 'B36', 'B38', 'B41', 'B46', 'B47')
STOP = 230_000_000_000
CAP = 250_000_000_000
THREADS = 8
SMALL_GUARDS = {
    B1 / 'summary.json': '97d80729969fb1f075d9a18c4134a8d9a94cd570d6f7643dbcfac6ef25c41b46',
    B1 / 'protocol.json': '2bc6fdf97604daece80259b79f0992c9b2ac81f1af94251ffa2fc5f3abcdc0f8',
    B1 / 'geneCounts.report.json': '42b1c7a563d7ee4c2cf36a088eb80536bdf585e001aedec8abce8b7c885c73a4',
    INDEX / 'summary.json': '3461a87af753436971a1340b04eeabdd9bd1ba87169b14869c219fef4a6214f8',
    INDEX / 'protocol.json': '1388892c71c7949bf667b8e0cdf756669a053126daa0c59af4fbc4f4203fec05',
    ROSTER: 'f91c00fb37387020a624a02a388ff00d440efc7fa5bdac5ef41da2a8b1ae7e78',
    HERE / 'codex_native_rna_pilot.py': '413b9617160e0f472682aecc21a629002fae0f54556e1244ffc3aac044e232c5',
    HERE / 'codex_build_native_rna_index.py': '8d6e64a8d270fc7133bdb5a07e56a2963bb305762b417212d01193dba8627c40',
    HERE / 'export_training_h3_target_transform.py': '3499fb927373a04c7fc7121bd89598709774deb38f652b26d5506275e09ee907',
    HERE / 'codex_native_tximport_gene_counts.R': '8aec882b16a9202348f0232c9f7d4365dd0c237386b2f98cea6cfb0219e86582',
}


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def guard(paths):
    for path, expected in paths.items():
        require(sha(path) == expected, 'Frozen/inherited input changed: ' + str(path))


def metadata_plan():
    guard(SMALL_GUARDS)
    first = json.loads((B1 / 'summary.json').read_text())
    protocol = json.loads((B1 / 'protocol.json').read_text())
    index_summary = json.loads((INDEX / 'summary.json').read_text())
    index_protocol = json.loads((INDEX / 'protocol.json').read_text())
    require(first['donor'] == 'B1' and first['native_quantification_completed']
            and first['frozen_predictions_completed'] and first['modeled_genes'] == 42163
            and first['native_genes'] == 60623 and not first['rows_fabricated']
            and not first['receiving_H3_read'] and not first['predictive_accuracy_evaluated'],
            'B1 is not the completed RNA-only native pilot')
    require(set(first['predictions']) == {'global', 'local'}
            and first['predictions']['global']['form'] == 'rrr'
            and first['predictions']['local']['form'] == 'rrr_offset_cis'
            and all(value['transport'] == 'raw' for value in first['predictions'].values()),
            'B1 frozen prediction forms differ')
    require(first['protocol_sha256'] == SMALL_GUARDS[B1 / 'protocol.json']
            and protocol['source_sha256'] == SMALL_GUARDS[HERE / 'codex_native_rna_pilot.py']
            and protocol['adapter_sha256'] == SMALL_GUARDS[HERE / 'codex_native_tximport_gene_counts.R'],
            'B1 source/protocol qualification differs')
    require(index_summary['index_built'] and index_summary['index_returncode'] == 0
            and index_summary['protocol_sha256'] == SMALL_GUARDS[INDEX / 'protocol.json']
            and index_protocol['source_sha256'] == SMALL_GUARDS[HERE / 'codex_build_native_rna_index.py']
            and index_protocol['transcripts'] == 227368 and index_protocol['decoys'] == 194
            and index_protocol['k'] == 31 and index_protocol['keep_duplicates']
            and not index_protocol['poly_a_clipping'], 'Qualified RNA reference recipe differs')
    with ROSTER.open(newline='') as handle:
        all_rows = list(csv.DictReader(handle, delimiter='\t'))
    planned = []
    for donor in DONORS:
        rows = [row for row in all_rows if row['donor'] == donor and row['assay'] == 'RNA']
        rows.sort(key=lambda row: (row['run_accession'], int(row['file_number'])))
        runs = sorted({row['run_accession'] for row in rows})
        require(len(runs) == 8 and len(rows) == 16, 'Fixed donor eight-run/sixteen-file RNA roster differs')
        require(len({row['fastq_url'] for row in rows}) == len(rows), 'Repeated RNA URL')
        for run in runs:
            pair = [row for row in rows if row['run_accession'] == run]
            require([row['file_number'] for row in pair] == ['1', '2']
                    and len({row['source_biosample'] for row in pair}) == 1
                    and len({row['experiment_accession'] for row in pair}) == 1,
                    'Paired mates differ in frozen RNA run/biological source')
        require(len({Path(row['fastq_url']).name for row in rows}) == len(rows), 'RNA basename collision')
        planned.append(dict(donor=donor, runs=runs, files=16,
                            expected_download_bytes=sum(int(row['fastq_bytes']) for row in rows), rows=rows))
    require(planned[0]['runs'] == protocol['technical_runs']
            and planned[0]['expected_download_bytes'] == protocol['expected_download_bytes'],
            'B1 roster differs from completed pilot')
    require(sum(item['expected_download_bytes'] for item in planned) == 79800518482,
            'All-seventeen frozen RNA bytes differ')
    require(sum(item['expected_download_bytes'] for item in planned[1:]) == 73538554197,
            'Remaining frozen RNA bytes differ')
    require(len({row['fastq_url'] for item in planned for row in item['rows']}) == 272,
            'RNA files are repeated across fixed donors')
    return first, protocol, index_summary, index_protocol, planned


def load_archived_module(name, original, data):
    require(name not in sys.modules, 'Numerical helper already imported: ' + name)
    module = types.ModuleType(name)
    # Preserve original relative-root and __file__ conventions, while executing
    # exactly the source bytes already archived and checked by this workflow.
    module.__file__ = str(original)
    sys.modules[name] = module
    exec(compile(data, str(original), 'exec'), module.__dict__)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true', help='small metadata/source checks only; no arrays, reads, subprocesses or downloads')
    args = parser.parse_args()
    first, inherited, index_summary, index_protocol, plan = metadata_plan()
    if args.inspect:
        print(json.dumps(dict(donor_order=list(DONORS), reused_B1=str(B1),
                              remaining=[{k: v for k, v in item.items() if k != 'rows'} for item in plan[1:]],
                              remaining_runs=128, remaining_files=256,
                              expected_remaining_download_bytes=73538554197,
                              sampled_stop_bytes=STOP, hard_cap_bytes=CAP), indent=2))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Native RNA work requires SLURM')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) >= THREADS, 'Requires eight allocated CPUs')
    output = REC / ('remaining_native_rna_' + os.environ['SLURM_JOB_ID'])
    require(not output.exists(), 'Refusing workflow output overwrite')
    fixed = {Path(path): expected for path, expected in inherited['input_sha256'].items()}
    fixed.update(SMALL_GUARDS)
    fixed.update({Path(path): expected for path, expected in index_protocol['input_sha256'].items()})
    fixed.update({INDEX / 'index' / name: expected for name, expected in index_summary['index_metadata_sha256'].items()})
    rscript = Path(inherited['quant_command'][0]).with_name('Rscript')
    fixed[rscript] = '750e8efcde0d04740fa754e5dc3cf22b66f5c176c7898a73d18e5b7d9635ad02'
    guard(fixed)
    b1_report = json.loads((B1 / 'geneCounts.report.json').read_text())
    require(b1_report['sample'] == 'B1' and b1_report['native_transcripts'] == 227368
            and not b1_report['rows_padded_or_removed'], 'B1 native aggregation receipt differs')
    reused = {B1 / 'geneCounts.tsv': b1_report['output_sha256']['gene_counts'],
              B1 / 'geneCounts.tsv.native_zero_transcripts.tsv': b1_report['output_sha256']['native_zero_transcripts'],
              B1 / 'quant/quant.sf': b1_report['inputs']['quant_sha256']}
    reused.update({Path(value['path']): value['sha256'] for value in first['predictions'].values()})
    require(all(path.resolve().is_relative_to(B1.resolve()) for path in reused), 'B1 output path escaped pilot')
    guard(reused)
    payload = {path: sha(path) for path in sorted((INDEX / 'index').rglob('*')) if path.is_file()}
    require(payload, 'Index payload is empty')
    scorer_original = Path(next(path for path in inherited['input_sha256'] if path.endswith('/score.py')))
    originals = [HERE / 'export_training_h3_target_transform.py', HERE / 'codex_build_native_rna_index.py',
                 HERE / 'codex_native_rna_pilot.py', HERE / 'codex_native_tximport_gene_counts.R', scorer_original]
    source_bytes = {path: path.read_bytes() for path in originals}
    for path, data in source_bytes.items():
        require(hashlib.sha256(data).hexdigest() == fixed[path], 'Source changed before archival')
    output.mkdir(exist_ok=False)
    archive = output / 'executed_sources'
    archive.mkdir()
    for path, data in source_bytes.items():
        (archive / path.name).write_bytes(data)
    (archive / Path(__file__).name).write_bytes(Path(__file__).read_bytes())
    load_archived_module('export_training_h3_target_transform', originals[0], source_bytes[originals[0]])
    builder = load_archived_module('codex_build_native_rna_index', originals[1], source_bytes[originals[1]])
    pilot = load_archived_module('codex_native_rna_pilot', originals[2], source_bytes[originals[2]])
    pilot.ADAPTER = archive / 'codex_native_tximport_gene_counts.R'
    pilot.SCORER = archive / 'score.py'
    pilot.CAP = STOP
    native_index = pilot.qualify_index_names_lengths()
    require(pilot.ROOT.resolve() == ROOT.resolve() and pilot.DONOR == 'B1' and pilot.THREADS == THREADS,
            'Inherited source/root/thread convention differs')
    require(subprocess.run([str(pilot.ENV / 'salmon'), '--version'], check=True,
                           capture_output=True, text=True).stdout.strip() == inherited['salmon_version'],
            'Salmon version differs from pilot')
    guard(fixed)
    phase, donor, completed = 'initialization', None, []

    def size_guard(stop=True):
        size = builder.output_bytes(output)
        require(size <= CAP, 'Hard remaining RNA cap exceeded; partial outputs preserved')
        if stop:
            require(size <= STOP, 'Sampled remaining RNA stop exceeded; partial outputs preserved')
        return size

    def write_json(path, value, fresh=False):
        data = (json.dumps(value, indent=2) + '\n').encode()
        require(not fresh or not path.exists(), 'Refusing existing result overwrite: ' + str(path))
        old = path.stat().st_size if path.exists() else 0
        require(size_guard() - old + len(data) <= STOP, 'Receipt/result would exceed sampled stop')
        with path.open('xb' if fresh else 'wb') as handle:
            handle.write(data)
        size_guard()

    protocol = dict(fixed_donor_order=list(DONORS), remaining_order=list(DONORS[1:]),
                    reused_B1_summary_sha256=SMALL_GUARDS[B1 / 'summary.json'],
                    source_inheritance='exact archived pilot/helper/index/scorer/adapter bytes; original __file__ semantics; DONOR set explicitly per fixed iteration',
                    source_sha256={str(path): hashlib.sha256(data).hexdigest() for path, data in source_bytes.items()},
                    fixed_input_sha256={str(path): value for path, value in fixed.items()},
                    index_payload_sha256={str(path): value for path, value in payload.items()},
                    B1_output_sha256={str(path): value for path, value in reused.items()},
                    native_index_metadata_check=native_index,
                    remaining_roster=[{k: v for k, v in item.items() if k != 'rows'} for item in plan[1:]],
                    expected_remaining_download_bytes=73538554197, expected_remaining_files=256,
                    sampled_stop_bytes=STOP, hard_cap_bytes=CAP, monitor_seconds=30,
                    threads=THREADS, requested_memory='64G', requested_wall_time='90:00:00',
                    allocated_cpus=os.environ.get('SLURM_CPUS_PER_TASK'),
                    allocated_memory_mb=os.environ.get('SLURM_MEM_PER_NODE'),
                    campaign_file_allowance_bytes=1500000000000,
                    declared_source_global_overrides=dict(DONOR='per fixed remaining donor', CAP=STOP,
                                                          ADAPTER=str(pilot.ADAPTER), SCORER=str(pilot.SCORER)),
                    library_type='A before mates', bias_flags=[], countsFromAbundance='no',
                    quantification_stochastic_seed=inherited['quantification_stochastic_seed'],
                    fitted_models={'global': 'rrr', 'local': 'rrr_offset_cis'}, transport='raw',
                    structural_native_zero_estimates_retained=True, private_PISCES_equivalence=False,
                    rows_fabricated=False, H3_read=False, models_fitted=False, predictive_accuracy_evaluated=False,
                    independently_validated=False, python=sys.version, numpy=pilot.np.__version__,
                    wrapper_sha256=sha(archive / Path(__file__).name),
                    launcher_sha256=sha(HERE / 'run_remaining_native_rna.sbatch'),
                    argv=sys.argv, slurm_job_id=os.environ['SLURM_JOB_ID'])
    write_json(output / 'protocol.json', protocol, fresh=True)
    try:
        for item in plan[1:]:
            donor, phase = item['donor'], 'donor_initialization'
            pilot.DONOR = donor
            guard(fixed)
            guard({archive / path.name: fixed[path] for path in originals})
            folder = output / donor
            folder.mkdir(exist_ok=False)
            raw = folder / 'raw'
            raw.mkdir()
            rows = item['rows']
            paths = {row['fastq_url']: raw / Path(row['fastq_url']).name for row in rows}
            command = [str(pilot.ENV / 'salmon'), 'quant', '-i', str(pilot.INDEX_DIR / 'index'), '-l', 'A',
                       '-1', *[str(paths[row['fastq_url']]) for row in rows if row['file_number'] == '1'],
                       '-2', *[str(paths[row['fastq_url']]) for row in rows if row['file_number'] == '2'],
                       '-p', str(THREADS), '-o', str(folder / 'quant')]
            adapter_command = [str(rscript), str(pilot.ADAPTER), '--quant', str(folder / 'quant/quant.sf'),
                               '--tx2gene', str(pilot.MAPPING), '--sample', donor,
                               '--out', str(folder / 'geneCounts.tsv'), '--report', str(folder / 'geneCounts.report.json')]
            write_json(folder / 'protocol.json', dict(donor=donor, technical_runs=item['runs'],
                       expected_download_bytes=item['expected_download_bytes'], rows=rows,
                       quant_command=command, adapter_command=adapter_command,
                       source_global_DONOR=pilot.DONOR, workflow_protocol_sha256=sha(output / 'protocol.json')), fresh=True)
            acquired, phase = [], 'acquisition'
            for row in rows:
                require(size_guard() + int(row['fastq_bytes']) <= STOP, 'Next raw file would exceed sampled stop')
                acquired.append(pilot.acquire(row, paths[row['fastq_url']]))
                write_json(folder / 'acquisition_receipts.json', acquired)
                print(json.dumps(dict(donor=donor, phase=phase, acquired_files=len(acquired), expected_files=16)), flush=True)
            phase = 'salmon_quantification'
            pilot.run_checked(command, folder / 'salmon_quant.log', output)
            phase = 'native_tximport'
            pilot.run_checked(adapter_command, folder / 'tximport.log', output)
            phase = 'frozen_predictions'
            prediction = pilot.predict_native(folder / 'geneCounts.tsv', folder)
            size_guard()
            guard(fixed)
            guard({archive / path.name: fixed[path] for path in originals})
            result = dict(**prediction, donor=donor, biological_n=1, acquired_files=16,
                          expected_download_bytes=item['expected_download_bytes'],
                          native_quantification_completed=True, frozen_predictions_completed=True,
                          protocol_sha256=sha(folder / 'protocol.json'), H3_read=False,
                          models_fitted=False, predictive_accuracy_evaluated=False, independently_validated=False)
            write_json(folder / 'summary.json', result, fresh=True)
            completed.append(dict(donor=donor, directory=str(folder), summary_sha256=sha(folder / 'summary.json'),
                                  predictions=prediction['predictions']))
            write_json(output / 'progress_summary.json', dict(completed_remaining=len(completed),
                       expected_remaining=16, completed=completed, saved_output_bytes=size_guard(),
                       H3_read=False, predictive_accuracy_evaluated=False))
            print(json.dumps(dict(donor=donor, phase='complete', completed_remaining=len(completed), expected_remaining=16)), flush=True)
        phase = 'final_hash_and_cap_checks'
        guard(fixed)
        guard(reused)
        guard(payload)
        require(set(payload) == {path for path in (INDEX / 'index').rglob('*') if path.is_file()},
                'Index payload inventory changed during remaining RNA calculation')
        guard({archive / path.name: fixed[path] for path in originals})
        require(len(completed) == 16 and tuple(item['donor'] for item in completed) == DONORS[1:],
                'Remaining fixed donor population incomplete')
        final = dict(completed_remaining=16, fixed_total_donors=17, B1_reused_directory=str(B1),
                     B1_summary_sha256=SMALL_GUARDS[B1 / 'summary.json'], completed=completed,
                     native_quantification_completed=True, frozen_predictions_completed=True,
                     saved_output_bytes_before_summary=size_guard(), protocol_sha256=sha(output / 'protocol.json'),
                     H3_read=False, predictive_accuracy_evaluated=False, models_fitted=False,
                     independently_validated=False, success=True)
        write_json(output / 'summary.json', final, fresh=True)
        require(builder.output_bytes(output) <= CAP, 'Final output including summary exceeds hard cap')
        print(json.dumps(dict(completed_remaining=16, fixed_total_donors=17, success=True)), flush=True)
    except BaseException:
        error = dict(donor=donor, phase=phase, completed_remaining=len(completed),
                     success=False, traceback=traceback.format_exc(), partial_outputs_preserved=True)
        try:
            write_json(output / 'failure.json', error, fresh=True)
        except Exception:
            print(json.dumps(error), file=sys.stderr, flush=True)
        raise


if __name__ == '__main__':
    main()
