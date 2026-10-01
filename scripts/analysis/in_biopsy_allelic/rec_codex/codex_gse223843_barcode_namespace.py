#!/usr/bin/env python3
"""Offline ARC codebook and processed-barcode identity qualification; no counts."""
import argparse
import collections
import contextlib
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time
import zlib

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ARC = Path('/nfs/sw/easybuild/software/cellranger-arc-2.0.2')
OLD = REC / 'codex_gse223843_offline_identity_21999912'
ANN = REC / 'codex_gse223843_rna_annotation_identity_22000649'
ACQ = REC / 'codex_gse223843_feature_barcode_ranges_21999459'
PERSONS = ('Control1', 'Control3', 'Fontan1', 'Fontan2', 'Fontan3', 'Fontan4')
AXIS_HASH = 'afbcc02143e68d1fbcda9631917f1344adfd3ac8efb752e2ebb7f044df934e60'
PINS = {
    REC / '20261001_codex_gse296875_donor_global_rna_atac_prespec.json': '88b9e4239b08cc1ccbe78b092e8eacd395778453a496fd8915b4acc9e2729ae1',
    REC / 'codex_gse223843_geometry_acquisition_proposal.json': '35107ad4b17a7463bf28831191d7327ea895e28776b652f85d58ec20360cd7e8',
    REC / 'codex_gse223843_identity_metadata/identity_metadata.json': '172aa3c42134a2a4f329cdf24db937d92c70da0a1300827f32c6cf13b3754018',
    ACQ / 'failure.json': 'ad5218429c3bd239a0cc46ae3e81b69b71480a901c1aa68fa87ee561fba2799a',
    ACQ / 'tar_identity_offsets.json': '8e0467bcb07978e943cbbe69046fc359e4860677b61511a39bba1569aaf4f878',
    OLD / 'summary.json': '0ef0b29278da9e38512204fa3ba05a9c5311fd6d71b7aed425e6692cb046268d',
    ANN / 'summary.json': 'da3f1310cfba1bc145ded47979a56e4134e2870baf09a25d553fc598024a7abe',
    ARC / 'lib/python/cellranger/barcodes/737K-arc-v1.txt.gz': '1279bb7ba28ddcd805e2f16ebd5d6933b7b1e7f128e18aabfe6ae8cd7b7bf710',
    ARC / 'lib/python/atac/barcodes/737K-arc-v1.txt.gz': '9a4c8f29a05f59b9a4d7525a2f44d1778ccfa540a8832ad9f5f84047f374cc26',
    ARC / 'lib/python/cellranger/barcodes/utils.py': 'e9c5c807d91af0b5ea9fc1573ab590fd48e0d0ade711057ad8aea36d78d070f9',
    ARC / 'lib/python/atac/barcodes/__init__.py': 'c896d07f2bde284d6dd36a618636e2c83d9b131a2195af627dc9eebd8c9069ba',
    ARC / 'lib/python/atac_rna/barcodes.py': '9f4e0c90eec48ac8bbadd0bd69377c8fca9b1a17cac59b216283dd50fc82f185',
    ARC / 'lib/python/cellranger/constants.py': '5783e51dfb65c5a6eb557f3aa4f22b6c9829944bbbb027168968804a2891b819',
    ARC / 'lib/python/cellranger/chemistry_defs.json': '5fb5b204d8900362a23124c72e3e66c2ad0feb28db8535995082be3d08f637a5',
    ARC / 'mro/atac_rna/sc_atac_gex_counter_cs.mro': '283bcc35bca83cd819d48a7b35d449bbed1bec2e3043f15089f6f4c50152c0ac',
    ARC / 'lib/python/atac_rna/matrix.py': '3795c9c494845e208c4f10ae9f53bb7e51199badec394221e17e3aa06caf47f8',
    ARC / 'mro/atac_rna/stages/processing/merge_atac_rna_matrices/__init__.py': 'eb1cbb59f50afb811110d790f558270f7e225aef30a07a2babace09dddc30ac3',
}
TOTAL_CAP = 128 * 1024**2
FILE_DECODE_CAP = 32 * 1024**2
MAX_LINE = 4096
MAX_ROWS = 1_000_000
TIME_CAP = 55 * 60
FAILURE_RESERVE = 8192
BARCODE = re.compile(rb'([ACGT]{16})-([1-9][0-9]{0,8})')
BASE = re.compile(rb'[ACGT]{16}')


def require(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def pin(path, expected, size=None):
    require(path.is_file(), 'source_file_missing')
    if size is not None:
        require(path.stat().st_size == size, 'source_size_changed')
    observed = digest(path)
    require(observed == expected, 'source_sha256_changed')
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': observed}


def read_json(path):
    require(path.stat().st_size <= 1024**2, 'metadata_file_too_large')
    return json.loads(path.read_text())


class Budget:
    def __init__(self, out, logs):
        self.out = out
        self.logs = logs
        self.start = time.monotonic()
        self.decoded = 0

    def saved(self):
        return sum(p.stat().st_size for p in self.out.iterdir() if p.is_file()) + sum(
            p.stat().st_size for p in self.logs if p.is_file())

    def check(self, extra=0):
        require(time.monotonic() - self.start < TIME_CAP, 'time_cap_exceeded')
        require(self.saved() + extra <= TOTAL_CAP - FAILURE_RESERVE, 'saved_byte_cap_exceeded')

    @contextlib.contextmanager
    def writer(self, name):
        # Reserve before each write, including bytes still in Python's buffer.
        base_bytes = self.saved()
        budget = self
        with (self.out / name).open('xb') as stream:
            class BoundedWriter:
                written = 0

                def write(self, body):
                    require(base_bytes + self.written + len(body) <= TOTAL_CAP - FAILURE_RESERVE,
                            'saved_byte_cap_exceeded')
                    self.written += len(body)
                    return stream.write(body)

                def flush(self):
                    stream.flush()
                    budget.check()
            yield BoundedWriter()
            stream.flush()
        self.check()

    def json(self, name, value):
        body = (json.dumps(value, indent=2, sort_keys=True) + '\n').encode()
        self.check(len(body))
        with (self.out / name).open('xb') as stream:
            stream.write(body)

    def lines(self, path, state):
        # Reading to EOF verifies every gzip member's CRC. No truncated-prefix success.
        with gzip.open(path, 'rb') as stream:
            while True:
                line = stream.readline(MAX_LINE + 1)
                if not line:
                    break
                require(len(line) <= MAX_LINE, 'identity_line_too_long')
                state['bytes'] += len(line)
                self.decoded += len(line)
                require(state['bytes'] <= FILE_DECODE_CAP, 'per_file_decode_cap_exceeded')
                require(self.decoded <= TOTAL_CAP, 'aggregate_decode_cap_exceeded')
                state['hash'].update(line)
                yield line
                if state['bytes'] % (1024**2) < MAX_LINE:
                    self.check()


def canonical(path, label, budget):
    state = {'bytes': 0, 'hash': hashlib.sha256()}
    values, index = [], {}
    comments = 0
    normalized = hashlib.sha256()
    with budget.writer(label + '_canonical_order.tsv') as out:
        for raw in budget.lines(path, state):
            # Identical canonical-loader rule; comment lines are not codebook slots.
            if b'#' in raw:
                comments += 1
                continue
            value = raw.strip()
            require(BASE.fullmatch(value) is not None, 'canonical_barcode_schema_invalid')
            require(value not in index, 'canonical_barcode_duplicate')
            require(len(values) < MAX_ROWS, 'canonical_row_cap_exceeded')
            index[value] = len(values)
            values.append(value)
            out.write(value + b'\n')
            normalized.update(value + b'\n')
        out.flush()
    budget.check()
    return values, index, {'rows': len(values), 'unique': True, 'gzip_CRC_verified': True,
        'raw_decoded_bytes': state['bytes'], 'raw_decoded_sha256': state['hash'].hexdigest(),
        'ordered_canonical_sha256': normalized.hexdigest(), 'comment_lines_ignored': comments,
        'order_codec': 'canonical loader strip rule, original noncomment order, barcode_LF',
        'orientation_choice': False, 'sorted_or_error_corrected': False}


def namespace(base, rna_index, atac_index):
    in_rna, in_atac = base in rna_index, base in atac_index
    return 'both' if in_rna and in_atac else 'RNA_only' if in_rna else 'ATAC_only' if in_atac else 'neither'


def processed(path, expected, rna_index, atac_index, budget):
    state = {'bytes': 0, 'hash': hashlib.sha256()}
    values, index = [], {}
    classes, gems = collections.Counter(), collections.Counter()
    with budget.writer(path.name + '.ordered_identity.tsv') as out:
        out.write(b'native_row_1based\tbarcode_original\tGEM_suffix_original\tnamespace_membership\tRNA_codebook_slot_1based\tATAC_codebook_slot_1based\n')
        for raw in budget.lines(path, state):
            value = raw.rstrip(b'\r\n')
            match = BARCODE.fullmatch(value)
            require(match is not None, 'processed_barcode_schema_invalid')
            require(value not in index, 'processed_barcode_duplicate')
            require(len(values) < MAX_ROWS, 'processed_row_cap_exceeded')
            index[value] = len(values)
            values.append(value)
            base, suffix = match.groups()
            category = namespace(base, rna_index, atac_index)
            classes[category] += 1
            gems[suffix.decode()] += 1
            rna_slot = str(rna_index[base] + 1) if base in rna_index else 'NA'
            atac_slot = str(atac_index[base] + 1) if base in atac_index else 'NA'
            out.write(b'\t'.join([str(len(values)).encode(), value, suffix, category.encode(),
                                  rna_slot.encode(), atac_slot.encode()]) + b'\n')
        out.flush()
    require(len(values) == expected['rows'], 'prior_barcode_row_count_changed')
    require(state['hash'].hexdigest() == expected['decoded_axis_sha256'], 'prior_decoded_barcode_axis_changed')
    require(state['bytes'] == expected['decoded_bytes'], 'prior_decoded_barcode_bytes_changed')
    budget.check()
    return values, index, {'rows': len(values), 'unique_within_source_sample_assay': True,
        'gzip_CRC_verified': True, 'decoded_bytes': state['bytes'],
        'decoded_axis_sha256': state['hash'].hexdigest(), 'namespace_membership_counts': dict(classes),
        'GEM_suffix_counts': dict(gems), 'all_unique_to_GEX_namespace': classes['RNA_only'] == len(values),
        'namespaces_ambiguous_or_unrecognized': classes['both'] + classes['neither'] > 0,
        'source_order_and_suffix_preserved': True, 'processed_ID_translation_performed': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-script-sha256', required=True)
    args = parser.parse_args()
    job = os.environ.get('SLURM_JOB_ID', '')
    require(re.fullmatch(r'[1-9][0-9]*', job) is not None, 'compute_SLURM_job_required')
    require(not os.environ.get('SLURM_ARRAY_TASK_ID'), 'array_not_allowed')
    source = Path(__file__).resolve()
    require(digest(source) == args.expected_script_sha256, 'executing_source_sha256_changed')
    launcher = source.with_name('run_codex_gse223843_barcode_namespace.sbatch')
    out = REC / ('codex_gse223843_barcode_namespace_' + job)
    logs = [REC / ('1592_barcode_namespace_' + job + suffix) for suffix in ('.out', '.err')]
    owned_output = False
    phase = 'source_guards'
    budget = None
    try:
        inventory = {str(path): pin(path, sha) for path, sha in PINS.items()}
        old = read_json(OLD / 'summary.json')
        ann = read_json(ANN / 'summary.json')
        identity = read_json(REC / 'codex_gse223843_identity_metadata/identity_metadata.json')
        proposal = read_json(REC / 'codex_gse223843_geometry_acquisition_proposal.json')
        failure = read_json(ACQ / 'failure.json')
        offsets = read_json(ACQ / 'tar_identity_offsets.json')['allowlisted_offsets']
        require(ann['all_six_full_annotation_qualified'] is True, 'six_annotations_not_qualified')
        require(tuple(identity['person_axis_by_exact_author_labels']) == PERSONS, 'source_person_order_changed')
        annotation_receipts = {}
        for person in PERSONS:
            name = person + '_annotation_receipt.json'
            expected = ann['output_artifacts'][name]
            path = ANN / name
            inventory[str(path)] = pin(path, expected['sha256'], expected['bytes'])
            receipt = read_json(path)
            require(receipt['person'] == person and receipt['full_annotation_schema_verified']
                    and receipt['gzip_CRC_passed'], 'annotation_receipt_not_qualified')
            require(receipt['validated_Gene_Expression_rows'] == 36601 and
                    receipt['validated_ordered_Gene_Expression_ID_symbol_sha256'] == AXIS_HASH,
                    'fixed_gene_axis_not_exact')
            require(receipt['source_feature_types']['Gene Expression'] == 36601 and
                    set(receipt['source_feature_types']) == {'Gene Expression', 'Peaks'}, 'unexpected_feature_types')
            annotation_receipts[person] = {
                'path': str(path), 'sha256': expected['sha256'],
                'source_feature_types': receipt['source_feature_types'],
                'ordered_Gene_Expression_ID_symbol_sha256': AXIS_HASH,
                'annotation_not_count_qualification': True}
        members = [row for row in proposal['members'] if row['role'] == 'barcodes']
        require(len(members) == 12, 'barcode_member_count_changed')
        by_person = {p: {} for p in PERSONS}
        for member in members:
            name = member['member_name']
            require(Path(name).name == name, 'member_path_not_basename')
            person, assay = member['source_person_label'], member['assay']
            require(person in by_person and assay in ('snRNA-seq', 'snATAC-seq') and
                    assay not in by_person[person], 'source_person_assay_duplicate_or_unknown')
            require(identity['person_axis_by_exact_author_labels'][person][assay] == member['accession'],
                    'source_accession_assignment_changed')
            require(name in identity['listed_member_filenames'], 'member_absent_from_official_identity')
            receipt = old['compressed_source_member_receipts'][name]
            size = member['declared_compressed_bytes']
            require(receipt['source_person_label'] == person and receipt['assay'] == assay and
                    receipt['accession'] == member['accession'] and receipt['role'] == 'barcodes' and
                    receipt['bytes'] == size and offsets[name]['bytes'] == size,
                    'compressed_source_identity_changed')
            source_path = ACQ / name
            require(Path(receipt['path']).resolve() == source_path.resolve(), 'source_member_path_changed')
            range_matches = [row for row in failure['HTTP_range_receipts']
                if row['start'] == offsets[name]['offset'] and row['end'] == offsets[name]['offset'] + size - 1]
            require(len(range_matches) == 1, 'original_range_receipt_missing_or_duplicate')
            original_range = range_matches[0]
            require(original_range['status'] == 206 and original_range['bytes_read'] == size and
                    original_range['bytes_expected'] == size and original_range['sha256'] == receipt['sha256'],
                    'original_range_not_exact')
            inventory[str(source_path)] = pin(source_path, receipt['sha256'], size)
            previous = old['per_file_results'][name]
            require(previous['status'] == 'identity_schema_and_CRC_verified' and previous['gzip_CRC_verified'],
                    'prior_barcode_not_qualified')
            by_person[person][assay] = (member, source_path, previous, original_range)
        require(all(set(v) == {'snRNA-seq', 'snATAC-seq'} for v in by_person.values()), 'source_assay_pair_missing')
        inventory[str(source)] = pin(source, args.expected_script_sha256)
        inventory[str(launcher)] = {'path': str(launcher), 'bytes': launcher.stat().st_size, 'sha256': digest(launcher)}
        out.mkdir(mode=0o700)
        owned_output = True
        budget = Budget(out, logs)
        shutil.copyfile(source, out / 'executed_source.py')
        shutil.copyfile(launcher, out / 'executed_launcher.sbatch')
        budget.json('source_contract.json', {'source_inventory': inventory,
            'annotation_receipts': annotation_receipts, 'author_person_order': PERSONS,
            'codebook_pairing': 'canonical original order zip, never sort-pairing',
            'source_same_nucleus_methods_pointer': identity['same_nucleus_pairing_source'],
            'barcodes_not_expression_or_genotype': True, 'matrix_read_authorized': False})
        phase = 'full_canonical_whitelists_CRC_and_bijection'
        rna, rna_index, rna_receipt = canonical(ARC / 'lib/python/cellranger/barcodes/737K-arc-v1.txt.gz', 'RNA', budget)
        atac, atac_index, atac_receipt = canonical(ARC / 'lib/python/atac/barcodes/737K-arc-v1.txt.gz', 'ATAC', budget)
        require(len(rna) == len(atac) and len(rna) > 0, 'canonical_lists_unequal_or_empty')
        same_slots = 0
        with budget.writer('canonical_ordered_bijection.tsv') as table:
            table.write(b'canonical_pair_slot_1based\tRNA_barcode_original\tATAC_barcode_original\n')
            for slot, (rna_base, atac_base) in enumerate(zip(rna, atac), 1):
                require(rna_index[rna_base] == slot - 1 and atac_index[atac_base] == slot - 1,
                        'canonical_bijection_inverse_failed')
                same_slots += rna_base == atac_base
                table.write(str(slot).encode() + b'\t' + rna_base + b'\t' + atac_base + b'\n')
                if slot % 65536 == 0:
                    table.flush()
                    budget.check()
            table.flush()
        codebook = {'RNA': rna_receipt, 'ATAC': atac_receipt, 'equal_lengths': True,
            'bijection_and_inverse_verified': True, 'pairing_rule': 'unaltered original ordered zip',
            'same_string_namespace_intersection': sum(base in atac_index for base in rna),
            'same_string_at_same_pair_slot': same_slots, 'overlap_used_to_choose_translation': False,
            'GEM_suffix_translation_rule': 'keep original suffix unchanged; no processed translation executed'}
        budget.json('canonical_codebook_receipt.json', codebook)
        files, rosters = {}, {}
        for person in PERSONS:
            axes = {}
            for assay in ('snRNA-seq', 'snATAC-seq'):
                member, path, previous, original_range = by_person[person][assay]
                phase = 'processed_barcode_CRC_namespace/' + member['member_name']
                values, index, receipt = processed(path, previous, rna_index, atac_index, budget)
                receipt.update({'source_person_label': person, 'assay': assay,
                    'accession': member['accession'], 'compressed_source': inventory[str(path)],
                    'original_HTTP_range_receipt': original_range})
                files[member['member_name']] = receipt
                axes[assay] = (values, index, receipt)
            phase = 'exact_author_person_roster/' + person
            rna_values, rna_rows, rna_file = axes['snRNA-seq']
            atac_values, atac_rows, atac_file = axes['snATAC-seq']
            shared = 0
            with budget.writer(person + '_exact_source_paired_barcode_roster.tsv') as roster:
                roster.write(b'author_person\tRNA_accession\tATAC_accession\tbarcode_original\tRNA_source_row_1based\tATAC_source_row_1based\tGEM_suffix_original\tnamespace_membership\n')
                for barcode in rna_values:
                    if barcode not in atac_rows:
                        continue
                    base, suffix = BARCODE.fullmatch(barcode).groups()
                    shared += 1
                    roster.write(b'\t'.join([person.encode(), by_person[person]['snRNA-seq'][0]['accession'].encode(),
                        by_person[person]['snATAC-seq'][0]['accession'].encode(), barcode,
                        str(rna_rows[barcode] + 1).encode(), str(atac_rows[barcode] + 1).encode(),
                        suffix, namespace(base, rna_index, atac_index).encode()]) + b'\n')
                roster.flush()
            expected = old['sample_barcode_joins'][person]
            observed = {'RNA_barcodes': len(rna_values), 'ATAC_barcodes': len(atac_values),
                'shared_exact_barcodes': shared, 'RNA_only_barcodes': len(rna_values) - shared,
                'ATAC_only_barcodes': len(atac_values) - shared}
            require(all(expected[key] == value for key, value in observed.items()), 'prior_exact_person_join_changed')
            require(expected['RNA_accession'] == by_person[person]['snRNA-seq'][0]['accession'] and
                    expected['ATAC_accession'] == by_person[person]['snATAC-seq'][0]['accession'], 'prior_join_accession_changed')
            rosters[person] = dict(observed, author_person_label=person,
                RNA_accession=expected['RNA_accession'], ATAC_accession=expected['ATAC_accession'],
                exact_existing_join_reproduced=True, roster_order='RNA native barcode file order',
                all_source_barcodes_unique_to_GEX_namespace=(rna_file['all_unique_to_GEX_namespace'] and
                                                            atac_file['all_unique_to_GEX_namespace']),
                no_count_or_outcome_based_eligibility=True, technical_identity_roster_only=True,
                source_methods_same_nucleus_pointer=identity['same_nucleus_pairing_source'],
                actual_well_identity_independently_verified=False,
                same_nucleus_identity_proven_by_barcode_coincidence=False,
                biological_person_disjointness_verified=False,
                translated_IDs_or_unknown_rows_dropped=False)
            budget.check()
        phase = 'final_source_accounting'
        after = {path: pin(Path(path), row['sha256'], row['bytes']) for path, row in inventory.items()}
        require(after == inventory, 'final_source_inventory_changed')
        artifacts = {p.name: {'bytes': p.stat().st_size, 'sha256': digest(p)}
                     for p in out.iterdir() if p.is_file()}
        summary = {'barcode_namespace_inspection_completed': True, 'SLURM_job_id': job,
            'all_12_processed_files_full_CRC_and_identity_verified': len(files) == 12,
            'canonical_codebook': codebook, 'per_file_results': files, 'sample_paired_rosters': rosters,
            'all_six_annotation_receipts': annotation_receipts, 'source_inventory': inventory,
            'source_inventory_before_after_identical': True, 'output_artifacts_before_summary': artifacts,
            'decoded_bytes': budget.decoded, 'decode_cap': TOTAL_CAP, 'per_file_decode_cap': FILE_DECODE_CAP,
            'saved_byte_cap_including_scheduler_logs': TOTAL_CAP,
            'saved_bytes_before_summary': budget.saved(), 'elapsed_seconds': time.monotonic() - budget.start,
            'Python_version': sys.version, 'zlib_runtime': zlib.ZLIB_RUNTIME_VERSION,
            'network_requests': 0, 'feature_payloads_decoded': False, 'matrix_or_counts_read': False,
            'raw_sequence_or_quality_read': False, 'clinical_or_genotype_or_allelic_values_read': False,
            'H3_or_outcomes_or_model_scores_read': False, 'lineage_assignment': False,
            'orientation_selection_or_error_correction': False, 'processed_barcode_translation': False,
            'source_roster_chosen_from_canonical_overlap': False, 'cross_person_deduplication': False,
            'model_fit_or_selection_or_admission': False, 'receiving_measurement_admission': False,
            'RNA_count_unit_or_denominator_verified_from_matrix': False,
            'full306706_ATAC_measurement_verified': False, 'original_jobs_and_files_preserved': True}
        budget.json('summary.json', summary)
        budget.check()
        print(json.dumps({'completed': True, 'phase': 'metadata_only', 'summary': str(out / 'summary.json'),
                          'summary_sha256': digest(out / 'summary.json'), 'saved_bytes': budget.saved()}))
        sys.stdout.flush()
        budget.check()
    except Exception as exc:
        # No source values or arbitrary exception text are emitted on schema failure.
        reason = str(exc) if isinstance(exc, RuntimeError) else type(exc).__name__
        report = {'completed': False, 'phase': phase, 'exception_type': type(exc).__name__,
            'categorical_reason': reason, 'partials_preserved': True, 'matrix_or_counts_read': False,
            'sequence_or_quality_read': False, 'network_requests': 0, 'model_or_measurement_admission': False}
        if owned_output:
            # Tiny categorical receipt even if normal-output cap was reached; never overwrite.
            try:
                with (out / 'failure.json').open('x') as stream:
                    json.dump(report, stream, indent=2)
                    stream.write('\n')
            except Exception:
                pass
        print(json.dumps(report), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
