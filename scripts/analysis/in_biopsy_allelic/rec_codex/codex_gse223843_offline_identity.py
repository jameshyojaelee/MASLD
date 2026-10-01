#!/usr/bin/env python3
"""Offline public feature/barcode identity checks; no network or count inputs."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import traceback
import zlib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ACQUIRED = REC / 'codex_gse223843_feature_barcode_ranges_21999459'
PROPOSAL = REC / 'codex_gse223843_geometry_acquisition_proposal.json'
IDENTITY = REC / 'codex_gse223843_identity_metadata/identity_metadata.json'
GEOMETRY = REC / 'codex_gse296875_metadata_census_21999159/native_custom_merged_atac_geometry.tsv'
CLASSIFICATION = REC / 'codex_gse223843_first_feature_schema_20261001.json'
GUARDS = {
    PROPOSAL: '35107ad4b17a7463bf28831191d7327ea895e28776b652f85d58ec20360cd7e8',
    IDENTITY: '172aa3c42134a2a4f329cdf24db937d92c70da0a1300827f32c6cf13b3754018',
    GEOMETRY: '50a52d097a5232edcf136dd53d6549ca2d1aa5dce9a57cf055ff10b7b2e0f3e3',
    CLASSIFICATION: 'ea7eae8f9d96ad1a0781dfd32bb9086f80d5cabe5a1121c0ec8aec5cc2e5f6f8',
    ACQUIRED / 'failure.json': 'ad5218429c3bd239a0cc46ae3e81b69b71480a901c1aa68fa87ee561fba2799a',
    ACQUIRED / 'tar_identity_offsets.json': '8e0467bcb07978e943cbbe69046fc359e4860677b61511a39bba1569aaf4f878',
    ACQUIRED / 'executed_source.py': '7f9845fd6ba4e756502b246fafd1500a55aa4b3987720e0569fc82c174224691',
    ACQUIRED / 'official_filelist.txt': 'd0e3fe3012e57aa3d0364378d90357e400dbdfb6346a2ebd892d6d13a74c295a',
}
FILE_CAP, TOTAL_CAP, SAVED_CAP = 64 * 1024**2, 512 * 1024**2, 64 * 1024**2
ROW_CAP, LINE_CAP = 2000000, 4096
FLAGS = dict(count_matrix_read=False, count_payload_inspected=False,
    clinical_values_read=False, genotype_values_read=False, model_scores_read=False,
    lineage_assignment=False, person_independence_verified=False,
    model_fit_or_selection=False, model_or_evaluation_admission=False,
    network_requests_made=False, zero_fill_or_intersection_target_created=False)


def require(test, message):
    if not test:
        raise ValueError(message)


class SchemaFailure(ValueError):
    def __init__(self, message, field_count):
        super().__init__(message)
        self.field_count = field_count


class DecodeLimit(ValueError):
    pass


class AggregateDecodeLimit(DecodeLimit):
    pass


def schema(test, message, field_count):
    if not test:
        raise SchemaFailure(message, field_count)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(65536), b''):
            digest.update(block)
    return digest.hexdigest()


def saved_bytes(output):
    return sum(path.stat().st_size for path in output.iterdir() if path.is_file())


def write_json(output, name, value):
    raw = (json.dumps(value, indent=2, allow_nan=False) + '\n').encode()
    require(saved_bytes(output) + len(raw) <= SAVED_CAP, 'Inclusive output cap exceeded')
    with (output / name).open('xb') as target:
        target.write(raw)


class Decoder:
    def __init__(self):
        self.total_bytes = 0
        self.file_bytes = 0
        self.file_rows = 0

    def rows(self, path):
        self.file_bytes, self.file_rows = 0, 0
        with path.open('rb') as source:
            schema(source.read(2) == b'\x1f\x8b', 'Expected gzip magic', None)
        with gzip.open(path, 'rb') as source:
            while True:
                raw = source.readline(LINE_CAP + 1)
                if not raw:
                    break  # Full successful iteration verifies gzip CRC/end.
                self.file_bytes += len(raw); self.total_bytes += len(raw); self.file_rows += 1
                if self.total_bytes > TOTAL_CAP:
                    raise AggregateDecodeLimit('Decoded aggregate cap exceeded')
                if len(raw) > LINE_CAP or self.file_bytes > FILE_CAP or self.file_rows > ROW_CAP:
                    raise DecodeLimit('Decoded line/file/row cap exceeded')
                try:
                    text = raw.rstrip(b'\r\n').decode('utf-8', errors='strict')
                except UnicodeDecodeError:
                    raise SchemaFailure('Non-UTF8 identity text; contents suppressed', None) from None
                schema(bool(text) and '\x00' not in text, 'Empty/NUL identity row', None)
                yield text, raw


def fixture_axis():
    ids, bed = set(), set()
    with GEOMETRY.open(newline='') as source:
        reader = csv.DictReader(source, delimiter='\t')
        require(reader.fieldnames == ['peak_id', 'chromosome', 'source_start_1based_closed',
            'source_end_1based_closed', 'bed_start_0based', 'bed_end_half_open'], 'Frozen geometry columns differ')
        for row in reader:
            a, b, c, d = [int(row[key]) for key in ('source_start_1based_closed',
                'source_end_1based_closed', 'bed_start_0based', 'bed_end_half_open')]
            coordinate = (row['chromosome'], c, d)
            require(a >= 1 and b >= a and c == a - 1 and d == b
                    and row['peak_id'] not in ids and coordinate not in bed, 'Invalid/duplicate frozen geometry')
            ids.add(row['peak_id']); bed.add(coordinate)
    require(len(ids) == len(bed) == 306706, 'Frozen full feature axis differs')
    return ids, bed


def chromosome(token):
    return re.fullmatch(r'chr[A-Za-z0-9_.]+', token) is not None


def narrow_identity(token):
    return re.fullmatch(r'[A-Za-z0-9_.-]+', token) is not None and any(c.isalpha() for c in token)


def validate_rna(fields):
    n = len(fields)
    schema(n in (1, 2, 3, 6), 'Unsupported RNA field count', n)
    schema(narrow_identity(fields[0]), 'Unsupported RNA identity syntax', n)
    schema(n == 1 or narrow_identity(fields[1]), 'Unsupported RNA symbol syntax', n)
    schema(n < 3 or fields[2] == 'Gene Expression', 'Unsupported RNA feature type; no rows filtered', n)
    if n == 6:
        schema(re.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?', fields[0]) is not None
            and chromosome(fields[3]) and re.fullmatch(r'[0-9]+', fields[4]) is not None
            and re.fullmatch(r'[0-9]+', fields[5]) is not None,
            'Unsupported six-field gene-coordinate annotation', n)
        schema(int(fields[5]) >= int(fields[4]), 'Reversed RNA gene annotation coordinates', n)


def rna_features(path, decoder):
    ids, genes, symbols, chromosomes = set(), set(), set(), Counter()
    digest, genes_digest, examples, types, width = hashlib.sha256(), hashlib.sha256(), [], Counter(), None
    for text, raw in decoder.rows(path):
        fields = text.split('\t')
        validate_rna(fields)  # No row examples before narrow schema validation.
        schema(width is None or width == len(fields), 'RNA schema width changes within file', len(fields))
        width = len(fields)
        schema(fields[0] not in ids, 'Duplicate RNA feature identity', width)
        ids.add(fields[0]); digest.update(raw)
        kind = fields[2] if width >= 3 else 'type_not_declared'
        types[kind] += 1
        if width >= 2:
            symbols.add(fields[1])
        if re.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?', fields[0]):
            genes.add(fields[0]); genes_digest.update(('\t'.join(fields[:3]) + '\n').encode())
        if width == 6:
            chromosomes[fields[3]] += 1
        if len(examples) < 3:
            examples.append(text)
    schema(bool(ids), 'Empty RNA feature axis', width)
    return dict(rows=len(ids), field_count=width, feature_types=dict(types),
        gene_ids_matching_ENSG_syntax=len(genes), symbol_identity_count=len(symbols),
        decoded_axis_sha256=digest.hexdigest(), hash_codec='original decoded bytes including line endings',
        ordered_gene_id_symbol_type_sha256=genes_digest.hexdigest(),
        ordered_gene_hash_codec='ENSG-row first up-to-three source fields separated by TAB followed by LF; no version stripping',
        first_three_validated_identity_examples=examples, gene_annotation_contigs=dict(chromosomes),
        gene_annotation_coordinate_origin_verified=False, feature_identity_does_not_measure_gene_expression=True)


def coordinate(fields):
    n = len(fields)
    schema(n in (1, 3), 'Unsupported ATAC field count', n)
    if n == 1:
        match = re.fullmatch(r'([^:\t]+):(\d+)-(\d+)', fields[0])
        if match is None:
            match = re.fullmatch(r'([^:\t]+)-(\d+)-(\d+)', fields[0])
        schema(match is not None, 'Unsupported ATAC interval string syntax', n)
        chrom, start, end = match.groups()
    else:
        chrom, start, end = fields
        schema(re.fullmatch(r'[0-9]+', start) is not None and re.fullmatch(r'[0-9]+', end) is not None,
               'ATAC coordinate fields are not unsigned integers', n)
    schema(re.fullmatch(r'[A-Za-z0-9_.]+', chrom) is not None, 'Unsupported ATAC chromosome syntax', n)
    a, b = int(start), int(end)
    schema(b >= a and (a >= 1 or b > a), 'ATAC interval invalid under both coordinate hypotheses', n)
    return chrom, a, b


def atac_features(path, decoder, fixed_ids, fixed_bed):
    coordinates, exact, zero, one = set(), set(), set(), set()
    digest, examples, width, rows, zero_invalid, one_invalid = hashlib.sha256(), [], None, 0, 0, 0
    for text, raw in decoder.rows(path):
        fields = text.split('\t')
        chrom, start, end = coordinate(fields)
        schema(width is None or width == len(fields), 'ATAC schema width changes within file', len(fields))
        width = len(fields); key = (chrom, start, end)
        schema(key not in coordinates, 'Duplicate ATAC coordinate identity', width)
        coordinates.add(key); rows += 1; digest.update(raw)
        if text in fixed_ids:
            exact.add(text)
        if end > start:
            if key in fixed_bed:
                zero.add(key)
        else:
            zero_invalid += 1
        if start >= 1 and end >= start:
            shifted = (chrom, start - 1, end)
            if shifted in fixed_bed:
                one.add(shifted)
        else:
            one_invalid += 1
        if len(examples) < 3:
            examples.append(text)
    schema(rows > 0, 'Empty ATAC feature axis', width)
    def hypothesis(matches, invalid):
        return dict(fixed_axis_rows=306706, exact_coordinate_matches=len(matches),
            missing_fixed_coordinates=306706 - len(matches), source_rows_invalid_under_hypothesis=invalid,
            source_rows_not_matching_fixed_axis=rows - invalid - len(matches),
            complete_coordinate_coverage_under_hypothesis=len(matches) == 306706,
            coordinate_origin_verified=False, count_compatibility_verified=False)
    return dict(rows=rows, field_count=width, decoded_axis_sha256=digest.hexdigest(),
        hash_codec='original decoded bytes including line endings', first_three_validated_geometry_examples=examples,
        exact_original_string_matches=len(exact), missing_fixed_original_strings=306706 - len(exact),
        hypotheses=dict(zero_based_half_open=hypothesis(zero, zero_invalid), one_based_closed=hypothesis(one, one_invalid)),
        coordinate_hypothesis_selected=False, missing_coordinates_not_observed_zeros=True,
        original_strings_preserved_in_source_compressed_member=True)


def barcode_file(path, decoder):
    values, digest, examples = set(), hashlib.sha256(), []
    for text, raw in decoder.rows(path):
        schema(re.fullmatch(r'[ACGT]{8,32}-[1-9]\d*', text) is not None,
               'Expected one canonical 10x barcode identity field', len(text.split('\t')))
        schema(text not in values, 'Duplicate within-assay barcode', 1)
        values.add(text); digest.update(raw)
        if len(examples) < 3:
            examples.append(text)
    schema(bool(values), 'Empty barcode axis', 1)
    return values, dict(rows=len(values), decoded_axis_sha256=digest.hexdigest(),
        hash_codec='original decoded bytes including line endings', first_three_validated_barcode_examples=examples,
        unique_within_source_sample_assay=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-script-sha256', required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Full offline identity checks require compute')
    require(sha(Path(__file__)) == args.expected_script_sha256, 'Reviewed offline helper source changed')
    output = REC / ('codex_gse223843_offline_identity_' + os.environ['SLURM_JOB_ID'])
    owned, checked, phase, decoder = False, {}, 'fresh_output', Decoder()
    try:
        output.mkdir(exist_ok=False); owned = True
        phase = 'frozen_receipt_and_fixture_hashes'
        for path, expected in GUARDS.items():
            checked[str(path)] = sha(path)
            require(checked[str(path)] == expected, 'Frozen receipt/identity/fixture source changed')
        proposal = json.loads(PROPOSAL.read_text()); identity = json.loads(IDENTITY.read_text())
        failure = json.loads((ACQUIRED / 'failure.json').read_text())
        offsets = json.loads((ACQUIRED / 'tar_identity_offsets.json').read_text())['allowlisted_offsets']
        members = proposal['members']
        require(len(members) == len({r['member_name'] for r in members}) == 24
            and sum(r['declared_compressed_bytes'] for r in members) == 16668757
            and len({(r['source_person_label'], r['assay'], r['role']) for r in members}) == 24,
            'Frozen member roster differs')
        require(failure['count_payload_inspected'] is False and failure['network_response_body_bytes'] == 16690996,
            'Original acquisition evidence differs')
        phase = 'every_acquired_member_size_and_sha'
        receipts = {}
        for row in members:
            name = row['member_name']; path = ACQUIRED / name; item = offsets[name]
            require(Path(name).name == name and row['role'] in ('features', 'barcodes')
                and row['assay'] in ('snRNA-seq', 'snATAC-seq')
                and row['accession'] == identity['person_axis_by_exact_author_labels'][row['source_person_label']][row['assay']]
                and name.startswith(row['accession'] + '_') and name in identity['listed_member_filenames'],
                'Unexpected acquired member role/person/path')
            ranges = [r for r in failure['HTTP_range_receipts']
                if r['purpose'] == 'allowlisted_identity_member' and r['start'] == item['offset']]
            require(len(ranges) == 1, 'Ambiguous original member range')
            receipt = ranges[0]
            require(receipt['status'] == 206 and receipt['bytes_read'] == receipt['bytes_expected']
                == item['bytes'] == row['declared_compressed_bytes'] == path.stat().st_size
                and receipt['end'] == receipt['start'] + receipt['bytes_read'] - 1,
                'Acquired member size/range evidence differs')
            digest = sha(path)
            require(digest == receipt['sha256'], 'Acquired member SHA differs from successful range receipt')
            receipts[name] = dict(path=str(path), bytes=path.stat().st_size, sha256=digest, **row)
        for name, path in [('executed_source.py', Path(__file__)),
            ('executed_launcher.sbatch', HERE / 'run_codex_gse223843_offline_identity.sbatch'),
            ('source_proposal.json', PROPOSAL), ('source_identity.json', IDENTITY),
            ('original_failure.json', ACQUIRED / 'failure.json'),
            ('original_tar_offsets.json', ACQUIRED / 'tar_identity_offsets.json')]:
            raw = path.read_bytes()
            require(saved_bytes(output) + len(raw) <= SAVED_CAP - 1048576, 'Source-copy output cap exceeded')
            with (output / name).open('xb') as target:
                target.write(raw)
        phase = 'full_fixed_geometry_identity'
        fixed_ids, fixed_bed = fixture_axis()
        results, barcode_sets = {}, {}
        for row in members:
            name = row['member_name']; path = ACQUIRED / name
            phase = 'per_file_identity/' + name
            try:
                if row['role'] == 'barcodes':
                    values, result = barcode_file(path, decoder)
                    barcode_sets[(row['source_person_label'], row['assay'])] = values
                elif row['assay'] == 'snRNA-seq':
                    result = rna_features(path, decoder)
                else:
                    result = atac_features(path, decoder, fixed_ids, fixed_bed)
                results[name] = dict(status='identity_schema_and_CRC_verified', gzip_CRC_verified=True,
                    decoded_bytes=decoder.file_bytes, **result)
            except (SchemaFailure, DecodeLimit, gzip.BadGzipFile, EOFError, zlib.error) as error:
                if isinstance(error, AggregateDecodeLimit):
                    raise
                # No rejected row values, valid prefix examples or partial axes are emitted.
                results[name] = dict(status='bounded_decoding_failure' if isinstance(error, DecodeLimit)
                    else 'unsupported_schema_or_compression_failure',
                    exception_type=type(error).__name__, reason=str(error) if isinstance(error, (SchemaFailure, DecodeLimit))
                    else 'Gzip integrity failure; decoded contents suppressed',
                    rejected_schema_field_count=error.field_count if isinstance(error, SchemaFailure) else None,
                    rejected_field_contents_saved=False, full_file_schema_verified=False,
                    gzip_CRC_verified=False, decoded_prefix_bytes=decoder.file_bytes,
                    partial_axes_not_reported_as_complete=True)
            # Preserve progress after each file, even if a later bound aborts the run.
            write_json(output, 'file_' + name + '.json', results[name])
        phase = 'sample_namespaced_barcode_joins'
        pairs = {}
        for person, assays in identity['person_axis_by_exact_author_labels'].items():
            rna = barcode_sets.get((person, 'snRNA-seq')); atac = barcode_sets.get((person, 'snATAC-seq'))
            pair = dict(author_person_label=person, RNA_accession=assays['snRNA-seq'], ATAC_accession=assays['snATAC-seq'],
                person_independence_verified=False, barcode_namespace=person,
                actual_well_identity_not_declared_by_this_source_metadata=True)
            if rna is None or atac is None:
                pair.update(join_status='unavailable_due_to_barcode_file_failure',
                    RNA_barcode_axis_available=rna is not None, ATAC_barcode_axis_available=atac is not None)
            else:
                pair.update(join_status='exact_sample_namespaced_identity_join', RNA_barcodes=len(rna),
                    ATAC_barcodes=len(atac), shared_exact_barcodes=len(rna & atac),
                    RNA_only_barcodes=len(rna - atac), ATAC_only_barcodes=len(atac - rna), full_sets_equal=rna == atac,
                    no_cross_person_barcode_suffix_matching=True, no_count_row_filtering_or_admission=True)
            pairs[person] = pair
        phase = 'completion_hashes_and_inclusive_output'
        require(all(sha(path) == expected for path, expected in GUARDS.items()), 'Frozen source changed during offline check')
        require(all(sha(Path(row['path'])) == row['sha256'] for row in receipts.values()), 'Acquired member changed during check')
        failed = [name for name, result in results.items() if result['status'] != 'identity_schema_and_CRC_verified']
        summary = dict(offline_inspection_completed=True, all_24_files_attempted=True,
            full_identity_qualification_completed=not failed, unsupported_or_failed_files=failed,
            qualified_files=24 - len(failed), source_acquisition_job='21999459', original_job_preserved=True,
            slurm_job_id=os.environ['SLURM_JOB_ID'], checked_source_sha256=checked,
            executed_source_sha256=sha(output / 'executed_source.py'),
            executed_launcher_sha256=sha(output / 'executed_launcher.sbatch'),
            compressed_source_member_receipts=receipts, per_file_results=results, sample_barcode_joins=pairs,
            fixed_target_axis_rows=306706, coordinate_origin_verified=False,
            coordinate_hypotheses_not_selected_by_best_overlap=True,
            compressed_gene_feature_members_acquired=True,
            full_gene_feature_schema_verified=all(results[r['member_name']]['status'] == 'identity_schema_and_CRC_verified'
                for r in members if r['assay'] == 'snRNA-seq' and r['role'] == 'features'),
            decoded_metadata_bytes=decoder.total_bytes, decoded_file_cap=FILE_CAP, decoded_total_cap=TOTAL_CAP,
            saved_output_cap=SAVED_CAP, saved_bytes_before_summary=saved_bytes(output),
            RNA_barcodes_are_identities_not_expression_values=True,
            gene_feature_identity_is_not_genotype_or_observed_gene_expression=True,
            source_population='six pediatric author-person labels; independence unverified',
            historical_source_exposure_or_cross_source_person_disjointness_verified=False,
            python=sys.version, gzip_module=gzip.__file__, zlib_runtime=zlib.ZLIB_RUNTIME_VERSION, **FLAGS)
        write_json(output, 'summary.json', summary)
        require(saved_bytes(output) <= SAVED_CAP, 'Final inclusive output cap exceeded')
        print(json.dumps(dict(offline_inspection_completed=True, qualified_files=24-len(failed),
            unsupported_or_failed_files=failed, output=str(output), saved_bytes=saved_bytes(output), **FLAGS)))
    except BaseException as error:
        if owned:
            failure = dict(offline_inspection_completed=False, phase=phase,
                exception_type=type(error).__name__, reason=str(error), traceback=traceback.format_exc(),
                checked_source_sha256=checked, partial_progress_receipts_preserved=True,
                decoded_metadata_bytes=decoder.total_bytes, **FLAGS)
            try:
                write_json(output, 'failure.json', failure)
            except Exception as receipt_error:
                print('Failure receipt write failed: ' + str(receipt_error), file=sys.stderr)
        raise


if __name__ == '__main__':
    main()
