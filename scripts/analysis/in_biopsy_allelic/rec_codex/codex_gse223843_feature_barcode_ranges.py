#!/usr/bin/env python3
"""Acquire only frozen public feature/barcode TAR members; never count payloads."""
import argparse
from collections import Counter
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import sys
import traceback
from urllib.request import HTTPRedirectHandler, Request, build_opener
import zlib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
PROPOSAL = REC / 'codex_gse223843_geometry_acquisition_proposal.json'
PROPOSAL_SHA = '35107ad4b17a7463bf28831191d7327ea895e28776b652f85d58ec20360cd7e8'
IDENTITY = REC / 'codex_gse223843_identity_metadata/identity_metadata.json'
IDENTITY_SHA = '172aa3c42134a2a4f329cdf24db937d92c70da0a1300827f32c6cf13b3754018'
GEOMETRY = REC / 'codex_gse296875_metadata_census_21999159/native_custom_merged_atac_geometry.tsv'
GEOMETRY_SHA = '50a52d097a5232edcf136dd53d6549ca2d1aa5dce9a57cf055ff10b7b2e0f3e3'
ARCHIVE_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223843/suppl/GSE223843_RAW.tar'
FILELIST_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223843/suppl/filelist.txt'
FILELIST_SHA = 'd0e3fe3012e57aa3d0364378d90357e400dbdfb6346a2ebd892d6d13a74c295a'
ARCHIVE_SIZE = 1066270720
NETWORK_CAP = 20000000
METADATA_CAP = 2000000
DECODED_FILE_CAP = 64 * 1024 * 1024
DECODED_TOTAL_CAP = 512 * 1024 * 1024
SAVED_CAP = 64 * 1024 * 1024
ROW_CAP = 2000000
LINE_CAP = 4096
FLAGS = dict(count_payload_inspected=False, clinical_values_inspected=False,
             lineage_assignment=False, person_independence_verified=False,
             genotype_values_inspected=False, model_scores_inspected=False,
             model_or_evaluation_admission=False)


def require(test, message):
    if not test:
        raise ValueError(message)


class SchemaFailure(ValueError):
    def __init__(self, message, field_count):
        super().__init__(message)
        self.field_count = field_count


def schema_require(test, message, field_count):
    if not test:
        raise SchemaFailure(message, field_count)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def saved_bytes(output):
    return sum(path.stat().st_size for path in output.rglob('*') if path.is_file())


def write_json(output, name, value):
    payload = (json.dumps(value, indent=2, allow_nan=False) + '\n').encode()
    require(saved_bytes(output) + len(payload) <= SAVED_CAP, 'Inclusive saved-output cap exceeded')
    with (output / name).open('xb') as target:
        target.write(payload)


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, request, fp, code, message, headers, new_url):
        raise ValueError('Redirect forbidden; only the exact frozen official URLs are allowed')


class Acquisition:
    def __init__(self, output):
        self.output = output
        self.opener = build_opener(NoRedirect())
        self.network_bytes = 0
        self.metadata_bytes = 0
        self.decoded_bytes = 0
        self.ranges = []

    def range(self, start, length, metadata=False, destination=None):
        require(start >= 0 and length > 0 and start + length <= ARCHIVE_SIZE,
                'Range outside frozen archive')
        require(self.network_bytes + length <= NETWORK_CAP, 'Network byte cap exceeded')
        if metadata:
            require(self.metadata_bytes + length <= METADATA_CAP, 'TAR metadata cap exceeded')
        end = start + length - 1
        request = Request(ARCHIVE_URL, headers={'Range': f'bytes={start}-{end}',
            'Accept-Encoding': 'identity', 'User-Agent': 'bounded-public-identity/1'})
        with self.opener.open(request, timeout=45) as response:
            # Validate before any response-body read. A 200 response is closed here.
            require(response.status == 206 and response.geturl() == ARCHIVE_URL,
                    'Server did not honor exact HTTP range; no archive fallback')
            require(response.headers.get('Content-Range') == f'bytes {start}-{end}/{ARCHIVE_SIZE}',
                    'Content-Range differs from exact request/frozen archive size')
            value = response.headers.get('Content-Length')
            require(value is None or value.isdecimal() and int(value) == length,
                    'Content-Length inconsistent with exact range')
            require(response.headers.get('Content-Encoding', 'identity').lower() == 'identity',
                    'Encoded HTTP response would change byte offsets')
            entry = dict(start=start, end=end, bytes_expected=length,
                         purpose='TAR_metadata' if metadata else 'allowlisted_identity_member',
                         bytes_read=0, status=response.status)
            self.ranges.append(entry)
            body = bytearray() if destination is None else None
            digest = hashlib.sha256()
            target = destination.open('xb') if destination is not None else None
            try:
                remaining = length
                while remaining:
                    block = response.read(min(65536, remaining))
                    require(block, 'Truncated exact-range response')
                    self.network_bytes += len(block)
                    if metadata:
                        self.metadata_bytes += len(block)
                    entry['bytes_read'] += len(block)
                    remaining -= len(block)
                    digest.update(block)
                    if target is None:
                        body.extend(block)
                    else:
                        require(saved_bytes(self.output) + len(block) <= SAVED_CAP - 1048576,
                                'Saved-output cap reached; partial retained')
                        target.write(block)
                # No extra-byte probe: that could cross into a forbidden matrix body.
                entry['sha256'] = digest.hexdigest()
            finally:
                if target is not None:
                    target.close()
            return bytes(body) if body is not None else dict(entry)

    def filelist(self):
        request = Request(FILELIST_URL, headers={'Accept-Encoding': 'identity'})
        with self.opener.open(request, timeout=45) as response:
            require(response.status == 200 and response.geturl() == FILELIST_URL,
                    'Unexpected official filelist response')
            value = response.headers.get('Content-Length')
            require(value is None or value == '2783', 'Filelist declared bytes changed')
            data = response.read(2784)
        self.network_bytes += len(data)
        require(self.network_bytes <= NETWORK_CAP and len(data) == 2783
                and hashlib.sha256(data).hexdigest() == FILELIST_SHA,
                'Official filelist hash/bytes changed')
        with (self.output / 'official_filelist.txt').open('xb') as target:
            target.write(data)
        reader = csv.DictReader(io.StringIO(data.decode('utf-8')), delimiter='\t')
        require(reader.fieldnames == ['#Archive/File', 'Name', 'Time', 'Size', 'Type'],
                'Filelist header changed')
        files, archive = {}, []
        for row in reader:
            require(None not in row and row['Size'].isdecimal(), 'Malformed filelist row')
            if row['#Archive/File'] == 'Archive':
                archive.append(row)
            else:
                require(row['#Archive/File'] == 'File' and row['Name'] not in files,
                        'Duplicate/unknown filelist entry')
                files[row['Name']] = int(row['Size'])
        require(len(archive) == 1 and archive[0]['Name'] == 'GSE223843_RAW.tar'
                and int(archive[0]['Size']) == ARCHIVE_SIZE, 'Frozen archive identity changed')
        return files


def tar_number(field):
    require(not field[0] & 128, 'Base-256 TAR number unsupported; stop closed')
    value = field.strip(b'\x00 ')
    require(not value or re.fullmatch(b'[0-7]+', value), 'Malformed octal TAR number')
    return int(value or b'0', 8)


def safe_path(raw):
    name = raw.decode('utf-8', errors='strict')
    while name.startswith('./'):
        name = name[2:]
    require(name and not name.startswith('/') and '\\' not in name
            and all(ord(c) >= 32 and ord(c) != 127 for c in name), 'Unsafe TAR path')
    require(all(part not in ('', '.', '..') for part in name.rstrip('/').split('/')),
            'Unsafe TAR path components')
    return name


def pax_fields(raw):
    fields, position = {}, 0
    while position < len(raw):
        space = raw.find(b' ', position)
        require(space > position and raw[position:space].isdigit(), 'Malformed PAX record length')
        length = int(raw[position:space])
        require(length > space - position + 2 and position + length <= len(raw), 'PAX length outside payload')
        record = raw[space + 1:position + length]
        require(record.endswith(b'\n') and b'=' in record, 'Malformed PAX key/value record')
        key, value = record[:-1].decode('utf-8').split('=', 1)
        require(key not in fields and key in {'path', 'size', 'mtime', 'atime', 'ctime',
                'uid', 'gid', 'uname', 'gname', 'charset'}, 'Unsupported/duplicate PAX key; stop closed')
        fields[key] = value
        position += length
    if 'path' in fields:
        safe_path(fields['path'].encode())
    if 'size' in fields:
        require(fields['size'].isdecimal(), 'Invalid PAX size')
    require(fields.get('charset', 'UTF-8') == 'UTF-8', 'Unsupported PAX charset')
    return fields


def scan_archive(acquisition, file_sizes, allowlist):
    offsets, seen, position, pending = {}, set(), 0, {}
    while position + 1024 <= ARCHIVE_SIZE:
        header = acquisition.range(position, 512, metadata=True)
        if header == bytes(512):
            require(not pending and acquisition.range(position + 512, 512, metadata=True) == bytes(512),
                    'Missing second TAR end block or unapplied PAX metadata')
            require((ARCHIVE_SIZE - position) % 512 == 0, 'Archive padding alignment differs')
            require(seen == set(file_sizes) and set(offsets) == allowlist,
                    'Archive names/allowlist incomplete or differ from official filelist')
            return offsets, dict(regular_members=len(seen), two_end_blocks_verified=True,
                trailing_padding_bytes_not_read=ARCHIVE_SIZE - position - 1024,
                matrix_bodies_read=False, matrix_members_skipped=sum(name.endswith('_matrix.mtx.gz') for name in seen))
        require(tar_number(header[148:156]) == sum(header[:148]) + 8 * 32 + sum(header[156:]),
                'TAR header checksum differs')
        require(header[257:263] in (b'ustar\x00', b'ustar ', bytes(6)), 'Unsupported TAR format')
        name = header[:100].split(b'\x00', 1)[0]
        if header[257:263] == b'ustar\x00':
            prefix = header[345:500].split(b'\x00', 1)[0]
            if prefix:
                name = prefix + b'/' + name
        name = safe_path(name)
        size, kind = tar_number(header[124:136]), header[156:157]
        payload_start = position + 512
        require(payload_start + size <= ARCHIVE_SIZE, 'TAR member body outside archive')
        if kind == b'x':
            require(not pending and 0 < size <= 65536, 'Chained/oversized PAX extension')
            pending = pax_fields(acquisition.range(payload_start, size, metadata=True))
        else:
            require(kind in (b'0', b'\x00', b'5'), 'Links/sparse/global/GNU extensions unsupported; stop closed')
            name = safe_path(pending.get('path', name).encode())
            if 'size' in pending:
                size = int(pending['size'])
            pending = {}
            require(payload_start + size <= ARCHIVE_SIZE, 'PAX member outside archive')
            if kind == b'5':
                require(size == 0, 'Nonempty directory member unsupported')
            else:
                require(name in file_sizes and name not in seen and size == file_sizes[name],
                        'TAR regular name/size differs from frozen official inventory')
                seen.add(name)
                if name in allowlist:
                    offsets[name] = dict(offset=payload_start, bytes=size)
                else:
                    require(name.endswith('_matrix.mtx.gz'), 'Nonallowlisted payload is not a declared matrix')
        # All regular bodies, including every matrix, are skipped without a read.
        position = payload_start + ((size + 511) // 512) * 512
    raise ValueError('Archive end missing within frozen size')


def lines(path, acquisition):
    rows, decoded = 0, 0
    with path.open('rb') as probe:
        require(probe.read(2) == b'\x1f\x8b', 'Member does not have gzip magic')
    with gzip.open(path, 'rb') as source:
        while True:
            raw = source.readline(LINE_CAP + 1)
            if not raw:
                break  # gzip verifies CRC/end while reaching EOF.
            decoded += len(raw)
            acquisition.decoded_bytes += len(raw)
            rows += 1
            require(len(raw) <= LINE_CAP and decoded <= DECODED_FILE_CAP
                    and acquisition.decoded_bytes <= DECODED_TOTAL_CAP and rows <= ROW_CAP,
                    'Decoded member/line/row/aggregate bound exceeded')
            try:
                text = raw.rstrip(b'\r\n').decode('utf-8', errors='strict')
            except UnicodeDecodeError:
                raise SchemaFailure('Identity text is not UTF-8; contents not recorded', None) from None
            require(text and '\x00' not in text, 'Empty/NUL identity row')
            yield text, raw


def barcodes(path, acquisition):
    values, digest, examples = set(), hashlib.sha256(), []
    for value, raw in lines(path, acquisition):
        schema_require(re.fullmatch(r'[ACGT]{8,32}-[1-9]\d*', value) is not None,
                'Barcode schema rejected; expected one canonical 10x identity field',
                len(value.split('\t')))
        require(value not in values, 'Duplicate barcode within source sample/assay')
        values.add(value); digest.update(raw)
        if len(examples) < 3:
            examples.append(value)
    require(values, 'Empty barcode axis')
    return values, dict(rows=len(values), decoded_axis_sha256=digest.hexdigest(),
        hash_codec='original UTF8 decompressed bytes including line endings', first_three=examples,
        duplicates=0, barcode_suffixes_not_used_for_cross_person_join=True)


def frozen_geometry():
    ids, bed = set(), set()
    with GEOMETRY.open(newline='') as source:
        reader = csv.DictReader(source, delimiter='\t')
        require(reader.fieldnames == ['peak_id', 'chromosome', 'source_start_1based_closed',
            'source_end_1based_closed', 'bed_start_0based', 'bed_end_half_open'], 'Frozen geometry header differs')
        for row in reader:
            a, b, c, d = [int(row[key]) for key in ['source_start_1based_closed',
                'source_end_1based_closed', 'bed_start_0based', 'bed_end_half_open']]
            key = (row['chromosome'], c, d)
            require(a >= 1 and b >= a and c == a - 1 and d == b and row['peak_id'] not in ids
                    and key not in bed, 'Frozen geometry is invalid/duplicate')
            ids.add(row['peak_id']); bed.add(key)
    require(len(ids) == len(bed) == 306706, 'Frozen geometry axis length differs')
    return ids, bed


def peak_coordinate(fields):
    if len(fields) == 1:
        match = re.fullmatch(r'([^:\t]+):(\d+)-(\d+)', fields[0])
        if match is None:
            match = re.fullmatch(r'([^:\t]+)-(\d+)-(\d+)', fields[0])
        require(match is not None, 'Unrecognized ATAC identity format; coordinate origin not inferred')
        chrom, start, end = match.groups()
        form = 'one string: chromosome:start-end or chromosome-start-end'
    else:
        require(len(fields) == 3 and fields[1].isdigit() and fields[2].isdigit(),
                'Unrecognized ATAC fields; coordinate origin not inferred')
        chrom, start, end = fields
        form = 'three tab-separated fields: chromosome,start,end'
    start, end = int(start), int(end)
    require(start >= 0 and end >= start and (start >= 1 or end > start)
            and re.fullmatch(r'[A-Za-z0-9_.]+', chrom),
            'Invalid/empty source ATAC interval or chromosome identity')
    return (chrom, start, end), form


def features(path, assay, acquisition, fixed_ids, fixed_bed):
    digest, examples, fields_seen, types = hashlib.sha256(), [], Counter(), Counter()
    identities, raw_coordinates, matches, zero_matches, one_matches = set(), set(), set(), set(), set()
    genes, gene_digest, formats, rows = set(), hashlib.sha256(), set(), 0
    for text, raw in lines(path, acquisition):
        fields = text.split('\t')
        schema_require(not fields_seen or len(fields) in fields_seen,
                'Feature field count changed; rejected schema, contents not recorded', len(fields))
        if assay == 'snRNA-seq':
            schema_require(len(fields) in (1, 2, 3), 'RNA schema rejected', len(fields))
            # No numeric expression-like row may be exposed as an identity example.
            schema_require(re.fullmatch(r'[A-Za-z0-9_.-]+', fields[0]) is not None
                    and any(c.isalpha() for c in fields[0]),
                    'RNA first field is not a narrow gene identity; contents not recorded', len(fields))
            schema_require(len(fields) == 1 or re.fullmatch(r'[A-Za-z0-9_.-]+', fields[1]) is not None
                    and any(c.isalpha() for c in fields[1]),
                    'RNA symbol field is not an identity; contents not recorded', len(fields))
            schema_require(len(fields) < 3 or fields[2] == 'Gene Expression',
                    'Unsupported RNA annotation/type; contents not recorded', len(fields))
            require(fields[0] not in identities, 'Duplicate source RNA feature identity')
            identities.add(fields[0])
            kind = fields[2] if len(fields) == 3 else 'type_not_declared'
            types[kind] += 1
            if re.fullmatch(r'ENSG\d+(?:\.\d+)?', fields[0]):
                require(kind in ('Gene Expression', 'type_not_declared'), 'Gene ID with unexpected feature type')
                genes.add(fields[0])
                gene_digest.update(('\t'.join(fields) + '\n').encode())
        else:
            # Three-field BED-like files repeat chromosome names; identity is the complete row.
            schema_require(len(fields) in (1, 3), 'ATAC schema rejected', len(fields))
            try:
                coordinate, form = peak_coordinate(fields)
            except ValueError as error:
                raise SchemaFailure(str(error), len(fields)) from None
            require(coordinate not in raw_coordinates, 'Duplicate source ATAC coordinates')
            identities.add(text)
            raw_coordinates.add(coordinate); formats.add(form)
            if text in fixed_ids:
                matches.add(text)
            chrom, start, end = coordinate
            if end > start and coordinate in fixed_bed:
                zero_matches.add(coordinate)
            shifted = (chrom, start - 1, end)
            if start >= 1 and end >= start and shifted in fixed_bed:
                one_matches.add(shifted)
        # Examples are eligible only after the narrow row schema was confirmed.
        rows += 1; digest.update(raw); fields_seen[len(fields)] += 1
        if len(examples) < 3:
            examples.append(text)
    require(rows > 0, 'Empty source feature axis')
    result = dict(rows=rows, decoded_axis_sha256=digest.hexdigest(),
        hash_codec='original UTF8 decompressed bytes including line endings', first_three=examples,
        field_counts=dict(fields_seen), original_strings_retained_in_compressed_member=True)
    if assay == 'snRNA-seq':
        result.update(feature_types=dict(types), gene_ids_matching_ENSG_syntax=len(genes),
            gene_id_symbol_type_ordered_sha256=gene_digest.hexdigest(),
            gene_hash_codec='original tab fields LF for each ENSG row; no version stripping',
            symbols_present=all(n >= 2 for n in fields_seen), feature_type_present=all(n == 3 for n in fields_seen),
            genetic_feature_identity_only=True, RNA_count_values_read=False)
    else:
        result.update(formats=sorted(formats), fixed_axis_rows=306706,
            exact_original_string_matches=len(matches), exact_original_string_missing=306706 - len(matches),
            coordinate_origin_verified=False,
            qualified_hypotheses=dict(source_zero_based_half_open_matches=len(zero_matches),
                source_one_based_closed_matches=len(one_matches)),
            coordinate_hypotheses_not_selected_by_best_overlap=True,
            missing_features_are_not_observed_zeros=True, raw_count_compatibility_claimed=False)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-script-sha256', required=True)
    args = parser.parse_args()
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Acquisition/identity parsing requires compute')
    require(sha(Path(__file__)) == args.expected_script_sha256, 'Reviewed helper source hash changed')
    guards = {str(PROPOSAL): PROPOSAL_SHA, str(IDENTITY): IDENTITY_SHA, str(GEOMETRY): GEOMETRY_SHA}
    output = REC / ('codex_gse223843_feature_barcode_ranges_' + os.environ['SLURM_JOB_ID'])
    owned, acquisition, phase, checked = False, None, 'fresh_output', {}
    try:
        output.mkdir(exist_ok=False); owned = True
        acquisition = Acquisition(output)
        phase = 'frozen_source_hashes'
        for name, expected in guards.items():
            path = Path(name)
            checked[name] = sha(path)
            require(checked[name] == expected, 'Frozen source hash changed: ' + name)
        proposal, identity = json.loads(PROPOSAL.read_text()), json.loads(IDENTITY.read_text())
        members = proposal['members']
        allowlist = {row['member_name'] for row in members}
        require(len(members) == len(allowlist) == 24
            and proposal['declared_archive_bytes'] == ARCHIVE_SIZE
            and proposal['official_archive_url'] == ARCHIVE_URL
            and proposal['official_filelist_sha256'] == FILELIST_SHA
            and sum(row['declared_compressed_bytes'] for row in members) == 16668757,
            'Frozen proposal acquisition identity differs')
        require(all(row['role'] in ('features', 'barcodes') and row['assay'] in ('snRNA-seq', 'snATAC-seq')
            and row['official_archive_url'] == ARCHIVE_URL
            and row['accession'] == identity['person_axis_by_exact_author_labels'][row['source_person_label']][row['assay']]
            and row['member_name'].startswith(row['accession'] + '_')
            and row['member_name'].endswith('_' + row['role'] + ('.tsv.gz' if row['assay'] == 'snRNA-seq' else '.txt.gz'))
            for row in members), 'Unexpected member role/assay/person/filename')
        for name, path in [('executed_source.py', Path(__file__)), ('executed_launcher.sbatch', HERE / 'run_codex_gse223843_feature_barcode_ranges.sbatch'),
                           ('source_proposal.json', PROPOSAL), ('source_identity.json', IDENTITY)]:
            with (output / name).open('xb') as target:
                target.write(path.read_bytes())
        phase = 'official_filelist_hash'
        file_sizes = acquisition.filelist()
        require(set(file_sizes) == set(identity['listed_member_filenames']) - {'GSE223843_RAW.tar'}
                and all(file_sizes[row['member_name']] == row['declared_compressed_bytes'] for row in members),
                'Live official filelist names/member bytes differ')
        phase = 'tar_metadata_ranges'
        offsets, tar_summary = scan_archive(acquisition, file_sizes, allowlist)
        write_json(output, 'tar_identity_offsets.json', dict(archive_url=ARCHIVE_URL,
            archive_bytes=ARCHIVE_SIZE, allowlisted_offsets=offsets, **tar_summary, **FLAGS))
        phase = 'allowlisted_gzip_members'
        receipts = {}
        for row in members:
            name, item = row['member_name'], offsets[row['member_name']]
            partial = output / (name + '.partial')
            receipt = acquisition.range(item['offset'], item['bytes'], destination=partial)
            require(partial.stat().st_size == row['declared_compressed_bytes']
                    and sha(partial) == receipt['sha256'], 'Member byte/hash receipt differs')
            final = output / name
            require(not final.exists(), 'Member output already exists')
            partial.rename(final)
            receipts[name] = dict(compressed_bytes=final.stat().st_size, sha256=receipt['sha256'],
                source_person_label=row['source_person_label'], assay=row['assay'], role=row['role'])
        phase = 'frozen_geometry_identity'
        fixed_ids, fixed_bed = frozen_geometry()
        phase = 'gzip_CRC_and_identity_formats'
        results, barcode_sets = {}, {}
        for row in members:
            path = output / row['member_name']
            phase = 'gzip_CRC_and_identity_formats/' + row['member_name']
            if row['role'] == 'barcodes':
                values, result = barcodes(path, acquisition)
                barcode_sets[(row['source_person_label'], row['assay'])] = values
            else:
                result = features(path, row['assay'], acquisition, fixed_ids, fixed_bed)
            results[row['member_name']] = dict(**result, gzip_CRC_verified=True)
        pairs = {}
        for person in identity['person_axis_by_exact_author_labels']:
            rna, atac = barcode_sets[(person, 'snRNA-seq')], barcode_sets[(person, 'snATAC-seq')]
            pairs[person] = dict(RNA_barcodes=len(rna), ATAC_barcodes=len(atac), shared_exact_barcodes=len(rna & atac),
                RNA_only_barcodes=len(rna - atac), ATAC_only_barcodes=len(atac - rna),
                exact_full_sets_equal=rna == atac, namespace=person,
                set_mismatch_does_not_admit_or_filter_count_rows=True)
        phase = 'completion_hashes_and_inclusive_bytes'
        require(all(sha(Path(name)) == expected for name, expected in guards.items()),
                'Frozen identity/fixture metadata changed during acquisition')
        summary = dict(completed=True, phase=phase, slurm_job_id=os.environ['SLURM_JOB_ID'],
            checked_source_sha256=checked, executed_source_sha256=sha(output / 'executed_source.py'),
            executed_launcher_sha256=sha(output / 'executed_launcher.sbatch'),
            archive_http_ranges_verified=True, archive_sha256_not_computed=True,
            member_receipts=receipts, identity_format_results=results, per_person_barcode_join=pairs,
            tar_scan=tar_summary, HTTP_range_receipts=acquisition.ranges,
            network_response_body_bytes=acquisition.network_bytes,
            tar_metadata_response_body_bytes=acquisition.metadata_bytes,
            decoded_identity_bytes=acquisition.decoded_bytes,
            network_cap=NETWORK_CAP, tar_metadata_cap=METADATA_CAP, saved_output_cap=SAVED_CAP,
            decoded_per_file_cap=DECODED_FILE_CAP, decoded_total_cap=DECODED_TOTAL_CAP,
            source_population='six pediatric author-person labels; independence unverified',
            fixed_axis=306706, coordinate_origin_unresolved=True,
            zero_fill_or_intersection_count_target_created=False, lineage_labels_not_read=True,
            source_exposure_and_cross_source_person_disjointness_verified=False,
            python=sys.version, gzip_module=gzip.__file__, zlib_runtime=zlib.ZLIB_RUNTIME_VERSION,
            saved_bytes_before_summary=saved_bytes(output), **FLAGS)
        write_json(output, 'summary.json', summary)
        require(saved_bytes(output) <= SAVED_CAP, 'Final inclusive output cap exceeded')
        print(json.dumps(dict(completed=True, output=str(output), saved_bytes=saved_bytes(output), **FLAGS)))
    except BaseException as error:
        if owned:
            failure = dict(completed=False, phase=phase, exception_type=type(error).__name__,
                message=str(error), traceback=traceback.format_exc(), checked_source_sha256=checked,
                partial_files_preserved=True, **FLAGS)
            if isinstance(error, SchemaFailure):
                failure['rejected_schema_field_count'] = error.field_count
                failure['rejected_schema_field_contents_saved'] = False
            if acquisition is not None:
                failure.update(network_response_body_bytes=acquisition.network_bytes,
                    tar_metadata_response_body_bytes=acquisition.metadata_bytes,
                    decoded_identity_bytes=acquisition.decoded_bytes, HTTP_range_receipts=acquisition.ranges)
            try:
                write_json(output, 'failure.json', failure)
            except Exception as receipt_error:
                print('Failure receipt could not be written: ' + str(receipt_error), file=sys.stderr)
        raise


if __name__ == '__main__':
    main()
