#!/usr/bin/env python3
"""Compute-only acquisition of six mixed GEX MTXs; qualify/aggregate GE rows only."""
import argparse
import csv
import gzip
import hashlib
import io
import json
import os
from pathlib import Path
import re
import shutil
import sys
import time
from urllib.request import HTTPRedirectHandler, Request, build_opener
import zlib

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = ROOT / 'scripts/analysis/in_biopsy_allelic/rec_codex'
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ANN = REC / 'codex_gse223843_rna_annotation_identity_22000649'
BC = REC / 'codex_gse223843_barcode_namespace_22000942'
ACQ = REC / 'codex_gse223843_feature_barcode_ranges_21999459'
IDENTITY = REC / 'codex_gse223843_identity_metadata/identity_metadata.json'
GENES = REC / 'codex_gse296875_aligned_native_pseudobulk_21999585/genes.tsv'
ARC_MATRIX = Path('/nfs/sw/easybuild/software/cellranger-arc-2.0.2/lib/python/atac_rna/matrix.py')
ARCHIVE_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223843/suppl/GSE223843_RAW.tar'
FILELIST_URL = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE223nnn/GSE223843/suppl/filelist.txt'
FILELIST_SHA = 'd0e3fe3012e57aa3d0364378d90357e400dbdfb6346a2ebd892d6d13a74c295a'
ARCHIVE_SIZE = 1066270720
PERSONS = ('Control1', 'Control3', 'Fontan1', 'Fontan2', 'Fontan3', 'Fontan4')
MATRICES = {
    'Control1': ('GSM6997741_GEX_Ctrl_matrix.mtx.gz', 55023947),
    'Control3': ('GSM6997744_GEX_Ctrl3_matrix.mtx.gz', 287164726),
    'Fontan1': ('GSM6997746_GEX_Fontan1_matrix.mtx.gz', 87519936),
    'Fontan2': ('GSM6997748_GEX_Fontan2_matrix.mtx.gz', 149823093),
    'Fontan3': ('GSM6997750_GEX_Fontan3_matrix.mtx.gz', 106057620),
    'Fontan4': ('GSM6997752_GEX_Fontan4_matrix.mtx.gz', 67537021),
}
AXIS_SHA = 'afbcc02143e68d1fbcda9631917f1344adfd3ac8efb752e2ebb7f044df934e60'
PINS = {
    IDENTITY: '172aa3c42134a2a4f329cdf24db937d92c70da0a1300827f32c6cf13b3754018',
    ANN / 'summary.json': 'da3f1310cfba1bc145ded47979a56e4134e2870baf09a25d553fc598024a7abe',
    BC / 'summary.json': '38561f9c1680926208bdd0e7c201fc4740131707ac52d109d4024236132ff768',
    ACQ / 'failure.json': 'ad5218429c3bd239a0cc46ae3e81b69b71480a901c1aa68fa87ee561fba2799a',
    ACQ / 'official_filelist.txt': FILELIST_SHA,
    ACQ / 'executed_source.py': '7f9845fd6ba4e756502b246fafd1500a55aa4b3987720e0569fc82c174224691',
    GENES: '2a33d4d339e2c6db695525ae2655af5156ccda008300bd992e306d48f01d38f0',
    ARC_MATRIX: '3795c9c494845e208c4f10ae9f53bb7e51199badec394221e17e3aa06caf47f8',
    REC / '20261001_codex_gse296875_donor_global_rna_atac_prespec.json':
        '88b9e4239b08cc1ccbe78b092e8eacd395778453a496fd8915b4acc9e2729ae1',
}
CAP = 1024**3                 # Network and inclusive saved bytes, independently.
DECODE_CAP = 8 * 1024**3       # All MTX decoded bytes, including discarded Peak tokens.
BITSET_CAP = 512 * 1024**2
RESERVE = 1024**2
LINE_CAP = 4096
TIME_CAP = 110 * 60
MAX_INT = 2**63 - 1
FLAGS = dict(ATAC_matrix_members_requested=False, ATAC_numeric_values_parsed=False,
    ATAC_counts_or_targets_or_scores_created=False, clinical_or_genotype_values_read=False,
    H3_or_allelic_or_protected_data_read=False, normalization_or_transform_applied=False,
    lineage_assignment=False, fitted_or_predicted=False, outcome_based_selection=False,
    source_person_disjointness_verified=False, receiving_measurement_admission=False,
    model_or_evaluation_admission=False)


class QualificationFailure(RuntimeError):
    """Fixed categorical reasons only: never include a source matrix line/value."""


def require(ok, reason):
    if not ok:
        raise QualificationFailure(reason)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    require(path.stat().st_size <= 1024**2, 'metadata_size_cap')
    return json.loads(path.read_text())


class Budget:
    def __init__(self, out, job):
        self.out = out
        self.logs = [REC / f'1592_gex_mtx_{job}.{s}' for s in ('out', 'err')]
        self.start = time.monotonic()
        self.network = 0
        self.decoded = 0
        self.phase = 'source_guards'
        self.ranges = []

    def saved(self):
        files = list(self.out.iterdir()) + self.logs
        require(not any(p.is_symlink() for p in files), 'output_symlink')
        return sum(p.stat().st_size for p in files if p.is_file())

    def check(self, extra=0):
        require(time.monotonic() - self.start < TIME_CAP, 'wall_clock_cap')
        require(self.saved() + extra <= CAP - RESERVE, 'inclusive_saved_cap')
        require(self.network <= CAP and self.decoded <= DECODE_CAP, 'network_or_decode_cap')

    def json(self, name, obj):
        body = (json.dumps(obj, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
        self.check(len(body))
        with (self.out / name).open('xb') as stream:
            stream.write(body)

    def progress(self, phase):
        self.phase = phase
        self.check()
        # Fixed phase names and byte accounting only, never a source entry.
        print(json.dumps({'phase': phase, 'network_bytes': self.network,
                          'decoded_bytes': self.decoded, 'elapsed_seconds': time.monotonic() - self.start}),
              flush=True)


class NoRedirect(HTTPRedirectHandler):
    def redirect_request(self, request, fp, code, message, headers, new_url):
        raise QualificationFailure('HTTP_redirect_forbidden')


class Acquisition:
    def __init__(self, budget, old_headers):
        self.budget = budget
        self.opener = build_opener(NoRedirect())
        self.old_headers = old_headers
        self.seen_headers = set()
        self.allowed_body_ranges = set()

    def filelist(self):
        request = Request(FILELIST_URL, headers={'Accept-Encoding': 'identity'})
        with self.opener.open(request, timeout=45) as response:
            require(response.status == 200 and response.geturl() == FILELIST_URL, 'filelist_HTTP_identity')
            require(response.headers.get('Content-Length') in (None, '2783'), 'filelist_length_changed')
            require(response.headers.get('Content-Encoding', 'identity').lower() == 'identity', 'filelist_encoding')
            raw = response.read(2784)
        self.budget.network += len(raw)
        require(len(raw) == 2783 and hashlib.sha256(raw).hexdigest() == FILELIST_SHA, 'filelist_SHA_changed')
        self.budget.check(len(raw))
        with (self.budget.out / 'official_filelist.txt').open('xb') as stream:
            stream.write(raw)
        reader = csv.DictReader(io.StringIO(raw.decode()), delimiter='\t')
        require(reader.fieldnames == ['#Archive/File', 'Name', 'Time', 'Size', 'Type'], 'filelist_schema')
        files, archives = {}, []
        for row in reader:
            require(None not in row and row['Size'].isdecimal(), 'filelist_row_schema')
            if row['#Archive/File'] == 'Archive':
                archives.append(row)
            else:
                require(row['#Archive/File'] == 'File' and row['Name'] not in files, 'filelist_duplicate_or_kind')
                files[row['Name']] = int(row['Size'])
        require(len(archives) == 1 and archives[0]['Name'] == 'GSE223843_RAW.tar'
                and int(archives[0]['Size']) == ARCHIVE_SIZE, 'archive_size_or_name_changed')
        return files

    def range(self, start, length, metadata=False, destination=None):
        key = (start, length)
        require(length > 0 and start >= 0 and start + length <= ARCHIVE_SIZE, 'range_outside_archive')
        require(key in (self.old_headers if metadata else self.allowed_body_ranges), 'range_not_allowlisted')
        require(metadata == (destination is None), 'range_destination_contract')
        self.budget.check(length if destination is not None else 0)
        require(self.budget.network + length <= CAP, 'network_cap_before_range')
        end = start + length - 1
        request = Request(ARCHIVE_URL, headers={'Range': f'bytes={start}-{end}',
            'Accept-Encoding': 'identity', 'User-Agent': 'bounded-GEX-qualification/1'})
        with self.opener.open(request, timeout=45) as response:
            # Check all response headers before touching a body; no HTTP200 fallback.
            require(response.status == 206 and response.geturl() == ARCHIVE_URL, 'exact_206_required')
            require(response.headers.get('Content-Range') == f'bytes {start}-{end}/{ARCHIVE_SIZE}', 'exact_content_range_required')
            require(response.headers.get('Content-Length') == str(length), 'exact_content_length_required')
            require(response.headers.get('Content-Encoding', 'identity').lower() == 'identity', 'range_encoding')
            entry = {'start': start, 'end': end, 'bytes_expected': length, 'bytes_read': 0,
                     'purpose': 'frozen_TAR_metadata' if metadata else 'allowlisted_mixed_GEX_MTX', 'status': 206}
            self.budget.ranges.append(entry)
            h, body = hashlib.sha256(), bytearray()
            target = destination.open('xb', buffering=0) if destination is not None else None
            try:
                remaining = length
                while remaining:
                    self.budget.check()
                    block = response.read(min(1024**2, remaining))
                    require(block and len(block) <= remaining, 'truncated_or_oversized_range')
                    remaining -= len(block)
                    self.budget.network += len(block)
                    entry['bytes_read'] += len(block)
                    h.update(block)
                    if target is None:
                        body.extend(block)
                    else:
                        self.budget.check(len(block))
                        target.write(block)
                # Declared Content-Length bounds this response. Never probe beyond the requested range.
                entry['sha256'] = h.hexdigest()
                if metadata:
                    require(entry['sha256'] == self.old_headers[key], 'frozen_TAR_header_SHA_changed')
                    self.seen_headers.add(key)
            finally:
                if target is not None:
                    target.close()
        return bytes(body) if metadata else entry


def tar_number(field):
    require(not field[0] & 128, 'TAR_base256_unsupported')
    value = field.strip(b'\x00 ')
    require(not value or re.fullmatch(b'[0-7]+', value), 'TAR_octal_schema')
    return int(value or b'0', 8)


def safe_path(raw):
    name = raw.decode('utf-8', errors='strict')
    while name.startswith('./'):
        name = name[2:]
    require(name and not name.startswith('/') and '\\' not in name
            and all(ord(c) >= 32 and ord(c) != 127 for c in name), 'TAR_path_schema')
    require(all(p not in ('', '.', '..') for p in name.rstrip('/').split('/')), 'TAR_path_components')
    return name


def pax_fields(raw):
    fields, position = {}, 0
    while position < len(raw):
        space = raw.find(b' ', position)
        require(space > position and raw[position:space].isdigit(), 'PAX_length_schema')
        length = int(raw[position:space])
        require(length > space - position + 2 and position + length <= len(raw), 'PAX_length_bound')
        record = raw[space + 1:position + length]
        require(record.endswith(b'\n') and b'=' in record, 'PAX_record_schema')
        key, value = record[:-1].decode('utf-8').split('=', 1)
        require(key not in fields and key in {'path', 'size', 'mtime', 'atime', 'ctime',
                'uid', 'gid', 'uname', 'gname', 'charset'}, 'PAX_key_unsupported')
        fields[key] = value
        position += length
    if 'path' in fields:
        safe_path(fields['path'].encode())
    require('size' not in fields or fields['size'].isdecimal(), 'PAX_size_schema')
    require(fields.get('charset', 'UTF-8') == 'UTF-8', 'PAX_charset')
    return fields


def scan_archive(acquisition, files):
    position, seen, offsets, pending = 0, set(), {}, {}
    allowlist = {name for name, _ in MATRICES.values()}
    while position + 1024 <= ARCHIVE_SIZE:
        header = acquisition.range(position, 512, metadata=True)
        if header == bytes(512):
            require(not pending and acquisition.range(position + 512, 512, metadata=True) == bytes(512), 'TAR_second_EOF')
            require(seen == set(files) and len(seen) == 36 and set(offsets) == allowlist, 'TAR_inventory_incomplete')
            require(acquisition.seen_headers == set(acquisition.old_headers), 'TAR_metadata_ranges_changed')
            require((ARCHIVE_SIZE - position) % 512 == 0, 'TAR_padding_alignment')
            acquisition.allowed_body_ranges = {(v['offset'], v['bytes']) for v in offsets.values()}
            return offsets
        require(tar_number(header[148:156]) == sum(header[:148]) + 8 * 32 + sum(header[156:]), 'TAR_checksum')
        require(header[257:263] in (b'ustar\x00', b'ustar ', bytes(6)), 'TAR_format_unsupported')
        name = header[:100].split(b'\x00', 1)[0]
        if header[257:263] == b'ustar\x00':
            prefix = header[345:500].split(b'\x00', 1)[0]
            if prefix:
                name = prefix + b'/' + name
        name = safe_path(name)
        size, kind, start = tar_number(header[124:136]), header[156:157], position + 512
        require(start + size <= ARCHIVE_SIZE, 'TAR_member_outside_archive')
        if kind == b'x':
            require(not pending and 0 < size <= 65536, 'PAX_chained_or_large')
            pending = pax_fields(acquisition.range(start, size, metadata=True))
        else:
            require(kind in (b'0', b'\x00', b'5'), 'TAR_link_sparse_extension_unsupported')
            name = safe_path(pending.get('path', name).encode())
            size = int(pending.get('size', size))
            pending = {}
            require(start + size <= ARCHIVE_SIZE, 'PAX_member_outside_archive')
            if kind == b'5':
                require(size == 0, 'TAR_nonempty_directory')
            else:
                require(name in files and name not in seen and size == files[name], 'TAR_member_name_or_size_changed')
                seen.add(name)
                if name in allowlist:
                    offsets[name] = {'offset': start, 'bytes': size}
        # Skip every regular body while discovering headers, including all six ATAC MTXs.
        position = start + ((size + 511) // 512) * 512
    raise QualificationFailure('TAR_EOF_missing')


def pin(path, expected, size=None):
    require(path.is_file() and not path.is_symlink(), 'pinned_source_missing_or_symlink')
    require(size is None or path.stat().st_size == size, 'pinned_source_bytes_changed')
    require(digest(path) == expected, 'pinned_source_SHA_changed')
    return {'path': str(path), 'bytes': path.stat().st_size, 'sha256': expected}


def fixed_genes():
    genes, axis = [], hashlib.sha256()
    with GENES.open(newline='') as stream:
        rows = csv.DictReader(stream, delimiter='\t')
        require(rows.fieldnames == ['gene_index', 'ensembl_id', 'native_gene_symbol'], 'native_gene_header')
        for i, row in enumerate(rows):
            require(None not in row and row['gene_index'] == str(i), 'native_gene_index')
            genes.append((row['ensembl_id'], row['native_gene_symbol']))
            axis.update((row['ensembl_id'] + '\0' + row['native_gene_symbol'] + '\n').encode())
    require(len(genes) == 36601 and axis.hexdigest() == AXIS_SHA, 'native_complete_gene_axis_changed')
    return genes


def contracts(inventory, identity, ann, bc, genes):
    result = {}
    require(tuple(identity['person_axis_by_exact_author_labels']) == PERSONS, 'author_person_order')
    require(ann['all_six_full_annotation_qualified'] and bc['barcode_namespace_inspection_completed']
            and bc['all_12_processed_files_full_CRC_and_identity_verified'], 'prior_annotation_or_barcode_failure')
    for person in PERSONS:
        annotation = ann['results'][person]
        require(annotation['validated_ordered_Gene_Expression_ID_symbol_sha256'] == AXIS_SHA
                and annotation['fixed36601_comparison']['ordered_exact_ID_and_name_axis']
                and annotation['gzip_CRC_passed'], 'prior_complete_gene_identity_unqualified')
        rows = annotation['full_source_feature_rows']
        require(rows == sum(annotation['source_feature_types'].values()) and
                set(annotation['source_feature_types']) == {'Gene Expression', 'Peaks'} and
                annotation['source_feature_types']['Gene Expression'] == 36601, 'mixed_feature_axis_schema')
        join_name = person + '_fixed36601_gene_identity_join.tsv'
        expected = ann['output_artifacts'][join_name]
        inventory[str(ANN / join_name)] = pin(ANN / join_name, expected['sha256'], expected['bytes'])
        row_to_gene, observed = [-1] * rows, 0
        with (ANN / join_name).open(newline='') as stream:
            table = csv.DictReader(stream, delimiter='\t')
            require(table.fieldnames == ['fixed_gene_index', 'fixed_gene_ID', 'fixed_native_name', 'ID_match_state',
                'source_Gene_Expression_indices', 'source_feature_indices', 'name_match_state'], 'gene_join_header')
            for i, row in enumerate(table):
                require(i < 36601 and None not in row and row['fixed_gene_index'] == str(i)
                        and (row['fixed_gene_ID'], row['fixed_native_name']) == genes[i]
                        and row['ID_match_state'] == row['name_match_state'] == 'exact'
                        and row['source_Gene_Expression_indices'] == str(i)
                        and row['source_feature_indices'].isdecimal(), 'gene_join_row')
                source_row = int(row['source_feature_indices'])
                require(0 <= source_row < rows and row_to_gene[source_row] == -1, 'duplicate_gene_join')
                row_to_gene[source_row] = i
                observed += 1
        require(observed == 36601, 'gene_join_incomplete')
        roster = bc['sample_paired_rosters'][person]
        require(roster['exact_existing_join_reproduced'] and roster['no_count_or_outcome_based_eligibility']
                and roster['roster_order'] == 'RNA native barcode file order', 'paired_roster_contract')
        columns = roster['RNA_barcodes']
        selected, previous = bytearray(columns), 0
        name = person + '_exact_source_paired_barcode_roster.tsv'
        expected = bc['output_artifacts_before_summary'][name]
        inventory[str(BC / name)] = pin(BC / name, expected['sha256'], expected['bytes'])
        n = 0
        with (BC / name).open(newline='') as stream:
            table = csv.DictReader(stream, delimiter='\t')
            require(table.fieldnames == ['author_person', 'RNA_accession', 'ATAC_accession', 'barcode_original',
                'RNA_source_row_1based', 'ATAC_source_row_1based', 'GEM_suffix_original', 'namespace_membership'], 'paired_roster_header')
            for row in table:
                require(None not in row and row['author_person'] == person
                        and row['RNA_accession'] == identity['person_axis_by_exact_author_labels'][person]['snRNA-seq']
                        and row['ATAC_accession'] == identity['person_axis_by_exact_author_labels'][person]['snATAC-seq']
                        and row['RNA_source_row_1based'].isdecimal(), 'paired_roster_row')
                col = int(row['RNA_source_row_1based'])
                require(previous < col <= columns and not selected[col - 1], 'paired_roster_order_or_duplicate')
                selected[col - 1] = 1
                previous = col
                n += 1
        require(n == roster['shared_exact_barcodes'] and n > 0, 'paired_roster_count_changed')
        result[person] = {'row_to_gene': row_to_gene, 'selected_columns': selected,
                          'feature_rows': rows, 'source_RNA_columns': columns, 'eligible_cells': n,
                          'RNA_accession': roster['RNA_accession']}
    return result


def matrix_lines(path, budget, state):
    with gzip.open(path, 'rb') as stream:
        while True:
            raw = stream.readline(LINE_CAP + 1)
            if not raw:
                break  # All gzip members, CRCs, trailers and full EOF must complete.
            budget.decoded += len(raw)
            state['decoded_bytes'] += len(raw)
            require(len(raw) <= LINE_CAP and budget.decoded <= DECODE_CAP, 'matrix_decode_or_line_cap')
            state['decoded_hash'].update(raw)
            yield raw


def qualify_matrix(path, contract, budget):
    state = {'decoded_bytes': 0, 'decoded_hash': hashlib.sha256()}
    lines = iter(matrix_lines(path, budget, state))
    header = next(lines, None)
    require(header is not None and header.rstrip(b'\r\n') == b'%%MatrixMarket matrix coordinate integer general', 'MTX_integer_general_required')
    dimension = None
    for raw in lines:
        if raw.startswith(b'%'):
            continue
        dimension = raw.split()
        break
    require(dimension is not None and len(dimension) == 3
            and all(re.fullmatch(rb'[0-9]{1,19}', f) for f in dimension), 'MTX_dimension_schema')
    rows, cols, declared = map(int, dimension)
    require(rows == contract['feature_rows'] and rows >= 1 and cols == contract['source_RNA_columns']
            and cols >= 1 and 0 < declared <= rows * cols, 'MTX_dimensions_or_entry_count')
    bit_bytes = (rows * cols + 7) // 8
    require(bit_bytes <= BITSET_CAP, 'duplicate_bitset_memory_cap')
    seen = bytearray(bit_bytes)
    counts, totals_by_column = [0] * 36601, [0] * cols
    observed, GE_entries, GE_explicit_zeros = 0, 0, 0
    library = 0
    row_map, selected = contract['row_to_gene'], contract['selected_columns']
    for raw in lines:
        # No Peak value is decoded, cast, regex-validated numerically, accumulated or printed.
        fields = raw.split()
        require(len(fields) == 3 and re.fullmatch(rb'[0-9]{1,19}', fields[0])
                and re.fullmatch(rb'[0-9]{1,19}', fields[1]), 'MTX_coordinate_schema')
        row, col = int(fields[0]), int(fields[1])
        require(1 <= row <= rows and 1 <= col <= cols, 'MTX_coordinate_bounds')
        linear = (col - 1) * rows + row - 1
        slot, mask = linear >> 3, 1 << (linear & 7)
        require(not seen[slot] & mask, 'MTX_duplicate_coordinate_strict_abort')
        seen[slot] |= mask
        observed += 1
        require(observed <= declared, 'MTX_extra_entries')
        gene = row_map[row - 1]
        if gene >= 0:
            require(re.fullmatch(rb'[0-9]{1,19}', fields[2]), 'GE_nonnegative_integer_schema')
            value = int(fields[2])
            require(value <= MAX_INT, 'GE_integer_bound')
            GE_entries += 1
            GE_explicit_zeros += value == 0
            if selected[col - 1]:
                counts[gene] += value
                totals_by_column[col - 1] += value
                library += value
                require(counts[gene] <= MAX_INT and totals_by_column[col - 1] <= MAX_INT
                        and library <= MAX_INT, 'GE_pseudobulk_integer_bound')
        # Peak fields[2] are discarded unchanged. Counts from unselected GE columns are checked only.
        if observed % 262144 == 0:
            budget.check()
    require(observed == declared, 'MTX_declared_entry_count_not_met')
    require(sum(counts) == library == sum(totals_by_column), 'GE_library_sum_invariant')
    budget.check()
    return counts, {'MTX_feature_rows': rows, 'MTX_RNA_columns': cols,
        'MTX_declared_entries': declared, 'MTX_observed_entries': observed,
        'MTX_structure_and_GE_integer_values_qualified': True,
        'Peak_numeric_schema_verified': False, 'Peak_value_tokens_discarded_without_numeric_parse': True,
        'mixed_GEX_compressed_and_decompressed_transport_contains_Peak_bytes': True,
        'full_gzip_CRC_and_EOF_verified': True, 'decoded_bytes': state['decoded_bytes'],
        'decoded_mixed_file_SHA256_not_a_Peak_count_summary': state['decoded_hash'].hexdigest(),
        'GE_stored_entries': GE_entries, 'GE_explicit_zero_entries_diagnostic_only': GE_explicit_zeros,
        'GE_explicit_zeros_did_not_change_eligibility': True,
        'duplicate_coordinate_policy': 'strict abort across every matrix coordinate; no additive fallback',
        'duplicate_bitset_bytes': bit_bytes, 'eligible_cells': contract['eligible_cells'],
        'eligible_cells_with_zero_GE_UMIs_diagnostic_only': sum(selected[i] and n == 0 for i, n in enumerate(totals_by_column)),
        'raw_UMI_library_sum_complete36601': library, 'positive_library_denominator': library > 0,
        'no_cell_or_gene_dropped_from_fixed_roster_or_axis': True,
        'integer_structure_is_not_independent_proof_of_native_UMI_measurement': True}


def save_counts(person, counts, genes, budget):
    raw = io.StringIO(newline='')
    writer = csv.writer(raw, delimiter='\t', lineterminator='\n')
    writer.writerow(['gene_index', 'ensembl_id', 'native_gene_symbol', 'raw_UMI_sum'])
    for i, (count, gene) in enumerate(zip(counts, genes)):
        writer.writerow([i, *gene, count])
    body = raw.getvalue().encode()
    budget.check(len(body))
    name = person + '_raw_complete36601_Gene_Expression_UMIs.tsv'
    with (budget.out / name).open('xb') as stream:
        stream.write(body)
    return {'path': str(budget.out / name), 'bytes': len(body), 'sha256': digest(budget.out / name)}


def focused_selfchecks(budget):
    """Tiny synthetic fixtures on compute only; exercise the same GE-only parser."""
    contract = {'row_to_gene': [0, 1, -1], 'selected_columns': bytearray([1, 0]),
                'feature_rows': 3, 'source_RNA_columns': 2, 'eligible_cells': 1}
    prefix = b'%%MatrixMarket matrix coordinate integer general\n% SYNTHETIC ONLY\n'
    bodies = {
        'selfcheck_synthetic_mixed.mtx.gz': prefix + b'3 2 4\n1 1 2\n2 1 3\n1 2 10\n3 1 not_numeric\n',
        'selfcheck_synthetic_explicit_zero.mtx.gz': prefix + b'3 2 4\n1 1 0\n2 1 3\n1 2 10\n3 1 not_numeric\n',
        'selfcheck_synthetic_duplicate.mtx.gz': prefix + b'3 2 2\n1 1 2\n1 1 2\n',
    }
    fixtures = {}
    for name, raw in bodies.items():
        compressed = gzip.compress(raw, mtime=0)
        budget.check(len(compressed))
        path = budget.out / name
        with path.open('xb') as stream:
            stream.write(compressed)
        fixtures[name] = {'bytes': len(compressed), 'sha256': digest(path), 'synthetic_only': True}
    truncated_name = 'selfcheck_synthetic_truncated.mtx.gz'
    truncated = gzip.compress(bodies['selfcheck_synthetic_mixed.mtx.gz'], mtime=0)[:-4]
    budget.check(len(truncated))
    with (budget.out / truncated_name).open('xb') as stream:
        stream.write(truncated)
    fixtures[truncated_name] = {'bytes': len(truncated), 'sha256': digest(budget.out / truncated_name),
                                'synthetic_only': True, 'gzip_trailer_deliberately_truncated': True}
    results = {}
    counts, receipt = qualify_matrix(budget.out / 'selfcheck_synthetic_mixed.mtx.gz', contract, budget)
    require(counts[:2] == [2, 3] and sum(counts) == 5
            and receipt['raw_UMI_library_sum_complete36601'] == 5
            and receipt['full_gzip_CRC_and_EOF_verified']
            and receipt['Peak_value_tokens_discarded_without_numeric_parse'], 'selfcheck_mixed_GE_selection_failed')
    results['mixed_GE_selection'] = {'passed': True, 'synthetic_first_two_GE_totals': [2, 3],
        'synthetic_complete_library_sum': 5, 'unselected_column_excluded': True,
        'Peak_value_not_numeric_parsed_or_printed': True, 'full_gzip_CRC_and_EOF_verified': True}
    counts, receipt = qualify_matrix(budget.out / 'selfcheck_synthetic_explicit_zero.mtx.gz', contract, budget)
    require(counts[:2] == [0, 3] and sum(counts) == 3
            and receipt['GE_explicit_zero_entries_diagnostic_only'] == 1
            and receipt['eligible_cells'] == 1
            and receipt['no_cell_or_gene_dropped_from_fixed_roster_or_axis'], 'selfcheck_explicit_zero_failed')
    results['explicit_GE_zero'] = {'passed': True, 'synthetic_eligible_cells_unchanged': 1,
                                    'synthetic_explicit_zero_records': 1}
    try:
        qualify_matrix(budget.out / 'selfcheck_synthetic_duplicate.mtx.gz', contract, budget)
    except QualificationFailure as exc:
        require(str(exc) == 'MTX_duplicate_coordinate_strict_abort', 'selfcheck_duplicate_wrong_failure')
        results['duplicate_coordinate'] = {'passed': True, 'categorical_refusal': str(exc)}
    else:
        raise QualificationFailure('selfcheck_duplicate_not_refused')
    try:
        qualify_matrix(budget.out / truncated_name, contract, budget)
    except (EOFError, gzip.BadGzipFile) as exc:
        results['truncated_gzip'] = {'passed': True, 'categorical_refusal': type(exc).__name__,
                                      'full_EOF_success_not_emitted': True}
    else:
        raise QualificationFailure('selfcheck_truncated_gzip_not_refused')
    report = {'all_four_focused_selfchecks_passed': True, 'synthetic_only': True,
        'completed_before_any_network_request': budget.network == 0,
        'fixtures': fixtures, 'results': results,
        'real_source_MTX_or_counts_read': False, 'fixture_and_receipt_byte_cap': 16 * 1024}
    body = (json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
    require(sum(v['bytes'] for v in fixtures.values()) + len(body) <= 16 * 1024, 'selfcheck_fixture_receipt_cap')
    require(budget.network == 0, 'selfchecks_must_precede_network')
    budget.json('focused_synthetic_selfchecks.json', report)
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-script-sha256', required=True)
    args = parser.parse_args()
    job = os.environ.get('SLURM_JOB_ID', '')
    require(re.fullmatch(r'[1-9][0-9]*', job) is not None, 'compute_SLURM_required')
    require(not os.environ.get('SLURM_ARRAY_TASK_ID'), 'array_not_allowed')
    source = Path(__file__).resolve()
    require(digest(source) == args.expected_script_sha256, 'reviewed_script_SHA_changed')
    out = REC / ('codex_gse223843_gex_matrix_qualification_' + job)
    owned, budget, inventory, results = False, None, {}, {}
    try:
        out.mkdir(mode=0o700, exist_ok=False)
        owned = True
        budget = Budget(out, job)
        inventory = {str(p): pin(p, h) for p, h in PINS.items()}
        inventory[str(source)] = pin(source, args.expected_script_sha256)
        launcher = HERE / 'run_codex_gse223843_gex_matrix_qualification.sbatch'
        inventory[str(launcher)] = pin(launcher, digest(launcher))
        identity, ann, bc = read_json(IDENTITY), read_json(ANN / 'summary.json'), read_json(BC / 'summary.json')
        genes = fixed_genes()
        fixed = contracts(inventory, identity, ann, bc, genes)
        old = read_json(ACQ / 'failure.json')
        headers = {}
        for item in old['HTTP_range_receipts']:
            if item['purpose'] != 'TAR_metadata':
                continue
            key = (item['start'], item['bytes_expected'])
            require(key not in headers and item['status'] == 206
                    and item['bytes_read'] == item['bytes_expected']
                    and item['end'] + 1 == item['start'] + item['bytes_expected'], 'frozen_TAR_receipt_schema')
            headers[key] = item['sha256']
        require(len(headers) == 38, 'frozen_36_TAR_headers_and_two_EOF_receipts')
        for name, p in [('executed_source.py', source), ('executed_launcher.sbatch', launcher)]:
            budget.check(p.stat().st_size)
            shutil.copyfile(p, out / name)
        budget.json('protocol.json', {'source_inventory': inventory, 'author_person_order': PERSONS,
            'allowlisted_mixed_GEX_MTX_members': MATRICES, 'compressed_member_bytes_total': 753126343,
            'official_archive_bytes': ARCHIVE_SIZE, 'official_archive_url': ARCHIVE_URL,
            'RNA_axis_rows': 36601, 'RNA_ordered_ID_name_sha256': AXIS_SHA,
            'eligibility': 'exact same-author-person common roster from 22000942, fixed before values',
            'duplicate_policy': 'strict abort; packed bitset per matrix, no order assumption',
            'GE_values': 'all GE columns nonnegative integers; aggregate only fixed eligible columns',
            'Peak_values': 'transport unavoidable in mixed GEX members; no numeric parse or aggregate',
            'RNA_units_method_source': 'https://www.10xgenomics.com/support/software/cell-ranger-arc/2.0/analysis/outputs/feature-barcode-matrices',
            'source_processing': 'https://pmc.ncbi.nlm.nih.gov/articles/PMC11103255/',
            'network_cap': CAP, 'inclusive_saved_cap': CAP, 'decoded_cap': DECODE_CAP,
            'duplicate_bitset_cap': BITSET_CAP, 'internal_seconds_cap': TIME_CAP, **FLAGS})
        budget.progress('focused_synthetic_selfchecks_before_network')
        selfchecks = focused_selfchecks(budget)
        acquisition = Acquisition(budget, headers)
        budget.progress('official_filelist_and_unchanged_TAR_headers')
        files = acquisition.filelist()
        require(set(files) == set(identity['listed_member_filenames']) - {'GSE223843_RAW.tar'}, 'official_member_names_changed')
        require(sum(n for _, n in MATRICES.values()) == 753126343, 'fixed_GEX_total_bytes')
        for person, (name, size) in MATRICES.items():
            require(name in files and files[name] == size and '_GEX_' in name and '_ATAC_' not in name
                    and name.startswith(fixed[person]['RNA_accession'] + '_'), 'GEX_member_identity')
        offsets = scan_archive(acquisition, files)
        budget.json('allowlisted_GEX_offsets.json', offsets)
        for person in PERSONS:
            name, size = MATRICES[person]
            budget.progress('acquire_mixed_GEX_MTX/' + person)
            partial = out / (name + '.partial')
            item = offsets[name]
            receipt = acquisition.range(item['offset'], item['bytes'], destination=partial)
            require(partial.stat().st_size == size and digest(partial) == receipt['sha256'], 'acquired_member_bytes_or_SHA')
            # Keep the original .partial name permanently; CRC qualification is in the receipt, not a rename.
            budget.json(person + '_compressed_range_receipt.json', receipt)
            budget.progress('full_MTX_structure_GE_values_and_raw_pseudobulk/' + person)
            counts, result = qualify_matrix(partial, fixed[person], budget)
            result.update({'author_person': person, 'RNA_accession': fixed[person]['RNA_accession'],
                           'compressed_source': receipt, 'raw_counts': save_counts(person, counts, genes, budget), **FLAGS})
            budget.json(person + '_raw_GE_qualification.json', result)
            results[person] = result
        budget.progress('final_source_identity_and_saved_accounting')
        after = {p: pin(Path(p), v['sha256'], v['bytes']) for p, v in inventory.items()}
        require(after == inventory, 'source_inventory_changed_during_job')
        require(len(results) == 6, 'all_six_required')
        budget.json('summary.json', {'all_six_raw_GE_qualified': True, 'SLURM_job_id': job,
            'focused_synthetic_selfchecks': selfchecks,
            'source_inventory_before_after_identical': True, 'source_inventory': inventory,
            'results': results, 'HTTP_range_receipts': budget.ranges,
            'network_response_body_bytes': budget.network, 'decoded_mixed_MTX_bytes': budget.decoded,
            'inclusive_saved_bytes_before_summary': budget.saved(),
            'network_cap': CAP, 'inclusive_saved_cap': CAP, 'decoded_cap': DECODE_CAP,
            'elapsed_seconds': time.monotonic() - budget.start, 'Python_version': sys.version,
            'zlib_runtime': zlib.ZLIB_RUNTIME_VERSION, 'biological_author_people': 6,
            'RNA_count_structural_qualification_not_measurement_or_evaluation_admission': True,
            'same_nucleus_and_person_independence_limits_remain': True, **FLAGS})
        budget.check()
        print(json.dumps({'completed': True, 'summary': str(out / 'summary.json'),
                          'summary_sha256': digest(out / 'summary.json'), 'inclusive_saved_bytes': budget.saved(), **FLAGS}), flush=True)
        budget.check()
    except Exception as exc:
        # Retain every partial, qualified earlier person's raw output and fixed categorical failure.
        failure = {'all_six_raw_GE_qualified': False, 'phase': budget.phase if budget else 'fresh_output',
            'exception_type': type(exc).__name__, 'categorical_reason': str(exc) if isinstance(exc, QualificationFailure) else type(exc).__name__,
            'partials_preserved': True, 'completed_people': list(results), **FLAGS}
        if budget:
            failure.update({'network_response_body_bytes': budget.network, 'decoded_mixed_MTX_bytes': budget.decoded,
                            'HTTP_range_receipts': budget.ranges})
        if owned:
            with (out / 'failure.json').open('x') as stream:
                json.dump(failure, stream, indent=2, allow_nan=False)
                stream.write('\n')
        print(json.dumps({'completed': False, 'phase': failure['phase'], 'categorical_reason': failure['categorical_reason'], **FLAGS}), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
