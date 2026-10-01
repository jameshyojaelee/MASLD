#!/usr/bin/env python3
"""Offline ARC six-field annotation identity only; no count or barcode inputs."""
import argparse
from collections import Counter, defaultdict
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import sys
import time
import zlib

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
ACQUIRED = REC / 'codex_gse223843_feature_barcode_ranges_21999459'
PROPOSAL = REC / 'codex_gse223843_geometry_acquisition_proposal.json'
IDENTITY = REC / 'codex_gse223843_identity_metadata/identity_metadata.json'
GENES = REC / 'codex_gse296875_aligned_native_pseudobulk_21999585/genes.tsv'
LOCK = ROOT / 'Analysis/SingleCell/continual_integration/reference/gse296875_source_lock_v44.json'
ARC_SOURCE = Path('/nfs/sw/easybuild/software/cellranger-arc-2.0.2/lib/python/atac_rna/feature_ref.py')
LAUNCHER = HERE / 'run_codex_gse223843_rna_annotation_identity.sbatch'
GUARDS = {
    PROPOSAL: '35107ad4b17a7463bf28831191d7327ea895e28776b652f85d58ec20360cd7e8',
    IDENTITY: '172aa3c42134a2a4f329cdf24db937d92c70da0a1300827f32c6cf13b3754018',
    ACQUIRED / 'failure.json': 'ad5218429c3bd239a0cc46ae3e81b69b71480a901c1aa68fa87ee561fba2799a',
    ACQUIRED / 'tar_identity_offsets.json': '8e0467bcb07978e943cbbe69046fc359e4860677b61511a39bba1569aaf4f878',
    ACQUIRED / 'executed_source.py': '7f9845fd6ba4e756502b246fafd1500a55aa4b3987720e0569fc82c174224691',
    GENES: '2a33d4d339e2c6db695525ae2655af5156ccda008300bd992e306d48f01d38f0',
    LOCK: '6be686b4187309f038097cf73d3f2a8dd34e72d666bc92733a53e99a1ccbb620',
    ARC_SOURCE: '1b80c8ce8675e1fe0c079f0167fc4209a84712f0bcf28162465251ac33c2ac51',
}
NATIVE_DIGEST = 'afbcc02143e68d1fbcda9631917f1344adfd3ac8efb752e2ebb7f044df934e60'
PERSONS = ('Control1', 'Control3', 'Fontan1', 'Fontan2', 'Fontan3', 'Fontan4')
PINNED_MEMBERS = {
    'GSM6997741_GEX_Ctrl_features.tsv.gz': (1034803, '96c0f53ba3d3bb979ca9580f269503128c27163387b2e24d94b3f0295dcb0c38'),
    'GSM6997744_GEX_Ctrl3_features.tsv.gz': (2317012, '1eecc3696db023bc8a014ef650ac73b25d5fb7ff9ecc85dff4ccdd49870cf430'),
    'GSM6997746_GEX_Fontan1_features.tsv.gz': (1453275, 'cf6e3ec34def8b956d103b163d6ba508c773e11dd826928614a7f0410e973a36'),
    'GSM6997748_GEX_Fontan2_features.tsv.gz': (1565546, '394a841689d1f36f160e8ed68db2ad8ced93fcdcbf4c5857c0bf9d2209ef3220'),
    'GSM6997750_GEX_Fontan3_features.tsv.gz': (1354139, '0d977c633b3ff1f846395dbf36db25c0f6ff59d86207b1fea4388f7414825f0c'),
    'GSM6997752_GEX_Fontan4_features.tsv.gz': (1097809, '400ac1358bf181a7701244ef624a4b317f13507aa1ffba03ed7d87b5304ed40d'),
}
CAP, FILE_DECODE_CAP, LINE_CAP, ROW_CAP = 256 * 1024**2, 64 * 1024**2, 4096, 2000000
FLAGS = dict(network_requests=False, counts_or_matrices_read=False, barcodes_read=False,
    raw_sequence_read=False, clinical_values_read=False, genotypes_read=False, H3_values_read=False,
    outcomes_or_model_scores_read=False, symbol_collapse=False, count_rows_created=False,
    normalization_redefined=False, orientation_or_coordinate_hypothesis_selected=False,
    feature_or_lineage_selection=False, model_or_evaluation_admission=False)


class QualificationFailure(Exception):
    pass


def require(condition, reason):
    if not condition:
        raise QualificationFailure(reason)


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(65536), b''):
            h.update(block)
    return h.hexdigest()


def stable_id(gene):
    return gene.split('.', 1)[0]


def annotation_name(value):
    # ARC writes source FeatureDef.name, not a punctuation-restricted symbol.
    return bool(value) and len(value.encode('utf-8')) <= 1024 and value.isprintable()


def typed_row(fields):
    """Return known annotation category or a fixed diagnostic, never field values."""
    if len(fields) != 6:
        return None, 'unsupported_field_count'
    identity, name, kind, chrom, start, end = fields
    if kind not in ('Gene Expression', 'Peaks'):
        return None, 'unknown_feature_type'
    if not annotation_name(name):
        return None, 'unsupported_annotation_name'
    if kind == 'Gene Expression' and not re.fullmatch(r'ENSG[0-9]{11}(?:\.[0-9]+)?', identity):
        return None, 'unsupported_human_Ensembl_gene_ID'
    if kind == 'Gene Expression' and (chrom, start, end) == ('', '-1', '-1'):
        return 'source_missing_TSS_annotation', None
    if not re.fullmatch(r'[A-Za-z0-9_.-]{1,128}', chrom):
        return None, 'unsupported_annotation_contig'
    if not re.fullmatch(r'[0-9]{1,19}', start) or not re.fullmatch(r'[0-9]{1,19}', end):
        return None, 'unsupported_BED_coordinate_lexical_type'
    a, b = int(start), int(end)
    if not (0 <= a < b <= 2**63 - 1):
        return None, 'invalid_BED_half_open_interval'
    if kind == 'Peaks':
        match = re.fullmatch(r'([A-Za-z0-9_.-]{1,128}):([0-9]{1,19})-([0-9]{1,19})', identity)
        if match is None or (match[1], int(match[2]), int(match[3])) != (chrom, a, b):
            return None, 'peak_ID_interval_annotation_mismatch'
        return 'source_peak_BED_interval', None
    return 'source_TSS_BED_interval', None


class Run:
    def __init__(self, output, job):
        self.output, self.job = output, job
        self.phase, self.started = 'source_guards', time.monotonic()
        self.decoded = 0
        self.results = {}
        self.written = 0
        self.next_filesystem_check = 1024**2

    def check_time(self):
        require(time.monotonic() - self.started < 3300, 'overall_time_cap')

    def sizes(self):
        self.check_time()
        paths = list(self.output.iterdir())
        require(not any(p.is_symlink() for p in paths), 'output_symlink')
        n = max(self.written, sum(p.stat().st_size for p in paths if p.is_file()))
        for suffix in ('out', 'err'):
            log = REC / f'1592_rna_annotation_{self.job}.{suffix}'
            require(not log.is_symlink(), 'log_symlink')
            if log.is_file():
                n += log.stat().st_size
        require(n <= CAP, 'inclusive_saved_byte_cap')
        return n

    def write(self, path, raw):
        self.check_time()
        # Reserve 1 MiB for tiny SLURM logs; tracked writes include buffered TSVs.
        require(self.written + len(raw) <= CAP - 1024**2, 'saved_cap_before_write')
        path.write(raw)
        self.written += len(raw)
        if self.written >= self.next_filesystem_check:
            self.sizes()
            self.next_filesystem_check = self.written + 1024**2

    def json(self, name, value):
        raw = (json.dumps(value, indent=2, allow_nan=False) + '\n').encode()
        with (self.output / name).open('xb') as target:
            self.write(target, raw)

    def line(self, target, values):
        # All fields are validated annotation strings or metadata integers.
        raw = ('\t'.join(str(x) for x in values) + '\n').encode('utf-8')
        self.write(target, raw)

    def source_guards(self, script_sha):
        inventory = {}
        for path, expected in {**GUARDS, Path(__file__).resolve(): script_sha}.items():
            self.check_time()
            require(not path.is_symlink() and digest(path) == expected, 'pinned_source_SHA_mismatch')
            inventory[str(path)] = dict(bytes=path.stat().st_size, sha256=expected)
        inventory[str(LAUNCHER)] = dict(bytes=LAUNCHER.stat().st_size, sha256=digest(LAUNCHER))
        return inventory

    def fixed_axis(self):
        lock = json.loads(LOCK.read_text())
        require(lock['raw_rna_registry']['rna_feature_id_and_symbol_sha256'] == NATIVE_DIGEST,
                'frozen_native_axis_lock_mismatch')
        rows, ordered = [], hashlib.sha256()
        with GENES.open('rb') as source:
            header = source.readline(LINE_CAP + 1)
            self.decoded += len(header)
            require(header == b'gene_index\tensembl_id\tnative_gene_symbol\r\n'
                    or header == b'gene_index\tensembl_id\tnative_gene_symbol\n', 'fixed_axis_header')
            for index, raw in enumerate(source):
                self.decoded += len(raw)
                require(len(raw) <= LINE_CAP and self.decoded <= CAP, 'fixed_axis_decode_cap')
                fields = raw.rstrip(b'\r\n').decode('utf-8', errors='strict').split('\t')
                require(len(fields) == 3 and fields[0] == str(index)
                        and re.fullmatch(r'ENSG[0-9]{11}(?:\.[0-9]+)?', fields[1])
                        and annotation_name(fields[2]), 'fixed_native_gene_identity_schema')
                rows.append((fields[1], fields[2]))
                ordered.update((fields[1] + '\0' + fields[2] + '\n').encode())
                self.check_time()
        require(len(rows) == len({g for g, _ in rows}) == 36601
                and ordered.hexdigest() == NATIVE_DIGEST, 'fixed_native_ordered_axis')
        return rows

    def compressed_roster(self):
        proposal = json.loads(PROPOSAL.read_text()); identity = json.loads(IDENTITY.read_text())
        failure = json.loads((ACQUIRED / 'failure.json').read_text())
        offsets = json.loads((ACQUIRED / 'tar_identity_offsets.json').read_text())['allowlisted_offsets']
        require(failure['count_payload_inspected'] is False, 'original_acquisition_scope')
        members = [r for r in proposal['members'] if r['assay'] == 'snRNA-seq' and r['role'] == 'features']
        require(len(members) == 6 and {r['member_name'] for r in members} == set(PINNED_MEMBERS)
                and {r['source_person_label'] for r in members} == set(PERSONS), 'six_gene_feature_members')
        selected = {}
        for row in members:
            name = row['member_name']; path = ACQUIRED / name
            size, sha = PINNED_MEMBERS[name]; offset = offsets[name]
            require(path.parent == ACQUIRED and Path(name).name == name and not path.is_symlink(),
                    'compressed_member_path')
            require(row['accession'] == identity['person_axis_by_exact_author_labels'][row['source_person_label']]['snRNA-seq']
                    and name in identity['listed_member_filenames'], 'source_person_assay_identity')
            ranges = [r for r in failure['HTTP_range_receipts'] if r['purpose'] == 'allowlisted_identity_member'
                      and r['start'] == offset['offset']]
            require(len(ranges) == 1, 'original_range_identity_not_unique')
            receipt = ranges[0]
            require(receipt['status'] == 206 and receipt['bytes_read'] == receipt['bytes_expected']
                    == offset['bytes'] == row['declared_compressed_bytes'] == path.stat().st_size == size
                    and receipt['end'] == receipt['start'] + size - 1
                    and receipt['sha256'] == sha == digest(path), 'compressed_member_range_size_SHA')
            selected[row['source_person_label']] = dict(**row, path=str(path), sha256=sha,
                original_HTTP_range_receipt=receipt)
        return selected

    def decode_member(self, person, member):
        """Stream all feature rows; unsupported schemas remain categorical, no contents."""
        path = Path(member['path']); self.phase = 'six_field_annotation/' + person
        source_genes, all_ids, categories = [], Counter(), Counter()
        failures, widths, annotation_states = Counter(), Counter(), Counter()
        unknown_hashes = Counter(); rows = decoded = 0; crc = False
        raw_digest, axis_digest, gene_digest = hashlib.sha256(), hashlib.sha256(), hashlib.sha256()
        axis_path = self.output / f'{person}_source_annotation_axis.tsv'
        bad_path = self.output / f'{person}_unsupported_rows.tsv'
        with axis_path.open('xb', buffering=1024**2) as axis, bad_path.open('xb', buffering=1024**2) as bad:
            self.line(axis, ['source_feature_index', 'feature_id', 'source_feature_name', 'feature_type',
                             'chromosome', 'source_BED_start', 'source_BED_end', 'annotation_state'])
            self.line(bad, ['source_feature_index', 'field_count', 'reason', 'feature_type_SHA256'])
            try:
                with gzip.open(path, 'rb') as source:
                    while True:
                        raw = source.readline(LINE_CAP + 1)
                        if not raw:
                            crc = True; break
                        rows += 1; decoded += len(raw); self.decoded += len(raw)
                        require(rows <= ROW_CAP and len(raw) <= LINE_CAP and decoded <= FILE_DECODE_CAP
                                and self.decoded <= CAP, 'annotation_decode_limit')
                        raw_digest.update(raw)
                        try:
                            fields = raw.rstrip(b'\r\n').decode('utf-8', errors='strict').split('\t')
                        except UnicodeDecodeError:
                            failures['non_UTF8_annotation'] += 1
                            self.line(bad, [rows - 1, 'unverified', 'non_UTF8_annotation', ''])
                            break  # Unsupported payload stays closed; no remainder scan.
                        widths[len(fields)] += 1
                        kind = fields[2] if len(fields) == 6 else None
                        if kind in ('Gene Expression', 'Peaks'):
                            categories[kind] += 1
                            type_sha = ''
                        elif kind is not None:
                            categories['unknown'] += 1
                            type_sha = hashlib.sha256(kind.encode('utf-8')).hexdigest()
                            unknown_hashes[type_sha] += 1
                            require(len(unknown_hashes) <= 1024, 'unknown_feature_type_category_cap')
                        else:
                            categories['width_unverified'] += 1; type_sha = ''
                        state, reason = typed_row(fields)
                        if reason:
                            failures[reason] += 1
                            self.line(bad, [rows - 1, len(fields), reason, type_sha])
                            break  # Never silently filter an unknown/malformed feature row.
                        identity, name, kind, chrom, start, end = fields
                        all_ids[identity] += 1; annotation_states[state] += 1
                        self.line(axis, [rows - 1, *fields, state])
                        axis_digest.update(('\t'.join(fields) + '\n').encode())
                        if kind == 'Gene Expression':
                            source_genes.append((identity, name, rows - 1))
                            gene_digest.update((identity + '\0' + name + '\n').encode())
                        self.check_time()
            except (gzip.BadGzipFile, EOFError, zlib.error):
                failures['gzip_CRC_or_stream_failure'] += 1
            except QualificationFailure:
                # Decode caps are aggregate hard stops, not permission to scan another file.
                raise
        if not rows:
            failures['empty_source_feature_axis'] += 1
        complete = crc and not failures
        gene_ids = Counter(g for g, _, _ in source_genes)
        stable_ids = Counter(stable_id(g) for g, _, _ in source_genes)
        names = Counter(n for _, n, _ in source_genes)
        return source_genes, dict(person=person, compressed_member=member,
            status='full_annotation_schema_and_CRC_verified' if complete else 'annotation_unqualified',
            full_annotation_schema_verified=complete, gzip_CRC_passed=crc,
            decoded_bytes=decoded, decoded_source_feature_rows=rows,
            full_source_feature_rows=rows if crc else None, field_counts=dict(widths),
            source_feature_types=dict(categories), unknown_type_sha256_counts=dict(unknown_hashes),
            unsupported_reason_counts=dict(failures), unsupported_field_contents_saved=False,
            unsupported_schema_stops_member_without_reading_remainder=True,
            raw_decoded_bytes_sha256=raw_digest.hexdigest(), validated_ordered_annotation_sha256=axis_digest.hexdigest(),
            validated_ordered_Gene_Expression_ID_symbol_sha256=gene_digest.hexdigest(),
            ID_symbol_hash_codec='UTF8_id_NUL_UTF8_symbol_LF_per_native_Gene_Expression_row',
            validated_Gene_Expression_rows=len(source_genes), annotation_states=dict(annotation_states),
            duplicate_validated_feature_ID_extra_rows=sum(n - 1 for n in all_ids.values()),
            duplicate_Gene_Expression_ID_extra_rows=sum(n - 1 for n in gene_ids.values()),
            duplicate_versionless_gene_ID_extra_rows=sum(n - 1 for n in stable_ids.values()),
            repeated_names_extra_rows=sum(n - 1 for n in names.values()),
            names_not_collapsed=True, TSS_bounds_are_not_gene_body_bounds=True,
            RNA_only_feature_file=complete and set(categories) == {'Gene Expression'},
            mixed_declared_gene_and_peak_file=bool(categories['Gene Expression'] and categories['Peaks']),
            source_missing_TSS_is_unavailable_annotation_not_zero_coordinate=True,
            annotation_axis_file=axis_path.name, unsupported_row_metadata_file=bad_path.name)

    def compare(self, person, genes, fixed, annotation):
        exact, bases, fixed_bases = defaultdict(list), defaultdict(list), defaultdict(list)
        for i, (identity, _, _) in enumerate(genes):
            exact[identity].append(i); bases[stable_id(identity)].append(i)
        for i, (identity, _) in enumerate(fixed):
            fixed_bases[stable_id(identity)].append(i)
        counts = Counter(); symbol_states = Counter(); used = set()
        joined = self.output / f'{person}_fixed36601_gene_identity_join.tsv'
        with joined.open('xb', buffering=1024**2) as target:
            self.line(target, ['fixed_gene_index', 'fixed_gene_ID', 'fixed_native_name', 'ID_match_state',
                               'source_Gene_Expression_indices', 'source_feature_indices', 'name_match_state'])
            for i, (identity, name) in enumerate(fixed):
                hits = exact.get(identity, [])
                if hits:
                    state = 'exact' if len(hits) == 1 else 'ambiguous_exact_duplicate'
                else:
                    hits = bases.get(stable_id(identity), [])
                    if not hits:
                        state = 'absent_from_validated_source_Gene_Expression_rows'
                    elif len(hits) == 1 and len(fixed_bases[stable_id(identity)]) == 1:
                        state = 'version_equivalent_ID'
                    else:
                        state = 'ambiguous_versionless_ID'
                names = [genes[j][1] for j in hits]
                name_state = ('absent' if not hits else 'exact' if all(n == name for n in names)
                              else 'different_source_name')
                counts[state] += 1; symbol_states[name_state] += 1; used.update(hits)
                self.line(target, [i, identity, name, state, ','.join(str(j) for j in hits),
                                   ','.join(str(genes[j][2]) for j in hits), name_state])
        extras = self.output / f'{person}_extra_source_gene_identities.tsv'
        with extras.open('xb', buffering=1024**2) as target:
            self.line(target, ['source_Gene_Expression_index', 'source_feature_index', 'source_gene_ID', 'source_name'])
            for i, (identity, name, feature_index) in enumerate(genes):
                if i not in used:
                    self.line(target, [i, feature_index, identity, name])
        complete = annotation['full_annotation_schema_verified']
        all_ids = [g for g, _, _ in genes]; target_ids = [g for g, _ in fixed]
        ordered_exact = complete and all_ids == target_ids
        ordered_pair_exact = complete and [(g, n) for g, n, _ in genes] == fixed
        stable_unique = len(bases) == len(genes) and len(fixed_bases) == len(fixed)
        ordered_version = complete and stable_unique and [stable_id(g) for g in all_ids] == [stable_id(g) for g in target_ids]
        return dict(fixed_rows=36601, source_validated_gene_rows=len(genes), ID_match_counts=dict(counts),
            source_name_match_counts=dict(symbol_states), source_gene_rows_without_any_fixed_ID_match=len(genes) - len(used),
            ordered_exact_ID_axis=ordered_exact, ordered_exact_ID_and_name_axis=ordered_pair_exact,
            ordered_versionless_ID_axis_equivalent=ordered_version,
            versionless_identity_one_to_one=stable_unique, full_gene_axis_comparison_qualified=complete,
            comparisons_are_partial_if_annotation_unqualified=not complete,
            fixed_join_file=joined.name, source_extra_gene_file=extras.name,
            mixed_rows_not_admitted_into_RNA_denominator=annotation['mixed_declared_gene_and_peak_file'],
            count_measurement_and_library_denominator_compatibility_verified=False,
            version_equivalence_is_metadata_classification_not_count_or_reference_equivalence=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-script-sha256', required=True)
    args = parser.parse_args()
    job = os.environ.get('SLURM_JOB_ID', '')
    require(job.isdigit() and re.fullmatch(r'[0-9a-f]{64}', args.expected_script_sha256), 'compute_allocation_and_hash_required')
    output = REC / f'codex_gse223843_rna_annotation_identity_{job}'
    own = False; run = Run(output, job)
    try:
        output.mkdir(exist_ok=False); own = True
        before = run.source_guards(args.expected_script_sha256)
        for name, source in [('executed_source.py', Path(__file__)), ('executed_launcher.sbatch', LAUNCHER)]:
            with (output / name).open('xb') as target:
                run.write(target, source.read_bytes())
        fixed = run.fixed_axis(); roster = run.compressed_roster()
        run.json('source_annotation_contract.json', dict(original_acquisition_job='21999459',
            original_jobs_and_sources_preserved=True, fixed_axis_file=str(GENES), fixed_axis_rows=36601,
            fixed_ordered_ID_name_digest=NATIVE_DIGEST, source_hashes=before, compressed_members=roster,
            source_schema_URL='https://www.10xgenomics.com/support/software/cell-ranger-arc/2.0/analysis/outputs/feature-barcode-matrices',
            installed_schema_source=str(ARC_SOURCE), source_BED_TSS_not_gene_body=True,
            explicit_source_missing_TSS_sentinel=['', '-1', '-1'], **FLAGS))
        for person in PERSONS:
            genes, result = run.decode_member(person, roster[person])
            result['fixed36601_comparison'] = run.compare(person, genes, fixed, result)
            run.results[person] = result
            run.json(person + '_annotation_receipt.json', result)
        run.phase = 'completion_source_hashes'
        require(run.source_guards(args.expected_script_sha256) == before, 'source_inventory_changed')
        require(all(digest(Path(r['path'])) == r['sha256'] for r in roster.values()), 'compressed_annotation_source_changed')
        files = {p.name: dict(bytes=p.stat().st_size, sha256=digest(p)) for p in output.iterdir() if p.is_file()}
        run.phase = 'completed_offline_annotation_attempts'
        summary = dict(annotation_inspection_completed=True, SLURM_job_id=job, all_six_files_attempted=True,
            all_six_full_annotation_qualified=all(r['full_annotation_schema_verified'] for r in run.results.values()),
            results=run.results, source_inventory_before_and_after_identical=True, source_hashes=before,
            output_artifacts=files, decoded_metadata_bytes=run.decoded, decode_cap=CAP,
            per_member_decode_cap=FILE_DECODE_CAP, saved_byte_cap=CAP,
            saved_bytes_before_summary=run.sizes(), elapsed_seconds=time.monotonic() - run.started,
            original_unsupported_parser_results_preserved_not_invalid_RNA_claim=True,
            annotation_identity_does_not_establish_counts_or_expression=True,
            Python_version=sys.version, zlib_runtime=zlib.ZLIB_RUNTIME_VERSION, **FLAGS)
        run.json('summary.json', summary); run.sizes()
        print(json.dumps(dict(annotation_inspection_completed=True, SLURM_job_id=job,
            qualified_files=sum(r['full_annotation_schema_verified'] for r in run.results.values()),
            output=str(output), **FLAGS)), flush=True)
    except Exception as error:
        failure = dict(annotation_inspection_completed=False, phase=run.phase, SLURM_job_id=job,
            reason=str(error) if isinstance(error, QualificationFailure) else 'unexpected_error_contents_suppressed',
            exception_type=type(error).__name__, partial_files_preserved=True,
            decoded_metadata_bytes=run.decoded, completed_member_results=run.results, **FLAGS)
        if own:
            try:
                run.json('failure.json', failure)
            except Exception:
                pass
        print(json.dumps(dict(annotation_inspection_completed=False, SLURM_job_id=job,
                             phase=run.phase, reason=failure['reason'])), flush=True)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
