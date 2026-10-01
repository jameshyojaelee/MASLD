#!/usr/bin/env python3
"""Compute-only elementwise comparison of frozen native/common aggregate counts."""
import argparse
import array
import ast
import csv
import hashlib
import json
import os
from pathlib import Path
import struct
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
HERE = ROOT / 'scripts/analysis/in_biopsy_allelic/rec_codex'
NATIVE = REC / 'codex_gse296875_aligned_native_pseudobulk_21999585'
COMMON = REC / 'codex_gse296875_common_fragment_counts_22000534'
LAUNCHER = HERE / 'run_codex_gse296875_compare_native_common_counts.sbatch'
SHAPE = (273, 306706)
PARTITIONS = ('full', 'halfA', 'halfB')
BLOCK_PEAKS = 1024
EXAMPLES = 20
CAP = 16 * 1024**2
GUARDS = {
    NATIVE / 'summary.json': 'dd1986b643e1e5ae2cec54d255fb8ed7619c6db668541194a4032277d8c00226',
    NATIVE / 'artifacts.json': '3c4102ecde36825237b7d31e6ff4c4395570dc97ee222e274429833cc000d8a5',
    NATIVE / 'units.tsv': '6ceeac2200ad1b452994e4df3fd6a4a60c51f812bcef4daa58116b84d4d69a1b',
    NATIVE / 'unit_identity.tsv': '5c0a1475cb5cbf71ae4e25f31652a822d4b773c4cd8f5f962bcc44b5bb80be68',
    NATIVE / 'ATAC_receipt.json': '4b540a1a0d605b9c1ffc4ec7a4cdc142163df6fc6811c651ecbcf66275f12e97',
    COMMON / 'summary.json': 'bd64347bc066cbd182695c466b957ddd3c6efc689b53f8db2d557aed34321f8e',
    COMMON / 'artifacts.json': '8735af6f215fc82be928f2349cb8182d3d76e834304cce8eb32e858f5880a960',
    COMMON / 'manifest.json': '36d516811ca1d80f67ad6fdabc11af9d35734d7afdc6b09661a41769ee51ea4e',
    COMMON / 'units.tsv': 'babe4cd3bc0f4420c6303526d9ebbe342ae69ca792237a2c05b2bf09aa780f2a',
    COMMON / 'unit_identity.tsv': '509dadfb442445c28eb4139d2f3a9699d6734c82a3ac6442b29e2d767417a964',
}
for _directory in (NATIVE, COMMON):
    GUARDS[_directory / 'cell_identity.tsv'] = '438f972c0883837fc254260a82f058caf95d634643d2a2b13156ba4089b848d4'
    GUARDS[_directory / 'native_peak_geometry.tsv'] = '50a52d097a5232edcf136dd53d6549ca2d1aa5dce9a57cf055ff10b7b2e0f3e3'


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(path):
    require(path.is_file() and not path.is_symlink(), 'Missing regular source or file symlink: ' + str(path))
    value = path.stat()
    return dict(bytes=value.st_size, inode=value.st_ino, device=value.st_dev,
                mtime_ns=value.st_mtime_ns, ctime_ns=value.st_ctime_ns)


def saved_bytes(output):
    return sum(p.stat().st_size for p in output.rglob('*') if p.is_file())


def write_json(output, name, value):
    path = output / name
    temporary = output / (name + '.partial')
    require(not path.exists() and not temporary.exists(), 'Refusing output overwrite')
    with temporary.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    require(saved_bytes(output) <= CAP, 'Inclusive comparison output exceeds 16 MiB')
    temporary.rename(path)


def table(path):
    identity = ('unit_index', 'donor_id', 'well_id', 'native_source_celltype', 'nuclei', 'halfA_nuclei', 'halfB_nuclei')
    totals = tuple('ATAC_' + p + '_total_native_peak_events' for p in PARTITIONS) if path.parent == NATIVE else tuple(p + '_summed_peak_counts' for p in PARTITIONS)
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(set(identity + totals) <= set(reader.fieldnames or ()), 'Unit metadata schema differs')
        rows = [{k: row[k] for k in identity + totals} for row in reader]
    require(len(rows) == SHAPE[0] and [int(row['unit_index']) for row in rows] == list(range(SHAPE[0])), 'Ordered unit dimension differs')
    return rows, identity


def header(path, fortran, artifact):
    with path.open('rb') as handle:
        require(handle.read(8) == b'\x93NUMPY\x01\x00', 'NPY magic/version differs')
        size = struct.unpack('<H', handle.read(2))[0]
        require(0 < size <= 65536, 'NPY header length differs')
        value = ast.literal_eval(handle.read(size).decode('latin1').strip())
        offset = handle.tell()
        handle.seek(0)
        digest = hashlib.sha256(handle.read(offset)).hexdigest()
    require(set(value) == {'descr', 'fortran_order', 'shape'} and value['descr'] == '<i4' and
            value['shape'] == SHAPE and value['fortran_order'] is fortran and offset == 128,
            'Count shape/dtype/order differs')
    require(path.stat().st_size == offset + SHAPE[0] * SHAPE[1] * 4 == artifact['bytes'], 'NPY payload byte dimension differs')
    return dict(dtype=value['descr'], shape=list(value['shape']), fortran_order=fortran,
                offset=offset, header_sha256=digest)


def integers(data, count):
    require(len(data) == count * 4, 'Truncated count block')
    values = array.array('i')
    values.frombytes(data)
    if sys.byteorder != 'little':
        values.byteswap()
    require(not values or min(values) >= 0, 'Negative source counts')
    return values


def compare_partition(partition, native_path, common_path, native_header, common_header, native_units, common_units):
    totals_native = [0] * SHAPE[0]
    totals_common = [0] * SHAPE[0]
    mismatches, maximum, examples, mismatch_units = 0, 0, [], set()
    began = time.monotonic()
    with native_path.open('rb') as native, common_path.open('rb') as common:
        for first in range(0, SHAPE[1], BLOCK_PEAKS):
            columns = min(BLOCK_PEAKS, SHAPE[1] - first)
            # Native Fortran storage has all 273 units contiguous for each peak.
            native.seek(native_header['offset'] + first * SHAPE[0] * 4)
            source = integers(native.read(columns * SHAPE[0] * 4), columns * SHAPE[0])
            for unit in range(SHAPE[0]):
                left = source[unit::SHAPE[0]]
                # Common C storage has all 306706 peaks contiguous for each unit.
                common.seek(common_header['offset'] + (unit * SHAPE[1] + first) * 4)
                right = integers(common.read(columns * 4), columns)
                require(len(left) == len(right) == columns, 'Fortran/C slice dimension differs')
                totals_native[unit] += sum(left)
                totals_common[unit] += sum(right)
                if left != right:
                    mismatch_units.add(unit)
                    for local, (a, b) in enumerate(zip(left, right)):
                        if a != b:
                            mismatches += 1
                            maximum = max(maximum, abs(a - b))
                            if len(examples) < EXAMPLES:
                                examples.append(dict(unit_index=unit, peak_index=first + local))
            if first // BLOCK_PEAKS % 64 == 0:
                print(json.dumps(dict(partition=partition, completed_peak_columns=first + columns,
                                      total_peak_columns=SHAPE[1], mismatched_elements=mismatches)), flush=True)
    require(totals_native == [int(u['ATAC_' + partition + '_total_native_peak_events']) for u in native_units],
            'Read native counts do not reconcile with immutable unit totals')
    require(totals_common == [int(u[partition + '_summed_peak_counts']) for u in common_units],
            'Read common counts do not reconcile with immutable unit totals')
    return dict(partition=partition, checked_elements=SHAPE[0] * SHAPE[1], mismatched_elements=mismatches,
                arrays_equal=mismatches == 0, maximum_absolute_difference=maximum,
                mismatched_units=len(mismatch_units), first_bounded_indices=examples,
                example_order='peak block, then unit, then peak within block; first 20 mismatches',
                native_grand_total=sum(totals_native), common_grand_total=sum(totals_common),
                all_273_persisted_row_totals_reconciled=True, elapsed_seconds=time.monotonic() - began)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true')
    parser.add_argument('--expected-helper-sha')
    args = parser.parse_args()
    if args.inspect:
        print(json.dumps(dict(role='data-only equality comparison', dimensions=list(SHAPE), partitions=PARTITIONS,
                              checked_elements=3 * SHAPE[0] * SHAPE[1], block_peaks=BLOCK_PEAKS,
                              count_block_bytes=SHAPE[0] * BLOCK_PEAKS * 4,
                              cpu=1, memory_GiB=2, hours=1, output_cap_bytes=CAP), indent=2))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdecimal() and os.environ.get('SLURM_CPUS_PER_TASK') == '1',
            'Comparison requires an individual one-CPU SLURM compute job')
    require(args.expected_helper_sha and sha(Path(__file__)) == args.expected_helper_sha, 'Reviewed comparison source drift')
    require(array.array('i').itemsize == 4, 'Installed stdlib int32 reader itemsize differs')
    output = REC / ('codex_gse296875_native_common_array_comparison_' + os.environ['SLURM_JOB_ID'])
    output.mkdir(mode=0o750)
    began = time.monotonic()
    try:
        archive = output / 'executed_sources'
        archive.mkdir()
        source_pins = {}
        for source in (Path(__file__), LAUNCHER):
            contents = source.read_bytes()
            (archive / source.name).write_bytes(contents)
            source_pins[str(source)] = hashlib.sha256(contents).hexdigest()
        require(source_pins[str(Path(__file__))] == args.expected_helper_sha, 'Source changed before archive')
        snapshots = {}
        for path, expected in GUARDS.items():
            before = fingerprint(path)
            require(sha(path) == expected and fingerprint(path) == before, 'Immutable source metadata differs: ' + str(path))
            snapshots[str(path)] = before
        native_units, identity = table(NATIVE / 'units.tsv')
        common_units, _ = table(COMMON / 'units.tsv')
        require([[u[k] for k in identity] for u in native_units] == [[u[k] for k in identity] for u in common_units],
                'Exact ordered native/common units differ')
        require(len({u['donor_id'] for u in native_units}) == 39 and sum(int(u['nuclei']) for u in native_units) == 68398,
                'Source person/nucleus axis differs')
        require(json.loads((NATIVE / 'summary.json').read_text())['aligned_native_count_aggregation_completed'] and
                json.loads((COMMON / 'summary.json').read_text())['common_fragment_overlap_counts_completed'],
                'Source producer completion markers differ')
        native_artifacts = json.loads((NATIVE / 'artifacts.json').read_text())
        common_artifacts = json.loads((COMMON / 'artifacts.json').read_text())
        counts = []
        for partition in PARTITIONS:
            native_name = 'ATAC_' + partition + '_counts_int32.npy'
            common_name = 'ATAC_fragment_overlap_' + partition + '_int32.npy'
            for path, artifact, order in ((NATIVE / native_name, native_artifacts['artifacts'][native_name], True),
                                          (COMMON / common_name, common_artifacts[common_name], False)):
                before = fingerprint(path)
                expected = artifact['sha256'] if order else artifact['computed_sha256']
                require(sha(path) == expected and fingerprint(path) == before, 'Count payload SHA/stat differs')
                parsed = header(path, order, artifact)
                recorded = native_artifacts['array_headers'][native_name] if order else artifact['header']
                require(parsed['header_sha256'] == recorded['header_sha256'], 'Recorded count header identity differs')
                counts.append(dict(path=str(path), partition=partition, producer='native' if order else 'common',
                                   computed_sha256=expected, source_stat_before=before, header=parsed))
        readers = {str(Path(sys.executable).resolve()): sha(Path(sys.executable).resolve()),
                   str(Path(array.__file__).resolve()): sha(Path(array.__file__).resolve())}
        write_json(output, 'manifest.json', dict(source_pins=source_pins, metadata_pins={str(p): v for p, v in GUARDS.items()},
                   source_stat_before=snapshots, input_counts=counts, argv=sys.argv, python_version=sys.version,
                   installed_stdlib_reader_hashes=readers, byteorder=sys.byteorder, resources=dict(cpus=1,memory_GiB=2,hours=1),
                   block_peaks=BLOCK_PEAKS, largest_count_block_bytes=SHAPE[0] * BLOCK_PEAKS * 4,
                   interpretation='Elementwise equality of aggregate unit-by-peak counts in this observed source only.',
                   no_normalization_fitting_scoring_or_filtering=True, no_RDS_fragment_RNA_or_external_payload_read=True,
                   no_count_arrays_written=True, output_cap_bytes_inclusive=CAP))
        results = []
        for partition in PARTITIONS:
            pair = [v for v in counts if v['partition'] == partition]
            result = compare_partition(partition, Path(pair[0]['path']), Path(pair[1]['path']),
                                       pair[0]['header'], pair[1]['header'], native_units, common_units)
            results.append(result)
            write_json(output, partition + '_comparison.json', result)
        # Rehash count inputs and immutable metadata after comparison; no source is edited.
        for source in counts:
            path = Path(source['path'])
            require(fingerprint(path) == source['source_stat_before'] and sha(path) == source['computed_sha256'],
                    'Count input changed during comparison')
            require(fingerprint(path) == source['source_stat_before'], 'Count input changed during final hash')
        for path, expected in GUARDS.items():
            require(fingerprint(path) == snapshots[str(path)] and sha(path) == expected, 'Metadata changed during comparison')
            require(fingerprint(path) == snapshots[str(path)], 'Metadata changed during final hash')
        for source, expected in source_pins.items():
            require(sha(Path(source)) == expected, 'Comparison helper/launcher changed')
        for source, expected in readers.items():
            require(sha(Path(source)) == expected, 'Installed stdlib reader changed')
        summary = dict(comparison_completed=True, checked_partitions=3,
                       checked_aggregate_elements=sum(r['checked_elements'] for r in results),
                       all_arrays_equal=all(r['arrays_equal'] for r in results),
                       mismatched_elements=sum(r['mismatched_elements'] for r in results),
                       maximum_absolute_difference=max(r['maximum_absolute_difference'] for r in results),
                       exact_273_unit_306706_peak_axis_verified=True, all_source_hashes_verified_before_and_after=True,
                       clinical_model_normalization_or_external_work_run=False,
                       qualification='Equality, if observed, concerns these stored aggregates; it does not prove the author original algorithm or compatibility with another source.',
                       elapsed_seconds=time.monotonic() - began, output_cap_bytes_inclusive=CAP,
                       saved_bytes_before_final_summary=saved_bytes(output), final_inclusive_cap_checked=True)
        write_json(output, 'summary.json', summary)
        print(json.dumps(summary), flush=True)
    except BaseException:
        write_json(output, 'failure.json', dict(completed=False, traceback=traceback.format_exc(), partial_outputs_preserved=True))
        raise


if __name__ == '__main__':
    main()
