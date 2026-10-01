#!/usr/bin/env python3
"""Compute-only, separately declared native fragment-body-overlap measurement."""
import argparse
import array
import ast
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import struct
import subprocess
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
BENCH = ROOT / 'Analysis/MASLD_Model_Benchmark'
REC = BENCH / 'executions/codex-rec-20260929T142434Z'
HERE = ROOT / 'scripts/analysis/in_biopsy_allelic/rec_codex'
ALIGNED = REC / 'codex_gse296875_aligned_native_pseudobulk_21999585'
SPLIT = BENCH / 'executions/gse296875-donor-lineage-fragments-21063889'
RAW = BENCH / 'executions/gse296875-fragments-direct-21063797'
MEMBERSHIP = BENCH / 'executions/gse296875-fragment-membership-21063829/cell_membership.tsv.gz'
LOCK = ROOT / 'Analysis/SingleCell/continual_integration/reference/gse296875_source_lock_v44.json'
ENV = Path('/gpfs/commons/home/jameslee/micromamba/envs/rnaseq')
CPP = HERE / 'codex_gse296875_common_fragment_counter.cpp'
LAUNCHER = HERE / 'run_codex_gse296875_common_fragment_counts.sbatch'
PARTITIONS = ('full', 'halfA', 'halfB')
LABELS = ('B cells', 'Cholangiocytes', 'Hepatocytes', 'Kupffer', 'LSEC', 'Mesenchymal', 'NK-T')
ALIASES = dict(zip(LABELS, ('b_cell', 'cholangiocyte', 'hepatocyte', 'macrophage',
                          'endothelial_cell', 'fibroblast', 't_cell')))
HALF_PREFIX = 'MASLD-NATIVE-HALVES-20261001\t'
CAP = 2 * 1024**3
SHAPE = (273, 306706)
RECORDS = 1091589686
READ_SUPPORT = 1441243289
INPUT_BYTES = 17012545862
GUARDS = {
    MEMBERSHIP: '81c43ee9af4a05784166c84a6b9ea0513f7b697ebbe33d9ff89982c52aa4e4c6',
    LOCK: '6be686b4187309f038097cf73d3f2a8dd34e72d666bc92733a53e99a1ccbb620',
    RAW / 'ARTIFACTS.json': '038bb09e1056d86bb1119137450fdef856f9129cb08cb8b2fdb27735d869a303',
    SPLIT / 'ARTIFACTS.json': 'ca8f157669d5a6953765525b62ab593db1d4a352e99a73ed1b6db353aab35ef4',
    SPLIT / 'contract.json': '7747291cf5cda2bb5a89e42105c729d6e241bce4e2b57ce0a328dfff575b4b43',
    SPLIT / 'fragment_manifest.tsv': '8dca7ca4d80302ef263e986e52fe3eb49cc3d565d0c7f411df9a2effedc5b5b6',
    SPLIT / 'validation/well_split_stats.json': 'ebe252d0b8fd5ecb89fbb88da94f9e6576a95499f72f62047a44818145d4007b',
    ALIGNED / 'summary.json': 'dd1986b643e1e5ae2cec54d255fb8ed7619c6db668541194a4032277d8c00226',
    ALIGNED / 'artifacts.json': '3c4102ecde36825237b7d31e6ff4c4395570dc97ee222e274429833cc000d8a5',
    ALIGNED / 'cell_identity.tsv': '438f972c0883837fc254260a82f058caf95d634643d2a2b13156ba4089b848d4',
    ALIGNED / 'unit_identity.tsv': '5c0a1475cb5cbf71ae4e25f31652a822d4b773c4cd8f5f962bcc44b5bb80be68',
    ALIGNED / 'native_peak_geometry.tsv': '50a52d097a5232edcf136dd53d6549ca2d1aa5dce9a57cf055ff10b7b2e0f3e3',
}
RAW_HEADER_SHA = (
    'c953b967cc476927932f08e745663ff417be9f925e119fb266e05f22364c9ba9',
    '2a3aa2f71d5c16e831ddf118b254835186d2e6dc58ab0a1c629289f70961d1da',
    '03b3c79596a8950ac563a4cda21602f61241dbfceed44583e7a56334aaf9ee31',
    'f8dcc145d1ff7ba9732d17712492622fb33e049fb26d4cc1806f6ac59c3bf35e',
    'beef39cc5b417b8c7c69da6f3750d70e024235aa67950f36ebfc5d9a1a699846',
    '3fa63806e8af37e9141b969fe45c84f1c94851f213b177314a3f5fefb17a6da0',
    '1a1fe89bf42ccf9b03ef3ea44bb11fbcab753af487f0538ff989331322719a40',
    '08b8cdfbc08556c237723f456bf2808cbbdadbd1e95d35306d1a3e475485e701',
)
DENOMINATORS = (
    'primary_fragment_records', 'autosomal_fragment_records', 'union_overlapped_fragments',
    'readSupport_diagnostic', 'summed_peak_counts',
)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def fingerprint(path):
    path = Path(path)
    require(path.is_file() and not path.is_symlink(), 'Missing regular source or file symlink: ' + str(path))
    info = path.stat()
    return dict(bytes=info.st_size, mtime_ns=info.st_mtime_ns, ctime_ns=info.st_ctime_ns,
                inode=info.st_ino, device=info.st_dev)


def saved_bytes(output):
    return sum(path.stat().st_size for path in output.rglob('*') if path.is_file())


def check_cap(output):
    require(saved_bytes(output) <= CAP, 'Inclusive saved output exceeds 2 GiB')


def write_json(path, value):
    partial = Path(str(path) + '.partial')
    require(not path.exists() and not partial.exists(), 'Refusing receipt overwrite')
    with partial.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    partial.rename(path)


def table(path, wanted, compressed=False):
    opener = gzip.open if compressed else open
    with opener(path, 'rt', encoding='utf-8', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(set(wanted) <= set(reader.fieldnames or ()), 'Identity schema differs: ' + str(path))
        return [{key: row[key] for key in wanted} for row in reader]


def write_tsv(path, fields, rows):
    with path.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def source_guards():
    snapshots = {}
    for path, expected in GUARDS.items():
        before = fingerprint(path)
        require(sha(path) == expected, 'Frozen source SHA256 differs: ' + str(path))
        require(fingerprint(path) == before, 'Source changed while hashing: ' + str(path))
        snapshots[str(path)] = before
    return snapshots


def identities():
    names = ('well_id', 'raw_barcode', 'cell_id', 'donor_id', 'source_label')
    original = table(MEMBERSHIP, names, compressed=True)
    frozen = table(ALIGNED / 'cell_identity.tsv', names + ('half_sort_sha256', 'nucleus_half'))
    require(len(original) == len(frozen) == 68398, 'Native cell count differs')
    require([tuple(row[n] for n in names) for row in original] ==
            [tuple(row[n] for n in names) for row in frozen], 'Aligned/native membership ordered join differs')
    require(len({row['cell_id'] for row in frozen}) == 68398, 'Duplicate native nucleus')
    grouped, donor_wells = {}, {}
    cell_hash = hashlib.sha256()
    for row in frozen:
        require(row['source_label'] in LABELS and row['donor_id'].isdecimal(), 'Native donor/label differs')
        require(re.fullmatch(r'well[1-8]', row['well_id']) is not None and
                re.fullmatch(r'[ACGT]+-[0-9]+', row['raw_barcode']) is not None and
                row['cell_id'] == row['well_id'] + '_' + row['raw_barcode'], 'Well namespace differs')
        require(donor_wells.setdefault(row['donor_id'], row['well_id']) == row['well_id'], 'Donor crosses wells')
        digest = hashlib.sha256((HALF_PREFIX + row['cell_id']).encode('utf-8')).hexdigest()
        require(digest == row['half_sort_sha256'], 'Frozen nucleus half hash differs')
        cell_hash.update(row['cell_id'].encode('utf-8') + b'\n')
        grouped.setdefault((row['donor_id'], row['source_label']), []).append(row)
    require(cell_hash.hexdigest() == 'a3d832f1128d2acea32a68d4c4287192a292695529722347aab1e4a07ce71217',
            'Native cell order SHA256 differs')
    for group in grouped.values():
        ordered = sorted(group, key=lambda row: (row['half_sort_sha256'], row['cell_id']))
        for position, row in enumerate(ordered):
            require(row['nucleus_half'] == ('A' if position < len(ordered) // 2 else 'B'), 'Half assignment differs')
    donors = sorted(donor_wells, key=int)
    require(len(donors) == 39 and set(donor_wells.values()) == {'well' + str(i) for i in range(1, 9)},
            'Donor/well axes differ')
    units = []
    for donor in donors:
        for label in LABELS:
            n = len(grouped.get((donor, label), ()))
            units.append(dict(unit_index=len(units), donor_id=donor, well_id=donor_wells[donor],
                              native_source_celltype=label, nuclei=n, halfA_nuclei=n // 2,
                              halfB_nuclei=n - n // 2, sampling_state='observed' if n else 'unsampled',
                              halfA_missing=n // 2 == 0, halfB_missing=n - n // 2 == 0))
    require([(u['donor_id'], u['native_source_celltype']) for u in units if not u['nuclei']] == [('733', 'B cells')],
            'Unsampled native stratum differs')
    fields = ('unit_index', 'donor_id', 'well_id', 'native_source_celltype', 'nuclei', 'halfA_nuclei', 'halfB_nuclei')
    aligned_units = table(ALIGNED / 'unit_identity.tsv', fields)
    require([[str(u[key]) for key in fields] for u in units] ==
            [[u[key] for key in fields] for u in aligned_units], 'Frozen ordered unit axis differs')
    return frozen, units, grouped


def fragment_roster(units):
    fields = ('donor_id', 'lineage_id', 'nuclei', 'wells', 'records', 'read_support', 'path', 'size_bytes', 'sha256')
    rows = table(SPLIT / 'fragment_manifest.tsv', fields)
    by_key = {(row['donor_id'], row['lineage_id']): row for row in rows}
    require(len(rows) == len(by_key) == 272, 'Selected fragment roster differs')
    manifest = json.loads((SPLIT / 'ARTIFACTS.json').read_text())
    artifacts = {row['path']: row for row in manifest['artifacts']}
    require(len(artifacts) == len(manifest['artifacts']), 'Duplicate selected artifact path')
    contract = json.loads((SPLIT / 'contract.json').read_text())
    require(contract['assay_signal_unit'] == 'one_unique_fragment_record_not_readSupport' and
            contract['source_coordinates'] == 'cellranger_arc_2_0_0_tn5_adjusted' and
            contract['barcode_policy'] == 'well_namespaced_cell_id', 'Fragment source unit contract differs')
    expected = set()
    for unit in units:
        key = (unit['donor_id'], ALIASES[unit['native_source_celltype']])
        if not unit['nuclei']:
            require(key not in by_key, 'Unexpected fragment source for unsampled unit')
            continue
        expected.add(key)
        row = by_key[key]
        require(int(row['nuclei']) == unit['nuclei'] and row['wells'] == unit['well_id'], 'Fragment cell/unit join differs')
        relative = f"fragments/donor_{unit['donor_id']}/{key[1]}.fragments.tsv.gz"
        require(row['path'] == relative and relative in artifacts, 'Fragment native alias/path differs')
        artifact = artifacts[relative]
        require(int(row['size_bytes']) == artifact['size_bytes'] and row['sha256'] == artifact['sha256'],
                'Selected fragment receipt identities differ')
        row['unit_index'] = unit['unit_index']
        row['source_path'] = str(SPLIT / relative)
        require(fingerprint(row['source_path'])['bytes'] == int(row['size_bytes']), 'Selected fragment size differs')
    require(set(by_key) == expected and sum(int(r['records']) for r in rows) == RECORDS and
            sum(int(r['read_support']) for r in rows) == READ_SUPPORT and
            sum(int(r['size_bytes']) for r in rows) == INPUT_BYTES, 'Full selected source totals differ')
    return sorted(rows, key=lambda row: row['unit_index'])


def raw_header_guards():
    raw_manifest = json.loads((RAW / 'ARTIFACTS.json').read_text())
    raw_artifacts = {row['path']: row for row in raw_manifest['artifacts']}
    split_stats = json.loads((SPLIT / 'validation/well_split_stats.json').read_text())
    require(len(split_stats) == 8 and {r['well_id'] for r in split_stats} == {'well' + str(i) for i in range(1, 9)},
            'Split well-header evidence differs')
    evidence = []
    for i in range(1, 9):
        relative = f'fragments/GSM{8979149 + 2 * i}_well{i}_atac_fragments.tsv.gz'
        path = RAW / relative
        before = fingerprint(path)
        require(before['bytes'] == raw_artifacts[relative]['size_bytes'], 'Original fragment size differs')
        header = bytearray()
        values = {}
        with gzip.open(path, 'rb') as handle:
            while True:
                marker = handle.read(1)
                if marker != b'#':
                    break  # No complete fragment record is retrieved.
                line = marker + handle.readline(8192)
                require(line.endswith(b'\n') and len(header) + len(line) <= 65536, 'Original header exceeds bound')
                header.extend(line)
                payload = line[1:].decode('utf-8').strip()
                if '=' in payload:
                    key, value = payload.split('=', 1)
                    values[key] = value
        require(hashlib.sha256(header).hexdigest() == RAW_HEADER_SHA[i - 1], 'Original fragment header SHA256 differs')
        require(values['pipeline_name'] == 'cellranger-arc' and values['pipeline_version'] == 'cellranger-arc-2.0.0'
                and values['reference_version'] == '2020-A', 'Original header pipeline/reference differs')
        require(fingerprint(path) == before, 'Original fragment changed while reading header')
        stat = next(row for row in split_stats if row['well_id'] == 'well' + str(i))
        require(stat['source_pipeline'] == values['pipeline_name'] and
                stat['source_pipeline_version'] == values['pipeline_version'] and
                stat['source_reference_version'] == values['reference_version'], 'Historical split header evidence differs')
        evidence.append(dict(well_id='well' + str(i), path=str(path), stat=before,
                             computed_header_sha256=RAW_HEADER_SHA[i - 1], header_fields=values,
                             recorded_full_compressed_sha256=raw_artifacts[relative]['sha256'],
                             full_raw_compressed_sha256_recomputed=False))
    return evidence


def npy_header(path):
    with path.open('rb') as handle:
        require(handle.read(6) == b'\x93NUMPY' and handle.read(2) == b'\x01\x00', 'NPY magic/version differs')
        size = struct.unpack('<H', handle.read(2))[0]
        require(0 < size <= 65536, 'NPY header length differs')
        header = ast.literal_eval(handle.read(size).decode('latin1').strip())
        require(set(header) == {'descr', 'fortran_order', 'shape'}, 'NPY header schema differs')
        offset = handle.tell()
        handle.seek(0)
        digest = hashlib.sha256(handle.read(offset)).hexdigest()
    return dict(shape=list(header['shape']), dtype=header['descr'], fortran_order=header['fortran_order'],
                offset=offset, header_sha256=digest)


def inherited_headers():
    artifact = json.loads((ALIGNED / 'artifacts.json').read_text())
    output = {}
    for modality, dtype, columns, order in (('RNA', '<u8', 36601, False), ('ATAC', '<i4', 306706, True)):
        for partition in PARTITIONS:
            name = f'{modality}_{partition}_counts_{"uint64" if modality == "RNA" else "int32"}.npy'
            path = ALIGNED / name
            header = npy_header(path)
            recorded = artifact['array_headers'][name]
            require(header['shape'] == [273, columns] and header['dtype'] == dtype and
                    header['fortran_order'] is order and header['header_sha256'] == recorded['header_sha256'] and
                    path.stat().st_size == artifact['artifacts'][name]['bytes'], 'Aligned source array header/size differs')
            output[name] = dict(path=str(path), header=header, stat=fingerprint(path),
                                recorded_payload_sha256=artifact['artifacts'][name]['sha256'], payload_read=False)
    return output


def make_npy(path):
    body = repr(dict(descr='<i4', fortran_order=False, shape=SHAPE)).encode('latin1')
    padding = (64 - (10 + len(body) + 1) % 64) % 64
    body += b' ' * padding + b'\n'
    with path.open('xb') as handle:
        handle.write(b'\x93NUMPY\x01\x00' + struct.pack('<H', len(body)) + body)
        offset = handle.tell()
        handle.truncate(offset + SHAPE[0] * SHAPE[1] * 4)
        handle.flush()
        os.fsync(handle.fileno())
    return offset


def run_unit(row, executable, output, offset, arrays, membership_path):
    path = Path(row['source_path'])
    before = fingerprint(path)
    require(before['bytes'] == int(row['size_bytes']) and sha(path) == row['sha256'], 'Selected compressed fragment hash differs')
    require(fingerprint(path) == before, 'Selected source changed during SHA256')
    counter_receipt = output / 'unit_receipts' / f"unit_{row['unit_index']:03d}_counter.json"
    log = output / 'unit_receipts' / f"unit_{row['unit_index']:03d}.log"
    command = [str(executable), str(path), str(output / 'native_peak_geometry.tsv'), str(membership_path)]
    command += [str(array_path) for array_path in arrays]
    command += [str(offset + row['unit_index'] * SHAPE[1] * 4), row['records'], str(counter_receipt)]
    with log.open('x') as handle:
        result = subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=False)
    require(result.returncode == 0, f"Unit {row['unit_index']} counter failed; see {log}")
    require(fingerprint(path) == before, 'Selected source changed while counting')
    receipt = json.loads(counter_receipt.read_text())
    require(receipt['records'] == int(row['records']) and receipt['readSupport_diagnostic'][0] == int(row['read_support'])
            and receipt['seen_cells'] == int(row['nuclei']) and receipt['duplicate_keys'] == 0
            and receipt['gzip_crc_and_eof_verified'] and receipt['full_equals_halves_per_peak'], 'Streamed unit/source invariants differ')
    for name in DENOMINATORS:
        values = receipt[name]
        require(len(values) == 3 and all(isinstance(v, int) and v >= 0 for v in values)
                and values[0] == values[1] + values[2], 'Unit denominator full=A+B failed')
    receipt.update(unit_index=row['unit_index'], source_path=str(path),
                   computed_compressed_sha256=row['sha256'], source_stat_before=before, source_stat_after=fingerprint(path),
                   count_argv=command)
    write_json(output / 'unit_receipts' / f"unit_{row['unit_index']:03d}_verified.json", receipt)
    return receipt


def verify_persisted(arrays, offset, receipts, units):
    require(array.array('i').itemsize == 4, 'Python int32 itemsize differs')
    digests = [hashlib.sha256() for _ in arrays]
    handles = [path.open('rb') for path in arrays]
    try:
        for i, handle in enumerate(handles):
            digests[i].update(handle.read(offset))
        for unit in units:
            blocks = [handle.read(SHAPE[1] * 4) for handle in handles]
            require(all(len(b) == SHAPE[1] * 4 for b in blocks), 'Persisted NPY row truncated')
            values = []
            for i, block in enumerate(blocks):
                digests[i].update(block)
                row = array.array('i')
                row.frombytes(block)
                if sys.byteorder != 'little':
                    row.byteswap()
                values.append(row)
            require(all(f >= 0 and a >= 0 and b >= 0 and f == a + b and f <= RECORDS
                        for f, a, b in zip(*values)), 'Persisted count full=A+B or int32 bound failed')
            sums = [sum(row) for row in values]
            expected = receipts[unit['unit_index']]['summed_peak_counts'] if unit['nuclei'] else [0, 0, 0]
            require(sums == expected, 'Persisted per-unit count totals differ')
        require(all(not handle.read(1) for handle in handles), 'NPY has trailing payload')
    finally:
        for handle in handles:
            handle.close()
    return [digest.hexdigest() for digest in digests]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true', help='Static source/resource inspection only')
    parser.add_argument('--expected-python-sha')
    parser.add_argument('--expected-cpp-sha')
    args = parser.parse_args()
    if args.inspect:
        print(json.dumps(dict(role='distinct native fragment-body-overlap counts, data only', shape=SHAPE,
                              workers=8, cpu=8, memory_GiB=32, hours=90, output_cap_bytes=CAP,
                              fragment_bytes=INPUT_BYTES, helper=str(Path(__file__)), counter=str(CPP)), indent=2))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdecimal(), 'Count execution requires a SLURM compute job')
    require(int(os.environ.get('SLURM_CPUS_PER_TASK', '0')) == 8, 'Requires exactly 8 allocated CPUs')
    require(args.expected_python_sha and args.expected_cpp_sha, 'Expected reviewed helper/counter hashes are required')
    require(sha(Path(__file__)) == args.expected_python_sha and sha(CPP) == args.expected_cpp_sha, 'Reviewed source SHA256 drift')
    output = REC / ('codex_gse296875_common_fragment_counts_' + os.environ['SLURM_JOB_ID'])
    output.mkdir(mode=0o750)
    began = time.monotonic()
    try:
        archive = output / 'executed_sources'
        archive.mkdir()
        source_hashes = {}
        for source in (Path(__file__), CPP, LAUNCHER):
            value = source.read_bytes()
            (archive / source.name).write_bytes(value)
            source_hashes[str(source)] = hashlib.sha256(value).hexdigest()
        require(source_hashes[str(Path(__file__))] == args.expected_python_sha and
                source_hashes[str(CPP)] == args.expected_cpp_sha, 'Source changed before archive')
        compiler = ENV / 'bin/g++'
        library = (ENV / 'lib/libz.so').resolve(strict=True)
        tools = {str(path): dict(sha256=sha(path), stat=fingerprint(path.resolve()))
                 for path in (compiler, library, ENV / 'include/zlib.h', Path(sys.executable))}
        executable = archive / 'common_fragment_counter'
        command = [str(compiler), '-O3', '-std=c++17', '-Wall', '-Wextra', '-Werror',
                   '-I' + str(ENV / 'include'), '-L' + str(ENV / 'lib'), '-Wl,-rpath,' + str(ENV / 'lib'),
                   str(archive / CPP.name), '-lz', '-o', str(executable)]
        with (output / 'compile.log').open('x') as handle:
            subprocess.run(command, stdout=handle, stderr=subprocess.STDOUT, check=True)
        self_check = output / 'self_checks'
        self_check.mkdir()
        subprocess.run([str(executable), '--self-check', str(self_check)], check=True)
        frozen_stats = source_guards()
        cells, units, groups = identities()
        rows = fragment_roster(units)
        raw_headers = raw_header_guards()
        aligned_headers = inherited_headers()
        require(json.loads((ALIGNED / 'summary.json').read_text())['aligned_native_count_aggregation_completed'],
                'Aligned source aggregation was not completed')
        shutil.copyfile(ALIGNED / 'native_peak_geometry.tsv', output / 'native_peak_geometry.tsv')
        require(sha(output / 'native_peak_geometry.tsv') == GUARDS[ALIGNED / 'native_peak_geometry.tsv'], 'Copied geometry differs')
        write_tsv(output / 'cell_identity.tsv', tuple(cells[0]), cells)
        write_tsv(output / 'unit_identity.tsv', tuple(units[0]), units)
        membership_dir = output / 'unit_membership'
        membership_dir.mkdir()
        (output / 'unit_receipts').mkdir()
        for unit in units:
            if unit['nuclei']:
                selected = groups[unit['donor_id'], unit['native_source_celltype']]
                write_tsv(membership_dir / f"unit_{unit['unit_index']:03d}.tsv", ('cell_id', 'nucleus_half'),
                          [{key: row[key] for key in ('cell_id', 'nucleus_half')} for row in selected])
        arrays = [output / f'ATAC_fragment_overlap_{partition}_int32.npy.partial' for partition in PARTITIONS]
        offsets = [make_npy(path) for path in arrays]
        require(len(set(offsets)) == 1, 'Output NPY header offsets differ')
        for path in arrays:
            header = npy_header(path)
            require(header['shape'] == list(SHAPE) and header['dtype'] == '<i4' and not header['fortran_order'],
                    'New array header differs')
        check_cap(output)
        manifest = dict(scientific_role='exposed public development, data-only distinct fragment-body-overlap endpoint',
                        argv=sys.argv, resources=dict(cpus=8, memory_GiB=32, time_hours=90),
                        source_hashes=source_hashes, tools=tools, compiler_argv=command,
                        compiler_version=subprocess.check_output([str(compiler), '--version'], text=True),
                        zlib_runtime_checked_by_counter=True, compiled_counter_sha256=sha(executable),
                        python_version=sys.version, source_guards={str(p): v for p, v in GUARDS.items()},
                        source_stat_before=frozen_stats, selected_fragment_inputs=rows,
                        original_comment_header_evidence=raw_headers, aligned_headers_only=aligned_headers,
                        exact_native_shape=list(SHAPE), nucleus_half_hash_prefix=HALF_PREFIX,
                        interval_rule='fragment.start < peak.BED_end and fragment.end > peak.BED_start',
                        unique_key=['well-prefixed cell_id', 'chromosome', 'fragment_start', 'fragment_end'],
                        duplicate_policy='abort on duplicate key within chromosome/start bucket; no collapsing',
                        readSupport_signal_weight=False, additional_Tn5_shift=False,
                        denominator_names=list(DENOMINATORS), partition_order=list(PARTITIONS),
                        matrix_sum_may_count_fragment_in_multiple_overlapping_peaks=True,
                        native_RDS_measurement_relabelled=False, author_RDS_payload_read=False,
                        original_fragment_assignment_to_author_peaks_verified=False,
                        source_coordinate_convention='author-object 1-based closed; fixed receipt BED [start-1,end)',
                        external_export_filter_assembly_and_label_compatibility_verified=False,
                        no_normalization_fitting_scoring_or_selection=True, no_old_masks_folds_or_clinical_fields=True,
                        RNA_referenced_without_payload_read_or_copy=True, output_cap_bytes_inclusive=CAP)
        write_json(output / 'manifest.json', manifest)
        receipts = {}
        print(json.dumps(dict(phase='counting', cells=len(cells), units=len(units), files=len(rows), workers=8)), flush=True)
        with ThreadPoolExecutor(max_workers=8) as pool:
            futures = {pool.submit(run_unit, row, executable, output, offsets[0], arrays,
                                   membership_dir / f"unit_{row['unit_index']:03d}.tsv"): row for row in rows}
            try:
                for future in as_completed(futures):
                    receipt = future.result()
                    receipts[receipt['unit_index']] = receipt
                    check_cap(output)
                    print(json.dumps(dict(completed_units=len(receipts), unit_index=receipt['unit_index'],
                                          records=receipt['records'], elapsed_seconds=receipt['elapsed_seconds'])), flush=True)
            except BaseException:
                for future in futures:
                    future.cancel()  # Already running workers finish; no partial files are removed.
                raise
        require(len(receipts) == 272 and sum(r['records'] for r in receipts.values()) == RECORDS and
                sum(r['readSupport_diagnostic'][0] for r in receipts.values()) == READ_SUPPORT, 'Global source totals differ')
        hashes = verify_persisted(arrays, offsets[0], receipts, units)
        for receipt in receipts.values():
            require(fingerprint(receipt['source_path']) == receipt['source_stat_after'],
                    'Selected fragment source changed before final completion')
        for path, before in frozen_stats.items():
            require(fingerprint(path) == before and sha(path) == GUARDS[Path(path)], 'Frozen metadata changed during counting')
        for source, expected in source_hashes.items():
            require(sha(source) == expected, 'Executed helper/launcher source changed during counting')
        for path, recorded in tools.items():
            require(sha(path) == recorded['sha256'] and fingerprint(Path(path).resolve()) == recorded['stat'], 'Compiler/reader changed')
        for evidence in raw_headers:
            require(fingerprint(evidence['path']) == evidence['stat'], 'Original source header authority changed')
        for evidence in aligned_headers.values():
            require(fingerprint(evidence['path']) == evidence['stat'] and
                    npy_header(Path(evidence['path'])) == evidence['header'], 'Referenced aligned array header/stat changed')
        artifact = {}
        for path, digest in zip(arrays, hashes):
            final = Path(str(path).removesuffix('.partial'))
            require(not final.exists(), 'Refusing final array overwrite')
            path.rename(final)
            artifact[final.name] = dict(bytes=final.stat().st_size, computed_sha256=digest, header=npy_header(final))
        unit_receipts = []
        for unit in units:
            row = dict(unit)
            receipt = receipts.get(unit['unit_index'])
            for name in DENOMINATORS:
                for partition_index, partition in enumerate(PARTITIONS):
                    row[partition + '_' + name] = receipt[name][partition_index] if receipt else 0
            unit_receipts.append(row)
        write_tsv(output / 'units.tsv', tuple(unit_receipts[0]), unit_receipts)
        for path in (output / 'units.tsv', output / 'cell_identity.tsv', output / 'unit_identity.tsv', output / 'native_peak_geometry.tsv'):
            artifact[path.name] = dict(bytes=path.stat().st_size, computed_sha256=sha(path))
        write_json(output / 'artifacts.json', artifact)
        totals = {name: [sum(r[name][i] for r in receipts.values()) for i in range(3)] for name in DENOMINATORS}
        summary = dict(common_fragment_overlap_counts_completed=True, distinct_native_author_counts_endpoint=True,
                       donors=39, cells=68398, strata=273, sampled_strata=272, peaks=306706,
                       nuclei_partitions=[68398, sum(u['halfA_nuclei'] for u in units), sum(u['halfB_nuclei'] for u in units)],
                       partition_order=list(PARTITIONS), global_denominators=totals,
                       full_equals_halves_per_peak_and_unit=True, all_selected_compressed_hashes_verified=True,
                       all_gzip_crc_and_source_record_totals_verified=True, duplicate_keys=0,
                       output_arrays_dtype='<i4', biological_replication_unit='donor; nuclei and halves are not independent people',
                       native_RDS_count_unit_equivalence_claimed=False, external_compatibility_verified=False,
                       normalization_model_scoring_selection_run=False, elapsed_seconds=time.monotonic() - began,
                       output_cap_bytes_inclusive=CAP, saved_bytes_before_summary=saved_bytes(output),
                       final_inclusive_cap_checked=True)
        # The final summary is the atomic completion marker, written only after all checks.
        temporary = output / 'summary.json.partial'
        with temporary.open('x') as handle:
            json.dump(summary, handle, indent=2, allow_nan=False)
            handle.write('\n')
            handle.flush()
            os.fsync(handle.fileno())
        check_cap(output)  # Includes summary.partial; renaming does not change saved bytes.
        temporary.rename(output / 'summary.json')
        print(json.dumps(summary), flush=True)
    except BaseException:
        failure = dict(completed=False, error=traceback.format_exc(), partials_preserved=True,
                       saved_bytes_including_partials=saved_bytes(output))
        write_json(output / 'failure.json', failure)
        raise


if __name__ == '__main__':
    main()
