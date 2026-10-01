#!/usr/bin/env python3
"""Compute-only raw donor sums over all seven frozen native lineage buckets."""
import argparse
import array
import ast
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import struct
import sys
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
HERE = ROOT / 'scripts/analysis/in_biopsy_allelic/rec_codex'
RNA = REC / 'codex_gse296875_aligned_native_pseudobulk_21999585'
ATAC = REC / 'codex_gse296875_common_fragment_counts_22000534'
PRESPEC = REC / '20261001_codex_gse296875_donor_global_rna_atac_prespec.json'
LAUNCHER = HERE / 'run_codex_gse296875_aggregate_donor_global.sbatch'
PARTITIONS = ('full', 'halfA', 'halfB')
LABELS = ('B cells', 'Cholangiocytes', 'Hepatocytes', 'Kupffer', 'LSEC', 'Mesenchymal', 'NK-T')
DENOMS = ('primary_fragment_records', 'autosomal_fragment_records', 'union_overlapped_fragments',
          'readSupport_diagnostic', 'summed_peak_counts')
CAP = 512 * 1024**2
MAX_U64 = 2**64 - 1
GUARDS = {
    PRESPEC: '88b9e4239b08cc1ccbe78b092e8eacd395778453a496fd8915b4acc9e2729ae1',
    RNA / 'summary.json': 'dd1986b643e1e5ae2cec54d255fb8ed7619c6db668541194a4032277d8c00226',
    RNA / 'artifacts.json': '3c4102ecde36825237b7d31e6ff4c4395570dc97ee222e274429833cc000d8a5',
    RNA / 'units.tsv': '6ceeac2200ad1b452994e4df3fd6a4a60c51f812bcef4daa58116b84d4d69a1b',
    RNA / 'unit_identity.tsv': '5c0a1475cb5cbf71ae4e25f31652a822d4b773c4cd8f5f962bcc44b5bb80be68',
    RNA / 'genes.tsv': '2a33d4d339e2c6db695525ae2655af5156ccda008300bd992e306d48f01d38f0',
    ATAC / 'summary.json': 'bd64347bc066cbd182695c466b957ddd3c6efc689b53f8db2d557aed34321f8e',
    ATAC / 'artifacts.json': '8735af6f215fc82be928f2349cb8182d3d76e834304cce8eb32e858f5880a960',
    ATAC / 'manifest.json': '36d516811ca1d80f67ad6fdabc11af9d35734d7afdc6b09661a41769ee51ea4e',
    ATAC / 'units.tsv': 'babe4cd3bc0f4420c6303526d9ebbe342ae69ca792237a2c05b2bf09aa780f2a',
    ATAC / 'unit_identity.tsv': '509dadfb442445c28eb4139d2f3a9699d6734c82a3ac6442b29e2d767417a964',
}
for _directory in (RNA, ATAC):
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
    require(path.is_file() and not path.is_symlink(), 'Missing regular source: ' + str(path))
    info = path.stat()
    return dict(bytes=info.st_size, inode=info.st_ino, device=info.st_dev,
                mtime_ns=info.st_mtime_ns, ctime_ns=info.st_ctime_ns)


def saved_bytes(output):
    return sum(p.stat().st_size for p in output.rglob('*') if p.is_file())


def check_cap(output):
    require(saved_bytes(output) <= CAP, 'Inclusive saved output exceeds 512 MiB')


def write_json(output, name, value):
    final, partial = output / name, output / (name + '.partial')
    require(not final.exists() and not partial.exists(), 'Refusing JSON overwrite')
    with partial.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    check_cap(output)
    partial.rename(final)


def read_table(path, wanted):
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(set(wanted) <= set(reader.fieldnames or ()), 'Metadata schema differs: ' + str(path))
        return [{k: row[k] for k in wanted} for row in reader]


def write_table(output, name, rows):
    with (output / name).open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=tuple(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    check_cap(output)


def npy_header(path):
    with path.open('rb') as handle:
        require(handle.read(8) == b'\x93NUMPY\x01\x00', 'NPY magic/version differs')
        length = struct.unpack('<H', handle.read(2))[0]
        require(0 < length <= 65536, 'NPY header exceeds bound')
        parsed = ast.literal_eval(handle.read(length).decode('latin1').strip())
        require(set(parsed) == {'descr', 'shape', 'fortran_order'}, 'NPY header schema differs')
        offset = handle.tell()
        handle.seek(0)
        digest = hashlib.sha256(handle.read(offset)).hexdigest()
    return dict(dtype=parsed['descr'], shape=list(parsed['shape']), fortran_order=parsed['fortran_order'],
                offset=offset, header_sha256=digest)


def start_npy(handle, columns):
    body = repr(dict(descr='<u8', fortran_order=False, shape=(39, columns))).encode('latin1')
    body += b' ' * ((64 - (10 + len(body) + 1) % 64) % 64) + b'\n'
    handle.write(b'\x93NUMPY\x01\x00' + struct.pack('<H', len(body)) + body)
    require(handle.tell() == 128, 'Output header dimension differs')


def read_values(handle, columns, typecode):
    width = 8 if typecode == 'Q' else 4
    data = handle.read(columns * width)
    require(len(data) == columns * width, 'Count row truncated')
    values = array.array(typecode)
    values.frombytes(data)
    if sys.byteorder != 'little':
        values.byteswap()
    require(not values or min(values) >= 0, 'Negative source raw counts')
    return values


def bytes_little(values):
    if sys.byteorder == 'little':
        return values.tobytes()
    copy = array.array('Q', values)
    copy.byteswap()
    return copy.tobytes()


def donor_metadata(recipe):
    identity = ('unit_index', 'donor_id', 'well_id', 'native_source_celltype', 'nuclei', 'halfA_nuclei', 'halfB_nuclei')
    rna_fields = identity + tuple('RNA_' + p + '_total_UMIs' for p in PARTITIONS)
    atac_fields = identity + tuple(p + '_' + name for p in PARTITIONS for name in DENOMS)
    rna_units = read_table(RNA / 'units.tsv', rna_fields)
    atac_units = read_table(ATAC / 'units.tsv', atac_fields)
    require(len(rna_units) == len(atac_units) == 273 and
            [[u[k] for k in identity] for u in rna_units] == [[u[k] for k in identity] for u in atac_units],
            'Native RNA/ATAC ordered unit join differs')
    donors = sorted({u['donor_id'] for u in rna_units}, key=int)
    require(len(donors) == 39 and [(u['donor_id'], u['native_source_celltype']) for u in rna_units] ==
            [(d, label) for d in donors for label in LABELS] and
            [int(u['unit_index']) for u in rna_units] == list(range(273)), 'Exact full Cartesian unit axis differs')
    require([(u['donor_id'], u['native_source_celltype']) for u in rna_units if int(u['nuclei']) == 0] == [('733', 'B cells')],
            'Empty native source bucket differs')
    require(sum(int(u['nuclei']) for u in rna_units) == 68398, 'Native nucleus roster differs')
    manifest = json.loads((ATAC / 'manifest.json').read_text())
    source_rows = {r['unit_index']: r for r in manifest['selected_fragment_inputs']}
    receipt_guards, summaries, fractions = {}, [], []
    for donor_index, donor in enumerate(donors):
        group = rna_units[donor_index * 7:(donor_index + 1) * 7]
        wells = {u['well_id'] for u in group}
        require(len(wells) == 1 and donor in recipe['source_identity']['well_donor_roster'][next(iter(wells))],
                'Donor/well identity differs')
        result = dict(donor_index=donor_index, donor_id=donor, well_id=next(iter(wells)),
                      stored_lineage_buckets=7, sampled_lineage_buckets=sum(int(u['nuclei']) > 0 for u in group))
        for partition, nucleus_field in zip(PARTITIONS, ('nuclei', 'halfA_nuclei', 'halfB_nuclei')):
            result[partition + '_nuclei'] = sum(int(u[nucleus_field]) for u in group)
            result[partition + '_RNA_total_UMIs'] = sum(int(u['RNA_' + partition + '_total_UMIs']) for u in group)
            for name in DENOMS:
                result[partition + '_' + name] = sum(int(atac_units[int(u['unit_index'])][partition + '_' + name]) for u in group)
            result[partition + '_RNA_input_defined'] = result[partition + '_RNA_total_UMIs'] > 0
            result[partition + '_ATAC_target_denominator_defined'] = result[partition + '_autosomal_fragment_records'] > 0
        for name in ('nuclei', 'RNA_total_UMIs') + DENOMS:
            require(result['full_' + name] == result['halfA_' + name] + result['halfB_' + name], 'Donor raw metadata full=A+B differs')
        for unit in group:
            i = int(unit['unit_index'])
            row = dict(donor_index=donor_index, donor_id=donor, well_id=result['well_id'],
                       unit_index=i, native_source_celltype=unit['native_source_celltype'],
                       empty_source_roster=int(unit['nuclei']) == 0, imputed_measurement=False)
            for partition, field in zip(PARTITIONS, ('nuclei', 'halfA_nuclei', 'halfB_nuclei')):
                n, total = int(unit[field]), result[partition + '_nuclei']
                row[partition + '_nuclei'] = n
                row[partition + '_nucleus_fraction'] = n / total if total else None
                if n == 0:
                    require(int(unit['RNA_' + partition + '_total_UMIs']) == 0 and
                            all(int(atac_units[i][partition + '_' + name]) == 0 for name in DENOMS),
                            'Empty source roster has nonempty measurement metadata')
            fractions.append(row)
            if int(unit['nuclei']) > 0:
                path = ATAC / 'unit_receipts' / f'unit_{i:03d}_verified.json'
                before = fingerprint(path)
                receipt = json.loads(path.read_text())
                require(receipt['unit_index'] == i and receipt['seen_cells'] == int(unit['nuclei']) and
                        receipt['computed_compressed_sha256'] == source_rows[i]['sha256'] and
                        receipt['gzip_crc_and_eof_verified'] and receipt['full_equals_halves_per_peak'] and receipt['duplicate_keys'] == 0,
                        'Sampled bucket source receipt missing or inconsistent')
                for name in DENOMS:
                    require(receipt[name] == [int(atac_units[i][p + '_' + name]) for p in PARTITIONS], 'Bucket denominator receipt differs')
                receipt_guards[str(path)] = dict(sha256=sha(path), stat=before)
                require(fingerprint(path) == before, 'Bucket receipt changed while checking')
        summaries.append(result)
    require(len(receipt_guards) == 272, 'Incomplete sampled bucket receipts')
    return rna_units, atac_units, summaries, fractions, receipt_guards


def aggregate(source, target, columns, typecode, partition, modality, units, donors):
    row_totals = []
    with source.open('rb') as input_handle, target.open('xb') as output_handle:
        input_handle.seek(128)
        start_npy(output_handle, columns)
        for donor in donors:
            accumulator = array.array('Q', [0]) * columns
            donor_index = donor['donor_index']
            for unit in units[donor_index * 7:(donor_index + 1) * 7]:
                values = read_values(input_handle, columns, typecode)
                expected = int(unit['RNA_' + partition + '_total_UMIs']) if modality == 'RNA' else int(unit[partition + '_summed_peak_counts'])
                require(sum(values) == expected and expected <= MAX_U64, 'Persisted source bucket row total differs')
                if int(unit['nuclei' if partition == 'full' else partition + '_nuclei']) == 0:
                    require(not any(values), 'Empty source roster has nonzero count payload')
                for feature, value in enumerate(values):
                    previous = accumulator[feature]
                    if value > MAX_U64 - previous:
                        raise OverflowError('Raw donor feature sum exceeds uint64')
                    accumulator[feature] = previous + value
            total = sum(accumulator)
            expected = donor[partition + '_RNA_total_UMIs'] if modality == 'RNA' else donor[partition + '_summed_peak_counts']
            require(total == expected and total <= MAX_U64, 'Raw donor feature/metadata sum differs')
            row_totals.append(total)
            output_handle.write(bytes_little(accumulator))
        require(not input_handle.read(1), 'Trailing source NPY data')
        output_handle.flush()
        os.fsync(output_handle.fileno())
    require(target.stat().st_size == 128 + 39 * columns * 8, 'Aggregate output byte dimension differs')
    return row_totals


def verify_outputs(paths, columns, donors, modality):
    handles = [p.open('rb') for p in paths]
    digests = [hashlib.sha256() for _ in paths]
    try:
        for handle, digest in zip(handles, digests):
            digest.update(handle.read(128))
        for donor in donors:
            rows = []
            for handle, digest in zip(handles, digests):
                data = handle.read(columns * 8)
                require(len(data) == columns * 8, 'Persisted donor row truncated')
                digest.update(data)
                values = array.array('Q')
                values.frombytes(data)
                if sys.byteorder != 'little':
                    values.byteswap()
                rows.append(values)
            require(all(f == a + b and a <= MAX_U64 - b for f, a, b in zip(*rows)), 'Persisted donor per-feature full=A+B differs')
            for partition, values in zip(PARTITIONS, rows):
                expected = donor[partition + '_RNA_total_UMIs'] if modality == 'RNA' else donor[partition + '_summed_peak_counts']
                require(sum(values) == expected, 'Persisted donor row/library totals differ')
        require(all(not h.read(1) for h in handles), 'Trailing aggregate array payload')
    finally:
        for handle in handles:
            handle.close()
    return [d.hexdigest() for d in digests]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true')
    parser.add_argument('--expected-helper-sha')
    args = parser.parse_args()
    if args.inspect:
        print(json.dumps(dict(role='raw donor sums only', donors=39, RNA_genes=36601, ATAC_peaks=306706,
                              largest_accumulator_bytes=306706 * 8, persisted_three_row_bytes=3 * 306706 * 8,
                              estimated_process_memory_GiB_upper=0.2, six_array_bytes_including_headers=321336120,
                              cpu=1, memory_GiB=2, hours=1, inclusive_cap_bytes=CAP), indent=2))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdecimal() and os.environ.get('SLURM_CPUS_PER_TASK') == '1', 'Requires one-CPU SLURM compute job')
    require(args.expected_helper_sha and sha(Path(__file__)) == args.expected_helper_sha, 'Reviewed aggregation helper drift')
    require(array.array('Q').itemsize == 8 and array.array('i').itemsize == 4, 'Installed integer reader widths differ')
    output = REC / ('codex_gse296875_donor_global_raw_' + os.environ['SLURM_JOB_ID'])
    output.mkdir(mode=0o750)
    began = time.monotonic()
    try:
        archive = output / 'executed_sources'
        archive.mkdir()
        source_pins = {}
        for source in (Path(__file__), LAUNCHER):
            value = source.read_bytes()
            (archive / source.name).write_bytes(value)
            source_pins[str(source)] = hashlib.sha256(value).hexdigest()
        require(source_pins[str(Path(__file__))] == args.expected_helper_sha, 'Source changed before archive')
        snapshots = {}
        for path, expected in GUARDS.items():
            before = fingerprint(path)
            require(sha(path) == expected and fingerprint(path) == before, 'Frozen metadata/source axis differs: ' + str(path))
            snapshots[str(path)] = before
        recipe = json.loads(PRESPEC.read_text())
        require(recipe['authorization']['aggregation_authorized'] and not recipe['authorization']['fitting_authorized'], 'Recipe aggregation authority differs')
        require(json.loads((RNA / 'summary.json').read_text())['aligned_native_count_aggregation_completed'] and
                json.loads((ATAC / 'summary.json').read_text())['common_fragment_overlap_counts_completed'], 'Source completion markers differ')
        rna_units, atac_units, donors, fractions, receipts = donor_metadata(recipe)
        genes = read_table(RNA / 'genes.tsv', ('gene_index', 'ensembl_id', 'native_gene_symbol'))
        require(len(genes) == len({g['ensembl_id'] for g in genes}) == 36601 and
                [int(g['gene_index']) for g in genes] == list(range(36601)), 'Complete native gene axis differs')
        gene_hash = hashlib.sha256()
        for gene in genes:
            gene_hash.update(gene['ensembl_id'].encode('utf-8') + b'\0' + gene['native_gene_symbol'].encode('utf-8') + b'\n')
        require(gene_hash.hexdigest() == recipe['source_identity']['rna_gene_id_symbol_order_sha256'], 'Native ordered gene identity/symbol codec differs')
        del genes
        rna_art = json.loads((RNA / 'artifacts.json').read_text())
        atac_art = json.loads((ATAC / 'artifacts.json').read_text())
        inputs = []
        for modality, directory, columns, dtype in (('RNA', RNA, 36601, '<u8'), ('ATAC', ATAC, 306706, '<i4')):
            for partition in PARTITIONS:
                name = f'RNA_{partition}_counts_uint64.npy' if modality == 'RNA' else f'ATAC_fragment_overlap_{partition}_int32.npy'
                path = directory / name
                artifact = rna_art['artifacts'][name] if modality == 'RNA' else atac_art[name]
                expected = artifact['sha256'] if modality == 'RNA' else artifact['computed_sha256']
                before = fingerprint(path)
                parsed = npy_header(path)
                require(parsed['shape'] == [273, columns] and parsed['dtype'] == dtype and parsed['fortran_order'] is False
                        and parsed['offset'] == 128 and before['bytes'] == artifact['bytes'] == 128 + 273 * columns * (8 if modality == 'RNA' else 4),
                        'Source raw count header/size differs')
                recorded = rna_art['array_headers'][name] if modality == 'RNA' else artifact['header']
                require(parsed['header_sha256'] == recorded['header_sha256'] and sha(path) == expected
                        and fingerprint(path) == before, 'Source raw array SHA/header/stat differs')
                inputs.append(dict(path=str(path), modality=modality, partition=partition, columns=columns,
                                   sha256=expected, stat_before=before, header=parsed))
        readers = {str(Path(sys.executable).resolve()): sha(Path(sys.executable).resolve()),
                   str(Path(array.__file__).resolve()): sha(Path(array.__file__).resolve())}
        write_json(output, 'manifest.json', dict(source_pins=source_pins, prespec_sha256=GUARDS[PRESPEC],
                   metadata_pins={str(p):v for p,v in GUARDS.items()}, input_arrays=inputs,
                   sampled_bucket_receipts=receipts, source_stat_before=snapshots, argv=sys.argv, python_version=sys.version,
                   stdlib_reader_hashes=readers, resources=dict(cpus=1,memory_GiB=2,hours=1),
                   aggregation='Raw integer sums across seven disjoint native buckets in numeric donor order, before all transforms.',
                   empty_source_roster='733/B cells contributes empty roster, not an imputed biological measurement.',
                   nucleus_fraction_metadata_only=True, no_filter_masks_equal_lineage_weights_or_clinical_fields=True,
                   no_transforms_model_fits_or_external_reads=True, output_cap_bytes_inclusive=CAP))
        outputs = {}
        for modality, columns, typecode, units in (('RNA', 36601, 'Q', rna_units), ('ATAC', 306706, 'i', atac_units)):
            partials = []
            for partition in PARTITIONS:
                source = next(row for row in inputs if row['modality'] == modality and row['partition'] == partition)
                partial = output / f'{modality}_donor_{partition}_counts_uint64.npy.partial'
                aggregate(Path(source['path']), partial, columns, typecode, partition, modality, units, donors)
                partials.append(partial)
                check_cap(output)
                print(json.dumps(dict(modality=modality, partition=partition, raw_donors_aggregated=39)), flush=True)
            digests = verify_outputs(partials, columns, donors, modality)
            for partial, digest in zip(partials, digests):
                parsed = npy_header(partial)
                require(parsed['shape'] == [39, columns] and parsed['dtype'] == '<u8' and parsed['fortran_order'] is False and parsed['offset'] == 128,
                        'Persisted donor array header differs')
                outputs[partial.name.removesuffix('.partial')] = dict(bytes=partial.stat().st_size, sha256=digest, header=parsed)
        for source in inputs:
            path = Path(source['path'])
            require(fingerprint(path) == source['stat_before'] and sha(path) == source['sha256']
                    and fingerprint(path) == source['stat_before'], 'Input count source changed while aggregating')
        for path, expected in GUARDS.items():
            require(fingerprint(path) == snapshots[str(path)] and sha(path) == expected and
                    fingerprint(path) == snapshots[str(path)], 'Frozen metadata changed while aggregating')
        for path, receipt in receipts.items():
            source = Path(path)
            require(fingerprint(source) == receipt['stat'] and sha(source) == receipt['sha256'], 'Sampled bucket receipt changed')
        for source, expected in {**source_pins, **readers}.items():
            require(sha(Path(source)) == expected, 'Executed helper or installed integer reader changed')
        for name in outputs:
            require(not (output / name).exists(), 'Refusing final array overwrite')
            (output / (name + '.partial')).rename(output / name)
        for source in (RNA / 'genes.tsv', RNA / 'native_peak_geometry.tsv'):
            require(not (output / source.name).exists(), 'Refusing feature identity overwrite')
            shutil.copyfile(source, output / source.name)
            require(sha(output / source.name) == GUARDS[source], 'Copied feature identity differs')
        write_table(output, 'donors.tsv', donors)
        write_table(output, 'lineage_nucleus_fractions.tsv', fractions)
        for name in ('genes.tsv', 'native_peak_geometry.tsv', 'donors.tsv', 'lineage_nucleus_fractions.tsv'):
            outputs[name] = dict(bytes=(output / name).stat().st_size, sha256=sha(output / name))
        write_json(output, 'artifacts.json', outputs)
        totals = {p: {name:sum(d[p + '_' + name] for d in donors) for name in ('nuclei','RNA_total_UMIs') + DENOMS} for p in PARTITIONS}
        require([totals[p]['nuclei'] for p in PARTITIONS] == [68398,34136,34262], 'All-donor nucleus partitions differ')
        require([totals[p]['RNA_total_UMIs'] for p in PARTITIONS] == [610383995,306227813,304156182], 'All-donor raw RNA conservation differs')
        require([totals[p]['summed_peak_counts'] for p in PARTITIONS] == [341091266,170071577,171019689], 'All-donor raw ATAC conservation differs')
        source_denoms = json.loads((ATAC / 'summary.json').read_text())['global_denominators']
        require(all(totals[p][name] == source_denoms[name][i] for i,p in enumerate(PARTITIONS) for name in DENOMS),
                'All-donor denominator families do not conserve source totals')
        write_json(output, 'summary.json', dict(raw_donor_global_aggregation_completed=True, donors=39, stored_source_buckets=273,
                   sampled_source_buckets=272, source_nuclei=68398, RNA_genes=36601, ATAC_peaks=306706,
                   count_dtype='<u8', partition_order=PARTITIONS, global_raw_metadata=totals,
                   persisted_per_donor_per_feature_full_equals_halves=True, source_and_output_row_global_totals_reconciled=True,
                   all_source_hashes_verified_before_after=True, no_imputation_or_source_bucket_selection=True,
                   full_donor_RNA_and_autosomal_denominators_defined=all(d['full_RNA_input_defined'] and d['full_ATAC_target_denominator_defined'] for d in donors),
                   normalization_model_fits_or_external_reads_run=False, elapsed_seconds=time.monotonic()-began,
                   output_cap_bytes_inclusive=CAP, saved_bytes_before_summary=saved_bytes(output), final_inclusive_cap_checked=True))
        print(json.dumps(dict(raw_donor_global_aggregation_completed=True, donors=39, saved_bytes=saved_bytes(output))), flush=True)
    except BaseException:
        write_json(output, 'failure.json', dict(completed=False, traceback=traceback.format_exc(), partials_preserved=True))
        raise


if __name__ == '__main__':
    main()
