#!/usr/bin/env python3
"""Compute-only native-author RNA/ATAC count aggregation; no modeling."""
import csv
import gc
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
BENCH = ROOT / 'Analysis/MASLD_Model_Benchmark'
REC = BENCH / 'executions/codex-rec-20260929T142434Z'
HERE = ROOT / 'scripts/analysis/in_biopsy_allelic/rec_codex'
LOCK = ROOT / 'Analysis/SingleCell/continual_integration/reference/gse296875_source_lock_v44.json'
MEMBERSHIP = BENCH / 'executions/gse296875-fragment-membership-21063829/cell_membership.tsv.gz'
LABELS = ('B cells', 'Cholangiocytes', 'Hepatocytes', 'Kupffer', 'LSEC', 'Mesenchymal', 'NK-T')
LOCK_SHA = '6be686b4187309f038097cf73d3f2a8dd34e72d666bc92733a53e99a1ccbb620'
MEMBERSHIP_SHA = '81c43ee9af4a05784166c84a6b9ea0513f7b697ebbe33d9ff89982c52aa4e4c6'
CELL_SHA = 'a3d832f1128d2acea32a68d4c4287192a292695529722347aab1e4a07ce71217'
GEOMETRY_SHA = '50a52d097a5232edcf136dd53d6549ca2d1aa5dce9a57cf055ff10b7b2e0f3e3'
RDS_SHA = '6af74b7dbdd3839841f90cb7f3e47a3184b4dd3dd02242f6f3bdbd14533cfcab'
CAP = 2 * 1024**3
ATAC_TOTAL = 341091266
HALF_PREFIX = 'MASLD-NATIVE-HALVES-20261001\t'
PARTITIONS = ('full', 'halfA', 'halfB')


def require(value, message):
    if not value:
        raise ValueError(message)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def fingerprint(path):
    s = Path(path).stat()
    require(Path(path).is_file(), 'Missing regular source: ' + str(path))
    return dict(bytes=s.st_size, mtime_ns=s.st_mtime_ns, ctime_ns=s.st_ctime_ns,
                inode=s.st_ino, device=s.st_dev)


def saved_bytes(output):
    return sum(p.stat().st_size for p in output.rglob('*') if p.is_file())


def check_cap(output):
    require(saved_bytes(output) <= CAP, 'Inclusive saved output exceeds 2 GiB')


def write_json(output, name, value):
    target = output / name
    temporary = output / (name + '.partial')
    require(not target.exists() and not temporary.exists(), 'Refusing metadata overwrite')
    with temporary.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    check_cap(output)
    temporary.rename(target)


def write_tsv(output, name, fields, rows):
    target = output / name
    require(not target.exists(), 'Refusing TSV overwrite')
    with target.open('x', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    check_cap(output)


def strings(values):
    return [v.decode('utf-8') if isinstance(v, bytes) else str(v) for v in values]


def identities():
    require(sha(MEMBERSHIP) == MEMBERSHIP_SHA, 'Frozen identity membership differs')
    fields = ('well_id', 'raw_barcode', 'cell_id', 'donor_id', 'source_label')
    with gzip.open(MEMBERSHIP, 'rt', newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(set(fields) <= set(reader.fieldnames or ()), 'Identity membership schema differs')
        rows = [{name: row[name] for name in fields} for row in reader]
    require(len(rows) == len({r['cell_id'] for r in rows}) == 68398, 'Identity cell axis differs')
    require({r['source_label'] for r in rows} == set(LABELS), 'Native author labels differ')
    cell_hash = hashlib.sha256()
    donor_wells = {}
    support = {}
    for row in rows:
        require(re.fullmatch(r'well[1-8]', row['well_id']) is not None
                and re.fullmatch(r'[ACGT]+-[0-9]+', row['raw_barcode']) is not None
                and row['cell_id'] == row['well_id'] + '_' + row['raw_barcode'],
                'Identity well namespace or barcode differs')
        donor = row['donor_id']
        require(donor.isdecimal(), 'Non-numeric source donor ID')
        require(donor_wells.setdefault(donor, row['well_id']) == row['well_id'], 'Donor crosses wells')
        key = (donor, row['source_label'])
        support[key] = support.get(key, 0) + 1
        cell_hash.update(row['cell_id'].encode('utf-8') + b'\n')
    require(cell_hash.hexdigest() == CELL_SHA, 'Original native cell order differs')
    donors = sorted(donor_wells, key=int)
    require(len(donors) == 39, 'Donor count differs')
    units = [dict(unit_index=index, donor_id=donor, well_id=donor_wells[donor],
                  native_source_celltype=label, nuclei=support.get((donor, label), 0),
                  sampling_state='observed' if support.get((donor, label), 0) else 'unsampled')
             for index, (donor, label) in enumerate((d, l) for d in donors for l in LABELS)]
    require([(u['donor_id'], u['native_source_celltype']) for u in units if not u['nuclei']]
            == [('733', 'B cells')], 'Unsampled native stratum differs')
    require(sum(u['nuclei'] for u in units) == 68398, 'Identity assignment is not exhaustive')
    grouped = {}
    for row in rows:
        row['half_sort_sha256'] = hashlib.sha256((HALF_PREFIX + row['cell_id']).encode('utf-8')).hexdigest()
        grouped.setdefault((row['donor_id'], row['source_label']), []).append(row)
    for key, group in grouped.items():
        ordered = sorted(group, key=lambda r: (r['half_sort_sha256'], r['cell_id']))
        for i, row in enumerate(ordered):
            row['nucleus_half'] = 'A' if i < len(ordered) // 2 else 'B'
    for unit in units:
        n = unit['nuclei']
        unit.update(halfA_nuclei=n // 2, halfB_nuclei=n - n // 2,
                    halfA_sampling_state='observed' if n // 2 else 'unsampled',
                    halfB_sampling_state='observed' if n - n // 2 else 'unsampled',
                    halfA_missing=n // 2 == 0, halfB_missing=n - n // 2 == 0,
                    split_diagnostic_available=n >= 2)
    return rows, units


def aggregate_rna(path, rows, unit_index, expected_axis, lock, np, h5py, sparse):
    with h5py.File(path, 'r') as handle:
        matrix = handle['matrix']
        feature_types = strings(matrix['features/feature_type'][:])
        gene_positions = np.asarray([i for i, v in enumerate(feature_types) if v == 'Gene Expression'])
        ids = strings(matrix['features/id'][:])
        symbols = strings(matrix['features/name'][:])
        require(len(ids) == len(symbols) == len(feature_types), 'RNA feature metadata dimensions differ')
        gene_axis = [(ids[i], symbols[i]) for i in gene_positions]
        require(len(gene_axis) == len({g for g, _ in gene_axis}) == 36601, 'Native gene axis differs')
        axis_hash = hashlib.sha256()
        for identity, symbol in gene_axis:
            axis_hash.update(identity.encode('utf-8') + b'\0' + symbol.encode('utf-8') + b'\n')
        require(axis_hash.hexdigest() == lock['raw_rna_registry']['rna_feature_id_and_symbol_sha256'],
                'Native RNA identity+symbol SHA256 differs')
        require(expected_axis is None or gene_axis == expected_axis, 'Ordered RNA axes differ across wells')
        barcodes = strings(matrix['barcodes'][:])
        position = {b: i for i, b in enumerate(barcodes)}
        require(len(position) == len(barcodes), 'Duplicate raw RNA barcode')
        selected = np.asarray([position[r['raw_barcode']] for r in rows], dtype=np.int64)
        require(len(set(selected.tolist())) == len(rows), 'RNA cell assigned more than once')
        full_groups = np.asarray([unit_index[r['donor_id'], r['source_label']] for r in rows], dtype=np.int64)
        half_groups = full_groups + np.asarray([273 if r['nucleus_half'] == 'A' else 546 for r in rows])
        shape = tuple(int(v) for v in matrix['shape'][:])
        require(shape == (len(feature_types), len(barcodes)), '10x sparse dimensions differ')
        data_ds, index_ds, pointer_ds = (matrix[name] for name in ('data', 'indices', 'indptr'))
        require(data_ds.dtype.kind in 'iu' and index_ds.dtype.kind in 'iu'
                and pointer_ds.dtype.kind in 'iu', '10x count/index storage must be integer')
        require(len(data_ds.shape) == len(index_ds.shape) == len(pointer_ds.shape) == 1,
                '10x sparse arrays must be one-dimensional')
        require(len(data_ds) == len(index_ds) < 2**31, '10x sparse entry dimension differs')
        data = np.empty(len(data_ds), dtype=np.uint64)
        for first in range(0, len(data_ds), 1000000):
            part = data_ds[first:first + 1000000]
            require(not np.any(part < 0) and (not len(part) or int(part.max()) < 2**31),
                    'RNA entries violate finite nonnegative bounded integer storage')
            data[first:first + len(part)] = part
        pointers = pointer_ds[:]
        indices = index_ds[:]
        require(len(pointers) == shape[1] + 1 and int(pointers[0]) == 0
                and int(pointers[-1]) == len(data)
                and np.all(pointers[1:] >= pointers[:-1]), '10x CSC pointer layout differs')
        require(not len(indices) or (int(indices.min()) >= 0 and int(indices.max()) < shape[0]),
                '10x CSC row index outside native feature axis')
        source = sparse.csc_matrix((data, indices, pointers), shape=shape)
        require(source.has_canonical_format, '10x CSC entries are not canonical')
        selected_gene_counts = source[gene_positions, :][:, selected]
        group = sparse.csc_matrix((np.ones(2 * len(rows), dtype=np.uint64),
                                   (np.tile(np.arange(len(rows)), 2), np.concatenate((full_groups, half_groups)))),
                                  shape=(len(rows), 819))
        result = np.asarray((selected_gene_counts @ group).toarray(), dtype=np.uint64).T
        source_feature_totals = np.asarray(selected_gene_counts.sum(axis=1), dtype=np.uint64).ravel()
        require(np.array_equal(result[:273], result[273:546] + result[546:]), 'RNA full=A+B failed')
        require(np.array_equal(result[:273].sum(axis=0, dtype=np.uint64), source_feature_totals),
                'RNA per-feature conservation failed')
        source_cell_totals = np.asarray(selected_gene_counts.sum(axis=0), dtype=np.uint64).ravel()
        expected_group_totals = np.zeros(819, dtype=np.uint64)
        np.add.at(expected_group_totals, full_groups, source_cell_totals)
        np.add.at(expected_group_totals, half_groups, source_cell_totals)
        require(np.array_equal(result.sum(axis=1, dtype=np.uint64), expected_group_totals),
                'RNA per-group and raw denominator conservation failed')
        exact_total = sum(int(v) for v in source_feature_totals)
        require(exact_total < 2**63, 'RNA grand sum exceeds exact chosen aggregation bound')
        return result, gene_axis, dict(nuclei=len(rows), raw_barcodes=len(barcodes),
            native_genes=36601, selected_total_umis=exact_total,
            feature_conservation=True, group_conservation=True, full_equals_halfA_plus_halfB=True,
            data_storage=str(data_ds.dtype),
            source_sparse_entries=len(data_ds), exact_integer_storage_bound='nnz<2^31 and value<2^31')


def main():
    job = os.environ.get('SLURM_JOB_ID', '')
    require(job.isdecimal() and os.environ.get('SLURM_JOB_NODELIST'), 'Requires SLURM compute allocation')
    require(len(sys.argv) == 1, 'Fixed-source builder accepts no overrides')
    output = REC / ('codex_gse296875_aligned_native_pseudobulk_' + job)
    output.mkdir(mode=0o750, exist_ok=False)
    phase = 'archive_executed_sources'
    try:
        files = [Path(__file__).resolve(), HERE / 'codex_gse296875_aggregate_native_atac.R',
                 HERE / 'run_codex_gse296875_aligned_native_pseudobulk.sbatch']
        code = {}
        archive = output / 'executed_sources'
        archive.mkdir(mode=0o750)
        for path in files:
            before = fingerprint(path)
            code[str(path)] = sha(path)
            target = archive / path.name
            require(not target.exists(), 'Refusing executed source overwrite')
            shutil.copyfile(path, target)
            require(sha(target) == code[str(path)] and fingerprint(path) == before, 'Source changed during archive')
        require(sha(LOCK) == LOCK_SHA, 'Raw source lock differs')
        lock_stat = fingerprint(LOCK)
        lock = json.loads(LOCK.read_text())
        require(set(lock['raw_rna_sources']) == {'well' + str(i) for i in range(1, 9)}, 'Raw well roster differs')
        require(lock['processed_source']['sha256'] == RDS_SHA
                and lock['processed_source']['bytes'] == 8534583098, 'Author ATAC source contract differs')
        rds_path = (LOCK.parent / lock['processed_source']['path']).resolve(strict=True)
        rds_stat = fingerprint(rds_path)
        require(rds_stat['bytes'] == 8534583098, 'Full author RDS size differs')
        phase = 'frozen_native_identity_partition'
        rows, units = identities()
        membership_stat = fingerprint(MEMBERSHIP)
        unit_index = {(r['donor_id'], r['native_source_celltype']): r['unit_index'] for r in units}
        write_tsv(output, 'unit_identity.tsv', tuple(units[0]), units)
        write_tsv(output, 'cell_identity.tsv', ('well_id', 'raw_barcode', 'cell_id', 'donor_id', 'source_label',
                                               'half_sort_sha256', 'nucleus_half'), rows)
        by_well = {well: [r for r in rows if r['well_id'] == well] for well in sorted(lock['raw_rna_sources'])}
        reader_paths = {name: Path(importlib.util.find_spec(name).origin) for name in ('numpy', 'scipy', 'h5py')}
        reader_hashes = {name: dict(path=str(path), sha256=sha(path)) for name, path in reader_paths.items()}
        import numpy as np
        import scipy
        import scipy.sparse as sparse
        import h5py
        rscript = Path(shutil.which('Rscript') or '')
        require(rscript.is_file(), 'Installed Rscript is unavailable')
        inventory = dict(python=sys.version, executable=sys.executable,
            executable_sha256=sha(sys.executable), numpy=np.__version__, scipy=scipy.__version__,
            h5py=h5py.__version__, hdf5=h5py.version.hdf5_version, reader_entrypoints=reader_hashes,
            Rscript=str(rscript), Rscript_sha256=sha(rscript))
        manifest = dict(schema_version='codex-gse296875-aligned-native-counts-v1', job_id=job,
            scientific_role='project-exposed public development; data-only native author measurement strata',
            biological_replication_unit='source donor; nuclei and disjoint halves are not independent people',
            source_lock=dict(path=str(LOCK), sha256=LOCK_SHA),
            ATAC_source=dict(path=str(rds_path), recorded_sha256=RDS_SHA,
                             source_stat_before=rds_stat, full_byte_hash_computed_by_R_helper=True),
            membership=dict(path=str(MEMBERSHIP), sha256=MEMBERSHIP_SHA),
            native_cell_axis_sha256=CELL_SHA, source_code_hashes=code, resources=dict(cpus=1,memory_GiB=64,hours=90),
            argv=sys.argv, inventory=inventory, output_cap_bytes_inclusive=CAP,
            deterministic_halves=dict(hash_prefix=HALF_PREFIX, encoding='UTF8',
                sorting='within native donor/author-lineage by SHA256, then cell_id; first floor(n/2) A, rest B',
                same_assignment_for_RNA_and_ATAC=True, no_half_reliability_selection=True),
            no_minimum_nuclei_mask=True, no_clinical_fields=True, no_legacy_roles_or_folds=True,
            no_normalization_or_model_scoring=True, whole_RDS_deserialization_limit=
            'ATAC readRDS deserializes the entire author object; analytical access only peaks counts and dimnames.')
        write_json(output, 'manifest.json', manifest)
        phase = 'fresh_eight_well_raw_RNA_aggregation'
        counts = np.zeros((819, 36601), dtype=np.uint64)
        axis = None
        raw_receipts = {}
        source_stats = {}
        for well, record in sorted(lock['raw_rna_sources'].items()):
            path = (LOCK.parent / record['path']).resolve(strict=True)
            before = fingerprint(path)
            require(before['bytes'] == record['bytes'] and sha(path) == record['sha256'],
                    'Raw RNA byte SHA or size differs: ' + well)
            source_stats[str(path)] = before
            block, observed_axis, receipt = aggregate_rna(path, by_well[well], unit_index, axis, lock, np, h5py, sparse)
            axis = observed_axis
            # A donor occupies one well, so accumulated rows never overlap.
            touched = sorted({unit_index[r['donor_id'], r['source_label']] for r in by_well[well]})
            touched_all = touched + [i + 273 for i in touched] + [i + 546 for i in touched]
            require(not np.any(counts[touched_all]), 'RNA unit was assigned in multiple wells')
            counts += block
            receipt.update(path=str(path), recorded_sha256=record['sha256'], computed_sha256=record['sha256'],
                           source_stat_before=before, source_stat_after=fingerprint(path))
            require(fingerprint(path) == before, 'Raw RNA changed while aggregating: ' + well)
            raw_receipts[well] = receipt
            write_json(output, well + '_RNA_receipt.json', receipt)
            print(json.dumps(dict(well=well, raw_RNA_aggregated=True, nuclei=len(by_well[well]))), flush=True)
            del block
            gc.collect()
        require(sum(r['nuclei'] for r in raw_receipts.values()) == 68398, 'RNA cell conservation differs')
        total_rna = sum(r['selected_total_umis'] for r in raw_receipts.values())
        require(total_rna < 2**63 and sum(int(v) for v in counts[:273].sum(axis=1, dtype=np.uint64)) == total_rna,
                'RNA all-well total conservation failed')
        require(np.array_equal(counts[:273], counts[273:546] + counts[546:]), 'All-well RNA full=A+B failed')
        for offset in (0, 273, 546):
            require(not np.any(counts[offset + unit_index['733', 'B cells']]), 'Unsampled RNA row is not empty')
        rna_libraries = {}
        for index, partition in enumerate(PARTITIONS):
            values = counts[index * 273:(index + 1) * 273]
            rna_path = output / ('RNA_' + partition + '_counts_uint64.npy')
            rna_partial = output / (rna_path.name + '.partial')
            with rna_partial.open('xb') as handle:
                np.save(handle, values, allow_pickle=False)
            check_cap(output)
            persisted = np.load(rna_partial, mmap_mode='r', allow_pickle=False)
            require(np.array_equal(persisted, values), 'Persisted RNA counts differ')
            del persisted
            rna_partial.rename(rna_path)
            rna_libraries[partition] = [int(v) for v in values.sum(axis=1, dtype=np.uint64)]
        write_tsv(output, 'genes.tsv', ('gene_index', 'ensembl_id', 'native_gene_symbol'),
                  (dict(gene_index=i, ensembl_id=g, native_gene_symbol=s) for i, (g, s) in enumerate(axis)))
        del counts, values
        gc.collect()
        phase = 'native_author_ATAC_aggregation'
        run = subprocess.run([str(rscript), '--vanilla', str(archive / 'codex_gse296875_aggregate_native_atac.R'),
                              str(output)], check=False)
        require(run.returncode == 0, 'Native ATAC helper failed; its partials/errors retained')
        atac_receipt = json.loads((output / 'ATAC_receipt.json').read_text())
        require(atac_receipt['completed'] is True and atac_receipt['total_native_peak_events'] == ATAC_TOTAL,
                'ATAC completion or grand-total receipt differs')
        require(sha(output / 'native_peak_geometry.tsv') == GEOMETRY_SHA, 'Native ATAC geometry digest differs')
        phase = 'persisted_ATAC_NPY_halves_and_group_conservation'
        atac_arrays = {p: np.load(output / ('ATAC_' + p + '_counts_int32.npy'), mmap_mode='r', allow_pickle=False)
                       for p in PARTITIONS}
        require(all(a.shape == (273,306706) and a.dtype.str == '<i4' for a in atac_arrays.values()),
                'Persisted ATAC NPY axes or dtype differ')
        conserved = {p: np.zeros(273, dtype=np.int64) for p in PARTITIONS}
        for first in range(0, 306706, 8192):
            blocks = {p:a[:, first:first + 8192] for p,a in atac_arrays.items()}
            require(all(not np.any(b < 0) for b in blocks.values()), 'Negative persisted native ATAC count')
            require(np.array_equal(blocks['full'], blocks['halfA'] + blocks['halfB']),
                    'Persisted ATAC full=A+B failed')
            for partition, block in blocks.items():
                conserved[partition] += block.sum(axis=1, dtype=np.int64)
        for partition in PARTITIONS:
            expected = np.asarray(atac_receipt['unit_totals_' + partition], dtype=np.int64)
            require(np.array_equal(conserved[partition], expected), 'Persisted ATAC group conservation failed')
        require(int(conserved['full'].sum()) == ATAC_TOTAL, 'Persisted ATAC grand sum differs')
        del atac_arrays, blocks, block
        check_cap(output)
        units_final = []
        for i,row in enumerate(units):
            value = dict(row)
            for partition in PARTITIONS:
                value['RNA_' + partition + '_total_UMIs'] = rna_libraries[partition][i]
                value['ATAC_' + partition + '_total_native_peak_events'] = int(conserved[partition][i])
            units_final.append(value)
        write_tsv(output, 'units.tsv', tuple(units_final[0]), units_final)
        phase = 'final_source_identity_and_saved_output_checks'
        for path, state in source_stats.items():
            require(fingerprint(path) == state, 'Raw RNA source drifted before completion')
        require(fingerprint(rds_path) == rds_stat, 'Author RDS drifted before completion')
        require(fingerprint(LOCK) == lock_stat and sha(LOCK) == LOCK_SHA
                and fingerprint(MEMBERSHIP) == membership_stat and sha(MEMBERSHIP) == MEMBERSHIP_SHA,
                'Identity source drifted before completion')
        for path, expected in code.items():
            require(sha(path) == expected, 'Producer code changed during execution')
        for value in reader_hashes.values():
            require(sha(value['path']) == value['sha256'], 'Python reader entrypoint changed')
        require(sha(rscript) == inventory['Rscript_sha256'], 'Rscript changed during execution')
        artifacts = {}
        array_names = [assay + '_' + p + '_counts_' + dtype + '.npy'
                       for assay,dtype in (('RNA','uint64'),('ATAC','int32')) for p in PARTITIONS]
        for name in array_names + ['genes.tsv','native_peak_geometry.tsv','unit_identity.tsv','cell_identity.tsv','units.tsv']:
            path = output / name
            artifacts[name] = dict(bytes=path.stat().st_size, sha256=sha(path))
        headers = {}
        for name in array_names:
            with (output / name).open('rb') as handle:
                version = np.lib.format.read_magic(handle)
                shape, fortran, dtype = np.lib.format._read_array_header(handle, version)
                headers[name] = dict(version=list(version), shape=list(shape), fortran_order=fortran, dtype=dtype.str)
                header_bytes = handle.tell()
                handle.seek(0)
                headers[name]['header_sha256'] = hashlib.sha256(handle.read(header_bytes)).hexdigest()
        for name, header in headers.items():
            require(header['shape'] == ([273,36601] if name.startswith('RNA') else [273,306706])
                    and header['dtype'] == ('<u8' if name.startswith('RNA') else '<i4')
                    and header['fortran_order'] == name.startswith('ATAC'), 'Final NPY metadata differs')
        write_json(output, 'artifacts.json', dict(artifacts=artifacts, array_headers=headers))
        write_json(output, 'summary.json', dict(aligned_native_count_aggregation_completed=True,
            scientific_role=manifest['scientific_role'], same_nuclei=68398, donors=39, native_author_lineages=list(LABELS),
            strata=273, sampled_strata=272, unsampled_stratum=dict(donor_id='733',native_source_celltype='B cells'),
            RNA_genes=36601, ATAC_native_peaks=306706, RNA_total_UMIs=total_rna,
            ATAC_total_native_peak_events=ATAC_TOTAL, per_feature_and_group_conservation=True,
            full_equals_halfA_plus_halfB=True, deterministic_halves=manifest['deterministic_halves'],
            biological_replication_unit=manifest['biological_replication_unit'],
            no_minimum_nuclei_mask=True, clinical_fields_accessed=False, legacy_counts_or_folds_read=False,
            model_scoring_or_normalization_run=False, count_unit_ATAC='native author peak-count events',
            unique_fragment_or_molecule_units_verified=False, independently_confirmed_source_coordinate_origin=False,
            source_code_hashes=code, output_cap_bytes_inclusive=CAP,
            saved_bytes_before_final_summary=saved_bytes(output), final_inclusive_cap_checked=True))
        print(json.dumps(dict(aligned_native_count_aggregation_completed=True, output=str(output),
                              saved_bytes=saved_bytes(output))), flush=True)
    except Exception as error:
        failure = dict(aligned_native_count_aggregation_completed=False, phase=phase,
                       error=repr(error), traceback=traceback.format_exc(), partial_files_preserved=True)
        try:
            write_json(output, 'failure.json', failure)
        except Exception:
            traceback.print_exc()
        raise


if __name__ == '__main__':
    main()
