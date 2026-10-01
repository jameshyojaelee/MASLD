#!/usr/bin/env python3
"""Full GSE296875 identity census; molecular datasets remain closed."""
import argparse
from collections import Counter, defaultdict
import csv
import hashlib
import itertools
import json
import os
from pathlib import Path
import re
import sys
import traceback

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
PREP = REC / 'codex_gse296875_metadata_census_preparation_v2'
SCRIPT = 'codex_gse296875_metadata_census.py'
LAUNCHER = 'run_codex_gse296875_metadata_census.sbatch'
OUTPUT_CAP = 100 * 1024 * 1024
IDENTITY_PATHS = ('matrix/features/id', 'matrix/features/name',
                  'matrix/features/feature_type', 'matrix/barcodes',
                  'rna/ensembl_id', 'rna/gene_name', 'atac/peak_id',
                  'atac/chromosome', 'atac/source_start_1based_closed',
                  'atac/source_end_1based_closed', 'atac/bed_start_0based',
                  'atac/bed_end_half_open')


def require(test, message):
    if not test:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            digest.update(chunk)
    return digest.hexdigest()


def json_write(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, allow_nan=False)
        handle.write('\n')


def write_tsv(path, fields, rows):
    with path.open('x', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t')
        writer.writerow(fields)
        writer.writerows(rows)


def exact_rows(path, columns):
    with Path(path).open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        require(tuple(reader.fieldnames or ()) == tuple(columns), 'Identity columns changed: ' + str(path))
        for row in reader:
            require(None not in row and all(row[name] is not None for name in columns),
                    'Malformed identity record: ' + str(path))
            yield row


def text(value):
    return value.decode('utf-8') if isinstance(value, bytes) else str(value)


def identity_dataset(handle, name):
    require(name in IDENTITY_PATHS, 'Non-identity HDF5 dataset forbidden: ' + name)
    require(name in handle, 'Missing identity dataset: ' + name)
    dataset = handle[name]
    require(len(dataset.shape) == 1, 'Identity axis must be one dimensional: ' + name)
    return dataset


def strings(handle, name, start=None, stop=None):
    dataset = identity_dataset(handle, name)
    require(dataset.dtype.kind in ('O', 'S', 'U'), 'Non-string identity dataset: ' + name)
    return [text(value) for value in dataset[slice(start, stop)]]


def native_rna_axis(handle):
    ds = [identity_dataset(handle, 'matrix/features/' + name)
          for name in ('id', 'name', 'feature_type')]
    require(len({value.shape for value in ds}) == 1, 'Native feature metadata lengths differ')
    ids, names, types = [], [], Counter()
    for start in range(0, len(ds[0]), 8192):
        part = [strings(handle, 'matrix/features/' + name, start, start + 8192)
                for name in ('id', 'name', 'feature_type')]
        for feature, name, kind in zip(*part):
            types[kind] += 1
            if kind == 'Gene Expression':
                ids.append(feature)
                names.append(name)
    require(len(ids) == len(set(ids)) == 36601, 'Full native RNA feature identity differs')
    require(all(re.fullmatch(r'ENSG[0-9]+(?:\.[0-9]+)?', gene) for gene in ids),
            'Native RNA Ensembl identity malformed')
    registry = hashlib.sha256()
    for identity, symbol in zip(ids, names):
        registry.update(identity.encode())
        registry.update(b'\0')
        registry.update(symbol.encode())
        registry.update(b'\n')
    digest = registry.hexdigest()
    return ids, names, digest, dict(types)


def census(output, mapping, labels):
    maps, labs = {}, {}
    mapping_ids, label_ids, pair_ids = Counter(), Counter(), Counter()
    mapping_rows = label_rows = 0
    empty_mapping = empty_labels = invalid_barcodes = cell_barcode_mismatch = 0
    for row in exact_rows(mapping, ('cell_id', 'donor_id', 'well_id', 'raw_barcode')):
        mapping_rows += 1
        mapping_ids[row['cell_id']] += 1
        pair_ids[(row['well_id'], row['raw_barcode'])] += 1
        empty_mapping += not all(row.values())
        maps.setdefault(row['cell_id'], row)
    for row in exact_rows(labels, ('cell_id', 'author_label')):
        label_rows += 1
        label_ids[row['cell_id']] += 1
        empty_labels += not all(row.values())
        labs.setdefault(row['cell_id'], row['author_label'])
    mismatches = dict(mapping_rows=mapping_rows, label_rows=label_rows,
        duplicate_mapping_cell_ids=sum(n - 1 for n in mapping_ids.values()),
        duplicate_label_cell_ids=sum(n - 1 for n in label_ids.values()),
        duplicate_well_barcode_keys=sum(n - 1 for n in pair_ids.values()),
        empty_mapping_rows=empty_mapping, empty_label_rows=empty_labels,
        mapping_without_label=len(set(maps) - set(labs)),
        labels_without_mapping=len(set(labs) - set(maps)))
    write_tsv(output / 'join_mismatch_cell_ids.tsv', ('cell_id', 'issue'),
        itertools.chain(((key, 'mapping_without_label') for key in sorted(set(maps) - set(labs))),
                        ((key, 'labels_without_mapping') for key in sorted(set(labs) - set(maps))),
                        ((key, 'duplicate_mapping_cell_id') for key, n in sorted(mapping_ids.items()) if n > 1),
                        ((key, 'duplicate_label_cell_id') for key, n in sorted(label_ids.items()) if n > 1)))
    joined, expected_barcodes = [], defaultdict(set)
    for key in sorted(set(maps) & set(labs)):
        row = maps[key]
        prefix = row['well_id'] + '_'
        encoded = row['raw_barcode']
        if row['cell_id'] != encoded:
            cell_barcode_mismatch += 1
            continue
        if not encoded.startswith(prefix) or re.fullmatch(r'[ACGT]+-[0-9]+', encoded[len(prefix):]) is None:
            invalid_barcodes += 1
            continue
        bare = encoded[len(prefix):]
        expected_barcodes[row['well_id']].add(bare)
        joined.append((key, row['donor_id'], row['well_id'], bare, labs[key]))
    mismatches['invalid_namespaced_barcode_rows'] = invalid_barcodes
    mismatches['cell_id_vs_namespaced_barcode_mismatch_rows'] = cell_barcode_mismatch
    json_write(output / 'join_diagnostics.json', mismatches)
    require(mapping_rows == label_rows == len(joined) == 68398,
            'Full source census differs from locked 68,398 nuclei; diagnostics retained')
    require(all(value == 0 for key, value in mismatches.items() if key not in ('mapping_rows', 'label_rows')),
            'Duplicate, missing or malformed source identity; diagnostics retained')
    donors = sorted({row[1] for row in joined})
    wells = sorted({row[2] for row in joined})
    native_labels = sorted({row[4] for row in joined})
    require(len(donors) == 39 and wells == ['well' + str(i) for i in range(1, 9)],
            'Full source donor/well identity differs')
    support = Counter((d, label, well) for _, d, well, _, label in joined)
    donor_wells = {d: sorted({w for _, donor, w, _, _ in joined if donor == d}) for d in donors}
    well_donors = {w: sorted({d for _, d, well, _, _ in joined if well == w}) for w in wells}
    bare_wells = defaultdict(set)
    for _, _, well, bare, _ in joined:
        bare_wells[bare].add(well)
    write_tsv(output / 'donor_native_celltype_well_support.tsv',
        ('donor_id', 'native_source_celltype', 'well_id', 'nuclei'),
        ((d, label, w, support[d, label, w]) for d in donors for label in native_labels for w in wells))
    write_tsv(output / 'donor_native_celltype_support.tsv',
        ('donor_id', 'native_source_celltype', 'nuclei', 'observed_wells'),
        ((d, label, sum(support[d, label, w] for w in wells),
          ';'.join(w for w in wells if support[d, label, w])) for d in donors for label in native_labels))
    write_tsv(output / 'well_pair_donor_overlap.tsv', ('well_a', 'well_b', 'shared_donors', 'donor_ids'),
        ((a, b, len(set(well_donors[a]) & set(well_donors[b])),
          ';'.join(sorted(set(well_donors[a]) & set(well_donors[b])))) for a in wells for b in wells))
    result = dict(nuclei=len(joined), biological_donors=39, donors=donors, wells=wells,
        native_source_celltypes=native_labels, donor_wells=donor_wells, well_donors=well_donors,
        donors_spanning_multiple_wells=[d for d in donors if len(donor_wells[d]) > 1],
        bare_barcodes_shared_across_wells=sum(len(v) > 1 for v in bare_wells.values()),
        donor_celltype_support_is_sampling_only=True, wells_are_technical_batches=True,
        borrowed_RNA_labels_are_not_ATAC_prediction_accuracy=True)
    return result, expected_barcodes


def axes(output, contract, lock, expected_barcodes):
    import h5py
    sources, base_axis, barcode_mismatch = {}, None, False
    for well in sorted(lock['raw_rna_sources']):
        receipt = lock['raw_rna_sources'][well]
        path = (Path(contract['source_lock']['path']).parent / receipt['path']).resolve()
        if not path.is_file():
            raise FileNotFoundError('Full raw source identity artifact missing: ' + str(path))
        require(path.stat().st_size == receipt['bytes'], 'Full raw H5 size differs: ' + well)
        before = fingerprint(path)
        with h5py.File(path, 'r') as handle:
            ids, names, digest, types = native_rna_axis(handle)
            require(digest == lock['raw_rna_registry']['rna_feature_id_and_symbol_sha256'],
                    'Native RNA feature identity hash differs: ' + well)
            if base_axis is None:
                base_axis = (ids, names)
                write_tsv(output / 'native_rna_gene_axis.tsv', ('ensembl_id', 'native_gene_symbol'), zip(ids, names))
            else:
                require((ids, names) == base_axis, 'Native RNA axes differ across wells')
            ds = identity_dataset(handle, 'matrix/barcodes')
            seen, found, duplicates = set(), set(), 0
            for start in range(0, len(ds), 65536):
                for barcode in strings(handle, 'matrix/barcodes', start, start + 65536):
                    duplicates += barcode in seen
                    seen.add(barcode)
                    if barcode in expected_barcodes[well]:
                        found.add(barcode)
            absent = sorted(expected_barcodes[well] - found)
            write_tsv(output / (well + '_missing_raw_barcode_identity.tsv'), ('raw_10x_barcode',), ((v,) for v in absent))
            sources[well] = dict(path=str(path), actual_bytes=path.stat().st_size,
                declared_source_sha256=receipt['sha256'], whole_H5_byte_hash_checked=False,
                RNA_gene_axis_sha256=digest, RNA_gene_axis_hash_codec='UTF8_id_NUL_UTF8_symbol_LF_per_row',
                native_features_by_type=types,
                raw_barcode_rows=len(ds), duplicate_raw_barcodes=duplicates,
                mapped_nuclei=len(expected_barcodes[well]), missing_mapped_barcodes=len(absent))
            barcode_mismatch |= bool(duplicates or absent)
        require(fingerprint(path) == before, 'Raw H5 changed during identity inspection')
        print(json.dumps(dict(well=well, identity_metadata_checked=True)), flush=True)
    json_write(output / 'native_raw_identity_sources.json', sources)
    require(not barcode_mismatch, 'Raw barcode duplicates or missing mapped cells; diagnostics retained')
    path = Path(contract['smoke_h5']['path'])
    if not path.is_file():
        raise FileNotFoundError('Custom merged ATAC axis artifact missing: ' + str(path))
    require(path.stat().st_size == contract['smoke_h5']['bytes'], 'Smoked full-axis H5 size differs')
    before = fingerprint(path)
    with h5py.File(path, 'r') as handle:
        require(strings(handle, 'rna/ensembl_id') == base_axis[0]
                and strings(handle, 'rna/gene_name') == base_axis[1], 'Smoke/native RNA identities differ')
        fields = ('peak_id', 'chromosome', 'source_start_1based_closed',
                  'source_end_1based_closed', 'bed_start_0based', 'bed_end_half_open')
        ds = [identity_dataset(handle, 'atac/' + name) for name in fields]
        require(all(value.shape == (306706,) for value in ds), 'Full custom merged peak axis differs')
        require(text(handle['atac'].attrs['source_coordinate_system']) == 'GRanges_1_based_closed'
                and text(handle['atac'].attrs['source_scale']) == 'raw_integer_peak_counts',
                'Source peak coordinate convention differs')
        peaks, contigs, invalid = set(), Counter(), 0
        primary_contigs = {'chr' + str(i) for i in range(1, 23)} | {'chrX', 'chrY'}
        with (output / 'native_custom_merged_atac_geometry.tsv').open('x', newline='') as target:
            writer = csv.writer(target, delimiter='\t'); writer.writerow(fields)
            for start in range(0, 306706, 8192):
                parts = [strings(handle, 'atac/' + name, start, start + 8192) for name in fields[:2]]
                for dataset in ds[2:]:
                    require(dataset.dtype.kind in ('i', 'u'), 'Peak geometry must be integer metadata')
                    parts.append([int(v) for v in dataset[start:start + 8192]])
                for peak, chrom, a, b, c, d in zip(*parts):
                    invalid += not (peak not in peaks and chrom in primary_contigs
                                    and a >= 1 and b >= a and c == a - 1 and d == b)
                    peaks.add(peak); contigs[chrom] += 1
                    writer.writerow((peak, chrom, a, b, c, d))
        require(invalid == 0, 'Duplicate or invalid custom peak geometry')
    require(fingerprint(path) == before, 'Custom full-axis H5 changed during identity inspection')
    return base_axis, dict(raw_sources=sources, h5py_version=h5py.__version__, custom_ATAC_peaks=len(peaks),
        emitted_RNA_gene_axis_sha256=sha(output / 'native_rna_gene_axis.tsv'),
        emitted_custom_ATAC_geometry_sha256=sha(output / 'native_custom_merged_atac_geometry.tsv'),
        custom_ATAC_contigs=dict(contigs), custom_axis_source=str(path),
        custom_axis_is_full_features_of_1000_cell_fixture=True,
        full_source_custom_ATAC_count_availability_verified=False,
        raw_Cell_Ranger_peak_axis_not_substituted=True, counts_datasets_opened=False)


def program_join(output, contract, axis):
    ids, names = axis
    symbol_ids = defaultdict(list)
    for identity, symbol in zip(ids, names):
        symbol_ids[symbol].append(identity)
    path = Path(contract['membership']['path'])
    with path.open(newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = ('program_uid', 'cell_type', 'module', 'source_gene', 'canonical_gene',
                  'source_weight', 'original_l1_weight', 'canonical_weight_text', 'membership_sha256')
        require(set(fields) <= set(reader.fieldnames or ()), 'Frozen membership identity schema differs')
        counts, programs, members = Counter(), defaultdict(Counter), 0
        with (output / 'frozen117_member_native_gene_join.tsv').open('x', newline='') as target:
            writer = csv.writer(target, delimiter='\t')
            writer.writerow((*fields, 'native_exact_source_symbol_match', 'native_ensembl_ids'))
            for row in reader:
                hits = symbol_ids.get(row['source_gene'], [])
                status = 'unique' if len(hits) == 1 else ('ambiguous' if hits else 'missing')
                writer.writerow((*[row[name] for name in fields], status, ';'.join(hits)))
                counts[status] += 1; programs[row['program_uid']][status] += 1; members += 1
        require(len(programs) == 117, 'Frozen 117 program identity differs')
        write_tsv(output / 'frozen117_native_gene_identity_coverage.tsv',
            ('program_uid', 'member_rows', 'unique_native_symbols', 'ambiguous_native_symbols', 'missing_native_symbols'),
            ((program, sum(n.values()), n['unique'], n['ambiguous'], n['missing']) for program, n in sorted(programs.items())))
    return dict(programs=117, member_rows=members, exact_source_symbol_matches=dict(counts),
        weights_preserved_as_original_text=True, imputation=False, alias_substitution=False,
        programs_are_read_only_posthoc_projection_identities=True,
        programs_as_model_inputs_or_training_targets_authorized=False)


def fingerprint(path):
    info = Path(path).stat()
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--inspect', action='store_true', help='small contract/stat/header inspection only')
    args = parser.parse_args()
    contract = json.loads((PREP / 'source_contract.json').read_text())
    if args.inspect:
        result = dict(contract=str(PREP / 'source_contract.json'), expected_nuclei=68398,
            metadata_inputs={name: dict(path=value['path'], exists=Path(value['path']).is_file(),
                actual_bytes=Path(value['path']).stat().st_size if Path(value['path']).is_file() else None)
                for name, value in contract.items() if isinstance(value, dict) and 'path' in value},
            HDF5_opened=False, full_source_tables_parsed=False)
        print(json.dumps(result, indent=2)); return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Full source metadata parsing requires compute')
    output = REC / ('codex_gse296875_metadata_census_' + os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    (output / 'executed_source.py').write_bytes(Path(__file__).read_bytes())
    (output / 'executed_launcher.sbatch').write_bytes((HERE / LAUNCHER).read_bytes())
    (output / 'source_contract.json').write_bytes((PREP / 'source_contract.json').read_bytes())
    phase, checked, fingerprints = 'source_identity_hashes', {}, {}
    try:
        for name in contract['hash_required']:
            item = contract[name]; path = Path(item['path'])
            if not path.is_file():
                raise FileNotFoundError('Required identity/source artifact missing: ' + str(path))
            fingerprints[str(path)] = fingerprint(path)
            checked[name] = sha(path)
            require(checked[name] == item['sha256'], 'Frozen metadata/source hash differs: ' + name)
            require(fingerprint(path) == fingerprints[str(path)], 'Metadata changed during hashing')
        lock = json.loads(Path(contract['source_lock']['path']).read_text())
        require(lock['schema_version'] == 'masld-cl-gse296875-source-lock-v44'
                and lock['published_expectations']['biological_donors_passing_qc'] == 39
                and set(lock['raw_rna_sources']) == {'well' + str(i) for i in range(1, 9)}, 'Source lock semantics differ')
        phase = 'full_mapping_and_native_labels'
        identity, barcodes = census(output, contract['mapping']['path'], contract['labels']['path'])
        json_write(output / 'full_source_identity_census.json', identity)
        phase = 'native_feature_identity_and_geometry'
        axis, feature_result = axes(output, contract, lock, barcodes)
        phase = 'frozen117_read_only_native_gene_join'
        projection = program_join(output, contract, axis)
        for path, previous in fingerprints.items():
            require(fingerprint(path) == previous, 'Frozen source changed before completion: ' + path)
        pre_summary_bytes = sum(path.stat().st_size for path in output.rglob('*') if path.is_file())
        summary = dict(metadata_census_completed=True,
            biological_replication_unit='source donor; not nucleus/well', identity=identity,
            native_feature_identity=feature_result, frozen117_read_only_join=projection,
            checked_sha256=checked,
            emitted_metadata_sha256={path.name: sha(path) for path in output.iterdir()
                if path.suffix == '.tsv'},
            molecular_values_read=False, clinical_columns_read=False,
            model_fit_or_selection=False, independent_validation=False,
            full_source_ATAC_measurement_readiness_claimed=False,
            source_is_project_exposed_development=True, python=sys.version,
            saved_bytes_before_summary=pre_summary_bytes, final_inclusive_output_cap_bytes=OUTPUT_CAP,
            slurm_job_id=os.environ['SLURM_JOB_ID'])
        summary_payload = json.dumps(summary, indent=2, allow_nan=False) + '\n'
        require(pre_summary_bytes + len(summary_payload.encode()) <= OUTPUT_CAP,
                'Census output including summary exceeds 100 MiB cap')
        with (output / 'summary.json').open('x') as handle:
            handle.write(summary_payload)
        require(sum(path.stat().st_size for path in output.rglob('*') if path.is_file()) <= OUTPUT_CAP,
                'Final census output including summary exceeds 100 MiB cap')
    except BaseException as error:
        json_write(output / 'failure.json', dict(metadata_census_completed=False, phase=phase,
            exception_type=type(error).__name__, missing_artifact=isinstance(error, FileNotFoundError),
            message=str(error), traceback=traceback.format_exc(), checked_sha256=checked,
            molecular_values_read=False, clinical_columns_read=False, model_fit_or_selection=False))
        raise


if __name__ == '__main__':
    main()
