#!/usr/bin/env python3
"""Independently reconstruct the frozen full lineage model's donor-equal losses.

Compute only. Reads full ATAC counts, identity metadata and clipped OOF profiles;
never imports the model, refits it or reads RNA, coefficients or external outcomes.
"""
import argparse
import ast
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import resource
import shutil
import sys
import time

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
MODEL = REC / 'codex_gse296875_lineage_model_22023750'
COMMON = REC / 'codex_gse296875_common_fragment_counts_22000534'
RECIPE = REC / '20261001_codex_gse296875_common_fragment_rna_atac_prespec.json'
RECIPE_SHA = '17ab84770eb3bfac46ea50d0a1400e52a70ffea1bbb7683e45b00a40e499403c'
PRODUCER_SHA = 'ccb6552d61b31f0d3fcf607f6b677c2b44909eaa6bdd6c1ae9777b2ec8e90430'
ARMS = ('lineage_mean', 'technical_only', 'donor_global_RNA', 'target_lineage_RNA')
ROWS, WINDOWS, CHUNK = 273, 306706, 8192
ATOL, RTOL = 1e-12, 1e-12  # Floating arithmetic agreement, not a scientific margin.
CAP, RESERVE, SECONDS = 4 * 1024**2, 65536, 55 * 60
PINS = {
    RECIPE: RECIPE_SHA,
    COMMON / 'summary.json': 'bd64347bc066cbd182695c466b957ddd3c6efc689b53f8db2d557aed34321f8e',
    COMMON / 'artifacts.json': '8735af6f215fc82be928f2349cb8182d3d76e834304cce8eb32e858f5880a960',
}


def require(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def file_state(path):
    require(path.is_file() and not path.is_symlink(), 'source_not_regular')
    s = path.stat()
    return {'bytes': s.st_size, 'device': s.st_dev, 'inode': s.st_ino,
            'mtime_ns': s.st_mtime_ns, 'ctime_ns': s.st_ctime_ns}


def header(path):
    with path.open('rb') as stream:
        prefix = stream.read(10)
        require(prefix[:8] == b'\x93NUMPY\x01\x00' and len(prefix) == 10, 'NPY_format_differs')
        length = int.from_bytes(prefix[8:], 'little')
        require(0 < length <= 65535, 'NPY_header_length')
        body = stream.read(length)
        require(len(body) == length, 'NPY_header_truncated')
    fields = ast.literal_eval(body.decode('latin1').strip())
    require(set(fields) == {'descr', 'fortran_order', 'shape'}, 'NPY_fields_differ')
    return {'dtype': fields['descr'], 'shape': list(fields['shape']),
            'fortran_order': fields['fortran_order'], 'offset': 10 + length,
            'header_sha256': hashlib.sha256(prefix + body).hexdigest()}


def json_read(path):
    require(path.stat().st_size <= 1024**2, 'metadata_size_cap')
    return json.loads(path.read_text())


def tsv_read(path):
    require(path.stat().st_size <= 1024**2, 'table_size_cap')
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


class Check:
    def __init__(self, directory, job):
        self.directory = directory
        self.logs = [REC / ('1593_verify_lineage_losses_' + job + '.' + ext) for ext in ('out', 'err')]
        self.start = time.monotonic()
        self.sources = {}
        self.comparisons = []
        self.phase = 'source_metadata'

    def budget(self, extra=0):
        files = list(self.directory.iterdir()) + self.logs
        require(not any(p.is_symlink() for p in files), 'output_symlink')
        used = sum(p.stat().st_size for p in files if p.is_file())
        require(used + extra <= CAP - RESERVE, 'inclusive_output_cap')
        require(time.monotonic() - self.start < SECONDS, 'verification_time_cap')
        return used

    def save_json(self, name, value):
        body = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
        self.budget(len(body))
        with (self.directory / name).open('xb') as stream:
            stream.write(body)

    def pin(self, path, digest, size=None, shape=None, dtype=None, expected_header=None):
        require(re.fullmatch('[0-9a-f]{64}', digest) is not None, 'invalid_source_SHA')
        before = file_state(path)
        require(size is None or before['bytes'] == size, 'source_bytes_differ')
        require(sha(path) == digest, 'source_SHA_differs')
        record = {'sha256': digest, 'stat': before}
        if shape is not None:
            h = header(path)
            require(h['shape'] == shape and h['dtype'] == dtype and h['fortran_order'] is False,
                    'array_header_differs')
            require(before['bytes'] == h['offset'] + math_product(shape) * np.dtype(dtype).itemsize,
                    'array_length_differs')
            require(expected_header is None or h == expected_header, 'declared_array_header_differs')
            record['header'] = h
        require(file_state(path) == before, 'source_changed_while_pinning')
        self.sources[str(path)] = record
        return record

    def artifact(self, base, name, records, shape=None, dtype=None):
        require(name in records, 'consumed_artifact_missing')
        item = records[name]
        digest = item.get('sha256', item.get('computed_sha256'))
        return self.pin(base / name, digest, item['bytes'], shape, dtype, item.get('header'))

    def compare(self, label, expected, reconstructed):
        a, b = np.asarray(expected, dtype=np.float64), np.asarray(reconstructed, dtype=np.float64)
        require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(), 'comparison_shape_or_finite')
        errors = np.abs(a - b)
        agrees = bool(np.all(errors <= ATOL + RTOL * np.abs(a)))
        self.comparisons.append({'comparison': label, 'elements': int(a.size),
            'max_absolute_difference': float(errors.max(initial=0)), 'atol': ATOL, 'rtol': RTOL,
            'all_agree': agrees})
        require(agrees, 'saved_loss_or_interval_does_not_reproduce')

    def after(self):
        self.phase = 'source_recheck'
        for name, item in self.sources.items():
            path = Path(name)
            require(file_state(path) == item['stat'] and sha(path) == item['sha256'], 'consumed_source_changed')
            if 'header' in item:
                require(header(path) == item['header'], 'consumed_header_changed')
            require(file_state(path) == item['stat'], 'source_changed_during_recheck')
            self.budget()


def math_product(values):
    result = 1
    for value in values:
        require(isinstance(value, int) and value > 0, 'array_dimension_invalid')
        result *= value
    return result


def column_values(rows, arm):
    return [float(row[arm + '_MSE']) for row in rows]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-script-sha256', required=True)
    parser.add_argument('--expected-summary-sha256', required=True)
    parser.add_argument('--expected-artifacts-sha256', required=True)
    args = parser.parse_args()
    job = os.environ.get('SLURM_JOB_ID', '')
    require(re.fullmatch('[0-9]+', job) is not None, 'compute_SLURM_job_required')
    require(all(os.environ.get(key) == '1' for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')),
            'single_thread_environment_required')
    global np
    import numpy as np
    source = Path(__file__).resolve()
    require(sha(source) == args.expected_script_sha256, 'verification_source_SHA_differs')
    directory = REC / ('codex_gse296875_verify_lineage_losses_' + job)
    directory.mkdir(mode=0o750)  # Refuse existing output; preserve every partial run.
    check = Check(directory, job)
    try:
        require(not (MODEL / 'failure.json').exists(), 'failed_producer')
        check.pin(source, args.expected_script_sha256)
        launcher = source.with_name('run_codex_gse296875_verify_lineage_losses.sbatch')
        check.pin(launcher, sha(launcher))
        for path, digest in PINS.items():
            check.pin(path, digest)
        check.pin(MODEL / 'summary.json', args.expected_summary_sha256)
        check.pin(MODEL / 'artifacts.json', args.expected_artifacts_sha256)
        summary = json_read(MODEL / 'summary.json')
        require(summary.get('development_execution_completed') is True and summary.get('complete_primary_available') is True and
                summary.get('final_scorers') == 28 and summary.get('source_inventory_before_after_identical') is True and
                summary.get('recipe_sha256') == RECIPE_SHA and summary.get('producer_sha256') == PRODUCER_SHA,
                'producer_completion_or_identity_differs')
        common_summary = json_read(COMMON / 'summary.json')
        require(common_summary.get('common_fragment_overlap_counts_completed') is True and
                common_summary.get('donors') == 39 and common_summary.get('sampled_strata') == 272 and
                common_summary.get('strata') == ROWS and common_summary.get('peaks') == WINDOWS,
                'common_count_identity_differs')
        model_artifacts = json_read(MODEL / 'artifacts.json')
        common_artifacts = json_read(COMMON / 'artifacts.json')
        tables = ('fixed_unit_identity.tsv', 'full/all273_unit_losses.tsv', 'full/donor_losses_and_effects.tsv',
                  'full/well_losses_and_effects.tsv', 'full/lineage_losses_and_effects.tsv', 'full/aggregate_summary.json')
        for name in tables:
            check.artifact(MODEL, name, model_artifacts)
        check.artifact(COMMON, 'units.tsv', common_artifacts)
        raw_name = 'ATAC_fragment_overlap_full_int32.npy'
        raw_record = check.artifact(COMMON, raw_name, common_artifacts, [ROWS, WINDOWS], '<i4')
        profiles = {}
        for arm in ARMS:
            name = 'full/' + arm + '_clipped_OOF.npy'
            profiles[arm] = check.artifact(MODEL, name, model_artifacts, [ROWS, WINDOWS], '<f8')
        draw_name = 'bootstrap_eight_well_indices.npy'
        check.artifact(MODEL, draw_name, model_artifacts, [10000, 8], '<i2')
        units = tsv_read(MODEL / 'fixed_unit_identity.tsv')
        measured = tsv_read(COMMON / 'units.tsv')
        recipe = json_read(RECIPE)
        lineages = recipe['source_identity']['ordered_author_lineages']
        roster = recipe['source_identity']['well_donor_roster']
        wells = tuple('well' + str(i) for i in range(1, 9))
        require(set(roster) == set(wells) and len(lineages) == 7 and len(set(lineages)) == 7 and
                len(units) == ROWS and len(measured) == ROWS, 'roster_dimensions_differ')
        donor_order = [d for well in wells for d in roster[well]]
        require(len(donor_order) == 39 and len(set(donor_order)) == 39, 'biological_donor_roster_differs')
        eligible = []
        keys = ('unit_index', 'donor_id', 'well_id', 'native_source_celltype', 'nuclei', 'sampling_state')
        for i, (unit, counted) in enumerate(zip(units, measured)):
            require(all(unit[key] == counted[key] for key in keys) and int(unit['unit_index']) == i and
                    unit['donor_id'] == donor_order[i // 7] and unit['native_source_celltype'] == lineages[i % 7] and
                    unit['donor_id'] in roster[unit['well_id']], 'RNA_ATAC_metadata_join_differs')
            sampled = int(unit['nuclei']) > 0
            require(sampled == (unit['sampling_state'] == 'observed'), 'sampling_state_differs')
            if sampled:
                require(0 < int(counted['full_autosomal_fragment_records']) < 2**53, 'target_denominator_undefined')
                eligible.append(i)
            else:
                require(unit['donor_id'] == '733' and unit['native_source_celltype'] == 'B cells' and
                        int(counted['full_autosomal_fragment_records']) == 0, 'unsampled_unit_differs')
        require(len(eligible) == 272, 'sampled_primary_roster_differs')
        aggregate = json_read(MODEL / 'full/aggregate_summary.json')
        require(aggregate == summary['results']['full'] and aggregate.get('complete_fixed_metadata_roster') is True and
                aggregate.get('sampled_units') == 272 and aggregate.get('biological_donors') == 39 and
                aggregate.get('well_clusters') == 8 and aggregate.get('bootstrap_draws') == 10000 and
                aggregate.get('bootstrap_seed') == 20261001 and aggregate.get('bootstrap_reweights_resampled_donor_count') is True,
                'primary_result_identity_differs')
        check.save_json('consumed_sources_before.json', check.sources)
        for p in (source, launcher):
            check.budget(p.stat().st_size)
            shutil.copyfile(p, directory / ('executed_' + p.name))
        check.phase = 'full_window_loss_reconstruction'
        losses = np.full((ROWS, 4), np.nan)
        streams = {arm: (MODEL / ('full/' + arm + '_clipped_OOF.npy')).open('rb') for arm in ARMS}
        try:
            with (COMMON / raw_name).open('rb') as counts_stream:
                for i in range(ROWS):
                    totals = np.zeros(4)
                    sampled = i in eligible
                    for begin in range(0, WINDOWS, CHUNK):
                        length = min(CHUNK, WINDOWS - begin)
                        if sampled:
                            counts_stream.seek(raw_record['header']['offset'] + (i * WINDOWS + begin) * 4)
                            count = np.fromfile(counts_stream, dtype='<i4', count=length)
                            require(count.size == length and np.all(count >= 0), 'ATAC_chunk_invalid')
                            # Independently ordered arithmetic: count/D, then per million, then log1p.
                            target = np.log1p((count.astype(np.float64) / int(measured[i]['full_autosomal_fragment_records'])) * 1e6)
                        for a, arm in enumerate(ARMS):
                            stream = streams[arm]
                            stream.seek(profiles[arm]['header']['offset'] + (i * WINDOWS + begin) * 8)
                            prediction = np.fromfile(stream, dtype='<f8', count=length)
                            require(prediction.size == length, 'prediction_chunk_truncated')
                            if sampled:
                                require(np.isfinite(prediction).all() and np.all(prediction >= 0), 'clipped_prediction_invalid')
                                error = prediction - target
                                totals[a] += float(np.dot(error, error))
                            else:
                                require(np.isnan(prediction).all(), 'unsampled_prediction_is_not_NaN')
                    if sampled:
                        losses[i] = totals / WINDOWS
                    check.budget()
        finally:
            for stream in streams.values():
                stream.close()
        saved_units = tsv_read(MODEL / 'full/all273_unit_losses.tsv')
        require(len(saved_units) == ROWS, 'saved_unit_roster_differs')
        for i, row in enumerate(saved_units):
            require(int(row['unit_index']) == i and row['donor_id'] == units[i]['donor_id'] and
                    row['well_id'] == units[i]['well_id'] and row['lineage'] == units[i]['native_source_celltype'] and
                    row['metadata_eligible'] == str(i in eligible), 'saved_unit_identity_differs')
            if i not in eligible:
                require(all(row[arm + '_MSE'] == '' for arm in ARMS), 'unsampled_loss_is_not_missing')
        for a, arm in enumerate(ARMS):
            check.compare('all272_unit_' + arm, [float(saved_units[i][arm + '_MSE']) for i in eligible], losses[eligible, a])
        check.phase = 'biological_donor_and_well_aggregation'
        donor_losses = []
        donor_wells = []
        donor_rows = tsv_read(MODEL / 'full/donor_losses_and_effects.tsv')
        require([row['donor_id'] for row in donor_rows] == donor_order, 'saved_donor_roster_differs')
        for donor, row in zip(donor_order, donor_rows):
            indices = [i for i in eligible if units[i]['donor_id'] == donor]
            well = units[indices[0]]['well_id']
            require(row['mode'] == 'full' and row['well_id'] == well and int(row['sampled_lineages']) == len(indices) and
                    len(indices) == (6 if donor == '733' else 7), 'saved_donor_replication_differs')
            donor_losses.append(np.mean(losses[indices], axis=0))
            donor_wells.append(well)
        donor_losses = np.asarray(donor_losses)
        groups = [('donor', donor_rows, donor_losses)]
        for a, arm in enumerate(ARMS):
            check.compare('all39_donor_' + arm, column_values(donor_rows, arm), donor_losses[:, a])
            check.compare('primary_equal39_' + arm, aggregate['arm_losses'][arm], np.mean(donor_losses[:, a]))
        for kind, key, order in (('well', 'well_id', wells), ('lineage', 'lineage', tuple(lineages))):
            rows = tsv_read(MODEL / ('full/' + kind + '_losses_and_effects.tsv'))
            require(tuple(row[key] for row in rows) == order, 'saved_group_roster_differs')
            group_losses = []
            for value, row in zip(order, rows):
                if kind == 'well':
                    indices = [j for j, well in enumerate(donor_wells) if well == value]
                    group_losses.append(np.mean(donor_losses[indices], axis=0))
                else:
                    indices = [i for i in eligible if units[i]['native_source_celltype'] == value]
                    group_losses.append(np.mean(losses[indices], axis=0))
                require(row['mode'] == 'full' and int(row['donor_n']) == len(indices), 'saved_group_donor_n_differs')
            group_losses = np.asarray(group_losses)
            for a, arm in enumerate(ARMS):
                check.compare('all_' + kind + '_' + arm, column_values(rows, arm), group_losses[:, a])
            groups.append((kind, rows, group_losses))
        for kind, rows, values in groups:
            values = donor_losses if kind == 'donor' else values
            for a, arm in enumerate(ARMS[:3]):
                check.compare('all_' + kind + '_lineage_minus_' + arm,
                              [float(row['lineage_minus_' + arm]) for row in rows], values[:, 3] - values[:, a])
        check.phase = 'expanded_repeated_donor_bootstrap'
        draws = np.load(MODEL / draw_name, allow_pickle=False)
        regenerated = np.random.Generator(np.random.PCG64(20261001)).integers(0, 8, size=(10000, 8), dtype=np.int16)
        require(np.array_equal(draws, regenerated), 'bootstrap_draws_do_not_reproduce_exactly')
        well_donors = [[j for j, value in enumerate(donor_wells) if value == well] for well in wells]
        require([len(js) for js in well_donors] == [len(roster[well]) for well in wells] and
                sorted(len(js) for js in well_donors) == [4, 5, 5, 5, 5, 5, 5, 5],
                'unequal_well_donor_counts_differ')
        boot = np.empty((10000, 4))
        for i, selected in enumerate(draws):
            # Expand each sampled well into all of its donors, retaining repetitions.
            repeated = [j for selected_well in selected for j in well_donors[int(selected_well)]]
            boot[i] = np.mean(donor_losses[repeated], axis=0)
        points = np.mean(donor_losses, axis=0)
        contrasts = {}
        for a, arm in enumerate(ARMS):
            check.compare('conditional_interval_' + arm, aggregate['arm_conditional_intervals'][arm],
                          np.percentile(boot[:, a], [2.5, 97.5], method='linear'))
            if a < 3:
                vector = boot[:, 3] - boot[:, a]
                saved = aggregate['contrasts'][arm]
                check.compare('primary_lineage_minus_' + arm, saved['lineage_minus_comparator'], points[3] - points[a])
                interval = np.percentile(vector, [2.5, 97.5], method='linear')
                check.compare('contrast_interval_lineage_minus_' + arm, saved['conditional_95_percentile_interval'], interval)
                check.budget(vector.nbytes + 4096)
                with (directory / ('reconstructed_bootstrap_lineage_minus_' + arm + '.npy')).open('xb') as stream:
                    np.save(stream, vector, allow_pickle=False)
                contrasts[arm] = {'lineage_minus_comparator': float(points[3] - points[a]),
                    'conditional_95_percentile_interval': interval.tolist()}
        require(aggregate['all_three_necessary_point_directions_negative'] == bool(np.all(points[3] - points[:3] < 0)),
                'necessary_direction_summary_differs')
        check.after()
        require(not (MODEL / 'failure.json').exists(), 'producer_failed_during_verification')
        check.save_json('consumed_sources_after.json', check.sources)
        check.save_json('numerical_comparisons.json', check.comparisons)
        check.save_json('summary.json', {'verification_completed': True, 'all_comparisons_agree': True,
            'biological_donors': 39, 'well_clusters': 8, 'sampled_units': 272, 'counted_windows': WINDOWS,
            'unsampled_733_B_cell_predictions_NaN': True, 'donor_equal_after_available_lineage_average': True,
            'bootstrap_expands_repeated_wells_into_repeated_donors': True, 'bootstrap_draws_exactly_reproduced': True,
            'contrasts': contrasts, 'atol': ATOL, 'rtol': RTOL, 'producer_sha256': PRODUCER_SHA,
            'producer_summary_sha256': args.expected_summary_sha256, 'producer_artifacts_sha256': args.expected_artifacts_sha256,
            'source_hashes_headers_and_stat_unchanged': True, 'source_model_code_imported_or_reused': False,
            'RNA_or_coefficients_or_raw_OOF_or_half_targets_or_external_outcomes_read': False,
            'model_refitted_or_frozen_recipe_changed': False, 'new_independent_validation_established': False,
            'conditional_fixed_loss_intervals_not_refit_or_population_calibrated': True,
            'Python': sys.version, 'NumPy': np.__version__, 'NumPy_path': np.__file__,
            'thread_environment': {key: os.environ.get(key) for key in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS')},
            'ru_maxrss_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'saved_output_cap_bytes_inclusive': CAP, 'elapsed_seconds': time.monotonic() - check.start})
        print(json.dumps({'verification_completed': True, 'summary_sha256': sha(directory / 'summary.json'),
                          'saved_bytes': check.budget()}), flush=True)
        check.budget()
    except Exception as exc:
        reason = str(exc) if isinstance(exc, RuntimeError) else type(exc).__name__
        failure = {'verification_completed': False, 'phase': check.phase, 'reason': reason,
                   'partial_outputs_preserved': True, 'no_model_or_recipe_change': True}
        body = (json.dumps(failure, sort_keys=True) + '\n').encode()
        with (directory / 'failure.json').open('xb') as stream:
            stream.write(body)
        print(json.dumps(failure), file=sys.stderr)
        return 1
    return 0


if __name__ == '__main__':
    sys.exit(main())
