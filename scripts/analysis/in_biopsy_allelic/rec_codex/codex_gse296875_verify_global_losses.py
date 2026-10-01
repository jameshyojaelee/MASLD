#!/usr/bin/env python3
"""Independent compute-only reconstruction of frozen source-global losses and CIs."""
import argparse
import ast
import contextlib
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import resource
import shutil
import sys
import time

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
HERE = ROOT / 'scripts/analysis/in_biopsy_allelic/rec_codex'
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
MODEL = REC / 'codex_gse296875_donor_global_model_22001071'
RAW = REC / 'codex_gse296875_donor_global_raw_22000941'
ARMS = ('donor_mean', 'technical_only', 'donor_global_RNA_PC', 'technical_source_fractions')
DIRECTIONS = (('full_RNA_to_full_ATAC', 'full'),
              ('halfA_RNA_to_halfB_ATAC', 'halfB'),
              ('halfB_RNA_to_halfA_ATAC', 'halfA'))
CONTRASTS = (('RNA_PC_minus_technical', 2, 1),
             ('RNA_PC_minus_mean', 2, 0),
             ('RNA_PC_minus_source_fractions', 2, 3))
WELLS = tuple('well' + str(i) for i in range(1, 9))
N, R, CHUNK = 39, 306706, 8192
DRAW_N, SEED = 10000, 20261001
ATOL, RTOL = 1e-12, 1e-12  # Floating arithmetic agreement only; no statistical margin.
CAP, RESERVE, SECONDS = 4 * 1024**2, 65536, 55 * 60
PINS = {
    MODEL / 'summary.json': 'd47a86e5948a1a5e29ac9aaf74dbb733133d751a3de37f6bc29ea2b6b66523d4',
    MODEL / 'artifacts.json': '38a85a0d922e8e65c5f2369ad2d9a7936cf4483b20d38037d00c2aebef6e4948',
    MODEL / 'manifest.json': 'ce2d1b06c5014f18864e39057799d2028dd2df12906eb5841281e841febe78fd',
    RAW / 'summary.json': '9684d83681b898ac507680e2df39a4781b42daccbd2abfac6342faa88031fb19',
    RAW / 'artifacts.json': '5629741635252c1869cc68da59c365833ad8912fa7f46d62fef25335beead111',
}
FLAGS = {'model_code_imported_or_reused': False, 'model_refitted': False,
    'RNA_SVD_or_coefficients_read': False, 'raw_OOF_profiles_read': False,
    'external_or_H3_or_clinical_or_genotype_or_allelic_values_read': False,
    'independent_validation_or_practical_utility_established': False,
    'bootstrap_is_conditional_descriptive_not_population_calibrated': True}


class VerificationFailure(RuntimeError):
    pass


def require(ok, reason):
    if not ok:
        raise VerificationFailure(reason)


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            h.update(block)
    return h.hexdigest()


def fingerprint(path):
    require(path.is_file() and not path.is_symlink(), 'source_not_regular')
    s = path.stat()
    return {'bytes': s.st_size, 'device': s.st_dev, 'inode': s.st_ino,
            'mtime_ns': s.st_mtime_ns, 'ctime_ns': s.st_ctime_ns}


def read_json(path):
    require(path.stat().st_size <= 1024**2, 'metadata_size_cap')
    return json.loads(path.read_text())


class Run:
    def __init__(self, out, job):
        self.out = out
        self.logs = [REC / f'1594_verify_global_losses_{job}.{s}' for s in ('out', 'err')]
        self.started = time.monotonic()
        self.phase = 'source_guards'
        self.comparisons = []

    def saved(self):
        paths = list(self.out.iterdir()) + self.logs
        require(not any(p.is_symlink() for p in paths), 'output_symlink')
        return sum(p.stat().st_size for p in paths if p.is_file())

    def check(self, extra=0):
        require(time.monotonic() - self.started < SECONDS, 'time_cap')
        require(self.saved() + extra <= CAP - RESERVE, 'inclusive_output_cap')

    def json(self, name, value):
        raw = (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode()
        self.check(len(raw))
        with (self.out / name).open('xb') as stream:
            stream.write(raw)

    def progress(self, phase):
        self.phase = phase
        self.check()
        print(json.dumps({'phase': phase, 'elapsed_seconds': time.monotonic() - self.started}), flush=True)

    def compare(self, label, reconstructed, persisted):
        a, b = np.asarray(reconstructed, dtype=np.float64), np.asarray(persisted, dtype=np.float64)
        require(a.shape == b.shape and np.isfinite(a).all() and np.isfinite(b).all(), 'comparison_shape_or_nonfinite')
        error = np.abs(a - b)
        allowed = ATOL + RTOL * np.abs(b)
        scale = np.maximum(np.abs(b), np.finfo(np.float64).tiny)
        row = {'comparison': label, 'shape': list(a.shape), 'elements': int(a.size),
            'max_absolute_difference': float(error.max(initial=0)),
            'max_relative_difference': float((error / scale).max(initial=0)),
            'atol': ATOL, 'rtol': RTOL, 'all_agree': bool(np.all(error <= allowed))}
        self.comparisons.append(row)
        return row


def npy_header(path):
    # Independent reader of the NPY 1.0 format, with no array-body read here.
    with path.open('rb') as stream:
        magic = stream.read(8)
        require(magic == b'\x93NUMPY\x01\x00', 'NPY_version_or_magic')
        length_bytes = stream.read(2)
        require(len(length_bytes) == 2, 'NPY_short_length')
        length = int.from_bytes(length_bytes, 'little')
        require(0 < length <= 65535, 'NPY_header_length')
        raw = stream.read(length)
        require(len(raw) == length, 'NPY_short_header')
        value = ast.literal_eval(raw.decode('latin1').strip())
        require(set(value) == {'descr', 'shape', 'fortran_order'}, 'NPY_header_fields')
    return {'dtype': value['descr'], 'shape': list(value['shape']),
            'fortran_order': value['fortran_order'], 'offset': 10 + length,
            'header_sha256': hashlib.sha256(magic + length_bytes + raw).hexdigest()}


def pin_source(path, descriptor, shape=None, dtype=None):
    before = fingerprint(path)
    require(before['bytes'] == descriptor['bytes'] and sha(path) == descriptor['sha256'], 'source_bytes_or_SHA_changed')
    row = {'sha256': descriptor['sha256'], 'stat': before}
    if shape is not None:
        header = npy_header(path)
        require(header == descriptor['header'] and header['shape'] == list(shape)
                and header['dtype'] == dtype and header['fortran_order'] is False, 'NPY_shape_dtype_header_changed')
        require(header['offset'] + math.prod(shape) * 8 == before['bytes'], 'NPY_payload_length')
        row['header'] = header
    require(fingerprint(path) == before, 'source_changed_during_initial_hash')
    return row


@contextlib.contextmanager
def readonly_map(path, shape, dtype):
    value = np.load(path, mmap_mode='r', allow_pickle=False)
    try:
        require(value.shape == shape and value.dtype.str == dtype and value.flags.c_contiguous
                and not value.flags.writeable, 'readonly_array_contract')
        yield value
    finally:
        value._mmap.close()


def table(path, expected_fields):
    with path.open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        require(reader.fieldnames == list(expected_fields), 'loss_table_header')
        rows = list(reader)
    require(all(None not in r for r in rows), 'loss_table_row_schema')
    return rows


def source_donors(manifest):
    # Read identity and the three ATAC D columns only; other donor metadata are unused strings.
    fields = ['donor_index', 'donor_id', 'well_id'] + [p + '_autosomal_fragment_records' for _, p in DIRECTIONS]
    with (RAW / 'donors.tsv').open(newline='') as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        require(set(fields) <= set(reader.fieldnames), 'donor_metadata_header')
        donors = [{k: r[k] for k in fields} for r in reader]
    require(len(donors) == N and [r['donor_index'] for r in donors] == list(map(str, range(N)))
            and len({r['donor_id'] for r in donors}) == N, 'donor_identity_axis')
    require([r['donor_id'] for r in donors] == manifest['ordered_donor_ids'], 'manifest_donor_order')
    groups = [np.array([i for i, r in enumerate(donors) if r['well_id'] == w], dtype=np.int64) for w in WELLS]
    require([len(g) for g in groups] == [5, 5, 4, 5, 5, 5, 5, 5], 'fixed_well_sizes')
    require(sorted(int(i) for g in groups for i in g) == list(range(N)), 'well_partition')
    require([f['well_id'] for f in manifest['folds']] == list(WELLS), 'manifest_well_order')
    for group, fold in zip(groups, manifest['folds']):
        require(group.tolist() == fold['held_indices'] and
                sorted(fold['train_indices']) == sorted(set(range(N)) - set(group.tolist())), 'manifest_held_partition')
    denoms = {}
    for _, partition in DIRECTIONS:
        raw = [d[partition + '_autosomal_fragment_records'] for d in donors]
        require(all(v.isdecimal() and 0 < int(v) <= 2**64 - 1 for v in raw), 'positive_integer_autosomal_D')
        denoms[partition] = np.array([int(v) for v in raw], dtype=np.uint64)
    return donors, groups, denoms


def reconstruct_losses(direction, partition, denoms, run):
    result = np.empty((N, len(ARMS)), dtype=np.float64)
    for arm_index, arm in enumerate(ARMS):
        run.progress('losses/' + direction + '/' + arm)
        contributions = [[] for _ in range(N)]
        # Reopen/close both maps for each head; never retain all OOF profiles or a full y matrix.
        with readonly_map(RAW / f'ATAC_donor_{partition}_counts_uint64.npy', (N, R), '<u8') as c, \
                readonly_map(MODEL / direction / (arm + '_OOF_clipped.npy'), (N, R), '<f8') as prediction:
            for start in range(0, R, CHUNK):
                stop = min(start + CHUNK, R)
                count_chunk = c[:, start:stop]
                require(np.all(count_chunk <= denoms[:, None]), 'interval_unique_count_exceeds_autosomal_total')
                target = np.array(count_chunk, dtype=np.float64)
                # Different evaluation order from the producer, same declared log1p FPM estimand.
                target *= 1e6
                target /= denoms.astype(np.float64)[:, None]
                np.log1p(target, out=target)
                clipped = prediction[:, start:stop]
                require(np.isfinite(target).all() and np.isfinite(clipped).all() and np.all(clipped >= 0), 'target_or_clipped_profile_domain')
                error = target - clipped
                np.square(error, out=error)
                sums = error.sum(axis=1, dtype=np.float64)
                for i, value in enumerate(sums):
                    contributions[i].append(float(value))
                del count_chunk, clipped, target, error, sums
                run.check()
        result[:, arm_index] = [math.fsum(parts) / R for parts in contributions]
    return result


def verify_tables(direction, losses, donors, groups, run):
    column_names = [a + '_MSE' for a in ARMS] + [name for name, _, _ in CONTRASTS]
    all_values = np.column_stack([losses] + [(losses[:, a] - losses[:, b])[:, None] for _, a, b in CONTRASTS])
    rows = table(MODEL / direction / 'donor_losses.tsv', ['donor_index', 'donor_id', 'well_id'] + column_names)
    require(len(rows) == N and all(all(row[k] == donor[k] for k in ('donor_index', 'donor_id', 'well_id'))
                for row, donor in zip(rows, donors)), 'persisted_loss_donor_identity')
    run.compare(direction + '/all_donor_arm_losses_and_contrasts', all_values,
                [[float(r[k]) for k in column_names] for r in rows])
    wells = table(MODEL / direction / 'well_losses.tsv', ['well_id', 'donors'] + column_names)
    require(len(wells) == 8 and [r['well_id'] for r in wells] == list(WELLS)
            and [int(r['donors']) for r in wells] == [len(g) for g in groups], 'persisted_well_loss_identity')
    well_values = np.array([[math.fsum(all_values[group, j].tolist()) / len(group)
                             for j in range(all_values.shape[1])] for group in groups])
    run.compare(direction + '/all_well_arm_losses_and_contrasts', well_values,
                [[float(r[k]) for k in column_names] for r in wells])
    return column_names, all_values, well_values


def verify_bootstrap(direction, losses, groups, draws, metrics, run):
    run.progress('donor_weighted_whole_well_bootstrap/' + direction)
    bootstrap = np.empty((DRAW_N, len(ARMS)), dtype=np.float64)
    donor_counts = np.empty(DRAW_N, dtype=np.int64)
    # Independently expand well multisets to repeated donor observations.
    # This does not reuse the producer's indexed sums of well totals.
    for i, draw in enumerate(draws):
        indices = [int(donor) for well in draw for donor in groups[int(well)]]
        donor_counts[i] = len(indices)
        for j in range(len(ARMS)):
            bootstrap[i, j] = math.fsum(float(losses[k, j]) for k in indices) / len(indices)
        if i % 1024 == 0:
            run.check()
    intervals = {}
    for name, a, b in CONTRASTS:
        values = bootstrap[:, a] - bootstrap[:, b]
        with readonly_map(MODEL / direction / (name + '_bootstrap.npy'), (DRAW_N,), '<f8') as saved:
            run.compare(direction + '/' + name + '/all_10000_bootstrap_contrasts', values, saved)
        ordered = np.sort(values)
        # Independent linear percentile calculation, not np.percentile.
        endpoints = []
        for q in (0.025, 0.975):
            position = (DRAW_N - 1) * q
            lo = math.floor(position)
            fraction = position - lo
            endpoints.append(float(ordered[lo] + fraction * (ordered[min(lo + 1, DRAW_N - 1)] - ordered[lo])))
        run.compare(direction + '/' + name + '/linear_95_percentiles', endpoints,
                    metrics['contrasts'][name]['descriptive_conditional_95_percentile_interval'])
        intervals[name] = endpoints
    return intervals, {'draws': DRAW_N, 'seed': SEED, 'bit_generator': 'PCG64',
        'resampled_donor_count_min': int(donor_counts.min()), 'resampled_donor_count_max': int(donor_counts.max()),
        'unequal_well_sizes_preserved': True, 'repeated_donors_from_repeated_wells_retained': True,
        'conditional_on_saved_overlapping_fits': True, 'calibrated_population_interval': False}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-script-sha256', required=True)
    args = parser.parse_args()
    job = os.environ.get('SLURM_JOB_ID', '')
    require(job.isdecimal() and int(job) > 0 and not os.environ.get('SLURM_ARRAY_TASK_ID'), 'individual_compute_job_required')
    require(sha(Path(__file__).resolve()) == args.expected_script_sha256, 'reviewed_verifier_SHA_changed')
    global np
    import numpy as np
    out = REC / ('codex_gse296875_verify_global_losses_' + job)
    owned, run, inventory, results = False, None, {}, {}
    try:
        out.mkdir(mode=0o700, exist_ok=False)
        owned, run = True, Run(out, job)
        for path, expected in PINS.items():
            inventory[str(path)] = pin_source(path, {'bytes': path.stat().st_size, 'sha256': expected})
        summary = read_json(MODEL / 'summary.json')
        require(summary['source_donor_global_model_completed'] and summary['source_development_only']
                and summary['donors'] == N and summary['wells'] == 8 and summary['windows'] == R
                and summary['complete_folds'] == 24 and not summary['external_evaluation_run'], 'model_not_complete_fixed_source_fit')
        raw_summary = read_json(RAW / 'summary.json')
        require(raw_summary['raw_donor_global_aggregation_completed'] and raw_summary['donors'] == N
                and raw_summary['ATAC_peaks'] == R, 'raw_aggregate_not_complete')
        model_artifacts, raw_artifacts = read_json(MODEL / 'artifacts.json'), read_json(RAW / 'artifacts.json')
        specifications = [(MODEL, 'donors.tsv', None, None), (RAW, 'donors.tsv', None, None),
                          (MODEL, 'bootstrap_well_indices.npy', (DRAW_N, 8), '<i8')]
        for direction, partition in DIRECTIONS:
            specifications.append((RAW, f'ATAC_donor_{partition}_counts_uint64.npy', (N, R), '<u8'))
            specifications += [(MODEL, direction + '/' + a + '_OOF_clipped.npy', (N, R), '<f8') for a in ARMS]
            specifications += [(MODEL, direction + '/' + name, None, None)
                               for name in ('donor_losses.tsv', 'well_losses.tsv', 'metrics.json')]
            specifications += [(MODEL, direction + '/' + c + '_bootstrap.npy', (DRAW_N,), '<f8') for c, _, _ in CONTRASTS]
        run.progress('all_consumed_array_SHA_stat_and_header_guards')
        for base, name, shape, dtype in specifications:
            path = base / name
            descriptor = (model_artifacts if base == MODEL else raw_artifacts)[name]
            inventory[str(path)] = pin_source(path, descriptor, shape, dtype)
            run.check()
        require(inventory[str(MODEL / 'donors.tsv')]['sha256'] == inventory[str(RAW / 'donors.tsv')]['sha256'], 'raw_model_donor_metadata_changed')
        source = Path(__file__).resolve()
        launcher = HERE / 'run_codex_gse296875_verify_global_losses.sbatch'
        for name, path in [('executed_source.py', source), ('executed_launcher.sbatch', launcher)]:
            run.check(path.stat().st_size)
            shutil.copyfile(path, out / name)
            inventory[str(path)] = pin_source(path, {'bytes': path.stat().st_size, 'sha256': sha(path)})
        run.json('protocol.json', {'source_inventory_before': inventory, 'source_model_job': '22001071',
            'source_raw_job': '22000941', 'target': 'ln(1 + 1e6 * unique_interval_fragment_count / unique_autosomal_fragment_total)',
            'units': 'natural log1p fragments per million autosomal unique fragment records',
            'donors': N, 'windows': R, 'chunk_windows': CHUNK, 'arms': ARMS,
            'directions_and_target_partitions': DIRECTIONS, 'bootstrap_seed': SEED,
            'bootstrap_well_draws': DRAW_N, 'numerical_atol': ATOL, 'numerical_rtol': RTOL,
            'statistical_margin_added': False, 'inclusive_saved_cap': CAP, 'internal_seconds_cap': SECONDS, **FLAGS})
        manifest = read_json(MODEL / 'manifest.json')
        donors, groups, denoms = source_donors(manifest)
        with readonly_map(MODEL / 'bootstrap_well_indices.npy', (DRAW_N, 8), '<i8') as saved:
            draws = np.random.Generator(np.random.PCG64(SEED)).integers(0, 8, size=(DRAW_N, 8))
            require(np.array_equal(draws, saved), 'exact_PCG64_well_draws_not_reproduced')
        for direction, partition in DIRECTIONS:
            losses = reconstruct_losses(direction, partition, denoms[partition], run)
            names, donor_values, well_values = verify_tables(direction, losses, donors, groups, run)
            metrics = read_json(MODEL / direction / 'metrics.json')
            embedded = summary['primary'] if partition == 'full' else summary['cross_half_diagnostics'][direction]
            require(metrics == embedded, 'summary_direction_metrics_differ_from_saved_metrics')
            require(metrics['donors'] == N and metrics['regions'] == R and metrics['wells'] == 8
                    and metrics['equal_donor_MSE'] and metrics['bootstrap']['draws'] == DRAW_N
                    and metrics['bootstrap']['seed'] == SEED and metrics['bootstrap']['bit_generator'] == 'PCG64'
                    and metrics['bootstrap']['percentile_method'] == 'linear', 'persisted_metric_contract')
            means = [math.fsum(losses[:, j].tolist()) / N for j in range(len(ARMS))]
            run.compare(direction + '/overall_donor_equal_MSE', means, [metrics['absolute_MSE'][a] for a in ARMS])
            for name, a, b in CONTRASTS:
                contrast = math.fsum((losses[:, a] - losses[:, b]).tolist()) / N
                run.compare(direction + '/' + name + '/overall_donor_equal_contrast', contrast,
                            metrics['contrasts'][name]['donor_equal_mean'])
            guard = means[2] < means[1] and means[2] < means[0]
            require(guard == metrics['necessary_descriptive_source_accuracy_guards_passed'], 'persisted_source_guard_differs')
            intervals, bootstrap_receipt = verify_bootstrap(direction, losses, groups, draws, metrics, run)
            results[direction] = {'target_partition': partition, 'column_order': names,
                'ordered_donor_identity': [{k: d[k] for k in ('donor_index', 'donor_id', 'well_id')} for d in donors],
                'donor_MSE_and_contrasts': donor_values.tolist(), 'well_order': WELLS,
                'well_donor_counts': [len(g) for g in groups], 'well_MSE_and_contrasts': well_values.tolist(),
                'overall_donor_equal_MSE': dict(zip(ARMS, means)), 'linear_descriptive_95_percentiles': intervals,
                'source_accuracy_guards_reproduced': True, 'source_accuracy_guards_passed': guard,
                'bootstrap': bootstrap_receipt, 'clipped_profiles_all_finite_nonnegative': True,
                'raw_to_clipped_transform_identity_not_rechecked': True}
            run.json(direction + '_verification.json', results[direction])
        run.progress('all_consumed_sources_final_SHA_stat_and_headers')
        for name, original in inventory.items():
            path = Path(name)
            require(fingerprint(path) == original['stat'] and sha(path) == original['sha256'], 'source_changed_before_final_hash')
            if 'header' in original:
                require(npy_header(path) == original['header'], 'source_header_changed')
            require(fingerprint(path) == original['stat'], 'source_changed_during_final_hash')
            run.check()
        passed = all(c['all_agree'] for c in run.comparisons)
        run.json('summary.json', {'verification_completed': True, 'all_numerical_comparisons_agree': passed,
            'SLURM_job_id': job, 'source_hash_stat_headers_before_after_identical': True,
            'source_inventory': inventory, 'comparisons': run.comparisons, 'results': results,
            'exact_10000_PCG64_well_multisets_reproduced': True, 'donor_count_weighting_verified': True,
            'biological_donors': N, 'well_clusters': 8, 'interval_limits': 'Eight clusters; fixed overlapping fits, no unconditional training uncertainty or calibrated population coverage.',
            'Python_version': sys.version, 'NumPy_version': np.__version__,
            'max_RSS_KiB_linux': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'elapsed_seconds': time.monotonic() - run.started, 'saved_bytes_before_summary': run.saved(), **FLAGS})
        run.check()
        print(json.dumps({'verification_completed': True, 'all_numerical_comparisons_agree': passed,
            'summary': str(out / 'summary.json'), 'summary_SHA256': sha(out / 'summary.json'),
            'inclusive_saved_bytes': run.saved(), **FLAGS}), flush=True)
        run.check()
        return 0 if passed else 1
    except Exception as exc:
        failure = {'verification_completed': False, 'phase': run.phase if run else 'fresh_output',
            'categorical_reason': str(exc) if isinstance(exc, VerificationFailure) else type(exc).__name__,
            'partial_outputs_preserved': True, 'completed_directions': list(results), **FLAGS}
        if owned:
            with (out / 'failure.json').open('x') as stream:
                json.dump(failure, stream, indent=2, allow_nan=False)
                stream.write('\n')
        print(json.dumps(failure), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
