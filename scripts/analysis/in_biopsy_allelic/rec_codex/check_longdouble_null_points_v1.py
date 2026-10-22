#!/usr/bin/env python3
"""Two fixed synthetic draw11 null points, without fitting or acceptance.

LongDoubleConditional keeps its float64 row return interface unchanged. Native
float64 model.predictors and chain-rule products are retained; only aggregation
of those complete row contributions uses longdouble. A float64 aggregate of
the same contributions isolates summation effects. This is not an evaluation
in extended precision throughout and not a mathematical accuracy certificate.
"""
import argparse
from dataclasses import replace
from decimal import Decimal, localcontext
import hashlib
import importlib
import importlib.util
import io
import contextlib
import json
import os
from pathlib import Path
import platform
import resource
import sys
import time

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
POINT = REC / 'repair_p3_gamma_draw_11_22081580'
PROBE = REC / 'repair_p3_draw11_stationarity_probe_22087829'
REFERENCE = REC / 'repair_p3_draw11_decimal_gradient_22087909'
HELPER = ROOT / 'scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_longdouble_conditional_v1.py'
PINNED = {
    HELPER: '66c269cf1679ded3461564537429988925670a05ff1beb7c3b94e8625009031d',
    REFERENCE / 'executed_checker.py': '8e7e7914db125007d7fe99d81144975f17d33d9d7f98c39bb71e61a7795819fd',
    REFERENCE / 'summary.json': '220f344cda959eee142d210c49b32a1f67032e15c816079ca3fb563bf7fc4406',
    REFERENCE / 'baseline_80.json': 'efbcbbc16424e0ce747ab883cf2e398809bebc3a8aa20c8b0ffc4bc6cf7295b0',
    REFERENCE / 'baseline_110.json': '89b9d48ad180c87bff98be515dc7d2e42d5b7725dcab2cc07226b4009b79a793',
    REFERENCE / 'full_original_Newton_probe_80.json': '20cd747bf444987b19fa8c43cb1b43b257ca5401dcc85e77eb57d13d80865beb',
    REFERENCE / 'full_original_Newton_probe_110.json': '586049d217bc8161bb5bca20f0c6da288dd0cfca9e322a44ce6584ea364949da',
}
CAP = 32 * 1024**2
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def load(path):
    require(Path(path).stat().st_size <= 2 * 1024**2, 'JSON metadata size cap exceeded')
    return json.loads(Path(path).read_text())


def module_at(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def metadata():
    # This archived checker exposes a stdlib-only metadata chain, which checks
    # the original protocol, saved source and synthetic input hashes.
    pins = {str(path): digest for path, digest in PINNED.items()}
    for absolute, digest in pins.items():
        require(sha(absolute) == digest, 'Fixed helper/reference pin changed')
    archived = module_at('fixed_draw11_decimal_gradient_metadata', REFERENCE / 'executed_checker.py')
    _, inherited = archived.metadata()
    for absolute, digest in inherited.items():
        require(absolute not in pins or pins[absolute] == digest, 'Conflicting input hash')
        pins[absolute] = digest
    for absolute, digest in pins.items():
        require(sha(absolute) == digest, 'Original synthetic source/input pin changed')
    return pins


class Output:
    def __init__(self, path):
        self.path = path
        self.saved = 0
        self.started = time.monotonic()

    def write(self, name, body):
        require(time.monotonic() - self.started < 3000, '50min diagnostic time cap exceeded')
        require(self.saved + len(body) < CAP - 65536, '32MiB diagnostic output cap exceeded')
        target = self.path / name
        target.parent.mkdir(parents=True, exist_ok=True)
        with target.open('xb') as stream:
            stream.write(body)
        self.saved += len(body)

    def json(self, name, value):
        self.write(name, (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode())

    def check(self):
        actual = sum(p.stat().st_size for p in self.path.rglob('*') if p.is_file())
        require(actual == self.saved and actual < CAP - 65536, 'Saved-output accounting differs')
        require(time.monotonic() - self.started < 3000, '50min diagnostic time cap exceeded')
        return actual


def stringify(np, values):
    values = np.asarray(values)
    if values.ndim:
        return [stringify(np, row) for row in values]
    text = np.format_float_scientific(np.longdouble(values), unique=True, trim='k')
    require(np.longdouble(text) == values, 'Extended aggregate string does not roundtrip')
    return text


def dec(value):
    return Decimal(value) if isinstance(value, str) else Decimal.from_float(float(value))


def compare(np, candidate, reference, include_hessian):
    with localcontext() as context:
        context.prec = 130
        gradient = [dec(x) - dec(y) for x, y in zip(candidate['gradient'], reference['gradient'])]
        answer = {'loglik_difference': str(dec(candidate['loglik']) - dec(reference['loglik'])),
                  'gradient_differences': list(map(str, gradient)),
                  'max_abs_gradient_difference': str(max(map(abs, gradient)))}
        if include_hessian and reference.get('hessian') is not None:
            require(np.shape(reference['hessian']) == (28, 28), 'Saved Hessian shape differs')
            differences = [[dec(x) - dec(y) for x, y in zip(a, b)]
                           for a, b in zip(candidate['hessian'], reference['hessian'])]
            answer['max_abs_Hessian_difference'] = str(max(abs(v) for row in differences for v in row))
        else:
            answer['Hessian_reference_available'] = False
        return answer


def projection(np, model, p, gradient, boundary_tolerance=1e-10):
    # Diagnostic application of the existing bound/information projection;
    # no parameter, optimizer convergence flag or original failed slot changes.
    bounds = [bound if model.informative()[j] else (p[j], p[j]) for j, bound in enumerate(model.bounds())]
    scores = list(map(dec, gradient))
    active = [((lo is not None and lo == hi) or
               (lo is not None and p[j] <= lo + boundary_tolerance and scores[j] < 0) or
               (hi is not None and p[j] >= hi - boundary_tolerance and scores[j] > 0))
              for j, (lo, hi) in enumerate(bounds)]
    projected = [Decimal(0) if a else s for a, s in zip(active, scores)]
    maximum = max(map(abs, projected))
    return {'bounds': bounds, 'active_coordinates': [j for j, a in enumerate(active) if a],
            'boundary_distance_tolerance': boundary_tolerance,
            'max_projected_gradient': str(maximum),
            'fixed_point_score_below_original_1e_6': maximum < Decimal('0.000001')}


def evaluate_fixed(np, model, p, helper_class, output):
    # Explicit evaluation injection: no module globals or production evaluator
    # are patched. Contexts are distinct for every gene and every fixed point.
    contexts, scope = {}, []
    for gene in range(model.G):
        rows = [row for row in model.rows if row.gene == gene]
        entries = sum(row.n + 1 for row in rows)
        scope.append({'gene': gene, 'rows': len(rows), 'count_axis_entries': entries})
        if not (0 < len(rows) <= 32 and entries <= 100000 and
                all(0 < row.n <= 5000 for row in rows)):
            return {'suitable': False, 'reason': 'Unchanged gene row/count-workspace cap exceeded',
                    'gene_scope': scope, 'caps_weakened': False}
    # Complete cap census precedes constructor allocation and row evaluation.
    for gene in range(model.G):
        contexts[gene] = helper_class(model, p, [row for row in model.rows if row.gene == gene])
    extended_ll = np.longdouble(0)
    extended_g = np.zeros(28, dtype=np.longdouble)
    extended_h = np.zeros((28, 28), dtype=np.longdouble)
    double_ll, double_g, double_h = 0., np.zeros(28), np.zeros((28, 28))
    row_records = []
    for index, row in enumerate(model.rows):
        require(time.monotonic() - output.started < 3000, '50min diagnostic time cap exceeded')
        eta, r, jac, eta_hess = model.predictors(row, p, 0., need_hessian=True)
        value, score, hessian = contexts[row.gene].terms(row, eta, r, need_hessian=True)
        require(isinstance(value, float) and jac.dtype == eta_hess.dtype == score.dtype == hessian.dtype == np.float64,
                'Native float64 predictor/row interface unexpectedly changed')
        require(jac.shape == (2, 28) and eta_hess.shape == (28, 28) and
                score.shape == (2,) and hessian.shape == (2, 2), 'Chain-rule shapes differ')
        # Exactly Model.evaluate's nonlinear eta chain rule, in float64.
        row_g = jac.T @ score
        row_h = jac.T @ hessian @ jac + score[0] * eta_hess
        extended_ll += np.longdouble(value)
        extended_g += np.asarray(row_g, dtype=np.longdouble)
        extended_h += np.asarray(row_h, dtype=np.longdouble)
        double_ll += value
        double_g += row_g
        double_h += row_h
        row_records.append({'row_index': index, 'gene': row.gene, 'individual': row.individual,
                            'row_loglik_float64': value, 'gradient_contribution_float64': row_g.tolist()})
    require(all(np.isfinite(a).all() for a in (extended_ll, extended_g, extended_h, double_ll, double_g, double_h)),
            'Nonfinite fixed-point result')
    return {'suitable': True, 'gene_scope': scope, 'row_order_preserved': True,
            'row_records': row_records,
            'longdouble_sum_of_float64_chain_contributions': {
                'loglik': stringify(np, extended_ll), 'gradient': stringify(np, extended_g),
                'hessian': stringify(np, extended_h)},
            'float64_sum_of_same_contributions': {'loglik': double_ll, 'gradient': double_g.tolist(),
                                                  'hessian': double_h.tolist()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true', help='Source/hash checks only, without numerical imports')
    parser.add_argument('--expected-script-sha256')
    args = parser.parse_args()
    sys.dont_write_bytecode = True
    pins = metadata()
    if args.inspect:
        require('numpy' not in sys.modules and 'scipy' not in sys.modules, 'Inspect imported numerical libraries')
        print(json.dumps({'pins_verified': len(pins), 'numerical_imports': False, 'fitting': False,
                          'helper_return_and_predictors_float64_preserved': True}))
        return 0
    job = os.environ.get('SLURM_JOB_ID', '')
    require(job.isdigit() and int(job) > 0 and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute SLURM allocation required')
    require(all(os.environ.get(key) == '1' for key in THREADS), 'One numerical thread required')
    require(os.environ.get('SLURM_CPUS_PER_TASK') == '1', 'One CPU allocation required')
    require(args.expected_script_sha256 and sha(Path(__file__)) == args.expected_script_sha256,
            'Reviewed checker source hash required')
    pins[str(Path(__file__).resolve())] = args.expected_script_sha256
    directory = REC / ('check_longdouble_null_points_v1_' + job)
    directory.mkdir(mode=0o750)
    output = Output(directory)
    phase = 'archive_sources'
    try:
        # Resource allocation is 1 CPU/2GB/1h; additionally bound CPU time and
        # decoded diagnostic outputs. Existing input/source files are read-only.
        resource.setrlimit(resource.RLIMIT_CPU, (3000, 3300))
        archived = set()
        for path, digest in pins.items():
            if Path(path).suffix == '.py':
                name = 'executed_sources/' + digest[:16] + '_' + Path(path).name
                if name not in archived:
                    output.write(name, Path(path).read_bytes())
                    archived.add(name)
        output.json('source_input_pins_before.json', pins)
        phase = 'recorded_environment_and_fixture'
        import numpy as np
        import scipy
        require((platform.python_version(), np.__version__, scipy.__version__) == ('3.10.18', '2.2.6', '1.15.2'),
                'Recorded synthetic environment required')
        require(np.finfo(np.longdouble).nmant >= 63, 'Longdouble significand must have at least64bits')
        config = io.StringIO()
        with contextlib.redirect_stdout(config):
            np.show_config()
        output.json('environment.json', {'python': platform.python_version(), 'numpy': np.__version__,
            'scipy': scipy.__version__, 'machine': platform.machine(), 'Python_executable': sys.executable,
            'Python_binary_sha256': sha(Path(sys.executable).resolve()), 'NumPy_config': config.getvalue(),
            'longdouble_nmant': np.finfo(np.longdouble).nmant, 'longdouble_bits': np.finfo(np.longdouble).bits,
            'thread_environment': {key: os.environ.get(key) for key in THREADS},
            'scheduler_memory_per_node_MB': os.environ.get('SLURM_MEM_PER_NODE'),
            'scheduler_CPUs_per_task': os.environ.get('SLURM_CPUS_PER_TASK'),
            'seed': None, 'seed_reason': 'Two deterministic saved points; no sampling, fitting or proposals',
            'limits': {'CPUs': 1, 'requested_memory_GB': 2, 'requested_wall_hours': 1,
                       'CPU_soft_seconds': 3000, 'elapsed_checked_seconds': 3000, 'saved_bytes': CAP,
                       'per_gene_rows': 32, 'per_gene_count_axis_entries': 100000, 'per_row_n': 5000}})
        archive = POINT / 'executed_sources'
        sys.path.insert(0, str(archive))
        loader = importlib.import_module('bfix_v1')
        importlib.import_module('tests_v1')
        for name in ('bfix_v1', 'design_v1', 'likelihood_v1', 'tests_v1'):
            require(Path(sys.modules[name].__file__).resolve() == (archive / (name + '.py')).resolve(),
                    'Synthetic source import escaped archive')
        design, settings = loader.load()
        template = design.model()
        counts = load(POINT / 'draw_counts.json')
        keys = [[row.individual, settings['gene_order'][row.gene]] for row in template.rows]
        require(len(template.rows) == len(counts['counts']) == 173 and len(design.participants) == 60 and
                keys == counts['full_row_keys'], 'Original synthetic row count/order differs')
        require(all(isinstance(a, int) and not isinstance(a, bool) for a in counts['counts']),
                'Original counts must be exact integers')
        require(hashlib.sha256(np.asarray(counts['counts'], dtype='<i8').tobytes()).hexdigest() == counts['count_vector_sha256'],
                'Original count vector hash differs')
        rows = [replace(row, a=a) for row, a in zip(template.rows, counts['counts'])]
        model = sys.modules['likelihood_v1'].Model(rows, gene_count=8)
        require((model.G, model.Z, model.size) == (8, 7, 28) and tuple(rows) == model.rows,
                'Model coordinates or original row order differ')
        for row in rows:
            require(all(int(v) == v for v in (row.n, row.depth, row.minor, row.percent)) and
                    row.orientation in (-1, 1), 'Original integer selection/orientation differs')
            lo = max(row.minor, (row.percent * row.n + 99) // 100)
            require(row.n >= row.depth and lo <= row.a <= row.n - lo, 'Original count fails inclusion predicate')
        helper_class = module_at('fixed_longdouble_conditional', HELPER).LongDoubleConditional
        saved_points = [('baseline', load(POINT / 'p3_011.json')['fit']['null']),
                        ('full_original_Newton_probe', load(PROBE / 'probe_01.json'))]
        results, references = {}, {}
        for label, saved in saved_points:
            phase = 'fixed_point/' + label
            p = np.asarray(saved['parameters'], dtype=np.float64)
            require(p.shape == (28,) and np.isfinite(p).all(), 'Saved parameter vector differs')
            refs = {digits: load(REFERENCE / (label + '_' + str(digits) + '.json')) for digits in (80, 110)}
            for digits, reference in refs.items():
                require(reference['parameters'] == saved['parameters'] and reference['point'] == label and
                        reference['digits'] == digits and len(reference['gradient']) == 28 and
                        [[r['individual'], r['gene']] for r in reference['rows']] ==
                        [[r.individual, r.gene] for r in model.rows], 'Decimal reference point/row identity differs')
            references[label] = refs
            evaluated = evaluate_fixed(np, model, p, helper_class, output)
            if not evaluated['suitable']:
                output.json('unsuitable.json', {'point': label, **evaluated, 'no_fitting': True,
                    'source_and_input_pins': pins,
                    'max_memory_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                    'elapsed_seconds': time.monotonic() - output.started})
                output.check()
                return 0
            ll, gradient, hessian = model.evaluate(p, need_hessian=True)
            native = {'loglik': ll, 'gradient': gradient.tolist(), 'hessian': hessian.tolist()}
            record = {'point': label, 'parameters': saved['parameters'], **evaluated,
                      'original_Model_evaluate_recomputed': native,
                      'saved_native_point': {key: saved[key] for key in ('loglik', 'gradient', 'hessian') if key in saved},
                      'comparisons': {}, 'projections': {}, 'original_null_fitter_projections': {}}
            for variant in ('longdouble_sum_of_float64_chain_contributions', 'float64_sum_of_same_contributions'):
                candidate = record[variant]
                record['comparisons'][variant] = {
                    'Decimal_references': {str(digits): compare(np, candidate, reference, False)
                                           for digits, reference in refs.items()},
                    'original_recomputed': compare(np, candidate, native, True),
                    'saved_native': compare(np, candidate, saved, True)}
                record['projections'][variant] = projection(np, model, p, candidate['gradient'])
                record['original_null_fitter_projections'][variant] = projection(
                    np, model, p, candidate['gradient'], boundary_tolerance=1e-8)
            record['original_recomputed_minus_saved'] = compare(np, native, saved, True)
            output.json(label + '.json', record)
            results[label] = record
            output.check()
        phase = 'point_difference_and_final_pins'
        with localcontext() as context:
            context.prec = 130
            reference_steps = {str(digits): dec(references['full_original_Newton_probe'][digits]['loglik']) -
                               dec(references['baseline'][digits]['loglik']) for digits in (80, 110)}
            steps = {}
            for variant in ('longdouble_sum_of_float64_chain_contributions', 'float64_sum_of_same_contributions',
                            'original_Model_evaluate_recomputed', 'saved_native_point'):
                delta = dec(results['full_original_Newton_probe'][variant]['loglik']) - dec(results['baseline'][variant]['loglik'])
                steps[variant] = {'trial_minus_baseline_LL': str(delta),
                    'same_positive_sign_as_both_Decimal_references': bool(delta > 0 and all(v > 0 for v in reference_steps.values())),
                    'difference_from_Decimal_steps': {digits: str(delta - value) for digits, value in reference_steps.items()}}
        for absolute, digest in pins.items():
            require(sha(absolute) == digest, 'Source/input changed during fixed-point diagnostic')
        output.json('source_input_pins_after.json', pins)
        output.json('summary.json', {'completed': True, 'suitable': True, 'SLURM_job_id': job,
            'points': list(results), 'rows': 173, 'source_individuals': 60, 'G': 8, 'Z': 7, 'parameters': 28,
            'gene_scope': results['baseline']['gene_scope'], 'Decimal_reference_steps': {k: str(v) for k, v in reference_steps.items()},
            'point_steps': steps, 'max_memory_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'elapsed_seconds': time.monotonic() - output.started, 'saved_bytes_before_summary': output.check(),
            'saved_cap_bytes': CAP, 'sources_before_after_identical': True,
            'precision_boundaries': {'native_predictors_and_jacobian': 'float64',
                'conditional_internal_arithmetic': 'unchanged LongDoubleConditional longdouble rising products and support sums',
                'conditional_row_return': 'float64 likelihood/score/Hessian',
                'chain_rule_products': 'float64, exactly native Model.evaluate expressions',
                'primary_aggregate': 'longdouble from already-rounded float64 row contributions',
                'control_aggregate': 'float64 from the same row contributions',
                'prior_weights_and_helper_literal_constants': 'unchanged original float64 sources',
                'JSON_extended_numbers': 'roundtrip unique longdouble decimal strings'},
            'saved_baseline_Hessian_available': 'hessian' in saved_points[0][1],
            'saved_probe_Hessian_available': 'hessian' in saved_points[1][1], 'independent_Decimal_Hessian_available': False,
            'fully_extended_precision_evaluation': False, 'mathematical_accuracy_certification': False,
            'fixed_Decimal_reference_projection_boundary_tolerance': 1e-10,
            'original_null_fitter_projection_boundary_tolerance': 1e-8,
            'no_fits_or_optimizer_changes': True, 'no_proposal_or_parameter_accepted': True,
            'original_failed_slot_and_warm_cold_rules_unchanged': True,
            'calibration_or_biological_claim': False, 'restricted_or_receiving_data_read': False})
        print(json.dumps({'completed': True, 'summary_sha256': sha(directory / 'summary.json'),
                          'saved_bytes': output.check(), 'no_fit': True}), flush=True)
        return 0
    except Exception as exc:
        failure = {'completed': False, 'phase': phase, 'exception_type': type(exc).__name__,
                   'reason': str(exc) if isinstance(exc, RuntimeError) else type(exc).__name__,
                   'max_memory_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                   'partials_preserved': True, 'original_failed_slot_unchanged': True}
        # Fresh owned output only; leave completed sources/receipts untouched.
        with (directory / 'failure.json').open('x') as stream:
            json.dump(failure, stream, indent=2)
            stream.write('\n')
        print(json.dumps(failure), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
