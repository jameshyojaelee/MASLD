#!/usr/bin/env python3
"""Separate synthetic draw11 null-only pilot; only evaluate arithmetic changes.

The original two-start fit, cold start, optimizer, constrained polish, acceptance,
bounds, warm/cold tie and strict cutoff are inherited. Unchanged long-double
conditional row terms and native float64 predictor/chain contributions are
summed in long double, then returned as float64 at the optimizer boundary.
No positive-variance fitting, revised acceptance rule or original-slot update.
"""
import argparse
import contextlib
from dataclasses import replace
from decimal import Decimal, localcontext
import hashlib
import importlib
import importlib.util
import io
import json
import os
from pathlib import Path
import platform
import resource
import sys
import time

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
FIXED = REC / 'check_longdouble_null_points_v1_22087992'
FIXED_CODE = FIXED / 'executed_sources/8dedd618015ec67a_check_longdouble_null_points_v1.py'
POINT = REC / 'repair_p3_gamma_draw_11_22081580'
REFERENCE = REC / 'repair_p3_draw11_decimal_gradient_22087909'
OBSERVED = REC / 'repair_p3_gamma_observed_native_22080548/observed_null.json'
HELPER = ROOT / 'scripts/analysis/in_biopsy_allelic/rec_codex/repair_p3_longdouble_conditional_v1.py'
PINS = {
    FIXED_CODE: '8dedd618015ec67a0d34f8941678e70472b0062d01dd9f7ad799b89fe0a56615',
    FIXED / 'summary.json': '39f15dfc27a527ba2b933b988794a90c5d160210e9fc77b3a203fb5401a452cc',
    FIXED / 'baseline.json': '35d6cb5929dd01456754382fd0330ff16d27ecedf551e5a03eb0a243a8b5062e',
    FIXED / 'full_original_Newton_probe.json': 'ee8fa215f09e74d27e0393fb83b7ae04fb34d87bd895b9a165c4bc78d1ec4ca6',
    FIXED / 'source_input_pins_before.json': 'e5b1cbb5aa7b5f06761635fa3b436a6a002663ac54b980ddee4522fcf471c580',
    FIXED / 'source_input_pins_after.json': 'e5b1cbb5aa7b5f06761635fa3b436a6a002663ac54b980ddee4522fcf471c580',
}
THREADS = ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'NUMEXPR_NUM_THREADS')
INHERITED = ('fit', 'fit_single', 'constrained_polish', 'cold_start', 'bounds', 'initial', 'informative', 'predictors')
CAP = 32 * 1024**2
SECONDS = 88 * 3600


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
    require(Path(path).stat().st_size <= 2 * 1024**2, 'JSON size cap exceeded')
    return json.loads(Path(path).read_text())


def module_at(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def metadata():
    pins = {str(path): digest for path, digest in PINS.items()}
    for path, digest in pins.items():
        require(sha(path) == digest, 'Completed fixed-point receipt/source differs')
    before = load(FIXED / 'source_input_pins_before.json')
    require(before == load(FIXED / 'source_input_pins_after.json'), 'Completed fixed-point input pins differ')
    for path, digest in before.items():
        require(path not in pins or pins[path] == digest, 'Conflicting original input pin')
        pins[path] = digest
    require(pins.get(str(HELPER)) == '66c269cf1679ded3461564537429988925670a05ff1beb7c3b94e8625009031d' and
            pins.get(str(OBSERVED)) == 'b980a5fbde75b6053ba7012a17df4cb7f7cd2c6d82d60a5cfbe918a5818e68d8',
            'Exact helper and original observed warm-start pins required')
    for path, digest in pins.items():
        require(sha(path) == digest, 'Original synthetic source/input changed')
    return pins


class Output:
    def __init__(self, directory):
        self.directory, self.saved, self.started = directory, 0, time.monotonic()

    def check(self):
        actual = sum(path.stat().st_size for path in self.directory.rglob('*') if path.is_file())
        require(actual == self.saved and actual < CAP - 65536, '32MiB output accounting/cap differs')
        require(time.monotonic() - self.started < SECONDS, '88h pilot elapsed cap exceeded')
        return actual

    def write(self, name, body):
        require(self.saved + len(body) < CAP - 65536, '32MiB output cap exceeded')
        require(time.monotonic() - self.started < SECONDS, '88h pilot elapsed cap exceeded')
        path = self.directory / name
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open('xb') as stream:
            stream.write(body)
        self.saved += len(body)

    def json(self, name, value):
        self.write(name, (json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + '\n').encode())


def plain(np, value):
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {key: plain(np, item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [plain(np, item) for item in value]
    return value


def record_evaluation(np, values):
    ll, gradient, hessian = values
    return {'loglik': float(ll), 'gradient': gradient.tolist(),
            'hessian': hessian.tolist() if hessian is not None else None}


def projection(model, p, gradient, start):
    # Exact original null fitter tolerance, rather than reference checker1e-10.
    bounds = [bound if model.informative()[j] else (start[j], start[j])
              for j, bound in enumerate(model.bounds())]
    scores = [Decimal(x) if isinstance(x, str) else Decimal.from_float(float(x)) for x in gradient]
    active = [((lo is not None and lo == hi) or
               (lo is not None and p[j] <= lo + 1e-8 and scores[j] < 0) or
               (hi is not None and p[j] >= hi - 1e-8 and scores[j] > 0))
              for j, (lo, hi) in enumerate(bounds)]
    projected = [Decimal(0) if flag else score for flag, score in zip(active, scores)]
    maximum = max(map(abs, projected))
    feasible = all((lo is None or p[j] >= lo) and (hi is None or p[j] <= hi)
                   for j, (lo, hi) in enumerate(bounds))
    return {'bounds': plain_tuple(bounds), 'active_coordinates': [j for j, flag in enumerate(active) if flag],
            'projected_gradient': list(map(str, projected)), 'max_projected_gradient': str(maximum),
            'strict_original_1e_6_score_met': maximum < Decimal('0.000001'),
            'original_boundary_tolerance': 1e-8, 'all_coordinates_within_bounds': bool(feasible)}


def plain_tuple(bounds):
    return [[None if v is None else float(v) for v in bound] for bound in bounds]


def model_class(np, original, helper_class, output, counters):
    original_functions = (original.Model.predictors, original.Model.evaluate, original.conditional_terms,
                          original.bb_terms, helper_class.__init__, helper_class.terms,
                          helper_class.rising, helper_class.sigmoid)

    class LongDoubleNullModel(original.Model):
        def evaluate(self, p, rows=None, u=0, need_hessian=True):
            require(rows is None and u == 0, 'Null-only full-row evaluator rejects subsets and nonzero random effects')
            require(time.monotonic() - output.started < SECONDS, '88h pilot elapsed cap exceeded')
            require((self.G, self.Z, self.size) == (8, 7, 28), 'Exact synthetic model scope required')
            require(all(getattr(type(self), name) is getattr(original.Model, name) for name in INHERITED),
                    'An inherited numerical fitting method changed')
            require(original_functions == (original.Model.predictors, original.Model.evaluate, original.conditional_terms,
                original.bb_terms, helper_class.__init__, helper_class.terms, helper_class.rising, helper_class.sigmoid),
                'Original or helper callback changed')
            require(self.predictors.__func__ is original.Model.predictors and
                    tuple(tuple(vars(row).items()) for row in self.rows) == self.fixed_row_signatures,
                    'Predictor callback or fixed rows changed')
            incoming = np.asarray(p)
            require(incoming.shape == (28,) and incoming.dtype == np.float64 and np.isfinite(incoming).all(),
                    'Finite native float64 parameter vector required')
            incoming_bytes = incoming.tobytes()
            point = incoming.copy()
            point.setflags(write=False)
            helpers = {gene: helper_class(self, point, [row for row in self.rows if row.gene == gene])
                       for gene in range(8)}
            ll = np.longdouble(0)
            gradient = np.zeros(28, dtype=np.longdouble)
            hessian = np.zeros((28, 28), dtype=np.longdouble) if need_hessian else None
            for row in self.rows:
                eta, r, jac, eta_hess = self.predictors(row, point, 0., need_hessian=need_hessian)
                val, score, row_h = helpers[row.gene].terms(row, eta, r, need_hessian=need_hessian)
                require(isinstance(val, float) and jac.dtype == score.dtype == np.float64,
                        'Native predictor/conditional float64 interface changed')
                ll += np.longdouble(val)
                gradient += np.asarray(jac.T @ score, dtype=np.longdouble)
                if need_hessian:
                    require(row_h.dtype == eta_hess.dtype == np.float64, 'Native Hessian interface changed')
                    hessian += np.asarray(jac.T @ row_h @ jac + score[0] * eta_hess, dtype=np.longdouble)
            require(incoming.tobytes() == incoming_bytes and not point.flags.writeable and
                    all(not helper.p.flags.writeable and helper.p.tobytes() == helper.parameter_bytes
                        for helper in helpers.values()), 'Fixed point or helper state changed')
            counters['calls'] += 1
            counters['Hessian_calls' if need_hessian else 'gradient_only_calls'] += 1
            counters['row_terms'] += len(self.rows)
            counters['gene_contexts'] += 8
            result = float(ll), np.asarray(gradient, dtype=np.float64), (
                np.asarray(hessian, dtype=np.float64) if need_hessian else None)
            require(all(a is None or np.isfinite(a).all() for a in result), 'Nonfinite revised evaluation')
            return result

    require(all(getattr(LongDoubleNullModel, name) is getattr(original.Model, name) for name in INHERITED),
            'Subclass must inherit every numerical fitting method unchanged')
    return LongDoubleNullModel


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true')
    parser.add_argument('--expected-source-sha', help='Reviewed runner SHA256, required for numeric execution')
    args = parser.parse_args()
    sys.dont_write_bytecode = True
    pins = metadata()
    if args.inspect:
        require('numpy' not in sys.modules and 'scipy' not in sys.modules, 'Inspect imported numerical libraries')
        print(json.dumps({'pins_verified': len(pins), 'numerical_imports': False, 'null_only': True,
                          'original_fitter_methods_inherited': list(INHERITED)}))
        return 0
    job = os.environ.get('SLURM_JOB_ID', '')
    require(job.isdigit() and int(job) > 0 and not os.environ.get('SLURM_ARRAY_TASK_ID') and
            os.environ.get('SLURM_CPUS_PER_TASK') == '1', 'Individual1CPU compute allocation required')
    require(all(os.environ.get(key) == '1' for key in THREADS), 'Four numerical thread variables must equal1')
    require(args.expected_source_sha and sha(Path(__file__)) == args.expected_source_sha, 'Reviewed runner SHA required')
    pins[str(Path(__file__).resolve())] = args.expected_source_sha
    directory = REC / ('continuation_longdouble_null_fit_v1_' + job)
    directory.mkdir(mode=0o750)
    output = Output(directory)
    phase, counters, terminals = 'archive_and_environment', dict(calls=0, Hessian_calls=0, gradient_only_calls=0,
                                                                row_terms=0, gene_contexts=0), []
    try:
        resource.setrlimit(resource.RLIMIT_CPU, (SECONDS, 89 * 3600))
        archived = set()
        for path, digest in pins.items():
            if Path(path).suffix == '.py':
                name = 'executed_sources/' + digest[:16] + '_' + Path(path).name
                if name not in archived:
                    output.write(name, Path(path).read_bytes())
                    archived.add(name)
        output.json('source_input_pins_before.json', pins)
        import numpy as np
        import scipy
        require((platform.python_version(), np.__version__, scipy.__version__) == ('3.10.18', '2.2.6', '1.15.2') and
                np.finfo(np.longdouble).nmant >= 63, 'Recorded environment and extended significand required')
        config = io.StringIO()
        with contextlib.redirect_stdout(config):
            np.show_config()
        output.json('environment.json', {'Python': platform.python_version(), 'NumPy': np.__version__,
            'SciPy': scipy.__version__, 'Python_binary_sha256': sha(Path(sys.executable).resolve()),
            'NumPy_config': config.getvalue(), 'longdouble_nmant': np.finfo(np.longdouble).nmant,
            'threads': {key: os.environ.get(key) for key in THREADS}, 'seed': None,
            'seed_reason': 'Fixed synthetic counts and deterministic original two-start optimizer',
            'limits': {'CPUs': 1, 'memory_GB': 2, 'wall_hours': 90, 'checked_elapsed_hours': 88,
                       'saved_bytes': CAP, 'per_gene_rows': 32, 'per_gene_count_entries': 100000, 'per_row_n': 5000},
            'scheduler_mem_per_node_MB': os.environ.get('SLURM_MEM_PER_NODE')})
        fixed = module_at('completed_longdouble_fixed_checker', FIXED_CODE)
        decimal_gradient = module_at('saved_draw11_decimal_gradient', REFERENCE / 'executed_checker.py')
        decimal_base = module_at('saved_draw11_decimal_probability', decimal_gradient.DECIMAL)
        archive = POINT / 'executed_sources'
        sys.path.insert(0, str(archive))
        loader = importlib.import_module('bfix_v1')
        importlib.import_module('tests_v1')
        for name in ('bfix_v1', 'design_v1', 'likelihood_v1', 'tests_v1'):
            require(Path(sys.modules[name].__file__).resolve() == (archive / (name + '.py')).resolve(), 'Import escaped archive')
        original = sys.modules['likelihood_v1']
        design, settings = loader.load()
        template, counts = design.model(), load(POINT / 'draw_counts.json')
        require(len(template.rows) == len(counts['counts']) == 173 and len(design.participants) == 60 and
                counts['full_row_keys'] == [[r.individual, settings['gene_order'][r.gene]] for r in template.rows],
                'Original synthetic row order/dimensions differ')
        require(all(isinstance(v, int) and not isinstance(v, bool) for v in counts['counts']) and
                hashlib.sha256(np.asarray(counts['counts'], dtype='<i8').tobytes()).hexdigest() == counts['count_vector_sha256'],
                'Exact original integer count vector differs')
        rows = [replace(row, a=a) for row, a in zip(template.rows, counts['counts'])]
        scope = [{'gene': g, 'rows': sum(r.gene == g for r in rows),
                  'count_entries': sum(r.n + 1 for r in rows if r.gene == g)} for g in range(8)]
        require(all(0 < s['rows'] <= 32 and s['count_entries'] <= 100000 for s in scope) and
                all(0 < r.n <= 5000 for r in rows), 'Unchanged helper row/workspace caps exceeded')
        helper_class = module_at('unchanged_longdouble_conditional', HELPER).LongDoubleConditional
        Revised = model_class(np, original, helper_class, output, counters)
        model = Revised(rows, gene_count=8)
        native = original.Model(rows, gene_count=8)
        require(tuple(rows) == model.rows == native.rows, 'Original row order changed')
        model.fixed_row_signatures = tuple(tuple(vars(row).items()) for row in model.rows)
        observed = np.asarray(load(OBSERVED)['fit']['parameters'], dtype=np.float64)
        require(observed.shape == (28,) and np.isfinite(observed).all(), 'Original observed warm start differs')
        output.json('protocol.json', {'draw': 11, 'count_vector_sha256': counts['count_vector_sha256'], 'rows': 173,
            'source_individuals': 60, 'gene_scope': scope, 'observed_warm_start': observed.tolist(),
            'method_source_file': str(archive / 'likelihood_v1.py'), 'inherited_methods': list(INHERITED),
            'null_only': True, 'subsets_and_nonzero_u_rejected': True,
            'no_parameter_or_row_cache': True, 'fresh_context_per_gene_per_evaluation': True,
            'native_predictor_and_row_return_and_chain_rule': 'float64', 'aggregation': 'longdouble',
            'optimizer_interface': 'float64', 'fully_extended_precision_evaluation': False,
            'cold_start_original_direct_conditional_path_unchanged': True,
            'instrumentation': 'Explicit instance fit_single wrapper delegates to saved original inherited bound method; no global profiling/hooks',
            'no_positive_variance_fit': True, 'original_failures_unchanged': True})
        phase = 'actual_revised_evaluator_prefit_checks'
        for label in ('baseline', 'full_original_Newton_probe'):
            receipt = load(FIXED / (label + '.json'))
            p = np.asarray(receipt['parameters'], dtype=np.float64)
            candidate = record_evaluation(np, model.evaluate(p))
            native_result = record_evaluation(np, native.evaluate(p))
            comparisons = {}
            for digits in (80, 110):
                reference = load(REFERENCE / (label + '_' + str(digits) + '.json'))
                require(reference['parameters'] == p.tolist(), 'Decimal fixed-point vector differs')
                comparison = fixed.compare(np, candidate, reference, False)
                comparison['LL_atol_1e_6_met'] = abs(Decimal(comparison['loglik_difference'])) <= Decimal('0.000001')
                comparison['gradient_atol_1e_8_rtol_1e_6_met'] = bool(np.allclose(
                    candidate['gradient'], np.asarray(reference['gradient'], dtype=np.float64), atol=1e-8, rtol=1e-6))
                comparisons[str(digits)] = comparison
            output.json('prefit_' + label + '.json', {'point': label, 'parameters': p.tolist(), 'actual_revised_evaluate': candidate,
                'Decimal_comparisons': comparisons, 'native_comparison': fixed.compare(np, candidate, native_result, True),
                'independent_Hessian': False, 'native_Hessian_check_descriptive_only': True})
            require(all(v['LL_atol_1e_6_met'] and v['gradient_atol_1e_8_rtol_1e_6_met'] for v in comparisons.values()),
                    'Actual revised evaluator failed fixed-point Decimal LL/gradient requirements')
        original_single = model.fit_single
        require(original_single.__func__ is original.Model.fit_single and
                all(getattr(Revised, name) is getattr(original.Model, name) for name in INHERITED),
                'Inherited original fitter identity differs before instrumentation')

        def recorded_single(start):
            require(len(terminals) < 2, 'Unexpected extra original fit_single call')
            label = ('warm', 'cold')[len(terminals)]
            exact_start = np.asarray(start).copy()
            output.json(label + '_start.json', {'parameters': exact_start.tolist(), 'calls_before': dict(counters),
                'source_bound_method': 'original.Model.fit_single', 'no_optimizer_or_acceptance_modification': True})
            fitted = original_single(start)
            require(np.array_equal(np.asarray(start), exact_start), 'Original fit_single mutated its input start')
            terminal = {'start': label, 'exact_start': exact_start.tolist(), 'fit': plain(np, fitted),
                        'calls_after': dict(counters)}
            output.json(label + '_terminal.json', terminal)
            terminals.append(terminal)
            output.check()
            return fitted

        model.fit_single = recorded_single
        phase = 'original_two_start_null_fit'
        warm_bytes = observed.tobytes()
        fit = model.fit(observed)
        require(len(terminals) == 2 and observed.tobytes() == warm_bytes and
                original_single.__func__ is original.Model.fit_single, 'Two-start instrumentation/observed input changed')
        output.json('null_fit.json', plain(np, fit))
        # Independently rederive the original float64 warm/cold tie, without
        # selecting an endpoint from Decimal or another arithmetic's likelihood.
        expected = 'cold' if terminals[1]['fit']['loglik'] - terminals[0]['fit']['loglik'] >= 1e-6 else 'warm'
        require(fit['kept_start'] == expected and np.array_equal(fit['parameters'], terminals[('warm', 'cold').index(expected)]['fit']['parameters']),
                'Original warm/cold tie or retained parameter identity differs')
        phase = 'native_and_independent_Decimal_endpoint_scores'
        endpoint_summaries = {}
        for terminal in terminals:
            label, saved = terminal['start'], terminal['fit']
            p = np.asarray(saved['parameters'], dtype=np.float64)
            start = terminal['exact_start']
            candidate = record_evaluation(np, model.evaluate(p))
            native_result = record_evaluation(np, native.evaluate(p))
            require(candidate == {key: saved[key] for key in ('loglik', 'gradient', 'hessian')},
                    'Revised terminal evaluation does not reproduce original fit_single receipt')
            projections = {'revised': projection(model, p, candidate['gradient'], start),
                           'native': projection(model, p, native_result['gradient'], start)}
            decimal_records, comparisons = {}, {}
            for digits in (80, 110):
                independent = decimal_gradient.gradient(decimal_base, model.rows, p.tolist(), digits)
                output.json(label + '_Decimal_' + str(digits) + '.json', independent)
                projections['Decimal_' + str(digits)] = projection(model, p, independent['gradient'], start)
                comparisons[str(digits)] = fixed.compare(np, candidate, independent, False)
                decimal_records[str(digits)] = {key: independent[key] for key in ('digits', 'loglik', 'gradient')}
            require(saved['converged'] == projections['revised']['strict_original_1e_6_score_met'] and
                    all(v['all_coordinates_within_bounds'] for v in projections.values()), 'Original endpoint cutoff/bounds differ')
            endpoint_summaries[label] = {'parameters': p.tolist(), 'native_recomputed': native_result,
                'actual_revised_evaluate': candidate, 'Decimal_reference_scores': decimal_records, 'projections': projections,
                'Decimal_comparisons': comparisons, 'native_comparison': fixed.compare(np, candidate, native_result, True),
                'independent_Hessian': False, 'fit_converged_original_rule': saved['converged']}
            output.json(label + '_endpoint_checks.json', endpoint_summaries[label])
        phase = 'final_source_and_output_checks'
        require(all(getattr(Revised, name) is getattr(original.Model, name) for name in INHERITED), 'Inherited fitter changed')
        for path, digest in pins.items():
            require(sha(path) == digest, 'Source/input changed during separate pilot')
        output.json('source_input_pins_after.json', pins)
        output.json('summary.json', {'completed': True, 'SLURM_job_id': job, 'draw': 11, 'null_only': True,
            'kept_start': fit['kept_start'], 'original_float64_tie_threshold': 1e-6,
            'start_logliks': plain(np, fit['start_logliks']), 'start_converged': plain(np, fit['start_converged']),
            'endpoint_checks': endpoint_summaries, 'call_counts': counters,
            'max_memory_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'elapsed_seconds': time.monotonic() - output.started, 'saved_bytes_before_summary': output.check(),
            'saved_cap_bytes': CAP, 'source_input_pins_before_after_identical': True,
            'original_inherited_methods_unchanged': list(INHERITED), 'only_numerical_change': 'evaluate arithmetic',
            'precision_boundaries': 'Unchanged internal longdouble conditional; native float64 predictor/row returns/chain products; longdouble aggregation; float64 optimizer return',
            'independent_Hessian_or_mathematical_certification': False, 'positive_variance_fitting': False,
            'original_failed_draw_LR_exceedance_and_warm_cold_rules_replaced': False,
            'calibration_or_biological_validation': False, 'restricted_or_receiving_values_read': False})
        print(json.dumps({'completed': True, 'summary_sha256': sha(directory / 'summary.json'),
                          'saved_bytes': output.check(), 'original_failure_unchanged': True}), flush=True)
        return 0
    except Exception as exc:
        failure = {'completed': False, 'phase': phase, 'exception_type': type(exc).__name__,
            'reason': str(exc) if isinstance(exc, RuntimeError) else type(exc).__name__,
            'max_memory_KiB': resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            'elapsed_seconds': time.monotonic() - output.started, 'call_counts': counters,
            'completed_start_receipts': [t['start'] for t in terminals], 'partials_preserved': True,
            'original_failed_draw_and_results_unchanged': True}
        with (directory / 'failure.json').open('x') as stream:
            json.dump(failure, stream, indent=2)
            stream.write('\n')
        print(json.dumps(failure), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
