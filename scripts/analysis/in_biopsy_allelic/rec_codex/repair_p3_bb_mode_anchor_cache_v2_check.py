#!/usr/bin/env python3
"""Synthetic full-support reuse comparisons, then the saved refused integral."""
import argparse
import ast
from decimal import Decimal, localcontext
import hashlib
import importlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import resource
import sys
import threading
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
POINT = REC/'repair_p3_observed_failure_anchor_22042612'
THETA = REC/'repair_p3_theta4_anchored_gl_gene2_22042347'
CASES = REC/'repair_p3_bb_mode_anchor_batch_fast_launch_20261005T172617Z/cases.json'
MODULES = ('likelihood_v1.py', 'design_v1.py', 'bfix_v1.py', 'tests_v1.py')
HELPERS = ('repair_p3_bb_mode_anchor.py', 'repair_p3_bb_mode_anchor_cache_v2.py',
           'repair_p3_composite_gl.py', 'repair_p3_integral.py')
KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
COUNT_SHA = '7009be4adeec6157fa4bbf3f46e48ec24aa398a0bee4d7906c78c231b7d3fbad'
CAP = 32*1024**2


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--expected-checker-sha', required=True)
    parser.add_argument('--pins', type=Path, required=True)
    parser.add_argument('--expected-pins-sha', required=True)
    parser.add_argument('--inspect', action='store_true')
    parser.add_argument('--component-only', action='store_true')
    args = parser.parse_args()
    require(sha(Path(__file__)) == args.expected_checker_sha and
            sha(args.pins) == args.expected_pins_sha, 'Reviewed checker/manifest differs')
    pins = {Path(p): h for p, h in json.loads(args.pins.read_text()).items()}
    for path, digest in pins.items():
        require(path.is_absolute() and path.resolve().is_relative_to(ROOT.resolve())
                and sha(path) == digest, 'Reviewed source/input differs: '+str(path))
        if path.suffix == '.py':
            ast.parse(path.read_text())
    required = [*(HERE/n for n in (*MODULES, *HELPERS)), CASES,
                *(POINT/n for n in ('protocol.json', 'summary.json', 'GL_candidate.json',
                                   'Decimal_reference_summary.json', 'GL_minus_Decimal.json')),
                THETA/'protocol.json']
    require(all(p in pins for p in required), 'Required source/input pin missing')
    cache_tree = ast.parse((HERE/'repair_p3_bb_mode_anchor_cache_v2.py').read_text())
    scalar_sha = next(ast.literal_eval(node.value) for node in cache_tree.body
                      if isinstance(node, ast.Assign) and any(isinstance(target, ast.Name)
                      and target.id == 'SCALAR_SHA' for target in node.targets))
    require(scalar_sha == pins[HERE/'repair_p3_bb_mode_anchor.py'],
            'Reuse helper scalar pin differs before allocation')
    parent = json.loads((POINT/'protocol.json').read_text())
    summary = json.loads((POINT/'summary.json').read_text())
    theta = json.loads((THETA/'protocol.json').read_text())
    packet = json.loads(CASES.read_text())
    require(summary['completed'] is True and summary['GL_candidate_refused'] is False
            and parent['gene'] == 2 and parent['call'] == 12 and parent['target_rows'] == 24
            and parent['fixed_v'] == .0009636462009520182
            and parent['count_vector_sha256'] == COUNT_SHA
            and theta['fixed_v'] == .0004 and len(packet['cases']) == 21,
            'Exact completed synthetic reference family differs')
    protocol = dict(source_and_input_sha256={str(p): h for p, h in pins.items()},
        selected_BB_points=21, conditional_rows=173, synthetic_participants=60,
        conditional_vectors=['theta4', 'failed_call12'], u_values=['0', '+sqrt(v)', '-sqrt(v)'],
        full_GL_gene=2, full_GL_with_hessian=False, original_GL_rules_and_caps=True,
        timing_rows=24, timing_u_points=32, timing_repeats=3,
        maximum_retained_support_bytes=8*1024**2, output_cap_bytes=CAP,
        no_fit=True, no_calibration=True, no_adoption=True)
    protocol['component_only'] = args.component_only
    protocol['full_GL_gene'] = None if args.component_only else 2
    if args.inspect:
        print(json.dumps(protocol))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and
            not os.environ.get('SLURM_ARRAY_TASK_ID'), 'Individual compute allocation required')
    output = REC/('repair_p3_bb_mode_anchor_cache_v2_check_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    archive = output/'executed_sources'
    archive.mkdir()
    for source in (Path(__file__), *(HERE/n for n in (*MODULES, *HELPERS))):
        target = archive/source.name
        target.write_bytes(source.read_bytes())
        require(sha(target) == (args.expected_checker_sha if source == Path(__file__)
                               else pins[source]), 'Executable archive differs')
    began = time.perf_counter()
    try:
        import numpy as np
        import scipy

        def serialize(value):
            if isinstance(value, np.ndarray):
                return value.tolist()
            if isinstance(value, np.generic):
                return value.item()
            raise TypeError(type(value).__name__)

        def save(name, value):
            body = (json.dumps(value, default=serialize, allow_nan=False, indent=2)+'\n').encode()
            require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body)
                    < CAP-65536, 'Output capacity exceeded; refusal reserve retained')
            with (output/(name+'.json')).open('xb') as handle:
                handle.write(body)
                handle.flush()
                os.fsync(handle.fileno())
            print(json.dumps(dict(saved=name, elapsed_seconds=time.perf_counter()-began)), flush=True)

        versions = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)
        require(versions == parent['environment'], 'Numerical environment differs')
        protocol['environment'] = versions
        save('protocol', protocol)
        sys.path.insert(0, str(archive))
        require(not any(n[:-3] in sys.modules for n in MODULES), 'Original modules already imported')
        loader = importlib.import_module('bfix_v1')
        importlib.import_module('tests_v1')
        for name in MODULES:
            module = sys.modules[name[:-3]]
            require(Path(module.__file__).resolve().parent == archive.resolve()
                    and sha(Path(module.__file__)) == pins[HERE/name], 'Original import escaped archive')
        private = load('cache_check_original', archive/'likelihood_v1.py')
        anchor = load('cache_check_scalar', archive/'repair_p3_bb_mode_anchor.py')
        reuse = load('cache_check_reuse', archive/'repair_p3_bb_mode_anchor_cache_v2.py')
        gl = load('cache_check_GL', archive/'repair_p3_composite_gl.py')
        cache = reuse.OnePointBBTerms(anchor.stable_bb_terms)
        design, settings = loader.load()
        model = private.Model(design.model().rows, gene_count=8)
        require(len(model.rows) == 173 and model.size == 28 and model.Z == 7
                and len(design.participants) == 60 and settings['seed'] == 20260930
                and hashlib.sha256(np.asarray([r.a for r in model.rows], dtype='<i8').tobytes()).hexdigest()
                == COUNT_SHA, 'Exact synthetic fixture differs')
        selected = []
        for case in packet['cases']:
            lp, g, h = cache([case['k']], case['n'], case['eta'], case['r'])
            ladder = case['reference_ladder']
            require([(z['digits'], z['step']) for z in ladder] ==
                    [(d, s) for d in (80, 110) for s in ('1e-15', '1e-20')], 'Independent ladder differs')
            with localcontext() as context:
                context.prec = 110
                flat = lambda z: [Decimal(z['logpmf']), *map(Decimal, z['score']),
                                  *(Decimal(v) for row in z['hessian'] for v in row)]
                final = flat(ladder[-1])
                stability = max(abs(a-b) for z in ladder for a, b in zip(flat(z), final))
                differences = [Decimal.from_float(float(a))-b
                               for a, b in zip([lp[0], *g[0], *h[0].ravel()], final)]
                passed = (abs(differences[0]) <= Decimal('1e-9') and
                          max(map(abs, differences[1:3])) <= Decimal('1e-8') and
                          max(map(abs, differences[3:])) <= Decimal('1e-5') and
                          stability <= Decimal('1e-25'))
            selected.append(dict(origin=case['origin'], n=case['n'], k=case['k'], eta=case['eta'], r=case['r'],
                raw_candidate_minus_Decimal=list(map(str, differences)), Decimal_stability=str(stability),
                selected_criterion_met=bool(passed)))
        save('selected_comparisons', selected)
        require(all(r['selected_criterion_met'] for r in selected), 'Selected independent BB criterion failed')

        def identical(a, b):
            return all((x is None and y is None) or (x is not None and y is not None
                       and np.asarray(x).dtype == np.asarray(y).dtype
                       and np.asarray(x).shape == np.asarray(y).shape
                       and np.asarray(x).tobytes() == np.asarray(y).tobytes())
                       for x, y in zip(a, b))

        # Hit/invalidation, n=0, empty/duplicate requests and returned-array mutation.
        cache = reuse.OnePointBBTerms(anchor.stable_bb_terms)
        for n, eta, r, need_hessian in [(5, 0., -15., True), (5, 0., -15., True),
                (5, 0., -15., False), (5, .1, -15., False), (5, .1, -15., True), (5, .1, -14., True),
                (6, .1, -14., True), (0, 0., 0., True), (5, 1200., -3., True),
                (5, -1200., -3., True), (5, -0., -3., True), (5, 0., -3., True), (5, 0., 2., True)]:
            counts = [n, 0, n//2, n//2]
            reference = anchor.stable_bb_terms(counts, n, eta, r, need_hessian)
            value = cache(counts, n, eta, r, need_hessian)
            require(identical(value, reference), 'Support reuse/invalidation differs')
            for array in value:
                if array is not None:
                    array.fill(123.)
            require(identical(cache(counts, n, eta, r, need_hessian), reference), 'Returned mutation poisoned support')
            require(identical(cache([], n, eta, r, need_hessian),
                              anchor.stable_bb_terms([], n, eta, r, need_hessian)), 'Empty request differs')
        for counts, n, eta, r in [([-1], 5, 0., 0.), ([6], 5, 0., 0.), ([.5], 5, 0., 0.),
                ([0], 5.5, 0., 0.), ([[0]], 5, 0., 0.), ([0], 5, math.nan, 0.)]:
            try:
                cache(counts, n, eta, r)
            except ValueError:
                continue
            raise ValueError('Invalid input accepted by reuse helper')
        for n, eta, r in [(-1, 0., 0.), (100001, 0., 0.), (5, 0., math.inf),
                          (5, np.longdouble(0.), 0.), (5, np.float32(0.), 0.),
                          (5, np.int64(0), 0.), (5, 2**53+1, 0.)]:
            try:
                cache([0], n, eta, r)
            except ValueError:
                continue
            raise ValueError('Unsupported depth/predictor accepted')
        thread_result = []
        def another_thread():
            try:
                cache([0], 5, 0., 0.)
            except ValueError as error:
                thread_result.append(str(error))
        worker = threading.Thread(target=another_thread)
        worker.start()
        worker.join()
        require(len(thread_result) == 1 and 'constructing process and thread' in thread_result[0],
                'A second thread could reuse the support')
        stats = cache.statistics()
        require(stats['calls'] == 39 and stats['misses'] == 12 and stats['hits'] == 27,
                'Expected actual reuse/invalidation counts differ')
        save('reuse_invariants', dict(passed=True, statistics=stats, thread_ownership_checked=True,
                                     unsupported_predictor_types_rejected=True))
        conditional = []
        for label, source in [('theta4', theta), ('failed_call12', parent)]:
            p, v = np.asarray(source['fixed_p']), source['fixed_v']
            for u in (0., math.sqrt(v), -math.sqrt(v)):
                for need_hessian in (False, True):
                    for row in model.rows:
                        eta, r = model.predictors(row, p, u, need_hessian=need_hessian)[:2]
                        private.bb_terms = anchor.stable_bb_terms
                        a = private.conditional_terms(row, eta, r, need_hessian)
                        private.bb_terms = cache
                        b = private.conditional_terms(row, eta, r, need_hessian)
                        require(identical(a, b), 'Native conditional LP/score/H changed')
                    conditional.append(dict(vector=label, u=u, Hessian=need_hessian, compared_rows=173,
                                            bitwise_identical=True))
        save('conditional_comparisons', dict(records=conditional, statistics=cache.statistics()))
        p, v = np.asarray(parent['fixed_p']), parent['fixed_v']
        rows = tuple(row for row in model.rows if row.gene == 2)
        require(len(rows) == 24, 'Saved failed gene row count differs')
        timing = []
        for need_hessian in (False, True):
            for repeat in range(3):
                samples = []
                durations = []
                helpers = (anchor.stable_bb_terms, cache) if (repeat+int(need_hessian))%2 == 0 else (cache, anchor.stable_bb_terms)
                for helper in helpers:
                    private.bb_terms = helper
                    start = time.perf_counter()
                    values = []
                    for u in np.linspace(-math.sqrt(v), math.sqrt(v), 32):
                        for row in rows:
                            eta, r = model.predictors(row, p, float(u), need_hessian=need_hessian)[:2]
                            values.append(private.conditional_terms(row, eta, r, need_hessian))
                    durations.append(time.perf_counter()-start)
                    samples.append(values)
                if helpers[0] is cache:
                    durations.reverse()
                require(all(identical(a, b) for a, b in zip(*samples)), 'Timed conditional results differ')
                timing.append(dict(Hessian=need_hessian, repeat=repeat, scalar_seconds=durations[0],
                    reuse_seconds=durations[1], scalar_over_reuse=durations[0]/durations[1],
                    scalar_first=helpers[0] is anchor.stable_bb_terms,
                    bitwise_identical=True, whole_fit_speed_claim=False))
        save('timing', dict(records=timing, statistics=cache.statistics(),
            measured_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        if args.component_only:
            for path, digest in pins.items():
                require(sha(path) == digest, 'Source/input changed during run')
            save('summary', dict(completed=True, selected21_criterion_met=True,
                conditional_comparisons=len(conditional)*173, conditional_LP_score_H_bitwise_identical=True,
                no_integration=True, statistics=cache.statistics(),
                measured_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                no_full_independent_Hessian=True, no_fit=True, no_calibration=True, no_adoption=True))
            return
        private.bb_terms = cache
        start = time.perf_counter()
        ll, g, h, report = gl.integrate_gene(model, p, v, rows, private.conditional_terms,
            with_hessian=False, expected_kernel_sha=KERNEL_SHA, kernel_path=archive/'repair_p3_integral.py')
        seconds = time.perf_counter()-start
        reference = json.loads((POINT/'GL_candidate.json').read_text())
        save('GL_candidate', dict(logmass=ll, gradient=g, Hessian=h, integration=report,
                                  elapsed_seconds=seconds, statistics=cache.statistics()))
        require(identical((ll, g, h), (reference['logmass'], reference['gradient'], None))
                and len(report['ladder']) == len(reference['integration']['ladder'])
                and report['likelihood_evaluations'] == reference['integration']['likelihood_evaluations'],
                'Full saved integral/value/score/rule count differs')
        for a, b in zip(report['ladder'], reference['integration']['ladder']):
            require(all(a[k] == b[k] for k in ('z_bound', 'grid_spacing', 'order', 'node_evaluations', 'logmass'))
                    and identical((a['gradient'],), (b['gradient'],)), 'Full refinement ladder differs')
        for path, digest in pins.items():
            require(sha(path) == digest, 'Source/input changed during run')
        save('summary', dict(completed=True, selected21_criterion_met=True, conditional_comparisons=len(conditional)*173,
            conditional_LP_score_H_bitwise_identical=True, full_saved_GL_LP_score_and_ladder_bitwise_identical=True,
            GL_elapsed_seconds=seconds, statistics=cache.statistics(),
            measured_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            no_full_independent_Hessian=True, no_fit=True, no_calibration=True, no_adoption=True))
    except BaseException:
        with (output/'failure.json').open('x') as handle:
            json.dump(dict(traceback=traceback.format_exc(), partials_preserved=True, no_fit=True,
                           no_calibration=True, no_adoption=True), handle, indent=2)
        raise


if __name__ == '__main__':
    main()
