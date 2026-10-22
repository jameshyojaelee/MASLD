#!/usr/bin/env python3
"""Selected independent checks of a direct gamma-difference native BB candidate."""
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
import time
import traceback

ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT/'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
HERE = ROOT/'scripts/analysis/in_biopsy_allelic/rec_codex'
POINT = REC/'repair_p3_observed_failure_anchor_22042612'
THETA = REC/'repair_p3_theta4_anchored_gl_gene2_22042347'
CASES = REC/'repair_p3_bb_mode_anchor_batch_fast_launch_20261005T172617Z/cases.json'
MODULES = ('likelihood_v1.py', 'design_v1.py', 'bfix_v1.py', 'tests_v1.py')
HELPERS = ('repair_p3_bb_mode_anchor.py','repair_p3_bb_gamma_difference_v3.py')
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
        no_GL=True,independent_rising_shapes=['-1200','log(1e-8)','0','log(31.999)','log(32)','log(32.001)','log(1000)','log(1e6)'],
        independent_rising_n=5000,independent_rising_log_atol='1e-9',independent_rising_derivative_atol='1e-10',
        timing_rows=24, timing_u_points=32, timing_repeats=3,
        maximum_n=5000,output_cap_bytes=CAP,bernoulli_terms=4,
        independent_rising_counts=[0,1,2,999,1000,1001,2500,5000],relative_derivative_tolerance='1e-10',
        no_fit=True, no_calibration=True, no_adoption=True)
    if args.inspect:
        print(json.dumps(protocol))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and
            not os.environ.get('SLURM_ARRAY_TASK_ID'), 'Individual compute allocation required')
    output = REC/('repair_p3_bb_gamma_difference_v3_check_'+os.environ['SLURM_JOB_ID'])
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
        candidate = load('gamma_difference_candidate', archive/'repair_p3_bb_gamma_difference_v3.py')
        require(private.conditional_terms.__globals__ is vars(private)
                and private.mixture_terms.__globals__ is vars(private)
                and private.het_inclusion.__globals__ is vars(private)
                and 'bb_terms' in private.mixture_terms.__code__.co_names
                and 'bb_terms' in private.het_inclusion.__code__.co_names,
                'Native conditional helper does not dynamically resolve the reviewed BB function')
        call_counts = dict(mode_ratio=0,gamma_difference=0)
        def counted(label,helper):
            def wrapper(*args,**kwargs):
                call_counts[label] += 1
                return helper(*args,**kwargs)
            return wrapper
        reference_helper = counted('mode_ratio',anchor.stable_bb_terms)
        cache = counted('gamma_difference',candidate.stable_bb_terms)
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

        rising_records = []
        for logshape in (-1200.,math.log(1e-8),0.,math.log(31.999),math.log(32),
                         math.log(32.001),math.log(1000),math.log(1e6)):
            counts = [0,1,2,999,1000,1001,2500,5000]
            actual = candidate.rising(counts,logshape)
            with localcontext() as context:
                context.prec = 80
                shape = Decimal.from_float(logshape).exp()
                logarithm = first = second = Decimal(0)
                expected = {0:(logarithm,first,second)}
                for j in range(5000):
                    term = shape+Decimal(j)
                    logarithm += term.ln()
                    first += shape/term
                    second += shape*Decimal(j)/(term*term)
                    if j+1 in counts:
                        expected[j+1] = (logarithm,first,second)
                differences = [[Decimal.from_float(float(actual[j][i]))-expected[k][j]
                                for j in range(3)] for i,k in enumerate(counts)]
                relative_differences = [[abs(d[j]/expected[k][j]) if abs(expected[k][j]) > Decimal('1e-250') else None
                    for j in (1,2)] for k,d in zip(counts,differences)]
                passed = all(abs(d[0]) <= Decimal('1e-9') and
                             max(map(abs,d[1:])) <= Decimal('1e-10') for d in differences)
                passed = passed and all(x is None or x <= Decimal('1e-10') for pair in relative_differences for x in pair)
            rising_records.append(dict(logshape=logshape,counts=counts,
                candidate_minus_Decimal=[list(map(str,d)) for d in differences],
                relative_derivative_differences=[[str(x) if x is not None else None for x in pair] for pair in relative_differences],
                evaluated_float_shape=float(np.exp(logshape)),passed=passed))
        save('independent_rising',rising_records)
        require(all(r['passed'] for r in rising_records),'Independent rising-product criterion failed')

        def identical(a,b):
            return all((x is None and y is None) or (x is not None and y is not None
                and np.allclose(x,y,atol=tol,rtol=0)) for x,y,tol in zip(a,b,(1e-9,1e-8,1e-5)))

        conditional = []
        for label, source in [('theta4', theta), ('failed_call12', parent)]:
            p, v = np.asarray(source['fixed_p']), source['fixed_v']
            for u in (0., math.sqrt(v), -math.sqrt(v)):
                for need_hessian in (False, True):
                    maxima = np.zeros(3)
                    before_counts = call_counts.copy()
                    for row in model.rows:
                        eta, r = model.predictors(row, p, u, need_hessian=need_hessian)[:2]
                        private.bb_terms = reference_helper
                        a = private.conditional_terms(row, eta, r, need_hessian)
                        private.bb_terms = cache
                        b = private.conditional_terms(row, eta, r, need_hessian)
                        require(identical(a, b), 'Native conditional LP/score/H exceeds unchanged nominal criterion')
                        maxima = np.maximum(maxima,[float(np.max(np.abs(x-y)))
                            if x is not None else 0. for x,y in zip(a,b)])
                    delta_counts = {key:call_counts[key]-before_counts[key] for key in call_counts}
                    require(all(n > 0 for n in delta_counts.values()),'Native conditional comparison did not call both helpers')
                    conditional.append(dict(helper_call_counts=delta_counts,vector=label, u=u, Hessian=need_hessian, compared_rows=173,
                                            within_nominal_criteria=True,maximum_absolute_differences=maxima))
        save('conditional_comparisons', dict(records=conditional))
        p, v = np.asarray(parent['fixed_p']), parent['fixed_v']
        rows = tuple(row for row in model.rows if row.gene == 2)
        require(len(rows) == 24, 'Saved failed gene row count differs')
        timing = []
        for need_hessian in (False, True):
            for repeat in range(3):
                samples,durations = {},{}
                helpers = [('mode_ratio',reference_helper),('gamma_difference',cache)]
                if repeat%2:
                    helpers.reverse()
                for label,helper in helpers:
                    private.bb_terms = helper
                    start = time.perf_counter()
                    values = []
                    for u in np.linspace(-math.sqrt(v), math.sqrt(v), 32):
                        for row in rows:
                            eta, r = model.predictors(row, p, float(u), need_hessian=need_hessian)[:2]
                            values.append(private.conditional_terms(row, eta, r, need_hessian))
                    durations[label] = time.perf_counter()-start
                    samples[label] = values
                require(all(identical(a, b) for a, b in zip(samples['mode_ratio'],samples['gamma_difference'])),
                        'Timed conditional results exceed nominal criteria')
                timing.append(dict(Hessian=need_hessian,repeat=repeat,
                    helper_order=[label for label,helper in helpers],durations=durations,
                    mode_ratio_over_gamma_difference=durations['mode_ratio']/durations['gamma_difference'],
                    within_nominal_criteria=True,whole_fit_speed_claim=False))
        save('timing', dict(records=timing,
            measured_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
        for path,digest in pins.items():
            require(sha(path) == digest,'Source/input changed during run')
        save('summary',dict(completed=True,selected21_criterion_met=True,
            independent_rising_criterion_met=True,conditional_comparisons=len(conditional)*173,
            conditional_nominal_LP_score_H_criteria_met=True,native_helper_call_counts_checked=True,
            measured_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            no_full_independent_Hessian=True,no_GL=True,no_fit=True,no_calibration=True,no_adoption=True))
    except BaseException:
        with (output/'failure.json').open('x') as handle:
            json.dump(dict(traceback=traceback.format_exc(), partials_preserved=True, no_fit=True,
                           no_calibration=True, no_adoption=True), handle, indent=2)
        raise


if __name__ == '__main__':
    main()

