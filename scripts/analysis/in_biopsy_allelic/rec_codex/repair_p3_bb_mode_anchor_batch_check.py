#!/usr/bin/env python3
"""Selected native RNA allelic BB comparisons and bounded batching timings."""
import argparse
import ast
from decimal import Decimal, localcontext
import hashlib
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
HERE = Path(__file__).resolve().parent
SCALAR_SHA = 'a31258aaac0fb39e91dda7aa634f0e810a6d0da144511b34d9c1b73a7028c266'
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
    parser.add_argument('--expected-batch-sha', required=True)
    parser.add_argument('--cases', type=Path, required=True)
    parser.add_argument('--expected-cases-sha', required=True)
    parser.add_argument('--inspect', action='store_true')
    args = parser.parse_args()
    pins = {Path(__file__).resolve(): args.expected_checker_sha,
            HERE/'repair_p3_bb_mode_anchor.py': SCALAR_SHA,
            HERE/'repair_p3_bb_mode_anchor_batch.py': args.expected_batch_sha,
            args.cases: args.expected_cases_sha}
    for path, digest in pins.items():
        require(sha(path) == digest, 'Reviewed source/case manifest changed')
        if path.suffix == '.py':
            ast.parse(path.read_text())
    packet = json.loads(args.cases.read_text())
    cases = packet['cases']
    require(len(cases) == 21 and packet['old_points'] == 18 and packet['observed_failure_points'] == 3,
            'Exact selected21-point family required')
    for relative, digest in packet['input_sha256'].items():
        path = REC/relative
        require(path.resolve().is_relative_to(REC.resolve()) and sha(path) == digest,
                'Selected independent reference receipt changed')
        pins[path] = digest
    for case in cases:
        source = json.loads((REC/case['receipt']).read_text())
        point = source.get('point', source)
        ladder = source['reference']['value'] if 'reference' in source else source['Decimal_reference']
        require(source['selected_family_criterion_met'] is True and
                all(case[k] == point[k] for k in ('n', 'k', 'eta', 'r')) and
                case['reference_ladder'] == ladder, 'Extracted selected point/reference differs')
    protocol = dict(selected_points=21, comparison_scope='21 selected BB points only',
        criteria=dict(logpmf='1e-9', score='1e-8', Hessian='1e-5', Decimal_stability='1e-25'),
        timing_case=dict(n=2466, batch=32, r=-15., repetitions=3, Hessian=True),
        deterministic_cases=True, output_cap_bytes=CAP, no_integration=True,
        no_fit=True, no_calibration=True, no_adoption=True,
        source_and_input_sha256={str(p): h for p, h in pins.items()})
    if args.inspect:
        print(json.dumps(protocol))
        return
    require(os.environ.get('SLURM_JOB_ID', '').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    output = REC/('repair_p3_bb_mode_anchor_batch_check_'+os.environ['SLURM_JOB_ID'])
    output.mkdir(exist_ok=False)
    archive = output/'executed_sources'
    archive.mkdir()
    for path in tuple(p for p in pins if p.suffix == '.py'):
        target = archive/path.name
        target.write_bytes(path.read_bytes())
        require(sha(target) == pins[path], 'Executable archive differs')
    (output/'input_cases.json').write_bytes(args.cases.read_bytes())
    try:
        import numpy as np
        import scipy
        versions = dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__)
        require(versions == dict(python='3.10.18', numpy='2.2.6', scipy='1.15.2'), 'Environment differs')
        scalar = load('batch_check_scalar', archive/'repair_p3_bb_mode_anchor.py')
        batch = load('batch_check_batch', archive/'repair_p3_bb_mode_anchor_batch.py')

        def serialize(value):
            if isinstance(value, np.ndarray):
                return value.tolist()
            if isinstance(value, np.generic):
                return value.item()
            raise TypeError(type(value).__name__)

        def save(name, value):
            body = (json.dumps(value, default=serialize, allow_nan=False, indent=2)+'\n').encode()
            require(sum(p.stat().st_size for p in output.rglob('*') if p.is_file())+len(body) < CAP,
                    'Output cap exceeded')
            with (output/(name+'.json')).open('xb') as handle:
                handle.write(body)

        protocol['environment'] = versions
        save('protocol', protocol)
        groups = {}
        for case in cases:
            groups.setdefault((case['n'], case['r']), []).append(case)
        results = []
        for index, ((n, r), group) in enumerate(groups.items()):
            ks = list(dict.fromkeys(case['k'] for case in group))
            etas = [case['eta'] for case in group]
            lp, g, h = batch.stable_bb_terms_batch(ks, n, etas, r)
            require(lp.shape == (len(group), len(ks)) and g.shape == lp.shape+(2,) and
                    h.shape == lp.shape+(2, 2) and all(np.isfinite(v).all() for v in (lp, g, h)),
                    'Finite batch shape required')
            lp0, g0, h0 = batch.stable_bb_terms_batch(ks, n, etas, r, False)
            require(h0 is None and np.array_equal(lp, lp0) and np.array_equal(g, g0),
                    'Score-only and Hessian calls disagree')
            for node, case in enumerate(group):
                column = ks.index(case['k'])
                ladder = case['reference_ladder']
                require([(z['digits'], z['step']) for z in ladder] ==
                        [(d, s) for d in (80, 110) for s in ('1e-15', '1e-20')], 'Reference ladder differs')
                with localcontext() as context:
                    context.prec = 110
                    flat = lambda z: [Decimal(z['logpmf']), *map(Decimal, z['score']),
                                      *(Decimal(v) for row in z['hessian'] for v in row)]
                    final = flat(ladder[-1])
                    stability = max(abs(a-b) for z in ladder for a, b in zip(flat(z), final))
                    values = [lp[node, column], *g[node, column], *h[node, column].ravel()]
                    differences = [Decimal.from_float(float(a))-b for a, b in zip(values, final)]
                    passed = (abs(differences[0]) <= Decimal('1e-9') and
                              max(map(abs, differences[1:3])) <= Decimal('1e-8') and
                              max(map(abs, differences[3:])) <= Decimal('1e-5') and
                              stability <= Decimal('1e-25'))
                result = dict(origin=case['origin'], n=n, k=case['k'], eta=case['eta'], r=r,
                    batch_group=index, batch_node=node, raw_candidate_minus_Decimal=list(map(str, differences)),
                    Decimal_stability=str(stability), selected_criterion_met=bool(passed))
                results.append(result)
            save('group_'+str(index), dict(n=n, r=r, etas=etas, counts=ks, logpmf=lp, score=g, Hessian=h))
        save('selected_comparisons', results)
        # Timing is component-level evidence only; the native scalar helper is
        # a comparator, and the21 independent reference points are kept separate.
        ks = [2466, 0, 822, 822]
        etas = np.linspace(-.3, .3, 32)
        timing = []
        for repeat in range(3):
            began = time.perf_counter()
            vector = batch.stable_bb_terms_batch(ks, 2466, etas, -15.)
            vector_seconds = time.perf_counter()-began
            began = time.perf_counter()
            rows = [scalar.stable_bb_terms(ks, 2466, float(eta), -15.) for eta in etas]
            scalar_seconds = time.perf_counter()-began
            serial = tuple(np.stack([row[k] for row in rows]) for k in range(3))
            differences = [float(np.max(np.abs(a-b))) for a, b in zip(vector, serial)]
            timing.append(dict(repeat=repeat, batch_seconds=vector_seconds, scalar_seconds=scalar_seconds,
                scalar_over_batch=scalar_seconds/vector_seconds, maximum_absolute_LP_score_H_difference=differences,
                nominal_comparison_within_1e_minus10=all(d <= 1e-10 for d in differences),
                independent_accuracy_or_fit_speed_claim=False))
        save('timing', dict(records=timing, memory_estimate=batch.memory_estimate(32, 2466, 4)))
        for path, digest in pins.items():
            require(sha(path) == digest, 'Source/input changed during execution')
        save('summary', dict(completed=True, selected_points=len(results),
            selected21_criterion_met=all(r['selected_criterion_met'] for r in results),
            batch_scalar_timing_comparison_within_1e_minus10=all(r['nominal_comparison_within_1e_minus10'] for r in timing),
            measured_process_peak_RSS_KiB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            no_integration=True, no_fit=True, no_calibration=True, no_adoption=True))
    except BaseException:
        with (output/'failure.json').open('x') as handle:
            json.dump(dict(traceback=traceback.format_exc(), partials_preserved=True, no_adoption=True), handle)
        raise


if __name__ == '__main__':
    main()
