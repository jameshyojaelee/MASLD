"""Compute-only equivalence check of discarded-Hessian removal; no fits."""
import argparse
import ast
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('campaign', type=Path)
    parser.add_argument('expected_sha256')
    args = parser.parse_args()
    source = Path(__file__).read_bytes()
    assert hashlib.sha256(source).hexdigest() == args.expected_sha256
    manifest = json.loads((args.campaign / 'campaign.json').read_text())
    for name, expected in manifest['files'].items():
        assert hashlib.sha256(Path(name).read_bytes()).hexdigest() == expected, name
    archive = args.campaign / 'executed_sources'
    original = (archive / 'likelihood_v1.py').read_text()
    replacements = [
        ('def hstats(u):', 'def hstats(u, need_curvature=True):'),
        ('value,score,hh=conditional_terms(row,eta,r)\n                    ll+=value;s+=score[0]*row.stage;h+=hh[0,0]*row.stage**2',
         'value,score,hh=conditional_terms(row,eta,r,need_hessian=need_curvature)\n                    ll+=value\n                    if need_curvature:\n                        s+=score[0]*row.stage;h+=hh[0,0]*row.stage**2'),
        ('if hstats(trial)[0]>=ll-1e-12:', 'if hstats(trial, need_curvature=False)[0]>=ll-1e-12:'),
        ('lp=bb_terms(k,row.n,eta-row.orientation*row.omega,linear_r)[0]',
         'lp=bb_terms(k,row.n,eta-row.orientation*row.omega,linear_r,need_hessian=False)[0]')]
    candidate = original
    for old, new in replacements:
        assert candidate.count(old) == 1, old
        candidate = candidate.replace(old, new)
    ast.parse(candidate)
    job = os.environ['SLURM_JOB_ID']
    output = args.campaign.parent / ('p3_cost_check_' + job)
    output.mkdir(exist_ok=False)
    (output / 'executed_check.py').write_bytes(source)
    candidate_path = output / 'likelihood_cost_candidate.py'
    candidate_path.write_text(candidate)
    sys.path.insert(0, str(archive))
    import numpy as np
    import scipy
    import likelihood_v1 as reference
    from bfix_v1 import load, FIXTURE
    spec = importlib.util.spec_from_file_location('likelihood_cost_candidate', candidate_path)
    fast = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = fast
    spec.loader.exec_module(fast)
    design, settings = load()
    rows = design.model().rows
    old_model = reference.Model(rows)
    fast_model = fast.Model(rows)
    uniforms = np.load(FIXTURE / 'draws/p3_uniforms.npy')
    null = json.loads((args.campaign / 'baseline/observed_null.json').read_text())
    points = [('observed_null', null['parameters'])]
    for b in [18, 19]:
        fit = json.loads((args.campaign / 'baseline' / f'p3_{b:03d}.json').read_text())['fit']
        points.append((f'draw_{b}_parameters', fit['parameters']))
    maximum = 0.
    def compare(a, b):
        nonlocal maximum
        x, y = np.asarray(a, dtype=float), np.asarray(b, dtype=float)
        assert x.shape == y.shape
        assert np.array_equal(np.isnan(x), np.isnan(y))
        assert np.array_equal(np.isposinf(x), np.isposinf(y))
        assert np.array_equal(np.isneginf(x), np.isneginf(y))
        good = np.isfinite(x) & np.isfinite(y)
        if good.any():
            delta = float(np.max(np.abs(x[good] - y[good])))
            maximum = max(maximum, delta)
            assert delta <= 1e-12, delta
    evaluations = []
    times = {'reference': 0., 'candidate': 0.}
    for label, p in points:
        for variance in [0., .00711885901713953, settings['tau_max_fixture'] ** 2]:
            for hessian in [False, True]:
                outcomes = []
                for name, model in [('reference', old_model), ('candidate', fast_model)]:
                    started = time.monotonic()
                    try:
                        value = model.marginal(np.asarray(p), variance, with_hessian=hessian)
                        outcome = ('value', value)
                    except (ValueError, FloatingPointError, np.linalg.LinAlgError) as error:
                        outcome = ('exception', type(error).__name__, str(error))
                    times[name] += time.monotonic() - started
                    outcomes.append(outcome)
                assert outcomes[0][0] == outcomes[1][0]
                if outcomes[0][0] == 'exception':
                    assert outcomes[0] == outcomes[1]
                else:
                    for a, b in zip(outcomes[0][1], outcomes[1][1]):
                        compare(a, b)
                record = {'point': label, 'variance': variance, 'hessian': hessian,
                          'outcome': outcomes[0][0], 'matched': True}
                evaluations.append(record)
                print(json.dumps(record), flush=True)
    sampling_checks = 0
    p = np.asarray(null['parameters'])
    for b in [0, 18, 19]:
        for row, (u1, u2) in zip(rows, uniforms[b]):
            eta, r, _, _ = old_model.predictors(row, p)
            a = reference.conditional_draw(row, eta, r, u1, u2)
            z = fast.conditional_draw(row, eta, r, u1, u2)
            assert a == z
            sampling_checks += 1
    result = {'status': 'pointwise_values_derivatives_failures_and_sampling_match',
              'absolute_tolerance': 1e-12, 'maximum_finite_difference': maximum,
              'evaluations': evaluations, 'sampling_checks': sampling_checks,
              'elapsed_seconds_by_version': times,
              'original_source_sha256': hashlib.sha256(original.encode()).hexdigest(),
              'candidate_source_sha256': hashlib.sha256(candidate.encode()).hexdigest(),
              'python': platform.python_version(), 'numpy': np.__version__, 'scipy': scipy.__version__,
              'fit_performed': False, 'active_code_changed': False,
              'limit': 'Pointwise equivalence is not fitted-result equivalence, calibration or measured end-to-end speedup.'}
    with (output / 'summary.json').open('x') as handle:
        json.dump(result, handle, indent=2)
    print(json.dumps(result), flush=True)


if __name__ == '__main__':
    main()
