"""Compute-only synthetic probe of exact log-shape beta-binomial arithmetic.

Does not edit/import live likelihood code, fit a model, or read fixture rows.
The comparison implementation is the immutable P3 executed-source snapshot.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import sys
import time

import numpy as np
import scipy
from scipy.special import expit, gammaln, logsumexp
from scipy.stats import betabinom


ROOT = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas')
REC = ROOT / 'Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z'
SNAPSHOT = REC / 'bfix_p3_21998519/executed_sources/likelihood_v1.py'
EPS = 1e-6


def rising_tables(n, logshape):
    """R_s, dR_s/dlogshape, d²R_s/dlogshape² for every integer s=0..n.

    The j=0 term is exactly logshape, with derivatives 1 and 0. No
    exp(logshape), reciprocal shape, or small-shape special function occurs.
    """
    if n == 0:
        return np.zeros(1), np.zeros(1), np.zeros(1)
    logs = np.empty(n)
    first = np.empty(n)
    second = np.empty(n)
    logs[0], first[0], second[0] = logshape, 1., 0.
    logj = np.log(np.arange(1, n, dtype=float))
    logs[1:] = np.logaddexp(logshape, logj)
    first[1:] = expit(logshape-logj)
    second[1:] = first[1:]*expit(logj-logshape)
    return tuple(np.r_[0., np.cumsum(x)] for x in (logs, first, second))


def stable_bb_terms(k, n, eta, linear_r):
    """Exact integer beta-binomial mass and eta/linear_r derivatives.

    This is an isolated candidate numerical representation, not a change to
    the production likelihood, optimizer, clipping, bounds, or test law.
    """
    k = np.atleast_1d(np.asarray(k))
    if (int(n) != n or n < 0 or np.any(~np.isfinite(k))
            or np.any(k != np.floor(k)) or np.any((k < 0) | (k > n))
            or not np.isfinite(eta) or not np.isfinite(linear_r)):
        raise ValueError('finite predictors and integer BB support required')
    n, k = int(n), k.astype(int)
    m, q = expit(eta), expit(-eta)
    v, vbar = expit(linear_r), expit(-linear_r)
    rho = EPS+(1-2*EPS)*v
    rp = (1-2*EPS)*v*vbar
    rpp = rp*(vbar-v)
    t = (1-rho)/rho
    tp = -rp/rho**2
    tpp = -rpp/rho**2+2*rp**2/rho**3
    logt = np.log(t)
    loga = logt-np.logaddexp(0., -eta)
    logb = logt-np.logaddexp(0., eta)
    ra, da, qa = rising_tables(n, loga)
    rb, db, qb = rising_tables(n, logb)
    rt, dt, qt = rising_tables(n, logt)
    da, db, dt = da[k], db[n-k], dt[n]
    qa, qb, qt = qa[k], qb[n-k], qt[n]
    choose = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)
    mass = choose+ra[k]+rb[n-k]-rt[n]
    u = tp/t
    w = tpp/t-u*u
    grad = np.column_stack((da*q-db*m, (da+db-dt)*u))
    hess = np.empty((len(k), 2, 2))
    hess[:, 0, 0] = qa*q*q+qb*m*m-(da+db)*m*q
    hess[:, 0, 1] = hess[:, 1, 0] = u*(qa*q-qb*m)
    hess[:, 1, 1] = u*u*(qa+qb-qt)+w*(da+db-dt)
    return mass, grad, hess


def conditional_terms(k, n, eta, linear_r):
    lp, g, h = stable_bb_terms(k, n, eta, linear_r)
    weights = np.exp(lp-logsumexp(lp))
    mean = weights@g
    norm_hess = np.einsum('i,ijk->jk', weights, h)
    norm_hess += np.einsum('i,ij,ik->jk', weights, g, g)-np.outer(mean, mean)
    return lp-logsumexp(lp), g-mean, h-norm_hess


def load_snapshot():
    spec = importlib.util.spec_from_file_location('codex_p3_frozen_bb_snapshot', SNAPSHOT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('probe must run on a compute allocation')
    started = time.monotonic()
    args.output.mkdir(exist_ok=False)
    source = Path(__file__).read_bytes()
    (args.output / Path(__file__).name).write_bytes(source)
    old = load_snapshot()
    maxima = {name: 0. for name in (
        'snapshot_mass', 'snapshot_gradient', 'snapshot_hessian', 'scipy_mass',
        'reflection_mass', 'reflection_gradient', 'reflection_hessian',
        'full_log_normalization', 'conditional_score_identity',
        'conditional_hessian_identity', 'fd_gradient', 'fd_hessian',
        'fd_conditional_gradient', 'fd_conditional_hessian')}
    cases = 0

    def close(name, actual, expected, atol, rtol=0.):
        error = float(np.max(np.abs(actual-expected)))
        maxima[name] = max(maxima[name], error)
        if not np.allclose(actual, expected, atol=atol, rtol=rtol):
            raise AssertionError(f'{name}: max absolute error {error}')

    # Moderate eta: all endpoint/interior counts, including the rho floor.
    with np.errstate(over='raise', divide='raise', invalid='raise'):
        for n in (20, 200):
            k = np.arange(n+1)
            for eta in (-6., 0., 6.):
                for r in (-15., -3., 5.):
                    new = stable_bb_terms(k, n, eta, r)
                    prior = old.bb_terms(k, n, eta, r)
                    for name, actual, expected, tolerance in zip(
                            ('snapshot_mass', 'snapshot_gradient', 'snapshot_hessian'),
                            new, prior, (1e-7, 1e-7, 1e-6)):
                        close(name, actual, expected, tolerance, 1e-7)
                    rho = EPS+(1-2*EPS)*expit(r)
                    t = (1-rho)/rho
                    reference = betabinom.logpmf(k, n, expit(eta)*t, expit(-eta)*t)
                    close('scipy_mass', new[0], reference, 1e-7, 1e-8)

        # Finite tails, full/truncated integer support, reflection and chain rule.
        # Percent-based inclusion is the actual integer predicate, not rounding.
        for n in (20, 200):
            full = np.arange(n+1)
            lo = max(3, (10*n+99)//100)
            selected = np.arange(lo, n-lo+1)
            for eta in (-1000., -40., -6., 0., 6., 40., 1000.):
                for r in (-1000., -15., -3., 5., 1000.):
                    values = stable_bb_terms(full, n, eta, r)
                    if not all(np.all(np.isfinite(a)) for a in values):
                        raise AssertionError('nonfinite BB value/gradient/Hessian')
                    close('full_log_normalization', np.array(logsumexp(values[0])), np.array(0.), 1e-8)
                    reflected = stable_bb_terms(n-full, n, -eta, r)
                    signs = np.array([-1., 1.])
                    close('reflection_mass', values[0], reflected[0], 1e-9)
                    close('reflection_gradient', values[1], reflected[1]*signs, 1e-9)
                    close('reflection_hessian', values[2], reflected[2]*np.outer(signs, signs), 1e-9)
                    for support, transform, prefix in (
                            (full, stable_bb_terms, 'fd_'),
                            (selected, conditional_terms, 'fd_conditional_')):
                        lp, g, h = transform(support, n, eta, r)
                        weights = np.exp(lp-logsumexp(lp))
                        if prefix == 'fd_conditional_':
                            close('conditional_score_identity', weights@g, np.zeros(2), 1e-9)
                            identity = np.einsum('i,ijk->jk', weights, h)
                            identity += np.einsum('i,ij,ik->jk', weights, g, g)
                            close('conditional_hessian_identity', identity, np.zeros((2, 2)), 1e-9)
                        step = 2e-4
                        fdg, fdh = np.empty_like(g), np.empty_like(h)
                        for axis in range(2):
                            plus, minus = np.array([eta, r]), np.array([eta, r])
                            plus[axis] += step
                            minus[axis] -= step
                            ap = transform(support, n, *plus)
                            am = transform(support, n, *minus)
                            fdg[:, axis] = (ap[0]-am[0])/(2*step)
                            fdh[:, :, axis] = (ap[1]-am[1])/(2*step)
                        close(prefix+'gradient', fdg, g, 2e-6, 2e-6)
                        close(prefix+'hessian', fdh, h, 2e-6, 2e-6)

                    # Inclusion-reweighted three-class CDF, with fixed native
                    # homozygous BB distributions, verifies generation arithmetic.
                    lp = stable_bb_terms(selected, n, eta, r)[0]
                    fixed_t = (1-.02)/.02
                    parts = np.stack((lp,
                        betabinom.logpmf(selected, n, .01*fixed_t, .99*fixed_t),
                        betabinom.logpmf(selected, n, .99*fixed_t, .01*fixed_t)))
                    classes = np.log(np.array([.35, .5, .15]))+logsumexp(parts, axis=1)
                    probabilities = np.exp(classes-logsumexp(classes))
                    distributions = np.exp(parts-logsumexp(parts, axis=1)[:, None])
                    for probabilities_one in (probabilities, *distributions):
                        cdf = np.cumsum(probabilities_one)
                        if (not np.all(np.isfinite(cdf)) or np.any(np.diff(cdf)<0)
                                or abs(cdf[-1]-1)>1e-10):
                            raise AssertionError('invalid class/count CDF')
                    cases += 1
    summary = {
        'passed': True, 'scope': 'synthetic BB arithmetic only; no fit or P3 size/power',
        'tail_grid_cases': cases, 'max_absolute_errors': maxima,
        'snapshot': str(SNAPSHOT), 'snapshot_sha256': hashlib.sha256(SNAPSHOT.read_bytes()).hexdigest(),
        'probe_sha256': hashlib.sha256(source).hexdigest(),
        'job_id': os.environ['SLURM_JOB_ID'], 'python': platform.python_version(),
        'numpy': np.__version__, 'scipy': scipy.__version__,
        'elapsed_seconds': time.monotonic()-started,
        'method': 'exact integer rising factorials in log-shape coordinates; no clipping',
        'finite_difference_step': 2e-4,
        'note': 'Arithmetic stability does not establish optimizer convergence or scientific calibration.'}
    (args.output / 'summary.json').write_text(json.dumps(summary, indent=2)+'\n')
    print(json.dumps({'passed': True, 'tail_grid_cases': cases, 'output': str(args.output)}))


if __name__ == '__main__':
    main()
