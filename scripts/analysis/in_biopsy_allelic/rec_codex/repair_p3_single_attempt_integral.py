"""Prespecified B15 single-run diagnostic; not an accepted integral candidate.

Integrates the normalized Gaussian-prior likelihood and derivative moments over
all located modes. This replaces the numerical 20-node rule, not the probability
law. Original outer optimizers, starts, bounds, tie and failure rules are inherited.
SciPy errors below are estimates, not certified bounds on floating arithmetic.
"""
from dataclasses import dataclass, replace
import os

import numpy as np
from scipy.integrate import quad_vec
from scipy.optimize import brentq, minimize_scalar
from scipy.special import expit, gammaln, logsumexp
from scipy.stats import norm


class IntegrationFailure(ValueError):
    """Preserves the original ValueError fitting-failure disposition."""


def require(ok, message):
    if not ok:
        raise IntegrationFailure(message)


@dataclass(frozen=True)
class Settings:
    logmass_tolerance: float = 1e-9
    gradient_tolerance: float = 1e-8
    hessian_tolerance: float = 1e-5
    quad_relative: float = 2e-12
    quad_absolute: float = 2e-12
    initial_z_bound: float = 8.0
    maximum_z_bound: float = 96.0
    grid_spacing: float = 0.25
    maximum_evaluations: int = 100000
    maximum_intervals: int = 6000


def rising_tables(n, logshape):
    """Exact integer rising factorials in log-shape coordinates; no shape floor."""
    if n == 0:
        return np.zeros(1), np.zeros(1), np.zeros(1)
    j = np.log(np.arange(1, n, dtype=float))
    value = np.r_[logshape, np.logaddexp(logshape, j)]
    score = np.r_[1., expit(logshape - j)]
    curvature = np.r_[0., expit(logshape - j) * expit(j - logshape)]
    return tuple(np.r_[0., np.cumsum(a)] for a in (value, score, curvature))


def stable_bb_terms(k, n, eta, linear_r, need_hessian=True):
    """Same exact BB law/derivatives as reviewed codex_stable_bb_probe.py.

    Derived from log(a)_k+log(b)_(n-k)-log(a+b)_n+log choose(n,k).
    Uses original rho=1e-6+(1-2e-6)expit(r), including its chain derivatives.
    """
    k = np.atleast_1d(np.asarray(k))
    require(int(n) == n and n >= 0 and np.isfinite(k).all() and np.all(k == np.floor(k)) and
            np.all((k >= 0) & (k <= n)) and np.isfinite(eta) and np.isfinite(linear_r), 'Invalid integer BB argument')
    n, k = int(n), k.astype(int)
    m, q, vr, qr = expit(eta), expit(-eta), expit(linear_r), expit(-linear_r)
    rho = 1e-6 + (1 - 2e-6) * vr
    rp = (1 - 2e-6) * vr * qr
    rpp = rp * (qr - vr)
    t = (1 - rho) / rho
    tp = -rp / rho**2
    tpp = -rpp / rho**2 + 2 * rp**2 / rho**3
    logt = np.log(t)
    ra, da, qa = rising_tables(n, logt - np.logaddexp(0., -eta))
    rb, db, qb = rising_tables(n, logt - np.logaddexp(0., eta))
    rt, dt, qt = rising_tables(n, logt)
    da, db, dt, qa, qb, qt = da[k], db[n-k], dt[n], qa[k], qb[n-k], qt[n]
    mass = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1) + ra[k] + rb[n-k] - rt[n]
    u, w = tp / t, tpp / t - (tp / t)**2
    score = np.column_stack((da * q - db * m, (da + db - dt) * u))
    if not need_hessian:
        return mass, score, None
    hess = np.empty((len(k), 2, 2))
    hess[:, 0, 0] = qa*q*q + qb*m*m - (da+db)*m*q
    hess[:, 0, 1] = hess[:, 1, 0] = u * (qa*q - qb*m)
    hess[:, 1, 1] = u*u*(qa+qb-qt) + w*(da+db-dt)
    return mass, score, hess


def derivative_bounds(model, p, rows):
    """Uniform-u score and conservative conditional-Hessian bounds.

    Each BB score <=n, BB Hessian entries <=4n (|dlogt/dr|<=1,
    |d2logt/dr2|<=2). Mixture/selection covariance adds bounded score
    products. 16n+8n^2 safely bounds each selected row Hessian entry.
    Predictor Jacobians and eta Hessians are independent of u here.
    """
    c = np.zeros(model.size)
    hh = np.zeros((model.size, model.size))
    for row in rows:
        _, _, jac, eta_hess = model.predictors(row, p, 0., need_hessian=True)
        j = np.abs(jac).sum(axis=0)
        c += 2 * row.n * j
        hh += (16*row.n + 8*row.n**2) * np.outer(j, j) + 2*row.n*np.abs(eta_hess)
    return c, hh


def tail_moments(bound, variance, logmass_lower, gradient, second, cp, hp, with_hessian):
    """Relative omitted mass/score/Hessian moment allowances, conditional on
    the numerical lower-mass estimate. Uses all-mode Gaussian envelope L<=1.
    z=u/sqrt(v); M2=2(a*phi(a)+sf(a)), M4=2((a^3+3a)phi(a)+3sf(a)).
    """
    logsf, logpdf, lv = norm.logsf(bound), norm.logpdf(bound), np.log(variance)
    logm0 = np.log(2.) + logsf
    logm2 = np.log(2.) + logsumexp([np.log(bound)+logpdf, logsf])
    logm4 = np.log(2.) + logsumexp([np.log(bound**3+3*bound)+logpdf, np.log(3.)+logsf])
    ratio = float(np.exp(logm0-logmass_lower))
    # Absolute variance score tail: (M2+M0)/(2v).
    vr = float(np.exp(logsumexp([logm2, logm0])-np.log(2.)-lv-logmass_lower))
    numerator = np.r_[cp*ratio, vr]
    gradient_error = numerator + np.abs(gradient)*ratio
    hessian_error = None
    if with_hessian:
        moment_tail = np.empty_like(second)
        moment_tail[:-1, :-1] = (hp + np.outer(cp, cp)) * ratio
        moment_tail[:-1, -1] = moment_tail[-1, :-1] = cp * vr
        # |Hvv+sv^2| <= (.75+1.5*z^2+.25*z^4)/v^2.
        moment_tail[-1, -1] = np.exp(logsumexp([np.log(.75)+logm0, np.log(1.5)+logm2,
                                               np.log(.25)+logm4])-2*lv-logmass_lower)
        hessian_error = moment_tail + np.abs(second)*ratio
        hessian_error += np.outer(np.abs(gradient), gradient_error) + np.outer(gradient_error, np.abs(gradient))
        hessian_error += np.outer(gradient_error, gradient_error)
    return dict(log_absolute_prior_mass_tail_bound=float(logm0), relative_mass_tail_allowance=ratio,
                relative_gradient_tail_allowance=gradient_error, relative_hessian_tail_allowance=hessian_error,
                error_allowances_conditional_on_estimated_lower_mass=True)


def integrate_gene(model, p, variance, rows, conditional, with_hessian=False, settings=Settings(), observer=None, diagnostic_precision=1.0):
    """Normalized integral and exact-integral score/Hessian moment construction.

    Modes set numerical breakpoints/scaling only. All grid intervals are integrated;
    a lower local mode never replaces the whole integral. Domain/refinement/cap
    failures are explicit. Grid/root search is not proof of all stationary points.
    """
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Numerical candidate requires compute')
    require(np.isfinite(variance) and variance > 0 and np.isfinite(p).all() and rows, 'Positive finite variance/parameters required')
    require(diagnostic_precision in (1.0, 0.5), 'Prespecified diagnostic precision required')
    width = model.size + 1
    cp, hp = derivative_bounds(model, p, rows)
    sqrtv = np.sqrt(variance)
    cache, evaluations = {}, 0
    current_run = {}
    run_number = 0
    attempt_context = {}

    def emit(event, **fields):
        if observer is not None:
            observer(dict(event=event, unique_evaluations=evaluations, maximum_evaluations=settings.maximum_evaluations,
                          **current_run, **fields))

    def stats(z):
        nonlocal evaluations
        z = float(z)
        if z in cache:
            return cache[z]
        evaluations += 1
        if evaluations > settings.maximum_evaluations:
            emit('evaluation_cap_refusal', quadrature_final_metrics_available=False)
        elif evaluations % 10000 == 0:
            emit('evaluation_progress', quadrature_final_metrics_available=False)
        require(evaluations <= settings.maximum_evaluations, 'Candidate integrand evaluation cap exhausted')
        u = sqrtv * z
        ll, su, hu = 0., -z, -1.
        gp = np.zeros(model.size)
        hh = np.zeros((width, width)) if with_hessian else None
        for row in rows:
            eta, dispersion, jac, eta_hess = model.predictors(row, p, u, need_hessian=with_hessian)
            value, score, row_hess = conditional(row, eta, dispersion, need_hessian=True)
            require(np.isfinite(value) and value <= 1e-8 and np.isfinite(score).all() and np.isfinite(row_hess).all(),
                    'Conditional likelihood/derivatives invalid')
            ll += value
            gp += jac.T @ score
            su += sqrtv * row.stage * score[0]
            hu += variance * row.stage**2 * row_hess[0, 0]
            if with_hessian:
                hh[:-1, :-1] += jac.T @ row_hess @ jac + score[0]*eta_hess
        # Prior scores differentiated at fixed u, not at moving quadrature nodes.
        gradient = np.r_[gp, (z*z-1)/(2*variance)]
        if with_hessian:
            hh[-1, -1] = (.5-z*z)/variance**2
        result = (ll-z*z/2, gradient, hh, su, hu)
        require(all(np.isfinite(a).all() for a in result if a is not None), 'Nonfinite normalized-integral statistics')
        cache[z] = result
        return result

    def locate(bound, spacing):
        grid = np.linspace(-bound, bound, int(np.ceil(2*bound/spacing))+1)
        values = [stats(z) for z in grid]
        candidates = [float(grid[int(np.argmax([a[0] for a in values]))])]
        stationary = []
        for left, right, sl, sr in zip(grid[:-1], grid[1:], values[:-1], values[1:]):
            if sl[3] * sr[3] < 0:
                root = brentq(lambda z: stats(z)[3], left, right, xtol=1e-11)
                stationary.append(root)
                candidates.append(root)
        # Explicitly bounded local searches cover sampled peaks with no score bracket.
        for i in range(1, len(grid)-1):
            if values[i][0] >= max(values[i-1][0], values[i+1][0]):
                peak = minimize_scalar(lambda z: -stats(z)[0], bounds=(grid[i-1], grid[i+1]), method='bounded',
                                       options={'xatol': 1e-10})
                require(peak.success, 'Bounded integration-center search failed')
                candidates.append(float(peak.x))
        center = max(candidates, key=lambda z: stats(z)[0])
        points = set(map(float, grid[1:-1]))
        peaks = []
        for root in stationary:
            ll, _, _, score, curvature = stats(root)
            peaks.append(dict(z=root, u=sqrtv*root, log_integrand_without_normalizer=ll,
                              z_score=score, z_curvature=curvature))
            if curvature < 0:
                scale = (-curvature)**-.5
                points.update(root + factor*scale for factor in (-8., -4., -1., 0., 1., 4., 8.) if -bound < root+factor*scale < bound)
        return center, sorted(points), peaks

    def run(bound, precision, spacing):
        nonlocal current_run, run_number
        run_number += 1
        current_run = dict(run_number=run_number, z_bound=bound, precision=precision, grid_spacing=spacing, **attempt_context)
        emit('run_entered', requested_epsabs=settings.quad_absolute*precision,
             requested_epsrel=settings.quad_relative*precision, with_hessian=with_hessian)
        try:
            center, points, peaks = locate(bound, spacing)
        except BaseException as error:
            emit('locate_refusal', exception_type=type(error).__name__, reason=str(error),
                 quadrature_final_metrics_available=False)
            raise
        emit('locate_completed', integration_center_z=center, quadrature_points=len(points), stationary_points=peaks)
        shift, gs, hs, _, _ = stats(center)
        gradient_scale = np.maximum(1., np.abs(gs))
        second_scale = np.maximum(1., np.abs(hs+np.outer(gs, gs))) if with_hessian else None
        scale = np.r_[1., gradient_scale, second_scale.ravel()] if with_hessian else np.r_[1., gradient_scale]
        def function(z):
            h, g, hh, _, _ = stats(z)
            value = np.r_[1., g, (hh+np.outer(g, g)).ravel()] if with_hessian else np.r_[1., g]
            return np.exp(h-shift) * value / scale
        try:
            vector, error, info = quad_vec(function, -bound, bound, points=points, norm='max',
                  quadrature='gk21', limit=settings.maximum_intervals, epsabs=settings.quad_absolute*precision,
                  epsrel=settings.quad_relative*precision, full_output=True)
        except BaseException as error:
            emit('quadrature_exception', exception_type=type(error).__name__, reason=str(error),
                 quadrature_final_metrics_available=False)
            raise
        emit('quadrature_returned', scaled_mass=float(vector[0]), error_estimate=float(error),
             quadrature_status=int(info.status), quadrature_message=info.message,
             quadrature_success=bool(info.success), quadrature_evaluations=int(info.neval),
             component_scaling=scale, quadrature_final_metrics_available=True,
             scaled_vector_maxnorm=float(np.max(np.abs(vector))),
             implied_maxnorm_tolerance=max(settings.quad_absolute*precision,
                                          settings.quad_relative*precision*float(np.max(np.abs(vector)))),
             scipy_global_error_target_divisor=8.)
        require(np.isfinite(vector).all() and vector[0] > error >= 0,
                'Composite integral failed: ' + info.message)
        mass = vector[0]
        value = vector * scale
        gradient = value[1:1+width] / mass
        logmass = shift + np.log(mass) - .5*np.log(2*np.pi)
        lower = shift + np.log(mass-error) - .5*np.log(2*np.pi)
        ge = error * (gradient_scale + np.abs(gradient)) / (mass-error)
        second = value[1+width:].reshape(width, width) / mass if with_hessian else None
        hessian = second - np.outer(gradient, gradient) if with_hessian else None
        he = None
        if with_hessian:
            he = error*(second_scale+np.abs(second))/(mass-error)
            he += np.outer(np.abs(gradient), ge) + np.outer(ge, np.abs(gradient)) + np.outer(ge, ge)
        quadrature_gradient_error = ge.copy()
        quadrature_hessian_error = he.copy() if with_hessian else None
        tail = tail_moments(bound, variance, lower, gradient, second, cp, hp, with_hessian)
        ge += tail['relative_gradient_tail_allowance']
        if with_hessian:
            he += tail['relative_hessian_tail_allowance']
        report = dict(z_bound=bound, u_bound=bound*sqrtv, integration_center_z=center, stationary_points=peaks,
                      requested_grid_spacing=spacing,
                      stationary_points_exhaustive=False, all_grid_intervals_integrated=True,
                      quadrature_points=len(points), quadrature_status=int(info.status), quadrature_message=info.message,
                      quadrature_evaluations=int(info.neval), error_estimate=float(error),
                      scaled_mass=float(mass), component_scaling=scale,
                      log_bounded_mass_lower_estimate=float(lower),
                      logmass_error_allowance=float(-np.log1p(-error/mass)+np.log1p(tail['relative_mass_tail_allowance'])),
                      gradient_error_allowance=ge, hessian_error_allowance=he,
                      quadrature_only_gradient_error_allowance=quadrature_gradient_error,
                      quadrature_only_hessian_error_allowance=quadrature_hessian_error,
                      quadrature_only_logmass_error_allowance=float(-np.log1p(-error/mass)), **tail)
        emit('run_numerical_allowances', logmass=logmass, gradient=gradient, hessian=hessian,
             report=report, quadrature_success=bool(info.success))
        require(info.success, 'Composite integral failed: ' + info.message)
        return logmass, gradient, hessian, report

    # One run only. No acceptance, retry, domain extension, or final refinement.
    attempt_context = dict(phase='prespecified_single_run_diagnostic', attempt=0)
    result = run(15., diagnostic_precision, settings.grid_spacing)
    report = result[3]
    budget_comparison = dict(
        logmass_budget=settings.logmass_tolerance/4,
        gradient_budget=settings.gradient_tolerance/4,
        hessian_budget=settings.hessian_tolerance/4,
        logmass_allowance_within_budget=report['logmass_error_allowance']<=settings.logmass_tolerance/4,
        gradient_allowance_within_budget=bool(np.max(report['gradient_error_allowance'])<=settings.gradient_tolerance/4),
        hessian_allowance_within_budget=bool(np.max(report['hessian_error_allowance'])<=settings.hessian_tolerance/4) if with_hessian else None,
        quadrature_only_gradient_max=float(np.max(report['quadrature_only_gradient_error_allowance'])),
        total_gradient_max=float(np.max(report['gradient_error_allowance'])),
        quadrature_only_hessian_max=float(np.max(report['quadrature_only_hessian_error_allowance'])) if with_hessian else None,
        total_hessian_max=float(np.max(report['hessian_error_allowance'])) if with_hessian else None)
    report.update(diagnostic_only=True, accepted_candidate_integral=False, full_refinement_completed=False,
        same_probability_estimand=True, diagnostic_precision=diagnostic_precision,
        unique_likelihood_evaluations=evaluations, budget_comparison=budget_comparison,
        quad_error_estimates_not_certified_bounds=True, floating_arithmetic_error_not_rigorously_bounded=True)
    emit('single_run_diagnostic_completed', report=report)
    return result

