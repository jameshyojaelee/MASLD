"""Isolated candidate numerical revision of the original selected likelihood.

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


def integrate_gene(model, p, variance, rows, conditional, with_hessian=False, settings=Settings()):
    """Normalized integral and exact-integral score/Hessian moment construction.

    Modes set numerical breakpoints/scaling only. All grid intervals are integrated;
    a lower local mode never replaces the whole integral. Domain/refinement/cap
    failures are explicit. Grid/root search is not proof of all stationary points.
    """
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Numerical candidate requires compute')
    require(np.isfinite(variance) and variance > 0 and np.isfinite(p).all() and rows, 'Positive finite variance/parameters required')
    width = model.size + 1
    cp, hp = derivative_bounds(model, p, rows)
    sqrtv = np.sqrt(variance)
    cache, evaluations = {}, 0

    def stats(z):
        nonlocal evaluations
        z = float(z)
        if z in cache:
            return cache[z]
        evaluations += 1
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

    def locate(bound):
        grid = np.linspace(-bound, bound, int(np.ceil(2*bound/settings.grid_spacing))+1)
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

    def run(bound, precision):
        center, points, peaks = locate(bound)
        shift, gs, hs, _, _ = stats(center)
        gradient_scale = np.maximum(1., np.abs(gs))
        second_scale = np.maximum(1., np.abs(hs+np.outer(gs, gs))) if with_hessian else None
        scale = np.r_[1., gradient_scale, second_scale.ravel()] if with_hessian else np.r_[1., gradient_scale]
        def function(z):
            h, g, hh, _, _ = stats(z)
            value = np.r_[1., g, (hh+np.outer(g, g)).ravel()] if with_hessian else np.r_[1., g]
            return np.exp(h-shift) * value / scale
        vector, error, info = quad_vec(function, -bound, bound, points=points, norm='max',
              quadrature='gk21', limit=settings.maximum_intervals, epsabs=settings.quad_absolute*precision,
              epsrel=settings.quad_relative*precision, full_output=True)
        require(info.success and np.isfinite(vector).all() and vector[0] > error >= 0,
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
        tail = tail_moments(bound, variance, lower, gradient, second, cp, hp, with_hessian)
        ge += tail['relative_gradient_tail_allowance']
        if with_hessian:
            he += tail['relative_hessian_tail_allowance']
        report = dict(z_bound=bound, u_bound=bound*sqrtv, integration_center_z=center, stationary_points=peaks,
                      stationary_points_exhaustive=False, all_grid_intervals_integrated=True,
                      quadrature_points=len(points), quadrature_status=int(info.status), quadrature_message=info.message,
                      quadrature_evaluations=int(info.neval), error_estimate=float(error),
                      scaled_mass=float(mass), component_scaling=scale,
                      log_bounded_mass_lower_estimate=float(lower),
                      logmass_error_allowance=float(-np.log1p(-error/mass)+np.log1p(tail['relative_mass_tail_allowance'])),
                      gradient_error_allowance=ge, hessian_error_allowance=he, **tail)
        return logmass, gradient, hessian, report

    bound = settings.initial_z_bound
    accepted = None
    for attempt in range(5):
        result = run(bound, 10.**-attempt)
        report = result[3]
        if report['relative_mass_tail_allowance'] > settings.logmass_tolerance / 8 or np.max(report['relative_gradient_tail_allowance']) > settings.gradient_tolerance / 8 or (
                with_hessian and np.max(report['relative_hessian_tail_allowance']) > settings.hessian_tolerance / 8):
            lower = report['log_bounded_mass_lower_estimate']
            requested = lower + np.log(min(settings.logmass_tolerance, settings.gradient_tolerance, settings.hessian_tolerance))
            requested -= np.log(16*(1+np.max(cp)+np.max(hp)+np.max(cp)**2))
            require(np.log(2.)+norm.logsf(settings.maximum_z_bound) < requested, 'Relative derivative-tail domain cap exhausted')
            desired = brentq(lambda z: np.log(2.)+norm.logsf(z)-requested, 0., settings.maximum_z_bound)
            bound = max(bound+2, float(np.ceil(desired)))
            require(bound+2 <= settings.maximum_z_bound, 'Integral domain refinement exceeds cap')
            continue
        if report['logmass_error_allowance'] <= settings.logmass_tolerance / 4 and np.max(report['gradient_error_allowance']) <= settings.gradient_tolerance / 4 and (
                not with_hessian or np.max(report['hessian_error_allowance']) <= settings.hessian_tolerance / 4):
            accepted = result
            break
    require(accepted is not None, 'Candidate numerical uncertainty remained above fixed precision allowance')
    refined = run(bound+2, 0.01)
    delta_ll = abs(refined[0]-accepted[0])
    delta_g = np.abs(refined[1]-accepted[1])
    delta_h = np.abs(refined[2]-accepted[2]) if with_hessian else None
    require(delta_ll+refined[3]['logmass_error_allowance'] <= settings.logmass_tolerance and
            np.max(delta_g+refined[3]['gradient_error_allowance']) <= settings.gradient_tolerance and
            (not with_hessian or np.max(delta_h+refined[3]['hessian_error_allowance']) <= settings.hessian_tolerance),
            'Expanded tighter integration did not agree within numerical allowances')
    report = dict(method='normalized-prior composite Gauss-Kronrod, all located modes and every grid interval',
                  numerical_revision_of_prescribed_20_node_approximation=True, same_probability_estimand=True,
                  first=accepted[3], refined=refined[3], logmass_refinement_difference=delta_ll,
                  gradient_refinement_difference=delta_g, hessian_refinement_difference=delta_h,
                  unique_likelihood_evaluations=evaluations, quad_error_estimates_not_certified_bounds=True,
                  floating_arithmetic_error_not_rigorously_bounded=True)
    return refined[0], refined[1], refined[2], report


def candidate_class(private_original, settings=Settings()):
    """Use only an independently loaded private copy of the hash-guarded source.

    The BB substitution reproduces the reviewed exact log-shape representation.
    No original module file is modified. Null/boundary/fit/statistical methods are
    inherited without source edits; failures remain ValueError subclasses.
    """
    require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Candidate construction requires compute')
    private_original.bb_terms = stable_bb_terms
    class Candidate(private_original.Model):
        def marginal(self, p, v, nodes=20, with_hessian=False, gene_terms=None):
            # nodes remains an inherited-call compatibility argument only.
            if v <= 0:
                self.last_integration_reports = []
                return super().marginal(p, v, nodes, with_hessian, gene_terms)
            total = 0.
            gradient = np.zeros(self.size+1)
            hessian = np.zeros((self.size+1, self.size+1)) if with_hessian else None
            reports = []
            present = sum(any(r.gene == gene for r in self.rows) for gene in range(self.G))
            gene_settings = replace(settings, logmass_tolerance=min(settings.logmass_tolerance, 1e-8/max(1, present)),
                                    gradient_tolerance=min(settings.gradient_tolerance, 1e-7/max(1, present)))
            for gene in range(self.G):
                rows = [r for r in self.rows if r.gene == gene]
                if not rows:
                    if gene_terms is not None:
                        gene_terms[gene] = 0.
                    continue
                ll, g, h, report = integrate_gene(self, p, v, rows, private_original.conditional_terms, with_hessian, gene_settings)
                total += ll
                gradient += g
                if with_hessian:
                    hessian += h
                if gene_terms is not None:
                    gene_terms[gene] = ll
                reports.append(dict(gene=gene, **report))
            self.last_integration_reports = reports
            self.last_total_error_allowance = dict(
                logmass=sum(r['refined']['logmass_error_allowance']+r['logmass_refinement_difference'] for r in reports),
                gradient=sum((r['refined']['gradient_error_allowance']+r['gradient_refinement_difference'] for r in reports),
                             np.zeros(self.size+1)))
            require(self.last_total_error_allowance['logmass'] <= 1e-8 and
                    np.max(self.last_total_error_allowance['gradient']) <= 1e-7,
                    'Accumulated gene numerical uncertainty exceeds full-model allowance')
            return (float(total), gradient, hessian) if with_hessian else (float(total), gradient)
    return Candidate
