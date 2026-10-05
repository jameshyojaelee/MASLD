"""Composite Gauss-Legendre diagnostic of the same normalized moments.

Integrates the normalized Gaussian-prior likelihood and derivative moments over
all located modes. This replaces the numerical 20-node rule, not the probability
law. Original outer optimizers, starts, bounds, tie and failure rules are inherited.
SciPy errors below are estimates, not certified bounds on floating arithmetic.
"""
from dataclasses import dataclass, replace
import os

import numpy as np
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


def integrate_gene(model, p, variance, rows, conditional, with_hessian=False, settings=Settings(), observer=None):
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

    geometry_cache = {}
    completed_results = []
    def run(bound, order, spacing, phase):
        nonlocal current_run, run_number
        run_number += 1
        current_run = dict(run_number=run_number,z_bound=bound,order=order,grid_spacing=spacing,phase=phase)
        emit('run_entered',with_hessian=with_hessian,diagnostic_only=True)
        key=(bound,spacing)
        if key not in geometry_cache:
            try:
                geometry_cache[key]=locate(bound,spacing)
            except BaseException as error:
                emit('locate_refusal',exception_type=type(error).__name__,reason=str(error))
                raise
        center,points,peaks=geometry_cache[key]
        emit('locate_completed',integration_center_z=center,quadrature_points=len(points),stationary_points=peaks,
             reused_geometry_across_orders=key==(15.,settings.grid_spacing))
        shift,gs,hs,_,_=stats(center)
        gradient_scale=np.maximum(1.,np.abs(gs))
        second_scale=np.maximum(1.,np.abs(hs+np.outer(gs,gs))) if with_hessian else None
        scale=np.r_[1.,gradient_scale,second_scale.ravel()] if with_hessian else np.r_[1.,gradient_scale]
        nodes,weights=np.polynomial.legendre.leggauss(order)
        require(np.all(weights>0) and abs(float(weights.sum())-2.)<1e-13,'GL node-weight normalization invalid')
        edges=np.asarray([-bound,*points,bound])
        require(np.all(np.diff(edges)>0),'GL interval geometry invalid')
        total=np.zeros(len(scale))
        compensation=np.zeros(len(scale))
        integrated_weight=0.
        try:
            for left,right in zip(edges[:-1],edges[1:]):
                half=(right-left)/2
                midpoint=(right+left)/2
                mapped=midpoint+half*nodes
                require(np.all(mapped>left) and np.all(mapped<right),'Mapped GL node outside interval')
                block=np.empty((order,len(scale)))
                for i,z in enumerate(mapped):
                    h,g,hh,_,_=stats(z)
                    value=np.r_[1.,g,(hh+np.outer(g,g)).ravel()] if with_hessian else np.r_[1.,g]
                    block[i]=half*weights[i]*np.exp(h-shift)*value/scale
                require(np.isfinite(block).all(),'Nonfinite GL integrand')
                # Pairwise node reduction, compensated interval accumulation.
                increment=np.sum(block,axis=0)-compensation
                updated=total+increment
                compensation=(updated-total)-increment
                total=updated
                integrated_weight+=half*float(weights.sum())
        except BaseException as error:
            emit('fixed_node_refusal',exception_type=type(error).__name__,reason=str(error),
                 final_integral_metrics_available=False)
            raise
        require(abs(integrated_weight-2*bound)<1e-10,'GL interval Jacobian normalization invalid')
        mass=float(total[0])
        require(np.isfinite(total).all() and mass>0,'Invalid GL integrated vector')
        values=total*scale
        gradient=values[1:1+width]/mass
        second=values[1+width:].reshape(width,width)/mass if with_hessian else None
        hessian=second-np.outer(gradient,gradient) if with_hessian else None
        logmass=shift+np.log(mass)-.5*np.log(2*np.pi)
        # No certified lower mass is available: retain the existing tail formula,
        # conditioning it explicitly on this numerical mass estimate.
        tail=tail_moments(bound,variance,logmass,gradient,second,cp,hp,with_hessian)
        report=dict(z_bound=bound,u_bound=bound*sqrtv,order=order,grid_spacing=spacing,
            integration_center_z=center,log_integrand_shift=shift,stationary_points=peaks,
            stationary_points_exhaustive=False,all_grid_intervals_integrated=True,
            intervals=len(edges)-1,fixed_nodes=(len(edges)-1)*order,
            unique_likelihood_evaluations=evaluations,scaled_mass=mass,component_scaling=scale,
            normalized_prior_constant=-.5*np.log(2*np.pi),integrated_unweighted_domain_length=integrated_weight,
            quadrature_error_estimate_available=False,mass_estimate_not_certified_lower_bound=True,
            floating_arithmetic_error_not_rigorously_bounded=True,
            diagnostic_only=True,accepted_candidate_integral=False,**tail)
        result=(float(logmass),gradient,hessian,report)
        completed_results.append(result)
        emit('order_result',logmass=result[0],gradient=gradient,hessian=hessian,report=report)
        return result

    orders=[run(15.,order,settings.grid_spacing,'order_ladder') for order in (16,32,64)]
    refined=run(17.,64,settings.grid_spacing/2,'expanded_half_grid_refinement')
    comparisons=[]
    for name,first,last,fraction in (
            ('order16_to32',orders[0],orders[1],.25),
            ('order32_to64',orders[1],orders[2],.25),
            ('order64_to_expanded_half_grid64',orders[2],refined,1.)):
        delta_ll=abs(last[0]-first[0])
        delta_g=np.abs(last[1]-first[1])
        delta_h=np.abs(last[2]-first[2]) if with_hessian else None
        mass_tail=last[3]['relative_mass_tail_allowance']
        ll_allowance=delta_ll+np.log1p(mass_tail)
        g_allowance=delta_g+last[3]['relative_gradient_tail_allowance']
        h_allowance=delta_h+last[3]['relative_hessian_tail_allowance'] if with_hessian else None
        budgets=dict(logmass=settings.logmass_tolerance*fraction,
                     gradient=settings.gradient_tolerance*fraction,hessian=settings.hessian_tolerance*fraction)
        comparison=dict(name=name,logmass_difference=delta_ll,gradient_difference=delta_g,hessian_difference=delta_h,
            logmass_difference_plus_tail=ll_allowance,gradient_difference_plus_tail=g_allowance,
            hessian_difference_plus_tail=h_allowance,budgets=budgets,
            descriptive_allowances_within_budgets=bool(ll_allowance<=budgets['logmass'] and
                np.max(g_allowance)<=budgets['gradient'] and
                (not with_hessian or np.max(h_allowance)<=budgets['hessian'])),
            differences_are_uncertified_estimates=True,unlocated_narrow_peak_may_be_missed=True)
        comparisons.append(comparison)
        emit('refinement_comparison',comparison=comparison)
    report=dict(method='Composite fixed Gauss-Legendre on all retained intervals',
        diagnostic_only=True,accepted_candidate_integral=False,no_optimizer_or_calibration=True,
        same_probability_estimand=True,orders=(16,32,64),base_z_bound=15.,refined_z_bound=17.,
        initial_grid_spacing=settings.grid_spacing,refined_grid_spacing=settings.grid_spacing/2,
        numerical_targets_and_evaluation_cap_unchanged=True,comparisons=comparisons,
        unique_likelihood_evaluations=evaluations,uncertified_refinement_estimates=True,
        unlocated_narrow_peak_may_be_missed=True,floating_arithmetic_error_not_rigorously_bounded=True)
    emit('composite_diagnostic_completed',report=report)
    return completed_results,report

