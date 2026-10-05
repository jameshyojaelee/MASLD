#!/usr/bin/env python3
"""Fixed-point composite GL numerical revision; no fits or calibration.

Every finite-grid interval is integrated. Order/grid/domain agreement is an
empirical precision check, not a certified quadrature or floating-error bound.
Independent scalar and finite-difference comparisons remain required.
"""
import argparse
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import sys


def inspect_metadata():
    return dict(method='composite Gauss-Legendre on every grid interval',
                orders=[16, 32, 64], initial_z_bound=8., maximum_z_bound=96.,
                initial_grid_spacing=.25, final_grid_spacing=.125,
                maximum_evaluations=100000, logmass_tolerance=1e-9,
                gradient_tolerance=1e-8, hessian_tolerance=1e-5,
                fixed_point_only=True, certified_error_bound=False)


def load_kernel(path, expected_sha):
    """Require an explicit reviewed source pin before numerical imports."""
    path = Path(path)
    if not expected_sha or hashlib.sha256(path.read_bytes()).hexdigest() != expected_sha:
        raise ValueError('Reviewed integration-kernel source pin differs')
    spec = importlib.util.spec_from_file_location('composite_gl_reviewed_kernel', path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def integrate_gene(model, p, variance, rows, conditional, with_hessian=False, *,
                   expected_kernel_sha, kernel_path=None, settings=None):
    """Same Gaussian-prior likelihood and fixed-u derivative moments.

    The caller supplies the isolated model and its matching conditional helper.
    This does not patch model methods or change optimizer/failure rules. Numerical
    refusal is the reviewed kernel's ValueError subclass, with partial reports.
    """
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise ValueError('Composite GL requires compute')
    kernel_path = Path(__file__).with_name('repair_p3_integral.py') if kernel_path is None else kernel_path
    kernel = load_kernel(kernel_path, expected_kernel_sha)
    import numpy as np
    from scipy.special import logsumexp
    from scipy.stats import norm
    settings = kernel.Settings() if settings is None else settings
    kernel.require(os.environ.get('SLURM_JOB_ID', '').isdigit(), 'Composite GL requires compute')
    kernel.require(variance > 0 and np.isfinite(variance) and np.isfinite(p).all() and rows,
                   'Positive finite variance/parameters required')
    cp, hp = kernel.derivative_bounds(model, p, rows)
    width = model.size + 1
    sqrtv = np.sqrt(variance)
    evaluations, reports = 0, []
    rules = {}

    def refuse(message):
        failure = kernel.IntegrationFailure(message)
        failure.integration_reports = reports
        failure.likelihood_evaluations = evaluations
        raise failure

    def stats(z):
        nonlocal evaluations
        evaluations += 1
        if evaluations > settings.maximum_evaluations:
            refuse('Composite GL evaluation cap exhausted')
        gp = np.zeros(model.size)
        hh = np.zeros((width, width)) if with_hessian else None
        ll = 0.
        for row in rows:
            eta, dispersion, jac, eta_hess = model.predictors(row, p, sqrtv*z,
                                                             need_hessian=with_hessian)
            value, score, row_hess = conditional(row, eta, dispersion,
                                                  need_hessian=with_hessian)
            kernel.require(np.isfinite(value) and value <= 1e-8 and np.isfinite(score).all(),
                           'Invalid conditional probability/score')
            ll += value
            gp += jac.T @ score
            if with_hessian:
                kernel.require(np.isfinite(row_hess).all(), 'Invalid conditional Hessian')
                hh[:-1, :-1] += jac.T @ row_hess @ jac + score[0]*eta_hess
        # Prior derivatives are at fixed u, despite using z for integration.
        g = np.r_[gp, (z*z-1)/(2*variance)]
        if with_hessian:
            hh[-1, -1] = (.5-z*z)/variance**2
        vector = np.r_[1., g, (hh+np.outer(g, g)).ravel()] if with_hessian else np.r_[1., g]
        kernel.require(np.isfinite(vector).all(), 'Nonfinite derivative moments')
        return ll-z*z/2-.5*np.log(2*np.pi), vector

    def run(bound, spacing, order):
        if order not in rules:
            rules[order] = np.polynomial.legendre.leggauss(order)
        nodes, weights = rules[order]
        panels = int(np.ceil(2*bound/spacing))
        edges = np.linspace(-bound, bound, panels+1)
        size = 1+width+width*width if with_hessian else 1+width
        total, correction = np.zeros(size), np.zeros(size)
        shift = -np.inf
        for left, right in zip(edges[:-1], edges[1:]):
            half, middle = (right-left)/2, (right+left)/2
            for node, weight in zip(nodes, weights):
                z = middle+half*node
                logdensity, vector = stats(float(z))
                exponent = logdensity+np.log(half*weight)
                if exponent > shift:
                    factor = np.exp(shift-exponent)
                    total *= factor
                    correction *= factor
                    shift = exponent
                # Kahan summation of signed moments; no per-node vector cache.
                increment = np.exp(exponent-shift)*vector-correction
                updated = total+increment
                correction = (updated-total)-increment
                total = updated
        kernel.require(total[0] > 0 and np.isfinite(total).all(), 'Invalid GL integral')
        ll = shift+np.log(total[0])
        g = total[1:1+width]/total[0]
        second = total[1+width:].reshape(width, width)/total[0] if with_hessian else None
        h = second-np.outer(g, g) if with_hessian else None
        # This is a numerical mass estimate, not a certified lower mass bound.
        tail = kernel.tail_moments(bound, variance, ll, g, second, cp, hp, with_hessian)
        report = dict(z_bound=bound, grid_spacing=spacing, panels=panels, order=order,
                      node_evaluations=panels*order, cumulative_evaluations=evaluations,
                      logmass=float(ll), gradient=g, hessian=h,
                      all_grid_intervals_integrated=True, stationary_points_exhaustive=False,
                      tail_allowances_use_estimated_mass_not_certified_lower_bound=True,
                      logmass_error_allowance=float(np.log1p(tail['relative_mass_tail_allowance'])),
                      gradient_error_allowance=tail['relative_gradient_tail_allowance'],
                      hessian_error_allowance=tail['relative_hessian_tail_allowance'],
                      **tail)
        reports.append(report)
        return ll, g, h, report

    def allowances(first, second):
        mass = abs(second[0]-first[0])+np.log1p(second[3]['relative_mass_tail_allowance'])
        grad = np.abs(second[1]-first[1])+second[3]['relative_gradient_tail_allowance']
        hess = (np.abs(second[2]-first[2])+second[3]['relative_hessian_tail_allowance']
                if with_hessian else None)
        return mass, grad, hess

    def within(errors, fraction):
        return (errors[0] <= settings.logmass_tolerance*fraction and
                np.max(errors[1]) <= settings.gradient_tolerance*fraction and
                (not with_hessian or np.max(errors[2]) <= settings.hessian_tolerance*fraction))

    bound, accepted = settings.initial_z_bound, None
    for attempt in range(5):
        ladder = [run(bound, settings.grid_spacing, order) for order in (16, 32, 64)]
        errors = allowances(ladder[1], ladder[2])
        ladder[-1][3]['order_refinement_16_32_allowances'] = allowances(ladder[0], ladder[1])
        ladder[-1][3]['order_refinement_32_64_allowances'] = errors
        tail = ladder[-1][3]
        tail_ok = (tail['relative_mass_tail_allowance'] <= settings.logmass_tolerance/8 and
                   np.max(tail['relative_gradient_tail_allowance']) <= settings.gradient_tolerance/8 and
                   (not with_hessian or np.max(tail['relative_hessian_tail_allowance']) <= settings.hessian_tolerance/8))
        if tail_ok:
            if not within(errors, .25):
                refuse('GL order refinement exceeds unchanged precision allowance')
            accepted = ladder[-1]
            break
        # Expand using the same probability envelope; variance moments must also
        # pass explicit tail checks after the new domain is evaluated.
        requested = ladder[-1][0]+np.log(min(settings.logmass_tolerance,
                    settings.gradient_tolerance, settings.hessian_tolerance))
        requested -= np.log(16*(1+np.max(cp)+np.max(hp)+np.max(cp)**2))
        if np.log(2.)+norm.logsf(settings.maximum_z_bound) >= requested:
            refuse('Composite GL relative-tail domain cap exhausted')
        from scipy.optimize import brentq
        desired = brentq(lambda z: np.log(2.)+norm.logsf(z)-requested,
                         0., settings.maximum_z_bound)
        bound = max(bound+2, float(np.ceil(desired)))
        if bound+2 > settings.maximum_z_bound:
            refuse('Composite GL expanded domain exceeds cap')
    if accepted is None:
        refuse('Composite GL tail allowances remain unresolved')
    if bound+2 > settings.maximum_z_bound:
        refuse('Composite GL final domain exceeds cap')
    final_ladder = [run(bound+2, settings.grid_spacing/2, order) for order in (16, 32, 64)]
    final = final_ladder[-1]
    final_order_errors = allowances(final_ladder[1], final)
    final[3]['order_refinement_16_32_allowances'] = allowances(final_ladder[0], final_ladder[1])
    final[3]['order_refinement_32_64_allowances'] = final_order_errors
    # Combine final order disagreement and domain/grid disagreement once each;
    # tail allowance is already included in final_order_errors.
    errors = (abs(final[0]-accepted[0])+final_order_errors[0],
              np.abs(final[1]-accepted[1])+final_order_errors[1],
              np.abs(final[2]-accepted[2])+final_order_errors[2] if with_hessian else None)
    final[3]['logmass_error_allowance'] = float(final_order_errors[0])
    final[3]['gradient_error_allowance'] = final_order_errors[1]
    final[3]['hessian_error_allowance'] = final_order_errors[2]
    if not within(errors, 1.):
        refuse('Composite GL expanded domain/grid refinement disagrees')
    report = dict(method='composite Gauss-Legendre on every finite-grid interval',
                  numerical_revision_of_prescribed_20_node_approximation=True,
                  same_probability_estimand=True, fixed_point_only=True,
                  first=accepted[3], refined=final[3],
                  logmass_refinement_difference=float(abs(final[0]-accepted[0])),
                  gradient_refinement_difference=np.abs(final[1]-accepted[1]),
                  hessian_refinement_difference=(np.abs(final[2]-accepted[2]) if with_hessian else None),
                  ladder=reports, empirical_refinement_plus_tail_allowances=errors,
                  likelihood_evaluations=evaluations, repeated_nodes_recomputed=True,
                  refinement_differences_not_certified_error_bounds=True,
                  floating_arithmetic_error_not_rigorously_bounded=True)
    return final[0], final[1], final[2], report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true', required=True)
    parser.parse_args()
    print(json.dumps(inspect_metadata(), sort_keys=True))
