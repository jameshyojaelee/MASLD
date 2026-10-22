#!/usr/bin/env python3
"""Gene-local GL with prespecified16/32 agreement before requesting64.

Every finite-grid interval is integrated. Order/grid/domain agreement is an
empirical precision check, not a certified quadrature or floating-error bound.
Independent scalar and finite-difference comparisons remain required.
An initial64 result may use the full budget before the unchanged final
domain/grid/order test. An initial32 result still needs the quarter budget.
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
                orders=[16, 32, 64], stop_at32_only_with_unchanged_allowances=True,
                initial64_budget_fraction=1.,final_combined_budget_fraction=1.,
                initial_z_bound=8., maximum_z_bound=96.,
                initial_grid_spacing=.25, final_grid_spacing=.125,
                maximum_evaluations=100000, logmass_tolerance=1e-9,
                gradient_tolerance=1e-8, hessian_tolerance=1e-5,
                requested_node_batches=[32,64], GL_workspace_estimate_cap_bytes=256*1024**2,
                probability_batch_callback_required=True, same_order_Kahan=True,
                fixed_point_only=True, certified_error_bound=False)


def load_kernel(path, expected_sha):
    """Require an explicit reviewed source pin before numerical imports."""
    path = Path(path)
    if not expected_sha or hashlib.sha256(path.read_bytes()).hexdigest() != expected_sha:
        raise ValueError('Reviewed integration-kernel source pin differs')
    spec = importlib.util.spec_from_file_location('local_gl_reviewed_'+path.stem, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def GL_workspace_estimate_bytes(model, rows, batch_size):
    """Conservative declared GL workspace capacity, excluding callback/RSS."""
    width = model.Z+8
    maximum_n = max(row.n for row in rows)
    return (1024*batch_size*(maximum_n+1)+64*batch_size*width**2
            +256*sum(row.n+1 for row in rows)+64*len(rows)*(width**2+width))


def integrate_gene(model, p, variance, rows, conditional, with_hessian=False, *,
                   expected_kernel_sha, kernel_path=None, settings=None, batch_size=32):
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
    local_path = Path(__file__).with_name('repair_p3_local_moments.py')
    local_sha = '95ec9f270d18bfbbc8f31f0857308f5f2acad053661a30b208049b7a3f59f31f'
    local_module = load_kernel(local_path, local_sha)
    local = local_module.LocalGeneMoments(model,p,rows,kernel)
    cp, hp = local.bounds()
    width = len(local.active)+1
    sqrtv = np.sqrt(variance)
    evaluations, reports = 0, []
    rules = {}
    kernel.require(batch_size in (32,64), 'Fixed requested node batch must be32 or64')
    # GL integration workspace estimate only, not callback or total process RSS.
    # Callback must separately bound its support arrays before allocation, and
    # the caller must record callback+GL estimates in the execution receipt.
    estimated_bytes = GL_workspace_estimate_bytes(model,rows,batch_size)
    kernel.require(estimated_bytes <= 256*1024**2, 'Batched GL estimated256MiB integration workspace exceeded')
    predictors = []
    for row in rows:
        base = local.predictors(row,p,0.,need_hessian=with_hessian)
        for u in (-.5,.5):
            actual = local.predictors(row,p,u,need_hessian=with_hessian)
            kernel.require(actual[0] == base[0]+u*row.stage and actual[1] == base[1]
                and np.array_equal(actual[2],base[2])
                and ((actual[3] is None and base[3] is None)
                     or (actual[3] is not None and base[3] is not None and np.array_equal(actual[3],base[3]))),
                'Constant predictor/Jacobian/Hessian invariant failed')
        predictors.append((row,*base))

    def refuse(message):
        failure = kernel.IntegrationFailure(message)
        failure.integration_reports = reports
        failure.likelihood_evaluations = evaluations
        raise failure

    def stats_batch(zs):
        nonlocal evaluations
        zs = np.asarray(zs,dtype=float)
        count = len(zs)
        evaluations += count
        if evaluations > settings.maximum_evaluations:
            refuse('Composite GL evaluation cap exhausted')
        gp = np.zeros((count,len(local.active)))
        hh = np.zeros((count,width,width)) if with_hessian else None
        ll = np.zeros(count)
        us = sqrtv*zs
        for row,eta0,dispersion,jac,eta_hess in predictors:
            # Same final base+u*stage operations as the frozen scalar predictor.
            etas = eta0+us*row.stage
            for ordinal in (0,count-1):
                actual = local.predictors(row,p,float(us[ordinal]),need_hessian=with_hessian)
                kernel.require(actual[0] == etas[ordinal] and actual[1] == dispersion
                    and np.array_equal(actual[2],jac)
                    and ((actual[3] is None and eta_hess is None)
                         or (actual[3] is not None and eta_hess is not None and np.array_equal(actual[3],eta_hess))),
                    'Actual-node predictor invariance failed')
            values,scores,row_hessians = conditional(row,etas,dispersion,need_hessian=with_hessian)
            kernel.require(values.shape == (count,) and scores.shape == (count,2)
                and np.isfinite(values).all() and np.all(values <= 1e-8)
                and np.isfinite(scores).all(), 'Invalid batched conditional probability/score')
            if with_hessian:
                kernel.require(row_hessians.shape == (count,2,2) and np.isfinite(row_hessians).all(),
                               'Invalid batched conditional Hessian')
            # Keep each node's row accumulation/matmul operation order unchanged.
            for ordinal in range(count):
                ll[ordinal] += values[ordinal]
                gp[ordinal] += jac.T @ scores[ordinal]
                if with_hessian:
                    hh[ordinal,:-1,:-1] += jac.T @ row_hessians[ordinal] @ jac + scores[ordinal,0]*eta_hess
        results = []
        for ordinal,z in enumerate(zs):
            z = float(z)
            g = np.r_[gp[ordinal],(z*z-1)/(2*variance)]
            h = hh[ordinal] if with_hessian else None
            if with_hessian:
                h[-1,-1] = (.5-z*z)/variance**2
            vector = np.r_[1.,g,(h+np.outer(g,g)).ravel()] if with_hessian else np.r_[1.,g]
            kernel.require(np.isfinite(vector).all(),'Nonfinite derivative moments')
            results.append((ll[ordinal]-z*z/2-.5*np.log(2*np.pi),vector))
        return results

    def run(bound, spacing, order):
        if order not in rules:
            rules[order] = np.polynomial.legendre.leggauss(order)
        nodes, weights = rules[order]
        panels = int(np.ceil(2*bound/spacing))
        edges = np.linspace(-bound, bound, panels+1)
        size = 1+width+width*width if with_hessian else 1+width
        total, correction = np.zeros(size), np.zeros(size)
        shift = -np.inf
        pending_z,pending_weights = [],[]
        def consume(zs,logweights):
            nonlocal total,correction,shift
            allowed = min(len(zs),settings.maximum_evaluations-evaluations)
            values = stats_batch(zs[:allowed]) if allowed else []
            for (logdensity,vector),logweight in zip(values,logweights[:allowed]):
                exponent = logdensity+logweight
                if exponent > shift:
                    factor = np.exp(shift-exponent)
                    total *= factor
                    correction *= factor
                    shift = exponent
                increment = np.exp(exponent-shift)*vector-correction
                updated = total+increment
                correction = (updated-total)-increment
                total = updated
            if allowed < len(zs):
                # Same first over-cap logical node and refusal counter as scalar.
                stats_batch(zs[allowed:allowed+1])
        for left,right in zip(edges[:-1],edges[1:]):
            half,middle = (right-left)/2,(right+left)/2
            for node,weight in zip(nodes,weights):
                pending_z.append(float(middle+half*node))
                pending_weights.append(np.log(half*weight))
                if len(pending_z) == batch_size:
                    consume(pending_z,pending_weights)
                    pending_z,pending_weights = [],[]
        if pending_z:
            consume(pending_z,pending_weights)
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

    def tail_passes(tail):
        return (tail['relative_mass_tail_allowance'] <= settings.logmass_tolerance/8 and
                np.max(tail['relative_gradient_tail_allowance']) <= settings.gradient_tolerance/8 and
                (not with_hessian or np.max(tail['relative_hessian_tail_allowance']) <= settings.hessian_tolerance/8))

    bound, accepted = settings.initial_z_bound, None
    for attempt in range(5):
        ladder = [run(bound, settings.grid_spacing, order) for order in (16, 32)]
        errors = allowances(ladder[0], ladder[1])
        ladder[-1][3]['order_refinement_16_32_allowances'] = errors
        ladder[-1][3]['order_refinement_32_64_allowances'] = None
        if tail_passes(ladder[-1][3]) and not within(errors, .25):
            ladder.append(run(bound, settings.grid_spacing, 64))
            errors = allowances(ladder[1], ladder[2])
            ladder[-1][3]['order_refinement_16_32_allowances'] = ladder[1][3]['order_refinement_16_32_allowances']
            ladder[-1][3]['order_refinement_32_64_allowances'] = errors
        if tail_passes(ladder[-1][3]):
            # Order64 may proceed to the unchanged final domain/grid check.
            if not within(errors, 1. if ladder[-1][3]['order']==64 else .25):
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
    final_ladder = [run(bound+2, settings.grid_spacing/2, order) for order in (16, 32)]
    final = final_ladder[-1]
    final_order_errors = allowances(final_ladder[0], final)
    final[3]['order_refinement_16_32_allowances'] = final_order_errors
    final[3]['order_refinement_32_64_allowances'] = None
    # Combine final order disagreement and domain/grid disagreement once each;
    # tail allowance is already included in final_order_errors.
    errors = (abs(final[0]-accepted[0])+final_order_errors[0],
              np.abs(final[1]-accepted[1])+final_order_errors[1],
              np.abs(final[2]-accepted[2])+final_order_errors[2] if with_hessian else None)
    if not within(final_order_errors, .25) or not within(errors, 1.):
        final_ladder.append(run(bound+2, settings.grid_spacing/2, 64))
        final = final_ladder[-1]
        final_order_errors = allowances(final_ladder[1], final)
        final[3]['order_refinement_16_32_allowances'] = final_ladder[1][3]['order_refinement_16_32_allowances']
        final[3]['order_refinement_32_64_allowances'] = final_order_errors
        errors = (abs(final[0]-accepted[0])+final_order_errors[0],
                  np.abs(final[1]-accepted[1])+final_order_errors[1],
                  np.abs(final[2]-accepted[2])+final_order_errors[2] if with_hessian else None)
    final[3]['logmass_error_allowance'] = float(final_order_errors[0])
    final[3]['gradient_error_allowance'] = final_order_errors[1]
    final[3]['hessian_error_allowance'] = final_order_errors[2]
    if not within(errors, 1.):
        refuse('Composite GL expanded domain/grid refinement disagrees')
    report = dict(method='adaptive16/32/64; initial64 full budget; unchanged final domain/grid precision',
                  selected_first_order=accepted[3]['order'],selected_refined_order=final[3]['order'],
                  accepted_order32_requires_16_32_and_total_allowances=True,
                  active_parameter_indices=local.active.tolist(),
                  augmented_global_indices=local.augmented.tolist(),
                  full_parameter_dimension=model.size, local_augmented_dimension=width,
                  ladder_arrays_are_local_coordinates=True, interface_arrays_embedded_only_at_return=True,
                  requested_node_batch_size=batch_size, estimated_GL_workspace_bytes=estimated_bytes,
                  GL_workspace_excludes_callback_and_process_RSS=True,
                  workspace_estimate_not_measured_peak=True, predictor_invariants_checked=True,
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
    # All integrations, covariance corrections, refinement and tail checks above
    # use active coordinates. Inactive moments are identically zero; max-budget
    # comparisons are therefore unchanged. Allocate the full interface only now.
    return final[0], local.embed_vector(final[1]), (local.embed_matrix(final[2])
        if with_hessian else None), report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect', action='store_true', required=True)
    parser.parse_args()
    print(json.dumps(inspect_metadata(), sort_keys=True))
