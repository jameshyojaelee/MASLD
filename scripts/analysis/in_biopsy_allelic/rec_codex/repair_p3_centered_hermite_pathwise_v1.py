"""Synthetic fixed-point importance Hermite quadrature; no fitting.

Integrate against a defensive mixture of the original normal and normal
proposals centered at sampled posterior peaks. Proposal parameters are held
fixed during moment evaluation. Pathwise derivatives concern the original
standard-normal variable, not the proposal. Empirical order agreement alone
does not certify missing modes, tails or floating point accuracy.
"""
import math
import os


def integrate(local, conditional, moments, variance, with_hessian=True):
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('Compute required')
    from scipy.optimize import minimize_scalar
    from scipy.special import logsumexp, roots_hermitenorm
    np = local.np
    ld = np.longdouble
    if not math.isfinite(variance) or variance <= 0:
        raise ValueError('Positive finite variance required')
    root = math.sqrt(variance)
    before = local.p.tobytes()
    probes = 0

    def posterior(z, curvature=False):
        nonlocal probes
        probes += 1
        values, second = [], -1.
        for row in local.rows:
            eta, r = local.predictors(row, local.p, root*z, curvature)[:2]
            value, _, hessian = conditional.terms(row, eta, r, need_hessian=curvature)
            values.append(value)
            if curvature:
                second += variance*row.stage**2*hessian[0, 0]
        return math.fsum(values)-z*z/2, float(second)

    grid = np.linspace(-16., 16., 65)
    values = [posterior(float(z))[0] for z in grid]
    candidates = [i for i in range(1, len(grid)-1)
                  if values[i] >= values[i-1] and values[i] >= values[i+1]]
    if values[0] >= values[1] or values[-1] >= values[-2] or not candidates or len(candidates) > 4:
        raise ValueError('Posterior peak scan unresolved')
    peaks = []
    for i in candidates:
        fit = minimize_scalar(lambda z: -posterior(z)[0], bounds=(grid[i-1], grid[i+1]),
                              method='bounded', options={'xatol': 1e-10, 'maxiter': 100})
        if not fit.success or not math.isfinite(fit.x):
            raise ValueError('Posterior peak refinement failed')
        height, curvature = posterior(float(fit.x), True)
        if not math.isfinite(curvature) or curvature >= 0:
            raise ValueError('Posterior peak has nonnegative curvature')
        sigma = min(1., max(1e-4, math.sqrt(-1/curvature)))
        peaks.append(dict(center=float(fit.x), scale=sigma, logheight=height, curvature=curvature))
    # The original prior has positive weight everywhere; each detected peak
    # also gets a wider proposal. None of these components is a new model.
    proposals = [(0., 1., .1)]
    for peak in peaks:
        proposals.extend([(peak['center'], peak['scale'], .45/len(peaks)),
                          (peak['center'], 2*peak['scale'], .45/len(peaks))])
    records, results = [], []
    for order in (32, 64, 128, 256):
        nodes, weights = roots_hermitenorm(order)
        weights = weights/math.sqrt(2*math.pi)
        if (not np.isfinite(nodes).all() or not np.all(weights > 0)
                or abs(float(weights.sum())-1) > 1e-12
                or abs(float(weights @ (nodes*nodes))-1) > 1e-10
                or abs(float(weights @ (nodes**4))-3) > 1e-10):
            raise ValueError('Normal rule failed')
        logs, gradients, seconds = [], [], []
        for center, scale, alpha in proposals:
            for node, weight in zip(nodes, weights):
                z = center+scale*float(node)
                row_values = []
                def callback(row, eta, r, need_hessian=True):
                    result = conditional.terms(row, eta, r, need_hessian=need_hessian)
                    row_values.append(result[0])
                    return result
                _, g, _, second = moments.point(local.p, z, variance, callback, with_hessian=with_hessian)
                # f(z)*phi(z)/q_mix(z), summed over alpha_j E_qj.
                logq = float(logsumexp([math.log(a)-math.log(s)-.5*math.log(2*math.pi)
                    -.5*((z-m)/s)**2 for m, s, a in proposals]))
                logs.append(math.fsum(row_values)-z*z/2-.5*math.log(2*math.pi)-logq
                            +math.log(alpha)+math.log(float(weight)))
                gradients.append(g)
                seconds.append(second)
        logs = np.asarray(logs, dtype=ld)
        shift = np.max(logs)
        weighted = np.exp(logs-shift)
        mass = weighted.sum(dtype=ld)
        probabilities = weighted/mass
        gradient = np.sum(probabilities[:, None]*np.asarray(gradients, dtype=ld), axis=0, dtype=ld)
        second = (np.sum(probabilities[:, None, None]*np.asarray(seconds, dtype=ld), axis=0, dtype=ld)
                  if with_hessian else None)
        hessian = second-np.outer(gradient, gradient) if with_hessian else None
        result = (float(shift+np.log(mass)), np.asarray(gradient, dtype=float),
                  np.asarray(hessian, dtype=float) if with_hessian else None)
        if not all(x is None or np.isfinite(x).all() for x in result):
            raise ValueError('Nonfinite importance moments')
        results.append(result)
        record = dict(order=order, logmass=result[0], gradient=local.embed_vector(result[1]),
                      hessian=local.embed_matrix(result[2]) if with_hessian else None,
                      nodes=order*len(proposals))
        if len(results) > 1:
            old = results[-2]
            differences = [abs(result[0]-old[0]), float(np.max(np.abs(result[1]-old[1]))),
                           float(np.max(np.abs(result[2]-old[2]))) if with_hessian else 0.]
            record['previous_order_differences'] = differences
            record['nominal_order_check_passed'] = all(a <= b for a, b in zip(differences, (1e-9, 1e-8, 1e-5)))
        records.append(record)
        if order >= 128 and record['nominal_order_check_passed']:
            break
    if local.p.tobytes() != before:
        raise ValueError('Fixed point changed')
    returned = records[-1]['nominal_order_check_passed']
    return dict(returned=returned, gene=local.g, logmass=result[0],
                gradient=local.embed_vector(result[1]), hessian=local.embed_matrix(result[2]) if with_hessian else None,
                orders=records, peaks=peaks, proposals=proposals, peak_probe_evaluations=probes,
                node_evaluations=sum(x['nodes'] for x in records),
                original_standard_normal_measure=True, proposal_derivatives_excluded=True,
                empirical_not_certified=True)
