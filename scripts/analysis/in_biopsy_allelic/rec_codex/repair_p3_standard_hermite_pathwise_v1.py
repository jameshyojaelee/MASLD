"""Fixed-measure Hermite diagnostic for saved synthetic G8 gene integrals.

This differs from the original adaptive20-node approximation. Standard normal
nodes stay fixed while differentiating u=sqrt(v)z. Order agreement is empirical,
not a quadrature certificate. No optimizer or fallback to another backend.
"""
import math
import os


def integrate(local, conditional, moments, variance):
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('Compute required')
    np = local.np
    ld = np.longdouble
    records = []
    results = []
    for order in (32, 64, 128, 256):
        zs, weights = np.polynomial.hermite_e.hermegauss(order)
        weights = weights/math.sqrt(2*math.pi)
        if (not np.isfinite(zs).all() or not np.isfinite(weights).all()
                or not np.all(weights > 0)
                or abs(float(weights.sum())-1) > 1e-12
                or abs(float(weights @ (zs*zs))-1) > 1e-10
                or abs(float(weights @ (zs**4))-3) > 1e-10):
            raise ValueError('Standard normal Hermite rule failed')
        logmasses, gradients, seconds = [], [], []
        for z in zs:
            row_values = []
            def callback(row, eta, r, need_hessian=True):
                result = conditional.terms(row, eta, r, need_hessian=need_hessian)
                row_values.append(result[0])
                return result
            _, g, _, second = moments.point(local.p, float(z), variance, callback)
            # Avoid subtracting and then adding the Gaussian log density.
            logmasses.append(math.fsum(row_values))
            gradients.append(g)
            seconds.append(second)
        logs = np.asarray(logmasses, dtype=ld)
        shift = np.max(logs)
        weighted = np.asarray(weights, dtype=ld)*np.exp(logs-shift)
        mass = weighted.sum(dtype=ld)
        if not np.isfinite(mass) or mass <= 0:
            raise ValueError('Invalid Hermite mass')
        probabilities = weighted/mass
        gradient = np.sum(probabilities[:, None]*np.asarray(gradients, dtype=ld), axis=0, dtype=ld)
        second = np.sum(probabilities[:, None, None]*np.asarray(seconds, dtype=ld), axis=0, dtype=ld)
        hessian = second-np.outer(gradient, gradient)
        result = (float(shift+np.log(mass)), np.asarray(gradient, dtype=float), np.asarray(hessian, dtype=float))
        if not all(np.isfinite(x).all() for x in result):
            raise ValueError('Nonfinite Hermite moments')
        results.append(result)
        record = dict(order=order, logmass=result[0], gradient=local.embed_vector(result[1]),
                      hessian=local.embed_matrix(result[2]), nodes=order)
        if len(results) > 1:
            old = results[-2]
            differences = [abs(result[0]-old[0]), float(np.max(np.abs(result[1]-old[1]))),
                           float(np.max(np.abs(result[2]-old[2])))]
            record['previous_order_differences'] = differences
            record['nominal_order_check_passed'] = all(a <= b for a, b in zip(differences, (1e-9, 1e-8, 1e-5)))
        records.append(record)
        if order >= 128 and record['nominal_order_check_passed']:
            return dict(returned=True, gene=local.g, logmass=result[0],
                        gradient=local.embed_vector(result[1]), hessian=local.embed_matrix(result[2]),
                        orders=records, node_evaluations=sum(x['nodes'] for x in records),
                        standard_normal_fixed_measure=True, empirical_not_certified=True)
    return dict(returned=False, gene=local.g, orders=records,
                node_evaluations=sum(x['nodes'] for x in records),
                refusal='Standard Hermite order refinement exceeds unchanged limits',
                empirical_not_certified=True)
