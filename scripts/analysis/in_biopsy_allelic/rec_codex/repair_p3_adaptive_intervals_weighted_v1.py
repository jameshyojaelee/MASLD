"""Synthetic adaptive moments with probability-weighted conditioning units.

Same integral, derivative-tail bounds and final acceptance as the original
worker; only the positive vector conditioning scales change. Their values
are sampled magnitudes, not certified bounds. No fitting or fallback.
"""
import math

WORK = None


def worker(gene):
    import numpy as np
    from scipy.integrate import quad_vec
    model, p, variance, peaks, modules = WORK
    local = modules['repair_p3_local_moments'].LocalGeneMoments(
        model, p, tuple(r for r in model.rows if r.gene == gene), None)
    conditional = modules['repair_p3_longdouble_conditional_v1'].LongDoubleConditional(model, p, local.rows)
    moments = modules['repair_p3_pathwise_moments_v1'].PathwiseGeneMoments(local)
    shift = max(x['logheight'] for x in peaks[gene])-.5*math.log(2*math.pi)
    # The failed importance rule contributes only locations, not any integral.
    points = sorted({x['center']+multiple*x['scale'] for x in peaks[gene]
                     for multiple in (-8., -4., -2., -1., 0., 1., 2., 4., 8.)})
    def raw(z):
        row_values = []
        def callback(row, eta, r, need_hessian=True):
            value = conditional.terms(row, eta, r, need_hessian=need_hessian)
            row_values.append(value[0])
            return value
        _, g, _, second = moments.point(p, float(z), variance, callback)
        logmass = math.fsum(row_values)-z*z/2-.5*math.log(2*math.pi)
        return logmass, np.r_[1., g, second.ravel()]
    width = len(local.active)+1
    scales = np.ones(1+width+width*width)
    for z in (0., -1., 1., *points):
        logmass, vector = raw(z)
        # Conditioning scales are arbitrary positive units, not tail bounds.
        scales = np.maximum(scales, math.exp(logmass-shift)*np.abs(vector))
    def integrand(z):
        logmass, vector = raw(z)
        return math.exp(logmass-shift)*vector/scales
    records = []
    for bound in (24., 26.):
        values, error, info = quad_vec(integrand, -bound, bound,
            points=[z for z in points if -bound < z < bound], epsabs=1e-12,
            epsrel=1e-12, norm='max', cache_size=4*1024**2,
            limit=2048, workers=1, full_output=True)
        record = dict(bound=bound, success=bool(info.success), status=int(info.status),
            message=info.message, evaluations=info.neval, error_estimate=float(error),
            scaled_values=values, scales=scales, log_shift=shift)
        records.append(record)
        if not info.success or not np.isfinite(values).all() or not math.isfinite(error) or error < 0:
            return dict(gene=gene, returned=False, records=records, refusal='Adaptive quadrature refused')
        values = values*scales
        mass = float(values[0])
        mt, gt, st = moments.shifted_tail_bounds(bound, variance, shift)
        lower = mass-scales[0]*error-mt
        if lower <= 0:
            return dict(gene=gene, returned=False, records=records, refusal='Mass allowance unresolved')
        gradient = values[1:1+width]/mass
        second = values[1+width:].reshape(width, width)/mass
        hessian = second-np.outer(gradient, gradient)
        mass_error = scales[0]*error+mt
        ga = (scales[1:1+width]*error+gt+np.abs(gradient)*mass_error)/lower
        sa = (scales[1+width:].reshape(width, width)*error+st+np.abs(second)*mass_error)/lower
        ha = sa+np.outer(np.abs(gradient), ga)+np.outer(ga, np.abs(gradient))+np.outer(ga, ga)
        record.update(logmass=shift+math.log(mass), gradient=local.embed_vector(gradient),
            hessian=local.embed_matrix(hessian), shifted_moments=values,
            estimated_logmass_allowance=mass_error/lower,
            estimated_gradient_allowance=local.embed_vector(ga),
            estimated_hessian_allowance=local.embed_matrix(ha),
            shifted_tail_bounds=[mt, gt, st])
    a, b = records
    deltas = [abs(a['logmass']-b['logmass']), np.abs(a['gradient']-b['gradient']),
              np.abs(a['hessian']-b['hessian'])]
    allowed = [deltas[0]+a['estimated_logmass_allowance']+b['estimated_logmass_allowance'],
        float(np.max(deltas[1]+a['estimated_gradient_allowance']+b['estimated_gradient_allowance'])),
        float(np.max(deltas[2]+a['estimated_hessian_allowance']+b['estimated_hessian_allowance']))]
    return dict(gene=gene, returned=all(x <= y for x, y in zip(allowed, (1e-9, 1e-8, 1e-5))),
        records=records, differences_plus_allowances=allowed, logmass=b['logmass'],
        gradient=b['gradient'], hessian=b['hessian'], empirical_not_certified=True)
