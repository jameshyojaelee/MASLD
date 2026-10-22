"""Original sampled/refined peak scan only; no importance quadrature.

Peaks supply adaptive interval breakpoints and a numerical log shift. They
do not restrict the likelihood integral to the discovered modes.
"""
import math
import os


class PeakRefusal(ValueError):
    pass


def scan(local, conditional, variance):
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('Compute required')
    from scipy.optimize import minimize_scalar
    np = local.np
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
        raise PeakRefusal('Posterior peak scan unresolved')
    peaks = []
    for i in candidates:
        fit = minimize_scalar(lambda z: -posterior(z)[0], bounds=(grid[i-1], grid[i+1]),
                              method='bounded', options={'xatol': 1e-10, 'maxiter': 100})
        if not fit.success or not math.isfinite(fit.x):
            raise PeakRefusal('Posterior peak refinement failed')
        height, curvature = posterior(float(fit.x), True)
        if not math.isfinite(curvature) or curvature >= 0:
            raise PeakRefusal('Posterior peak has nonnegative curvature')
        sigma = min(1., max(1e-4, math.sqrt(-1/curvature)))
        peaks.append(dict(center=float(fit.x), scale=sigma, logheight=height, curvature=curvature))
    if local.p.tobytes() != before:
        raise RuntimeError('Fixed peak scan point changed')
    return dict(peaks=peaks, probe_evaluations=probes)
