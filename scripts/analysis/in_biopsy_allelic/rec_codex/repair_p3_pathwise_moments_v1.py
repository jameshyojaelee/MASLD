"""Fixed synthetic derivative diagnostic under u=sqrt(v)*z; no fitter.

The standard normal measure is fixed. Derivatives are with respect to v,
not tau. Conditional eta is affine in u and its nuisance Jacobian does not
depend on u. Tail envelopes apply to the selected RNA-count likelihood.
"""
import math
import os


class PathwiseGeneMoments:
    def __init__(self, local):
        if not os.environ.get('SLURM_JOB_ID', '').isdigit():
            raise RuntimeError('Compute allocation required')
        self.local = local
        self.np = local.np

    def point(self, p, z, variance, conditional, with_hessian=True):
        local, np = self.local, self.np
        local.validate(p)
        if not math.isfinite(z) or not math.isfinite(variance) or variance <= 0:
            raise ValueError('Finite z and positive variance required')
        root = math.sqrt(variance)
        width = len(local.active)
        ll, gu = 0., 0.
        gp = np.zeros(width)
        h = np.zeros((width+1, width+1)) if with_hessian else None
        for row in local.rows:
            eta, r, jac, eh = local.predictors(row, p, root*z, with_hessian)
            value, score, rh = conditional(row, eta, r, need_hessian=with_hessian)
            ll += value
            gp += jac.T @ score
            gu += row.stage*score[0]
            if with_hessian:
                h[:-1, :-1] += jac.T @ rh @ jac + score[0]*eh
                h[:-1, -1] += z/(2*root)*row.stage*(rh[0] @ jac)
                h[-1, -1] += z*z/(4*variance)*row.stage**2*rh[0, 0]
        gradient = np.r_[gp, z/(2*root)*gu]
        if with_hessian:
            h[-1, :-1] = h[:-1, -1]
            h[-1, -1] -= z/(4*variance*root)*gu
        second = h + np.outer(gradient, gradient) if with_hessian else None
        return ll-z*z/2-.5*math.log(2*math.pi), gradient, h, second

    def envelopes(self):
        """Absolute conditional nuisance/u score and curvature bounds."""
        local, np = self.local, self.np
        cp, hp = local.bounds()
        cu, huu = 0., 0.
        hpu = np.zeros(len(cp))
        for row in local.rows:
            _, _, jac, _ = local.predictors(row, local.p, 0., False)
            b = 16*row.n + 8*row.n**2
            cu += 2*row.n*abs(row.stage)
            huu += b*row.stage**2
            hpu += b*abs(row.stage)*np.abs(jac).sum(axis=0)
        return cp, hp, cu, huu, hpu

    def shifted_tail_bounds(self, bound, variance, log_shift, with_hessian=True):
        """Absolute numerator tails after division by exp(log_shift).

        L<=1. M0, M1, M2 integrate 1, |z|, z^2 against the two normal tails.
        These bounds do not cover floating point or quadrature errors.
        """
        np = self.np
        if bound <= 0 or variance <= 0:
            raise ValueError('Positive bound/variance required')
        cp, hp, cu, huu, hpu = self.envelopes()
        phi = math.exp(-.5*bound*bound-.5*math.log(2*math.pi)-log_shift)
        m0 = math.exp(math.log(math.erfc(bound/math.sqrt(2)))-log_shift)
        m1 = 2*phi
        m2 = 2*bound*phi + m0
        root = math.sqrt(variance)
        gt = np.r_[cp*m0, cu/(2*root)*m1]
        st = None
        if with_hessian:
            st = np.zeros((len(cp)+1, len(cp)+1))
            st[:-1, :-1] = (hp+np.outer(cp, cp))*m0
            st[:-1, -1] = (hpu+cp*cu)/(2*root)*m1
            st[-1, :-1] = st[:-1, -1]
            st[-1, -1] = (huu+cu*cu)/(4*variance)*m2 + cu/(4*variance*root)*m1
        return m0, gt, st
