"""Bounded synthetic native selected-count probabilities in long double.

Same oriented three-genotype mixture and integer inclusion support. Rising
products use their integer factors, with analytic eta/r first and second
derivatives. The inclusion probability is summed directly on its full support.
No fitting or production model replacement is provided.
"""
import os


class LongDoubleConditional:
    def __init__(self, model, p, rows):
        if not os.environ.get('SLURM_JOB_ID', '').isdigit():
            raise RuntimeError('Compute allocation required')
        import numpy as np
        self.np = np
        self.dtype = np.longdouble
        if np.finfo(np.longdouble).nmant < 63:
            raise RuntimeError('At least 64-bit significand required')
        self.model = model
        self.p = np.asarray(p).copy()
        self.p.setflags(write=False)
        self.parameter_bytes = self.p.tobytes()
        self.rows = tuple(rows)
        if (model.G, model.Z, model.size) != (8, 7, 28) or not 0 < len(self.rows) <= 32:
            raise ValueError('Saved synthetic G8/Z7/P28 scope required')
        if len({r.gene for r in self.rows}) != 1 or sum(r.n+1 for r in self.rows) > 100000:
            raise ValueError('One gene and bounded count workspace required')
        self.signatures = tuple(tuple(vars(r).items()) for r in self.rows)
        self.prepared = {}
        ld = self.dtype
        for row in self.rows:
            n = row.n
            lo = max(row.minor, (row.percent*n+99)//100)
            if int(n) != n or not 0 < n <= 5000 or n < row.depth or not lo <= row.a <= n-lo:
                raise ValueError('Original integer inclusion required')
            support = np.arange(lo, n-lo+1)
            k = np.arange(1, n+1, dtype=ld)
            choose = np.r_[ld(0), np.cumsum(np.log((ld(n)+1-k)/k), dtype=ld)]
            r = model.predictors(row, p, 0., need_hessian=False)[1]
            s, q = self.sigmoid(ld(r)), self.sigmoid(-ld(r))
            rho = ld(1e-6)+ld(1-2e-6)*s
            rp = ld(1-2e-6)*s*q
            rpp = rp*(q-s)
            t = (1-rho)/rho
            tp = -rp/rho**2
            tpp = -rpp/rho**2+2*rp**2/rho**3
            rt, dt, qt = self.rising(n, np.log(t))
            hom = []
            for mean in (row.e, 1-row.e):
                m = ld(mean)
                if row.phi == 0:
                    logs = choose[support]+support*np.log(m)+(n-support)*np.log1p(-m)
                else:
                    ht = (1-ld(row.phi))/ld(row.phi)
                    ha = self.rising(n, np.log(ht*m))[0]
                    hb = self.rising(n, np.log(ht*(1-m)))[0]
                    hd = self.rising(n, np.log(ht))[0][-1]
                    logs = choose[support]+ha[support]+hb[n-support]-hd
                hom.append(logs)
            with np.errstate(divide='ignore'):
                prior = np.log(np.asarray(row.prior, dtype=ld))
            self.prepared[id(row)] = dict(row=row, r=r, n=n, support=support,
                observed=row.a-lo, choose=choose[support], logt=np.log(t),
                rt=rt[-1], dt=dt[-1], qt=qt[-1], u=tp/t, w=tpp/t-(tp/t)**2,
                logprior=prior, hom=hom)

    def sigmoid(self, x):
        np = self.np
        e = np.exp(-np.abs(x))
        return np.where(x >= 0, 1/(1+e), e/(1+e))

    def rising(self, n, logshape):
        np, ld = self.np, self.dtype
        j = np.arange(1, n, dtype=ld)
        logj = np.log(j)
        values = np.r_[ld(logshape), np.logaddexp(logshape, logj)]
        d = self.sigmoid(logshape-logj)
        scores = np.r_[ld(1), d]
        curvatures = np.r_[ld(0), d*self.sigmoid(logj-logshape)]
        return tuple(np.r_[ld(0), np.cumsum(a, dtype=ld)] for a in (values, scores, curvatures))

    def terms(self, row, eta, r, need_hessian=True):
        np, ld = self.np, self.dtype
        if (self.p.tobytes() != self.parameter_bytes
                or tuple(tuple(vars(x).items()) for x in self.rows) != self.signatures):
            raise ValueError('Fixed context changed')
        d = self.prepared.get(id(row))
        if d is None or d['row'] is not row or r != d['r'] or not np.isfinite(eta):
            raise ValueError('Row/dispersion context differs')
        eta = ld(eta)-ld(row.orientation)*ld(row.omega)
        m, q = self.sigmoid(eta), self.sigmoid(-eta)
        ra, da, qa = self.rising(d['n'], d['logt']-np.logaddexp(ld(0), -eta))
        rb, db, qb = self.rising(d['n'], d['logt']-np.logaddexp(ld(0), eta))
        k, nk = d['support'], d['n']-d['support']
        da, db, qa, qb = da[k], db[nk], qa[k], qb[nk]
        lp = d['choose']+ra[k]+rb[nk]-d['rt']
        score = np.column_stack((da*q-db*m, (da+db-d['dt'])*d['u']))
        h = None
        if need_hessian:
            h = np.empty((len(k), 2, 2), dtype=ld)
            h[:, 0, 0] = qa*q*q+qb*m*m-(da+db)*m*q
            h[:, 0, 1] = h[:, 1, 0] = d['u']*(qa*q-qb*m)
            h[:, 1, 1] = d['u']**2*(qa+qb-d['qt'])+d['w']*(da+db-d['dt'])
        parts = np.column_stack((lp+d['logprior'][0], d['hom'][0]+d['logprior'][1],
                                 d['hom'][1]+d['logprior'][2]))
        pivot = np.max(parts, axis=1)
        weights = np.exp(parts-pivot[:, None])
        total = weights.sum(axis=1, dtype=ld)
        logmass = pivot+np.log(total)
        wh = weights[:, 0]/total
        gs = wh[:, None]*score
        hs = (wh[:, None, None]*h+(wh*(1-wh))[:, None, None]*np.einsum('ni,nj->nij', score, score)
              if need_hessian else None)
        pivot = np.max(logmass)
        den = np.exp(logmass-pivot)
        den /= den.sum(dtype=ld)
        dg = np.einsum('n,ni->i', den, gs)
        dh = (np.einsum('n,nij->ij', den, hs)
              +np.einsum('n,ni,nj->ij', den, gs, gs)-np.outer(dg, dg)
              if need_hessian else None)
        obs = d['observed']
        value = logmass[obs]-pivot-np.log(np.exp(logmass-pivot).sum(dtype=ld))
        answer = value, gs[obs]-dg, hs[obs]-dh if need_hessian else None
        if not all(a is None or np.isfinite(a).all() for a in answer):
            raise ValueError('Nonfinite selected-count result')
        # Conversion occurs only at the scalar diagnostic interface.
        return float(answer[0]), np.asarray(answer[1], dtype=float), np.asarray(answer[2], dtype=float) if need_hessian else None
