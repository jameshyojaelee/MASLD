"""v2: per-call eta prefixes shared across supports; no persistent u cache."""
import os
import sys
import time


class PreparedConditional:
    def __init__(self, model, p, rows):
        if not os.environ.get('SLURM_JOB_ID', '').isdigit():
            raise RuntimeError('Compute allocation required')
        import numpy as np
        from scipy.special import expit, gammaln
        self.np, self.expit = np, expit
        self.original = sys.modules[model.evaluate.__func__.__module__]
        self.model, self.p = model, np.asarray(p).copy()
        if self.p.ndim != 1 or not np.isfinite(self.p).all():
            raise ValueError('Finite fixed parameter vector required')
        self.rows = {}
        for row in rows:
            _, r, _, _ = model.predictors(row, self.p, 0., need_hessian=False)
            n = row.n
            lo = max(row.minor, (row.percent*n+99)//100)
            if n < row.depth or not lo <= row.a <= n-lo:
                raise ValueError('Row fails original integer inclusion rule')
            vr, qr = expit(r), expit(-r)
            rho = 1e-6 + (1-2e-6)*vr
            rp = (1-2e-6)*vr*qr
            rpp = rp*(qr-vr)
            t = (1-rho)/rho
            tp = -rp/rho**2
            tpp = -rpp/rho**2+2*rp**2/rho**3
            logt = np.log(t)
            logj = np.log(np.arange(1,n,dtype=float))
            rt, dt, qt = self.rising(n, logt, logj)
            tails = np.r_[np.arange(lo),np.arange(n-lo+1,n+1)]
            supports = {}
            for k in (np.array([row.a]), tails, np.arange(lo,n-lo+1)):
                supports[tuple(k)] = (k, gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1))
            with np.errstate(divide='ignore'):
                logprior = np.log(row.prior)
            hom1 = self.original.fixed_bb(np.array([row.a]),n,row.e,row.phi)
            hom2 = self.original.fixed_bb(np.array([row.a]),n,1-row.e,row.phi)
            self.rows[id(row)] = dict(row=row, signature=tuple(vars(row).items()), r=r,
                n=n, lo=lo, logt=logt, logj=logj, rt=rt[n], dt=dt[n], qt=qt[n],
                u=tp/t, w=tpp/t-(tp/t)**2, tails=tails, supports=supports,
                logprior=logprior, hom1=hom1, hom2=hom2,
                hom=self.original.hom_log_inclusion(n,lo,row.e,row.phi))

    def rising(self, n, logshape, logj):
        np, expit = self.np, self.expit
        if n == 0:
            return np.zeros(1),np.zeros(1),np.zeros(1)
        value = np.r_[logshape,np.logaddexp(logshape,logj)]
        score = np.r_[1.,expit(logshape-logj)]
        curvature = np.r_[0.,expit(logshape-logj)*expit(logj-logshape)]
        return tuple(np.r_[0.,np.cumsum(a)] for a in (value,score,curvature))

    def validate(self, p):
        np = self.np
        if np.asarray(p).shape != self.p.shape or not np.array_equal(p,self.p):
            raise ValueError('Prepared context cannot span changed parameters')

    def eta_prefixes(self, eta, d):
        np, expit = self.np, self.expit
        m,q = expit(eta),expit(-eta)
        ra,da,qa = self.rising(d['n'],d['logt']-np.logaddexp(0.,-eta),d['logj'])
        rb,db,qb = self.rising(d['n'],d['logt']-np.logaddexp(0.,eta),d['logj'])
        return m,q,ra,da,qa,rb,db,qb

    def bb(self, k, d, prefixes, need_hessian):
        np = self.np
        m,q,ra,da,qa,rb,db,qb = prefixes
        k,choose = d['supports'][tuple(k)]
        nk = d['n']-k
        da,db,qa,qb = da[k],db[nk],qa[k],qb[nk]
        mass = choose+ra[k]+rb[nk]-d['rt']
        u,w = d['u'],d['w']
        score = np.column_stack((da*q-db*m,(da+db-d['dt'])*u))
        if not need_hessian:
            return mass,score,None
        hess = np.empty((len(k),2,2))
        hess[:,0,0] = qa*q*q+qb*m*m-(da+db)*m*q
        hess[:,0,1] = hess[:,1,0] = u*(qa*q-qb*m)
        hess[:,1,1] = u*u*(qa+qb-d['qt'])+w*(da+db-d['dt'])
        return mass,score,hess

    def terms(self, row, eta, linear_r, p, need_hessian=True):
        self.validate(p)
        np, orig = self.np,self.original
        d = self.rows.get(id(row))
        if d is None or d['row'] is not row or d['signature'] != tuple(vars(row).items()):
            raise ValueError('Prepared row changed')
        if linear_r != d['r'] or not np.isfinite(eta):
            raise ValueError('Prepared dispersion changed or eta nonfinite')
        eta = eta-row.orientation*row.omega
        prefixes = self.eta_prefixes(eta,d)
        if not len(d['tails']):
            logc,score,hess = 0.,np.zeros(2),np.zeros((2,2)) if need_hessian else None
        else:
            logtail,score,hess = orig.sum_terms(*self.bb(d['tails'],d,prefixes,need_hessian))
            if logtail < np.log(.5):
                logc = np.log1p(-np.exp(logtail))
                ratio = np.exp(logtail-logc)
                c_score = -ratio*score
                hess = (-ratio*(hess+np.outer(score,score))-np.outer(c_score,c_score)
                        if need_hessian else None)
                score = c_score
            else:
                logc,score,hess = orig.sum_terms(*self.bb(np.arange(d['lo'],d['n']-d['lo']+1),d,prefixes,need_hessian))
        logs = d['logprior']+np.array([logc,d['hom'],d['hom']])
        logden,ds,dh = orig.sum_terms(logs,np.vstack((score,np.zeros((2,2)))),
            np.concatenate((hess[None,:,:],np.zeros((2,2,2)))) if need_hessian else None)
        lp,grad,hess = self.bb(np.array([row.a]),d,prefixes,need_hessian)
        parts = np.stack((lp+d['logprior'][0],d['hom1']+d['logprior'][1],d['hom2']+d['logprior'][2]))
        from scipy.special import logsumexp
        mass = logsumexp(parts,axis=0)
        wh = np.exp(parts[0]-mass)
        ns = wh[:,None]*grad
        nh = (wh[:,None,None]*hess+(wh*(1-wh))[:,None,None]*np.einsum('ni,nj->nij',grad,grad)
              if need_hessian else None)
        return float(mass[0]-logden),ns[0]-ds,nh[0]-dh if need_hessian else None


def compare_cached_predictors(model, p, rows):
    """Compatibility hook for the immutable synthetic reconstruction checker.

    Compares conditional outputs, not predictors. Never overrides model methods.
    """
    if not os.environ.get('SLURM_JOB_ID', '').isdigit():
        raise RuntimeError('Compute allocation required')
    import numpy as np
    rows = tuple(rows)
    started = time.perf_counter()
    prepared = PreparedConditional(model,p,rows)
    preparation_seconds = time.perf_counter()-started
    cache_layout = (tuple(prepared.rows), tuple(tuple(d['supports']) for d in prepared.rows.values()))
    original = prepared.original.conditional_terms
    cases,timings = [],[]
    for need_hessian in (False,True):
        inputs = [(row,*model.predictors(row,p,u,need_hessian=False)[:2])
                  for u in (-.5,0.,.5) for row in rows]
        started = time.perf_counter()
        references = [original(row,eta,r,need_hessian) for row,eta,r in inputs]
        unprepared_seconds = time.perf_counter()-started
        started = time.perf_counter()
        actuals = [prepared.terms(row,eta,r,p,need_hessian) for row,eta,r in inputs]
        prepared_seconds = time.perf_counter()-started
        for i,(ref,actual) in enumerate(zip(references,actuals)):
            errors = [None if a is None else float(np.max(np.abs(a-b))) for a,b in zip(ref,actual)]
            exact = [a is b if a is None else bool(np.array_equal(a,b)) for a,b in zip(ref,actual)]
            cases.append(dict(ordinal=i,need_hessian=need_hessian,errors=errors,exact_equal=exact))
            if any(e is not None and e > 1e-10 for e in errors):
                raise ValueError('Prepared conditional differs beyond fixed1e-10 check allowance')
        timings.append(dict(need_hessian=need_hessian,unprepared_seconds=unprepared_seconds,
                            prepared_seconds=prepared_seconds,preparation_seconds=preparation_seconds,
                            rows_times_u=len(inputs),single_run_not_end_to_end_speed=True))
    changed = np.asarray(p).copy()
    changed[0] += .01
    rejected = False
    try:
        row,eta,r = inputs[0]
        prepared.terms(row,eta,r,changed)
    except ValueError:
        rejected = True
    if not rejected:
        raise ValueError('Parameter mutation was not refused')
    # Valid context remains usable after refusal; no original globals were changed.
    prepared.terms(row,eta,r,p)
    if prepared.original.conditional_terms is not original:
        raise ValueError('Original conditional callable changed')
    if cache_layout != (tuple(prepared.rows), tuple(tuple(d['supports']) for d in prepared.rows.values())):
        raise ValueError('Per-call cache grew during eta/support evaluations')
    return dict(conditional_comparisons=cases,timings=timings,preparation_seconds=preparation_seconds,
        changed_parameters_rejected=rejected,valid_context_after_refusal=True,original_callable_unchanged=True,
        persistent_cache_layout_unchanged=True,eta_prefix_storage='local per-call only',
        rows=len(rows),fixed_u_points=[-.5,0.,.5],no_integration=True,no_fits=True,no_calibration=True)
