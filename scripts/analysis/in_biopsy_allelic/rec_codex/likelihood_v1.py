"""Independent Model A v1 arithmetic. Synthetic fixtures only until Gate B.

Parameter order: alpha[G], r[G], kappa, delta[Z], omega_S, rho_S,b,d.
Variance v is a separate final coordinate for the marginal model.
"""
from dataclasses import dataclass, replace
from functools import lru_cache

import numpy as np
from scipy.optimize import minimize
from scipy.special import betaln, digamma, expit, gammaln, logsumexp, polygamma
from scipy.stats import binom

EPS = 1e-6


def bb_logmass(k, n, aa, bb):
    """Full constants at each run's first count, recurrence within runs."""
    k = np.atleast_1d(np.asarray(k, dtype=float))
    mass = np.empty(len(k))
    if not len(k):
        return mass
    breaks = np.r_[0, np.flatnonzero(np.diff(k) != 1)+1, len(k)]
    for start, stop in zip(breaks[:-1], breaks[1:]):
        first = k[start]
        choose=gammaln(n+1)-gammaln(first+1)-gammaln(n-first+1)
        t=aa+bb
        if t>1e4:
            # Exact rising-factorial identity avoids subtracting large lbeta
            # values near the dispersion floor. This remains beta-binomial.
            left=int(first);right=int(n-first)
            mass[start]=(choose+left*np.log(aa/t)+right*np.log(bb/t)
                         +np.log1p(np.arange(left)/aa).sum()
                         +np.log1p(np.arange(right)/bb).sum()
                         -np.log1p(np.arange(n)/t).sum())
        else:
            mass[start]=choose+betaln(first+aa,n-first+bb)-betaln(aa,bb)
        counts = k[start:stop-1]
        increments = (np.log(n-counts)+np.log(counts+aa)
                      -np.log(counts+1)-np.log(n-counts-1+bb))
        mass[start+1:stop] = mass[start]+np.cumsum(increments)
    return mass


def allowed_counts(n, depth, minor, percent):
    if any(int(x) != x for x in (n, depth, minor, percent)) or not 0 <= percent <= 50:
        raise ValueError('integer het predicate required')
    lo = max(minor, (percent * n + 99) // 100)
    return np.arange(lo, n-lo+1, dtype=int) if n >= depth else np.array([], dtype=int)


@dataclass(frozen=True)
class Row:
    gene: int
    individual: str
    cohort: str
    n: int
    a: int
    p_ref: float
    p_alt: float
    p_x: float
    e: float
    phi: float
    omega: float
    orientation: int
    stage: float
    z: tuple
    background: float
    duplication: float
    depth: int
    minor: int
    percent: int

    @property
    def allowed(self):
        return allowed_counts(self.n, self.depth, self.minor, self.percent)

    @property
    def prior(self):
        r, a, x = self.p_ref, self.p_alt, self.p_x
        prior = np.array([2*r*a, r*r+2*r*x, a*a+2*a*x])
        return prior if self.orientation == 1 else prior[[0, 2, 1]]


def bb_terms(k, n, eta, linear_r, need_hessian=True):
    """Full BB log mass, analytic gradient and Hessian in eta and linear_r."""
    k = np.atleast_1d(np.asarray(k, dtype=float))
    m = expit(eta)
    one_m = expit(-eta)
    v = expit(linear_r)
    rho = EPS + (1-2*EPS)*v
    rp = (1-2*EPS)*v*(1-v)
    rpp = rp*(1-2*v)
    t = (1-rho)/rho
    tp = -rp/rho**2
    tpp = -rpp/rho**2 + 2*rp**2/rho**3
    aa, bb = m*t, one_m*t
    lp = bb_logmass(k, n, aa, bb)
    ga = digamma(k+aa)-digamma(n+t)-digamma(aa)+digamma(t)
    gb = digamma(n-k+bb)-digamma(n+t)-digamma(bb)+digamma(t)
    mp = m*one_m
    ja = np.array([t*mp, m*tp])
    jb = np.array([-t*mp, one_m*tp])
    grad = ga[:,None]*ja + gb[:,None]*jb
    if not need_hessian:
        return lp, grad, None
    hab = -polygamma(1,n+t)+polygamma(1,t)
    haa = polygamma(1,k+aa)-polygamma(1,aa)+hab
    hbb = polygamma(1,n-k+bb)-polygamma(1,bb)+hab
    hess = (haa[:,None,None]*np.outer(ja,ja)+hbb[:,None,None]*np.outer(jb,jb)
            + hab*(np.outer(ja,jb)+np.outer(jb,ja)))
    a2 = np.array([[t*mp*(1-2*m), mp*tp],[mp*tp,m*tpp]])
    b2 = np.array([[-t*mp*(1-2*m),-mp*tp],[-mp*tp,one_m*tpp]])
    hess += ga[:,None,None]*a2 + gb[:,None,None]*b2
    return lp, grad, hess


def fixed_bb(k, n, mean, rho):
    if rho == 0:
        return binom.logpmf(k,n,mean)
    t = (1-rho)/rho
    return bb_logmass(k, n, mean*t, (1-mean)*t)


def mixture_terms(row, k, eta, linear_r, need_hessian=True):
    k = np.asarray(k)
    lp, grad, hess = bb_terms(k,row.n,eta-row.orientation*row.omega,linear_r,need_hessian)
    with np.errstate(divide='ignore'):
        logprior = np.log(row.prior)
    parts = np.stack((lp+logprior[0],
                     fixed_bb(k,row.n,row.e,row.phi)+logprior[1],
                     fixed_bb(k,row.n,1-row.e,row.phi)+logprior[2]))
    mass = logsumexp(parts,axis=0)
    wh = np.exp(parts[0]-mass)
    outgrad = wh[:,None]*grad
    if not need_hessian:
        return mass,outgrad,None
    outhess = wh[:,None,None]*hess + (wh*(1-wh))[:,None,None]*np.einsum('ni,nj->nij',grad,grad)
    return mass,outgrad,outhess


def sum_terms(lp, score, hess):
    logden = logsumexp(lp)
    prob = np.exp(lp-logden)
    mean_score = prob @ score
    if hess is None:
        return float(logden), mean_score, None
    den_hess = (np.einsum('n,nij->ij',prob,hess)
                + np.einsum('n,ni,nj->ij',prob,score,score)-np.outer(mean_score,mean_score))
    return float(logden), mean_score, den_hess


@lru_cache(maxsize=8192)
def hom_log_inclusion(n, lo, e, phi):
    # Symmetric K gives exactly the same mass for the two hom classes.
    return float(logsumexp(fixed_bb(np.arange(lo, n-lo+1), n, e, phi)))


def het_inclusion(n, lo, eta, linear_r, need_hessian=True):
    tails = np.r_[np.arange(lo), np.arange(n-lo+1, n+1)]
    if not len(tails):
        return 0., np.zeros(2), np.zeros((2,2)) if need_hessian else None
    logtail, score, hess = sum_terms(*bb_terms(tails, n, eta, linear_r,need_hessian))
    if logtail < np.log(.5):
        logc = np.log1p(-np.exp(logtail))
        ratio = np.exp(logtail-logc)
        c_score = -ratio*score
        if not need_hessian:
            return float(logc), c_score, None
        c_hess = -ratio*(hess+np.outer(score,score))-np.outer(c_score,c_score)
        return float(logc), c_score, c_hess
    return sum_terms(*bb_terms(np.arange(lo, n-lo+1), n, eta, linear_r,need_hessian))


def conditional_terms(row, eta, linear_r, need_hessian=True):
    lo = max(row.minor, (row.percent*row.n+99)//100)
    if row.n < row.depth or not lo <= row.a <= row.n-lo:
        raise ValueError('row fails integer het predicate')
    logc, score, hess = het_inclusion(row.n, lo, eta-row.orientation*row.omega, linear_r,need_hessian)
    hom = hom_log_inclusion(row.n, lo, row.e, row.phi)
    with np.errstate(divide='ignore'):
        logs = np.log(row.prior)+np.array([logc,hom,hom])
    logden, ds, dh = sum_terms(logs, np.vstack((score,np.zeros((2,2)))),
                              np.concatenate((hess[None,:,:],np.zeros((2,2,2)))) if need_hessian else None)
    lp, ns, nh = mixture_terms(row,np.array([row.a]),eta,linear_r,need_hessian)
    return float(lp[0]-logden), ns[0]-ds, nh[0]-dh if need_hessian else None


class Model:
    def __init__(self, rows, gene_count=None):
        self.rows = tuple(sorted(rows,key=lambda r:(r.individual,r.gene)))
        self.G = gene_count if gene_count is not None else max(r.gene for r in rows)+1
        self.Z = len(rows[0].z)
        if any(r.gene not in range(self.G) for r in rows):
            raise ValueError('dense gene index required')
        self.size = 2*self.G+self.Z+5

    def initial(self):
        p = np.zeros(self.size)
        p[self.G:2*self.G] = np.log(.02/.98)
        return p

    def bounds(self):
        return [(-6,6)]*self.G+[(-15,5)]*self.G+[(-.9,.9)]+[(None,None)]*self.Z+[(-3,3)]*4

    def predictors(self,row,p,u=0,need_hessian=True):
        g, G, Z = row.gene,self.G,self.Z
        k = 2*G
        factor = 1+p[k]*row.stage+np.dot(p[k+1:k+1+Z],row.z)
        eta = p[g]*factor-row.orientation*p[k+1+Z]*row.stage+u*row.stage
        r = p[G+g]+np.dot(p[-3:],(row.stage,row.background,row.duplication))
        jac = np.zeros((2,self.size))
        jac[0,g]=factor
        jac[0,k]=p[g]*row.stage
        jac[0,k+1:k+1+Z]=p[g]*np.asarray(row.z)
        jac[0,k+1+Z]=-row.orientation*row.stage
        jac[1,G+g]=1
        jac[1,-3:]=[row.stage,row.background,row.duplication]
        eta_hess=None
        if need_hessian:
            eta_hess=np.zeros((self.size,self.size))
            eta_hess[g,k]=eta_hess[k,g]=row.stage
            eta_hess[g,k+1:k+1+Z]=row.z
            eta_hess[k+1:k+1+Z,g]=row.z
        return eta,r,jac,eta_hess

    def evaluate(self,p,rows=None,u=0,need_hessian=True):
        ll=0.; grad=np.zeros(self.size); hess=np.zeros((self.size,self.size)) if need_hessian else None
        for row in self.rows if rows is None else rows:
            eta,r,jac,eta_hess=self.predictors(row,p,u,need_hessian=need_hessian)
            val,s,h=conditional_terms(row,eta,r,need_hessian)
            ll+=val;grad+=jac.T@s
            if need_hessian:hess+=jac.T@h@jac+s[0]*eta_hess
        return ll,grad,hess

    def cold_start(self):
        start=self.initial()
        for gene in range(self.G):
            rs=[r for r in self.rows if r.gene==gene]
            if not rs:continue
            # The fixed cohort offset belongs in the starting likelihood.
            def intercept_objective(v):
                values=[conditional_terms(r,v[0],v[1],False) for r in rs]
                return -sum(t[0] for t in values),-np.sum([t[1] for t in values],axis=0)
            candidates=[minimize(intercept_objective,[a,np.log(.02/.98)],jac=True,
                                 method='L-BFGS-B',bounds=[(-6,6),(-15,5)]) for a in (-1.,0.,1.)]
            best=min(candidates,key=lambda f:f.fun)
            start[gene],start[self.G+gene]=best.x
        return start

    def informative(self):
        keep=np.ones(self.size,dtype=bool)
        present={r.gene for r in self.rows}
        for gene in set(range(self.G))-present:keep[gene]=keep[self.G+gene]=False
        for j in range(self.Z):
            if all(r.z[j]==0 for r in self.rows):keep[2*self.G+1+j]=False
        return keep

    def fit(self,start=None):
        cold=self.cold_start()
        if start is None:
            fit=self.fit_single(cold)
            return {**fit,'kept_start':'cold','start_logliks':{'cold':fit['loglik']}}
        warm=np.asarray(start,float)
        cold[~self.informative()]=warm[~self.informative()]
        fits={'warm':self.fit_single(warm),'cold':self.fit_single(cold)}
        kept='cold' if fits['cold']['loglik']-fits['warm']['loglik']>=1e-6 else 'warm'
        return {**fits[kept],'kept_start':kept,
                'start_logliks':{k:f['loglik'] for k,f in fits.items()},
                'start_converged':{k:f['converged'] for k,f in fits.items()}}

    def fit_single(self,start):
        start=np.asarray(start,float)
        informative=self.informative()
        bounds=[b if informative[j] else (start[j],start[j]) for j,b in enumerate(self.bounds())]
        def objective(p):
            ll,g,_=self.evaluate(p,need_hessian=False)
            return -ll,-g
        candidates=[]
        for shift in (0.,-.1,.1):
            p=start.copy()
            p[2*self.G]=np.clip(p[2*self.G]+shift,-.9,.9)
            fit=minimize(objective,p,jac=True,method='L-BFGS-B',bounds=bounds,
                         options={'ftol':1e-13,'gtol':1e-7,'maxiter':1000})
            candidates.append(fit)
        fit=min(candidates,key=lambda f:f.fun)
        p=fit.x.copy()
        # Damped Newton polish on interior coordinates.
        for _ in range(12):
            ll,g,h=self.evaluate(p)
            active=np.array([(lo is not None and lo==hi) or (lo is not None and p[j]<=lo+1e-8 and g[j]<0) or
                             (hi is not None and p[j]>=hi-1e-8 and g[j]>0)
                             for j,(lo,hi) in enumerate(bounds)])
            free=~active
            if np.max(np.abs(g[free]),initial=0)<1e-6:break
            step=np.zeros(self.size)
            step[free]=np.linalg.lstsq(-h[np.ix_(free,free)],g[free],rcond=1e-10)[0]
            improved=False
            for damping in 2.**-np.arange(15):
                trial=p+damping*step
                for j,(lo,hi) in enumerate(bounds):trial[j]=np.clip(trial[j],-np.inf if lo is None else lo,np.inf if hi is None else hi)
                if self.evaluate(trial,need_hessian=False)[0]>ll:p=trial;improved=True;break
            if not improved:break
        ll,g,h=self.evaluate(p)
        projected=g.copy()
        for j,(lo,hi) in enumerate(bounds):
            if (lo is not None and lo==hi) or (lo is not None and p[j]<=lo+1e-8 and g[j]<0) or (hi is not None and p[j]>=hi-1e-8 and g[j]>0):projected[j]=0
        if np.max(np.abs(projected))>=1e-6:
            p,_=self.constrained_polish(p,bounds)
            ll,g,h=self.evaluate(p)
            projected=g.copy()
            for j,(lo,hi) in enumerate(bounds):
                if (lo is not None and lo==hi) or (lo is not None and p[j]<=lo+1e-8 and g[j]<0) or (hi is not None and p[j]>=hi-1e-8 and g[j]>0):projected[j]=0
        return {'parameters':p,'loglik':ll,'gradient':g,'hessian':h,
                'converged':bool(np.max(np.abs(projected))<1e-6)}

    def constrained_polish(self,start,bounds,maxiter=100):
        """Solve a bound-constrained Newton subproblem when clipping stalls.

        Direction regularization changes only the proposal. Acceptance uses
        the exact likelihood, original bounds and original convergence rule.
        """
        p=np.asarray(start,float).copy();history=[]
        for iteration in range(maxiter):
            ll,g,h=self.evaluate(p)
            active=np.array([(lo is not None and lo==hi) or
                             (lo is not None and p[j]<=lo+1e-8 and g[j]<0) or
                             (hi is not None and p[j]>=hi-1e-8 and g[j]>0)
                             for j,(lo,hi) in enumerate(bounds)])
            free=np.flatnonzero(~active)
            maximum=float(np.max(np.abs(g[free]),initial=0))
            history.append({'iteration':iteration,'loglik':ll,'max_projected_gradient':maximum})
            if maximum<1e-6:break
            curvature=-h[np.ix_(free,free)]
            values,vectors=np.linalg.eigh((curvature+curvature.T)/2)
            floor=max(1e-9,float(np.max(np.abs(values)))*1e-10)
            curvature=(vectors*np.maximum(values,floor))@vectors.T
            score=g[free]
            limits=[(None if bounds[j][0] is None else bounds[j][0]-p[j],
                     None if bounds[j][1] is None else bounds[j][1]-p[j]) for j in free]
            def quadratic(step):
                return .5*step@curvature@step-score@step,curvature@step-score
            direction=minimize(quadratic,np.zeros(len(free)),jac=True,method='SLSQP',bounds=limits,
                               options={'ftol':1e-14,'maxiter':1000})
            step=np.zeros(self.size);step[free]=direction.x
            improved=False
            for damping in 2.**-np.arange(24):
                trial=p+damping*step
                for j,(lo,hi) in enumerate(bounds):
                    trial[j]=np.clip(trial[j],-np.inf if lo is None else lo,np.inf if hi is None else hi)
                if self.evaluate(trial,need_hessian=False)[0]>ll:
                    p=trial;improved=True;break
            if not improved:break
        return p,history

    def sandwich(self,p,clusters=None):
        _,g,h=self.evaluate(p)
        scores={}
        for row in self.rows:
            cluster=clusters.get(row.individual,row.individual) if clusters else row.individual
            scores.setdefault(cluster,np.zeros(self.size))
            scores[cluster]+=self.evaluate(p,[row],need_hessian=False)[1]
        free=np.array([not ((lo is not None and p[j]<=lo+1e-8) or
                           (hi is not None and p[j]>=hi-1e-8)) for j,(lo,hi) in enumerate(self.bounds())])
        free &= self.informative()
        k=2*self.G
        free[k]=True
        ek=np.zeros(self.size);ek[k]=1
        direction=np.zeros(self.size)
        try:
            direction[free]=np.linalg.solve(-h[np.ix_(free,free)],ek[free])
        except np.linalg.LinAlgError:
            return float('nan')
        return float(np.sqrt(sum((direction@s)**2 for s in scores.values())))

    def marginal(self,p,v,nodes=20,with_hessian=False,gene_terms=None):
        if v<0:raise ValueError('negative variance')
        if v==0:
            ll,g,hp=self.evaluate(p,need_hessian=with_hessian);score_v=0.
            for gene in range(self.G):
                su=0.;hu=0.
                if gene_terms is not None:
                    gene_terms[gene]=self.evaluate(p,[r for r in self.rows if r.gene==gene],need_hessian=False)[0]
                for row in self.rows:
                    if row.gene!=gene:continue
                    eta,r,_,_=self.predictors(row,p,need_hessian=False)
                    _,s,h=conditional_terms(row,eta,r)
                    su+=s[0]*row.stage;hu+=h[0,0]*row.stage**2
                score_v+=.5*(su**2+hu)
            if with_hessian:
                # At a KKT variance boundary Newton uses nuisance coordinates
                # only. The v Hessian would require higher-order row derivatives.
                h=np.full((self.size+1,self.size+1),np.nan)
                h[:-1,:-1]=hp
                return ll,np.r_[g,score_v],h
            return ll,np.r_[g,score_v]
        roots,weights=np.polynomial.hermite.hermgauss(nodes)
        total=0.;gradient=np.zeros(self.size+1)
        hessian=np.zeros((self.size+1,self.size+1)) if with_hessian else None
        for gene in range(self.G):
            rs=[r for r in self.rows if r.gene==gene]
            if not rs:
                if gene_terms is not None:gene_terms[gene]=0.
                continue
            def hstats(u):
                ll=-u*u/(2*v);s=-u/v;h=-1/v
                for row in rs:
                    eta,r,_,_=self.predictors(row,p,u,need_hessian=False)
                    value,score,hh=conditional_terms(row,eta,r)
                    ll+=value;s+=score[0]*row.stage;h+=hh[0,0]*row.stage**2
                return ll,s,h
            mode=0.
            for _ in range(100):
                ll,s,h=hstats(mode)
                step=s/max(-h,1e-8)
                accepted=False
                for damping in 2.**-np.arange(30):
                    trial=mode+damping*step
                    if hstats(trial)[0]>=ll-1e-12:
                        mode=trial;accepted=True;break
                if not accepted or abs(damping*step)<1e-10:break
            curvature=hstats(mode)[2]
            if curvature>=0:raise ValueError('nonconcave mode in adaptive quadrature')
            scale=(-curvature)**-.5
            terms=[];scores=[];hessians=[]
            for x,w in zip(roots,weights):
                u=mode+np.sqrt(2)*scale*x
                ll,g,h=self.evaluate(p,rs,u,need_hessian=with_hessian)
                terms.append(np.log(w)+x*x+ll-u*u/(2*v))
                scores.append(np.r_[g,(u*u-v)/(2*v*v)])
                if with_hessian:
                    hh=np.zeros((self.size+1,self.size+1))
                    hh[:-1,:-1]=h
                    hh[-1,-1]=1/(2*v*v)-u*u/(v*v*v)
                    hessians.append(hh)
            logsum=logsumexp(terms)
            contribution=np.log(np.sqrt(2)*scale)+logsum-.5*np.log(2*np.pi*v)
            total+=contribution
            if gene_terms is not None:
                gene_terms[gene]=float(contribution)
            probability=np.exp(np.asarray(terms)-logsum)
            scores=np.asarray(scores)
            gs=probability@scores
            gradient+=gs
            if with_hessian:
                hessian+=(np.einsum('n,nij->ij',probability,np.asarray(hessians))
                          +np.einsum('n,ni,nj->ij',probability,scores,scores)-np.outer(gs,gs))
        if with_hessian:
            return float(total),gradient,hessian
        return float(total),gradient

    def fit_marginal(self,tau_max,start=None,fixed_v=None,null_fit=None,profile_start=None,
                     section4_from_null=False):
        null=self.fit(start) if null_fit is None else null_fit
        p0=np.asarray(null['parameters'],float)
        if fixed_v is not None:
            warm=np.r_[p0 if profile_start is None else profile_start,fixed_v]
            cold=np.r_[p0,fixed_v]
        else:
            warm=np.r_[p0,0.]
            cold=np.r_[p0 if section4_from_null else self.cold_start(),tau_max**2]
        cold[:-1][~self.informative()]=p0[~self.informative()]
        fits={'warm':self.fit_marginal_single(warm,fixed_v),
              'cold':self.fit_marginal_single(cold,fixed_v)}
        kept='cold' if fits['cold']['loglik']-fits['warm']['loglik']>=1e-6 else 'warm'
        fit=fits[kept]
        if fixed_v is None and fit['v']==0. and fit['gradient'][-1]>1e-6:
            restart=self.fit_marginal_single(np.r_[fit['parameters'],tau_max**2],None)
            if restart['loglik']>fit['loglik']:fit=restart;kept='boundary_restart'
        if fixed_v is None and fit['v']==0. and null['loglik']>=fit['loglik']-1e-6:
            ll,g=self.marginal(p0,0.)
            fit={**fit,'parameters':p0,'loglik':ll,'gradient':g,
                 'converged':bool(null['converged'] and g[-1]<=1e-6)}
        return {**fit,'null':null,'kept_start':kept,
                'start_logliks':{k:f['loglik'] for k,f in fits.items()},
                'start_converged':{k:f['converged'] for k,f in fits.items()}}

    def fit_marginal_single(self,start,fixed_v):
        informative=self.informative()
        bounds=[b if informative[j] else (start[j],start[j]) for j,b in enumerate(self.bounds())]
        bounds += [(0,None) if fixed_v is None else (fixed_v,fixed_v)]
        candidates=[]
        def objective(p):
            ll,g=self.marginal(p[:-1],p[-1])
            return -ll,-g
        for shift in (0.,-.1,.1):
            p=np.asarray(start,float).copy()
            p[2*self.G]=np.clip(p[2*self.G]+shift,-.9,.9)
            fit=minimize(objective,p,jac=True,method='L-BFGS-B',bounds=bounds,
                         options={'ftol':1e-13,'gtol':1e-7,'maxiter':500})
            candidates.append(fit)
        fit=min(candidates,key=lambda f:f.fun)
        p=fit.x.copy()
        for _ in range(20):
            ll,g,h=self.marginal(p[:-1],p[-1],with_hessian=True)
            active=np.array([(lo is not None and lo==hi) or
                             (lo is not None and p[j]<=lo+1e-10 and g[j]<0) or
                             (hi is not None and p[j]>=hi-1e-10 and g[j]>0)
                             for j,(lo,hi) in enumerate(bounds)])
            free=~active
            if np.max(np.abs(g[free]),initial=0)<1e-6:
                break
            hf=h[np.ix_(free,free)]
            if not np.isfinite(hf).all():
                # Positive boundary score needs an interior optimizer step.
                break
            step=np.zeros(len(p))
            step[free]=np.linalg.lstsq(-hf,g[free],rcond=1e-12)[0]
            if np.dot(step,g)<=0:
                break
            improved=False
            for damping in 2.**-np.arange(20):
                trial=p+damping*step
                for j,(lo,hi) in enumerate(bounds):
                    trial[j]=np.clip(trial[j],-np.inf if lo is None else lo,np.inf if hi is None else hi)
                if self.marginal(trial[:-1],trial[-1])[0]>ll:
                    p=trial;improved=True;break
            if not improved:
                break
        ll,g=self.marginal(p[:-1],p[-1])
        projected=g.copy()
        for j,(lo,hi) in enumerate(bounds):
            if (lo is not None and lo==hi) or (lo is not None and p[j]<=lo+1e-10 and g[j]<0) or (hi is not None and p[j]>=hi-1e-10 and g[j]>0):projected[j]=0
        return {'parameters':p[:-1],'v':p[-1],'loglik':ll,'gradient':g,
                'converged':bool(np.max(np.abs(projected))<1e-6)}


def conditional_draw(row,eta,linear_r,uniform_class,uniform_count):
    """P3's two uniforms: inclusion-reweighted genotype, then truncated count."""
    k=row.allowed
    lp=bb_terms(k,row.n,eta-row.orientation*row.omega,linear_r)[0]
    parts=np.stack((lp,fixed_bb(k,row.n,row.e,row.phi),fixed_bb(k,row.n,1-row.e,row.phi)))
    with np.errstate(divide='ignore'): class_logs=np.log(row.prior)+logsumexp(parts,axis=1)
    class_prob=np.exp(class_logs-logsumexp(class_logs))
    g=min(2,int(np.searchsorted(np.cumsum(class_prob),uniform_class,side='left')))
    p=np.exp(parts[g]-logsumexp(parts[g]))
    chosen=min(len(k)-1,int(np.searchsorted(np.cumsum(p),uniform_count,side='left')))
    return replace(row,a=int(k[chosen]))
