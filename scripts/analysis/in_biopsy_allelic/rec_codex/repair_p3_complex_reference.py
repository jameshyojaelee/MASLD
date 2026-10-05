"""Independent analytic complex probability reference; compute-only fixed points.

No original conditional evaluator, BB derivatives, candidate scores or Hessians.
Complex perturbations extend the exact integer rising-factorial probability law.
Quadrature estimates and floating arithmetic are not certified error bounds.
"""
import os
import warnings

if not os.environ.get('SLURM_JOB_ID', '').isdigit() or os.environ.get('SLURM_ARRAY_TASK_ID'):
    raise RuntimeError('Individual compute allocation required')

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize_scalar
from scipy.special import gammaln, logsumexp
from scipy.stats import norm

STEPS = (1e-10, 1e-15, 1e-20)


class ReferenceFailure(ValueError):
    def __init__(self, message, record=None):
        super().__init__(message)
        self.reference_record = record


def require(ok, message):
    if not ok:
        raise ReferenceFailure(message)


def logistic(x):
    x = complex(x)
    if x.real >= 0:
        return 1 / (1 + np.exp(-x))
    e = np.exp(x)
    return e / (1 + e)


def log_logistic(x):
    x = complex(x)
    return -np.log1p(np.exp(-x)) if x.real >= 0 else x - np.log1p(np.exp(x))


def rising_prefix(logshape, n):
    result = np.empty(n + 1, dtype=complex)
    result[0] = 0.
    if n:
        terms = np.empty(n, dtype=complex)
        terms[0] = logshape
        logj = np.log(np.arange(1, n, dtype=float))
        shift = np.maximum(complex(logshape).real, logj)
        terms[1:] = shift + np.log(np.exp(logshape-shift) + np.exp(logj-shift))
        np.cumsum(terms, out=result[1:])
    return result


def complex_logsumexp(values, axis=None):
    # A real scale is only a numerical identity, not a transformation of the law.
    shift = np.max(values.real, axis=axis, keepdims=True)
    result = shift + np.log(np.sum(np.exp(values-shift), axis=axis, keepdims=True))
    return result.squeeze() if axis is None else np.squeeze(result, axis=axis)


class BB:
    def __init__(self, n, linear_r=None, *, mean=None, phi=None):
        require(isinstance(n, (int, np.integer)) and 0 <= n <= 100000, 'Invalid integer support')
        self.n = int(n)
        k = np.arange(n + 1)
        self.choose = gammaln(n+1) - gammaln(k+1) - gammaln(n-k+1)
        if linear_r is not None:
            rho = 1e-6 + (1-2e-6)*logistic(linear_r)
            self.logt = np.log((1-rho)/rho)
            self.fixed = None
        else:
            require(0 < mean < 1 and 0 <= phi < 1, 'Invalid fixed homozygous law')
            self.fixed = (mean,phi)
            self.logt = np.log((1-phi)/phi) if phi else None
        self.denominator = rising_prefix(self.logt,self.n)[self.n] if self.logt is not None else None

    def logpmf(self, eta=None):
        if self.fixed is not None:
            mean,phi = self.fixed
            if phi == 0:
                k = np.arange(self.n+1)
                return (self.choose+k*np.log(mean)+(self.n-k)*np.log1p(-mean)).astype(complex)
            la,lb = self.logt+np.log(mean),self.logt+np.log1p(-mean)
        else:
            la,lb = self.logt+log_logistic(eta),self.logt+log_logistic(-eta)
        return self.choose+rising_prefix(la,self.n)+rising_prefix(lb,self.n)[::-1]-self.denominator


class Reference:
    def __init__(self, rows, gene_count, z_count, maximum_evaluations=100000):
        self.rows, self.G, self.Z = list(rows),gene_count,z_count
        self.size = 2*gene_count+z_count+5
        self.maximum_evaluations = maximum_evaluations
        self.evaluations = 0
        self.parts = []
        for row in self.rows:
            lo = max(row.minor,(row.percent*row.n+99)//100)
            support = np.arange(lo,row.n-lo+1)
            require(row.n >= row.depth and row.a in support and len(row.z)==z_count, 'Selected row identity differs')
            prior = np.array([2*row.p_ref*row.p_alt,row.p_ref**2+2*row.p_ref*row.p_x,row.p_alt**2+2*row.p_alt*row.p_x])
            if row.orientation == -1:
                prior = prior[[0,2,1]]
            require(np.all(prior>=0) and prior.sum()>0, 'Invalid raw three-class weights')
            lp = np.full(3,-np.inf,dtype=complex)
            nonzero = prior>0
            lp[nonzero] = np.log(prior[nonzero])
            self.parts.append((support,lp,BB(row.n,mean=row.e,phi=row.phi).logpmf()[support],
                               BB(row.n,mean=1-row.e,phi=row.phi).logpmf()[support],int(row.a-lo)))

    def prepare(self, q):
        require(len(q)==self.size and np.isfinite(q).all(), 'Invalid parameter vector')
        return [BB(row.n,q[self.G+row.gene]+np.dot(q[-3:],(row.stage,row.background,row.duplication))) for row in self.rows]

    def loglik(self,q,u,prepared):
        self.evaluations += 1
        require(self.evaluations<=self.maximum_evaluations, 'Independent complex probability evaluation cap exhausted')
        total = 0j
        for row,part,bb in zip(self.rows,self.parts,prepared):
            eta = q[row.gene]*(1+q[2*self.G]*row.stage+np.dot(q[2*self.G+1:2*self.G+1+self.Z],row.z))
            eta += u*row.stage-row.orientation*q[2*self.G+1+self.Z]*row.stage-row.orientation*row.omega
            support,prior,homref,homalt,at = part
            mixture = complex_logsumexp(np.stack((bb.logpmf(eta)[support],homref,homalt))+prior[:,None],axis=0)
            total += mixture[at]-complex_logsumexp(mixture)
        require(np.isfinite(total), 'Nonfinite complex selected probability')
        return total

    def score_bounds(self,p):
        # Explicit independent-u predictor Jacobians; no original derivative call.
        cp = np.zeros(self.size)
        for row in self.rows:
            eta_j,r_j = np.zeros(self.size),np.zeros(self.size)
            eta_j[row.gene] = 1+p[2*self.G]*row.stage+np.dot(p[2*self.G+1:2*self.G+1+self.Z],row.z)
            eta_j[2*self.G] = p[row.gene]*row.stage
            eta_j[2*self.G+1:2*self.G+1+self.Z] = p[row.gene]*np.asarray(row.z)
            eta_j[2*self.G+1+self.Z] = -row.orientation*row.stage
            r_j[self.G+row.gene] = 1
            r_j[-3:] = (row.stage,row.background,row.duplication)
            cp += 2*row.n*(np.abs(eta_j)+np.abs(r_j))
        return cp

    def geometry(self,p,variance,bound,spacing):
        prepared = self.prepare(p)
        cache = {}
        def h(z):
            if float(z) not in cache:
                value = self.loglik(p,np.sqrt(variance)*z,prepared)
                require(abs(value.imag)<1e-12, 'Real baseline produced imaginary probability')
                cache[float(z)] = value.real-z*z/2-.5*np.log(2*np.pi)
            return cache[float(z)]
        grid = np.linspace(-bound,bound,int(np.ceil(2*bound/spacing))+1)
        values = np.array([h(z) for z in grid])
        points = set(map(float,grid[1:-1]))
        centers = [float(grid[values.argmax()])]
        for i in range(1,len(grid)-1):
            if values[i]>=max(values[i-1],values[i+1]):
                peak = minimize_scalar(lambda z:-h(z),bounds=(grid[i-1],grid[i+1]),method='bounded',options={'xatol':1e-10})
                require(peak.success,'Independent real integration-center search failed')
                centers.append(float(peak.x))
                points.update(peak.x+a for a in (-.5,-.1,0.,.1,.5) if -bound<peak.x+a<bound)
        return dict(bound=bound,spacing=spacing,points=sorted(points),shift=max(h(c) for c in centers),
                    centers=centers,mode_search_not_exhaustive=True)

    def integrate(self,function,geometry,kind):
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter('always')
            answer = quad(function,-geometry['bound'],geometry['bound'],points=geometry['points'],limit=4000,
                          epsabs=2e-12,epsrel=2e-12,full_output=True)
        value,error,info = answer[:3]
        record = dict(kind=kind,value=float(value),estimated_error=float(error),evaluations=int(info['neval']),
                      warnings=[str(w.message) for w in caught],message=answer[3] if len(answer)>3 else None,
                      reported_quad_error_not_certified=True)
        if len(answer)!=3 or not np.isfinite(value) or not np.isfinite(error) or error<0:
            raise ReferenceFailure('Independent real scalar quadrature refused',record)
        return record

    def mass(self,p,variance,geometry):
        prepared = self.prepare(p)
        def function(z):
            return float(np.exp(self.loglik(p,np.sqrt(variance)*z,prepared).real-z*z/2-.5*np.log(2*np.pi)-geometry['shift']))
        record = self.integrate(function,geometry,'real_mass')
        require(record['value']>record['estimated_error'], 'Independent mass lower estimate invalid')
        record.update(logmass=geometry['shift']+np.log(record['value']),
                      logmass_lower_estimate=geometry['shift']+np.log(record['value']-record['estimated_error']),
                      arithmetic_roundoff_not_bounded=True,geometry=geometry)
        return record

    def derivative(self,p,variance,axis,step,geometry,mass):
        require(axis in range(self.size+1) and step in STEPS,'Undeclared complex derivative')
        q = np.asarray(p,dtype=complex).copy()
        cv = complex(variance)
        if axis==self.size:
            cv += 1j*step
        else:
            q[axis] += 1j*step
        prepared = self.prepare(q)
        rootv = np.sqrt(cv)
        def function(z):
            exponent = self.loglik(q,rootv*z,prepared)-z*z/2-.5*np.log(2*np.pi)-geometry['shift']
            # Scale before quadrature so its absolute tolerance applies to the derivative.
            return float(np.exp(exponent).imag/step)
        numerator = self.integrate(function,geometry,'imaginary_numerator_divided_by_step')
        m,me = mass['value'],mass['estimated_error']
        gradient = numerator['value']/m
        quadrature_error = (numerator['estimated_error']+abs(gradient)*me)/(m-me)
        b,lower = geometry['bound'],mass['logmass_lower_estimate']
        logm0 = np.log(2.)+norm.logsf(b)
        ratio = float(np.exp(logm0-lower))
        if axis==self.size:
            logm2 = np.log(2.)+logsumexp([np.log(b)+norm.logpdf(b),norm.logsf(b)])
            abs_tail = float(np.exp(logsumexp([logm2,logm0])-np.log(2*variance)-lower))
            endpoint = float(np.exp(np.log(b)+norm.logpdf(b)-np.log(variance)-lower))
        else:
            abs_tail = float(self.score_bounds(p)[axis]*ratio)
            endpoint = 0.
        tail = abs_tail+abs(gradient)*ratio
        return dict(axis=axis,step=step,gradient=float(gradient),numerator=numerator,real_mass=mass['value'],
            estimated_numerator_denominator_gradient_error=float(quadrature_error),
            relative_prior_mass_tail_allowance=ratio,relative_derivative_tail_allowance=tail,
            finite_domain_variance_endpoint_allowance=endpoint,endpoint_term_not_subtracted=True,
            estimated_combined_gradient_uncertainty=float(quadrature_error+tail+endpoint),
            complex_step_truncation_and_floating_arithmetic_not_bounded=True,
            conditional_on_estimated_lower_mass=True,uncertified_estimate=True,
            complex_sqrt_variance_pathwise_derivative=axis==self.size)
