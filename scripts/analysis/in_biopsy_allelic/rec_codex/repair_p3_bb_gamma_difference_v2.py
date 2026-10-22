"""Isolated native BB arithmetic using stable gamma differences.

Positive shapes >=32 use eight Bernoulli terms of the Stirling expansion.
Smaller shapes separate the j=0 term before special-function evaluation.
This avoids building all n+1 masses when a subset is requested. Mathematical
series remainders are not a floating arithmetic accuracy certificate.
Sources: https://dlmf.nist.gov/5.11 and https://dlmf.nist.gov/5.15
No integration, fitter, bound, genotype mixture or selection rule is changed.
"""
import os

BERNOULLI = (1/6,-1/30,1/42,-1/30,5/66,-691/2730,7/6,-3617/510)
MAXIMUM_N = 5000


def require(ok,message):
    if not ok:
        raise ValueError(message)


def rising(counts,logshape,need_second=True):
    """log rising product and its first two log-shape derivatives."""
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    import numpy as np
    from scipy.special import digamma,gammaln,polygamma
    counts = np.atleast_1d(np.asarray(counts,dtype=float))
    require(counts.ndim == 1 and np.isfinite(counts).all() and (counts >= 0).all()
            and (counts <= MAXIMUM_N).all() and (counts == np.floor(counts)).all()
            and np.isfinite(logshape) and logshape <= 16,'Invalid bounded rising-product inputs')
    with np.errstate(under='ignore'):
        shape = float(np.exp(logshape))
    if shape >= 32:
        x = counts/shape
        l = np.log1p(x)
        f = (1+x)*l-x
        g = l-x/(1+x) if need_second else None
        f_series = np.zeros_like(x)
        g_series = np.zeros_like(x) if need_second else None
        for j in range(8,1,-1):
            f_series = f_series*x+(-1)**j/(j*(j-1))
            if need_second:
                g_series = g_series*x+(-1)**j*(j-1)/j
        f = np.where(x <= 1e-3,x*x*f_series,f)
        if need_second:
            g = np.where(x <= 1e-3,x*x*g_series,g)
        correction = shape*f-.5*l
        first = shape*l+.5*x/(1+x)
        second = shape*g+.5*(x/(1+x)+np.expm1(-2*l)) if need_second else None
        for j,bernoulli in enumerate(BERNOULLI,1):
            power = 2*j-1
            correction += bernoulli/(2*j*power)*shape**(-power)*np.expm1(-power*l)
            first += bernoulli/(2*j)*shape**(1-2*j)*(-np.expm1(-2*j*l))
            if need_second:
                second += bernoulli*shape**(1-2*j)*(
                    -np.expm1(-2*j*l)/(2*j)+np.expm1(-(2*j+1)*l))
        logs = counts*logshape+correction
    else:
        logs,first = np.zeros_like(counts),np.zeros_like(counts)
        second = np.zeros_like(counts) if need_second else None
        present = counts > 0
        s = counts[present]
        # This remains defined if exp(logshape) underflows: j=0 is retained
        # as logshape, with exact log-shape derivatives 1 and 0.
        logs[present] = logshape+gammaln(shape+s)-gammaln(shape+1)
        rest = shape*(digamma(shape+s)-digamma(shape+1))
        first[present] = 1+rest
        if need_second:
            second[present] = rest-shape*shape*(polygamma(1,shape+1)-polygamma(1,shape+s))
    logs = np.where(counts == 0,0.,np.where(counts == 1,logshape,logs))
    first = np.where(counts == 0,0.,np.where(counts == 1,1.,first))
    if need_second:
        second = np.where(counts <= 1,0.,second)
    require(all(v is None or np.isfinite(v).all() for v in (logs,first,second)), 'Nonfinite rising product')
    return logs,first,second


def stable_bb_terms(k,n,eta,linear_r,need_hessian=True):
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    import numpy as np
    from scipy.special import expit,gammaln
    k = np.atleast_1d(np.asarray(k))
    require(np.isfinite(n) and int(n) == n and 0 <= n <= MAXIMUM_N and k.ndim == 1
            and np.isfinite(k).all() and (k == np.floor(k)).all() and ((k >= 0)&(k <= n)).all()
            and np.isfinite(eta) and np.isfinite(linear_r),'Invalid bounded BB inputs')
    n,k = int(n),k.astype(int)
    if n == 0:
        return np.zeros(len(k)),np.zeros((len(k),2)),np.zeros((len(k),2,2)) if need_hessian else None
    m,q = expit(eta),expit(-eta)
    vr,qr = expit(linear_r),expit(-linear_r)
    rho = 1e-6+(1-2e-6)*vr
    rp = (1-2e-6)*vr*qr
    rpp = rp*(qr-vr)
    logt = np.log1p(-rho)-np.log(rho)
    loga,logb = logt-np.logaddexp(0.,-eta),logt-np.logaddexp(0.,eta)
    ra,da,qa = rising(k,loga,need_hessian)
    rb,db,qb = rising(n-k,logb,need_hessian)
    rt_array,dt_array,qt_array = rising([n],logt,need_hessian)
    rt,dt = rt_array[0],dt_array[0]
    qt = qt_array[0] if need_hessian else None
    choose = gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)
    lp = choose+ra+rb-rt
    denominator = rho*(1-rho)
    u = -rp/denominator
    w = -rpp/denominator+rp*rp*(1-2*rho)/(denominator*denominator)
    score = np.column_stack((da*q-db*m,(da+db-dt)*u))
    hessian = None
    if need_hessian:
        hessian = np.empty((len(k),2,2))
        hessian[:,0,0] = qa*q*q+qb*m*m-(da+db)*m*q
        hessian[:,0,1] = hessian[:,1,0] = u*(qa*q-qb*m)
        hessian[:,1,1] = u*u*(qa+qb-qt)+w*(da+db-dt)
    require(np.isfinite(lp).all() and np.isfinite(score).all()
            and (hessian is None or np.isfinite(hessian).all()),'Nonfinite candidate BB result')
    return lp,score,hessian

