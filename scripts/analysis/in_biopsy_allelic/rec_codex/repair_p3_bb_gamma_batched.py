"""Bounded node-by-count version of the checked Gamma v3 BB arithmetic.

No conditional-mixture assembly, integration, optimization or adoption. Shapes
and integer counts use the same formulas and branches as Gamma v3. Comparison
with its scalar outputs and the independent saved references is required.
"""
import hashlib
import os
from pathlib import Path

SCALAR_SHA='546e7ea2f1db21e6fcdd6aed72204d3bbe48c2d818996ab21de5f5ab5215bcf3'
MAXIMUM_ENTRIES=100000
MAXIMUM_NODES=64
MAXIMUM_N=5000
BERNOULLI=(1/6,-1/30,1/42,-1/30)


def require(ok,message):
    if not ok:
        raise RuntimeError(message)


def compute_guard():
    require(os.environ.get('SLURM_JOB_ID','').isdigit()
            and not os.environ.get('SLURM_ARRAY_TASK_ID'),'Individual compute allocation required')


def rising(counts,logshapes,need_second=True):
    compute_guard()
    import numpy as np
    from scipy.special import digamma,gammaln,polygamma
    counts=np.atleast_1d(np.asarray(counts,dtype=float))
    logshapes=np.atleast_1d(np.asarray(logshapes,dtype=float))
    require(counts.ndim==1 and logshapes.ndim==1 and 0<len(logshapes)<=MAXIMUM_NODES
            and len(counts)*len(logshapes)<=MAXIMUM_ENTRIES
            and np.isfinite(counts).all() and np.all((counts>=0)&(counts<=MAXIMUM_N))
            and np.all(counts==np.floor(counts)) and np.isfinite(logshapes).all()
            and np.all(logshapes<=16),'Invalid bounded rising-product batch')
    with np.errstate(under='ignore'):
        shapes=np.exp(logshapes)
    logs=np.zeros((len(shapes),len(counts)))
    first=np.zeros_like(logs)
    second=np.zeros_like(logs) if need_second else None
    large=np.flatnonzero(shapes>=32)
    if len(large):
        shape=shapes[large,None]
        x=counts[None,:]/shape
        l=np.log1p(x)
        f=(1+x)*l-x
        g=l-x/(1+x) if need_second else None
        fs=np.zeros_like(x);gs=np.zeros_like(x) if need_second else None
        for j in range(8,1,-1):
            fs=fs*x+(-1)**j/(j*(j-1))
            if need_second:gs=gs*x+(-1)**j*(j-1)/j
        f=np.where(x<=1e-3,x*x*fs,f)
        if need_second:g=np.where(x<=1e-3,x*x*gs,g)
        correction=shape*f-.5*l
        d=shape*l+.5*x/(1+x)
        q=shape*g+.5*(x/(1+x)+np.expm1(-2*l)) if need_second else None
        for j,b in enumerate(BERNOULLI,1):
            power=2*j-1
            correction+=b/(2*j*power)*shape**(-power)*np.expm1(-power*l)
            d+=b/(2*j)*shape**(1-2*j)*(-np.expm1(-2*j*l))
            if need_second:q+=b*shape**(1-2*j)*(-np.expm1(-2*j*l)/(2*j)+np.expm1(-(2*j+1)*l))
        logs[large,:]=counts[None,:]*logshapes[large,None]+correction
        first[large,:]=d
        if need_second:second[large,:]=q
    small=np.flatnonzero(shapes<32)
    positive=np.flatnonzero(counts>0)
    if len(small) and len(positive):
        a=shapes[small,None];s=counts[None,positive]
        index=np.ix_(small,positive)
        logs[index]=logshapes[small,None]+gammaln(a+s)-gammaln(a+1)
        rest=a*(digamma(a+s)-digamma(a+1))
        first[index]=1+rest
        if need_second:second[index]=rest-a*a*(polygamma(1,a+1)-polygamma(1,a+s))
    zero=counts==0;one=counts==1
    logs[:,zero]=0.;logs[:,one]=logshapes[:,None]
    first[:,zero]=0.;first[:,one]=1.
    if need_second:second[:,counts<=1]=0.
    require(all(x is None or np.isfinite(x).all() for x in (logs,first,second)),
            'Nonfinite rising-product batch')
    return logs,first,second


def bb_terms(k,n,etas,linear_r,scalar,need_hessian=True):
    compute_guard()
    import numpy as np
    from scipy.special import expit,gammaln
    require(hashlib.sha256(Path(scalar.__file__).read_bytes()).hexdigest()==SCALAR_SHA,
            'Checked scalar Gamma source differs')
    k=np.atleast_1d(np.asarray(k));eta=np.atleast_1d(np.asarray(etas,dtype=float))
    require(np.isfinite(n) and int(n)==n and 0<=n<=MAXIMUM_N and k.ndim==1
            and eta.ndim==1 and 0<len(eta)<=MAXIMUM_NODES
            and len(k)*len(eta)<=MAXIMUM_ENTRIES and np.isfinite(k).all()
            and np.all(k==np.floor(k)) and np.all((k>=0)&(k<=n))
            and np.isfinite(eta).all() and np.asarray(linear_r).shape==()
            and np.isfinite(linear_r),'Invalid bounded BB batch')
    n,k=int(n),k.astype(int)
    if n==0:
        shape=(len(eta),len(k))
        return np.zeros(shape),np.zeros((*shape,2)),np.zeros((*shape,2,2)) if need_hessian else None
    m,q=expit(eta)[:,None],expit(-eta)[:,None]
    vr,qr=expit(linear_r),expit(-linear_r)
    rho=1e-6+(1-2e-6)*vr;rp=(1-2e-6)*vr*qr;rpp=rp*(qr-vr)
    logt=np.log1p(-rho)-np.log(rho)
    ra,da,qa=rising(k,logt-np.logaddexp(0.,-eta),need_hessian)
    rb,db,qb=rising(n-k,logt-np.logaddexp(0.,eta),need_hessian)
    rt,dt,qt=scalar.rising([n],logt,need_hessian)
    choose=gammaln(n+1)-gammaln(k+1)-gammaln(n-k+1)
    lp=choose[None,:]+ra+rb-rt[0]
    denominator=rho*(1-rho);u=-rp/denominator
    w=-rpp/denominator+rp*rp*(1-2*rho)/(denominator*denominator)
    score=np.stack((da*q-db*m,(da+db-dt[0])*u),axis=2)
    hessian=None
    if need_hessian:
        hessian=np.empty((*lp.shape,2,2))
        hessian[:,:,0,0]=qa*q*q+qb*m*m-(da+db)*m*q
        hessian[:,:,0,1]=hessian[:,:,1,0]=u*(qa*q-qb*m)
        hessian[:,:,1,1]=u*u*(qa+qb-qt[0])+w*(da+db-dt[0])
    require(np.isfinite(lp).all() and np.isfinite(score).all()
            and (hessian is None or np.isfinite(hessian).all()),'Nonfinite BB batch')
    return lp,score,hessian
