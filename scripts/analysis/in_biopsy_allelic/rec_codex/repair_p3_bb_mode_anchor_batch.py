#!/usr/bin/env python3
"""Bounded batched anchored native BB prototype; no integration or accuracy claim.

One n/r and one requested count support are shared by a batch of eta values.
The scalar source's mode-relative ratios, shape derivatives and compensated
normalization moments are retained. Vectorization changes array operations,
requiring independent point checks before use in an integrator.
"""
import hashlib
import math
import os
from pathlib import Path

SCALAR_SHA = 'a31258aaac0fb39e91dda7aa634f0e810a6d0da144511b34d9c1b73a7028c266'
SCALAR = Path(__file__).with_name('repair_p3_bb_mode_anchor.py')
MAXIMUM_N = 100000
MAXIMUM_BATCH = 32
MAXIMUM_SUPPORT_ENTRIES = 100000
MAXIMUM_REQUESTED_ENTRIES = 100000
MAXIMUM_ESTIMATED_BYTES = 256*1024**2


def require(ok,message):
    if not ok:
        raise ValueError(message)


def memory_estimate(batch,n,requested,need_hessian=True):
    """Conservative simultaneous array capacity estimate, not measured RSS.

1024 bytes per full-support entry covers ratio/g/H, their mode-relative
counterparts, chain intermediates, masks, cumsum scratch and normalization
arrays; 256 bytes per requested entry covers all returned arrays and scratch.
Interpreter/library/caller allocations are excluded. The same conservative
bound is used when Hessians are omitted.
"""
    support = batch*(n+1)
    output = batch*requested
    return dict(batch=batch,n=n,requested_counts=requested,support_entries=support,
        requested_entries=output,estimated_array_bytes=1024*support+256*output+1048576,
        maximum_estimated_array_bytes=MAXIMUM_ESTIMATED_BYTES,
        includes_derivative_and_cumulative_sum_workspaces=True,
        excludes_interpreter_libraries_and_caller_memory=True,estimate_not_measured_peak=True)


def stable_bb_terms_batch(k,n,etas,linear_r,need_hessian=True):
    """Return B×K logPMF, B×K×2 score and optional B×K×2×2 Hessian.

Each node has its own anchor. Directional sums start there and never subtract
global prefixes. The arbitrary anchor offset cancels through full-support
normalization, including at changes of mode. Shape ratios stay in log space
when eta is large enough to underflow a beta shape.
"""
    require(os.environ.get('SLURM_JOB_ID','').isdigit() and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    require(hashlib.sha256(SCALAR.read_bytes()).hexdigest()==SCALAR_SHA,
            'Reviewed scalar anchored BB source differs')
    import numpy as np
    from scipy.special import expit,logsumexp
    require(isinstance(n,(int,float,np.integer,np.floating)) and np.isfinite(n) and
            int(n)==n and 0<=n<=MAXIMUM_N and np.ndim(linear_r)==0 and np.isfinite(linear_r),
            'Invalid bounded n or shared dispersion')
    require(isinstance(need_hessian,(bool,np.bool_)),'Hessian flag must be boolean')
    require(hasattr(etas,'__len__') and 1<=len(etas)<=MAXIMUM_BATCH,'Eta batch size outside declared scope')
    etas = np.asarray(etas,dtype=float)
    counts = np.atleast_1d(np.asarray(k))
    require(etas.ndim==1 and np.isfinite(etas).all() and counts.ndim==1 and np.isfinite(counts).all() and
            np.all(counts==np.floor(counts)) and np.all((counts>=0)&(counts<=n)),
            'Finite 1D eta batch and bounded integer counts required')
    require(len(etas)*len(counts)<=MAXIMUM_REQUESTED_ENTRIES,
            'Requested output entry count exceeds declared scope before integer-copy allocation')
    n,counts = int(n),counts.astype(int)
    batch = len(etas)
    estimate = memory_estimate(batch,n,len(counts),need_hessian)
    require(estimate['support_entries']<=MAXIMUM_SUPPORT_ENTRIES and
            estimate['requested_entries']<=MAXIMUM_REQUESTED_ENTRIES and
            estimate['estimated_array_bytes']<=MAXIMUM_ESTIMATED_BYTES,
            'Declared batch support/output/array-memory capacity exceeded')
    if n==0:
        return (np.zeros((batch,len(counts))),np.zeros((batch,len(counts),2)),
                np.zeros((batch,len(counts),2,2)) if need_hessian else None)
    m,q,vr,qr = expit(etas),expit(-etas),expit(linear_r),expit(-linear_r)
    rho = 1e-6+(1-2e-6)*vr
    rp = (1-2e-6)*vr*qr
    rpp = rp*(qr-vr)
    t = (1-rho)/rho
    tp = -rp/rho**2
    tpp = -rpp/rho**2+2*rp**2/rho**3
    u,w = tp/t,tpp/t-(tp/t)**2
    loga = np.log(t)-np.logaddexp(0.,-etas)
    logb = np.log(t)-np.logaddexp(0.,etas)
    j = np.arange(n,dtype=float)
    logj = np.full(n,-np.inf)
    logj[1:] = np.log(j[1:])
    logreverse = logj[::-1]
    ratio = (np.log(n-j)[None,:]-np.log(j+1)[None,:]+np.logaddexp(loga[:,None],logj[None,:])
             -np.logaddexp(logb[:,None],logreverse[None,:]))
    A,B = expit(loga[:,None]-logj[None,:]),expit(logb[:,None]-logreverse[None,:])
    sa,sb = np.column_stack((q,np.full(batch,u))),np.column_stack((-m,np.full(batch,u)))
    ratio_g = A[:,:,None]*sa[:,None,:]-B[:,:,None]*sb[:,None,:]
    ratio_h = None
    if need_hessian:
        shape_h = np.zeros((batch,2,2))
        shape_h[:,0,0],shape_h[:,1,1] = -m*q,w
        ratio_h = ((A-B)[:,:,None,None]*shape_h[:,None,:,:]
            +(A*expit(logj[None,:]-loga[:,None]))[:,:,None,None]*(sa[:,:,None]*sa[:,None,:])[:,None,:,:]
            -(B*expit(logreverse[None,:]-logb[:,None]))[:,:,None,None]*(sb[:,:,None]*sb[:,None,:])[:,None,:,:])
    anchors = np.zeros(batch,dtype=int)
    for node in range(batch):
        if loga[node]>0 and logb[node]>0:
            am1,bm1 = math.expm1(float(loga[node])),math.expm1(float(logb[node]))
            anchor = min(n,max(0,math.floor((n+1)*am1/(am1+bm1))))
            if anchor<n and ratio[node,anchor]>0:
                anchor += 1
            if anchor>0 and ratio[node,anchor-1]<0:
                anchor -= 1
            require((anchor==0 or ratio[node,anchor-1]>=-1e-12) and
                    (anchor==n or ratio[node,anchor]<=1e-12),'Interior node mode unresolved')
        elif loga[node]<=0 and logb[node]>=0:
            anchor = 0
        elif loga[node]>=0 and logb[node]<=0:
            anchor = n
        else:
            anchor = n if math.fsum(ratio[node].tolist())>0 else 0
        anchors[node] = anchor
    right = np.arange(n)[None,:]>=anchors[:,None]
    left = np.arange(n)[None,:]<anchors[:,None]

    def from_anchors(increments):
        trailing = (1,)*(increments.ndim-2)
        answer = np.zeros((batch,n+1)+increments.shape[2:])
        # Zeros outside each directional support do not involve cancellation
        # against another prefix. Addition order within a direction follows
        # the scalar helper's forward/right or reversed/left count order.
        right_values = np.where(right.reshape(right.shape+trailing),increments,0.)
        answer[:,1:] = np.cumsum(right_values,axis=1)
        del right_values
        left_values = np.where(left.reshape(left.shape+trailing),increments,0.)
        answer[:,:n] -= np.cumsum(left_values[:,::-1],axis=1)[:,::-1]
        return answer

    raw = from_anchors(ratio)
    raw_g = from_anchors(ratio_g)
    logpmf = raw-logsumexp(raw,axis=1)[:,None]
    probability = np.exp(logpmf)
    require(np.isfinite(raw).all() and np.isfinite(raw_g).all() and np.isfinite(logpmf).all(),
            'Nonfinite anchored batch mass/score')
    # C-level conversion preserves the float64 values and their summation order
    # while avoiding a Python conversion for every support entry.
    mean_g = np.array([[math.fsum((probability[node]*raw_g[node,:,d]).tolist())
                       for d in range(2)] for node in range(batch)])
    centered_g = raw_g-mean_g[:,None,:]
    hessian = None
    if need_hessian:
        raw_h = from_anchors(ratio_h)
        require(np.isfinite(raw_h).all(),'Nonfinite anchored batch Hessian')
        correction = np.empty((batch,2,2))
        for node in range(batch):
            for a in range(2):
                for b in range(2):
                    mean_h = math.fsum((probability[node]*raw_h[node,:,a,b]).tolist())
                    covariance = math.fsum((probability[node]*centered_g[node,:,a]*centered_g[node,:,b]).tolist())
                    correction[node,a,b] = mean_h+covariance
        hessian = raw_h[:,counts]-correction[:,None,:,:]
    score = centered_g[:,counts]
    require(np.isfinite(score).all() and (hessian is None or np.isfinite(hessian).all()),
            'Nonfinite normalized batch derivatives')
    return logpmf[:,counts],score,hessian
