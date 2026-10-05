#!/usr/bin/env python3
"""Isolated full-support beta-binomial ratio prototype; compute-only.

Same probability law and eta/r derivatives as the pinned stable kernel.
Direct left/right sums start at a mode, then normalize the complete support.
This source has no integration/fitter hook and supplies no acceptance claim.
"""
import argparse
import ast
import hashlib
import math
import os
from pathlib import Path

KERNEL_SHA = '461d325a91cd5ba728f91255311bb264f2dd392f5702a42a90231d905d0a7ac3'
KERNEL = Path('/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/MASLD_Model_Benchmark/executions/codex-rec-20260929T142434Z/repair_p3_launch_20261005T121740Z/repair_p3_integral.py')
MAXIMUM_N = 100000


def require(ok, message):
    if not ok:
        raise ValueError(message)


def stable_bb_terms(k, n, eta, linear_r, need_hessian=True):
    """Return logPMF, two-coordinate score, optional2x2 Hessian for requested k.

rho=1e-6+(1-2e-6)expit(r), exactly as in the pinned kernel. The anchor is
held fixed during raw derivative construction; normalization cancels its
arbitrary parameter-dependent offset, including at mode changes.
"""
    require(os.environ.get('SLURM_JOB_ID','').isdigit()
            and not os.environ.get('SLURM_ARRAY_TASK_ID'),
            'Individual compute allocation required before numerical imports')
    import numpy as np
    from scipy.special import expit, logsumexp
    k = np.atleast_1d(np.asarray(k))
    require(isinstance(n,(int,float,np.integer,np.floating)) and np.isfinite(n)
            and int(n)==n and 0<=n<=MAXIMUM_N and k.ndim==1 and np.isfinite(k).all()
            and np.all(k==np.floor(k)) and np.all((k>=0)&(k<=n))
            and np.isfinite(eta) and np.isfinite(linear_r),'Invalid bounded integer BB input')
    n, k = int(n), k.astype(int)
    if n==0:
        return (np.zeros(len(k)),np.zeros((len(k),2)),
                np.zeros((len(k),2,2)) if need_hessian else None)
    m,q,vr,qr = expit(eta),expit(-eta),expit(linear_r),expit(-linear_r)
    rho = 1e-6+(1-2e-6)*vr
    rp = (1-2e-6)*vr*qr
    rpp = rp*(qr-vr)
    t = (1-rho)/rho
    tp = -rp/rho**2
    tpp = -rpp/rho**2+2*rp**2/rho**3
    u,w = tp/t,tpp/t-(tp/t)**2
    loga = np.log(t)-np.logaddexp(0.,-eta)
    logb = np.log(t)-np.logaddexp(0.,eta)
    j = np.arange(n,dtype=float)
    logj = np.full(n,-np.inf)
    logj[1:] = np.log(j[1:])
    logreverse = logj[::-1]
    # No exponentiation of an underflowing shape is needed for these ratios.
    ratio = (np.log(n-j)-np.log(j+1)+np.logaddexp(loga,logj)
             -np.logaddexp(logb,logreverse))
    A,B = expit(loga-logj),expit(logb-logreverse)
    sa,sb = np.array([q,u]),np.array([-m,u])
    ratio_g = A[:,None]*sa-B[:,None]*sb
    ratio_h = None
    if need_hessian:
        shape_h = np.diag([-m*q,w])
        ratio_h = ((A-B)[:,None,None]*shape_h
            +(A*expit(logj-loga))[:,None,None]*np.outer(sa,sa)
            -(B*expit(logreverse-logb))[:,None,None]*np.outer(sb,sb))
    if loga>0 and logb>0:
        am1,bm1 = math.expm1(float(loga)),math.expm1(float(logb))
        anchor = min(n,max(0,math.floor((n+1)*am1/(am1+bm1))))
        # A roundoff-level tied mode may choose either neighboring maximizer.
        if anchor<n and ratio[anchor]>0: anchor += 1
        if anchor>0 and ratio[anchor-1]<0: anchor -= 1
        require((anchor==0 or ratio[anchor-1]>=-1e-12)
                and (anchor==n or ratio[anchor]<=1e-12),'Interior mode unresolved')
    elif loga<=0 and logb>=0:
        anchor = 0
    elif loga>=0 and logb<=0:
        anchor = n
    else:
        # Both shapes <=1: two endpoint maxima. fsum is used only to choose
        # the higher endpoint, never to subtract global prefixes for PMFs.
        anchor = n if math.fsum(map(float,ratio))>0 else 0

    def from_anchor(increments):
        answer = np.zeros((n+1,)+increments.shape[1:],dtype=float)
        if anchor<n: answer[anchor+1:] = np.cumsum(increments[anchor:],axis=0)
        if anchor>0: answer[:anchor] = -np.cumsum(increments[:anchor][::-1],axis=0)[::-1]
        return answer

    raw = from_anchor(ratio)
    raw_g = from_anchor(ratio_g)
    logpmf = raw-logsumexp(raw)
    probability = np.exp(logpmf)
    require(np.isfinite(raw).all() and np.isfinite(raw_g).all()
            and np.isfinite(logpmf).all(),'Nonfinite anchored mass/score')
    # fsum reduces cancellation in the normalization moments. Computing all
    # support entries is necessary even when only a subset k was requested.
    mean_g = np.array([math.fsum(map(float,probability*raw_g[:,d])) for d in range(2)])
    centered_g = raw_g-mean_g
    hessian = None
    if need_hessian:
        raw_h = from_anchor(ratio_h)
        require(np.isfinite(raw_h).all(),'Nonfinite anchored Hessian')
        correction = np.empty((2,2))
        for a in range(2):
            for b in range(2):
                mean_h = math.fsum(map(float,probability*raw_h[:,a,b]))
                covariance = math.fsum(map(float,probability*centered_g[:,a]*centered_g[:,b]))
                correction[a,b] = mean_h+covariance
        hessian = raw_h[k]-correction
    score = centered_g[k]
    require(np.isfinite(score).all() and (hessian is None or np.isfinite(hessian).all()),
            'Nonfinite normalized derivatives')
    return logpmf[k],score,hessian


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inspect',action='store_true',required=True)
    args = parser.parse_args()
    require(hashlib.sha256(KERNEL.read_bytes()).hexdigest()==KERNEL_SHA,'Pinned kernel source differs')
    ast.parse(KERNEL.read_text())
    print('Source-only prototype; pinned kernel checked; maximum n=100000; numerical checks not run.')


if __name__=='__main__': main()
