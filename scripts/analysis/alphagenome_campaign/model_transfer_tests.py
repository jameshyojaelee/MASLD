#!/usr/bin/env python3
"""Small scientific invariants for the separate-unit shared/private objective."""
import json
import numpy as np
from model_transfer_fit import delta,variant_weights,fit_prediction,compress


def main():
    rng=np.random.default_rng(1103)
    alleles=rng.normal(size=(5,4,7))
    np.testing.assert_allclose(delta(alleles[:,[1,0,3,2]]),-delta(alleles),atol=0)
    identical=alleles.copy();identical[:,1]=identical[:,0];identical[:,3]=identical[:,2]
    np.testing.assert_array_equal(delta(identical),np.zeros((5,7)))
    n=40;v=rng.normal(size=(n,3));t=rng.normal(size=(n,2));c=rng.normal(size=(n,2))
    pi=np.arange(n);train=np.arange(n)<30;test=~train;wp=np.ones(n)
    yp=v@np.array([.3,-.5,.8])+t@np.array([.1,-.2])
    yr=v@np.array([.3,-.5,.8])+c@np.array([-.1,.2])
    p,r,_,_=fit_prediction(v,t,c,pi,yp,yr,wp,train,train,test,test,.01,1)
    # Independent assay rescaling must rescale that assay's native predictions
    # without changing the other output or normalized joint optimization.
    ps,rs,_,_=fit_prediction(v,t,c,pi,7*yp,.25*yr,wp,train,train,test,test,.01,1)
    np.testing.assert_allclose(ps,7*p,rtol=1e-10,atol=1e-10)
    np.testing.assert_allclose(rs,.25*r,rtol=1e-10,atol=1e-10)
    # Duplicate one target measurement with reciprocal weights: assay loss and
    # coefficients must remain unchanged, rather than counting it twice.
    dup=np.r_[np.arange(n),0];td=t[dup];pid=pi[dup];yd=yp[dup]
    trd=train[dup];ted=np.r_[test,False];weights=variant_weights(pid)
    dp,dr,_,_=fit_prediction(v,td,c,pid,yd,yr,weights,trd,train,ted,test,.01,1)
    np.testing.assert_allclose(dp,p,rtol=1e-10,atol=1e-10)
    np.testing.assert_allclose(dr,r,rtol=1e-10,atol=1e-10)
    # Held examples cannot change the fitted representation of training rows.
    x=rng.normal(size=(n,12));a=compress(x,train,4)
    x[test]*=1000;b=compress(x,train,4)
    np.testing.assert_array_equal(a[train],b[train])
    print(json.dumps({"identical_and_swapped_alleles":"passed","independent_assay_units":"passed",
        "variant_balanced_duplicate_target":"passed","training_only_projection":"passed","seed":1103},indent=2))


if __name__ == "__main__":main()
