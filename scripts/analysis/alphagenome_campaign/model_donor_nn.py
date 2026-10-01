#!/usr/bin/env python3
"""Small FiLM/additive context heads; sequence features remain fixed."""
import time
import jax,jax.numpy as jnp
import numpy as np
import optax
from model_donor_common import SEED


def initialize(seed=SEED):
    keys=jax.random.split(jax.random.PRNGKey(seed),6)
    return {"sequence":jax.random.normal(keys[0],(32,32))*np.sqrt(1/32),
        "context":jax.random.normal(keys[1],(32,32))*np.sqrt(1/32),
        "gamma":jax.random.normal(keys[2],(32,32))*np.sqrt(1/32),
        "beta":jax.random.normal(keys[3],(32,32))*np.sqrt(1/32),
        "output":jax.random.normal(keys[4],(32,))*np.sqrt(1/32),"bias":jnp.zeros(32)}


def forward(params,rna,sequence,interaction):
    s=jax.nn.gelu(sequence@params["sequence"]+params["bias"])
    c=jax.nn.gelu(rna@params["context"])
    if interaction:
        h=jax.nn.gelu(s[None]*(1+(c@params["gamma"])[:,None])+ (c@params["beta"])[:,None])
        return jnp.einsum("drk,k->dr",h,params["output"])
    # Exact additive control: separate nonlinear RNA and sequence terms.
    return (s@params["output"])[None]+((c@params["beta"])@params["output"])[:,None]


def balanced_partners(donors,epoch):
    rng=np.random.default_rng(SEED+epoch);ring=rng.permutation(donors)
    shifts=rng.choice(np.arange(1,len(ring)),2,replace=False)
    a=np.tile(ring,2);b=np.concatenate([np.roll(ring,int(k)) for k in shifts])
    aa=np.zeros(198,np.int32);bb=np.zeros(198,np.int32);mask=np.zeros(198,np.float32)
    aa[:len(a)]=a;bb[:len(b)]=b;mask[:len(a)]=1
    return aa,bb,mask


OPTIMIZER=optax.chain(optax.clip_by_global_norm(1.),optax.adamw(1e-3,weight_decay=1e-3))


def loss(params,rna,seq,target,mask,donor_mask,scale,level,a,b,pair_mask,objective,interaction):
    q=forward(params,rna,seq,interaction)
    if objective == "profile":prediction=level+scale*q;residual=(prediction-target)*mask
    elif objective == "residual":
        centered=q-jnp.sum(q*donor_mask[:,None],axis=0,keepdims=True)/jnp.sum(donor_mask)
        residual=(scale*centered-target)*mask
    elif objective == "pairwise":
        residual=(scale*(q[a]-q[b])-(target[a]-target[b]))*pair_mask
        return jnp.sum(residual**2)/jnp.sum(pair_mask)
    else:raise ValueError(objective)
    return jnp.sum(residual**2)/jnp.sum(mask)


def _step(params,state,rna,seq,target,mask,donor_mask,scale,level,a,b,pair_mask,*,objective,interaction):
    value,gradient=jax.value_and_grad(loss)(params,rna,seq,target,mask,donor_mask,scale,level,a,b,pair_mask,objective,interaction)
    updates,state=OPTIMIZER.update(gradient,state,params)
    return optax.apply_updates(params,updates),state,value


step=jax.jit(_step,static_argnames=("objective","interaction"))
predict=jax.jit(forward,static_argnames=("interaction",))


def train(rna,sequence,y,donors,regions,baseline,objective,interaction,checkpoints,progress=None):
    """Only explicitly indexed training targets enter optimizer or scales."""
    params=initialize();state=OPTIMIZER.init(params);n,nr=y.shape
    mask=np.zeros_like(y,np.float32);mask[np.ix_(donors,regions)]=1
    donor_mask=np.zeros(n,np.float32);donor_mask[donors]=1
    observed=y[np.ix_(donors,regions)];mu=observed.mean(0)
    target=np.zeros_like(y,np.float32)
    if objective == "profile":
        target[np.ix_(donors,regions)]=observed;level=float(observed.mean());scale=max(float(observed.std()),1e-6)
    else:
        target[np.ix_(donors,regions)]=observed-mu;level=0.;scale=max(float((observed-mu).std()),1e-6)
    arrays=[jnp.asarray(v) for v in (rna,sequence,target,mask,donor_mask,np.float32(scale),np.float32(level))]
    outputs={};start=time.monotonic();timings=[]
    for epoch in range(1,max(checkpoints)+1):
        tick=time.monotonic();a,b,valid=balanced_partners(donors,epoch)
        pair_mask=valid[:,None]*mask[donors[0]][None]
        params,state,value=step(params,state,*arrays,jnp.asarray(a),jnp.asarray(b),jnp.asarray(pair_mask),objective=objective,interaction=interaction)
        number=float(value)
        if not np.isfinite(number):raise ValueError("Nonfinite donor objective")
        timings.append(time.monotonic()-tick)
        if epoch == 100 and progress is not None:
            first100=time.monotonic()-start;inference_tick=time.monotonic()
            np.asarray(predict(params,jnp.asarray(rna),jnp.asarray(sequence),interaction=interaction))
            progress({"first100_seconds":first100,"warm_iteration_seconds":float(np.median(timings[10:])),
                "warm_iteration_p90_seconds":float(np.quantile(timings[10:],.9)),
                "full_panel_inference_seconds":time.monotonic()-inference_tick,"steps":100})
        if epoch in checkpoints:
            q=np.asarray(predict(params,jnp.asarray(rna),jnp.asarray(sequence),interaction=interaction),dtype=float)
            outputs[epoch]=level+scale*q if objective == "profile" else baseline[None]+scale*(q-q[donors].mean(0))
    return outputs,{"seconds":time.monotonic()-start,"steps":max(checkpoints),"objective":objective,"interaction":interaction,
        "training_donors":len(donors),"training_regions":len(regions),"normalization":"training_targets_only",
        "pair_partners_per_donor_per_epoch":2,"global_level":level,"target_scale":scale}
