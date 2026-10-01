#!/usr/bin/env python3
"""Shared numerical definitions for the bounded 99-donor comparison."""
import hashlib
import numpy as np
import pandas as pd

SEED=20260915
RIDGE_C=np.logspace(-5,1,13)
PCR_K=(4,8,16,32,64)
PLS_K=(2,4,8,10,16)


def stable(value):return int(hashlib.sha256(str(value).encode()).hexdigest()[:15],16)


def region_splits(reg):
    """Keep every tile and any overlapping 16-kb input windows together."""
    tiles=sorted(reg.tile_index.unique());parent={t:t for t in tiles}
    def root(t):
        while parent[t] != t:t=parent[t]
        return t
    intervals=[]
    for row in reg.itertuples():
        mid=(int(row.start0)+int(row.end0))//2
        intervals.append((row.chrom,mid-8192,mid+8192,row.tile_index))
    for i,a in enumerate(intervals):
        for b in intervals[i+1:]:
            if a[0] == b[0] and a[1] < b[2] and b[1] < a[2]:parent[root(b[3])]=root(a[3])
    components=np.array([root(t) for t in reg.tile_index])
    outer=reg.region_role.to_numpy()
    for block in set(components):
        if len(set(outer[components == block])) != 1:raise ValueError("Actual input windows span outer region split")
    training_components=sorted(set(components[outer == "train"]),key=lambda c:stable("donor-inner-region|"+str(c)))
    if len(training_components)<4:raise ValueError("Too few independent inner region components")
    inner_valid=set(training_components[:len(training_components)//2])
    inner=np.array(["outer_held" if o == "held" else "inner_valid" if c in inner_valid else "inner_train" for c,o in zip(components,outer)])
    return components,inner


def kernel_projection(x,donors,keep):
    """Full-feature standardized RNA kernel and its exact training sample space."""
    donors=np.asarray(donors);raw=np.asarray(x[:,keep],float)
    mean=raw[donors].mean(0);sd=raw[donors].std(0)
    admitted=sd>1e-10
    if admitted.sum()<2:raise ValueError("Too few variable input genes")
    z=(raw[:,admitted]-mean[admitted])/sd[admitted]
    kernel=z@z[donors].T/z.shape[1]
    values,vectors=np.linalg.eigh(kernel[donors]);order=np.argsort(values)[::-1]
    values=values[order];vectors=vectors[:,order]
    valid=values>max(values[0]*1e-10,1e-12);values=values[valid];vectors=vectors[:,valid]
    scores=kernel@vectors/np.sqrt(values)[None]
    nn=np.zeros((len(x),32),dtype=np.float32);rank=min(32,scores.shape[1])
    nn[:,:rank]=(scores[:,:rank]/np.maximum(scores[donors,:rank].std(0),1e-10)).astype(np.float32)
    return {"kernel":kernel,"values":values,"vectors":vectors,"scores":scores,"nn":nn,
            "genes":np.flatnonzero(keep)[admitted],"training_donors":donors}


def sequence_projection(raw,regions):
    regions=np.asarray(regions);raw=np.asarray(raw,float)
    mean=raw[regions].mean(0);sd=raw[regions].std(0);keep=sd>1e-10
    z=(raw[:,keep]-mean[keep])/sd[keep]
    kernel=z@z[regions].T/max(z.shape[1],1)
    values,vectors=np.linalg.eigh(kernel[regions]);order=np.argsort(values)[::-1]
    values=values[order];vectors=vectors[:,order];valid=values>max(values[0]*1e-10,1e-12)
    values=values[valid];vectors=vectors[:,valid]
    scores=kernel@vectors/np.sqrt(values)[None]
    nn=np.zeros((len(raw),32),np.float32);rank=min(32,scores.shape[1])
    nn[:,:rank]=(scores[:,:rank]/np.maximum(scores[regions,:rank].std(0),1e-10)).astype(np.float32)
    return {"kernel":kernel,"nn":nn,"training_regions":regions}


def ridge_kernel(kernel,train,y,c):
    mean=y.mean(0);k=kernel[train]
    penalty=float(c)*max(float(np.linalg.eigvalsh(k).max()),1e-12)
    return kernel@np.linalg.solve(k+penalty*np.eye(len(train)),y-mean)+mean


def rna_candidates(projection,y):
    """Separate native-scale predictions; PLS loadings fitted to training outputs."""
    tr=projection["training_donors"];scores=projection["scores"]
    out={};k=projection["kernel"]
    def path(kernel,prefix,grid):
        eigen,vectors=np.linalg.eigh(kernel[tr]);eigen=np.maximum(eigen,0)
        cross=kernel@vectors;mean=y.mean(0);target=vectors.T@(y-mean)
        for c in grid:
            out[f"{prefix}:c{c:g}"]=cross@(target/(eigen[:,None]+float(c)*max(eigen.max(),1e-12)))+mean
    path(k,"full_rna_ridge",RIDGE_C)
    for requested in PCR_K:
        rank=min(requested,scores.shape[1]);p=scores[:,:rank];kernel=p@p[tr].T
        path(kernel,f"pcr:k{requested}",RIDGE_C)
    yc=(y-y.mean(0))/np.maximum(y.std(0),1e-8)
    directions,_,_=np.linalg.svd(scores[tr].T@yc,full_matrices=False)
    for requested in PLS_K:
        rank=min(requested,directions.shape[1]);p=scores@directions[:,:rank];kernel=p@p[tr].T
        path(kernel,f"pls_svd_ridge:k{requested}",(0.001,0.01,0.1,1.))
    return out


def sequence_mean(projection,means,c):
    return ridge_kernel(projection["kernel"],projection["training_regions"],means[:,None],c)[:,0]


def product_prediction(rna,sequence,y,c,baseline,additive=False):
    """Full-RNA linear effect with global or sequence-dependent coefficients."""
    d=rna["training_donors"];r=sequence["training_regions"]
    residual=y-y.mean(0)
    if additive:
        donor=ridge_kernel(rna["kernel"],d,residual.mean(1,keepdims=True),c)[:,0]
        return baseline[None]+donor[:,None]
    kx=rna["kernel"];ks=sequence["kernel"]+1.
    lx,ux=np.linalg.eigh(kx[d]);ls,us=np.linalg.eigh(ks[r])
    lx=np.maximum(lx,0);ls=np.maximum(ls,0)
    penalty=float(c)*max(float(lx.max()*ls.max()),1e-12)
    coefficients=ux@((ux.T@residual@us)/(lx[:,None]*ls[None,:]+penalty))@us.T
    return baseline[None]+kx@coefficients@ks.T


def groups(reg,levels):
    out=[("all_rna","all",np.arange(len(reg))),("annotated_all_rna","all",np.arange(len(reg)))]
    for level in levels:
        if level == "chromosome":
            for chrom in sorted(reg.chrom.unique()):out.extend([(level,chrom,np.flatnonzero(reg.chrom == chrom)),("random_"+level,chrom,np.flatnonzero(reg.chrom == chrom))])
        elif level in ("overlap","within_1mb"):
            for j,row in enumerate(reg.itertuples()):out.extend([(level,row.region_key,np.array([j])),("random_"+level,row.region_key,np.array([j]))])
        else:raise ValueError(level)
    return out


def gene_mask(x,donors,gene_chrom,gene_start,gene_end,reg,arm,group,targets):
    if arm == "all_rna":return np.ones(x.shape[1],bool),np.zeros(x.shape[1],bool)
    variable=x[donors].var(0)>1e-10;known=gene_chrom != "unknown"
    if arm == "annotated_all_rna":return known,np.zeros(x.shape[1],bool)
    level=arm.removeprefix("random_")
    if level == "chromosome":remove=(gene_chrom == group)&variable
    else:
        row=reg.iloc[int(targets[0])];distance=1_000_000 if level == "within_1mb" else 0
        remove=(gene_chrom == row.chrom)&(gene_start<row.end0+distance)&(gene_end>row.start0-distance)&variable
    if arm.startswith("random_"):
        rng=np.random.default_rng(SEED+stable(arm+"|"+group));count=int(remove.sum())
        remove=np.zeros(x.shape[1],bool);remove[rng.choice(np.flatnonzero(known&variable),count,replace=False)]=True
    return known&~remove,remove
