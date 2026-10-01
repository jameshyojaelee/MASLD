#!/usr/bin/env python3
"""Twelve cached frozen comparisons using installed JAX and SciPy only.

Ridge minimizes ||Xw-y||² + alpha*||w||² via scipy.sparse.linalg.lsqr
with damp=sqrt(alpha), matching the zero-intercept ridge recipe. This avoids
requiring scikit-learn in the isolated local AlphaGenome environment.
"""
import argparse
import json
from pathlib import Path
import time

import numpy as np
import pandas as pd
from scipy.sparse.linalg import lsqr

from model_frozen import C, POOLINGS, fit_shared


def fit(args):
    out=args.out/"comparisons"
    out.mkdir(parents=True,exist_ok=False)
    labels=pd.read_csv(args.labels,sep="\t").set_index("lead_variant_id")
    manifests={l:pd.read_csv(args.out/str(l)/"manifest.tsv",sep="\t") for l in (2048,16384)}
    if not np.array_equal(manifests[2048].key,manifests[16384].key):
        raise ValueError("Length manifests differ")
    matched=np.ones(len(manifests[2048]),dtype=bool)
    coverage=[]
    for length in (2048,16384):
        meta=json.loads((args.out/str(length)/"complete.json").read_text())
        if not meta["complete"] or meta["card_tag"] != "l40s":
            raise ValueError("Extraction incomplete or hardware differs")
        with np.load(args.out/str(length)/"coverage.npz",allow_pickle=False) as archive:
            for pool in POOLINGS:
                matched &= archive[pool]
                coverage.append({"length":length,"pooling":pool,"eligible_rows":int(archive[pool].sum()),"total_rows":len(matched)})
    manifest=manifests[2048].loc[matched].copy()
    y=labels.loc["chr"+manifest.key,"beta_alt"].to_numpy(dtype=np.float32)
    folds=manifest.heldout_fold.to_numpy()
    if not len(y) or len(set(folds)) != 5:
        raise ValueError("Matched rows cannot support all five folds")
    pd.DataFrame(coverage).to_csv(out/"length_specific_coverage.tsv",sep="\t",index=False)
    recipes=[json.loads(p.read_text()) for p in sorted(args.configs.glob("frozen_*.json"))]
    if len(recipes) != 12:
        raise ValueError("Expected twelve frozen recipes")
    predictions=manifest[["key","block_1mb","heldout_fold"]].copy()
    predictions["beta_alt"]=y
    metrics,failures=[],[]
    for recipe in recipes:
        tick=time.monotonic()
        array=np.load(args.out/str(recipe["length"])/(recipe["pooling"]+".npy"),mmap_mode="r")
        features=np.asarray(array[matched],dtype=np.float32)
        pred=np.full(len(y),np.nan,dtype=np.float32)
        for fold in range(5):
            train,test=folds != fold,folds == fold
            try:
                if recipe["head"] == "ridge":
                    delta=features[:,1]-features[:,0]
                    scale=delta[train].std(axis=0);scale[scale < 1e-6]=1
                    result=lsqr(delta[train]/scale,y[train],damp=np.sqrt(recipe["alpha"]),atol=1e-6,btol=1e-6,iter_lim=10000)
                    if result[1] not in (0,1,2,4,5):
                        raise ValueError(f"Ridge LSQR did not converge: istop={result[1]}")
                    pred[test]=(delta[test]/scale) @ result[0]
                elif recipe["head"] == "shared_mlp64":
                    pred[test]=fit_shared(features[train,0],features[train,1],y[train],features[test,0],features[test,1],recipe)
                else:
                    raise ValueError("Unsupported head")
            except Exception as exc:
                failures.append({"recipe":recipe["id"],"fold":fold,"error":repr(exc)})
        predictions[recipe["id"]]=pred
        okay=np.isfinite(pred).all()
        metrics.append({"recipe":recipe["id"],"matched_rows":len(y),"blocks":manifest.block_1mb.nunique(),
            "macro_spearman":C.macro_over_folds(y,pred,folds) if okay else np.nan,
            "RMSE_beta_units":float(np.sqrt(np.mean((pred-y)**2))) if okay else np.nan,
            "seconds":time.monotonic()-tick,"failed_folds":sum(r["recipe"] == recipe["id"] for r in failures),
            "evidence":"fixed_recipe_single_seed_development_cross_validation"})
        predictions.to_csv(out/"predictions.tsv.gz",sep="\t",index=False)
        pd.DataFrame(metrics).to_csv(out/"performance.tsv",sep="\t",index=False)
        pd.DataFrame(failures,columns=["recipe","fold","error"]).to_csv(out/"failed_fits.tsv",sep="\t",index=False)
        print(json.dumps(metrics[-1]),flush=True)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out",type=Path,required=True)
    parser.add_argument("--labels",type=Path,default=C.C2_LABELS)
    parser.add_argument("--configs",type=Path,required=True)
    fit(parser.parse_args())
