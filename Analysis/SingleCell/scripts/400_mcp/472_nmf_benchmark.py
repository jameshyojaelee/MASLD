#!/usr/bin/env python
"""
472_nmf_benchmark.py — Plain NMF (sklearn) vs cNMF consensus at matched k.

Produces cophenetic correlation, program-Jaccard, and runtime comparison for one name/k.
"""
from __future__ import annotations

import argparse
import os
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import NMF
from scipy.cluster.hierarchy import linkage, cophenet

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
INPUTS = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/inputs"
CNMF_ROOT = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/cnmf_runs"
BENCH = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/mcp/benchmarks"
BENCH.mkdir(parents=True, exist_ok=True)


def main(args):
    t0 = time.time()
    f = INPUTS / f"atlas_cnmf_{args.name[4:] if args.name.startswith('pct_') else 'global'}.h5ad"
    a = ad.read_h5ad(f)
    X = a.layers["counts"]
    # Top-HVG filter via variance
    if sp.issparse(X):
        var = np.asarray(X.power(2).mean(axis=0) - np.asarray(X.mean(axis=0)) ** 2).ravel()
    else:
        var = X.var(axis=0)
    order = np.argsort(var)[::-1][: args.n_hvg]
    Xf = X[:, order]
    if sp.issparse(Xf):
        Xf = Xf.toarray()
    # Normalize by total counts
    row_sum = Xf.sum(axis=1, keepdims=True).clip(min=1)
    Xn = Xf / row_sum * 1e4
    Xn = np.log1p(Xn)
    Xn = np.clip(Xn, 0, None)
    print(f"[472] Xn shape {Xn.shape}")

    # n_replicates plain NMFs
    basis_stack = []
    for rep in range(args.n_rep):
        m = NMF(n_components=args.k, init="random", random_state=rep, max_iter=500, tol=1e-4)
        H = m.fit_transform(Xn)
        W = m.components_
        basis_stack.append(W)
    t1 = time.time()
    # Cophenetic correlation via column-cosine of stacked basis
    B = np.concatenate(basis_stack, axis=0)  # (n_rep*k) x genes
    D = 1 - (B @ B.T) / (np.linalg.norm(B, axis=1)[:, None] * np.linalg.norm(B, axis=1)[None, :] + 1e-8)
    Z = linkage(D[np.triu_indices(D.shape[0], k=1)], method="average")
    c, _ = cophenet(Z, D[np.triu_indices(D.shape[0], k=1)])

    # cNMF comparison (if available)
    dt_tag = "0_03"
    full_f = CNMF_ROOT / args.name / args.name / f"{args.name}.gene_spectra_score.k_{args.k}.dt_{dt_tag}.txt"
    cnmf_c = np.nan
    jac_vals = []
    if full_f.exists():
        full = pd.read_csv(full_f, sep="\t", index_col=0)
        # Plain NMF program top genes (ensembl_base via adata var)
        gene_names = a.var_names[order]
        rep_tops = []
        for W in basis_stack:
            # W: k x genes_hvg
            tops = [set(gene_names[np.argsort(-W[c, :])[:100]]) for c in range(W.shape[0])]
            rep_tops.append(tops)
        # Jaccard of cNMF top-100 against best plain-NMF program
        for p in full.index:
            sf = set(full.loc[p].nlargest(100).index)
            best = 0.0
            for reptop in rep_tops:
                for s in reptop:
                    j = len(sf & s) / max(1, len(sf | s))
                    if j > best:
                        best = j
            jac_vals.append({"program": p, "best_jaccard_vs_plain_nmf_rep": best})
    pd.DataFrame(jac_vals).to_csv(BENCH / f"plain_nmf_vs_cnmf_{args.name}_k{args.k}.tsv", sep="\t", index=False)
    pd.DataFrame([{
        "name": args.name,
        "k": args.k,
        "n_rep": args.n_rep,
        "n_hvg": args.n_hvg,
        "plain_nmf_cophenetic": float(c),
        "plain_nmf_time_sec": t1 - t0,
    }]).to_csv(BENCH / f"plain_nmf_stability_{args.name}_k{args.k}.tsv", sep="\t", index=False)
    print(f"[472] plain NMF cophenetic={c:.3f}; wrote benchmarks.")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--name", required=True)
    ap.add_argument("--k", type=int, required=True)
    ap.add_argument("--n-rep", type=int, default=20)
    ap.add_argument("--n-hvg", type=int, default=2500)
    main(ap.parse_args())
