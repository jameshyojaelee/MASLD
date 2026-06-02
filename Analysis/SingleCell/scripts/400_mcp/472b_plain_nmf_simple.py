#!/usr/bin/env python
"""
472b_plain_nmf_simple.py — Minimal plain-NMF benchmark at k=16.
Fit sklearn NMF 20 times with different seeds, compute cophenetic correlation
of stacked program components. Compare to cNMF's cophenetic.
"""
from __future__ import annotations
import os, time
from pathlib import Path
import anndata as ad, numpy as np, pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import NMF
from scipy.cluster.hierarchy import linkage, cophenet
from scipy.spatial.distance import squareform

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
K = 16
N_REP = 20
N_HVG = 2500

# Use the cnmf normalized counts for direct comparison
a = ad.read_h5ad(MCP / "cnmf_runs/global/cnmf_tmp/global.norm_counts.h5ad")
print(f"[472b] norm_counts shape: {a.shape}")
X = a.X
if sp.issparse(X):
    X = X.toarray()
# Variance-sort for HVG
v = X.var(axis=0)
order = np.argsort(v)[::-1][:N_HVG]
Xf = X[:, order]
print(f"[472b] HVG-filtered: {Xf.shape}")

# Run plain NMFs
t0 = time.time()
basis_stack = []
for rep in range(N_REP):
    m = NMF(n_components=K, init="random", random_state=rep, max_iter=300, tol=1e-4)
    _ = m.fit_transform(Xf)
    W = m.components_  # K x HVG
    # Unit-normalize
    W = W / (np.linalg.norm(W, axis=1, keepdims=True) + 1e-10)
    basis_stack.append(W)
    print(f"[472b] rep {rep}: done at {time.time()-t0:.0f}s")
B = np.concatenate(basis_stack, axis=0)  # (N_REP*K) x HVG

# Cophenetic via cosine distances
D = 1 - B @ B.T
np.fill_diagonal(D, 0.0)
D = np.clip(D, 0, 2)
try:
    Z = linkage(squareform(D, checks=False), method="average")
    coph, _ = cophenet(Z, squareform(D, checks=False))
except Exception as e:
    coph = np.nan

out_dir = MCP / "benchmarks"
out_dir.mkdir(exist_ok=True)
pd.DataFrame([{
    "method": "plain_NMF_sklearn", "k": K, "n_rep": N_REP, "n_hvg": N_HVG,
    "cophenetic": float(coph), "runtime_s": float(time.time() - t0),
}]).to_csv(out_dir / "plain_nmf_vs_cnmf.tsv", sep="\t", index=False)
print(f"[472b] plain NMF cophenetic: {coph:.3f} (cNMF was 0.951 at k=16)")
