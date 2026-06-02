#!/usr/bin/env python
"""
501_plain_nmf_converged.py — Plain sklearn NMF benchmark at matched configuration.

Runs 50 replicates of plain NMF at k in {10, 16} using the same HVG gene pool
and data matrix as cNMF, enabling a fair stability comparison.

Configuration (matched to cNMF factorize defaults):
    max_iter = 1000, tol = 1e-4, beta_loss = "frobenius"
    solver = "cd", init = "random", seeds 0..49

Metrics computed across replicates:
    - Cophenetic correlation: clustering agreement across replicate pairs
      (Brunet 2004; values near 1.0 indicate stable factorization)
    - Dispersion coefficient: spread of consensus clustering entries (Kim & Park 2007)
    - Silhouette width: within-cluster vs between-cluster distance on the consensus
    - ARI: adjusted Rand index across seed pairs (factor-assignment stability)
    - Mean top-100-gene Jaccard vs cNMF reference programs (program-content overlap)

Input : cnmf_runs/global/cnmf_tmp/global.norm_counts.h5ad  (895,542 × 2,500)
Gene list: cnmf_runs/global/global.overdispersed_genes.txt  (2,499; we intersect)
Output: benchmarks/factorization_showdown_plain_nmf.tsv     (long-format rows)
        benchmarks/plain_nmf_components_k{K}.npz            (stacked components)
"""
from __future__ import annotations

import os
import time
import json
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from sklearn.decomposition import NMF
from sklearn.metrics import adjusted_rand_score, silhouette_score
from scipy.cluster.hierarchy import linkage, cophenet, fcluster
from scipy.spatial.distance import squareform

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
NORM_H5AD = MCP / "cnmf_runs/global/cnmf_tmp/global.norm_counts.h5ad"
HVG_FILE = MCP / "cnmf_runs/global/global.overdispersed_genes.txt"
OUT_DIR = MCP / "benchmarks"
OUT_DIR.mkdir(exist_ok=True, parents=True)

# cNMF reference spectra (consensus.txt) — tab-separated: first col = program idx,
# remaining cols = gene names (the var_names used in cnmf run).
CNMF_CONSENSUS = {
    10: MCP / "cnmf_runs/global/global.spectra.k_10.dt_0_03.consensus.txt",
    16: MCP / "cnmf_runs/global/global.spectra.k_16.dt_0_03.consensus.txt",
}

N_REPLICATES = int(os.environ.get("PLAIN_N_REP", 50))
K_LIST = [int(x) for x in os.environ.get("PLAIN_K_LIST", "10,16").split(",")]
MAX_ITER = int(os.environ.get("PLAIN_MAX_ITER", 1000))
TOL = float(os.environ.get("PLAIN_TOL", 1e-4))
N_CELLS_SUB = int(os.environ.get("PLAIN_N_CELLS", 0))  # 0 = full
CONSENSUS_SUB = int(os.environ.get("PLAIN_CONSENSUS_N_CELLS", 30000))
TOP_N = 100


def load_input():
    print(f"[501] loading {NORM_H5AD}")
    a = ad.read_h5ad(NORM_H5AD)
    print(f"[501] shape: {a.shape}  dtype: {a.X.dtype}")
    genes = pd.Index(a.var_names)
    X = a.X
    if sp.issparse(X):
        X = X.toarray()
    X = np.ascontiguousarray(X, dtype=np.float64)  # NMF cd solver expects float64
    # Ensure non-negative (should already be, but guard)
    X = np.clip(X, 0.0, None)
    if N_CELLS_SUB and X.shape[0] > N_CELLS_SUB:
        rng = np.random.default_rng(0)
        idx = rng.choice(X.shape[0], N_CELLS_SUB, replace=False)
        X = X[idx]
        print(f"[501] subsampled to {X.shape[0]} cells (PLAIN_N_CELLS={N_CELLS_SUB})")
    return X, genes


def load_cnmf_reference(k):
    path = CNMF_CONSENSUS[k]
    df = pd.read_csv(path, sep="\t", index_col=0)
    # Rows = programs (1..k), cols = genes
    return df  # DataFrame: rows programs, cols genes


def top_genes(components, gene_names, n=TOP_N):
    """Return list[list[str]] of top-n genes per component (row-wise)."""
    out = []
    for row in components:
        order = np.argsort(row)[::-1][:n]
        out.append(set(gene_names[order]))
    return out


def jaccard(a: set, b: set):
    if not a and not b:
        return 0.0
    return len(a & b) / len(a | b)


def consensus_co_cluster(assignments_list):
    """Build n_samples x n_samples consensus co-clustering probability matrix.
    `assignments_list`: list of 1D arrays (length n_samples), one per replicate.
    Uses float32 to save memory.
    """
    n_samples = len(assignments_list[0])
    C = np.zeros((n_samples, n_samples), dtype=np.float32)
    for a in assignments_list:
        # Co-cluster matrix 1 if same cluster, else 0
        # Efficient construction via equality broadcast in blocks (to keep memory tame)
        a = a.astype(np.int32)
        # block size
        block = 2000
        for i0 in range(0, n_samples, block):
            i1 = min(i0 + block, n_samples)
            C[i0:i1] += (a[i0:i1, None] == a[None, :]).astype(np.float32)
    C /= len(assignments_list)
    return C


def cophenetic_from_consensus(C):
    D = 1.0 - C
    np.fill_diagonal(D, 0.0)
    D = np.clip(D, 0.0, 1.0)
    # condensed
    condensed = squareform(D, checks=False)
    Z = linkage(condensed, method="average")
    coph, _ = cophenet(Z, condensed)
    return float(coph), Z, D


def dispersion_coefficient(C):
    """Kim & Park 2007:  rho = (1/n^2) * sum_ij 4 * (C_ij - 0.5)^2 ."""
    n = C.shape[0]
    return float((4.0 * (C - 0.5) ** 2).sum() / (n * n))


def silhouette_from_consensus(C, Z, k):
    """Cluster consensus matrix into k clusters, compute silhouette width on 1-C distance."""
    D = 1.0 - C
    np.fill_diagonal(D, 0.0)
    D = np.clip(D, 0.0, 1.0)
    labels = fcluster(Z, t=k, criterion="maxclust")
    # silhouette on precomputed distance
    if len(np.unique(labels)) < 2:
        return float("nan")
    # silhouette_score supports precomputed distances
    return float(silhouette_score(D, labels, metric="precomputed"))


def mean_pairwise_ari(assignments_list):
    aris = []
    for i in range(len(assignments_list)):
        for j in range(i + 1, len(assignments_list)):
            aris.append(adjusted_rand_score(assignments_list[i], assignments_list[j]))
    return float(np.mean(aris)) if aris else float("nan")


def mean_program_jaccard_to_cnmf(components_list, gene_names, ref_df):
    """For each replicate's k programs, compute top-100 Jaccard against cNMF reference
    programs; take best match per rep program; then average over programs and reps."""
    ref_genes = ref_df.columns
    # top-100 per cNMF reference program
    ref_top = []
    for pidx in ref_df.index:
        vals = ref_df.loc[pidx].values
        order = np.argsort(vals)[::-1][:TOP_N]
        ref_top.append(set(ref_genes[order]))
    jaccs = []
    for comp in components_list:
        rep_top = top_genes(comp, gene_names, n=TOP_N)
        # For each rep program, best Jaccard to any ref program
        for rs in rep_top:
            best = max(jaccard(rs, r) for r in ref_top)
            jaccs.append(best)
    return float(np.mean(jaccs)) if jaccs else float("nan")


def run_one_k(X, gene_names, k, ref_df, log_prefix=""):
    t0 = time.time()
    n_samples = X.shape[0]
    components_list = []   # list of (k, n_genes) components
    assignments_consensus = []   # downsampled assignment vectors for consensus
    rng_sub = np.random.default_rng(12345)
    sub_n = min(CONSENSUS_SUB, n_samples)
    sub_idx = rng_sub.choice(n_samples, sub_n, replace=False)
    sub_idx.sort()
    reconstruction_errs = []
    n_iters = []
    for rep in range(N_REPLICATES):
        trep = time.time()
        m = NMF(
            n_components=k,
            init="random",
            random_state=rep,
            max_iter=MAX_ITER,
            tol=TOL,
            beta_loss="frobenius",
            solver="cd",
        )
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            W = m.fit_transform(X)  # (n_samples, k)
            H = m.components_        # (k, n_genes)
            conv_warn = any("did not converge" in str(w.message).lower() for w in caught)
        # Unit-normalize H rows for downstream program similarity
        H_norm = H / (np.linalg.norm(H, axis=1, keepdims=True) + 1e-10)
        components_list.append(H_norm)
        # Cell->program assignments (argmax W) on the downsampled subset
        assignments_consensus.append(W[sub_idx].argmax(axis=1))
        reconstruction_errs.append(float(m.reconstruction_err_))
        n_iters.append(int(m.n_iter_))
        print(f"{log_prefix}[k={k}][rep {rep}] iter={m.n_iter_} err={m.reconstruction_err_:.4f} "
              f"t={time.time()-trep:.1f}s conv_warn={conv_warn}")
    # Stack components for persistence
    stacked = np.concatenate(components_list, axis=0)  # (N_REP*k, n_genes)
    np.savez_compressed(
        OUT_DIR / f"plain_nmf_components_k{k}.npz",
        components=stacked, genes=gene_names.to_numpy(),
        n_rep=N_REPLICATES, k=k,
    )

    # Metrics on consensus matrix (downsampled)
    print(f"{log_prefix}[k={k}] building consensus co-clustering matrix ({sub_n}×{sub_n}) ...")
    C = consensus_co_cluster(assignments_consensus)
    coph, Z, _ = cophenetic_from_consensus(C)
    disp = dispersion_coefficient(C)
    sil = silhouette_from_consensus(C, Z, k)
    ari = mean_pairwise_ari(assignments_consensus)
    jac = mean_program_jaccard_to_cnmf(components_list, gene_names, ref_df) if ref_df is not None else float("nan")
    runtime_s = time.time() - t0
    # Long-format rows: one row per replicate for reconstruction_err/iters + one summary row for metrics
    rows = []
    for rep in range(N_REPLICATES):
        rows.append({
            "method": "plain_NMF_sklearn",
            "k": k,
            "replicate": rep,
            "metric": "reconstruction_err",
            "value": reconstruction_errs[rep],
            "note": f"n_iter={n_iters[rep]}",
        })
        rows.append({
            "method": "plain_NMF_sklearn",
            "k": k,
            "replicate": rep,
            "metric": "n_iter",
            "value": float(n_iters[rep]),
            "note": "max_iter=%d tol=%g" % (MAX_ITER, TOL),
        })
    summary_common = dict(
        method="plain_NMF_sklearn", k=k, replicate=-1,
        note=(f"n_rep={N_REPLICATES} consensus_n={sub_n} "
              f"max_iter={MAX_ITER} tol={TOL} solver=cd loss=frobenius init=random"),
    )
    rows += [
        dict(summary_common, metric="cophenetic", value=coph),
        dict(summary_common, metric="dispersion", value=disp),
        dict(summary_common, metric="silhouette", value=sil),
        dict(summary_common, metric="ari_mean", value=ari),
        dict(summary_common, metric="jaccard_top100_vs_cnmf", value=jac),
        dict(summary_common, metric="runtime_s", value=float(runtime_s)),
    ]
    return rows


def main():
    X, gene_names = load_input()
    # Align gene_names to cNMF reference genes (which use the same norm_counts var_names).
    all_rows = []
    for k in K_LIST:
        ref_path = CNMF_CONSENSUS.get(k)
        ref_df = None
        if ref_path is not None and ref_path.exists():
            ref_df = load_cnmf_reference(k)
            # Keep only genes present in both
            common = gene_names.intersection(ref_df.columns)
            ref_df = ref_df.reindex(columns=gene_names).fillna(0.0)
            print(f"[501] cNMF ref k={k} genes: {ref_df.shape}, common with norm_counts: {len(common)}")
        else:
            print(f"[501] WARN: no cNMF reference consensus file for k={k} — Jaccard will be NaN")
        rows = run_one_k(X, gene_names, k, ref_df, log_prefix="  ")
        all_rows.extend(rows)

    out_tsv = OUT_DIR / "factorization_showdown_plain_nmf.tsv"
    df = pd.DataFrame(all_rows)
    df.to_csv(out_tsv, sep="\t", index=False)
    print(f"[501] wrote {out_tsv}  ({len(df)} rows)")

    # Also write a compact summary JSON
    summary = {}
    for k in K_LIST:
        sub = df[(df["k"] == k) & (df["replicate"] == -1)].set_index("metric")["value"].to_dict()
        summary[f"k={k}"] = sub
    (OUT_DIR / "factorization_showdown_plain_nmf_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    print("[501] DONE.")


if __name__ == "__main__":
    main()
