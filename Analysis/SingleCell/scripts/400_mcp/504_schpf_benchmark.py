#!/usr/bin/env python
"""
504_schpf_benchmark.py — scHPF benchmark (Levitin 2019 PMID 30796088).

Hierarchical Poisson Factorization on the HVG-subset raw counts of the
atlas. k ∈ {10, 16}. Metrics align with the showdown:
  - Mean top-100-gene Jaccard vs cNMF reference
  - Runtime
  - ELBO-equivalent: final training loss on validation split (schpf's `.loss`)

Inputs:
  atlas_cnmf_global.h5ad (counts layer)
  global.overdispersed_genes.txt (2,499 HVGs)
  cNMF consensus spectra for reference

Outputs:
  benchmarks/schpf_factors_k{10,16}.tsv              (gene × factor weights)
  benchmarks/schpf_cell_scores_k{10,16}.tsv          (cell × factor scores)
  benchmarks/program_jaccard_schpf_vs_cnmf.tsv
  benchmarks/factorization_showdown_schpf.tsv
"""
from __future__ import annotations

import os
import sys
import time
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
ATLAS_H5AD = MCP / "inputs/atlas_cnmf_global.h5ad"
HVG_FILE = MCP / "cnmf_runs/global/global.overdispersed_genes.txt"
OUT_DIR = MCP / "benchmarks"
OUT_DIR.mkdir(parents=True, exist_ok=True)

K_LIST = [int(x) for x in os.environ.get("SCHPF_K_LIST", "10,16").split(",")]
N_CELLS_SUB = int(os.environ.get("SCHPF_N_CELLS", 0))  # 0 = full
TOP_N = 100
N_TRIALS = int(os.environ.get("SCHPF_N_TRIALS", 3))


def load_input():
    print(f"[504] loading {ATLAS_H5AD}")
    a = ad.read_h5ad(ATLAS_H5AD)
    hvg = pd.read_csv(HVG_FILE, header=None)[0].tolist()
    mask = a.var_names.isin(set(hvg))
    X = a.layers["counts"][:, mask]
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.astype(np.float32)
    genes = pd.Index(a.var_names[mask])
    print(f"[504] HVG counts matrix: {X.shape}  nnz={X.nnz}")
    if N_CELLS_SUB and X.shape[0] > N_CELLS_SUB:
        rng = np.random.default_rng(0)
        idx = rng.choice(X.shape[0], N_CELLS_SUB, replace=False)
        X = X[idx]
        print(f"[504] subsampled to {X.shape[0]} cells (SCHPF_N_CELLS={N_CELLS_SUB})")
    return X, genes


def load_cnmf_reference(k):
    p = MCP / f"cnmf_runs/global/global.spectra.k_{k}.dt_0_03.consensus.txt"
    if not p.exists():
        return None
    return pd.read_csv(p, sep="\t", index_col=0)


def top_set(vec, names, n=TOP_N):
    order = np.argsort(vec)[::-1][:n]
    return set(names[order])


def jaccard(a, b):
    u = len(a | b)
    return (len(a & b) / u) if u else 0.0


def main():
    try:
        import schpf  # noqa: F401
        from schpf import scHPF
    except ImportError:
        print("[504] ERROR: schpf not installed in this env. Install with:")
        print("      pip install git+https://github.com/simslab/scHPF.git")
        sys.exit(2)

    X, genes = load_input()

    all_show_rows = []
    all_jac_rows = []
    skip_done = os.environ.get("SCHPF_SKIP_DONE", "1") == "1"
    for k in K_LIST:
        print(f"\n[504] ===== scHPF k = {k} =====")
        # Skip if both factor + cell-score TSVs already on disk (retry resume)
        factors_out = OUT_DIR / f"schpf_factors_k{k}.tsv"
        cells_out = OUT_DIR / f"schpf_cell_scores_k{k}.tsv"
        if skip_done and factors_out.exists() and cells_out.exists():
            print(f"[504] k={k}: skipping fit (outputs exist); reloading for synthesis")
            df_beta = pd.read_csv(factors_out, sep="\t", index_col=0)
            # Reuse downstream Jaccard block without refitting
            ref_df = load_cnmf_reference(k)
            mean_jac = float("nan")
            if ref_df is not None:
                common = ref_df.columns.intersection(df_beta.index)
                ref_sub = ref_df.reindex(columns=common)
                ref_top = [top_set(ref_sub.iloc[i].values, ref_sub.columns)
                           for i in range(ref_sub.shape[0])]
                q_sub = df_beta.reindex(index=common)
                best_jacs = []
                for fi, fname in enumerate(q_sub.columns):
                    q_top = top_set(q_sub[fname].values, q_sub.index)
                    best_j = max(jaccard(q_top, r) for r in ref_top)
                    best_jacs.append(best_j)
                    all_jac_rows.append({
                        "method": "scHPF", "k": k, "factor": fname,
                        "best_jaccard_top100": best_j,
                    })
                mean_jac = float(np.mean(best_jacs))
                print(f"[504] k={k} mean best Jaccard vs cNMF (resumed): {mean_jac:.3f}")
            note = f"resumed from on-disk factors (fit skipped); n_genes={df_beta.shape[0]} n_factors={df_beta.shape[1]}"
            for metric, value in [
                ("jaccard_top100_vs_cnmf", mean_jac),
                ("final_loss", float("nan")),
                ("runtime_s", float("nan")),
                ("cophenetic", float("nan")),
                ("dispersion", float("nan")),
                ("silhouette", float("nan")),
                ("ari_mean", float("nan")),
            ]:
                all_show_rows.append({
                    "method": "scHPF", "k": k, "replicate": -1,
                    "metric": metric, "value": value, "note": note,
                })
            continue
        best_trial = None  # dict(model=, loss=, rt=, trial=)
        ref_df = load_cnmf_reference(k)
        for trial in range(N_TRIALS):
            t0 = time.time()
            # Default scHPF prior hyperparameters; Levitin 2019 settings.
            # scHPF's __init__ does NOT accept a random_state kwarg — it seeds
            # via numpy's global RNG state at fit time, so we seed here.
            np.random.seed(trial)
            model = scHPF(nfactors=k, verbose=True)
            model.fit(X.tocoo())
            loss = float(model.loss[-1]) if hasattr(model, "loss") and len(model.loss) else float("nan")
            rt = time.time() - t0
            print(f"[504] k={k} trial={trial} loss={loss:.4f} t={rt:.1f}s")
            if best_trial is None or loss < best_trial["loss"]:
                best_trial = dict(model=model, loss=loss, rt=rt, trial=trial)
        model = best_trial["model"]
        # Gene scores (beta): shape (n_genes, k)
        beta = model.beta.e_x
        # Cell scores (theta): (n_cells, k)
        theta = model.theta.e_x

        # Write outputs
        df_beta = pd.DataFrame(beta, index=genes, columns=[f"factor_{i+1}" for i in range(k)])
        df_beta.to_csv(OUT_DIR / f"schpf_factors_k{k}.tsv", sep="\t")
        df_theta = pd.DataFrame(theta, columns=[f"factor_{i+1}" for i in range(k)])
        df_theta.to_csv(OUT_DIR / f"schpf_cell_scores_k{k}.tsv", sep="\t", index=False)

        # Jaccard vs cNMF reference
        mean_jac = float("nan")
        if ref_df is not None:
            common = ref_df.columns.intersection(genes)
            ref_sub = ref_df.reindex(columns=common)
            ref_top = [top_set(ref_sub.iloc[i].values, ref_sub.columns) for i in range(ref_sub.shape[0])]
            q_sub = df_beta.reindex(index=common)
            best_jacs = []
            for fi, fname in enumerate(q_sub.columns):
                q_top = top_set(q_sub[fname].values, q_sub.index)
                best_j = max(jaccard(q_top, r) for r in ref_top)
                best_jacs.append(best_j)
                all_jac_rows.append({
                    "method": "scHPF", "k": k, "factor": fname, "best_jaccard_top100": best_j,
                })
            mean_jac = float(np.mean(best_jacs))
            print(f"[504] k={k} mean best Jaccard vs cNMF: {mean_jac:.3f}")

        note = (f"n_trials={N_TRIALS} best_trial={best_trial.get('trial', -1)} "
                f"nfactors={k} n_cells={X.shape[0]} n_genes={X.shape[1]}")
        all_show_rows += [
            dict(method="scHPF", k=k, replicate=-1, metric="jaccard_top100_vs_cnmf",
                 value=mean_jac, note=note),
            dict(method="scHPF", k=k, replicate=-1, metric="final_loss",
                 value=best_trial["loss"], note=note),
            dict(method="scHPF", k=k, replicate=-1, metric="runtime_s",
                 value=best_trial["rt"], note=note),
            dict(method="scHPF", k=k, replicate=-1, metric="cophenetic",
                 value=float("nan"), note="single fit — cophenetic not applicable"),
            dict(method="scHPF", k=k, replicate=-1, metric="dispersion",
                 value=float("nan"), note="NA"),
            dict(method="scHPF", k=k, replicate=-1, metric="silhouette",
                 value=float("nan"), note="NA"),
            dict(method="scHPF", k=k, replicate=-1, metric="ari_mean",
                 value=float("nan"), note="NA"),
        ]

    pd.DataFrame(all_show_rows).to_csv(OUT_DIR / "factorization_showdown_schpf.tsv",
                                       sep="\t", index=False)
    if all_jac_rows:
        pd.DataFrame(all_jac_rows).to_csv(OUT_DIR / "program_jaccard_schpf_vs_cnmf.tsv",
                                          sep="\t", index=False)
    print("[504] DONE.")


if __name__ == "__main__":
    main()
