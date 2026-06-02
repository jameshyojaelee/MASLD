#!/usr/bin/env python
"""
505_pyliger_inmf.py — pyliger iNMF benchmark (Welch 2019 PMID 31178122),
Python drop-in for Mission A1 Task 5. Replaces the 505_liger_inmf.R run that
could not compile `rliger` (hdf5r / RcppPlanc deps failed in celltype_bio).

Integrative NMF across the 7 datasets of atlas_cnmf_global.h5ad, with per-
factor dataset-specificity computed from the V matrices. k ∈ {10, 16}.

Outputs (in results_gpu_v2/mcp/benchmarks/):
  liger_factors_k{10,16}.tsv              — gene × factor shared W
  liger_dataset_factors_k{10,16}.tsv      — gene × factor × dataset (V long)
  liger_dataset_specificity_k{10,16}.tsv  — factor × dataset specificity
  program_jaccard_liger_vs_cnmf.tsv       — factor best top-100 Jaccard
  factorization_showdown_liger.tsv        — long-format showdown rows
"""
from __future__ import annotations

import gc
import json
import os
import sys
import time
import types
from pathlib import Path

# pyliger's package-level __init__ eagerly imports umap.plot, which
# transitively imports datashader → cudf. On CPU-only bigmem nodes the
# CUDA driver probe in cudf raises `CUDARuntimeError` and the whole
# import tree crashes. We only need pyliger's numpy iNMF, not its UMAP
# plotting helpers, so stub umap.plot with no-ops BEFORE importing pyliger.
_up = types.ModuleType("umap.plot")
_up.points = lambda *a, **k: None
_up.connectivity = lambda *a, **k: None
_up.diagnostic = lambda *a, **k: None
sys.modules["umap.plot"] = _up

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmread

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
MCP = ROOT / "Analysis/SingleCell/results_gpu_v2/mcp"
STAGING_DIR = MCP / "liger_staging"
OUT_DIR = MCP / "benchmarks"
OUT_DIR.mkdir(parents=True, exist_ok=True)
CNMF_DIR = MCP / "cnmf_runs/global"

K_LIST = [int(x) for x in os.environ.get("LIGER_K_LIST", "10,16").split(",")]
TOP_N = 100
N_ITERS = int(os.environ.get("LIGER_N_ITERS", 30))


def load_datasets():
    """Build an AnnData per dataset from the staged .mtx + barcode/gene files."""
    manifest = json.loads((STAGING_DIR / "manifest.json").read_text())
    adatas = {}
    for d in manifest["datasets"]:
        mtx = mmread(d["mtx_path"]).tocsr()   # genes × cells in staging file
        genes = [g.strip() for g in open(d["genes_path"])]
        # Prefix barcodes with dataset name to guarantee cross-dataset uniqueness
        # (pyliger raises if any cell name collides across datasets)
        raw_barcodes = [b.strip() for b in open(d["barcodes_path"])]
        barcodes = [f"{d['name']}::{b}" for b in raw_barcodes]
        # pyliger expects cells × genes
        X = mtx.T.tocsr().astype(np.float32)
        a = ad.AnnData(X=X)
        a.var_names = pd.Index(genes)
        a.obs_names = pd.Index(barcodes)
        # pyliger's _create_liger_matrix() checks `adata.obs.index.name` and
        # `adata.var.index.name` (not the Index values) — defaults to None,
        # so we must assign them explicitly or createLiger raises
        # "Raw data must have both row (cell) and column (gene) names."
        a.obs.index.name = "cell"
        a.var.index.name = "gene"
        a.obs_names_make_unique()
        a.var_names_make_unique()
        a.uns["sample_name"] = d["name"]
        adatas[d["name"]] = a
        print(f"[505py]   {d['name']}: {X.shape}")
    return adatas


def top_set(vec, names, n=TOP_N):
    order = np.argsort(vec)[::-1][:n]
    return set(names[order])


def jaccard(a, b):
    u = len(a | b)
    return (len(a & b) / u) if u else 0.0


def load_cnmf_ref(k):
    p = CNMF_DIR / f"global.spectra.k_{k}.dt_0_03.consensus.txt"
    if not p.exists():
        return None
    return pd.read_csv(p, sep="\t", index_col=0)


def main():
    try:
        import pyliger  # noqa: F401
    except ImportError:
        print("[505py] pyliger not installed — run `pip install pyliger` first.")
        sys.exit(2)

    adatas = load_datasets()
    print(f"[505py] {len(adatas)} datasets loaded")

    all_show_rows = []
    all_jac_rows = []

    for k in K_LIST:
        print(f"\n[505py] ===== pyliger iNMF k = {k} =====")
        # Rebuild the Liger object each iteration (pyliger mutates in place)
        adatas_copy = {n: a.copy() for n, a in adatas.items()}
        lig = pyliger.create_liger(list(adatas_copy.values()), remove_missing=False)
        pyliger.normalize(lig)
        # Use union of per-dataset varnames as HVG (we already HVG-subset upstream).
        # Seed pyliger's selected genes with the staging HVG list.
        union_hvg = sorted(set().union(*(set(a.var_names.tolist()) for a in adatas_copy.values())))
        try:
            pyliger.select_genes(lig, var_thresh=0.1, do_plot=False)
        except Exception as e:
            print(f"[505py] select_genes fell back to union HVG ({e})")
        # Ensure all HVG are in the selected set
        lig.var_genes = union_hvg
        pyliger.scale_not_center(lig)

        t0 = time.time()
        # pyliger 0.2.x: optimize_ALS(ligerex, k, ...)
        pyliger.optimize_ALS(lig, k=k, nrep=1, thresh=1e-6, max_iters=N_ITERS,
                             rand_seed=42, value_lambda=5.0)
        rt = time.time() - t0
        print(f"[505py] k={k} trained in {rt:.1f} s")

        # pyliger stores W/V on each sub-AnnData's varm:
        #   lig.adata_list[i].varm["W"]  = gene × factor (shared, duplicated)
        #   lig.adata_list[i].varm["V"]  = gene × factor (dataset-specific)
        #   lig.adata_list[i].obsm["H"]  = cell × factor
        factor_names = [f"factor_{i+1}" for i in range(k)]
        # Use first dataset's varm["W"] as the canonical shared W (all copies are equal).
        first = lig.adata_list[0]
        W_mat = np.asarray(first.varm["W"])
        genes0 = list(first.var_names)
        if W_mat.shape[1] != k and W_mat.shape[0] == k:
            W_mat = W_mat.T
        W_df = pd.DataFrame(W_mat, index=genes0, columns=factor_names)
        W_df.to_csv(OUT_DIR / f"liger_factors_k{k}.tsv", sep="\t")
        w_norms = np.linalg.norm(W_df.values, axis=0)  # per-factor ||W||

        spec_rows = []
        v_long_rows = []
        for sub in lig.adata_list:
            dname = sub.uns.get("sample_name") if "sample_name" in sub.uns else None
            if dname is None:
                # sample_names is a list aligned with adata_list order
                dname = lig.sample_names[lig.adata_list.index(sub)]
            V = np.asarray(sub.varm.get("V"))
            if V is None or V.size == 0:
                print(f"[505py] WARN: {dname} has no V matrix (skipped)")
                continue
            if V.shape[1] != k and V.shape[0] == k:
                V = V.T
            V_df = pd.DataFrame(V, index=list(sub.var_names), columns=factor_names)
            v_norms = np.linalg.norm(V_df.values, axis=0)
            for fi, fname in enumerate(factor_names):
                spec = v_norms[fi] / (w_norms[fi] + v_norms[fi] + 1e-12)
                spec_rows.append({
                    "k": k, "factor": fname, "dataset": dname,
                    "w_norm": float(w_norms[fi]),
                    "v_norm": float(v_norms[fi]),
                    "dataset_specificity": float(spec),
                })
            for fi, fname in enumerate(factor_names):
                for gi, g in enumerate(V_df.index):
                    v_long_rows.append({
                        "dataset": dname, "factor": fname,
                        "gene": g, "dataset_weight": float(V_df.iloc[gi, fi]),
                    })

        spec_df = pd.DataFrame(spec_rows)
        spec_df.to_csv(OUT_DIR / f"liger_dataset_specificity_k{k}.tsv",
                       sep="\t", index=False)

        # Write the long-format V table in chunks to save memory
        v_long_df = pd.DataFrame(v_long_rows)
        v_long_df.to_csv(OUT_DIR / f"liger_dataset_factors_k{k}.tsv",
                         sep="\t", index=False)
        del v_long_df, v_long_rows
        gc.collect()

        # Per-factor max specificity → dataset-loading flag
        max_spec = spec_df.groupby("factor")["dataset_specificity"].max().reset_index()
        max_spec["dataset_loading_flag"] = max_spec["dataset_specificity"] > 0.5
        print("[505py] per-factor max dataset-specificity:")
        print(max_spec.to_string(index=False))
        n_dataset_loading = int(max_spec["dataset_loading_flag"].sum())

        # Jaccard vs cNMF
        ref = load_cnmf_ref(k)
        mean_jac = float("nan")
        if ref is not None:
            common = ref.columns.intersection(W_df.index)
            ref_sub = ref.reindex(columns=common)
            ref_top = [top_set(ref_sub.iloc[i].values, ref_sub.columns)
                       for i in range(ref_sub.shape[0])]
            q_sub = W_df.reindex(index=common)
            for fi, fname in enumerate(factor_names):
                q_top = top_set(q_sub[fname].values, q_sub.index)
                best = max(jaccard(q_top, r) for r in ref_top)
                spec_row = max_spec.loc[max_spec["factor"] == fname].iloc[0]
                all_jac_rows.append({
                    "method": "LIGER(pyliger)", "k": k, "factor": fname,
                    "best_jaccard_top100": best,
                    "max_dataset_specificity": float(spec_row["dataset_specificity"]),
                    "dataset_loading_flag": bool(spec_row["dataset_loading_flag"]),
                })
            mean_jac = float(np.mean([r["best_jaccard_top100"] for r in all_jac_rows if r["method"] == "LIGER(pyliger)" and r["k"] == k]))
            print(f"[505py] k={k} mean best Jaccard (LIGER W vs cNMF): {mean_jac:.3f}")

        note = (f"pyliger optimize_ALS; n_datasets={len(lig.sample_names)}; "
                f"k={k}; max_iters={N_ITERS}; seed=42; lambda=5.0; rt={rt:.0f}s")
        for metric, value in [
            ("cophenetic", float("nan")),
            ("dispersion", float("nan")),
            ("silhouette", float("nan")),
            ("ari_mean", float("nan")),
            ("jaccard_top100_vs_cnmf", mean_jac),
            ("runtime_s", rt),
            ("n_dataset_loading_factors", float(n_dataset_loading)),
        ]:
            all_show_rows.append({
                "method": "LIGER(pyliger)", "k": k, "replicate": -1,
                "metric": metric, "value": value, "note": note,
            })

    pd.DataFrame(all_show_rows).to_csv(OUT_DIR / "factorization_showdown_liger.tsv",
                                       sep="\t", index=False)
    if all_jac_rows:
        pd.DataFrame(all_jac_rows).to_csv(
            OUT_DIR / "program_jaccard_liger_vs_cnmf.tsv", sep="\t", index=False)
    print("[505py] DONE.")


if __name__ == "__main__":
    main()
