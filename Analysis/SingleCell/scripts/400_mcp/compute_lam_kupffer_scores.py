#!/usr/bin/env python
# ============================================================================
# compute_lam_kupffer_scores.py
# Compute per-macrophage Kupffer and LAM identity scores from
# Macrophages_subset.h5ad. Replaces the noisy "sign-flip consensus pseudotime"
# task: consensus pseudotime captures macrophage maturation (both Kupffer AND
# LAM markers negative with pt), not the Kupffer->LAM identity axis. We need a
# direct identity score.
#
# Inputs:
#   Analysis/SingleCell/results_gpu_v2/pseudotime/Macrophages_subset.h5ad
#   Analysis/SingleCell/results_gpu_v2/macrophage_trajectory/
#       macrophage_marker_panels_pseudotime.csv (canonical panel definitions)
#
# Outputs:
#   Analysis/SingleCell/results_gpu_v2/macrophage_trajectory/
#       macrophage_KL_axis_per_cell.csv  (cell, kupffer, lam, kl_axis)
#       macrophage_KL_axis_summary.csv   (per-stage means + Spearman vs stage)
#       macrophage_KL_receptor_rho.csv   (recomputed receptor rho vs kl_axis
#           for hep->mac LR pairs, replaces noisy receptor_rho-vs-pseudotime)
# ============================================================================

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scipy.stats import spearmanr

BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
PT_DIR    = BASE / "Analysis/SingleCell/results_gpu_v2/pseudotime"
TRAJ_DIR  = BASE / "Analysis/SingleCell/results_gpu_v2/macrophage_trajectory"
H5AD      = PT_DIR / "Macrophages_subset.h5ad"
PANELS    = TRAJ_DIR / "macrophage_marker_panels_pseudotime.csv"
LR_FILE   = TRAJ_DIR / "hep_to_mac_ligands_x_pseudotime_receptors.csv"
OUT_CELLS = TRAJ_DIR / "macrophage_KL_axis_per_cell.csv"
OUT_SUMMARY = TRAJ_DIR / "macrophage_KL_axis_summary.csv"
OUT_LR    = TRAJ_DIR / "macrophage_KL_receptor_rho.csv"

# Curated marker panels (literature-anchored: Ramachandran 2019, Guilliams 2022,
# Andrews 2022, Remmerie 2020, Jaitin 2019). The CSV-derived rho table mixes
# weak and strong markers; for an identity score we use the strongest, most
# specific markers per panel.
KUPFFER = ["MARCO", "CD5L", "VSIG4", "TIMD4", "LYVE1", "FCN1", "CETP", "CLEC4F"]
LAM     = ["TREM2", "GPNMB", "CD9", "SPP1", "LIPA", "LGALS3", "CTSD", "PLIN2",
           "FABP5", "CTSL"]


def score_panel(adata, gene_list, score_name, n_ctrl=50):
    """sc.tl.score_genes wrapper that handles missing genes gracefully."""
    keep = [g for g in gene_list if g in adata.var_names]
    missing = sorted(set(gene_list) - set(keep))
    if missing:
        print(f"  [{score_name}] missing from var_names: {missing}")
    if not keep:
        adata.obs[score_name] = 0.0
        return []
    sc.tl.score_genes(
        adata, gene_list=keep, score_name=score_name,
        ctrl_size=n_ctrl, n_bins=25, random_state=42,
        copy=False, use_raw=False)
    return keep


def main():
    print(f"Loading {H5AD} ...")
    adata = ad.read_h5ad(H5AD)
    print(f"  shape = {adata.shape} (cells x genes)")
    print(f"  obs cols = {list(adata.obs.columns)[:30]}")
    print(f"  X dtype  = {adata.X.dtype}")
    print(f"  layers   = {list(adata.layers.keys())}")

    # If counts layer exists, prefer it for scoring; else use X
    if "counts" in adata.layers:
        # Move counts to .X temporarily for scoring on log-normalized data; if
        # X is already log1p-normalized, scoring proceeds directly.
        pass

    print("\nScoring panels ...")
    used_kc = score_panel(adata, KUPFFER, "kupffer_score")
    used_lam = score_panel(adata, LAM, "lam_score")
    print(f"  Used Kupffer markers: {used_kc}")
    print(f"  Used LAM markers:     {used_lam}")

    # Composite identity axis: positive = LAM-like, negative = Kupffer-like
    adata.obs["kl_axis"] = adata.obs["lam_score"] - adata.obs["kupffer_score"]

    # ---- Per-cell output -----------------------------------------------------
    out = adata.obs[[
        "sample", "dataset", "condition",
        "kupffer_score", "lam_score", "kl_axis"]].copy()
    if "disease_stage_coarse" in adata.obs.columns:
        out["disease_stage_coarse"] = adata.obs["disease_stage_coarse"]
    if "fibrosis_stage" in adata.obs.columns:
        out["fibrosis_stage"] = adata.obs["fibrosis_stage"]
    out.index.name = "cell"
    out.to_csv(OUT_CELLS)
    print(f"\nWrote per-cell scores -> {OUT_CELLS}  (n={len(out)})")

    # ---- Summary by stage ----------------------------------------------------
    stage_col = ("disease_stage_coarse"
                 if "disease_stage_coarse" in out.columns
                 else "condition")
    summary = (out.groupby(stage_col)
                  .agg(n_cells=("kl_axis", "size"),
                       kupffer_mean=("kupffer_score", "mean"),
                       lam_mean=("lam_score", "mean"),
                       kl_axis_mean=("kl_axis", "mean"),
                       kl_axis_sd=("kl_axis", "std"))
                  .reset_index())
    summary.to_csv(OUT_SUMMARY, index=False)
    print(f"Wrote stage summary -> {OUT_SUMMARY}")
    print(summary)

    # ---- Recompute receptor rho vs kl_axis ----------------------------------
    if not LR_FILE.exists():
        print(f"WARN: {LR_FILE} not found; skipping receptor-rho recomputation.")
        return

    lr = pd.read_csv(LR_FILE)
    receptors = sorted(set(lr["receptor"].dropna()))
    # Some receptors are heteromers like "B2M_FCGRT"; split on '_' and keep
    # the primary subunit (we'll record both in case)
    receptor_genes = []
    for r in receptors:
        for sub in str(r).split("_"):
            if sub in adata.var_names:
                receptor_genes.append((r, sub))
                break

    print(f"\nReceptors with mappable genes: {len(receptor_genes)} / {len(receptors)}")

    # Pull receptor expression (cells x receptor_genes)
    rec_genes_unique = sorted(set(g for _, g in receptor_genes))
    expr = adata[:, rec_genes_unique].X
    if hasattr(expr, "toarray"):
        expr = expr.toarray()

    kl_vec = adata.obs["kl_axis"].values
    rows = []
    for label, gene in receptor_genes:
        idx = rec_genes_unique.index(gene)
        x = expr[:, idx]
        # Skip if essentially constant
        if np.nanstd(x) < 1e-9:
            rows.append({"receptor": label, "gene": gene,
                         "rho_kl_axis": np.nan, "p_kl_axis": np.nan,
                         "rho_disease": np.nan, "p_disease": np.nan,
                         "n_cells": int(len(x))})
            continue
        rho_kl, p_kl = spearmanr(x, kl_vec, nan_policy="omit")
        # Disease vs healthy as binary contrast (fallback if stage not numeric)
        cond = adata.obs["condition"].astype(str).values
        binary = np.where(np.isin(cond, ["MASLD", "MASH", "Cirrhotic"]), 1, 0)
        rho_dis, p_dis = spearmanr(x, binary, nan_policy="omit")
        rows.append({"receptor": label, "gene": gene,
                     "rho_kl_axis": float(rho_kl), "p_kl_axis": float(p_kl),
                     "rho_disease": float(rho_dis), "p_disease": float(p_dis),
                     "n_cells": int(len(x))})

    out_rho = pd.DataFrame(rows)
    out_rho.to_csv(OUT_LR, index=False)
    print(f"Wrote receptor-vs-kl_axis -> {OUT_LR}  (n={len(out_rho)})")
    print(f"  positive rho (LAM-side): {(out_rho['rho_kl_axis'] > 0).sum()}")
    print(f"  negative rho (KC-side):  {(out_rho['rho_kl_axis'] < 0).sum()}")


if __name__ == "__main__":
    main()
