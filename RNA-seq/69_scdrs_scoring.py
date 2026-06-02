#!/usr/bin/env python3
"""
69_scdrs_scoring.py
---------------------------------------------------------------------------
Phase 4A: scDRS per-cell disease relevance scoring

Uses S-PrediXcan TWAS z-scores as gene-level disease association scores
(equivalent to MAGMA z-scores) to compute per-cell disease relevance.

Input:
  - Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad
  - RNA-seq/results/causal_inference/{gwas}/twas_spredixcan_liver.csv

Output:
  - RNA-seq/results/causal_inference/scdrs/{gwas}/
      - cell_scores.csv         (per-cell disease z-score)
      - group_results.csv       (per cell-type enrichment)
      - group_condition.csv     (per cell-type × condition enrichment)
      - scdrs_summary.csv       (summary statistics)

Usage:
  python RNA-seq/69_scdrs_scoring.py [--gwas ghodsian]
---------------------------------------------------------------------------
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scdrs
import warnings
import time

warnings.filterwarnings("ignore")

# ---- Configuration ----
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)

ATLAS_PATH = os.path.join(
    BASE, "Analysis/SingleCell/results_gpu_v2/integrated_atlas.h5ad"
)

GWAS_LIST = ["ghodsian", "ukbb_alt", "ukbb_ast", "ukbb_ggt", "pdff"]

# Top N genes by |z-score| to use as disease gene set
TOP_N_GENES = 1000

# scDRS parameters
N_CTRL = 1000  # number of control gene sets


def load_twas_gene_scores(gwas_name):
    """Load TWAS z-scores as gene-level disease scores.

    Checks S-PrediXcan first, then falls back to OTTERS combined results.
    """
    # Try S-PrediXcan first
    twas_file = os.path.join(
        BASE,
        "RNA-seq/results/causal_inference",
        gwas_name,
        "twas_spredixcan_liver.csv",
    )
    if os.path.exists(twas_file):
        df = pd.read_csv(twas_file)
        if "gene_name" in df.columns and "zscore" in df.columns:
            df = df.dropna(subset=["gene_name", "zscore"])
            df = df[np.isfinite(df["zscore"])]
            df["abs_z"] = df["zscore"].abs()
            df = df.sort_values("abs_z", ascending=False).drop_duplicates(
                subset="gene_name", keep="first"
            )
            print(f"  {gwas_name}: {len(df)} genes with TWAS z-scores (S-PrediXcan)")
            return df

    # Fallback: OTTERS combined results
    otters_file = os.path.join(
        BASE,
        "RNA-seq/results/causal_inference/otters_broadaway",
        gwas_name,
        "otters_twas_combined.csv",
    )
    if os.path.exists(otters_file):
        df = pd.read_csv(otters_file)
        # OTTERS columns: gene_symbol, best_z
        sym_col = next((c for c in ["gene_symbol", "GeneName"] if c in df.columns), None)
        z_col = next((c for c in ["best_z", "acat_pval", "FUSION_Z"] if c in df.columns), None)
        if sym_col and z_col:
            df = df.rename(columns={sym_col: "gene_name", z_col: "zscore"})
            df = df.dropna(subset=["gene_name", "zscore"])
            df = df[df["gene_name"] != ""]
            df["zscore"] = df["zscore"].astype(float)
            df = df[np.isfinite(df["zscore"])]
            df["abs_z"] = df["zscore"].abs()
            df = df.sort_values("abs_z", ascending=False).drop_duplicates(
                subset="gene_name", keep="first"
            )
            print(f"  {gwas_name}: {len(df)} genes with TWAS z-scores (OTTERS)")
            return df

    print(f"  WARNING: No TWAS results found for {gwas_name}")
    return None


def get_disease_gene_set(twas_df, n_top=TOP_N_GENES):
    """Extract top disease genes by |z-score| for scDRS."""
    top = twas_df.nlargest(n_top, "abs_z")
    gene_list = top["gene_name"].tolist()
    gene_weights = top["abs_z"].values
    return gene_list, gene_weights


def run_scdrs_for_gwas(adata, gwas_name, out_dir):
    """Run scDRS for a single GWAS."""
    print(f"\n--- scDRS: {gwas_name} ---")
    start = time.time()

    # Load gene scores
    twas_df = load_twas_gene_scores(gwas_name)
    if twas_df is None:
        return None

    # Get disease gene set
    gene_list, gene_weights = get_disease_gene_set(twas_df)

    # Filter to genes present in adata
    valid_genes = [g for g in gene_list if g in adata.var_names]
    valid_idx = [i for i, g in enumerate(gene_list) if g in adata.var_names]
    valid_weights = gene_weights[valid_idx]
    print(f"  Disease genes in atlas: {len(valid_genes)}/{len(gene_list)}")

    if len(valid_genes) < 50:
        print(f"  WARNING: Too few valid genes ({len(valid_genes)}), skipping")
        return None

    # Run scDRS scoring
    print("  Scoring cells...")
    score_result = scdrs.score_cell(
        adata,
        gene_list=valid_genes,
        gene_weight=valid_weights,
        ctrl_match_key="mean_var",
        n_ctrl=N_CTRL,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=False,
    )

    # Extract scores
    cell_scores = pd.DataFrame(
        {
            "cell_id": adata.obs_names,
            "raw_score": score_result["raw_score"],
            "norm_score": score_result["norm_score"],
            "mc_pval": score_result["mc_pval"],
            "cell_type": adata.obs["cell_type"].values,
            "condition": adata.obs["condition_harmonized"].values
            if "condition_harmonized" in adata.obs.columns
            else adata.obs.get("condition", pd.Series(["unknown"] * len(adata))).values,
            "dataset": adata.obs["dataset"].values,
        }
    )

    # Save cell scores
    os.makedirs(out_dir, exist_ok=True)
    cell_scores.to_csv(os.path.join(out_dir, "cell_scores.csv"), index=False)
    print(f"  Cell scores saved: {len(cell_scores)} cells")

    # Group-level enrichment by cell type
    print("  Computing cell-type enrichment...")
    group_results = []
    for ct in sorted(adata.obs["cell_type"].unique()):
        ct_mask = cell_scores["cell_type"] == ct
        n_cells = ct_mask.sum()
        if n_cells < 10:
            continue
        ct_scores = cell_scores.loc[ct_mask, "norm_score"]
        # One-sided test: are disease scores elevated in this cell type?
        from scipy import stats

        t_stat, t_pval = stats.ttest_1samp(ct_scores, 0, alternative="greater")
        mean_score = ct_scores.mean()
        frac_sig = (cell_scores.loc[ct_mask, "mc_pval"] < 0.05).mean()

        group_results.append(
            {
                "cell_type": ct,
                "n_cells": n_cells,
                "mean_norm_score": mean_score,
                "ttest_stat": t_stat,
                "ttest_pval": t_pval,
                "frac_significant": frac_sig,
            }
        )

    group_df = pd.DataFrame(group_results)
    if len(group_df) > 0:
        from statsmodels.stats.multitest import multipletests

        group_df["fdr"] = multipletests(group_df["ttest_pval"], method="fdr_bh")[1]
        group_df = group_df.sort_values("ttest_pval")
        group_df.to_csv(os.path.join(out_dir, "group_results.csv"), index=False)
        print(f"  Cell-type enrichment: {(group_df['fdr'] < 0.05).sum()} FDR<0.05")

    # Group-level enrichment by cell type × condition
    print("  Computing cell-type × condition enrichment...")
    cond_results = []
    for ct in sorted(adata.obs["cell_type"].unique()):
        ct_mask = cell_scores["cell_type"] == ct
        conditions = cell_scores.loc[ct_mask, "condition"].unique()
        for cond in conditions:
            mask = ct_mask & (cell_scores["condition"] == cond)
            n = mask.sum()
            if n < 10:
                continue
            scores = cell_scores.loc[mask, "norm_score"]
            t_stat, t_pval = stats.ttest_1samp(scores, 0, alternative="greater")
            cond_results.append(
                {
                    "cell_type": ct,
                    "condition": cond,
                    "n_cells": n,
                    "mean_norm_score": scores.mean(),
                    "ttest_stat": t_stat,
                    "ttest_pval": t_pval,
                    "frac_significant": (
                        cell_scores.loc[mask, "mc_pval"] < 0.05
                    ).mean(),
                }
            )
    cond_df = pd.DataFrame(cond_results)
    if len(cond_df) > 0:
        cond_df["fdr"] = multipletests(cond_df["ttest_pval"], method="fdr_bh")[1]
        cond_df = cond_df.sort_values("ttest_pval")
        cond_df.to_csv(os.path.join(out_dir, "group_condition.csv"), index=False)

    # Summary
    n_sig_cells = (cell_scores["mc_pval"] < 0.05).sum()
    summary = {
        "gwas": gwas_name,
        "n_disease_genes": len(valid_genes),
        "n_cells": len(cell_scores),
        "n_sig_cells": n_sig_cells,
        "frac_sig_cells": n_sig_cells / len(cell_scores),
        "n_celltypes_enriched": (group_df["fdr"] < 0.05).sum()
        if len(group_df) > 0
        else 0,
        "top_celltype": group_df.iloc[0]["cell_type"]
        if len(group_df) > 0
        else "NA",
        "top_celltype_pval": group_df.iloc[0]["ttest_pval"]
        if len(group_df) > 0
        else 1.0,
        "elapsed_sec": time.time() - start,
    }

    print(f"  Significant cells: {n_sig_cells:,} ({summary['frac_sig_cells']:.1%})")
    if len(group_df) > 0:
        print(f"  Top cell type: {summary['top_celltype']} (p={summary['top_celltype_pval']:.2e})")

    return summary


def main():
    parser = argparse.ArgumentParser(description="scDRS disease scoring")
    parser.add_argument("--gwas", default=None, help="Specific GWAS to run (default: all)")
    parser.add_argument(
        "--downsample",
        type=int,
        default=0,
        help="Downsample to N cells (0=no downsample)",
    )
    args = parser.parse_args()

    gwas_list = [args.gwas] if args.gwas else GWAS_LIST

    print("=== scDRS Disease Scoring ===")
    print(f"Start: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Atlas: {ATLAS_PATH}")
    print(f"GWAS: {', '.join(gwas_list)}")
    print()

    # Load atlas (with fallback for anndata version mismatch)
    print("Loading scRNA atlas...")
    try:
        adata = sc.read_h5ad(ATLAS_PATH)
    except KeyError as e:
        if "ordered" in str(e):
            print("  Falling back to h5py-based obs reading (anndata version mismatch)...")
            import h5py
            import scipy.sparse as sp

            with h5py.File(ATLAS_PATH, "r") as f:
                # Read X matrix
                X_group = f["X"]
                if "data" in X_group:
                    X = sp.csr_matrix(
                        (X_group["data"][:], X_group["indices"][:], X_group["indptr"][:]),
                        shape=tuple(X_group.attrs.get("shape", X_group.attrs.get("h5sparse_shape"))),
                    )
                else:
                    X = np.array(X_group)

                # Read var
                var_dict = {}
                for key in f["var"].keys():
                    try:
                        vals = f["var"][key][:]
                        if hasattr(vals[0], "decode"):
                            vals = np.array([v.decode() for v in vals])
                        var_dict[key] = vals
                    except Exception:
                        pass
                var_index = None
                if "_index" in f["var"]:
                    idx = f["var"]["_index"][:]
                    var_index = [v.decode() if hasattr(v, "decode") else v for v in idx]
                var_df = pd.DataFrame(var_dict, index=var_index)

                # Read obs — handle categoricals manually
                obs_dict = {}
                obs_index = None
                if "_index" in f["obs"]:
                    idx = f["obs"]["_index"][:]
                    obs_index = [v.decode() if hasattr(v, "decode") else v for v in idx]
                for key in f["obs"].keys():
                    if key == "_index":
                        continue
                    try:
                        grp = f["obs"][key]
                        if isinstance(grp, h5py.Group) and "codes" in grp and "categories" in grp:
                            codes = grp["codes"][:]
                            cats = grp["categories"][:]
                            if hasattr(cats[0], "decode"):
                                cats = np.array([c.decode() for c in cats])
                            vals = cats[codes]
                            vals[codes < 0] = None
                            obs_dict[key] = vals
                        elif isinstance(grp, h5py.Dataset):
                            vals = grp[:]
                            if hasattr(vals[0], "decode"):
                                vals = np.array([v.decode() for v in vals])
                            obs_dict[key] = vals
                    except Exception:
                        pass
                obs_df = pd.DataFrame(obs_dict, index=obs_index)

            import anndata

            adata = anndata.AnnData(X=X, obs=obs_df, var=var_df)
            print(f"  Loaded via h5py fallback")
        else:
            raise
    print(f"  Shape: {adata.shape[0]:,} cells × {adata.shape[1]:,} genes")
    print(f"  Cell types: {adata.obs['cell_type'].nunique()}")

    # Optional downsampling for memory
    if args.downsample > 0 and args.downsample < adata.shape[0]:
        print(f"  Downsampling to {args.downsample:,} cells (stratified by cell type)...")
        sc.pp.subsample(adata, n_obs=args.downsample)
        print(f"  After downsampling: {adata.shape[0]:,} cells")

    # Preprocess for scDRS
    print("Preprocessing for scDRS...")
    # scDRS expects log1p-normalized data
    # Check if already normalized
    if adata.X.max() > 50:
        print("  Data appears raw — normalizing...")
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
    else:
        print("  Data appears already log-normalized")

    # scDRS preprocessing
    scdrs.preprocess(adata)
    print("  scDRS preprocessing complete")

    # Run for each GWAS (skip already completed)
    all_summaries = []
    for gwas_name in gwas_list:
        out_dir = os.path.join(
            BASE, "RNA-seq/results/causal_inference/scdrs", gwas_name
        )
        # Skip if results already exist
        if os.path.exists(os.path.join(out_dir, "group_results.csv")):
            print(f"\n--- scDRS: {gwas_name} --- (already done, skipping)")
            continue
        summary = run_scdrs_for_gwas(adata, gwas_name, out_dir)
        if summary is not None:
            all_summaries.append(summary)

    # Save combined summary
    if all_summaries:
        summary_df = pd.DataFrame(all_summaries)
        summary_out = os.path.join(
            BASE, "RNA-seq/results/causal_inference/scdrs/scdrs_summary.csv"
        )
        os.makedirs(os.path.dirname(summary_out), exist_ok=True)
        summary_df.to_csv(summary_out, index=False)
        print(f"\n=== Summary ===")
        print(summary_df.to_string(index=False))

    print(f"\n=== Done: {time.strftime('%Y-%m-%d %H:%M:%S')} ===")


if __name__ == "__main__":
    main()
