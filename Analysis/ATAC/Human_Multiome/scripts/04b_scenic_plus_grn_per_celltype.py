#!/usr/bin/env python3
"""
04b_scenic_plus_grn_per_celltype.py — Per-cell-type gene-activity-based GRN

Generalizes Script 04b (`04b_scenic_grn_from_activity.py`, hepatocyte only) to
arbitrary cell-type subsets so that disease regulons can be built for
Macrophage+Kupffer, Stellate, Endothelial+LSEC, and Cholangiocyte populations
in parallel. The pipeline is otherwise identical:

  1. Load SnapATAC2 processed h5ad
  2. Subset to cells whose `cell_type` is in --target-labels (comma list,
     merged into one analysis subset, e.g. "Macrophage,Kupffer_Cell")
  3. Stratified subsample to --n-cells across `condition` levels
  4. Compute gene activity from the AnnDataSet (positional indices)
  5. Parse GENCODE GTF for TSS, identify expressed TFs
  6. Vectorized TF-enhancer and enhancer-gene Spearman correlations
  7. Assemble regulons, score AUCell-like activity, MASLD vs Normal Wilcoxon
  8. Export 4 CSVs named by --celltype-name:
       {NAME}_regulons.csv
       {NAME}_enhancer_gene_links.csv
       {NAME}_regulon_activity_scores.csv
       {NAME}_disease_regulons.csv

The original `04b_scenic_grn_from_activity.py` is left untouched. This script
imports its helper functions (TSS parsing, rank-matrix correlations, regulon
assembly, scoring, export) and only overrides the cell-type subset logic and
output naming.

Environment:
  micromamba activate snapatac2
  export PYTHONNOUSERSITE=1

Usage (single CT or merged group):
  python scripts/04b_scenic_plus_grn_per_celltype.py \
      --target-labels "Macrophage,Kupffer_Cell" \
      --celltype-name Macrophage

Sbatch array:
  sbatch scripts/run_04b_scenic_per_ct.sbatch

Author: MASLD-Atlas pipeline (B2 of the ATAC improvement plan)
"""

from __future__ import annotations

import argparse
import gc
import logging
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

# Reuse all helpers from Script 04b. Importing the module keeps a single source
# of truth for correlation math, TSS parsing, and regulon assembly.
_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

import importlib.util
_M_NAME = "scenic_grn_from_activity"
_M_PATH = _SCRIPT_DIR / "04b_scenic_grn_from_activity.py"
_spec = importlib.util.spec_from_file_location(_M_NAME, _M_PATH)
_grn = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_grn)

KNOWN_TFS = _grn.KNOWN_TFS
parse_tss_from_gtf = _grn.parse_tss_from_gtf
parse_peak_coords = _grn.parse_peak_coords
build_peak_gene_proximity = _grn.build_peak_gene_proximity
_rank_matrix = _grn._rank_matrix
compute_tf_enhancer_correlations = _grn.compute_tf_enhancer_correlations
compute_enhancer_gene_correlations = _grn.compute_enhancer_gene_correlations
assemble_regulons = _grn.assemble_regulons
score_regulon_activity = _grn.score_regulon_activity


log = logging.getLogger("scenic_grn_per_ct")


# ---------------------------------------------------------------------------
# Cell-type-aware data loaders (analogue of `load_and_get_hepatocyte_indices`
# and `compute_gene_activity` from Script 04b)
# ---------------------------------------------------------------------------
def load_and_get_celltype_indices(
    h5ad_path: str,
    target_labels: list[str],
    celltype_name: str,
    n_cells: int = 15000,
    seed: int = 42,
) -> tuple[ad.AnnData, np.ndarray]:
    """Load SnapATAC2 processed h5ad and subset to cells whose cell_type is in
    `target_labels`. Returns (subsampled_adata, global_indices).
    """
    log.info("Loading ATAC data: %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    log.info("  Full data: %d cells x %d peaks", adata.n_obs, adata.n_vars)

    if "cell_type" not in adata.obs.columns:
        raise ValueError("ATAC h5ad missing required 'cell_type' obs column")

    obs_ct = adata.obs["cell_type"].astype(str)
    present = sorted(set(obs_ct.unique()))
    missing = [lbl for lbl in target_labels if lbl not in present]
    if missing:
        log.warning(
            "Requested labels not present in h5ad: %s (present: %s)",
            missing, present,
        )

    mask = obs_ct.isin(target_labels).to_numpy()
    ct_global_idx = np.where(mask)[0]
    log.info(
        "  Cell-type subset '%s' (labels=%s): %d cells",
        celltype_name, target_labels, len(ct_global_idx),
    )
    if len(ct_global_idx) < 100:
        raise ValueError(
            f"Too few cells ({len(ct_global_idx)}) for subset '{celltype_name}'"
        )

    if len(ct_global_idx) > n_cells:
        rng = np.random.RandomState(seed)
        conditions = adata.obs["condition"].astype(str).to_numpy()[ct_global_idx]
        unique_conds = np.unique(conditions)
        total = len(ct_global_idx)
        selected_local = []
        for c in unique_conds:
            c_local_idx = np.where(conditions == c)[0]
            n_sample = max(1, int(n_cells * len(c_local_idx) / total))
            n_sample = min(n_sample, len(c_local_idx))
            selected_local.extend(rng.choice(c_local_idx, n_sample, replace=False))
        selected_local = sorted(selected_local)[:n_cells]
        ct_global_idx = ct_global_idx[selected_local]
        log.info("  Subsampled to %d cells (stratified by condition)",
                 len(ct_global_idx))

    adata_sub = adata[ct_global_idx].copy()
    log.info(
        "  Condition distribution: %s",
        dict(adata_sub.obs["condition"].astype(str).value_counts()),
    )
    log.info(
        "  Per-label cell counts: %s",
        dict(adata_sub.obs["cell_type"].astype(str).value_counts()),
    )
    return adata_sub, ct_global_idx


def compute_gene_activity(anndataset_path: str, global_indices: np.ndarray) -> ad.AnnData:
    """Compute gene activity from the AnnDataSet and subset by positional
    indices. Mirrors `04b_scenic_grn_from_activity.compute_gene_activity` but
    relabels logging from 'hepatocyte' to generic.
    """
    import snapatac2 as snap

    log.info("Computing gene activity from AnnDataSet: %s", anndataset_path)
    t0 = time.time()
    dataset = snap.read_dataset(anndataset_path)
    log.info("  AnnDataSet: %d cells x %d features", dataset.n_obs, dataset.n_vars)
    gene_mat = snap.pp.make_gene_matrix(dataset, gene_anno=snap.genome.hg38)
    log.info("  Gene activity matrix: %d cells x %d genes",
             gene_mat.n_obs, gene_mat.n_vars)
    if hasattr(gene_mat, "to_memory"):
        gene_mat = gene_mat.to_memory()

    log.info("  Subsetting gene activity to %d cells by position",
             len(global_indices))
    gene_mat_sub = gene_mat[global_indices].copy()
    gene_mat_sub.obs_names = [f"cell_{i}" for i in range(gene_mat_sub.n_obs)]

    del gene_mat, dataset
    gc.collect()

    log.info("  Gene activity subset: %d cells x %d genes in %.1fs",
             gene_mat_sub.n_obs, gene_mat_sub.n_vars, time.time() - t0)
    return gene_mat_sub


# ---------------------------------------------------------------------------
# Output writer (cell-type-prefixed schema, mirrors Script 04b columns)
# ---------------------------------------------------------------------------
def export_results(
    out_dir: str,
    celltype_name: str,
    regulon_df: pd.DataFrame,
    enhancer_df: pd.DataFrame,
    activity_df: pd.DataFrame,
    disease_df: pd.DataFrame,
):
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    tag = celltype_name.lower()

    path = out / f"{tag}_regulons.csv"
    regulon_df.to_csv(path, index=False)
    log.info("  Saved: %s (%d rows)", path, len(regulon_df))

    path = out / f"{tag}_enhancer_gene_links.csv"
    if enhancer_df is not None and not enhancer_df.empty:
        enhancer_df.to_csv(path, index=False)
    else:
        pd.DataFrame(columns=[
            "enhancer_chr", "enhancer_start", "enhancer_end",
            "target_gene", "tf_name", "correlation_rna_atac",
            "tf_enhancer_corr", "distance_to_tss",
        ]).to_csv(path, index=False)
    log.info("  Saved: %s (%d links)", path,
             len(enhancer_df) if enhancer_df is not None else 0)

    path = out / f"{tag}_regulon_activity_scores.csv"
    activity_df.to_csv(path)
    log.info("  Saved: %s (%s)", path, activity_df.shape)

    path = out / f"{tag}_disease_regulons.csv"
    disease_df.to_csv(path, index=False)
    log.info("  Saved: %s (%d disease regulons)", path, len(disease_df))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description=(
            "Per-cell-type gene-activity-based GRN (generalizes Script 04b "
            "from hepatocytes to any merged cell-type subset)."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--atac-h5ad",
        default=(
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
            "Analysis/ATAC/Human_Multiome/results/snapatac2/snapatac2_processed.h5ad"
        ),
        help="SnapATAC2 processed h5ad with cell_type + condition obs columns",
    )
    parser.add_argument(
        "--anndataset",
        default=(
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
            "Analysis/ATAC/Human_Multiome/results/snapatac2/combined.h5ads"
        ),
        help="AnnDataSet path for gene activity computation",
    )
    parser.add_argument(
        "--gtf",
        default=(
            "/gpfs/commons/home/jameslee/reference_genome/gencode_v49/"
            "gencode.v49.chr_patch_hapl_scaff.annotation.gtf.gz"
        ),
        help="GENCODE GTF (gzipped) for TSS annotation",
    )
    parser.add_argument(
        "--target-labels", required=True,
        help=(
            "Comma-separated list of cell_type labels to MERGE into the "
            "analysis subset, e.g. 'Macrophage,Kupffer_Cell'"
        ),
    )
    parser.add_argument(
        "--celltype-name", required=True,
        help="Name used for output filenames, e.g. 'Macrophage'",
    )
    parser.add_argument(
        "--output-dir",
        default=(
            "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
            "Analysis/ATAC/Human_Multiome/scenic_plus"
        ),
        help="Output directory for the 4 CSV files",
    )
    parser.add_argument(
        "--n-cells", type=int, default=15000,
        help="Max cells to subsample (balanced across conditions)",
    )
    parser.add_argument(
        "--window", type=int, default=500000,
        help="Peak-gene proximity window in bp",
    )
    parser.add_argument(
        "--min-corr", type=float, default=0.1,
        help="Min |Spearman rho| for significant correlations",
    )
    parser.add_argument(
        "--top-peaks", type=int, default=50000,
        help="Top accessible peaks to use for TF-enhancer correlations",
    )
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    target_labels = [lbl.strip() for lbl in args.target_labels.split(",") if lbl.strip()]
    if not target_labels:
        log.error("--target-labels must be a non-empty comma-separated list")
        sys.exit(1)

    log.info("=" * 70)
    log.info("Module 3 (per-CT): GRN for %s (labels=%s)",
             args.celltype_name, target_labels)
    log.info("=" * 70)
    t_start = time.time()

    # -----------------------------------------------------------------------
    # Step 1: Cell-type subset
    # -----------------------------------------------------------------------
    atac_adata, ct_global_idx = load_and_get_celltype_indices(
        args.atac_h5ad,
        target_labels=target_labels,
        celltype_name=args.celltype_name,
        n_cells=args.n_cells,
        seed=args.seed,
    )

    # -----------------------------------------------------------------------
    # Step 2: Gene activity matrix from AnnDataSet
    # -----------------------------------------------------------------------
    activity_adata = compute_gene_activity(args.anndataset, ct_global_idx)

    assert atac_adata.n_obs == activity_adata.n_obs, (
        f"Cell count mismatch: ATAC={atac_adata.n_obs}, "
        f"activity={activity_adata.n_obs}"
    )
    activity_adata.obs_names = atac_adata.obs_names.copy()
    activity_adata.obs["condition"] = atac_adata.obs["condition"].values
    # Transfer donor_id so score_regulon_activity's differential aggregates per
    # DONOR (n=18), not per cell. Without this it silently falls back to grouping
    # by condition -> pseudoreplication. (Fix 2026-05-30; mirrors hepatocyte main.)
    if "donor_id" in atac_adata.obs:
        activity_adata.obs["donor_id"] = atac_adata.obs["donor_id"].values
    else:
        log.warning("atac_adata.obs lacks 'donor_id'; donor-level test will fall "
                    "back to condition grouping (pseudoreplication risk).")
    log.info("Aligned %d cells between ATAC and gene activity", atac_adata.n_obs)

    # -----------------------------------------------------------------------
    # Step 3: TSS from GTF + TFs present in activity matrix
    # -----------------------------------------------------------------------
    gene_tss = parse_tss_from_gtf(args.gtf)

    activity_genes = set(activity_adata.var_names)
    tfs_present = sorted(KNOWN_TFS & activity_genes)
    tf_to_idx = {tf: list(activity_adata.var_names).index(tf) for tf in tfs_present}
    log.info("TFs present in gene activity: %d / %d",
             len(tfs_present), len(KNOWN_TFS))
    log.info("  TFs: %s", ", ".join(tfs_present))
    if not tfs_present:
        log.error("No known TFs found in gene activity matrix. Exiting.")
        sys.exit(1)

    # -----------------------------------------------------------------------
    # Step 4: Select top accessible peaks and densify
    # -----------------------------------------------------------------------
    log.info("Selecting top %d accessible peaks from sparse ATAC matrix...",
             args.top_peaks)
    t_rank = time.time()
    atac_sparse = atac_adata.X
    if not sp.issparse(atac_sparse):
        atac_sparse = sp.csr_matrix(atac_sparse)
    else:
        atac_sparse = sp.csr_matrix(atac_sparse)

    peak_accessibility = np.asarray(atac_sparse.mean(axis=0)).ravel()
    top_peak_idx = np.argsort(peak_accessibility)[-args.top_peaks:]
    top_peak_idx = np.sort(top_peak_idx)

    all_peak_names = list(atac_adata.var_names)
    top_peak_names = [all_peak_names[i] for i in top_peak_idx]
    log.info(
        "  Selected %d peaks (mean accessibility range: %.4f - %.4f)",
        len(top_peak_idx),
        peak_accessibility[top_peak_idx[0]],
        peak_accessibility[top_peak_idx[-1]],
    )

    atac_X = np.asarray(atac_sparse[:, top_peak_idx].toarray(), dtype=np.float32)
    log.info("  Dense ATAC subset: %d x %d (%.1f GiB)",
             *atac_X.shape, atac_X.nbytes / 1e9)
    del atac_sparse
    gc.collect()

    act_X = activity_adata.X
    if sp.issparse(act_X):
        act_X = act_X.toarray()
    act_X = np.asarray(act_X, dtype=np.float32)

    log.info("  Ranking ATAC matrix (%d x %d)...", *atac_X.shape)
    atac_X_rank = _rank_matrix(atac_X)
    log.info("  Ranking activity matrix (%d x %d)...", *act_X.shape)
    activity_X_rank = _rank_matrix(act_X)

    del atac_X, act_X
    gc.collect()
    log.info("  Rank matrices computed in %.1fs", time.time() - t_rank)

    # -----------------------------------------------------------------------
    # Step 5: Peak-gene proximity + TF-enhancer + enhancer-gene correlations
    # -----------------------------------------------------------------------
    peak_df = parse_peak_coords(top_peak_names)
    proximity = build_peak_gene_proximity(
        peak_df, gene_tss,
        gene_names=list(activity_adata.var_names),
        window=args.window,
    )

    tf_enhancer_df = compute_tf_enhancer_correlations(
        atac_X_rank=atac_X_rank,
        activity_X_rank=activity_X_rank,
        tf_indices=tf_to_idx,
        peak_names=top_peak_names,
        top_n_peaks=len(top_peak_names),
        peak_accessibility=None,
        min_corr=args.min_corr,
    )

    enhancer_gene_df = compute_enhancer_gene_correlations(
        atac_X_rank=atac_X_rank,
        activity_X_rank=activity_X_rank,
        tf_enhancer_df=tf_enhancer_df,
        proximity=proximity,
        peak_df=peak_df,
        min_corr=args.min_corr,
    )

    # -----------------------------------------------------------------------
    # Step 6: Regulons + differential activity (MASLD vs Normal)
    # -----------------------------------------------------------------------
    regulons = assemble_regulons(tf_enhancer_df, enhancer_gene_df)
    if len(regulons) < 10:
        log.warning(
            "Only %d TFs with regulons for %s (threshold: 10). "
            "Consider lowering --min-corr.",
            len(regulons), args.celltype_name,
        )

    # score_regulon_activity now returns 4 frames (added the TF-own-gene-excluded
    # differential for the circularity check). (Fix 2026-05-30; was 3-tuple.)
    activity_df, disease_df, regulon_summary_df, disease_df_excl = score_regulon_activity(
        activity_adata, regulons,
    )

    # -----------------------------------------------------------------------
    # Step 7: Export
    # -----------------------------------------------------------------------
    export_results(
        args.output_dir,
        celltype_name=args.celltype_name,
        regulon_df=regulon_summary_df,
        enhancer_df=enhancer_gene_df,
        activity_df=activity_df,
        disease_df=disease_df,
    )
    # Self-excluded (TF-own-gene-removed) donor-level differential, for the
    # regulon-circularity robustness check (parallels the hepatocyte output).
    if disease_df_excl is not None and len(disease_df_excl):
        excl_path = Path(args.output_dir) / f"{args.celltype_name}_disease_regulons_excl_self.csv"
        disease_df_excl.to_csv(excl_path, index=False)
        log.info("  Wrote %s (%d rows)", excl_path, len(disease_df_excl))

    elapsed = time.time() - t_start
    log.info("=" * 70)
    log.info("Per-CT GRN COMPLETE for %s in %.1f minutes",
             args.celltype_name, elapsed / 60)
    log.info("  Cells analyzed:        %d", atac_adata.n_obs)
    log.info("  TFs analyzed:          %d", len(tfs_present))
    log.info("  TFs with regulons:     %d", len(regulons))
    log.info("  TF-enhancer links:     %d", len(tf_enhancer_df))
    log.info("  Enhancer-gene links:   %d", len(enhancer_gene_df))
    log.info(
        "  Unique target genes:   %d",
        enhancer_gene_df["target_gene"].nunique()
        if not enhancer_gene_df.empty else 0,
    )
    log.info("  Disease regulons:      %d (padj<0.05)", len(disease_df))
    log.info("  Output directory:      %s", args.output_dir)
    log.info("=" * 70)


if __name__ == "__main__":
    main()
