#!/usr/bin/env python3
"""
05_annotation_sensitivity.py — Sensitivity comparison between gene-activity-based
and label-transfer-based scATAC cell type annotations.

Evaluates impact on downstream analyses:
  1. Cell-level concordance rate and per-cell-type precision/recall
  2. chromVAR TF activity: Spearman rho between old and new
  3. SCENIC+ regulon overlap: recovery of disease regulons
  4. Atlas S5 column changes

Pipeline position:
  02b (label transfer) → 02c (peaks) → 03 (chromVAR) → 04 (SCENIC+)
                                                        → **05_annotation_sensitivity.py**

Inputs:
  - results/label_transfer/snapatac2_label_transferred.h5ad
  - results/chromvar/chromvar_tf_activity.csv (original)
  - results/chromvar_v2/chromvar_tf_activity.csv (re-run)
  - results/scenic_plus/hepatocyte_regulons.csv (original)
  - results/scenic_plus_v2/hepatocyte_regulons.csv (re-run)

Outputs (in results/label_transfer/sensitivity/):
  - sensitivity_report.txt       Full comparison report
  - umap_old_vs_new.pdf          UMAP colored by old vs new labels
  - tf_activity_scatter.pdf      chromVAR TF activity scatter (old vs new)
  - regulon_overlap.csv          SCENIC+ regulon recovery

Environment:
  micromamba activate snapatac2

Usage:
  cd Analysis/ATAC/Human_Multiome
  python scripts/05_annotation_sensitivity.py \
      --label_transferred results/label_transfer/snapatac2_label_transferred.h5ad \
      --chromvar_old results/chromvar/chromvar_tf_activity.csv \
      --chromvar_new results/chromvar_v2/chromvar_tf_activity.csv \
      --scenic_old results/scenic_plus/hepatocyte_regulons.csv \
      --scenic_new results/scenic_plus_v2/hepatocyte_regulons.csv \
      --output_dir results/label_transfer/sensitivity
"""

from __future__ import annotations

import argparse
import gc
import json
import logging
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

warnings.filterwarnings("ignore", category=FutureWarning)
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger(__name__)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import seaborn as sns  # noqa: E402


def parse_args():
    p = argparse.ArgumentParser(
        description="Sensitivity analysis: gene-activity vs label-transfer annotations"
    )
    p.add_argument(
        "--label_transferred",
        default="results/label_transfer/snapatac2_label_transferred.h5ad",
    )
    p.add_argument(
        "--chromvar_old", default="results/chromvar/chromvar_tf_activity.csv"
    )
    p.add_argument(
        "--chromvar_new", default="results/chromvar_v2/chromvar_tf_activity.csv"
    )
    p.add_argument(
        "--scenic_old", default="results/scenic_plus/hepatocyte_regulons.csv"
    )
    p.add_argument(
        "--scenic_new", default="results/scenic_plus_v2/hepatocyte_regulons.csv"
    )
    p.add_argument(
        "--output_dir", default="results/label_transfer/sensitivity"
    )
    return p.parse_args()


def harmonize_labels(labels: pd.Series) -> pd.Series:
    """Map cell type labels to a shared major-type vocabulary.

    Gene-activity labels (e.g., "Hepatocyte", "Stellate_Cell") and CellTypist
    labels (e.g., "Hepatocytes", "Fibroblasts") use different naming. This
    maps both to a common set for meaningful concordance comparison.
    """
    mapping = {
        # Gene-activity -> major type
        "Hepatocyte": "Hepatocyte",
        "Cholangiocyte": "Cholangiocyte",
        "Stellate_Cell": "Mesenchymal",
        "Endothelial": "Endothelial",
        "LSEC": "Endothelial",
        "Kupffer_Cell": "Myeloid",
        "Macrophage": "Myeloid",
        "NK_T_Cell": "Lymphoid",
        "B_Cell": "Lymphoid",
        "Plasma_Cell": "Plasma_cell",
        # CellTypist -> major type
        "Hepatocytes": "Hepatocyte",
        "Cholangiocytes": "Cholangiocyte",
        "Fibroblasts": "Mesenchymal",
        "Endothelial cells": "Endothelial",
        "Endothelial_cells": "Endothelial",
        "Macrophages": "Myeloid",
        "Mono+mono derived cells": "Myeloid",
        "Mono+mono_derived_cells": "Myeloid",
        "Kupffer cells": "Myeloid",
        "T cells": "Lymphoid",
        "T_cells": "Lymphoid",
        "Circulating NK/NKT": "Lymphoid",
        "Circulating_NK/NKT": "Lymphoid",
        "Circulating_NK_NKT": "Lymphoid",
        "Resident NK": "Lymphoid",
        "Resident_NK": "Lymphoid",
        "B cells": "Lymphoid",
        "B_cells": "Lymphoid",
        "Plasma cells": "Plasma_cell",
        "Plasma_cells": "Plasma_cell",
        "cDC1s": "Myeloid",
        "cDC2s": "Myeloid",
        "pDCs": "Myeloid",
        "Neutrophils": "Myeloid",
        "Basophils": "Myeloid",
    }
    return labels.map(lambda x: mapping.get(x, x))


def annotation_concordance(adata, output_dir: Path) -> dict:
    """Compute cell-level concordance, per-cell-type precision/recall.

    Uses both raw labels (for full confusion matrix) and harmonized major-type
    labels (for meaningful concordance across different naming conventions).
    """
    old = adata.obs["cell_type_gene_activity"].astype(str)
    new = adata.obs["cell_type_transferred"].astype(str)

    # Exclude low-confidence / unassigned for metrics
    valid = ~new.isin(["Low_confidence", "Unassigned"])
    old_v = old[valid]
    new_v = new[valid]

    concordance = (old_v == new_v).mean()
    n_changed = (old_v != new_v).sum()

    # Per-cell-type overlap metrics (no ground truth assumed — report symmetric Jaccard)
    all_types = sorted(set(old_v.unique()) | set(new_v.unique()))
    ct_metrics = []
    for ct in all_types:
        tp = ((old_v == ct) & (new_v == ct)).sum()
        fp = ((old_v != ct) & (new_v == ct)).sum()  # in new but not old
        fn = ((old_v == ct) & (new_v != ct)).sum()  # in old but not new
        precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
        recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
        f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0.0
        old_count = (old_v == ct).sum()
        new_count = (new_v == ct).sum()

        # Jaccard
        old_set = set(adata.obs_names[valid][old_v == ct])
        new_set = set(adata.obs_names[valid][new_v == ct])
        union = len(old_set | new_set)
        jaccard = len(old_set & new_set) / union if union > 0 else 0.0

        ct_metrics.append({
            "cell_type": ct,
            "old_count": int(old_count),
            "new_count": int(new_count),
            "precision": precision,
            "recall": recall,
            "f1": f1,
            "jaccard": jaccard,
        })

    ct_df = pd.DataFrame(ct_metrics)
    ct_df.to_csv(output_dir / "per_celltype_metrics.csv", index=False)

    # --- Harmonized major-type concordance ---
    old_harm = harmonize_labels(old_v)
    new_harm = harmonize_labels(new_v)
    concordance_harmonized = (old_harm == new_harm).mean()
    n_changed_harmonized = (old_harm != new_harm).sum()
    log.info(
        "Harmonized concordance (major types): %.1f%% (%d changed)",
        100 * concordance_harmonized, n_changed_harmonized,
    )

    # Per-major-type Jaccard
    major_types = sorted(set(old_harm.unique()) | set(new_harm.unique()))
    major_metrics = []
    for mt in major_types:
        old_set = set(adata.obs_names[valid][old_harm == mt])
        new_set = set(adata.obs_names[valid][new_harm == mt])
        union = len(old_set | new_set)
        jaccard = len(old_set & new_set) / union if union > 0 else 0.0
        major_metrics.append({
            "major_type": mt,
            "old_count": int((old_harm == mt).sum()),
            "new_count": int((new_harm == mt).sum()),
            "jaccard": jaccard,
        })

    major_df = pd.DataFrame(major_metrics)
    major_df.to_csv(output_dir / "per_majortype_metrics.csv", index=False)

    return {
        "concordance_raw": float(concordance),
        "concordance_harmonized": float(concordance_harmonized),
        "n_valid": int(valid.sum()),
        "n_changed_raw": int(n_changed),
        "n_changed_harmonized": int(n_changed_harmonized),
        "n_low_confidence": int((~valid).sum()),
        "ct_metrics": ct_df,
        "major_metrics": major_df,
    }


def umap_comparison(adata, output_dir: Path):
    """Side-by-side UMAP colored by old vs new labels."""
    if "X_umap" not in adata.obsm:
        log.warning("No UMAP embedding found, skipping UMAP comparison plot")
        return

    fig, axes = plt.subplots(1, 2, figsize=(24, 10))

    # Old labels
    sc.pl.umap(adata, color="cell_type_gene_activity", ax=axes[0], show=False,
               title="Gene Activity Scoring (original)", legend_loc="right margin",
               frameon=False)

    # New labels
    sc.pl.umap(adata, color="cell_type_transferred", ax=axes[1], show=False,
               title="Label Transfer (snRNA reference)", legend_loc="right margin",
               frameon=False)

    plt.tight_layout()
    fig.savefig(output_dir / "umap_old_vs_new.pdf", dpi=150, bbox_inches="tight")
    plt.close(fig)
    log.info("Saved UMAP comparison: %s", output_dir / "umap_old_vs_new.pdf")


def chromvar_comparison(old_path: str, new_path: str, output_dir: Path) -> dict:
    """Compare chromVAR TF activity between old and new cell type annotations."""
    old_p = Path(old_path)
    new_p = Path(new_path)

    if not old_p.exists():
        log.warning("Old chromVAR results not found: %s", old_p)
        return {"status": "old_missing"}
    if not new_p.exists():
        log.warning("New chromVAR results not found: %s", new_p)
        return {"status": "new_missing"}

    old_df = pd.read_csv(old_p)
    new_df = pd.read_csv(new_p)

    log.info("chromVAR old: %d rows, new: %d rows", len(old_df), len(new_df))

    # Harmonize cell type names in both dataframes for merge
    if "cell_type" in old_df.columns:
        old_df["cell_type_harmonized"] = harmonize_labels(old_df["cell_type"])
    if "cell_type" in new_df.columns:
        new_df["cell_type_harmonized"] = harmonize_labels(new_df["cell_type"])

    # Find common TFs and cell types for comparison
    merge_cols = []
    for col in ["tf_name", "tf", "motif"]:
        if col in old_df.columns and col in new_df.columns:
            merge_cols.append(col)
            break

    # Use harmonized cell type for merge (original names differ between runs)
    ct_merge = "cell_type_harmonized" if "cell_type_harmonized" in old_df.columns else "cell_type"
    if ct_merge in old_df.columns and ct_merge in new_df.columns:
        merge_cols.append(ct_merge)

    if len(merge_cols) < 2:
        log.warning("Cannot merge chromVAR results: missing key columns")
        return {"status": "incompatible_columns"}

    merged = old_df.merge(new_df, on=merge_cols, suffixes=("_old", "_new"))
    log.info("Merged rows: %d", len(merged))

    # Compare mean deviation scores
    dev_col_old = None
    dev_col_new = None
    for prefix in ["mean_deviation", "logFC_deviation", "mean_diff"]:
        if f"{prefix}_old" in merged.columns:
            dev_col_old = f"{prefix}_old"
            dev_col_new = f"{prefix}_new"
            break

    results = {"status": "compared", "n_merged": len(merged)}

    if dev_col_old and dev_col_new:
        valid = merged[[dev_col_old, dev_col_new]].dropna()
        if len(valid) > 2:
            rho, pval = spearmanr(valid[dev_col_old], valid[dev_col_new])
            results["spearman_rho"] = float(rho)
            results["spearman_pval"] = float(pval)

            # Direction concordance
            same_dir = (np.sign(valid[dev_col_old]) == np.sign(valid[dev_col_new])).mean()
            results["direction_concordance"] = float(same_dir)

            log.info("chromVAR Spearman rho: %.3f (p=%.2e)", rho, pval)
            log.info("Direction concordance: %.1f%%", 100 * same_dir)

            # Scatter plot
            fig, ax = plt.subplots(figsize=(8, 8))
            ax.scatter(valid[dev_col_old], valid[dev_col_new], alpha=0.3, s=10)
            lim = max(abs(valid[dev_col_old]).max(), abs(valid[dev_col_new]).max()) * 1.1
            ax.plot([-lim, lim], [-lim, lim], "r--", alpha=0.5)
            ax.set_xlabel(f"Old annotation ({dev_col_old.replace('_old', '')})")
            ax.set_ylabel(f"New annotation ({dev_col_new.replace('_new', '')})")
            ax.set_title(f"chromVAR TF Activity: rho={rho:.3f}, dir={same_dir:.1%}")
            ax.set_xlim(-lim, lim)
            ax.set_ylim(-lim, lim)
            ax.set_aspect("equal")
            plt.tight_layout()
            fig.savefig(output_dir / "tf_activity_scatter.pdf", dpi=150, bbox_inches="tight")
            plt.close(fig)

    return results


def scenic_comparison(old_path: str, new_path: str, output_dir: Path) -> dict:
    """Compare SCENIC+ regulons between old and new annotations."""
    old_p = Path(old_path)
    new_p = Path(new_path)

    # Guard against self-comparison (review 2026-05-30, F019): if both paths
    # resolve to the same file the test trivially reports 100% TF/regulon
    # recovery, which is a tautology, not evidence of robustness.
    if old_p.resolve() == new_p.resolve():
        raise ValueError(
            f"scenic_comparison: --scenic_old and --scenic_new resolve to the same "
            f"file ({old_p}); a self-comparison is meaningless (guaranteed 100% "
            f"recovery). Point them at the two distinct annotation runs (e.g. the "
            f"gene-activity backup vs the label-transfer v2 SCENIC+ regulons)."
        )

    if not old_p.exists():
        log.warning("Old SCENIC+ results not found: %s", old_p)
        return {"status": "old_missing"}
    if not new_p.exists():
        log.warning("New SCENIC+ results not found: %s", new_p)
        return {"status": "new_missing"}

    old_df = pd.read_csv(old_p)
    new_df = pd.read_csv(new_p)

    log.info("SCENIC+ old: %d regulons, new: %d regulons", len(old_df), len(new_df))

    # TF-level comparison
    tf_col = "tf_name" if "tf_name" in old_df.columns else old_df.columns[0]
    old_tfs = set(old_df[tf_col].unique())
    new_tfs = set(new_df[tf_col].unique())

    shared_tfs = old_tfs & new_tfs
    old_only = old_tfs - new_tfs
    new_only = new_tfs - old_tfs

    results = {
        "status": "compared",
        "old_n_regulons": len(old_df),
        "new_n_regulons": len(new_df),
        "old_n_tfs": len(old_tfs),
        "new_n_tfs": len(new_tfs),
        "shared_tfs": len(shared_tfs),
        "recovery_rate": len(shared_tfs) / len(old_tfs) if old_tfs else 0.0,
        "old_only_tfs": sorted(old_only),
        "new_only_tfs": sorted(new_only),
    }

    log.info(
        "SCENIC+ TF recovery: %d/%d (%.1f%%)",
        len(shared_tfs), len(old_tfs), 100 * results["recovery_rate"],
    )
    if old_only:
        log.info("  Lost TFs: %s", ", ".join(sorted(old_only)[:10]))
    if new_only:
        log.info("  New TFs: %s", ", ".join(sorted(new_only)[:10]))

    # Check key hepatocyte TFs
    key_tfs = ["HNF4A", "THRB", "NR1H4", "PPARA", "RXRA", "CEBPA", "CEBPB"]
    key_recovery = {tf: tf in new_tfs for tf in key_tfs}
    results["key_tf_recovery"] = key_recovery
    log.info("Key TF recovery: %s", key_recovery)

    # Save overlap table
    overlap_df = pd.DataFrame({
        "tf": sorted(old_tfs | new_tfs),
        "in_old": [tf in old_tfs for tf in sorted(old_tfs | new_tfs)],
        "in_new": [tf in new_tfs for tf in sorted(old_tfs | new_tfs)],
    })
    overlap_df.to_csv(output_dir / "regulon_overlap.csv", index=False)

    # Disease regulon comparison
    disease_old_path = old_p.parent / "disease_regulons.csv"
    disease_new_path = new_p.parent / "disease_regulons.csv"

    if disease_old_path.exists() and disease_new_path.exists():
        dis_old = pd.read_csv(disease_old_path)
        dis_new = pd.read_csv(disease_new_path)
        old_dis_tfs = set(dis_old[tf_col].unique()) if tf_col in dis_old.columns else set()
        new_dis_tfs = set(dis_new[tf_col].unique()) if tf_col in dis_new.columns else set()
        results["disease_regulon_recovery"] = {
            "old_count": len(old_dis_tfs),
            "new_count": len(new_dis_tfs),
            "shared": len(old_dis_tfs & new_dis_tfs),
            "recovery_rate": len(old_dis_tfs & new_dis_tfs) / len(old_dis_tfs) if old_dis_tfs else 0.0,
        }
        log.info(
            "Disease regulon recovery: %d/%d (%.1f%%)",
            len(old_dis_tfs & new_dis_tfs), len(old_dis_tfs),
            100 * results["disease_regulon_recovery"]["recovery_rate"],
        )

    return results


def main():
    args = parse_args()
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fh = logging.FileHandler(output_dir / "sensitivity.log")
    fh.setFormatter(logging.Formatter("%(asctime)s [%(levelname)s] %(message)s"))
    log.addHandler(fh)

    log.info("=" * 60)
    log.info("Annotation Sensitivity Analysis")
    log.info("=" * 60)

    # --- 1. Cell-level concordance ---
    log.info("\n=== 1. Cell-Level Concordance ===")
    adata = sc.read_h5ad(args.label_transferred)

    # Report on valid cells (excluding Low_confidence/Unassigned)
    ann_results = annotation_concordance(adata, output_dir)
    log.info("Raw concordance: %.1f%%", 100 * ann_results["concordance_raw"])
    log.info("Harmonized concordance: %.1f%%", 100 * ann_results["concordance_harmonized"])
    log.info("Changed (harmonized): %d / %d valid cells", ann_results["n_changed_harmonized"], ann_results["n_valid"])
    log.info(
        "NOTE: %d low-confidence cells (%.1f%%) excluded from concordance. "
        "Including them would lower concordance.",
        ann_results["n_low_confidence"],
        100 * ann_results["n_low_confidence"] / adata.n_obs,
    )

    # --- 2. UMAP comparison ---
    log.info("\n=== 2. UMAP Comparison ===")
    umap_comparison(adata, output_dir)

    del adata
    gc.collect()

    # --- 3. chromVAR comparison ---
    log.info("\n=== 3. chromVAR TF Activity Comparison ===")
    cv_results = chromvar_comparison(args.chromvar_old, args.chromvar_new, output_dir)

    # --- 4. SCENIC+ comparison ---
    log.info("\n=== 4. SCENIC+ Regulon Comparison ===")
    sc_results = scenic_comparison(args.scenic_old, args.scenic_new, output_dir)

    # --- Write report ---
    report_path = output_dir / "sensitivity_report.txt"
    with open(report_path, "w") as f:
        f.write("Annotation Sensitivity Analysis Report\n")
        f.write("=" * 60 + "\n\n")

        f.write("1. CELL-LEVEL CONCORDANCE\n")
        f.write(f"   Raw concordance:        {ann_results['concordance_raw']:.3f} (different naming conventions)\n")
        f.write(f"   Harmonized concordance: {ann_results['concordance_harmonized']:.3f} (major cell types)\n")
        f.write(f"   Valid cells:            {ann_results['n_valid']}\n")
        f.write(f"   Changed (harmonized):   {ann_results['n_changed_harmonized']}\n")
        f.write(f"   Low confidence:         {ann_results['n_low_confidence']}\n\n")

        f.write("   Per-major-type metrics (harmonized):\n")
        for _, row in ann_results["major_metrics"].iterrows():
            f.write(
                f"     {row['major_type']:20s}  J={row['jaccard']:.3f}  "
                f"(old={row['old_count']}, new={row['new_count']})\n"
            )
        f.write("\n   Per-cell-type metrics (raw labels):\n")
        for _, row in ann_results["ct_metrics"].iterrows():
            f.write(
                f"     {row['cell_type']:25s}  P={row['precision']:.3f}  "
                f"R={row['recall']:.3f}  F1={row['f1']:.3f}  J={row['jaccard']:.3f}  "
                f"(old={row['old_count']}, new={row['new_count']})\n"
            )

        f.write("\n2. CHROMVAR TF ACTIVITY\n")
        if cv_results.get("spearman_rho") is not None:
            f.write(f"   Spearman rho:      {cv_results['spearman_rho']:.3f}\n")
            f.write(f"   Direction conc.:   {cv_results['direction_concordance']:.3f}\n")
        else:
            f.write(f"   Status: {cv_results.get('status', 'unknown')}\n")

        f.write("\n3. SCENIC+ REGULONS\n")
        if sc_results.get("status") == "compared":
            f.write(f"   Old regulons:      {sc_results['old_n_regulons']}\n")
            f.write(f"   New regulons:      {sc_results['new_n_regulons']}\n")
            f.write(f"   TF recovery:       {sc_results['shared_tfs']}/{sc_results['old_n_tfs']} "
                    f"({100*sc_results['recovery_rate']:.1f}%)\n")
            if "key_tf_recovery" in sc_results:
                f.write("   Key TFs:\n")
                for tf, found in sc_results["key_tf_recovery"].items():
                    f.write(f"     {tf:10s}  {'RECOVERED' if found else 'LOST'}\n")
            if "disease_regulon_recovery" in sc_results:
                dr = sc_results["disease_regulon_recovery"]
                f.write(f"   Disease regulons:  {dr['shared']}/{dr['old_count']} "
                        f"({100*dr['recovery_rate']:.1f}%)\n")
        else:
            f.write(f"   Status: {sc_results.get('status', 'unknown')}\n")

    log.info("\nSensitivity report saved: %s", report_path)

    # Save results as JSON
    json_results = {
        "annotation": {
            "concordance_raw": ann_results["concordance_raw"],
            "concordance_harmonized": ann_results["concordance_harmonized"],
            "n_valid": ann_results["n_valid"],
            "n_changed_harmonized": ann_results["n_changed_harmonized"],
        },
        "chromvar": {
            k: v for k, v in cv_results.items()
            if isinstance(v, (int, float, str, bool))
        },
        "scenic": {
            k: v for k, v in sc_results.items()
            if isinstance(v, (int, float, str, bool, list, dict))
        },
    }
    with open(output_dir / "sensitivity_results.json", "w") as f:
        json.dump(json_results, f, indent=2, default=str)

    log.info("Done!")


if __name__ == "__main__":
    main()
