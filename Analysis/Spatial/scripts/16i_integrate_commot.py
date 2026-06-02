#!/usr/bin/env python3
"""
16i_integrate_commot.py — Integrate COMMOT results into multi-evidence atlas.

Adds COMMOT communication columns to the multi-evidence atlas: signaling
regulation status, top pathway, sender/receiver scores, and disease delta.

SLURM: --partition=cpu --cpus=4 --mem=32G --time=1:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_multi_evidence_atlas,
    save_csv, print_header, print_step,
)


def load_signaling_summary():
    """Load signaling-regulated gene summary from 16e."""
    path = RESULTS_DIR / "commot" / "signaling_regulated_genes_summary.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, index_col=0)


def load_differential_communication():
    """Load differential communication from 16f."""
    path = RESULTS_DIR / "commot" / "differential_communication.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, index_col=0)


def load_sender_receiver_shift():
    """Load sender/receiver cell type shift from 16f."""
    path = RESULTS_DIR / "commot" / "sender_receiver_shift.csv"
    if not path.exists():
        return pd.DataFrame()
    return pd.read_csv(path, index_col=0)


def load_per_condition_signaling():
    """Load per-condition signaling results to compute per-gene scores."""
    results = {}
    for condition in ["Healthy", "Steatotic"]:
        cond_label = condition.replace(" ", "_").replace("/", "_")
        path = RESULTS_DIR / "commot" / f"signaling_regulated_genes_{cond_label}.csv"
        if path.exists():
            results[condition] = pd.read_csv(path, index_col=0)
    return results


def build_atlas_columns(atlas, sig_summary, diff_comm, sr_shift, per_cond):
    """Build COMMOT atlas columns for each gene."""
    sym_col = "human_symbol" if "human_symbol" in atlas.columns else "symbol"

    # Initialize columns
    atlas["commot_is_signaling_regulated"] = False
    atlas["commot_top_pathway"] = np.nan
    atlas["commot_sender_score"] = np.nan
    atlas["commot_receiver_score"] = np.nan
    atlas["commot_disease_delta"] = np.nan

    # ── Signaling-regulated status and top pathway ──
    if len(sig_summary) > 0:
        gene_col = "gene" if "gene" in sig_summary.columns else sig_summary.index.name

        if gene_col and gene_col in sig_summary.columns:
            sig_map = sig_summary.set_index(gene_col)
        else:
            sig_map = sig_summary

        sig_gene_set = set(sig_map.index)

        # Map to atlas
        is_regulated = atlas[sym_col].isin(sig_gene_set)
        atlas.loc[is_regulated, "commot_is_signaling_regulated"] = True

        # Top pathway
        if "top_pathway" in sig_map.columns:
            top_pw_map = sig_map["top_pathway"].to_dict()
            atlas.loc[is_regulated, "commot_top_pathway"] = \
                atlas.loc[is_regulated, sym_col].map(top_pw_map)
        elif "pathway" in sig_map.columns:
            # If no top_pathway, use the first pathway entry
            top_pw_map = sig_map.groupby(sig_map.index)["pathway"].first().to_dict()
            atlas.loc[is_regulated, "commot_top_pathway"] = \
                atlas.loc[is_regulated, sym_col].map(top_pw_map)

        print(f"  Signaling-regulated genes mapped: {is_regulated.sum()}")

    # ── Sender/receiver scores ──
    if len(sr_shift) > 0:
        # Compute per-gene mean sender/receiver scores across pathways
        # sr_shift is per cell_type x pathway x direction x condition
        # We aggregate across all cell types for per-pathway scores

        for direction in ["sender", "receiver"]:
            dir_df = sr_shift[sr_shift["direction"] == direction] \
                if "direction" in sr_shift.columns else pd.DataFrame()

            if len(dir_df) == 0:
                continue

            # For each pathway, compute the mean score across conditions
            pw_scores = dir_df.groupby("pathway")["mean_score"].mean().to_dict()

            # Map genes to their pathway scores via sig_summary
            if len(sig_summary) > 0:
                gene_col = "gene" if "gene" in sig_summary.columns else sig_summary.index.name
                if gene_col and gene_col in sig_summary.columns:
                    for _, row in sig_summary.iterrows():
                        gene = row[gene_col] if gene_col in row.index else row.name
                        pw = row.get("top_pathway", row.get("pathway", ""))
                        if pw in pw_scores:
                            mask = atlas[sym_col] == gene
                            atlas.loc[mask, f"commot_{direction}_score"] = pw_scores[pw]

    # ── Disease delta (fold change in signaling involvement) ──
    if per_cond:
        healthy_genes = {}
        disease_genes = {}

        for condition, df in per_cond.items():
            gene_col = "gene" if "gene" in df.columns else df.index.name
            score_col = "impact_score" if "impact_score" in df.columns else \
                ("mean_impact" if "mean_impact" in df.columns else None)

            if gene_col is None or score_col is None:
                continue

            if gene_col in df.columns:
                gene_scores = df.groupby(gene_col)[score_col].mean().to_dict()
            else:
                gene_scores = df[score_col].to_dict()

            if "healthy" in condition.lower():
                healthy_genes = gene_scores
            else:
                disease_genes = gene_scores

        # Compute delta for each gene
        all_genes = set(healthy_genes.keys()) | set(disease_genes.keys())
        for gene in all_genes:
            score_h = healthy_genes.get(gene, 0.0)
            score_d = disease_genes.get(gene, 0.0)
            if score_h > 0:
                delta = score_d / score_h
            elif score_d > 0:
                delta = float("inf")
            else:
                delta = 1.0

            mask = atlas[sym_col] == gene
            if mask.any():
                atlas.loc[mask, "commot_disease_delta"] = delta

    return atlas


def main():
    print_header("16i: Integrate COMMOT into Multi-Evidence Atlas")

    config = load_config()
    output_dir = RESULTS_DIR / "commot"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Load multi-evidence atlas ──
    print_step("Loading multi-evidence atlas", 1, 4)
    atlas = load_multi_evidence_atlas()
    n_cols_before = len(atlas.columns)
    print(f"  Atlas: {len(atlas)} genes × {n_cols_before} columns")

    # ── Load COMMOT results ──
    print_step("Loading COMMOT results", 2, 4)
    sig_summary = load_signaling_summary()
    diff_comm = load_differential_communication()
    sr_shift = load_sender_receiver_shift()
    per_cond = load_per_condition_signaling()

    print(f"  Signaling summary: {len(sig_summary)} genes")
    print(f"  Differential communication: {len(diff_comm)} pathways")
    print(f"  Sender/receiver shift: {len(sr_shift)} records")
    print(f"  Per-condition results: {list(per_cond.keys())}")

    # ── Build atlas columns ──
    print_step("Building atlas columns", 3, 4)
    atlas = build_atlas_columns(atlas, sig_summary, diff_comm, sr_shift, per_cond)

    n_new = len(atlas.columns) - n_cols_before
    print(f"  Added {n_new} COMMOT columns")

    # Report coverage
    commot_cols = [c for c in atlas.columns if c.startswith("commot_")]
    for col in commot_cols:
        if atlas[col].dtype == bool:
            n_true = atlas[col].sum()
            print(f"    {col}: {n_true} genes (True)")
        else:
            n_annotated = atlas[col].notna().sum()
            print(f"    {col}: {n_annotated} genes annotated")

    # ── Save outputs ──
    print_step("Saving results", 4, 4)
    sym_col = "human_symbol" if "human_symbol" in atlas.columns else "symbol"
    # Include the gene-symbol column so the atlas-column CSV carries a join key.
    # Previously only commot_* were saved → bare numeric index, no gene id, and
    # the 17a merge silently annotated 0 genes.
    save_csv(atlas[[sym_col] + commot_cols], "commot_atlas_columns.csv", subdir="commot")

    # Summary statistics
    sym_col = "human_symbol" if "human_symbol" in atlas.columns else "symbol"
    n_regulated = atlas["commot_is_signaling_regulated"].sum()
    n_with_pathway = atlas["commot_top_pathway"].notna().sum()
    n_with_delta = atlas["commot_disease_delta"].notna().sum()

    print(f"\n  Summary:")
    print(f"    Genes in atlas: {len(atlas)}")
    print(f"    Signaling-regulated: {n_regulated} ({100*n_regulated/len(atlas):.1f}%)")
    print(f"    With top pathway: {n_with_pathway}")
    print(f"    With disease delta: {n_with_delta}")

    # Top pathway distribution
    if n_with_pathway > 0:
        pw_dist = atlas["commot_top_pathway"].value_counts().head(10)
        print(f"\n  Top pathway distribution:")
        for pw, n in pw_dist.items():
            print(f"    {pw}: {n} genes")

    print_header("16i: Complete")


if __name__ == "__main__":
    main()
