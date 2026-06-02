#!/usr/bin/env python3
"""
16e_commot_regulated_genes.py — Identify signaling-regulated genes via COMMOT.

Uses COMMOT's communication_impact (gradient boosted trees) to identify genes
whose expression is explained by incoming signaling per pathway. Restricts
analysis to dream DEGs for computational efficiency.

SLURM: --partition=cpu --cpus=16 --mem=128G --time=8:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_dream_degs,
    save_csv, print_header, print_step,
)
from spatial_stats import get_commot_pathways, match_key_pathways


def load_commot_adata(condition):
    """Load COMMOT results h5ad for a condition."""
    cond_label = condition.replace(" ", "_").replace("/", "_")
    path = RESULTS_DIR / "commot" / f"adata_commot_{cond_label}.h5ad"
    if not path.exists():
        print(f"  WARNING: {path} not found")
        return None
    adata = sc.read_h5ad(path)
    print(f"  Loaded {condition}: {adata.n_obs} spots")
    return adata


def get_available_pathways(adata):
    """Extract pathway names available in COMMOT results.

    Fixed (F080 + coverage gap): real pathway-level matrices are obsp keys
    ``commot-CellChat-<PATHWAY>`` (no 2nd '-'). The old obsm/'sum' parse only
    matched ``commot-CellChat-sum-sender`` and returned ['sum'], which then got
    passed to communication_impact as pathway_name='sum' (meaningless, and a
    contributor to the to_adata failure). Delegate to canonical discovery.
    """
    return get_commot_pathways(adata, db="CellChat")


def get_deg_genes_in_adata(adata, dream_degs):
    """Get intersection of dream DEG symbols with adata var_names."""
    gene_col = "symbol" if "symbol" in dream_degs.columns else dream_degs.columns[0]
    deg_symbols = set(dream_degs[gene_col].dropna())
    overlap = sorted(deg_symbols & set(adata.var_names))
    print(f"    Dream DEGs in adata: {len(overlap)}/{len(deg_symbols)}")
    return overlap


def run_communication_impact(adata, pathway, target_genes, condition):
    """Run COMMOT communication_impact for a pathway."""
    import commot as ct

    try:
        ct.tl.communication_impact(
            adata,
            database_name="CellChat",
            pathway_name=pathway,
            ds_genes=target_genes,
        )
    except Exception as e:
        print(f"      WARNING: communication_impact failed for "
              f"{pathway}/{condition}: {e}")
        return None

    # Extract impact scores from adata.var or uns
    # COMMOT stores impact results in adata.var with columns like
    # 'commot_impact-CellChat-{pathway}-{gene}'
    impact_cols = [c for c in adata.var.columns
                   if f"commot_impact" in c and pathway in c]

    if not impact_cols:
        # Try alternative storage in uns
        uns_key = f"commot_impact-CellChat-{pathway}"
        if uns_key in adata.uns:
            impact_data = adata.uns[uns_key]
            if isinstance(impact_data, pd.DataFrame):
                return impact_data
            elif isinstance(impact_data, dict):
                return pd.DataFrame([impact_data])

        print(f"      No impact results found for {pathway}")
        return None

    # Build impact summary
    records = []
    for col in impact_cols:
        # Parse gene name from column
        gene = col.replace(f"commot_impact-CellChat-{pathway}-", "")
        score = adata.var[col].values
        if hasattr(score, "mean"):
            records.append({
                "gene": gene,
                "pathway": pathway,
                "condition": condition,
                "impact_score": float(np.nanmean(score)),
            })

    return pd.DataFrame(records) if records else None


def main():
    print_header("16e: COMMOT Signaling-Regulated Genes")

    config = load_config()
    commot_config = config["commot"]
    key_pathways = commot_config.get("key_pathways", [])
    output_dir = RESULTS_DIR / "commot"
    output_dir.mkdir(parents=True, exist_ok=True)

    # ── Load dream DEGs for target gene restriction ──
    print_step("Loading dream DEGs", 1, 4)
    dream_degs = load_dream_degs(padj_thresh=0.1, lfc_thresh=0.0)

    # ── Process each condition ──
    print_step("Loading COMMOT results", 2, 4)
    conditions = ["Healthy", "Steatotic"]
    all_impact = []

    for cond_idx, condition in enumerate(conditions):
        print(f"\n  Condition: {condition}")
        adata = load_commot_adata(condition)
        if adata is None:
            continue

        available = get_available_pathways(adata)
        # Case-insensitive + WNT-family match (F080/F082); fall back to ALL real
        # pathways (never the obsm 'sum' aggregate) if config names don't match.
        pathways_to_process = match_key_pathways(key_pathways, available)
        if not pathways_to_process:
            if not available:
                print(f"    ERROR: no real CellChat pathways found in {condition} "
                      f"obsp — COMMOT result is malformed.")
                continue
            pathways_to_process = available
            print(f"    No config key_pathways matched (config={key_pathways}, "
                  f"available={available}); using all {len(pathways_to_process)} "
                  f"available pathways")

        # Get DEG genes present in this adata
        target_genes = get_deg_genes_in_adata(adata, dream_degs)
        if not target_genes:
            print(f"    WARNING: No DEG genes found in adata, using top variable genes")
            target_genes = adata.var_names[:500].tolist()

        # ── Run impact per pathway ──
        print_step(f"Running communication_impact ({condition})", 3, 4)
        condition_results = []

        for i, pathway in enumerate(pathways_to_process):
            print(f"    [{i+1}/{len(pathways_to_process)}] {pathway}...")
            impact_df = run_communication_impact(
                adata, pathway, target_genes, condition,
            )
            if impact_df is not None and len(impact_df) > 0:
                condition_results.append(impact_df)
                print(f"      Genes with impact scores: {len(impact_df)}")

        if condition_results:
            cond_df = pd.concat(condition_results, ignore_index=True)
            cond_label = condition.replace(" ", "_").replace("/", "_")
            save_csv(cond_df, f"signaling_regulated_genes_{cond_label}.csv",
                     subdir="commot")
            all_impact.append(cond_df)

            # Report top regulated genes per pathway
            for pathway in cond_df["pathway"].unique():
                pw_df = cond_df[cond_df["pathway"] == pathway]
                pw_df = pw_df.sort_values("impact_score", ascending=False)
                top3 = pw_df.head(3)
                if len(top3) > 0:
                    top_str = ", ".join(
                        f"{r['gene']}({r['impact_score']:.3f})"
                        for _, r in top3.iterrows()
                    )
                    print(f"      {pathway} top: {top_str}")

    # ── Merge across conditions ──
    print_step("Merging results across conditions", 4, 4)
    if all_impact:
        combined = pd.concat(all_impact, ignore_index=True)

        # Identify top signaling-regulated genes (highest impact across pathways)
        gene_summary = combined.groupby("gene").agg(
            mean_impact=("impact_score", "mean"),
            max_impact=("impact_score", "max"),
            n_pathways=("pathway", "nunique"),
            n_conditions=("condition", "nunique"),
            top_pathway=("impact_score", lambda x: combined.loc[x.idxmax(), "pathway"]
                         if len(x) > 0 else ""),
        ).reset_index()
        gene_summary = gene_summary.sort_values("mean_impact", ascending=False)
        save_csv(gene_summary, "signaling_regulated_genes_summary.csv",
                 subdir="commot")

        print(f"\n  Summary:")
        print(f"    Total gene-pathway-condition entries: {len(combined)}")
        print(f"    Unique regulated genes: {combined['gene'].nunique()}")
        print(f"    Pathways with results: {combined['pathway'].nunique()}")
        print(f"\n  Top 10 signaling-regulated genes:")
        for _, row in gene_summary.head(10).iterrows():
            print(f"    {row['gene']}: impact={row['mean_impact']:.3f}, "
                  f"pathways={int(row['n_pathways'])}, "
                  f"top={row['top_pathway']}")
    else:
        # Guard (coverage-gap P0): communication_impact produced ZERO rows for
        # every pathway/condition. Previously 16e silently 'completed' with no
        # CSV written, so 16g/16i enriched against an EMPTY gene set. Fail loudly
        # instead so the pipeline halts rather than emitting void downstream
        # numbers. (No signaling_regulated_genes_summary.csv is written here.)
        print("\n  ERROR: communication_impact returned NO signaling-regulated "
              "genes for any pathway/condition. Not writing an empty summary; "
              "downstream 16g/16i must not run against an empty gene set. "
              "Check the COMMOT communication_impact API/version (the "
              "'NoneType has no attribute to_adata' error in the run log) and "
              "that pathway discovery returns real CellChat pathways.")
        print_header("16e: FAILED (empty impact set)")
        sys.exit(1)

    print_header("16e: Complete")


if __name__ == "__main__":
    main()
