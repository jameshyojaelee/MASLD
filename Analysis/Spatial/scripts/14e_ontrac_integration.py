#!/usr/bin/env python3
"""
14e_ontrac_integration.py — Integrate the DPT niche-trajectory results
(ONTraC fallback; ONTraC CLI unavailable) into the multi-evidence atlas.

PROVENANCE NOTE (F090): the "niche trajectory" (NT) score is a DPT
diffusion-pseudotime (rooted at the highest-hepatocyte spot, computed in 14b),
NOT a native ONTraC trajectory — the ONTraC CLI is unavailable, so 14b ran its
DPT fallback (parameter_sweep.csv method=diffusion_pseudotime). The atlas COLUMN
names (spatial_nt_score_mean, spatial_niche_cluster, spatial_is_tag,
spatial_tag_direction) and the ``ontrac/`` result paths are kept verbatim for
downstream/17a join compatibility, but the underlying axis is DPT pseudotime —
see the 14c provenance note.

Maps TAGs to atlas genes, adds spatial niche trajectory columns, and runs
enrichment tests (TAGs vs dream DEGs, Conserved, drug targets, GWAS loci).

F092: the Fisher universe is the genes that were ELIGIBLE to be TAGs — i.e. the
spatially-tested gene set (var_names of the deconvolved adata, ~10,738),
intersected with atlas symbols — NOT the full ~14k atlas. Genes never measured
spatially can never be TAGs; padding the universe with them inflates odds ratios.

SLURM: --partition=cpu --cpus=8 --mem=64G --time=2:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    RESULTS_DIR, load_config, load_multi_evidence_atlas,
    load_deconvolved_adata, load_dream_degs, load_conserved,
    load_drug_targets, load_coloc_results, fisher_test,
    save_csv, print_header, print_step,
)


def load_nt_scores():
    """Load niche trajectory scores from 14b."""
    path = RESULTS_DIR / "ontrac" / "niche_trajectory_scores.csv"
    if not path.exists():
        print("  ERROR: Run 14b_ontrac_run.py first")
        sys.exit(1)
    return pd.read_csv(path, index_col=0)


def load_tags():
    """Load trajectory-associated genes from 14c."""
    path = RESULTS_DIR / "ontrac" / "trajectory_associated_genes.csv"
    if not path.exists():
        print("  WARNING: No TAG results from 14c")
        return pd.DataFrame()
    return pd.read_csv(path, index_col=0)


def compute_per_gene_nt_stats(adata, nt_df):
    """Compute mean NT score in spots where each gene is expressed.

    Also determines dominant niche cluster per gene.
    Returns DataFrame indexed by gene symbol.
    """
    shared = adata.obs_names.intersection(nt_df.index)
    adata_sub = adata[shared].copy()
    nt_scores = nt_df.loc[shared, "nt_score"].values
    niche_clusters = nt_df.loc[shared, "niche_cluster"].values

    records = []
    n_genes = adata_sub.n_vars

    print(f"  Computing per-gene NT stats for {n_genes} genes...")
    for i in range(n_genes):
        if issparse(adata_sub.X):
            expr = np.asarray(adata_sub.X[:, i].todense()).flatten()
        else:
            expr = adata_sub.X[:, i].flatten()

        # Spots where gene is expressed (>0)
        expressed = expr > 0
        if expressed.sum() < 10:
            continue

        nt_expressed = nt_scores[expressed]
        nc_expressed = niche_clusters[expressed]

        # Dominant niche cluster (mode)
        unique, counts = np.unique(nc_expressed, return_counts=True)
        dominant_niche = unique[np.argmax(counts)]

        records.append({
            "gene": adata_sub.var_names[i],
            "spatial_nt_score_mean": np.nanmean(nt_expressed),
            "spatial_niche_cluster": int(dominant_niche),
            "n_spots_expressed": int(expressed.sum()),
        })

        if (i + 1) % 5000 == 0:
            print(f"    Processed {i + 1}/{n_genes} genes...")

    return pd.DataFrame(records).set_index("gene")


def add_atlas_columns(atlas, tag_df, gene_stats):
    """Add DPT niche-trajectory columns (ONTraC fallback) to the multi-evidence atlas.

    Column names (spatial_nt_score_mean, spatial_niche_cluster, spatial_is_tag,
    spatial_tag_direction) are kept verbatim for downstream/17a join
    compatibility; the values are DPT diffusion-pseudotime, not ONTraC (F090).
    """
    sym_col = "human_symbol" if "human_symbol" in atlas.columns else (
        "symbol" if "symbol" in atlas.columns else atlas.columns[0]
    )

    # NT score mean per gene
    if len(gene_stats) > 0:
        atlas["spatial_nt_score_mean"] = atlas[sym_col].map(
            gene_stats["spatial_nt_score_mean"].to_dict()
        )
        atlas["spatial_niche_cluster"] = atlas[sym_col].map(
            gene_stats["spatial_niche_cluster"].to_dict()
        )

    # TAG status and direction
    if len(tag_df) > 0:
        gene_col = "gene" if "gene" in tag_df.columns else tag_df.index.name
        if gene_col and gene_col in tag_df.columns:
            tag_map = tag_df.set_index(gene_col) if gene_col in tag_df.columns else tag_df
        else:
            tag_map = tag_df

        if "is_tag" in tag_map.columns:
            atlas["spatial_is_tag"] = atlas[sym_col].map(
                tag_map["is_tag"].to_dict()
            ).fillna(False).astype(bool)

        if "direction" in tag_map.columns:
            atlas["spatial_tag_direction"] = atlas[sym_col].map(
                tag_map["direction"].to_dict()
            )

    return atlas


def _tag_gene_series(tag_df):
    """Return (all_tested_genes, is_tag_genes) sets from the 14c TAG table.

    The TAG table has one row per spatially-tested gene; that full set is the
    Fisher universe (genes eligible to be a TAG). ``is_tag`` marks the hits.
    """
    if len(tag_df) == 0:
        return set(), set()
    gene_col = "gene" if "gene" in tag_df.columns else tag_df.index.name
    if gene_col and gene_col in tag_df.columns:
        tested = set(tag_df[gene_col].dropna())
        hits = set(tag_df[tag_df["is_tag"]][gene_col]) if "is_tag" in tag_df.columns else set()
    else:
        tested = set(tag_df.index)
        hits = set(tag_df[tag_df["is_tag"]].index) if "is_tag" in tag_df.columns else set()
    return tested, hits


def run_enrichment_tests(tag_df, atlas):
    """Run Fisher's exact enrichment tests for TAGs vs key gene sets.

    F092: the universe is the intersection of atlas symbols AND the
    spatially-tested gene set (every gene that was eligible to be a TAG, taken
    from the 14c TAG table which has one row per tested gene), NOT the full
    atlas. Every reference set is also clipped to this universe before the test.
    """
    sym_col = "human_symbol" if "human_symbol" in atlas.columns else (
        "symbol" if "symbol" in atlas.columns else atlas.columns[0]
    )
    atlas_genes = set(atlas[sym_col].dropna())

    tested_genes, tag_genes = _tag_gene_series(tag_df)

    if len(tag_genes) == 0:
        print("  No TAGs for enrichment testing")
        return pd.DataFrame()

    # F092: restrict universe to spatially-tested ∩ atlas (TAG-eligible) genes.
    if tested_genes:
        all_genes = atlas_genes & tested_genes
    else:
        # Defensive fallback: if 14c didn't carry the full tested set, the best
        # available eligibility proxy is the atlas (old behaviour) — warn.
        print("  WARNING (F092): TAG table lacks the full spatially-tested gene "
              "list; falling back to atlas universe (odds ratios may be inflated).")
        all_genes = atlas_genes

    tag_genes = tag_genes & all_genes
    print(f"  Testing {len(tag_genes)} TAGs against reference gene sets "
          f"(universe={len(all_genes)} TAG-eligible genes)")

    tests = []

    # 1. Dream DEGs
    try:
        dream = load_dream_degs(padj_thresh=0.1, lfc_thresh=0.0)
        gene_col = "symbol" if "symbol" in dream.columns else dream.columns[0]
        deg_set = set(dream[gene_col].dropna()) & all_genes
        odds, pval, n_overlap = fisher_test(tag_genes, deg_set, all_genes)
        tests.append({
            "gene_set": "Dream_DEGs",
            "n_tags": len(tag_genes), "n_reference": len(deg_set),
            "n_overlap": n_overlap, "odds_ratio": odds, "fisher_pval": pval,
        })
        sig = "*" if pval < 0.05 else "ns"
        print(f"    TAGs vs Dream DEGs: OR={odds:.2f}, p={pval:.2e}, "
              f"overlap={n_overlap} {sig}")
    except Exception as e:
        print(f"    WARNING: Dream DEGs test failed: {e}")

    # 2. Conserved
    cc = load_conserved()
    if cc:
        cc_set = set(cc) & all_genes
        odds, pval, n_overlap = fisher_test(tag_genes, cc_set, all_genes)
        tests.append({
            "gene_set": "Conserved",
            "n_tags": len(tag_genes), "n_reference": len(cc_set),
            "n_overlap": n_overlap, "odds_ratio": odds, "fisher_pval": pval,
        })
        sig = "*" if pval < 0.05 else "ns"
        print(f"    TAGs vs Conserved: OR={odds:.2f}, p={pval:.2e}, "
              f"overlap={n_overlap} {sig}")

    # 3. Drug targets
    drug_df = load_drug_targets()
    if len(drug_df) > 0:
        dt_col = next(
            (c for c in drug_df.columns if "gene" in c.lower() or "target" in c.lower()),
            None,
        )
        if dt_col:
            dt_set = set(drug_df[dt_col].dropna()) & all_genes
            odds, pval, n_overlap = fisher_test(tag_genes, dt_set, all_genes)
            tests.append({
                "gene_set": "Drug_targets",
                "n_tags": len(tag_genes), "n_reference": len(dt_set),
                "n_overlap": n_overlap, "odds_ratio": odds, "fisher_pval": pval,
            })
            sig = "*" if pval < 0.05 else "ns"
            print(f"    TAGs vs Drug targets: OR={odds:.2f}, p={pval:.2e}, "
                  f"overlap={n_overlap} {sig}")

    # 4. COLOC (GWAS loci)
    coloc_df = load_coloc_results(pp4_threshold=0.5)
    if len(coloc_df) > 0:
        gene_cols = ["gene_symbol", "symbol", "gene_name", "gene"]
        coloc_gene_col = next((c for c in gene_cols if c in coloc_df.columns), None)
        if coloc_gene_col:
            coloc_set = set(coloc_df[coloc_gene_col].dropna()) & all_genes
            odds, pval, n_overlap = fisher_test(tag_genes, coloc_set, all_genes)
            tests.append({
                "gene_set": "COLOC_GWAS_loci",
                "n_tags": len(tag_genes), "n_reference": len(coloc_set),
                "n_overlap": n_overlap, "odds_ratio": odds, "fisher_pval": pval,
            })
            sig = "*" if pval < 0.05 else "ns"
            print(f"    TAGs vs COLOC loci: OR={odds:.2f}, p={pval:.2e}, "
                  f"overlap={n_overlap} {sig}")

    # 5. Positive TAGs vs negative TAGs (test for directional bias)
    if "direction" in tag_df.columns and "is_tag" in tag_df.columns:
        gene_col_tag = "gene" if "gene" in tag_df.columns else tag_df.index.name
        if gene_col_tag and gene_col_tag in tag_df.columns:
            pos_tags = set(tag_df[(tag_df["is_tag"]) & (tag_df["direction"] == "positive")][gene_col_tag]) & all_genes
            neg_tags = set(tag_df[(tag_df["is_tag"]) & (tag_df["direction"] == "negative")][gene_col_tag]) & all_genes
        else:
            pos_tags = set(tag_df[(tag_df["is_tag"]) & (tag_df["direction"] == "positive")].index) & all_genes
            neg_tags = set(tag_df[(tag_df["is_tag"]) & (tag_df["direction"] == "negative")].index) & all_genes

        # Test positive TAGs vs dream up-DEGs
        try:
            dream_up = load_dream_degs(padj_thresh=0.1, lfc_thresh=0.0)
            gene_col_d = "symbol" if "symbol" in dream_up.columns else dream_up.columns[0]
            up_set = set(dream_up[dream_up["logFC"] > 0][gene_col_d].dropna()) & all_genes
            if pos_tags:
                odds, pval, n_overlap = fisher_test(pos_tags, up_set, all_genes)
                tests.append({
                    "gene_set": "Positive_TAGs_vs_UpDEGs",
                    "n_tags": len(pos_tags), "n_reference": len(up_set),
                    "n_overlap": n_overlap, "odds_ratio": odds, "fisher_pval": pval,
                })
        except Exception:
            pass

    tests_df = pd.DataFrame(tests)
    if len(tests_df) > 0:
        from statsmodels.stats.multitest import multipletests
        _, padj, _, _ = multipletests(tests_df["fisher_pval"].fillna(1.0).values, method="fdr_bh")
        tests_df["fisher_padj_bh"] = padj

    return tests_df


def main():
    print_header("14e: DPT Niche-Trajectory Integration with Multi-Evidence Atlas "
                 "(ONTraC fallback; CLI unavailable)")

    config = load_config()
    output_dir = RESULTS_DIR / "ontrac"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load DPT niche-trajectory results (ONTraC fallback; CLI unavailable)
    print_step("Loading DPT niche-trajectory results (ONTraC fallback)")
    nt_df = load_nt_scores()
    tag_df = load_tags()
    print(f"  NT scores: {len(nt_df)} spots")
    if len(tag_df) > 0:
        n_tags = tag_df["is_tag"].sum() if "is_tag" in tag_df.columns else len(tag_df)
        print(f"  TAGs: {n_tags}")

    # Load multi-evidence atlas
    print_step("Loading multi-evidence atlas")
    atlas = load_multi_evidence_atlas()
    n_cols_before = len(atlas.columns)
    print(f"  Atlas: {len(atlas)} genes x {n_cols_before} columns")

    # Compute per-gene NT statistics from spatial expression data
    print_step("Computing per-gene NT statistics")
    try:
        adata = load_deconvolved_adata()
        gene_stats = compute_per_gene_nt_stats(adata, nt_df)
        print(f"  Gene stats computed for {len(gene_stats)} genes")
    except Exception as e:
        print(f"  WARNING: Could not compute per-gene stats: {e}")
        gene_stats = pd.DataFrame()

    # Add columns to atlas (DPT-derived; column names kept for 17a compatibility)
    print_step("Adding DPT niche-trajectory columns to atlas (ONTraC fallback)")
    atlas = add_atlas_columns(atlas, tag_df, gene_stats)
    n_new = len(atlas.columns) - n_cols_before
    print(f"  Added {n_new} columns")

    # Report coverage
    new_cols = [c for c in atlas.columns if c.startswith("spatial_") and c not in
                ["spatial_morans_i", "spatial_is_svg", "spatial_zonation_class",
                 "spatial_spearman_rho", "spatial_svg_category",
                 "spatial_hep_spatial_fc", "spatial_hep_validated",
                 "spatial_hep_wilcoxon_padj_bh", "spatial_coexpr_module",
                 "spatial_zone_de_class", "spatial_zone_de_n_sig_bins",
                 "spatial_coloc_zone", "spatial_niche_domain", "spatial_niche_fc"]]
    ontrac_cols = [c for c in atlas.columns if "nt_score" in c or "niche_cluster" in c
                   or "is_tag" in c or "tag_direction" in c]
    for col in ontrac_cols:
        n_annotated = atlas[col].notna().sum()
        if col == "spatial_is_tag":
            n_annotated = atlas[col].sum()
        print(f"    {col}: {n_annotated} genes annotated")

    # Save atlas columns (not full atlas — let 06_integration.py merge)
    sym_col = "human_symbol" if "human_symbol" in atlas.columns else (
        "symbol" if "symbol" in atlas.columns else atlas.columns[0]
    )
    atlas_export = atlas[[sym_col] + ontrac_cols].copy()
    save_csv(atlas_export, "ontrac_atlas_columns.csv", subdir="ontrac")

    # Enrichment tests
    print_step("Running enrichment tests")
    enrichment = run_enrichment_tests(tag_df, atlas)
    if len(enrichment) > 0:
        save_csv(enrichment, "enrichment_tests.csv", subdir="ontrac")
        print(f"\n  Enrichment summary:")
        for _, row in enrichment.iterrows():
            padj_str = f" (BH: {row['fisher_padj_bh']:.2e})" if "fisher_padj_bh" in row.index else ""
            sig = "***" if row["fisher_pval"] < 0.001 else "**" if row["fisher_pval"] < 0.01 else "*" if row["fisher_pval"] < 0.05 else "ns"
            print(f"    {row['gene_set']}: OR={row['odds_ratio']:.2f}, "
                  f"p={row['fisher_pval']:.2e}{padj_str} {sig} "
                  f"({row['n_overlap']}/{row['n_tags']})")

    print_header("14e: Complete (DPT diffusion-pseudotime; ONTraC fallback)")


if __name__ == "__main__":
    main()
