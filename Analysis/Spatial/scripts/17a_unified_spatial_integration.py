#!/usr/bin/env python3
"""
17a_unified_spatial_integration.py — Merge all 4 enhanced spatial method outputs
into the master multi-evidence atlas and compute cross-method convergence.

Methods integrated:
  - MultiSP: spatial domains (multisp_domains.csv, validation_*.csv)
  - ONTraC: niche trajectory-associated genes (ontrac_atlas_columns.csv)
  - gsMap: GWAS-zonation enrichment (gsmap_atlas_columns.csv)
  - COMMOT: signaling-regulated genes (commot_atlas_columns.csv)

Outputs:
  - Updated multi_evidence_atlas.csv (in-place update of S6 spatial columns)
  - results/spatial_integration/cross_method_convergence.csv
  - results/spatial_integration/method_coverage_summary.csv

SLURM: --partition=cpu --cpus=8 --mem=64G --time=2:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config, load_multi_evidence_atlas,
    load_conserved, load_dream_degs, load_drug_targets,
    save_csv, print_header, print_step, fisher_test,
)


# ── Method Loaders ──────────────────────────────────────────────────────────

def _drop_pandas_index_leak(df):
    """T0.3 (2026-04-22): drop leading `Unnamed: 0` column written by save_csv(index=True)
    in spatial_utils. Without this, downstream atlas merge suffixes it with the method
    prefix (e.g. `ontrac_Unnamed: 0`) and leaks a spurious integer column into the atlas.
    """
    if df is None:
        return None
    leak_cols = [c for c in df.columns if c.startswith("Unnamed:")]
    if leak_cols:
        df = df.drop(columns=leak_cols)
    return df


def load_ontrac_atlas_columns():
    """Load ONTraC per-gene atlas columns."""
    path = RESULTS_DIR / "ontrac" / "ontrac_atlas_columns.csv"
    if not path.exists():
        print(f"  WARNING: {path} not found")
        return None
    df = _drop_pandas_index_leak(pd.read_csv(path))
    print(f"  ONTraC: {len(df)} genes loaded from {path.name}")
    return df


def load_gsmap_atlas_columns():
    """Load gsMap per-gene atlas columns."""
    path = RESULTS_DIR / "gsmap" / "gsmap_atlas_columns.csv"
    if not path.exists():
        print(f"  WARNING: {path} not found")
        return None
    df = _drop_pandas_index_leak(pd.read_csv(path))
    print(f"  gsMap: {len(df)} genes loaded from {path.name}")
    return df


def load_commot_atlas_columns():
    """Load COMMOT per-gene atlas columns."""
    path = RESULTS_DIR / "commot" / "commot_atlas_columns.csv"
    if not path.exists():
        print(f"  WARNING: {path} not found")
        return None
    df = _drop_pandas_index_leak(pd.read_csv(path))
    print(f"  COMMOT: {len(df)} genes loaded from {path.name}")
    return df


def load_multisp_results():
    """Load MultiSP domain classification results."""
    path = RESULTS_DIR / "multisp" / "multisp_domains.csv"
    if not path.exists():
        print(f"  WARNING: {path} not found")
        return None
    df = pd.read_csv(path)
    print(f"  MultiSP: {len(df)} entries loaded from {path.name}")
    return df


# ── Integration Logic ───────────────────────────────────────────────────────

def identify_gene_column(df):
    """Find the gene symbol column in a DataFrame."""
    for col in ["human_symbol", "symbol", "gene", "gene_symbol"]:
        if col in df.columns:
            return col
    # Check if index looks like gene symbols
    if df.index.dtype == object and not df.index.str.startswith("ENS").any():
        return None  # index is gene symbol
    return df.columns[0]


def merge_method_columns(atlas, method_df, method_name, gene_col_atlas):
    """Left-join per-method atlas columns into master atlas."""
    if method_df is None:
        return atlas

    gene_col_method = identify_gene_column(method_df)
    if gene_col_method is None:
        # Gene symbols are in the index
        method_df = method_df.copy()
        method_df["_merge_key"] = method_df.index
        gene_col_method = "_merge_key"

    # Identify columns to merge (exclude the gene key itself)
    merge_cols = [c for c in method_df.columns if c != gene_col_method]
    if not merge_cols:
        print(f"  WARNING: No data columns in {method_name}")
        return atlas

    # Rename columns to avoid collisions (prefix with method name if needed)
    # Per-method atlas columns should already be prefixed (e.g., ontrac_*, gsmap_*)
    # Check if already prefixed; if not, add prefix
    prefix = method_name.lower() + "_"
    rename_map = {}
    for col in merge_cols:
        if not col.startswith(prefix) and not col.startswith("spatial_"):
            rename_map[col] = prefix + col

    method_df = method_df.rename(columns=rename_map)
    merge_cols = [rename_map.get(c, c) for c in merge_cols]

    # Build mapping: gene symbol -> row values
    merge_df = method_df[[gene_col_method] + merge_cols].copy()
    merge_df[gene_col_method] = merge_df[gene_col_method].astype(str)
    merge_df = merge_df.drop_duplicates(subset=[gene_col_method], keep="first")

    # Left join
    n_before = len(atlas.columns)
    atlas = atlas.merge(
        merge_df,
        left_on=gene_col_atlas,
        right_on=gene_col_method,
        how="left",
        suffixes=("", f"_{method_name}"),
    )
    # Drop duplicate merge key if created
    if gene_col_method != gene_col_atlas and gene_col_method in atlas.columns:
        atlas = atlas.drop(columns=[gene_col_method])

    n_added = len(atlas.columns) - n_before
    n_annotated = atlas[merge_cols[0]].notna().sum() if merge_cols else 0
    print(f"  {method_name}: {n_added} columns added, {n_annotated} genes annotated")

    return atlas


def build_multisp_gene_flags(atlas, multisp_df, gene_col):
    """Map MultiSP domain-level results to per-gene flags.

    Genes are flagged as 'in disease domain' if they are differentially
    expressed in a disease-enriched MultiSP domain.
    """
    if multisp_df is None:
        return atlas

    # Look for per-gene validation results first
    val_path = RESULTS_DIR / "multisp" / "validation_gene_flags.csv"
    if val_path.exists():
        flags = pd.read_csv(val_path)
        flag_gene_col = identify_gene_column(flags)
        if flag_gene_col is None:
            flags["_key"] = flags.index
            flag_gene_col = "_key"
        for col in flags.columns:
            if col != flag_gene_col and col.startswith("multisp_"):
                atlas[col] = atlas[gene_col].map(
                    dict(zip(flags[flag_gene_col], flags[col]))
                )
                n_ann = atlas[col].notna().sum()
                print(f"  MultiSP (validation): {col} -> {n_ann} genes")
        return atlas

    # Fallback: look for domain_gene_markers or domain_marker_genes
    marker_path = RESULTS_DIR / "multisp" / "domain_marker_genes.csv"
    if not marker_path.exists():
        marker_path = RESULTS_DIR / "multisp" / "multisp_domain_markers.csv"
    if marker_path.exists():
        markers = pd.read_csv(marker_path)
        marker_gene_col = identify_gene_column(markers)
        if marker_gene_col is None:
            markers["_key"] = markers.index
            marker_gene_col = "_key"
        # Flag genes that are markers for any disease-associated domain
        if "domain" in markers.columns:
            gene_domain_map = dict(zip(markers[marker_gene_col], markers["domain"]))
            atlas["multisp_domain"] = atlas[gene_col].map(gene_domain_map)
            atlas["multisp_in_disease_domain"] = atlas["multisp_domain"].notna()
            n_ann = atlas["multisp_in_disease_domain"].sum()
            print(f"  MultiSP (domain markers): {n_ann} genes in disease domain")
    else:
        # Minimal fallback: just flag presence in MultiSP output
        mp_gene_col = identify_gene_column(multisp_df)
        if mp_gene_col:
            genes_in_multisp = set(multisp_df[mp_gene_col].dropna())
        else:
            genes_in_multisp = set(multisp_df.index)
        atlas["multisp_in_disease_domain"] = atlas[gene_col].isin(genes_in_multisp)
        n_ann = atlas["multisp_in_disease_domain"].sum()
        print(f"  MultiSP (presence): {n_ann} genes flagged")

    return atlas


# ── Convergence Scoring ─────────────────────────────────────────────────────

def compute_convergence(atlas, gene_col):
    """Compute per-gene spatial method convergence score.

    For each gene, count how many of the 4 methods provide evidence:
      - ONTraC: trajectory-associated gene (TAG)
      - gsMap: risk-localized gene
      - COMMOT: signaling-regulated gene
      - MultiSP: in disease-enriched domain

    Produces:
      - spatial_n_methods: integer count 0-4
      - spatial_convergence_score: weighted combination
    """
    convergence = pd.DataFrame(index=atlas.index)
    convergence[gene_col] = atlas[gene_col]

    # 1. ONTraC TAG indicator
    ontrac_cols = [c for c in atlas.columns if c.startswith("ontrac_")]
    if ontrac_cols:
        # Look for a binary TAG indicator or padj column
        tag_col = None
        for candidate in ["ontrac_is_tag", "ontrac_tag", "ontrac_trajectory_associated"]:
            if candidate in atlas.columns:
                tag_col = candidate
                break
        if tag_col:
            convergence["ontrac_evidence"] = atlas[tag_col].fillna(False).astype(bool)
        else:
            # Try padj-based: significant if ontrac padj < threshold
            padj_candidates = [c for c in ontrac_cols if "padj" in c or "pval" in c or "fdr" in c]
            if padj_candidates:
                convergence["ontrac_evidence"] = atlas[padj_candidates[0]].fillna(1.0) < 0.05
            else:
                # Any non-null ontrac data = evidence
                convergence["ontrac_evidence"] = atlas[ontrac_cols[0]].notna()
    else:
        convergence["ontrac_evidence"] = False

    # 2. gsMap risk-localized indicator
    gsmap_cols = [c for c in atlas.columns if c.startswith("gsmap_")]
    if gsmap_cols:
        risk_col = None
        for candidate in ["gsmap_is_risk_localized", "gsmap_risk_localized",
                          "gsmap_significant"]:
            if candidate in atlas.columns:
                risk_col = candidate
                break
        if risk_col:
            convergence["gsmap_evidence"] = atlas[risk_col].fillna(False).astype(bool)
        else:
            padj_candidates = [c for c in gsmap_cols if "padj" in c or "pval" in c or "fdr" in c]
            if padj_candidates:
                convergence["gsmap_evidence"] = atlas[padj_candidates[0]].fillna(1.0) < 0.05
            else:
                convergence["gsmap_evidence"] = atlas[gsmap_cols[0]].notna()
    else:
        convergence["gsmap_evidence"] = False

    # 3. COMMOT signaling-regulated indicator
    commot_cols = [c for c in atlas.columns if c.startswith("commot_")]
    if commot_cols:
        reg_col = None
        for candidate in ["commot_is_signaling_regulated", "commot_signaling_regulated",
                          "commot_significant"]:
            if candidate in atlas.columns:
                reg_col = candidate
                break
        if reg_col:
            convergence["commot_evidence"] = atlas[reg_col].fillna(False).astype(bool)
        else:
            padj_candidates = [c for c in commot_cols if "padj" in c or "pval" in c or "fdr" in c]
            if padj_candidates:
                convergence["commot_evidence"] = atlas[padj_candidates[0]].fillna(1.0) < 0.05
            else:
                convergence["commot_evidence"] = atlas[commot_cols[0]].notna()
    else:
        convergence["commot_evidence"] = False

    # 4. MultiSP disease domain indicator
    if "multisp_in_disease_domain" in atlas.columns:
        convergence["multisp_evidence"] = atlas["multisp_in_disease_domain"].fillna(False).astype(bool)
    else:
        convergence["multisp_evidence"] = False

    # Count methods with evidence
    evidence_cols = ["ontrac_evidence", "gsmap_evidence", "commot_evidence", "multisp_evidence"]
    convergence["spatial_n_methods"] = convergence[evidence_cols].sum(axis=1).astype(int)

    # Weighted convergence score: equal weights (0.25 each)
    convergence["spatial_convergence_score"] = convergence[evidence_cols].astype(float).mean(axis=1)

    return convergence


# ── Enrichment Tests ────────────────────────────────────────────────────────

def run_enrichment_tests(atlas, convergence, gene_col):
    """Fisher's exact tests: genes with >=3/4 spatial methods vs key gene sets."""
    converged_genes = set(
        convergence.loc[convergence["spatial_n_methods"] >= 3, gene_col].dropna()
    )
    if len(converged_genes) < 5:
        print(f"  Only {len(converged_genes)} genes with >=3 methods, "
              f"trying >=2 threshold")
        converged_genes = set(
            convergence.loc[convergence["spatial_n_methods"] >= 2, gene_col].dropna()
        )
    if len(converged_genes) < 5:
        print("  Insufficient converged genes for enrichment tests")
        return pd.DataFrame()

    all_genes = set(atlas[gene_col].dropna())
    tests = []

    # 1. Conserved
    cc = set(load_conserved()) & all_genes
    if cc:
        odds, pval, overlap = fisher_test(converged_genes, cc, all_genes)
        tests.append({
            "test_set": "Conserved",
            "n_converged": len(converged_genes),
            "n_test_set": len(cc),
            "n_overlap": overlap,
            "odds_ratio": odds,
            "fisher_pval": pval,
        })

    # 2. Dream DEGs
    if "dream_padj" in atlas.columns:
        dream_genes = set(atlas.loc[atlas["dream_padj"] < 0.1, gene_col].dropna()) & all_genes
        if dream_genes:
            odds, pval, overlap = fisher_test(converged_genes, dream_genes, all_genes)
            tests.append({
                "test_set": "Dream_DEGs",
                "n_converged": len(converged_genes),
                "n_test_set": len(dream_genes),
                "n_overlap": overlap,
                "odds_ratio": odds,
                "fisher_pval": pval,
            })

    # 3. Drug targets
    drug_df = load_drug_targets()
    if len(drug_df) > 0:
        drug_gene_col = identify_gene_column(drug_df)
        if drug_gene_col:
            drug_genes = set(drug_df[drug_gene_col].dropna()) & all_genes
        else:
            drug_genes = set(drug_df.index) & all_genes
        if drug_genes:
            odds, pval, overlap = fisher_test(converged_genes, drug_genes, all_genes)
            tests.append({
                "test_set": "Drug_targets",
                "n_converged": len(converged_genes),
                "n_test_set": len(drug_genes),
                "n_overlap": overlap,
                "odds_ratio": odds,
                "fisher_pval": pval,
            })

    # 4. DGIdb druggable
    if "dgidb_druggable" in atlas.columns:
        druggable = set(atlas.loc[atlas["dgidb_druggable"] == True, gene_col].dropna()) & all_genes
        if druggable:
            odds, pval, overlap = fisher_test(converged_genes, druggable, all_genes)
            tests.append({
                "test_set": "Druggable_genes",
                "n_converged": len(converged_genes),
                "n_test_set": len(druggable),
                "n_overlap": overlap,
                "odds_ratio": odds,
                "fisher_pval": pval,
            })

    # 5. Top 5% multi-evidence
    if "layers_active" in atlas.columns:
        threshold = atlas["layers_active"].quantile(0.95)
        top_me = set(atlas.loc[atlas["layers_active"] >= threshold, gene_col].dropna()) & all_genes
        if top_me:
            odds, pval, overlap = fisher_test(converged_genes, top_me, all_genes)
            tests.append({
                "test_set": "Top5pct_multi_evidence",
                "n_converged": len(converged_genes),
                "n_test_set": len(top_me),
                "n_overlap": overlap,
                "odds_ratio": odds,
                "fisher_pval": pval,
            })

    tests_df = pd.DataFrame(tests)
    if len(tests_df) > 0:
        _, padj, _, _ = multipletests(tests_df["fisher_pval"].values, method="fdr_bh")
        tests_df["fisher_padj_bh"] = padj

    return tests_df


# ── Coverage Summary ────────────────────────────────────────────────────────

def compute_coverage_summary(atlas, convergence, gene_col):
    """Compute per-method and combined coverage statistics."""
    n_total = len(atlas)
    rows = []

    # Per-method coverage
    for method, col in [
        ("ONTraC", "ontrac_evidence"),
        ("gsMap", "gsmap_evidence"),
        ("COMMOT", "commot_evidence"),
        ("MultiSP", "multisp_evidence"),
    ]:
        if col in convergence.columns:
            n_with = convergence[col].sum()
            rows.append({
                "method": method,
                "n_genes_with_evidence": int(n_with),
                "pct_coverage": n_with / n_total * 100,
            })

    # Combined coverage (any method)
    evidence_cols = [c for c in ["ontrac_evidence", "gsmap_evidence",
                                 "commot_evidence", "multisp_evidence"]
                     if c in convergence.columns]
    if evidence_cols:
        any_evidence = convergence[evidence_cols].any(axis=1).sum()
        rows.append({
            "method": "Any_method",
            "n_genes_with_evidence": int(any_evidence),
            "pct_coverage": any_evidence / n_total * 100,
        })

    # Per n_methods breakdown
    for k in range(5):
        n_k = (convergence["spatial_n_methods"] == k).sum()
        rows.append({
            "method": f"n_methods={k}",
            "n_genes_with_evidence": int(n_k),
            "pct_coverage": n_k / n_total * 100,
        })

    return pd.DataFrame(rows)


# ── Main ────────────────────────────────────────────────────────────────────

def main():
    print_header("17a: Unified Spatial Integration")

    config = load_config()

    # ── 1. Load current multi-evidence atlas ─────────────────────────────
    print_step("Loading multi-evidence atlas")
    atlas = load_multi_evidence_atlas()
    gene_col = "human_symbol" if "human_symbol" in atlas.columns else "symbol"
    n_cols_before = len(atlas.columns)
    print(f"  Atlas: {len(atlas)} genes x {n_cols_before} columns")

    # Record old S6 coverage
    old_spatial_cols = [c for c in atlas.columns if c.startswith("spatial_")]
    old_s6_coverage = 0
    if old_spatial_cols:
        old_s6_coverage = atlas[old_spatial_cols].notna().any(axis=1).sum()
    print(f"  Old S6 coverage: {old_s6_coverage} genes ({old_s6_coverage / len(atlas) * 100:.1f}%)")

    # ── 2. Load per-method atlas columns ─────────────────────────────────
    print_step("Loading per-method results")
    ontrac_df = load_ontrac_atlas_columns()
    gsmap_df = load_gsmap_atlas_columns()
    commot_df = load_commot_atlas_columns()
    multisp_df = load_multisp_results()

    n_methods_found = sum(x is not None for x in [ontrac_df, gsmap_df, commot_df, multisp_df])
    print(f"\n  Methods with data: {n_methods_found}/4")

    if n_methods_found == 0:
        print("  ERROR: No method results found. Exiting.")
        print_header("17a: Complete (no data)")
        return

    # ── 3. Merge all columns into atlas ──────────────────────────────────
    print_step("Merging method columns into atlas")
    atlas = merge_method_columns(atlas, ontrac_df, "ontrac", gene_col)
    atlas = merge_method_columns(atlas, gsmap_df, "gsmap", gene_col)
    atlas = merge_method_columns(atlas, commot_df, "commot", gene_col)
    atlas = build_multisp_gene_flags(atlas, multisp_df, gene_col)

    # ── 4. Compute cross-method convergence ──────────────────────────────
    print_step("Computing cross-method convergence")
    convergence = compute_convergence(atlas, gene_col)

    # Add convergence columns to atlas
    atlas["spatial_n_methods"] = convergence["spatial_n_methods"].values
    atlas["spatial_convergence_score"] = convergence["spatial_convergence_score"].values

    n_cols_after = len(atlas.columns)
    print(f"  Added {n_cols_after - n_cols_before} new columns (total: {n_cols_after})")

    # ── 5. Update S6 spatial coverage statistics ─────────────────────────
    print_step("Computing updated S6 coverage")
    # Coverage = genes with >=1 method of ACTUAL evidence. Do NOT count the
    # never-NaN summary columns (spatial_n_methods / spatial_convergence_score),
    # which previously inflated coverage to a meaningless 100% (17a coverage bug).
    summary_cols = {"spatial_n_methods", "spatial_convergence_score"}
    new_spatial_cols = [c for c in atlas.columns
                        if (c.startswith("spatial_") or c.startswith("ontrac_")
                            or c.startswith("gsmap_") or c.startswith("commot_")
                            or c.startswith("multisp_"))
                        and c not in summary_cols]
    new_s6_coverage = int((convergence["spatial_n_methods"] >= 1).sum())
    print(f"  New S6 coverage: {new_s6_coverage} genes ({new_s6_coverage / len(atlas) * 100:.1f}%)")
    print(f"  Change: {old_s6_coverage} -> {new_s6_coverage} "
          f"(+{new_s6_coverage - old_s6_coverage})")

    # ── 6. Enrichment tests ──────────────────────────────────────────────
    print_step("Running enrichment tests (convergent genes vs key gene sets)")
    enrichment = run_enrichment_tests(atlas, convergence, gene_col)
    if len(enrichment) > 0:
        for _, row in enrichment.iterrows():
            sig = ("***" if row["fisher_pval"] < 0.001
                   else "**" if row["fisher_pval"] < 0.01
                   else "*" if row["fisher_pval"] < 0.05
                   else "ns")
            padj_str = ""
            if "fisher_padj_bh" in row.index:
                padj_str = f" (BH: {row['fisher_padj_bh']:.2e})"
            print(f"    {row['test_set']}: OR={row['odds_ratio']:.2f}, "
                  f"p={row['fisher_pval']:.2e}{padj_str} {sig} "
                  f"({row['n_overlap']}/{row['n_converged']})")

    # ── 7. Save outputs ──────────────────────────────────────────────────
    print_step("Saving outputs")

    # Save updated master atlas
    atlas_path = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
    atlas.to_csv(atlas_path, index=False)
    print(f"  Saved atlas: {atlas_path} ({len(atlas)} genes x {len(atlas.columns)} cols)")

    # Save convergence table
    outdir = RESULTS_DIR / "spatial_integration"
    outdir.mkdir(parents=True, exist_ok=True)
    save_csv(convergence, "cross_method_convergence.csv", subdir="spatial_integration")

    # Save coverage summary
    coverage = compute_coverage_summary(atlas, convergence, gene_col)
    save_csv(coverage, "method_coverage_summary.csv", subdir="spatial_integration")

    # Save enrichment tests
    if len(enrichment) > 0:
        save_csv(enrichment, "convergence_enrichment_tests.csv", subdir="spatial_integration")

    # ── 8. Print summary ─────────────────────────────────────────────────
    print_step("Summary")
    print(f"\n  Atlas dimensions: {len(atlas)} genes x {len(atlas.columns)} columns")
    print(f"  Columns added: {n_cols_after - n_cols_before}")
    print(f"  S6 coverage: {old_s6_coverage} -> {new_s6_coverage} genes "
          f"({new_s6_coverage / len(atlas) * 100:.1f}%)")
    print(f"\n  Convergence distribution:")
    for k in range(5):
        n_k = (convergence["spatial_n_methods"] == k).sum()
        print(f"    {k} methods: {n_k} genes ({n_k / len(atlas) * 100:.1f}%)")

    print_header("17a: Complete")


if __name__ == "__main__":
    main()
