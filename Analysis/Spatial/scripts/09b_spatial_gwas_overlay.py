#!/usr/bin/env python3
"""
09b_spatial_gwas_overlay.py — Map COLOC hits onto zonation axis.

Maps canonical SuSiE-COLOC hits (PP.H4.susie>0.5 on the PolyFun 28-GWAS atlas;
~538 genes) onto the periportal-pericentral zonation gradient and spatial
domains. Tests zonation enrichment with Fisher's exact tests.

COLOC source repointed 2026-06-01 (F001) from the deprecated 4-GWAS ABF panel
(RNA-seq/results/causal_inference/broadaway_*/coloc_results.csv) to the
canonical GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv, matching
Script 213. MR arm removed 2026-06-01 (MR permanently ditched 2026-04-22).

SLURM: --partition=cpu --cpus=8 --mem=64G --time=2:00:00
"""

import pathlib
import sys
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config, load_multi_evidence_atlas,
    save_csv, print_header, print_step,
)


def load_coloc_results():
    """Load canonical SuSiE-COLOC gene-level results (F001 repoint 2026-06-01).

    Reads GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv (the
    PolyFun-canonical 28-GWAS SuSiE-COLOC table), filtering on
    coloc_best_susie_pp4 > 0.5 with ABF (coloc_best_pp4) fallback — exactly as
    Script 213 does (213:40-49). Replaces the deprecated 4-GWAS ABF panel
    (broadaway_*/coloc_results.csv, mtime Mar 2026) that predated both the
    2026-04-21 SuSiE-canonical rebuild and the 2026-05-06 PolyFun swap.

    Returns one row per colocalizing gene with a 'gene' symbol column, the
    effective 'coloc_best_pp4' (SuSiE-preferred), and a 'gwas' column derived
    from coloc_best_susie_gwas (ABF gwas fallback) so the downstream per-GWAS
    grouping in main() still works.
    """
    coloc_path = (PROJECT_ROOT /
                  "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
    if not coloc_path.exists():
        print(f"  WARNING: {coloc_path} not found")
        return pd.DataFrame()

    df = pd.read_csv(coloc_path)
    df = df[df["gene"].notna() & (df["gene"].astype(str) != "")].copy()

    # Prefer SuSiE PP4 over legacy ABF (matches 213 T0.4 logic).
    if "coloc_best_susie_pp4" in df.columns:
        df["coloc_best_abf_pp4"] = df["coloc_best_pp4"]
        df["coloc_best_pp4"] = df["coloc_best_susie_pp4"].where(
            df["coloc_best_susie_pp4"].notna(), df["coloc_best_pp4"])

    # Best GWAS label: SuSiE best gwas with ABF best-gwas fallback.
    if "coloc_best_susie_gwas" in df.columns:
        df["gwas"] = df["coloc_best_susie_gwas"].where(
            df["coloc_best_susie_gwas"].notna() &
            (df["coloc_best_susie_gwas"].astype(str) != ""),
            df.get("coloc_best_gwas"))
    else:
        df["gwas"] = df.get("coloc_best_gwas")
    df["gwas"] = df["gwas"].fillna("unknown")

    # Canonical colocalization threshold: PP.H4 > 0.5 on the effective posterior.
    df = df[df["coloc_best_pp4"] > 0.5].copy()
    print(f"  SuSiE-COLOC (PP.H4>0.5, SuSiE-preferred): {len(df)} genes "
          f"across {df['gwas'].nunique()} best-GWAS labels")
    return df.reset_index(drop=True)


def map_genes_to_zonation(genes, deg_zon_df):
    """Map gene list to zonation classifications from 04b results."""
    if len(deg_zon_df) == 0:
        return pd.DataFrame()

    gene_col = "gene" if "gene" in deg_zon_df.columns else deg_zon_df.columns[0]
    mapped = deg_zon_df[deg_zon_df[gene_col].isin(genes)].copy()
    return mapped


def fisher_enrichment(genes_of_interest, all_genes_df, zone_class, gene_col="gene"):
    """Fisher's exact test: are genes enriched in a specific zone class?"""
    in_zone = set(all_genes_df[all_genes_df["zonation_class"] == zone_class][gene_col])
    not_in_zone = set(all_genes_df[all_genes_df["zonation_class"] != zone_class][gene_col])
    interest_set = set(genes_of_interest)

    a = len(interest_set & in_zone)
    b = len(interest_set & not_in_zone)
    c = len(in_zone - interest_set)
    d = len(not_in_zone - interest_set)

    if a + b == 0 or a + c == 0:
        return np.nan, np.nan, 0, 0

    odds_ratio, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
    return odds_ratio, pval, a, a + b


def main():
    print_header("09b: Spatial GWAS Integration")

    config = load_config()
    output_dir = RESULTS_DIR / "gwas_spatial"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load zonation classification from 04b
    zon_path = RESULTS_DIR / "zonation" / "deg_zonation_classification.csv"
    if not zon_path.exists():
        print("  ERROR: Run 04b_deg_zonation_mapping.py first")
        sys.exit(1)
    deg_zon = pd.read_csv(zon_path)
    # Resolve the gene-symbol column once: the regenerated (2026-05-29)
    # deg_zonation_classification.csv uses 'human_symbol', older versions used
    # 'gene'. Resolve defensively so the loader repoint does not crash on the
    # current schema.
    zon_gene_col = next((c for c in ["gene", "human_symbol", "symbol"]
                         if c in deg_zon.columns), deg_zon.columns[0])
    print(f"  Loaded zonation classifications: {len(deg_zon)} genes "
          f"(gene col '{zon_gene_col}')")

    # Load COLOC results (canonical SuSiE-COLOC; MR removed 2026-06-01)
    print("\n  Loading SuSiE-COLOC results...")
    coloc_df = load_coloc_results()

    # Map COLOC genes to zonation
    if len(coloc_df) > 0:
        # Get gene symbols from COLOC (may be in different columns)
        gene_cols = ["gene_symbol", "symbol", "gene_name", "gene"]
        coloc_gene_col = next((c for c in gene_cols if c in coloc_df.columns), None)
        if coloc_gene_col is None:
            print("  WARNING: Cannot find gene symbol column in COLOC results")
            coloc_genes = []
        else:
            coloc_genes = coloc_df[coloc_gene_col].dropna().unique().tolist()
            print(f"\n  Total unique COLOC genes: {len(coloc_genes)}")

            # Map to zonation
            coloc_zon = map_genes_to_zonation(coloc_genes, deg_zon)
            print(f"  COLOC genes with zonation: {len(coloc_zon)}")

            if len(coloc_zon) > 0:
                print("\n  COLOC genes by zonation class:")
                for cls, n in coloc_zon["zonation_class"].value_counts().items():
                    print(f"    {cls}: {n}")

            # Per-GWAS zonation mapping
            gwas_zon_records = []
            for gwas in coloc_df["gwas"].unique():
                gwas_genes = coloc_df[coloc_df["gwas"] == gwas][coloc_gene_col].dropna().unique()
                gwas_zon = map_genes_to_zonation(gwas_genes, deg_zon)

                for _, row in gwas_zon.iterrows():
                    gwas_zon_records.append({
                        "gene": row[zon_gene_col],
                        "gwas": gwas,
                        "zonation_class": row["zonation_class"],
                        "spearman_rho": row.get("spearman_rho", np.nan),
                        "kruskal_pval": row.get("kruskal_pval", np.nan),
                    })

            if gwas_zon_records:
                gwas_zon_df = pd.DataFrame(gwas_zon_records)
                save_csv(gwas_zon_df, "coloc_zonation_mapping.csv", subdir="gwas_spatial")

            # Fisher's exact enrichment tests
            print("\n  Zonation enrichment tests (Fisher's exact):")
            enrichment_records = []
            # Derive zone classes from the actual table rather than hardcoding
            # label strings: the 2026-05-29 deg_zonation_classification.csv uses
            # {Periportal, Pericentral, Non-zoned}; older versions used
            # {Periportal-enriched, Pericentral-enriched, Pan-lobular}. Exclude
            # the non-zoned / pan-lobular catch-all from the tested zones.
            present_classes = [c for c in deg_zon["zonation_class"].dropna().unique()
                               if str(c) not in ("Non-zoned", "Pan-lobular")]
            zone_classes = sorted(present_classes)

            # Overall COLOC enrichment
            for zone in zone_classes:
                OR, pval, n_hit, n_total = fisher_enrichment(
                    coloc_genes, deg_zon, zone, gene_col=zon_gene_col
                )
                enrichment_records.append({
                    "gene_set": "All_COLOC", "zone": zone,
                    "odds_ratio": OR, "pval": pval,
                    "n_in_zone": n_hit, "n_total": n_total,
                })
                sig = "*" if pval < 0.05 else ""
                print(f"    All COLOC → {zone}: OR={OR:.2f}, p={pval:.3f} "
                      f"({n_hit}/{n_total}) {sig}")

            # Per-GWAS enrichment
            for gwas in coloc_df["gwas"].unique():
                gwas_genes = coloc_df[coloc_df["gwas"] == gwas][coloc_gene_col].dropna().unique()
                for zone in zone_classes:
                    OR, pval, n_hit, n_total = fisher_enrichment(
                        gwas_genes, deg_zon, zone, gene_col=zon_gene_col
                    )
                    enrichment_records.append({
                        "gene_set": f"COLOC_{gwas}", "zone": zone,
                        "odds_ratio": OR, "pval": pval,
                        "n_in_zone": n_hit, "n_total": n_total,
                    })

            enrichment_df = pd.DataFrame(enrichment_records)
            # BH correction
            pvals = enrichment_df["pval"].fillna(1.0).values
            _, padj, _, _ = multipletests(pvals, method="fdr_bh")
            enrichment_df["padj_bh"] = padj
            save_csv(enrichment_df, "gwas_spatial_enrichment.csv", subdir="gwas_spatial")

    # Cross-reference with disease-emergent SVGs
    svg_diff_path = RESULTS_DIR / "svg" / "differential_svgs.csv"
    if svg_diff_path.exists() and len(coloc_df) > 0 and coloc_gene_col:
        print("\n  Cross-referencing with disease-emergent SVGs...")
        svg_diff = pd.read_csv(svg_diff_path, index_col=0)
        emergent_svgs = set(svg_diff[svg_diff["category"] == "disease_emergent_SVG"].index)
        coloc_set = set(coloc_genes)
        overlap = coloc_set & emergent_svgs
        print(f"    COLOC ∩ disease-emergent SVGs: {len(overlap)}")
        if overlap:
            print(f"    Genes: {sorted(overlap)[:20]}")

        # Hypergeometric test
        from scipy.stats import hypergeom
        M = len(svg_diff)  # total genes
        n = len(emergent_svgs)  # emergent SVGs
        N = len(coloc_set & set(svg_diff.index))  # COLOC in SVG universe
        k = len(overlap)  # COLOC ∩ emergent
        pval = hypergeom.sf(k - 1, M, n, N)
        print(f"    Hypergeometric p = {pval:.3e}")

    # NOTE (2026-06-01): The MR-genes-on-zonation block was removed here.
    # Mendelian Randomization was permanently ditched from the paper on
    # 2026-04-22 (memory/feedback_no_mr.md). TWAS+COLOC+INTACT is the sole
    # sanctioned causal framework, so only the SuSiE-COLOC overlay above runs.

    print_header("09b: Complete")


if __name__ == "__main__":
    main()
