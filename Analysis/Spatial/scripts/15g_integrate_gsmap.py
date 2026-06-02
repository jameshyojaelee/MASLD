#!/usr/bin/env python3
"""
15g_integrate_gsmap.py — Integrate gsMap results into multi-evidence atlas.

Adds gsMap GWAS-spatial risk-localization columns to the multi-evidence atlas.

THE FIX (findings F101 + F104):
  - F101 (path): gsMap quick_mode writes spatial-LDSC scores to
    ``spatial_ldsc/`` (NOT ``ldsc/``).
  - F104 (gene mapping): the spatial_ldsc files are PER-SPOT (columns
    ``spot, beta, se, z, p``) with NO ``gene`` column, so the old
    ``groupby("gene")`` skipped everything and produced no output. Gene-level
    risk localization must instead be derived by combining the per-spot LDSC
    p-values with the per-spot-per-gene marker scores from gsMap's
    ``latent_to_gene/{sample}_gene_marker_score.feather`` (genes x spots).

    A gene is "risk-localized" for a trait if it is a top marker of the spots
    that carry that trait's LDSC risk signal. This logic lives in
    ``15f_gsmap_analysis.gene_risk_localization`` and is reused here (15f
    writes ``gene_risk_localization.csv``; 15g recomputes it directly if that
    CSV is absent, so 15g is runnable standalone).

Output ``gsmap_atlas_columns.csv`` columns (joined by 17a on ``human_symbol``):
  - human_symbol            : gene symbol (merge key for 17a)
  - gsmap_is_risk_localized : bool, gene is risk-localized for >= 1 trait
  - gsmap_n_traits          : int, number of distinct traits localized
  - gsmap_top_trait         : str, trait with the strongest marker enrichment
  - gsmap_min_pval          : float, best (smallest) LDSC spot p-value among the
                              risk spots where the gene is a top marker
                              (< 0.05 by construction for localized genes)

SLURM: --partition=cpu --cpus=4 --mem=32G --time=1:00:00
"""

import pathlib
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).parent))
from spatial_utils import (
    PROJECT_ROOT, RESULTS_DIR, load_config, load_multi_evidence_atlas,
    save_csv, print_header, print_step,
)

# Reuse the spot->gene risk-localization machinery from 15f so the two scripts
# can never diverge on schema/logic.
from importlib import import_module
_m15f = import_module("15f_gsmap_analysis")
load_gsmap_results = _m15f.load_gsmap_results
load_gsmap_marker_scores = _m15f.load_gsmap_marker_scores
gene_risk_localization = _m15f.gene_risk_localization

GSMAP_DIR = RESULTS_DIR / "gsmap"

# Significance threshold (spot-level LDSC p) defining a risk-localized gene
SPATIAL_SIG_THRESHOLD = 0.05


def load_gene_risk_localization():
    """Return the long gene x trait risk-localization table.

    Prefers ``gene_risk_localization.csv`` written by 15f; falls back to
    recomputing it directly from the on-disk gsMap output so 15g can run even
    if 15f was not (the logic is identical — same 15f functions).

    Columns: dataset, gene, trait, enrichment, min_pval, n_samples_flagged.
    """
    cached = GSMAP_DIR / "gene_risk_localization.csv"
    if cached.exists():
        df = pd.read_csv(cached)
        # save_csv writes index=True -> drop the leaked unnamed index column
        df = df.loc[:, ~df.columns.str.match(r"^Unnamed")]
        if {"gene", "trait", "enrichment", "min_pval"}.issubset(df.columns):
            print(f"    Loaded cached gene risk localization: {len(df)} rows "
                  f"from {cached.name}")
            return df
        print(f"    WARNING: {cached.name} missing expected columns — recomputing")

    print("    Recomputing gene risk localization from gsMap output...")
    gsmap_results = load_gsmap_results()
    marker_scores = load_gsmap_marker_scores()
    return gene_risk_localization(gsmap_results, marker_scores)


# Discriminative threshold: a gene is "risk-localized" only if it localizes
# across >= this many of the 7 liver GWAS traits. A >=1-trait flag covered 74%
# of tested genes (near-universal, useless as a convergence vote); >=5 (majority
# of traits) gives a ~19% discriminative set.
GSMAP_MIN_TRAITS = 5


def build_atlas_columns(gene_risk_df):
    """Collapse the long gene x trait table to one row per gene.

    For each gene, aggregates across traits (and the datasets in which it was
    flagged):
      - gsmap_is_risk_localized = True (every gene in this table is localized
        for >= 1 trait)
      - gsmap_n_traits          = number of distinct traits localized
      - gsmap_top_trait         = trait with the maximum marker enrichment
      - gsmap_min_pval          = smallest LDSC spot p across its risk spots
    """
    records = []
    for gene, sub in gene_risk_df.groupby("gene"):
        # Strongest trait by marker enrichment (max over datasets/samples)
        trait_enrich = sub.groupby("trait")["enrichment"].max()
        top_trait = trait_enrich.idxmax()
        n_traits = int(sub["trait"].nunique())
        records.append({
            "human_symbol": gene,
            "gsmap_is_risk_localized": n_traits >= GSMAP_MIN_TRAITS,
            "gsmap_n_traits": n_traits,
            "gsmap_top_trait": str(top_trait),
            "gsmap_min_pval": float(sub["min_pval"].min()),
        })
    out = pd.DataFrame(records)
    if len(out) > 0:
        out = out.sort_values(
            ["gsmap_n_traits", "gsmap_min_pval"], ascending=[False, True]
        ).reset_index(drop=True)
    return out


def main():
    print_header("15g: Integrate gsMap into Multi-Evidence Atlas")

    config = load_config()

    # Load multi-evidence atlas
    atlas = load_multi_evidence_atlas()
    symbol_col = ("human_symbol" if "human_symbol" in atlas.columns
                  else "symbol" if "symbol" in atlas.columns
                  else atlas.columns[0])
    n_cols_before = len(atlas.columns)
    print(f"  Atlas: {len(atlas)} genes x {n_cols_before} columns")

    # Derive gene-level GWAS risk localization (spot -> gene; F104 fix)
    print("\n  Loading gsMap gene-level risk localization...")
    gene_risk_df = load_gene_risk_localization()

    if gene_risk_df is None or len(gene_risk_df) == 0:
        print("  ERROR: No gsMap gene-level risk localization. Run 15e + 15f first.")
        sys.exit(1)

    print(f"  Gene x trait rows: {len(gene_risk_df)} "
          f"({gene_risk_df['gene'].nunique()} genes, "
          f"{gene_risk_df['trait'].nunique()} traits)")

    # Collapse to one row per gene -> atlas columns
    gsmap_cols_df = build_atlas_columns(gene_risk_df)
    print(f"  Risk-localized genes: {len(gsmap_cols_df)}")

    gsmap_value_cols = [
        "gsmap_is_risk_localized", "gsmap_n_traits",
        "gsmap_top_trait", "gsmap_min_pval",
    ]

    # Standalone atlas-columns file (merge key = human_symbol; 17a joins on it)
    save_csv(
        gsmap_cols_df.set_index("human_symbol"),
        "gsmap_atlas_columns.csv", subdir="gsmap",
    )

    # Also write a full atlas annotated with the gsMap columns (left-join;
    # genes absent from the gsMap table get is_risk_localized=False / NaN).
    merged = atlas.merge(
        gsmap_cols_df, how="left",
        left_on=symbol_col, right_on="human_symbol",
    )
    if "human_symbol" != symbol_col and "human_symbol" in merged.columns:
        merged = merged.drop(columns=["human_symbol"])
    merged["gsmap_is_risk_localized"] = (
        merged["gsmap_is_risk_localized"].fillna(False).astype(bool)
    )
    merged["gsmap_n_traits"] = merged["gsmap_n_traits"].fillna(0).astype(int)

    n_localized = int(merged["gsmap_is_risk_localized"].sum())
    n_cols_after = len(merged.columns)
    print(f"\n  Atlas genes flagged risk-localized: {n_localized}")
    print(f"  Added {n_cols_after - n_cols_before} gsMap columns to atlas")

    atlas_out = GSMAP_DIR / "multi_evidence_atlas_with_gsmap.csv"
    merged.to_csv(atlas_out, index=False)
    print(f"  Full atlas saved: {atlas_out} "
          f"({len(merged)} genes x {n_cols_after} columns)")

    print_header("15g: Complete")


if __name__ == "__main__":
    main()
