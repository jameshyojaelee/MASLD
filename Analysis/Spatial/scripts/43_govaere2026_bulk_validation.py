#!/usr/bin/env python3
"""
43_govaere2026_bulk_validation.py — Correlate Wave 3 GeoMx + CosMx DE results
against the bulk dream Disease-vs-Control DEG signature and quantify concordance.

Inputs:
  Analysis/Spatial/results/govaere2026/
    geomx_de_sh_vs_pt.csv                       — CD68+ SH vs PT (lmer)
    geomx_de_sh_vs_ls.csv                       — CD68+ SH vs LS (lmer)
    geomx_de_panck_vs_cd68.csv                  — panCK vs CD68 (segment ID)
    geomx_de_cd45_vs_cd68.csv                   — CD45 vs CD68 (segment ID)
    cosmx_de_Hepatocyte_MASH_vs_noMASH.csv      — Hepatocyte (wilcoxon)
    cosmx_de_KC_MASH_vs_noMASH.csv              — KC cluster (incl. MetMac)

  RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/
    dream_results_ashr.csv                      — canonical bulk DEG file
      (cols: gene [ENSG], logFC, t, P.Value, padj, shrunk_logFC, lfsr, symbol)

Outputs (Analysis/Spatial/results/govaere2026/):
  bulk_validation_geomx_sh_vs_pt.csv
  bulk_validation_geomx_sh_vs_ls.csv
  bulk_validation_geomx_panck_vs_cd68.csv
  bulk_validation_geomx_cd45_vs_cd68.csv
  bulk_validation_cosmx_hepatocyte.csv
  bulk_validation_cosmx_kc.csv
  bulk_validation_summary.tsv
  bulk_validation_top10_concordant.tsv

Method (per contrast):
  1. Standardize Govaere DE to (gene_symbol, logFC, padj).
  2. Inner-join with bulk dream on gene_symbol.
  3. Compute:
     - Spearman ρ on logFC across all overlapping genes (nominal).
     - Spearman ρ on logFC restricted to genes with padj<0.05 in BOTH.
     - Sign concordance %: fraction of overlapping genes where signs match.
     - Jaccard at padj<0.05: |intersect|/|union| of significant gene sets.
     - n_overlap (total genes after join).

Environment:
  micromamba activate rnaseq   (or `spatial`; both ship pandas+scipy)

SLURM:
  --partition=io --qos=interactive --cpus-per-task=4 --mem=16G --time=2:00:00
"""

from __future__ import annotations

import os
import pathlib
import sys
import time
from typing import Dict, Tuple

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = pathlib.Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
GOVAERE_DIR = PROJECT_ROOT / "Analysis" / "Spatial" / "results" / "govaere2026"
BULK_DREAM_CSV = (
    PROJECT_ROOT
    / "RNA-seq" / "Human" / "Patient_Cohorts" / "analysis" / "integration"
    / "results" / "integration" / "dream_results_ashr.csv"
)

# Per-contrast inputs and output names
CONTRASTS: Dict[str, Dict[str, str]] = {
    "geomx_sh_vs_pt": {
        "input": "geomx_de_sh_vs_pt.csv",
        "platform": "GeoMx",
        "description": "CD68+ SH vs PT (lmer)",
        "out": "bulk_validation_geomx_sh_vs_pt.csv",
    },
    "geomx_sh_vs_ls": {
        "input": "geomx_de_sh_vs_ls.csv",
        "platform": "GeoMx",
        "description": "CD68+ SH vs LS (lmer)",
        "out": "bulk_validation_geomx_sh_vs_ls.csv",
    },
    "geomx_panck_vs_cd68": {
        "input": "geomx_de_panck_vs_cd68.csv",
        "platform": "GeoMx",
        "description": "panCK vs CD68 (segment identity)",
        "out": "bulk_validation_geomx_panck_vs_cd68.csv",
    },
    "geomx_cd45_vs_cd68": {
        "input": "geomx_de_cd45_vs_cd68.csv",
        "platform": "GeoMx",
        "description": "CD45 vs CD68 (segment identity)",
        "out": "bulk_validation_geomx_cd45_vs_cd68.csv",
    },
    "cosmx_hepatocyte": {
        "input": "cosmx_de_Hepatocyte_MASH_vs_noMASH.csv",
        "platform": "CosMx",
        "description": "Hepatocyte MASH vs no_MASH (wilcoxon)",
        "out": "bulk_validation_cosmx_hepatocyte.csv",
    },
    "cosmx_kc": {
        "input": "cosmx_de_KC_MASH_vs_noMASH.csv",
        "platform": "CosMx",
        "description": "KC (macrophage incl. MetMac) MASH vs no_MASH (wilcoxon)",
        "out": "bulk_validation_cosmx_kc.csv",
    },
}

SIG_THRESH = 0.05  # padj cutoff used for both Govaere and bulk dream

# ---------------------------------------------------------------------------
# I/O helpers
# ---------------------------------------------------------------------------


def _log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def load_bulk_dream() -> pd.DataFrame:
    """Return dream DEG table with columns: gene_symbol, bulk_dream_logFC,
    bulk_dream_padj. Collapses duplicate symbols by retaining the row with the
    smallest padj (most significant)."""
    df = pd.read_csv(BULK_DREAM_CSV)
    needed = {"symbol", "logFC", "padj"}
    missing = needed - set(df.columns)
    if missing:
        raise ValueError(f"bulk dream missing columns: {missing}")
    df = df.rename(
        columns={
            "symbol": "gene_symbol",
            "logFC": "bulk_dream_logFC",
            "padj": "bulk_dream_padj",
        }
    )[["gene_symbol", "bulk_dream_logFC", "bulk_dream_padj"]]
    df = df.dropna(subset=["gene_symbol"])
    df["gene_symbol"] = df["gene_symbol"].astype(str)
    # Collapse duplicate symbols (rare) keeping most significant
    df = df.sort_values("bulk_dream_padj", na_position="last").drop_duplicates(
        subset=["gene_symbol"], keep="first"
    )
    return df


def load_govaere_de(contrast_key: str) -> pd.DataFrame:
    """Standardize Govaere DE CSV to (gene_symbol, govaere_logFC, govaere_padj).
    GeoMx schema: gene_symbol, logFC, ..., padj_bh
    CosMx  schema: gene, logfoldchange, pval, pval_adj, ...
    """
    meta = CONTRASTS[contrast_key]
    path = GOVAERE_DIR / meta["input"]
    df = pd.read_csv(path)
    if meta["platform"] == "GeoMx":
        df = df.rename(
            columns={
                "gene_symbol": "gene_symbol",
                "logFC": "govaere_logFC",
                "padj_bh": "govaere_padj",
            }
        )[["gene_symbol", "govaere_logFC", "govaere_padj"]]
    else:  # CosMx
        df = df.rename(
            columns={
                "gene": "gene_symbol",
                "logfoldchange": "govaere_logFC",
                "pval_adj": "govaere_padj",
            }
        )[["gene_symbol", "govaere_logFC", "govaere_padj"]]
    df = df.dropna(subset=["gene_symbol"])
    df["gene_symbol"] = df["gene_symbol"].astype(str)
    # Collapse duplicate symbols (rare; e.g. CosMx panel uses unique gene names)
    df = df.sort_values("govaere_padj", na_position="last").drop_duplicates(
        subset=["gene_symbol"], keep="first"
    )
    return df


# ---------------------------------------------------------------------------
# Concordance metrics
# ---------------------------------------------------------------------------


def _safe_spearman(x: pd.Series, y: pd.Series) -> Tuple[float, float, int]:
    mask = x.notna() & y.notna()
    if mask.sum() < 3:
        return (np.nan, np.nan, int(mask.sum()))
    rho, pval = spearmanr(x[mask], y[mask])
    return (float(rho), float(pval), int(mask.sum()))


def concordance_for_contrast(
    contrast_key: str, bulk: pd.DataFrame
) -> Tuple[pd.DataFrame, Dict[str, float]]:
    """Inner-join Govaere DE with bulk dream, compute concordance metrics."""
    meta = CONTRASTS[contrast_key]
    gov = load_govaere_de(contrast_key)
    merged = gov.merge(bulk, on="gene_symbol", how="inner").copy()

    n_overlap = len(merged)
    # Spearman on logFC (nominal, all overlapping genes)
    rho_nom, pval_nom, n_used_nom = _safe_spearman(
        merged["govaere_logFC"], merged["bulk_dream_logFC"]
    )

    # Restrict to genes significant in BOTH at padj<SIG_THRESH
    both_sig = (
        (merged["govaere_padj"] < SIG_THRESH)
        & (merged["bulk_dream_padj"] < SIG_THRESH)
    )
    rho_sig, pval_sig, n_sig_both = _safe_spearman(
        merged.loc[both_sig, "govaere_logFC"],
        merged.loc[both_sig, "bulk_dream_logFC"],
    )

    # Sign concordance (over all overlapping genes with non-NA logFC on both sides)
    mask_signable = (
        merged["govaere_logFC"].notna()
        & merged["bulk_dream_logFC"].notna()
        & (merged["govaere_logFC"] != 0)
        & (merged["bulk_dream_logFC"] != 0)
    )
    sign_match = np.sign(merged.loc[mask_signable, "govaere_logFC"]) == np.sign(
        merged.loc[mask_signable, "bulk_dream_logFC"]
    )
    merged["sign_concordant"] = pd.NA
    merged.loc[mask_signable, "sign_concordant"] = sign_match.astype(int)
    sign_concord_pct = (
        100.0 * float(sign_match.sum()) / float(mask_signable.sum())
        if mask_signable.sum() > 0
        else np.nan
    )

    # Jaccard at padj<0.05 (sig union/intersection within the inner-join universe)
    gov_sig = set(
        merged.loc[merged["govaere_padj"] < SIG_THRESH, "gene_symbol"]
    )
    bulk_sig = set(
        merged.loc[merged["bulk_dream_padj"] < SIG_THRESH, "gene_symbol"]
    )
    union = gov_sig | bulk_sig
    intersect = gov_sig & bulk_sig
    jaccard = (
        float(len(intersect)) / float(len(union)) if len(union) > 0 else np.nan
    )

    summary = {
        "contrast": contrast_key,
        "platform": meta["platform"],
        "description": meta["description"],
        "n_overlap": int(n_overlap),
        "spearman_rho_logfc": rho_nom,
        "spearman_pval_logfc": pval_nom,
        "n_for_nominal_rho": int(n_used_nom),
        "spearman_rho_logfc_padj05_in_both": rho_sig,
        "spearman_pval_logfc_padj05_in_both": pval_sig,
        "n_sig_both": int(n_sig_both),
        "sign_concordance_pct": sign_concord_pct,
        "n_signable": int(mask_signable.sum()),
        "n_govaere_sig": int(len(gov_sig)),
        "n_bulk_sig_in_overlap": int(len(bulk_sig)),
        "n_sig_intersect": int(len(intersect)),
        "n_sig_union": int(len(union)),
        "jaccard_padj05": jaccard,
    }
    return merged, summary


# ---------------------------------------------------------------------------
# Top concordant table
# ---------------------------------------------------------------------------


def top_concordant_for_contrast(
    contrast_key: str, merged: pd.DataFrame, top_n: int = 10
) -> pd.DataFrame:
    """Return top-N genes that are padj<0.05 in BOTH Govaere and bulk dream,
    ordered by combined evidence (govaere_padj × bulk_dream_padj, smallest first).
    """
    both_sig = (
        (merged["govaere_padj"] < SIG_THRESH)
        & (merged["bulk_dream_padj"] < SIG_THRESH)
    )
    sub = merged.loc[both_sig].copy()
    if len(sub) == 0:
        return pd.DataFrame()
    # Combined ordering: rank by sum of -log10(padj) on both sides (robust to 0s)
    eps = 1e-300
    sub["combined_neglog10_padj"] = -(
        np.log10(sub["govaere_padj"].clip(lower=eps))
        + np.log10(sub["bulk_dream_padj"].clip(lower=eps))
    )
    sub = sub.sort_values("combined_neglog10_padj", ascending=False).head(top_n)
    sub.insert(0, "contrast", contrast_key)
    return sub[
        [
            "contrast",
            "gene_symbol",
            "govaere_logFC",
            "govaere_padj",
            "bulk_dream_logFC",
            "bulk_dream_padj",
            "sign_concordant",
            "combined_neglog10_padj",
        ]
    ]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> int:
    GOVAERE_DIR.mkdir(parents=True, exist_ok=True)

    _log(f"Loading bulk dream DEG file: {BULK_DREAM_CSV}")
    bulk = load_bulk_dream()
    _log(
        f"  bulk dream rows: {len(bulk):,} unique gene symbols; "
        f"sig padj<{SIG_THRESH}: "
        f"{int((bulk['bulk_dream_padj'] < SIG_THRESH).sum()):,}"
    )

    summaries = []
    top_tables = []
    for key in CONTRASTS:
        meta = CONTRASTS[key]
        _log(f"--- Contrast: {key} ({meta['description']}) ---")
        merged, summary = concordance_for_contrast(key, bulk)
        out_csv = GOVAERE_DIR / meta["out"]
        merged.to_csv(out_csv, index=False)
        _log(f"  wrote {out_csv.name}  (rows={len(merged):,})")
        _log(
            f"  n_overlap={summary['n_overlap']:,}  "
            f"rho_nom={summary['spearman_rho_logfc']:.3f}  "
            f"rho_sig_both={summary['spearman_rho_logfc_padj05_in_both']:.3f} "
            f"(n={summary['n_sig_both']})  "
            f"sign_concord={summary['sign_concordance_pct']:.1f}%  "
            f"jaccard={summary['jaccard_padj05']:.3f}"
        )
        summaries.append(summary)
        top_tables.append(top_concordant_for_contrast(key, merged, top_n=10))

    # Write summary TSV
    summary_df = pd.DataFrame(summaries)
    summary_path = GOVAERE_DIR / "bulk_validation_summary.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)
    _log(f"Wrote {summary_path.name}")

    # Write top-10 concordant per contrast
    top_concat = (
        pd.concat([t for t in top_tables if len(t) > 0], ignore_index=True)
        if any(len(t) > 0 for t in top_tables)
        else pd.DataFrame()
    )
    top_path = GOVAERE_DIR / "bulk_validation_top10_concordant.tsv"
    top_concat.to_csv(top_path, sep="\t", index=False)
    _log(f"Wrote {top_path.name}  (rows={len(top_concat):,})")

    # Final inline report
    _log("===== SUMMARY =====")
    cols_show = [
        "contrast",
        "n_overlap",
        "spearman_rho_logfc",
        "spearman_rho_logfc_padj05_in_both",
        "n_sig_both",
        "sign_concordance_pct",
        "jaccard_padj05",
    ]
    with pd.option_context("display.max_columns", None, "display.width", 200):
        print(summary_df[cols_show].to_string(index=False), flush=True)

    _log("DONE")
    return 0


if __name__ == "__main__":
    sys.exit(main())
