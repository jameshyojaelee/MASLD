#!/usr/bin/env python3
"""
Generate the ancestry-COLOC web-data pair for the MASLD atlas portal (data
contract §5.3):

  - coloc_ancestry_summary.json  (per-ancestry counts, confidence tiers,
                                  cross-ancestry-replicated genes, headline
                                  effector-gene counts, Broadaway overlap)
  - coloc_by_gwas.parquet        (long-format, one row per gene x GWAS x method
                                  colocalization; PP.H4.abf > 0.01; sorted by
                                  `symbol`)

This single parquet is consumed BOTH by the gene page (WHERE symbol=?) and by
/genetics (per contract §8 item 5, it supersedes the separate
`gene_coloc_by_gwas.parquet` of §3.3 — materialized once, here).

Headline effector-gene counts are the MAIN Tier-1/2 liver-specific portfolio
(35 GWAS; NAFLD/NASH/PDFF + ALT/AST/GGT), read from the pre-tiered
`gene_level_coloc_tier12.csv`:
    SuSiE (PP.H4.susie > 0.5)                = 473
    union (SuSiE ∪ ABF, PP.H4 > 0.5)         = 1031
The full 50-GWAS portfolio (incl. Tier-3/4 cirrhosis/ChronLiver/albumin/
platelet) drives the per-ancestry breadth, confidence tiers and the parquet.

Env: micromamba run -n spatial (NumPy2/pandas ABI). Honors MASLD_PROJECT_ROOT
and --output-dir. Compact JSON, NaN/Inf -> None.

Usage:
  micromamba run -n spatial python scripts/portal/generate_ancestry_coloc_data.py \
      [--output-dir masld-atlas-v2/public/data]
"""

from __future__ import annotations

import argparse
import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
SUSIE_DIR = PROJECT_ROOT / "GWAS/finemapping/results/susie_coloc"
GENE_LEVEL = SUSIE_DIR / "gene_level_coloc.csv"              # full 50-GWAS aggregate
GENE_LEVEL_T12 = SUSIE_DIR / "gene_level_coloc_tier12.csv"   # Tier-1/2 headline set
ALL_GWAS = SUSIE_DIR / "susie_coloc_all_gwas.csv"            # long, 927k rows
BROADAWAY = SUSIE_DIR / "broadaway_coloc_comparison.csv"

# Fixed ancestry order (contract §5.3).
ANCESTRIES = ["EUR", "AFR", "EAS", "AMR", "SAS"]

# Long-file columns we actually need (keep the read light).
LONG_USECOLS = [
    "gwas_name",
    "gene",
    "ensembl",
    "PP.H4.abf",
    "PP.H4.susie",
    "method",
    "ancestry",
    "ld_panel",
    "ld_reliability",
]


def _is_true(x) -> bool:
    """R-style 'TRUE'/'FALSE' (or real booleans) -> Python bool."""
    if isinstance(x, bool):
        return x
    return str(x).strip().upper() == "TRUE"


def _num_or_none(x):
    """Float, mapping NaN/Inf -> None (valid JSON)."""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(v):
        return None
    return v


def _round_or_none(x, ndp):
    v = _num_or_none(x)
    return None if v is None else round(v, ndp)


# ---------------------------------------------------------------------------
# Summary pieces
# ---------------------------------------------------------------------------
def headline_effector_counts() -> dict:
    """Tier-1/2 (main) headline effector-gene counts from the pre-tiered file."""
    t12 = pd.read_csv(GENE_LEVEL_T12)
    susie = t12["coloc_best_susie_pp4"].astype(float) > 0.5
    abf = t12["coloc_best_abf_pp4"].astype(float) > 0.5
    return {
        "n_effector_genes_susie": int(susie.sum()),
        "n_effector_genes_union": int((susie | abf).sum()),
    }


def confidence_tiers() -> list[dict]:
    """SuSiE confidence-tier counts (canonical SuSiE-primary column) from the
    full aggregate. Excludes 'none'."""
    gl = pd.read_csv(GENE_LEVEL, usecols=["coloc_susie_conf_tier"])
    counts = gl["coloc_susie_conf_tier"].astype(str).str.lower().value_counts()
    out = []
    for raw, label in [("high", "High"), ("suggestive", "Suggestive"), ("nominal", "Nominal")]:
        out.append({"tier": label, "n": int(counts.get(raw, 0))})
    return out


def broadaway_overlap() -> dict:
    bd = pd.read_csv(BROADAWAY)
    return {
        "n_in_broadaway_747": int(bd["in_broadaway_747"].map(_is_true).sum()),
        "n_novel_to_us": int(bd["novel_to_us"].map(_is_true).sum()),
    }


def cross_ancestry_replicated(gene_ancestries: dict[str, list[str]]) -> list[dict]:
    """Genes flagged colocalizing in >=2 ancestries. best_pp4 = max(ABF, SuSiE)
    from the aggregate; ancestries = fixed-order list from the long-file
    membership map (best row-PP4 > 0.5)."""
    gl = pd.read_csv(
        GENE_LEVEL,
        usecols=[
            "gene",
            "coloc_best_pp4",
            "coloc_best_susie_pp4",
            "coloc_cross_ancestry_replicated",
            "coloc_best_ancestry",
        ],
    )
    gl = gl[gl["coloc_cross_ancestry_replicated"].map(_is_true)].copy()
    rows = []
    for _, r in gl.iterrows():
        sym = r["gene"]
        if not isinstance(sym, str) or not sym or sym.lower() == "nan":
            continue
        best = max(
            _num_or_none(r["coloc_best_pp4"]) or float("-inf"),
            _num_or_none(r["coloc_best_susie_pp4"]) or float("-inf"),
        )
        best = None if not math.isfinite(best) else round(best, 4)
        ancs = gene_ancestries.get(sym, [])
        if not ancs and isinstance(r["coloc_best_ancestry"], str):
            ancs = [r["coloc_best_ancestry"]]
        rows.append({"gene": sym, "best_pp4": best, "ancestries": ancs})
    rows.sort(key=lambda d: (-(d["best_pp4"] or 0.0), d["gene"]))
    return rows


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--output-dir",
        default="masld-atlas-v2/public/data",
        help="Output dir (relative to MASLD_PROJECT_ROOT unless absolute).",
    )
    args = ap.parse_args()

    out_dir = Path(args.output_dir)
    if not out_dir.is_absolute():
        out_dir = PROJECT_ROOT / out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    for p in (GENE_LEVEL, GENE_LEVEL_T12, ALL_GWAS, BROADAWAY):
        if not p.exists():
            raise SystemExit(f"Missing source: {p}")

    # ---- long file (full portfolio) -------------------------------------
    print(f"Reading {ALL_GWAS.name} (usecols={len(LONG_USECOLS)}) ...")
    long = pd.read_csv(
        ALL_GWAS,
        usecols=LONG_USECOLS,
        dtype={
            "gwas_name": "string",
            "gene": "string",
            "ensembl": "string",
            "method": "string",
            "ancestry": "string",
            "ld_panel": "string",
            "ld_reliability": "string",
        },
    )
    long["PP.H4.abf"] = pd.to_numeric(long["PP.H4.abf"], errors="coerce")
    long["PP.H4.susie"] = pd.to_numeric(long["PP.H4.susie"], errors="coerce")
    print(f"  {len(long):,} rows; {long['gwas_name'].nunique()} GWAS; "
          f"ancestries={sorted(long['ancestry'].dropna().unique())}")

    # Row-level best PP4 (union of ABF/SuSiE) for per-ancestry gene counts.
    long["_best"] = long[["PP.H4.abf", "PP.H4.susie"]].max(axis=1, skipna=True)

    # by_ancestry: distinct genes >0.5 / >0.8 + #GWAS per ancestry.
    by_anc = []
    for a in ANCESTRIES:
        sub = long[long["ancestry"] == a]
        if sub.empty:
            continue
        by_anc.append(
            {
                "ancestry": a,
                "n_gwas": int(sub["gwas_name"].nunique()),
                "n_genes_h4_05": int(sub.loc[sub["_best"] > 0.5, "gene"].nunique()),
                "n_genes_h4_08": int(sub.loc[sub["_best"] > 0.8, "gene"].nunique()),
            }
        )

    # per-gene ancestry membership (best row-PP4 > 0.5), fixed ancestry order.
    hit = long[(long["_best"] > 0.5) & long["gene"].notna()]
    anc_by_gene_raw = hit.groupby("gene")["ancestry"].apply(lambda s: set(s.dropna()))
    gene_ancestries = {
        g: [a for a in ANCESTRIES if a in s] for g, s in anc_by_gene_raw.items()
    }

    # ---- parquet (filtered PP.H4.abf > 0.01) ----------------------------
    filt = long[long["PP.H4.abf"] > 0.01].copy()
    parquet = pd.DataFrame(
        {
            "symbol": filt["gene"].astype("string").fillna(""),
            "ensembl": filt["ensembl"].astype("string").fillna(""),
            "gwas": filt["gwas_name"].astype("string"),
            "ancestry": filt["ancestry"].astype("string"),
            "method": filt["method"].astype("string"),
            "pp4": filt["PP.H4.abf"].round(4),
            "pp4_susie": filt["PP.H4.susie"].round(4),
            "ld_panel": filt["ld_panel"].astype("string"),
            "ld_reliability": filt["ld_reliability"].astype("string"),
        }
    )
    # ±Inf -> NaN before write (contract §0.4); sort by symbol (row-group skip).
    parquet = parquet.replace([np.inf, -np.inf], np.nan)
    parquet = parquet.sort_values("symbol", kind="stable").reset_index(drop=True)
    pq_path = out_dir / "coloc_by_gwas.parquet"
    parquet.to_parquet(pq_path, engine="pyarrow", index=False)
    print(f"  -> {pq_path.name}: {len(parquet):,} rows, {parquet.shape[1]} cols")

    # ---- summary json ---------------------------------------------------
    summary = {
        **headline_effector_counts(),
        "n_gwas": int(long["gwas_name"].nunique()),
        "ancestries": ANCESTRIES,
        "by_ancestry": by_anc,
        "conf_tier_counts": confidence_tiers(),
        "cross_ancestry_replicated": cross_ancestry_replicated(gene_ancestries),
        "broadaway_overlap": broadaway_overlap(),
    }

    json_path = out_dir / "coloc_ancestry_summary.json"
    with open(json_path, "w") as fh:
        json.dump(summary, fh, separators=(",", ":"))
    print(f"  -> {json_path.name}: {json_path.stat().st_size/1024:.1f} KB")

    # ---- console verification ------------------------------------------
    print("\n=== coloc_ancestry_summary.json ===")
    print(f"  n_effector_genes_susie = {summary['n_effector_genes_susie']} (target 473)")
    print(f"  n_effector_genes_union = {summary['n_effector_genes_union']} (target 1031)")
    print(f"  n_gwas                 = {summary['n_gwas']}")
    print(f"  by_ancestry            = {summary['by_ancestry']}")
    print(f"  conf_tier_counts       = {summary['conf_tier_counts']}")
    print(f"  cross_ancestry_replicated = {len(summary['cross_ancestry_replicated'])} genes")
    print(f"  broadaway_overlap      = {summary['broadaway_overlap']}")
    print("\n=== coloc_by_gwas.parquet ===")
    print(f"  cols     = {list(parquet.columns)}")
    print(f"  rows     = {len(parquet):,}")
    print(f"  ancestry = {sorted(parquet['ancestry'].dropna().unique())}")
    print(f"  n_gwas   = {parquet['gwas'].nunique()}")


if __name__ == "__main__":
    main()
