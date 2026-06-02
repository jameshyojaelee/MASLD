#!/usr/bin/env python3
"""Generate proteomics_summary.json for the portal /proteomics page.

Inputs:
  - Analysis/Proteomics/results/protein_differential_results_v3.csv
  - Analysis/Proteomics/results/mrna_protein_concordance_stratified.csv
  - Analysis/Proteomics/results/protein_validated_targets.csv
  - RNA-seq/results/multi_evidence/multi_evidence_atlas.csv  (for proteomics-significant count)

Output:
  masld-atlas-v2/public/data/proteomics_summary.json

Run:
  micromamba run -n spatial python scripts/portal/generate_proteomics_data.py
"""

from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PROT_DE = PROJECT_ROOT / "Analysis/Proteomics/results/protein_differential_results_v3.csv"
CONCORDANCE = PROJECT_ROOT / "Analysis/Proteomics/results/mrna_protein_concordance_stratified.csv"
TARGETS = PROJECT_ROOT / "Analysis/Proteomics/results/protein_validated_targets.csv"
ATLAS = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
OUT = PROJECT_ROOT / "masld-atlas-v2/public/data/proteomics_summary.json"


def _r(v, n=4):
    if v is None:
        return None
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    if math.isnan(v) or math.isinf(v):
        return None
    return round(v, n)


def main():
    print(f"Loading {PROT_DE}")
    de = pd.read_csv(PROT_DE)
    de = de[de["gene"].notna() & de["padj"].notna()].copy()
    # Best-per-gene (smallest padj)
    de = de.sort_values("padj").drop_duplicates("gene", keep="first")
    n_proteins = len(de)
    print(f"  {n_proteins} unique proteins")

    # Top up/down by padj, logFC sign
    sig = de[de["padj"] < 0.05].copy()
    up = sig[sig["logFC"] > 0].sort_values("padj").head(50)
    down = sig[sig["logFC"] < 0].sort_values("padj").head(50)

    def _row(r):
        return {
            "gene": str(r["gene"]),
            "logfc": _r(r.get("logFC"), 3),
            "padj": _r(r.get("padj"), 6),
            "tstat": _r(r.get("t"), 2),
            "dataset": str(r.get("dataset") or ""),
        }

    top_up = [_row(r) for _, r in up.iterrows()]
    top_down = [_row(r) for _, r in down.iterrows()]

    # Concordance stats
    print(f"Loading {CONCORDANCE}")
    conc = pd.read_csv(CONCORDANCE)
    conc_rows = []
    for _, r in conc.iterrows():
        conc_rows.append({
            "stratum": str(r["stratum"]),
            "n": int(r["n"]),
            "rho": _r(r["rho"], 3),
            "direction_pct": _r(r["direction_pct"], 1),
        })

    overall_rho = None
    core_rho = None
    for row in conc_rows:
        if row["stratum"] == "All genes":
            overall_rho = row["rho"]
        if row["stratum"] == "Conserved":  # stratum key emitted by mrna_protein_concordance.R
            core_rho = row["rho"]

    # Validated targets
    print(f"Loading {TARGETS}")
    tgt = pd.read_csv(TARGETS)
    _detected = tgt["detected"].astype(str).str.upper().eq("TRUE")
    # count distinct genes, not gene×dataset rows (a gene detected in >1 dataset must not double-count)
    n_validated = int(tgt.loc[_detected, "gene"].nunique()) if "gene" in tgt.columns else int(_detected.sum())

    # Significant proteins from atlas column
    print(f"Loading atlas for best_protein_padj count")
    atlas = pd.read_csv(ATLAS, usecols=["human_symbol", "best_protein_logFC", "best_protein_padj"], low_memory=False)
    n_tested_atlas = int(atlas["best_protein_padj"].notna().sum())
    n_sig_atlas = int((pd.to_numeric(atlas["best_protein_padj"], errors="coerce") < 0.05).sum())

    summary = {
        "n_proteins_tested": n_proteins,
        "n_proteins_significant": int((de["padj"] < 0.05).sum()),
        "n_atlas_tested": n_tested_atlas,
        "n_atlas_significant": n_sig_atlas,
        "overall_rho": overall_rho,
        "conserved_rho": core_rho,
        "n_validated_targets": n_validated,
        "concordance_strata": conc_rows,
        "top_upregulated": top_up,
        "top_downregulated": top_down,
        "datasets": sorted({str(d) for d in de["dataset"].dropna().unique()}),
        "provenance": {
            "pxd052937": "PRIDE PXD052937 (DIA-MS plasma, 72 samples)",
            "pxd051911": "PRIDE PXD051911 (DIA-MS liver tissue, 58 samples)",
            "gse276114": "GSE276114 (liver fibrosis cohort)",
            "olink": "Olink Explore 1536 plasma panel (Yang et al. 2025; supervised results withdrawn — see review 2026-06-01)",
        },
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(summary, fh, separators=(",", ":"))
    size_kb = OUT.stat().st_size / 1024
    print(f"  -> {OUT.name}: {size_kb:.1f} KB")
    print(f"  n_proteins_tested={n_proteins}, atlas_sig={n_sig_atlas}, overall_rho={overall_rho}, core_rho={core_rho}, validated={n_validated}")


if __name__ == "__main__":
    main()
