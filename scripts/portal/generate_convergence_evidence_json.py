#!/usr/bin/env python3
"""
Generate convergence_ranking.json for the portal from
RNA-seq/results/multi_evidence/convergence_evidence.csv (output of 46d).

NOTE: The score is an evidence-weighted convergence ranking (NOT a Bayesian
posterior — no prior or likelihood is specified). Renamed from
bayesian_evidence_ranking / bayesian_ranking on 2026-05-19.
See docs/refactor/2026-05-19-bayes-to-convergence-rename.md.

Usage: python scripts/portal/generate_convergence_evidence_json.py
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
EV_CSV = PROJECT_ROOT / "RNA-seq/results/multi_evidence/convergence_evidence.csv"
# Canonical convergence ranking. Replaces the legacy bayesian_ranking.json
# (retired Script 46b) and the previous bayesian_evidence_ranking.json (46d).
OUT = PROJECT_ROOT / "masld-atlas-v2/public/data/convergence_ranking.json"


def main() -> None:
    if not EV_CSV.exists():
        raise SystemExit(f"Missing {EV_CSV}. Run 46d first.")

    print(f"Reading {EV_CSV}")
    # keep_default_na=False ensures "None" string (from R fifelse) isn't coerced
    # to NaN — preserves the literal druggability tier value.
    ev = pd.read_csv(
        EV_CSV,
        keep_default_na=False,
        na_values=["", "NA"],
        low_memory=False,
    )
    # Exclude the confounder-stripped rows and those without rank
    ev = ev[~ev["excluded_from_ranking"].astype(bool) & ev["convergence_rank"].notna()].copy()
    ev["convergence_rank"] = ev["convergence_rank"].astype(int)
    ev = ev.sort_values("convergence_rank")

    n = len(ev)
    print(f"  {n} genes in rank table")

    # Optional atlas-derived convenience flags (is_deg, is_coloc, is_druggable)
    # to match convergence_ranking.json's schema for the portal UI.
    is_deg = ev["log_BF_S1"].abs() > 1.1  # log_BF_S1 > log(3) ≈ active DE evidence
    is_coloc = ev["coloc_best_pp4_S2a"].fillna(0.0) >= 0.5
    is_druggable = ev["druggability_tier"].fillna("None") != "None"
    pct = 100.0 * (1.0 - (ev["convergence_rank"].values - 1) / max(n, 1))

    # 2026-05-21 Govaere2026 signature extension columns (optional — gracefully
    # absent for older runs).
    has_govaere = "signature_govaere2026_n_evidence" in ev.columns
    if has_govaere:
        print("  including Govaere2026 signature columns in portal JSON")

    genes = []
    for i, row in ev.iterrows():
        symbol = row["human_symbol"]
        if not isinstance(symbol, str) or not symbol or symbol.lower() == "nan":
            continue
        gene_entry = {
            "rank": int(row["convergence_rank"]),
            "symbol": symbol,
            "convergence_score": round(float(row["convergence_score"]), 5),
            "log_convergence_odds": round(float(row["log_convergence_odds"]), 3),
            "evidence_total": round(float(row["evidence_total"]), 3),
            "tier": row["tier"],
            "concordance_state": row["concordance_state"],
            "dominant_stage_S1": row["dominant_stage_S1"]
            if pd.notna(row["dominant_stage_S1"])
            else None,
            # cross_species_concordant field REMOVED 2026-07-05 — no cross-species in atlas paper
            "druggability_tier": row["druggability_tier"],
            "coloc_best_gwas": row["coloc_best_gwas_S2a"]
            if pd.notna(row["coloc_best_gwas_S2a"])
            else None,
            "coloc_best_pp4": round(float(row["coloc_best_pp4_S2a"]), 3)
            if pd.notna(row["coloc_best_pp4_S2a"])
            else None,
            "n_modalities_active": int(row["n_modalities_active"]),
            "percentile": round(float(pct[list(ev.index).index(i)]), 2),
            "is_deg": bool(is_deg.loc[i]),
            "is_coloc": bool(is_coloc.loc[i]),
            "is_druggable": bool(is_druggable.loc[i]),
        }
        if has_govaere:
            gene_entry["govaere2026_evidence"] = bool(
                row.get("signature_govaere2026_evidence", False)
            )
            gene_entry["govaere2026_n_evidence"] = int(
                row.get("signature_govaere2026_n_evidence", 0) or 0
            )
            gv_score = row.get("convergence_score_with_govaere2026")
            gv_rank = row.get("convergence_rank_with_govaere2026")
            gene_entry["convergence_score_with_govaere2026"] = (
                round(float(gv_score), 5) if pd.notna(gv_score) else None
            )
            gene_entry["convergence_rank_with_govaere2026"] = (
                int(gv_rank) if pd.notna(gv_rank) else None
            )
        genes.append(gene_entry)

    payload = {
        "genes": genes,
        "metadata": {
            "n_genes": len(genes),
            "source": f"atlas file: {EV_CSV.relative_to(PROJECT_ROOT)}",
            "score_definition": (
                "convergence_score = evidence-weighted ranking across 6 active human evidence modalities. "
                "NOT a Bayesian posterior (no prior or likelihood specified). "
                "Mouse/cross-species evidence is excluded from the canonical score. "
                "Unsupervised heuristic; calibrated empirically via permutation FDR. "
                "Renamed from posterior_prob / bayesian_evidence on 2026-05-19."
            ),
            "tier_definitions": {
                "1_Genetic_validated": "PP4 > 0.7 OR |TWAS z| > 4, with ≥2 active modalities",
                "2_Convergent": "≥5 active modalities, convergence_score > 0.9, concordant direction",
                "3_Suggestive": "3-4 active modalities, convergence_score > 0.7",
                "4_Weak": "background",
            },
            "concordance_states": [
                "Concordant-up",
                "Concordant-down",
                "Protective-LOF",
                "Conflicted",
                "Insufficient-evidence",
            ],
        },
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(payload, fh, separators=(",", ":"))
    size_kb = OUT.stat().st_size / 1024
    print(f"  -> {OUT.relative_to(PROJECT_ROOT)}: {size_kb:.1f} KB, {len(genes)} genes")


if __name__ == "__main__":
    main()
