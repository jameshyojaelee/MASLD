#!/usr/bin/env python
"""Post-hoc S8 layer — decoupled combined ranking on top of 46d output.

Per plan §16 + user decision 2026-05-21 Option C: 46d (convergence-evidence
ranking) remains UNCHANGED — it's a sophisticated, evolving Wakefield-ABF +
Brown's correction + INTACT framework that is unsafe to monkey-patch with a
naive additive S8 term.

This script consumes 46d output (convergence_evidence.csv) AND the atlas
S8 columns (perturb_validation_tier, perturb_n_arms_passing, etc. — added by
Script 27a when run with INCLUDE_S8=TRUE) and produces a side-by-side
combined ranking — WITHOUT blending the numerical scores into one.

Output schema (combined_ranking.csv):
    human_symbol, ensembl_id,
    # 46d axis (observational)
    convergence_score, convergence_tier, concordance_state,
    convergence_rank,
    # S8 axis (interventional / counterfactual)
    perturb_validation_tier, perturb_n_arms_applicable, perturb_n_arms_passing,
    perturb_rank,
    # combined decision label (4 cases — NOT a blended score)
    decision_label ∈ {convergent_only, perturbation_only, both, neither}
    avg_percentile  # auxiliary, lightweight: mean of percentile-ranks across the 2 axes
                    # presented as supplementary, NOT canonical

Reviewer-defensible framing: S8 is interventional evidence that's categorically
different from S1–S7 observational sources. Combining them via Bayes-factor
math would require calibrating S8 against held-out experiments (which we don't
have for in-silico perturbation predictions). Until that's done, the safe and
honest framing is: "high-confidence targets are those where convergent
observational evidence AND counterfactual perturbation prediction both agree."

Usage:
    # default paths
    python post_46d_s8_layer.py

    # explicit paths
    python post_46d_s8_layer.py \
        --convergence-csv RNA-seq/results/multi_evidence/convergence_evidence.csv \
        --atlas-csv RNA-seq/results/multi_evidence/multi_evidence_atlas.csv \
        --out-csv RNA-seq/results/multi_evidence/combined_ranking.csv

Environment: rnaseq (any env with pandas + numpy works).
"""
from __future__ import annotations

import argparse
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
ME = PROJECT_ROOT / "RNA-seq/results/multi_evidence"

DEFAULT_CONVERGENCE = ME / "convergence_evidence.csv"
DEFAULT_ATLAS = ME / "multi_evidence_atlas.csv"
DEFAULT_OUT = ME / "combined_ranking.csv"

# Convergence-tier names from 46d (per docs/refactor/2026-05-19-bayes-to-convergence-rename.md)
CONVERGENCE_TIER_ORDER = ("1_Genetic_validated", "2_Strong", "3_Moderate", "4_Weak")
PERTURB_TIER_ORDER = ("strong", "validated", "partial", "—")


def _parse_args():
    p = argparse.ArgumentParser(__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--convergence-csv", default=str(DEFAULT_CONVERGENCE))
    p.add_argument("--atlas-csv", default=str(DEFAULT_ATLAS))
    p.add_argument("--out-csv", default=str(DEFAULT_OUT))
    p.add_argument(
        "--no-perturb-required",
        action="store_true",
        help="If set, do NOT require atlas to have perturb_* columns. "
        "Run on convergence alone with empty S8 columns.",
    )
    return p.parse_args()


def _classify(row) -> str:
    """Return one of {convergent_only, perturbation_only, both, neither}.

    'High' = convergence_tier in {1_Genetic_validated, 2_Strong} OR
            perturb_validation_tier in {strong, validated}.
    """
    conv_high = str(row.get("convergence_tier", "")) in (
        "1_Genetic_validated",
        "2_Strong",
    )
    perturb_tier = str(row.get("perturb_validation_tier", ""))
    perturb_high = perturb_tier in ("strong", "validated")
    if conv_high and perturb_high:
        return "both"
    if conv_high and not perturb_high:
        return "convergent_only"
    if perturb_high and not conv_high:
        return "perturbation_only"
    return "neither"


def main() -> None:
    args = _parse_args()
    convergence_path = Path(args.convergence_csv)
    atlas_path = Path(args.atlas_csv)
    out_path = Path(args.out_csv)

    if not convergence_path.exists():
        raise FileNotFoundError(
            f"Convergence file not found at {convergence_path}. "
            f"Run RNA-seq/46d_convergence_evidence.R first."
        )

    print(f"[post_46d_s8] loading convergence: {convergence_path.name}")
    conv = pd.read_csv(convergence_path)
    print(f"  shape: {conv.shape}")
    print(f"  tiers: {dict(conv['convergence_tier'].value_counts(dropna=False)) if 'convergence_tier' in conv else 'no convergence_tier col'}")

    # Atlas read — defensive; the S8 columns may not exist yet
    atlas_cols_needed = [
        "human_symbol",
        "ensembl_id",
        "perturb_validation_tier",
        "perturb_n_arms_applicable",
        "perturb_n_arms_passing",
    ]
    if not atlas_path.exists():
        if not args.no_perturb_required:
            raise FileNotFoundError(
                f"Atlas file not found at {atlas_path}. "
                f"Pass --no-perturb-required to run on convergence alone."
            )
        print(f"[post_46d_s8] atlas missing; running with empty S8 axis")
        atlas = pd.DataFrame(columns=atlas_cols_needed)
    else:
        print(f"[post_46d_s8] loading atlas: {atlas_path.name}")
        atlas = pd.read_csv(atlas_path, usecols=lambda c: c in atlas_cols_needed)
        missing_perturb = [c for c in atlas_cols_needed if c not in atlas.columns]
        if missing_perturb:
            print(f"  WARN atlas missing cols {missing_perturb} (Script 27a INCLUDE_S8 not yet run)")
            for c in missing_perturb:
                atlas[c] = pd.NA

    # Join on human_symbol (canonical atlas key)
    if "human_symbol" not in conv.columns:
        # 46d may use `gene` instead
        if "gene" in conv.columns:
            conv = conv.rename(columns={"gene": "human_symbol"})
        else:
            raise ValueError(f"Convergence has no human_symbol/gene column. Cols: {list(conv.columns[:10])}")

    merged = conv.merge(atlas, on="human_symbol", how="left", suffixes=("", "__atlas"))

    # Compute ranks per axis (1 = best)
    if "convergence_score" in merged.columns:
        merged["convergence_rank"] = merged["convergence_score"].rank(ascending=False, method="min").astype("Int64")
    else:
        merged["convergence_rank"] = pd.NA

    # Perturb rank: derive from n_arms_passing (higher = better); tier-stratified
    if "perturb_n_arms_passing" in merged.columns:
        # Sort by tier then n_arms_passing then convergence_rank as tiebreaker
        tier_order_map = {t: i for i, t in enumerate(PERTURB_TIER_ORDER)}
        merged["_tier_ord"] = merged["perturb_validation_tier"].map(tier_order_map).fillna(99)
        merged["_pass"] = merged["perturb_n_arms_passing"].fillna(-1)
        merged["perturb_rank"] = (
            merged.sort_values(["_tier_ord", "_pass"], ascending=[True, False])
                  .reset_index()
                  .reset_index()
                  .set_index("index")["level_0"] + 1
        )
        merged["perturb_rank"] = merged["perturb_rank"].astype("Int64")
        merged = merged.drop(columns=["_tier_ord", "_pass"])
    else:
        merged["perturb_rank"] = pd.NA

    # Decision label (4 cases — NOT a blended score)
    merged["decision_label"] = merged.apply(_classify, axis=1)

    # Auxiliary blended percentile (presented as supplementary)
    n = len(merged)
    conv_pct = (merged["convergence_rank"].fillna(n + 1) / n).clip(upper=1.0)
    perturb_pct = (merged["perturb_rank"].fillna(n + 1) / n).clip(upper=1.0)
    merged["avg_percentile"] = ((conv_pct + perturb_pct) / 2).round(4)

    # Reorder columns
    out_cols_priority = [
        "human_symbol",
        "ensembl_id",
        "convergence_score",
        "convergence_tier",
        "concordance_state",
        "convergence_rank",
        "perturb_validation_tier",
        "perturb_n_arms_applicable",
        "perturb_n_arms_passing",
        "perturb_rank",
        "decision_label",
        "avg_percentile",
    ]
    keep_cols = [c for c in out_cols_priority if c in merged.columns]
    extras = [c for c in merged.columns if c not in keep_cols and not c.endswith("__atlas")]
    out_df = merged[keep_cols + extras]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(out_path, index=False)
    print(f"[post_46d_s8] wrote {out_path} ({len(out_df):,} rows × {out_df.shape[1]} cols)")

    # Audit
    print(f"\n=== Audit ===")
    print(f"Decision-label distribution:")
    print(out_df["decision_label"].value_counts())
    if "perturb_validation_tier" in out_df.columns:
        print(f"\nPerturbation tier distribution:")
        print(out_df["perturb_validation_tier"].value_counts(dropna=False))
    print(f"\nConvergence tier × Perturb tier cross-tab (top 6 categories):")
    if "convergence_tier" in out_df.columns and "perturb_validation_tier" in out_df.columns:
        ct = pd.crosstab(
            out_df["convergence_tier"].fillna("—"),
            out_df["perturb_validation_tier"].fillna("—"),
        )
        print(ct)


if __name__ == "__main__":
    main()
