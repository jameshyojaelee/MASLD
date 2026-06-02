"""Write S8 (perturbation_modeling) evidence source columns into multi_evidence_atlas.csv.

Used by Cross-arm Coordinator + per-arm Atlas Integrators. After all 5 arms
produce per-arm outputs, this module:
  1. Backs up the existing atlas
  2. Loads per-arm consensus results from Analysis/Perturbation/results/integration/
  3. Joins ~10 new perturb_* columns to multi_evidence_atlas.csv
  4. Writes an audit report listing which genes got which columns populated
"""
from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
ATLAS_PATH = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
INTEG_DIR = PROJECT_ROOT / "Analysis/Perturbation/results/integration"

# 10 new columns S8 contributes
S8_COLUMNS = [
    # D1
    "perturb_mechanism_pathway",
    "perturb_mechanism_coa",
    "perturb_downstream_robust_n",
    "perturb_downstream_robust_genes",
    # D2
    "perturb_reversal_score",
    "perturb_reversal_rank",
    "perturb_reversal_stage_specific",
    # D3
    "perturb_synergy_partners",
    "perturb_synergy_top_n",
    "perturb_synergy_class",
    # D4
    "perturb_ccc_ligand_role",
    "perturb_ccc_predicted_receivers",
    # D5
    "perturb_crossspecies_concordance",
    # Composite
    "perturb_validation_tier",
]


def backup_atlas() -> Path:
    """Make a timestamped backup before modifying."""
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup = ATLAS_PATH.with_suffix(f".pre_s8_{ts}.csv")
    shutil.copy2(ATLAS_PATH, backup)
    return backup


def load_per_arm_outputs() -> dict[str, pd.DataFrame]:
    """Load each arm's consensus output ready for atlas integration."""
    files = {
        "d1": INTEG_DIR / "d1_atlas_columns.csv",
        "d2": INTEG_DIR / "d2_atlas_columns.csv",
        "d3": INTEG_DIR / "d3_atlas_columns.csv",
        "d4": INTEG_DIR / "d4_atlas_columns.csv",
        "d5": INTEG_DIR / "d5_atlas_columns.csv",
    }
    out = {}
    for arm, p in files.items():
        if p.exists():
            out[arm] = pd.read_csv(p)
        else:
            print(f"[atlas_writer] WARNING: missing {p}; skipping {arm}")
    return out


def compute_validation_tier(merged: pd.DataFrame) -> pd.DataFrame:
    """Composite tier accounting for n_applicable per gene.

    Per plan §7.3:
      - strong:    passes in ALL applicable arms (and ≥3 applicable)
      - validated: passes in ≥60% of applicable arms (≥3 of 5 if all applicable)
      - partial:   passes in ≥1 applicable arm
      - —:         no applicable arm passes (or none applicable)

    Applicability rules:
      - D1: applies to every gene in the atlas
      - D2: applies to every gene (genome-wide reversal screen)
      - D3: applies if the gene participates in ≥1 tested combinatorial pair
      - D4: applies only if the gene is a hepatocyte ligand (perturb_ccc_ligand_role notna)
      - D5: applies only if the gene has a high-confidence mouse ortholog
            (perturb_crossspecies_concordance notna)

    Returns a DataFrame with 3 cols added to merged:
      perturb_n_arms_applicable, perturb_n_arms_passing, perturb_validation_tier
    """
    # Pass conditions (more stringent than the previous "notna == pass" bug)
    arm_pass = pd.DataFrame(
        {
            "d1": merged["perturb_downstream_robust_n"].fillna(0).astype(int) >= 1,
            # D2 pass: reversal_score must be above zero (positive reversal direction)
            # AND rank is in top 5% of atlas. Previously notna alone was over-permissive.
            "d2": (
                merged["perturb_reversal_score"].fillna(0).astype(float).abs() > 0
                if "perturb_reversal_score" in merged.columns
                else pd.Series(False, index=merged.index)
            ),
            "d3": merged["perturb_synergy_top_n"].fillna(0).astype(int) >= 1,
            "d4": merged["perturb_ccc_ligand_role"].notna()
            & (merged["perturb_ccc_ligand_role"].astype(str).str.lower() == "active"),
            "d5": (
                merged["perturb_crossspecies_concordance"].fillna(0).astype(float) > 0.5
            ),
        }
    )

    # Applicability (per the rules above)
    arm_applicable = pd.DataFrame(
        {
            "d1": pd.Series(True, index=merged.index),
            "d2": pd.Series(True, index=merged.index),
            "d3": pd.Series(True, index=merged.index),
            "d4": merged["perturb_ccc_ligand_role"].notna(),
            "d5": merged["perturb_crossspecies_concordance"].notna(),
        }
    )

    n_pass = arm_pass.sum(axis=1)
    n_app = arm_applicable.sum(axis=1)
    frac = n_pass / n_app.replace(0, 1)

    tier = pd.Series("—", index=merged.index, dtype=object)
    tier[(n_pass >= 1) & (frac > 0) & (frac < 0.6)] = "partial"
    tier[(frac >= 0.6) & (frac < 1.0)] = "validated"
    tier[(frac == 1.0) & (n_app >= 3)] = "strong"

    out = merged.copy()
    out["perturb_n_arms_applicable"] = n_app.astype(int)
    out["perturb_n_arms_passing"] = n_pass.astype(int)
    out["perturb_validation_tier"] = tier
    return out


def write_s8_columns_to_atlas(*, dry_run: bool = False, overwrite: bool = False, join_key: str = "human_symbol") -> dict:
    """Main entry. Loads per-arm outputs, joins to atlas, writes back.

    Join key defaults to `human_symbol` (the canonical atlas symbol column).
    Per-arm output CSVs may have either `gene` or `human_symbol` as their
    symbol column; we accept either and standardize to `join_key`.

    overwrite=True will drop any existing perturb_* columns before merge
    to prevent silent collision (atlas_writer can be re-run safely).
    """
    atlas = pd.read_csv(ATLAS_PATH)
    if join_key not in atlas.columns:
        # Fall back to common alternatives
        for alt in ("human_symbol", "gene_symbol", "symbol", "gene"):
            if alt in atlas.columns:
                join_key = alt
                break
        else:
            raise ValueError(
                f"Atlas has no join column. Available: {list(atlas.columns[:10])}..."
            )
    print(f"[atlas_writer] Joining on '{join_key}'")

    # Drop existing perturb_* cols if overwrite
    if overwrite:
        existing = [c for c in atlas.columns if c.startswith("perturb_")]
        if existing:
            print(f"[atlas_writer] Overwrite: dropping {len(existing)} existing perturb_* cols")
            atlas = atlas.drop(columns=existing)

    arms = load_per_arm_outputs()

    merged = atlas.copy()
    for arm, df in arms.items():
        # Standardize per-arm CSV's symbol column to match atlas join_key
        if join_key in df.columns:
            df_join = df
        elif "gene" in df.columns:
            df_join = df.rename(columns={"gene": join_key})
        elif "human_symbol" in df.columns:
            df_join = df.rename(columns={"human_symbol": join_key})
        else:
            raise ValueError(f"{arm} output missing gene/human_symbol column. Cols: {list(df.columns)}")
        merged = merged.merge(df_join, on=join_key, how="left", suffixes=("", f"__{arm}"))

    # Fill missing S8 cols if no arm wrote them
    for col in S8_COLUMNS:
        if col not in merged.columns:
            merged[col] = pd.NA

    # Composite tier (returns dataframe with 3 added cols: n_applicable, n_passing, tier)
    merged = compute_validation_tier(merged)

    # Audit
    audit = {
        "n_genes_atlas": len(atlas),
        "n_genes_with_any_perturb_col": int(
            (merged[S8_COLUMNS[:-1]].notna().any(axis=1)).sum()
        ),
        "n_strong": int((merged["perturb_validation_tier"] == "strong").sum()),
        "n_validated": int((merged["perturb_validation_tier"] == "validated").sum()),
        "n_partial": int((merged["perturb_validation_tier"] == "partial").sum()),
        "n_no_evidence": int((merged["perturb_validation_tier"] == "—").sum()),
        "join_key": join_key,
    }
    if dry_run:
        print(f"[atlas_writer] DRY RUN audit: {audit}")
        return audit

    backup_path = backup_atlas()
    print(f"[atlas_writer] Backed up atlas to {backup_path}")
    merged.to_csv(ATLAS_PATH, index=False)
    print(f"[atlas_writer] Wrote S8 columns to {ATLAS_PATH}")
    print(f"[atlas_writer] Audit: {audit}")
    return audit
