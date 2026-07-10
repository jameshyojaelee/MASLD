#!/usr/bin/env python3
"""generate_drug_targets_data.py
================================
Emit `drug_targets.parquet` for the /drugs page from the CANONICAL, full
drug-target-universe classification (`data/external/drug_targets/
drug_target_classification.tsv`, 27,943 genes keyed by symbol).

Why a separate file (not atlas_core): the /drugs pipeline bars are about the
full drug universe, so a live GROUP BY here hits the literal canonical headline
(masld_approved=2 / masld_clinical=208 / masld_preclinical=1504) EXACTLY. The
atlas_core drug_dev_status column is restricted to the 27,187 atlas genes and
therefore lands at the atlas-coverage ceiling (1/186/1345) — correct for the
gene page, but not the pipeline headline.

Env: `spatial` (NumPy2/pandas ABI). Honors MASLD_PROJECT_ROOT and --output-dir.

Run:
  micromamba run -n spatial python scripts/portal/generate_drug_targets_data.py \
    --output-dir masld-atlas-v2/public/data
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
DRUG_CLASS_TSV = (
    PROJECT_ROOT / "data/external/drug_targets/drug_target_classification.tsv"
)

# Lean but useful column set (all cheap, present in the TSV).
KEEP_COLS = [
    "symbol",
    "ensembl_id",
    "in_atlas",
    "drug_dev_status",
    "max_phase_masld",
    "max_phase_any",
    "has_masld_trial",
    "n_masld_trials",
    "masld_preclinical_evidence",
    "pharos_tdl",
    "dgidb_n_drugs",
    "dgidb_n_approved",
]

# Collapse duplicate symbols keeping the most-advanced status.
STATUS_RANK = {
    "masld_approved": 0, "masld_clinical": 1, "masld_preclinical": 2,
    "masld_discontinued": 3, "drugged_other_indication": 4,
    "discovery": 5, "undetermined": 6,
}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--output-dir",
        default=str(PROJECT_ROOT / "masld-atlas-v2/public/data"),
        help="Output directory for drug_targets.parquet",
    )
    args = ap.parse_args()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if not DRUG_CLASS_TSV.exists():
        raise SystemExit(f"ERROR: drug class TSV not found: {DRUG_CLASS_TSV}")

    print(f"Loading {DRUG_CLASS_TSV} ...")
    df = pd.read_csv(DRUG_CLASS_TSV, sep="\t", low_memory=False)
    print(f"  TSV shape: {df.shape}")

    cols = [c for c in KEEP_COLS if c in df.columns]
    missing = [c for c in KEEP_COLS if c not in df.columns]
    if missing:
        print(f"  [warn] TSV missing requested cols (skipped): {missing}")
    out = df[cols].copy()

    # Drop rows with no symbol; collapse dup symbols by most-advanced status.
    out = out.dropna(subset=["symbol"])
    out["_r"] = out["drug_dev_status"].map(STATUS_RANK).fillna(9)
    out = (
        out.sort_values("_r")
        .drop_duplicates("symbol", keep="first")
        .drop(columns="_r")
        .sort_values("symbol")
        .reset_index(drop=True)
    )

    # Clean dtypes: ±Inf -> NaN; in_atlas -> bool.
    num_cols = out.select_dtypes(include=[np.number]).columns
    out[num_cols] = out[num_cols].replace([np.inf, -np.inf], np.nan)
    if "in_atlas" in out.columns:
        out["in_atlas"] = out["in_atlas"].map(
            lambda x: str(x).strip().upper() in ("TRUE", "1", "1.0")
        )

    out_path = out_dir / "drug_targets.parquet"
    out.to_parquet(out_path, index=False, engine="pyarrow")
    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(f"  drug_targets.parquet: {size_mb:.2f} MB, {out.shape[0]} rows x {out.shape[1]} cols")

    vc = out["drug_dev_status"].value_counts()
    print("  drug_dev_status GROUP BY (full universe):")
    for k in ("masld_approved", "masld_clinical", "masld_preclinical",
              "masld_discontinued", "drugged_other_indication",
              "discovery", "undetermined"):
        print(f"    {k}: {int(vc.get(k, 0))}")
    print("  columns:", list(out.columns))


if __name__ == "__main__":
    main()
