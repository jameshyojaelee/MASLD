#!/usr/bin/env python
"""D5 atlas integrator -> 1 column: perturb_crossspecies_concordance (0-1)."""
from __future__ import annotations

import sys
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

import pandas as pd  # noqa: E402

from arm_driver import RESULTS_ROOT  # noqa: E402
import concordance_scorer  # noqa: E402

ARM_DIR = RESULTS_ROOT / "d5_mouse"
CONCORDANCE_CSV = ARM_DIR / "concordance_scores.csv"
OUT_CSV = RESULTS_ROOT / "integration" / "d5_atlas_columns.csv"


def main() -> None:
    if not CONCORDANCE_CSV.exists():
        concordance_scorer.score(CONCORDANCE_CSV)
    df = pd.read_csv(CONCORDANCE_CSV)
    if df.empty:
        OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(columns=["gene", "perturb_crossspecies_concordance"]).to_csv(OUT_CSV, index=False)
        print(f"[D5/atlas] empty concordance -> {OUT_CSV}")
        return

    out = (
        df.groupby("human_gene")["jaccard_top_k"].mean().reset_index()
        .rename(columns={"human_gene": "gene", "jaccard_top_k": "perturb_crossspecies_concordance"})
    )
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT_CSV, index=False)
    print(f"[D5/atlas] wrote {len(out)} rows -> {OUT_CSV}")


if __name__ == "__main__":
    main()
