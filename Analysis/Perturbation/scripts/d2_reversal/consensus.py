#!/usr/bin/env python
"""D2 reversal consensus: gene must rank top-5% in >= 2 models."""
from __future__ import annotations

import sys
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

from arm_driver import RESULTS_ROOT, load_runs  # noqa: E402
from consensus_aggregator import aggregate_reversal_consensus  # noqa: E402
from output_schema import ReversalPrediction  # noqa: E402

ARM_DIR = RESULTS_ROOT / "d2_reversal"
OUT_CSV = ARM_DIR / "consensus_predictions.csv"


def main() -> None:
    runs = load_runs(ARM_DIR)
    print(f"[D2/consensus] loaded {len(runs)} runs")
    for run in runs:
        run.predictions = [
            p if not isinstance(p, dict) else ReversalPrediction(**p)
            for p in run.predictions
        ]
    df = aggregate_reversal_consensus(runs=runs, top_pct=0.05, min_models=2)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    if df.empty:
        df = df.assign()  # ensure header present
    df.to_csv(OUT_CSV, index=False)
    n_pass = int(df["pass_consensus"].sum()) if "pass_consensus" in df else 0
    print(f"[D2/consensus] wrote {len(df)} rows ({n_pass} pass) -> {OUT_CSV}")


if __name__ == "__main__":
    main()
