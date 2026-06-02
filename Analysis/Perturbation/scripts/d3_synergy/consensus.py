#!/usr/bin/env python
"""D3 synergy consensus: synergy_magnitude > 2sigma above additive in >=2 models."""
from __future__ import annotations

import sys
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

from arm_driver import RESULTS_ROOT, load_runs  # noqa: E402
from consensus_aggregator import aggregate_synergy_consensus  # noqa: E402
from output_schema import SynergyPrediction  # noqa: E402

ARM_DIR = RESULTS_ROOT / "d3_synergy"
OUT_CSV = ARM_DIR / "consensus_predictions.csv"


def main() -> None:
    runs = load_runs(ARM_DIR)
    print(f"[D3/consensus] loaded {len(runs)} runs")
    for run in runs:
        run.predictions = [
            p if not isinstance(p, dict) else SynergyPrediction(**p)
            for p in run.predictions
        ]
    df = aggregate_synergy_consensus(runs=runs, sigma_threshold=2.0, min_models=2)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    n_pass = int(df["pass_consensus"].sum()) if "pass_consensus" in df else 0
    print(f"[D3/consensus] wrote {len(df)} ({n_pass} pass) -> {OUT_CSV}")


if __name__ == "__main__":
    main()
