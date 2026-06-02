#!/usr/bin/env python
"""D4 circuit consensus: 2-of-3 propagation methods agree on receiver response."""
from __future__ import annotations

import sys
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

from arm_driver import RESULTS_ROOT, load_runs  # noqa: E402
from consensus_aggregator import aggregate_circuit_consensus  # noqa: E402
from output_schema import CircuitPrediction  # noqa: E402

ARM_DIR = RESULTS_ROOT / "d4_circuit"
OUT_CSV = ARM_DIR / "consensus_predictions.csv"


def main() -> None:
    runs = load_runs(ARM_DIR)
    print(f"[D4/consensus] loaded {len(runs)} runs")
    for run in runs:
        run.predictions = [
            p if not isinstance(p, dict) else CircuitPrediction(**p)
            for p in run.predictions
        ]
    df = aggregate_circuit_consensus(runs=runs, min_methods=2)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    n_pass = int(df["pass_consensus"].sum()) if "pass_consensus" in df else 0
    print(f"[D4/consensus] wrote {len(df)} ({n_pass} pass) -> {OUT_CSV}")


if __name__ == "__main__":
    main()
