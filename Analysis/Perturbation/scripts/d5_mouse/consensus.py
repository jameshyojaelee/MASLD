#!/usr/bin/env python
"""D5 mouse consensus: same downstream-overlap rule as D1, applied on mouse-side."""
from __future__ import annotations

import sys
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

import pandas as pd  # noqa: E402

from arm_driver import RESULTS_ROOT, load_runs  # noqa: E402
from consensus_aggregator import aggregate_downstream_consensus  # noqa: E402
from output_schema import GenePrediction  # noqa: E402

ARM = "D5"
ARM_DIR = RESULTS_ROOT / "d5_mouse"
OUT_CSV = ARM_DIR / "consensus_predictions.csv"


def main() -> None:
    runs = load_runs(ARM_DIR)
    print(f"[D5/consensus] loaded {len(runs)} runs")
    for run in runs:
        run.predictions = [
            p if not isinstance(p, dict) else GenePrediction(**p)
            for p in run.predictions
        ]
    consensus = aggregate_downstream_consensus(
        arm=ARM, runs=runs, min_models=2, direction_must_agree=True
    )
    rows = []
    for c in consensus:
        target, downstream = c.gene.split("->", 1)
        rows.append(dict(
            mouse_target_gene=target,
            mouse_downstream_gene=downstream,
            n_models_agree=c.n_models_agree,
            consensus_direction=c.consensus_direction,
            consensus_magnitude=c.consensus_magnitude,
            pass_consensus=c.pass_consensus,
            contributing_models=";".join(c.contributing_models),
        ))
    df = pd.DataFrame(rows)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"[D5/consensus] wrote {len(df)} rows -> {OUT_CSV}")


if __name__ == "__main__":
    main()
