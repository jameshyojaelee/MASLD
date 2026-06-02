#!/usr/bin/env python
"""D1 mechanism consensus: 2-of-4 downstream-gene overlap.

Reads every ModelRunOutput JSON in results/d1_mechanism/ (state, scgpt,
geneformer, tahoe — across modalities + contexts) and produces a
consensus_predictions.csv listing target->downstream pairs that pass.
"""
from __future__ import annotations

import sys
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

import pandas as pd  # noqa: E402

from arm_driver import RESULTS_ROOT, load_runs  # noqa: E402
from consensus_aggregator import aggregate_downstream_consensus  # noqa: E402

ARM = "D1"
ARM_DIR = RESULTS_ROOT / "d1_mechanism"
OUT_CSV = ARM_DIR / "consensus_predictions.csv"


def main() -> None:
    runs = load_runs(ARM_DIR)
    print(f"[D1/consensus] loaded {len(runs)} ModelRunOutputs from {ARM_DIR}")
    if not runs:
        print("[D1/consensus] no runs found; nothing to aggregate")
        OUT_CSV.write_text("target_gene,downstream_gene,n_models_agree,"
                           "consensus_direction,consensus_magnitude,pass_consensus,"
                           "contributing_models\n")
        return

    # Patch predictions back into pydantic objects so aggregator can use attributes.
    from output_schema import GenePrediction
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
        rows.append(
            dict(
                target_gene=target,
                downstream_gene=downstream,
                n_models_agree=c.n_models_agree,
                consensus_direction=c.consensus_direction,
                consensus_magnitude=c.consensus_magnitude,
                pass_consensus=c.pass_consensus,
                contributing_models=";".join(c.contributing_models),
            )
        )
    COLS = [
        "target_gene", "downstream_gene", "n_models_agree",
        "consensus_direction", "consensus_magnitude", "pass_consensus",
        "contributing_models",
    ]
    if rows:
        df = pd.DataFrame(rows).sort_values(
            ["pass_consensus", "n_models_agree", "consensus_magnitude"],
            ascending=[False, False, False],
        )
    else:
        df = pd.DataFrame(columns=COLS)
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    n_pass = int(df["pass_consensus"].sum()) if not df.empty else 0
    print(f"[D1/consensus] wrote {len(df)} rows -> {OUT_CSV} "
          f"({n_pass} pass)")


if __name__ == "__main__":
    main()
