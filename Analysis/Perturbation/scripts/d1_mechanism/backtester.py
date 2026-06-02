#!/usr/bin/env python
"""D1 mechanism backtester.

For each model run, ranks the *downstream* genes by predicted magnitude
(aggregated across all target hits) and computes recovery on the canonical
backtest set. A model whose fold_enrichment < 2x at top-10% gets recommended
for ejection from the consensus pool.
"""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

import pandas as pd  # noqa: E402

from arm_driver import (  # noqa: E402
    RESULTS_ROOT,
    load_runs,
    plot_recovery_curves,
    recommend_ejection,
    run_backtest_for_runs,
)

ARM = "D1"
ARM_DIR = RESULTS_ROOT / "d1_mechanism"
OUT_CSV = ARM_DIR / "recovery_curves.csv"
OUT_PER_SOURCE = ARM_DIR / "recovery_per_source.csv"
OUT_PDF = ARM_DIR / "recovery_curves.pdf"
EJECT_TXT = ARM_DIR / "ejected_models.txt"


def extract_ranked_genes(run):
    """For D1: aggregate predictions, rank downstream genes by mean |logFC|."""
    # Coerce dicts back to objects for attribute access
    score: dict[str, float] = defaultdict(float)
    count: dict[str, int] = defaultdict(int)
    for p in run.predictions:
        g = p["downstream_gene"] if isinstance(p, dict) else p.downstream_gene
        lfc = p["logFC_predicted"] if isinstance(p, dict) else p.logFC_predicted
        score[g] += abs(lfc)
        count[g] += 1
    if not score:
        return []
    df = pd.DataFrame(
        {"gene": list(score), "score": [score[g] / max(1, count[g]) for g in score]}
    ).sort_values("score", ascending=False)
    return df["gene"].tolist()


def main() -> None:
    runs = load_runs(ARM_DIR)
    print(f"[D1/backtester] loaded {len(runs)} ModelRunOutputs")
    curves, per_source = run_backtest_for_runs(
        runs=runs, arm=ARM, extract_ranked_genes=extract_ranked_genes
    )
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    curves.to_csv(OUT_CSV, index=False)
    per_source.to_csv(OUT_PER_SOURCE, index=False)
    plot_recovery_curves(curves, OUT_PDF, title="D1 Mechanism — recovery curves")
    losers = recommend_ejection(curves, fold_threshold=2.0, at_frac=0.10)
    EJECT_TXT.write_text("\n".join(losers))
    print(f"[D1/backtester] wrote {OUT_CSV}, {OUT_PDF}; ejection candidates: {losers}")


if __name__ == "__main__":
    main()
