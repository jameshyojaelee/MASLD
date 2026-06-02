#!/usr/bin/env python
"""D3 synergy backtester.

Rank genes by total synergy contribution (summed across all pairs they appear
in) and recover canonical backtest set.
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

ARM_DIR = RESULTS_ROOT / "d3_synergy"
OUT_CSV = ARM_DIR / "recovery_curves.csv"
OUT_PER_SOURCE = ARM_DIR / "recovery_per_source.csv"
OUT_PDF = ARM_DIR / "recovery_curves.pdf"
EJECT_TXT = ARM_DIR / "ejected_models.txt"


def extract_ranked_genes(run):
    score: dict[str, float] = defaultdict(float)
    for p in run.predictions:
        sigma = p["sigma_above_additive"] if isinstance(p, dict) else p.sigma_above_additive
        for f in ("gene1", "gene2", "gene3", "gene4"):
            g = p.get(f) if isinstance(p, dict) else getattr(p, f, None)
            if g:
                score[g] += abs(sigma)
    if not score:
        return []
    return [g for g, _ in sorted(score.items(), key=lambda kv: kv[1], reverse=True)]


def main() -> None:
    runs = load_runs(ARM_DIR)
    curves, per_source = run_backtest_for_runs(
        runs=runs, arm="D3", extract_ranked_genes=extract_ranked_genes
    )
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    curves.to_csv(OUT_CSV, index=False)
    per_source.to_csv(OUT_PER_SOURCE, index=False)
    plot_recovery_curves(curves, OUT_PDF, title="D3 Synergy — recovery curves")
    losers = recommend_ejection(curves, fold_threshold=2.0, at_frac=0.10)
    EJECT_TXT.write_text("\n".join(losers))
    print(f"[D3/backtester] losers={losers}")


if __name__ == "__main__":
    main()
