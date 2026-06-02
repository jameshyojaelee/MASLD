#!/usr/bin/env python
"""D4 circuit backtester. Rank sender ligands by response magnitude per run."""
from __future__ import annotations

import sys
from collections import defaultdict
from pathlib import Path

SHARED = Path(__file__).resolve().parents[1] / "shared"
sys.path.insert(0, str(SHARED))

from arm_driver import (  # noqa: E402
    RESULTS_ROOT,
    load_runs,
    plot_recovery_curves,
    recommend_ejection,
    run_backtest_for_runs,
)

ARM_DIR = RESULTS_ROOT / "d4_circuit"
OUT_CSV = ARM_DIR / "recovery_curves.csv"
OUT_PER_SOURCE = ARM_DIR / "recovery_per_source.csv"
OUT_PDF = ARM_DIR / "recovery_curves.pdf"
EJECT_TXT = ARM_DIR / "ejected_models.txt"


def extract_ranked_genes(run):
    score: dict[str, float] = defaultdict(float)
    for p in run.predictions:
        sender = p["sender_gene"] if isinstance(p, dict) else p.sender_gene
        mag = p["response_magnitude"] if isinstance(p, dict) else p.response_magnitude
        score[sender] += abs(mag)
    return [g for g, _ in sorted(score.items(), key=lambda kv: kv[1], reverse=True)]


def main() -> None:
    runs = load_runs(ARM_DIR)
    curves, per_source = run_backtest_for_runs(
        runs=runs, arm="D4", extract_ranked_genes=extract_ranked_genes
    )
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    curves.to_csv(OUT_CSV, index=False)
    per_source.to_csv(OUT_PER_SOURCE, index=False)
    plot_recovery_curves(curves, OUT_PDF, title="D4 Circuit — recovery curves")
    losers = recommend_ejection(curves, fold_threshold=2.0, at_frac=0.10)
    EJECT_TXT.write_text("\n".join(losers))
    print(f"[D4/backtester] losers={losers}")


if __name__ == "__main__":
    main()
