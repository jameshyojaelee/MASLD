#!/usr/bin/env python
"""D5 mouse backtester.

For D5 the canonical backtest set is the *human* Conserved_Core + COLOC list,
so we map mouse predictions back to human symbols before recovery.
"""
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
from ortholog_mapper import mouse_to_human  # noqa: E402

ARM_DIR = RESULTS_ROOT / "d5_mouse"
OUT_CSV = ARM_DIR / "recovery_curves.csv"
OUT_PER_SOURCE = ARM_DIR / "recovery_per_source.csv"
OUT_PDF = ARM_DIR / "recovery_curves.pdf"
EJECT_TXT = ARM_DIR / "ejected_models.txt"


def extract_ranked_genes(run):
    score: dict[str, float] = defaultdict(float)
    for p in run.predictions:
        g = p["downstream_gene"] if isinstance(p, dict) else p.downstream_gene
        lfc = p["logFC_predicted"] if isinstance(p, dict) else p.logFC_predicted
        # Map mouse downstream symbol back to human
        for h in mouse_to_human(g) or [g.upper()]:
            score[h] += abs(lfc)
    return [g for g, _ in sorted(score.items(), key=lambda kv: kv[1], reverse=True)]


def main() -> None:
    runs = load_runs(ARM_DIR)
    curves, per_source = run_backtest_for_runs(
        runs=runs, arm="D5", extract_ranked_genes=extract_ranked_genes
    )
    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    curves.to_csv(OUT_CSV, index=False)
    per_source.to_csv(OUT_PER_SOURCE, index=False)
    plot_recovery_curves(curves, OUT_PDF, title="D5 Mouse — recovery curves")
    losers = recommend_ejection(curves, fold_threshold=2.0, at_frac=0.10)
    EJECT_TXT.write_text("\n".join(losers))
    print(f"[D5/backtester] losers={losers}")


if __name__ == "__main__":
    main()
