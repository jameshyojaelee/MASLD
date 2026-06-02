"""Phase 1 backtest harness.

Runs all candidate models zero-shot on validation set, produces recovery
curves, recommends model eject decisions. Called by Phase 1 orchestrator
once envs + checkpoints are ready.

Used by Cross-arm Coordinator. Each arm's Backtester subagent calls this for
their specific arm via the `arm` parameter.
"""
from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

from backtest_framework import (
    load_backtest_set,
    per_source_recovery,
    recovery_curve,
)

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
OUT_DIR = PROJECT_ROOT / "Analysis/Perturbation/results/backtest"
OUT_DIR.mkdir(parents=True, exist_ok=True)


def run_phase1_backtest(
    *,
    arm: str,
    model_rankings: dict[str, list[str]],  # model_name -> ranked gene list
    n_atlas_total: int = 33943,
    gate_fold: float = 2.0,
) -> dict:
    """Run Phase 1 backtest for an arm.

    Returns: dict with per-model decisions (keep/eject), recovery dataframes,
    output paths.
    """
    backtest = load_backtest_set()
    bt_genes = backtest["gene"].tolist()

    # Per-source recovery for each model
    per_src_dfs = []
    curve_dfs = []
    eject_decisions = {}

    for model_name, ranked in model_rankings.items():
        per_src = per_source_recovery(
            model=model_name,
            arm=arm,
            ranked_genes=ranked,
            n_atlas_total=n_atlas_total,
        )
        per_src_dfs.append(per_src)

        curve = recovery_curve(
            model=model_name,
            arm=arm,
            ranked_genes=ranked,
            backtest_genes=bt_genes,
            n_atlas_total=n_atlas_total,
        )
        curve_dfs.append(curve)

        # Eject decision: fail if fold_enrichment < gate_fold at top 10% on
        # the FULL backtest set (not slice).
        top10 = curve[curve["top_frac"] == 0.10]
        fold: float | None = None
        if len(top10) == 0:
            decision = "no_top10_data"
        else:
            fold = float(top10["fold_enrichment"].iloc[0])
            decision = "KEEP" if fold >= gate_fold else "EJECT"
        eject_decisions[model_name] = {
            "decision": decision,
            "fold_enrichment_top10": fold,
        }

    # Combine
    per_src_all = pd.concat(per_src_dfs, ignore_index=True)
    curves_all = pd.concat(curve_dfs, ignore_index=True)

    per_src_out = OUT_DIR / f"phase1_{arm}_per_source_recovery.csv"
    curves_out = OUT_DIR / f"phase1_{arm}_recovery_curves.csv"
    per_src_all.to_csv(per_src_out, index=False)
    curves_all.to_csv(curves_out, index=False)

    # Plot recovery curves
    fig, ax = plt.subplots(1, 1, figsize=(7, 5))
    for model_name, ranked in model_rankings.items():
        sub = curves_all[curves_all["model"] == model_name]
        ax.plot(
            sub["top_frac"] * 100,
            sub["fold_enrichment"],
            marker="o",
            label=model_name,
        )
    ax.axhline(gate_fold, color="red", ls="--", alpha=0.5, label=f"Eject gate ({gate_fold}x)")
    ax.set_xscale("log")
    ax.set_xlabel("Top % of ranked list")
    ax.set_ylabel("Fold-enrichment vs chance")
    ax.set_title(f"Phase 1 backtest recovery — arm {arm}")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig_out = OUT_DIR / f"phase1_{arm}_recovery_curves.pdf"
    fig.savefig(fig_out)
    plt.close(fig)

    # Eject report
    report = {
        "arm": arm,
        "timestamp": datetime.now().isoformat(),
        "n_atlas_total": n_atlas_total,
        "gate_fold": gate_fold,
        "n_models_evaluated": len(model_rankings),
        "decisions": eject_decisions,
        "n_kept": sum(1 for d in eject_decisions.values() if d["decision"] == "KEEP"),
        "n_ejected": sum(1 for d in eject_decisions.values() if d["decision"] == "EJECT"),
        "per_source_csv": str(per_src_out),
        "curves_csv": str(curves_out),
        "curves_pdf": str(fig_out),
    }
    report_out = OUT_DIR / f"phase1_{arm}_decisions.json"
    with open(report_out, "w") as f:
        json.dump(report, f, indent=2)

    return report


def aggregate_phase1_report(arm_reports: Sequence[dict]) -> Path:
    """Cross-arm Coordinator: aggregate per-arm Phase 1 reports into
    single decisions document used by Phase 2 (fine-tuning)."""
    rows = []
    for r in arm_reports:
        for model, d in r["decisions"].items():
            rows.append(
                dict(
                    arm=r["arm"],
                    model=model,
                    decision=d["decision"],
                    fold_enrichment_top10=d.get("fold_enrichment_top10"),
                )
            )
    df = pd.DataFrame(rows)
    out = OUT_DIR / "phase1_backtest_report.csv"
    df.to_csv(out, index=False)
    return out
