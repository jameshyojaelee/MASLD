"""Backtest framework: positive-control recovery curves for any perturbation model.

Used by every arm's Backtester subagent. Backtester is BLINDED to consensus output
before computing recovery — this preserves consensus statistical validity.

A model whose recovery on the canonical backtest set is at chance gets EJECTED from
the consensus pool for that arm.

Recovery metric: % of validation genes appearing in top-N of model's ranking.
Chance level: baseline_frac = (N_in_validation / N_total_in_atlas).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
BACKTEST_SET = (
    PROJECT_ROOT
    / "Analysis/Perturbation/results/validation_sets/backtest_set.csv"
)


@dataclass
class RecoveryResult:
    model: str
    arm: str
    source: str  # validation source slice (e.g., "conserved_core")
    n_validation: int
    n_atlas_total: int
    top_n: int
    n_recovered: int
    recovery_pct: float
    chance_pct: float
    fold_enrichment: float
    pass_gate: bool  # >= 2× chance counts as "above chance"


def load_backtest_set(sources: list[str] | None = None) -> pd.DataFrame:
    """Load canonical backtest gene set. Optional source filter."""
    df = pd.read_csv(BACKTEST_SET)
    if sources is not None:
        keep = df["source"].apply(
            lambda s: any(src in str(s).split(";") for src in sources)
        )
        df = df[keep]
    return df


def compute_recovery(
    *,
    model: str,
    arm: str,
    ranked_genes: list[str],
    backtest_genes: list[str],
    n_atlas_total: int,
    top_n: int = None,
    gate_fold: float = 2.0,
    source_label: str = "all",
) -> RecoveryResult:
    """Compute recovery at top-N (default top-10% of ranked list).

    fold_enrichment = recovery / chance ; pass_gate iff >= gate_fold (default 2x).
    """
    if top_n is None:
        top_n = max(1, len(ranked_genes) // 10)
    ranked_top = set(ranked_genes[:top_n])
    bt = set(backtest_genes)
    n_recov = len(ranked_top & bt)
    recovery_pct = n_recov / max(1, len(bt)) * 100.0
    chance_pct = top_n / max(1, n_atlas_total) * 100.0
    fold = recovery_pct / max(0.001, chance_pct)
    return RecoveryResult(
        model=model,
        arm=arm,
        source=source_label,
        n_validation=len(bt),
        n_atlas_total=n_atlas_total,
        top_n=top_n,
        n_recovered=n_recov,
        recovery_pct=recovery_pct,
        chance_pct=chance_pct,
        fold_enrichment=fold,
        pass_gate=fold >= gate_fold,
    )


def recovery_curve(
    *,
    model: str,
    arm: str,
    ranked_genes: list[str],
    backtest_genes: list[str],
    n_atlas_total: int,
    fractions: tuple[float, ...] = (0.01, 0.05, 0.10, 0.20, 0.50),
) -> pd.DataFrame:
    """Recovery curve at multiple top-fraction cutoffs."""
    rows = []
    for f in fractions:
        n = max(1, int(len(ranked_genes) * f))
        r = compute_recovery(
            model=model,
            arm=arm,
            ranked_genes=ranked_genes,
            backtest_genes=backtest_genes,
            n_atlas_total=n_atlas_total,
            top_n=n,
        )
        rows.append(
            dict(
                model=r.model,
                arm=r.arm,
                top_frac=f,
                top_n=r.top_n,
                n_recovered=r.n_recovered,
                recovery_pct=r.recovery_pct,
                chance_pct=r.chance_pct,
                fold_enrichment=r.fold_enrichment,
                pass_gate=r.pass_gate,
            )
        )
    return pd.DataFrame(rows)


def per_source_recovery(
    *,
    model: str,
    arm: str,
    ranked_genes: list[str],
    n_atlas_total: int,
    top_n: int | None = None,
) -> pd.DataFrame:
    """Recovery sliced by validation source (conserved_core, coloc, drug_target, cas13).

    Reports per-slice fold-enrichment so we can see whether a model is biased
    (e.g., recovers Conserved_Core but not COLOC).
    """
    bt = load_backtest_set()
    rows = []
    for source in [
        "conserved_core",
        "coloc_pp4_high",
        "drug_target_approved",
        "drug_target_trial",
        "drug_target_genetic_risk",
        "cas13_pos_control",
    ]:
        bt_src = bt[bt["source"].str.contains(source, na=False)]["gene"].tolist()
        if len(bt_src) < 3:
            continue
        r = compute_recovery(
            model=model,
            arm=arm,
            ranked_genes=ranked_genes,
            backtest_genes=bt_src,
            n_atlas_total=n_atlas_total,
            top_n=top_n,
            source_label=source,
        )
        rows.append(
            dict(
                model=r.model,
                arm=r.arm,
                source=source,
                n_validation=r.n_validation,
                top_n=r.top_n,
                n_recovered=r.n_recovered,
                recovery_pct=r.recovery_pct,
                chance_pct=r.chance_pct,
                fold_enrichment=r.fold_enrichment,
                pass_gate=r.pass_gate,
            )
        )
    return pd.DataFrame(rows)
