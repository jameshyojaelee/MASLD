"""Per-arm driver helpers: load all ModelRunOutput JSONs, aggregate, plot.

Used by every arm's consensus.py + backtester.py + atlas_integrator.py to keep
arm-specific scripts thin.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Iterable, Sequence

# Make `shared/` importable from arm subdirs.
SHARED_DIR = Path(__file__).resolve().parent
if str(SHARED_DIR) not in sys.path:
    sys.path.insert(0, str(SHARED_DIR))

import pandas as pd  # noqa: E402

import backtest_framework as btf  # noqa: E402
from output_schema import ModelRunOutput  # noqa: E402

PROJECT_ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
PERTURB_ROOT = PROJECT_ROOT / "Analysis/Perturbation"
RESULTS_ROOT = PERTURB_ROOT / "results"


def list_run_jsons(arm_dir: Path, modality: str | None = None) -> list[Path]:
    """Return all ModelRunOutput JSONs for an arm."""
    arm_dir = Path(arm_dir)
    files = sorted(arm_dir.glob("*.json"))
    if modality is not None:
        files = [f for f in files if f"_{modality}_" in f.name]
    return files


def load_run(path: Path) -> ModelRunOutput:
    with open(path) as fh:
        raw = json.load(fh)
    # pydantic v2 will revalidate; predictions are kept as dicts and validated
    # only against ModelRunOutput's `list` typing (no per-item discriminator).
    return ModelRunOutput.model_validate(raw)


def load_runs(arm_dir: Path, modality: str | None = None) -> list[ModelRunOutput]:
    return [load_run(p) for p in list_run_jsons(arm_dir, modality)]


# ---------------------------------------------------------------------------
# Backtest harness
# ---------------------------------------------------------------------------


def run_backtest_for_runs(
    *,
    runs: Sequence[ModelRunOutput],
    arm: str,
    extract_ranked_genes,
    n_atlas_total: int = 33943,
    fractions: tuple[float, ...] = (0.01, 0.05, 0.10, 0.20, 0.50),
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run recovery + per-source recovery for every model in `runs`.

    `extract_ranked_genes(run)` is an arm-specific lambda that returns a
    ranked list of gene symbols from a ModelRunOutput.
    """
    if not runs:
        return pd.DataFrame(), pd.DataFrame()

    bt = btf.load_backtest_set()
    backtest_genes = bt["gene"].tolist()

    curve_frames = []
    per_source_frames = []

    for run in runs:
        ranked = extract_ranked_genes(run)
        if not ranked:
            continue
        curve = btf.recovery_curve(
            model=f"{run.model}_{run.modality}",
            arm=arm,
            ranked_genes=ranked,
            backtest_genes=backtest_genes,
            n_atlas_total=n_atlas_total,
            fractions=fractions,
        )
        curve_frames.append(curve)

        psrc = btf.per_source_recovery(
            model=f"{run.model}_{run.modality}",
            arm=arm,
            ranked_genes=ranked,
            n_atlas_total=n_atlas_total,
        )
        per_source_frames.append(psrc)

    return (
        pd.concat(curve_frames, ignore_index=True) if curve_frames else pd.DataFrame(),
        pd.concat(per_source_frames, ignore_index=True) if per_source_frames else pd.DataFrame(),
    )


def plot_recovery_curves(
    curves: pd.DataFrame, out_pdf: Path, title: str
) -> Path:
    """Matplotlib line plot: top_frac (x) vs fold_enrichment (y) per model."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: E402

    fig, ax = plt.subplots(figsize=(6, 4))
    if curves.empty:
        ax.text(0.5, 0.5, "no runs", ha="center", va="center")
    else:
        for model, sub in curves.groupby("model"):
            sub = sub.sort_values("top_frac")
            ax.plot(sub["top_frac"], sub["fold_enrichment"], marker="o", label=model)
        ax.axhline(2.0, ls="--", color="grey", lw=0.8, label="2x gate")
        ax.set_xscale("log")
        ax.set_xlabel("top fraction of ranked list")
        ax.set_ylabel("fold enrichment over chance")
        ax.set_title(title)
        ax.legend(fontsize=7, loc="best")
        ax.grid(alpha=0.3)
    fig.tight_layout()
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, format="pdf")
    plt.close(fig)
    return out_pdf


def recommend_ejection(
    curves: pd.DataFrame, fold_threshold: float = 2.0, at_frac: float = 0.10
) -> list[str]:
    """Return list of models whose fold_enrichment at `at_frac` < threshold."""
    if curves.empty:
        return []
    sub = curves[curves["top_frac"] == at_frac]
    losers = sub[sub["fold_enrichment"] < fold_threshold]["model"].tolist()
    return losers
