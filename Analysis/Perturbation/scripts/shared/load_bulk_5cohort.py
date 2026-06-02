"""Load 5-cohort canonical bulk pseudobulk for perturbation models that consume bulk.

Returns a (gene × sample) matrix + sample metadata (cohort, sex, stage_coarse).
Used by Data Engineer subagents in D1, D2, D5 where models accept bulk.

5-cohort canonical (2026-05-01): Suppli, Hoang, Govaere, Bril, Chen — n=847 samples,
control-bearing. Gerhard/PRJNA512027 dropped from paper presentation 2026-05-15 (L0/S0
batch confound); yaml retained for Script 05d provenance.

CSV SUBSTRATES (P0 fix 2026-05-21): Canonical .rds files at integration/results/
must be pre-converted to CSV via `convert_bulk_rds_to_csv.R` (run once via SLURM).
Outputs land at Analysis/Perturbation/data/pseudobulks/:
  bulk_counts_raw.csv         (raw counts, gene × sample)
  bulk_logcpm.csv             (batch-corrected log-CPM, gene × sample)
  bulk_dge_metadata.csv       (sample-level DGEList metadata)
"""
from __future__ import annotations

import os
from pathlib import Path

import pandas as pd

PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
INTEG_ROOT = (
    PROJECT_ROOT
    / "RNA-seq/Human/Patient_Cohorts/analysis/integration"
)
PSEUDOBULK_DIR = PROJECT_ROOT / "Analysis/Perturbation/data/pseudobulks"

# Canonical 5-cohort substrates. Prefer pre-converted CSVs at PSEUDOBULK_DIR.
# Fallback to direct integration/ paths (in case they're added later).
COUNTS_CANDIDATES = [
    PSEUDOBULK_DIR / "bulk_counts_raw.csv",
    PSEUDOBULK_DIR / "bulk_logcpm.csv",
    INTEG_ROOT / "results/integration/dream_input_counts.csv",
    INTEG_ROOT / "results/integration/counts_canonical_5cohort.csv",
]
LOGCPM_CANDIDATES = [
    PSEUDOBULK_DIR / "bulk_logcpm.csv",
]
META_CANDIDATES = [
    PSEUDOBULK_DIR / "bulk_dge_metadata.csv",
    INTEG_ROOT / "metadata/unified_metadata.csv",
]

CANONICAL_COHORTS = ("Suppli", "Hoang", "Govaere", "Bril", "Chen")
CONVERSION_SCRIPT = (
    PROJECT_ROOT / "Analysis/Perturbation/scripts/shared/convert_bulk_rds_to_csv.R"
)


def _first_existing(paths: list[Path]) -> Path:
    for p in paths:
        if p.exists():
            return p
    raise FileNotFoundError(
        f"None of these substrates exist:\n"
        + "\n".join(f"  {p}" for p in paths)
        + f"\n\nThe upstream .rds files live at:\n  {INTEG_ROOT / 'results/integration/'}\n"
        + f"Convert them to CSV via SLURM:\n"
        + f"  sbatch --partition=gpu --gres=gpu:0 --mem=40G --time=1:00:00 "
        + f"--wrap='module load R && Rscript {CONVERSION_SCRIPT}'\n"
    )


def load_bulk_counts(prefer: str = "raw") -> pd.DataFrame:
    """Return gene × sample counts matrix.

    prefer ∈ {"raw", "logcpm"}: which substrate to load by default.
    Returns DataFrame indexed by gene.
    """
    cands = LOGCPM_CANDIDATES + COUNTS_CANDIDATES if prefer == "logcpm" else COUNTS_CANDIDATES
    path = _first_existing(cands)
    df = pd.read_csv(path)
    # Standardize: first col is gene symbol
    if "gene" in df.columns:
        df = df.set_index("gene")
    else:
        df = df.set_index(df.columns[0])
    return df


def load_bulk_metadata() -> pd.DataFrame:
    path = _first_existing(META_CANDIDATES)
    return pd.read_csv(path)


def load_5cohort_bulk(
    filter_to_canonical: bool = True,
    prefer: str = "logcpm",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Return (counts_df gene × sample, meta_df sample-level).

    prefer="logcpm" (default for perturbation models that consume normalized expression)
    or "raw" (for models that want raw counts).

    If filter_to_canonical, restrict to CANONICAL_COHORTS (matches mega-analysis 847).
    """
    counts = load_bulk_counts(prefer=prefer)
    meta = load_bulk_metadata()
    if filter_to_canonical and "cohort" in meta.columns:
        meta = meta[meta["cohort"].isin(CANONICAL_COHORTS)].copy()
        sample_col = (
            "sample_id" if "sample_id" in meta.columns
            else ("sample" if "sample" in meta.columns else meta.columns[0])
        )
        keep = [c for c in counts.columns if c in set(meta[sample_col])]
        counts = counts[keep]
        meta = meta.set_index(sample_col).loc[keep].reset_index()
    return counts, meta
