"""Load + stratify the 657k hepatocyte atlas for perturbation-model fine-tuning.

Provides:
- load_hepatocyte_adata()         : load annotated atlas (with anndata IORegistryError workaround)
- load_with_joined_metadata()     : auto-join hepatocyte_subtype_metadata.csv into obs
- stratify_by_sex_stage_subtype() : iterator over (context, view) pairs
- pseudobulk_by_context()         : pseudobulk per context for bulk-compat models
- donor_holdout_split()           : reproducible 90/10 donor-level hold-out

Used by every arm's Data Engineer subagent.

CRITICAL DEFAULTS:
- The actual h5ad obs lacks `inferred_sex`, `meta_subtype`, `donor` columns
  (shared_utils_review.md P0-H1). Defaults set to the actual on-disk schema:
    subtype_col = "hepatocyte_subtype_label"  (e.g., "Hep_0_Healthy")
    sex_col     = None                         (until task #15 backfill lands)
    donor_col   = "sample"                     (cell-ID; donor inferred via _parse_donor)
- ad.read_h5ad(path, backed='r') raises IORegistryError on /uns/log1p/base
  null encoding (shared_utils_review.md P0-H2). Workaround in load_hepatocyte_adata.
"""
from __future__ import annotations

import hashlib
import os
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Iterator

import anndata as ad
import numpy as np
import pandas as pd

PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
HEP_ATLAS_PATH = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad"
)
HEP_META_PATH = (
    PROJECT_ROOT
    / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv"
)
HEP_META_WITH_SEX_PATH = HEP_META_PATH.with_name(
    "hepatocyte_subtype_metadata_with_sex.csv"
)

# Actual subtype labels from `hepatocyte_subtype_label` column (e.g., Hep_0_Healthy).
# Meta-subtypes (Progressor/Moderate/Stable/Neutral/Healthy) live in a separate
# mapping at meta_subtype_mapping.csv — join via hepatocyte_subtype_label.
DEFAULT_SUBTYPE_LABELS = None  # None = use all subtypes present
DEFAULT_META_SUBTYPES = ("Progressor", "Moderate", "Stable", "Neutral", "Healthy")
DEFAULT_SEXES = ("F", "M")
COARSE_STAGE = ("Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis")
# Cirrhosis is sparse (127 cells); recommend excluding from fine-tuning per
# substrate_inventory finding #9 (use F0-F3 = Healthy/Steatosis/Steatohepatitis).
FINETUNE_STAGES = ("Healthy", "Steatosis", "Steatohepatitis")


@dataclass
class StratifiedContext:
    sex: str | None  # None when sex-pool (sex backfill not yet available)
    stage: str
    subtype: str | None
    n_cells: int
    cell_index: np.ndarray  # cells in this context


def _parse_donor(sample_id: str) -> str:
    """Best-effort donor extraction from sample-id string.

    Sample IDs from the hep atlas obs `sample` column come from many datasets
    with varying conventions:
      - "GSE244832_S1_AAAA-1"   → donor = "GSE244832_S1"
      - "Liver_Atlas_HD1_cell"  → donor = "Liver_Atlas_HD1"
    Strip the last `_<cellbarcode>` segment.
    """
    if not isinstance(sample_id, str):
        return str(sample_id)
    # Strip trailing -N suffix
    s = re.sub(r"-\d+$", "", sample_id)
    parts = s.split("_")
    # Heuristic: donor is everything up to the last segment if it looks like a cell barcode
    if len(parts) > 1 and re.match(r"^[ACGT]{6,}$", parts[-1]):
        return "_".join(parts[:-1])
    return s


def load_hepatocyte_adata(backed: str | None = None) -> ad.AnnData:
    """Load hepatocyte atlas with anndata IORegistryError workaround.

    backed='r' fails on /uns/log1p/base null encoding in the current anndata
    version. Workaround: read the h5ad file with backed=None (full RAM load)
    and emit a warning if the user requested backed mode. For 657k cells
    this needs ~30-50GB RAM — use bigmem or gpu partition.

    If you genuinely need backed mode, pre-process the h5ad to strip
    /uns/log1p/base first (h5repack / scanpy delete).
    """
    if backed == "r":
        try:
            return ad.read_h5ad(HEP_ATLAS_PATH, backed="r")
        except Exception as e:
            warnings.warn(
                f"backed='r' read failed ({type(e).__name__}: {str(e)[:120]}); "
                f"falling back to full RAM load (~30-50GB needed for 657k cells). "
                f"Use --partition=bigmem or gpu --mem=80G.",
                stacklevel=2,
            )
    return ad.read_h5ad(HEP_ATLAS_PATH)


def load_with_joined_metadata(
    adata: ad.AnnData | None = None,
    *,
    prefer_sex_backfilled: bool = True,
) -> ad.AnnData:
    """Load atlas + auto-join hepatocyte_subtype_metadata CSV into obs.

    Joins on `sample` column. If prefer_sex_backfilled and the
    `_with_sex.csv` exists (task #15 deliverable), uses it; otherwise
    falls back to the plain metadata CSV (no sex column).

    Adds these obs columns when available:
      meta_subtype          (Progressor / Moderate / Stable / Neutral / Healthy)
      inferred_sex_final    (F / M / U)  — only if backfilled metadata exists
      donor_id              (parsed from sample via _parse_donor)
    """
    if adata is None:
        adata = load_hepatocyte_adata()

    # Pick metadata file
    meta_path = HEP_META_PATH
    if prefer_sex_backfilled and HEP_META_WITH_SEX_PATH.exists():
        meta_path = HEP_META_WITH_SEX_PATH

    if not meta_path.exists():
        warnings.warn(f"No metadata CSV at {meta_path}; obs will be unenriched")
        return adata

    meta = pd.read_csv(meta_path)
    join_col = "sample" if "sample" in meta.columns and "sample" in adata.obs.columns else None
    if join_col is None:
        warnings.warn("No `sample` join column found; skipping metadata join")
        return adata

    # Take only new columns from meta (avoid overwriting existing obs)
    new_cols = [c for c in meta.columns if c not in adata.obs.columns and c != join_col]
    if new_cols:
        adata.obs = adata.obs.merge(
            meta[[join_col] + new_cols], on=join_col, how="left"
        ).set_index(adata.obs.index)

    # Add parsed donor_id if not present
    if "donor_id" not in adata.obs.columns:
        adata.obs["donor_id"] = adata.obs["sample"].map(_parse_donor)
    return adata


def stratify_by_sex_stage_subtype(
    adata: ad.AnnData,
    *,
    sexes=DEFAULT_SEXES,
    stages=FINETUNE_STAGES,
    subtypes=DEFAULT_META_SUBTYPES,
    sex_col: str | None = "inferred_sex_final",
    stage_col: str = "disease_stage_coarse",
    subtype_col: str = "meta_subtype",
    min_cells: int = 100,
) -> Iterator[StratifiedContext]:
    """Iterator over stratified contexts with ≥ min_cells.

    DEFAULTS use the post-join schema from load_with_joined_metadata():
      sex_col = "inferred_sex_final" (from task #15 backfill; falls back to None pool)
      stage_col = "disease_stage_coarse" (Healthy/Steatosis/Steatohepatitis/Cirrhosis)
      subtype_col = "meta_subtype" (Progressor/Moderate/Stable/Neutral/Healthy)

    Cirrhosis stage EXCLUDED by default (only 127 cells across the whole atlas).
    Pass stages=COARSE_STAGE explicitly to include it.

    If sex_col is None or absent in obs, yields sex='pool' contexts (no sex
    stratification). This is the safe path until task #15 sex backfill lands.

    If subtype_col is absent, falls back to "hepatocyte_subtype_label" if present.
    """
    obs = adata.obs

    # Subtype column fallback
    if subtype_col not in obs.columns:
        if "hepatocyte_subtype_label" in obs.columns:
            subtype_col = "hepatocyte_subtype_label"
            warnings.warn(
                f"meta_subtype not in obs; using hepatocyte_subtype_label "
                f"(load via load_with_joined_metadata() for proper meta_subtype)"
            )
        else:
            raise KeyError(
                f"Neither meta_subtype nor hepatocyte_subtype_label in obs. "
                f"Available: {list(obs.columns)}"
            )

    # Sex column fallback: pool across sex if not available
    sex_iter = sexes if (sex_col is not None and sex_col in obs.columns) else (None,)
    if sex_iter == (None,):
        warnings.warn(
            f"sex_col='{sex_col}' not in obs (task #15 backfill pending); "
            f"yielding 'pool' contexts (no sex stratification)"
        )

    for sex in sex_iter:
        for stage in stages:
            for subtype in subtypes:
                mask = (
                    (obs[stage_col] == stage) & (obs[subtype_col] == subtype)
                )
                if sex is not None:
                    mask = mask & (obs[sex_col] == sex)
                idx = np.where(mask.values)[0]
                if len(idx) < min_cells:
                    continue
                yield StratifiedContext(
                    sex=sex if sex is not None else "pool",
                    stage=stage,
                    subtype=subtype,
                    n_cells=len(idx),
                    cell_index=idx,
                )


def pseudobulk_by_context(
    adata: ad.AnnData,
    contexts: list[StratifiedContext],
    *,
    layer: str | None = None,
    method: str = "sum",
) -> pd.DataFrame:
    """Pseudobulk per context. Returns (gene × context) matrix.

    For bulk-compat models (scGen sometimes, CPA, Tahoe drug arm).
    Method = 'sum' for raw counts, 'mean' for log-normalised.
    """
    cols = []
    for ctx in contexts:
        X = adata.X if layer is None else adata.layers[layer]
        sub = X[ctx.cell_index, :]
        if method == "sum":
            v = np.asarray(sub.sum(axis=0)).ravel()
        elif method == "mean":
            v = np.asarray(sub.mean(axis=0)).ravel()
        else:
            raise ValueError(f"Unknown method {method}")
        name = f"{ctx.sex}__{ctx.stage}__{ctx.subtype or 'any'}"
        cols.append((name, v))
    df = pd.DataFrame({k: v for k, v in cols}, index=adata.var_names)
    return df


def donor_holdout_split(
    adata: ad.AnnData,
    *,
    donor_col: str = "donor_id",
    holdout_frac: float = 0.10,
    seed: int = 42,
) -> tuple[np.ndarray, np.ndarray]:
    """Reproducible donor-level 90/10 split.

    Returns (train_mask, holdout_mask) over cells.
    Donor assignment is deterministic via sha256(donor || seed).

    Default donor_col is "donor_id" (added by load_with_joined_metadata via
    _parse_donor on the `sample` column). If donor_col is absent, falls back
    to sample-level split (a per-cell split, NOT per-donor — less ideal).
    """
    if donor_col not in adata.obs.columns:
        warnings.warn(
            f"donor_col='{donor_col}' not in obs; falling back to per-sample split. "
            f"Use load_with_joined_metadata() to get donor_id."
        )
        donor_col = "sample" if "sample" in adata.obs.columns else adata.obs.columns[0]

    donors = adata.obs[donor_col].unique()
    keyed = sorted(
        donors,
        key=lambda d: hashlib.sha256(f"{d}_{seed}".encode()).hexdigest(),
    )
    n_holdout = max(1, int(round(holdout_frac * len(keyed))))
    holdout_donors = set(keyed[:n_holdout])
    train_donors = set(keyed[n_holdout:])
    is_holdout = adata.obs[donor_col].isin(holdout_donors).values
    is_train = adata.obs[donor_col].isin(train_donors).values
    return is_train, is_holdout
