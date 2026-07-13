#!/usr/bin/env python3
"""Shared IO helpers for the atlas-web generator scripts (data contract §0).

Conventions enforced here (so every generator agrees):
  - JSON compact (separators=(",",":")), NaN/+-Inf -> None recursively,
    R-style "TRUE"/"FALSE" strings -> Python bool.
  - Parquet via pyarrow, index=False, +-Inf -> NaN before write; symbol-keyed
    long parquets sorted by `symbol`.
  - Project root from MASLD_PROJECT_ROOT env.
"""

from __future__ import annotations

import json
import math
import os
from pathlib import Path

import numpy as np
import pandas as pd

DEFAULT_ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"


def project_root() -> Path:
    return Path(os.environ.get("MASLD_PROJECT_ROOT", DEFAULT_ROOT))


def r(v, n: int = 4):
    """Round a scalar to n dp; NaN/Inf/None -> None."""
    if v is None:
        return None
    try:
        v = float(v)
    except (TypeError, ValueError):
        return None
    if math.isnan(v) or math.isinf(v):
        return None
    return round(v, n)


def _to_bool(v):
    if isinstance(v, str):
        u = v.strip().upper()
        if u == "TRUE":
            return True
        if u == "FALSE":
            return False
    return v


def clean(obj):
    """Recursively sanitize for JSON: NaN/+-Inf -> None, 'TRUE'/'FALSE' -> bool,
    numpy scalars -> python scalars."""
    if isinstance(obj, dict):
        return {k: clean(val) for k, val in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v) for v in obj]
    if isinstance(obj, (np.floating,)):
        obj = float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.bool_,)):
        return bool(obj)
    if isinstance(obj, float):
        if math.isnan(obj) or math.isinf(obj):
            return None
        return obj
    if isinstance(obj, str):
        return _to_bool(obj)
    return obj


def dump_json(obj, path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(clean(obj), fh, separators=(",", ":"))
    return path.stat().st_size / 1024.0


def write_parquet(df: pd.DataFrame, path: Path, sort_col: str = "symbol"):
    """Write a long parquet: +-Inf -> NaN, optional sort, index=False, pyarrow."""
    df = df.copy()
    num = df.select_dtypes(include=[np.number]).columns
    if len(num):
        df[num] = df[num].replace([np.inf, -np.inf], np.nan)
    if sort_col and sort_col in df.columns:
        df = df.sort_values(sort_col, kind="stable").reset_index(drop=True)
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(path, engine="pyarrow", index=False)
    return df


def strip_version(s):
    """ENSG...12 -> ENSG..."""
    if s is None:
        return None
    s = str(s)
    return s.split(".")[0] if s.startswith("ENSG") else s


def load_ensembl_symbol_map(root: Path) -> dict:
    """version-stripped ensembl_id -> human_symbol from the multi-evidence atlas."""
    atlas = root / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
    a = pd.read_csv(atlas, usecols=["human_symbol", "ensembl_id"], low_memory=False)
    a = a.dropna(subset=["ensembl_id", "human_symbol"])
    a["ensembl_id"] = a["ensembl_id"].map(strip_version)
    return dict(zip(a["ensembl_id"], a["human_symbol"]))
