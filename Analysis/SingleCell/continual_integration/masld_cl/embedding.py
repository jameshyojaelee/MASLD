"""Immutable embedding bundles shared by training and benchmark code."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

from .contracts import ContractError, sha256_path


REQUIRED_CELL_COLUMNS = {
    "row_index", "cell_id", "library_id", "donor_id", "dataset",
    "audit_cell_type", "predicted_cell_type", "preparation", "technical_batch",
    "strict_reference", "primary_query", "query_control", "analysis_eligible",
}


def load_embedding(manifest_path: str | Path) -> tuple[dict[str, Any], np.ndarray, pd.DataFrame]:
    manifest_path = Path(manifest_path).resolve()
    with manifest_path.open() as handle:
        manifest = json.load(handle)
    if manifest.get("schema_version") != "masld-cl-embedding-v1":
        raise ContractError("unsupported embedding manifest")
    latent_path = manifest_path.parent / manifest["latent_file"]
    cells_path = manifest_path.parent / manifest["cells_file"]
    for path, key in ((latent_path, "latent_sha256"), (cells_path, "cells_sha256")):
        if not path.is_file() or sha256_path(path) != manifest.get(key):
            raise ContractError(f"embedding artifact changed or is missing: {path.name}")
    latent = np.load(latent_path, mmap_mode="r")
    cells = pd.read_csv(cells_path, sep="\t", compression="gzip", low_memory=False)
    if not REQUIRED_CELL_COLUMNS.issubset(cells.columns):
        missing = sorted(REQUIRED_CELL_COLUMNS - set(cells.columns))
        raise ContractError(f"embedding cell metadata is incomplete: {missing}")
    if latent.ndim != 2 or len(cells) != len(latent) or len(cells) != manifest.get("n_cells"):
        raise ContractError("embedding latent and cell roster have inconsistent dimensions")
    if cells["cell_id"].duplicated().any():
        raise ContractError("embedding contains duplicate cell IDs")
    expected_rows = np.arange(len(cells), dtype=np.int64)
    if not np.array_equal(cells["row_index"].to_numpy(), expected_rows):
        raise ContractError("embedding row indices are not exact and ordered")
    if not np.isfinite(latent).all():
        raise ContractError("embedding contains non-finite values")
    for key in ("strict_reference", "primary_query", "query_control", "analysis_eligible"):
        if cells[key].dtype != bool:
            normalized = cells[key].astype(str).str.lower()
            if not normalized.isin(["true", "false"]).all():
                raise ContractError(f"embedding contains invalid booleans: {key}")
            cells[key] = normalized == "true"
    return manifest, latent, cells


def matched_rows(left: pd.DataFrame, right: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    left_index = {value: index for index, value in enumerate(left["cell_id"].astype(str))}
    right_index = {value: index for index, value in enumerate(right["cell_id"].astype(str))}
    shared = sorted(set(left_index) & set(right_index))
    if not shared:
        raise ContractError("embedding bundles have no shared cells")
    return (
        np.asarray([left_index[x] for x in shared], dtype=np.int64),
        np.asarray([right_index[x] for x in shared], dtype=np.int64),
    )
