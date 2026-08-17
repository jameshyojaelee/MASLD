#!/usr/bin/env python3
"""Inspect a pinned external H5AD without loading its count matrix into memory."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any

import anndata as ad
import numpy as np
from scipy import sparse

from masld_cl.config import write_json_exclusive


def sha256_path(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def matrix_summary(matrix: Any) -> dict[str, Any]:
    return {
        "class": type(matrix).__name__,
        "dtype": str(matrix.dtype),
        "shape": [int(value) for value in matrix.shape],
    }


def sampled_count_check(matrix: Any) -> dict[str, Any]:
    n_rows = int(matrix.shape[0])
    indices = sorted({0, n_rows // 2, n_rows - 1})
    sample = matrix[indices]
    if sparse.issparse(sample):
        values = sample.data
    else:
        values = np.asarray(sample).ravel()
    finite = bool(np.isfinite(values).all())
    nonnegative = bool((values >= 0).all()) if values.size else True
    integer = bool(np.equal(values, np.floor(values)).all()) if values.size else True
    return {
        "row_indices": indices,
        "stored_values": int(values.size),
        "finite": finite,
        "nonnegative": nonnegative,
        "integer": integer,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--expected-bytes", required=True, type=int)
    parser.add_argument("--expected-dataset-version", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    path = Path(args.h5ad)
    if path.stat().st_size != args.expected_bytes:
        raise ValueError(
            f"external H5AD size mismatch: {path.stat().st_size} != {args.expected_bytes}"
        )
    value = ad.read_h5ad(path, backed="r")
    try:
        obs_fields: dict[str, Any] = {}
        priority = {
            "STUDY", "study", "donor", "donor_id", "donor_uuid", "sample_id",
            "library_uuid", "library_alias", "suspension_type", "assay",
            "cell_type", "lineage", "author_cell_type", "disease", "tissue",
        }
        for column in value.obs.columns:
            if column in priority or any(
                token in column.lower()
                for token in ("study", "donor", "sample", "library", "cell_type", "lineage")
            ):
                counts = value.obs[column].astype(str).value_counts(dropna=False)
                obs_fields[column] = {
                    "categories": int(len(counts)),
                    "top": {str(key): int(count) for key, count in counts.head(200).items()},
                }
        layers = {
            name: matrix_summary(value.layers[name]) for name in value.layers.keys()
        }
        raw = None
        if value.raw is not None:
            raw = {
                "matrix": matrix_summary(value.raw.X),
                "var_columns": list(value.raw.var.columns),
                "sampled_count_check": sampled_count_check(value.raw.X),
            }
        result = {
            "schema_version": "masld-external-reference-inspection-v1",
            "dataset_version_id": args.expected_dataset_version,
            "h5ad": str(path.resolve()),
            "h5ad_bytes": path.stat().st_size,
            "h5ad_sha256": sha256_path(path),
            "shape": [int(value.n_obs), int(value.n_vars)],
            "x": matrix_summary(value.X),
            "x_sampled_count_check": sampled_count_check(value.X),
            "layers": layers,
            "raw": raw,
            "obs_columns": list(value.obs.columns),
            "var_columns": list(value.var.columns),
            "selected_obs_fields": obs_fields,
        }
    finally:
        value.file.close()
    write_json_exclusive(args.output, result)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
