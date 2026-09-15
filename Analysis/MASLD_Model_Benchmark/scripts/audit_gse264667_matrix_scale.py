#!/usr/bin/env python3
"""Scan GSE264667 HepG2 matrix scale and row-sum provenance."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import h5py
import numpy as np


class MatrixScaleError(RuntimeError):
    """Raised when the bounded numeric scale audit differs."""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not args.source.is_file() or args.output.exists():
        raise MatrixScaleError("matrix scale request differs")
    block_rows = 256
    total_values = 0
    zero_values = 0
    noninteger_values = 0
    nonfinite_values = 0
    negative_values = 0
    observed_min = math.inf
    observed_max = -math.inf
    maximum_row_sum_difference = 0.0
    minimum_signed_row_sum_difference = math.inf
    maximum_signed_row_sum_difference = -math.inf
    row_sum_equal_count = 0
    row_sum_less_than_or_equal_count = 0
    with h5py.File(args.source, mode="r", locking=True) as handle:
        matrix = handle["X"]
        umi_counts = handle["obs/UMI_count"]
        if matrix.shape != (145473, 9624) or umi_counts.shape != (145473,):
            raise MatrixScaleError("matrix or UMI-count axes differ")
        for start in range(0, matrix.shape[0], block_rows):
            stop = min(start + block_rows, matrix.shape[0])
            block = np.asarray(matrix[start:stop, :], dtype=np.float32)
            declared = np.asarray(umi_counts[start:stop], dtype=np.float64)
            total_values += int(block.size)
            zero_values += int(np.count_nonzero(block == 0))
            noninteger_values += int(np.count_nonzero(block != np.floor(block)))
            nonfinite_values += int(np.count_nonzero(~np.isfinite(block)))
            negative_values += int(np.count_nonzero(block < 0))
            observed_min = min(observed_min, float(np.min(block)))
            observed_max = max(observed_max, float(np.max(block)))
            row_sums = np.sum(block, axis=1, dtype=np.float64)
            signed_differences = row_sums - declared
            maximum_row_sum_difference = max(
                maximum_row_sum_difference,
                float(np.max(np.abs(signed_differences))),
            )
            minimum_signed_row_sum_difference = min(
                minimum_signed_row_sum_difference,
                float(np.min(signed_differences)),
            )
            maximum_signed_row_sum_difference = max(
                maximum_signed_row_sum_difference,
                float(np.max(signed_differences)),
            )
            row_sum_equal_count += int(np.count_nonzero(signed_differences == 0))
            row_sum_less_than_or_equal_count += int(np.count_nonzero(signed_differences <= 0))
    result = {
        "schema_version": "masld-bench-gse264667-hepg2-matrix-scale-v2",
        "dataset_id": "gse264667",
        "cell_line": "HepG2",
        "matrix_shape": [145473, 9624],
        "matrix_dtype": "float32",
        "block_rows": block_rows,
        "total_values_scanned": total_values,
        "zero_value_count": zero_values,
        "nonzero_value_count": total_values - zero_values,
        "noninteger_value_count": noninteger_values,
        "nonfinite_value_count": nonfinite_values,
        "negative_value_count": negative_values,
        "minimum": observed_min,
        "maximum": observed_max,
        "maximum_absolute_X_row_sum_minus_obs_UMI_count": maximum_row_sum_difference,
        "minimum_signed_X_row_sum_minus_obs_UMI_count": minimum_signed_row_sum_difference,
        "maximum_signed_X_row_sum_minus_obs_UMI_count": maximum_signed_row_sum_difference,
        "X_row_sum_equals_obs_UMI_count_count": row_sum_equal_count,
        "X_row_sum_less_than_or_equal_obs_UMI_count_count": row_sum_less_than_or_equal_count,
        "integer_nonnegative_finite_matrix_supported": bool(
            noninteger_values == 0 and nonfinite_values == 0 and negative_values == 0
        ),
        "raw_count_scale_supported": bool(
            noninteger_values == 0
            and nonfinite_values == 0
            and negative_values == 0
            and maximum_row_sum_difference == 0
        ),
        "normalization_performed": False,
        "feature_selection_performed": False,
        "model_fit_performed": False,
        "biological_unit_count": 1,
        "gem_groups_as_biological_replicates": False,
        "gse313774_accessed": False,
    }
    args.output.mkdir(parents=True, exist_ok=False)
    (args.output / "matrix_scale.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
