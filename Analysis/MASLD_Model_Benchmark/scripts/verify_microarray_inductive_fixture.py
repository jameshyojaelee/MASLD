#!/usr/bin/env python3
"""Verify that frozen microarray transforms never learn from held arrays."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path

import numpy as np


class InductiveFixtureError(RuntimeError):
    """Raised when a held-array isolation invariant fails."""


def fit_quantile_target(training: np.ndarray) -> np.ndarray:
    """Fit one quantile target from feature-by-training-array values only."""
    values = np.asarray(training, dtype=np.float64)
    if values.ndim != 2 or values.shape[0] < 2 or values.shape[1] < 2:
        raise InductiveFixtureError("training matrix must contain features and arrays")
    if not np.isfinite(values).all():
        raise InductiveFixtureError("training matrix contains nonfinite values")
    return np.sort(values, axis=0).mean(axis=1)


def apply_quantile_target(array: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Apply a frozen target to one array without consulting another array."""
    values = np.asarray(array, dtype=np.float64)
    reference = np.asarray(target, dtype=np.float64)
    if values.ndim != 1 or reference.ndim != 1 or values.shape != reference.shape:
        raise InductiveFixtureError("array and frozen target axes differ")
    if not np.isfinite(values).all() or not np.isfinite(reference).all():
        raise InductiveFixtureError("array or frozen target contains nonfinite values")
    order = np.argsort(values, kind="mergesort")
    transformed = np.empty_like(values)
    sorted_values = values[order]
    start = 0
    while start < values.size:
        stop = start + 1
        while stop < values.size and sorted_values[stop] == sorted_values[start]:
            stop += 1
        transformed[order[start:stop]] = reference[start:stop].mean()
        start = stop
    return transformed


def per_array_fractional_rank(array: np.ndarray) -> np.ndarray:
    """Rank one already summarized expression vector independently."""
    values = np.asarray(array, dtype=np.float64)
    if values.ndim != 1 or values.size < 2 or not np.isfinite(values).all():
        raise InductiveFixtureError("rank baseline requires one finite expression vector")
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=np.float64)
    sorted_values = values[order]
    start = 0
    while start < values.size:
        stop = start + 1
        while stop < values.size and sorted_values[stop] == sorted_values[start]:
            stop += 1
        ranks[order[start:stop]] = np.arange(start, stop, dtype=np.float64).mean()
        start = stop
    return ranks / float(values.size - 1)


def sha256_array(value: np.ndarray) -> str:
    array = np.asarray(value, dtype="<f8", order="C")
    payload = (
        b"masld-bench-frozen-array-v1\0"
        + json.dumps(array.shape).encode("ascii")
        + b"\0"
        + array.tobytes(order="C")
    )
    return sha256(payload).hexdigest()


def run_fixture(output: Path) -> dict[str, object]:
    training = np.asarray(
        [
            [10.0, 12.0, 11.0],
            [2.0, 1.0, 3.0],
            [7.0, 8.0, 9.0],
            [4.0, 6.0, 5.0],
            [20.0, 18.0, 19.0],
            [13.0, 14.0, 15.0],
        ],
        dtype=np.float64,
    )
    held_a = np.asarray([900.0, 2.0, 50.0, 8.0, 12.0, 100.0], dtype=np.float64)
    held_b = np.asarray([-50.0, 1_000.0, 0.0, 90.0, 3.0, 4.0], dtype=np.float64)
    altered_held_b = held_b * 1_000.0 + 77.0
    held_ties = np.asarray([4.0, 1.0, 4.0, 9.0, 1.0, 3.0], dtype=np.float64)

    target = fit_quantile_target(training)
    transformed_a_alone = apply_quantile_target(held_a, target)
    _ = apply_quantile_target(held_b, target)
    transformed_a_with_b = apply_quantile_target(held_a, target)
    _ = apply_quantile_target(altered_held_b, target)
    transformed_a_with_altered_b = apply_quantile_target(held_a, target)
    if not (
        np.array_equal(transformed_a_alone, transformed_a_with_b)
        and np.array_equal(transformed_a_alone, transformed_a_with_altered_b)
    ):
        raise InductiveFixtureError("held-array transform depends on another held array")

    shifted_training = training.copy()
    shifted_training[0, 0] += 100.0
    shifted_target = fit_quantile_target(shifted_training)
    if np.array_equal(target, shifted_target):
        raise InductiveFixtureError("training perturbation did not change fitted target")

    rank_a_alone = per_array_fractional_rank(transformed_a_alone)
    _ = per_array_fractional_rank(apply_quantile_target(held_b, target))
    rank_a_after_b = per_array_fractional_rank(transformed_a_alone)
    if not np.array_equal(rank_a_alone, rank_a_after_b):
        raise InductiveFixtureError("per-array rank depends on another held array")

    permutation = np.asarray([2, 4, 0, 5, 1, 3], dtype=np.int64)
    inverse = np.argsort(permutation)
    ties_transformed = apply_quantile_target(held_ties, target)
    permuted_transformed = apply_quantile_target(held_ties[permutation], target)[inverse]
    ties_rank = per_array_fractional_rank(held_ties)
    permuted_rank = per_array_fractional_rank(held_ties[permutation])[inverse]
    if not np.array_equal(ties_transformed, permuted_transformed):
        raise InductiveFixtureError("tie-aware quantile application is feature-order dependent")
    if not np.array_equal(ties_rank, permuted_rank):
        raise InductiveFixtureError("tie-aware rank is feature-order dependent")

    output.mkdir(parents=True, exist_ok=False)
    summary = {
        "schema_version": "masld-bench-microarray-inductive-fixture-v1",
        "status": "pass_synthetic_held_array_isolation",
        "fixture_role": "algorithmic_firewall_only_not_biological_preprocessing",
        "input_axis": "already_summarized_feature_by_array_expression",
        "quantile_target_fit_scope": "outer_training_arrays_only",
        "held_application_scope": "one_array_at_a_time",
        "held_A_invariant_to_held_B_presence": True,
        "held_A_invariant_to_held_B_distribution": True,
        "per_array_rank_invariant_to_other_held_arrays": True,
        "exact_ties_receive_average_target": True,
        "exact_ties_receive_average_fractional_rank": True,
        "tie_transform_invariant_to_feature_permutation": True,
        "tie_rank_invariant_to_feature_permutation": True,
        "training_perturbation_changes_frozen_object": True,
        "training_matrix_sha256": sha256_array(training),
        "frozen_target_sha256": sha256_array(target),
        "held_A_input_sha256": sha256_array(held_a),
        "held_A_transformed_sha256": sha256_array(transformed_a_alone),
        "held_A_rank_sha256": sha256_array(rank_a_alone),
        "raw_probe_grid_rank_as_gene_expression_allowed": False,
        "biological_CEL_summarization_activated": False,
        "model_training_activated": False,
        "labels_read": False,
    }
    (output / "inductive_fixture_audit.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    run_fixture(arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
