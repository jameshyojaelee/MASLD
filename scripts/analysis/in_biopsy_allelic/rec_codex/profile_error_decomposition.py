"""Descriptive split of profile error into regional means and donor variation.

This changes neither model predictions nor the primary absolute-profile loss.
Centering here defines a diagnostic; it is not a transported prediction rule.
"""
from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys

import numpy as np


def profile_error_decomposition(target, prediction):
    target, prediction = np.asarray(target, dtype=np.float64), np.asarray(prediction, dtype=np.float64)
    if target.ndim != 2 or prediction.shape != target.shape or min(target.shape) < 1:
        raise ValueError("Require matching nonempty donor-by-region arrays")
    if not np.isfinite(target).all() or not np.isfinite(prediction).all():
        raise ValueError("Nonfinite target or prediction")
    error = target - prediction
    donor_loss = np.mean(error * error, axis=1)
    region_mean_error = error.mean(axis=0)
    mean_component = float(np.mean(region_mean_error * region_mean_error))
    within_component = float(np.mean((error - region_mean_error) ** 2))
    primary_mse = float(donor_loss.mean())
    difference = primary_mse - mean_component - within_component
    if abs(difference) > 1e-12 + 1e-10 * primary_mse:
        raise ValueError("Exact profile-error decomposition failed")
    return dict(biological_n=target.shape[0], n_regions=target.shape[1],
                donor_mse=donor_loss.tolist(), primary_mse=primary_mse,
                regional_mean_error_mse=mean_component,
                within_region_across_donor_error_mse=within_component,
                decomposition_roundoff=difference,
                variance_divisor="n, preserving the exact descriptive MSE identity",
                receiving_prediction_transform_fitted=False,
                independent_pair_count_asserted=False)


def check_invariants():
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Execute numerical checks on a compute node")
    # Independently calculated: mean errors [2,3], mean squared=13/2;
    # deviations [-2,0,2] at each region give 8/3; total=55/6.
    target = np.array([[0., 1.], [2., 3.], [4., 5.]])
    zero = np.zeros_like(target)
    original = profile_error_decomposition(target, zero)
    assert np.isclose(original["primary_mse"], 55 / 6, atol=1e-12, rtol=0)
    assert np.isclose(original["regional_mean_error_mse"], 13 / 2, atol=1e-12, rtol=0)
    assert np.isclose(original["within_region_across_donor_error_mse"], 8 / 3, atol=1e-12, rtol=0)
    shifted = profile_error_decomposition(target + [7., -4.], zero)
    assert np.isclose(shifted["regional_mean_error_mse"], 41., atol=1e-12, rtol=0)
    assert np.isclose(shifted["within_region_across_donor_error_mse"], 8 / 3, atol=1e-12, rtol=0)
    assert profile_error_decomposition(target, target)["primary_mse"] == 0
    offset_prediction = profile_error_decomposition(target, target + [7., -4.])
    assert offset_prediction["within_region_across_donor_error_mse"] == 0
    assert offset_prediction["primary_mse"] == 65 / 2
    permuted = profile_error_decomposition(target[[2, 0, 1]], zero)
    repeated_regions = profile_error_decomposition(np.tile(target, (1, 3)), np.tile(zero, (1, 3)))
    for result in (permuted, repeated_regions):
        for key in ("primary_mse", "regional_mean_error_mse", "within_region_across_donor_error_mse"):
            assert np.isclose(result[key], original[key], atol=1e-12, rtol=0)
    for bad_target, bad_prediction in ((target, zero[:2]), (target * np.nan, zero)):
        try:
            profile_error_decomposition(bad_target, bad_prediction)
        except ValueError:
            pass
        else:
            raise AssertionError("Invalid input was accepted")
    report = dict(hand_calculated_identity=True, regional_offset_invariance_of_within_term=True,
                  exact_prediction_and_offset_cases=True, donor_permutation_invariance=True,
                  repeated_region_weight_invariance=True, invalid_inputs_refused=True,
                  receiving_data_read=False, biological_example_n=3, example_regions=2,
                  original=original, shifted=shifted,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  python=sys.version, numpy=np.__version__, slurm_job_id=os.environ["SLURM_JOB_ID"])
    output = Path(os.environ["CODEX_REC_OUTPUT"]) / ("profile_error_decomposition_" + os.environ["SLURM_JOB_ID"])
    output.mkdir(exist_ok=False)
    (output / "summary.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({k:v for k,v in report.items() if isinstance(v, bool)}, sort_keys=True))


if __name__ == "__main__":
    check_invariants()
