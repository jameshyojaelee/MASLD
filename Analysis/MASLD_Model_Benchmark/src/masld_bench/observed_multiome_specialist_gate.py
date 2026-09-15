"""Fail-closed comparison check for GSE296875 observed-multiome specialists.

Native model outputs remain in their native lanes.  A specialist enters this
common comparison only through a prediction-level donor-lineage target-profile
projection frozen before outcome scoring.  The check consumes unit metrics only
after that prediction separation has been satisfied.
"""

from __future__ import annotations

from dataclasses import dataclass
from statistics import fmean
from typing import Any, Mapping, Sequence

import numpy as np

from .evaluators.stats import two_way_donor_block_bootstrap


SEEDS = (20260824, 20260825, 20260826, 20260827, 20260828)
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
GENOMIC_BLOCKS = (0, 1, 2, 3, 4)
BASELINE_MODEL_ID = "observed_atac_glm"
COMPARISON_OUTPUT_FAMILY = "donor_lineage_target_profile"
N_DONORS = 39
N_UNIT_ROWS = N_DONORS * len(LINEAGES) * len(GENOMIC_BLOCKS)
BOOTSTRAP_REPLICATES = 10_000
BOOTSTRAP_SEED = 20260825
MINIMUM_RELATIVE_DEVIANCE_GAIN = 0.05
MINIMUM_IMPROVED_LINEAGES = 4
MAXIMUM_LINEAGE_DEGRADATION = 0.02
MINIMUM_POSITIVE_GAIN_SEEDS = 4

NATIVE_OUTPUT_FAMILIES = {
    "base_resolution_accessibility_track",
    "peak_accessibility_rate",
    COMPARISON_OUTPUT_FAMILY,
}


class SpecialistGateError(ValueError):
    """Raised when a specialist comparison is not eligible."""


@dataclass(frozen=True, slots=True)
class UnitAxis:
    """The frozen donor-lineage-by-genomic-block evaluation axis."""

    donor_hash: tuple[str, ...]
    lineage: tuple[str, ...]
    genomic_block: tuple[int, ...]


def validate_unit_axis(axis: UnitAxis) -> None:
    if not (
        len(axis.donor_hash)
        == len(axis.lineage)
        == len(axis.genomic_block)
        == N_UNIT_ROWS
    ):
        raise SpecialistGateError("the exact 975-row unit axis is required")
    if any(not value for value in axis.donor_hash):
        raise SpecialistGateError("donor hashes must be non-empty")
    keys = tuple(zip(axis.donor_hash, axis.lineage, axis.genomic_block, strict=True))
    if len(set(keys)) != N_UNIT_ROWS:
        raise SpecialistGateError("donor-lineage-block rows must be unique")
    donors = tuple(dict.fromkeys(axis.donor_hash))
    if len(donors) != N_DONORS:
        raise SpecialistGateError("exactly 39 biological donors are required")
    if set(axis.lineage) != set(LINEAGES) or set(axis.genomic_block) != set(GENOMIC_BLOCKS):
        raise SpecialistGateError("lineage or genomic-block roster differs")
    observed = set(keys)
    expected = {
        (donor, lineage, block)
        for donor in donors
        for lineage in LINEAGES
        for block in GENOMIC_BLOCKS
    }
    if observed != expected:
        raise SpecialistGateError("the donor-by-lineage-by-block grid must be complete")


def require_exact_axis(candidate: UnitAxis, baseline: UnitAxis) -> None:
    validate_unit_axis(baseline)
    validate_unit_axis(candidate)
    if candidate != baseline:
        raise SpecialistGateError("candidate unit identifiers or order differ from the baseline")


def validate_candidate_metadata(metadata: Mapping[str, Any]) -> None:
    required = {
        "schema_version",
        "dataset_id",
        "stage",
        "model_id",
        "native_output_family",
        "comparison_output_family",
        "common_profile_projection",
        "fixed_seeds",
        "seed_ensemble",
        "prediction_committed_before_outcome",
        "partial_rectangle_used",
        "supports_external_claim",
        "supports_champion_claim",
        "sealed_outcomes_read",
    }
    if set(metadata) != required:
        raise SpecialistGateError("candidate metadata field roster differs")
    if metadata["schema_version"] != "masld-bench-observed-multiome-specialist-candidate-v1":
        raise SpecialistGateError("candidate schema differs")
    if metadata["dataset_id"] != "gse296875" or metadata["stage"] != "development":
        raise SpecialistGateError("candidate dataset or stage differs")
    model_id = metadata["model_id"]
    if not isinstance(model_id, str) or not model_id or model_id == BASELINE_MODEL_ID:
        raise SpecialistGateError("candidate model identity differs")
    native_family = metadata["native_output_family"]
    if native_family not in NATIVE_OUTPUT_FAMILIES:
        raise SpecialistGateError("native output family is not registered")
    if metadata["comparison_output_family"] != COMPARISON_OUTPUT_FAMILY:
        raise SpecialistGateError("cross-family scores cannot enter the common-profile gate")
    projection = metadata["common_profile_projection"]
    if not isinstance(projection, dict) or set(projection) != {
        "projection_id",
        "source_sha256",
        "fit_scope",
        "frozen_before_scoring",
        "held_outcome_access",
        "sealed_outcome_access",
        "native_lane_preserved",
    }:
        raise SpecialistGateError("common-profile projection record differs")
    if not isinstance(projection["projection_id"], str) or not projection["projection_id"]:
        raise SpecialistGateError("a named common-profile projection is required")
    digest = projection["source_sha256"]
    if not isinstance(digest, str) or len(digest) != 64 or any(value not in "0123456789abcdef" for value in digest):
        raise SpecialistGateError("projection source must have an exact SHA-256")
    if projection["fit_scope"] not in {"parameter_free", "outer_training_only"}:
        raise SpecialistGateError("projection may be parameter-free or fit in outer training only")
    if projection["frozen_before_scoring"] is not True or projection["held_outcome_access"] is not False or projection["sealed_outcome_access"] is not False or projection["native_lane_preserved"] is not True:
        raise SpecialistGateError("projection firewall or native-lane separation differs")
    if tuple(metadata["fixed_seeds"]) != SEEDS:
        raise SpecialistGateError("the exact five fixed seeds are required")
    if metadata["seed_ensemble"] != "arithmetic_mean_prediction_before_scoring":
        raise SpecialistGateError("the prediction-level five-seed ensemble differs")
    if metadata["prediction_committed_before_outcome"] is not True:
        raise SpecialistGateError("predictions must be committed before outcome access")
    for field in (
        "partial_rectangle_used",
        "supports_external_claim",
        "supports_champion_claim",
        "sealed_outcomes_read",
    ):
        if metadata[field] is not False:
            raise SpecialistGateError(f"candidate disposition differs: {field}")


def _skills(values: Sequence[float] | np.ndarray, shape: tuple[int, ...], label: str) -> np.ndarray:
    result = np.asarray(values, dtype=np.float64)
    if result.shape != shape or not np.isfinite(result).all():
        raise SpecialistGateError(f"{label} shape or finiteness differs")
    if np.any(result > 1.0 + 1e-12):
        raise SpecialistGateError(f"{label} exceeds the maximum profile-deviance skill")
    return result


def relative_residual_deviance_gain(candidate_skill: np.ndarray, baseline_skill: np.ndarray) -> np.ndarray:
    candidate = np.asarray(candidate_skill, dtype=np.float64)
    baseline = np.asarray(baseline_skill, dtype=np.float64)
    if candidate.shape != baseline.shape:
        raise SpecialistGateError("candidate and baseline metric shapes differ")
    denominator = 1.0 - baseline
    if np.any(denominator <= 0.0):
        raise SpecialistGateError("baseline residual-deviance denominator is not positive")
    return (candidate - baseline) / denominator


def compare_specialist(
    *,
    metadata: Mapping[str, Any],
    candidate_axis: UnitAxis,
    baseline_axis: UnitAxis,
    candidate_ensemble_skill: Sequence[float] | np.ndarray,
    candidate_seed_skill: Sequence[Sequence[float]] | np.ndarray,
    baseline_ensemble_skill: Sequence[float] | np.ndarray,
    baseline_seed_skill: Sequence[Sequence[float]] | np.ndarray,
    n_resamples: int = BOOTSTRAP_REPLICATES,
    bootstrap_seed: int = BOOTSTRAP_SEED,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    """Compare one frozen common-profile projection with ``observed_atac_glm``."""

    validate_candidate_metadata(metadata)
    require_exact_axis(candidate_axis, baseline_axis)
    candidate_ensemble = _skills(candidate_ensemble_skill, (N_UNIT_ROWS,), "candidate ensemble")
    baseline_ensemble = _skills(baseline_ensemble_skill, (N_UNIT_ROWS,), "baseline ensemble")
    candidate_seeds = _skills(candidate_seed_skill, (len(SEEDS), N_UNIT_ROWS), "candidate seed")
    baseline_seeds = _skills(baseline_seed_skill, (len(SEEDS), N_UNIT_ROWS), "baseline seed")

    ensemble_gain = relative_residual_deviance_gain(candidate_ensemble, baseline_ensemble)
    seed_gain = relative_residual_deviance_gain(candidate_seeds, baseline_seeds)
    zeros = np.zeros(N_UNIT_ROWS, dtype=np.float64)
    interval = two_way_donor_block_bootstrap(
        ensemble_gain,
        zeros,
        baseline_axis.donor_hash,
        baseline_axis.genomic_block,
        n_resamples=n_resamples,
        seed=bootstrap_seed,
        expected_n_donors=N_DONORS,
        expected_n_genomic_blocks=len(GENOMIC_BLOCKS),
    )
    lineage_gain = {
        lineage: fmean(
            float(value)
            for value, observed_lineage in zip(ensemble_gain, baseline_axis.lineage, strict=True)
            if observed_lineage == lineage
        )
        for lineage in LINEAGES
    }
    seed_points = {
        str(seed): fmean(float(value) for value in seed_gain[index])
        for index, seed in enumerate(SEEDS)
    }
    improved_lineages = sum(value > 0.0 for value in lineage_gain.values())
    positive_gain_seeds = sum(value > 0.0 for value in seed_points.values())
    checks = {
        "minimum_relative_deviance_gain": interval.estimate >= MINIMUM_RELATIVE_DEVIANCE_GAIN,
        "paired_bootstrap_lower_bound_above_zero": interval.lower > 0.0,
        "minimum_four_improved_lineages": improved_lineages >= MINIMUM_IMPROVED_LINEAGES,
        "no_lineage_worse_by_more_than_0_02": min(lineage_gain.values()) >= -MAXIMUM_LINEAGE_DEGRADATION,
        "gain_positive_in_at_least_four_seeds": positive_gain_seeds >= MINIMUM_POSITIVE_GAIN_SEEDS,
    }
    passed = all(checks.values())
    receipt = {
        "schema_version": "masld-bench-observed-multiome-specialist-comparison-v1",
        "dataset_id": "gse296875",
        "stage": "development",
        "model_id": metadata["model_id"],
        "native_output_family": metadata["native_output_family"],
        "comparison_output_family": COMPARISON_OUTPUT_FAMILY,
        "baseline_model_id": BASELINE_MODEL_ID,
        "biological_donors": N_DONORS,
        "lineages": list(LINEAGES),
        "genomic_blocks": len(GENOMIC_BLOCKS),
        "unit_rows": N_UNIT_ROWS,
        "five_seed_ensemble_relative_residual_deviance_gain": interval.estimate,
        "paired_two_way_bootstrap": interval.to_dict(),
        "lineage_relative_residual_deviance_gain": lineage_gain,
        "improved_lineages": improved_lineages,
        "seed_relative_residual_deviance_gain": seed_points,
        "positive_gain_seeds": positive_gain_seeds,
        "development_gate_checks": checks,
        "development_gate_passed": passed,
        "partial_rectangle_used": False,
        "seeds_are_biological_replicates": False,
        "native_lane_preserved": True,
        "supports_external_claim": False,
        "supports_champion_claim": False,
        "sealed_outcomes_read": False,
    }
    arrays = {
        "unit_hash": np.asarray(baseline_axis.donor_hash),
        "lineage": np.asarray(baseline_axis.lineage),
        "genomic_fold": np.asarray(baseline_axis.genomic_block, dtype=np.int8),
        "ensemble_relative_residual_deviance_gain": ensemble_gain,
        "seed_relative_residual_deviance_gain": seed_gain,
    }
    return receipt, arrays
