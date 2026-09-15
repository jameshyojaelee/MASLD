"""Corrected family-native check for observed-multiome specialists."""

from __future__ import annotations

from statistics import fmean
from typing import Any, Mapping, Sequence

import numpy as np

from .evaluators.stats import two_way_donor_block_bootstrap
from .observed_multiome_specialist_gate import (
    BASELINE_MODEL_ID,
    BOOTSTRAP_REPLICATES,
    BOOTSTRAP_SEED,
    COMPARISON_OUTPUT_FAMILY,
    GENOMIC_BLOCKS,
    LINEAGES,
    MAXIMUM_LINEAGE_DEGRADATION,
    MINIMUM_IMPROVED_LINEAGES,
    MINIMUM_POSITIVE_GAIN_SEEDS,
    MINIMUM_RELATIVE_DEVIANCE_GAIN,
    N_DONORS,
    N_UNIT_ROWS,
    SEEDS,
    SpecialistGateError,
    UnitAxis,
    _skills,
    relative_residual_deviance_gain,
    require_exact_axis,
    validate_unit_axis,
)


MODEL_NATIVE_OUTPUT_FAMILY = {
    "epibert": "functional_track",
    "epcotv2": "functional_track",
    "scooby_epicardioids": "functional_track",
    "scooby_neurips": "functional_track",
    "multivi": "joint_representation",
    "peakvi": "joint_representation",
}
NATIVE_OUTPUT_FAMILIES = {
    "functional_track",
    "joint_representation",
    COMPARISON_OUTPUT_FAMILY,
}


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
    if metadata["schema_version"] != "masld-bench-observed-multiome-specialist-candidate-v2":
        raise SpecialistGateError("candidate schema differs")
    if metadata["dataset_id"] != "gse296875" or metadata["stage"] != "development":
        raise SpecialistGateError("candidate dataset or stage differs")
    model_id = metadata["model_id"]
    if model_id not in MODEL_NATIVE_OUTPUT_FAMILY:
        raise SpecialistGateError("candidate model is absent from the frozen specialist roster")
    native_family = metadata["native_output_family"]
    if native_family not in NATIVE_OUTPUT_FAMILIES or native_family != MODEL_NATIVE_OUTPUT_FAMILY[model_id]:
        raise SpecialistGateError("native output family differs from the frozen census")
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
    for field in ("partial_rectangle_used", "supports_external_claim", "supports_champion_claim", "sealed_outcomes_read"):
        if metadata[field] is not False:
            raise SpecialistGateError(f"candidate disposition differs: {field}")


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
    """Compare a frozen common-profile projection with ``observed_atac_glm``."""

    validate_candidate_metadata(metadata)
    require_exact_axis(candidate_axis, baseline_axis)
    candidate_ensemble = _skills(candidate_ensemble_skill, (N_UNIT_ROWS,), "candidate ensemble")
    baseline_ensemble = _skills(baseline_ensemble_skill, (N_UNIT_ROWS,), "baseline ensemble")
    candidate_seeds = _skills(candidate_seed_skill, (len(SEEDS), N_UNIT_ROWS), "candidate seed")
    baseline_seeds = _skills(baseline_seed_skill, (len(SEEDS), N_UNIT_ROWS), "baseline seed")
    ensemble_gain = relative_residual_deviance_gain(candidate_ensemble, baseline_ensemble)
    seed_gain = relative_residual_deviance_gain(candidate_seeds, baseline_seeds)
    interval = two_way_donor_block_bootstrap(
        ensemble_gain,
        np.zeros(N_UNIT_ROWS, dtype=np.float64),
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
    receipt = {
        "schema_version": "masld-bench-observed-multiome-specialist-comparison-v2",
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
        "development_gate_passed": all(checks.values()),
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
