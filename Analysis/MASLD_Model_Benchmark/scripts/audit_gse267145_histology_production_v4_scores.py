#!/usr/bin/env python3
"""Independently audit frozen GSE267145 v4 development score output files."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.contracts import ArtifactRef, MissingState, PredictionBundle
from masld_bench.hashing import canonical_sha256, sha256_file

import scripts.evaluate_gse267145_histology_predictions as reference


TASK_ID = "paired_bulk_histology_state"
DATASET_IDS = ("gse267145_znf469_human_liver",)
SPLIT_ID = "gse267145_participant_outer_v1"
UNIT_NAMESPACE = "gse267145_participant_id"
SCORING_REVISION_ID = "model-scoring-078-production-v4"
PRODUCTION_REVISION_ID = "model-training-069-production-v4"
EXPECTED_PARTICIPANTS = 99
EXPECTED_PREDICTION_BUNDLES = 55
EXPECTED_STAGE5_DISPOSITIONS = 11
EXPECTED_METRICS_PER_BUNDLE = 9
EXPECTED_METRIC_ROWS = 495
EXPECTED_COMPUTED_METRIC_ROWS = 99
EXPECTED_NOT_APPLICABLE_METRIC_ROWS = 396
EXPECTED_SEX_AUDIT_ROWS = 22
BOOTSTRAP_REPLICATES = 10_000
BOOTSTRAP_SEED = 20260824
REFERENCE_EVALUATOR_SHA256 = (
    "1d753b8e23fda8f9582e4d8cc2ce9c90eeae65d370ce342bf9d32ffdd2ba5920"
)
TASK_CONTRACT_SHA256 = (
    "53a85a843e75b1c0fa088af86cb2f2593976d28039807dfa17d1ec571603778c"
)
SURFACE_SHA256 = "23367a28bd16266c79adde854534b124cfad7d2062788625839e78f3267a8c0b"
OUTCOMES_ARTIFACTS_SHA256 = (
    "98a74b97f6a9712a109d2516c3803a4acf15030ed607f861e1ce8c22b4f9753a"
)
FOLDS_ARTIFACTS_SHA256 = (
    "9f5b96e69f217ba643b4b2fd4f3e165049b0cc2bfe9d65c99a7789816b4d7ee3"
)
PRODUCTION_AGGREGATE_ARTIFACTS_SHA256 = (
    "fa6b0ca7512c7ec1bd615d7f981b658d1d0c9bab2053ac6254b55cf8eb0458ac"
)
MULTIPLICITY_FAMILY = "gse267145_histology_development_secondary"
MULTIPLICITY_SCOPE = "within_model_comparison_secondary_endpoint_family"
MULTIPLICITY_METHOD = "bh_q0.05_within_prespecified_secondary_endpoint_family"
NO_CONFIRMATORY_STATE = "not_calculated_development_no_confirmatory_lock"
CLAIM_MODE = "single_cohort_participant_held_development_only"
SOURCE_STAGE5_DISPOSITION = "not_applicable_no_registered_prediction"
SEED_AGGREGATION = "mean_predictions_across_all_five_fixed_model_seeds_before_scoring"
PREDICTED_ENDPOINTS = (
    "stage3",
    "fibrosis_group3",
    "fibrosis_cumulative",
    "fibrosis_regression",
    "nash_crn_component_sum",
)
ENDPOINT_VALUE_FIELDS = {
    "stage3": (
        "probability_NOR",
        "probability_NAFL",
        "probability_NASH",
        "predicted_stage3",
    ),
    "fibrosis_group3": (
        "probability_fibrosis_F0",
        "probability_fibrosis_F1",
        "probability_fibrosis_F2_3",
        "predicted_fibrosis_group3",
    ),
    "fibrosis_cumulative": ("predicted_fibrosis_cumulative_expected",),
    "fibrosis_regression": ("predicted_fibrosis_regression",),
    "nash_crn_component_sum": ("predicted_nash_crn_component_sum",),
}
ENDPOINT_METRICS = {
    "stage3": {"stage3_macro_f1", "stage3_multiclass_brier"},
    "fibrosis_group3": {"fibrosis_group3_macro_f1"},
    "fibrosis_cumulative": {
        "fibrosis_cumulative_ordinal_mae",
        "fibrosis_cumulative_spearman",
    },
    "fibrosis_regression": {
        "fibrosis_regression_mae",
        "fibrosis_regression_spearman",
    },
    "nash_crn_component_sum": {
        "nash_crn_component_sum_spearman",
        "nash_crn_component_sum_mae",
    },
}
NUMERIC_PREDICTION_FIELDS = (
    "probability_NOR",
    "probability_NAFL",
    "probability_NASH",
    "predicted_nash_crn_component_sum",
    "predicted_fibrosis_cumulative_expected",
    "predicted_fibrosis_regression",
    "probability_fibrosis_F0",
    "probability_fibrosis_F1",
    "probability_fibrosis_F2_3",
)
STANDARDIZED_BASE_FIELDS = (
    "prediction_row_id",
    "participant_id",
    "outer_fold",
    "endpoint_id",
    "prediction_state",
)
ROW_ID_FIELDS = ("prediction_row_id", "participant_id")
METRIC_FIELDS = (
    "task_id",
    "model_id",
    "endpoint_id",
    "prediction_aggregation",
    "prediction_set_id",
    "metric_id",
    "endpoint_role",
    "analysis_role",
    "applicability_state",
    "applicability_reason",
    "direction",
    "estimate",
    "ci95_low",
    "ci95_high",
    "valid_bootstrap_replicates",
    "biological_resampling_unit",
    "multiplicity_family",
    "multiplicity_scope",
    "multiplicity_method",
    "p_value",
    "bh_adjusted_q_value",
    "confirmatory_inference_allowed",
    "claim_mode",
)
METRIC_SPECS = (
    ("stage3_macro_f1", "primary", "maximize"),
    ("stage3_multiclass_brier", "calibration", "minimize"),
    ("fibrosis_group3_macro_f1", "secondary", "maximize"),
    ("fibrosis_cumulative_ordinal_mae", "secondary", "minimize"),
    ("fibrosis_cumulative_spearman", "secondary", "maximize"),
    ("fibrosis_regression_mae", "secondary", "minimize"),
    ("fibrosis_regression_spearman", "secondary", "maximize"),
    ("nash_crn_component_sum_spearman", "secondary", "maximize"),
    ("nash_crn_component_sum_mae", "secondary", "minimize"),
)


class HistologyProductionScoreAuditError(RuntimeError):
    """Raised when score output files do not rederive from frozen inputs."""


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise HistologyProductionScoreAuditError(f"invalid JSON: {path}") from error
    if not isinstance(value, dict):
        raise HistologyProductionScoreAuditError(f"JSON is not an object: {path}")
    return value


def _read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise HistologyProductionScoreAuditError(f"TSV has no header: {path}")
            return tuple(reader.fieldnames), [dict(row) for row in reader]
    except OSError as error:
        raise HistologyProductionScoreAuditError(f"cannot read TSV: {path}") from error


def _verify_file(path: Path, expected: str, label: str) -> None:
    if not path.is_file() or path.is_symlink() or sha256_file(path) != expected:
        raise HistologyProductionScoreAuditError(f"{label} SHA-256 differs")


def _prediction_sets() -> list[tuple[str, str]]:
    return [
        (model_id, endpoint_id)
        for model_id in reference.MODEL_IDS
        for endpoint_id in PREDICTED_ENDPOINTS
    ]


def _prediction_row_id(
    production_sha256: str,
    model_id: str,
    endpoint_id: str,
    participant_id: str,
    outer_fold: str,
) -> str:
    return canonical_sha256(
        {
            "production_artifacts_sha256": production_sha256,
            "model_id": model_id,
            "endpoint_id": endpoint_id,
            "participant_id": participant_id,
            "outer_fold": outer_fold,
        }
    )


def _load_outcomes(
    outcomes: Path, folds: Path
) -> tuple[list[str], list[str], dict[str, np.ndarray]]:
    fold_fields, fold_rows = _read_tsv(folds / "participant_outer_folds.tsv")
    endpoint_fields, endpoint_rows = _read_tsv(outcomes / "participant_endpoints.tsv")
    participants = [row["participant_id"] for row in fold_rows]
    outer_folds = [row["outer_fold"] for row in fold_rows]
    if (
        fold_fields != ("participant_id", "outer_fold")
        or len(fold_rows) != EXPECTED_PARTICIPANTS
        or len(set(participants)) != EXPECTED_PARTICIPANTS
        or [row["participant_id"] for row in endpoint_rows] != participants
        or [row["outer_fold"] for row in endpoint_rows] != outer_folds
        or not {
            "stage3",
            "fibrosis",
            "nash_crn_component_sum",
            "recorded_sex",
        }
        <= set(endpoint_fields)
    ):
        raise HistologyProductionScoreAuditError("participant outcome/fold join differs")
    fibrosis = np.asarray([int(row["fibrosis"]) for row in endpoint_rows], dtype=float)
    stage3 = np.asarray([row["stage3"] for row in endpoint_rows], dtype=str)
    sex = np.asarray([row["recorded_sex"] for row in endpoint_rows], dtype=str)
    if (
        Counter(outer_folds) != {"0": 21, "1": 21, "2": 21, "3": 19, "4": 17}
        or Counter(sex) != {"F": 85, "M": 14}
        or set(stage3[sex == "M"]) != {"NASH"}
        or Counter(fibrosis)[3.0] != 4
    ):
        raise HistologyProductionScoreAuditError("participant census differs")
    return participants, outer_folds, {
        "stage3": stage3,
        "fibrosis": fibrosis,
        "fibrosis_group3": np.asarray(
            ["F0" if value == 0 else "F1" if value == 1 else "F2_3" for value in fibrosis],
            dtype=str,
        ),
        "nas": np.asarray(
            [int(row["nash_crn_component_sum"]) for row in endpoint_rows],
            dtype=float,
        ),
        "sex": sex,
    }


def _prediction_arrays(rows: Sequence[Mapping[str, str]]) -> dict[str, np.ndarray]:
    stage_probabilities = np.asarray(
        [
            [float(row[f"probability_{label}"]) for label in reference.STAGE3]
            for row in rows
        ],
        dtype=float,
    )
    fibrosis_probabilities = np.asarray(
        [
            [
                float(row[f"probability_fibrosis_{label}"])
                for label in reference.FIBROSIS_GROUP3
            ]
            for row in rows
        ],
        dtype=float,
    )
    result = {
        "stage_probabilities": stage_probabilities,
        "stage3": np.asarray([row["predicted_stage3"] for row in rows], dtype=str),
        "fibrosis_group3": np.asarray(
            [row["predicted_fibrosis_group3"] for row in rows], dtype=str
        ),
        "fibrosis_cumulative": np.asarray(
            [float(row["predicted_fibrosis_cumulative_expected"]) for row in rows]
        ),
        "fibrosis_regression": np.asarray(
            [float(row["predicted_fibrosis_regression"]) for row in rows]
        ),
        "nas": np.asarray(
            [float(row["predicted_nash_crn_component_sum"]) for row in rows]
        ),
    }
    if (
        not np.all(np.isfinite(stage_probabilities))
        or not np.all(np.isfinite(fibrosis_probabilities))
        or np.any(stage_probabilities < 0.0)
        or np.any(fibrosis_probabilities < 0.0)
        or not np.allclose(stage_probabilities.sum(axis=1), 1.0, atol=1e-6, rtol=0.0)
        or not np.allclose(
            fibrosis_probabilities.sum(axis=1), 1.0, atol=1e-6, rtol=0.0
        )
        or not set(result["stage3"]) <= set(reference.STAGE3)
        or not set(result["fibrosis_group3"]) <= set(reference.FIBROSIS_GROUP3)
        or any(
            not np.all(np.isfinite(result[key]))
            for key in ("fibrosis_cumulative", "fibrosis_regression", "nas")
        )
    ):
        raise HistologyProductionScoreAuditError("prediction values differ")
    return result


def _frozen_bundle_prediction_arrays(
    *,
    prediction_collection: Path,
    model_id: str,
    participants: Sequence[str],
    outer_folds: Sequence[str],
) -> dict[str, np.ndarray]:
    """Assemble one participant-once prediction axis from five endpoint bundles."""
    combined = [
        {"participant_id": participant, "outer_fold": outer_fold}
        for participant, outer_fold in zip(participants, outer_folds, strict=True)
    ]
    for endpoint_id in PREDICTED_ENDPOINTS:
        fields, rows = _read_tsv(
            prediction_collection
            / f"{model_id}--endpoint-{endpoint_id}/predictions.tsv"
        )
        expected_fields = STANDARDIZED_BASE_FIELDS + ENDPOINT_VALUE_FIELDS[endpoint_id]
        if (
            fields != expected_fields
            or len(rows) != EXPECTED_PARTICIPANTS
            or [row["participant_id"] for row in rows] != list(participants)
            or [row["outer_fold"] for row in rows] != list(outer_folds)
            or any(row["endpoint_id"] != endpoint_id for row in rows)
            or any(row["prediction_state"] != "observed" for row in rows)
        ):
            raise HistologyProductionScoreAuditError(
                "frozen endpoint PredictionBundle scoring axis differs"
            )
        for target, row in zip(combined, rows, strict=True):
            target.update(
                {field: row[field] for field in ENDPOINT_VALUE_FIELDS[endpoint_id]}
            )
    return _prediction_arrays(combined)


def _verify_five_seed_ensembles(
    *,
    prediction_root: Path,
    participants: Sequence[str],
    outer_folds: Sequence[str],
) -> tuple[dict[str, list[dict[str, str]]], dict[str, dict[str, Any]]]:
    expected_files = {
        f"{model_id}.tsv" for model_id in reference.MODEL_IDS
    } | {
        f"{model_id}--seed-{seed}.tsv"
        for model_id in reference.MODEL_IDS
        for seed in reference.SEEDS
    }
    if {path.name for path in prediction_root.glob("*.tsv")} != expected_files:
        raise HistologyProductionScoreAuditError("production prediction roster differs")
    ensembles: dict[str, list[dict[str, str]]] = {}
    lineage: dict[str, dict[str, Any]] = {}
    for model_id in reference.MODEL_IDS:
        seed_rows: list[list[dict[str, str]]] = []
        seed_hashes: dict[str, str] = {}
        for seed in reference.SEEDS:
            source = prediction_root / f"{model_id}--seed-{seed}.tsv"
            fields, rows = _read_tsv(source)
            if (
                fields != reference.PREDICTION_FIELDS
                or len(rows) != EXPECTED_PARTICIPANTS
                or [row["participant_id"] for row in rows] != list(participants)
                or [row["outer_fold"] for row in rows] != list(outer_folds)
            ):
                raise HistologyProductionScoreAuditError("seed prediction axis differs")
            _prediction_arrays(rows)
            seed_rows.append(rows)
            seed_hashes[str(seed)] = sha256_file(source)
        ensemble_path = prediction_root / f"{model_id}.tsv"
        fields, ensemble_rows = _read_tsv(ensemble_path)
        if fields != reference.PREDICTION_FIELDS or len(ensemble_rows) != EXPECTED_PARTICIPANTS:
            raise HistologyProductionScoreAuditError("ensemble schema differs")
        for row_index, observed in enumerate(ensemble_rows):
            members = [rows[row_index] for rows in seed_rows]
            averaged = {
                field: math.fsum(float(row[field]) for row in members)
                / len(reference.SEEDS)
                for field in NUMERIC_PREDICTION_FIELDS
            }
            stage = [averaged[f"probability_{label}"] for label in reference.STAGE3]
            fibrosis = [
                averaged[f"probability_fibrosis_{label}"]
                for label in reference.FIBROSIS_GROUP3
            ]
            expected = {
                "participant_id": participants[row_index],
                "outer_fold": outer_folds[row_index],
                **{
                    field: format(averaged[field], ".17g")
                    for field in NUMERIC_PREDICTION_FIELDS[:3]
                },
                "predicted_stage3": reference.STAGE3[
                    max(range(3), key=stage.__getitem__)
                ],
                **{
                    field: format(averaged[field], ".17g")
                    for field in NUMERIC_PREDICTION_FIELDS[3:6]
                },
                **{
                    field: format(averaged[field], ".17g")
                    for field in NUMERIC_PREDICTION_FIELDS[6:]
                },
                "predicted_fibrosis_group3": reference.FIBROSIS_GROUP3[
                    max(range(3), key=fibrosis.__getitem__)
                ],
            }
            if observed != expected:
                raise HistologyProductionScoreAuditError(
                    "ensemble is not the exact five-seed mean"
                )
        ensembles[model_id] = ensemble_rows
        lineage[model_id] = {
            "ensemble_prediction_sha256": sha256_file(ensemble_path),
            "seed_prediction_sha256": seed_hashes,
            "five_seed_membership_sha256": canonical_sha256(seed_hashes),
        }
    return ensembles, lineage


def _reference_points(
    outcomes: Mapping[str, np.ndarray], predictions: Mapping[str, np.ndarray]
) -> dict[str, str]:
    functions = reference._metric_functions(outcomes, predictions)
    axis = np.arange(EXPECTED_PARTICIPANTS)
    result: dict[str, str] = {}
    for metric_id, _, _ in METRIC_SPECS:
        value = float(functions[metric_id](axis))
        result[metric_id] = (
            "not_estimable" if not math.isfinite(value) else format(value, ".17g")
        )
    return result


def _metric_policy(role: str) -> tuple[str, str, str, str, str]:
    if role == "primary":
        return (
            "not_applicable_primary_endpoint",
            "not_applicable_primary_endpoint",
            "primary_effect_and_interval_without_confirmatory_p_value",
            "not_calculated_primary_development_only",
            "not_applicable_primary_development_only",
        )
    if role == "calibration":
        return (
            "not_applicable_calibration_endpoint",
            "not_applicable_calibration_endpoint",
            "not_applicable_calibration_endpoint",
            "not_calculated_calibration_development_only",
            "not_applicable_calibration_development_only",
        )
    return (
        MULTIPLICITY_FAMILY,
        MULTIPLICITY_SCOPE,
        MULTIPLICITY_METHOD,
        NO_CONFIRMATORY_STATE,
        NO_CONFIRMATORY_STATE,
    )


def _validate_metric_rows(
    *,
    rows: Sequence[Mapping[str, str]],
    model_id: str,
    endpoint_id: str,
    set_id: str,
    points: Mapping[str, str],
) -> tuple[int, int]:
    if len(rows) != EXPECTED_METRICS_PER_BUNDLE or [
        row["metric_id"] for row in rows
    ] != [item[0] for item in METRIC_SPECS]:
        raise HistologyProductionScoreAuditError("metric roster differs")
    computed = 0
    not_applicable = 0
    for row, (metric_id, role, direction) in zip(rows, METRIC_SPECS, strict=True):
        family, scope, method, p_value, q_value = _metric_policy(role)
        native = metric_id in ENDPOINT_METRICS[endpoint_id]
        expected_estimate = points[metric_id] if native else "not_applicable"
        expected_state = "observed" if native else "not_applicable"
        expected_reason = (
            "metric_native_to_endpoint_bundle"
            if native
            else "metric_not_native_to_endpoint_bundle"
        )
        if (
            row.get("task_id") != TASK_ID
            or row.get("model_id") != model_id
            or row.get("endpoint_id") != endpoint_id
            or row.get("prediction_aggregation") != SEED_AGGREGATION
            or row.get("prediction_set_id") != set_id
            or row.get("endpoint_role") != role
            or row.get("analysis_role") != "ensemble_development_endpoint"
            or row.get("applicability_state") != expected_state
            or row.get("applicability_reason") != expected_reason
            or row.get("direction") != direction
            or row.get("estimate") != expected_estimate
            or row.get("biological_resampling_unit") != "participant"
            or row.get("multiplicity_family") != family
            or row.get("multiplicity_scope") != scope
            or row.get("multiplicity_method") != method
            or row.get("p_value") != p_value
            or row.get("bh_adjusted_q_value") != q_value
            or row.get("confirmatory_inference_allowed") != "false"
            or row.get("claim_mode") != CLAIM_MODE
        ):
            raise HistologyProductionScoreAuditError(
                "metric applicability, policy, or point estimate differs"
            )
        try:
            valid = int(row["valid_bootstrap_replicates"])
        except ValueError as error:
            raise HistologyProductionScoreAuditError("bootstrap count differs") from error
        if native:
            computed += 1
            if expected_estimate == "not_estimable":
                if (
                    valid != 0
                    or row["ci95_low"] != "not_estimable"
                    or row["ci95_high"] != "not_estimable"
                ):
                    raise HistologyProductionScoreAuditError(
                        "endpoint-native non-estimable metric state differs"
                    )
            else:
                if not 0 < valid <= BOOTSTRAP_REPLICATES:
                    raise HistologyProductionScoreAuditError("bootstrap count differs")
                try:
                    low = float(row["ci95_low"])
                    high = float(row["ci95_high"])
                except ValueError as error:
                    raise HistologyProductionScoreAuditError(
                        "bootstrap interval differs"
                    ) from error
                if not math.isfinite(low) or not math.isfinite(high) or low > high:
                    raise HistologyProductionScoreAuditError("bootstrap interval differs")
        else:
            not_applicable += 1
            if (
                valid != 0
                or row["ci95_low"] != "not_applicable"
                or row["ci95_high"] != "not_applicable"
            ):
                raise HistologyProductionScoreAuditError(
                    "non-native metric state differs"
                )
    return computed, not_applicable


def _validate_contract(
    contract: Mapping[str, Any], production_artifacts_sha256: str
) -> None:
    production = contract.get("production_campaign", {})
    task = contract.get("task", {})
    inputs = contract.get("evaluator_inputs", {})
    predictions = contract.get("prediction_artifacts", {})
    metrics = contract.get("metrics", {})
    claims = contract.get("claim_boundary", {})
    validation = contract.get("independent_source_validation", {})
    submission = contract.get("submission_gate", {})
    if (
        contract.get("schema_version")
        != "masld-bench-gse267145-histology-production-v4-scoring-contract-v1"
        or contract.get("status") != "bound_root_verified_production_artifacts"
        or contract.get("scoring_revision_id") != SCORING_REVISION_ID
        or production.get("artifacts_sha256") != production_artifacts_sha256
        or production.get("slurm_job_id") != 21089724
        or production.get("path") != "executions/model-training-069-production-v4"
        or production.get("required_status") != "passed_unscored"
        or production.get("aggregate_artifacts_sha256")
        != PRODUCTION_AGGREGATE_ARTIFACTS_SHA256
        or not isinstance(validation.get("path"), str)
        or not validation["path"].startswith("executions/model-check-079-")
        or not isinstance(validation.get("artifacts_sha256"), str)
        or len(validation["artifacts_sha256"]) != 64
        or any(character not in "0123456789abcdef" for character in validation["artifacts_sha256"])
        or validation.get("status")
        != "passed_unscored_source_and_unit_validation"
        or validation.get("production_scores_calculated") is not False
        or task.get("task_id") != TASK_ID
        or task.get("dataset_ids") != list(DATASET_IDS)
        or task.get("split_id") != SPLIT_ID
        or task.get("biological_unit") != "participant"
        or task.get("participants") != EXPECTED_PARTICIPANTS
        or task.get("outer_folds") != 5
        or task.get("model_seeds") != list(reference.SEEDS)
        or task.get("seeds_are_biological_replicates") is not False
        or inputs.get("reference_evaluator_sha256") != REFERENCE_EVALUATOR_SHA256
        or task.get("task_contract_sha256") != TASK_CONTRACT_SHA256
        or task.get("surface_sha256") != SURFACE_SHA256
        or inputs.get("outcomes_artifacts_sha256") != OUTCOMES_ARTIFACTS_SHA256
        or inputs.get("folds_artifacts_sha256") != FOLDS_ARTIFACTS_SHA256
        or inputs.get("model_or_fit_source_import_allowed") is not False
        or predictions.get("standard") != "masld-bench-prediction-bundle-v1"
        or predictions.get("prediction_bundles") != EXPECTED_PREDICTION_BUNDLES
        or predictions.get("predicted_endpoints_per_model") != list(PREDICTED_ENDPOINTS)
        or predictions.get("source_stage5_endpoint_dispositions")
        != EXPECTED_STAGE5_DISPOSITIONS
        or predictions.get("source_stage5_disposition") != SOURCE_STAGE5_DISPOSITION
        or predictions.get("rows_per_set") != EXPECTED_PARTICIPANTS
        or predictions.get("fold_aggregation")
        != "each_participant_appears_once_from_its_held_outer_fold"
        or predictions.get("seed_aggregation") != SEED_AGGREGATION
        or predictions.get("seed_or_fold_results_emitted") is not False
        or predictions.get("frozen_before_outcomes_read") is not True
        or predictions.get("mutation_allowed") is not False
        or metrics.get("bootstrap_replicates") != BOOTSTRAP_REPLICATES
        or metrics.get("bootstrap_seed") != BOOTSTRAP_SEED
        or metrics.get("secondary_multiplicity_family") != MULTIPLICITY_FAMILY
        or metrics.get("secondary_multiplicity_scope") != MULTIPLICITY_SCOPE
        or metrics.get("secondary_multiplicity_method") != MULTIPLICITY_METHOD
        or metrics.get("secondary_p_and_q_state") != NO_CONFIRMATORY_STATE
        or metrics.get("source_stage5")
        != "not_scored_no_separate_registered_prediction_artifact"
        or metrics.get("uncertainty_method")
        != "paired_participant_cluster_bootstrap"
        or metrics.get("primary_inference")
        != "effect_and_interval_without_confirmatory_p_value"
        or metrics.get("real_prediction_bundle_metric_rows") != EXPECTED_METRIC_ROWS
        or metrics.get("endpoint_native_computed_metric_rows")
        != EXPECTED_COMPUTED_METRIC_ROWS
        or metrics.get("non_native_not_applicable_metric_rows")
        != EXPECTED_NOT_APPLICABLE_METRIC_ROWS
        or claims.get("claim_mode") != CLAIM_MODE
        or claims.get("model_ranking_allowed") is not False
        or claims.get("champion_claim_allowed") is not False
        or claims.get("external_transfer_claim_allowed") is not False
        or claims.get("diagnostic_or_prognostic_claim_allowed") is not False
        or claims.get("confirmatory_inference_allowed") is not False
        or submission.get("root_must_verify_production_artifacts_sha256") is not True
        or submission.get("scoring_wrapper_submission_allowed_before_binding") is not False
        or submission.get("independent_compute_node_validation_required") is not True
    ):
        raise HistologyProductionScoreAuditError("scoring contract differs")


def audit(
    *,
    production: Path,
    production_artifacts_sha256: str,
    scores: Path,
    scores_artifacts_sha256: str,
    scorer_source: Path,
    scorer_source_sha256: str,
    outcomes: Path,
    folds: Path,
    reference_evaluator: Path,
    task_contract: Path,
    surface: Path,
    scoring_contract: Path,
    scoring_contract_sha256: str,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise HistologyProductionScoreAuditError(
            f"refusing to overwrite audit output: {output}"
        )
    for path, expected, label in (
        (production / "ARTIFACTS.json", production_artifacts_sha256, "production"),
        (
            production / "aggregate/ARTIFACTS.json",
            PRODUCTION_AGGREGATE_ARTIFACTS_SHA256,
            "production aggregate",
        ),
        (scores / "ARTIFACTS.json", scores_artifacts_sha256, "scores"),
        (scorer_source, scorer_source_sha256, "scorer source"),
        (outcomes / "ARTIFACTS.json", OUTCOMES_ARTIFACTS_SHA256, "outcomes"),
        (folds / "ARTIFACTS.json", FOLDS_ARTIFACTS_SHA256, "folds"),
        (reference_evaluator, REFERENCE_EVALUATOR_SHA256, "reference evaluator"),
        (task_contract, TASK_CONTRACT_SHA256, "task contract"),
        (surface, SURFACE_SHA256, "surface"),
        (scoring_contract, scoring_contract_sha256, "scoring contract"),
    ):
        _verify_file(path, expected, label)
    contract = _read_json(scoring_contract)
    _validate_contract(contract, production_artifacts_sha256)
    source_validation_artifacts_sha256 = contract["independent_source_validation"][
        "artifacts_sha256"
    ]
    production_manifest = verify_frozen_tree(production)
    score_manifest = verify_frozen_tree(scores)
    verify_frozen_tree(outcomes)
    verify_frozen_tree(folds)
    production_metadata = production_manifest.get("metadata", {})
    score_metadata = score_manifest.get("metadata", {})
    if (
        production_metadata.get("artifact_class")
        != "gse267145_histology_production_prediction_bundle"
        or production_metadata.get("status") != "passed_unscored"
        or production_metadata.get("outcomes_read") is not False
        or production_metadata.get("metrics_calculated") is not False
        or score_metadata.get("artifact_class")
        != "gse267145_histology_production_v4_development_scores"
        or score_metadata.get("production_artifacts_sha256")
        != production_artifacts_sha256
        or score_metadata.get("scoring_contract_sha256") != scoring_contract_sha256
        or score_metadata.get("scorer_source_sha256") != scorer_source_sha256
        or score_metadata.get("prediction_bundles") != EXPECTED_PREDICTION_BUNDLES
        or score_metadata.get("source_stage5_endpoint_dispositions")
        != EXPECTED_STAGE5_DISPOSITIONS
        or score_metadata.get("metric_artifacts") != EXPECTED_PREDICTION_BUNDLES
        or score_metadata.get("metric_rows") != EXPECTED_METRIC_ROWS
        or score_metadata.get("endpoint_native_computed_metric_rows")
        != EXPECTED_COMPUTED_METRIC_ROWS
        or score_metadata.get("non_native_not_applicable_metric_rows")
        != EXPECTED_NOT_APPLICABLE_METRIC_ROWS
        or score_metadata.get("biological_resampling_unit") != "participant"
        or score_metadata.get("p_values_calculated") is not False
        or score_metadata.get("bh_adjustment_calculated") is not False
        or score_metadata.get("model_ranking_performed") is not False
        or score_metadata.get("champion_claim_allowed") is not False
        or score_metadata.get("external_claim_allowed") is not False
        or score_metadata.get("diagnostic_or_prognostic_claim_allowed") is not False
    ):
        raise HistologyProductionScoreAuditError("production or score metadata differs")

    receipt = _read_json(scores / "scoring_receipt.json")
    if (
        receipt.get("status") != "passed_development_metrics_no_confirmatory_claim"
        or receipt.get("scoring_revision_id") != SCORING_REVISION_ID
        or receipt.get("production_revision_id") != PRODUCTION_REVISION_ID
        or receipt.get("production_artifacts_sha256") != production_artifacts_sha256
        or receipt.get("scoring_contract_sha256") != scoring_contract_sha256
        or receipt.get("scorer_source_sha256") != scorer_source_sha256
        or receipt.get("aggregate_artifacts_sha256")
        != PRODUCTION_AGGREGATE_ARTIFACTS_SHA256
        or receipt.get("source_validation_artifacts_sha256")
        != source_validation_artifacts_sha256
        or receipt.get("standardized_prediction_bundles")
        != EXPECTED_PREDICTION_BUNDLES
        or receipt.get("source_stage5_endpoint_dispositions")
        != EXPECTED_STAGE5_DISPOSITIONS
        or receipt.get("source_stage5_disposition") != SOURCE_STAGE5_DISPOSITION
        or receipt.get("standardized_metric_artifacts")
        != EXPECTED_PREDICTION_BUNDLES
        or receipt.get("standardized_metric_rows") != EXPECTED_METRIC_ROWS
        or receipt.get("endpoint_native_computed_metric_rows")
        != EXPECTED_COMPUTED_METRIC_ROWS
        or receipt.get("non_native_not_applicable_metric_rows")
        != EXPECTED_NOT_APPLICABLE_METRIC_ROWS
        or receipt.get("predicted_endpoints_per_model") != list(PREDICTED_ENDPOINTS)
        or receipt.get("prediction_aggregation") != SEED_AGGREGATION
        or receipt.get("model_seed_roster") != list(reference.SEEDS)
        or receipt.get("outer_fold_aggregation")
        != "each_participant_once_from_held_fold"
        or receipt.get("seed_or_fold_results_emitted") is not False
        or receipt.get("participants_per_prediction_set") != EXPECTED_PARTICIPANTS
        or receipt.get("bootstrap_replicates") != BOOTSTRAP_REPLICATES
        or receipt.get("bootstrap_seed") != BOOTSTRAP_SEED
        or receipt.get("biological_resampling_unit") != "participant"
        or receipt.get("seeds_are_biological_replicates") is not False
        or receipt.get("prediction_bundles_frozen_before_outcomes_read") is not True
        or receipt.get("five_seed_ensemble_rederived_before_outcomes_read") is not True
        or receipt.get("raw_production_prediction_files_read_after_outcomes") is not False
        or receipt.get("metric_scoring_prediction_source")
        != "55_frozen_endpoint_prediction_bundles"
        or receipt.get("prediction_bundles_mutated") is not False
        or receipt.get("model_or_fit_source_imported") is not False
        or receipt.get("baseline_comparison_performed") is not False
        or receipt.get("model_ranking_performed") is not False
        or receipt.get("p_values_calculated") is not False
        or receipt.get("bh_adjustment_calculated") is not False
        or receipt.get("secondary_multiplicity_family_retained")
        != MULTIPLICITY_FAMILY
        or receipt.get("secondary_multiplicity_scope_retained") != MULTIPLICITY_SCOPE
        or receipt.get("secondary_multiplicity_method_retained")
        != MULTIPLICITY_METHOD
        or receipt.get("secondary_p_and_q_state") != NO_CONFIRMATORY_STATE
        or receipt.get("champion_claim_allowed") is not False
        or receipt.get("external_transfer_claim_allowed") is not False
        or receipt.get("diagnostic_or_prognostic_claim_allowed") is not False
        or receipt.get("confirmatory_inference_allowed") is not False
        or receipt.get("claim_mode") != CLAIM_MODE
    ):
        raise HistologyProductionScoreAuditError("scoring receipt differs")

    participants, outer_folds, outcomes_map = _load_outcomes(outcomes, folds)
    prediction_root = production / "aggregate/predictions"
    ensemble_rows, ensemble_lineage = _verify_five_seed_ensembles(
        prediction_root=prediction_root,
        participants=participants,
        outer_folds=outer_folds,
    )
    prediction_collection = scores / "prediction_bundles"
    disposition_collection = scores / "endpoint_dispositions"
    metric_collection = scores / "metric_artifacts"
    prediction_collection_manifest = verify_frozen_tree(prediction_collection)
    disposition_collection_manifest = verify_frozen_tree(disposition_collection)
    metric_collection_manifest = verify_frozen_tree(metric_collection)
    if (
        sha256_file(prediction_collection / "ARTIFACTS.json")
        != receipt.get("prediction_collection_artifacts_sha256")
        or sha256_file(disposition_collection / "ARTIFACTS.json")
        != receipt.get("endpoint_disposition_collection_artifacts_sha256")
        or sha256_file(metric_collection / "ARTIFACTS.json")
        != receipt.get("metric_collection_artifacts_sha256")
        or prediction_collection_manifest.get("metadata", {}).get("prediction_bundles")
        != EXPECTED_PREDICTION_BUNDLES
        or prediction_collection_manifest.get("metadata", {}).get(
            "seed_or_fold_results_emitted"
        )
        is not False
        or disposition_collection_manifest.get("metadata", {}).get(
            "endpoint_dispositions"
        )
        != EXPECTED_STAGE5_DISPOSITIONS
        or disposition_collection_manifest.get("metadata", {}).get(
            "prediction_bundles_created"
        )
        != 0
        or metric_collection_manifest.get("metadata", {}).get("metric_rows")
        != EXPECTED_METRIC_ROWS
        or metric_collection_manifest.get("metadata", {}).get(
            "endpoint_native_computed_metric_rows"
        )
        != EXPECTED_COMPUTED_METRIC_ROWS
        or metric_collection_manifest.get("metadata", {}).get(
            "non_native_not_applicable_metric_rows"
        )
        != EXPECTED_NOT_APPLICABLE_METRIC_ROWS
        or metric_collection_manifest.get("metadata", {}).get(
            "diagnostic_or_prognostic_claim_allowed"
        )
        is not False
    ):
        raise HistologyProductionScoreAuditError("nested collection receipt differs")
    points_by_model = {
        model_id: _reference_points(
            outcomes_map,
            _frozen_bundle_prediction_arrays(
                prediction_collection=prediction_collection,
                model_id=model_id,
                participants=participants,
                outer_folds=outer_folds,
            ),
        )
        for model_id in reference.MODEL_IDS
    }

    prediction_index_fields, prediction_index_rows = _read_tsv(
        prediction_collection / "prediction_bundle_index.tsv"
    )
    metric_index_fields, metric_index_rows = _read_tsv(
        metric_collection / "metric_artifact_index.tsv"
    )
    disposition_index_fields, disposition_index_rows = _read_tsv(
        disposition_collection / "endpoint_disposition_index.tsv"
    )
    expected_sets = {
        f"{model_id}--endpoint-{endpoint_id}"
        for model_id, endpoint_id in _prediction_sets()
    }
    prediction_index = {row["prediction_set_id"]: row for row in prediction_index_rows}
    metric_index = {row["prediction_set_id"]: row for row in metric_index_rows}
    if (
        len(prediction_index_rows) != EXPECTED_PREDICTION_BUNDLES
        or len(metric_index_rows) != EXPECTED_PREDICTION_BUNDLES
        or set(prediction_index) != expected_sets
        or set(metric_index) != expected_sets
        or "endpoint_id" not in prediction_index_fields
        or "endpoint_id" not in metric_index_fields
        or len(disposition_index_rows) != EXPECTED_STAGE5_DISPOSITIONS
        or "disposition_id" not in disposition_index_fields
    ):
        raise HistologyProductionScoreAuditError("artifact index roster differs")

    disposition_models: set[str] = set()
    for record in disposition_index_rows:
        model_id = record["model_id"]
        target = disposition_collection / model_id
        manifest = verify_frozen_tree(target)
        document_path = target / "endpoint_disposition.json"
        document = _read_json(document_path)
        claimed_id = document.pop("disposition_id", None)
        if (
            model_id not in reference.MODEL_IDS
            or model_id in disposition_models
            or record.get("endpoint_id") != "source_stage5"
            or record.get("disposition") != SOURCE_STAGE5_DISPOSITION
            or record.get("disposition_document_sha256") != sha256_file(document_path)
            or record.get("disposition_artifacts_sha256")
            != sha256_file(target / "ARTIFACTS.json")
            or claimed_id != record.get("disposition_id")
            or claimed_id != canonical_sha256(document)
            or document.get("schema_version")
            != "masld-bench-endpoint-prediction-disposition-v1"
            or document.get("model_id") != model_id
            or document.get("endpoint_id") != "source_stage5"
            or document.get("production_artifacts_sha256")
            != production_artifacts_sha256
            or document.get("prediction_aggregation") != SEED_AGGREGATION
            or document.get("model_seed_roster") != list(reference.SEEDS)
            or document.get("disposition") != SOURCE_STAGE5_DISPOSITION
            or document.get("prediction_bundle_created") is not False
            or document.get("metric_calculated") is not False
            or document.get("diagnostic_or_prognostic_claim_allowed") is not False
            or manifest.get("metadata", {}).get("prediction_bundle_created") is not False
            or manifest.get("metadata", {}).get("metric_calculated") is not False
            or manifest.get("metadata", {}).get(
                "diagnostic_or_prognostic_claim_allowed"
            )
            is not False
            or (target / "prediction_bundle.json").exists()
        ):
            raise HistologyProductionScoreAuditError("source_stage5 disposition differs")
        disposition_models.add(model_id)
    if disposition_models != set(reference.MODEL_IDS):
        raise HistologyProductionScoreAuditError("source_stage5 disposition roster differs")

    all_metric_rows: list[dict[str, str]] = []
    computed_rows = 0
    not_applicable_rows = 0
    for model_id, endpoint_id in _prediction_sets():
        set_id = f"{model_id}--endpoint-{endpoint_id}"
        prediction_record = prediction_index[set_id]
        metric_record = metric_index[set_id]
        bundle_root = prediction_collection / set_id
        bundle_manifest = verify_frozen_tree(bundle_root)
        bundle_document = bundle_root / "prediction_bundle.json"
        bundle = PredictionBundle.load_json(bundle_document)
        bundle.validate_artifacts(bundle_root)
        standardized_fields = STANDARDIZED_BASE_FIELDS + ENDPOINT_VALUE_FIELDS[endpoint_id]
        table_schema = canonical_sha256(
            {"format": "tsv", "fields": list(standardized_fields)}
        )
        lineage = ensemble_lineage[model_id]
        join_key = canonical_sha256(
            {
                "task_id": TASK_ID,
                "dataset_ids": list(DATASET_IDS),
                "split_id": SPLIT_ID,
                "row_id_field": "prediction_row_id",
                "unit_id_field": "participant_id",
                "unit_id_namespace": UNIT_NAMESPACE,
                "biological_unit": "participant",
            }
        )
        bundle_identity = {
            "run_id": production_artifacts_sha256,
            "task_id": TASK_ID,
            "model_id": model_id,
            "endpoint_id": endpoint_id,
            "prediction_aggregation": SEED_AGGREGATION,
            "dataset_ids": list(DATASET_IDS),
            "split_id": SPLIT_ID,
            "standardized_table_sha256": bundle.standardized_table.sha256,
            "row_ids_sha256": bundle.row_ids.sha256,
            "n_predictions": EXPECTED_PARTICIPANTS,
            "source_join_key_sha256": join_key,
        }
        if (
            prediction_record.get("model_id") != model_id
            or prediction_record.get("endpoint_id") != endpoint_id
            or prediction_record.get("prediction_aggregation") != SEED_AGGREGATION
            or prediction_record.get("source_filename") != f"{model_id}.tsv"
            or prediction_record.get("source_prediction_sha256")
            != lineage["ensemble_prediction_sha256"]
            or prediction_record.get("bundle_id") != bundle.bundle_id
            or prediction_record.get("bundle_document_sha256")
            != sha256_file(bundle_document)
            or prediction_record.get("bundle_artifacts_sha256")
            != sha256_file(bundle_root / "ARTIFACTS.json")
            or bundle_manifest.get("metadata", {}).get("outcomes_read") is not False
            or bundle_manifest.get("metadata", {}).get("metrics_calculated") is not False
            or bundle.run_id != production_artifacts_sha256
            or bundle.bundle_id != canonical_sha256(bundle_identity)
            or bundle.task_id != TASK_ID
            or bundle.model_id != model_id
            or tuple(bundle.dataset_ids) != DATASET_IDS
            or bundle.split_id != SPLIT_ID
            or bundle.n_predictions != EXPECTED_PARTICIPANTS
            or bundle.row_id_field != "prediction_row_id"
            or bundle.unit_id_field != "participant_id"
            or bundle.unit_id_namespace != UNIT_NAMESPACE
            or bundle.biological_unit != "participant"
            or bundle.table_schema_sha256 != table_schema
            or bundle.source_join_key_sha256 != join_key
            or bundle.missing_state is not MissingState.OBSERVED
            or bundle.metadata.get("endpoint_id") != endpoint_id
            or bundle.metadata.get("prediction_aggregation") != SEED_AGGREGATION
            or tuple(bundle.metadata.get("model_seed_roster", ())) != reference.SEEDS
            or bundle.metadata.get("outer_fold_aggregation")
            != "each_participant_once_from_held_fold"
            or bundle.metadata.get("source_seed_prediction_sha256")
            != lineage["seed_prediction_sha256"]
            or bundle.metadata.get("five_seed_membership_sha256")
            != lineage["five_seed_membership_sha256"]
            or bundle.metadata.get("prediction_mutated") is not False
            or bundle.metadata.get("frozen_before_outcomes_read") is not True
            or bundle.metadata.get("seeds_are_biological_replicates") is not False
            or bundle.metadata.get("champion_claim_allowed") is not False
            or bundle.metadata.get("external_claim_allowed") is not False
            or bundle.metadata.get("diagnostic_or_prognostic_claim_allowed") is not False
        ):
            raise HistologyProductionScoreAuditError("PredictionBundle contract differs")
        raw_rows = ensemble_rows[model_id]
        table_fields, table_rows = _read_tsv(
            bundle.standardized_table.validate(bundle_root, require_relative=True)
        )
        row_fields, row_rows = _read_tsv(
            bundle.row_ids.validate(bundle_root, require_relative=True)
        )
        if (
            table_fields != standardized_fields
            or row_fields != ROW_ID_FIELDS
            or len(table_rows) != EXPECTED_PARTICIPANTS
            or len(row_rows) != EXPECTED_PARTICIPANTS
        ):
            raise HistologyProductionScoreAuditError("standardized prediction schema differs")
        for raw, standardized, row_id_record, participant, outer_fold in zip(
            raw_rows, table_rows, row_rows, participants, outer_folds, strict=True
        ):
            expected_row_id = _prediction_row_id(
                production_artifacts_sha256,
                model_id,
                endpoint_id,
                participant,
                outer_fold,
            )
            expected_row = {
                "prediction_row_id": expected_row_id,
                "participant_id": participant,
                "outer_fold": outer_fold,
                "endpoint_id": endpoint_id,
                "prediction_state": "observed",
                **{field: raw[field] for field in ENDPOINT_VALUE_FIELDS[endpoint_id]},
            }
            if (
                standardized != expected_row
                or row_id_record
                != {"prediction_row_id": expected_row_id, "participant_id": participant}
            ):
                raise HistologyProductionScoreAuditError("standardized prediction row differs")

        metric_root = metric_collection / set_id
        metric_manifest = verify_frozen_tree(metric_root)
        metric_document_path = metric_root / "development_metric_artifact.json"
        metric_document = _read_json(metric_document_path)
        claimed_id = metric_document.pop("metric_artifact_id", None)
        if (
            metric_record.get("model_id") != model_id
            or metric_record.get("endpoint_id") != endpoint_id
            or metric_record.get("prediction_aggregation") != SEED_AGGREGATION
            or metric_record.get("metric_artifact_id") != claimed_id
            or metric_record.get("metric_artifact_document_sha256")
            != sha256_file(metric_document_path)
            or metric_record.get("metric_artifact_artifacts_sha256")
            != sha256_file(metric_root / "ARTIFACTS.json")
            or metric_record.get("prediction_bundle_id") != bundle.bundle_id
            or metric_record.get("prediction_bundle_artifacts_sha256")
            != prediction_record["bundle_artifacts_sha256"]
            or claimed_id != canonical_sha256(metric_document)
            or metric_document.get("scoring_revision_id") != SCORING_REVISION_ID
            or metric_document.get("model_id") != model_id
            or metric_document.get("endpoint_id") != endpoint_id
            or metric_document.get("prediction_aggregation") != SEED_AGGREGATION
            or metric_document.get("model_seed_roster") != list(reference.SEEDS)
            or metric_document.get("outer_fold_aggregation")
            != "each_participant_once_from_held_fold"
            or metric_document.get("production_artifacts_sha256")
            != production_artifacts_sha256
            or metric_document.get("prediction_bundle_id") != bundle.bundle_id
            or metric_document.get("prediction_bundle_artifacts_sha256")
            != prediction_record["bundle_artifacts_sha256"]
            or metric_document.get("scoring_contract_sha256")
            != scoring_contract_sha256
            or metric_document.get("scorer_source_sha256") != scorer_source_sha256
            or metric_document.get("outcomes_artifacts_sha256")
            != OUTCOMES_ARTIFACTS_SHA256
            or metric_document.get("folds_artifacts_sha256") != FOLDS_ARTIFACTS_SHA256
            or metric_document.get("secondary_p_and_q_state") != NO_CONFIRMATORY_STATE
            or metric_document.get("secondary_multiplicity_scope") != MULTIPLICITY_SCOPE
            or metric_document.get("model_ranking_performed") is not False
            or metric_document.get("champion_claim_allowed") is not False
            or metric_document.get("external_claim_allowed") is not False
            or metric_document.get("diagnostic_or_prognostic_claim_allowed") is not False
            or metric_document.get("confirmatory_inference_allowed") is not False
            or metric_manifest.get("metadata", {}).get("metric_rows")
            != EXPECTED_METRICS_PER_BUNDLE
            or metric_manifest.get("metadata", {}).get("p_values_calculated") is not False
            or metric_manifest.get("metadata", {}).get("bh_adjustment_calculated") is not False
            or metric_manifest.get("metadata", {}).get("champion_claim_allowed") is not False
            or metric_manifest.get("metadata", {}).get("external_claim_allowed") is not False
            or metric_manifest.get("metadata", {}).get(
                "diagnostic_or_prognostic_claim_allowed"
            )
            is not False
        ):
            raise HistologyProductionScoreAuditError("metric artifact identity differs")
        metric_ref = ArtifactRef.from_dict(metric_document["metrics_table"])
        metric_fields, metric_rows = _read_tsv(
            metric_ref.validate(metric_root, require_relative=True)
        )
        if metric_fields != METRIC_FIELDS:
            raise HistologyProductionScoreAuditError("metric table schema differs")
        native, non_native = _validate_metric_rows(
            rows=metric_rows,
            model_id=model_id,
            endpoint_id=endpoint_id,
            set_id=set_id,
            points=points_by_model[model_id],
        )
        computed_rows += native
        not_applicable_rows += non_native
        all_metric_rows.extend(metric_rows)

    aggregate_fields, aggregate_rows = _read_tsv(scores / "standardized_metrics.tsv")
    sex_fields, sex_rows = _read_tsv(scores / "sex_error_audit.tsv")
    expected_stage3_sets = {
        f"{model_id}--endpoint-stage3" for model_id in reference.MODEL_IDS
    }
    if (
        aggregate_fields != METRIC_FIELDS
        or aggregate_rows != all_metric_rows
        or len(aggregate_rows) != EXPECTED_METRIC_ROWS
        or computed_rows != EXPECTED_COMPUTED_METRIC_ROWS
        or not_applicable_rows != EXPECTED_NOT_APPLICABLE_METRIC_ROWS
        or len(sex_rows) != EXPECTED_SEX_AUDIT_ROWS
        or "between_sex_global_gap_claim_allowed" not in sex_fields
        or {row["prediction_set_id"] for row in sex_rows} != expected_stage3_sets
        or any(row.get("endpoint_id") != "stage3" for row in sex_rows)
        or any(row.get("prediction_aggregation") != SEED_AGGREGATION for row in sex_rows)
        or any(row["between_sex_global_gap_claim_allowed"] != "false" for row in sex_rows)
        or any(
            row["stage3_macro_f1"] != "not_applicable_missing_NOR_and_NAFL"
            for row in sex_rows
            if row["recorded_sex"] == "M"
        )
    ):
        raise HistologyProductionScoreAuditError("aggregate metric or sex audit differs")

    output.mkdir(parents=True)
    audit_receipt = {
        "schema_version": "masld-bench-gse267145-histology-production-v4-score-audit-v1",
        "status": "passed_independent_development_score_audit",
        "production_artifacts_sha256": production_artifacts_sha256,
        "scores_artifacts_sha256": scores_artifacts_sha256,
        "scoring_contract_sha256": scoring_contract_sha256,
        "scorer_source_sha256": scorer_source_sha256,
        "source_validation_artifacts_sha256": source_validation_artifacts_sha256,
        "auditor_source_sha256": sha256_file(Path(__file__).resolve()),
        "prediction_bundles_verified": EXPECTED_PREDICTION_BUNDLES,
        "source_stage5_endpoint_dispositions_verified": EXPECTED_STAGE5_DISPOSITIONS,
        "metric_artifacts_verified": EXPECTED_PREDICTION_BUNDLES,
        "metric_rows_verified": EXPECTED_METRIC_ROWS,
        "endpoint_native_computed_metric_rows_verified": EXPECTED_COMPUTED_METRIC_ROWS,
        "non_native_not_applicable_metric_rows_verified": (
            EXPECTED_NOT_APPLICABLE_METRIC_ROWS
        ),
        "ensemble_membership_rederived_from_five_seed_predictions": True,
        "participant_held_fold_axis_rederived": True,
        "point_estimates_rederived_with_frozen_reference_evaluator": True,
        "point_estimate_prediction_source": "55_frozen_endpoint_prediction_bundles",
        "raw_production_predictions_used_only_for_pre_outcome_lineage_audit": True,
        "bootstrap_contract_verified": True,
        "bootstrap_interval_values_rederived": False,
        "bootstrap_interval_state": "structural_and_range_audit_only",
        "prediction_aggregation": SEED_AGGREGATION,
        "outer_fold_aggregation": "each_participant_once_from_held_fold",
        "seed_or_fold_results_emitted": False,
        "biological_resampling_unit": "participant",
        "seeds_are_biological_replicates": False,
        "secondary_multiplicity_family_retained": MULTIPLICITY_FAMILY,
        "secondary_multiplicity_method_retained": MULTIPLICITY_METHOD,
        "secondary_p_and_q_state": NO_CONFIRMATORY_STATE,
        "p_values_calculated": False,
        "bh_adjustment_calculated": False,
        "baseline_comparison_performed": False,
        "model_ranking_performed": False,
        "champion_claim_allowed": False,
        "external_transfer_claim_allowed": False,
        "diagnostic_or_prognostic_claim_allowed": False,
        "confirmatory_inference_allowed": False,
        "claim_mode": CLAIM_MODE,
    }
    write_json_exclusive(output / "audit_receipt.json", audit_receipt, mode=0o440)
    freeze_tree(
        output,
        {
            "artifact_class": "gse267145_histology_production_v4_score_audit",
            "production_artifacts_sha256": production_artifacts_sha256,
            "scores_artifacts_sha256": scores_artifacts_sha256,
            "scoring_contract_sha256": scoring_contract_sha256,
            "scorer_source_sha256": scorer_source_sha256,
            "prediction_bundles_verified": EXPECTED_PREDICTION_BUNDLES,
            "source_stage5_endpoint_dispositions_verified": (
                EXPECTED_STAGE5_DISPOSITIONS
            ),
            "metric_artifacts_verified": EXPECTED_PREDICTION_BUNDLES,
            "metric_rows_verified": EXPECTED_METRIC_ROWS,
            "endpoint_native_computed_metric_rows_verified": (
                EXPECTED_COMPUTED_METRIC_ROWS
            ),
            "non_native_not_applicable_metric_rows_verified": (
                EXPECTED_NOT_APPLICABLE_METRIC_ROWS
            ),
            "biological_resampling_unit": "participant",
            "p_values_calculated": False,
            "bh_adjustment_calculated": False,
            "model_ranking_performed": False,
            "champion_claim_allowed": False,
            "external_claim_allowed": False,
            "diagnostic_or_prognostic_claim_allowed": False,
            "status": "passed_independent_audit",
        },
    )
    verify_frozen_tree(output)
    return audit_receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production", required=True, type=Path)
    parser.add_argument("--production-artifacts-sha256", required=True)
    parser.add_argument("--scores", required=True, type=Path)
    parser.add_argument("--scores-artifacts-sha256", required=True)
    parser.add_argument("--scorer-source", required=True, type=Path)
    parser.add_argument("--scorer-source-sha256", required=True)
    parser.add_argument("--outcomes", required=True, type=Path)
    parser.add_argument("--folds", required=True, type=Path)
    parser.add_argument("--reference-evaluator", required=True, type=Path)
    parser.add_argument("--task-contract", required=True, type=Path)
    parser.add_argument("--surface", required=True, type=Path)
    parser.add_argument("--scoring-contract", required=True, type=Path)
    parser.add_argument("--scoring-contract-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    result = audit(
        production=arguments.production,
        production_artifacts_sha256=arguments.production_artifacts_sha256,
        scores=arguments.scores,
        scores_artifacts_sha256=arguments.scores_artifacts_sha256,
        scorer_source=arguments.scorer_source,
        scorer_source_sha256=arguments.scorer_source_sha256,
        outcomes=arguments.outcomes,
        folds=arguments.folds,
        reference_evaluator=arguments.reference_evaluator,
        task_contract=arguments.task_contract,
        surface=arguments.surface,
        scoring_contract=arguments.scoring_contract,
        scoring_contract_sha256=arguments.scoring_contract_sha256,
        output=arguments.output,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
