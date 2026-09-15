#!/usr/bin/env python3
"""Freeze and score the development-only GSE267145 v4 prediction campaign."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
import math
from pathlib import Path
import shutil
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
from scipy.stats import rankdata, spearmanr

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

PREDICTED_ENDPOINTS = (
    "stage3",
    "fibrosis_group3",
    "fibrosis_cumulative",
    "fibrosis_regression",
    "nash_crn_component_sum",
)
SOURCE_STAGE5_DISPOSITION = "not_applicable_no_registered_prediction"
SEED_AGGREGATION = "mean_predictions_across_all_five_fixed_model_seeds_before_scoring"
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
SEX_AUDIT_FIELDS = (
    "model_id",
    "endpoint_id",
    "prediction_aggregation",
    "prediction_set_id",
    "recorded_sex",
    "participants",
    "stage3_census",
    "stage3_brier",
    "stage3_macro_f1",
    "nash_recall",
    "mean_nash_true_class_probability",
    "between_sex_global_gap_claim_allowed",
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


class HistologyProductionScoringError(RuntimeError):
    """Raised when production scoring crosses a frozen requirement boundary."""


def _read_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise HistologyProductionScoringError(f"invalid JSON: {path}") from error
    if not isinstance(value, dict):
        raise HistologyProductionScoringError(f"JSON is not an object: {path}")
    return value


def _read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    try:
        with path.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if reader.fieldnames is None:
                raise HistologyProductionScoringError(f"TSV has no header: {path}")
            return tuple(reader.fieldnames), [dict(row) for row in reader]
    except OSError as error:
        raise HistologyProductionScoringError(f"cannot read TSV: {path}") from error


def _write_tsv(
    path: Path, fields: Sequence[str], rows: Iterable[Mapping[str, Any]]
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def _production_prediction_files() -> set[str]:
    result: set[str] = set()
    for model_id in reference.MODEL_IDS:
        result.add(f"{model_id}.tsv")
        result.update(f"{model_id}--seed-{seed}.tsv" for seed in reference.SEEDS)
    return result


def _endpoint_bundles() -> list[tuple[str, str, str]]:
    return [
        (model_id, endpoint_id, f"{model_id}.tsv")
        for model_id in reference.MODEL_IDS
        for endpoint_id in PREDICTED_ENDPOINTS
    ]


def _metric_policy(endpoint_role: str) -> tuple[str, str, str, str, str]:
    if endpoint_role == "primary":
        return (
            "not_applicable_primary_endpoint",
            "not_applicable_primary_endpoint",
            "primary_effect_and_interval_without_confirmatory_p_value",
            "not_calculated_primary_development_only",
            "not_applicable_primary_development_only",
        )
    if endpoint_role == "calibration":
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


def _validate_contract(
    contract: Mapping[str, Any], *, production_artifacts_sha256: str
) -> None:
    task = contract.get("task", {})
    inputs = contract.get("evaluator_inputs", {})
    prediction = contract.get("prediction_artifacts", {})
    metrics = contract.get("metrics", {})
    claims = contract.get("claim_boundary", {})
    production = contract.get("production_campaign", {})
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
        or prediction.get("standard") != "masld-bench-prediction-bundle-v1"
        or prediction.get("prediction_bundles") != 55
        or prediction.get("predicted_endpoints_per_model") != list(PREDICTED_ENDPOINTS)
        or prediction.get("source_stage5_endpoint_dispositions") != 11
        or prediction.get("source_stage5_disposition") != SOURCE_STAGE5_DISPOSITION
        or prediction.get("rows_per_set") != EXPECTED_PARTICIPANTS
        or prediction.get("fold_aggregation")
        != "each_participant_appears_once_from_its_held_outer_fold"
        or prediction.get("seed_aggregation")
        != "mean_predictions_across_all_five_fixed_model_seeds_before_scoring"
        or prediction.get("seed_or_fold_results_emitted") is not False
        or prediction.get("frozen_before_outcomes_read") is not True
        or prediction.get("mutation_allowed") is not False
        or metrics.get("primary") != ["stage3_macro_f1"]
        or metrics.get("calibration") != ["stage3_multiclass_brier"]
        or metrics.get("secondary")
        != [metric_id for metric_id, role, _ in METRIC_SPECS if role == "secondary"]
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
        or metrics.get("rectangular_metric_rows_per_real_prediction_bundle") != 9
        or metrics.get("real_prediction_bundle_metric_rows") != 495
        or metrics.get("endpoint_native_computed_metric_rows") != 99
        or metrics.get("non_native_not_applicable_metric_rows") != 396
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
        raise HistologyProductionScoringError("scoring contract differs")


def _verify_file(path: Path, expected: str, label: str) -> None:
    if not path.is_file() or path.is_symlink() or sha256_file(path) != expected:
        raise HistologyProductionScoringError(f"{label} SHA-256 differs")


def _verify_production(production: Path, expected_sha256: str) -> dict[str, Any]:
    _verify_file(production / "ARTIFACTS.json", expected_sha256, "production campaign")
    manifest = verify_frozen_tree(production)
    metadata = manifest.get("metadata", {})
    receipt = _read_json(production / "campaign_receipt.json")
    specification = _read_json(production / "campaign_spec/campaign_spec.json")
    aggregate = production / "aggregate"
    _verify_file(
        aggregate / "ARTIFACTS.json",
        PRODUCTION_AGGREGATE_ARTIFACTS_SHA256,
        "production aggregate",
    )
    aggregate_manifest = verify_frozen_tree(aggregate)
    aggregate_receipt = _read_json(aggregate / "receipt.json")
    aggregate_sha256 = sha256_file(aggregate / "ARTIFACTS.json")
    if (
        metadata.get("artifact_class")
        != "gse267145_histology_production_prediction_bundle"
        or metadata.get("logical_outer_seed_units") != 25
        or metadata.get("model_family_unit_fits") != 275
        or metadata.get("prediction_files") != 66
        or metadata.get("outcomes_read") is not False
        or metadata.get("metrics_calculated") is not False
        or metadata.get("status") != "passed_unscored"
        or receipt.get("schema_version")
        != "masld-bench-gse267145-production-bundle-v1"
        or receipt.get("status") != "passed_predictions_unscored"
        or receipt.get("campaign_id") != PRODUCTION_REVISION_ID
        or receipt.get("aggregate_artifacts_sha256") != aggregate_sha256
        or aggregate_sha256 != PRODUCTION_AGGREGATE_ARTIFACTS_SHA256
        or receipt.get("outcomes_read") is not False
        or receipt.get("metrics_calculated") is not False
        or receipt.get("scorer_called") is not False
        or receipt.get("seeds_are_biological_replicates") is not False
        or specification.get("revision_id") != PRODUCTION_REVISION_ID
        or specification.get("outcomes_read") is not False
        or specification.get("metrics_calculated") is not False
        or specification.get("scorer_called") is not False
        or aggregate_manifest.get("metadata", {}).get("outcomes_read") is not False
        or aggregate_manifest.get("metadata", {}).get("metrics_calculated") is not False
        or aggregate_receipt.get("participants") != EXPECTED_PARTICIPANTS
        or aggregate_receipt.get("model_ids") != list(reference.MODEL_IDS)
        or aggregate_receipt.get("model_seeds") != list(reference.SEEDS)
        or aggregate_receipt.get("prediction_files") != 66
        or aggregate_receipt.get("outcomes_read") is not False
        or aggregate_receipt.get("metrics_calculated") is not False
        or aggregate_receipt.get("scorer_called") is not False
    ):
        raise HistologyProductionScoringError("production campaign contract differs")
    expected_files = _production_prediction_files()
    prediction_root = aggregate / "predictions"
    if {path.name for path in prediction_root.glob("*.tsv")} != expected_files:
        raise HistologyProductionScoringError("production prediction roster differs")
    return {
        "aggregate_path": aggregate,
        "aggregate_artifacts_sha256": aggregate_sha256,
        "prediction_root": prediction_root,
    }


def _participant_axis(folds: Path) -> tuple[list[str], list[str]]:
    fields, rows = _read_tsv(folds / "participant_outer_folds.tsv")
    if fields != ("participant_id", "outer_fold") or len(rows) != EXPECTED_PARTICIPANTS:
        raise HistologyProductionScoringError("participant fold axis differs")
    participants = [row["participant_id"] for row in rows]
    outer_folds = [row["outer_fold"] for row in rows]
    if len(set(participants)) != EXPECTED_PARTICIPANTS or Counter(outer_folds) != {
        "0": 21,
        "1": 21,
        "2": 21,
        "3": 19,
        "4": 17,
    }:
        raise HistologyProductionScoringError("participant fold census differs")
    return participants, outer_folds


def _prediction_row_id(
    *,
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


def _verify_five_seed_ensembles(
    *,
    source_root: Path,
    participants: Sequence[str],
    outer_folds: Sequence[str],
) -> dict[str, dict[str, Any]]:
    """Prove each scored ensemble is the exact frozen five-seed mean."""
    lineage: dict[str, dict[str, Any]] = {}
    for model_id in reference.MODEL_IDS:
        seed_rows: list[list[dict[str, str]]] = []
        seed_hashes: dict[str, str] = {}
        for seed in reference.SEEDS:
            source = source_root / f"{model_id}--seed-{seed}.tsv"
            fields, rows = _read_tsv(source)
            if (
                fields != reference.PREDICTION_FIELDS
                or len(rows) != EXPECTED_PARTICIPANTS
                or [row["participant_id"] for row in rows] != list(participants)
                or [row["outer_fold"] for row in rows] != list(outer_folds)
            ):
                raise HistologyProductionScoringError(
                    f"seed prediction axis differs: {source.name}"
                )
            _load_predictions(rows)
            seed_rows.append(rows)
            seed_hashes[str(seed)] = sha256_file(source)
        ensemble = source_root / f"{model_id}.tsv"
        ensemble_fields, ensemble_rows = _read_tsv(ensemble)
        if (
            ensemble_fields != reference.PREDICTION_FIELDS
            or len(ensemble_rows) != EXPECTED_PARTICIPANTS
        ):
            raise HistologyProductionScoringError(
                f"ensemble prediction schema differs: {ensemble.name}"
            )
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
                raise HistologyProductionScoringError(
                    f"ensemble is not the exact five-seed mean: {ensemble.name}"
                )
        lineage[model_id] = {
            "ensemble_prediction_sha256": sha256_file(ensemble),
            "seed_prediction_sha256": seed_hashes,
            "five_seed_membership_sha256": canonical_sha256(seed_hashes),
        }
    return lineage


def _freeze_prediction_bundles(
    *,
    source_root: Path,
    output: Path,
    production_sha256: str,
    participants: Sequence[str],
    outer_folds: Sequence[str],
) -> tuple[dict[str, dict[str, Any]], str]:
    output.mkdir()
    records: dict[str, dict[str, Any]] = {}
    index_rows: list[dict[str, Any]] = []
    ensemble_lineage = _verify_five_seed_ensembles(
        source_root=source_root,
        participants=participants,
        outer_folds=outer_folds,
    )
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
    for model_id, endpoint_id, filename in _endpoint_bundles():
        set_id = f"{model_id}--endpoint-{endpoint_id}"
        target = output / set_id
        target.mkdir()
        source = source_root / filename
        fields, rows = _read_tsv(source)
        if (
            fields != reference.PREDICTION_FIELDS
            or len(rows) != EXPECTED_PARTICIPANTS
            or [row["participant_id"] for row in rows] != list(participants)
            or [row["outer_fold"] for row in rows] != list(outer_folds)
        ):
            raise HistologyProductionScoringError(
                f"production prediction axis differs: {filename}"
            )
        standardized_rows: list[dict[str, str]] = []
        row_id_rows: list[dict[str, str]] = []
        standardized_fields = STANDARDIZED_BASE_FIELDS + ENDPOINT_VALUE_FIELDS[
            endpoint_id
        ]
        for row in rows:
            row_id = _prediction_row_id(
                production_sha256=production_sha256,
                model_id=model_id,
                endpoint_id=endpoint_id,
                participant_id=row["participant_id"],
                outer_fold=row["outer_fold"],
            )
            standardized_rows.append(
                {
                    "prediction_row_id": row_id,
                    "participant_id": row["participant_id"],
                    "outer_fold": row["outer_fold"],
                    "endpoint_id": endpoint_id,
                    "prediction_state": "observed",
                    **{field: row[field] for field in ENDPOINT_VALUE_FIELDS[endpoint_id]},
                }
            )
            row_id_rows.append(
                {"prediction_row_id": row_id, "participant_id": row["participant_id"]}
            )
        table = target / "predictions.tsv"
        row_ids = target / "row_ids.tsv"
        _write_tsv(table, standardized_fields, standardized_rows)
        _write_tsv(row_ids, ROW_ID_FIELDS, row_id_rows)
        table_ref = ArtifactRef.from_path(
            table,
            relative_to=target,
            media_type="text/tab-separated-values",
            role=f"standardized_prediction_table:{TASK_ID}",
        )
        row_ref = ArtifactRef.from_path(
            row_ids,
            relative_to=target,
            media_type="text/tab-separated-values",
            role=f"prediction_row_ids:{TASK_ID}",
        )
        identity = {
            "run_id": production_sha256,
            "task_id": TASK_ID,
            "model_id": model_id,
            "endpoint_id": endpoint_id,
            "prediction_aggregation": SEED_AGGREGATION,
            "dataset_ids": list(DATASET_IDS),
            "split_id": SPLIT_ID,
            "standardized_table_sha256": table_ref.sha256,
            "row_ids_sha256": row_ref.sha256,
            "n_predictions": EXPECTED_PARTICIPANTS,
            "source_join_key_sha256": join_key,
        }
        bundle = PredictionBundle(
            schema_version="masld-bench-prediction-bundle-v1",
            bundle_id=canonical_sha256(identity),
            run_id=production_sha256,
            task_id=TASK_ID,
            model_id=model_id,
            dataset_ids=DATASET_IDS,
            split_id=SPLIT_ID,
            artifacts=(table_ref, row_ref),
            standardized_table=table_ref,
            row_ids=row_ref,
            n_predictions=EXPECTED_PARTICIPANTS,
            row_id_field="prediction_row_id",
            unit_id_field="participant_id",
            unit_id_namespace=UNIT_NAMESPACE,
            biological_unit="participant",
            table_schema_sha256=canonical_sha256(
                {"format": "tsv", "fields": list(standardized_fields)}
            ),
            source_join_key_sha256=join_key,
            format_version="gse267145-histology-endpoint-prediction-tsv-v1",
            missing_state=MissingState.OBSERVED,
            metadata={
                "prediction_set_id": set_id,
                "endpoint_id": endpoint_id,
                "prediction_role": "five_seed_ensemble_development_endpoint",
                "prediction_aggregation": SEED_AGGREGATION,
                "model_seed_roster": list(reference.SEEDS),
                "outer_fold_aggregation": "each_participant_once_from_held_fold",
                "source_prediction_filename": filename,
                "source_prediction_sha256": sha256_file(source),
                "source_seed_prediction_sha256": ensemble_lineage[model_id][
                    "seed_prediction_sha256"
                ],
                "five_seed_membership_sha256": ensemble_lineage[model_id][
                    "five_seed_membership_sha256"
                ],
                "production_artifacts_sha256": production_sha256,
                "frozen_before_outcomes_read": True,
                "prediction_mutated": False,
                "seeds_are_biological_replicates": False,
                "champion_claim_allowed": False,
                "external_claim_allowed": False,
                "diagnostic_or_prognostic_claim_allowed": False,
            },
        )
        document = target / "prediction_bundle.json"
        write_json_exclusive(document, bundle.to_dict(), mode=0o440)
        loaded = PredictionBundle.load_json(document)
        loaded.validate_artifacts(target)
        artifacts_sha256 = freeze_tree(
            target,
            {
                "artifact_class": "gse267145_histology_standardized_prediction_bundle",
                "bundle_id": bundle.bundle_id,
                "model_id": model_id,
                "endpoint_id": endpoint_id,
                "prediction_aggregation": SEED_AGGREGATION,
                "participants": EXPECTED_PARTICIPANTS,
                "biological_unit": "participant",
                "outcomes_read": False,
                "metrics_calculated": False,
                "status": "passed_unscored",
            },
        )
        record = {
            "prediction_set_id": set_id,
            "model_id": model_id,
            "endpoint_id": endpoint_id,
            "prediction_aggregation": SEED_AGGREGATION,
            "source_filename": filename,
            "source_prediction_sha256": sha256_file(source),
            "bundle_path": f"{set_id}/prediction_bundle.json",
            "bundle_id": bundle.bundle_id,
            "bundle_document_sha256": sha256_file(document),
            "bundle_artifacts_sha256": artifacts_sha256,
        }
        records[set_id] = record
        index_rows.append(record)
    if len(records) != 55 or len(index_rows) != 55:
        raise HistologyProductionScoringError("PredictionBundle count differs")
    _write_tsv(
        output / "prediction_bundle_index.tsv",
        (
            "prediction_set_id",
            "model_id",
            "endpoint_id",
            "prediction_aggregation",
            "source_filename",
            "source_prediction_sha256",
            "bundle_path",
            "bundle_id",
            "bundle_document_sha256",
            "bundle_artifacts_sha256",
        ),
        index_rows,
    )
    artifacts_sha256 = freeze_tree(
        output,
        {
            "artifact_class": "gse267145_histology_prediction_bundle_collection",
            "prediction_bundles": len(records),
            "models": len(reference.MODEL_IDS),
            "predicted_endpoints_per_model": list(PREDICTED_ENDPOINTS),
            "participants_per_bundle": EXPECTED_PARTICIPANTS,
            "prediction_aggregation": SEED_AGGREGATION,
            "seed_or_fold_results_emitted": False,
            "production_artifacts_sha256": production_sha256,
            "outcomes_read": False,
            "metrics_calculated": False,
            "status": "passed_unscored",
        },
    )
    return records, artifacts_sha256


def _freeze_source_stage5_dispositions(
    *, output: Path, production_sha256: str
) -> tuple[list[dict[str, str]], str]:
    output.mkdir()
    index_rows: list[dict[str, str]] = []
    for model_id in reference.MODEL_IDS:
        target = output / model_id
        target.mkdir()
        identity = {
            "schema_version": "masld-bench-endpoint-prediction-disposition-v1",
            "task_id": TASK_ID,
            "model_id": model_id,
            "endpoint_id": "source_stage5",
            "production_artifacts_sha256": production_sha256,
            "prediction_aggregation": SEED_AGGREGATION,
            "model_seed_roster": list(reference.SEEDS),
            "disposition": SOURCE_STAGE5_DISPOSITION,
            "prediction_bundle_created": False,
            "metric_calculated": False,
            "reason": "source_stage5 requires a separate registered prediction artifact and none exists",
            "champion_claim_allowed": False,
            "external_claim_allowed": False,
            "diagnostic_or_prognostic_claim_allowed": False,
            "claim_mode": CLAIM_MODE,
        }
        payload = {"disposition_id": canonical_sha256(identity), **identity}
        document = target / "endpoint_disposition.json"
        write_json_exclusive(document, payload, mode=0o440)
        artifacts_sha256 = freeze_tree(
            target,
            {
                "artifact_class": "gse267145_histology_endpoint_prediction_disposition",
                "disposition_id": payload["disposition_id"],
                "model_id": model_id,
                "endpoint_id": "source_stage5",
                "disposition": SOURCE_STAGE5_DISPOSITION,
                "prediction_bundle_created": False,
                "metric_calculated": False,
                "diagnostic_or_prognostic_claim_allowed": False,
                "status": "not_applicable",
            },
        )
        index_rows.append(
            {
                "model_id": model_id,
                "endpoint_id": "source_stage5",
                "disposition": SOURCE_STAGE5_DISPOSITION,
                "disposition_path": f"{model_id}/endpoint_disposition.json",
                "disposition_id": payload["disposition_id"],
                "disposition_document_sha256": sha256_file(document),
                "disposition_artifacts_sha256": artifacts_sha256,
            }
        )
    if len(index_rows) != 11:
        raise HistologyProductionScoringError("source_stage5 disposition count differs")
    _write_tsv(
        output / "endpoint_disposition_index.tsv",
        (
            "model_id",
            "endpoint_id",
            "disposition",
            "disposition_path",
            "disposition_id",
            "disposition_document_sha256",
            "disposition_artifacts_sha256",
        ),
        index_rows,
    )
    collection_sha256 = freeze_tree(
        output,
        {
            "artifact_class": "gse267145_histology_endpoint_disposition_collection",
            "endpoint_dispositions": len(index_rows),
            "endpoint_id": "source_stage5",
            "disposition": SOURCE_STAGE5_DISPOSITION,
            "prediction_bundles_created": 0,
            "metrics_calculated": False,
            "diagnostic_or_prognostic_claim_allowed": False,
            "production_artifacts_sha256": production_sha256,
            "status": "not_applicable",
        },
    )
    return index_rows, collection_sha256


def _load_outcomes(outcomes: Path, participants: Sequence[str], folds: Sequence[str]) -> dict[str, np.ndarray]:
    fields, rows = _read_tsv(outcomes / "participant_endpoints.tsv")
    required = {
        "participant_id",
        "outer_fold",
        "stage3",
        "nash_crn_component_sum",
        "fibrosis",
        "recorded_sex",
    }
    if (
        not required <= set(fields)
        or len(rows) != EXPECTED_PARTICIPANTS
        or [row["participant_id"] for row in rows] != list(participants)
        or [row["outer_fold"] for row in rows] != list(folds)
    ):
        raise HistologyProductionScoringError("outcome participant axis differs")
    stage = np.asarray([row["stage3"] for row in rows], dtype=str)
    fibrosis = np.asarray([int(row["fibrosis"]) for row in rows], dtype=np.float64)
    sex = np.asarray([row["recorded_sex"] for row in rows], dtype=str)
    if (
        not set(stage) <= set(reference.STAGE3)
        or Counter(sex) != {"F": 85, "M": 14}
        or set(stage[sex == "M"]) != {"NASH"}
        or Counter(fibrosis)[3.0] != 4
    ):
        raise HistologyProductionScoringError("outcome census differs")
    return {
        "stage3": stage,
        "fibrosis": fibrosis,
        "fibrosis_group3": np.asarray(
            ["F0" if value == 0 else "F1" if value == 1 else "F2_3" for value in fibrosis],
            dtype=str,
        ),
        "nas": np.asarray(
            [int(row["nash_crn_component_sum"]) for row in rows], dtype=np.float64
        ),
        "sex": sex,
    }


def _probabilities(
    rows: Sequence[Mapping[str, str]], prefix: str, roster: Sequence[str]
) -> np.ndarray:
    values = np.asarray(
        [[float(row[f"{prefix}{label}"]) for label in roster] for row in rows],
        dtype=np.float64,
    )
    if (
        not np.all(np.isfinite(values))
        or np.any(values < 0.0)
        or not np.allclose(values.sum(axis=1), 1.0, atol=1.0e-6, rtol=0.0)
    ):
        raise HistologyProductionScoringError("prediction probabilities differ")
    return values


def _load_predictions(rows: Sequence[Mapping[str, str]]) -> dict[str, np.ndarray]:
    stage_probabilities = _probabilities(rows, "probability_", reference.STAGE3)
    fibrosis_probabilities = _probabilities(
        rows, "probability_fibrosis_", reference.FIBROSIS_GROUP3
    )
    result = {
        "stage_probabilities": stage_probabilities,
        "stage3": np.asarray([row["predicted_stage3"] for row in rows], dtype=str),
        "fibrosis_probabilities": fibrosis_probabilities,
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
    return _validate_prediction_arrays(result)


def _validate_prediction_arrays(
    result: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    if (
        set(result)
        != {
            "stage_probabilities",
            "stage3",
            "fibrosis_probabilities",
            "fibrosis_group3",
            "fibrosis_cumulative",
            "fibrosis_regression",
            "nas",
        }
        or result["stage_probabilities"].shape
        != (EXPECTED_PARTICIPANTS, len(reference.STAGE3))
        or result["fibrosis_probabilities"].shape
        != (EXPECTED_PARTICIPANTS, len(reference.FIBROSIS_GROUP3))
        or not np.all(np.isfinite(result["stage_probabilities"]))
        or not np.all(np.isfinite(result["fibrosis_probabilities"]))
        or np.any(result["stage_probabilities"] < 0.0)
        or np.any(result["fibrosis_probabilities"] < 0.0)
        or not np.allclose(
            result["stage_probabilities"].sum(axis=1),
            1.0,
            atol=1.0e-6,
            rtol=0.0,
        )
        or not np.allclose(
            result["fibrosis_probabilities"].sum(axis=1),
            1.0,
            atol=1.0e-6,
            rtol=0.0,
        )
        or not set(result["stage3"]) <= set(reference.STAGE3)
        or not set(result["fibrosis_group3"]) <= set(reference.FIBROSIS_GROUP3)
        or not np.array_equal(
            result["stage3"],
            np.asarray(reference.STAGE3)[
                np.argmax(result["stage_probabilities"], axis=1)
            ],
        )
        or not np.array_equal(
            result["fibrosis_group3"],
            np.asarray(reference.FIBROSIS_GROUP3)[
                np.argmax(result["fibrosis_probabilities"], axis=1)
            ],
        )
        or any(
            not np.all(np.isfinite(result[key]))
            for key in ("fibrosis_cumulative", "fibrosis_regression", "nas")
        )
    ):
        raise HistologyProductionScoringError("prediction values differ")
    return result


def _load_frozen_endpoint_predictions(
    *,
    prediction_root: Path,
    prediction_records: Mapping[str, Mapping[str, Any]],
    production_sha256: str,
    participants: Sequence[str],
    outer_folds: Sequence[str],
) -> dict[str, dict[str, np.ndarray]]:
    """Load scoring inputs only from the 55 frozen endpoint bundles."""
    verify_frozen_tree(prediction_root)
    expected_sets = {
        f"{model_id}--endpoint-{endpoint_id}"
        for model_id, endpoint_id, _ in _endpoint_bundles()
    }
    if set(prediction_records) != expected_sets or len(prediction_records) != 55:
        raise HistologyProductionScoringError("frozen PredictionBundle roster differs")
    result: dict[str, dict[str, np.ndarray]] = {}
    for model_id in reference.MODEL_IDS:
        values: dict[str, np.ndarray] = {}
        for endpoint_id in PREDICTED_ENDPOINTS:
            set_id = f"{model_id}--endpoint-{endpoint_id}"
            target = prediction_root / set_id
            manifest = verify_frozen_tree(target)
            record = prediction_records[set_id]
            bundle = PredictionBundle.load_json(target / "prediction_bundle.json")
            bundle.validate_artifacts(target)
            fields, rows = _read_tsv(
                bundle.standardized_table.validate(target, require_relative=True)
            )
            expected_fields = STANDARDIZED_BASE_FIELDS + ENDPOINT_VALUE_FIELDS[
                endpoint_id
            ]
            if (
                fields != expected_fields
                or len(rows) != EXPECTED_PARTICIPANTS
                or [row["participant_id"] for row in rows] != list(participants)
                or [row["outer_fold"] for row in rows] != list(outer_folds)
                or any(row["endpoint_id"] != endpoint_id for row in rows)
                or any(row["prediction_state"] != "observed" for row in rows)
                or bundle.bundle_id != record.get("bundle_id")
                or sha256_file(target / "ARTIFACTS.json")
                != record.get("bundle_artifacts_sha256")
                or manifest.get("metadata", {}).get("outcomes_read") is not False
                or manifest.get("metadata", {}).get("metrics_calculated") is not False
            ):
                raise HistologyProductionScoringError(
                    "frozen PredictionBundle scoring input differs"
                )
            for row, participant, outer_fold in zip(
                rows, participants, outer_folds, strict=True
            ):
                if row["prediction_row_id"] != _prediction_row_id(
                    production_sha256=production_sha256,
                    model_id=model_id,
                    endpoint_id=endpoint_id,
                    participant_id=participant,
                    outer_fold=outer_fold,
                ):
                    raise HistologyProductionScoringError(
                        "frozen PredictionBundle row identity differs"
                    )
            if endpoint_id == "stage3":
                values["stage_probabilities"] = np.asarray(
                    [
                        [float(row[f"probability_{label}"]) for label in reference.STAGE3]
                        for row in rows
                    ],
                    dtype=np.float64,
                )
                values["stage3"] = np.asarray(
                    [row["predicted_stage3"] for row in rows], dtype=str
                )
            elif endpoint_id == "fibrosis_group3":
                values["fibrosis_probabilities"] = np.asarray(
                    [
                        [
                            float(row[f"probability_fibrosis_{label}"])
                            for label in reference.FIBROSIS_GROUP3
                        ]
                        for row in rows
                    ],
                    dtype=np.float64,
                )
                values["fibrosis_group3"] = np.asarray(
                    [row["predicted_fibrosis_group3"] for row in rows], dtype=str
                )
            elif endpoint_id == "fibrosis_cumulative":
                values["fibrosis_cumulative"] = np.asarray(
                    [
                        float(row["predicted_fibrosis_cumulative_expected"])
                        for row in rows
                    ],
                    dtype=np.float64,
                )
            elif endpoint_id == "fibrosis_regression":
                values["fibrosis_regression"] = np.asarray(
                    [float(row["predicted_fibrosis_regression"]) for row in rows],
                    dtype=np.float64,
                )
            else:
                values["nas"] = np.asarray(
                    [float(row["predicted_nash_crn_component_sum"]) for row in rows],
                    dtype=np.float64,
                )
        result[model_id] = _validate_prediction_arrays(values)
    if len(result) != 11:
        raise HistologyProductionScoringError("frozen PredictionBundle model roster differs")
    return result


def _macro_f1(observed: np.ndarray, predicted: np.ndarray, roster: Sequence[str]) -> float:
    values = []
    for label in roster:
        tp = int(np.sum((observed == label) & (predicted == label)))
        fp = int(np.sum((observed != label) & (predicted == label)))
        fn = int(np.sum((observed == label) & (predicted != label)))
        denominator = 2 * tp + fp + fn
        values.append(0.0 if denominator == 0 else 2.0 * tp / denominator)
    return float(np.mean(values))


def _point_metrics(
    outcomes: Mapping[str, np.ndarray], predictions: Mapping[str, np.ndarray]
) -> dict[str, float]:
    stage_one_hot = np.asarray(
        [[float(value == label) for label in reference.STAGE3] for value in outcomes["stage3"]]
    )
    values = {
        "stage3_macro_f1": _macro_f1(
            outcomes["stage3"], predictions["stage3"], reference.STAGE3
        ),
        "stage3_multiclass_brier": float(
            np.mean(np.sum((predictions["stage_probabilities"] - stage_one_hot) ** 2, axis=1))
        ),
        "fibrosis_group3_macro_f1": _macro_f1(
            outcomes["fibrosis_group3"],
            predictions["fibrosis_group3"],
            reference.FIBROSIS_GROUP3,
        ),
        "fibrosis_cumulative_ordinal_mae": float(
            np.mean(np.abs(outcomes["fibrosis"] - predictions["fibrosis_cumulative"]))
        ),
        "fibrosis_cumulative_spearman": float(
            spearmanr(outcomes["fibrosis"], predictions["fibrosis_cumulative"]).statistic
        ),
        "fibrosis_regression_mae": float(
            np.mean(np.abs(outcomes["fibrosis"] - predictions["fibrosis_regression"]))
        ),
        "fibrosis_regression_spearman": float(
            spearmanr(outcomes["fibrosis"], predictions["fibrosis_regression"]).statistic
        ),
        "nash_crn_component_sum_spearman": float(
            spearmanr(outcomes["nas"], predictions["nas"]).statistic
        ),
        "nash_crn_component_sum_mae": float(
            np.mean(np.abs(outcomes["nas"] - predictions["nas"]))
        ),
    }
    return values


def _bootstrap_macro_f1(
    observed: np.ndarray,
    predicted: np.ndarray,
    roster: Sequence[str],
    indices: np.ndarray,
) -> np.ndarray:
    observed_sample = observed[indices]
    predicted_sample = predicted[indices]
    result = np.zeros(len(indices), dtype=np.float64)
    for label in roster:
        tp = np.sum((observed_sample == label) & (predicted_sample == label), axis=1)
        fp = np.sum((observed_sample != label) & (predicted_sample == label), axis=1)
        fn = np.sum((observed_sample == label) & (predicted_sample != label), axis=1)
        denominator = 2 * tp + fp + fn
        result += np.divide(
            2.0 * tp,
            denominator,
            out=np.zeros_like(result),
            where=denominator != 0,
        )
    return result / len(roster)


def _bootstrap_spearman(
    observed: np.ndarray, predicted: np.ndarray, indices: np.ndarray
) -> np.ndarray:
    observed_rank = rankdata(observed[indices], axis=1, method="average")
    predicted_rank = rankdata(predicted[indices], axis=1, method="average")
    observed_rank -= observed_rank.mean(axis=1, keepdims=True)
    predicted_rank -= predicted_rank.mean(axis=1, keepdims=True)
    numerator = np.sum(observed_rank * predicted_rank, axis=1)
    denominator = np.sqrt(
        np.sum(observed_rank**2, axis=1) * np.sum(predicted_rank**2, axis=1)
    )
    return np.divide(
        numerator,
        denominator,
        out=np.full(len(indices), np.nan, dtype=np.float64),
        where=denominator > 0.0,
    )


def _bootstrap_metrics(
    outcomes: Mapping[str, np.ndarray],
    predictions: Mapping[str, np.ndarray],
    indices: np.ndarray,
) -> dict[str, np.ndarray]:
    stage_one_hot = np.asarray(
        [[float(value == label) for label in reference.STAGE3] for value in outcomes["stage3"]]
    )
    stage_brier = np.sum(
        (predictions["stage_probabilities"] - stage_one_hot) ** 2, axis=1
    )
    return {
        "stage3_macro_f1": _bootstrap_macro_f1(
            outcomes["stage3"], predictions["stage3"], reference.STAGE3, indices
        ),
        "stage3_multiclass_brier": np.mean(stage_brier[indices], axis=1),
        "fibrosis_group3_macro_f1": _bootstrap_macro_f1(
            outcomes["fibrosis_group3"],
            predictions["fibrosis_group3"],
            reference.FIBROSIS_GROUP3,
            indices,
        ),
        "fibrosis_cumulative_ordinal_mae": np.mean(
            np.abs(outcomes["fibrosis"] - predictions["fibrosis_cumulative"])[indices],
            axis=1,
        ),
        "fibrosis_cumulative_spearman": _bootstrap_spearman(
            outcomes["fibrosis"], predictions["fibrosis_cumulative"], indices
        ),
        "fibrosis_regression_mae": np.mean(
            np.abs(outcomes["fibrosis"] - predictions["fibrosis_regression"])[indices],
            axis=1,
        ),
        "fibrosis_regression_spearman": _bootstrap_spearman(
            outcomes["fibrosis"], predictions["fibrosis_regression"], indices
        ),
        "nash_crn_component_sum_spearman": _bootstrap_spearman(
            outcomes["nas"], predictions["nas"], indices
        ),
        "nash_crn_component_sum_mae": np.mean(
            np.abs(outcomes["nas"] - predictions["nas"])[indices], axis=1
        ),
    }


def _format_metric(value: float) -> str:
    return "not_estimable" if not math.isfinite(value) else format(value, ".17g")


def _metric_rows(
    *,
    model_id: str,
    endpoint_id: str,
    set_id: str,
    points: Mapping[str, float],
    distributions: Mapping[str, np.ndarray],
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for metric_id, endpoint_role, direction in METRIC_SPECS:
        applicable = metric_id in ENDPOINT_METRICS[endpoint_id]
        if applicable:
            distribution = distributions[metric_id]
            valid = distribution[np.isfinite(distribution)]
            low = float(np.quantile(valid, 0.025)) if len(valid) else math.nan
            high = float(np.quantile(valid, 0.975)) if len(valid) else math.nan
            estimate = _format_metric(points[metric_id])
            ci95_low = _format_metric(low)
            ci95_high = _format_metric(high)
            valid_replicates = len(valid)
            applicability_state = "observed"
            applicability_reason = "metric_native_to_endpoint_bundle"
        else:
            estimate = "not_applicable"
            ci95_low = "not_applicable"
            ci95_high = "not_applicable"
            valid_replicates = 0
            applicability_state = "not_applicable"
            applicability_reason = "metric_not_native_to_endpoint_bundle"
        family, scope, method, p_value, q_value = _metric_policy(endpoint_role)
        rows.append(
            {
                "task_id": TASK_ID,
                "model_id": model_id,
                "endpoint_id": endpoint_id,
                "prediction_aggregation": SEED_AGGREGATION,
                "prediction_set_id": set_id,
                "metric_id": metric_id,
                "endpoint_role": endpoint_role,
                "analysis_role": "ensemble_development_endpoint",
                "applicability_state": applicability_state,
                "applicability_reason": applicability_reason,
                "direction": direction,
                "estimate": estimate,
                "ci95_low": ci95_low,
                "ci95_high": ci95_high,
                "valid_bootstrap_replicates": valid_replicates,
                "biological_resampling_unit": "participant",
                "multiplicity_family": family,
                "multiplicity_scope": scope,
                "multiplicity_method": method,
                "p_value": p_value,
                "bh_adjusted_q_value": q_value,
                "confirmatory_inference_allowed": "false",
                "claim_mode": CLAIM_MODE,
            }
        )
    return rows


def _sex_rows(
    *,
    model_id: str,
    set_id: str,
    outcomes: Mapping[str, np.ndarray],
    predictions: Mapping[str, np.ndarray],
) -> list[dict[str, Any]]:
    result: list[dict[str, Any]] = []
    for recorded_sex in ("F", "M"):
        mask = outcomes["sex"] == recorded_sex
        nash_mask = mask & (outcomes["stage3"] == "NASH")
        stage_one_hot = np.asarray(
            [
                [float(value == label) for label in reference.STAGE3]
                for value in outcomes["stage3"][mask]
            ]
        )
        brier = float(
            np.mean(
                np.sum(
                    (predictions["stage_probabilities"][mask] - stage_one_hot) ** 2,
                    axis=1,
                )
            )
        )
        result.append(
            {
                "model_id": model_id,
                "endpoint_id": "stage3",
                "prediction_aggregation": SEED_AGGREGATION,
                "prediction_set_id": set_id,
                "recorded_sex": recorded_sex,
                "participants": int(mask.sum()),
                "stage3_census": json.dumps(
                    dict(sorted(Counter(outcomes["stage3"][mask]).items())),
                    sort_keys=True,
                    separators=(",", ":"),
                ),
                "stage3_brier": format(brier, ".17g"),
                "stage3_macro_f1": (
                    format(
                        _macro_f1(
                            outcomes["stage3"][mask],
                            predictions["stage3"][mask],
                            reference.STAGE3,
                        ),
                        ".17g",
                    )
                    if recorded_sex == "F"
                    else "not_applicable_missing_NOR_and_NAFL"
                ),
                "nash_recall": format(
                    float(np.mean(predictions["stage3"][nash_mask] == "NASH")),
                    ".17g",
                ),
                "mean_nash_true_class_probability": format(
                    float(
                        np.mean(
                            predictions["stage_probabilities"][
                                nash_mask, reference.STAGE3.index("NASH")
                            ]
                        )
                    ),
                    ".17g",
                ),
                "between_sex_global_gap_claim_allowed": "false",
            }
        )
    return result


def score(
    *,
    production: Path,
    production_artifacts_sha256: str,
    outcomes: Path,
    folds: Path,
    reference_evaluator: Path,
    task_contract: Path,
    surface: Path,
    scoring_contract: Path,
    scoring_contract_sha256: str,
    output: Path,
    bootstrap_replicates: int = BOOTSTRAP_REPLICATES,
    bootstrap_seed: int = BOOTSTRAP_SEED,
) -> dict[str, Any]:
    if output.exists():
        raise HistologyProductionScoringError(f"refusing to overwrite scoring output: {output}")
    if bootstrap_replicates != BOOTSTRAP_REPLICATES or bootstrap_seed != BOOTSTRAP_SEED:
        raise HistologyProductionScoringError("bootstrap contract differs")
    scorer_source_sha256 = sha256_file(Path(__file__).resolve())
    for path, expected, label in (
        (reference_evaluator, REFERENCE_EVALUATOR_SHA256, "reference evaluator"),
        (task_contract, TASK_CONTRACT_SHA256, "task contract"),
        (surface, SURFACE_SHA256, "benchmark surface"),
        (scoring_contract, scoring_contract_sha256, "scoring contract"),
        (folds / "ARTIFACTS.json", FOLDS_ARTIFACTS_SHA256, "folds"),
    ):
        _verify_file(path, expected, label)
    contract = _read_json(scoring_contract)
    _validate_contract(contract, production_artifacts_sha256=production_artifacts_sha256)
    source_validation_artifacts_sha256 = contract["independent_source_validation"][
        "artifacts_sha256"
    ]
    production_state = _verify_production(production, production_artifacts_sha256)
    verify_frozen_tree(folds)
    participants, outer_folds = _participant_axis(folds)
    output.mkdir(parents=True)
    prediction_records, prediction_collection_sha256 = _freeze_prediction_bundles(
        source_root=production_state["prediction_root"],
        output=output / "prediction_bundles",
        production_sha256=production_artifacts_sha256,
        participants=participants,
        outer_folds=outer_folds,
    )
    disposition_records, disposition_collection_sha256 = (
        _freeze_source_stage5_dispositions(
            output=output / "endpoint_dispositions",
            production_sha256=production_artifacts_sha256,
        )
    )

    _verify_file(outcomes / "ARTIFACTS.json", OUTCOMES_ARTIFACTS_SHA256, "outcomes")
    verify_frozen_tree(outcomes)
    outcome_values = _load_outcomes(outcomes, participants, outer_folds)
    scoring_predictions = _load_frozen_endpoint_predictions(
        prediction_root=output / "prediction_bundles",
        prediction_records=prediction_records,
        production_sha256=production_artifacts_sha256,
        participants=participants,
        outer_folds=outer_folds,
    )
    bootstrap_indices = np.random.default_rng(bootstrap_seed).integers(
        0,
        EXPECTED_PARTICIPANTS,
        size=(bootstrap_replicates, EXPECTED_PARTICIPANTS),
        endpoint=False,
    )
    bootstrap_indices_sha256 = canonical_sha256(bootstrap_indices.tolist())
    metrics_root = output / "metric_artifacts"
    metrics_root.mkdir()
    all_metric_rows: list[dict[str, Any]] = []
    all_sex_rows: list[dict[str, Any]] = []
    metric_index_rows: list[dict[str, Any]] = []
    model_results: dict[
        str, tuple[dict[str, float], dict[str, np.ndarray], dict[str, np.ndarray]]
    ] = {}
    for model_id in reference.MODEL_IDS:
        predictions = scoring_predictions[model_id]
        model_results[model_id] = (
            _point_metrics(outcome_values, predictions),
            _bootstrap_metrics(outcome_values, predictions, bootstrap_indices),
            predictions,
        )
        all_sex_rows.extend(
            _sex_rows(
                model_id=model_id,
                set_id=f"{model_id}--endpoint-stage3",
                outcomes=outcome_values,
                predictions=predictions,
            )
        )
    for model_id, endpoint_id, filename in _endpoint_bundles():
        set_id = f"{model_id}--endpoint-{endpoint_id}"
        points, distributions, _ = model_results[model_id]
        rows = _metric_rows(
            model_id=model_id,
            endpoint_id=endpoint_id,
            set_id=set_id,
            points=points,
            distributions=distributions,
        )
        all_metric_rows.extend(rows)
        target = metrics_root / set_id
        target.mkdir()
        metric_table = target / "metrics.tsv"
        _write_tsv(metric_table, METRIC_FIELDS, rows)
        metric_ref = ArtifactRef.from_path(
            metric_table,
            relative_to=target,
            media_type="text/tab-separated-values",
            role=f"development_metrics:{TASK_ID}:{set_id}",
        )
        prediction_binding = prediction_records[set_id]
        identity = {
            "schema_version": "masld-bench-development-metric-artifact-v1",
            "scoring_revision_id": SCORING_REVISION_ID,
            "task_id": TASK_ID,
            "model_id": model_id,
            "endpoint_id": endpoint_id,
            "prediction_aggregation": SEED_AGGREGATION,
            "model_seed_roster": list(reference.SEEDS),
            "outer_fold_aggregation": "each_participant_once_from_held_fold",
            "prediction_set_id": set_id,
            "production_artifacts_sha256": production_artifacts_sha256,
            "prediction_bundle_id": prediction_binding["bundle_id"],
            "prediction_bundle_artifacts_sha256": prediction_binding[
                "bundle_artifacts_sha256"
            ],
            "prediction_collection_artifacts_sha256": prediction_collection_sha256,
            "outcomes_artifacts_sha256": OUTCOMES_ARTIFACTS_SHA256,
            "folds_artifacts_sha256": FOLDS_ARTIFACTS_SHA256,
            "task_contract_sha256": TASK_CONTRACT_SHA256,
            "surface_sha256": SURFACE_SHA256,
            "reference_evaluator_sha256": REFERENCE_EVALUATOR_SHA256,
            "scoring_contract_sha256": scoring_contract_sha256,
            "scorer_source_sha256": scorer_source_sha256,
            "metrics_table": metric_ref.to_dict(),
            "metric_ids": [metric_id for metric_id, _, _ in METRIC_SPECS],
            "endpoint_native_metric_ids": sorted(ENDPOINT_METRICS[endpoint_id]),
            "endpoint_native_computed_metric_rows": len(
                ENDPOINT_METRICS[endpoint_id]
            ),
            "non_native_not_applicable_metric_rows": len(METRIC_SPECS)
            - len(ENDPOINT_METRICS[endpoint_id]),
            "bootstrap_replicates": bootstrap_replicates,
            "bootstrap_seed": bootstrap_seed,
            "bootstrap_indices_sha256": bootstrap_indices_sha256,
            "biological_resampling_unit": "participant",
            "secondary_multiplicity_family": MULTIPLICITY_FAMILY,
            "secondary_multiplicity_scope": MULTIPLICITY_SCOPE,
            "secondary_multiplicity_method": MULTIPLICITY_METHOD,
            "secondary_p_and_q_state": NO_CONFIRMATORY_STATE,
            "model_ranking_performed": False,
            "champion_claim_allowed": False,
            "external_claim_allowed": False,
            "diagnostic_or_prognostic_claim_allowed": False,
            "confirmatory_inference_allowed": False,
            "claim_mode": CLAIM_MODE,
        }
        payload = {"metric_artifact_id": canonical_sha256(identity), **identity}
        document = target / "development_metric_artifact.json"
        write_json_exclusive(document, payload, mode=0o440)
        metric_artifacts_sha256 = freeze_tree(
            target,
            {
                "artifact_class": "gse267145_histology_development_metric_artifact",
                "metric_artifact_id": payload["metric_artifact_id"],
                "model_id": model_id,
                "endpoint_id": endpoint_id,
                "prediction_aggregation": SEED_AGGREGATION,
                "metric_rows": len(METRIC_SPECS),
                "endpoint_native_computed_metric_rows": len(
                    ENDPOINT_METRICS[endpoint_id]
                ),
                "non_native_not_applicable_metric_rows": len(METRIC_SPECS)
                - len(ENDPOINT_METRICS[endpoint_id]),
                "biological_resampling_unit": "participant",
                "p_values_calculated": False,
                "bh_adjustment_calculated": False,
                "champion_claim_allowed": False,
                "external_claim_allowed": False,
                "diagnostic_or_prognostic_claim_allowed": False,
                "status": "passed_development_metrics",
            },
        )
        metric_index_rows.append(
            {
                "prediction_set_id": set_id,
                "model_id": model_id,
                "endpoint_id": endpoint_id,
                "prediction_aggregation": SEED_AGGREGATION,
                "metric_artifact_path": f"{set_id}/development_metric_artifact.json",
                "metric_artifact_id": payload["metric_artifact_id"],
                "metric_artifact_document_sha256": sha256_file(document),
                "metric_artifact_artifacts_sha256": metric_artifacts_sha256,
                "prediction_bundle_id": prediction_binding["bundle_id"],
                "prediction_bundle_artifacts_sha256": prediction_binding[
                    "bundle_artifacts_sha256"
                ],
            }
        )
    computed_rows = sum(
        row["applicability_state"] == "observed" for row in all_metric_rows
    )
    not_applicable_rows = sum(
        row["applicability_state"] == "not_applicable" for row in all_metric_rows
    )
    if (
        len(prediction_records) != 55
        or len(disposition_records) != 11
        or len(metric_index_rows) != 55
        or len(all_metric_rows) != 495
        or computed_rows != 99
        or not_applicable_rows != 396
        or len(all_sex_rows) != 22
    ):
        raise HistologyProductionScoringError("scoring artifact count mapping differs")
    _write_tsv(output / "standardized_metrics.tsv", METRIC_FIELDS, all_metric_rows)
    _write_tsv(output / "sex_error_audit.tsv", SEX_AUDIT_FIELDS, all_sex_rows)
    _write_tsv(
        metrics_root / "metric_artifact_index.tsv",
        (
            "prediction_set_id",
            "model_id",
            "endpoint_id",
            "prediction_aggregation",
            "metric_artifact_path",
            "metric_artifact_id",
            "metric_artifact_document_sha256",
            "metric_artifact_artifacts_sha256",
            "prediction_bundle_id",
            "prediction_bundle_artifacts_sha256",
        ),
        metric_index_rows,
    )
    metric_collection_sha256 = freeze_tree(
        metrics_root,
        {
            "artifact_class": "gse267145_histology_development_metric_collection",
            "metric_artifacts": len(metric_index_rows),
            "metric_rows": len(all_metric_rows),
            "endpoint_native_computed_metric_rows": computed_rows,
            "non_native_not_applicable_metric_rows": not_applicable_rows,
            "biological_resampling_unit": "participant",
            "p_values_calculated": False,
            "bh_adjustment_calculated": False,
            "model_ranking_performed": False,
            "champion_claim_allowed": False,
            "external_claim_allowed": False,
            "diagnostic_or_prognostic_claim_allowed": False,
            "status": "passed_development_metrics",
        },
    )
    receipt = {
        "schema_version": "masld-bench-gse267145-histology-production-v4-scoring-v1",
        "status": "passed_development_metrics_no_confirmatory_claim",
        "scoring_revision_id": SCORING_REVISION_ID,
        "production_revision_id": PRODUCTION_REVISION_ID,
        "production_artifacts_sha256": production_artifacts_sha256,
        "scoring_contract_sha256": scoring_contract_sha256,
        "scorer_source_sha256": scorer_source_sha256,
        "reference_evaluator_sha256": REFERENCE_EVALUATOR_SHA256,
        "task_contract_sha256": TASK_CONTRACT_SHA256,
        "surface_sha256": SURFACE_SHA256,
        "outcomes_artifacts_sha256": OUTCOMES_ARTIFACTS_SHA256,
        "folds_artifacts_sha256": FOLDS_ARTIFACTS_SHA256,
        "aggregate_artifacts_sha256": production_state["aggregate_artifacts_sha256"],
        "source_validation_artifacts_sha256": source_validation_artifacts_sha256,
        "prediction_collection_artifacts_sha256": prediction_collection_sha256,
        "endpoint_disposition_collection_artifacts_sha256": disposition_collection_sha256,
        "metric_collection_artifacts_sha256": metric_collection_sha256,
        "standardized_prediction_bundles": len(prediction_records),
        "source_stage5_endpoint_dispositions": len(disposition_records),
        "source_stage5_disposition": SOURCE_STAGE5_DISPOSITION,
        "standardized_metric_artifacts": len(metric_index_rows),
        "standardized_metric_rows": len(all_metric_rows),
        "endpoint_native_computed_metric_rows": computed_rows,
        "non_native_not_applicable_metric_rows": not_applicable_rows,
        "models": len(reference.MODEL_IDS),
        "predicted_endpoints_per_model": list(PREDICTED_ENDPOINTS),
        "real_prediction_bundles_per_model": len(PREDICTED_ENDPOINTS),
        "prediction_aggregation": SEED_AGGREGATION,
        "model_seed_roster": list(reference.SEEDS),
        "outer_fold_aggregation": "each_participant_once_from_held_fold",
        "seed_or_fold_results_emitted": False,
        "inference_mapping": (
            "participant bootstrap over 99 participants after held-fold assembly "
            "and fixed five-seed prediction averaging"
        ),
        "participants_per_prediction_set": EXPECTED_PARTICIPANTS,
        "primary_metrics": 1,
        "calibration_metrics": 1,
        "secondary_metrics": 7,
        "source_stage5_state": SOURCE_STAGE5_DISPOSITION,
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_seed": bootstrap_seed,
        "bootstrap_indices_sha256": bootstrap_indices_sha256,
        "biological_resampling_unit": "participant",
        "seeds_are_biological_replicates": False,
        "prediction_bundles_frozen_before_outcomes_read": True,
        "five_seed_ensemble_rederived_before_outcomes_read": True,
        "raw_production_prediction_files_read_after_outcomes": False,
        "metric_scoring_prediction_source": "55_frozen_endpoint_prediction_bundles",
        "prediction_bundles_mutated": False,
        "model_or_fit_source_imported": False,
        "baseline_comparison_performed": False,
        "model_ranking_performed": False,
        "p_values_calculated": False,
        "bh_adjustment_calculated": False,
        "secondary_multiplicity_family_retained": MULTIPLICITY_FAMILY,
        "secondary_multiplicity_scope_retained": MULTIPLICITY_SCOPE,
        "secondary_multiplicity_method_retained": MULTIPLICITY_METHOD,
        "secondary_p_and_q_state": NO_CONFIRMATORY_STATE,
        "champion_claim_allowed": False,
        "external_transfer_claim_allowed": False,
        "diagnostic_or_prognostic_claim_allowed": False,
        "confirmatory_inference_allowed": False,
        "claim_mode": CLAIM_MODE,
    }
    write_json_exclusive(output / "scoring_receipt.json", receipt, mode=0o440)
    freeze_tree(
        output,
        {
            "artifact_class": "gse267145_histology_production_v4_development_scores",
            "production_artifacts_sha256": production_artifacts_sha256,
            "production_aggregate_artifacts_sha256": (
                PRODUCTION_AGGREGATE_ARTIFACTS_SHA256
            ),
            "source_validation_artifacts_sha256": (
                source_validation_artifacts_sha256
            ),
            "scoring_contract_sha256": scoring_contract_sha256,
            "scorer_source_sha256": scorer_source_sha256,
            "prediction_bundles": len(prediction_records),
            "source_stage5_endpoint_dispositions": len(disposition_records),
            "metric_artifacts": len(metric_index_rows),
            "metric_rows": len(all_metric_rows),
            "endpoint_native_computed_metric_rows": computed_rows,
            "non_native_not_applicable_metric_rows": not_applicable_rows,
            "biological_resampling_unit": "participant",
            "p_values_calculated": False,
            "bh_adjustment_calculated": False,
            "model_ranking_performed": False,
            "champion_claim_allowed": False,
            "external_claim_allowed": False,
            "diagnostic_or_prognostic_claim_allowed": False,
            "status": "passed_development_metrics",
        },
    )
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--production", required=True, type=Path)
    parser.add_argument("--production-artifacts-sha256", required=True)
    parser.add_argument("--outcomes", required=True, type=Path)
    parser.add_argument("--folds", required=True, type=Path)
    parser.add_argument("--reference-evaluator", required=True, type=Path)
    parser.add_argument("--task-contract", required=True, type=Path)
    parser.add_argument("--surface", required=True, type=Path)
    parser.add_argument("--scoring-contract", required=True, type=Path)
    parser.add_argument("--scoring-contract-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    receipt = score(
        production=arguments.production,
        production_artifacts_sha256=arguments.production_artifacts_sha256,
        outcomes=arguments.outcomes,
        folds=arguments.folds,
        reference_evaluator=arguments.reference_evaluator,
        task_contract=arguments.task_contract,
        surface=arguments.surface,
        scoring_contract=arguments.scoring_contract,
        scoring_contract_sha256=arguments.scoring_contract_sha256,
        output=arguments.output,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
