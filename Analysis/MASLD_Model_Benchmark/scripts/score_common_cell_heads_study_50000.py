#!/usr/bin/env python3
"""Independently score audited 50,000-cell common-head predictions."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import re
from typing import Any

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import canonical_sha256, sha256_file
from scripts import score_cell_baselines_study_50000 as metric_contract


DATASET_VIEW_ID = metric_contract.DATASET_VIEW_ID
SPLIT_ID = metric_contract.SPLIT_ID
ROSTER = metric_contract.ROSTER
HEADS = ("linear", "two_layer_mlp")
SCREEN_SEEDS = (1103, 1201, 1301)
EXPECTED_ROWS = metric_contract.EXPECTED_ROWS
EXPECTED_DONORS = metric_contract.EXPECTED_DONORS
EXPECTED_STUDIES = metric_contract.EXPECTED_STUDIES
BOOTSTRAP_RESAMPLES = metric_contract.BOOTSTRAP_RESAMPLES
BOOTSTRAP_SEED = metric_contract.BOOTSTRAP_SEED
MODEL_ID_PATTERN = re.compile(r"^[a-z0-9][a-z0-9._-]{2,79}$")
EXPOSURE_STATUSES = frozenset(
    {
        "clean_declared",
        "target_label_unexposed",
        "encoder_seen",
        "continual_seen",
        "reference_only",
        "downstream_demo",
        "unknown",
    }
)


class CommonHeadScoreError(ValueError):
    """Raised when development predictions cannot be independently evaluated."""


def require_model_identity(
    model_id: str, exposure_status: str, sealed_champion_eligible: bool
) -> tuple[str, str, bool]:
    if MODEL_ID_PATTERN.fullmatch(model_id) is None:
        raise CommonHeadScoreError("expected model_id is invalid")
    if exposure_status not in EXPOSURE_STATUSES:
        raise CommonHeadScoreError("expected exposure status is invalid")
    if type(sealed_champion_eligible) is not bool:
        raise CommonHeadScoreError("expected champion eligibility is invalid")
    if sealed_champion_eligible and exposure_status not in {
        "clean_declared",
        "target_label_unexposed",
    }:
        raise CommonHeadScoreError("champion eligibility contradicts exposure")
    return model_id, exposure_status, sealed_champion_eligible


def mean_seed_probabilities(values: list[np.ndarray]) -> np.ndarray:
    if len(values) != len(SCREEN_SEEDS) or any(
        value.shape != (EXPECTED_ROWS, len(ROSTER)) for value in values
    ):
        raise CommonHeadScoreError("seed prediction tensors differ")
    result = np.mean(np.stack(values, axis=0), axis=0)
    if not np.all(np.isfinite(result)) or not np.allclose(
        result.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-12
    ):
        raise CommonHeadScoreError("seed-mean probability ensemble differs")
    return result


def run(
    *,
    source: Path,
    split: Path,
    predictions: Path,
    output: Path,
    expected_source_artifacts_sha256: str,
    expected_split_artifacts_sha256: str,
    expected_prediction_artifacts_sha256: str,
    expected_model_id: str,
    expected_exposure_status: str,
    expected_sealed_champion_eligible: bool,
) -> None:
    if output.exists():
        raise CommonHeadScoreError("refusing to overwrite common-head evaluation")
    (
        expected_model_id,
        expected_exposure_status,
        expected_sealed_champion_eligible,
    ) = require_model_identity(
        expected_model_id,
        expected_exposure_status,
        expected_sealed_champion_eligible,
    )
    for root, expected in (
        (source, expected_source_artifacts_sha256),
        (split, expected_split_artifacts_sha256),
        (predictions, expected_prediction_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise CommonHeadScoreError("input ARTIFACTS SHA-256 differs")
        verify_frozen_tree(root)

    source_metadata = verify_frozen_tree(source)["metadata"]
    split_metadata = verify_frozen_tree(split)["metadata"]
    prediction_metadata = verify_frozen_tree(predictions)["metadata"]
    if (
        source_metadata.get("subset_id") != DATASET_VIEW_ID
        or source_metadata.get("row_count") != EXPECTED_ROWS
        or source_metadata.get("sealed_outcomes_read") is not False
        or source_metadata.get("histology_read") is not False
        or split_metadata.get("split_id") != SPLIT_ID
        or split_metadata.get("rows") != EXPECTED_ROWS
        or split_metadata.get("donors") != EXPECTED_DONORS
        or split_metadata.get("studies") != EXPECTED_STUDIES
        or split_metadata.get("target_labels_used_for_assignment") is not False
        or split_metadata.get("sealed_outcomes_read") is not False
        or prediction_metadata.get("artifact_class")
        != "cell_foundation_common_head_audited_predictions"
        or prediction_metadata.get("model_id") != expected_model_id
        or prediction_metadata.get("dataset_view_id") != DATASET_VIEW_ID
        or prediction_metadata.get("split_id") != SPLIT_ID
        or prediction_metadata.get("rows") != EXPECTED_ROWS
        or prediction_metadata.get("metrics_calculated") is not False
        or prediction_metadata.get("sealed_outcomes_read") is not False
    ):
        raise CommonHeadScoreError("frozen evaluator input metadata differs")

    _, source_rows = metric_contract.read_tsv(source / "selection.tsv")
    split_fields, split_rows = metric_contract.read_tsv(split / "row_outer_folds.tsv")
    if split_fields != ["row_id", "donor_id", "dataset", "outer_fold"]:
        raise CommonHeadScoreError("split table schema differs")
    if len(source_rows) != EXPECTED_ROWS or len(split_rows) != EXPECTED_ROWS:
        raise CommonHeadScoreError("source or split row count differs")
    row_ids = np.asarray([row["row_id"] for row in source_rows], dtype=str)
    donors = np.asarray([row["donor_id"] for row in source_rows], dtype=str)
    studies = np.asarray([row["dataset"] for row in source_rows], dtype=str)
    labels = np.asarray([row["broad_label"] for row in source_rows], dtype=str)
    folds = np.asarray([int(row["outer_fold"]) for row in split_rows], dtype=np.int8)
    if (
        len(set(row_ids.tolist())) != EXPECTED_ROWS
        or len(set(donors.tolist())) != EXPECTED_DONORS
        or len(set(studies.tolist())) != EXPECTED_STUDIES
        or set(labels.tolist()) != set(ROSTER)
        or set(folds.tolist()) != set(range(5))
        or row_ids.tolist() != [row["row_id"] for row in split_rows]
        or donors.tolist() != [row["donor_id"] for row in split_rows]
        or studies.tolist() != [row["dataset"] for row in split_rows]
        or any(len(set(folds[donors == donor])) != 1 for donor in set(donors))
        or any(len(set(folds[studies == study])) != 1 for study in set(studies))
    ):
        raise CommonHeadScoreError("source labels or donor-safe split alignment differs")
    label_index = {label: index for index, label in enumerate(ROSTER)}
    truth = np.asarray([label_index[label] for label in labels], dtype=np.int8)

    receipt = json.loads(
        (predictions / "prediction_receipt.json").read_text(encoding="utf-8")
    )
    if (
        receipt.get("status") != "pass_development_predictions_audited"
        or receipt.get("model_id") != expected_model_id
        or receipt.get("heads") != list(HEADS)
        or receipt.get("screen_seeds") != list(SCREEN_SEEDS)
        or receipt.get("rows") != EXPECTED_ROWS
        or receipt.get("donors") != EXPECTED_DONORS
        or receipt.get("studies") != EXPECTED_STUDIES
        or receipt.get("head_seed_fold_census_complete") is not True
        or receipt.get("donor_safe_outer_folds_verified") is not True
        or receipt.get("whole_study_outer_folds_verified") is not True
        or receipt.get("probability_schema_verified") is not True
        or receipt.get("exposure_status") != expected_exposure_status
        or receipt.get("sealed_champion_eligible")
        is not expected_sealed_champion_eligible
        or receipt.get("metrics_calculated") is not False
        or receipt.get("development_labels_read_by_aggregator") != []
        or receipt.get("prediction_tables_contain_observed_labels") is not False
        or receipt.get("sealed_outcomes_read") is not False
        or receipt.get("source_artifacts_sha256_from_shards")
        != expected_source_artifacts_sha256
        or receipt.get("split_artifacts_sha256")
        != expected_split_artifacts_sha256
    ):
        raise CommonHeadScoreError("audited prediction receipt differs")

    prediction_fields = [
        "row_id",
        "donor_id",
        "dataset",
        "outer_fold",
        "predicted_class",
        *(f"probability::{label}" for label in ROSTER),
    ]

    def read_aligned_prediction(path: Path) -> np.ndarray:
        fields, rows = metric_contract.read_tsv(path)
        if fields != prediction_fields or len(rows) != EXPECTED_ROWS:
            raise CommonHeadScoreError("standardized prediction table schema differs")
        if any("label" in field.lower() or "histolog" in field.lower() for field in fields):
            raise CommonHeadScoreError("prediction table contains evaluator-only fields")
        if (
            [row["row_id"] for row in rows] != row_ids.tolist()
            or [row["donor_id"] for row in rows] != donors.tolist()
            or [row["dataset"] for row in rows] != studies.tolist()
            or [int(row["outer_fold"]) for row in rows] != folds.tolist()
        ):
            raise CommonHeadScoreError("prediction rows differ from source and split")
        values = metric_contract.probability_matrix(rows)
        expected_hard = [ROSTER[index] for index in np.argmax(values, axis=1)]
        if [row["predicted_class"] for row in rows] != expected_hard:
            raise CommonHeadScoreError("hard predictions differ from probabilities")
        return values

    head_probabilities: dict[str, np.ndarray] = {}
    seed_metrics_rows: list[dict[str, Any]] = []
    unique_studies = sorted(set(studies.tolist()))
    for head in HEADS:
        seed_values: list[np.ndarray] = []
        for seed in SCREEN_SEEDS:
            values = read_aligned_prediction(
                predictions / "predictions" / f"{head}__seed{seed}.tsv"
            )
            seed_values.append(values)
            overall = metric_contract.score_subset(truth, donors, values)
            study_macro = float(
                np.mean(
                    [
                        metric_contract.score_subset(
                            truth[studies == study],
                            donors[studies == study],
                            values[studies == study],
                        )["macro_f1"]
                        for study in unique_studies
                    ]
                )
            )
            seed_metrics_rows.append(
                {
                    "head_id": head,
                    "seed": seed,
                    "donor_class_balanced_macro_f1": overall["macro_f1"],
                    "study_balanced_macro_f1": study_macro,
                    "multiclass_brier": overall["multiclass_brier"],
                    "top_label_ece": overall["calibration"]["top_label_ece"],
                }
            )
        head_probabilities[head] = mean_seed_probabilities(seed_values)

    unique_donors = sorted(set(donors.tolist()))
    donor_to_index = {donor: index for index, donor in enumerate(unique_donors)}
    donor_indices = np.asarray([donor_to_index[donor] for donor in donors], dtype=np.int16)
    donor_studies = np.asarray(
        [next(iter(set(studies[donors == donor].tolist()))) for donor in unique_donors],
        dtype=str,
    )
    global_present = np.zeros((EXPECTED_DONORS, len(ROSTER)), dtype=bool)
    for donor in range(EXPECTED_DONORS):
        global_present[donor, np.unique(truth[donor_indices == donor])] = True
    multiplicities = metric_contract.build_multiplicities(
        donor_studies,
        global_present,
        n_resamples=BOOTSTRAP_RESAMPLES,
        seed=BOOTSTRAP_SEED,
    )
    if (
        multiplicities.shape != (BOOTSTRAP_RESAMPLES, EXPECTED_DONORS)
        or np.any(multiplicities.sum(axis=1) != EXPECTED_DONORS)
    ):
        raise CommonHeadScoreError("donor bootstrap multiplicities differ")

    output.mkdir(mode=0o750)
    metrics_rows: list[dict[str, Any]] = []
    class_rows: list[dict[str, Any]] = []
    donor_rows: list[dict[str, Any]] = []
    study_rows: list[dict[str, Any]] = []
    calibration_rows: list[dict[str, Any]] = []
    interval_rows: list[dict[str, Any]] = []
    distributions: dict[str, np.ndarray] = {}
    all_donor_indices = np.arange(EXPECTED_DONORS)
    all_classes = np.ones(len(ROSTER), dtype=bool)

    for head, probabilities in head_probabilities.items():
        overall = metric_contract.score_subset(truth, donors, probabilities)
        study_scores: dict[str, dict[str, Any]] = {}
        for study in unique_studies:
            selected = studies == study
            value = metric_contract.score_subset(
                truth[selected], donors[selected], probabilities[selected]
            )
            study_scores[study] = value
            study_rows.append(
                {
                    "head_id": head,
                    "study": study,
                    "rows": value["rows"],
                    "donors": value["donors"],
                    "evaluable_classes": ",".join(value["classes"]),
                    "donor_class_balanced_macro_f1": value["macro_f1"],
                    "multiclass_brier": value["multiclass_brier"],
                    "top_label_ece": value["calibration"]["top_label_ece"],
                    "mean_classwise_ece": value["calibration"]["mean_classwise_ece"],
                }
            )
        study_macro = float(np.mean([value["macro_f1"] for value in study_scores.values()]))
        point = {
            "donor_class_balanced_macro_f1": overall["macro_f1"],
            "study_balanced_macro_f1": study_macro,
            "multiclass_brier": overall["multiclass_brier"],
            "top_label_ece": overall["calibration"]["top_label_ece"],
            "mean_classwise_ece": overall["calibration"]["mean_classwise_ece"],
        }
        metrics_rows.append({"head_id": head, **point})
        for label in ROSTER:
            study_values = [
                value["per_class_f1"][label]
                for value in study_scores.values()
                if label in value["per_class_f1"]
            ]
            class_rows.append(
                {
                    "head_id": head,
                    "class": label,
                    "donor_class_balanced_f1": overall["per_class_f1"][label],
                    "study_balanced_f1": float(np.mean(study_values)),
                    "evaluable_studies": len(study_values),
                }
            )
        for record in overall["calibration_bins"]:
            calibration_rows.append({"head_id": head, "scope": "all_studies", **record})
        for donor in unique_donors:
            selected = donors == donor
            value = metric_contract.score_subset(
                truth[selected], donors[selected], probabilities[selected]
            )
            donor_rows.append(
                {
                    "head_id": head,
                    "donor_id": donor,
                    "study": studies[selected][0],
                    "rows": value["rows"],
                    "evaluable_classes": ",".join(value["classes"]),
                    "class_balanced_macro_f1": value["macro_f1"],
                    "multiclass_brier": value["multiclass_brier"],
                    "top_label_ece": value["calibration"]["top_label_ece"],
                }
            )

        stats = metric_contract.sufficient_statistics(
            truth,
            np.argmax(probabilities, axis=1),
            probabilities,
            donor_indices,
            donors=EXPECTED_DONORS,
        )
        identity = metric_contract.endpoints_from_multiplicities(
            stats,
            np.ones((1, EXPECTED_DONORS), dtype=np.int16),
            all_donor_indices,
            all_classes,
        )
        identity_study_parts = []
        for study in unique_studies:
            selected_donors = np.flatnonzero(donor_studies == study)
            study_classes = stats["present"][selected_donors].any(axis=0)
            identity_study_parts.append(
                metric_contract.endpoints_from_multiplicities(
                    stats,
                    np.ones((1, EXPECTED_DONORS), dtype=np.int16),
                    selected_donors,
                    study_classes,
                )["macro_f1"][0]
            )
        if (
            not math.isclose(
                float(identity["macro_f1"][0]),
                overall["macro_f1"],
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            or not math.isclose(
                float(identity["brier"][0]),
                overall["multiclass_brier"],
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            or not math.isclose(
                float(identity["top_label_ece"][0]),
                overall["calibration"]["top_label_ece"],
                rel_tol=0.0,
                abs_tol=1e-12,
            )
            or not math.isclose(
                float(np.mean(identity_study_parts)),
                study_macro,
                rel_tol=0.0,
                abs_tol=1e-12,
            )
        ):
            raise CommonHeadScoreError("independent sufficient-statistic rederivation differs")
        donor_bootstrap = metric_contract.endpoints_from_multiplicities(
            stats, multiplicities, all_donor_indices, all_classes
        )
        study_bootstrap_parts = []
        for study in unique_studies:
            selected_donors = np.flatnonzero(donor_studies == study)
            study_classes = stats["present"][selected_donors].any(axis=0)
            study_bootstrap_parts.append(
                metric_contract.endpoints_from_multiplicities(
                    stats, multiplicities, selected_donors, study_classes
                )["macro_f1"]
            )
        metric_distributions = {
            "donor_class_balanced_macro_f1": donor_bootstrap["macro_f1"],
            "study_balanced_macro_f1": np.mean(
                np.stack(study_bootstrap_parts, axis=1), axis=1
            ),
            "multiclass_brier": donor_bootstrap["brier"],
            "top_label_ece": donor_bootstrap["top_label_ece"],
        }
        for index, label in enumerate(ROSTER):
            metric_distributions[f"class_f1::{label}"] = donor_bootstrap["class_f1"][:, index]
        for metric, values in metric_distributions.items():
            distributions[f"{head}::{metric}"] = values.astype(np.float64, copy=False)
            estimate = (
                overall["per_class_f1"][metric.split("::", 1)[1]]
                if metric.startswith("class_f1::")
                else point[metric]
            )
            interval_rows.append(
                {"head_id": head, "metric": metric, "estimate": estimate, **metric_contract.interval(values)}
            )

    metric_contract.write_tsv(output / "head_metrics.tsv", list(metrics_rows[0]), metrics_rows)
    metric_contract.write_tsv(output / "seed_metrics.tsv", list(seed_metrics_rows[0]), seed_metrics_rows)
    metric_contract.write_tsv(output / "per_class_f1.tsv", list(class_rows[0]), class_rows)
    metric_contract.write_tsv(output / "donor_summaries.tsv", list(donor_rows[0]), donor_rows)
    metric_contract.write_tsv(output / "study_summaries.tsv", list(study_rows[0]), study_rows)
    metric_contract.write_tsv(output / "calibration_bins.tsv", list(calibration_rows[0]), calibration_rows)
    metric_contract.write_tsv(output / "bootstrap_intervals.tsv", list(interval_rows[0]), interval_rows)
    with (output / "bootstrap_distributions.npz").open("xb") as handle:
        np.savez_compressed(handle, **distributions)

    evaluation_receipt = {
        "schema_version": "masld-bench-common-cell-head-study-evaluator-v1",
        "status": "pass_development_evaluation_only",
        "evaluator_id": "common_cell_heads_study_50000_evaluation_only_v1",
        "model_id": expected_model_id,
        "checkpoint_sha256": receipt["checkpoint_sha256"],
        "exposure_status": receipt["exposure_status"],
        "sealed_champion_eligible": receipt["sealed_champion_eligible"],
        "dataset_view_id": DATASET_VIEW_ID,
        "split_id": SPLIT_ID,
        "rows": EXPECTED_ROWS,
        "donors": EXPECTED_DONORS,
        "studies": EXPECTED_STUDIES,
        "classes": list(ROSTER),
        "heads": list(HEADS),
        "screen_seeds": list(SCREEN_SEEDS),
        "final_prediction": "unweighted_mean_probability_ensemble_within_head",
        "metric_contract": "score_cell_baselines_study_50000_v1",
        "metrics_computed_by_independent_evaluator": True,
        "training_adapter_imported": False,
        "source_artifacts_sha256": expected_source_artifacts_sha256,
        "split_artifacts_sha256": expected_split_artifacts_sha256,
        "prediction_artifacts_sha256": expected_prediction_artifacts_sha256,
        "source_labels_read_inside_evaluator": ["broad_label"],
        "source_cell_type_used_for_scoring": False,
        "prediction_tables_contain_observed_labels": False,
        "models_fit": False,
        "preprocessing_fit": False,
        "calibration_fit": False,
        "head_selection_performed": False,
        "comparison_to_baseline_performed": False,
        "champion_promotion_performed": False,
        "project_sealed_external_evaluation_performed": False,
        "bootstrap": {
            "method": "study_stratified_class_coverage_preserving_donor_cluster_bootstrap",
            "resamples": BOOTSTRAP_RESAMPLES,
            "seed": BOOTSTRAP_SEED,
            "donors_are_independent_units": True,
            "cells_are_independent_units": False,
            "study_donor_counts_preserved": True,
            "distribution_sha256": canonical_sha256(
                {key: values.tolist() for key, values in sorted(distributions.items())}
            ),
        },
        "histology_read": False,
        "sealed_outcomes_read": False,
        "clinical_claim_allowed": False,
    }
    write_json_exclusive(output / "evaluation_receipt.json", evaluation_receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "cell_foundation_common_head_study_evaluation",
            "model_id": expected_model_id,
            "dataset_view_id": DATASET_VIEW_ID,
            "split_id": SPLIT_ID,
            "rows": EXPECTED_ROWS,
            "donors": EXPECTED_DONORS,
            "models_fit": False,
            "head_selection_performed": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )


def parser() -> argparse.ArgumentParser:
    value = argparse.ArgumentParser(description=__doc__)
    value.add_argument("--source", required=True, type=Path)
    value.add_argument("--split", required=True, type=Path)
    value.add_argument("--predictions", required=True, type=Path)
    value.add_argument("--output", required=True, type=Path)
    value.add_argument("--expected-source-artifacts-sha256", required=True)
    value.add_argument("--expected-split-artifacts-sha256", required=True)
    value.add_argument("--expected-prediction-artifacts-sha256", required=True)
    value.add_argument("--expected-model-id", required=True)
    value.add_argument(
        "--expected-exposure-status", required=True, choices=sorted(EXPOSURE_STATUSES)
    )
    value.add_argument(
        "--expected-sealed-champion-eligible",
        required=True,
        choices=("true", "false"),
    )
    return value


def main() -> int:
    arguments = parser().parse_args()
    run(
        source=arguments.source,
        split=arguments.split,
        predictions=arguments.predictions,
        output=arguments.output,
        expected_source_artifacts_sha256=arguments.expected_source_artifacts_sha256,
        expected_split_artifacts_sha256=arguments.expected_split_artifacts_sha256,
        expected_prediction_artifacts_sha256=arguments.expected_prediction_artifacts_sha256,
        expected_model_id=arguments.expected_model_id,
        expected_exposure_status=arguments.expected_exposure_status,
        expected_sealed_champion_eligible=(
            arguments.expected_sealed_champion_eligible == "true"
        ),
    )
    print(json.dumps({"output": arguments.output.resolve().as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
