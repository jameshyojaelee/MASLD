#!/usr/bin/env python3
"""Score frozen study-held-out cell predictions without fitting any model."""

from __future__ import annotations

import argparse
from collections import defaultdict
import csv
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.evaluators.metrics import (
    donor_class_balanced_weights,
    multiclass_brier_score,
    weighted_macro_f1,
)
from masld_bench.hashing import canonical_sha256, sha256_file


DATASET_VIEW_ID = "resource_atlas_frozen_screen_50000_v1"
SPLIT_ID = "resource_atlas_study_outer_5fold_v1"
ROSTER = (
    "cholangiocyte",
    "endothelial",
    "hepatocyte",
    "immune",
    "mesenchymal_stromal",
)
MODEL_IDS = (
    "hvg_pca_nearest_centroid",
    "hvg_pca_knn",
    "hvg_pca_logistic",
    "hvg_pca_elastic_net",
    "hvg_pca_linear_svm",
)
SCREEN_SEEDS = (20260824, 20260825, 20260826)
EXPECTED_ROWS = 50_000
EXPECTED_DONORS = 102
EXPECTED_STUDIES = 7
CALIBRATION_BINS = 10
BOOTSTRAP_RESAMPLES = 10_000
BOOTSTRAP_SEED = 20260824


class CellStudyScoreError(ValueError):
    """Raised when frozen evaluator inputs or metrics does not meet the requirements."""


def read_tsv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise CellStudyScoreError(f"TSV lacks a header: {path}")
        return list(reader.fieldnames), list(reader)


def write_tsv(
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


def probability_matrix(rows: Sequence[Mapping[str, str]]) -> np.ndarray:
    fields = [f"probability::{label}" for label in ROSTER]
    try:
        values = np.asarray(
            [[float(row[field]) for field in fields] for row in rows],
            dtype=np.float64,
        )
    except (KeyError, TypeError, ValueError) as error:
        raise CellStudyScoreError("prediction probability fields differ") from error
    if (
        values.shape != (len(rows), len(ROSTER))
        or not np.all(np.isfinite(values))
        or np.any(values < 0.0)
        or not np.allclose(values.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-7)
    ):
        raise CellStudyScoreError("prediction probabilities differ")
    return values


def calibration_metrics(
    truth: np.ndarray,
    probabilities: np.ndarray,
    weights: np.ndarray,
    *,
    bins: int = CALIBRATION_BINS,
) -> tuple[dict[str, Any], list[dict[str, Any]]]:
    if (
        truth.ndim != 1
        or probabilities.shape != (len(truth), len(ROSTER))
        or weights.shape != (len(truth),)
        or bins < 2
    ):
        raise CellStudyScoreError("calibration arrays differ")
    total_weight = float(weights.sum())
    if total_weight <= 0.0 or not np.all(np.isfinite(weights)) or np.any(weights < 0.0):
        raise CellStudyScoreError("calibration weights differ")
    normalized = weights / total_weight
    predicted = np.argmax(probabilities, axis=1)
    confidence = probabilities[np.arange(len(truth)), predicted]
    correct = (predicted == truth).astype(np.float64)
    bin_ids = np.minimum((confidence * bins).astype(np.int64), bins - 1)
    records: list[dict[str, Any]] = []
    ece = 0.0
    mce = 0.0
    for bin_id in range(bins):
        selected = bin_ids == bin_id
        mass = float(normalized[selected].sum())
        mean_confidence = (
            float(np.average(confidence[selected], weights=normalized[selected]))
            if mass > 0.0
            else None
        )
        accuracy = (
            float(np.average(correct[selected], weights=normalized[selected]))
            if mass > 0.0
            else None
        )
        gap = (
            abs(float(accuracy) - float(mean_confidence))
            if mass > 0.0
            else None
        )
        if gap is not None:
            ece += mass * gap
            mce = max(mce, gap)
        records.append(
            {
                "bin": bin_id,
                "lower": bin_id / bins,
                "upper": (bin_id + 1) / bins,
                "weighted_mass": mass,
                "mean_confidence": mean_confidence,
                "accuracy": accuracy,
                "absolute_gap": gap,
            }
        )

    classwise: dict[str, float] = {}
    for class_index, label in enumerate(ROSTER):
        scores = probabilities[:, class_index]
        observed = (truth == class_index).astype(np.float64)
        class_bins = np.minimum((scores * bins).astype(np.int64), bins - 1)
        value = 0.0
        for bin_id in range(bins):
            selected = class_bins == bin_id
            mass = float(normalized[selected].sum())
            if mass <= 0.0:
                continue
            predicted_rate = float(np.average(scores[selected], weights=normalized[selected]))
            observed_rate = float(np.average(observed[selected], weights=normalized[selected]))
            value += mass * abs(observed_rate - predicted_rate)
        classwise[label] = value
    return {
        "top_label_ece": ece,
        "top_label_mce": mce,
        "mean_classwise_ece": float(np.mean(list(classwise.values()))),
        "classwise_ece": classwise,
        "bins": bins,
    }, records


def score_subset(
    truth_labels: np.ndarray,
    donors: np.ndarray,
    probabilities: np.ndarray,
) -> dict[str, Any]:
    labels = np.asarray([ROSTER[index] for index in truth_labels], dtype=str)
    predicted_indices = np.argmax(probabilities, axis=1)
    predicted = np.asarray([ROSTER[index] for index in predicted_indices], dtype=str)
    weights = np.asarray(
        donor_class_balanced_weights(donors.tolist(), labels.tolist()),
        dtype=np.float64,
    )
    macro_f1, per_class = weighted_macro_f1(
        labels.tolist(), predicted.tolist(), weights.tolist()
    )
    probability_rows = [
        {label: float(probabilities[row, index]) for index, label in enumerate(ROSTER)}
        for row in range(len(probabilities))
    ]
    brier = multiclass_brier_score(
        labels.tolist(), probability_rows, weights.tolist()
    )
    calibration, bins = calibration_metrics(
        truth_labels, probabilities, weights
    )
    return {
        "macro_f1": float(macro_f1),
        "per_class_f1": {str(key): float(value) for key, value in per_class.items()},
        "multiclass_brier": float(brier),
        "calibration": calibration,
        "calibration_bins": bins,
        "rows": len(truth_labels),
        "donors": len(set(donors.tolist())),
        "classes": sorted(set(labels.tolist())),
    }


def choose_strongest(metrics: Mapping[str, Mapping[str, float]]) -> str:
    if set(metrics) != set(MODEL_IDS):
        raise CellStudyScoreError("strongest-baseline roster differs")
    return min(
        MODEL_IDS,
        key=lambda model_id: (
            -float(metrics[model_id]["donor_class_balanced_macro_f1"]),
            -float(metrics[model_id]["study_balanced_macro_f1"]),
            float(metrics[model_id]["multiclass_brier"]),
            model_id,
        ),
    )


def build_multiplicities(
    donor_studies: np.ndarray,
    present: np.ndarray,
    *,
    n_resamples: int,
    seed: int,
) -> np.ndarray:
    """Study-stratified donor bootstrap with fixed evaluable-class coverage."""

    if donor_studies.ndim != 1 or present.shape[0] != len(donor_studies):
        raise CellStudyScoreError("bootstrap donor strata differ")
    rng = np.random.default_rng(seed)
    result = np.zeros((n_resamples, len(donor_studies)), dtype=np.int16)
    for study in sorted(set(donor_studies.tolist())):
        donor_indices = np.flatnonzero(donor_studies == study)
        class_roster = present[donor_indices].any(axis=0)
        remaining = np.arange(n_resamples)
        attempts = 0
        while len(remaining):
            attempts += 1
            if attempts > 10_000:
                raise CellStudyScoreError("class-preserving donor bootstrap stalled")
            draws = rng.integers(
                0, len(donor_indices), size=(len(remaining), len(donor_indices))
            )
            covered = present[donor_indices[draws]].any(axis=1)
            accepted = np.all(covered[:, class_roster], axis=1)
            accepted_rows = remaining[accepted]
            accepted_draws = draws[accepted]
            if len(accepted_rows):
                row_axis = np.repeat(accepted_rows, len(donor_indices))
                column_axis = donor_indices[accepted_draws.reshape(-1)]
                np.add.at(result, (row_axis, column_axis), 1)
            remaining = remaining[~accepted]
    return result


def sufficient_statistics(
    truth: np.ndarray,
    predicted: np.ndarray,
    probabilities: np.ndarray,
    donor_index: np.ndarray,
    *,
    donors: int,
) -> dict[str, np.ndarray]:
    classes = len(ROSTER)
    present = np.zeros((donors, classes), dtype=bool)
    confusion = np.zeros((donors, classes, classes), dtype=np.float64)
    brier = np.zeros((donors, classes), dtype=np.float64)
    bin_mass = np.zeros((donors, classes, CALIBRATION_BINS), dtype=np.float64)
    bin_confidence = np.zeros_like(bin_mass)
    bin_correct = np.zeros_like(bin_mass)
    confidence = probabilities[np.arange(len(truth)), predicted]
    correct = (truth == predicted).astype(np.float64)
    bin_ids = np.minimum(
        (confidence * CALIBRATION_BINS).astype(np.int64), CALIBRATION_BINS - 1
    )
    losses = np.sum(
        (probabilities - np.eye(classes, dtype=np.float64)[truth]) ** 2,
        axis=1,
    )
    for donor in range(donors):
        for truth_class in range(classes):
            selected = (donor_index == donor) & (truth == truth_class)
            count = int(selected.sum())
            if not count:
                continue
            present[donor, truth_class] = True
            confusion[donor, truth_class] = (
                np.bincount(predicted[selected], minlength=classes) / count
            )
            brier[donor, truth_class] = float(losses[selected].mean())
            for bin_id in range(CALIBRATION_BINS):
                in_bin = selected & (bin_ids == bin_id)
                bin_mass[donor, truth_class, bin_id] = in_bin.sum() / count
                bin_confidence[donor, truth_class, bin_id] = confidence[in_bin].sum() / count
                bin_correct[donor, truth_class, bin_id] = correct[in_bin].sum() / count
    return {
        "present": present,
        "confusion": confusion,
        "brier": brier,
        "bin_mass": bin_mass,
        "bin_confidence": bin_confidence,
        "bin_correct": bin_correct,
    }


def endpoints_from_multiplicities(
    stats: Mapping[str, np.ndarray],
    multiplicities: np.ndarray,
    donor_indices: np.ndarray,
    class_roster: np.ndarray,
) -> dict[str, np.ndarray]:
    mult = multiplicities[:, donor_indices].astype(np.float64)
    present = stats["present"][donor_indices][:, class_roster]
    denominators = mult @ present.astype(np.float64)
    if np.any(denominators <= 0.0):
        raise CellStudyScoreError("bootstrap lost a frozen evaluation class")
    confusion_sum = np.einsum(
        "bd,dcp->bcp",
        mult,
        stats["confusion"][donor_indices][:, class_roster, :],
        optimize=True,
    )
    weighted_confusion = (
        confusion_sum / denominators[:, :, np.newaxis] / int(class_roster.sum())
    )
    target_positions = np.flatnonzero(class_roster)
    true_positive = np.stack(
        [weighted_confusion[:, position, label] for position, label in enumerate(target_positions)],
        axis=1,
    )
    false_positive = np.stack(
        [weighted_confusion[:, :, label].sum(axis=1) for label in target_positions],
        axis=1,
    ) - true_positive
    false_negative = weighted_confusion.sum(axis=2) - true_positive
    denominator = 2.0 * true_positive + false_positive + false_negative
    class_f1 = np.divide(
        2.0 * true_positive,
        denominator,
        out=np.zeros_like(true_positive),
        where=denominator > 0.0,
    )
    brier_sum = mult @ stats["brier"][donor_indices][:, class_roster]
    brier = np.mean(brier_sum / denominators, axis=1)
    confidence_sum = np.einsum(
        "bd,dck->bck",
        mult,
        stats["bin_confidence"][donor_indices][:, class_roster, :],
        optimize=True,
    )
    correct_sum = np.einsum(
        "bd,dck->bck",
        mult,
        stats["bin_correct"][donor_indices][:, class_roster, :],
        optimize=True,
    )
    weighted_confidence = (
        confidence_sum / denominators[:, :, np.newaxis] / int(class_roster.sum())
    ).sum(axis=1)
    weighted_correct = (
        correct_sum / denominators[:, :, np.newaxis] / int(class_roster.sum())
    ).sum(axis=1)
    ece = np.abs(weighted_correct - weighted_confidence).sum(axis=1)
    return {
        "macro_f1": class_f1.mean(axis=1),
        "class_f1": class_f1,
        "brier": brier,
        "top_label_ece": ece,
    }


def interval(values: np.ndarray) -> dict[str, Any]:
    if values.shape != (BOOTSTRAP_RESAMPLES,) or not np.all(np.isfinite(values)):
        raise CellStudyScoreError("bootstrap distribution differs")
    return {
        "lower": float(np.quantile(values, 0.025, method="linear")),
        "median": float(np.quantile(values, 0.5, method="linear")),
        "upper": float(np.quantile(values, 0.975, method="linear")),
        "n_resamples": BOOTSTRAP_RESAMPLES,
        "confidence_level": 0.95,
        "seed": BOOTSTRAP_SEED,
    }


def run(
    *,
    source: Path,
    split: Path,
    predictions: Path,
    output: Path,
    expected_source_artifacts_sha256: str,
    expected_split_artifacts_sha256: str,
    expected_prediction_artifacts_sha256: str,
) -> None:
    if output.exists():
        raise CellStudyScoreError("refusing to overwrite scorer output")
    for root, expected in (
        (source, expected_source_artifacts_sha256),
        (split, expected_split_artifacts_sha256),
        (predictions, expected_prediction_artifacts_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise CellStudyScoreError("input ARTIFACTS SHA-256 differs")
        verify_frozen_tree(root)

    source_manifest = verify_frozen_tree(source)["metadata"]
    split_manifest = verify_frozen_tree(split)["metadata"]
    prediction_manifest = verify_frozen_tree(predictions)["metadata"]
    if (
        source_manifest.get("subset_id") != DATASET_VIEW_ID
        or source_manifest.get("row_count") != EXPECTED_ROWS
        or source_manifest.get("sealed_outcomes_read") is not False
        or source_manifest.get("histology_read") is not False
        or split_manifest.get("split_id") != SPLIT_ID
        or split_manifest.get("rows") != EXPECTED_ROWS
        or split_manifest.get("donors") != EXPECTED_DONORS
        or split_manifest.get("studies") != EXPECTED_STUDIES
        or split_manifest.get("target_labels_used_for_assignment") is not False
        or split_manifest.get("sealed_outcomes_read") is not False
        or split_manifest.get("histology_read") is not False
        or prediction_manifest.get("dataset_view_id") != DATASET_VIEW_ID
        or prediction_manifest.get("split_id") != SPLIT_ID
        or prediction_manifest.get("rows") != EXPECTED_ROWS
        or prediction_manifest.get("metrics_calculated") is not False
        or prediction_manifest.get("sealed_outcomes_read") is not False
    ):
        raise CellStudyScoreError("frozen input metadata differs")

    _, source_rows = read_tsv(source / "selection.tsv")
    split_fields, split_rows = read_tsv(split / "row_outer_folds.tsv")
    if split_fields != ["row_id", "donor_id", "dataset", "outer_fold"]:
        raise CellStudyScoreError("split table schema differs")
    if len(source_rows) != EXPECTED_ROWS or len(split_rows) != EXPECTED_ROWS:
        raise CellStudyScoreError("source or split row count differs")
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
    ):
        raise CellStudyScoreError("source labels or donor-safe split alignment differs")
    label_index = {label: index for index, label in enumerate(ROSTER)}
    truth = np.asarray([label_index[label] for label in labels], dtype=np.int8)

    prediction_receipt = json.loads(
        (predictions / "prediction_receipt.json").read_text(encoding="utf-8")
    )
    if (
        prediction_receipt.get("status") != "pass_development_predictions"
        or prediction_receipt.get("models") != list(MODEL_IDS)
        or prediction_receipt.get("rows") != EXPECTED_ROWS
        or prediction_receipt.get("donors") != EXPECTED_DONORS
        or prediction_receipt.get("studies") != EXPECTED_STUDIES
        or prediction_receipt.get("parameters", {}).get("seeds")
        != list(SCREEN_SEEDS)
        or prediction_receipt.get("parameters", {}).get("final_prediction")
        != "unweighted_mean_probability_ensemble"
        or prediction_receipt.get("metrics_calculated") is not False
        or prediction_receipt.get("prediction_tables_contain_observed_labels") is not False
        or prediction_receipt.get("histology_read") is not False
        or prediction_receipt.get("sealed_outcomes_read") is not False
    ):
        raise CellStudyScoreError("prediction receipt differs")

    prediction_fields = [
        "row_id",
        "donor_id",
        "dataset",
        "outer_fold",
        "predicted_class",
        *(f"probability::{label}" for label in ROSTER),
    ]
    def read_aligned_prediction(path: Path) -> np.ndarray:
        fields, rows = read_tsv(path)
        if fields != prediction_fields or len(rows) != EXPECTED_ROWS:
            raise CellStudyScoreError("standardized prediction table schema differs")
        if any("label" in field.lower() or "histolog" in field.lower() for field in fields):
            raise CellStudyScoreError("prediction table contains evaluator-only fields")
        if (
            [row["row_id"] for row in rows] != row_ids.tolist()
            or [row["donor_id"] for row in rows] != donors.tolist()
            or [row["dataset"] for row in rows] != studies.tolist()
            or [int(row["outer_fold"]) for row in rows] != folds.tolist()
        ):
            raise CellStudyScoreError("prediction rows differ from frozen source and split")
        values = probability_matrix(rows)
        expected_hard = [ROSTER[index] for index in np.argmax(values, axis=1)]
        if [row["predicted_class"] for row in rows] != expected_hard:
            raise CellStudyScoreError("hard predictions differ from frozen probabilities")
        return values

    model_probabilities: dict[str, np.ndarray] = {}
    seed_metrics_rows: list[dict[str, Any]] = []
    for model_id in MODEL_IDS:
        ensemble = read_aligned_prediction(
            predictions / "predictions" / f"{model_id}.tsv"
        )
        seed_values = []
        for seed in SCREEN_SEEDS:
            values = read_aligned_prediction(
                predictions / "predictions" / f"{model_id}--seed-{seed}.tsv"
            )
            seed_values.append(values)
            overall = score_subset(truth, donors, values)
            study_macro = float(
                np.mean(
                    [
                        score_subset(
                            truth[studies == study],
                            donors[studies == study],
                            values[studies == study],
                        )["macro_f1"]
                        for study in sorted(set(studies.tolist()))
                    ]
                )
            )
            seed_metrics_rows.append(
                {
                    "model_id": model_id,
                    "seed": seed,
                    "donor_class_balanced_macro_f1": overall["macro_f1"],
                    "study_balanced_macro_f1": study_macro,
                    "multiclass_brier": overall["multiclass_brier"],
                    "top_label_ece": overall["calibration"]["top_label_ece"],
                }
            )
        if not np.allclose(
            ensemble,
            np.mean(np.stack(seed_values, axis=0), axis=0),
            rtol=0.0,
            atol=1.0e-15,
        ):
            raise CellStudyScoreError("ensemble prediction differs from frozen seed mean")
        model_probabilities[model_id] = ensemble

    unique_donors = sorted(set(donors.tolist()))
    unique_studies = sorted(set(studies.tolist()))
    donor_to_index = {donor: index for index, donor in enumerate(unique_donors)}
    donor_indices = np.asarray([donor_to_index[donor] for donor in donors], dtype=np.int16)
    donor_studies = np.asarray(
        [next(iter(set(studies[donors == donor].tolist()))) for donor in unique_donors],
        dtype=str,
    )
    global_present = np.zeros((EXPECTED_DONORS, len(ROSTER)), dtype=bool)
    for donor in range(EXPECTED_DONORS):
        global_present[donor, np.unique(truth[donor_indices == donor])] = True
    multiplicities = build_multiplicities(
        donor_studies,
        global_present,
        n_resamples=BOOTSTRAP_RESAMPLES,
        seed=BOOTSTRAP_SEED,
    )
    if (
        multiplicities.shape != (BOOTSTRAP_RESAMPLES, EXPECTED_DONORS)
        or np.any(multiplicities.sum(axis=1) != EXPECTED_DONORS)
    ):
        raise CellStudyScoreError("donor bootstrap multiplicities differ")

    output.mkdir(mode=0o750)
    metrics_rows: list[dict[str, Any]] = []
    class_rows: list[dict[str, Any]] = []
    donor_rows: list[dict[str, Any]] = []
    study_rows: list[dict[str, Any]] = []
    calibration_rows: list[dict[str, Any]] = []
    interval_rows: list[dict[str, Any]] = []
    point_metrics: dict[str, dict[str, float]] = {}
    distributions: dict[str, np.ndarray] = {}

    all_donor_indices = np.arange(EXPECTED_DONORS)
    all_classes = np.ones(len(ROSTER), dtype=bool)
    for model_id, probabilities in model_probabilities.items():
        overall = score_subset(truth, donors, probabilities)
        study_scores: dict[str, dict[str, Any]] = {}
        for study in unique_studies:
            selected = studies == study
            value = score_subset(truth[selected], donors[selected], probabilities[selected])
            study_scores[study] = value
            study_rows.append(
                {
                    "model_id": model_id,
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
        study_balanced_macro = float(
            np.mean([value["macro_f1"] for value in study_scores.values()])
        )
        point_metrics[model_id] = {
            "donor_class_balanced_macro_f1": overall["macro_f1"],
            "study_balanced_macro_f1": study_balanced_macro,
            "multiclass_brier": overall["multiclass_brier"],
            "top_label_ece": overall["calibration"]["top_label_ece"],
            "mean_classwise_ece": overall["calibration"]["mean_classwise_ece"],
        }
        metrics_rows.append({"model_id": model_id, **point_metrics[model_id]})
        for label in ROSTER:
            study_values = [
                value["per_class_f1"][label]
                for value in study_scores.values()
                if label in value["per_class_f1"]
            ]
            class_rows.append(
                {
                    "model_id": model_id,
                    "class": label,
                    "donor_class_balanced_f1": overall["per_class_f1"][label],
                    "study_balanced_f1": float(np.mean(study_values)),
                    "evaluable_studies": len(study_values),
                }
            )
        for record in overall["calibration_bins"]:
            calibration_rows.append(
                {"model_id": model_id, "scope": "all_studies", **record}
            )
        for donor in unique_donors:
            selected = donors == donor
            value = score_subset(truth[selected], donors[selected], probabilities[selected])
            donor_rows.append(
                {
                    "model_id": model_id,
                    "donor_id": donor,
                    "study": studies[selected][0],
                    "rows": value["rows"],
                    "evaluable_classes": ",".join(value["classes"]),
                    "class_balanced_macro_f1": value["macro_f1"],
                    "multiclass_brier": value["multiclass_brier"],
                    "top_label_ece": value["calibration"]["top_label_ece"],
                }
            )

        stats = sufficient_statistics(
            truth,
            np.argmax(probabilities, axis=1),
            probabilities,
            donor_indices,
            donors=EXPECTED_DONORS,
        )
        identity_multiplicities = np.ones((1, EXPECTED_DONORS), dtype=np.int16)
        identity = endpoints_from_multiplicities(
            stats, identity_multiplicities, all_donor_indices, all_classes
        )
        identity_study_parts = []
        for study in unique_studies:
            selected_donors = np.flatnonzero(donor_studies == study)
            study_classes = stats["present"][selected_donors].any(axis=0)
            identity_study_parts.append(
                endpoints_from_multiplicities(
                    stats,
                    identity_multiplicities,
                    selected_donors,
                    study_classes,
                )["macro_f1"][0]
            )
        if (
            not math.isclose(
                float(identity["macro_f1"][0]),
                overall["macro_f1"],
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            or not math.isclose(
                float(identity["brier"][0]),
                overall["multiclass_brier"],
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            or not math.isclose(
                float(identity["top_label_ece"][0]),
                overall["calibration"]["top_label_ece"],
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
            or not math.isclose(
                float(np.mean(identity_study_parts)),
                study_balanced_macro,
                rel_tol=0.0,
                abs_tol=1.0e-12,
            )
        ):
            raise CellStudyScoreError(
                "sufficient-statistic endpoint rederivation differs"
            )
        donor_bootstrap = endpoints_from_multiplicities(
            stats, multiplicities, all_donor_indices, all_classes
        )
        study_macro_parts = []
        for study in unique_studies:
            selected_donors = np.flatnonzero(donor_studies == study)
            study_classes = stats["present"][selected_donors].any(axis=0)
            study_macro_parts.append(
                endpoints_from_multiplicities(
                    stats, multiplicities, selected_donors, study_classes
                )["macro_f1"]
            )
        study_bootstrap = np.mean(np.stack(study_macro_parts, axis=1), axis=1)
        metric_distributions = {
            "donor_class_balanced_macro_f1": donor_bootstrap["macro_f1"],
            "study_balanced_macro_f1": study_bootstrap,
            "multiclass_brier": donor_bootstrap["brier"],
            "top_label_ece": donor_bootstrap["top_label_ece"],
        }
        for class_index, label in enumerate(ROSTER):
            metric_distributions[f"class_f1::{label}"] = donor_bootstrap["class_f1"][:, class_index]
        for metric, values in metric_distributions.items():
            key = f"{model_id}::{metric}"
            distributions[key] = values.astype(np.float64, copy=False)
            observed = (
                overall["per_class_f1"][metric.split("::", 1)[1]]
                if metric.startswith("class_f1::")
                else point_metrics[model_id][metric]
            )
            interval_rows.append(
                {"model_id": model_id, "metric": metric, "estimate": observed, **interval(values)}
            )

    strongest = choose_strongest(point_metrics)
    delta_rows: list[dict[str, Any]] = []
    for comparator in MODEL_IDS:
        if comparator == strongest:
            continue
        for metric, higher_is_better in (
            ("donor_class_balanced_macro_f1", True),
            ("study_balanced_macro_f1", True),
            ("multiclass_brier", False),
        ):
            left = distributions[f"{strongest}::{metric}"]
            right = distributions[f"{comparator}::{metric}"]
            values = left - right if higher_is_better else right - left
            observed = (
                point_metrics[strongest][metric] - point_metrics[comparator][metric]
                if higher_is_better
                else point_metrics[comparator][metric] - point_metrics[strongest][metric]
            )
            delta_rows.append(
                {
                    "strongest_model_id": strongest,
                    "comparator_model_id": comparator,
                    "metric": metric,
                    "oriented_estimate": observed,
                    **interval(values),
                    "probability_improvement": float(np.mean(values > 0.0)),
                    "selection_conditioned_descriptive_only": True,
                }
            )

    write_tsv(output / "model_metrics.tsv", list(metrics_rows[0]), metrics_rows)
    write_tsv(
        output / "seed_metrics.tsv", list(seed_metrics_rows[0]), seed_metrics_rows
    )
    write_tsv(output / "per_class_f1.tsv", list(class_rows[0]), class_rows)
    write_tsv(output / "donor_summaries.tsv", list(donor_rows[0]), donor_rows)
    write_tsv(output / "study_summaries.tsv", list(study_rows[0]), study_rows)
    write_tsv(output / "calibration_bins.tsv", list(calibration_rows[0]), calibration_rows)
    write_tsv(output / "bootstrap_intervals.tsv", list(interval_rows[0]), interval_rows)
    write_tsv(output / "strongest_vs_baselines.tsv", list(delta_rows[0]), delta_rows)
    with (output / "bootstrap_distributions.npz").open("xb") as handle:
        np.savez_compressed(handle, **distributions)

    strongest_value = {
        "strongest_task_native_baseline": strongest,
        "selection_rule": "maximum donor-class-balanced macro-F1, then maximum study-balanced macro-F1, then minimum multiclass Brier, then ascending model ID",
        "metrics": point_metrics[strongest],
        "selection_conditioned_uncertainty_is_descriptive": True,
    }
    write_json_exclusive(output / "strongest_baseline.json", strongest_value)
    receipt = {
        "schema_version": "masld-bench-cell-study-heldout-evaluator-v1",
        "status": "pass_development_evaluation_only",
        "evaluator_id": "cell_study_50000_evaluation_only_v1",
        "dataset_view_id": DATASET_VIEW_ID,
        "split_id": SPLIT_ID,
        "rows": EXPECTED_ROWS,
        "donors": EXPECTED_DONORS,
        "studies": EXPECTED_STUDIES,
        "classes": list(ROSTER),
        "models": list(MODEL_IDS),
        "screen_seeds": list(SCREEN_SEEDS),
        "final_prediction": "unweighted_mean_probability_ensemble",
        "strongest_task_native_baseline": strongest,
        "source_artifacts_sha256": expected_source_artifacts_sha256,
        "split_artifacts_sha256": expected_split_artifacts_sha256,
        "prediction_artifacts_sha256": expected_prediction_artifacts_sha256,
        "source_labels_read_inside_evaluator": ["source_cell_type", "broad_label"],
        "source_cell_type_used_for_scoring": False,
        "prediction_tables_contain_observed_labels": False,
        "models_fit": False,
        "preprocessing_fit": False,
        "calibration_fit": False,
        "calibration_evaluation": "fixed_10_bin_top_label_and_one_vs_rest_classwise_ece",
        "metric_definitions": {
            "donor_class_balanced_macro_f1": "Every class receives equal mass, every donor containing that class receives equal within-class mass, and cells divide their donor-class share.",
            "study_balanced_macro_f1": "Arithmetic mean of seven study-specific donor-class-balanced macro-F1 values; each study receives equal mass.",
            "per_class_f1": "Global donor-class-balanced class F1 plus the equal-study mean among studies where the source class is observed.",
            "multiclass_brier": "Donor-class-balanced sum of squared error over the frozen five-class probability vector.",
            "calibration": "Fixed ten-bin top-label ECE/MCE and one-vs-rest classwise ECE; no calibrator is fit.",
            "strongest_baseline": "Maximum donor-class-balanced macro-F1, then maximum study-balanced macro-F1, then minimum Brier, then ascending model ID.",
        },
        "bootstrap": {
            "method": "study_stratified_class_coverage_preserving_donor_cluster_bootstrap",
            "resamples": BOOTSTRAP_RESAMPLES,
            "seed": BOOTSTRAP_SEED,
            "donors_are_independent_units": True,
            "cells_are_independent_units": False,
            "study_donor_counts_preserved": True,
            "study_evaluable_class_rosters_preserved": True,
            "distribution_sha256": canonical_sha256(
                {key: values.tolist() for key, values in sorted(distributions.items())}
            ),
        },
        "histology_read": False,
        "sealed_outcomes_read": False,
        "clinical_claim_allowed": False,
    }
    write_json_exclusive(output / "evaluation_receipt.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "cell_state_study_heldout_evaluation",
            "dataset_view_id": DATASET_VIEW_ID,
            "split_id": SPLIT_ID,
            "rows": EXPECTED_ROWS,
            "donors": EXPECTED_DONORS,
            "studies": EXPECTED_STUDIES,
            "models_fit": False,
            "histology_read": False,
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
    )
    print(json.dumps({"output": arguments.output.resolve().as_posix(), "status": "passed"}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
