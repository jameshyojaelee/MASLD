#!/usr/bin/env python3
"""Independently score frozen GSE260666 three-state external predictions."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np


CLASSES = ("NOR", "NAFL", "NASH")
EXPECTED_SOURCE_TO_EVALUATION = {
    "healthy_control": "NOR",
    "non-alcoholic_fatty_liver_disease_(nafld)": "NAFL",
    "non-alcoholic_steatohepatitis_(nash)": "NASH",
}
EXPECTED_SOURCE_COUNTS = {
    "healthy_control": 6,
    "non-alcoholic_fatty_liver_disease_(nafld)": 6,
    "non-alcoholic_steatohepatitis_(nash)": 4,
}
PROBABILITY_FIELDS = tuple(f"probability_{value}" for value in CLASSES)
MANDATORY_BASELINES = (
    "training_stage_distribution",
    "rna_hvg_pca_nearest_centroid",
    "rna_hvg_pca_elastic_net",
    "rna_hvg_pca_linear_svm",
)


class ExternalTransferEvaluatorError(RuntimeError):
    """Raised when predictions or outcomes does not meet the frozen evaluator requirements."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ExternalTransferEvaluatorError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def write_tsv(
    path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]
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


def _class_statistics(observed: np.ndarray, predicted: np.ndarray) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for class_id, label in enumerate(CLASSES):
        true = observed == class_id
        called = predicted == class_id
        tp = int(np.sum(true & called))
        fp = int(np.sum(~true & called))
        fn = int(np.sum(true & ~called))
        precision = 0.0 if tp + fp == 0 else tp / (tp + fp)
        recall = 0.0 if tp + fn == 0 else tp / (tp + fn)
        f1 = 0.0 if 2 * tp + fp + fn == 0 else 2 * tp / (2 * tp + fp + fn)
        rows.append(
            {
                "class": label,
                "support": int(np.sum(true)),
                "precision": precision,
                "recall": recall,
                "f1": f1,
                "tp": tp,
                "fp": fp,
                "fn": fn,
            }
        )
    return rows


def _ece(observed: np.ndarray, predicted: np.ndarray, probabilities: np.ndarray) -> float:
    confidence = probabilities.max(axis=1)
    correct = (observed == predicted).astype(np.float64)
    edges = np.linspace(0.0, 1.0, 5)
    total = 0.0
    for index in range(4):
        member = (confidence >= edges[index]) & (
            confidence <= edges[index + 1] if index == 3 else confidence < edges[index + 1]
        )
        if np.any(member):
            total += float(np.mean(member)) * abs(
                float(np.mean(confidence[member]) - np.mean(correct[member]))
            )
    return total


def metrics(observed: np.ndarray, probabilities: np.ndarray) -> dict[str, float]:
    predicted = np.argmax(probabilities, axis=1)
    per_class = _class_statistics(observed, predicted)
    one_hot = np.eye(len(CLASSES), dtype=np.float64)[observed]
    confidence = probabilities.max(axis=1)
    correct = (observed == predicted).astype(np.float64)
    clipped = np.clip(probabilities[np.arange(len(observed)), observed], 1e-15, 1.0)
    return {
        "macro_f1": float(np.mean([row["f1"] for row in per_class])),
        "accuracy": float(np.mean(correct)),
        "multiclass_brier": float(np.mean(np.sum((probabilities - one_hot) ** 2, axis=1))),
        "log_loss": float(-np.mean(np.log(clipped))),
        "confidence_calibration_gap": float(np.mean(confidence) - np.mean(correct)),
        "ece4_descriptive": _ece(observed, predicted, probabilities),
    }


def _load_predictions(path: Path, row_ids: Sequence[str]) -> tuple[str, np.ndarray]:
    fields, rows = read_tsv(path)
    required = {"row_id", "model_id", "predicted_stage3", *PROBABILITY_FIELDS}
    if not required <= set(fields) or len(rows) != len(row_ids):
        raise ExternalTransferEvaluatorError("prediction schema or row count differs")
    model_ids = {row["model_id"] for row in rows}
    if len(model_ids) != 1:
        raise ExternalTransferEvaluatorError("prediction bundle contains multiple models")
    by_id = {row["row_id"]: row for row in rows}
    if len(by_id) != len(rows) or set(by_id) != set(row_ids):
        raise ExternalTransferEvaluatorError("prediction row IDs differ")
    ordered = [by_id[value] for value in row_ids]
    probabilities = np.asarray(
        [[float(row[field]) for field in PROBABILITY_FIELDS] for row in ordered],
        dtype=np.float64,
    )
    if (
        not np.all(np.isfinite(probabilities))
        or np.any(probabilities < 0.0)
        or np.any(probabilities > 1.0)
        or not np.allclose(probabilities.sum(axis=1), 1.0, rtol=0.0, atol=1e-8)
    ):
        raise ExternalTransferEvaluatorError("prediction probabilities are invalid")
    expected = [CLASSES[index] for index in np.argmax(probabilities, axis=1)]
    if expected != [row["predicted_stage3"] for row in ordered]:
        raise ExternalTransferEvaluatorError("predicted labels differ from fixed-order argmax")
    return next(iter(model_ids)), probabilities


def evaluate(
    *,
    labels_path: Path,
    prediction_paths: Sequence[Path],
    output: Path,
    bootstrap_replicates: int = 10000,
    bootstrap_seed: int = 20260824,
    mandatory_baselines: Sequence[str] = MANDATORY_BASELINES,
    prediction_bundle_artifacts_sha256: Sequence[str] | None = None,
) -> dict[str, Any]:
    if output.exists():
        raise ExternalTransferEvaluatorError(f"refusing to overwrite evaluation: {output}")
    label_fields, labels = read_tsv(labels_path)
    required_labels = {"row_id", "source_label", "evaluation_stage3"}
    if not required_labels <= set(label_fields) or len(labels) != 16:
        raise ExternalTransferEvaluatorError("evaluator label schema or participant count differs")
    row_ids = [row["row_id"] for row in labels]
    if len(row_ids) != len(set(row_ids)):
        raise ExternalTransferEvaluatorError("evaluator row ID is duplicated")
    if any(row["evaluation_stage3"] not in CLASSES for row in labels):
        raise ExternalTransferEvaluatorError("evaluation class differs")
    source_to_evaluation = {}
    for row in labels:
        prior = source_to_evaluation.setdefault(row["source_label"], row["evaluation_stage3"])
        if prior != row["evaluation_stage3"]:
            raise ExternalTransferEvaluatorError("source mapping is not one-to-one")
    if source_to_evaluation != EXPECTED_SOURCE_TO_EVALUATION:
        raise ExternalTransferEvaluatorError("external source labels were pooled or relabeled")
    if Counter(row["source_label"] for row in labels) != Counter(EXPECTED_SOURCE_COUNTS):
        raise ExternalTransferEvaluatorError("external source-label counts differ from 6/6/4")
    observed = np.asarray([CLASSES.index(row["evaluation_stage3"]) for row in labels], dtype=int)
    model_probabilities: dict[str, np.ndarray] = {}
    if prediction_bundle_artifacts_sha256 is not None and len(
        prediction_bundle_artifacts_sha256
    ) != len(prediction_paths):
        raise ExternalTransferEvaluatorError("prediction artifact SHA roster differs")
    prediction_file_sha256: dict[str, str] = {}
    prediction_artifacts_sha256: dict[str, str] = {}
    for path_index, path in enumerate(prediction_paths):
        model_id, probabilities = _load_predictions(path, row_ids)
        if model_id in model_probabilities:
            raise ExternalTransferEvaluatorError("model prediction bundle is duplicated")
        model_probabilities[model_id] = probabilities
        prediction_file_sha256[model_id] = sha256_file(path)
        if prediction_bundle_artifacts_sha256 is not None:
            prediction_artifacts_sha256[model_id] = prediction_bundle_artifacts_sha256[
                path_index
            ]
    missing = sorted(set(mandatory_baselines) - set(model_probabilities))
    if missing:
        raise ExternalTransferEvaluatorError(f"mandatory external baselines are missing: {missing}")
    output.mkdir(parents=True)
    metric_rows: list[dict[str, Any]] = []
    class_rows: list[dict[str, Any]] = []
    point_metrics: dict[str, dict[str, float]] = {}
    point_per_class: dict[str, dict[str, float]] = {}
    confusion: dict[str, Any] = {}
    for model_id, probabilities in sorted(model_probabilities.items()):
        values = metrics(observed, probabilities)
        point_metrics[model_id] = values
        metric_rows.extend(
            {"model_id": model_id, "metric": key, "estimate": value}
            for key, value in values.items()
        )
        predicted = np.argmax(probabilities, axis=1)
        model_class = _class_statistics(observed, predicted)
        point_per_class[model_id] = {row["class"]: float(row["f1"]) for row in model_class}
        class_rows.extend({"model_id": model_id, **row} for row in model_class)
        confusion[model_id] = {
            observed_label: {
                predicted_label: int(
                    np.sum((observed == observed_index) & (predicted == predicted_index))
                )
                for predicted_index, predicted_label in enumerate(CLASSES)
            }
            for observed_index, observed_label in enumerate(CLASSES)
        }
    strongest_baseline = min(
        mandatory_baselines,
        key=lambda value: (
            -point_metrics[value]["macro_f1"],
            point_metrics[value]["multiclass_brier"],
            value,
        ),
    )
    rng = np.random.default_rng(bootstrap_seed)
    by_class = [np.flatnonzero(observed == class_id) for class_id in range(len(CLASSES))]
    bootstrap_values = {
        model_id: {metric_id: [] for metric_id in point_metrics[model_id]}
        for model_id in model_probabilities
    }
    bootstrap_class_f1 = {
        model_id: {class_id: [] for class_id in CLASSES}
        for model_id in model_probabilities
    }
    gain_values = {
        model_id: []
        for model_id in model_probabilities
        if model_id not in set(mandatory_baselines)
    }
    for _ in range(bootstrap_replicates):
        sampled = np.concatenate(
            [rng.choice(indices, size=len(indices), replace=True) for indices in by_class]
        )
        baseline_metric = metrics(observed[sampled], model_probabilities[strongest_baseline][sampled])
        for model_id, probabilities in model_probabilities.items():
            value = metrics(observed[sampled], probabilities[sampled])
            for metric_id, estimate in value.items():
                bootstrap_values[model_id][metric_id].append(estimate)
            sampled_predicted = np.argmax(probabilities[sampled], axis=1)
            for row in _class_statistics(observed[sampled], sampled_predicted):
                bootstrap_class_f1[model_id][row["class"]].append(float(row["f1"]))
            if model_id in gain_values:
                gain_values[model_id].append(value["macro_f1"] - baseline_metric["macro_f1"])
    bootstrap_rows: list[dict[str, Any]] = []
    for model_id, metrics_by_id in sorted(bootstrap_values.items()):
        for metric_id, values in metrics_by_id.items():
            vector = np.asarray(values, dtype=np.float64)
            bootstrap_rows.append(
                {
                    "model_id": model_id,
                    "metric": metric_id,
                    "estimate": point_metrics[model_id][metric_id],
                    "ci_low": float(np.quantile(vector, 0.025)),
                    "ci_high": float(np.quantile(vector, 0.975)),
                }
            )
    gain_rows = []
    for model_id, values in sorted(gain_values.items()):
        vector = np.asarray(values, dtype=np.float64)
        point = point_metrics[model_id]["macro_f1"] - point_metrics[strongest_baseline]["macro_f1"]
        gain_rows.append(
            {
                "model_id": model_id,
                "baseline_model_id": strongest_baseline,
                "macro_f1_gain": point,
                "ci_low": float(np.quantile(vector, 0.025)),
                "ci_high": float(np.quantile(vector, 0.975)),
            }
        )
    class_bootstrap_rows: list[dict[str, Any]] = []
    for model_id, by_class_id in sorted(bootstrap_class_f1.items()):
        for class_id, values in by_class_id.items():
            vector = np.asarray(values, dtype=np.float64)
            class_bootstrap_rows.append(
                {
                    "model_id": model_id,
                    "class": class_id,
                    "f1": point_per_class[model_id][class_id],
                    "ci_low": float(np.quantile(vector, 0.025)),
                    "ci_high": float(np.quantile(vector, 0.975)),
                }
            )
    write_tsv(output / "metrics.tsv", ("model_id", "metric", "estimate"), metric_rows)
    write_tsv(
        output / "per_class.tsv",
        ("model_id", "class", "support", "precision", "recall", "f1", "tp", "fp", "fn"),
        class_rows,
    )
    write_tsv(
        output / "bootstrap_intervals.tsv",
        ("model_id", "metric", "estimate", "ci_low", "ci_high"),
        bootstrap_rows,
    )
    write_tsv(
        output / "per_class_f1_bootstrap_intervals.tsv",
        ("model_id", "class", "f1", "ci_low", "ci_high"),
        class_bootstrap_rows,
    )
    if gain_rows:
        write_tsv(
            output / "candidate_vs_strongest_baseline.tsv",
            ("model_id", "baseline_model_id", "macro_f1_gain", "ci_low", "ci_high"),
            gain_rows,
        )
    (output / "confusion_matrices.json").write_text(
        json.dumps(confusion, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    receipt = {
        "schema_version": "masld-bench-gse260666-external-evaluation-v1",
        "status": "passed_external_development_scoring",
        "participants": 16,
        "biological_unit": "participant",
        "models": sorted(model_probabilities),
        "mandatory_baselines": list(mandatory_baselines),
        "strongest_baseline": strongest_baseline,
        "strongest_baseline_role": "post_hoc_external_development_reporting_only_not_source_model_selection",
        "primary_metric": "class_balanced_macro_f1",
        "bootstrap_method": "source_class_stratified_participant_bootstrap",
        "bootstrap_replicates": bootstrap_replicates,
        "bootstrap_seed": bootstrap_seed,
        "source_labels_retained": sorted(source_to_evaluation),
        "source_to_evaluation_mapping": source_to_evaluation,
        "nafl_nash_pooled": False,
        "query_fit_or_calibration_performed": False,
        "model_selected_from_external_outcomes": False,
        "external_outcomes_used_for_model_repair": False,
        "nominal_or_confirmatory_p_values_computed": False,
        "prediction_file_sha256": dict(sorted(prediction_file_sha256.items())),
        "prediction_bundle_artifacts_sha256": dict(
            sorted(prediction_artifacts_sha256.items())
        ),
        "project_sealed": False,
        "external_development_only": True,
        "champion_claim_eligible": False,
        "clinical_claim_eligible": False,
        "diagnostic_or_prognostic_claim_eligible": False,
        "confirmatory_claim_eligible": False,
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evaluator-only", required=True, type=Path)
    parser.add_argument("--evaluator-only-artifacts-sha256", required=True)
    parser.add_argument(
        "--prediction",
        action="append",
        required=True,
        help="Frozen prediction root; must contain ARTIFACTS.json and predictions.tsv",
    )
    parser.add_argument("--prediction-artifacts-sha256", action="append", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    if sha256_file(arguments.evaluator_only / "ARTIFACTS.json") != arguments.evaluator_only_artifacts_sha256:
        raise ExternalTransferEvaluatorError("evaluator-only ARTIFACTS SHA differs")
    prediction_roots = [Path(value) for value in arguments.prediction]
    if len(prediction_roots) != len(arguments.prediction_artifacts_sha256):
        raise ExternalTransferEvaluatorError("prediction roots and SHA rosters differ")
    for root, expected in zip(
        prediction_roots, arguments.prediction_artifacts_sha256, strict=True
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise ExternalTransferEvaluatorError(f"prediction ARTIFACTS SHA differs: {root}")
    result = evaluate(
        labels_path=arguments.evaluator_only / "labels.tsv",
        prediction_paths=[root / "predictions.tsv" for root in prediction_roots],
        output=arguments.output,
        prediction_bundle_artifacts_sha256=arguments.prediction_artifacts_sha256,
    )
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
