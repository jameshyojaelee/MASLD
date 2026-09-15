#!/usr/bin/env python3
"""Evaluate committed GSE296875 cell-state predictions by biological donor."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
from pathlib import Path
from typing import Any

import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import sha256_file
from scripts.score_cell_baselines_study_50000 import (
    calibration_metrics,
    endpoints_from_multiplicities,
    score_subset,
    sufficient_statistics,
)


VIEW_ID = "gse296875_rna_cell_state_external_development_7500_v1"
ROSTER = (
    "cholangiocyte",
    "endothelial",
    "hepatocyte",
    "immune",
    "mesenchymal_stromal",
)
BASELINE_MODELS = (
    "hvg_pca_nearest_centroid",
    "hvg_pca_knn",
    "hvg_pca_logistic",
    "hvg_pca_elastic_net",
    "hvg_pca_linear_svm",
)
TF_MODELS = (
    "transcriptformer_tf_sapiens__linear",
    "transcriptformer_tf_sapiens__two_layer_mlp",
)
BOOTSTRAP_RESAMPLES = 10_000
BOOTSTRAP_SEED = 20260825


class TransferEvaluationError(ValueError):
    """Raised when the isolated evaluator requirement differs."""


def _read_tsv(path: Path, *, compressed: bool = False) -> tuple[list[str], list[dict[str, str]]]:
    opener = gzip.open if compressed else open
    with opener(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return list(reader.fieldnames or ()), list(reader)


def _prediction(path: Path, expected_sha256: str, row_ids: list[str]) -> np.ndarray:
    if sha256_file(path) != expected_sha256:
        raise TransferEvaluationError("locked prediction file SHA-256 differs")
    fields, rows = _read_tsv(path)
    expected = ["row_id", "predicted_class", *(f"probability::{label}" for label in ROSTER)]
    if fields != expected or [row["row_id"] for row in rows] != row_ids:
        raise TransferEvaluationError("locked prediction row contract differs")
    values = np.asarray(
        [[float(row[f"probability::{label}"]) for label in ROSTER] for row in rows],
        dtype=np.float64,
    )
    if (
        values.shape != (7_500, 5)
        or not np.all(np.isfinite(values))
        or np.any(values < 0.0)
        or not np.allclose(values.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-7)
        or [row["predicted_class"] for row in rows]
        != [ROSTER[index] for index in np.argmax(values, axis=1)]
    ):
        raise TransferEvaluationError("locked prediction probabilities differ")
    return values


def _interval(values: np.ndarray) -> dict[str, Any]:
    if values.shape != (BOOTSTRAP_RESAMPLES,) or not np.all(np.isfinite(values)):
        raise TransferEvaluationError("bootstrap distribution differs")
    return {
        "lower": float(np.quantile(values, 0.025, method="linear")),
        "median": float(np.quantile(values, 0.5, method="linear")),
        "upper": float(np.quantile(values, 0.975, method="linear")),
        "n_resamples": BOOTSTRAP_RESAMPLES,
        "confidence_level": 0.95,
        "seed": BOOTSTRAP_SEED,
    }


def run(args: argparse.Namespace) -> dict[str, Any]:
    if args.output.exists():
        raise TransferEvaluationError("refusing to overwrite transfer evaluation")
    lock = args.prediction_lock.resolve(strict=True)
    evaluator = args.evaluator.resolve(strict=True)
    for root, expected in (
        (lock, args.expected_prediction_lock_sha256),
        (evaluator, args.expected_evaluator_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise TransferEvaluationError(f"input ARTIFACTS SHA-256 differs: {root}")
        verify_frozen_tree(root)
    lock_meta = verify_frozen_tree(lock)["metadata"]
    evaluator_meta = verify_frozen_tree(evaluator)["metadata"]
    if (
        lock_meta.get("artifact_class") != "gse296875_cell_state_transfer_prediction_lock"
        or lock_meta.get("view_id") != VIEW_ID
        or lock_meta.get("evaluator_labels_read") is not False
        or lock_meta.get("metrics_calculated") is not False
        or evaluator_meta.get("artifact_class")
        != "gse296875_rna_cell_state_external_development_evaluator"
        or evaluator_meta.get("view_id") != VIEW_ID
        or evaluator_meta.get("rows") != 7_500
        or evaluator_meta.get("donors") != 39
        or evaluator_meta.get("champion_claim_allowed") is not False
    ):
        raise TransferEvaluationError("lock or evaluator authority differs")
    selection = json.loads((lock / "selection_lock.json").read_text())
    expected_roster = [*BASELINE_MODELS, *TF_MODELS]
    if (
        selection.get("status") != "locked_before_evaluator_join"
        or selection.get("candidate_roster") != expected_roster
        or selection.get("evaluator_labels_read") is not False
        or selection.get("external_or_champion_claim_allowed") is not False
    ):
        raise TransferEvaluationError("selection lock differs")

    fields, labels = _read_tsv(evaluator / "labels.tsv.gz", compressed=True)
    if fields != ["row_id", "donor_id", "broad_label"] or len(labels) != 7_500:
        raise TransferEvaluationError("evaluator label schema differs")
    row_ids = [row["row_id"] for row in labels]
    donors = np.asarray([row["donor_id"] for row in labels], dtype=str)
    observed_labels = np.asarray([row["broad_label"] for row in labels], dtype=str)
    if (
        len(set(row_ids)) != 7_500
        or len(set(donors)) != 39
        or set(observed_labels) != set(ROSTER)
        or any(set(observed_labels[donors == donor]) != set(ROSTER) for donor in set(donors))
    ):
        raise TransferEvaluationError("evaluator donor/class contract differs")
    truth = np.asarray([ROSTER.index(label) for label in observed_labels], dtype=np.int8)

    probabilities: dict[str, np.ndarray] = {}
    for model_id in expected_roster:
        file_lock = selection["prediction_files"][model_id]
        probabilities[model_id] = _prediction(
            Path(file_lock["path"]), file_lock["sha256"], row_ids
        )

    unique_donors = sorted(set(donors.tolist()))
    donor_to_index = {donor: index for index, donor in enumerate(unique_donors)}
    donor_index = np.asarray([donor_to_index[donor] for donor in donors], dtype=np.int16)
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    draws = rng.integers(0, len(unique_donors), size=(BOOTSTRAP_RESAMPLES, len(unique_donors)))
    multiplicities = np.zeros((BOOTSTRAP_RESAMPLES, len(unique_donors)), dtype=np.int16)
    np.add.at(
        multiplicities,
        (np.repeat(np.arange(BOOTSTRAP_RESAMPLES), len(unique_donors)), draws.reshape(-1)),
        1,
    )
    if np.any(multiplicities.sum(axis=1) != 39):
        raise TransferEvaluationError("donor bootstrap multiplicity differs")

    point: dict[str, dict[str, Any]] = {}
    distributions: dict[str, dict[str, np.ndarray]] = {}
    for model_id, values in probabilities.items():
        result = score_subset(truth, donors, values)
        point[model_id] = {
            "donor_class_balanced_macro_f1": result["macro_f1"],
            "per_class_f1": result["per_class_f1"],
            "multiclass_brier": result["multiclass_brier"],
            "top_label_ece": result["calibration"]["top_label_ece"],
            "mean_classwise_ece": result["calibration"]["mean_classwise_ece"],
        }
        stats = sufficient_statistics(
            truth, np.argmax(values, axis=1), values, donor_index, donors=39
        )
        distributions[model_id] = endpoints_from_multiplicities(
            stats,
            multiplicities,
            np.arange(39),
            np.ones(len(ROSTER), dtype=bool),
        )
    strongest = min(
        BASELINE_MODELS,
        key=lambda model_id: (
            -float(point[model_id]["donor_class_balanced_macro_f1"]),
            float(point[model_id]["multiclass_brier"]),
            model_id,
        ),
    )
    comparisons = {}
    for model_id in TF_MODELS:
        delta_f1 = distributions[model_id]["macro_f1"] - distributions[strongest]["macro_f1"]
        delta_brier = distributions[model_id]["brier"] - distributions[strongest]["brier"]
        delta_ece = distributions[model_id]["top_label_ece"] - distributions[strongest]["top_label_ece"]
        comparison = {
            "model_id": model_id,
            "strongest_task_native_baseline": strongest,
            "selection_rule": selection["strongest_baseline_selection_rule"],
            "macro_f1_gain": float(
                point[model_id]["donor_class_balanced_macro_f1"]
                - point[strongest]["donor_class_balanced_macro_f1"]
            ),
            "macro_f1_gain_interval": _interval(delta_f1),
            "brier_difference": float(
                point[model_id]["multiclass_brier"] - point[strongest]["multiclass_brier"]
            ),
            "brier_difference_interval": _interval(delta_brier),
            "top_label_ece_difference": float(
                point[model_id]["top_label_ece"] - point[strongest]["top_label_ece"]
            ),
            "top_label_ece_difference_interval": _interval(delta_ece),
            "development_transfer_gate": {
                "macro_f1_at_least_0_70": point[model_id]["donor_class_balanced_macro_f1"] >= 0.70,
                "macro_f1_gain_at_least_0_02": float(
                    point[model_id]["donor_class_balanced_macro_f1"]
                    - point[strongest]["donor_class_balanced_macro_f1"]
                ) >= 0.02,
                "paired_95pct_lower_bound_above_zero": _interval(delta_f1)["lower"] > 0.0,
                "every_class_f1_at_least_0_50": min(point[model_id]["per_class_f1"].values()) >= 0.50,
                "brier_not_worse_by_more_than_0_01": float(
                    point[model_id]["multiclass_brier"] - point[strongest]["multiclass_brier"]
                ) <= 0.01,
            },
            "confirmatory_or_champion_gate": False,
        }
        comparison["development_transfer_gate"]["all_components_pass"] = all(
            comparison["development_transfer_gate"].values()
        )
        comparisons[model_id] = comparison

    args.output.mkdir(mode=0o750)
    write_json_exclusive(args.output / "point_metrics.json", point)
    write_json_exclusive(args.output / "comparisons.json", comparisons)
    with (args.output / "bootstrap_distributions.npz").open("xb") as handle:
        np.savez_compressed(
            handle,
            multiplicities=multiplicities,
            **{
                f"{model_id}::{endpoint}": values
                for model_id, endpoints in distributions.items()
                for endpoint, values in endpoints.items()
            },
        )
    receipt = {
        "schema_version": "masld-bench-gse296875-cell-state-evaluation-v1",
        "status": "pass_project_exposed_development_evaluation",
        "view_id": VIEW_ID,
        "rows": 7_500,
        "donors": 39,
        "classes": list(ROSTER),
        "biological_replication_unit": "donor",
        "primary_endpoint": "donor_class_balanced_macro_f1",
        "uncertainty": "10000_paired_donor_cluster_bootstrap_resamples",
        "candidate_roster": expected_roster,
        "strongest_task_native_baseline": strongest,
        "prediction_lock_artifacts_sha256": args.expected_prediction_lock_sha256,
        "evaluator_artifacts_sha256": args.expected_evaluator_sha256,
        "evaluator_labels_read_after_prediction_lock": True,
        "target_adaptation_performed": False,
        "multiplicity_family": "development_selection_descriptive_no_confirmatory_hypothesis",
        "project_exposed_development_only": True,
        "external_or_champion_claim_allowed": False,
        "donor_fingerprint_overlap_status": "unresolved",
        "sealed_outcomes_read": False,
    }
    write_json_exclusive(args.output / "evaluation_receipt.json", receipt)
    freeze_tree(
        args.output,
        {
            "artifact_class": "gse296875_cell_state_transfer_development_evaluation",
            "view_id": VIEW_ID,
            "rows": 7_500,
            "donors": 39,
            "candidate_roster": expected_roster,
            "biological_replication_unit": "donor",
            "project_exposed_development_only": True,
            "external_or_champion_claim_allowed": False,
            "sealed_outcomes_read": False,
            "status": "passed",
        },
    )
    return {**receipt, "comparisons": comparisons}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--prediction-lock", type=Path, required=True)
    parser.add_argument("--expected-prediction-lock-sha256", required=True)
    parser.add_argument("--evaluator", type=Path, required=True)
    parser.add_argument("--expected-evaluator-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), sort_keys=True))
