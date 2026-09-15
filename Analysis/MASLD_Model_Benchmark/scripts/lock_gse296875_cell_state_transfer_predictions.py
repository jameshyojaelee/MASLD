#!/usr/bin/env python3
"""Commit outcome-blind GSE296875 predictions before evaluator labels are opened."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import h5py
import numpy as np

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive
from masld_bench.hashing import canonical_sha256, sha256_file


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


class PredictionLockError(ValueError):
    """Raised when prediction files cannot be committed safely."""


def _decode(values: Any) -> list[str]:
    return [value.decode() if isinstance(value, bytes) else str(value) for value in values]


def _validate_prediction(path: Path, row_ids: list[str]) -> dict[str, Any]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = list(reader.fieldnames or ())
        rows = list(reader)
    expected = ["row_id", "predicted_class", *(f"probability::{label}" for label in ROSTER)]
    if fields != expected or len(rows) != 7_500:
        raise PredictionLockError(f"prediction schema differs: {path}")
    if any("donor" in field.lower() or "label" in field.lower() or "histolog" in field.lower() for field in fields):
        raise PredictionLockError("prediction contains evaluator-only fields")
    if [row["row_id"] for row in rows] != row_ids:
        raise PredictionLockError("prediction row order differs")
    try:
        probabilities = np.asarray(
            [[float(row[f"probability::{label}"]) for label in ROSTER] for row in rows],
            dtype=np.float64,
        )
    except (KeyError, TypeError, ValueError) as error:
        raise PredictionLockError("prediction probabilities differ") from error
    if (
        not np.all(np.isfinite(probabilities))
        or np.any(probabilities < 0.0)
        or not np.allclose(probabilities.sum(axis=1), 1.0, rtol=0.0, atol=1.0e-7)
        or [row["predicted_class"] for row in rows]
        != [ROSTER[index] for index in np.argmax(probabilities, axis=1)]
    ):
        raise PredictionLockError("prediction probability or hard-call contract differs")
    return {"path": path.resolve().as_posix(), "sha256": sha256_file(path), "rows": len(rows)}


def run(args: argparse.Namespace) -> dict[str, Any]:
    if args.output.exists():
        raise PredictionLockError("refusing to overwrite prediction lock")
    features = args.features.resolve(strict=True)
    baseline = args.baseline_branch.resolve(strict=True)
    transcriptformer = args.transcriptformer_branch.resolve(strict=True)
    for root, expected in (
        (features, args.expected_features_sha256),
        (baseline, args.expected_baseline_sha256),
        (transcriptformer, args.expected_transcriptformer_sha256),
    ):
        if sha256_file(root / "ARTIFACTS.json") != expected:
            raise PredictionLockError(f"input ARTIFACTS SHA-256 differs: {root}")
        verify_frozen_tree(root)
    feature_meta = verify_frozen_tree(features)["metadata"]
    baseline_meta = verify_frozen_tree(baseline)["metadata"]
    tf_meta = verify_frozen_tree(transcriptformer)["metadata"]
    if (
        feature_meta.get("artifact_class")
        != "gse296875_rna_cell_state_external_development_features"
        or feature_meta.get("view_id") != VIEW_ID
        or feature_meta.get("labels_present") is not False
        or feature_meta.get("donor_ids_present") is not False
        or baseline_meta.get("artifact_class")
        != "gse296875_cell_state_outcome_blind_prediction_branch"
        or baseline_meta.get("branch_id") != "task_native_outer_ensemble"
        or baseline_meta.get("models") != list(BASELINE_MODELS)
        or baseline_meta.get("source_fit_census_per_model") != 15
        or tf_meta.get("artifact_class")
        != "gse296875_cell_state_outcome_blind_prediction_branch"
        or tf_meta.get("branch_id") != "transcriptformer_outer_ensemble"
        or tf_meta.get("models") != list(TF_MODELS)
        or tf_meta.get("source_fit_census_per_model") != 15
        or tf_meta.get("source_fit_census_total") != 30
    ):
        raise PredictionLockError("prediction branch authority differs")
    for metadata in (baseline_meta, tf_meta):
        if any(
            metadata.get(key) is not False
            for key in (
                "query_labels_read", "query_donor_ids_read", "query_atac_read",
                "target_adaptation_performed", "metrics_calculated", "sealed_outcomes_read",
            )
        ):
            raise PredictionLockError("prediction branch outcome firewall differs")
    with h5py.File(features / "rna_query.h5", "r") as handle:
        row_ids = _decode(handle["obs/row_id"][:])
    if len(row_ids) != 7_500 or len(set(row_ids)) != 7_500:
        raise PredictionLockError("feature row contract differs")

    files: dict[str, dict[str, Any]] = {}
    for model_id in BASELINE_MODELS:
        files[model_id] = _validate_prediction(
            baseline / "predictions" / f"{model_id}.tsv", row_ids
        )
    for model_id in TF_MODELS:
        files[model_id] = _validate_prediction(
            transcriptformer / "predictions" / f"{model_id}.tsv", row_ids
        )
    args.output.mkdir(mode=0o750)
    receipt = {
        "schema_version": "masld-bench-gse296875-cell-state-prediction-lock-v1",
        "status": "locked_before_evaluator_join",
        "view_id": VIEW_ID,
        "candidate_roster": [*BASELINE_MODELS, *TF_MODELS],
        "strongest_baseline_selection_rule": "maximum_evaluator_donor_class_balanced_macro_f1_then_minimum_brier_then_model_id",
        "rows": 7_500,
        "row_ids_sha256": canonical_sha256(row_ids),
        "branches": {
            "task_native": {
                "path": baseline.as_posix(),
                "artifacts_sha256": args.expected_baseline_sha256,
            },
            "transcriptformer": {
                "path": transcriptformer.as_posix(),
                "artifacts_sha256": args.expected_transcriptformer_sha256,
            },
        },
        "prediction_files": files,
        "evaluator_path_resolved": False,
        "evaluator_labels_read": False,
        "evaluator_donor_ids_read": False,
        "target_adaptation_performed": False,
        "metrics_calculated": False,
        "project_exposed_development_only": True,
        "external_or_champion_claim_allowed": False,
        "sealed_outcomes_read": False,
    }
    write_json_exclusive(args.output / "selection_lock.json", receipt)
    freeze_tree(
        args.output,
        {
            "artifact_class": "gse296875_cell_state_transfer_prediction_lock",
            "view_id": VIEW_ID,
            "candidate_roster": [*BASELINE_MODELS, *TF_MODELS],
            "rows": 7_500,
            "evaluator_labels_read": False,
            "evaluator_donor_ids_read": False,
            "target_adaptation_performed": False,
            "metrics_calculated": False,
            "external_or_champion_claim_allowed": False,
            "sealed_outcomes_read": False,
            "status": "locked",
        },
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--expected-features-sha256", required=True)
    parser.add_argument("--baseline-branch", type=Path, required=True)
    parser.add_argument("--expected-baseline-sha256", required=True)
    parser.add_argument("--transcriptformer-branch", type=Path, required=True)
    parser.add_argument("--expected-transcriptformer-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


if __name__ == "__main__":
    print(json.dumps(run(parse_args()), sort_keys=True))
