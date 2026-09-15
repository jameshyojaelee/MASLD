#!/usr/bin/env python3
"""Commit label-free GSE281367 ATAC transport predictions before scoring."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-atac-transport-prediction-bundle-v1"
COMMIT_SCHEMA = "masld-bench-atac-transport-prediction-commit-v1"
ROTATIONS = {
    "valid_context_test_target": ("valid", "test"),
    "test_context_valid_target": ("test", "valid"),
}
OUTPUT_FAMILIES = {
    "masked_accessibility_profile",
    "masked_accessibility_count",
    "functional_track_profile",
    "functional_track_count",
}
INPUT_REGIMES = {"sequence_only", "observed_atac"}
LAYOUTS = {"donor_lineage", "lineage_broadcast"}


class ATACPredictionCommitError(ValueError):
    """Raised when a prediction bundle cannot be committed outcome-blind."""


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ATACPredictionCommitError(f"JSON object required: {path}")
    return value


def validate_evidence_contract(manifest: Mapping[str, Any]) -> None:
    evidence = manifest.get("evidence_contract")
    if not isinstance(evidence, dict):
        raise ATACPredictionCommitError("prediction evidence contract is absent")
    for field in (
        "condition_labels_consumed",
        "phenotype_values_consumed",
        "development_outcomes_consumed",
        "sealed_data_consumed",
        "target_role_atac_consumed",
        "target_role_contigs_entered_query_transform",
        "missing_as_zero",
    ):
        if evidence.get(field) is not False:
            raise ATACPredictionCommitError(f"prediction firewall opened: {field}")
    regime = manifest.get("input_regime")
    context_role, target_role = ROTATIONS.get(str(manifest.get("rotation_id")), (None, None))
    if context_role is None:
        raise ATACPredictionCommitError("rotation ID differs")
    if evidence.get("scored_target_role") != target_role:
        raise ATACPredictionCommitError("scored target role differs")
    if regime == "sequence_only":
        if (
            evidence.get("observed_atac_context_role") != "none"
            or evidence.get("query_atac_consumed") is not False
            or manifest.get("prediction_layout") != "lineage_broadcast"
        ):
            raise ATACPredictionCommitError("sequence-only view consumed or encoded query ATAC")
    elif regime == "observed_atac":
        if (
            evidence.get("observed_atac_context_role") != context_role
            or evidence.get("query_atac_consumed") is not True
            or manifest.get("prediction_layout") != "donor_lineage"
        ):
            raise ATACPredictionCommitError("observed-ATAC view lacks rotation-specific query evidence")
    else:
        raise ATACPredictionCommitError("input regime differs")
    receptive_field = evidence.get("receptive_field_bp")
    if not isinstance(receptive_field, int) or not 1 <= receptive_field <= 524288:
        raise ATACPredictionCommitError("receptive field exceeds the exchange contract")
    if evidence.get("query_transform_scope") not in {
        "source_training_only",
        "per_donor_context_only",
        "target_context_transductive_separate_lane",
    }:
        raise ATACPredictionCommitError("query-transform scope differs")


def validate_array_contract(
    values: Any,
    missing: Any,
    *,
    output_family: str,
    layout: str,
    minimum_coverage: float = 0.95,
) -> dict[str, Any]:
    import numpy as np

    score = np.asarray(values)
    states = np.asarray(missing)
    profile = output_family.endswith("profile")
    expected = (12, 4, 16000, 20) if layout == "donor_lineage" and profile else None
    if layout == "donor_lineage" and not profile:
        expected = (12, 4, 16000)
    if layout == "lineage_broadcast" and profile:
        expected = (4, 16000, 20)
    if layout == "lineage_broadcast" and not profile:
        expected = (4, 16000)
    expected_missing = expected[:-1] if profile else expected
    if tuple(score.shape) != expected or tuple(states.shape) != expected_missing:
        raise ATACPredictionCommitError("prediction or missing-state shape differs")
    if states.dtype.kind not in "ui" or np.any(states > 7):
        raise ATACPredictionCommitError("missing-state codes differ")
    observed = states == 0
    expanded = np.broadcast_to(observed[..., None], score.shape) if profile else observed
    if np.any(~np.isfinite(score[expanded])) or np.any(score[expanded] < 0):
        raise ATACPredictionCommitError("observed predictions are not finite and nonnegative")
    if np.any(~np.isnan(score[~expanded])):
        raise ATACPredictionCommitError("missing predictions must be NaN, never zero")
    if layout == "donor_lineage":
        coverage = observed.reshape(48, 16000).mean(axis=1)
    else:
        coverage = observed.reshape(4, 16000).mean(axis=1)
    if np.any(coverage < minimum_coverage):
        raise ATACPredictionCommitError("prediction coverage is below the frozen threshold")
    return {
        "minimum_coverage": float(coverage.min()),
        "maximum_coverage": float(coverage.max()),
        "observed_predictions": int(observed.sum()),
        "missing_predictions": int((~observed).sum()),
    }


def validate_bundle(
    prediction_root: Path,
    prediction_artifacts_sha256: str,
    axis_root: Path,
    axis_artifacts_sha256: str,
) -> dict[str, Any]:
    import numpy as np

    prediction_root = prediction_root.resolve(strict=True)
    axis_root = axis_root.resolve(strict=True)
    prediction_artifacts = verify_frozen_tree(prediction_root)
    axis_artifacts = verify_frozen_tree(axis_root)
    if _digest(prediction_root / "ARTIFACTS.json") != prediction_artifacts_sha256:
        raise ATACPredictionCommitError("prediction tree hash differs")
    if (
        _digest(axis_root / "ARTIFACTS.json") != axis_artifacts_sha256
        or axis_artifacts.get("metadata", {}).get("artifact_class")
        != "gse281367_label_free_atac_exchange_axis"
    ):
        raise ATACPredictionCommitError("exchange-axis authority differs")
    manifest = _load_json(prediction_root / "prediction_bundle.json")
    if (
        manifest.get("schema_version") != SCHEMA
        or manifest.get("input_regime") not in INPUT_REGIMES
        or manifest.get("output_family") not in OUTPUT_FAMILIES
        or manifest.get("prediction_layout") not in LAYOUTS
        or manifest.get("axis_artifacts_sha256") != axis_artifacts_sha256
        or manifest.get("target_dataset_id") != "gse281367"
        or manifest.get("source_dataset_id") != "gse296875"
    ):
        raise ATACPredictionCommitError("prediction manifest identity differs")
    validate_evidence_contract(manifest)
    model_manifest = manifest.get("model_manifest")
    if not isinstance(model_manifest, dict):
        raise ATACPredictionCommitError("model manifest binding is absent")
    model_path = Path(str(model_manifest.get("path", ""))).resolve(strict=True)
    if model_path.is_symlink() or not model_path.is_file() or _digest(model_path) != model_manifest.get("sha256"):
        raise ATACPredictionCommitError("model manifest binding differs")
    arrays = manifest.get("arrays")
    if not isinstance(arrays, dict) or set(arrays) != {"predictions", "missing_state"}:
        raise ATACPredictionCommitError("prediction array roster differs")
    paths: dict[str, Path] = {}
    for name, expected_dtype in (("predictions", "float32"), ("missing_state", "uint8")):
        record = arrays[name]
        path = (prediction_root / str(record.get("path", ""))).resolve(strict=True)
        path.relative_to(prediction_root)
        if path.is_symlink() or _digest(path) != record.get("sha256") or record.get("dtype") != expected_dtype:
            raise ATACPredictionCommitError(f"{name} array binding differs")
        paths[name] = path
    predictions = np.load(paths["predictions"], mmap_mode="r", allow_pickle=False)
    missing = np.load(paths["missing_state"], mmap_mode="r", allow_pickle=False)
    array_receipt = validate_array_contract(
        predictions,
        missing,
        output_family=str(manifest["output_family"]),
        layout=str(manifest["prediction_layout"]),
    )
    return {
        "schema_version": COMMIT_SCHEMA,
        "status": "committed_before_outcome_open",
        "bundle_id": manifest.get("bundle_id"),
        "model_id": manifest.get("model_id"),
        "input_regime": manifest.get("input_regime"),
        "output_family": manifest.get("output_family"),
        "prediction_layout": manifest.get("prediction_layout"),
        "rotation_id": manifest.get("rotation_id"),
        "prediction_artifacts_sha256": prediction_artifacts_sha256,
        "axis_artifacts_sha256": axis_artifacts_sha256,
        "model_manifest_sha256": model_manifest["sha256"],
        "array_validation": array_receipt,
        "prediction_values_read_for_structure_validation": True,
        "condition_labels_read": False,
        "phenotype_values_read": False,
        "development_outcomes_read": False,
        "metrics_calculated": False,
        "sealed_data_read": False,
        "thresholds_or_calibration_changed": False,
    }


def commit(
    prediction_root: Path,
    prediction_artifacts_sha256: str,
    axis_root: Path,
    axis_artifacts_sha256: str,
    output: Path,
) -> dict[str, Any]:
    if output.exists():
        raise ATACPredictionCommitError("refusing to overwrite prediction commit")
    receipt = validate_bundle(
        prediction_root,
        prediction_artifacts_sha256,
        axis_root,
        axis_artifacts_sha256,
    )
    output.mkdir(parents=True, mode=0o750)
    write_json_exclusive(output / "commit.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "gse281367_atac_transport_prediction_commit",
            "bundle_id": receipt["bundle_id"],
            "model_id": receipt["model_id"],
            "prediction_artifacts_sha256": prediction_artifacts_sha256,
            "axis_artifacts_sha256": axis_artifacts_sha256,
            "development_outcomes_read": False,
            "metrics_calculated": False,
            "status": "passed",
        },
    )
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prediction-root", required=True, type=Path)
    parser.add_argument("--prediction-artifacts-sha256", required=True)
    parser.add_argument("--axis-root", required=True, type=Path)
    parser.add_argument("--axis-artifacts-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    result = commit(
        arguments.prediction_root,
        arguments.prediction_artifacts_sha256,
        arguments.axis_root,
        arguments.axis_artifacts_sha256,
        arguments.output,
    )
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
