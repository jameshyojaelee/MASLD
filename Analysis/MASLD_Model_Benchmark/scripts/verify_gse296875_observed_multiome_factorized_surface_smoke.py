#!/usr/bin/env python3
"""Independently verify the first biological surface and unranked metrics."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

import h5py
import numpy as np
from scipy import sparse

from masld_bench.artifacts import freeze_tree, reject_symlink_components, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-observed-multiome-factorized-surface-smoke-verification-v1"
CAMPAIGN_ID = "gse296875_observed_multiome_factorized_surface_smoke_verification_20260825"


class SurfaceSmokeVerificationError(ValueError):
    """Raised when process timing, binding, joins, or metrics differ."""


def _digest(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SurfaceSmokeVerificationError("JSON object required")
    return value


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    path = reject_symlink_components(root / str(record.get("path", "")), label=label).resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise SurfaceSmokeVerificationError(f"{label} drifted")
    return path


def validate_config(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    if config.get("schema_version") != SCHEMA or config.get("campaign_id") != CAMPAIGN_ID or config.get("dataset_id") != "gse296875" or config.get("stage") != "smoke":
        raise SurfaceSmokeVerificationError("surface verification identity differs")
    expected = config.get("expected")
    if not isinstance(expected, dict) or expected.get("held_genomic_fold") != 0 or expected.get("held_donor_fold") != 0 or expected.get("seed") != 20260825 or expected.get("biological_donors") != 11 or expected.get("donor_lineage_units") != 55 or expected.get("held_targets") != 1000 or expected.get("absolute_metric_tolerance") != 1e-12 or len(expected.get("models", [])) != 5:
        raise SurfaceSmokeVerificationError("surface verification census differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise SurfaceSmokeVerificationError("surface verification firewall is open")
    parents = config.get("parents")
    if not isinstance(parents, dict) or set(parents) != {"surface_smoke", "evaluator_outcome"}:
        raise SurfaceSmokeVerificationError("surface verification parents differ")
    return {key: _tree(root, record, key) for key, record in parents.items()}


def _decode(values: Any) -> list[str]:
    return [value.decode("utf-8") if isinstance(value, bytes) else str(value) for value in values]


def _read_csr(group: Any) -> sparse.csr_matrix:
    return sparse.csr_matrix((np.asarray(group["data"][:]), np.asarray(group["indices"][:], dtype=np.int64), np.asarray(group["indptr"][:], dtype=np.int64)), shape=tuple(int(value) for value in group["shape"][:]))


def _skills(observed: np.ndarray, predicted: np.ndarray, pseudocount: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    totals = observed.sum(axis=1)
    probabilities = (predicted + pseudocount) / (predicted.sum(axis=1, keepdims=True) + pseudocount * predicted.shape[1])
    model = np.zeros(observed.shape[0], dtype=float)
    null = np.zeros(observed.shape[0], dtype=float)
    for index in range(observed.shape[0]):
        positive = observed[index] > 0
        if totals[index] > 0:
            model[index] = 2 * np.sum(observed[index, positive] * np.log(observed[index, positive] / (totals[index] * probabilities[index, positive])))
            null[index] = 2 * np.sum(observed[index, positive] * np.log(observed[index, positive] / (totals[index] / observed.shape[1])))
    valid = (totals > 0) & (null > 0)
    skill = np.full(observed.shape[0], np.nan)
    skill[valid] = 1 - model[valid] / null[valid]
    return skill, model, null


def run(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    resolved = validate_config(root, config)
    smoke = resolved["surface_smoke"]
    receipt = _json(smoke / "receipt.json")
    if receipt.get("prediction_committed_before_evaluator_launch") is not True or receipt.get("evaluator_artifact_read_before_prediction_commit") is not False or receipt.get("partial_ranking_authorized") is not False:
        raise SurfaceSmokeVerificationError("surface process receipt differs")
    fit_binding = _json(smoke / "bindings/fit.json")
    predict_binding = _json(smoke / "bindings/predict.json")
    evaluate_binding = _json(smoke / "bindings/evaluate.json")
    if "evaluator_outcome" in fit_binding or "held_target_features" in fit_binding or "training_label" in predict_binding or "evaluator_outcome" in predict_binding or "frozen_fit_state" in evaluate_binding or "model_input" in evaluate_binding or "training_label" in evaluate_binding:
        raise SurfaceSmokeVerificationError("surface process firewall differs")
    commit = _json(smoke / "prediction_commit/commit.json")
    evaluator_mtime = (smoke / "evaluator/receipt.json").stat().st_mtime_ns
    if int(commit["committed_unix_time_ns"]) >= evaluator_mtime:
        raise SurfaceSmokeVerificationError("evaluator artifact predates prediction commit")
    predictions = smoke / "predictions"
    verify_frozen_tree(predictions)
    if _digest(predictions / "ARTIFACTS.json") != commit["prediction_bundle_artifacts_sha256"]:
        raise SurfaceSmokeVerificationError("surface prediction commit differs")
    identifiers = _json(predictions / "identifiers.json")
    with h5py.File(resolved["evaluator_outcome"] / "evaluator_outcomes.h5", "r") as evaluator:
        folds = np.asarray(evaluator["rows/donor_fold"][:], dtype=int)
        query = np.flatnonzero(folds == 0)
        rows = _decode(evaluator["rows/row_hash"][:])
        units = _decode(evaluator["rows/unit_hash"][:])
        lineages = _decode(evaluator["rows/lineage"][:])
        targets = _decode(evaluator["target_atac/target_hash"][:])
        counts = _read_csr(evaluator["target_atac/counts_csr"])[query].toarray().astype(float)
    query_rows = [rows[index] for index in query]
    query_units = [units[index] for index in query]
    query_lineages = [lineages[index] for index in query]
    if identifiers != {"row_hash": query_rows, "unit_hash": query_units, "lineage": query_lineages, "target_hash": targets} or len(set(query_units)) != 11 or counts.shape != (55, 1000):
        raise SurfaceSmokeVerificationError("surface evaluator join differs")
    reported = _json(smoke / "evaluator/receipt.json")
    tolerance = float(config["expected"]["absolute_metric_tolerance"])
    maximum_difference = 0.0
    for model_id in config["expected"]["models"]:
        with (predictions / f"{model_id}.npy").open("rb") as handle:
            estimate = np.load(handle, allow_pickle=False)
        skill, model_deviance, null_deviance = _skills(counts, estimate, float(evaluate_binding["evaluation"]["pseudocount"]))
        valid = np.isfinite(skill)
        donor_means = [np.mean([skill[index] for index, unit_value in enumerate(query_units) if unit_value == unit and valid[index]]) for unit in sorted(set(query_units))]
        observed = reported["metrics_unranked"][model_id]
        derived = {"donor_macro_profile_deviance_skill": float(np.mean(donor_means)), "mean_model_deviance": float(np.mean(model_deviance[valid])), "mean_uniform_null_deviance": float(np.mean(null_deviance[valid]))}
        for key, value in derived.items():
            difference = abs(value - float(observed[key]))
            maximum_difference = max(maximum_difference, difference)
            if difference > tolerance:
                raise SurfaceSmokeVerificationError("independent surface metric differs")
    return {"schema_version": "masld-bench-observed-multiome-factorized-surface-smoke-verification-receipt-v1", "campaign_id": CAMPAIGN_ID, "dataset_id": "gse296875", "held_genomic_fold": 0, "held_donor_fold": 0, "seed": 20260825, "models_verified": 5, "biological_donors": 11, "donor_lineage_units": 55, "held_targets": 1000, "prediction_commit_verified": True, "process_firewall_verified": True, "independent_metrics_verified": True, "maximum_absolute_metric_difference": maximum_difference, "fit_state_bound": False, "model_input_bound": False, "training_label_bound": False, "partial_ranking_authorized": False, "promotion_authorized": False, "sealed_outcomes_read": False, "full_rectangle_planning_authorized": True, "full_rectangle_execution_authorized": False, "next_gate": "prespecified_full_rectangle_contract"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    output = reject_symlink_components(args.output, label="surface smoke verification output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    receipt = run(root, _json(config_path))
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "tests/unit/test_verify_observed_multiome_factorized_surface_smoke.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_factorized_surface_smoke_verification", "campaign_id": CAMPAIGN_ID, "full_rectangle_planning_authorized": True, "full_rectangle_execution_authorized": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
