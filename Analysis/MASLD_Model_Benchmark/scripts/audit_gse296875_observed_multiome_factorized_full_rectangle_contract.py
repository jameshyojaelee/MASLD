#!/usr/bin/env python3
"""Freeze the complete five-seed observed-multiome baseline rectangle."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import freeze_tree, reject_symlink_components, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-observed-multiome-factorized-full-rectangle-contract-v1"
CONTRACT_ID = "gse296875_observed_multiome_factorized_full_rectangle_20260825"


class FullRectangleContractError(ValueError):
    """Raised when the prospective rectangle or read-only parents differ."""


def _digest(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise FullRectangleContractError("JSON object required")
    return value


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    path = reject_symlink_components(root / str(record.get("path", "")), label=label).resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise FullRectangleContractError(f"{label} drifted")
    return path


def validate(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    if config.get("schema_version") != SCHEMA or config.get("contract_id") != CONTRACT_ID or config.get("dataset_id") != "gse296875" or config.get("stage") != "development":
        raise FullRectangleContractError("full rectangle identity differs")
    rectangle = config.get("rectangle")
    if rectangle != {"genomic_folds": [0, 1, 2, 3, 4], "donor_folds": [0, 1, 2, 3, 4], "seeds": [20260824, 20260825, 20260826, 20260827, 20260828], "surfaces_per_seed": 25, "surface_seed_runs": 125, "models_per_run": 5, "prediction_matrices": 625, "partial_rectangle_usable": False}:
        raise FullRectangleContractError("full rectangle census differs")
    if config.get("bundles") != {"strategy": "one_non_array_cpu_bundle_per_seed_with_25_sequential_surfaces", "jobs": 5, "partition": "cpu", "account": "nslab", "qos": "nslab", "cpus_per_job": 8, "memory_gb_per_job": 64, "wall_hours_per_job": 4, "maximum_concurrent_cpus": 40, "maximum_concurrent_memory_gb": 320, "total_cpu_hour_ceiling": 160, "total_memory_gb_hour_ceiling": 1280, "gpu_hours": 0}:
        raise FullRectangleContractError("full rectangle resource plan differs")
    preprocessing = config.get("preprocessing")
    if not isinstance(preprocessing, dict) or preprocessing.get("components") != 8 or preprocessing.get("fit_on_outer_training_rows_only") is not True or preprocessing.get("fit_on_training_targets_only") is not True:
        raise FullRectangleContractError("full rectangle preprocessing differs")
    evaluation = config.get("evaluation")
    if not isinstance(evaluation, dict) or evaluation.get("unit_metrics_required") is not True or evaluation.get("bootstrap_replicates_after_complete_rectangle") != 10000 or evaluation.get("ranking_before_complete_rectangle") is not False:
        raise FullRectangleContractError("full rectangle evaluation differs")
    failure = config.get("failure_policy")
    if not isinstance(failure, dict) or set(failure.values()) != {True}:
        raise FullRectangleContractError("full rectangle failure policy differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise FullRectangleContractError("full rectangle firewall is open")
    sources = config.get("source_bindings")
    if not isinstance(sources, dict) or set(sources) != {"surface_runner", "surface_preprocessing", "factorized_adapter"}:
        raise FullRectangleContractError("full rectangle source roster differs")
    for label, record in sources.items():
        path = (root / str(record.get("path", ""))).resolve(strict=True)
        path.relative_to(root)
        if _digest(path) != record.get("sha256"):
            raise FullRectangleContractError(f"{label} source drifted")
    parents = config.get("parents")
    expected = {"surface_smoke_verification", "model_inputs", "training_labels", "target_features", "evaluator_outcomes", "verified_model_inputs", "verified_training_labels", "verified_target_features", "verified_evaluator_outcomes"}
    if not isinstance(parents, dict) or set(parents) != expected:
        raise FullRectangleContractError("full rectangle parent roster differs")
    resolved = {}
    for key, record in parents.items():
        if key == "evaluator_outcomes":
            path = reject_symlink_components(root / str(record.get("path", "")), label=key).resolve(strict=True)
            path.relative_to(root)
            resolved[key] = path
        else:
            resolved[key] = _tree(root, record, key)
    if _json(resolved["surface_smoke_verification"] / "receipt.json").get("full_rectangle_planning_authorized") is not True:
        raise FullRectangleContractError("surface smoke did not authorize rectangle planning")
    for key in ("verified_model_inputs", "verified_training_labels", "verified_target_features", "verified_evaluator_outcomes"):
        if _json(resolved[key] / "receipt.json").get("promotion_gate_passed") is not True:
            raise FullRectangleContractError(f"{key} is not promoted")
    model_receipt = _json(resolved["model_inputs"] / "receipt.json")
    label_receipt = _json(resolved["training_labels"] / "receipt.json")
    target_receipt = _json(resolved["target_features"] / "receipt.json")
    evaluator_receipt = _json(resolved["evaluator_outcomes"] / "receipt.json")
    if model_receipt.get("child_artifact_count") != 5 or label_receipt.get("child_artifact_count") != 25 or target_receipt.get("child_count") != 10 or evaluator_receipt.get("child_artifact_count") != 5:
        raise FullRectangleContractError("full rectangle input child census differs")
    labels = {(int(row["held_genomic_fold"]), int(row["held_donor_fold"])) for row in label_receipt["child_artifacts"]}
    targets = {(row["role"], int(row["held_genomic_fold"])) for row in target_receipt["children"]}
    if labels != {(genomic, donor) for genomic in range(5) for donor in range(5)} or targets != {(role, genomic) for role in ("training", "held") for genomic in range(5)}:
        raise FullRectangleContractError("full rectangle input coverage differs")
    return {"schema_version": "masld-bench-observed-multiome-factorized-full-rectangle-contract-receipt-v1", "contract_id": CONTRACT_ID, "dataset_id": "gse296875", "genomic_folds": 5, "donor_folds": 5, "seeds": 5, "surface_seed_runs": 125, "prediction_matrices_planned": 625, "cpu_bundle_jobs": 5, "total_cpu_hour_ceiling": 160, "total_memory_gb_hour_ceiling": 1280, "gpu_hours": 0, "biological_matrix_values_deserialized": False, "evaluator_artifact_bytes_read": False, "prediction_bundle_read": False, "development_metric_calculated": False, "partial_ranking_authorized": False, "promotion_authorized": False, "sealed_outcomes_read": False, "full_rectangle_execution_authorized": True, "next_gate": "five_seed_cpu_bundle_execution"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    output = reject_symlink_components(args.output, label="full rectangle contract output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    receipt = validate(root, _json(config_path))
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "tests/unit/test_observed_multiome_factorized_full_rectangle_contract.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_factorized_full_rectangle_contract", "contract_id": CONTRACT_ID, "full_rectangle_execution_authorized": True, "partial_ranking_authorized": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
