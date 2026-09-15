#!/usr/bin/env python3
"""Authorize exact five-seed observed-multiome CPU bundle execution."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import freeze_tree, reject_symlink_components, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-observed-multiome-factorized-seed-bundle-admission-v1"
ADMISSION_ID = "gse296875_observed_multiome_factorized_seed_bundle_admission_20260825"
SEEDS = [20260824, 20260825, 20260826, 20260827, 20260828]


class SeedBundleAdmissionError(ValueError):
    """Raised when the exact production execution surface differs."""


def _digest(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SeedBundleAdmissionError("JSON object required")
    return value


def _bound_file(root: Path, record: Mapping[str, Any], label: str) -> Path:
    path = reject_symlink_components(root / str(record.get("path", "")), label=label).resolve(strict=True)
    path.relative_to(root)
    if _digest(path) != record.get("sha256"):
        raise SeedBundleAdmissionError(f"{label} drifted")
    return path


def validate(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    if config.get("schema_version") != SCHEMA or config.get("admission_id") != ADMISSION_ID:
        raise SeedBundleAdmissionError("admission identity differs")
    if config.get("dataset_id") != "gse296875" or config.get("stage") != "development":
        raise SeedBundleAdmissionError("dataset or stage differs")

    parent = config.get("parent_contract")
    if not isinstance(parent, dict):
        raise SeedBundleAdmissionError("parent contract is missing")
    contract = reject_symlink_components(root / str(parent.get("path", "")), label="parent contract").resolve(strict=True)
    contract.relative_to(root)
    verify_frozen_tree(contract)
    if _digest(contract / "ARTIFACTS.json") != parent.get("artifacts_sha256"):
        raise SeedBundleAdmissionError("parent contract drifted")
    contract_receipt = _json(contract / "receipt.json")
    if contract_receipt.get("full_rectangle_execution_authorized") is not True or contract_receipt.get("partial_ranking_authorized") is not False:
        raise SeedBundleAdmissionError("parent contract disposition differs")

    bindings = config.get("source_bindings")
    expected_bindings = {"master_contract_config", "seed_bundle_driver", "seed_bundle_unit_test", "seed_bundle_sbatch", "surface_runner"}
    if not isinstance(bindings, dict) or set(bindings) != expected_bindings:
        raise SeedBundleAdmissionError("source binding roster differs")
    paths = {label: _bound_file(root, record, label) for label, record in bindings.items()}

    execution = config.get("execution")
    expected_commands = [f"sbatch --export=ALL,SEED={seed} slurm/run_gse296875_observed_multiome_factorized_seed_bundle_cpu.sbatch" for seed in SEEDS]
    if execution != {"seeds": SEEDS, "jobs": 5, "surfaces_per_job": 25, "surface_seed_runs": 125, "models_per_surface": 5, "prediction_matrices": 625, "submission_commands": expected_commands}:
        raise SeedBundleAdmissionError("execution census or commands differ")
    resources = config.get("resources")
    expected_resources = {"partition": "cpu", "account": "nslab", "qos": "nslab", "nodes_per_job": 1, "tasks_per_job": 1, "cpus_per_job": 8, "memory_gb_per_job": 64, "wall_hours_per_job": 4, "maximum_concurrent_cpus": 40, "maximum_concurrent_memory_gb": 320, "total_cpu_hour_ceiling": 160, "total_memory_gb_hour_ceiling": 1280, "gpu_hours": 0}
    if resources != expected_resources:
        raise SeedBundleAdmissionError("resource plan differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or firewall != {"prediction_commit_precedes_evaluator": True, "unit_metrics_written_per_surface": True, "evaluator_values_read_by_admission": False, "prediction_values_read_by_admission": False, "development_metrics_calculated_by_admission": False, "partial_ranking_authorized": False, "promotion_authorized": False, "sealed_outcome_read": False}:
        raise SeedBundleAdmissionError("firewall differs")
    failure = config.get("failure_policy")
    if not isinstance(failure, dict) or set(failure.values()) != {True}:
        raise SeedBundleAdmissionError("failure policy differs")

    driver_text = paths["seed_bundle_driver"].read_text(encoding="utf-8")
    for token in ("prediction_commit", "committed_unix_time_ns", '"evaluate"', "unit_metrics_written", "partial_ranking_authorized"):
        if token not in driver_text:
            raise SeedBundleAdmissionError(f"driver firewall token is missing: {token}")
    sbatch_text = paths["seed_bundle_sbatch"].read_text(encoding="utf-8")
    for token in ("#SBATCH --partition=cpu", "#SBATCH --account=nslab", "#SBATCH --qos=nslab", "#SBATCH --cpus-per-task=8", "#SBATCH --mem=64G", "#SBATCH --time=04:00:00", "SEED must be exported", "model-check-321-21099900"):
        if token not in sbatch_text:
            raise SeedBundleAdmissionError(f"sbatch token is missing: {token}")

    return {
        "schema_version": "masld-bench-observed-multiome-factorized-seed-bundle-admission-receipt-v1",
        "admission_id": ADMISSION_ID,
        "dataset_id": "gse296875",
        "seed_bundle_jobs_authorized": 5,
        "surface_seed_runs_authorized": 125,
        "prediction_matrices_authorized": 625,
        "submission_commands": expected_commands,
        "total_cpu_hour_ceiling": 160,
        "total_memory_gb_hour_ceiling": 1280,
        "gpu_hours": 0,
        "source_bindings_verified": True,
        "parent_contract_verified": True,
        "evaluator_values_read": False,
        "prediction_values_read": False,
        "development_metrics_calculated": False,
        "partial_ranking_authorized": False,
        "promotion_authorized": False,
        "sealed_outcomes_read": False,
        "seed_bundle_execution_authorized": True,
        "next_gate": "five_seed_cpu_bundle_execution_then_complete_rectangle_aggregation"
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    output = reject_symlink_components(args.output, label="seed bundle admission output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    receipt = validate(root, _json(config_path))
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "tests/unit/test_observed_multiome_factorized_seed_bundle_admission.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_factorized_seed_bundle_admission", "admission_id": ADMISSION_ID, "seed_bundle_execution_authorized": True, "partial_ranking_authorized": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
