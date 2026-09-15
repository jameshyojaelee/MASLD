#!/usr/bin/env python3
"""Run one frozen outcome-blind task-native sequence bundle sequentially."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import subprocess
from typing import Any, Mapping

from masld_bench.artifacts import verify_frozen_tree


MODELS = {
    "bpnet",
    "chrombpnet",
    "sequence_cnn_control",
    "sequence_transformer_control",
}
REQUIRED_OUTPUTS = {
    "bpnet": {
        "model/bpnet.inference.h5",
        "predictions/profile_probabilities.h5",
        "predictions/regional_counts.tsv",
    },
    "chrombpnet": {
        "model/chrombpnet.h5",
        "model/chrombpnet_nobias.h5",
        "predictions/full_model/profile_probabilities.h5",
        "predictions/full_model/regional_counts.tsv",
        "predictions/nobias_model/profile_probabilities.h5",
        "predictions/nobias_model/regional_counts.tsv",
    },
    "sequence_cnn_control": {
        "model/sequence_cnn_control.h5",
        "predictions/fixed_ccre/profile_probabilities.h5",
        "predictions/fixed_ccre/regional_counts.tsv",
        "predictions/fixed_ccre/standardized_predictions.tsv",
    },
    "sequence_transformer_control": {
        "model/sequence_transformer_control.h5",
        "predictions/fixed_ccre/profile_probabilities.h5",
        "predictions/fixed_ccre/regional_counts.tsv",
        "predictions/fixed_ccre/standardized_predictions.tsv",
    },
}
FALSE_FLAGS = {
    "bpnet": ("benchmark_metrics_calculated", "evaluator_outcomes_exposed"),
    "chrombpnet": ("benchmark_metrics_calculated", "test_outcomes_used"),
    "sequence_cnn_control": ("held_donor_atac_exposed", "test_outcomes_used"),
    "sequence_transformer_control": ("held_donor_atac_exposed", "test_outcomes_used"),
}
EXPORT_KEYS = {"OUTER_FOLD", "LINEAGE", "SEED", "INPUT_JOB_ID", "INPUT_ARTIFACTS_SHA256"}
PASS_STATES = {"passed", "recovered_passed"}


class BundleError(RuntimeError):
    """Raised when bundle execution would violate the frozen plan."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--bundle-id", required=True)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--stage", type=Path, required=True)
    return parser.parse_args()


def now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def inside(root: Path, value: str, label: str, *, strict: bool = True) -> Path:
    path = Path(value).resolve(strict=strict)
    if root not in path.parents:
        raise BundleError(f"{label} escapes benchmark root")
    return path


def parse_exports(value: str) -> dict[str, str]:
    exports: dict[str, str] = {}
    for binding in value.split(";"):
        key, separator, item = binding.partition("=")
        if separator != "=" or key not in EXPORT_KEYS or not re.fullmatch(r"[A-Za-z0-9_.:-]+", item):
            raise BundleError("task export contract differs")
        exports[key] = item
    if set(exports) != EXPORT_KEYS:
        raise BundleError("task export key roster differs")
    return exports


def expected_artifact(root: Path, task: Mapping[str, str], job_id: str) -> Path:
    relative = task["expected_artifact_template"].replace("{SLURM_JOB_ID}", job_id)
    path = (root / relative).resolve(strict=False)
    if root / "executions" not in path.parents or "{SLURM_JOB_ID}" in relative:
        raise BundleError("expected artifact path differs")
    return path


def validate_fit(path: Path, task: Mapping[str, str]) -> str:
    verify_frozen_tree(path)
    artifact_path = path / "ARTIFACTS.json"
    payload = json.loads(artifact_path.read_text(encoding="utf-8"))
    metadata = payload.get("metadata", {})
    model = task["model_id"]
    if (
        model not in MODELS
        or metadata.get("model_id") != model
        or metadata.get("dataset_id") != "gse296875"
        or metadata.get("lineage_id") != task["lineage_id"]
        or metadata.get("split_id") != task["split_id"]
        or metadata.get("seed") != int(task["seed"])
        or metadata.get("status") != "passed"
        or any(metadata.get(flag) is not False for flag in FALSE_FLAGS[model])
    ):
        raise BundleError("fit metadata or outcome firewall differs")
    paths = {str(row.get("path")) for row in payload.get("artifacts", ())}
    if not REQUIRED_OUTPUTS[model].issubset(paths):
        raise BundleError("fit output roster differs")
    if any("outcome" in item.lower() or "metric" in item.lower() for item in paths):
        raise BundleError("fit artifact contains outcome or metric path")
    return digest(artifact_path)


def main() -> None:
    args = parse_args()
    root = args.root.resolve(strict=True)
    plan_root = args.plan_root.resolve(strict=True)
    stage = args.stage.resolve(strict=True)
    if not re.fullmatch(r"[0-9]+", args.job_id):
        raise BundleError("Slurm job ID differs")
    verify_frozen_tree(plan_root)
    summary = json.loads((plan_root / "plan/summary.json").read_text(encoding="utf-8"))
    if (
        summary.get("status") != "pass"
        or summary.get("expected_fits") != 500
        or summary.get("expected_prediction_views") != 625
        or summary.get("prediction_or_metric_values_read") is not False
        or summary.get("outcome_paths_available") is not False
        or summary.get("manual_gpu_submission_allowed") is not False
    ):
        raise BundleError("frozen plan safety summary differs")
    allocations = read_tsv(plan_root / "plan/bundle_allocations.tsv")
    matches = [row for row in allocations if row["bundle_id"] == args.bundle_id]
    if len(matches) != 1:
        raise BundleError("bundle identity differs")
    allocation = matches[0]
    if (
        allocation["partition"] != "gpu"
        or allocation["account"] != "nslab"
        or allocation["qos"] != "nslab"
        or allocation["gres"] != "gpu:l40s:1"
        or allocation["cpus"] != "8"
        or allocation["memory"] != "32G"
        or allocation["time_limit"] != "48:00:00"
        or allocation["submission_authority"] != "dispatcher_only_no_manual_gpu_sbatch"
    ):
        raise BundleError("bundle resources or authority differ")
    tasks = sorted(
        (row for row in read_tsv(plan_root / "plan/bundle_tasks.tsv") if row["bundle_id"] == args.bundle_id),
        key=lambda row: int(row["task_order"]),
    )
    if len(tasks) != int(allocation["logical_fits"]) or not 8 <= len(tasks) <= 16:
        raise BundleError("bundle fit cardinality differs")
    logical_ids = [row["logical_work_id"] for row in tasks]
    if sha256("\n".join(logical_ids).encode()).hexdigest() != allocation["logical_work_ids_sha256"]:
        raise BundleError("bundle logical-work hash differs")
    if [int(row["task_order"]) for row in tasks] != list(range(1, len(tasks) + 1)):
        raise BundleError("bundle task order differs")

    for task in tasks:
        if (
            task["outcome_paths_available_to_fit"] != "false"
            or task["sealed_paths_available_to_fit"] != "false"
            or task["seed"] not in {"20260824", "20260825", "20260826", "20260827", "20260828"}
        ):
            raise BundleError("task outcome, seal, or seed contract differs")
        wrapper = (root / task["task_wrapper"]).resolve(strict=True)
        if root not in wrapper.parents or wrapper.is_symlink() or digest(wrapper) != task["task_wrapper_sha256"]:
            raise BundleError("task wrapper source differs")
        input_path = inside(root, task["input_path"], "task input")
        verify_frozen_tree(input_path)
        if digest(input_path / "ARTIFACTS.json") != task["input_artifacts_sha256"]:
            raise BundleError("task input artifact hash differs")
        parse_exports(task["export_contract"])
        expected_artifact(root, task, args.job_id)
    atomic_json(stage / "preflight.json", {
        "status": "pass",
        "bundle_id": args.bundle_id,
        "logical_fits": len(tasks),
        "outcome_paths_available_to_fits": False,
        "sealed_paths_available_to_fits": False,
        "prediction_or_metric_values_read": False,
    })

    receipts: dict[str, dict[str, Any]] = {}
    for task in tasks:
        logical_id = task["logical_work_id"]
        receipt_path = stage / "task_receipts" / f"{logical_id}.json"
        artifact = expected_artifact(root, task, args.job_id)
        started = now()
        status = "failed"
        return_code = 1
        artifact_hash = "not_available"
        recovered = False
        if (artifact / "ARTIFACTS.json").is_file():
            artifact_hash = validate_fit(artifact, task)
            status = "recovered_passed"
            return_code = 0
            recovered = True
        elif Path(f"{artifact}.failed").exists():
            status = "existing_failed"
        else:
            environment = os.environ.copy()
            environment.update(parse_exports(task["export_contract"]))
            wrapper = root / task["task_wrapper"]
            with (stage / "logs" / f"{logical_id}.stdout.txt").open("x", encoding="utf-8") as stdout, (
                stage / "logs" / f"{logical_id}.stderr.txt"
            ).open("x", encoding="utf-8") as stderr:
                result = subprocess.run(
                    ["bash", str(wrapper)], cwd=root, env=environment,
                    stdout=stdout, stderr=stderr, check=False,
                )
            return_code = result.returncode
            if result.returncode == 0 and (artifact / "ARTIFACTS.json").is_file():
                try:
                    artifact_hash = validate_fit(artifact, task)
                except Exception:
                    status = "failed_invalid_artifact"
                    return_code = 96
                else:
                    status = "passed"
            elif result.returncode == 0:
                status = "failed_missing_artifact"
                return_code = 97
        receipt = {
            "logical_work_id": logical_id,
            "model_id": task["model_id"],
            "lineage_id": task["lineage_id"],
            "outer_fold": int(task["outer_fold"]),
            "seed": int(task["seed"]),
            "status": status,
            "return_code": return_code,
            "expected_artifact": str(artifact),
            "artifacts_sha256": artifact_hash,
            "started_at": started,
            "ended_at": now(),
            "recovered_from_existing_frozen_artifact": recovered,
        }
        atomic_json(receipt_path, receipt)
        receipts[logical_id] = receipt
        print(f"{logical_id}: {status}", flush=True)

    counts = Counter(row["status"] for row in receipts.values())
    terminal_status = "passed" if set(counts) <= PASS_STATES else "completed_with_fit_failures"
    atomic_json(stage / "bundle_summary.json", {
        "schema_version": "masld-bench-sequence-task-native-five-seed-bundle-receipt-v1",
        "bundle_id": args.bundle_id,
        "slurm_job_id": int(args.job_id),
        "logical_fits": len(tasks),
        "status_counts": dict(sorted(counts.items())),
        "terminal_status": terminal_status,
        "checkpoint_granularity": "logical_fit",
        "continue_after_isolated_fit_failure": True,
        "qos": "nslab",
        "outcome_paths_available_to_fits": False,
        "sealed_paths_available_to_fits": False,
        "prediction_or_metric_values_read": False,
    })


if __name__ == "__main__":
    main()
