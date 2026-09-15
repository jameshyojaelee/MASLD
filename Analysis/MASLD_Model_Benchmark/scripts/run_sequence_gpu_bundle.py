#!/usr/bin/env python3
"""Execute one frozen sequence-model bundle with task-level checkpoints."""

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
import sys
from typing import Any

from masld_bench.artifacts import verify_frozen_tree


ALLOWED_EXPORTS = {
    "OUTER_FOLD",
    "LINEAGE",
    "SEED",
    "INPUT_JOB_ID",
    "INPUT_ARTIFACTS_SHA256",
    "MODEL_JOB_ID",
    "MODEL_ARTIFACTS_SHA256",
    "RUN_ID",
}
PASS_STATES = {"passed", "recovered_passed"}
TERMINAL_STATES = PASS_STATES | {
    "blocked_prerequisite",
    "existing_failed",
    "failed",
    "failed_missing_artifact",
    "failed_invalid_artifact",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--bundle-id", required=True)
    parser.add_argument("--requested-time", required=True)
    parser.add_argument("--job-id", required=True)
    parser.add_argument("--stage", type=Path, required=True)
    parser.add_argument("--preflight-only", action="store_true")
    return parser.parse_args()


def utc_now() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_json(path: Path, payload: dict[str, Any]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    temporary.replace(path)


def file_sha256(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def parse_exports(contract: str, job_id: str) -> dict[str, str]:
    exports: dict[str, str] = {}
    for binding in contract.split(";"):
        key, separator, value = binding.partition("=")
        if separator != "=" or key not in ALLOWED_EXPORTS or not re.fullmatch(r"[A-Z0-9_]+", key):
            raise ValueError(f"unsafe or unknown export binding: {binding}")
        if value == "FROM_PREREQUISITE":
            value = job_id
        if not re.fullmatch(r"[A-Za-z0-9_.:-]+", value):
            raise ValueError(f"unsafe export value for {key}")
        exports[key] = value
    return exports


def input_root(root: Path, task: dict[str, str]) -> Path:
    if task["model_id"] == "scbasset":
        name = f"scbasset-{task['split_id']}-inputs-{task['input_job_id']}"
    else:
        name = (
            f"chrombpnet-{task['lineage_id']}-{task['split_id']}"
            f"-inputs-{task['input_job_id']}"
        )
        if (
            task["lineage_id"] == "hepatocyte"
            and task["outer_fold"] == "0"
            and task["input_job_id"] == "21064244"
        ):
            name = "chrombpnet-hepatocyte-donor0-genomic0-inputs-21064244"
    return root / "executions" / name


def verify_task_input(root: Path, task: dict[str, str]) -> str:
    artifact = input_root(root, task)
    verify_frozen_tree(artifact)
    actual = file_sha256(artifact / "ARTIFACTS.json")
    expected = task["input_artifacts_sha256"]
    if expected != "AUTO" and actual != expected:
        raise ValueError(f"input artifact hash differs for {task['logical_work_id']}")
    return actual


def expected_artifact(root: Path, task: dict[str, str], job_id: str) -> Path:
    relative = task["expected_artifact_template"].replace("{SLURM_JOB_ID}", job_id)
    path = (root / relative).resolve()
    execution_root = (root / "executions").resolve()
    if execution_root not in path.parents or "{SLURM_JOB_ID}" in relative:
        raise ValueError("unsafe expected artifact path")
    return path


def load_receipt(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("status") not in TERMINAL_STATES:
        raise ValueError(f"nonterminal checkpoint state in {path}")
    return payload


def recovered_receipt(task: dict[str, str], artifact: Path) -> dict[str, Any]:
    verify_frozen_tree(artifact)
    now = utc_now()
    return {
        "logical_work_id": task["logical_work_id"],
        "status": "recovered_passed",
        "return_code": 0,
        "expected_artifact": str(artifact),
        "artifacts_sha256": file_sha256(artifact / "ARTIFACTS.json"),
        "started_at": now,
        "ended_at": now,
        "recovered_from_existing_frozen_artifact": True,
    }


def main() -> None:
    args = parse_args()
    root = args.root.resolve(strict=True)
    plan_root = args.plan_root.resolve(strict=True)
    stage = args.stage.resolve(strict=True)
    if not re.fullmatch(r"[0-9]+", args.job_id):
        raise SystemExit("invalid Slurm job ID")
    allocations = read_tsv(plan_root / "plan/bundle_allocations.tsv")
    matches = [row for row in allocations if row["bundle_id"] == args.bundle_id]
    if len(matches) != 1:
        raise SystemExit("bundle allocation identity differs")
    allocation = matches[0]
    if (
        allocation["requested_time"] != args.requested_time
        or allocation["qos"] != "nslab"
        or allocation["partition"] != "gpu"
        or allocation["gres"] != "gpu:l40s:1"
        or allocation["cpus"] != "8"
        or allocation["memory"] != "32G"
        or allocation["submission_authority"] != "dispatcher_only_no_manual_sbatch"
    ):
        raise SystemExit("bundle resource or authority contract differs")
    all_tasks = read_tsv(plan_root / "plan/bundle_tasks.tsv")
    tasks = sorted(
        (row for row in all_tasks if row["bundle_id"] == args.bundle_id),
        key=lambda row: int(row["task_order"]),
    )
    if len(tasks) != int(allocation["logical_tasks"]) or not 4 <= len(tasks) <= 8:
        raise SystemExit("bundle logical-task cardinality differs")
    if [int(task["task_order"]) for task in tasks] != list(range(1, len(tasks) + 1)):
        raise SystemExit("bundle task order differs")
    logical_ids = [task["logical_work_id"] for task in tasks]
    if sha256("\n".join(logical_ids).encode("utf-8")).hexdigest() != allocation["logical_work_ids_sha256"]:
        raise SystemExit("bundle logical-task hash differs")
    if any(
        task["submission_authority"] != "dispatcher_only_no_manual_sbatch"
        or task["failure_scope"] != "isolated_continue"
        or "sealed" in task["export_contract"].lower()
        or "outcome" in task["export_contract"].lower()
        for task in tasks
    ):
        raise SystemExit("bundle task safety contract differs")

    resolved_input_hashes = {}
    for task in tasks:
        script = (root / task["script"]).resolve(strict=True)
        if root not in script.parents or script.is_symlink():
            raise SystemExit("unsafe task wrapper path")
        parse_exports(task["export_contract"], args.job_id)
        resolved_input_hashes[task["logical_work_id"]] = verify_task_input(root, task)
        expected_artifact(root, task, args.job_id)
    atomic_json(stage / "preflight.json", {
        "status": "pass",
        "bundle_id": args.bundle_id,
        "logical_tasks": len(tasks),
        "resolved_input_artifacts_sha256": resolved_input_hashes,
        "outcome_paths_available_to_models": False,
        "sealed_paths_available_to_models": False,
    })
    if args.preflight_only:
        return

    receipt_dir = stage / "task_receipts"
    log_dir = stage / "logs"
    receipts: dict[str, dict[str, Any]] = {}
    for task in tasks:
        logical_id = task["logical_work_id"]
        receipt_path = receipt_dir / f"{logical_id}.json"
        existing_receipt = load_receipt(receipt_path)
        if existing_receipt is not None:
            receipts[logical_id] = existing_receipt
            continue
        prerequisite = task["prerequisite_logical_id"]
        if prerequisite != "none":
            prerequisite_receipt = receipts.get(prerequisite) or load_receipt(
                receipt_dir / f"{prerequisite}.json"
            )
            if prerequisite_receipt is None or prerequisite_receipt["status"] not in PASS_STATES:
                now = utc_now()
                receipt = {
                    "logical_work_id": logical_id,
                    "status": "blocked_prerequisite",
                    "return_code": 98,
                    "prerequisite_logical_work_id": prerequisite,
                    "expected_artifact": str(expected_artifact(root, task, args.job_id)),
                    "artifacts_sha256": "not_available",
                    "started_at": now,
                    "ended_at": now,
                }
                atomic_json(receipt_path, receipt)
                receipts[logical_id] = receipt
                continue
        artifact = expected_artifact(root, task, args.job_id)
        if (artifact / "ARTIFACTS.json").is_file():
            receipt = recovered_receipt(task, artifact)
            atomic_json(receipt_path, receipt)
            receipts[logical_id] = receipt
            continue
        if Path(f"{artifact}.failed").is_dir():
            now = utc_now()
            receipt = {
                "logical_work_id": logical_id,
                "status": "existing_failed",
                "return_code": 1,
                "expected_artifact": str(artifact),
                "artifacts_sha256": "not_available",
                "started_at": now,
                "ended_at": now,
            }
            atomic_json(receipt_path, receipt)
            receipts[logical_id] = receipt
            continue

        environment = os.environ.copy()
        environment.update(parse_exports(task["export_contract"], args.job_id))
        environment["INPUT_ARTIFACTS_SHA256"] = resolved_input_hashes[logical_id]
        started_at = utc_now()
        with (log_dir / f"{logical_id}.stdout.txt").open("x", encoding="utf-8") as stdout, (
            log_dir / f"{logical_id}.stderr.txt"
        ).open("x", encoding="utf-8") as stderr:
            result = subprocess.run(
                ["bash", str(root / task["script"])],
                cwd=root,
                env=environment,
                stdout=stdout,
                stderr=stderr,
                check=False,
            )
        status = "failed"
        artifact_hash = "not_available"
        return_code = result.returncode
        if result.returncode == 0 and (artifact / "ARTIFACTS.json").is_file():
            try:
                verify_frozen_tree(artifact)
            except Exception:
                status = "failed_invalid_artifact"
                return_code = 96
            else:
                status = "passed"
                artifact_hash = file_sha256(artifact / "ARTIFACTS.json")
        elif result.returncode == 0:
            status = "failed_missing_artifact"
            return_code = 97
        receipt = {
            "logical_work_id": logical_id,
            "model_id": task["model_id"],
            "stage": task["stage"],
            "lineage_id": task["lineage_id"],
            "outer_fold": int(task["outer_fold"]),
            "seed": int(task["seed"]),
            "status": status,
            "return_code": return_code,
            "expected_artifact": str(artifact),
            "artifacts_sha256": artifact_hash,
            "started_at": started_at,
            "ended_at": utc_now(),
            "recovered_from_existing_frozen_artifact": False,
        }
        atomic_json(receipt_path, receipt)
        receipts[logical_id] = receipt
        print(f"{logical_id}: {status}", flush=True)

    status_counts: dict[str, int] = {}
    for receipt in receipts.values():
        status = str(receipt["status"])
        status_counts[status] = status_counts.get(status, 0) + 1
    terminal_status = "passed" if set(status_counts) <= PASS_STATES else "completed_with_task_failures"
    atomic_json(stage / "bundle_summary.json", {
        "schema_version": "masld-bench-sequence-gpu-bundle-receipt-v2",
        "bundle_id": args.bundle_id,
        "slurm_job_id": int(args.job_id),
        "logical_tasks": len(tasks),
        "status_counts": dict(sorted(status_counts.items())),
        "terminal_status": terminal_status,
        "continue_after_isolated_failure": True,
        "checkpoint_granularity": "logical_task",
        "requested_time": args.requested_time,
        "qos": "nslab",
    })


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"bundle controller failed: {error}", file=sys.stderr)
        raise
