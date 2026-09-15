#!/usr/bin/env python3
"""Submit one frozen ready-queue snapshot through four dependency chains."""

from __future__ import annotations

import argparse
from collections.abc import Mapping, Sequence
import fcntl
from hashlib import sha256
import json
from pathlib import Path
import re
import subprocess
import time
from typing import Any

from masld_bench.artifacts import freeze_tree, reject_symlink_components, write_json_exclusive
from scripts import gpu_bundle_dispatcher as dispatcher


SCHEMA = "masld-bench-gpu-queue-prime-all-ready-v1"
CAMPAIGN_ID = "gpu_queue_prime_all_ready_20260825"


class QueuePrimeError(RuntimeError):
    """Raised when the frozen queue snapshot cannot be submitted safely."""


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        while block := handle.read(1024 * 1024):
            value.update(block)
    return value.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise QueuePrimeError("JSON object required")
    return value


def _header(text: str, key: str) -> str:
    match = re.search(rf"^#SBATCH --{re.escape(key)}(?:=|\s+)([^\s]+)", text, re.MULTILINE)
    if match is None:
        raise QueuePrimeError(f"missing sbatch header: {key}")
    return match.group(1)


def _hours(value: str) -> float:
    base = 0.0
    if "-" in value:
        days, value = value.split("-", 1)
        base = 24.0 * int(days)
    fields = [int(field) for field in value.split(":")]
    if len(fields) != 3:
        raise QueuePrimeError("wall time must use HH:MM:SS or D-HH:MM:SS")
    return base + fields[0] + fields[1] / 60.0 + fields[2] / 3600.0


def _memory_gb(value: str) -> float:
    if value.endswith("G"):
        return float(value[:-1])
    if value.endswith("M"):
        return float(value[:-1]) / 1024.0
    raise QueuePrimeError("memory must use G or M units")


def load_selected(root: Path, queue: Path, state: Path, config: Mapping[str, Any]) -> list[dict[str, Any]]:
    if config.get("schema_version") != SCHEMA or config.get("campaign_id") != CAMPAIGN_ID:
        raise QueuePrimeError("queue-prime identity differs")
    if config.get("account") != "nslab" or config.get("qos") != "nslab" or config.get("partition") != "gpu" or config.get("gres") != "gpu:l40s:1":
        raise QueuePrimeError("GPU routing differs")
    if config.get("maximum_running_campaign_jobs") != 4 or config.get("pending_job_limit") is not None or config.get("dependency_type") != "afterany":
        raise QueuePrimeError("queue concurrency policy differs")
    bindings = config.get("source_bindings")
    if not isinstance(bindings, dict) or set(bindings) != {"submitter", "unit_test", "sbatch"}:
        raise QueuePrimeError("queue-prime source roster differs")
    for label, record in bindings.items():
        if not isinstance(record, dict):
            raise QueuePrimeError(f"invalid source binding: {label}")
        path = (root / str(record.get("path", ""))).resolve(strict=True)
        path.relative_to(root)
        if _digest(path) != record.get("sha256"):
            raise QueuePrimeError(f"queue-prime source drifted: {label}")
    bundle_ids = config.get("bundle_ids")
    if not isinstance(bundle_ids, list) or len(bundle_ids) != len(set(bundle_ids)) or any(not isinstance(value, str) for value in bundle_ids):
        raise QueuePrimeError("bundle roster differs")
    by_id: dict[str, dict[str, Any]] = {}
    for path in sorted(queue.glob("*.json")):
        item = dispatcher.load_item(path, root)
        if item["bundle_id"] in bundle_ids:
            if item["bundle_id"] in by_id:
                raise QueuePrimeError("bundle identity is duplicated")
            item["_queue_path"] = path
            by_id[item["bundle_id"]] = item
    if set(by_id) != set(bundle_ids):
        raise QueuePrimeError("frozen bundle is absent from the live queue")
    claims = state / "claims"
    selected = []
    for bundle_id in bundle_ids:
        item = by_id[bundle_id]
        if item["enabled"] is not True or item["_ready"] is not True or (claims / f"{bundle_id}.json").exists():
            raise QueuePrimeError(f"bundle is not ready and unclaimed: {bundle_id}")
        text = item["_wrapper"].read_text(encoding="utf-8")
        if "#SBATCH --array" in text or "innovation" in text.lower() or "masld" in _header(text, "job-name").lower():
            raise QueuePrimeError(f"wrapper routing or generic name differs: {bundle_id}")
        item["_job_name"] = _header(text, "job-name")
        item["_cpus"] = int(_header(text, "cpus-per-task"))
        item["_memory_gb"] = _memory_gb(_header(text, "mem"))
        item["_hours"] = _hours(_header(text, "time"))
        selected.append(item)
    expected = config.get("expected")
    observed = {
        "new_jobs": len(selected),
        "total_visible_gpu_jobs_after_submission": len(selected) + 1,
        "logical_tasks": sum(int(item["logical_tasks"]) for item in selected),
        "requested_gpu_hours_ceiling": round(sum(float(item["_hours"]) for item in selected)),
        "requested_cpu_hours_ceiling": round(sum(float(item["_hours"]) * int(item["_cpus"]) for item in selected)),
        "requested_memory_gb_hours_ceiling": round(sum(float(item["_hours"]) * float(item["_memory_gb"]) for item in selected)),
        "generic_job_name_violations": 0,
        "innovation_jobs": 0,
        "array_jobs": 0,
    }
    if observed != expected:
        raise QueuePrimeError(f"resource census differs: {observed}")
    return selected


def plan_chains(items: Sequence[Mapping[str, Any]], root_job_id: str, *, lanes: int = 4, root_hours: float = 120.0) -> list[dict[str, Any]]:
    if lanes != 4 or not root_job_id.isdigit() or not items:
        raise QueuePrimeError("dependency-chain inputs differ")
    parents: list[str | None] = [root_job_id, None, None, None]
    loads = [float(root_hours), 0.0, 0.0, 0.0]
    plan: list[dict[str, Any]] = []
    for index, item in enumerate(items):
        if index < lanes - 1:
            lane = index + 1
        else:
            lane = min(range(lanes), key=lambda value: (loads[value], value))
        dependency = parents[lane]
        token = f"planned:{item['bundle_id']}"
        plan.append({"bundle_id": item["bundle_id"], "lane": lane, "dependency_job_id": dependency})
        parents[lane] = token
        loads[lane] += float(item["_hours"])
    return plan


def _gpu_snapshot(user: str) -> list[dict[str, str]]:
    result = subprocess.run(["squeue", "-u", user, "-p", "gpu", "-h", "-o", "%i|%T"], check=True, text=True, capture_output=True)
    rows = []
    for line in result.stdout.splitlines():
        if line.strip():
            job_id, state = line.strip().split("|", 1)
            rows.append({"job_id": job_id, "state": state})
    return rows


def _submit(item: Mapping[str, Any], dependency: str | None) -> str:
    exports = ["ALL", *(f"{key}={value}" for key, value in sorted(item["exports"].items()))]
    command = ["sbatch", "--parsable", "--account=nslab", "--qos=nslab"]
    if dependency is not None:
        command.append(f"--dependency=afterany:{dependency}")
    command.extend([f"--export={','.join(exports)}", str(item["_wrapper"])])
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    job_id = result.stdout.strip().split(";", 1)[0]
    if not job_id.isdigit():
        raise QueuePrimeError("sbatch did not return a numeric job identifier")
    return job_id


def submit_all(root: Path, queue: Path, state: Path, config: Mapping[str, Any], user: str) -> list[dict[str, Any]]:
    root_job_id = str(config["existing_root_job_id"])
    state.mkdir(parents=True, exist_ok=True)
    claims = state / "claims"
    claims.mkdir(exist_ok=True)
    submitted: list[dict[str, Any]] = []
    with (state / "dispatcher.lock").open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        snapshot = _gpu_snapshot(user)
        if len(snapshot) != 1 or snapshot[0]["job_id"] != root_job_id or snapshot[0]["state"] not in {"PENDING", "RUNNING"}:
            raise QueuePrimeError(f"GPU scheduler snapshot differs: {snapshot}")
        items = load_selected(root, queue, state, config)
        planned = plan_chains(items, root_job_id)
        actual_parents: dict[str, str] = {}
        for item, row in zip(items, planned, strict=True):
            dependency = row["dependency_job_id"]
            if isinstance(dependency, str) and dependency.startswith("planned:"):
                dependency = actual_parents[dependency.split(":", 1)[1]]
            marker = claims / f"{item['bundle_id']}.json"
            claim = {
                "schema_version": "masld-bench-gpu-bundle-dispatch-claim-v1",
                "bundle_id": item["bundle_id"],
                "queue_sha256": item["_queue_sha256"],
                "wrapper_sha256": item["wrapper_sha256"],
                "dispatcher_job_id": str(__import__("os").environ.get("SLURM_JOB_ID", "queue-prime")),
                "claimed_at_unix": int(time.time()),
                "status": "claimed_before_submission",
                "queue_prime_campaign_id": CAMPAIGN_ID,
                "dependency_job_id": dependency,
                "lane": row["lane"],
            }
            dispatcher.atomic_json(marker, claim)
            try:
                job_id = _submit(item, dependency)
            except Exception as error:
                claim.update(status="submission_failed", error=f"{type(error).__name__}: {error}")
                dispatcher.atomic_json(marker, claim)
                raise
            actual_parents[item["bundle_id"]] = job_id
            claim.update(status="submitted", job_id=job_id, submitted_at_unix=int(time.time()))
            dispatcher.atomic_json(marker, claim)
            dispatcher.append_ledger(state / "submissions.tsv", {"timestamp_unix": int(time.time()), "bundle_id": item["bundle_id"], "family": item["family"], "logical_tasks": item["logical_tasks"], "queue_sha256": item["_queue_sha256"], "wrapper_sha256": item["wrapper_sha256"], "job_id": job_id, "account": "nslab", "qos": "nslab", "status": "submitted_queue_prime_all_ready"})
            submitted.append({"bundle_id": item["bundle_id"], "family": item["family"], "logical_tasks": item["logical_tasks"], "job_id": job_id, "lane": row["lane"], "dependency_job_id": dependency, "queue_sha256": item["_queue_sha256"], "wrapper_sha256": item["wrapper_sha256"], "requested_gpu_hours_ceiling": item["_hours"]})
    return submitted


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--queue", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--user", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    queue = args.queue.resolve(strict=True)
    state = args.state.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    for path in (queue, state, config_path):
        path.relative_to(root)
    output = reject_symlink_components(args.output, label="GPU queue prime output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    submitted = submit_all(root, queue, state, _json(config_path), args.user)
    snapshot = _gpu_snapshot(args.user)
    if len(snapshot) != len(submitted) + 1:
        raise QueuePrimeError("post-submission GPU census differs")
    receipt = {"schema_version": "masld-bench-gpu-queue-prime-all-ready-receipt-v1", "campaign_id": CAMPAIGN_ID, "existing_root_job_id": "21083052", "new_jobs_submitted": len(submitted), "total_visible_gpu_jobs": len(snapshot), "maximum_running_campaign_jobs": 4, "pending_job_limit": None, "submissions": submitted, "scheduler_snapshot": snapshot, "innovation_used": False, "array_jobs_used": False, "sealed_outcomes_read": False}
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "tests/unit/test_prime_all_ready_gpu_queue.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gpu_queue_prime_all_ready", "campaign_id": CAMPAIGN_ID, "new_jobs_submitted": len(submitted), "maximum_running_campaign_jobs": 4, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest, "new_jobs_submitted": len(submitted)}, sort_keys=True))


if __name__ == "__main__":
    main()
