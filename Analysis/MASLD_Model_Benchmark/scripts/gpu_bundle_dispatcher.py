#!/usr/bin/env python3
"""Keep at most five GPU allocations by dispatching frozen sequential bundles."""

from __future__ import annotations

import argparse
import csv
import fcntl
from hashlib import sha256
import json
import os
from pathlib import Path
import re
import subprocess
import time
from typing import Any, Mapping, Sequence


SCHEMA = "masld-bench-gpu-bundle-queue-item-v1"
BUNDLE_ID = re.compile(r"^[a-z0-9][a-z0-9._-]{2,79}$")
EXPORT_KEY = re.compile(r"^[A-Z][A-Z0-9_]{0,79}$")


class GPUDispatcherError(RuntimeError):
    """Raised when queue authority or scheduler state is unsafe."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def atomic_json(path: Path, value: Mapping[str, Any]) -> None:
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temporary.open("x", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(temporary, path)


def resolve_inside(root: Path, value: str, label: str) -> Path:
    candidate = (root / value).resolve(strict=True)
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise GPUDispatcherError(f"{label} escapes benchmark root") from error
    return candidate


def resolve_maybe_inside(root: Path, value: str, label: str) -> Path:
    candidate = (root / value).resolve(strict=False)
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise GPUDispatcherError(f"{label} escapes benchmark root") from error
    return candidate


def load_item(path: Path, root: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GPUDispatcherError(f"invalid queue item {path.name}") from error
    required = {
        "schema_version",
        "bundle_id",
        "priority",
        "enabled",
        "wrapper_path",
        "wrapper_sha256",
        "exports",
        "required_paths",
        "logical_tasks",
        "family",
    }
    if not isinstance(value, dict) or set(value) != required or value["schema_version"] != SCHEMA:
        raise GPUDispatcherError(f"queue item schema differs: {path.name}")
    if not isinstance(value["bundle_id"], str) or not BUNDLE_ID.fullmatch(value["bundle_id"]):
        raise GPUDispatcherError("bundle identifier differs")
    if not isinstance(value["priority"], int) or value["priority"] < 0:
        raise GPUDispatcherError("bundle priority differs")
    if not isinstance(value["enabled"], bool):
        raise GPUDispatcherError("bundle enabled flag differs")
    if not isinstance(value["logical_tasks"], int) or value["logical_tasks"] < 1:
        raise GPUDispatcherError("bundle logical-task count differs")
    if not isinstance(value["family"], str) or not value["family"]:
        raise GPUDispatcherError("bundle family differs")
    wrapper = resolve_inside(root, value["wrapper_path"], "bundle wrapper")
    if not wrapper.is_file() or digest(wrapper) != value["wrapper_sha256"]:
        raise GPUDispatcherError("bundle wrapper checksum differs")
    text = wrapper.read_text(encoding="utf-8")
    if (
        "#SBATCH --partition=gpu" not in text
        or "#SBATCH --qos=nslab" not in text
        or "gpu:l40s:1" not in text
        or "#SBATCH --array" in text
        or "--qos=innovation" in text
    ):
        raise GPUDispatcherError("bundle wrapper GPU/QOS/array contract differs")
    exports = value["exports"]
    if not isinstance(exports, dict):
        raise GPUDispatcherError("bundle exports must be an object")
    for key, item in exports.items():
        if not isinstance(key, str) or not EXPORT_KEY.fullmatch(key):
            raise GPUDispatcherError("bundle export key differs")
        if not isinstance(item, str) or not item or any(token in item for token in ("\n", "\r", ",")):
            raise GPUDispatcherError("bundle export value differs")
    required_paths = value["required_paths"]
    if not isinstance(required_paths, list) or any(not isinstance(item, str) for item in required_paths):
        raise GPUDispatcherError("bundle required paths differ")
    resolved_required_paths = [
        resolve_maybe_inside(root, item, "bundle required path") for item in required_paths
    ]
    value["_wrapper"] = wrapper
    value["_queue_sha256"] = digest(path)
    value["_ready"] = all(item.is_file() for item in resolved_required_paths)
    return value


def gpu_job_states(user: str) -> list[str]:
    result = subprocess.run(
        ["squeue", "-u", user, "-p", "gpu", "-h", "-o", "%T"],
        check=True,
        text=True,
        capture_output=True,
    )
    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def submission_allowed(
    states: Sequence[str], *, cap: int, max_running: int, max_pending: int
) -> bool:
    if not 1 <= cap <= 5 or not 1 <= max_running <= cap or max_pending != 1:
        raise GPUDispatcherError("GPU running/pending policy differs")
    pending = sum(state == "PENDING" for state in states)
    active = len(states) - pending
    # The submitted job can start immediately.  Require an open running slot
    # before submission instead of assuming the new allocation will pend.
    return len(states) < cap and active < max_running and pending < max_pending


def append_ledger(path: Path, row: Mapping[str, object]) -> None:
    fields = (
        "timestamp_unix",
        "bundle_id",
        "family",
        "logical_tasks",
        "queue_sha256",
        "wrapper_sha256",
        "job_id",
        "account",
        "qos",
        "status",
    )
    exists = path.exists()
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        if not exists:
            writer.writeheader()
        writer.writerow(row)
        handle.flush()
        os.fsync(handle.fileno())


def submit(item: Mapping[str, Any]) -> str:
    exports = ["ALL", *(f"{key}={value}" for key, value in sorted(item["exports"].items()))]
    command = [
        "sbatch",
        "--parsable",
        "--account=nslab",
        "--qos=nslab",
        f"--export={','.join(exports)}",
        str(item["_wrapper"]),
    ]
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    job_id = result.stdout.strip().split(";", 1)[0]
    if not job_id.isdigit():
        raise GPUDispatcherError("sbatch did not return a numeric job identifier")
    return job_id


def dispatch_once(
    *,
    root: Path,
    queue: Path,
    state: Path,
    user: str,
    cap: int,
    max_running: int,
    max_pending: int,
) -> list[str]:
    root = root.resolve(strict=True)
    queue.mkdir(parents=True, exist_ok=True)
    state.mkdir(parents=True, exist_ok=True)
    claims = state / "claims"
    claims.mkdir(exist_ok=True)
    submitted: list[str] = []
    with (state / "dispatcher.lock").open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        current_states = gpu_job_states(user)
        if not submission_allowed(
            current_states,
            cap=cap,
            max_running=max_running,
            max_pending=max_pending,
        ):
            return submitted
        items = [load_item(path, root) for path in sorted(queue.glob("*.json"))]
        items.sort(key=lambda item: (item["priority"], item["bundle_id"]))
        for item in items:
            claim_retry_of = None
            if not submission_allowed(
                current_states,
                cap=cap,
                max_running=max_running,
                max_pending=max_pending,
            ):
                break
            marker = claims / f"{item['bundle_id']}.json"
            if not item["enabled"] or not item["_ready"]:
                continue
            # Recovery reads the marker's STATUS, not merely its existence. A
            # hard kill between writing the claim and submitting leaves the
            # marker at claimed_before_submission; treating that as "already
            # handled" retires the bundle silently and forever. The
            # claim-before-submit ordering itself is a write-ahead intent
            # record and is correct -- do not reorder it to fix this.
            if marker.exists():
                try:
                    prior = json.loads(marker.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError):
                    prior = {"status": "unreadable"}
                if prior.get("status") != "claimed_before_submission":
                    continue
                # An interrupted claim is retried, and the retry is recorded
                # rather than overwriting the evidence that it happened.
                claim_retry_of = prior
            claim = {
                "schema_version": "masld-bench-gpu-bundle-dispatch-claim-v1",
                "retry_of_interrupted_claim": locals().get("claim_retry_of") is not None,
                "bundle_id": item["bundle_id"],
                "queue_sha256": item["_queue_sha256"],
                "wrapper_sha256": item["wrapper_sha256"],
                "dispatcher_job_id": os.environ.get("SLURM_JOB_ID", "interactive"),
                "claimed_at_unix": int(time.time()),
                "status": "claimed_before_submission",
            }
            atomic_json(marker, claim)
            try:
                job_id = submit(item)
            except Exception as error:
                claim.update(status="submission_failed", error=f"{type(error).__name__}: {error}")
                atomic_json(marker, claim)
                append_ledger(
                    state / "submissions.tsv",
                    {
                        "timestamp_unix": int(time.time()),
                        "bundle_id": item["bundle_id"],
                        "family": item["family"],
                        "logical_tasks": item["logical_tasks"],
                        "queue_sha256": item["_queue_sha256"],
                        "wrapper_sha256": item["wrapper_sha256"],
                        "job_id": "",
                        "account": "nslab",
                        "qos": "nslab",
                        "status": "submission_failed",
                    },
                )
                continue
            claim.update(
                status="submitted",
                job_id=job_id,
                submitted_at_unix=int(time.time()),
                account="nslab",
                qos="nslab",
            )
            atomic_json(marker, claim)
            append_ledger(
                state / "submissions.tsv",
                {
                    "timestamp_unix": int(time.time()),
                    "bundle_id": item["bundle_id"],
                    "family": item["family"],
                    "logical_tasks": item["logical_tasks"],
                    "queue_sha256": item["_queue_sha256"],
                    "wrapper_sha256": item["wrapper_sha256"],
                    "job_id": job_id,
                    "account": "nslab",
                    "qos": "nslab",
                    "status": "submitted",
                },
            )
            submitted.append(job_id)
            current_states.append("PENDING")
    return submitted


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--queue", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--user", required=True)
    parser.add_argument("--cap", type=int, default=5)
    parser.add_argument("--max-running", type=int, default=4)
    parser.add_argument("--max-pending", type=int, default=1)
    parser.add_argument("--poll-seconds", type=int, default=30)
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()
    if not 10 <= args.poll_seconds <= 300:
        raise GPUDispatcherError("poll interval must be between 10 and 300 seconds")
    while True:
        jobs = dispatch_once(
            root=args.root,
            queue=args.queue,
            state=args.state,
            user=args.user,
            cap=args.cap,
            max_running=args.max_running,
            max_pending=args.max_pending,
        )
        if jobs:
            print(json.dumps({"submitted": jobs, "timestamp_unix": int(time.time())}), flush=True)
        if args.once:
            break
        time.sleep(args.poll_seconds)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
