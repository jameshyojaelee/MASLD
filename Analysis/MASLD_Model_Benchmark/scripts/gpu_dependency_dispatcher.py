#!/usr/bin/env python3
"""Persistently append ready GPU bundles to four dependency lanes.

This is the successor to ``gpu_bundle_dispatcher.py``.  It intentionally has
no visible-job or pending-job ceiling.  Safety comes from serial ``afterany``
chains: every submitted campaign job belongs to one of exactly four lanes.
"""

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

from scripts import gpu_bundle_dispatcher as legacy


POLICY_SCHEMA = "masld-bench-gpu-dependency-dispatch-policy-v1"
STATE_SCHEMA = "masld-bench-gpu-dependency-lanes-v1"
CLAIM_SCHEMA = "masld-bench-gpu-dependency-dispatch-claim-v1"
GENERIC_JOB_NAME = re.compile(r"^model-(?:training|work|probe)(?:-[0-9]{3})?$")
ACTIVE_STATES = {
    "PENDING",
    "RUNNING",
    "CONFIGURING",
    "COMPLETING",
    "SUSPENDED",
    "RESIZING",
}
TERMINAL_STATES = {
    "BOOT_FAIL",
    "CANCELLED",
    "COMPLETED",
    "DEADLINE",
    "FAILED",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "PREEMPTED",
    "REVOKED",
    "SPECIAL_EXIT",
    "TIMEOUT",
}


class GPUDependencyDispatcherError(RuntimeError):
    """Raised when dependency-lane dispatch cannot be proven safe."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GPUDependencyDispatcherError(f"invalid JSON: {path}") from error
    if not isinstance(value, dict):
        raise GPUDependencyDispatcherError(f"JSON object required: {path}")
    return value


def resolve_inside(root: Path, value: str, label: str) -> Path:
    candidate = (root / value).resolve(strict=True)
    try:
        candidate.relative_to(root)
    except ValueError as error:
        raise GPUDependencyDispatcherError(f"{label} escapes benchmark root") from error
    return candidate


def sbatch_header(text: str, key: str) -> str:
    match = re.search(
        rf"^#SBATCH\s+--{re.escape(key)}(?:=|\s+)(\S+)", text, re.MULTILINE
    )
    if match is None:
        raise GPUDependencyDispatcherError(f"missing sbatch header: {key}")
    return match.group(1)


def wall_hours(value: str) -> float:
    total = 0.0
    if "-" in value:
        days, value = value.split("-", 1)
        total = 24.0 * int(days)
    fields = value.split(":")
    if len(fields) != 3:
        raise GPUDependencyDispatcherError("wall time must be HH:MM:SS or D-HH:MM:SS")
    hours, minutes, seconds = (int(field) for field in fields)
    if minutes > 59 or seconds > 59:
        raise GPUDependencyDispatcherError("invalid wall time")
    return total + hours + minutes / 60.0 + seconds / 3600.0


def validate_item_route(item: dict[str, Any]) -> dict[str, Any]:
    text = item["_wrapper"].read_text(encoding="utf-8")
    expected = {
        "partition": "gpu",
        "account": "nslab",
        "qos": "nslab",
        "gres": "gpu:l40s:1",
    }
    observed = {key: sbatch_header(text, key) for key in expected}
    if observed != expected:
        raise GPUDependencyDispatcherError(
            f"GPU wrapper routing differs for {item['bundle_id']}: {observed}"
        )
    job_name = sbatch_header(text, "job-name")
    if GENERIC_JOB_NAME.fullmatch(job_name) is None:
        raise GPUDependencyDispatcherError(
            f"GPU job name is not generic for {item['bundle_id']}: {job_name}"
        )
    if "#SBATCH --array" in text or "innovation" in text.lower():
        raise GPUDependencyDispatcherError(
            f"GPU wrapper array/QOS contract differs for {item['bundle_id']}"
        )
    item["_job_name"] = job_name
    item["_hours"] = wall_hours(sbatch_header(text, "time"))
    return item


def _validate_source_bindings(root: Path, bindings: Any) -> None:
    required = {
        "dependency_dispatcher",
        "legacy_queue_parser",
        "unit_test",
        "validation_wrapper",
        "runtime_wrapper",
    }
    if not isinstance(bindings, dict) or set(bindings) != required:
        raise GPUDependencyDispatcherError("policy source roster differs")
    for label, record in bindings.items():
        if not isinstance(record, dict) or set(record) != {"path", "sha256"}:
            raise GPUDependencyDispatcherError(f"invalid source binding: {label}")
        path = resolve_inside(root, str(record["path"]), f"source binding {label}")
        if digest(path) != record["sha256"]:
            raise GPUDependencyDispatcherError(f"policy source drifted: {label}")


def load_policy(path: Path, root: Path) -> dict[str, Any]:
    policy = load_json(path)
    expected_keys = {
        "schema_version",
        "campaign_id",
        "routing",
        "lane_count",
        "dependency_type",
        "pending_job_limit",
        "submission_batch_limit",
        "poll_seconds",
        "legacy_dispatcher_job_id",
        "bootstrap",
        "source_bindings",
    }
    if set(policy) != expected_keys or policy.get("schema_version") != POLICY_SCHEMA:
        raise GPUDependencyDispatcherError("dependency-dispatch policy schema differs")
    if policy.get("campaign_id") != "gpu_dependency_dispatcher_20260825":
        raise GPUDependencyDispatcherError("dependency-dispatch campaign differs")
    if policy.get("routing") != {
        "account": "nslab",
        "qos": "nslab",
        "partition": "gpu",
        "gres": "gpu:l40s:1",
        "job_name_pattern": GENERIC_JOB_NAME.pattern,
        "arrays_allowed": False,
        "innovation_allowed": False,
    }:
        raise GPUDependencyDispatcherError("dependency-dispatch routing differs")
    if (
        policy.get("lane_count") != 4
        or policy.get("dependency_type") != "afterany"
        or policy.get("pending_job_limit") is not None
        or policy.get("submission_batch_limit") is not None
    ):
        raise GPUDependencyDispatcherError("dependency-dispatch concurrency differs")
    if not isinstance(policy.get("poll_seconds"), int) or not 10 <= policy["poll_seconds"] <= 300:
        raise GPUDependencyDispatcherError("dependency-dispatch poll interval differs")
    if not str(policy.get("legacy_dispatcher_job_id", "")).isdigit():
        raise GPUDependencyDispatcherError("legacy dispatcher identifier differs")
    _validate_source_bindings(root, policy["source_bindings"])
    _validate_bootstrap(policy, root)
    policy["_sha256"] = digest(path)
    return policy


def _validate_bootstrap(policy: Mapping[str, Any], root: Path) -> None:
    bootstrap = policy.get("bootstrap")
    required = {
        "artifact_path",
        "artifact_sha256",
        "receipt_path",
        "receipt_sha256",
        "root_job_id",
        "root_requested_gpu_hours",
        "lane_tails",
    }
    if not isinstance(bootstrap, dict) or set(bootstrap) != required:
        raise GPUDependencyDispatcherError("bootstrap schema differs")
    artifact = resolve_inside(root, str(bootstrap["artifact_path"]), "bootstrap artifact")
    receipt_path = resolve_inside(root, str(bootstrap["receipt_path"]), "bootstrap receipt")
    if digest(artifact) != bootstrap["artifact_sha256"] or digest(receipt_path) != bootstrap["receipt_sha256"]:
        raise GPUDependencyDispatcherError("bootstrap artifact drifted")
    receipt = load_json(receipt_path)
    if (
        receipt.get("maximum_running_campaign_jobs") != 4
        or receipt.get("pending_job_limit") is not None
        or receipt.get("innovation_used") is not False
        or receipt.get("array_jobs_used") is not False
        or receipt.get("sealed_outcomes_read") is not False
    ):
        raise GPUDependencyDispatcherError("bootstrap receipt policy differs")
    root_job_id = str(bootstrap["root_job_id"])
    if root_job_id != str(receipt.get("existing_root_job_id")) or not root_job_id.isdigit():
        raise GPUDependencyDispatcherError("bootstrap root differs")
    submissions = receipt.get("submissions")
    if not isinstance(submissions, list) or not submissions:
        raise GPUDependencyDispatcherError("bootstrap submissions are absent")
    parents: list[str | None] = [root_job_id, None, None, None]
    root_count = 1
    seen_jobs = {root_job_id}
    loads = [float(bootstrap["root_requested_gpu_hours"]), 0.0, 0.0, 0.0]
    bundle_tails: list[str | None] = ["model-work-208", None, None, None]
    for row in submissions:
        if not isinstance(row, dict) or not isinstance(row.get("lane"), int):
            raise GPUDependencyDispatcherError("bootstrap submission row differs")
        lane = row["lane"]
        if lane not in range(4):
            raise GPUDependencyDispatcherError("bootstrap lane differs")
        dependency = row.get("dependency_job_id")
        if dependency != parents[lane]:
            raise GPUDependencyDispatcherError("bootstrap dependency chain is discontinuous")
        job_id = str(row.get("job_id", ""))
        if not job_id.isdigit() or job_id in seen_jobs:
            raise GPUDependencyDispatcherError("bootstrap job identity differs")
        if dependency is None:
            root_count += 1
        parents[lane] = job_id
        seen_jobs.add(job_id)
        loads[lane] += float(row.get("requested_gpu_hours_ceiling", 0.0))
        bundle_tails[lane] = str(row.get("bundle_id", ""))
    if root_count != 4:
        raise GPUDependencyDispatcherError("bootstrap does not have exactly four roots")
    expected_tails = [
        {
            "lane": lane,
            "tail_job_id": parents[lane],
            "tail_bundle_id": bundle_tails[lane],
            "requested_gpu_hours": loads[lane],
        }
        for lane in range(4)
    ]
    if bootstrap["lane_tails"] != expected_tails:
        raise GPUDependencyDispatcherError("bootstrap lane tails differ")


def bootstrap_lane_state(policy: Mapping[str, Any]) -> dict[str, Any]:
    return {
        "schema_version": STATE_SCHEMA,
        "campaign_id": policy["campaign_id"],
        "policy_sha256": policy["_sha256"],
        "revision": 0,
        "updated_at_unix": None,
        "lanes": [dict(row) for row in policy["bootstrap"]["lane_tails"]],
    }


def load_lane_state(state_dir: Path, policy: Mapping[str, Any]) -> tuple[dict[str, Any], bool]:
    path = state_dir / "dependency_lanes.json"
    if not path.exists():
        return bootstrap_lane_state(policy), False
    state = load_json(path)
    if (
        state.get("schema_version") != STATE_SCHEMA
        or state.get("campaign_id") != policy["campaign_id"]
        or state.get("policy_sha256") != policy["_sha256"]
        or not isinstance(state.get("revision"), int)
        or not isinstance(state.get("lanes"), list)
        or len(state["lanes"]) != 4
    ):
        raise GPUDependencyDispatcherError("persisted dependency-lane state differs")
    for lane, row in enumerate(state["lanes"]):
        if (
            not isinstance(row, dict)
            or row.get("lane") != lane
            or not str(row.get("tail_job_id", "")).isdigit()
            or not isinstance(row.get("tail_bundle_id"), str)
            or not isinstance(row.get("requested_gpu_hours"), (int, float))
        ):
            raise GPUDependencyDispatcherError("persisted dependency lane differs")
    return state, True


def scheduler_snapshot(user: str) -> list[dict[str, str]]:
    result = subprocess.run(
        [
            "squeue",
            "-u",
            user,
            "-p",
            "gpu",
            "-h",
            "-o",
            "%i|%u|%T|%a|%q|%P|%j|%b|%R|%E",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    rows: list[dict[str, str]] = []
    keys = ("job_id", "user", "state", "account", "qos", "partition", "job_name", "gres", "reason", "dependency")
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        fields = line.rstrip().split("|", 9)
        if len(fields) != len(keys):
            raise GPUDependencyDispatcherError("unexpected squeue row")
        rows.append(dict(zip(keys, fields, strict=True)))
    return rows


def legacy_dispatcher_active(job_id: str) -> bool:
    result = subprocess.run(
        ["squeue", "-j", job_id, "-h", "-o", "%T"],
        check=False,
        text=True,
        capture_output=True,
    )
    if result.returncode == 1 and not result.stdout.strip() and "Invalid job id specified" in result.stderr:
        return False
    if result.returncode != 0:
        raise GPUDependencyDispatcherError(
            f"legacy-dispatcher scheduler query failed: {result.stderr.strip()}"
        )
    return any(state.strip() in ACTIVE_STATES for state in result.stdout.splitlines())


def _tracked_lane_job_ids(state_dir: Path, policy: Mapping[str, Any]) -> set[str]:
    identifiers = {str(policy["bootstrap"]["root_job_id"])}
    receipt = load_json(resolve_inside(Path(policy["_root"]), policy["bootstrap"]["receipt_path"], "bootstrap receipt"))
    identifiers.update(str(row["job_id"]) for row in receipt["submissions"])
    claims = state_dir / "claims"
    if claims.is_dir():
        for path in claims.glob("*.json"):
            value = load_json(path)
            if value.get("schema_version") == CLAIM_SCHEMA and value.get("status") == "submitted":
                identifiers.add(str(value["job_id"]))
    return identifiers


def _campaign_job_ids(state_dir: Path, tracked_lane_job_ids: set[str]) -> set[str]:
    identifiers = set(tracked_lane_job_ids)
    ledger = state_dir / "submissions.tsv"
    if ledger.is_file():
        with ledger.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                job_id = str(row.get("job_id", ""))
                if job_id.isdigit():
                    identifiers.add(job_id)
    return identifiers


def audit_live_campaign(
    rows: Sequence[Mapping[str, str]],
    campaign_job_ids: set[str],
    tracked_lane_job_ids: set[str] | None = None,
) -> dict[str, int]:
    campaign = [row for row in rows if row.get("job_id") in campaign_job_ids]
    if tracked_lane_job_ids is not None:
        untracked = sorted(
            str(row["job_id"])
            for row in campaign
            if row.get("job_id") not in tracked_lane_job_ids
        )
        if untracked:
            raise GPUDependencyDispatcherError(
                "live campaign jobs are outside the four dependency lanes: "
                + ",".join(untracked)
            )
    for row in campaign:
        if (
            row.get("account") != "nslab"
            or row.get("qos") != "nslab"
            or row.get("partition") != "gpu"
            or "l40s" not in str(row.get("gres", "")).lower()
            or GENERIC_JOB_NAME.fullmatch(str(row.get("job_name", ""))) is None
        ):
            raise GPUDependencyDispatcherError(
                f"live campaign GPU routing differs: {row.get('job_id')}"
            )
    def dependency_pending(row: Mapping[str, str]) -> bool:
        reason = str(row.get("reason", "")).strip().strip("()")
        return row.get("state") == "PENDING" and reason == "Dependency"

    runnable = sum(
        row.get("state") in ACTIVE_STATES and not dependency_pending(row)
        for row in campaign
    )
    if runnable > 4:
        raise GPUDependencyDispatcherError(
            f"more than four campaign jobs are runnable/running: {runnable}"
        )
    return {
        "visible_campaign_jobs": len(campaign),
        "runnable_or_running_campaign_jobs": runnable,
        "dependency_pending_campaign_jobs": sum(
            dependency_pending(row) for row in campaign
        ),
    }


def assert_no_incomplete_dependency_claims(state_dir: Path) -> None:
    claims = state_dir / "claims"
    if not claims.is_dir():
        return
    incomplete = []
    for path in sorted(claims.glob("*.json")):
        value = load_json(path)
        if value.get("schema_version") == CLAIM_SCHEMA and value.get("status") != "submitted":
            incomplete.append(path.name)
    if incomplete:
        raise GPUDependencyDispatcherError(
            "incomplete dependency-dispatch claims require reconciliation: "
            + ",".join(incomplete)
        )


def _terminal_tail_state(job_id: str, user: str) -> str:
    result = subprocess.run(
        [
            "sacct",
            "-X",
            "-j",
            job_id,
            "--noheader",
            "--parsable2",
            "-o",
            "JobIDRaw,User,State",
        ],
        check=True,
        text=True,
        capture_output=True,
    )
    matches = []
    for line in result.stdout.splitlines():
        fields = line.strip().split("|")
        if len(fields) >= 3 and fields[0] == job_id and fields[1] == user:
            matches.append(fields[2].split("+", 1)[0])
    if len(matches) != 1 or matches[0] not in TERMINAL_STATES:
        raise GPUDependencyDispatcherError(
            f"tail job is neither live nor uniquely terminal: {job_id}"
        )
    return matches[0]


def tail_activity(
    lane_state: Mapping[str, Any], rows: Sequence[Mapping[str, str]], user: str
) -> dict[int, bool]:
    live = {str(row["job_id"]): str(row["state"]) for row in rows}
    activity: dict[int, bool] = {}
    for lane in lane_state["lanes"]:
        job_id = str(lane["tail_job_id"])
        if job_id in live:
            if live[job_id] not in ACTIVE_STATES:
                raise GPUDependencyDispatcherError(f"unexpected live tail state: {job_id}")
            activity[int(lane["lane"])] = True
        else:
            _terminal_tail_state(job_id, user)
            activity[int(lane["lane"])] = False
    return activity


def plan_items(
    items: Sequence[Mapping[str, Any]],
    lane_state: Mapping[str, Any],
    active_tails: Mapping[int, bool],
) -> list[dict[str, Any]]:
    if set(active_tails) != set(range(4)):
        raise GPUDependencyDispatcherError("tail activity roster differs")
    parents: list[str | None] = []
    loads: list[float] = []
    for lane in lane_state["lanes"]:
        index = int(lane["lane"])
        parents.append(str(lane["tail_job_id"]) if active_tails[index] else None)
        loads.append(float(lane["requested_gpu_hours"]) if active_tails[index] else 0.0)
    planned: list[dict[str, Any]] = []
    for item in items:
        lane = min(range(4), key=lambda index: (loads[index], index))
        dependency = parents[lane]
        planned.append(
            {
                "bundle_id": item["bundle_id"],
                "lane": lane,
                "dependency_job_id": dependency,
                "reset_lane_load": dependency is None,
            }
        )
        parents[lane] = f"planned:{item['bundle_id']}"
        loads[lane] += float(item["_hours"])
    return planned


def ready_unclaimed_items(root: Path, queue: Path, state_dir: Path) -> list[dict[str, Any]]:
    claims = state_dir / "claims"
    values = []
    for path in sorted(queue.glob("*.json")):
        summary = load_json(path)
        bundle_id = summary.get("bundle_id")
        if not isinstance(bundle_id, str):
            raise GPUDependencyDispatcherError(f"queue bundle identity differs: {path.name}")
        # Claimed work is read-only history.  Do not re-admit its wrapper after
        # source evolution; the shared claim is the duplicate-prevention check.
        if (claims / f"{bundle_id}.json").exists() or summary.get("enabled") is False:
            continue
        item = legacy.load_item(path, root)
        if (
            item["enabled"] is True
            and item["_ready"] is True
        ):
            values.append(validate_item_route(item))
    values.sort(key=lambda item: (item["priority"], item["bundle_id"]))
    return values


def _submit(item: Mapping[str, Any], dependency: str | None) -> str:
    exports = ["ALL", *(f"{key}={value}" for key, value in sorted(item["exports"].items()))]
    command = [
        "sbatch",
        "--parsable",
        "--partition=gpu",
        "--account=nslab",
        "--qos=nslab",
        "--gres=gpu:l40s:1",
        f"--job-name={item['_job_name']}",
        f"--comment=mbgpu-{item['_queue_sha256'][:24]}",
    ]
    if dependency is not None:
        command.append(f"--dependency=afterany:{dependency}")
    command.extend([f"--export={','.join(exports)}", str(item["_wrapper"])])
    result = subprocess.run(command, check=True, text=True, capture_output=True)
    job_id = result.stdout.strip().split(";", 1)[0]
    if not job_id.isdigit():
        raise GPUDependencyDispatcherError("sbatch did not return a numeric job identifier")
    return job_id


def _persist_lane_state(
    state_dir: Path,
    state: dict[str, Any],
    lane: int,
    item: Mapping[str, Any],
    job_id: str,
    *,
    reset_lane_load: bool,
) -> None:
    previous = 0.0 if reset_lane_load else float(
        state["lanes"][lane]["requested_gpu_hours"]
    )
    state["lanes"][lane] = {
        "lane": lane,
        "tail_job_id": job_id,
        "tail_bundle_id": item["bundle_id"],
        "requested_gpu_hours": previous + float(item["_hours"]),
    }
    state["revision"] += 1
    state["updated_at_unix"] = int(time.time())
    legacy.atomic_json(state_dir / "dependency_lanes.json", state)


def dispatch_once(
    *,
    root: Path,
    queue: Path,
    state_dir: Path,
    policy: Mapping[str, Any],
    user: str,
    audit_only: bool,
) -> dict[str, Any]:
    root = root.resolve(strict=True)
    state_dir.mkdir(parents=True, exist_ok=True)
    claims = state_dir / "claims"
    claims.mkdir(exist_ok=True)
    with (state_dir / "dispatcher.lock").open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        assert_no_incomplete_dependency_claims(state_dir)
        rows = scheduler_snapshot(user)
        tracked_lane_jobs = _tracked_lane_job_ids(state_dir, policy)
        census = audit_live_campaign(
            rows,
            _campaign_job_ids(state_dir, tracked_lane_jobs),
            tracked_lane_jobs,
        )
        lane_state, persisted = load_lane_state(state_dir, policy)
        active_tails = tail_activity(lane_state, rows, user)
        items = ready_unclaimed_items(root, queue, state_dir)
        plan = plan_items(items, lane_state, active_tails)
        legacy_active = legacy_dispatcher_active(str(policy["legacy_dispatcher_job_id"]))
        if audit_only:
            return {
                "audit_only": True,
                "legacy_dispatcher_active": legacy_active,
                "lane_state_persisted": persisted,
                "ready_unclaimed_jobs": len(items),
                "planned_jobs": len(plan),
                **census,
            }
        if legacy_active:
            raise GPUDependencyDispatcherError(
                "legacy dispatcher remains active; refusing dual-dispatch activation"
            )
        actual_parents: dict[str, str] = {}
        submitted = []
        for item, row in zip(items, plan, strict=True):
            dependency = row["dependency_job_id"]
            if isinstance(dependency, str) and dependency.startswith("planned:"):
                dependency = actual_parents[dependency.split(":", 1)[1]]
            marker = claims / f"{item['bundle_id']}.json"
            claim = {
                "schema_version": CLAIM_SCHEMA,
                "campaign_id": policy["campaign_id"],
                "bundle_id": item["bundle_id"],
                "queue_sha256": item["_queue_sha256"],
                "wrapper_sha256": item["wrapper_sha256"],
                "dispatcher_job_id": os.environ.get("SLURM_JOB_ID", "interactive"),
                "claimed_at_unix": int(time.time()),
                "status": "claimed_before_submission",
                "lane": row["lane"],
                "dependency_job_id": dependency,
                "submission_comment": f"mbgpu-{item['_queue_sha256'][:24]}",
            }
            legacy.atomic_json(marker, claim)
            try:
                job_id = _submit(item, dependency)
            except Exception as error:
                claim.update(
                    status="submission_indeterminate",
                    error=f"{type(error).__name__}: {error}",
                )
                legacy.atomic_json(marker, claim)
                raise
            _persist_lane_state(
                state_dir,
                lane_state,
                int(row["lane"]),
                item,
                job_id,
                reset_lane_load=bool(row["reset_lane_load"]),
            )
            claim.update(status="submitted", job_id=job_id, submitted_at_unix=int(time.time()))
            legacy.atomic_json(marker, claim)
            legacy.append_ledger(
                state_dir / "submissions.tsv",
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
                    "status": "submitted_dependency_lane",
                },
            )
            actual_parents[item["bundle_id"]] = job_id
            submitted.append(
                {
                    "bundle_id": item["bundle_id"],
                    "job_id": job_id,
                    "lane": row["lane"],
                    "dependency_job_id": dependency,
                }
            )
        return {"audit_only": False, "submitted": submitted, **census}


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--queue", type=Path, required=True)
    parser.add_argument("--state", type=Path, required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--user", required=True)
    parser.add_argument("--audit-only", action="store_true")
    parser.add_argument("--once", action="store_true")
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    queue = args.queue.resolve(strict=True)
    state_dir = args.state.resolve(strict=True)
    policy_path = args.policy.resolve(strict=True)
    for path in (queue, state_dir, policy_path):
        path.relative_to(root)
    policy = load_policy(policy_path, root)
    policy["_root"] = str(root)
    while True:
        receipt = dispatch_once(
            root=root,
            queue=queue,
            state_dir=state_dir,
            policy=policy,
            user=args.user,
            audit_only=args.audit_only,
        )
        print(json.dumps(receipt, sort_keys=True), flush=True)
        if args.once or args.audit_only:
            break
        time.sleep(int(policy["poll_seconds"]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
