#!/usr/bin/env python3
"""Reconcile corrected five-seed queue sources with the live dispatcher."""

from __future__ import annotations

import argparse
import csv
import fcntl
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
import subprocess
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import verify_frozen_tree


SCHEMA = "masld-bench-sequence-five-seed-live-dispatch-reconciliation-v1"
CLAIM_STATES = {"claimed_before_submission", "submitted", "submission_failed"}


class ReconciliationError(RuntimeError):
    """Raised when corrected queue or dispatcher safety differs."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--correction-root", type=Path, required=True)
    parser.add_argument("--queue-admission-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--user", required=True)
    return parser.parse_args()


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def write_tsv(
    path: Path,
    fields: Sequence[str],
    rows: Iterable[Mapping[str, object]],
) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def load_contract(path: Path) -> dict[str, Any]:
    contract = json.loads(path.read_text(encoding="utf-8"))
    if contract.get("schema_version") != SCHEMA:
        raise ReconciliationError("live dispatch reconciliation schema differs")
    queue = contract.get("queue", {})
    v1 = queue.get("v1_dispatcher_ids", ())
    v2 = queue.get("v2_dispatcher_ids", ())
    if (
        v1 != [f"model-training-{value}" for value in range(701, 731)]
        or v2 != [f"model-training-{value}" for value in range(801, 831)]
        or set(v1) & set(v2)
        or queue.get("logical_fits") != 438
        or queue.get("v1_required_enabled_state") is not False
        or queue.get("v2_required_enabled_state") is not True
        or queue.get("v1_claims_allowed") is not False
        or queue.get("v1_ledger_submissions_allowed") is not False
    ):
        raise ReconciliationError("v1/v2 identity or disposition contract differs")
    policy = contract.get("gpu_policy", {})
    if (
        policy.get("submission_authority") != "central_dispatcher_only"
        or policy.get("account") != "nslab"
        or policy.get("qos") != "nslab"
        or policy.get("gres") != "gpu:l40s:1"
        or policy.get("gpu_job_ceiling") != 5
        or policy.get("maximum_running") != 4
        or policy.get("maximum_pending") != 1
        or policy.get("manual_gpu_submission_allowed") is not False
        or policy.get("innovation_allowed") is not False
        or policy.get("array_allowed") is not False
    ):
        raise ReconciliationError("GPU policy differs")
    return contract


def load_module(path: Path, name: str) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ReconciliationError(f"cannot load source: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_controller(root: Path, contract: Mapping[str, Any]) -> dict[str, object]:
    sources = contract["frozen_sources"]
    paths = {
        "v1_controller": root / "scripts/run_sequence_task_native_five_seed_bundle.py",
        "v2_controller": root / "scripts/run_sequence_task_native_five_seed_bundle_v2.py",
        "v1_wrapper": root / "slurm/run_sequence_task_native_five_seed_bundle.sbatch",
        "v2_wrapper": root / "slurm/run_sequence_task_native_five_seed_bundle_v2.sbatch",
        "dispatcher": root / "scripts/gpu_bundle_dispatcher.py",
    }
    for name, path in paths.items():
        if digest(path) != sources[f"{name}_sha256"]:
            raise ReconciliationError(f"frozen source differs: {name}")
    controller = load_module(paths["v2_controller"], "sequence_five_seed_controller_v2")
    if controller.Counter(["passed", "passed"]) != {"passed": 2}:
        raise ReconciliationError("v2 controller Counter import is incomplete")
    passed_counts, passed_status = controller.terminal_summary({
        "one": {"status": "passed"},
        "two": {"status": "recovered_passed"},
    })
    failed_counts, failed_status = controller.terminal_summary({
        "one": {"status": "passed"},
        "two": {"status": "failed_missing_artifact"},
    })
    if (
        passed_counts != {"passed": 1, "recovered_passed": 1}
        or passed_status != "passed"
        or failed_counts != {"passed": 1, "failed_missing_artifact": 1}
        or failed_status != "completed_with_fit_failures"
    ):
        raise ReconciliationError("v2 terminal summary is incomplete")
    wrapper_text = paths["v2_wrapper"].read_text(encoding="utf-8")
    if (
        "#SBATCH --job-name=model-training-700" not in wrapper_text
        or "#SBATCH --partition=gpu" not in wrapper_text
        or "#SBATCH --account=nslab" not in wrapper_text
        or "#SBATCH --qos=nslab" not in wrapper_text
        or "#SBATCH --gres=gpu:l40s:1" not in wrapper_text
        or "#SBATCH --array" in wrapper_text
        or "--qos=innovation" in wrapper_text
    ):
        raise ReconciliationError("v2 wrapper generic-name or GPU policy differs")
    return {
        "counter_import_complete": True,
        "pass_terminal_summary_executed": True,
        "failure_terminal_summary_executed": True,
        "v2_wrapper_generic_nslab_l40s": True,
    }


def validate_frozen_roots(
    *,
    root: Path,
    contract: Mapping[str, Any],
    plan_root: Path,
    correction_root: Path,
    queue_admission_root: Path,
) -> None:
    sources = contract["frozen_sources"]
    bindings = (
        (plan_root, sources["rectangle_plan_artifacts_sha256"]),
        (correction_root, sources["execution_correction_artifacts_sha256"]),
        (queue_admission_root, sources["corrected_queue_artifacts_sha256"]),
    )
    for path, expected in bindings:
        if root not in path.parents:
            raise ReconciliationError("frozen source escapes benchmark root")
        if digest(path / "ARTIFACTS.json") != expected:
            raise ReconciliationError(f"frozen artifact binding differs: {path.name}")
        verify_frozen_tree(path)


def load_claim(path: Path) -> dict[str, Any]:
    claim = json.loads(path.read_text(encoding="utf-8"))
    if (
        claim.get("schema_version") != "masld-bench-gpu-bundle-dispatch-claim-v1"
        or claim.get("status") not in CLAIM_STATES
    ):
        raise ReconciliationError(f"dispatcher claim schema differs: {path.name}")
    return claim


def gpu_scheduler_snapshot(user: str) -> list[dict[str, str]]:
    result = subprocess.run(
        ["squeue", "-u", user, "-p", "gpu", "-h", "-o", "%i|%T|%j|%a|%q"],
        check=True,
        text=True,
        capture_output=True,
    )
    rows: list[dict[str, str]] = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        values = line.strip().split("|")
        if len(values) != 5:
            raise ReconciliationError("GPU scheduler snapshot schema differs")
        rows.append(dict(zip(("job_id", "state", "job_name", "account", "qos"), values, strict=True)))
    return rows


def dispatcher_snapshot(user: str) -> list[dict[str, str]]:
    result = subprocess.run(
        ["squeue", "-u", user, "-h", "-o", "%i|%T|%j|%P"],
        check=True,
        text=True,
        capture_output=True,
    )
    rows: list[dict[str, str]] = []
    for line in result.stdout.splitlines():
        if not line.strip():
            continue
        values = line.strip().split("|")
        if len(values) != 4:
            raise ReconciliationError("dispatcher scheduler snapshot schema differs")
        rows.append(dict(zip(("job_id", "state", "job_name", "partition"), values, strict=True)))
    return rows


def reconcile(
    *,
    root: Path,
    contract: Mapping[str, Any],
    queue_admission_root: Path,
    user: str,
) -> tuple[dict[str, Any], list[dict[str, object]], list[dict[str, object]], list[dict[str, str]], list[dict[str, str]]]:
    dispatcher = load_module(root / "scripts/gpu_bundle_dispatcher.py", "gpu_bundle_dispatcher_reconciliation")
    live_queue = root / "config/campaigns/gpu_bundle_queue"
    claims_root = root / "executions/gpu-bundle-dispatch-state/claims"
    state_root = claims_root.parent
    frozen_queue = queue_admission_root / "queue_items"
    manifest = json.loads((frozen_queue / "manifest.json").read_text(encoding="utf-8"))
    v1_ids = set(contract["queue"]["v1_dispatcher_ids"])
    v2_ids = set(contract["queue"]["v2_dispatcher_ids"])

    all_raw = []
    for path in sorted(live_queue.glob("*.json")):
        value = json.loads(path.read_text(encoding="utf-8"))
        all_raw.append((path, value))
    all_ids = [value.get("bundle_id") for _path, value in all_raw]
    if any(not isinstance(value, str) for value in all_ids) or len(all_ids) != len(set(all_ids)):
        raise ReconciliationError("live dispatcher bundle identities collide")

    v1_rows: list[dict[str, object]] = []
    live_v1 = [(path, value) for path, value in all_raw if value.get("bundle_id") in v1_ids]
    if len(live_v1) != 30 or {value["bundle_id"] for _path, value in live_v1} != v1_ids:
        raise ReconciliationError("live v1 queue roster differs")
    for path, value in sorted(live_v1):
        item = dispatcher.load_item(path, root)
        claim_path = claims_root / f"{item['bundle_id']}.json"
        if item["enabled"] or claim_path.exists():
            raise ReconciliationError("v1 item is enabled or claimed")
        v1_rows.append({
            "bundle_id": item["bundle_id"],
            "queue_path": str(path.resolve()),
            "queue_sha256": item["_queue_sha256"],
            "enabled": item["enabled"],
            "ready": item["_ready"],
            "claim_exists": False,
        })

    ledger = read_tsv(state_root / "submissions.tsv")
    if any(row["bundle_id"] in v1_ids for row in ledger):
        raise ReconciliationError("v1 dispatcher identity appears in submission ledger")

    bindings = {row["dispatcher_bundle_id"]: row for row in manifest["bindings"]}
    v2_rows: list[dict[str, object]] = []
    claim_rows: list[dict[str, object]] = []
    live_v2 = [(path, value) for path, value in all_raw if value.get("bundle_id") in v2_ids]
    if len(live_v2) != 30 or {value["bundle_id"] for _path, value in live_v2} != v2_ids:
        raise ReconciliationError("live v2 queue roster differs")
    for path, value in sorted(live_v2):
        frozen_matches = list(frozen_queue.glob(f"*-{value['bundle_id']}-v2.json"))
        if len(frozen_matches) != 1 or path.read_bytes() != frozen_matches[0].read_bytes():
            raise ReconciliationError(f"live v2 queue differs from frozen admission: {path.name}")
        item = dispatcher.load_item(path, root)
        binding = bindings.get(item["bundle_id"])
        if (
            not item["enabled"]
            or not item["_ready"]
            or binding is None
            or item["exports"].get("RECTANGLE_BUNDLE_ID") != binding["task_bundle_id"]
            or item["logical_tasks"] != binding["logical_fits"]
            or item["priority"] != binding["priority"]
            or item["wrapper_sha256"] != contract["frozen_sources"]["v2_wrapper_sha256"]
        ):
            raise ReconciliationError("v2 readiness, binding, or wrapper differs")
        claim_path = claims_root / f"{item['bundle_id']}.json"
        claimed = claim_path.exists()
        if claimed:
            claim = load_claim(claim_path)
            if (
                claim.get("bundle_id") != item["bundle_id"]
                or claim.get("queue_sha256") != item["_queue_sha256"]
                or claim.get("wrapper_sha256") != item["wrapper_sha256"]
                or claim.get("account", "nslab") != "nslab"
                or claim.get("qos", "nslab") != "nslab"
            ):
                raise ReconciliationError("v2 claim binding differs")
            claim_rows.append({
                "bundle_id": item["bundle_id"],
                "status": claim["status"],
                "job_id": claim.get("job_id", "not_available"),
                "claim_path": str(claim_path.resolve()),
            })
        v2_rows.append({
            "dispatcher_bundle_id": item["bundle_id"],
            "task_bundle_id": binding["task_bundle_id"],
            "priority": item["priority"],
            "logical_fits": item["logical_tasks"],
            "queue_path": str(path.resolve()),
            "queue_sha256": item["_queue_sha256"],
            "enabled": item["enabled"],
            "ready": item["_ready"],
            "claimed": claimed,
        })
    if (
        set(bindings) != v2_ids
        or sum(int(row["logical_fits"]) for row in v2_rows) != 438
        or {str(row["task_bundle_id"]) for row in v2_rows} != v1_ids
    ):
        raise ReconciliationError("v2 dispatcher-to-task mapping differs")

    with (state_root / "dispatcher.lock").open("a+") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_SH)
        gpu_rows = gpu_scheduler_snapshot(user)
        dispatchers = dispatcher_snapshot(user)
        states = [row["state"] for row in gpu_rows]
        submission_allowed = dispatcher.submission_allowed(
            states,
            cap=5,
            max_running=4,
            max_pending=1,
        )
    pending = sum(row["state"] == "PENDING" for row in gpu_rows)
    active = len(gpu_rows) - pending
    if (
        len(gpu_rows) > 5
        or active > 4
        or pending > 1
        or any(row["account"] != "nslab" or row["qos"] != "nslab" for row in gpu_rows)
        or any("masld" in row["job_name"].lower() for row in gpu_rows)
    ):
        raise ReconciliationError("live GPU cap, QOS, or generic-name policy differs")
    central = [
        row for row in dispatchers
        if row["job_name"] == "model-dispatcher" and row["state"] == "RUNNING"
    ]
    if not central:
        raise ReconciliationError("central dispatcher is not running")

    if claim_rows:
        readiness = "v2_dispatch_in_progress"
    elif submission_allowed:
        readiness = "v2_ready_for_next_dispatcher_poll"
    else:
        readiness = "v2_ready_waiting_for_gpu_cap"
    summary = {
        "schema_version": "masld-bench-sequence-five-seed-live-dispatch-readiness-v1",
        "status": "pass",
        "readiness": readiness,
        "v1_queue_items": len(v1_rows),
        "v1_enabled": 0,
        "v1_claimed": 0,
        "v1_ledger_submissions": 0,
        "v2_queue_items": len(v2_rows),
        "v2_enabled_ready": sum(bool(row["enabled"]) and bool(row["ready"]) for row in v2_rows),
        "v2_claimed": len(claim_rows),
        "logical_fits": sum(int(row["logical_fits"]) for row in v2_rows),
        "dispatcher_bundle_identities_distinct": True,
        "dispatcher_bundle_ids_globally_unique": True,
        "live_v2_byte_identical_to_frozen_admission": True,
        "central_dispatcher_running": True,
        "central_dispatcher_job_ids": sorted(row["job_id"] for row in central),
        "gpu_jobs": len(gpu_rows),
        "gpu_running_or_other_active": active,
        "gpu_pending": pending,
        "dispatcher_submission_currently_allowed": submission_allowed,
        "gpu_job_ceiling": 5,
        "maximum_running": 4,
        "maximum_pending": 1,
        "submission_authority": "central_dispatcher_only",
        "manual_gpu_submission_performed": False,
        "account": "nslab",
        "qos": "nslab",
        "innovation_used": False,
        "prediction_or_metric_values_read": False,
        "outcomes_read": False,
        "sealed_data_read": False,
        "cross_model_ranking_performed": False,
    }
    return summary, v1_rows, v2_rows, claim_rows, gpu_rows


def main() -> None:
    args = parse_args()
    root = args.root.resolve(strict=True)
    contract_path = args.contract.resolve(strict=True)
    plan_root = args.plan_root.resolve(strict=True)
    correction_root = args.correction_root.resolve(strict=True)
    queue_admission_root = args.queue_admission_root.resolve(strict=True)
    output = args.output.resolve(strict=False)
    if output.exists():
        raise ReconciliationError("output already exists")
    for path in (contract_path, plan_root, correction_root, queue_admission_root, output):
        if root not in path.parents:
            raise ReconciliationError("input or output escapes benchmark root")
    contract = load_contract(contract_path)
    validate_frozen_roots(
        root=root,
        contract=contract,
        plan_root=plan_root,
        correction_root=correction_root,
        queue_admission_root=queue_admission_root,
    )
    controller = validate_controller(root, contract)
    summary, v1_rows, v2_rows, claim_rows, gpu_rows = reconcile(
        root=root,
        contract=contract,
        queue_admission_root=queue_admission_root,
        user=args.user,
    )
    summary.update(controller)
    summary["contract_sha256"] = digest(contract_path)
    output.mkdir(mode=0o750)
    (output / "contract_snapshot.json").write_bytes(contract_path.read_bytes())
    write_tsv(output / "v1_disposition.tsv", (
        "bundle_id", "queue_path", "queue_sha256", "enabled", "ready", "claim_exists",
    ), v1_rows)
    write_tsv(output / "v2_readiness.tsv", (
        "dispatcher_bundle_id", "task_bundle_id", "priority", "logical_fits",
        "queue_path", "queue_sha256", "enabled", "ready", "claimed",
    ), v2_rows)
    write_tsv(output / "v2_claims.tsv", (
        "bundle_id", "status", "job_id", "claim_path",
    ), claim_rows)
    write_tsv(output / "gpu_scheduler_snapshot.tsv", (
        "job_id", "state", "job_name", "account", "qos",
    ), gpu_rows)
    (output / "summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, sort_keys=True))


if __name__ == "__main__":
    main()
