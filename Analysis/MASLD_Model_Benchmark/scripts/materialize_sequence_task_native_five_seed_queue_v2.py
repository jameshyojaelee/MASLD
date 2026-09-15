#!/usr/bin/env python3
"""Materialize corrected, noncolliding dispatcher identities for the rectangle."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path

from masld_bench.artifacts import verify_frozen_tree


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--plan-root", type=Path, required=True)
    parser.add_argument("--correction-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--priority-start", type=int, default=70)
    return parser.parse_args()


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    args = parse_args()
    root = args.root.resolve(strict=True)
    plan_root = args.plan_root.resolve(strict=True)
    correction_root = args.correction_root.resolve(strict=True)
    output = args.output.resolve(strict=False)
    if (
        root not in plan_root.parents
        or root not in correction_root.parents
        or output.exists()
        or args.priority_start < 1
    ):
        raise SystemExit("unsafe corrected queue materialization path or priority")
    verify_frozen_tree(plan_root)
    verify_frozen_tree(correction_root)
    plan_hash = digest(plan_root / "ARTIFACTS.json")
    correction_hash = digest(correction_root / "ARTIFACTS.json")
    plan_summary = json.loads((plan_root / "plan/summary.json").read_text(encoding="utf-8"))
    correction_metadata = json.loads((correction_root / "ARTIFACTS.json").read_text(encoding="utf-8"))["metadata"]
    if (
        plan_summary.get("status") != "pass"
        or plan_summary.get("manual_gpu_submission_allowed") is not False
        or correction_metadata.get("status") != "passed"
        or correction_metadata.get("v1_queue_items_claimed") != 0
        or correction_metadata.get("v1_queue_items_submitted") != 0
        or correction_metadata.get("scientific_task_contract_changed") is not False
        or correction_metadata.get("manual_gpu_submission_allowed") is not False
    ):
        raise SystemExit("plan or correction safety authority differs")
    with (plan_root / "plan/bundle_allocations.tsv").open(encoding="utf-8", newline="") as handle:
        allocations = list(csv.DictReader(handle, delimiter="\t"))
    if len(allocations) != 30:
        raise SystemExit("corrected queue requires exactly 30 frozen allocations")
    wrapper = root / "slurm/run_sequence_task_native_five_seed_bundle_v2.sbatch"
    wrapper_hash = digest(wrapper)
    if wrapper_hash != "7cfd70d94ad9fb57c4ae6337ed3696e8c4a381bc5d5a96b6800e8ef744c13e3d":
        raise SystemExit("corrected bundle wrapper differs")
    output.mkdir(parents=True)
    bindings: list[dict[str, object]] = []
    for allocation in allocations:
        rank = int(allocation["bundle_rank"])
        dispatcher_id = f"model-training-{800 + rank:03d}"
        task_bundle_id = allocation["bundle_id"]
        item = {
            "schema_version": "masld-bench-gpu-bundle-queue-item-v1",
            "bundle_id": dispatcher_id,
            "priority": args.priority_start + rank - 1,
            "enabled": True,
            "wrapper_path": str(wrapper.relative_to(root)),
            "wrapper_sha256": wrapper_hash,
            "exports": {
                "MASLD_GPU_DISPATCHER_AUTHORITY": "sequence_five_seed_rectangle_v1",
                "RECTANGLE_BUNDLE_ID": task_bundle_id,
                "RECTANGLE_PLAN_ARTIFACTS_SHA256": plan_hash,
                "RECTANGLE_PLAN_ROOT": str(plan_root),
            },
            "required_paths": [
                str((plan_root / "ARTIFACTS.json").relative_to(root)),
                str((correction_root / "ARTIFACTS.json").relative_to(root)),
            ],
            "logical_tasks": int(allocation["logical_fits"]),
            "family": "model_training",
        }
        path = output / f"{item['priority']:03d}-{dispatcher_id}-v2.json"
        path.write_text(json.dumps(item, sort_keys=True, indent=2) + "\n", encoding="utf-8")
        bindings.append({
            "dispatcher_bundle_id": dispatcher_id,
            "task_bundle_id": task_bundle_id,
            "priority": item["priority"],
            "logical_fits": item["logical_tasks"],
        })
    manifest = {
        "schema_version": "masld-bench-sequence-five-seed-corrected-queue-materialization-v1",
        "status": "pass",
        "plan_root": str(plan_root),
        "plan_artifacts_sha256": plan_hash,
        "correction_root": str(correction_root),
        "correction_artifacts_sha256": correction_hash,
        "queue_items": len(bindings),
        "logical_fits": sum(int(row["logical_fits"]) for row in bindings),
        "priority_range": [args.priority_start, args.priority_start + len(bindings) - 1],
        "dispatcher_identity_range": [bindings[0]["dispatcher_bundle_id"], bindings[-1]["dispatcher_bundle_id"]],
        "task_identity_range": [bindings[0]["task_bundle_id"], bindings[-1]["task_bundle_id"]],
        "bindings": bindings,
        "v1_queue_action": "disabled_preserved",
        "submission_authority": "dispatcher_only_no_manual_gpu_sbatch",
        "qos": "nslab",
        "gpu_squeue_ceiling": 5,
        "manual_gpu_sbatch_allowed": False,
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
