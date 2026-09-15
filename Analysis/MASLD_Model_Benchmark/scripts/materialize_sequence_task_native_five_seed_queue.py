#!/usr/bin/env python3
"""Materialize dispatcher queue items from a frozen five-seed rectangle plan."""

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
    output = args.output.resolve(strict=False)
    if root not in plan_root.parents or output.exists() or args.priority_start < 1:
        raise SystemExit("unsafe plan, output, or priority")
    verify_frozen_tree(plan_root)
    plan_hash = digest(plan_root / "ARTIFACTS.json")
    summary = json.loads((plan_root / "plan/summary.json").read_text(encoding="utf-8"))
    if (
        summary.get("status") != "pass"
        or summary.get("manual_gpu_submission_allowed") is not False
        or summary.get("prediction_or_metric_values_read") is not False
        or summary.get("outcome_paths_available") is not False
    ):
        raise SystemExit("frozen plan safety contract differs")
    with (plan_root / "plan/bundle_allocations.tsv").open(encoding="utf-8", newline="") as handle:
        allocations = list(csv.DictReader(handle, delimiter="\t"))
    if len(allocations) != summary.get("bundle_allocations"):
        raise SystemExit("bundle allocation count differs")
    output.mkdir(parents=True)
    for allocation in allocations:
        wrapper = (root / allocation["bundle_wrapper"]).resolve(strict=True)
        if root not in wrapper.parents or digest(wrapper) != allocation["bundle_wrapper_sha256"]:
            raise SystemExit("bundle wrapper differs")
        rank = int(allocation["bundle_rank"])
        item = {
            "schema_version": "masld-bench-gpu-bundle-queue-item-v1",
            "bundle_id": allocation["bundle_id"],
            "priority": args.priority_start + rank - 1,
            "enabled": True,
            "wrapper_path": allocation["bundle_wrapper"],
            "wrapper_sha256": allocation["bundle_wrapper_sha256"],
            "exports": {
                "MASLD_GPU_DISPATCHER_AUTHORITY": "sequence_five_seed_rectangle_v1",
                "RECTANGLE_BUNDLE_ID": allocation["bundle_id"],
                "RECTANGLE_PLAN_ARTIFACTS_SHA256": plan_hash,
                "RECTANGLE_PLAN_ROOT": str(plan_root),
            },
            "required_paths": [
                str((plan_root / "ARTIFACTS.json").relative_to(root)),
            ],
            "logical_tasks": int(allocation["logical_fits"]),
            "family": "model_training",
        }
        path = output / f"{item['priority']:03d}-{item['bundle_id']}.json"
        path.write_text(json.dumps(item, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    manifest = {
        "schema_version": "masld-bench-sequence-five-seed-queue-materialization-v1",
        "status": "pass",
        "plan_root": str(plan_root),
        "plan_artifacts_sha256": plan_hash,
        "queue_items": len(allocations),
        "priority_range": [args.priority_start, args.priority_start + len(allocations) - 1],
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
