#!/usr/bin/env python3
"""Build a deterministic, dispatcher-only GPU bundle plan from a frozen backlog."""

from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from hashlib import sha256
import json
import math
from pathlib import Path
from typing import Iterable


BACKLOG_STATES = {"backlog_ready", "backlog_waiting_input", "backlog_waiting_model"}
ACTIVE_OR_COMPLETE = {"active_pending", "active_running", "completed"}
SEQUENCE_MODELS = (
    "bpnet",
    "chrombpnet",
    "sequence_cnn_control",
    "sequence_transformer_control",
)
ANCHOR_TASKS = (
    "seq__bpnet__cholangiocyte__d0g0__s20260824",
    "seq__chrombpnet__cholangiocyte__d0g0__s20260824",
    "seq__bpnet__fibroblast__d0g0__s20260824",
    "seq__chrombpnet__fibroblast__d0g0__s20260824",
)
ANCHOR_BUNDLE_002 = (
    "seq__bpnet__t_cell__d0g0__s20260824",
    "seq__chrombpnet__t_cell__d0g0__s20260824",
    "seq__bpnet__fibroblast__d1g1__s20260824",
    "seq__chrombpnet__fibroblast__d1g1__s20260824",
    "seq__bpnet__macrophage__d1g1__s20260824",
    "seq__chrombpnet__macrophage__d1g1__s20260824",
    "seq__bpnet__cholangiocyte__d2g2__s20260824",
    "seq__chrombpnet__cholangiocyte__d2g2__s20260824",
)
ANCHOR_BUNDLE_003 = (
    "seq__bpnet__macrophage__d2g2__s20260824",
    "seq__chrombpnet__macrophage__d2g2__s20260824",
    "seq__bpnet__t_cell__d2g2__s20260824",
    "seq__chrombpnet__t_cell__d2g2__s20260824",
    "seq__bpnet__cholangiocyte__d3g3__s20260824",
    "seq__chrombpnet__cholangiocyte__d3g3__s20260824",
    "seq__bpnet__fibroblast__d3g3__s20260824",
    "seq__chrombpnet__fibroblast__d3g3__s20260824",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--backlog", type=Path, required=True)
    parser.add_argument("--policy", type=Path, required=True)
    parser.add_argument("--bundle-001", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: Iterable[dict[str, object]], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def balanced_chunks(rows: list[dict[str, str]], maximum: int = 8, minimum: int = 4) -> list[list[dict[str, str]]]:
    if not rows:
        return []
    groups = math.ceil(len(rows) / maximum)
    if len(rows) // groups < minimum:
        raise ValueError(f"cannot divide {len(rows)} tasks into {minimum}-{maximum}-task bundles")
    base, remainder = divmod(len(rows), groups)
    sizes = [base + (index < remainder) for index in range(groups)]
    chunks = []
    offset = 0
    for size in sizes:
        chunks.append(rows[offset : offset + size])
        offset += size
    return chunks


def digest_text(values: Iterable[str]) -> str:
    return sha256("\n".join(values).encode("utf-8")).hexdigest()


def expected_artifact(row: dict[str, str]) -> str:
    split_id = row["split_id"]
    seed = row["seed"]
    if row["model_id"] == "scbasset":
        suffix = f"scbasset-{split_id}-seed{seed}"
        if row["stage"] == "predict_valid":
            suffix += "-valid"
    else:
        suffix = f"{row['model_id']}-{row['lineage_id']}-{split_id}-seed{seed}"
    return f"executions/{suffix}-{{SLURM_JOB_ID}}"


def allocation(
    rank: int,
    label: str,
    rows: list[dict[str, str]],
    readiness: str,
    estimated_hours: int,
    requested_hours: int,
    family: str,
) -> dict[str, object]:
    logical_ids = [row["logical_work_id"] for row in rows]
    short_hash = digest_text(logical_ids)[:10]
    bundle_id = f"sequence_gpu_bundle_{rank:03d}_{label}_{short_hash}"
    models = ",".join(sorted({row["model_id"] for row in rows}))
    return {
        "bundle_rank": rank,
        "bundle_id": bundle_id,
        "bundle_family": family,
        "runtime_class": "tf2_chromatin_sequence_l40s",
        "models": models,
        "logical_tasks": len(rows),
        "initial_readiness": readiness,
        "estimated_gpu_hours": estimated_hours,
        "requested_time": f"{requested_hours:02d}:00:00",
        "requested_gpu_hours": requested_hours,
        "gres": "gpu:l40s:1",
        "cpus": 8,
        "memory": "32G",
        "partition": "gpu",
        "qos": "nslab",
        "failure_policy": "continue_after_isolated_task_failure",
        "checkpoint_policy": "atomic_receipt_after_each_logical_task",
        "submission_authority": "dispatcher_only_no_manual_sbatch",
        "wrapper": "slurm/run_sequence_gpu_bundle.sbatch",
        "logical_work_ids_sha256": digest_text(logical_ids),
        "rows": rows,
    }


def main() -> None:
    args = parse_args()
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    backlog = read_tsv(args.backlog)
    old_policy = json.loads(args.policy.read_text(encoding="utf-8"))
    bundle_001 = read_tsv(args.bundle_001)
    by_id = {row["logical_work_id"]: row for row in backlog}
    if len(backlog) != 250 or len(by_id) != 250:
        raise SystemExit("frozen backlog cardinality differs")
    if tuple(row["logical_work_id"] for row in bundle_001) != ANCHOR_TASKS:
        raise SystemExit("bundle 001 anchor differs")
    for task in bundle_001:
        original = by_id.get(task["logical_work_id"])
        if original is None or original["current_state"] != "backlog_ready":
            raise SystemExit("bundle 001 is no longer a frozen ready anchor")

    future = [row for row in backlog if row["current_state"] in BACKLOG_STATES]
    excluded = [row for row in backlog if row["current_state"] in ACTIVE_OR_COMPLETE]
    if len(future) != 240 or len(excluded) != 10:
        raise SystemExit("future/excluded task census differs")
    if set(ANCHOR_TASKS) - {row["logical_work_id"] for row in future}:
        raise SystemExit("bundle 001 anchor is absent from future backlog")

    bundles: list[dict[str, object]] = []
    anchor_specs = (
        (1, "sequence_gpu_bundle_001", ANCHOR_TASKS, 16, 20, "sequence_balanced"),
        (2, "sequence_gpu_bundle_002", ANCHOR_BUNDLE_002, 32, 40, "bpnet_chrombpnet_breadth"),
        (3, "sequence_gpu_bundle_003", ANCHOR_BUNDLE_003, 32, 40, "bpnet_chrombpnet_breadth"),
    )
    consumed: set[str] = set()
    for bundle_rank, bundle_id, logical_ids, estimated, requested, family in anchor_specs:
        anchor_rows = [by_id[logical_id] for logical_id in logical_ids]
        if any(row["current_state"] != "backlog_ready" for row in anchor_rows):
            raise SystemExit(f"{bundle_id} contains a non-ready anchor task")
        anchor_bundle = allocation(
            bundle_rank, "balanced", anchor_rows, "ready", estimated, requested, family
        )
        anchor_bundle["bundle_id"] = bundle_id
        if bundle_rank == 1:
            anchor_bundle["wrapper"] = "slurm/sequence_gpu_bundle_001.sbatch"
        bundles.append(anchor_bundle)
        consumed.update(logical_ids)
    rank = 4
    for model_id in SEQUENCE_MODELS:
        for current_state, readiness in (
            ("backlog_ready", "ready"),
            ("backlog_waiting_input", "waiting_input"),
        ):
            rows = sorted(
                (
                    row for row in future
                    if row["logical_work_id"] not in consumed
                    and row["model_id"] == model_id
                    and row["current_state"] == current_state
                ),
                key=lambda row: int(row["queue_rank"]),
            )
            for chunk in balanced_chunks(rows):
                requested = math.ceil(len(chunk) * 4 * 1.25)
                bundles.append(allocation(
                    rank, f"{model_id}-{readiness}", chunk, readiness,
                    len(chunk) * 4, requested, "sequence_specialist",
                ))
                consumed.update(row["logical_work_id"] for row in chunk)
                rank += 1

    scbasset_rows = [row for row in future if row["model_id"] == "scbasset"]
    scbasset_by_id = {row["logical_work_id"]: row for row in scbasset_rows}
    train_rows = [row for row in scbasset_rows if row["stage"] == "train"]
    for train_state, readiness in (
        ("backlog_ready", "ready_with_internal_dependencies"),
        ("backlog_waiting_input", "waiting_input_with_internal_dependencies"),
    ):
        pairs: list[list[dict[str, str]]] = []
        for train in sorted(
            (row for row in train_rows if row["current_state"] == train_state),
            key=lambda row: (int(row["seed"]), int(row["outer_fold"])),
        ):
            prediction_id = train["logical_work_id"].replace("__train__", "__predict_valid__")
            prediction = scbasset_by_id.get(prediction_id)
            if prediction is None or prediction["prerequisite_logical_id"] != train["logical_work_id"]:
                raise SystemExit("scBasset train/predict pairing differs")
            pairs.append([train, prediction])
        pair_chunks = balanced_chunks(
            [{"pair_index": str(index)} for index in range(len(pairs))],
            maximum=4,
            minimum=2,
        )
        offset = 0
        for pair_chunk in pair_chunks:
            selected_pairs = pairs[offset : offset + len(pair_chunk)]
            offset += len(pair_chunk)
            chunk = [task for pair in selected_pairs for task in pair]
            pair_count = len(selected_pairs)
            bundles.append(allocation(
                rank, f"scbasset-{readiness}", chunk, readiness,
                pair_count * 6, math.ceil(pair_count * 6 * 1.25), "scbasset_train_predict",
            ))
            consumed.update(row["logical_work_id"] for row in chunk)
            rank += 1

    future_ids = {row["logical_work_id"] for row in future}
    if consumed != future_ids:
        missing = sorted(future_ids - consumed)
        extra = sorted(consumed - future_ids)
        raise SystemExit(f"bundle coverage differs; missing={missing}, extra={extra}")
    if any(not 4 <= int(bundle["logical_tasks"]) <= 8 for bundle in bundles):
        raise SystemExit("bundle size differs from 4-8 contract")
    if any(int(bundle["requested_gpu_hours"]) > 48 for bundle in bundles):
        raise SystemExit("bundle exceeds 48-hour allocation cap")

    allocation_fields = [
        "bundle_rank", "bundle_id", "bundle_family", "runtime_class", "models",
        "logical_tasks", "initial_readiness", "estimated_gpu_hours", "requested_time",
        "requested_gpu_hours", "gres", "cpus", "memory", "partition", "qos",
        "failure_policy", "checkpoint_policy", "submission_authority", "wrapper",
        "logical_work_ids_sha256",
    ]
    write_tsv(
        output / "bundle_allocations.tsv",
        ({key: bundle[key] for key in allocation_fields} for bundle in bundles),
        allocation_fields,
    )
    task_fields = [
        "bundle_rank", "bundle_id", "task_order", "logical_work_id", "model_id",
        "stage", "lineage_id", "outer_fold", "split_id", "seed", "current_state",
        "input_job_id", "input_state", "input_artifacts_sha256",
        "prerequisite_logical_id", "script", "export_contract",
        "expected_artifact_template", "failure_scope", "submission_authority",
    ]
    task_output = []
    for bundle in bundles:
        for task_order, row in enumerate(bundle["rows"], 1):
            task_output.append({
                "bundle_rank": bundle["bundle_rank"],
                "bundle_id": bundle["bundle_id"],
                "task_order": task_order,
                "logical_work_id": row["logical_work_id"],
                "model_id": row["model_id"],
                "stage": row["stage"],
                "lineage_id": row["lineage_id"],
                "outer_fold": row["outer_fold"],
                "split_id": row["split_id"],
                "seed": row["seed"],
                "current_state": row["current_state"],
                "input_job_id": row["input_job_id"],
                "input_state": row["input_state"],
                "input_artifacts_sha256": row["input_artifacts_sha256"],
                "prerequisite_logical_id": row["prerequisite_logical_id"],
                "script": row["script"],
                "export_contract": row["export_contract"],
                "expected_artifact_template": expected_artifact(row),
                "failure_scope": "isolated_continue",
                "submission_authority": "dispatcher_only_no_manual_sbatch",
            })
    write_tsv(output / "bundle_tasks.tsv", task_output, task_fields)
    excluded_fields = [
        "logical_work_id", "model_id", "stage", "lineage_id", "outer_fold", "seed",
        "current_state", "active_job_id", "attempt_job_ids", "attempt_state_summary",
        "exclusion_reason",
    ]
    write_tsv(
        output / "excluded_active_or_completed_tasks.tsv",
        ({
            **{key: row[key] for key in excluded_fields if key != "exclusion_reason"},
            "exclusion_reason": "not_future_backlog_at_snapshot",
        } for row in sorted(excluded, key=lambda row: int(row["queue_rank"]))),
        excluded_fields,
    )
    summary = {
        "schema_version": "masld-bench-sequence-gpu-bundle-plan-v2",
        "status": "pass",
        "source_logical_tasks": len(backlog),
        "future_logical_tasks": len(future),
        "excluded_active_or_completed_tasks": len(excluded),
        "bundle_allocations": len(bundles),
        "bundle_task_count_range": [
            min(int(bundle["logical_tasks"]) for bundle in bundles),
            max(int(bundle["logical_tasks"]) for bundle in bundles),
        ],
        "estimated_gpu_hours": sum(int(bundle["estimated_gpu_hours"]) for bundle in bundles),
        "requested_gpu_hours": sum(int(bundle["requested_gpu_hours"]) for bundle in bundles),
        "requested_cpu_hours": sum(
            int(bundle["requested_gpu_hours"]) * int(bundle["cpus"]) for bundle in bundles
        ),
        "readiness_counts": dict(sorted(Counter(str(bundle["initial_readiness"]) for bundle in bundles).items())),
        "model_task_counts": dict(sorted(Counter(row["model_id"] for row in task_output).items())),
        "manual_gpu_sbatch_prohibited": True,
        "gpu_squeue_ceiling": old_policy["gpu_squeue_ceiling"],
        "max_running_gpu_jobs": old_policy["max_running_gpu_jobs"],
        "max_pending_gpu_jobs": old_policy["max_pending_gpu_jobs"],
        "dispatcher_submission_batch_size": 1,
        "current_running_jobs_untouched": old_policy["active_snapshot"]["running"],
        "current_pending_jobs_untouched": old_policy["active_snapshot"]["pending"],
        "priority_bundle_anchors": {
            "sequence_gpu_bundle_001": list(ANCHOR_TASKS),
            "sequence_gpu_bundle_002": list(ANCHOR_BUNDLE_002),
            "sequence_gpu_bundle_003": list(ANCHOR_BUNDLE_003),
        },
    }
    (output / "bundle_plan_summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
