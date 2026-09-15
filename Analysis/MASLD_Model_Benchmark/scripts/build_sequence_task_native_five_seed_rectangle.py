#!/usr/bin/env python3
"""Freeze an outcome-blind five-seed task-native sequence fit census and plan."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from datetime import datetime
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import statistics
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import verify_frozen_tree


SCHEMA = "masld-bench-sequence-task-native-five-seed-rectangle-v1"
MODELS = (
    "bpnet",
    "chrombpnet",
    "sequence_cnn_control",
    "sequence_transformer_control",
)
TASK_WRAPPERS = {
    "bpnet": "slurm/bpnet_train_predict_five_seed_v2.sbatch",
    "chrombpnet": "slurm/chrombpnet_train_predict_five_seed_v2.sbatch",
    "sequence_cnn_control": "slurm/sequence_control_train_cnn_five_seed_v2.sbatch",
    "sequence_transformer_control": "slurm/sequence_control_train_transformer_five_seed_v2.sbatch",
}
REQUIRED_OUTPUTS = {
    "bpnet": {
        "model/bpnet.inference.h5",
        "predictions/profile_probabilities.h5",
        "predictions/regional_counts.tsv",
        "predictions/summary.json",
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
REQUIRED_FALSE = {
    "bpnet": ("benchmark_metrics_calculated", "evaluator_outcomes_exposed"),
    "chrombpnet": ("benchmark_metrics_calculated", "test_outcomes_used"),
    "sequence_cnn_control": ("held_donor_atac_exposed", "test_outcomes_used"),
    "sequence_transformer_control": ("held_donor_atac_exposed", "test_outcomes_used"),
}
EXPECTED_CLASSES = {
    "bpnet": "bpnet_training_prediction_campaign",
    "chrombpnet": "chrombpnet_full_depth_training_and_predictions",
    "sequence_cnn_control": "sequence_control_production",
    "sequence_transformer_control": "sequence_control_production",
}
TERMINAL_JOB = re.compile(r"(?<![0-9])([0-9]{8})(?![0-9])")


class RectangleError(RuntimeError):
    """Raised when the prospective rectangle cannot be included."""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--frozen-root", type=Path, required=True)
    parser.add_argument("--bundle-wrapper", type=Path, required=True)
    return parser.parse_args()


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def write_tsv(path: Path, rows: Iterable[Mapping[str, object]], fields: Sequence[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def load_contract(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if value.get("schema_version") != SCHEMA:
        raise RectangleError("rectangle contract schema differs")
    if tuple(value.get("fit_families", ())) != MODELS:
        raise RectangleError("fit-family roster differs")
    if len(value.get("lineages", ())) != 5 or len(set(value["lineages"])) != 5:
        raise RectangleError("lineage roster differs")
    if value.get("diagonal_outer_folds") != list(range(5)):
        raise RectangleError("diagonal-fold roster differs")
    if value.get("fixed_seeds") != [20260824, 20260825, 20260826, 20260827, 20260828]:
        raise RectangleError("fixed-seed roster differs")
    if value["fit_count_contract"].get("total_fits") != 500:
        raise RectangleError("fit-count contract differs")
    if value["fit_count_contract"].get("total_prediction_views") != 625:
        raise RectangleError("prediction-view contract differs")
    execution = value["execution"]
    if (
        execution.get("submission_authority") != "dispatcher_only_no_manual_gpu_sbatch"
        or execution.get("qos") != "nslab"
        or execution.get("account") != "nslab"
        or execution.get("gres") != "gpu:l40s:1"
        or execution.get("gpu_squeue_ceiling") != 5
        or execution.get("maximum_pending_gpu_jobs") != 1
    ):
        raise RectangleError("GPU routing or ceiling contract differs")
    if any(value["claims"].get(key) is not False for key in (
        "champion_claim_allowed", "universal_claim_allowed", "cross_model_ranking_allowed",
        "post_hoc_unit_selection_allowed",
    )):
        raise RectangleError("outcome-blind claim contract differs")
    return value


def terminal_job_id(path: Path) -> int:
    matches = TERMINAL_JOB.findall(path.name)
    return int(matches[-1]) if matches else 10**12


def canonical_choice_key(row: Mapping[str, Any]) -> tuple[int, str]:
    """Choose duplicates without using predictions, metrics, or runtime."""
    path = Path(str(row["artifact_path"]))
    return terminal_job_id(path), str(path)


def balanced_sizes(total: int, maximum: int, minimum: int) -> list[int]:
    if total < 1:
        return []
    groups = math.ceil(total / maximum)
    base, remainder = divmod(total, groups)
    if base < minimum:
        raise RectangleError(f"cannot bundle {total} fits within {minimum}-{maximum}")
    return [base + (index < remainder) for index in range(groups)]


def artifact_candidate(
    path: Path,
    contract: Mapping[str, Any],
) -> tuple[dict[str, Any] | None, str]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None, "invalid_artifacts_json"
    metadata = payload.get("metadata", {})
    model = metadata.get("model_id")
    if model not in MODELS or metadata.get("dataset_id") != contract["dataset_id"]:
        return None, "outside_model_dataset_scope"
    lineage = metadata.get("lineage_id")
    split_id = metadata.get("split_id")
    seed = metadata.get("seed")
    if lineage not in contract["lineages"] or seed not in contract["fixed_seeds"]:
        return None, "outside_rectangle"
    match = re.fullmatch(r"donor([0-4])_genomic([0-4])", str(split_id))
    if match is None or match.group(1) != match.group(2):
        return None, "not_diagonal_fold"
    fold = int(match.group(1))
    if metadata.get("status") != "passed":
        return None, "status_not_passed"
    if metadata.get("artifact_class") != EXPECTED_CLASSES[model]:
        return None, "artifact_class_differs"
    if any(metadata.get(flag) is not False for flag in REQUIRED_FALSE[model]):
        return None, "explicit_outcome_firewall_flag_missing_or_true"
    if any(metadata.get(flag) is True for flag in contract["artifact_admission"]["forbid_true_metadata_flags"]):
        return None, "forbidden_metadata_flag_true"
    artifact_paths = {str(row.get("path")) for row in payload.get("artifacts", ())}
    if not REQUIRED_OUTPUTS[model].issubset(artifact_paths):
        return None, "required_fit_or_prediction_artifact_missing"
    if any("outcome" in item.lower() or "metric" in item.lower() for item in artifact_paths):
        return None, "outcome_or_metric_artifact_present"
    try:
        verify_frozen_tree(path.parent)
    except Exception:
        return None, "frozen_tree_verification_failed"
    return {
        "model_id": model,
        "lineage_id": lineage,
        "outer_fold": fold,
        "split_id": split_id,
        "seed": int(seed),
        "artifact_path": str(path.parent.resolve()),
        "artifacts_sha256": digest(path),
        "artifact_bytes": sum(int(row["size_bytes"]) for row in payload["artifacts"]),
        "terminal_slurm_job_id": terminal_job_id(path.parent),
    }, "admitted"


def scan_existing(root: Path, contract: Mapping[str, Any]) -> tuple[list[dict[str, Any]], list[dict[str, Any]], list[dict[str, str]]]:
    grouped: dict[tuple[str, str, int, int], list[dict[str, Any]]] = defaultdict(list)
    rejected: list[dict[str, str]] = []
    for path in sorted((root / "executions").glob("*/ARTIFACTS.json")):
        row, reason = artifact_candidate(path, contract)
        if row is None:
            if reason not in {"outside_model_dataset_scope", "outside_rectangle"}:
                rejected.append({"artifact_path": str(path.parent.resolve()), "reason": reason})
            continue
        key = (row["model_id"], row["lineage_id"], row["outer_fold"], row["seed"])
        grouped[key].append(row)
    canonical: list[dict[str, Any]] = []
    duplicates: list[dict[str, Any]] = []
    for key in sorted(grouped):
        values = sorted(grouped[key], key=canonical_choice_key)
        canonical.append(values[0])
        for alias in values[1:]:
            duplicates.append({
                **{name: alias[name] for name in (
                    "model_id", "lineage_id", "outer_fold", "split_id", "seed",
                    "artifact_path", "artifacts_sha256", "terminal_slurm_job_id",
                )},
                "canonical_artifact_path": values[0]["artifact_path"],
                "canonical_rule": "lowest_terminal_slurm_job_id_then_lexicographic_path",
            })
    return canonical, duplicates, rejected


def scan_inputs(root: Path, contract: Mapping[str, Any]) -> list[dict[str, Any]]:
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    for path in sorted((root / "executions").glob("chrombpnet-*-inputs-*/ARTIFACTS.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        metadata = payload.get("metadata", {})
        if metadata.get("artifact_class") != "chrombpnet_training_inputs":
            continue
        if metadata.get("dataset_id") != contract["dataset_id"]:
            continue
        lineage = metadata.get("lineage_id")
        match = re.fullmatch(r"donor([0-4])_genomic([0-4])", str(metadata.get("split_id")))
        if lineage not in contract["lineages"] or match is None or match.group(1) != match.group(2):
            continue
        if (
            metadata.get("status") != "passed"
            or metadata.get("held_atac_used_for_peak_discovery") is not False
            or metadata.get("test_outcomes_available_to_training") is not False
        ):
            raise RectangleError(f"input outcome firewall differs: {path.parent.name}")
        verify_frozen_tree(path.parent)
        fold = int(match.group(1))
        grouped[(str(lineage), fold)].append({
            "lineage_id": lineage,
            "outer_fold": fold,
            "split_id": metadata["split_id"],
            "input_path": str(path.parent.resolve()),
            "input_job_id": terminal_job_id(path.parent),
            "input_artifacts_sha256": digest(path),
        })
    expected = {(lineage, fold) for lineage in contract["lineages"] for fold in range(5)}
    if set(grouped) != expected or any(len(values) != 1 for values in grouped.values()):
        raise RectangleError("exact 25-cell-state/fold input census differs")
    return [grouped[key][0] for key in sorted(grouped)]


def receipt_hours(root: Path, admitted_paths: set[str]) -> dict[str, list[float]]:
    values: dict[str, list[float]] = defaultdict(list)
    for path in sorted((root / "executions").glob("sequence_gpu_bundle*/task_receipts/*.json")):
        try:
            receipt = json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        if (
            receipt.get("status") != "passed"
            or receipt.get("recovered_from_existing_frozen_artifact") is not False
            or receipt.get("expected_artifact") not in admitted_paths
            or receipt.get("model_id") not in MODELS
        ):
            continue
        try:
            started = datetime.fromisoformat(receipt["started_at"].replace("Z", "+00:00"))
            ended = datetime.fromisoformat(receipt["ended_at"].replace("Z", "+00:00"))
        except (KeyError, TypeError, ValueError):
            continue
        hours = (ended - started).total_seconds() / 3600
        if hours > 0:
            values[receipt["model_id"]].append(hours)
    return values


def percentile90(values: Sequence[float]) -> float:
    ordered = sorted(values)
    return ordered[max(0, math.ceil(0.9 * len(ordered)) - 1)]


def main() -> None:
    args = parse_args()
    root = args.root.resolve(strict=True)
    contract_path = args.contract.resolve(strict=True)
    output = args.output.resolve(strict=False)
    frozen_root = args.frozen_root.resolve(strict=False)
    bundle_wrapper = args.bundle_wrapper.resolve(strict=True)
    if output.exists():
        raise RectangleError("output already exists")
    if root not in contract_path.parents or root not in bundle_wrapper.parents:
        raise RectangleError("source path escapes benchmark root")
    if root not in frozen_root.parents:
        raise RectangleError("frozen plan root escapes benchmark root")
    contract = load_contract(contract_path)
    output.mkdir(parents=True)
    (output / "contract_snapshot.json").write_bytes(contract_path.read_bytes())

    canonical, duplicates, rejected = scan_existing(root, contract)
    inputs = scan_inputs(root, contract)
    input_index = {(row["lineage_id"], row["outer_fold"]): row for row in inputs}
    existing_index = {
        (row["model_id"], row["lineage_id"], row["outer_fold"], row["seed"]): row
        for row in canonical
    }
    expected_keys = [
        (model, lineage, fold, seed)
        for model in MODELS
        for lineage in contract["lineages"]
        for fold in contract["diagonal_outer_folds"]
        for seed in contract["fixed_seeds"]
    ]
    if len(expected_keys) != 500 or len(set(expected_keys)) != 500:
        raise RectangleError("expected rectangle is not exactly 500 unique fits")

    matrix: list[dict[str, Any]] = []
    missing: list[dict[str, Any]] = []
    for model, lineage, fold, seed in expected_keys:
        existing = existing_index.get((model, lineage, fold, seed))
        input_row = input_index[(lineage, fold)]
        row = {
            "model_id": model,
            "lineage_id": lineage,
            "outer_fold": fold,
            "split_id": f"donor{fold}_genomic{fold}",
            "seed": seed,
            "prediction_views": ",".join(contract["prediction_views"][model]),
            "fit_state": "complete_reusable" if existing else "missing",
            "canonical_artifact_path": existing["artifact_path"] if existing else "not_available",
            "canonical_artifacts_sha256": existing["artifacts_sha256"] if existing else "not_available",
            "input_path": input_row["input_path"],
            "input_artifacts_sha256": input_row["input_artifacts_sha256"],
        }
        matrix.append(row)
        if existing is None:
            task_wrapper = (root / TASK_WRAPPERS[model]).resolve(strict=True)
            missing.append({
                **row,
                "logical_work_id": f"rect__{model}__{lineage}__d{fold}g{fold}__s{seed}",
                "input_job_id": input_row["input_job_id"],
                "task_wrapper": str(task_wrapper.relative_to(root)),
                "task_wrapper_sha256": digest(task_wrapper),
                "export_contract": (
                    f"OUTER_FOLD={fold};LINEAGE={lineage};SEED={seed};"
                    f"INPUT_JOB_ID={input_row['input_job_id']};"
                    f"INPUT_ARTIFACTS_SHA256={input_row['input_artifacts_sha256']}"
                ),
            })

    admitted_paths = {row["artifact_path"] for row in canonical}
    durations = receipt_hours(root, admitted_paths)
    runtime_rows: list[dict[str, Any]] = []
    runtime_p90: dict[str, float] = {}
    storage_median: dict[str, float] = {}
    for model in MODELS:
        model_durations = durations[model]
        if len(model_durations) < 5:
            raise RectangleError(f"fewer than five observed task runtimes for {model}")
        sizes = [row["artifact_bytes"] for row in canonical if row["model_id"] == model]
        runtime_p90[model] = percentile90(model_durations)
        storage_median[model] = statistics.median(sizes)
        runtime_rows.append({
            "model_id": model,
            "observed_runtime_tasks": len(model_durations),
            "median_gpu_hours": f"{statistics.median(model_durations):.6f}",
            "p90_gpu_hours": f"{runtime_p90[model]:.6f}",
            "maximum_gpu_hours": f"{max(model_durations):.6f}",
            "observed_artifacts": len(sizes),
            "median_storage_bytes": int(storage_median[model]),
            "maximum_storage_bytes": max(sizes),
            "source": "frozen_task_receipt_timestamps_and_artifact_manifests_only",
        })

    maximum = int(contract["execution"]["maximum_logical_fits_per_bundle"])
    minimum = int(contract["execution"]["minimum_logical_fits_per_bundle"])
    allocations: list[dict[str, Any]] = []
    bundled_tasks: list[dict[str, Any]] = []
    rank = 0
    for model in MODELS:
        model_missing = [row for row in missing if row["model_id"] == model]
        sizes = balanced_sizes(len(model_missing), maximum, minimum)
        offset = 0
        for size in sizes:
            rank += 1
            tasks = model_missing[offset : offset + size]
            offset += size
            bundle_id = f"model-training-{700 + rank:03d}"
            logical_ids = [row["logical_work_id"] for row in tasks]
            allocations.append({
                "bundle_rank": rank,
                "bundle_id": bundle_id,
                "model_id": model,
                "logical_fits": len(tasks),
                "logical_work_ids_sha256": sha256("\n".join(logical_ids).encode()).hexdigest(),
                "estimated_p90_gpu_hours": f"{sum(runtime_p90[model] for _ in tasks):.6f}",
                "time_limit": contract["execution"]["time_limit"],
                "partition": "gpu",
                "account": "nslab",
                "qos": "nslab",
                "gres": "gpu:l40s:1",
                "cpus": 8,
                "memory": "32G",
                "failure_policy": "continue_after_isolated_fit_failure",
                "submission_authority": "dispatcher_only_no_manual_gpu_sbatch",
                "bundle_wrapper": str(bundle_wrapper.relative_to(root)),
                "bundle_wrapper_sha256": digest(bundle_wrapper),
            })
            for order, task in enumerate(tasks, 1):
                bundled_tasks.append({
                    "bundle_rank": rank,
                    "bundle_id": bundle_id,
                    "task_order": order,
                    **{key: task[key] for key in (
                        "logical_work_id", "model_id", "lineage_id", "outer_fold", "split_id",
                        "seed", "input_path", "input_job_id", "input_artifacts_sha256",
                        "task_wrapper", "task_wrapper_sha256", "export_contract",
                    )},
                    "expected_artifact_template": (
                        f"executions/{model}-{task['lineage_id']}-{task['split_id']}"
                        f"-seed{task['seed']}-{{SLURM_JOB_ID}}"
                    ),
                    "outcome_paths_available_to_fit": "false",
                    "sealed_paths_available_to_fit": "false",
                })
    if len(bundled_tasks) != len(missing) or not 1 <= len(allocations) <= 40:
        raise RectangleError("missing-fit bundle coverage differs")

    matrix_fields = [
        "model_id", "lineage_id", "outer_fold", "split_id", "seed", "prediction_views",
        "fit_state", "canonical_artifact_path", "canonical_artifacts_sha256", "input_path",
        "input_artifacts_sha256",
    ]
    write_tsv(output / "expected_fit_matrix.tsv", matrix, matrix_fields)
    write_tsv(output / "canonical_existing_fits.tsv", canonical, [
        "model_id", "lineage_id", "outer_fold", "split_id", "seed", "artifact_path",
        "artifacts_sha256", "artifact_bytes", "terminal_slurm_job_id",
    ])
    write_tsv(output / "duplicate_artifact_aliases.tsv", duplicates, [
        "model_id", "lineage_id", "outer_fold", "split_id", "seed", "artifact_path",
        "artifacts_sha256", "terminal_slurm_job_id", "canonical_artifact_path", "canonical_rule",
    ])
    write_tsv(output / "rejected_candidate_artifacts.tsv", rejected, ["artifact_path", "reason"])
    write_tsv(output / "input_artifacts.tsv", inputs, [
        "lineage_id", "outer_fold", "split_id", "input_path", "input_job_id",
        "input_artifacts_sha256",
    ])
    write_tsv(output / "missing_fit_tasks.tsv", missing, [
        *matrix_fields, "logical_work_id", "input_job_id", "task_wrapper",
        "task_wrapper_sha256", "export_contract",
    ])
    write_tsv(output / "runtime_storage_priors.tsv", runtime_rows, list(runtime_rows[0]))
    write_tsv(output / "bundle_allocations.tsv", allocations, list(allocations[0]))
    write_tsv(output / "bundle_tasks.tsv", bundled_tasks, list(bundled_tasks[0]))

    counts = Counter(row["model_id"] for row in canonical)
    missing_counts = Counter(row["model_id"] for row in missing)
    complete_views = sum(len(contract["prediction_views"][row["model_id"]]) for row in canonical)
    missing_views = 625 - complete_views
    estimates = {
        "schema_version": "masld-bench-sequence-five-seed-resource-estimate-v1",
        "basis": "frozen_task_receipt_timestamps_and_artifact_manifest_sizes_only",
        "missing_fit_count": len(missing),
        "missing_prediction_view_count": missing_views,
        "median_estimated_gpu_hours": sum(
            statistics.median(durations[model]) * count for model, count in missing_counts.items()
        ),
        "p90_estimated_gpu_hours": sum(runtime_p90[model] * count for model, count in missing_counts.items()),
        "estimated_storage_bytes": int(sum(storage_median[model] * count for model, count in missing_counts.items())),
        "bundle_allocations": len(allocations),
        "requested_gpu_hours_ceiling": len(allocations) * 48,
        "requested_cpu_hours_ceiling": len(allocations) * 48 * 8,
        "concurrency_policy": "at_most_four_running_plus_one_pending_GPU_jobs_user_wide",
    }
    (output / "resource_estimate.json").write_text(
        json.dumps(estimates, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    summary = {
        "schema_version": "masld-bench-sequence-task-native-five-seed-census-v1",
        "status": "pass",
        "contract_sha256": digest(contract_path),
        "fit_families": list(MODELS),
        "lineages": contract["lineages"],
        "diagonal_outer_folds": contract["diagonal_outer_folds"],
        "fixed_seeds": contract["fixed_seeds"],
        "expected_fits": 500,
        "expected_prediction_views": 625,
        "complete_reusable_fits": len(canonical),
        "complete_reusable_fit_counts": dict(sorted(counts.items())),
        "complete_reusable_prediction_views": complete_views,
        "missing_fits": len(missing),
        "missing_fit_counts": dict(sorted(missing_counts.items())),
        "missing_prediction_views": missing_views,
        "duplicate_aliases": len(duplicates),
        "rejected_candidates": len(rejected),
        "input_artifacts": len(inputs),
        "bundle_allocations": len(allocations),
        "bundle_fit_count_range": [
            min(row["logical_fits"] for row in allocations),
            max(row["logical_fits"] for row in allocations),
        ],
        "full_rectangle_retained": True,
        "chrombpnet_views_share_one_fit": True,
        "prediction_or_metric_values_read": False,
        "outcome_paths_available": False,
        "sealed_paths_available": False,
        "cross_model_ranking_performed": False,
        "post_hoc_unit_selection_performed": False,
        "manual_gpu_submission_allowed": False,
        "frozen_plan_root": str(frozen_root),
        "bundle_wrapper": str(bundle_wrapper.relative_to(root)),
        "bundle_wrapper_sha256": digest(bundle_wrapper),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
