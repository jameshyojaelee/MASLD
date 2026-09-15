#!/usr/bin/env python3
"""Audit the complete ChromBPNet two-view rectangle dispatch lane."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import ArtifactError, verify_frozen_tree
from scripts.audit_chrombpnet_full_head_readiness import (
    audit_readiness as audit_head_semantics,
)


SCHEMA_VERSION = "masld-bench-chrombpnet-full-head-rectangle-execution-readiness-v1"
STATUS = "AUTHORIZED_EXISTING_CENTRAL_QUEUE_WAITING_FOR_FULL_RECTANGLE_LICENSE_AND_EXTERNAL_EVALUATION"
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
FOLDS = tuple(range(5))
SEEDS = (20260824, 20260825, 20260826, 20260827, 20260828)
REQUIRED_OUTPUTS = {
    "model/chrombpnet.h5",
    "model/chrombpnet_nobias.h5",
    "predictions/full_model/profile_probabilities.h5",
    "predictions/full_model/regional_counts.tsv",
    "predictions/nobias_model/profile_probabilities.h5",
    "predictions/nobias_model/regional_counts.tsv",
}


class ChromBPNetRectangleError(ValueError):
    """Raised when the full-head rectangle or dispatch lane differs."""


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ChromBPNetRectangleError(f"JSON authority differs: {path}") from error
    if not isinstance(value, dict):
        raise ChromBPNetRectangleError(f"JSON authority is not an object: {path}")
    return value


def _resolve_file(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ChromBPNetRectangleError("authority path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ChromBPNetRectangleError(f"authority is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise ChromBPNetRectangleError(f"authority escapes or is absent: {relative}") from error
    if not resolved.is_file():
        raise ChromBPNetRectangleError(f"authority is not a file: {relative}")
    return resolved


def _resolve_directory(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ChromBPNetRectangleError("artifact path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ChromBPNetRectangleError(f"artifact is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise ChromBPNetRectangleError(f"artifact escapes or is absent: {relative}") from error
    if not resolved.is_dir():
        raise ChromBPNetRectangleError(f"artifact is not a directory: {relative}")
    return resolved


def _require_file(root: Path, record: Mapping[str, Any]) -> Path:
    path = _resolve_file(root, str(record.get("path", "")))
    if _digest(path) != record.get("sha256"):
        raise ChromBPNetRectangleError(f"file authority differs: {path}")
    return path


def _require_frozen(
    root: Path,
    record: Mapping[str, Any],
    artifact_class: str,
) -> tuple[Path, dict[str, Any]]:
    artifact_root = _resolve_directory(root, str(record.get("path", "")))
    try:
        manifest = verify_frozen_tree(artifact_root)
    except ArtifactError as error:
        raise ChromBPNetRectangleError(f"frozen artifact differs: {artifact_root}") from error
    if (
        _digest(artifact_root / "ARTIFACTS.json") != record.get("artifacts_sha256")
        or manifest.get("metadata", {}).get("artifact_class") != artifact_class
        or manifest.get("metadata", {}).get("status") != "passed"
    ):
        raise ChromBPNetRectangleError(f"frozen artifact identity differs: {artifact_root}")
    return artifact_root, dict(manifest)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _assert_false(mapping: Mapping[str, Any], fields: tuple[str, ...]) -> None:
    for field in fields:
        if mapping.get(field) is not False:
            raise ChromBPNetRectangleError(f"fail-closed field opened: {field}")


def _audit_frozen_authorities(
    root: Path, contract: Mapping[str, Any]
) -> tuple[Path, Path, dict[str, Any]]:
    authorities = contract.get("frozen_authorities", {})
    head_contract = _require_file(root, authorities.get("head_semantics_contract", {}))
    _require_file(root, authorities.get("head_semantics_audit", {}))
    _require_file(root, authorities.get("campaign_contract", {}))
    _require_file(root, authorities.get("production_wrapper", {}))
    _require_file(root, authorities.get("bundle_wrapper", {}))
    _require_file(root, authorities.get("bundle_controller", {}))
    _, head_admission = _require_frozen(
        root,
        authorities.get("head_semantics_admission", {}),
        "chrombpnet_full_head_production_readiness",
    )
    plan_root, plan = _require_frozen(
        root,
        authorities.get("rectangle_admission", {}),
        "sequence_task_native_five_seed_rectangle_admission",
    )
    _, evaluator = _require_frozen(
        root,
        authorities.get("evaluator_readiness", {}),
        "sequence_task_native_five_seed_evaluator_readiness",
    )
    _require_frozen(
        root,
        authorities.get("execution_correction", {}),
        "sequence_task_native_five_seed_execution_correction",
    )
    queue_root, queue = _require_frozen(
        root,
        authorities.get("queue_admission", {}),
        "sequence_task_native_five_seed_corrected_queue_admission",
    )
    _, live = _require_frozen(
        root,
        authorities.get("live_dispatch_readiness", {}),
        "sequence_task_native_five_seed_live_dispatch_readiness",
    )

    base = audit_head_semantics(root, head_contract)
    if (
        base.get("status") != "pass_full_head_reconciled_waiting_for_rectangle_and_license"
        or base.get("full_and_nobias_share_one_fit") is not True
        or base.get("primary_biological_view") != "chrombpnet_nobias"
        or base.get("assay_qc_view") != "chrombpnet_full"
        or base.get("chrombpnet_complete_fits_at_freeze") != 19
        or base.get("chrombpnet_missing_fits_at_freeze") != 106
        or base.get("license_authority_contradiction_present") is not True
        or base.get("weight_redistribution_authorized") is not False
        or base.get("prediction_values_read") is not False
        or base.get("metric_values_read") is not False
    ):
        raise ChromBPNetRectangleError("base full-head semantics admission differs")
    if (
        head_admission.get("metadata", {}).get("contract_sha256") != _digest(head_contract)
        or head_admission.get("metadata", {}).get("complete_fits_at_freeze") != 19
        or head_admission.get("metadata", {}).get("missing_fits_at_freeze") != 106
        or head_admission.get("metadata", {}).get("prediction_values_read") is not False
        or plan.get("metadata", {}).get("prediction_or_metric_values_read") is not False
        or evaluator.get("metadata", {}).get("gate_status") != "waiting_for_full_rectangle"
        or evaluator.get("metadata", {}).get("outcome_access_authorized") is not False
        or evaluator.get("metadata", {}).get("partial_cross_model_ranking_authorized") is not False
        or queue.get("metadata", {}).get("submission_authority")
        != "dispatcher_only_no_manual_gpu_sbatch"
        or queue.get("metadata", {}).get("manual_gpu_sbatch_allowed") is not False
        or live.get("metadata", {}).get("submission_authority") != "central_dispatcher_only"
        or live.get("metadata", {}).get("maximum_running") != 4
        or live.get("metadata", {}).get("maximum_pending") != 1
        or live.get("metadata", {}).get("gpu_job_ceiling") != 5
        or live.get("metadata", {}).get("qos") != "nslab"
    ):
        raise ChromBPNetRectangleError("frozen evaluator or dispatcher authority differs")
    return plan_root, queue_root, base


def _audit_rectangle(
    root: Path, plan_root: Path, contract: Mapping[str, Any]
) -> tuple[int, int, dict[str, dict[str, int]], list[dict[str, str]]]:
    rows = [
        row
        for row in _read_tsv(plan_root / "plan/expected_fit_matrix.tsv")
        if row.get("model_id") == "chrombpnet"
    ]
    expected = {(lineage, fold, seed) for lineage in LINEAGES for fold in FOLDS for seed in SEEDS}
    observed = {(row["lineage_id"], int(row["outer_fold"]), int(row["seed"])) for row in rows}
    if (
        len(rows) != 125
        or observed != expected
        or any(row.get("prediction_views") != "chrombpnet_full,chrombpnet_nobias" for row in rows)
        or any(row.get("fit_state") not in {"complete_reusable", "missing"} for row in rows)
    ):
        raise ChromBPNetRectangleError("ChromBPNet prospective rectangle differs")
    complete = [row for row in rows if row["fit_state"] == "complete_reusable"]
    missing = [row for row in rows if row["fit_state"] == "missing"]
    for row in complete:
        artifact_root = Path(row["canonical_artifact_path"])
        try:
            artifact_root.resolve(strict=True).relative_to(root / "executions")
        except (OSError, ValueError) as error:
            raise ChromBPNetRectangleError("canonical ChromBPNet artifact escapes executions") from error
        manifest_path = artifact_root / "ARTIFACTS.json"
        manifest = _load_json(manifest_path)
        metadata = manifest.get("metadata", {})
        inventory = {item.get("path") for item in manifest.get("artifacts", ()) if isinstance(item, dict)}
        if (
            _digest(manifest_path) != row["canonical_artifacts_sha256"]
            or metadata.get("artifact_class") != "chrombpnet_full_depth_training_and_predictions"
            or metadata.get("model_id") != "chrombpnet"
            or metadata.get("dataset_id") != "gse296875"
            or metadata.get("lineage_id") != row["lineage_id"]
            or metadata.get("split_id") != row["split_id"]
            or metadata.get("seed") != int(row["seed"])
            or metadata.get("status") != "passed"
            or metadata.get("primary_biological_score") != "nobias_model"
            or metadata.get("assay_qc_score") != "full_model"
            or metadata.get("benchmark_metrics_calculated") is not False
            or metadata.get("test_outcomes_used") is not False
            or not REQUIRED_OUTPUTS.issubset(inventory)
            or any("metric" in str(value).lower() or "outcome" in str(value).lower() for value in inventory)
        ):
            raise ChromBPNetRectangleError("canonical ChromBPNet manifest differs")

    counts = Counter(row["lineage_id"] for row in complete)
    lineage_counts = {
        lineage: {"complete": counts[lineage], "missing": 25 - counts[lineage]}
        for lineage in LINEAGES
    }
    rectangle = contract.get("rectangle", {})
    if (
        len(complete) != rectangle.get("complete_fits_at_freeze")
        or len(missing) != rectangle.get("missing_fits_at_freeze")
        or 2 * len(complete) != rectangle.get("complete_prediction_views_at_freeze")
        or 2 * len(missing) != rectangle.get("missing_prediction_views_at_freeze")
        or lineage_counts != rectangle.get("lineage_fit_counts")
        or rectangle.get("full_compatible_rectangle") is not False
        or rectangle.get("partial_evaluator_access_allowed") is not False
        or rectangle.get("partial_results_may_be_ranked") is not False
    ):
        raise ChromBPNetRectangleError("ChromBPNet rectangle count contract differs")
    return len(complete), len(missing), lineage_counts, missing


def _audit_dispatch(
    root: Path,
    plan_root: Path,
    queue_root: Path,
    contract: Mapping[str, Any],
    missing_rows: list[dict[str, str]],
) -> None:
    dispatch = contract.get("central_dispatch_authorization", {})
    task_rows = [
        row
        for row in _read_tsv(plan_root / "plan/bundle_tasks.tsv")
        if row.get("model_id") == "chrombpnet"
    ]
    allocations = {
        row["bundle_id"]: row
        for row in _read_tsv(plan_root / "plan/bundle_allocations.tsv")
        if row.get("model_id") == "chrombpnet"
    }
    task_contexts = {(row["lineage_id"], int(row["outer_fold"]), int(row["seed"])) for row in task_rows}
    missing_contexts = {(row["lineage_id"], int(row["outer_fold"]), int(row["seed"])) for row in missing_rows}
    if len(task_rows) != 106 or len(task_contexts) != 106 or task_contexts != missing_contexts:
        raise ChromBPNetRectangleError("ChromBPNet task bundles do not exactly cover missing fits")

    live_root = _resolve_directory(root, str(dispatch.get("queue_live_root", "")))
    wrapper = contract["frozen_authorities"]["bundle_wrapper"]
    plan_sha = contract["frozen_authorities"]["rectangle_admission"]["artifacts_sha256"]
    seen_tasks: set[str] = set()
    seen_dispatchers: set[str] = set()
    queued = 0
    for bundle in dispatch.get("bundles", ()):
        task_id = str(bundle.get("task_bundle_id", ""))
        dispatcher_id = str(bundle.get("dispatcher_bundle_id", ""))
        filename = str(bundle.get("queue_filename", ""))
        frozen_path = queue_root / "queue_items" / filename
        live_path = live_root / filename
        if (
            task_id in seen_tasks
            or dispatcher_id in seen_dispatchers
            or not frozen_path.is_file()
            or not live_path.is_file()
            or frozen_path.read_bytes() != live_path.read_bytes()
            or _digest(frozen_path) != bundle.get("queue_sha256")
        ):
            raise ChromBPNetRectangleError("ChromBPNet queue byte lock differs")
        seen_tasks.add(task_id)
        seen_dispatchers.add(dispatcher_id)
        item = _load_json(live_path)
        tasks = [row for row in task_rows if row["bundle_id"] == task_id]
        allocation = allocations.get(task_id, {})
        if (
            item.get("bundle_id") != dispatcher_id
            or item.get("enabled") is not True
            or item.get("logical_tasks") != bundle.get("logical_fits")
            or item.get("priority") != bundle.get("priority")
            or item.get("wrapper_path") != wrapper["path"]
            or item.get("wrapper_sha256") != wrapper["sha256"]
            or item.get("exports", {}).get("RECTANGLE_BUNDLE_ID") != task_id
            or item.get("exports", {}).get("RECTANGLE_PLAN_ARTIFACTS_SHA256") != plan_sha
            or len(tasks) != bundle.get("logical_fits")
            or allocation.get("submission_authority") != "dispatcher_only_no_manual_gpu_sbatch"
            or allocation.get("partition") != "gpu"
            or allocation.get("qos") != "nslab"
            or any(row.get("task_wrapper") != "slurm/chrombpnet_train_predict_five_seed_v2.sbatch" for row in tasks)
            or any(row.get("outcome_paths_available_to_fit") != "false" for row in tasks)
            or any(row.get("sealed_paths_available_to_fit") != "false" for row in tasks)
        ):
            raise ChromBPNetRectangleError("ChromBPNet queue item differs")
        queued += len(tasks)
    if len(seen_tasks) != 7 or len(seen_dispatchers) != 7 or queued != 106:
        raise ChromBPNetRectangleError("ChromBPNet bundle coverage differs")


def audit_readiness(root: Path, contract_path: Path) -> dict[str, Any]:
    """Verify full-head lane metadata without opening predictions or outcomes."""

    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    try:
        contract_path.relative_to(root)
    except ValueError as error:
        raise ChromBPNetRectangleError("contract escapes benchmark root") from error
    contract = _load_json(contract_path)
    if (
        contract.get("schema_version") != SCHEMA_VERSION
        or contract.get("status") != STATUS
        or contract.get("model_id") != "chrombpnet"
        or contract.get("lane_id") != "chrombpnet_full_head_two_view_five_seed_rectangle"
    ):
        raise ChromBPNetRectangleError("ChromBPNet lane contract identity differs")
    head = contract.get("head_contract", {})
    if (
        head.get("one_fit_two_views") is not True
        or head.get("primary_biological_view") != "chrombpnet_nobias"
        or head.get("assay_qc_view") != "chrombpnet_full"
        or head.get("full_view_is_independent_architecture") is not False
        or head.get("full_and_nobias_are_distinct_from_bpnet") is not True
        or head.get("expected_views_per_fit") != 2
    ):
        raise ChromBPNetRectangleError("ChromBPNet head lane semantics differ")

    plan_root, queue_root, base = _audit_frozen_authorities(root, contract)
    complete, missing, lineage_counts, missing_rows = _audit_rectangle(root, plan_root, contract)
    _audit_dispatch(root, plan_root, queue_root, contract, missing_rows)

    release = contract.get("license_and_release_gate", {})
    dispatch = contract.get("central_dispatch_authorization", {})
    firewall = contract.get("claim_firewall", {})
    if (
        release.get("license_authority_contradiction_present") is not True
        or release.get("weight_redistribution_authorized") is not False
        or release.get("open_champion_eligible_now") is not False
        or release.get("restricted_development_comparator_execution_allowed") is not True
        or dispatch.get("new_queue_registration_needed") is not False
        or dispatch.get("direct_gpu_submission_authorized") is not False
    ):
        raise ChromBPNetRectangleError("license, release, or submission gate differs")
    _assert_false(
        firewall,
        (
            "partial_results_ranked",
            "prediction_values_read",
            "metric_values_read",
            "development_outcomes_read",
            "sealed_features_read",
            "sealed_labels_or_outcomes_read",
            "thresholds_changed",
            "gpu_job_submitted_by_this_authorization",
            "champion_claim_allowed",
            "external_transfer_claim_allowed",
            "universal_claim_allowed",
        ),
    )
    return {
        "schema_version": "masld-bench-chrombpnet-full-head-rectangle-readiness-receipt-v1",
        "status": "pass_chrombpnet_full_head_lane_existing_central_queue_authorized",
        "model_id": "chrombpnet",
        "lane_id": "chrombpnet_full_head_two_view_five_seed_rectangle",
        "one_fit_two_views": True,
        "primary_biological_view": "chrombpnet_nobias",
        "assay_qc_view": "chrombpnet_full",
        "full_view_is_independent_architecture": False,
        "expected_fits": 125,
        "complete_fits_at_freeze": complete,
        "missing_fits_at_freeze": missing,
        "expected_prediction_views": 250,
        "complete_prediction_views_at_freeze": 2 * complete,
        "missing_prediction_views_at_freeze": 2 * missing,
        "lineage_fit_counts": lineage_counts,
        "full_compatible_rectangle": False,
        "partial_evaluator_access_authorized": False,
        "missing_fits_already_centrally_registered": missing,
        "central_task_bundles": 7,
        "central_dispatcher_bundle_ids": [row["dispatcher_bundle_id"] for row in dispatch["bundles"]],
        "new_gpu_queue_registration_needed": False,
        "direct_gpu_submission_authorized": False,
        "gpu_squeue_ceiling": 5,
        "maximum_running_gpu_jobs": 4,
        "maximum_pending_gpu_jobs": 1,
        "qos": "nslab",
        "license_authority_contradiction_present": base["license_authority_contradiction_present"],
        "weight_redistribution_authorized": False,
        "open_champion_eligible_now": False,
        "partial_results_ranked": False,
        "prediction_values_read": False,
        "metric_values_read": False,
        "development_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_labels_or_outcomes_read": False,
        "thresholds_changed": False,
        "gpu_job_submitted": False,
        "contract_sha256": _digest(contract_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = audit_readiness(args.project_root, args.contract)
    rendered = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
