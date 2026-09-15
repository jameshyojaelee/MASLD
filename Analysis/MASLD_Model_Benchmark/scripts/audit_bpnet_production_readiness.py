#!/usr/bin/env python3
"""Audit BPNet identity, frozen rectangle, and central dispatch authorization."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


SCHEMA_VERSION = "masld-bench-bpnet-production-readiness-v1"
STATUS = "AUTHORIZED_EXISTING_CENTRAL_QUEUE_WAITING_FOR_FULL_RECTANGLE_AND_EXTERNAL_EVALUATION"
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
FOLDS = tuple(range(5))
SEEDS = (20260824, 20260825, 20260826, 20260827, 20260828)
LEGACY_SEEDS = (11, 29, 47, 71, 101)
REQUIRED_OUTPUTS = {
    "model/bpnet.inference.h5",
    "predictions/profile_probabilities.h5",
    "predictions/regional_counts.tsv",
}


class BPNetReadinessError(ValueError):
    """Raised when a frozen BPNet authority differs."""


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
        raise BPNetReadinessError(f"JSON authority differs: {path}") from error
    if not isinstance(value, dict):
        raise BPNetReadinessError(f"JSON authority is not an object: {path}")
    return value


def _resolve_file(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise BPNetReadinessError("authority path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise BPNetReadinessError(f"authority is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise BPNetReadinessError(f"authority escapes or is absent: {relative}") from error
    if not resolved.is_file():
        raise BPNetReadinessError(f"authority is not a file: {relative}")
    return resolved


def _resolve_directory(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise BPNetReadinessError("artifact path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise BPNetReadinessError(f"artifact is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise BPNetReadinessError(f"artifact escapes or is absent: {relative}") from error
    if not resolved.is_dir():
        raise BPNetReadinessError(f"artifact is not a directory: {relative}")
    return resolved


def _require_file(root: Path, record: Mapping[str, Any]) -> Path:
    path = _resolve_file(root, str(record.get("path", "")))
    if _digest(path) != record.get("sha256"):
        raise BPNetReadinessError(f"file authority differs: {path}")
    return path


def _require_frozen_artifact(
    root: Path,
    relative: str,
    expected_sha256: str,
    expected_class: str,
) -> tuple[Path, dict[str, Any]]:
    artifact_root = _resolve_directory(root, relative)
    try:
        manifest = verify_frozen_tree(artifact_root)
    except ArtifactError as error:
        raise BPNetReadinessError(f"frozen artifact differs: {relative}") from error
    if (
        _digest(artifact_root / "ARTIFACTS.json") != expected_sha256
        or manifest.get("metadata", {}).get("artifact_class") != expected_class
        or manifest.get("metadata", {}).get("status") != "passed"
    ):
        raise BPNetReadinessError(f"frozen artifact identity differs: {relative}")
    return artifact_root, dict(manifest)


def _require_manifest_only(
    root: Path,
    relative: str,
    expected_sha256: str,
    expected_class: str,
) -> tuple[Path, dict[str, Any]]:
    artifact_root = _resolve_directory(root, relative)
    manifest_path = artifact_root / "ARTIFACTS.json"
    manifest = _load_json(manifest_path)
    if (
        _digest(manifest_path) != expected_sha256
        or manifest.get("metadata", {}).get("artifact_class") != expected_class
    ):
        raise BPNetReadinessError(f"artifact manifest identity differs: {relative}")
    return artifact_root, manifest


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _assert_false(mapping: Mapping[str, Any], fields: tuple[str, ...]) -> None:
    for field in fields:
        if mapping.get(field) is not False:
            raise BPNetReadinessError(f"fail-closed field opened: {field}")


def _audit_identity(root: Path, contract: Mapping[str, Any]) -> None:
    authorities = contract.get("frozen_authorities", {})
    checkpoints = _load_json(_require_file(root, authorities.get("checkpoint_authority", {})))
    crosswalk = _load_json(_require_file(root, authorities.get("development_crosswalk", {})))
    exposure = _load_json(_require_file(root, authorities.get("exposure_audit", {})))
    licenses = [
        _load_json(_require_file(root, authorities.get(name, {})))
        for name in ("code_license", "weights_license", "derivative_weights_license")
    ]
    bpnet_wrapper = _require_file(root, authorities.get("bpnet_production_wrapper", {}))
    chrombpnet_wrapper = _require_file(root, authorities.get("chrombpnet_production_wrapper", {}))
    chrombpnet = _load_json(_require_file(root, authorities.get("chrombpnet_head_readiness", {})))
    _require_file(root, authorities.get("bundle_wrapper", {}))
    _require_file(root, authorities.get("bundle_controller", {}))

    if (
        checkpoints.get("implementation", {}).get("implementation_type") != "train_locally"
        or checkpoints.get("implementation", {}).get("repository_revision")
        != "f4593eedca51741f25b8c5cc0c0d647faf8c8a3c"
        or checkpoints.get("architecture_contract", {}).get("input_length_bp") != 2114
        or checkpoints.get("architecture_contract", {}).get("output_length_bp") != 1000
        or checkpoints.get("input_contract", {}).get("control_policy", "").startswith(
            "No control or bias bigWig"
        )
        is not True
        or checkpoints.get("inference_contract", {}).get("atac_available_at_inference") is not False
        or checkpoints.get("inference_contract", {}).get("bias_inputs")
        != "none_for_the_prespecified_sequence_only_baseline"
        or checkpoints.get("weight_license") != "NOT_APPLICABLE_PROJECT_TRAINED"
    ):
        raise BPNetReadinessError("BPNet checkpoint or architecture identity differs")
    if (
        tuple(crosswalk.get("local_model_roster", {}).get("cell_states", ())) != LINEAGES
        or tuple(crosswalk.get("local_model_roster", {}).get("seeds", ())) != LEGACY_SEEDS
        or crosswalk.get("findings", {}).get("gse296875", {}).get("exposure_state")
        != "continual_seen"
        or exposure.get("checkpoint_findings", {}).get("bpnet", {}).get("upstream_pretraining")
        != "none"
        or exposure.get("checkpoint_findings", {}).get("bpnet", {}).get("exposure_state")
        != "target_label_unexposed"
    ):
        raise BPNetReadinessError("BPNet crosswalk or exposure boundary differs")
    if (
        licenses[0].get("declared_license") != "MIT"
        or any(row.get("use_allowed") is not True for row in licenses)
        or any(row.get("redistribution_allowed") is not True for row in licenses)
        or any(row.get("derivative_redistribution_allowed") is not True for row in licenses)
    ):
        raise BPNetReadinessError("BPNet license authority differs")

    acquisition = authorities.get("source_acquisition", {})
    _, source = _require_manifest_only(
        root,
        str(acquisition.get("path", "")),
        str(acquisition.get("artifacts_sha256", "")),
        "bpnet_source_acquisition",
    )
    if (
        source.get("metadata", {}).get("revision")
        != "f4593eedca51741f25b8c5cc0c0d647faf8c8a3c"
        or source.get("metadata", {}).get("upstream_checkpoint") is not None
    ):
        raise BPNetReadinessError("BPNet source acquisition differs")

    bpnet_text = bpnet_wrapper.read_text(encoding="utf-8")
    chrombpnet_text = chrombpnet_wrapper.read_text(encoding="utf-8")
    if (
        bpnet_wrapper == chrombpnet_wrapper
        or 'training["architecture"]["parameter_count"] != 109124' not in bpnet_text
        or 'training["architecture"]["bias_or_control_input"]' not in bpnet_text
        or '"artifact_class": "bpnet_training_prediction_campaign"' not in bpnet_text
        or '"model_id": "bpnet"' not in bpnet_text
        or 'BIAS_MODEL="${BIAS_ROOT}/models/bias_models/ATAC/' not in chrombpnet_text
        or '"primary_biological_score": "nobias_model"' not in chrombpnet_text
        or '"assay_qc_score": "full_model"' not in chrombpnet_text
    ):
        raise BPNetReadinessError("BPNet and ChromBPNet execution identities collapsed")
    head = chrombpnet.get("scientific_head_contract", {})
    identity = contract.get("scientific_identity", {})
    if (
        identity.get("bpnet_is_chrombpnet_nobias_view") is not False
        or identity.get("chrombpnet_nobias_is_a_bpnet_checkpoint") is not False
        or head.get("one_fit_two_views") is not True
        or head.get("primary_biological_view") != "chrombpnet_nobias"
        or head.get("full_view_is_independent_architecture") is not False
    ):
        raise BPNetReadinessError("BPNet/ChromBPNet scientific separation differs")


def _audit_rectangle(root: Path, contract: Mapping[str, Any]) -> tuple[int, int, dict[str, dict[str, int]]]:
    rectangle = contract.get("five_seed_rectangle", {})
    campaign_path = _resolve_file(root, str(rectangle.get("contract_path", "")))
    if _digest(campaign_path) != rectangle.get("contract_sha256"):
        raise BPNetReadinessError("five-seed campaign contract differs")
    campaign = _load_json(campaign_path)
    if (
        tuple(campaign.get("lineages", ())) != LINEAGES
        or tuple(campaign.get("diagonal_outer_folds", ())) != FOLDS
        or tuple(campaign.get("fixed_seeds", ())) != SEEDS
        or campaign.get("prediction_views", {}).get("bpnet") != ["bpnet"]
        or campaign.get("fit_count_contract", {}).get("fits_per_family") != 125
        or campaign.get("execution", {}).get("submission_authority")
        != "dispatcher_only_no_manual_gpu_sbatch"
    ):
        raise BPNetReadinessError("five-seed campaign semantics differ")
    _assert_false(
        campaign.get("claims", {}),
        ("champion_claim_allowed", "cross_model_ranking_allowed", "universal_claim_allowed"),
    )

    plan_root, plan = _require_frozen_artifact(
        root,
        str(rectangle.get("admission_artifact_root", "")),
        str(rectangle.get("admission_artifacts_sha256", "")),
        "sequence_task_native_five_seed_rectangle_admission",
    )
    _, evaluator = _require_frozen_artifact(
        root,
        str(rectangle.get("evaluator_readiness_artifact_root", "")),
        str(rectangle.get("evaluator_readiness_artifacts_sha256", "")),
        "sequence_task_native_five_seed_evaluator_readiness",
    )
    _assert_false(plan.get("metadata", {}), ("prediction_or_metric_values_read", "outcome_paths_available", "sealed_paths_available"))
    if (
        evaluator.get("metadata", {}).get("gate_status") != "waiting_for_full_rectangle"
        or evaluator.get("metadata", {}).get("full_compatible_rectangle") is not False
        or evaluator.get("metadata", {}).get("outcome_access_authorized") is not False
        or evaluator.get("metadata", {}).get("partial_cross_model_ranking_authorized") is not False
    ):
        raise BPNetReadinessError("evaluator firewall differs")

    rows = [
        row
        for row in _read_tsv(plan_root / "plan/expected_fit_matrix.tsv")
        if row.get("model_id") == "bpnet"
    ]
    expected = {(lineage, fold, seed) for lineage in LINEAGES for fold in FOLDS for seed in SEEDS}
    observed = {(row["lineage_id"], int(row["outer_fold"]), int(row["seed"])) for row in rows}
    if len(rows) != 125 or observed != expected or any(row["prediction_views"] != "bpnet" for row in rows):
        raise BPNetReadinessError("BPNet rectangle is not the exact prospective Cartesian product")
    complete = [row for row in rows if row["fit_state"] == "complete_reusable"]
    missing = [row for row in rows if row["fit_state"] == "missing"]
    if len(complete) + len(missing) != len(rows):
        raise BPNetReadinessError("BPNet fit disposition differs")
    expected_complete = {tuple(value) for value in rectangle.get("complete_contexts", ())}
    actual_complete = {(row["lineage_id"], int(row["outer_fold"]), int(row["seed"])) for row in complete}
    if actual_complete != expected_complete:
        raise BPNetReadinessError("BPNet complete-fit contexts differ")

    for row in complete:
        artifact_root = Path(row["canonical_artifact_path"])
        try:
            artifact_root.resolve(strict=True).relative_to(root / "executions")
        except (OSError, ValueError) as error:
            raise BPNetReadinessError("canonical BPNet artifact escapes executions") from error
        manifest_path = artifact_root / "ARTIFACTS.json"
        manifest = _load_json(manifest_path)
        metadata = manifest.get("metadata", {})
        inventory = {item.get("path") for item in manifest.get("artifacts", ()) if isinstance(item, dict)}
        if (
            _digest(manifest_path) != row["canonical_artifacts_sha256"]
            or metadata.get("artifact_class") != "bpnet_training_prediction_campaign"
            or metadata.get("model_id") != "bpnet"
            or metadata.get("dataset_id") != "gse296875"
            or metadata.get("lineage_id") != row["lineage_id"]
            or metadata.get("split_id") != row["split_id"]
            or metadata.get("seed") != int(row["seed"])
            or metadata.get("status") != "passed"
            or metadata.get("bias_or_control_input") is not False
            or metadata.get("benchmark_metrics_calculated") is not False
            or metadata.get("evaluator_outcomes_exposed") is not False
            or not REQUIRED_OUTPUTS.issubset(inventory)
            or any("metric" in str(value).lower() or "outcome" in str(value).lower() for value in inventory)
        ):
            raise BPNetReadinessError("canonical BPNet artifact manifest differs")

    counts = Counter(row["lineage_id"] for row in complete)
    lineage_counts = {
        lineage: {"complete": counts[lineage], "missing": 25 - counts[lineage]}
        for lineage in LINEAGES
    }
    if (
        len(complete) != rectangle.get("complete_fits_at_freeze")
        or len(missing) != rectangle.get("missing_fits_at_freeze")
        or lineage_counts != rectangle.get("lineage_fit_counts")
    ):
        raise BPNetReadinessError("BPNet rectangle count contract differs")
    return len(complete), len(missing), lineage_counts


def _audit_dispatch(root: Path, contract: Mapping[str, Any], missing_count: int) -> None:
    dispatch = contract.get("central_dispatch_authorization", {})
    queue_root, queue_admission = _require_frozen_artifact(
        root,
        str(dispatch.get("queue_admission_artifact_root", "")),
        str(dispatch.get("queue_admission_artifacts_sha256", "")),
        "sequence_task_native_five_seed_corrected_queue_admission",
    )
    _, correction = _require_frozen_artifact(
        root,
        str(dispatch.get("execution_correction_artifact_root", "")),
        str(dispatch.get("execution_correction_artifacts_sha256", "")),
        "sequence_task_native_five_seed_execution_correction",
    )
    _, live = _require_frozen_artifact(
        root,
        str(dispatch.get("live_readiness_artifact_root", "")),
        str(dispatch.get("live_readiness_artifacts_sha256", "")),
        "sequence_task_native_five_seed_live_dispatch_readiness",
    )
    queue_metadata = queue_admission.get("metadata", {})
    if (
        queue_metadata.get("submission_authority") != "dispatcher_only_no_manual_gpu_sbatch"
        or queue_metadata.get("qos") != "nslab"
        or queue_metadata.get("manual_gpu_sbatch_allowed") is not False
        or correction.get("metadata", {}).get("manual_gpu_submission_allowed") is not False
        or live.get("metadata", {}).get("maximum_running") != 4
        or live.get("metadata", {}).get("maximum_pending") != 1
        or live.get("metadata", {}).get("gpu_job_ceiling") != 5
        or live.get("metadata", {}).get("submission_authority") != "central_dispatcher_only"
        or live.get("metadata", {}).get("qos") != "nslab"
    ):
        raise BPNetReadinessError("central GPU governance differs")
    _assert_false(queue_metadata, ("outcome_paths_available", "sealed_paths_available", "prediction_or_metric_values_read"))

    rectangle = contract["five_seed_rectangle"]
    plan_root = _resolve_directory(root, str(rectangle["admission_artifact_root"]))
    task_rows = [
        row
        for row in _read_tsv(plan_root / "plan/bundle_tasks.tsv")
        if row.get("model_id") == "bpnet"
    ]
    allocation_rows = {
        row["bundle_id"]: row
        for row in _read_tsv(plan_root / "plan/bundle_allocations.tsv")
        if row.get("model_id") == "bpnet"
    }
    matrix_rows = [
        row
        for row in _read_tsv(plan_root / "plan/expected_fit_matrix.tsv")
        if row.get("model_id") == "bpnet" and row.get("fit_state") == "missing"
    ]
    task_contexts = {(row["lineage_id"], int(row["outer_fold"]), int(row["seed"])) for row in task_rows}
    missing_contexts = {(row["lineage_id"], int(row["outer_fold"]), int(row["seed"])) for row in matrix_rows}
    if len(task_rows) != missing_count or len(task_contexts) != missing_count or task_contexts != missing_contexts:
        raise BPNetReadinessError("central BPNet tasks do not exactly cover missing fits")

    live_root = _resolve_directory(root, str(dispatch.get("queue_live_root", "")))
    seen_tasks: set[str] = set()
    seen_dispatchers: set[str] = set()
    queued_fits = 0
    wrapper_record = contract["frozen_authorities"]["bundle_wrapper"]
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
            raise BPNetReadinessError("BPNet queue identity or byte lock differs")
        seen_tasks.add(task_id)
        seen_dispatchers.add(dispatcher_id)
        item = _load_json(live_path)
        tasks = [row for row in task_rows if row["bundle_id"] == task_id]
        allocation = allocation_rows.get(task_id, {})
        if (
            item.get("bundle_id") != dispatcher_id
            or item.get("enabled") is not True
            or item.get("logical_tasks") != bundle.get("logical_fits")
            or item.get("priority") != bundle.get("priority")
            or item.get("wrapper_path") != wrapper_record["path"]
            or item.get("wrapper_sha256") != wrapper_record["sha256"]
            or item.get("exports", {}).get("RECTANGLE_BUNDLE_ID") != task_id
            or item.get("exports", {}).get("RECTANGLE_PLAN_ARTIFACTS_SHA256")
            != rectangle["admission_artifacts_sha256"]
            or len(tasks) != bundle.get("logical_fits")
            or allocation.get("submission_authority") != "dispatcher_only_no_manual_gpu_sbatch"
            or allocation.get("partition") != "gpu"
            or allocation.get("qos") != "nslab"
            or any(row.get("task_wrapper") != "slurm/bpnet_train_predict_five_seed_v2.sbatch" for row in tasks)
            or any(row.get("outcome_paths_available_to_fit") != "false" for row in tasks)
            or any(row.get("sealed_paths_available_to_fit") != "false" for row in tasks)
        ):
            raise BPNetReadinessError("BPNet central queue item differs")
        queued_fits += len(tasks)
    if len(seen_tasks) != 7 or len(seen_dispatchers) != 7 or queued_fits != missing_count:
        raise BPNetReadinessError("BPNet bundle authorization count differs")


def audit_readiness(root: Path, contract_path: Path) -> dict[str, Any]:
    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    try:
        contract_path.relative_to(root)
    except ValueError as error:
        raise BPNetReadinessError("contract escapes project root") from error
    contract = _load_json(contract_path)
    if (
        contract.get("schema_version") != SCHEMA_VERSION
        or contract.get("status") != STATUS
        or contract.get("model_id") != "bpnet"
    ):
        raise BPNetReadinessError("BPNet readiness contract identity differs")
    _audit_identity(root, contract)
    complete, missing, lineage_counts = _audit_rectangle(root, contract)
    _audit_dispatch(root, contract, missing)

    release = contract.get("license_and_release_gate", {})
    firewall = contract.get("claim_firewall", {})
    dispatch = contract.get("central_dispatch_authorization", {})
    if (
        release.get("license_blocks_open_champion") is not False
        or release.get("open_champion_eligible_now") is not False
        or dispatch.get("new_queue_registration_needed") is not False
        or dispatch.get("direct_gpu_submission_authorized") is not False
    ):
        raise BPNetReadinessError("release or submission gate differs")
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
        "schema_version": "masld-bench-bpnet-production-readiness-receipt-v1",
        "status": "pass_bpnet_reconciled_existing_central_queue_authorized",
        "model_id": "bpnet",
        "bpnet_is_separate_architecture_from_chrombpnet": True,
        "bpnet_is_chrombpnet_nobias_view": False,
        "prediction_views_per_bpnet_fit": 1,
        "expected_fits": 125,
        "complete_fits_at_freeze": complete,
        "missing_fits_at_freeze": missing,
        "lineage_fit_counts": lineage_counts,
        "full_compatible_rectangle": False,
        "partial_evaluator_access_authorized": False,
        "missing_fits_already_centrally_registered": missing,
        "central_task_bundles": 7,
        "central_dispatcher_bundle_ids": [
            row["dispatcher_bundle_id"] for row in dispatch["bundles"]
        ],
        "new_gpu_queue_registration_needed": False,
        "direct_gpu_submission_authorized": False,
        "gpu_squeue_ceiling": 5,
        "maximum_running_gpu_jobs": 4,
        "maximum_pending_gpu_jobs": 1,
        "qos": "nslab",
        "code_and_project_weight_redistribution_authorized": True,
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
    text = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(text, encoding="utf-8")
    print(text, end="")


if __name__ == "__main__":
    main()
