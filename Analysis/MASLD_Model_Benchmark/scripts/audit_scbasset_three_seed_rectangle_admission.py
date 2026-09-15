#!/usr/bin/env python3
"""Audit the outcome-blind scBasset three-seed production inclusion."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any, Mapping

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


SCHEMA_VERSION = "masld-bench-scbasset-three-seed-rectangle-v1"
STATUS = "PROSPECTIVE_DEVELOPMENT_RECTANGLE_WAITING_FOR_CENTRAL_DISPATCH"
FOLDS = tuple(range(5))
SEEDS = (20260824, 20260825, 20260826)
LEGACY_SEEDS = (11, 29, 47, 71, 101)
COMPLETE_CONTEXTS = {
    (1, 20260824),
    (2, 20260824),
    (3, 20260824),
    (4, 20260824),
}
EXPECTED_CONTEXTS = {(fold, seed) for fold in FOLDS for seed in SEEDS}
MISSING_CONTEXTS = EXPECTED_CONTEXTS - COMPLETE_CONTEXTS


class ScBassetRectangleAdmissionError(ValueError):
    """Raised when a production authority or fail-closed check differs."""


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
        raise ScBassetRectangleAdmissionError(f"JSON authority differs: {path}") from error
    if not isinstance(value, dict):
        raise ScBassetRectangleAdmissionError(f"JSON authority is not an object: {path}")
    return value


def _resolve_file(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ScBassetRectangleAdmissionError("authority path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ScBassetRectangleAdmissionError(f"authority is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise ScBassetRectangleAdmissionError(
            f"authority escapes or is absent: {relative}"
        ) from error
    if not resolved.is_file():
        raise ScBassetRectangleAdmissionError(f"authority is not a file: {relative}")
    return resolved


def _resolve_directory(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ScBassetRectangleAdmissionError("artifact path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ScBassetRectangleAdmissionError(f"artifact is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise ScBassetRectangleAdmissionError(
            f"artifact escapes or is absent: {relative}"
        ) from error
    if not resolved.is_dir():
        raise ScBassetRectangleAdmissionError(f"artifact is not a directory: {relative}")
    return resolved


def _manifest_only(
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
        raise ScBassetRectangleAdmissionError(f"artifact manifest differs: {relative}")
    return artifact_root, manifest


def _frozen_artifact(
    root: Path,
    relative: str,
    expected_sha256: str,
    expected_class: str,
) -> tuple[Path, dict[str, Any]]:
    artifact_root = _resolve_directory(root, relative)
    try:
        manifest = verify_frozen_tree(artifact_root)
    except ArtifactError as error:
        raise ScBassetRectangleAdmissionError(f"frozen artifact differs: {relative}") from error
    if (
        _digest(artifact_root / "ARTIFACTS.json") != expected_sha256
        or manifest.get("metadata", {}).get("artifact_class") != expected_class
    ):
        raise ScBassetRectangleAdmissionError(f"frozen identity differs: {relative}")
    return artifact_root, dict(manifest)


def _inventory(manifest: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    rows = manifest.get("artifacts", [])
    if not isinstance(rows, list):
        raise ScBassetRectangleAdmissionError("artifact inventory is not a list")
    indexed = {str(row.get("path")): dict(row) for row in rows if isinstance(row, dict)}
    if len(indexed) != len(rows):
        raise ScBassetRectangleAdmissionError("artifact inventory has duplicate or invalid rows")
    return indexed


def _assert_false(mapping: Mapping[str, Any], fields: tuple[str, ...]) -> None:
    for field in fields:
        if mapping.get(field) is not False:
            raise ScBassetRectangleAdmissionError(f"fail-closed field opened: {field}")


def _audit_geometry(contract: Mapping[str, Any]) -> None:
    complete = {tuple(row) for row in contract.get("complete_contexts", [])}
    missing = {tuple(row) for row in contract.get("missing_contexts", [])}
    if (
        contract.get("model_id") != "scbasset"
        or contract.get("dataset_id") != "gse296875"
        or contract.get("task_role")
        != "sequence_only_profile_baseline_for_rna_conditioned_atac_development"
        or contract.get("native_observed_atac_role") != "outer_training_supervision_only"
        or contract.get("observed_query_atac_consumed") is not False
        or contract.get("rna_context_consumed") is not False
        or contract.get("biological_outer_unit") != "donor"
        or contract.get("genomic_outer_unit") != "whole_chromosome_group"
        or tuple(contract.get("outer_folds", ())) != FOLDS
        or tuple(contract.get("fixed_screening_seeds", ())) != SEEDS
        or tuple(contract.get("legacy_seeds_excluded", ())) != LEGACY_SEEDS
        or contract.get("expected_fits") != 15
        or contract.get("expected_valid_prediction_views") != 15
        or contract.get("complete_compatible_fits_at_freeze") != 4
        or contract.get("missing_fits_at_freeze") != 11
        or complete != COMPLETE_CONTEXTS
        or missing != MISSING_CONTEXTS
        or complete & missing
        or complete | missing != EXPECTED_CONTEXTS
        or any(seed in LEGACY_SEEDS for _, seed in complete | missing)
    ):
        raise ScBassetRectangleAdmissionError("prospective rectangle geometry differs")


def _audit_readiness(root: Path, contract: Mapping[str, Any]) -> None:
    authority = contract.get("readiness_authority", {})
    _, manifest = _frozen_artifact(
        root,
        str(authority.get("artifact_root", "")),
        str(authority.get("artifacts_sha256", "")),
        str(authority.get("artifact_class", "")),
    )
    receipt = _load_json(
        _resolve_file(root, f"{authority.get('artifact_root', '')}/readiness/receipt.json")
    )
    metadata = manifest.get("metadata", {})
    if (
        metadata.get("status") != "passed"
        or metadata.get("fixed_seed_rectangle_complete_fits") != 4
        or metadata.get("fixed_seed_rectangle_missing_fits") != 11
        or metadata.get("prediction_values_read") is not False
        or metadata.get("metric_values_read") is not False
        or metadata.get("sealed_labels_or_outcomes_read") is not False
        or receipt.get("status") != authority.get("previous_status")
        or receipt.get("valid_next_centrally_dispatched_gpu_unit_exists") is not False
        or receipt.get("partial_ranking_authorized") is not False
        or authority.get("supersession_scope")
        != "exact_11_fit_fixed_seed_development_rectangle_only"
        or authority.get("observed_multiome_eligibility_changed") is not False
        or authority.get("rna_conditioned_model_status_changed") is not False
        or authority.get("champion_eligibility_changed") is not False
    ):
        raise ScBassetRectangleAdmissionError("readiness authority or supersession differs")


def _audit_sources_and_inputs(root: Path, contract: Mapping[str, Any]) -> None:
    sources = contract.get("source_authorities", {})
    for record in sources.values():
        path = _resolve_file(root, str(record.get("path", "")))
        if _digest(path) != record.get("sha256"):
            raise ScBassetRectangleAdmissionError(f"source differs: {path}")
    training_wrapper = _resolve_file(root, sources["training_wrapper"]["path"])
    prediction_wrapper = _resolve_file(root, sources["prediction_wrapper"]["path"])
    for path in (training_wrapper, prediction_wrapper):
        text = path.read_text(encoding="utf-8")
        if (
            ': "${INPUT_ROOT_OVERRIDE:=}"' not in text
            or '*) printf \'Invalid INPUT_ROOT_OVERRIDE\\n\'' not in text
            or "20260824|20260825|20260826" not in text
            or "--qos=innovation" in text
        ):
            raise ScBassetRectangleAdmissionError(f"input override or seed gate differs: {path}")

    records = contract.get("input_artifacts", [])
    if not isinstance(records, list) or len(records) != 5:
        raise ScBassetRectangleAdmissionError("input artifact census differs")
    observed = set()
    for record in records:
        fold = int(record.get("outer_fold", -1))
        split_id = str(record.get("split_id", ""))
        artifact_root, manifest = _manifest_only(
            root,
            str(record.get("artifact_root", "")),
            str(record.get("artifacts_sha256", "")),
            "scbasset_training_inputs",
        )
        metadata = manifest.get("metadata", {})
        inventory = _inventory(manifest)
        summary_path = artifact_root / "inputs/summary.json"
        if _digest(summary_path) != inventory.get("inputs/summary.json", {}).get("sha256"):
            raise ScBassetRectangleAdmissionError(f"input summary differs: fold {fold}")
        summary = _load_json(summary_path)
        if (
            fold not in FOLDS
            or split_id != f"donor{fold}_genomic{fold}"
            or metadata.get("status") != "passed"
            or metadata.get("model_id") != "scbasset"
            or metadata.get("dataset_id") != "gse296875"
            or metadata.get("split_id") != split_id
            or metadata.get("held_donor_atac_used") is not False
            or metadata.get("genomic_test_atac_used") is not False
            or summary.get("donor_test_fold") != fold
            or summary.get("donor_valid_fold") != (fold + 1) % 5
            or summary.get("donor_train_folds")
            != sorted(set(FOLDS).difference({fold, (fold + 1) % 5}))
            or summary.get("held_donor_atac_read") is not False
            or summary.get("genomic_test_atac_used_for_matrix") is not False
            or summary.get("missing_evidence_encoded_as_zero") is not False
        ):
            raise ScBassetRectangleAdmissionError(f"input split isolation differs: fold {fold}")
        observed.add(fold)
    if observed != set(FOLDS):
        raise ScBassetRectangleAdmissionError("input folds differ")


def _audit_complete_artifacts(root: Path, contract: Mapping[str, Any]) -> None:
    records = contract.get("complete_artifacts", [])
    if not isinstance(records, list) or len(records) != 4:
        raise ScBassetRectangleAdmissionError("complete artifact census differs")
    observed = set()
    for record in records:
        fold = int(record.get("outer_fold", -1))
        seed = int(record.get("seed", -1))
        _, model = _manifest_only(
            root,
            str(record.get("model_root", "")),
            str(record.get("model_artifacts_sha256", "")),
            "scbasset_fold_model",
        )
        _, prediction = _manifest_only(
            root,
            str(record.get("prediction_root", "")),
            str(record.get("prediction_artifacts_sha256", "")),
            "scbasset_fold_predictions",
        )
        model_meta = model.get("metadata", {})
        prediction_meta = prediction.get("metadata", {})
        split_id = f"donor{fold}_genomic{fold}"
        if (
            (fold, seed) not in COMPLETE_CONTEXTS
            or model_meta.get("status") != "passed"
            or model_meta.get("model_id") != "scbasset"
            or model_meta.get("split_id") != split_id
            or model_meta.get("seed") != seed
            or model_meta.get("held_donor_atac_used") is not False
            or model_meta.get("genomic_test_atac_used") is not False
            or prediction_meta.get("status") != "passed"
            or prediction_meta.get("model_id") != "scbasset"
            or prediction_meta.get("split_id") != split_id
            or prediction_meta.get("seed") != seed
            or prediction_meta.get("evaluation_role") != "valid"
            or prediction_meta.get("held_donor_atac_used") is not False
            or prediction_meta.get("held_cell_embedding_available") is not False
        ):
            raise ScBassetRectangleAdmissionError(f"complete fit differs: fold {fold}")
        observed.add((fold, seed))
    if observed != COMPLETE_CONTEXTS:
        raise ScBassetRectangleAdmissionError("complete contexts differ")

    filesystem_complete = set()
    pattern = re.compile(
        r"^scbasset-donor([0-4])_genomic\1-seed(20260824|20260825|20260826)-([0-9]+)$"
    )
    for path in sorted((root / "executions").iterdir()):
        match = pattern.fullmatch(path.name)
        if match is None or not (path / "ARTIFACTS.json").is_file():
            continue
        manifest = _load_json(path / "ARTIFACTS.json")
        metadata = manifest.get("metadata", {})
        prediction = path.with_name(
            f"scbasset-donor{match.group(1)}_genomic{match.group(1)}-"
            f"seed{match.group(2)}-valid-{match.group(3)}"
        )
        if (
            metadata.get("artifact_class") == "scbasset_fold_model"
            and metadata.get("status") == "passed"
            and (prediction / "ARTIFACTS.json").is_file()
            and _load_json(prediction / "ARTIFACTS.json").get("metadata", {}).get("status")
            == "passed"
        ):
            filesystem_complete.add((int(match.group(1)), int(match.group(2))))
    if filesystem_complete != COMPLETE_CONTEXTS:
        raise ScBassetRectangleAdmissionError("filesystem rectangle changed after freeze")


def _audit_runtime_and_bundle(root: Path, contract: Mapping[str, Any]) -> None:
    runtime = contract.get("runtime_authorities", {})
    for prefix, artifact_class in (
        ("environment", "chrombpnet_native_runtime_acquisition"),
        ("core_overlay", "chrombpnet_native_core_overlay"),
        ("source", "scbasset_source_acquisition"),
        ("runtime_probe", "scbasset_native_runtime_probe"),
        ("training_probe", "scbasset_training_script_probe"),
        ("prediction_probe", "scbasset_prediction_script_probe"),
    ):
        _, manifest = _manifest_only(
            root,
            str(runtime.get(f"{prefix}_root", "")),
            str(runtime.get(f"{prefix}_artifacts_sha256", "")),
            artifact_class,
        )
        if manifest.get("metadata", {}).get("status") not in {None, "passed"}:
            raise ScBassetRectangleAdmissionError(f"runtime authority failed: {prefix}")

    bundle = contract.get("production_bundle", {})
    if (
        bundle.get("bundle_id") != "model-training-604"
        or bundle.get("queue_priority") != 63
        or bundle.get("job_name") != "model-training"
        or bundle.get("partition") != "gpu"
        or bundle.get("account") != "nslab"
        or bundle.get("qos") != "nslab"
        or bundle.get("gres") != "gpu:l40s:1"
        or bundle.get("cpus") != 8
        or bundle.get("memory") != "64G"
        or bundle.get("time") != "08:00:00"
        or bundle.get("array") is not False
        or bundle.get("fits") != 11
        or bundle.get("logical_tasks") != 22
        or bundle.get("submission_authority") != "central_dispatcher_only"
        or bundle.get("direct_gpu_submission_authorized") is not False
        or bundle.get("manual_gpu_sbatch_allowed") is not False
        or bundle.get("gpu_job_ceiling") != 5
        or bundle.get("maximum_running_gpu_jobs") != 4
        or bundle.get("maximum_pending_gpu_jobs") != 1
    ):
        raise ScBassetRectangleAdmissionError("bundle resource or dispatch contract differs")
    wrapper = _resolve_file(root, str(bundle.get("wrapper_path", "")))
    text = wrapper.read_text(encoding="utf-8")
    required = (
        "#SBATCH --job-name=model-training",
        "#SBATCH --partition=gpu",
        "#SBATCH --account=nslab",
        "#SBATCH --qos=nslab",
        "#SBATCH --gres=gpu:l40s:1",
        "#SBATCH --cpus-per-task=8",
        "#SBATCH --mem=64G",
        "#SBATCH --time=08:00:00",
        "SCBASSET_RECTANGLE_AUTHORITY",
        "SCBASSET_PREFLIGHT_ONLY",
        "INPUT_ROOT_OVERRIDE",
        "continue_after_isolated_fit_failure",
    )
    if (
        any(token not in text for token in required)
        or "#SBATCH --array" in text
        or "--qos=innovation" in text
        or re.search(r"^\s*sbatch\s", text, flags=re.MULTILINE)
    ):
        raise ScBassetRectangleAdmissionError("production wrapper semantics differ")


def _audit_queue(root: Path, contract: Mapping[str, Any]) -> None:
    bundle = contract["production_bundle"]
    queue_path = _resolve_file(
        root, "config/campaigns/gpu_bundle_queue/063-model-training-604.json"
    )
    queue = _load_json(queue_path)
    wrapper = _resolve_file(root, str(bundle["wrapper_path"]))
    contract_path = _resolve_file(
        root, "config/campaigns/scbasset_three_seed_rectangle_20260825.json"
    )
    required_paths = {
        "config/campaigns/scbasset_three_seed_rectangle_20260825.json",
        f"{bundle['admission_artifact_root']}/ARTIFACTS.json",
        f"{contract['readiness_authority']['artifact_root']}/ARTIFACTS.json",
        *(f"{record['artifact_root']}/ARTIFACTS.json" for record in contract["input_artifacts"]),
    }
    if (
        set(queue) != {
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
        or queue.get("schema_version") != "masld-bench-gpu-bundle-queue-item-v1"
        or queue.get("bundle_id") != bundle.get("bundle_id")
        or queue.get("priority") != bundle.get("queue_priority")
        or queue.get("enabled") is not True
        or queue.get("wrapper_path") != bundle.get("wrapper_path")
        or queue.get("wrapper_sha256") != _digest(wrapper)
        or queue.get("logical_tasks") != 22
        or queue.get("family") != "single_cell_chromatin"
        or set(queue.get("required_paths", [])) != required_paths
        or queue.get("exports", {}).get("SCBASSET_RECTANGLE_AUTHORITY")
        != "scbasset_three_seed_rectangle_v1"
        or queue.get("exports", {}).get("SCBASSET_RECTANGLE_CONTRACT")
        != str(contract_path)
        or queue.get("exports", {}).get("SCBASSET_RECTANGLE_CONTRACT_SHA256")
        != _digest(contract_path)
    ):
        raise ScBassetRectangleAdmissionError("central queue item differs")
    gate = root / f"{bundle['admission_artifact_root']}/ARTIFACTS.json"
    claim = root / "executions/gpu-bundle-dispatch-state/claims/model-training-604.json"
    other_required = required_paths - {f"{bundle['admission_artifact_root']}/ARTIFACTS.json"}
    if any(not (root / item).is_file() for item in other_required):
        raise ScBassetRectangleAdmissionError("non-gate queue dependency is absent")
    if not gate.exists() and claim.exists():
        raise ScBassetRectangleAdmissionError("queue item was claimed before admission")
    if gate.exists():
        admission = verify_frozen_tree(gate.parent)
        if (
            admission.get("metadata", {}).get("artifact_class")
            != "scbasset_three_seed_rectangle_admission"
            or admission.get("metadata", {}).get("contract_sha256")
            != _digest(contract_path)
            or admission.get("metadata", {}).get("status") != "passed"
        ):
            raise ScBassetRectangleAdmissionError("existing admission gate differs")


def audit_admission(root: Path, contract_path: Path) -> dict[str, Any]:
    """Verify the prospective rectangle without opening predictions or outcomes."""

    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    try:
        contract_path.relative_to(root)
    except ValueError as error:
        raise ScBassetRectangleAdmissionError("contract escapes benchmark root") from error
    contract = _load_json(contract_path)
    if contract.get("schema_version") != SCHEMA_VERSION or contract.get("status") != STATUS:
        raise ScBassetRectangleAdmissionError("rectangle schema or status differs")
    _audit_geometry(contract)
    _audit_readiness(root, contract)
    _audit_sources_and_inputs(root, contract)
    _audit_complete_artifacts(root, contract)
    _audit_runtime_and_bundle(root, contract)
    _audit_queue(root, contract)
    firewall = contract.get("evaluation_firewall", {})
    _assert_false(
        firewall,
        (
            "test_role_predictions_allowed",
            "prediction_values_may_be_read_during_admission",
            "metric_values_may_be_read_during_admission",
            "raw_outcomes_may_be_read_during_admission",
            "sealed_features_may_be_read",
            "sealed_labels_or_outcomes_may_be_read",
            "partial_results_may_be_ranked",
            "thresholds_may_be_changed",
            "development_evaluator_opens_before_full_rectangle",
            "champion_claim_allowed",
            "universal_claim_allowed",
        ),
    )
    return {
        "schema_version": "masld-bench-scbasset-three-seed-admission-receipt-v1",
        "status": "pass_exact_11_fit_rectangle_admitted_for_central_dispatch",
        "contract_sha256": _digest(contract_path),
        "expected_fits": 15,
        "complete_fits_at_freeze": 4,
        "missing_fits_admitted": 11,
        "fixed_seeds": list(SEEDS),
        "legacy_seeds_excluded": list(LEGACY_SEEDS),
        "central_bundle_id": "model-training-604",
        "central_bundle_count": 1,
        "logical_tasks": 22,
        "requested_gpu_allocations": 1,
        "requested_gpu_hours": 8,
        "central_queue_registration_present": True,
        "central_queue_unclaimed_during_validation": not (
            root
            / "executions/gpu-bundle-dispatch-state/claims/model-training-604.json"
        ).exists(),
        "admission_gate_present_during_validation": (
            root
            / f"{contract['production_bundle']['admission_artifact_root']}/ARTIFACTS.json"
        ).is_file(),
        "direct_gpu_submission_authorized": False,
        "manual_gpu_sbatch_allowed": False,
        "observed_multiome_eligibility_changed": False,
        "rna_conditioned_model_status_changed": False,
        "partial_ranking_authorized": False,
        "development_evaluator_open": False,
        "champion_claim_allowed": False,
        "prediction_values_read": False,
        "metric_values_read": False,
        "raw_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_labels_or_outcomes_read": False,
        "gpu_job_submitted_by_admission": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    receipt = audit_admission(arguments.project_root, arguments.contract)
    rendered = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(rendered, end="")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
