#!/usr/bin/env python3
"""Audit the outcome-blind scBasset five-seed rectangle extension."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any, Mapping

from masld_bench.artifacts import ArtifactError, verify_frozen_tree
from scripts.audit_scbasset_three_seed_rectangle_admission import (
    audit_admission as audit_three_seed_admission,
)
from scripts.validate_scbasset_broad_rg_prediction_line_exposure import (
    validate_incident,
)


SCHEMA_VERSION = "masld-bench-scbasset-five-seed-rectangle-v1"
STATUS = "PROSPECTIVE_FAMILY_NATIVE_FIVE_SEED_RECTANGLE_WAITING_FOR_CENTRAL_DISPATCH"
FOLDS = tuple(range(5))
SEEDS = (20260824, 20260825, 20260826, 20260827, 20260828)
LEGACY_SEEDS = (11, 29, 47, 71, 101)
COMPLETE_CONTEXTS = {
    (1, 20260824),
    (2, 20260824),
    (3, 20260824),
    (4, 20260824),
}
EXPECTED_CONTEXTS = {(fold, seed) for fold in FOLDS for seed in SEEDS}
MISSING_CONTEXTS = EXPECTED_CONTEXTS - COMPLETE_CONTEXTS
EXTENSION_CONTEXTS = {(fold, seed) for fold in FOLDS for seed in SEEDS[-2:]}


class ScBassetFiveSeedAdmissionError(ValueError):
    """Raised when a five-seed authority or fail-closed check differs."""


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
        raise ScBassetFiveSeedAdmissionError(f"JSON authority differs: {path}") from error
    if not isinstance(value, dict):
        raise ScBassetFiveSeedAdmissionError(f"JSON authority is not an object: {path}")
    return value


def _resolve_file(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ScBassetFiveSeedAdmissionError("authority path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ScBassetFiveSeedAdmissionError(f"authority is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise ScBassetFiveSeedAdmissionError(f"authority escapes or is absent: {relative}") from error
    if not resolved.is_file():
        raise ScBassetFiveSeedAdmissionError(f"authority is not a file: {relative}")
    return resolved


def _resolve_directory(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ScBassetFiveSeedAdmissionError("artifact path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ScBassetFiveSeedAdmissionError(f"artifact is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise ScBassetFiveSeedAdmissionError(f"artifact escapes or is absent: {relative}") from error
    if not resolved.is_dir():
        raise ScBassetFiveSeedAdmissionError(f"artifact is not a directory: {relative}")
    return resolved


def _require_file(root: Path, relative: str, expected_sha256: str) -> Path:
    path = _resolve_file(root, relative)
    if _digest(path) != expected_sha256:
        raise ScBassetFiveSeedAdmissionError(f"file authority differs: {relative}")
    return path


def _require_frozen(
    root: Path,
    relative: str,
    expected_sha256: str,
    artifact_class: str,
) -> dict[str, Any]:
    artifact_root = _resolve_directory(root, relative)
    try:
        manifest = verify_frozen_tree(artifact_root)
    except ArtifactError as error:
        raise ScBassetFiveSeedAdmissionError(f"frozen artifact differs: {relative}") from error
    if (
        _digest(artifact_root / "ARTIFACTS.json") != expected_sha256
        or manifest.get("metadata", {}).get("artifact_class") != artifact_class
        or manifest.get("metadata", {}).get("status") != "passed"
    ):
        raise ScBassetFiveSeedAdmissionError(f"frozen identity differs: {relative}")
    return dict(manifest)


def _assert_false(mapping: Mapping[str, Any], fields: tuple[str, ...]) -> None:
    for field in fields:
        if mapping.get(field) is not False:
            raise ScBassetFiveSeedAdmissionError(f"fail-closed field opened: {field}")


def _audit_geometry_and_views(contract: Mapping[str, Any]) -> None:
    complete = {tuple(row) for row in contract.get("complete_contexts", ())}
    missing = {tuple(row) for row in contract.get("missing_contexts", ())}
    old = {tuple(row) for row in contract.get("frozen_three_seed_authority", {}).get("covered_missing_contexts", ())}
    extension = {tuple(row) for row in contract.get("five_seed_extension", {}).get("covered_missing_contexts", ())}
    if (
        contract.get("model_id") != "scbasset"
        or contract.get("dataset_id") != "gse296875"
        or contract.get("task_role") != "outer_training_atac_supervised_sequence_profile_baseline"
        or tuple(contract.get("outer_folds", ())) != FOLDS
        or tuple(contract.get("fixed_seeds", ())) != SEEDS
        or tuple(contract.get("legacy_seeds_excluded", ())) != LEGACY_SEEDS
        or contract.get("expected_fits") != 25
        or contract.get("complete_compatible_fits_at_freeze") != 4
        or contract.get("missing_fits_at_freeze") != 21
        or complete != COMPLETE_CONTEXTS
        or missing != MISSING_CONTEXTS
        or complete & missing
        or complete | missing != EXPECTED_CONTEXTS
        or len(old) != 11
        or len(extension) != 10
        or extension != EXTENSION_CONTEXTS
        or old & extension
        or old | extension != MISSING_CONTEXTS
        or any(seed in LEGACY_SEEDS for _, seed in EXPECTED_CONTEXTS)
    ):
        raise ScBassetFiveSeedAdmissionError("five-seed rectangle geometry differs")

    view = contract.get("prediction_view_contract", {})
    if (
        view.get("view_id") != "valid_training_cell_state_mean_sequence_profile"
        or view.get("views_per_fit") != 1
        or view.get("expected_valid_prediction_views") != 25
        or view.get("complete_valid_prediction_views_at_freeze") != 4
        or view.get("missing_valid_prediction_views_at_freeze") != 21
        or view.get("donor_context") != "none_sequence_only"
        or view.get("held_donor_profiles_identical_within_lineage") is not True
        or view.get("observed_query_atac_consumed") is not False
        or view.get("rna_context_consumed") is not False
        or view.get("held_cell_embedding_available") is not False
        or view.get("native_training_cell_embedding_is_transductive_only") is not True
        or view.get("observed_multiome_inductive_claim_allowed") is not False
        or view.get("rna_conditioned_claim_allowed") is not False
        or view.get("standalone_signed_variant_effect_claim_allowed") is not False
    ):
        raise ScBassetFiveSeedAdmissionError("scBasset prediction-view boundary differs")


def _audit_three_seed_authority(root: Path, contract: Mapping[str, Any]) -> None:
    authority = contract.get("frozen_three_seed_authority", {})
    old_contract = _require_file(root, str(authority.get("contract_path", "")), str(authority.get("contract_sha256", "")))
    _require_file(root, str(authority.get("audit_path", "")), str(authority.get("audit_sha256", "")))
    _require_file(root, str(authority.get("queue_path", "")), str(authority.get("queue_sha256", "")))
    admission = _require_frozen(
        root,
        str(authority.get("admission_artifact_root", "")),
        str(authority.get("admission_artifacts_sha256", "")),
        "scbasset_three_seed_rectangle_admission",
    )
    old_receipt = audit_three_seed_admission(root, old_contract)
    if (
        old_receipt.get("status") != "pass_exact_11_fit_rectangle_admitted_for_central_dispatch"
        or old_receipt.get("complete_fits_at_freeze") != 4
        or old_receipt.get("missing_fits_admitted") != 11
        or old_receipt.get("central_bundle_id") != "model-training-604"
        or old_receipt.get("central_queue_unclaimed_during_validation") is not True
        or old_receipt.get("prediction_values_read") is not False
        or old_receipt.get("metric_values_read") is not False
        or admission.get("metadata", {}).get("logical_tasks") != 22
    ):
        raise ScBassetFiveSeedAdmissionError("three-seed coverage authority differs")


def _audit_reusable_fits(root: Path) -> None:
    model_pattern = re.compile(r"^scbasset-donor([0-4])_genomic\1-seed(2026082[4-8])-([0-9]+)$")
    prediction_pattern = re.compile(
        r"^scbasset-donor([0-4])_genomic\1-seed(2026082[4-8])-valid-([0-9]+)$"
    )
    models: dict[tuple[int, int, str], dict[str, Any]] = {}
    predictions: dict[tuple[int, int, str], dict[str, Any]] = {}
    for path in (root / "executions").iterdir():
        match = model_pattern.fullmatch(path.name)
        if match is not None and (path / "ARTIFACTS.json").is_file():
            models[(int(match.group(1)), int(match.group(2)), match.group(3))] = _load_json(
                path / "ARTIFACTS.json"
            ).get("metadata", {})
        match = prediction_pattern.fullmatch(path.name)
        if match is not None and (path / "ARTIFACTS.json").is_file():
            predictions[(int(match.group(1)), int(match.group(2)), match.group(3))] = _load_json(
                path / "ARTIFACTS.json"
            ).get("metadata", {})
    complete = set()
    for key, model in models.items():
        prediction = predictions.get(key, {})
        if (
            model.get("artifact_class") == "scbasset_fold_model"
            and model.get("status") == "passed"
            and model.get("held_donor_atac_used") is False
            and model.get("genomic_test_atac_used") is False
            and prediction.get("artifact_class") == "scbasset_fold_predictions"
            and prediction.get("status") == "passed"
            and prediction.get("evaluation_role") == "valid"
            and prediction.get("held_donor_atac_used") is False
            and prediction.get("held_cell_embedding_available") is False
        ):
            complete.add(key[:2])
    if complete != COMPLETE_CONTEXTS:
        raise ScBassetFiveSeedAdmissionError("reusable scBasset fit census changed")


def _audit_derivation_and_bundle(root: Path, contract: Mapping[str, Any]) -> None:
    derivation = contract.get("source_wrapper_derivation", {})
    for role, path_key, hash_key in (
        ("training", "training_wrapper_path", "training_wrapper_sha256"),
        ("prediction", "prediction_wrapper_path", "prediction_wrapper_sha256"),
    ):
        path = _require_file(root, str(derivation.get(path_key, "")), str(derivation.get(hash_key, "")))
        text = path.read_text(encoding="utf-8")
        old_gate = str(derivation.get("old_seed_gate", ""))
        new_gate = str(derivation.get("new_seed_gate", ""))
        if (
            text.count(old_gate) != derivation.get("expected_replacements_per_wrapper")
            or old_gate not in text
            or new_gate in text
            or "--qos=innovation" in text
        ):
            raise ScBassetFiveSeedAdmissionError(f"{role} wrapper derivation differs")
        derived = text.replace(old_gate, new_gate)
        if derived.count(new_gate) != derivation.get("expected_replacements_per_wrapper"):
            raise ScBassetFiveSeedAdmissionError(f"{role} derived wrapper differs")
    if (
        derivation.get("other_source_changes_allowed") is not False
        or derivation.get("derived_wrappers_frozen_with_bundle") is not True
    ):
        raise ScBassetFiveSeedAdmissionError("derived-wrapper provenance differs")

    bundle = contract.get("five_seed_extension", {})
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
        "#SBATCH --time=12:00:00",
        "SCBASSET_FIVE_SEED_AUTHORITY",
        "SCBASSET_PREFLIGHT_ONLY",
        "continue_after_isolated_fit_failure",
        "derived_wrappers",
    )
    if (
        bundle.get("bundle_id") != "model-training-605"
        or bundle.get("queue_priority") != 64
        or bundle.get("fits") != 10
        or bundle.get("logical_tasks") != 20
        or bundle.get("array") is not False
        or bundle.get("submission_authority") != "central_dispatcher_only"
        or bundle.get("direct_gpu_submission_authorized") is not False
        or bundle.get("manual_gpu_sbatch_allowed") is not False
        or bundle.get("gpu_job_ceiling") != 5
        or bundle.get("maximum_running_gpu_jobs") != 4
        or bundle.get("maximum_pending_gpu_jobs") != 1
        or any(token not in text for token in required)
        or "#SBATCH --array" in text
        or "--qos=innovation" in text
        or re.search(r"^\s*sbatch\s", text, flags=re.MULTILINE)
    ):
        raise ScBassetFiveSeedAdmissionError("five-seed bundle resource or authority differs")


def _audit_queue(root: Path, contract_path: Path, contract: Mapping[str, Any]) -> None:
    bundle = contract["five_seed_extension"]
    queue_path = _resolve_file(root, "config/campaigns/gpu_bundle_queue/064-model-training-605.json")
    queue = _load_json(queue_path)
    wrapper = _resolve_file(root, str(bundle["wrapper_path"]))
    old_contract = _load_json(
        _resolve_file(root, str(contract["frozen_three_seed_authority"]["contract_path"]))
    )
    required_paths = {
        "config/campaigns/scbasset_five_seed_rectangle_20260825.json",
        "config/artifacts/incidents/scbasset_broad_rg_prediction_line_exposure_20260825.json",
        f"{bundle['admission_artifact_root']}/ARTIFACTS.json",
        f"{contract['frozen_three_seed_authority']['admission_artifact_root']}/ARTIFACTS.json",
        f"{old_contract['readiness_authority']['artifact_root']}/ARTIFACTS.json",
        *(f"{record['artifact_root']}/ARTIFACTS.json" for record in old_contract["input_artifacts"]),
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
        or queue.get("bundle_id") != "model-training-605"
        or queue.get("priority") != 64
        or queue.get("enabled") is not True
        or queue.get("wrapper_path") != bundle.get("wrapper_path")
        or queue.get("wrapper_sha256") != _digest(wrapper)
        or queue.get("logical_tasks") != 20
        or queue.get("family") != "single_cell_chromatin"
        or set(queue.get("required_paths", ())) != required_paths
        or queue.get("exports", {}).get("SCBASSET_FIVE_SEED_AUTHORITY")
        != "scbasset_five_seed_rectangle_v1"
        or queue.get("exports", {}).get("SCBASSET_FIVE_SEED_CONTRACT") != str(contract_path)
        or queue.get("exports", {}).get("SCBASSET_FIVE_SEED_CONTRACT_SHA256")
        != _digest(contract_path)
    ):
        raise ScBassetFiveSeedAdmissionError("five-seed central queue item differs")
    gate = root / f"{bundle['admission_artifact_root']}/ARTIFACTS.json"
    claim = root / "executions/gpu-bundle-dispatch-state/claims/model-training-605.json"
    other_required = required_paths - {f"{bundle['admission_artifact_root']}/ARTIFACTS.json"}
    if any(not (root / relative).is_file() for relative in other_required):
        raise ScBassetFiveSeedAdmissionError("non-gate queue dependency is absent")
    if claim.exists():
        raise ScBassetFiveSeedAdmissionError("five-seed queue item was claimed before admission")
    if gate.exists():
        admission = verify_frozen_tree(gate.parent)
        if (
            admission.get("metadata", {}).get("artifact_class")
            != "scbasset_five_seed_rectangle_admission"
            or admission.get("metadata", {}).get("contract_sha256") != _digest(contract_path)
            or admission.get("metadata", {}).get("status") != "passed"
        ):
            raise ScBassetFiveSeedAdmissionError("existing five-seed gate differs")


def audit_admission(root: Path, contract_path: Path) -> dict[str, Any]:
    """Verify the five-seed campaign using manifests and structural receipts only."""

    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    try:
        contract_path.relative_to(root)
    except ValueError as error:
        raise ScBassetFiveSeedAdmissionError("contract escapes benchmark root") from error
    contract = _load_json(contract_path)
    if contract.get("schema_version") != SCHEMA_VERSION or contract.get("status") != STATUS:
        raise ScBassetFiveSeedAdmissionError("five-seed schema or status differs")
    _audit_geometry_and_views(contract)
    _audit_three_seed_authority(root, contract)
    _audit_reusable_fits(root)
    _audit_derivation_and_bundle(root, contract)
    _audit_queue(root, contract_path, contract)

    firewall = contract.get("evaluation_firewall", {})
    _assert_false(
        firewall,
        (
            "test_role_predictions_allowed",
            "development_evaluator_opens_before_full_rectangle",
            "partial_results_may_be_ranked",
            "scbasset_prediction_values_may_be_read_during_admission",
            "metric_values_may_be_read_during_admission",
            "raw_outcomes_may_be_read_during_admission",
            "sealed_features_may_be_read",
            "sealed_labels_or_outcomes_may_be_read",
            "thresholds_may_be_changed",
            "champion_claim_allowed",
            "universal_claim_allowed",
        ),
    )
    incident = contract.get("development_search_incident", {})
    incident_path = _require_file(
        root,
        str(incident.get("incident_path", "")),
        str(incident.get("incident_sha256", "")),
    )
    incident_receipt = validate_incident(incident_path)
    if (
        incident.get("incident_id")
        != "scbasset_broad_rg_prediction_line_exposure_20260825"
        or
        incident.get("unrelated_prediction_lines_accidentally_emitted") is not True
        or incident.get("scbasset_prediction_values_emitted") is not False
        or incident.get("metric_values_emitted") is not False
        or incident.get("raw_outcomes_emitted") is not False
        or incident.get("sealed_data_emitted") is not False
        or incident.get("values_used_for_rectangle_or_execution_design") is not False
        or incident.get("incident_disclosed_to_parent_agent") is not True
        or incident_receipt.get("status")
        != "pass_structure_only_incident_validation"
        or incident_receipt.get("incident_values_used") is not False
        or incident_receipt.get("referenced_prediction_files_reopened_by_validator")
        is not False
    ):
        raise ScBassetFiveSeedAdmissionError("development-search incident record differs")
    return {
        "schema_version": "masld-bench-scbasset-five-seed-admission-receipt-v1",
        "status": "pass_exact_21_missing_fit_rectangle_covered_by_two_central_bundles",
        "contract_sha256": _digest(contract_path),
        "expected_fits": 25,
        "complete_fits_at_freeze": 4,
        "missing_fits_admitted": 21,
        "expected_valid_prediction_views": 25,
        "complete_valid_prediction_views_at_freeze": 4,
        "missing_valid_prediction_views_admitted": 21,
        "prediction_view": "valid_training_cell_state_mean_sequence_profile",
        "donor_context": "none_sequence_only",
        "observed_query_atac_consumed": False,
        "rna_context_consumed": False,
        "held_cell_embedding_available": False,
        "observed_multiome_inductive_claim_allowed": False,
        "rna_conditioned_claim_allowed": False,
        "existing_bundle_id": "model-training-604",
        "existing_bundle_missing_fits": 11,
        "extension_bundle_id": "model-training-605",
        "extension_bundle_missing_fits": 10,
        "central_bundle_count": 2,
        "logical_tasks": 42,
        "requested_new_gpu_allocations": 1,
        "requested_new_gpu_hours": 12,
        "central_queue_registration_present": True,
        "central_queue_unclaimed_during_validation": True,
        "direct_gpu_submission_authorized": False,
        "manual_gpu_sbatch_allowed": False,
        "gpu_squeue_ceiling": 5,
        "maximum_running_gpu_jobs": 4,
        "maximum_pending_gpu_jobs": 1,
        "partial_ranking_authorized": False,
        "development_evaluator_open": False,
        "champion_claim_allowed": False,
        "admission_audit_scbasset_prediction_values_read": False,
        "metric_values_read": False,
        "raw_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_labels_or_outcomes_read": False,
        "development_search_incident_recorded": True,
        "development_search_incident_id": incident_receipt["incident_id"],
        "development_search_incident_sha256": incident_receipt["incident_sha256"],
        "development_search_incident_structure_validated": True,
        "unrelated_prediction_lines_accidentally_emitted": True,
        "incident_values_used": False,
        "gpu_job_submitted_by_admission": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = audit_admission(args.project_root, args.contract)
    rendered = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is not None:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
