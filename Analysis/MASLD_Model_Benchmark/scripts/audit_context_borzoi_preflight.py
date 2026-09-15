#!/usr/bin/env python3
"""Audit the outcome-blind preflight for the generic conditional context slot."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any


SCHEMA_VERSION = "masld-bench-context-borzoi-outcome-blind-preflight-v1"
BLOCKED_STATUS = "PREFLIGHT_ONLY_BLOCKED_TRIGGER_AND_COMPONENT_ADMISSION"
ABLATIONS = [
    "sequence_only",
    "trans_only",
    "permuted_context",
    "reverse_complement",
    "matched_shuffle",
]
MISSING_STATES = [
    "observed",
    "structurally_missing",
    "not_applicable",
    "below_qc",
    "unavailable_permission",
    "join_unresolved",
    "withheld_sealed",
    "derivable_not_processed",
]


class ContextBorzoiPreflightError(ValueError):
    """Raised when the conditional architecture preflight drifts or opens a check."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _canonical_sha256(value: Any) -> str:
    payload = json.dumps(
        value,
        sort_keys=True,
        separators=(",", ":"),
        ensure_ascii=False,
        allow_nan=False,
    ).encode("utf-8")
    return sha256(payload).hexdigest()


def _resolve_member(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ContextBorzoiPreflightError("preflight member path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ContextBorzoiPreflightError(f"preflight member is a symlink: {relative}")
    member = candidate.resolve(strict=True)
    try:
        member.relative_to(root)
    except ValueError as error:
        raise ContextBorzoiPreflightError("preflight member escapes benchmark root") from error
    if not member.is_file():
        raise ContextBorzoiPreflightError(f"preflight member is not a file: {relative}")
    return member


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ContextBorzoiPreflightError(f"JSON authority is not an object: {path}")
    return value


def _load_toml(path: Path) -> dict[str, Any]:
    with path.open("rb") as handle:
        value = tomllib.load(handle)
    if not isinstance(value, dict):
        raise ContextBorzoiPreflightError(f"TOML authority is not an object: {path}")
    return value


def _artifact_index(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise ContextBorzoiPreflightError("artifact inventory differs")
    result: dict[str, dict[str, Any]] = {}
    for record in artifacts:
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            raise ContextBorzoiPreflightError("artifact record differs")
        if record["path"] in result:
            raise ContextBorzoiPreflightError("artifact inventory contains duplicates")
        result[record["path"]] = record
    return result


def audit_preflight(root: Path, authority_path: Path) -> dict[str, Any]:
    """Verify only code/config authorities and outcome-blind fixture metadata."""

    root = root.resolve(strict=True)
    authority_path = authority_path.resolve(strict=True)
    try:
        authority_path.relative_to(root)
    except ValueError as error:
        raise ContextBorzoiPreflightError("preflight authority escapes benchmark root") from error
    authority = _load_json(authority_path)
    if authority.get("schema_version") != SCHEMA_VERSION:
        raise ContextBorzoiPreflightError("conditional preflight schema differs")
    if authority.get("status") != BLOCKED_STATUS:
        raise ContextBorzoiPreflightError("conditional preflight status is not blocked")

    authority_hashes = authority.get("authority_hashes")
    implementation_hashes = authority.get("implementation_hashes")
    manifest_hashes = authority.get("artifact_manifest_hashes")
    if not isinstance(authority_hashes, dict) or len(authority_hashes) != 11:
        raise ContextBorzoiPreflightError("conditional authority hash family differs")
    if not isinstance(implementation_hashes, dict) or len(implementation_hashes) != 9:
        raise ContextBorzoiPreflightError("conditional implementation hash family differs")
    if not isinstance(manifest_hashes, dict) or len(manifest_hashes) != 2:
        raise ContextBorzoiPreflightError("conditional fixture-manifest family differs")

    verified: dict[str, str] = {}
    for family in (authority_hashes, implementation_hashes):
        for relative, expected in sorted(family.items()):
            if not isinstance(relative, str) or not isinstance(expected, str):
                raise ContextBorzoiPreflightError("conditional hash record differs")
            member = _resolve_member(root, relative)
            observed = _digest(member)
            if observed != expected:
                raise ContextBorzoiPreflightError(f"conditional authority drifted: {relative}")
            verified[relative] = observed

    manifests: dict[str, dict[str, Any]] = {}
    for relative, expected in sorted(manifest_hashes.items()):
        if (
            not isinstance(relative, str)
            or not relative.startswith("executions/")
            or not relative.endswith("/ARTIFACTS.json")
            or not isinstance(expected, str)
        ):
            raise ContextBorzoiPreflightError("fixture manifest allowlist differs")
        member = _resolve_member(root, relative)
        observed = _digest(member)
        if observed != expected:
            raise ContextBorzoiPreflightError(f"conditional fixture drifted: {relative}")
        manifest = _load_json(member)
        if manifest.get("schema_version") != "masld-bench-artifacts-v1":
            raise ContextBorzoiPreflightError("conditional fixture schema differs")
        manifests[relative] = manifest
        verified[relative] = observed

    campaign = _load_toml(
        _resolve_member(root, "config/campaigns/v1_conditional_model.toml")
    )
    if (
        campaign.get("campaign_id") != "v1-conditional-context-model"
        or campaign.get("wave") != "conditional_model"
        or campaign.get("status") != "blocked_trigger"
        or campaign.get("submit_enabled") is not False
        or campaign.get("conditional_model_spec_path") != "UNRESOLVED"
        or campaign.get("conditional_model_spec_sha256") != "UNRESOLVED"
        or campaign.get("task_ids")
        != ["variant_to_regulation", "rna_conditioned_atac"]
        or campaign.get("execution", {}).get("architecture_count_cap") != 1
    ):
        raise ContextBorzoiPreflightError("conditional campaign is not fail-closed")

    context_checkpoint = _load_json(
        _resolve_member(root, "config/artifacts/models/context_borzoi/checkpoints.json")
    )
    context_crosswalk = _load_json(
        _resolve_member(
            root, "config/artifacts/models/context_borzoi/development_crosswalk.json"
        )
    )
    context_exposure = _load_json(
        _resolve_member(root, "config/artifacts/models/context_borzoi/exposure_audit.json")
    )
    borzoi_checkpoint = _load_json(
        _resolve_member(root, "config/artifacts/models/borzoi/checkpoints.json")
    )
    borzoi_exposure = _load_json(
        _resolve_member(root, "config/artifacts/models/borzoi/exposure_audit.json")
    )
    corgi_checkpoint = _load_json(
        _resolve_member(root, "config/artifacts/models/corgi/checkpoints.json")
    )
    corgi_exposure = _load_json(
        _resolve_member(root, "config/artifacts/models/corgi/exposure_audit.json")
    )

    identity = authority.get("identity", {})
    authorization = context_checkpoint.get("authorization_contract", {})
    checkpoint = context_checkpoint.get("checkpoint", {})
    components = context_checkpoint.get("component_selection", {})
    architecture = context_checkpoint.get("architecture_authority", {})
    noninheritance = context_checkpoint.get("noninheritance_contract", {})
    if (
        identity.get("model_id") != "context_borzoi"
        or identity.get("generic_conditional_slot") is not True
        or identity.get("borzoi_component_selected") is not False
        or identity.get("corgi_component_selected") is not False
        or identity.get("selected_sequence_backbone") is not None
        or identity.get("selected_cell_encoder") is not None
        or identity.get("trained_checkpoint") is not None
        or identity.get("architecture_count_cap") != 1
        or identity.get("zero_shot_lane_eligible") is not False
        or checkpoint.get("artifact_exists") is not False
        or checkpoint.get("checkpoint_sha256") != "NOT_APPLICABLE"
        or checkpoint.get("weight_downloaded") is not False
        or components.get("selection_lock") != "ABSENT"
        or not str(components.get("base_sequence_backbone", "")).startswith(
            "UNSELECTED"
        )
        or not str(components.get("cell_state_encoder", "")).startswith("UNSELECTED")
        or noninheritance.get("borzoi_checkpoint_bytes_inherited") is not False
    ):
        raise ContextBorzoiPreflightError("conditional model identity is no longer empty")

    gate = authority.get("complementarity_gate", {})
    required_evidence = authorization.get("required_evidence", {})
    if (
        gate.get("state") != "NOT_EVALUATED"
        or gate.get("build_authorized") is not False
        or authorization.get("trigger_result") != "NOT_EVALUATED"
        or authorization.get("build_authorized") is not False
        or gate.get("residual_correlation_strictly_below")
        != required_evidence.get("residual_correlation_strictly_below")
        or gate.get("qualifying_relative_deviance_gain_at_least")
        != required_evidence.get("qualifying_relative_deviance_gain")
        or gate.get("qualifying_absolute_correlation_or_f1_gain_at_least")
        != required_evidence.get("qualifying_absolute_correlation_or_f1_gain")
        or gate.get("minimum_qualifying_mechanistic_endpoints")
        != required_evidence.get("minimum_qualifying_mechanistic_endpoints")
        or gate.get("evaluated_seed_count")
        != required_evidence.get("evaluated_seed_count")
        or gate.get("minimum_qualifying_seed_count")
        != required_evidence.get("gain_persists_in_at_least_seeds")
        or gate.get("minimum_study_count") != required_evidence.get("minimum_studies")
        or gate.get("cross_fitted") is not True
        or gate.get("nonnegative_stack") is not True
        or gate.get("best_open_models_compared") is not True
        or gate.get("qualifying_endpoint_keys_must_include")
        != ["rna_conditioned_atac", "variant_to_regulation"]
    ):
        raise ContextBorzoiPreflightError("complementarity trigger contract differs")

    if (
        architecture.get("architecture_count_cap") != 1
        or architecture.get("conditioning") != ["gated_film", "lora_adapters"]
        or architecture.get("heads")
        != ["atac_profile", "atac_count", "rna_coverage"]
        or architecture.get("forbidden_heads") != ["masld_histone"]
        or architecture.get("required_ablations") != ABLATIONS
        or architecture.get("missing_modality_masks") is not True
        or architecture.get("status") != "specified_but_not_implemented"
    ):
        raise ContextBorzoiPreflightError("conditional architecture roster differs")

    component_gate = authority.get("component_admission_gate", {})
    if (
        component_gate.get("state") != "BLOCKED_NO_SELECTED_COMPONENTS"
        or component_gate.get("component_findings_inherited_by_composite") is not False
        or component_gate.get("component_weights_may_be_opened_now") is not False
        or component_gate.get("local_derivative_may_be_built_now") is not False
        or borzoi_checkpoint.get("safe_loading_contract", {}).get(
            "checkpoint_content_downloaded"
        )
        is not False
        or borzoi_checkpoint.get("license_audit", {}).get("weight_use_status")
        != "BLOCKED_PENDING_EXPLICIT_WEIGHT_TERMS_OR_INSTITUTIONAL_LEGAL_DETERMINATION"
        or borzoi_exposure.get("sealed_champion_eligibility")
        != "exposure_eligible_but_license_and_runtime_blocked"
        or corgi_checkpoint.get("verified_artifacts", {}).get(
            "native_numeric_parity_performed"
        )
        is not False
        or not str(corgi_checkpoint.get("admission_status", "")).startswith("BLOCKED")
        or corgi_exposure.get("license_disposition", {}).get("open_champion_eligibility")
        != "ineligible_until_code_and_derivative_weight_terms_are_cleared"
    ):
        raise ContextBorzoiPreflightError("component admission gate differs")

    expected_missing = set(MISSING_STATES)
    missing_contract = authority.get("input_and_missingness_contract", {})
    crosswalk_missing = context_crosswalk.get("missingness_and_pairing_contract", {})
    if (
        set(missing_contract.get("missingness_states", [])) != expected_missing
        or set(crosswalk_missing.get("allowed_missingness_states", []))
        != expected_missing
        or missing_contract.get("missing_feature")
        != {"mask": False, "value": None, "numeric_zero_allowed": False}
        or crosswalk_missing.get("missing_as_zero") != "forbidden"
        or context_checkpoint.get("input_contract", {}).get("atac_inference_input")
        != "forbidden"
        or context_checkpoint.get("input_contract", {}).get("sealed_outcome_inputs")
        != "forbidden"
        or missing_contract.get("missing_context_film_lora_behavior")
        != "context gate exactly zero and neutral residual identity"
    ):
        raise ContextBorzoiPreflightError("conditional missingness contract differs")

    fixture_manifest_path = (
        "executions/context-conditioning-precondition-21064508/ARTIFACTS.json"
    )
    fixture_manifest = manifests[fixture_manifest_path]
    fixture_metadata = fixture_manifest.get("metadata")
    if fixture_metadata != {
        "artifact_class": "context_conditioning_precondition_fixture",
        "status": "pass",
        "conditional_model_authorized": False,
        "architecture_built": False,
        "model_training_started": False,
        "sealed_outcomes_exposed": False,
    }:
        raise ContextBorzoiPreflightError("context precondition metadata differs")
    fixture_index = _artifact_index(fixture_manifest)
    fixture_path = str(Path(fixture_manifest_path).parent / "fixture.json")
    fixture_member = _resolve_member(root, fixture_path)
    if _digest(fixture_member) != fixture_index.get("fixture.json", {}).get("sha256"):
        raise ContextBorzoiPreflightError("context precondition fixture drifted")
    fixture = _load_json(fixture_member)
    if (
        fixture.get("status") != "pass"
        or fixture.get("donor_overlap") is not False
        or fixture.get("genomic_block_overlap") is not False
        or fixture.get("observed_atac_input_exposed") is not False
        or fixture.get("sealed_outcomes_exposed") is not False
        or fixture.get("missing_context_numeric_zero_allowed") is not False
        or fixture.get("neutral_film_identity") is not True
        or fixture.get("zero_initialized_lora_identity") is not True
        or fixture.get("architecture_built") is not False
        or fixture.get("model_training_started") is not False
        or fixture.get("conditional_model_authorized") is not False
    ):
        raise ContextBorzoiPreflightError("context precondition safety fixture differs")

    ablation_manifest_path = (
        "executions/conditional-ablation-preconditions-21064541/ARTIFACTS.json"
    )
    ablation_manifest = manifests[ablation_manifest_path]
    ablation_metadata = ablation_manifest.get("metadata")
    if (
        not isinstance(ablation_metadata, dict)
        or ablation_metadata.get("artifact_class")
        != "conditional_ablation_preconditions"
        or ablation_metadata.get("status") != "precondition_only_parent_unresolved"
        or ablation_metadata.get("ablation_count") != 5
        or ablation_metadata.get("architecture_built") is not False
        or ablation_metadata.get("model_training_started") is not False
        or ablation_metadata.get("conditional_model_authorized") is not False
        or ablation_metadata.get("sealed_outcomes_exposed") is not False
    ):
        raise ContextBorzoiPreflightError("conditional ablation metadata differs")
    ablation_index = _artifact_index(ablation_manifest)
    ablation_path = str(
        Path(ablation_manifest_path).parent / "conditional_ablation_preconditions.json"
    )
    ablation_member = _resolve_member(root, ablation_path)
    if _digest(ablation_member) != ablation_index.get(
        "conditional_ablation_preconditions.json", {}
    ).get("sha256"):
        raise ContextBorzoiPreflightError("conditional ablation document drifted")
    ablation = _load_json(ablation_member)
    identity_payload = dict(ablation)
    claimed_identity = identity_payload.pop("precondition_spec_id", None)
    if (
        claimed_identity != ablation_metadata.get("precondition_spec_id")
        or _canonical_sha256(identity_payload) != claimed_identity
        or ablation.get("ablation_names") != ABLATIONS
        or set(ablation.get("ablations", {})) != set(ABLATIONS)
        or authority.get("ablation_roster", {}).get("names") != ABLATIONS
    ):
        raise ContextBorzoiPreflightError("conditional ablation identity differs")
    for name in ABLATIONS:
        spec = ablation["ablations"][name]
        if (
            spec.get("status") != "precondition_only_parent_unresolved"
            or spec.get("parent_binding_required") is not True
            or spec.get("same_parent_rows_output_axis_splits_seeds_and_evaluator_required")
            is not True
            or spec.get("outcome_blind") is not True
            or spec.get("sealed_outcomes_allowed") is not False
            or spec.get("observed_query_atac_allowed") is not False
            or spec.get("independent_champion_identity") is not False
            or spec.get("executable") is not False
            or spec.get("model_training_allowed") is not False
            or spec.get("future_realization_requires_new_frozen_parent_bound_spec")
            is not True
        ):
            raise ContextBorzoiPreflightError(f"conditional ablation {name} differs")

    variant_capability = _load_toml(
        _resolve_member(
            root, "config/evaluation/variant_to_regulation_capabilities.toml"
        )
    )
    records = [
        record
        for record in variant_capability.get("models", [])
        if record.get("model_id") == "context_borzoi"
    ]
    if len(records) != 1:
        raise ContextBorzoiPreflightError("context variant capability record differs")
    variant_record = records[0]
    if (
        variant_record.get("role") != "conditional_only"
        or variant_record.get("native_outputs")
        != ["accessibility_delta", "gene_expression_delta"]
        or variant_record.get("requires_fitted_head") is not True
        or variant_record.get("requires_observed_target_context") is not False
    ):
        raise ContextBorzoiPreflightError("context variant native boundary differs")
    task_boundaries = authority.get("family_native_task_boundaries", {})
    if (
        task_boundaries.get("rna_conditioned_atac", {}).get("profile_fixture_passed")
        is not False
        or task_boundaries.get("rna_conditioned_atac", {}).get("scoring_authorized")
        is not False
        or task_boundaries.get("variant_to_regulation", {}).get("scoring_authorized")
        is not False
        or not str(
            task_boundaries.get("variant_to_regulation", {}).get(
                "signed_cell_type_eqtl_effect", ""
            )
        ).startswith("BLOCKED")
    ):
        raise ContextBorzoiPreflightError("conditional task boundary opened")

    drift = authority.get("known_authority_drift_blockers")
    if not isinstance(drift, list) or len(drift) != 4:
        raise ContextBorzoiPreflightError("known authority-drift blocker family differs")
    for finding in drift:
        if not isinstance(finding, dict):
            raise ContextBorzoiPreflightError("authority-drift finding differs")
        observed = finding.get("observed_sha256")
        path = finding.get("path")
        if not isinstance(path, str) or not isinstance(observed, str):
            raise ContextBorzoiPreflightError("authority-drift binding differs")
        if _digest(_resolve_member(root, path)) != observed:
            raise ContextBorzoiPreflightError(f"known authority drift changed: {path}")
        embedded = finding.get("embedded_sha256")
        if embedded is None:
            embedded = finding.get("embedded_in_exposure_audit_sha256")
        if embedded == observed:
            raise ContextBorzoiPreflightError(f"authority drift unexpectedly closed: {path}")

    if (
        context_exposure.get("conditional_status", {}).get("trigger_result")
        != "NOT_EVALUATED"
        or context_exposure.get("conditional_status", {}).get("build_authorized")
        is not False
        or context_exposure.get("checkpoint_findings", {})
        .get("context_borzoi", {})
        .get("exposure_state")
        != "unknown"
        or any(
            finding.get("exposure_state") != "unknown"
            for finding in context_crosswalk.get(
                "current_model_exposure_findings", {}
            ).values()
        )
    ):
        raise ContextBorzoiPreflightError("conditional exposure boundary differs")

    disposition = authority.get("execution_disposition", {})
    if disposition.get("cpu_metadata_preflight_authorized") is not True or any(
        disposition.get(field) is not False
        for field in (
            "architecture_import_authorized",
            "checkpoint_download_or_open_authorized",
            "mapper_fit_authorized",
            "model_build_authorized",
            "training_authorized",
            "prediction_authorized",
            "development_scoring_authorized",
            "sealed_scoring_authorized",
            "gpu_sbatch_authorized",
        )
    ):
        raise ContextBorzoiPreflightError("conditional execution disposition opened")

    return {
        "schema_version": "masld-bench-context-borzoi-preflight-audit-receipt-v1",
        "status": "pass_outcome_blind_preflight_blocked",
        "authority_sha256": _digest(authority_path),
        "verified_metadata_files": len(verified),
        "complementarity_gate_passed": False,
        "component_admission_passed": False,
        "architecture_executable": False,
        "training_authorized": False,
        "gpu_sbatch_authorized": False,
        "ablation_count": 5,
        "missing_as_zero_allowed": False,
        "checkpoint_opened": False,
        "development_outcomes_opened": False,
        "sealed_data_opened": False,
        "evaluator_outputs_opened": False,
        "borzoi_weight_admitted": False,
        "corgi_weight_admitted_as_component": False,
    }


def main() -> None:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument(
        "--authority",
        type=Path,
        default=default_root
        / "config/artifacts/models/context_borzoi/outcome_blind_architecture_preflight_20260824.json",
    )
    parser.add_argument(
        "--require-executable",
        action="store_true",
        help="Fail after audit because both activation gates are currently closed.",
    )
    arguments = parser.parse_args()
    receipt = audit_preflight(arguments.root, arguments.authority)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if arguments.require_executable:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
