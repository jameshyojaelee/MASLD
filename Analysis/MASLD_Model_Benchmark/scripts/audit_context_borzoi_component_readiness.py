#!/usr/bin/env python3
"""Audit the outcome-blind context-Borzoi component interface requirements."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


SCHEMA_VERSION = "masld-bench-context-borzoi-component-interface-readiness-v1"
STATUS = "INTERFACE_CONTRACT_READY_BUILD_BLOCKED"
MISSING_STATES = {
    "observed",
    "structurally_missing",
    "not_applicable",
    "below_qc",
    "unavailable_permission",
    "join_unresolved",
    "withheld_sealed",
    "derivable_not_processed",
}
EXECUTION_FLAGS = (
    "checkpoint_download_or_open_authorized",
    "component_selection_authorized",
    "architecture_build_authorized",
    "mapper_fit_authorized",
    "training_authorized",
    "prediction_authorized",
    "development_scoring_authorized",
    "sealed_scoring_authorized",
    "gpu_bundle_registration_authorized",
    "gpu_job_submission_authorized",
)


class ContextBorzoiReadinessError(ValueError):
    """Raised when a readiness authority drifts or opens a closed check."""


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
        raise ContextBorzoiReadinessError(f"JSON authority differs: {path}") from error
    if not isinstance(value, dict):
        raise ContextBorzoiReadinessError(f"JSON authority is not an object: {path}")
    return value


def _resolve_file(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ContextBorzoiReadinessError("authority path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ContextBorzoiReadinessError(f"authority is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise ContextBorzoiReadinessError(
            f"authority escapes or is absent: {relative}"
        ) from error
    if not resolved.is_file():
        raise ContextBorzoiReadinessError(f"authority is not a file: {relative}")
    return resolved


def _require_frozen_artifact(
    root: Path,
    relative: str,
    expected_sha256: str,
    *,
    artifact_class: str,
) -> dict[str, Any]:
    artifact_root = root / relative
    if artifact_root.is_symlink():
        raise ContextBorzoiReadinessError(f"artifact is a symlink: {relative}")
    try:
        resolved = artifact_root.resolve(strict=True)
        resolved.relative_to(root)
        manifest = verify_frozen_tree(resolved)
    except (OSError, ValueError, ArtifactError) as error:
        raise ContextBorzoiReadinessError(
            f"frozen artifact differs: {relative}"
        ) from error
    if (
        _digest(resolved / "ARTIFACTS.json") != expected_sha256
        or manifest.get("metadata", {}).get("artifact_class") != artifact_class
    ):
        raise ContextBorzoiReadinessError(
            f"frozen artifact identity differs: {relative}"
        )
    return dict(manifest)


def _assert_false(mapping: Mapping[str, Any], names: tuple[str, ...]) -> None:
    for name in names:
        if mapping.get(name) is not False:
            raise ContextBorzoiReadinessError(f"fail-closed gate opened: {name}")


def audit_readiness(root: Path, contract_path: Path) -> dict[str, Any]:
    """Verify only frozen metadata, interfaces, and closed execution checks."""

    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    try:
        contract_path.relative_to(root)
    except ValueError as error:
        raise ContextBorzoiReadinessError("contract escapes benchmark root") from error
    contract = _load_json(contract_path)
    if contract.get("schema_version") != SCHEMA_VERSION or contract.get("status") != STATUS:
        raise ContextBorzoiReadinessError("readiness schema or status differs")

    slot = contract.get("conditional_slot", {})
    preflight_path = _resolve_file(root, slot.get("existing_preflight_path", ""))
    preflight = _load_json(preflight_path)
    if (
        _digest(preflight_path) != slot.get("existing_preflight_sha256")
        or preflight.get("status")
        != "PREFLIGHT_ONLY_BLOCKED_TRIGGER_AND_COMPONENT_ADMISSION"
        or preflight.get("complementarity_gate", {}).get("state") != "NOT_EVALUATED"
        or preflight.get("complementarity_gate", {}).get("build_authorized") is not False
        or slot.get("model_id") != "context_borzoi"
        or slot.get("architecture_count_cap") != 1
        or slot.get("component_selection_lock") != "ABSENT"
        or slot.get("selected_sequence_backbone") is not None
        or slot.get("selected_cell_encoder") is not None
        or slot.get("architecture_built") is not False
        or slot.get("training_started") is not False
    ):
        raise ContextBorzoiReadinessError("conditional slot is no longer empty and blocked")

    sequence = contract.get("provisional_sequence_latent_interface", {})
    if (
        sequence.get("status") != "SCHEMA_READY_COMPONENT_UNSELECTED"
        or sequence.get("input_axis_order") != ["batch", "nucleotide", "base"]
        or sequence.get("input_shape") != [None, 4, 524288]
        or sequence.get("latent_axis_order")
        != ["batch", "channel", "central_bin"]
        or sequence.get("latent_shape") != [None, 1920, 6144]
        or sequence.get("central_bin_resolution_bp") != 32
        or sequence.get("sequence_only_identity_required") is not True
        or sequence.get("geometry_change_requires_new_parent_bound_spec") is not True
    ):
        raise ContextBorzoiReadinessError("provisional sequence-latent interface differs")

    context = contract.get("context_interface", {})
    if (
        context.get("status") != "SCHEMA_READY_ENCODER_UNSELECTED"
        or "observed_query_atac" not in context.get("forbidden_query_sources", [])
        or "held_atac_target" not in context.get("forbidden_query_sources", [])
        or "sealed_label_or_outcome" not in context.get("forbidden_query_sources", [])
        or set(context.get("missingness_states", [])) != MISSING_STATES
        or context.get("missing_context_value") is not None
        or context.get("missing_context_mask") is not False
        or context.get("missing_as_numeric_zero_allowed") is not False
        or context.get("missing_context_effect") != "exact_neutral_identity"
        or context.get("mapper_fit_scope") != "outer_training_partition_only"
    ):
        raise ContextBorzoiReadinessError("context or missingness interface differs")

    conditioning = contract.get("conditioning_interface", {})
    if (
        conditioning.get("film", {}).get("gate")
        != "explicit_context_observed_mask"
        or conditioning.get("film", {}).get("missing_context_behavior")
        != "latent_unchanged"
        or conditioning.get("lora", {}).get("exact_module_names")
        != "UNRESOLVED_UNTIL_PARENT_SELECTION"
        or conditioning.get("lora", {}).get("neutral_initialization") != "zero_delta"
        or conditioning.get("lora", {}).get("missing_context_behavior")
        != "delta_exactly_zero"
    ):
        raise ContextBorzoiReadinessError("conditioning interface differs")

    heads = contract.get("head_interface", {})
    tasks = contract.get("family_native_task_boundaries", {})
    if (
        heads.get("allowed_heads") != ["atac_profile", "atac_count", "rna_coverage"]
        or heads.get("forbidden_heads") != ["masld_histone"]
        or heads.get("profile_and_count_factorized") is not True
        or tasks.get("rna_conditioned_atac", {}).get("query_atac_allowed") is not False
        or tasks.get("rna_conditioned_atac", {}).get("biological_unit") != "donor"
        or tasks.get("observed_multiome", {}).get("eligible_interface") is not False
        or tasks.get("observed_multiome", {}).get("results_may_support_universal_claim")
        is not False
        or tasks.get("variant_to_regulation", {}).get("native_outputs")
        != ["accessibility_delta", "gene_expression_delta"]
    ):
        raise ContextBorzoiReadinessError("head or family-native task boundary differs")

    components = contract.get("component_reconciliation", {})
    borzoi = components.get("borzoi", {})
    borzoi_authority = _resolve_file(root, borzoi.get("official_authority_path", ""))
    if (
        _digest(borzoi_authority) != borzoi.get("official_authority_sha256")
        or borzoi.get("selected_as_parent") is not False
        or borzoi.get("weights_reused_by_this_contract") is not False
        or "native TensorFlow numeric parity" not in borzoi.get("not_proven", [])
        or "official weight terms" not in borzoi.get("not_proven", [])
    ):
        raise ContextBorzoiReadinessError("Borzoi authority or boundary differs")
    borzoi_artifacts = {
        "static_audit": _require_frozen_artifact(
            root,
            borzoi["static_converted_port_audit_artifact"],
            borzoi["static_converted_port_audit_artifacts_sha256"],
            artifact_class="borzoi_grelu_converted_port_static_audit",
        ),
        "semantic_mapping": _require_frozen_artifact(
            root,
            borzoi["semantic_mapping_artifact"],
            borzoi["semantic_mapping_artifacts_sha256"],
            artifact_class="borzoi_grelu_to_local_port_semantic_mapping_probe",
        ),
        "semantic_diagnostic": _require_frozen_artifact(
            root,
            borzoi["semantic_diagnostic_artifact"],
            borzoi["semantic_diagnostic_artifacts_sha256"],
            artifact_class="borzoi_grelu_local_semantics_diagnostic",
        ),
        "disabled_gpu_queue": _require_frozen_artifact(
            root,
            borzoi["disabled_full_window_gpu_queue_artifact"],
            borzoi["disabled_full_window_gpu_queue_artifacts_sha256"],
            artifact_class="model_training_603_disabled_gpu_queue_admission",
        ),
    }
    borzoi_meta = {name: value["metadata"] for name, value in borzoi_artifacts.items()}
    if (
        borzoi_meta["static_audit"].get("model_forward_executed") is not False
        or borzoi_meta["static_audit"].get("native_borzoi_substitute") is not False
        or borzoi_meta["semantic_mapping"].get("all_four_strict_loads_passed")
        is not True
        or borzoi_meta["semantic_mapping"].get("native_borzoi_substitute") is not False
        or borzoi_meta["semantic_diagnostic"].get("native_numeric_parity_established")
        is not False
        or borzoi_meta["semantic_diagnostic"].get("local_port_execution_gate_opened")
        is not False
        or borzoi_meta["disabled_gpu_queue"].get("queue_enabled") is not False
        or borzoi_meta["disabled_gpu_queue"].get("gpu_job_submitted") is not False
    ):
        raise ContextBorzoiReadinessError("Borzoi evidence opened a prohibited gate")

    corgi = components.get("corgi", {})
    corgi_artifacts = {
        "film_contract": _require_frozen_artifact(
            root,
            corgi["film_contract_artifact"],
            corgi["film_contract_artifacts_sha256"],
            artifact_class="corgi_regular_film_plus_head_prospective_contract",
        ),
        "adapter_unit": _require_frozen_artifact(
            root,
            corgi["adapter_unit_artifact"],
            corgi["adapter_unit_artifacts_sha256"],
            artifact_class="corgi_regular_film_plus_head_adapter_unit_validation_v2",
        ),
        "outcome_blind_fixture": _require_frozen_artifact(
            root,
            corgi["outcome_blind_fixture_artifact"],
            corgi["outcome_blind_fixture_artifacts_sha256"],
            artifact_class="corgi_regular_film_plus_head_outcome_blind_fixture_v2",
        ),
    }
    corgi_meta = {name: value["metadata"] for name, value in corgi_artifacts.items()}
    if (
        corgi_meta["film_contract"].get("training_execution_authorized") is not False
        or corgi_meta["film_contract"].get("development_outcomes_read") is not False
        or corgi_meta["adapter_unit"].get("checkpoint_loaded") is not False
        or corgi_meta["adapter_unit"].get("training_performed") is not False
        or corgi_meta["outcome_blind_fixture"].get("development_atac_outcomes_read")
        is not False
        or corgi_meta["outcome_blind_fixture"].get("model_fit_performed") is not False
        or corgi.get("corgi_checkpoint_selected_or_reused") is not False
        or corgi.get("corgi_film_weights_reused") is not False
    ):
        raise ContextBorzoiReadinessError("Corgi interface evidence differs")

    scooby = components.get("scooby", {})
    scooby_authority = _load_json(
        _resolve_file(root, scooby.get("checkpoint_authority_path", ""))
    )
    if (
        _digest(_resolve_file(root, scooby["checkpoint_authority_path"]))
        != scooby.get("checkpoint_authority_sha256")
        or scooby_authority.get("architecture_contract", {}).get(
            "sequence_embedding_dim"
        )
        != 1920
        or scooby_authority.get("architecture_contract", {}).get(
            "output_central_bins"
        )
        != 6144
        or scooby.get("checkpoint_or_code_reused") is not False
        or scooby.get("registered_rna_conditioned_atac_models") != []
    ):
        raise ContextBorzoiReadinessError("Scooby authority or interface differs")
    scooby_artifact = _require_frozen_artifact(
        root,
        scooby["released_lane_audit_artifact"],
        scooby["released_lane_audit_artifacts_sha256"],
        artifact_class="scooby_released_lane_audit",
    )
    scooby_meta = scooby_artifact["metadata"]
    if (
        scooby_meta.get("checkpoint_forward_allowed") != []
        or scooby_meta.get("rna_conditioned_atac_models") != []
        or scooby_meta.get("checkpoint_payloads_downloaded") is not False
        or scooby_meta.get("project_data_read") is not False
        or scooby_meta.get("sealed_outcomes_read") is not False
    ):
        raise ContextBorzoiReadinessError("Scooby released-lane boundary differs")

    execution = contract.get("execution_disposition", {})
    if execution.get("cpu_metadata_and_unit_validation_authorized") is not True:
        raise ContextBorzoiReadinessError("CPU readiness validation is not authorized")
    _assert_false(execution, EXECUTION_FLAGS)
    firewall = contract.get("outcome_firewall", {})
    if any(value is not False for value in firewall.values()):
        raise ContextBorzoiReadinessError("outcome firewall differs")
    if len(contract.get("remaining_activation_gates", [])) != 8:
        raise ContextBorzoiReadinessError("activation-gate roster differs")

    return {
        "schema_version": "masld-bench-context-borzoi-component-readiness-receipt-v1",
        "status": "pass_interface_contract_build_blocked",
        "contract_sha256": _digest(contract_path),
        "verified_frozen_artifacts": 8,
        "provisional_input_shape": [None, 4, 524288],
        "provisional_latent_shape": [None, 1920, 6144],
        "context_missing_as_zero_allowed": False,
        "query_atac_allowed_for_rna_conditioned_task": False,
        "borzoi_selected": False,
        "corgi_selected_or_reused": False,
        "scooby_selected_or_reused": False,
        "component_interface_schema_ready": True,
        "complementarity_trigger_passed": False,
        "component_selection_locked": False,
        "architecture_built": False,
        "checkpoint_opened_by_audit": False,
        "biological_data_read": False,
        "development_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_labels_or_outcomes_read": False,
        "training_authorized": False,
        "prediction_authorized": False,
        "gpu_bundle_registration_authorized": False,
        "gpu_job_submission_authorized": False,
    }


def main() -> None:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_root)
    parser.add_argument(
        "--contract",
        type=Path,
        default=default_root
        / "config/artifacts/models/context_borzoi/component_interface_readiness_20260825.json",
    )
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    receipt = audit_readiness(arguments.project_root, arguments.contract)
    payload = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if arguments.output is not None:
        if arguments.output.exists() or arguments.output.is_symlink():
            raise ContextBorzoiReadinessError("refusing to overwrite readiness receipt")
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(payload, encoding="utf-8")
    print(payload, end="")


if __name__ == "__main__":
    main()
