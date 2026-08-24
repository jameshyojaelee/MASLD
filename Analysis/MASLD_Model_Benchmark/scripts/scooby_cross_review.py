#!/usr/bin/env python3
"""Fail-closed cross-review for the three released scooby checkpoints."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any


class ScoobyCrossReviewError(ValueError):
    """Raised when a scooby authority or task boundary differs."""


MODEL_CONTRACTS: dict[str, dict[str, Any]] = {
    "scooby_onek1k": {
        "finding_id": "onek1k",
        "context_dim": 10,
        "tracks": ["rna_plus", "rna_minus"],
        "topology": "rna_only_pbmc_scpoli",
        "observed_atac_in_context": False,
        "requires_observed_target_context": False,
        "domain": "blood_PBMC",
        "license": "declared_MIT_with_unresolved_inherited_Borzoi_authority",
        "checkpoint_download_allowed": True,
        "checkpoint_execution_allowed": False,
        "terminal_disposition": (
            "metadata_and_bounded_safetensors_header_allowed; forward_blocked_"
            "missing_exact_10D_encoder_runtime_reference_parity_and_inherited_terms"
        ),
    },
    "scooby_epicardioids": {
        "finding_id": "epicardioids",
        "context_dim": 50,
        "tracks": ["rna_plus", "rna_minus", "atac_insertion"],
        "topology": "same_study_unpaired_RNA_ATAC_scGLUE_pseudomatches",
        "observed_atac_in_context": True,
        "requires_observed_target_context": True,
        "domain": "cardiac_epicardioid",
        "license": "declared_MIT_with_unresolved_inherited_Borzoi_authority",
        "checkpoint_download_allowed": True,
        "checkpoint_execution_allowed": False,
        "terminal_disposition": (
            "metadata_and_bounded_safetensors_header_allowed; forward_blocked_"
            "missing_exact_50D_encoder_runtime_reference_parity_and_inherited_terms"
        ),
    },
    "scooby_neurips": {
        "finding_id": "neurips",
        "context_dim": 14,
        "tracks": ["rna_plus", "rna_minus", "atac_insertion"],
        "topology": "same_nucleus_RNA_ATAC_Poisson_MultiVI",
        "observed_atac_in_context": True,
        "requires_observed_target_context": True,
        "domain": "bone_marrow_hematopoietic",
        "license": "undeclared_weight_terms_and_unresolved_inherited_Borzoi_authority",
        "checkpoint_download_allowed": False,
        "checkpoint_execution_allowed": False,
        "terminal_disposition": (
            "metadata_and_pointer_only; checkpoint_object_and_forward_blocked_"
            "undeclared_terms_missing_exact_14D_encoder_and_runtime_reference_parity"
        ),
    },
}


FORBIDDEN_OUTCOME_TOKENS = (
    "gse289173_label_values",
    "gse289173_outcome_values",
    "gse289173_eqtl_table",
    "gse289173_ieqtl_table",
    "sealed_truth",
    "observed_y",
)


def _sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 1_000_000:
        raise ScoobyCrossReviewError(f"authority path differs: {path}")
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ScoobyCrossReviewError(f"authority JSON differs: {path}") from exc
    if not isinstance(payload, dict):
        raise ScoobyCrossReviewError(f"authority JSON is not an object: {path}")
    return payload


def _model_blocks(path: Path) -> dict[str, str]:
    if path.is_symlink() or not path.is_file() or path.stat().st_size > 1_000_000:
        raise ScoobyCrossReviewError(f"model registry path differs: {path}")
    text = path.read_text(encoding="utf-8")
    blocks: dict[str, str] = {}
    for block in text.split("[[models]]")[1:]:
        match = re.search(r'^model_id\s*=\s*"([^"]+)"\s*$', block, re.MULTILINE)
        if match:
            blocks[match.group(1)] = block
    return blocks


def _quoted_value(block: str, key: str) -> str:
    match = re.search(
        rf'^{re.escape(key)}\s*=\s*"([^"]+)"\s*$', block, re.MULTILINE
    )
    if not match:
        raise ScoobyCrossReviewError(f"missing registry field: {key}")
    return match.group(1)


def _bool_value(block: str, key: str) -> bool:
    match = re.search(
        rf"^{re.escape(key)}\s*=\s*(true|false)\s*$", block, re.MULTILINE
    )
    if not match:
        raise ScoobyCrossReviewError(f"missing registry boolean: {key}")
    return match.group(1) == "true"


def _string_list(block: str, key: str) -> list[str]:
    match = re.search(
        rf"^{re.escape(key)}\s*=\s*(\[[^\n]*\])\s*$", block, re.MULTILINE
    )
    if not match:
        raise ScoobyCrossReviewError(f"missing registry list: {key}")
    try:
        value = json.loads(match.group(1))
    except json.JSONDecodeError as exc:
        raise ScoobyCrossReviewError(f"registry list differs: {key}") from exc
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ScoobyCrossReviewError(f"registry list differs: {key}")
    return value


def _context_dim(cell_decoder: str) -> int:
    match = re.search(r"\b(\d+)-dimensional\b", cell_decoder)
    if not match:
        raise ScoobyCrossReviewError("cell-context dimension is missing")
    return int(match.group(1))


def _review_model(
    model_id: str,
    checkpoint: dict[str, Any],
    crosswalk: dict[str, Any],
    exposure: dict[str, Any],
    registry_block: str,
    capability_block: str,
) -> dict[str, Any]:
    expected = MODEL_CONTRACTS[model_id]
    finding = exposure.get("checkpoint_findings", {}).get(expected["finding_id"], {})
    if (
        checkpoint.get("schema_version")
        != "masld-bench-upstream-checkpoint-preflight-v1"
        or crosswalk.get("schema_version")
        != "masld-bench-development-crosswalk-v1"
        or exposure.get("schema_version")
        != "masld-bench-checkpoint-exposure-audit-v1"
        or finding.get("exposure_state") != "target_label_unexposed"
        or crosswalk.get("findings", {}).get("gse289173", {}).get("exposure_state")
        != "target_label_unexposed"
    ):
        raise ScoobyCrossReviewError(f"exposure contract differs: {model_id}")

    architecture = checkpoint.get("architecture_contract", {})
    output = checkpoint.get("output_contract", {})
    inputs = checkpoint.get("input_contract", {})
    if (
        _context_dim(str(architecture.get("cell_decoder", "")))
        != expected["context_dim"]
        or architecture.get("decoder_tracks") != len(expected["tracks"])
        or output.get("profiles") != expected["tracks"]
        or "NOT_BUNDLED" not in str(inputs.get("cell_context_encoder_status", ""))
        or "524288" not in str(inputs.get("sequence", ""))
    ):
        raise ScoobyCrossReviewError(f"input/output contract differs: {model_id}")

    observed_atac_text = " ".join(
        str(value)
        for key, value in inputs.items()
        if key in {"cell_context", "topology", "warning"}
    ).lower()
    observed_atac_detected = "atac" in observed_atac_text
    if observed_atac_detected is not expected["observed_atac_in_context"]:
        raise ScoobyCrossReviewError(f"observed-ATAC topology differs: {model_id}")
    rna_atac_fit = str(
        checkpoint.get("task_fit", {}).get("rna_conditioned_atac", "")
    ).lower()
    if rna_atac_fit.startswith("supported"):
        raise ScoobyCrossReviewError(f"RNA-conditioned ATAC gate differs: {model_id}")

    if (
        _string_list(registry_block, "supported_tasks") != ["variant_to_regulation"]
        or not _bool_value(registry_block, "admission_blocking")
        or _quoted_value(registry_block, "exposure_status")
        != "target_label_unexposed"
        or _bool_value(capability_block, "primary_eligible")
        or _bool_value(capability_block, "requires_observed_target_context")
        is not expected["requires_observed_target_context"]
    ):
        raise ScoobyCrossReviewError(f"task registry differs: {model_id}")

    weight_status = str(checkpoint.get("license_audit", {}).get("weight_status", ""))
    registry_license = _quoted_value(registry_block, "license_status")
    if model_id == "scooby_neurips":
        if "UNDECLARED" not in weight_status or "UNDECLARED" not in registry_license:
            raise ScoobyCrossReviewError("NeurIPS weight terms are no longer undeclared")
    elif (
        "UNRESOLVED_INHERITED_BORZOI" not in weight_status
        or "UNRESOLVED" not in registry_license
    ):
        raise ScoobyCrossReviewError(f"inherited weight terms differ: {model_id}")

    sealed = str(exposure.get("sealed_champion_eligibility", ""))
    if "ineligible" not in sealed.lower():
        raise ScoobyCrossReviewError(f"sealed champion gate differs: {model_id}")

    return {
        "model_id": model_id,
        "context_dim": expected["context_dim"],
        "decoder_tracks": expected["tracks"],
        "native_context_topology": expected["topology"],
        "native_domain": expected["domain"],
        "observed_atac_in_checkpoint_context": expected["observed_atac_in_context"],
        "observed_query_atac_comparator": expected["observed_atac_in_context"],
        "rna_conditioned_atac_eligible": False,
        "registered_tasks": ["variant_to_regulation"],
        "gse289173_target_label_exposure": "target_label_unexposed",
        "sealed_masld_champion_eligible": False,
        "license_disposition": expected["license"],
        "checkpoint_download_allowed": expected["checkpoint_download_allowed"],
        "checkpoint_deserialization_allowed": False,
        "checkpoint_execution_allowed": expected["checkpoint_execution_allowed"],
        "terminal_disposition": expected["terminal_disposition"],
    }


def _validate_outcome_firewall(firewall: dict[str, Any]) -> None:
    lowered = json.dumps(firewall, sort_keys=True).lower()
    if any(token in lowered for token in FORBIDDEN_OUTCOME_TOKENS):
        raise ScoobyCrossReviewError("outcome firewall contains a forbidden key")
    expected = {
        "review_inputs": "authority_metadata_only",
        "project_data_read": False,
        "sealed_features_read": False,
        "sealed_labels_or_outcomes_read": False,
        "gse289173_eqtl_or_ieqtl_read": False,
        "gse296875_outcomes_read": False,
        "frozen_117_hotspot_programs_read": False,
        "checkpoint_bytes_read": False,
        "observed_query_atac_results_must_be_separate": True,
        "held_atac_forbidden_in_rna_conditioned_track": True,
        "prediction_hash_committed_before_sealed_label_join": True,
    }
    if firewall != expected:
        raise ScoobyCrossReviewError("outcome firewall differs")


def review(project_root: Path, output: Path) -> dict[str, Any]:
    root = project_root.resolve(strict=True)
    if output.exists() or output.is_symlink():
        raise ScoobyCrossReviewError("output already exists")
    config = root / "config"
    registry_path = config / "models" / "regulatory_sequence.toml"
    capability_path = config / "evaluation" / "variant_to_regulation_capabilities.toml"
    registry_blocks = _model_blocks(registry_path)
    capability_blocks = _model_blocks(capability_path)
    receipts: dict[str, dict[str, Any]] = {}
    authority_files: dict[str, dict[str, Any]] = {
        "model_registry": {
            "path": str(registry_path.relative_to(root)),
            "sha256": _sha256(registry_path),
        },
        "variant_capabilities": {
            "path": str(capability_path.relative_to(root)),
            "sha256": _sha256(capability_path),
        },
    }
    for model_id, expected in MODEL_CONTRACTS.items():
        bundle = config / "artifacts" / "models" / model_id
        paths = {
            name: bundle / name
            for name in (
                "checkpoints.json",
                "development_crosswalk.json",
                "exposure_audit.json",
            )
        }
        checkpoint = _load_json(paths["checkpoints.json"])
        crosswalk = _load_json(paths["development_crosswalk.json"])
        exposure = _load_json(paths["exposure_audit.json"])
        if exposure.get("development_crosswalk_record", {}).get("sha256") != _sha256(
            paths["development_crosswalk.json"]
        ):
            raise ScoobyCrossReviewError(f"crosswalk binding differs: {model_id}")
        if model_id not in registry_blocks or model_id not in capability_blocks:
            raise ScoobyCrossReviewError(f"registry binding is missing: {model_id}")
        receipts[model_id] = _review_model(
            model_id,
            checkpoint,
            crosswalk,
            exposure,
            registry_blocks[model_id],
            capability_blocks[model_id],
        )
        for name, path in paths.items():
            authority_files[f"{model_id}/{name}"] = {
                "path": str(path.relative_to(root)),
                "sha256": _sha256(path),
            }

    outcome_firewall = {
        "review_inputs": "authority_metadata_only",
        "project_data_read": False,
        "sealed_features_read": False,
        "sealed_labels_or_outcomes_read": False,
        "gse289173_eqtl_or_ieqtl_read": False,
        "gse296875_outcomes_read": False,
        "frozen_117_hotspot_programs_read": False,
        "checkpoint_bytes_read": False,
        "observed_query_atac_results_must_be_separate": True,
        "held_atac_forbidden_in_rna_conditioned_track": True,
        "prediction_hash_committed_before_sealed_label_join": True,
    }
    _validate_outcome_firewall(outcome_firewall)
    receipt = {
        "schema_version": "masld-bench-scooby-cross-review-v1",
        "status": "pass",
        "review_order": [
            "task_topology",
            "license",
            "checkpoint_exposure",
            "outcome_firewall",
            "production_disposition",
        ],
        "authority_files": authority_files,
        "models": receipts,
        "outcome_firewall": outcome_firewall,
        "rna_conditioned_atac_models": [],
        "observed_query_atac_models": [
            "scooby_epicardioids",
            "scooby_neurips",
        ],
        "production_allowed": {
            "scooby_onek1k": "metadata_and_bounded_header_only",
            "scooby_epicardioids": "metadata_and_bounded_header_only",
            "scooby_neurips": "metadata_and_pointer_only",
        },
        "checkpoint_forward_allowed": [],
        "sealed_inference_allowed": [],
        "open_champion_eligible": [],
        "terminal_disposition": (
            "cross_review_passed; bounded_admission_only; all_checkpoint_forwards_"
            "blocked_pending_model_specific_terms_context_runtime_and_reference_gates"
        ),
    }
    output.mkdir(parents=True, mode=0o750)
    (output / "scooby_cross_review_receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    review(arguments.project_root, arguments.output)


if __name__ == "__main__":
    main()
