#!/usr/bin/env python3
"""Build the biology-free context-Borzoi RNA adapter fixture on CPU only."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping

from masld_bench.adapters.context_borzoi_fixture import (
    AdapterFixtureSpec,
    ContextBorzoiFixtureError,
    REQUIRED_ARMS,
    RNAContextMask,
    SyntheticContextRow,
    derange_context_assignments,
    export_fixture_bundle,
    run_adapter_arm,
)
from masld_bench.artifacts import ArtifactError, verify_frozen_tree
from masld_bench.context_preconditions import require_precondition_only_campaign
from masld_bench.topology import MISSING_STATES


CONTRACT_SCHEMA = "masld-bench-context-borzoi-rna-adapter-fixture-contract-v1"
CONTRACT_STATUS = "SYNTHETIC_FIXTURE_ALLOWED_ARCHITECTURE_BLOCKED"


class FixtureBuildError(RuntimeError):
    """Raised when a frozen authority or fixture boundary differs."""


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise FixtureBuildError(f"JSON authority is not an object: {path}")
    return value


def _resolve(root: Path, relative: str, *, directory: bool = False) -> Path:
    if not relative or Path(relative).is_absolute():
        raise FixtureBuildError("authority path must be relative")
    path = root / relative
    if path.is_symlink():
        raise FixtureBuildError(f"authority cannot be a symlink: {relative}")
    try:
        resolved = path.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise FixtureBuildError(f"authority escapes or is absent: {relative}") from error
    if directory != resolved.is_dir():
        raise FixtureBuildError(f"authority type differs: {relative}")
    return resolved


def validate_contract(root: Path, contract: Mapping[str, object]) -> AdapterFixtureSpec:
    if contract.get("schema_version") != CONTRACT_SCHEMA or contract.get("status") != CONTRACT_STATUS:
        raise FixtureBuildError("fixture contract schema or status differs")
    component = contract.get("component_authority")
    if not isinstance(component, dict):
        raise FixtureBuildError("component authority missing")
    artifact = _resolve(root, str(component["artifact_path"]), directory=True)
    try:
        manifest = verify_frozen_tree(artifact)
    except ArtifactError as error:
        raise FixtureBuildError("component authority is not a valid frozen tree") from error
    if (
        _sha256(artifact / "ARTIFACTS.json") != component["artifacts_sha256"]
        or manifest.get("metadata", {}).get("artifact_class") != component["artifact_class"]
    ):
        raise FixtureBuildError("component artifact identity differs")
    receipt_path = artifact / str(component["receipt_path"])
    if _sha256(receipt_path) != component["receipt_sha256"]:
        raise FixtureBuildError("component receipt identity differs")
    receipt = _load_json(receipt_path)
    if (
        receipt.get("status") != component["required_status"]
        or receipt.get("complementarity_trigger_passed") is not False
        or receipt.get("component_selection_locked") is not False
        or receipt.get("architecture_built") is not False
        or receipt.get("checkpoint_opened_by_audit") is not False
        or receipt.get("training_authorized") is not False
        or receipt.get("prediction_authorized") is not False
        or receipt.get("gpu_job_submission_authorized") is not False
    ):
        raise FixtureBuildError("component receipt opened a closed gate")

    interface_path = _resolve(root, str(component["interface_path"]))
    if _sha256(interface_path) != component["interface_sha256"]:
        raise FixtureBuildError("component interface identity differs")
    interface = _load_json(interface_path)
    slot = interface.get("conditional_slot", {})
    context = interface.get("context_interface", {})
    conditioning = interface.get("conditioning_interface", {})
    tasks = interface.get("family_native_task_boundaries", {})
    if (
        interface.get("status") != "INTERFACE_CONTRACT_READY_BUILD_BLOCKED"
        or slot.get("complementarity_trigger_state") != "NOT_EVALUATED"
        or slot.get("component_selection_lock") != "ABSENT"
        or slot.get("architecture_built") is not False
        or context.get("missing_as_numeric_zero_allowed") is not False
        or context.get("missing_context_effect") != "exact_neutral_identity"
        or "observed_query_atac" not in context.get("forbidden_query_sources", [])
        or conditioning.get("film", {}).get("gate") != "explicit_context_observed_mask"
        or conditioning.get("lora", {}).get("exact_module_names")
        != "UNRESOLVED_UNTIL_PARENT_SELECTION"
        or tasks.get("observed_multiome", {}).get("eligible_interface") is not False
    ):
        raise FixtureBuildError("component interface boundary differs")

    campaign_spec = contract.get("campaign_authority")
    if not isinstance(campaign_spec, dict):
        raise FixtureBuildError("campaign authority missing")
    campaign_path = _resolve(root, str(campaign_spec["path"]))
    if _sha256(campaign_path) != campaign_spec["sha256"]:
        raise FixtureBuildError("campaign identity differs")
    with campaign_path.open("rb") as handle:
        campaign = tomllib.load(handle)
    require_precondition_only_campaign(campaign)
    if (
        campaign.get("status") != campaign_spec["required_status"]
        or campaign.get("conditional_model_spec_path")
        != campaign_spec["conditional_model_spec_path"]
        or campaign.get("submit_enabled") is not False
    ):
        raise FixtureBuildError("conditional campaign is no longer blocked")

    fixture = contract.get("fixture_interface")
    if not isinstance(fixture, dict):
        raise FixtureBuildError("fixture interface missing")
    if (
        fixture.get("source_kind") != "synthetic_fixture_only"
        or fixture.get("task_id") != "rna_conditioned_atac"
        or fixture.get("arms") != list(REQUIRED_ARMS)
        or set(fixture.get("missingness_states", [])) != set(MISSING_STATES)
        or fixture.get("missing_value") is not None
        or fixture.get("missing_mask") is not False
        or fixture.get("missing_as_numeric_zero_allowed") is not False
        or fixture.get("missing_context_effect") != "exact_sequence_only_identity"
        or fixture.get("remaining_parent_bound_ablations")
        != ["reverse_complement", "matched_shuffle"]
    ):
        raise FixtureBuildError("fixture mask or ablation contract differs")

    family = contract.get("family_native_boundary")
    if not isinstance(family, dict) or any(
        family.get(key) is not False
        for key in (
            "observed_atac_input_allowed",
            "held_atac_target_read",
            "observed_multiome_task_modified",
            "histology_or_disease_label_allowed",
            "sealed_input_allowed",
            "biological_unit_claimed",
        )
    ):
        raise FixtureBuildError("family-native task boundary differs")
    export = contract.get("export_contract")
    if (
        not isinstance(export, dict)
        or export.get("format") != "canonical_json_utf8"
        or any(export.get(key) is not False for key in ("pickle_allowed", "tensor_state_allowed", "model_weights_allowed", "checkpoint_allowed"))
    ):
        raise FixtureBuildError("non-pickle export contract differs")
    execution = contract.get("execution_disposition")
    if not isinstance(execution, dict) or execution.get("cpu_unit_fixture_authorized") is not True:
        raise FixtureBuildError("CPU fixture is not authorized")
    if any(value is not False for key, value in execution.items() if key != "cpu_unit_fixture_authorized"):
        raise FixtureBuildError("fixture execution gate opened")
    gate = contract.get("conditional_gate")
    if (
        not isinstance(gate, dict)
        or gate.get("state") != "NOT_EVALUATED"
        or any(value is not False for key, value in gate.items() if key.endswith("passed") or key.endswith("complete") or key.endswith("frozen") or key.endswith("authorized"))
    ):
        raise FixtureBuildError("conditional trigger is no longer closed")

    spec = AdapterFixtureSpec(
        feature_width=int(fixture["feature_width"]),
        lora_rank=int(fixture["lora_rank"]),
        lora_alpha=float(fixture["lora_alpha"]),
        permutation_seed=int(fixture["permutation_seed"]),
    )
    spec.validate()
    return spec


def _context(values: tuple[float | None, ...]) -> RNAContextMask:
    return RNAContextMask(
        values=values,
        observed_mask=tuple(value is not None for value in values),
        states=tuple(
            "observed" if value is not None else "structurally_missing"
            for value in values
        ),
    )


def _row(row_id: str, values: tuple[float | None, ...]) -> SyntheticContextRow:
    return SyntheticContextRow(
        row_id=row_id,
        outer_partition="synthetic-outer-0",
        dataset_id="synthetic-dataset",
        assay="synthetic-rna",
        lineage="synthetic-lineage",
        topology="synthetic-unpaired",
        rna=_context(values),
    )


def build_fixture(root: Path, contract: Mapping[str, object], output: Path) -> dict[str, object]:
    spec = validate_contract(root, contract)
    row_a = _row("synthetic-row-a", (1.0, 2.0, 3.0, 4.0))
    row_b = _row("synthetic-row-b", (4.0, 3.0, 2.0, 1.0))
    row_missing = _row("synthetic-row-missing", (1.0, None, 3.0, 4.0))
    assignments = derange_context_assignments(
        (row_a, row_b), width=spec.feature_width, seed=spec.permutation_seed
    )
    sequence = (
        (0.25, -0.50, 0.75, 1.00),
        (1.25, 1.50, -1.75, 2.00),
        (2.25, -2.50, 2.75, -3.00),
    )
    outputs = {
        "conditioned": run_adapter_arm(spec, sequence, row_a, arm="conditioned"),
        "sequence_only": run_adapter_arm(spec, sequence, row_a, arm="sequence_only"),
        "trans_only": run_adapter_arm(spec, sequence, row_a, arm="trans_only"),
        "permuted_context": run_adapter_arm(
            spec,
            sequence,
            row_a,
            arm="permuted_context",
            permuted_source=assignments[row_a.row_id],
        ),
        "missing_context": run_adapter_arm(
            spec, sequence, row_missing, arm="conditioned"
        ),
    }
    if outputs["sequence_only"].values != sequence:
        raise FixtureBuildError("sequence-only arm is not exact identity")
    if outputs["missing_context"].values != outputs["sequence_only"].values:
        raise FixtureBuildError("missing context did not route to sequence-only identity")
    if outputs["conditioned"].values == outputs["permuted_context"].values:
        raise FixtureBuildError("permuted context did not change the synthetic fixture")
    alternate_sequence = tuple(tuple(value + 9.0 for value in row) for row in sequence)
    alternate_trans = run_adapter_arm(
        spec, alternate_sequence, row_a, arm="trans_only"
    )
    if outputs["trans_only"].values != alternate_trans.values:
        raise FixtureBuildError("trans-only arm depends on sequence values")
    export_manifest = export_fixture_bundle(
        output / "fixture_export", spec, outputs, assignments
    )
    receipt = {
        "schema_version": "masld-bench-context-borzoi-rna-adapter-fixture-receipt-v1",
        "status": "pass_synthetic_fixture_architecture_blocked",
        "component_authority_verified": True,
        "complementarity_trigger_passed": False,
        "component_selection_locked": False,
        "parent_component_selected": False,
        "architecture_built": False,
        "checkpoint_opened": False,
        "model_weights_present": False,
        "training_performed": False,
        "biological_prediction_performed": False,
        "synthetic_forward_fixture_executed": True,
        "biological_data_read": False,
        "development_outcomes_read": False,
        "sealed_features_or_outcomes_read": False,
        "observed_atac_read": False,
        "observed_multiome_task_modified": False,
        "gpu_executed": False,
        "arms_executed": sorted(outputs),
        "sequence_only_exact_identity": True,
        "missing_context_exact_sequence_only_identity": True,
        "trans_only_sequence_invariant": True,
        "permuted_context_deranged_within_stratum": True,
        "explicit_missing_masks": True,
        "export_format": "canonical_json_utf8",
        "pickle_present": export_manifest["pickle_present"],
        "tensor_state_present": export_manifest["tensor_state_present"],
        "conditional_model_authorized": False,
        "gpu_job_submission_authorized": False,
        "remaining_parent_bound_ablations": ["reverse_complement", "matched_shuffle"],
    }
    (output / "receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve(strict=True)
    contract_path = args.contract.resolve(strict=True)
    try:
        contract_path.relative_to(root)
    except ValueError as error:
        raise FixtureBuildError("fixture contract escapes benchmark root") from error
    contract = _load_json(contract_path)
    try:
        receipt = build_fixture(root, contract, args.output)
    except ContextBorzoiFixtureError as error:
        raise FixtureBuildError(str(error)) from error
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
