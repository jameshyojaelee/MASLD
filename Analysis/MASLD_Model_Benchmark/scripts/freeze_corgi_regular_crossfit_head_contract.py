#!/usr/bin/env python3
"""Validate a prospective, outcome-free Regular Corgi profile-head requirement."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping, Sequence


SCHEMA_VERSION = "masld-bench-corgi-crossfit-profile-head-contract-v1"
RECEIPT_SCHEMA_VERSION = "masld-bench-corgi-crossfit-profile-head-preflight-v1"
FOLDS = tuple(range(5))
CONTEXT_CONTROLS = {
    "training_lineage_mean_released_rank",
    "nearest_training_released_rank",
    "shuffled_valid_released_rank",
}


class CorgiHeadContractError(RuntimeError):
    """Raised when the prospective head or its read-only inputs drift."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def _load_json(path: Path) -> Mapping[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise CorgiHeadContractError(f"JSON authority differs: {path}")
    return value


def _artifact_metadata(root: Path) -> Mapping[str, Any]:
    value = _load_json(root / "ARTIFACTS.json")
    if value.get("schema_version") != "masld-bench-artifacts-v1":
        raise CorgiHeadContractError(f"artifact schema differs: {root}")
    metadata = value.get("metadata")
    if not isinstance(metadata, dict):
        raise CorgiHeadContractError(f"artifact metadata differs: {root}")
    return metadata


def _check_source(root: Path, relative: str, expected: str) -> Path:
    if not isinstance(relative, str) or not isinstance(expected, str) or len(expected) != 64:
        raise CorgiHeadContractError("source binding differs")
    path = (root / relative).resolve(strict=True)
    try:
        path.relative_to(root.resolve(strict=True))
    except ValueError as error:
        raise CorgiHeadContractError("source binding escapes project root") from error
    target = path / "ARTIFACTS.json" if path.is_dir() else path
    if digest(target) != expected:
        raise CorgiHeadContractError(f"source hash differs: {relative}")
    return path


def validate_fold_plan(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    if len(rows) != len(FOLDS):
        raise CorgiHeadContractError("cross-fit fold count differs")
    normalized: list[dict[str, Any]] = []
    seen: set[int] = set()
    for row in rows:
        evaluation_fold = int(row.get("evaluation_fold", -1))
        evaluation_base = int(row.get("evaluation_base_outer_fold", -1))
        fit_valid = [int(value) for value in row.get("fit_valid_folds", [])]
        fit_base = [int(value) for value in row.get("fit_base_outer_folds", [])]
        if evaluation_fold not in FOLDS or evaluation_fold in seen:
            raise CorgiHeadContractError("evaluation fold identity differs")
        seen.add(evaluation_fold)
        expected_valid = [fold for fold in FOLDS if fold != evaluation_fold]
        expected_base = [(fold - 1) % len(FOLDS) for fold in expected_valid]
        if (
            evaluation_base != (evaluation_fold - 1) % len(FOLDS)
            or fit_valid != expected_valid
            or fit_base != expected_base
            or evaluation_fold in fit_valid
            or evaluation_base in fit_base
        ):
            raise CorgiHeadContractError("donor/genomic cross-fit firewall differs")
        normalized.append(
            {
                "evaluation_fold": evaluation_fold,
                "evaluation_base_outer_fold": evaluation_base,
                "fit_valid_folds": fit_valid,
                "fit_base_outer_folds": fit_base,
            }
        )
    if seen != set(FOLDS):
        raise CorgiHeadContractError("cross-fit evaluation folds are incomplete")
    return sorted(normalized, key=lambda item: item["evaluation_fold"])


def validate_contract(value: Mapping[str, Any]) -> list[dict[str, Any]]:
    if (
        value.get("schema_version") != SCHEMA_VERSION
        or value.get("status") != "frozen_development_preflight_only"
        or value.get("model_id") != "corgi_regular"
        or value.get("task_id") != "rna_conditioned_atac"
        or value.get("dataset_id") != "gse296875"
        or value.get("lineage_id") != "hepatocyte"
        or value.get("adaptation_rung") != "new_assay_head"
        or value.get("biological_unit") != "donor"
        or value.get("genomic_block_unit") != "chromosome"
        or value.get("same_nucleus_topology_required") is not True
        or value.get("gse296875_masld_diagnosis_claimed") is not False
        or value.get("histology_or_disease_endpoint_used") is not False
        or value.get("test_or_sealed_features_used") is not False
        or value.get("test_or_sealed_outcomes_used") is not False
        or value.get("champion_or_external_transfer_claim_allowed") is not False
        or value.get("global_census_revision") is not False
        or value.get("historical_audit_hashes_modified") is not False
    ):
        raise CorgiHeadContractError("top-level scientific boundary differs")
    base = value.get("base_representation", {})
    if (
        base.get("score_name") != "strand_tta_orientation_mean_softplus_regional_sum"
        or base.get("primary_context_arm") != "actual_released_rank_masked"
        or set(base.get("context_controls", [])) != CONTEXT_CONTROLS
        or base.get("context_mapper_selected_in_this_contract") is not False
        or base.get("model_weights_updated") is not False
        or base.get("base_predictions_are_out_of_fold") is not True
    ):
        raise CorgiHeadContractError("base representation contract differs")
    head = value.get("profile_head", {})
    alpha = [float(item) for item in head.get("alpha_grid", [])]
    if (
        head.get("identity") != "one_parameter_geometric_profile_interpolation"
        or alpha != [0.0, 0.25, 0.5, 0.75, 1.0]
        or float(head.get("baseline_alpha", -1)) != 0.0
        or float(head.get("raw_corgi_alpha", -1)) != 1.0
        or float(head.get("epsilon_total_mass", -1)) != 1.0e-6
        or head.get("fit_loss") != "multinomial_deviance_per_insertion"
        or head.get("fit_aggregation") != "equal_weight_donor_by_chromosome_units"
        or head.get("selection_rule")
        != "one_standard_error_from_minimum_mean_fit_loss_then_smallest_alpha"
        or head.get("randomness") != "none"
        or head.get("fit_context_arm") != "actual_released_rank_masked_only"
        or head.get("apply_same_selected_alpha_to_context_controls") is not True
        or head.get("context_arm_specific_fitting_forbidden") is not True
        or head.get("held_fold_outcome_during_fit_forbidden") is not True
    ):
        raise CorgiHeadContractError("profile head contract differs")
    evaluation = value.get("evaluation", {})
    if (
        evaluation.get("primary_endpoint")
        != "conditional_multinomial_profile_deviance_skill_over_training_only_mean"
        or evaluation.get("secondary_endpoint") != "regional_count_spearman"
        or evaluation.get("macro_averaging_unit") != "donor_by_chromosome"
        or evaluation.get("promotion_gate_evaluated") is not False
        or evaluation.get("development_diagnostic_only") is not True
        or evaluation.get("full_family_ranking_allowed") is not False
    ):
        raise CorgiHeadContractError("evaluation boundary differs")
    terminal = value.get("terminal_boundaries", {})
    if any(
        terminal.get(field) is not False
        for field in (
            "native_numeric_parity_established",
            "open_champion_eligible",
            "external_validation_present",
            "full_donor_lineage_screen_complete",
        )
    ):
        raise CorgiHeadContractError("terminal boundary differs")
    return validate_fold_plan(value.get("crossfit_folds", []))


def freeze_contract(project_root: Path, contract_path: Path, output: Path) -> Mapping[str, Any]:
    if output.exists():
        raise CorgiHeadContractError("refusing to overwrite contract receipt")
    project_root = project_root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    value = tomllib.loads(contract_path.read_text(encoding="utf-8"))
    folds = validate_contract(value)
    sources = value.get("source_artifacts", {})
    bindings = (
        ("adaptation_contract_path", "adaptation_contract_sha256"),
        ("bounded_smoke_disposition_path", "bounded_smoke_disposition_sha256"),
        ("dataset_manifest_path", "dataset_manifest_sha256"),
        ("base_predictions_path", "base_predictions_artifacts_sha256"),
        ("donor_atac_path", "donor_atac_artifacts_sha256"),
        ("training_mean_atac_path", "training_mean_atac_artifacts_sha256"),
        (
            "design_exposed_smoke_evaluation_path",
            "design_exposed_smoke_evaluation_artifacts_sha256",
        ),
    )
    resolved = {
        path_key: _check_source(project_root, sources.get(path_key), sources.get(hash_key))
        for path_key, hash_key in bindings
    }
    prediction_metadata = _artifact_metadata(resolved["base_predictions_path"])
    donor_metadata = _artifact_metadata(resolved["donor_atac_path"])
    training_metadata = _artifact_metadata(resolved["training_mean_atac_path"])
    exposed_metadata = _artifact_metadata(resolved["design_exposed_smoke_evaluation_path"])
    if (
        prediction_metadata.get("artifact_class")
        != "Corgi_outcome_aligned_tile_smoke_bundle_v3"
        or prediction_metadata.get("outcomes_read") is not False
        or prediction_metadata.get("outcome_role") != "valid"
        or prediction_metadata.get("logical_tasks") != 5
        or donor_metadata.get("artifact_class") != "gse296875_deduplicated_tn5_bigwigs"
        or donor_metadata.get("biological_unit") != "donor"
        or training_metadata.get("artifact_class") != "gse296875_training_fold_pseudobulk"
        or training_metadata.get("biological_outer_unit") != "donor"
        or exposed_metadata.get("artifact_class")
        != "corgi_regional_development_smoke_evaluation"
        or exposed_metadata.get("held_test_or_sealed_outcomes_read") is not False
        or exposed_metadata.get("histology_or_disease_labels_read") is not False
        or exposed_metadata.get("promotion_gate_evaluated") is not False
        or exposed_metadata.get("champion_claim_allowed") is not False
    ):
        raise CorgiHeadContractError("source artifact firewall differs")
    for outer_fold in FOLDS:
        receipt = _load_json(resolved["base_predictions_path"] / f"fold{outer_fold}/receipt.json")
        valid_fold = (outer_fold + 1) % len(FOLDS)
        if (
            receipt.get("schema_version") != "masld-bench-corgi-outcome-aligned-prediction-v3"
            or receipt.get("status") != "pass_outcome_free_prediction"
            or receipt.get("mapper_outer_fold") != outer_fold
            or receipt.get("valid_donor_fold") != valid_fold
            or receipt.get("valid_genomic_fold") != valid_fold
            or receipt.get("outcome_role") != "valid"
            or receipt.get("held_ATAC_or_other_outcomes_used") is not False
            or receipt.get("test_or_sealed_features_or_labels_read") is not False
            or receipt.get("model_fitted_or_adapted") is not False
        ):
            raise CorgiHeadContractError("base prediction fold firewall differs")
    output.mkdir(parents=True, exist_ok=False)
    receipt = {
        "schema_version": RECEIPT_SCHEMA_VERSION,
        "status": "pass_outcome_free_crossfit_head_contract",
        "contract_id": value["contract_id"],
        "contract_sha256": digest(contract_path),
        "model_id": "corgi_regular",
        "task_id": "rna_conditioned_atac",
        "dataset_id": "gse296875",
        "lineage_id": "hepatocyte",
        "adaptation_rung": "new_assay_head",
        "folds": folds,
        "fit_units_per_head": 4,
        "evaluation_units_per_head": 1,
        "donor_and_genomic_fold_crossfit": True,
        "source_signal_arrays_read": False,
        "development_outcomes_read": False,
        "histology_or_disease_labels_read": False,
        "test_or_sealed_features_or_outcomes_read": False,
        "model_fit_performed": False,
        "evaluation_performed": False,
        "promotion_gate_evaluated": False,
        "champion_or_external_claim_allowed": False,
        "global_census_modified": False,
        "next_action": value["terminal_boundaries"]["next_action"],
    }
    (output / "contract_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    receipt = freeze_contract(args.project_root, args.contract, args.output)
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
