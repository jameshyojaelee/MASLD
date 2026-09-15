#!/usr/bin/env python3
"""Freeze the additive supervised-training separation for observed multiome."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import (
    freeze_tree,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)


SCHEMA = "masld-bench-observed-multiome-supervised-training-contract-v1"
CONTRACT_ID = "gse296875_observed_multiome_supervised_training_20260825"


class SupervisedTrainingContractError(ValueError):
    """Raised when held donor or genomic outcomes can enter supervised fit."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SupervisedTrainingContractError("JSON object required")
    return value


def _tree(root: Path, record: Mapping[str, Any], label: str) -> Path:
    lexical = reject_symlink_components(root / str(record.get("path", "")), label=label)
    path = lexical.resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise SupervisedTrainingContractError(f"{label} drifted")
    return path


def validate(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    if config.get("schema_version") != SCHEMA or config.get("contract_id") != CONTRACT_ID or config.get("dataset_id") != "gse296875" or config.get("stage") != "smoke":
        raise SupervisedTrainingContractError("contract identity differs")
    parents = config.get("parents")
    expected_parents = {"target_axis_contract", "factorized_fixture", "materialization_contract_v1", "verified_model_inputs", "verified_evaluator_outcomes"}
    if not isinstance(parents, dict) or set(parents) != expected_parents:
        raise SupervisedTrainingContractError("parent roster differs")
    resolved = {key: _tree(root, record, key) for key, record in parents.items()}
    model_receipt = _json(resolved["verified_model_inputs"] / "receipt.json")
    evaluator_receipt = _json(resolved["verified_evaluator_outcomes"] / "receipt.json")
    if model_receipt.get("promotion_gate_passed") is not True or evaluator_receipt.get("promotion_gate_passed") is not True:
        raise SupervisedTrainingContractError("materialized sibling is not verified")
    gap = config.get("resolved_gap")
    if not isinstance(gap, dict) or gap != {"materialization_contract_v1_model_and_evaluator_separation_remains_binding": True, "materialization_contract_v1_fit_binding_is_insufficient_for_supervised_models": True, "reason": "model_inputs_contain_no_held_target_values_and_fit_may_not_bind_evaluator_outcomes", "resolution": "add_distinct_outer_training_label_artifacts_that_never_contain_held_donor_or_held_genomic_outcomes", "v1_model_and_evaluator_artifacts_rebuilt": False, "supervised_fit_authorized_now": False}:
        raise SupervisedTrainingContractError("resolved gap differs")
    surfaces = config.get("crossed_outer_surfaces")
    if surfaces != {"donor_folds": 5, "genomic_folds": 5, "surfaces": 25, "biological_unit": "donor_by_lineage", "cells_or_nuclei_as_replicates": False}:
        raise SupervisedTrainingContractError("outer surface contract differs")
    target = config.get("training_target_plan")
    if not isinstance(target, dict) or target != {"targets_per_nonheld_source_genomic_fold": 250, "training_targets_per_held_genomic_fold": 1000, "selection_uses_count_values": False, "selection_order": "sha256_of_training_target_namespace_held_fold_and_peak_id", "must_be_outside_held_genomic_fold": True, "must_be_disjoint_from_observed_atac_input_peaks": True, "must_be_disjoint_from_evaluator_target_peaks": True, "minimum_distance_from_observed_atac_input_peak_on_same_chromosome_bp": 524288, "insufficient_candidates": "fail_closed_without_reducing_buffer_or_target_count"}:
        raise SupervisedTrainingContractError("training target plan differs")
    labels = config.get("training_label_artifacts")
    if (
        not isinstance(labels, dict)
        or labels.get("children") != 25
        or labels.get("one_separately_frozen_child_per_held_genomic_by_held_donor_fold") is not True
        or labels.get("held_donor_rows") != "forbidden"
        or labels.get("held_genomic_target_values") != "forbidden"
        or labels.get("rna_values") != "forbidden"
        or labels.get("observed_atac_input_values") != "forbidden"
        or labels.get("evaluator_artifact_values") != "forbidden"
        or labels.get("raw_donor_ids") != "forbidden"
    ):
        raise SupervisedTrainingContractError("training label artifact differs")
    binding = config.get("run_binding_v2")
    if (
        not isinstance(binding, dict)
        or "one_matching_verified_training_label_artifact_for_held_genomic_and_held_donor_fold" not in binding.get("fit_may_bind", [])
        or "any_evaluator_outcome_artifact_or_verifier" not in binding.get("fit_may_not_bind", [])
        or "held_donor_training_target_values" not in binding.get("fit_may_not_bind", [])
        or "training_label_artifact" not in binding.get("predict_may_not_bind", [])
        or "prediction_bundle_committed_before_join" not in binding.get("evaluator_may_bind", [])
    ):
        raise SupervisedTrainingContractError("run binding differs")
    if config.get("promotion_sequence") != ["identifier_only_training_target_mask_plan", "training_label_materialization", "independent_training_label_verification", "prediction_only_synthetic_broker_fixture", "single_surface_single_seed_biological_smoke", "twenty_five_surface_five_seed_development_screen"]:
        raise SupervisedTrainingContractError("promotion sequence differs")
    firewall = config.get("firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise SupervisedTrainingContractError("contract firewall is open")
    return {
        "schema_version": "masld-bench-observed-multiome-supervised-training-contract-receipt-v1",
        "contract_id": CONTRACT_ID,
        "crossed_outer_surfaces": 25,
        "training_label_artifacts_planned": 25,
        "training_targets_per_held_genomic_fold": 1000,
        "held_donor_rows_in_training_labels": False,
        "held_genomic_values_in_training_labels": False,
        "model_fit_may_bind_evaluator_outcomes": False,
        "model_fit_may_bind_verified_fold_local_training_labels": True,
        "model_and_evaluator_artifacts_reused_without_rebuild": True,
        "biological_count_values_read": False,
        "development_metric_calculated": False,
        "model_fit": False,
        "sealed_outcomes_read": False,
        "supervised_fit_authorized": False,
        "next_gate": "identifier_only_training_target_mask_plan"
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    output = reject_symlink_components(args.output, label="contract output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    receipt = validate(root, _json(config_path))
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "tests/unit/test_observed_multiome_supervised_training_contract.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_supervised_training_contract", "contract_id": CONTRACT_ID, "supervised_fit_authorized": False, "sealed_outcomes_accessed": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
