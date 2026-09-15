#!/usr/bin/env python3
"""Freeze the production target-axis boundary for observed-multiome baselines."""

from __future__ import annotations

import argparse
import ast
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

from masld_bench.artifacts import freeze_tree, verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-observed-multiome-target-axis-contract-v1"
MODELS = {
    "masked_modality",
    "observed_atac_glm",
    "observed_atac_only",
    "rna_only",
    "shuffled_modality",
}


class TargetAxisContractError(ValueError):
    """Raised when the target-axis, crossed split, or outcome separation drifts."""


def _load(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise TargetAxisContractError(f"JSON object required: {path}")
    return value


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _verify_parent(root: Path, record: dict[str, Any], label: str) -> dict[str, Any]:
    path = (root / str(record.get("path", ""))).resolve(strict=True)
    path.relative_to(root)
    manifest = verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise TargetAxisContractError(f"{label} artifact drifted")
    return manifest


def validate(root: Path, config: dict[str, Any]) -> dict[str, Any]:
    if config.get("schema_version") != SCHEMA or config.get("contract_id") != "observed_multiome_target_axis_20260825":
        raise TargetAxisContractError("target-axis contract identity differs")
    parents = config.get("parent_artifacts")
    if not isinstance(parents, dict) or set(parents) != {"census", "synthetic_fixture"}:
        raise TargetAxisContractError("parent artifact roster differs")
    _verify_parent(root, parents["census"], "census")
    synthetic_manifest = _verify_parent(root, parents["synthetic_fixture"], "synthetic fixture")
    if (
        parents["synthetic_fixture"].get("disposition")
        != "mechanics_only_not_production_target_axis_compatible"
        or synthetic_manifest.get("metadata", {}).get("artifact_class")
        != "observed_multiome_baseline_fixture"
    ):
        raise TargetAxisContractError("synthetic fixture disposition differs")

    adapter = config.get("current_synthetic_adapter")
    if not isinstance(adapter, dict) or adapter.get("production_eligible") is not False:
        raise TargetAxisContractError("current adapter must remain production-ineligible")
    adapter_path = (root / str(adapter.get("path", ""))).resolve(strict=True)
    adapter_path.relative_to(root)
    if _digest(adapter_path) != adapter.get("sha256"):
        raise TargetAxisContractError("current synthetic adapter drifted")
    tree = ast.parse(adapter_path.read_text(encoding="utf-8"))
    functions = {node.name: node for node in tree.body if isinstance(node, ast.FunctionDef)}
    if not {"fit", "predict", "export"}.issubset(functions):
        raise TargetAxisContractError("synthetic adapter interface differs")
    prediction_args = {
        argument.arg
        for argument in (
            functions["predict"].args.args + functions["predict"].args.kwonlyargs
        )
    }
    if {"target", "target_counts"}.intersection(prediction_args):
        raise TargetAxisContractError("prediction interface accepts a target")

    axes = config.get("production_axes")
    expected_axes = {
        "biological_row_axis": "donor_by_lineage",
        "context_feature_axis": "outer_training_fitted_row_features",
        "genomic_target_axis": "prespecified_peak_or_profile_bin",
        "sequence_feature_axis": "genomic_target_by_sequence_feature",
        "target_matrix_axis": "donor_by_lineage_by_genomic_target",
        "cells_as_replicates": False,
    }
    if axes != expected_axes:
        raise TargetAxisContractError("production axes differ")

    split = config.get("crossed_split")
    if not isinstance(split, dict):
        raise TargetAxisContractError("crossed split is absent")
    if (
        split.get("donor_folds") != 5
        or split.get("genomic_folds") != 5
        or split.get("seeds") != [20260824, 20260825, 20260826, 20260827, 20260828]
        or split.get("fit_pair_rule") != "donor_not_in_held_fold_and_genomic_target_not_in_held_block"
        or split.get("prediction_pair_rule") != "donor_in_held_fold_and_genomic_target_in_held_block"
        or split.get("largest_receptive_field_buffer_bp") != 524288
    ):
        raise TargetAxisContractError("crossed split identity differs")
    false_fields = {
        "held_target_or_buffer_allowed_in_query_transform",
        "held_target_outcome_allowed_in_fit",
        "held_donor_allowed_in_fit",
        "target_specific_parameter_learned_from_held_block",
    }
    if any(split.get(field) is not False for field in false_fields):
        raise TargetAxisContractError("crossed split opens held evidence")

    baselines = config.get("baseline_production_contracts")
    if not isinstance(baselines, dict) or set(baselines) != MODELS:
        raise TargetAxisContractError("production baseline roster differs")
    for model_id in MODELS.difference({"shuffled_modality"}):
        if baselines[model_id].get("required_design") not in {
            "factorized_or_bilinear_target_generalizing",
            "factorized_target_generalizing",
        }:
            raise TargetAxisContractError(f"{model_id} lacks a target-generalizing design")
    shuffled = baselines["shuffled_modality"]
    if (
        shuffled.get("parent_bound") is not True
        or shuffled.get("target_outcome_used_for_permutation") is not False
        or shuffled.get("held_query_composition_used_for_permutation") is not False
    ):
        raise TargetAxisContractError("shuffled-modality firewall differs")

    prediction = config.get("prediction_contract")
    if (
        not isinstance(prediction, dict)
        or set(prediction.get("required_identifiers", []))
        != {"row_hash", "unit_hash", "block_hash", "model_id", "fold", "seed"}
        or prediction.get("target_argument_allowed_at_prediction") is not False
        or prediction.get("target_specific_intercept_from_held_block_allowed") is not False
        or prediction.get("partial_rectangle_ranking_allowed") is not False
        or prediction.get("benchmark_metric_calculated_by_adapter") is not False
    ):
        raise TargetAxisContractError("prediction contract differs")
    firewall = config.get("outcome_firewall")
    if not isinstance(firewall, dict) or set(firewall.values()) != {False}:
        raise TargetAxisContractError("outcome firewall is open")

    return {
        "schema_version": "masld-bench-observed-multiome-target-axis-audit-v1",
        "contract_id": config["contract_id"],
        "model_count": len(MODELS),
        "synthetic_fixture_mechanics_passed": True,
        "synthetic_fixture_production_eligible": False,
        "production_target_axis_frozen": True,
        "production_implementation_ready": False,
        "held_donor_or_block_evidence_opened": False,
        "biological_matrix_read": False,
        "development_outcome_read": False,
        "sealed_outcome_read": False,
        "metric_calculated": False,
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
    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    receipt = validate(root, _load(config_path))
    write_json_exclusive(output / "receipt.json", receipt)
    source_paths = [
        config_path,
        Path(__file__).resolve(strict=True),
        root / "tests/unit/test_observed_multiome_target_axis.py",
    ]
    (output / "source.sha256").write_text(
        "".join(f"{_digest(path)}  {path}\n" for path in source_paths),
        encoding="utf-8",
    )
    digest = freeze_tree(
        output,
        metadata={
            "artifact_class": "observed_multiome_target_axis_contract",
            "contract_id": receipt["contract_id"],
            "production_implementation_ready": False,
            "sealed_outcomes_accessed": False,
        },
    )
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
