#!/usr/bin/env python3
"""Validate the outcome-blind observed-multiome specialist check."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

import numpy as np

from masld_bench.artifacts import freeze_tree, reject_symlink_components, verify_frozen_tree, write_json_exclusive
from masld_bench.observed_multiome_specialist_gate import (
    BASELINE_MODEL_ID,
    BOOTSTRAP_REPLICATES,
    BOOTSTRAP_SEED,
    COMPARISON_OUTPUT_FAMILY,
    GENOMIC_BLOCKS,
    LINEAGES,
    N_DONORS,
    N_UNIT_ROWS,
    SEEDS,
    UnitAxis,
    validate_unit_axis,
)


SCHEMA = "masld-bench-observed-multiome-specialist-gate-readiness-v1"


class SpecialistGateValidationError(ValueError):
    """Raised when the frozen check requirements or baseline authority differs."""


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            value.update(block)
    return value.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SpecialistGateValidationError("JSON object required")
    return value


def _bound_tree(root: Path, record: dict[str, Any], label: str) -> Path:
    path = reject_symlink_components(root / str(record.get("path", "")), label=label).resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise SpecialistGateValidationError(f"{label} drifted")
    return path


def validate(root: Path, config: dict[str, Any]) -> dict[str, Any]:
    if config.get("schema_version") != SCHEMA or config.get("dataset_id") != "gse296875" or config.get("stage") != "development":
        raise SpecialistGateValidationError("gate identity differs")
    sources = config.get("source_bindings")
    if not isinstance(sources, dict) or set(sources) != {"module", "validator", "unit_test", "sbatch"}:
        raise SpecialistGateValidationError("source binding roster differs")
    for label, record in sources.items():
        path = reject_symlink_components(root / str(record.get("path", "")), label=label).resolve(strict=True)
        path.relative_to(root)
        if _digest(path) != record.get("sha256"):
            raise SpecialistGateValidationError(f"source drifted: {label}")

    production = _bound_tree(root, config.get("baseline_production", {}), "baseline production")
    verification = _bound_tree(root, config.get("independent_verification", {}), "independent verification")
    verified = _json(verification / "receipt.json")
    expected_verification = {
        "independent_rederivation_passed": True,
        "strongest_development_control_verified": BASELINE_MODEL_ID,
        "production_artifacts_sha256": config["baseline_production"]["artifacts_sha256"],
        "biological_donors_verified": N_DONORS,
        "lineages_verified": list(LINEAGES),
        "genomic_blocks_verified": len(GENOMIC_BLOCKS),
        "surface_seed_runs_verified": 125,
        "partial_rectangle_used": False,
        "supports_external_claim": False,
        "supports_champion_claim": False,
        "sealed_outcomes_read": False,
    }
    if any(verified.get(key) != value for key, value in expected_verification.items()):
        raise SpecialistGateValidationError("independent baseline verification differs")

    with np.load(production / "unit_results.npz", allow_pickle=False) as archive:
        if set(archive.files) != {
            "unit_hash",
            "lineage",
            "genomic_fold",
            "model_id",
            "ensemble_profile_deviance_skill",
            "seed_profile_deviance_skill",
        }:
            raise SpecialistGateValidationError("baseline unit artifact field roster differs")
        models = tuple(archive["model_id"].astype(str))
        axis = UnitAxis(
            tuple(archive["unit_hash"].astype(str)),
            tuple(archive["lineage"].astype(str)),
            tuple(int(value) for value in archive["genomic_fold"]),
        )
    validate_unit_axis(axis)
    if BASELINE_MODEL_ID not in models:
        raise SpecialistGateValidationError("verified baseline model is absent")

    candidates = config.get("candidate_registry")
    if not isinstance(candidates, list) or {row.get("model_id") for row in candidates} != {"epibert", "epcotv2", "scooby", "multivi", "peakvi"}:
        raise SpecialistGateValidationError("specialist candidate roster differs")
    expected_native = {
        "epibert": "base_resolution_accessibility_track",
        "epcotv2": "base_resolution_accessibility_track",
        "scooby": "base_resolution_accessibility_track",
        "multivi": "peak_accessibility_rate",
        "peakvi": "peak_accessibility_rate",
    }
    for record in candidates:
        model_id = record.get("model_id")
        if record != {
            "model_id": model_id,
            "native_output_family": expected_native.get(model_id),
            "native_lane_preserved": True,
            "common_profile_entry": "requires_prediction_level_projection_frozen_before_outcome_scoring",
            "direct_cross_family_comparison_allowed": False,
        }:
            raise SpecialistGateValidationError(f"candidate family contract differs: {model_id}")

    rules = config.get("comparison_rules")
    expected_rules = {
        "baseline_model_id": BASELINE_MODEL_ID,
        "comparison_output_family": COMPARISON_OUTPUT_FAMILY,
        "fixed_seeds": list(SEEDS),
        "seed_ensemble": "arithmetic_mean_prediction_before_scoring",
        "biological_donors": N_DONORS,
        "lineages": list(LINEAGES),
        "genomic_blocks": list(GENOMIC_BLOCKS),
        "unit_rows": N_UNIT_ROWS,
        "bootstrap": "paired_two_way_donor_by_genomic_block_with_lineage_macro_within_cell",
        "bootstrap_replicates": BOOTSTRAP_REPLICATES,
        "bootstrap_seed": BOOTSTRAP_SEED,
        "minimum_relative_residual_deviance_gain": 0.05,
        "minimum_improved_lineages": 4,
        "maximum_lineage_degradation": 0.02,
        "minimum_positive_gain_seeds": 4,
        "paired_bootstrap_lower_bound_must_exceed_zero": True,
        "partial_rectangle_allowed": False,
        "external_claim_allowed": False,
        "champion_claim_allowed": False,
    }
    if rules != expected_rules:
        raise SpecialistGateValidationError("comparison rules differ")
    return {
        "schema_version": "masld-bench-observed-multiome-specialist-gate-readiness-receipt-v1",
        "dataset_id": "gse296875",
        "stage": "development",
        "baseline_production_artifacts_sha256": config["baseline_production"]["artifacts_sha256"],
        "independent_verification_artifacts_sha256": config["independent_verification"]["artifacts_sha256"],
        "baseline_model_id": BASELINE_MODEL_ID,
        "baseline_unit_axis_verified": True,
        "candidate_models_registered": [record["model_id"] for record in candidates],
        "native_output_families_preserved": True,
        "direct_cross_family_comparison_allowed": False,
        "common_profile_projection_required_before_scoring": True,
        "candidate_value_fields_loaded": [],
        "baseline_value_fields_loaded": [],
        "new_outcomes_read": False,
        "sealed_outcomes_read": False,
        "gate_ready": True,
        "supports_external_claim": False,
        "supports_champion_claim": False,
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
    output = reject_symlink_components(args.output, label="specialist gate validation output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    receipt = validate(root, _json(config_path))
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True), root / "src/masld_bench/observed_multiome_specialist_gate.py"]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_specialist_gate_readiness", "outcome_blind": True, "sealed_outcomes_accessed": False, "supports_external_claim": False, "supports_champion_claim": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
