#!/usr/bin/env python3
"""Run the import-corrected v4 GSE267145 production bundle."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import verify_frozen_tree

import scripts.run_gse267145_histology_production_bundle_v3 as v3


REVISION_ID = "model-training-069-production-v4"
DELEGATED_ALGORITHM_REVISION_ID = "model-training-069-production-v3"
SOFTWARE_CONTRACT_PATH = (
    Path(__file__).resolve().parents[1]
    / "config/evaluation/gse267145_histology_production_v4.json"
)
SOFTWARE_CONTRACT_SHA256 = (
    "df1f77f82c3657e6786be9a67100b3ce138384f8b28eb6a22b48bc3e0cf97de2"
)
V4_FITTER_SHA256 = "1e0f9ec881a3d04601f0e2c4c375458ef7664c72cd5fb733e400ab08e51f3cff"
DIAGNOSTIC_PATH = (
    Path(__file__).resolve().parents[1] / "executions/model-check-075-21088989"
)
DIAGNOSTIC_ARTIFACTS_SHA256 = (
    "b068e3d9af7952ce6bfe49d7482b4453c9bbb15e6d49cbfd3c912dba9b4fd60d"
)


def _validate_diagnostic_receipt(receipt: Mapping[str, Any]) -> None:
    units = receipt.get("units")
    control = units.get("1/1701") if isinstance(units, Mapping) else None
    failure = units.get("1/1721") if isinstance(units, Mapping) else None
    if (
        receipt.get("schema_version")
        != "masld-bench-gse267145-production-v4-validation-v1"
        or receipt.get("status") != "passed_unscored"
        or receipt.get("revision_id") != REVISION_ID
        or receipt.get("delegated_algorithm_revision_id")
        != DELEGATED_ALGORITHM_REVISION_ID
        or receipt.get("actual_child_environment_self_test_passed") is not True
        or receipt.get("unit_tests_passed") is not True
        or receipt.get("v2_outer_1_seed_1701_predictions_byte_identical") is not True
        or receipt.get("primary_endpoint_code_path_changed") is not False
        or receipt.get("secondary_endpoint_code_paths_changed") is not False
        or receipt.get("outer_test_outcomes_read") is not False
        or receipt.get("metrics_calculated") is not False
        or receipt.get("scorer_called") is not False
        or receipt.get("production_submission_authorized") is not False
        or not isinstance(control, Mapping)
        or control.get("activation_count") != 0
        or control.get("activated_model_ids") != []
        or not isinstance(failure, Mapping)
        or failure.get("activation_count") != 1
        or failure.get("activated_model_ids") != ["h3_variance_pca_knn"]
    ):
        raise v3.v2.ProductionBundleError("v4 parity/failure validation differs")


def _verify_v4_diagnostic(path: Path) -> None:
    manifest = verify_frozen_tree(path)
    metadata = manifest.get("metadata", {})
    if (
        metadata.get("artifact_class")
        != "gse267145_histology_production_v4_validation"
        or metadata.get("actual_child_environment_self_test_passed") is not True
        or metadata.get("status") != "passed_unscored"
        or metadata.get("outcomes_read") is not False
        or metadata.get("metrics_calculated") is not False
        or metadata.get("scorer_called") is not False
        or metadata.get("production_submission_authorized") is not False
    ):
        raise v3.v2.ProductionBundleError("v4 validation metadata differs")
    _validate_diagnostic_receipt(
        v3.v2._read_json_object(path / "validation/receipt.json")
    )


def _source_hash_manifest(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise v3.v2.ProductionBundleError("v4 bundle source manifest is unreadable") from error
    for line in lines:
        digest, separator, source = line.partition("  ")
        if (
            separator != "  "
            or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
            or not source
            or source in records
        ):
            raise v3.v2.ProductionBundleError("v4 bundle source manifest differs")
        records[source] = digest
    return records


def _verify_bundle_validation(path: Path, expected_sha256: str) -> None:
    if (
        len(expected_sha256) != 64
        or any(character not in "0123456789abcdef" for character in expected_sha256)
    ):
        raise v3.v2.ProductionBundleError("v4 bundle validation SHA-256 is invalid")
    v3.v2._check_hash(
        path / "ARTIFACTS.json", expected_sha256, label="v4 bundle validation"
    )
    manifest = verify_frozen_tree(path)
    metadata = manifest.get("metadata", {})
    receipt = v3.v2._read_json_object(path / "validation/receipt.json")
    production = receipt.get("production_slurm", {})
    if (
        metadata.get("artifact_class")
        != "gse267145_histology_production_v4_bundle_validation"
        or metadata.get("reused_v2_units_byte_identical") != 22
        or metadata.get("v4_units_to_compute") != 3
        or metadata.get("status") != "passed_unscored"
        or metadata.get("outcomes_read") is not False
        or metadata.get("metrics_calculated") is not False
        or metadata.get("scorer_called") is not False
        or metadata.get("production_wrapper_submitted") is not False
        or receipt.get("schema_version")
        != "masld-bench-gse267145-production-v4-bundle-validation-v1"
        or receipt.get("status")
        != "passed_unscored_production_wrapper_not_submitted"
        or receipt.get("revision_id") != REVISION_ID
        or receipt.get("parity_validation_artifacts_sha256")
        != DIAGNOSTIC_ARTIFACTS_SHA256
        or receipt.get("reused_v2_units_byte_identical") != 22
        or receipt.get("v4_units_to_compute")
        != ["1/1721", "1/1723", "1/1733"]
        or receipt.get("failed_v3_attempts_reused") is not False
        or receipt.get("participant_outer_folds_artifacts_sha256")
        != v3.v2.FOLDS_ARTIFACTS_SHA256
        or production.get("job_name") != "model-training-077"
        or production.get("partition") != "cpu"
        or production.get("cpus") != 16
        or production.get("memory") != "32G"
        or production.get("wall_time") != "04:00:00"
        or production.get("fold_workers") != 5
        or production.get("blas_threads_per_active_fit") != 3
        or production.get("maximum_concurrent_new_unit_fits_after_reuse") != 1
        or receipt.get("outer_test_outcomes_read") is not False
        or receipt.get("metrics_calculated") is not False
        or receipt.get("scorer_called") is not False
        or receipt.get("production_wrapper_submitted") is not False
    ):
        raise v3.v2.ProductionBundleError("v4 bundle validation differs")

    records = _source_hash_manifest(path / "source.sha256")
    root = Path(__file__).resolve().parents[1]
    required = (
        "scripts/run_gse267145_histology_production_bundle.py",
        "scripts/run_gse267145_histology_production_bundle_v3.py",
        "scripts/run_gse267145_histology_production_bundle_v4.py",
        "scripts/fit_gse267145_histology_baselines.py",
        "scripts/fit_gse267145_histology_baselines_v3.py",
        "scripts/fit_gse267145_histology_baselines_v4.py",
        "scripts/aggregate_gse267145_histology_production_predictions.py",
        "config/evaluation/gse267145_histology_production_v3.json",
        "config/evaluation/gse267145_histology_production_v4.json",
        "config/evaluation/gse267145_histology_production_v3_reuse_manifest.json",
        "slurm/run_gse267145_histology_production_bundle_v4_cpu.sbatch",
    )
    for relative in required:
        if records.get(relative) != v3.v2.sha256_file(root / relative):
            raise v3.v2.ProductionBundleError(
                f"v4 bundle does not pin current production source: {relative}"
            )
    diagnostic_key = (DIAGNOSTIC_PATH / "ARTIFACTS.json").as_posix()
    if records.get(diagnostic_key) != DIAGNOSTIC_ARTIFACTS_SHA256:
        raise v3.v2.ProductionBundleError("v4 bundle does not pin parity validation")


def _bundle_validation_environment() -> tuple[Path, str]:
    root_value = os.environ.get("BUNDLE_VALIDATION_ROOT")
    expected = os.environ.get("BUNDLE_VALIDATION_SHA256")
    if not root_value or not expected:
        raise v3.v2.ProductionBundleError(
            "BUNDLE_VALIDATION_ROOT and BUNDLE_VALIDATION_SHA256 are required"
        )
    root = Path(root_value)
    if not root.is_absolute():
        raise v3.v2.ProductionBundleError("BUNDLE_VALIDATION_ROOT must be absolute")
    _verify_bundle_validation(root, expected)
    return root, expected


def _input_contract_v4(arguments: argparse.Namespace) -> None:
    v3._input_contract_v2(arguments)
    v3.v2._check_hash(
        SOFTWARE_CONTRACT_PATH,
        SOFTWARE_CONTRACT_SHA256,
        label="v4 software contract",
    )
    v3.v2._check_hash(
        v3.REUSE_MANIFEST_PATH,
        v3.REUSE_MANIFEST_SHA256,
        label="v3 reuse manifest retained by v4",
    )
    if arguments.campaign.name != REVISION_ID:
        raise v3.v2.ProductionBundleError("v4 campaign path has the wrong revision ID")
    if arguments.diagnostic != DIAGNOSTIC_PATH:
        raise v3.v2.ProductionBundleError("v4 diagnostic path differs")
    contract = v3.v2._read_json_object(SOFTWARE_CONTRACT_PATH)
    correction = contract.get("software_correction", {})
    if (
        contract.get("revision_id") != REVISION_ID
        or contract.get("parent_revision") != DELEGATED_ALGORITHM_REVISION_ID
        or correction.get("v3_fitter_sha256") != v3.V3_FITTER_SHA256
        or correction.get("modeling_code_changed") is not False
        or contract.get("scoring_authorized") is not False
    ):
        raise v3.v2.ProductionBundleError("v4 software contract differs")
    reused, _ = v3._reuse_manifest()
    for (outer_fold, seed), expected in sorted(reused.items()):
        v3._verify_reused_unit(
            v3._v2_unit(outer_fold, seed), outer_fold, seed, expected
        )
    _verify_v4_diagnostic(arguments.diagnostic)
    _bundle_validation_environment()


def _campaign_spec_v4(arguments: argparse.Namespace) -> dict[str, Any]:
    specification = v3._campaign_spec_v2(arguments)
    specification["schema_version"] = (
        "masld-bench-gse267145-production-campaign-spec-v4"
    )
    specification["revision_id"] = REVISION_ID
    specification["delegated_algorithm_revision_id"] = (
        DELEGATED_ALGORITHM_REVISION_ID
    )
    specification["software_contract"] = {
        "path": SOFTWARE_CONTRACT_PATH.as_posix(),
        "sha256": SOFTWARE_CONTRACT_SHA256,
        "modeling_code_changed": False,
    }
    specification["fallback_validation"] = {
        "artifacts_sha256": DIAGNOSTIC_ARTIFACTS_SHA256,
        "actual_child_environment_self_test_passed": True,
        "control_predictions_byte_identical_to_v2": True,
        "observed_failure_unit_fallback_activated": True,
    }
    reused, recomputed = v3._reuse_manifest()
    specification["unit_execution"] = {
        "reuse_manifest_path": v3.REUSE_MANIFEST_PATH.as_posix(),
        "reuse_manifest_sha256": v3.REUSE_MANIFEST_SHA256,
        "byte_identical_frozen_v2_units": len(reused),
        "v4_units_to_compute": [
            f"{outer}/{seed}" for outer, seed in sorted(recomputed)
        ],
        "failed_v3_attempts_reused": False,
    }
    validation_root, validation_sha256 = _bundle_validation_environment()
    specification["bundle_validation"] = {
        "path": validation_root.as_posix(),
        "artifacts_sha256": validation_sha256,
        "current_production_sources_pinned": True,
    }
    return specification


def main() -> int:
    v3.v2.FITTER_SHA256 = V4_FITTER_SHA256
    v3.v2.DIAGNOSTIC_ARTIFACTS_SHA256 = DIAGNOSTIC_ARTIFACTS_SHA256
    v3.v2._input_contract = _input_contract_v4
    v3.v2._campaign_spec = _campaign_spec_v4
    v3.v2._prepare_spec = v3._prepare_spec_v3
    v3.v2._verify_unit = v3._verify_unit_v3
    return v3.v2.main()


if __name__ == "__main__":
    raise SystemExit(main())
