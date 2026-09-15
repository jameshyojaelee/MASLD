#!/usr/bin/env python3
"""Run the v3 GSE267145 production bundle after frozen fallback validation."""

from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
from typing import Any, Mapping

from masld_bench.artifacts import verify_frozen_tree

import scripts.run_gse267145_histology_production_bundle as v2


REVISION_ID = "model-training-069-production-v3"
PARENT_REVISION = "model-training-069-production-v2"
FALLBACK_ID = "training_only_fibrosis_group3_class_prior_v1"
CONVERGENCE_FAILURE = (
    "HistologyBaselineError:elastic-net logistic candidate did not converge"
)
CONTRACT_PATH = (
    Path(__file__).resolve().parents[1]
    / "config/evaluation/gse267145_histology_production_v3.json"
)
REUSE_MANIFEST_PATH = (
    Path(__file__).resolve().parents[1]
    / "config/evaluation/gse267145_histology_production_v3_reuse_manifest.json"
)
V2_CAMPAIGN_PATH = (
    Path(__file__).resolve().parents[1]
    / "executions/model-training-069-production-v2"
)
CONTRACT_SHA256 = "8be11a742def0f80465bca92fd02baf626d2d207ee30a87cdb25bcf59e47e29f"
REUSE_MANIFEST_SHA256 = "8616379cdf1f892fd994b754758028c38bf3d74f5b4dc20657848d3b7a8aedef"
BASE_RUNNER_SHA256 = "5e2caef9396a040173842df59d9c01e83cd6ac35031fe7050045e6209ad3ca8e"
V3_FITTER_SHA256 = "0cf28366d86040ce18fec73e3a5cde862ab33dad71cf0719bdd1d15a14508d3e"
DIAGNOSTIC_ARTIFACTS_SHA256 = "e65a54a3afefd386ae4fa239de37665ad45bf7b8ce34f7ad8f0a15157653972a"


_input_contract_v2 = v2._input_contract
_campaign_spec_v2 = v2._campaign_spec
_prepare_spec_v2 = v2._prepare_spec
_verify_unit_v2 = v2._verify_unit


def _reuse_manifest() -> tuple[dict[tuple[int, int], str], set[tuple[int, int]]]:
    manifest = v2._read_json_object(REUSE_MANIFEST_PATH)
    reused_raw = manifest.get("reused_unit_artifacts_sha256")
    recomputed_raw = manifest.get("recomputed_units")
    if not isinstance(reused_raw, Mapping) or not isinstance(recomputed_raw, list):
        raise v2.ProductionBundleError("v3 reuse manifest roster is missing")
    try:
        reused = {
            tuple(int(value) for value in str(key).split("/")): str(digest)
            for key, digest in reused_raw.items()
        }
        recomputed = {
            tuple(int(value) for value in str(key).split("/"))
            for key in recomputed_raw
        }
    except ValueError as error:
        raise v2.ProductionBundleError("v3 reuse manifest unit key differs") from error
    complete = {(outer, seed) for outer in v2.OUTER_FOLDS for seed in v2.SEEDS}
    if (
        manifest.get("schema_version")
        != "masld-bench-gse267145-histology-production-v3-reuse-v1"
        or manifest.get("revision_id") != REVISION_ID
        or manifest.get("parent_revision") != PARENT_REVISION
        or manifest.get("parent_campaign")
        != "executions/model-training-069-production-v2"
        or len(reused) != 22
        or manifest.get("reused_unit_count") != 22
        or recomputed != {(1, 1721), (1, 1723), (1, 1733)}
        or manifest.get("recomputed_unit_count") != 3
        or set(reused).intersection(recomputed)
        or set(reused).union(recomputed) != complete
        or manifest.get("participant_outer_folds_artifacts_sha256")
        != v2.FOLDS_ARTIFACTS_SHA256
        or manifest.get("outer_test_outcomes_read") is not False
        or manifest.get("metrics_calculated") is not False
        or manifest.get("scorer_called") is not False
    ):
        raise v2.ProductionBundleError("v3 reuse manifest differs")
    if any(len(key) != 2 or len(digest) != 64 for key, digest in reused.items()):
        raise v2.ProductionBundleError("v3 reuse manifest hash roster differs")
    return reused, recomputed


def _v2_unit(outer_fold: int, seed: int) -> Path:
    return V2_CAMPAIGN_PATH / "units" / f"outer_{outer_fold}" / f"seed_{seed}"


def _verify_reused_unit(path: Path, outer_fold: int, seed: int, expected: str) -> None:
    _verify_unit_v2(path, outer_fold, seed)
    if v2.sha256_file(path / "ARTIFACTS.json") != expected:
        raise v2.ProductionBundleError("reused v2 unit ARTIFACTS SHA-256 differs")


def _populate_reused_units(arguments: argparse.Namespace, job_id: str) -> None:
    reused, _ = _reuse_manifest()
    for (outer_fold, seed), expected in sorted(reused.items()):
        source = _v2_unit(outer_fold, seed)
        _verify_reused_unit(source, outer_fold, seed, expected)
        target = (
            arguments.campaign / "units" / f"outer_{outer_fold}" / f"seed_{seed}"
        )
        if target.exists():
            _verify_reused_unit(target, outer_fold, seed, expected)
            continue
        target.parent.mkdir(parents=True, exist_ok=True)
        stage = target.with_name(f"{target.name}.reuse-staging-{job_id}")
        if stage.exists():
            raise v2.ProductionBundleError(f"refusing to reuse v3 copy stage: {stage}")
        shutil.copytree(source, stage, copy_function=shutil.copy2, symlinks=False)
        _verify_reused_unit(stage, outer_fold, seed, expected)
        v2.publish_directory_noreplace(stage, target)
        _verify_reused_unit(target, outer_fold, seed, expected)


def _validate_diagnostic_receipt(receipt: Mapping[str, Any]) -> None:
    units = receipt.get("units")
    control = units.get("1/1701") if isinstance(units, Mapping) else None
    failed_unit = units.get("1/1721") if isinstance(units, Mapping) else None
    if (
        receipt.get("schema_version")
        != "masld-bench-gse267145-production-v3-validation-v1"
        or receipt.get("status") != "passed_unscored"
        or receipt.get("revision_id") != REVISION_ID
        or receipt.get("unit_tests_passed") is not True
        or receipt.get("v2_outer_1_seed_1701_predictions_byte_identical") is not True
        or receipt.get("primary_endpoint_code_path_changed") is not False
        or receipt.get("other_secondary_endpoint_code_paths_changed") is not False
        or receipt.get("outer_test_outcomes_read") is not False
        or receipt.get("metrics_calculated") is not False
        or receipt.get("scorer_called") is not False
        or receipt.get("production_submission_authorized") is not False
        or not isinstance(control, Mapping)
        or control.get("activation_count") != 0
        or control.get("activated_model_ids") != []
        or not isinstance(failed_unit, Mapping)
        or not isinstance(failed_unit.get("activation_count"), int)
        or failed_unit["activation_count"] < 1
        or len(failed_unit.get("activated_model_ids", []))
        != failed_unit["activation_count"]
    ):
        raise v2.ProductionBundleError("v3 fallback validation receipt differs")


def _verify_v3_diagnostic(path: Path) -> None:
    manifest = verify_frozen_tree(path)
    metadata = manifest.get("metadata", {})
    if (
        metadata.get("artifact_class")
        != "gse267145_histology_production_v3_validation"
        or metadata.get("status") != "passed_unscored"
        or metadata.get("outcomes_read") is not False
        or metadata.get("metrics_calculated") is not False
        or metadata.get("scorer_called") is not False
        or metadata.get("production_submission_authorized") is not False
    ):
        raise v2.ProductionBundleError("v3 fallback validation metadata differs")
    _validate_diagnostic_receipt(
        v2._read_json_object(path / "validation/receipt.json")
    )


def _source_hash_manifest(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        raise v2.ProductionBundleError("bundle validation source manifest is unreadable") from error
    for line in lines:
        digest, separator, source = line.partition("  ")
        if (
            separator != "  "
            or len(digest) != 64
            or any(character not in "0123456789abcdef" for character in digest)
            or not source
            or source in records
        ):
            raise v2.ProductionBundleError("bundle validation source manifest differs")
        records[source] = digest
    return records


def _verify_bundle_validation(path: Path, expected_sha256: str) -> None:
    if (
        len(expected_sha256) != 64
        or any(character not in "0123456789abcdef" for character in expected_sha256)
    ):
        raise v2.ProductionBundleError("bundle validation SHA-256 is invalid")
    v2._check_hash(
        path / "ARTIFACTS.json",
        expected_sha256,
        label="v3 production bundle validation",
    )
    manifest = verify_frozen_tree(path)
    metadata = manifest.get("metadata", {})
    receipt = v2._read_json_object(path / "validation/receipt.json")
    production = receipt.get("production_slurm", {})
    if (
        metadata.get("artifact_class")
        != "gse267145_histology_production_v3_bundle_validation"
        or metadata.get("reused_v2_units_byte_identical") != 22
        or metadata.get("v3_units_to_compute") != 3
        or metadata.get("outcomes_read") is not False
        or metadata.get("metrics_calculated") is not False
        or metadata.get("scorer_called") is not False
        or metadata.get("production_wrapper_submitted") is not False
        or metadata.get("status") != "passed_unscored"
        or receipt.get("schema_version")
        != "masld-bench-gse267145-production-v3-bundle-validation-v1"
        or receipt.get("status")
        != "passed_unscored_production_wrapper_not_submitted"
        or receipt.get("revision_id") != REVISION_ID
        or receipt.get("fallback_validation_artifacts_sha256")
        != DIAGNOSTIC_ARTIFACTS_SHA256
        or receipt.get("reused_v2_units_byte_identical") != 22
        or receipt.get("v3_units_to_compute")
        != ["1/1721", "1/1723", "1/1733"]
        or receipt.get("failed_v2_attempt_reused") is not False
        or receipt.get("participant_outer_folds_artifacts_sha256")
        != v2.FOLDS_ARTIFACTS_SHA256
        or production.get("job_name") != "model-training-073"
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
        raise v2.ProductionBundleError("v3 production bundle validation differs")

    records = _source_hash_manifest(path / "source.sha256")
    root = Path(__file__).resolve().parents[1]
    required = (
        "scripts/run_gse267145_histology_production_bundle.py",
        "scripts/run_gse267145_histology_production_bundle_v3.py",
        "scripts/fit_gse267145_histology_baselines.py",
        "scripts/fit_gse267145_histology_baselines_v3.py",
        "scripts/aggregate_gse267145_histology_production_predictions.py",
        "config/evaluation/gse267145_histology_production_v3.json",
        "config/evaluation/gse267145_histology_production_v3_reuse_manifest.json",
        "slurm/run_gse267145_histology_production_bundle_v3_cpu.sbatch",
    )
    for relative in required:
        source = root / relative
        if records.get(relative) != v2.sha256_file(source):
            raise v2.ProductionBundleError(
                f"bundle validation does not pin current production source: {relative}"
            )
    diagnostic_key = (
        root / "executions/model-check-072-21088262/ARTIFACTS.json"
    ).as_posix()
    if records.get(diagnostic_key) != DIAGNOSTIC_ARTIFACTS_SHA256:
        raise v2.ProductionBundleError(
            "bundle validation does not pin the fallback validation"
        )


def _bundle_validation_environment() -> tuple[Path, str]:
    root_value = os.environ.get("BUNDLE_VALIDATION_ROOT")
    expected = os.environ.get("BUNDLE_VALIDATION_SHA256")
    if not root_value or not expected:
        raise v2.ProductionBundleError(
            "BUNDLE_VALIDATION_ROOT and BUNDLE_VALIDATION_SHA256 are required"
        )
    root = Path(root_value)
    if not root.is_absolute():
        raise v2.ProductionBundleError("BUNDLE_VALIDATION_ROOT must be absolute")
    _verify_bundle_validation(root, expected)
    return root, expected


def _input_contract_v3(arguments: argparse.Namespace) -> None:
    _input_contract_v2(arguments)
    v2._check_hash(CONTRACT_PATH, CONTRACT_SHA256, label="v3 fallback contract")
    v2._check_hash(
        REUSE_MANIFEST_PATH,
        REUSE_MANIFEST_SHA256,
        label="v3 reuse manifest",
    )
    v2._check_hash(
        Path(v2.__file__).resolve(), BASE_RUNNER_SHA256, label="v2 production runner"
    )
    if arguments.campaign.name != REVISION_ID:
        raise v2.ProductionBundleError("v3 campaign path has the wrong revision ID")
    contract = v2._read_json_object(CONTRACT_PATH)
    fallback = contract.get("fallback_contract", {})
    if (
        contract.get("revision_id") != REVISION_ID
        or contract.get("parent_revision") != PARENT_REVISION
        or fallback.get("fallback_id") != FALLBACK_ID
        or fallback.get("eligible_endpoint") != "fibrosis_group3"
        or contract.get("scoring_authorized") is not False
    ):
        raise v2.ProductionBundleError("v3 fallback contract differs")
    reused, _ = _reuse_manifest()
    for (outer_fold, seed), expected in sorted(reused.items()):
        _verify_reused_unit(_v2_unit(outer_fold, seed), outer_fold, seed, expected)
    _verify_v3_diagnostic(arguments.diagnostic)
    _bundle_validation_environment()


def _campaign_spec_v3(arguments: argparse.Namespace) -> dict[str, Any]:
    specification = _campaign_spec_v2(arguments)
    specification["schema_version"] = (
        "masld-bench-gse267145-production-campaign-spec-v3"
    )
    specification["revision_id"] = REVISION_ID
    specification["parent_revision"] = PARENT_REVISION
    specification["fallback_contract"] = {
        "contract_path": CONTRACT_PATH.as_posix(),
        "contract_sha256": CONTRACT_SHA256,
        "fallback_id": FALLBACK_ID,
        "eligible_endpoint": "fibrosis_group3",
        "eligible_endpoint_role": "secondary_only",
        "activation_requires_zero_valid_exact_grid_convergence_only_failures": True,
        "outer_test_outcomes_used": False,
        "metrics_calculated": False,
    }
    specification["validation"] = {
        "artifacts_sha256": DIAGNOSTIC_ARTIFACTS_SHA256,
        "control_predictions_byte_identical_to_v2": True,
        "observed_failure_unit_fallback_activated": True,
        "production_submission_authorized_by_validation_alone": False,
    }
    reused, recomputed = _reuse_manifest()
    specification["unit_execution"] = {
        "reuse_manifest_path": REUSE_MANIFEST_PATH.as_posix(),
        "reuse_manifest_sha256": REUSE_MANIFEST_SHA256,
        "byte_identical_frozen_v2_units": len(reused),
        "v3_units_to_compute": [
            f"{outer}/{seed}" for outer, seed in sorted(recomputed)
        ],
        "failed_v2_attempt_reused": False,
    }
    validation_root, validation_sha256 = _bundle_validation_environment()
    specification["bundle_validation"] = {
        "path": validation_root.as_posix(),
        "artifacts_sha256": validation_sha256,
        "current_production_sources_pinned": True,
    }
    return specification


def _prepare_spec_v3(arguments: argparse.Namespace, job_id: str) -> str:
    specification_sha256 = _prepare_spec_v2(arguments, job_id)
    _populate_reused_units(arguments, job_id)
    return specification_sha256


def _validate_v3_fit_receipt(
    receipt: Mapping[str, Any], selections: Mapping[str, Any]
) -> None:
    fallback = receipt.get("fibrosis_group3_fallback")
    if not isinstance(fallback, Mapping):
        raise v2.ProductionBundleError("v3 fit fallback receipt is missing")
    activated = fallback.get("activated_model_ids")
    if (
        receipt.get("schema_version")
        != "masld-bench-gse267145-histology-baseline-fit-preflight-v3"
        or receipt.get("revision_id") != REVISION_ID
        or receipt.get("status") != "passed_timing_predictions_unscored"
        or receipt.get("primary_endpoint_code_path_changed") is not False
        or receipt.get("other_secondary_endpoint_code_paths_changed") is not False
        or receipt.get("outer_test_outcomes_read") is not False
        or receipt.get("outer_test_metrics_calculated") is not False
        or fallback.get("fallback_id") != FALLBACK_ID
        or fallback.get("eligible_endpoint") != "fibrosis_group3"
        or fallback.get("eligible_endpoint_role") != "secondary_only"
        or fallback.get("outer_test_outcomes_used") is not False
        or fallback.get("metrics_calculated") is not False
        or not isinstance(activated, list)
        or activated != sorted(set(activated))
        or fallback.get("activation_count") != len(activated)
    ):
        raise v2.ProductionBundleError("v3 fit fallback receipt differs")

    observed: list[str] = []
    for model_id, model_receipt in selections.items():
        if not isinstance(model_receipt, Mapping):
            continue
        for endpoint in (
            "stage3",
            "nash_crn_component_sum",
            "fibrosis_exact_regression",
            "fibrosis_cumulative",
        ):
            selected = model_receipt.get(endpoint)
            if isinstance(selected, Mapping) and selected.get("fallback_triggered") is True:
                raise v2.ProductionBundleError("v3 fallback escaped fibrosis_group3")
        group = model_receipt.get("fibrosis_group3")
        if not isinstance(group, Mapping) or group.get("fallback_triggered") is not True:
            continue
        if (
            group.get("fallback_id") != FALLBACK_ID
            or group.get("parameters") != {"fallback_id": FALLBACK_ID}
            or group.get("valid_candidate_count") != 0
            or group.get("candidate_count") != 12
            or group.get("candidate_failure_count") != 12
            or group.get("candidate_failure_reasons")
            != {CONVERGENCE_FAILURE: 12}
            or group.get("one_standard_error_applied") is not False
            or group.get(
                "one_standard_error_not_applicable_no_fallback_hyperparameters"
            )
            is not True
            or group.get("outer_test_features_used") is not False
            or group.get("outer_test_outcomes_used") is not False
            or group.get("recorded_sex_used") is not False
            or group.get("source_stage5_used") is not False
            or group.get("randomness_used") is not False
        ):
            raise v2.ProductionBundleError("v3 activated fallback selection differs")
        observed.append(str(model_id))
    if sorted(observed) != activated:
        raise v2.ProductionBundleError("v3 fallback activation roster differs")


def _verify_unit_v3(unit: Path, outer_fold: int, seed: int) -> None:
    reused, recomputed = _reuse_manifest()
    key = (outer_fold, seed)
    if key in reused:
        _verify_reused_unit(unit, outer_fold, seed, reused[key])
        return
    if key not in recomputed:
        raise v2.ProductionBundleError("v3 unit is outside the frozen execution roster")
    _verify_unit_v2(unit, outer_fold, seed)
    receipt = v2._read_json_object(unit / "fit/receipt.json")
    selections = v2._read_json_object(
        unit
        / "fit/fit_receipts"
        / f"outer_{outer_fold}"
        / f"seed_{seed}"
        / "selection_receipts.json"
    )
    _validate_v3_fit_receipt(receipt, selections)


def main() -> int:
    v2.FITTER_SHA256 = V3_FITTER_SHA256
    v2.DIAGNOSTIC_ARTIFACTS_SHA256 = DIAGNOSTIC_ARTIFACTS_SHA256
    v2._input_contract = _input_contract_v3
    v2._campaign_spec = _campaign_spec_v3
    v2._prepare_spec = _prepare_spec_v3
    v2._verify_unit = _verify_unit_v3
    return v2.main()


if __name__ == "__main__":
    raise SystemExit(main())
