#!/usr/bin/env python3
"""Audit scBasset readiness without opening predictions, metrics, or outcomes."""

from __future__ import annotations

import argparse
from collections import Counter
import csv
from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


SCHEMA_VERSION = "masld-bench-scbasset-production-readiness-v1"
STATUS = "RECONCILED_DEVELOPMENT_ONLY_NO_GPU_UNIT"
FOLDS = tuple(range(5))
SCREENING_SEEDS = (20260824, 20260825, 20260826)
LEGACY_SEEDS = (11, 29, 47, 71, 101)
COMPLETE_CONTEXTS = ((1, 20260824), (2, 20260824), (3, 20260824), (4, 20260824))


class ScBassetReadinessError(ValueError):
    """Raised when an authority drifts or a closed check opens."""


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
        raise ScBassetReadinessError(f"JSON authority differs: {path}") from error
    if not isinstance(value, dict):
        raise ScBassetReadinessError(f"JSON authority is not an object: {path}")
    return value


def _load_toml(path: Path) -> dict[str, Any]:
    try:
        with path.open("rb") as handle:
            value = tomllib.load(handle)
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise ScBassetReadinessError(f"TOML authority differs: {path}") from error
    if not isinstance(value, dict):
        raise ScBassetReadinessError(f"TOML authority is not an object: {path}")
    return value


def _resolve_file(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ScBassetReadinessError("authority path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ScBassetReadinessError(f"authority is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise ScBassetReadinessError(
            f"authority escapes or is absent: {relative}"
        ) from error
    if not resolved.is_file():
        raise ScBassetReadinessError(f"authority is not a file: {relative}")
    return resolved


def _resolve_directory(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ScBassetReadinessError("artifact path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ScBassetReadinessError(f"artifact is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise ScBassetReadinessError(
            f"artifact escapes or is absent: {relative}"
        ) from error
    if not resolved.is_dir():
        raise ScBassetReadinessError(f"artifact is not a directory: {relative}")
    return resolved


def _require_file(root: Path, record: Mapping[str, Any]) -> Path:
    path = _resolve_file(root, str(record.get("path", "")))
    if _digest(path) != record.get("sha256"):
        raise ScBassetReadinessError(f"file authority differs: {path}")
    return path


def _require_frozen_artifact(
    root: Path,
    relative: str,
    expected_sha256: str,
    expected_class: str,
) -> tuple[Path, dict[str, Any]]:
    artifact_root = _resolve_directory(root, relative)
    try:
        manifest = verify_frozen_tree(artifact_root)
    except ArtifactError as error:
        raise ScBassetReadinessError(
            f"frozen artifact differs: {relative}"
        ) from error
    if (
        _digest(artifact_root / "ARTIFACTS.json") != expected_sha256
        or manifest.get("metadata", {}).get("artifact_class") != expected_class
    ):
        raise ScBassetReadinessError(f"frozen artifact identity differs: {relative}")
    return artifact_root, dict(manifest)


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _assert_false(mapping: Mapping[str, Any], fields: tuple[str, ...]) -> None:
    for field in fields:
        if mapping.get(field) is not False:
            raise ScBassetReadinessError(f"fail-closed field opened: {field}")


def _audit_semantics(root: Path, contract: Mapping[str, Any]) -> None:
    authorities = contract.get("source_authorities", {})
    checkpoint_path = _require_file(root, authorities.get("checkpoint_authority", {}))
    crosswalk_path = _require_file(root, authorities.get("development_crosswalk", {}))
    exposure_path = _require_file(root, authorities.get("exposure_audit", {}))
    registry_path = _require_file(root, authorities.get("model_registry", {}))
    rna_capability_path = _require_file(
        root, authorities.get("rna_conditioned_capabilities", {})
    )
    observed_task_path = _require_file(root, authorities.get("observed_multiome_task", {}))
    observed_capability_path = _require_file(
        root, authorities.get("observed_multiome_capabilities", {})
    )
    native_path = _require_file(root, authorities.get("native_adapter", {}))
    profile_path = _require_file(root, authorities.get("sequence_profile_adapter", {}))
    input_path = _require_file(root, authorities.get("input_builder", {}))
    training_path = _require_file(root, authorities.get("training_source", {}))
    prediction_path = _require_file(root, authorities.get("prediction_source", {}))

    checkpoint = _load_json(checkpoint_path)
    crosswalk = _load_json(crosswalk_path)
    exposure = _load_json(exposure_path)
    if (
        checkpoint.get("code_license") != "Apache-2.0"
        or checkpoint.get("architecture_contract", {}).get("input_length_bp") != 1344
        or checkpoint.get("architecture_contract", {}).get("output_tasks")
        != "one_sigmoid_output_per_training_cell"
        or checkpoint.get("inference_contract", {}).get("atac_available_at_inference")
        is not False
        or checkpoint.get("inference_contract", {}).get("held_cell_embedding")
        != "structurally_unavailable"
        or checkpoint.get("terms_audit", {}).get("result")
        != "local_training_allowed_but_public_tutorial_weights_not_admitted"
    ):
        raise ScBassetReadinessError("checkpoint semantics or terms differ")
    if (
        tuple(crosswalk.get("local_model_roster", {}).get("seeds", ()))
        != LEGACY_SEEDS
        or crosswalk.get("findings", {}).get("gse296875", {}).get("exposure_state")
        != "continual_seen"
        or crosswalk.get("findings", {}).get("gse289173", {}).get("exposure_state")
        != "target_label_unexposed"
        or exposure.get("checkpoint_findings", {}).get("scbasset_local", {}).get(
            "upstream_pretraining"
        )
        != "none"
        or exposure.get("checkpoint_findings", {}).get(
            "scbasset_tutorial_objects", {}
        ).get("exposure_state")
        != "reference_only"
    ):
        raise ScBassetReadinessError("crosswalk or exposure semantics differ")

    native = native_path.read_text(encoding="utf-8")
    profile = profile_path.read_text(encoding="utf-8")
    input_source = input_path.read_text(encoding="utf-8")
    training = training_path.read_text(encoding="utf-8")
    prediction = prediction_path.read_text(encoding="utf-8")
    if any(
        token not in native
        for token in (
            "SEQUENCE_LENGTH = 1344",
            "BOTTLENECK_DIMENSION = 32",
            'current, units=n_training_cells, activation="sigmoid"',
            ") != (None, n_training_cells):",
        )
    ):
        raise ScBassetReadinessError("native adapter semantics differ")
    if any(
        token not in profile
        for token in (
            "No held-donor ATAC or RNA is accepted at inference.",
            'raise ScBassetProfileError("scBasset has no inductive held-cell output columns")',
            "profile = predictions[:, indices].mean(axis=1, dtype=np.float64)",
            '"donor_context": "none_sequence_only"',
            '"held_donor_profiles_identical_within_lineage": True',
            '"observed_atac_exported": False',
            '"rna_input_exposed": False',
        )
    ):
        raise ScBassetReadinessError("sequence-profile semantics differ")
    if any(
        token not in training
        for token in (
            'summary.get("held_donor_atac_read") is not False',
            'summary.get("genomic_test_atac_used_for_matrix") is not False',
            '"training_initialization": "from_scratch"',
            '"public_tutorial_weights_used": False',
            '"held_donor_atac_used": False',
            '"genomic_test_atac_used": False',
        )
    ) or any(
        token not in prediction
        for token in (
            '"training_cell_output_columns_only": True',
            '"held_cell_embedding_available": False',
            '"held_donor_atac_used": False',
            '"rna_used": False',
        )
    ):
        raise ScBassetReadinessError("training or prediction firewall differs")
    if any(
        token not in input_source
        for token in (
            '"genomic_test_atac_used_for_matrix": False',
            '"model_training_input_contains_genomic_test_atac": False',
            '"missing_evidence_encoded_as_zero": False',
            '"held_donor_atac_exported": False',
        )
    ):
        raise ScBassetReadinessError("input-builder firewall differs")

    registry = _load_toml(registry_path)
    models = [row for row in registry.get("models", []) if row.get("model_id") == "scbasset"]
    if (
        len(models) != 1
        or models[0].get("license_status") != "Apache-2.0"
        or models[0].get("exposure_status") != "target_label_unexposed"
        or models[0].get("status") != "candidate"
        or models[0].get("admission_blocking") is not True
        or models[0].get("supported_tasks") != ["rna_conditioned_atac"]
    ):
        raise ScBassetReadinessError("global model registry differs")
    rna_capability = _load_toml(rna_capability_path)
    observed_task = _load_toml(observed_task_path)
    observed_capability = _load_toml(observed_capability_path)
    observed_models = set()
    for field in (
        "observed_atac_conditioned_track_candidates",
        "sequence_plus_rna_track_candidates",
        "joint_representation_candidates",
        "registered_baseline_models",
    ):
        observed_models.update(observed_capability.get(field, []))
    if (
        "scbasset" not in rna_capability.get("profile_baseline_models", [])
        or "scbasset" in rna_capability.get("profile_fixture_passed_models", [])
        or observed_task.get("task_id") != "observed_multiome"
        or observed_task.get("datasets_sealed") != []
        or "scbasset" in observed_models
        or observed_capability.get("execution_authorized") is not False
        or observed_capability.get("primary_eligible_models") != []
    ):
        raise ScBassetReadinessError("task capability boundary differs")

    semantics = contract.get("native_semantics", {})
    if (
        semantics.get("observed_atac_role") != "outer_training_supervision_only"
        or semantics.get("observed_query_atac_consumed") is not False
        or semantics.get("rna_context_consumed") is not False
        or semantics.get("held_cell_embedding_available") is not False
        or semantics.get("donor_context") != "none_sequence_only"
        or semantics.get("held_donor_profiles_identical_within_lineage") is not True
        or semantics.get("inductive_observed_multiome_model") is not False
        or semantics.get("allowed_observed_multiome_role")
        != "none_under_current_inductive_masked_query_task"
    ):
        raise ScBassetReadinessError("reconciled native semantics differ")


def _audit_smoke(root: Path, contract: Mapping[str, Any]) -> None:
    smoke = contract.get("smoke_and_runtime", {})
    _, runtime = _require_frozen_artifact(
        root,
        str(smoke.get("runtime_artifact_root", "")),
        str(smoke.get("runtime_artifacts_sha256", "")),
        "scbasset_native_runtime_probe",
    )
    _, training = _require_frozen_artifact(
        root,
        str(smoke.get("training_probe_artifact_root", "")),
        str(smoke.get("training_probe_artifacts_sha256", "")),
        "scbasset_training_script_probe",
    )
    _, prediction = _require_frozen_artifact(
        root,
        str(smoke.get("prediction_probe_artifact_root", "")),
        str(smoke.get("prediction_probe_artifacts_sha256", "")),
        "scbasset_prediction_script_probe",
    )
    runtime_meta = runtime.get("metadata", {})
    training_meta = training.get("metadata", {})
    prediction_meta = prediction.get("metadata", {})
    if (
        runtime_meta.get("device_class") != "NVIDIA_L40S_sm89"
        or runtime_meta.get("forward_backward_passed") is not True
        or runtime_meta.get("hdf5_weights_only_load_passed") is not True
        or runtime_meta.get("checkpoint_resume_step_passed") is not True
        or runtime_meta.get("public_tutorial_weights_used") is not False
        or training_meta.get("held_donor_atac_used") is not False
        or training_meta.get("public_tutorial_weights_used") is not False
        or prediction_meta.get("profile_fixture_passed") is not True
        or prediction_meta.get("reverse_complement_tta") is not True
        or prediction_meta.get("held_cell_embedding_available") is not False
        or prediction_meta.get("held_donor_atac_used") is not False
    ):
        raise ScBassetReadinessError("smoke or runtime evidence differs")


def _audit_coverage(root: Path, contract: Mapping[str, Any]) -> None:
    coverage = contract.get("five_fold_valid_coverage", {})
    audit_root, audit_manifest = _require_frozen_artifact(
        root,
        str(coverage.get("campaign_audit_artifact_root", "")),
        str(coverage.get("campaign_audit_artifacts_sha256", "")),
        "scbasset_all5_campaign_audit",
    )
    metadata = audit_manifest.get("metadata", {})
    if (
        metadata.get("outer_folds") != list(FOLDS)
        or metadata.get("donor_safe_artifact_campaign_complete") is not True
        or metadata.get("fixed_seed_three_seed_screen_complete") is not False
        or metadata.get("development_gate_passed") is not False
        or metadata.get("gpu_production_authorized") is not False
        or metadata.get("sealed_data_read") is not False
    ):
        raise ScBassetReadinessError("five-fold campaign metadata differs")
    receipt = _load_json(audit_root / "audit/campaign_audit.json")
    if (
        receipt.get("outer_folds_present") != list(FOLDS)
        or receipt.get("donors_partitioned_exactly_once_as_test") is not True
        or receipt.get("donors_partitioned_exactly_once_as_validation") is not True
        or receipt.get("uniform_prediction_lock_five_fold_campaign_complete")
        is not False
        or receipt.get("fixed_seed_three_seed_screen_complete") is not False
        or receipt.get("observed_models") != 5
        or receipt.get("prespecified_screen_models") != 15
        or receipt.get("seed_by_outer_fold")
        != {"0": 11, "1": 20260824, "2": 20260824, "3": 20260824, "4": 20260824}
        or receipt.get("raw_outcome_authority_opened") is not False
        or receipt.get("outcomes_used_for_fitting") is not False
        or receipt.get("champion_claim_allowed") is not False
        or receipt.get("development_gate_passed") is not False
        or receipt.get("gpu_production_authorized") is not False
    ):
        raise ScBassetReadinessError("five-fold structural receipt differs")
    if (
        coverage.get("outer_folds") != list(FOLDS)
        or coverage.get("donors") != 39
        or coverage.get("prediction_roles_present") != ["valid"]
        or coverage.get("test_role_predictions_present") is not False
        or coverage.get("mixed_seed_meta_aggregate") is not True
        or coverage.get("uniform_prediction_lock_folds") != [1, 2, 3, 4]
        or coverage.get("legacy_unlocked_fold") != 0
        or coverage.get("all_five_means_all_outer_folds_not_all_seeds_or_test_blocks")
        is not True
    ):
        raise ScBassetReadinessError("five-fold coverage interpretation differs")


def _audit_rectangle_and_dispatch(root: Path, contract: Mapping[str, Any]) -> None:
    rectangle = contract.get("fixed_seed_rectangle", {})
    _, matrix_gate = _require_frozen_artifact(
        root,
        str(rectangle.get("matrix_gate_artifact_root", "")),
        str(rectangle.get("matrix_gate_artifacts_sha256", "")),
        "scbasset_diagonal_matrix_gate",
    )
    plan_root, plan = _require_frozen_artifact(
        root,
        str(rectangle.get("historical_bundle_plan_artifact_root", "")),
        str(rectangle.get("historical_bundle_plan_artifacts_sha256", "")),
        "sequence_gpu_bundle_plan",
    )
    gate_meta = matrix_gate.get("metadata", {})
    if (
        gate_meta.get("outer_folds") != list(FOLDS)
        or gate_meta.get("models_per_fold") != 3
        or gate_meta.get("screening_seeds") != list(SCREENING_SEEDS)
        or gate_meta.get("prediction_roles") != ["valid"]
        or plan.get("metadata", {}).get("manual_gpu_sbatch_prohibited") is not True
        or plan.get("metadata", {}).get("future_qos") != "nslab"
    ):
        raise ScBassetReadinessError("rectangle gate or historical plan differs")

    tasks = [
        row
        for row in _read_tsv(plan_root / "plan/bundle_tasks.tsv")
        if row.get("model_id") == "scbasset"
    ]
    train = [row for row in tasks if row.get("stage") == "train"]
    predict = [row for row in tasks if row.get("stage") == "predict_valid"]
    expected = {(fold, seed) for fold in FOLDS for seed in SCREENING_SEEDS}
    if (
        len(tasks) != 30
        or len(train) != 15
        or len(predict) != 15
        or {(int(row["outer_fold"]), int(row["seed"])) for row in train} != expected
        or {(int(row["outer_fold"]), int(row["seed"])) for row in predict} != expected
        or any(row.get("prerequisite_logical_id") in {"", "none"} for row in predict)
        or any(
            row.get("submission_authority") != "dispatcher_only_no_manual_sbatch"
            for row in tasks
        )
    ):
        raise ScBassetReadinessError("fixed-seed task rectangle differs")

    complete = tuple(sorted(tuple(row) for row in rectangle.get("complete_contexts", [])))
    if (
        tuple(rectangle.get("outer_folds", ())) != FOLDS
        or tuple(rectangle.get("screening_seeds", ())) != SCREENING_SEEDS
        or rectangle.get("expected_fits") != 15
        or rectangle.get("complete_compatible_fits") != 4
        or rectangle.get("missing_fits") != 11
        or complete != COMPLETE_CONTEXTS
        or rectangle.get("legacy_fold0_seed11_is_rectangle_member") is not False
        or tuple(rectangle.get("legacy_crosswalk_seeds", ())) != LEGACY_SEEDS
        or rectangle.get("legacy_crosswalk_seeds_may_be_mixed_with_screening_rectangle")
        is not False
        or rectangle.get("fixed_seed_rectangle_complete") is not False
        or rectangle.get("partial_ranking_allowed") is not False
    ):
        raise ScBassetReadinessError("fixed-seed coverage differs")

    dispatcher = contract.get("dispatcher_reconciliation", {})
    queue_path = _resolve_file(root, str(dispatcher.get("historical_queue_item_path", "")))
    claim_path = _resolve_file(root, str(dispatcher.get("historical_claim_path", "")))
    if (
        _digest(queue_path) != dispatcher.get("historical_queue_item_sha256")
        or _digest(claim_path) != dispatcher.get("historical_claim_sha256")
    ):
        raise ScBassetReadinessError("historical dispatch identity differs")
    queue = _load_json(queue_path)
    claim = _load_json(claim_path)
    candidate = str(dispatcher.get("historical_plan_candidate_bundle_id", ""))
    allocations = {
        row["bundle_id"]: row
        for row in _read_tsv(plan_root / "plan/bundle_allocations.tsv")
    }
    candidate_row = allocations.get(candidate, {})
    if (
        queue.get("bundle_id") != "scbasset-fold1to4-chain-v2-generic1"
        or claim.get("bundle_id") != queue.get("bundle_id")
        or claim.get("status") != "submitted"
        or str(claim.get("job_id")) != "21069969"
        or candidate_row.get("models") != "scbasset"
        or candidate_row.get("logical_tasks") != "6"
        or candidate_row.get("initial_readiness") != "ready_with_internal_dependencies"
        or candidate_row.get("submission_authority")
        != "dispatcher_only_no_manual_sbatch"
    ):
        raise ScBassetReadinessError("historical or candidate dispatch semantics differ")
    registered_candidate_paths = []
    for path in sorted((root / "config/campaigns/gpu_bundle_queue").glob("*.json")):
        value = _load_json(path)
        if value.get("bundle_id") == candidate or value.get("exports", {}).get("BUNDLE_ID") == candidate:
            registered_candidate_paths.append(path.name)
    if registered_candidate_paths:
        raise ScBassetReadinessError("blocked historical plan candidate is registered")
    if (
        dispatcher.get("historical_item_already_claimed") is not True
        or dispatcher.get("historical_item_is_next_gpu_unit") is not False
        or dispatcher.get("historical_plan_predates_later_terminal_development_audit")
        is not True
        or dispatcher.get("later_development_gate_supersedes_gpu_candidate") is not True
        or dispatcher.get("valid_next_centrally_dispatched_gpu_unit_exists") is not False
        or dispatcher.get("new_queue_registration_authorized") is not False
        or dispatcher.get("direct_gpu_submission_authorized") is not False
    ):
        raise ScBassetReadinessError("dispatch gate differs")


def audit_readiness(root: Path, contract_path: Path) -> dict[str, Any]:
    """Verify only source, manifest, and structural receipt authorities."""

    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    try:
        contract_path.relative_to(root)
    except ValueError as error:
        raise ScBassetReadinessError("contract escapes benchmark root") from error
    contract = _load_json(contract_path)
    if contract.get("schema_version") != SCHEMA_VERSION or contract.get("status") != STATUS:
        raise ScBassetReadinessError("readiness schema or status differs")
    _audit_semantics(root, contract)
    _audit_smoke(root, contract)
    _audit_coverage(root, contract)
    _audit_rectangle_and_dispatch(root, contract)

    promotion = contract.get("license_exposure_and_promotion", {})
    if (
        promotion.get("code_license") != "Apache-2.0"
        or promotion.get("project_weights") != "locally_trained_from_scratch"
        or promotion.get("project_weight_redistribution_decision") != "not_yet_frozen"
        or promotion.get("tutorial_weight_terms") != "undeclared"
        or promotion.get("tutorial_weights_admitted") is not False
        or promotion.get("rna_conditioned_profile_fixture_artifact_exists") is not True
        or promotion.get("rna_conditioned_global_capability_registry_lists_fixture_passed")
        is not False
        or promotion.get("observed_multiome_capability_registry_lists_scbasset")
        is not False
        or promotion.get("observed_multiome_execution_authorized") is not False
        or promotion.get("target_label_unexposed_equals_champion_eligible") is not False
        or promotion.get("open_champion_eligible_now") is not False
        or promotion.get("universal_claim_allowed") is not False
    ):
        raise ScBassetReadinessError("license, exposure, or promotion gate differs")

    firewall = contract.get("firewall", {})
    _assert_false(
        firewall,
        (
            "prediction_values_may_be_read",
            "metric_values_may_be_read",
            "raw_outcomes_may_be_read",
            "sealed_features_may_be_read",
            "sealed_labels_or_outcomes_may_be_read",
            "partial_results_may_be_ranked",
            "thresholds_may_be_changed",
            "gpu_queue_may_be_registered_by_this_audit",
            "gpu_job_may_be_submitted_by_this_audit",
        ),
    )

    return {
        "schema_version": "masld-bench-scbasset-production-readiness-receipt-v1",
        "status": "pass_reconciled_development_only_no_gpu_unit",
        "contract_sha256": _digest(contract_path),
        "native_observed_atac_role": "outer_training_supervision_only",
        "observed_query_atac_consumed": False,
        "inductive_observed_multiome_model": False,
        "rna_conditioned_model": False,
        "donor_context": "none_sequence_only",
        "held_cell_embedding_available": False,
        "five_outer_fold_valid_campaign_present": True,
        "test_role_predictions_present": False,
        "mixed_seed_meta_aggregate": True,
        "uniform_five_fold_prediction_lock_complete": False,
        "fixed_seed_rectangle_expected_fits": 15,
        "fixed_seed_rectangle_complete_fits": 4,
        "fixed_seed_rectangle_missing_fits": 11,
        "fixed_seed_rectangle_complete": False,
        "partial_ranking_authorized": False,
        "development_gate_passed": False,
        "profile_fixture_artifact_exists": True,
        "global_profile_fixture_admission_complete": False,
        "observed_multiome_task_eligible": False,
        "tutorial_weights_admitted": False,
        "open_champion_eligible_now": False,
        "historical_gpu_item_already_claimed": True,
        "blocked_plan_candidate_registered": False,
        "valid_next_centrally_dispatched_gpu_unit_exists": False,
        "new_gpu_queue_registration_authorized": False,
        "direct_gpu_submission_authorized": False,
        "prediction_values_read": False,
        "metric_values_read": False,
        "raw_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_labels_or_outcomes_read": False,
        "partial_results_ranked": False,
        "thresholds_changed": False,
        "gpu_queue_registered": False,
        "gpu_job_submitted": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    receipt = audit_readiness(args.project_root, args.contract)
    rendered = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
