#!/usr/bin/env python3
"""Audit the frozen regular-Corgi bounded prediction-smoke disposition."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any


SCHEMA_VERSION = "masld-bench-corgi-regular-bounded-prediction-smoke-v1"
COMPLETE_STATUS = "BOUNDED_OUTCOME_FREE_PREDICTION_SMOKE_COMPLETE"
ADMISSION_CLASS = "restricted_development_comparator"
CHECKPOINT_SHA256 = (
    "cd51539f01de1e66a3a4f9b89e772bdeae07df2f0178d5b692f2368db2f4d2a4"
)
CONTEXT_ARMS = [
    "actual_released_rank_masked",
    "actual_length_adjusted_tpm_rank_masked",
    "training_lineage_mean_released_rank",
    "nearest_training_released_rank",
    "shuffled_valid_released_rank",
]
SCORE_NAMES = [
    "forward_raw_regional_mean",
    "reverse_complement_raw_regional_mean",
    "strand_tta_raw_regional_mean",
    "strand_tta_orientation_mean_softplus_regional_sum",
]


class CorgiRegularSmokeError(ValueError):
    """Raised when the bounded regular-Corgi smoke requirements drifts."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _resolve_member(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise CorgiRegularSmokeError("authority member path is not relative")
    member = (root / relative).resolve(strict=True)
    try:
        member.relative_to(root)
    except ValueError as error:
        raise CorgiRegularSmokeError("authority member escapes benchmark root") from error
    if not member.is_file() or member.is_symlink():
        raise CorgiRegularSmokeError(f"authority member is not a regular file: {relative}")
    return member


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise CorgiRegularSmokeError(f"JSON authority is not an object: {path}")
    return value


def _artifact_index(manifest: dict[str, Any]) -> dict[str, dict[str, Any]]:
    artifacts = manifest.get("artifacts")
    if not isinstance(artifacts, list):
        raise CorgiRegularSmokeError("artifact manifest inventory differs")
    indexed: dict[str, dict[str, Any]] = {}
    for record in artifacts:
        if not isinstance(record, dict) or not isinstance(record.get("path"), str):
            raise CorgiRegularSmokeError("artifact manifest record differs")
        path = record["path"]
        if path in indexed:
            raise CorgiRegularSmokeError(f"duplicate artifact manifest member: {path}")
        indexed[path] = record
    return indexed


def _require_metadata(
    manifest: dict[str, Any], expected: dict[str, Any], label: str
) -> None:
    if manifest.get("schema_version") != "masld-bench-artifacts-v1":
        raise CorgiRegularSmokeError(f"{label} artifact schema differs")
    metadata = manifest.get("metadata")
    if not isinstance(metadata, dict):
        raise CorgiRegularSmokeError(f"{label} artifact metadata differs")
    for key, value in expected.items():
        if metadata.get(key) != value:
            raise CorgiRegularSmokeError(f"{label} metadata differs at {key}")


def audit_smoke_plan(root: Path, authority_path: Path) -> dict[str, Any]:
    """Verify only frozen authorities, ARTIFACTS metadata, and JSON receipts.

    This function intentionally does not verify frozen trees because doing so would
    read prediction arrays. It does not traverse evaluator directories.
    """

    root = root.resolve(strict=True)
    authority_path = authority_path.resolve(strict=True)
    try:
        authority_path.relative_to(root)
    except ValueError as error:
        raise CorgiRegularSmokeError("disposition authority escapes benchmark root") from error
    if authority_path.is_symlink() or not authority_path.is_file():
        raise CorgiRegularSmokeError("disposition authority is not a regular file")
    authority = _load_json(authority_path)
    if authority.get("schema_version") != SCHEMA_VERSION:
        raise CorgiRegularSmokeError("regular-Corgi smoke schema differs")
    if authority.get("status") != COMPLETE_STATUS:
        raise CorgiRegularSmokeError("regular-Corgi smoke completion status differs")
    if authority.get("admission_class") != ADMISSION_CLASS:
        raise CorgiRegularSmokeError("regular-Corgi admission class differs")

    authority_hashes = authority.get("authority_hashes")
    implementation_hashes = authority.get("implementation_hashes")
    manifest_hashes = authority.get("artifact_manifest_hashes")
    if not isinstance(authority_hashes, dict) or len(authority_hashes) != 6:
        raise CorgiRegularSmokeError("regular-Corgi authority hash family differs")
    if not isinstance(implementation_hashes, dict) or len(implementation_hashes) != 3:
        raise CorgiRegularSmokeError("regular-Corgi implementation hash family differs")
    if not isinstance(manifest_hashes, dict) or len(manifest_hashes) != 16:
        raise CorgiRegularSmokeError("regular-Corgi artifact-manifest family differs")

    verified: dict[str, str] = {}
    for hash_family in (authority_hashes, implementation_hashes):
        for relative, expected in sorted(hash_family.items()):
            if not isinstance(relative, str) or not isinstance(expected, str):
                raise CorgiRegularSmokeError("regular-Corgi hash record differs")
            member = _resolve_member(root, relative)
            observed = _digest(member)
            if observed != expected:
                raise CorgiRegularSmokeError(f"regular-Corgi authority drifted: {relative}")
            verified[relative] = observed

    manifests: dict[str, dict[str, Any]] = {}
    for relative, expected in sorted(manifest_hashes.items()):
        if (
            not isinstance(relative, str)
            or not relative.startswith("executions/")
            or not relative.endswith("/ARTIFACTS.json")
            or not isinstance(expected, str)
        ):
            raise CorgiRegularSmokeError("artifact-manifest allowlist differs")
        member = _resolve_member(root, relative)
        observed = _digest(member)
        if observed != expected:
            raise CorgiRegularSmokeError(f"regular-Corgi artifact manifest drifted: {relative}")
        manifest = _load_json(member)
        if manifest.get("schema_version") != "masld-bench-artifacts-v1":
            raise CorgiRegularSmokeError(f"artifact manifest schema differs: {relative}")
        manifests[relative] = manifest
        verified[relative] = observed

    checkpoint_authority = _load_json(
        _resolve_member(root, "config/artifacts/models/corgi/checkpoints.json")
    )
    release_authority = _load_json(
        _resolve_member(root, "config/artifacts/models/corgi/release.json")
    )
    exposure_authority = _load_json(
        _resolve_member(root, "config/artifacts/models/corgi/exposure_audit.json")
    )
    crosswalk_authority = _load_json(
        _resolve_member(root, "config/artifacts/models/corgi/development_crosswalk.json")
    )
    plus_authority = _load_json(
        _resolve_member(
            root,
            "config/artifacts/models/corgi_plus/release_disposition_20260824.json",
        )
    )

    identity = authority.get("model_identity", {})
    exact = authority.get("exact_contract", {})
    checkpoint = checkpoint_authority.get("checkpoint", {})
    architecture = checkpoint_authority.get("architecture_contract", {})
    release_weight = release_authority.get("weight", {})
    if (
        identity.get("model_id") != "corgi_regular"
        or identity.get("checkpoint_sha256") != CHECKPOINT_SHA256
        or checkpoint.get("sha256") != CHECKPOINT_SHA256
        or release_weight.get("sha256") != CHECKPOINT_SHA256
        or identity.get("checkpoint_size_bytes") != checkpoint.get("size_bytes")
        or identity.get("checkpoint_state_key_count")
        != checkpoint.get("state_key_count")
        or identity.get("checkpoint_state_numel")
        != architecture.get("checkpoint_state_numel")
        or identity.get("code_commit") != release_authority.get("code_commit")
        or identity.get("corgi_plus_identity_claimed") is not False
        or identity.get("relabeling_as_corgi_plus_prohibited") is not True
    ):
        raise CorgiRegularSmokeError("regular-Corgi checkpoint or identity differs")
    if (
        exact.get("sequence_input_bp") != architecture.get("input_length_bp")
        or exact.get("trans_context_features")
        != architecture.get("trans_context_features")
        or exact.get("output_channels")
        != checkpoint_authority.get("output_contract", {}).get("count")
        or exact.get("output_bins") != architecture.get("output_central_bins")
        or exact.get("output_resolution_bp")
        != architecture.get("output_resolution_bp")
        or exact.get("output_central_bp") != architecture.get("output_length_bp")
        or exact.get("atac_channel_index") != 1
        or exact.get("hgnc_gene_order_sha256")
        != checkpoint_authority.get("input_contract", {})
        .get("context_gene_order", {})
        .get("hgnc_list_sha256")
        or exact.get("ensembl_gene_order_sha256")
        != checkpoint_authority.get("input_contract", {})
        .get("context_gene_order", {})
        .get("ensembl_list_sha256")
        or exact.get("released_reference_sha256")
        != checkpoint_authority.get("input_contract", {}).get(
            "released_reference_sha256"
        )
        or exact.get("output_activation")
        != checkpoint_authority.get("output_contract", {}).get("activation")
    ):
        raise CorgiRegularSmokeError("regular-Corgi input or output contract differs")
    if (
        checkpoint_authority.get("verified_artifacts", {}).get(
            "native_numeric_parity_performed"
        )
        is not False
        or checkpoint_authority.get("license_audit", {}).get("code_status")
        != "CONFLICT_UNRESOLVED"
        or exposure_authority.get("sealed_champion_eligibility")
        == "eligible"
        or crosswalk_authority.get("findings", {})
        .get("gse296875", {})
        .get("exposure_state")
        != "clean_declared"
    ):
        raise CorgiRegularSmokeError("regular-Corgi admission boundary differs")

    plus_separation = authority.get("corgi_plus_separation", {})
    if (
        plus_authority.get("status") != "TERMINAL_BLOCKED_CURRENT_UPSTREAM_RELEASE"
        or plus_authority.get("model_identity", {}).get("relabeling_prohibited")
        is not True
        or plus_separation.get("regular_corgi_is_not_a_corgi_plus_substitute")
        is not True
        or plus_separation.get("relabeling_prohibited") is not True
    ):
        raise CorgiRegularSmokeError("Corgi+ terminal separation differs")

    runtime_path = "executions/corgi-native-runtime-l40s-r8-21065022/ARTIFACTS.json"
    _require_metadata(
        manifests[runtime_path],
        {
            "artifact_class": "corgi_native_runtime_l40s_probe",
            "regular_strict_restore": True,
            "full_524288_bp_forward": True,
            "shortened_window_backward": True,
            "native_numeric_parity_performed": False,
            "evaluator_outcomes_exposed": False,
            "status": "pass",
        },
        "L40S compatibility probe",
    )
    context_path = "executions/corgi-gse296875-context-counts-21066278/ARTIFACTS.json"
    _require_metadata(
        manifests[context_path],
        {
            "artifact_class": "masked_donor_lineage_rna_context_counts",
            "biological_unit": "donor_by_lineage_pseudobulk",
            "dataset_id": "gse296875",
            "donors": 39,
            "model_id": "corgi_regular",
            "model_training_activated": False,
            "outcomes_read": False,
            "status": "passed",
        },
        "development context",
    )
    tile_path = "executions/corgi-outcome-aligned-tiles-21068311/ARTIFACTS.json"
    _require_metadata(
        manifests[tile_path],
        {
            "artifact_class": "corgi_outcome_aligned_tile_contract",
            "dataset_id": "gse296875",
            "tiles": 6635,
            "scoreable_windows": 79956,
            "outcomes_read": False,
            "predictions_read": False,
            "status": "passed",
        },
        "development tile contract",
    )
    subset_path = "executions/corgi-tile-smoke-subset-21068454/ARTIFACTS.json"
    _require_metadata(
        manifests[subset_path],
        {
            "artifact_class": "corgi_outcome_aligned_tile_smoke_subset",
            "dataset_id": "gse296875",
            "tiles": 640,
            "windows": 7863,
            "outcomes_read": False,
            "predictions_read": False,
            "status": "passed",
        },
        "development smoke subset",
    )
    predictor_test_path = (
        "executions/corgi-outcome-aligned-predictor-test-21068997/ARTIFACTS.json"
    )
    _require_metadata(
        manifests[predictor_test_path],
        {
            "artifact_class": "corgi_outcome_aligned_predictor_cpu_test",
            "gpu_used": False,
            "model_id": "corgi_regular",
            "outcomes_read": False,
            "status": "passed",
        },
        "predictor CPU test",
    )
    preflight_path = (
        "executions/corgi-hepatocyte-tile-bundle-v3-preflight-21069014/ARTIFACTS.json"
    )
    _require_metadata(
        manifests[preflight_path],
        {
            "artifact_class": "Corgi_outcome_aligned_tile_smoke_bundle_v3_preflight",
            "gpu_used": False,
            "outcomes_read": False,
            "wrapper_sha256": implementation_hashes[
                "slurm/run_corgi_hepatocyte_tile_smoke_bundle.sbatch"
            ],
            "status": "passed",
        },
        "bundle CPU preflight",
    )
    final_path = (
        "executions/corgi-hepatocyte-tile-smoke-bundle-v3-21069656/ARTIFACTS.json"
    )
    _require_metadata(
        manifests[final_path],
        {
            "artifact_class": "Corgi_outcome_aligned_tile_smoke_bundle_v3",
            "model_id": "corgi_regular",
            "dataset_id": "gse296875",
            "lineage_id": "hepatocyte",
            "logical_tasks": 5,
            "mapper_outer_folds": list(range(5)),
            "outcome_role": "valid",
            "crossed_split_ids": [f"donor{fold}_genomic{fold}" for fold in range(5)],
            "outcomes_read": False,
            "open_champion_eligible": False,
            "status": "passed",
        },
        "completed prediction smoke",
    )

    task = authority.get("development_task", {})
    execution = authority.get("completed_execution", {})
    if (
        task.get("dataset_id") != "gse296875"
        or task.get("lineage_id") != "hepatocyte"
        or task.get("outcome_role") != "valid"
        or task.get("outer_folds") != 5
        or task.get("crossed_split_ids")
        != [f"donor{fold}_genomic{fold}" for fold in range(5)]
        or task.get("tiles_per_fold") != 128
        or task.get("total_tiles") != 640
        or task.get("context_arms") != CONTEXT_ARMS
        or task.get("regional_score_names") != SCORE_NAMES
        or task.get("model_fitted_or_adapted") is not False
        or task.get("held_atac_or_other_outcomes_used") is not False
        or task.get("test_or_sealed_features_or_labels_read") is not False
        or task.get("outcomes_read") is not False
        or task.get("evaluator_outputs_opened_for_this_audit") is not False
        or execution.get("resubmission_authorized") is not False
        or execution.get("resubmission_needed") is not False
        or execution.get("evaluation_authorized_by_this_disposition") is not False
        or execution.get("bundle_artifacts_sha256") != manifest_hashes[final_path]
        or execution.get("status") != "pass_outcome_free_prediction"
        or execution.get("logical_tasks") != 5
        or execution.get("requested_resources")
        != {
            "partition": "gpu",
            "account": "nslab",
            "qos": "nslab",
            "gpus": "l40s:1",
            "cpus": 8,
            "memory_gb": 32,
            "wall_time": "08:00:00",
        }
        or execution.get("observed_peak_memory_bytes") != 14399623168
        or execution.get("repeat_max_abs") != 0.0
        or execution.get("repeat_tolerance") != 0.001
    ):
        raise CorgiRegularSmokeError("bounded development-task disposition differs")

    fold_receipts = authority.get("fold_receipts")
    if not isinstance(fold_receipts, list) or len(fold_receipts) != 5:
        raise CorgiRegularSmokeError("prediction receipt family differs")
    final_index = _artifact_index(manifests[final_path])
    mapper_paths = [
        "executions/corgi-gse296875-mapper-fold0-21066334/ARTIFACTS.json",
        "executions/corgi-gse296875-mapper-fold1-21066335/ARTIFACTS.json",
        "executions/corgi-gse296875-mapper-fold2-21066337/ARTIFACTS.json",
        "executions/corgi-gse296875-mapper-fold3-21066336/ARTIFACTS.json",
        "executions/corgi-gse296875-mapper-fold4-21066338/ARTIFACTS.json",
    ]
    total_windows = 0
    total_valid_donors = 0
    for fold, record in enumerate(fold_receipts):
        if (
            not isinstance(record, dict)
            or record.get("fold") != fold
            or record.get("valid_fold") != (fold + 1) % 5
        ):
            raise CorgiRegularSmokeError(f"prediction fold {fold} authority differs")
        mapper_path = mapper_paths[fold]
        mapper_manifest = manifests[mapper_path]
        _require_metadata(
            mapper_manifest,
            {
                "artifact_class": "fold_fitted_masked_Corgi_context_mapper",
                "dataset_id": "gse296875",
                "model_id": "corgi_regular",
                "outcomes_read": False,
                "outer_fold": fold,
                "valid_fold": (fold + 1) % 5,
                "status": "passed",
            },
            f"mapper fold {fold}",
        )
        mapper_index = _artifact_index(mapper_manifest)
        mapper_receipt_record = mapper_index.get("receipt.json", {})
        mapper_receipt_path = str(Path(mapper_path).parent / "receipt.json")
        mapper_receipt_member = _resolve_member(root, mapper_receipt_path)
        mapper_receipt_sha = _digest(mapper_receipt_member)
        if (
            mapper_receipt_sha != mapper_receipt_record.get("sha256")
            or record.get("mapper_artifacts_sha256") != manifest_hashes[mapper_path]
        ):
            raise CorgiRegularSmokeError(f"mapper fold {fold} receipt linkage differs")
        mapper_receipt = _load_json(mapper_receipt_member)
        if (
            mapper_receipt.get("schema_version")
            != "masld-bench-corgi-masked-context-mapper-v1"
            or mapper_receipt.get("status") != "pass_outcome_free_mapper"
            or mapper_receipt.get("input_artifacts_sha256") != manifest_hashes[context_path]
            or mapper_receipt.get("outer_fold") != fold
            or mapper_receipt.get("valid_fold") != (fold + 1) % 5
            or mapper_receipt.get("held_RNA_used_for_context") is not True
            or mapper_receipt.get("held_ATAC_or_other_outcomes_used") is not False
            or mapper_receipt.get("test_or_sealed_outcomes_read") is not False
            or mapper_receipt.get("native_missing_as_zero_mapper_run") is not False
        ):
            raise CorgiRegularSmokeError(f"mapper fold {fold} safety receipt differs")

        receipt_relative = record.get("path")
        expected_receipt_sha = record.get("sha256")
        if not isinstance(receipt_relative, str) or not receipt_relative.endswith(
            f"/fold{fold}/receipt.json"
        ):
            raise CorgiRegularSmokeError(f"prediction fold {fold} receipt path differs")
        receipt_member = _resolve_member(root, receipt_relative)
        observed_receipt_sha = _digest(receipt_member)
        manifest_receipt = final_index.get(f"fold{fold}/receipt.json", {})
        if (
            observed_receipt_sha != expected_receipt_sha
            or observed_receipt_sha != manifest_receipt.get("sha256")
        ):
            raise CorgiRegularSmokeError(f"prediction fold {fold} receipt drifted")
        receipt = _load_json(receipt_member)
        if (
            receipt.get("schema_version")
            != "masld-bench-corgi-outcome-aligned-prediction-v3"
            or receipt.get("status") != "pass_outcome_free_prediction"
            or receipt.get("terminal_disposition")
            != "restricted_development_comparator_pending_license_reference_and_native_parity"
            or receipt.get("model_id") != "corgi_regular"
            or receipt.get("dataset_id") != "gse296875"
            or receipt.get("lineage_id") != "hepatocyte"
            or receipt.get("outcome_role") != "valid"
            or receipt.get("split_id") != f"donor{fold}_genomic{fold}"
            or receipt.get("mapper_outer_fold") != fold
            or receipt.get("donor_test_fold") != fold
            or receipt.get("genomic_test_fold") != fold
            or receipt.get("valid_donor_fold") != (fold + 1) % 5
            or receipt.get("valid_genomic_fold") != (fold + 1) % 5
            or receipt.get("context_arms") != CONTEXT_ARMS
            or receipt.get("regional_score_names") != SCORE_NAMES
            or receipt.get("checkpoint_sha256") != CHECKPOINT_SHA256
            or receipt.get("gpu_name") != "NVIDIA L40S"
            or receipt.get("tiles") != 128
            or receipt.get("held_RNA_used") is not True
            or receipt.get("held_ATAC_or_other_outcomes_used") is not False
            or receipt.get("test_or_sealed_features_or_labels_read") is not False
            or receipt.get("model_fitted_or_adapted") is not False
            or receipt.get("native_numeric_parity_established") is not False
            or receipt.get("open_champion_eligible") is not False
            or receipt.get("repeat_max_abs") != 0.0
            or receipt.get("count_like_tta_order")
            != "orientation_wise_softplus_then_average"
            or receipt.get("tile_contract_artifacts_sha256") != manifest_hashes[tile_path]
            or receipt.get("tile_roster_artifacts_sha256") != manifest_hashes[subset_path]
            or receipt.get("mapper_artifacts_sha256") != manifest_hashes[mapper_path]
            or receipt.get("valid_donors") != record.get("valid_donors")
            or receipt.get("scoreable_windows") != record.get("scoreable_windows")
        ):
            raise CorgiRegularSmokeError(f"prediction fold {fold} safety receipt differs")
        prediction_record = final_index.get(f"fold{fold}/regional_predictions.npz", {})
        if receipt.get("predictions_sha256") != prediction_record.get("sha256"):
            raise CorgiRegularSmokeError(f"prediction fold {fold} manifest linkage differs")
        total_windows += receipt["scoreable_windows"]
        total_valid_donors += receipt["valid_donors"]

    if total_windows != task.get("total_scoreable_windows") or total_valid_donors != 39:
        raise CorgiRegularSmokeError("bounded smoke fold totals differ")

    tournament = authority.get("tournament_disposition", {})
    if (
        tournament.get("bounded_development_prediction_complete") is not True
        or any(
            tournament.get(field) is not False
            for field in (
                "native_lane_admitted",
                "open_champion_eligible",
                "comparative_claim_eligible",
                "production_resubmission_authorized",
                "sealed_evaluation_authorized",
                "development_evaluation_authorized",
            )
        )
        or len(authority.get("remaining_blockers", [])) != 5
    ):
        raise CorgiRegularSmokeError("regular-Corgi tournament gate differs")

    return {
        "schema_version": "masld-bench-corgi-regular-smoke-audit-receipt-v1",
        "status": "pass_bounded_outcome_free_prediction_smoke",
        "authority_sha256": _digest(authority_path),
        "verified_metadata_files": len(verified),
        "admission_class": ADMISSION_CLASS,
        "native_lane_admitted": False,
        "open_champion_eligible": False,
        "resubmission_needed": False,
        "prediction_arrays_opened": False,
        "evaluator_outputs_opened": False,
        "outcomes_opened": False,
        "sealed_features_or_labels_opened": False,
        "corgi_plus_status": "TERMINAL_BLOCKED_CURRENT_UPSTREAM_RELEASE",
        "regular_corgi_relabelled_as_corgi_plus": False,
    }


def main() -> None:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=default_root)
    parser.add_argument(
        "--authority",
        type=Path,
        default=default_root
        / "config/artifacts/models/corgi/bounded_prediction_smoke_disposition_20260824.json",
    )
    parser.add_argument(
        "--require-native-admission",
        action="store_true",
        help="Fail after auditing because this disposition does not establish native admission.",
    )
    arguments = parser.parse_args()
    receipt = audit_smoke_plan(arguments.root, arguments.authority)
    print(json.dumps(receipt, indent=2, sort_keys=True))
    if arguments.require_native_admission:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
