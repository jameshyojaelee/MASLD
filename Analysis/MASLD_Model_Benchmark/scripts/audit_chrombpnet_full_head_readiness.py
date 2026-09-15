#!/usr/bin/env python3
"""Audit ChromBPNet full/no-bias semantics and readiness."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


SCHEMA_VERSION = "masld-bench-chrombpnet-full-head-production-readiness-v1"
STATUS = "RECONCILED_WAITING_FOR_FULL_RECTANGLE_LICENSE_AND_EXTERNAL_EVALUATION"
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
FOLDS = tuple(range(5))
SEEDS = (20260824, 20260825, 20260826, 20260827, 20260828)
LEGACY_SEEDS = (11, 29, 47, 71, 101)
REQUIRED_OUTPUTS = {
    "model/chrombpnet.h5",
    "model/chrombpnet_nobias.h5",
    "predictions/full_model/profile_probabilities.h5",
    "predictions/full_model/regional_counts.tsv",
    "predictions/nobias_model/profile_probabilities.h5",
    "predictions/nobias_model/regional_counts.tsv",
}
EXPECTED_LINEAGE_COUNTS = {
    "cholangiocyte": {"complete": 5, "missing": 20},
    "fibroblast": {"complete": 4, "missing": 21},
    "hepatocyte": {"complete": 1, "missing": 24},
    "macrophage": {"complete": 4, "missing": 21},
    "t_cell": {"complete": 5, "missing": 20},
}


class ChromBPNetReadinessError(ValueError):
    """Raised when a frozen authority or a fail-closed check differs."""


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
        raise ChromBPNetReadinessError(f"JSON authority differs: {path}") from error
    if not isinstance(value, dict):
        raise ChromBPNetReadinessError(f"JSON authority is not an object: {path}")
    return value


def _resolve_file(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ChromBPNetReadinessError("authority path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ChromBPNetReadinessError(f"authority is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise ChromBPNetReadinessError(
            f"authority escapes or is absent: {relative}"
        ) from error
    if not resolved.is_file():
        raise ChromBPNetReadinessError(f"authority is not a file: {relative}")
    return resolved


def _resolve_directory(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise ChromBPNetReadinessError("artifact path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise ChromBPNetReadinessError(f"artifact is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise ChromBPNetReadinessError(
            f"artifact escapes or is absent: {relative}"
        ) from error
    if not resolved.is_dir():
        raise ChromBPNetReadinessError(f"artifact is not a directory: {relative}")
    return resolved


def _require_file(root: Path, record: Mapping[str, Any]) -> Path:
    path = _resolve_file(root, str(record.get("path", "")))
    if _digest(path) != record.get("sha256"):
        raise ChromBPNetReadinessError(f"file authority differs: {path}")
    return path


def _require_frozen_artifact(
    root: Path,
    relative: str,
    expected_sha256: str,
    expected_class: str,
) -> dict[str, Any]:
    artifact_root = _resolve_directory(root, relative)
    try:
        manifest = verify_frozen_tree(artifact_root)
    except ArtifactError as error:
        raise ChromBPNetReadinessError(
            f"frozen artifact differs: {relative}"
        ) from error
    if (
        _digest(artifact_root / "ARTIFACTS.json") != expected_sha256
        or manifest.get("metadata", {}).get("artifact_class") != expected_class
        or manifest.get("metadata", {}).get("status") != "passed"
    ):
        raise ChromBPNetReadinessError(f"frozen artifact identity differs: {relative}")
    return dict(manifest)


def _require_manifest_only(
    root: Path,
    relative: str,
    expected_sha256: str,
    expected_class: str,
) -> tuple[Path, dict[str, Any]]:
    artifact_root = _resolve_directory(root, relative)
    manifest_path = artifact_root / "ARTIFACTS.json"
    manifest = _load_json(manifest_path)
    if (
        _digest(manifest_path) != expected_sha256
        or manifest.get("metadata", {}).get("artifact_class") != expected_class
    ):
        raise ChromBPNetReadinessError(f"artifact manifest identity differs: {relative}")
    return artifact_root, manifest


def _artifact_rows(manifest: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    rows = manifest.get("artifacts", [])
    if not isinstance(rows, list):
        raise ChromBPNetReadinessError("artifact inventory is not a list")
    indexed = {str(row.get("path")): dict(row) for row in rows if isinstance(row, dict)}
    if len(indexed) != len(rows):
        raise ChromBPNetReadinessError("artifact inventory has duplicate or invalid rows")
    return indexed


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _assert_false(mapping: Mapping[str, Any], fields: tuple[str, ...]) -> None:
    for field in fields:
        if mapping.get(field) is not False:
            raise ChromBPNetReadinessError(f"fail-closed field opened: {field}")


def _audit_head_semantics(root: Path, contract: Mapping[str, Any]) -> None:
    authorities = contract.get("frozen_authorities", {})
    checkpoint_path = _require_file(root, authorities.get("checkpoint_authority", {}))
    crosswalk_path = _require_file(root, authorities.get("development_crosswalk", {}))
    exposure_path = _require_file(root, authorities.get("exposure_audit", {}))
    upstream_path = _require_file(root, authorities.get("upstream_full_head_source", {}))
    prediction_path = _require_file(root, authorities.get("project_prediction_source", {}))
    _require_file(root, authorities.get("production_wrapper", {}))

    checkpoints = _load_json(checkpoint_path)
    crosswalk = _load_json(crosswalk_path)
    exposure = _load_json(exposure_path)
    if (
        checkpoints.get("code_license") != "MIT"
        or checkpoints.get("architecture_contract", {}).get("outputs")
        != ["profile_logits_1000bp", "log_total_count"]
        or checkpoints.get("architecture_contract", {}).get("assay_model", {}).get(
            "input_length_bp"
        )
        != 2114
        or checkpoints.get("architecture_contract", {}).get("assay_model", {}).get(
            "output_length_bp"
        )
        != 1000
        or checkpoints.get("inference_contract", {}).get("atac_available_at_inference")
        is not False
        or "bias-free TF model" not in checkpoints.get("inference_contract", {}).get(
            "profile_prediction", ""
        )
    ):
        raise ChromBPNetReadinessError("checkpoint or inference semantics differ")
    if (
        tuple(crosswalk.get("local_model_roster", {}).get("cell_states", ()))
        != LINEAGES
        or tuple(crosswalk.get("local_model_roster", {}).get("seeds", ()))
        != LEGACY_SEEDS
        or crosswalk.get("findings", {}).get("gse296875", {}).get("exposure_state")
        != "continual_seen"
        or crosswalk.get("findings", {}).get("gse289173", {}).get("exposure_state")
        != "target_label_unexposed"
        or exposure.get("sealed_champion_eligibility")
        != "conditional_after_runtime_data_and_local_artifact_locks_as_a_sequence_accessibility_component_not_a_standalone_signed_gene_effect_model"
    ):
        raise ChromBPNetReadinessError("crosswalk or exposure boundary differs")

    upstream = upstream_path.read_text(encoding="utf-8")
    upstream_tokens = (
        'name="model_wo_bias"',
        'profile_out = Add(name="logits_profile_predictions")([output_wo_bias[0],bias_output[0]])',
        "Concatenate(axis=-1)([output_wo_bias[1], bias_output[1]])",
        "tf.math.reduce_logsumexp(x, axis=-1, keepdims=True)",
        'model.get_layer("model_wo_bias").output',
        'model_without_bias.save(output_prefix+"_nobias.h5")',
    )
    if any(token not in upstream for token in upstream_tokens):
        raise ChromBPNetReadinessError("upstream full/no-bias head semantics differ")

    prediction = prediction_path.read_text(encoding="utf-8")
    prediction_tokens = (
        'if transform == "log1p_absolute":',
        "return np.maximum(np.expm1(logcounts), 0.0)",
        'if transform == "log_component":',
        "return np.exp(logcounts)",
        "reverse_probability = softmax(reverse_profile_raw)[:, ::-1]",
        "0.5 * (forward_probability + reverse_probability)",
        "regional_mass_from_logscore(forward_logcount, count_transform)",
        "regional_mass_from_logscore(reverse_logcount, count_transform)",
    )
    if any(token not in prediction for token in prediction_tokens):
        raise ChromBPNetReadinessError("project inference semantics differ")

    head = contract.get("scientific_head_contract", {})
    if (
        head.get("input_shape") != [None, 2114, 4]
        or head.get("profile_output_shape") != [None, 1000]
        or head.get("count_output_shape") != [None, 1]
        or head.get("one_fit_two_views") is not True
        or head.get("primary_biological_view") != "chrombpnet_nobias"
        or head.get("assay_qc_view") != "chrombpnet_full"
        or head.get("full_view_is_independent_architecture") is not False
        or head.get("full_count_inverse_transform") != "expm1_then_nonnegative_clip"
        or head.get("nobias_count_inverse_transform") != "exp"
        or head.get("standalone_signed_gene_effect_claim_allowed") is not False
    ):
        raise ChromBPNetReadinessError("scientific head contract differs")


def _audit_reference_fit(root: Path, contract: Mapping[str, Any]) -> None:
    reference = contract.get("reference_fit", {})
    artifact_root, manifest = _require_manifest_only(
        root,
        str(reference.get("artifact_root", "")),
        str(reference.get("artifacts_sha256", "")),
        str(reference.get("expected_artifact_class", "")),
    )
    metadata = manifest.get("metadata", {})
    if (
        metadata.get("status") != "passed"
        or metadata.get("dataset_id") != reference.get("dataset_id")
        or metadata.get("lineage_id") != reference.get("lineage_id")
        or metadata.get("split_id") != reference.get("split_id")
        or metadata.get("seed") != reference.get("seed")
        or metadata.get("primary_biological_score") != "nobias_model"
        or metadata.get("assay_qc_score") != "full_model"
        or metadata.get("source_bias_role") != "restricted_smoke_comparator"
        or metadata.get("champion_eligible") is not False
    ):
        raise ChromBPNetReadinessError("reference fit metadata differs")
    _assert_false(metadata, ("benchmark_metrics_calculated", "test_outcomes_used"))

    inventory = _artifact_rows(manifest)
    if not REQUIRED_OUTPUTS.issubset(inventory):
        raise ChromBPNetReadinessError("reference fit lacks full/no-bias outputs")
    if inventory["model/chrombpnet.h5"].get("sha256") != reference.get(
        "full_checkpoint_sha256"
    ) or inventory["model/chrombpnet_nobias.h5"].get("sha256") != reference.get(
        "nobias_checkpoint_sha256"
    ):
        raise ChromBPNetReadinessError("reference checkpoint identity differs")
    if set(reference.get("required_artifact_paths", ())) != REQUIRED_OUTPUTS:
        raise ChromBPNetReadinessError("reference output contract differs")

    validation_row = inventory.get("validation/summary.json", {})
    validation_path = artifact_root / "validation/summary.json"
    if (
        not validation_path.is_file()
        or _digest(validation_path) != validation_row.get("sha256")
    ):
        raise ChromBPNetReadinessError("reference validation summary differs")
    validation = _load_json(validation_path)
    if (
        validation.get("observed_held_donor_atac_used_for_inference") is not False
        or validation.get("primary_biological_score") != "nobias_model"
        or validation.get("assay_qc_score") != "full_model"
        or validation.get("champion_eligible") is not False
    ):
        raise ChromBPNetReadinessError("reference validation boundary differs")
    _assert_false(validation, ("benchmark_metrics_calculated", "test_outcomes_used"))


def _audit_rectangle(root: Path, contract: Mapping[str, Any]) -> None:
    rectangle = contract.get("five_seed_rectangle", {})
    contract_path = _resolve_file(root, str(rectangle.get("contract_path", "")))
    if _digest(contract_path) != rectangle.get("contract_sha256"):
        raise ChromBPNetReadinessError("five-seed contract differs")
    campaign = _load_json(contract_path)
    if (
        tuple(campaign.get("lineages", ())) != LINEAGES
        or tuple(campaign.get("diagonal_outer_folds", ())) != FOLDS
        or tuple(campaign.get("fixed_seeds", ())) != SEEDS
        or campaign.get("prediction_views", {}).get("chrombpnet")
        != ["chrombpnet_full", "chrombpnet_nobias"]
        or campaign.get("fit_count_contract", {}).get(
            "chrombpnet_full_and_nobias_share_one_fit"
        )
        is not True
        or campaign.get("execution", {}).get("submission_authority")
        != "dispatcher_only_no_manual_gpu_sbatch"
    ):
        raise ChromBPNetReadinessError("five-seed campaign semantics differ")
    _assert_false(
        campaign.get("claims", {}),
        (
            "champion_claim_allowed",
            "universal_claim_allowed",
            "cross_model_ranking_allowed",
            "post_hoc_unit_selection_allowed",
        ),
    )

    plan = _require_frozen_artifact(
        root,
        str(rectangle.get("admission_artifact_root", "")),
        str(rectangle.get("admission_artifacts_sha256", "")),
        "sequence_task_native_five_seed_rectangle_admission",
    )
    evaluator = _require_frozen_artifact(
        root,
        str(rectangle.get("evaluator_readiness_artifact_root", "")),
        str(rectangle.get("evaluator_readiness_artifacts_sha256", "")),
        "sequence_task_native_five_seed_evaluator_readiness",
    )
    if (
        plan.get("metadata", {}).get("prediction_or_metric_values_read") is not False
        or evaluator.get("metadata", {}).get("gate_status")
        != "waiting_for_full_rectangle"
        or evaluator.get("metadata", {}).get("full_compatible_rectangle") is not False
        or evaluator.get("metadata", {}).get("outcome_access_authorized") is not False
        or evaluator.get("metadata", {}).get("partial_cross_model_ranking_authorized")
        is not False
    ):
        raise ChromBPNetReadinessError("rectangle or evaluator firewall differs")

    plan_root = _resolve_directory(root, str(rectangle.get("admission_artifact_root", "")))
    rows = [
        row
        for row in _read_tsv(plan_root / "plan/expected_fit_matrix.tsv")
        if row.get("model_id") == "chrombpnet"
    ]
    expected_keys = {
        (lineage, fold, seed)
        for lineage in LINEAGES
        for fold in FOLDS
        for seed in SEEDS
    }
    actual_keys = {
        (row["lineage_id"], int(row["outer_fold"]), int(row["seed"]))
        for row in rows
    }
    if len(rows) != 125 or actual_keys != expected_keys:
        raise ChromBPNetReadinessError("ChromBPNet rectangle is not 5x5x5")
    if any(
        row.get("split_id") != f"donor{row['outer_fold']}_genomic{row['outer_fold']}"
        or row.get("prediction_views") != "chrombpnet_full,chrombpnet_nobias"
        or row.get("fit_state") not in {"complete_reusable", "missing"}
        for row in rows
    ):
        raise ChromBPNetReadinessError("ChromBPNet row semantics differ")

    states = Counter(row["fit_state"] for row in rows)
    lineage_states: dict[str, Counter[str]] = defaultdict(Counter)
    for row in rows:
        lineage_states[row["lineage_id"]][row["fit_state"]] += 1
    observed_lineage_counts = {
        lineage: {
            "complete": lineage_states[lineage]["complete_reusable"],
            "missing": lineage_states[lineage]["missing"],
        }
        for lineage in LINEAGES
    }
    complete_contexts = sorted(
        [
            [row["lineage_id"], int(row["outer_fold"]), int(row["seed"])]
            for row in rows
            if row["fit_state"] == "complete_reusable"
        ]
    )
    if (
        states != Counter({"missing": 106, "complete_reusable": 19})
        or observed_lineage_counts != EXPECTED_LINEAGE_COUNTS
        or observed_lineage_counts != rectangle.get("lineage_fit_counts")
        or complete_contexts != sorted(rectangle.get("complete_contexts", []))
        or rectangle.get("expected_fits") != 125
        or rectangle.get("complete_prediction_views_at_freeze") != 38
        or rectangle.get("missing_prediction_views_at_freeze") != 212
        or rectangle.get("full_compatible_rectangle") is not False
        or rectangle.get("partial_evaluator_access_allowed") is not False
    ):
        raise ChromBPNetReadinessError("ChromBPNet frozen rectangle counts differ")

    history = contract.get("authority_history", {})
    if (
        tuple(history.get("legacy_crosswalk_seeds", ())) != LEGACY_SEEDS
        or tuple(history.get("production_rectangle_seeds", ())) != SEEDS
        or history.get("production_seed_authority")
        != "later_frozen_five_seed_rectangle_contract"
        or history.get("legacy_crosswalk_rewritten_by_this_reconciliation") is not False
        or history.get("seed_sets_may_be_mixed") is not False
    ):
        raise ChromBPNetReadinessError("seed-authority history differs")


def _audit_license_and_dispatch(root: Path, contract: Mapping[str, Any]) -> None:
    license_gate = contract.get("license_and_release_gate", {})
    bias_root, bias = _require_manifest_only(
        root,
        str(license_gate.get("bias_artifact_root", "")),
        str(license_gate.get("bias_artifacts_sha256", "")),
        "chrombpnet_pretrained_bias_models",
    )
    bias_metadata = bias.get("metadata", {})
    inventory = _artifact_rows(bias)
    validation_path = bias_root / "validation.json"
    record_path = bias_root / "zenodo-record-7443683.json"
    if (
        _digest(validation_path) != inventory.get("validation.json", {}).get("sha256")
        or _digest(record_path)
        != inventory.get("zenodo-record-7443683.json", {}).get("sha256")
    ):
        raise ChromBPNetReadinessError("bias license evidence differs")
    validation = _load_json(validation_path)
    record = _load_json(record_path)
    if (
        bias_metadata.get("weight_license_status") != "unknown_not_redistributable"
        or bias_metadata.get("allowed_role") != "restricted_smoke_comparator"
        or validation.get("zenodo_license_field", {}).get("id") != "cc-by-4.0"
        or record.get("metadata", {}).get("license", {}).get("id") != "cc-by-4.0"
        or license_gate.get("authority_contradiction_present") is not True
        or license_gate.get(
            "formal_terms_and_derivative_attribution_reconciliation_required"
        )
        is not True
        or license_gate.get("automatic_authority_rewrite_allowed") is not False
        or license_gate.get("full_view_weight_redistribution_allowed") is not False
        or license_gate.get("nobias_view_weight_redistribution_allowed") is not False
        or license_gate.get("open_champion_eligible_now") is not False
    ):
        raise ChromBPNetReadinessError("license and release gate differs")

    dispatcher = contract.get("dispatcher_gate", {})
    _require_frozen_artifact(
        root,
        str(dispatcher.get("execution_correction_artifact_root", "")),
        str(dispatcher.get("execution_correction_artifacts_sha256", "")),
        "sequence_task_native_five_seed_execution_correction",
    )
    queue = _require_frozen_artifact(
        root,
        str(dispatcher.get("corrected_queue_artifact_root", "")),
        str(dispatcher.get("corrected_queue_artifacts_sha256", "")),
        "sequence_task_native_five_seed_corrected_queue_admission",
    )
    live = _require_frozen_artifact(
        root,
        str(dispatcher.get("live_readiness_artifact_root", "")),
        str(dispatcher.get("live_readiness_artifacts_sha256", "")),
        "sequence_task_native_five_seed_live_dispatch_readiness",
    )
    next_item_path = _resolve_file(root, str(dispatcher.get("next_corrected_queue_item_path", "")))
    disabled_item_path = _resolve_file(root, str(dispatcher.get("disabled_v1_queue_item_path", "")))
    if (
        _digest(next_item_path) != dispatcher.get("next_corrected_queue_item_sha256")
        or _digest(disabled_item_path) != dispatcher.get("disabled_v1_queue_item_sha256")
    ):
        raise ChromBPNetReadinessError("queue item identity differs")
    next_item = _load_json(next_item_path)
    disabled_item = _load_json(disabled_item_path)
    queue_metadata = queue.get("metadata", {})
    live_metadata = live.get("metadata", {})
    if (
        next_item.get("enabled") is not True
        or next_item.get("bundle_id") != "model-training-801"
        or next_item.get("exports", {}).get("RECTANGLE_BUNDLE_ID") != "model-training-701"
        or next_item.get("logical_tasks") != 15
        or disabled_item.get("enabled") is not False
        or queue_metadata.get("submission_authority")
        != "dispatcher_only_no_manual_gpu_sbatch"
        or queue_metadata.get("manual_gpu_sbatch_allowed") is not False
        or queue_metadata.get("qos") != "nslab"
        or live_metadata.get("submission_authority") != "central_dispatcher_only"
        or live_metadata.get("readiness") != "v2_ready_waiting_for_gpu_cap"
        or live_metadata.get("gpu_pending") != 1
        or live_metadata.get("innovation_used") is not False
        or dispatcher.get("valid_next_gpu_unit_already_registered") is not True
        or dispatcher.get("new_gpu_unit_registration_needed") is not False
        or dispatcher.get("direct_gpu_submission_authorized") is not False
    ):
        raise ChromBPNetReadinessError("dispatcher production gate differs")


def audit_readiness(root: Path, contract_path: Path) -> dict[str, Any]:
    """Verify outcome-blind metadata without opening prediction or outcome values."""

    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    try:
        contract_path.relative_to(root)
    except ValueError as error:
        raise ChromBPNetReadinessError("contract escapes benchmark root") from error
    contract = _load_json(contract_path)
    if contract.get("schema_version") != SCHEMA_VERSION or contract.get("status") != STATUS:
        raise ChromBPNetReadinessError("readiness schema or status differs")
    _audit_head_semantics(root, contract)
    _audit_reference_fit(root, contract)
    _audit_rectangle(root, contract)
    _audit_license_and_dispatch(root, contract)

    claims = contract.get("exposure_and_claim_gate", {})
    firewall = contract.get("firewall", {})
    _assert_false(
        claims,
        (
            "target_label_unexposed_equals_champion_eligible",
            "external_evaluation_complete",
            "sealed_evaluation_complete",
            "development_champion_claim_allowed",
            "external_champion_claim_allowed",
            "universal_claim_allowed",
        ),
    )
    _assert_false(
        firewall,
        (
            "prediction_values_may_be_read",
            "metric_values_may_be_read",
            "development_outcomes_may_be_read",
            "sealed_features_may_be_read",
            "sealed_labels_or_outcomes_may_be_read",
            "thresholds_may_be_changed",
            "gpu_job_may_be_submitted_by_this_audit",
        ),
    )

    return {
        "schema_version": "masld-bench-chrombpnet-full-head-readiness-receipt-v1",
        "status": "pass_full_head_reconciled_waiting_for_rectangle_and_license",
        "contract_sha256": _digest(contract_path),
        "full_and_nobias_share_one_fit": True,
        "primary_biological_view": "chrombpnet_nobias",
        "assay_qc_view": "chrombpnet_full",
        "full_view_is_independent_architecture": False,
        "chrombpnet_expected_fits": 125,
        "chrombpnet_complete_fits_at_freeze": 19,
        "chrombpnet_missing_fits_at_freeze": 106,
        "chrombpnet_expected_prediction_views": 250,
        "chrombpnet_complete_prediction_views_at_freeze": 38,
        "chrombpnet_missing_prediction_views_at_freeze": 212,
        "full_compatible_rectangle": False,
        "partial_evaluator_access_authorized": False,
        "legacy_seed_authority_superseded_for_rectangle": True,
        "license_authority_contradiction_present": True,
        "weight_redistribution_authorized": False,
        "open_champion_eligible_now": False,
        "target_label_unexposed_equals_champion_eligible": False,
        "valid_next_gpu_unit_already_registered": True,
        "next_dispatcher_bundle_id": "model-training-801",
        "new_gpu_unit_registration_needed": False,
        "direct_gpu_submission_authorized": False,
        "prediction_values_read": False,
        "metric_values_read": False,
        "development_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_labels_or_outcomes_read": False,
        "thresholds_changed": False,
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
