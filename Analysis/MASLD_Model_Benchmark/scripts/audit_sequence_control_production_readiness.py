#!/usr/bin/env python3
"""Audit outcome-blind CNN and transformer readiness."""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any, Mapping

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


SCHEMA_VERSION = "masld-bench-sequence-control-production-readiness-v1"
STATUS = (
    "RECONCILED_WAITING_FOR_FULL_RECTANGLES_PROJECT_LICENSE_AND_EXTERNAL_EVALUATION"
)
MODEL_IDS = ("sequence_cnn_control", "sequence_transformer_control")
LINEAGES = ("cholangiocyte", "fibroblast", "hepatocyte", "macrophage", "t_cell")
FOLDS = tuple(range(5))
SEEDS = (20260824, 20260825, 20260826, 20260827, 20260828)
LEGACY_SEEDS = (1103, 2909, 4721, 6673, 8111)
EXPECTED_COUNTS = {
    "sequence_cnn_control": {
        "complete": 12,
        "missing": 113,
        "lineage": {
            "cholangiocyte": 4,
            "fibroblast": 2,
            "hepatocyte": 1,
            "macrophage": 2,
            "t_cell": 3,
        },
        "fold": {0: 4, 1: 3, 2: 3, 3: 1, 4: 1},
        "seed": {20260824: 12, 20260825: 0, 20260826: 0, 20260827: 0, 20260828: 0},
    },
    "sequence_transformer_control": {
        "complete": 11,
        "missing": 114,
        "lineage": {
            "cholangiocyte": 3,
            "fibroblast": 2,
            "hepatocyte": 1,
            "macrophage": 3,
            "t_cell": 2,
        },
        "fold": {0: 4, 1: 3, 2: 2, 3: 1, 4: 1},
        "seed": {20260824: 10, 20260825: 1, 20260826: 0, 20260827: 0, 20260828: 0},
    },
}


class SequenceControlReadinessError(ValueError):
    """Raised when a frozen sequence-control authority differs."""


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
        raise SequenceControlReadinessError(f"JSON authority differs: {path}") from error
    if not isinstance(value, dict):
        raise SequenceControlReadinessError(f"JSON authority is not an object: {path}")
    return value


def _load_toml(path: Path) -> dict[str, Any]:
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise SequenceControlReadinessError(f"TOML authority differs: {path}") from error
    if not isinstance(value, dict):
        raise SequenceControlReadinessError(f"TOML authority is not a table: {path}")
    return value


def _resolve_file(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise SequenceControlReadinessError("authority path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise SequenceControlReadinessError(f"authority is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise SequenceControlReadinessError(
            f"authority escapes or is absent: {relative}"
        ) from error
    if not resolved.is_file():
        raise SequenceControlReadinessError(f"authority is not a file: {relative}")
    return resolved


def _resolve_directory(root: Path, relative: str) -> Path:
    if not relative or Path(relative).is_absolute():
        raise SequenceControlReadinessError("artifact path is not relative")
    candidate = root / relative
    if candidate.is_symlink():
        raise SequenceControlReadinessError(f"artifact is a symlink: {relative}")
    try:
        resolved = candidate.resolve(strict=True)
        resolved.relative_to(root)
    except (OSError, ValueError) as error:
        raise SequenceControlReadinessError(
            f"artifact escapes or is absent: {relative}"
        ) from error
    if not resolved.is_dir():
        raise SequenceControlReadinessError(f"artifact is not a directory: {relative}")
    return resolved


def _require_file(root: Path, record: Mapping[str, Any]) -> Path:
    path = _resolve_file(root, str(record.get("path", "")))
    if _digest(path) != record.get("sha256"):
        raise SequenceControlReadinessError(f"file authority differs: {path}")
    return path


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
        raise SequenceControlReadinessError(f"artifact manifest differs: {relative}")
    return artifact_root, manifest


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
        raise SequenceControlReadinessError(f"frozen artifact differs: {relative}") from error
    if (
        _digest(artifact_root / "ARTIFACTS.json") != expected_sha256
        or manifest.get("metadata", {}).get("artifact_class") != expected_class
    ):
        raise SequenceControlReadinessError(f"frozen artifact identity differs: {relative}")
    return artifact_root, dict(manifest)


def _inventory(manifest: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
    rows = manifest.get("artifacts", [])
    if not isinstance(rows, list):
        raise SequenceControlReadinessError("artifact inventory is not a list")
    indexed = {str(row.get("path")): dict(row) for row in rows if isinstance(row, dict)}
    if len(indexed) != len(rows):
        raise SequenceControlReadinessError("artifact inventory has invalid or duplicate rows")
    return indexed


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _assert_false(mapping: Mapping[str, Any], fields: tuple[str, ...]) -> None:
    for field in fields:
        if mapping.get(field) is not False:
            raise SequenceControlReadinessError(f"fail-closed field opened: {field}")


def _audit_architectures(root: Path, contract: Mapping[str, Any]) -> None:
    authorities = contract.get("frozen_source_authorities", {})
    config_path = _require_file(root, authorities.get("architecture_config", {}))
    architecture_path = _require_file(root, authorities.get("architecture_source", {}))
    _require_file(root, authorities.get("training_input_source", {}))
    training_path = _require_file(root, authorities.get("training_source", {}))
    prediction_path = _require_file(root, authorities.get("prediction_source", {}))
    driver_path = _require_file(root, authorities.get("five_seed_driver", {}))
    _require_file(root, authorities.get("cnn_wrapper", {}))
    _require_file(root, authorities.get("transformer_wrapper", {}))
    _require_file(root, authorities.get("corrected_bundle_wrapper", {}))

    config = _load_json(config_path)
    shared = config.get("shared_contract", {})
    architectures = config.get("architectures", {})
    if (
        set(architectures) != set(MODEL_IDS)
        or shared.get("input_length") != 2114
        or shared.get("profile_output_length") != 1000
        or shared.get("from_scratch") is not True
        or shared.get("held_donor_atac_available_to_model") is not False
        or shared.get("test_outcomes_used") is not False
        or tuple(shared.get("seed_authority", {}).get("finalist", ())) != LEGACY_SEEDS
        or config.get("parameter_budget", {}).get("maximum_pairwise_ratio") != 1.25
    ):
        raise SequenceControlReadinessError("architecture config differs")
    if architectures["sequence_cnn_control"] != {
        "blocks": 4,
        "channels": 256,
        "dropout": 0.1,
        "kernel_size": 7,
        "stem_kernel_size": 21,
    } or architectures["sequence_transformer_control"] != {
        "attention_dropout": 0.1,
        "blocks": 8,
        "dropout": 0.1,
        "ffn_width": 768,
        "heads": 8,
        "profile_kernel_size": 7,
        "stem_kernel_size": 15,
        "token_stride": 4,
        "width": 192,
    }:
        raise SequenceControlReadinessError("CNN or transformer architecture differs")

    architecture = architecture_path.read_text(encoding="utf-8")
    required_architecture_tokens = (
        "def _build_cnn(",
        "def _build_transformer(",
        'name="sequence_cnn_control"',
        'name="sequence_transformer_control"',
        "tf.keras.layers.MultiHeadAttention(",
        "tf.keras.layers.GlobalAveragePooling1D(name=\"cnn_global_average\")",
        "loss=[multinomial_nll, \"mse\"]",
    )
    if any(token not in architecture for token in required_architecture_tokens):
        raise SequenceControlReadinessError("architecture source semantics differ")
    training = training_path.read_text(encoding="utf-8")
    prediction = prediction_path.read_text(encoding="utf-8")
    driver = driver_path.read_text(encoding="utf-8")
    if any(
        token not in training
        for token in (
            'choices=MODEL_IDS',
            '"from_scratch": True',
            '"held_donor_atac_exposed": False',
            '"test_outcomes_used": False',
        )
    ) or any(
        token not in prediction
        for token in (
            "reverse_complement_one_hot",
            "reverse_probability = softmax(reverse_profile_raw)[:, ::-1]",
            "0.5 * (forward_probability + reverse_probability)",
            "mass_from_log1p",
        )
    ) or any(
        token not in driver
        for token in (
            "sequence_cnn_control",
            "sequence_transformer_control",
            "20260824|20260825|20260826|20260827|20260828",
        )
    ):
        raise SequenceControlReadinessError("training, prediction, or seed semantics differ")

    identity = contract.get("architecture_identity", {})
    cnn = identity.get("sequence_cnn_control", {})
    transformer = identity.get("sequence_transformer_control", {})
    shared_identity = identity.get("shared", {})
    ratio = max(cnn.get("reference_parameter_count", 0), transformer.get("reference_parameter_count", 0)) / min(
        cnn.get("reference_parameter_count", 0), transformer.get("reference_parameter_count", 0)
    )
    if (
        cnn.get("profile_receptive_field_bp") != 69
        or cnn.get("count_receptive_field_bp") != 2114
        or transformer.get("profile_receptive_field_bp") != 2114
        or transformer.get("count_receptive_field_bp") != 2114
        or ratio > identity.get("maximum_pairwise_parameter_ratio", 0)
        or abs(ratio - identity.get("observed_pairwise_parameter_ratio", 0)) > 1.0e-9
        or shared_identity.get("same_architecture_family") is not False
        or shared_identity.get("checkpoint_sharing_allowed") is not False
        or shared_identity.get("separate_fit_and_prediction_identity_required") is not True
        or identity.get("controls_may_be_pooled_as_one_architecture") is not False
    ):
        raise SequenceControlReadinessError("reconciled architecture identity differs")


def _audit_reference_artifacts(root: Path, contract: Mapping[str, Any]) -> None:
    for model_id in MODEL_IDS:
        record = contract.get("reference_artifacts", {}).get(model_id, {})
        artifact_root, manifest = _require_manifest_only(
            root,
            str(record.get("artifact_root", "")),
            str(record.get("artifacts_sha256", "")),
            str(record.get("expected_artifact_class", "")),
        )
        metadata = manifest.get("metadata", {})
        if (
            metadata.get("status") != "passed"
            or metadata.get("model_id") != model_id
            or metadata.get("dataset_id") != record.get("dataset_id")
            or metadata.get("lineage_id") != record.get("lineage_id")
            or metadata.get("split_id") != record.get("split_id")
            or metadata.get("seed") != record.get("seed")
            or metadata.get("held_donor_atac_exposed") is not False
            or metadata.get("test_outcomes_used") is not False
        ):
            raise SequenceControlReadinessError(f"reference metadata differs: {model_id}")
        inventory = _inventory(manifest)
        checkpoint_member = str(record.get("checkpoint_member", ""))
        summary_member = str(record.get("validation_summary_member", ""))
        if (
            inventory.get(checkpoint_member, {}).get("sha256")
            != record.get("checkpoint_sha256")
            or inventory.get(summary_member, {}).get("sha256")
            != record.get("validation_summary_sha256")
        ):
            raise SequenceControlReadinessError(f"reference inventory differs: {model_id}")
        summary_path = artifact_root / summary_member
        if _digest(summary_path) != record.get("validation_summary_sha256"):
            raise SequenceControlReadinessError(f"reference summary bytes differ: {model_id}")
        summary = _load_json(summary_path)
        expected_parameters = contract["architecture_identity"][model_id][
            "reference_parameter_count"
        ]
        if (
            summary.get("model_id") != model_id
            or summary.get("input_length") != 2114
            or summary.get("output_length") != 1000
            or summary.get("parameter_count") != expected_parameters
            or summary.get("reverse_complement_tta") is not True
            or summary.get("ready_for_independent_evaluator_expansion") is not True
            or summary.get("benchmark_metrics_calculated") is not False
            or summary.get("held_donor_atac_exposed") is not False
            or summary.get("test_outcomes_used") is not False
        ):
            raise SequenceControlReadinessError(f"reference validation differs: {model_id}")


def _audit_authority_history(root: Path, contract: Mapping[str, Any]) -> None:
    authorities = contract.get("frozen_source_authorities", {})
    checkpoints = _load_json(_require_file(root, authorities.get("legacy_checkpoint_authority", {})))
    crosswalk = _load_json(_require_file(root, authorities.get("development_crosswalk", {})))
    exposure = _load_json(_require_file(root, authorities.get("exposure_audit", {})))
    cnn = checkpoints.get("models", {}).get("sequence_cnn_control", {})
    transformer = checkpoints.get("models", {}).get("sequence_transformer_control", {})
    history = contract.get("authority_history", {})
    if (
        cnn.get("implementation_status")
        != history.get("legacy_registry_CNN_state")
        or transformer.get("implementation_status")
        != history.get("legacy_registry_transformer_state")
        or crosswalk.get("crosswalk_scope", {}).get("status")
        != "blocked_pending_implementation_runtime_and_fixed_fixtures"
        or exposure.get("license_and_release_disposition", {}).get(
            "project_implementation_license"
        )
        != "UNRESOLVED_for_sequence_cnn_control_despite_exact_code_artifact"
        or tuple(history.get("legacy_architecture_config_finalist_seeds", ()))
        != LEGACY_SEEDS
        or tuple(history.get("production_rectangle_seeds", ())) != SEEDS
        or history.get("production_seed_authority")
        != "later_frozen_five_seed_rectangle_contract"
        or history.get("legacy_registry_rewritten_by_this_reconciliation") is not False
        or history.get("later_artifacts_upgrade_open_release_or_champion_eligibility")
        is not False
        or history.get("seed_sets_may_be_mixed") is not False
    ):
        raise SequenceControlReadinessError("authority history differs")


def _audit_split_isolation(root: Path, contract: Mapping[str, Any]) -> None:
    isolation = contract.get("donor_and_genomic_isolation", {})
    membership_root, membership = _require_manifest_only(
        root,
        str(isolation.get("fragment_membership_artifact_root", "")),
        str(isolation.get("fragment_membership_artifacts_sha256", "")),
        "gse296875_fragment_membership",
    )
    membership_meta = membership.get("metadata", {})
    membership_inventory = _inventory(membership)
    donor_path = membership_root / "donor_folds.tsv"
    if (
        membership_meta.get("donors") != 39
        or membership_meta.get("cells") != 68398
        or membership_meta.get("primary_lineages") != 5
        or _digest(donor_path)
        != membership_inventory.get("donor_folds.tsv", {}).get("sha256")
    ):
        raise SequenceControlReadinessError("donor membership authority differs")
    donor_rows = _read_tsv(donor_path)
    donor_ids = [row.get("donor_id") for row in donor_rows]
    donor_folds = [int(row["outer_fold"]) for row in donor_rows]
    if (
        len(donor_rows) != 39
        or len(set(donor_ids)) != 39
        or Counter(donor_folds)
        != Counter(dict(enumerate(isolation.get("donor_fold_counts", ()))))
        or isolation.get("all_nuclei_from_one_donor_in_one_outer_fold") is not True
    ):
        raise SequenceControlReadinessError("donor fold isolation differs")

    pseudobulk_root, pseudobulk = _require_manifest_only(
        root,
        str(isolation.get("training_pseudobulk_artifact_root", "")),
        str(isolation.get("training_pseudobulk_artifacts_sha256", "")),
        "gse296875_training_fold_pseudobulk",
    )
    pseudobulk_inventory = _inventory(pseudobulk)
    pseudobulk_contract_path = pseudobulk_root / "contract.json"
    if _digest(pseudobulk_contract_path) != pseudobulk_inventory.get(
        "contract.json", {}
    ).get("sha256"):
        raise SequenceControlReadinessError("pseudobulk contract bytes differ")
    pseudobulk_contract = _load_json(pseudobulk_contract_path)
    if (
        pseudobulk.get("metadata", {}).get("biological_outer_unit") != "donor"
        or pseudobulk_contract.get("donor_folds") != 5
        or pseudobulk_contract.get("training_folds_per_outer_split") != 3
        or pseudobulk_contract.get("donor_valid_fold_rule") != "test_plus_1_mod_5"
        or pseudobulk_contract.get("held_donor_test_used") is not False
        or pseudobulk_contract.get("held_donor_validation_used") is not False
    ):
        raise SequenceControlReadinessError("training pseudobulk isolation differs")

    split_root, split_manifest = _require_manifest_only(
        root,
        str(isolation.get("sequence_split_artifact_root", "")),
        str(isolation.get("sequence_split_artifacts_sha256", "")),
        "sequence_split_contract",
    )
    split_inventory = _inventory(split_manifest)
    split_contract_path = split_root / "contract.json"
    crossed_path = split_root / "crossed_outer_splits.tsv"
    genomic_path = split_root / "genomic_folds.tsv"
    for member, path in (
        ("contract.json", split_contract_path),
        ("crossed_outer_splits.tsv", crossed_path),
        ("genomic_folds.tsv", genomic_path),
    ):
        if _digest(path) != split_inventory.get(member, {}).get("sha256"):
            raise SequenceControlReadinessError(f"split member differs: {member}")
    split_contract = _load_json(split_contract_path)
    crossed = _read_tsv(crossed_path)
    genomic = _read_tsv(genomic_path)
    diagonal = [
        row
        for row in crossed
        if row["split_id"] == f"donor{row['donor_test_fold']}_genomic{row['genomic_test_fold']}"
        and row["donor_test_fold"] == row["genomic_test_fold"]
    ]
    contigs = [contig for row in genomic for contig in row["contigs"].split(",")]
    primary = {f"chr{index}" for index in range(1, 23)} | {"chrX", "chrY"}
    if (
        len(crossed) != 25
        or len(diagonal) != 5
        or len(genomic) != 5
        or len(contigs) != len(set(contigs))
        or set(contigs) != primary
        or split_contract.get("genomic_outer_unit") != "whole_chromosome_group"
        or split_contract.get("boundary_buffer_required") is not False
        or split_contract.get("boundary_buffer_reason")
        != "whole_contigs_do_not_share_linear_boundaries"
        or split_contract.get("held_atac_used_for_window_selection") is not False
        or split_contract.get("ccre_windows") != 80000
        or isolation.get("largest_control_receptive_field_bp") != 2114
    ):
        raise SequenceControlReadinessError("genomic split or buffer isolation differs")


def _audit_rectangle_and_evaluator(root: Path, contract: Mapping[str, Any]) -> None:
    rectangle = contract.get("five_seed_rectangle", {})
    campaign_path = _resolve_file(root, str(rectangle.get("contract_path", "")))
    if _digest(campaign_path) != rectangle.get("contract_sha256"):
        raise SequenceControlReadinessError("rectangle contract differs")
    campaign = _load_json(campaign_path)
    if (
        tuple(campaign.get("lineages", ())) != LINEAGES
        or tuple(campaign.get("diagonal_outer_folds", ())) != FOLDS
        or tuple(campaign.get("fixed_seeds", ())) != SEEDS
        or campaign.get("prediction_views", {}).get("sequence_cnn_control")
        != ["sequence_cnn_control"]
        or campaign.get("prediction_views", {}).get("sequence_transformer_control")
        != ["sequence_transformer_control"]
        or campaign.get("execution", {}).get("submission_authority")
        != "dispatcher_only_no_manual_gpu_sbatch"
    ):
        raise SequenceControlReadinessError("rectangle campaign semantics differ")
    _assert_false(
        campaign.get("claims", {}),
        (
            "champion_claim_allowed",
            "universal_claim_allowed",
            "cross_model_ranking_allowed",
            "post_hoc_unit_selection_allowed",
        ),
    )

    plan_root, plan = _require_frozen_artifact(
        root,
        str(rectangle.get("admission_artifact_root", "")),
        str(rectangle.get("admission_artifacts_sha256", "")),
        "sequence_task_native_five_seed_rectangle_admission",
    )
    _, evaluator = _require_frozen_artifact(
        root,
        str(rectangle.get("evaluator_readiness_artifact_root", "")),
        str(rectangle.get("evaluator_readiness_artifacts_sha256", "")),
        "sequence_task_native_five_seed_evaluator_readiness",
    )
    plan_meta = plan.get("metadata", {})
    evaluator_meta = evaluator.get("metadata", {})
    if (
        plan_meta.get("complete_reusable_fits") != 62
        or plan_meta.get("expected_fits") != 500
        or plan_meta.get("expected_prediction_views") != 625
        or plan_meta.get("prediction_or_metric_values_read") is not False
        or evaluator_meta.get("gate_status") != "waiting_for_full_rectangle"
        or evaluator_meta.get("full_compatible_rectangle") is not False
        or evaluator_meta.get("outcome_access_authorized") is not False
        or evaluator_meta.get("partial_cross_model_ranking_authorized") is not False
        or evaluator_meta.get("development_cross_model_ranking_authorized") is not False
        or evaluator_meta.get("prediction_or_metric_values_read") is not False
    ):
        raise SequenceControlReadinessError("rectangle or evaluator gate differs")

    rows = _read_tsv(plan_root / "plan/expected_fit_matrix.tsv")
    expected_keys = {
        (lineage, fold, seed)
        for lineage in LINEAGES
        for fold in FOLDS
        for seed in SEEDS
    }
    for model_id in MODEL_IDS:
        model_rows = [row for row in rows if row.get("model_id") == model_id]
        actual_keys = {
            (row["lineage_id"], int(row["outer_fold"]), int(row["seed"]))
            for row in model_rows
        }
        if (
            len(model_rows) != 125
            or actual_keys != expected_keys
            or any(
                row.get("split_id")
                != f"donor{row['outer_fold']}_genomic{row['outer_fold']}"
                or row.get("prediction_views") != model_id
                or row.get("fit_state") not in {"complete_reusable", "missing"}
                for row in model_rows
            )
        ):
            raise SequenceControlReadinessError(f"rectangle geometry differs: {model_id}")
        complete = [row for row in model_rows if row["fit_state"] == "complete_reusable"]
        observed = {
            "complete": len(complete),
            "missing": len(model_rows) - len(complete),
            "lineage": dict(Counter(row["lineage_id"] for row in complete)),
            "fold": dict(Counter(int(row["outer_fold"]) for row in complete)),
            "seed": {
                seed: Counter(int(row["seed"]) for row in complete)[seed]
                for seed in SEEDS
            },
        }
        expected = EXPECTED_COUNTS[model_id]
        recorded = rectangle.get(model_id, {})
        contexts = sorted(
            [
                [row["lineage_id"], int(row["outer_fold"]), int(row["seed"])]
                for row in complete
            ]
        )
        if (
            observed != expected
            or recorded.get("complete_fits_at_freeze") != expected["complete"]
            or recorded.get("missing_fits_at_freeze") != expected["missing"]
            or recorded.get("complete_by_lineage") != expected["lineage"]
            or {int(k): v for k, v in recorded.get("complete_by_fold", {}).items()}
            != expected["fold"]
            or {int(k): v for k, v in recorded.get("complete_by_seed", {}).items()}
            != expected["seed"]
            or contexts != sorted(recorded.get("complete_contexts", []))
            or recorded.get("full_rectangle") is not False
        ):
            raise SequenceControlReadinessError(f"rectangle counts differ: {model_id}")


def _audit_partial_diagnostics(root: Path, contract: Mapping[str, Any]) -> None:
    for model_id in MODEL_IDS:
        record = contract.get("partial_diagnostics", {}).get(model_id, {})
        _, manifest = _require_manifest_only(
            root,
            str(record.get("artifact_root", "")),
            str(record.get("artifacts_sha256", "")),
            str(record.get("artifact_class", "")),
        )
        metadata = manifest.get("metadata", {})
        if (
            metadata.get("dataset_id") != "gse296875"
            or metadata.get("evaluation_role") != "valid"
            or metadata.get("test_outcomes_read") is not False
            or metadata.get("sealed_data_read") is not False
            or metadata.get("cross_model_ranking_calculated") is not False
            or metadata.get("promotion_gate_evaluated") is not False
        ):
            raise SequenceControlReadinessError(f"partial diagnostic differs: {model_id}")
        units = metadata.get("units", metadata.get("evaluation_units"))
        if units != 8:
            raise SequenceControlReadinessError(f"partial diagnostic units differ: {model_id}")
    partial = contract.get("partial_diagnostics", {})
    if (
        partial.get("metric_or_prediction_members_opened_by_this_reconciliation")
        is not False
        or partial.get("partial_results_may_select_or_rank_models") is not False
    ):
        raise SequenceControlReadinessError("partial-result firewall differs")


def _audit_license_queue_and_task_isolation(
    root: Path, contract: Mapping[str, Any]
) -> None:
    license_gate = contract.get("license_and_exposure", {})
    upstream_root, upstream = _require_manifest_only(
        root,
        str(license_gate.get("upstream_acquisition_artifact_root", "")),
        str(license_gate.get("upstream_acquisition_artifacts_sha256", "")),
        "chrombpnet_source_acquisition",
    )
    license_path = _resolve_file(root, str(license_gate.get("upstream_license_member", "")))
    try:
        license_member = str(license_path.relative_to(upstream_root))
    except ValueError as error:
        raise SequenceControlReadinessError(
            "upstream license escapes acquisition artifact"
        ) from error
    upstream_inventory = _inventory(upstream)
    if (
        _digest(license_path) != license_gate.get("upstream_license_sha256")
        or upstream_inventory.get(license_member, {}).get("sha256")
        != license_gate.get("upstream_license_sha256")
        or not license_path.read_text(encoding="utf-8").startswith("MIT License\n")
        or license_gate.get("upstream_code_license") != "MIT"
    ):
        raise SequenceControlReadinessError("upstream license evidence differs")
    repository_root = root.parents[1]
    local_license_files = list(root.glob("LICENSE*")) + list(repository_root.glob("LICENSE*"))
    if (
        local_license_files
        or license_gate.get("project_repository_license")
        != "UNRESOLVED_no_repository_level_license_detected"
        or license_gate.get("source_redistribution_allowed") is not False
        or license_gate.get("derivative_weight_redistribution_allowed") is not False
        or license_gate.get("checkpoint_corpus_contamination")
        != "not_applicable_random_initialization_no_pretraining"
        or license_gate.get("GSE296875_fit_exposure")
        != "continual_seen_for_exact_outer_training_donors_loci_sequences_and_targets"
        or license_gate.get("external_evaluation_complete") is not False
        or license_gate.get("sealed_evaluation_complete") is not False
        or license_gate.get("open_champion_eligible_now") is not False
    ):
        raise SequenceControlReadinessError("project license or exposure gate differs")

    queue = contract.get("central_queue", {})
    _require_frozen_artifact(
        root,
        str(queue.get("execution_correction_artifact_root", "")),
        str(queue.get("execution_correction_artifacts_sha256", "")),
        "sequence_task_native_five_seed_execution_correction",
    )
    queue_root, queue_manifest = _require_frozen_artifact(
        root,
        str(queue.get("queue_admission_artifact_root", "")),
        str(queue.get("queue_admission_artifacts_sha256", "")),
        "sequence_task_native_five_seed_corrected_queue_admission",
    )
    _, live = _require_frozen_artifact(
        root,
        str(queue.get("live_readiness_artifact_root", "")),
        str(queue.get("live_readiness_artifacts_sha256", "")),
        "sequence_task_native_five_seed_live_dispatch_readiness",
    )
    queue_meta = queue_manifest.get("metadata", {})
    live_meta = live.get("metadata", {})
    if (
        queue_meta.get("logical_fits") != 438
        or queue_meta.get("queue_items") != 30
        or queue_meta.get("submission_authority")
        != "dispatcher_only_no_manual_gpu_sbatch"
        or queue_meta.get("manual_gpu_sbatch_allowed") is not False
        or queue_meta.get("qos") != "nslab"
        or live_meta.get("central_dispatcher_running") is not True
        or live_meta.get("live_v2_byte_identical_to_frozen_admission") is not True
        or live_meta.get("submission_authority") != "central_dispatcher_only"
        or live_meta.get("innovation_used") is not False
    ):
        raise SequenceControlReadinessError("central queue authority differs")

    plan_root = _resolve_directory(
        root, str(contract["five_seed_rectangle"]["admission_artifact_root"])
    )
    tasks = _read_tsv(plan_root / "plan/bundle_tasks.tsv")
    for model_id in MODEL_IDS:
        record = queue.get(model_id, {})
        model_tasks = [row for row in tasks if row.get("model_id") == model_id]
        task_bundles = sorted({row["bundle_id"] for row in model_tasks})
        if (
            len(model_tasks) != record.get("missing_fits_registered")
            or task_bundles != record.get("task_bundle_ids")
        ):
            raise SequenceControlReadinessError(f"queue fit census differs: {model_id}")
        dispatcher_bundles = record.get("dispatcher_bundle_ids", [])
        if len(dispatcher_bundles) != len(task_bundles):
            raise SequenceControlReadinessError(f"dispatcher census differs: {model_id}")
        for task_bundle, dispatcher_bundle in zip(task_bundles, dispatcher_bundles, strict=True):
            if int(dispatcher_bundle.rsplit("-", 1)[1]) != int(task_bundle.rsplit("-", 1)[1]) + 100:
                raise SequenceControlReadinessError(f"dispatcher mapping differs: {model_id}")
            frozen_candidates = list(
                (queue_root / "queue_items").glob(f"*-{dispatcher_bundle}-v2.json")
            )
            active_candidates = list(
                (root / "config/campaigns/gpu_bundle_queue").glob(
                    f"*-{dispatcher_bundle}-v2.json"
                )
            )
            if len(frozen_candidates) != 1 or len(active_candidates) != 1:
                raise SequenceControlReadinessError(f"queue item multiplicity differs: {model_id}")
            if active_candidates[0].read_bytes() != frozen_candidates[0].read_bytes():
                raise SequenceControlReadinessError(f"active queue item drifted: {model_id}")
            item = _load_json(active_candidates[0])
            if (
                item.get("enabled") is not True
                or item.get("bundle_id") != dispatcher_bundle
                or item.get("exports", {}).get("RECTANGLE_BUNDLE_ID") != task_bundle
                or item.get("wrapper_path")
                != "slurm/run_sequence_task_native_five_seed_bundle_v2.sbatch"
            ):
                raise SequenceControlReadinessError(f"queue item semantics differ: {model_id}")
        first_path = _resolve_file(root, str(record.get("first_active_queue_item_path", "")))
        if _digest(first_path) != record.get("first_active_queue_item_sha256"):
            raise SequenceControlReadinessError(f"first queue item differs: {model_id}")
    if (
        queue.get("all_missing_control_fits_already_registered") is not True
        or queue.get("new_GPU_queue_registration_needed") is not False
        or queue.get("direct_GPU_submission_authorized") is not False
        or queue.get("manual_GPU_sbatch_allowed") is not False
    ):
        raise SequenceControlReadinessError("queue eligibility disposition differs")

    task = contract.get("task_isolation", {})
    control_path = _resolve_file(root, str(task.get("GSE281364_control_contract_path", "")))
    production_path = _resolve_file(
        root, str(task.get("GSE281364_production_contract_path", ""))
    )
    if (
        _digest(control_path) != task.get("GSE281364_control_contract_sha256")
        or _digest(production_path) != task.get("GSE281364_production_contract_sha256")
    ):
        raise SequenceControlReadinessError("GSE281364 task contract differs")
    gse281364 = _load_toml(control_path)
    gse281364_production = _load_toml(production_path)
    _, preflight = _require_manifest_only(
        root,
        str(task.get("GSE281364_preflight_artifact_root", "")),
        str(task.get("GSE281364_preflight_artifacts_sha256", "")),
        "gse281364_task_native_sequence_control_production_preflight",
    )
    if (
        gse281364.get("dataset_id") != "gse281364"
        or gse281364.get("sequence_contract", {}).get("input_length_bp") != 230
        or tuple(gse281364.get("row_contract", {}).get("fixed_seeds", ()))
        != LEGACY_SEEDS
        or gse281364_production.get("execution", {}).get("fit_count") != 50
        or preflight.get("metadata", {}).get("gpu_jobs_submitted") != 0
        or task.get("GSE281364_fits_may_fill_GSE296875_rectangle") is not False
    ):
        raise SequenceControlReadinessError("task isolation differs")


def audit_readiness(root: Path, contract_path: Path) -> dict[str, Any]:
    """Verify structural authorities without opening predictions, metrics, or outcomes."""

    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    try:
        contract_path.relative_to(root)
    except ValueError as error:
        raise SequenceControlReadinessError("contract escapes benchmark root") from error
    contract = _load_json(contract_path)
    if contract.get("schema_version") != SCHEMA_VERSION or contract.get("status") != STATUS:
        raise SequenceControlReadinessError("readiness schema or status differs")

    _audit_architectures(root, contract)
    _audit_reference_artifacts(root, contract)
    _audit_authority_history(root, contract)
    _audit_split_isolation(root, contract)
    _audit_rectangle_and_evaluator(root, contract)
    _audit_partial_diagnostics(root, contract)
    _audit_license_queue_and_task_isolation(root, contract)

    evaluator = contract.get("evaluator_closure", {})
    _assert_false(
        evaluator,
        (
            "outcome_access_authorized",
            "partial_cross_model_ranking_authorized",
            "development_cross_model_ranking_authorized",
            "external_evaluation_complete",
            "sealed_evaluation_complete",
            "champion_claim_allowed",
            "universal_claim_allowed",
        ),
    )
    firewall = contract.get("firewall", {})
    _assert_false(
        firewall,
        (
            "prediction_values_may_be_read",
            "metric_values_may_be_read",
            "development_outcomes_may_be_read",
            "sealed_features_may_be_read",
            "sealed_labels_or_outcomes_may_be_read",
            "partial_results_may_be_ranked",
            "thresholds_may_be_changed",
            "GPU_job_may_be_submitted_by_this_audit",
        ),
    )

    return {
        "schema_version": "masld-bench-sequence-control-production-readiness-receipt-v1",
        "status": "pass_controls_reconciled_waiting_for_full_rectangles_license_and_external_evaluation",
        "contract_sha256": _digest(contract_path),
        "model_ids": list(MODEL_IDS),
        "controls_are_separate_architectures": True,
        "controls_share_checkpoint": False,
        "cnn_parameter_count": 3696386,
        "transformer_parameter_count": 3683330,
        "parameter_budget_matched": True,
        "biological_donors": 39,
        "donor_fold_counts": [11, 9, 7, 4, 8],
        "genomic_outer_unit": "whole_chromosome_group",
        "linear_boundary_buffer_required": False,
        "whole_contigs_shared_across_genomic_folds": False,
        "cnn_expected_fits": 125,
        "cnn_complete_fits_at_freeze": 12,
        "cnn_missing_fits_at_freeze": 113,
        "transformer_expected_fits": 125,
        "transformer_complete_fits_at_freeze": 11,
        "transformer_missing_fits_at_freeze": 114,
        "full_compatible_rectangle": False,
        "partial_evaluator_access_authorized": False,
        "legacy_transformer_registry_state_superseded_for_development": True,
        "legacy_seed_authority_superseded_for_rectangle": True,
        "upstream_code_license": "MIT",
        "project_repository_license_resolved": False,
        "source_redistribution_authorized": False,
        "derivative_weight_redistribution_authorized": False,
        "open_champion_eligible_now": False,
        "cnn_missing_fits_already_centrally_registered": 113,
        "transformer_missing_fits_already_centrally_registered": 114,
        "cnn_next_dispatcher_bundle_id": "model-training-815",
        "transformer_next_dispatcher_bundle_id": "model-training-823",
        "new_gpu_queue_registration_needed": False,
        "direct_gpu_submission_authorized": False,
        "gse281364_fits_fill_gse296875_rectangle": False,
        "prediction_values_read": False,
        "metric_values_read": False,
        "development_outcomes_read": False,
        "sealed_features_read": False,
        "sealed_labels_or_outcomes_read": False,
        "partial_results_ranked": False,
        "thresholds_changed": False,
        "gpu_job_submitted": False,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--contract", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    receipt = audit_readiness(arguments.project_root, arguments.contract)
    rendered = json.dumps(receipt, indent=2, sort_keys=True) + "\n"
    if arguments.output is None:
        print(rendered, end="")
    else:
        arguments.output.parent.mkdir(parents=True, exist_ok=True)
        arguments.output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
