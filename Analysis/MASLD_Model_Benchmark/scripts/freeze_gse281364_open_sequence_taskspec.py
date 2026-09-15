#!/usr/bin/env python3
"""Validate and materialize the prespecified GSE281364 head-campaign TaskSpec."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
import tomllib
from typing import Any, Mapping

from masld_bench.artifacts import verify_frozen_tree, write_json_exclusive


SCHEMA = "masld-bench-gse281364-open-sequence-head-campaign-v1"
STATUS = "prespecified_exposed_development_head_campaign"
SEEDS = (1103, 2909, 4721, 6673, 8111)
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
MODELS = ("caduceus", "dnabert2", "hyenadna")
HEADS = ("delta_ridge", "full_ridge", "two_layer_gelu")
BASELINES = ("zero", "training_mean", "allele_identity_ridge")
MISSING_BASELINES = (
    "deltaSVM",
    "gkm-SVM",
    "sequence_CNN",
    "sequence_transformer",
    "MPRALegNet_signed_MPRA_rekeyed_to_1033_elements_239_blocks",
)
PREDICTION_FIELDS = (
    "seed",
    "row_hash",
    "unit_hash",
    "block_hash",
    "stratum",
    "outer_fold",
    "study_id",
    "assay_context_id",
    "element_id",
    "source_locus_group_id",
    "long_range_block_id",
    "model_id",
    "head_id",
    "observed",
    "prediction",
    "experimental_replicates",
    "biological_donors",
    "outcome_role",
)
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class TaskSpecError(RuntimeError):
    """Raised when the campaign is not exactly prespecified or outcome-safe."""


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_config(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise TaskSpecError("campaign config is not a regular file")
    try:
        config = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise TaskSpecError(f"invalid campaign config: {error}") from error
    if not isinstance(config, dict):
        raise TaskSpecError("campaign config must be a table")
    return config


def safe_path(root: Path, relative_value: object, *, label: str) -> Path:
    if not isinstance(relative_value, str) or not relative_value:
        raise TaskSpecError(f"{label} path is missing")
    relative = Path(relative_value)
    if relative.is_absolute() or ".." in relative.parts:
        raise TaskSpecError(f"unsafe {label} path")
    try:
        path = (root / relative).resolve(strict=True)
        path.relative_to(root)
    except (OSError, ValueError) as error:
        raise TaskSpecError(f"{label} path is missing or escapes root") from error
    return path


def verify_tree(
    root: Path, relative_value: object, expected_sha256: object, *, label: str
) -> Path:
    if not isinstance(expected_sha256, str) or not SHA256.fullmatch(expected_sha256):
        raise TaskSpecError(f"{label} hash differs")
    path = safe_path(root, relative_value, label=label)
    if file_sha256(path / "ARTIFACTS.json") != expected_sha256:
        raise TaskSpecError(f"{label} ARTIFACTS hash differs")
    verify_frozen_tree(path)
    return path


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != STATUS
        or config.get("dataset_id") != "gse281364"
        or config.get("task_id") != "variant_to_regulation"
        or config.get("endpoint_id") != "mpra_allelic_direction"
        or config.get("evaluator_id")
        != "gse281364_long_range_open_sequence_development_v1"
        or config.get("outcome_access_authorized") is not True
        or config.get("prediction_generation_authorized") is not True
        or config.get("model_head_fit_authorized") is not True
        or config.get("base_checkpoint_tuning_authorized") is not False
        or config.get("stack_fit_authorized") is not False
        or config.get("residual_correlation_authorized") is not False
        or config.get("conditional_model_build_authorized") is not False
        or config.get("sealed_asset_access_authorized") is not False
    ):
        raise TaskSpecError("campaign authorization boundary differs")
    task = config.get("task_spec", {})
    if (
        task.get("claim_role") != "exposed_development_MPRA_proxy_only"
        or task.get("external_evaluation") is not False
        or task.get("champion_claim_authorized") is not False
        or task.get("clinical_claim_authorized") is not False
        or task.get("biological_donors") != 0
        or task.get("measurement_replicates_per_context") != 4
        or task.get("measurement_replicates_are_biological_replicates") is not False
        or task.get("allele_sign") != "ALT_minus_REF"
        or tuple(task.get("contexts", ())) != CONTEXTS
        or task.get("selected_elements") != 1033
        or task.get("rows_per_seed") != 2066
        or tuple(task.get("fixed_seeds", ())) != SEEDS
        or task.get("seeded_rows_per_candidate") != 10330
        or task.get("outer_folds") != 5
        or task.get("long_range_blocks") != 239
        or task.get("largest_receptive_field_buffer_bp") != 524288
        or task.get("seeds_are_independent_biological_replicates") is not False
        or task.get("minimum_positive_gain_seeds_for_partial_one_se_screen") != 4
    ):
        raise TaskSpecError("TaskSpec replication, row, or split contract differs")
    completion = config.get("completion_firewall", {})
    if (
        completion.get("campaign_scope") != "partial_common_head_campaign"
        or completion.get("mandatory_baselines_complete") is not False
        or tuple(completion.get("mandatory_baselines_missing", ()))
        != MISSING_BASELINES
        or completion.get("strongest_available_control_only") is not True
        or completion.get("allele_identity_ridge_role")
        != "16_substitution_class_control_not_strongest_task_native_baseline"
        or completion.get("shortlist_blocked") is not True
        or completion.get("finalist_claim_blocked") is not True
        or completion.get("complementarity_blocked") is not True
        or completion.get("conditional_trigger_blocked") is not True
    ):
        raise TaskSpecError("partial-campaign completion firewall differs")
    uncertainty = config.get("uncertainty", {})
    if (
        uncertainty.get("method") != "paired_long_range_block_bootstrap"
        or uncertainty.get("resamples") != 10000
        or uncertainty.get("seed") != 20260825
        or uncertainty.get("resampling_unit") != "long_range_block_id"
        or uncertainty.get("contexts_resampled_together") is not True
        or uncertainty.get("source_locus_groups_within_sampled_block_kept_together")
        is not True
        or uncertainty.get("confidence_level") != 0.95
    ):
        raise TaskSpecError("uncertainty contract differs")
    multiplicity = config.get("multiplicity", {})
    if (
        multiplicity.get("confirmatory_family") != "none_exposed_development_screen"
        or multiplicity.get("confirmatory_p_values_reported") is not False
        or multiplicity.get("familywise_claim_authorized") is not False
        or multiplicity.get("development_pairwise_tests") is not False
    ):
        raise TaskSpecError("multiplicity boundary differs")
    projection = config.get("projection_contract", {})
    if (
        projection.get("outer_projection_fit") != "training_long_range_blocks_only"
        or projection.get("projection_width") != 256
        or projection.get("outer_head_feature_width") != 1024
        or tuple(projection.get("feature_blocks", ()))
        != ("REF", "ALT", "ALT_minus_REF", "absolute_ALT_minus_REF")
        or projection.get("delta_feature_slice_start0") != 512
        or projection.get("delta_feature_slice_end0") != 768
        or projection.get("held_fold_outcomes_read_during_fit") is not False
        or projection.get("held_fold_features_used_for_fit_or_tuning") is not False
    ):
        raise TaskSpecError("nested projection boundary differs")
    ridge = config.get("ridge_recipe", {})
    if (
        ridge.get("solver") != "sklearn_Ridge_lsqr"
        or tuple(ridge.get("alpha_grid", ()))
        != (0.001, 0.01, 0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)
        or ridge.get("selection_metric") != "inner_grouped_RMSE"
    ):
        raise TaskSpecError("ridge recipe differs")
    mlp = config.get("two_layer_recipe", {})
    if (
        mlp.get("included") is not True
        or mlp.get("input_width") != 1024
        or mlp.get("hidden_width") != 32
        or mlp.get("activation") != "GELU"
        or mlp.get("dropout") != 0.0
        or mlp.get("batch_normalization") is not False
        or mlp.get("optimizer") != "AdamW"
        or mlp.get("learning_rate") != 0.0005
        or mlp.get("weight_decay") != 0.01
        or mlp.get("gradient_clip_norm") != 1.0
        or mlp.get("maximum_epochs") != 100
        or mlp.get("patience") != 12
        or mlp.get("minimum_delta") != 0.000001
        or mlp.get("epoch_selection_seed") != 1103
        or mlp.get("selected_epoch_reused_across_final_seeds") is not True
        or mlp.get("deterministic_algorithms") is not True
        or mlp.get("cpu_only") is not True
    ):
        raise TaskSpecError("small-n two-layer recipe differs")
    prediction = config.get("prediction_schema", {})
    if (
        tuple(prediction.get("fields", ())) != PREDICTION_FIELDS
        or prediction.get("missing_allowed") is not False
        or prediction.get("missing_as_zero_allowed") is not False
        or prediction.get("same_rows_all_candidates") is not True
        or prediction.get("same_rows_all_seeds") is not True
        or prediction.get("same_outer_fold_all_candidates") is not True
    ):
        raise TaskSpecError("standardized prediction schema differs")
    candidates = config.get("candidate", ())
    expected = {
        "available_simple_controls": BASELINES,
        **{model: HEADS for model in MODELS},
    }
    if len(candidates) != 4 or {row.get("model_id") for row in candidates} != set(expected):
        raise TaskSpecError("candidate roster differs")
    for row in candidates:
        if tuple(row.get("head_ids", ())) != expected[row["model_id"]]:
            raise TaskSpecError(f"head roster differs for {row['model_id']}")
        if row.get("open_release_eligible") is not True:
            raise TaskSpecError("restricted model entered open head campaign")


def validate_authorities(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    authorities: dict[str, Path] = {}
    for authority_id in ("row_authority", "outcome_authority", "projection_authority"):
        authority = config[authority_id]
        tree = verify_tree(
            root,
            authority["tree_path"],
            authority["artifacts_sha256"],
            label=authority_id,
        )
        authorities[authority_id] = tree
    if file_sha256(authorities["row_authority"] / config["row_authority"]["row_member"]) != config["row_authority"]["row_sha256"]:
        raise TaskSpecError("row-universe member hash differs")
    if file_sha256(authorities["outcome_authority"] / config["outcome_authority"]["outcome_member"]) != config["outcome_authority"]["outcome_sha256"]:
        raise TaskSpecError("development outcome member hash differs")
    if file_sha256(authorities["projection_authority"] / config["projection_authority"]["projection_member"]) != config["projection_authority"]["projection_receipt_sha256"]:
        raise TaskSpecError("projection receipt hash differs")
    raw_models = set()
    for index, raw in enumerate(config.get("raw_authority", ())):
        raw_tree = verify_tree(
            root,
            raw["tree_path"],
            raw["artifacts_sha256"],
            label=f"raw_authority_{index}",
        )
        fixture_tree = verify_tree(
            root,
            raw["fixture_tree_path"],
            raw["fixture_artifacts_sha256"],
            label=f"raw_fixture_{index}",
        )
        authorities[f"raw_authority_{index}"] = raw_tree
        authorities[f"raw_fixture_{index}"] = fixture_tree
        raw_models.update(raw["model_ids"])
    if raw_models != set(MODELS):
        raise TaskSpecError("raw model authority roster differs")
    for authority_id, authority in config.get("file_authorities", {}).items():
        path = safe_path(root, authority["path"], label=authority_id)
        if path.is_symlink() or not path.is_file() or file_sha256(path) != authority["sha256"]:
            raise TaskSpecError(f"file authority differs: {authority_id}")
        authorities[authority_id] = path
    return authorities


def freeze_spec(*, root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise TaskSpecError("TaskSpec output exists")
    config = load_config(config_path)
    validate_config(config)
    validate_authorities(root, config)
    output.mkdir(parents=True, mode=0o750)
    task_spec = {
        "schema_version": SCHEMA,
        "status": STATUS,
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "endpoint_id": "mpra_allelic_direction",
        "evaluator_id": config["evaluator_id"],
        "task_spec": config["task_spec"],
        "completion_firewall": config["completion_firewall"],
        "uncertainty": config["uncertainty"],
        "multiplicity": config["multiplicity"],
        "projection_contract": config["projection_contract"],
        "ridge_recipe": config["ridge_recipe"],
        "two_layer_recipe": config["two_layer_recipe"],
        "prediction_schema": config["prediction_schema"],
        "candidates": config["candidate"],
        "authority_bindings": {
            "row": config["row_authority"]["artifacts_sha256"],
            "outcomes": config["outcome_authority"]["artifacts_sha256"],
            "projection": config["projection_authority"]["artifacts_sha256"],
        },
    }
    write_json_exclusive(output / "task_spec.json", task_spec)
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_prespecified_taskspec",
        "candidate_count": 12,
        "model_count": 3,
        "available_simple_control_count": 3,
        "fixed_seeds": list(SEEDS),
        "outer_folds": 5,
        "long_range_blocks": 239,
        "rows_per_candidate": 10330,
        "biological_donors": 0,
        "measurement_replicates_per_context": 4,
        "measurement_replicates_are_biological_replicates": False,
        "bootstrap_resamples": 10000,
        "confirmatory_family": "none_exposed_development_screen",
        "campaign_scope": "partial_common_head_campaign",
        "mandatory_baselines_complete": False,
        "mandatory_baselines_missing": list(MISSING_BASELINES),
        "shortlist_blocked": True,
        "finalist_claim_blocked": True,
        "complementarity_blocked": True,
        "conditional_trigger_blocked": True,
        "outcome_values_read": False,
        "model_fit": False,
        "predictions_generated": False,
        "stack_fit": False,
        "residual_correlation_calculated": False,
        "conditional_model_built_or_fit": False,
        "sealed_assets_read": False,
        "config_sha256": file_sha256(config_path),
    }
    write_json_exclusive(output / "receipt.json", receipt)
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    freeze_spec(
        root=arguments.root.resolve(strict=True),
        config_path=arguments.config.resolve(strict=True),
        output=arguments.output,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
