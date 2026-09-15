#!/usr/bin/env python3
"""Validate and freeze the prospective Regular Corgi FiLM-plus-head requirements."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import tomllib
from typing import Any


SCHEMA_VERSION = "masld-bench-corgi-film-plus-head-contract-v1"
EXPECTED_ARMS = [
    "actual_released_rank_masked",
    "actual_length_adjusted_tpm_rank_masked",
    "training_lineage_mean_released_rank",
    "nearest_training_released_rank",
    "shuffled_valid_released_rank",
]
EXECUTION_FLAGS = (
    "training_execution_authorized",
    "prediction_execution_authorized",
    "evaluation_execution_authorized",
)
CLAIM_FLAGS = (
    "gse296875_masld_diagnosis_claimed",
    "histology_or_disease_endpoint_used",
    "test_or_sealed_features_used",
    "test_or_sealed_outcomes_used",
    "champion_or_external_transfer_claim_allowed",
    "global_census_revision",
    "global_promotion_gate_revision",
    "historical_audit_hashes_modified",
)
BUFFER_SUFFIXES = (".running_mean", ".running_var", ".num_batches_tracked")


class CorgiFilmHeadContractError(ValueError):
    """Raised when the prospective FiLM-plus-head requirements are unsafe or drifts."""


def _digest(path: Path) -> str:
    return sha256(path.read_bytes()).hexdigest()


def _load_toml(path: Path) -> dict[str, Any]:
    value = tomllib.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise CorgiFilmHeadContractError(f"TOML authority is not a table: {path}")
    return value


def _load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise CorgiFilmHeadContractError(f"JSON authority is not an object: {path}")
    return value


def _resolve_member(root: Path, relative: str) -> Path:
    member = (root / relative).resolve(strict=True)
    try:
        member.relative_to(root)
    except ValueError as error:
        raise CorgiFilmHeadContractError("contract member escapes benchmark root") from error
    if not member.is_file() or member.is_symlink():
        raise CorgiFilmHeadContractError(f"contract member is not a regular file: {relative}")
    return member


def _verify_path_hash_table(
    root: Path, table: dict[str, Any], *, skip: set[str] | None = None
) -> dict[str, str]:
    skip = skip or set()
    verified: dict[str, str] = {}
    for key, relative in sorted(table.items()):
        if not key.endswith("_path") or key in skip:
            continue
        hash_key = f"{key[:-5]}_sha256"
        expected = table.get(hash_key)
        if not isinstance(relative, str) or not isinstance(expected, str):
            raise CorgiFilmHeadContractError(f"path/hash pair is incomplete: {key}")
        member = _resolve_member(root, relative)
        observed = _digest(member)
        if observed != expected:
            raise CorgiFilmHeadContractError(f"contract member drifted: {relative}")
        verified[relative] = observed
    return verified


def _assert_false(mapping: dict[str, Any], fields: tuple[str, ...]) -> None:
    for field in fields:
        if mapping.get(field) is not False:
            raise CorgiFilmHeadContractError(f"gate is not fail-closed: {field}")


def _audit_folds(contract: dict[str, Any], prior: dict[str, Any]) -> list[dict[str, Any]]:
    folds = contract.get("crossfit_folds", [])
    prior_folds = prior.get("crossfit_folds", [])
    if len(folds) != 5 or len(prior_folds) != 5:
        raise CorgiFilmHeadContractError("five crossed folds are required")
    shared = (
        "evaluation_fold",
        "evaluation_base_outer_fold",
        "fit_valid_folds",
        "fit_base_outer_folds",
    )
    for expected_fold, (fold, old) in enumerate(zip(folds, prior_folds, strict=True)):
        if any(fold.get(field) != old.get(field) for field in shared):
            raise CorgiFilmHeadContractError(
                f"FiLM rung does not preserve prior crossed fold {expected_fold}"
            )
        if fold.get("evaluation_fold") != expected_fold:
            raise CorgiFilmHeadContractError("evaluation fold order differs")
        fit = set(fold.get("fit_valid_folds", []))
        inner_valid = fold.get("inner_validation_fold")
        inner_train = set(fold.get("inner_training_folds", []))
        if (
            fit != set(range(5)) - {expected_fold}
            or inner_valid not in fit
            or inner_train != fit - {inner_valid}
            or expected_fold in inner_train
            or expected_fold == inner_valid
        ):
            raise CorgiFilmHeadContractError(
                f"inner selection leaks or omits fit folds for outer fold {expected_fold}"
            )
    return folds


def _audit_checkpoint_schema(
    contract: dict[str, Any], schema: dict[str, Any]
) -> dict[str, int]:
    architecture = contract.get("architecture", {})
    source = contract.get("source_artifacts", {})
    if (
        schema.get("checkpoint_sha256") != source.get("checkpoint_sha256")
        or schema.get("checkpoint") != "corgi_model.pt"
        or schema.get("schema_version") != "masld-bench-corgi-checkpoint-schema-v1"
    ):
        raise CorgiFilmHeadContractError("checkpoint schema identity differs")
    state = schema.get("state", [])
    if not isinstance(state, list):
        raise CorgiFilmHeadContractError("checkpoint state inventory differs")
    counts: dict[str, int] = {}
    trainable = 0
    for prefix in ("film_mlp_conv.", "film_mlp_transformer.", "output_head."):
        records = [row for row in state if row.get("key", "").startswith(prefix)]
        counts[f"{prefix}entries"] = len(records)
        counts[f"{prefix}state_numel"] = sum(int(row["numel"]) for row in records)
        if prefix.startswith("film_mlp"):
            trainable += sum(
                int(row["numel"])
                for row in records
                if not str(row["key"]).endswith(BUFFER_SUFFIXES)
            )
    if (
        counts["film_mlp_conv.entries"]
        != architecture.get("checkpoint_film_mlp_conv_state_entries")
        or counts["film_mlp_transformer.entries"]
        != architecture.get("checkpoint_film_mlp_transformer_state_entries")
        or counts["output_head.entries"]
        != architecture.get("checkpoint_output_head_state_entries")
        or counts["film_mlp_conv.state_numel"]
        != architecture.get("checkpoint_film_mlp_conv_state_numel_including_buffers")
        or counts["film_mlp_transformer.state_numel"]
        != architecture.get("checkpoint_film_mlp_transformer_state_numel_including_buffers")
    ):
        raise CorgiFilmHeadContractError("FiLM or output-head checkpoint surface differs")
    head_numel = 1920 + 1
    if (
        architecture.get("expected_replacement_head_parameter_numel") != head_numel
        or architecture.get("expected_trainable_parameter_numel") != trainable + head_numel
        or architecture.get("trainable_parameter_prefixes")
        != ["film_mlp_conv.", "film_mlp_transformer.", "masld_atac_head."]
    ):
        raise CorgiFilmHeadContractError("prospective trainable parameter surface differs")
    counts["film_trainable_parameter_numel"] = trainable
    counts["replacement_head_parameter_numel"] = head_numel
    counts["total_trainable_parameter_numel"] = trainable + head_numel
    return counts


def validate_contract(root: Path, contract_path: Path) -> dict[str, Any]:
    root = root.resolve(strict=True)
    contract_path = contract_path.resolve(strict=True)
    contract = _load_toml(contract_path)
    if contract.get("schema_version") != SCHEMA_VERSION:
        raise CorgiFilmHeadContractError("FiLM-plus-head schema differs")
    if contract.get("status") != "frozen_prospective_preflight_execution_not_yet_authorized":
        raise CorgiFilmHeadContractError("prospective status differs")
    if (
        contract.get("model_id") != "corgi_regular"
        or contract.get("task_id") != "rna_conditioned_atac"
        or contract.get("dataset_id") != "gse296875"
        or contract.get("lineage_id") != "hepatocyte"
        or contract.get("adaptation_rung") != "film_plus_head"
        or contract.get("observed_atac_task") is not True
    ):
        raise CorgiFilmHeadContractError("prospective task identity differs")
    _assert_false(contract, EXECUTION_FLAGS + CLAIM_FLAGS)

    source = contract.get("source_artifacts", {})
    released = contract.get("released_source", {})
    verified = _verify_path_hash_table(root, source)
    verified.update(_verify_path_hash_table(root, released, skip={"source_root"}))
    prior = _load_toml(_resolve_member(root, source["new_head_contract_path"]))
    negative = _load_json(_resolve_member(root, source["negative_disposition_path"]))
    schema = _load_json(_resolve_member(root, source["checkpoint_schema_path"]))
    tile_summary = _load_json(_resolve_member(root, source["tile_subset_summary_path"]))
    folds = _audit_folds(contract, prior)
    counts = _audit_checkpoint_schema(contract, schema)

    if (
        negative.get("status")
        != "NEGATIVE_DEVELOPMENT_SMOKE_NO_INCREMENTAL_CORGI_PROFILE_CONTRIBUTION"
        or negative.get("next_rung", {}).get("name") != "film_plus_head"
        or negative.get("next_rung", {}).get(
            "training_or_prediction_authorized_by_this_disposition"
        )
        is not False
    ):
        raise CorgiFilmHeadContractError("negative prior-rung boundary differs")
    if (
        tile_summary.get("status") != "pass_outcome_free_smoke_subset"
        or tile_summary.get("tiles_per_fold") != 128
        or tile_summary.get("tiles") != 640
        or tile_summary.get("outcomes_or_labels_read") is not False
        or tile_summary.get("model_predictions_read") is not False
    ):
        raise CorgiFilmHeadContractError("outcome-free tile subset differs")

    architecture = contract.get("architecture", {})
    state = contract.get("module_state", {})
    objective = contract.get("objective", {})
    optimization = contract.get("optimization", {})
    sampling = contract.get("sampling", {})
    ablations = contract.get("context_ablations", {})
    evaluation = contract.get("evaluation", {})
    gates = contract.get("required_pre_execution_gates", {})
    terminal = contract.get("terminal_boundaries", {})
    if (
        architecture.get("replacement_head_initialization")
        != "copy_pretrained_output_head_channel_1_weight_and_bias"
        or architecture.get("replacement_head_output_transform") != "softplus"
        or architecture.get("final_conv_trainable") is not False
        or architecture.get("transformer_blocks_trainable") is not False
        or architecture.get("convolution_backbone_trainable") is not False
        or architecture.get("film_layers_have_trainable_parameters") is not False
        or released.get("upstream_finetune_helper_reused") is not False
    ):
        raise CorgiFilmHeadContractError("architecture or source-reuse boundary differs")
    if (
        state.get("whole_model_mode_during_optimization") != "eval"
        or state.get("film_dropout_active") is not False
        or state.get("all_batchnorm_running_statistics_frozen") is not True
        or state.get("film_batchnorm_affine_parameters_trainable") is not True
        or state.get("checkpoint_resume_required_before_campaign") is not True
    ):
        raise CorgiFilmHeadContractError("module-state contract differs")
    if (
        objective.get("training_context_arm") != EXPECTED_ARMS[0]
        or objective.get("context_arm_specific_training_forbidden") is not True
        or objective.get("held_fold_atac_outcomes_during_fit_forbidden") is not True
        or objective.get("histology_or_disease_labels_forbidden") is not True
        or objective.get("sealed_or_test_data_forbidden") is not True
        or objective.get("loss") != "poisson_multinomial"
    ):
        raise CorgiFilmHeadContractError("training objective firewall differs")
    if (
        optimization.get("seed") != 20260825
        or optimization.get("batch_size") != 1
        or optimization.get("gradient_accumulation_steps") != 4
        or optimization.get("refit_from_original_checkpoint") is not True
        or optimization.get("deterministic_algorithms") is not True
    ):
        raise CorgiFilmHeadContractError("optimization contract differs")
    if (
        sampling.get("tile_universe")
        != "frozen_outcome_free_128_tiles_per_genomic_fold"
        or sampling.get("sampling_uses_atac_signal") is not False
        or sampling.get("every_training_donor_seen_each_epoch") is not True
        or sampling.get("held_evaluation_uses_all_128_tiles") is not True
        or sampling.get("boundary_buffer_bp") != 524288
    ):
        raise CorgiFilmHeadContractError("sampling or boundary contract differs")
    if (
        ablations.get("held_prediction_arms") != EXPECTED_ARMS
        or ablations.get("fit_context_arm") != EXPECTED_ARMS[0]
        or ablations.get("same_fitted_weights_applied_to_every_arm") is not True
        or ablations.get("arm_specific_early_stopping_forbidden") is not True
        or ablations.get("arm_specific_calibration_forbidden") is not True
    ):
        raise CorgiFilmHeadContractError("context ablations differ from prior rung")
    if (
        evaluation.get("independent_evaluator_required") is not True
        or evaluation.get("fit_predict_stage_computes_benchmark_metrics") is not False
        or evaluation.get("same_unit_and_formula_as_new_assay_head_rung") is not True
        or evaluation.get("numeric_tolerance_cannot_be_interpreted_as_gain") is not True
        or evaluation.get("promotion_gate_evaluated") is not False
        or evaluation.get("development_diagnostic_only") is not True
        or evaluation.get("full_family_ranking_allowed") is not False
    ):
        raise CorgiFilmHeadContractError("independent evaluation firewall differs")
    gate_names = (
        "adapter_implementation_present",
        "parameter_mask_unit_test_passed",
        "single_step_l40s_test_passed",
        "checkpoint_resume_test_passed",
        "training_fixture_leakage_test_passed",
        "held_prediction_export_test_passed",
        "separate_evaluator_test_passed",
        "gpu_dispatch_bundle_registered",
    )
    _assert_false(gates, gate_names)
    if gates.get("all_gates_required_before_training") is not True:
        raise CorgiFilmHeadContractError("pre-execution gates are not conjunctive")
    if (
        terminal.get("native_numeric_parity_established") is not False
        or terminal.get("open_champion_eligible") is not False
        or terminal.get("external_validation_present") is not False
        or terminal.get("negative_new_head_result_overridden") is not False
    ):
        raise CorgiFilmHeadContractError("terminal boundary differs")

    return {
        "schema_version": "masld-bench-corgi-film-plus-head-contract-receipt-v1",
        "status": "pass_prospective_contract_execution_blocked",
        "contract_sha256": _digest(contract_path),
        "verified_source_hashes": verified,
        "crossfit_folds": len(folds),
        "context_arms": EXPECTED_ARMS,
        "trainable_parameter_numel": counts["total_trainable_parameter_numel"],
        "film_trainable_parameter_numel": counts["film_trainable_parameter_numel"],
        "replacement_head_parameter_numel": counts[
            "replacement_head_parameter_numel"
        ],
        "development_outcomes_read": False,
        "histology_or_disease_labels_read": False,
        "test_or_sealed_features_or_outcomes_read": False,
        "training_execution_authorized": False,
        "prediction_execution_authorized": False,
        "evaluation_execution_authorized": False,
        "global_census_modified": False,
        "global_promotion_gate_modified": False,
    }


def freeze_contract(root: Path, contract_path: Path, output: Path) -> dict[str, Any]:
    receipt = validate_contract(root, contract_path)
    if output.exists():
        raise CorgiFilmHeadContractError(f"refusing to overwrite output: {output}")
    output.mkdir(mode=0o750)
    (output / "contract_receipt.json").write_text(
        json.dumps(receipt, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return receipt


def main() -> None:
    default_root = Path(__file__).resolve().parents[1]
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=default_root)
    parser.add_argument(
        "--contract",
        type=Path,
        default=default_root
        / "config/evaluation/corgi_regular_film_plus_head_smoke.toml",
    )
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()
    if arguments.output is None:
        receipt = validate_contract(arguments.project_root, arguments.contract)
    else:
        receipt = freeze_contract(
            arguments.project_root, arguments.contract, arguments.output
        )
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
