#!/usr/bin/env python3
"""Reconcile four DNA language encoders on one development MPRA head requirements."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import tomllib
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import (
    ArtifactError,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
)


SCHEMA = "masld-bench-dna-lm-family-reconciliation-v1"
MODELS = (
    "dnabert2",
    "nucleotide_transformer",
    "hyenadna",
    "caduceus",
)
SHA256 = re.compile(r"^[0-9a-f]{64}$")
COMMON_MODELS = {"dnabert2", "nucleotide_transformer", "hyenadna"}
EXPECTED_HEADS = ("linear_ridge", "two_layer_gelu")
EXPECTED_BASELINES = ("zero", "training_mean", "allele_identity_ridge")
EXPECTED_CONTEXTS = ("HepG2_control", "HepG2_PAOA")
EXPECTED_FILE_AUTHORITIES = {
    "registry",
    "common_task_contract",
    "caduceus_task_contract",
    "common_fixture_contract",
    "common_evaluator",
    "caduceus_evaluator",
    *{
        f"{model}_{suffix}"
        for model in MODELS
        for suffix in ("checkpoints", "exposure", "crosswalk")
    },
}
EXPECTED_TREE_AUTHORITIES = {
    "source_admission",
    "project_reference",
    "common_fixture",
    "caduceus_fixture",
    "token_coordinate_contract",
    "dnabert2_runtime_probe",
    "nucleotide_transformer_runtime_probe",
    "hyenadna_runtime_probe",
    "caduceus_common_runtime_probe",
    "caduceus_131k_runtime_probe",
    "common_embeddings",
    "caduceus_embeddings",
    "development_outcomes",
    "common_evaluation",
    "caduceus_evaluation",
}
EXPECTED_CLAIM_BOUNDARIES = {
    "sequence_embeds_are_biologically_interpretable_without_head",
    "native_likelihood_is_cell_type_effect_direction",
    "supervised_eQTL_or_ieQTL_head_fit",
    "cell_state_context_supported",
    "external_transfer_claim_supported",
    "universal_model_claim_supported",
    "clinical_claim_supported",
    "sealed_assets_read_by_reconciliation",
}
PROMOTION_FIELDS = (
    "model_id",
    "selected_head",
    "mean_two_context_spearman",
    "allele_baseline_mean_spearman",
    "spearman_gain",
    "positive_spearman_both_contexts",
    "rmse_better_than_training_mean_both_contexts",
    "three_seed_screening_recommended",
    "restricted_comparator",
    "open_champion_eligible_after_task_gates",
)


class DnaLmReconciliationError(RuntimeError):
    """Raised when a frozen family authority or comparison differs."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise DnaLmReconciliationError(f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise DnaLmReconciliationError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise DnaLmReconciliationError(f"{label} must be a JSON object")
    return value


def safe_project_path(root: Path, relative: str) -> Path:
    value = Path(relative)
    if value.is_absolute() or not value.parts or ".." in value.parts:
        raise DnaLmReconciliationError(f"unsafe project path: {relative}")
    path = reject_symlink_components(root / value, label="DNA LM authority")
    try:
        path.resolve(strict=True).relative_to(root)
    except (OSError, ValueError) as error:
        raise DnaLmReconciliationError(
            f"DNA LM authority is missing or escapes root: {relative}"
        ) from error
    return path


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status")
        != "pass_reconciled_development_smoke_no_model_promoted"
        or tuple(config.get("model_order", ())) != MODELS
        or set(config.get("models", {})) != set(MODELS)
        or config.get("family_id") != "dna_language"
    ):
        raise DnaLmReconciliationError("family reconciliation identity differs")
    task = config.get("task_contract", {})
    head = task.get("head_contract", {})
    if (
        task.get("dataset_id") != "gse281364"
        or task.get("biological_donors") != 0
        or task.get("experimental_replicates_per_context") != 4
        or task.get("replicates_used_as_independent_donors") is not False
        or task.get("outcome_role") != "exposed_development_MPRA_only"
        or tuple(task.get("contexts", ())) != EXPECTED_CONTEXTS
        or tuple(head.get("heads", ())) != EXPECTED_HEADS
        or head.get("input_width") != 1024
        or head.get("projection_fit_unit")
        != "outer_training_locus_sequence_groups_only"
        or task.get("external_evaluation") is not False
        or task.get("champion_eligible") is not False
    ):
        raise DnaLmReconciliationError("task or head contract differs")
    models = config["models"]
    if (
        models["nucleotide_transformer"]["restricted_comparator"] is not True
        or models["nucleotide_transformer"][
            "open_champion_eligible_after_task_and_external_gates"
        ]
        is not False
        or any(
            models[model]["restricted_comparator"]
            for model in MODELS
            if model != "nucleotide_transformer"
        )
        or any(
            models[model]["task_input_bp"] != 4096 for model in COMMON_MODELS
        )
        or models["caduceus"]["task_input_bp"] != 131072
    ):
        raise DnaLmReconciliationError("license or family-native input boundary differs")
    boundaries = config.get("claim_boundaries", {})
    if set(boundaries) != EXPECTED_CLAIM_BOUNDARIES or any(boundaries.values()):
        raise DnaLmReconciliationError("claim boundary differs")


def validate_authorities(
    root: Path, config: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    results: dict[str, dict[str, Any]] = {}
    files = config.get("frozen_file_authorities")
    trees = config.get("frozen_tree_authorities")
    if not isinstance(files, dict) or set(files) != EXPECTED_FILE_AUTHORITIES:
        raise DnaLmReconciliationError("frozen file census differs")
    if not isinstance(trees, dict) or set(trees) != EXPECTED_TREE_AUTHORITIES:
        raise DnaLmReconciliationError("frozen tree census differs")
    for name, binding in files.items():
        path = safe_project_path(root, str(binding.get("path", "")))
        expected = str(binding.get("sha256", ""))
        if not SHA256.fullmatch(expected) or sha256_file(path) != expected:
            raise DnaLmReconciliationError(f"frozen file changed: {name}")
        results[name] = {
            "kind": "file",
            "path": str(path.relative_to(root)),
            "sha256": expected,
        }
    for name, binding in trees.items():
        path = safe_project_path(root, str(binding.get("path", "")))
        expected = str(binding.get("artifacts_sha256", ""))
        if not SHA256.fullmatch(expected):
            raise DnaLmReconciliationError(f"frozen tree hash differs: {name}")
        try:
            verify_frozen_tree(path)
        except ArtifactError as error:
            raise DnaLmReconciliationError(
                f"frozen tree verification failed: {name}: {error}"
            ) from error
        if sha256_file(path / "ARTIFACTS.json") != expected:
            raise DnaLmReconciliationError(f"frozen tree changed: {name}")
        results[name] = {
            "kind": "tree",
            "path": str(path.relative_to(root)),
            "artifacts_sha256": expected,
        }
    return results


def registry_models(path: Path) -> dict[str, Mapping[str, Any]]:
    raw = tomllib.loads(path.read_text(encoding="utf-8"))
    rows = {
        row["model_id"]: row
        for row in raw.get("models", [])
        if row.get("model_id") in MODELS
    }
    if set(rows) != set(MODELS):
        raise DnaLmReconciliationError("DNA LM registry model census differs")
    return rows


def validate_model_authorities(
    root: Path, config: Mapping[str, Any]
) -> list[dict[str, Any]]:
    file_bindings = config["frozen_file_authorities"]
    tree_bindings = config["frozen_tree_authorities"]
    registry = registry_models(root / file_bindings["registry"]["path"])
    rows: list[dict[str, Any]] = []
    for model in MODELS:
        expected = config["models"][model]
        checkpoint = load_json(
            root / file_bindings[f"{model}_checkpoints"]["path"],
            label=f"{model} checkpoint audit",
        )
        exposure = load_json(
            root / file_bindings[f"{model}_exposure"]["path"],
            label=f"{model} exposure audit",
        )
        crosswalk = load_json(
            root / file_bindings[f"{model}_crosswalk"]["path"],
            label=f"{model} development crosswalk",
        )
        checkpoint_record = checkpoint.get("checkpoint", {})
        exposure_finding = exposure.get("checkpoint_finding", {})
        effect_finding = crosswalk.get("findings", {}).get(
            "GSE289173_eQTL_ieQTL_effects", {}
        )
        registry_row = registry[model]
        if (
            checkpoint_record.get("repository") != expected["checkpoint_repository"]
            or checkpoint_record.get("revision") != expected["checkpoint_revision"]
            or checkpoint_record.get("sha256")
            != expected["source_checkpoint_sha256"]
            or checkpoint.get("weight_license")
            not in {
                expected["weight_license"],
                expected["weight_license"].removesuffix(
                    "_checkpoint_repository_LICENSE"
                ),
            }
            or exposure_finding.get("exposure_state")
            != expected["exposure_state"]
            or effect_finding.get("exposure_state") != "target_label_unexposed"
            or registry_row.get("checkpoint_sha256")
            != expected["source_checkpoint_sha256"]
            or bool(registry_row.get("status") == "restricted_comparator")
            != expected["restricted_comparator"]
        ):
            raise DnaLmReconciliationError(f"{model} source authority differs")
        runtime_root = root / tree_bindings[expected["runtime_probe_authority"]]["path"]
        runtime = load_json(
            runtime_root / expected["runtime_receipt"],
            label=f"{model} runtime receipt",
        )
        parameter_field = (
            "learned_parameter_count" if model == "caduceus" else "model_parameter_count"
        )
        eligible = runtime.get(
            "champion_eligible_after_task_and_external_evaluation_gates",
            runtime.get(
                "open_champion_eligible_after_task_gates",
                runtime.get(
                    "open_champion_eligible",
                    expected[
                        "open_champion_eligible_after_task_and_external_gates"
                    ],
                ),
            ),
        )
        expected_feature_width = expected.get(
            "aligned_feature_width", expected["raw_feature_width"]
        )
        if (
            runtime.get("status") != "pass"
            or runtime.get("checkpoint_sha256")
            != expected["executable_checkpoint_sha256"]
            or runtime.get(parameter_field) != expected["learned_parameters"]
            or runtime.get("feature_width") != expected_feature_width
            or runtime.get("observed_outcomes_loaded") is not False
            or runtime.get("sealed_outcomes_loaded") is not False
            or runtime.get("head_fit") is not False
            or bool(eligible)
            != expected["open_champion_eligible_after_task_and_external_gates"]
        ):
            raise DnaLmReconciliationError(f"{model} runtime receipt differs")
        rows.append(
            {
                "model_id": model,
                "checkpoint_repository": expected["checkpoint_repository"],
                "checkpoint_revision": expected["checkpoint_revision"],
                "source_checkpoint_sha256": expected["source_checkpoint_sha256"],
                "executable_checkpoint_sha256": expected[
                    "executable_checkpoint_sha256"
                ],
                "weight_license": expected["weight_license"],
                "restricted_comparator": expected["restricted_comparator"],
                "open_champion_eligible_after_task_and_external_gates": expected[
                    "open_champion_eligible_after_task_and_external_gates"
                ],
                "exposure_state": expected["exposure_state"],
                "reference_sequence_exposure": expected[
                    "reference_sequence_exposure"
                ],
                "task_input_bp": expected["task_input_bp"],
                "tokenization": expected["tokenization"],
                "aligned_feature_width": expected_feature_width,
                "learned_parameters": expected["learned_parameters"],
            }
        )
    return rows


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def indexed_rows(
    rows: Sequence[Mapping[str, str]], fields: Sequence[str]
) -> dict[tuple[str, ...], dict[str, str]]:
    result: dict[tuple[str, ...], dict[str, str]] = {}
    for row in rows:
        key = tuple(row[field] for field in fields)
        if key in result:
            raise DnaLmReconciliationError(f"duplicate TSV row: {key}")
        result[key] = dict(row)
    return result


def validate_shared_evaluator_code(root: Path, config: Mapping[str, Any]) -> None:
    paths = config["frozen_file_authorities"]
    source = (root / paths["caduceus_evaluator"]["path"]).read_text(encoding="utf-8")
    required = (
        "from scripts.evaluate_gse281364_dna_lm_common_lane import (",
        "fit_mlp_outer,",
        "fit_ridge_outer,",
        "load_outcomes,",
        "metric_values,",
        "bootstrap_metrics,",
        "CONTEXTS,",
        "SEED,",
    )
    if any(token not in source for token in required):
        raise DnaLmReconciliationError("Caduceus evaluator does not reuse common heads")


def validate_embedding_contracts(
    root: Path, config: Mapping[str, Any]
) -> dict[str, Any]:
    files = config["frozen_file_authorities"]
    trees = config["frozen_tree_authorities"]
    common_task = load_json(
        root / files["common_task_contract"]["path"],
        label="common DNA LM task contract",
    )
    caduceus_task = load_json(
        root / files["caduceus_task_contract"]["path"],
        label="Caduceus task contract",
    )
    common_fixture_contract = load_json(
        root / files["common_fixture_contract"]["path"],
        label="common fixture contract",
    )
    if (
        set(common_task.get("models", {})) != COMMON_MODELS
        or common_task.get("fixture", {}).get("input_length_bp") != 4096
        or set(common_fixture_contract.get("common_lane", {}).get("models", ()))
        != COMMON_MODELS
        or common_fixture_contract.get("common_lane", {}).get("input_length_bp")
        != 4096
        or caduceus_task.get("model_id") != "caduceus"
        or caduceus_task.get("fixture", {}).get("input_length_bp") != 131072
        or common_task.get("fixture", {}).get("elements") != 1033
        or caduceus_task.get("fixture", {}).get("elements") != 1033
    ):
        raise DnaLmReconciliationError("family-native task fixture contract differs")
    common = root / trees["common_embeddings"]["path"]
    caduceus = root / trees["caduceus_embeddings"]["path"]
    common_projection = load_json(
        common / "projected/receipt.json", label="common projection receipt"
    )
    caduceus_projection = load_json(
        caduceus / "projected/receipt.json", label="Caduceus projection receipt"
    )
    shared_fields = (
        "status",
        "dataset_id",
        "elements",
        "outer_locus_sequence_groups",
        "outer_folds",
        "projection_width",
        "head_feature_blocks",
        "head_input_width",
        "linear_head",
        "two_layer_head",
        "projection_fit_on_held_out_fold",
        "reporter_counts_read",
        "reporter_outcomes_read",
        "sealed_labels_read",
        "downstream_head_fit",
    )
    for field in shared_fields:
        if common_projection.get(field) != caduceus_projection.get(field):
            raise DnaLmReconciliationError(
                f"projected feature contract differs between lanes: {field}"
            )
    if (
        set(common_projection.get("models", {})) != COMMON_MODELS
        or set(caduceus_projection.get("models", {})) != {"caduceus"}
        or common_projection.get("projection_width") != 256
        or common_projection.get("head_input_width") != 1024
        or common_projection.get("projection_fit_on_held_out_fold") is not False
        or common_projection.get("reporter_outcomes_read") is not False
        or common_projection.get("sealed_labels_read") is not False
        or common_projection.get("downstream_head_fit") is not False
    ):
        raise DnaLmReconciliationError("projected feature firewall differs")

    common_config = config["models"]
    for model in MODELS:
        embedding_root = caduceus if model == "caduceus" else common
        raw = load_json(
            embedding_root / f"raw/{model}/receipt.json",
            label=f"{model} task embedding receipt",
        )
        expected_width = common_config[model].get(
            "aligned_feature_width", common_config[model]["raw_feature_width"]
        )
        eligibility = raw.get(
            "open_champion_eligible_after_task_and_external_evaluation_gates",
            raw.get("open_champion_eligible_after_task_gates"),
        )
        if (
            raw.get("status") != "pass_outcome_blind_embedding_extraction"
            or raw.get("model_id") != model
            or raw.get("checkpoint_sha256")
            != common_config[model]["executable_checkpoint_sha256"]
            or raw.get("elements") != 1033
            or raw.get("outer_locus_sequence_groups") != 1033
            or raw.get("outer_folds") != 5
            or raw.get("hidden_width") != expected_width
            or raw.get("restricted_comparator")
            != common_config[model]["restricted_comparator"]
            or eligibility
            != common_config[model][
                "open_champion_eligible_after_task_and_external_gates"
            ]
            or raw.get("reporter_counts_read") is not False
            or raw.get("reporter_outcomes_read") is not False
            or raw.get("sealed_labels_read") is not False
            or raw.get("downstream_head_fit") is not False
        ):
            raise DnaLmReconciliationError(f"{model} task embedding receipt differs")
        if model == "caduceus" and raw.get("input_length_bp") != 131072:
            raise DnaLmReconciliationError("Caduceus native input length differs")
        projected_model = (
            caduceus_projection if model == "caduceus" else common_projection
        )["models"][model]
        folds = projected_model.get("folds", ())
        if (
            len(folds) != 5
            or {fold.get("held_out_fold") for fold in folds} != set(range(5))
            or any(fold.get("projection_width") != 256 for fold in folds)
            or any(fold.get("head_input_width") != 1024 for fold in folds)
        ):
            raise DnaLmReconciliationError(f"{model} projected fold contract differs")
    return {
        "outcome_blind_embedding_models": list(MODELS),
        "projection_width": 256,
        "head_input_width": 1024,
        "feature_blocks": [
            "REF",
            "ALT",
            "ALT_minus_REF",
            "absolute_ALT_minus_REF",
        ],
        "projection_fit_on_held_out_fold": False,
        "reporter_outcomes_read_during_embedding_or_projection": False,
        "sealed_labels_read": False,
        "downstream_head_fit_during_embedding_or_projection": False,
    }


def validate_evaluations(
    root: Path, config: Mapping[str, Any]
) -> tuple[list[dict[str, str]], dict[str, Any]]:
    trees = config["frozen_tree_authorities"]
    common = root / trees["common_evaluation"]["path"] / "evaluation"
    caduceus = root / trees["caduceus_evaluation"]["path"] / "evaluation"
    common_receipt = load_json(common / "receipt.json", label="common evaluation receipt")
    caduceus_receipt = load_json(
        caduceus / "receipt.json", label="Caduceus evaluation receipt"
    )
    shared_fields = (
        "status",
        "dataset_id",
        "elements",
        "outer_locus_sequence_groups",
        "outer_folds",
        "contexts",
        "endpoint",
        "heads",
        "baselines",
        "head_fit_scope",
        "projection_fit_scope",
        "hyperparameter_selection_scope",
        "seed",
        "bootstrap_replicates",
        "bootstrap_unit",
        "donor_count",
        "experimental_replicates_per_context",
        "experimental_replicates_used_as_independent_donors",
        "outcome_role",
        "external_evaluation",
        "champion_eligible",
        "confirmatory_inference",
        "clinical_claim_supported",
        "multiple_testing",
        "sealed_outcomes_read",
    )
    for field in shared_fields:
        if common_receipt.get(field) != caduceus_receipt.get(field):
            raise DnaLmReconciliationError(
                f"evaluation receipt differs between lanes: {field}"
            )
    if (
        set(common_receipt.get("models", ())) != COMMON_MODELS
        or caduceus_receipt.get("models") != ["caduceus"]
        or tuple(common_receipt.get("heads", ())) != EXPECTED_HEADS
        or tuple(common_receipt.get("baselines", ())) != EXPECTED_BASELINES
        or tuple(common_receipt.get("contexts", ())) != EXPECTED_CONTEXTS
        or common_receipt.get("donor_count") != 0
        or common_receipt.get("experimental_replicates_used_as_independent_donors")
        is not False
        or common_receipt.get("champion_eligible") is not False
    ):
        raise DnaLmReconciliationError("evaluation topology or head census differs")

    common_metrics = read_tsv(common / "aggregate_metrics.tsv")
    caduceus_metrics = read_tsv(caduceus / "aggregate_metrics.tsv")
    key_fields = ("model_id", "head_id", "context_id")
    common_index = indexed_rows(common_metrics, key_fields)
    caduceus_index = indexed_rows(caduceus_metrics, key_fields)
    baseline_keys = {
        key for key in common_index if key[0] == "task_native_baseline"
    }
    if (
        baseline_keys
        != {key for key in caduceus_index if key[0] == "task_native_baseline"}
        or any(common_index[key] != caduceus_index[key] for key in baseline_keys)
    ):
        raise DnaLmReconciliationError("task-native baselines differ between lanes")

    common_promotion = read_tsv(common / "promotion.tsv")
    caduceus_promotion = read_tsv(caduceus / "promotion.tsv")
    promotion = indexed_rows(common_promotion + caduceus_promotion, ("model_id",))
    expected = config["expected_development_disposition"]
    if set(key[0] for key in promotion) != set(MODELS):
        raise DnaLmReconciliationError("promotion model census differs")
    output_rows: list[dict[str, str]] = []
    for model in MODELS:
        row = promotion[(model,)]
        model_expected = expected["models"][model]
        for field in (
            "selected_head",
            "mean_two_context_spearman",
            "spearman_gain",
            "three_seed_screening_recommended",
        ):
            if row.get(field) != model_expected[field]:
                raise DnaLmReconciliationError(
                    f"{model} development disposition differs: {field}"
                )
        if (
            row.get("allele_baseline_mean_spearman")
            != expected["allele_identity_ridge_mean_two_context_spearman"]
            or row.get("three_seed_screening_recommended") != "false"
            or float(row["spearman_gain"]) >= 0.0
        ):
            raise DnaLmReconciliationError(f"{model} promotion gate differs")
        output_rows.append({field: row[field] for field in PROMOTION_FIELDS})
    ranked = sorted(
        output_rows,
        key=lambda row: float(row["mean_two_context_spearman"]),
        reverse=True,
    )
    if (
        ranked[0]["model_id"] != expected["best_model_by_mean_spearman"]
        or ranked[0]["restricted_comparator"] != "true"
        or expected["family_promotion"] != "none"
    ):
        raise DnaLmReconciliationError("family ranking or restriction differs")
    return output_rows, {
        "baseline_rows_identical_between_evaluators": True,
        "identical_head_functions_reused": True,
        "head_and_hyperparameter_contract_identical": True,
        "contexts_identical": True,
        "outer_locus_sequence_groups": 1033,
        "outer_folds": 5,
        "biological_donors": 0,
        "experimental_replicates_per_context": 4,
        "best_model_by_mean_spearman": ranked[0]["model_id"],
        "best_model_is_restricted": True,
        "all_model_gains_over_allele_identity_ridge_are_negative": True,
        "three_seed_screening_recommended": [],
        "family_promotion": "none",
    }


def preflight(*, root: Path, config_path: Path) -> dict[str, Any]:
    root = reject_symlink_components(root, label="benchmark root").resolve(strict=True)
    config_path = reject_symlink_components(
        config_path, label="DNA LM reconciliation config"
    ).resolve(strict=True)
    try:
        config_path.relative_to(root)
    except ValueError as error:
        raise DnaLmReconciliationError("config escapes benchmark root") from error
    config = load_json(config_path, label="DNA LM reconciliation config")
    validate_config(config)
    authorities = validate_authorities(root, config)
    validate_shared_evaluator_code(root, config)
    model_rows = validate_model_authorities(root, config)
    embedding_audit = validate_embedding_contracts(root, config)
    promotion_rows, evaluation_audit = validate_evaluations(root, config)
    return {
        "schema_version": "masld-bench-dna-lm-family-reconciliation-receipt-v1",
        "status": config["status"],
        "family_id": "dna_language",
        "config_sha256": sha256_file(config_path),
        "authorities": authorities,
        "model_census": model_rows,
        "promotion_summary": promotion_rows,
        "task_id": config["task_contract"]["task_id"],
        "task_role": "exposed_development_MPRA_only",
        "embedding_audit": embedding_audit,
        "evaluation_audit": evaluation_audit,
        "family_native_input_lengths_bp": {
            model: config["models"][model]["task_input_bp"] for model in MODELS
        },
        "representation_context_matched": False,
        "downstream_head_contract_matched": True,
        "development_summary_metrics_read": True,
        "raw_reporter_outcomes_parsed": False,
        "sealed_assets_read": False,
        "supervised_eQTL_or_ieQTL_head_fit": False,
        "external_evaluation": False,
        "champion_eligible": False,
        "universal_claim_supported": False,
        "clinical_claim_supported": False,
        "next_action": config["expected_development_disposition"]["next_action"],
    }


def write_outputs(output: Path, receipt: Mapping[str, Any]) -> None:
    output = reject_symlink_components(output, label="DNA LM reconciliation output")
    if output.exists() or output.is_symlink():
        raise DnaLmReconciliationError(f"refusing to overwrite {output}")
    output.mkdir(parents=True, mode=0o750)
    for name, rows in (
        ("model_census.tsv", receipt["model_census"]),
        ("promotion_summary.tsv", receipt["promotion_summary"]),
    ):
        if not rows:
            raise DnaLmReconciliationError(f"empty output rows: {name}")
        with (output / name).open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
    frozen = dict(receipt)
    frozen["model_census_rows"] = len(receipt["model_census"])
    frozen["promotion_summary_rows"] = len(receipt["promotion_summary"])
    del frozen["model_census"]
    del frozen["promotion_summary"]
    frozen["model_census_sha256"] = sha256_file(output / "model_census.tsv")
    frozen["promotion_summary_sha256"] = sha256_file(
        output / "promotion_summary.tsv"
    )
    write_json_exclusive(output / "receipt.json", frozen, mode=0o640)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    arguments = parser.parse_args()
    receipt = preflight(root=arguments.root, config_path=arguments.config)
    write_outputs(arguments.output, receipt)
    printable = dict(receipt)
    printable["model_census"] = len(receipt["model_census"])
    printable["promotion_summary"] = len(receipt["promotion_summary"])
    print(json.dumps(printable, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
