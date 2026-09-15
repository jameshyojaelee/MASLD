#!/usr/bin/env python3
"""Freeze DNA-language feature coverage and the seeded MPRA head TaskSpec."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np

from masld_bench.artifacts import ArtifactError, verify_frozen_tree


SCHEMA = "masld-bench-gse281364-dna-language-seeded-head-campaign-v1"
MODELS = ("caduceus", "dnabert2", "hyenadna", "nucleotide_transformer")
OPEN_MODELS = ("caduceus", "dnabert2", "hyenadna")
RESTRICTED_MODELS = ("nucleotide_transformer",)
SEEDS = (1103, 2909, 4721, 6673, 8111)
INNER_BOOTSTRAP_SALT_OFFSET = 10000
OUTER_BOOTSTRAP_SALT_OFFSET = 20000
TRAINING_MEAN_BOOTSTRAP_SALT_OFFSET = 30000
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
HEADS = ("delta_ridge", "full_ridge")
CONTROLS = ("zero", "outer_training_mean", "allele_identity_ridge")
CANDIDATES = tuple(
    [("available_simple_controls", head) for head in CONTROLS]
    + [(model, head) for model in MODELS for head in HEADS]
)
MISSING_BASELINES = (
    "deltaSVM",
    "gkm-SVM",
    "sequence_CNN",
    "sequence_transformer",
    "MPRALegNet_signed_rekeyed_to_1033_elements_239_blocks",
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
    "prediction",
    "experimental_replicates",
    "biological_donors",
    "outcome_role",
)


class DNAFeatureReconciliationError(RuntimeError):
    """Raised when feature, checkpoint, runtime, or exposure status differs."""


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise DNAFeatureReconciliationError(f"invalid JSON {path}: {error}") from error
    if not isinstance(value, dict):
        raise DNAFeatureReconciliationError(f"JSON object required: {path}")
    return value


def validate_config(config: Mapping[str, Any]) -> None:
    authorization = config.get("authorization", {})
    fit = config.get("fit", {})
    firewall = config.get("completion_firewall", {})
    shared_control = config.get("shared_control_reference", {})
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != "prespecified_open_and_restricted_development_head_campaign"
        or config.get("dataset_id") != "gse281364"
        or config.get("task_id") != "variant_to_regulation"
        or authorization.get("outcome_blind_reprojection") is not True
        or authorization.get("outcome_access_after_taskspec_freeze") is not True
        or authorization.get("static_head_fit") is not True
        or any(
            authorization.get(field) is not False
            for field in (
                "base_checkpoint_tuning",
                "sealed_asset_access",
                "stack_fit",
                "global_census_mutation",
            )
        )
        or tuple(row.get("model_id") for row in config.get("models", ())) != MODELS
        or tuple(config.get("head_ids", ())) != HEADS
        or tuple(config.get("control_head_ids", ())) != CONTROLS
        or tuple(config.get("task", {}).get("seeds", ())) != SEEDS
        or tuple(config.get("task", {}).get("contexts", ())) != CONTEXTS
        or config.get("task", {}).get("elements") != 1033
        or config.get("task", {}).get("long_range_blocks") != 239
        or fit.get("five_seed_strategy") != "outer_training_long_range_block_bootstrap"
        or fit.get("bootstrap_draw_contract") != "exact_frozen_enformer_sei_offsets"
        or fit.get("shared_control_prediction_policy") != "reuse_bit_exact_frozen_enformer_sei_oof"
        or fit.get("inner_bootstrap_salt_offset") != INNER_BOOTSTRAP_SALT_OFFSET
        or fit.get("outer_bootstrap_salt_offset") != OUTER_BOOTSTRAP_SALT_OFFSET
        or fit.get("training_mean_bootstrap_salt_offset") != TRAINING_MEAN_BOOTSTRAP_SALT_OFFSET
        or fit.get("inner_projection") != "raw_embeddings_fit_on_inner_training_folds_only"
        or fit.get("outer_projection") != "raw_embeddings_fit_on_outer_training_folds_only"
        or fit.get("projection_width") != 256
        or fit.get("full_feature_width") != 1024
        or fit.get("held_fold_outcomes_used_for_fit_or_tuning") is not False
        or fit.get("held_fold_features_used_for_preprocessing_selection_or_tuning") is not False
        or firewall.get("mandatory_task_native_baselines_complete") is not False
        or tuple(firewall.get("missing_baselines", ())) != MISSING_BASELINES
        or firewall.get("models_reported_separately") is not True
        or shared_control.get("tree_path") != "executions/model-cpu-train-604-21097189/fit"
        or shared_control.get("artifacts_sha256") != "9306d6352163a4ccabf06406709d3545b0a73cca1b47e07ab817bcae874a829e"
        or shared_control.get("member") != "oof_predictions.tsv.gz"
        or shared_control.get("member_sha256") != "7b96aced402064e6a5d06a86e77d59cbce7be5377c3a47e30aadb6a96ad36c06"
        or shared_control.get("model_id") != "available_simple_controls"
        or shared_control.get("head_id") != "allele_identity_ridge"
        or shared_control.get("expected_rows") != 10330
        or any(
            firewall.get(field) is not False
            for field in (
                "cross_family_ranking",
                "shortlist_authorized",
                "champion_claim_authorized",
                "sealed_claim_authorized",
                "clinical_claim_authorized",
            )
        )
    ):
        raise DNAFeatureReconciliationError("campaign contract differs")


def load_config(path: Path) -> dict[str, Any]:
    config = load_json(path)
    validate_config(config)
    return config


def safe_tree(root: Path, binding: Mapping[str, Any], *, label: str) -> Path:
    value = binding.get("tree_path")
    if not isinstance(value, str) or not value:
        raise DNAFeatureReconciliationError(f"{label} tree path missing")
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise DNAFeatureReconciliationError(f"unsafe {label} tree path")
    try:
        tree = (root / relative).resolve(strict=True)
        tree.relative_to(root)
    except (OSError, ValueError) as error:
        raise DNAFeatureReconciliationError(f"{label} tree missing or escapes root") from error
    if file_sha256(tree / "ARTIFACTS.json") != binding.get("artifacts_sha256"):
        raise DNAFeatureReconciliationError(f"{label} ARTIFACTS identity differs")
    try:
        verify_frozen_tree(tree)
    except ArtifactError as error:
        raise DNAFeatureReconciliationError(f"{label} frozen tree differs: {error}") from error
    return tree


def write_tsv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise DNAFeatureReconciliationError(f"cannot write empty TSV: {path}")
    fields: list[str] = []
    for row in rows:
        for field in row:
            if field not in fields:
                fields.append(field)
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def raw_binding(config: Mapping[str, Any], model: str) -> tuple[Mapping[str, Any], Mapping[str, Any]]:
    source_name = next(row["raw_source"] for row in config["models"] if row["model_id"] == model)
    source = config[source_name]
    return source, source["models"][model]


def verify_raw_model(
    root: Path, config: Mapping[str, Any], model_record: Mapping[str, Any]
) -> tuple[dict[str, Any], dict[str, Any]]:
    model = model_record["model_id"]
    source, member = raw_binding(config, model)
    tree = safe_tree(root, source, label=f"{model}_raw")
    receipt_path = tree / f"raw/{model}/receipt.json"
    matrix_path = tree / f"raw/{model}/allele_embeddings.npz"
    if file_sha256(receipt_path) != member["receipt_sha256"] or file_sha256(matrix_path) != member["matrix_sha256"]:
        raise DNAFeatureReconciliationError(f"{model} raw member identity differs")
    receipt = load_json(receipt_path)
    if (
        receipt.get("status") != "pass_outcome_blind_embedding_extraction"
        or receipt.get("model_id") != model
        or receipt.get("elements") != 1033
        or receipt.get("outer_locus_sequence_groups") != 1033
        or receipt.get("checkpoint_sha256") != model_record["checkpoint_sha256"]
        or (
            model != "caduceus"
            and receipt.get("successful_probe_artifacts_sha256")
            != model_record["runtime_probe_artifacts_sha256"]
        )
        or receipt.get("license") != model_record["license"]
        or bool(receipt.get("restricted_comparator")) != bool(model_record["restricted"])
        or receipt.get("checkpoint_loaded_with_safetensors") is not True
        or receipt.get("runtime_network_allowed") is not False
        or receipt.get("reporter_counts_read") is not False
        or receipt.get("reporter_outcomes_read") is not False
        or receipt.get("sealed_labels_read") is not False
        or receipt.get("downstream_head_fit") is not False
    ):
        raise DNAFeatureReconciliationError(f"{model} raw receipt differs")
    with np.load(matrix_path, allow_pickle=False) as data:
        ids = data["fixture_ids"].astype(str)
        groups = data["outer_locus_sequence_group_ids"].astype(str)
        folds = data["outer_folds"].astype(np.int64)
        allele_order = tuple(data["allele_order"].astype(str).tolist())
        values = data["embeddings"]
    expected_shape = (1033, 4, int(model_record["hidden_width"]))
    if (
        ids.shape != (1033,)
        or len(set(ids.tolist())) != 1033
        or groups.shape != (1033,)
        or len(set(groups.tolist())) != 1033
        or folds.shape != (1033,)
        or set(folds.tolist()) != set(range(5))
        or allele_order != ("REF", "ALT", "REF_RC", "ALT_RC")
        or values.shape != expected_shape
        or not np.isfinite(values).all()
    ):
        raise DNAFeatureReconciliationError(f"{model} raw feature coverage differs")
    return receipt, {
        "model_id": model,
        "feature_stage": "raw_allele_embedding",
        "rows": 1033,
        "alleles": 4,
        "feature_width": values.shape[2],
        "outer_folds": 5,
        "long_range_blocks": "not_applicable_before_reprojection",
        "complete": "true",
        "member_sha256": member["matrix_sha256"],
        "missing_or_imputed_features": "false",
    }


def verify_projection_model(
    projection_root: Path, model: str, row_elements: set[str]
) -> list[dict[str, Any]]:
    receipt = load_json(projection_root / f"{model}/receipt.json")
    if (
        receipt.get("status") != "pass_outcome_blind_long_range_projection"
        or receipt.get("model_id") != model
        or receipt.get("elements") != 1033
        or receipt.get("long_range_blocks") != 239
        or receipt.get("outer_folds") != 5
        or receipt.get("source_groups_reassigned_from_original_fold_map") != 830
        or receipt.get("projection_width") != 256
        or receipt.get("head_input_width") != 1024
        or receipt.get("outcomes_read")
        or receipt.get("prediction_values_read")
        or receipt.get("downstream_head_fit")
        or receipt.get("sealed_assets_read")
    ):
        raise DNAFeatureReconciliationError(f"{model} projection receipt differs")
    rows = []
    for held_fold in range(5):
        path = projection_root / f"{model}/heldout_fold{held_fold}/head_features.npz"
        with np.load(path, allow_pickle=False) as data:
            elements = data["element_ids"].astype(str)
            blocks = data["long_range_block_ids"].astype(str)
            folds = data["outer_folds"].astype(np.int64)
            order = tuple(data["feature_block_order"].astype(str).tolist())
            features = data["features"]
        if (
            elements.shape != (1033,)
            or set(elements.tolist()) != row_elements
            or len(set(blocks.tolist())) != 239
            or set(folds.tolist()) != set(range(5))
            or order != ("REF", "ALT", "ALT_minus_REF", "absolute_ALT_minus_REF")
            or features.shape != (1033, 1024)
            or not np.isfinite(features).all()
        ):
            raise DNAFeatureReconciliationError(f"{model} projected feature coverage differs")
        rows.append({
            "model_id": model,
            "feature_stage": "long_range_outer_projection",
            "held_out_fold": held_fold,
            "rows": 1033,
            "feature_width": 1024,
            "outer_folds": 5,
            "long_range_blocks": 239,
            "complete": "true",
            "member_sha256": file_sha256(path),
            "missing_or_imputed_features": "false",
        })
    return rows


def reconcile(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.root.resolve(strict=True)
    config_path = arguments.config.resolve(strict=True)
    config = load_config(config_path)
    if arguments.output.exists():
        raise DNAFeatureReconciliationError("reconciliation output exists")
    row_tree = safe_tree(root, config["row_authority"], label="row_authority")
    row_path = row_tree / config["row_authority"]["member"]
    if file_sha256(row_path) != config["row_authority"]["member_sha256"]:
        raise DNAFeatureReconciliationError("row-universe member differs")
    with row_path.open("r", encoding="utf-8", newline="") as handle:
        row_rows = [dict(row) for row in csv.DictReader(handle, delimiter="\t")]
    row_elements = {row["element_id"] for row in row_rows}
    row_blocks = {row["long_range_block_id"] for row in row_rows}
    if len(row_rows) != 10330 or len(row_elements) != 1033 or len(row_blocks) != 239:
        raise DNAFeatureReconciliationError("row-universe denominator differs")

    raw_receipts = {}
    coverage_rows = []
    status_rows = []
    for model_record in config["models"]:
        model = model_record["model_id"]
        receipt, coverage = verify_raw_model(root, config, model_record)
        raw_receipts[model] = receipt
        coverage_rows.append(coverage)
        runtime_tree = safe_tree(
            root,
            {
                "tree_path": model_record["runtime_probe_tree_path"],
                "artifacts_sha256": model_record["runtime_probe_artifacts_sha256"],
            },
            label=f"{model}_runtime_probe",
        )
        checkpoint_path = root / f"config/artifacts/models/{model}/checkpoints.json"
        exposure_path = root / f"config/artifacts/models/{model}/exposure_audit.json"
        if (
            file_sha256(checkpoint_path) != model_record["checkpoint_manifest_sha256"]
            or file_sha256(exposure_path) != model_record["exposure_manifest_sha256"]
        ):
            raise DNAFeatureReconciliationError(f"{model} checkpoint/exposure manifest differs")
        exposure = load_json(exposure_path)["checkpoint_finding"]
        if exposure.get("exposure_state") != model_record["exposure"] or exposure.get("target_effect_label_exposure") != "unexposed":
            raise DNAFeatureReconciliationError(f"{model} exposure status differs")
        status_rows.append({
            "model_id": model,
            "family_native_input": model_record["family_native_input"],
            "checkpoint_sha256": model_record["checkpoint_sha256"],
            "checkpoint_manifest_sha256": model_record["checkpoint_manifest_sha256"],
            "runtime_probe_artifacts_sha256": model_record["runtime_probe_artifacts_sha256"],
            "runtime_probe_tree_path": runtime_tree.relative_to(root).as_posix(),
            "runtime_status": "passed_exact_checkpoint_offline_embedding_extraction",
            "exposure": model_record["exposure"],
            "reference_sequence_exposure": exposure.get("reference_sequence_exposure"),
            "target_effect_label_exposure": exposure.get("target_effect_label_exposure"),
            "license": model_record["license"],
            "restricted": str(bool(model_record["restricted"])).lower(),
            "open_champion_eligible_after_external_gates": str(bool(model_record["open_champion_eligible_after_external_gates"])).lower(),
            "raw_feature_complete": "true",
            "long_range_projection_complete": "true",
            "blocked_feature_reason": "none",
        })

    open_tree = safe_tree(root, config["open_projection_authority"], label="open_projection_authority")
    open_root = open_tree / "open_sequence_projection"
    open_receipt = load_json(open_root / "receipt.json")
    if tuple(open_receipt.get("models", ())) != OPEN_MODELS:
        raise DNAFeatureReconciliationError("open projection roster differs")
    nt_tree = arguments.nt_projection_tree.resolve(strict=True)
    if file_sha256(nt_tree / "ARTIFACTS.json") != arguments.nt_projection_sha256:
        raise DNAFeatureReconciliationError("Nucleotide Transformer projection identity differs")
    verify_frozen_tree(nt_tree)
    nt_receipt = load_json(nt_tree / "receipt.json")
    if tuple(nt_receipt.get("models", ())) != RESTRICTED_MODELS:
        raise DNAFeatureReconciliationError("restricted projection roster differs")
    for model in OPEN_MODELS:
        coverage_rows.extend(verify_projection_model(open_root, model, row_elements))
    coverage_rows.extend(verify_projection_model(nt_tree, "nucleotide_transformer", row_elements))

    blockers = [
        {"scope": "campaign", "blocker": value, "feature_incomplete": "false"}
        for value in MISSING_BASELINES
    ]
    blockers.extend(
        [
            {"scope": "nucleotide_transformer", "blocker": "CC-BY-NC-SA-4.0_restricted_comparator_only", "feature_incomplete": "false"},
            {"scope": "all_models", "blocker": "no_external_or_sealed_evaluation_in_this_development_campaign", "feature_incomplete": "false"},
            {"scope": "all_models", "blocker": "no_cross_family_shortlist_until_mandatory_baselines_finish", "feature_incomplete": "false"},
        ]
    )
    arguments.output.mkdir(parents=True, mode=0o750)
    write_tsv(arguments.output / "model_status.tsv", status_rows)
    write_tsv(arguments.output / "feature_coverage.tsv", coverage_rows)
    write_tsv(arguments.output / "blockers.tsv", blockers)
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_complete_feature_reconciliation_and_prespecified_taskspec",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "models": list(MODELS),
        "open_models": list(OPEN_MODELS),
        "restricted_models": list(RESTRICTED_MODELS),
        "candidate_count": len(CANDIDATES),
        "elements": 1033,
        "long_range_blocks": 239,
        "outer_folds": 5,
        "fixed_seeds": list(SEEDS),
        "feature_coverage_rows": len(coverage_rows),
        "feature_complete_models": list(MODELS),
        "feature_incomplete_models": [],
        "features_imputed": False,
        "obsolete_6kb_folds_reused": False,
        "source_groups_reassigned_from_obsolete_folds": 830,
        "config_sha256": file_sha256(config_path),
        "open_projection_artifacts_sha256": config["open_projection_authority"]["artifacts_sha256"],
        "nucleotide_transformer_projection_artifacts_sha256": arguments.nt_projection_sha256,
        "outcomes_read": False,
        "prediction_values_read": False,
        "model_fit": False,
        "sealed_assets_read": False,
        "mandatory_task_native_baselines_complete": False,
        "missing_baselines": list(MISSING_BASELINES),
        "models_ranked": False,
        "shortlist_created": False,
        "champion_claim": False,
    }
    (arguments.output / "receipt.json").write_text(json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--nt-projection-tree", type=Path, required=True)
    parser.add_argument("--nt-projection-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    reconcile(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
