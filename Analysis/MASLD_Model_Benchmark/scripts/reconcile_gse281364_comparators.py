#!/usr/bin/env python3
"""Re-key static GSE281364 comparators without reading outcomes or metrics."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import json
import math
from pathlib import Path
import re
import tomllib
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import verify_frozen_tree


SCHEMA = "masld-bench-gse281364-comparator-reconciliation-v1"
STATUS = "prespecified_outcome_blind_join_coverage_disposition"
SEEDS = (1103, 2909, 4721, 6673, 8111)
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
MODELS = (
    "mpralegnet",
    "sequence_cnn_control",
    "sequence_transformer_control",
    "gkmsvm",
    "deltasvm",
)
ROW_FIELDS = (
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
)
MPRA_FIELDS = (
    "element_id",
    "outer_locus_sequence_group_id",
    "outer_fold",
    "contig",
    "variant_pos0",
    "ref",
    "alt",
    "cell_context",
    "input_lane",
    "reference_forward_score",
    "alternative_forward_score",
    "reference_reverse_complement_score",
    "alternative_reverse_complement_score",
    "reference_orientation_mean_score",
    "alternative_orientation_mean_score",
    "alt_minus_ref_score",
)
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class ReconciliationError(RuntimeError):
    """Raised when a join, authority, or outcome separation differs."""


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ReconciliationError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise ReconciliationError(f"{label} must be an object")
    return value


def load_config(path: Path) -> dict[str, Any]:
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise ReconciliationError(f"invalid reconciliation config: {error}") from error
    if not isinstance(value, dict):
        raise ReconciliationError("reconciliation config must be a table")
    return value


def read_tsv(path: Path) -> list[dict[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            return [dict(row) for row in reader]
    except (OSError, UnicodeDecodeError, csv.Error) as error:
        raise ReconciliationError(f"cannot read TSV {path}: {error}") from error


def write_tsv(path: Path, fields: Sequence[str], rows: Sequence[Mapping[str, Any]]) -> None:
    if not rows:
        raise ReconciliationError(f"cannot write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def safe_path(root: Path, value: object, *, label: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ReconciliationError(f"{label} path is missing")
    relative = Path(value)
    if relative.is_absolute() or ".." in relative.parts:
        raise ReconciliationError(f"unsafe {label} path")
    try:
        path = (root / relative).resolve(strict=True)
        path.relative_to(root)
    except (OSError, ValueError) as error:
        raise ReconciliationError(f"{label} path is unavailable") from error
    return path


def verify_tree(root: Path, authority: Mapping[str, Any], *, label: str) -> Path:
    expected = authority.get("artifacts_sha256")
    if not isinstance(expected, str) or SHA256.fullmatch(expected) is None:
        raise ReconciliationError(f"{label} ARTIFACTS hash is invalid")
    tree = safe_path(root, authority.get("tree_path"), label=label)
    if file_sha256(tree / "ARTIFACTS.json") != expected:
        raise ReconciliationError(f"{label} ARTIFACTS hash differs")
    verify_frozen_tree(tree)
    return tree


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != STATUS
        or config.get("dataset_id") != "gse281364"
        or config.get("task_id") != "variant_to_regulation"
        or config.get("outcome_access_authorized") is not False
        or config.get("model_fit_authorized") is not False
        or config.get("prediction_generation_authorized") is not False
        or config.get("prediction_rekey_authorized") is not True
        or config.get("metric_calculation_authorized") is not False
        or config.get("shortlist_authorized") is not False
        or config.get("complementarity_authorized") is not False
        or config.get("conditional_trigger_authorized") is not False
        or config.get("sealed_asset_access_authorized") is not False
    ):
        raise ReconciliationError("authorization firewall differs")
    row = config.get("row_contract", {})
    if (
        row.get("selected_elements") != 1033
        or row.get("source_locus_groups") != 1033
        or row.get("long_range_blocks") != 239
        or row.get("outer_folds") != 5
        or tuple(row.get("contexts", ())) != CONTEXTS
        or tuple(row.get("fixed_schema_seeds", ())) != SEEDS
        or row.get("row_keys") != 10330
        or row.get("biological_donors") != 0
    ):
        raise ReconciliationError("row contract differs")
    completion = config.get("completion_firewall", {})
    if (
        completion.get("campaign_scope")
        != "outcome_blind_comparator_reconciliation_only"
        or completion.get("mandatory_comparator_evaluation_complete") is not False
        or completion.get("candidate_summary_metric_values_read") is not False
        or completion.get("deterministic_schema_seed_repeats_are_independent_fits")
        is not False
        or completion.get("stochastic_stability_claim_blocked_for_deterministic_heads")
        is not True
        or completion.get("shortlist_blocked") is not True
        or completion.get("finalist_claim_blocked") is not True
        or completion.get("complementarity_blocked") is not True
        or completion.get("conditional_trigger_blocked") is not True
        or completion.get("champion_claim_blocked") is not True
    ):
        raise ReconciliationError("completion firewall differs")
    comparators = config.get("comparator", ())
    if len(comparators) != 5 or tuple(row.get("model_id") for row in comparators) != MODELS:
        raise ReconciliationError("comparator roster differs")


def validate_file_authorities(root: Path, config: Mapping[str, Any]) -> dict[str, Path]:
    paths: dict[str, Path] = {}
    for label, authority in config.get("file_authorities", {}).items():
        path = safe_path(root, authority.get("path"), label=label)
        if path.is_symlink() or not path.is_file() or file_sha256(path) != authority.get("sha256"):
            raise ReconciliationError(f"file authority differs: {label}")
        paths[label] = path
    if len(paths) != 12:
        raise ReconciliationError("file-authority roster differs")
    return paths


def validate_registry_semantics(paths: Mapping[str, Path]) -> dict[str, Any]:
    mpra_checkpoint = load_json(paths["mpralegnet_checkpoints"], label="MPRALegNet checkpoint")
    mpra_exposure = load_json(paths["mpralegnet_exposure"], label="MPRALegNet exposure")
    controls = load_json(paths["sequence_controls_checkpoints"], label="sequence controls")
    control_exposure = load_json(paths["sequence_controls_exposure"], label="sequence control exposure")
    gkm = load_json(paths["gkmsvm_checkpoints"], label="gkm-SVM checkpoint")
    delta = load_json(paths["deltasvm_checkpoints"], label="deltaSVM checkpoint")
    gkm_exposure = load_json(paths["gkmsvm_exposure"], label="gkm-SVM exposure")
    delta_exposure = load_json(paths["deltasvm_exposure"], label="deltaSVM exposure")
    mpra_finding = mpra_exposure.get("checkpoint_finding", {})
    mpra_terms = mpra_checkpoint.get("terms_audit", {})
    control_models = controls.get("model_specifications", controls.get("models", {}))
    if not isinstance(control_models, dict):
        raise ReconciliationError("sequence-control model registry differs")
    cnn = control_models.get("sequence_cnn_control", {})
    transformer = control_models.get("sequence_transformer_control", {})
    if (
        mpra_checkpoint.get("code", {}).get("license") != "MIT"
        or mpra_checkpoint.get("weight_license")
        != "MIT_declared_by_checkpoint_model_card"
        or mpra_terms.get("result")
        != "conditional_open_artifact_pending_provenance_equivalence_or_explicit_acceptance_of_third_party_weights"
        or mpra_finding.get("overall_exposure_state") != "target_label_unexposed"
        or mpra_finding.get("gse281364_target_exposure") != "clean_declared"
        or cnn.get("implementation_status")
        != "terminal_development_production_bundle_registry_blocked"
        or cnn.get("registry_license_status") != "UNRESOLVED"
        or cnn.get("registry_exposure_status") != "unknown"
        or transformer.get("implementation_status")
        != "not_implemented_architecture_unresolved"
        or control_exposure.get("license_and_release_disposition", {}).get(
            "project_implementation_license"
        )
        != "UNRESOLVED_for_sequence_cnn_control_despite_exact_code_artifact"
        or gkm.get("runtime_contract", {}).get("environment_status")
        != "UNRESOLVED_CPU_SLURM_BUILD_NUMERIC_AND_RESUME_PROBES_REQUIRED"
        or delta.get("runtime_contract", {}).get("environment_status")
        != "UNRESOLVED_CPU_SLURM_BUILD_NUMERIC_AND_RESUME_PROBES_REQUIRED"
        or "not_yet_created"
        not in gkm_exposure.get("checkpoint_finding", {}).get("checkpoint_identity", "")
        or "not_yet_created"
        not in delta_exposure.get("checkpoint_finding", {}).get("checkpoint_identity", "")
    ):
        raise ReconciliationError("license, exposure, or implementation semantics differ")
    return {
        "mpralegnet_code_license": "MIT",
        "mpralegnet_weight_license": "MIT_declared_by_checkpoint_model_card",
        "mpralegnet_exposure": "target_label_unexposed",
        "mpralegnet_gse281364_exposure": "clean_declared",
        "sequence_cnn_license": "UNRESOLVED",
        "sequence_cnn_exposure": "unknown",
        "sequence_transformer_implementation": "not_implemented_architecture_unresolved",
        "gkmsvm_runtime": "unresolved_no_project_model",
        "deltasvm_runtime": "unresolved_no_shared_project_model_or_wrapper",
    }


def read_candidate_summary_header(path: Path) -> tuple[str, ...]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            line = handle.readline().rstrip("\n")
    except (OSError, UnicodeDecodeError) as error:
        raise ReconciliationError(f"cannot read candidate-summary header: {error}") from error
    fields = tuple(line.split("\t"))
    required = {
        "model_id",
        "head_id",
        "positive_gain_seeds",
        "development_shortlist",
        "shortlist_blocked",
        "finalist_claim",
    }
    if not required.issubset(fields):
        raise ReconciliationError("candidate-summary interpretation fields differ")
    return fields


def seed_interpretation_audit(
    production_root: Path, task_spec: Mapping[str, Any]
) -> list[dict[str, Any]]:
    artifacts = load_json(production_root / "ARTIFACTS.json", label="partial campaign artifacts")
    records = {row["path"]: row["sha256"] for row in artifacts.get("artifacts", ())}
    declared = {
        (candidate["model_id"], head): candidate["seed_behavior"]
        for candidate in task_spec.get("candidates", ())
        for head in candidate["head_ids"]
    }
    expected = {
        **{("available_simple_controls", head): "deterministic_repeated_for_schema_only" for head in ("zero", "training_mean", "allele_identity_ridge")},
        **{
            (model, head): "ridge_deterministic_two_layer_seed_specific"
            for model in ("caduceus", "dnabert2", "hyenadna")
            for head in ("delta_ridge", "full_ridge", "two_layer_gelu")
        },
    }
    if declared != expected:
        raise ReconciliationError("partial-campaign declared seed behavior differs")
    rows: list[dict[str, Any]] = []
    for model_id, head_id in expected:
        stochastic = head_id == "two_layer_gelu"
        hash_audit = "contract_only_no_saved_state"
        if head_id not in {"zero", "training_mean"}:
            for context in CONTEXTS:
                for fold in range(5):
                    hashes = []
                    for seed in SEEDS:
                        path = f"fit/heads/{model_id}/{head_id}/{context}/seed{seed}_fold{fold}.npz"
                        if path not in records:
                            raise ReconciliationError(f"partial-campaign head state is missing: {path}")
                        hashes.append(records[path])
                    expected_unique = 5 if stochastic else 1
                    if len(set(hashes)) != expected_unique:
                        raise ReconciliationError("seed-specific state-hash behavior differs")
            hash_audit = (
                "five_distinct_seed_states_per_context_fold"
                if stochastic
                else "identical_state_hash_across_five_schema_seeds_per_context_fold"
            )
        rows.append(
            {
                "model_id": model_id,
                "head_id": head_id,
                "declared_seed_behavior": expected[(model_id, head_id)],
                "schema_seed_rows": 5,
                "independent_fitted_predictions": 5 if stochastic else 1,
                "state_hash_audit": hash_audit,
                "positive_gain_seeds_is_stochastic_stability_evidence": str(stochastic).lower(),
                "stability_claim_blocked": str(not stochastic).lower(),
                "interpretation": (
                    "five_seed_specific_two_layer_fits_but_all_shortlist_and_finalist_claims_remain_blocked"
                    if stochastic
                    else "five_schema_repeats_of_one_deterministic_fit_not_four_of_five_stability_evidence"
                ),
            }
        )
    return rows


def reconcile(arguments: argparse.Namespace) -> dict[str, Any]:
    root = arguments.root.resolve(strict=True)
    config_path = arguments.config.resolve(strict=True)
    config = load_config(config_path)
    validate_config(config)
    if arguments.output.exists():
        raise ReconciliationError("reconciliation output exists")
    files = validate_file_authorities(root, config)
    registry = validate_registry_semantics(files)
    trees = {
        label: verify_tree(root, config[label], label=label)
        for label in (
            "row_authority",
            "mpralegnet_prediction_authority",
            "sequence_cnn_authority",
            "sequence_transformer_authority",
            "partial_campaign_authority",
        )
    }
    for label, member_key, sha_key in (
        ("row_authority", "member", "member_sha256"),
        ("mpralegnet_prediction_authority", "score_member", "score_member_sha256"),
        ("mpralegnet_prediction_authority", "receipt_member", "receipt_member_sha256"),
        ("sequence_cnn_authority", "receipt_member", "receipt_member_sha256"),
        ("sequence_transformer_authority", "receipt_member", "receipt_member_sha256"),
        ("partial_campaign_authority", "task_spec_member", "task_spec_member_sha256"),
        ("partial_campaign_authority", "candidate_summary_member", "candidate_summary_member_sha256"),
    ):
        authority = config[label]
        member = trees[label] / authority[member_key]
        if file_sha256(member) != authority[sha_key]:
            raise ReconciliationError(f"authority member differs: {label}/{member_key}")

    row_rows = read_tsv(trees["row_authority"] / config["row_authority"]["member"])
    if not row_rows or tuple(row_rows[0]) != ROW_FIELDS or len(row_rows) != 10330:
        raise ReconciliationError("row-universe schema or denominator differs")
    element_metadata: dict[str, dict[str, str]] = {}
    row_identities = set()
    for row in row_rows:
        identity = (int(row["seed"]), row["row_hash"])
        if (
            int(row["seed"]) not in SEEDS
            or row["assay_context_id"] not in CONTEXTS
            or identity in row_identities
        ):
            raise ReconciliationError("row-universe identity differs")
        row_identities.add(identity)
        metadata = {
            "source_locus_group_id": row["source_locus_group_id"],
            "long_range_block_id": row["long_range_block_id"],
            "outer_fold": row["outer_fold"],
        }
        if element_metadata.setdefault(row["element_id"], metadata) != metadata:
            raise ReconciliationError("element metadata differs across schema rows")
    if len(element_metadata) != 1033 or len({row["long_range_block_id"] for row in element_metadata.values()}) != 239:
        raise ReconciliationError("selected-element or long-range-block census differs")

    mpra_rows = read_tsv(
        trees["mpralegnet_prediction_authority"]
        / config["mpralegnet_prediction_authority"]["score_member"]
    )
    mpra_receipt = load_json(
        trees["mpralegnet_prediction_authority"]
        / config["mpralegnet_prediction_authority"]["receipt_member"],
        label="MPRALegNet prediction receipt",
    )
    if (
        not mpra_rows
        or tuple(mpra_rows[0]) != MPRA_FIELDS
        or len(mpra_rows) != 4359
        or mpra_receipt.get("status") != "pass_outcome_blind_prediction"
        or mpra_receipt.get("elements") != 4359
        or mpra_receipt.get("outer_locus_sequence_groups") != 1033
        or mpra_receipt.get("allele_effect_sign") != "ALT_minus_REF"
        or mpra_receipt.get("outcomes_read") is not False
        or mpra_receipt.get("reporter_counts_read") is not False
        or mpra_receipt.get("sealed_outcomes_read") is not False
        or mpra_receipt.get("standalone_champion_eligible") is not False
    ):
        raise ReconciliationError("MPRALegNet source receipt or denominator differs")
    mpra_by_element: dict[str, dict[str, str]] = {}
    for row in mpra_rows:
        if row["element_id"] in mpra_by_element:
            raise ReconciliationError("duplicate MPRALegNet element")
        values = [
            float(row[field])
            for field in (
                "reference_forward_score",
                "alternative_forward_score",
                "reference_reverse_complement_score",
                "alternative_reverse_complement_score",
                "reference_orientation_mean_score",
                "alternative_orientation_mean_score",
                "alt_minus_ref_score",
            )
        ]
        if (
            not all(math.isfinite(value) for value in values)
            or row["cell_context"] != "HepG2"
            or row["input_lane"] != "GRCh38p14_230bp_genomic_window_transfer"
            or row["ref"] not in "ACGT"
            or row["alt"] not in "ACGT"
            or row["ref"] == row["alt"]
            or abs(values[4] - (values[0] + values[2]) / 2.0) > 1.0e-6
            or abs(values[5] - (values[1] + values[3]) / 2.0) > 1.0e-6
            or abs(values[6] - (values[5] - values[4])) > 1.0e-6
        ):
            raise ReconciliationError("MPRALegNet sign, strand, allele, or score semantics differ")
        mpra_by_element[row["element_id"]] = row
    if set(element_metadata) - set(mpra_by_element):
        raise ReconciliationError("selected row universe is not covered by MPRALegNet")

    selected_rows = []
    reassigned = 0
    for element in sorted(element_metadata):
        metadata = element_metadata[element]
        source = mpra_by_element[element]
        if source["outer_locus_sequence_group_id"] != metadata["source_locus_group_id"]:
            raise ReconciliationError("MPRALegNet source-locus group differs")
        old_fold = int(source["outer_fold"])
        new_fold = int(metadata["outer_fold"].removeprefix("fold-"))
        reassigned += old_fold != new_fold
        selected_rows.append(
            {
                "element_id": element,
                "source_locus_group_id": metadata["source_locus_group_id"],
                "long_range_block_id": metadata["long_range_block_id"],
                "outer_fold": metadata["outer_fold"],
                "original_6kb_outer_fold": old_fold,
                "model_id": "mpralegnet_hepg2_test1_val2",
                "native_score_field": "mean_orientation_ALT_minus_REF_uncalibrated_lentiMPRA_reporter_expression",
                "allele_effect_sign": "ALT_minus_REF",
                "alt_minus_ref_score": format(float(source["alt_minus_ref_score"]), ".17g"),
                "cell_context": "HepG2",
                "input_lane": source["input_lane"],
                "condition_specific": "false",
                "seed_specific": "false",
                "calibrated_to_gse281364_outcome": "false",
                "checkpoint_exposure": registry["mpralegnet_exposure"],
                "code_license": registry["mpralegnet_code_license"],
                "weight_license": registry["mpralegnet_weight_license"],
                "standalone_champion_eligible": "false",
            }
        )
    if reassigned != 830:
        raise ReconciliationError("MPRALegNet long-range fold reassignment census differs")
    selected_lookup = {row["element_id"]: row for row in selected_rows}
    joined_rows = []
    for row in row_rows:
        selected = selected_lookup[row["element_id"]]
        joined_rows.append(
            {
                **row,
                "model_id": selected["model_id"],
                "native_alt_minus_ref_score": selected["alt_minus_ref_score"],
                "native_score_repeated_across_contexts": "true",
                "native_score_repeated_across_schema_seeds": "true",
                "independent_fitted_predictions": 1,
                "condition_interaction_modeled": "false",
                "outcome_calibration_state": "not_calibrated",
                "future_role": "fixed_secondary_native_reporter_activity_comparator_only",
            }
        )

    cnn_receipt = load_json(
        trees["sequence_cnn_authority"] / config["sequence_cnn_authority"]["receipt_member"],
        label="sequence CNN bundle",
    )
    transformer_receipt = load_json(
        trees["sequence_transformer_authority"]
        / config["sequence_transformer_authority"]["receipt_member"],
        label="sequence transformer probe",
    )
    if (
        cnn_receipt.get("terminal_status") != "passed"
        or cnn_receipt.get("logical_tasks") != 8
        or transformer_receipt.get("status") != "pass"
        or transformer_receipt.get("models")
        != ["sequence_cnn_control", "sequence_transformer_control"]
        or transformer_receipt.get("test_outcomes_used") is not False
        or transformer_receipt.get("benchmark_metrics_calculated") is not False
    ):
        raise ReconciliationError("sequence-control non-GSE281364 evidence differs")

    production = trees["partial_campaign_authority"]
    task_spec = load_json(
        production / config["partial_campaign_authority"]["task_spec_member"],
        label="partial campaign TaskSpec",
    )
    candidate_header = read_candidate_summary_header(
        production / config["partial_campaign_authority"]["candidate_summary_member"]
    )
    seed_audit_rows = seed_interpretation_audit(production, task_spec)

    dispositions = []
    for comparator in config["comparator"]:
        dispositions.append(
            {
                "model_id": comparator["model_id"],
                "disposition": comparator["expected_disposition"],
                "joined_row_keys": comparator["expected_joined_row_keys"],
                "unique_native_scores": comparator["expected_unique_native_scores"],
                "valid_for_future_native_descriptive_scoring": str(
                    comparator["valid_for_future_native_descriptive_scoring"]
                ).lower(),
                "valid_for_current_shortlist": "false",
                "metrics_calculated": "false",
                "outcomes_read": "false",
            }
        )

    arguments.output.mkdir(parents=True, mode=0o750)
    write_tsv(arguments.output / "mpralegnet_selected_elements.tsv", tuple(selected_rows[0]), selected_rows)
    write_tsv(arguments.output / "mpralegnet_row_key_join.tsv", tuple(joined_rows[0]), joined_rows)
    write_tsv(arguments.output / "comparator_disposition.tsv", tuple(dispositions[0]), dispositions)
    write_tsv(arguments.output / "terminal_seed_interpretation_audit.tsv", tuple(seed_audit_rows[0]), seed_audit_rows)
    interpretation = {
        "schema_version": SCHEMA,
        "status": "pass_terminal_seed_interpretation_audit",
        "partial_campaign_artifacts_sha256": config["partial_campaign_authority"]["artifacts_sha256"],
        "candidate_summary_member_sha256": config["partial_campaign_authority"]["candidate_summary_member_sha256"],
        "candidate_summary_header_fields": list(candidate_header),
        "candidate_summary_header_read": True,
        "candidate_summary_metric_values_read": False,
        "deterministic_schema_seed_repeats_are_independent_fits": False,
        "deterministic_heads_independent_fitted_predictions": 1,
        "positive_gain_seeds_for_deterministic_heads_is_stability_evidence": False,
        "four_of_five_stability_gate_satisfied_by_deterministic_repeats": False,
        "provisional_partial_campaign_leader_may_be_used_for_finalist_claim": False,
        "mandatory_comparator_evaluation_complete": False,
        "shortlist_blocked": True,
        "finalist_claim_blocked": True,
        "complementarity_blocked": True,
    }
    (arguments.output / "terminal_interpretation_audit.json").write_text(
        json.dumps(interpretation, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    receipt = {
        "schema_version": SCHEMA,
        "status": "pass_outcome_blind_join_coverage_disposition",
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "selected_elements": 1033,
        "source_locus_groups": 1033,
        "long_range_blocks": 239,
        "row_keys": 10330,
        "mpralegnet_source_elements": 4359,
        "mpralegnet_selected_elements": 1033,
        "mpralegnet_join_coverage": 1.0,
        "mpralegnet_reassigned_from_original_fold": reassigned,
        "mpralegnet_context_specific": False,
        "mpralegnet_seed_specific": False,
        "mpralegnet_calibrated_to_gse281364_outcome": False,
        "mpralegnet_exposure": registry["mpralegnet_exposure"],
        "mpralegnet_code_license": registry["mpralegnet_code_license"],
        "mpralegnet_weight_license": registry["mpralegnet_weight_license"],
        "mpralegnet_checkpoint_provenance_equivalence_resolved": False,
        "sequence_cnn_gse281364_prediction_rows": 0,
        "sequence_transformer_gse281364_prediction_rows": 0,
        "gkmsvm_gse281364_prediction_rows": 0,
        "deltasvm_gse281364_prediction_rows": 0,
        "gkmsvm_and_deltasvm_independent_evidence": False,
        "candidate_summary_header_read": True,
        "candidate_summary_metric_values_read": False,
        "metrics_calculated": False,
        "outcomes_read": False,
        "reporter_counts_read": False,
        "model_fit": False,
        "new_predictions_generated": False,
        "sealed_assets_read": False,
        "shortlist_blocked": True,
        "finalist_claim_blocked": True,
        "complementarity_blocked": True,
        "conditional_trigger_blocked": True,
        "champion_claim": False,
        "mandatory_comparator_evaluation_complete": False,
        "config_sha256": file_sha256(config_path),
    }
    (arguments.output / "receipt.json").write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    reconcile(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
