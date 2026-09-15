#!/usr/bin/env python3
"""Freeze the outcome-blind GSE281364 sequence/context refit DAG."""

from __future__ import annotations

import argparse
import csv
from hashlib import sha256
import io
import json
from pathlib import Path
import re
import tomllib
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import write_json_exclusive, write_text_exclusive


SCHEMA = "masld-bench-gse281364-sequence-context-refit-dag-v1"
STATUS = "prepared_open_sequence_refit_context_restricted_no_trigger_pair"
SEEDS = (1103, 2909, 4721, 6673, 8111)
CONTEXTS = ("HepG2_control", "HepG2_PAOA")
COMPONENT_FIELDS = (
    "seed",
    "row_hash",
    "unit_hash",
    "block_hash",
    "stratum",
    "outer_fold",
    "study_id",
    "observed",
    "sequence_prediction",
    "context_prediction",
)
MODEL_IDS = {
    "borzoi_ensemble",
    "caduceus",
    "corgi_regular",
    "dnabert2",
    "enformer_crested_restricted_port",
    "hyenadna",
    "mpralegnet",
    "nucleotide_transformer",
    "sei",
    "sequence_cnn_control",
    "sequence_transformer_control",
    "zero_training_mean_allele_identity_ridge",
}
READY_STATES = {
    "complete_after_this_validator",
    "ready_after_projection_entrypoint_validation",
    "ready_after_context_roster_entrypoint_validation",
}
INVENTORY_FIELDS = (
    "artifact_id",
    "model_id",
    "family",
    "source_role",
    "tree_path",
    "artifacts_sha256",
    "receipt_member",
    "current_rows",
    "current_seeds",
    "current_fold_dependency",
    "action",
    "license",
    "exposure",
    "open_after_refit_and_external_gates",
    "diagnostic_only",
)
STAGE_FIELDS = (
    "order",
    "stage_id",
    "lane",
    "state",
    "dependencies",
    "entrypoint",
    "command",
    "outcomes_required",
    "gpu_required",
    "cpus",
    "memory_gb",
    "hours",
    "output_role",
    "blockers",
)
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class RefitDagError(RuntimeError):
    """Raised when the outcome separation, inventory, or refit DAG differs."""


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise RefitDagError(f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise RefitDagError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise RefitDagError(f"{label} must be a JSON object")
    return value


def load_config(path: Path) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise RefitDagError("config is not a regular file")
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise RefitDagError(f"invalid config: {error}") from error
    if not isinstance(value, dict):
        raise RefitDagError("config must be a TOML table")
    return value


def safe_path(root: Path, relative_value: object, *, label: str) -> Path:
    if not isinstance(relative_value, str) or not relative_value:
        raise RefitDagError(f"{label} path is missing")
    relative = Path(relative_value)
    if relative.is_absolute() or ".." in relative.parts:
        raise RefitDagError(f"unsafe {label} path")
    try:
        path = (root / relative).resolve(strict=True)
        path.relative_to(root)
    except (OSError, ValueError) as error:
        raise RefitDagError(f"{label} path is missing or escapes root") from error
    return path


def manifest_authority(
    root: Path, relative_path: object, expected_sha256: object, *, label: str
) -> tuple[Path, dict[str, Any], dict[str, Mapping[str, Any]]]:
    if not isinstance(expected_sha256, str) or not SHA256.fullmatch(expected_sha256):
        raise RefitDagError(f"{label} artifacts hash differs")
    tree = safe_path(root, relative_path, label=label)
    if tree.is_symlink() or not tree.is_dir():
        raise RefitDagError(f"{label} is not a frozen directory")
    manifest_path = tree / "ARTIFACTS.json"
    complete_path = tree / "COMPLETE"
    if file_sha256(manifest_path) != expected_sha256:
        raise RefitDagError(f"{label} ARTIFACTS hash differs")
    manifest = load_json(manifest_path, label=f"{label} ARTIFACTS")
    complete = load_json(complete_path, label=f"{label} COMPLETE")
    artifacts = manifest.get("artifacts")
    if (
        manifest.get("schema_version") != "masld-bench-artifacts-v1"
        or not isinstance(artifacts, list)
        or complete
        != {
            "artifact_count": len(artifacts),
            "manifest_sha256": expected_sha256,
            "schema_version": "masld-bench-complete-v1",
        }
    ):
        raise RefitDagError(f"{label} frozen-tree contract differs")
    members: dict[str, Mapping[str, Any]] = {}
    for item in artifacts:
        if not isinstance(item, dict) or not {"path", "sha256", "size_bytes"} <= set(item):
            raise RefitDagError(f"{label} artifact member schema differs")
        member = str(item["path"])
        if member in members or Path(member).is_absolute() or ".." in Path(member).parts:
            raise RefitDagError(f"{label} artifact member path differs")
        members[member] = item
    return tree, manifest, members


def bound_member(
    tree: Path,
    members: Mapping[str, Mapping[str, Any]],
    member: object,
    *,
    label: str,
) -> Path:
    if not isinstance(member, str) or member not in members:
        raise RefitDagError(f"{label} is not bound by ARTIFACTS")
    path = tree / member
    authority = members[member]
    if (
        path.is_symlink()
        or not path.is_file()
        or path.stat().st_size != int(authority["size_bytes"])
        or file_sha256(path) != authority["sha256"]
    ):
        raise RefitDagError(f"{label} member bytes differ")
    return path


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != STATUS
        or config.get("dataset_id") != "gse281364"
        or config.get("task_id") != "variant_to_regulation"
        or config.get("endpoint_id") != "mpra_allelic_direction"
        or config.get("planning_authorized") is not True
        or any(
            config.get(field) is not False
            for field in (
                "outcome_access_authorized",
                "prediction_generation_authorized",
                "model_fit_authorized",
                "stack_fit_authorized",
                "residual_correlation_authorized",
                "conditional_model_build_authorized",
                "sealed_asset_access_authorized",
            )
        )
    ):
        raise RefitDagError("top-level outcome firewall differs")
    row = config.get("row_binding", {})
    if (
        row.get("source_locus_groups") != 1033
        or row.get("long_range_blocks") != 239
        or row.get("base_rows_per_seed") != 2066
        or row.get("seeded_rows") != 10330
        or tuple(row.get("fixed_seeds", ())) != SEEDS
        or tuple(row.get("assay_contexts", ())) != CONTEXTS
        or row.get("largest_receptive_field_buffer_bp") != 524288
        or row.get("original_6kb_predictions_reusable") is not False
    ):
        raise RefitDagError("row or long-range-fold binding differs")
    schema = config.get("common_prediction_schema", {})
    if (
        tuple(schema.get("fields", ())) != COMPONENT_FIELDS
        or tuple(schema.get("fixed_seeds", ())) != SEEDS
        or schema.get("required_study_id") != "gse281364"
        or schema.get("required_stratum") != "all"
        or schema.get("allele_sign") != "ALT_minus_REF"
        or schema.get("missing_allowed") is not False
        or schema.get("missing_as_zero_allowed") is not False
        or not all(
            schema.get(field) is True
            for field in (
                "same_rows_all_models",
                "same_rows_all_seeds",
                "same_outer_fold_all_models",
            )
        )
    ):
        raise RefitDagError("common prediction schema differs")
    eligibility = config.get("eligibility", {})
    if (
        tuple(eligibility.get("allowed_checkpoint_exposure", ()))
        != ("clean_declared", "target_label_unexposed")
        or eligibility.get("restricted_or_unknown_may_enter_open_trigger") is not False
        or eligibility.get("cross_cohort_hepatocyte_context_is_condition_matched") is not False
        or eligibility.get("conditional_architecture_may_supply_its_own_trigger") is not False
        or tuple(eligibility.get("current_open_sequence_paths_after_refit", ()))
        != ("caduceus", "dnabert2", "hyenadna")
        or eligibility.get("preferred_open_sequence_path") != "caduceus"
        or eligibility.get("current_open_context_paths") != []
        or eligibility.get("current_trigger_compatible_pair_count") != 0
    ):
        raise RefitDagError("license, exposure, or trigger eligibility differs")
    context = config.get("context_transfer_contract", {})
    if (
        context.get("source_dataset") != "gse296875"
        or context.get("context_key") != "training_lineage_mean_released_rank"
        or context.get("lineage") != "hepatocyte"
        or context.get("lineage_index") != 2
        or context.get("contexts") != 5
        or context.get("context_width") != 2891
        or context.get("query_dataset") != "gse281364"
        or context.get("condition_matched") is not False
    ):
        raise RefitDagError("cross-cohort Corgi context semantics differ")


def validate_authorities(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    file_authorities = config.get("file_authorities", {})
    if set(file_authorities) != {
        "corgi_checkpoint",
        "corgi_exposure",
        "dataset_manifest",
        "dna_language_registry",
        "mandatory_baseline_registry",
        "regulatory_local_registry",
        "regulatory_sequence_registry",
        "stacking_implementation",
    }:
        raise RefitDagError("file authority roster differs")
    files = {}
    for authority_id, authority in file_authorities.items():
        path = safe_path(root, authority.get("path"), label=authority_id)
        if path.is_symlink() or not path.is_file() or file_sha256(path) != authority.get("sha256"):
            raise RefitDagError(f"{authority_id} authority differs")
        files[authority_id] = path

    row = config["row_binding"]
    row_tree, _, row_members = manifest_authority(
        root, row["tree_path"], row["artifacts_sha256"], label="row universe"
    )
    row_receipt_path = bound_member(
        row_tree, row_members, row["receipt_member"], label="row-universe receipt"
    )
    row_tsv_path = bound_member(
        row_tree, row_members, row["row_member"], label="row-universe table"
    )
    row_receipt = load_json(row_receipt_path, label="row-universe receipt")
    if (
        file_sha256(row_tsv_path) != row["row_tsv_sha256"]
        or row_receipt.get("base_row_set_sha256") != row["base_row_set_sha256"]
        or row_receipt.get("source_locus_groups") != 1033
        or row_receipt.get("long_range_blocks") != 239
        or row_receipt.get("seeded_rows") != 10330
        or row_receipt.get("source_groups_reassigned_from_original_6kb_fold_map") != 830
        or row_receipt.get("outcomes_read")
        or row_receipt.get("prediction_values_read")
        or row_receipt.get("stack_fit")
        or row_receipt.get("conditional_model_built_or_fit")
    ):
        raise RefitDagError("row-universe receipt differs")

    outcome = config["outcome_authority"]
    outcome_tree, outcome_manifest, outcome_members = manifest_authority(
        root,
        outcome["tree_path"],
        outcome["artifacts_sha256"],
        label="development outcome authority",
    )
    outcome_member = outcome_members.get(outcome["member"])
    if (
        outcome_member is None
        or outcome_member.get("sha256") != outcome["member_sha256"]
        or outcome_manifest.get("metadata", {}).get("artifact_class")
        != "gse281364_replicate_safe_outcomes"
        or outcome_manifest.get("metadata", {}).get("donor_count") != 0
        or outcome.get("read_by_reconciliation") is not False
    ):
        raise RefitDagError("development outcome authority differs")
    # Deliberately do not open or hash the outcome member here.

    reference = config["reference_authority"]
    manifest_authority(
        root,
        reference["tree_path"],
        reference["artifacts_sha256"],
        label="project reference",
    )

    fixtures = {}
    for fixture in config.get("projection_fixture", ()):
        tree, _, members = manifest_authority(
            root,
            fixture["tree_path"],
            fixture["artifacts_sha256"],
            label=f"projection fixture {fixture['fixture_id']}",
        )
        bound_member(
            tree,
            members,
            fixture["manifest_member"],
            label=f"projection fixture {fixture['fixture_id']} manifest",
        )
        fixtures[fixture["fixture_id"]] = tree
    if set(fixtures) != {"caduceus_131072", "dna_lm_common_4096"}:
        raise RefitDagError("projection fixture roster differs")

    return {
        "files": files,
        "row_tree": row_tree,
        "row_receipt": row_receipt,
        "row_tsv": row_tsv_path,
        "outcome_tree": outcome_tree,
        "fixtures": fixtures,
    }


def validate_registry_states(
    authorities: Mapping[str, Any]
) -> dict[str, str]:
    expected = {
        "borzoi_ensemble": ("blocked_terms", "target_label_unexposed"),
        "caduceus": ("candidate", "target_label_unexposed"),
        "corgi_regular": ("candidate", "target_label_unexposed"),
        "dnabert2": ("candidate", "target_label_unexposed"),
        "enformer": ("restricted_comparator", "target_label_unexposed"),
        "hyenadna": ("candidate", "target_label_unexposed"),
        "mpralegnet": ("candidate", "target_label_unexposed"),
        "nucleotide_transformer": (
            "restricted_comparator",
            "target_label_unexposed",
        ),
        "sei": ("restricted_comparator", "target_label_unexposed"),
        "sequence_cnn_control": ("candidate", "unknown"),
        "sequence_transformer_control": ("candidate", "unknown"),
    }
    records: dict[str, Mapping[str, Any]] = {}
    for authority_id in (
        "dna_language_registry",
        "regulatory_sequence_registry",
        "regulatory_local_registry",
        "mandatory_baseline_registry",
    ):
        path = authorities["files"][authority_id]
        try:
            source = tomllib.loads(path.read_text(encoding="utf-8"))
        except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
            raise RefitDagError(f"cannot parse {authority_id}: {error}") from error
        for row in source.get("models", ()):
            model_id = row.get("model_id")
            if model_id in expected:
                if model_id in records:
                    raise RefitDagError(f"duplicate registry model {model_id}")
                records[model_id] = row
    if set(records) != set(expected):
        raise RefitDagError("registry model census differs")
    states = {}
    for model_id, (status, exposure) in expected.items():
        row = records[model_id]
        if row.get("status") != status or row.get("exposure_status") != exposure:
            raise RefitDagError(f"registry state differs for {model_id}")
        states[model_id] = f"{status}|{exposure}"
    return states


def validate_inventory(root: Path, config: Mapping[str, Any]) -> list[dict[str, Any]]:
    inventory = config.get("artifact_inventory", ())
    if len(inventory) != 12 or {item.get("model_id") for item in inventory} != MODEL_IDS:
        raise RefitDagError("model inventory roster differs")
    seen: set[str] = set()
    output = []
    for item in inventory:
        if set(item) != set(INVENTORY_FIELDS) or item["artifact_id"] in seen:
            raise RefitDagError("artifact inventory schema or identity differs")
        seen.add(item["artifact_id"])
        tree, _, members = manifest_authority(
            root,
            item["tree_path"],
            item["artifacts_sha256"],
            label=f"inventory {item['artifact_id']}",
        )
        receipt_path = bound_member(
            tree,
            members,
            item["receipt_member"],
            label=f"inventory {item['artifact_id']} receipt",
        )
        receipt = load_json(receipt_path, label=f"inventory {item['artifact_id']} receipt")
        model = item["model_id"]
        if model in {"dnabert2", "hyenadna", "nucleotide_transformer"}:
            if (
                receipt.get("status") != "pass_outcome_blind_embedding_extraction"
                or receipt.get("elements") != 1033
                or receipt.get("downstream_head_fit")
                or receipt.get("reporter_outcomes_read")
            ):
                raise RefitDagError(f"{model} raw receipt differs")
        elif model == "caduceus":
            if (
                receipt.get("status") != "pass_outcome_blind_embedding_extraction"
                or receipt.get("elements") != 1033
                or receipt.get("downstream_head_fit")
                or receipt.get("reporter_outcomes_read")
            ):
                raise RefitDagError("Caduceus raw receipt differs")
        elif model == "mpralegnet" and (
            receipt.get("elements") != 4359 or receipt.get("outcomes_read")
        ):
            raise RefitDagError("MPRALegNet denominator differs")
        elif model == "sei" and (
            receipt.get("elements") != 1033 or receipt.get("outcomes_read")
        ):
            raise RefitDagError("Sei static-score receipt differs")
        elif model == "enformer_crested_restricted_port" and (
            receipt.get("elements") != 1033 or receipt.get("outcomes_read")
        ):
            raise RefitDagError("Enformer static-score receipt differs")
        elif model == "borzoi_ensemble" and (
            receipt.get("model_forward_executed") or receipt.get("outcomes_read")
        ):
            raise RefitDagError("Borzoi fixture-only state differs")
        elif model == "corgi_regular":
            checkpoints = receipt.get("checkpoints", ())
            regular = [row for row in checkpoints if row.get("filename") == "corgi_model.pt"]
            if (
                len(regular) != 1
                or regular[0].get("sha256")
                != "cd51539f01de1e66a3a4f9b89e772bdeae07df2f0178d5b692f2368db2f4d2a4"
            ):
                raise RefitDagError("Corgi checkpoint identity differs")
        output.append(dict(item))
    return output


def validate_context_mappers(
    root: Path, config: Mapping[str, Any]
) -> list[dict[str, Any]]:
    mappers = config.get("context_mapper", ())
    if len(mappers) != 5 or {item.get("mapper_fold") for item in mappers} != set(range(5)):
        raise RefitDagError("Corgi mapper roster differs")
    output = []
    for item in sorted(mappers, key=lambda value: value["mapper_fold"]):
        tree, _, members = manifest_authority(
            root,
            item["tree_path"],
            item["artifacts_sha256"],
            label=f"Corgi mapper fold {item['mapper_fold']}",
        )
        receipt_path = bound_member(
            tree,
            members,
            item["receipt_member"],
            label=f"Corgi mapper fold {item['mapper_fold']} receipt",
        )
        context_path = bound_member(
            tree,
            members,
            item["context_member"],
            label=f"Corgi mapper fold {item['mapper_fold']} contexts",
        )
        receipt = load_json(receipt_path, label="Corgi mapper receipt")
        if (
            receipt.get("schema_version") != "masld-bench-corgi-masked-context-mapper-v1"
            or receipt.get("status") != "pass_outcome_free_mapper"
            or receipt.get("outer_fold") != item["mapper_fold"]
            or receipt.get("held_ATAC_or_other_outcomes_used")
            or receipt.get("test_or_sealed_outcomes_read")
            or "lineage_mean_contexts" not in receipt.get("training_only_objects", ())
        ):
            raise RefitDagError("Corgi mapper receipt differs")
        output.append({**dict(item), "tree": tree, "context_path": context_path})
    return output


def validate_stages(root: Path, config: Mapping[str, Any]) -> list[dict[str, Any]]:
    stages = config.get("stages", ())
    if len(stages) != 12:
        raise RefitDagError("refit stage count differs")
    identifiers = {stage.get("stage_id") for stage in stages}
    if len(identifiers) != len(stages):
        raise RefitDagError("duplicate refit stage identity")
    ordered = sorted(stages, key=lambda stage: stage.get("order", -1))
    if [stage.get("order") for stage in ordered] != list(range(12)):
        raise RefitDagError("refit stage order differs")
    prior: set[str] = set()
    for stage in ordered:
        if set(stage) != set(STAGE_FIELDS):
            raise RefitDagError(f"stage {stage.get('stage_id')} schema differs")
        dependencies = set(stage["dependencies"])
        if not dependencies <= prior:
            raise RefitDagError(f"stage {stage['stage_id']} has a forward dependency")
        prior.add(stage["stage_id"])
        if stage["state"] in READY_STATES:
            entrypoint = safe_path(root, stage["entrypoint"], label=stage["stage_id"])
            if not entrypoint.is_file() or stage["command"] == "BLOCKED" or stage["blockers"]:
                raise RefitDagError(f"ready stage {stage['stage_id']} is not executable")
        else:
            if stage["command"] != "BLOCKED" or not stage["blockers"]:
                raise RefitDagError(f"blocked stage {stage['stage_id']} is not fail-closed")
        if stage["outcomes_required"] and stage["state"] in READY_STATES:
            raise RefitDagError("outcome-requiring stage is prematurely ready")
    if any(
        token in stage["stage_id"]
        for stage in ordered
        for token in ("residual", "stack_fit", "conditional_model", "novel_model")
    ):
        raise RefitDagError("prohibited modeling stage entered refit DAG")
    return [dict(stage) for stage in ordered]


def tsv_text(fields: Sequence[str], rows: Sequence[Mapping[str, object]]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow(
            {
                field: ";".join(str(value) for value in row[field])
                if isinstance(row[field], list)
                else str(row[field]).lower()
                if isinstance(row[field], bool)
                else row[field]
                for field in fields
            }
        )
    return output.getvalue()


def materialize_contexts(
    config: Mapping[str, Any], mappers: Sequence[Mapping[str, Any]], output: Path
) -> dict[str, Any]:
    import numpy as np

    if output.exists():
        raise RefitDagError("Corgi context-roster output exists")
    key = config["context_transfer_contract"]["context_key"]
    lineage_index = int(config["context_transfer_contract"]["lineage_index"])
    values = []
    roster = []
    for item in mappers:
        with np.load(item["context_path"], allow_pickle=False) as source:
            if key not in source.files:
                raise RefitDagError("Corgi mapper lacks training lineage means")
            means = source[key]
        if means.shape != (5, 2891) or not np.isfinite(means).all():
            raise RefitDagError("Corgi training lineage-mean geometry differs")
        values.append(np.asarray(means[lineage_index], dtype=np.float32))
        roster.append(
            {
                "context_index": len(roster),
                "context_id": f"gse296875_mapper_fold{item['mapper_fold']}_training_hepatocyte_mean",
                "mapper_fold": item["mapper_fold"],
                "source_dataset": "gse296875",
                "lineage": "hepatocyte",
                "query_dataset": "gse281364",
                "condition_matched": "false",
                "claim_role": "cross_cohort_MASLD_hepatocyte_context_transfer_diagnostic_only",
            }
        )
    matrix = np.vstack(values)
    output.mkdir(parents=True, mode=0o750)
    np.savez_compressed(
        output / "contexts.npz",
        context_ids=np.asarray([row["context_id"] for row in roster]),
        mapper_folds=np.asarray([row["mapper_fold"] for row in roster], dtype=np.int64),
        context_values=matrix,
    )
    roster_fields = tuple(roster[0])
    write_text_exclusive(output / "context_roster.tsv", tsv_text(roster_fields, roster))
    receipt = {
        "schema_version": "masld-bench-gse281364-corgi-cross-cohort-context-roster-v1",
        "status": "pass_outcome_blind_restricted_context_roster",
        "source_dataset": "gse296875",
        "query_dataset": "gse281364",
        "contexts": 5,
        "context_width": 2891,
        "lineage": "hepatocyte",
        "context_key": key,
        "condition_matched": False,
        "claim_role": "cross_cohort_MASLD_hepatocyte_context_transfer_diagnostic_only",
        "outcomes_read": False,
        "prediction_values_read": False,
        "model_forward_executed": False,
        "model_fit": False,
        "stack_fit": False,
        "conditional_model_built_or_fit": False,
    }
    write_json_exclusive(output / "receipt.json", receipt)
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def prepare(root: Path, config_path: Path, output: Path) -> dict[str, Any]:
    if output.exists():
        raise RefitDagError("refit-DAG output exists")
    config = load_config(config_path)
    validate_config(config)
    authorities = validate_authorities(root, config)
    registry_states = validate_registry_states(authorities)
    inventory = validate_inventory(root, config)
    mappers = validate_context_mappers(root, config)
    stages = validate_stages(root, config)
    output.mkdir(parents=True, mode=0o750)
    write_text_exclusive(output / "artifact_inventory.tsv", tsv_text(INVENTORY_FIELDS, inventory))
    write_text_exclusive(output / "refit_stages.tsv", tsv_text(STAGE_FIELDS, stages))
    schema_rows = [
        {
            "position": index,
            "field": field,
            "nullable": "false",
            "rule": {
                "seed": "one_of_1103_2909_4721_6673_8111",
                "row_hash": "exact_frozen_row_universe_identity",
                "unit_hash": "source_locus_group_hash",
                "block_hash": "524288bp_long_range_component_hash",
                "stratum": "all",
                "outer_fold": "frozen_long_range_fold",
                "study_id": "gse281364",
                "observed": "development_MPRA_mean_signed_ALT_minus_REF_activity",
                "sequence_prediction": "outer_block_out_of_fold_assay_native_unit",
                "context_prediction": "outer_block_out_of_fold_assay_native_unit",
            }[field],
        }
        for index, field in enumerate(COMPONENT_FIELDS)
    ]
    write_text_exclusive(
        output / "common_prediction_schema.tsv",
        tsv_text(("position", "field", "nullable", "rule"), schema_rows),
    )
    dag = {
        "schema_version": SCHEMA,
        "status": STATUS,
        "nodes": stages,
        "edges": [
            {"from": dependency, "to": stage["stage_id"]}
            for stage in stages
            for dependency in stage["dependencies"]
        ],
        "open_sequence_priority": ["caduceus", "dnabert2", "hyenadna"],
        "restricted_sequence_diagnostics": [
            "nucleotide_transformer",
            "sei",
            "enformer_crested_restricted_port",
            "mpralegnet",
        ],
        "context_priority": "corgi_regular_cross_cohort_gse296875_training_hepatocyte_mean",
        "context_condition_matched": False,
        "current_trigger_compatible_pair_count": 0,
    }
    write_json_exclusive(output / "refit_dag.json", dag)
    receipt = {
        "schema_version": SCHEMA,
        "status": STATUS,
        "dataset_id": "gse281364",
        "task_id": "variant_to_regulation",
        "inventory_rows": len(inventory),
        "stage_rows": len(stages),
        "ready_or_complete_stages": sum(stage["state"] in READY_STATES for stage in stages),
        "blocked_stages": sum(stage["state"] not in READY_STATES for stage in stages),
        "open_sequence_paths_after_refit": ["caduceus", "dnabert2", "hyenadna"],
        "open_context_paths": [],
        "restricted_context_roster_ready": True,
        "current_trigger_compatible_pair_count": 0,
        "row_universe_artifacts_sha256": config["row_binding"]["artifacts_sha256"],
        "row_universe_tsv_sha256": config["row_binding"]["row_tsv_sha256"],
        "base_row_set_sha256": authorities["row_receipt"]["base_row_set_sha256"],
        "fixed_seeds": list(SEEDS),
        "base_rows_per_seed": 2066,
        "seeded_rows": 10330,
        "long_range_blocks": 239,
        "source_groups_reassigned_from_original_6kb_fold_map": 830,
        "outcome_manifest_bound": True,
        "outcome_member_opened": False,
        "outcomes_read": False,
        "prediction_artifact_payloads_parsed": False,
        "prediction_values_read": False,
        "predictions_generated": False,
        "model_fit": False,
        "stack_fit": False,
        "residual_correlation_calculated": False,
        "conditional_model_built_or_fit": False,
        "sealed_assets_read": False,
        "context_condition_matched": False,
        "context_claim_role": "cross_cohort_MASLD_hepatocyte_context_transfer_diagnostic_only",
        "registry_states": registry_states,
    }
    write_json_exclusive(output / "receipt.json", receipt)
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--materialize-corgi-contexts", action="store_true")
    arguments = parser.parse_args()
    root = arguments.root.resolve(strict=True)
    config = load_config(arguments.config)
    validate_config(config)
    if arguments.materialize_corgi_contexts:
        validate_authorities(root, config)
        mappers = validate_context_mappers(root, config)
        materialize_contexts(config, mappers, arguments.output)
    else:
        prepare(root, arguments.config, arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
