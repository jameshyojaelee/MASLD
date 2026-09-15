#!/usr/bin/env python3
"""Audit whether current development predictions can support the stack trigger."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import tomllib
from typing import Any, Iterable, Mapping, Sequence

from masld_bench.artifacts import (
    ArtifactError,
    reject_symlink_components,
    verify_frozen_tree,
    write_json_exclusive,
    write_text_exclusive,
)


SCHEMA = "masld-bench-complementarity-readiness-audit-v1"
STATUS = "blocked_no_trigger_compatible_prediction_pair"
TASKS = ("variant_to_regulation", "rna_conditioned_atac")
SEEDS = (1103, 2909, 4721, 6673, 8111)
LANES = (
    "gse296875_local_sequence_profile",
    "gse296875_corgi_context_regional",
    "gse296875_midas_inductive_context",
    "gse296875_scbasset_sequence",
    "gse281364_dna_lm_common_mpra",
    "gse281364_native_sequence_transfer",
    "gse281364_borzoi_native",
    "gse244832_rna_atac_development",
)
PAIRS = (
    "local_sequence_family_internal",
    "local_sequence_vs_corgi",
    "scbasset_vs_midas",
    "dna_lm_common_internal",
    "native_mpra_cross_family",
    "any_current_variant_vs_borzoi",
)
SOURCE_FIELDS = (
    "variant_development_shortlist",
    "rna_atac_development_shortlist",
    "sequence_development_residual_bundle",
    "context_development_residual_bundle",
    "variant_development_stack_gain_bundle",
    "rna_atac_development_stack_gain_bundle",
)
FILE_AUTHORITIES = {
    "stacking_source",
    "context_preflight",
    "context_checkpoints",
    "rna_atac_task",
    "variant_task",
    "error_audit_campaign",
    "conditional_campaign",
    "family_native_tournament",
    "local_variant_config",
}
TREE_AUTHORITIES = {
    "local_variant_reconciliation",
    "local_profile_incomplete_matrix",
    "corgi_hepatocyte_tile_predictions",
    "midas_inductive_predictions",
    "scbasset_fold0_predictions",
    "scbasset_fold1to4_prediction_lock",
    "dna_lm_family_reconciliation",
    "native_mpra_evaluation",
    "enformer_mpra_evaluation",
    "borzoi_native_execution_preflight",
}
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class ComplementarityReadinessError(RuntimeError):
    """Raised when a source, row identity, or trigger boundary differs."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ComplementarityReadinessError(f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise ComplementarityReadinessError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise ComplementarityReadinessError(f"{label} must be a JSON object")
    return value


def load_toml(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise ComplementarityReadinessError(f"{label} is not a regular file")
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, tomllib.TOMLDecodeError) as error:
        raise ComplementarityReadinessError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise ComplementarityReadinessError(f"{label} must be a TOML table")
    return value


def safe_path(root: Path, configured: str, *, label: str) -> Path:
    value = Path(configured)
    if not value.is_absolute() and (not value.parts or ".." in value.parts):
        raise ComplementarityReadinessError(f"unsafe {label}: {configured}")
    lexical = value if value.is_absolute() else root / value
    try:
        path = reject_symlink_components(lexical, label=label).resolve(strict=True)
        path.relative_to(root)
    except (ArtifactError, OSError, ValueError) as error:
        raise ComplementarityReadinessError(
            f"{label} is missing or escapes project root: {configured}"
        ) from error
    return path


def manifest_metadata(path: Path, *, label: str) -> Mapping[str, Any]:
    manifest = load_json(path / "ARTIFACTS.json", label=f"{label} manifest")
    metadata = manifest.get("metadata", {})
    if not isinstance(metadata, Mapping):
        raise ComplementarityReadinessError(f"{label} metadata differs")
    return metadata


def validate_config(config: Mapping[str, Any]) -> None:
    lanes = config.get("prediction_lanes", ())
    pairs = config.get("pair_audits", ())
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != STATUS
        or config.get("audit_id")
        != "current_development_prediction_complementarity_readiness"
        or config.get("build_authorized") is not False
        or config.get("stack_fit_or_residual_correlation_authorized") is not False
        or config.get("sealed_assets_read") is not False
        or config.get("raw_or_summary_outcomes_read") is not False
        or config.get("prediction_values_read") is not False
        or config.get("row_identifiers_may_be_read") is not True
        or tuple(row.get("lane_id") for row in lanes) != LANES
        or tuple(row.get("pair_id") for row in pairs) != PAIRS
    ):
        raise ComplementarityReadinessError("audit identity or roster differs")

    trigger = config.get("trigger_contract", {})
    if (
        tuple(trigger.get("required_task_ids", ())) != TASKS
        or tuple(trigger.get("fixed_seeds", ())) != SEEDS
        or trigger.get("variant_endpoint_id") != "mpra_allelic_direction"
        or trigger.get("variant_study_ids") != ["gse281364"]
        or trigger.get("variant_minimum_gain") != 0.02
        or trigger.get("rna_atac_endpoint_id") != "rna_atac_regulatory_profile"
        or trigger.get("rna_atac_study_ids") != ["gse244832", "gse296875"]
        or trigger.get("rna_atac_minimum_gain") != 0.05
        or trigger.get("minimum_qualifying_seed_count") != 4
        or trigger.get("residual_correlation_strictly_below") != 0.8
        or trigger.get("minimum_study_count") != 2
        or trigger.get("cross_fitted") is not True
        or trigger.get("nonnegative_stack") is not True
        or trigger.get("best_open_models_compared") is not True
        or trigger.get("distinct_sequence_and_context_source_families") is not True
        or trigger.get("required_source_counts")
        != {
            "development_shortlist": 2,
            "development_residual_bundle": 2,
            "development_stack_gain_bundle": 2,
        }
    ):
        raise ComplementarityReadinessError("trigger contract differs")

    sources = config.get("required_stacking_sources", {})
    if set(sources) != set(SOURCE_FIELDS) or any(
        sources[field] != "UNBOUND" for field in SOURCE_FIELDS
    ):
        raise ComplementarityReadinessError("stacking sources must remain unbound")

    for row in lanes:
        if (
            row.get("residual_ready") is not False
            or row.get("stack_gain_ready") is not False
            or row.get("open_component_ready") is not False
            or not row.get("blocking_findings")
        ):
            raise ComplementarityReadinessError(
                f"prediction lane readiness differs: {row.get('lane_id')}"
            )
    if any(row.get("trigger_ready") is not False for row in pairs):
        raise ComplementarityReadinessError("a pair cannot be trigger-ready")
    if (
        next(row for row in lanes if row["lane_id"] == "gse244832_rna_atac_development")[
            "authority_ids"
        ]
        != []
        or next(row for row in lanes if row["lane_id"] == "gse281364_borzoi_native")[
            "prediction_schema"
        ]
        != "none_fixture_only"
    ):
        raise ComplementarityReadinessError("absent prediction authority differs")

    ask = config.get("bundled_cpu_validation_ask", {})
    if (
        ask.get("job_name") != "model-check-219"
        or ask.get("partition") != "cpu"
        or ask.get("qos") != "nslab"
        or ask.get("cpus_per_task") != 2
        or ask.get("memory_gb") != 24
        or ask.get("walltime") != "04:00:00"
        or ask.get("array") is not False
        or ask.get("job_count") != 1
        or ask.get("total_cpu_hours") != 8
        or ask.get("total_gpu_hours") != 0
        or ask.get("submitted") is not False
    ):
        raise ComplementarityReadinessError("CPU validation ask differs")


def validate_authorities(
    root: Path, config: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    files = config.get("frozen_file_authorities", {})
    trees = config.get("frozen_tree_authorities", {})
    if set(files) != FILE_AUTHORITIES or set(trees) != TREE_AUTHORITIES:
        raise ComplementarityReadinessError("authority census differs")
    result: dict[str, dict[str, Any]] = {}
    for name, binding in files.items():
        path = safe_path(root, str(binding.get("path", "")), label=f"file {name}")
        expected = str(binding.get("sha256", ""))
        if not SHA256.fullmatch(expected) or sha256_file(path) != expected:
            raise ComplementarityReadinessError(f"frozen file changed: {name}")
        result[name] = {"kind": "file", "path": path.relative_to(root).as_posix(), "sha256": expected}
    for name, binding in trees.items():
        path = safe_path(root, str(binding.get("path", "")), label=f"tree {name}")
        expected = str(binding.get("artifacts_sha256", ""))
        if not SHA256.fullmatch(expected):
            raise ComplementarityReadinessError(f"tree hash differs: {name}")
        try:
            verify_frozen_tree(path)
        except ArtifactError as error:
            raise ComplementarityReadinessError(
                f"frozen tree verification failed: {name}: {error}"
            ) from error
        if sha256_file(path / "ARTIFACTS.json") != expected:
            raise ComplementarityReadinessError(f"frozen tree changed: {name}")
        result[name] = {"kind": "tree", "path": path.relative_to(root).as_posix(), "artifacts_sha256": expected}
    return result


def validate_trigger_authorities(root: Path, config: Mapping[str, Any]) -> None:
    files = config["frozen_file_authorities"]
    rna_task = load_toml(root / files["rna_atac_task"]["path"], label="RNA-ATAC TaskSpec")
    variant_task = load_toml(root / files["variant_task"]["path"], label="variant TaskSpec")
    error_campaign = load_toml(
        root / files["error_audit_campaign"]["path"], label="error-audit campaign"
    )
    conditional_campaign = load_toml(
        root / files["conditional_campaign"]["path"], label="conditional campaign"
    )
    context = load_json(root / files["context_preflight"]["path"], label="context preflight")
    gate = context.get("complementarity_gate", {})
    if (
        rna_task.get("task_id") != "rna_conditioned_atac"
        or rna_task.get("datasets_train") != ["gse296875"]
        or rna_task.get("datasets_development") != ["gse244832"]
        or rna_task.get("primary_evaluator_id")
        != "rna_atac_two_way_deviance_reduction_v1"
        or variant_task.get("task_id") != "variant_to_regulation"
        or variant_task.get("unit_of_inference") != "independent_locus"
        or variant_task.get("datasets_train")[-2:] != ["gse281364", "gse296875"]
        or error_campaign.get("status") != "blocked_prerequisites"
        or error_campaign.get("prediction_bundle_set_sha256") != "UNRESOLVED"
        or conditional_campaign.get("status") != "blocked_trigger"
        or conditional_campaign.get("conditional_model_spec_path") != "UNRESOLVED"
        or gate.get("state") != "NOT_EVALUATED"
        or gate.get("build_authorized") is not False
        or gate.get("source_binding_counts")
        != {
            "development_shortlist": 2,
            "development_residual_bundle": 2,
            "development_stack_gain_bundle": 2,
        }
        or gate.get("evaluated_seed_count") != 5
        or gate.get("minimum_qualifying_seed_count") != 4
        or gate.get("minimum_study_count") != 2
    ):
        raise ComplementarityReadinessError("campaign or trigger authority differs")


def validate_prediction_metadata(root: Path, config: Mapping[str, Any]) -> dict[str, Any]:
    trees = config["frozen_tree_authorities"]
    local_root = root / trees["local_variant_reconciliation"]["path"]
    local = load_json(local_root / "reconciliation/receipt.json", label="local family receipt")
    matrix_root = root / trees["local_profile_incomplete_matrix"]["path"]
    matrix = load_json(matrix_root / "matrix_summary.json", label="local profile matrix")
    corgi_root = root / trees["corgi_hepatocyte_tile_predictions"]["path"]
    corgi_meta = manifest_metadata(corgi_root, label="Corgi tile predictions")
    corgi_folds = [load_json(corgi_root / f"fold{fold}/receipt.json", label=f"Corgi fold {fold}") for fold in range(5)]
    midas_root = root / trees["midas_inductive_predictions"]["path"]
    midas = load_json(midas_root / "predictions/receipt.json", label="MIDAS receipt")
    scb_lock_root = root / trees["scbasset_fold1to4_prediction_lock"]["path"]
    scb_lock = load_json(scb_lock_root / "prediction_lock.json", label="scBasset prediction lock")
    scb0_root = root / trees["scbasset_fold0_predictions"]["path"]
    scb0_meta = manifest_metadata(scb0_root, label="scBasset fold0 predictions")
    dna_meta = manifest_metadata(
        root / trees["dna_lm_family_reconciliation"]["path"], label="DNA LM reconciliation"
    )
    native_meta = manifest_metadata(
        root / trees["native_mpra_evaluation"]["path"], label="native MPRA evaluation"
    )
    enformer_meta = manifest_metadata(
        root / trees["enformer_mpra_evaluation"]["path"], label="Enformer MPRA evaluation"
    )
    borzoi = load_json(
        root / trees["borzoi_native_execution_preflight"]["path"] / "receipt.json",
        label="Borzoi execution preflight",
    )
    if (
        local.get("status") != "pass_reconciled_native_tasks_no_model_promoted"
        or local.get("active_campaign_outputs_bound_or_read") is not False
        or local.get("additional_development_screens_allowed_now") != []
        or matrix.get("evaluation_groups") != 7
        or matrix.get("cross_fold_candidate_ranking_calculated") is not False
        or corgi_meta.get("lineage_id") != "hepatocyte"
        or corgi_meta.get("logical_tasks") != 5
        or corgi_meta.get("outcomes_read") is not False
        or any(row.get("held_ATAC_or_other_outcomes_used") is not False for row in corgi_folds)
        or any(row.get("open_champion_eligible") is not False for row in corgi_folds)
        or any(row.get("native_numeric_parity_established") is not False for row in corgi_folds)
        or midas.get("fold_count") != 5
        or midas.get("base_seed") != 2711
        or midas.get("held_atac_input_exposed") is not False
        or midas.get("outcomes_read") is not False
        or scb_lock.get("seed") != 20260824
        or scb_lock.get("folds") != [1, 2, 3, 4]
        or scb_lock.get("outcomes_read") is not False
        or scb0_meta.get("seed") != 11
        or scb0_meta.get("evaluation_role") != "valid"
        or dna_meta.get("family_promotion") != "none"
        or dna_meta.get("model_count") != 4
        or native_meta.get("models_forced_to_common_endpoint") is not False
        or native_meta.get("models_ranked") is not False
        or enformer_meta.get("models_ranked") is not False
        or enformer_meta.get("restricted_comparator") is not True
        or borzoi.get("executable") is not False
        or borzoi.get("model_forward_executed") is not False
        or borzoi.get("reporter_outcomes_read") is not False
        or borzoi.get("sealed_assets_read") is not False
    ):
        raise ComplementarityReadinessError("prediction metadata disposition differs")
    return {
        "local_profile_evaluation_groups": matrix["evaluation_groups"],
        "corgi_lineages": [corgi_meta["lineage_id"]],
        "corgi_outer_folds": len(corgi_folds),
        "midas_seed_ids": [midas["base_seed"]],
        "scbasset_seed_ids": [scb0_meta["seed"], scb_lock["seed"]],
        "dna_lm_family_promotion": dna_meta["family_promotion"],
        "borzoi_prediction_available": False,
        "gse244832_prediction_authority_bound": False,
    }


def _row_id_records(path: Path) -> Iterable[tuple[str, str]]:
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if tuple(reader.fieldnames or ()) != ("row_hash", "donor_hash"):
                raise ComplementarityReadinessError(
                    f"row identifier table has wrong header: {path}"
                )
            for index, row in enumerate(reader):
                row_hash = row["row_hash"]
                donor_hash = row["donor_hash"]
                if not SHA256.fullmatch(row_hash) or not SHA256.fullmatch(donor_hash):
                    raise ComplementarityReadinessError(
                        f"invalid row identifier at {path}:{index + 2}"
                    )
                yield row_hash, donor_hash
    except (OSError, UnicodeDecodeError, csv.Error) as error:
        raise ComplementarityReadinessError(
            f"cannot read row identifiers {path}: {error}"
        ) from error


def _bundle_row_ids(bundle_root: Path, bundle: Mapping[str, Any]) -> Path:
    row_ids = bundle.get("row_ids", {})
    if (
        bundle.get("schema_version") != "masld-bench-prediction-bundle-v1"
        or bundle.get("task_id") != "rna_conditioned_atac"
        or bundle.get("biological_unit") != "donor"
        or bundle.get("row_id_field") != "row_hash"
        or bundle.get("unit_id_field") != "donor_hash"
        or bundle.get("table_schema_sha256")
        != "a622b95c2aec79a86a448179960301815f31f0ced8735d3b11e5c2294dff5c01"
        or not isinstance(row_ids, Mapping)
        or row_ids.get("role") != "prediction_row_ids:rna_conditioned_atac"
    ):
        raise ComplementarityReadinessError("PredictionBundle row contract differs")
    path = safe_path(bundle_root, str(row_ids.get("path", "")), label="row identifiers")
    expected = str(row_ids.get("sha256", ""))
    if not SHA256.fullmatch(expected) or sha256_file(path) != expected:
        raise ComplementarityReadinessError("row identifier artifact changed")
    return path


def _scbasset_bundles(root: Path, config: Mapping[str, Any]) -> dict[int, tuple[Path, dict[str, Any]]]:
    trees = config["frozen_tree_authorities"]
    result: dict[int, tuple[Path, dict[str, Any]]] = {}
    fold0_root = root / trees["scbasset_fold0_predictions"]["path"]
    fold0_bundle_root = fold0_root / "predictions/bundle"
    fold0 = load_json(fold0_bundle_root / "prediction_bundle.json", label="scBasset fold0 bundle")
    valid_fold = int(fold0.get("metadata", {}).get("donor_valid_fold", -1))
    result[valid_fold] = (fold0_bundle_root, fold0)

    lock_root = root / trees["scbasset_fold1to4_prediction_lock"]["path"]
    lock = load_json(lock_root / "prediction_lock.json", label="scBasset lock")
    for record in lock.get("records", ()):
        prediction_root = safe_path(
            root, str(record.get("prediction_root", "")), label="scBasset prediction root"
        )
        manifest = load_json(
            prediction_root / "ARTIFACTS.json", label="scBasset prediction manifest"
        )
        if (
            not (prediction_root / "COMPLETE").is_file()
            or sha256_file(prediction_root / "ARTIFACTS.json")
            != record.get("prediction_artifacts_sha256")
            or manifest.get("metadata", {}).get("artifact_class")
            != "scbasset_fold_predictions"
        ):
            raise ComplementarityReadinessError("scBasset lock binding changed")
        bundle_root = prediction_root / "predictions/bundle"
        bundle = load_json(bundle_root / "prediction_bundle.json", label="scBasset bundle")
        fold = int(record.get("valid_fold", -1))
        if fold in result:
            raise ComplementarityReadinessError("duplicate scBasset valid fold")
        result[fold] = (bundle_root, bundle)
    if set(result) != set(range(5)):
        raise ComplementarityReadinessError("scBasset valid-fold coverage differs")
    return result


def compare_scbasset_midas_rows(
    root: Path, config: Mapping[str, Any]
) -> list[dict[str, Any]]:
    trees = config["frozen_tree_authorities"]
    midas_root = root / trees["midas_inductive_predictions"]["path"] / "predictions"
    scbasset = _scbasset_bundles(root, config)
    rows: list[dict[str, Any]] = []
    for held_fold in range(5):
        midas_bundle_root = midas_root / f"fold_{held_fold}"
        midas_bundle = load_json(
            midas_bundle_root / "prediction_bundle.json", label=f"MIDAS fold {held_fold} bundle"
        )
        scb_root, scb_bundle = scbasset[held_fold]
        if (
            midas_bundle.get("metadata", {}).get("held_out_fold") != held_fold
            or scb_bundle.get("metadata", {}).get("donor_valid_fold") != held_fold
        ):
            raise ComplementarityReadinessError("held donor fold mapping differs")
        midas_path = _bundle_row_ids(midas_bundle_root, midas_bundle)
        scb_path = _bundle_row_ids(scb_root, scb_bundle)
        midas_rows: dict[str, str] = {}
        for row_hash, donor_hash in _row_id_records(midas_path):
            if row_hash in midas_rows:
                raise ComplementarityReadinessError("duplicate MIDAS row hash")
            midas_rows[row_hash] = donor_hash
        scb_count = 0
        scb_seen: set[str] = set()
        intersection = 0
        donor_mismatches = 0
        for row_hash, donor_hash in _row_id_records(scb_path):
            if row_hash in scb_seen:
                raise ComplementarityReadinessError("duplicate scBasset row hash")
            scb_seen.add(row_hash)
            scb_count += 1
            if row_hash in midas_rows:
                intersection += 1
                donor_mismatches += midas_rows[row_hash] != donor_hash
        midas_count = len(midas_rows)
        rows.append(
            {
                "held_donor_fold": held_fold,
                "midas_rows": midas_count,
                "scbasset_rows": scb_count,
                "intersection_rows": intersection,
                "midas_coverage_fraction": intersection / midas_count,
                "scbasset_coverage_fraction": intersection / scb_count,
                "identical_row_set": midas_count == scb_count == intersection,
                "midas_subset_of_scbasset": intersection == midas_count,
                "donor_hash_mismatches": donor_mismatches,
                "same_table_schema_sha256": midas_bundle.get("table_schema_sha256")
                == scb_bundle.get("table_schema_sha256"),
                "same_unit_id_namespace": midas_bundle.get("unit_id_namespace")
                == scb_bundle.get("unit_id_namespace"),
                "same_prediction_scale": midas_bundle.get("metadata", {}).get(
                    "prediction_scale"
                )
                == scb_bundle.get("metadata", {}).get("prediction_scale"),
            }
        )
    return rows


def tsv(rows: Sequence[Mapping[str, Any]], fields: Sequence[str]) -> str:
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        buffer, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        writer.writerow({field: row.get(field, "") for field in fields})
    return buffer.getvalue()


def preflight(
    *, root: Path, config_path: Path, output: Path | None = None
) -> dict[str, Any]:
    root = root.resolve(strict=True)
    config = load_toml(config_path, label="complementarity readiness config")
    validate_config(config)
    authorities = validate_authorities(root, config)
    validate_trigger_authorities(root, config)
    metadata = validate_prediction_metadata(root, config)
    overlap_rows = compare_scbasset_midas_rows(root, config)
    lane_rows = [dict(row) for row in config["prediction_lanes"]]
    pair_rows = [dict(row) for row in config["pair_audits"]]
    receipt = {
        "schema_version": "masld-bench-complementarity-readiness-audit-receipt-v1",
        "status": STATUS,
        "audit_id": config["audit_id"],
        "authority_rows": len(authorities),
        "prediction_lane_rows": len(lane_rows),
        "pair_audit_rows": len(pair_rows),
        "required_stacking_source_rows": len(SOURCE_FIELDS),
        "bound_stacking_source_rows": 0,
        "trigger_ready_pair_rows": 0,
        "prediction_metadata": metadata,
        "scbasset_midas_row_overlap": overlap_rows,
        "scbasset_midas_any_identical_row_set": any(
            row["identical_row_set"] for row in overlap_rows
        ),
        "scbasset_midas_all_rows_share_prediction_scale": all(
            row["same_prediction_scale"] for row in overlap_rows
        ),
        "residual_correlation_calculated": False,
        "stack_weights_fit": False,
        "prediction_values_read": False,
        "row_identifiers_read": True,
        "raw_or_summary_outcomes_read": False,
        "sealed_assets_read": False,
        "novel_model_built_or_fit": False,
        "conditional_model_build_authorized": False,
        "next_required_artifact": "complete_task_specific_open_development_shortlists_before_residual_or_stack_production",
        "bundled_cpu_validation_ask": config["bundled_cpu_validation_ask"],
    }
    if output is not None:
        output.mkdir(parents=True, exist_ok=False)
        write_json_exclusive(output / "receipt.json", receipt)
        write_text_exclusive(
            output / "prediction_lane_inventory.tsv",
            tsv(
                lane_rows,
                (
                    "lane_id",
                    "task_id",
                    "dataset_ids",
                    "model_ids",
                    "source_family_role",
                    "prediction_schema",
                    "outer_folds_observed",
                    "seed_ids_observed",
                    "row_identity_ready",
                    "residual_ready",
                    "stack_gain_ready",
                    "blocking_findings",
                ),
            ),
        )
        write_text_exclusive(
            output / "pair_audit.tsv",
            tsv(
                pair_rows,
                (
                    "pair_id",
                    "left_lane",
                    "right_lane",
                    "same_task",
                    "same_prediction_schema",
                    "distinct_sequence_context_families",
                    "trigger_ready",
                    "disposition",
                ),
            ),
        )
        write_text_exclusive(
            output / "scbasset_midas_row_overlap.tsv",
            tsv(
                overlap_rows,
                (
                    "held_donor_fold",
                    "midas_rows",
                    "scbasset_rows",
                    "intersection_rows",
                    "midas_coverage_fraction",
                    "scbasset_coverage_fraction",
                    "identical_row_set",
                    "midas_subset_of_scbasset",
                    "donor_hash_mismatches",
                    "same_table_schema_sha256",
                    "same_unit_id_namespace",
                    "same_prediction_scale",
                ),
            ),
        )
    return receipt


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    receipt = preflight(root=args.root, config_path=args.config, output=args.output)
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
