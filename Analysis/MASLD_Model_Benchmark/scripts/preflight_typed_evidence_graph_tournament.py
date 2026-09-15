#!/usr/bin/env python3
"""Fail closed unless a leakage-safe typed-evidence graph task is executable."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path
import re
import tomllib
from typing import Any, Mapping, Sequence

from masld_bench.artifacts import (
    ArtifactError,
    reject_symlink_components,
    write_json_exclusive,
    write_text_exclusive,
)


SCHEMA = "masld-bench-typed-evidence-graph-tournament-preflight-v1"
BLOCKED_STATUS = "blocked_no_executable_wholly_removed_source_task"
MODEL_ORDER = (
    "evidence_count",
    "nearest_gene",
    "graph_degree",
    "graph_source_count",
    "regularized_linear_evidence",
    "node2vec",
    "bionic",
    "graphsage",
    "rgcn",
    "hgt",
    "typed_transformer",
    "late_fusion",
)
MANDATORY_BASELINES = MODEL_ORDER[:5]
EVIDENCE_GRAPH_MODELS = MODEL_ORDER[5:]
ARTIFACT_GROUPS = (
    "evidence_count",
    "nearest_gene",
    "native_baselines",
    "node2vec",
    "bionic",
    "graphsage",
    "rgcn",
    "hgt",
    "typed_transformer",
    "late_fusion",
)
BASE_AUTHORITIES = {
    "evidence_graph_registry",
    "mandatory_baseline_registry",
    "task",
    "split",
    "promotion_gates",
    "current_coloc_gwas_dataset",
    "gse281160_dataset",
    "gse281364_dataset",
}
EXPECTED_AUTHORITIES = BASE_AUTHORITIES | {
    f"{group}_{suffix}"
    for group in ARTIFACT_GROUPS
    for suffix in ("checkpoints", "crosswalk", "exposure")
}
EXPECTED_GATES = (
    "promoted_corrected_coloc_release",
    "immutable_graph_snapshot",
    "frozen_evidence_ontology_and_direction",
    "source_family_alias_and_derivative_map",
    "source_independent_candidate_universe",
    "wholly_removed_source_identity",
    "explicit_supported_and_tested_negative_denominator",
    "held_source_removal_receipt",
    "evidence_provenance_applicability_observability_and_missing_states",
    "ancestry_aware_native_build_ld_blocks",
    "gene_blocks_and_held_locus_incident_edge_removal",
    "degree_source_count_and_observability_matched_controls",
    "split_safe_negative_sampling",
    "source_blind_calibration_contract",
    "implementation_runtime_and_tensor_only_artifact_receipts",
    "model_audit_source_record_reconciliation",
    "wholly_sealed_external_graph_source",
)
EXPECTED_BINDINGS = {
    "graph_snapshot",
    "evidence_ontology",
    "source_family_map",
    "candidate_universe",
    "held_source_manifest",
    "explicit_outcome_denominator",
    "observability_manifest",
    "ld_block_manifest",
    "gene_block_manifest",
    "matched_control_manifest",
    "source_removal_fixture",
    "negative_sampling_fixture",
}
EXPECTED_CANDIDATES = {"current_coloc_gwas", "gse281160", "gse281364"}
SHA256 = re.compile(r"^[0-9a-f]{64}$")


class GraphTournamentPreflightError(RuntimeError):
    """Raised when an outcome-blind graph authority differs from the selection record."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise GraphTournamentPreflightError(f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise GraphTournamentPreflightError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise GraphTournamentPreflightError(f"{label} must be a JSON object")
    return value


def load_toml(path: Path, *, label: str) -> dict[str, Any]:
    if path.is_symlink() or not path.is_file():
        raise GraphTournamentPreflightError(f"{label} is not a regular file")
    try:
        value = tomllib.loads(path.read_text(encoding="utf-8"))
    except (OSError, tomllib.TOMLDecodeError) as error:
        raise GraphTournamentPreflightError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise GraphTournamentPreflightError(f"{label} must be a TOML table")
    return value


def safe_project_path(root: Path, relative: str) -> Path:
    value = Path(relative)
    if value.is_absolute() or not value.parts or ".." in value.parts:
        raise GraphTournamentPreflightError(f"unsafe project path: {relative}")
    try:
        path = reject_symlink_components(root / value, label="graph authority")
        path.resolve(strict=True).relative_to(root)
    except (ArtifactError, OSError, ValueError) as error:
        raise GraphTournamentPreflightError(
            f"graph authority is missing or escapes root: {relative}"
        ) from error
    return path


def validate_config(config: Mapping[str, Any]) -> None:
    if (
        config.get("schema_version") != SCHEMA
        or config.get("status") != BLOCKED_STATUS
        or config.get("task_id") != "typed_evidence_graph"
        or config.get("strongest_currently_executable_task") is not None
        or tuple(config.get("tournament_model_order", ())) != MODEL_ORDER
        or set(config.get("models", {})) != set(MODEL_ORDER)
    ):
        raise GraphTournamentPreflightError("tournament identity or roster differs")

    candidate = config.get("first_activation_candidate", {})
    if (
        candidate.get("candidate_id")
        != "current_coloc_gwas_leave_one_study_family_out"
        or candidate.get("status") != "blocked"
        or candidate.get("held_source_unit")
        != "one_frozen_independent_study_source_family"
        or "explicit tested-negative" not in candidate.get("reason", "")
        or "Absence from a study is never a negative"
        not in candidate.get("activation_boundary", "")
    ):
        raise GraphTournamentPreflightError("activation candidate differs")

    models = config["models"]
    if (
        tuple(model for model in MODEL_ORDER if models[model].get("mandatory"))
        != MANDATORY_BASELINES
        or models["node2vec"].get("evaluation_support")
        != "common_support_only_source_only_and_isolated_nodes_untestable"
        or models["bionic"].get("evaluation_support")
        != "common_support_only_unsupported_nodes_untestable"
        or "source_independent_root_features"
        not in models["graphsage"].get("evaluation_support", "")
        or models["rgcn"].get("evaluation_support")
        != "cannot_directly_complete_an_unseen_relation_type"
        or models["hgt"].get("evaluation_support")
        != "cannot_score_untrained_unseen_node_or_edge_types"
        or not models["typed_transformer"].get("current_status", "").startswith(
            "blocked_no_implementation_identity"
        )
        or "identical_source_removal"
        not in models["late_fusion"].get("evaluation_support", "")
    ):
        raise GraphTournamentPreflightError("model topology boundary differs")

    source_audit = config.get("candidate_source_audit", {})
    if set(source_audit) != EXPECTED_CANDIDATES or any(
        row.get("currently_executable") is not False
        for row in source_audit.values()
    ):
        raise GraphTournamentPreflightError("candidate source disposition differs")
    if (
        source_audit["current_coloc_gwas"].get("expected_biological_units") != 50
        or source_audit["gse281160"].get("expected_biological_units")
        != "UNRESOLVED"
        or source_audit["gse281364"].get("expected_biological_units") != 16
        or "MPRA_activity_is_not_an_endogenous_enhancer_gene_link"
        not in source_audit["gse281364"].get("blocking_findings", ())
    ):
        raise GraphTournamentPreflightError("candidate source evidence differs")

    lineage = config.get("authority_lineage_audit", {})
    if (
        lineage.get("mandatory_baseline_registry_current_sha256")
        != "cd7e67aa068c2a041fbf842564fc20d82be008020bb5226d70bc36c937e4a74d"
        or lineage.get("graph_baseline_audits_recorded_registry_sha256")
        != "8e89cc25c298627617f29927604dc7aa5288a89eb3c1bd3b87b3d5bc0ef6915e"
        or lineage.get("current_registry_matches_recorded_registry") is not False
        or lineage.get("graph_model_rows_revalidated_from_current_registry")
        is not True
        or not lineage.get("status", "").startswith("stale_whole_registry_hash")
    ):
        raise GraphTournamentPreflightError("model-audit authority lineage differs")

    gates = config.get("activation_gates", ())
    if (
        tuple(row.get("gate_id") for row in gates) != EXPECTED_GATES
        or any(row.get("required") is not True for row in gates[:-1])
        or any(row.get("status") != "unresolved" for row in gates[:-1])
        or gates[-1]
        != {
            "gate_id": "wholly_sealed_external_graph_source",
            "status": "absent",
            "required": False,
        }
    ):
        raise GraphTournamentPreflightError("activation gate census differs")

    bindings = config.get("required_artifact_bindings", {})
    if set(bindings) != EXPECTED_BINDINGS or any(
        value is not None for value in bindings.values()
    ):
        raise GraphTournamentPreflightError(
            "unresolved artifact bindings must remain explicit nulls"
        )

    evaluation = config.get("evaluation_contract", {})
    if (
        evaluation.get("primary_metric") != "held_source_auprc"
        or tuple(evaluation.get("mandatory_baselines", ())) != MANDATORY_BASELINES
        or "1Mb_or_r2_at_least_0.8" not in evaluation.get("grouping", "")
        or not evaluation.get("uncertainty", "").startswith("10000_paired_two_way")
        or evaluation.get("promotion")
        != "none_without_a_wholly_sealed_evidence_source"
        or evaluation.get("causal_claim") is not False
        or evaluation.get("clinical_claim") is not False
    ):
        raise GraphTournamentPreflightError("evaluation contract differs")

    firewall = config.get("outcome_firewall", {})
    if not firewall or any(value is not False for value in firewall.values()):
        raise GraphTournamentPreflightError("outcome firewall differs")

    ask = config.get("bundled_cpu_validation_ask", {})
    if (
        ask.get("job_name") != "model-check-217"
        or ask.get("partition") != "cpu"
        or ask.get("qos") != "nslab"
        or ask.get("cpus_per_task") != 1
        or ask.get("memory_gb") != 8
        or ask.get("walltime") != "02:00:00"
        or ask.get("array") is not False
        or ask.get("job_count") != 1
        or ask.get("total_cpu_hours") != 2
        or ask.get("total_gpu_hours") != 0
        or ask.get("submitted") is not False
    ):
        raise GraphTournamentPreflightError("bundled CPU validation ask differs")


def validate_authorities(
    root: Path, config: Mapping[str, Any]
) -> dict[str, dict[str, Any]]:
    authorities = config.get("frozen_file_authorities", {})
    if set(authorities) != EXPECTED_AUTHORITIES:
        raise GraphTournamentPreflightError("frozen authority census differs")
    result: dict[str, dict[str, Any]] = {}
    for name, binding in authorities.items():
        if not isinstance(binding, Mapping):
            raise GraphTournamentPreflightError(f"authority binding differs: {name}")
        path = safe_project_path(root, str(binding.get("path", "")))
        expected = str(binding.get("sha256", ""))
        if not SHA256.fullmatch(expected) or sha256_file(path) != expected:
            raise GraphTournamentPreflightError(f"frozen authority changed: {name}")
        result[name] = {
            "path": path.relative_to(root).as_posix(),
            "sha256": expected,
            "size_bytes": path.stat(follow_symlinks=False).st_size,
        }
    return result


def model_rows_from_registry(path: Path) -> dict[str, Mapping[str, Any]]:
    registry = load_toml(path, label=f"model registry {path.name}")
    result: dict[str, Mapping[str, Any]] = {}
    for row in registry.get("models", ()):  # type: ignore[union-attr]
        model_id = row.get("model_id")
        if model_id in MODEL_ORDER:
            if model_id in result:
                raise GraphTournamentPreflightError(
                    f"duplicate graph model registry row: {model_id}"
                )
            result[model_id] = row
    return result


def validate_model_registries(root: Path, config: Mapping[str, Any]) -> list[dict[str, Any]]:
    paths = config["frozen_file_authorities"]
    baseline_rows = model_rows_from_registry(
        root / paths["mandatory_baseline_registry"]["path"]
    )
    graph_rows = model_rows_from_registry(
        root / paths["evidence_graph_registry"]["path"]
    )
    if set(baseline_rows) != set(MANDATORY_BASELINES):
        raise GraphTournamentPreflightError("mandatory graph baseline roster differs")
    if set(graph_rows) != set(EVIDENCE_GRAPH_MODELS):
        raise GraphTournamentPreflightError("graph architecture roster differs")

    output: list[dict[str, Any]] = []
    for model_id in MODEL_ORDER:
        expected = config["models"][model_id]
        row = baseline_rows.get(model_id, graph_rows.get(model_id))
        if row is None:
            raise GraphTournamentPreflightError(f"model is absent: {model_id}")
        if (
            row.get("status") != "candidate"
            or row.get("admission_blocking") is not True
            or "typed_evidence_graph" not in row.get("supported_tasks", ())
            or not row.get("blockers")
            or expected.get("current_status", "").split("_", 1)[0] != "blocked"
        ):
            raise GraphTournamentPreflightError(
                f"model admission state differs: {model_id}"
            )
        if model_id == "typed_transformer" and (
            row.get("checkpoint_revision") != "UNRESOLVED"
            or row.get("license_status") != "UNRESOLVED"
        ):
            raise GraphTournamentPreflightError(
                "typed transformer must remain unresolved until separately admitted"
            )
        output.append(
            {
                "model_id": model_id,
                "registry": expected["registry"],
                "model_class": expected["model_class"],
                "evaluation_support": expected["evaluation_support"],
                "mandatory": expected["mandatory"],
                "registry_status": row["status"],
                "admission_blocking": row["admission_blocking"],
                "implementation_type": row["implementation_type"],
                "checkpoint_revision": row["checkpoint_revision"],
                "license_status": row["license_status"],
                "exposure_status": row["exposure_status"],
                "current_status": expected["current_status"],
            }
        )
    return output


def validate_task_split_and_promotion(root: Path, config: Mapping[str, Any]) -> None:
    paths = config["frozen_file_authorities"]
    task = load_toml(root / paths["task"]["path"], label="graph TaskSpec")
    split = load_toml(root / paths["split"]["path"], label="graph split")
    promotion = load_toml(
        root / paths["promotion_gates"]["path"], label="promotion gates"
    ).get("tasks", {}).get("typed_evidence_graph", {})
    if (
        task.get("task_id") != "typed_evidence_graph"
        or task.get("status") != "blocked"
        or task.get("unit_of_inference") != "independent_locus"
        or task.get("datasets_train") != ["current_coloc_gwas"]
        or task.get("datasets_development") != ["gse281364", "gse281160"]
        or task.get("datasets_sealed") != []
        or tuple(task.get("baseline_model_ids", ()))
        != (
            "evidence_count",
            "nearest_gene",
            "graph_degree",
            "graph_source_count",
            "regularized_linear_evidence",
        )
        or task.get("bootstrap_replicates") != 10000
    ):
        raise GraphTournamentPreflightError("typed-evidence TaskSpec differs")
    if (
        split.get("split_id") != "evidence_graph_outer"
        or split.get("entity") != "independent_locus"
        or split.get("group_keys")
        != ["graph_snapshot_id", "native_build", "ld_block_id", "locus_id"]
        or "Held-locus labels and every edge derived from them are removed"
        not in split.get("label_policy", "")
        or "not sealed graph labels" not in split.get("exposure_policy", "")
    ):
        raise GraphTournamentPreflightError("typed-evidence split differs")
    if (
        promotion.get("allowed_promotion_modes")
        != ["blocked_missing_sealed_source", "conditional_sealed"]
        or promotion.get("sealed_holdout_dataset_ids") != []
        or promotion.get("delta_metric") != "held_source_auprc_gain"
        or promotion.get("delta_min") != 0.03
        or promotion.get("brier_delta_max") != 0.01
    ):
        raise GraphTournamentPreflightError("typed-evidence promotion gate differs")


def validate_dataset_dispositions(root: Path, config: Mapping[str, Any]) -> None:
    paths = config["frozen_file_authorities"]
    coloc = load_toml(
        root / paths["current_coloc_gwas_dataset"]["path"],
        label="current COLOC DatasetManifest",
    )
    crop = load_toml(
        root / paths["gse281160_dataset"]["path"], label="GSE281160 DatasetManifest"
    )
    mpra = load_toml(
        root / paths["gse281364_dataset"]["path"], label="GSE281364 DatasetManifest"
    )
    if (
        coloc.get("expected_biological_units") != 50
        or coloc.get("biological_unit") != "study"
        or coloc.get("admission_blocking") is not True
        or coloc.get("exposure_status") != "downstream_demo"
        or "Corrected 50-study COLOC rerun must be promoted"
        not in " ".join(coloc.get("blockers", ()))
    ):
        raise GraphTournamentPreflightError("current COLOC disposition differs")
    if (
        crop.get("expected_biological_units") != "UNRESOLVED"
        or crop.get("biological_unit") != "experimental_replicate"
        or crop.get("pairing_levels") != ["same_cell"]
        or crop.get("admission_blocking") is not True
        or crop.get("exposure_status") != "downstream_demo"
    ):
        raise GraphTournamentPreflightError("GSE281160 disposition differs")
    if (
        mpra.get("expected_biological_units") != 16
        or mpra.get("biological_unit") != "experimental_replicate"
        or mpra.get("pairing_levels") != ["same_sample_different_aliquot"]
        or mpra.get("admission_blocking") is not True
        or mpra.get("exposure_status") != "downstream_demo"
    ):
        raise GraphTournamentPreflightError("GSE281364 disposition differs")


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
    config = load_json(config_path, label="graph tournament preflight config")
    validate_config(config)
    authorities = validate_authorities(root, config)
    validate_task_split_and_promotion(root, config)
    validate_dataset_dispositions(root, config)
    model_rows = validate_model_registries(root, config)
    gate_rows = [dict(row) for row in config["activation_gates"]]
    source_rows = [
        {
            "dataset_id": dataset_id,
            "currently_executable": row["currently_executable"],
            "role": row["role"],
            "unit": row["unit"],
            "expected_biological_units": row["expected_biological_units"],
            "blocking_findings": "|".join(row["blocking_findings"]),
        }
        for dataset_id, row in config["candidate_source_audit"].items()
    ]
    receipt = {
        "schema_version": "masld-bench-typed-evidence-graph-tournament-preflight-receipt-v1",
        "status": BLOCKED_STATUS,
        "task_id": "typed_evidence_graph",
        "strongest_currently_executable_task": None,
        "first_activation_candidate": config["first_activation_candidate"],
        "authority_lineage_audit": config["authority_lineage_audit"],
        "model_census_rows": len(model_rows),
        "mandatory_baseline_rows": sum(bool(row["mandatory"]) for row in model_rows),
        "candidate_source_rows": len(source_rows),
        "activation_gate_rows": len(gate_rows),
        "required_gates_unresolved": sum(
            row["required"] and row["status"] != "pass" for row in gate_rows
        ),
        "sealed_graph_source_registered": False,
        "graph_snapshot_bound": False,
        "wholly_removed_source_bound": False,
        "explicit_tested_negative_denominator_bound": False,
        "training_or_scoring_allowed": False,
        "champion_claim_allowed": False,
        "causal_claim_allowed": False,
        "clinical_claim_allowed": False,
        "raw_or_summary_labels_read": False,
        "sealed_assets_read": False,
        "model_import_fit_predict_or_score": False,
        "network_download_or_install": False,
        "frozen_authorities": authorities,
        "bundled_cpu_validation_ask": config["bundled_cpu_validation_ask"],
    }
    if output is not None:
        try:
            output = reject_symlink_components(output, label="graph preflight output")
        except ArtifactError as error:
            raise GraphTournamentPreflightError(str(error)) from error
        if output.exists():
            raise GraphTournamentPreflightError(
                f"refusing to overwrite graph preflight output: {output}"
            )
        if not output.parent.is_dir():
            raise GraphTournamentPreflightError(
                f"graph preflight output parent is missing: {output.parent}"
            )
        output.mkdir(mode=0o750)
        write_text_exclusive(
            output / "model_census.tsv",
            tsv(
                model_rows,
                (
                    "model_id",
                    "registry",
                    "model_class",
                    "evaluation_support",
                    "mandatory",
                    "registry_status",
                    "admission_blocking",
                    "implementation_type",
                    "checkpoint_revision",
                    "license_status",
                    "exposure_status",
                    "current_status",
                ),
            ),
        )
        write_text_exclusive(
            output / "candidate_source_audit.tsv",
            tsv(
                source_rows,
                (
                    "dataset_id",
                    "currently_executable",
                    "role",
                    "unit",
                    "expected_biological_units",
                    "blocking_findings",
                ),
            ),
        )
        write_text_exclusive(
            output / "activation_gates.tsv",
            tsv(gate_rows, ("gate_id", "status", "required")),
        )
        write_json_exclusive(output / "receipt.json", receipt)
    return receipt


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", type=Path)
    return parser.parse_args()


def main() -> None:
    arguments = parse_args()
    receipt = preflight(
        root=arguments.root,
        config_path=arguments.config,
        output=arguments.output,
    )
    print(json.dumps(receipt, sort_keys=True))


if __name__ == "__main__":
    main()
