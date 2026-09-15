#!/usr/bin/env python3
"""Independently freeze outcome-blind readiness for external cohort inputs.

This auditor deliberately does not import either cohort's preparation code. It
reads frozen manifests, no-fit receipts, topology axes, and missingness masks.
Frozen-tree verification hashes every registered output file byte, but the auditor
never deserializes expression matrices or reads phenotype/outcome tables, fits,
predictions, or evaluator results. A passing receipt grants only the narrow
development scope stated in that receipt.
"""

from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import re
from typing import Any, Iterable, Mapping

import numpy as np

from masld_bench.artifacts import (
    freeze_tree,
    verify_frozen_tree,
    write_json_exclusive,
)
from masld_bench.hashing import sha256_file


GSE105127_PLAN_FIELDS = (
    "bundle_id",
    "row_id",
    "participant_group_id",
    "zone",
    "pairing_topology",
)
GSE105127_FORBIDDEN_FIELDS = frozenset(
    {
        "phenotype",
        "label",
        "outcome",
        "disease",
        "fibrosis",
        "nas",
        "sex",
        "age",
        "bmi",
    }
)
GSE105127_MODEL_INPUT_FIELDS = (
    "row_id",
    "participant_group_id",
    "zone",
    "pairing_topology",
    "rna_row_index",
    "rrbs_collapsed_path",
    "rrbs_collapsed_artifacts_sha256",
    "cpg_crosswalk_path",
    "rna_matrix_root",
)
GSE105127_FOLD_FIELDS = (
    "row_id",
    "participant_group_id",
    "zone",
    "participant_outer_fold",
    "fold_assignment_uses_outcomes",
    "available_to_model_as_feature",
)
GSE105127_ROW_AXIS_FIELDS = (
    "row_index",
    "row_id",
    "participant_group_id",
    "zone",
    "pairing_topology",
)
GSE268273_PLAN_FIELDS = (
    "bundle_id",
    "row_id",
    "experiment_accession",
    "biosample_accession",
    "technical_runs",
    "fastq_files",
    "fastq_bytes",
    "effective_library_layout",
    "aggregation_rule",
)
GSE268273_PARTICIPANT_FIELDS = (
    "row_index",
    "row_id",
    "cohort_family_id",
    "rna_observation_state",
    "raw_technical_run_count",
    "raw_fastq_bytes",
    "quantification_measurement",
)


class ExternalCohortReadinessError(ValueError):
    """Raised when a cohort cannot receive the narrow readiness disposition."""


def read_tsv(path: Path) -> tuple[tuple[str, ...], list[dict[str, str]]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ExternalCohortReadinessError(f"TSV has no header: {path}")
        return tuple(reader.fieldnames), list(reader)


def load_json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ExternalCohortReadinessError(f"JSON root is not an object: {path}")
    return value


def require_frozen(
    root: Path,
    *,
    artifact_class: str,
    expected_sha256: str | None = None,
) -> dict[str, Any]:
    if expected_sha256 is not None and sha256_file(root / "ARTIFACTS.json") != expected_sha256:
        raise ExternalCohortReadinessError(
            f"pinned ARTIFACTS hash differs: {root}"
        )
    manifest = verify_frozen_tree(root)
    if manifest.get("metadata", {}).get("artifact_class") != artifact_class:
        raise ExternalCohortReadinessError(
            f"artifact class differs for {root}: {artifact_class}"
        )
    return manifest


def require_false(receipt: Mapping[str, Any], fields: Iterable[str], label: str) -> None:
    for field in fields:
        if receipt.get(field) is not False:
            raise ExternalCohortReadinessError(f"{label} firewall differs: {field}")


def validate_gse105127_plan_rows(
    rna_fields: tuple[str, ...],
    rna_rows: list[dict[str, str]],
    rrbs_fields: tuple[str, ...],
    rrbs_rows: list[dict[str, str]],
) -> dict[str, Any]:
    if GSE105127_FORBIDDEN_FIELDS & (set(rna_fields) | set(rrbs_fields)):
        raise ExternalCohortReadinessError("GSE105127 plan exposes forbidden fields")
    if not set(GSE105127_PLAN_FIELDS) <= set(rna_fields) or not set(
        GSE105127_PLAN_FIELDS
    ) <= set(rrbs_fields):
        raise ExternalCohortReadinessError("GSE105127 plan topology fields differ")
    if len(rna_rows) != 57 or len(rrbs_rows) != 57:
        raise ExternalCohortReadinessError("GSE105127 row census differs")
    rna_by_id = {row["row_id"]: row for row in rna_rows}
    rrbs_by_id = {row["row_id"]: row for row in rrbs_rows}
    if len(rna_by_id) != 57 or set(rna_by_id) != set(rrbs_by_id):
        raise ExternalCohortReadinessError("GSE105127 RNA/RRBS row join differs")
    zones: dict[str, set[str]] = defaultdict(set)
    bundles: dict[str, set[str]] = defaultdict(set)
    for row_id in sorted(rna_by_id):
        rna = rna_by_id[row_id]
        rrbs = rrbs_by_id[row_id]
        identity = ("participant_group_id", "zone", "pairing_topology")
        if any(rna[field] != rrbs[field] for field in identity):
            raise ExternalCohortReadinessError("GSE105127 modality identity differs")
        if rna["pairing_topology"] != "adjacent_section":
            raise ExternalCohortReadinessError("GSE105127 false pairing detected")
        participant = rna["participant_group_id"]
        zones[participant].add(rna["zone"])
        bundles[participant].update((rna["bundle_id"], rrbs["bundle_id"]))
    if len(zones) != 19 or any(value != {"CV", "IZ", "PP"} for value in zones.values()):
        raise ExternalCohortReadinessError("GSE105127 participant-zone topology differs")
    if any(len(value) != 1 for value in bundles.values()):
        raise ExternalCohortReadinessError("GSE105127 participant split across bundles")
    if len({row.get("run_accession", "") for row in rna_rows}) != 57:
        raise ExternalCohortReadinessError("GSE105127 RNA run identity differs")
    if any(
        row.get("library_layout") != "SINGLE"
        or row.get("read_length") != "76"
        or int(row.get("read_count", "0")) <= 0
        or int(row.get("fastq_bytes", "0")) <= 0
        or re.fullmatch(r"[0-9a-f]{32}", row.get("fastq_md5", "")) is None
        for row in rna_rows
    ):
        raise ExternalCohortReadinessError("GSE105127 RNA run mechanics differ")
    return {
        "participants": 19,
        "participant_zone_rows": 57,
        "rna_runs": 57,
        "rrbs_libraries": 57,
        "zones_per_participant": 3,
        "pairing_topology": "adjacent_section",
    }


def validate_gse105127_fold_rows(
    fields: tuple[str, ...], rows: list[dict[str, str]]
) -> dict[str, int]:
    if fields != GSE105127_FOLD_FIELDS or len(rows) != 57:
        raise ExternalCohortReadinessError("GSE105127 fold authority schema differs")
    folds: dict[str, set[int]] = defaultdict(set)
    for row in rows:
        if (
            row["fold_assignment_uses_outcomes"].lower() != "false"
            or row["available_to_model_as_feature"].lower() != "false"
        ):
            raise ExternalCohortReadinessError("GSE105127 fold firewall differs")
        fold = int(row["participant_outer_fold"])
        if fold not in range(5):
            raise ExternalCohortReadinessError("GSE105127 outer fold differs")
        folds[row["participant_group_id"]].add(fold)
    if len(folds) != 19 or any(len(value) != 1 for value in folds.values()):
        raise ExternalCohortReadinessError("GSE105127 participant crosses folds")
    counts = Counter(next(iter(value)) for value in folds.values())
    if sorted(counts.values()) != [3, 4, 4, 4, 4]:
        raise ExternalCohortReadinessError("GSE105127 fold balance differs")
    return {str(fold): counts[fold] for fold in range(5)}


def validate_gse268273_topology(
    fields: tuple[str, ...], rows: list[dict[str, str]]
) -> dict[str, Any]:
    if fields != GSE268273_PLAN_FIELDS or len(rows) != 109:
        raise ExternalCohortReadinessError("GSE268273 plan schema or census differs")
    if (
        len({row["row_id"] for row in rows}) != 109
        or len({row["experiment_accession"] for row in rows}) != 109
        or len({row["biosample_accession"] for row in rows}) != 109
    ):
        raise ExternalCohortReadinessError("GSE268273 participant identity differs")
    layouts = Counter(row["effective_library_layout"] for row in rows)
    run_distribution = Counter(int(row["technical_runs"]) for row in rows)
    if layouts != Counter({"single_end": 73, "paired_end": 36}):
        raise ExternalCohortReadinessError("GSE268273 effective layout differs")
    if run_distribution != Counter({8: 49, 4: 36, 12: 24}):
        raise ExternalCohortReadinessError("GSE268273 run distribution differs")
    for row in rows:
        runs = int(row["technical_runs"])
        files = int(row["fastq_files"])
        if row["effective_library_layout"] == "single_end":
            expected_files = runs
            expected_rule = "comma_ordered_single_end_technical_runs_then_one_participant_RSEM_library"
        else:
            expected_files = 2 * runs
            expected_rule = "matewise_comma_ordered_technical_runs_then_one_participant_RSEM_library"
        if files != expected_files or row["aggregation_rule"] != expected_rule:
            raise ExternalCohortReadinessError("GSE268273 run aggregation differs")
        if int(row["bundle_id"]) not in range(8) or int(row["fastq_bytes"]) <= 0:
            raise ExternalCohortReadinessError("GSE268273 bundle or byte count differs")
    totals = {
        "participants": len(rows),
        "technical_runs": sum(int(row["technical_runs"]) for row in rows),
        "fastq_files": sum(int(row["fastq_files"]) for row in rows),
        "fastq_bytes": sum(int(row["fastq_bytes"]) for row in rows),
        "effective_single_end_participants": layouts["single_end"],
        "effective_paired_end_participants": layouts["paired_end"],
    }
    if totals != {
        "participants": 109,
        "technical_runs": 824,
        "fastq_files": 968,
        "fastq_bytes": 512881099727,
        "effective_single_end_participants": 73,
        "effective_paired_end_participants": 36,
    }:
        raise ExternalCohortReadinessError("GSE268273 raw totals differ")
    return totals


def audit_gse105127(
    *,
    plan_root: Path,
    campaign_root: Path,
    source_audit_root: Path,
    output: Path,
    expected_plan_artifacts_sha256: str,
    expected_source_audit_artifacts_sha256: str,
) -> dict[str, Any]:
    if output.exists():
        raise ExternalCohortReadinessError("refusing to overwrite GSE105127 disposition")
    require_frozen(
        plan_root,
        artifact_class="gse105127_production_plan",
        expected_sha256=expected_plan_artifacts_sha256,
    )
    require_frozen(
        source_audit_root,
        artifact_class="gse105127_independent_source_dag_audit_execution",
        expected_sha256=expected_source_audit_artifacts_sha256,
    )
    source_audit_manifest = require_frozen(
        source_audit_root / "audit",
        artifact_class="gse105127_independent_source_dag_audit",
    )
    control_root = campaign_root / "control"
    activation_root = campaign_root / "activation"
    terminal_root = campaign_root / "terminal"
    crosswalk_root = campaign_root / "cpg_crosswalk"
    rna_root = campaign_root / "rna_consolidated"
    source_reference_root = plan_root.parent / "reference"
    rsem_reference_root = campaign_root / "rsem_reference"
    require_frozen(
        control_root, artifact_class="gse105127_gzip_preflight_revision_v2_control"
    )
    require_frozen(
        source_reference_root, artifact_class="gse105127_reference_bundle"
    )
    require_frozen(
        rsem_reference_root, artifact_class="gse105127_rsem_star_reference"
    )
    require_frozen(activation_root, artifact_class="gse105127_production_activation")
    require_frozen(
        terminal_root,
        artifact_class="gse105127_gzip_preflight_revision_v2_terminal",
    )
    require_frozen(
        crosswalk_root, artifact_class="gse105127_failure_aware_cpg_crosswalk"
    )
    require_frozen(rna_root, artifact_class="gse105127_rna_continuous_matrices")

    plan = load_json(plan_root / "plan.json")
    source_audit = load_json(
        source_audit_root / "audit/independent_source_dag_audit.json"
    )
    control = load_json(control_root / "control.json")
    source_reference = load_json(source_reference_root / "receipt.json")
    rsem_reference = load_json(rsem_reference_root / "receipt.json")
    activation = load_json(activation_root / "receipt.json")
    terminal = load_json(terminal_root / "receipt.json")
    crosswalk = load_json(crosswalk_root / "receipt.json")
    rna = load_json(rna_root / "receipt.json")
    require_false(plan, ("labels_accessed", "fit_or_score_performed"), "GSE105127 plan")
    require_false(
        source_audit,
        (
            "labels_read",
            "outcomes_read",
            "predictions_read",
            "fit_or_score_performed",
            "modeling_or_scoring_authorized_by_this_receipt",
            "sealed_evaluation_authorized_by_this_receipt",
            "clean_or_sealed_champion_eligible_by_this_receipt",
        ),
        "GSE105127 source audit",
    )
    require_false(
        control,
        (
            "labels_accessed",
            "outcomes_accessed",
            "fit_or_score_performed",
            "molecular_values_opened_during_control_preparation",
        ),
        "GSE105127 revision control",
    )
    require_false(
        source_reference,
        ("labels_accessed", "fit_or_score_performed"),
        "GSE105127 source reference",
    )
    require_false(
        rsem_reference,
        ("labels_accessed", "fit_or_score_performed"),
        "GSE105127 RSEM reference",
    )
    require_false(
        activation,
        ("labels_accessed", "fit_or_score_performed", "model_fit_run", "model_score_run"),
        "GSE105127 activation",
    )
    require_false(
        terminal,
        ("labels_accessed", "outcomes_accessed", "fit_or_score_performed"),
        "GSE105127 terminal",
    )
    require_false(crosswalk, ("labels_accessed", "fit_or_score_performed"), "GSE105127 crosswalk")
    require_false(rna, ("labels_accessed", "fit_or_score_performed"), "GSE105127 RNA")
    if (
        plan.get("status") != "planned_no_fit"
        or source_audit_manifest.get("metadata", {}).get("status")
        != "passed_audit_source_mechanics_verified_promotion_blocked"
        or source_audit.get("status")
        != "passed_audit_source_mechanics_verified_promotion_blocked"
        or source_audit.get("source_mechanics_independently_rederived") is not True
        or control.get("status")
        != "source_audited_revision_locked_before_submission"
        or control.get("revision_id") != "4f96d314e84de086"
        or control.get("independent_source_audit_artifacts_sha256")
        != sha256_file(source_audit_root / "ARTIFACTS.json")
        or source_reference.get("status") != "passed_exact_legacy_stream"
        or source_reference.get("legacy_stream_admission", {}).get("sha256")
        != "2f9cd9e853a9284c53884e6a551b1c7284795dd053f255d630aeeb114d1fa81f"
        or rsem_reference.get("status") != "passed"
        or rsem_reference.get("assembly") != "GRCh38.p14"
        or rsem_reference.get("annotation") != "GENCODE_v49"
        or rsem_reference.get("fasta_sha256")
        != "9489780d014865df158afc650e1f0ccc204b3a9878e147c1aae2ac982df26215"
        or rsem_reference.get("gtf_sha256")
        != "73bbbbd6eb2f114d1536f7cbf2339a653e1edc889b0cbe78992e92560b5aaba9"
        or activation.get("status") != "assays_ready_no_fit_evaluator_implementation_pending"
        or terminal.get("status") != "assays_ready_no_fit_gzip_contract_repaired"
        or terminal.get("clean_or_sealed_champion_eligible_by_this_receipt") is not False
        or crosswalk.get("status") != "passed_failure_aware_roundtrip"
        or crosswalk.get("source_reference_status") != "passed_exact_legacy_stream"
        or crosswalk.get("all_source_intervals_retained") is not True
        or rna.get("status") != "passed_label_free_raw_scale"
        or rna.get("post_quantification_normalization_or_feature_filter_applied") is not False
    ):
        raise ExternalCohortReadinessError("GSE105127 no-fit status differs")
    activation_sources = activation.get("source_artifacts")
    if activation_sources != {
        "plan": sha256_file(plan_root / "ARTIFACTS.json"),
        "reference": sha256_file(source_reference_root / "ARTIFACTS.json"),
        "crosswalk": sha256_file(crosswalk_root / "ARTIFACTS.json"),
        "rna": sha256_file(rna_root / "ARTIFACTS.json"),
    }:
        raise ExternalCohortReadinessError("GSE105127 activation source hashes differ")
    if (
        terminal.get("control_artifacts_sha256")
        != sha256_file(control_root / "ARTIFACTS.json")
        or terminal.get("independent_source_audit_artifacts_sha256")
        != sha256_file(source_audit_root / "ARTIFACTS.json")
        or terminal.get("activation_artifacts_sha256")
        != sha256_file(activation_root / "ARTIFACTS.json")
    ):
        raise ExternalCohortReadinessError("GSE105127 terminal source hashes differ")

    rna_fields, rna_rows = read_tsv(plan_root / "rna_rows.tsv")
    rrbs_fields, rrbs_rows = read_tsv(plan_root / "rrbs_rows.tsv")
    topology = validate_gse105127_plan_rows(
        rna_fields, rna_rows, rrbs_fields, rrbs_rows
    )
    row_axis_fields, row_axis = read_tsv(rna_root / "row_axis.tsv")
    if (
        row_axis_fields != GSE105127_ROW_AXIS_FIELDS
        or len(row_axis) != 57
        or [int(row["row_index"]) for row in row_axis] != list(range(57))
        or any(row["pairing_topology"] != "adjacent_section" for row in row_axis)
    ):
        raise ExternalCohortReadinessError("GSE105127 RNA row axis differs")
    model_fields, model_rows = read_tsv(activation_root / "model_input_manifest.tsv")
    if (
        model_fields != GSE105127_MODEL_INPUT_FIELDS
        or len(model_rows) != 57
        or {row["row_id"] for row in model_rows} != {row["row_id"] for row in rna_rows}
        or any(row["pairing_topology"] != "adjacent_section" for row in model_rows)
    ):
        raise ExternalCohortReadinessError("GSE105127 model-input topology differs")
    fold_fields, fold_rows = read_tsv(activation_root / "participant_folds.tsv")
    fold_counts = validate_gse105127_fold_rows(fold_fields, fold_rows)

    raw_rna_sources = json.loads(
        (rna_root / "source_receipts.json").read_text(encoding="utf-8")
    )
    if not isinstance(raw_rna_sources, list) or len(raw_rna_sources) != 57:
        raise ExternalCohortReadinessError("GSE105127 RNA source receipts differ")
    rna_source_by_id = {
        str(item.get("row_id")): item
        for item in raw_rna_sources
        if isinstance(item, dict)
    }
    if len(rna_source_by_id) != 57:
        raise ExternalCohortReadinessError("GSE105127 RNA source identity differs")
    for planned in rna_rows:
        row_id = planned["row_id"]
        root = campaign_root / "rna_quant" / "participants" / row_id
        require_frozen(root, artifact_class="gse105127_rna_rsem_v2")
        receipt = load_json(root / "receipt.json")
        require_false(
            receipt,
            ("labels_accessed", "fit_or_score_performed"),
            "GSE105127 RNA participant",
        )
        source = rna_source_by_id.get(row_id, {})
        if (
            receipt.get("status") != "passed_raw_scale"
            or receipt.get("participant_group_id") != planned["participant_group_id"]
            or receipt.get("zone") != planned["zone"]
            or receipt.get("pairing_topology") != "adjacent_section"
            or receipt.get("post_quantification_normalization_or_feature_filter_applied")
            is not False
            or source.get("source_artifacts_sha256")
            != sha256_file(root / "ARTIFACTS.json")
            or source.get("source_gene_results_sha256")
            != receipt.get("gene_results_sha256")
        ):
            raise ExternalCohortReadinessError("GSE105127 RNA participant differs")

    rrbs_by_id = {row["row_id"]: row for row in rrbs_rows}
    missingness = Counter()
    for row_id in sorted(rrbs_by_id):
        root = campaign_root / "rrbs_collapsed" / "participants" / row_id
        require_frozen(root, artifact_class="gse105127_collapsed_rrbs")
        receipt = load_json(root / "receipt.json")
        require_false(receipt, ("labels_accessed", "fit_or_score_performed", "missing_as_zero"), "GSE105127 RRBS")
        planned = rrbs_by_id[row_id]
        if (
            receipt.get("status") != "passed"
            or receipt.get("participant_group_id") != planned["participant_group_id"]
            or receipt.get("zone") != planned["zone"]
            or receipt.get("pairing_topology") != "adjacent_section"
            or int(receipt.get("coverage", 0)) <= 0
            or int(receipt.get("canonical_cpgs", 0)) <= 0
        ):
            raise ExternalCohortReadinessError("GSE105127 RRBS receipt differs")
        missingness["single_strand_cpgs"] += int(receipt["single_strand_cpgs"])
        missingness["paired_cpgs"] += int(receipt["paired_cpgs"])
    mapping_states = crosswalk.get("mapping_states")
    if (
        not isinstance(mapping_states, dict)
        or sum(int(value) for value in mapping_states.values())
        != int(crosswalk.get("unique_source_cpg_intervals", -1))
        or int(mapping_states.get("mapped_unique_cpg", 0)) <= 0
    ):
        raise ExternalCohortReadinessError("GSE105127 crosswalk missingness differs")

    expected_links = {"rsem_reference", "cpg_crosswalk"}
    expected_links.update(f"rrbs_collapse_bundle_{bundle:02d}" for bundle in range(8))
    expected_links.update(f"rna_quantification_bundle_{bundle:02d}" for bundle in range(8))
    if set(terminal.get("stage_links", {})) != expected_links:
        raise ExternalCohortReadinessError("GSE105127 terminal stage links differ")
    for name in sorted(expected_links):
        root = campaign_root / "stage_links" / name
        require_frozen(root, artifact_class="gse105127_gzip_stage_link_v2")
        link = load_json(root / "receipt.json")
        require_false(
            link,
            (
                "molecular_values_opened_by_linker",
                "labels_accessed",
                "outcomes_accessed",
                "fit_or_score_performed",
            ),
            "GSE105127 stage link",
        )
        if link.get("status") != "job_local_preflight_bound_to_completed_stage":
            raise ExternalCohortReadinessError("GSE105127 stage-link status differs")
        terminal_link = terminal["stage_links"][name]
        if (
            terminal_link.get("artifacts_sha256")
            != sha256_file(root / "ARTIFACTS.json")
            or terminal_link.get("slurm_job_id") != link.get("slurm_job_id")
            or int(terminal_link.get("preflight_file_count", -1))
            != int(link.get("preflight_file_count", -2))
        ):
            raise ExternalCohortReadinessError("GSE105127 terminal link hash differs")

    receipt = {
        "schema_version": "masld-bench-external-cohort-development-readiness-v1",
        "cohort": "GSE105127",
        "status": "assay_inputs_ready_development_modeling_blocked",
        "eligible_scope": "development_only_adjacent_section_rna_rrbs_input_mechanics",
        "topology": topology,
        "participant_outer_fold_counts": fold_counts,
        "explicit_missingness": {
            "rna_rows_observed": 57,
            "rrbs_missing_cpgs_masked_never_zero": True,
            "rrbs_single_strand_cpg_observations_retained": missingness[
                "single_strand_cpgs"
            ],
            "rrbs_paired_cpg_observations_retained": missingness["paired_cpgs"],
            "crosswalk_mapping_states": mapping_states,
        },
        "immutable_inputs": {
            "plan_artifacts_sha256": sha256_file(plan_root / "ARTIFACTS.json"),
            "source_audit_artifacts_sha256": sha256_file(source_audit_root / "ARTIFACTS.json"),
            "control_artifacts_sha256": sha256_file(control_root / "ARTIFACTS.json"),
            "source_reference_artifacts_sha256": sha256_file(
                source_reference_root / "ARTIFACTS.json"
            ),
            "rsem_reference_artifacts_sha256": sha256_file(
                rsem_reference_root / "ARTIFACTS.json"
            ),
            "activation_artifacts_sha256": sha256_file(activation_root / "ARTIFACTS.json"),
            "terminal_artifacts_sha256": sha256_file(terminal_root / "ARTIFACTS.json"),
            "crosswalk_artifacts_sha256": sha256_file(crosswalk_root / "ARTIFACTS.json"),
            "rna_artifacts_sha256": sha256_file(rna_root / "ARTIFACTS.json"),
        },
        "labels_read": False,
        "outcomes_read": False,
        "molecular_artifact_bytes_hashed_for_integrity": True,
        "expression_values_deserialized_by_readiness_auditor": False,
        "fit_or_score_performed": False,
        "modeling_authorized_by_this_receipt": False,
        "sealed_evaluation_authorized_by_this_receipt": False,
        "clean_or_sealed_champion_eligible_by_this_receipt": False,
        "remaining_blockers": [
            "implement_and_fixture_test_the_frozen_coverage_aware_evaluator",
            "freeze_model_receptive_field_specific_genomic_blocks_and_buffers",
            "complete_model_specific_derivative_weight_review",
        ],
    }
    output.mkdir(mode=0o750)
    write_json_exclusive(output / "development_readiness.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "external_cohort_development_readiness_disposition",
            "cohort": "GSE105127",
            "eligible_scope": receipt["eligible_scope"],
            "modeling_authorized": False,
            "clean_or_sealed_champion_eligible": False,
            "status": receipt["status"],
        },
    )
    verify_frozen_tree(output)
    return receipt


def audit_gse268273(
    *,
    plan_root: Path,
    prior_model_input_root: Path,
    consolidation_root: Path,
    independent_audit_root: Path,
    output: Path,
    expected_plan_artifacts_sha256: str,
    expected_prior_model_input_artifacts_sha256: str,
) -> dict[str, Any]:
    if output.exists():
        raise ExternalCohortReadinessError("refusing to overwrite GSE268273 disposition")
    require_frozen(
        plan_root,
        artifact_class="gse268273_raw_campaign_plan",
        expected_sha256=expected_plan_artifacts_sha256,
    )
    require_frozen(
        prior_model_input_root,
        artifact_class="gse268273_fibrosis_ood_model_input",
        expected_sha256=expected_prior_model_input_artifacts_sha256,
    )
    require_frozen(
        consolidation_root,
        artifact_class="gse268273_raw_rsem_consolidation_execution",
    )
    source_root = consolidation_root / "model_input"
    require_frozen(
        source_root, artifact_class="gse268273_outcome_free_raw_rsem_model_input"
    )
    require_frozen(
        independent_audit_root,
        artifact_class="gse268273_independent_raw_source_integrity_audit_execution",
    )
    audit_root = independent_audit_root / "audit"
    require_frozen(
        audit_root, artifact_class="gse268273_independent_raw_source_integrity_audit"
    )

    plan_fields, plan_rows = read_tsv(plan_root / "participants.tsv")
    topology = validate_gse268273_topology(plan_fields, plan_rows)
    prior_fields, prior_rows = read_tsv(prior_model_input_root / "participant_axis.tsv")
    source_fields, source_rows = read_tsv(source_root / "participant_axis.tsv")
    if (
        prior_fields != GSE268273_PARTICIPANT_FIELDS[:-1]
        or source_fields != GSE268273_PARTICIPANT_FIELDS
        or len(prior_rows) != 109
        or len(source_rows) != 109
        or [row["row_id"] for row in prior_rows] != [row["row_id"] for row in source_rows]
    ):
        raise ExternalCohortReadinessError("GSE268273 participant axis differs")
    plan_by_id = {row["row_id"]: row for row in plan_rows}
    if set(plan_by_id) != {row["row_id"] for row in source_rows}:
        raise ExternalCohortReadinessError("GSE268273 plan-to-source join differs")
    for prior, source in zip(prior_rows, source_rows, strict=True):
        expected = dict(prior)
        expected["rna_observation_state"] = "observed"
        expected["quantification_measurement"] = "RSEM_expected_count_raw_count_scale"
        if source != expected:
            raise ExternalCohortReadinessError("GSE268273 source state mutation differs")
        planned = plan_by_id[source["row_id"]]
        if (
            source["cohort_family_id"] != "gse268273_imid_masld"
            or int(source["raw_technical_run_count"])
            != int(planned["technical_runs"])
            or int(source["raw_fastq_bytes"]) != int(planned["fastq_bytes"])
        ):
            raise ExternalCohortReadinessError("GSE268273 participant-run join differs")

    source_receipt = load_json(source_root / "receipt.json")
    independent = load_json(
        audit_root / "independent_source_integrity_audit.json"
    )
    require_false(
        source_receipt,
        (
            "normalization_applied",
            "missing_encoded_as_zero",
            "labels_included",
            "clinical_covariates_included",
            "outcomes_accessed",
            "fit_or_score_performed",
            "source_global_voom_used",
            "processed_differential_expression_used",
        ),
        "GSE268273 source",
    )
    require_false(
        independent,
        (
            "endpoint_decision_encoded",
            "model_decision_encoded",
            "transform_decision_encoded",
            "normalization_applied",
            "feature_selection_performed",
            "labels_read",
            "clinical_covariates_read",
            "outcomes_read",
            "predictions_read",
            "fit_or_score_performed",
            "modeling_or_scoring_authorized_by_this_receipt",
            "sealed_evaluation_authorized_by_this_receipt",
            "clean_or_sealed_champion_eligible_by_this_receipt",
        ),
        "GSE268273 independent audit",
    )
    if (
        source_receipt.get("status")
        != "passed_outcome_free_quantification_training_transform_blocked"
        or source_receipt.get("rna_observed_mask_all_true") is not True
        or source_receipt.get("measured_zeros_are_values") is not True
        or independent.get("status") != "passed_source_integrity_mechanics_only"
        or independent.get("source_integrity_mechanics_independently_rederived") is not True
        or independent.get("source_integrity_mechanics_eligible") is not True
        or independent.get("eligible_scope") != "source_integrity_mechanics_only"
    ):
        raise ExternalCohortReadinessError("GSE268273 source disposition differs")
    observed_mask = np.load(source_root / "rna_observed_mask.npy", allow_pickle=False)
    if (
        observed_mask.shape != (109,)
        or observed_mask.dtype != np.bool_
        or not bool(np.all(observed_mask))
    ):
        raise ExternalCohortReadinessError("GSE268273 missingness mask differs")

    receipt = {
        "schema_version": "masld-bench-external-cohort-development-readiness-v1",
        "cohort": "GSE268273",
        "status": "raw_source_mechanics_ready_development_training_blocked",
        "eligible_scope": "development_only_outcome_free_raw_rsem_source_mechanics",
        "topology": topology,
        "explicit_missingness": {
            "participant_rna_rows_observed": 109,
            "rna_observed_mask_all_true": True,
            "measured_zero_is_value": True,
            "missing_encoded_as_zero": False,
            "unavailable_clinical_or_histology_fields_not_read": True,
        },
        "immutable_inputs": {
            "plan_artifacts_sha256": sha256_file(plan_root / "ARTIFACTS.json"),
            "prior_model_input_artifacts_sha256": sha256_file(
                prior_model_input_root / "ARTIFACTS.json"
            ),
            "consolidation_artifacts_sha256": sha256_file(
                consolidation_root / "ARTIFACTS.json"
            ),
            "source_artifacts_sha256": sha256_file(source_root / "ARTIFACTS.json"),
            "independent_audit_execution_artifacts_sha256": sha256_file(
                independent_audit_root / "ARTIFACTS.json"
            ),
            "independent_audit_artifacts_sha256": sha256_file(
                audit_root / "ARTIFACTS.json"
            ),
        },
        "labels_read": False,
        "outcomes_read": False,
        "molecular_artifact_bytes_hashed_for_integrity": True,
        "expression_values_deserialized_by_readiness_auditor": False,
        "fit_or_score_performed": False,
        "modeling_authorized_by_this_receipt": False,
        "sealed_evaluation_authorized_by_this_receipt": False,
        "clean_or_sealed_champion_eligible_by_this_receipt": False,
        "remaining_blockers": [
            "freeze_training_cohort_only_feature_filter_normalization_and_transform",
            "freeze_model_family_preprocessing_calibration_thresholds_and_ensemble_without_GSE268273",
            "freeze_separate_outcome_evaluator_and_prediction_commit_firewall",
            "retain_project_exposed_external_development_not_sealed_status",
        ],
    }
    output.mkdir(mode=0o750)
    write_json_exclusive(output / "development_readiness.json", receipt)
    freeze_tree(
        output,
        {
            "artifact_class": "external_cohort_development_readiness_disposition",
            "cohort": "GSE268273",
            "eligible_scope": receipt["eligible_scope"],
            "modeling_authorized": False,
            "clean_or_sealed_champion_eligible": False,
            "status": receipt["status"],
        },
    )
    verify_frozen_tree(output)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="cohort", required=True)
    gse105127 = subparsers.add_parser("gse105127")
    gse105127.add_argument("--plan-root", required=True, type=Path)
    gse105127.add_argument("--campaign-root", required=True, type=Path)
    gse105127.add_argument("--source-audit-root", required=True, type=Path)
    gse105127.add_argument("--output", required=True, type=Path)
    gse105127.add_argument("--expected-plan-artifacts-sha256", required=True)
    gse105127.add_argument("--expected-source-audit-artifacts-sha256", required=True)
    gse268273 = subparsers.add_parser("gse268273")
    gse268273.add_argument("--plan-root", required=True, type=Path)
    gse268273.add_argument("--prior-model-input-root", required=True, type=Path)
    gse268273.add_argument("--consolidation-root", required=True, type=Path)
    gse268273.add_argument("--independent-audit-root", required=True, type=Path)
    gse268273.add_argument("--output", required=True, type=Path)
    gse268273.add_argument("--expected-plan-artifacts-sha256", required=True)
    gse268273.add_argument(
        "--expected-prior-model-input-artifacts-sha256", required=True
    )
    arguments = parser.parse_args()
    values = vars(arguments)
    cohort = values.pop("cohort")
    if cohort == "gse105127":
        receipt = audit_gse105127(**values)
    else:
        receipt = audit_gse268273(**values)
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
