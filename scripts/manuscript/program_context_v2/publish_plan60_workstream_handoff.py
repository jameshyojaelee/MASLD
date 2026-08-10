#!/usr/bin/env python3
"""Publish exact, signed, candidate-only Plan 60 workstream handoffs.

The publisher never copies, edits, or promotes scientific artifacts. It hashes
the exact terminal artifact set required by the release coordinator, writes one
owner manifest inside the workstream root, and records a detached non-promotion
attestation. Existing handoffs are never overwritten.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import tempfile
from pathlib import Path

from release_common import (
    ARTIFACT_FIELDS,
    CANDIDATE_ID,
    ReleaseContractError,
    allowed_workstream_roots,
    atomic_write_json,
    atomic_write_tsv,
    clean_relative_path,
    parse_nonnegative_int,
    project_relative,
    read_tsv_exact,
    require_sha256,
    require_within_any,
    resolve_project_path,
    sha256_file,
    validate_artifact_manifest,
)


HANDOFF_SIGNATURE_VERSION = "plan60_workstream_handoff_v1"
WORKSTREAM_LANES = {
    "PLAN13": "SP-INT",
    "PLAN20": "HS-V2",
    "PLAN30": "GEN",
    "PLAN40": "HLF",
    "PLAN50": "PASS",
}
PLAN50_PAYLOAD_MANIFEST_FIELDS = (
    "relative_path",
    "sha256",
    "bytes",
    "producer",
    "environment",
    "producer_sha256",
    "upstream_artifact_sha256",
    "release_status",
)
PLAN50_TERMINAL_CHAIN = (
    ("passport_terminal", "PASS06_VALIDATED", "terminal_verdict"),
    (
        "passport_release_manifest",
        "passport_release_manifest.tsv",
        "release_manifest",
    ),
    (
        "passport_validation",
        "passport_validation_report.tsv",
        "workstream_validation",
    ),
    (
        "passport_terminal_provenance",
        "passport_terminal_provenance.tsv",
        "terminal_provenance",
    ),
    (
        "passport_inner_plan60_handoff",
        "passport_plan60_handoff.tsv",
        "inner_plan60_handoff",
    ),
)


def plan50_full_bundle_contract(
    root: Path,
) -> tuple[tuple[str, str, str, str], ...]:
    """Return every sealed Plan50 payload plus its five-file terminal chain."""

    manifest = root / "passport_release_manifest.tsv"
    payload_rows = read_tsv_exact(manifest, PLAN50_PAYLOAD_MANIFEST_FIELDS)
    if not payload_rows:
        raise ReleaseContractError("Plan50 sealed payload manifest is empty")
    specifications = [
        (artifact_id, relative, role, relative)
        for artifact_id, relative, role in PLAN50_TERMINAL_CHAIN
    ]
    seen_paths = {row[1] for row in specifications}
    seen_ids = {row[0] for row in specifications}
    for row in sorted(payload_rows, key=lambda value: value["relative_path"]):
        relative = clean_relative_path(
            row["relative_path"], "Plan50 payload relative_path"
        )
        parts = Path(relative).parts
        if (
            relative in seen_paths
            or "fixtures" in parts
            or "superseded_bundles" in parts
            or relative.startswith("plan60_terminal_artifacts")
        ):
            raise ReleaseContractError(
                f"unsafe or duplicate Plan50 payload path: {relative}"
            )
        source = root / relative
        expected_hash = require_sha256(
            row["sha256"], f"Plan50 payload {relative} SHA256"
        )
        expected_bytes = parse_nonnegative_int(
            row["bytes"], f"Plan50 payload {relative} bytes"
        )
        if (
            not source.is_file()
            or source.is_symlink()
            or source.stat().st_size != expected_bytes
            or sha256_file(source) != expected_hash
        ):
            raise ReleaseContractError(
                f"Plan50 payload manifest/file drift: {relative}"
            )
        if row["release_status"] != "validated_candidate_handoff":
            raise ReleaseContractError(f"Plan50 payload is not terminal: {relative}")
        artifact_id = (
            "passport_payload_"
            + hashlib.sha256(relative.encode("utf-8")).hexdigest()[:16]
        )
        if artifact_id in seen_ids:
            raise ReleaseContractError(
                f"Plan50 payload artifact ID collision: {relative}"
            )
        seen_ids.add(artifact_id)
        seen_paths.add(relative)
        specifications.append(
            (artifact_id, relative, "sealed_passport_payload", relative)
        )
    return tuple(specifications)


def validate_handoff_rows_in_memory(
    project: Path,
    workstream_id: str,
    rows: list[dict[str, object]],
) -> None:
    """Apply the manifest contract before either handoff file is created."""

    if not rows:
        raise ReleaseContractError(f"empty {workstream_id} handoff row set")
    lane = WORKSTREAM_LANES[workstream_id]
    roots = allowed_workstream_roots(project, workstream_id)
    artifact_ids: set[str] = set()
    snapshots: set[str] = set()
    terminal_count = 0
    for row in rows:
        if tuple(row) != ARTIFACT_FIELDS:
            raise ReleaseContractError(
                f"{workstream_id} in-memory handoff schema drift"
            )
        artifact_id = str(row["artifact_id"]).strip()
        if (
            not artifact_id
            or artifact_id in artifact_ids
            or artifact_id == "terminal_artifact_manifest"
        ):
            raise ReleaseContractError(
                f"blank, duplicate, or reserved {workstream_id} artifact_id: {artifact_id!r}"
            )
        artifact_ids.add(artifact_id)
        source = resolve_project_path(
            project, str(row["source_path"]), f"{workstream_id} {artifact_id} source"
        )
        require_within_any(source, roots, f"{workstream_id} {artifact_id} source")
        if not source.is_file() or source.is_symlink():
            raise ReleaseContractError(
                f"artifact source must be a regular non-symlink file: {source}"
            )
        digest = require_sha256(
            str(row["source_sha256"]),
            f"{workstream_id} {artifact_id} source_sha256",
        )
        size = parse_nonnegative_int(
            str(row["source_bytes"]), f"{workstream_id} {artifact_id} source_bytes"
        )
        if source.stat().st_size != size or sha256_file(source) != digest:
            raise ReleaseContractError(
                f"artifact source hash/byte mismatch: {workstream_id} {artifact_id}"
            )
        snapshot = clean_relative_path(
            str(row["snapshot_relpath"]),
            f"{workstream_id} {artifact_id} snapshot_relpath",
        )
        if not snapshot.startswith(f"inputs/{lane}/") or snapshot in snapshots:
            raise ReleaseContractError(
                f"invalid or duplicate {workstream_id} snapshot path: {snapshot}"
            )
        snapshots.add(snapshot)
        downstream = clean_relative_path(
            str(row["downstream_read_path"]),
            f"{workstream_id} {artifact_id} downstream_read_path",
        )
        if downstream != snapshot:
            raise ReleaseContractError(
                f"{workstream_id} downstream path must equal its snapshot path"
            )
        role = str(row["artifact_role"]).strip()
        if not role:
            raise ReleaseContractError(
                f"blank artifact role: {workstream_id} {artifact_id}"
            )
        terminal_count += role == "terminal_verdict"
    if terminal_count != 1:
        raise ReleaseContractError(
            f"{workstream_id} requires exactly one terminal_verdict artifact; observed {terminal_count}"
        )


def workstream_contract(
    project: Path,
) -> dict[str, tuple[Path, tuple[tuple[str, str, str, str], ...]]]:
    multimodal = (
        project / "Analysis/Multimodal_Program_Projection/candidates" / CANDIDATE_ID
    )
    plan13 = multimodal / "spatial_context_semantic_v2_2026-08-08"
    plan20 = multimodal / "hotspot"
    plan30 = multimodal / "genetics_context"
    plan40 = multimodal / "myojin_hlf"
    plan50 = project / "RNA-seq/results/evidence_passports/candidates" / CANDIDATE_ID
    contracts = {
        "PLAN13": (
            plan13,
            (
                (
                    "plan13_terminal",
                    "SEMANTIC_V2_READY",
                    "terminal_verdict",
                    "SEMANTIC_V2_READY",
                ),
                (
                    "plan13_release_manifest",
                    "semantic_v2_release_manifest.tsv",
                    "release_manifest",
                    "semantic_v2_release_manifest.tsv",
                ),
                (
                    "plan13_figure4_matrix",
                    "final_integration/figure4_program_matrix.tsv",
                    "figure_source",
                    "final_integration/figure4_program_matrix.tsv",
                ),
                (
                    "plan13_dataset_verdict",
                    "final_integration/figure4_dataset_verdict.tsv",
                    "dataset_verdict",
                    "final_integration/figure4_dataset_verdict.tsv",
                ),
                (
                    "plan13_validation",
                    "validation_report.tsv",
                    "workstream_validation",
                    "validation_report.tsv",
                ),
                (
                    "plan13_producers",
                    "final_integration_code_freeze.tsv",
                    "producer_manifest",
                    "final_integration_code_freeze.tsv",
                ),
                (
                    "plan13_final_bundle",
                    "final_integration/bundle_manifest.tsv",
                    "final_bundle_manifest",
                    "final_integration/bundle_manifest.tsv",
                ),
                (
                    "plan13_final_validation",
                    "final_integration/validation_report.tsv",
                    "final_validation",
                    "final_integration/validation_report.tsv",
                ),
                (
                    "plan13_input_seal",
                    "final_integration/input_seal_manifest.tsv",
                    "input_seal_manifest",
                    "final_integration/input_seal_manifest.tsv",
                ),
                (
                    "plan13_source_baseline",
                    "semantic_v2_source_baseline.tsv",
                    "source_baseline",
                    "semantic_v2_source_baseline.tsv",
                ),
                (
                    "plan13_v1_preservation",
                    "v1_preservation_final.tsv",
                    "v1_preservation",
                    "v1_preservation_final.tsv",
                ),
                (
                    "plan13_adapter_bundle",
                    "real_adapter_bundle_manifest.tsv",
                    "adapter_bundle_manifest",
                    "real_adapter_bundle_manifest.tsv",
                ),
                (
                    "plan13_adapter_validation",
                    "real_adapter_validation.tsv",
                    "adapter_validation",
                    "real_adapter_validation.tsv",
                ),
            ),
        ),
        "PLAN20": (
            plan20,
            (
                (
                    "hotspot_semantic_terminal",
                    "SEMANTIC_ADJUDICATION_READY",
                    "terminal_verdict",
                    "SEMANTIC_ADJUDICATION_READY",
                ),
                ("hotspot_ready", "READY", "ready_seal", "READY"),
                (
                    "hotspot_registry",
                    "program_registry_v2.tsv",
                    "program_registry",
                    "program_registry_v2.tsv",
                ),
                (
                    "hotspot_membership",
                    "program_membership_v2.tsv",
                    "program_membership",
                    "program_membership_v2.tsv",
                ),
                (
                    "hotspot_semantics",
                    "program_registry_v2_semantic_adjudication.tsv",
                    "semantic_adjudication",
                    "program_registry_v2_semantic_adjudication.tsv",
                ),
                (
                    "hotspot_figure2",
                    "fig2_program_source.tsv",
                    "figure_source",
                    "fig2_program_source.tsv",
                ),
                (
                    "composition_ready",
                    "COMPOSITION_QC_READY",
                    "composition_ready_seal",
                    "COMPOSITION_QC_READY",
                ),
                (
                    "composition_results",
                    "composition_sample_qc_v2.tsv",
                    "sample_composition",
                    "composition_sample_qc_v2.tsv",
                ),
                (
                    "composition_testability",
                    "composition_sample_qc_testability.tsv",
                    "sample_composition_testability",
                    "composition_sample_qc_testability.tsv",
                ),
                (
                    "composition_audit",
                    "composition_sample_qc_audit.tsv",
                    "sample_composition_audit",
                    "composition_sample_qc_audit.tsv",
                ),
                (
                    "composition_validation",
                    "composition_sample_qc_validation.tsv",
                    "sample_composition_validation",
                    "composition_sample_qc_validation.tsv",
                ),
                (
                    "composition_inputs",
                    "composition_sample_qc_input_manifest.tsv",
                    "sample_composition_input_manifest",
                    "composition_sample_qc_input_manifest.tsv",
                ),
                (
                    "composition_sensitivity",
                    "composition_sample_qc_sensitivity.tsv",
                    "sample_composition_sensitivity",
                    "composition_sample_qc_sensitivity.tsv",
                ),
                (
                    "hotspot_analysis_spec",
                    "analysis_specification.tsv",
                    "analysis_specification",
                    "analysis_specification.tsv",
                ),
                (
                    "hotspot_cohort_lodo",
                    "cohort_and_lodo_effects.tsv",
                    "cohort_lodo_effects",
                    "cohort_and_lodo_effects.tsv",
                ),
                (
                    "hotspot_tested_universe",
                    "tested_universe.tsv",
                    "tested_universe",
                    "tested_universe.tsv",
                ),
                (
                    "hotspot_release_manifest",
                    "release_manifest.tsv",
                    "workstream_release_manifest",
                    "release_manifest.tsv",
                ),
                (
                    "hotspot_validation",
                    "validation_status.tsv",
                    "workstream_validation",
                    "validation_status.tsv",
                ),
                (
                    "hotspot_stage_design",
                    "stage_dataset_design_audit.tsv",
                    "stage_design_audit",
                    "stage_dataset_design_audit.tsv",
                ),
                (
                    "hotspot_environment",
                    "environment_record.tsv",
                    "environment_record",
                    "environment_record.tsv",
                ),
                (
                    "hotspot_inputs",
                    "input_manifest.tsv",
                    "workstream_input_manifest",
                    "input_manifest.tsv",
                ),
                (
                    "hotspot_v1_preservation",
                    "v1_preservation.tsv",
                    "v1_preservation",
                    "v1_preservation.tsv",
                ),
                (
                    "hotspot_semantic_manifest",
                    "semantic_adjudication_manifest.tsv",
                    "semantic_adjudication_manifest",
                    "semantic_adjudication_manifest.tsv",
                ),
                (
                    "hotspot_semantic_validation",
                    "semantic_adjudication_validation.tsv",
                    "semantic_adjudication_validation",
                    "semantic_adjudication_validation.tsv",
                ),
                (
                    "hotspot_documented_fstage",
                    "documented_fstage_sensitivity.tsv",
                    "documented_fstage_sensitivity",
                    "documented_fstage_sensitivity.tsv",
                ),
                (
                    "nmf_continuous_ready",
                    "nmf_continuous_supplement/NMF_CONTINUOUS_SUPPLEMENT_READY",
                    "supplement_ready_seal",
                    "nmf_continuous/NMF_CONTINUOUS_SUPPLEMENT_READY",
                ),
                (
                    "nmf_continuous_supplement",
                    "nmf_continuous_supplement/nmf_continuous_supplement.tsv",
                    "figure_source_adapter",
                    "nmf_continuous_supplement.tsv",
                ),
                (
                    "nmf_continuous_loadings",
                    "nmf_continuous_supplement/nmf_continuous_loadings.tsv",
                    "supplement_loadings",
                    "nmf_continuous/nmf_continuous_loadings.tsv",
                ),
                (
                    "nmf_continuous_source_manifest",
                    "nmf_continuous_supplement/source_manifest.tsv",
                    "supplement_source_manifest",
                    "nmf_continuous/source_manifest.tsv",
                ),
                (
                    "nmf_continuous_validation",
                    "nmf_continuous_supplement/validation_report.tsv",
                    "supplement_validation",
                    "nmf_continuous/validation_report.tsv",
                ),
                (
                    "nmf_continuous_producers",
                    "nmf_continuous_supplement/producer_manifest.tsv",
                    "supplement_producer_manifest",
                    "nmf_continuous/producer_manifest.tsv",
                ),
                (
                    "nmf_k4_loadings_source",
                    "nmf_continuous_supplement/source_inputs/nmf_assignments_k4_pre_k6restore.csv",
                    "frozen_source_input",
                    "nmf_continuous/source_inputs/nmf_assignments_k4_pre_k6restore.csv",
                ),
                (
                    "nmf_k6_loadings_source",
                    "nmf_continuous_supplement/source_inputs/nmf_assignments_k6.csv",
                    "frozen_source_input",
                    "nmf_continuous/source_inputs/nmf_assignments_k6.csv",
                ),
                (
                    "nmf_k4_labels_source",
                    "nmf_continuous_supplement/source_inputs/program_labels_k4_pre_k6restore.csv",
                    "frozen_source_input",
                    "nmf_continuous/source_inputs/program_labels_k4_pre_k6restore.csv",
                ),
                (
                    "nmf_k6_labels_source",
                    "nmf_continuous_supplement/source_inputs/program_labels_k6.csv",
                    "frozen_source_input",
                    "nmf_continuous/source_inputs/program_labels_k6.csv",
                ),
                (
                    "nmf_three_seed_metrics_source",
                    "nmf_continuous_supplement/source_inputs/per_seed_metric_summary.csv",
                    "frozen_source_input",
                    "nmf_continuous/source_inputs/per_seed_metric_summary.csv",
                ),
                (
                    "nmf_three_seed_stability_source",
                    "nmf_continuous_supplement/source_inputs/cross_seed_stability_summary.csv",
                    "frozen_source_input",
                    "nmf_continuous/source_inputs/cross_seed_stability_summary.csv",
                ),
            ),
        ),
        "PLAN30": (
            plan30,
            (
                (
                    "gen_terminal",
                    "GEN_TERMINAL_CLOSURE_READY",
                    "terminal_verdict",
                    "GEN_TERMINAL_CLOSURE_READY",
                ),
                (
                    "gen_closure",
                    "terminal_closure.tsv",
                    "terminal_closure",
                    "terminal_closure.tsv",
                ),
                (
                    "gen_classes",
                    "frozen_evidence_classes.tsv",
                    "frozen_evidence_classes",
                    "frozen_evidence_classes.tsv",
                ),
                (
                    "gen_phenotypes",
                    "phenotype_registry.tsv",
                    "phenotype_registry",
                    "phenotype_registry.tsv",
                ),
                (
                    "gen_power_interface",
                    "power_stratified_interface.tsv",
                    "power_interface",
                    "power_stratified_interface.tsv",
                ),
                (
                    "gen_inputs",
                    "input_manifest.tsv",
                    "workstream_input_manifest",
                    "input_manifest.tsv",
                ),
                (
                    "gen_producers",
                    "producer_manifest.tsv",
                    "producer_manifest",
                    "producer_manifest.tsv",
                ),
                (
                    "gen_validation",
                    "validation_report.tsv",
                    "workstream_validation",
                    "validation_report.tsv",
                ),
                (
                    "gen_source_provenance",
                    "source_provenance.tsv",
                    "source_provenance",
                    "source_provenance.tsv",
                ),
                (
                    "gen_observability",
                    "gene_observability.tsv",
                    "gene_observability",
                    "gene_observability.tsv",
                ),
                (
                    "gen_tested_universe",
                    "source_tested_gene_universe_status.tsv",
                    "tested_universe_status",
                    "source_tested_gene_universe_status.tsv",
                ),
                (
                    "gen_power_balance",
                    "power_match_balance.tsv",
                    "power_match_balance",
                    "power_match_balance.tsv",
                ),
                (
                    "gen_sensitivity",
                    "sensitivity.tsv",
                    "sensitivity",
                    "sensitivity.tsv",
                ),
                (
                    "gen_closure_manifest",
                    "terminal_closure_manifest.tsv",
                    "terminal_closure_manifest",
                    "terminal_closure_manifest.tsv",
                ),
                (
                    "gen_closure_validation",
                    "terminal_closure_validation.tsv",
                    "terminal_closure_validation",
                    "terminal_closure_validation.tsv",
                ),
            ),
        ),
        "PLAN40": (
            plan40,
            (
                (
                    "myojin_terminal",
                    "PHASE_C_VALIDATED",
                    "terminal_verdict",
                    "PHASE_C_VALIDATED",
                ),
                (
                    "myojin_verdict",
                    "fig5_verdict.tsv",
                    "figure_verdict",
                    "fig5_verdict.tsv",
                ),
                (
                    "myojin_classes",
                    "class_effects.tsv",
                    "class_effects",
                    "class_effects.tsv",
                ),
                (
                    "myojin_programs",
                    "program_effects.tsv",
                    "program_effects",
                    "program_effects.tsv",
                ),
                (
                    "myojin_blind_spec",
                    "blind_analysis_spec.yaml",
                    "blind_analysis_specification",
                    "blind_analysis_spec.yaml",
                ),
                (
                    "myojin_predictions",
                    "prediction_manifest.tsv",
                    "outcome_blind_prediction_manifest",
                    "prediction_manifest.tsv",
                ),
                (
                    "myojin_spec_hash",
                    "specification_sha256.txt",
                    "specification_hash",
                    "specification_sha256.txt",
                ),
                (
                    "myojin_code_environment",
                    "phase_c_code_manifest.tsv",
                    "producer_environment_manifest",
                    "phase_c_code_manifest.tsv",
                ),
                (
                    "myojin_release_manifest",
                    "phase_c_release_manifest.tsv",
                    "workstream_release_manifest",
                    "phase_c_release_manifest.tsv",
                ),
                (
                    "myojin_validation",
                    "phase_c_validation_status.tsv",
                    "workstream_validation",
                    "phase_c_validation_status.tsv",
                ),
                (
                    "myojin_preflight",
                    "phase_c_preflight_audit.tsv",
                    "preflight_audit",
                    "phase_c_preflight_audit.tsv",
                ),
                (
                    "myojin_testability",
                    "program_testability.tsv",
                    "program_testability",
                    "program_testability.tsv",
                ),
                (
                    "myojin_tested_universe",
                    "screen_gene_universe.tsv",
                    "tested_universe",
                    "screen_gene_universe.tsv",
                ),
                (
                    "myojin_screen_results",
                    "screen_results.tsv",
                    "screen_results",
                    "screen_results.tsv",
                ),
                (
                    "myojin_sensitivity",
                    "sensitivity.tsv",
                    "sensitivity",
                    "sensitivity.tsv",
                ),
                (
                    "myojin_permutation_audit",
                    "permutation_audit.tsv",
                    "permutation_audit",
                    "permutation_audit.tsv",
                ),
                (
                    "myojin_matched_null_audit",
                    "matched_null_audit.tsv",
                    "matched_null_audit",
                    "matched_null_audit.tsv",
                ),
                (
                    "myojin_gene_mapping",
                    "gene_mapping_audit.tsv",
                    "gene_mapping_audit",
                    "gene_mapping_audit.tsv",
                ),
                (
                    "myojin_covariate_coverage",
                    "covariate_coverage_audit.tsv",
                    "covariate_coverage_audit",
                    "covariate_coverage_audit.tsv",
                ),
                (
                    "myojin_known_hit_exclusion",
                    "known_hit_exclusion.tsv",
                    "known_hit_exclusion",
                    "known_hit_exclusion.tsv",
                ),
                (
                    "myojin_source_manifest",
                    "source_manifest.tsv",
                    "workstream_source_manifest",
                    "source_manifest.tsv",
                ),
            ),
        ),
        "PLAN50": (
            plan50,
            (
                (
                    "passport_terminal",
                    "PASS06_VALIDATED",
                    "terminal_verdict",
                    "PASS06_VALIDATED",
                ),
                (
                    "passport_release_manifest",
                    "passport_release_manifest.tsv",
                    "release_manifest",
                    "passport_release_manifest.tsv",
                ),
                (
                    "passport_gene_index",
                    "passport_gene_index.parquet",
                    "gene_index",
                    "passport_gene_index.parquet",
                ),
                (
                    "passport_gene_evidence",
                    "passport_evidence_long.parquet",
                    "gene_evidence",
                    "passport_evidence_long.parquet",
                ),
                (
                    "passport_gene_coverage",
                    "passport_coverage_long.parquet",
                    "gene_coverage",
                    "passport_coverage_long.parquet",
                ),
                (
                    "passport_domain_coverage",
                    "passport_domain_coverage.tsv",
                    "domain_coverage",
                    "passport_domain_coverage.tsv",
                ),
                (
                    "passport_assay_status",
                    "passport_assay_status.tsv",
                    "native_grain_assay_status",
                    "passport_assay_status.tsv",
                ),
                (
                    "passport_program_context",
                    "passport_program_context.parquet",
                    "program_context",
                    "passport_program_context.parquet",
                ),
                (
                    "passport_next_experiment",
                    "passport_next_experiment.tsv",
                    "next_experiment",
                    "passport_next_experiment.tsv",
                ),
                (
                    "passport_input_selection",
                    "passport_input_selection.tsv",
                    "input_selection",
                    "passport_input_selection.tsv",
                ),
                (
                    "passport_input_selection_signature",
                    "passport_input_selection.signature.json",
                    "input_selection_signature",
                    "passport_input_selection.signature.json",
                ),
                (
                    "passport_manual_acceptance",
                    "passport_manual_acceptance.tsv",
                    "manual_acceptance",
                    "passport_manual_acceptance.tsv",
                ),
                (
                    "passport_manual_acceptance_signature",
                    "passport_manual_acceptance.signature.json",
                    "manual_acceptance_signature",
                    "passport_manual_acceptance.signature.json",
                ),
                (
                    "passport_validation",
                    "passport_validation_report.tsv",
                    "workstream_validation",
                    "passport_validation_report.tsv",
                ),
                (
                    "passport_terminal_provenance",
                    "passport_terminal_provenance.tsv",
                    "terminal_provenance",
                    "passport_terminal_provenance.tsv",
                ),
                (
                    "passport_inner_plan60_handoff",
                    "passport_plan60_handoff.tsv",
                    "inner_plan60_handoff",
                    "passport_plan60_handoff.tsv",
                ),
                (
                    "passport_source_nodes",
                    "passport_source_nodes.tsv",
                    "source_nodes",
                    "passport_source_nodes.tsv",
                ),
                (
                    "passport_source_edges",
                    "passport_source_edges.tsv",
                    "source_edges",
                    "passport_source_edges.tsv",
                ),
                (
                    "passport_rulebook",
                    "passport_rulebook.tsv",
                    "rulebook",
                    "passport_rulebook.tsv",
                ),
                (
                    "passport_vocabularies",
                    "passport_controlled_vocabularies.tsv",
                    "controlled_vocabularies",
                    "passport_controlled_vocabularies.tsv",
                ),
                (
                    "passport_dictionary",
                    "passport_data_dictionary.tsv",
                    "data_dictionary",
                    "passport_data_dictionary.tsv",
                ),
                (
                    "passport_gate",
                    "passport_gate_status.tsv",
                    "gate_status",
                    "passport_gate_status.tsv",
                ),
                (
                    "passport_program_index",
                    "passport_program_index.parquet",
                    "program_index",
                    "passport_program_index.parquet",
                ),
                (
                    "passport_program_membership",
                    "passport_program_membership.parquet",
                    "program_membership",
                    "passport_program_membership.parquet",
                ),
                (
                    "passport_accepted_input_audit",
                    "accepted_input_audit.tsv",
                    "accepted_input_audit",
                    "accepted_input_audit.tsv",
                ),
                (
                    "passport_adapter_audit",
                    "passport_adapter_audit.tsv",
                    "adapter_audit",
                    "passport_adapter_audit.tsv",
                ),
                (
                    "passport_ui_index",
                    "portal_candidate/index.html",
                    "local_review_ui",
                    "portal_candidate/index.html",
                ),
                (
                    "passport_ui_contract",
                    "portal_candidate/ui_contract.json",
                    "ui_contract",
                    "portal_candidate/ui_contract.json",
                ),
                (
                    "passport_ui_source_manifest",
                    "portal_candidate/ui_source_manifest.tsv",
                    "ui_source_manifest",
                    "portal_candidate/ui_source_manifest.tsv",
                ),
                (
                    "passport_ui_review_overview",
                    "portal_candidate/review/overview.png",
                    "ui_review_image",
                    "portal_candidate/review/overview.png",
                ),
                (
                    "passport_ui_review_thrb",
                    "portal_candidate/review/THRB.png",
                    "ui_review_image",
                    "portal_candidate/review/THRB.png",
                ),
                (
                    "passport_ui_review_hkdc1",
                    "portal_candidate/review/HKDC1.png",
                    "ui_review_image",
                    "portal_candidate/review/HKDC1.png",
                ),
                (
                    "passport_ui_review_glp1r",
                    "portal_candidate/review/GLP1R.png",
                    "ui_review_image",
                    "portal_candidate/review/GLP1R.png",
                ),
                (
                    "passport_ui_review_mtarc1",
                    "portal_candidate/review/MTARC1.png",
                    "ui_review_image",
                    "portal_candidate/review/MTARC1.png",
                ),
                (
                    "passport_ui_review_manifest",
                    "portal_candidate/review/review_manifest.tsv",
                    "ui_review_manifest",
                    "portal_candidate/review/review_manifest.tsv",
                ),
            ),
        ),
    }
    # The static subset above is retained only as readable role history.  The
    # executable Plan50 contract is the complete sealed payload inventory plus
    # its terminal chain, so a relocated snapshot can pass read-only PASS06
    # verification without resolving the live workstream root.
    contracts["PLAN50"] = (plan50, plan50_full_bundle_contract(plan50))
    return contracts


def publish(
    project_root: Path,
    workstream_id: str,
    signed_by: str,
    signed_at_utc: str,
) -> dict[str, object]:
    project = project_root.resolve()
    contracts = workstream_contract(project)
    root, artifacts = contracts[workstream_id]
    manifest = root / "plan60_terminal_artifacts.tsv"
    signature = root / "plan60_terminal_artifacts.signature.json"
    if (
        manifest.exists()
        or manifest.is_symlink()
        or signature.exists()
        or signature.is_symlink()
    ):
        raise ReleaseContractError(
            f"refusing to overwrite an existing {workstream_id} Plan60 handoff"
        )
    if not signed_by.strip() or not signed_at_utc.strip():
        raise ReleaseContractError("handoff signer and UTC timestamp are required")
    lane = WORKSTREAM_LANES[workstream_id]
    rows = []
    for artifact_id, relative, role, snapshot_relative in artifacts:
        source = root / relative
        if not source.is_file() or source.is_symlink():
            raise ReleaseContractError(
                f"{workstream_id} terminal artifact is missing or symlinked: {source}"
            )
        snapshot = f"inputs/{lane}/{snapshot_relative}"
        rows.append(
            {
                "artifact_id": artifact_id,
                "source_path": project_relative(project, source),
                "source_sha256": sha256_file(source),
                "source_bytes": source.stat().st_size,
                "snapshot_relpath": snapshot,
                "artifact_role": role,
                "downstream_read_path": snapshot,
            }
        )
    validate_handoff_rows_in_memory(project, workstream_id, rows)
    # Build and validate the pair in a private directory.  Final publication
    # uses no-overwrite hard links; any failure rolls back only links created by
    # this invocation, so a missing signature/post-write failure cannot leave
    # an apparently complete owner handoff.
    with tempfile.TemporaryDirectory(prefix=".plan60_handoff.", dir=root) as temporary:
        staging = Path(temporary)
        staged_manifest = staging / manifest.name
        staged_signature = staging / signature.name
        atomic_write_tsv(staged_manifest, rows, ARTIFACT_FIELDS)
        validated = validate_artifact_manifest(project, workstream_id, staged_manifest)
        if len(validated) != len(rows):
            raise ReleaseContractError("staged handoff row count drift")
        manifest_sha256 = sha256_file(staged_manifest)
        signature_payload = {
            "attestation_version": HANDOFF_SIGNATURE_VERSION,
            "candidate_id": CANDIDATE_ID,
            "workstream_id": workstream_id,
            "manifest_sha256": manifest_sha256,
            "signed_by": signed_by.strip(),
            "signed_at_utc": signed_at_utc.strip(),
            "canonical_promotion_authorized": False,
        }
        atomic_write_json(staged_signature, signature_payload)
        if (
            json.loads(staged_signature.read_text(encoding="utf-8"))
            != signature_payload
        ):
            raise ReleaseContractError("staged handoff signature does not reproduce")

        created: list[tuple[Path, int, int]] = []
        try:
            for staged, final in (
                (staged_manifest, manifest),
                (staged_signature, signature),
            ):
                if final.exists() or final.is_symlink():
                    raise ReleaseContractError(
                        f"handoff target appeared during publication: {final}"
                    )
                os.link(staged, final)
                final_stat = final.stat()
                created.append((final, final_stat.st_dev, final_stat.st_ino))
            final_validated = validate_artifact_manifest(
                project, workstream_id, manifest
            )
            if len(final_validated) != len(rows):
                raise ReleaseContractError("published handoff row count drift")
            if (
                sha256_file(manifest) != manifest_sha256
                or json.loads(signature.read_text(encoding="utf-8"))
                != signature_payload
            ):
                raise ReleaseContractError("published handoff pair drift")
        except Exception:
            for created_path, device, inode in reversed(created):
                if created_path.is_file() and not created_path.is_symlink():
                    current = created_path.stat()
                    if (current.st_dev, current.st_ino) == (device, inode):
                        created_path.unlink()
            raise
    return {
        "workstream_id": workstream_id,
        "manifest": str(manifest),
        "manifest_sha256": manifest_sha256,
        "n_artifacts": len(rows),
        "canonical_promotion_authorized": False,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--workstream", choices=tuple(WORKSTREAM_LANES), required=True)
    parser.add_argument("--signed-by", required=True)
    parser.add_argument("--signed-at-utc", required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    try:
        result = publish(
            args.project_root,
            args.workstream,
            args.signed_by,
            args.signed_at_utc,
        )
    except ReleaseContractError as error:
        raise SystemExit(f"HANDOFF_BLOCKED: {error}") from error
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
