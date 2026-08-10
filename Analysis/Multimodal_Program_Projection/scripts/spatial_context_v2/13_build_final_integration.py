#!/usr/bin/env python3
"""Build SP-INT-06/07 Figure 4 source tables in an isolated candidate root.

``--mode fixture`` is outcome-free infrastructure testing.  ``--mode real``
hard-fails until both native Visium READY seals exist and then performs a
deterministic union of already-validated assay-native evidence.
"""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping

from contract_lib import (
    ADAPTER_REGISTRY_COLUMNS,
    ContractError,
    RELEASE_ID,
    SCHEMAS,
    parse_bool,
    read_tsv,
    sha256_file,
    validate_candidate_contract,
    write_tsv,
)
from final_integration_lib import (
    FIGURE4_MATRIX_COLUMNS,
    HOTFIX_EXPECTED_SHA256,
    INPUT_SEAL_COLUMNS,
    INTEGRATED_EFFECT_COLUMNS,
    SOURCE_TABLE_MANIFEST_COLUMNS,
    VERDICT_COLUMNS,
    adapter_row,
    count_data_rows,
    expected_programs,
    native_spatial_rows,
    project_relative,
    require_complete_groups,
    require_no_cross_assay_construct,
    seal_rows,
    verify_hotfix,
    write_effect_outputs,
)


STATUS_COLUMNS = (
    "release_id",
    "status",
    "mode",
    "n_programs",
    "n_dataset_assay_groups",
    "n_integrated_program_rows",
    "n_figure4_source_tables",
    "input_seal_manifest_sha256",
    "integrated_effects_sha256",
    "figure4_matrix_sha256",
    "figure4_verdict_sha256",
    "source_table_manifest_sha256",
    "hotfix_manifest_sha256",
    "cross_assay_score_constructed",
    "pdf_rendered",
    "canonical_write_authorized",
    "plan13_complete",
    "built_utc",
)


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def one_row(path: Path, required: tuple[str, ...]) -> dict[str, str]:
    _, rows = read_tsv(path, required)
    if len(rows) != 1:
        raise ContractError(f"expected one row in {path}, found {len(rows)}")
    return rows[0]


def ready_validation_report(root: Path, ready: Mapping[str, str]) -> Path:
    """Resolve the READY-linked validation report without assuming its filename."""
    relative = ready.get("validation_report_path", "validation_report.tsv").strip()
    if not relative or Path(relative).is_absolute() or ".." in Path(relative).parts:
        raise ContractError(f"invalid READY validation_report_path: {relative!r}")
    report = root / relative
    if not report.is_file() or ready["validation_report_sha256"] != sha256_file(report):
        raise ContractError(f"native spatial validation link drift: {report}")
    return report


def verify_output_manifest(root: Path) -> None:
    _, rows = read_tsv(root / "output_manifest.tsv", ("relative_path", "bytes", "sha256"))
    if not rows:
        raise ContractError(f"empty native output manifest: {root}")
    for row in rows:
        path = root / row["relative_path"]
        if not path.is_file() or path.stat().st_size != int(row["bytes"]) or sha256_file(path) != row["sha256"]:
            raise ContractError(f"native output-manifest drift: {path}")


def verify_code_freeze(project_root: Path, candidate_root: Path) -> Path:
    path = candidate_root / "final_integration_code_freeze.tsv"
    _, rows = read_tsv(path, ("relative_path", "bytes", "sha256", "role"))
    if not rows:
        raise ContractError("final integration code freeze is empty")
    for row in rows:
        source = project_root / row["relative_path"]
        if not source.is_file() or source.stat().st_size != int(row["bytes"]) or sha256_file(source) != row["sha256"]:
            raise ContractError(f"final integration code-freeze drift: {source}")
    return path


def verify_real_seals(
    project_root: Path,
    candidate_root: Path,
    hotspot_root: Path,
    yak_native_root: Path,
    yak_adapter_root: Path,
) -> dict[str, Path]:
    real_ready = candidate_root / "REAL_ADAPTERS_READY"
    real = one_row(
        real_ready,
        (
            "status", "real_adapter_registry_sha256", "bundle_manifest_sha256",
            "validation_report_sha256", "n_real_adapters", "plan13_complete",
        ),
    )
    if real["status"] != "ready_partial_real_adapter_seam_not_plan13_complete" or real["n_real_adapters"] != "3":
        raise ContractError("real-adapter READY status/count drift")
    if parse_bool(real["plan13_complete"], "real.plan13"):
        raise ContractError("partial real-adapter READY overstates Plan 13 completion")
    for field, path in (
        ("real_adapter_registry_sha256", candidate_root / "real_adapter_registry.tsv"),
        ("bundle_manifest_sha256", candidate_root / "real_adapter_bundle_manifest.tsv"),
        ("validation_report_sha256", candidate_root / "real_adapter_validation.tsv"),
    ):
        if real[field] != sha256_file(path):
            raise ContractError(f"real-adapter READY link drift: {field}")

    protein_root = candidate_root / "protein_atac"
    protein_ready = protein_root / "READY"
    protein = one_row(
        protein_ready,
        (
            "status", "adapter_registry_sha256", "bundle_manifest_sha256",
            "validation_report_sha256", "n_adapters", "plan13_complete",
        ),
    )
    if protein["status"] != "ready_sp_int_05_candidate_only_not_plan13_complete" or protein["n_adapters"] != "3":
        raise ContractError("protein/ATAC READY status/count drift")
    if parse_bool(protein["plan13_complete"], "protein.plan13"):
        raise ContractError("SP-INT-05 READY overstates Plan 13 completion")
    for field, path in (
        ("adapter_registry_sha256", protein_root / "adapter_registry.tsv"),
        ("bundle_manifest_sha256", protein_root / "bundle_manifest.tsv"),
        ("validation_report_sha256", protein_root / "validation.tsv"),
    ):
        if protein[field] != sha256_file(path):
            raise ContractError(f"protein/ATAC READY link drift: {field}")

    plan12_gate = project_root / "Analysis/Spatial/candidates" / RELEASE_ID / "gse287826/gate_status.tsv"
    gate = one_row(plan12_gate, ("status", "inference_authorized", "outcomes_read", "reason"))
    if gate["status"] != "skipped_no_donor_key":
        raise ContractError("GSE287826 terminal gate no longer records skipped_no_donor_key")
    if parse_bool(gate["inference_authorized"], "gse.inference") or parse_bool(gate["outcomes_read"], "gse.outcomes"):
        raise ContractError("GSE287826 skip gate improperly authorizes or reads inference")

    yak_root = yak_native_root
    yak_ready = yak_root / "READY"
    native_yak = one_row(
        yak_ready,
        ("status", "candidate_release_manifest_sha256", "n_frozen_family", "n_robust_programs", "canonical_promotion_authorized"),
    )
    if native_yak["status"] != "ready_for_plan13_candidate_integration" or native_yak["n_frozen_family"] != "2":
        raise ContractError("Plan 11 native READY status/family drift")
    if parse_bool(native_yak["canonical_promotion_authorized"], "yak.promotion"):
        raise ContractError("Plan 11 candidate improperly authorizes canonical promotion")
    if native_yak["candidate_release_manifest_sha256"] != sha256_file(yak_root / "candidate_release_manifest.tsv"):
        raise ContractError("Plan 11 native READY candidate-release link drift")
    verify_hotfix(yak_root / "validator_hotfix_manifest.tsv")
    yak_adapter = yak_adapter_root
    yak_adapter_ready = yak_adapter / "PLAN13_ADAPTER_READY"
    yak = one_row(
        yak_adapter_ready,
        (
            "status", "native_ready_sha256", "adapter_registry_sha256",
            "validation_report_sha256", "canonical_promotion_authorized",
        ),
    )
    if yak["status"] != "validated_plan13_adapter_candidate":
        raise ContractError("Plan 11 -> Plan 13 adapter READY status drift")
    if parse_bool(yak["canonical_promotion_authorized"], "yak.adapter.promotion"):
        raise ContractError("Plan 11 adapter improperly authorizes canonical promotion")
    if yak["native_ready_sha256"] != sha256_file(yak_ready):
        raise ContractError("Plan 11 adapter does not pin native READY")
    if yak["adapter_registry_sha256"] != sha256_file(yak_adapter / "adapter_registry.tsv"):
        raise ContractError("Plan 11 adapter registry link drift")
    if yak["validation_report_sha256"] != sha256_file(yak_adapter / "validation_report.tsv"):
        raise ContractError("Plan 11 adapter validation link drift")

    native_parent = candidate_root / "native_spatial"
    v1 = native_parent / "v1_regression"
    v2 = native_parent / "v2_candidate"
    v1_ready = one_row(v1 / "READY", ("registry_version", "status", "validation_report_sha256"))
    v2_ready = one_row(v2 / "READY", ("registry_version", "status", "validation_report_sha256"))
    if v1_ready["registry_version"] != "v1" or v1_ready["status"] != "pass_v1_regression":
        raise ContractError("native spatial v1 READY is absent or invalid")
    if v2_ready["registry_version"] != "v2" or v2_ready["status"] != "pass_v2_candidate":
        raise ContractError("native spatial v2 READY is absent or invalid")
    for root, ready in ((v1, v1_ready), (v2, v2_ready)):
        ready_validation_report(root, ready)
        verify_output_manifest(root)
    execution = {row["parameter"]: row["value"] for row in read_tsv(v2 / "execution_manifest.tsv")[1]}
    if execution.get("v1_regression_ready_sha256") != sha256_file(v1 / "READY"):
        raise ContractError("native spatial v2 does not pin the passing v1 READY")
    if execution.get("canonical_write_allowed") != "FALSE" or execution.get("program_count") != "2":
        raise ContractError("native spatial v2 execution contract drift")

    validate_candidate_contract(candidate_root, candidate_root / "real_adapter_registry.tsv", hotspot_root, project_root)
    validate_candidate_contract(candidate_root, protein_root / "adapter_registry.tsv", hotspot_root, project_root)
    validate_candidate_contract(yak_adapter, yak_adapter / "adapter_registry.tsv", hotspot_root, project_root)
    return {
        "real_ready": real_ready,
        "protein_ready": protein_ready,
        "plan12_gate": plan12_gate,
        "yak_ready": yak_ready,
        "yak_adapter_ready": yak_adapter_ready,
        "yak_hotfix": yak_root / "validator_hotfix_manifest.tsv",
        "native_v1_ready": v1 / "READY",
        "native_v2_ready": v2 / "READY",
    }


def semantics_for(row: Mapping[str, str]) -> str:
    state = row["evidence_state"]
    if state in {"untestable", "not_applicable", "skipped"}:
        return "not_applicable_terminal_state"
    if row["interval_type"] == "donor_effect_range":
        return "descriptive_donor_effect_range_not_confidence_interval"
    if row["dataset"] in {"GSE281367", "GSE244832"}:
        return "Welch_mean_difference_SE_with_exact_label_permutation_primary_pvalue"
    return "source_native_adapter_uncertainty"


ROLE_MAP = {
    "gse287826_skip": ("source_gate_boundary", "independent disease replication not estimable without donor key", "retain_transparent_skip_rows"),
    "govaere2026_geomx_dropped": ("transparent_assay_exclusion", "excluded local GeoMx reanalysis", "retain_exclusion_boundary"),
    "govaere2026_cosmx_il32_context": ("cell_proximity_context", "source-native IL32 proximity only; whole-program inference not applicable", "retain_source_native_context_only"),
    "pxd051911_diams_program": ("protein_observability", "frozen-program DIA observability boundary", "retain_program_untestable_rows_and_fixed25_descriptive"),
    "gse281367_snatac_dynamic": ("chromatin_context", "independent donor-level promoter accessibility", "retain_complete_independent_ATAC_panel"),
    "gse244832_snatac_dynamic": ("chromatin_context", "source-dependent donor-level promoter accessibility", "retain_complete_source_dependent_ATAC_panel"),
    "yakubovsky2026_binary_lipid": (
        "local_lipid_context",
        "tissue-spot projection of a hepatocyte-derived signature without hepatocyte filtering; source-defined lipid-zone association after zonation adjustment; cell composition remains a rival explanation",
        "retain_complete_indeterminate_lipid_panel_with_separate_descriptive_and_inferential_directions",
    ),
}


def import_registered_adapters(
    project_root: Path,
    candidate_root: Path,
    registry_path: Path,
    registry_base: Path,
    ready_path: Path,
) -> tuple[list[dict[str, object]], list[tuple[Path, str, str]]]:
    _, adapters = read_tsv(registry_path, ADAPTER_REGISTRY_COLUMNS)
    ready_rel = project_relative(project_root, ready_path)
    ready_sha = sha256_file(ready_path)
    rows = []
    sources: list[tuple[Path, str, str]] = [(ready_path, "validated_adapter_ready", "sealed")]
    for adapter in adapters:
        adapter_id = adapter["adapter_id"]
        if adapter_id not in ROLE_MAP:
            raise ContractError(f"final integration lacks an explicit role for adapter {adapter_id}")
        role, claim, interpretation = ROLE_MAP[adapter_id]
        effect_path = registry_base / adapter["adapter_root"] / "program_effects.tsv"
        _, effects = read_tsv(effect_path, SCHEMAS["program_effects.tsv"])
        effect_rel = project_relative(project_root, effect_path)
        effect_sha = sha256_file(effect_path)
        sources.append((effect_path, f"{adapter_id}_program_effects", "validated"))
        for effect in effects:
            rows.append(
                adapter_row(
                    effect,
                    adapter_id=adapter_id,
                    evidence_role=role,
                    claim_scope=claim,
                    uncertainty_semantics=semantics_for(effect),
                    ready_path=ready_rel,
                    ready_sha256=ready_sha,
                    effect_path=effect_rel,
                    effect_sha256=effect_sha,
                    interpretation=interpretation,
                )
            )
    sources.append((registry_path, "validated_adapter_registry", "sealed"))
    return rows, sources


def copy_source_native(project_root: Path, candidate_root: Path, staging: Path) -> list[dict[str, object]]:
    source_root = staging / "source_native"
    source_root.mkdir(parents=True)
    copies = (
        (
            candidate_root / "real_adapters/govaere2026_cosmx_il32_context/source_native_context.tsv",
            source_root / "govaere_cosmx_il32_context.tsv",
            "govaere_cosmx_source_native",
            "cell_proximity_context",
            "array-level source-native IL32 proximity; no whole-program refit",
        ),
        (
            candidate_root / "protein_atac/adapters/pxd051911_diams_program/fixed_25_protein_context.tsv",
            source_root / "fixed_25_protein_context.tsv",
            "fixed_25_protein_context",
            "protein_context",
            "selection-conditioned descriptive same-cohort protein display",
        ),
        (
            candidate_root / "protein_atac/source_native/atac_static_context/static_atac_context.tsv",
            source_root / "static_atac_context.tsv",
            "static_atac_context",
            "static_chromatin_context",
            "descriptive fine-mapped variant/open-promoter program context",
        ),
        (
            candidate_root / "protein_atac/source_native/atac_static_context/static_atac_hits.tsv",
            source_root / "static_atac_hits.tsv",
            "static_atac_hits",
            "static_chromatin_context",
            "source-native gene-level fine-mapped/open-promoter links",
        ),
        (
            candidate_root / "protein_atac/source_native/atac_static_context/static_study_coverage.tsv",
            source_root / "static_atac_study_coverage.tsv",
            "static_atac_study_coverage",
            "static_chromatin_provenance",
            "prespecified study-coverage context; not a program effect",
        ),
    )
    rows = []
    for source, target, table_id, role, interpretation in copies:
        if not source.is_file():
            raise ContractError(f"missing source-native Figure 4 table: {source}")
        shutil.copyfile(source, target)
        if sha256_file(source) != sha256_file(target):
            raise ContractError(f"source-native copy drift: {source}")
        rows.append(
            {
                "release_id": RELEASE_ID,
                "table_id": table_id,
                "relative_path": target.relative_to(staging).as_posix(),
                "data_rows": count_data_rows(target),
                "bytes": target.stat().st_size,
                "sha256": sha256_file(target),
                "figure4_role": role,
                "interpretation": interpretation,
            }
        )
    return rows


def real_input_paths(
    project_root: Path,
    candidate_root: Path,
    seals: Mapping[str, Path],
    code_freeze: Path,
    adapter_sources: list[tuple[Path, str, str]],
    yak_native_root: Path,
    yak_adapter_root: Path,
) -> list[tuple[Path, str, str]]:
    paths: list[tuple[Path, str, str]] = [(code_freeze, "outcome_blind_final_assembly_code_freeze", "sealed")]
    paths.extend((path, role, "sealed") for role, path in seals.items())
    paths.extend(adapter_sources)
    for relative, role in (
        ("real_adapter_bundle_manifest.tsv", "real_adapter_bundle_manifest"),
        ("real_adapter_validation.tsv", "real_adapter_validation"),
        ("protein_atac/bundle_manifest.tsv", "protein_atac_bundle_manifest"),
        ("protein_atac/validation.tsv", "protein_atac_validation"),
        ("native_spatial/v1_regression/output_manifest.tsv", "native_spatial_v1_output_manifest"),
        ("native_spatial/v2_candidate/output_manifest.tsv", "native_spatial_v2_output_manifest"),
        ("native_spatial/v2_candidate/execution_manifest.tsv", "native_spatial_v2_execution_manifest"),
        ("native_spatial/v2_candidate/source_manifest.tsv", "native_spatial_v2_source_manifest"),
        ("native_spatial/v2_candidate/tested_universe.tsv", "native_spatial_v2_tested_universe"),
        ("native_spatial/v2_candidate/spatial_program_results.tsv", "native_spatial_v2_program_effects"),
        ("native_spatial/v2_candidate/native_design_audit.tsv", "native_spatial_v2_design"),
        ("native_spatial/v2_candidate/spatial_testability_audit.tsv", "native_spatial_v2_testability"),
        ("native_spatial/v2_candidate/spatial_sensitivity_audit.tsv", "native_spatial_v2_sensitivity"),
        ("native_spatial/v2_candidate/spatial_null_summary.tsv", "native_spatial_v2_null_summary"),
    ):
        paths.append((candidate_root / relative, role, "validated"))
    for relative, role in (
        ("native_spatial/v1_regression", "native_spatial_v1_validation_report"),
        ("native_spatial/v2_candidate", "native_spatial_v2_validation_report"),
    ):
        native_root = candidate_root / relative
        native_ready = one_row(
            native_root / "READY", ("validation_report_sha256",)
        )
        paths.append((ready_validation_report(native_root, native_ready), role, "validated"))
    paths.append(
        (
            yak_native_root / "candidate_release_manifest.tsv",
            "plan11_candidate_release_manifest",
            "sealed",
        )
    )
    for relative, role in (
        ("source_manifest.tsv", "plan11_adapter_source_manifest"),
        ("execution_manifest.tsv", "plan11_adapter_execution_manifest"),
        ("validation_report.tsv", "plan11_adapter_validation"),
    ):
        paths.append((yak_adapter_root / relative, role, "sealed"))
    gse_root = project_root / "Analysis/Spatial/candidates" / RELEASE_ID / "gse287826"
    paths.append((gse_root / "execution_manifest.tsv", "plan12_skip_execution_manifest", "sealed"))
    return paths


def build_real(
    project_root: Path,
    candidate_root: Path,
    staging: Path,
    yak_native_root: Path,
    yak_adapter_root: Path,
) -> None:
    hotspot_root = candidate_root.parent / "hotspot"
    selected = expected_programs(hotspot_root)
    code_freeze = verify_code_freeze(project_root, candidate_root)
    seals = verify_real_seals(
        project_root, candidate_root, hotspot_root, yak_native_root, yak_adapter_root
    )
    rows = []
    adapter_sources: list[tuple[Path, str, str]] = []
    for registry_path, registry_base, ready in (
        (candidate_root / "real_adapter_registry.tsv", candidate_root, seals["real_ready"]),
        (candidate_root / "protein_atac/adapter_registry.tsv", candidate_root, seals["protein_ready"]),
    ):
        imported, sources = import_registered_adapters(
            project_root, candidate_root, registry_path, registry_base, ready
        )
        rows.extend(imported)
        adapter_sources.extend(sources)
    imported, sources = import_registered_adapters(
        project_root,
        candidate_root,
        yak_adapter_root / "adapter_registry.tsv",
        yak_adapter_root,
        seals["yak_adapter_ready"],
    )
    rows.extend(imported)
    adapter_sources.extend(sources)
    rows.extend(native_spatial_rows(project_root, candidate_root / "native_spatial/v2_candidate", hotspot_root))
    require_complete_groups(rows, selected)
    if len(rows) != 18:
        raise ContractError(f"final integration expected 18 program rows across nine groups, found {len(rows)}")
    require_no_cross_assay_construct(rows, INTEGRATED_EFFECT_COLUMNS)

    input_rows = seal_rows(
        project_root,
        real_input_paths(
            project_root,
            candidate_root,
            seals,
            code_freeze,
            adapter_sources,
            yak_native_root,
            yak_adapter_root,
        ),
    )
    write_tsv(staging / "input_seal_manifest.tsv", INPUT_SEAL_COLUMNS, input_rows)
    write_effect_outputs(staging, rows)
    source_tables = copy_source_native(project_root, candidate_root, staging)
    for table_id, filename, role, interpretation in (
        ("integrated_program_effects", "integrated_program_effects.tsv", "complete_program_evidence", "assay-native long contract; no cross-assay score"),
        ("figure4_program_matrix", "figure4_program_matrix.tsv", "program_panel_source", "all preselected rows including null and terminal states"),
        ("figure4_dataset_verdict", "figure4_dataset_verdict.tsv", "dataset_adjudication", "panel inclusion boundary without positivity requirement"),
    ):
        path = staging / filename
        source_tables.append(
            {
                "release_id": RELEASE_ID,
                "table_id": table_id,
                "relative_path": filename,
                "data_rows": count_data_rows(path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "figure4_role": role,
                "interpretation": interpretation,
            }
        )
    write_tsv(staging / "figure4_source_table_manifest.tsv", SOURCE_TABLE_MANIFEST_COLUMNS, source_tables)
    write_tsv(
        staging / "assembly_status.tsv",
        STATUS_COLUMNS,
        [
            {
                "release_id": RELEASE_ID,
                "status": "assembled_pending_independent_validation",
                "mode": "real",
                "n_programs": 2,
                "n_dataset_assay_groups": 9,
                "n_integrated_program_rows": len(rows),
                "n_figure4_source_tables": len(source_tables),
                "input_seal_manifest_sha256": sha256_file(staging / "input_seal_manifest.tsv"),
                "integrated_effects_sha256": sha256_file(staging / "integrated_program_effects.tsv"),
                "figure4_matrix_sha256": sha256_file(staging / "figure4_program_matrix.tsv"),
                "figure4_verdict_sha256": sha256_file(staging / "figure4_dataset_verdict.tsv"),
                "source_table_manifest_sha256": sha256_file(staging / "figure4_source_table_manifest.tsv"),
                "hotfix_manifest_sha256": HOTFIX_EXPECTED_SHA256,
                "cross_assay_score_constructed": False,
                "pdf_rendered": False,
                "canonical_write_authorized": False,
                "plan13_complete": False,
                "built_utc": utc_now(),
            }
        ],
    )


def synthetic_effect(reg: Mapping[str, str], dataset: str, state: str, index: int) -> dict[str, object]:
    tested = state in {"robust", "indeterminate", "tested_negative"}
    estimate = 0.2 + index / 10 if tested else ""
    base = {
        "release_id": RELEASE_ID,
        "registry_sha256": "fixture_registry_sha256",
        "membership_sha256": reg["membership_sha256"],
        "program_uid": reg["program_uid"],
        "legacy_program_id": f"{reg['cell_type']}:{reg['module']}",
        "program_label": reg["module_name"],
        "cell_type": reg["cell_type"],
        "dataset": dataset,
        "assay": "synthetic_contract_fixture",
        "source_publication": "synthetic_no_real_outcome",
        "source_dependence": "independent" if dataset in {"FixtureSpatial", "FixtureProtein"} else "source_dependent",
        "analysis_set_id": "synthetic_null_fixture",
        "biological_unit": "donor" if dataset != "FixtureSkip" else "unknown_public_biological_unit",
        "biological_unit_resolution": "resolved" if dataset != "FixtureSkip" else "unresolved",
        "n_biological": 6 if dataset != "FixtureSkip" else "",
        "technical_unit": "synthetic_unit",
        "n_technical": 6 if dataset != "FixtureSkip" else 0,
        "contrast_or_exposure": "synthetic_exposure",
        "effect_unit": "synthetic_assay_native_unit",
        "estimate": estimate,
        "std_error": 0.1 if tested else "",
        "matched_null_sd": "",
        "interval_low": "",
        "interval_high": "",
        "interval_type": "none" if tested else "not_applicable",
        "pvalue": 0.02 + index / 100 if tested else "",
        "padj": 0.04 + index / 100 if tested else "",
        "pvalue_method": "synthetic_exact_null" if tested else "not_applicable",
        "multiplicity_family": "synthetic_complete_two_program_family",
        "n_genes_measured": 12 if tested else 0,
        "retained_l1_weight": 0.5 if tested else 0.0,
        "testable": tested,
        "testability_reason": "synthetic_fixture",
        "direction_expected": "positive",
        "direction_observed": "positive" if tested else "",
        "descriptive_effect_direction": "positive" if tested else "",
        "inferential_test_direction": "positive" if tested else "",
        "direction_agreement": tested,
        "heterogeneity_statistic": "",
        "heterogeneity_df": "",
        "heterogeneity_pvalue": "",
        "sensitivity_sign_agree": tested,
        "robustness_pass": state == "robust",
        "negative_call_rule_id": (
            "fixture_adequate_negative_margin_rule_v1" if state == "tested_negative" else ""
        ),
        "cross_assay_comparable": False,
        "evidence_state": state,
        "producer": "synthetic_fixture",
        "producer_sha256": "synthetic_fixture",
        "source_manifest_sha256": "synthetic_fixture",
        "source_adapter_id": f"{dataset.lower()}_fixture",
        "evidence_role": "synthetic_contract_role",
        "claim_scope": "synthetic_only",
        "uncertainty_semantics": "synthetic_only",
        "import_method": "synthetic_fixture",
        "source_ready_path": "synthetic_fixture",
        "source_ready_sha256": "synthetic_fixture",
        "source_effect_path": "synthetic_fixture",
        "source_effect_sha256": "synthetic_fixture",
        "source_effect_row_sha256": "synthetic_fixture",
        "figure4_include": True,
        "figure4_interpretation": "synthetic_only",
    }
    return base


def build_fixture(project_root: Path, candidate_root: Path, staging: Path) -> None:
    hotspot_root = candidate_root.parent / "hotspot"
    selected = expected_programs(hotspot_root)
    states = {
        "FixtureSpatial": ("robust", "indeterminate"),
        "FixtureLipid": ("indeterminate", "tested_negative"),
        "FixtureProtein": ("untestable", "untestable"),
        "FixtureSkip": ("skipped", "skipped"),
    }
    rows = []
    for dataset, dataset_states in states.items():
        for index, uid in enumerate(sorted(selected)):
            rows.append(synthetic_effect(selected[uid], dataset, dataset_states[index], index))
    require_complete_groups(rows, selected)
    require_no_cross_assay_construct(rows, INTEGRATED_EFFECT_COLUMNS)
    write_tsv(
        staging / "input_seal_manifest.tsv",
        INPUT_SEAL_COLUMNS,
        seal_rows(
            project_root,
            [
                (hotspot_root / "READY", "fixture_registry_ready_only_no_external_outcome", "sealed"),
                (hotspot_root / "program_registry_v2.tsv", "fixture_registry_only_no_external_outcome", "sealed"),
            ],
        ),
    )
    write_effect_outputs(staging, rows)
    source_tables = []
    for table_id, filename in (
        ("integrated_program_effects", "integrated_program_effects.tsv"),
        ("figure4_program_matrix", "figure4_program_matrix.tsv"),
        ("figure4_dataset_verdict", "figure4_dataset_verdict.tsv"),
    ):
        path = staging / filename
        source_tables.append(
            {
                "release_id": RELEASE_ID,
                "table_id": table_id,
                "relative_path": filename,
                "data_rows": count_data_rows(path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "figure4_role": "synthetic_contract_fixture",
                "interpretation": "synthetic/null infrastructure test; no real outcome",
            }
        )
    write_tsv(staging / "figure4_source_table_manifest.tsv", SOURCE_TABLE_MANIFEST_COLUMNS, source_tables)
    write_tsv(
        staging / "assembly_status.tsv",
        STATUS_COLUMNS,
        [
            {
                "release_id": RELEASE_ID,
                "status": "synthetic_fixture_assembled_pending_validation",
                "mode": "fixture",
                "n_programs": 2,
                "n_dataset_assay_groups": 4,
                "n_integrated_program_rows": 8,
                "n_figure4_source_tables": 3,
                "input_seal_manifest_sha256": sha256_file(staging / "input_seal_manifest.tsv"),
                "integrated_effects_sha256": sha256_file(staging / "integrated_program_effects.tsv"),
                "figure4_matrix_sha256": sha256_file(staging / "figure4_program_matrix.tsv"),
                "figure4_verdict_sha256": sha256_file(staging / "figure4_dataset_verdict.tsv"),
                "source_table_manifest_sha256": sha256_file(staging / "figure4_source_table_manifest.tsv"),
                "hotfix_manifest_sha256": "not_applicable_fixture",
                "cross_assay_score_constructed": False,
                "pdf_rendered": False,
                "canonical_write_authorized": False,
                "plan13_complete": False,
                "built_utc": utc_now(),
            }
        ],
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=("fixture", "real"), required=True)
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--candidate-root", type=Path, default=None)
    parser.add_argument("--yak-native-root", type=Path, default=None)
    parser.add_argument("--yak-adapter-root", type=Path, default=None)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    candidate_root = args.candidate_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "spatial_context"
    )
    yak_native_root = args.yak_native_root or (
        project_root / "Analysis/Spatial/candidates" / RELEASE_ID / "yakubovsky2026"
    )
    yak_adapter_root = args.yak_adapter_root or (yak_native_root / "plan13_adapter")
    target = candidate_root / ("final_integration_fixture" if args.mode == "fixture" else "final_integration")
    staging = candidate_root / f".{target.name}.incomplete.{os.environ.get('SLURM_JOB_ID', os.getpid())}"
    try:
        if target.exists():
            raise ContractError(f"refusing to overwrite final integration target: {target}")
        if staging.exists():
            raise ContractError(f"staging path already exists: {staging}")
        staging.mkdir(parents=True)
        if args.mode == "fixture":
            build_fixture(project_root, candidate_root, staging)
        else:
            build_real(
                project_root,
                candidate_root,
                staging,
                yak_native_root.resolve(),
                yak_adapter_root.resolve(),
            )
        os.replace(staging, target)
        print(f"PASS: assembled {args.mode} final integration pending independent validation")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
