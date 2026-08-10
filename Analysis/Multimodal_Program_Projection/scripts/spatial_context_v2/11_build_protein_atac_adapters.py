#!/usr/bin/env python3
"""Assemble candidate-only protein and ATAC adapters for SP-INT-05.

No protein outcome is refit: both frozen v2 programs fail the prespecified
PXD051911 observability threshold, so their program rows are terminal
``untestable``.  The protected 25-protein display is copied byte-for-byte and
remains explicitly selection-conditioned/descriptive.

The dynamic snATAC rows consume the bounded v2 producer output from
``10_run_atac_v2_projection.R``.  Static fine-mapped/open-promoter links are
re-keyed source-native rows, not inferential program tests.
"""

from __future__ import annotations

import argparse
import csv
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from contract_lib import (
    ADAPTER_REGISTRY_COLUMNS,
    ContractError,
    RELEASE_ID,
    SCHEMAS,
    canonical_id_hash,
    parse_bool,
    parse_float,
    parse_int,
    program_universe,
    read_tsv,
    sha256_file,
    validate_hotspot_seal,
    write_tsv,
)


STATIC_CONTEXT_COLUMNS = (
    "release_id",
    "registry_sha256",
    "program_uid",
    "membership_sha256",
    "legacy_program_id",
    "program_label",
    "cell_type",
    "dataset",
    "assay",
    "source_dependence",
    "trait_scope",
    "n_source_loci",
    "n_study_loci",
    "n_variants",
    "n_genes",
    "n_studies",
    "n_ancestries",
    "n_registry_studies",
    "n_studies_with_credible_sets",
    "effect_unit",
    "interpretation",
)

STATIC_HIT_COLUMNS = (
    "release_id",
    "registry_sha256",
    "program_uid",
    "membership_sha256",
    "legacy_program_id",
    "cell_type",
    "gene_symbol",
    "study",
    "ancestry",
    "trait",
    "trait_scope",
    "source_locus",
    "chromosome",
    "locus",
    "variant_id",
    "recommended_pip",
    "chr_hg38",
    "pos_hg38",
    "source_dependence",
    "interpretation",
)


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _project_relative(project_root: Path, path: Path) -> str:
    try:
        return path.resolve().relative_to(project_root.resolve()).as_posix()
    except ValueError as exc:
        raise ContractError(f"path outside project root: {path}") from exc


def _read_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise ContractError(f"missing CSV: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            raise ContractError(f"missing CSV header: {path}")
        return [dict(row) for row in reader]


def _source_manifest(
    project_root: Path,
    adapter_root: Path,
    inputs: Sequence[tuple[Path, str, bool, str]],
) -> str:
    rows = []
    seen: set[Path] = set()
    for path, role, public_access, retrieved in inputs:
        resolved = path.resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        if not path.is_file():
            raise ContractError(f"missing source-manifest input: {path}")
        rows.append(
            {
                "path_scope": "project_relative",
                "relative_path": _project_relative(project_root, path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "source_role": role,
                "public_access": public_access,
                "retrieved_utc": retrieved,
            }
        )
    write_tsv(adapter_root / "source_manifest.tsv", SCHEMAS["source_manifest.tsv"], rows)
    return sha256_file(adapter_root / "source_manifest.tsv")


def _execution_manifest(adapter_root: Path, filenames: Iterable[str], created: str) -> str:
    rows = []
    for filename in filenames:
        path = adapter_root / filename
        if not path.is_file():
            raise ContractError(f"cannot seal missing adapter artifact: {path}")
        rows.append(
            {
                "role": "sp_int_05_candidate_artifact",
                "relative_path": filename,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "created_utc": created,
            }
        )
    write_tsv(adapter_root / "execution_manifest.tsv", SCHEMAS["execution_manifest.tsv"], rows)
    return sha256_file(adapter_root / "execution_manifest.tsv")


def _write_gate(
    root: Path,
    spec: Mapping[str, object],
    seal,
    source_hash: str,
    execution_hash: str,
    status: str,
    reason: str,
    outcomes_tested: bool,
) -> None:
    write_tsv(
        root / "gate_status.tsv",
        SCHEMAS["gate_status.tsv"],
        [
            {
                "release_id": RELEASE_ID,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "gate_status": status,
                "gate_reason": reason,
                "outcomes_tested": outcomes_tested,
                "registry_sha256": seal.registry_sha256,
                "source_manifest_sha256": source_hash,
                "execution_manifest_sha256": execution_hash,
            }
        ],
    )


def _mapped_membership(hotspot_root: Path, selected: Mapping[str, Mapping[str, str]]) -> dict[str, list[dict[str, str]]]:
    _, rows = read_tsv(
        hotspot_root / "program_membership_v2.tsv",
        (
            "program_uid",
            "membership_sha256",
            "mapped_symbol",
            "mapped_symbol_status",
            "original_l1_weight",
        ),
    )
    result = {uid: [] for uid in selected}
    for row in rows:
        uid = row["program_uid"]
        if uid not in result:
            continue
        symbol = row["mapped_symbol"].strip()
        if not symbol or symbol.upper() == "NA" or row["mapped_symbol_status"] == "unmapped_or_ambiguous":
            continue
        result[uid].append(row)
    if any(not rows for rows in result.values()):
        raise ContractError("selected v2 program lacks mapped membership")
    return result


def _legacy_id(row: Mapping[str, str]) -> str:
    """Return the v2 registry's canonical legacy display identifier."""
    return f"{row['cell_type']}:{row['module']}"


def _v1_program_id(row: Mapping[str, str]) -> str:
    """Return the historical v1 program_id used by accepted source tables."""
    return f"{row['cell_type']}::{row['module']}"


def _protein_spec(seal) -> dict[str, object]:
    return {
        "adapter_id": "pxd051911_diams_program",
        "adapter_root": "protein_atac/adapters/pxd051911_diams_program",
        "dataset": "PXD051911",
        "assay": "liver_DIA_MS",
        "source_publication": "PXD051911_ProteomeXchange_public_record",
        "source_dependence": "independent",
        "biological_unit": "patient",
        "biological_unit_resolution": "resolved",
        "technical_unit": "liver_DIA_sample",
        "analysis_set_id": "MASLD_vs_Control_adjusted_program_projection",
        "contrast_or_exposure": "MASLD_minus_Control",
        "effect_unit": "adjusted_patient_program_score_difference",
        "cross_assay_comparable": False,
        "program_universe": "external_test_eligible",
        "min_genes_testable": 8,
        "min_retained_l1_weight": 0.2,
        "gate_expectation": "untestable",
        "registry_sha256": seal.registry_sha256,
        "ready_sha256": seal.ready_sha256,
        "status": "real_v2_program_unobservable_fixed25_descriptive_only",
    }


def _atac_spec(seal, cohort: str) -> dict[str, object]:
    if cohort == "GSE281367":
        dependence = "independent"
        status = "real_external_donor_dynamic_program_projection"
    elif cohort == "GSE244832":
        dependence = "source_dependent"
        status = "real_internal_crossmodal_donor_dynamic_program_projection"
    else:
        raise ContractError(f"unsupported ATAC cohort: {cohort}")
    return {
        "adapter_id": f"{cohort.lower()}_snatac_dynamic",
        "adapter_root": f"protein_atac/adapters/{cohort.lower()}_snatac_dynamic",
        "dataset": cohort,
        "assay": "snATAC_donor_pseudobulk_promoter_accessibility",
        "source_publication": f"{cohort}_NCBI_GEO_public_record",
        "source_dependence": dependence,
        "biological_unit": "donor",
        "biological_unit_resolution": "resolved",
        "technical_unit": "donor_pseudobulk",
        "analysis_set_id": "hepatocyte_MASH_vs_NORMAL_program_accessibility",
        "contrast_or_exposure": "MASH_minus_NORMAL",
        "effect_unit": "donor_program_score_mean_difference",
        "cross_assay_comparable": False,
        "program_universe": "external_test_eligible",
        "min_genes_testable": 8,
        "min_retained_l1_weight": 0.2,
        "gate_expectation": "indeterminate",
        "registry_sha256": seal.registry_sha256,
        "ready_sha256": seal.ready_sha256,
        "status": status,
    }


def _base_program_identity(reg: Mapping[str, str], spec: Mapping[str, object], source_hash: str, producer: str, producer_hash: str) -> dict[str, object]:
    return {
        "release_id": RELEASE_ID,
        "registry_sha256": spec["registry_sha256"],
        "membership_sha256": reg["membership_sha256"],
        "program_uid": reg["program_uid"],
        "legacy_program_id": _legacy_id(reg),
        "program_label": reg["module_name"],
        "cell_type": reg["cell_type"],
        "dataset": spec["dataset"],
        "assay": spec["assay"],
        "source_publication": spec["source_publication"],
        "source_dependence": spec["source_dependence"],
        "analysis_set_id": spec["analysis_set_id"],
        "biological_unit": spec["biological_unit"],
        "biological_unit_resolution": spec["biological_unit_resolution"],
        "technical_unit": spec["technical_unit"],
        "contrast_or_exposure": spec["contrast_or_exposure"],
        "effect_unit": spec["effect_unit"],
        "matched_null_sd": "",
        "descriptive_effect_direction": "",
        "inferential_test_direction": "",
        "direction_agreement": False,
        "heterogeneity_statistic": "",
        "heterogeneity_df": "",
        "heterogeneity_pvalue": "",
        "negative_call_rule_id": "",
        "cross_assay_comparable": False,
        "producer": producer,
        "producer_sha256": producer_hash,
        "source_manifest_sha256": source_hash,
    }


def _build_protein(
    project_root: Path,
    root: Path,
    spec: Mapping[str, object],
    seal,
    selected: Mapping[str, Mapping[str, str]],
    membership: Mapping[str, Sequence[Mapping[str, str]]],
    producer_relative: str,
    producer_hash: str,
    created: str,
) -> None:
    v1 = project_root / "Analysis/Multimodal_Program_Projection/results/proteomics"
    hotspot = project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "hotspot"
    protein_de = v1 / "protein_de_adjusted.tsv"
    _, de_rows = read_tsv(protein_de, ("gene",))
    measured_genes = {row["gene"] for row in de_rows}
    _, accepted_rows = read_tsv(
        v1 / "module_protein_results.tsv",
        ("program_id", "n_measured", "retained_l1_weight", "testable"),
    )
    accepted = {row["program_id"]: row for row in accepted_rows}
    _, meta = read_tsv(v1 / "panel4c_metadata.tsv", ("sample_id", "group", "acquisition_batch"))
    if len(meta) != 58 or len({row["sample_id"] for row in meta}) != 58:
        raise ContractError("PXD051911 accepted metadata must contain 58 unique patient samples")

    root.mkdir(parents=True)
    source_hash = _source_manifest(
        project_root,
        root,
        (
            (project_root / "data/PXD051911/liver_protein_quant.txt", "public_DIA_protein_group_intensity", True, "not_recorded_legacy_input"),
            (project_root / "data/PXD051911/meta_data.txt", "public_patient_metadata", True, "not_recorded_legacy_input"),
            (protein_de, "accepted_group_blind_adjusted_protein_universe", False, "not_recorded_legacy_output"),
            (v1 / "module_protein_results.tsv", "accepted_v1_program_projection_for_compatibility", False, "not_recorded_legacy_output"),
            (v1 / "panel4c_metadata.tsv", "accepted_58_patient_design", False, "not_recorded_legacy_output"),
            (v1 / "panel4c_mrna_protein.tsv", "protected_fixed_25_descriptive_rows", False, "not_recorded_legacy_output"),
            (project_root / "Analysis/Multimodal_Program_Projection/config/panel4c_fixed_rows.tsv", "versioned_fixed_25_row_contract", False, "not_recorded_legacy_input"),
            (project_root / "Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R", "validated_v1_assay_native_producer", False, "not_applicable_code"),
            (hotspot / "READY", "sealed_v2_registry_ready", False, "not_applicable_seal"),
            (hotspot / "program_registry_v2.tsv", "sealed_v2_registry", False, "not_applicable_seal"),
            (hotspot / "program_membership_v2.tsv", "sealed_v2_membership", False, "not_applicable_seal"),
        ),
    )

    sample_rows = [
        {
            "release_id": RELEASE_ID,
            "dataset": spec["dataset"],
            "analysis_set_id": spec["analysis_set_id"],
            "biological_id": row["sample_id"],
            "technical_id": row["sample_id"],
            "biological_unit": spec["biological_unit"],
            "technical_unit": spec["technical_unit"],
            "include_primary": True,
            "gate_state": f"included_design_{row['group']}",
        }
        for row in meta
    ]
    ids = [row["sample_id"] for row in meta]
    design = {
        "release_id": RELEASE_ID,
        "dataset": spec["dataset"],
        "analysis_set_id": spec["analysis_set_id"],
        "contrast_or_exposure": spec["contrast_or_exposure"],
        "biological_unit": spec["biological_unit"],
        "biological_unit_resolution": spec["biological_unit_resolution"],
        "n_biological": 58,
        "technical_unit": spec["technical_unit"],
        "n_technical": 58,
        "design_status": "resolved_58_patients_but_selected_v2_programs_fail_observability",
        "biological_ids_sha256": canonical_id_hash(ids),
        "technical_ids_sha256": canonical_id_hash(ids),
    }
    mapping_rows = []
    testability_rows = []
    effect_rows = []
    sensitivity_rows = []
    for uid in sorted(selected):
        reg = selected[uid]
        by_symbol: dict[str, float] = {}
        for row in membership[uid]:
            symbol = row["mapped_symbol"]
            by_symbol[symbol] = by_symbol.get(symbol, 0.0) + float(row["original_l1_weight"])
        measured = sorted(set(by_symbol) & measured_genes)
        retained = sum(by_symbol[symbol] for symbol in measured)
        source_program_id = _v1_program_id(reg)
        if source_program_id not in accepted:
            raise ContractError(f"accepted protein row absent for {source_program_id}")
        old = accepted[source_program_id]
        if parse_int(old["n_measured"], f"protein.n_measured[{source_program_id}]") != len(measured):
            raise ContractError(f"v2/v1 measured protein subset differs for {source_program_id}")
        old_weight = parse_float(old["retained_l1_weight"], f"protein.weight[{source_program_id}]")
        if old_weight is None or abs(old_weight - retained) > 1e-12:
            raise ContractError(f"v2/v1 retained protein weight differs for {source_program_id}")
        if parse_bool(old["testable"], f"protein.testable[{source_program_id}]"):
            raise ContractError(f"selected protein program unexpectedly testable: {source_program_id}")
        reason = (
            f"prespecified_observability_fail;n_genes={len(measured)}<8_or_"
            f"retained_l1_weight={retained:.17g}<0.20;no_program_outcome_test"
        )
        mapping_rows.append(
            {
                "release_id": RELEASE_ID,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "program_uid": uid,
                "membership_sha256": reg["membership_sha256"],
                "n_source_genes": reg["n_source_genes"],
                "n_genes_measured": len(measured),
                "retained_l1_weight": retained,
                "mapping_status": "group_blind_DIA_unambiguous_single_gene_coverage_audited",
            }
        )
        testability_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_sha256": seal.registry_sha256,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "program_uid": uid,
                "membership_sha256": reg["membership_sha256"],
                "testable": False,
                "n_genes_measured": len(measured),
                "retained_l1_weight": retained,
                "testability_reason": reason,
                "evidence_state": "untestable",
            }
        )
        effect = _base_program_identity(reg, spec, source_hash, producer_relative, producer_hash)
        effect.update(
            {
                "n_biological": 58,
                "n_technical": 58,
                "estimate": "",
                "std_error": "",
                "interval_low": "",
                "interval_high": "",
                "interval_type": "not_applicable",
                "pvalue": "",
                "padj": "",
                "pvalue_method": "not_applicable",
                "multiplicity_family": "complete_two_program_family_zero_testable",
                "n_genes_measured": len(measured),
                "retained_l1_weight": retained,
                "testable": False,
                "testability_reason": reason,
                "direction_expected": reg["primary_direction"],
                "direction_observed": "",
                "sensitivity_sign_agree": False,
                "robustness_pass": False,
                "evidence_state": "untestable",
            }
        )
        effect_rows.append(effect)
        sensitivity_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_sha256": seal.registry_sha256,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "program_uid": uid,
                "membership_sha256": reg["membership_sha256"],
                "sensitivity_id": "not_run_program_unobservable",
                "estimate": "",
                "effect_unit": spec["effect_unit"],
                "direction_observed": "",
                "sign_agree": False,
                "status": "untestable",
            }
        )

    write_tsv(root / "sample_manifest.tsv", SCHEMAS["sample_manifest.tsv"], sample_rows)
    write_tsv(root / "gene_mapping_audit.tsv", SCHEMAS["gene_mapping_audit.tsv"], mapping_rows)
    write_tsv(root / "design_audit.tsv", SCHEMAS["design_audit.tsv"], [design])
    write_tsv(root / "program_testability.tsv", SCHEMAS["program_testability.tsv"], testability_rows)
    write_tsv(root / "per_sample_program_scores.tsv", SCHEMAS["per_sample_program_scores.tsv"], [])
    write_tsv(root / "program_effects.tsv", SCHEMAS["program_effects.tsv"], effect_rows)
    write_tsv(root / "sensitivity.tsv", SCHEMAS["sensitivity.tsv"], sensitivity_rows)
    shutil.copyfile(v1 / "panel4c_mrna_protein.tsv", root / "fixed_25_protein_context.tsv")
    execution_hash = _execution_manifest(
        root,
        (
            "sample_manifest.tsv",
            "gene_mapping_audit.tsv",
            "design_audit.tsv",
            "program_testability.tsv",
            "per_sample_program_scores.tsv",
            "program_effects.tsv",
            "sensitivity.tsv",
            "fixed_25_protein_context.tsv",
            "source_manifest.tsv",
        ),
        created,
    )
    _write_gate(
        root,
        spec,
        seal,
        source_hash,
        execution_hash,
        "untestable",
        "both_prespecified_v2_programs_fail_DIA_observability;fixed_25_rows_retained_descriptive_selection_conditioned",
        False,
    )


def _build_atac(
    project_root: Path,
    root: Path,
    spec: Mapping[str, object],
    seal,
    selected: Mapping[str, Mapping[str, str]],
    native_root: Path,
    producer_relative: str,
    producer_hash: str,
    created: str,
) -> None:
    cohort = str(spec["dataset"])
    _, native_effects_all = read_tsv(
        native_root / "dynamic_program_accessibility.tsv",
        (
            "release_id", "registry_sha256", "program_uid", "membership_sha256",
            "cohort", "n_measured", "retained_l1_weight", "testable", "effect",
            "std_error", "pvalue", "equal_effect", "leave_top_effect", "padj",
            "sensitivity_sign_agree", "robust", "direction_expected", "n_normal",
            "n_mash", "pvalue_method", "multiplicity_family",
        ),
    )
    native_effects = [row for row in native_effects_all if row["cohort"] == cohort]
    native_by_uid = {row["program_uid"]: row for row in native_effects}
    if len(native_by_uid) != 2 or set(native_by_uid) != set(selected):
        raise ContractError(f"{cohort} native ATAC output is not the complete two-program family")
    _, native_scores_all = read_tsv(
        native_root / "dynamic_program_scores.tsv",
        (
            "release_id", "registry_sha256", "program_uid", "membership_sha256",
            "cohort", "sample_id", "condition", "score", "equal_score", "leave_top_score",
        ),
    )
    native_scores = [row for row in native_scores_all if row["cohort"] == cohort]
    sample_condition: dict[str, str] = {}
    for row in native_scores:
        sample = row["sample_id"]
        condition = row["condition"]
        prior = sample_condition.setdefault(sample, condition)
        if prior != condition:
            raise ContractError(f"{cohort} donor condition drift for {sample}")
    expected_n = 12 if cohort == "GSE281367" else 14
    if len(sample_condition) != expected_n:
        raise ContractError(f"{cohort} expected {expected_n} donors, found {len(sample_condition)}")

    hotspot = project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "hotspot"
    if cohort == "GSE281367":
        counts = project_root / "Analysis/ATAC/Human_External/pseudobulk/hep_pseudobulk_counts_GSE281367.tsv.gz"
        metadata = project_root / "Analysis/ATAC/Human_External/pseudobulk/hep_pseudobulk_coldata_GSE281367.tsv"
    else:
        counts = project_root / "Analysis/ATAC/Human_Multiome/results/snapatac2/hep_pseudobulk_counts.tsv.gz"
        metadata = project_root / "Analysis/ATAC/Human_Multiome/metadata/donor_metadata_curated.tsv"
    root.mkdir(parents=True)
    source_hash = _source_manifest(
        project_root,
        root,
        (
            (counts, "public_accession_derived_donor_pseudobulk_peak_counts", True, "not_recorded_legacy_input"),
            (metadata, "public_accession_derived_donor_condition_metadata", True, "not_recorded_legacy_input"),
            (native_root / "input_manifest.tsv", "exact_bounded_producer_input_manifest_including_external_GTF", False, "not_applicable_candidate_output"),
            (native_root / "dynamic_program_accessibility.tsv", "bounded_v2_assay_native_effects", False, "not_applicable_candidate_output"),
            (native_root / "dynamic_program_scores.tsv", "bounded_v2_donor_scores_and_sensitivities", False, "not_applicable_candidate_output"),
            (native_root / "v1_selected_regression.tsv", "candidate_engine_v1_numeric_regression", False, "not_applicable_candidate_output"),
            (native_root / "producer_status.tsv", "bounded_v2_producer_status", False, "not_applicable_candidate_output"),
            (project_root / "Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2/10_run_atac_v2_projection.R", "bounded_candidate_ATAC_producer", False, "not_applicable_code"),
            (project_root / "Analysis/Multimodal_Program_Projection/scripts/03_atac_projection.R", "validated_v1_native_logic_source", False, "not_applicable_code"),
            (hotspot / "READY", "sealed_v2_registry_ready", False, "not_applicable_seal"),
            (hotspot / "program_registry_v2.tsv", "sealed_v2_registry", False, "not_applicable_seal"),
            (hotspot / "program_membership_v2.tsv", "sealed_v2_membership", False, "not_applicable_seal"),
        ),
    )

    sample_rows = [
        {
            "release_id": RELEASE_ID,
            "dataset": cohort,
            "analysis_set_id": spec["analysis_set_id"],
            "biological_id": sample,
            "technical_id": sample,
            "biological_unit": spec["biological_unit"],
            "technical_unit": spec["technical_unit"],
            "include_primary": True,
            "gate_state": f"included_{condition}",
        }
        for sample, condition in sorted(sample_condition.items())
    ]
    ids = sorted(sample_condition)
    design = {
        "release_id": RELEASE_ID,
        "dataset": cohort,
        "analysis_set_id": spec["analysis_set_id"],
        "contrast_or_exposure": spec["contrast_or_exposure"],
        "biological_unit": spec["biological_unit"],
        "biological_unit_resolution": spec["biological_unit_resolution"],
        "n_biological": expected_n,
        "technical_unit": spec["technical_unit"],
        "n_technical": expected_n,
        "design_status": "complete_unpaired_donor_label_permutation_design",
        "biological_ids_sha256": canonical_id_hash(ids),
        "technical_ids_sha256": canonical_id_hash(ids),
    }
    mapping_rows = []
    testability_rows = []
    effect_rows = []
    sensitivity_rows = []
    score_rows = []
    for uid in sorted(selected):
        reg = selected[uid]
        row = native_by_uid[uid]
        if row["release_id"] != RELEASE_ID or row["registry_sha256"] != seal.registry_sha256:
            raise ContractError(f"{cohort} native registry linkage drift for {uid}")
        if row["membership_sha256"] != reg["membership_sha256"]:
            raise ContractError(f"{cohort} native membership drift for {uid}")
        n_measured = parse_int(row["n_measured"], f"{cohort}.n_measured[{uid}]")
        retained = parse_float(row["retained_l1_weight"], f"{cohort}.weight[{uid}]")
        testable = parse_bool(row["testable"], f"{cohort}.testable[{uid}]")
        if not testable or n_measured is None or retained is None:
            raise ContractError(f"{cohort} selected v2 ATAC program is not testable: {uid}")
        robust = parse_bool(row["robust"], f"{cohort}.robust[{uid}]")
        sign_agree = parse_bool(row["sensitivity_sign_agree"], f"{cohort}.sign_agree[{uid}]")
        estimate = parse_float(row["effect"], f"{cohort}.effect[{uid}]")
        std_error = parse_float(row["std_error"], f"{cohort}.se[{uid}]")
        pvalue = parse_float(row["pvalue"], f"{cohort}.p[{uid}]")
        padj = parse_float(row["padj"], f"{cohort}.q[{uid}]")
        equal_effect = parse_float(row["equal_effect"], f"{cohort}.equal[{uid}]")
        leave_effect = parse_float(row["leave_top_effect"], f"{cohort}.leave[{uid}]")
        assert None not in (estimate, std_error, pvalue, padj, equal_effect, leave_effect)
        evidence_state = "robust" if robust else "indeterminate"
        direction = "positive" if estimate > 0 else "negative" if estimate < 0 else "zero"
        reason = "passed_8_gene_20pct_L1_observability"
        mapping_rows.append(
            {
                "release_id": RELEASE_ID,
                "dataset": cohort,
                "analysis_set_id": spec["analysis_set_id"],
                "program_uid": uid,
                "membership_sha256": reg["membership_sha256"],
                "n_source_genes": reg["n_source_genes"],
                "n_genes_measured": n_measured,
                "retained_l1_weight": retained,
                "mapping_status": "GENCODE_v49_plus_group_blind_promoter_activity_filter",
            }
        )
        testability_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_sha256": seal.registry_sha256,
                "dataset": cohort,
                "analysis_set_id": spec["analysis_set_id"],
                "program_uid": uid,
                "membership_sha256": reg["membership_sha256"],
                "testable": True,
                "n_genes_measured": n_measured,
                "retained_l1_weight": retained,
                "testability_reason": reason,
                "evidence_state": evidence_state,
            }
        )
        effect = _base_program_identity(reg, spec, source_hash, producer_relative, producer_hash)
        effect.update(
            {
                "n_biological": expected_n,
                "n_technical": expected_n,
                "estimate": estimate,
                "std_error": std_error,
                "interval_low": "",
                "interval_high": "",
                "interval_type": "none",
                "pvalue": pvalue,
                "padj": padj,
                "pvalue_method": row["pvalue_method"],
                "multiplicity_family": row["multiplicity_family"],
                "n_genes_measured": n_measured,
                "retained_l1_weight": retained,
                "testable": True,
                "testability_reason": reason,
                "direction_expected": row["direction_expected"],
                "direction_observed": direction,
                "descriptive_effect_direction": direction,
                "inferential_test_direction": direction,
                "direction_agreement": True,
                "sensitivity_sign_agree": sign_agree,
                "robustness_pass": robust,
                "evidence_state": evidence_state,
            }
        )
        effect_rows.append(effect)
        for sensitivity_id, sensitivity_estimate in (
            ("equal_weight", equal_effect),
            ("leave_highest_weight_gene_out", leave_effect),
        ):
            sensitivity_direction = (
                "positive" if sensitivity_estimate > 0 else "negative" if sensitivity_estimate < 0 else "zero"
            )
            sensitivity_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "registry_sha256": seal.registry_sha256,
                    "dataset": cohort,
                    "analysis_set_id": spec["analysis_set_id"],
                    "program_uid": uid,
                    "membership_sha256": reg["membership_sha256"],
                    "sensitivity_id": sensitivity_id,
                    "estimate": sensitivity_estimate,
                    "effect_unit": spec["effect_unit"],
                    "direction_observed": sensitivity_direction,
                    "sign_agree": sensitivity_direction == direction,
                    "status": "tested",
                }
            )
        uid_scores = [score for score in native_scores if score["program_uid"] == uid]
        if len(uid_scores) != expected_n:
            raise ContractError(f"{cohort} incomplete donor scores for {uid}")
        for score in uid_scores:
            score_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "registry_sha256": seal.registry_sha256,
                    "dataset": cohort,
                    "analysis_set_id": spec["analysis_set_id"],
                    "program_uid": uid,
                    "membership_sha256": reg["membership_sha256"],
                    "biological_id": score["sample_id"],
                    "technical_unit_count": 1,
                    "program_score": score["score"],
                    "score_unit": "weighted_z_standardized_donor_promoter_accessibility",
                }
            )

    write_tsv(root / "sample_manifest.tsv", SCHEMAS["sample_manifest.tsv"], sample_rows)
    write_tsv(root / "gene_mapping_audit.tsv", SCHEMAS["gene_mapping_audit.tsv"], mapping_rows)
    write_tsv(root / "design_audit.tsv", SCHEMAS["design_audit.tsv"], [design])
    write_tsv(root / "program_testability.tsv", SCHEMAS["program_testability.tsv"], testability_rows)
    write_tsv(root / "per_sample_program_scores.tsv", SCHEMAS["per_sample_program_scores.tsv"], score_rows)
    write_tsv(root / "program_effects.tsv", SCHEMAS["program_effects.tsv"], effect_rows)
    write_tsv(root / "sensitivity.tsv", SCHEMAS["sensitivity.tsv"], sensitivity_rows)
    write_tsv(
        root / "scoring_sensitivity_scores.tsv",
        (
            "release_id", "registry_sha256", "program_uid", "membership_sha256",
            "cohort", "sample_id", "condition", "score", "equal_score", "leave_top_score",
        ),
        native_scores,
    )
    execution_hash = _execution_manifest(
        root,
        (
            "sample_manifest.tsv",
            "gene_mapping_audit.tsv",
            "design_audit.tsv",
            "program_testability.tsv",
            "per_sample_program_scores.tsv",
            "program_effects.tsv",
            "sensitivity.tsv",
            "scoring_sensitivity_scores.tsv",
            "source_manifest.tsv",
        ),
        created,
    )
    _write_gate(
        root,
        spec,
        seal,
        source_hash,
        execution_hash,
        "indeterminate",
        "complete_two_program_donor_projection_with_exact_label_permutation_and_complete_BH_family",
        True,
    )


def _build_static_context(
    project_root: Path,
    root: Path,
    seal,
    selected: Mapping[str, Mapping[str, str]],
    membership: Mapping[str, Sequence[Mapping[str, str]]],
    created: str,
) -> None:
    v1 = project_root / "Analysis/Multimodal_Program_Projection/results/atac"
    _, hit_rows = read_tsv(
        v1 / "static_gwas_atac_hits.tsv",
        (
            "program_id", "cell_type", "gene_symbol", "study", "ancestry", "trait",
            "trait_scope", "source_locus", "chromosome", "locus", "variant_id",
            "recommended_pip", "chr_hg38", "pos_hg38",
        ),
    )
    _, accepted_summary = read_tsv(
        v1 / "static_gwas_atac_summary.tsv",
        (
            "trait_scope", "program_id", "n_source_loci", "n_study_loci", "n_variants",
            "n_genes", "n_studies", "n_ancestries", "n_registry_studies",
            "n_studies_with_credible_sets",
        ),
    )
    legacy_to_uid = {_v1_program_id(reg): uid for uid, reg in selected.items()}
    filtered = [row for row in hit_rows if row["program_id"] in legacy_to_uid]
    mapped_symbols = {
        uid: {row["mapped_symbol"] for row in membership[uid]}
        for uid in selected
    }
    output_hits = []
    for row in filtered:
        uid = legacy_to_uid[row["program_id"]]
        if row["gene_symbol"] not in mapped_symbols[uid]:
            raise ContractError(
                f"static ATAC accepted hit gene is not mapped in v2 membership: {uid}/{row['gene_symbol']}"
            )
        reg = selected[uid]
        output_hits.append(
            {
                "release_id": RELEASE_ID,
                "registry_sha256": seal.registry_sha256,
                "program_uid": uid,
                "membership_sha256": reg["membership_sha256"],
                "legacy_program_id": _legacy_id(reg),
                "cell_type": row["cell_type"],
                "gene_symbol": row["gene_symbol"],
                "study": row["study"],
                "ancestry": row["ancestry"],
                "trait": row["trait"],
                "trait_scope": row["trait_scope"],
                "source_locus": row["source_locus"],
                "chromosome": row["chromosome"],
                "locus": row["locus"],
                "variant_id": row["variant_id"],
                "recommended_pip": row["recommended_pip"],
                "chr_hg38": row["chr_hg38"],
                "pos_hg38": row["pos_hg38"],
                "source_dependence": "descriptive",
                "interpretation": "static_gene_level_finemapped_variant_in_lineage_open_promoter_not_program_inference",
            }
        )

    accepted_by_key = {
        (row["program_id"], row["trait_scope"]): row
        for row in accepted_summary
        if row["program_id"] in legacy_to_uid
    }
    context_rows = []
    for uid in sorted(selected):
        reg = selected[uid]
        legacy = _legacy_id(reg)
        source_program_id = _v1_program_id(reg)
        for scope in ("direct_disease_PDFF", "liver_enzyme"):
            rows = [row for row in output_hits if row["program_uid"] == uid and row["trait_scope"] == scope]
            counts = {
                "n_source_loci": len({row["source_locus"] for row in rows}),
                "n_study_loci": len({(row["study"], row["source_locus"]) for row in rows}),
                "n_variants": len({row["variant_id"] for row in rows}),
                "n_genes": len({row["gene_symbol"] for row in rows}),
                "n_studies": len({row["study"] for row in rows}),
                "n_ancestries": len({row["ancestry"] for row in rows}),
            }
            accepted = accepted_by_key.get((source_program_id, scope))
            if accepted is None:
                raise ContractError(f"accepted static ATAC summary missing {source_program_id}/{scope}")
            for field, value in counts.items():
                if parse_int(accepted[field], f"static.{field}[{source_program_id}/{scope}]") != value:
                    raise ContractError(f"static ATAC {field} does not rederive for {source_program_id}/{scope}")
            context_rows.append(
                {
                    "release_id": RELEASE_ID,
                    "registry_sha256": seal.registry_sha256,
                    "program_uid": uid,
                    "membership_sha256": reg["membership_sha256"],
                    "legacy_program_id": legacy,
                    "program_label": reg["module_name"],
                    "cell_type": reg["cell_type"],
                    "dataset": "35_study_Tier1_Tier2_GWAS_registry_plus_lineage_snATAC",
                    "assay": "fine_mapped_variant_in_lineage_open_promoter",
                    "source_dependence": "descriptive",
                    "trait_scope": scope,
                    **counts,
                    "n_registry_studies": accepted["n_registry_studies"],
                    "n_studies_with_credible_sets": accepted["n_studies_with_credible_sets"],
                    "effect_unit": "unique_program_genes_with_open_promoter_finemapped_variant",
                    "interpretation": "static_gene_link_context_not_dynamic_accessibility_or_program_association",
                }
            )

    root.mkdir(parents=True)
    source_hash = _source_manifest(
        project_root,
        root,
        (
            (v1 / "static_gwas_atac_hits.tsv", "accepted_source_native_static_gene_links", False, "not_recorded_legacy_output"),
            (v1 / "static_gwas_atac_summary.tsv", "accepted_source_native_static_program_summary", False, "not_recorded_legacy_output"),
            (v1 / "static_study_coverage.tsv", "prespecified_35_study_registry_coverage", False, "not_recorded_legacy_output"),
            (project_root / "Analysis/Multimodal_Program_Projection/scripts/03_atac_projection.R", "validated_v1_static_link_producer", False, "not_applicable_code"),
            (project_root / "GWAS/finemapping/results/credible_sets.csv", "fine_mapped_credible_set_input", False, "not_recorded_legacy_input"),
            (project_root / "GWAS/finemapping/config/gwas_trait_tier.tsv", "prespecified_trait_tier_registry", False, "not_recorded_legacy_input"),
            (project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "hotspot/program_registry_v2.tsv", "sealed_v2_registry", False, "not_applicable_seal"),
            (project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "hotspot/program_membership_v2.tsv", "sealed_v2_membership", False, "not_applicable_seal"),
        ),
    )
    write_tsv(root / "static_atac_context.tsv", STATIC_CONTEXT_COLUMNS, context_rows)
    write_tsv(root / "static_atac_hits.tsv", STATIC_HIT_COLUMNS, output_hits)
    shutil.copyfile(v1 / "static_study_coverage.tsv", root / "static_study_coverage.tsv")
    write_tsv(
        root / "source_native_status.tsv",
        ("release_id", "status", "n_programs", "n_context_rows", "n_hit_rows", "source_manifest_sha256"),
        [
            {
                "release_id": RELEASE_ID,
                "status": "descriptive_source_native_static_context_no_program_test",
                "n_programs": 2,
                "n_context_rows": len(context_rows),
                "n_hit_rows": len(output_hits),
                "source_manifest_sha256": source_hash,
            }
        ],
    )
    _execution_manifest(
        root,
        (
            "static_atac_context.tsv",
            "static_atac_hits.tsv",
            "static_study_coverage.tsv",
            "source_native_status.tsv",
            "source_manifest.tsv",
        ),
        created,
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--candidate-root", type=Path, default=None)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    candidate_root = args.candidate_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "spatial_context"
    )
    hotspot_root = project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "hotspot"
    native_root = candidate_root / "protein_atac_native/atac_v2"
    target = candidate_root / "protein_atac"
    try:
        if target.exists():
            raise ContractError(f"refusing to overwrite SP-INT-05 candidate: {target}")
        if not native_root.is_dir():
            raise ContractError("bounded v2 ATAC native output must exist before adapter assembly")
        seal = validate_hotspot_seal(hotspot_root)
        selected = program_universe(seal, "external_test_eligible")
        if len(selected) != 2:
            raise ContractError(f"expected two external-test programs, found {len(selected)}")
        membership = _mapped_membership(hotspot_root, selected)
        producer_path = Path(__file__).resolve()
        producer_relative = _project_relative(project_root, producer_path)
        producer_hash = sha256_file(producer_path)
        created = utc_now()

        staging = candidate_root / f".protein_atac.incomplete.{os.getpid()}"
        if staging.exists():
            raise ContractError(f"staging path already exists: {staging}")
        staging.mkdir(parents=True)
        specs = [_protein_spec(seal), _atac_spec(seal, "GSE281367"), _atac_spec(seal, "GSE244832")]
        _build_protein(
            project_root, staging / "adapters/pxd051911_diams_program", specs[0], seal,
            selected, membership, producer_relative, producer_hash, created,
        )
        _build_atac(
            project_root, staging / "adapters/gse281367_snatac_dynamic", specs[1], seal,
            selected, native_root, producer_relative, producer_hash, created,
        )
        _build_atac(
            project_root, staging / "adapters/gse244832_snatac_dynamic", specs[2], seal,
            selected, native_root, producer_relative, producer_hash, created,
        )
        _build_static_context(
            project_root, staging / "source_native/atac_static_context", seal,
            selected, membership, created,
        )
        write_tsv(staging / "adapter_registry.tsv", ADAPTER_REGISTRY_COLUMNS, specs)
        os.replace(staging, target)
        print("PASS: built isolated protein/ATAC candidate adapters and descriptive static context")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
