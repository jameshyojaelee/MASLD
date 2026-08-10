#!/usr/bin/env python3
"""Executable PASS-02--PASS-06 synthetic fixture suite."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any

import build_evidence_passport_bundle as builder
import adjudicate_gen_ensembl_identity as gen_identity
import generate_evidence_passports as contract
import validate_evidence_passports as validator


SCRIPT_PATH = Path(__file__).resolve()
FIXTURE_ROOT = contract.DEFAULT_FIXTURE_ROOT / "pass02_06"
INPUT_ROOT = FIXTURE_ROOT / "inputs"
POSITIVE_ROOT = FIXTURE_ROOT / "positive"
REPEAT_ROOT = FIXTURE_ROOT / "positive_repeat"


GEN_IDENTITY_PREFLIGHT_PATHS = {
    "GEN_IDENTITY_ADJUDICATION_READY": "GEN_IDENTITY_ADJUDICATION_READY",
    "gen_frozen_classes_ensembl_adjudicated.tsv": "gen_frozen_classes_ensembl_adjudicated.tsv",
    "gen_identity_adjudication_audit.tsv": "gen_identity_adjudication_audit.tsv",
    "gen_identity_quarantine.tsv": "gen_identity_quarantine.tsv",
}


def stage_identity_preflight(input_root: Path, output_root: Path) -> None:
    source_root = input_root / "gen_identity_adjudication"
    destination = output_root / "preflight" / "gen_identity_adjudication_v1"
    destination.mkdir(parents=True, exist_ok=True)
    for source_name, destination_name in GEN_IDENTITY_PREFLIGHT_PATHS.items():
        shutil.copy2(source_root / source_name, destination / destination_name)


def write_table(path: Path, rows: list[dict[str, Any]]) -> Path:
    if not rows:
        raise RuntimeError(f"fixture table cannot be empty: {path}")
    contract.atomic_write_tsv(path, rows, list(rows[0]))
    return path


def fixture_payloads(root: Path, variant: str = "positive") -> tuple[Path, Path]:
    root.mkdir(parents=True, exist_ok=True)
    genes = [
        ("ENSG00000000001", "GENE1"),
        ("ENSG00000000002", "GENE2"),
        ("ENSG00000000003", "GENE3"),
        ("ENSG00000000004", "GENE4"),
        ("ENSG00000000005", "GENE5"),
    ]
    identity = write_table(
        root / "gencode_fixture.tsv",
        [{"ensembl_base": ensembl, "gene_name": symbol} for ensembl, symbol in genes],
    )
    closure = write_table(
        root / "gen_terminal_closure.tsv",
        [
            {
                "release_id": "fixture-pass02-06-v1",
                "terminal_status": "coverage_limited_terminal",
                "context_rescue_authorized": "false",
                "negative_claim_authorized": "false",
            }
        ],
    )
    closure_hash = "0" * 64 if variant == "stale_gen" else contract.sha256_file(closure)
    gen_ready = write_table(
        root / "GEN_TERMINAL_CLOSURE_READY",
        [
            {
                "release_id": "fixture-pass02-06-v1",
                "status": "coverage_limited_terminal_validated",
                "closure_sha256": closure_hash,
                "context_rescue_authorized": "false",
                "negative_claim_authorized": "false",
            }
        ],
    )
    static_specs = [
        (
            "GENE1",
            genes[0][0],
            "true",
            "true",
            "true",
            "true",
            "convergent",
            "0.8",
            "0.001",
            "0.91",
        ),
        (
            "GENE2",
            genes[1][0],
            "true",
            "true",
            "false",
            "true",
            "disease_state_only",
            "0.5",
            "0.02",
            "0.10",
        ),
        (
            "GENE3",
            genes[2][0],
            "true",
            "true",
            "true",
            "false",
            "genetic_only",
            "0.1",
            "0.8",
            "0.85",
        ),
        (
            "GENE4",
            genes[3][0],
            "false",
            "false",
            "false",
            "false",
            "indeterminate_not_jointly_testable",
            "",
            "",
            "",
        ),
    ]
    gen_rows = []
    for (
        symbol,
        ensembl,
        bulk_tested,
        genetic_tested,
        genetic,
        state,
        static_class,
        lfc,
        fdr,
        pp4,
    ) in static_specs:
        gen_rows.append(
            {
                "gene_symbol": symbol,
                "ensembl_bulk": ensembl,
                "ensembl_genetic": ensembl if genetic_tested == "true" else "",
                "bulk_tested": bulk_tested,
                "primary_genetic_map_tested": genetic_tested,
                "joint_testable": str(
                    bulk_tested == "true" and genetic_tested == "true"
                ).lower(),
                "primary_genetic": genetic,
                "established_state_associated": state,
                "static_class": static_class,
                "bulk_row_count": "1" if bulk_tested == "true" else "0",
                "bulk_logFC": lfc,
                "bulk_t": "2.0" if lfc else "",
                "bulk_AveExpr": "5.0" if lfc else "",
                "bulk_treat_fdr": fdr,
                "genetic_row_count": "1" if genetic_tested == "true" else "0",
                "coloc_best_susie_pp4": pp4,
                "coloc_best_abf_pp4": "",
                "driving_gwas": "FIXTURE_MASLD" if genetic_tested == "true" else "",
                "driving_trait": "MASLD diagnosis" if genetic_tested == "true" else "",
                "context_annotation": "not_evaluated_preflight",
            }
        )
    gen_classes = write_table(root / "gen_frozen_classes.tsv", gen_rows)
    gen_identity_products = gen_identity.adjudicate(
        gen_classes,
        identity,
        root / "gen_identity_adjudication",
        "fixture-pass02-06-v1",
        True,
    )
    accepted_ensembl = (
        f"{genes[4][0]};{genes[3][0]}" if variant == "multi_ensembl" else genes[4][0]
    )
    accepted = write_table(
        root / "accepted_gene_evidence.tsv",
        [
            {
                "source_input_row_id": "ACCEPTED001",
                "ensembl_id": accepted_ensembl,
                "symbol": "GENE5",
                "primary_evidence_class": "unresolved",
                "evidence_domain": "chromatin"
                if variant == "gene_domain_scope"
                else "transcriptomics",
                "assay": "fixture_scRNA_pseudobulk",
                "dataset_id": "FIXTURE_SCRNA",
                "phenotype": "MASH",
                "context": "hepatocyte transcript abundance",
                "biological_unit": "donor",
                "contrast_or_exposure": "MASH minus control",
                "effect_unit": "log2 fold change",
                "estimate": "0.1",
                "standard_error": "0.2",
                "ci_lower": "-0.3",
                "ci_upper": "0.5",
                "p_value": "0.6",
                "q_value": "0.7",
                "direction": "unknown",
                "testability_state": "testable",
                "testability_reason": "testable",
                "call_state": "indeterminate",
                "negative_call_rule_id": "",
                "negative_decision_boundary": "",
                "negative_margin": "",
                "negative_call_passed": "",
                "gate_id": "FIXTURE_SCRNA_GATE",
                "n_biological_units": "5",
                "n_technical_units": "10",
                "allowed_wording": "Complete synthetic transcriptomic test with no positive or negative call.",
                "limitation": "Synthetic PASS-02 fixture only.",
            }
        ],
    )
    membership_hashes = {"P_HEP_1": "1" * 64, "P_HEP_2": "2" * 64}
    semantic_rows = []
    for index, uid in enumerate(membership_hashes, start=1):
        semantic_rows.append(
            {
                "release_id": "fixture-pass02-06-v1",
                "program_uid": uid,
                "cell_type": "hepatocytes",
                "module": str(index),
                "module_name": f"Fixture hepatocyte program {index}",
                "membership_sha256": membership_hashes[uid],
                "registry_sha256": "3" * 64,
                "primary_estimable": "TRUE",
                "primary_qvalue": "0.2",
                "primary_hc3_qvalue": "0.3",
                "stability_median": "0.8",
                "primary_selected": "FALSE",
                "selected_unstable": "FALSE",
                "selected_hc3_fragile": "FALSE",
                "robust_display": "FALSE",
                "external_test_eligible": "TRUE",
                "legacy_tested_negative": "TRUE",
                "legacy_registry_state": "tested_nonsignificant",
                "legacy_negative_semantics_deprecated": "true",
                "adjudicated_state": "indeterminate_nonconfirmatory",
                "selection_tier": "not_selected_q_nonsignificant",
                "tested_negative_authorized": "false",
                "evidence_scope": "internal_donor_level_stage_association_only",
                "external_outcome_state": "not_adjudicated_here",
                "adjudication_reason": "nonsignificance is nonconfirmatory",
                "semantic_contract": "fixture_semantics_v1",
            }
        )
    semantics = write_table(root / "hotspot_semantics.tsv", semantic_rows)
    hotspot_ready = write_table(
        root / "SEMANTIC_ADJUDICATION_READY",
        [
            {
                "release_id": "fixture-pass02-06-v1",
                "status": "semantic_adjudication_validated",
                "adjudication_sha256": contract.sha256_file(semantics),
                "n_tested_negative_authorized": "0",
            }
        ],
    )
    memberships = write_table(
        root / "hotspot_membership.tsv",
        [
            {
                "cell_type": "hepatocytes",
                "module": "1",
                "source_gene": "GENE1",
                "canonical_gene": "GENE1",
                "source_weight": "2.0",
                "mapped_symbol": "GENE1",
                "mapped_symbol_status": "gencode_v49_unique_symbol_confirmed",
                "original_l1_weight": "0.6",
                "canonical_weight_text": "2.0",
                "membership_sha256": membership_hashes["P_HEP_1"],
                "program_uid": "P_HEP_1",
            },
            {
                "cell_type": "hepatocytes",
                "module": "2",
                "source_gene": "GENE2",
                "canonical_gene": "GENE2",
                "source_weight": "1.0",
                "mapped_symbol": "GENE2",
                "mapped_symbol_status": "gencode_v49_unique_symbol_confirmed",
                "original_l1_weight": "0.4",
                "canonical_weight_text": "1.0",
                "membership_sha256": membership_hashes["P_HEP_2"],
                "program_uid": "P_HEP_2",
            },
            {
                "cell_type": "hepatocytes",
                "module": "1",
                "source_gene": genes[2][0],
                "canonical_gene": genes[2][0],
                "source_weight": "0.5",
                "mapped_symbol": genes[2][0],
                "mapped_symbol_status": "gencode_v49_unambiguous_ensembl_to_symbol",
                "original_l1_weight": "0.1",
                "canonical_weight_text": "0.5",
                "membership_sha256": membership_hashes["P_HEP_1"],
                "program_uid": "P_HEP_1",
            },
        ],
    )
    plan13_rows = []
    for index, (uid, state, estimate, padj, dependence) in enumerate(
        [
            ("P_HEP_1", "robust", "0.4", "0.01", "independent"),
            ("P_HEP_2", "source_dependent", "0.2", "0.20", "source_dependent"),
        ],
        start=1,
    ):
        plan13_rows.append(
            {
                "release_id": "fixture-pass02-06-v1",
                "registry_sha256": "3" * 64,
                "membership_sha256": membership_hashes[uid],
                "program_uid": uid,
                "legacy_program_id": f"hepatocytes:{index}",
                "program_label": f"Fixture hepatocyte program {index}",
                "cell_type": "hepatocytes",
                "dataset": "YakubovskyFixture",
                "assay": "synthetic_spatial",
                "source_publication": "Fixture publication",
                "source_dependence": dependence,
                "analysis_set_id": "fixture_lipid",
                "biological_unit": "donor",
                "biological_unit_resolution": "resolved",
                "n_biological": "4",
                "technical_unit": "spot",
                "n_technical": "400",
                "contrast_or_exposure": "local lipid",
                "effect_unit": "standardized donor slope",
                "estimate": estimate,
                "std_error": "0.1",
                "interval_low": "",
                "interval_high": "",
                "interval_type": "none",
                "pvalue": padj,
                "padj": padj,
                "pvalue_method": "fixture_block_permutation",
                "multiplicity_family": "two_programs",
                "n_genes_measured": "2",
                "retained_l1_weight": "0.8",
                "testable": "TRUE",
                "testability_reason": "testable",
                "direction_expected": "positive",
                "direction_observed": "positive",
                "sensitivity_sign_agree": "TRUE",
                "robustness_pass": "TRUE" if state == "robust" else "FALSE",
                "cross_assay_comparable": "FALSE",
                "evidence_state": state,
                "producer": "fixture",
                "producer_sha256": "4" * 64,
                "source_manifest_sha256": "5" * 64,
                "source_adapter_id": "fixture",
                "evidence_role": "fixture",
                "claim_scope": "program_only",
                "uncertainty_semantics": "fixture",
                "import_method": "fixture",
                "source_ready_path": "fixture",
                "source_ready_sha256": "6" * 64,
                "source_effect_path": "fixture",
                "source_effect_sha256": "7" * 64,
                "source_effect_row_sha256": contract.stable_sha256({"uid": uid}),
                "figure4_include": "TRUE",
                "figure4_interpretation": "program context only",
            }
        )
    if variant == "semantic_spoof":
        plan13_rows[1]["evidence_state"] = "tested_negative"
    plan13_effects = write_table(root / "integrated_program_effects.tsv", plan13_rows)
    dataset_verdict = write_table(
        root / "figure4_dataset_verdict.tsv",
        [
            {
                "release_id": "fixture-pass02-06-v1",
                "dataset": "YakubovskyFixture",
                "assay": "synthetic_spatial",
                "analysis_set_id": "fixture_lipid",
                "source_dependence": "independent",
                "biological_unit": "donor",
                "n_biological": "4",
                "technical_unit": "spot",
                "n_technical": "400",
                "effect_unit": "standardized donor slope",
                "n_programs": "2",
                "n_robust": "1",
                "n_tested_negative": "0",
                "n_untestable": "0",
                "n_not_applicable": "0",
                "n_skipped": "0",
                "dataset_gate": "pass",
                "figure4_verdict": "retain_complete_program_panel",
                "figure4_role": "lipid_context",
                "interpretation": "Synthetic program-grain context.",
            },
            {
                "release_id": "fixture-pass02-06-v1",
                "dataset": "GSE287826",
                "assay": "GeoMx",
                "analysis_set_id": "donor_replication",
                "source_dependence": "independent",
                "biological_unit": "donor",
                "n_biological": "",
                "technical_unit": "AOI",
                "n_technical": "30",
                "effect_unit": "not applicable",
                "n_programs": "0",
                "n_robust": "0",
                "n_tested_negative": "0",
                "n_untestable": "0",
                "n_not_applicable": "0",
                "n_skipped": "1",
                "dataset_gate": "skipped_no_donor_key",
                "figure4_verdict": "skip",
                "figure4_role": "conditional_replication",
                "interpretation": "No authoritative public donor key; dataset-level skip only.",
            },
        ],
    )
    plan13_ready = write_table(
        root / "PLAN13_READY",
        [
            {
                "release_id": "fixture-pass02-06-v1",
                "status": "validated_pass02_06_fixture",
                "mode": "fixture_pass02_06",
                "plan13_complete": "FALSE"
                if variant == "incomplete_plan13"
                else "TRUE",
                "integrated_effects_sha256": contract.sha256_file(plan13_effects),
                "figure4_verdict_sha256": contract.sha256_file(dataset_verdict),
            }
        ],
    )
    myojin_ready = write_table(
        root / "PHASE_C_VALIDATED",
        [
            {
                "release_id": "fixture-pass02-06-v1",
                "status": "phase_c_release_validated",
                "external_outcomes_read": "TRUE",
            }
        ],
    )
    myojin_verdict = write_table(
        root / "myojin_verdict.tsv",
        [
            {
                "main_figure_eligible": "FALSE",
                "source_gate_pass": "TRUE",
                "verdict_reason_codes": "class_branch_failed;program_branch_failed",
            }
        ],
    )
    myojin_class = write_table(
        root / "myojin_class.tsv",
        [
            {
                "contrast": "omnibus_evidence_class",
                "estimate": "",
                "HC3_SE": "",
                "CI95_low": "",
                "CI95_high": "",
                "permutation_p": "0.5",
                "BH_q": "",
            }
        ],
    )
    myojin_program = write_table(
        root / "myojin_program.tsv",
        [
            {
                "program_uid": "P_HEP_1",
                "testable": "TRUE",
                "testability_reason": "testable",
                "signed_effect": "0.05",
                "empirical_p": "0.8",
                "BH_q": "0.8",
            }
        ],
    )
    specs = [
        (
            "IDENTITY",
            "PASS01",
            "gencode_identity",
            "gencode_v49_identity_v1",
            "reference",
            identity,
            "identity_only",
            "",
            "independent",
            "FIXTURE_GENCODE",
        ),
        (
            "GEN_CLOSURE",
            "GEN",
            "gen_terminal_closure",
            "gen_terminal_closure_v1",
            "gate",
            closure,
            "gate_only",
            "",
            "independent",
            "FIXTURE_GEN",
        ),
        (
            "GEN_READY",
            "GEN",
            "gen_terminal_ready",
            "gen_terminal_ready_v1",
            "gate",
            gen_ready,
            "gate_only",
            "",
            "independent",
            "FIXTURE_GEN",
        ),
        (
            "GEN_RAW",
            "GEN",
            "gen_frozen_classes_raw_identity_source",
            "gen_raw_identity_source_v1",
            "provenance",
            gen_classes,
            "source_provenance_only",
            "GEN_READY",
            "independent",
            "FIXTURE_GEN",
        ),
        (
            "GEN_ID_READY",
            "PASS03",
            "gen_identity_adjudication_ready",
            "gen_identity_adjudication_ready_v1",
            "gate",
            gen_identity_products["ready"],
            "gate_only",
            "GEN_READY",
            "independent",
            "FIXTURE_GEN",
        ),
        (
            "GEN_ID_AUDIT",
            "PASS03",
            "gen_identity_adjudication_audit",
            "gen_identity_adjudication_audit_v1",
            "provenance",
            gen_identity_products["audit"],
            "source_provenance_only",
            "GEN_ID_READY",
            "independent",
            "FIXTURE_GEN",
        ),
        (
            "GEN_ID_QUAR",
            "PASS03",
            "gen_identity_adjudication_quarantine",
            "gen_identity_adjudication_quarantine_v1",
            "provenance",
            gen_identity_products["quarantine"],
            "source_provenance_only",
            "GEN_ID_READY",
            "independent",
            "FIXTURE_GEN",
        ),
        (
            "GEN_CLASSES",
            "GEN",
            "gen_frozen_classes_ensembl_adjudicated",
            "gen_adjudicated_classes_v1",
            "gene",
            gen_identity_products["adjudicated"],
            "gene_evidence",
            "GEN_ID_READY",
            "independent",
            "FIXTURE_GEN",
        ),
        (
            "ACCEPTED_GENE",
            "PASS03",
            "accepted_gene_evidence",
            "accepted_gene_evidence_v1",
            "gene",
            accepted,
            "gene_evidence",
            "GEN_READY",
            "independent",
            "FIXTURE_SCRNA",
        ),
        (
            "HOTSPOT_READY",
            "PLAN20",
            "hotspot_semantic_ready",
            "hotspot_semantic_ready_v1",
            "gate",
            hotspot_ready,
            "gate_only",
            "",
            "independent",
            "FIXTURE_HOTSPOT",
        ),
        (
            "HOTSPOT_SEM",
            "PLAN20",
            "hotspot_program_semantics",
            "hotspot_program_semantics_v1",
            "program",
            semantics,
            "program_context",
            "HOTSPOT_READY",
            "partially_dependent",
            "FIXTURE_HOTSPOT",
        ),
        (
            "HOTSPOT_MEM",
            "PLAN20",
            "hotspot_program_membership",
            "hotspot_program_membership_v1",
            "program_membership",
            memberships,
            "program_membership_only",
            "HOTSPOT_READY",
            "partially_dependent",
            "FIXTURE_HOTSPOT",
        ),
        (
            "PLAN13_READY",
            "PLAN13",
            "plan13_terminal_ready",
            "plan13_terminal_ready_v1",
            "gate",
            plan13_ready,
            "gate_only",
            "",
            "independent",
            "YakubovskyFixture;GSE287826",
        ),
        (
            "PLAN13_EFFECTS",
            "PLAN13",
            "plan13_program_context",
            "plan13_program_context_v1",
            "program",
            plan13_effects,
            "program_context",
            "PLAN13_READY",
            "partially_dependent",
            "YakubovskyFixture;GSE287826",
        ),
        (
            "PLAN13_DATASETS",
            "PLAN13",
            "plan13_dataset_status",
            "dataset_status_v1",
            "dataset",
            dataset_verdict,
            "dataset_status_only",
            "PLAN13_READY",
            "partially_dependent",
            "YakubovskyFixture;GSE287826",
        ),
        (
            "MYOJIN_READY",
            "PLAN40",
            "myojin_terminal_ready",
            "myojin_terminal_ready_v1",
            "gate",
            myojin_ready,
            "gate_only",
            "",
            "independent",
            "MYOJIN_FIXTURE",
        ),
        (
            "MYOJIN_ASSAY",
            "PLAN40",
            "myojin_assay_non_support",
            "myojin_assay_verdict_v1",
            "assay",
            myojin_verdict,
            "assay_status_only",
            "MYOJIN_READY",
            "independent",
            "MYOJIN_FIXTURE",
        ),
        (
            "MYOJIN_CLASS",
            "PLAN40",
            "myojin_class_non_support",
            "myojin_class_nonsupport_v1",
            "assay",
            myojin_class,
            "assay_status_only",
            "MYOJIN_READY",
            "independent",
            "MYOJIN_FIXTURE",
        ),
        (
            "MYOJIN_PROGRAM",
            "PLAN40",
            "myojin_program_non_support",
            "myojin_program_nonsupport_v1",
            "assay",
            myojin_program,
            "assay_status_only",
            "MYOJIN_READY",
            "independent",
            "MYOJIN_FIXTURE",
        ),
    ]
    selection_rows = []
    for (
        input_id,
        workstream,
        role,
        adapter,
        grain,
        path,
        scope,
        terminal,
        provenance,
        datasets,
    ) in specs:
        producer_path = (
            gen_identity.SCRIPT_PATH
            if adapter
            in {
                "gen_identity_adjudication_ready_v1",
                "gen_identity_adjudication_audit_v1",
                "gen_identity_adjudication_quarantine_v1",
                "gen_adjudicated_classes_v1",
            }
            else SCRIPT_PATH
        )
        parent_ids = terminal
        if input_id == "GEN_ID_READY":
            parent_ids = "GEN_READY;GEN_RAW;IDENTITY"
        selection_rows.append(
            {
                "input_id": input_id,
                "workstream_id": workstream,
                "artifact_role": role,
                "adapter_id": adapter,
                "artifact_grain": grain,
                "artifact_path": str(path.resolve()),
                "artifact_sha256": contract.sha256_file(path),
                "artifact_bytes": path.stat().st_size,
                "producer_script": str(
                    producer_path.relative_to(contract.PROJECT_ROOT)
                ),
                "producer_sha256": contract.sha256_file(producer_path),
                "source_release_id": "fixture-pass02-06-v1",
                "gate_verdict": "include",
                "terminal_gate_input_id": terminal,
                "allowed_claim_wording": "Synthetic source-faithful fixture result.",
                "limitation": "Synthetic PASS-02--PASS-06 fixture only.",
                "source_datasets": datasets,
                "source_cohorts": "FIXTURE_COHORT",
                "source_publication": "Synthetic PASS-02--PASS-06 fixture",
                "source_url": "fixture://pass02-06",
                "biological_unit": "source-defined biological unit",
                "n_biological_units": "4",
                "n_technical_units": "4",
                "provenance_state": provenance,
                "claim_scope": scope,
                "inclusion_destination": "passport_only",
                "selection_decision": "accepted",
                "parent_input_ids": parent_ids,
                "fixture_only": "true",
            }
        )
    selection = root / "passport_input_selection.tsv"
    contract.atomic_write_tsv(
        selection, selection_rows, contract.SELECTION_REQUIRED_COLUMNS
    )
    signature = root / "passport_input_selection.signature.json"
    contract.atomic_write_json(
        signature,
        {
            "attestation_version": contract.SELECTION_ATTESTATION_VERSION,
            "selection_sha256": contract.sha256_file(selection),
            "signed_by": "PASS02_06_SYNTHETIC_FIXTURE",
            "signed_at_utc": "2026-08-08T00:00:00Z",
            "decision_register_id": "FIXTURE_PASS02_06",
            "authority_document": "docs/plans/2026-08-07_paper_program/50_EVIDENCE_PASSPORTS_AND_PORTAL.md",
            "fixture_only": True,
            "coordinator_attested": False,
            "analysis_release_id": "fixture-pass02-06-v1",
            "gen_terminal_closure_input_id": "GEN_CLOSURE",
            "gen_terminal_ready_input_id": "GEN_READY",
            "plan13_terminal_ready_input_id": "PLAN13_READY",
        },
    )
    return selection, signature


def compare_tree(first: Path, second: Path) -> list[str]:
    excluded = {
        "passport_validation_report.tsv",
        "passport_plan60_handoff.tsv",
        validator.TERMINAL_PROVENANCE_FILENAME,
        "PASS06_VALIDATED",
        "FIXTURE_PASS06_VALIDATED",
    }

    def included(path: Path, root: Path) -> bool:
        relative = path.relative_to(root)
        return (
            path.is_file()
            and relative.as_posix() not in excluded
            and "superseded_bundles" not in relative.parts
        )

    first_files = sorted(
        path.relative_to(first).as_posix()
        for path in first.rglob("*")
        if included(path, first)
    )
    second_files = sorted(
        path.relative_to(second).as_posix()
        for path in second.rglob("*")
        if included(path, second)
    )
    if first_files != second_files:
        return ["inventory"]
    return [
        name
        for name in first_files
        if contract.sha256_file(first / name) != contract.sha256_file(second / name)
    ]


def write_manual_v2_fixture(bundle: Path, output: Path) -> tuple[Path, Path]:
    """Create a production-shaped manual-v2 record bound to fixture UI bytes."""

    selection = bundle / "passport_input_selection.tsv"
    selection_signature = json.loads(
        (bundle / "passport_input_selection.signature.json").read_text(
            encoding="utf-8"
        )
    )
    ui_paths = {
        "ui_index_sha256": bundle / "portal_candidate/index.html",
        "ui_contract_sha256": bundle / "portal_candidate/ui_contract.json",
        "ui_source_manifest_sha256": (
            bundle / "portal_candidate/ui_source_manifest.tsv"
        ),
        "ui_review_manifest_sha256": (
            bundle / "portal_candidate/review/review_manifest.tsv"
        ),
    }
    ui_hashes = {
        field: contract.sha256_file(path) for field, path in ui_paths.items()
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    row = {
        "review_id": "FIXTURE_MANUAL_V2",
        "analysis_release_id": selection_signature["analysis_release_id"],
        "selection_sha256": contract.sha256_file(selection),
        "reviewer": "FIXTURE_REVIEWER",
        "reviewed_at_utc": "2026-08-08T00:00:00Z",
        **ui_hashes,
        "call_states_reviewed": "true",
        "provenance_states_reviewed": "true",
        "hero_genes_reviewed": "THRB;HKDC1;GLP1R;MTARC1",
        "visible_boundary_pass": "true",
        "visible_testability_pass": "true",
        "visible_source_dependence_pass": "true",
        "visible_falsifier_pass": "true",
        "decision": "accepted",
    }
    contract.atomic_write_tsv(output, [row], contract.MANUAL_ACCEPTANCE_COLUMNS)
    signature = output.with_name("passport_manual_acceptance.signature.json")
    contract.atomic_write_json(
        signature,
        {
            "attestation_version": contract.MANUAL_ACCEPTANCE_ATTESTATION_VERSION,
            "acceptance_sha256": contract.sha256_file(output),
            **ui_hashes,
            "signed_by": row["reviewer"],
            "signed_at_utc": row["reviewed_at_utc"],
            "decision_register_id": "FIXTURE_MANUAL_V2",
            "fixture_only": False,
        },
    )
    return output, signature


def mutate_embedded_native_record(
    bundle: Path, payload_key: str, field: str = "wording"
) -> None:
    """Make payload and lower-level hashes self-consistent but scientifically stale."""

    index = bundle / "portal_candidate/index.html"
    text = index.read_text(encoding="utf-8")
    match = re.search(
        r'(<script id="passport-data" type="application/json">)(.*?)(</script>)',
        text,
        flags=re.DOTALL,
    )
    if match is None:
        raise RuntimeError("fixture embedded payload is absent")
    payload = json.loads(match.group(2))
    records = payload[payload_key]
    if not records:
        raise RuntimeError(f"fixture {payload_key} is empty")
    records[0][field] += " Tampered but internally rehashed."
    embedded = validator.passport_ui.stable_json(payload).replace("<", "\\u003c")
    index.write_text(
        text[: match.start(2)] + embedded + text[match.end(2) :],
        encoding="utf-8",
    )
    contract_path = bundle / "portal_candidate/ui_contract.json"
    ui_contract = json.loads(contract_path.read_text(encoding="utf-8"))
    ui_contract["embedded_payload_sha256"] = hashlib.sha256(
        embedded.encode("utf-8")
    ).hexdigest()
    contract.atomic_write_json(contract_path, ui_contract)
    review_path = bundle / "portal_candidate/review/review_manifest.tsv"
    review_rows = validator.read_tsv(review_path)
    for row in review_rows:
        row["source_index_sha256"] = contract.sha256_file(index)
        row["ui_contract_sha256"] = contract.sha256_file(contract_path)
    contract.atomic_write_tsv(
        review_path,
        review_rows,
        list(validator.passport_ui_renderer.REVIEW_MANIFEST_COLUMNS),
    )


def run_suite() -> Path:
    FIXTURE_ROOT.mkdir(parents=True, exist_ok=True)
    gen_identity.run_self_test()
    selection, _ = fixture_payloads(INPUT_ROOT / "positive")
    stage_identity_preflight(INPUT_ROOT / "positive", POSITIVE_ROOT)
    stage_identity_preflight(INPUT_ROOT / "positive", REPEAT_ROOT)
    builder.build_from_selection(selection, POSITIVE_ROOT, fixture_mode=True)
    builder.build_from_selection(selection, REPEAT_ROOT, fixture_mode=True)
    differences = compare_tree(POSITIVE_ROOT, REPEAT_ROOT)
    if differences:
        raise contract.PassportContractError(
            "NONDETERMINISTIC_PRODUCTION_FIXTURE", str(differences)
        )
    report, assembly, _, signature = validator.validate_production_bundle(
        POSITIVE_ROOT, fixture_mode=True
    )
    superseded_fixture = POSITIVE_ROOT / "superseded_bundles/prior/PASS06_VALIDATED"
    superseded_fixture.parent.mkdir(parents=True, exist_ok=True)
    superseded_fixture.write_text("archived fixture seal\n", encoding="utf-8")
    validator.finalize_pass06(
        POSITIVE_ROOT,
        report,
        assembly,
        signature,
        fixture_mode=True,
        manual_acceptance=None,
    )
    source_nodes = validator.read_tsv(POSITIVE_ROOT / "passport_source_nodes.tsv")
    local_uris = [
        row["uri"]
        for row in source_nodes
        if row["node_type"] in {"derived_object", "frozen_release"}
    ]
    if not local_uris or any(
        uri.startswith("/") or "/gpfs/" in uri for uri in local_uris
    ):
        raise contract.PassportContractError(
            "NONRELOCATABLE_SOURCE_GRAPH", str(local_uris[:5])
        )
    manifest_rows = {
        row["relative_path"]: row
        for row in validator.read_tsv(POSITIVE_ROOT / "passport_release_manifest.tsv")
    }
    if any(path.startswith("superseded_bundles/") for path in manifest_rows):
        raise contract.PassportContractError(
            "SUPERSEDED_BUNDLE_IN_ACTIVE_MANIFEST", "archival subtree"
        )
    expected_producers = {
        "passport_input_selection.tsv": (
            builder.SELECTION_LOGICAL_PRODUCER_ID, builder.SELECTION_PRODUCER_PATH
        ),
        "passport_input_selection.signature.json": (
            builder.ATTESTATION_LOGICAL_PRODUCER_ID,
            builder.ATTESTATION_PRODUCER_PATH,
        ),
        "passport_gate_status.tsv": (
            builder.VALIDATOR_LOGICAL_PRODUCER_ID, builder.VALIDATOR_PRODUCER_PATH
        ),
        "passport_build_status.tsv": (
            builder.VALIDATOR_LOGICAL_PRODUCER_ID, builder.VALIDATOR_PRODUCER_PATH
        ),
        "preflight/gen_identity_adjudication_v1/GEN_IDENTITY_ADJUDICATION_READY": (
            builder.GEN_IDENTITY_LOGICAL_PRODUCER_ID,
            builder.GEN_IDENTITY_PRODUCER_PATH,
        ),
        "preflight/gen_identity_adjudication_v1/gen_frozen_classes_ensembl_adjudicated.tsv": (
            builder.GEN_IDENTITY_LOGICAL_PRODUCER_ID,
            builder.GEN_IDENTITY_PRODUCER_PATH,
        ),
        "preflight/gen_identity_adjudication_v1/gen_identity_adjudication_audit.tsv": (
            builder.GEN_IDENTITY_LOGICAL_PRODUCER_ID,
            builder.GEN_IDENTITY_PRODUCER_PATH,
        ),
        "preflight/gen_identity_adjudication_v1/gen_identity_quarantine.tsv": (
            builder.GEN_IDENTITY_LOGICAL_PRODUCER_ID,
            builder.GEN_IDENTITY_PRODUCER_PATH,
        ),
        "portal_candidate/index.html": (
            builder.UI_GENERATOR_LOGICAL_PRODUCER_ID,
            builder.UI_GENERATOR_PRODUCER_PATH,
        ),
        "portal_candidate/ui_contract.json": (
            builder.UI_GENERATOR_LOGICAL_PRODUCER_ID,
            builder.UI_GENERATOR_PRODUCER_PATH,
        ),
        "portal_candidate/ui_source_manifest.tsv": (
            builder.UI_GENERATOR_LOGICAL_PRODUCER_ID,
            builder.UI_GENERATOR_PRODUCER_PATH,
        ),
        "portal_candidate/review/overview.png": (
            builder.UI_RENDERER_LOGICAL_PRODUCER_ID,
            builder.UI_RENDERER_PRODUCER_PATH,
        ),
        "portal_candidate/review/THRB.png": (
            builder.UI_RENDERER_LOGICAL_PRODUCER_ID,
            builder.UI_RENDERER_PRODUCER_PATH,
        ),
        "portal_candidate/review/HKDC1.png": (
            builder.UI_RENDERER_LOGICAL_PRODUCER_ID,
            builder.UI_RENDERER_PRODUCER_PATH,
        ),
        "portal_candidate/review/GLP1R.png": (
            builder.UI_RENDERER_LOGICAL_PRODUCER_ID,
            builder.UI_RENDERER_PRODUCER_PATH,
        ),
        "portal_candidate/review/MTARC1.png": (
            builder.UI_RENDERER_LOGICAL_PRODUCER_ID,
            builder.UI_RENDERER_PRODUCER_PATH,
        ),
        "portal_candidate/review/review_manifest.tsv": (
            builder.UI_RENDERER_LOGICAL_PRODUCER_ID,
            builder.UI_RENDERER_PRODUCER_PATH,
        ),
    }
    for relative_path, (producer_id, producer_path) in expected_producers.items():
        row = manifest_rows.get(relative_path)
        if row is None or row["producer"] != producer_id:
            raise contract.PassportContractError(
                "MANIFEST_PRODUCER_ATTRIBUTION",
                f"{relative_path}: expected {producer_id}; observed {row}",
            )
        if row["producer_sha256"] != contract.sha256_file(producer_path):
            raise contract.PassportContractError(
                "MANIFEST_PRODUCER_HASH", relative_path
            )
    renderer_environments = {
        row["environment"]
        for relative, row in manifest_rows.items()
        if relative.startswith("portal_candidate/review/")
    }
    if len(renderer_environments) != 1 or not all(
        token in next(iter(renderer_environments))
        for token in (
            "renderer:pinned_puppeteer_chrome",
            "browser_sha256=",
            f"mesa={validator.passport_ui_renderer.EXPECTED_MESA_MODULE}",
            "fonts=system_fonts_unfrozen",
        )
    ):
        raise contract.PassportContractError(
            "MANIFEST_RENDERER_ENVIRONMENT", repr(renderer_environments)
        )

    def bundle_state(root: Path) -> dict[str, tuple[str, int, int]]:
        return {
            path.relative_to(root).as_posix(): (
                contract.sha256_file(path),
                path.stat().st_size,
                path.stat().st_mtime_ns,
            )
            for path in root.rglob("*")
            if path.is_file()
        }

    state_before = bundle_state(POSITIVE_ROOT)
    validator.verify_sealed_bundle(POSITIVE_ROOT)
    state_after = bundle_state(POSITIVE_ROOT)
    if state_before != state_after:
        raise contract.PassportContractError(
            "POST_SEAL_VERIFIER_MUTATED_BUNDLE", str(POSITIVE_ROOT)
        )

    terminal_failures = []
    with tempfile.TemporaryDirectory(
        prefix="relocated-passport-", dir=FIXTURE_ROOT
    ) as temporary:
        relocation_root = Path(temporary) / "relocated"
        shutil.copytree(POSITIVE_ROOT, relocation_root)
        validator.verify_sealed_bundle(relocation_root)

        ui_failures = []

        def expect_ui_failure(name: str, mutate, expected_code: str) -> None:
            target = Path(temporary) / f"ui_{name}"
            shutil.copytree(POSITIVE_ROOT, target)
            mutate(target)
            try:
                validator.check_candidate_ui(
                    target, fixture_mode=True, sealed=True
                )
            except validator.PassportValidationError as exc:
                if exc.code != expected_code:
                    raise contract.PassportContractError(
                        "UI_NEGATIVE_FIXTURE_WRONG_FAILURE",
                        f"{name}: expected {expected_code}; observed {exc.code}: {exc.message}",
                    ) from exc
                ui_failures.append((name, expected_code))
            else:
                raise contract.PassportContractError(
                    "UI_NEGATIVE_FIXTURE_FALSE_PASS", name
                )

        expect_ui_failure(
            "missing_index",
            lambda path: (path / "portal_candidate/index.html").unlink(),
            "UI_FILE_MISSING",
        )

        def drift_embedded_payload(path: Path) -> None:
            index = path / "portal_candidate/index.html"
            content = index.read_text(encoding="utf-8")
            token = '"alphabeticalDefault":true'
            if token not in content:
                raise RuntimeError("fixture embedded payload token is absent")
            index.write_text(
                content.replace(token, '"alphabeticalDefault":false', 1),
                encoding="utf-8",
            )

        expect_ui_failure(
            "embedded_payload_drift",
            drift_embedded_payload,
            "UI_EMBEDDED_PAYLOAD_HASH",
        )

        def drift_ui_source_hash(path: Path) -> None:
            source_manifest = path / "portal_candidate/ui_source_manifest.tsv"
            rows = validator.read_tsv(source_manifest)
            rows[0]["sha256"] = "0" * 64
            contract.atomic_write_tsv(
                source_manifest,
                rows,
                list(validator.passport_ui.SOURCE_MANIFEST_COLUMNS),
            )

        expect_ui_failure(
            "source_table_hash_drift",
            drift_ui_source_hash,
            "UI_SOURCE_HASH",
        )
        def drift_review_image(path: Path) -> None:
            with (path / "portal_candidate/review/THRB.png").open("ab") as handle:
                handle.write(b"tamper")

        expect_ui_failure(
            "review_image_drift",
            drift_review_image,
            "UI_REVIEW_IMAGE_HASH",
        )

        def drift_hero_review_row(path: Path) -> None:
            review = path / "portal_candidate/review/review_manifest.tsv"
            rows = validator.read_tsv(review)
            for row in rows:
                if row["review_id"] == "passport_ui_review:GLP1R":
                    row["gene_symbol"] = "THRB"
            contract.atomic_write_tsv(
                review,
                rows,
                list(validator.passport_ui_renderer.REVIEW_MANIFEST_COLUMNS),
            )

        expect_ui_failure(
            "hero_review_identity_drift",
            drift_hero_review_row,
            "UI_REVIEW_MANIFEST_SEMANTICS",
        )

        expect_ui_failure(
            "gene_self_consistent_payload_drift",
            lambda path: mutate_embedded_native_record(path, "genes", "claim"),
            "UI_GENE_SEMANTIC_DRIFT",
        )
        expect_ui_failure(
            "program_context_self_consistent_payload_drift",
            lambda path: mutate_embedded_native_record(
                path, "programContextRecords"
            ),
            "UI_PROGRAM_CONTEXT_SEMANTIC_DRIFT",
        )
        expect_ui_failure(
            "assay_grain_self_consistent_payload_drift",
            lambda path: mutate_embedded_native_record(path, "assayGrainRecords"),
            "UI_ASSAY_GRAIN_SEMANTIC_DRIFT",
        )

        chrome = validator.passport_ui_renderer._chrome_path(None)  # noqa: SLF001
        native_route_results = []
        for route, expected_view, expected_heading in (
            ("program-context", "program-context", "Program Context"),
            ("assay-records", "assay-records", "Assay and class status"),
        ):
            command = [
                *validator.passport_ui_renderer._browser_command(  # noqa: SLF001
                    chrome, 1440, 1200
                ),
                "--dump-dom",
                (POSITIVE_ROOT / "portal_candidate/index.html").as_uri()
                + f"?view={route}",
            ]
            rendered = subprocess.run(
                command,
                check=True,
                capture_output=True,
                text=True,
                timeout=60,
            ).stdout
            rendered_main = re.search(
                r'<main id="app"[^>]*>(.*?)</main>', rendered, flags=re.DOTALL
            )
            if (
                f'data-view="{expected_view}"' not in rendered
                or rendered_main is None
                or expected_heading not in rendered_main.group(1)
                or "Gene passport filters" in rendered_main.group(1)
            ):
                raise contract.PassportContractError(
                    "UI_NATIVE_ROUTE_FALSE_PASS", route
                )
            native_route_results.append(route)

        manual_failures = []
        manual_root = Path(temporary) / "manual_v2"
        manual, _ = write_manual_v2_fixture(
            POSITIVE_ROOT, manual_root / "passport_manual_acceptance.tsv"
        )
        analysis_release_id = json.loads(
            (POSITIVE_ROOT / "passport_input_selection.signature.json").read_text(
                encoding="utf-8"
            )
        )["analysis_release_id"]
        selection_sha = contract.sha256_file(
            POSITIVE_ROOT / "passport_input_selection.tsv"
        )
        validator.validate_manual_acceptance(
            manual, analysis_release_id, selection_sha, POSITIVE_ROOT
        )

        def expect_manual_failure(
            name: str, manual_path: Path, bundle: Path, expected_code: str
        ) -> None:
            try:
                validator.validate_manual_acceptance(
                    manual_path, analysis_release_id, selection_sha, bundle
                )
            except validator.PassportValidationError as exc:
                if exc.code != expected_code:
                    raise contract.PassportContractError(
                        "MANUAL_V2_NEGATIVE_WRONG_FAILURE",
                        f"{name}: expected {expected_code}; observed {exc.code}",
                    ) from exc
                manual_failures.append((name, expected_code))
            else:
                raise contract.PassportContractError(
                    "MANUAL_V2_NEGATIVE_FALSE_PASS", name
                )

        ui_drift_bundle = Path(temporary) / "manual_ui_drift"
        shutil.copytree(POSITIVE_ROOT, ui_drift_bundle)
        with (ui_drift_bundle / "portal_candidate/index.html").open(
            "a", encoding="utf-8"
        ) as handle:
            handle.write("\n")
        expect_manual_failure(
            "same_selection_changed_ui",
            manual,
            ui_drift_bundle,
            "MANUAL_ACCEPTANCE_UI_HASH",
        )

        schema_manual_root = Path(temporary) / "manual_schema"
        schema_manual, schema_signature = write_manual_v2_fixture(
            POSITIVE_ROOT,
            schema_manual_root / "passport_manual_acceptance.tsv",
        )
        schema_payload = json.loads(schema_signature.read_text(encoding="utf-8"))
        schema_payload["unexpected_field"] = "must fail"
        contract.atomic_write_json(schema_signature, schema_payload)
        expect_manual_failure(
            "signature_exact_key_mismatch",
            schema_manual,
            POSITIVE_ROOT,
            "MANUAL_ACCEPTANCE_SIGNATURE_SCHEMA",
        )

        coverage_bundle = Path(temporary) / "manual_hero_coverage"
        shutil.copytree(POSITIVE_ROOT, coverage_bundle)
        coverage_genes = validator.pd.read_parquet(
            coverage_bundle / "passport_gene_index.parquet", engine="pyarrow"
        ).sort_values(["symbol", "ensembl_id"], kind="stable")
        nonhero_id = str(coverage_genes.iloc[-1]["passport_id"])
        coverage_evidence = validator.pd.read_parquet(
            coverage_bundle / "passport_evidence_long.parquet", engine="pyarrow"
        )
        nonhero_rows = coverage_evidence.index[
            coverage_evidence["passport_id"].astype(str) == nonhero_id
        ]
        if len(nonhero_rows) == 0:
            raise contract.PassportContractError(
                "MANUAL_HERO_COVERAGE_FIXTURE", nonhero_id
            )
        coverage_evidence.loc[nonhero_rows[0], "call_state"] = "tested_negative"
        contract.atomic_write_parquet(
            coverage_bundle / "passport_evidence_long.parquet",
            contract.coerce_dataframe(
                coverage_evidence,
                contract.PRODUCTION_PARQUET_SCHEMAS[
                    "passport_evidence_long.parquet"
                ],
            ),
            contract.PRODUCTION_PARQUET_SCHEMAS["passport_evidence_long.parquet"],
        )
        expect_manual_failure(
            "hero_state_coverage",
            manual,
            coverage_bundle,
            "MANUAL_HERO_STATE_COVERAGE",
        )

        no_overwrite_bundle = Path(temporary) / "manual_no_overwrite"
        shutil.copytree(REPEAT_ROOT, no_overwrite_bundle)
        no_overwrite_report, no_overwrite_assembly, _, no_overwrite_signature = (
            validator.validate_production_bundle(
                no_overwrite_bundle, fixture_mode=True
            )
        )
        no_overwrite_manual, _ = write_manual_v2_fixture(
            no_overwrite_bundle,
            Path(temporary) / "manual_no_overwrite_stage/passport_manual_acceptance.tsv",
        )
        (no_overwrite_bundle / "passport_manual_acceptance.tsv").write_text(
            "stale partial decision\n", encoding="utf-8"
        )
        try:
            validator.finalize_pass06(
                no_overwrite_bundle,
                no_overwrite_report,
                no_overwrite_assembly,
                no_overwrite_signature,
                fixture_mode=False,
                manual_acceptance=no_overwrite_manual,
            )
        except validator.PassportValidationError as exc:
            if exc.code != "MANUAL_ACCEPTANCE_DESTINATION_EXISTS":
                raise
            manual_failures.append(
                ("finalize_no_overwrite", "MANUAL_ACCEPTANCE_DESTINATION_EXISTS")
            )
        else:
            raise contract.PassportContractError(
                "MANUAL_ACCEPTANCE_OVERWRITE_FALSE_PASS", str(no_overwrite_bundle)
            )

        relocated_tool_root = Path(temporary) / "relocated_portal_tools"
        relocated_tool_root.mkdir()
        relocated_sources = {
            *builder.manifest_producer_paths().values(),
            contract.SCRIPT_PATH,
        }
        for source in relocated_sources:
            shutil.copy2(source, relocated_tool_root / source.name)
        relocated_build = Path(temporary) / "relocated_producer_build"
        stage_identity_preflight(INPUT_ROOT / "positive", relocated_build)
        relocated_environment = os.environ.copy()
        relocated_environment["MASLD_PROJECT_ROOT"] = str(contract.PROJECT_ROOT)
        relocated_environment["PYTHONDONTWRITEBYTECODE"] = "1"
        subprocess.run(
            [
                sys.executable,
                "-B",
                str(relocated_tool_root / "build_evidence_passport_bundle.py"),
                "--selection",
                str(selection),
                "--output-root",
                str(relocated_build),
                "--fixture-mode",
            ],
            check=True,
            capture_output=True,
            text=True,
            timeout=600,
            env=relocated_environment,
        )
        subprocess.run(
            [
                sys.executable,
                "-B",
                str(relocated_tool_root / "validate_evidence_passports.py"),
                "--bundle",
                str(relocated_build),
                "--pass02-06-fixture",
                "--finalize-pass06",
            ],
            check=True,
            capture_output=True,
            text=True,
            timeout=600,
            env=relocated_environment,
        )
        subprocess.run(
            [
                sys.executable,
                "-B",
                str(relocated_tool_root / "validate_evidence_passports.py"),
                "--bundle",
                str(relocated_build),
                "--verify-sealed",
            ],
            check=True,
            capture_output=True,
            text=True,
            timeout=600,
            env=relocated_environment,
        )
        logical_paths = {
            logical: relocated_tool_root / Path(logical).name
            for logical in builder.manifest_producer_paths()
        }
        for row in validator.read_tsv(
            relocated_build / "passport_release_manifest.tsv"
        ):
            physical = logical_paths.get(row["producer"])
            if physical is None or row["producer_sha256"] != contract.sha256_file(
                physical
            ):
                raise contract.PassportContractError(
                    "RELOCATED_LOGICAL_PRODUCER_FALSE_PASS",
                    f"{row['relative_path']}:{row['producer']}",
                )
        for row in validator.read_tsv(
            relocated_build / validator.TERMINAL_PROVENANCE_FILENAME
        ):
            physical = logical_paths.get(row["producer"])
            if physical is None or row["producer_sha256"] != contract.sha256_file(
                physical
            ):
                raise contract.PassportContractError(
                    "RELOCATED_TERMINAL_PRODUCER_FALSE_PASS",
                    f"{row['artifact_role']}:{row['producer']}",
                )

        def expect_terminal_failure(name: str, mutate, expected_code: str) -> None:
            target = Path(temporary) / name
            shutil.copytree(POSITIVE_ROOT, target)
            mutate(target)
            try:
                validator.verify_sealed_bundle(target)
            except validator.PassportValidationError as exc:
                if exc.code != expected_code:
                    raise contract.PassportContractError(
                        "TERMINAL_NEGATIVE_FIXTURE_WRONG_FAILURE",
                        f"{name}: expected {expected_code}; observed {exc.code}: {exc.message}",
                    ) from exc
                terminal_failures.append((name, expected_code))
            else:
                raise contract.PassportContractError(
                    "TERMINAL_NEGATIVE_FIXTURE_FALSE_PASS", name
                )

        def append_text(path: Path, relative: str) -> None:
            with (path / relative).open("a", encoding="utf-8") as handle:
                handle.write("#tamper\n")

        external_root = Path(temporary) / "external_same_byte_artifacts"
        external_root.mkdir()

        def replace_with_external_symlink(path: Path, relative: str) -> None:
            source = path / relative
            external = external_root / f"{path.name}-{Path(relative).name}"
            shutil.copy2(source, external)
            source.unlink()
            source.symlink_to(external)

        expect_terminal_failure(
            "symlinked_payload",
            lambda path: replace_with_external_symlink(
                path, "passport_gene_index.parquet"
            ),
            "TERMINAL_PAYLOAD_SYMLINK",
        )
        expect_terminal_failure(
            "symlinked_ui_review",
            lambda path: replace_with_external_symlink(
                path, "portal_candidate/review/HKDC1.png"
            ),
            "TERMINAL_PAYLOAD_SYMLINK",
        )
        expect_terminal_failure(
            "symlinked_terminal_artifact",
            lambda path: replace_with_external_symlink(
                path, "passport_plan60_handoff.tsv"
            ),
            "TERMINAL_ARTIFACT_SYMLINK",
        )

        internal_root = external_root / "portal_export"
        shutil.copytree(POSITIVE_ROOT / "portal_export", internal_root)

        def replace_internal_directory_with_symlink(path: Path) -> None:
            internal = path / "portal_export"
            shutil.rmtree(internal)
            internal.symlink_to(internal_root, target_is_directory=True)

        expect_terminal_failure(
            "symlinked_internal_component",
            replace_internal_directory_with_symlink,
            "TERMINAL_INTERNAL_SYMLINK",
        )

        symlinked_root = Path(temporary) / "symlinked_bundle_root"
        symlinked_root.symlink_to(POSITIVE_ROOT, target_is_directory=True)
        try:
            validator.verify_sealed_bundle(symlinked_root)
        except validator.PassportValidationError as exc:
            if exc.code != "TERMINAL_BUNDLE_ROOT_SYMLINK":
                raise contract.PassportContractError(
                    "TERMINAL_NEGATIVE_FIXTURE_WRONG_FAILURE",
                    "symlinked_bundle_root: expected TERMINAL_BUNDLE_ROOT_SYMLINK; "
                    f"observed {exc.code}: {exc.message}",
                ) from exc
            terminal_failures.append(
                ("symlinked_bundle_root", "TERMINAL_BUNDLE_ROOT_SYMLINK")
            )
        else:
            raise contract.PassportContractError(
                "TERMINAL_NEGATIVE_FIXTURE_FALSE_PASS", "symlinked_bundle_root"
            )

        expect_terminal_failure(
            "tampered_handoff",
            lambda path: append_text(path, "passport_plan60_handoff.tsv"),
            "TERMINAL_PROVENANCE_ARTIFACT_HASH",
        )
        expect_terminal_failure(
            "tampered_terminal_provenance",
            lambda path: append_text(path, validator.TERMINAL_PROVENANCE_FILENAME),
            "TERMINAL_PROVENANCE_HASH",
        )
        expect_terminal_failure(
            "tampered_validation_report",
            lambda path: append_text(path, "passport_validation_report.tsv"),
            "TERMINAL_SEAL_HASH",
        )
        expect_terminal_failure(
            "atomic_temp_residue",
            lambda path: (path / ".passport_source_edges.tsv.fixture.tmp").write_text(
                "incomplete\n", encoding="utf-8"
            ),
            "TERMINAL_TEMP_RESIDUE",
        )

        def rewrite_absolute_handoff(path: Path) -> None:
            handoff_path = path / "passport_plan60_handoff.tsv"
            handoff_rows = validator.read_tsv(handoff_path)
            handoff_rows[0]["bundle_path"] = "/nonrelocatable/live/path"
            contract.atomic_write_tsv(
                handoff_path,
                handoff_rows,
                contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_plan60_handoff.tsv"],
            )
            terminal_path = path / validator.TERMINAL_PROVENANCE_FILENAME
            terminal_rows = validator.read_tsv(terminal_path)
            for row in terminal_rows:
                if row["artifact_role"] == "plan60_handoff":
                    row["sha256"] = contract.sha256_file(handoff_path)
                    row["bytes"] = str(handoff_path.stat().st_size)
            contract.atomic_write_tsv(
                terminal_path,
                terminal_rows,
                contract.PASS06_TERMINAL_PROVENANCE_COLUMNS,
            )
            report_path = path / "passport_validation_report.tsv"
            report_rows = validator.read_tsv(report_path)
            for row in report_rows:
                if row["check_id"] == validator.TERMINAL_PROVENANCE_CHECK_ID:
                    detail = json.loads(row["detail"])
                    detail["terminal_provenance_sha256"] = contract.sha256_file(
                        terminal_path
                    )
                    row["detail"] = json.dumps(
                        detail, sort_keys=True, separators=(",", ":")
                    )
            contract.atomic_write_tsv(
                report_path,
                report_rows,
                contract.PRODUCTION_ONLY_TSV_SCHEMAS["passport_validation_report.tsv"],
            )
            seal_path = path / "FIXTURE_PASS06_VALIDATED"
            seal_rows = validator.read_tsv(seal_path)
            seal_rows[0]["validation_report_sha256"] = contract.sha256_file(report_path)
            contract.atomic_write_tsv(
                seal_path, seal_rows, contract.PASS06_TERMINAL_COLUMNS
            )

        expect_terminal_failure(
            "absolute_handoff_uri",
            rewrite_absolute_handoff,
            "TERMINAL_HANDOFF_SEMANTICS",
        )
    expectations = {
        "semantic_spoof": "TESTED_NEGATIVE_RULE_REQUIRED",
        "multi_ensembl": "IDENTITY_QUARANTINE_NONEMPTY",
        "incomplete_plan13": "PLAN13_INCOMPLETE",
        "stale_gen": "GEN_TERMINAL_STALE",
        "gene_domain_scope": "GENE_DOMAIN_OUT_OF_CANDIDATE_SCOPE",
    }
    results = [
        {
            "check": "gen_identity_outcome_blind_adjudication",
            "status": "pass",
            "detail": "unique v49 identity accepted; irreducible multi-ID, absent-v49, and discordant duplicate calls quarantined",
        },
        {
            "check": "positive_pass02_06",
            "status": "pass",
            "detail": "all adapters and PASS06 fixture validation passed",
        },
        {
            "check": "deterministic_rebuild",
            "status": "pass",
            "detail": "byte-identical preterminal bundles",
        },
        {
            "check": "zero_tested_negative",
            "status": "pass",
            "detail": "complete bundle accepted with zero tested-negative rows",
        },
        {
            "check": "release_relative_source_graph",
            "status": "pass",
            "detail": "derived-object and frozen-release URIs contain no live absolute GPFS path",
        },
        {
            "check": "file_level_producer_provenance",
            "status": "pass",
            "detail": "selection, attestation, and validator-mutated files bind their actual producer hashes",
        },
        {
            "check": "gen_preflight_producer_provenance",
            "status": "pass",
            "detail": "all four GEN identity preflight products bind adjudicate_gen_ensembl_identity.py and its exact hash",
        },
        {
            "check": "post_seal_relocatable_read_only",
            "status": "pass",
            "detail": "relocated sealed bundle verifies without changing any byte, size, or mtime",
        },
        {
            "check": "terminal_chain_negative_fixtures",
            "status": "pass",
            "detail": ";".join(f"{name}:{code}" for name, code in terminal_failures),
        },
        {
            "check": "local_candidate_ui",
            "status": "pass",
            "detail": (
                "self-contained alphabetical UI, identity-safe navigation, source-table "
                "semantics, five review renders, and adversarial failures passed: "
                + ";".join(f"{name}:{code}" for name, code in ui_failures)
            ),
        },
        {
            "check": "ui_native_grain_routes",
            "status": "pass",
            "detail": (
                "browser-executed routes reached native program and assay views: "
                + ";".join(native_route_results)
            ),
        },
        {
            "check": "manual_acceptance_v2",
            "status": "pass",
            "detail": (
                "valid UI-byte-bound attestation passed; adversarial failures: "
                + ";".join(f"{name}:{code}" for name, code in manual_failures)
            ),
        },
        {
            "check": "relocated_logical_producers",
            "status": "pass",
            "detail": (
                "relocated frozen portal tools built and validated with stable logical "
                "producer IDs bound to the executing copied bytes"
            ),
        },
        {
            "check": "superseded_bundle_isolation",
            "status": "pass",
            "detail": "archival subtrees never enter the active release manifest",
        },
    ]
    for variant, expected in expectations.items():
        variant_selection, _ = fixture_payloads(INPUT_ROOT / variant, variant=variant)
        try:
            builder.build_from_selection(
                variant_selection,
                FIXTURE_ROOT / f"negative_{variant}",
                fixture_mode=True,
            )
        except contract.PassportContractError as exc:
            if exc.code != expected:
                raise contract.PassportContractError(
                    "NEGATIVE_FIXTURE_WRONG_FAILURE",
                    f"{variant}: expected {expected}; observed {exc.code}: {exc.message}",
                ) from exc
            results.append(
                {
                    "check": f"negative_{variant}",
                    "status": "pass",
                    "detail": f"rejected with {exc.code}",
                }
            )
        else:
            raise contract.PassportContractError("NEGATIVE_FIXTURE_FALSE_PASS", variant)
    unsigned_dir = INPUT_ROOT / "unsigned_selection"
    unsigned_selection, unsigned_signature = fixture_payloads(unsigned_dir)
    rows = contract.read_tsv_rows(unsigned_selection)
    for row in rows:
        row["fixture_only"] = "false"
    contract.atomic_write_tsv(
        unsigned_selection, rows, contract.SELECTION_REQUIRED_COLUMNS
    )
    signature_payload = json.loads(unsigned_signature.read_text(encoding="utf-8"))
    signature_payload.update(
        selection_sha256=contract.sha256_file(unsigned_selection),
        fixture_only=False,
        coordinator_attested=False,
    )
    contract.atomic_write_json(unsigned_signature, signature_payload)
    try:
        contract.validate_signed_selection(unsigned_selection, allow_fixture=False)
    except contract.PassportContractError as exc:
        if exc.code != "UNSIGNED_COORDINATOR_SELECTION":
            raise
        results.append(
            {
                "check": "negative_unsigned_selection",
                "status": "pass",
                "detail": "rejected with UNSIGNED_COORDINATOR_SELECTION",
            }
        )
    else:
        raise contract.PassportContractError(
            "UNSIGNED_SELECTION_FALSE_PASS", str(unsigned_selection)
        )
    report_path = FIXTURE_ROOT / "pass02_06_self_test_report.tsv"
    contract.atomic_write_tsv(report_path, results, ["check", "status", "detail"])
    return report_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    report = run_suite()
    print(f"PASS: PASS-02--PASS-06 fixtures validated; report={report}")


if __name__ == "__main__":
    try:
        main()
    except contract.PassportContractError as exc:
        raise SystemExit(str(exc)) from exc
