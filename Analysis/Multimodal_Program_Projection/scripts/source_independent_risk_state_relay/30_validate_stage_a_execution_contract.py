#!/usr/bin/env python3
"""Independently validate the outcome-blind Stage-A execution contract."""

from __future__ import annotations

import csv
import json

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


EXPECTED_SCHEMA_FIELDS = {
    "pilot_power_audit_template.tsv": {
        "platform_id", "pilot_id", "blinded_pilot", "scientific_outcomes_opened",
        "transfer_route", "challenge_id", "challenge_exposure_mode",
        "challenge_duration_value", "challenge_duration_unit", "early_time_value", "intermediate_time_value",
        "late_time_value", "stage_a_backgrounds",
        "stage_a_differentiations_per_background_condition", "familywise_alpha",
        "target_power", "power_scope", "power_variance_estimate", "power_effect_size",
        "viability_threshold", "maturation_threshold",
        "composition_deviation_threshold", "reviewer_1", "reviewer_2",
        "review_concordant",
    },
    "experimental_unit_manifest_template.tsv": {
        "biological_unit_id", "background_id", "differentiation_id", "platform_id",
        "source_lineage", "independent_background", "independent_differentiation",
    },
    "assay_registry_template.tsv": {
        "assay_id", "assay_role", "modality", "model_family", "normalization",
        "required_arm_ids", "required_time_roles", "required_challenge_roles",
        "biological_unit_definition", "technical_unit", "outcome_file_path",
        "outcomes_present", "frozen_before_unblinding", "reviewer_1", "reviewer_2",
        "review_concordant",
    },
    "phenotype_registry_template.tsv": {
        "phenotype_id", "assay_id", "phenotype_class", "expected_direction",
        "model_family", "multiplicity_family", "selected_before_outcomes",
        "reviewer_1", "reviewer_2", "review_concordant",
    },
    "randomization_manifest_template.tsv": {
        "sample_id", "blinded_label", "biological_unit_id", "background_id",
        "differentiation_id", "target_uid", "guide_id", "arm_id", "time_role",
        "challenge_role", "challenge_id", "assay_id", "randomized_allocation",
        "unblinding_status", "exclusion_status",
    },
    "execution_source_manifest_template.tsv": {
        "source_id", "source_path", "size_bytes", "sha256", "source_role",
        "scientific_outcome",
    },
    "dual_review_signoff_template.tsv": {
        "review_scope", "reviewer_1", "reviewer_2", "review_concordant",
        "scientific_outcomes_opened", "signed_utc",
    },
}


def header(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter="\t")
        return next(reader)


def main() -> None:
    seal_path = CANDIDATE_ROOT / "STAGE_A_EXECUTION_CONTRACT_SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_outcome_blind_stage_a_execution_contract":
        raise RuntimeError("Invalid Stage-A execution-contract status")
    for field in [
        "experimental_targets_frozen", "actual_stage_a_execution_frozen",
        "scientific_outcomes_inspected", "stage_b_design_frozen",
    ]:
        if seal.get(field) is not False:
            raise RuntimeError(f"Outcome-blind execution contract improperly sets {field}")
    if seal.get("governing_plan_integrity_bound") is not False:
        raise RuntimeError("Mutable plan prose must not be integrity-bound")
    if seal.get("required_primary_phenotypes") != 2:
        raise RuntimeError("Orthogonal phenotype count drift")
    if seal.get("required_min_backgrounds") != 3:
        raise RuntimeError("Minimum background count drift")
    if seal.get("required_min_differentiations_per_background_condition") != 2:
        raise RuntimeError("Minimum differentiation count drift")

    for name, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Execution-contract hash mismatch: {name}")

    for name, required in EXPECTED_SCHEMA_FIELDS.items():
        path = CANDIDATE_ROOT / name
        if read_tsv(path):
            raise RuntimeError(f"Outcome-blind schema contains rows: {name}")
        observed = set(header(path))
        if not required <= observed:
            raise RuntimeError(f"Execution schema lacks required fields: {name}")

    policies = read_tsv(CANDIDATE_ROOT / "stage_a_execution_policy.tsv")
    if [row["rule_id"] for row in policies] != [f"EX{i:02d}" for i in range(1, 13)]:
        raise RuntimeError("Stage-A execution policy universe drift")
    policy_text = "\n".join(row["frozen_rule"] for row in policies).lower()
    for required in ["three independent backgrounds", "two independent differentiations", "exactly two non-rna", "valid null"]:
        if required not in policy_text:
            raise RuntimeError(f"Execution policy omitted: {required}")

    assays = read_tsv(CANDIDATE_ROOT / "stage_a_required_assay_roles.tsv")
    roles = {row["assay_role"] for row in assays}
    if roles != {
        "CIS_TARGET", "RECIPIENT_RNA", "SECRETOME_DISCOVERY", "VIABILITY",
        "COMPOSITION", "PHENOTYPE",
    }:
        raise RuntimeError("Required assay-role universe drift")
    recipient = next(row for row in assays if row["assay_role"] == "RECIPIENT_RNA")
    if not {"early_cis", "late_relay"} <= set(recipient["required_time_roles"].split(";")):
        raise RuntimeError("Recipient RNA lacks time-ordering requirements")

    estimands = read_tsv(CANDIDATE_ROOT / "stage_a_analysis_estimand_contract.tsv")
    if [row["estimand_id"] for row in estimands] != [
        "CIS01", "RELAY01", "ORIGIN01", "ROUTE01", "COMP01", "PHENO01"
    ]:
        raise RuntimeError("Execution estimand hierarchy drift")
    if any(row["correction"] != "BH" for row in estimands):
        raise RuntimeError("Multiplicity correction is not explicitly frozen")

    manifest = read_tsv(CANDIDATE_ROOT / "stage_a_execution_contract_input_manifest.tsv")
    expected_roles = {
        "stage_a_execution_contract_producer", "stage_a_execution_contract_validator",
        "stage_a_target_freeze_producer", "stage_a_target_freeze_validator",
        "stage_a_execution_freeze_producer", "stage_a_execution_freeze_validator",
        "stage_a_template_seal",
    }
    if {row["role"] for row in manifest} != expected_roles:
        raise RuntimeError("Execution-contract input universe drift")
    for row in manifest:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Execution-contract source drift: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Execution-contract source hash drift: {path}")

    forbidden = {"ACSL5", "CENPQ", "ERCC2", "IL18R1", "NECAB2", "PNPLA6", "RANBP17", "SHMT1", "ZBTB41"}
    text = "\n".join(
        (CANDIDATE_ROOT / name).read_text(encoding="utf-8")
        for name in seal["output_sha256"]
    )
    if any(gene in text for gene in forbidden):
        raise RuntimeError("Outcome-blind execution contract contains a legacy target")
    print(
        "STAGE_A_EXECUTION_CONTRACT_VALIDATION_PASS "
        "rules=12 assays=6 estimands=6 targets_frozen=false outcomes_opened=false"
    )


if __name__ == "__main__":
    main()
