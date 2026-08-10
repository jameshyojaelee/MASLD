#!/usr/bin/env python3
"""Independently validate the Stage-A terminal/Stage-B freeze contract."""

from __future__ import annotations

import csv
import json

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def header(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.reader(handle, delimiter="\t"))


REQUIRED = {
    "stage_a_primary_estimand_results_template.tsv": {"estimand_id", "target_uid", "guide_id", "endpoint_id", "effect", "p_value", "q_value", "result_status"},
    "stage_a_sensitivity_results_template.tsv": {"sensitivity_id", "target_uid", "guide_id", "recipient_lineage", "effect", "result_status"},
    "stage_a_mediator_candidates_template.tsv": {"target_uid", "recipient_lineage", "mediator_id", "p_value", "q_value", "blocking_reagent", "upstream_cis_preservation_assay", "off_target_control"},
    "stage_a_outcome_signoff_template.tsv": {"unblinding_seal_path", "unblinding_seal_sha256", "primary_results_path", "primary_results_sha256", "complete_outcome_universe"},
    "stage_b_exact_design_template.tsv": {"target_uid", "recipient_lineage", "mediator_id", "variant_id", "snp_hg19", "orientation_uid", "exact_design_uid", "risk_allele", "editing_method", "pegrna_id_1", "pegrna_id_2", "rescue_method", "blocking_reagent", "phenotype_id_1", "phenotype_id_2"},
    "stage_b_power_units_template.tsv": {"target_uid", "background_id", "allele_role", "clone_id", "differentiation_id", "independent_clone", "independent_differentiation"},
    "stage_b_randomization_template.tsv": {"sample_id", "blinded_label", "target_uid", "background_id", "allele_role", "clone_id", "differentiation_id", "arm_id", "assay_id", "randomized_allocation", "unblinding_status"},
    "stage_b_freeze_signoff_template.tsv": {"stage_b_design_path", "stage_b_design_sha256", "power_units_path", "power_units_sha256", "randomization_path", "randomization_sha256", "stage_b_outcomes_inspected"},
}


def main() -> None:
    seal = json.loads((CANDIDATE_ROOT / "STAGE_A_TERMINAL_CONTRACT_SEALED.json").read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_outcome_blind_stage_a_terminal_and_stage_b_freeze_contract":
        raise RuntimeError("Invalid Stage-A terminal contract status")
    for field in [
        "scientific_outcomes_inspected", "stage_a_terminal_verdict_frozen",
        "stage_b_design_frozen", "stage_b_promotion_authorized",
    ]:
        if seal.get(field) is not False:
            raise RuntimeError(f"Outcome-blind terminal contract improperly sets {field}")
    expected_constants = {
        "maximum_stage_b_nominations": 2, "primary_q_threshold": 0.05,
        "mediator_q_threshold": 0.05, "minimum_stage_a_backgrounds": 2,
        "minimum_stage_b_backgrounds": 3,
        "minimum_clones_per_allele_background": 2,
        "minimum_differentiations_per_clone": 2, "n_policy_rules": 10,
        "n_input_schemas": 8,
    }
    for field, expected in expected_constants.items():
        if seal.get(field) != expected:
            raise RuntimeError(f"Terminal contract constant drift: {field}")
    for name, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Terminal contract hash mismatch: {name}")
    for name, required in REQUIRED.items():
        path = CANDIDATE_ROOT / name
        if read_tsv(path) or not required <= set(header(path)):
            raise RuntimeError(f"Terminal contract schema drift: {name}")
    policies = read_tsv(CANDIDATE_ROOT / "stage_a_terminal_policy.tsv")
    if [row["rule_id"] for row in policies] != [f"ADJ{i:02d}" for i in range(1, 11)]:
        raise RuntimeError("Terminal policy universe drift")
    text = "\n".join(row["frozen_rule"] for row in policies).lower()
    for phrase in ["both guides", "at most two loci", "fixed-cell-count", "replacement-locus", "valid null"]:
        if phrase not in text:
            raise RuntimeError(f"Terminal contract policy omission: {phrase}")
    manifest = read_tsv(CANDIDATE_ROOT / "stage_a_terminal_contract_input_manifest.tsv")
    if {row["role"] for row in manifest} != {
        "terminal_contract_producer", "terminal_contract_validator",
        "terminal_adjudication_producer", "terminal_adjudication_validator",
        "execution_contract_seal", "qc_unblinding_contract_seal",
    }:
        raise RuntimeError("Terminal contract source universe drift")
    for row in manifest:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Terminal contract source drift: {path}")
    print("STAGE_A_TERMINAL_CONTRACT_VALIDATION_PASS rules=10 schemas=8 outcomes_opened=false stage_b_frozen=false")


if __name__ == "__main__":
    main()
