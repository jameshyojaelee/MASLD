#!/usr/bin/env python3
"""Independently validate the Stage-A blinded-QC/unblinding contract."""

from __future__ import annotations

import csv
import json

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def header(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.reader(handle, delimiter="\t"))


REQUIRED_FIELDS = {
    "blinded_sample_file_manifest_template.tsv": {
        "sample_id", "blinded_label", "assay_id", "raw_source_path", "size_bytes",
        "sha256", "file_readable", "scientific_condition_inspected",
    },
    "blinded_assay_qc_thresholds_template.tsv": {
        "assay_id", "metric_name", "operator", "threshold", "frozen_utc",
        "raw_qc_accessed_before_freeze", "reviewer_1", "reviewer_2",
        "review_concordant",
    },
    "threshold_freeze_signoff_template.tsv": {
        "thresholds_path", "thresholds_sha256", "reviewer_1", "reviewer_2",
        "review_concordant", "raw_qc_accessed_before_freeze",
        "scientific_condition_inspected", "signed_utc",
    },
    "blinded_sample_qc_metrics_template.tsv": {
        "sample_id", "blinded_label", "assay_id", "metric_name", "metric_value",
        "scientific_condition_inspected",
    },
    "target_independent_guide_response_template.tsv": {
        "gene_id", "contrast", "effect", "p_value", "q_value", "assay_id",
        "n_backgrounds", "exclusion_gene", "condition_labels_opened",
    },
    "blinded_qc_signoff_template.tsv": {
        "review_scope", "reviewer_1", "reviewer_2", "review_concordant",
        "scientific_condition_inspected", "signed_utc",
    },
    "unblinding_key_template.tsv": {
        "sample_id", "blinded_label", "biological_unit_id", "background_id",
        "differentiation_id", "target_uid", "guide_id", "arm_id", "time_role",
        "challenge_role", "assay_id", "key_generated_before_qc",
        "key_opened_after_qc_seal",
    },
    "unblinding_signoff_template.tsv": {
        "qc_seal_path", "qc_seal_sha256", "unblinding_key_path",
        "unblinding_key_sha256", "data_manager", "analysis_lead",
        "key_opened_after_qc_seal", "scientific_outcomes_inspected",
    },
}


def main() -> None:
    seal_path = CANDIDATE_ROOT / "STAGE_A_QC_UNBLINDING_CONTRACT_SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_outcome_blind_stage_a_qc_unblinding_contract":
        raise RuntimeError("Invalid QC/unblinding contract status")
    for field in [
        "scientific_outcomes_inspected", "scientific_condition_labels_opened",
        "sample_exclusions_frozen", "unblinding_authorized", "stage_b_design_frozen",
    ]:
        if seal.get(field) is not False:
            raise RuntimeError(f"Outcome-blind QC contract improperly sets {field}")
    if seal.get("governing_plan_integrity_bound") is not False:
        raise RuntimeError("Mutable plan prose must not be integrity-bound")
    if seal.get("guide_response_q_threshold") != 0.05 or seal.get("guide_response_abs_effect_threshold") != 0.25:
        raise RuntimeError("Guide-response exclusion threshold drift")
    if seal.get("minimum_backgrounds_per_claim_bearing_comparison") != 2:
        raise RuntimeError("Post-QC comparison background threshold drift")
    if seal.get("n_policy_rules") != 10 or seal.get("n_input_schemas") != 8:
        raise RuntimeError("QC/unblinding contract family-size drift")
    for name, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"QC/unblinding contract hash mismatch: {name}")
    for name, required in REQUIRED_FIELDS.items():
        path = CANDIDATE_ROOT / name
        if read_tsv(path):
            raise RuntimeError(f"QC/unblinding schema contains rows: {name}")
        if not required <= set(header(path)):
            raise RuntimeError(f"QC/unblinding schema lacks required fields: {name}")

    policies = read_tsv(CANDIDATE_ROOT / "stage_a_qc_unblinding_policy.tsv")
    if [row["rule_id"] for row in policies] != [f"QC{i:02d}" for i in range(1, 11)]:
        raise RuntimeError("QC/unblinding policy universe drift")
    text = "\n".join(row["frozen_rule"] for row in policies).lower()
    for phrase in [
        "opaque labels", "after the blinded-qc seal", "cannot add or replace a sample",
        "at least two backgrounds", "valid null",
    ]:
        if phrase not in text:
            raise RuntimeError(f"QC/unblinding policy omitted: {phrase}")

    manifest = read_tsv(CANDIDATE_ROOT / "stage_a_qc_unblinding_contract_input_manifest.tsv")
    expected_roles = {
        "qc_contract_producer", "qc_contract_validator", "qc_threshold_producer",
        "qc_threshold_validator", "blinded_qc_producer", "blinded_qc_validator",
        "unblinding_producer", "unblinding_validator", "execution_contract_seal",
    }
    if {row["role"] for row in manifest} != expected_roles:
        raise RuntimeError("QC/unblinding contract source universe drift")
    for row in manifest:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"QC/unblinding contract source drift: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"QC/unblinding contract source hash drift: {path}")
    print(
        "STAGE_A_QC_UNBLINDING_CONTRACT_VALIDATION_PASS "
        "rules=10 schemas=8 conditions_opened=false outcomes_opened=false"
    )


if __name__ == "__main__":
    main()
