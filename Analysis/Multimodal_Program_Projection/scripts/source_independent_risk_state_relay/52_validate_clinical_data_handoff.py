#!/usr/bin/env python3
"""Independently validate the outcome-free clinical collaboration packet."""

from __future__ import annotations

import json

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def main() -> None:
    seal = json.loads((CANDIDATE_ROOT / "CLINICAL_DATA_HANDOFF_SEALED.json").read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_outcome_free_clinical_data_collaborator_handoff":
        raise RuntimeError("Invalid clinical-data handoff status")
    for field in [
        "participant_data_accessed", "clinical_outcomes_accessed",
        "molecular_outcomes_accessed", "target_identity_present",
        "permission_granted", "paper_promotion_authorized",
    ]:
        if seal.get(field) is not False:
            raise RuntimeError(f"Clinical handoff improperly sets {field}")
    if seal.get("required_independent_human_cohorts") != 2:
        raise RuntimeError("Clinical handoff weakens the two-cohort requirement")
    for relative, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Clinical handoff output drift: {relative}")

    source_rows = read_tsv(CANDIDATE_ROOT / "clinical_handoff_source_manifest.tsv")
    if {row["role"] for row in source_rows} != {"goal_audit", "producer", "validator"}:
        raise RuntimeError("Clinical handoff source universe drift")
    for row in source_rows:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Clinical handoff source drift: {path}")
    audit_path = PROJECT_ROOT / next(row["source_path"] for row in source_rows if row["role"] == "goal_audit")
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if audit.get("status") != "goal_not_achieved_new_data_required" or audit.get("full_risk_to_state_architecture_achieved") is not False:
        raise RuntimeError("Clinical handoff is not bound to the incomplete full-goal audit")

    questions = read_tsv(CANDIDATE_ROOT / "clinical_cohort_feasibility_questions.tsv")
    if [row["question_id"] for row in questions] != [f"CF{i:02d}" for i in range(1, 17)]:
        raise RuntimeError("Clinical feasibility universe drift")
    if any(
        row["cohort_role_candidate"] or row["answer"] or row["evidence_path"]
        or row["owner"] or row["signed_utc"]
        for row in questions
    ):
        raise RuntimeError("Clinical feasibility fields were prefilled")
    impacts = {row["impact"] for row in questions}
    if impacts != {"cohort_blocking", "replication_blocking", "orthogonal_blocking", "governance_blocking"}:
        raise RuntimeError("Clinical feasibility impact vocabulary drift")

    dictionary = read_tsv(CANDIDATE_ROOT / "minimum_data_dictionary.tsv")
    required_fields = {
        "participant_id", "cohort_id", "sample_id", "timepoint",
        "intervention", "nas", "steatohepatitis_status", "fibrosis_stage",
        "responder_definition_source", "rna_assay_id", "processing_qc",
    }
    if not required_fields <= {row["field"] for row in dictionary}:
        raise RuntimeError("Clinical minimum data dictionary is incomplete")
    if any(row["field"] in {"name", "date_of_birth", "medical_record_number"} for row in dictionary):
        raise RuntimeError("Clinical handoff requests direct identifiers")
    assays = read_tsv(CANDIDATE_ROOT / "assay_availability_matrix.tsv")
    if {row["assay_id"] for row in assays} != {
        "bulk_rna", "histology", "tissue_proteomics", "secreted_proteomics",
        "spatial_transcriptomics_or_proteomics",
    }:
        raise RuntimeError("Clinical assay universe drift")
    if any(row["available"] or row["n_baseline"] or row["n_followup"] or row["n_authoritative_pairs"] or row["evidence_path"] for row in assays):
        raise RuntimeError("Clinical assay availability was prefilled")
    roles = read_tsv(CANDIDATE_ROOT / "clinical_role_firewall.tsv")
    if {row["role"] for row in roles} != {
        "clinical_custodian", "histopathology_lead", "data_manager",
        "assay_lead", "analysis_lead", "independent_validator",
    }:
        raise RuntimeError("Clinical role universe drift")
    if any(row["assigned_person"] or row["accepted_utc"] for row in roles):
        raise RuntimeError("Clinical roles were preassigned")
    brief = (CANDIDATE_ROOT / "CLINICAL_DATA_COLLABORATOR_BRIEF.md").read_text(encoding="utf-8").lower()
    for phrase in [
        "two genuinely independent", "metadata and availability only",
        "both histologic improvers and non-improvers", "participants, never biopsies",
        "no protected health information",
    ]:
        if phrase not in brief:
            raise RuntimeError(f"Clinical brief omitted: {phrase}")
    if seal.get("n_feasibility_questions") != len(questions) or seal.get("n_data_dictionary_fields") != len(dictionary) or seal.get("n_assay_roles") != len(assays) or seal.get("n_firewall_roles") != len(roles):
        raise RuntimeError("Clinical handoff count drift")
    print(
        "CLINICAL_DATA_HANDOFF_VALIDATION_PASS questions=16 cohorts_required=2 "
        "participant_data=false outcomes=false promotion=false"
    )


if __name__ == "__main__":
    main()
