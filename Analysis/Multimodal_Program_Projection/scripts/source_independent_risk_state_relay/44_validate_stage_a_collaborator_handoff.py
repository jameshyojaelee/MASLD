#!/usr/bin/env python3
"""Independently validate the target-independent collaborator handoff."""

from __future__ import annotations

import json

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def main() -> None:
    seal = json.loads((CANDIDATE_ROOT / "STAGE_A_COLLABORATOR_HANDOFF_SEALED.json").read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_target_independent_stage_a_collaborator_handoff":
        raise RuntimeError("Invalid collaborator-handoff status")
    for field in [
        "experimental_targets_frozen", "scientific_conditions_opened",
        "scientific_outcomes_inspected", "wetlab_samples_exist",
        "stage_b_design_frozen", "paper_promotion_authorized",
    ]:
        if seal.get(field) is not False:
            raise RuntimeError(f"Collaborator handoff improperly sets {field}")
    if seal.get("bundle_read_only") is not True:
        raise RuntimeError("Collaborator handoff does not preserve immutable response separation")
    constants = {
        "n_input_contracts": 5, "n_feasibility_questions": 12,
        "n_deliverables": 14, "n_roles": 7,
    }
    for field, expected in constants.items():
        if seal.get(field) != expected:
            raise RuntimeError(f"Collaborator handoff count drift: {field}")
    if seal.get("primary_recipient_lineages") != [
        "Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Macrophages"
    ]:
        raise RuntimeError("Collaborator recipient-axis universe drift")
    for relative, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Collaborator output hash mismatch: {relative}")

    questions = read_tsv(CANDIDATE_ROOT / "platform_feasibility_questions.tsv")
    if [row["question_id"] for row in questions] != [f"FQ{i:02d}" for i in range(1, 13)]:
        raise RuntimeError("Collaborator feasibility-question universe drift")
    if any(row["answer"] or row["evidence_path"] or row["owner"] or row["signed_utc"] for row in questions):
        raise RuntimeError("Target-independent feasibility answers were prefilled")
    if {row["impact"] for row in questions} != {"blocking", "stage_b_blocking"}:
        raise RuntimeError("Collaborator feasibility impact classes drift")
    deliverables = read_tsv(CANDIDATE_ROOT / "collaborator_deliverable_register.tsv")
    if [row["deliverable_id"] for row in deliverables] != [f"D{i:02d}" for i in range(1, 15)]:
        raise RuntimeError("Collaborator deliverable universe drift")
    if any(row["status"] != "not_started" for row in deliverables):
        raise RuntimeError("Collaborator deliverables were prematurely started")
    roles = read_tsv(CANDIDATE_ROOT / "role_and_blinding_firewall.tsv")
    required_roles = {
        "platform_lead", "wetlab_lead", "data_manager", "qc_reviewer_1",
        "qc_reviewer_2", "analysis_lead", "independent_validator",
    }
    if {row["role"] for row in roles} != required_roles or any(row["assigned_person"] or row["accepted_utc"] for row in roles):
        raise RuntimeError("Collaborator role universe/assignment drift")

    bundle = read_tsv(CANDIDATE_ROOT / "contract_bundle_manifest.tsv")
    if {row["contract_role"] for row in bundle} != {
        "stage_a_template", "target_adjudication", "execution",
        "qc_unblinding", "terminal",
    }:
        raise RuntimeError("Collaborator contract-bundle universe drift")
    for row in bundle:
        source = PROJECT_ROOT / row["source_path"]
        copy = CANDIDATE_ROOT / row["bundle_path"]
        if (
            not source.is_file() or not copy.is_file()
            or source.stat().st_size != int(row["size_bytes"])
            or sha256_file(source) != row["sha256"]
            or sha256_file(copy) != row["sha256"]
        ):
            raise RuntimeError(f"Collaborator contract copy/source drift: {copy}")
    manifest = read_tsv(CANDIDATE_ROOT / "handoff_source_manifest.tsv")
    if {row["role"] for row in manifest} != {
        "handoff_producer", "handoff_validator", "stage_a_template_seal",
        "target_adjudication_seal", "execution_seal", "qc_unblinding_seal",
        "terminal_seal",
    }:
        raise RuntimeError("Collaborator source-manifest universe drift")
    for row in manifest:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Collaborator source drift: {path}")
    handoff = (CANDIDATE_ROOT / "COLLABORATOR_HANDOFF.md").read_text(encoding="utf-8").lower()
    for phrase in [
        "contains no target", "at least three backgrounds", "two independent differentiations",
        "cells, wells, organoids, and lanes never become biological replicates",
        "exact endogenous", "failed primary gate",
    ]:
        if phrase not in handoff:
            raise RuntimeError(f"Collaborator handoff omitted: {phrase}")
    response = (CANDIDATE_ROOT / "RESPONSE_INSTRUCTIONS.md").read_text(encoding="utf-8").lower()
    for phrase in [
        "immutable and read-only", "do not edit", "outside this candidate root",
        "parent candidate id", "sha256", "never write the response back",
    ]:
        if phrase not in response:
            raise RuntimeError(f"Collaborator response separation omitted: {phrase}")
    print(
        "STAGE_A_COLLABORATOR_HANDOFF_VALIDATION_PASS "
        "contracts=5 questions=12 deliverables=14 targets_frozen=false outcomes_opened=false"
    )


if __name__ == "__main__":
    main()
