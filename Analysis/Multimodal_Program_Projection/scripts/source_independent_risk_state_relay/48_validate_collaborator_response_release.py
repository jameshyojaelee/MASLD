#!/usr/bin/env python3
"""Independently rederive a sealed collaborator-feasibility adjudication."""

from __future__ import annotations

import json

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def truth(value: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise RuntimeError(f"Invalid serialized boolean: {value}")


def expected_verdict(answers: dict[str, str]) -> tuple[str, bool, bool]:
    stage_a = [answers[f"FQ{i:02d}"] for i in range(1, 11)]
    stage_b = [answers[f"FQ{i:02d}"] for i in range(11, 13)]
    stage_a_ok = all(value == "yes" for value in stage_a)
    full_ok = stage_a_ok and all(value == "yes" for value in stage_b)
    if "no" in stage_a:
        verdict = "complete_stage_a_infeasible"
    elif "unresolved" in stage_a:
        verdict = "pending_stage_a_feasibility"
    elif "no" in stage_b:
        verdict = "complete_stage_a_only_not_full_mechanism"
    elif "unresolved" in stage_b:
        verdict = "pending_stage_b_feasibility"
    else:
        verdict = "complete_full_mechanism_feasible"
    return verdict, stage_a_ok, full_ok


def main() -> None:
    seal = json.loads((CANDIDATE_ROOT / "COLLABORATOR_RESPONSE_SEALED.json").read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_collaborator_response_adjudication":
        raise RuntimeError("Invalid collaborator-response release status")
    for field in [
        "target_identity_present", "target_frozen", "scientific_outcomes_present",
        "stage_a_outcome_generation_permitted", "stage_b_frozen",
        "paper_promotion_authorized",
    ]:
        if seal.get(field) is not False:
            raise RuntimeError(f"Collaborator response improperly sets {field}")
    for relative, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Collaborator-response output drift: {relative}")

    source_rows = read_tsv(CANDIDATE_ROOT / "collaborator_response_source_manifest.tsv")
    expected_roles = {"response_producer", "response_validator", "handoff_seal", "response_draft_marker"}
    if {row["role"] for row in source_rows} != expected_roles:
        raise RuntimeError("Collaborator-response source universe drift")
    source_by_role = {row["role"]: row for row in source_rows}
    for row in source_rows:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Collaborator-response source drift: {path}")
    handoff_seal_path = PROJECT_ROOT / source_by_role["handoff_seal"]["source_path"]
    handoff_root = handoff_seal_path.parent
    if handoff_root.name != seal.get("parent_handoff_candidate_id") or sha256_file(handoff_seal_path) != seal.get("parent_handoff_seal_sha256"):
        raise RuntimeError("Collaborator-response parent binding drift")
    handoff_seal = json.loads(handoff_seal_path.read_text(encoding="utf-8"))
    if handoff_seal.get("bundle_read_only") is not True:
        raise RuntimeError("Collaborator-response parent is not read-only")

    questions = read_tsv(CANDIDATE_ROOT / "response/platform_feasibility_questions.tsv")
    frozen_questions = read_tsv(handoff_root / "platform_feasibility_questions.tsv")
    static_q = ["question_id", "domain", "frozen_question", "impact"]
    if len(questions) != 12 or any(
        any(left[field] != right[field] for field in static_q)
        for left, right in zip(questions, frozen_questions)
    ):
        raise RuntimeError("Collaborator-response question drift")
    answers = {row["question_id"]: row["answer"].strip().lower() for row in questions}
    if set(answers) != {f"FQ{i:02d}" for i in range(1, 13)} or not set(answers.values()) <= {"yes", "no", "unresolved"}:
        raise RuntimeError("Collaborator-response answer universe drift")

    roles = read_tsv(CANDIDATE_ROOT / "response/role_and_blinding_firewall.tsv")
    frozen_roles = read_tsv(handoff_root / "role_and_blinding_firewall.tsv")
    static_r = ["role", "responsibility", "target_access", "outcome_access", "prohibition"]
    if len(roles) != 7 or any(
        any(left[field] != right[field] for field in static_r)
        for left, right in zip(roles, frozen_roles)
    ):
        raise RuntimeError("Collaborator-response role drift")
    assignments = {row["role"]: row["assigned_person"].strip().casefold() for row in roles}
    critical = [assignments[role] for role in ["data_manager", "qc_reviewer_1", "qc_reviewer_2", "analysis_lead"]]
    role_gate = all(critical) and len(set(critical)) == 4 and assignments["independent_validator"] != assignments["analysis_lead"]

    deliverables = read_tsv(CANDIDATE_ROOT / "response/collaborator_deliverable_register.tsv")
    d01_complete = {row["deliverable_id"]: row["status"] for row in deliverables}.get("D01") == "complete"
    verdict, stage_a_ok, full_ok = expected_verdict(answers)
    adjudication_rows = read_tsv(CANDIDATE_ROOT / "collaborator_response_adjudication.tsv")
    if len(adjudication_rows) != 1:
        raise RuntimeError("Collaborator-response adjudication cardinality drift")
    adjudication = adjudication_rows[0]
    expected = {
        "feasibility_verdict": verdict,
        "stage_a_platform_feasible": stage_a_ok,
        "full_mechanism_platform_feasible": full_ok,
        "role_firewall_accepted": bool(role_gate),
        "d01_complete": d01_complete,
        "stage_a_planning_permitted": stage_a_ok and bool(role_gate) and d01_complete,
        "target_disclosure_permitted": full_ok and bool(role_gate) and d01_complete,
    }
    for field, value in expected.items():
        observed = adjudication[field] if isinstance(value, str) else truth(adjudication[field])
        if observed != value or seal.get(field) != value:
            raise RuntimeError(f"Collaborator-response verdict drift: {field}")
    counts = {
        "n_yes": sum(value == "yes" for value in answers.values()),
        "n_no": sum(value == "no" for value in answers.values()),
        "n_unresolved": sum(value == "unresolved" for value in answers.values()),
    }
    if any(int(adjudication[field]) != value for field, value in counts.items()):
        raise RuntimeError("Collaborator-response answer-count drift")

    evidence = read_tsv(CANDIDATE_ROOT / "response/evidence_file_register.tsv")
    registered = set()
    for row in evidence:
        path = CANDIDATE_ROOT / "response" / row["relative_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Collaborator-response evidence drift: {path}")
        registered.add(row["relative_path"])
    referenced = {
        part.strip()
        for row in questions
        for part in row["evidence_path"].split(";")
        if part.strip()
    }
    if registered != referenced or len(registered) != int(adjudication["n_evidence_files"]) or len(registered) != seal.get("n_evidence_files"):
        raise RuntimeError("Collaborator-response evidence-universe drift")
    print(
        "COLLABORATOR_RESPONSE_RELEASE_VALIDATION_PASS "
        f"verdict={verdict} target_disclosure={str(expected['target_disclosure_permitted']).lower()} "
        "target_frozen=false outcomes_present=false"
    )


if __name__ == "__main__":
    main()
