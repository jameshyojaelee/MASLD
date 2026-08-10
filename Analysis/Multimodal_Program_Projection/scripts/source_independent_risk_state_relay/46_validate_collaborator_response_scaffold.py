#!/usr/bin/env python3
"""Validate the initially blank, editable collaborator-response scaffold."""

from __future__ import annotations

import json

from relay_common import CANDIDATE_ROOT, SCRIPT_ROOT, read_tsv, sha256_file


def main() -> None:
    marker_path = CANDIDATE_ROOT / "DRAFT_EDITABLE.json"
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    if marker.get("status") != "editable_collaborator_response_scaffold":
        raise RuntimeError("Invalid collaborator-response scaffold status")
    if marker.get("mutable_response_workspace") is not True:
        raise RuntimeError("Response scaffold is not explicitly editable")
    for field in [
        "target_identity_present", "scientific_outcomes_present",
        "target_disclosure_permitted",
    ]:
        if marker.get(field) is not False:
            raise RuntimeError(f"Response scaffold improperly sets {field}")
    if marker.get("parent_bundle_read_only") is not True:
        raise RuntimeError("Response scaffold does not preserve the parent handoff")
    if marker.get("producer_sha256") != sha256_file(
        SCRIPT_ROOT / "45_build_collaborator_response_scaffold.py"
    ):
        raise RuntimeError("Response-scaffold producer drift")
    for relative, expected in marker["initial_template_sha256"].items():
        path = CANDIDATE_ROOT / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Initial response scaffold drift: {relative}")

    questions = read_tsv(CANDIDATE_ROOT / "platform_feasibility_questions.tsv")
    if [row["question_id"] for row in questions] != [f"FQ{i:02d}" for i in range(1, 13)]:
        raise RuntimeError("Response-scaffold feasibility universe drift")
    if any(row["answer"] or row["evidence_path"] or row["owner"] or row["signed_utc"] for row in questions):
        raise RuntimeError("Response-scaffold feasibility fields were prefilled")
    roles = read_tsv(CANDIDATE_ROOT / "role_and_blinding_firewall.tsv")
    if len(roles) != 7 or any(row["assigned_person"] or row["accepted_utc"] for row in roles):
        raise RuntimeError("Response-scaffold role fields were prefilled")
    deliverables = read_tsv(CANDIDATE_ROOT / "collaborator_deliverable_register.tsv")
    if [row["deliverable_id"] for row in deliverables] != [f"D{i:02d}" for i in range(1, 15)]:
        raise RuntimeError("Response-scaffold deliverable universe drift")
    if any(row["status"] != "not_started" for row in deliverables):
        raise RuntimeError("Response-scaffold deliverables were prefilled")
    identity = read_tsv(CANDIDATE_ROOT / "response_identity.tsv")
    if len(identity) != 1 or identity[0]["response_status"] != "editable_draft":
        raise RuntimeError("Response-scaffold identity drift")
    if any(identity[0][field] for field in ["organization", "prepared_by", "created_utc"]):
        raise RuntimeError("Response-scaffold identity was prefilled")
    if read_tsv(CANDIDATE_ROOT / "evidence_file_register.tsv"):
        raise RuntimeError("Response-scaffold evidence register was prefilled")
    print(
        "COLLABORATOR_RESPONSE_SCAFFOLD_VALIDATION_PASS "
        "questions=12 roles=7 deliverables=14 mutable=true target_present=false"
    )


if __name__ == "__main__":
    main()
