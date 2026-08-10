#!/usr/bin/env python3
"""Independently validate metadata-only Plan 46A clinical outreach."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[4]
ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-clinical-outreach-v1-2026-08-10"
PARENT_ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-clinical-data-handoff-v1-2026-08-10"
PARENT_SEAL = PARENT_ROOT / "CLINICAL_DATA_HANDOFF_SEALED.json"
EXPECTED = {"nash_crn_as116": 8, "maestro_nash": 6, "essence": 6}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    seal_path = ROOT / "CLINICAL_OUTREACH_SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_metadata_only_clinical_outreach":
        raise RuntimeError("Invalid clinical outreach status")
    for field, expected in {
        "n_routes": 3,
        "required_independent_human_cohorts": 2,
        "n_candidate_cohort_A_routes": 1,
        "n_candidate_cohort_B_routes": 2,
        "n_route_questions": 20,
        "bundle_read_only": True,
        "message_sent": False,
        "participant_data_accessed": False,
        "clinical_outcomes_accessed": False,
        "molecular_outcomes_accessed": False,
        "permission_granted": False,
        "analysis_authorized": False,
        "paper_promotion_authorized": False,
    }.items():
        if seal.get(field) != expected:
            raise RuntimeError(f"Clinical outreach seal drift: {field}")
    if seal.get("parent_handoff_seal_sha256") != sha256_file(PARENT_SEAL):
        raise RuntimeError("Clinical parent handoff drift")
    if (ROOT / "CLINICAL_OUTREACH_SEAL_SHA256.txt").read_text(encoding="utf-8").strip() != sha256_file(seal_path):
        raise RuntimeError("Clinical outreach seal sidecar mismatch")
    for relative, expected in seal["output_sha256"].items():
        path = ROOT / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Clinical outreach payload drift: {relative}")
    manifest = read_tsv(ROOT / "file_manifest.tsv")
    for row in manifest:
        path = ROOT / row["relative_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Clinical file-manifest drift: {path}")
    source = read_tsv(ROOT / "source_manifest.tsv")
    if {row["role"] for row in source} != {"plan46a", "parent_clinical_handoff_seal", "producer", "validator"}:
        raise RuntimeError("Clinical source-manifest role drift")
    for row in source:
        path = PROJECT_ROOT / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Clinical source drift: {path}")

    routes = read_tsv(ROOT / "route_registry.tsv")
    if {row["route_id"] for row in routes} != set(EXPECTED):
        raise RuntimeError("Clinical route universe drift")
    if [sum(row["candidate_cohort_role"] == role for row in routes) for role in ("candidate_cohort_A", "candidate_cohort_B")] != [1, 2]:
        raise RuntimeError("Clinical candidate-cohort role drift")
    if len({row["independence_group"] for row in routes}) != 3:
        raise RuntimeError("Clinical routes are not institutionally distinct")
    for row in routes:
        for field in ("message_sent", "participant_data_received", "outcome_received", "permission_granted", "response_received"):
            if row[field] != "false":
                raise RuntimeError(f"Clinical route prematurely advanced: {row['route_id']} {field}")

    parent_files = {
        name: (PARENT_ROOT / name).read_bytes()
        for name in ("clinical_cohort_feasibility_questions.tsv", "assay_availability_matrix.tsv", "clinical_role_firewall.tsv", "minimum_data_dictionary.tsv")
    }
    for route_id, n_questions in EXPECTED.items():
        route_root = ROOT / "routes" / route_id
        for name, expected in parent_files.items():
            if (route_root / name).read_bytes() != expected:
                raise RuntimeError(f"Clinical parent table drift in {route_id}: {name}")
        questions = read_tsv(route_root / "clinical_cohort_feasibility_questions.tsv")
        if [row["question_id"] for row in questions] != [f"CF{i:02d}" for i in range(1, 17)]:
            raise RuntimeError(f"Clinical frozen question drift in {route_id}")
        if any(row["cohort_role_candidate"] or row["answer"] or row["evidence_path"] or row["owner"] or row["signed_utc"] for row in questions):
            raise RuntimeError(f"Clinical frozen response prefilled in {route_id}")
        assays = read_tsv(route_root / "assay_availability_matrix.tsv")
        if any(row["available"] or row["n_baseline"] or row["n_followup"] or row["n_authoritative_pairs"] or row["evidence_path"] for row in assays):
            raise RuntimeError(f"Clinical assay response prefilled in {route_id}")
        addendum = read_tsv(route_root / "route_specific_questions.tsv")
        if len(addendum) != n_questions or any(row["answer"] or row["evidence_path"] or row["owner"] or row["signed_utc"] for row in addendum):
            raise RuntimeError(f"Clinical route addendum drift in {route_id}")
        inventory = read_tsv(route_root / "public_inventory_audit.tsv")
        if {row["evidence_state"] for row in inventory} != {"public_aggregate", "unverified_route_gate"}:
            raise RuntimeError(f"Clinical public/unverified boundary drift in {route_id}")
        if any(row["custodian_confirmation"] for row in inventory):
            raise RuntimeError(f"Clinical custodian confirmation prefilled in {route_id}")
        email = " ".join((route_root / "EMAIL_DRAFT.md").read_text(encoding="utf-8").lower().split())
        if "before requesting" not in email and "at this stage we request no participant-level data and no outcome" not in email:
            raise RuntimeError(f"Clinical outcome-free boundary omitted in {route_id}")
        response = " ".join((route_root / "RESPONSE_INSTRUCTIONS.md").read_text(encoding="utf-8").lower().split())
        for phrase in ("metadata-only", "do not send participant-level", "does not grant permission", "14 calendar days"):
            if phrase not in response:
                raise RuntimeError(f"Clinical response firewall omitted in {route_id}: {phrase}")

    send_log = read_tsv(ROOT / "send_log_template.tsv")
    if {row["route_id"] for row in send_log} != set(EXPECTED) or any(row["status"] != "not_sent" for row in send_log):
        raise RuntimeError("Clinical send log drift")
    firewall = read_tsv(ROOT / "outcome_firewall.tsv")
    if len(firewall) != 1 or firewall[0] != {
        "public_aggregate_inventory_only": "true",
        "participant_rows": "0",
        "direct_identifiers": "0",
        "responder_linked_molecular_outcomes": "0",
        "permissions_assumed": "0",
        "messages_sent": "0",
        "verdict": "pass",
    }:
        raise RuntimeError("Clinical outcome firewall drift")
    for path in ROOT.rglob("*"):
        if path.is_file() and (path.stat().st_mode & 0o222):
            raise RuntimeError(f"Clinical outreach file remains writable: {path}")
    if ROOT.stat().st_mode & 0o222:
        raise RuntimeError("Clinical outreach root remains writable")
    print(
        "CLINICAL_OUTREACH_VALIDATION_PASS "
        f"routes=3 cohort_A=1 cohort_B=2 route_questions=20 participant_data=false outcomes=false seal_sha256={sha256_file(seal_path)}"
    )


if __name__ == "__main__":
    main()
