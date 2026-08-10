#!/usr/bin/env python3
"""Independently validate the immutable Plan 45A outreach packets."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[4]
SCRIPT_ROOT = Path(__file__).resolve().parent
ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-platform-outreach-v3-2026-08-10"
PARENT_ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-stage-a-collaborator-handoff-v3-2026-08-10"
PARENT_SEAL = PARENT_ROOT / "STAGE_A_COLLABORATOR_HANDOFF_SEALED.json"
EXPECTED_ROUTES = {
    "du_tsinghua": 7,
    "takebe_cincinnati": 6,
    "ebrahimkhani_pittsburgh": 6,
    "ge_iorgantech": 6,
}


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
    seal_path = ROOT / "OUTREACH_PACKET_SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_target_blind_platform_outreach":
        raise RuntimeError("Invalid outreach status")
    for field, expected in {
        "n_routes": 4,
        "n_primary_parallel_routes": 3,
        "n_supporting_routes": 1,
        "bundle_read_only": True,
        "message_sent": False,
        "target_disclosed": False,
        "guide_disclosed": False,
        "condition_allocation_disclosed": False,
        "scientific_outcome_disclosed": False,
        "target_freeze_authorized": False,
        "stage_a_authorized": False,
        "stage_b_authorized": False,
        "paper_promotion_authorized": False,
    }.items():
        if seal.get(field) != expected:
            raise RuntimeError(f"Outreach seal drift: {field}")
    if seal.get("parent_handoff_seal_sha256") != sha256_file(PARENT_SEAL):
        raise RuntimeError("Parent handoff identity drift")
    if (ROOT / "PACKET_SEAL_SHA256.txt").read_text(encoding="utf-8").strip() != sha256_file(seal_path):
        raise RuntimeError("Outreach seal sidecar mismatch")
    for relative, expected in seal["output_sha256"].items():
        path = ROOT / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Sealed outreach output drift: {relative}")

    manifest = read_tsv(ROOT / "file_manifest.tsv")
    if len({row["relative_path"] for row in manifest}) != len(manifest):
        raise RuntimeError("Duplicate file-manifest path")
    for row in manifest:
        path = ROOT / row["relative_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"File-manifest mismatch: {path}")

    source_manifest = read_tsv(ROOT / "source_manifest.tsv")
    if {row["role"] for row in source_manifest} != {"plan45a", "parent_handoff_seal", "producer", "validator"}:
        raise RuntimeError("Source-manifest role drift")
    for row in source_manifest:
        path = PROJECT_ROOT / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Source drift: {path}")

    routes = read_tsv(ROOT / "route_registry.tsv")
    if {row["route_id"] for row in routes} != set(EXPECTED_ROUTES):
        raise RuntimeError("Route universe drift")
    if sum(row["send_priority"] == "primary_parallel" for row in routes) != 3:
        raise RuntimeError("Primary send universe drift")
    if any(row["message_sent"] != "false" or row["target_disclosed"] != "false" or row["response_received"] != "false" for row in routes):
        raise RuntimeError("Outreach route was prematurely advanced")
    if any(row["contact_status"] == "requires_current_confirmation" and row["contact_value"] for row in routes):
        raise RuntimeError("Unverified contact was prefilled")

    parent_questions = (PARENT_ROOT / "platform_feasibility_questions.tsv").read_bytes()
    parent_roles = (PARENT_ROOT / "role_and_blinding_firewall.tsv").read_bytes()
    parent_deliverables = (PARENT_ROOT / "collaborator_deliverable_register.tsv").read_bytes()
    for route_id, expected_questions in EXPECTED_ROUTES.items():
        route_root = ROOT / "routes" / route_id
        for name, expected in {
            "platform_feasibility_questions.tsv": parent_questions,
            "role_and_blinding_firewall.tsv": parent_roles,
            "collaborator_deliverable_register.tsv": parent_deliverables,
        }.items():
            if (route_root / name).read_bytes() != expected:
                raise RuntimeError(f"Frozen parent table drift in {route_id}: {name}")
        frozen = read_tsv(route_root / "platform_feasibility_questions.tsv")
        if [row["question_id"] for row in frozen] != [f"FQ{i:02d}" for i in range(1, 13)]:
            raise RuntimeError(f"FQ universe drift in {route_id}")
        if any(row["answer"] or row["evidence_path"] or row["owner"] or row["signed_utc"] for row in frozen):
            raise RuntimeError(f"Frozen FQ response prefilled in {route_id}")
        addendum = read_tsv(route_root / "route_specific_questions.tsv")
        if len(addendum) != expected_questions or len({row["route_question_id"] for row in addendum}) != expected_questions:
            raise RuntimeError(f"Route-question universe drift in {route_id}")
        if any(row["answer"] or row["evidence_path"] or row["owner"] or row["signed_utc"] for row in addendum):
            raise RuntimeError(f"Route-specific answer prefilled in {route_id}")
        audit = read_tsv(route_root / "public_capability_audit.tsv")
        states = {row["public_evidence_state"] for row in audit}
        if states != {"publicly_supported", "not_publicly_demonstrated"}:
            raise RuntimeError(f"Capability evidence-state boundary drift in {route_id}")
        email = (route_root / "EMAIL_DRAFT.md").read_text(encoding="utf-8").lower()
        for phrase in ("before", "target"):
            if phrase not in email:
                raise RuntimeError(f"Target-blind email language omitted in {route_id}: {phrase}")
        target_blind_boundaries = (
            "disclose no target",
            "without target disclosure",
            "no biological outcome or target",
            "before any target is disclosed",
            "target-independent",
        )
        if not any(phrase in email for phrase in target_blind_boundaries):
            raise RuntimeError(f"Target-disclosure boundary omitted in {route_id}")
        response = (route_root / "RESPONSE_INSTRUCTIONS.md").read_text(encoding="utf-8").lower()
        for phrase in ("read-only", "outside this candidate", "14 calendar days", "cannot freeze a target"):
            if phrase not in response:
                raise RuntimeError(f"Response firewall omitted in {route_id}: {phrase}")

    send_log = read_tsv(ROOT / "send_log_template.tsv")
    if {row["route_id"] for row in send_log} != set(EXPECTED_ROUTES) or any(row["status"] != "not_sent" for row in send_log):
        raise RuntimeError("Send-log state drift")
    if any(row["recipient_used"] or row["sender"] or row["sent_utc"] or row["message_id"] for row in send_log):
        raise RuntimeError("Send log was prematurely populated")
    firewall = read_tsv(ROOT / "target_firewall_scan.tsv")
    if len(firewall) != 1 or firewall[0]["verdict"] != "pass" or firewall[0]["forbidden_identifier_hits"] != "0":
        raise RuntimeError("Target firewall did not pass")

    for path in ROOT.rglob("*"):
        if path.is_file() and (path.stat().st_mode & 0o222):
            raise RuntimeError(f"Outreach file remains writable: {path}")
    if ROOT.stat().st_mode & 0o222:
        raise RuntimeError("Outreach root remains writable")
    print(
        "PLATFORM_OUTREACH_VALIDATION_PASS "
        f"routes=4 primary=3 supporting=1 generic_questions=12 route_questions=25 "
        f"messages_sent=false targets_disclosed=false seal_sha256={sha256_file(seal_path)}"
    )


if __name__ == "__main__":
    main()
