#!/usr/bin/env python3
"""Seal a completed collaborator response and adjudicate feasibility only."""

from __future__ import annotations

import datetime as dt
import json
import os
import shutil
import tempfile
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_json, read_tsv,
    require_within, sha256_file, write_tsv,
)


TABLES = [
    "platform_feasibility_questions.tsv",
    "role_and_blinding_firewall.tsv",
    "collaborator_deliverable_register.tsv",
    "response_identity.tsv",
    "evidence_file_register.tsv",
]
STATIC_FIELDS = {
    "platform_feasibility_questions.tsv": ["question_id", "domain", "frozen_question", "impact"],
    "role_and_blinding_firewall.tsv": ["role", "responsibility", "target_access", "outcome_access", "prohibition"],
    "collaborator_deliverable_register.tsv": [
        "deliverable_id", "deliverable", "owner_role", "required_before", "immutable_template_paths",
    ],
}
ALLOWED_DELIVERABLE_STATES = {"not_started", "in_progress", "complete", "not_applicable"}
CRITICAL_DISTINCT_ROLES = ["data_manager", "qc_reviewer_1", "qc_reviewer_2", "analysis_lead"]


def candidate_source(variable: str) -> Path:
    raw = os.environ.get(variable, "").strip()
    if not raw:
        raise RuntimeError(f"{variable} is required")
    path = Path(raw)
    path = path if path.is_absolute() else PROJECT_ROOT / path
    path = path.resolve()
    allowed = (PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates").resolve()
    if allowed not in path.parents:
        raise RuntimeError(f"{variable} escapes the candidate root: {path}")
    return path


def parse_utc(value: str, label: str) -> None:
    try:
        parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise RuntimeError(f"Invalid UTC timestamp for {label}: {value}") from error
    if parsed.tzinfo is None or parsed.utcoffset() != dt.timedelta(0):
        raise RuntimeError(f"Timestamp is not UTC for {label}: {value}")


def copy_binary(source: Path, destination: Path) -> None:
    require_within(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{destination.name}.", dir=destination.parent)
    try:
        with source.open("rb") as src, os.fdopen(fd, "wb") as dst:
            shutil.copyfileobj(src, dst)
        os.replace(temporary, destination)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def validate_handoff(root: Path) -> tuple[Path, dict[str, object]]:
    seal_path = root / "STAGE_A_COLLABORATOR_HANDOFF_SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if (
        seal.get("status") != "sealed_target_independent_stage_a_collaborator_handoff"
        or seal.get("bundle_read_only") is not True
        or seal.get("experimental_targets_frozen") is not False
        or seal.get("scientific_outcomes_inspected") is not False
    ):
        raise RuntimeError("Invalid collaborator handoff")
    for relative, expected in seal["output_sha256"].items():  # type: ignore[union-attr]
        path = root / str(relative)
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Collaborator handoff drift: {relative}")
    return seal_path, seal


def exact_static_rows(response: Path, handoff: Path, name: str) -> list[dict[str, str]]:
    observed = read_tsv(response / name)
    frozen = read_tsv(handoff / name)
    fields = STATIC_FIELDS[name]
    if len(observed) != len(frozen):
        raise RuntimeError(f"Response row-count drift: {name}")
    for index, (left, right) in enumerate(zip(observed, frozen), start=1):
        if any(left.get(field) != right.get(field) for field in fields):
            raise RuntimeError(f"Frozen response field drift: {name} row {index}")
        if set(left) != set(right):
            raise RuntimeError(f"Response schema drift: {name}")
    return observed


def classify(answers: dict[str, str], role_gate: bool, d01_complete: bool) -> dict[str, object]:
    stage_a = [answers[f"FQ{i:02d}"] for i in range(1, 11)]
    stage_b = [answers[f"FQ{i:02d}"] for i in range(11, 13)]
    stage_a_feasible = all(value == "yes" for value in stage_a)
    full_mechanism_feasible = stage_a_feasible and all(value == "yes" for value in stage_b)
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
    stage_a_planning = stage_a_feasible and role_gate and d01_complete
    return {
        "feasibility_verdict": verdict,
        "stage_a_platform_feasible": stage_a_feasible,
        "full_mechanism_platform_feasible": full_mechanism_feasible,
        "role_firewall_accepted": role_gate,
        "d01_complete": d01_complete,
        "stage_a_planning_permitted": stage_a_planning,
        "target_disclosure_permitted": full_mechanism_feasible and role_gate and d01_complete,
    }


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite response release: {CANDIDATE_ROOT}")
    handoff_root = candidate_source("PLAN45_COLLABORATOR_HANDOFF_ROOT")
    response_root = candidate_source("PLAN45_COLLABORATOR_RESPONSE_DRAFT_ROOT")
    if response_root == handoff_root:
        raise RuntimeError("Response draft cannot be the sealed handoff")
    handoff_seal_path, _ = validate_handoff(handoff_root)
    draft_marker_path = response_root / "DRAFT_EDITABLE.json"
    marker = json.loads(draft_marker_path.read_text(encoding="utf-8"))
    if (
        marker.get("status") != "editable_collaborator_response_scaffold"
        or marker.get("parent_candidate_id") != handoff_root.name
        or marker.get("parent_seal_sha256") != sha256_file(handoff_seal_path)
        or marker.get("mutable_response_workspace") is not True
        or marker.get("target_identity_present") is not False
        or marker.get("scientific_outcomes_present") is not False
    ):
        raise RuntimeError("Invalid collaborator-response draft provenance")

    allowed_root_files = set(TABLES) | {"DRAFT_EDITABLE.json", "EDITABLE_RESPONSE_README.md"}
    actual_root_files = {path.name for path in response_root.iterdir() if path.is_file()}
    if actual_root_files != allowed_root_files:
        raise RuntimeError(f"Unexpected response-root files: {sorted(actual_root_files ^ allowed_root_files)}")
    unexpected_dirs = {path.name for path in response_root.iterdir() if path.is_dir()} - {"evidence"}
    if unexpected_dirs:
        raise RuntimeError(f"Unexpected response directories: {sorted(unexpected_dirs)}")

    identity = read_tsv(response_root / "response_identity.tsv")
    if len(identity) != 1:
        raise RuntimeError("Response identity must contain one row")
    identity_row = identity[0]
    expected_identity = {
        "response_package_id", "parent_candidate_id", "parent_seal_sha256",
        "organization", "prepared_by", "created_utc", "response_status",
    }
    if set(identity_row) != expected_identity:
        raise RuntimeError("Response identity schema drift")
    if (
        identity_row["response_package_id"] != response_root.name
        or identity_row["parent_candidate_id"] != handoff_root.name
        or identity_row["parent_seal_sha256"] != sha256_file(handoff_seal_path)
        or identity_row["response_status"] != "submitted"
        or not identity_row["organization"].strip()
        or not identity_row["prepared_by"].strip()
    ):
        raise RuntimeError("Incomplete or inconsistent response identity")
    parse_utc(identity_row["created_utc"], "response identity")

    questions = exact_static_rows(response_root, handoff_root, TABLES[0])
    answers: dict[str, str] = {}
    referenced_evidence: set[str] = set()
    for row in questions:
        answer = row["answer"].strip().lower()
        if answer not in {"yes", "no", "unresolved"}:
            raise RuntimeError(f"Invalid feasibility answer: {row['question_id']}")
        if not row["owner"].strip() or not row["signed_utc"].strip() or not row["evidence_path"].strip():
            raise RuntimeError(f"Incomplete feasibility attestation: {row['question_id']}")
        parse_utc(row["signed_utc"], row["question_id"])
        answers[row["question_id"]] = answer
        referenced_evidence.update(part.strip() for part in row["evidence_path"].split(";") if part.strip())

    roles = exact_static_rows(response_root, handoff_root, TABLES[1])
    assignments = {}
    for row in roles:
        person = row["assigned_person"].strip()
        if not person or not row["accepted_utc"].strip():
            raise RuntimeError(f"Incomplete role acceptance: {row['role']}")
        parse_utc(row["accepted_utc"], row["role"])
        assignments[row["role"]] = person.casefold()
    critical_people = [assignments[role] for role in CRITICAL_DISTINCT_ROLES]
    role_gate = len(set(critical_people)) == len(critical_people)
    if assignments["independent_validator"] == assignments["analysis_lead"]:
        role_gate = False

    deliverables = exact_static_rows(response_root, handoff_root, TABLES[2])
    if any(row["status"] not in ALLOWED_DELIVERABLE_STATES for row in deliverables):
        raise RuntimeError("Invalid collaborator-deliverable status")
    deliverable_state = {row["deliverable_id"]: row["status"] for row in deliverables}
    d01_complete = deliverable_state.get("D01") == "complete"

    evidence_rows = read_tsv(response_root / "evidence_file_register.tsv")
    expected_evidence_fields = {"evidence_id", "relative_path", "description", "sha256", "size_bytes"}
    evidence_paths: set[str] = set()
    evidence_ids: set[str] = set()
    for row in evidence_rows:
        if set(row) != expected_evidence_fields:
            raise RuntimeError("Evidence register schema drift")
        relative = Path(row["relative_path"])
        if relative.is_absolute() or ".." in relative.parts or not relative.parts or relative.parts[0] != "evidence":
            raise RuntimeError(f"Unsafe evidence path: {relative}")
        if row["evidence_id"] in evidence_ids or str(relative) in evidence_paths:
            raise RuntimeError("Duplicate evidence ID or path")
        evidence_ids.add(row["evidence_id"])
        evidence_paths.add(str(relative))
        source = response_root / relative
        if source.is_symlink() or not source.is_file():
            raise RuntimeError(f"Missing or linked evidence: {relative}")
        if source.stat().st_size != int(row["size_bytes"]) or sha256_file(source) != row["sha256"]:
            raise RuntimeError(f"Evidence identity drift: {relative}")
        if not row["description"].strip():
            raise RuntimeError(f"Evidence description missing: {relative}")
    actual_evidence = {
        str(path.relative_to(response_root))
        for path in (response_root / "evidence").rglob("*")
        if path.is_file()
    } if (response_root / "evidence").is_dir() else set()
    if actual_evidence != evidence_paths or not referenced_evidence <= evidence_paths:
        raise RuntimeError("Evidence register/reference universe mismatch")
    if referenced_evidence != evidence_paths:
        raise RuntimeError("Every evidence object must support at least one feasibility answer")

    verdict = classify(answers, role_gate, d01_complete)
    CANDIDATE_ROOT.mkdir(parents=True)
    outputs: list[Path] = []
    for name in TABLES:
        destination = CANDIDATE_ROOT / "response" / name
        copy_binary(response_root / name, destination)
        outputs.append(destination)
    for relative in sorted(evidence_paths):
        destination = CANDIDATE_ROOT / "response" / relative
        copy_binary(response_root / relative, destination)
        outputs.append(destination)
    marker_copy = CANDIDATE_ROOT / "response" / "source_DRAFT_EDITABLE.json"
    copy_binary(draft_marker_path, marker_copy)
    outputs.append(marker_copy)

    adjudication_path = CANDIDATE_ROOT / "collaborator_response_adjudication.tsv"
    adjudication = {
        **verdict,
        "n_yes": sum(value == "yes" for value in answers.values()),
        "n_no": sum(value == "no" for value in answers.values()),
        "n_unresolved": sum(value == "unresolved" for value in answers.values()),
        "n_evidence_files": len(evidence_paths),
        "target_frozen": False,
        "stage_a_outcome_generation_permitted": False,
        "stage_b_frozen": False,
        "paper_promotion_authorized": False,
    }
    write_tsv(adjudication_path, [adjudication], list(adjudication))
    outputs.append(adjudication_path)

    source_manifest_path = CANDIDATE_ROOT / "collaborator_response_source_manifest.tsv"
    source_rows = []
    for role, path in [
        ("response_producer", SCRIPT_ROOT / "47_seal_and_adjudicate_collaborator_response.py"),
        ("response_validator", SCRIPT_ROOT / "48_validate_collaborator_response_release.py"),
        ("handoff_seal", handoff_seal_path),
        ("response_draft_marker", draft_marker_path),
    ]:
        source_rows.append({
            "role": role,
            "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    write_tsv(source_manifest_path, source_rows, ["role", "source_path", "size_bytes", "sha256"])
    outputs.append(source_manifest_path)

    seal = {
        "status": "sealed_collaborator_response_adjudication",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "parent_handoff_candidate_id": handoff_root.name,
        "parent_handoff_seal_sha256": sha256_file(handoff_seal_path),
        "source_response_draft_id": response_root.name,
        **verdict,
        "n_feasibility_questions": len(answers),
        "n_evidence_files": len(evidence_paths),
        "target_identity_present": False,
        "target_frozen": False,
        "scientific_outcomes_present": False,
        "stage_a_outcome_generation_permitted": False,
        "stage_b_frozen": False,
        "paper_promotion_authorized": False,
        "output_sha256": {
            str(path.relative_to(CANDIDATE_ROOT)): sha256_file(path) for path in outputs
        },
    }
    atomic_write_json(CANDIDATE_ROOT / "COLLABORATOR_RESPONSE_SEALED.json", seal)
    print(json.dumps(seal, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
