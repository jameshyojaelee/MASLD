#!/usr/bin/env python3
"""Seal and outcome-blindly adjudicate a completed Plan 46 clinical response."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import tempfile
from pathlib import Path

from clinical_response_common import (
    ASSAY_IDS, CENSUS_METRICS, DEFAULT_OUTREACH_ROOT, PROJECT_ROOT, ROUTE_ROLES,
    SCRIPT_ROOT, candidate_root, classify_response, nonnegative_integer,
    parse_utc, read_tsv, resolve_candidate, sha256_file, validate_outreach,
    write_json, write_tsv,
)


TABLES = (
    "clinical_cohort_feasibility_questions.tsv",
    "assay_availability_matrix.tsv",
    "clinical_role_firewall.tsv",
    "minimum_data_dictionary.tsv",
    "route_specific_questions.tsv",
    "public_inventory_audit.tsv",
    "clinical_response_identity.tsv",
    "blinded_cohort_census.tsv",
    "evidence_file_register.tsv",
)
STATIC_FIELDS = {
    "clinical_cohort_feasibility_questions.tsv": (
        "question_id", "domain", "frozen_question", "impact",
    ),
    "assay_availability_matrix.tsv": (
        "assay_id", "biospecimen", "requirement", "scientific_role",
    ),
    "clinical_role_firewall.tsv": ("role", "responsibility", "prohibition"),
    "route_specific_questions.tsv": (
        "route_question_id", "frozen_cf_links", "metadata_only_question",
    ),
    "public_inventory_audit.tsv": (
        "inventory_item", "public_value", "evidence_state", "source_url",
    ),
}


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--draft-root", required=True)
    parser.add_argument("--candidate-id", required=True)
    parser.add_argument("--outreach-root", default=str(DEFAULT_OUTREACH_ROOT))
    return parser.parse_args()


def copy_binary(source: Path, destination: Path) -> None:
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


def exact_static(response: Path, frozen: Path, name: str) -> list[dict[str, str]]:
    observed = read_tsv(response / name)
    source = read_tsv(frozen / name)
    if len(observed) != len(source):
        raise RuntimeError(f"Clinical response row-count drift: {name}")
    fields = STATIC_FIELDS[name]
    for index, (left, right) in enumerate(zip(observed, source), start=1):
        if set(left) != set(right) or any(left[field] != right[field] for field in fields):
            raise RuntimeError(f"Frozen clinical response field drift: {name} row {index}")
    return observed


def evidence_references(rows: list[dict[str, str]]) -> set[str]:
    return {
        part.strip()
        for row in rows
        for part in row.get("evidence_path", "").split(";")
        if part.strip()
    }


def make_read_only(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o555 if path.is_dir() else 0o444)
    root.chmod(0o555)


def main() -> None:
    args = arguments()
    output_root = candidate_root(args.candidate_id)
    if output_root.exists():
        raise RuntimeError(f"Refusing to overwrite clinical response release: {output_root}")
    draft_root = resolve_candidate(args.draft_root)
    outreach_root = Path(args.outreach_root).resolve()
    outreach_seal_path, _, routes = validate_outreach(outreach_root)
    marker_path = draft_root / "DRAFT_EDITABLE.json"
    marker = json.loads(marker_path.read_text(encoding="utf-8"))
    if (
        marker.get("status") != "editable_metadata_only_clinical_response"
        or marker.get("parent_outreach_candidate_id") != outreach_root.name
        or marker.get("parent_outreach_seal_sha256") != sha256_file(outreach_seal_path)
        or marker.get("mutable_response_workspace") is not True
    ):
        raise RuntimeError("Invalid clinical response draft provenance")
    for field in (
        "participant_rows_present", "clinical_outcomes_present",
        "molecular_outcomes_present", "permission_assumed",
        "analysis_authorized", "paper_promotion_authorized",
    ):
        if marker.get(field) is not False:
            raise RuntimeError(f"Clinical response draft firewall crossed: {field}")
    route_id = str(marker["route_id"])
    if route_id not in routes:
        raise RuntimeError("Clinical response route is not frozen")
    frozen_root = outreach_root / "routes" / route_id

    expected_root_files = set(TABLES) | {
        "DRAFT_EDITABLE.json", "EDITABLE_CLINICAL_RESPONSE_README.md",
    }
    actual_root_files = {path.name for path in draft_root.iterdir() if path.is_file()}
    if actual_root_files != expected_root_files:
        raise RuntimeError(f"Unexpected clinical response files: {sorted(actual_root_files ^ expected_root_files)}")
    unexpected_dirs = {path.name for path in draft_root.iterdir() if path.is_dir()} - {"evidence"}
    if unexpected_dirs:
        raise RuntimeError(f"Unexpected clinical response directories: {sorted(unexpected_dirs)}")

    identity_rows = read_tsv(draft_root / "clinical_response_identity.tsv")
    if len(identity_rows) != 1:
        raise RuntimeError("Clinical response identity must contain one row")
    identity = identity_rows[0]
    fixed_identity = {
        "response_package_id": draft_root.name,
        "parent_outreach_candidate_id": outreach_root.name,
        "parent_outreach_seal_sha256": sha256_file(outreach_seal_path),
        "route_id": route_id,
        "candidate_cohort_role": routes[route_id]["candidate_cohort_role"],
        "independence_group": routes[route_id]["independence_group"],
        "organization": routes[route_id]["organization"],
    }
    if any(identity.get(field) != value for field, value in fixed_identity.items()):
        raise RuntimeError("Clinical response identity provenance drift")
    required_identity = (
        "cohort_uid", "institution_uid", "trial_uid", "named_custodian",
        "governance_route", "credible_completion_utc", "prepared_by", "submitted_utc",
    )
    if identity.get("response_status") != "submitted" or any(not identity.get(field, "").strip() for field in required_identity):
        raise RuntimeError("Clinical response identity is incomplete")
    parse_utc(identity["credible_completion_utc"], "credible completion")
    parse_utc(identity["submitted_utc"], "response submission")
    recoverable = identity["single_missing_gate_recoverable"].strip().lower()
    if recoverable not in {"yes", "no"}:
        raise RuntimeError("Invalid single-missing-gate recoverability")
    if recoverable == "yes":
        if not identity["missing_gate_id"].strip() or not identity["recovery_due_utc"].strip():
            raise RuntimeError("Recoverable gate lacks ID or deadline")
        parse_utc(identity["recovery_due_utc"], "gate recovery deadline")
    elif identity["missing_gate_id"].strip() or identity["recovery_due_utc"].strip():
        raise RuntimeError("Non-recoverable response improperly names a recovery gate")

    questions = exact_static(draft_root, frozen_root, "clinical_cohort_feasibility_questions.tsv")
    answers: dict[str, str] = {}
    for row in questions:
        answer = row["answer"].strip().lower()
        if answer not in {"yes", "no", "unresolved"}:
            raise RuntimeError(f"Invalid clinical answer: {row['question_id']}")
        if row["cohort_role_candidate"] != routes[route_id]["candidate_cohort_role"]:
            raise RuntimeError(f"Clinical cohort role not attested: {row['question_id']}")
        if not row["evidence_path"].strip() or not row["owner"].strip() or not row["signed_utc"].strip():
            raise RuntimeError(f"Incomplete clinical attestation: {row['question_id']}")
        parse_utc(row["signed_utc"], row["question_id"])
        answers[row["question_id"]] = answer

    assays_rows = exact_static(draft_root, frozen_root, "assay_availability_matrix.tsv")
    assays: dict[str, str] = {}
    for row in assays_rows:
        available = row["available"].strip().lower()
        if available not in {"yes", "no", "unresolved"} or not row["evidence_path"].strip():
            raise RuntimeError(f"Incomplete clinical assay response: {row['assay_id']}")
        baseline = nonnegative_integer(row["n_baseline"], f"{row['assay_id']} baseline")
        followup = nonnegative_integer(row["n_followup"], f"{row['assay_id']} followup")
        pairs = nonnegative_integer(row["n_authoritative_pairs"], f"{row['assay_id']} pairs")
        if pairs > min(baseline, followup) or (available == "yes" and pairs == 0):
            raise RuntimeError(f"Inconsistent clinical assay counts: {row['assay_id']}")
        assays[row["assay_id"]] = available
    if set(assays) != set(ASSAY_IDS):
        raise RuntimeError("Clinical assay universe drift")

    roles = exact_static(draft_root, frozen_root, "clinical_role_firewall.tsv")
    assignments: dict[str, str] = {}
    for row in roles:
        person = row["assigned_person"].strip()
        timestamp = row["accepted_utc"].strip()
        if bool(person) != bool(timestamp):
            raise RuntimeError(f"Partial clinical role acceptance: {row['role']}")
        if timestamp:
            parse_utc(timestamp, row["role"])
        assignments[row["role"]] = person.casefold()
    if not assignments.get("clinical_custodian"):
        raise RuntimeError("A named clinical custodian is required")
    if assignments.get("data_manager") and assignments.get("data_manager") == assignments.get("analysis_lead"):
        raise RuntimeError("Data manager and analysis lead must differ")
    if assignments.get("independent_validator") and assignments.get("independent_validator") == assignments.get("analysis_lead"):
        raise RuntimeError("Independent validator and analysis lead must differ")

    if (draft_root / "minimum_data_dictionary.tsv").read_bytes() != (frozen_root / "minimum_data_dictionary.tsv").read_bytes():
        raise RuntimeError("Clinical minimum data dictionary drift")
    route_questions = exact_static(draft_root, frozen_root, "route_specific_questions.tsv")
    for row in route_questions:
        if row["answer"].strip().lower() not in {"yes", "no", "unresolved"}:
            raise RuntimeError(f"Invalid route answer: {row['route_question_id']}")
        if not row["evidence_path"].strip() or not row["owner"].strip() or not row["signed_utc"].strip():
            raise RuntimeError(f"Incomplete route attestation: {row['route_question_id']}")
        parse_utc(row["signed_utc"], row["route_question_id"])
    inventory = exact_static(draft_root, frozen_root, "public_inventory_audit.tsv")
    if any(not row["custodian_confirmation"].strip() for row in inventory):
        raise RuntimeError("Every public/unverified inventory row needs custodian confirmation")

    census_rows = read_tsv(draft_root / "blinded_cohort_census.tsv")
    if [row["metric"] for row in census_rows] != list(CENSUS_METRICS):
        raise RuntimeError("Clinical blinded census universe drift")
    census = {}
    for row in census_rows:
        if not row["evidence_path"].strip():
            raise RuntimeError(f"Clinical census evidence missing: {row['metric']}")
        census[row["metric"]] = nonnegative_integer(row["value"], row["metric"])
    pairs = census["authoritative_baseline_followup_pairs"]
    if (
        pairs > census["unique_participants"]
        or census["complete_histology_pairs"] > pairs
        or census["paired_genomewide_rna_pairs"] > pairs
        or census["histologic_improvers"] + census["histologic_non_improvers"] > census["complete_histology_pairs"]
        or census["five_assay_complete_pairs"] > min(
            pairs, census["complete_histology_pairs"], census["paired_genomewide_rna_pairs"]
        )
    ):
        raise RuntimeError("Clinical blinded census is internally inconsistent")

    evidence_rows = read_tsv(draft_root / "evidence_file_register.tsv")
    registered: set[str] = set()
    for row in evidence_rows:
        relative = Path(row["relative_path"])
        if relative.is_absolute() or ".." in relative.parts or not relative.parts or relative.parts[0] != "evidence":
            raise RuntimeError(f"Unsafe clinical evidence path: {relative}")
        relative_text = relative.as_posix()
        if relative_text in registered or not row["evidence_id"].strip() or not row["description"].strip():
            raise RuntimeError("Duplicate or incomplete clinical evidence register")
        path = draft_root / relative
        if path.is_symlink() or not path.is_file():
            raise RuntimeError(f"Missing or linked clinical evidence: {relative}")
        if path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Clinical evidence identity drift: {relative}")
        registered.add(relative_text)
    referenced = set().union(
        evidence_references(questions), evidence_references(assays_rows),
        evidence_references(route_questions), evidence_references(census_rows),
    )
    actual = {
        path.relative_to(draft_root).as_posix()
        for path in (draft_root / "evidence").rglob("*") if path.is_file()
    } if (draft_root / "evidence").is_dir() else set()
    if not registered or registered != referenced or registered != actual:
        raise RuntimeError("Clinical evidence register/reference universe mismatch")

    verdict = classify_response(
        answers=answers,
        assays=assays,
        census=census,
        identity_complete=all(bool(identity[field].strip()) for field in required_identity),
        custodian_named=bool(assignments["clinical_custodian"]),
        missing_gate_recoverable=recoverable == "yes",
        missing_gate_id=identity["missing_gate_id"].strip(),
    )
    output_root.mkdir(parents=True)
    outputs: list[Path] = []
    for name in TABLES:
        destination = output_root / "response" / name
        copy_binary(draft_root / name, destination)
        outputs.append(destination)
    for relative in sorted(registered):
        destination = output_root / "response" / relative
        copy_binary(draft_root / relative, destination)
        outputs.append(destination)
    marker_copy = output_root / "response" / "source_DRAFT_EDITABLE.json"
    copy_binary(marker_path, marker_copy)
    outputs.append(marker_copy)

    adjudication = {
        "route_id": route_id,
        "candidate_cohort_role": routes[route_id]["candidate_cohort_role"],
        "independence_group": routes[route_id]["independence_group"],
        "cohort_uid": identity["cohort_uid"],
        "institution_uid": identity["institution_uid"],
        "trial_uid": identity["trial_uid"],
        **verdict,
        "n_authoritative_pairs": pairs,
        "n_complete_histology_pairs": census["complete_histology_pairs"],
        "n_paired_genomewide_rna_pairs": census["paired_genomewide_rna_pairs"],
        "n_histologic_improvers": census["histologic_improvers"],
        "n_histologic_non_improvers": census["histologic_non_improvers"],
        "n_five_assay_complete_pairs": census["five_assay_complete_pairs"],
        "n_evidence_files": len(registered),
    }
    adjudication_path = output_root / "clinical_response_adjudication.tsv"
    write_tsv(adjudication_path, [adjudication], list(adjudication))
    outputs.append(adjudication_path)
    manifest_path = output_root / "clinical_response_source_manifest.tsv"
    source_paths = (
        ("producer", SCRIPT_ROOT / "59_seal_and_adjudicate_clinical_response.py"),
        ("validator", SCRIPT_ROOT / "60_validate_clinical_response_release.py"),
        ("common_rules", SCRIPT_ROOT / "clinical_response_common.py"),
        ("outreach_seal", outreach_seal_path),
        ("draft_marker", marker_path),
    )
    write_tsv(manifest_path, [{
        "role": role,
        "path": path.relative_to(PROJECT_ROOT),
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    } for role, path in source_paths], ["role", "path", "sha256", "size_bytes"])
    outputs.append(manifest_path)
    seal = {
        "status": "sealed_metadata_only_clinical_response_adjudication",
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "parent_outreach_candidate_id": outreach_root.name,
        "parent_outreach_seal_sha256": sha256_file(outreach_seal_path),
        "source_response_draft_id": draft_root.name,
        **adjudication,
        "participant_rows_present": False,
        "clinical_outcomes_present": False,
        "molecular_outcomes_present": False,
        "output_sha256": {
            path.relative_to(output_root).as_posix(): sha256_file(path)
            for path in outputs
        },
    }
    seal_path = output_root / "CLINICAL_RESPONSE_SEALED.json"
    write_json(seal_path, seal)
    (output_root / "CLINICAL_RESPONSE_SEAL_SHA256.txt").write_text(
        sha256_file(seal_path) + "\n", encoding="utf-8"
    )
    make_read_only(output_root)
    print(json.dumps(seal, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

