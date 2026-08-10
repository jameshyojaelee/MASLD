#!/usr/bin/env python3
"""Independently rederive a sealed metadata-only clinical response verdict."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from clinical_response_common import (
    CENSUS_METRICS, PROJECT_ROOT, classify_response, read_tsv, sha256_file,
)


def truth(value: str) -> bool:
    if value == "True":
        return True
    if value == "False":
        return False
    raise RuntimeError(f"Invalid serialized boolean: {value}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True)
    root = Path(parser.parse_args().root).resolve()
    seal_path = root / "CLINICAL_RESPONSE_SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_metadata_only_clinical_response_adjudication":
        raise RuntimeError("Invalid sealed clinical response")
    for field in (
        "participant_rows_present", "clinical_outcomes_present",
        "molecular_outcomes_present", "participant_data_request_authorized",
        "molecular_outcome_access_authorized", "analysis_authorized",
        "paper_promotion_authorized",
    ):
        if seal.get(field) is not False:
            raise RuntimeError(f"Clinical response release firewall crossed: {field}")
    if (root / "CLINICAL_RESPONSE_SEAL_SHA256.txt").read_text(encoding="utf-8").strip() != sha256_file(seal_path):
        raise RuntimeError("Clinical response seal sidecar mismatch")
    for relative, expected in seal["output_sha256"].items():
        path = root / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Clinical response payload drift: {relative}")
    source = read_tsv(root / "clinical_response_source_manifest.tsv")
    if {row["role"] for row in source} != {"producer", "validator", "common_rules", "outreach_seal", "draft_marker"}:
        raise RuntimeError("Clinical response source universe drift")
    for row in source:
        path = PROJECT_ROOT / row["path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Clinical response source drift: {path}")

    response = root / "response"
    questions = read_tsv(response / "clinical_cohort_feasibility_questions.tsv")
    answers = {row["question_id"]: row["answer"].strip().lower() for row in questions}
    assays = {
        row["assay_id"]: row["available"].strip().lower()
        for row in read_tsv(response / "assay_availability_matrix.tsv")
    }
    census = {
        row["metric"]: int(row["value"])
        for row in read_tsv(response / "blinded_cohort_census.tsv")
    }
    if set(census) != set(CENSUS_METRICS):
        raise RuntimeError("Clinical response census drift")
    identity_rows = read_tsv(response / "clinical_response_identity.tsv")
    if len(identity_rows) != 1:
        raise RuntimeError("Clinical response identity cardinality drift")
    identity = identity_rows[0]
    required_identity = (
        "cohort_uid", "institution_uid", "trial_uid", "named_custodian",
        "governance_route", "credible_completion_utc", "prepared_by", "submitted_utc",
    )
    roles = {
        row["role"]: row["assigned_person"].strip()
        for row in read_tsv(response / "clinical_role_firewall.tsv")
    }
    expected = classify_response(
        answers=answers,
        assays=assays,
        census=census,
        identity_complete=all(bool(identity[field].strip()) for field in required_identity),
        custodian_named=bool(roles.get("clinical_custodian")),
        missing_gate_recoverable=identity["single_missing_gate_recoverable"].lower() == "yes",
        missing_gate_id=identity["missing_gate_id"].strip(),
    )
    adjudication_rows = read_tsv(root / "clinical_response_adjudication.tsv")
    if len(adjudication_rows) != 1:
        raise RuntimeError("Clinical response adjudication cardinality drift")
    adjudication = adjudication_rows[0]
    for field, value in expected.items():
        if isinstance(value, str):
            observed = adjudication[field]
        elif isinstance(value, bool):
            observed = truth(adjudication[field])
        else:
            observed = int(adjudication[field])
        if observed != value or seal.get(field) != value:
            raise RuntimeError(f"Clinical response adjudication drift: {field}")
    evidence = read_tsv(response / "evidence_file_register.tsv")
    registered = set()
    for row in evidence:
        path = response / row["relative_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Clinical response evidence drift: {path}")
        registered.add(row["relative_path"])
    referenced = {
        part.strip()
        for name in (
            "clinical_cohort_feasibility_questions.tsv", "assay_availability_matrix.tsv",
            "route_specific_questions.tsv", "blinded_cohort_census.tsv",
        )
        for row in read_tsv(response / name)
        for part in row.get("evidence_path", "").split(";")
        if part.strip()
    }
    if registered != referenced or len(registered) != int(adjudication["n_evidence_files"]):
        raise RuntimeError("Clinical response evidence universe drift")
    for path in root.rglob("*"):
        if path.is_file() and path.stat().st_mode & 0o222:
            raise RuntimeError(f"Clinical response release remains writable: {path}")
    print(
        "CLINICAL_RESPONSE_RELEASE_VALIDATION_PASS "
        f"route={seal['route_id']} verdict={expected['feasibility_verdict']} "
        f"full_nested={str(expected['full_goal_nested_candidate']).lower()} outcomes=false"
    )


if __name__ == "__main__":
    main()
