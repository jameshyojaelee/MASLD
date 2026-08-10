#!/usr/bin/env python3
"""Independently rederive a sealed, still outcome-locked clinical portfolio."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from clinical_response_common import PROJECT_ROOT, classify_portfolio, read_tsv, sha256_file


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
    seal_path = root / "CLINICAL_PORTFOLIO_SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_outcome_locked_two_cohort_clinical_portfolio":
        raise RuntimeError("Invalid clinical portfolio release")
    for field in (
        "participant_rows_present", "clinical_outcomes_present",
        "molecular_outcomes_present", "participant_data_request_authorized",
        "molecular_outcome_access_authorized", "analysis_authorized",
        "paper_promotion_authorized",
    ):
        if seal.get(field) is not False:
            raise RuntimeError(f"Clinical portfolio firewall crossed: {field}")
    if (root / "CLINICAL_PORTFOLIO_SEAL_SHA256.txt").read_text(encoding="utf-8").strip() != sha256_file(seal_path):
        raise RuntimeError("Clinical portfolio seal sidecar mismatch")
    for relative, expected in seal["output_sha256"].items():
        path = root / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Clinical portfolio payload drift: {relative}")

    response_rows = read_tsv(root / "selected_clinical_responses.tsv")
    if [row["portfolio_role"] for row in response_rows] != ["cohort_A", "cohort_B"]:
        raise RuntimeError("Clinical portfolio response-role drift")
    loaded = []
    for row in response_rows:
        response_seal_path = PROJECT_ROOT / row["seal_path"]
        if not response_seal_path.is_file() or sha256_file(response_seal_path) != row["seal_sha256"]:
            raise RuntimeError(f"Clinical portfolio response seal drift: {response_seal_path}")
        response_seal = json.loads(response_seal_path.read_text(encoding="utf-8"))
        loaded.append({
            "candidate_cohort_role": response_seal["candidate_cohort_role"],
            "route_id": response_seal["route_id"],
            "independence_group": response_seal["independence_group"],
            "institution_uid": response_seal["institution_uid"],
            "trial_uid": response_seal["trial_uid"],
            "pass_to_blinded_census": response_seal["pass_to_blinded_census"],
            "full_goal_nested_candidate": response_seal["full_goal_nested_candidate"],
        })
    audit_rows = read_tsv(root / "pairwise_independence_audit.tsv")
    if len(audit_rows) != 1:
        raise RuntimeError("Clinical portfolio audit cardinality drift")
    audit = audit_rows[0]
    independence = {
        field: audit[field].lower() == "yes"
        for field in (
            "institutionally_independent", "trial_independent", "participant_independent",
            "overlap_resolved", "cohort_b_locked_before_cohort_a_outcomes",
        )
    }
    expected = classify_portfolio(loaded[0], loaded[1], independence)
    adjudication_rows = read_tsv(root / "clinical_portfolio_adjudication.tsv")
    if len(adjudication_rows) != 1:
        raise RuntimeError("Clinical portfolio adjudication cardinality drift")
    adjudication = adjudication_rows[0]
    for field, value in expected.items():
        observed = adjudication[field] if isinstance(value, str) else truth(adjudication[field])
        if observed != value or seal.get(field) != value:
            raise RuntimeError(f"Clinical portfolio adjudication drift: {field}")
    source = read_tsv(root / "clinical_portfolio_source_manifest.tsv")
    if {row["role"] for row in source} != {
        "producer", "validator", "common_rules", "cohort_a_seal",
        "cohort_b_seal", "independence_audit",
    }:
        raise RuntimeError("Clinical portfolio source universe drift")
    for row in source:
        path = Path(row["path"])
        path = path if path.is_absolute() else PROJECT_ROOT / path
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Clinical portfolio source drift: {path}")
    for path in root.rglob("*"):
        if path.is_file() and path.stat().st_mode & 0o222:
            raise RuntimeError(f"Clinical portfolio remains writable: {path}")
    print(
        "CLINICAL_PORTFOLIO_VALIDATION_PASS "
        f"verdict={expected['portfolio_verdict']} "
        f"orthogonal_nested={str(expected['orthogonal_nested_source_available']).lower()} outcomes=false"
    )


if __name__ == "__main__":
    main()

