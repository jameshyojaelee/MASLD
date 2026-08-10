#!/usr/bin/env python3
"""Adjudicate two sealed responses into a still outcome-locked Plan 46 portfolio."""

from __future__ import annotations

import argparse
import datetime as dt
import json
import os
import shutil
import tempfile
from pathlib import Path

from clinical_response_common import (
    PROJECT_ROOT, SCRIPT_ROOT, candidate_root, classify_portfolio, parse_utc,
    read_tsv, resolve_candidate, sha256_file, write_json, write_tsv,
)


AUDIT_FIELDS = (
    "cohort_a_candidate_id", "cohort_b_candidate_id",
    "institutionally_independent", "trial_independent",
    "participant_independent", "overlap_resolved",
    "cohort_b_locked_before_cohort_a_outcomes", "evidence_path", "owner",
    "signed_utc",
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--cohort-a-root", required=True)
    parser.add_argument("--cohort-b-root", required=True)
    parser.add_argument("--audit-root", required=True)
    parser.add_argument("--candidate-id", required=True)
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


def load_response(root: Path) -> tuple[Path, dict[str, object]]:
    seal_path = root / "CLINICAL_RESPONSE_SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_metadata_only_clinical_response_adjudication":
        raise RuntimeError(f"Invalid clinical response release: {root}")
    for field in (
        "participant_rows_present", "clinical_outcomes_present",
        "molecular_outcomes_present", "participant_data_request_authorized",
        "molecular_outcome_access_authorized", "analysis_authorized",
        "paper_promotion_authorized",
    ):
        if seal.get(field) is not False:
            raise RuntimeError(f"Clinical response firewall crossed in {root}: {field}")
    for relative, expected in seal["output_sha256"].items():
        path = root / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Clinical response drift in {root}: {relative}")
    return seal_path, seal


def make_read_only(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o555 if path.is_dir() else 0o444)
    root.chmod(0o555)


def main() -> None:
    args = arguments()
    cohort_a_root = resolve_candidate(args.cohort_a_root)
    cohort_b_root = resolve_candidate(args.cohort_b_root)
    if cohort_a_root == cohort_b_root:
        raise RuntimeError("Clinical portfolio cohorts must differ")
    audit_root = Path(args.audit_root).resolve()
    output_root = candidate_root(args.candidate_id)
    if output_root.exists():
        raise RuntimeError(f"Refusing to overwrite clinical portfolio: {output_root}")
    seal_a_path, seal_a = load_response(cohort_a_root)
    seal_b_path, seal_b = load_response(cohort_b_root)
    audit_path = audit_root / "pairwise_independence_audit.tsv"
    audit_rows = read_tsv(audit_path)
    if len(audit_rows) != 1 or tuple(audit_rows[0]) != AUDIT_FIELDS:
        raise RuntimeError("Clinical pairwise independence audit schema drift")
    audit = audit_rows[0]
    if (
        audit["cohort_a_candidate_id"] != cohort_a_root.name
        or audit["cohort_b_candidate_id"] != cohort_b_root.name
        or not audit["owner"].strip()
        or not audit["signed_utc"].strip()
    ):
        raise RuntimeError("Clinical pairwise independence audit identity drift")
    parse_utc(audit["signed_utc"], "pairwise independence audit")
    independence = {}
    for field in (
        "institutionally_independent", "trial_independent", "participant_independent",
        "overlap_resolved", "cohort_b_locked_before_cohort_a_outcomes",
    ):
        if audit[field].strip().lower() not in {"yes", "no"}:
            raise RuntimeError(f"Invalid clinical independence answer: {field}")
        independence[field] = audit[field].strip().lower() == "yes"

    evidence_register_path = audit_root / "evidence_file_register.tsv"
    evidence_rows = read_tsv(evidence_register_path)
    registered = set()
    for row in evidence_rows:
        relative = Path(row["relative_path"])
        if relative.is_absolute() or ".." in relative.parts or not relative.parts or relative.parts[0] != "evidence":
            raise RuntimeError(f"Unsafe portfolio evidence path: {relative}")
        path = audit_root / relative
        if path.is_symlink() or not path.is_file():
            raise RuntimeError(f"Missing portfolio evidence: {relative}")
        if path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Clinical portfolio evidence drift: {relative}")
        registered.add(relative.as_posix())
    referenced = {part.strip() for part in audit["evidence_path"].split(";") if part.strip()}
    actual = {
        path.relative_to(audit_root).as_posix()
        for path in (audit_root / "evidence").rglob("*") if path.is_file()
    } if (audit_root / "evidence").is_dir() else set()
    if not registered or registered != referenced or registered != actual:
        raise RuntimeError("Clinical portfolio evidence universe mismatch")

    cohort_a = {
        "candidate_cohort_role": seal_a["candidate_cohort_role"],
        "route_id": seal_a["route_id"],
        "independence_group": seal_a["independence_group"],
        "institution_uid": seal_a["institution_uid"],
        "trial_uid": seal_a["trial_uid"],
        "pass_to_blinded_census": seal_a["pass_to_blinded_census"],
        "full_goal_nested_candidate": seal_a["full_goal_nested_candidate"],
    }
    cohort_b = {
        "candidate_cohort_role": seal_b["candidate_cohort_role"],
        "route_id": seal_b["route_id"],
        "independence_group": seal_b["independence_group"],
        "institution_uid": seal_b["institution_uid"],
        "trial_uid": seal_b["trial_uid"],
        "pass_to_blinded_census": seal_b["pass_to_blinded_census"],
        "full_goal_nested_candidate": seal_b["full_goal_nested_candidate"],
    }
    verdict = classify_portfolio(cohort_a, cohort_b, independence)
    output_root.mkdir(parents=True)
    outputs: list[Path] = []
    for source, destination in (
        (audit_path, output_root / "pairwise_independence_audit.tsv"),
        (evidence_register_path, output_root / "evidence_file_register.tsv"),
    ):
        copy_binary(source, destination)
        outputs.append(destination)
    for relative in sorted(registered):
        destination = output_root / relative
        copy_binary(audit_root / relative, destination)
        outputs.append(destination)
    response_rows = [
        {
            "portfolio_role": role,
            "candidate_id": root.name,
            "seal_path": seal_path.relative_to(PROJECT_ROOT),
            "seal_sha256": sha256_file(seal_path),
            "route_id": seal["route_id"],
            "candidate_cohort_role": seal["candidate_cohort_role"],
            "independence_group": seal["independence_group"],
            "institution_uid": seal["institution_uid"],
            "trial_uid": seal["trial_uid"],
            "feasibility_verdict": seal["feasibility_verdict"],
            "full_goal_nested_candidate": seal["full_goal_nested_candidate"],
        }
        for role, root, seal_path, seal in (
            ("cohort_A", cohort_a_root, seal_a_path, seal_a),
            ("cohort_B", cohort_b_root, seal_b_path, seal_b),
        )
    ]
    responses_path = output_root / "selected_clinical_responses.tsv"
    write_tsv(responses_path, response_rows, list(response_rows[0]))
    outputs.append(responses_path)
    adjudication = {
        **verdict,
        "cohort_a_candidate_id": cohort_a_root.name,
        "cohort_b_candidate_id": cohort_b_root.name,
        "n_independent_cohorts": 2 if verdict["two_cohort_source_gate_pass"] else 0,
        "n_nested_full_goal_candidates": sum(
            bool(seal["full_goal_nested_candidate"]) for seal in (seal_a, seal_b)
        ),
    }
    adjudication_path = output_root / "clinical_portfolio_adjudication.tsv"
    write_tsv(adjudication_path, [adjudication], list(adjudication))
    outputs.append(adjudication_path)
    manifest_path = output_root / "clinical_portfolio_source_manifest.tsv"
    source_paths = (
        ("producer", SCRIPT_ROOT / "61_adjudicate_clinical_response_portfolio.py"),
        ("validator", SCRIPT_ROOT / "62_validate_clinical_response_portfolio.py"),
        ("common_rules", SCRIPT_ROOT / "clinical_response_common.py"),
        ("cohort_a_seal", seal_a_path),
        ("cohort_b_seal", seal_b_path),
        ("independence_audit", audit_path),
    )
    write_tsv(manifest_path, [{
        "role": role,
        "path": path.relative_to(PROJECT_ROOT) if PROJECT_ROOT in path.parents else path,
        "sha256": sha256_file(path),
        "size_bytes": path.stat().st_size,
    } for role, path in source_paths], ["role", "path", "sha256", "size_bytes"])
    outputs.append(manifest_path)
    seal = {
        "status": "sealed_outcome_locked_two_cohort_clinical_portfolio",
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        **adjudication,
        "participant_rows_present": False,
        "clinical_outcomes_present": False,
        "molecular_outcomes_present": False,
        "output_sha256": {
            path.relative_to(output_root).as_posix(): sha256_file(path)
            for path in outputs
        },
    }
    seal_path = output_root / "CLINICAL_PORTFOLIO_SEALED.json"
    write_json(seal_path, seal)
    (output_root / "CLINICAL_PORTFOLIO_SEAL_SHA256.txt").write_text(
        sha256_file(seal_path) + "\n", encoding="utf-8"
    )
    make_read_only(output_root)
    print(json.dumps(seal, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

