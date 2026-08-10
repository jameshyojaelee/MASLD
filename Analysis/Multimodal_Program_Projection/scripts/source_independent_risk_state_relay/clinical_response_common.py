#!/usr/bin/env python3
"""Shared fail-closed rules for metadata-only Plan 46 clinical responses."""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
import re
import tempfile
from pathlib import Path
from typing import Iterable, Mapping, Sequence


PROJECT_ROOT = Path(__file__).resolve().parents[4]
SCRIPT_ROOT = Path(__file__).resolve().parent
CANDIDATE_PARENT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates"
DEFAULT_OUTREACH_ROOT = CANDIDATE_PARENT / "source-independent-risk-state-relay-clinical-outreach-v1-2026-08-10"
CORE_QUESTIONS = tuple(
    [f"CF{i:02d}" for i in range(1, 9)]
    + [f"CF{i:02d}" for i in range(13, 17)]
)
ORTHOGONAL_QUESTIONS = tuple(f"CF{i:02d}" for i in range(9, 13))
ASSAY_IDS = (
    "bulk_rna",
    "histology",
    "tissue_proteomics",
    "secreted_proteomics",
    "spatial_transcriptomics_or_proteomics",
)
CENSUS_METRICS = (
    "unique_participants",
    "authoritative_baseline_followup_pairs",
    "complete_histology_pairs",
    "paired_genomewide_rna_pairs",
    "histologic_improvers",
    "histologic_non_improvers",
    "five_assay_complete_pairs",
)
ROUTE_ROLES = {
    "nash_crn_as116": "candidate_cohort_A",
    "maestro_nash": "candidate_cohort_B",
    "essence": "candidate_cohort_B",
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


def write_tsv(path: Path, rows: Iterable[Mapping[str, object]], fields: Sequence[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in rows:
                writer.writerow({field: row.get(field, "") for field in fields})
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def write_json(path: Path, value: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            json.dump(value, handle, indent=2, sort_keys=True)
            handle.write("\n")
        os.replace(temporary, path)
    except Exception:
        try:
            os.unlink(temporary)
        except FileNotFoundError:
            pass
        raise


def parse_utc(value: str, label: str) -> dt.datetime:
    try:
        parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as error:
        raise RuntimeError(f"Invalid UTC timestamp for {label}: {value}") from error
    if parsed.tzinfo is None or parsed.utcoffset() != dt.timedelta(0):
        raise RuntimeError(f"Timestamp is not UTC for {label}: {value}")
    return parsed


def candidate_root(candidate_id: str) -> Path:
    if not re.fullmatch(r"source-independent-risk-state-relay-clinical-[A-Za-z0-9._-]+", candidate_id):
        raise RuntimeError(f"Invalid clinical candidate ID: {candidate_id!r}")
    return CANDIDATE_PARENT / candidate_id


def resolve_candidate(path: str | Path) -> Path:
    resolved = Path(path)
    resolved = resolved if resolved.is_absolute() else PROJECT_ROOT / resolved
    resolved = resolved.resolve()
    if CANDIDATE_PARENT.resolve() not in resolved.parents:
        raise RuntimeError(f"Path escapes candidate root: {resolved}")
    return resolved


def validate_outreach(root: Path) -> tuple[Path, dict[str, object], dict[str, dict[str, str]]]:
    seal_path = root / "CLINICAL_OUTREACH_SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_metadata_only_clinical_outreach":
        raise RuntimeError("Invalid clinical outreach parent")
    for field in (
        "message_sent", "participant_data_accessed", "clinical_outcomes_accessed",
        "molecular_outcomes_accessed", "permission_granted", "analysis_authorized",
        "paper_promotion_authorized",
    ):
        if seal.get(field) is not False:
            raise RuntimeError(f"Clinical outreach firewall crossed: {field}")
    for relative, expected in seal["output_sha256"].items():
        path = root / relative
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Clinical outreach drift: {relative}")
    routes = {row["route_id"]: row for row in read_tsv(root / "route_registry.tsv")}
    if set(routes) != set(ROUTE_ROLES):
        raise RuntimeError("Clinical outreach route universe drift")
    return seal_path, seal, routes


def nonnegative_integer(value: str, label: str) -> int:
    if not re.fullmatch(r"0|[1-9][0-9]*", value.strip()):
        raise RuntimeError(f"Expected non-negative integer for {label}: {value!r}")
    return int(value)


def classify_response(
    answers: Mapping[str, str],
    assays: Mapping[str, str],
    census: Mapping[str, int],
    identity_complete: bool,
    custodian_named: bool,
    missing_gate_recoverable: bool = False,
    missing_gate_id: str = "",
) -> dict[str, object]:
    """Adjudicate availability only; never authorize outcome access."""
    if set(answers) != {f"CF{i:02d}" for i in range(1, 17)}:
        raise RuntimeError("Clinical feasibility question universe drift")
    if not set(answers.values()) <= {"yes", "no", "unresolved"}:
        raise RuntimeError("Invalid clinical feasibility answer")
    if set(assays) != set(ASSAY_IDS) or not set(assays.values()) <= {"yes", "no", "unresolved"}:
        raise RuntimeError("Clinical assay universe drift")
    if set(census) != set(CENSUS_METRICS) or any(value < 0 for value in census.values()):
        raise RuntimeError("Clinical census universe drift")

    core_no = [qid for qid in CORE_QUESTIONS if answers[qid] == "no"]
    core_unresolved = [qid for qid in CORE_QUESTIONS if answers[qid] == "unresolved"]
    counts_consistent = all(
        census[field] > 0
        for field in (
            "authoritative_baseline_followup_pairs", "complete_histology_pairs",
            "paired_genomewide_rna_pairs", "histologic_improvers",
            "histologic_non_improvers",
        )
    )
    required_assays = assays["bulk_rna"] == "yes" and assays["histology"] == "yes"
    structural_gate = identity_complete and custodian_named and counts_consistent and required_assays
    all_core_yes = not core_no and not core_unresolved
    pass_to_census = all_core_yes and structural_gate

    one_recoverable = (
        not core_no
        and len(core_unresolved) == 1
        and missing_gate_recoverable
        and missing_gate_id == core_unresolved[0]
        and identity_complete
        and custodian_named
    )
    orthogonal_available = any(
        assays[assay] == "yes"
        for assay in (
            "tissue_proteomics", "secreted_proteomics",
            "spatial_transcriptomics_or_proteomics",
        )
    )
    supportive_only = (
        answers["CF05"] == "no"
        and all(answers[qid] == "yes" for qid in ("CF01", "CF02", "CF03", "CF04", "CF14"))
        and orthogonal_available
    )

    if pass_to_census:
        verdict = "pass_to_blinded_census"
    elif one_recoverable:
        verdict = "conditional_missing_one_gate"
    elif supportive_only:
        verdict = "supportive_only"
    else:
        verdict = "fail_inaccessible_or_unpaired"

    full_nested = (
        pass_to_census
        and all(answers[qid] == "yes" for qid in ORTHOGONAL_QUESTIONS)
        and all(assays[assay] == "yes" for assay in ASSAY_IDS)
        and census["five_assay_complete_pairs"] > 0
    )
    return {
        "feasibility_verdict": verdict,
        "pass_to_blinded_census": pass_to_census,
        "full_goal_nested_candidate": full_nested,
        "n_core_no": len(core_no),
        "n_core_unresolved": len(core_unresolved),
        "structural_gate_pass": structural_gate,
        "advance_to_overlap_power_audit": pass_to_census,
        "participant_data_request_authorized": False,
        "molecular_outcome_access_authorized": False,
        "analysis_authorized": False,
        "paper_promotion_authorized": False,
    }


def classify_portfolio(
    cohort_a: Mapping[str, object],
    cohort_b: Mapping[str, object],
    independence: Mapping[str, bool],
) -> dict[str, object]:
    required_independence = (
        "institutionally_independent", "trial_independent", "participant_independent",
        "overlap_resolved", "cohort_b_locked_before_cohort_a_outcomes",
    )
    if set(independence) != set(required_independence):
        raise RuntimeError("Clinical portfolio independence universe drift")
    roles_ok = (
        cohort_a.get("candidate_cohort_role") == "candidate_cohort_A"
        and cohort_b.get("candidate_cohort_role") == "candidate_cohort_B"
    )
    routes_distinct = (
        cohort_a.get("route_id") != cohort_b.get("route_id")
        and cohort_a.get("independence_group") != cohort_b.get("independence_group")
        and cohort_a.get("institution_uid") != cohort_b.get("institution_uid")
        and cohort_a.get("trial_uid") != cohort_b.get("trial_uid")
    )
    both_pass = bool(cohort_a.get("pass_to_blinded_census")) and bool(cohort_b.get("pass_to_blinded_census"))
    independent = all(independence.values())
    source_gate = roles_ok and routes_distinct and both_pass and independent
    full_nested = source_gate and (
        bool(cohort_a.get("full_goal_nested_candidate"))
        or bool(cohort_b.get("full_goal_nested_candidate"))
    )
    return {
        "portfolio_verdict": (
            "pass_two_cohort_blinded_preflight" if source_gate
            else "fail_two_cohort_source_gate"
        ),
        "two_cohort_source_gate_pass": source_gate,
        "orthogonal_nested_source_available": full_nested,
        "advance_to_blinded_overlap_power_model_freeze": source_gate,
        "participant_data_request_authorized": False,
        "molecular_outcome_access_authorized": False,
        "analysis_authorized": False,
        "paper_promotion_authorized": False,
    }

