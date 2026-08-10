#!/usr/bin/env python3
"""Independently validate the Hotspot v2 post-freeze semantic adjudication."""

from __future__ import annotations

import argparse
import csv
import hashlib
import math
import os
import tempfile
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable, Mapping, Sequence


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SCRIPT_DIR.parents[3]
RELEASE_ID = "program-context-v2-candidate-2026-08-07"
CANDIDATE_REL = Path(
    "Analysis/Multimodal_Program_Projection/candidates"
) / RELEASE_ID / "hotspot"
EXPECTED_REGISTRY_SHA256 = (
    "b134b5f8e46a716e857af22271aa54e1294432b9f9ab735559d7f603cdd29fe7"
)
EXPECTED_STATE_COUNTS = {
    "indeterminate_nonconfirmatory": 110,
    "indeterminate_stability_failed": 3,
    "indeterminate_hc3_fragile": 2,
    "supported_internal_stage_association": 2,
}
EXPECTED_TIER_COUNTS = {
    "not_selected_q_nonsignificant": 110,
    "q_significant_stability_failed": 3,
    "primary_selected_hc3_fragile": 2,
    "primary_selected_internal_robust": 2,
}


class ContractError(RuntimeError):
    """Raised when the independent semantic audit fails."""


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None:
            raise ContractError(f"missing TSV header: {path}")
        return [dict(row) for row in reader]


def atomic_write_tsv(
    path: Path,
    rows: Iterable[Mapping[str, object]],
    fieldnames: Sequence[str],
) -> None:
    materialized = list(rows)
    with tempfile.NamedTemporaryFile(
        mode="w",
        encoding="utf-8",
        newline="",
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(materialized)
    os.replace(temporary, path)


def parse_bool(value: str, field: str) -> bool:
    if value in {"TRUE", "true"}:
        return True
    if value in {"FALSE", "false"}:
        return False
    raise ContractError(f"{field} is not boolean: {value!r}")


def parse_float(value: str, field: str) -> float:
    try:
        result = float(value)
    except ValueError as error:
        raise ContractError(f"{field} is not numeric: {value!r}") from error
    if not math.isfinite(result):
        raise ContractError(f"{field} is not finite: {value!r}")
    return result


def expected_state(row: Mapping[str, str]) -> tuple[str, str]:
    estimable = parse_bool(row["primary_estimable"], "primary_estimable")
    if not estimable:
        return "untestable", "not_selected_untestable"
    q_value = parse_float(row["primary_qvalue"], "primary_qvalue")
    selected = parse_bool(row["primary_selected"], "primary_selected")
    unstable = parse_bool(row["selected_unstable"], "selected_unstable")
    fragile = parse_bool(row["selected_hc3_fragile"], "selected_hc3_fragile")
    robust = parse_bool(row["robust_display"], "robust_display")
    external = parse_bool(row["external_test_eligible"], "external_test_eligible")
    stability = parse_float(row["stability_median"], "stability_median")
    direction_agree = parse_bool(
        row["scoring_direction_agree"], "scoring_direction_agree"
    )
    hc3_supported = parse_bool(row["hc3_supported"], "hc3_supported")

    independently_selected = q_value < 0.05 and stability >= 0.50
    independently_robust = independently_selected and direction_agree and hc3_supported
    if selected != independently_selected:
        raise ContractError(f"primary selection does not rederive: {row['program_uid']}")
    if unstable != (q_value < 0.05 and stability < 0.50):
        raise ContractError(f"stability-failure flag does not rederive: {row['program_uid']}")
    if robust != independently_robust:
        raise ContractError(f"robust-display flag does not rederive: {row['program_uid']}")
    if external != (independently_robust and row["cell_type"] == "hepatocytes"):
        raise ContractError(f"external-test flag does not rederive: {row['program_uid']}")

    if q_value >= 0.05:
        return "indeterminate_nonconfirmatory", "not_selected_q_nonsignificant"
    if unstable:
        return "indeterminate_stability_failed", "q_significant_stability_failed"
    if robust:
        return (
            "supported_internal_stage_association",
            "primary_selected_internal_robust",
        )
    if fragile:
        return "indeterminate_hc3_fragile", "primary_selected_hc3_fragile"
    raise ContractError(f"unclassified q-significant row: {row['program_uid']}")


def add_check(
    checks: list[dict[str, str]],
    check_id: str,
    observed: object,
    expected: object,
    note: str,
) -> None:
    if observed != expected:
        raise ContractError(f"{check_id}: observed {observed!r}; expected {expected!r}")
    checks.append(
        {
            "check_id": check_id,
            "status": "pass",
            "observed": str(observed),
            "expected": str(expected),
            "note": note,
        }
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument(
        "--candidate-root", type=Path, default=PROJECT_ROOT / CANDIDATE_REL
    )
    args = parser.parse_args()

    project_root = args.project_root.resolve()
    candidate = args.candidate_root.resolve()
    expected_candidate = (project_root / CANDIDATE_REL).resolve()
    if project_root != PROJECT_ROOT.resolve() or candidate != expected_candidate:
        raise ContractError("project or candidate root drift")

    registry_path = candidate / "program_registry_v2.tsv"
    adjudication_path = candidate / "program_registry_v2_semantic_adjudication.tsv"
    manifest_path = candidate / "semantic_adjudication_manifest.tsv"
    registry_sha_before = sha256_file(registry_path)
    if registry_sha_before != EXPECTED_REGISTRY_SHA256:
        raise ContractError("frozen registry differs from its sealed SHA256")

    registry = read_tsv(registry_path)
    adjudication = read_tsv(adjudication_path)
    if len(registry) != 117 or len(adjudication) != 117:
        raise ContractError("registry/adjudication row count is not 117")
    registry_by_uid = {row["program_uid"]: row for row in registry}
    adjudication_by_uid = {row["program_uid"]: row for row in adjudication}
    if len(registry_by_uid) != 117 or set(registry_by_uid) != set(adjudication_by_uid):
        raise ContractError("registry/adjudication program_uid set mismatch")

    independently_expected: dict[str, tuple[str, str]] = {}
    for program_uid, row in registry_by_uid.items():
        state, tier = expected_state(row)
        independently_expected[program_uid] = state, tier
        observed = adjudication_by_uid[program_uid]
        if observed["registry_sha256"] != registry_sha_before:
            raise ContractError(f"unanchored adjudication row: {program_uid}")
        if observed["membership_sha256"] != row["membership_sha256"]:
            raise ContractError(f"membership hash mismatch: {program_uid}")
        if observed["adjudicated_state"] != state or observed["selection_tier"] != tier:
            raise ContractError(f"semantic rederivation mismatch: {program_uid}")
        if parse_bool(observed["tested_negative_authorized"], "tested_negative_authorized"):
            raise ContractError(f"negative evidence was authorized: {program_uid}")
        legacy_negative = parse_bool(row["tested_negative"], "tested_negative")
        if parse_bool(observed["legacy_tested_negative"], "legacy_tested_negative") != legacy_negative:
            raise ContractError(f"legacy negative flag mismatch: {program_uid}")
        if parse_bool(
            observed["legacy_negative_semantics_deprecated"],
            "legacy_negative_semantics_deprecated",
        ) != legacy_negative:
            raise ContractError(f"legacy deprecation flag mismatch: {program_uid}")

    state_counts = Counter(row["adjudicated_state"] for row in adjudication)
    tier_counts = Counter(row["selection_tier"] for row in adjudication)
    q_significant = sum(
        parse_float(row["primary_qvalue"], "primary_qvalue") < 0.05
        for row in registry
    )
    primary_selected = sum(
        parse_bool(row["primary_selected"], "primary_selected") for row in registry
    )
    robust_display = sum(
        parse_bool(row["robust_display"], "robust_display")
        for row in registry
    )
    legacy_negative = sum(
        parse_bool(row["tested_negative"], "tested_negative") for row in registry
    )

    checks: list[dict[str, str]] = []
    add_check(
        checks,
        "frozen_registry_sha256",
        registry_sha_before,
        EXPECTED_REGISTRY_SHA256,
        "the outcome-blind frozen registry is byte-identical",
    )
    add_check(checks, "program_uid_family", len(registry_by_uid), 117, "complete family")
    add_check(checks, "primary_q_significant", q_significant, 7, "117-family OLS BH")
    add_check(
        checks,
        "q_significant_stability_failed",
        tier_counts["q_significant_stability_failed"],
        3,
        "q-significant rows failing the frozen stability gate remain unselected",
    )
    add_check(
        checks,
        "primary_selected",
        primary_selected,
        4,
        "selection is unchanged by semantic adjudication",
    )
    add_check(
        checks,
        "robust_display",
        robust_display,
        2,
        "the sealed internally robust display subset is unchanged",
    )
    add_check(
        checks,
        "legacy_q_fail_negative_flags",
        legacy_negative,
        110,
        "historical producer flags are preserved but semantically deprecated",
    )
    add_check(
        checks,
        "adjudicated_state_counts",
        dict(state_counts),
        EXPECTED_STATE_COUNTS,
        "q failure is indeterminate, not tested-negative",
    )
    add_check(
        checks,
        "selection_tier_counts",
        dict(tier_counts),
        EXPECTED_TIER_COUNTS,
        "all seven q-significant rows and four selected rows are explicitly separated",
    )
    add_check(
        checks,
        "tested_negative_authorizations",
        sum(
            parse_bool(row["tested_negative_authorized"], "tested_negative_authorized")
            for row in adjudication
        ),
        0,
        "no equivalence or minimum-informative-effect gate exists in Plan 20",
    )

    manifest = read_tsv(manifest_path)
    if not manifest:
        raise ContractError("semantic manifest is empty")
    for row in manifest:
        path = project_root / row["relative_path"]
        if not path.is_file():
            raise ContractError(f"manifested artifact is missing: {path}")
        if sha256_file(path) != row["sha256"] or path.stat().st_size != int(row["bytes"]):
            raise ContractError(f"manifest hash/size mismatch: {path}")
    add_check(
        checks,
        "semantic_manifest_entries",
        len(manifest),
        6,
        "frozen sources, derived table, code, validator, and review-only wrapper are pinned",
    )

    registry_sha_after = sha256_file(registry_path)
    add_check(
        checks,
        "frozen_registry_postvalidation_sha256",
        registry_sha_after,
        registry_sha_before,
        "validation performed no in-place registry mutation",
    )

    report_path = candidate / "semantic_adjudication_validation.tsv"
    atomic_write_tsv(report_path, checks, list(checks[0]))
    ready_path = candidate / "SEMANTIC_ADJUDICATION_READY"
    ready_rows = [
        {
            "release_id": RELEASE_ID,
            "status": "semantic_adjudication_validated",
            "registry_sha256": registry_sha_after,
            "adjudication_sha256": sha256_file(adjudication_path),
            "manifest_sha256": sha256_file(manifest_path),
            "validation_sha256": sha256_file(report_path),
            "producer_sha256": sha256_file(
                SCRIPT_DIR / "515_adjudicate_hotspot_v2_semantics.py"
            ),
            "validator_sha256": sha256_file(Path(__file__).resolve()),
            "n_programs": 117,
            "n_q_significant": 7,
            "n_q_significant_stability_failed": 3,
            "n_primary_selected": 4,
            "n_robust_display": 2,
            "n_indeterminate_nonconfirmatory": 110,
            "n_tested_negative_authorized": 0,
            "validated_at_utc": datetime.now(timezone.utc)
            .replace(microsecond=0)
            .isoformat(),
        }
    ]
    atomic_write_tsv(ready_path, ready_rows, list(ready_rows[0]))
    print(
        "PASS: Hotspot v2 semantic adjudication independently validated; "
        "frozen registry unchanged"
    )


if __name__ == "__main__":
    main()
