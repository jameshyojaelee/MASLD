#!/usr/bin/env python3
"""Independently validate the GEN terminal closure and identifier units."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from pathlib import Path

from genetics_common import (
    ContractError,
    DEFAULT_CANDIDATE_ROOT,
    PROJECT_ROOT,
    assert_candidate_root,
    atomic_write_tsv,
    clean,
    ensembl_base,
    read_table,
    sha256_file,
)


RELEASE_ID = "program-context-v2-candidate-2026-08-07"
EXPECTED_PREFLIGHT_SHA256 = (
    "b639664c3fc9199a54719a75b94bd60a4b11e57a303d91e65134894f98cf8aeb"
)
EXPECTED_HONG_GATE_SHA256 = (
    "1e9a2a48b0a6488008267c06a3cbf47861eb1b0e51c28a571f6a406fc60c2d25"
)
EXPECTED_HONG_VALIDATION_SHA256 = (
    "fe3856c9fa73166cde55a0320b8f8e1cb5d5d83897fbf0be8f98b01925cf49e8"
)
EXPECTED_MISSING_ENSGS = {"ENSG00000137411", "ENSG00000204610"}


def fail(condition: bool, message: str) -> None:
    if condition:
        raise ContractError(message)


def add_check(
    checks: list[dict[str, str]],
    check_id: str,
    observed: object,
    expected: object,
    note: str,
) -> None:
    fail(observed != expected, f"{check_id}: observed {observed!r}; expected {expected!r}")
    checks.append(
        {
            "check_id": check_id,
            "status": "pass",
            "observed": str(observed),
            "expected": str(expected),
            "claim_authorized": "false",
            "note": note,
        }
    )


def source_primary_ensgs(path: Path) -> set[str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fail(
            reader.fieldnames is None
            or not {"Ensembl", "Signal"}.issubset(reader.fieldnames),
            "source eQTL snapshot schema drift",
        )
        identifiers = {
            ensembl_base(row["Ensembl"])
            for row in reader
            if clean(row["Signal"]).startswith("1:")
        }
    fail("" in identifiers, "blank primary source Ensembl identifier")
    return identifiers


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--candidate-root", type=Path, default=DEFAULT_CANDIDATE_ROOT)
    args = parser.parse_args()

    fail(args.project_root.resolve() != PROJECT_ROOT.resolve(), "project-root drift")
    project_root = args.project_root.resolve()
    candidate = assert_candidate_root(project_root, args.candidate_root)

    closure_path = candidate / "terminal_closure.tsv"
    manifest_path = candidate / "terminal_closure_manifest.tsv"
    preflight_path = candidate / "gate_status.tsv"
    hong_gate_path = candidate / "hong_source_gate_status.tsv"
    hong_validation_path = candidate / "hong_validation_report.tsv"
    for path in (
        closure_path,
        manifest_path,
        preflight_path,
        hong_gate_path,
        hong_validation_path,
    ):
        fail(not path.is_file(), f"required closure artifact missing: {path}")

    checks: list[dict[str, str]] = []
    add_check(
        checks,
        "historical_preflight_sha256",
        sha256_file(preflight_path),
        EXPECTED_PREFLIGHT_SHA256,
        "historical GEN00-02 preflight is preserved rather than overwritten",
    )
    add_check(
        checks,
        "terminal_hong_gate_sha256",
        sha256_file(hong_gate_path),
        EXPECTED_HONG_GATE_SHA256,
        "terminal GEN03 authority is pinned",
    )
    add_check(
        checks,
        "terminal_hong_validation_sha256",
        sha256_file(hong_validation_path),
        EXPECTED_HONG_VALIDATION_SHA256,
        "independent GEN03 terminal validation is pinned",
    )

    closure_rows = read_table(closure_path, "tsv")
    fail(len(closure_rows) != 1, "terminal closure must contain exactly one row")
    closure = closure_rows[0]
    expected_fields = {
        "release_id": RELEASE_ID,
        "closure_id": "GEN_TERMINAL_CLOSURE_V1",
        "terminal_status": "coverage_limited_terminal",
        "downstream_selection_authority": "terminal_closure",
        "authoritative_gate_path": "hong_source_gate_status.tsv",
        "authoritative_gate_sha256": EXPECTED_HONG_GATE_SHA256,
        "independent_validation_path": "hong_validation_report.tsv",
        "independent_validation_sha256": EXPECTED_HONG_VALIDATION_SHA256,
        "historical_preflight_path": "gate_status.tsv",
        "historical_preflight_sha256": EXPECTED_PREFLIGHT_SHA256,
        "superseded_preflight_gate_ids": "GEN03_HONG_SOURCE_AUDIT;PREFLIGHT_OVERALL",
        "supersession_scope": "downstream_gate_selection_only",
        "historical_preflight_preserved": "true",
        "context_rescue_authorized": "false",
        "negative_claim_authorized": "false",
        "canonical_promotion_status": "not_promoted",
    }
    for field, expected in expected_fields.items():
        add_check(
            checks,
            f"closure_field_{field}",
            closure.get(field),
            expected,
            "terminal selection semantics are explicit and fail closed",
        )

    preflight = {
        row["gate_id"]: row for row in read_table(preflight_path, "tsv")
    }
    add_check(
        checks,
        "historical_gen03_stale_row",
        preflight["GEN03_HONG_SOURCE_AUDIT"]["status"],
        "not_started_outside_preflight_scope",
        "this preserved row is superseded only for downstream selection",
    )
    add_check(
        checks,
        "historical_preflight_overall_stale_row",
        preflight["PREFLIGHT_OVERALL"]["status"],
        "preflight_complete_claim_integration_blocked",
        "this preserved row is superseded only for downstream selection",
    )
    hong_gate = {
        row["gate_id"]: row for row in read_table(hong_gate_path, "tsv")
    }
    add_check(
        checks,
        "terminal_gen03_overall",
        hong_gate["GEN03_OVERALL"]["status"],
        "coverage_limited_terminal",
        "this is the downstream GEN03 selection verdict",
    )
    add_check(
        checks,
        "terminal_context_rescue_authorization",
        hong_gate["GEN03_CONTEXT_RESCUE_AUTHORIZATION"]["status"],
        "prohibited",
        "no context-rescue panel or rescued fraction is authorized",
    )
    add_check(
        checks,
        "terminal_negative_claim_authorization",
        hong_gate["GEN03_NEGATIVE_CLAIM_AUTHORIZATION"]["status"],
        "prohibited",
        "coverage limitation is not a genetic or context null",
    )

    source_path = (
        candidate
        / "work/input_snapshots/Liver_eQTL_Meta_Leads_ST3_20240530.tsv"
    )
    source_ensgs = source_primary_ensgs(source_path)
    observability = read_table(candidate / "gene_observability.tsv", "tsv")
    mapped_rows = [
        row for row in observability if clean(row["source_defined_egene"]) == "true"
    ]
    assignments = [
        ensembl_base(identifier)
        for row in mapped_rows
        for identifier in clean(row["source_eqtl_ensembl"]).split(";")
        if clean(identifier)
    ]
    represented = set(assignments)
    missing = source_ensgs - represented
    fail(not represented.issubset(source_ensgs), "mapped rows include a non-source eGene")

    source_reproduction = {
        row["metric"]: row
        for row in read_table(candidate / "source_eqtl_reproduction.tsv", "tsv")
    }
    mapping_audit = {
        row["mapping_basis"]: row
        for row in read_table(candidate / "gene_mapping_audit.tsv", "tsv")
    }
    add_check(
        checks,
        "source_wide_unique_source_defined_egenes",
        len(source_ensgs),
        6_564,
        "source-level unit: unique primary-signal Ensembl eGenes",
    )
    add_check(
        checks,
        "source_reproduction_table_egenes",
        int(source_reproduction["source_defined_egenes"]["observed"]),
        6_564,
        "published source-wide count independently agrees",
    )
    add_check(
        checks,
        "mapped_symbol_annotation_rows",
        len(mapped_rows),
        6_583,
        "atlas mapping unit: symbol annotation rows, not unique eGenes",
    )
    add_check(
        checks,
        "mapping_audit_symbol_rows",
        int(mapping_audit["ensembl"]["gene_count"]),
        6_583,
        "producer mapping audit agrees with the independent row count",
    )
    add_check(
        checks,
        "mapped_source_ensg_assignments",
        len(assignments),
        6_584,
        "one mapped symbol row carries two source ENSG assignments",
    )
    add_check(
        checks,
        "represented_unique_source_ensgs",
        len(represented),
        6_562,
        "identifier-level unit represented among mapped symbol rows",
    )
    add_check(
        checks,
        "unrepresented_source_ensgs",
        ";".join(sorted(missing)),
        ";".join(sorted(EXPECTED_MISSING_ENSGS)),
        "two source-wide eGenes have no mapped symbol annotation row",
    )
    add_check(
        checks,
        "closure_source_wide_count",
        int(closure["source_wide_unique_source_defined_egenes"]),
        6_564,
        "closure stores the source-wide unit explicitly",
    )
    add_check(
        checks,
        "closure_mapped_symbol_rows",
        int(closure["mapped_symbol_annotation_rows"]),
        6_583,
        "closure stores the mapped-row unit explicitly",
    )
    add_check(
        checks,
        "closure_represented_source_ensgs",
        int(closure["represented_unique_source_ensgs"]),
        6_562,
        "closure stores the represented-identifier unit explicitly",
    )

    manifest = read_table(manifest_path, "tsv")
    add_check(
        checks,
        "terminal_closure_manifest_entries",
        len(manifest),
        11,
        "authority, historical state, source identifiers, output, code, and wrapper are pinned",
    )
    for row in manifest:
        path = project_root / row["relative_path"]
        fail(not path.is_file(), f"manifested closure artifact missing: {path}")
        fail(
            sha256_file(path) != row["sha256"]
            or path.stat().st_size != int(row["bytes"]),
            f"closure manifest hash/size mismatch: {path}",
        )

    report_path = candidate / "terminal_closure_validation.tsv"
    atomic_write_tsv(report_path, checks, list(checks[0]))
    ready_path = candidate / "GEN_TERMINAL_CLOSURE_READY"
    ready = {
        "release_id": RELEASE_ID,
        "status": "coverage_limited_terminal_validated",
        "closure_sha256": sha256_file(closure_path),
        "terminal_gate_sha256": sha256_file(hong_gate_path),
        "terminal_validation_sha256": sha256_file(hong_validation_path),
        "historical_preflight_sha256": sha256_file(preflight_path),
        "closure_manifest_sha256": sha256_file(manifest_path),
        "closure_validation_sha256": sha256_file(report_path),
        "validator_sha256": sha256_file(Path(__file__).resolve()),
        "source_wide_unique_source_defined_egenes": 6_564,
        "mapped_symbol_annotation_rows": 6_583,
        "represented_unique_source_ensgs": 6_562,
        "context_rescue_authorized": "false",
        "negative_claim_authorized": "false",
        "canonical_promotion_status": "not_promoted",
        "validated_at_utc": datetime.now(timezone.utc)
        .replace(microsecond=0)
        .isoformat(),
    }
    atomic_write_tsv(ready_path, [ready], list(ready))
    print(
        "PASS: GEN terminal closure validated; 6564 source eGenes, 6583 mapped "
        "symbol rows, and 6562 represented source ENSGs remain distinct"
    )


if __name__ == "__main__":
    main()
