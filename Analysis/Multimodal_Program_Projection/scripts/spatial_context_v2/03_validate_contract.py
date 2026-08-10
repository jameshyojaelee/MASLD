#!/usr/bin/env python3
"""Validate the common adapter contract and run negative synthetic tests."""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from contract_lib import (
    ADAPTER_REGISTRY_COLUMNS,
    ALLOWED_EVIDENCE_STATES,
    ALLOWED_INTERVAL_TYPES,
    ContractError,
    RELEASE_ID,
    SCHEMAS,
    Check,
    read_tsv,
    sha256_file,
    validate_candidate_contract,
    validate_hotspot_seal,
    write_tsv,
)


REPORT_COLUMNS = ("scope", "check_id", "status", "detail")
SELFTEST_COLUMNS = ("test_id", "expected", "observed", "status", "detail")


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def validate_metadata(candidate_root: Path) -> list[Check]:
    _, schema_rows = read_tsv(
        candidate_root / "contract_schema.tsv",
        ("artifact", "column", "required", "logical_type", "semantics"),
    )
    expected = {(artifact, column) for artifact, columns in SCHEMAS.items() for column in columns}
    observed = {(row["artifact"], row["column"]) for row in schema_rows}
    if observed != expected or len(schema_rows) != len(expected):
        raise ContractError("contract_schema.tsv does not exactly encode the required adapter columns")
    if any(row["required"] != "TRUE" for row in schema_rows):
        raise ContractError("every base-contract column must be marked required")

    _, state_rows = read_tsv(
        candidate_root / "allowed_evidence_states.tsv",
        ("evidence_state", "ordered", "numeric_encoding_permitted", "semantics"),
    )
    if {row["evidence_state"] for row in state_rows} != ALLOWED_EVIDENCE_STATES:
        raise ContractError("allowed_evidence_states.tsv vocabulary drift")
    for row in state_rows:
        if row["ordered"] != "FALSE" or row["numeric_encoding_permitted"] != "FALSE":
            raise ContractError("evidence states must remain unordered, nonnumeric categories")

    _, interval_rows = read_tsv(
        candidate_root / "allowed_interval_types.tsv",
        ("interval_type", "requires_paired_bounds", "implies_confidence_coverage", "semantics"),
    )
    if {row["interval_type"] for row in interval_rows} != ALLOWED_INTERVAL_TYPES:
        raise ContractError("allowed_interval_types.tsv vocabulary drift")
    interval_by_type = {row["interval_type"]: row for row in interval_rows}
    for interval_type, row in interval_by_type.items():
        expected_bounds = interval_type in {"confidence_interval", "donor_effect_range"}
        if row["requires_paired_bounds"] != ("TRUE" if expected_bounds else "FALSE"):
            raise ContractError(f"interval bound requirement drift for {interval_type}")
    if interval_by_type["donor_effect_range"]["implies_confidence_coverage"] != "FALSE":
        raise ContractError("donor_effect_range must never imply confidence coverage")
    if interval_by_type["confidence_interval"]["implies_confidence_coverage"] != "TRUE":
        raise ContractError("confidence_interval metadata must declare confidence coverage")

    _, adapters = read_tsv(candidate_root / "adapter_registry.tsv", ADAPTER_REGISTRY_COLUMNS)
    units = {row["effect_unit"] for row in adapters if row["gate_expectation"] in {"pass", "valid_null"}}
    if len(units) < 2:
        raise ContractError("synthetic fixture must exercise at least two distinct assay-native effect units")
    return [
        Check("contract_schema", "pass", f"{len(expected)} required artifact columns"),
        Check("state_vocabulary", "pass", f"{len(ALLOWED_EVIDENCE_STATES)} unordered categorical states"),
        Check("interval_vocabulary", "pass", f"{len(ALLOWED_INTERVAL_TYPES)} explicit interval types"),
        Check("distinct_native_units", "pass", f"fixture exercises {len(units)} non-comparable native units"),
    ]


def copy_fixture_candidate(source: Path, target: Path) -> None:
    target.mkdir(parents=True)
    shutil.copy2(source / "adapter_registry.tsv", target / "adapter_registry.tsv")
    shutil.copytree(source / "synthetic_fixture", target / "synthetic_fixture")


def mutate_tsv(path: Path, mutate) -> None:
    header, rows = read_tsv(path)
    mutate(header, rows)
    write_tsv(path, header, rows)


def reseal_adapter(adapter_root: Path, changed_relative: str) -> None:
    manifest_path = adapter_root / "execution_manifest.tsv"
    header, rows = read_tsv(manifest_path, SCHEMAS["execution_manifest.tsv"])
    found = False
    changed_path = adapter_root / changed_relative
    for row in rows:
        if row["relative_path"] == changed_relative:
            row["bytes"] = str(changed_path.stat().st_size)
            row["sha256"] = sha256_file(changed_path)
            found = True
    if not found:
        raise RuntimeError(f"execution manifest lacks {changed_relative}")
    write_tsv(manifest_path, header, rows)
    gate_path = adapter_root / "gate_status.tsv"
    gate_header, gate_rows = read_tsv(gate_path, SCHEMAS["gate_status.tsv"])
    gate_rows[0]["execution_manifest_sha256"] = sha256_file(manifest_path)
    write_tsv(gate_path, gate_header, gate_rows)


def expect_rejection(
    test_id: str,
    candidate_root: Path,
    hotspot_root: Path,
    project_root: Path,
    mutator,
) -> dict[str, str]:
    with tempfile.TemporaryDirectory(prefix="spatial_contract_") as tmp:
        copied = Path(tmp) / "candidate"
        copy_fixture_candidate(candidate_root, copied)
        mutator(copied)
        try:
            validate_candidate_contract(copied, copied / "adapter_registry.tsv", hotspot_root, project_root)
        except ContractError as exc:
            return {
                "test_id": test_id,
                "expected": "reject",
                "observed": "reject",
                "status": "pass",
                "detail": str(exc),
            }
        return {
            "test_id": test_id,
            "expected": "reject",
            "observed": "accept",
            "status": "fail",
            "detail": "invalid fixture was accepted",
        }


def run_selftests(candidate_root: Path, hotspot_root: Path, project_root: Path) -> list[dict[str, str]]:
    results = []
    try:
        validate_candidate_contract(candidate_root, candidate_root / "adapter_registry.tsv", hotspot_root, project_root)
        results.append(
            {
                "test_id": "valid_contract_accepts",
                "expected": "accept",
                "observed": "accept",
                "status": "pass",
                "detail": "all valid synthetic adapters accepted",
            }
        )
    except ContractError as exc:
        results.append(
            {
                "test_id": "valid_contract_accepts",
                "expected": "accept",
                "observed": "reject",
                "status": "fail",
                "detail": str(exc),
            }
        )

    def mismatch_count(root: Path) -> None:
        adapter = root / "synthetic_fixture/valid/synthetic_visium_pass"

        def mutate(_header, rows):
            rows[0]["n_biological"] = "99"

        mutate_tsv(adapter / "program_effects.tsv", mutate)
        reseal_adapter(adapter, "program_effects.tsv")

    results.append(expect_rejection("biological_count_mismatch_rejected", candidate_root, hotspot_root, project_root, mismatch_count))

    def universal_score(root: Path) -> None:
        path = root / "synthetic_fixture/valid/synthetic_visium_pass/program_effects.tsv"

        def mutate(header, rows):
            header.append("universal_score")
            for row in rows:
                row["universal_score"] = "1"

        mutate_tsv(path, mutate)

    results.append(expect_rejection("universal_score_column_rejected", candidate_root, hotspot_root, project_root, universal_score))

    def invalid_state(root: Path) -> None:
        adapter = root / "synthetic_fixture/valid/synthetic_visium_pass"

        def mutate(_header, rows):
            rows[0]["evidence_state"] = "supported_plus"

        mutate_tsv(adapter / "program_effects.tsv", mutate)
        reseal_adapter(adapter, "program_effects.tsv")

    results.append(expect_rejection("invalid_evidence_state_rejected", candidate_root, hotspot_root, project_root, invalid_state))

    def negative_without_rule(root: Path) -> None:
        adapter = root / "synthetic_fixture/valid/synthetic_visium_pass"

        def mutate_effect(_header, rows):
            rows[1]["evidence_state"] = "tested_negative"
            rows[1]["negative_call_rule_id"] = ""

        def mutate_testability(_header, rows):
            rows[1]["evidence_state"] = "tested_negative"

        mutate_tsv(adapter / "program_effects.tsv", mutate_effect)
        mutate_tsv(adapter / "program_testability.tsv", mutate_testability)
        reseal_adapter(adapter, "program_effects.tsv")
        reseal_adapter(adapter, "program_testability.tsv")

    results.append(
        expect_rejection(
            "tested_negative_without_adequate_negative_rule_rejected",
            candidate_root,
            hotspot_root,
            project_root,
            negative_without_rule,
        )
    )

    def double_labeled_dispersion(root: Path) -> None:
        adapter = root / "synthetic_fixture/valid/synthetic_visium_pass"

        def mutate(_header, rows):
            rows[0]["matched_null_sd"] = rows[0]["std_error"]

        mutate_tsv(adapter / "program_effects.tsv", mutate)
        reseal_adapter(adapter, "program_effects.tsv")

    results.append(
        expect_rejection(
            "sampling_se_and_matched_null_sd_double_label_rejected",
            candidate_root,
            hotspot_root,
            project_root,
            double_labeled_dispersion,
        )
    )

    def registry_drift(root: Path) -> None:
        adapter = root / "synthetic_fixture/valid/synthetic_visium_pass"

        def mutate(_header, rows):
            rows[0]["registry_sha256"] = "0" * 64

        mutate_tsv(adapter / "program_effects.tsv", mutate)
        reseal_adapter(adapter, "program_effects.tsv")

    results.append(expect_rejection("registry_hash_drift_rejected", candidate_root, hotspot_root, project_root, registry_drift))

    def missing_uncertainty(root: Path) -> None:
        adapter = root / "synthetic_fixture/valid/synthetic_visium_pass"

        def mutate(_header, rows):
            rows[0]["std_error"] = ""
            rows[0]["interval_low"] = ""
            rows[0]["interval_high"] = ""
            rows[0]["interval_type"] = "none"

        mutate_tsv(adapter / "program_effects.tsv", mutate)
        reseal_adapter(adapter, "program_effects.tsv")

    results.append(expect_rejection("missing_se_and_interval_rejected", candidate_root, hotspot_root, project_root, missing_uncertainty))

    def half_interval(root: Path) -> None:
        adapter = root / "synthetic_fixture/valid/synthetic_visium_pass"

        def mutate(_header, rows):
            rows[0]["interval_high"] = ""

        mutate_tsv(adapter / "program_effects.tsv", mutate)
        reseal_adapter(adapter, "program_effects.tsv")

    results.append(expect_rejection("half_interval_rejected", candidate_root, hotspot_root, project_root, half_interval))

    def invalid_interval_type(root: Path) -> None:
        adapter = root / "synthetic_fixture/valid/synthetic_visium_pass"

        def mutate(_header, rows):
            rows[0]["interval_type"] = "credible_interval"

        mutate_tsv(adapter / "program_effects.tsv", mutate)
        reseal_adapter(adapter, "program_effects.tsv")

    results.append(expect_rejection("invalid_interval_type_rejected", candidate_root, hotspot_root, project_root, invalid_interval_type))

    with tempfile.TemporaryDirectory(prefix="spatial_ready_") as tmp:
        copied_hotspot = Path(tmp) / "hotspot"
        copied_hotspot.mkdir()
        for filename in (
            "READY",
            "validation_status.tsv",
            "release_manifest.tsv",
            "program_registry_v2.tsv",
            "program_membership_v2.tsv",
        ):
            shutil.copy2(hotspot_root / filename, copied_hotspot / filename)

        def mutate_ready(_header, rows):
            rows[0]["registry_sha256"] = "f" * 64

        mutate_tsv(copied_hotspot / "READY", mutate_ready)
        try:
            validate_hotspot_seal(copied_hotspot)
        except ContractError as exc:
            results.append(
                {
                    "test_id": "ready_registry_linkage_rejected",
                    "expected": "reject",
                    "observed": "reject",
                    "status": "pass",
                    "detail": str(exc),
                }
            )
        else:
            results.append(
                {
                    "test_id": "ready_registry_linkage_rejected",
                    "expected": "reject",
                    "observed": "accept",
                    "status": "fail",
                    "detail": "corrupted READY registry hash was accepted",
                }
            )
    return results


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--candidate-root", type=Path, default=None)
    parser.add_argument("--hotspot-root", type=Path, default=None)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    candidate_root = args.candidate_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "spatial_context"
    )
    hotspot_root = args.hotspot_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "hotspot"
    )
    try:
        checks = validate_metadata(candidate_root)
        checks.extend(
            validate_candidate_contract(
                candidate_root,
                candidate_root / "adapter_registry.tsv",
                hotspot_root,
                project_root,
            )
        )
        report_rows = [
            {"scope": "contract_fixture", "check_id": row.check_id, "status": row.status, "detail": row.detail}
            for row in checks
        ]
        write_tsv(candidate_root / "validation_report.tsv", REPORT_COLUMNS, report_rows)
        selftests = run_selftests(candidate_root, hotspot_root, project_root)
        write_tsv(candidate_root / "synthetic_fixture/selftest_report.tsv", SELFTEST_COLUMNS, selftests)
        if any(row["status"] != "pass" for row in selftests):
            raise ContractError("one or more synthetic negative tests failed")
        write_tsv(
            candidate_root / "sp_int_01_02_status.tsv",
            ("release_id", "status", "real_outcomes_imported", "plan13_complete", "validated_at_utc"),
            [
                {
                    "release_id": RELEASE_ID,
                    "status": "sp_int_01_02_complete_contract_only",
                    "real_outcomes_imported": False,
                    "plan13_complete": False,
                    "validated_at_utc": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                }
            ],
        )
        print(f"PASS: {len(checks)} contract checks and {len(selftests)} synthetic acceptance/rejection tests")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
