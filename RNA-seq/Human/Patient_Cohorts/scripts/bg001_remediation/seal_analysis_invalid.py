#!/usr/bin/env python3
"""Publish/verify an immutable INVALID terminal state for hard scientific gates."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from safe_io import publish_new_bytes, require_regular_file


SCHEMA = "bg001-analysis-invalid-v1"
COMMON_BOUND = (
    "contract/run_contract.json",
    "contract/source_manifest.tsv",
    "BASELINE_FROZEN.json",
    "manifests/baseline_frozen_files.sha256",
    "manifests/baseline_manifests.sha256",
    "RECOUNT_COMPLETE",
    "manifests/recount_artifacts.sha256",
    "comparisons/structural_validation.json",
)
R0_BOUND = (
    "comparisons/R0_reproduction.json",
    "arms/R0/results/integration/deg_results.csv",
    "arms/R0/results/integration/model_design.tsv",
    "arms/R0/provenance/arm_validation.tsv",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bound_hashes(root: Path) -> dict[str, str]:
    relatives = COMMON_BOUND + R0_BOUND
    result: dict[str, str] = {}
    for relative in relatives:
        path = root / relative
        require_regular_file(path)
        result[relative] = sha256(path)
    return dict(sorted(result.items()))


def expected_payload(root: Path) -> dict[str, object]:
    report_path = root / "comparisons/R0_reproduction.json"
    require_regular_file(report_path)
    report = json.loads(report_path.read_text())
    if report.get("status") != "FAIL" or not report.get("failures"):
        raise SystemExit("R0 INVALID state requires a nonempty FAIL reproduction report")
    return {
        "schema": SCHEMA,
        "run_id": root.name,
        "status": "INVALID",
        "hard_gate": "R0_REPRODUCTION",
        "failure_report": str(report_path.relative_to(root)),
        "failures": report["failures"],
        "bound_inputs": bound_hashes(root),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("seal", "verify"))
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--kind", required=True, choices=("r0",))
    args = parser.parse_args()
    root = args.run_root.resolve(strict=True)
    marker = root / "ANALYSIS_INVALID.json"
    complete = root / "ANALYSIS_COMPLETE.json"
    if complete.exists() or complete.is_symlink():
        raise SystemExit("INVALID and COMPLETE analysis states are mutually exclusive")
    expected = expected_payload(root)
    if args.action == "seal":
        if marker.exists() or marker.is_symlink():
            raise SystemExit("Refusing an existing ANALYSIS_INVALID marker")
        publish_new_bytes(
            marker,
            (json.dumps(expected, indent=2, sort_keys=True) + "\n").encode(),
            root,
        )
    else:
        require_regular_file(marker)
        if json.loads(marker.read_text()) != expected:
            raise SystemExit("ANALYSIS_INVALID marker or bound hard-gate inputs drifted")
    print(f"PASS analysis {args.action} INVALID: {expected['hard_gate']}")


if __name__ == "__main__":
    main()
