#!/usr/bin/env python3
"""Seal the completion-audit remediation code before outcome regeneration.

This freeze is method-only.  It cannot change the frozen programs, hypotheses,
contrasts, evidence classes, multiplicity families, or promotion thresholds.
"""

from __future__ import annotations

import datetime as dt
import json

from public_functional_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    SCRIPT_ROOT,
    atomic_write_text,
    sha256_file,
    stable_json_sha256,
    write_tsv,
)


EXPECTED_AMENDMENT_SHA256 = "2de5f76c6b171d517ae1889ba8c6d73ed19faab54037417f001b0e7824b51c75"


def main() -> None:
    destination = CANDIDATE_ROOT / "METHOD_COMPLIANCE_EXECUTION_READY.json"
    if destination.exists():
        raise RuntimeError(f"Method-compliance execution is already frozen: {destination}")
    amendment = CANDIDATE_ROOT / "COMPLETION_AUDIT_AMENDMENT_01.json"
    if not amendment.is_file() or sha256_file(amendment) != EXPECTED_AMENDMENT_SHA256:
        raise RuntimeError("Completion-audit amendment is absent or has drifted")

    files = [
        "03_run_source_gates.py",
        "04_run_bulk_assays.R",
        "05_build_schlo_pseudobulk.py",
        "06_run_schlo_programs.R",
        "07_run_evidence_classes.R",
        "11_run_null_permutations.py",
        "12_build_final_integration.py",
        "13_validate_final.py",
        "16_capture_environment.py",
        "20_freeze_completion_audit_amendment.py",
        "21_freeze_method_compliance_execution.py",
        "functional_analysis_common.R",
        "public_functional_common.py",
        "test_public_functional.R",
        "run_source_gates.sbatch",
        "run_bulk_assays.sbatch",
        "run_schlo.sbatch",
        "run_evidence_classes.sbatch",
        "run_tests.sbatch",
        "run_final.sbatch",
    ]
    rows = []
    for name in files:
        path = SCRIPT_ROOT / name
        if not path.is_file():
            raise RuntimeError(f"Missing method-compliance producer: {path}")
        rows.append(
            {
                "path": str(path.relative_to(PROJECT_ROOT)),
                "sha256": sha256_file(path),
                "size_bytes": path.stat().st_size,
            }
        )

    contract = {
        "freeze_id": "public-functional-completion-audit-method-compliance-v1",
        "frozen_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "completion_audit_amendment_sha256": EXPECTED_AMENDMENT_SHA256,
        "files": rows,
        "allowed_regeneration": [
            "source gates and schema audit",
            "GSE200418 program projections",
            "GSE207889 pseudobulk program projections",
            "GSE106737 program projections",
            "GSE106737 evidence-class tests",
            "terminal null, integration, environment, manifest, and validation products",
        ],
        "scientific_contract_change": "none",
        "biological_positivity_required_for_acceptance": False,
    }
    contract["code_bundle_sha256"] = stable_json_sha256(contract)
    write_tsv(
        CANDIDATE_ROOT / "method_compliance_code_manifest.tsv",
        rows,
        ["path", "sha256", "size_bytes"],
    )
    atomic_write_text(destination, json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
