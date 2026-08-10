#!/usr/bin/env python3
"""Freeze the complete-family output correction before rerunning projections."""

from __future__ import annotations

import datetime as dt
import json

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_text, sha256_file, stable_json_sha256, write_tsv


def main() -> None:
    destination = CANDIDATE_ROOT / "PROGRAM_COMPLETENESS_AMENDMENT_READY.json"
    if destination.exists():
        raise RuntimeError(f"Completeness amendment already exists: {destination}")
    prior = CANDIDATE_ROOT / "PROGRAM_EXECUTION_CODE_READY.json"
    biopsy = CANDIDATE_ROOT / "BIOPSY_EXECUTION_AMENDMENT_READY.json"
    if not prior.is_file() or not biopsy.is_file():
        raise RuntimeError("Prior execution and biopsy-remediation freezes are required")
    code = [
        SCRIPT_ROOT / "functional_analysis_common.R", SCRIPT_ROOT / "04_run_bulk_assays.R",
        SCRIPT_ROOT / "05_build_schlo_pseudobulk.py", SCRIPT_ROOT / "06_run_schlo_programs.R",
        SCRIPT_ROOT / "13_validate_final.py",
    ]
    rows = [{"path": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size} for path in code]
    contract = {
        "freeze_id": "public-functional-complete-117-family-v1",
        "frozen_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "prior_program_freeze_sha256": sha256_file(prior),
        "biopsy_remediation_freeze_sha256": sha256_file(biopsy),
        "files": rows,
        "reason": "Pre-adjudication schema audit found that untestable programs were omitted from effect rows, silently shrinking output and BH families.",
        "correction": "emit exactly 117 rows per dataset/contrast/scoring/lineage family; label missing statistics untestable; use fixed BH n=117 and confirmatory n=2; no model, contrast, direction, weight, or selection change",
    }
    contract["code_bundle_sha256"] = stable_json_sha256(contract)
    write_tsv(CANDIDATE_ROOT / "program_completeness_code_manifest.tsv", rows, ["path", "sha256", "size_bytes"])
    atomic_write_text(destination, json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
