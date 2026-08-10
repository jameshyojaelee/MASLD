#!/usr/bin/env python3
"""Freeze null, integration, validation, and current evidence-class code before final jobs."""

from __future__ import annotations

import datetime as dt
import json

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_text, sha256_file, stable_json_sha256, write_tsv


def main() -> None:
    destination = CANDIDATE_ROOT / "FINAL_EXECUTION_CODE_READY.json"
    if destination.exists():
        raise RuntimeError(f"Final execution already frozen: {destination}")
    prerequisites = [
        "SEALED.json", "EVIDENCE_CLASS_CODE_READY.json", "BIOPSY_EXECUTION_AMENDMENT_READY.json",
        "PROGRAM_COMPLETENESS_AMENDMENT_READY.json",
    ]
    for name in prerequisites:
        if not (CANDIDATE_ROOT / name).is_file():
            raise FileNotFoundError(CANDIDATE_ROOT / name)
    files = [
        "public_functional_common.py", "functional_analysis_common.R", "07_run_evidence_classes.R",
        "11_run_null_permutations.py", "12_build_final_integration.py", "13_validate_final.py",
        "run_evidence_classes.sbatch", "run_null_permutations.sbatch", "run_final.sbatch",
    ]
    rows = []
    for name in files:
        path = SCRIPT_ROOT / name
        rows.append({"path": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size})
    contract = {
        "freeze_id": "public-functional-final-execution-v1",
        "frozen_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "prerequisite_hashes": {name: sha256_file(CANDIDATE_ROOT / name) for name in prerequisites},
        "files": rows,
        "scope": "evidence-class execution after biopsy annotation remediation; exact biological-unit nulls; complete-family collation; mechanical six-gate verdict; independent validation",
        "biological_positivity_required_for_software_acceptance": False,
    }
    contract["code_bundle_sha256"] = stable_json_sha256(contract)
    write_tsv(CANDIDATE_ROOT / "final_execution_code_manifest.tsv", rows, ["path", "sha256", "size_bytes"])
    atomic_write_text(destination, json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
