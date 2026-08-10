#!/usr/bin/env python3
"""Freeze the terminal execution after the biopsy mapping v2 remediation."""

from __future__ import annotations

import datetime as dt
import json

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_text, sha256_file, stable_json_sha256, write_tsv


def main() -> None:
    destination = CANDIDATE_ROOT / "FINAL_EXECUTION_V3_READY.json"
    if destination.exists():
        raise RuntimeError(f"Final v3 already frozen: {destination}")
    prerequisite = CANDIDATE_ROOT / "BIOPSY_MAPPING_V2_READY.json"
    prior = CANDIDATE_ROOT / "FINAL_EXECUTION_CODE_AMENDMENT_01_READY.json"
    if not prerequisite.is_file() or not prior.is_file():
        raise RuntimeError("Biopsy mapping v2 and prior final provenance freeze are required")
    files = [
        "functional_analysis_common.R", "04_run_bulk_assays.R", "07_run_evidence_classes.R",
        "11_run_null_permutations.py", "12_build_final_integration.py", "13_validate_final.py",
        "16_capture_environment.py", "run_final.sbatch",
    ]
    rows = []
    for name in files:
        path = SCRIPT_ROOT / name
        rows.append({"path": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size})
    contract = {
        "freeze_id": "public-functional-terminal-v3",
        "frozen_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "supersedes_sha256": sha256_file(prior),
        "biopsy_mapping_v2_sha256": sha256_file(prerequisite),
        "files": rows,
        "scope": "complete 117-family native effects, biopsy source remediation, evidence-class tests, nulls, source/environment/execution manifests, mechanical verdict, and validation",
        "scientific_change_from_v2": "none; only biopsy mapping implementation changed after a reproduced software incompatibility",
        "biological_positivity_required_for_software_acceptance": False,
    }
    contract["code_bundle_sha256"] = stable_json_sha256(contract)
    write_tsv(CANDIDATE_ROOT / "final_v3_code_manifest.tsv", rows, ["path", "sha256", "size_bytes"])
    atomic_write_text(destination, json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
