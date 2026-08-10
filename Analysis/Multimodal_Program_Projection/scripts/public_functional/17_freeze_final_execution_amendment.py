#!/usr/bin/env python3
"""Freeze provenance-complete final execution after the pre-run manifest audit."""

from __future__ import annotations

import datetime as dt
import json

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_text, sha256_file, stable_json_sha256, write_tsv


def main() -> None:
    destination = CANDIDATE_ROOT / "FINAL_EXECUTION_CODE_AMENDMENT_01_READY.json"
    if destination.exists():
        raise RuntimeError(f"Final amendment already frozen: {destination}")
    prior = CANDIDATE_ROOT / "FINAL_EXECUTION_CODE_READY.json"
    if not prior.is_file():
        raise FileNotFoundError(prior)
    files = ["11_run_null_permutations.py", "12_build_final_integration.py", "13_validate_final.py", "16_capture_environment.py", "run_final.sbatch"]
    rows = []
    for name in files:
        path = SCRIPT_ROOT / name
        rows.append({"path": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size})
    contract = {
        "freeze_id": "public-functional-final-execution-provenance-v2",
        "frozen_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "supersedes_sha256": sha256_file(prior),
        "files": rows,
        "reason": "Pre-run provenance audit required same-study OS-H5AD and GPL16686 annotation amendment manifests plus explicit Python/R/execution identities in the root release contract.",
        "scientific_change": "none",
        "biological_positivity_required_for_software_acceptance": False,
    }
    contract["code_bundle_sha256"] = stable_json_sha256(contract)
    write_tsv(CANDIDATE_ROOT / "final_execution_amendment_manifest.tsv", rows, ["path", "sha256", "size_bytes"])
    atomic_write_text(destination, json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
