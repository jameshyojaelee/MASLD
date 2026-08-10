#!/usr/bin/env python3
"""Seal the log-excluding terminal validator before revalidation."""

from __future__ import annotations

import datetime as dt
import json

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_text, sha256_file, stable_json_sha256, write_tsv


EXPECTED_AMENDMENT_SHA256 = "8a342aa7d12ba12babcc8e4cd0f3fdaabdc44a10b150e167b49e7a0b070d6624"


def main() -> None:
    destination = CANDIDATE_ROOT / "TERMINAL_REVALIDATION_READY.json"
    if destination.exists():
        raise RuntimeError(f"Terminal revalidation is already frozen: {destination}")
    amendment = CANDIDATE_ROOT / "RELEASE_MANIFEST_LOG_AMENDMENT_01.json"
    if not amendment.is_file() or sha256_file(amendment) != EXPECTED_AMENDMENT_SHA256:
        raise RuntimeError("Release-manifest log amendment is absent or has drifted")
    files = ["13_validate_final.py", "23_freeze_terminal_revalidation.py", "run_terminal_revalidation.sbatch"]
    rows = []
    for name in files:
        path = SCRIPT_ROOT / name
        rows.append({"path": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size})
    contract = {
        "freeze_id": "public-functional-terminal-log-exclusion-v1",
        "frozen_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "release_manifest_log_amendment_sha256": EXPECTED_AMENDMENT_SHA256,
        "files": rows,
        "execution_scope": "terminal validation and manifest regeneration only",
        "scientific_analysis_rerun": False,
        "scientific_contract_change": "none",
    }
    contract["code_bundle_sha256"] = stable_json_sha256(contract)
    write_tsv(CANDIDATE_ROOT / "terminal_revalidation_code_manifest.tsv", rows, ["path", "sha256", "size_bytes"])
    atomic_write_text(destination, json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
