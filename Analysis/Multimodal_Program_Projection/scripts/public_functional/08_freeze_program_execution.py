#!/usr/bin/env python3
"""Freeze the exact native program-projection code before result inspection."""

from __future__ import annotations

import datetime as dt
import json

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_text, require_sealed, sha256_file, stable_json_sha256, write_tsv


FILES = [
    "public_functional_common.py",
    "functional_analysis_common.R",
    "03_run_source_gates.py",
    "04_run_bulk_assays.R",
    "05_build_schlo_pseudobulk.py",
    "06_run_schlo_programs.R",
    "run_bulk_assays.sbatch",
    "run_schlo.sbatch",
]


def main() -> None:
    seal = require_sealed()
    destination = CANDIDATE_ROOT / "PROGRAM_EXECUTION_CODE_READY.json"
    if destination.exists():
        raise RuntimeError(f"Program execution code already frozen: {destination}")
    rows = []
    for name in FILES:
        path = SCRIPT_ROOT / name
        rows.append({"path": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size})
    bundle = {
        "freeze_id": "public-functional-program-execution-v1",
        "frozen_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "parent_specification_sha256": seal["specification_sha256"],
        "files": rows,
        "scope": "program scoring, assay-native models, sensitivity scoring, and source-replicate pseudobulk; evidence-class and final-integration phases freeze separately",
    }
    bundle["code_bundle_sha256"] = stable_json_sha256(bundle)
    write_tsv(CANDIDATE_ROOT / "program_execution_code_manifest.tsv", rows, ["path", "sha256", "size_bytes"])
    atomic_write_text(destination, json.dumps(bundle, indent=2, sort_keys=True) + "\n")
    print(json.dumps(bundle, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
