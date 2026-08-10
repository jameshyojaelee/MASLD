#!/usr/bin/env python3
"""Freeze the biopsy-only annotation remediation code before rerunning it."""

from __future__ import annotations

import datetime as dt
import json

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_text, sha256_file, stable_json_sha256, write_tsv


def main() -> None:
    destination = CANDIDATE_ROOT / "BIOPSY_EXECUTION_AMENDMENT_READY.json"
    if destination.exists():
        raise RuntimeError(f"Biopsy code amendment already frozen: {destination}")
    required = [
        CANDIDATE_ROOT / "PROGRAM_EXECUTION_CODE_READY.json",
        CANDIDATE_ROOT / "EVIDENCE_CLASS_CODE_READY.json",
        CANDIDATE_ROOT / "SOURCE_AMENDMENT_02.json",
        CANDIDATE_ROOT / "sources/GSE106737/annotation/annotation_source_manifest.tsv",
    ]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    code = [SCRIPT_ROOT / "04_run_bulk_assays.R", SCRIPT_ROOT / "07_run_evidence_classes.R"]
    rows = [{"path": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size} for path in code]
    contract = {
        "freeze_id": "public-functional-biopsy-annotation-remediation-v1",
        "frozen_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "supersedes_for_GSE106737_only": [sha256_file(required[0]), sha256_file(required[1])],
        "source_amendment_sha256": json.loads(required[2].read_text(encoding="utf-8"))["amendment_sha256"],
        "annotation_manifest_sha256": sha256_file(required[3]),
        "files": rows,
        "change_boundary": "replace impossible GEO-symbol lookup with checksum-pinned AnnotationDbi PROBEID-to-SYMBOL lookup; discard zero/multi-symbol probes; all models and gates unchanged",
    }
    contract["code_bundle_sha256"] = stable_json_sha256(contract)
    write_tsv(CANDIDATE_ROOT / "biopsy_execution_amendment_manifest.tsv", rows, ["path", "sha256", "size_bytes"])
    atomic_write_text(destination, json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
