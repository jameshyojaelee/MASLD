#!/usr/bin/env python3
"""Freeze transparent SQLite→Entrez→symbol remediation after loadDb failure."""

from __future__ import annotations

import datetime as dt
import json

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_text, sha256_file, stable_json_sha256, write_tsv


def main() -> None:
    destination = CANDIDATE_ROOT / "BIOPSY_MAPPING_V2_READY.json"
    if destination.exists():
        raise RuntimeError(f"Biopsy mapping v2 already frozen: {destination}")
    error_logs = [CANDIDATE_ROOT / "logs/bulk_19655578_2.err", CANDIDATE_ROOT / "logs/class_19655528_2.err"]
    expected = "invalid data of mode 'character' (too short)"
    if not all(path.is_file() and expected in path.read_text(encoding="utf-8") for path in error_logs):
        raise RuntimeError("Expected AnnotationDbi::loadDb incompatibility was not reproduced in both branches")
    code = [SCRIPT_ROOT / "04_run_bulk_assays.R", SCRIPT_ROOT / "07_run_evidence_classes.R"]
    rows = [{"path": str(path.relative_to(PROJECT_ROOT)), "sha256": sha256_file(path), "size_bytes": path.stat().st_size} for path in code]
    contract = {
        "freeze_id": "public-functional-biopsy-mapping-v2",
        "frozen_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "files": rows,
        "failure_class": "current AnnotationDbi loadDb cannot instantiate the legacy ChipDb package SQLite despite valid schema and checksum",
        "mapping": "RSQLite probes(probe_id,gene_id,is_multiple=0) then official org.Hs.eg.db ENTREZID-to-SYMBOL; discard ambiguous/missing symbols",
        "scientific_change": "none; same outcome-blind mapping and probe-selection rules",
        "sqlite_sha256": "fa6cd1ec2c75cf7af0579676dcdcff6f58488900e1f91d4f2f505063e64ce34f",
    }
    contract["code_bundle_sha256"] = stable_json_sha256(contract)
    write_tsv(CANDIDATE_ROOT / "biopsy_mapping_v2_code_manifest.tsv", rows, ["path", "sha256", "size_bytes"])
    atomic_write_text(destination, json.dumps(contract, indent=2, sort_keys=True) + "\n")
    print(json.dumps(contract, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
