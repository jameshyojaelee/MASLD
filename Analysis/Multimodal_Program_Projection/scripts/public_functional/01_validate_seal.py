#!/usr/bin/env python3
"""Independently validate the outcome-blind public functional specification."""

from __future__ import annotations

import json

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def main() -> None:
    seal_path = CANDIDATE_ROOT / "SEALED.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    checks: list[tuple[str, bool, str]] = []
    checks.append(("candidate_id", seal["candidate_id"] == "public-functional-map-2026-08-09", str(seal["candidate_id"])))
    checks.append(("registry_rows", seal["registry_rows"] == 117, str(seal["registry_rows"])))
    checks.append(("primary_program_count", len(seal["primary_program_uids"]) == 2, str(len(seal["primary_program_uids"]))))
    checks.append(("primary_membership_rows", seal["primary_membership_rows"] == 71, str(seal["primary_membership_rows"])))
    for row in read_tsv(CANDIDATE_ROOT / "frozen_input_manifest.tsv"):
        source = PROJECT_ROOT / row["source_path"]
        snapshot = PROJECT_ROOT / row["snapshot_path"]
        checks.append((f"source_hash:{row['input_id']}", source.is_file() and sha256_file(source) == row["sha256"], row["source_path"]))
        checks.append((f"snapshot_hash:{row['input_id']}", snapshot.is_file() and sha256_file(snapshot) == row["sha256"], row["snapshot_path"]))
    for name in ["public_assays.tsv", "source_files.tsv", "contrasts.tsv", "frozen_inputs.tsv"]:
        checks.append((f"frozen_config:{name}", (CANDIDATE_ROOT / "frozen_spec" / name).is_file(), name))
    failures = [check for check in checks if not check[1]]
    for check_id, passed, detail in checks:
        print(f"{check_id}\t{'PASS' if passed else 'FAIL'}\t{detail}")
    if failures:
        raise SystemExit(f"Seal validation failed: {len(failures)} checks")
    print(f"SEAL_VALIDATED\t{len(checks)}/{len(checks)}")


if __name__ == "__main__":
    main()
