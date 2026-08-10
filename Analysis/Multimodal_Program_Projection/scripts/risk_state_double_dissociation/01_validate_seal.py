#!/usr/bin/env python3
"""Independently validate the Plan 42 seal before source outcomes are acquired."""

from __future__ import annotations

import json
from collections import Counter

from risk_state_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_text,
    read_tsv,
    sha256_file,
    stable_json_sha256,
)


EXPECTED_CLASSES = {
    "genetic_only": 413,
    "disease_state_only": 1227,
    "convergent": 34,
    "neither": 13257,
    "indeterminate_not_jointly_testable": 16298,
}


def main() -> None:
    seal_path = CANDIDATE_ROOT / "SEALED.json"
    if not seal_path.is_file():
        raise FileNotFoundError(seal_path)
    if (CANDIDATE_ROOT / "source").exists():
        raise RuntimeError("Source directory exists before seal validation; outcome firewall violated")
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    specification_sha = seal.pop("specification_sha256")
    if stable_json_sha256(seal) != specification_sha:
        raise RuntimeError("Specification digest cannot be independently rederived")

    checks: list[dict[str, object]] = []
    for row in read_tsv(CANDIDATE_ROOT / "frozen_input_manifest.tsv"):
        source = PROJECT_ROOT / row["source_path"]
        snapshot = PROJECT_ROOT / row["snapshot_path"]
        observed_source = sha256_file(source)
        observed_snapshot = sha256_file(snapshot)
        passed = observed_source == row["sha256"] == observed_snapshot
        checks.append({"check": f"frozen_input:{row['input_id']}", "passed": passed, "detail": observed_source})
    classes = read_tsv(CANDIDATE_ROOT / "frozen_inputs/evidence_classes__frozen_evidence_classes.tsv")
    counts = Counter(row["static_class"] for row in classes)
    checks.append({"check": "evidence_class_counts", "passed": dict(counts) == EXPECTED_CLASSES, "detail": json.dumps(dict(counts), sort_keys=True)})
    loci = read_tsv(CANDIDATE_ROOT / "genetic_locus_registry.tsv")
    representatives = Counter(row["locus_uid"] for row in loci if row["is_representative"] == "true")
    locus_ids = {row["locus_uid"] for row in loci}
    checks.append({"check": "genetic_rows_413", "passed": len(loci) == 413, "detail": len(loci)})
    checks.append({"check": "one_representative_per_locus", "passed": set(representatives) == locus_ids and all(value == 1 for value in representatives.values()), "detail": f"{len(representatives)}/{len(locus_ids)}"})
    checks.append({"check": "locus_registry_digest", "passed": sha256_file(CANDIDATE_ROOT / "genetic_locus_registry.tsv") == seal["locus_registry_sha256"], "detail": seal["locus_registry_sha256"]})
    for path in CANDIDATE_ROOT.joinpath("frozen_spec").glob("*.tsv"):
        checks.append({"check": f"nonempty_config:{path.name}", "passed": path.stat().st_size > 0, "detail": path.stat().st_size})

    failures = [row for row in checks if not row["passed"]]
    validation = {
        "candidate_id": seal["candidate_id"],
        "status": "passed" if not failures else "failed",
        "seal_sha256": sha256_file(seal_path),
        "specification_sha256": specification_sha,
        "n_checks": len(checks),
        "n_failed": len(failures),
        "checks": checks,
    }
    if failures:
        raise RuntimeError(json.dumps(validation, indent=2, sort_keys=True))
    atomic_write_text(
        CANDIDATE_ROOT / "SEAL_VALIDATED.json",
        json.dumps(validation, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()

