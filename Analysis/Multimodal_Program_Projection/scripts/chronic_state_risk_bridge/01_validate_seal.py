#!/usr/bin/env python3
"""Independently validate the Plan 43 outcome-blind seal."""

from __future__ import annotations

import json
import math
from collections import Counter

from bridge_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_text,
    read_tsv,
    sha256_file,
)


def main() -> None:
    seal_path = CANDIDATE_ROOT / "SEALED.json"
    if not seal_path.is_file():
        raise FileNotFoundError(seal_path)
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    checks: list[dict[str, object]] = []

    def add(name: str, passed: bool, detail: object) -> None:
        checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})

    inputs = read_tsv(CANDIDATE_ROOT / "frozen_input_manifest.tsv")
    for row in inputs:
        path = PROJECT_ROOT / row["snapshot_path"]
        add(f"input_hash:{row['input_id']}", path.is_file() and sha256_file(path) == row["sha256"], row["snapshot_path"])
    coloc = read_tsv(CANDIDATE_ROOT / "tier12_coloc_manifest.tsv")
    add("tier12_study_count", len(coloc) == 35, len(coloc))
    for row in coloc:
        path = PROJECT_ROOT / row["snapshot_path"]
        add(f"tier12_hash:{row['study_name']}", path.is_file() and sha256_file(path) == row["sha256"], row["snapshot_path"])

    loci = read_tsv(CANDIDATE_ROOT / "corrected_genetic_locus_registry.tsv")
    locus_ids = {row["tier12_locus_uid"] for row in loci}
    representatives = Counter(row["tier12_locus_uid"] for row in loci if row["is_representative"] == "true")
    add("genetic_row_count", len(loci) == 413, len(loci))
    add("corrected_locus_count", len(locus_ids) == 326, len(locus_ids))
    add("one_representative_per_locus", set(representatives) == locus_ids and all(value == 1 for value in representatives.values()), len(representatives))
    add("driver_correction_count", sum(row["driver_corrected"] == "true" for row in loci) == 80, sum(row["driver_corrected"] == "true" for row in loci))
    add("tier12_only", all(row["tier12_driving_tier"] in {"1", "2"} for row in loci), Counter(row["tier12_driving_tier"] for row in loci))
    add("lead_snps_resolved", all(row["tier12_top_snp"] and row["position"] for row in loci), sum(not row["position"] for row in loci))

    classes = read_tsv(CANDIDATE_ROOT / "corrected_evidence_classes.tsv")
    expected = {"genetic_only": 413, "disease_state_only": 1227, "convergent": 34, "neither": 13257, "indeterminate_not_jointly_testable": 16298}
    add("class_counts_unchanged", dict(Counter(row["static_class"] for row in classes)) == expected, Counter(row["static_class"] for row in classes))

    axis = read_tsv(CANDIDATE_ROOT / "frozen_state_axis.tsv")
    primary = [float(row["primary_loading"]) for row in axis]
    lineage = [float(row["lineage_balanced_loading"]) for row in axis]
    add("state_axis_117", len(axis) == 117, len(axis))
    add("state_axis_primary_l1", math.isclose(sum(abs(value) for value in primary), 1.0, rel_tol=0, abs_tol=1e-12), sum(abs(value) for value in primary))
    add("state_axis_lineage_l1", math.isclose(sum(abs(value) for value in lineage), 1.0, rel_tol=0, abs_tol=1e-12), sum(abs(value) for value in lineage))
    add("state_axis_no_selection", all(row["stage_beta"] != "" for row in axis), "all 117 betas present")

    hash_fields = {
        "frozen_input_manifest_sha256": "frozen_input_manifest.tsv",
        "tier12_coloc_manifest_sha256": "tier12_coloc_manifest.tsv",
        "corrected_locus_registry_sha256": "corrected_genetic_locus_registry.tsv",
        "corrected_evidence_classes_sha256": "corrected_evidence_classes.tsv",
        "frozen_state_axis_sha256": "frozen_state_axis.tsv",
        "config_manifest_sha256": "config_manifest.tsv",
    }
    for field, name in hash_fields.items():
        add(f"seal_hash:{field}", seal.get(field) == sha256_file(CANDIDATE_ROOT / name), name)
    add("outcomes_closed", seal.get("source_outcomes_loaded_before_seal") is False, seal.get("source_outcomes_loaded_before_seal"))

    failures = [row for row in checks if not row["passed"]]
    validation = {
        "status": "pass" if not failures else "fail",
        "n_checks": len(checks),
        "n_passed": len(checks) - len(failures),
        "n_failed": len(failures),
        "seal_sha256": sha256_file(seal_path),
        "specification_sha256": seal["specification_sha256"],
        "checks": checks,
    }
    atomic_write_text(
        CANDIDATE_ROOT / "SEAL_VALIDATED.json",
        json.dumps(validation, indent=2, sort_keys=True) + "\n",
    )
    if failures:
        raise RuntimeError(f"Plan 43 seal validation failed: {failures}")
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
