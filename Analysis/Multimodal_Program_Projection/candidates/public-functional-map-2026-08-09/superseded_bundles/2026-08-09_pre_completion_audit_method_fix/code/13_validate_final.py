#!/usr/bin/env python3
"""Independently validate multiplicity, biological units, hashes, and verdict mechanics."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from public_functional_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, require_sealed, sha256_file, write_tsv


def bh(values: np.ndarray, n_total: int) -> np.ndarray:
    output = np.full(len(values), np.nan)
    valid = np.isfinite(values)
    p = values[valid]
    order = np.argsort(p)
    ranked = p[order]
    adjusted = np.minimum.accumulate((ranked * n_total / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    restored = np.empty(len(ranked)); restored[order] = np.minimum(adjusted, 1)
    output[valid] = restored
    return output


def main() -> None:
    seal = require_sealed()
    checks = []
    def check(check_id: str, passed: bool, detail: str) -> None:
        checks.append({"check_id": check_id, "status": "PASS" if passed else "FAIL", "detail": detail})

    required = ["public_assay_registry.tsv", "sample_manifest.tsv", "source_gate_status.tsv", "source_manifest.tsv", "environment_manifest.tsv", "execution_manifest.tsv", "program_testability.tsv", "program_effects.tsv", "evidence_class_effects.tsv", "sensitivity.tsv", "cross_assay_support.tsv", "figure5_verdict.tsv", "null_permutation_audit.tsv"]
    for name in required:
        check(f"required:{name}", (CANDIDATE_ROOT / name).is_file(), name)
    programs = pd.read_csv(CANDIDATE_ROOT / "program_effects.tsv", sep="\t")
    primary = {row["program_uid"] for row in read_tsv(CANDIDATE_ROOT / "frozen_inputs/external_test_programs.tsv")}
    for keys, frame in programs.groupby(["dataset_id", "contrast_id", "scoring_scheme", "lineage"], dropna=False):
        landscape_ok = len(frame) == 117 and np.allclose(bh(frame.p.to_numpy(float), 117), frame.landscape_117_q.to_numpy(float), equal_nan=True, atol=1e-12)
        check("BH117:" + ":".join(map(str, keys)), landscape_ok, f"n={len(frame)}")
        primary_frame = frame[frame.program_uid.isin(primary)]
        primary_ok = len(primary_frame) == 2 and np.allclose(bh(primary_frame.p.to_numpy(float), 2), primary_frame.primary_family_q.to_numpy(float), equal_nan=True, atol=1e-12)
        check("BH2:" + ":".join(map(str, keys)), primary_ok, f"n={len(primary_frame)}")
    gates = pd.DataFrame(read_tsv(CANDIDATE_ROOT / "source_gate_status.tsv"))
    expected_gates = {"GSE200418": "pass", "GSE207889": "pass_with_source_limit", "GSE253380": "pass_model_specific_no_line_key", "GSE106737": "pass", "CLCC1": "skipped_source_unavailable", "MYOJIN_HLF": "pass_reuse_only"}
    observed = dict(zip(gates.dataset_id, gates.status))
    check("source_gate_states", observed == expected_gates, str(observed))
    source_manifest = pd.read_csv(CANDIDATE_ROOT / "source_manifest.tsv", sep="\t")
    check("OS_H5AD_provenance", bool((source_manifest.file_id == "os_h5ad_source_amendment").any()), "same-study lineage remediation must be visible")
    check("GPL_annotation_provenance", bool((source_manifest.file_id == "GPL16686_annotation_sqlite").any()), "official array annotation remediation must be visible")
    samples = pd.read_csv(CANDIDATE_ROOT / "sample_manifest.tsv", sep="\t")
    check("PCLS_donors", samples[(samples.dataset_id == "GSE200418") & samples.include_in_inference].biological_unit_id.nunique() == 7, "must be seven donors")
    check("biopsy_participant_rows", len(samples[(samples.dataset_id == "GSE106737") & samples.include_in_inference]) == 82, "41 paired participants x 2")
    check("schlo_source_units", samples[(samples.dataset_id == "GSE207889") & samples.include_in_inference].replicate.nunique() == 2, "cells are not biological n")
    nulls = pd.read_csv(CANDIDATE_ROOT / "null_permutation_audit.tsv", sep="\t")
    check("centered_nulls", bool(nulls.null_centered.all()), f"{nulls.null_centered.sum()}/{len(nulls)}")
    frozen_manifest = pd.read_csv(CANDIDATE_ROOT / "frozen_input_manifest.tsv", sep="\t")
    frozen_ok = all(sha256_file(PROJECT_ROOT / row.snapshot_path) == row.sha256 for row in frozen_manifest.itertuples())
    check("frozen_snapshot_hashes", frozen_ok, f"n={len(frozen_manifest)}")
    verdict = pd.read_csv(CANDIDATE_ROOT / "figure5_verdict.tsv", sep="\t")
    workstream = verdict[(verdict.object_type == "workstream") & (verdict.object_id == "public_functional_map")]
    computed = bool((verdict[verdict.object_type == "program"].main_figure_eligible == True).any())
    check("mechanical_verdict", len(workstream) == 1 and bool(workstream.main_figure_eligible.iloc[0]) == computed, f"computed={computed}")
    check("no_secondary_promotion", set(verdict[verdict.object_type == "program"].object_id) == primary, "only two externally frozen programs adjudicated")
    failures = [row for row in checks if row["status"] == "FAIL"]
    write_tsv(CANDIDATE_ROOT / "validation_report.tsv", checks, ["check_id", "status", "detail"])
    if failures:
        raise SystemExit(f"FINAL_VALIDATION_FAILED: {len(failures)}/{len(checks)}")
    manifest_rows = []
    for path in sorted(CANDIDATE_ROOT.rglob("*")):
        if path.is_file() and "sources" not in path.parts and path.name != "release_manifest.tsv":
            manifest_rows.append({"path": str(path.relative_to(PROJECT_ROOT)), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)})
    write_tsv(CANDIDATE_ROOT / "release_manifest.tsv", manifest_rows, ["path", "size_bytes", "sha256"])
    ready = {"status": "validated", "checks_passed": len(checks), "checks_failed": 0, "main_figure_eligible": computed, "specification_sha256": seal["specification_sha256"], "release_manifest_sha256": sha256_file(CANDIDATE_ROOT / "release_manifest.tsv"), "canonical_promotion_authorized": False}
    from public_functional_common import atomic_write_text
    atomic_write_text(CANDIDATE_ROOT / "VALIDATED.json", json.dumps(ready, indent=2, sort_keys=True) + "\n")
    print(json.dumps(ready, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
