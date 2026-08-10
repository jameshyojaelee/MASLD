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


def as_bool(series: pd.Series) -> pd.Series:
    """Parse serialized booleans fail-closed instead of relying on pandas inference."""
    normalized = series.astype(str).str.strip().str.lower()
    if not normalized.isin({"true", "false"}).all():
        invalid = sorted(normalized[~normalized.isin({"true", "false"})].unique())
        raise RuntimeError(f"Invalid serialized boolean values: {invalid}")
    return normalized.eq("true")


def main() -> None:
    seal = require_sealed()
    checks = []
    def check(check_id: str, passed: bool, detail: str) -> None:
        checks.append({"check_id": check_id, "status": "PASS" if passed else "FAIL", "detail": detail})

    required = ["public_assay_registry.tsv", "sample_manifest.tsv", "source_gate_status.tsv", "source_schema_audit.tsv", "source_manifest.tsv", "environment_manifest.tsv", "execution_manifest.tsv", "program_testability.tsv", "program_effects.tsv", "evidence_class_effects.tsv", "evidence_class_sensitivity.tsv", "sensitivity.tsv", "cross_assay_support.tsv", "figure5_verdict.tsv", "null_permutation_audit.tsv", "COMPLETION_AUDIT_AMENDMENT_01.json", "METHOD_COMPLIANCE_EXECUTION_READY.json", "RELEASE_MANIFEST_LOG_AMENDMENT_01.json", "TERMINAL_REVALIDATION_READY.json"]
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
    schema = pd.read_csv(CANDIDATE_ROOT / "source_schema_audit.tsv", sep="\t")
    accounting = schema[(schema.dataset_id == "GSE200418") & schema.item.isin(["source_attempted_samples", "source_qc_excluded_samples", "source_deposited_after_qc"])]
    accounting_observed = dict(zip(accounting.item, accounting.value.astype(int)))
    check("PCLS_source_accounting", accounting_observed == {"source_attempted_samples": 112, "source_qc_excluded_samples": 13, "source_deposited_after_qc": 99} and bool(as_bool(accounting["pass"]).all()), str(accounting_observed))
    source_manifest = pd.read_csv(CANDIDATE_ROOT / "source_manifest.tsv", sep="\t")
    check("OS_H5AD_provenance", bool((source_manifest.file_id == "os_h5ad_source_amendment").any()), "same-study lineage remediation must be visible")
    check("GPL_annotation_provenance", bool((source_manifest.file_id == "GPL16686_annotation_sqlite").any()), "official array annotation remediation must be visible")
    samples = pd.read_csv(CANDIDATE_ROOT / "sample_manifest.tsv", sep="\t")
    included = as_bool(samples.include_in_inference)
    check("PCLS_donors", samples[(samples.dataset_id == "GSE200418") & included].biological_unit_id.nunique() == 7, "must be seven donors")
    check("biopsy_participant_rows", len(samples[(samples.dataset_id == "GSE106737") & included]) == 82, "41 paired participants x 2")
    check("schlo_source_units", samples[(samples.dataset_id == "GSE207889") & included].replicate.nunique() == 2, "cells are not biological n")
    schlo = programs[programs.dataset_id == "GSE207889"]
    check("schlo_assay_native_method", set(schlo.analysis_method.dropna()) == {"edgeR_TMM_voom_then_paired_limma_empirical_bayes"} and len(schlo) == 4212, "all 4 lineages x 3 contrasts x 3 schemes x 117 programs must use the frozen method")
    sensitivities = pd.read_csv(CANDIDATE_ROOT / "sensitivity.tsv", sep="\t", low_memory=False)
    primary_lodo = sensitivities[(sensitivities.dataset_id == "GSE200418") & sensitivities.program_uid.isin(primary) & (sensitivities.contrast_id == "PCLS_GFIPO_GFI_48h") & (sensitivities.sensitivity == "leave_one_donor_out")]
    check("PCLS_LODO_contributing_donors_only", len(primary_lodo) == 12 and primary_lodo.groupby("program_uid").held_out_unit.nunique().eq(6).all(), "two primary programs x six contributing donors")
    probe_sensitivity = sensitivities[(sensitivities.dataset_id == "GSE106737") & (sensitivities.sensitivity == "median_across_probes")]
    check("biopsy_program_probe_sensitivity", len(probe_sensitivity) == 234 and set(probe_sensitivity.contrast_id) == {"LSI_RESPONSE_DIFF", "RYGB_RESPONSE"} and probe_sensitivity.groupby("contrast_id").size().eq(117).all(), "117 programs x two biopsy contrasts")
    class_sensitivity = pd.read_csv(CANDIDATE_ROOT / "evidence_class_sensitivity.tsv", sep="\t")
    median_camera = class_sensitivity[(class_sensitivity.dataset_id == "GSE106737") & (class_sensitivity.method == "camera_median_across_unambiguous_probes")]
    check("biopsy_class_probe_sensitivity", len(median_camera) == 3 and set(median_camera.evidence_class) == {"disease_state_only", "genetic_only", "convergent"}, "three frozen evidence classes")
    nulls = pd.read_csv(CANDIDATE_ROOT / "null_permutation_audit.tsv", sep="\t")
    null_centered = as_bool(nulls.null_centered)
    check("centered_nulls", bool(null_centered.all()), f"{int(null_centered.sum())}/{len(nulls)}")
    frozen_manifest = pd.read_csv(CANDIDATE_ROOT / "frozen_input_manifest.tsv", sep="\t")
    frozen_ok = all(sha256_file(PROJECT_ROOT / row.snapshot_path) == row.sha256 for row in frozen_manifest.itertuples())
    check("frozen_snapshot_hashes", frozen_ok, f"n={len(frozen_manifest)}")
    method_manifest = pd.read_csv(CANDIDATE_ROOT / "method_compliance_code_manifest.tsv", sep="\t")
    superseded_validator = "Analysis/Multimodal_Program_Projection/scripts/public_functional/13_validate_final.py"
    unchanged_method_rows = method_manifest[method_manifest.path != superseded_validator]
    method_ok = all(sha256_file(PROJECT_ROOT / row.path) == row.sha256 for row in unchanged_method_rows.itertuples())
    check("method_compliance_code_hashes", method_ok and len(method_manifest) == 20 and len(unchanged_method_rows) == 19, f"unchanged={len(unchanged_method_rows)}/19; validator superseded by terminal amendment")
    terminal_manifest = pd.read_csv(CANDIDATE_ROOT / "terminal_revalidation_code_manifest.tsv", sep="\t")
    terminal_ok = all(sha256_file(PROJECT_ROOT / row.path) == row.sha256 for row in terminal_manifest.itertuples())
    check("terminal_revalidation_code_hashes", terminal_ok and len(terminal_manifest) == 3, f"n={len(terminal_manifest)}")
    archive_manifest = pd.read_csv(CANDIDATE_ROOT / "superseded_bundles/2026-08-09_pre_completion_audit_method_fix/archive_manifest.tsv", sep="\t")
    archive_ok = all(sha256_file(PROJECT_ROOT / row.archive) == row.sha256 for row in archive_manifest.itertuples())
    check("pre_remediation_archive_hashes", archive_ok, f"n={len(archive_manifest)}")
    verdict = pd.read_csv(CANDIDATE_ROOT / "figure5_verdict.tsv", sep="\t")
    workstream = verdict[(verdict.object_type == "workstream") & (verdict.object_id == "public_functional_map")]
    program_eligible = as_bool(verdict.loc[verdict.object_type == "program", "main_figure_eligible"])
    computed = bool(program_eligible.any())
    workstream_eligible = as_bool(workstream.main_figure_eligible) if len(workstream) else pd.Series(dtype=bool)
    check("mechanical_verdict", len(workstream) == 1 and bool(workstream_eligible.iloc[0]) == computed, f"computed={computed}")
    check("no_secondary_promotion", set(verdict[verdict.object_type == "program"].object_id) == primary, "only two externally frozen programs adjudicated")
    excluded_release_parts = {"sources", "superseded_attempts", "superseded_bundles", "logs"}
    # A file cannot truthfully contain its own checksum.  VALIDATED.json carries
    # the terminal manifest hash, so both terminal marker files are deliberately
    # outside the checksummed payload rather than being recursively self-hashed.
    excluded_release_names = {"release_manifest.tsv", "VALIDATED.json"}
    manifest_rows = []
    for path in sorted(CANDIDATE_ROOT.rglob("*")):
        if path.is_file() and not excluded_release_parts.intersection(path.parts) and path.name not in excluded_release_names:
            manifest_rows.append({"path": str(path.relative_to(PROJECT_ROOT)), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)})
    check("active_manifest_excludes_superseded", all("superseded_attempts" not in row["path"] and "superseded_bundles" not in row["path"] for row in manifest_rows), f"n={len(manifest_rows)}")
    check("active_manifest_excludes_mutable_logs", all("logs" not in Path(row["path"]).parts for row in manifest_rows), f"n={len(manifest_rows)}")
    check("terminal_markers_not_self_hashed", all(Path(row["path"]).name not in excluded_release_names for row in manifest_rows), "release_manifest.tsv and VALIDATED.json are intentionally outside the checksummed payload")
    failures = [row for row in checks if row["status"] == "FAIL"]
    write_tsv(CANDIDATE_ROOT / "validation_report.tsv", checks, ["check_id", "status", "detail"])
    if failures:
        raise SystemExit(f"FINAL_VALIDATION_FAILED: {len(failures)}/{len(checks)}")
    manifest_rows = []
    for path in sorted(CANDIDATE_ROOT.rglob("*")):
        if path.is_file() and not excluded_release_parts.intersection(path.parts) and path.name not in excluded_release_names:
            manifest_rows.append({"path": str(path.relative_to(PROJECT_ROOT)), "size_bytes": path.stat().st_size, "sha256": sha256_file(path)})
    write_tsv(CANDIDATE_ROOT / "release_manifest.tsv", manifest_rows, ["path", "size_bytes", "sha256"])
    ready = {"status": "validated", "checks_passed": len(checks), "checks_failed": 0, "main_figure_eligible": computed, "specification_sha256": seal["specification_sha256"], "release_manifest_sha256": sha256_file(CANDIDATE_ROOT / "release_manifest.tsv"), "canonical_promotion_authorized": False}
    from public_functional_common import atomic_write_text
    atomic_write_text(CANDIDATE_ROOT / "VALIDATED.json", json.dumps(ready, indent=2, sort_keys=True) + "\n")
    print(json.dumps(ready, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
