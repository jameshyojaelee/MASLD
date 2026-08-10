#!/usr/bin/env python3
"""Independently validate Plan 43 outputs and freeze the terminal release."""

from __future__ import annotations

import datetime as dt
import json
import math
from pathlib import Path

from bridge_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_text,
    read_tsv,
    require_validated_seal,
    sha256_file,
    write_tsv,
)


def number(value: str) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def truth(value: str) -> bool:
    return str(value).strip().lower() == "true"


def bh(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    output = [math.nan] * len(values)
    running = 1.0
    for reverse_rank in range(len(values) - 1, -1, -1):
        index = order[reverse_rank]
        running = min(running, values[index] * len(values) / (reverse_rank + 1))
        output[index] = min(running, 1.0)
    return output


def main() -> None:
    seal = require_validated_seal()
    checks: list[dict[str, object]] = []

    def add(check_id: str, passed: bool, detail: object) -> None:
        checks.append({"check_id": check_id, "passed": str(bool(passed)).lower(), "detail": detail})

    for row in read_tsv(CANDIDATE_ROOT / "frozen_input_manifest.tsv"):
        source = PROJECT_ROOT / row["source_path"]
        snapshot = PROJECT_ROOT / row["snapshot_path"]
        passed = source.is_file() and snapshot.is_file() and sha256_file(source) == row["sha256"] == sha256_file(snapshot)
        add(f"immutable_input:{row['input_id']}", passed, row["sha256"])

    source_manifest = read_tsv(CANDIDATE_ROOT / "source_manifest.tsv")
    source_hashes_ok = True
    for row in source_manifest:
        if not truth(row["acquired"]):
            source_hashes_ok = False
            continue
        source = CANDIDATE_ROOT / row["relative_path"]
        source_hashes_ok &= source.is_file() and source.stat().st_size == int(row["size_bytes"]) and sha256_file(source) == row["sha256"]
    add("source_payload_checksums", source_hashes_ok and len(source_manifest) > 0, len(source_manifest))

    registry = read_tsv(CANDIDATE_ROOT / "frozen_inputs/program_registry__program_registry_v2.tsv")
    programs = {row["program_uid"] for row in registry}
    membership = read_tsv(CANDIDATE_ROOT / "frozen_inputs/program_membership__program_membership_v2.tsv")
    add("complete_117_program_registry", len(registry) == 117 and len(programs) == 117, len(registry))
    add("membership_only_frozen_programs", {row["program_uid"] for row in membership} == programs, len(membership))
    axis = read_tsv(CANDIDATE_ROOT / "frozen_state_axis.tsv")
    add("state_axis_all_117_no_significance_filter", len(axis) == 117 and {row["program_uid"] for row in axis} == programs, len(axis))
    add("state_axis_absolute_loading_sums_one", abs(sum(abs(number(row["primary_loading"])) for row in axis) - 1) < 1e-10, "sum_abs")

    loci = read_tsv(CANDIDATE_ROOT / "corrected_genetic_locus_registry.tsv")
    corrected_classes = read_tsv(CANDIDATE_ROOT / "integrity_correction_v2/corrected_evidence_classes_v2.tsv")
    add(
        "corrected_tier12_locus_count",
        len(loci) == 413
        and len({row["tier12_locus_uid"] for row in loci}) == 326
        and sum(truth(row["is_representative"]) for row in loci) == 326,
        f"gene_rows={len(loci)}; loci={len({row['tier12_locus_uid'] for row in loci})}",
    )
    add("genetic_class_row_count", sum(row["static_class"] == "genetic_only" for row in corrected_classes) == 413, "expected=413")
    add("integrity_manifest_promoted_nonoverwriting", (CANDIDATE_ROOT / "integrity_correction_manifest.tsv").is_file(), "root copy")

    human = read_tsv(CANDIDATE_ROOT / "human_reversal_effects.tsv")
    primary = [row for row in human if row["score_scheme"] == "weighted__primary"]
    observed_q = [number(row["q"]) for row in primary]
    expected_q = bh([number(row["p"]) for row in primary])
    add("human_primary_three_datasets", len(primary) == 3, len(primary))
    add("human_primary_BH_rederived", all(abs(x - y) < 1e-10 for x, y in zip(observed_q, expected_q)), f"expected={expected_q}; observed={observed_q}")
    direction = {row["dataset_id"]: row["expected_direction"] for row in primary}
    add("human_direction_convention", direction == {"GSE106737": "negative", "GSE83452": "negative", "GSE48452": "positive"}, direction)
    exact_counts = {row["dataset_id"]: int(float(row["n_permutations"])) for row in primary}
    add("human_exact_permutations", exact_counts == {"GSE106737": 184756, "GSE83452": 74613, "GSE48452": 10810800}, exact_counts)
    meta = read_tsv(CANDIDATE_ROOT / "human_reversal/human_reversal_meta_analysis.tsv")[0]
    add("overlap_collapsed_to_two_cohorts", int(meta["n_effective_cohorts"]) == 2 and not truth(meta["three_cohort_preregistered_gate_evaluable"]), meta["reason"])
    fingerprints = read_tsv(CANDIDATE_ROOT / "source_gates/sample_fingerprint_matches.tsv")
    reciprocal = [row for row in fingerprints if truth(row.get("reciprocal_best", "false")) and number(row.get("pearson_r", "")) >= 0.9999]
    add("cohort_overlap_expression_confirmed", len(reciprocal) >= 78, len(reciprocal))

    crop = read_tsv(CANDIDATE_ROOT / "cropseq_risk_bridge.tsv")[0]
    add("cropseq_fail_closed_before_outcome", crop["status"] == "skipped_insufficient_oriented_loci" and crop["estimate"] == "" and int(crop["n_eligible_loci"]) == 5, crop["status"])

    protein = read_tsv(CANDIDATE_ROOT / "protein_transportability.tsv")
    protein_map = {row["dataset_id"]: row for row in protein}
    add("protein_unseen_source_count", int(protein_map["PXD052787"]["n_pairs"]) == 88, protein_map["PXD052787"]["n_pairs"])
    add("protein_unseen_not_significant", number(protein_map["PXD052787"]["q"]) >= 0.05, protein_map["PXD052787"]["q"])
    add("protein_matching_gate_failed_explicitly", not truth(protein_map["PXD052787"]["matching_gate_pass"]) and not truth(protein_map["PXD051911"]["matching_gate_pass"]), "both false")

    context = read_tsv(CANDIDATE_ROOT / "context_topology.tsv")
    pcls = [row for row in context if row["dataset_id"] == "GSE200418" and row["score_scheme"] == "weighted__primary"]
    schlo = [row for row in context if row["dataset_id"] == "GSE207889" and row["score_scheme"] == "weighted__primary"]
    add("PCLS_donor_units_not_samples", all(int(row["n_biological_units"]) <= 6 for row in pcls), [row["n_biological_units"] for row in pcls])
    add("scHLO_replicates_not_cells", all(int(row["n_biological_units"]) == 2 for row in schlo), len(schlo))

    spatial = read_tsv(CANDIDATE_ROOT / "spatial_lineage_support.tsv")
    add("four_frozen_program_spatial_rows", sum(row["analysis_id"] == "accepted_program_spatial_support" for row in spatial) == 4, len(spatial))
    add("spatial_class_comparison_fail_closed", any(row["status"].startswith("skipped_missing_prespecified") for row in spatial), "missing covariates")
    identity = read_tsv(CANDIDATE_ROOT / "program_identity_audit.tsv")
    add("M8_M20_identity_complete", len(identity) == 2 and all(int(row["n_donors"]) >= 200 for row in identity), [(row["program_uid"], row["n_donors"]) for row in identity])
    add("doublet_limitation_explicit", all(not truth(row["authoritative_doublet_field_available"]) for row in identity), "no source doublet field")

    clcc = read_tsv(CANDIDATE_ROOT / "clcc1_integrity_correction/crispr_class_effects.tsv")
    add("clcc_corrected_primary_family", len(clcc) == 3 and all(number(row["q"]) >= 0.05 for row in clcc), [(row["comparison"], row["q"]) for row in clcc])
    add("clcc_correction_nonconfirmatory", read_tsv(CANDIDATE_ROOT / "clcc1_integrity_correction/crispr_verdict.tsv")[0]["status"] == "valid_null_or_nonconfirmatory", "supplementary")

    verdict = {row["component_id"]: row for row in read_tsv(CANDIDATE_ROOT / "promotion_verdict.tsv")}
    component_ids = {"human_histologic_reversal", "chronic_context_topology", "regulatory_risk_cropseq_bridge", "protein_translation", "spatial_lineage_remodeling", "complete_risk_to_state_architecture"}
    add("mechanical_component_verdicts_present", component_ids.issubset(verdict), sorted(verdict))
    add("full_escalation_false", not truth(verdict["complete_risk_to_state_architecture"]["gate_pass"]), verdict["complete_risk_to_state_architecture"]["evidence"])

    mandatory = [
        "integrity_correction_manifest.tsv", "source_gate_status.tsv", "frozen_state_axis.tsv",
        "paired_sample_manifest.tsv", "human_reversal_effects.tsv", "context_topology.tsv",
        "cropseq_risk_bridge.tsv", "protein_transportability.tsv", "spatial_lineage_support.tsv",
        "program_identity_audit.tsv", "promotion_verdict.tsv", "source_manifest.tsv", "environment_manifest.tsv", "execution_manifest.tsv",
    ]
    for relative in mandatory:
        path = CANDIDATE_ROOT / relative
        add(f"mandatory_output:{relative}", path.is_file() and path.stat().st_size > 0, path.stat().st_size if path.is_file() else 0)

    failures = [row for row in checks if row["passed"] != "true"]
    write_tsv(CANDIDATE_ROOT / "validation_report.tsv", checks, ["check_id", "passed", "detail"])
    if failures:
        raise RuntimeError(json.dumps(failures, indent=2, sort_keys=True))

    excluded_parts = {"logs", "sources", "__pycache__"}
    excluded_names = {"release_manifest.tsv", "VALIDATED.json", "COMPLETE"}
    release = []
    for path in sorted(CANDIDATE_ROOT.rglob("*")):
        if not path.is_file() or excluded_parts.intersection(path.parts) or path.name in excluded_names:
            continue
        release.append({
            "relative_path": str(path.relative_to(CANDIDATE_ROOT)),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    write_tsv(CANDIDATE_ROOT / "release_manifest.tsv", release, ["relative_path", "bytes", "sha256"])
    validation = {
        "candidate_id": seal["candidate_id"],
        "status": "complete_no_promotion",
        "validated_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "checks_passed": len(checks), "checks_failed": 0,
        "release_payloads": len(release),
        "release_manifest_sha256": sha256_file(CANDIDATE_ROOT / "release_manifest.tsv"),
        "figure5_promotion": False,
        "canonical_promotion_authorized": False,
    }
    atomic_write_text(CANDIDATE_ROOT / "VALIDATED.json", json.dumps(validation, indent=2, sort_keys=True) + "\n")
    atomic_write_text(CANDIDATE_ROOT / "COMPLETE", validation["release_manifest_sha256"] + "\n")
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
