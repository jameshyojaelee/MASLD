#!/usr/bin/env python3
"""Independently validate the spatial Resource candidate and release products."""

from __future__ import annotations

import argparse
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from spatial_resource_lib import (
    EXPECTED_N_CONFIRMATORY,
    EXPECTED_N_PROGRAMS,
    EXPECTED_REGISTRY_SHA256,
    PROGRAM_RELEASE_ID,
    RESOURCE_RELEASE_ID,
    SpatialResourceError,
    load_dataset_registry,
    load_frozen_programs,
    parse_bool,
    read_tsv,
    require_complete_coverage,
    sha256_file,
    write_tsv,
)


def check(condition: bool, check_id: str, observed: object, expected: object, rows: list[dict[str, object]]) -> None:
    rows.append({
        "check_id": check_id,
        "status": "PASS" if condition else "FAIL",
        "observed": str(observed),
        "expected": str(expected),
    })
    if not condition:
        raise SpatialResourceError(f"{check_id}: observed={observed!r}, expected={expected!r}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--candidate-root", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    candidate = args.candidate_root.resolve()
    validation = candidate / "validation"
    if (validation / "READY").exists():
        raise SpatialResourceError(f"immutable validation candidate already exists: {validation}")
    checks: list[dict[str, object]] = []

    config = Path(__file__).with_name("spatial_resource_datasets.tsv")
    registry = load_dataset_registry(config)
    registry_by_id = {row["dataset_id"]: row for row in registry}
    registry_json = json.loads((candidate / "registry/spatial_dataset_registry.json").read_text(encoding="utf-8"))
    check(registry_json["release_id"] == RESOURCE_RELEASE_ID, "registry_release_id", registry_json["release_id"], RESOURCE_RELEASE_ID, checks)
    check(len(registry_json["datasets"]) == len(registry), "registry_complete", len(registry_json["datasets"]), len(registry), checks)

    hotspot = project / "Analysis/Multimodal_Program_Projection/candidates" / PROGRAM_RELEASE_ID / "hotspot"
    programs, _ = load_frozen_programs(hotspot)
    selected = {uid for uid, row in programs.items() if row["external_test_eligible"].upper() == "TRUE"}
    check(len(programs) == EXPECTED_N_PROGRAMS, "program_family_117", len(programs), EXPECTED_N_PROGRAMS, checks)
    check(len(selected) == EXPECTED_N_CONFIRMATORY, "confirmatory_family_two", len(selected), EXPECTED_N_CONFIRMATORY, checks)

    _, coverage = read_tsv(candidate / "coverage/spatial_program_coverage.tsv")
    require_complete_coverage(coverage, registry)
    required_assays = sum(parse_bool(row["coverage_required"], row["dataset_id"]) for row in registry)
    check(len(coverage) == EXPECTED_N_PROGRAMS * required_assays, "coverage_cartesian_complete", len(coverage), EXPECTED_N_PROGRAMS * required_assays, checks)
    check(not any(key in coverage[0] for key in ("pvalue", "qvalue", "validation", "winner")), "coverage_has_no_outcome_fields", sorted(coverage[0]), "no outcome fields", checks)
    check(all(row["registry_sha256"] == EXPECTED_REGISTRY_SHA256 for row in coverage), "coverage_registry_hash", "all_match", EXPECTED_REGISTRY_SHA256, checks)

    _, effects = read_tsv(candidate / "effects/spatial_program_effects.tsv")
    groups = defaultdict(list)
    for row in effects:
        groups[(row["dataset_id"], row["assay_id"])].append(row)
    check(all(len(rows) == 2 and {row["program_uid"] for row in rows} == selected for rows in groups.values()), "effects_complete_two_program_groups", len(groups), "all groups exactly two", checks)
    check(not any(row["pvalue"] or row["qvalue"] for row in groups[("HRA007511_HMSMA", "visium_rna")]), "hmsma_no_population_pvalues", "none", "none", checks)
    check(all(row["biological_unit_resolution"] == "unresolved" and not row["n_biological"] for row in groups[("HRA007511_HMSMA", "visium_rna")]), "hmsma_unresolved_unit", "unresolved/no donor count", "unresolved/no donor count", checks)
    check(all(row["source_dependence"] == "source_dependent" and row["biological_unit_resolution"] == "unresolved" for row in groups[("Vu_et_al_2025", "visium_rna")]), "vu_source_dependent", "source_dependent/unresolved", "source_dependent/unresolved", checks)
    check(all(row["biological_unit_resolution"] == "resolved" and row["n_biological"] == "4" for row in groups[("GSE192741", "visium_rna")]), "gse192741_donor_first", "4 resolved donors", "4 resolved donors", checks)
    coverage_by_key = {(row["dataset_id"], row["program_uid"]): row for row in coverage}
    compatible_match = True
    for dataset_id in ("GSE192741", "Vu_et_al_2025"):
        for effect in groups[(dataset_id, "visium_rna")]:
            observed = coverage_by_key[(dataset_id, effect["program_uid"])]
            compatible_match &= observed["n_genes_measured"] == effect["n_genes_measured"]
            compatible_match &= math.isclose(
                float(observed["retained_l1_weight"]),
                float(effect["retained_l1_weight"]),
                rel_tol=0,
                abs_tol=1e-12,
            )
    check(compatible_match, "coverage_matches_confirmatory_detection_gate", compatible_match, True, checks)

    gse_states = {row["program_uid"]: row["evidence_state"] for row in groups[("GSE192741", "visium_rna")]}
    igfbp7 = "hotspot_hepatocytes_f05c535ae5bbc0b9"
    bicc1 = "hotspot_hepatocytes_48f39dd4d817a10e"
    check(gse_states == {igfbp7: "supported", bicc1: "indeterminate"}, "asymmetric_primary_result", gse_states, "IGFBP7 supported; BICC1 indeterminate", checks)
    vu_rows = {row["program_uid"]: row for row in groups[("Vu_et_al_2025", "visium_rna")]}
    check(vu_rows[igfbp7]["within_source_result"] == "supported" and vu_rows[igfbp7]["evidence_state"] == "source_dependent", "vu_igfbp7_qualified_support", f"{vu_rows[igfbp7]['within_source_result']}/{vu_rows[igfbp7]['evidence_state']}", "supported/source_dependent", checks)
    check(vu_rows[bicc1]["evidence_state"] == "indeterminate", "vu_bicc1_indeterminate", vu_rows[bicc1]["evidence_state"], "indeterminate", checks)
    check(all(row["matched_null_sd"] == "" or row["uncertainty_semantics"] != "sampling_standard_error" for row in effects), "null_sd_not_sampling_se", "qualified", "never sampling SE", checks)
    check(registry_by_id["Govaere2026_GeoMx"]["dataset_gate"] == "dropped", "geomx_dropped", registry_by_id["Govaere2026_GeoMx"]["dataset_gate"], "dropped", checks)
    check(registry_by_id["GSE287826"]["biological_unit_resolution"] == "unresolved", "gse287826_no_donor_inference", registry_by_id["GSE287826"]["biological_unit_resolution"], "unresolved", checks)

    _, old_ready = read_tsv(
        project / "Analysis/Multimodal_Program_Projection/candidates" / PROGRAM_RELEASE_ID
        / "spatial_context_semantic_v2_2026-08-08/final_integration/READY"
    )
    check(old_ready[0]["v1_preservation_unchanged"] == "TRUE", "v1_numerical_regression_preserved", old_ready[0]["v1_preservation_unchanged"], "TRUE", checks)
    check(old_ready[0]["n_integrated_program_rows"] == "18", "sealed_semantic_source_complete", old_ready[0]["n_integrated_program_rows"], "18", checks)

    portal_manifest_fields, portal_manifest = read_tsv(candidate / "portal/product_manifest.tsv")
    check(all(row["release_id"] == RESOURCE_RELEASE_ID for row in portal_manifest), "portal_release_id_consistent", "all_match", RESOURCE_RELEASE_ID, checks)
    for row in portal_manifest:
        product = candidate / "portal" / row["product"]
        check(product.is_file() and sha256_file(product) == row["sha256"], f"portal_checksum_{row['product']}", row["sha256"], "file hash matches", checks)
    semantics = json.loads((candidate / "portal/spatial_semantics.json").read_text(encoding="utf-8"))
    check(semantics["display_contract"]["zero_is_negative_only_with_adequate_negative_rule"], "portal_zero_not_automatic_negative", True, True, checks)
    check(semantics["legacy_compatibility"]["excluded_from_rankings"], "legacy_s5_removed_from_ranking", True, True, checks)
    check(semantics["legacy_compatibility"]["excluded_from_active_layer_counts"], "legacy_s5_removed_from_active_layers", True, True, checks)

    program_parquet = pd.read_parquet(candidate / "portal/spatial_program_coverage.parquet")
    effect_parquet = pd.read_parquet(candidate / "portal/spatial_program_effects.parquet")
    check(set(program_parquet["release_id"]) == {RESOURCE_RELEASE_ID}, "coverage_parquet_release", set(program_parquet["release_id"]), {RESOURCE_RELEASE_ID}, checks)
    check(set(effect_parquet["release_id"]) == {RESOURCE_RELEASE_ID}, "effects_parquet_release", set(effect_parquet["release_id"]), {RESOURCE_RELEASE_ID}, checks)

    validation.mkdir(parents=True, exist_ok=True)
    write_tsv(validation / "validation_report.tsv", ("check_id", "status", "observed", "expected"), checks)
    write_tsv(
        validation / "READY",
        ("release_id", "status", "n_checks", "n_failed", "validation_sha256", "canonical_promotion_authorized"),
        [{
            "release_id": RESOURCE_RELEASE_ID,
            "status": "validated_candidate_awaiting_figure_and_portal_UI_adjudication",
            "n_checks": len(checks),
            "n_failed": Counter(row["status"] for row in checks)["FAIL"],
            "validation_sha256": sha256_file(validation / "validation_report.tsv"),
            "canonical_promotion_authorized": "FALSE",
        }],
    )
    print(f"PASS {len(checks)} spatial Resource candidate checks")


if __name__ == "__main__":
    main()
