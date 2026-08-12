#!/usr/bin/env python3
"""Integrate the sealed two-program family and HMSMA descriptive branch."""

from __future__ import annotations

import argparse
import statistics
from collections import defaultdict
from pathlib import Path

import pandas as pd

from spatial_resource_lib import (
    EFFECT_COLUMNS,
    EXPECTED_N_CONFIRMATORY,
    EXPECTED_REGISTRY_SHA256,
    PROGRAM_RELEASE_ID,
    RESOURCE_RELEASE_ID,
    SpatialResourceError,
    canonical_row_sha256,
    load_dataset_registry,
    load_frozen_programs,
    read_tsv,
    sha256_file,
    write_tsv,
)


DATASET_ALIASES = {"Yakubovsky_2026": "Yakubovsky2026"}
ASSAY_IDS = {
    "GSE192741": "visium_rna",
    "Vu_et_al_2025": "visium_rna",
    "Yakubovsky2026": "spatial_transcriptomics_rna",
    "Govaere2026_CosMx": "cosmx_targeted_rna",
    "Govaere2026_GeoMx": "geomx_targeted_rna",
    "GSE287826": "geomx_targeted_rna",
    "GSE244832": "snatac_promoter_accessibility",
    "GSE281367": "snatac_promoter_accessibility",
    "PXD051911": "liver_dia_ms",
}


def blank_if_na(value: str) -> str:
    return "" if value.strip().lower() in {"", "na", "nan", "none"} else value


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--candidate-root", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    candidate = args.candidate_root.resolve()
    output = candidate / "effects"
    if (output / "READY").exists():
        raise SpatialResourceError(f"immutable effects candidate already exists: {output}")

    config = Path(__file__).with_name("spatial_resource_datasets.tsv")
    registry_rows = load_dataset_registry(config)
    registry = {(row["dataset_id"], row["assay_id"]): row for row in registry_rows}
    hotspot = project / "Analysis/Multimodal_Program_Projection/candidates" / PROGRAM_RELEASE_ID / "hotspot"
    programs, _ = load_frozen_programs(hotspot)
    selected = {uid: row for uid, row in programs.items() if row["external_test_eligible"].upper() == "TRUE"}

    semantic_root = (
        project / "Analysis/Multimodal_Program_Projection/candidates" / PROGRAM_RELEASE_ID
        / "spatial_context_semantic_v2_2026-08-08"
    )
    source_path = semantic_root / "final_integration/integrated_program_effects.tsv"
    source_fields, source_rows = read_tsv(source_path, ("program_uid", "dataset", "evidence_state"))
    source_rows = [row for row in source_rows if row["program_uid"] in selected]

    _, sensitivity_rows = read_tsv(
        semantic_root / "native_spatial/v2_candidate/spatial_sensitivity_audit.tsv",
        ("dataset", "program_id", "sensitivity_id", "sign_agree_with_primary"),
    )
    sensitivities = defaultdict(dict)
    for row in sensitivity_rows:
        sensitivities[(row["dataset"], row["program_id"])][row["sensitivity_id"]] = row["sign_agree_with_primary"]

    effects = []
    for source in source_rows:
        dataset_id = DATASET_ALIASES.get(source["dataset"], source["dataset"])
        assay_id = ASSAY_IDS[dataset_id]
        contract = registry[(dataset_id, assay_id)]
        underlying = "supported" if source["evidence_state"] == "robust" else source["evidence_state"]
        if contract["dataset_gate"] == "dropped":
            state = "dropped"
        elif dataset_id == "GSE287826":
            state = "untestable"
        elif underlying == "supported" and contract["source_dependence"] == "source_dependent":
            state = "source_dependent"
        elif underlying in {"skipped", "untestable", "not_applicable", "indeterminate", "supported"}:
            state = underlying
        else:
            raise SpatialResourceError(f"unmapped evidence state: {source['evidence_state']}")
        sensitivity = sensitivities.get((source["dataset"], source["program_uid"]), {})
        effects.append({
            "release_id": RESOURCE_RELEASE_ID,
            "program_release_id": PROGRAM_RELEASE_ID,
            "registry_sha256": EXPECTED_REGISTRY_SHA256,
            "dataset_id": dataset_id,
            "assay_id": assay_id,
            "program_uid": source["program_uid"],
            "program_label": source["program_label"],
            "membership_sha256": source["membership_sha256"],
            "biological_unit": contract["biological_unit"],
            "biological_unit_resolution": contract["biological_unit_resolution"],
            "n_biological": contract["n_biological"],
            "technical_unit": contract["technical_unit"],
            "n_technical": contract["n_technical"],
            "source_dependence": contract["source_dependence"],
            "dataset_gate": contract["dataset_gate"],
            "estimand": source["contrast_or_exposure"],
            "effect_unit": source["effect_unit"],
            "estimate": blank_if_na(source["estimate"]),
            "matched_null_sd": blank_if_na(source["matched_null_sd"]),
            "pvalue": blank_if_na(source["pvalue"]),
            "qvalue": blank_if_na(source["padj"]),
            "multiplicity_family": source["multiplicity_family"],
            "n_null_draws": "9999" if "9999" in source["pvalue_method"] else "",
            "n_genes_measured": source["n_genes_measured"],
            "retained_l1_weight": source["retained_l1_weight"],
            "equal_weight_sign_agree": sensitivity.get("equal_weight", source["sensitivity_sign_agree"]),
            "leave_top_weighted_gene_sign_agree": sensitivity.get("leave_highest_weight_gene_out", source["sensitivity_sign_agree"]),
            "within_source_result": underlying,
            "evidence_state": state,
            "testability_reason": source["testability_reason"],
            "uncertainty_semantics": source["uncertainty_semantics"],
            "source_release_id": source["release_id"],
            "source_row_sha256": canonical_row_sha256(source, source_fields),
        })

    _, hmsma_rows = read_tsv(
        candidate / "hmsma_label_blind/per_array_program_organization.tsv",
        ("program_uid", "residual_moran_i", "n_genes_measured", "retained_l1_weight"),
    )
    by_program = defaultdict(list)
    for row in hmsma_rows:
        by_program[row["program_uid"]].append(row)
    hmsma_contract = registry[("HRA007511_HMSMA", "visium_rna")]
    _, coverage_rows = read_tsv(
        candidate / "coverage/spatial_program_coverage.tsv",
        ("dataset_id", "program_uid", "n_genes_measured", "retained_l1_weight"),
    )
    hmsma_coverage = {
        row["program_uid"]: row
        for row in coverage_rows
        if row["dataset_id"] == "HRA007511_HMSMA"
    }
    for uid, program in selected.items():
        part = by_program[uid]
        if len(part) != 35:
            raise SpatialResourceError(f"HMSMA {uid} expected 35 descriptive array rows, found {len(part)}")
        values = [float(row["residual_moran_i"]) for row in part]
        effects.append({
            "release_id": RESOURCE_RELEASE_ID,
            "program_release_id": PROGRAM_RELEASE_ID,
            "registry_sha256": EXPECTED_REGISTRY_SHA256,
            "dataset_id": "HRA007511_HMSMA",
            "assay_id": "visium_rna",
            "program_uid": uid,
            "program_label": program["module_name"],
            "membership_sha256": program["membership_sha256"],
            "biological_unit": hmsma_contract["biological_unit"],
            "biological_unit_resolution": "unresolved",
            "n_biological": "",
            "technical_unit": "physical_array",
            "n_technical": "35",
            "source_dependence": "source_dependent",
            "dataset_gate": "metadata_pending",
            "estimand": "median_within_array_library_and_detection_adjusted_spatial_organization_descriptive_only",
            "effect_unit": "median_array_residual_moran_i",
            "estimate": f"{statistics.median(values):.17g}",
            "matched_null_sd": "",
            "pvalue": "",
            "qvalue": "",
            "multiplicity_family": "not_applicable_no_population_inference",
            "n_null_draws": "",
            "n_genes_measured": hmsma_coverage[uid]["n_genes_measured"],
            "retained_l1_weight": hmsma_coverage[uid]["retained_l1_weight"],
            "equal_weight_sign_agree": "",
            "leave_top_weighted_gene_sign_agree": "",
            "within_source_result": "descriptive_only",
            "evidence_state": "source_dependent",
            "testability_reason": "array_level_descriptive_score_available_but_donor_identity_unresolved",
            "uncertainty_semantics": "array_distribution_is_not_a_sampling_distribution_over_independent_donors",
            "source_release_id": RESOURCE_RELEASE_ID,
            "source_row_sha256": sha256_file(candidate / "hmsma_label_blind/per_array_program_organization.tsv"),
        })

    groups = defaultdict(list)
    for row in effects:
        groups[(row["dataset_id"], row["assay_id"])].append(row["program_uid"])
    for key, uids in groups.items():
        if len(uids) != EXPECTED_N_CONFIRMATORY or set(uids) != set(selected):
            raise SpatialResourceError(f"incomplete two-program family for {key}: {uids}")
    output.mkdir(parents=True, exist_ok=True)
    write_tsv(output / "spatial_program_effects.tsv", EFFECT_COLUMNS, effects)
    pd.DataFrame(effects, columns=EFFECT_COLUMNS).to_parquet(output / "spatial_program_effects.parquet", index=False)
    write_tsv(
        output / "READY",
        ("release_id", "status", "n_programs", "n_dataset_assay_groups", "n_effect_rows", "effects_sha256", "source_effects_sha256"),
        [{
            "release_id": RESOURCE_RELEASE_ID,
            "status": "complete_two_program_family_candidate_no_canonical_promotion",
            "n_programs": len(selected),
            "n_dataset_assay_groups": len(groups),
            "n_effect_rows": len(effects),
            "effects_sha256": sha256_file(output / "spatial_program_effects.parquet"),
            "source_effects_sha256": sha256_file(source_path),
        }],
    )
    print(f"WROTE {len(effects)} effects across {len(groups)} dataset/assay groups")


if __name__ == "__main__":
    main()
