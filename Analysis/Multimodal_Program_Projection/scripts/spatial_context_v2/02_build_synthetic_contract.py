#!/usr/bin/env python3
"""Materialize outcome-free fixtures for the common spatial-context contract."""

from __future__ import annotations

import argparse
import sys
from datetime import datetime, timezone
from pathlib import Path

from contract_lib import (
    ADAPTER_REGISTRY_COLUMNS,
    ALLOWED_EVIDENCE_STATES,
    ALLOWED_INTERVAL_TYPES,
    RELEASE_ID,
    SCHEMAS,
    canonical_id_hash,
    program_universe,
    sha256_file,
    validate_hotspot_seal,
    write_tsv,
)


SCHEMA_COLUMNS = ("artifact", "column", "required", "logical_type", "semantics")
STATE_COLUMNS = ("evidence_state", "ordered", "numeric_encoding_permitted", "semantics")
INTERVAL_COLUMNS = (
    "interval_type",
    "requires_paired_bounds",
    "implies_confidence_coverage",
    "semantics",
)


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[4]


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def logical_type(column: str) -> str:
    if column.startswith("n_") or column in {"bytes", "technical_unit_count"}:
        return "integer_or_empty"
    if column in {
        "estimate",
        "std_error",
        "matched_null_sd",
        "interval_low",
        "interval_high",
        "pvalue",
        "padj",
        "heterogeneity_statistic",
        "heterogeneity_df",
        "heterogeneity_pvalue",
        "retained_l1_weight",
        "program_score",
    }:
        return "finite_number_or_empty"
    if column in {
        "testable",
        "include_primary",
        "outcomes_tested",
        "sensitivity_sign_agree",
        "robustness_pass",
        "cross_assay_comparable",
        "direction_agreement",
        "sign_agree",
        "public_access",
    }:
        return "TRUE_or_FALSE"
    if column.endswith("sha256"):
        return "sha256_hex"
    return "string"


def semantics(column: str) -> str:
    special = {
        "estimate": "assay-native estimate; never cross-assay standardized",
        "effect_unit": "explicit assay-native effect unit",
        "std_error": "finite positive standard error when available; blank when no valid combined SE exists",
        "matched_null_sd": "matched-null dispersion only; never a sampling standard error or confidence interval",
        "interval_low": "lower bound whose meaning is declared by interval_type; never assumed to be a confidence limit",
        "interval_high": "upper bound whose meaning is declared by interval_type; never assumed to be a confidence limit",
        "interval_type": "explicit interval semantics from the closed allowed_interval_types vocabulary",
        "pvalue_method": "explicit native method used to derive pvalue",
        "evidence_state": "categorical state; never ordered, scored, or summed",
        "descriptive_effect_direction": "direction of the assay-native descriptive estimate",
        "inferential_test_direction": "direction of the statistic that generated the inferential P value",
        "direction_agreement": "whether descriptive and inferential directions agree",
        "negative_call_rule_id": "required only for tested_negative; identifies a prespecified adequate-negative rule",
        "n_biological": "unique inferential biological units in sample_manifest",
        "n_technical": "unique included technical units in sample_manifest",
        "registry_sha256": "file hash anchored by the Hotspot READY seal",
        "membership_sha256": "content hash for the exact frozen program definition",
        "cross_assay_comparable": "must be FALSE for this contract",
    }
    return special.get(column, "required adapter provenance or assay-native field")


def write_contract_metadata(candidate_root: Path) -> None:
    schema_rows = []
    for artifact, columns in SCHEMAS.items():
        for column in columns:
            schema_rows.append(
                {
                    "artifact": artifact,
                    "column": column,
                    "required": True,
                    "logical_type": logical_type(column),
                    "semantics": semantics(column),
                }
            )
    write_tsv(candidate_root / "contract_schema.tsv", SCHEMA_COLUMNS, schema_rows)
    state_semantics = {
        "robust": "native statistical rule and all prespecified sensitivity gates pass",
        "indeterminate": "valid native test performed, but neither the robust rule nor an explicit adequate-negative rule passes",
        "tested_negative": "a prespecified adequate-negative criterion passes and negative_call_rule_id identifies that rule",
        "untestable": "coverage or design prevents a valid test",
        "not_applicable": "assay cannot observe the lineage or scientific question",
        "skipped": "public-source or donor gate terminates the dataset before outcomes",
        "source_dependent": "informative context that is not independent validation",
    }
    write_tsv(
        candidate_root / "allowed_evidence_states.tsv",
        STATE_COLUMNS,
        [
            {
                "evidence_state": state,
                "ordered": False,
                "numeric_encoding_permitted": False,
                "semantics": state_semantics[state],
            }
            for state in sorted(ALLOWED_EVIDENCE_STATES)
        ],
    )
    interval_semantics = {
        "confidence_interval": "paired finite statistical confidence limits supplied by the native producer",
        "donor_effect_range": "paired finite minimum-to-maximum donor effects; descriptive heterogeneity only, with no confidence-coverage implication",
        "none": "no interval is supplied; a tested effect requires either a finite positive sampling standard error or an explicitly separate matched_null_sd",
        "not_applicable": "terminal untested state for which interval uncertainty is not applicable",
    }
    write_tsv(
        candidate_root / "allowed_interval_types.tsv",
        INTERVAL_COLUMNS,
        [
            {
                "interval_type": interval_type,
                "requires_paired_bounds": interval_type in {"confidence_interval", "donor_effect_range"},
                "implies_confidence_coverage": interval_type == "confidence_interval",
                "semantics": interval_semantics[interval_type],
            }
            for interval_type in sorted(ALLOWED_INTERVAL_TYPES)
        ],
    )


def source_manifest(adapter_root: Path, retrieved: str) -> str:
    source_path = adapter_root / "synthetic_source.txt"
    source_path.write_text(
        "SYNTHETIC CONTRACT FIXTURE ONLY. No biological or external outcome is represented.\n",
        encoding="utf-8",
    )
    write_tsv(
        adapter_root / "source_manifest.tsv",
        SCHEMAS["source_manifest.tsv"],
        [
            {
                "path_scope": "adapter_relative",
                "relative_path": "synthetic_source.txt",
                "bytes": source_path.stat().st_size,
                "sha256": sha256_file(source_path),
                "source_role": "synthetic_contract_fixture",
                "public_access": False,
                "retrieved_utc": retrieved,
            }
        ],
    )
    return sha256_file(adapter_root / "source_manifest.tsv")


def execution_manifest(adapter_root: Path, created: str) -> str:
    filenames = (
        "sample_manifest.tsv",
        "gene_mapping_audit.tsv",
        "design_audit.tsv",
        "program_testability.tsv",
        "per_sample_program_scores.tsv",
        "program_effects.tsv",
        "sensitivity.tsv",
        "source_manifest.tsv",
    )
    rows = []
    for filename in filenames:
        path = adapter_root / filename
        rows.append(
            {
                "role": "synthetic_adapter_artifact",
                "relative_path": filename,
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "created_utc": created,
            }
        )
    write_tsv(adapter_root / "execution_manifest.tsv", SCHEMAS["execution_manifest.tsv"], rows)
    return sha256_file(adapter_root / "execution_manifest.tsv")


def adapter_specifications(seal) -> list[dict[str, object]]:
    common = {
        "source_publication": "SYNTHETIC_FIXTURE_NOT_A_PUBLICATION",
        "source_dependence": "independent",
        "biological_unit": "donor",
        "biological_unit_resolution": "resolved",
        "analysis_set_id": "synthetic_primary",
        "cross_assay_comparable": False,
        "program_universe": "robust_display",
        "min_genes_testable": 8,
        "min_retained_l1_weight": 0.2,
        "registry_sha256": seal.registry_sha256,
        "ready_sha256": seal.ready_sha256,
        "status": "synthetic_fixture_only",
    }
    return [
        {
            **common,
            "adapter_id": "synthetic_visium_pass",
            "adapter_root": "synthetic_fixture/valid/synthetic_visium_pass",
            "dataset": "SYNTHETIC_VISIUM",
            "assay": "Visium",
            "technical_unit": "section",
            "contrast_or_exposure": "synthetic_spatial_organization",
            "effect_unit": "matched_null_residual_moran_z",
            "gate_expectation": "pass",
        },
        {
            **common,
            "adapter_id": "synthetic_lipid_pass",
            "adapter_root": "synthetic_fixture/valid/synthetic_lipid_pass",
            "dataset": "SYNTHETIC_LIPID",
            "assay": "Visium_lipid_context",
            "technical_unit": "section",
            "contrast_or_exposure": "source_defined_lipid_zone_vs_non_lipid_zone",
            "effect_unit": "source_lipid_zone_minus_non_lipid_zone_score",
            "gate_expectation": "pass",
        },
        {
            **common,
            "adapter_id": "synthetic_geomx_skip",
            "adapter_root": "synthetic_fixture/valid/synthetic_geomx_skip",
            "dataset": "SYNTHETIC_GEOMX_SKIPPED",
            "assay": "GeoMx",
            "technical_unit": "AOI",
            "contrast_or_exposure": "synthetic_MASH_minus_healthy",
            "effect_unit": "donor_score_difference",
            "gate_expectation": "skipped",
        },
    ]


def sample_rows(spec: dict[str, object]) -> list[dict[str, object]]:
    if spec["gate_expectation"] == "skipped":
        return []
    if spec["adapter_id"] == "synthetic_visium_pass":
        ids = (("V_D1", "V_S1"), ("V_D1", "V_S2"), ("V_D2", "V_S3"))
    else:
        ids = (("L_D1", "L_S1"), ("L_D2", "L_S2"), ("L_D3", "L_S3"))
    return [
        {
            "release_id": RELEASE_ID,
            "dataset": spec["dataset"],
            "analysis_set_id": spec["analysis_set_id"],
            "biological_id": biological,
            "technical_id": technical,
            "biological_unit": spec["biological_unit"],
            "technical_unit": spec["technical_unit"],
            "include_primary": True,
            "gate_state": "included",
        }
        for biological, technical in ids
    ]


def build_adapter(project_root: Path, candidate_root: Path, spec: dict[str, object], seal, created: str) -> None:
    adapter_root = candidate_root / str(spec["adapter_root"])
    adapter_root.mkdir(parents=True, exist_ok=True)
    source_manifest_hash = source_manifest(adapter_root, created)
    universe = program_universe(seal, str(spec["program_universe"]))
    programs = [universe[uid] for uid in sorted(universe)]
    samples = sample_rows(spec)
    biological_ids = sorted({str(row["biological_id"]) for row in samples})
    technical_ids = sorted(str(row["technical_id"]) for row in samples)
    write_tsv(adapter_root / "sample_manifest.tsv", SCHEMAS["sample_manifest.tsv"], samples)

    is_skip = spec["gate_expectation"] == "skipped"
    design = {
        "release_id": RELEASE_ID,
        "dataset": spec["dataset"],
        "analysis_set_id": spec["analysis_set_id"],
        "contrast_or_exposure": spec["contrast_or_exposure"],
        "biological_unit": spec["biological_unit"],
        "biological_unit_resolution": spec["biological_unit_resolution"],
        "n_biological": len(biological_ids),
        "technical_unit": spec["technical_unit"],
        "n_technical": len(technical_ids),
        "design_status": "skipped_source_gate" if is_skip else "passed_synthetic_design",
        "biological_ids_sha256": canonical_id_hash(biological_ids),
        "technical_ids_sha256": canonical_id_hash(technical_ids),
    }
    write_tsv(adapter_root / "design_audit.tsv", SCHEMAS["design_audit.tsv"], [design])

    mapping_rows = []
    testability_rows = []
    score_rows = []
    effect_rows = []
    sensitivity_rows = []
    producer_relative = Path(__file__).resolve().relative_to(project_root).as_posix()
    producer_hash = sha256_file(Path(__file__).resolve())
    technical_count_by_donor = {
        biological: sum(1 for row in samples if row["biological_id"] == biological)
        for biological in biological_ids
    }

    for index, reg in enumerate(programs):
        uid = reg["program_uid"]
        testable = not is_skip
        n_genes = 0 if is_skip else 12 + index
        retained = 0.0 if is_skip else 0.35 + 0.05 * index
        state = "skipped" if is_skip else ("robust" if index == 0 else "indeterminate")
        mapping_rows.append(
            {
                "release_id": RELEASE_ID,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "program_uid": uid,
                "membership_sha256": reg["membership_sha256"],
                "n_source_genes": reg["n_source_genes"],
                "n_genes_measured": n_genes,
                "retained_l1_weight": retained,
                "mapping_status": "skipped_source_gate" if is_skip else "passed_synthetic_coverage",
            }
        )
        testability_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_sha256": seal.registry_sha256,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "program_uid": uid,
                "membership_sha256": reg["membership_sha256"],
                "testable": testable,
                "n_genes_measured": n_genes,
                "retained_l1_weight": retained,
                "testability_reason": "skipped_source_gate" if is_skip else "passed_synthetic_coverage",
                "evidence_state": state,
            }
        )
        if testable:
            for donor_index, biological_id in enumerate(biological_ids):
                sign = 1.0 if index == 0 else -1.0
                score_rows.append(
                    {
                        "release_id": RELEASE_ID,
                        "registry_sha256": seal.registry_sha256,
                        "dataset": spec["dataset"],
                        "analysis_set_id": spec["analysis_set_id"],
                        "program_uid": uid,
                        "membership_sha256": reg["membership_sha256"],
                        "biological_id": biological_id,
                        "technical_unit_count": technical_count_by_donor[biological_id],
                        "program_score": f"{sign * (donor_index + 1) / 10:.6f}",
                        "score_unit": "frozen_weighted_program_score",
                    }
                )

        if is_skip:
            estimate = se = interval_low = interval_high = pvalue = padj = ""
            interval_type = "not_applicable"
            pvalue_method = "not_applicable"
            direction = ""
            robustness = False
            sign_agree = False
        elif index == 0:
            estimate, interval_low, interval_high = 0.5, 0.2, 0.8
            pvalue, padj = 0.01, 0.02
            if spec["adapter_id"] == "synthetic_lipid_pass":
                se = ""
                interval_type = "donor_effect_range"
                pvalue_method = "equal_donor_signed_stouffer"
            else:
                se = 0.15
                interval_type = "confidence_interval"
                pvalue_method = "synthetic_native_spatial_test"
            direction = "positive"
            robustness = True
            sign_agree = True
        else:
            estimate, interval_low, interval_high = -0.1, -0.6, 0.4
            pvalue, padj = 0.6, 0.6
            if spec["adapter_id"] == "synthetic_lipid_pass":
                se = ""
                interval_type = "donor_effect_range"
                pvalue_method = "equal_donor_signed_stouffer"
            else:
                se = 0.25
                interval_type = "confidence_interval"
                pvalue_method = "synthetic_native_spatial_test"
            direction = "negative"
            robustness = False
            sign_agree = True
        effect_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_sha256": seal.registry_sha256,
                "membership_sha256": reg["membership_sha256"],
                "program_uid": uid,
                "legacy_program_id": f"{reg['cell_type']}:{reg['module']}",
                "program_label": reg["module_name"],
                "cell_type": reg["cell_type"],
                "dataset": spec["dataset"],
                "assay": spec["assay"],
                "source_publication": spec["source_publication"],
                "source_dependence": spec["source_dependence"],
                "analysis_set_id": spec["analysis_set_id"],
                "biological_unit": spec["biological_unit"],
                "biological_unit_resolution": spec["biological_unit_resolution"],
                "n_biological": len(biological_ids),
                "technical_unit": spec["technical_unit"],
                "n_technical": len(technical_ids),
                "contrast_or_exposure": spec["contrast_or_exposure"],
                "effect_unit": spec["effect_unit"],
                "estimate": estimate,
                "std_error": se,
                "matched_null_sd": "",
                "interval_low": interval_low,
                "interval_high": interval_high,
                "interval_type": interval_type,
                "pvalue": pvalue,
                "padj": padj,
                "pvalue_method": pvalue_method,
                "multiplicity_family": "not_tested" if is_skip else "synthetic_complete_family_n2",
                "n_genes_measured": n_genes,
                "retained_l1_weight": retained,
                "testable": testable,
                "testability_reason": "skipped_source_gate" if is_skip else "passed_synthetic_coverage",
                "direction_expected": reg.get("primary_direction", "not_prespecified"),
                "direction_observed": direction,
                "descriptive_effect_direction": direction,
                "inferential_test_direction": direction,
                "direction_agreement": False if is_skip else True,
                "heterogeneity_statistic": "",
                "heterogeneity_df": "",
                "heterogeneity_pvalue": "",
                "sensitivity_sign_agree": sign_agree,
                "robustness_pass": robustness,
                "negative_call_rule_id": "",
                "cross_assay_comparable": False,
                "evidence_state": state,
                "producer": producer_relative,
                "producer_sha256": producer_hash,
                "source_manifest_sha256": source_manifest_hash,
            }
        )
        sensitivity_rows.append(
            {
                "release_id": RELEASE_ID,
                "registry_sha256": seal.registry_sha256,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "program_uid": uid,
                "membership_sha256": reg["membership_sha256"],
                "sensitivity_id": "synthetic_equal_weight",
                "estimate": "" if is_skip else estimate,
                "effect_unit": spec["effect_unit"],
                "direction_observed": direction,
                "sign_agree": False if is_skip else True,
                "status": "skipped" if is_skip else "synthetic_pass",
            }
        )

    write_tsv(adapter_root / "gene_mapping_audit.tsv", SCHEMAS["gene_mapping_audit.tsv"], mapping_rows)
    write_tsv(adapter_root / "program_testability.tsv", SCHEMAS["program_testability.tsv"], testability_rows)
    write_tsv(adapter_root / "per_sample_program_scores.tsv", SCHEMAS["per_sample_program_scores.tsv"], score_rows)
    write_tsv(adapter_root / "program_effects.tsv", SCHEMAS["program_effects.tsv"], effect_rows)
    write_tsv(adapter_root / "sensitivity.tsv", SCHEMAS["sensitivity.tsv"], sensitivity_rows)
    execution_hash = execution_manifest(adapter_root, created)
    write_tsv(
        adapter_root / "gate_status.tsv",
        SCHEMAS["gate_status.tsv"],
        [
            {
                "release_id": RELEASE_ID,
                "dataset": spec["dataset"],
                "analysis_set_id": spec["analysis_set_id"],
                "gate_status": spec["gate_expectation"],
                "gate_reason": "synthetic_fixture_skip" if is_skip else "synthetic_fixture_pass",
                "outcomes_tested": not is_skip,
                "registry_sha256": seal.registry_sha256,
                "source_manifest_sha256": source_manifest_hash,
                "execution_manifest_sha256": execution_hash,
            }
        ],
    )


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=project_root_from_script())
    parser.add_argument("--candidate-root", type=Path, default=None)
    parser.add_argument("--hotspot-root", type=Path, default=None)
    args = parser.parse_args()
    project_root = args.project_root.resolve()
    candidate_root = args.candidate_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "spatial_context"
    )
    hotspot_root = args.hotspot_root or (
        project_root / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID / "hotspot"
    )
    try:
        seal = validate_hotspot_seal(hotspot_root)
        candidate_root.mkdir(parents=True, exist_ok=True)
        write_contract_metadata(candidate_root)
        specs = adapter_specifications(seal)
        write_tsv(candidate_root / "adapter_registry.tsv", ADAPTER_REGISTRY_COLUMNS, specs)
        created = utc_now()
        for spec in specs:
            build_adapter(project_root, candidate_root, spec, seal, created)
        print(
            f"built {len(specs)} synthetic-only adapters for "
            f"{sum(1 for row in seal.registry_rows if row['robust_display'] == 'TRUE')} sealed programs"
        )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
