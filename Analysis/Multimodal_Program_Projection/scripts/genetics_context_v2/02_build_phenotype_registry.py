#!/usr/bin/env python3
"""Build the prespecified 35-GWAS phenotype and source-provenance registry."""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

from genetics_common import (
    ContractError,
    DEFAULT_CANDIDATE_ROOT,
    PROJECT_ROOT,
    SCRIPT_DIR,
    assert_candidate_root,
    atomic_write_tsv,
    clean,
    load_config_tsv,
    read_table,
    write_stage_seal,
)


ALLOWED_STRATA = {
    "direct_masld_mash_diagnosis",
    "mri_pdff_or_histologic_steatosis",
    "ct1_fibrosis_or_longitudinal_progression",
    "alt_ast_or_ggt",
    "systemic_metabolic_trait",
}
EXPECTED_COUNTS = {
    "direct_masld_mash_diagnosis": 12,
    "mri_pdff_or_histologic_steatosis": 3,
    "alt_ast_or_ggt": 20,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--candidate-root", type=Path, default=DEFAULT_CANDIDATE_ROOT)
    parser.add_argument(
        "--catalog",
        type=Path,
        default=SCRIPT_DIR / "config" / "phenotype_source_catalog_v1.tsv",
    )
    return parser.parse_args()


def input_from_manifest(candidate: Path, input_id: str) -> Path:
    rows = read_table(candidate / "input_manifest.tsv", "tsv")
    matching = [row for row in rows if row["input_id"] == input_id]
    if len(matching) != 1:
        raise ContractError(f"expected one input-manifest row for {input_id}")
    row = matching[0]
    return candidate / row["snapshot_path"] if row["snapshot_path"] else Path(row["source_path"])


def main() -> None:
    args = parse_args()
    project = args.project_root.resolve()
    candidate = assert_candidate_root(project, args.candidate_root)
    stage1 = candidate / "work" / "stage_seals" / "01_freeze_and_rederive.json"
    if not stage1.is_file():
        raise ContractError("GEN-00 seal is missing; run 01_freeze_and_rederive.py first")

    registry = read_table(input_from_manifest(candidate, "gwas_registry"), "tsv")
    tier = read_table(input_from_manifest(candidate, "gwas_trait_tier"), "tsv")
    catalog = load_config_tsv(args.catalog.resolve())

    registry_by_name = {row["study_name"]: row for row in registry}
    if len(registry_by_name) != len(registry):
        raise ContractError("duplicate study_name in GWAS registry")
    tier_main = [row for row in tier if row["placement"] == "main"]
    tier_by_name = {row["study_name"]: row for row in tier_main}
    catalog_by_name = {row["study_name"]: row for row in catalog}
    if len(tier_main) != 35 or len(tier_by_name) != 35:
        raise ContractError(f"expected 35 unique primary studies; observed {len(tier_by_name)}")
    if len(catalog) != 35 or len(catalog_by_name) != 35:
        raise ContractError(f"phenotype catalog must contain exactly 35 unique studies")
    if set(tier_by_name) != set(catalog_by_name):
        raise ContractError(
            "catalog/main-registry mismatch: "
            f"missing={sorted(set(tier_by_name) - set(catalog_by_name))}; "
            f"extra={sorted(set(catalog_by_name) - set(tier_by_name))}"
        )
    if not set(tier_by_name).issubset(registry_by_name):
        raise ContractError("one or more primary studies are absent from gwas_registry.tsv")

    output_rows: list[dict[str, object]] = []
    for study_name in sorted(tier_by_name):
        reg = registry_by_name[study_name]
        trait = tier_by_name[study_name]
        source = catalog_by_name[study_name]
        stratum = source["phenotype_stratum"]
        if stratum not in ALLOWED_STRATA:
            raise ContractError(f"invalid phenotype stratum for {study_name}: {stratum}")
        if not all(clean(source[column]) for column in ("phenotype_definition", "source_citation", "source_url")):
            raise ContractError(f"incomplete source provenance for {study_name}")
        if source["source_ancestry"] != reg["ancestry"]:
            raise ContractError(
                f"source/registry ancestry mismatch for {study_name}: "
                f"{source['source_ancestry']} vs {reg['ancestry']}"
            )
        design_ok = (
            reg["trait_type"] == "binary" and source["source_design"] == "case_control"
        ) or (
            reg["trait_type"] == "quantitative" and source["source_design"] == "quantitative"
        )
        if not design_ok:
            raise ContractError(f"source/registry design mismatch for {study_name}")
        source_n_text = clean(source["source_sample_size"])
        registry_n = int(reg["N_tot"])
        if source_n_text.startswith("~"):
            source_n = int(source_n_text[1:])
            n_status = (
                "rounded_source_consistent_with_registry"
                if abs(source_n - registry_n) / max(source_n, 1) < 0.01
                else "rounded_source_registry_mismatch"
            )
        elif source_n_text.isdigit():
            source_n = int(source_n_text)
            n_status = (
                "source_consistent_with_registry"
                if source_n == registry_n
                else "source_registry_mismatch"
            )
        else:
            n_status = "source_sample_size_not_resolved"
        eqtl_ancestry_status = (
            "ancestry_matched_eur"
            if reg["ancestry"] == "EUR"
            else "cross_ancestry_eqtl_limited"
        )
        output_rows.append(
            {
                "study_name": study_name,
                "trait": trait["trait"],
                "tier": trait["tier"],
                "placement": trait["placement"],
                "phenotype_stratum": stratum,
                "phenotype_definition": source["phenotype_definition"],
                "direct_or_proxy": source["direct_or_proxy"],
                "design": reg["trait_type"],
                "ancestry": reg["ancestry"],
                "registry_N_total": reg["N_tot"],
                "registry_N_cases": reg["N_cases"],
                "source_sample_size": source["source_sample_size"],
                "sample_size_audit": n_status,
                "source_genome_build": source["source_build"],
                "analysis_genome_build": "hg19",
                "eqtl_panel": "Broadaway_liver_meta_N1183_EUR_hg19",
                "ld_panel": reg["ld_panel"],
                "ld_panel_n": reg["ld_panel_n"],
                "gwas_ld_ancestry_status": "registry_ancestry_matched",
                "regulatory_ancestry_status": eqtl_ancestry_status,
                "source_citation": source["source_citation"],
                "source_url": source["source_url"],
                "source_verification": source["source_verification"],
                "source_note": source["source_note"],
                "source_dependence": "partially_dependent_shared_broadaway_eqtl_panel",
                "bulk_discovery_overlap_status": "not_resolved_preflight",
            }
        )

    observed_counts = Counter(row["phenotype_stratum"] for row in output_rows)
    if dict(observed_counts) != EXPECTED_COUNTS:
        raise ContractError(
            f"primary phenotype composition drift: expected {EXPECTED_COUNTS}, "
            f"observed {dict(observed_counts)}"
        )

    registry_path = candidate / "phenotype_registry.tsv"
    registry_fields = [
        "study_name",
        "trait",
        "tier",
        "placement",
        "phenotype_stratum",
        "phenotype_definition",
        "direct_or_proxy",
        "design",
        "ancestry",
        "registry_N_total",
        "registry_N_cases",
        "source_sample_size",
        "sample_size_audit",
        "source_genome_build",
        "analysis_genome_build",
        "eqtl_panel",
        "ld_panel",
        "ld_panel_n",
        "gwas_ld_ancestry_status",
        "regulatory_ancestry_status",
        "source_citation",
        "source_url",
        "source_verification",
        "source_note",
        "source_dependence",
        "bulk_discovery_overlap_status",
    ]
    atomic_write_tsv(registry_path, output_rows, registry_fields)

    provenance_rows = [
        {
            "result_id": row["study_name"],
            "result_layer": "primary_gwas_coloc",
            "phenotype_or_panel": row["phenotype_definition"],
            "tissue_or_context": "GWAS trait joined to static bulk liver eQTL",
            "ancestry": row["ancestry"],
            "cohort_or_panel": row["source_citation"],
            "tested_universe": "canonical gene-level COLOC aggregation; not a source eQTL tested-gene denominator",
            "source_publication": row["source_citation"],
            "source_url": row["source_url"],
            "source_dependence": row["source_dependence"],
            "reuse_note": row["source_note"],
        }
        for row in output_rows
    ]
    provenance_rows.append(
        {
            "result_id": "Broadaway_liver_eQTL_meta",
            "result_layer": "static_liver_eqtl",
            "phenotype_or_panel": "liver cis-eQTL meta-analysis",
            "tissue_or_context": "bulk liver; static cis regulation",
            "ancestry": "EUR-prioritized",
            "cohort_or_panel": "four liver microarray cohorts; up to N=1183",
            "tested_universe": "18,322 expression genes reported in paper; complete deposited tested-gene table not established",
            "source_publication": "Broadaway et al. AJHG 2024; doi:10.1016/j.ajhg.2024.07.017",
            "source_url": "https://doi.org/10.1016/j.ajhg.2024.07.017",
            "source_dependence": "reused_source_derived_for_every_static_coloc",
            "reuse_note": "The same eQTL panel is reused across all 35 GWAS; GWAS rows are not independent eQTL experiments.",
        }
    )
    provenance_path = candidate / "source_provenance.tsv"
    atomic_write_tsv(
        provenance_path,
        provenance_rows,
        [
            "result_id",
            "result_layer",
            "phenotype_or_panel",
            "tissue_or_context",
            "ancestry",
            "cohort_or_panel",
            "tested_universe",
            "source_publication",
            "source_url",
            "source_dependence",
            "reuse_note",
        ],
    )

    audit_rows = [
        {
            "audit": "primary_study_count",
            "expected": 35,
            "observed": len(output_rows),
            "status": "pass",
            "note": "exactly one phenotype/source row per primary GWAS",
        },
        *[
            {
                "audit": f"stratum_count:{stratum}",
                "expected": expected,
                "observed": observed_counts[stratum],
                "status": "pass",
                "note": "prespecified before phenotype outcome comparison",
            }
            for stratum, expected in EXPECTED_COUNTS.items()
        ],
        {
            "audit": "non_eur_eqtl_limitation",
            "expected": sum(row["ancestry"] != "EUR" for row in output_rows),
            "observed": sum(
                row["regulatory_ancestry_status"] == "cross_ancestry_eqtl_limited"
                for row in output_rows
            ),
            "status": "pass",
            "note": "European Broadaway eQTL support is not ancestry-matched regulatory validation for non-EUR GWAS",
        },
        {
            "audit": "source_registry_sample_size_discrepancies",
            "expected": "report_not_repair",
            "observed": sum(row["sample_size_audit"] != "exact" for row in output_rows),
            "status": "warning",
            "note": "Source/registry differences are retained rather than repaired, including rounded 2023 cohort sizes and exact BBJ source sizes",
        },
    ]
    audit_path = candidate / "phenotype_registry_audit.tsv"
    atomic_write_tsv(
        audit_path,
        audit_rows,
        ["audit", "expected", "observed", "status", "note"],
    )
    write_stage_seal(
        candidate,
        "02_build_phenotype_registry",
        [registry_path, provenance_path, audit_path],
        [stage1, args.catalog.resolve()],
    )
    print(f"PASS: phenotype registry contains {len(output_rows)} primary GWAS")


if __name__ == "__main__":
    main()
