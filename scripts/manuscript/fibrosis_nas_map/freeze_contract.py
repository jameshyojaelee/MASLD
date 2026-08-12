#!/usr/bin/env python3
"""Freeze the outcome-independent fibrosis-NAS analysis contract."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path


ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
RELEASE_ID = "resource-f-five-coloc-v6-candidate-2026-08-10"
WORKSTREAM_ID = "BULK-PROGRAM-MAP-v8"
CANDIDATE = ROOT / "RNA-seq/results/manuscript_release/candidates" / RELEASE_ID
WORKSTREAM = CANDIDATE / "workstreams" / WORKSTREAM_ID
SCRIPT_ROOT = ROOT / "scripts/manuscript/fibrosis_nas_map"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


INPUTS = {
    "precontract_nonholdout_dge": CANDIDATE / "workstreams/BULK-PROGRAM-MAP-v2/sealed_inputs/nonholdout_dge.rds",
    "precontract_nonholdout_meta": CANDIDATE / "workstreams/BULK-PROGRAM-MAP-v2/sealed_inputs/nonholdout_meta.rds",
    "precontract_testability_audit": CANDIDATE / "workstreams/BULK-PROGRAM-MAP-v4/precontract/testability_by_program_cohort.tsv.gz",
    "corrected_nine_cohort_dge": CANDIDATE / "inputs/BG001-DECISION/arms/F_legacy/results/integration/merged_dge.rds",
    "corrected_nine_cohort_meta": CANDIDATE / "inputs/BG001-DECISION/arms/F_legacy/results/integration/meta_matched.rds",
    "gencode_v49_gene_metadata": CANDIDATE / "inputs/BULK-F-FIVE/frozen_model_inputs/gencode_v49_gene_metadata.tsv.gz",
    "gse193066_crosswalk": CANDIDATE / "workstreams/BULK-STAGE/PREFLIGHT-v1/audits/gse193066_participant_crosswalk.tsv",
    "program_registry": ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_registry_v2.tsv",
    "program_membership": ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/program_membership_v2.tsv",
    "music_proportions": ROOT / "RNA-seq/results/celltype_attribution/persample_celltype_proportions.csv",
    "composition_testability": ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/composition_sample_qc_testability.tsv",
    "composition_ready": ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/hotspot/COMPOSITION_QC_READY",
    "hallmark_gmt": ROOT / "Analysis/downstream_analysis/pathway_analysis/data/genesets/hallmark.gmt",
    "govaere_panel": ROOT / "data/published_gene_panels/govaere_2020_panel.tsv",
    "pantano_panel": ROOT / "data/published_gene_panels/pantano_2021_panel.tsv",
    "moylan_panel": ROOT / "data/published_gene_panels/moylan_2014_panel.tsv",
    "arendt_panel": ROOT / "data/published_gene_panels/arendt_2015_panel.tsv",
}


def main() -> None:
    missing = [str(path) for path in INPUTS.values() if not path.is_file()]
    if missing:
        raise SystemExit("Missing contract inputs:\n" + "\n".join(missing))

    scripts = sorted(
        path for path in SCRIPT_ROOT.iterdir()
        if path.is_file() and path.suffix in {".R", ".py", ".sbatch"}
    )
    contract = {
        "schema_version": "1.0.0",
        "release_id": RELEASE_ID,
        "workstream_id": WORKSTREAM_ID,
        "state": "outcome_independent_contract_frozen",
        "supersedes_prior_workstreams": [
            "BULK-PROGRAM-MAP-v1 used the superseded 27638-gene/1260-sample locked object and failed before partitioning",
            "BULK-PROGRAM-MAP-v2 stopped before phenotype modeling when a full-denominator coverage audit rederived 113 rather than the provisional 114 testable programs",
            "BULK-PROGRAM-MAP-v3 was aborted before contract freeze because its outcome-blind coverage audit loaded the monolithic object before subsetting, violating the literal holdout-access boundary; no holdout statistic was printed or used",
            "BULK-PROGRAM-MAP-v4 froze successfully but failed before partition publication or phenotype modeling because a manifest row duplicated its grouping column",
            "BULK-PROGRAM-MAP-v5 was an outcome-blind smoke run: all models completed, but discovery stopped before publication when an unsupported data.table pronoun prevented LOCO summary assembly; no coefficient or p-value table was inspected",
            "BULK-PROGRAM-MAP-v6 sealed discovery without opening the holdout, but a pre-holdout code audit found that whole-sample composition missingness was not removed before transport permutations; no v6 coefficient or p-value table was inspected",
            "BULK-PROGRAM-MAP-v7 sealed discovery without opening the holdout; its primary and LOCO outputs passed independent numeric validation, but its descriptive raw-mean linearity rule called every program nonmonotonic, so the audit was replaced by a complete-family HC3 adjacent-contrast test before holdout access",
        ],
        "claim": "conditional cross-sectional fibrosis- and NAS-associated molecular remodeling",
        "forbidden_methods": [
            "transcriptomic ordering", "pseudotime", "sliding windows",
            "outcome-selected genes", "classifier tuning", "hero swapping",
        ],
        "cohort_roles": {
            "discovery_joint_fibrosis_nas": {
                "GSE130970": 76, "GSE135251": 214,
                "GSE162694": 86, "GSE174478": 93,
            },
            "locked_paired_holdout": {"GSE193066": {"pairs": 54, "changed": 46, "unchanged": 8, "discordant": 11}},
            "fibrosis_only_transport": {"GSE240729": 64},
            "diagnosis_transport": {"GSE126848": 53, "GSE167523": 96, "GSE213621": 361},
        },
        "features": {
            "primary": "117 frozen Hotspot programs",
            "program_score": "L1-weighted mean of within-cohort gene z scores",
            "standardization_population": "all QC expression samples within each assay cohort; phenotype eligibility is applied only at model fitting",
            "program_sensitivities": ["unweighted mean-z", "weighted within-sample rank"],
            "program_weight_coverage_minimum": 0.80,
            "program_coverage_denominator": "sum of original_l1_weight over every frozen source membership row, including unmapped members",
            "outcome_blind_testable_programs": 113,
            "outcome_blind_untestable_programs": [
                "hotspot_cholangiocytes_fe9b79d294a58bf5",
                "hotspot_hepatocytes_c749bbf208c0ff37",
                "hotspot_hepatocytes_00250b1e5f818798",
                "hotspot_tcells_eaa9c262b58d9323",
            ],
            "comparators": ["50 Hallmark sets", "four verbatim MASLD signatures", "rebuilt cell composition"],
            "gene_set_minimum_genes": 10,
            "gene_set_coverage_minimum": 0.80,
            "composition_pseudocount": 1e-6,
        },
        "precontract_observability_audit": {
            "source": "BULK-PROGRAM-MAP-v2 sealed non-holdout partition only",
            "holdout_expression_accessed": False,
            "result": "113 programs pass 80% nonconstant frozen weight in each of four discovery cohorts",
        },
        "models": {
            "primary": "score ~ sex_final + centered_fibrosis + centered_NAS",
            "cohort_standard_errors": "HC3",
            "meta_analysis": "REML random effects with Knapp-Hartung inference",
            "primary_multiple_testing": "BH across every testable program by both axes",
            "interaction": "centered_fibrosis * centered_NAS; separate BH family",
            "linearity_audit": "categorical adjusted marginal means with HC3 adjacent-level contrasts; BH across all program-axis-cohort contrasts; reversal requires a BH-supported contrast opposing the conditional linear meta-effect",
            "age_sensitivity_cohorts": ["GSE130970", "GSE162694", "GSE174478"],
            "composition_conditioned": "score ~ sex + centered_fibrosis + centered_NAS + four non-hepatocyte-to-hepatocyte log-ratios",
        },
        "validation": {
            "loco": "four held-out cohort effect-vector Spearman correlations per axis",
            "loco_global": "mean Fisher-z across folds; 10000 joint fibrosis/NAS label permutations; Holm across two axes",
            "paired_primary": "median participant cosine between predicted and observed program deltas",
            "paired_resampling": "10000 participant bootstraps and 10000 paired-delta permutations",
            "empirical_p": "(1 + number at least as extreme) / (B + 1)",
        },
        "promotion": {
            "both_loco_axes_holm_p_below": 0.05,
            "minimum_positive_folds_per_axis": 3,
            "paired_empirical_p_below": 0.05,
            "paired_bootstrap_ci_lower_above": 0.0,
            "failed_paired_claim": "cross-sectional only; no longitudinal-change claim",
        },
        "seed": 20260811,
        "resamples": 10000,
        "input_manifest": {
            name: {"path": str(path.relative_to(ROOT)), "sha256": sha256(path)}
            for name, path in INPUTS.items()
        },
        "code_manifest": {
            str(path.relative_to(ROOT)): sha256(path) for path in scripts
        },
        "composition_rebuild": {
            "source": str(INPUTS["music_proportions"].relative_to(ROOT)),
            "census": "filter deposited MuSiC estimates to the corrected 1257-sample fragment-count object; preserve assay-missing samples as missing",
            "lineages": "all 22 deposited columns retained; 16 prespecified non-structural lineages tested",
            "transform": "replace zeros by 1e-6, renormalize, then centered log-ratio",
            "available_corrected_samples": 1219,
            "assay_missing_corrected_samples": 38,
        },
    }

    destination = WORKSTREAM / "contract" / "analysis_contract.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    encoded = json.dumps(contract, indent=2, sort_keys=True) + "\n"
    if destination.exists():
        if destination.read_text() != encoded:
            raise SystemExit(f"Frozen contract differs; refusing overwrite: {destination}")
        print(f"Contract already frozen and identical: {destination}")
        return
    with destination.open("x") as handle:
        handle.write(encoded)
    os.chmod(destination, 0o440)
    print(f"Frozen contract: {destination}")
    print(f"SHA256: {sha256(destination)}")


if __name__ == "__main__":
    main()
