#!/usr/bin/env python3
"""Import accepted spatial evidence and fail closed on unavailable class covariates."""

from __future__ import annotations

import csv
from pathlib import Path

from bridge_common import CANDIDATE_ROOT, PROJECT_ROOT, require_validated_seal, sha256_file, write_tsv


SOURCE = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/spatial_context_semantic_v2_2026-08-08/final_integration/integrated_program_effects.tsv"
PROGRAMS = {"hotspot_hepatocytes_f05c535ae5bbc0b9", "hotspot_hepatocytes_48f39dd4d817a10e"}


def main() -> None:
    require_validated_seal()
    with SOURCE.open(newline="", encoding="utf-8") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["program_uid"] in PROGRAMS and row["evidence_role"] == "spatial_organization"]
    if len(rows) != 4:
        raise RuntimeError(f"Expected four accepted spatial-organization rows for M8/M20, found {len(rows)}")
    output = [{
        "analysis_id": "accepted_program_spatial_support",
        "dataset_id": row["dataset"], "program_uid": row["program_uid"], "program_label": row["program_label"],
        "biological_unit": row["biological_unit"], "n_biological": row["n_biological"],
        "effect_unit": row["effect_unit"], "estimate": row["estimate"], "p": row["pvalue"], "q": row["padj"],
        "robustness_pass": row["robustness_pass"], "evidence_state": row["evidence_state"],
        "source_dependence": row["source_dependence"], "claim_scope": row["claim_scope"],
        "source_sha256": sha256_file(SOURCE), "status": "accepted_frozen_evidence",
    } for row in rows]
    output.append({
        "analysis_id": "evidence_class_spatial_comparison", "dataset_id": "accepted_gene_level_spatial_tables",
        "program_uid": "", "program_label": "", "biological_unit": "gene", "n_biological": "",
        "effect_unit": "disease_state_only_vs_corrected_genetic_only", "estimate": "", "p": "", "q": "",
        "robustness_pass": "false", "evidence_state": "skipped", "source_dependence": "not_applicable",
        "claim_scope": "none", "source_sha256": sha256_file(SOURCE),
        "status": "skipped_missing_prespecified_gene_level_spatial_variability_and_coverage_covariates",
    })
    write_tsv(CANDIDATE_ROOT / "spatial_lineage_support.tsv", output, list(output[0]))
    print("SPATIAL_SUPPORT_COMPLETE")


if __name__ == "__main__": main()
