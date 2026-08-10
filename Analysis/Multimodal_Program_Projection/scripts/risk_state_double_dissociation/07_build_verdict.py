#!/usr/bin/env python3
"""Integrate assay-native Plan 42 arm verdicts without inventing a universal score."""

from __future__ import annotations

import json
import math
import platform
import subprocess

from risk_state_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    SCRIPT_ROOT,
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


def main() -> None:
    seal = require_validated_seal()
    hypotheses = read_tsv(CANDIDATE_ROOT / "frozen_spec/frozen_hypotheses.tsv")
    write_tsv(
        CANDIDATE_ROOT / "frozen_hypotheses.tsv",
        hypotheses,
        list(hypotheses[0]),
    )
    required = [
        "crispr_class_effects.tsv",
        "crispr_verdict.tsv",
        "stress_response_effects.tsv",
        "stress_response_verdict.tsv",
        "source_gate_status.tsv",
    ]
    for filename in required:
        if not (CANDIDATE_ROOT / filename).is_file():
            raise FileNotFoundError(CANDIDATE_ROOT / filename)

    crispr = read_tsv(CANDIDATE_ROOT / "crispr_class_effects.tsv")
    crispr_map = {row["comparison"]: row for row in crispr}
    t_crispr = number(crispr_map["genetic_vs_state"]["estimate"])
    t_crispr_q = number(crispr_map["genetic_vs_state"]["q"])
    direct_verdict = read_tsv(CANDIDATE_ROOT / "crispr_verdict.tsv")[0]
    stress_verdict = read_tsv(CANDIDATE_ROOT / "stress_response_verdict.tsv")[0]
    gates = read_tsv(CANDIDATE_ROOT / "source_gate_status.tsv")
    gate_map = {(row["source_id"], row["gate"]): row for row in gates}
    maldi = gate_map[("HMSMA_MALDI", "public_registered_product")]
    clinical = gate_map[("HMSMA_VISIUM", "GSA_to_clinical_sample_key")]

    spatial_rows = [
        {
            "dataset_id": "HMSMA",
            "assay": "paired_Visium_MALDI_MSI",
            "biological_unit": "registered_adjacent_section",
            "contrast_or_exposure": "source_metabolic_module_5",
            "effect_unit": "standardized_donor_or_section_spatial_slope",
            "estimate": "",
            "se": "",
            "p": "",
            "q": "",
            "n_biological": 0,
            "n_technical": 0,
            "status": "skipped_no_public_registered_maldi",
            "source_dependence": "OMIX009098_controlled_access; portal_requires_author_contact",
            "detail": maldi["detail"],
            "release_id": seal["candidate_id"],
        },
        {
            "dataset_id": "HMSMA",
            "assay": "Visium_histology",
            "biological_unit": "section",
            "contrast_or_exposure": "fibrosis_conditional_on_steatosis",
            "effect_unit": "donor_level_effect",
            "estimate": "",
            "se": "",
            "p": "",
            "q": "",
            "n_biological": 0,
            "n_technical": 0,
            "status": "skipped_no_authoritative_sample_key",
            "source_dependence": "GSA_HRA_identifiers_not_joinable_to_clinical_identifiers",
            "detail": clinical["detail"],
            "release_id": seal["candidate_id"],
        },
    ]
    write_tsv(
        CANDIDATE_ROOT / "spatial_niche_effects.tsv",
        spatial_rows,
        [
            "dataset_id", "assay", "biological_unit", "contrast_or_exposure", "effect_unit",
            "estimate", "se", "p", "q", "n_biological", "n_technical", "status",
            "source_dependence", "detail", "release_id",
        ],
    )

    components = [
        {
            "statistic_id": "T_CRISPR",
            "definition": "median matched residual abs(CRISPR Z) in genetic-only minus disease-state-only",
            "estimate": t_crispr,
            "component_q": t_crispr_q,
            "expected_direction": "positive",
            "status": "valid_nonconfirmatory",
            "gate_pass": "false",
            "not_computed_reason": "",
        },
        {
            "statistic_id": "T_SPATIAL",
            "definition": "median matched residual spatial-niche Z in disease-state-only minus genetic-only",
            "estimate": "",
            "component_q": "",
            "expected_direction": "positive",
            "status": "not_computed_source_gate_failed",
            "gate_pass": "false",
            "not_computed_reason": "no public no-permission registered MALDI product and no authoritative GSA-to-clinical sample key",
        },
        {
            "statistic_id": "D",
            "definition": "rank_normalize(T_CRISPR) plus rank_normalize(T_SPATIAL)",
            "estimate": "",
            "component_q": "",
            "expected_direction": "positive; permutation p below 0.01",
            "status": "not_computed_prerequisite_failed",
            "gate_pass": "false",
            "not_computed_reason": "T_CRISPR failed corrected significance and T_SPATIAL is source-gated",
        },
    ]
    for row in components:
        row["release_id"] = seal["candidate_id"]
        row["minimum_permutations_if_computable"] = 99999
    write_tsv(
        CANDIDATE_ROOT / "double_dissociation.tsv",
        components,
        [
            "statistic_id", "definition", "estimate", "component_q", "expected_direction",
            "status", "gate_pass", "not_computed_reason", "minimum_permutations_if_computable",
            "release_id",
        ],
    )

    promotion_gates = [
        {"gate_order": 1, "gate_id": "CRISPR_COMPONENT", "passed": "false", "detail": f"T_CRISPR={t_crispr:.6g}; q={t_crispr_q:.6g}; direct_arm={direct_verdict['status']}"},
        {"gate_order": 2, "gate_id": "SPATIAL_COMPONENT", "passed": "false", "detail": f"{maldi['status']}; {clinical['status']}"},
        {"gate_order": 3, "gate_id": "INTERACTION", "passed": "false", "detail": "not computed because both prerequisite components did not pass"},
        {"gate_order": 4, "gate_id": "OBSERVABILITY", "passed": "false", "detail": "CLCC1 adjustment completed, but the required effect was nonconfirmatory; spatial adjustment not testable"},
        {"gate_order": 5, "gate_id": "PHENOTYPE_PROVENANCE", "passed": "false", "detail": "proxy-excluded CLCC1 sensitivity retained only eight matched loci and was non-significant"},
        {"gate_order": 6, "gate_id": "PROVENANCE", "passed": "true", "detail": "no source-owned headline gene or inferred sample/registration key was used"},
    ]
    for row in promotion_gates:
        row["release_id"] = seal["candidate_id"]
    write_tsv(
        CANDIDATE_ROOT / "promotion_gate_status.tsv",
        promotion_gates,
        ["gate_order", "gate_id", "passed", "detail", "release_id"],
    )

    verdict = {
        "candidate_id": seal["candidate_id"],
        "status": "complete_no_promotion",
        "functional_partition_claim_pass": False,
        "stress_contingent_language_pass": stress_verdict["promotion_gate_pass"].lower() == "true",
        "figure5_promotion": False,
        "paper_action": "retain_five_figure_Cell_Genomics_Resource_posture",
        "main_text_use": "none_as_new_mechanistic_or_functional_architecture_claim",
        "supplement_use": "complete_valid_and_null_CLCC1_PNPLA3_and_source_gate_package",
        "replacement_search_allowed": False,
        "reason": "The direct-lipid component was nonconfirmatory, the PNPLA3 stress-alignment boundary was nonconfirmatory, and the HMSMA registered MALDI/histology branch is unavailable without permission or an authoritative sample key.",
    }
    write_tsv(
        CANDIDATE_ROOT / "promotion_verdict.tsv",
        [verdict],
        list(verdict),
    )

    code_rows = []
    for path in sorted(SCRIPT_ROOT.glob("*")):
        if path.is_file() and path.suffix in {".py", ".R", ".sbatch", ".md"}:
            code_rows.append(
                {
                    "relative_path": str(path.relative_to(PROJECT_ROOT)),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                }
            )
    write_tsv(CANDIDATE_ROOT / "execution_manifest.tsv", code_rows, ["relative_path", "bytes", "sha256"])

    r_version = subprocess.run(
        ["micromamba", "run", "-n", "rnaseq", "Rscript", "--version"],
        cwd=PROJECT_ROOT, check=True, capture_output=True, text=True,
    )
    environment = [
        {"component": "python", "version": platform.python_version(), "detail": platform.platform()},
        {"component": "R_rnaseq", "version": (r_version.stderr or r_version.stdout).strip(), "detail": "micromamba environment rnaseq"},
        {"component": "openpyxl", "version": __import__("openpyxl").__version__, "detail": "CLCC1 workbook schema gate"},
    ]
    write_tsv(CANDIDATE_ROOT / "environment_manifest.tsv", environment, ["component", "version", "detail"])
    atomic_write_text(
        CANDIDATE_ROOT / "INTEGRATED.json",
        json.dumps(verdict, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(verdict, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
