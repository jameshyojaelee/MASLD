#!/usr/bin/env python3
"""Independent structural validator for the Myojin outcome-blind firewall."""

from __future__ import annotations

import argparse
import collections
import csv
import math
from pathlib import Path

from myojin_firewall_common import (
    CANDIDATE_ROOT,
    HLF_MODEL_ID,
    PROJECT_ROOT,
    RELEASE_ID,
    atomic_write_text,
    md5_file,
    parse_bool,
    read_tsv,
    require_within,
    sha256_file,
    stable_bundle_sha256,
    write_tsv,
)


EXPECTED_DEPMAP_FILES = {
    "Model.csv": (645696, "675210d17675f3517b0ce39a3c274f16"),
    "CRISPRGeneEffect.csv": (428678699, "6edf7ade09b9b34199210b559d4745d3"),
    "OmicsExpressionProteinCodingGenesTPMLogp1.csv": (
        506628654,
        "71794802b750ce77c422dad0720a40af",
    ),
}
EXPECTED_DEPMAP_ENTREZ_CONFLICTS = {
    "MEF2B": ("4207", "100271849"),
    "RLN2": ("6013", "6019"),
}
EXPECTED_CLASS_COVERAGE = {
    "neither": (9896, 9635, 7056),
    "genetic_only": (307, 297, 228),
    "disease_state_only": (956, 939, 463),
    "convergent": (28, 25, 15),
}


MASKED_ALLOWLIST = {
    "source_row",
    "source_gene_id",
    "guide_count_d21_d0",
    "guide_count_pa_vehicle",
    "d21_d0_score_present",
    "d21_d0_pvalue_present",
    "d21_d0_fdr_present",
    "d21_d0_lfc_present",
    "pa_vehicle_score_present",
    "pa_vehicle_pvalue_present",
    "pa_vehicle_fdr_present",
    "pa_vehicle_lfc_present",
    "pa_vehicle_primary_fields_complete",
    "source_gene_id_duplicate",
    "outcome_values_exported",
}
FORBIDDEN_VALUE_COLUMNS = {
    "pos|score",
    "pos|p-value",
    "pos|fdr",
    "pos|rank",
    "pos|goodsgrna",
    "pos|lfc",
    "pa_vs_vehicle_pos_lfc",
    "pa_vs_vehicle_pos_p",
    "pa_vs_vehicle_pos_fdr",
    "protective_hit",
}


def header(path: Path) -> list[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def add_check(checks: list[dict[str, str]], name: str, passed: bool, detail: str) -> None:
    checks.append(
        {
            "check": name,
            "passed": str(passed).upper(),
            "detail": detail,
            "release_id": RELEASE_ID,
        }
    )
    if not passed:
        raise AssertionError(f"{name}: {detail}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=CANDIDATE_ROOT)
    args = parser.parse_args()
    outdir = require_within(args.outdir, CANDIDATE_ROOT)
    required = [
        "PHASE_A_READY",
        "source_manifest.tsv",
        "source_schema_audit.tsv",
        "masked_screen_schema.tsv",
        "blind_analysis_spec.yaml",
        "specification_sha256.txt",
        "screen_gene_universe.tsv",
        "gene_mapping_audit.tsv",
        "covariate_coverage_audit.tsv",
        "prediction_manifest.tsv",
        "program_testability.tsv",
        "known_hit_exclusion.tsv",
        "gate_status.tsv",
        "firewall_input_manifest.tsv",
        "release_manifest.tsv",
        "local_depmap_preflight.tsv",
        "deviations.tsv",
        "unseal_record.tsv",
    ]
    if (outdir / "DEPMAP_READY").is_file():
        required.extend(
            [
                "depmap_source_manifest.tsv",
                "depmap_gene_mapping_audit.tsv",
                "hlf_identity_audit.tsv",
                "detectability_covariates.tsv",
            ]
        )
    missing = [name for name in required if not (outdir / name).is_file()]
    if missing:
        raise FileNotFoundError(f"Missing blind-phase artifacts: {missing}")
    if not ((outdir / "SEALED").is_file() ^ (outdir / "BLOCKED").is_file()):
        raise AssertionError("Exactly one of SEALED or BLOCKED must exist")

    checks: list[dict[str, str]] = []
    masked_path = outdir / "masked_screen_schema.tsv"
    masked_header = set(header(masked_path))
    add_check(
        checks,
        "masked_schema_exact_allowlist",
        masked_header == MASKED_ALLOWLIST,
        f"observed={sorted(masked_header)}",
    )
    add_check(
        checks,
        "masked_schema_no_outcome_value_columns",
        not (masked_header & FORBIDDEN_VALUE_COLUMNS),
        f"intersection={sorted(masked_header & FORBIDDEN_VALUE_COLUMNS)}",
    )
    masked = read_tsv(masked_path)
    add_check(checks, "source_record_count", len(masked) == 18343, f"n={len(masked)}")
    add_check(
        checks,
        "outcome_export_marker_false",
        all(row["outcome_values_exported"] == "FALSE" for row in masked),
        "all masked rows must be FALSE",
    )
    add_check(
        checks,
        "masked_guide_counts_nonnegative",
        all(
            row[field] == "" or int(row[field]) >= 0
            for row in masked
            for field in ("guide_count_d21_d0", "guide_count_pa_vehicle")
        ),
        "both neutral guide-count fields audited",
    )

    source_manifest = read_tsv(outdir / "source_manifest.tsv")
    add_check(
        checks,
        "raw_workbook_checksum",
        len(source_manifest) == 1
        and source_manifest[0]["sha256"]
        == "2663a811676a1abe2778375e6d28f6397eaded1f895d274f3c241733669c1c66",
        source_manifest[0]["sha256"] if source_manifest else "missing",
    )
    raw_workbook = outdir / source_manifest[0]["relative_path"]
    add_check(
        checks,
        "raw_workbook_file_authenticated",
        raw_workbook.is_file()
        and raw_workbook.stat().st_size == int(source_manifest[0]["bytes"])
        and sha256_file(raw_workbook) == source_manifest[0]["sha256"]
        and raw_workbook.stat().st_mode & 0o222 == 0,
        "on-disk bytes/hash match manifest and source is read-only",
    )
    add_check(
        checks,
        "raw_workbook_shape",
        source_manifest[0]["sheet_count"] == "1"
        and source_manifest[0]["used_range"] == "A1:Q18346"
        and source_manifest[0]["gene_records"] == "18343",
        "one sheet; A1:Q18346; 18,343 records",
    )

    prediction = read_tsv(outdir / "prediction_manifest.tsv")
    add_check(
        checks,
        "prediction_outcomes_unread",
        prediction
        and all(row["external_outcomes_read"] == "FALSE" for row in prediction),
        f"rows={len(prediction)}",
    )
    class_prediction = [row for row in prediction if row["object_type"] == "evidence_class_gene"]
    observed_classes = {row["primary_evidence_class"] for row in class_prediction}
    add_check(
        checks,
        "four_classes_mutually_labeled",
        observed_classes == {"neither", "genetic_only", "disease_state_only", "convergent"},
        f"classes={sorted(observed_classes)}",
    )
    symbol_labels: dict[str, set[str]] = collections.defaultdict(set)
    for row in class_prediction:
        symbol_labels[row["gene_symbol"]].add(row["primary_evidence_class"])
    add_check(
        checks,
        "class_membership_mutually_exclusive",
        all(len(labels) == 1 for labels in symbol_labels.values()),
        f"symbols={len(symbol_labels)}",
    )
    add_check(
        checks,
        "genetics_only_direction_unspecified",
        all(
            row["expected_direction"] == "not_specified"
            for row in class_prediction
            if row["primary_evidence_class"] == "genetic_only"
        ),
        "no allele-direction invention",
    )

    program_prediction = [
        row for row in prediction if row["object_type"] == "hepatocyte_program_gene"
    ]
    program_weights: dict[str, float] = collections.defaultdict(float)
    for row in program_prediction:
        weight = float(row["frozen_weight"])
        if not math.isfinite(weight) or weight <= 0:
            raise AssertionError("Program weights must be finite and positive")
        program_weights[row["object_id"]] += weight
    add_check(
        checks,
        "program_weight_sums",
        len(program_weights) == 2
        and all(math.isclose(value, 1.0, rel_tol=0, abs_tol=1e-12) for value in program_weights.values()),
        f"sums={dict(program_weights)}",
    )
    add_check(
        checks,
        "program_registry_hash",
        all(
            row["registry_sha256"]
            == "b134b5f8e46a716e857af22271aa54e1294432b9f9ab735559d7f603cdd29fe7"
            for row in program_prediction
        ),
        "HS-V2 registry seal",
    )

    known = read_tsv(outdir / "known_hit_exclusion.tsv")
    add_check(
        checks,
        "known_hit_set",
        {row["canonical_symbol"] for row in known}
        == {"ACSL3", "INSIG1", "NF2", "CASP8", "GPAT3", "RNF213"},
        f"symbols={sorted(row['canonical_symbol'] for row in known)}",
    )
    add_check(
        checks,
        "known_hit_primary_and_sensitivity",
        all(
            row["included_in_primary"] == "TRUE"
            and row["remove_all_known_hits_sensitivity"] == "TRUE"
            for row in known
        ),
        "primary remains genome-wide; exclusion is mandatory sensitivity",
    )

    universe = read_tsv(outdir / "screen_gene_universe.tsv")
    add_check(
        checks,
        "universe_outcomes_unread",
        len(universe) == 18343 and all(row["outcome_values_read"] == "FALSE" for row in universe),
        f"rows={len(universe)}",
    )
    add_check(
        checks,
        "actual_vehicle_fitness_not_observed",
        all(row["actual_vehicle_fitness_observed"] == "FALSE" for row in universe),
        "Chronos is generic HLF fitness only",
    )

    local_preflight = read_tsv(outdir / "local_depmap_preflight.tsv")
    local_preflight_pass = len(local_preflight) == 3 and all(
        row["accepted_for_confirmatory_use"] == "FALSE" for row in local_preflight
    )
    for row in local_preflight:
        if row["path"] == "missing":
            local_preflight_pass &= row["sha256"] == "NA"
            continue
        source = PROJECT_ROOT / row["path"]
        local_preflight_pass &= source.is_file() and sha256_file(source) == row["sha256"]
    add_check(
        checks,
        "unversioned_local_depmap_rejected",
        local_preflight_pass,
        "local Model/Chronos hashes verified; missing matched expression; all rejected",
    )
    depmap_ready_path = outdir / "DEPMAP_READY"
    if depmap_ready_path.is_file():
        depmap_manifest = read_tsv(outdir / "depmap_source_manifest.tsv")
        manifest_by_name = {row["filename"]: row for row in depmap_manifest}
        source_root = outdir / "source/depmap_24q4"
        depmap_sources_pass = set(manifest_by_name) == set(EXPECTED_DEPMAP_FILES)
        for filename, (expected_bytes, expected_md5) in EXPECTED_DEPMAP_FILES.items():
            row = manifest_by_name.get(filename, {})
            source = source_root / filename
            depmap_sources_pass &= (
                source.is_file()
                and source.stat().st_size == expected_bytes
                and md5_file(source) == expected_md5
                and row.get("depmap_release") == "DepMap Public 24Q4"
                and row.get("doi") == "10.25452/figshare.plus.27993248.v1"
                and row.get("bytes") == str(expected_bytes)
                and row.get("published_md5") == expected_md5
                and row.get("observed_md5") == expected_md5
                and row.get("observed_sha256") == sha256_file(source)
                and row.get("license") == "CC_BY_4.0"
            )
        add_check(
            checks,
            "official_depmap_24q4_triplet_authenticated",
            depmap_sources_pass,
            "three official Figshare files match bytes, MD5, SHA256, DOI, release, and license",
        )

        identity = read_tsv(outdir / "hlf_identity_audit.tsv")
        identity_pass = (
            len(identity) == 1
            and identity[0]["model_id"] == HLF_MODEL_ID
            and identity[0]["cell_line_name"] == "HLF"
            and identity[0]["stripped_cell_line_name"] == "HLF"
            and identity[0]["oncotree_lineage"] == "Liver"
            and identity[0]["hlf_unique"] == "TRUE"
            and identity[0]["hlfa_conflated"] == "FALSE"
            and identity[0]["actual_vehicle_fitness_observed"] == "FALSE"
        )
        add_check(
            checks,
            "depmap_hlf_identity",
            identity_pass,
            "ACH-000393 is unique HLF liver/HCC and is not HLF-a",
        )

        covariates = read_tsv(outdir / "detectability_covariates.tsv")
        covariates_by_symbol = {row["gene_symbol"]: row for row in covariates}
        mapping_audit = read_tsv(outdir / "depmap_gene_mapping_audit.tsv")
        conflicts = {
            row["gene_symbol"]: row
            for row in mapping_audit
            if row["depmap_entrez_conflict"] == "TRUE"
        }
        observed_conflicts = {
            symbol: (row["expression_entrez_id"], row["chronos_entrez_id"])
            for symbol, row in conflicts.items()
        }
        conflict_pass = observed_conflicts == EXPECTED_DEPMAP_ENTREZ_CONFLICTS
        for symbol, row in conflicts.items():
            covariate = covariates_by_symbol.get(symbol, {})
            conflict_pass &= (
                row["depmap_mapping_state"] == "entrez_conflict_values_blanked"
                and row["covariate_values_blanked"] == "TRUE"
                and row["eligible_for_primary_mapping"] == "FALSE"
                and covariate.get("depmap_entrez_conflict") == "TRUE"
                and covariate.get("HLF_TPM") == ""
                and covariate.get("log1p_HLF_TPM") == ""
                and covariate.get("HLF_Chronos") == ""
            )
        conflict_universe = [
            row for row in universe if row["gene_symbol"] in observed_conflicts
        ]
        conflict_pass &= bool(conflict_universe) and all(
            row["depmap_entrez_conflict"] == "TRUE"
            and row["primary_eligible"] == "FALSE"
            and "depmap_entrez_conflict_values_blanked" in row["exclusion_reasons"]
            for row in conflict_universe
        )
        add_check(
            checks,
            "depmap_entrez_conflicts_fail_closed",
            conflict_pass,
            f"conflicts={observed_conflicts}; no symbol-only merge; covariates blanked",
        )

        complete_covariates = sum(
            row["HLF_TPM"] != "" and row["HLF_Chronos"] != ""
            for row in covariates
        )
        add_check(
            checks,
            "hlf_expression_and_chronos_covariates",
            len(covariates) == len(mapping_audit)
            and len(covariates_by_symbol) == len(covariates)
            and complete_covariates > 0
            and all(row["model_id"] == HLF_MODEL_ID for row in covariates),
            f"genes={len(covariates)} complete_both={complete_covariates}",
        )

        ready = read_tsv(depmap_ready_path)
        depmap_bundle_hash, _ = stable_bundle_sha256(
            [
                outdir / "detectability_covariates.tsv",
                outdir / "depmap_gene_mapping_audit.tsv",
                outdir / "hlf_identity_audit.tsv",
                outdir / "depmap_source_manifest.tsv",
            ],
            outdir,
        )
        add_check(
            checks,
            "depmap_ready_bundle_rederived",
            len(ready) == 1
            and ready[0]["outcomes_read"] == "FALSE"
            and ready[0]["depmap_release"] == "DepMap Public 24Q4"
            and ready[0]["entrez_conflicts_excluded"] == "2"
            and ready[0]["bundle_sha256"] == depmap_bundle_hash,
            f"bundle_sha256={depmap_bundle_hash}",
        )

    coverage_audit = read_tsv(outdir / "covariate_coverage_audit.tsv")
    coverage_by_class = {
        row["primary_evidence_class"]: row for row in coverage_audit
    }
    coverage_pass = set(coverage_by_class) == set(EXPECTED_CLASS_COVERAGE)
    rederived_coverage: dict[str, tuple[int, int, int]] = {}
    for class_name in EXPECTED_CLASS_COVERAGE:
        class_rows = [
            row
            for row in universe
            if row["primary_evidence_class"] == class_name
            and row["neutral_source_base"] == "TRUE"
        ]
        complete = [
            row
            for row in class_rows
            if row["primary_covariates_numeric_complete"] == "TRUE"
        ]
        detectable = [
            row
            for row in complete
            if row["HLF_expression_floor_pass"] == "TRUE"
        ]
        observed = (len(class_rows), len(complete), len(detectable))
        rederived_coverage[class_name] = observed
        audit = coverage_by_class.get(class_name, {})
        coverage_pass &= (
            observed == EXPECTED_CLASS_COVERAGE[class_name]
            and audit.get("neutral_source_base_n") == str(observed[0])
            and audit.get("numeric_TPM_and_Chronos_complete_n") == str(observed[1])
            and audit.get("TPM_ge_1_primary_eligible_n") == str(observed[2])
            and audit.get("coverage_gate_uses_expression_floor") == "FALSE"
            and audit.get("external_outcomes_read") == "FALSE"
        )
    add_check(
        checks,
        "covariate_completeness_separate_from_expression_floor",
        coverage_pass,
        f"denominator/complete/TPM>=1={rederived_coverage}",
    )
    add_check(
        checks,
        "nonreference_numeric_coverage_gate",
        all(
            coverage_by_class[class_name]["coverage_gate_pass"] == "TRUE"
            and EXPECTED_CLASS_COVERAGE[class_name][1]
            / EXPECTED_CLASS_COVERAGE[class_name][0]
            >= 0.80
            for class_name in ("genetic_only", "disease_state_only", "convergent")
        ),
        "numeric TPM+Chronos coverage passes; TPM>=1 remains eligibility-only",
    )
    add_check(
        checks,
        "eligibility_mechanical",
        all(
            not parse_bool(row["primary_eligible"])
            or (
                row["mapping_state"]
                in {"unique_protein_coding_exact", "unique_protein_coding_alias"}
                and row["primary_evidence_class"]
                in {"neither", "genetic_only", "disease_state_only", "convergent"}
                and row["depmap_entrez_conflict"] == "FALSE"
                and row["neutral_source_base"] == "TRUE"
                and row["primary_covariates_numeric_complete"] == "TRUE"
                and row["HLF_expression_floor_pass"] == "TRUE"
                and int(row["guide_count_min"]) >= 4
                and row["direct_primary_fields_complete"] == "TRUE"
                and float(row["HLF_TPM"]) >= 1.0
                and row["HLF_Chronos"] != ""
            )
            for row in universe
        ),
        f"eligible={sum(parse_bool(row['primary_eligible']) for row in universe)}",
    )

    program_testability = read_tsv(outdir / "program_testability.tsv")
    add_check(
        checks,
        "program_testability_mechanical",
        len(program_testability) == 2
        and all(
            row["testable_if_unsealed_now"] == "FALSE"
            or (
                int(row["n_primary_eligible_genes"]) >= 8
                and float(row["retained_original_L1_weight"]) >= 0.20
            )
            for row in program_testability
        ),
        ";".join(
            f"{row['program_uid']}:{row['n_primary_eligible_genes']}/{row['retained_original_L1_weight']}"
            for row in program_testability
        ),
    )

    gates = read_tsv(outdir / "gate_status.tsv")
    gate_map = {row["gate"]: parse_bool(row["passed"]) for row in gates}
    all_gates = all(gate_map.values())
    is_sealed = (outdir / "SEALED").is_file()
    add_check(
        checks,
        "seal_matches_gates",
        is_sealed == all_gates,
        f"sealed={is_sealed} all_gates={all_gates}",
    )
    add_check(
        checks,
        "no_unseal_record",
        len(read_tsv(outdir / "unseal_record.tsv")) == 0 and not (outdir / "UNSEALED").exists(),
        "phase C has not begun",
    )
    add_check(
        checks,
        "no_deviations",
        len(read_tsv(outdir / "deviations.tsv")) == 0,
        "pre-unseal specification has no deviations",
    )

    first_line = (outdir / "specification_sha256.txt").read_text(encoding="utf-8").splitlines()[0]
    label, recorded_spec_hash = first_line.split("\t")
    add_check(
        checks,
        "specification_hash_label",
        label == "specification_bundle_sha256" and len(recorded_spec_hash) == 64,
        first_line,
    )
    spec_text = (outdir / "blind_analysis_spec.yaml").read_text(encoding="utf-8")
    exact_analysis_rules = (
        "cross_matrix_entrez_conflict_policy: blank_both_covariates_mark_mapping_ambiguous_and_untestable",
        "covariate_coverage_definition: finite_numeric_TPM_and_Chronos_before_expression_floor",
        "expression_floor_role: primary_eligibility_and_testability_only_not_missingness",
        "coefficient_confidence_interval: HC3_robust_Wald_95pct",
        "standardized_pairwise_effect: adjusted_class_coefficient_divided_by_full_model_residual_SD",
        "risk_difference_estimator: shared_strata_size_weighted_class_minus_neither_difference_in_proportions",
        "draw_rule: one_control_per_program_member_without_replacement_within_draw",
        "matched_control_weight: inherit_the_corresponding_program_member_frozen_weight",
        "insufficient_exact_cell_policy: program_untestable_no_relaxation",
    )
    add_check(
        checks,
        "analysis_estimands_fully_frozen",
        all(rule in spec_text for rule in exact_analysis_rules),
        "class scale/CI, binary estimator, and matched-null construction are explicit",
    )
    component_paths = [
        outdir / "blind_analysis_spec.yaml",
        outdir / "prediction_manifest.tsv",
        outdir / "known_hit_exclusion.tsv",
        outdir / "screen_gene_universe.tsv",
        outdir / "gene_mapping_audit.tsv",
        outdir / "covariate_coverage_audit.tsv",
        outdir / "program_testability.tsv",
        outdir / "gate_status.tsv",
        outdir / "firewall_input_manifest.tsv",
    ]
    rederived_hash, _ = stable_bundle_sha256(component_paths, outdir)
    add_check(
        checks,
        "specification_hash_rederived",
        rederived_hash == recorded_spec_hash,
        f"recorded={recorded_spec_hash} rederived={rederived_hash}",
    )

    manifest = read_tsv(outdir / "firewall_input_manifest.tsv")
    producer_rows = [row for row in manifest if row["role"] == "producer"]
    producer_hashes_match = True
    for row in producer_rows:
        path = PROJECT_ROOT / row["path"]
        producer_hashes_match &= path.is_file() and sha256_file(path) == row["sha256"]
    add_check(
        checks,
        "producer_hashes_match",
        len(producer_rows) == 5 and producer_hashes_match,
        f"producers={len(producer_rows)}",
    )

    release_manifest = read_tsv(outdir / "release_manifest.tsv")
    release_manifest_pass = bool(release_manifest)
    for row in release_manifest:
        path = outdir / row["relative_path"]
        release_manifest_pass &= (
            path.is_file()
            and sha256_file(path) == row["sha256"]
            and row["external_outcomes_read"] == "FALSE"
            and row["canonical"] == "FALSE"
        )
    add_check(
        checks,
        "release_manifest_rederived",
        release_manifest_pass,
        f"artifacts={len(release_manifest)}",
    )

    validation_path = outdir / "validation_status.tsv"
    write_tsv(validation_path, checks, ["check", "passed", "detail", "release_id"])
    validation_hash = sha256_file(validation_path)
    if is_sealed:
        atomic_write_text(
            outdir / "FIREWALL_READY",
            "release_id\tstatus\tspecification_sha256\tvalidation_sha256\texternal_outcomes_read\n"
            f"{RELEASE_ID}\tready_for_one_pass_unseal\t{recorded_spec_hash}\t{validation_hash}\tFALSE\n",
        )
        print(
            f"FIREWALL_VALIDATION_PASS checks={len(checks)} sealed=TRUE "
            f"specification_sha256={recorded_spec_hash} outcomes_read=FALSE"
        )
    else:
        ready_path = outdir / "FIREWALL_READY"
        if ready_path.exists():
            ready_path.unlink()
        print(
            f"FIREWALL_BLOCKED_VALIDATION_PASS checks={len(checks)} sealed=FALSE "
            f"specification_sha256={recorded_spec_hash} outcomes_read=FALSE"
        )


if __name__ == "__main__":
    main()
    md5_file,
