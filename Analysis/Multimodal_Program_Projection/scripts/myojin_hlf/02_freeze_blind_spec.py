#!/usr/bin/env python3
"""Freeze outcome-blind HLF predictions, universe, and analysis decisions."""

from __future__ import annotations

import argparse
import collections
import csv
import datetime as dt
import gzip
import math
import re
from pathlib import Path

from myojin_firewall_common import (
    CANDIDATE_ROOT,
    HLF_MODEL_ID,
    PROJECT_ROOT,
    RELEASE_ID,
    atomic_write_text,
    parse_bool,
    read_tsv,
    require_within,
    sha256_file,
    split_gene_header,
    stable_bundle_sha256,
    write_tsv,
)


CLASS_PATH = PROJECT_ROOT / "RNA-seq/results/manuscript_release/2026-07-15-r2/evidence_class_table.tsv"
HOTSPOT_ROOT = (
    PROJECT_ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / RELEASE_ID
    / "hotspot"
)
GENCODE_PATH = PROJECT_ROOT / "data/gencode_v49_gene_metadata.tsv.gz"
EXPECTED_HASHES = {
    "evidence_class_table.tsv": "427659ff0884ee913855b41508f3b63a3fa85f61ec0fb7caaa7fc76be626b247",
    "READY": "cf109fcb2518470600764cbfb88dde59681e59c0dbb61dc3e0f3d7dfb2a45966",
    "external_test_programs.tsv": "c84f6fc47b6d274024862e698fb571714a0c88060cf877ad836940af4d87ca17",
    "program_membership_v2.tsv": "feec3fc9ccaaffaa9abc1ac5d593cc06845ed8140460c8873a51bc8652d7336b",
    "gencode_v49_gene_metadata.tsv.gz": "49dcaae77323ffb20f5fd8399bb2899cd553e3634755a37b60a133d62549c189",
}
FOUR_CLASSES = ("neither", "genetic_only", "disease_state_only", "convergent")
KNOWN_HITS = (
    ("ACSL3", "ACSL3", "paper_text_and_figure_1"),
    ("INSIG1", "INSIG1", "paper_text_and_figure_1"),
    ("NF2", "NF2", "paper_text_and_figure_1"),
    ("CASP8", "CASP8", "paper_figure_1"),
    ("AGPAT9/GPAT3", "GPAT3", "paper_figure_1;GENCODE_v49_symbol_GPAT3"),
    ("RNF213", "RNF213", "paper_figure_1"),
)
ALIASES = {"AGPAT9": "GPAT3"}


def verify_hash(path: Path, key: str) -> None:
    observed = sha256_file(path)
    expected = EXPECTED_HASHES[key]
    if observed != expected:
        raise ValueError(f"Frozen hash mismatch for {path}: {observed} != {expected}")


def load_gencode() -> dict[str, list[dict[str, str]]]:
    by_symbol: dict[str, list[dict[str, str]]] = collections.defaultdict(list)
    with gzip.open(GENCODE_PATH, "rt", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            by_symbol[row["gene_name"]].append(row)
    return by_symbol


def load_classes() -> tuple[dict[str, dict[str, str]], list[dict[str, str]]]:
    rows = read_tsv(CLASS_PATH)
    by_symbol: dict[str, dict[str, str]] = {}
    for row in rows:
        symbol = row["symbol"]
        if not symbol:
            continue
        if symbol in by_symbol:
            raise ValueError(f"Duplicate evidence-class symbol: {symbol}")
        by_symbol[symbol] = row
    return by_symbol, rows


def load_programs() -> tuple[list[dict[str, str]], dict[str, list[dict[str, object]]]]:
    programs = read_tsv(HOTSPOT_ROOT / "external_test_programs.tsv")
    if not programs:
        raise ValueError("No external-test programs")
    if any(row["cell_type"] != "hepatocytes" for row in programs):
        raise ValueError("Non-hepatocyte program entered the HLF external-test universe")
    program_ids = {row["program_uid"] for row in programs}
    if len(program_ids) != len(programs):
        raise ValueError("Duplicate external-test program IDs")
    aggregated: dict[str, dict[str, float]] = {
        program_uid: collections.defaultdict(float) for program_uid in program_ids
    }
    for row in read_tsv(HOTSPOT_ROOT / "program_membership_v2.tsv"):
        program_uid = row["program_uid"]
        if program_uid not in program_ids:
            continue
        symbol = row["mapped_symbol"]
        weight = float(row["original_l1_weight"])
        if not symbol or not math.isfinite(weight) or weight <= 0:
            raise ValueError(f"Invalid frozen program member: {program_uid} {symbol} {weight}")
        aggregated[program_uid][symbol] += weight
    membership: dict[str, list[dict[str, object]]] = {}
    for program_uid, weights in aggregated.items():
        total = sum(weights.values())
        if not math.isclose(total, 1.0, rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"Program weights do not sum to one: {program_uid} {total}")
        membership[program_uid] = [
            {"gene_symbol": symbol, "original_l1_weight": weight}
            for symbol, weight in sorted(weights.items())
        ]
    return programs, membership


def expected_gene_direction(class_row: dict[str, str]) -> str:
    evidence_class = class_row["primary_evidence_class"]
    if evidence_class not in {"disease_state_only", "convergent"}:
        return "not_specified"
    effect = float(class_row["bulk_logFC"])
    if effect > 0:
        return "positive_KO_LFC_counter_state_compatible"
    if effect < 0:
        return "negative_KO_LFC_counter_state_compatible"
    return "zero_bulk_effect_no_direction"


def load_covariates(outdir: Path) -> tuple[dict[str, dict[str, str]], bool]:
    ready = outdir / "DEPMAP_READY"
    covariate_path = outdir / "detectability_covariates.tsv"
    identity_path = outdir / "hlf_identity_audit.tsv"
    manifest_path = outdir / "depmap_source_manifest.tsv"
    mapping_audit_path = outdir / "depmap_gene_mapping_audit.tsv"
    if not all(
        path.is_file()
        for path in (
            ready,
            covariate_path,
            identity_path,
            manifest_path,
            mapping_audit_path,
        )
    ):
        return {}, False
    ready_rows = read_tsv(ready)
    if (
        len(ready_rows) != 1
        or ready_rows[0]["outcomes_read"] != "FALSE"
        or int(ready_rows[0]["entrez_conflicts_excluded"]) < 0
    ):
        raise ValueError("Invalid DEPMAP_READY outcome-blind marker")
    identity_rows = read_tsv(identity_path)
    if (
        len(identity_rows) != 1
        or identity_rows[0]["model_id"] != HLF_MODEL_ID
        or identity_rows[0]["hlf_unique"] != "TRUE"
        or identity_rows[0]["hlfa_conflated"] != "FALSE"
        or identity_rows[0]["actual_vehicle_fitness_observed"] != "FALSE"
    ):
        raise ValueError("HLF identity gate failed")
    rows = read_tsv(covariate_path)
    covariates = {row["gene_symbol"]: row for row in rows}
    if len(covariates) != len(rows):
        raise ValueError("Duplicate DepMap covariate symbols")
    mapping_audit = read_tsv(mapping_audit_path)
    conflicts = [row for row in mapping_audit if row["depmap_entrez_conflict"] == "TRUE"]
    if len(conflicts) != int(ready_rows[0]["entrez_conflicts_excluded"]):
        raise ValueError("DEPMAP_READY conflict count does not match mapping audit")
    for conflict in conflicts:
        covariate = covariates.get(conflict["gene_symbol"])
        if (
            covariate is None
            or covariate["depmap_entrez_conflict"] != "TRUE"
            or covariate["HLF_TPM"] != ""
            or covariate["log1p_HLF_TPM"] != ""
            or covariate["HLF_Chronos"] != ""
            or conflict["covariate_values_blanked"] != "TRUE"
            or conflict["eligible_for_primary_mapping"] != "FALSE"
        ):
            raise ValueError(f"DepMap conflict was not fail-closed: {conflict['gene_symbol']}")
    return covariates, True


def map_source_gene(
    raw: str, gencode: dict[str, list[dict[str, str]]]
) -> tuple[str, str, str, str]:
    parsed, _ = split_gene_header(raw)
    canonical = ALIASES.get(parsed, parsed)
    records = gencode.get(canonical, [])
    if not records:
        return parsed, canonical, "unmapped", ""
    if len(records) != 1:
        return parsed, canonical, "ambiguous_gencode_symbol", ""
    record = records[0]
    if record["gene_biotype"] != "protein_coding":
        return parsed, canonical, "unique_non_protein_coding", record["ensembl_base"]
    state = "unique_protein_coding_alias" if canonical != parsed else "unique_protein_coding_exact"
    return parsed, canonical, state, record["ensembl_base"]


def depmap_covariate_state(covariate: dict[str, str]) -> dict[str, object]:
    """Separate numeric availability from the prespecified TPM eligibility floor."""

    tpm = (
        float(covariate["HLF_TPM"])
        if covariate.get("HLF_TPM", "") != ""
        else None
    )
    chronos = (
        float(covariate["HLF_Chronos"])
        if covariate.get("HLF_Chronos", "") != ""
        else None
    )
    tpm_available = tpm is not None and math.isfinite(tpm)
    chronos_available = chronos is not None and math.isfinite(chronos)
    entrez_conflict = parse_bool(covariate.get("depmap_entrez_conflict", "FALSE"))
    complete = tpm_available and chronos_available and not entrez_conflict
    return {
        "tpm": tpm,
        "chronos": chronos,
        "tpm_available": tpm_available,
        "chronos_available": chronos_available,
        "entrez_conflict": entrez_conflict,
        "complete": complete,
        "expression_floor_pass": tpm_available and tpm >= 1.0,
    }


def blind_spec_yaml(sealed: bool, depmap_ready: bool) -> str:
    status = "sealed_ready_for_unseal" if sealed else "specification_fixed_data_gate_blocked"
    return f"""# Outcome-blind specification; no Myojin genome-wide outcome value was read.
release_id: {RELEASE_ID}
specification_status: {status}
external_outcomes_read: false
question: "Do frozen evidence classes or hepatocyte programs differ in HLF knockout survival under 200 uM palmitate relative to vehicle?"
interpretation: phenotype_specific_compatibility_or_falsification_not_generic_functional_validation
source:
  paper: "Myojin et al., Hepatology Communications 2026, e0884"
  supplemental_table_url: "https://links.lww.com/HC9/C218"
  primary_contrast: "Palmitic acid vs Vehicle"
  direct_vehicle_contrast_used: true
  actual_vehicle_fitness_observed: false
  raw_T0_counts_publicly_recovered: false
  T0_status: T0_coverage_source_QC_only
depmap:
  required_release: "DepMap Public 24Q4"
  doi: "10.25452/figshare.plus.27993248.v1"
  model_id: {HLF_MODEL_ID}
  model_identity: HLF_hepatocellular_carcinoma_not_HLF_a
  complete_version_matched_triplet_ready: {str(depmap_ready).lower()}
  expression_file: OmicsExpressionProteinCodingGenesTPMLogp1.csv
  expression_unit: log2_TPM_plus_1
  chronos_role: generic_HLF_fitness_covariate_not_Myojin_vehicle_effect
  cross_matrix_entrez_conflict_policy: blank_both_covariates_mark_mapping_ambiguous_and_untestable
  pinned_cross_matrix_entrez_conflicts:
    MEF2B: {{expression: 4207, Chronos: 100271849}}
    RLN2: {{expression: 6013, Chronos: 6019}}
eligible_universe:
  identifier: unique_unambiguous_GENCODE_v49_symbol_and_nonconflicting_DepMap_Entrez
  biotype: protein_coding
  classes: [neither, genetic_only, disease_state_only, convergent]
  primary_result_presence: direct_PA_vs_vehicle_score_p_fdr_lfc_all_nonmissing
  guide_count_definition: min_num_across_D21_vs_D0_and_PA_vs_vehicle
  primary_guide_floor: 4
  sensitivity_guide_floor: 5
  expression_floor_TPM: 1.0
  primary_covariates_complete: [log1p_HLF_TPM, guide_count_min, HLF_Chronos]
  nonreference_class_covariate_coverage_gate: 0.80
  covariate_coverage_definition: finite_numeric_TPM_and_Chronos_before_expression_floor
  expression_floor_role: primary_eligibility_and_testability_only_not_missingness
  prohibited_eligibility_fields: [goodsgrna, MAGeCK_p, MAGeCK_FDR, MAGeCK_rank, PA_vs_D0_LFC, PA_vs_vehicle_LFC]
evidence_classes:
  source_release: 2026-07-15-r2
  reference: neither
  tested: [genetic_only, disease_state_only, convergent]
  hypothesis: omnibus_difference_not_ordered_superiority
programs:
  source: Plan20_external_test_eligible_hepatocyte_programs
  minimum_unique_genes: 8
  minimum_retained_original_L1_weight: 0.20
  primary_weight: frozen_positive_original_L1_weight_renormalized_over_eligible_genes
  sensitivities: [equal_weight, leave_highest_original_weight_out]
primary_outcomes:
  continuous: PA_vs_vehicle_positive_selection_LFC
  binary: "PA_vs_vehicle_positive_selection_FDR < 0.05 AND LFC >= 0.25"
  negative_lfc_policy: continuous_direction_only_never_source_significant_sensitizing_hit
class_model:
  formula: "PA_vs_vehicle_pos_lfc ~ evidence_class + log1p_HLF_TPM + guide_count_min + HLF_Chronos"
  reduced_formula_for_permutation: "PA_vs_vehicle_pos_lfc ~ log1p_HLF_TPM + guide_count_min + HLF_Chronos"
  primary_population: complete_case
  coefficient_estimator: OLS
  coefficient_confidence_interval: HC3_robust_Wald_95pct
  omnibus_statistic: partial_F_for_three_class_indicators
  standardized_pairwise_effect: adjusted_class_coefficient_divided_by_full_model_residual_SD
  residual_SD_definition: sqrt_full_model_RSS_divided_by_n_minus_full_model_rank
  design_rank_gate: full_and_reduced_designs_must_be_full_column_rank_else_block
  pairwise_contrasts: [genetic_only_vs_neither, disease_state_only_vs_neither, convergent_vs_neither]
  permutation: Freedman_Lane_residual
  permutation_pvalue: "(1 + count(null_statistic >= observed_statistic)) / (accepted_draws + 1)"
  accepted_draws: 100000
  seed: 2026080740
  strata:
    guide_count: [4, 5, 6, 7, 8_or_more]
    expression: frozen_quintiles_on_primary_eligible_universe
    chronos_missingness: preserved
  quintile_assignment: rank_first_after_sorting_by_value_then_gene_symbol
  singleton_strata: fixed_not_moved
  invalid_draw_policy: deterministic_next_seeded_draw_until_accepted_count
binary_model:
  estimates: [risk_difference, odds_ratio]
  risk_difference_estimator: shared_strata_size_weighted_class_minus_neither_difference_in_proportions
  odds_ratio_estimator: Mantel_Haenszel_over_strata_containing_both_groups
  inference: stratified_label_permutation
  confirmatory_statistic: absolute_risk_difference
  permutation_unit: evidence_class_label_within_frozen_stratum_preserving_class_counts
  zero_cell_policy: no_class_deletion_or_hit_continuity_filter_permutation_p_remains_primary
  permutation_pvalue: "(1 + count(abs(null_RD) >= abs(observed_RD))) / (accepted_draws + 1)"
  accepted_draws: 100000
  seed: 2026080741
program_model:
  gene_outcome_transform: zscore_within_primary_eligible_universe
  zscore_definition: subtract_universe_mean_divide_by_universe_sample_SD
  observed_statistic: frozen_stage_direction_times_weighted_mean_gene_zscore
  effect: observed_signed_statistic_minus_mean_matched_null_signed_statistic
  null_draws: 10000
  seed: 2026080742
  null_pool: primary_eligible_genes_excluding_members_of_the_program_being_tested
  exact_match_cells:
    guide_count: [4, 5, 6, 7, 8_or_more]
    expression: frozen_quintiles
    HLF_Chronos: frozen_quintiles
    chronos_missingness: exact
  draw_rule: one_control_per_program_member_without_replacement_within_draw
  matched_control_weight: inherit_the_corresponding_program_member_frozen_weight
  control_reuse_between_draws: allowed
  empirical_pvalue: "(1 + count(null_signed_statistic >= observed_signed_statistic)) / (null_draws + 1)"
  insufficient_exact_cell_policy: program_untestable_no_relaxation
multiplicity:
  class_omnibus: one_primary_test_no_adjustment
  class_continuous_pairwise: BH_across_three
  class_binary_pairwise: BH_across_three
  programs: BH_across_complete_testable_external_hepatocyte_family
known_hit_sensitivity:
  canonical_symbols: [ACSL3, INSIG1, NF2, CASP8, GPAT3, RNF213]
  alias_note: AGPAT9_is_GPAT3_in_GENCODE_v49
  primary_includes_known_hits: true
  mandatory_remove_all_known_hits: true
  mandatory_remove_ACSL3_only: true
figure_gate:
  alpha: 0.05
  class_minimum_standardized_pairwise_effect: 0.20
  class_magnitude_metric: adjusted_class_coefficient_divided_by_full_model_residual_SD
  program_minimum_signed_standardized_effect: 0.20
  program_magnitude_metric: observed_signed_statistic_minus_mean_matched_null_signed_statistic
  qualifying_direction_class: positive_class_minus_neither_KO_survival_compatibility
  qualifying_direction_program: positive_after_multiplying_by_frozen_stage_direction
  require_known_hit_exclusion_survival: true
  require_program_weight_sensitivity_survival_when_program_branch_used: true
  main_figure_default_if_any_component_fails: false
mandatory_sensitivities:
  - guide_floor_5
  - missing_Chronos_indicator_plus_within_missingness_stratified_permutation
  - remove_all_publication_known_hits
  - remove_ACSL3_only
  - PA_D21_vs_D0_combined_survival_selection
  - absolute_LFC_nondirectional_exploratory
  - positive_selection_FDR_hit_without_LFC_floor
missing_chronos_sensitivity:
  population: otherwise_primary_eligible_including_missing_HLF_Chronos
  imputation: median_HLF_Chronos_from_complete_primary_universe
  added_model_term: HLF_Chronos_missing_indicator
  permutation: within_frozen_guide_expression_and_Chronos_missingness_strata
  interpretation: sensitivity_only_primary_remains_complete_case
prohibited:
  - outcome_driven_gene_universe_tuning
  - use_goodsgrna_rank_p_or_FDR_as_detectability_covariates
  - reconstruct_vehicle_fitness_by_subtracting_MAGeCK_contrasts
  - analyze_CRISPR_effects_as_expression_profiles
  - call_negative_LFC_source_significant_sensitizing_without_negative_selection_statistics
  - modify_classes_programs_thresholds_or_significance_families_after_unseal
"""


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--outdir", type=Path, default=CANDIDATE_ROOT)
    parser.add_argument("--sealed-at-utc", default="")
    args = parser.parse_args()
    outdir = require_within(args.outdir, CANDIDATE_ROOT)
    if (outdir / "UNSEALED").exists():
        raise RuntimeError("Refusing to alter a specification after UNSEALED")
    unseal_record = outdir / "unseal_record.tsv"
    if unseal_record.exists() and len(read_tsv(unseal_record)) > 0:
        raise RuntimeError("Refusing to alter a specification with a populated unseal record")

    verify_hash(CLASS_PATH, "evidence_class_table.tsv")
    verify_hash(HOTSPOT_ROOT / "READY", "READY")
    verify_hash(HOTSPOT_ROOT / "external_test_programs.tsv", "external_test_programs.tsv")
    verify_hash(HOTSPOT_ROOT / "program_membership_v2.tsv", "program_membership_v2.tsv")
    verify_hash(GENCODE_PATH, "gencode_v49_gene_metadata.tsv.gz")

    phase_a_ready = outdir / "PHASE_A_READY"
    masked_path = outdir / "masked_screen_schema.tsv"
    if not phase_a_ready.is_file() or not masked_path.is_file():
        raise FileNotFoundError("Phase-A masked source artifacts are required")
    phase_a = read_tsv(phase_a_ready)
    if len(phase_a) != 1 or phase_a[0]["outcome_values_exported"] != "FALSE":
        raise ValueError("Phase-A outcome-isolation marker failed")

    classes, class_rows = load_classes()
    programs, membership = load_programs()
    gencode = load_gencode()
    covariates, depmap_ready = load_covariates(outdir)
    covariate_path = outdir / "detectability_covariates.tsv"
    if not depmap_ready and not covariate_path.exists():
        write_tsv(
            covariate_path,
            [],
            [
                "gene_symbol",
                "depmap_entrez_id",
                "expression_entrez_id",
                "chronos_entrez_id",
                "depmap_mapping_state",
                "depmap_entrez_conflict",
                "expression_source_header",
                "chronos_source_header",
                "HLF_expression_log2_tpm_plus_1",
                "HLF_TPM",
                "log1p_HLF_TPM",
                "HLF_Chronos",
                "expression_mapped",
                "chronos_mapped",
                "depmap_release",
                "model_id",
            ],
        )
    masked = read_tsv(masked_path)
    raw_counts = collections.Counter(row["source_gene_id"] for row in masked)

    mapping_rows = []
    universe_rows = []
    coverage_denominator = collections.Counter()
    coverage_numerator = collections.Counter()
    coverage_expression_floor = collections.Counter()
    eligible_symbols: set[str] = set()
    mapped_seen = collections.Counter()
    staged = []
    for row in masked:
        parsed, canonical, mapping_state, ensembl = map_source_gene(
            row["source_gene_id"], gencode
        )
        mapped_seen[canonical] += 1
        class_row = classes.get(canonical)
        evidence_class = "" if class_row is None else class_row["primary_evidence_class"]
        class_valid = evidence_class in FOUR_CLASSES
        protein_coding = mapping_state in {
            "unique_protein_coding_exact",
            "unique_protein_coding_alias",
        }
        guide_values = []
        for field in ("guide_count_d21_d0", "guide_count_pa_vehicle"):
            if row[field] != "":
                guide_values.append(int(row[field]))
        guide_min = min(guide_values) if len(guide_values) == 2 else None
        primary_fields = parse_bool(row["pa_vehicle_primary_fields_complete"])
        neutral_source_base = (
            protein_coding
            and class_valid
            and raw_counts[row["source_gene_id"]] == 1
            and mapped_seen[canonical] == 1
            and guide_min is not None
            and guide_min >= 4
            and primary_fields
        )
        staged.append(
            (
                row,
                parsed,
                canonical,
                mapping_state,
                ensembl,
                class_row,
                evidence_class,
                guide_min,
                neutral_source_base,
            )
        )

    # A second pass makes canonical-symbol duplication fail closed regardless of
    # row order.
    for (
        row,
        parsed,
        canonical,
        mapping_state,
        ensembl,
        class_row,
        evidence_class,
        guide_min,
        neutral_source_base,
    ) in staged:
        duplicate_mapping = mapped_seen[canonical] > 1
        cov = covariates.get(canonical, {})
        depmap_mapping_state = cov.get(
            "depmap_mapping_state", "missing_from_version_matched_depmap"
        )
        covariate_state = depmap_covariate_state(cov)
        depmap_entrez_conflict = bool(covariate_state["entrez_conflict"])
        tpm = covariate_state["tpm"]
        chronos = covariate_state["chronos"]
        tpm_available = bool(covariate_state["tpm_available"])
        chronos_available = bool(covariate_state["chronos_available"])
        covariates_complete = bool(covariate_state["complete"])
        expression_pass = bool(covariate_state["expression_floor_pass"])
        base = neutral_source_base and not duplicate_mapping
        if base and evidence_class in FOUR_CLASSES:
            coverage_denominator[evidence_class] += 1
            if covariates_complete:
                coverage_numerator[evidence_class] += 1
            if covariates_complete and expression_pass:
                coverage_expression_floor[evidence_class] += 1
        primary_eligible = (
            base
            and covariates_complete
            and expression_pass
            and depmap_ready
        )
        if primary_eligible:
            eligible_symbols.add(canonical)
        reasons = []
        if mapping_state not in {"unique_protein_coding_exact", "unique_protein_coding_alias"}:
            reasons.append(mapping_state)
        if duplicate_mapping:
            reasons.append("duplicate_canonical_screen_mapping")
        if class_row is None or evidence_class not in FOUR_CLASSES:
            reasons.append("not_in_primary_four_class_universe")
        if guide_min is None:
            reasons.append("missing_guide_count")
        elif guide_min < 4:
            reasons.append("guide_count_below_4")
        if not parse_bool(row["pa_vehicle_primary_fields_complete"]):
            reasons.append("direct_result_fields_incomplete")
        if not depmap_ready:
            reasons.append("version_matched_depmap_triplet_missing")
        else:
            if depmap_entrez_conflict:
                reasons.append("depmap_entrez_conflict_values_blanked")
            if not tpm_available:
                reasons.append("HLF_TPM_missing")
            elif not expression_pass:
                reasons.append("HLF_TPM_below_1")
            if not chronos_available:
                reasons.append("HLF_Chronos_missing")
        mapping_rows.append(
            {
                "source_row": row["source_row"],
                "source_gene_id": row["source_gene_id"],
                "parsed_symbol": parsed,
                "canonical_symbol": canonical,
                "mapping_state": mapping_state,
                "gencode_ensembl_base": ensembl,
                "canonical_screen_mapping_count": mapped_seen[canonical],
                "primary_evidence_class": evidence_class,
                "depmap_entrez_id": cov.get("depmap_entrez_id", ""),
                "expression_entrez_id": cov.get("expression_entrez_id", ""),
                "chronos_entrez_id": cov.get("chronos_entrez_id", ""),
                "depmap_mapping_state": depmap_mapping_state,
                "depmap_entrez_conflict": str(depmap_entrez_conflict).upper(),
                "mapping_eligible": str(
                    mapping_state
                    in {"unique_protein_coding_exact", "unique_protein_coding_alias"}
                    and not duplicate_mapping
                    and class_row is not None
                    and evidence_class in FOUR_CLASSES
                    and not depmap_entrez_conflict
                ).upper(),
            }
        )
        universe_rows.append(
            {
                "source_row": row["source_row"],
                "source_gene_id": row["source_gene_id"],
                "gene_symbol": canonical,
                "gencode_ensembl_base": ensembl,
                "mapping_state": mapping_state,
                "gene_biotype": "" if class_row is None else class_row["gene_biotype"],
                "primary_evidence_class": evidence_class,
                "depmap_entrez_id": cov.get("depmap_entrez_id", ""),
                "expression_entrez_id": cov.get("expression_entrez_id", ""),
                "chronos_entrez_id": cov.get("chronos_entrez_id", ""),
                "depmap_mapping_state": depmap_mapping_state,
                "depmap_entrez_conflict": str(depmap_entrez_conflict).upper(),
                "neutral_source_base": str(base).upper(),
                "primary_covariates_numeric_complete": str(covariates_complete).upper(),
                "HLF_expression_floor_pass": str(expression_pass).upper(),
                "guide_count_d21_d0": row["guide_count_d21_d0"],
                "guide_count_pa_vehicle": row["guide_count_pa_vehicle"],
                "guide_count_min": "" if guide_min is None else guide_min,
                "direct_primary_fields_complete": row["pa_vehicle_primary_fields_complete"],
                "HLF_TPM": "" if tpm is None else f"{tpm:.17g}",
                "log1p_HLF_TPM": cov.get("log1p_HLF_TPM", ""),
                "HLF_Chronos": "" if chronos is None else f"{chronos:.17g}",
                "actual_vehicle_fitness_observed": "FALSE",
                "primary_eligible": str(primary_eligible).upper(),
                "guide5_sensitivity_eligible": str(
                    primary_eligible and guide_min is not None and guide_min >= 5
                ).upper(),
                "exclusion_reasons": ";".join(reasons) if reasons else "eligible",
                "outcome_values_read": "FALSE",
            }
        )

    mapping_path = outdir / "gene_mapping_audit.tsv"
    write_tsv(
        mapping_path,
        mapping_rows,
        [
            "source_row",
            "source_gene_id",
            "parsed_symbol",
            "canonical_symbol",
            "mapping_state",
            "gencode_ensembl_base",
            "canonical_screen_mapping_count",
            "primary_evidence_class",
            "depmap_entrez_id",
            "expression_entrez_id",
            "chronos_entrez_id",
            "depmap_mapping_state",
            "depmap_entrez_conflict",
            "mapping_eligible",
        ],
    )
    universe_path = outdir / "screen_gene_universe.tsv"
    write_tsv(
        universe_path,
        universe_rows,
        [
            "source_row",
            "source_gene_id",
            "gene_symbol",
            "gencode_ensembl_base",
            "mapping_state",
            "gene_biotype",
            "primary_evidence_class",
            "depmap_entrez_id",
            "expression_entrez_id",
            "chronos_entrez_id",
            "depmap_mapping_state",
            "depmap_entrez_conflict",
            "neutral_source_base",
            "primary_covariates_numeric_complete",
            "HLF_expression_floor_pass",
            "guide_count_d21_d0",
            "guide_count_pa_vehicle",
            "guide_count_min",
            "direct_primary_fields_complete",
            "HLF_TPM",
            "log1p_HLF_TPM",
            "HLF_Chronos",
            "actual_vehicle_fitness_observed",
            "primary_eligible",
            "guide5_sensitivity_eligible",
            "exclusion_reasons",
            "outcome_values_read",
        ],
    )

    prediction_rows = []
    known_symbols = {canonical for _, canonical, _ in KNOWN_HITS}
    for class_row in sorted(class_rows, key=lambda row: row["symbol"]):
        if (
            class_row["gene_biotype"] != "protein_coding"
            or class_row["primary_evidence_class"] not in FOUR_CLASSES
            or not class_row["symbol"]
        ):
            continue
        prediction_rows.append(
            {
                "object_type": "evidence_class_gene",
                "object_id": class_row["primary_evidence_class"],
                "gene_symbol": class_row["symbol"],
                "primary_evidence_class": class_row["primary_evidence_class"],
                "expected_direction": expected_gene_direction(class_row),
                "frozen_weight": "1",
                "frozen_source_effect": class_row["bulk_logFC"],
                "source_effect_unit": "canonical_bulk_log2FC",
                "screen_primary_eligible": str(class_row["symbol"] in eligible_symbols).upper(),
                "publication_known_hit": str(class_row["symbol"] in known_symbols).upper(),
                "registry_sha256": "not_applicable",
                "membership_sha256": "not_applicable",
                "external_outcomes_read": "FALSE",
            }
        )
    program_by_id = {row["program_uid"]: row for row in programs}
    for program_uid in sorted(membership):
        program = program_by_id[program_uid]
        for member in membership[program_uid]:
            symbol = str(member["gene_symbol"])
            prediction_rows.append(
                {
                    "object_type": "hepatocyte_program_gene",
                    "object_id": program_uid,
                    "gene_symbol": symbol,
                    "primary_evidence_class": classes.get(symbol, {}).get(
                        "primary_evidence_class", "not_applicable"
                    ),
                    "expected_direction": (
                        "positive_KO_LFC_counter_state_compatible"
                        if float(program["primary_beta"]) > 0
                        else "negative_KO_LFC_counter_state_compatible"
                    ),
                    "frozen_weight": f"{float(member['original_l1_weight']):.17g}",
                    "frozen_source_effect": program["primary_beta"],
                    "source_effect_unit": "donor_stage_ordinal_beta",
                    "screen_primary_eligible": str(symbol in eligible_symbols).upper(),
                    "publication_known_hit": str(symbol in known_symbols).upper(),
                    "registry_sha256": program["registry_sha256"],
                    "membership_sha256": program["membership_sha256"],
                    "external_outcomes_read": "FALSE",
                }
            )
    prediction_path = outdir / "prediction_manifest.tsv"
    write_tsv(
        prediction_path,
        prediction_rows,
        [
            "object_type",
            "object_id",
            "gene_symbol",
            "primary_evidence_class",
            "expected_direction",
            "frozen_weight",
            "frozen_source_effect",
            "source_effect_unit",
            "screen_primary_eligible",
            "publication_known_hit",
            "registry_sha256",
            "membership_sha256",
            "external_outcomes_read",
        ],
    )

    known_path = outdir / "known_hit_exclusion.tsv"
    write_tsv(
        known_path,
        [
            {
                "source_name": source,
                "canonical_symbol": canonical,
                "freeze_basis": basis,
                "included_in_primary": "TRUE",
                "remove_all_known_hits_sensitivity": "TRUE",
                "ACSL3_only_sensitivity": str(canonical == "ACSL3").upper(),
                "external_outcomes_read": "FALSE",
            }
            for source, canonical, basis in KNOWN_HITS
        ],
        [
            "source_name",
            "canonical_symbol",
            "freeze_basis",
            "included_in_primary",
            "remove_all_known_hits_sensitivity",
            "ACSL3_only_sensitivity",
            "external_outcomes_read",
        ],
    )

    program_testability_rows = []
    for program in programs:
        members = membership[program["program_uid"]]
        eligible_members = [
            member for member in members if str(member["gene_symbol"]) in eligible_symbols
        ]
        retained = sum(float(member["original_l1_weight"]) for member in eligible_members)
        testable = depmap_ready and len(eligible_members) >= 8 and retained >= 0.20
        reasons = []
        if not depmap_ready:
            reasons.append("depmap_source_gate_blocked")
        if len(eligible_members) < 8:
            reasons.append("fewer_than_8_primary_eligible_genes")
        if retained < 0.20:
            reasons.append("retained_L1_weight_below_0.20")
        program_testability_rows.append(
            {
                "program_uid": program["program_uid"],
                "module_name": program["module_name"],
                "registry_sha256": program["registry_sha256"],
                "membership_sha256": program["membership_sha256"],
                "n_frozen_unique_genes": len(members),
                "n_primary_eligible_genes": len(eligible_members),
                "retained_original_L1_weight": f"{retained:.17g}",
                "testable_if_unsealed_now": str(testable).upper(),
                "testability_reason": "testable" if testable else ";".join(reasons),
                "external_outcomes_read": "FALSE",
            }
        )
    testability_path = outdir / "program_testability.tsv"
    write_tsv(
        testability_path,
        program_testability_rows,
        [
            "program_uid",
            "module_name",
            "registry_sha256",
            "membership_sha256",
            "n_frozen_unique_genes",
            "n_primary_eligible_genes",
            "retained_original_L1_weight",
            "testable_if_unsealed_now",
            "testability_reason",
            "external_outcomes_read",
        ],
    )

    coverage_audit_path = outdir / "covariate_coverage_audit.tsv"
    coverage_audit_rows = []
    for class_name in FOUR_CLASSES:
        denominator = coverage_denominator[class_name]
        complete = coverage_numerator[class_name]
        detectable = coverage_expression_floor[class_name]
        coverage_audit_rows.append(
            {
                "primary_evidence_class": class_name,
                "neutral_source_base_n": denominator,
                "numeric_TPM_and_Chronos_complete_n": complete,
                "numeric_covariate_coverage_fraction": (
                    "" if denominator == 0 else f"{complete / denominator:.17g}"
                ),
                "TPM_ge_1_primary_eligible_n": detectable,
                "TPM_ge_1_fraction_of_neutral_base": (
                    "" if denominator == 0 else f"{detectable / denominator:.17g}"
                ),
                "expression_floor_attrition_n": complete - detectable,
                "coverage_gate_threshold": "0.80",
                "coverage_gate_pass": str(
                    denominator > 0 and complete / denominator >= 0.80
                ).upper(),
                "coverage_gate_uses_expression_floor": "FALSE",
                "external_outcomes_read": "FALSE",
                "release_id": RELEASE_ID,
            }
        )
    write_tsv(
        coverage_audit_path,
        coverage_audit_rows,
        [
            "primary_evidence_class",
            "neutral_source_base_n",
            "numeric_TPM_and_Chronos_complete_n",
            "numeric_covariate_coverage_fraction",
            "TPM_ge_1_primary_eligible_n",
            "TPM_ge_1_fraction_of_neutral_base",
            "expression_floor_attrition_n",
            "coverage_gate_threshold",
            "coverage_gate_pass",
            "coverage_gate_uses_expression_floor",
            "external_outcomes_read",
            "release_id",
        ],
    )

    nonreference_coverage_pass = all(
        coverage_denominator[class_name] > 0
        and coverage_numerator[class_name] / coverage_denominator[class_name] >= 0.80
        for class_name in ("genetic_only", "disease_state_only", "convergent")
    )
    mapping_gate = all(mapped_seen[symbol] == 1 for symbol in eligible_symbols)
    depmap_conflict_symbols = {
        symbol
        for symbol, row in covariates.items()
        if parse_bool(row.get("depmap_entrez_conflict", "FALSE"))
    }
    depmap_conflicts_excluded = depmap_ready and all(
        symbol not in eligible_symbols for symbol in depmap_conflict_symbols
    )
    expression_floor_audit_pass = depmap_ready and all(
        coverage_denominator[class_name] > 0
        and coverage_expression_floor[class_name] <= coverage_numerator[class_name]
        for class_name in FOUR_CLASSES
    )
    gate_rows = []
    gate_values = {
        "masked_source_ready": True,
        "external_outcomes_unread": True,
        "frozen_evidence_class_hash": True,
        "frozen_hotspot_registry_hashes": True,
        "gencode_v49_hash": True,
        "version_matched_depmap_triplet": depmap_ready,
        "HLF_identity_unique_ACH_000393": depmap_ready,
        "HLF_expression_available": depmap_ready,
        "depmap_entrez_conflicts_explicitly_excluded": depmap_conflicts_excluded,
        "mapping_unambiguous_for_primary_universe": mapping_gate,
        "HLF_TPM_ge_1_eligibility_attrition_reported": expression_floor_audit_pass,
        "nonreference_class_complete_covariates_at_least_80pct": (
            depmap_ready and nonreference_coverage_pass
        ),
    }
    for gate, passed in gate_values.items():
        detail = "pass"
        if gate == "version_matched_depmap_triplet" and not passed:
            detail = "awaiting_DepMap_24Q4_Model_Chronos_expression_triplet"
        elif gate == "HLF_identity_unique_ACH_000393" and not passed:
            detail = "identity_known_locally_but_not_accepted_without_version_matched_triplet"
        elif gate == "HLF_expression_available" and not passed:
            detail = "no_version_matched_public_HLF_expression_matrix_present"
        elif gate == "depmap_entrez_conflicts_explicitly_excluded":
            detail = (
                "excluded_symbols=" + ",".join(sorted(depmap_conflict_symbols))
                if passed
                else "DepMap_conflicts_missing_or_not_fail_closed"
            )
        elif gate == "HLF_TPM_ge_1_eligibility_attrition_reported":
            detail = ";".join(
                f"{name}:{coverage_expression_floor[name]}/{coverage_denominator[name]}"
                for name in FOUR_CLASSES
            )
        elif gate == "nonreference_class_complete_covariates_at_least_80pct":
            detail = "complete=" + ";".join(
                f"{name}:{coverage_numerator[name]}/{coverage_denominator[name]}"
                for name in ("genetic_only", "disease_state_only", "convergent")
            ) + "|TPM_ge_1=" + ";".join(
                f"{name}:{coverage_expression_floor[name]}/{coverage_denominator[name]}"
                for name in ("genetic_only", "disease_state_only", "convergent")
            )
        gate_rows.append(
            {
                "gate": gate,
                "passed": str(passed).upper(),
                "detail": detail,
                "external_outcomes_read": "FALSE",
                "release_id": RELEASE_ID,
            }
        )
    sealed = all(gate_values.values())
    gate_path = outdir / "gate_status.tsv"
    write_tsv(
        gate_path,
        gate_rows,
        ["gate", "passed", "detail", "external_outcomes_read", "release_id"],
    )

    spec_path = outdir / "blind_analysis_spec.yaml"
    atomic_write_text(spec_path, blind_spec_yaml(sealed, depmap_ready))
    deviations_path = outdir / "deviations.tsv"
    write_tsv(
        deviations_path,
        [],
        [
            "deviation_id",
            "recorded_at_utc",
            "category",
            "original_rule",
            "deviation",
            "reason",
            "confirmatory_status_invalidated",
            "new_release_id_required",
        ],
    )
    write_tsv(
        unseal_record,
        [],
        [
            "unsealed_at_utc",
            "specification_sha256",
            "executor",
            "raw_outcome_source_sha256",
            "one_pass_release_id",
        ],
    )

    producer_paths = [
        Path(__file__).resolve(),
        Path(__file__).resolve().parent / "00_mask_myojin_workbook.py",
        Path(__file__).resolve().parent / "01_extract_depmap_hlf.py",
        Path(__file__).resolve().parent / "03_validate_blind_spec.py",
        Path(__file__).resolve().parent / "myojin_firewall_common.py",
    ]
    input_paths = [
        CLASS_PATH,
        HOTSPOT_ROOT / "READY",
        HOTSPOT_ROOT / "external_test_programs.tsv",
        HOTSPOT_ROOT / "program_membership_v2.tsv",
        GENCODE_PATH,
        masked_path,
        outdir / "source_manifest.tsv",
        outdir / "source_schema_audit.tsv",
    ]
    local_preflight_path = outdir / "local_depmap_preflight.tsv"
    if local_preflight_path.is_file():
        input_paths.append(local_preflight_path)
    if depmap_ready:
        input_paths.extend(
            [
                outdir / "detectability_covariates.tsv",
                outdir / "depmap_gene_mapping_audit.tsv",
                outdir / "depmap_source_manifest.tsv",
                outdir / "hlf_identity_audit.tsv",
                outdir / "DEPMAP_READY",
            ]
        )
    input_manifest_path = outdir / "firewall_input_manifest.tsv"
    input_rows = []
    for role, paths in (("frozen_input", input_paths), ("producer", producer_paths)):
        for path in paths:
            input_rows.append(
                {
                    "role": role,
                    "path": str(path.relative_to(PROJECT_ROOT)),
                    "bytes": path.stat().st_size,
                    "sha256": sha256_file(path),
                    "external_outcomes_read": "FALSE",
                    "release_id": RELEASE_ID,
                }
            )
    write_tsv(
        input_manifest_path,
        input_rows,
        ["role", "path", "bytes", "sha256", "external_outcomes_read", "release_id"],
    )

    specification_components = [
        spec_path,
        prediction_path,
        known_path,
        universe_path,
        mapping_path,
        coverage_audit_path,
        testability_path,
        gate_path,
        input_manifest_path,
    ]
    specification_hash, component_rows = stable_bundle_sha256(
        specification_components, outdir
    )
    specification_hash_path = outdir / "specification_sha256.txt"
    text = f"specification_bundle_sha256\t{specification_hash}\n"
    text += "".join(
        f"component\t{row['relative_path']}\t{row['bytes']}\t{row['sha256']}\n"
        for row in component_rows
    )
    atomic_write_text(specification_hash_path, text)

    sealed_at = args.sealed_at_utc or dt.datetime.now(dt.timezone.utc).strftime(
        "%Y-%m-%dT%H:%M:%SZ"
    )
    status_path = outdir / ("SEALED" if sealed else "BLOCKED")
    other_status = outdir / ("BLOCKED" if sealed else "SEALED")
    if other_status.exists():
        other_status.unlink()
    atomic_write_text(
        status_path,
        "release_id\tstatus\tspecification_sha256\tsealed_at_utc\texternal_outcomes_read\n"
        f"{RELEASE_ID}\t{'sealed_ready_for_unseal' if sealed else 'blocked_source_or_covariate_gate'}\t"
        f"{specification_hash}\t{sealed_at}\tFALSE\n",
    )
    release_artifacts = [
        phase_a_ready,
        outdir / "source_manifest.tsv",
        outdir / "source_schema_audit.tsv",
        masked_path,
        covariate_path,
        mapping_path,
        coverage_audit_path,
        universe_path,
        prediction_path,
        known_path,
        testability_path,
        spec_path,
        specification_hash_path,
        deviations_path,
        unseal_record,
        gate_path,
        input_manifest_path,
        status_path,
    ]
    if depmap_ready:
        release_artifacts.extend(
            [
                outdir / "depmap_gene_mapping_audit.tsv",
                outdir / "depmap_source_manifest.tsv",
                outdir / "hlf_identity_audit.tsv",
                outdir / "DEPMAP_READY",
            ]
        )
    if local_preflight_path.is_file():
        release_artifacts.append(local_preflight_path)
    release_manifest_path = outdir / "release_manifest.tsv"
    release_rows = []
    for path in release_artifacts:
        with path.open(encoding="utf-8", errors="replace") as handle:
            line_count = sum(1 for _ in handle)
        release_rows.append(
            {
                "role": "blocked_specification_artifact" if not sealed else "sealed_specification_artifact",
                "relative_path": str(path.relative_to(outdir)),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
                "line_count": line_count,
                "external_outcomes_read": "FALSE",
                "canonical": "FALSE",
                "release_id": RELEASE_ID,
            }
        )
    write_tsv(
        release_manifest_path,
        release_rows,
        [
            "role",
            "relative_path",
            "bytes",
            "sha256",
            "line_count",
            "external_outcomes_read",
            "canonical",
            "release_id",
        ],
    )
    print(
        f"SPECIFICATION_{'SEALED' if sealed else 'BLOCKED'} "
        f"masked_rows={len(masked)} primary_eligible={len(eligible_symbols)} "
        f"programs={len(programs)} depmap_ready={depmap_ready} "
        f"specification_sha256={specification_hash} external_outcomes_read=FALSE"
    )


if __name__ == "__main__":
    main()
