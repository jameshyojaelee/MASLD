#!/usr/bin/env python3
"""Seal the outcome-free Plan 46 human analysis contract before cohort access."""

from __future__ import annotations

import csv
import datetime as dt
import hashlib
import json
import os
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parents[4]
SCRIPT_ROOT = Path(__file__).resolve().parent
PLAN43 = PROJECT_ROOT / "docs/plans/2026-08-07_paper_program/43_CHRONIC_STATE_REVERSAL_AND_REGULATORY_RISK_BRIDGE.md"
PLAN46 = PROJECT_ROOT / "docs/plans/2026-08-07_paper_program/46_FULL_GOAL_COMPLETION_NEW_DATA_PROTOCOL.md"
PLAN46B = PROJECT_ROOT / "docs/plans/2026-08-07_paper_program/46B_CLINICAL_RESPONSE_INTAKE_AND_TWO_COHORT_GATE.md"
PLAN46C = PROJECT_ROOT / "docs/plans/2026-08-07_paper_program/46C_BLINDED_CLINICAL_PREFLIGHT_AND_MODEL_FREEZE.md"
STATE_ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/chronic-state-risk-bridge-2026-08-09"
STATE_AXIS = STATE_ROOT / "frozen_state_axis.tsv"
STATE_SPEC = STATE_ROOT / "frozen_spec/state_axis_spec.tsv"
ROOT = PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates/source-independent-risk-state-relay-clinical-preflight-template-v1-2026-08-10"

EXPECTED_AXIS_SHA256 = "03242574b4a4aee810000cd2304abb0ba9f3c960ee7b2e0f584c3e3c8a180e33"
EXPECTED_SPEC_SHA256 = "beaab7cc62a413dc5532159a757dd17ca3e6315f8bfa39b4188f8adccc1a2575"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def make_read_only(root: Path) -> None:
    for path in sorted(root.rglob("*"), reverse=True):
        path.chmod(0o555 if path.is_dir() else 0o444)
    root.chmod(0o555)


RESPONDER_ROWS = [
    {
        "endpoint_id": "primary_histologic_improvement",
        "population": "baseline_biopsy_defined_MASH_with_complete_paired_central_histology",
        "improver_definition": "NAS_decrease_ge_2_with_ge_1_point_ballooning_or_lobular_inflammation_improvement_and_no_fibrosis_worsening",
        "nonimprover_definition": "complete_eligible_pair_not_meeting_primary_improver_definition",
        "direction_for_delta_state": "improver_minus_nonimprover_less_than_0",
        "role": "primary_confirmatory",
        "substitution_allowed": "false",
    },
    {
        "endpoint_id": "mash_resolution_no_fibrosis_worsening",
        "population": "baseline_biopsy_defined_MASH_with_complete_paired_central_histology",
        "improver_definition": "followup_ballooning_0_and_lobular_inflammation_0_or_1_and_no_fibrosis_worsening",
        "nonimprover_definition": "complete_eligible_pair_not_meeting_resolution_definition",
        "direction_for_delta_state": "improver_minus_nonimprover_less_than_0",
        "role": "prespecified_secondary",
        "substitution_allowed": "false",
    },
    {
        "endpoint_id": "fibrosis_improvement",
        "population": "complete_paired_central_histology",
        "improver_definition": "fibrosis_stage_decrease_ge_1_without_worsening_ballooning_or_lobular_inflammation",
        "nonimprover_definition": "complete_eligible_pair_not_meeting_fibrosis_improvement_definition",
        "direction_for_delta_state": "improver_minus_nonimprover_less_than_0",
        "role": "prespecified_secondary",
        "substitution_allowed": "false",
    },
]


ESTIMANDS = [
    {
        "estimand_id": "HUM_PRIMARY",
        "analysis_unit": "participant",
        "outcome": "delta_state_axis_followup_minus_baseline",
        "contrast": "mean_delta_in_primary_improvers_minus_mean_delta_in_nonimprovers",
        "model": "unadjusted_studentized_difference_with_responder_labels_permuted_within_frozen_trial_and_intervention_strata",
        "expected_direction": "negative",
        "inference": "exact_enumeration_when_feasible_else_999999_seeded_permutations",
        "multiplicity": "one_primary_test_per_cohort",
        "promotion_role": "each_cohort_must_independently_pass_two_sided_p_lt_0.05_in_expected_direction",
    },
    {
        "estimand_id": "HUM_COVARIATE_SENSITIVITY",
        "analysis_unit": "participant",
        "outcome": "delta_state_axis_followup_minus_baseline",
        "contrast": "primary_improver_coefficient",
        "model": "robust_regression_with_frozen_intervention_delta_BMI_interval_baseline_NAS_baseline_fibrosis_age_sex_covariates_when_deposited",
        "expected_direction": "negative",
        "inference": "heteroskedasticity_robust_interval_and_Freedman_Lane_residual_permutation",
        "multiplicity": "sensitivity_not_primary_family",
        "promotion_role": "direction_and_material_effect_stability_required",
    },
    {
        "estimand_id": "HUM_TWO_COHORT_META",
        "analysis_unit": "cohort",
        "outcome": "signed_primary_cohort_Z",
        "contrast": "equal_cohort_signed_Stouffer",
        "model": "two_locked_cohorts_equal_weight",
        "expected_direction": "negative",
        "inference": "two_sided_p_lt_0.01_and_both_individual_p_lt_0.05",
        "multiplicity": "single_frozen_meta_test",
        "promotion_role": "cannot_rescue_failed_or_retuned_cohort",
    },
]


POWER_POLICIES = [
    ("biological_unit", "participant", "biopsies_arrays_sections_and_assays_never_inflate_n"),
    ("minimum_class_floor", "at_least_5_primary_improvers_and_5_nonimprovers", "software_floor_not_a_power_claim"),
    ("exact_test_resolution", "conservative_two_sided_minimum_p_le_0.05", "derive_from_frozen_class_counts_before_molecular_outcomes"),
    ("target_power", "at_least_0.80", "requires_recorded_alpha_effect_and_variance_source"),
    ("variance_source", "external_or_response_blinded_only", "frozen_state_scores_and_responder_linked_molecular_values_prohibited"),
    ("attrition", "predeclared_from_assay_QC_not_state_score", "no_outcome_selected_exclusion"),
    ("underpowered_action", "fail_source_gate_or_later_manuscript", "no_responder_or_cohort_substitution"),
]


QC_POLICIES = [
    ("QC01", "authoritative_pair_key", "required", "never_infer_pairs_from_sample_order_title_geometry_or_similarity"),
    ("QC02", "one_profile_per_biopsy", "required", "technical_replicates_collapse_before_participant_inference"),
    ("QC03", "RNA_quality_and_library_QC", "thresholds_frozen_before_response_key_release", "no_program_or_state_score_based_exclusion"),
    ("QC04", "baseline_standardization", "baseline_samples_only_within_cohort", "followup_samples_never_change_gene_mean_or_SD"),
    ("QC05", "program_testability", "ge_8_mapped_genes_and_ge_0.20_retained_original_L1_weight", "frozen_Plan43_rule"),
    ("QC06", "dataset_testability", "ge_90_programs_and_ge_0.80_absolute_axis_loading", "failure_is_untestable_not_negative"),
    ("QC07", "participant_exclusion", "dual_reviewed_source_or_technical_reason_only", "reason_and_timestamp_precede_state_projection"),
    ("QC08", "cohort_B_lock", "complete_protocol_hash_precedes_cohort_A_outcome_release", "no_post_discovery_change"),
]


ORTHOGONAL_ROWS = [
    {
        "assay_id": "tissue_proteomics",
        "biological_unit": "participant",
        "primary_object": "paired_observable_protein_state_projection",
        "mandatory_adjustments": "abundance_coverage_and_continuous_RNA_effect",
        "expected_direction": "improvement_associated_reversal",
        "native_null": "participant_label_or_residual_permutation",
        "substitution_allowed": "false",
    },
    {
        "assay_id": "liver_linked_secreted_proteomics",
        "biological_unit": "participant",
        "primary_object": "paired_liver_attributed_secreted_protein_projection",
        "mandatory_adjustments": "systemic_source_annotation_processing_batch_and_abundance",
        "expected_direction": "improvement_associated_reversal",
        "native_null": "participant_label_or_residual_permutation",
        "substitution_allowed": "false",
    },
    {
        "assay_id": "physical_niche",
        "biological_unit": "participant",
        "primary_object": "paired_donor_collapsed_spatial_state_organization",
        "mandatory_adjustments": "zonation_composition_coverage_and_spatial_blocking",
        "expected_direction": "improvement_associated_reversal",
        "native_null": "spatially_blocked_within_participant_then_participant_level_inference",
        "substitution_allowed": "false",
    },
]


def payload_files() -> list[Path]:
    return sorted(
        path for path in ROOT.rglob("*")
        if path.is_file() and path.name not in {
            "CLINICAL_PREFLIGHT_TEMPLATE_SEALED.json",
            "CLINICAL_PREFLIGHT_TEMPLATE_SEAL_SHA256.txt",
        }
    )


def main() -> None:
    if ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite clinical preflight template: {ROOT}")
    for path in (PLAN43, PLAN46, PLAN46B, PLAN46C, STATE_AXIS, STATE_SPEC, SCRIPT_ROOT / "64_validate_blinded_clinical_preflight_template.py"):
        if not path.is_file():
            raise RuntimeError(f"Missing clinical preflight source: {path}")
    if sha256_file(STATE_AXIS) != EXPECTED_AXIS_SHA256 or sha256_file(STATE_SPEC) != EXPECTED_SPEC_SHA256:
        raise RuntimeError("Frozen Plan 43 state geometry drift")
    axis = read_tsv(STATE_AXIS)
    if len(axis) != 117 or len({row["program_uid"] for row in axis}) != 117:
        raise RuntimeError("Frozen state axis is not the complete 117-program universe")
    if abs(sum(abs(float(row["primary_loading"])) for row in axis) - 1.0) > 1e-10:
        raise RuntimeError("Frozen state-axis absolute loading does not sum to one")

    ROOT.mkdir(parents=True)
    write_tsv(ROOT / "state_axis_binding.tsv", [{
        "state_axis_path": STATE_AXIS.relative_to(PROJECT_ROOT),
        "state_axis_sha256": EXPECTED_AXIS_SHA256,
        "state_spec_path": STATE_SPEC.relative_to(PROJECT_ROOT),
        "state_spec_sha256": EXPECTED_SPEC_SHA256,
        "n_programs": 117,
        "score_formula": "sum_L1_normalized_stage_beta_times_standardized_program_score_over_testable_absolute_loading",
        "reference_standardization": "baseline_only_within_cohort",
        "minimum_testable_programs": 90,
        "minimum_retained_axis_mass": 0.80,
        "new_program_selection_allowed": "false",
    }], [
        "state_axis_path", "state_axis_sha256", "state_spec_path",
        "state_spec_sha256", "n_programs", "score_formula",
        "reference_standardization", "minimum_testable_programs",
        "minimum_retained_axis_mass", "new_program_selection_allowed",
    ])
    write_tsv(ROOT / "histologic_responder_definitions.tsv", RESPONDER_ROWS, list(RESPONDER_ROWS[0]))
    write_tsv(ROOT / "prospective_human_estimands.tsv", ESTIMANDS, list(ESTIMANDS[0]))
    write_tsv(ROOT / "power_and_resolution_policy.tsv", [
        {"policy_id": f"PW{i:02d}", "parameter": parameter, "frozen_value": value, "rationale": rationale}
        for i, (parameter, value, rationale) in enumerate(POWER_POLICIES, start=1)
    ], ["policy_id", "parameter", "frozen_value", "rationale"])
    write_tsv(ROOT / "qc_and_exclusion_policy.tsv", [
        {"policy_id": pid, "domain": domain, "frozen_rule": rule, "prohibition": prohibition}
        for pid, domain, rule, prohibition in QC_POLICIES
    ], ["policy_id", "domain", "frozen_rule", "prohibition"])
    write_tsv(ROOT / "cohort_protocol_template.tsv", [{
        "portfolio_role": role,
        "clinical_response_candidate_id": "",
        "cohort_uid": "",
        "protocol_sha256": "",
        "responder_endpoint_id": "primary_histologic_improvement",
        "primary_estimand_id": "HUM_PRIMARY",
        "state_axis_sha256": EXPECTED_AXIS_SHA256,
        "protocol_locked_utc": "",
        "outcome_key_released_utc": "",
        "molecular_outcome_accessed": "false",
        "status": "blank_requires_passing_plan46b_portfolio",
    } for role in ("cohort_A", "cohort_B")], [
        "portfolio_role", "clinical_response_candidate_id", "cohort_uid",
        "protocol_sha256", "responder_endpoint_id", "primary_estimand_id",
        "state_axis_sha256", "protocol_locked_utc", "outcome_key_released_utc",
        "molecular_outcome_accessed", "status",
    ])
    write_tsv(ROOT / "blinded_power_audit_template.tsv", [{
        "portfolio_role": role,
        "n_authoritative_pairs": "",
        "n_primary_improvers": "",
        "n_nonimprovers": "",
        "minimum_two_sided_exact_p": "",
        "target_alpha": "0.05",
        "target_power": "0.80",
        "target_effect_definition": "",
        "variance_source": "",
        "required_n": "",
        "power_gate": "pending",
        "state_axis_or_responder_linked_molecular_outcome_inspected": "false",
        "evidence_path": "",
    } for role in ("cohort_A", "cohort_B")], [
        "portfolio_role", "n_authoritative_pairs", "n_primary_improvers",
        "n_nonimprovers", "minimum_two_sided_exact_p", "target_alpha",
        "target_power", "target_effect_definition", "variance_source",
        "required_n", "power_gate",
        "state_axis_or_responder_linked_molecular_outcome_inspected",
        "evidence_path",
    ])
    write_tsv(ROOT / "outcome_lock_firewall.tsv", [
        {
            "role": "clinical_data_manager",
            "responsibility": "hold_participant_response_and_treatment_keys_until_protocol_and_QC_freeze",
            "assigned_person": "",
            "accepted_utc": "",
            "prohibition": "cannot_be_analysis_lead",
        },
        {
            "role": "analysis_lead",
            "responsibility": "freeze_code_models_covariates_multiplicity_and_output_contract",
            "assigned_person": "",
            "accepted_utc": "",
            "prohibition": "cannot_receive_response_key_before_freeze",
        },
        {
            "role": "independent_validator",
            "responsibility": "rederive_pairing_power_protocol_hashes_and_release_chronology",
            "assigned_person": "",
            "accepted_utc": "",
            "prohibition": "cannot_be_analysis_lead",
        },
        {
            "role": "histopathology_lead",
            "responsibility": "freeze_histology_provenance_and_responder_labels_before_state_projection",
            "assigned_person": "",
            "accepted_utc": "",
            "prohibition": "cannot_change_responder_rule_after_molecular_access",
        },
    ], ["role", "responsibility", "assigned_person", "accepted_utc", "prohibition"])
    write_tsv(ROOT / "orthogonal_translation_protocol.tsv", ORTHOGONAL_ROWS, list(ORTHOGONAL_ROWS[0]))
    write_tsv(ROOT / "required_preflight_inputs.tsv", [
        {"input_id": input_id, "required_object": object_name, "contains_participant_rows": contains, "contains_molecular_outcomes": outcomes, "status": "pending"}
        for input_id, object_name, contains, outcomes in (
            ("IN01", "sealed_passing_Plan46B_two_cohort_portfolio", "false", "false"),
            ("IN02", "deidentified_authoritative_pair_census_and_overlap_hashes", "aggregate_or_hashed_only", "false"),
            ("IN03", "central_histology_definition_and_blinded_class_counts", "false", "false"),
            ("IN04", "outcome_blind_RNA_and_orthogonal_assay_QC_inventory", "false", "false"),
            ("IN05", "external_or_response_blinded_power_inputs", "false", "false"),
            ("IN06", "signed_role_firewall_and_Cohort_B_lock_chronology", "false", "false"),
        )
    ], ["input_id", "required_object", "contains_participant_rows", "contains_molecular_outcomes", "status"])
    (ROOT / "README.md").write_text(
        """# Outcome-free paired-human preflight template

This immutable template freezes the Plan 46 human scoring, responder,
estimand, sensitivity, QC, power, orthogonal-assay, and replication-lock rules
before any new cohort is selected or molecular outcome is opened. It is not an
analysis authorization. It contains no participant row, partner response,
clinical outcome, expression value, state score, permission assertion, or
paper-promotion decision.

Instantiation requires a sealed passing Plan 46B two-cohort portfolio. Cohort
A and Cohort B must be frozen together; Cohort B's protocol hash must predate
Cohort A outcome release. Failure of blinded resolution or power closes or
defers the cohort and never changes the responder definition.
""",
        encoding="utf-8",
    )
    write_tsv(ROOT / "source_manifest.tsv", [
        {"role": role, "path": path.relative_to(PROJECT_ROOT), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
        for role, path in (
            ("plan43", PLAN43), ("plan46", PLAN46), ("plan46b", PLAN46B),
            ("plan46c", PLAN46C),
            ("state_axis", STATE_AXIS), ("state_spec", STATE_SPEC),
            ("producer", Path(__file__).resolve()),
            ("validator", SCRIPT_ROOT / "64_validate_blinded_clinical_preflight_template.py"),
        )
    ], ["role", "path", "sha256", "size_bytes"])
    files = payload_files()
    write_tsv(ROOT / "file_manifest.tsv", [
        {"relative_path": path.relative_to(ROOT), "sha256": sha256_file(path), "size_bytes": path.stat().st_size}
        for path in files
    ], ["relative_path", "sha256", "size_bytes"])
    outputs = payload_files()
    seal = {
        "candidate_id": ROOT.name,
        "status": "sealed_outcome_free_two_cohort_clinical_preflight_template",
        "generated_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "n_frozen_programs": 117,
        "n_cohort_roles": 2,
        "n_responder_endpoints": len(RESPONDER_ROWS),
        "n_estimands": len(ESTIMANDS),
        "n_orthogonal_assays": len(ORTHOGONAL_ROWS),
        "bundle_read_only": True,
        "real_partner_response_present": False,
        "participant_rows_present": False,
        "new_clinical_outcomes_present": False,
        "new_molecular_outcomes_present": False,
        "cohort_selected": False,
        "outcome_key_release_authorized": False,
        "analysis_authorized": False,
        "paper_promotion_authorized": False,
        "state_axis_sha256": EXPECTED_AXIS_SHA256,
        "output_sha256": {
            path.relative_to(ROOT).as_posix(): sha256_file(path) for path in outputs
        },
    }
    seal_path = ROOT / "CLINICAL_PREFLIGHT_TEMPLATE_SEALED.json"
    seal_path.write_text(json.dumps(seal, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    (ROOT / "CLINICAL_PREFLIGHT_TEMPLATE_SEAL_SHA256.txt").write_text(
        sha256_file(seal_path) + "\n", encoding="utf-8"
    )
    make_read_only(ROOT)
    print(json.dumps(seal, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
