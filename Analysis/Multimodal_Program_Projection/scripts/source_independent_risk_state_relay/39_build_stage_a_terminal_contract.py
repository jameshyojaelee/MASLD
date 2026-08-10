#!/usr/bin/env python3
"""Seal the Plan 45 Stage-A terminal adjudication and Stage-B freeze contract."""

from __future__ import annotations

import datetime as dt
import json
import os
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT, PROJECT_ROOT, SCRIPT_ROOT, atomic_write_json, sha256_file,
    write_tsv,
)


def candidate_source(variable: str) -> Path:
    raw = os.environ.get(variable, "").strip()
    if not raw:
        raise RuntimeError(f"{variable} is required")
    path = Path(raw)
    path = path if path.is_absolute() else PROJECT_ROOT / path
    path = path.resolve()
    allowed = (
        PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates"
    ).resolve()
    if allowed not in path.parents:
        raise RuntimeError(f"{variable} escapes the candidate root: {path}")
    return path


SCHEMAS = {
    "stage_a_primary_estimand_results_template.tsv": [
        "estimand_id", "target_uid", "guide_id", "endpoint_id", "effect",
        "se", "p_value", "q_value", "expected_direction_agreement",
        "n_backgrounds", "n_biological_units", "assay_native_model",
        "full_family_size", "result_status",
    ],
    "stage_a_sensitivity_results_template.tsv": [
        "sensitivity_id", "target_uid", "guide_id", "recipient_lineage",
        "effect", "expected_direction_agreement", "n_backgrounds",
        "result_status",
    ],
    "stage_a_mediator_candidates_template.tsv": [
        "target_uid", "recipient_lineage", "mediator_id", "effect", "se",
        "p_value", "q_value", "expected_direction_agreement",
        "intermediate_time_detected", "source_secretome_detected",
        "recipient_receptor_detected", "n_backgrounds", "blocking_reagent",
        "blocking_dose", "blocking_dose_unit", "blocking_feasibility_pass",
        "upstream_cis_preservation_assay", "off_target_control", "result_status",
    ],
    "stage_a_outcome_signoff_template.tsv": [
        "unblinding_seal_path", "unblinding_seal_sha256",
        "primary_results_path", "primary_results_sha256",
        "sensitivity_results_path", "sensitivity_results_sha256",
        "mediator_results_path", "mediator_results_sha256", "analyst_1",
        "analyst_2", "complete_outcome_universe", "signed_utc",
    ],
    "stage_b_exact_design_template.tsv": [
        "target_uid", "recipient_lineage", "mediator_id", "variant_id",
        "snp_hg19", "orientation_uid", "exact_design_uid", "genome_build",
        "protective_allele", "risk_allele", "editing_method", "pegrna_id_1",
        "pegrna_id_2", "edit_design_sha256", "source_lineage", "challenge_id",
        "early_time_value", "early_time_unit", "intermediate_time_value",
        "intermediate_time_unit", "late_time_value", "late_time_unit", "rescue_method",
        "rescue_direction", "blocking_reagent", "blocking_dose",
        "blocking_dose_unit", "phenotype_id_1", "phenotype_id_2",
        "target_removed_from_axes", "programs_with_target_removed",
        "outcomes_inspected_before_design_freeze",
    ],
    "stage_b_power_units_template.tsv": [
        "target_uid", "background_id", "allele_role", "clone_id",
        "differentiation_id", "power_family_alpha", "target_power",
        "pilot_variance_source", "maximum_feasible_expansion",
        "independent_clone", "independent_differentiation",
    ],
    "stage_b_randomization_template.tsv": [
        "sample_id", "blinded_label", "target_uid", "background_id",
        "allele_role", "clone_id", "differentiation_id", "arm_id",
        "time_role", "challenge_role", "assay_id", "batch_id", "plate_id",
        "well_id", "randomized_allocation", "unblinding_status",
        "preoutcome_exclusion_status",
    ],
    "stage_b_freeze_signoff_template.tsv": [
        "stage_b_design_path", "stage_b_design_sha256", "power_units_path",
        "power_units_sha256", "randomization_path", "randomization_sha256",
        "reviewer_1", "reviewer_2", "review_concordant",
        "stage_b_outcomes_inspected", "signed_utc",
    ],
}


POLICIES = [
    ("ADJ01", "complete_family", "all frozen target-guide estimands are reported; missing or invalid rows cannot nominate a locus"),
    ("ADJ02", "cis_gate", "both guides pass CIS01 at BH q<0.05 in the expected direction with at least two backgrounds"),
    ("ADJ03", "relay_gate", "both guides pass RELAY01 and ORIGIN01 for the same prespecified recipient lineage after CIS01"),
    ("ADJ04", "composition_gate", "RELAY01 direction survives both substate-adjusted and fixed-cell-count sensitivities for both guides"),
    ("ADJ05", "route_completeness", "the frozen transfer-versus-contact ROUTE01 is estimable for the advancing recipient without forcing significance"),
    ("ADJ06", "nomination_cap", "at most two loci advance by deterministic worst-gate q ranking with target UID as tie-breaker"),
    ("ADJ07", "mediator_freeze", "one family-corrected intermediate mediator with source secretion, recipient receptor, feasible blocker, upstream-cis assay, and off-target control freezes per advancing locus"),
    ("ADJ08", "stage_b_design", "exact endogenous edit, target rescue, mediator blockade, two non-RNA phenotypes, three backgrounds, two clones per allele/background, and two differentiations per clone freeze before Stage-B outcomes"),
    ("ADJ09", "one_way_gate", "Stage-B target, mediator, endpoint, phenotype, or replacement-locus substitution after this freeze is prohibited"),
    ("ADJ10", "valid_terminal", "valid null, cell-autonomous boundary, composition-mediated effect, no mediator, or QC-indeterminate outcome terminates escalation without replacement search"),
]


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite Stage-A terminal contract: {CANDIDATE_ROOT}")
    execution_contract_root = candidate_source("PLAN45_STAGE_A_EXECUTION_CONTRACT_ROOT")
    qc_contract_root = candidate_source("PLAN45_STAGE_A_QC_CONTRACT_ROOT")
    execution_seal_path = execution_contract_root / "STAGE_A_EXECUTION_CONTRACT_SEALED.json"
    qc_seal_path = qc_contract_root / "STAGE_A_QC_UNBLINDING_CONTRACT_SEALED.json"
    execution_seal = json.loads(execution_seal_path.read_text(encoding="utf-8"))
    qc_seal = json.loads(qc_seal_path.read_text(encoding="utf-8"))
    if execution_seal.get("status") != "sealed_outcome_blind_stage_a_execution_contract":
        raise RuntimeError("Invalid Stage-A execution-contract dependency")
    if qc_seal.get("status") != "sealed_outcome_blind_stage_a_qc_unblinding_contract":
        raise RuntimeError("Invalid Stage-A QC-contract dependency")
    if execution_seal.get("scientific_outcomes_inspected") is not False or qc_seal.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("Terminal contract dependency reports outcome access")
    for root, seal in [(execution_contract_root, execution_seal), (qc_contract_root, qc_seal)]:
        for name, expected in seal["output_sha256"].items():
            path = root / name
            if not path.is_file() or sha256_file(path) != expected:
                raise RuntimeError(f"Terminal-contract dependency drift: {path}")

    CANDIDATE_ROOT.mkdir(parents=True)
    outputs: list[Path] = []
    for name, fields in SCHEMAS.items():
        path = CANDIDATE_ROOT / name
        write_tsv(path, [], fields)
        outputs.append(path)
    policy_path = CANDIDATE_ROOT / "stage_a_terminal_policy.tsv"
    write_tsv(
        policy_path,
        [{"rule_id": rule, "domain": domain, "frozen_rule": text} for rule, domain, text in POLICIES],
        ["rule_id", "domain", "frozen_rule"],
    )
    outputs.append(policy_path)
    manifest_rows = []
    for role, path in [
        ("terminal_contract_producer", SCRIPT_ROOT / "39_build_stage_a_terminal_contract.py"),
        ("terminal_contract_validator", SCRIPT_ROOT / "40_validate_stage_a_terminal_contract.py"),
        ("terminal_adjudication_producer", SCRIPT_ROOT / "41_adjudicate_stage_a_and_freeze_stage_b.py"),
        ("terminal_adjudication_validator", SCRIPT_ROOT / "42_validate_stage_a_terminal_adjudication.py"),
        ("execution_contract_seal", execution_seal_path),
        ("qc_unblinding_contract_seal", qc_seal_path),
    ]:
        if not path.is_file():
            raise RuntimeError(f"Terminal-contract dependency absent: {path}")
        manifest_rows.append({
            "role": role, "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size, "sha256": sha256_file(path),
        })
    manifest_path = CANDIDATE_ROOT / "stage_a_terminal_contract_input_manifest.tsv"
    write_tsv(manifest_path, manifest_rows, ["role", "source_path", "size_bytes", "sha256"])
    outputs.append(manifest_path)
    payload = {
        "status": "sealed_outcome_blind_stage_a_terminal_and_stage_b_freeze_contract",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "scientific_outcomes_inspected": False,
        "stage_a_terminal_verdict_frozen": False,
        "stage_b_design_frozen": False,
        "maximum_stage_b_nominations": 2,
        "primary_q_threshold": 0.05,
        "mediator_q_threshold": 0.05,
        "minimum_stage_a_backgrounds": 2,
        "minimum_stage_b_backgrounds": 3,
        "minimum_clones_per_allele_background": 2,
        "minimum_differentiations_per_clone": 2,
        "n_policy_rules": len(POLICIES),
        "n_input_schemas": len(SCHEMAS),
        "stage_b_promotion_authorized": False,
        "output_sha256": {path.name: sha256_file(path) for path in outputs},
    }
    atomic_write_json(CANDIDATE_ROOT / "STAGE_A_TERMINAL_CONTRACT_SEALED.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
