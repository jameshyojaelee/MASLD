#!/usr/bin/env python3
"""Seal the outcome-blind contract for importing guides and freezing Stage A.

This contract is created before any external guide result, technical pilot, or
experimental disease-state outcome is opened.  It separates Stage-A promoter
CRISPRa/i designs from Stage-B exact regulatory-allele designs and fixes the
target-selection algorithm.  It contains schemas only and cannot name a locus.
"""

from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_json,
    sha256_file,
    write_tsv,
)


STAGE_A_TEMPLATE_ROOT = (
    PROJECT_ROOT
    / "Analysis/Multimodal_Program_Projection/candidates/"
    "source-independent-risk-state-relay-stage-a-template-v2-2026-08-10"
)
PLAN = (
    PROJECT_ROOT
    / "docs/archive/plans/2026-08-07_paper_program/"
    "45_SOURCE_INDEPENDENT_RISK_TO_STATE_RELAY.md"
)
SCRIPT_ROOT = Path(__file__).resolve().parent

STAGE_A_AUDIT_FIELDS = [
    "promoter_request_uid",
    "orientation_uid",
    "audit_row_type",
    "design_status",
    "guide_id",
    "guide_sequence",
    "perturbation_mode",
    "source_tool",
    "source_tool_version",
    "source_design_rank",
    "source_recommended",
    "sequence_qc_pass",
    "off_target_review_pass",
    "reviewer_1",
    "reviewer_2",
    "review_concordant",
    "raw_export_path",
    "raw_export_sha256",
]
STAGE_B_AUDIT_FIELDS = [
    "design_uid",
    "orientation_uid",
    "audit_row_type",
    "design_status",
    "pegrna_id",
    "spacer_sequence",
    "pbs_sequence",
    "rtt_sequence",
    "ngrna_sequence",
    "source_tool",
    "source_tool_version",
    "source_design_rank",
    "source_recommended",
    "complete_design",
    "intended_edit_matches_worklist",
    "off_target_review_pass",
    "reviewer_1",
    "reviewer_2",
    "review_concordant",
    "raw_export_path",
    "raw_export_sha256",
]
PROTOCOL_FIELDS = [
    "platform_id",
    "pilot_blinded",
    "scientific_outcomes_opened",
    "n_independent_backgrounds",
    "n_independent_differentiations_per_background_condition",
    "viability_gate_pass",
    "maturation_gate_pass",
    "source_lineage_gate_pass",
    "lineage_mix_gate_pass",
    "chronic_challenge_gate_pass",
    "early_intermediate_late_timepoints_locked",
    "randomization_locked",
    "blinding_locked",
    "missingness_locked",
    "power_inputs_locked",
    "protocol_path",
    "protocol_sha256",
    "signed_by_1",
    "signed_by_2",
]


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(
            f"Refusing to overwrite target-adjudication contract: {CANDIDATE_ROOT}"
        )
    guide_script = SCRIPT_ROOT / "14_prepare_guide_design_worklist.py"
    guide_validator = SCRIPT_ROOT / "15_validate_guide_design_worklist.py"
    template_seal = STAGE_A_TEMPLATE_ROOT / "STAGE_A_TEMPLATE_SEALED.json"
    # PLAN is the governing rationale but is intentionally not integrity-bound:
    # it must receive append-only execution status and seal hashes. Hashing it
    # here would create a self-invalidating contract. The policy table emitted
    # below is the frozen experimental-selection specification.
    required = [guide_script, guide_validator, template_seal]
    for path in required:
        if not path.is_file() or path.stat().st_size == 0:
            raise RuntimeError(f"Missing target-adjudication dependency: {path}")
    template = json.loads(template_seal.read_text(encoding="utf-8"))
    if template.get("status") != "sealed_target_independent_stage_a_template":
        raise RuntimeError("Stage-A source is not the authoritative sealed template")
    if template.get("target_freeze_status") != "prohibited":
        raise RuntimeError("Stage-A template unexpectedly opened target selection")

    CANDIDATE_ROOT.mkdir(parents=True)
    stage_a_path = CANDIDATE_ROOT / "stage_a_crispick_audit_template.tsv"
    stage_b_path = CANDIDATE_ROOT / "stage_b_primedesign_audit_template.tsv"
    protocol_path = CANDIDATE_ROOT / "target_independent_protocol_gate_template.tsv"
    write_tsv(stage_a_path, [], STAGE_A_AUDIT_FIELDS)
    write_tsv(stage_b_path, [], STAGE_B_AUDIT_FIELDS)
    write_tsv(protocol_path, [], PROTOCOL_FIELDS)

    policy_rows = [
        {
            "rule_id": "TG01",
            "scope": "stage_a_guideability",
            "rule": (
                "at least two nonidentical CRISPick-recommended promoter guides; "
                "CRISPRa for risk-increases-expression and CRISPRi for "
                "risk-decreases-expression; sequence QC, off-target review, and "
                "dual-review concordance must pass"
            ),
        },
        {
            "rule_id": "TG02",
            "scope": "stage_b_exact_editability",
            "rule": (
                "one accessibility-passing shared-signal variant must have at least "
                "two nonidentical complete PrimeDesign pegRNAs that install the exact "
                "worklist allele and pass off-target and dual-review gates"
            ),
        },
        {
            "rule_id": "TG03",
            "scope": "source_lineage",
            "rule": (
                "primary architecture requires hepatocyte source plus independent "
                "two-cohort hepatocyte accessibility; a nonhepatocyte target cannot "
                "freeze from the current single-cohort accessibility reference"
            ),
        },
        {
            "rule_id": "TG04",
            "scope": "cross_study_orientation",
            "rule": (
                "all eligible direct-trait orientation records for the same physical "
                "locus and target gene must agree; multiple records count as support, "
                "not automatically as statistically independent replication"
            ),
        },
        {
            "rule_id": "TG05",
            "scope": "screen_or_single_locus",
            "rule": (
                "freeze four to six Stage-A targets only when at least four pass and "
                "the eligible universe contains both risk-expression directions and "
                "both direct phenotype strata; otherwise freeze only the top eligible "
                "target as a single-locus mechanism study"
            ),
        },
        {
            "rule_id": "TG06",
            "scope": "deterministic_priority",
            "rule": (
                "maximize uncovered direction/phenotype categories, then order by "
                "number of concordant participant-overlap-collapsed GWAS evidence "
                "families (excluding overlap-unresolved meta-analysis rows from "
                "breadth), orientation consensus, SuSiE PP.H4, exact variant shared "
                "posterior, physical locus ID, and Ensembl ID"
            ),
        },
        {
            "rule_id": "TG07",
            "scope": "protocol",
            "rule": (
                "a target-independent blinded platform pilot must pass viability, "
                "maturation, lineage, challenge, timing, randomization, blinding, "
                "missingness, and power-input gates in at least three backgrounds and "
                "two differentiations per background/condition"
            ),
        },
        {
            "rule_id": "TG08",
            "scope": "outcome_firewall",
            "rule": (
                "guide and protocol audits may contain design/QC facts only; cis, "
                "recipient-axis, mediator, phenotype, and disease-alignment outcomes "
                "are prohibited before the target registry and Stage-A design seal"
            ),
        },
    ]
    policy_path = CANDIDATE_ROOT / "target_selection_policy.tsv"
    write_tsv(policy_path, policy_rows, ["rule_id", "scope", "rule"])

    manifest_rows = []
    for role, path in [
        ("guide_worklist_producer", guide_script),
        ("guide_worklist_validator", guide_validator),
        ("stage_a_template_seal", template_seal),
    ]:
        manifest_rows.append(
            {
                "role": role,
                "source_path": str(path.relative_to(PROJECT_ROOT)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    manifest_path = CANDIDATE_ROOT / "target_adjudication_input_manifest.tsv"
    write_tsv(
        manifest_path,
        manifest_rows,
        ["role", "source_path", "size_bytes", "sha256"],
    )

    payload = {
        "status": "sealed_outcome_blind_target_adjudication_contract",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "candidate_id": CANDIDATE_ROOT.name,
        "governing_plan_path": str(PLAN.relative_to(PROJECT_ROOT)),
        "governing_plan_integrity_bound": False,
        "governing_plan_integrity_note": (
            "append-only status updates would self-invalidate the contract; the "
            "emitted target_selection_policy.tsv is the frozen selection authority"
        ),
        "scientific_outcomes_inspected": False,
        "experimental_targets_frozen": False,
        "stage_a_action": (
            "promoter CRISPRa for risk-increases-expression; promoter CRISPRi for "
            "risk-decreases-expression"
        ),
        "stage_b_action": "exact endogenous accessible shared-signal allele edit",
        "primary_source_lineage": "Hepatocytes",
        "primary_recipient_lineages": template["primary_recipient_lineages"],
        "screen_target_min": 4,
        "screen_target_max": 6,
        "single_locus_fallback_count": 1,
        "required_stage_a_guides": 2,
        "required_stage_b_pegrnas_per_exact_variant": 2,
        "output_sha256": {
            stage_a_path.name: sha256_file(stage_a_path),
            stage_b_path.name: sha256_file(stage_b_path),
            protocol_path.name: sha256_file(protocol_path),
            policy_path.name: sha256_file(policy_path),
            manifest_path.name: sha256_file(manifest_path),
        },
    }
    atomic_write_json(
        CANDIDATE_ROOT / "TARGET_ADJUDICATION_CONTRACT_SEALED.json", payload
    )
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
