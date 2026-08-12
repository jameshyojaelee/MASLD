#!/usr/bin/env python3
"""Seal the outcome-blind contract for a target-bound Stage-A execution freeze."""

from __future__ import annotations

import datetime as dt
import json
import os
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    SCRIPT_ROOT,
    atomic_write_json,
    read_tsv,
    sha256_file,
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


SCHEMAS: dict[str, list[str]] = {
    "pilot_power_audit_template.tsv": [
        "platform_id", "pilot_id", "blinded_pilot", "scientific_outcomes_opened",
        "transfer_route", "challenge_id", "challenge_concentration", "challenge_unit",
        "challenge_exposure_mode", "challenge_duration_value", "challenge_duration_unit",
        "early_time_value", "early_time_unit", "intermediate_time_value",
        "intermediate_time_unit", "late_time_value", "late_time_unit",
        "stage_a_backgrounds", "stage_a_differentiations_per_background_condition",
        "familywise_alpha", "target_power", "power_scope", "power_variance_estimate",
        "power_effect_size", "max_backgrounds",
        "max_differentiations_per_background_condition", "viability_threshold",
        "maturation_threshold", "composition_deviation_threshold", "reviewer_1",
        "reviewer_2", "review_concordant",
    ],
    "experimental_unit_manifest_template.tsv": [
        "biological_unit_id", "background_id", "differentiation_id", "clone_id",
        "platform_id", "source_lineage", "independent_background",
        "independent_differentiation", "lineage_mix_qc_planned",
    ],
    "assay_registry_template.tsv": [
        "assay_id", "assay_role", "modality", "model_family", "normalization",
        "outcome_unit", "effect_unit", "required_arm_ids", "required_time_roles",
        "required_challenge_roles", "biological_unit_definition", "technical_unit",
        "primary_status", "outcome_file_path", "outcomes_present",
        "frozen_before_unblinding", "reviewer_1", "reviewer_2",
        "review_concordant",
    ],
    "phenotype_registry_template.tsv": [
        "phenotype_id", "assay_id", "phenotype_class", "expected_direction",
        "transform", "model_family", "primary_effect_unit", "multiplicity_family",
        "selected_before_outcomes", "source_supported", "reviewer_1", "reviewer_2",
        "review_concordant",
    ],
    "randomization_manifest_template.tsv": [
        "sample_id", "blinded_label", "biological_unit_id", "background_id",
        "differentiation_id", "batch_id", "plate_id", "well_id", "target_uid",
        "guide_id", "arm_id", "time_role", "time_value", "time_unit",
        "challenge_role", "challenge_id", "assay_id", "technical_replicate_id",
        "randomized_allocation", "unblinding_status", "exclusion_status",
        "exclusion_reason",
    ],
    "execution_source_manifest_template.tsv": [
        "source_id", "source_path", "size_bytes", "sha256", "source_role",
        "scientific_outcome",
    ],
    "dual_review_signoff_template.tsv": [
        "review_scope", "reviewer_1", "reviewer_2", "review_concordant",
        "scientific_outcomes_opened", "signed_utc",
    ],
}


POLICIES = [
    ("EX01", "target_gate", "consume one independently validated Stage-A target freeze; stop when it freezes zero targets"),
    ("EX02", "outcome_firewall", "all scientific outcomes remain unopened until this execution freeze validates"),
    ("EX03", "biological_units", "at least three independent backgrounds and two independent differentiations per background/condition; wells/cells/lanes are technical"),
    ("EX04", "pilot_lock", "blinded target-independent pilot fixes a continuous-or-repeated chronic challenge, 24-48-hour cis time, 3-8-day mediator time, 10-21-day relay time, transfer route, viability/maturation/composition limits, and Stage-B planning power inputs"),
    ("EX05", "factorial_completeness", "every frozen target-guide, biological unit, prespecified arm/time, challenge, and required assay cell must be randomized before outcome access"),
    ("EX06", "blinding", "analysis labels are opaque and exclusions are assigned before unblinding; outcome-based removal is prohibited"),
    ("EX07", "assay_native", "count RNA, protein, imaging, target-expression, composition, and viability assays retain modality-native units and models"),
    ("EX08", "orthogonal_phenotypes", "exactly two non-RNA phenotypes from distinct biological classes are frozen before outcome access"),
    ("EX09", "mediator_discovery", "Stage A may discover a mediator, but its identity, direction, blockade reagent, dose, and Stage-B test freeze before Stage-B unblinding"),
    ("EX10", "multiplicity", "CIS01 and RELAY01 use complete prespecified BH families; both guides must agree and no secondary program can replace the primary axes"),
    ("EX11", "stage_boundary", "this freeze authorizes Stage A only; exact-edit Stage B, rescue, and blockade remain unfrozen"),
    ("EX12", "valid_null", "a valid null experiment must pass software and QC acceptance; biological positivity is never a release criterion"),
]


ASSAY_REQUIREMENTS = [
    {
        "assay_role": "CIS_TARGET",
        "minimum_modality_class": "target_expression",
        "required_arm_ids": "SRC_CTRL;SRC_RISK;MOS_CTRL;MOS_RISK",
        "required_time_roles": "early_cis",
        "required_challenge_roles": "basal;primary_chronic",
        "primary_use": "CIS01",
    },
    {
        "assay_role": "RECIPIENT_RNA",
        "minimum_modality_class": "rna_counts",
        "required_arm_ids": "MOS_CTRL;MOS_RISK;RECIP_CTRL;RECIP_RISK;{TRANSFER_CTRL};{TRANSFER_RISK}",
        "required_time_roles": "early_cis;late_relay",
        "required_challenge_roles": "basal;primary_chronic",
        "primary_use": "RELAY01;ORIGIN01;ROUTE01;COMP01",
    },
    {
        "assay_role": "SECRETOME_DISCOVERY",
        "minimum_modality_class": "secreted_protein",
        "required_arm_ids": "SRC_CTRL;SRC_RISK;MOS_CTRL;MOS_RISK",
        "required_time_roles": "intermediate_mediator",
        "required_challenge_roles": "basal;primary_chronic",
        "primary_use": "Stage-A mediator nomination only",
    },
    {
        "assay_role": "VIABILITY",
        "minimum_modality_class": "viability",
        "required_arm_ids": "{ACTIVE_ARMS}",
        "required_time_roles": "{ACTIVE_TIMES}",
        "required_challenge_roles": "basal;primary_chronic",
        "primary_use": "release gate",
    },
    {
        "assay_role": "COMPOSITION",
        "minimum_modality_class": "lineage_composition",
        "required_arm_ids": "MOS_CTRL;MOS_RISK;RECIP_CTRL;RECIP_RISK",
        "required_time_roles": "early_cis;late_relay",
        "required_challenge_roles": "basal;primary_chronic",
        "primary_use": "COMP01 and release gate",
    },
    {
        "assay_role": "PHENOTYPE",
        "minimum_modality_class": "non_rna_phenotype",
        "required_arm_ids": "MOS_CTRL;MOS_RISK",
        "required_time_roles": "late_relay",
        "required_challenge_roles": "basal;primary_chronic",
        "primary_use": "PHENO01; exactly two registered assays",
    },
]


ESTIMANDS = [
    ("CIS01", 1, "all frozen target x guide tests", "BH", "expected target-expression direction for each guide; q<0.05; direction in at least two backgrounds"),
    ("RELAY01", 2, "CIS-pass target x guide x four recipient axes", "BH", "positive late-minus-early direct-mosaic effect; q<0.05 for both guides; at least two backgrounds"),
    ("ORIGIN01", 3, "CIS-and-RELAY-pass target x guide x recipient", "BH", "source perturbation exceeds reciprocal-recipient perturbation"),
    ("ROUTE01", 4, "complete passing target x guide x recipient", "BH", "classify transfer versus direct-contact route without forcing a direction"),
    ("COMP01", 5, "same complete family as RELAY01", "BH", "RELAY01 direction survives substate adjustment and fixed-cell-count pseudobulk"),
    ("PHENO01", 6, "advancing target x guide x two frozen phenotype assays", "BH", "expected direction in both distinct non-RNA phenotype classes; both guides agree"),
]


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite execution contract: {CANDIDATE_ROOT}")
    template_root = candidate_source("PLAN45_STAGE_A_TEMPLATE_ROOT")
    template_seal_path = template_root / "STAGE_A_TEMPLATE_SEALED.json"
    template_seal = json.loads(template_seal_path.read_text(encoding="utf-8"))
    if template_seal.get("status") != "sealed_target_independent_stage_a_template":
        raise RuntimeError("Invalid Stage-A template dependency")
    if template_seal.get("experimental_outcomes_inspected") is not False:
        raise RuntimeError("Stage-A template reports outcome access")
    for stem, expected in template_seal["output_sha256"].items():
        path = template_root / f"{stem}.tsv"
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Stage-A template source drift: {path}")

    CANDIDATE_ROOT.mkdir(parents=True)
    output_paths: list[Path] = []
    for name, fields in SCHEMAS.items():
        path = CANDIDATE_ROOT / name
        write_tsv(path, [], fields)
        output_paths.append(path)

    policy_path = CANDIDATE_ROOT / "stage_a_execution_policy.tsv"
    write_tsv(
        policy_path,
        [
            {"rule_id": rule_id, "domain": domain, "frozen_rule": rule}
            for rule_id, domain, rule in POLICIES
        ],
        ["rule_id", "domain", "frozen_rule"],
    )
    output_paths.append(policy_path)

    assay_path = CANDIDATE_ROOT / "stage_a_required_assay_roles.tsv"
    write_tsv(assay_path, ASSAY_REQUIREMENTS, list(ASSAY_REQUIREMENTS[0]))
    output_paths.append(assay_path)

    estimand_path = CANDIDATE_ROOT / "stage_a_analysis_estimand_contract.tsv"
    write_tsv(
        estimand_path,
        [
            {
                "estimand_id": estimand_id,
                "hierarchy": hierarchy,
                "multiplicity_family": family,
                "correction": correction,
                "pass_rule": pass_rule,
            }
            for estimand_id, hierarchy, family, correction, pass_rule in ESTIMANDS
        ],
        ["estimand_id", "hierarchy", "multiplicity_family", "correction", "pass_rule"],
    )
    output_paths.append(estimand_path)

    manifest_rows = []
    for role, path in [
        ("stage_a_execution_contract_producer", SCRIPT_ROOT / "29_build_stage_a_execution_contract.py"),
        ("stage_a_execution_contract_validator", SCRIPT_ROOT / "30_validate_stage_a_execution_contract.py"),
        ("stage_a_target_freeze_producer", SCRIPT_ROOT / "27_adjudicate_and_freeze_stage_a_targets.py"),
        ("stage_a_target_freeze_validator", SCRIPT_ROOT / "28_validate_stage_a_target_freeze.py"),
        ("stage_a_execution_freeze_producer", SCRIPT_ROOT / "31_freeze_stage_a_execution.py"),
        ("stage_a_execution_freeze_validator", SCRIPT_ROOT / "32_validate_stage_a_execution_freeze.py"),
        ("stage_a_template_seal", template_seal_path),
    ]:
        if not path.is_file():
            raise RuntimeError(f"Execution-contract dependency absent: {path}")
        manifest_rows.append(
            {
                "role": role,
                "source_path": str(path.relative_to(PROJECT_ROOT)),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    manifest_path = CANDIDATE_ROOT / "stage_a_execution_contract_input_manifest.tsv"
    write_tsv(manifest_path, manifest_rows, ["role", "source_path", "size_bytes", "sha256"])
    output_paths.append(manifest_path)

    payload = {
        "status": "sealed_outcome_blind_stage_a_execution_contract",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "experimental_targets_frozen": False,
        "actual_stage_a_execution_frozen": False,
        "scientific_outcomes_inspected": False,
        "stage_b_design_frozen": False,
        "n_execution_policy_rules": len(POLICIES),
        "n_required_assay_roles": len(ASSAY_REQUIREMENTS),
        "n_estimands": len(ESTIMANDS),
        "required_primary_phenotypes": 2,
        "required_min_backgrounds": 3,
        "required_min_differentiations_per_background_condition": 2,
        "governing_plan_path": "docs/archive/plans/2026-08-07_paper_program/45_SOURCE_INDEPENDENT_RISK_TO_STATE_RELAY.md",
        "governing_plan_integrity_bound": False,
        "output_sha256": {path.name: sha256_file(path) for path in output_paths},
    }
    atomic_write_json(CANDIDATE_ROOT / "STAGE_A_EXECUTION_CONTRACT_SEALED.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
