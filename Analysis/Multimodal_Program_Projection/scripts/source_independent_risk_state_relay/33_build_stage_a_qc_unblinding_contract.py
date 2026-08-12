#!/usr/bin/env python3
"""Seal the Plan 45 blinded-QC and unblinding-authorization contract."""

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
    "blinded_sample_file_manifest_template.tsv": [
        "sample_id", "blinded_label", "assay_id", "raw_source_path",
        "size_bytes", "sha256", "file_readable", "scientific_condition_inspected",
        "reviewer_1", "reviewer_2", "review_concordant",
    ],
    "blinded_assay_qc_thresholds_template.tsv": [
        "assay_id", "metric_name", "operator", "threshold", "metric_unit",
        "threshold_source", "frozen_utc", "raw_qc_accessed_before_freeze",
        "reviewer_1", "reviewer_2", "review_concordant",
    ],
    "threshold_freeze_signoff_template.tsv": [
        "thresholds_path", "thresholds_sha256", "reviewer_1", "reviewer_2",
        "review_concordant", "raw_qc_accessed_before_freeze",
        "scientific_condition_inspected", "signed_utc",
    ],
    "blinded_sample_qc_metrics_template.tsv": [
        "sample_id", "blinded_label", "assay_id", "metric_name", "metric_value",
        "metric_unit", "metric_source_path", "metric_source_sha256",
        "scientific_condition_inspected", "reviewer_1", "reviewer_2",
        "review_concordant",
    ],
    "target_independent_guide_response_template.tsv": [
        "gene_id", "contrast", "effect", "p_value", "q_value", "assay_id",
        "n_backgrounds", "exclusion_gene", "condition_labels_opened",
        "reviewer_1", "reviewer_2", "review_concordant",
    ],
    "blinded_qc_signoff_template.tsv": [
        "review_scope", "reviewer_1", "reviewer_2", "review_concordant",
        "scientific_condition_inspected", "signed_utc",
    ],
    "unblinding_key_template.tsv": [
        "sample_id", "blinded_label", "biological_unit_id", "background_id",
        "differentiation_id", "target_uid", "guide_id", "arm_id", "time_role",
        "challenge_role", "challenge_id", "assay_id", "key_generated_before_qc",
        "key_opened_after_qc_seal", "reviewer_1", "reviewer_2",
        "review_concordant",
    ],
    "unblinding_signoff_template.tsv": [
        "qc_seal_path", "qc_seal_sha256", "unblinding_key_path",
        "unblinding_key_sha256", "data_manager", "analysis_lead",
        "key_opened_after_qc_seal", "scientific_outcomes_inspected",
        "signed_utc",
    ],
}


POLICIES = [
    ("QC01", "execution_prerequisite", "consume exactly one validated target-bound Stage-A execution freeze"),
    ("QC02", "complete_sample_universe", "every randomized sample receives one blinded file audit and every assay-specific threshold metric"),
    ("QC03", "threshold_precommitment", "assay thresholds are dual-reviewed and cryptographically sealed in a separate release before raw QC values are accessed"),
    ("QC04", "blinded_exclusions", "all technical sample exclusions freeze under opaque labels before the condition key opens"),
    ("QC05", "guide_response", "target-independent NTC-versus-untransduced guide-response genes are excluded at BH q<0.05 and absolute effect at least 0.25 over the full observable family"),
    ("QC06", "condition_firewall", "blinded QC may inspect technical metrics and the target-independent guide response but no risk/control or target condition"),
    ("QC07", "key_bijection", "the data-manager key must reproduce every frozen randomization identity exactly once and cannot add or replace a sample"),
    ("QC08", "temporal_order", "the unblinding key opens only after the blinded-QC seal and exclusion-gene hash are final"),
    ("QC09", "post_qc_pairing", "only prespecified risk/control pairs retained within a biological unit enter analysis; at least two backgrounds are required for a claim-bearing Stage-A comparison"),
    ("QC10", "valid_null", "a valid null or QC-incomplete comparison remains reportable; biological positivity is never required for release"),
]


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite QC/unblinding contract: {CANDIDATE_ROOT}")
    execution_contract_root = candidate_source("PLAN45_STAGE_A_EXECUTION_CONTRACT_ROOT")
    execution_seal_path = execution_contract_root / "STAGE_A_EXECUTION_CONTRACT_SEALED.json"
    execution_seal = json.loads(execution_seal_path.read_text(encoding="utf-8"))
    if execution_seal.get("status") != "sealed_outcome_blind_stage_a_execution_contract":
        raise RuntimeError("Invalid Stage-A execution-contract dependency")
    if execution_seal.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("Execution-contract dependency reports outcome access")
    for name, expected in execution_seal["output_sha256"].items():
        path = execution_contract_root / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Execution-contract source drift: {path}")

    CANDIDATE_ROOT.mkdir(parents=True)
    outputs: list[Path] = []
    for name, fields in SCHEMAS.items():
        path = CANDIDATE_ROOT / name
        write_tsv(path, [], fields)
        outputs.append(path)
    policy_path = CANDIDATE_ROOT / "stage_a_qc_unblinding_policy.tsv"
    write_tsv(
        policy_path,
        [
            {"rule_id": rule_id, "domain": domain, "frozen_rule": rule}
            for rule_id, domain, rule in POLICIES
        ],
        ["rule_id", "domain", "frozen_rule"],
    )
    outputs.append(policy_path)

    manifest_rows = []
    for role, path in [
        ("qc_contract_producer", SCRIPT_ROOT / "33_build_stage_a_qc_unblinding_contract.py"),
        ("qc_contract_validator", SCRIPT_ROOT / "34_validate_stage_a_qc_unblinding_contract.py"),
        ("qc_threshold_producer", SCRIPT_ROOT / "35a_freeze_stage_a_qc_thresholds.py"),
        ("qc_threshold_validator", SCRIPT_ROOT / "35b_validate_stage_a_qc_thresholds.py"),
        ("blinded_qc_producer", SCRIPT_ROOT / "35_seal_blinded_stage_a_qc.py"),
        ("blinded_qc_validator", SCRIPT_ROOT / "36_validate_blinded_stage_a_qc.py"),
        ("unblinding_producer", SCRIPT_ROOT / "37_authorize_stage_a_unblinding.py"),
        ("unblinding_validator", SCRIPT_ROOT / "38_validate_stage_a_unblinding.py"),
        ("execution_contract_seal", execution_seal_path),
    ]:
        if not path.is_file():
            raise RuntimeError(f"QC/unblinding dependency absent: {path}")
        manifest_rows.append({
            "role": role,
            "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    manifest_path = CANDIDATE_ROOT / "stage_a_qc_unblinding_contract_input_manifest.tsv"
    write_tsv(manifest_path, manifest_rows, ["role", "source_path", "size_bytes", "sha256"])
    outputs.append(manifest_path)

    payload = {
        "status": "sealed_outcome_blind_stage_a_qc_unblinding_contract",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "scientific_outcomes_inspected": False,
        "scientific_condition_labels_opened": False,
        "sample_exclusions_frozen": False,
        "unblinding_authorized": False,
        "stage_b_design_frozen": False,
        "n_policy_rules": len(POLICIES),
        "n_input_schemas": len(SCHEMAS),
        "guide_response_q_threshold": 0.05,
        "guide_response_abs_effect_threshold": 0.25,
        "minimum_backgrounds_per_claim_bearing_comparison": 2,
        "governing_plan_path": "docs/archive/plans/2026-08-07_paper_program/45_SOURCE_INDEPENDENT_RISK_TO_STATE_RELAY.md",
        "governing_plan_integrity_bound": False,
        "output_sha256": {path.name: sha256_file(path) for path in outputs},
    }
    atomic_write_json(CANDIDATE_ROOT / "STAGE_A_QC_UNBLINDING_CONTRACT_SEALED.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
