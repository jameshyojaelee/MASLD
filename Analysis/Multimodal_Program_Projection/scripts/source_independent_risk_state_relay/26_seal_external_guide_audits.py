#!/usr/bin/env python3
"""Validate and seal standardized external guide/protocol audits.

Raw CRISPick/PrimeDesign exports remain source artifacts. Two reviewers convert
them into the exact contract schemas; this script checks full-universe coverage,
binds every row back to its raw export by SHA256, and copies only design/QC facts
into an immutable candidate. Scientific outcomes are structurally impossible
because extra columns are rejected.
"""

from __future__ import annotations

import datetime as dt
import json
import os
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_json,
    read_tsv,
    sha256_file,
    write_tsv,
)


STAGE_A_FIELDS = [
    "promoter_request_uid", "orientation_uid", "audit_row_type", "design_status",
    "guide_id", "guide_sequence", "perturbation_mode", "source_tool",
    "source_tool_version", "source_design_rank", "source_recommended",
    "sequence_qc_pass", "off_target_review_pass", "reviewer_1", "reviewer_2",
    "review_concordant", "raw_export_path", "raw_export_sha256",
]
STAGE_B_FIELDS = [
    "design_uid", "orientation_uid", "audit_row_type", "design_status",
    "pegrna_id", "spacer_sequence", "pbs_sequence", "rtt_sequence",
    "ngrna_sequence", "source_tool", "source_tool_version", "source_design_rank",
    "source_recommended", "complete_design", "intended_edit_matches_worklist",
    "off_target_review_pass", "reviewer_1", "reviewer_2", "review_concordant",
    "raw_export_path", "raw_export_sha256",
]
PROTOCOL_FIELDS = [
    "platform_id", "pilot_blinded", "scientific_outcomes_opened",
    "n_independent_backgrounds", "n_independent_differentiations_per_background_condition",
    "viability_gate_pass", "maturation_gate_pass", "source_lineage_gate_pass",
    "lineage_mix_gate_pass", "chronic_challenge_gate_pass",
    "early_intermediate_late_timepoints_locked", "randomization_locked",
    "blinding_locked", "missingness_locked", "power_inputs_locked",
    "protocol_path", "protocol_sha256", "signed_by_1", "signed_by_2",
]


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


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
        raise RuntimeError(f"{variable} escapes the Plan 45 candidate root: {path}")
    return path


def read_exact(path: Path, fields: list[str]) -> list[dict[str, str]]:
    import csv

    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != fields:
            raise RuntimeError(
                f"Audit schema mismatch for {path.name}: {reader.fieldnames}"
            )
        rows = list(reader)
    if not rows:
        raise RuntimeError(f"Audit table has no rows: {path}")
    return rows


def valid_dna(value: str, minimum: int = 1, maximum: int = 200) -> bool:
    sequence = value.strip().upper()
    return minimum <= len(sequence) <= maximum and set(sequence) <= set("ACGT")


def verify_bound_file(path_value: str, expected_sha: str, role: str) -> Path:
    path = Path(path_value)
    path = path if path.is_absolute() else PROJECT_ROOT / path
    path = path.resolve()
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Missing {role} source file: {path}")
    observed = sha256_file(path)
    if observed != expected_sha:
        raise RuntimeError(f"{role} SHA256 mismatch: {path}")
    return path


def dual_review_pass(row: dict[str, str]) -> bool:
    return (
        bool(row["reviewer_1"].strip())
        and bool(row["reviewer_2"].strip())
        and row["reviewer_1"].strip() != row["reviewer_2"].strip()
        and yes(row["review_concordant"])
    )


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite guide-audit candidate: {CANDIDATE_ROOT}")
    guide_root = candidate_source("PLAN45_GUIDE_ROOT")
    contract_root = candidate_source("PLAN45_ADJUDICATION_CONTRACT_ROOT")
    input_root = candidate_source("PLAN45_EXTERNAL_AUDIT_INPUT_ROOT")

    guide_seal_path = guide_root / "GUIDE_DESIGN_WORKLIST_SEALED.json"
    contract_seal_path = contract_root / "TARGET_ADJUDICATION_CONTRACT_SEALED.json"
    guide_seal = json.loads(guide_seal_path.read_text(encoding="utf-8"))
    contract_seal = json.loads(contract_seal_path.read_text(encoding="utf-8"))
    if guide_seal.get("status") != "separate_stage_a_promoter_and_stage_b_exact_designs_pending":
        raise RuntimeError("Guide source predates the corrected Stage-A/Stage-B separation")
    if contract_seal.get("status") != "sealed_outcome_blind_target_adjudication_contract":
        raise RuntimeError("Invalid target-adjudication contract source")
    if guide_seal.get("experimental_targets_frozen") is not False:
        raise RuntimeError("Guide source already claims frozen targets")

    stage_a_work_path = guide_root / "stage_a_promoter_design_worklist.tsv"
    stage_b_work_path = guide_root / "stage_b_exact_edit_worklist.tsv"
    for path in [stage_a_work_path, stage_b_work_path]:
        expected = guide_seal["output_sha256"].get(path.name)
        if not expected or sha256_file(path) != expected:
            raise RuntimeError(f"Guide worklist hash mismatch: {path}")
    stage_a_work = {row["promoter_request_uid"]: row for row in read_tsv(stage_a_work_path)}
    stage_b_work = {row["design_uid"]: row for row in read_tsv(stage_b_work_path)}

    source_a = input_root / "stage_a_crispick_design_audit.tsv"
    source_b = input_root / "stage_b_primedesign_design_audit.tsv"
    source_protocol = input_root / "target_independent_protocol_gate.tsv"
    stage_a = read_exact(source_a, STAGE_A_FIELDS)
    stage_b = read_exact(source_b, STAGE_B_FIELDS)
    protocol = read_exact(source_protocol, PROTOCOL_FIELDS)

    if {row["promoter_request_uid"] for row in stage_a} != set(stage_a_work):
        raise RuntimeError("Stage-A audit does not cover the complete promoter-request universe")
    if {row["design_uid"] for row in stage_b} != set(stage_b_work):
        raise RuntimeError("Stage-B audit does not cover the complete exact-variant universe")

    bound_files: dict[Path, str] = {}
    seen_stage_a: set[tuple[str, str]] = set()
    for row in stage_a:
        work = stage_a_work[row["promoter_request_uid"]]
        if row["orientation_uid"] != work["orientation_uid"]:
            raise RuntimeError("Stage-A audit/worklist orientation mismatch")
        if row["perturbation_mode"] != work["stage_a_perturbation_mode"]:
            raise RuntimeError("Stage-A audit risk-direction mode mismatch")
        if row["source_tool"] != "CRISPick" or not row["source_tool_version"].strip():
            raise RuntimeError("Stage-A source tool/version is incomplete")
        if not dual_review_pass(row):
            raise RuntimeError("Stage-A audit lacks concordant independent review")
        raw = verify_bound_file(row["raw_export_path"], row["raw_export_sha256"], "CRISPick")
        bound_files[raw] = "CRISPick_raw_export"
        if row["audit_row_type"] == "guide":
            if row["design_status"] != "design_available":
                raise RuntimeError("Stage-A guide row has invalid design status")
            if not row["guide_id"].strip() or not valid_dna(row["guide_sequence"], 18, 30):
                raise RuntimeError("Stage-A guide identity or sequence is invalid")
            try:
                if int(row["source_design_rank"]) < 1:
                    raise ValueError
            except ValueError as error:
                raise RuntimeError("Stage-A source rank must be a positive integer") from error
            key = (row["promoter_request_uid"], row["guide_id"])
            if key in seen_stage_a:
                raise RuntimeError(f"Duplicate Stage-A guide: {key}")
            seen_stage_a.add(key)
        elif row["audit_row_type"] == "no_design":
            if row["design_status"] not in {"no_design_returned", "source_export_failed"}:
                raise RuntimeError("Invalid Stage-A no-design status")
            if row["guide_id"].strip() or row["guide_sequence"].strip():
                raise RuntimeError("Stage-A no-design row contains a guide")
        else:
            raise RuntimeError("Invalid Stage-A audit row type")

    seen_stage_b: set[tuple[str, str]] = set()
    for row in stage_b:
        work = stage_b_work[row["design_uid"]]
        if row["orientation_uid"] != work["orientation_uid"]:
            raise RuntimeError("Stage-B audit/worklist orientation mismatch")
        if row["source_tool"] != "PrimeDesign" or not row["source_tool_version"].strip():
            raise RuntimeError("Stage-B source tool/version is incomplete")
        if not dual_review_pass(row):
            raise RuntimeError("Stage-B audit lacks concordant independent review")
        raw = verify_bound_file(row["raw_export_path"], row["raw_export_sha256"], "PrimeDesign")
        bound_files[raw] = "PrimeDesign_raw_export"
        if row["audit_row_type"] == "design":
            if row["design_status"] != "design_available":
                raise RuntimeError("Stage-B design row has invalid design status")
            sequences = [
                row["spacer_sequence"], row["pbs_sequence"], row["rtt_sequence"],
                row["ngrna_sequence"],
            ]
            if not row["pegrna_id"].strip() or not all(valid_dna(value) for value in sequences):
                raise RuntimeError("Stage-B pegRNA component is missing or invalid")
            try:
                if int(row["source_design_rank"]) < 1:
                    raise ValueError
            except ValueError as error:
                raise RuntimeError("Stage-B source rank must be a positive integer") from error
            key = (row["design_uid"], row["pegrna_id"])
            if key in seen_stage_b:
                raise RuntimeError(f"Duplicate Stage-B pegRNA: {key}")
            seen_stage_b.add(key)
        elif row["audit_row_type"] == "no_design":
            if row["design_status"] not in {"no_design_returned", "source_export_failed"}:
                raise RuntimeError("Invalid Stage-B no-design status")
            if row["pegrna_id"].strip():
                raise RuntimeError("Stage-B no-design row contains a pegRNA")
        else:
            raise RuntimeError("Invalid Stage-B audit row type")

    for row in protocol:
        if not row["platform_id"].strip():
            raise RuntimeError("Protocol audit lacks platform ID")
        if yes(row["scientific_outcomes_opened"]):
            raise RuntimeError("Protocol audit reports scientific outcome access")
        if not row["signed_by_1"].strip() or not row["signed_by_2"].strip():
            raise RuntimeError("Protocol audit lacks two signers")
        if row["signed_by_1"].strip() == row["signed_by_2"].strip():
            raise RuntimeError("Protocol audit requires two distinct signers")
        raw = verify_bound_file(row["protocol_path"], row["protocol_sha256"], "protocol")
        bound_files[raw] = "target_independent_protocol"

    CANDIDATE_ROOT.mkdir(parents=True)
    output_a = CANDIDATE_ROOT / "stage_a_crispick_design_audit.tsv"
    output_b = CANDIDATE_ROOT / "stage_b_primedesign_design_audit.tsv"
    output_protocol = CANDIDATE_ROOT / "target_independent_protocol_gate.tsv"
    write_tsv(output_a, stage_a, STAGE_A_FIELDS)
    write_tsv(output_b, stage_b, STAGE_B_FIELDS)
    write_tsv(output_protocol, protocol, PROTOCOL_FIELDS)

    manifests = []
    for role, path in [
        ("guide_release", guide_seal_path),
        ("target_adjudication_contract", contract_seal_path),
        ("stage_a_standardized_input", source_a),
        ("stage_b_standardized_input", source_b),
        ("protocol_standardized_input", source_protocol),
    ] + [(role, path) for path, role in sorted(bound_files.items(), key=lambda x: str(x[0]))]:
        manifests.append(
            {
                "role": role,
                "source_path": (
                    str(path.relative_to(PROJECT_ROOT))
                    if PROJECT_ROOT in path.parents else str(path)
                ),
                "size_bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            }
        )
    manifest_path = CANDIDATE_ROOT / "external_guide_audit_input_manifest.tsv"
    write_tsv(manifest_path, manifests, ["role", "source_path", "size_bytes", "sha256"])

    payload = {
        "status": "external_guide_and_target_independent_protocol_audits_sealed",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "source_guide_root": str(guide_root.relative_to(PROJECT_ROOT)),
        "source_contract_root": str(contract_root.relative_to(PROJECT_ROOT)),
        "n_stage_a_audit_rows": len(stage_a),
        "n_stage_b_audit_rows": len(stage_b),
        "n_protocol_rows": len(protocol),
        "scientific_outcomes_inspected": False,
        "experimental_targets_frozen": False,
        "output_sha256": {
            output_a.name: sha256_file(output_a),
            output_b.name: sha256_file(output_b),
            output_protocol.name: sha256_file(output_protocol),
            manifest_path.name: sha256_file(manifest_path),
        },
    }
    atomic_write_json(CANDIDATE_ROOT / "EXTERNAL_GUIDE_AUDIT_SEALED.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
