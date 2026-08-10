#!/usr/bin/env python3
"""Independently validate the outcome-blind Plan 45 target contract."""

from __future__ import annotations

import json

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


def main() -> None:
    seal_path = CANDIDATE_ROOT / "TARGET_ADJUDICATION_CONTRACT_SEALED.json"
    if not seal_path.is_file():
        raise RuntimeError("Missing target-adjudication contract seal")
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "sealed_outcome_blind_target_adjudication_contract":
        raise RuntimeError("Invalid target-adjudication contract status")
    if seal.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("Contract claims scientific outcome access")
    if seal.get("experimental_targets_frozen") is not False:
        raise RuntimeError("Schema contract improperly froze a target")
    if seal.get("governing_plan_integrity_bound") is not False:
        raise RuntimeError("Mutable governing plan must not be self-referentially hash-bound")
    if seal.get("screen_target_min") != 4 or seal.get("screen_target_max") != 6:
        raise RuntimeError("Stage-A screen-size rule drift")
    if seal.get("required_stage_a_guides") != 2:
        raise RuntimeError("Stage-A guide-count rule drift")
    if seal.get("required_stage_b_pegrnas_per_exact_variant") != 2:
        raise RuntimeError("Stage-B pegRNA-count rule drift")

    for name, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Target-contract output hash mismatch: {name}")

    policies = read_tsv(CANDIDATE_ROOT / "target_selection_policy.tsv")
    if [row["rule_id"] for row in policies] != [f"TG{i:02d}" for i in range(1, 9)]:
        raise RuntimeError("Target-selection policy universe drift")
    if any(read_tsv(CANDIDATE_ROOT / name) for name in [
        "stage_a_crispick_audit_template.tsv",
        "stage_b_primedesign_audit_template.tsv",
        "target_independent_protocol_gate_template.tsv",
    ]):
        raise RuntimeError("An outcome-blind schema template contains data rows")

    manifests = read_tsv(CANDIDATE_ROOT / "target_adjudication_input_manifest.tsv")
    if {row["role"] for row in manifests} != {
        "guide_worklist_producer",
        "guide_worklist_validator",
        "stage_a_template_seal",
    }:
        raise RuntimeError("Target-contract manifest universe drift")
    for row in manifests:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Target-contract source drift: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Target-contract source hash drift: {path}")

    print(
        "TARGET_ADJUDICATION_CONTRACT_VALIDATION_PASS "
        "rules=8 targets_frozen=false outcomes_opened=false"
    )


if __name__ == "__main__":
    main()
