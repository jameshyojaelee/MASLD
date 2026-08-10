#!/usr/bin/env python3
"""Independently rederive the target-bound Stage-A execution freeze."""

from __future__ import annotations

import json
from collections import Counter, defaultdict

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


PRIMARY_RECIPIENTS = 4


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def values(value: str) -> set[str]:
    return {item.strip() for item in value.split(";") if item.strip()}


def to_hours(value: str, unit: str) -> float:
    number = float(value)
    if unit == "hours":
        return number
    if unit == "days":
        return number * 24.0
    raise RuntimeError(f"Unsupported frozen time unit: {unit}")


def main() -> None:
    seal_path = CANDIDATE_ROOT / "FROZEN_STAGE_A_EXECUTION.json"
    seal = json.loads(seal_path.read_text(encoding="utf-8"))
    if seal.get("status") != "stage_a_execution_frozen_ready_for_blinded_outcome_generation":
        raise RuntimeError("Invalid Stage-A execution-freeze status")
    if seal.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("Execution freeze reports scientific outcome access")
    if seal.get("actual_stage_a_execution_frozen") is not True:
        raise RuntimeError("Execution freeze is not marked frozen")
    if seal.get("stage_b_design_frozen") is not False:
        raise RuntimeError("Execution freeze improperly authorizes Stage B")
    for name, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Stage-A execution output hash mismatch: {name}")

    manifest = read_tsv(CANDIDATE_ROOT / "stage_a_execution_input_manifest.tsv")
    roles = {row["role"] for row in manifest}
    required_roles = {
        "pilot_power_audit", "experimental_unit_manifest", "assay_registry",
        "phenotype_registry", "randomization_manifest", "execution_source_manifest",
        "dual_review_signoff", "stage_a_target_freeze_seal",
        "stage_a_execution_contract_seal",
    }
    if roles != required_roles:
        raise RuntimeError("Stage-A execution input-manifest universe drift")
    for row in manifest:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]):
            raise RuntimeError(f"Stage-A execution source drift: {path}")
        if sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Stage-A execution source hash drift: {path}")

    target_manifest_row = next(row for row in manifest if row["role"] == "stage_a_target_freeze_seal")
    target_seal_path = PROJECT_ROOT / target_manifest_row["source_path"]
    target_root = target_seal_path.parent
    target_seal = json.loads(target_seal_path.read_text(encoding="utf-8"))
    if not target_seal.get("experimental_targets_frozen") or target_seal.get("stage_b_design_frozen") is not False:
        raise RuntimeError("Bound target freeze is not Stage-A-only")
    for name, expected in target_seal["output_sha256"].items():
        path = target_root / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Bound target-freeze output drift: {path}")

    contract_manifest_row = next(row for row in manifest if row["role"] == "stage_a_execution_contract_seal")
    contract_seal_path = PROJECT_ROOT / contract_manifest_row["source_path"]
    contract_root = contract_seal_path.parent
    contract_seal = json.loads(contract_seal_path.read_text(encoding="utf-8"))
    if contract_seal.get("status") != "sealed_outcome_blind_stage_a_execution_contract":
        raise RuntimeError("Bound execution contract is invalid")
    for name, expected in contract_seal["output_sha256"].items():
        path = contract_root / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Bound execution-contract output drift: {path}")

    pilots = read_tsv(CANDIDATE_ROOT / "frozen_pilot_power_audit.tsv")
    if len(pilots) != 1:
        raise RuntimeError("Execution freeze must contain one pilot row")
    pilot = pilots[0]
    if not yes(pilot["blinded_pilot"]) or yes(pilot["scientific_outcomes_opened"]):
        raise RuntimeError("Frozen pilot is not blinded/outcome-free")
    if pilot["challenge_exposure_mode"] != "continuous_or_repeated_chronic":
        raise RuntimeError("Frozen challenge is not chronic")
    schedule = {
        "early": to_hours(pilot["early_time_value"], pilot["early_time_unit"]),
        "intermediate": to_hours(
            pilot["intermediate_time_value"], pilot["intermediate_time_unit"]
        ),
        "late": to_hours(pilot["late_time_value"], pilot["late_time_unit"]),
    }
    duration = to_hours(
        pilot["challenge_duration_value"], pilot["challenge_duration_unit"]
    )
    if not (
        24 <= schedule["early"] <= 48
        and 72 <= schedule["intermediate"] <= 192
        and 240 <= schedule["late"] <= 504
        and duration >= schedule["late"]
    ):
        raise RuntimeError("Frozen chronic schedule is outside the preregistered windows")
    if pilot["power_scope"] != "stage_b_planning_input_only":
        raise RuntimeError("Frozen power inputs overclaim Stage-A power")

    units = read_tsv(CANDIDATE_ROOT / "frozen_experimental_unit_manifest.tsv")
    backgrounds = {row["background_id"] for row in units}
    if len(backgrounds) != int(pilot["stage_a_backgrounds"]):
        raise RuntimeError("Execution background count drift")
    by_background = Counter(row["background_id"] for row in units)
    if set(by_background.values()) != {int(pilot["stage_a_differentiations_per_background_condition"])}:
        raise RuntimeError("Execution differentiation count drift")
    if len({row["biological_unit_id"] for row in units}) != len(units):
        raise RuntimeError("Execution biological units are duplicated")

    target_registry = read_tsv(target_root / "frozen_target_registry.tsv")
    selected_targets = {
        row["target_uid"] for row in target_registry
        if row["selection_status"] == "selected_for_stage_a"
    }
    if selected_targets != set(seal["frozen_target_uids"]):
        raise RuntimeError("Execution target universe differs from the target seal")
    guides = read_tsv(target_root / "frozen_stage_a_guides.tsv")
    target_guides = {(row["target_uid"], row["guide_id"]) for row in guides}
    if len(target_guides) != 2 * len(selected_targets):
        raise RuntimeError("Execution target-guide universe drift")

    assays = read_tsv(CANDIDATE_ROOT / "frozen_assay_registry.tsv")
    role_counts = Counter(row["assay_role"] for row in assays)
    if role_counts != Counter({
        "CIS_TARGET": 1, "RECIPIENT_RNA": 1, "SECRETOME_DISCOVERY": 1,
        "VIABILITY": 1, "COMPOSITION": 1, "PHENOTYPE": 2,
    }):
        raise RuntimeError("Execution assay-role universe drift")
    if any(yes(row["outcomes_present"]) or not yes(row["frozen_before_unblinding"]) for row in assays):
        raise RuntimeError("Frozen assay registry indicates outcome access")
    expected_modalities = {
        "CIS_TARGET": "target_expression",
        "RECIPIENT_RNA": "rna_counts",
        "SECRETOME_DISCOVERY": "secreted_protein",
        "VIABILITY": "viability",
        "COMPOSITION": "lineage_composition",
        "PHENOTYPE": "non_rna_phenotype",
    }
    if any(row["modality"] != expected_modalities[row["assay_role"]] for row in assays):
        raise RuntimeError("Frozen assay role/modality mismatch")
    if any(
        row["primary_status"]
        != ("stage_a_discovery" if row["assay_role"] == "SECRETOME_DISCOVERY" else "primary_gate")
        for row in assays
    ):
        raise RuntimeError("Frozen assay inferential role drift")
    phenotypes = read_tsv(CANDIDATE_ROOT / "frozen_phenotype_registry.tsv")
    if len(phenotypes) != 2 or len({row["phenotype_class"] for row in phenotypes}) != 2:
        raise RuntimeError("Execution freeze lacks two distinct phenotype classes")
    if {row["assay_id"] for row in phenotypes} != {
        row["assay_id"] for row in assays if row["assay_role"] == "PHENOTYPE"
    }:
        raise RuntimeError("Phenotype assays do not match the assay registry")

    target_design = read_tsv(target_root / "stage_a_design.tsv")
    arm_times: dict[str, set[str]] = defaultdict(set)
    for row in target_design:
        arm_times[row["arm_id"]] |= values(row["required_time_roles"])
    randomization = read_tsv(CANDIDATE_ROOT / "frozen_randomization_manifest.tsv")
    observed = {
        (
            row["target_uid"], row["guide_id"], row["biological_unit_id"],
            row["background_id"], row["differentiation_id"], row["arm_id"],
            row["time_role"], row["challenge_role"], row["assay_id"],
        )
        for row in randomization
    }
    expected = set()
    unit_keys = {
        (row["biological_unit_id"], row["background_id"], row["differentiation_id"])
        for row in units
    }
    for assay in assays:
        for target_uid, guide_id in target_guides:
            for biological_unit_id, background_id, differentiation_id in unit_keys:
                for arm_id in values(assay["required_arm_ids"]):
                    for time_role in values(assay["required_time_roles"]) & arm_times[arm_id]:
                        for challenge_role in values(assay["required_challenge_roles"]):
                            expected.add(
                                (
                                    target_uid, guide_id, biological_unit_id,
                                    background_id, differentiation_id, arm_id,
                                    time_role, challenge_role, assay["assay_id"],
                                )
                            )
    if expected != observed:
        raise RuntimeError(
            f"Independent design-cell rederivation failed: missing={len(expected-observed)} extra={len(observed-expected)}"
        )
    if len({row["sample_id"] for row in randomization}) != len(randomization):
        raise RuntimeError("Frozen randomized sample IDs are duplicated")
    if any(
        not yes(row["randomized_allocation"])
        or row["unblinding_status"] != "sealed"
        or row["exclusion_status"] != "included_preoutcome"
        for row in randomization
    ):
        raise RuntimeError("Frozen randomization/blinding state drift")
    transfer_pair = (
        ("TRW_CTRL", "TRW_RISK")
        if pilot["transfer_route"] == "transwell"
        else ("CM_CTRL", "CM_RISK")
    )
    pairs = {
        "source": ("SRC_CTRL", "SRC_RISK"),
        "mosaic": ("MOS_CTRL", "MOS_RISK"),
        "reciprocal": ("RECIP_CTRL", "RECIP_RISK"),
        "transfer": transfer_pair,
    }
    arm_pair = {arm: pair_id for pair_id, pair in pairs.items() for arm in pair}
    blocked: dict[tuple[str, ...], dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    for row in randomization:
        pair_id = arm_pair[row["arm_id"]]
        key = (
            row["target_uid"], row["guide_id"], row["biological_unit_id"],
            row["assay_id"], row["time_role"], row["challenge_role"], pair_id,
        )
        blocked[key][row["arm_id"]].add(row["batch_id"])
    for key, batches in blocked.items():
        pair = pairs[key[-1]]
        if set(batches) != set(pair) or not (batches[pair[0]] & batches[pair[1]]):
            raise RuntimeError("Independent within-unit batch-pairing audit failed")

    source_rows = read_tsv(CANDIDATE_ROOT / "frozen_execution_source_manifest.tsv")
    if any(yes(row["scientific_outcome"]) for row in source_rows):
        raise RuntimeError("Frozen execution sources include a scientific outcome")
    signoffs = read_tsv(CANDIDATE_ROOT / "frozen_dual_review_signoff.tsv")
    required_scopes = {"pilot", "biological_units", "assays", "phenotypes", "randomization", "sources"}
    if {row["review_scope"] for row in signoffs} != required_scopes:
        raise RuntimeError("Frozen dual-review scope drift")
    if any(
        row["reviewer_1"] == row["reviewer_2"]
        or not yes(row["review_concordant"])
        or yes(row["scientific_outcomes_opened"])
        for row in signoffs
    ):
        raise RuntimeError("Frozen dual-review signoff is invalid")

    completeness = read_tsv(CANDIDATE_ROOT / "stage_a_design_completeness.tsv")
    expected_counts = Counter(key[-1] for key in expected)
    if {
        row["assay_id"]: int(row["n_expected_design_cells"])
        for row in completeness
    } != dict(expected_counts):
        raise RuntimeError("Design completeness table does not rederive")
    if any(not yes(row["design_complete"]) for row in completeness):
        raise RuntimeError("Execution freeze contains an incomplete assay")

    families = {
        row["estimand_id"]: int(row["n_primary_tests"])
        for row in read_tsv(CANDIDATE_ROOT / "stage_a_analysis_family_sizes.tsv")
    }
    n_guides = len(target_guides)
    expected_families = {
        "CIS01": n_guides,
        "RELAY01": n_guides * PRIMARY_RECIPIENTS,
        "ORIGIN01": n_guides * PRIMARY_RECIPIENTS,
        "ROUTE01": n_guides * PRIMARY_RECIPIENTS,
        "COMP01": n_guides * PRIMARY_RECIPIENTS,
        "PHENO01": n_guides * 2,
    }
    if families != expected_families:
        raise RuntimeError("Stage-A multiplicity-family sizes do not rederive")

    gates = read_tsv(CANDIDATE_ROOT / "stage_a_execution_gate_status.tsv")
    if len(gates) != 1 or gates[0]["execution_freeze_status"] != "ready_for_blinded_stage_a_outcome_generation":
        raise RuntimeError("Execution gate status drift")
    if yes(gates[0]["scientific_outcomes_opened"]) or yes(gates[0]["stage_b_frozen"]):
        raise RuntimeError("Execution gate improperly opens outcomes or Stage B")
    print(
        "STAGE_A_EXECUTION_FREEZE_VALIDATION_PASS "
        f"targets={len(selected_targets)} guides={n_guides} backgrounds={len(backgrounds)} "
        f"design_cells={len(expected)} outcomes_opened=false stage_b_frozen=false"
    )


if __name__ == "__main__":
    main()
