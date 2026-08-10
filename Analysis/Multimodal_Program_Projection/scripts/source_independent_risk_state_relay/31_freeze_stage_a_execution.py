#!/usr/bin/env python3
"""Freeze the actual, target-bound Stage-A execution before outcome access."""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
import os
import re
from collections import Counter, defaultdict
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT,
    PROJECT_ROOT,
    atomic_write_json,
    read_tsv,
    sha256_file,
    write_tsv,
)


INPUT_FILES = [
    "pilot_power_audit.tsv",
    "experimental_unit_manifest.tsv",
    "assay_registry.tsv",
    "phenotype_registry.tsv",
    "randomization_manifest.tsv",
    "execution_source_manifest.tsv",
    "dual_review_signoff.tsv",
]
PRIMARY_RECIPIENTS = {
    "Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Macrophages"
}
BASE_ACTIVE_ARMS = {
    "SRC_CTRL", "SRC_RISK", "MOS_CTRL", "MOS_RISK", "RECIP_CTRL", "RECIP_RISK"
}
ALLOWED_MODELS = {
    "target_expression": {"linear_mixed_model"},
    "rna_counts": {"edgeR_QL_pseudobulk", "DESeq2_pseudobulk", "limma_voom_pseudobulk"},
    "secreted_protein": {"limma_mixed_model", "linear_mixed_model"},
    "viability": {"binomial_glmm", "beta_mixed_model", "linear_mixed_model"},
    "lineage_composition": {"logratio_linear_mixed_model", "dirichlet_multinomial"},
    "non_rna_phenotype": {
        "linear_mixed_model", "beta_mixed_model", "negative_binomial_glmm",
        "ordinal_mixed_model",
    },
}


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def finite_number(value: object, field: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"Non-finite {field}: {value!r}")
    return result


def split_values(value: str) -> set[str]:
    return {item.strip() for item in value.split(";") if item.strip()}


def hours(value: object, unit: str) -> float:
    result = finite_number(value, "time")
    if unit == "hours":
        return result
    if unit == "days":
        return result * 24.0
    raise RuntimeError(f"Unsupported time unit: {unit}")


def validate_chronic_schedule(time_hours: dict[str, float], duration_hours: float) -> None:
    if not 24 <= time_hours["early_cis"] <= 48:
        raise RuntimeError("Early cis readout must be 24-48 hours")
    if not 72 <= time_hours["intermediate_mediator"] <= 192:
        raise RuntimeError("Intermediate mediator readout must be 3-8 days")
    if not 240 <= time_hours["late_relay"] <= 504:
        raise RuntimeError("Late relay readout must be 10-21 days")
    if duration_hours < time_hours["late_relay"]:
        raise RuntimeError("Chronic challenge duration ends before the late relay readout")


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


def header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def validate_seal_outputs(root: Path, seal: dict[str, object], suffix: str = "") -> None:
    for name, expected in seal["output_sha256"].items():  # type: ignore[union-attr]
        path = root / f"{name}{suffix}"
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Sealed source drift: {path}")


def exact_schema_rows(input_root: Path, contract_root: Path, basename: str) -> list[dict[str, str]]:
    path = input_root / basename
    template = contract_root / basename.replace(".tsv", "_template.tsv")
    if not path.is_file() or not template.is_file():
        raise RuntimeError(f"Stage-A execution input/schema absent: {path} / {template}")
    observed = header(path)
    expected = header(template)
    if observed != expected:
        raise RuntimeError(f"Stage-A execution schema drift: {basename}")
    rows = read_tsv(path)
    if not rows:
        raise RuntimeError(f"Stage-A execution input has no rows: {basename}")
    return rows


def transfer_arms(route: str) -> tuple[str, str]:
    if route == "transwell":
        return "TRW_CTRL", "TRW_RISK"
    if route == "conditioned_medium":
        return "CM_CTRL", "CM_RISK"
    raise RuntimeError("transfer_route must be transwell or conditioned_medium")


def validate_blocked_pairing(rows: list[dict[str, str]], route: str) -> None:
    transfer_ctrl, transfer_risk = transfer_arms(route)
    pairs = {
        "source": ("SRC_CTRL", "SRC_RISK"),
        "mosaic": ("MOS_CTRL", "MOS_RISK"),
        "reciprocal": ("RECIP_CTRL", "RECIP_RISK"),
        "transfer": (transfer_ctrl, transfer_risk),
    }
    by_arm: dict[tuple[str, ...], dict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )
    arm_to_pair = {
        arm: pair_id for pair_id, pair in pairs.items() for arm in pair
    }
    for row in rows:
        if row["arm_id"] not in arm_to_pair:
            continue
        key = (
            row["target_uid"], row["guide_id"], row["biological_unit_id"],
            row["assay_id"], row["time_role"], row["challenge_role"],
            arm_to_pair[row["arm_id"]],
        )
        by_arm[key][row["arm_id"]].add(row["batch_id"])
    for key, batches in by_arm.items():
        pair = pairs[key[-1]]
        if set(batches) != set(pair) or not (batches[pair[0]] & batches[pair[1]]):
            raise RuntimeError(
                "Risk/control arms are not paired within biological unit and batch: "
                + "|".join(key)
            )


def resolve_contract_values(value: str, active_arms: set[str], active_times: set[str], route: str) -> set[str]:
    transfer_ctrl, transfer_risk = transfer_arms(route)
    if value == "{ACTIVE_ARMS}":
        return set(active_arms)
    if value == "{ACTIVE_TIMES}":
        return set(active_times)
    return {
        item.replace("{TRANSFER_CTRL}", transfer_ctrl).replace("{TRANSFER_RISK}", transfer_risk)
        for item in split_values(value)
    }


def expected_design_keys(
    targets_guides: set[tuple[str, str]],
    units: set[tuple[str, str, str]],
    assays: list[dict[str, str]],
    allowed_arm_times: dict[str, set[str]],
) -> set[tuple[str, ...]]:
    """Expand the exact inferential design; technical replicates never add units."""
    expected: set[tuple[str, ...]] = set()
    for assay in assays:
        arms = split_values(assay["required_arm_ids"])
        times = split_values(assay["required_time_roles"])
        challenges = split_values(assay["required_challenge_roles"])
        for target_uid, guide_id in targets_guides:
            for biological_unit_id, background_id, differentiation_id in units:
                for arm_id in arms:
                    for time_role in times & allowed_arm_times.get(arm_id, set()):
                        for challenge_role in challenges:
                            expected.add(
                                (
                                    target_uid, guide_id, biological_unit_id,
                                    background_id, differentiation_id, arm_id,
                                    time_role, challenge_role, assay["assay_id"],
                                )
                            )
    return expected


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite Stage-A execution freeze: {CANDIDATE_ROOT}")
    target_root = candidate_source("PLAN45_TARGET_FREEZE_ROOT")
    contract_root = candidate_source("PLAN45_STAGE_A_EXECUTION_CONTRACT_ROOT")
    input_root = candidate_source("PLAN45_STAGE_A_EXECUTION_INPUT_ROOT")

    target_seal_path = target_root / "FROZEN_STAGE_A_TARGETS.json"
    contract_seal_path = contract_root / "STAGE_A_EXECUTION_CONTRACT_SEALED.json"
    target_seal = json.loads(target_seal_path.read_text(encoding="utf-8"))
    contract_seal = json.loads(contract_seal_path.read_text(encoding="utf-8"))
    if target_seal.get("status") not in {
        "stage_a_screen_targets_frozen",
        "single_locus_stage_a_target_frozen_architecture_generalization_prohibited",
    }:
        raise RuntimeError("Stage-A execution requires at least one frozen target")
    if not target_seal.get("experimental_targets_frozen"):
        raise RuntimeError("Stage-A target seal freezes no target")
    if target_seal.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("Target-freeze dependency reports outcome access")
    if target_seal.get("stage_b_design_frozen") is not False:
        raise RuntimeError("Target-freeze dependency improperly froze Stage B")
    if contract_seal.get("status") != "sealed_outcome_blind_stage_a_execution_contract":
        raise RuntimeError("Invalid Stage-A execution-contract dependency")
    if contract_seal.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("Execution contract reports outcome access")
    validate_seal_outputs(target_root, target_seal)
    validate_seal_outputs(contract_root, contract_seal)

    input_rows = {
        name: exact_schema_rows(input_root, contract_root, name)
        for name in INPUT_FILES
    }
    pilot_rows = input_rows["pilot_power_audit.tsv"]
    if len(pilot_rows) != 1:
        raise RuntimeError("Pilot/power audit must contain exactly one row")
    pilot = pilot_rows[0]
    if not yes(pilot["blinded_pilot"]) or yes(pilot["scientific_outcomes_opened"]):
        raise RuntimeError("Pilot is not target-independent and outcome-blinded")
    if not yes(pilot["review_concordant"]) or not pilot["reviewer_1"] or not pilot["reviewer_2"]:
        raise RuntimeError("Pilot lacks dual concordant review")
    if pilot["reviewer_1"] == pilot["reviewer_2"]:
        raise RuntimeError("Pilot dual review uses one identity twice")
    time_hours = {
        role: hours(pilot[f"{prefix}_time_value"], pilot[f"{prefix}_time_unit"])
        for role, prefix in [
            ("early_cis", "early"),
            ("intermediate_mediator", "intermediate"),
            ("late_relay", "late"),
        ]
    }
    if not 0 < time_hours["early_cis"] < time_hours["intermediate_mediator"] < time_hours["late_relay"]:
        raise RuntimeError("Pilot times are not strictly ordered")
    if pilot["challenge_exposure_mode"] != "continuous_or_repeated_chronic":
        raise RuntimeError("An acute-pulse challenge cannot satisfy the chronic-state design")
    challenge_duration_hours = hours(
        pilot["challenge_duration_value"], pilot["challenge_duration_unit"]
    )
    validate_chronic_schedule(time_hours, challenge_duration_hours)
    n_backgrounds = int(pilot["stage_a_backgrounds"])
    n_differentiations = int(pilot["stage_a_differentiations_per_background_condition"])
    if n_backgrounds < int(contract_seal["required_min_backgrounds"]):
        raise RuntimeError("Stage-A background count below contract")
    if n_differentiations < int(contract_seal["required_min_differentiations_per_background_condition"]):
        raise RuntimeError("Stage-A differentiation count below contract")
    if int(pilot["max_backgrounds"]) < n_backgrounds or int(pilot["max_differentiations_per_background_condition"]) < n_differentiations:
        raise RuntimeError("Maximum feasible expansion is below the frozen Stage-A design")
    alpha = finite_number(pilot["familywise_alpha"], "familywise_alpha")
    power = finite_number(pilot["target_power"], "target_power")
    if pilot["power_scope"] != "stage_b_planning_input_only":
        raise RuntimeError("Pilot power inputs cannot be relabeled as Stage-A confirmatory power")
    if not 0 < alpha <= 0.05 or not 0.80 <= power < 1:
        raise RuntimeError("Pilot alpha/power contract is invalid")
    for field in [
        "power_variance_estimate", "power_effect_size", "viability_threshold",
        "maturation_threshold", "composition_deviation_threshold",
    ]:
        if finite_number(pilot[field], field) <= 0:
            raise RuntimeError(f"Pilot field must be positive: {field}")
    if finite_number(pilot["challenge_concentration"], "challenge_concentration") <= 0 or not pilot["challenge_unit"] or pilot["challenge_id"] == "basal":
        raise RuntimeError("Chronic challenge identity/concentration is not frozen")

    unit_rows = input_rows["experimental_unit_manifest.tsv"]
    if len({row["biological_unit_id"] for row in unit_rows}) != len(unit_rows):
        raise RuntimeError("Biological-unit identifiers are not unique")
    if any(
        not yes(row["independent_background"])
        or not yes(row["independent_differentiation"])
        or not yes(row["lineage_mix_qc_planned"])
        or row["clone_id"] != "not_applicable_stage_a"
        or row["platform_id"] != pilot["platform_id"]
        or row["source_lineage"] != "Hepatocytes"
        for row in unit_rows
    ):
        raise RuntimeError("Biological-unit independence/source-lineage gate failed")
    backgrounds = {row["background_id"] for row in unit_rows}
    if len(backgrounds) != n_backgrounds:
        raise RuntimeError("Biological-unit background count differs from the pilot freeze")
    by_background = Counter(row["background_id"] for row in unit_rows)
    if set(by_background.values()) != {n_differentiations}:
        raise RuntimeError("Each background must have the exact frozen differentiation count")
    if any(
        len({row["differentiation_id"] for row in unit_rows if row["background_id"] == background})
        != n_differentiations
        for background in backgrounds
    ):
        raise RuntimeError("Differentiation identifiers are not unique within background")

    registry = read_tsv(target_root / "frozen_target_registry.tsv")
    selected = [row for row in registry if row["selection_status"] == "selected_for_stage_a"]
    target_uids = {row["target_uid"] for row in selected}
    if target_uids != set(target_seal["frozen_target_uids"]):
        raise RuntimeError("Target universe differs from the validated target seal")
    target_symbols = {row["target_uid"]: row["gene_symbol"] for row in selected}
    guides = read_tsv(target_root / "frozen_stage_a_guides.tsv")
    targets_guides = {(row["target_uid"], row["guide_id"]) for row in guides}
    if len(targets_guides) != 2 * len(target_uids):
        raise RuntimeError("Every Stage-A target must retain exactly two guides")

    target_design = read_tsv(target_root / "stage_a_design.tsv")
    allowed_arm_times: dict[str, set[str]] = defaultdict(set)
    for row in target_design:
        allowed_arm_times[row["arm_id"]] |= split_values(row["required_time_roles"])
    transfer_ctrl, transfer_risk = transfer_arms(pilot["transfer_route"])
    active_arms = BASE_ACTIVE_ARMS | {transfer_ctrl, transfer_risk}
    if not active_arms <= set(allowed_arm_times):
        raise RuntimeError("Active arm universe is absent from the target-bound design")
    active_times = set().union(*(allowed_arm_times[arm] for arm in active_arms))

    assay_rows = input_rows["assay_registry.tsv"]
    if len({row["assay_id"] for row in assay_rows}) != len(assay_rows):
        raise RuntimeError("Assay identifiers are not unique")
    role_counts = Counter(row["assay_role"] for row in assay_rows)
    if role_counts != Counter({
        "CIS_TARGET": 1, "RECIPIENT_RNA": 1, "SECRETOME_DISCOVERY": 1,
        "VIABILITY": 1, "COMPOSITION": 1, "PHENOTYPE": 2,
    }):
        raise RuntimeError("Actual assay registry must contain five singleton roles and two phenotype assays")
    contract_assays = {
        row["assay_role"]: row
        for row in read_tsv(contract_root / "stage_a_required_assay_roles.tsv")
    }
    for row in assay_rows:
        if yes(row["outcomes_present"]) or not yes(row["frozen_before_unblinding"]):
            raise RuntimeError("Assay registry indicates outcome access or post-unblinding selection")
        if not yes(row["review_concordant"]) or not row["reviewer_1"] or not row["reviewer_2"] or row["reviewer_1"] == row["reviewer_2"]:
            raise RuntimeError("Assay registry lacks independent dual review")
        if row["modality"] not in ALLOWED_MODELS or row["model_family"] not in ALLOWED_MODELS[row["modality"]]:
            raise RuntimeError(f"Assay-native model mismatch: {row['assay_id']}")
        if row["biological_unit_definition"] != "background_x_independent_differentiation":
            raise RuntimeError("Assay registry changes the inferential unit")
        if row["technical_unit"] in {"biological_unit", "background", "differentiation"}:
            raise RuntimeError("Assay technical unit is mislabeled as biological")
        required = contract_assays[row["assay_role"]]
        if row["modality"] != required["minimum_modality_class"]:
            raise RuntimeError(f"Assay role is paired with the wrong modality: {row['assay_id']}")
        expected_primary_status = (
            "stage_a_discovery" if row["assay_role"] == "SECRETOME_DISCOVERY"
            else "primary_gate"
        )
        if row["primary_status"] != expected_primary_status:
            raise RuntimeError(f"Assay inferential role drift: {row['assay_id']}")
        if not row["normalization"] or not row["outcome_unit"] or not row["effect_unit"] or not row["outcome_file_path"]:
            raise RuntimeError(f"Assay registry has an incomplete outcome contract: {row['assay_id']}")
        minimum_arms = resolve_contract_values(required["required_arm_ids"], active_arms, active_times, pilot["transfer_route"])
        minimum_times = resolve_contract_values(required["required_time_roles"], active_arms, active_times, pilot["transfer_route"])
        minimum_challenges = split_values(required["required_challenge_roles"])
        declared_arms = split_values(row["required_arm_ids"])
        declared_times = split_values(row["required_time_roles"])
        declared_challenges = split_values(row["required_challenge_roles"])
        if not minimum_arms <= declared_arms or not minimum_times <= declared_times or not minimum_challenges <= declared_challenges:
            raise RuntimeError(f"Assay coverage below the frozen minimum: {row['assay_id']}")
        if not declared_arms <= active_arms or not declared_times <= active_times or declared_challenges != {"basal", "primary_chronic"}:
            raise RuntimeError(f"Assay declares an unauthorized design cell: {row['assay_id']}")

    phenotype_rows = input_rows["phenotype_registry.tsv"]
    if len(phenotype_rows) != int(contract_seal["required_primary_phenotypes"]):
        raise RuntimeError("Exactly two non-RNA phenotypes must be frozen")
    if len({row["phenotype_id"] for row in phenotype_rows}) != 2 or len({row["phenotype_class"] for row in phenotype_rows}) != 2:
        raise RuntimeError("Phenotypes must have unique IDs and distinct biological classes")
    phenotype_assays = {row["assay_id"] for row in assay_rows if row["assay_role"] == "PHENOTYPE"}
    if {row["assay_id"] for row in phenotype_rows} != phenotype_assays:
        raise RuntimeError("Phenotype registry does not cover exactly the two phenotype assays")
    for row in phenotype_rows:
        if row["expected_direction"] not in {"increase", "decrease"}:
            raise RuntimeError("Phenotype direction must be frozen")
        if row["multiplicity_family"] != "PHENO01_all_targets_two_assays":
            raise RuntimeError("Phenotype multiplicity family drift")
        if not yes(row["selected_before_outcomes"]) or not yes(row["source_supported"]) or not yes(row["review_concordant"]):
            raise RuntimeError("Phenotype was not source-supported and frozen before outcomes")
        if not row["reviewer_1"] or not row["reviewer_2"] or row["reviewer_1"] == row["reviewer_2"]:
            raise RuntimeError("Phenotype dual review failed")
        assay = next(value for value in assay_rows if value["assay_id"] == row["assay_id"])
        if row["model_family"] != assay["model_family"] or row["primary_effect_unit"] != assay["effect_unit"]:
            raise RuntimeError("Phenotype model/effect unit differs from its assay registry")

    expected_units = {
        (row["biological_unit_id"], row["background_id"], row["differentiation_id"])
        for row in unit_rows
    }
    expected = expected_design_keys(targets_guides, expected_units, assay_rows, allowed_arm_times)
    randomization = input_rows["randomization_manifest.tsv"]
    if len({row["sample_id"] for row in randomization}) != len(randomization):
        raise RuntimeError("Randomization sample IDs are not unique")
    if len({row["blinded_label"] for row in randomization}) != len(randomization):
        raise RuntimeError("Blinded labels are not unique")
    observed: set[tuple[str, ...]] = set()
    unit_lookup = {
        row["biological_unit_id"]: (row["background_id"], row["differentiation_id"])
        for row in unit_rows
    }
    assay_ids = {row["assay_id"] for row in assay_rows}
    time_values = {
        "early_cis": (pilot["early_time_value"], pilot["early_time_unit"]),
        "intermediate_mediator": (pilot["intermediate_time_value"], pilot["intermediate_time_unit"]),
        "late_relay": (pilot["late_time_value"], pilot["late_time_unit"]),
    }
    for row in randomization:
        if (row["target_uid"], row["guide_id"]) not in targets_guides:
            raise RuntimeError("Randomization contains an unfrozen target-guide pair")
        if row["biological_unit_id"] not in unit_lookup or unit_lookup[row["biological_unit_id"]] != (row["background_id"], row["differentiation_id"]):
            raise RuntimeError("Randomization biological-unit mapping drift")
        if row["assay_id"] not in assay_ids or row["arm_id"] not in active_arms:
            raise RuntimeError("Randomization contains an unauthorized assay/arm")
        if row["time_role"] not in allowed_arm_times[row["arm_id"]]:
            raise RuntimeError("Randomization contains an arm/time mismatch")
        if (row["time_value"], row["time_unit"]) != time_values[row["time_role"]]:
            raise RuntimeError("Randomization time differs from the blinded pilot freeze")
        expected_challenge = "basal" if row["challenge_role"] == "basal" else pilot["challenge_id"]
        if row["challenge_role"] not in {"basal", "primary_chronic"} or row["challenge_id"] != expected_challenge:
            raise RuntimeError("Randomization challenge identity drift")
        if not yes(row["randomized_allocation"]) or row["unblinding_status"] != "sealed" or row["exclusion_status"] != "included_preoutcome" or row["exclusion_reason"]:
            raise RuntimeError("Randomization/blinding/pre-outcome exclusion firewall failed")
        if not all(row[field] for field in ["batch_id", "plate_id", "well_id", "technical_replicate_id"]):
            raise RuntimeError("Randomization lacks a technical allocation identifier")
        label = row["blinded_label"]
        forbidden_tokens = {
            row["arm_id"].lower(), row["target_uid"].lower(),
            target_symbols[row["target_uid"]].lower(), "risk", "control",
        }
        if not re.fullmatch(r"BLIND_[A-Z0-9]{8,}", label) or any(token and token in label.lower() for token in forbidden_tokens):
            raise RuntimeError("Analysis-facing label is informative or malformed")
        observed.add(
            (
                row["target_uid"], row["guide_id"], row["biological_unit_id"],
                row["background_id"], row["differentiation_id"], row["arm_id"],
                row["time_role"], row["challenge_role"], row["assay_id"],
            )
        )
    if observed != expected:
        missing = len(expected - observed)
        extra = len(observed - expected)
        raise RuntimeError(f"Randomized design is incomplete or expanded: missing={missing}, extra={extra}")
    validate_blocked_pairing(randomization, pilot["transfer_route"])

    source_rows = input_rows["execution_source_manifest.tsv"]
    if len({row["source_id"] for row in source_rows}) != len(source_rows):
        raise RuntimeError("Execution source IDs are not unique")
    for row in source_rows:
        if yes(row["scientific_outcome"]):
            raise RuntimeError("Execution freeze source manifest contains a scientific outcome")
        source_path = (PROJECT_ROOT / row["source_path"]).resolve()
        if PROJECT_ROOT.resolve() not in source_path.parents or not source_path.is_file():
            raise RuntimeError(f"Execution source absent/outside project: {source_path}")
        if source_path.stat().st_size != int(row["size_bytes"]) or sha256_file(source_path) != row["sha256"]:
            raise RuntimeError(f"Execution source provenance drift: {source_path}")

    signoffs = input_rows["dual_review_signoff.tsv"]
    required_scopes = {"pilot", "biological_units", "assays", "phenotypes", "randomization", "sources"}
    if {row["review_scope"] for row in signoffs} != required_scopes or len(signoffs) != len(required_scopes):
        raise RuntimeError("Dual-review signoff universe drift")
    for row in signoffs:
        if not row["reviewer_1"] or not row["reviewer_2"] or row["reviewer_1"] == row["reviewer_2"] or not yes(row["review_concordant"]) or yes(row["scientific_outcomes_opened"]):
            raise RuntimeError("Execution signoff is not independently reviewed before outcomes")

    input_manifest_rows = []
    for name in INPUT_FILES:
        path = input_root / name
        input_manifest_rows.append({
            "role": name.removesuffix(".tsv"),
            "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })
    for role, path in [
        ("stage_a_target_freeze_seal", target_seal_path),
        ("stage_a_execution_contract_seal", contract_seal_path),
    ]:
        input_manifest_rows.append({
            "role": role,
            "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        })

    completeness_rows = []
    observed_counts = Counter(key[-1] for key in observed)
    expected_counts = Counter(key[-1] for key in expected)
    for assay_id in sorted(assay_ids):
        completeness_rows.append({
            "assay_id": assay_id,
            "n_expected_design_cells": expected_counts[assay_id],
            "n_observed_design_cells": observed_counts[assay_id],
            "design_complete": str(expected_counts[assay_id] == observed_counts[assay_id]).lower(),
        })
    family_rows = [
        {"estimand_id": "CIS01", "n_primary_tests": len(targets_guides), "family_status": "frozen_complete_family"},
        {"estimand_id": "RELAY01", "n_primary_tests": len(targets_guides) * len(PRIMARY_RECIPIENTS), "family_status": "frozen_complete_family_after_CIS01_gate"},
        {"estimand_id": "ORIGIN01", "n_primary_tests": len(targets_guides) * len(PRIMARY_RECIPIENTS), "family_status": "hierarchical_upper_bound_frozen"},
        {"estimand_id": "ROUTE01", "n_primary_tests": len(targets_guides) * len(PRIMARY_RECIPIENTS), "family_status": "hierarchical_upper_bound_frozen"},
        {"estimand_id": "COMP01", "n_primary_tests": len(targets_guides) * len(PRIMARY_RECIPIENTS), "family_status": "same_family_as_RELAY01"},
        {"estimand_id": "PHENO01", "n_primary_tests": len(targets_guides) * 2, "family_status": "frozen_complete_family"},
    ]
    gate_rows = [{
        "execution_freeze_status": "ready_for_blinded_stage_a_outcome_generation",
        "n_targets": len(target_uids),
        "n_guides": len(targets_guides),
        "n_backgrounds": n_backgrounds,
        "n_differentiations_per_background": n_differentiations,
        "n_assays": len(assay_rows),
        "n_primary_phenotypes": len(phenotype_rows),
        "n_randomized_samples": len(randomization),
        "n_expected_design_cells": len(expected),
        "scientific_outcomes_opened": "false",
        "stage_b_frozen": "false",
    }]

    outputs: list[Path] = []
    for name in INPUT_FILES:
        destination = CANDIDATE_ROOT / f"frozen_{name}"
        rows = input_rows[name]
        write_tsv(destination, rows, header(input_root / name))
        outputs.append(destination)
    manifest_path = CANDIDATE_ROOT / "stage_a_execution_input_manifest.tsv"
    write_tsv(manifest_path, input_manifest_rows, ["role", "source_path", "size_bytes", "sha256"])
    outputs.append(manifest_path)
    completeness_path = CANDIDATE_ROOT / "stage_a_design_completeness.tsv"
    write_tsv(completeness_path, completeness_rows, list(completeness_rows[0]))
    outputs.append(completeness_path)
    family_path = CANDIDATE_ROOT / "stage_a_analysis_family_sizes.tsv"
    write_tsv(family_path, family_rows, list(family_rows[0]))
    outputs.append(family_path)
    gate_path = CANDIDATE_ROOT / "stage_a_execution_gate_status.tsv"
    write_tsv(gate_path, gate_rows, list(gate_rows[0]))
    outputs.append(gate_path)

    payload = {
        "status": "stage_a_execution_frozen_ready_for_blinded_outcome_generation",
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "frozen_target_uids": sorted(target_uids),
        "n_targets": len(target_uids),
        "n_guides": len(targets_guides),
        "n_backgrounds": n_backgrounds,
        "n_differentiations_per_background": n_differentiations,
        "n_assays": len(assay_rows),
        "n_primary_phenotypes": len(phenotype_rows),
        "n_randomized_samples": len(randomization),
        "scientific_outcomes_inspected": False,
        "actual_stage_a_execution_frozen": True,
        "stage_b_design_frozen": False,
        "architecture_generalization_permitted": target_seal.get("architecture_generalization_permitted", False),
        "next_gate": "generate blinded Stage-A outcomes, pass QC, then analyze the frozen CIS01-to-PHENO01 hierarchy",
        "output_sha256": {path.name: sha256_file(path) for path in outputs},
    }
    atomic_write_json(CANDIDATE_ROOT / "FROZEN_STAGE_A_EXECUTION.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
