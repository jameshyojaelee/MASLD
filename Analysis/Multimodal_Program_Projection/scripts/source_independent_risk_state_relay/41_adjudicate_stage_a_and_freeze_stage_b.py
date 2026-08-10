#!/usr/bin/env python3
"""Adjudicate the complete Stage-A family and freeze any eligible Stage-B design."""

from __future__ import annotations

import csv
import datetime as dt
import json
import math
import os
import re
from collections import defaultdict
from pathlib import Path

from relay_common import (
    CANDIDATE_ROOT, PROJECT_ROOT, atomic_write_json, read_tsv, sha256_file,
    write_tsv,
)


RECIPIENTS = ("Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Macrophages")
SENSITIVITIES = ("substate_adjusted", "fixed_cell_count")
STAGE_B_ARMS = {
    "protective_reference", "source_risk", "reciprocal_risk", "all_lineage_risk",
    "source_risk_target_rescue", "source_risk_mediator_blockade", "sham", "parental",
}


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def number(value: object, field: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise RuntimeError(f"Non-finite {field}: {value!r}")
    return result


def utc(value: str) -> dt.datetime:
    parsed = dt.datetime.fromisoformat(value.replace("Z", "+00:00"))
    if parsed.tzinfo is None:
        raise RuntimeError(f"Timestamp lacks timezone: {value!r}")
    return parsed.astimezone(dt.timezone.utc)


def bh_adjust(p_values: list[float]) -> list[float]:
    if any(not 0 <= value <= 1 for value in p_values):
        raise RuntimeError("BH input p-value outside [0,1]")
    n = len(p_values)
    order = sorted(range(n), key=lambda index: (p_values[index], index))
    result = [1.0] * n
    running = 1.0
    for position in range(n - 1, -1, -1):
        index = order[position]
        running = min(running, p_values[index] * n / (position + 1))
        result[index] = min(1.0, running)
    return result


def candidate_source(variable: str) -> Path:
    raw = os.environ.get(variable, "").strip()
    if not raw:
        raise RuntimeError(f"{variable} is required")
    path = Path(raw)
    path = path if path.is_absolute() else PROJECT_ROOT / path
    path = path.resolve()
    allowed = (PROJECT_ROOT / "Analysis/Multimodal_Program_Projection/candidates").resolve()
    if allowed not in path.parents:
        raise RuntimeError(f"{variable} escapes the candidate root: {path}")
    return path


def header(path: Path) -> list[str]:
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.reader(handle, delimiter="\t"))


def exact_rows(
    input_root: Path, contract_root: Path, basename: str, allow_empty: bool = False
) -> list[dict[str, str]]:
    path = input_root / basename
    template = contract_root / basename.replace(".tsv", "_template.tsv")
    if not path.is_file() or not template.is_file() or header(path) != header(template):
        raise RuntimeError(f"Stage-A terminal schema/source failure: {basename}")
    rows = read_tsv(path)
    if not rows and not allow_empty:
        raise RuntimeError(f"Stage-A terminal input is empty: {basename}")
    return rows


def validate_seal(root: Path, basename: str, status: str) -> tuple[Path, dict[str, object]]:
    path = root / basename
    seal = json.loads(path.read_text(encoding="utf-8"))
    if seal.get("status") != status:
        raise RuntimeError(f"Invalid dependency status: {path}")
    for name, expected in seal["output_sha256"].items():  # type: ignore[union-attr]
        output = root / str(name)
        if not output.is_file() or sha256_file(output) != expected:
            raise RuntimeError(f"Terminal-adjudication dependency drift: {output}")
    return path, seal


def validate_q_family(rows: list[dict[str, str]], family: str) -> None:
    adjusted = bh_adjust([number(row["p_value"], f"{family} p") for row in rows])
    for row, q_value in zip(rows, adjusted):
        if abs(number(row["q_value"], f"{family} q") - q_value) > 1e-8:
            raise RuntimeError(f"{family} BH q-value does not rederive")
        if int(row["full_family_size"]) != len(rows):
            raise RuntimeError(f"{family} declared family size drift")


def row_pass(row: dict[str, str], q_threshold: float, min_backgrounds: int) -> bool:
    return (
        row["result_status"] in {"valid_estimate", "tested_null"}
        and number(row["q_value"], "primary q") < q_threshold
        and yes(row["expected_direction_agreement"])
        and int(row["n_backgrounds"]) >= min_backgrounds
    )


def main() -> None:
    if CANDIDATE_ROOT.exists():
        raise RuntimeError(f"Refusing to overwrite terminal adjudication: {CANDIDATE_ROOT}")
    execution_root = candidate_source("PLAN45_STAGE_A_EXECUTION_ROOT")
    unblinding_root = candidate_source("PLAN45_STAGE_A_UNBLINDING_ROOT")
    contract_root = candidate_source("PLAN45_STAGE_A_TERMINAL_CONTRACT_ROOT")
    input_root = candidate_source("PLAN45_STAGE_A_TERMINAL_INPUT_ROOT")
    execution_path, execution = validate_seal(
        execution_root, "FROZEN_STAGE_A_EXECUTION.json",
        "stage_a_execution_frozen_ready_for_blinded_outcome_generation",
    )
    unblinding_path, unblinding = validate_seal(
        unblinding_root, "STAGE_A_UNBLINDING_AUTHORIZED.json",
        "stage_a_unblinding_authorized_for_prespecified_post_qc_pairs",
    )
    contract_path, contract = validate_seal(
        contract_root, "STAGE_A_TERMINAL_CONTRACT_SEALED.json",
        "sealed_outcome_blind_stage_a_terminal_and_stage_b_freeze_contract",
    )
    if unblinding.get("scientific_outcomes_inspected") is not False:
        raise RuntimeError("Stage-A unblinding dependency reports prior outcome access")
    unblinding_manifest = read_tsv(unblinding_root / "unblinding_input_manifest.tsv")
    bound_execution = next(row for row in unblinding_manifest if row["role"] == "stage_a_execution_seal")
    if (
        bound_execution["source_path"] != str(execution_path.relative_to(PROJECT_ROOT))
        or bound_execution["sha256"] != sha256_file(execution_path)
    ):
        raise RuntimeError("Unblinding release is not bound to the supplied execution freeze")

    primary = exact_rows(input_root, contract_root, "stage_a_primary_estimand_results.tsv")
    sensitivity = exact_rows(input_root, contract_root, "stage_a_sensitivity_results.tsv")
    mediators = exact_rows(input_root, contract_root, "stage_a_mediator_candidates.tsv", allow_empty=True)
    outcome_signoffs = exact_rows(input_root, contract_root, "stage_a_outcome_signoff.tsv")
    stage_b_design = exact_rows(input_root, contract_root, "stage_b_exact_design.tsv", allow_empty=True)
    stage_b_units = exact_rows(input_root, contract_root, "stage_b_power_units.tsv", allow_empty=True)
    stage_b_randomization = exact_rows(input_root, contract_root, "stage_b_randomization.tsv", allow_empty=True)
    stage_b_signoffs = exact_rows(input_root, contract_root, "stage_b_freeze_signoff.tsv", allow_empty=True)

    if len(outcome_signoffs) != 1:
        raise RuntimeError("Stage-A outcome signoff must contain exactly one row")
    signoff = outcome_signoffs[0]
    file_bindings = [
        ("unblinding_seal_path", "unblinding_seal_sha256", unblinding_path),
        ("primary_results_path", "primary_results_sha256", input_root / "stage_a_primary_estimand_results.tsv"),
        ("sensitivity_results_path", "sensitivity_results_sha256", input_root / "stage_a_sensitivity_results.tsv"),
        ("mediator_results_path", "mediator_results_sha256", input_root / "stage_a_mediator_candidates.tsv"),
    ]
    for path_field, hash_field, path in file_bindings:
        if (
            Path(signoff[path_field]).as_posix() != path.relative_to(PROJECT_ROOT).as_posix()
            or signoff[hash_field] != sha256_file(path)
        ):
            raise RuntimeError(f"Stage-A outcome signoff binding failed: {path_field}")
    if (
        not yes(signoff["complete_outcome_universe"])
        or not signoff["analyst_1"] or not signoff["analyst_2"]
        or signoff["analyst_1"] == signoff["analyst_2"]
        or utc(signoff["signed_utc"]) <= utc(str(unblinding["created_utc"]))
    ):
        raise RuntimeError("Stage-A outcome signoff is incomplete or chronologically invalid")

    randomization = read_tsv(execution_root / "frozen_randomization_manifest.tsv")
    target_guides = sorted({(row["target_uid"], row["guide_id"]) for row in randomization})
    targets = sorted({target for target, _ in target_guides})
    guides_by_target: dict[str, set[str]] = defaultdict(set)
    for target, guide in target_guides:
        guides_by_target[target].add(guide)
    if any(len(guides) != 2 for guides in guides_by_target.values()):
        raise RuntimeError("Terminal adjudication requires exactly two frozen guides per target")
    phenotype_rows = read_tsv(execution_root / "frozen_phenotype_registry.tsv")
    phenotype_ids = sorted(row["phenotype_id"] for row in phenotype_rows)
    if len(phenotype_ids) != 2 or len(set(phenotype_ids)) != 2:
        raise RuntimeError("Stage-A execution does not contain exactly two phenotypes")
    execution_manifest = read_tsv(execution_root / "stage_a_execution_input_manifest.tsv")
    target_seal_row = next(
        row for row in execution_manifest if row["role"] == "stage_a_target_freeze_seal"
    )
    target_root = (PROJECT_ROOT / target_seal_row["source_path"]).parent
    target_registry = {
        row["target_uid"]: row
        for row in read_tsv(target_root / "frozen_target_registry.tsv")
        if row["selection_status"] == "selected_for_stage_a"
    }
    exact_editability: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in read_tsv(target_root / "stage_b_exact_editability_registry.tsv"):
        exact_editability[row["target_uid"]].append(row)
    if set(target_registry) != set(targets):
        raise RuntimeError("Execution target universe differs from the bound target freeze")
    pilot_rows = read_tsv(execution_root / "frozen_pilot_power_audit.tsv")
    if len(pilot_rows) != 1:
        raise RuntimeError("Bound Stage-A execution must contain one pilot schedule")
    pilot = pilot_rows[0]

    expected: dict[str, set[tuple[str, str, str]]] = {
        "CIS01": {(target, guide, "SOURCE_TARGET") for target, guide in target_guides},
        "RELAY01": {(target, guide, recipient) for target, guide in target_guides for recipient in RECIPIENTS},
        "ORIGIN01": {(target, guide, recipient) for target, guide in target_guides for recipient in RECIPIENTS},
        "ROUTE01": {(target, guide, recipient) for target, guide in target_guides for recipient in RECIPIENTS},
        "COMP01": {(target, guide, recipient) for target, guide in target_guides for recipient in RECIPIENTS},
        "PHENO01": {(target, guide, phenotype) for target, guide in target_guides for phenotype in phenotype_ids},
    }
    rows_by_estimand: dict[str, list[dict[str, str]]] = defaultdict(list)
    observed: dict[str, set[tuple[str, str, str]]] = defaultdict(set)
    for row in primary:
        estimand = row["estimand_id"]
        if estimand not in expected:
            raise RuntimeError(f"Unexpected Stage-A estimand: {estimand}")
        key = (row["target_uid"], row["guide_id"], row["endpoint_id"])
        if key in observed[estimand]:
            raise RuntimeError(f"Duplicate Stage-A estimand row: {estimand} {key}")
        observed[estimand].add(key)
        rows_by_estimand[estimand].append(row)
        number(row["effect"], "primary effect")
        number(row["se"], "primary SE")
        if int(row["n_biological_units"]) < int(row["n_backgrounds"]):
            raise RuntimeError("Stage-A biological-unit count is smaller than backgrounds")
        if not row["assay_native_model"] or row["result_status"] not in {"valid_estimate", "tested_null", "untestable_qc"}:
            raise RuntimeError("Stage-A result lacks an assay-native model/status")
    if observed != expected:
        raise RuntimeError("Stage-A primary estimand universe is incomplete or expanded")
    for estimand, rows in rows_by_estimand.items():
        validate_q_family(rows, estimand)

    sensitivity_observed: dict[tuple[str, str, str, str], dict[str, str]] = {}
    expected_sensitivity = {
        (sensitivity_id, target, guide, recipient)
        for sensitivity_id in SENSITIVITIES
        for target, guide in target_guides
        for recipient in RECIPIENTS
    }
    for row in sensitivity:
        key = (row["sensitivity_id"], row["target_uid"], row["guide_id"], row["recipient_lineage"])
        if key in sensitivity_observed:
            raise RuntimeError("Duplicate Stage-A composition sensitivity")
        sensitivity_observed[key] = row
        number(row["effect"], "sensitivity effect")
        if row["result_status"] not in {"valid_estimate", "tested_null", "untestable_qc"}:
            raise RuntimeError("Stage-A sensitivity status is invalid")
    if set(sensitivity_observed) != expected_sensitivity:
        raise RuntimeError("Stage-A sensitivity universe is incomplete or expanded")

    primary_index = {
        (row["estimand_id"], row["target_uid"], row["guide_id"], row["endpoint_id"]): row
        for row in primary
    }
    q_threshold = float(contract["primary_q_threshold"])
    min_backgrounds = int(contract["minimum_stage_a_backgrounds"])
    gate_rows = []
    eligible: list[dict[str, object]] = []
    initial_status_by_target: dict[str, str] = {}
    for target in targets:
        guides = sorted(guides_by_target[target])
        cis_pass = all(row_pass(primary_index[("CIS01", target, guide, "SOURCE_TARGET")], q_threshold, min_backgrounds) for guide in guides)
        recipient_evidence = []
        for recipient in RECIPIENTS:
            relay_rows = [primary_index[("RELAY01", target, guide, recipient)] for guide in guides]
            origin_rows = [primary_index[("ORIGIN01", target, guide, recipient)] for guide in guides]
            route_rows = [primary_index[("ROUTE01", target, guide, recipient)] for guide in guides]
            relay_pass = all(row_pass(row, q_threshold, min_backgrounds) for row in relay_rows)
            origin_pass = all(row_pass(row, q_threshold, min_backgrounds) for row in origin_rows)
            route_complete = all(row["result_status"] in {"valid_estimate", "tested_null"} for row in route_rows)
            composition_pass = all(
                sensitivity_observed[(sensitivity_id, target, guide, recipient)]["result_status"] in {"valid_estimate", "tested_null"}
                and yes(sensitivity_observed[(sensitivity_id, target, guide, recipient)]["expected_direction_agreement"])
                and int(sensitivity_observed[(sensitivity_id, target, guide, recipient)]["n_backgrounds"]) >= min_backgrounds
                for sensitivity_id in SENSITIVITIES for guide in guides
            )
            passed = cis_pass and relay_pass and origin_pass and composition_pass and route_complete
            worst_q = max(number(row["q_value"], "gate q") for row in relay_rows + origin_rows)
            recipient_evidence.append((recipient, passed, relay_pass, origin_pass, composition_pass, route_complete, worst_q))
            gate_rows.append({
                "target_uid": target, "recipient_lineage": recipient,
                "cis_both_guides_pass": str(cis_pass).lower(),
                "relay_both_guides_pass": str(relay_pass).lower(),
                "origin_both_guides_pass": str(origin_pass).lower(),
                "composition_sensitivities_pass": str(composition_pass).lower(),
                "route_estimable": str(route_complete).lower(),
                "stage_a_nomination_gate_pass": str(passed).lower(),
                "worst_claim_gate_q": f"{worst_q:.17g}",
            })
        passing = [item for item in recipient_evidence if item[1]]
        if passing:
            best = min(passing, key=lambda item: (item[6], item[0]))
            eligible.append({"target_uid": target, "recipient_lineage": best[0], "worst_gate_q": best[6]})
            initial_status_by_target[target] = "stage_a_gate_pass_pending_mediator"
        elif not cis_pass:
            initial_status_by_target[target] = "stage_a_cis_gate_failed"
        elif not any(item[2] for item in recipient_evidence):
            initial_status_by_target[target] = "stage_a_cell_autonomous_boundary"
        elif not any(item[2] and item[3] for item in recipient_evidence):
            initial_status_by_target[target] = "stage_a_propagation_not_source_specific"
        elif not any(item[2] and item[3] and item[4] for item in recipient_evidence):
            initial_status_by_target[target] = "stage_a_composition_mediated_boundary"
        else:
            initial_status_by_target[target] = "stage_a_route_indeterminate"

    if mediators:
        mediator_q = bh_adjust([number(row["p_value"], "mediator p") for row in mediators])
        for row, q_value in zip(mediators, mediator_q):
            if abs(number(row["q_value"], "mediator q") - q_value) > 1e-8:
                raise RuntimeError("Mediator BH q-value does not rederive")
            if row["result_status"] not in {"valid_estimate", "tested_null", "untestable_qc"}:
                raise RuntimeError("Mediator result status is invalid")
    mediator_candidates: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in mediators:
        passes = (
            row["result_status"] in {"valid_estimate", "tested_null"}
            and number(row["q_value"], "mediator q") < float(contract["mediator_q_threshold"])
            and yes(row["expected_direction_agreement"])
            and yes(row["intermediate_time_detected"])
            and yes(row["source_secretome_detected"])
            and yes(row["recipient_receptor_detected"])
            and int(row["n_backgrounds"]) >= min_backgrounds
            and yes(row["blocking_feasibility_pass"])
            and all(row[field] for field in [
                "blocking_reagent", "blocking_dose", "blocking_dose_unit",
                "upstream_cis_preservation_assay", "off_target_control",
            ])
        )
        if passes:
            mediator_candidates[(row["target_uid"], row["recipient_lineage"])].append(row)

    maximum = int(contract["maximum_stage_b_nominations"])
    ranked = sorted(
        eligible,
        key=lambda row: (float(row["worst_gate_q"]), str(row["target_uid"])),
    )[:maximum]
    nominations = []
    for row in ranked:
        candidates = mediator_candidates.get((str(row["target_uid"]), str(row["recipient_lineage"])), [])
        if not candidates:
            continue
        mediator = min(candidates, key=lambda item: (number(item["q_value"], "mediator q"), item["mediator_id"]))
        nominations.append({
            "nomination_rank": len(nominations) + 1,
            "target_uid": row["target_uid"],
            "recipient_lineage": row["recipient_lineage"],
            "worst_claim_gate_q": f"{float(row['worst_gate_q']):.17g}",
            "mediator_id": mediator["mediator_id"],
            "mediator_q": mediator["q_value"],
            "stage_b_nomination_status": "eligible_pending_exact_design_validation",
        })

    nominated_targets = {str(row["target_uid"]) for row in nominations}
    stage_b_files_have_rows = any([stage_b_design, stage_b_units, stage_b_randomization, stage_b_signoffs])
    if not nominated_targets and stage_b_files_have_rows:
        raise RuntimeError("Stage-B inputs are populated despite zero eligible nominations")
    stage_b_frozen = False
    if nominated_targets:
        if not all([stage_b_design, stage_b_units, stage_b_randomization]) or len(stage_b_signoffs) != 1:
            raise RuntimeError("Eligible Stage-B nominations require a complete exact-design freeze")
        if {row["target_uid"] for row in stage_b_design} != nominated_targets or len(stage_b_design) != len(nominated_targets):
            raise RuntimeError("Stage-B exact-design target universe differs from nominations")
        nomination_by_target = {str(row["target_uid"]): row for row in nominations}
        phenotype_class = {row["phenotype_id"]: row["phenotype_class"] for row in phenotype_rows}
        for row in stage_b_design:
            nomination = nomination_by_target[row["target_uid"]]
            if row["recipient_lineage"] != nomination["recipient_lineage"] or row["mediator_id"] != nomination["mediator_id"]:
                raise RuntimeError("Stage-B design substitutes a nominated recipient/mediator")
            frozen_target = target_registry[row["target_uid"]]
            role = frozen_target["selected_exact_reference_to_alternate_role"]
            reference = frozen_target["selected_exact_reference_allele"]
            alternate = frozen_target["selected_exact_alternate_allele"]
            if role == "install_risk_from_protective_reference":
                expected_protective, expected_risk = reference, alternate
            elif role == "install_protective_from_risk_reference":
                expected_protective, expected_risk = alternate, reference
            else:
                raise RuntimeError("Frozen target has unresolved exact-edit allele orientation")
            if (
                row["snp_hg19"] != frozen_target["selected_exact_snp_hg19"]
                or row["orientation_uid"] != frozen_target["best_orientation_uid"]
                or row["exact_design_uid"] != frozen_target["selected_exact_design_uid"]
                or row["protective_allele"] != expected_protective
                or row["risk_allele"] != expected_risk
            ):
                raise RuntimeError("Stage-B design substitutes the frozen variant/orientation/alleles")
            if (
                row["source_lineage"] != frozen_target["dominant_celltype"]
                or row["challenge_id"] != pilot["challenge_id"]
                or any(
                    row[field] != pilot[field]
                    for field in [
                        "early_time_value", "early_time_unit",
                        "intermediate_time_value", "intermediate_time_unit",
                        "late_time_value", "late_time_unit",
                    ]
                )
            ):
                raise RuntimeError("Stage-B design changes the frozen lineage/challenge/time schedule")
            eligible_pegrnas = {
                item["pegrna_id"] for item in exact_editability[row["target_uid"]]
                if item["design_uid"] == row["exact_design_uid"]
            }
            if (
                row["editing_method"] != "prime_editing"
                or row["pegrna_id_1"] == row["pegrna_id_2"]
                or not {row["pegrna_id_1"], row["pegrna_id_2"]} <= eligible_pegrnas
                or not re.fullmatch(r"[0-9a-f]{64}", row["edit_design_sha256"])
            ):
                raise RuntimeError("Stage-B exact-edit design is not bound to two preverified pegRNAs")
            if row["protective_allele"] == row["risk_allele"] or not all(row[field] for field in [
                "variant_id", "genome_build", "editing_method", "edit_design_sha256",
                "source_lineage", "challenge_id", "early_time_value", "early_time_unit",
                "intermediate_time_value", "intermediate_time_unit", "late_time_value",
                "late_time_unit", "rescue_method", "rescue_direction", "blocking_reagent",
                "blocking_dose", "blocking_dose_unit",
            ]):
                raise RuntimeError("Stage-B exact allele/rescue/blockade design is incomplete")
            if row["rescue_direction"] != "restore_protective_target_dosage":
                raise RuntimeError("Stage-B rescue does not restore protective target dosage")
            selected_mediator = min(
                mediator_candidates[(row["target_uid"], row["recipient_lineage"])],
                key=lambda item: (number(item["q_value"], "mediator q"), item["mediator_id"]),
            )
            if any(row[field] != selected_mediator[field] for field in ["blocking_reagent", "blocking_dose", "blocking_dose_unit"]):
                raise RuntimeError("Stage-B blocker differs from the frozen mediator candidate")
            phenotypes = [row["phenotype_id_1"], row["phenotype_id_2"]]
            if set(phenotypes) != set(phenotype_ids) or phenotype_class[phenotypes[0]] == phenotype_class[phenotypes[1]]:
                raise RuntimeError("Stage-B design does not preserve both distinct frozen phenotypes")
            if not yes(row["target_removed_from_axes"]) or not yes(row["programs_with_target_removed"]) or yes(row["outcomes_inspected_before_design_freeze"]):
                raise RuntimeError("Stage-B design violates score-removal/outcome firewall")

        units_by_target: dict[str, list[dict[str, str]]] = defaultdict(list)
        for row in stage_b_units:
            if row["target_uid"] not in nominated_targets or not yes(row["independent_clone"]) or not yes(row["independent_differentiation"]):
                raise RuntimeError("Stage-B power/unit table contains unauthorized or dependent units")
            if not (0 < number(row["power_family_alpha"], "Stage-B alpha") <= 0.05) or number(row["target_power"], "Stage-B power") < 0.8:
                raise RuntimeError("Stage-B power audit is below the frozen confirmatory requirement")
            if not row["pilot_variance_source"] or int(row["maximum_feasible_expansion"]) < 0:
                raise RuntimeError("Stage-B power audit lacks variance/expansion bounds")
            units_by_target[row["target_uid"]].append(row)
        min_backgrounds_b = int(contract["minimum_stage_b_backgrounds"])
        min_clones = int(contract["minimum_clones_per_allele_background"])
        min_differentiations = int(contract["minimum_differentiations_per_clone"])
        for target in nominated_targets:
            rows = units_by_target[target]
            backgrounds = {row["background_id"] for row in rows}
            if len(backgrounds) < min_backgrounds_b:
                raise RuntimeError("Stage-B unit freeze has fewer than three backgrounds")
            for background in backgrounds:
                for allele in ("protective", "risk"):
                    allele_rows = [row for row in rows if row["background_id"] == background and row["allele_role"] == allele]
                    clones = {row["clone_id"] for row in allele_rows}
                    if len(clones) < min_clones:
                        raise RuntimeError("Stage-B unit freeze lacks two clones per allele/background")
                    for clone in clones:
                        if len({row["differentiation_id"] for row in allele_rows if row["clone_id"] == clone}) < min_differentiations:
                            raise RuntimeError("Stage-B unit freeze lacks two differentiations per clone")
                if not {"parental", "sham"} <= {row["allele_role"] for row in rows if row["background_id"] == background}:
                    raise RuntimeError("Stage-B unit freeze lacks parental/sham controls")

        random_ids = {row["sample_id"] for row in stage_b_randomization}
        if len(random_ids) != len(stage_b_randomization) or {row["target_uid"] for row in stage_b_randomization} != nominated_targets:
            raise RuntimeError("Stage-B randomization sample/target universe drift")
        randomized_units = {
            (row["target_uid"], row["background_id"], row["allele_role"], row["clone_id"], row["differentiation_id"])
            for row in stage_b_randomization
        }
        frozen_units = {
            (row["target_uid"], row["background_id"], row["allele_role"], row["clone_id"], row["differentiation_id"])
            for row in stage_b_units
        }
        if not frozen_units <= randomized_units:
            raise RuntimeError("Stage-B randomization does not cover every frozen biological unit")
        for target in nominated_targets:
            if not STAGE_B_ARMS <= {row["arm_id"] for row in stage_b_randomization if row["target_uid"] == target}:
                raise RuntimeError("Stage-B randomization lacks a mandatory rescue/blockade/control arm")
        for row in stage_b_randomization:
            if (
                not re.fullmatch(r"BLIND_[A-Z0-9]{8,}", row["blinded_label"])
                or not yes(row["randomized_allocation"])
                or row["unblinding_status"] != "sealed"
                or row["preoutcome_exclusion_status"] != "included_preoutcome"
            ):
                raise RuntimeError("Stage-B randomization/blinding firewall failed")

        freeze_signoff = stage_b_signoffs[0]
        for path_field, hash_field, path in [
            ("stage_b_design_path", "stage_b_design_sha256", input_root / "stage_b_exact_design.tsv"),
            ("power_units_path", "power_units_sha256", input_root / "stage_b_power_units.tsv"),
            ("randomization_path", "randomization_sha256", input_root / "stage_b_randomization.tsv"),
        ]:
            if Path(freeze_signoff[path_field]).as_posix() != path.relative_to(PROJECT_ROOT).as_posix() or freeze_signoff[hash_field] != sha256_file(path):
                raise RuntimeError("Stage-B freeze signoff is not bound to exact inputs")
        if (
            not freeze_signoff["reviewer_1"] or not freeze_signoff["reviewer_2"]
            or freeze_signoff["reviewer_1"] == freeze_signoff["reviewer_2"]
            or not yes(freeze_signoff["review_concordant"])
            or yes(freeze_signoff["stage_b_outcomes_inspected"])
            or utc(freeze_signoff["signed_utc"]) <= utc(signoff["signed_utc"])
        ):
            raise RuntimeError("Stage-B freeze signoff is invalid or predates Stage-A adjudication")
        stage_b_frozen = True

    for row in nominations:
        row["stage_b_nomination_status"] = (
            "frozen_ready_for_stage_b_execution" if stage_b_frozen
            else "not_frozen_incomplete_stage_b_design"
        )
    status_by_target = dict(initial_status_by_target)
    for row in eligible:
        status_by_target[str(row["target_uid"])] = "stage_a_pass_no_freezable_mediator"
    for row in nominations:
        status_by_target[str(row["target_uid"])] = row["stage_b_nomination_status"]
    target_verdict_rows = [
        {"target_uid": target, "terminal_status": status_by_target[target]}
        for target in targets
    ]

    CANDIDATE_ROOT.mkdir(parents=True)
    outputs: list[Path] = []
    gate_path = CANDIDATE_ROOT / "stage_a_gate_matrix.tsv"
    write_tsv(gate_path, gate_rows, list(gate_rows[0]))
    outputs.append(gate_path)
    target_path = CANDIDATE_ROOT / "stage_a_target_terminal_verdict.tsv"
    write_tsv(target_path, target_verdict_rows, list(target_verdict_rows[0]))
    outputs.append(target_path)
    nomination_path = CANDIDATE_ROOT / "stage_a_stage_b_nominations.tsv"
    nomination_fields = ["nomination_rank", "target_uid", "recipient_lineage", "worst_claim_gate_q", "mediator_id", "mediator_q", "stage_b_nomination_status"]
    write_tsv(nomination_path, nominations, nomination_fields)
    outputs.append(nomination_path)
    if stage_b_frozen:
        for basename, rows in [
            ("frozen_stage_b_exact_design.tsv", stage_b_design),
            ("frozen_stage_b_power_units.tsv", stage_b_units),
            ("frozen_stage_b_randomization.tsv", stage_b_randomization),
            ("frozen_stage_b_signoff.tsv", stage_b_signoffs),
        ]:
            path = CANDIDATE_ROOT / basename
            source_name = basename.removeprefix("frozen_")
            write_tsv(path, rows, header(input_root / source_name))
            outputs.append(path)
    input_manifest = []
    for role, path in [
        ("stage_a_execution_seal", execution_path),
        ("stage_a_unblinding_seal", unblinding_path),
        ("stage_a_terminal_contract_seal", contract_path),
    ] + [(name.removesuffix(".tsv"), input_root / name) for name in [
        "stage_a_primary_estimand_results.tsv", "stage_a_sensitivity_results.tsv",
        "stage_a_mediator_candidates.tsv", "stage_a_outcome_signoff.tsv",
        "stage_b_exact_design.tsv", "stage_b_power_units.tsv",
        "stage_b_randomization.tsv", "stage_b_freeze_signoff.tsv",
    ]]:
        input_manifest.append({
            "role": role, "source_path": str(path.relative_to(PROJECT_ROOT)),
            "size_bytes": path.stat().st_size, "sha256": sha256_file(path),
        })
    manifest_path = CANDIDATE_ROOT / "stage_a_terminal_input_manifest.tsv"
    write_tsv(manifest_path, input_manifest, ["role", "source_path", "size_bytes", "sha256"])
    outputs.append(manifest_path)
    status = (
        "stage_a_terminal_stage_b_nomination_frozen" if stage_b_frozen
        else "stage_a_terminal_pass_no_freezable_mediator"
        if eligible else "stage_a_terminal_no_stage_b_nomination"
    )
    payload = {
        "status": status,
        "created_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "n_stage_a_targets": len(targets),
        "n_stage_a_gate_eligible_targets": len(eligible),
        "n_stage_b_nominations": len(nominations),
        "stage_a_terminal_verdict_frozen": True,
        "stage_b_design_frozen": stage_b_frozen,
        "stage_b_promotion_authorized": False,
        "scientific_outcomes_inspected": True,
        "replacement_target_search_permitted": False,
        "next_gate": (
            "independent exact-edit Stage-B execution under the frozen design"
            if stage_b_frozen else "terminal stop; report complete Stage-A result without replacement search"
        ),
        "output_sha256": {path.name: sha256_file(path) for path in outputs},
    }
    atomic_write_json(CANDIDATE_ROOT / "STAGE_A_TERMINAL_VERDICT.json", payload)
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
