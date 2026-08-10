#!/usr/bin/env python3
"""Independently rederive the Stage-A terminal and Stage-B nomination verdict."""

from __future__ import annotations

import json
from collections import defaultdict

from relay_common import CANDIDATE_ROOT, PROJECT_ROOT, read_tsv, sha256_file


RECIPIENTS = ("Cholangiocytes", "Endothelial_cells", "Fibroblasts", "Macrophages")
SENSITIVITIES = ("substate_adjusted", "fixed_cell_count")


def yes(value: object) -> bool:
    return str(value).strip().lower() in {"true", "t", "1", "yes"}


def bh(p_values: list[float]) -> list[float]:
    n = len(p_values)
    order = sorted(range(n), key=lambda index: (p_values[index], index))
    result = [1.0] * n
    running = 1.0
    for position in range(n - 1, -1, -1):
        index = order[position]
        running = min(running, p_values[index] * n / (position + 1))
        result[index] = min(1.0, running)
    return result


def passes(row: dict[str, str], q_threshold: float, minimum: int) -> bool:
    return (
        row["result_status"] in {"valid_estimate", "tested_null"}
        and float(row["q_value"]) < q_threshold
        and yes(row["expected_direction_agreement"])
        and int(row["n_backgrounds"]) >= minimum
    )


def main() -> None:
    seal = json.loads((CANDIDATE_ROOT / "STAGE_A_TERMINAL_VERDICT.json").read_text(encoding="utf-8"))
    if seal.get("status") not in {
        "stage_a_terminal_stage_b_nomination_frozen",
        "stage_a_terminal_pass_no_freezable_mediator",
        "stage_a_terminal_no_stage_b_nomination",
    }:
        raise RuntimeError("Invalid Stage-A terminal status")
    if seal.get("stage_a_terminal_verdict_frozen") is not True:
        raise RuntimeError("Stage-A terminal verdict is not frozen")
    if seal.get("stage_b_promotion_authorized") is not False or seal.get("replacement_target_search_permitted") is not False:
        raise RuntimeError("Terminal release improperly authorizes promotion/replacement search")
    for name, expected in seal["output_sha256"].items():
        path = CANDIDATE_ROOT / name
        if not path.is_file() or sha256_file(path) != expected:
            raise RuntimeError(f"Stage-A terminal output hash mismatch: {name}")

    manifest = read_tsv(CANDIDATE_ROOT / "stage_a_terminal_input_manifest.tsv")
    expected_roles = {
        "stage_a_execution_seal", "stage_a_unblinding_seal",
        "stage_a_terminal_contract_seal", "stage_a_primary_estimand_results",
        "stage_a_sensitivity_results", "stage_a_mediator_candidates",
        "stage_a_outcome_signoff", "stage_b_exact_design", "stage_b_power_units",
        "stage_b_randomization", "stage_b_freeze_signoff",
    }
    if {row["role"] for row in manifest} != expected_roles:
        raise RuntimeError("Stage-A terminal input-manifest universe drift")
    for row in manifest:
        path = PROJECT_ROOT / row["source_path"]
        if not path.is_file() or path.stat().st_size != int(row["size_bytes"]) or sha256_file(path) != row["sha256"]:
            raise RuntimeError(f"Stage-A terminal source drift: {path}")
    paths = {row["role"]: PROJECT_ROOT / row["source_path"] for row in manifest}
    contract = json.loads(paths["stage_a_terminal_contract_seal"].read_text(encoding="utf-8"))
    execution_root = paths["stage_a_execution_seal"].parent
    randomization = read_tsv(execution_root / "frozen_randomization_manifest.tsv")
    target_guides = sorted({(row["target_uid"], row["guide_id"]) for row in randomization})
    targets = sorted({target for target, _ in target_guides})
    guides_by_target: dict[str, set[str]] = defaultdict(set)
    for target, guide in target_guides:
        guides_by_target[target].add(guide)
    phenotype_ids = sorted(
        row["phenotype_id"]
        for row in read_tsv(execution_root / "frozen_phenotype_registry.tsv")
    )
    expected_sizes = {
        "CIS01": len(target_guides),
        "RELAY01": len(target_guides) * len(RECIPIENTS),
        "ORIGIN01": len(target_guides) * len(RECIPIENTS),
        "ROUTE01": len(target_guides) * len(RECIPIENTS),
        "COMP01": len(target_guides) * len(RECIPIENTS),
        "PHENO01": len(target_guides) * len(phenotype_ids),
    }
    primary = read_tsv(paths["stage_a_primary_estimand_results"])
    by_family: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in primary:
        by_family[row["estimand_id"]].append(row)
    if set(by_family) != set(expected_sizes):
        raise RuntimeError("Independent primary-estimand family universe drift")
    for family, rows in by_family.items():
        if len(rows) != expected_sizes[family]:
            raise RuntimeError(f"Independent {family} family size drift")
        adjusted = bh([float(row["p_value"]) for row in rows])
        if any(abs(float(row["q_value"]) - q_value) > 1e-8 for row, q_value in zip(rows, adjusted)):
            raise RuntimeError(f"Independent {family} BH rederivation failed")
    primary_index = {
        (row["estimand_id"], row["target_uid"], row["guide_id"], row["endpoint_id"]): row
        for row in primary
    }
    sensitivities = {
        (row["sensitivity_id"], row["target_uid"], row["guide_id"], row["recipient_lineage"]): row
        for row in read_tsv(paths["stage_a_sensitivity_results"])
    }
    expected_sensitivity = {
        (sensitivity_id, target, guide, recipient)
        for sensitivity_id in SENSITIVITIES
        for target, guide in target_guides for recipient in RECIPIENTS
    }
    if set(sensitivities) != expected_sensitivity:
        raise RuntimeError("Independent sensitivity universe drift")
    threshold = float(contract["primary_q_threshold"])
    minimum = int(contract["minimum_stage_a_backgrounds"])
    expected_gate = {}
    eligible = []
    for target in targets:
        guides = sorted(guides_by_target[target])
        cis = all(passes(primary_index[("CIS01", target, guide, "SOURCE_TARGET")], threshold, minimum) for guide in guides)
        target_candidates = []
        for recipient in RECIPIENTS:
            relay_rows = [primary_index[("RELAY01", target, guide, recipient)] for guide in guides]
            origin_rows = [primary_index[("ORIGIN01", target, guide, recipient)] for guide in guides]
            relay = all(passes(row, threshold, minimum) for row in relay_rows)
            origin = all(passes(row, threshold, minimum) for row in origin_rows)
            route = all(primary_index[("ROUTE01", target, guide, recipient)]["result_status"] in {"valid_estimate", "tested_null"} for guide in guides)
            composition = all(
                sensitivities[(sensitivity_id, target, guide, recipient)]["result_status"] in {"valid_estimate", "tested_null"}
                and yes(sensitivities[(sensitivity_id, target, guide, recipient)]["expected_direction_agreement"])
                and int(sensitivities[(sensitivity_id, target, guide, recipient)]["n_backgrounds"]) >= minimum
                for sensitivity_id in SENSITIVITIES for guide in guides
            )
            passed = cis and relay and origin and route and composition
            worst_q = max(float(row["q_value"]) for row in relay_rows + origin_rows)
            expected_gate[(target, recipient)] = (cis, relay, origin, composition, route, passed, worst_q)
            if passed:
                target_candidates.append((worst_q, recipient))
        if target_candidates:
            worst_q, recipient = min(target_candidates)
            eligible.append((worst_q, target, recipient))
    observed_gate = {
        (row["target_uid"], row["recipient_lineage"]): row
        for row in read_tsv(CANDIDATE_ROOT / "stage_a_gate_matrix.tsv")
    }
    if set(observed_gate) != set(expected_gate):
        raise RuntimeError("Independent gate-matrix universe drift")
    for key, values in expected_gate.items():
        row = observed_gate[key]
        booleans = [
            "cis_both_guides_pass", "relay_both_guides_pass",
            "origin_both_guides_pass", "composition_sensitivities_pass",
            "route_estimable", "stage_a_nomination_gate_pass",
        ]
        if [yes(row[field]) for field in booleans] != list(values[:6]) or abs(float(row["worst_claim_gate_q"]) - values[6]) > 1e-8:
            raise RuntimeError("Independent gate-matrix rederivation failed")

    mediator_rows = read_tsv(paths["stage_a_mediator_candidates"])
    if mediator_rows:
        adjusted = bh([float(row["p_value"]) for row in mediator_rows])
        if any(abs(float(row["q_value"]) - q_value) > 1e-8 for row, q_value in zip(mediator_rows, adjusted)):
            raise RuntimeError("Independent mediator BH rederivation failed")
    mediator_pass: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in mediator_rows:
        if (
            row["result_status"] in {"valid_estimate", "tested_null"}
            and float(row["q_value"]) < float(contract["mediator_q_threshold"])
            and yes(row["expected_direction_agreement"])
            and yes(row["intermediate_time_detected"])
            and yes(row["source_secretome_detected"])
            and yes(row["recipient_receptor_detected"])
            and int(row["n_backgrounds"]) >= minimum
            and yes(row["blocking_feasibility_pass"])
            and all(row[field] for field in ["blocking_reagent", "blocking_dose", "blocking_dose_unit", "upstream_cis_preservation_assay", "off_target_control"])
        ):
            mediator_pass[(row["target_uid"], row["recipient_lineage"])].append(row)
    top = sorted(eligible)[: int(contract["maximum_stage_b_nominations"])]
    expected_nominations = []
    for worst_q, target, recipient in top:
        candidates = mediator_pass.get((target, recipient), [])
        if candidates:
            mediator = min(candidates, key=lambda row: (float(row["q_value"]), row["mediator_id"]))
            expected_nominations.append((target, recipient, mediator["mediator_id"]))
    observed_nominations = [
        (row["target_uid"], row["recipient_lineage"], row["mediator_id"])
        for row in read_tsv(CANDIDATE_ROOT / "stage_a_stage_b_nominations.tsv")
    ]
    if observed_nominations != expected_nominations:
        raise RuntimeError("Independent Stage-B nomination rederivation failed")
    if int(seal["n_stage_a_gate_eligible_targets"]) != len(eligible) or int(seal["n_stage_b_nominations"]) != len(expected_nominations):
        raise RuntimeError("Terminal verdict count drift")
    stage_b_frozen = seal.get("stage_b_design_frozen") is True
    if stage_b_frozen != (seal["status"] == "stage_a_terminal_stage_b_nomination_frozen"):
        raise RuntimeError("Stage-B frozen/status inconsistency")
    for basename in [
        "frozen_stage_b_exact_design.tsv", "frozen_stage_b_power_units.tsv",
        "frozen_stage_b_randomization.tsv", "frozen_stage_b_signoff.tsv",
    ]:
        if (CANDIDATE_ROOT / basename).exists() != stage_b_frozen:
            raise RuntimeError("Frozen Stage-B output presence/status mismatch")
    print(
        "STAGE_A_TERMINAL_VALIDATION_PASS "
        f"targets={len(targets)} eligible={len(eligible)} nominations={len(expected_nominations)} "
        f"stage_b_frozen={str(stage_b_frozen).lower()} promotion=false"
    )


if __name__ == "__main__":
    main()
