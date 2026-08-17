"""Terminal, fail-closed adoption decision for the MASLD continual integration."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path


POLICY_SCHEMA = "masld-cl-final-adoption-policy-v42"
EXPECTED_CONTRACT = {
    "strict_reference": {"cells": 216957, "donors": 7, "libraries": 29},
    "primary_query": {
        "cells": 687559, "donors": 64, "libraries": 184,
        "controls": 9, "cases": 55,
    },
    "analyzed": {"cells": 1232285, "donors": 102, "libraries": 273},
    "descriptive": {"cells": 1232318, "donors": 104, "libraries": 275},
}


def _gate(
    gate_id: str, category: str, observed: Any, threshold: Any, passed: bool,
    reason: str,
) -> dict[str, Any]:
    return {
        "gate_id": gate_id,
        "category": category,
        "observed": observed,
        "threshold": threshold,
        "passed": bool(passed),
        "reason": reason,
    }


def _load_policy(config: dict[str, Any], value: str | Path):
    path = Path(value).resolve()
    expected = (
        Path(config["_config_path"]).resolve().parent
        / "reference" / "final_adoption_policy_v42.json"
    )
    if path != expected:
        raise ContractError("final adoption requires its source-controlled V42 policy")
    with path.open() as handle:
        policy = json.load(handle)
    if (
        policy.get("schema_version") != POLICY_SCHEMA
        or policy.get("config_sha256") != config["_config_sha256"]
        or policy.get("decision_rule")
        != "primary promotion requires every required gate; any integrity or preservation failure leaves Harmony primary"
        or policy.get("external_failure_may_trigger_revision_on_gse212837") is not False
    ):
        raise ContractError("final adoption policy identity differs")
    artifacts = {}
    for name, spec in policy.get("sources", {}).items():
        source = (path.parent / spec["path"]).resolve()
        if sha256_path(source) != spec["sha256"]:
            raise ContractError(f"final adoption source changed: {name}")
        with source.open() as handle:
            artifact = json.load(handle)
        if artifact.get("schema_version") != spec["schema_version"]:
            raise ContractError(f"final adoption source schema differs: {name}")
        if artifact.get("config_sha256") != config["_config_sha256"]:
            raise ContractError(f"final adoption config hash differs: {name}")
        artifacts[name] = (source, artifact)
    if set(artifacts) != set(policy["sources"]):
        raise ContractError("final adoption source roster differs")
    return path, policy, artifacts


def parse_unittest_log(path_value: str | Path, policy: dict[str, Any]) -> dict[str, Any]:
    path = Path(path_value).resolve()
    text = path.read_text()
    match = re.search(r"Ran\s+(\d+)\s+tests?\s+in\s+[0-9.]+s", text)
    terminal = re.search(r"^OK(?:\s+\(skipped=(\d+)\))?\s*$", text, re.MULTILINE)
    if match is None or terminal is None:
        raise ContractError("software-test log is incomplete or not successful")
    tests = int(match.group(1))
    skipped = int(terminal.group(1) or 0)
    requirements = policy["software_tests"]
    return {
        "path": str(path),
        "sha256": sha256_path(path),
        "tests_run": tests,
        "skipped": skipped,
        "failures": 0,
        "errors": 0,
        "terminal_status": "OK",
        "passed": bool(
            tests >= requirements["minimum_tests"]
            and skipped <= requirements["maximum_skips"]
        ),
    }


def _build_gates(artifacts: dict[str, tuple[Path, dict[str, Any]]], tests):
    get = lambda name: artifacts[name][1]
    contract = get("contract")
    confirmation = get("internal_confirmation")
    selection = get("selection_lock")
    outcomes = get("outcome_decision")
    reference = get("reference_decision")
    projection = get("full_atlas_projection")
    production = get("production_expansion")
    stress = get("secondary_stress")
    firewall = get("program_firewall")
    external = get("independent_external")
    gpu = get("gpu_tolerance")

    routing = external["routing_audit"]
    external_gates = external["gates"]
    gpu_reproducible = all([
        gpu["same_source_identity"], gpu["same_input_artifacts"],
        gpu["same_arguments_except_output"], gpu["same_cell_roster_and_order"],
        gpu["full_embedding"]["within_absolute_1e_6"],
        gpu["training_embedding"]["within_absolute_1e_6"],
        gpu["checkpoint"]["bitwise_identical_tensors"],
    ])
    production_roster = {
        "all_lineage", "hepatocytes", "macrophages", "fibroblasts",
        "cholangiocytes", "t_cells",
    }
    return [
        _gate(
            "atlas_contract", "integrity", contract["observed_contract"],
            EXPECTED_CONTRACT, contract["production_ready"]
            and contract["observed_contract"] == EXPECTED_CONTRACT,
            "The immutable census and integer raw-count contract must match exactly.",
        ),
        _gate(
            "control_only_selection_lock", "integrity",
            {"selection_frozen": selection["selection_frozen"],
             "outcomes_unlocked": selection["outcomes_unlocked"]},
            {"selection_frozen": True, "outcomes_unlocked": True},
            selection["selection_frozen"] and selection["outcomes_unlocked"],
            "The frozen control-only selection lock must explicitly authorize outcome evaluation.",
        ),
        _gate(
            "internal_confirmation", "preservation_and_alignment",
            confirmation["confirmation_pass"], True,
            confirmation["confirmation_pass"],
            "The five-seed, held-study, label, disease, and order matrix passed.",
        ),
        _gate(
            "outcome_lock", "disease_preservation",
            {"disease_preservation_pass": outcomes["disease_preservation_pass"],
             "decision": outcomes["decision"]},
            {"disease_preservation_pass": True,
             "decision": "continue_remaining_gates"},
            outcomes["disease_preservation_pass"]
            and outcomes["decision"] == "continue_remaining_gates",
            "Locked donor-level raw-count-PCA comparators preserve cross-sectional disease variation.",
        ),
        _gate(
            "healthy_reference_choice", "reference",
            reference["decision"], "retain_common_strict7",
            reference["decision"] == "retain_common_strict7",
            "The strict reference already spans two datasets; the tested 26-donor expansion failed preservation.",
        ),
        _gate(
            "full_atlas_projection", "production",
            {"contract": projection["contract"],
             "v32_coordinates_bitwise_identical": projection["v32_coordinates_bitwise_identical"],
             "secondary_cells_used_for_fit": projection["secondary_and_descriptive_cells_used_for_fit"]},
            "1,232,318 descriptive cells; primary coordinates unchanged; secondary cells projection-only",
            projection["contract"]["descriptive_cells"] == 1232318
            and projection["contract"]["analyzed_cells"] == 1232285
            and projection["v32_coordinates_bitwise_identical"]
            and not projection["secondary_and_descriptive_cells_used_for_fit"],
            "The complete frozen census was projected without refitting on secondary or descriptive-only cells.",
        ),
        _gate(
            "six_model_production_expansion", "production",
            {"models": sorted(production["model_roster"]),
             "conditions_available_to_fit": production["conditions_available_to_fit"]},
            "all-lineage plus five lineages; no conditions available to fit",
            set(production["model_roster"]) == production_roster
            and set(production["models"]) == production_roster
            and not production["conditions_available_to_fit"],
            "The stage-blind production bundle contains every prespecified model.",
        ),
        _gate(
            "secondary_technical_stress", "robustness",
            {"technical_stress_pass": stress["technical_stress_pass"],
             "independent_confirmation_pass": stress["independent_secondary_confirmation_pass"]},
            "technical stress passes; independence reported separately",
            stress["technical_stress_pass"],
            "GSE136103 passes technically but is not independent; GSE174748 is underpowered.",
        ),
        _gate(
            "frozen_program_firewall", "program_firewall",
            {"passed": firewall["passed"], "programs": firewall["registry"]["programs"]},
            {"passed": True, "programs": 117},
            firewall["passed"] and firewall["registry"]["programs"] == 117,
            "Program identities, weights, scores, and inferential families remain frozen.",
        ),
        _gate(
            "gpu_reproducibility_recorded", "reproducibility",
            {"maximum_absolute_difference": gpu["observed_gpu_tolerance"],
             "checkpoint_tensors_bitwise_identical": gpu["checkpoint"]["bitwise_identical_tensors"]},
            "measured and <= 1e-6 on independent supported GPUs",
            gpu_reproducible,
            "Two separate L40S executions of the exact paper setting were bit-identical.",
        ),
        _gate(
            "independent_external_alignment", "external_validation",
            external_gates["all_required_control_alignment_gates"], True,
            external_gates["all_required_control_alignment_gates"],
            "All four evaluable GSE212837 models passed locked control-alignment gates.",
        ),
        _gate(
            "independent_external_disease_preservation", "external_validation",
            external_gates["all_required_disease_preservation_gates"], True,
            external_gates["all_required_disease_preservation_gates"],
            "All four evaluable GSE212837 models preserved donor-level disease geometry.",
        ),
        _gate(
            "independent_external_routing", "preservation",
            routing["mean_donor_macro_f1"], routing["minimum"], routing["pass"],
            "Untouched GSE212837 broad-label routing failed the locked donor-balanced macro-F1 floor.",
        ),
        _gate(
            "software_regression_suite", "integrity", tests,
            "successful suite meeting the V42 minimum test and skip limits",
            tests["passed"], "The final source tree must pass its compute-node regression suite.",
        ),
    ]


def write_final_adoption_decision(
    config: dict[str, Any], policy_value: str | Path, test_log: str | Path,
    output_value: str | Path,
) -> dict[str, Any]:
    policy_path, policy, artifacts = _load_policy(config, policy_value)
    tests = parse_unittest_log(test_log, policy)
    gates = _build_gates(artifacts, tests)
    failed = [gate for gate in gates if not gate["passed"]]
    implementation_gate_ids = {
        gate["gate_id"] for gate in gates
        if gate["gate_id"] != "independent_external_routing"
    }
    implementation_complete = all(
        gate["passed"] for gate in gates if gate["gate_id"] in implementation_gate_ids
    )
    external = artifacts["independent_external"][1]
    reference = artifacts["reference_decision"][1]
    decision = {
        "schema_version": "masld-cl-final-adoption-decision-v42",
        "config_sha256": config["_config_sha256"],
        "decision": "reject_primary_promotion" if failed else "eligible_pending_human_approval",
        "implementation_complete": implementation_complete,
        "promotion_eligible": not failed,
        "all_required_promotion_gates_pass": not failed,
        "harmony_remains_primary": bool(failed),
        "candidate_role": "supplementary_sensitivity_only" if failed else "promotion_candidate",
        "paper_claim_allowed": not failed,
        "method_identity": policy["method_identity"],
        "failure_reasons": [gate["reason"] for gate in failed],
        "failed_gate_ids": [gate["gate_id"] for gate in failed],
        "gates": gates,
        "independent_external_interpretation": {
            "cohort": external["cohort"],
            "alignment_pass": external["gates"]["all_required_control_alignment_gates"],
            "disease_preservation_pass": external["gates"]["all_required_disease_preservation_gates"],
            "routing_macro_f1": external["routing_audit"]["mean_donor_macro_f1"],
            "routing_threshold": external["routing_audit"]["minimum"],
            "method_revision_on_this_cohort_allowed": policy["external_failure_may_trigger_revision_on_gse212837"],
        },
        "healthy_reference_decision": {
            "selected": "strict7",
            "datasets": 2,
            "donors": 7,
            "external26_expansion_decision": reference["decision"],
        },
        "release_actions": {
            "canonical_h5ad_changed": False,
            "paper_authorities_changed": False,
            "manuscript_changed": False,
            "promotion_approval_requested": False,
        },
        "policy": {"path": str(policy_path), "sha256": sha256_path(policy_path)},
        "software_tests": tests,
        "sources": {
            name: {"path": str(path), "sha256": sha256_path(path)}
            for name, (path, _) in artifacts.items()
        },
    }
    output = Path(output_value).resolve()
    output.mkdir(parents=True, exist_ok=False)
    write_json_exclusive(output / "promotion_decision.json", decision)
    return decision
