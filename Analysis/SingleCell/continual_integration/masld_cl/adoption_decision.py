"""Fail-closed adoption decision for the replay-plus-EWC and V9 pilots."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .firewall import load_retargeted_bridge_selection_lock, validate_program_firewall


def _read(path: str | Path) -> tuple[Path, dict[str, Any]]:
    resolved = Path(path).resolve()
    with resolved.open() as handle:
        return resolved, json.load(handle)


def _gate(
    gate_id: str, category: str, observed: Any, threshold: str, passed: bool,
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


def _minimum_held_protocol_harmony(held: dict[str, Any]) -> dict[str, Any]:
    rows = []
    for study, result in held["held_studies"].items():
        for row in result["held_control_rows"]:
            if row["testable"]:
                rows.append({
                    "study": study,
                    "model_kind": row["model_kind"],
                    "improvement": float(row["protocol_harmony"]["improvement"]),
                })
    if not rows:
        raise ContractError("held-study artifact has no evaluable protocol/Harmony rows")
    return min(rows, key=lambda row: row["improvement"])


def _check_v9_identity(
    artifact: dict[str, Any], path: Path, config_sha256: str,
    selection: dict[str, Any], selection_path: Path,
) -> None:
    if artifact.get("config_sha256") != config_sha256:
        raise ContractError(f"config hash differs in {path}")
    if artifact.get("selection_lock_sha256") != selection["lock_sha256"]:
        raise ContractError(f"selection lock hash differs in {path}")
    if artifact.get("selection_lock_file_sha256") != sha256_path(selection_path):
        raise ContractError(f"selection lock file hash differs in {path}")


def write_adoption_decision(
    config: dict[str, Any], rescue_decision: str | Path,
    selection_lock: str | Path, outcome_decision: str | Path,
    seed_confirmation: str | Path, order_invariant: str | Path,
    held_study: str | Path, label_benchmark: str | Path,
    output: str | Path,
) -> dict[str, Any]:
    """Verify all early-stage evidence and write the terminal adoption decision."""
    config_sha256 = config["_config_sha256"]
    selection_path = Path(selection_lock).resolve()
    selection = load_retargeted_bridge_selection_lock(selection_path, config)
    program_firewall_sha256 = validate_program_firewall(config)

    rescue_path, rescue = _read(rescue_decision)
    if (
        rescue.get("schema_version") != "masld-cl-control-only-rescue-result-v2"
        or rescue.get("config_sha256") != config_sha256
        or not rescue.get("control_only")
        or not rescue.get("complete_grid")
    ):
        raise ContractError("replay-plus-EWC rescue artifact is incomplete or incompatible")

    artifacts = {}
    for name, value in {
        "outcome_decision": outcome_decision,
        "seed_confirmation": seed_confirmation,
        "order_invariant": order_invariant,
        "held_study": held_study,
        "label_benchmark": label_benchmark,
    }.items():
        path, artifact = _read(value)
        _check_v9_identity(artifact, path, config_sha256, selection, selection_path)
        artifacts[name] = (path, artifact)

    outcome = artifacts["outcome_decision"][1]
    seeds = artifacts["seed_confirmation"][1]
    order = artifacts["order_invariant"][1]
    held = artifacts["held_study"][1]
    labels = artifacts["label_benchmark"][1]
    selected_summary = selection["summaries"][
        str(selection["selected_global_reference_weight"])
    ]

    settings = rescue["settings"]
    best_reference_f1 = max(
        float(row["metrics"]["reference_macro_f1_change_ci_low"]) for row in settings
    )
    best_jaccard_loss = min(
        float(row["metrics"]["reference_neighborhood_jaccard_loss"]) for row in settings
    )
    best_harmony_improvement = max(
        float(row["metrics"]["alignment_improvement_vs_harmony"]) for row in settings
    )
    worst_held = _minimum_held_protocol_harmony(held)
    min_lineage_change = min(labels["major_lineage_f1_changes_vs_best_non_cl"].values())

    gates = [
        _gate(
            "replay_ewc_control_only_survivor", "method",
            len(rescue["eligible_lambdas"]), ">= 1", bool(rescue["eligible_lambdas"]),
            "No replay-plus-EWC setting passed the reference-retention gates.",
        ),
        _gate(
            "replay_ewc_reference_macro_f1", "preservation", best_reference_f1,
            f"> {-config['gates']['reference_macro_f1_margin']}",
            best_reference_f1 > -config["gates"]["reference_macro_f1_margin"],
            "Best lower confidence bound across the rescue grid.",
        ),
        _gate(
            "replay_ewc_reference_neighborhood", "preservation", best_jaccard_loss,
            f"<= {config['gates']['reference_neighborhood_jaccard_loss']}",
            best_jaccard_loss <= config["gates"]["reference_neighborhood_jaccard_loss"],
            "Lowest reference k=30 neighborhood-Jaccard loss across the rescue grid.",
        ),
        _gate(
            "replay_ewc_alignment_vs_harmony", "alignment", best_harmony_improvement,
            f">= {config['gates']['minimum_control_alignment_improvement']}",
            best_harmony_improvement >= config["gates"]["minimum_control_alignment_improvement"],
            "Best Harmony improvement across the rescue grid.",
        ),
        _gate(
            "v9_is_replay_plus_ewc", "method", "retargeted_orthogonal_bridge",
            "published replay-plus-EWC objective", False,
            "V9 is a closed-form geometry-preserving adapter and cannot be relabeled as continual learning.",
        ),
        _gate(
            "v9_control_only_alignment", "alignment", selected_summary,
            "selected control-only setting eligible", selected_summary["eligible"],
            "Control alignment passed before outcome access.",
        ),
        _gate(
            "v9_disease_preservation", "disease", outcome["disease_preservation_pass"],
            "true", outcome["disease_preservation_pass"],
            "Frozen outcome analysis used donor-level raw-count-PCA comparators.",
        ),
        _gate(
            "v9_seed_confirmation", "robustness", seeds["all_five_seeds_pass"],
            "all five seeds pass", seeds["all_five_seeds_pass"],
            "Seeds resampled the donor-balanced closed-form bridge fit.",
        ),
        _gate(
            "v9_acquisition_order", "robustness",
            {"passing_orders": order["passing_order_count"],
             "minimum_distance_spearman": order["minimum_order_distance_spearman"]},
            "6/6 orders; distance Spearman >= 0.90", order["order_gate_pass"],
            "All six order applications were coordinate-equivalent.",
        ),
        _gate(
            "v9_query_macro_f1", "preservation",
            labels["query_macro_f1_change_vs_best_non_cl"],
            f">= {-config['gates']['query_macro_f1_margin']}",
            labels["query_macro_f1_change_vs_best_non_cl"]
            >= -config["gates"]["query_macro_f1_margin"],
            f"Best non-continual comparator was {labels['best_non_cl_method']}.",
        ),
        _gate(
            "v9_major_lineage_f1", "preservation", min_lineage_change,
            f">= {-config['gates']['major_lineage_f1_margin']}",
            min_lineage_change >= -config["gates"]["major_lineage_f1_margin"],
            "Worst major-lineage change versus the best non-continual comparator.",
        ),
        _gate(
            "v9_leave_query_study_out", "robustness", worst_held,
            f"all held studies pass and no protocol stratum < {-config['gates']['maximum_stratum_worsening']}",
            held["leave_query_study_out_pass"],
            "The worst held-study protocol/Harmony stratum failed the prespecified bound.",
        ),
        _gate(
            "frozen_program_firewall", "program_firewall", program_firewall_sha256,
            "all configured program hashes match", True,
            "The 117-program registry and donor-score authorities were validated read-only.",
        ),
    ]

    sources = {
        "replay_ewc_rescue": {"path": str(rescue_path), "sha256": sha256_path(rescue_path)},
        "selection_lock": {"path": str(selection_path), "sha256": sha256_path(selection_path)},
    }
    for name, (path, _) in artifacts.items():
        sources[name] = {"path": str(path), "sha256": sha256_path(path)}

    decision = {
        "schema_version": "masld-cl-adoption-decision-v1",
        "config_sha256": config_sha256,
        "selection_lock_sha256": selection["lock_sha256"],
        "decision": "reject",
        "harmony_remains_primary": True,
        "claim_allowed": False,
        "failure_reasons": [
            "The replay-plus-EWC grid had no setting that preserved the reference.",
            "The V9 bridge is not replay-plus-EWC, failed query-label noninferiority, and failed one held-study protocol stratum.",
        ],
        "gates": gates,
        "all_required_gates_pass": all(gate["passed"] for gate in gates),
        "not_run_after_early_stop": [
            "joint production update over every analyzed non-reference cell",
            "secondary GSE174748/GSE136103 stress tests",
            "GSE189600 descriptive projection",
            "BI-guided replay sensitivity",
            "paper or canonical-atlas promotion",
        ],
        "method_identity": {
            "replay_ewc": "corrected scvi-tools 1.3.3 replay-plus-EWC",
            "v9": "retargeted orthogonal bridge; diagnostic fallback only",
        },
        "program_firewall_sha256": program_firewall_sha256,
        "sources": sources,
    }
    write_json_exclusive(output, decision)
    return decision
