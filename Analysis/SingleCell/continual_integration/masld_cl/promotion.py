"""Prespecified main-figure promotion decision from long-form metrics."""

from __future__ import annotations

import csv
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable

from .config import write_json_exclusive
from .firewall import load_selection_lock, validate_program_firewall


@dataclass(frozen=True)
class GateResult:
    gate_id: str
    category: str
    passed: bool
    observed: Any
    threshold: str
    reason: str


class PromotionError(RuntimeError):
    pass


def read_metrics(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required = {"metric", "scope", "value"}
    if not rows or not required.issubset(rows[0]):
        raise PromotionError(f"evaluation metrics require columns {sorted(required)}")
    return rows


def _values(rows: list[dict[str, str]], metric: str) -> list[float]:
    found = [float(row["value"]) for row in rows if row["metric"] == metric]
    if not found:
        raise PromotionError(f"required promotion metric is missing: {metric}")
    return found


def _gate(
    rows: list[dict[str, str]], metric: str, category: str, reducer: Callable[[list[float]], float],
    predicate: Callable[[float], bool], threshold: str,
) -> GateResult:
    try:
        observed = reducer(_values(rows, metric))
        passed = bool(predicate(observed))
        reason = "pass" if passed else f"observed {observed:g} does not satisfy {threshold}"
    except PromotionError as exc:
        observed = None
        passed = False
        reason = str(exc)
    return GateResult(metric, category, passed, observed, threshold, reason)


def _decision_from_gates(results: list[GateResult]) -> str:
    if all(result.passed for result in results):
        return "promote_main"
    hard_categories = {
        "integrity", "preservation", "disease", "robustness", "program_firewall",
    }
    hard_pass = all(
        result.passed for result in results if result.category in hard_categories
    )
    noninferior_alignment = next(
        (
            result.passed for result in results
            if result.gate_id == "maximum_protocol_stratum_worsening"
        ),
        False,
    )
    return "supplement_only" if hard_pass and noninferior_alignment else "reject"


def evaluate_promotion(
    config: dict[str, Any], rows: list[dict[str, str]], selection_lock: str | Path
) -> dict[str, Any]:
    selection = load_selection_lock(selection_lock, config)
    firewall = validate_program_firewall(config)
    gates = config["gates"]
    results = [
        _gate(rows, "contract_valid", "integrity", min, lambda x: x == 1, "all == 1"),
        _gate(rows, "runs_complete", "integrity", min, lambda x: x == 1, "all == 1"),
        _gate(rows, "input_hashes_valid", "integrity", min, lambda x: x == 1, "all == 1"),
        _gate(rows, "reference_macro_f1_change_ci_low", "preservation", min,
              lambda x: x > -gates["reference_macro_f1_margin"],
              f"> {-gates['reference_macro_f1_margin']}"),
        _gate(rows, "reference_neighborhood_jaccard_loss", "preservation", max,
              lambda x: x <= gates["reference_neighborhood_jaccard_loss"],
              f"<= {gates['reference_neighborhood_jaccard_loss']}"),
        _gate(rows, "query_macro_f1_change", "preservation", min,
              lambda x: x >= -gates["query_macro_f1_margin"],
              f">= {-gates['query_macro_f1_margin']}"),
        _gate(rows, "major_lineage_f1_change", "preservation", min,
              lambda x: x >= -gates["major_lineage_f1_margin"],
              f">= {-gates['major_lineage_f1_margin']}"),
        _gate(rows, "alignment_improvement_vs_harmony", "alignment", min,
              lambda x: x >= gates["minimum_control_alignment_improvement"],
              f">= {gates['minimum_control_alignment_improvement']}"),
        _gate(rows, "alignment_improvement_vs_harmony_ci_low", "alignment", min,
              lambda x: x > 0, "> 0"),
        _gate(rows, "alignment_improvement_vs_architecture_surgery", "alignment", min,
              lambda x: x >= gates["minimum_control_alignment_improvement"],
              f">= {gates['minimum_control_alignment_improvement']}"),
        _gate(rows, "alignment_improvement_vs_architecture_surgery_ci_low", "alignment", min,
              lambda x: x > 0, "> 0"),
        _gate(rows, "improved_lineage_count", "alignment", min,
              lambda x: x >= gates["minimum_improved_lineages"],
              f">= {gates['minimum_improved_lineages']}"),
        _gate(rows, "all_lineage_alignment_pass", "alignment", min, lambda x: x == 1, "== 1"),
        _gate(rows, "maximum_protocol_stratum_worsening", "alignment", max,
              lambda x: x <= gates["maximum_stratum_worsening"],
              f"<= {gates['maximum_stratum_worsening']}"),
        _gate(rows, "pooled_disease_retention", "disease", min,
              lambda x: x >= gates["pooled_disease_retention"],
              f">= {gates['pooled_disease_retention']}"),
        _gate(rows, "per_study_disease_retention", "disease", min,
              lambda x: x >= gates["per_study_disease_retention"],
              f">= {gates['per_study_disease_retention']}"),
        _gate(rows, "within_study_distance_spearman", "disease", min,
              lambda x: x >= gates["within_study_distance_spearman"],
              f">= {gates['within_study_distance_spearman']}"),
        _gate(rows, "all_orders_reference_and_disease_pass", "robustness", min,
              lambda x: x == 1, "== 1"),
        _gate(rows, "passing_order_count", "robustness", min,
              lambda x: x >= gates["minimum_passing_orders"],
              f">= {gates['minimum_passing_orders']}"),
        _gate(rows, "worst_order_degradation", "robustness", max,
              lambda x: x <= gates["maximum_order_degradation"],
              f"<= {gates['maximum_order_degradation']}"),
        _gate(rows, "order_distance_spearman", "robustness", min,
              lambda x: x >= gates["order_distance_spearman"],
              f">= {gates['order_distance_spearman']}"),
        _gate(rows, "order_to_seed_variability_ratio", "robustness", max,
              lambda x: x <= gates["maximum_order_to_seed_variability_ratio"],
              f"<= {gates['maximum_order_to_seed_variability_ratio']}"),
        _gate(rows, "selected_seed_reference_and_disease_pass", "robustness", min,
              lambda x: x == 1, "== 1"),
        _gate(rows, "selected_seed_alignment_pass_count", "robustness", min,
              lambda x: x == len(config["screen"]["confirmation_seeds"]),
              f"== {len(config['screen']['confirmation_seeds'])}"),
        _gate(rows, "held_study_reference_and_disease_pass", "robustness", min,
              lambda x: x == 1, "== 1"),
        _gate(rows, "held_study_alignment_pass", "robustness", min,
              lambda x: x == 1, "== 1"),
        _gate(rows, "secondary_stress_pass", "robustness", min,
              lambda x: x == 1, "== 1"),
        _gate(rows, "gpu_tolerance_measured", "integrity", min,
              lambda x: x == 1, "== 1"),
        _gate(rows, "descriptive_projection_only", "integrity", min,
              lambda x: x == 1, "== 1"),
        _gate(rows, "frozen_program_scores_unchanged", "program_firewall", min,
              lambda x: x == 1, "== 1"),
        _gate(rows, "weighted_bh_family_size", "program_firewall", min,
              lambda x: x == config["program_inference_contract"]["weighted_bh_family_size"],
              f"== {config['program_inference_contract']['weighted_bh_family_size']}"),
        _gate(rows, "unweighted_bh_family_size", "program_firewall", min,
              lambda x: x == config["program_inference_contract"]["unweighted_bh_family_size"],
              f"== {config['program_inference_contract']['unweighted_bh_family_size']}"),
    ]
    decision = _decision_from_gates(results)
    return {
        "schema_version": "masld-cl-promotion-v1",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "config_sha256": config["_config_sha256"],
        "selection_lock_sha256": selection["lock_sha256"],
        "selected_setting": selection["selected"],
        "program_firewall_sha256": firewall,
        "decision": decision,
        "all_gates_pass": all(x.passed for x in results),
        "gates": [asdict(x) for x in results],
        "claim_allowed": decision == "promote_main",
        "claim_text": (
            "Continual expansion improved alignment of protocol-compatible control donors "
            "while preserving donor-level MASLD-associated variation under held-study, seed, "
            "and acquisition-order perturbations."
            if decision == "promote_main" else None
        ),
    }


def write_promotion_decision(
    config: dict[str, Any], metrics_path: str | Path, metrics_lock: str | Path,
    selection_lock: str | Path, output: str | Path
) -> dict[str, Any]:
    from .metric_bundle import verify_metric_bundle
    selection = load_selection_lock(selection_lock, config)
    verified = verify_metric_bundle(config, metrics_path, metrics_lock, selection)
    decision = evaluate_promotion(config, read_metrics(metrics_path), selection_lock)
    decision["metric_bundle_lock_sha256"] = verified["lock_sha256"]
    write_json_exclusive(output, decision)
    return decision
