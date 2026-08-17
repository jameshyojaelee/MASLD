"""Fail-closed disease-preservation decision for the selected V9 bridge."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .firewall import load_retargeted_bridge_selection_lock, validate_program_firewall


def decide_orthogonal_bridge_outcomes(
    config: dict[str, Any], selection_lock: str | Path,
    metric_values: list[str | Path], output: str | Path,
) -> dict[str, Any]:
    selection_path = Path(selection_lock).resolve()
    selection = load_retargeted_bridge_selection_lock(selection_path, config)
    validate_program_firewall(config)
    expected = {"all_lineage", *config["lineages"]}
    if len(metric_values) != len(expected):
        raise ContractError("V9 outcome decision requires the exact six-model roster")
    thresholds = {
        "pooled_disease_retention": float(config["gates"]["pooled_disease_retention"]),
        "per_study_disease_retention": float(config["gates"]["per_study_disease_retention"]),
        "within_study_distance_spearman": float(config["gates"]["within_study_distance_spearman"]),
    }
    models, sources = {}, []
    for value in metric_values:
        path = Path(value).resolve()
        details_path = path.with_suffix(path.suffix + ".details.json")
        with details_path.open() as handle:
            details = json.load(handle)
        kind = details.get("model_kind")
        if (
            details.get("schema_version") != "masld-cl-orthogonal-bridge-matched-outcome-v9"
            or details.get("config_sha256") != config["_config_sha256"]
            or details.get("selection_lock_sha256") != selection["lock_sha256"]
            or details.get("metrics_sha256") != sha256_path(path)
            or float(details.get("selected_global_reference_weight", -1))
            != float(selection["selected_global_reference_weight"])
            or kind not in expected or kind in models
        ):
            raise ContractError("invalid or duplicate V9 matched outcome result")
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        gates = []
        for row in rows:
            metric = row["metric"]
            if metric not in thresholds:
                raise ContractError(f"unexpected V9 outcome metric: {metric}")
            observed, threshold = float(row["value"]), thresholds[metric]
            gates.append({
                "metric": metric, "scope": row["scope"], "threshold": threshold,
                "comparison": ">=", "observed": observed, "pass": bool(observed >= threshold),
            })
        if not gates:
            raise ContractError("V9 outcome result has no testable gates")
        models[kind] = {"pass": bool(all(row["pass"] for row in gates)), "gates": gates}
        sources.append({
            "metrics": str(path), "metrics_sha256": sha256_path(path),
            "details": str(details_path.resolve()), "details_sha256": sha256_path(details_path),
        })
    if set(models) != expected:
        raise ContractError("V9 outcome decision model roster is incomplete")
    passed = bool(all(row["pass"] for row in models.values()))
    decision = {
        "schema_version": "masld-cl-orthogonal-bridge-outcome-decision-v9",
        "config_sha256": config["_config_sha256"],
        "selection_lock": str(selection_path), "selection_lock_file_sha256": sha256_path(selection_path),
        "selection_lock_sha256": selection["lock_sha256"],
        "selected_control_offset_weight": selection["selected_control_offset_weight"],
        "selected_global_reference_weight": selection["selected_global_reference_weight"],
        "thresholds": thresholds, "models": models, "sources": sources,
        "disease_preservation_pass": passed, "promotion_eligible": False,
        "decision": "continue_remaining_gates" if passed else "reject",
        "failure_consequence": "Harmony remains primary if this or any remaining gate fails.",
        "language_constraint": "cross-sectional stage-associated remodeling",
    }
    write_json_exclusive(output, decision)
    return decision
