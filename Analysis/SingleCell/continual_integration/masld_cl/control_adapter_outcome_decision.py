"""Fail-closed disease-preservation decision for the selected v7 adapter."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .firewall import load_control_adapter_selection_lock, validate_program_firewall


def decide_control_adapter_outcomes(
    config: dict[str, Any], selection_lock: str | Path,
    metric_values: list[str | Path], output: str | Path,
) -> dict[str, Any]:
    selection_path = Path(selection_lock).resolve()
    selection = load_control_adapter_selection_lock(selection_path, config)
    validate_program_firewall(config)
    expected = {"all_lineage", *config["lineages"]}
    if len(metric_values) != len(expected):
        raise ContractError("outcome decision requires the exact six-model roster")
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
            details.get("schema_version") != "masld-cl-control-adapter-outcome-v7"
            or details.get("config_sha256") != config["_config_sha256"]
            or details.get("selection_lock_sha256") != selection["lock_sha256"]
            or details.get("metrics_sha256") != sha256_path(path)
            or details.get("global_reference_weight")
            != selection["selected_global_reference_weight"]
            or kind not in expected or kind in models
        ):
            raise ContractError("invalid or duplicate control-adapter outcome result")
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        gates = []
        for row in rows:
            metric = row["metric"]
            if metric not in thresholds:
                raise ContractError(f"unexpected outcome metric: {metric}")
            observed = float(row["value"])
            threshold = thresholds[metric]
            gates.append({
                "metric": metric, "scope": row["scope"], "threshold": threshold,
                "comparison": ">=", "observed": observed,
                "pass": bool(observed >= threshold),
            })
        if not gates:
            raise ContractError("outcome result has no testable gates")
        models[kind] = {
            "pass": bool(all(row["pass"] for row in gates)),
            "gates": gates,
        }
        sources.append({
            "metrics": str(path), "metrics_sha256": sha256_path(path),
            "details": str(details_path.resolve()),
            "details_sha256": sha256_path(details_path),
        })
    if set(models) != expected:
        raise ContractError("outcome decision model roster is incomplete")
    decision = {
        "schema_version": "masld-cl-control-adapter-outcome-decision-v7",
        "config_sha256": config["_config_sha256"],
        "selection_lock": str(selection_path),
        "selection_lock_file_sha256": sha256_path(selection_path),
        "selection_lock_sha256": selection["lock_sha256"],
        "selected_global_reference_weight": selection["selected_global_reference_weight"],
        "thresholds": thresholds,
        "models": models,
        "sources": sources,
        "disease_preservation_pass": bool(all(row["pass"] for row in models.values())),
        "promotion_eligible": False,
        "decision": "reject" if not all(row["pass"] for row in models.values()) else "continue_remaining_gates",
        "failure_consequence": "Harmony remains primary; do not run confirmatory seeds or robustness jobs.",
        "language_constraint": "cross-sectional stage-associated remodeling",
    }
    write_json_exclusive(output, decision)
    return decision
