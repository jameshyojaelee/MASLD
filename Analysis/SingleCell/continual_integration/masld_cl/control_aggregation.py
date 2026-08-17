"""Aggregate six model-specific control evaluations for one global selection."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np

from .config import write_json_exclusive
from .contracts import ContractError, sha256_path
from .execution import require_execution_ownership, verify_execution_record
from .firewall import assert_control_only_metric_names


EXPECTED_METRICS = {
    "reference_macro_f1_change_ci_low", "reference_neighborhood_jaccard_loss",
    "shift_control", "shift_control_standard_error",
}


def aggregate_control_metrics(
    config: dict[str, Any], inputs: list[str | Path], output: str | Path,
) -> list[dict[str, Any]]:
    expected_models = {"all_lineage", *config["lineages"]}
    records: dict[str, dict[str, dict[str, float]]] = {}
    parameters: dict[str, tuple[float, float]] = {}
    sources = []
    execution_records: dict[str, dict[str, str]] = {}
    pipeline_root = Path(config["_config_path"]).resolve().parent
    for value in inputs:
        path = Path(value).resolve()
        details_path = path.with_suffix(path.suffix + ".details.json")
        if not details_path.is_file():
            raise ContractError(f"control metric source lacks provenance sidecar: {path}")
        with details_path.open() as handle:
            details = json.load(handle)
        if (
            details.get("schema_version") != "masld-cl-control-evaluation-v1"
            or details.get("config_sha256") != config["_config_sha256"]
            or details.get("metrics_realpath") != str(path)
            or details.get("metrics_sha256") != sha256_path(path)
            or details.get("candidate_seed") != config["screen"]["seed"]
        ):
            raise ContractError(f"control metric provenance mismatch: {path}")
        for prefix in ("reference_run_manifest", "candidate_run_manifest"):
            artifact = Path(details[prefix]).resolve()
            if sha256_path(artifact) != details[f"{prefix}_sha256"]:
                raise ContractError(f"control metric training source changed: {artifact}")
        for role in ("reference", "candidate"):
            record_path = Path(details[f"{role}_execution_record"]).resolve()
            if sha256_path(record_path) != details[f"{role}_execution_record_sha256"]:
                raise ContractError(f"control {role} execution record changed")
            record = verify_execution_record(
                record_path, pipeline_root, config["_config_sha256"]
            )
            require_execution_ownership(
                record,
                [details[f"{role}_run_manifest"], details[f"{role}_embedding"]],
                role=f"control {role}",
            )
            execution_records[str(record_path)] = {
                "path": str(record_path), "sha256": sha256_path(record_path)
            }
        sources.append({
            "metrics": str(path), "metrics_sha256": sha256_path(path),
            "details": str(details_path.resolve()),
            "details_sha256": sha256_path(details_path),
        })
        with path.open(newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        required = {"setting_id", "ewc_lambda", "replay_fraction", "model_kind", "metric", "value"}
        if not rows or not required.issubset(rows[0]):
            raise ContractError(f"model-specific control metrics are malformed: {value}")
        assert_control_only_metric_names([row["metric"] for row in rows])
        for row in rows:
            setting, model, metric = row["setting_id"], row["model_kind"], row["metric"]
            if model not in expected_models or metric not in EXPECTED_METRICS:
                raise ContractError(f"unexpected control metric scope: {model}|{metric}")
            model_values = records.setdefault(setting, {}).setdefault(model, {})
            if metric in model_values:
                raise ContractError(f"duplicate model-specific control metric: {setting}|{model}|{metric}")
            model_values[metric] = float(row["value"])
            pair = (float(row["ewc_lambda"]), float(row["replay_fraction"]))
            if setting in parameters and parameters[setting] != pair:
                raise ContractError(f"setting parameters vary across model kinds: {setting}")
            parameters[setting] = pair
    output_rows = []
    for setting, models in sorted(records.items()):
        if set(models) != expected_models:
            raise ContractError(f"setting lacks all six model kinds: {setting}")
        if any(set(values) != EXPECTED_METRICS for values in models.values()):
            raise ContractError(f"setting has an incomplete control metric family: {setting}")
        aggregate = {
            "reference_macro_f1_change_ci_low": min(x["reference_macro_f1_change_ci_low"] for x in models.values()),
            "reference_neighborhood_jaccard_loss": max(x["reference_neighborhood_jaccard_loss"] for x in models.values()),
            "shift_control": float(np.mean([x["shift_control"] for x in models.values()])),
            "shift_control_standard_error": max(x["shift_control_standard_error"] for x in models.values()),
        }
        ewc_lambda, replay = parameters[setting]
        for metric, value in aggregate.items():
            output_rows.append({
                "setting_id": setting, "ewc_lambda": ewc_lambda,
                "replay_fraction": replay, "model_kind": "six_model_aggregate",
                "metric": metric, "scope": "prespecified_six_model_aggregate", "value": value,
            })
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    write_json_exclusive(output.with_suffix(output.suffix + ".details.json"), {
        "schema_version": "masld-cl-control-aggregate-v1",
        "config_sha256": config["_config_sha256"],
        "metrics_realpath": str(output.resolve()),
        "metrics_sha256": sha256_path(output),
        "expected_model_kinds": sorted(expected_models),
        "sources": sources,
        "execution_records": [execution_records[key] for key in sorted(execution_records)],
    })
    return output_rows
