"""Control-only hyperparameter selection and immutable selection lock."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from .config import canonical_json_bytes, write_json_exclusive
from .contracts import ContractError, sha256_path
from .execution import verify_execution_record
from .firewall import assert_control_only_metric_names


class SelectionError(RuntimeError):
    pass


REQUIRED_METRICS = {
    "reference_macro_f1_change_ci_low",
    "reference_neighborhood_jaccard_loss",
    "shift_control",
    "shift_control_standard_error",
}


def read_control_metrics(path: str | Path) -> list[dict[str, str]]:
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    required_columns = {"setting_id", "ewc_lambda", "replay_fraction", "metric", "value"}
    if not rows or not required_columns.issubset(rows[0]):
        raise SelectionError(f"control metrics require columns {sorted(required_columns)}")
    assert_control_only_metric_names([row["metric"] for row in rows])
    for row in rows:
        if not math.isfinite(float(row["value"])):
            raise SelectionError("control metrics contain a non-finite value")
    return rows


def _setting_id(ewc_lambda: float, replay_fraction: float) -> str:
    return f"lambda_{ewc_lambda:g}__replay_{replay_fraction:g}"


def _validate_exact_grid(
    config: dict[str, Any], rows: list[dict[str, str]], *, best_lambda: float | None,
) -> None:
    expected = {
        _setting_id(float(value), 0.20): (float(value), 0.20)
        for value in config["screen"]["lambda_values"]
    }
    if best_lambda is not None:
        for replay in config["screen"]["replay_values"]:
            expected[_setting_id(float(best_lambda), float(replay))] = (
                float(best_lambda), float(replay)
            )
    expected_count = 6 if best_lambda is None else 9
    if len(expected) != expected_count:
        raise SelectionError(f"prespecified grid must contain {expected_count} unique settings")
    observed: dict[str, tuple[float, float]] = {}
    metric_count: dict[str, int] = {}
    for row in rows:
        setting = row["setting_id"]
        pair = (float(row["ewc_lambda"]), float(row["replay_fraction"]))
        if setting in observed and observed[setting] != pair:
            raise SelectionError(f"setting parameters vary within grid: {setting}")
        observed[setting] = pair
        metric_count[setting] = metric_count.get(setting, 0) + 1
        if row.get("model_kind", "six_model_aggregate") != "six_model_aggregate":
            raise SelectionError("selection requires the prespecified six-model aggregate")
    if observed != expected:
        raise SelectionError(
            f"control grid differs from the exact prespecified settings; "
            f"expected={sorted(expected)}, observed={sorted(observed)}"
        )
    if any(metric_count[x] != len(REQUIRED_METRICS) for x in expected):
        raise SelectionError("every grid setting must contain exactly four control metrics")


def _verify_aggregate_provenance(
    config: dict[str, Any], metric_path: str | Path, expected_settings: int,
) -> dict[str, Any]:
    path = Path(metric_path).resolve()
    details_path = path.with_suffix(path.suffix + ".details.json")
    if not details_path.is_file():
        raise SelectionError("aggregated control metrics lack a provenance sidecar")
    with details_path.open() as handle:
        details = json.load(handle)
    expected_models = sorted({"all_lineage", *config["lineages"]})
    if (
        details.get("schema_version") != "masld-cl-control-aggregate-v1"
        or details.get("config_sha256") != config["_config_sha256"]
        or details.get("metrics_realpath") != str(path)
        or details.get("metrics_sha256") != sha256_path(path)
        or details.get("expected_model_kinds") != expected_models
        or len(details.get("sources", [])) != expected_settings * len(expected_models)
    ):
        raise SelectionError("aggregated control provenance is incomplete or inconsistent")
    for source in details["sources"]:
        for path_key, hash_key in (
            ("metrics", "metrics_sha256"), ("details", "details_sha256")
        ):
            if sha256_path(source[path_key]) != source[hash_key]:
                raise SelectionError("a control-evaluation source changed after aggregation")
    records = details.get("execution_records", [])
    expected_records = expected_settings * len(expected_models) + len(expected_models)
    if len(records) != expected_records or len({row.get("path") for row in records}) != expected_records:
        raise SelectionError(
            f"control aggregate requires exactly {expected_records} unique GPU records"
        )
    pipeline_root = Path(config["_config_path"]).resolve().parent
    for source in records:
        if sha256_path(source["path"]) != source["sha256"]:
            raise SelectionError("a control GPU record changed after aggregation")
        verify_execution_record(
            source["path"], pipeline_root, config["_config_sha256"]
        )
    return {
        "aggregate_details_realpath": str(details_path.resolve()),
        "aggregate_details_sha256": sha256_path(details_path),
        "execution_records": records,
    }


def select_setting(config: dict[str, Any], rows: list[dict[str, str]]) -> dict[str, Any]:
    by_setting: dict[str, dict[str, float]] = {}
    parameters: dict[str, tuple[float, float]] = {}
    for row in rows:
        setting = row["setting_id"]
        if row["metric"] in by_setting.get(setting, {}):
            raise SelectionError(
                f"setting {setting} contains duplicate metric {row['metric']}; "
                "aggregate the six model kinds first"
            )
        by_setting.setdefault(setting, {})[row["metric"]] = float(row["value"])
        current = (float(row["ewc_lambda"]), float(row["replay_fraction"]))
        if setting in parameters and parameters[setting] != current:
            raise SelectionError(f"setting has inconsistent parameters: {setting}")
        parameters[setting] = current
    candidates = []
    for setting, metrics in by_setting.items():
        missing = REQUIRED_METRICS - metrics.keys()
        if missing:
            raise SelectionError(f"setting {setting} missing metrics: {sorted(missing)}")
        eligible = (
            parameters[setting][0] > 0
            and parameters[setting][1] > 0
            and
            metrics["reference_macro_f1_change_ci_low"]
            > -config["gates"]["reference_macro_f1_margin"]
            and metrics["reference_neighborhood_jaccard_loss"]
            <= config["gates"]["reference_neighborhood_jaccard_loss"]
        )
        if eligible:
            candidates.append((setting, metrics, parameters[setting]))
    if not candidates:
        raise SelectionError("no positive-EWC, positive-replay setting passes reference-retention gates")
    best_shift = min(x[1]["shift_control"] for x in candidates)
    best_se = next(
        x[1]["shift_control_standard_error"]
        for x in candidates if x[1]["shift_control"] == best_shift
    )
    one_se = [x for x in candidates if x[1]["shift_control"] <= best_shift + best_se]
    one_se.sort(key=lambda x: (x[2][0], x[2][1], x[1]["shift_control"], x[0]))
    selected = one_se[0]
    runner_pool = sorted(
        (x for x in candidates if x[0] != selected[0]),
        key=lambda x: (x[1]["shift_control"], x[2][0], x[2][1], x[0]),
    )
    return {
        "selected": {
            "setting_id": selected[0],
            "ewc_lambda": selected[2][0],
            "replay_fraction": selected[2][1],
            "control_metrics": selected[1],
        },
        "pareto_runner_up": None if not runner_pool else {
            "setting_id": runner_pool[0][0],
            "ewc_lambda": runner_pool[0][2][0],
            "replay_fraction": runner_pool[0][2][1],
            "control_metrics": runner_pool[0][1],
        },
        "eligible_setting_ids": sorted(x[0] for x in candidates),
        "one_standard_error_threshold": best_shift + best_se,
    }


def select_provisional_lambda(config: dict[str, Any], rows: list[dict[str, str]]) -> dict[str, Any]:
    """Choose the lambda used to instantiate the replay screen without freezing outcomes."""
    _validate_exact_grid(config, rows, best_lambda=None)
    selection = select_setting(config, rows)
    selected = selection["selected"]
    if float(selected["replay_fraction"]) != 0.20:
        raise SelectionError("provisional lambda input must contain only replay=0.20 settings")
    return {
        "schema_version": "masld-cl-provisional-lambda-v1",
        "config_sha256": config["_config_sha256"],
        "outcomes_unlocked": False,
        "best_lambda": float(selected["ewc_lambda"]),
        "setting_id": selected["setting_id"],
        "selection_basis": "reference_retention_and_query_control_shift_only",
    }


def write_provisional_lambda_lock(
    config: dict[str, Any], metric_path: str | Path, output: str | Path,
) -> dict[str, Any]:
    rows = read_control_metrics(metric_path)
    result = select_provisional_lambda(config, rows)
    provenance = _verify_aggregate_provenance(config, metric_path, expected_settings=6)
    result.update({
        "control_metrics_realpath": str(Path(metric_path).resolve()),
        "control_metrics_sha256": sha256_path(metric_path),
        **provenance,
    })
    result["lock_sha256"] = hashlib.sha256(canonical_json_bytes(result)).hexdigest()
    write_json_exclusive(output, result)
    return result


def load_provisional_lambda_lock(
    config: dict[str, Any], path: str | Path,
) -> dict[str, Any]:
    with Path(path).open() as handle:
        lock = json.load(handle)
    payload = {key: value for key, value in lock.items() if key != "lock_sha256"}
    if (
        lock.get("schema_version") != "masld-cl-provisional-lambda-v1"
        or lock.get("config_sha256") != config["_config_sha256"]
        or hashlib.sha256(canonical_json_bytes(payload)).hexdigest()
        != lock.get("lock_sha256")
        or sha256_path(lock.get("control_metrics_realpath", ""))
        != lock.get("control_metrics_sha256")
        or sha256_path(lock.get("aggregate_details_realpath", ""))
        != lock.get("aggregate_details_sha256")
    ):
        raise SelectionError("provisional lambda lock is invalid or its sources changed")
    records = lock.get("execution_records", [])
    expected_records = 6 * (len(config["lineages"]) + 1) + len(config["lineages"]) + 1
    if len(records) != expected_records:
        raise SelectionError("provisional lambda lock lacks its exact GPU record roster")
    pipeline_root = Path(config["_config_path"]).resolve().parent
    for source in records:
        if sha256_path(source["path"]) != source["sha256"]:
            raise SelectionError("provisional lambda GPU record changed")
        verify_execution_record(
            source["path"], pipeline_root, config["_config_sha256"]
        )
    return lock


def write_selection_lock(
    config: dict[str, Any], metric_path: str | Path,
    provisional_lambda_lock: str | Path, output: str | Path,
) -> dict[str, Any]:
    rows = read_control_metrics(metric_path)
    provisional = load_provisional_lambda_lock(config, provisional_lambda_lock)
    best_lambda = float(provisional["best_lambda"])
    _validate_exact_grid(config, rows, best_lambda=best_lambda)
    provenance = _verify_aggregate_provenance(config, metric_path, expected_settings=9)
    selection = select_setting(config, rows)
    metric_hash = hashlib.sha256(Path(metric_path).read_bytes()).hexdigest()
    lock = {
        "schema_version": "masld-cl-selection-v1",
        "config_sha256": config["_config_sha256"],
        "control_metrics_sha256": metric_hash,
        "control_metrics_realpath": str(Path(metric_path).resolve()),
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "selection_frozen": True,
        "provisional_lambda_lock_realpath": str(Path(provisional_lambda_lock).resolve()),
        "provisional_lambda_lock_sha256": sha256_path(provisional_lambda_lock),
        **provenance,
        **selection,
    }
    lock["lock_sha256"] = hashlib.sha256(canonical_json_bytes(lock)).hexdigest()
    write_json_exclusive(output, lock)
    return lock
