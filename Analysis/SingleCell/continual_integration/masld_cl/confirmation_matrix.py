"""Lock the selected, runner-up, and paper-setting five-seed GPU matrix."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from .config import canonical_json_bytes, write_json_exclusive
from .contracts import ContractError, sha256_path
from .execution import verify_execution_record
from .firewall import load_selection_lock


def confirmation_settings(selection: dict[str, Any]) -> dict[tuple[float, float], list[str]]:
    by_pair: dict[tuple[float, float], list[str]] = {}

    def add(role: str, value: dict[str, Any]) -> None:
        pair = (float(value["ewc_lambda"]), float(value["replay_fraction"]))
        if pair[0] <= 0 or pair[1] <= 0:
            raise ContractError(f"confirmation role is not replay-plus-EWC: {role}")
        by_pair.setdefault(pair, []).append(role)

    add("selected", selection["selected"])
    if selection.get("pareto_runner_up") is not None:
        add("pareto_runner_up", selection["pareto_runner_up"])
    add("paper", {"ewc_lambda": 100.0, "replay_fraction": 0.20})
    return {pair: sorted(roles) for pair, roles in sorted(by_pair.items())}


def _expected_keys(
    config: dict[str, Any], settings: dict[tuple[float, float], list[str]],
) -> set[tuple[float, float, int, str]]:
    return {
        (pair[0], pair[1], int(seed), model_kind)
        for pair in settings
        for seed in config["screen"]["confirmation_seeds"]
        for model_kind in ["all_lineage", *config["lineages"]]
    }


def _manifest_key(
    manifest: dict[str, Any], config: dict[str, Any], selection: dict[str, Any],
) -> tuple[float, float, int, str]:
    if (
        manifest.get("schema_version") != "masld-cl-update-v1"
        or manifest.get("config_sha256") != config["_config_sha256"]
        or manifest.get("selection_lock_sha256") != selection["lock_sha256"]
        or manifest.get("method") != "continual_learning"
        or manifest.get("replay_mode") != "random"
        or manifest.get("sensitivity_only") is not False
        or manifest.get("production") is not False
        or manifest.get("held_out_datasets")
        or set(manifest.get("query_datasets", []))
        != set(config["evaluation"]["powered_query_studies"])
        or int(manifest.get("stopping_epoch", 0)) <= 0
    ):
        raise ContractError("confirmation execution is not a qualifying development fit")
    return (
        float(manifest["ewc_lambda"]), float(manifest["replay_fraction"]),
        int(manifest["seed"]), str(manifest["model_kind"]),
    )


def write_confirmation_matrix_lock(
    config: dict[str, Any], selection_lock: str | Path,
    execution_records: list[str | Path], output: str | Path,
) -> dict[str, Any]:
    selection = load_selection_lock(selection_lock, config)
    settings = confirmation_settings(selection)
    expected = _expected_keys(config, settings)
    paths = [str(Path(value).resolve()) for value in execution_records]
    if len(paths) != len(set(paths)) or len(paths) != len(expected):
        raise ContractError(
            f"confirmation matrix requires exactly {len(expected)} unique GPU records"
        )
    pipeline_root = Path(config["_config_path"]).resolve().parent
    observed: dict[tuple[float, float, int, str], dict[str, Any]] = {}
    sources = []
    for path_value in sorted(paths):
        record = verify_execution_record(
            path_value, pipeline_root, config["_config_sha256"]
        )
        manifest_path = Path(record["result_identity"]["result_manifest_realpath"])
        with manifest_path.open() as handle:
            manifest = json.load(handle)
        key = _manifest_key(manifest, config, selection)
        if key not in expected or key in observed:
            raise ContractError(f"unexpected or duplicated confirmation fit: {key}")
        observed[key] = manifest
        sources.append({
            "path": path_value, "sha256": sha256_path(path_value),
            "result_manifest": str(manifest_path.resolve()),
            "result_manifest_sha256": sha256_path(manifest_path),
            "ewc_lambda": key[0], "replay_fraction": key[1],
            "seed": key[2], "model_kind": key[3],
            "roles": settings[(key[0], key[1])],
        })
    if set(observed) != expected:
        raise ContractError("confirmation GPU records do not contain the exact required matrix")
    lock = {
        "schema_version": "masld-cl-confirmation-matrix-v1",
        "config_sha256": config["_config_sha256"],
        "selection_lock_sha256": selection["lock_sha256"],
        "settings": [
            {"ewc_lambda": pair[0], "replay_fraction": pair[1], "roles": roles}
            for pair, roles in settings.items()
        ],
        "expected_fit_count": len(expected),
        "execution_records": sources,
    }
    lock["lock_sha256"] = hashlib.sha256(canonical_json_bytes(lock)).hexdigest()
    write_json_exclusive(output, lock)
    return lock


def load_confirmation_matrix_lock(
    path: str | Path, config: dict[str, Any], selection: dict[str, Any],
) -> dict[str, Any]:
    with Path(path).open() as handle:
        lock = json.load(handle)
    payload = {key: value for key, value in lock.items() if key != "lock_sha256"}
    if (
        lock.get("schema_version") != "masld-cl-confirmation-matrix-v1"
        or lock.get("config_sha256") != config["_config_sha256"]
        or lock.get("selection_lock_sha256") != selection["lock_sha256"]
        or hashlib.sha256(canonical_json_bytes(payload)).hexdigest()
        != lock.get("lock_sha256")
    ):
        raise ContractError("confirmation matrix lock is invalid")
    settings = confirmation_settings(selection)
    expected = _expected_keys(config, settings)
    sources = lock.get("execution_records", [])
    if lock.get("expected_fit_count") != len(expected) or len(sources) != len(expected):
        raise ContractError("confirmation matrix lock has the wrong fit count")
    pipeline_root = Path(config["_config_path"]).resolve().parent
    observed = set()
    for source in sources:
        if (
            sha256_path(source["path"]) != source["sha256"]
            or sha256_path(source["result_manifest"])
            != source["result_manifest_sha256"]
        ):
            raise ContractError("confirmation matrix source changed")
        record = verify_execution_record(
            source["path"], pipeline_root, config["_config_sha256"]
        )
        if (
            str(Path(record["result_identity"]["result_manifest_realpath"]).resolve())
            != str(Path(source["result_manifest"]).resolve())
        ):
            raise ContractError("confirmation record owns a different result manifest")
        with Path(source["result_manifest"]).open() as handle:
            key = _manifest_key(json.load(handle), config, selection)
        if key in observed or key not in expected:
            raise ContractError("confirmation matrix source is duplicated or unexpected")
        if source.get("roles") != settings[(key[0], key[1])]:
            raise ContractError("confirmation matrix role mapping changed")
        observed.add(key)
    if observed != expected:
        raise ContractError("confirmation matrix lock lacks required fits")
    return lock
