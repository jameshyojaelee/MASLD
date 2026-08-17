"""Freeze median stopping epochs from verified selected-setting GPU runs."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

from .config import canonical_json_bytes, write_json_exclusive
from .contracts import sha256_path
from .execution import verify_execution_record
from .firewall import load_selection_lock


class RefitError(RuntimeError):
    pass


def write_refit_lock(
    config: dict[str, Any], selection_lock: str | Path,
    execution_records: list[str | Path], output: str | Path,
) -> dict[str, Any]:
    selection = load_selection_lock(selection_lock, config)
    expected_count = (len(config["lineages"]) + 1) * len(
        config["screen"]["confirmation_seeds"]
    )
    resolved_records = [str(Path(value).resolve()) for value in execution_records]
    if len(resolved_records) != expected_count or len(set(resolved_records)) != expected_count:
        raise RefitError(f"refit lock requires exactly {expected_count} unique GPU records")
    pipeline_root = Path(config["_config_path"]).resolve().parent
    selected = selection["selected"]
    rows = []
    record_sources = []
    for value in resolved_records:
        path = Path(value)
        record = verify_execution_record(
            path, pipeline_root, config["_config_sha256"]
        )
        identity = record["result_identity"]
        root = Path(identity["output_realpath"])
        manifest_paths = [
            root / item["path"] for item in identity["output_files"]
            if item["path"].endswith("update_manifest.json")
        ]
        matched = []
        for manifest_path in manifest_paths:
            with manifest_path.open() as handle:
                manifest = json.load(handle)
            if (
                manifest.get("schema_version") == "masld-cl-update-v1"
                and manifest.get("method") == "continual_learning"
                and manifest.get("replay_mode") == "random"
                and manifest.get("sensitivity_only") is False
                and manifest.get("production") is False
                and manifest.get("selection_lock_sha256") == selection["lock_sha256"]
                and float(manifest.get("ewc_lambda")) == float(selected["ewc_lambda"])
                and float(manifest.get("replay_fraction")) == float(selected["replay_fraction"])
                and not manifest.get("held_out_datasets")
                and set(manifest.get("query_datasets", []))
                == set(config["evaluation"]["powered_query_studies"])
            ):
                matched.append(manifest)
        if len(matched) != 1:
            raise RefitError(
                f"each refit GPU record must own exactly one selected development fit: {path}"
            )
        rows.extend(matched)
        record_sources.append({"path": str(path), "sha256": sha256_path(path)})
    expected_seeds = set(config["screen"]["confirmation_seeds"])
    expected_models = {"all_lineage", *config["lineages"]}
    observed_pairs = [(row.get("model_kind"), row.get("seed")) for row in rows]
    expected_pairs = {(model, seed) for model in expected_models for seed in expected_seeds}
    if set(observed_pairs) != expected_pairs or len(observed_pairs) != len(expected_pairs):
        raise RefitError("verified runs do not contain the exact 30 selected-setting model/seed fits")
    epochs = {}
    for model_kind in sorted(expected_models):
        values = np.asarray([
            int(row["stopping_epoch"]) for row in rows if row["model_kind"] == model_kind
        ], dtype=int)
        if np.any(values <= 0):
            raise RefitError("stopping epochs must be positive")
        epochs[f"continual_learning|{model_kind}"] = int(np.median(values))
    lock = {
        "schema_version": "masld-cl-refit-v2",
        "config_sha256": config["_config_sha256"],
        "selection_lock_sha256": selection["lock_sha256"],
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "median_stopping_epochs": epochs,
        "execution_records": record_sources,
    }
    lock["lock_sha256"] = hashlib.sha256(canonical_json_bytes(lock)).hexdigest()
    write_json_exclusive(output, lock)
    return lock


def load_refit_lock(
    path: str | Path, config: dict[str, Any], selection: dict[str, Any]
) -> dict[str, Any]:
    with Path(path).open() as handle:
        lock = json.load(handle)
    if lock.get("schema_version") != "masld-cl-refit-v2":
        raise RefitError("unsupported refit lock")
    payload = {key: value for key, value in lock.items() if key != "lock_sha256"}
    if hashlib.sha256(canonical_json_bytes(payload)).hexdigest() != lock.get("lock_sha256"):
        raise RefitError("refit lock content hash mismatch")
    if lock.get("config_sha256") != config["_config_sha256"]:
        raise RefitError("refit lock config hash mismatch")
    if lock.get("selection_lock_sha256") != selection["lock_sha256"]:
        raise RefitError("refit lock selection hash mismatch")
    for source in lock.get("execution_records", []):
        if sha256_path(source["path"]) != source["sha256"]:
            raise RefitError("refit execution record changed after locking")
        verify_execution_record(
            source["path"], Path(config["_config_path"]).resolve().parent,
            config["_config_sha256"],
        )
    expected_count = (len(config["lineages"]) + 1) * len(
        config["screen"]["confirmation_seeds"]
    )
    if len(lock.get("execution_records", [])) != expected_count:
        raise RefitError("refit lock lacks the exact selected-run roster")
    expected = {
        f"continual_learning|{model_kind}" for model_kind in ["all_lineage", *config["lineages"]]
    }
    if set(lock.get("median_stopping_epochs", {})) != expected:
        raise RefitError("refit lock lacks the exact six selected model kinds")
    return lock
