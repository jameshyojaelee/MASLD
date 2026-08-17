"""Selection-lock and frozen-program firewalls."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

from .config import canonical_json_bytes, repo_path
from .contracts import sha256_path


class FirewallError(RuntimeError):
    pass


FORBIDDEN_BEFORE_SELECTION = {
    "case", "disease", "stage", "fibrosis", "program", "hero", "cas13", "perturbation"
}


def sha256_file(path: str | Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(chunk_size), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_program_firewall(config: dict[str, Any]) -> dict[str, str]:
    observed = {}
    for relpath, expected in config.get("program_firewall", {}).items():
        path = repo_path(config, relpath)
        value = sha256_file(path)
        observed[relpath] = value
        if value != expected:
            raise FirewallError(f"frozen program artifact changed: {relpath}")
    if not observed:
        raise FirewallError("program firewall is empty")
    return observed


def assert_control_only_metric_names(names: list[str]) -> None:
    bad = sorted(
        name for name in names
        if any(token in name.lower() for token in FORBIDDEN_BEFORE_SELECTION)
    )
    if bad:
        raise FirewallError(f"outcome-bearing metrics were exposed before selection: {bad}")


def load_selection_lock(path: str | Path, config: dict[str, Any]) -> dict[str, Any]:
    with Path(path).open() as handle:
        lock = json.load(handle)
    if lock.get("schema_version") != "masld-cl-selection-v1":
        raise FirewallError("unsupported selection lock")
    if lock.get("config_sha256") != config["_config_sha256"]:
        raise FirewallError("selection lock config hash mismatch")
    payload = {k: v for k, v in lock.items() if k != "lock_sha256"}
    observed = hashlib.sha256(canonical_json_bytes(payload)).hexdigest()
    if observed != lock.get("lock_sha256"):
        raise FirewallError("selection lock content hash mismatch")
    if not lock.get("selection_frozen"):
        raise FirewallError("selection lock is not frozen")
    for path_key, hash_key in (
        ("control_metrics_realpath", "control_metrics_sha256"),
        ("aggregate_details_realpath", "aggregate_details_sha256"),
        ("provisional_lambda_lock_realpath", "provisional_lambda_lock_sha256"),
    ):
        if sha256_path(lock.get(path_key, "")) != lock.get(hash_key):
            raise FirewallError(f"selection source changed: {path_key}")
    records = lock.get("execution_records", [])
    expected = 9 * (len(config["lineages"]) + 1) + len(config["lineages"]) + 1
    if len(records) != expected or len({row.get("path") for row in records}) != expected:
        raise FirewallError("selection lock lacks the exact 60 GPU record roster")
    from .execution import verify_execution_record
    pipeline_root = Path(config["_config_path"]).resolve().parent
    for source in records:
        if sha256_path(source["path"]) != source["sha256"]:
            raise FirewallError("selection GPU record changed after locking")
        verify_execution_record(
            source["path"], pipeline_root, config["_config_sha256"]
        )
    return lock


def load_control_adapter_selection_lock(
    path: str | Path, config: dict[str, Any]
) -> dict[str, Any]:
    """Verify the immutable control-only v7 selection without requiring current source identity."""
    with Path(path).open() as handle:
        lock = json.load(handle)
    if lock.get("schema_version") != "masld-cl-control-adapter-selection-v7":
        raise FirewallError("unsupported control-adapter selection lock")
    if lock.get("config_sha256") != config["_config_sha256"]:
        raise FirewallError("control-adapter selection config hash mismatch")
    payload = {k: v for k, v in lock.items() if k != "lock_sha256"}
    observed = hashlib.sha256(canonical_json_bytes(payload)).hexdigest()
    if observed != lock.get("lock_sha256"):
        raise FirewallError("control-adapter selection content hash mismatch")
    if lock.get("selection_frozen") is not True or lock.get("outcomes_unlocked") is not True:
        raise FirewallError("control-adapter outcome firewall is not unlocked")
    if sha256_path(lock.get("policy_realpath", "")) != lock.get("policy_sha256"):
        raise FirewallError("control-adapter selection policy changed")
    sources = lock.get("control_results", [])
    expected = 3 * (len(config["lineages"]) + 1)
    if len(sources) != expected or len({row.get("path") for row in sources}) != expected:
        raise FirewallError("control-adapter lock lacks the exact 18-result grid")
    keys = set()
    frozen_source_identity = lock.get("source_identity")
    for source in sources:
        if sha256_path(source.get("path", "")) != source.get("sha256"):
            raise FirewallError("control-adapter result changed after locking")
        with Path(source["path"]).open() as handle:
            result = json.load(handle)
        if (
            result.get("schema_version") != "masld-cl-control-adapter-pilot-v7"
            or result.get("config_sha256") != config["_config_sha256"]
            or result.get("control_only") is not True
            or result.get("outcomes_unlocked") is not False
        ):
            raise FirewallError("invalid control-adapter result in selection lock")
        manifest_path = result.get("sources", {}).get("adapter_manifest", "")
        if sha256_path(manifest_path) != result.get("sources", {}).get("adapter_manifest_sha256"):
            raise FirewallError("control-adapter manifest changed after locking")
        with Path(manifest_path).open() as handle:
            manifest = json.load(handle)
        if manifest.get("adapter_source_identity") != frozen_source_identity:
            raise FirewallError("control-adapter result has a different frozen source identity")
        keys.add((float(result["global_reference_weight"]), result["model_kind"]))
    expected_keys = {
        (weight, kind)
        for weight in (0.25, 0.5, 0.75)
        for kind in ("all_lineage", *config["lineages"])
    }
    if keys != expected_keys:
        raise FirewallError("control-adapter result grid is incomplete")
    if float(lock.get("selected_global_reference_weight", -1)) not in (0.25, 0.5, 0.75):
        raise FirewallError("control-adapter selected weight is invalid")
    return lock


def load_retargeted_bridge_selection_lock(
    path: str | Path, config: dict[str, Any]
) -> dict[str, Any]:
    """Verify the immutable V9 control-only selection and its complete grid."""
    from .execution import verify_source_identity_payload

    with Path(path).open() as handle:
        lock = json.load(handle)
    if lock.get("schema_version") != "masld-cl-orthogonal-bridge-retarget-selection-v9":
        raise FirewallError("unsupported retargeted-bridge selection lock")
    if lock.get("config_sha256") != config["_config_sha256"]:
        raise FirewallError("retargeted-bridge selection config hash mismatch")
    payload = {key: value for key, value in lock.items() if key != "lock_sha256"}
    observed = hashlib.sha256(canonical_json_bytes(payload)).hexdigest()
    if observed != lock.get("lock_sha256"):
        raise FirewallError("retargeted-bridge selection content hash mismatch")
    if lock.get("selection_frozen") is not True or lock.get("outcomes_unlocked") is not True:
        raise FirewallError("retargeted-bridge outcome firewall is not unlocked")
    if sha256_path(lock.get("policy_realpath", "")) != lock.get("policy_sha256"):
        raise FirewallError("retargeted-bridge selection policy changed")
    frozen_source = verify_source_identity_payload(lock.get("source_identity"))
    sources = lock.get("control_results", [])
    grid = (0.0, 0.25, 0.5, 0.75, 1.0)
    expected_count = len(grid) * (len(config["lineages"]) + 1)
    if len(sources) != expected_count or len({row.get("path") for row in sources}) != expected_count:
        raise FirewallError("retargeted-bridge lock lacks the exact 30-result grid")
    keys = set()
    for source in sources:
        if sha256_path(source.get("path", "")) != source.get("sha256"):
            raise FirewallError("retargeted-bridge result changed after locking")
        with Path(source["path"]).open() as handle:
            result = json.load(handle)
        if (
            result.get("schema_version") != "masld-cl-orthogonal-bridge-retarget-pilot-v9"
            or result.get("config_sha256") != config["_config_sha256"]
            or result.get("control_only") is not True
            or result.get("outcomes_unlocked") is not False
        ):
            raise FirewallError("invalid retargeted-bridge result in selection lock")
        manifest_path = result.get("sources", {}).get("bridge_manifest", "")
        if sha256_path(manifest_path) != result.get("sources", {}).get("bridge_manifest_sha256"):
            raise FirewallError("retargeted-bridge manifest changed after locking")
        with Path(manifest_path).open() as handle:
            manifest = json.load(handle)
        if (
            manifest.get("schema_version") != "masld-cl-orthogonal-bridge-retarget-v9"
            or verify_source_identity_payload(manifest.get("source_identity")) != frozen_source
        ):
            raise FirewallError("retargeted-bridge result has a different frozen source identity")
        keys.add((float(result["global_reference_weight"]), result["model_kind"]))
    expected_keys = {
        (weight, kind)
        for weight in grid
        for kind in ("all_lineage", *config["lineages"])
    }
    if keys != expected_keys:
        raise FirewallError("retargeted-bridge result grid is incomplete")
    if float(lock.get("selected_control_offset_weight", -1)) != 1.0:
        raise FirewallError("retargeted-bridge control offset differs from the frozen setting")
    if float(lock.get("selected_global_reference_weight", -1)) not in grid:
        raise FirewallError("retargeted-bridge selected target weight is invalid")
    return lock


def load_outcome_selection_lock(path: str | Path, config: dict[str, Any]) -> dict[str, Any]:
    """Dispatch to the matching immutable selection-lock verifier."""
    with Path(path).open() as handle:
        schema = json.load(handle).get("schema_version")
    if schema == "masld-cl-selection-v1":
        return load_selection_lock(path, config)
    if schema == "masld-cl-control-adapter-selection-v7":
        return load_control_adapter_selection_lock(path, config)
    if schema == "masld-cl-orthogonal-bridge-retarget-selection-v9":
        return load_retargeted_bridge_selection_lock(path, config)
    raise FirewallError("unsupported outcome selection lock")
