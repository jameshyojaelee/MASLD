"""Configuration loading and immutable run-directory helpers."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class ConfigError(ValueError):
    """Raised when the versioned pipeline configuration is invalid."""


def canonical_json_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def sha256_json(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def load_config(path: str | Path) -> dict[str, Any]:
    path = Path(path).resolve()
    with path.open() as handle:
        config = json.load(handle)
    if config.get("schema_version") != "masld-cl-v1":
        raise ConfigError("config schema_version must be 'masld-cl-v1'")
    required = {
        "input", "expected_contract", "roles", "features", "architecture",
        "sampling", "screen", "bootstrap", "gates",
    }
    missing = sorted(required - config.keys())
    if missing:
        raise ConfigError(f"config missing required sections: {missing}")
    features = config["features"]
    if features.get("selection_batch_key") != "technical_batch":
        raise ConfigError("the only permitted conditioning key is technical_batch")
    for key in ("categorical_covariate_keys", "continuous_covariate_keys"):
        if features.get(key, []) not in ([], None):
            raise ConfigError(f"additional model covariates are forbidden: {key}")
    forbidden = set(config["input"].get("forbidden_conditioning_fields", []))
    if features["selection_batch_key"] in forbidden:
        raise ConfigError("the configured batch key is a forbidden biological covariate")
    config["_config_path"] = str(path)
    public = {k: v for k, v in config.items() if not k.startswith("_")}
    config["_config_sha256"] = sha256_json(public)
    return config


def repo_path(config: dict[str, Any], value: str | Path) -> Path:
    path = Path(value)
    if path.is_absolute():
        return path
    config_path = Path(config["_config_path"])
    # config is Analysis/SingleCell/continual_integration/config_v1.json
    root = config_path.parents[3]
    return (root / path).resolve()


def new_candidate_dir(base: str | Path) -> Path:
    """Create a new UTC-stamped directory without overwriting an existing run."""
    base = Path(base)
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H%M%SZ")
    candidate = base / stamp
    candidate.mkdir(parents=True, exist_ok=False)
    return candidate


def write_json_exclusive(path: str | Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")
