#!/usr/bin/env python3
"""Verify the frozen validation bundle required by GSE281364 production v2."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import re
from typing import Any

from masld_bench.artifacts import verify_frozen_tree


SCHEMA = "masld-bench-gse281364-open-sequence-production-binding-v2"
VALIDATION_ARTIFACT_CLASS = (
    "gse281364_open_sequence_head_campaign_validation_v2"
)
MISSING_BASELINES = (
    "deltaSVM",
    "gkm-SVM",
    "sequence_CNN",
    "sequence_transformer",
    "MPRALegNet_signed_MPRA_rekeyed_to_1033_elements_239_blocks",
)
REQUIRED_SOURCES = (
    "config/gse281364_open_sequence_head_campaign.toml",
    "scripts/freeze_gse281364_open_sequence_taskspec.py",
    "scripts/fit_gse281364_open_sequence_heads.py",
    "scripts/evaluate_gse281364_open_sequence_heads.py",
    "scripts/verify_gse281364_open_sequence_production_binding.py",
    "tests/unit/test_gse281364_open_sequence_head_campaign.py",
    "tests/unit/test_verify_gse281364_open_sequence_production_binding.py",
    "slurm/run_gse281364_open_sequence_head_campaign_v2_cpu.sbatch",
    "slurm/validate_gse281364_open_sequence_head_campaign_v2_cpu.sbatch",
)
PRODUCTION_WRAPPER = (
    "slurm/run_gse281364_open_sequence_head_campaign_v2_cpu.sbatch"
)
CONFIG = "config/gse281364_open_sequence_head_campaign.toml"
SHA256 = re.compile(r"^[0-9a-f]{64}$")
PRODUCTION_VALIDATION_NAME = re.compile(r"^model-check-224-[0-9]+$")


class BindingError(RuntimeError):
    """Raised when production is not bound to its exact validation bundle."""


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError) as error:
        raise BindingError(f"invalid {label}: {error}") from error
    if not isinstance(value, dict):
        raise BindingError(f"{label} must be a JSON object")
    return value


def load_source_manifest(path: Path) -> dict[str, str]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeDecodeError) as error:
        raise BindingError(f"cannot read validation source manifest: {error}") from error
    records: dict[str, str] = {}
    for line in lines:
        match = re.fullmatch(r"([0-9a-f]{64})  ([^\n]+)", line)
        if match is None:
            raise BindingError("validation source-manifest syntax differs")
        digest, relative_value = match.groups()
        relative = Path(relative_value)
        if (
            relative.is_absolute()
            or ".." in relative.parts
            or relative_value in records
        ):
            raise BindingError("validation source-manifest path differs")
        records[relative_value] = digest
    if set(records) != set(REQUIRED_SOURCES):
        raise BindingError("validation source-manifest roster differs")
    return records


def verify_binding(
    *,
    root: Path,
    validation_root: Path,
    validation_sha256: str,
    require_production_location: bool = True,
) -> dict[str, Any]:
    if not SHA256.fullmatch(validation_sha256):
        raise BindingError("VALIDATION_SHA is not a lowercase SHA-256")
    root = root.resolve(strict=True)
    executions = (root / "executions").resolve(strict=True)
    if not validation_root.is_absolute() or validation_root.is_symlink():
        raise BindingError("VALIDATION_ROOT must be an absolute non-symlink path")
    try:
        validation = validation_root.resolve(strict=True)
        validation.relative_to(executions)
    except (OSError, ValueError) as error:
        raise BindingError("VALIDATION_ROOT is outside the benchmark executions") from error
    if require_production_location and (
        validation.parent != executions
        or PRODUCTION_VALIDATION_NAME.fullmatch(validation.name) is None
    ):
        raise BindingError("VALIDATION_ROOT is not the frozen v2 validation attempt")
    artifacts_path = validation / "ARTIFACTS.json"
    if artifacts_path.is_symlink() or file_sha256(artifacts_path) != validation_sha256:
        raise BindingError("VALIDATION_SHA does not bind validation ARTIFACTS.json")
    verify_frozen_tree(validation)
    artifacts = load_json(artifacts_path, label="validation ARTIFACTS")
    metadata = artifacts.get("metadata", {})
    receipt = load_json(validation / "validation/receipt.json", label="validation receipt")
    source_manifest = validation / "source.sha256"
    records = load_source_manifest(source_manifest)
    source_manifest_sha256 = file_sha256(source_manifest)
    if (
        metadata.get("artifact_class") != VALIDATION_ARTIFACT_CLASS
        or metadata.get("status") != "passed_bound_production_wrapper"
        or metadata.get("campaign_scope") != "partial_common_head_campaign"
        or metadata.get("mandatory_baselines_complete") is not False
        or tuple(metadata.get("mandatory_baselines_missing", ()))
        != MISSING_BASELINES
        or metadata.get("shortlist_blocked") is not True
        or metadata.get("finalist_claim_blocked") is not True
        or metadata.get("complementarity_blocked") is not True
        or metadata.get("conditional_trigger_blocked") is not True
        or metadata.get("outcomes_read") is not False
        or metadata.get("model_fit") is not False
        or metadata.get("predictions_generated") is not False
        or metadata.get("production_wrapper_submitted") is not False
        or metadata.get("source_manifest_sha256") != source_manifest_sha256
        or receipt.get("schema_version")
        != "masld-bench-gse281364-open-sequence-validation-v2"
        or receipt.get("status") != "passed_bound_production_wrapper"
        or tuple(receipt.get("mandatory_baselines_missing", ()))
        != MISSING_BASELINES
        or receipt.get("production_wrapper_submitted") is not False
        or receipt.get("outcomes_read") is not False
        or receipt.get("model_fit") is not False
        or receipt.get("predictions_generated") is not False
        or receipt.get("source_manifest_sha256") != source_manifest_sha256
    ):
        raise BindingError("v2 validation metadata or firewall differs")
    for relative_value in REQUIRED_SOURCES:
        path = root / relative_value
        if path.is_symlink() or not path.is_file():
            raise BindingError(f"validated source is unavailable: {relative_value}")
        if file_sha256(path) != records[relative_value]:
            raise BindingError(f"validated source drifted: {relative_value}")
    wrapper_sha256 = file_sha256(root / PRODUCTION_WRAPPER)
    config_sha256 = file_sha256(root / CONFIG)
    if (
        metadata.get("production_wrapper_sha256") != wrapper_sha256
        or receipt.get("production_wrapper_sha256") != wrapper_sha256
        or records[PRODUCTION_WRAPPER] != wrapper_sha256
        or metadata.get("config_sha256") != config_sha256
        or receipt.get("config_sha256") != config_sha256
        or records[CONFIG] != config_sha256
    ):
        raise BindingError("v2 production wrapper or config binding differs")
    result = {
        "schema_version": SCHEMA,
        "status": "pass_exact_frozen_validation_binding",
        "validation_root": validation.as_posix(),
        "validation_artifacts_sha256": validation_sha256,
        "source_manifest_sha256": source_manifest_sha256,
        "validated_source_count": len(records),
        "production_wrapper_sha256": wrapper_sha256,
        "config_sha256": config_sha256,
        "campaign_scope": "partial_common_head_campaign",
        "mandatory_baselines_complete": False,
        "mandatory_baselines_missing": list(MISSING_BASELINES),
        "shortlist_blocked": True,
        "finalist_claim_blocked": True,
        "complementarity_blocked": True,
        "conditional_trigger_blocked": True,
        "outcomes_read": False,
    }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--validation-root", type=Path, required=True)
    parser.add_argument("--validation-sha", required=True)
    arguments = parser.parse_args()
    receipt = verify_binding(
        root=arguments.root,
        validation_root=arguments.validation_root,
        validation_sha256=arguments.validation_sha,
    )
    print(json.dumps(receipt, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
