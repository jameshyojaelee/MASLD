#!/usr/bin/env python3
"""Freeze a source-complete wrapper around the validated specialist check."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any

from masld_bench.artifacts import freeze_tree, reject_symlink_components, verify_frozen_tree, write_json_exclusive
from scripts.validate_gse296875_observed_multiome_specialist_gate_v2 import validate as validate_revision2


SCHEMA = "masld-bench-observed-multiome-specialist-gate-readiness-v3"


class SpecialistGateValidationError(ValueError):
    """Raised when the source-complete wrapper differs."""


def _digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        while block := handle.read(8 * 1024 * 1024):
            value.update(block)
    return value.hexdigest()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise SpecialistGateValidationError("JSON object required")
    return value


def _bound_file(root: Path, record: dict[str, Any], label: str) -> Path:
    path = reject_symlink_components(root / str(record.get("path", "")), label=label).resolve(strict=True)
    path.relative_to(root)
    if _digest(path) != record.get("sha256"):
        raise SpecialistGateValidationError(f"{label} drifted")
    return path


def _bound_tree(root: Path, record: dict[str, Any], label: str) -> Path:
    path = reject_symlink_components(root / str(record.get("path", "")), label=label).resolve(strict=True)
    path.relative_to(root)
    verify_frozen_tree(path)
    if _digest(path / "ARTIFACTS.json") != record.get("artifacts_sha256"):
        raise SpecialistGateValidationError(f"{label} drifted")
    return path


def validate(root: Path, config: dict[str, Any]) -> dict[str, Any]:
    if config.get("schema_version") != SCHEMA or config.get("dataset_id") != "gse296875" or config.get("stage") != "development" or config.get("revision") != 3:
        raise SpecialistGateValidationError("wrapper identity differs")
    incident_path = _bound_file(root, config.get("incident", {}), "revision-2 incident")
    if _json(incident_path).get("terminal_disposition") != "superseded_incomplete_source_lock":
        raise SpecialistGateValidationError("revision-2 incident disposition differs")
    sources = config.get("source_bindings")
    expected_sources = {"base_module", "stats_module", "gate_module", "revision2_validator", "wrapper", "unit_test", "sbatch"}
    if not isinstance(sources, dict) or set(sources) != expected_sources:
        raise SpecialistGateValidationError("transitive source binding roster differs")
    for label, record in sources.items():
        _bound_file(root, record, label)

    revision2_config = _bound_file(root, config.get("revision2_config", {}), "revision-2 config")
    recalculated = validate_revision2(root, _json(revision2_config))
    revision2_artifact = _bound_tree(root, config.get("revision2_artifact", {}), "revision-2 artifact")
    frozen_receipt = _json(revision2_artifact / "receipt.json")
    if recalculated != frozen_receipt or frozen_receipt.get("gate_ready") is not True:
        raise SpecialistGateValidationError("revision-2 receipt does not reproduce exactly")
    return {
        **frozen_receipt,
        "schema_version": "masld-bench-observed-multiome-specialist-gate-readiness-receipt-v3",
        "revision": 3,
        "supersedes_incomplete_source_lock_job_id": 21100612,
        "revision2_artifacts_sha256": config["revision2_artifact"]["artifacts_sha256"],
        "transitive_source_bindings_complete": True,
        "revision2_receipt_exactly_reproduced": True,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    config_path = args.config.resolve(strict=True)
    config_path.relative_to(root)
    output = reject_symlink_components(args.output, label="specialist gate source-complete output")
    output.resolve(strict=False).parent.relative_to(root / "executions")
    output.mkdir(parents=True, exist_ok=False)
    receipt = validate(root, _json(config_path))
    write_json_exclusive(output / "receipt.json", receipt)
    sources = [config_path, Path(__file__).resolve(strict=True)] + [root / record["path"] for record in _json(config_path)["source_bindings"].values()]
    (output / "source.sha256").write_text("".join(f"{_digest(path)}  {path}\n" for path in sources), encoding="utf-8")
    digest = freeze_tree(output, metadata={"artifact_class": "gse296875_observed_multiome_specialist_gate_readiness", "revision": 3, "source_complete": True, "outcome_blind": True, "sealed_outcomes_accessed": False, "supports_external_claim": False, "supports_champion_claim": False})
    print(json.dumps({"output": str(output), "artifacts_sha256": digest}, sort_keys=True))


if __name__ == "__main__":
    main()
