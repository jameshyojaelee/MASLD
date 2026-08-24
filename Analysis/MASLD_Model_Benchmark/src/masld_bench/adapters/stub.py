"""Deterministic fixture adapter for end-to-end control-plane tests only."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
from typing import Any, Mapping

from ..artifacts import write_json_exclusive
from .base import AdapterAction, AdapterError


def _prediction(run_id: str, row_id: str) -> float:
    value = int.from_bytes(sha256(f"{run_id}\0{row_id}".encode()).digest()[:8], "big")
    return value / float((1 << 64) - 1)


def run_fixture(
    action: AdapterAction,
    request: Mapping[str, Any],
    output_dir: str | Path,
) -> Path:
    """Create protocol-valid artifacts; never calculate evaluation metrics."""

    output = Path(output_dir)
    if output.exists():
        raise AdapterError(f"fixture adapter output already exists: {output}")
    output.mkdir(parents=True)
    run_id = str(request.get("run_id", ""))
    if len(run_id) != 64 or any(
        character not in "0123456789abcdef" for character in run_id
    ):
        raise AdapterError("fixture request must contain a full lowercase SHA-256 run_id")
    run_spec = request.get("run_spec")
    if not isinstance(run_spec, Mapping):
        raise AdapterError("fixture request must contain a RunSpec object")
    inputs = run_spec.get("inputs")
    if not isinstance(inputs, list):
        raise AdapterError("fixture RunSpec inputs must be a list")
    environment_inputs = [
        item
        for item in inputs
        if isinstance(item, Mapping)
        and str(item.get("role", "")).startswith("environment:")
    ]
    if len(environment_inputs) != 1:
        raise AdapterError("fixture RunSpec must bind exactly one environment artifact")
    environment_sha256 = str(environment_inputs[0].get("sha256", ""))
    runtime_id = str(run_spec.get("runtime_id", ""))
    fit_dataset_ids = request.get("fit_dataset_ids")
    if not isinstance(fit_dataset_ids, list) or any(
        not isinstance(item, str) for item in fit_dataset_ids
    ):
        raise AdapterError(
            "fixture request must bind the frozen fit_dataset_ids array"
        )
    payload: dict[str, Any]
    artifact_name: str
    if action == AdapterAction.PREDICT:
        row_ids = request.get("row_ids", [])
        if not isinstance(row_ids, list) or not row_ids:
            raise AdapterError("fixture prediction requires a non-empty row_ids list")
        payload = {
            "schema_version": "masld-bench-fixture-predictions-v1",
            "run_id": run_id,
            "rows": [
                {"row_id": str(row_id), "prediction": _prediction(run_id, str(row_id))}
                for row_id in row_ids
            ],
        }
        artifact_name = "predictions.json"
    else:
        payload = {
            "schema_version": "masld-bench-fixture-artifact-v1",
            "run_id": run_id,
            "action": action.value,
            "fixture_only": True,
        }
        artifact_name = f"{action.value}.json"
    record = write_json_exclusive(output / artifact_name, payload)
    receipt = {
        "schema_version": "masld-bench-adapter-receipt-v1",
        "action": action.value,
        "run_id": run_id,
        "status": "complete",
        "artifacts": [
            {
                "path": artifact_name,
                "sha256": record.sha256,
                "size_bytes": record.size_bytes,
            }
        ],
        "metadata": {
            "adapter": "deterministic_fixture",
            "not_for_scientific_use": True,
            "environment_artifact_sha256": environment_sha256,
            "runtime_id": runtime_id,
            "fit_dataset_ids": list(fit_dataset_ids),
        },
    }
    write_json_exclusive(output / "adapter_receipt.json", receipt)
    return output


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--action", required=True, choices=[item.value for item in AdapterAction])
    parser.add_argument("--request", required=True)
    parser.add_argument("--output", required=True)
    arguments = parser.parse_args(argv)
    request = json.loads(Path(arguments.request).read_text(encoding="utf-8"))
    run_fixture(AdapterAction(arguments.action), request, arguments.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
