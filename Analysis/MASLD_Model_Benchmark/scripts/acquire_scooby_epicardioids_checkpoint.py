#!/usr/bin/env python3
"""Verify and freeze the exact Scooby Epicardioids safetensors payload."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from typing import Any, Mapping

from masld_bench.artifacts import ArtifactError, verify_frozen_tree, write_json_exclusive
from scripts.scooby_remote_safetensors_header import inspect_reader


SCHEMA = "masld-bench-scooby-epicardioids-checkpoint-acquisition-v1"


class ScoobyCheckpointAcquisitionError(RuntimeError):
    """Raised when the downloaded payload differs from its frozen requirements."""


class LocalRangeReader:
    """Read only the bounded byte ranges requested by the header validator."""

    def __init__(self, path: Path) -> None:
        self.path = path
        self.size = path.stat().st_size
        self.request_count = 0
        self.bytes_downloaded = 0

    def read(self, start: int, length: int) -> bytes:
        if start < 0 or length <= 0 or start + length > self.size:
            raise ScoobyCheckpointAcquisitionError("local range request is outside the checkpoint")
        with self.path.open("rb") as handle:
            handle.seek(start)
            payload = handle.read(length)
        require(len(payload) == length, "local range read was not exact")
        self.request_count += 1
        self.bytes_downloaded += len(payload)
        return payload


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ScoobyCheckpointAcquisitionError(message)


def digest(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def load_json(path: Path, *, label: str) -> dict[str, Any]:
    require(path.is_file() and not path.is_symlink(), f"{label} is not a regular file")
    try:
        value = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as error:
        raise ScoobyCheckpointAcquisitionError(f"invalid {label}: {error}") from error
    require(isinstance(value, dict), f"{label} is not an object")
    return value


def verify_authority(root: Path, binding: Mapping[str, Any], *, label: str) -> Path:
    path = root / str(binding["path"])
    try:
        verify_frozen_tree(path)
    except ArtifactError as error:
        raise ScoobyCheckpointAcquisitionError(f"{label} artifact differs: {error}") from error
    require(digest(path / "ARTIFACTS.json") == binding["artifacts_sha256"], f"{label} identity differs")
    return path


def validate_task_boundary(config: Mapping[str, Any]) -> None:
    task = config["task_scope"]
    policy = config["gate_policy"]
    require(task["observed_atac_in_context"], "observed-ATAC provenance was lost")
    require(not task["rna_conditioned_atac_eligible"], "observed-ATAC model cannot enter RNA-conditioned lane")
    require(not task["liver_query_encoder_available"], "unreleased liver query encoder was asserted")
    require(not task["full_sequence_forward_allowed"], "acquisition cannot execute a full sequence forward")
    require(not task["sealed_inference_allowed"], "acquisition cannot authorize sealed inference")
    require(not task["universal_or_champion_claim_allowed"], "acquisition cannot authorize a champion claim")
    require(not policy["checkpoint_deserialization_allowed"], "acquisition cannot deserialize weights")
    require(not policy["model_forward_allowed"], "acquisition cannot execute a model forward")
    require(not policy["project_data_read_allowed"] and not policy["sealed_data_read_allowed"], "data firewall differs")


def compare_header(observed: Mapping[str, Any], frozen: Mapping[str, Any]) -> None:
    keys = (
        "file_size_bytes",
        "header_size_bytes",
        "tensor_payload_size_bytes",
        "tensor_count",
        "parameter_count",
        "dtypes",
        "metadata",
        "tensors",
    )
    require(all(observed.get(key) == frozen.get(key) for key in keys), "local safetensors header differs from frozen bounded header")


def inspect_local_header(checkpoint: Path, expected_parameters: int) -> dict[str, Any]:
    require(checkpoint.is_file() and not checkpoint.is_symlink(), "checkpoint is not a regular file")
    return inspect_reader(LocalRangeReader(checkpoint), expected_parameters)


def audit(root: Path, config_path: Path, checkpoint: Path, output: Path) -> dict[str, Any]:
    root = root.resolve(strict=True)
    config = load_json(config_path.resolve(strict=True), label="acquisition config")
    require(config.get("schema_version") == SCHEMA, "acquisition config identity differs")
    require(not output.exists() and not output.is_symlink(), "refusing to overwrite audit output")
    validate_task_boundary(config)
    schema_artifact = verify_authority(root, config["authorities"]["runtime_schema_artifact"], label="runtime schema")
    verify_authority(root, config["authorities"]["bounded_header_artifact"], label="bounded header")
    schema_receipt = load_json(schema_artifact / "audit/receipt.json", label="runtime schema receipt")
    require(schema_receipt["status"] == "pass_checkpoint_download_and_strict_restore_probe_allowed", "runtime schema did not authorize acquisition")
    for name in ("bounded_header", "checkpoint_contract"):
        binding = config["authorities"][name]
        require(digest(root / binding["path"]) == binding["sha256"], f"{name} authority differs")

    contract = config["checkpoint"]
    checkpoint = checkpoint.resolve(strict=True)
    require(checkpoint.stat().st_size == contract["size_bytes"], "checkpoint byte size differs")
    observed_sha256 = digest(checkpoint)
    require(observed_sha256 == contract["sha256"], "checkpoint SHA-256 differs")
    observed_header = inspect_local_header(checkpoint, contract["state_dict_elements_including_buffers"])
    frozen_header = load_json(root / config["authorities"]["bounded_header"]["path"], label="bounded header")
    compare_header(observed_header, frozen_header)
    require(observed_header["header_size_bytes"] == contract["header_size_bytes"], "checkpoint header size differs")
    require(observed_header["tensor_count"] == contract["state_dict_tensors"], "checkpoint tensor count differs")

    output.mkdir(parents=True, mode=0o750)
    receipt = {
        "schema_version": "masld-bench-scooby-epicardioids-checkpoint-acquisition-receipt-v1",
        "status": "pass_exact_checkpoint_acquired_safe_restore_allowed",
        "config_sha256": digest(config_path),
        "checkpoint": {
            "repository": contract["repository"],
            "revision": contract["revision"],
            "filename": contract["filename"],
            "format": contract["format"],
            "sha256": observed_sha256,
            "size_bytes": checkpoint.stat().st_size,
            "header_size_bytes": observed_header["header_size_bytes"],
            "state_dict_tensors": observed_header["tensor_count"],
            "state_dict_elements_including_buffers": observed_header["parameter_count"],
            "full_payload_hash_verified": True,
        },
        "license": config["license"],
        "task_scope": config["task_scope"],
        "checkpoint_payload_downloaded": True,
        "checkpoint_deserialized": False,
        "model_forward_executed": False,
        "full_sequence_forward_executed": False,
        "project_data_read": False,
        "sealed_data_read": False,
        "next_gate": config["gate_policy"]["next_if_pass"],
    }
    write_json_exclusive(output / "receipt.json", receipt, mode=0o640)
    return receipt


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--checkpoint", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.root, args.config, args.checkpoint, args.output), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
