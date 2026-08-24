#!/usr/bin/env python3
"""Safely inventory a tensor-only Corgi checkpoint without importing Corgi."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import torch


def _hash_file(path: Path, algorithm: str) -> str:
    digest = hashlib.new(algorithm)
    with path.open("rb") as handle:
        while chunk := handle.read(16 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _tensor_mapping_candidates(
    value: Any,
    path: tuple[str, ...] = (),
) -> list[tuple[tuple[str, ...], Mapping[str, torch.Tensor]]]:
    candidates: list[tuple[tuple[str, ...], Mapping[str, torch.Tensor]]] = []
    if not isinstance(value, Mapping):
        return candidates
    if value and all(
        isinstance(key, str) and isinstance(item, torch.Tensor)
        for key, item in value.items()
    ):
        candidates.append((path, value))
    for key, item in value.items():
        if isinstance(key, str) and isinstance(item, Mapping):
            candidates.extend(_tensor_mapping_candidates(item, path + (key,)))
    return candidates


def inventory_checkpoint(
    path: Path,
    *,
    expected_size: int,
    expected_md5: str,
    expected_sha256: str,
) -> dict[str, object]:
    path = path.resolve(strict=True)
    if not path.is_file() or path.is_symlink():
        raise ValueError("checkpoint must be a regular non-symlink file")
    if path.stat().st_size != expected_size:
        raise ValueError("checkpoint size differs")
    if _hash_file(path, "md5") != expected_md5:
        raise ValueError("checkpoint MD5 differs")
    if _hash_file(path, "sha256") != expected_sha256:
        raise ValueError("checkpoint SHA-256 differs")

    loaded = torch.load(
        path,
        map_location="cpu",
        weights_only=True,
        mmap=True,
    )
    candidates = _tensor_mapping_candidates(loaded)
    if len(candidates) != 1:
        paths = [".".join(candidate_path) or "<root>" for candidate_path, _ in candidates]
        raise ValueError(f"checkpoint has {len(candidates)} tensor mappings: {paths}")
    wrapper_path, state = candidates[0]
    records = []
    total_numel = 0
    total_bytes = 0
    for key in sorted(state):
        tensor = state[key]
        numel = tensor.numel()
        storage_bytes = numel * tensor.element_size()
        total_numel += numel
        total_bytes += storage_bytes
        records.append({
            "key": key,
            "shape": list(tensor.shape),
            "dtype": str(tensor.dtype),
            "numel": numel,
            "storage_bytes": storage_bytes,
        })
    aux_records = [record for record in records if "aux" in str(record["key"]).lower()]
    return {
        "schema_version": "masld-bench-corgi-checkpoint-schema-v1",
        "checkpoint": path.name,
        "checkpoint_size_bytes": expected_size,
        "checkpoint_md5": expected_md5,
        "checkpoint_sha256": expected_sha256,
        "safe_loader": {
            "torch_load_weights_only": True,
            "map_location": "cpu",
            "mmap": True,
            "corgi_imported": False,
            "tensor_values_exported": False,
        },
        "wrapper_path": list(wrapper_path),
        "state_key_count": len(records),
        "total_numel": total_numel,
        "total_storage_bytes": total_bytes,
        "auxiliary_key_records": aux_records,
        "state": records,
    }


def _write_json_exclusive(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o640)
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--expected-size", type=int, required=True)
    parser.add_argument("--expected-md5", required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    payload = inventory_checkpoint(
        args.checkpoint,
        expected_size=args.expected_size,
        expected_md5=args.expected_md5,
        expected_sha256=args.expected_sha256,
    )
    _write_json_exclusive(args.output, payload)


if __name__ == "__main__":
    main()
