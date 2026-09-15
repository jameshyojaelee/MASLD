#!/usr/bin/env python3
"""Verify and inventory a safetensors file without loading tensor data."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import struct
from typing import Mapping


DTYPE_BYTES = {
    "BOOL": 1,
    "U8": 1,
    "I8": 1,
    "I16": 2,
    "U16": 2,
    "F16": 2,
    "BF16": 2,
    "I32": 4,
    "U32": 4,
    "F32": 4,
    "F8_E4M3": 1,
    "F8_E5M2": 1,
    "F8_E8M0": 1,
    "I64": 8,
    "U64": 8,
    "F64": 8,
    "C64": 8,
    "C128": 16,
}


class SafetensorsInventoryError(ValueError):
    """Raised when a safetensors file does not meet its no-load requirements."""


def _hash(path: Path) -> tuple[str, int]:
    value = sha256()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            value.update(block)
            size += len(block)
    return value.hexdigest(), size


def _elements(shape: list[int]) -> int:
    value = 1
    for dimension in shape:
        if not isinstance(dimension, int) or dimension < 0:
            raise SafetensorsInventoryError("tensor shape differs")
        value *= dimension
    return value


def inventory(
    path: Path,
    output: Path,
    *,
    expected_sha256: str,
    expected_size: int,
    maximum_header_bytes: int = 100_000_000,
    maximum_tensors: int = 1_000_000,
) -> dict[str, object]:
    if output.exists() or path.is_symlink() or not path.is_file():
        raise SafetensorsInventoryError("safetensors request differs")
    digest, size = _hash(path)
    if digest != expected_sha256 or size != expected_size:
        raise SafetensorsInventoryError("safetensors bytes differ")
    with path.open("rb") as handle:
        prefix = handle.read(8)
        if len(prefix) != 8:
            raise SafetensorsInventoryError("safetensors header prefix differs")
        header_length = struct.unpack("<Q", prefix)[0]
        if header_length < 2 or header_length > maximum_header_bytes:
            raise SafetensorsInventoryError("safetensors header length differs")
        header_bytes = handle.read(header_length)
        if len(header_bytes) != header_length:
            raise SafetensorsInventoryError("safetensors header is truncated")
    try:
        header = json.loads(header_bytes.decode("utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SafetensorsInventoryError("safetensors header JSON differs") from exc
    if not isinstance(header, dict):
        raise SafetensorsInventoryError("safetensors header root differs")
    metadata = header.pop("__metadata__", {})
    if not isinstance(metadata, dict) or len(header) > maximum_tensors:
        raise SafetensorsInventoryError("safetensors metadata or tensor count differs")
    data_size = size - 8 - header_length
    rows: list[dict[str, object]] = []
    ranges: list[tuple[int, int, str]] = []
    dtype_counts: dict[str, int] = {}
    parameter_count = 0
    for name, raw in header.items():
        if not isinstance(name, str) or not name or not isinstance(raw, Mapping):
            raise SafetensorsInventoryError("safetensors tensor record differs")
        dtype = raw.get("dtype")
        shape = raw.get("shape")
        offsets = raw.get("data_offsets")
        if dtype not in DTYPE_BYTES or not isinstance(shape, list) or not isinstance(offsets, list) or len(offsets) != 2:
            raise SafetensorsInventoryError(f"safetensors tensor schema differs: {name}")
        start, end = offsets
        if (
            not isinstance(start, int)
            or not isinstance(end, int)
            or start < 0
            or end < start
            or end > data_size
        ):
            raise SafetensorsInventoryError(f"safetensors tensor offsets differ: {name}")
        elements = _elements(shape)
        expected_bytes = elements * DTYPE_BYTES[str(dtype)]
        if end - start != expected_bytes:
            raise SafetensorsInventoryError(f"safetensors tensor byte span differs: {name}")
        ranges.append((start, end, name))
        parameter_count += elements
        dtype_counts[str(dtype)] = dtype_counts.get(str(dtype), 0) + 1
        rows.append(
            {
                "name": name,
                "dtype": dtype,
                "shape": shape,
                "elements": elements,
                "data_offsets": [start, end],
                "size_bytes": expected_bytes,
            }
        )
    ranges.sort()
    for left, right in zip(ranges, ranges[1:]):
        if left[1] > right[0]:
            raise SafetensorsInventoryError(
                f"safetensors tensor ranges overlap: {left[2]} and {right[2]}"
            )
    maximum_end = max((row[1] for row in ranges), default=0)
    if maximum_end != data_size:
        raise SafetensorsInventoryError("safetensors data section has trailing bytes")
    result = {
        "schema_version": "masld-bench-safe-safetensors-inventory-v1",
        "status": "pass",
        "filename": path.name,
        "sha256": digest,
        "size_bytes": size,
        "header_length_bytes": header_length,
        "data_length_bytes": data_size,
        "metadata": metadata,
        "tensor_count": len(rows),
        "parameter_and_buffer_elements": parameter_count,
        "dtype_tensor_counts": dict(sorted(dtype_counts.items())),
        "tensor_data_loaded": False,
        "checkpoint_deserialized": False,
        "tensors": rows,
    }
    output.write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--expected-size", type=int, required=True)
    args = parser.parse_args()
    result = inventory(
        args.input,
        args.output,
        expected_sha256=args.expected_sha256,
        expected_size=args.expected_size,
    )
    print(
        json.dumps(
            {key: value for key, value in result.items() if key != "tensors"},
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
