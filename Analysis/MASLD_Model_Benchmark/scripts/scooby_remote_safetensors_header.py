#!/usr/bin/env python3
"""Validate a remote safetensors header without downloading tensor payloads."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import struct
from typing import Any, Protocol
import urllib.error
import urllib.request


MAX_HEADER_BYTES = 16 * 1024 * 1024
DTYPE_BYTES = {
    "BOOL": 1,
    "U8": 1,
    "I8": 1,
    "F8_E4M3": 1,
    "F8_E5M2": 1,
    "I16": 2,
    "U16": 2,
    "F16": 2,
    "BF16": 2,
    "I32": 4,
    "U32": 4,
    "F32": 4,
    "F64": 8,
    "I64": 8,
    "U64": 8,
}


class SafetensorsHeaderError(ValueError):
    """Raised when a remote safetensors header differs."""


class RangeReader(Protocol):
    size: int
    request_count: int
    bytes_downloaded: int

    def read(self, start: int, length: int) -> bytes:
        """Return an exact byte range."""


class HTTPRangeReader:
    def __init__(self, url: str, expected_size: int) -> None:
        if not url.startswith("https://") or expected_size <= 8:
            raise SafetensorsHeaderError("remote safetensors request differs")
        self.url = url
        self.size = expected_size
        self.request_count = 0
        self.bytes_downloaded = 0

    def read(self, start: int, length: int) -> bytes:
        if start < 0 or length <= 0 or start + length > self.size:
            raise SafetensorsHeaderError("range request is outside the object")
        end = start + length - 1
        request = urllib.request.Request(
            self.url,
            headers={
                "Accept-Encoding": "identity",
                "Range": f"bytes={start}-{end}",
                "User-Agent": "masld-bench-scooby-header/1",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                status = int(getattr(response, "status", 0))
                content_range = response.headers.get("Content-Range", "")
                content_encoding = response.headers.get("Content-Encoding", "identity")
                payload = response.read(length + 1)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            raise SafetensorsHeaderError("remote range request failed") from exc
        if (
            status != 206
            or content_range != f"bytes {start}-{end}/{self.size}"
            or content_encoding not in {"", "identity", None}
            or len(payload) != length
        ):
            raise SafetensorsHeaderError("remote server did not honor the exact range")
        self.request_count += 1
        self.bytes_downloaded += len(payload)
        return payload


class BytesRangeReader:
    def __init__(self, payload: bytes) -> None:
        self.payload = payload
        self.size = len(payload)
        self.request_count = 0
        self.bytes_downloaded = 0

    def read(self, start: int, length: int) -> bytes:
        if start < 0 or length <= 0 or start + length > self.size:
            raise SafetensorsHeaderError("range request is outside the object")
        self.request_count += 1
        self.bytes_downloaded += length
        return self.payload[start : start + length]


def _unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in pairs:
        if key in result:
            raise SafetensorsHeaderError("safetensors header has duplicate keys")
        result[key] = value
    return result


def _tensor_nbytes(dtype: str, shape: list[int]) -> int:
    if dtype not in DTYPE_BYTES or not isinstance(shape, list):
        raise SafetensorsHeaderError("tensor dtype or shape differs")
    elements = 1
    for dimension in shape:
        if not isinstance(dimension, int) or dimension < 0:
            raise SafetensorsHeaderError("tensor shape differs")
        elements *= dimension
        if elements > 10**12:
            raise SafetensorsHeaderError("tensor shape exceeds the bounded contract")
    return elements * DTYPE_BYTES[dtype]


def inspect_reader(reader: RangeReader, expected_parameters: int) -> dict[str, object]:
    if expected_parameters <= 0:
        raise SafetensorsHeaderError("expected parameter count differs")
    prefix = reader.read(0, 8)
    header_length = int(struct.unpack("<Q", prefix)[0])
    if header_length <= 2 or header_length > MAX_HEADER_BYTES or 8 + header_length >= reader.size:
        raise SafetensorsHeaderError("safetensors header length differs")
    raw_header = reader.read(8, header_length)
    try:
        header = json.loads(
            raw_header.decode("utf-8"), object_pairs_hook=_unique_object
        )
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise SafetensorsHeaderError("safetensors header JSON differs") from exc
    if not isinstance(header, dict):
        raise SafetensorsHeaderError("safetensors header is not an object")
    metadata = header.pop("__metadata__", {})
    if not isinstance(metadata, dict):
        raise SafetensorsHeaderError("safetensors metadata differs")
    tensors = []
    total_parameters = 0
    previous_end = 0
    data_bytes = reader.size - 8 - header_length
    for name, spec in sorted(
        header.items(), key=lambda item: item[1].get("data_offsets", [-1, -1])[0]
        if isinstance(item[1], dict)
        else -1
    ):
        if not name or not isinstance(spec, dict) or set(spec) != {"dtype", "shape", "data_offsets"}:
            raise SafetensorsHeaderError("tensor schema differs")
        dtype = spec["dtype"]
        shape = spec["shape"]
        offsets = spec["data_offsets"]
        if (
            not isinstance(dtype, str)
            or not isinstance(offsets, list)
            or len(offsets) != 2
            or not all(isinstance(value, int) for value in offsets)
        ):
            raise SafetensorsHeaderError("tensor offsets differ")
        start, end = offsets
        expected_bytes = _tensor_nbytes(dtype, shape)
        if start != previous_end or end < start or end - start != expected_bytes or end > data_bytes:
            raise SafetensorsHeaderError("tensor byte spans are not contiguous and exact")
        elements = expected_bytes // DTYPE_BYTES[dtype]
        total_parameters += elements
        tensors.append(
            {
                "name": name,
                "dtype": dtype,
                "shape": shape,
                "data_offsets": offsets,
                "parameter_count": elements,
            }
        )
        previous_end = end
    if previous_end != data_bytes or total_parameters != expected_parameters:
        raise SafetensorsHeaderError("payload span or parameter count differs")
    return {
        "schema_version": "masld-bench-remote-safetensors-header-v1",
        "status": "pass",
        "file_size_bytes": reader.size,
        "header_size_bytes": header_length,
        "tensor_payload_size_bytes": data_bytes,
        "tensor_count": len(tensors),
        "parameter_count": total_parameters,
        "dtypes": sorted({tensor["dtype"] for tensor in tensors}),
        "metadata": metadata,
        "tensors": tensors,
        "checkpoint_payload_downloaded": False,
        "checkpoint_deserialized": False,
        "model_forward_executed": False,
    }


def inspect_remote(
    url: str,
    expected_size: int,
    expected_parameters: int,
    reported_parameters: int,
    expected_sha256: str,
    output: Path,
) -> dict[str, object]:
    if (
        output.exists()
        or len(expected_sha256) != 64
        or reported_parameters <= 0
    ):
        raise SafetensorsHeaderError("output or declared SHA-256 differs")
    reader = HTTPRangeReader(url, expected_size)
    receipt = inspect_reader(reader, expected_parameters)
    receipt.update(
        {
            "remote_url": url,
            "declared_file_sha256": expected_sha256,
            "declared_file_sha256_verified_from_full_payload": False,
            "registry_reported_parameter_count": reported_parameters,
            "registry_to_header_element_difference": (
                receipt["parameter_count"] - reported_parameters
            ),
            "http_range_request_count": reader.request_count,
            "http_bytes_downloaded": reader.bytes_downloaded,
        }
    )
    output.parent.mkdir(parents=True, exist_ok=True, mode=0o750)
    output.write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--url", required=True)
    parser.add_argument("--expected-size", required=True, type=int)
    parser.add_argument("--expected-parameters", required=True, type=int)
    parser.add_argument("--reported-parameters", required=True, type=int)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    inspect_remote(
        arguments.url,
        arguments.expected_size,
        arguments.expected_parameters,
        arguments.reported_parameters,
        arguments.expected_sha256,
        arguments.output,
    )


if __name__ == "__main__":
    main()
