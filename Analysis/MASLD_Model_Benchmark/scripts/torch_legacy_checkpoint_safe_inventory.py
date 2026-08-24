#!/usr/bin/env python3
"""Inventory a legacy PyTorch pickle stream without executing pickle opcodes."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import pickletools
import zipfile


class TorchLegacyInventoryError(ValueError):
    """Raised when a legacy checkpoint differs from its admitted structure."""


def _sha256_file(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def inventory(
    path: Path,
    output: Path,
    expected_sha256: str,
    expected_size: int,
) -> dict[str, object]:
    if output.exists() or path.is_symlink() or not path.is_file():
        raise TorchLegacyInventoryError("legacy checkpoint request differs")
    if path.stat().st_size != expected_size or _sha256_file(path) != expected_sha256:
        raise TorchLegacyInventoryError("legacy checkpoint bytes differ")
    if zipfile.is_zipfile(path):
        raise TorchLegacyInventoryError("checkpoint unexpectedly uses ZIP serialization")

    streams: list[dict[str, object]] = []
    globals_seen: set[str] = set()
    with path.open("rb") as handle:
        for stream_index in range(5):
            start = handle.tell()
            opcode_count = 0
            stopped = False
            for opcode, argument, _ in pickletools.genops(handle):
                opcode_count += 1
                if opcode.name in {"GLOBAL", "STACK_GLOBAL"}:
                    globals_seen.add(str(argument))
                if opcode.name == "STOP":
                    stopped = True
            if not stopped or handle.tell() <= start:
                raise TorchLegacyInventoryError("legacy pickle stream framing differs")
            streams.append(
                {
                    "stream_index": stream_index,
                    "start_offset": start,
                    "end_offset": handle.tell(),
                    "opcode_count": opcode_count,
                }
            )
        storage_start = handle.tell()
    if storage_start >= expected_size or len(streams) != 5:
        raise TorchLegacyInventoryError("legacy tensor-storage boundary differs")
    receipt = {
        "schema_version": "masld-bench-safe-torch-legacy-inventory-v1",
        "status": "pass",
        "filename": path.name,
        "sha256": expected_sha256,
        "size_bytes": expected_size,
        "container_type": "PyTorch_legacy_pickle_stream",
        "pickle_streams": streams,
        "pickle_globals": sorted(globals_seen),
        "tensor_storage_start_offset": storage_start,
        "tensor_storage_bytes": expected_size - storage_start,
        "pickle_executed": False,
        "torch_imported": False,
        "checkpoint_deserialized": False,
    }
    output.write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(receipt, sort_keys=True))
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-sha256", required=True)
    parser.add_argument("--expected-size", type=int, required=True)
    arguments = parser.parse_args()
    inventory(
        arguments.input,
        arguments.output,
        arguments.expected_sha256,
        arguments.expected_size,
    )


if __name__ == "__main__":
    main()
