#!/usr/bin/env python3
"""Inventory a PyTorch ZIP checkpoint without importing pickle or torch."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import stat
from typing import BinaryIO
import zipfile


class TorchCheckpointInventoryError(ValueError):
    """Raised when a PyTorch checkpoint does not meet its no-pickle requirements."""


def _hash_stream(handle: BinaryIO) -> tuple[str, int]:
    value = sha256()
    size = 0
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        value.update(block)
        size += len(block)
    return value.hexdigest(), size


def _hash_file(path: Path) -> tuple[str, int]:
    with path.open("rb") as handle:
        return _hash_stream(handle)


def _safe_name(value: str) -> str:
    if not value or value.startswith("/") or "\\" in value or "\x00" in value:
        raise TorchCheckpointInventoryError("checkpoint member has an unsafe name")
    parts = PurePosixPath(value).parts
    if not parts or any(part in {"", ".", ".."} for part in parts):
        raise TorchCheckpointInventoryError("checkpoint member path traverses")
    return PurePosixPath(*parts).as_posix()


def inventory(
    path: Path,
    output: Path,
    *,
    expected_sha256: str,
    expected_size: int,
    maximum_members: int = 1_000_000,
    maximum_uncompressed_bytes: int = 34_359_738_368,
) -> dict[str, object]:
    if output.exists() or path.is_symlink() or not path.is_file():
        raise TorchCheckpointInventoryError("checkpoint request differs")
    digest, size = _hash_file(path)
    if digest != expected_sha256 or size != expected_size:
        raise TorchCheckpointInventoryError("checkpoint bytes differ")
    if not zipfile.is_zipfile(path):
        raise TorchCheckpointInventoryError(
            "checkpoint is not a ZIP-format torch container; raw pickle inspection is blocked"
        )
    rows: list[dict[str, object]] = []
    seen: set[str] = set()
    uncompressed_total = 0
    with zipfile.ZipFile(path, mode="r") as archive:
        infos = archive.infolist()
        if not infos or len(infos) > maximum_members:
            raise TorchCheckpointInventoryError("checkpoint member count differs")
        for info in infos:
            name = _safe_name(info.filename)
            if name in seen:
                raise TorchCheckpointInventoryError("checkpoint has duplicate members")
            seen.add(name)
            mode = (info.external_attr >> 16) & 0xFFFF
            if stat.S_ISLNK(mode):
                raise TorchCheckpointInventoryError("checkpoint has a symbolic link member")
            is_directory = info.is_dir()
            if is_directory and info.file_size:
                raise TorchCheckpointInventoryError("checkpoint directory has data")
            record: dict[str, object] = {
                "path": name,
                "type": "directory" if is_directory else "regular_file",
                "size_bytes": int(info.file_size),
                "compressed_size_bytes": int(info.compress_size),
                "crc32": f"{info.CRC:08x}",
                "compression_type": int(info.compress_type),
            }
            if not is_directory:
                uncompressed_total += int(info.file_size)
                if uncompressed_total > maximum_uncompressed_bytes:
                    raise TorchCheckpointInventoryError(
                        "checkpoint uncompressed-byte ceiling exceeded"
                    )
                with archive.open(info, mode="r") as handle:
                    member_sha256, observed_size = _hash_stream(handle)
                if observed_size != info.file_size:
                    raise TorchCheckpointInventoryError("checkpoint member size differs")
                record["sha256"] = member_sha256
            rows.append(record)
    pickle_members = [
        row["path"]
        for row in rows
        if row["type"] == "regular_file" and str(row["path"]).endswith(".pkl")
    ]
    if len(pickle_members) != 1 or not str(pickle_members[0]).endswith("data.pkl"):
        raise TorchCheckpointInventoryError("torch pickle metadata member differs")
    storage_members = [
        row["path"]
        for row in rows
        if row["type"] == "regular_file" and "/data/" in str(row["path"])
    ]
    if not storage_members:
        raise TorchCheckpointInventoryError("torch tensor storage members are absent")
    result = {
        "schema_version": "masld-bench-safe-torch-zip-inventory-v1",
        "status": "pass",
        "filename": path.name,
        "sha256": digest,
        "size_bytes": size,
        "container_type": "PyTorch_ZIP_with_pickle_metadata",
        "member_count": len(rows),
        "regular_file_count": sum(row["type"] == "regular_file" for row in rows),
        "uncompressed_regular_bytes": uncompressed_total,
        "pickle_metadata_members": pickle_members,
        "tensor_storage_member_count": len(storage_members),
        "archive_extracted": False,
        "pickle_imported": False,
        "torch_imported": False,
        "checkpoint_deserialized": False,
        "members": rows,
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
    parser.add_argument("--maximum-uncompressed-bytes", type=int, default=34_359_738_368)
    args = parser.parse_args()
    result = inventory(
        args.input,
        args.output,
        expected_sha256=args.expected_sha256,
        expected_size=args.expected_size,
        maximum_uncompressed_bytes=args.maximum_uncompressed_bytes,
    )
    print(
        json.dumps(
            {key: value for key, value in result.items() if key != "members"},
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
