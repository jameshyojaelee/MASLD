#!/usr/bin/env python3
"""Inventory a remote ZIP central directory using bounded HTTP range requests."""

from __future__ import annotations

import argparse
import json
from pathlib import Path, PurePosixPath
import re
import stat
import struct
from typing import Protocol
import urllib.error
import urllib.request


EOCD_SIGNATURE = b"PK\x05\x06"
ZIP64_LOCATOR_SIGNATURE = b"PK\x06\x07"
ZIP64_EOCD_SIGNATURE = b"PK\x06\x06"
CENTRAL_SIGNATURE = b"PK\x01\x02"
MAX_EOCD_SEARCH = 22 + 65_535 + 20 + 64
MAX_CENTRAL_DIRECTORY_BYTES = 256 * 1024 * 1024
MAX_ENTRIES = 1_000_000


class RemoteZipInventoryError(ValueError):
    """Raised when the remote object or ZIP structure differs."""


class RangeReader(Protocol):
    size: int

    def read(self, start: int, length: int) -> bytes:
        """Return exactly ``length`` bytes beginning at ``start``."""


class HTTPRangeReader:
    """Strict immutable-object reader backed by HTTP byte ranges."""

    def __init__(self, url: str, expected_size: int) -> None:
        if not url.startswith("https://") or expected_size <= 0:
            raise RemoteZipInventoryError("remote ZIP request differs")
        self.url = url
        self.size = expected_size
        self.request_count = 0
        self.bytes_downloaded = 0

    def read(self, start: int, length: int) -> bytes:
        if start < 0 or length <= 0 or start + length > self.size:
            raise RemoteZipInventoryError("range request is outside the object")
        end = start + length - 1
        request = urllib.request.Request(
            self.url,
            headers={
                "Accept-Encoding": "identity",
                "Range": f"bytes={start}-{end}",
                "User-Agent": "masld-bench-scooby-inventory/1",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=120) as response:
                status = int(getattr(response, "status", 0))
                content_range = response.headers.get("Content-Range", "")
                content_encoding = response.headers.get("Content-Encoding", "identity")
                payload = response.read(length + 1)
        except (urllib.error.URLError, TimeoutError, OSError) as exc:
            raise RemoteZipInventoryError("remote range request failed") from exc
        expected_content_range = f"bytes {start}-{end}/{self.size}"
        if (
            status != 206
            or content_range != expected_content_range
            or content_encoding not in {"", "identity", None}
            or len(payload) != length
        ):
            raise RemoteZipInventoryError(
                "remote server did not honor the exact identity range"
            )
        self.request_count += 1
        self.bytes_downloaded += len(payload)
        return payload


class BytesRangeReader:
    """In-memory strict range reader used by unit fixtures."""

    def __init__(self, payload: bytes) -> None:
        self.payload = payload
        self.size = len(payload)
        self.request_count = 0
        self.bytes_downloaded = 0

    def read(self, start: int, length: int) -> bytes:
        if start < 0 or length <= 0 or start + length > self.size:
            raise RemoteZipInventoryError("range request is outside the object")
        self.request_count += 1
        self.bytes_downloaded += length
        return self.payload[start : start + length]


def _find_eocd(reader: RangeReader) -> tuple[int, bytes]:
    tail_length = min(reader.size, MAX_EOCD_SEARCH)
    tail_start = reader.size - tail_length
    tail = reader.read(tail_start, tail_length)
    search_end = len(tail)
    while True:
        index = tail.rfind(EOCD_SIGNATURE, 0, search_end)
        if index < 0:
            raise RemoteZipInventoryError("ZIP end-of-central-directory is missing")
        if index + 22 <= len(tail):
            comment_length = struct.unpack_from("<H", tail, index + 20)[0]
            absolute = tail_start + index
            if absolute + 22 + comment_length == reader.size:
                return absolute, tail[index : index + 22 + comment_length]
        search_end = index


def _central_directory_contract(reader: RangeReader) -> dict[str, int | bool]:
    eocd_offset, eocd = _find_eocd(reader)
    (
        _signature,
        disk_number,
        central_disk,
        entries_on_disk,
        entries_total,
        central_size,
        central_offset,
        _comment_length,
    ) = struct.unpack_from("<4s4H2IH", eocd, 0)
    if disk_number != 0 or central_disk != 0 or entries_on_disk != entries_total:
        raise RemoteZipInventoryError("multi-disk ZIP archives are unsupported")
    zip64 = (
        entries_total == 0xFFFF
        or central_size == 0xFFFFFFFF
        or central_offset == 0xFFFFFFFF
    )
    if zip64:
        if eocd_offset < 20:
            raise RemoteZipInventoryError("ZIP64 locator is missing")
        locator = reader.read(eocd_offset - 20, 20)
        signature, zip64_disk, zip64_offset, total_disks = struct.unpack(
            "<4sIQI", locator
        )
        if (
            signature != ZIP64_LOCATOR_SIGNATURE
            or zip64_disk != 0
            or total_disks != 1
            or zip64_offset + 56 > reader.size
        ):
            raise RemoteZipInventoryError("ZIP64 locator differs")
        fixed = reader.read(zip64_offset, 56)
        (
            signature,
            record_size,
            _made_by,
            _needed,
            disk_number,
            central_disk,
            entries_on_disk,
            entries_total,
            central_size,
            central_offset,
        ) = struct.unpack("<4sQ2H2I4Q", fixed)
        if (
            signature != ZIP64_EOCD_SIGNATURE
            or record_size < 44
            or disk_number != 0
            or central_disk != 0
            or entries_on_disk != entries_total
        ):
            raise RemoteZipInventoryError("ZIP64 central-directory record differs")
    if (
        entries_total > MAX_ENTRIES
        or central_size > MAX_CENTRAL_DIRECTORY_BYTES
        or central_offset + central_size > reader.size
    ):
        raise RemoteZipInventoryError("central directory exceeds the bounded contract")
    return {
        "zip64": zip64,
        "entries_total": int(entries_total),
        "central_size": int(central_size),
        "central_offset": int(central_offset),
        "eocd_offset": eocd_offset,
    }


def _zip64_values(
    extra: bytes,
    compressed_size: int,
    uncompressed_size: int,
    local_offset: int,
    disk_start: int,
) -> tuple[int, int, int, int]:
    fields: dict[int, bytes] = {}
    position = 0
    while position < len(extra):
        if position + 4 > len(extra):
            raise RemoteZipInventoryError("central extra field is truncated")
        field_id, field_length = struct.unpack_from("<HH", extra, position)
        position += 4
        if position + field_length > len(extra) or field_id in fields:
            raise RemoteZipInventoryError("central extra field differs")
        fields[field_id] = extra[position : position + field_length]
        position += field_length
    needs_zip64 = (
        uncompressed_size == 0xFFFFFFFF,
        compressed_size == 0xFFFFFFFF,
        local_offset == 0xFFFFFFFF,
        disk_start == 0xFFFF,
    )
    if not any(needs_zip64):
        return compressed_size, uncompressed_size, local_offset, disk_start
    payload = fields.get(0x0001)
    if payload is None:
        raise RemoteZipInventoryError("ZIP64 member extra field is missing")
    cursor = 0

    def take(width: int) -> int:
        nonlocal cursor
        if cursor + width > len(payload):
            raise RemoteZipInventoryError("ZIP64 member extra field is truncated")
        fmt = "<Q" if width == 8 else "<I"
        value = int(struct.unpack_from(fmt, payload, cursor)[0])
        cursor += width
        return value

    if needs_zip64[0]:
        uncompressed_size = take(8)
    if needs_zip64[1]:
        compressed_size = take(8)
    if needs_zip64[2]:
        local_offset = take(8)
    if needs_zip64[3]:
        disk_start = take(4)
    return compressed_size, uncompressed_size, local_offset, disk_start


def _safe_member_name(raw_name: bytes, flags: int) -> str:
    encoding = "utf-8" if flags & 0x0800 else "cp437"
    try:
        name = raw_name.decode(encoding)
    except UnicodeDecodeError as exc:
        raise RemoteZipInventoryError("member filename encoding differs") from exc
    path = PurePosixPath(name)
    if (
        not name
        or "\x00" in name
        or "\\" in name
        or path.is_absolute()
        or ".." in path.parts
    ):
        raise RemoteZipInventoryError("unsafe ZIP member path")
    return name


def _context_candidate(name: str) -> bool:
    normalized = name.lower()
    tokens = (
        "embedding",
        "latent",
        "scpoli",
        "scglue",
        "multivi",
        "multi_vi",
        "poisson",
        "query_model",
        "reference_model",
    )
    suffixes = (".h5ad", ".h5mu", ".h5", ".npz", ".parquet", ".pkl", ".pt", ".ckpt")
    return any(token in normalized for token in tokens) or normalized.endswith(suffixes)


def inventory_reader(reader: RangeReader) -> dict[str, object]:
    contract = _central_directory_contract(reader)
    central = reader.read(
        int(contract["central_offset"]), int(contract["central_size"])
    )
    entries = []
    position = 0
    for index in range(int(contract["entries_total"])):
        if position + 46 > len(central):
            raise RemoteZipInventoryError("central directory is truncated")
        fields = struct.unpack_from("<4s6H3I5H2I", central, position)
        (
            signature,
            version_made,
            _version_needed,
            flags,
            compression,
            _mtime,
            _mdate,
            crc32,
            compressed_size,
            uncompressed_size,
            filename_length,
            extra_length,
            comment_length,
            disk_start,
            _internal_attributes,
            external_attributes,
            local_offset,
        ) = fields
        if signature != CENTRAL_SIGNATURE or flags & 0x0001:
            raise RemoteZipInventoryError("central member signature or encryption differs")
        payload_start = position + 46
        payload_end = payload_start + filename_length + extra_length + comment_length
        if payload_end > len(central):
            raise RemoteZipInventoryError("central member payload is truncated")
        raw_name = central[payload_start : payload_start + filename_length]
        extra = central[
            payload_start + filename_length : payload_start + filename_length + extra_length
        ]
        name = _safe_member_name(raw_name, flags)
        compressed_size, uncompressed_size, local_offset, disk_start = _zip64_values(
            extra, compressed_size, uncompressed_size, local_offset, disk_start
        )
        if disk_start != 0 or local_offset >= reader.size:
            raise RemoteZipInventoryError("central member offset differs")
        unix_mode = (external_attributes >> 16) & 0xFFFF
        if (version_made >> 8) == 3 and stat.S_IFMT(unix_mode) == stat.S_IFLNK:
            raise RemoteZipInventoryError("symlink member is not admitted")
        entries.append(
            {
                "index": index,
                "path": name,
                "compressed_size": compressed_size,
                "uncompressed_size": uncompressed_size,
                "compression": compression,
                "crc32": f"{crc32:08x}",
                "local_header_offset": local_offset,
                "directory": name.endswith("/"),
                "context_candidate": _context_candidate(name),
            }
        )
        position = payload_end
    if position != len(central):
        raise RemoteZipInventoryError("central directory contains trailing bytes")
    candidates = [entry for entry in entries if entry["context_candidate"]]
    return {
        "schema_version": "masld-bench-remote-zip-inventory-v1",
        "status": "pass",
        "archive_size_bytes": reader.size,
        "zip64": contract["zip64"],
        "central_directory_offset": contract["central_offset"],
        "central_directory_size_bytes": contract["central_size"],
        "entry_count": len(entries),
        "members": entries,
        "context_candidate_count": len(candidates),
        "context_candidates": candidates,
        "archive_payload_downloaded": False,
        "archive_extracted": False,
    }


def inventory_remote(
    url: str, expected_size: int, expected_md5: str, output: Path
) -> dict[str, object]:
    if output.exists() or not re.fullmatch(r"[0-9a-f]{32}", expected_md5):
        raise RemoteZipInventoryError("output or declared MD5 contract differs")
    reader = HTTPRangeReader(url, expected_size)
    receipt = inventory_reader(reader)
    receipt.update(
        {
            "remote_url": url,
            "declared_archive_md5": expected_md5,
            "declared_archive_md5_verified_from_full_payload": False,
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
    parser.add_argument("--expected-md5", required=True)
    parser.add_argument("--output", required=True, type=Path)
    arguments = parser.parse_args()
    inventory_remote(
        arguments.url,
        arguments.expected_size,
        arguments.expected_md5,
        arguments.output,
    )


if __name__ == "__main__":
    main()
