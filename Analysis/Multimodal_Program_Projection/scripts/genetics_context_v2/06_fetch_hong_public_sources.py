#!/usr/bin/env python3
"""Fetch the bounded public Hong et al. context-eQTL source contract.

The Zenodo deposit is a 6.6 GB ZIP with HTTP byte-range support.  GEN-03 only
needs the 23 source interaction tables, the authors' significant-call table,
and the reported quartet table.  This downloader parses a frozen central
directory, retrieves those entries by byte range, and verifies each ZIP CRC32.
It never reads repository-local sc-eQTL outcome products.
"""

from __future__ import annotations

import argparse
import binascii
import os
import struct
import tempfile
import time
import urllib.error
import urllib.request
import zlib
from dataclasses import dataclass
from pathlib import Path

from genetics_common import (
    ContractError,
    DEFAULT_CANDIDATE_ROOT,
    PROJECT_ROOT,
    assert_candidate_root,
    assert_owned_output,
    atomic_write_tsv,
    sha256_file,
)


ZENODO_RECORD_ID = "14586466"
ZENODO_ARCHIVE_NAME = "Zenodo_250212.zip"
ZENODO_ARCHIVE_BYTES = 6_647_735_352
ZENODO_ARCHIVE_MD5 = "1fd1a014829a8d6a0f2e97fcb07e69fb"
ZENODO_URL = (
    "https://zenodo.org/api/records/14586466/files/Zenodo_250212.zip/content"
)
CENTRAL_DIRECTORY_REL = Path("work/public_sources/zenodo_zip_central_directory.bin")
OUTPUT_REL = Path("work/public_sources/zenodo_context")

EXPECTED_CONTEXTS = {
    "hepatocyte": {"disease_group3", "Hep-M3", "Hep-M4", "Hep-M6", "Hep-M8", "Hep-M12", "Hep-M14"},
    "cholangiocyte": {"disease_group3", "Chol-M1", "Chol-M2", "Chol-M3", "Chol-M4", "Chol-M6", "Chol-M7"},
    "stellate_cell": {"disease_group3", "HSC-M1", "HSC-M2", "HSC-M3"},
    "endothelial_cell": {"disease_group3", "endo-M1", "endo-M4", "endo-M7", "endo-M9"},
}

ANNOTATION_ENTRY = "Zenodo/sceQTL_output/significant_sc_eQTLs_with_annotations.txt.gz"
QUARTET_ENTRY = "Zenodo/quartets/significant_quartets_n601.txt.gz"
CONTEXT_PREFIX = "Zenodo/sceQTL_output/ieQTLs/interaction."


@dataclass(frozen=True)
class ZipEntry:
    name: str
    uncompressed_bytes: int
    compressed_bytes: int
    compression: int
    local_header_offset: int
    crc32: int
    flags: int


def parse_zip64_extra(
    extra: bytes,
    uncompressed_bytes: int,
    compressed_bytes: int,
    local_offset: int,
    disk_start: int,
) -> tuple[int, int, int, int]:
    position = 0
    while position + 4 <= len(extra):
        tag, size = struct.unpack_from("<HH", extra, position)
        payload = extra[position + 4 : position + 4 + size]
        position += 4 + size
        if tag != 0x0001:
            continue
        cursor = 0
        if uncompressed_bytes == 0xFFFFFFFF:
            uncompressed_bytes = struct.unpack_from("<Q", payload, cursor)[0]
            cursor += 8
        if compressed_bytes == 0xFFFFFFFF:
            compressed_bytes = struct.unpack_from("<Q", payload, cursor)[0]
            cursor += 8
        if local_offset == 0xFFFFFFFF:
            local_offset = struct.unpack_from("<Q", payload, cursor)[0]
            cursor += 8
        if disk_start == 0xFFFF:
            disk_start = struct.unpack_from("<L", payload, cursor)[0]
        break
    return uncompressed_bytes, compressed_bytes, local_offset, disk_start


def parse_central_directory(path: Path) -> dict[str, ZipEntry]:
    payload = path.read_bytes()
    position = 0
    entries: dict[str, ZipEntry] = {}
    while position < len(payload):
        if payload[position : position + 4] != b"PK\x01\x02":
            raise ContractError(f"invalid ZIP central-directory signature at byte {position}")
        fields = struct.unpack_from("<4s6H3L5H2L", payload, position)
        (
            _,
            _version_made,
            _version_needed,
            flags,
            compression,
            _mod_time,
            _mod_date,
            crc32,
            compressed_bytes,
            uncompressed_bytes,
            name_length,
            extra_length,
            comment_length,
            disk_start,
            _internal_attributes,
            _external_attributes,
            local_offset,
        ) = fields
        name_start = position + 46
        encoding = "utf-8" if flags & 0x800 else "cp437"
        name = payload[name_start : name_start + name_length].decode(encoding)
        extra_start = name_start + name_length
        extra = payload[extra_start : extra_start + extra_length]
        uncompressed_bytes, compressed_bytes, local_offset, disk_start = parse_zip64_extra(
            extra,
            uncompressed_bytes,
            compressed_bytes,
            local_offset,
            disk_start,
        )
        if disk_start != 0:
            raise ContractError(f"multi-disk ZIP entry is unsupported: {name}")
        if name in entries:
            raise ContractError(f"duplicate ZIP entry: {name}")
        entries[name] = ZipEntry(
            name=name,
            uncompressed_bytes=uncompressed_bytes,
            compressed_bytes=compressed_bytes,
            compression=compression,
            local_header_offset=local_offset,
            crc32=crc32,
            flags=flags,
        )
        position = extra_start + extra_length + comment_length
    if position != len(payload):
        raise ContractError("central-directory parser did not consume the complete payload")
    return entries


def context_from_entry(name: str) -> tuple[str, str] | None:
    if not name.startswith(CONTEXT_PREFIX) or not name.endswith(".txt.gz"):
        return None
    stem = name[len(CONTEXT_PREFIX) : -len(".txt.gz")]
    for cell_type in sorted(EXPECTED_CONTEXTS, key=len, reverse=True):
        prefix = f"{cell_type}_"
        if stem.startswith(prefix):
            return cell_type, stem[len(prefix) :]
    raise ContractError(f"unrecognized Hong context filename: {name}")


def select_entries(entries: dict[str, ZipEntry]) -> list[tuple[ZipEntry, str, str, str]]:
    selected: list[tuple[ZipEntry, str, str, str]] = []
    observed: dict[str, set[str]] = {cell_type: set() for cell_type in EXPECTED_CONTEXTS}
    for name, entry in entries.items():
        parsed = context_from_entry(name)
        if parsed is None:
            continue
        cell_type, context = parsed
        observed[cell_type].add(context)
        selected.append((entry, "interaction_summary", cell_type, context))
    if observed != EXPECTED_CONTEXTS:
        raise ContractError(
            "Hong context inventory drift: "
            + "; ".join(
                f"{cell_type} expected={sorted(EXPECTED_CONTEXTS[cell_type])} "
                f"observed={sorted(observed[cell_type])}"
                for cell_type in EXPECTED_CONTEXTS
                if observed[cell_type] != EXPECTED_CONTEXTS[cell_type]
            )
        )
    for name, role in ((ANNOTATION_ENTRY, "significant_call_annotation"), (QUARTET_ENTRY, "quartet_table")):
        if name not in entries:
            raise ContractError(f"required Hong source entry is absent: {name}")
        selected.append((entries[name], role, "", ""))
    return sorted(selected, key=lambda item: item[0].name)


def request_range(start: int, end: int, attempts: int = 6) -> bytes:
    expected = end - start + 1
    last_error: Exception | None = None
    for attempt in range(attempts):
        request = urllib.request.Request(
            ZENODO_URL,
            headers={
                "Range": f"bytes={start}-{end}",
                "User-Agent": "MASLD-Hong-source-audit/1.0",
            },
        )
        try:
            with urllib.request.urlopen(request, timeout=180) as response:
                if response.status != 206:
                    raise ContractError(f"range request returned HTTP {response.status}")
                content_range = response.headers.get("Content-Range", "")
                expected_range = f"bytes {start}-{end}/{ZENODO_ARCHIVE_BYTES}"
                if content_range != expected_range:
                    raise ContractError(
                        f"unexpected Content-Range {content_range!r}; expected {expected_range!r}"
                    )
                payload = response.read()
            if len(payload) != expected:
                raise ContractError(
                    f"short range response for {start}-{end}: {len(payload)} != {expected}"
                )
            return payload
        except (OSError, urllib.error.URLError, ContractError) as error:
            last_error = error
            if attempt + 1 == attempts:
                break
            time.sleep(min(2**attempt, 30))
    raise ContractError(f"failed byte-range request {start}-{end}: {last_error}")


def extract_entry(entry: ZipEntry, output: Path) -> None:
    header = request_range(entry.local_header_offset, entry.local_header_offset + 65_535)
    if header[:4] != b"PK\x03\x04":
        raise ContractError(f"invalid local header for {entry.name}")
    local_fields = struct.unpack_from("<4s5H3L2H", header, 0)
    name_length, extra_length = local_fields[-2:]
    name_encoding = "utf-8" if entry.flags & 0x800 else "cp437"
    local_name = header[30 : 30 + name_length].decode(name_encoding)
    if local_name != entry.name:
        raise ContractError(f"central/local name mismatch: {entry.name} != {local_name}")
    data_start = entry.local_header_offset + 30 + name_length + extra_length
    compressed = request_range(data_start, data_start + entry.compressed_bytes - 1)
    if entry.compression == 8:
        payload = zlib.decompress(compressed, -15)
    elif entry.compression == 0:
        payload = compressed
    else:
        raise ContractError(f"unsupported ZIP compression method {entry.compression}: {entry.name}")
    if len(payload) != entry.uncompressed_bytes:
        raise ContractError(
            f"uncompressed-size mismatch for {entry.name}: "
            f"{len(payload)} != {entry.uncompressed_bytes}"
        )
    observed_crc = binascii.crc32(payload) & 0xFFFFFFFF
    if observed_crc != entry.crc32:
        raise ContractError(
            f"CRC mismatch for {entry.name}: {observed_crc:08x} != {entry.crc32:08x}"
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        mode="wb",
        dir=output.parent,
        prefix=f".{output.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
        handle.write(payload)
    os.replace(temporary, output)


def file_crc32(path: Path, chunk_bytes: int = 8 * 1024 * 1024) -> int:
    checksum = 0
    with path.open("rb") as handle:
        while True:
            block = handle.read(chunk_bytes)
            if not block:
                break
            checksum = binascii.crc32(block, checksum)
    return checksum & 0xFFFFFFFF


def output_is_valid(path: Path, entry: ZipEntry) -> bool:
    return (
        path.is_file()
        and path.stat().st_size == entry.uncompressed_bytes
        and file_crc32(path) == entry.crc32
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--candidate-root", type=Path, default=DEFAULT_CANDIDATE_ROOT)
    parser.add_argument(
        "--inventory-only",
        action="store_true",
        help="Validate the frozen central directory and emit the bounded inventory without downloading entries.",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    project_root = args.project_root.resolve()
    if project_root != PROJECT_ROOT.resolve():
        raise ContractError(f"project root must be exactly {PROJECT_ROOT.resolve()}")
    candidate = assert_candidate_root(project_root, args.candidate_root)
    central_path = candidate / CENTRAL_DIRECTORY_REL
    output_dir = candidate / OUTPUT_REL
    assert_owned_output(candidate, central_path)
    assert_owned_output(candidate, output_dir)
    if not central_path.is_file():
        raise ContractError(f"frozen Zenodo central directory is missing: {central_path}")

    entries = parse_central_directory(central_path)
    selected = select_entries(entries)
    rows = []
    for entry, role, cell_type, context in selected:
        output = output_dir / Path(entry.name).name
        if not args.inventory_only and not output_is_valid(output, entry):
            extract_entry(entry, output)
        available = output_is_valid(output, entry)
        rows.append(
            {
                "archive_entry": entry.name,
                "role": role,
                "cell_type": cell_type,
                "context": context,
                "archive_uncompressed_bytes": entry.uncompressed_bytes,
                "archive_compressed_bytes": entry.compressed_bytes,
                "archive_crc32": f"{entry.crc32:08x}",
                "local_header_offset": entry.local_header_offset,
                "local_path": str(output.relative_to(candidate)) if available else "",
                "local_sha256": sha256_file(output) if available else "",
                "retrieval_status": "crc_verified" if available else "not_retrieved",
                "source_url": ZENODO_URL,
                "zenodo_record_id": ZENODO_RECORD_ID,
                "archive_name": ZENODO_ARCHIVE_NAME,
                "archive_bytes": ZENODO_ARCHIVE_BYTES,
                "published_archive_md5": ZENODO_ARCHIVE_MD5,
            }
        )
    inventory = candidate / "hong_public_archive_inventory.tsv"
    assert_owned_output(candidate, inventory)
    atomic_write_tsv(inventory, rows, list(rows[0]))
    complete = all(row["retrieval_status"] == "crc_verified" for row in rows)
    print(
        f"{'PASS' if complete else 'INVENTORY_ONLY'}: Hong public source inventory "
        f"contains {len(rows)} bounded entries; retrieved={sum(row['retrieval_status'] == 'crc_verified' for row in rows)}"
    )


if __name__ == "__main__":
    main()
