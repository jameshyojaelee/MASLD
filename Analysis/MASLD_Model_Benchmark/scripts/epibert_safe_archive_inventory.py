#!/usr/bin/env python3
"""Hash and safely inventory a tar archive without extracting any member."""

from __future__ import annotations

import argparse
from hashlib import md5, sha256
import json
from pathlib import Path, PurePosixPath
import tarfile
from typing import BinaryIO


class ArchiveInventoryError(ValueError):
    """Raised when an upstream archive does not meet the safe-inventory requirements."""


def _hash_stream(handle: BinaryIO) -> tuple[str, int]:
    digest = sha256()
    size = 0
    for block in iter(lambda: handle.read(1024 * 1024), b""):
        digest.update(block)
        size += len(block)
    return digest.hexdigest(), size


def _hash_file(path: Path) -> tuple[str, str, int]:
    sha_digest = sha256()
    md5_digest = md5()
    size = 0
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            sha_digest.update(block)
            md5_digest.update(block)
            size += len(block)
    return sha_digest.hexdigest(), md5_digest.hexdigest(), size


def _safe_name(value: str) -> str:
    if not value or "\x00" in value or "\\" in value or value.startswith("/"):
        raise ArchiveInventoryError("archive member has an unsafe name")
    parts = PurePosixPath(value).parts
    if not parts or any(part in {"", ".", ".."} for part in parts):
        raise ArchiveInventoryError("archive member has unsafe path traversal")
    return PurePosixPath(*parts).as_posix()


def inventory(
    archive: Path,
    output: Path,
    *,
    expected_size: int | None,
    expected_md5: str | None,
    maximum_members: int,
    maximum_uncompressed_bytes: int,
) -> dict[str, object]:
    if output.exists() or archive.is_symlink() or not archive.is_file():
        raise ArchiveInventoryError("archive or output contract is invalid")
    archive_sha256, archive_md5, archive_size = _hash_file(archive)
    if expected_size is not None and archive_size != expected_size:
        raise ArchiveInventoryError("archive size differs")
    if expected_md5 is not None and archive_md5 != expected_md5:
        raise ArchiveInventoryError("archive MD5 differs")
    seen: set[str] = set()
    members: list[dict[str, object]] = []
    uncompressed_bytes = 0
    with tarfile.open(archive, mode="r|gz") as handle:
        for member in handle:
            name = _safe_name(member.name)
            if name in seen:
                raise ArchiveInventoryError("archive has duplicate normalized members")
            seen.add(name)
            if len(seen) > maximum_members:
                raise ArchiveInventoryError("archive member ceiling exceeded")
            if member.issym() or member.islnk() or member.isdev() or member.isfifo():
                raise ArchiveInventoryError("archive has a linked or special member")
            if not (member.isdir() or member.isfile()):
                raise ArchiveInventoryError("archive has an unsupported member type")
            record: dict[str, object] = {
                "path": name,
                "type": "directory" if member.isdir() else "regular_file",
                "size_bytes": int(member.size),
                "mode": int(member.mode),
            }
            if member.isfile():
                uncompressed_bytes += int(member.size)
                if uncompressed_bytes > maximum_uncompressed_bytes:
                    raise ArchiveInventoryError("archive uncompressed-byte ceiling exceeded")
                extracted = handle.extractfile(member)
                if extracted is None:
                    raise ArchiveInventoryError("regular archive member is unreadable")
                member_sha256, observed_size = _hash_stream(extracted)
                if observed_size != member.size:
                    raise ArchiveInventoryError("archive member size differs while streaming")
                record["sha256"] = member_sha256
            members.append(record)
    if not members:
        raise ArchiveInventoryError("archive is empty")
    result = {
        "schema_version": "masld-bench-safe-tar-inventory-v1",
        "status": "pass",
        "archive": archive.name,
        "archive_sha256": archive_sha256,
        "archive_md5": archive_md5,
        "archive_size_bytes": archive_size,
        "archive_extracted": False,
        "safe_inventory_passed": True,
        "member_count": len(members),
        "regular_file_count": sum(row["type"] == "regular_file" for row in members),
        "uncompressed_regular_bytes": uncompressed_bytes,
        "checkpoint_index_members": [
            row["path"] for row in members if str(row["path"]).endswith(".index")
        ],
        "checkpoint_data_members": [
            row["path"] for row in members if ".data-" in str(row["path"])
        ],
        "members": members,
    }
    output.write_text(
        json.dumps(result, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--expected-size", type=int)
    parser.add_argument("--expected-md5")
    parser.add_argument("--maximum-members", type=int, default=100_000)
    parser.add_argument(
        "--maximum-uncompressed-bytes", type=int, default=68_719_476_736
    )
    arguments = parser.parse_args()
    result = inventory(
        arguments.archive,
        arguments.output,
        expected_size=arguments.expected_size,
        expected_md5=arguments.expected_md5,
        maximum_members=arguments.maximum_members,
        maximum_uncompressed_bytes=arguments.maximum_uncompressed_bytes,
    )
    print(json.dumps({key: value for key, value in result.items() if key != "members"}, sort_keys=True))


if __name__ == "__main__":
    main()
