#!/usr/bin/env python3
"""Materialize a previously inventoried tar archive without using tar.extract."""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path, PurePosixPath
import tarfile


class SafeMaterializationError(ValueError):
    """Raised when an archive differs from its frozen safe inventory."""


def _sha256_file(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _safe_name(value: str) -> str:
    if not value or value.startswith("/") or "\\" in value or "\x00" in value:
        raise SafeMaterializationError("archive member has an unsafe name")
    parts = PurePosixPath(value).parts
    if not parts or any(part in {"", ".", ".."} for part in parts):
        raise SafeMaterializationError("archive member traverses its destination")
    return PurePosixPath(*parts).as_posix()


def materialize(
    archive: Path,
    inventory_path: Path,
    output: Path,
    *,
    expected_archive_sha256: str,
) -> dict[str, object]:
    if output.exists() or archive.is_symlink() or not archive.is_file():
        raise SafeMaterializationError("materialization request differs")
    if _sha256_file(archive) != expected_archive_sha256:
        raise SafeMaterializationError("archive SHA-256 differs")
    inventory = json.loads(inventory_path.read_text(encoding="utf-8"))
    if (
        inventory.get("schema_version") != "masld-bench-safe-tar-inventory-v1"
        or inventory.get("status") != "pass"
        or inventory.get("safe_inventory_passed") is not True
        or inventory.get("archive_sha256") != expected_archive_sha256
        or inventory.get("archive_extracted") is not False
    ):
        raise SafeMaterializationError("safe archive inventory is not admissible")
    expected_rows = inventory.get("members")
    if not isinstance(expected_rows, list) or not expected_rows:
        raise SafeMaterializationError("safe inventory member roster is missing")
    expected = {str(row["path"]): row for row in expected_rows}
    if len(expected) != len(expected_rows):
        raise SafeMaterializationError("safe inventory contains duplicate members")

    output.mkdir(mode=0o750)
    observed: set[str] = set()
    with tarfile.open(archive, mode="r|gz") as handle:
        for member in handle:
            name = _safe_name(member.name)
            if name in observed or name not in expected:
                raise SafeMaterializationError("archive member roster differs")
            observed.add(name)
            expected_row = expected[name]
            destination = output.joinpath(*PurePosixPath(name).parts)
            if member.isdir():
                if expected_row.get("type") != "directory" or member.size != 0:
                    raise SafeMaterializationError("archive directory differs")
                destination.mkdir(mode=0o750, parents=True, exist_ok=False)
                continue
            if not member.isfile() or expected_row.get("type") != "regular_file":
                raise SafeMaterializationError("archive member type differs")
            if member.size != expected_row.get("size_bytes"):
                raise SafeMaterializationError("archive member size differs")
            destination.parent.mkdir(mode=0o750, parents=True, exist_ok=True)
            source = handle.extractfile(member)
            if source is None:
                raise SafeMaterializationError("archive member is unreadable")
            digest = sha256()
            size = 0
            with destination.open("xb") as target:
                for block in iter(lambda: source.read(1024 * 1024), b""):
                    digest.update(block)
                    size += len(block)
                    target.write(block)
            if size != member.size or digest.hexdigest() != expected_row.get("sha256"):
                raise SafeMaterializationError("archive member digest differs")
            destination.chmod(0o440)
    if observed != set(expected):
        raise SafeMaterializationError("archive member roster is incomplete")
    receipt = {
        "schema_version": "masld-bench-safe-tar-materialization-v1",
        "status": "pass",
        "archive_sha256": expected_archive_sha256,
        "inventory_sha256": _sha256_file(inventory_path),
        "member_count": len(observed),
        "regular_file_count": sum(row.get("type") == "regular_file" for row in expected_rows),
        "materialized_root": output.name,
        "members": expected_rows,
    }
    return receipt


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive", type=Path, required=True)
    parser.add_argument("--inventory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--receipt", type=Path, required=True)
    parser.add_argument("--expected-archive-sha256", required=True)
    arguments = parser.parse_args()
    if arguments.receipt.exists():
        raise SafeMaterializationError("receipt already exists")
    receipt = materialize(
        arguments.archive,
        arguments.inventory,
        arguments.output,
        expected_archive_sha256=arguments.expected_archive_sha256,
    )
    arguments.receipt.write_text(
        json.dumps(receipt, sort_keys=True, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps({key: value for key, value in receipt.items() if key != "members"}, sort_keys=True))


if __name__ == "__main__":
    main()
