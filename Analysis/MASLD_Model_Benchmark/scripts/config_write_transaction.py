#!/usr/bin/env python3
"""Make a config-writing run an actual transaction.

A run that edits config and then dies leaves the edits in place. The execution
directory gets renamed to ``.failed``, so anyone reading job status concludes
nothing landed while the tree has in fact moved. The guard fires and the damage
escapes through the failure path rather than through the check.

``snapshot`` copies every file a run intends to write, with its pre-write
digest, into the staging tree. ``restore`` puts them back and proves each one
re-derives to the digest recorded before the run touched it, then writes a
receipt saying positively what was undone.

Neither mode reasons about time. A restore is decided by the caller's exit
status and verified by digest, never by comparing modification times.
"""

from __future__ import annotations

import argparse
from hashlib import sha256
import json
from pathlib import Path
import shutil


class TransactionError(RuntimeError):
    """Raised when a snapshot or a restore cannot be trusted."""


def digest(path: Path) -> str:
    value = sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def snapshot(root: Path, targets: list[str], into: Path) -> dict:
    """Copy each target and record the digest it had before any write."""

    if into.exists():
        raise TransactionError(f"snapshot directory exists: {into}")
    into.mkdir(parents=True)
    entries = []
    for relative in targets:
        source = root / relative
        if not source.is_file():
            raise TransactionError(f"cannot snapshot a missing file: {source}")
        stored = into / relative.replace("/", "__")
        shutil.copy2(source, stored)
        entries.append({
            "path": relative,
            "stored_as": stored.name,
            "sha256_before": digest(source),
            "size_before": source.stat().st_size,
        })
    manifest = {
        "schema_version": "masld-bench-config-write-transaction-v1",
        "files": entries,
        "restored": False,
    }
    (into / "MANIFEST.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def restore(root: Path, into: Path) -> dict:
    """Put every snapshotted file back and prove it matches its pre-write digest."""

    manifest_path = into / "MANIFEST.json"
    if not manifest_path.is_file():
        raise TransactionError(f"no snapshot manifest at {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    results = []
    failures = []
    for entry in manifest["files"]:
        target = root / entry["path"]
        stored = into / entry["stored_as"]
        if not stored.is_file():
            failures.append(f"{entry['path']}: snapshot copy missing")
            continue
        changed = target.is_file() and digest(target) != entry["sha256_before"]
        shutil.copy2(stored, target)
        actual = digest(target)
        if actual != entry["sha256_before"]:
            failures.append(
                f"{entry['path']}: restored to {actual}, expected {entry['sha256_before']}"
            )
        results.append({
            "path": entry["path"],
            "was_modified_by_the_run": changed,
            "restored_to_sha256": actual,
            "verified": actual == entry["sha256_before"],
        })
    manifest["restored"] = True
    manifest["restore_results"] = results
    manifest["files_the_run_had_modified"] = sum(
        1 for r in results if r["was_modified_by_the_run"]
    )
    manifest["all_verified"] = not failures
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if failures:
        raise TransactionError("restore could not be verified:\n" + "\n".join(failures))
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("snapshot", "restore"))
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--into", type=Path, required=True)
    parser.add_argument("--target", action="append", default=[])
    args = parser.parse_args()
    if args.mode == "snapshot":
        result = snapshot(args.root, args.target, args.into)
        print(f"snapshotted {len(result['files'])} files into {args.into}")
    else:
        result = restore(args.root, args.into)
        print(
            f"restored {len(result['restore_results'])} files; "
            f"{result['files_the_run_had_modified']} had been modified by the run; "
            f"all verified: {result['all_verified']}"
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
