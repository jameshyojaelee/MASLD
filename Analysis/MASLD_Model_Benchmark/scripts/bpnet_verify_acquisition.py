#!/usr/bin/env python3
"""Verify the immutable legacy BPNet acquisition receipt without deserializing it."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path, PurePosixPath


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def verify_acquisition(root: Path) -> dict[str, object]:
    root = root.resolve(strict=True)
    manifest_path = root / "ARTIFACTS.json"
    complete_path = root / "COMPLETE"
    if manifest_path.is_symlink() or complete_path.is_symlink():
        raise ValueError("BPNet acquisition control files may not be symlinks")
    if not manifest_path.is_file() or not complete_path.is_file():
        raise ValueError("BPNet acquisition control files are missing")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if set(manifest) != {"schema_version", "metadata", "artifacts"}:
        raise ValueError("BPNet acquisition manifest fields differ")
    if manifest["schema_version"] != "masld-bench-artifacts-v1":
        raise ValueError("BPNet acquisition manifest schema differs")
    if manifest["metadata"].get("artifact_class") != "bpnet_source_acquisition":
        raise ValueError("BPNet acquisition artifact class differs")
    if manifest["metadata"].get("model_id") != "bpnet":
        raise ValueError("BPNet acquisition model identity differs")

    complete: dict[str, str] = {}
    for line in complete_path.read_text(encoding="utf-8").splitlines():
        key, value = line.split("\t", 1)
        if key in complete:
            raise ValueError("duplicate BPNet acquisition completion field")
        complete[key] = value
    manifest_sha256 = _sha256(manifest_path)
    if complete != {
        "status": "complete",
        "artifact_class": "bpnet_source_acquisition",
        "artifacts_sha256": manifest_sha256,
    }:
        raise ValueError("BPNet acquisition completion marker differs")

    expected: set[str] = set()
    for item in manifest["artifacts"]:
        if set(item) != {"path", "sha256", "size_bytes"}:
            raise ValueError("BPNet acquisition artifact fields differ")
        relative = PurePosixPath(str(item["path"]))
        if relative.is_absolute() or not relative.parts or ".." in relative.parts:
            raise ValueError(f"unsafe BPNet acquisition artifact path: {relative}")
        relative_text = relative.as_posix()
        if relative_text in expected:
            raise ValueError(f"duplicate BPNet acquisition artifact: {relative}")
        expected.add(relative_text)
        path = root.joinpath(*relative.parts)
        if not path.is_file() or path.is_symlink():
            raise ValueError(f"missing or unsafe BPNet acquisition artifact: {relative}")
        if path.stat().st_size != item["size_bytes"]:
            raise ValueError(f"BPNet acquisition artifact size differs: {relative}")
        if _sha256(path) != item["sha256"]:
            raise ValueError(f"BPNet acquisition artifact digest differs: {relative}")

    observed = {
        path.relative_to(root).as_posix()
        for path in root.rglob("*")
        if path.is_file() and path.name not in {"ARTIFACTS.json", "COMPLETE"}
    }
    if observed != expected:
        raise ValueError("BPNet acquisition artifact roster differs")
    return {
        "status": "pass",
        "artifact_count": len(expected),
        "manifest_sha256": manifest_sha256,
    }


def _write_json_exclusive(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o640)
    with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--acquisition", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    _write_json_exclusive(args.output, verify_acquisition(args.acquisition))


if __name__ == "__main__":
    main()
