#!/usr/bin/env python3
"""Seal/verify exact arm and comparison artifact file sets."""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import os
import re
import stat
from pathlib import Path

from safe_io import publish_new_bytes, require_regular_file


FIELDS = ("relative_path", "size_bytes", "sha256")
SCHEMA = "bg001-analysis-artifact-seal-v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def checked_inside(path: Path, root: Path) -> Path:
    try:
        resolved_parent = path.parent.resolve(strict=True)
        resolved_parent.relative_to(root)
    except (FileNotFoundError, ValueError) as exc:
        raise SystemExit(f"Path escapes run root: {path}") from exc
    return path


def enumerate_artifacts(artifact_root: Path, excluded: set[Path]) -> list[dict[str, str | int]]:
    rows: list[dict[str, str | int]] = []
    for directory, directories, filenames in os.walk(artifact_root, topdown=True, followlinks=False):
        directory_path = Path(directory)
        for name in directories:
            path = directory_path / name
            mode = os.lstat(path).st_mode
            if stat.S_ISLNK(mode) or not stat.S_ISDIR(mode):
                raise SystemExit(f"Artifact tree contains symlink/special directory: {path}")
        for name in filenames:
            path = directory_path / name
            if path in excluded:
                continue
            mode = os.lstat(path).st_mode
            if stat.S_ISLNK(mode) or not stat.S_ISREG(mode):
                raise SystemExit(f"Artifact tree contains symlink/special file: {path}")
            if name.endswith(".tmp") or ".tmp." in name:
                raise SystemExit(f"Artifact tree contains temporary file: {path}")
            rows.append(
                {
                    "relative_path": str(path.relative_to(artifact_root)),
                    "size_bytes": path.stat().st_size,
                    "sha256": sha256(path),
                }
            )
    rows.sort(key=lambda row: str(row["relative_path"]))
    if not rows:
        raise SystemExit(f"Refusing to seal empty artifact tree: {artifact_root}")
    return rows


def encode_manifest(rows: list[dict[str, str | int]]) -> bytes:
    output = io.StringIO(newline="")
    writer = csv.DictWriter(output, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue().encode()


def read_manifest(path: Path) -> list[dict[str, str]]:
    require_regular_file(path)
    payload = path.read_bytes()
    if b"\r" in payload:
        raise SystemExit("Artifact manifest must be UTF-8 LF-only")
    reader = csv.DictReader(io.StringIO(payload.decode("utf-8", errors="strict")), delimiter="\t")
    if tuple(reader.fieldnames or ()) != FIELDS:
        raise SystemExit("Artifact manifest schema/order differs")
    rows = list(reader)
    if not rows or len({row["relative_path"] for row in rows}) != len(rows):
        raise SystemExit("Artifact manifest is empty or contains duplicate paths")
    for row in rows:
        if not re.fullmatch(r"[0-9a-f]{64}", row["sha256"]) or not row["size_bytes"].isdigit():
            raise SystemExit("Artifact manifest contains an invalid size/hash")
        rel = Path(row["relative_path"])
        if rel.is_absolute() or ".." in rel.parts or rel.as_posix() != row["relative_path"]:
            raise SystemExit("Artifact manifest contains unsafe relative path")
    return rows


def bound_hashes(paths: list[Path], run_root: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    for path in paths:
        checked_inside(path, run_root)
        require_regular_file(path)
        relative = str(path.relative_to(run_root))
        if relative in result:
            raise SystemExit(f"Duplicate bound input: {relative}")
        result[relative] = sha256(path)
    return dict(sorted(result.items()))


def seal(args) -> None:
    root = args.run_root.resolve(strict=True)
    if args.artifact_root.is_symlink():
        raise SystemExit(f"Artifact root is symlinked: {args.artifact_root}")
    artifact_root = args.artifact_root.resolve(strict=True)
    artifact_root.relative_to(root)
    manifest = checked_inside(args.manifest, root)
    marker = checked_inside(args.marker, root)
    if manifest.exists() or manifest.is_symlink() or marker.exists() or marker.is_symlink():
        raise SystemExit("Refusing an existing artifact manifest or completion marker")
    excluded = {manifest, marker}
    rows = enumerate_artifacts(artifact_root, excluded)
    publish_new_bytes(manifest, encode_manifest(rows), root)
    marker_payload = {
        "schema": SCHEMA,
        "run_id": root.name,
        "label": args.label,
        "artifact_root": str(artifact_root.relative_to(root)),
        "manifest": str(manifest.relative_to(root)),
        "manifest_sha256": sha256(manifest),
        "artifact_count": len(rows),
        "bound_inputs": bound_hashes(args.bound_input, root),
    }
    publish_new_bytes(
        marker,
        (json.dumps(marker_payload, indent=2, sort_keys=True) + "\n").encode(),
        root,
    )
    print(f"SEALED {args.label}: {len(rows)} artifacts")


def verify(args) -> None:
    root = args.run_root.resolve(strict=True)
    if args.artifact_root.is_symlink():
        raise SystemExit(f"Artifact root is symlinked: {args.artifact_root}")
    manifest = checked_inside(args.manifest, root)
    marker = checked_inside(args.marker, root)
    require_regular_file(marker)
    marker_payload = json.loads(marker.read_text())
    if marker_payload.get("schema") != SCHEMA or marker_payload.get("run_id") != root.name:
        raise SystemExit("Artifact completion marker schema/run differs")
    if marker_payload.get("label") != args.label:
        raise SystemExit("Artifact completion marker label differs")
    artifact_root = root / marker_payload.get("artifact_root", "")
    if artifact_root.resolve(strict=True) != args.artifact_root.resolve(strict=True):
        raise SystemExit("Artifact completion marker root differs")
    if marker_payload.get("manifest") != str(manifest.relative_to(root)):
        raise SystemExit("Artifact completion marker manifest path differs")
    require_regular_file(manifest)
    if marker_payload.get("manifest_sha256") != sha256(manifest):
        raise SystemExit("Artifact manifest hash differs from completion marker")
    expected_rows = read_manifest(manifest)
    actual_rows = enumerate_artifacts(artifact_root, {manifest, marker})
    if actual_rows != [
        {"relative_path": row["relative_path"], "size_bytes": int(row["size_bytes"]), "sha256": row["sha256"]}
        for row in expected_rows
    ]:
        raise SystemExit("Artifact tree file set, size, or hash differs from sealed manifest")
    if marker_payload.get("artifact_count") != len(expected_rows):
        raise SystemExit("Artifact completion marker count differs")
    expected_bound = marker_payload.get("bound_inputs")
    if not isinstance(expected_bound, dict) or not expected_bound:
        raise SystemExit("Artifact completion marker lacks bound inputs")
    actual_bound = bound_hashes([root / relative for relative in expected_bound], root)
    if actual_bound != expected_bound:
        raise SystemExit("Artifact bound input drift")
    print(f"PASS sealed {args.label}: {len(expected_rows)} artifacts")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=("seal", "verify"))
    parser.add_argument("--run-root", required=True, type=Path)
    parser.add_argument("--artifact-root", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--marker", required=True, type=Path)
    parser.add_argument("--label", required=True)
    parser.add_argument("--bound-input", action="append", default=[], type=Path)
    args = parser.parse_args()
    if args.action == "seal" and not args.bound_input:
        raise SystemExit("At least one --bound-input is required")
    if args.action == "seal":
        seal(args)
    else:
        verify(args)


if __name__ == "__main__":
    main()
