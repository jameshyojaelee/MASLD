#!/usr/bin/env python3
"""Publish the fresh baseline manifests and immutable transaction marker."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import stat
from pathlib import Path

from safe_io import publish_new_bytes, require_regular_file
from verify_analysis_contract import EXPECTED_BASELINE_INPUTS
from bam_manifest_finalization_guard import GuardError, verify_artifacts


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def regular_tree_files(root: Path) -> list[Path]:
    files: list[Path] = []
    for directory, directories, filenames in os.walk(root, topdown=True, followlinks=False):
        directory_path = Path(directory)
        for name in directories:
            path = directory_path / name
            mode = os.lstat(path).st_mode
            if stat.S_ISLNK(mode) or not stat.S_ISDIR(mode):
                raise SystemExit(f"Baseline tree contains symlink/special directory: {path}")
        for name in filenames:
            path = directory_path / name
            mode = os.lstat(path).st_mode
            if stat.S_ISLNK(mode) or not stat.S_ISREG(mode):
                raise SystemExit(f"Baseline tree contains symlink/special file: {path}")
            if name.endswith(".tmp") or ".tmp." in name:
                raise SystemExit(f"Baseline tree contains a temporary artifact: {path}")
            files.append(path)
    return sorted(files)


def sha_lines(paths: list[Path]) -> bytes:
    return "".join(f"{sha256(path)}  {path}\n" for path in paths).encode()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=Path)
    args = parser.parse_args()
    root = args.run_root.resolve(strict=True)
    frozen_root = root / "frozen_sets"
    if not frozen_root.is_dir():
        raise SystemExit("Missing frozen_sets directory")
    files_manifest = root / "manifests/baseline_frozen_files.sha256"
    inputs_manifest = root / "manifests/baseline_manifests.sha256"
    marker_path = root / "BASELINE_FROZEN.json"
    compatibility_path = root / "BASELINE_FROZEN"
    for target in (files_manifest, inputs_manifest, marker_path, compatibility_path):
        if target.exists() or target.is_symlink():
            raise SystemExit(f"Refusing existing baseline transaction target: {target}")

    try:
        verify_artifacts(root)
    except GuardError as exc:
        raise SystemExit(f"Cannot seal an invalid BAM-manifest guard transaction: {exc}") from exc

    frozen_files = regular_tree_files(frozen_root)
    if not frozen_files:
        raise SystemExit("Cannot seal an empty frozen baseline")
    input_files = [root / relative for relative in EXPECTED_BASELINE_INPUTS]
    for path in input_files:
        require_regular_file(path)
    publish_new_bytes(files_manifest, sha_lines(frozen_files), root)
    publish_new_bytes(inputs_manifest, sha_lines(input_files), root)
    marker = {
        "schema": "bg001-baseline-seal-v1",
        "run_id": root.name,
        "baseline_frozen_files_manifest": str(files_manifest.relative_to(root)),
        "baseline_frozen_files_manifest_sha256": sha256(files_manifest),
        "baseline_manifests_manifest": str(inputs_manifest.relative_to(root)),
        "baseline_manifests_manifest_sha256": sha256(inputs_manifest),
        "frozen_file_count": len(frozen_files),
        "bound_input_count": len(input_files),
        "run_contract_sha256": sha256(root / "contract/run_contract.json"),
        "source_manifest_sha256": sha256(root / "contract/source_manifest.tsv"),
        "compatibility_marker": str(compatibility_path.relative_to(root)),
        "compatibility_marker_contract": "lowercase_sha256(BASELINE_FROZEN.json)+LF",
    }
    publish_new_bytes(
        marker_path,
        (json.dumps(marker, indent=2, sort_keys=True) + "\n").encode(),
        root,
    )
    # Frozen recount launchers predate the structured transaction and only test
    # for BASELINE_FROZEN.  Publish a nonempty compatibility marker whose bytes
    # cryptographically bind the structured marker; never recreate the former
    # unauthenticated zero-byte sentinel for a fresh run.
    publish_new_bytes(
        compatibility_path,
        (sha256(marker_path) + "\n").encode(),
        root,
    )
    print(f"SEALED baseline: {len(frozen_files)} files, {len(input_files)} bound inputs")


if __name__ == "__main__":
    main()
