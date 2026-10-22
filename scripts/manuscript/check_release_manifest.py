#!/usr/bin/env python3
"""Check a named release against its manifest.

Every manifest row must exist with its recorded sha256, and every main-figure
panel in the release must have its source table in the release.
"""

from __future__ import annotations

import hashlib
import os
import sys
from pathlib import Path

import pandas as pd


ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas",
))
RELEASE_ID = os.environ.get("MANUSCRIPT_RELEASE_ID", "2026-10-21-r9")
OUT = ROOT / "RNA-seq/results/manuscript_release" / RELEASE_ID
PANEL_PLAN = ROOT / "figures/candidates/panel-deep-fixes-20261017T035633Z/INSTALL_PLAN.tsv"
SOURCE_SUBDIR = "data/panel_deep_20261017"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    errors: list[str] = []
    manifest = pd.read_csv(OUT / "release_manifest.tsv", sep="\t", dtype=str)
    if set(manifest["release_id"]) != {RELEASE_ID}:
        errors.append(f"manifest names releases {sorted(set(manifest['release_id']))}")

    for row in manifest.itertuples():
        path = ROOT / row.path
        if not path.exists():
            errors.append(f"missing: {row.path}")
        elif sha256(path) != row.sha256:
            errors.append(f"changed since release: {row.path}")

    listed = set(manifest["path"])
    plan = pd.read_csv(PANEL_PLAN, sep="\t", dtype=str, keep_default_na=False)
    redrawn = plan[plan["status"] == "redrawn"]
    for row in redrawn.itertuples():
        if row.install_target not in listed:
            errors.append(f"panel not in release: {row.panel} {row.install_target}")
        table = row.source_table.split(" (+")[0]
        expected = str(Path(row.install_target).parent / SOURCE_SUBDIR / Path(table).name)
        if row.panel.startswith("7") and not (ROOT / expected).exists():
            expected = str(Path(row.install_target).parent / SOURCE_SUBDIR / f"fig{row.panel}.tsv")
        if expected not in listed:
            errors.append(f"source table not in release: {row.panel} {expected}")

    print(f"RELEASE\t{RELEASE_ID}")
    print(f"MANIFEST_ROWS\t{len(manifest)}")
    print(f"MAIN_PANELS_CHECKED\t{len(redrawn)}")
    if errors:
        print("RELEASE_MANIFEST_CHECK\tFAIL")
        for message in errors:
            print(f"ERROR\t{message}")
        sys.exit(1)
    print("RELEASE_MANIFEST_CHECK\tPASS")


if __name__ == "__main__":
    main()
