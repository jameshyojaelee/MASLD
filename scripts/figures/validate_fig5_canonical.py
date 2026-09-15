#!/usr/bin/env python3
"""Read-only final checks for the promoted canonical Figure 5 package."""

from __future__ import annotations

import csv
import hashlib
import subprocess
from pathlib import Path

from PyPDF2 import PdfReader


ROOT = Path(__file__).resolve().parents[2]
FIGURE = ROOT / "figures/main/fig5_molecular_context"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def sha256(path: Path) -> str:
    value = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            value.update(block)
    return value.hexdigest()


def main() -> None:
    manifest = read_tsv(FIGURE / "CANONICAL_MAIN_PANELS.tsv")
    if [row["callout"] for row in manifest] != [f"5{letter}" for letter in "ABCDEF"]:
        raise RuntimeError("canonical manifest is not Figure 5A-F")
    for row in manifest:
        path = FIGURE / row["filename"]
        if sha256(path) != row["sha256"] or path.stat().st_size != int(row["bytes"]):
            raise RuntimeError(f"canonical hash or size mismatch: {row['callout']}")
        if len(PdfReader(str(path)).pages) != 1 or b"/Subtype /Type3" in path.read_bytes():
            raise RuntimeError(f"canonical PDF structural failure: {row['callout']}")
    assembled_manifest = read_tsv(FIGURE / "ASSEMBLED_FIGURE.tsv")
    if len(assembled_manifest) != 1 or not assembled_manifest[0]["status"].startswith("promoted_and_postcopy_validated"):
        raise RuntimeError("assembled Figure 5 manifest is not promoted")
    assembled = FIGURE / assembled_manifest[0]["filename"]
    if sha256(assembled) != assembled_manifest[0]["sha256"] or assembled.stat().st_size != int(assembled_manifest[0]["bytes"]):
        raise RuntimeError("assembled Figure 5 hash or size mismatch")
    if len(PdfReader(str(assembled)).pages) != 1 or b"/Subtype /Type3" in assembled.read_bytes():
        raise RuntimeError("assembled Figure 5 structural failure")
    source_d = read_tsv(FIGURE / "panels/data/fig5d_variant_consistent_accessibility.tsv")
    if len(source_d) != 50 or sum(int(row["n_pairs"]) for row in source_d if row["lineage"] == "hepatocyte") != 816:
        raise RuntimeError("canonical Figure 5D source grid is invalid")
    if (FIGURE / "data/fig5d_snatac_accessibility_source.tsv").exists():
        raise RuntimeError("blocked pre-variant-consistent Figure 5D source remains active")
    retired_panel_root = ROOT / "figures/main" / ("fig4_" + "validation") / "panels"
    if retired_panel_root.exists():
        raise RuntimeError("retired fig4_validation panel path was recreated")
    if not (FIGURE / "_legacy/2026-08-24_pre_fig5_completion").is_dir():
        raise RuntimeError("pre-promotion Figure 5 archive is absent")
    active_text = "\n".join(
        path.read_text(errors="replace")
        for path in (
            ROOT / "docs/STATUS.md",
            ROOT / "docs/RESULTS.md",
            ROOT / "docs/ROADMAP.md",
            FIGURE / "README.md",
        )
    )
    stale = (
        "5D remains blocked",
        "Figure 5D remain unbuilt",
        "replay itself has not been run",
        "Panel 5A is scientifically accepted, but",
    )
    if any(term in active_text for term in stale):
        raise RuntimeError("active authority retains a stale Figure 5 status")
    subprocess.run(
        ["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", str(assembled)],
        check=True,
    )
    print("FIGURE_5_CANONICAL_VALIDATION\tPASS")


if __name__ == "__main__":
    main()
