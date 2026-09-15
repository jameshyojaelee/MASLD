#!/usr/bin/env python3
"""Promote the hash-pinned minimal Figure 4F trajectory-only revision."""

from __future__ import annotations

import csv
import hashlib
import os
import shutil
import tempfile
from pathlib import Path


ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", ""))
if not ROOT.is_dir():
    raise SystemExit("MASLD_PROJECT_ROOT is unset or invalid")

FIG4 = ROOT / "figures/main/fig4_singlecell_programs"
CANDIDATE = ROOT / (
    "figures/candidates/histology-continuum-fig4f-minimal-20260825T131700Z"
)
SOURCE_PANEL = CANDIDATE / "panels/fig4f_continuum_program_trajectories.pdf"
DESTINATION = FIG4 / "panels/fig4f_continuum_program_trajectories.pdf"
SOURCE_MIRROR = FIG4 / "source_tables/current_candidate/fig4f_minimal_20260825"
MANIFEST = FIG4 / "manifests/fig4f_minimal_promotion_20260825.tsv"
OLD_SHA256 = "51f8d5ad4d82f112f36a224d93af7710e6914934bdb1e1248cff03efcf3eba30"
NEW_SHA256 = "3813f4db40f469cb7ae4fab70273ae477c913a9e1366d843e861d24ff7736d45"
MIRROR_FILES = (
    "READY_FOR_REVIEW.ok",
    "README.md",
    "figure_manifest.tsv",
    "source_tables/focal_hotspot_five_cohort_windows.tsv",
    "provenance/input_manifest.tsv",
    "provenance/sessionInfo.txt",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def copy_atomic(source: Path, destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=destination.parent, prefix=f".{destination.name}.", suffix=".tmp", delete=False
    ) as handle:
        temporary = Path(handle.name)
    try:
        shutil.copy2(source, temporary)
        if sha256(temporary) != sha256(source):
            raise RuntimeError(f"Temporary copy hash mismatch: {temporary}")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> None:
    if not (CANDIDATE / "READY_FOR_REVIEW.ok").is_file():
        raise SystemExit("Candidate does not have READY_FOR_REVIEW.ok")
    if sha256(SOURCE_PANEL) != NEW_SHA256:
        raise SystemExit("Candidate Figure 4F hash mismatch")
    if sha256(DESTINATION) != OLD_SHA256:
        raise SystemExit("Active Figure 4F preflight hash mismatch")
    if MANIFEST.exists() or SOURCE_MIRROR.exists():
        raise SystemExit("Refusing to overwrite a Figure 4F promotion record")

    copy_atomic(SOURCE_PANEL, DESTINATION)
    if sha256(DESTINATION) != NEW_SHA256:
        raise SystemExit("Active Figure 4F post-promotion hash mismatch")

    copied: list[dict[str, str]] = []
    for relative in MIRROR_FILES:
        source = CANDIDATE / relative
        target = SOURCE_MIRROR / relative
        if not source.is_file():
            raise SystemExit(f"Missing candidate provenance file: {source}")
        copy_atomic(source, target)
        copied.append({"relative_path": relative, "sha256": sha256(target)})
    with (SOURCE_MIRROR / "promoted_source_checksums.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=("relative_path", "sha256"), delimiter="\t")
        writer.writeheader()
        writer.writerows(copied)

    row = {
        "callout": "4F",
        "destination": str(DESTINATION.relative_to(ROOT)),
        "source": str(SOURCE_PANEL.relative_to(ROOT)),
        "old_sha256": OLD_SHA256,
        "new_sha256": NEW_SHA256,
        "trajectory_cohorts_displayed": "5",
        "forest_panel_displayed": "FALSE",
        "embedded_explanatory_text": "FALSE",
        "participant_size_legend_displayed": "FALSE",
        "source_overlap_wording_in_art": "FALSE",
        "window_values_changed": "FALSE",
        "statistical_model_changed": "FALSE",
    }
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys(), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
