#!/usr/bin/env python3
"""Promote the hash-pinned main Figure 4H systems comparison."""

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

CANDIDATE = ROOT / (
    "figures/candidates/histology-continuum-fig4h-systems-20260825T140800Z"
)
SOURCE_PANEL = CANDIDATE / "panels/fig4h_molecular_systems_stage_vs_continuum.pdf"
FIG4 = ROOT / "figures/main/fig4_singlecell_programs"
DESTINATION = FIG4 / "panels/fig4h_molecular_systems_stage_vs_continuum.pdf"
SOURCE_MIRROR = FIG4 / (
    "source_tables/current_candidate/fig4h_system_stage_continuum_20260825"
)
MANIFEST = FIG4 / "manifests/fig4h_system_stage_continuum_promotion_20260825.tsv"
NEW_SHA256 = "c6098c550f496d8afe427ed7eb0a76d73be3f3f0ca268b4da0f683c2df777615"
MIRROR_FILES = (
    "READY_FOR_REVIEW.ok",
    "README.md",
    "figure_manifest.tsv",
    "source_tables/fig4h_system_stage_continuum.tsv",
    "source_tables/fig4h_fixed_label_roster.tsv",
    "provenance/input_manifest.tsv",
    "provenance/audit.tsv",
    "provenance/output_checksums.tsv",
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
        dir=destination.parent, prefix=f".{destination.name}.", suffix=".tmp",
        delete=False
    ) as handle:
        temporary = Path(handle.name)
    try:
        shutil.copy2(source, temporary)
        if sha256(temporary) != sha256(source):
            raise RuntimeError(f"Temporary copy hash mismatch: {temporary}")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main() -> None:
    if not (CANDIDATE / "READY_FOR_REVIEW.ok").is_file():
        raise SystemExit("Candidate does not have READY_FOR_REVIEW.ok")
    if sha256(SOURCE_PANEL) != NEW_SHA256:
        raise SystemExit("Candidate Figure 4H hash mismatch")
    if DESTINATION.exists() or MANIFEST.exists() or SOURCE_MIRROR.exists():
        raise SystemExit("Refusing to overwrite a Figure 4H authority")

    figures = read_tsv(CANDIDATE / "figure_manifest.tsv")
    audits = read_tsv(CANDIDATE / "provenance/audit.tsv")
    if not (
        len(figures) == 1
        and figures[0]["callout"] == "4H"
        and figures[0]["complete_system_family"] == "43/43"
        and figures[0]["strict_maxt_supported"] == "18/43"
        and figures[0]["same_direction"] == "36/43"
        and figures[0]["fixed_outcome_blind_labels"] == "15/43"
        and figures[0]["window_values_used"] == "FALSE"
        and all(row["passed"] == "TRUE" for row in audits)
    ):
        raise SystemExit("Candidate Figure 4H contract failed")

    copy_atomic(SOURCE_PANEL, DESTINATION)
    if sha256(DESTINATION) != NEW_SHA256:
        raise SystemExit("Promoted Figure 4H hash mismatch")

    copied: list[dict[str, str]] = []
    for relative in MIRROR_FILES:
        source = CANDIDATE / relative
        target = SOURCE_MIRROR / relative
        if not source.is_file():
            raise SystemExit(f"Missing candidate provenance file: {source}")
        copy_atomic(source, target)
        copied.append({"relative_path": relative, "sha256": sha256(target)})
    with (SOURCE_MIRROR / "promoted_source_checksums.tsv").open(
        "w", newline=""
    ) as handle:
        writer = csv.DictWriter(
            handle, fieldnames=("relative_path", "sha256"), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(copied)

    row = {
        "callout": "4H",
        "destination": str(DESTINATION.relative_to(ROOT)),
        "source": str(SOURCE_PANEL.relative_to(ROOT)),
        "new_sha256": NEW_SHA256,
        "complete_system_family": "43/43",
        "strict_maxt_supported": "18/43",
        "same_direction": "36/43",
        "spearman_rho": figures[0]["spearman_rho"],
        "fixed_outcome_blind_labels": "15/43",
        "window_values_used": "FALSE",
        "system_geometry_used": "FALSE",
        "statistical_model_changed": "FALSE",
    }
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys(), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)


if __name__ == "__main__":
    main()
