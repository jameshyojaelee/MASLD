#!/usr/bin/env python3
"""Promote the hash-verified equal-cohort pooled molecular-systems panel."""

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

SOURCE = ROOT / (
    "figures/candidates/"
    "histology-continuum-molecular-systems-pooled-20260824T202331Z"
)
FIG4 = ROOT / "figures/main/fig4_singlecell_programs"
DESTINATION = (
    FIG4
    / "panels/supplementary/figs4u_molecular_systems_pooled_continuum_windows.pdf"
)
SOURCE_DESTINATION = (
    FIG4
    / "source_tables/current_candidate/molecular_systems_pooled_20260824"
)
PROMOTION_MANIFEST = (
    FIG4 / "manifests/molecular_systems_pooled_promotion_20260824.tsv"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def atomic_copy(source: Path, destination: Path) -> None:
    if destination.exists():
        raise SystemExit(f"Refusing to overwrite promoted artifact: {destination}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
    try:
        shutil.copy2(source, temporary)
        if sha256(temporary) != sha256(source):
            raise SystemExit(f"Temporary copy hash mismatch: {source}")
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)


def main() -> None:
    if not (SOURCE / "READY_FOR_REVIEW.ok").is_file():
        raise SystemExit("Pooled molecular-systems candidate is not review-ready")
    audit_rows = read_tsv(SOURCE / "provenance/audit.tsv")
    if not audit_rows or any(row["passed"].upper() != "TRUE" for row in audit_rows):
        raise SystemExit("Pooled molecular-systems audit did not pass")

    figure_rows = read_tsv(SOURCE / "figure_manifest.tsv")
    if len(figure_rows) != 1 or figure_rows[0]["figure_id"] != "pooled_continuum_windows":
        raise SystemExit("Pooled figure manifest drift")
    source_panel = Path(figure_rows[0]["path"])
    if not source_panel.is_file() or sha256(source_panel) != figure_rows[0]["sha256"]:
        raise SystemExit("Pooled panel hash mismatch")
    if figure_rows[0]["pooling_weight"] != "equal_cohort":
        raise SystemExit("Pooled panel does not use equal cohort weights")

    output_rows = read_tsv(SOURCE / "provenance/output_checksums.tsv")
    for row in output_rows:
        path = Path(row["path"])
        if not path.is_file() or sha256(path) != row["sha256"]:
            raise SystemExit(f"Candidate output hash drift: {path}")

    if SOURCE_DESTINATION.exists():
        raise SystemExit(f"Source-table namespace exists: {SOURCE_DESTINATION}")
    if PROMOTION_MANIFEST.exists():
        raise SystemExit(f"Promotion manifest exists: {PROMOTION_MANIFEST}")

    atomic_copy(source_panel, DESTINATION)
    SOURCE_DESTINATION.mkdir(parents=True)
    copied: list[Path] = []
    for relative in (
        "figure_manifest.tsv",
        "READY_FOR_REVIEW.ok",
        "source_tables",
        "provenance",
    ):
        source = SOURCE / relative
        destination = SOURCE_DESTINATION / relative
        if source.is_dir():
            shutil.copytree(source, destination)
            copied.extend(path for path in destination.rglob("*") if path.is_file())
        else:
            shutil.copy2(source, destination)
            copied.append(destination)

    for destination in copied:
        source = SOURCE / destination.relative_to(SOURCE_DESTINATION)
        if sha256(destination) != sha256(source):
            raise SystemExit(f"Promoted source hash mismatch: {destination}")

    row = {
        "callout": "S4U",
        "destination": str(DESTINATION.relative_to(ROOT)),
        "source": str(source_panel.relative_to(ROOT)),
        "sha256": sha256(DESTINATION),
        "role": "equal-cohort pooled molecular-system continuum windows",
        "scope": "five_cohort_descriptive_equal_weight",
    }
    PROMOTION_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with PROMOTION_MANIFEST.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=row.keys(), delimiter="\t")
        writer.writeheader()
        writer.writerow(row)

    checksum_path = SOURCE_DESTINATION / "promoted_source_checksums.tsv"
    with checksum_path.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("path", "sha256"))
        for path in sorted(copied):
            writer.writerow((str(path.relative_to(ROOT)), sha256(path)))


if __name__ == "__main__":
    main()
