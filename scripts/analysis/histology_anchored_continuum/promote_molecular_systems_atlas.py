#!/usr/bin/env python3
"""Promote the corrected, hash-verified molecular-systems atlas panels."""

from __future__ import annotations

import csv
import hashlib
import os
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(os.environ.get("MASLD_PROJECT_ROOT", ""))
if not ROOT.is_dir():
    raise SystemExit("MASLD_PROJECT_ROOT is unset or invalid")

SOURCE = ROOT / (
    "figures/candidates/"
    "histology-continuum-molecular-systems-atlas-20260824T195148Z"
)
FIG4 = ROOT / "figures/main/fig4_singlecell_programs"
SOURCE_DESTINATION = (
    FIG4 / "source_tables/current_candidate/molecular_systems_atlas_20260824"
)
PROMOTION_MANIFEST = (
    FIG4 / "manifests/molecular_systems_atlas_promotion_20260824.tsv"
)


@dataclass(frozen=True)
class Panel:
    callout: str
    figure_id: str
    source_name: str
    destination: Path
    role: str


PANELS = (
    Panel(
        "4G", "continuum_effect_map",
        "s4_molecular_systems_continuum_effect_map.pdf",
        FIG4 / "panels/fig4g_molecular_systems_continuum_effect_map.pdf",
        "stage-adjusted continuum effects on outcome-blind membership systems",
    ),
    Panel(
        "S4Q", "membership_atlas",
        "s4_molecular_systems_membership_atlas.pdf",
        FIG4 / "panels/supplementary/figs4q_molecular_systems_membership_atlas.pdf",
        "signature-excluded outcome-blind membership geometry",
    ),
    Panel(
        "S4R", "fibrosis_stage_map",
        "s4_molecular_systems_fibrosis_stage_effect_map.pdf",
        FIG4 / "panels/supplementary/figs4r_molecular_systems_fibrosis_stage_map.pdf",
        "fibrosis-stage effects on the same membership geometry",
    ),
    Panel(
        "S4S", "system_stage_continuum_summary",
        "s4_molecular_systems_stage_vs_continuum_summary.pdf",
        FIG4 / "panels/supplementary/figs4s_molecular_systems_stage_vs_continuum.pdf",
        "descriptive system-level stage-versus-continuum alignment",
    ),
    Panel(
        "S4T", "five_cohort_sliding_windows",
        "s4_molecular_systems_five_cohort_sliding_windows.pdf",
        FIG4 / "panels/supplementary/figs4t_molecular_systems_continuum_windows.pdf",
        "five-cohort descriptive molecular-system sliding windows",
    ),
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


def require_passed_audit(path: Path) -> None:
    rows = read_tsv(path)
    if not rows or any(row.get("passed", "").upper() != "TRUE" for row in rows):
        raise SystemExit(f"Audit did not pass: {path}")


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
    if SOURCE.name != "histology-continuum-molecular-systems-atlas-20260824T195148Z":
        raise SystemExit("Unexpected candidate identity")
    if not (SOURCE / "READY_FOR_REVIEW.ok").is_file():
        raise SystemExit("Corrected candidate has not passed its review gates")
    require_passed_audit(SOURCE / "provenance/outcome_blind_audit.tsv")
    require_passed_audit(SOURCE / "provenance/overlay_audit.tsv")

    figure_rows = {row["figure_id"]: row for row in read_tsv(SOURCE / "figure_manifest.tsv")}
    for item in PANELS:
        row = figure_rows.get(item.figure_id)
        source = SOURCE / "panels" / item.source_name
        if row is None or Path(row["path"]).resolve() != source.resolve():
            raise SystemExit(f"Figure-manifest path mismatch: {item.figure_id}")
        if row.get("promoted", "").upper() != "FALSE":
            raise SystemExit(f"Candidate promotion-state mismatch: {item.figure_id}")
        if not source.is_file() or sha256(source) != row["sha256"]:
            raise SystemExit(f"Candidate panel hash mismatch: {source}")
        if item.destination.exists():
            raise SystemExit(f"Destination already exists: {item.destination}")

    if SOURCE_DESTINATION.exists():
        raise SystemExit(f"Source-table namespace already exists: {SOURCE_DESTINATION}")
    if PROMOTION_MANIFEST.exists():
        raise SystemExit(f"Promotion manifest already exists: {PROMOTION_MANIFEST}")

    panel_rows: list[dict[str, str]] = []
    for item in PANELS:
        source = SOURCE / "panels" / item.source_name
        atomic_copy(source, item.destination)
        panel_rows.append({
            "callout": item.callout,
            "destination": str(item.destination.relative_to(ROOT)),
            "source": str(source.relative_to(ROOT)),
            "sha256": sha256(item.destination),
            "role": item.role,
            "inference_scope": (
                "GSE162694;GSE213621_stage_sex_adjusted"
                if item.callout in {"4G", "S4R", "S4S"}
                else "outcome_blind_or_five_cohort_descriptive"
            ),
        })

    SOURCE_DESTINATION.mkdir(parents=True)
    copied_source_files: list[Path] = []
    for relative in (
        "figure_manifest.tsv",
        "READY_FOR_REVIEW.ok",
        "geometry",
        "overlays",
        "provenance",
        "source_tables",
    ):
        source = SOURCE / relative
        destination = SOURCE_DESTINATION / relative
        if source.is_dir():
            shutil.copytree(source, destination)
            copied_source_files.extend(path for path in destination.rglob("*") if path.is_file())
        else:
            shutil.copy2(source, destination)
            copied_source_files.append(destination)

    for destination in copied_source_files:
        source = SOURCE / destination.relative_to(SOURCE_DESTINATION)
        if sha256(destination) != sha256(source):
            raise SystemExit(f"Promoted source-table hash mismatch: {destination}")

    PROMOTION_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with PROMOTION_MANIFEST.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=panel_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(panel_rows)

    source_checksums = SOURCE_DESTINATION / "promoted_source_checksums.tsv"
    with source_checksums.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("path", "sha256"))
        for path in sorted(copied_source_files):
            writer.writerow((str(path.relative_to(ROOT)), sha256(path)))


if __name__ == "__main__":
    main()
