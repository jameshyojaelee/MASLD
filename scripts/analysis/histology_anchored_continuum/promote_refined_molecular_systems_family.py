#!/usr/bin/env python3
"""Promote the visually reviewed molecular-systems figure family."""

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
    "histology-continuum-molecular-systems-refined-20260824T224359Z"
)
FIG4 = ROOT / "figures/main/fig4_singlecell_programs"
SOURCE_DESTINATION = FIG4 / (
    "source_tables/current_candidate/"
    "molecular_systems_visual_revision_20260824"
)
PROMOTION_MANIFEST = FIG4 / (
    "manifests/molecular_systems_visual_revision_promotion_20260824.tsv"
)


@dataclass(frozen=True)
class Panel:
    callout: str
    figure_id: str
    source_name: str
    destination: Path
    old_sha256: str
    role: str


PANELS = (
    Panel(
        "4G", "continuum_effect_map",
        "s4_molecular_systems_continuum_effect_map.pdf",
        FIG4 / "panels/fig4g_molecular_systems_continuum_effect_map.pdf",
        "bdf6bfa560cbca91d8094315005994c3f69630e545518a58d0c3d6f0268b1a5d",
        "stage-adjusted continuum effects on fixed membership geometry",
    ),
    Panel(
        "S4Q", "membership_atlas",
        "s4_molecular_systems_membership_atlas.pdf",
        FIG4 / "panels/supplementary/figs4q_molecular_systems_membership_atlas.pdf",
        "c2746b840971da24f59f88d3588c7df5a81973cf68adc8586f5fd86b0283e7b1",
        "signature-excluded outcome-blind membership geometry",
    ),
    Panel(
        "S4R", "fibrosis_stage_map",
        "s4_molecular_systems_fibrosis_stage_effect_map.pdf",
        FIG4 / "panels/supplementary/figs4r_molecular_systems_fibrosis_stage_map.pdf",
        "e23419f15ac6c8d62c189120bf6b21a4e3b830ab420a5a6f3e72d4058c0bbfa5",
        "fibrosis-stage effects on fixed membership geometry",
    ),
    Panel(
        "S4S", "system_stage_continuum_summary",
        "s4_molecular_systems_stage_vs_continuum_summary.pdf",
        FIG4 / "panels/supplementary/figs4s_molecular_systems_stage_vs_continuum.pdf",
        "2596b986dd0c2dcee4409cbe9a039260710227a67936f653390e07500a0e92a6",
        "descriptive system-level stage-versus-continuum alignment",
    ),
    Panel(
        "S4T", "five_cohort_sliding_windows",
        "s4_molecular_systems_five_cohort_sliding_windows.pdf",
        FIG4 / "panels/supplementary/figs4t_molecular_systems_continuum_windows.pdf",
        "7f966e7d96b38e67904e6e6e0d49311ce623aca7b1eaee195225b002cde3260c",
        "five-cohort descriptive windows in pooled late-early order",
    ),
    Panel(
        "S4U", "pooled_continuum_windows",
        "s4_molecular_systems_pooled_continuum_windows.pdf",
        FIG4 / "panels/supplementary/figs4u_molecular_systems_pooled_continuum_windows.pdf",
        "854e7744c88a38dea7ae943ecbb46653e3cd54fd24658090a2e31a04fdc612c6",
        "equal-cohort pooled windows in descriptive late-early order",
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


def stage_copy(source: Path, destination: Path) -> Path:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(
        dir=destination.parent,
        prefix=f".{destination.name}.",
        suffix=".tmp",
        delete=False,
    ) as handle:
        temporary = Path(handle.name)
    shutil.copy2(source, temporary)
    if sha256(temporary) != sha256(source):
        temporary.unlink(missing_ok=True)
        raise SystemExit(f"Staged panel hash mismatch: {source}")
    return temporary


def main() -> None:
    if not (SOURCE / "READY_FOR_REVIEW.ok").is_file():
        raise SystemExit("Refined molecular-systems candidate is not review-ready")
    audit_rows = read_tsv(SOURCE / "provenance/audit.tsv")
    if not audit_rows or any(row["passed"].upper() != "TRUE" for row in audit_rows):
        raise SystemExit("Refined molecular-systems candidate audit did not pass")
    figure_rows = {
        row["figure_id"]: row for row in read_tsv(SOURCE / "figure_manifest.tsv")
    }
    if set(figure_rows) != {item.figure_id for item in PANELS}:
        raise SystemExit("Refined six-panel figure manifest drift")
    for row in read_tsv(SOURCE / "provenance/output_checksums.tsv"):
        path = Path(row["path"])
        if not path.is_file() or sha256(path) != row["sha256"]:
            raise SystemExit(f"Refined candidate output hash drift: {path}")

    if SOURCE_DESTINATION.exists() or PROMOTION_MANIFEST.exists():
        raise SystemExit("Refined molecular-systems promotion namespace already exists")

    staged: list[tuple[Path, Path]] = []
    try:
        for item in PANELS:
            row = figure_rows[item.figure_id]
            source = SOURCE / "panels" / item.source_name
            if Path(row["path"]).resolve() != source.resolve():
                raise SystemExit(f"Candidate manifest path mismatch: {item.figure_id}")
            if row["promoted"].upper() != "FALSE":
                raise SystemExit(f"Candidate promotion-state drift: {item.figure_id}")
            if not source.is_file() or sha256(source) != row["sha256"]:
                raise SystemExit(f"Candidate panel hash mismatch: {source}")
            if not item.destination.is_file() or sha256(item.destination) != item.old_sha256:
                raise SystemExit(f"Active panel hash drift: {item.destination}")
            staged.append((stage_copy(source, item.destination), item.destination))

        for temporary, destination in staged:
            os.replace(temporary, destination)
    finally:
        for temporary, _ in staged:
            temporary.unlink(missing_ok=True)

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
            raise SystemExit(f"Promoted visual source hash mismatch: {destination}")

    panel_rows: list[dict[str, str]] = []
    for item in PANELS:
        row = figure_rows[item.figure_id]
        panel_rows.append({
            "callout": item.callout,
            "destination": str(item.destination.relative_to(ROOT)),
            "source": str((SOURCE / "panels" / item.source_name).relative_to(ROOT)),
            "old_sha256": item.old_sha256,
            "new_sha256": sha256(item.destination),
            "role": item.role,
            "display_order_rule": row["display_order_rule"],
            "result_adaptive_display_order": row["result_adaptive_display_order"],
            "data_values_changed": row["data_values_changed"],
            "geometry_changed": row["geometry_changed"],
            "revision": "compact_layout_consistent_labels_semantic_direction_colors",
        })

    PROMOTION_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with PROMOTION_MANIFEST.open("x", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=panel_rows[0].keys(), delimiter="\t"
        )
        writer.writeheader()
        writer.writerows(panel_rows)

    checksum_path = SOURCE_DESTINATION / "promoted_source_checksums.tsv"
    with checksum_path.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("path", "sha256"))
        for path in sorted(copied):
            writer.writerow((str(path.relative_to(ROOT)), sha256(path)))
    (SOURCE_DESTINATION / "PROMOTED.ok").write_text(
        "visually_reviewed_six_panel_revision\n"
    )


if __name__ == "__main__":
    main()
