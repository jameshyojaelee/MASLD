#!/usr/bin/env python3
"""Promote the aspect-matched narrow pooled molecular-systems heatmap."""

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
    "histology-continuum-molecular-systems-refined-20260824T225846Z"
)
FIG4 = ROOT / "figures/main/fig4_singlecell_programs"
DESTINATION = FIG4 / (
    "panels/supplementary/figs4u_molecular_systems_pooled_continuum_windows.pdf"
)
SOURCE_DESTINATION = FIG4 / (
    "source_tables/current_candidate/molecular_systems_pooled_narrow_20260824"
)
PROMOTION_MANIFEST = FIG4 / (
    "manifests/molecular_systems_pooled_narrow_promotion_20260824.tsv"
)
EXPECTED_OLD_SHA256 = (
    "46ba62c3d5d58f01e404debab4a5acc1ce4eeea3e82720e552e41efc09dfad50"
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


def replace_verified(source: Path, destination: Path) -> None:
    if not destination.is_file() or sha256(destination) != EXPECTED_OLD_SHA256:
        raise SystemExit(f"Existing S4U does not match the visual-revision hash: {destination}")
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
        raise SystemExit("Narrow pooled candidate is not review-ready")
    audit_rows = read_tsv(SOURCE / "provenance/audit.tsv")
    if not audit_rows or any(row["passed"].upper() != "TRUE" for row in audit_rows):
        raise SystemExit("Narrow pooled candidate audit did not pass")
    figure_rows = {
        row["figure_id"]: row for row in read_tsv(SOURCE / "figure_manifest.tsv")
    }
    row = figure_rows.get("pooled_continuum_windows")
    source_panel = SOURCE / "panels/s4_molecular_systems_pooled_continuum_windows.pdf"
    if row is None or Path(row["path"]).resolve() != source_panel.resolve():
        raise SystemExit("Narrow pooled figure-manifest path drift")
    if row["promoted"].upper() != "FALSE":
        raise SystemExit("Narrow pooled candidate promotion-state drift")
    if not source_panel.is_file() or sha256(source_panel) != row["sha256"]:
        raise SystemExit("Narrow pooled candidate panel hash mismatch")
    for output in read_tsv(SOURCE / "provenance/output_checksums.tsv"):
        path = Path(output["path"])
        if not path.is_file() or sha256(path) != output["sha256"]:
            raise SystemExit(f"Narrow pooled candidate output hash drift: {path}")
    if SOURCE_DESTINATION.exists() or PROMOTION_MANIFEST.exists():
        raise SystemExit("Narrow pooled promotion namespace already exists")

    replace_verified(source_panel, DESTINATION)
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
            raise SystemExit(f"Promoted narrow source hash mismatch: {destination}")

    promotion_row = {
        "callout": "S4U",
        "destination": str(DESTINATION.relative_to(ROOT)),
        "source": str(source_panel.relative_to(ROOT)),
        "old_sha256": EXPECTED_OLD_SHA256,
        "new_sha256": sha256(DESTINATION),
        "role": "equal-cohort pooled molecular-system continuum windows",
        "display_order_rule": row["display_order_rule"],
        "data_values_changed": row["data_values_changed"],
        "geometry_changed": row["geometry_changed"],
        "revision": "coord_fixed_ratio_1_canvas_2.25_inches_match_s4t_tile_aspect",
    }
    PROMOTION_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    with PROMOTION_MANIFEST.open("x", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=promotion_row.keys(), delimiter="\t"
        )
        writer.writeheader()
        writer.writerow(promotion_row)

    checksum_path = SOURCE_DESTINATION / "promoted_source_checksums.tsv"
    with checksum_path.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("path", "sha256"))
        for path in sorted(copied):
            writer.writerow((str(path.relative_to(ROOT)), sha256(path)))
    (SOURCE_DESTINATION / "PROMOTED.ok").write_text(
        "aspect_matched_narrow_s4u\n"
    )


if __name__ == "__main__":
    main()
