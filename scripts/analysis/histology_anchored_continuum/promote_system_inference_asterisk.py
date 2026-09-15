#!/usr/bin/env python3
"""Promote the reviewed strict-support asterisk revision of Figure S4U."""

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
    "histology-continuum-molecular-systems-significance-asterisk-20260825T015450Z"
)
FIG4 = ROOT / "figures/main/fig4_singlecell_programs"
SOURCE_PANEL = SOURCE / (
    "panels/s4_molecular_systems_pooled_continuum_windows_with_inference.pdf"
)
DESTINATION_PANEL = FIG4 / (
    "panels/supplementary/figs4u_molecular_systems_pooled_continuum_windows.pdf"
)
SOURCE_DESTINATION = FIG4 / (
    "source_tables/current_candidate/molecular_systems_inference_asterisk_20260825"
)
PROMOTION_MANIFEST = FIG4 / (
    "manifests/molecular_systems_inference_asterisk_promotion_20260825.tsv"
)
EXPECTED_OLD_SHA256 = (
    "7e80527d15b28208a731112797156a47bcd5df108925ed1ea27d6bb3d3cbc360"
)
EXPECTED_NEW_SHA256 = (
    "f75fab4cf487bd0c78184790058bf0e799dff2c58c12273583d569a7297bdc4a"
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


def staged_copy(source: Path, destination: Path) -> Path:
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
        raise SystemExit("Staged S4U panel hash mismatch")
    return temporary


def main() -> None:
    if not (SOURCE / "READY_FOR_REVIEW.ok").is_file():
        raise SystemExit("Strict-support S4U candidate is not review-ready")
    audit = read_tsv(SOURCE / "provenance/audit.tsv")
    if not audit or any(row["passed"].upper() != "TRUE" for row in audit):
        raise SystemExit("Strict-support S4U audit did not pass")
    figures = read_tsv(SOURCE / "figure_manifest.tsv")
    if len(figures) != 1 or figures[0]["promoted"].upper() != "FALSE":
        raise SystemExit("Strict-support S4U figure manifest drift")
    if Path(figures[0]["path"]).resolve() != SOURCE_PANEL.resolve():
        raise SystemExit("Strict-support S4U figure path drift")
    if sha256(SOURCE_PANEL) != EXPECTED_NEW_SHA256:
        raise SystemExit("Reviewed strict-support S4U hash drift")
    for row in read_tsv(SOURCE / "provenance/output_checksums.tsv"):
        path = Path(row["path"])
        if not path.is_file() or sha256(path) != row["sha256"]:
            raise SystemExit(f"Strict-support candidate output drift: {path}")
    strip = read_tsv(SOURCE / "source_tables/system_inference_strip.tsv")
    if len(strip) != 86:
        raise SystemExit("Strict-support source family drift")
    rendered = [row for row in strip if row["rendered_in_panel"].upper() == "TRUE"]
    if len(rendered) != 18 or any(
        row["inference_column"] != "maxT" or row["render_symbol"] != "*"
        for row in rendered
    ):
        raise SystemExit("Strict-support marker registry drift")
    if not DESTINATION_PANEL.is_file() or sha256(DESTINATION_PANEL) != EXPECTED_OLD_SHA256:
        raise SystemExit("Active S4U hash drift")
    if SOURCE_DESTINATION.exists() or PROMOTION_MANIFEST.exists():
        raise SystemExit("Strict-support S4U promotion namespace already exists")

    temporary = staged_copy(SOURCE_PANEL, DESTINATION_PANEL)
    try:
        os.replace(temporary, DESTINATION_PANEL)
    finally:
        temporary.unlink(missing_ok=True)
    if sha256(DESTINATION_PANEL) != EXPECTED_NEW_SHA256:
        raise SystemExit("Promoted strict-support S4U hash mismatch")

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

    with PROMOTION_MANIFEST.open("x", newline="") as handle:
        fields = (
            "callout", "destination", "source", "old_sha256", "new_sha256",
            "complete_system_family", "primary_bh_supported",
            "strict_maxT_supported", "rendered_strict_markers",
            "window_values_changed", "display_order_changed",
            "windows_inferentially_tested", "revision",
        )
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({
            "callout": "S4U",
            "destination": str(DESTINATION_PANEL.relative_to(ROOT)),
            "source": str(SOURCE_PANEL.relative_to(ROOT)),
            "old_sha256": EXPECTED_OLD_SHA256,
            "new_sha256": EXPECTED_NEW_SHA256,
            "complete_system_family": "43/43",
            "primary_bh_supported": "33/43",
            "strict_maxT_supported": "18/43",
            "rendered_strict_markers": "18/43",
            "window_values_changed": "FALSE",
            "display_order_changed": "FALSE",
            "windows_inferentially_tested": "FALSE",
            "revision": "strict_maxt_asterisk_no_support_columns",
        })

    checksum_path = SOURCE_DESTINATION / "promoted_source_checksums.tsv"
    with checksum_path.open("x", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("path", "sha256"))
        for path in sorted(copied):
            writer.writerow((str(path.relative_to(ROOT)), sha256(path)))
    (SOURCE_DESTINATION / "PROMOTED.ok").write_text(
        "user_approved_s4u_strict_support_asterisk_revision\n"
    )


if __name__ == "__main__":
    main()
