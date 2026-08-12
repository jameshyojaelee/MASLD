#!/usr/bin/env python3
"""Validate the six-panel candidate Figure 4 and exact source-table hashes."""

from __future__ import annotations

import argparse
import subprocess
from pathlib import Path

from spatial_resource_lib import RESOURCE_RELEASE_ID, SpatialResourceError, sha256_file, write_tsv


PANELS = {
    "4A": ("fig4a_overview_cascade.pdf", "resource_spatial_firewall"),
    "4B": ("fig4b_protein_triage.pdf", "accepted_protein_triage_regenerated"),
    "4C": ("fig4c_mrna_protein_composite.pdf", "accepted_protein_cornerstone_regenerated"),
    "4D": ("fig4d_snatac_accessibility.pdf", "accepted_donor_aware_chromatin_context_regenerated"),
    "4E": ("fig4e_spatial_program_maps.pdf", "two_predeclared_program_maps"),
    "4F": ("fig4f_multimodal_program_summary.pdf", "two_program_source_state_matrix"),
}

SOURCE_TABLES = {
    "4A": "main/fig4_validation/data/fig4a_spatial_firewall.tsv",
    "4B": "main/fig4_validation/panels/data/fig4b_protein_triage_summary.tsv",
    "4C": "main/fig4_validation/data/composite_mrna_protein_corrected.csv",
    "4D": "main/fig4_validation/data/fig4d_snatac_accessibility_source.tsv",
    "4E": "main/fig4_validation/data/fig4e_spatial_program_maps_selection.tsv",
    "4F": "main/fig4_validation/data/fig4f_two_program_source_matrix.tsv",
}


def pdf_pages(path: Path) -> int:
    result = subprocess.run(
        ["pdfinfo", str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    for line in result.stdout.splitlines():
        if line.startswith("Pages:"):
            return int(line.split(":", 1)[1].strip())
    raise SpatialResourceError(f"pdfinfo did not report a page count: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-root", type=Path, required=True)
    args = parser.parse_args()
    candidate = args.candidate_root.resolve()
    figure_root = candidate / "figures"
    panels = figure_root / "main/fig4_validation/panels"
    if (figure_root / "READY").exists():
        raise SpatialResourceError("immutable candidate Figure 4 is already validated")

    panel_rows = []
    source_rows = []
    for callout, (filename, role) in PANELS.items():
        path = panels / filename
        if not path.is_file() or path.stat().st_size == 0:
            raise SpatialResourceError(f"missing candidate panel {callout}: {path}")
        pages = pdf_pages(path)
        if pages != 1:
            raise SpatialResourceError(f"candidate panel {callout} has {pages} PDF pages")
        panel_rows.append({
            "release_id": RESOURCE_RELEASE_ID,
            "callout": callout,
            "relative_path": path.relative_to(candidate).as_posix(),
            "role": role,
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "pdf_pages": pages,
            "status": "candidate_not_promoted",
        })
        source = figure_root / SOURCE_TABLES[callout]
        if not source.is_file() or source.stat().st_size == 0:
            raise SpatialResourceError(f"missing source table for {callout}: {source}")
        source_rows.append({
            "release_id": RESOURCE_RELEASE_ID,
            "callout": callout,
            "relative_path": source.relative_to(candidate).as_posix(),
            "bytes": source.stat().st_size,
            "sha256": sha256_file(source),
        })
    write_tsv(
        figure_root / "figure4_panel_manifest.tsv",
        ("release_id", "callout", "relative_path", "role", "bytes", "sha256", "pdf_pages", "status"),
        panel_rows,
    )
    write_tsv(
        figure_root / "figure4_source_table_manifest.tsv",
        ("release_id", "callout", "relative_path", "bytes", "sha256"),
        source_rows,
    )
    write_tsv(
        figure_root / "READY",
        ("release_id", "status", "n_panels", "n_source_tables", "panel_manifest_sha256", "source_manifest_sha256", "canonical_figure_written"),
        [{
            "release_id": RESOURCE_RELEASE_ID,
            "status": "six_panel_candidate_figure4_validated_awaiting_adjudication",
            "n_panels": len(panel_rows),
            "n_source_tables": len(source_rows),
            "panel_manifest_sha256": sha256_file(figure_root / "figure4_panel_manifest.tsv"),
            "source_manifest_sha256": sha256_file(figure_root / "figure4_source_table_manifest.tsv"),
            "canonical_figure_written": "FALSE",
        }],
    )
    print("PASS six candidate Figure 4 PDFs and source tables")


if __name__ == "__main__":
    main()
