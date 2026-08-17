#!/usr/bin/env python3
"""Validate the indexed working-main Figure 2 A-K panel sequence."""

from __future__ import annotations

import csv
import hashlib
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
PANELS = ROOT / "figures/main/fig2_genetics/panels"
SIZE_SPEC = ROOT / "figures/layout_specs/figure2_panel_sizes.tsv"
LEGEND = ROOT / "docs/manuscript/05_figure_legends.md"
REINDEX_MANIFEST = PANELS / "Figure2_reindex_manifest.tsv"
PROHIBITED = ("ABF-fallback", "exploratory", "causal gene", "ancestry-specific")
TOLERANCE_IN = 0.015


def pdf_info(path: Path) -> tuple[int, float, float]:
    text = subprocess.check_output(["pdfinfo", str(path)], text=True)
    pages = int(re.search(r"^Pages:\s+(\d+)", text, re.M).group(1))
    size = re.search(r"^Page size:\s+([\d.]+) x ([\d.]+) pts", text, re.M)
    return pages, float(size.group(1)) / 72, float(size.group(2)) / 72


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    errors: list[str] = []
    with SIZE_SPEC.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = {
        Path(row["pdf"]).name: (float(row["width_in"]), float(row["height_in"]))
        for row in rows
    }
    indexed = sorted(expected)
    if len(indexed) != 11:
        errors.append(f"size index has {len(indexed)} panels, expected 11")

    for name, (want_w, want_h) in expected.items():
        path = PANELS / name
        if not path.exists():
            errors.append(f"missing indexed panel {name}")
            continue
        pages, width, height = pdf_info(path)
        if pages != 1:
            errors.append(f"{name}: {pages} pages")
        if abs(width - want_w) > TOLERANCE_IN or abs(height - want_h) > TOLERANCE_IN:
            errors.append(
                f"{name}: {width:.3f}x{height:.3f} in, expected {want_w:.2f}x{want_h:.2f}"
            )
        artwork = subprocess.check_output(["pdftotext", str(path), "-"], text=True)
        for term in PROHIBITED:
            if term.lower() in artwork.lower():
                errors.append(f"{name}: prohibited artwork label {term!r}")
        if b"/FontFile" not in path.read_bytes():
            errors.append(f"{name}: no embedded font object detected")
        if subprocess.run(
            ["gs", "-q", "-dNOPAUSE", "-dBATCH", "-sDEVICE=nullpage", str(path)],
            check=False,
        ).returncode:
            errors.append(f"{name}: Ghostscript parse failed")

    required_sources = (
        "Fig2B_PIP_vs_SuSiE-coloc_source.csv",
        "Fig2C_coding_noncoding_source.csv",
        "Fig2D_pip_architecture_by_trait_directness_source.tsv",
        "Fig2E_ancestry_unique_coloc_GWS_source.tsv",
        "Fig2F_crossancestry_coloc_source.tsv",
        "Fig2G_RORA_NAFLD_source.csv",
        "Fig2H_FABP1_ALT_source.csv",
        "Fig2J_eqtl_observability_source.csv",
        "Fig2K_phenotype_provenance_source.csv",
    )
    for name in required_sources:
        if not (PANELS / name).exists():
            errors.append(f"missing source sidecar {name}")

    fig2f_source = PANELS / "Fig2F_crossancestry_coloc_source.tsv"
    fig2f_pdf = PANELS / "Fig2F_crossancestry_coloc.pdf"
    if fig2f_source.exists() and fig2f_pdf.exists():
        with fig2f_source.open(newline="") as handle:
            fig2f_rows = list(csv.DictReader(handle, delimiter="\t"))
        observed_states = {row["evidence_state"] for row in fig2f_rows}
        fig2f_text = subprocess.check_output(["pdftotext", str(fig2f_pdf), "-"], text=True)
        state_labels = {
            "multi_signal": "Multi-signal COLOC",
            "single_signal_only": "Single-signal COLOC only",
            "evaluated_no_support": "Evaluated, PP.H4 ≤ 0.5",
            "not_evaluable": "Not evaluable",
        }
        for state, label in state_labels.items():
            if state not in observed_states and label.lower() in fig2f_text.lower():
                errors.append(f"Fig2F legend includes unused evidence state {state}")

    if not REINDEX_MANIFEST.exists():
        errors.append("missing Figure2_reindex_manifest.tsv")
    else:
        with REINDEX_MANIFEST.open(newline="") as handle:
            manifest_rows = list(csv.DictReader(handle, delimiter="\t"))
        if [row["panel"] for row in manifest_rows] != list("ABCDEFGHIJK"):
            errors.append("reindex manifest does not contain ordered panels A-K")
        for row in manifest_rows:
            path = PANELS / row["canonical_pdf"]
            if path.exists() and sha256(path) != row["sha256"]:
                errors.append(f"{row['canonical_pdf']}: hash disagrees with reindex manifest")

    stale = (
        "Fig2D_ancestry_unique_coloc_GWS.pdf",
        "Fig2E_crossancestry_coloc.pdf",
        "Fig2F_RORA_NAFLD.pdf",
        "Fig2G_FABP1_ALT.pdf",
        "Fig2H_locus_legend.pdf",
        "fig2I_eqtl_observability.pdf",
        "Fig2J_pip_architecture_by_trait_directness.pdf",
        "fig2K_phenotype_provenance.pdf",
    )
    for name in stale:
        if (PANELS / name).exists():
            errors.append(f"stale pre-reindex canonical file remains: {name}")

    legend_text = LEGEND.read_text()
    for term in PROHIBITED:
        if term.lower() in legend_text.lower():
            errors.append(f"legend contains prohibited label {term!r}")

    if errors:
        print("FAIL")
        print("\n".join(f"- {error}" for error in errors))
        return 1
    print("PASS")
    print(f"indexed_panels\t{len(indexed)}")
    print("one_page_exact_dimensions\tTRUE")
    print("embedded_fonts\tTRUE")
    print("prohibited_labels_absent\tTRUE")
    print("source_sidecars_present\tTRUE")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
