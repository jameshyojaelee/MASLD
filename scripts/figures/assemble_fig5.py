#!/usr/bin/env python3
"""Assemble the six validated Figure 5 panel PDFs on one vector PDF page."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("pdf")
import matplotlib.pyplot as plt
from PyPDF2 import PageObject, PdfReader, PdfWriter, Transformation


PAGE_WIDTH_IN = 7.2
PAGE_HEIGHT_IN = 9.2


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate-root", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = arguments()
    figure_root = args.candidate_root.resolve() / "figure5"
    panel_root = figure_root / "panels"
    output = figure_root / "fig5_molecular_context.pdf"
    if output.exists():
        raise RuntimeError(f"refusing to overwrite assembled Figure 5: {output}")

    panels = {
        "A": ("fig5a_input_firewall.pdf", 0.400, 0.200),
        "B": ("fig5b_protein_triage.pdf", 3.600, 0.200),
        "C": ("fig5c_mrna_protein_composite.pdf", 0.285, 2.450),
        "D": ("fig5d_snatac_accessibility.pdf", 0.228, 5.410),
        "E": ("fig5e_multimodal_program_summary.pdf", 3.348, 5.410),
        "F": ("fig5f_spatial_program_calibration.pdf", 1.209, 7.420),
    }
    missing = [name for name, (filename, _, _) in panels.items() if not (panel_root / filename).is_file()]
    if missing:
        raise RuntimeError(f"missing Figure 5 panels: {missing}")

    matplotlib.rcParams.update({
        "pdf.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Liberation Sans"],
        "font.size": 6,
    })
    label_layer = figure_root / ".fig5_labels.pdf"
    fig = plt.figure(figsize=(PAGE_WIDTH_IN, PAGE_HEIGHT_IN))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_axis_off()
    for label, (_, x_in, top_in) in panels.items():
        ax.text(
            x_in / PAGE_WIDTH_IN,
            1 - (top_in - 0.015) / PAGE_HEIGHT_IN,
            label,
            ha="left",
            va="bottom",
            fontsize=6,
            fontweight="normal",
            color="#202124",
            transform=ax.transAxes,
        )
    fig.savefig(label_layer, format="pdf", facecolor="white", edgecolor="none")
    plt.close(fig)

    base_reader = PdfReader(str(label_layer))
    page = base_reader.pages[0]
    for _, (filename, x_in, top_in) in panels.items():
        source = PdfReader(str(panel_root / filename)).pages[0]
        width_in = float(source.mediabox.width) / 72.0
        height_in = float(source.mediabox.height) / 72.0
        if x_in + width_in > PAGE_WIDTH_IN + 1e-6 or top_in + height_in > PAGE_HEIGHT_IN + 1e-6:
            raise RuntimeError(f"panel {filename} exceeds the assembly page")
        y_in = PAGE_HEIGHT_IN - top_in - height_in
        ctm = Transformation().translate(
            tx=x_in * 72.0, ty=y_in * 72.0
        ).ctm
        # PyPDF2 3.0.1's documented add_transformation()+merge_page() path
        # leaves the source trim box at the origin and clips translated panels.
        # This is the transformation-aware internal path used by the library's
        # retired mergeTransformedPage implementation; it transforms the clip
        # and content together while preserving vector panel content.
        page._merge_page(
            source,
            lambda content, source=source, ctm=ctm: PageObject._add_transformation_matrix(
                content, source.pdf, ctm
            ),
            ctm,
            False,
        )

    writer = PdfWriter()
    writer.add_page(page)
    with output.open("wb") as handle:
        writer.write(handle)
    label_layer.unlink()
    print(output)


if __name__ == "__main__":
    main()
