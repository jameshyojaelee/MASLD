#!/usr/bin/env python3
"""Render candidate-only Figure 4A, 4E, and 4F from sealed source tables."""

from __future__ import annotations

import argparse
import gzip
from pathlib import Path

import matplotlib

matplotlib.use("pdf")
matplotlib.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 6,
    "axes.titlesize": 6,
    "axes.labelsize": 6,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
})
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import FancyBboxPatch, Rectangle
import pandas as pd

from spatial_resource_lib import RESOURCE_RELEASE_ID, SpatialResourceError, sha256_file, write_tsv


CYAN = "#007C91"
CYAN_MID = "#68B7C7"
CYAN_LIGHT = "#D9EEF2"
GRAY = "#9E9E9E"
INK = "#222222"
WHITE = "#FFFFFF"


def panel_a(output: Path) -> Path:
    rows = [
        {"stage": "frozen_input", "n_programs": 117, "role": "immutable disease-state programs"},
        {"stage": "coverage", "n_programs": 117, "role": "complete outcome-free observability by assay"},
        {"stage": "confirmatory_family", "n_programs": 2, "role": "predeclared IGFBP7 and BICC1 programs"},
        {"stage": "assay_native_inference", "n_programs": 2, "role": "matched-null inference only where unit contracts permit"},
    ]
    write_tsv(output.parent / "data/fig4a_spatial_firewall.tsv", ("stage", "n_programs", "role"), rows)
    fig, ax = plt.subplots(figsize=(5.5, 1.55))
    ax.set_xlim(0, 5.5)
    ax.set_ylim(0, 1.55)
    ax.axis("off")
    boxes = [
        (0.10, 0.48, 1.08, 0.62, "117 frozen\nprograms", CYAN_LIGHT),
        (1.48, 0.48, 1.20, 0.62, "117-program\ncoverage", CYAN_LIGHT),
        (3.00, 0.48, 1.00, 0.62, "2 predeclared\nprograms", CYAN_MID),
        (4.30, 0.48, 1.10, 0.62, "assay-native\ninference", CYAN),
    ]
    for x, y, width, height, label, color in boxes:
        patch = FancyBboxPatch((x, y), width, height, boxstyle="round,pad=0.02,rounding_size=0.04", facecolor=color, edgecolor=INK, linewidth=0.35)
        ax.add_patch(patch)
        ax.text(x + width / 2, y + height / 2, label, ha="center", va="center", color=WHITE if color == CYAN else INK)
    for start, stop in ((1.18, 1.48), (2.68, 3.00), (4.00, 4.30)):
        ax.annotate("", xy=(stop, 0.79), xytext=(start, 0.79), arrowprops={"arrowstyle": "-|>", "lw": 0.5, "color": GRAY})
    ax.text(2.74, 1.30, "No spatial outcome used for coverage or selection", ha="center", color=INK)
    ax.text(2.74, 0.20, "Unresolved units or joins: source-dependent / metadata-pending / untestable", ha="center", color="#555555")
    path = output / "fig4a_overview_cascade.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def panel_e(project: Path, output: Path) -> Path:
    source = (
        project / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07"
        / "spatial_context/map_source_v2/per_spot_program_map.tsv.gz"
    )
    selection = source.parent / "map_selection.tsv"
    with gzip.open(source, "rt", encoding="utf-8") as handle:
        maps = pd.read_csv(handle, sep="\t")
    if set(maps["program_uid"]) != {
        "hotspot_hepatocytes_f05c535ae5bbc0b9", "hotspot_hepatocytes_48f39dd4d817a10e"
    }:
        raise SpatialResourceError("Figure 4E source is not the exact two-program family")
    selected = pd.read_csv(selection, sep="\t")
    selected.to_csv(output.parent / "data/fig4e_spatial_program_maps_selection.tsv", sep="\t", index=False)
    program_order = [
        ("hotspot_hepatocytes_f05c535ae5bbc0b9", "Stromal ECM (IGFBP7)"),
        ("hotspot_hepatocytes_48f39dd4d817a10e", "Ductular injury (BICC1)"),
    ]
    dataset_order = [("GSE192741", "GSE192741 · 4 donors / 5 sections"), ("Vu_et_al_2025", "Vu et al. · 10 arrays; donor key unresolved")]
    cmap = LinearSegmentedColormap.from_list("spatial_cyan", ["#E7F3F6", CYAN_MID, CYAN])
    norm = Normalize(vmin=-2, vmax=2, clip=True)
    fig, axes = plt.subplots(2, 2, figsize=(4.30, 3.15))
    for row_index, (dataset, dataset_label) in enumerate(dataset_order):
        for column_index, (uid, program_label) in enumerate(program_order):
            ax = axes[row_index, column_index]
            part = maps[(maps["dataset"] == dataset) & (maps["program_uid"] == uid)]
            ax.scatter(part["x"], -part["y"], c=part["residual_program_score_z"], cmap=cmap, norm=norm, s=2.0, linewidths=0, rasterized=True)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if row_index == 0:
                ax.set_title(program_label, pad=2)
            if column_index == 0:
                ax.text(-0.04, 0.5, dataset_label, transform=ax.transAxes, ha="right", va="center", rotation=90)
            unit = str(part["reporting_unit_id"].iloc[0])
            ax.text(0.01, 0.01, unit, transform=ax.transAxes, ha="left", va="bottom", color="#555555")
    fig.subplots_adjust(left=0.12, right=0.99, top=0.92, bottom=0.14, wspace=0.04, hspace=0.08)
    cax = fig.add_axes([0.34, 0.055, 0.34, 0.018])
    colorbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    colorbar.set_ticks([-2, 0, 2])
    colorbar.ax.tick_params(width=0.35, length=1.5, pad=1)
    colorbar.outline.set_linewidth(0.35)
    colorbar.set_label("abundance-adjusted program score (z)", labelpad=1)
    path = output / "fig4e_spatial_program_maps.pdf"
    fig.savefig(path, bbox_inches="tight", dpi=400)
    plt.close(fig)
    return path


def panel_f(candidate: Path, output: Path) -> Path:
    effects = pd.read_parquet(candidate / "effects/spatial_program_effects.parquet")
    dataset_order = [
        "GSE192741", "Vu_et_al_2025", "HRA007511_HMSMA", "Yakubovsky2026",
        "Govaere2026_CosMx", "Govaere2026_GeoMx", "GSE287826",
        "GSE244832", "GSE281367", "PXD051911",
    ]
    program_order = [
        "hotspot_hepatocytes_f05c535ae5bbc0b9",
        "hotspot_hepatocytes_48f39dd4d817a10e",
    ]
    program_labels = {
        program_order[0]: "Stromal ECM\n(IGFBP7)",
        program_order[1]: "Ductular injury\n(BICC1)",
    }
    states = ["supported", "source_dependent", "indeterminate", "untestable", "not_applicable", "skipped", "dropped"]
    colors = {
        "supported": CYAN,
        "source_dependent": CYAN_LIGHT,
        "indeterminate": "#D8D8D8",
        "untestable": "#9E9E9E",
        "not_applicable": WHITE,
        "skipped": "#F1F1F1",
        "dropped": "#4D4D4D",
    }
    symbols = {"supported": "●", "source_dependent": "S", "indeterminate": "?", "untestable": "×", "not_applicable": "–", "skipped": "skip", "dropped": "drop"}
    source_table = effects[effects["dataset_id"].isin(dataset_order)].copy()
    source_table.to_csv(output.parent / "data/fig4f_two_program_source_matrix.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(4.85, 3.55))
    ax.set_xlim(0, 2)
    ax.set_ylim(0, len(dataset_order))
    ax.invert_yaxis()
    for row_index, dataset in enumerate(dataset_order):
        for column_index, uid in enumerate(program_order):
            row = source_table[(source_table["dataset_id"] == dataset) & (source_table["program_uid"] == uid)]
            if len(row) != 1:
                raise SpatialResourceError(f"Figure 4F requires one row for {dataset}/{uid}, found {len(row)}")
            state = str(row.iloc[0]["evidence_state"])
            ax.add_patch(Rectangle((column_index, row_index), 1, 1, facecolor=colors[state], edgecolor=WHITE, linewidth=1.0))
            text_color = WHITE if state == "dropped" else INK
            ax.text(column_index + 0.5, row_index + 0.52, symbols[state], ha="center", va="center", color=text_color)
    ax.set_xticks([0.5, 1.5], [program_labels[uid] for uid in program_order])
    ax.xaxis.tick_top()
    ax.set_yticks([index + 0.5 for index in range(len(dataset_order))], [label.replace("_", " ") for label in dataset_order])
    ax.tick_params(length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    handles = [Rectangle((0, 0), 1, 1, facecolor=colors[state], edgecolor=GRAY, linewidth=0.35) for state in states]
    ax.legend(handles, [state.replace("_", " ") for state in states], ncol=2, frameon=False, loc="upper left", bbox_to_anchor=(1.02, 1.01), handlelength=1.0, handleheight=0.8, columnspacing=0.8, borderaxespad=0)
    ax.text(0, len(dataset_order) + 0.40, "S = within-source support, but source-dependent unit or join", ha="left", va="top", color="#555555")
    path = output / "fig4f_multimodal_program_summary.pdf"
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--candidate-root", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    candidate = args.candidate_root.resolve()
    output = candidate / "figures/main/fig4_validation/panels"
    if (candidate / "figures/READY").exists():
        raise SpatialResourceError("immutable candidate figures already exist")
    output.mkdir(parents=True, exist_ok=True)
    (output.parent / "data").mkdir(parents=True, exist_ok=True)
    paths = [panel_a(output), panel_e(project, output), panel_f(candidate, output)]
    write_tsv(
        candidate / "figures/python_panel_manifest.tsv",
        ("release_id", "panel", "relative_path", "bytes", "sha256", "status"),
        [{
            "release_id": RESOURCE_RELEASE_ID,
            "panel": path.stem.split("_")[0].replace("fig", "").upper(),
            "relative_path": path.relative_to(candidate).as_posix(),
            "bytes": path.stat().st_size,
            "sha256": sha256_file(path),
            "status": "candidate_not_promoted",
        } for path in paths],
    )
    print(f"WROTE candidate panels: {', '.join(path.name for path in paths)}")


if __name__ == "__main__":
    main()
