#!/usr/bin/env python3
# KEY MESSAGE: spatial figures should show one biological decision per panel;
# provenance and claim boundaries remain in source tables rather than the artwork.
"""Render a compact, publication-facing spatial figure family from sealed results."""

from __future__ import annotations

import argparse
import os
import shutil
import sys
from pathlib import Path

import matplotlib

matplotlib.use("pdf")
matplotlib.rcParams.update(
    {
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 6,
        "axes.titlesize": 6,
        "axes.labelsize": 6,
        "xtick.labelsize": 6,
        "ytick.labelsize": 6,
        "legend.fontsize": 6,
        "figure.titlesize": 6,
        "axes.linewidth": 0.45,
        "xtick.major.width": 0.45,
        "ytick.major.width": 0.45,
    }
)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch, Patch, Rectangle
import numpy as np
import pandas as pd

from spatial_resource_lib import SpatialResourceError, sha256_file


RELEASE_ID = "spatial-publication-figures-candidate-2026-08-11-r4"
CYAN = "#007C91"
CYAN_MID = "#68B7C7"
CYAN_LIGHT = "#CDE7EE"
GRAY = "#9E9E9E"
GRAY_LIGHT = "#E6E6E5"
DARK_GRAY = "#4D4D4D"
INK = "#222222"
WHITE = "#FFFFFF"
MAGENTA = "#C9265E"

ROOTS = {
    "resource": "Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11",
    "downstream": "Analysis/Multimodal_Program_Projection/candidates/spatial-figure-integration-candidate-2026-08-11-r2",
    "composition": "Analysis/Multimodal_Program_Projection/candidates/spatial-full-composition-figure-candidate-2026-08-11-r4",
    "cell_context": "Analysis/Multimodal_Program_Projection/candidates/spatial-cell-context-attribution-candidate-2026-08-11-r2",
    "unit_null": "Analysis/Multimodal_Program_Projection/candidates/spatial-unit-null-diagnostics-candidate-2026-08-11-r2",
    "mask": "Analysis/Multimodal_Program_Projection/candidates/hmsma-mask-sensitivity-candidate-2026-08-11-r2",
}
READY_RELATIVE = {"resource": "validation/READY"}

PANEL_NAMES = (
    "fig1_spatial_assay_observability.pdf",
    "fig4a_spatial_firewall.pdf",
    "fig4e_spatial_program_maps.pdf",
    "fig4f_multimodal_program_summary.pdf",
    "fig5_spatial_passport_examples.pdf",
    "figS_spatial_program_coverage.pdf",
    "figS_hmsma_label_blind_organization.pdf",
    "figS_full_composition_sensitivity.pdf",
    "figS_cell_context_association.pdf",
    "figS_spatial_unit_heterogeneity.pdf",
    "figS_spatial_matched_null_calibration.pdf",
    "figS_spatial_weight_sensitivity.pdf",
    "figS_hmsma_proxy_mask_sensitivity.pdf",
)


def save_pdf(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, format="pdf", facecolor=WHITE, bbox_inches="tight", dpi=400)
    plt.close(fig)


def clean_axis(ax: plt.Axes, left: bool = True, bottom: bool = True) -> None:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(left)
    ax.spines["bottom"].set_visible(bottom)
    ax.tick_params(length=2, pad=2)


def state_from_join(value: object) -> str:
    value = str(value)
    if value == "resolved":
        return "resolved"
    if value == "not_applicable":
        return "not_applicable"
    return "unresolved"


def unit_abbreviation(row: pd.Series) -> str:
    unit_codes = {"section": "S", "physical_array": "A", "AOI": "AOI"}
    technical_count = str(row["n_technical"]).strip()
    tech = "—" if technical_count in {"", "nan", "None"} else str(int(float(technical_count)))
    tech = f"{tech}{unit_codes.get(str(row['technical_unit']), '')}"
    if str(row["biological_unit_resolution"]) != "resolved" or pd.isna(row["n_biological"]):
        bio = "?D"
    else:
        bio = f"{int(float(row['n_biological']))}D"
    return f"{bio} / {tech}"


def render_fig1(source: pd.DataFrame, path: Path) -> pd.DataFrame:
    # The main inset keeps only decision-bearing fields. Coordinates are resolved
    # throughout and therefore remain in the complete source table, not the art.
    display = source.copy()
    display["unit_state"] = display["biological_unit_resolution"].map(state_from_join)
    display["clinical_state"] = np.where(
        (display["donor_join"] == "resolved") & (display["phenotype_join"].isin(["resolved", "not_applicable"])),
        "resolved",
        "unresolved",
    )
    display["histology_state"] = display["histology_join"].map(state_from_join)
    observable = pd.to_numeric(display["n_observable"], errors="raise")
    display["coverage_state"] = np.select(
        [observable >= 100, observable > 0],
        ["resolved", "partial"],
        default="unresolved",
    )
    display["unit_label"] = display.apply(unit_abbreviation, axis=1)
    display["gate_mark"] = display["dataset_gate"].map(
        {"pass": "", "source_dependent": "*", "metadata_pending": "†", "dropped": "×"}
    ).fillna("")

    columns = ["unit_state", "clinical_state", "histology_state", "coverage_state"]
    labels = ["unit", "clinical", "H&E", "programs"]
    colors = {
        "resolved": CYAN,
        "partial": CYAN_LIGHT,
        "unresolved": GRAY,
        "not_applicable": WHITE,
    }
    fig, ax = plt.subplots(figsize=(4.55, 2.52))
    ax.set_xlim(0, 5.45)
    ax.set_ylim(0, len(display))
    ax.invert_yaxis()
    for i, row in display.iterrows():
        for j, column in enumerate(columns):
            state = row[column]
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=colors[state], edgecolor=WHITE, linewidth=0.9))
            if column == "coverage_state":
                text = f"{int(row['n_observable'])}/117"
                color = WHITE if state == "resolved" else INK
                ax.text(j + 0.5, i + 0.52, text, ha="center", va="center", color=color)
            elif state == "resolved":
                ax.text(j + 0.5, i + 0.52, "•", ha="center", va="center", color=WHITE)
            elif state == "not_applicable":
                ax.text(j + 0.5, i + 0.52, "–", ha="center", va="center", color=INK)
        ax.text(4.18, i + 0.52, row["unit_label"], ha="left", va="center", color=INK)
    ylabels = [f"{row.display_label}{row.gate_mark}" for row in display.itertuples(index=False)]
    ax.set_xticks([0.5, 1.5, 2.5, 3.5, 4.52], labels + ["bio / tech"])
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(display)) + 0.5, ylabels)
    ax.tick_params(length=0, pad=2)
    for spine in ax.spines.values():
        spine.set_visible(False)
    handles = [
        Patch(facecolor=CYAN, edgecolor=GRAY, linewidth=0.3, label="resolved"),
        Patch(facecolor=CYAN_LIGHT, edgecolor=GRAY, linewidth=0.3, label="partial"),
        Patch(facecolor=GRAY, edgecolor=GRAY, linewidth=0.3, label="unresolved"),
        Patch(facecolor=WHITE, edgecolor=GRAY, linewidth=0.3, label="n/a"),
    ]
    ax.legend(handles=handles, frameon=False, ncol=4, loc="upper left", bbox_to_anchor=(0, -0.035), borderaxespad=0, handlelength=0.9, columnspacing=0.8)
    ax.text(0, -0.18, "* source-dependent   † metadata pending", transform=ax.transAxes, ha="left", va="top", color=DARK_GRAY)
    fig.subplots_adjust(left=0.20, right=0.99, top=0.88, bottom=0.30)
    save_pdf(fig, path)
    return display


def render_fig4a(path: Path) -> pd.DataFrame:
    source = pd.DataFrame(
        [
            {"path": "coverage", "n_programs": 117, "action": "report observability in every assay"},
            {"path": "inference", "n_programs": 2, "action": "test only where biological units resolve"},
        ]
    )
    fig, ax = plt.subplots(figsize=(4.15, 1.42))
    ax.set_xlim(0, 4.15)
    ax.set_ylim(0, 1.42)
    ax.axis("off")
    left = FancyBboxPatch((0.05, 0.43), 0.95, 0.55, boxstyle="round,pad=0.02,rounding_size=0.035", facecolor=CYAN_LIGHT, edgecolor=GRAY, linewidth=0.45)
    ax.add_patch(left)
    ax.text(0.525, 0.705, "117 frozen\nprograms", ha="center", va="center")
    ax.annotate("", xy=(1.34, 0.71), xytext=(1.02, 0.71), arrowprops={"arrowstyle": "-|>", "lw": 0.5, "color": GRAY})
    coverage = FancyBboxPatch((1.36, 0.78), 1.25, 0.48, boxstyle="round,pad=0.02,rounding_size=0.035", facecolor=CYAN_LIGHT, edgecolor=GRAY, linewidth=0.45)
    inference = FancyBboxPatch((1.36, 0.16), 1.25, 0.48, boxstyle="round,pad=0.02,rounding_size=0.035", facecolor=CYAN_MID, edgecolor=GRAY, linewidth=0.45)
    ax.add_patch(coverage)
    ax.add_patch(inference)
    ax.text(1.985, 1.02, "coverage: all 117", ha="center", va="center")
    ax.text(1.985, 0.40, "inference: 2", ha="center", va="center")
    ax.annotate("", xy=(2.95, 1.02), xytext=(2.63, 1.02), arrowprops={"arrowstyle": "-|>", "lw": 0.5, "color": GRAY})
    ax.annotate("", xy=(2.95, 0.40), xytext=(2.63, 0.40), arrowprops={"arrowstyle": "-|>", "lw": 0.5, "color": GRAY})
    ax.text(3.02, 1.02, "observability\nin every assay", ha="left", va="center")
    ax.text(3.02, 0.40, "only where units\npermit inference", ha="left", va="center")
    ax.plot([1.15, 1.15], [0.18, 1.24], color=DARK_GRAY, linewidth=0.55)
    ax.text(1.15, 1.34, "firewall", ha="center", va="bottom", color=DARK_GRAY)
    save_pdf(fig, path)
    return source


def render_fig4f(effects: pd.DataFrame, path: Path) -> pd.DataFrame:
    dataset_order = [
        "GSE192741", "Vu_et_al_2025", "HRA007511_HMSMA", "Yakubovsky2026",
        "Govaere2026_CosMx", "Govaere2026_GeoMx", "GSE287826",
        "GSE244832", "GSE281367", "PXD051911",
    ]
    labels = {
        "GSE192741": "GSE192741", "Vu_et_al_2025": "Vu et al.",
        "HRA007511_HMSMA": "HMSMA", "Yakubovsky2026": "Yakubovsky",
        "Govaere2026_CosMx": "CosMx", "Govaere2026_GeoMx": "GeoMx",
        "GSE287826": "GSE287826", "GSE244832": "GSE244832",
        "GSE281367": "GSE281367", "PXD051911": "PXD051911",
    }
    programs = [
        "hotspot_hepatocytes_f05c535ae5bbc0b9",
        "hotspot_hepatocytes_48f39dd4d817a10e",
    ]
    program_labels = ["Stromal ECM\n(IGFBP7)", "Ductular injury\n(BICC1)"]
    colors = {
        "supported": CYAN,
        "source_dependent": CYAN_LIGHT,
        "indeterminate": GRAY_LIGHT,
        "untestable": GRAY,
        "not_applicable": WHITE,
        "skipped": "#F5F5F5",
        "dropped": DARK_GRAY,
    }
    state_text = {
        "supported": "supported",
        "source_dependent": "source dep.",
        "indeterminate": "indeterminate",
        "untestable": "untestable",
        "not_applicable": "n/a",
        "skipped": "skipped",
        "dropped": "dropped",
    }
    source = effects[effects["dataset_id"].isin(dataset_order) & effects["program_uid"].isin(programs)].copy()
    if len(source) != 20:
        raise SpatialResourceError(f"Figure 4F requires 20 rows, found {len(source)}")
    fig, ax = plt.subplots(figsize=(3.95, 3.08))
    ax.set_xlim(0, 2)
    ax.set_ylim(0, len(dataset_order))
    ax.invert_yaxis()
    for i, dataset in enumerate(dataset_order):
        for j, uid in enumerate(programs):
            row = source[(source["dataset_id"] == dataset) & (source["program_uid"] == uid)]
            if len(row) != 1:
                raise SpatialResourceError(f"missing Figure 4F row: {dataset}/{uid}")
            state = str(row.iloc[0]["evidence_state"])
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=colors[state], edgecolor=WHITE, linewidth=0.9))
            color = WHITE if state in {"supported", "dropped"} else INK
            ax.text(j + 0.5, i + 0.52, state_text[state], ha="center", va="center", color=color)
    ax.set_xticks([0.5, 1.5], program_labels)
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(dataset_order)) + 0.5, [labels[x] for x in dataset_order])
    ax.tick_params(length=0, pad=2)
    for y in (7, 9):
        ax.axhline(y, color=INK, linewidth=0.5)
    ax.text(-0.58, 3.5, "spatial RNA", rotation=90, ha="center", va="center", color=DARK_GRAY, clip_on=False)
    ax.text(-0.58, 8.0, "chromatin", rotation=90, ha="center", va="center", color=DARK_GRAY, clip_on=False)
    ax.text(-0.58, 9.5, "protein", rotation=90, ha="center", va="center", color=DARK_GRAY, clip_on=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.subplots_adjust(left=0.28, right=0.99, top=0.87, bottom=0.03)
    save_pdf(fig, path)
    return source


def render_fig5(source: pd.DataFrame, path: Path) -> pd.DataFrame:
    display = source.copy()
    messages = {
        "igfbp7_program": ("reproducible organization\n(composition-linked)", "donor-resolved\nRNA + protein"),
        "bicc1_program": ("adequate test; no\norganization call", "targeted assay in\nductular regions"),
        "hmsma_clinical_gate": ("35 arrays; no donor\nor clinical key", "authoritative join\n+ registered H&E"),
    }
    display["spatial_message"] = display["example_id"].map(lambda x: messages[x][0])
    display["short_next_test"] = display["example_id"].map(lambda x: messages[x][1])
    state_colors = {"supported": CYAN, "indeterminate": GRAY_LIGHT, "untestable": GRAY}
    fig, ax = plt.subplots(figsize=(4.95, 1.88))
    ax.set_xlim(0, 6.35)
    ax.set_ylim(0, 3.45)
    ax.invert_yaxis()
    ax.axis("off")
    headers = [(0.05, "example"), (1.40, "state"), (2.55, "what spatial adds"), (4.60, "next test")]
    for x, label in headers:
        ax.text(x, 0.22, label, ha="left", va="center")
    for i, row in display.iterrows():
        y = 0.52 + i * 0.93
        ax.axhline(y + 0.82, color=GRAY_LIGHT, linewidth=0.45, xmin=0.0, xmax=1.0)
        ax.text(0.05, y + 0.41, str(row["display_label"]), ha="left", va="center")
        state = str(row["primary_evidence_state"])
        box = FancyBboxPatch((1.40, y + 0.18), 0.92, 0.44, boxstyle="round,pad=0.02,rounding_size=0.04", facecolor=state_colors[state], edgecolor=GRAY, linewidth=0.35)
        ax.add_patch(box)
        ax.text(1.86, y + 0.40, state, ha="center", va="center", color=WHITE if state in {"supported", "untestable"} else INK)
        ax.text(2.55, y + 0.41, row["spatial_message"], ha="left", va="center")
        ax.annotate("", xy=(4.47, y + 0.41), xytext=(4.15, y + 0.41), arrowprops={"arrowstyle": "-|>", "lw": 0.45, "color": GRAY})
        ax.text(4.60, y + 0.41, row["short_next_test"], ha="left", va="center")
    save_pdf(fig, path)
    return display


def render_hmsma_label_blind(source: pd.DataFrame, path: Path) -> pd.DataFrame:
    pivot = source.pivot(index="array_id", columns="program_label", values="residual_moran_i")
    ig = pivot["Stromal ECM (IGFBP7)"].to_numpy(float)
    bi = pivot["Ductular injury (BICC1)"].to_numpy(float)
    x = np.array([0, 1])
    fig, ax = plt.subplots(figsize=(2.75, 2.45))
    for left, right in zip(ig, bi, strict=True):
        ax.plot(x, [left, right], color=GRAY_LIGHT, linewidth=0.5, zorder=1)
    ax.scatter(np.zeros(len(ig)), ig, s=9, color=CYAN, linewidth=0, zorder=2)
    ax.scatter(np.ones(len(bi)), bi, s=9, color=CYAN_MID, linewidth=0, zorder=2)
    medians = [float(np.median(ig)), float(np.median(bi))]
    ax.hlines(medians, x - 0.14, x + 0.14, color=INK, linewidth=1.0, zorder=3)
    for xx, median in zip(x, medians, strict=True):
        ax.text(xx, median + 0.045, f"median {median:.3f}", ha="center", va="bottom")
    ax.axhline(0, color=GRAY, linewidth=0.45, linestyle="--")
    ax.set_xticks(x, ["Stromal ECM\n(IGFBP7)", "Ductular injury\n(BICC1)"])
    ax.set_ylabel("Residual Moran's I")
    ax.set_xlim(-0.35, 1.35)
    clean_axis(ax, left=True, bottom=True)
    fig.subplots_adjust(left=0.21, right=0.98, top=0.95, bottom=0.22)
    save_pdf(fig, path)
    return source


def render_composition(source: pd.DataFrame, path: Path) -> pd.DataFrame:
    source = source.sort_values("display_order").copy()
    y = np.arange(len(source))[::-1]
    before = source["primary_centered_moran_i"].to_numpy(float)
    after = source["full_composition_centered_moran_i"].to_numpy(float)
    fig, ax = plt.subplots(figsize=(4.25, 2.02))
    ax.axvline(0, color=GRAY_LIGHT, linewidth=0.5)
    for yy, row, x0, x1 in zip(y, source.itertuples(index=False), before, after, strict=True):
        open_symbol = row.dataset == "Vu_et_al_2025"
        ax.annotate("", xy=(x1, yy), xytext=(x0, yy), arrowprops={"arrowstyle": "-|>", "lw": 0.65, "color": GRAY})
        ax.scatter(x0, yy, s=18, marker="o", facecolor=WHITE if open_symbol else GRAY, edgecolor=GRAY, linewidth=0.65, zorder=3)
        ax.scatter(x1, yy, s=20, marker="D", facecolor=WHITE if open_symbol else CYAN, edgecolor=CYAN, linewidth=0.65, zorder=4)
        ax.text(x1, yy - 0.24, f"{x1:.3f}", ha="center", va="top", color=CYAN)
    short_labels = []
    for row in source.itertuples(index=False):
        program = "IGFBP7" if "IGFBP7" in row.program_label else "BICC1"
        dataset = "GSE" if row.dataset == "GSE192741" else "Vu*"
        short_labels.append(f"{program}  |  {dataset}")
    ax.set_yticks(y, short_labels)
    ax.set_xlabel("Centered Moran's I relative to matched-gene null")
    ax.set_xlim(-0.005, max(before.max(), after.max()) + 0.018)
    ax.set_ylim(-0.55, len(source) - 0.35)
    clean_axis(ax, left=False, bottom=True)
    ax.tick_params(axis="y", length=0)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markeredgecolor=GRAY, markersize=3.6, label="hepatocyte-adjusted"),
            Line2D([0], [0], marker="D", color="none", markerfacecolor=CYAN, markeredgecolor=CYAN, markersize=3.6, label="all-cell adjusted"),
        ],
        frameon=False,
        ncol=2,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.05),
        borderaxespad=0,
        handletextpad=0.35,
        columnspacing=0.8,
    )
    fig.subplots_adjust(left=0.25, right=0.98, top=0.87, bottom=0.24)
    save_pdf(fig, path)
    return source


def render_cell_context(source: pd.DataFrame, path: Path) -> pd.DataFrame:
    requested = [
        ("hotspot_hepatocytes_f05c535ae5bbc0b9", "Fibroblasts"),
        ("hotspot_hepatocytes_48f39dd4d817a10e", "Cholangiocytes"),
    ]
    parts = []
    for uid, factor in requested:
        part = source[(source["program_uid"] == uid) & (source["cell2location_factor"] == factor)].copy()
        if len(part) != 2:
            raise SpatialResourceError(f"cell-context focus requires two source rows: {uid}/{factor}")
        parts.append(part)
    focused = pd.concat(parts, ignore_index=True)
    dataset_order = ["GSE192741", "Vu_et_al_2025"]
    dataset_labels = ["GSE\n4 donors", "Vu*\n10 arrays"]
    fig, axes = plt.subplots(1, 2, figsize=(4.15, 1.82), sharex=True)
    for ax, (uid, factor), title in zip(
        axes,
        requested,
        ["IGFBP7 / fibroblasts", "BICC1 / cholangiocytes"],
        strict=True,
    ):
        part = focused[(focused["program_uid"] == uid) & (focused["cell2location_factor"] == factor)].set_index("dataset")
        x = np.arange(2)
        med = np.array([part.loc[d, "median_spearman_rho"] for d in dataset_order], dtype=float)
        lo = np.array([part.loc[d, "q25_spearman_rho"] for d in dataset_order], dtype=float)
        hi = np.array([part.loc[d, "q75_spearman_rho"] for d in dataset_order], dtype=float)
        ax.errorbar(x, med, yerr=[med - lo, hi - med], fmt="none", ecolor=GRAY, elinewidth=0.7, capsize=2, capthick=0.7, zorder=1)
        ax.scatter(x[0], med[0], s=25, marker="D", color=CYAN, linewidth=0, zorder=2)
        ax.scatter(x[1], med[1], s=28, marker="^", facecolor=WHITE, edgecolor=CYAN, linewidth=0.8, zorder=2)
        for xx, value in zip(x, med, strict=True):
            x_text = xx + 0.055 if xx == 0 else xx - 0.055
            alignment = "left" if xx == 0 else "right"
            ax.text(x_text, value + 0.035, f"{value:.2f}", ha=alignment, va="bottom")
        ax.axhline(0, color=GRAY_LIGHT, linewidth=0.5)
        ax.set_xticks(x, dataset_labels)
        ax.set_title(title, pad=3)
        ax.set_ylim(-0.05, 0.63)
        clean_axis(ax, left=True, bottom=True)
    axes[0].set_ylabel("Median within-unit Spearman rho")
    axes[1].set_ylabel("")
    fig.subplots_adjust(left=0.15, right=0.99, top=0.86, bottom=0.25, wspace=0.32)
    save_pdf(fig, path)
    return focused


def donor_means(unit: pd.DataFrame) -> pd.DataFrame:
    gse = unit[unit["dataset"] == "GSE192741"].copy()
    return (
        gse.groupby(["program_id", "program_label", "aggregation_group"], observed=True, as_index=False)
        .agg(residual_moran_i=("residual_moran_i", "mean"), n_sections=("sample_id", "nunique"))
        .rename(columns={"aggregation_group": "unit_id"})
    )


def render_unit_heterogeneity(source: pd.DataFrame, path: Path) -> pd.DataFrame:
    donors = donor_means(source)
    vu = source[source["dataset"] == "Vu_et_al_2025"].copy()
    programs = [
        ("hotspot_hepatocytes_f05c535ae5bbc0b9", "IGFBP7"),
        ("hotspot_hepatocytes_48f39dd4d817a10e", "BICC1"),
    ]
    rng = np.random.default_rng(410)
    fig, ax = plt.subplots(figsize=(4.35, 2.25))
    positions = []
    labels = []
    position = 0
    for uid, label in programs:
        for dataset, dataset_label in (("GSE192741", "GSE"), ("Vu_et_al_2025", "Vu*")):
            if dataset == "GSE192741":
                values = donors[donors["program_id"] == uid]["residual_moran_i"].to_numpy(float)
                marker, face, edge = "D", CYAN, CYAN
            else:
                values = vu[vu["program_id"] == uid]["residual_moran_i"].to_numpy(float)
                marker, face, edge = "^", WHITE, CYAN
            jitter = rng.uniform(-0.10, 0.10, len(values))
            ax.scatter(np.full(len(values), position) + jitter, values, s=18, marker=marker, facecolor=face, edgecolor=edge, linewidth=0.7, zorder=2)
            median = float(np.median(values))
            ax.hlines(median, position - 0.18, position + 0.18, color=INK, linewidth=1.0, zorder=3)
            positions.append(position)
            labels.append(f"{label}\n{dataset_label}")
            position += 1
        position += 0.45
    # Preserve the only repeated-section example without labeling every unit.
    h35 = source[(source["dataset"] == "GSE192741") & (source["aggregation_group"] == "H35") & (source["program_id"] == programs[0][0])]
    ax.scatter(np.full(len(h35), positions[0] - 0.24), h35["residual_moran_i"], s=15, facecolor=WHITE, edgecolor=CYAN, linewidth=0.65, zorder=3)
    ax.plot(np.full(len(h35), positions[0] - 0.24), h35["residual_moran_i"], color=GRAY, linewidth=0.5, zorder=1)
    ax.text(positions[0] - 0.24, float(h35["residual_moran_i"].max()) + 0.035, "H35\nsections", ha="center", va="bottom", color=DARK_GRAY)
    ax.axhline(0, color=GRAY, linewidth=0.45, linestyle="--")
    ax.set_xticks(positions, labels)
    ax.set_ylabel("Residual Moran's I")
    ax.set_xlim(-0.55, positions[-1] + 0.45)
    clean_axis(ax, left=True, bottom=True)
    fig.subplots_adjust(left=0.14, right=0.99, top=0.91, bottom=0.24)
    save_pdf(fig, path)
    return source


def render_null(source: pd.DataFrame, path: Path) -> pd.DataFrame:
    source = source.copy()
    order = [
        ("hotspot_hepatocytes_f05c535ae5bbc0b9", "GSE192741"),
        ("hotspot_hepatocytes_f05c535ae5bbc0b9", "Vu_et_al_2025"),
        ("hotspot_hepatocytes_48f39dd4d817a10e", "GSE192741"),
        ("hotspot_hepatocytes_48f39dd4d817a10e", "Vu_et_al_2025"),
    ]
    labels = ["IGFBP7 | GSE", "IGFBP7 | Vu*", "BICC1 | GSE", "BICC1 | Vu*"]
    states = ["supported", "supported*", "indeterminate", "indeterminate*"]
    fig, ax = plt.subplots(figsize=(4.35, 2.15))
    y = np.arange(4)[::-1]
    for yy, key, state in zip(y, order, states, strict=True):
        row = source[(source["program_id"] == key[0]) & (source["dataset"] == key[1])].iloc[0]
        ax.plot([row.q050, row.q950], [yy, yy], color=GRAY, linewidth=3.0, solid_capstyle="butt", zorder=1)
        ax.scatter(row.q500, yy, s=18, marker="D", color=DARK_GRAY, linewidth=0, zorder=2)
        marker = "^" if key[1] == "Vu_et_al_2025" else "o"
        ax.scatter(row.residual_moran_i, yy, s=24, marker=marker, facecolor=WHITE if marker == "^" else CYAN, edgecolor=CYAN, linewidth=0.75, zorder=3)
        ax.text(0.405, yy, state, ha="right", va="center", color=CYAN if "supported" in state else DARK_GRAY)
    ax.axvline(0, color=GRAY, linewidth=0.45, linestyle="--")
    ax.set_yticks(y, labels)
    ax.set_xlabel("Residual Moran's I")
    ax.set_xlim(-0.02, 0.42)
    ax.set_ylim(-0.55, 3.75)
    clean_axis(ax, left=False, bottom=True)
    ax.tick_params(axis="y", length=0)
    ax.legend(
        handles=[
            Line2D([0], [0], color=GRAY, linewidth=3, label="matched-null 5th–95th"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor=CYAN, markeredgecolor=CYAN, markersize=3.8, label="observed"),
        ],
        frameon=False,
        ncol=2,
        loc="upper left",
        bbox_to_anchor=(0, 1.01),
        borderaxespad=0,
        handlelength=1.5,
        columnspacing=0.9,
    )
    fig.subplots_adjust(left=0.24, right=0.98, top=0.82, bottom=0.23)
    save_pdf(fig, path)
    return source


def render_weight(source: pd.DataFrame, path: Path) -> pd.DataFrame:
    order = [
        ("hotspot_hepatocytes_f05c535ae5bbc0b9", "GSE192741"),
        ("hotspot_hepatocytes_f05c535ae5bbc0b9", "Vu_et_al_2025"),
        ("hotspot_hepatocytes_48f39dd4d817a10e", "GSE192741"),
        ("hotspot_hepatocytes_48f39dd4d817a10e", "Vu_et_al_2025"),
    ]
    labels = ["IGFBP7 | GSE", "IGFBP7 | Vu*", "BICC1 | GSE", "BICC1 | Vu*"]
    sensitivity = [
        ("primary_original_l1_weight", "frozen", "o", CYAN),
        ("equal_weight", "equal", "s", WHITE),
        ("leave_highest_weight_gene_out", "leave-top", "^", WHITE),
    ]
    fig, ax = plt.subplots(figsize=(3.95, 2.08))
    y = np.arange(4)[::-1]
    for yy, key in zip(y, order, strict=True):
        part = source[(source["program_id"] == key[0]) & (source["dataset"] == key[1])].set_index("sensitivity_id")
        values = [float(part.loc[sid, "centered_residual_moran_i"]) for sid, _, _, _ in sensitivity]
        ax.plot([min(values), max(values)], [yy, yy], color=GRAY, linewidth=0.6, zorder=1)
        for value, (_, _, marker, face) in zip(values, sensitivity, strict=True):
            ax.scatter(value, yy, s=22, marker=marker, facecolor=face, edgecolor=CYAN, linewidth=0.75, zorder=2)
    ax.axvline(0, color=GRAY, linewidth=0.45, linestyle="--")
    ax.set_yticks(y, labels)
    ax.set_xlabel("Centered Moran's I")
    ax.set_xlim(-0.006, 0.158)
    clean_axis(ax, left=False, bottom=True)
    ax.tick_params(axis="y", length=0)
    handles = [Line2D([0], [0], marker=marker, color="none", markerfacecolor=face, markeredgecolor=CYAN, markersize=3.8, label=label) for _, label, marker, face in sensitivity]
    ax.legend(handles=handles, frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(0, 1.04), borderaxespad=0, handletextpad=0.35, columnspacing=0.8)
    fig.subplots_adjust(left=0.26, right=0.98, top=0.83, bottom=0.24)
    save_pdf(fig, path)
    return source


def render_mask(summary: pd.DataFrame, path: Path) -> pd.DataFrame:
    plot = summary[summary["min_umi"].isin([1000, 2000])].copy()
    programs = [
        ("hotspot_hepatocytes_f05c535ae5bbc0b9", "IGFBP7", CYAN, "o"),
        ("hotspot_hepatocytes_48f39dd4d817a10e", "BICC1", CYAN_MID, "s"),
    ]
    fig, axes = plt.subplots(1, 2, figsize=(4.15, 1.82), sharex=True)
    metrics = [
        ("spearman_rho_with_500", "Rank correlation\nwith 500-UMI mask"),
        ("sign_agreement_with_500_fraction", "Sign agreement\nwith 500-UMI mask"),
    ]
    for ax, (metric, ylabel) in zip(axes, metrics, strict=True):
        for uid, label, color, marker in programs:
            part = plot[plot["program_uid"] == uid].sort_values("min_umi")
            x = part["min_umi"].to_numpy(int)
            y = part[metric].to_numpy(float)
            ax.plot(x, y, color=color, linewidth=0.7, zorder=1)
            ax.scatter(x, y, s=25, marker=marker, facecolor=WHITE if marker == "s" else color, edgecolor=color, linewidth=0.8, label=label, zorder=2)
            for xx, yy in zip(x, y, strict=True):
                ax.text(xx, yy + 0.025, f"{yy:.2f}", ha="center", va="bottom", color=color)
        ax.set_xticks([1000, 2000], ["1,000", "2,000"])
        ax.set_xlabel("Minimum UMI")
        ax.set_ylabel(ylabel)
        ax.set_ylim(0.68, 1.04)
        clean_axis(ax, left=True, bottom=True)
    axes[0].legend(frameon=False, loc="lower left", borderaxespad=0, handletextpad=0.35)
    fig.subplots_adjust(left=0.14, right=0.99, top=0.95, bottom=0.26, wspace=0.34)
    save_pdf(fig, path)
    return summary


def require_ready(project: Path) -> dict[str, Path]:
    roots = {name: project / relative for name, relative in ROOTS.items()}
    for name, root in roots.items():
        ready = root / READY_RELATIVE.get(name, "READY")
        if not ready.is_file():
            raise SpatialResourceError(f"sealed input missing READY: {name}: {root}")
    return roots


def write_table(table: pd.DataFrame, path: Path) -> None:
    table.to_csv(path, sep="\t", index=False)


def build(project: Path, output: Path) -> Path:
    roots = require_ready(project)
    expected_parent = project / "Analysis/Multimodal_Program_Projection/candidates"
    if output.parent.resolve() != expected_parent.resolve() or output.name != RELEASE_ID:
        raise SpatialResourceError("output must be the named isolated publication-figure candidate")
    if output.exists():
        raise SpatialResourceError(f"refusing to overwrite immutable candidate: {output}")
    staging = output.parent / f".{output.name}.incomplete.{os.environ.get('SLURM_JOB_ID', os.getpid())}"
    if staging.exists():
        raise SpatialResourceError(f"staging path already exists: {staging}")
    panel_dir = staging / "panels"
    data_dir = staging / "data"
    panel_dir.mkdir(parents=True)
    data_dir.mkdir()

    fig1 = pd.read_csv(
        roots["downstream"] / "data/fig1_spatial_assay_observability.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    fig5 = pd.read_csv(
        roots["downstream"] / "data/fig5_spatial_passport_examples.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    hmsma = pd.read_csv(roots["downstream"] / "data/figS_hmsma_label_blind_organization.tsv", sep="\t")
    effects = pd.read_parquet(roots["resource"] / "effects/spatial_program_effects.parquet")
    composition = pd.read_csv(roots["composition"] / "figS_full_composition_sensitivity.tsv", sep="\t")
    cell = pd.read_csv(
        roots["cell_context"] / "data/figS_cell_context_association.tsv",
        sep="\t",
        dtype=str,
        keep_default_na=False,
    )
    unit = pd.read_csv(roots["unit_null"] / "data/figS_spatial_unit_heterogeneity.tsv", sep="\t")
    null = pd.read_csv(roots["unit_null"] / "data/figS_spatial_matched_null_calibration.tsv", sep="\t")
    weight = pd.read_csv(roots["unit_null"] / "data/figS_spatial_weight_sensitivity.tsv", sep="\t")
    mask = pd.read_csv(roots["mask"] / "data/hmsma_proxy_mask_summary.tsv", sep="\t")

    outputs: list[tuple[str, pd.DataFrame]] = []
    outputs.append(("fig1_spatial_assay_observability", render_fig1(fig1, panel_dir / PANEL_NAMES[0])))
    outputs.append(("fig4a_spatial_firewall", render_fig4a(panel_dir / PANEL_NAMES[1])))

    # These two panels already match the main-figure grammar; retain them byte-for-byte.
    retained = [
        (
            roots["resource"] / "figures/main/fig4_validation/panels/fig4e_spatial_program_maps.pdf",
            panel_dir / "fig4e_spatial_program_maps.pdf",
        ),
        (
            roots["downstream"] / "panels/figS_spatial_program_coverage.pdf",
            panel_dir / "figS_spatial_program_coverage.pdf",
        ),
    ]
    for source_path, target_path in retained:
        shutil.copy2(source_path, target_path)

    outputs.append(("fig4f_multimodal_program_summary", render_fig4f(effects, panel_dir / PANEL_NAMES[3])))
    outputs.append(("fig5_spatial_passport_examples", render_fig5(fig5, panel_dir / PANEL_NAMES[4])))
    outputs.append(("figS_hmsma_label_blind_organization", render_hmsma_label_blind(hmsma, panel_dir / PANEL_NAMES[6])))
    outputs.append(("figS_full_composition_sensitivity", render_composition(composition, panel_dir / PANEL_NAMES[7])))
    outputs.append(("figS_cell_context_association", render_cell_context(cell, panel_dir / PANEL_NAMES[8])))
    outputs.append(("figS_spatial_unit_heterogeneity", render_unit_heterogeneity(unit, panel_dir / PANEL_NAMES[9])))
    outputs.append(("figS_spatial_matched_null_calibration", render_null(null, panel_dir / PANEL_NAMES[10])))
    outputs.append(("figS_spatial_weight_sensitivity", render_weight(weight, panel_dir / PANEL_NAMES[11])))
    outputs.append(("figS_hmsma_proxy_mask_sensitivity", render_mask(mask, panel_dir / PANEL_NAMES[12])))

    for stem, table in outputs:
        write_table(table, data_dir / f"{stem}.tsv")
    # Whole-family source tables are copied byte-for-byte after rendering so
    # CSV type inference cannot change integer spelling, booleans, or precision.
    exact_table_copies = {
        "figS_hmsma_label_blind_organization.tsv": roots["downstream"] / "data/figS_hmsma_label_blind_organization.tsv",
        "figS_full_composition_sensitivity.tsv": roots["composition"] / "figS_full_composition_sensitivity.tsv",
        "figS_spatial_unit_heterogeneity.tsv": roots["unit_null"] / "data/figS_spatial_unit_heterogeneity.tsv",
        "figS_spatial_matched_null_calibration.tsv": roots["unit_null"] / "data/figS_spatial_matched_null_calibration.tsv",
        "figS_spatial_weight_sensitivity.tsv": roots["unit_null"] / "data/figS_spatial_weight_sensitivity.tsv",
        "figS_hmsma_proxy_mask_sensitivity.tsv": roots["mask"] / "data/hmsma_proxy_mask_summary.tsv",
    }
    for name, source_path in exact_table_copies.items():
        shutil.copy2(source_path, data_dir / name)
    shutil.copy2(
        roots["downstream"] / "data/figS_spatial_program_coverage.tsv",
        data_dir / "figS_spatial_program_coverage.tsv",
    )
    shutil.copy2(
        roots["resource"] / "figures/main/fig4_validation/data/fig4e_spatial_program_maps_selection.tsv",
        data_dir / "fig4e_spatial_program_maps_selection.tsv",
    )

    dispositions = pd.DataFrame(
        [
            ("fig1_spatial_assay_observability", "simplified", "four decision-bearing fields; full joins retained in source table"),
            ("fig4a_spatial_firewall", "simplified", "parallel coverage and inference paths replace prose-heavy cascade"),
            ("fig4e_spatial_program_maps", "retained", "maps already provide one clear assay-native message"),
            ("fig4f_multimodal_program_summary", "simplified", "states written directly in cells; legend removed"),
            ("fig5_spatial_passport_examples", "simplified", "numbers and long caveats moved to source table and caption"),
            ("figS_spatial_program_coverage", "retained", "complete 117-program heatmap already compact"),
            ("figS_hmsma_label_blind_organization", "simplified", "paired distributions retained; footer removed"),
            ("figS_full_composition_sensitivity", "simplified", "short labels and one adjusted value per row"),
            ("figS_cell_context_association", "simplified", "focus on prespecified fibroblast and cholangiocyte contexts"),
            ("figS_spatial_unit_heterogeneity", "simplified", "unit distributions replace four labeled small multiples"),
            ("figS_spatial_matched_null_calibration", "simplified", "direct evidence states replace q-value prose"),
            ("figS_spatial_weight_sensitivity", "simplified", "compact sensitivity comparison without claim-boundary prose"),
            ("figS_hmsma_proxy_mask_sensitivity", "simplified", "two summary stability metrics replace four array scatters"),
        ],
        columns=["panel", "decision", "reason"],
    )
    dispositions.insert(0, "figure_release_id", RELEASE_ID)
    write_table(dispositions, staging / "review_disposition.tsv")

    source_rows = []
    for name, root in roots.items():
        ready = root / READY_RELATIVE.get(name, "READY")
        source_rows.append(
            {
                "figure_release_id": RELEASE_ID,
                "source_role": f"sealed_{name}_ready",
                "relative_path": ready.relative_to(project).as_posix(),
                "bytes": ready.stat().st_size,
                "sha256": sha256_file(ready),
            }
        )
    exact_inputs = [
        (roots["downstream"] / "data/fig1_spatial_assay_observability.tsv", "figure1_exact_source"),
        (roots["downstream"] / "data/fig5_spatial_passport_examples.tsv", "figure5_exact_source"),
        (roots["downstream"] / "data/figS_hmsma_label_blind_organization.tsv", "hmsma_label_blind_exact_source"),
        (roots["downstream"] / "data/figS_spatial_program_coverage.tsv", "coverage_exact_source"),
        (roots["resource"] / "effects/spatial_program_effects.parquet", "figure4f_exact_source"),
        (roots["composition"] / "figS_full_composition_sensitivity.tsv", "composition_exact_source"),
        (roots["cell_context"] / "data/figS_cell_context_association.tsv", "cell_context_complete_source"),
        (roots["unit_null"] / "data/figS_spatial_unit_heterogeneity.tsv", "unit_exact_source"),
        (roots["unit_null"] / "data/figS_spatial_matched_null_calibration.tsv", "matched_null_exact_source"),
        (roots["unit_null"] / "data/figS_spatial_weight_sensitivity.tsv", "weight_sensitivity_exact_source"),
        (roots["mask"] / "data/hmsma_proxy_mask_summary.tsv", "mask_summary_exact_source"),
        (roots["resource"] / "figures/main/fig4_validation/panels/fig4e_spatial_program_maps.pdf", "retained_figure4e_panel"),
        (roots["downstream"] / "panels/figS_spatial_program_coverage.pdf", "retained_coverage_panel"),
        (roots["resource"] / "figures/main/fig4_validation/data/fig4e_spatial_program_maps_selection.tsv", "figure4e_selection_source"),
        (
            project / "Analysis/Multimodal_Program_Projection/candidates/program-context-v2-candidate-2026-08-07/spatial_context/map_source_v2/per_spot_program_map.tsv.gz",
            "figure4e_per_spot_source",
        ),
    ]
    code_dependencies = [
        (Path(__file__).resolve(), "producer"),
        (Path(__file__).resolve().parent / "48_validate_spatial_publication_figures.py", "validator"),
        (Path(__file__).resolve().parent / "run_spatial_publication_figures.sbatch", "execution_wrapper"),
    ]
    for source_path, role in exact_inputs + code_dependencies:
        if not source_path.is_file():
            raise SpatialResourceError(f"manifested input is absent: {source_path}")
        source_rows.append(
            {
                "figure_release_id": RELEASE_ID,
                "source_role": role,
                "relative_path": source_path.relative_to(project).as_posix(),
                "bytes": source_path.stat().st_size,
                "sha256": sha256_file(source_path),
            }
        )
    write_table(pd.DataFrame(source_rows), staging / "source_manifest.tsv")

    panel_rows = []
    for name in PANEL_NAMES:
        panel = panel_dir / name
        panel_rows.append(
            {
                "figure_release_id": RELEASE_ID,
                "panel": Path(name).stem,
                "relative_path": panel.relative_to(staging).as_posix(),
                "bytes": panel.stat().st_size,
                "sha256": sha256_file(panel),
                "canonical_written": "FALSE",
            }
        )
    write_table(pd.DataFrame(panel_rows), staging / "panel_manifest.tsv")
    write_table(
        pd.DataFrame(
            [
                {
                    "figure_release_id": RELEASE_ID,
                    "status": "rendered_awaiting_independent_validation",
                    "n_panels": len(PANEL_NAMES),
                    "n_simplified": int((dispositions["decision"] == "simplified").sum()),
                    "n_retained": int((dispositions["decision"] == "retained").sum()),
                    "canonical_written": "FALSE",
                }
            ]
        ),
        staging / "BUILD_COMPLETE",
    )
    os.replace(staging, output)
    print(f"rendered {len(PANEL_NAMES)} publication-facing spatial panels: {output}", flush=True)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=None)
    args = parser.parse_args()
    project = args.project_root.resolve()
    output = (args.output_root or project / "Analysis/Multimodal_Program_Projection/candidates" / RELEASE_ID).resolve()
    try:
        build(project, output)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
