#!/usr/bin/env python3
# KEY MESSAGE: IGFBP7 spatial organization attenuates but remains supported after all-cell-type adjustment, whereas BICC1 remains indeterminate.
"""Render visually corrected r3 of the full-composition sensitivity panel."""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path


FIGURE_RELEASE_ID = "spatial-full-composition-figure-candidate-2026-08-11-r3"


def import_base(script_dir: Path):
    path = script_dir / "39_render_full_composition_figure.py"
    spec = importlib.util.spec_from_file_location("full_composition_figure_r2_base", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import r2 renderer: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.FIGURE_RELEASE_ID = FIGURE_RELEASE_ID
    return module


BASE = import_base(Path(__file__).resolve().parent)
INPUT_RELEASE_ID = BASE.INPUT_RELEASE_ID
INPUT_READY_SHA256 = BASE.INPUT_READY_SHA256
PANEL_NAME = BASE.PANEL_NAME
SOURCE_NAME = BASE.SOURCE_NAME
sha256_file = BASE.sha256_file
verify_input = BASE.verify_input
build_source = BASE.build_source


def default_input(project_root: Path) -> Path:
    return BASE.default_input(project_root)


def default_output(project_root: Path) -> Path:
    return project_root / "Analysis/Multimodal_Program_Projection/candidates" / FIGURE_RELEASE_ID


def render(source, path: Path) -> None:
    import matplotlib as mpl
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
            "font.size": 6,
            "axes.titlesize": 6,
            "axes.labelsize": 6,
            "xtick.labelsize": 6,
            "ytick.labelsize": 6,
            "legend.fontsize": 6,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, ax = plt.subplots(figsize=(5.4, 2.35))
    y = list(range(len(source) - 1, -1, -1))
    primary = source["primary_centered_moran_i"].to_numpy(float)
    full = source["full_composition_centered_moran_i"].to_numpy(float)
    ax.axvline(0, color="#D7D7D7", linewidth=0.55, zorder=0)
    for position, row, before, after in zip(y, source.itertuples(index=False), primary, full, strict=True):
        source_dependent = row.dataset == "Vu_et_al_2025"
        ax.plot([before, after], [position, position], color="#C9C9C9", linewidth=0.7, zorder=1)
        ax.scatter(
            before,
            position,
            s=20,
            marker="o",
            facecolor="white" if source_dependent else BASE.GRAY,
            edgecolor=BASE.GRAY,
            linewidth=0.7,
            zorder=2,
        )
        ax.scatter(
            after,
            position,
            s=23,
            marker="D",
            facecolor="white" if source_dependent else BASE.CYAN,
            edgecolor=BASE.CYAN,
            linewidth=0.7,
            zorder=3,
        )
        ax.text(before, position + 0.20, f"{before:.3f}", ha="center", va="bottom", color="#666666")
        ax.text(after, position - 0.22, f"{after:.3f}", ha="center", va="top", color=BASE.CYAN)

    labels = [f"{row.program_label}\n{row.dataset_label}" for row in source.itertuples(index=False)]
    ax.set_yticks(y, labels)
    upper = max(primary.max(), full.max()) + 0.022
    ax.set_xlim(-0.005, upper)
    ax.set_xlabel("Centered residual Moran's I relative to matched-gene null")
    ax.set_ylim(-0.65, len(source) - 0.35)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.tick_params(axis="x", width=0.5, length=2)
    ax.legend(
        handles=[
            Line2D([0], [0], marker="o", color="none", markerfacecolor=BASE.GRAY, markeredgecolor=BASE.GRAY, markersize=3.8, label="Hepatocyte + QC adjustment"),
            Line2D([0], [0], marker="D", color="none", markerfacecolor=BASE.CYAN, markeredgecolor=BASE.CYAN, markersize=3.8, label="All 16 cell types + QC"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=BASE.DARK, markersize=3.8, label="Open symbols: source-dependent"),
        ],
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        frameon=False,
        handletextpad=0.4,
        borderaxespad=0,
        borderpad=0,
        labelspacing=0.3,
    )
    fig.text(
        0.99,
        0.015,
        "Attenuation indicates composition-linked context; it is not evidence of mediation or causality.",
        ha="right",
        va="bottom",
        color="#555555",
    )
    # Reserve the right 28% for the legend so no legend glyph or label can
    # enter the data field; retain a separate bottom band for the claim boundary.
    fig.subplots_adjust(left=0.30, right=0.70, top=0.96, bottom=0.24)
    fig.savefig(path, format="pdf", facecolor="white")
    plt.close(fig)


def run(project_root: Path, input_root: Path | None, output_root: Path | None) -> Path:
    import pandas as pd

    root = project_root.resolve()
    input_path = (input_root or default_input(root)).resolve()
    output_path = (output_root or default_output(root)).resolve()
    candidate_parent = (root / "Analysis/Multimodal_Program_Projection/candidates").resolve()
    if output_path.parent != candidate_parent or output_path.name != FIGURE_RELEASE_ID:
        raise RuntimeError("figure output must be the named isolated r3 candidate root")
    if output_path.exists():
        raise RuntimeError(f"refusing to overwrite figure candidate: {output_path}")
    verify_input(input_path)
    staging = output_path.parent / f".{output_path.name}.incomplete.{os.environ.get('SLURM_JOB_ID', os.getpid())}"
    if staging.exists():
        raise RuntimeError(f"staging path already exists: {staging}")
    staging.mkdir()
    source = build_source(input_path)
    columns = [
        "figure_release_id", "release_id", "program_id", "program_label", "dataset",
        "dataset_label", "display_order", "primary_adjustment", "full_adjustment",
        "primary_centered_moran_i", "full_composition_centered_moran_i", "centered_change",
        "primary_residual_pvalue", "primary_residual_padj", "residual_pvalue",
        "residual_padj", "primary_robust", "robust", "evidence_state", "interpretation",
        "matched_set_sha256", "matched_set_hash_identical", "claim_boundary",
    ]
    source[columns].to_csv(staging / SOURCE_NAME, sep="\t", index=False)
    render(source, staging / PANEL_NAME)
    pd.DataFrame(
        [{
            "figure_release_id": FIGURE_RELEASE_ID,
            "panel": PANEL_NAME,
            "panel_sha256": sha256_file(staging / PANEL_NAME),
            "source_table": SOURCE_NAME,
            "source_sha256": sha256_file(staging / SOURCE_NAME),
            "input_release_id": INPUT_RELEASE_ID,
            "input_ready_sha256": INPUT_READY_SHA256,
            "producer_sha256": sha256_file(Path(__file__).resolve()),
            "supersedes_visual_release": "spatial-full-composition-figure-candidate-2026-08-11-r2",
            "supersession_reason": "legend_removed_from_data_field_after_200dpi_visual_review",
            "canonical_written": "FALSE",
        }]
    ).to_csv(staging / "panel_manifest.tsv", sep="\t", index=False)
    os.replace(staging, output_path)
    print(f"full-composition supplementary figure r3 rendered: {output_path}", flush=True)
    return output_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--input-root", type=Path, default=None)
    parser.add_argument("--output-root", type=Path, default=None)
    args = parser.parse_args()
    try:
        run(args.project_root, args.input_root, args.output_root)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
