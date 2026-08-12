#!/usr/bin/env python3
# KEY MESSAGE: IGFBP7 spatial organization attenuates but remains supported after all-cell-type adjustment, whereas BICC1 remains indeterminate.
"""Render the candidate-only supplementary full-composition sensitivity panel."""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path


INPUT_RELEASE_ID = "spatial-full-composition-sensitivity-candidate-2026-08-11-r2"
INPUT_READY_SHA256 = "e15a92d129f26b913556016246bbb8dc40aef7045f318fadfb9f1756dcd83ae4"
FIGURE_RELEASE_ID = "spatial-full-composition-figure-candidate-2026-08-11-r2"
PANEL_NAME = "figS_full_composition_sensitivity.pdf"
SOURCE_NAME = "figS_full_composition_sensitivity.tsv"
PROGRAM_LABELS = {
    "hotspot_hepatocytes_f05c535ae5bbc0b9": "Stromal ECM (IGFBP7)",
    "hotspot_hepatocytes_48f39dd4d817a10e": "Ductular injury (BICC1)",
}
PROGRAM_ORDER = list(PROGRAM_LABELS)
DATASET_ORDER = ["GSE192741", "Vu_et_al_2025"]
GRAY = "#9E9E9E"
CYAN = "#007C91"
DARK = "#252525"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def default_input(project_root: Path) -> Path:
    return project_root / "Analysis/Multimodal_Program_Projection/candidates" / INPUT_RELEASE_ID


def default_output(project_root: Path) -> Path:
    return project_root / "Analysis/Multimodal_Program_Projection/candidates" / FIGURE_RELEASE_ID


def verify_input(root: Path) -> None:
    ready = root / "READY"
    if not ready.is_file() or sha256_file(ready) != INPUT_READY_SHA256:
        raise RuntimeError("validated r2 full-composition READY is missing or drifted")


def build_source(input_root: Path):
    import pandas as pd

    comparison = pd.read_csv(input_root / "comparison_to_primary.tsv", sep="\t")
    if len(comparison) != 4 or comparison.duplicated(["dataset", "program_id"]).any():
        raise RuntimeError("expected four unique dataset/program sensitivity comparisons")
    expected = {(dataset, program) for program in PROGRAM_ORDER for dataset in DATASET_ORDER}
    if set(zip(comparison["dataset"], comparison["program_id"], strict=True)) != expected:
        raise RuntimeError("sensitivity comparison is not the complete two-program/two-dataset family")
    if not comparison["matched_set_hash_identical"].astype(str).str.lower().isin({"true", "1"}).all():
        raise RuntimeError("figure input does not preserve all four matched-gene sets")
    comparison["figure_release_id"] = FIGURE_RELEASE_ID
    comparison["program_label"] = comparison["program_id"].map(PROGRAM_LABELS)
    comparison["dataset_label"] = comparison["dataset"].map(
        {"GSE192741": "GSE192741 (4 donors)", "Vu_et_al_2025": "Vu (10 arrays; source-dependent)"}
    )
    comparison["primary_adjustment"] = "hepatocyte_q05_plus_log_library_and_detected_genes"
    comparison["full_adjustment"] = "all_16_cell2location_q05_plus_log_library_and_detected_genes"
    comparison["interpretation"] = comparison.apply(
        lambda row: (
            "supported_after_full_composition_adjustment"
            if row["dataset"] == "GSE192741" and bool(row["robust"])
            else "within_source_supported_but_source_dependent"
            if row["dataset"] == "Vu_et_al_2025" and bool(row["robust"])
            else "indeterminate"
            if row["dataset"] == "GSE192741"
            else "indeterminate_within_source_and_source_dependent"
        ),
        axis=1,
    )
    comparison["claim_boundary"] = (
        "attenuation_is_composition_linked_context_not_mediation_causality_or_cell_intrinsic_proof"
    )
    order = {(program, dataset): i for i, (program, dataset) in enumerate(
        (pair for program in PROGRAM_ORDER for pair in ((program, "GSE192741"), (program, "Vu_et_al_2025")))
    )}
    comparison["display_order"] = [order[(program, dataset)] for program, dataset in zip(comparison["program_id"], comparison["dataset"], strict=True)]
    return comparison.sort_values("display_order").reset_index(drop=True)


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
    fig, ax = plt.subplots(figsize=(4.7, 2.35))
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
            facecolor="white" if source_dependent else GRAY,
            edgecolor=GRAY,
            linewidth=0.7,
            zorder=2,
        )
        ax.scatter(
            after,
            position,
            s=23,
            marker="D",
            facecolor="white" if source_dependent else CYAN,
            edgecolor=CYAN,
            linewidth=0.7,
            zorder=3,
        )
        ax.text(before, position + 0.20, f"{before:.3f}", ha="center", va="bottom", color="#666666")
        ax.text(after, position - 0.22, f"{after:.3f}", ha="center", va="top", color=CYAN)

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
            Line2D([0], [0], marker="o", color="none", markerfacecolor=GRAY, markeredgecolor=GRAY, markersize=3.8, label="Hepatocyte + QC adjustment"),
            Line2D([0], [0], marker="D", color="none", markerfacecolor=CYAN, markeredgecolor=CYAN, markersize=3.8, label="All 16 cell types + QC"),
            Line2D([0], [0], marker="o", color="none", markerfacecolor="white", markeredgecolor=DARK, markersize=3.8, label="Open symbols: source-dependent"),
        ],
        loc="upper right",
        frameon=False,
        handletextpad=0.4,
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
    fig.subplots_adjust(left=0.34, right=0.98, top=0.96, bottom=0.24)
    fig.savefig(path, format="pdf", facecolor="white")
    plt.close(fig)


def run(project_root: Path, input_root: Path | None, output_root: Path | None) -> Path:
    root = project_root.resolve()
    input_path = (input_root or default_input(root)).resolve()
    output_path = (output_root or default_output(root)).resolve()
    candidate_parent = (root / "Analysis/Multimodal_Program_Projection/candidates").resolve()
    if output_path.parent != candidate_parent or output_path.name != FIGURE_RELEASE_ID:
        raise RuntimeError("figure output must be the named isolated candidate root")
    if output_path.exists():
        raise RuntimeError(f"refusing to overwrite figure candidate: {output_path}")
    verify_input(input_path)
    staging = output_path.parent / f".{output_path.name}.incomplete.{os.environ.get('SLURM_JOB_ID', os.getpid())}"
    if staging.exists():
        raise RuntimeError(f"staging path already exists: {staging}")
    staging.mkdir()
    try:
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
        manifest = [
            {
                "figure_release_id": FIGURE_RELEASE_ID,
                "panel": PANEL_NAME,
                "panel_sha256": sha256_file(staging / PANEL_NAME),
                "source_table": SOURCE_NAME,
                "source_sha256": sha256_file(staging / SOURCE_NAME),
                "input_release_id": INPUT_RELEASE_ID,
                "input_ready_sha256": INPUT_READY_SHA256,
                "producer_sha256": sha256_file(Path(__file__).resolve()),
                "canonical_written": "FALSE",
            }
        ]
        import pandas as pd

        pd.DataFrame(manifest).to_csv(staging / "panel_manifest.tsv", sep="\t", index=False)
        os.replace(staging, output_path)
    except Exception:
        raise
    print(f"full-composition supplementary figure rendered: {output_path}", flush=True)
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
