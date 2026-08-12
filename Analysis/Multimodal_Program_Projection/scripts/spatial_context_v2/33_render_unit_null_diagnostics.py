#!/usr/bin/env python3
# KEY MESSAGE: positive unit-level spatial organization is not sufficient evidence;
# donor/array aggregation and a matched-gene null distinguish support from indeterminacy.
"""Render candidate-only spatial unit, matched-null, and weight-sensitivity diagnostics."""

from __future__ import annotations

import argparse
import os
import shutil
import tempfile
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
    "legend.fontsize": 6,
    "figure.titlesize": 6,
})
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from spatial_resource_lib import PROGRAM_RELEASE_ID, SpatialResourceError, sha256_file, write_tsv


DIAGNOSTIC_RELEASE_ID = "spatial-unit-null-diagnostics-candidate-2026-08-11-r2"
CYAN = "#007C91"
CYAN_MID = "#68B7C7"
CYAN_LIGHT = "#CDE7EE"
GRAY = "#9E9E9E"
GRAY_LIGHT = "#E6E6E5"
DARK_GRAY = "#4D4D4D"
INK = "#222222"
WHITE = "#FFFFFF"

DATASETS = ["GSE192741", "Vu_et_al_2025"]
PROGRAMS = [
    "hotspot_hepatocytes_f05c535ae5bbc0b9",
    "hotspot_hepatocytes_48f39dd4d817a10e",
]
PROGRAM_LABELS = {
    PROGRAMS[0]: "Stromal ECM (IGFBP7)",
    PROGRAMS[1]: "Ductular injury (BICC1)",
}
GSE_SECTION_ORDER = ["JBO014", "JBO015", "JBO018", "JBO019", "JBO022"]
VU_ARRAY_ORDER = [
    "VLP115_A", "VLP115_D", "VLP116_A", "VLP116_D", "VLP119_A",
    "VLP119_D", "VLP120_A", "VLP120_D", "VLP121_A", "VLP121_D",
]
SENSITIVITY_ORDER = [
    "primary_original_l1_weight",
    "equal_weight",
    "leave_highest_weight_gene_out",
]
SENSITIVITY_LABELS = {
    "primary_original_l1_weight": "frozen weights",
    "equal_weight": "equal weights",
    "leave_highest_weight_gene_out": "leave top-weighted gene out",
}


def require_upstream(upstream: Path) -> None:
    ready = upstream / "READY"
    if not ready.is_file():
        raise SpatialResourceError(f"sealed v2 candidate is missing READY: {ready}")
    seal = pd.read_csv(ready, sep="\t", dtype=str, keep_default_na=False)
    if len(seal) != 1 or seal.iloc[0].get("status") != "pass_v2_candidate":
        raise SpatialResourceError("unit/null diagnostics require the independently validated v2 candidate")
    if seal.iloc[0].get("release_id") != PROGRAM_RELEASE_ID:
        raise SpatialResourceError("upstream program release ID does not match the frozen spatial contract")
    required = (
        "spatial_section_results.tsv",
        "spatial_program_results.tsv",
        "spatial_null_summary.tsv",
        "spatial_sensitivity_audit.tsv",
    )
    for filename in required:
        if not (upstream / filename).is_file():
            raise SpatialResourceError(f"upstream v2 candidate is incomplete: missing {filename}")


def save_pdf(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight", dpi=400)
    plt.close(fig)


def load_inputs(upstream: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    section = pd.read_csv(upstream / "spatial_section_results.tsv", sep="\t")
    program = pd.read_csv(upstream / "spatial_program_results.tsv", sep="\t")
    null = pd.read_csv(upstream / "spatial_null_summary.tsv", sep="\t")
    sensitivity = pd.read_csv(upstream / "spatial_sensitivity_audit.tsv", sep="\t")
    expected_pairs = {(dataset, uid) for uid in PROGRAMS for dataset in DATASETS}
    for name, table in (("section", section), ("program", program), ("sensitivity", sensitivity)):
        pairs = set(zip(table["dataset"], table["program_id"]))
        if pairs != expected_pairs:
            raise SpatialResourceError(f"{name} table does not contain the complete two-program by two-dataset family")
    residual_null = null[null["statistic"] == "residual_moran_i"]
    if set(zip(residual_null["dataset"], residual_null["program_id"])) != expected_pairs:
        raise SpatialResourceError("residual matched-null table does not contain the complete family")
    expected_units = {
        "GSE192741": set(GSE_SECTION_ORDER),
        "Vu_et_al_2025": set(VU_ARRAY_ORDER),
    }
    for dataset, unit_ids in expected_units.items():
        part = section[section["dataset"] == dataset]
        if set(part["sample_id"]) != unit_ids or len(part) != len(unit_ids) * len(PROGRAMS):
            raise SpatialResourceError(f"unexpected physical-unit family for {dataset}")
    return section, program, null, sensitivity


def build_unit_source(section: pd.DataFrame) -> pd.DataFrame:
    source = section.copy()
    source["program_label"] = source["program_id"].map(PROGRAM_LABELS)
    source["physical_unit_type"] = np.where(source["dataset"] == "GSE192741", "Visium_section", "Visium_array")
    source["aggregation_group"] = np.where(source["dataset"] == "GSE192741", source["individual"], source["sample_id"])
    source["aggregation_group_resolution"] = np.where(source["dataset"] == "GSE192741", "resolved_donor", "unresolved_physical_unit")
    source["source_dependent"] = np.where(source["dataset"] == "Vu_et_al_2025", "TRUE", "FALSE")
    source["inferential_pvalue_authorized"] = "FALSE"
    source["collapsed_residual_moran_i"] = source.groupby(
        ["dataset", "program_id", "aggregation_group"], observed=True
    )["residual_moran_i"].transform("mean")
    source["collapse_rule"] = np.where(
        source["dataset"] == "GSE192741",
        "mean_sections_within_donor_then_mean_four_donors",
        "mean_ten_physical_arrays_donor_mapping_unresolved",
    )
    source.insert(0, "upstream_release_id", PROGRAM_RELEASE_ID)
    source.insert(0, "diagnostic_release_id", DIAGNOSTIC_RELEASE_ID)
    order = {sample: index for index, sample in enumerate(GSE_SECTION_ORDER + VU_ARRAY_ORDER)}
    source["_program_order"] = source["program_id"].map({uid: index for index, uid in enumerate(PROGRAMS)})
    source["_dataset_order"] = source["dataset"].map({dataset: index for index, dataset in enumerate(DATASETS)})
    source["_unit_order"] = source["sample_id"].map(order)
    source = source.sort_values(["_program_order", "_dataset_order", "_unit_order"]).drop(
        columns=["_program_order", "_dataset_order", "_unit_order"]
    )
    return source


def panel_unit_heterogeneity(source: pd.DataFrame, panel_dir: Path, data_dir: Path) -> Path:
    source.to_csv(data_dir / "figS_spatial_unit_heterogeneity.tsv", sep="\t", index=False)
    fig, axes = plt.subplots(2, 2, figsize=(5.5, 3.75), sharey=True)
    ymax = max(0.55, float(source["residual_moran_i"].max()) + 0.03)
    ymin = min(-0.035, float(source["residual_moran_i"].min()) - 0.02)
    for row_index, uid in enumerate(PROGRAMS):
        for column_index, dataset in enumerate(DATASETS):
            ax = axes[row_index, column_index]
            part = source[(source["program_id"] == uid) & (source["dataset"] == dataset)].copy()
            if dataset == "GSE192741":
                part = part.set_index("sample_id").loc[GSE_SECTION_ORDER].reset_index()
                x = np.arange(len(part), dtype=float)
                ax.scatter(x, part["residual_moran_i"], s=17, color=CYAN_LIGHT, edgecolor=CYAN, linewidth=0.45, zorder=3)
                donor_centers: list[float] = []
                donor_values: list[float] = []
                for donor, donor_part in part.groupby("aggregation_group", sort=False):
                    indices = np.flatnonzero(part["aggregation_group"].to_numpy() == donor).astype(float)
                    values = donor_part["residual_moran_i"].to_numpy(dtype=float)
                    if len(indices) > 1:
                        ax.plot(indices, values, color=GRAY, linewidth=0.55, zorder=2)
                    donor_centers.append(float(indices.mean()))
                    donor_values.append(float(donor_part["collapsed_residual_moran_i"].iloc[0]))
                ax.scatter(donor_centers, donor_values, s=21, marker="D", color=CYAN, edgecolor=WHITE, linewidth=0.35, zorder=4)
                labels = ["H35\nJBO014", "H35\nJBO015", "H36\nJBO018", "H37\nJBO019", "H38\nJBO022"]
                ax.set_xticks(x, labels, rotation=55, ha="right", rotation_mode="anchor")
                h35_mean = float(part.loc[part["aggregation_group"] == "H35", "collapsed_residual_moran_i"].iloc[0])
                ax.annotate(
                    "H35 donor mean", xy=(0.5, h35_mean), xytext=(0, -9 if row_index == 0 else 8),
                    textcoords="offset points", ha="center", va="center", color=DARK_GRAY,
                )
                dataset_label = "GSE192741 | 5 sections / 4 donors"
            else:
                part = part.set_index("sample_id").loc[VU_ARRAY_ORDER].reset_index()
                x = np.arange(len(part), dtype=float)
                ax.scatter(x, part["residual_moran_i"], s=17, marker="^", color=CYAN, edgecolor=WHITE, linewidth=0.35, zorder=3)
                ax.set_xticks(x, [value.replace("VLP", "") for value in VU_ARRAY_ORDER], rotation=55, ha="right", rotation_mode="anchor")
                dataset_label = "Vu et al. | 10 physical arrays (source-dependent)"
            ax.axhline(0, color=GRAY, linewidth=0.45, linestyle="--", zorder=0)
            ax.set_ylim(ymin, ymax)
            ax.set_title(dataset_label if row_index == 0 else "", pad=3, fontweight="normal")
            if column_index == 0:
                ax.set_ylabel(f"{PROGRAM_LABELS[uid]}\nresidual Moran's I")
            else:
                ax.set_ylabel("")
            ax.spines["top"].set_visible(False)
            ax.spines["right"].set_visible(False)
            ax.spines["left"].set_linewidth(0.35)
            ax.spines["bottom"].set_linewidth(0.35)
            ax.tick_params(width=0.35, length=2, pad=1.5)
    handles = [
        Line2D([], [], marker="o", linestyle="none", markerfacecolor=CYAN_LIGHT, markeredgecolor=CYAN, markersize=4, label="section"),
        Line2D([], [], marker="D", linestyle="none", markerfacecolor=CYAN, markeredgecolor=WHITE, markersize=4, label="donor-collapsed value"),
        Line2D([], [], marker="^", linestyle="none", markerfacecolor=CYAN, markeredgecolor=WHITE, markersize=4, label="physical array; donor key unresolved"),
    ]
    fig.legend(handles=handles, frameon=False, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 0.075), columnspacing=0.9, handletextpad=0.35)
    fig.text(0.01, 0.015, "Physical-unit estimates are descriptive; no section- or array-level P values.", ha="left", va="bottom", color=DARK_GRAY)
    fig.subplots_adjust(left=0.14, right=0.99, top=0.93, bottom=0.30, wspace=0.18, hspace=0.28)
    path = panel_dir / "figS_spatial_unit_heterogeneity.pdf"
    save_pdf(fig, path)
    return path


def build_null_source(program: pd.DataFrame, null: pd.DataFrame, sensitivity: pd.DataFrame) -> pd.DataFrame:
    residual_null = null[null["statistic"] == "residual_moran_i"].copy()
    program_columns = [
        "program_id", "dataset", "program_name", "residual_moran_i", "residual_pvalue",
        "residual_padj", "residual_null_mean", "residual_null_sd", "matched_set_sha256",
    ]
    source = program[program_columns].merge(
        residual_null[[
            "program_id", "dataset", "n_null", "matched_set_sha256", "null_mean", "null_sd",
            "q010", "q050", "q500", "q950", "q990",
        ]],
        on=["program_id", "dataset", "matched_set_sha256"],
        validate="one_to_one",
    )
    primary = sensitivity[sensitivity["sensitivity_id"] == "primary_original_l1_weight"][[
        "program_id", "dataset", "centered_residual_moran_i"
    ]].rename(columns={"centered_residual_moran_i": "sensitivity_primary_centered_residual_moran_i"})
    source = source.merge(primary, on=["program_id", "dataset"], validate="one_to_one")
    source["centered_residual_moran_i"] = source["residual_moran_i"] - source["null_mean"]
    source["aggregation_level"] = np.where(source["dataset"] == "GSE192741", "four_donor_mean", "ten_physical_array_mean")
    source["n_aggregation_units"] = np.where(source["dataset"] == "GSE192741", 4, 10)
    source["n_technical_units"] = np.where(source["dataset"] == "GSE192741", 5, 10)
    source["source_dependent"] = np.where(source["dataset"] == "Vu_et_al_2025", "TRUE", "FALSE")
    source["inferential_level"] = np.where(source["dataset"] == "GSE192741", "dataset_donor_first", "dataset_physical_array_source_dependent")
    source["uncertainty_semantics"] = "matched_gene_null_quantiles_not_sampling_confidence_interval"
    source.insert(0, "upstream_release_id", PROGRAM_RELEASE_ID)
    source.insert(0, "diagnostic_release_id", DIAGNOSTIC_RELEASE_ID)
    source["_program_order"] = source["program_id"].map({uid: index for index, uid in enumerate(PROGRAMS)})
    source["_dataset_order"] = source["dataset"].map({dataset: index for index, dataset in enumerate(DATASETS)})
    source = source.sort_values(["_program_order", "_dataset_order"]).drop(columns=["_program_order", "_dataset_order"])
    return source


def panel_matched_null(source: pd.DataFrame, panel_dir: Path, data_dir: Path) -> Path:
    source.to_csv(data_dir / "figS_spatial_matched_null_calibration.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(5.25, 2.45))
    y = np.arange(len(source), dtype=float)
    for row_index, (_, row) in enumerate(source.iterrows()):
        ax.plot([row["q010"], row["q990"]], [row_index, row_index], color=GRAY_LIGHT, linewidth=1.1, solid_capstyle="round", zorder=1)
        ax.plot([row["q050"], row["q950"]], [row_index, row_index], color=GRAY, linewidth=3.1, solid_capstyle="round", zorder=2)
        ax.scatter(row["q500"], row_index, s=17, marker="D", color=DARK_GRAY, edgecolor=WHITE, linewidth=0.3, zorder=3)
        marker = "o" if row["dataset"] == "GSE192741" else "^"
        ax.scatter(row["residual_moran_i"], row_index, s=27, marker=marker, color=CYAN, edgecolor=WHITE, linewidth=0.35, zorder=4)
        q_text = f"delta={row['centered_residual_moran_i']:.3f}; q={row['residual_padj']:.4f}"
        ax.text(0.305, row_index, q_text, ha="left", va="center", color=INK)
    labels = [
        f"{PROGRAM_LABELS[row.program_id]} | {'GSE' if row.dataset == 'GSE192741' else 'Vu (S)'}"
        for row in source.itertuples(index=False)
    ]
    ax.set_yticks(y, labels)
    ax.set_ylim(3.45, -0.62)
    ax.set_xlim(-0.02, 0.42)
    ax.set_xlabel("abundance-adjusted residual Moran's I")
    ax.axvline(0, color=GRAY, linewidth=0.45, linestyle="--", zorder=0)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.35)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.tick_params(axis="x", width=0.35, length=2)
    handles = [
        Line2D([], [], color=GRAY, linewidth=3.1, label="matched-null 5th-95th percentile"),
        Line2D([], [], marker="D", linestyle="none", color=DARK_GRAY, markersize=4, label="matched-null median"),
        Line2D([], [], marker="o", linestyle="none", color=CYAN, markersize=4, label="observed donor/array-collapsed statistic"),
    ]
    ax.legend(handles=handles, frameon=False, ncol=1, loc="lower left", bbox_to_anchor=(0, -0.47), borderaxespad=0, handlelength=1.7)
    ax.text(0.42, -0.50, "S = source-dependent. Null quantiles are not sampling confidence intervals.", ha="right", va="center", color=DARK_GRAY)
    fig.subplots_adjust(left=0.31, right=0.99, top=0.98, bottom=0.32)
    path = panel_dir / "figS_spatial_matched_null_calibration.pdf"
    save_pdf(fig, path)
    return path


def build_sensitivity_source(sensitivity: pd.DataFrame) -> pd.DataFrame:
    source = sensitivity.copy()
    source["program_label"] = source["program_id"].map(PROGRAM_LABELS)
    source["sensitivity_label"] = source["sensitivity_id"].map(SENSITIVITY_LABELS)
    source["source_dependent"] = np.where(source["dataset"] == "Vu_et_al_2025", "TRUE", "FALSE")
    source["inferential_pvalue_authorized"] = "FALSE"
    source.insert(0, "diagnostic_release_id", DIAGNOSTIC_RELEASE_ID)
    source["_program_order"] = source["program_id"].map({uid: index for index, uid in enumerate(PROGRAMS)})
    source["_dataset_order"] = source["dataset"].map({dataset: index for index, dataset in enumerate(DATASETS)})
    source["_sensitivity_order"] = source["sensitivity_id"].map({value: index for index, value in enumerate(SENSITIVITY_ORDER)})
    if source["_sensitivity_order"].isna().any():
        raise SpatialResourceError("unexpected sensitivity analysis entered the diagnostic family")
    source = source.sort_values(["_program_order", "_dataset_order", "_sensitivity_order"]).drop(
        columns=["_program_order", "_dataset_order", "_sensitivity_order"]
    )
    return source


def panel_sensitivity(source: pd.DataFrame, panel_dir: Path, data_dir: Path) -> Path:
    source.to_csv(data_dir / "figS_spatial_weight_sensitivity.tsv", sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(4.7, 2.45))
    combinations = [(uid, dataset) for uid in PROGRAMS for dataset in DATASETS]
    shapes = {SENSITIVITY_ORDER[0]: "o", SENSITIVITY_ORDER[1]: "s", SENSITIVITY_ORDER[2]: "^"}
    fills = {SENSITIVITY_ORDER[0]: CYAN, SENSITIVITY_ORDER[1]: WHITE, SENSITIVITY_ORDER[2]: CYAN_LIGHT}
    offsets = {SENSITIVITY_ORDER[0]: -0.14, SENSITIVITY_ORDER[1]: 0.0, SENSITIVITY_ORDER[2]: 0.14}
    for row_index, (uid, dataset) in enumerate(combinations):
        part = source[(source["program_id"] == uid) & (source["dataset"] == dataset)].set_index("sensitivity_id").loc[SENSITIVITY_ORDER]
        values = part["centered_residual_moran_i"].to_numpy(dtype=float)
        ax.plot(values, row_index + np.array([offsets[value] for value in SENSITIVITY_ORDER]), color=GRAY, linewidth=0.55, zorder=1)
        for sensitivity_id, value in zip(SENSITIVITY_ORDER, values):
            ax.scatter(
                value, row_index + offsets[sensitivity_id], s=24, marker=shapes[sensitivity_id],
                facecolor=fills[sensitivity_id], edgecolor=CYAN, linewidth=0.55, zorder=2,
            )
    labels = [
        f"{PROGRAM_LABELS[uid]} | {'GSE' if dataset == 'GSE192741' else 'Vu (S)'}"
        for uid, dataset in combinations
    ]
    ax.set_yticks(range(len(combinations)), labels)
    ax.set_ylim(3.45, -0.62)
    ax.axvline(0, color=GRAY, linewidth=0.5, linestyle="--", zorder=0)
    ax.set_xlim(-0.008, 0.155)
    ax.set_xlabel("observed - matched-null mean residual Moran's I")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_linewidth(0.35)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.tick_params(axis="x", width=0.35, length=2)
    handles = [
        Line2D([], [], marker=shapes[value], linestyle="none", markerfacecolor=fills[value], markeredgecolor=CYAN, markersize=4, label=SENSITIVITY_LABELS[value])
        for value in SENSITIVITY_ORDER
    ]
    ax.legend(handles=handles, frameon=False, ncol=1, loc="lower left", bbox_to_anchor=(0, -0.48), borderaxespad=0, handletextpad=0.4)
    ax.text(0.155, -0.50, "S = source-dependent. Sensitivities test direction, not a new inferential family.", ha="right", va="center", color=DARK_GRAY)
    fig.subplots_adjust(left=0.34, right=0.99, top=0.98, bottom=0.32)
    path = panel_dir / "figS_spatial_weight_sensitivity.pdf"
    save_pdf(fig, path)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    upstream = args.upstream_root.resolve()
    output = args.output_root.resolve()
    if output.exists():
        raise SpatialResourceError(f"immutable diagnostic candidate already exists: {output}")
    require_upstream(upstream)
    section, program, null, sensitivity = load_inputs(upstream)
    temp = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        panel_dir = temp / "panels"
        data_dir = temp / "data"
        panel_dir.mkdir(parents=True)
        data_dir.mkdir(parents=True)
        unit_source = build_unit_source(section)
        null_source = build_null_source(program, null, sensitivity)
        sensitivity_source = build_sensitivity_source(sensitivity)
        panels = [
            panel_unit_heterogeneity(unit_source, panel_dir, data_dir),
            panel_matched_null(null_source, panel_dir, data_dir),
            panel_sensitivity(sensitivity_source, panel_dir, data_dir),
        ]
        upstream_rows = []
        for filename in (
            "READY", "spatial_section_results.tsv", "spatial_program_results.tsv",
            "spatial_null_summary.tsv", "spatial_sensitivity_audit.tsv",
        ):
            path = upstream / filename
            upstream_rows.append({
                "diagnostic_release_id": DIAGNOSTIC_RELEASE_ID,
                "upstream_release_id": PROGRAM_RELEASE_ID,
                "source_path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            })
        write_tsv(temp / "upstream_source_manifest.tsv", tuple(upstream_rows[0]), upstream_rows)
        write_tsv(
            temp / "BUILD_COMPLETE",
            (
                "diagnostic_release_id", "upstream_release_id", "status", "n_panels",
                "upstream_ready_sha256", "upstream_source_manifest_sha256", "canonical_figure_written",
            ),
            [{
                "diagnostic_release_id": DIAGNOSTIC_RELEASE_ID,
                "upstream_release_id": PROGRAM_RELEASE_ID,
                "status": "rendered_candidate_awaiting_independent_validation",
                "n_panels": len(panels),
                "upstream_ready_sha256": sha256_file(upstream / "READY"),
                "upstream_source_manifest_sha256": sha256_file(temp / "upstream_source_manifest.tsv"),
                "canonical_figure_written": "FALSE",
            }],
        )
        os.replace(temp, output)
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise
    print(f"WROTE {len(panels)} candidate-only spatial unit/null diagnostic panels to {output}")


if __name__ == "__main__":
    main()
