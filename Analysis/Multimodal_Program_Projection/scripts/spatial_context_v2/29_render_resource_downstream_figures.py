#!/usr/bin/env python3
# KEY MESSAGE: spatial assays differ in what they can observe and infer; technical
# coverage must remain distinct from biological replication and negative evidence.
"""Render candidate-only Figure 1, Figure 5, and supplement spatial panels."""

from __future__ import annotations

import argparse
import json
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
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np
import pandas as pd

from spatial_resource_lib import RESOURCE_RELEASE_ID, SpatialResourceError, sha256_file, write_tsv


INTEGRATION_RELEASE_ID = "spatial-figure-integration-candidate-2026-08-11-r2"
CYAN = "#007C91"
CYAN_MID = "#68B7C7"
CYAN_LIGHT = "#CDE7EE"
GRAY = "#9E9E9E"
GRAY_LIGHT = "#E6E6E5"
INK = "#222222"
WHITE = "#FFFFFF"
DARK_GRAY = "#4D4D4D"

COVERAGE_DATASETS = [
    "GSE192741",
    "Vu_et_al_2025",
    "HRA007511_HMSMA",
    "Yakubovsky2026",
    "Govaere2026_CosMx",
    "Govaere2026_GeoMx",
    "GSE287826",
]
DATASET_LABELS = {
    "GSE192741": "GSE192741",
    "Vu_et_al_2025": "Vu et al.",
    "HRA007511_HMSMA": "HMSMA",
    "Yakubovsky2026": "Yakubovsky",
    "Govaere2026_CosMx": "CosMx",
    "Govaere2026_GeoMx": "GeoMx",
    "GSE287826": "GSE287826",
}
CONFIRMATORY_UIDS = {
    "hotspot_hepatocytes_f05c535ae5bbc0b9": "IGFBP7",
    "hotspot_hepatocytes_48f39dd4d817a10e": "BICC1",
}


def require_ready(upstream: Path) -> None:
    for relative in ("registry/READY", "coverage/READY", "effects/READY", "validation/READY", "figures/READY"):
        if not (upstream / relative).is_file():
            raise SpatialResourceError(f"upstream candidate is not sealed: missing {relative}")
    validation = pd.read_csv(upstream / "validation/READY", sep="\t", dtype=str)
    if len(validation) != 1 or validation.iloc[0]["canonical_promotion_authorized"] != "FALSE":
        raise SpatialResourceError("downstream figures require the validated, unpromoted spatial candidate")
    if validation.iloc[0]["release_id"] != RESOURCE_RELEASE_ID:
        raise SpatialResourceError("upstream release ID does not match the spatial contract")


def save_pdf(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, bbox_inches="tight", dpi=400)
    plt.close(fig)


def display_count(value: object) -> str:
    if value is None or str(value).strip() in {"", "nan", "None"}:
        return "—"
    return str(value).strip()


def unit_phrase(count: object, unit: object, resolution: object, biological: bool) -> str:
    value = display_count(count)
    unit_text = str(unit)
    labels = {
        "donor": "donors",
        "section": "sections",
        "physical_array": "arrays",
        "AOI": "AOIs",
        "unknown_public_biological_unit": "biological units",
    }
    label = labels.get(unit_text, unit_text.replace("_", " "))
    if biological and str(resolution) != "resolved":
        return f"— unresolved {label}"
    return f"{value} {label}"


def join_state(value: str) -> str:
    if value == "resolved":
        return "resolved"
    if value == "not_applicable":
        return "not_applicable"
    return "unresolved"


def panel_fig1(registry: list[dict[str, object]], coverage: pd.DataFrame, panel_dir: Path, data_dir: Path) -> Path:
    registry_by_dataset = {str(row["dataset_id"]): row for row in registry}
    rows = []
    for dataset in COVERAGE_DATASETS:
        if dataset not in registry_by_dataset:
            raise SpatialResourceError(f"Figure 1 source missing from registry: {dataset}")
        source = registry_by_dataset[dataset]
        part = coverage[coverage["dataset_id"] == dataset]
        counts = part["coverage_status"].value_counts().to_dict()
        if len(part) != 117:
            raise SpatialResourceError(f"Figure 1 requires 117 coverage rows for {dataset}")
        rows.append({
            "integration_release_id": INTEGRATION_RELEASE_ID,
            "spatial_release_id": RESOURCE_RELEASE_ID,
            "dataset_id": dataset,
            "display_label": DATASET_LABELS[dataset],
            "dataset_gate": source["dataset_gate"],
            "biological_unit": source["biological_unit"],
            "biological_unit_resolution": source["biological_unit_resolution"],
            "donor_join": source["donor_join"],
            "phenotype_join": source["phenotype_join"],
            "histology_join": source["histology_join"],
            "spatial_coordinate_join": source["spatial_coordinate_join"],
            "n_biological": source["n_biological"],
            "technical_unit": source["technical_unit"],
            "n_technical": source["n_technical"],
            "n_programs": len(part),
            "n_observable": int(counts.get("observable", 0)),
            "n_partial": int(counts.get("partial", 0)),
            "n_untestable": int(counts.get("untestable", 0)),
        })
    table = pd.DataFrame(rows)
    table.to_csv(data_dir / "fig1_spatial_assay_observability.tsv", sep="\t", index=False)

    columns = ["gate", "unit", "donor", "phenotype", "histology", "coordinates", "coverage"]
    state_colors = {
        "resolved": CYAN,
        "source_dependent": CYAN_LIGHT,
        "metadata_pending": GRAY,
        "unresolved": GRAY,
        "not_applicable": WHITE,
        "dropped": DARK_GRAY,
        "near_complete": CYAN,
        "limited": CYAN_LIGHT,
        "untestable": GRAY,
    }
    state_symbols = {
        "resolved": "●",
        "source_dependent": "S",
        "metadata_pending": "?",
        "unresolved": "?",
        "not_applicable": "–",
        "dropped": "×",
    }
    fig, ax = plt.subplots(figsize=(5.5, 3.05))
    ax.set_xlim(0, 9.4)
    ax.set_ylim(0, len(table))
    ax.invert_yaxis()
    ax.set_xticks([index + 0.5 for index in range(len(columns))], columns)
    ax.xaxis.tick_top()
    ax.set_yticks([index + 0.5 for index in range(len(table))], table["display_label"])
    ax.tick_params(length=0, pad=2)
    for row_index, row in table.iterrows():
        gate = str(row["dataset_gate"])
        cell_states = [
            "resolved" if gate == "pass" else gate,
            join_state(str(row["biological_unit_resolution"])),
            join_state(str(row["donor_join"])),
            join_state(str(row["phenotype_join"])),
            join_state(str(row["histology_join"])),
            join_state(str(row["spatial_coordinate_join"])),
            "near_complete" if row["n_observable"] >= 100 else ("limited" if row["n_observable"] > 0 else "untestable"),
        ]
        for column_index, state in enumerate(cell_states):
            ax.add_patch(Rectangle((column_index, row_index), 1, 1, facecolor=state_colors[state], edgecolor=WHITE, linewidth=1.0))
            if column_index == len(columns) - 1:
                symbol = f"{row['n_observable']}/117"
            else:
                symbol = state_symbols.get(state, "●")
            text_color = WHITE if state in {"resolved", "near_complete", "dropped"} else INK
            ax.text(column_index + 0.5, row_index + 0.52, symbol, ha="center", va="center", color=text_color)
        biological = unit_phrase(row["n_biological"], row["biological_unit"], row["biological_unit_resolution"], True)
        technical = unit_phrase(row["n_technical"], row["technical_unit"], "resolved", False)
        ax.text(7.18, row_index + 0.52, f"{biological}  |  {technical}", ha="left", va="center", color=INK)
    ax.text(7.18, -0.22, "biological units  |  technical units", ha="left", va="bottom", color=INK)
    for spine in ax.spines.values():
        spine.set_visible(False)
    legend = [
        Patch(facecolor=CYAN, edgecolor=GRAY, linewidth=0.35, label="resolved / observable"),
        Patch(facecolor=CYAN_LIGHT, edgecolor=GRAY, linewidth=0.35, label="source-dependent / limited"),
        Patch(facecolor=GRAY, edgecolor=GRAY, linewidth=0.35, label="unresolved / untestable"),
        Patch(facecolor=WHITE, edgecolor=GRAY, linewidth=0.35, label="not applicable"),
        Patch(facecolor=DARK_GRAY, edgecolor=GRAY, linewidth=0.35, label="dropped"),
    ]
    ax.text(0, -0.085, "Observable counts are programs, not samples; S = source-dependent.", transform=ax.transAxes, ha="left", va="top", color=DARK_GRAY)
    ax.legend(handles=legend, frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.16), borderaxespad=0, columnspacing=0.9, handlelength=1.0)
    path = panel_dir / "fig1_spatial_assay_observability.pdf"
    save_pdf(fig, path)
    return path


def effect_row(effects: pd.DataFrame, dataset: str, uid: str) -> pd.Series:
    row = effects[(effects["dataset_id"] == dataset) & (effects["program_uid"] == uid)]
    if len(row) != 1:
        raise SpatialResourceError(f"expected one effect row for {dataset}/{uid}, found {len(row)}")
    return row.iloc[0]


def panel_fig5(effects: pd.DataFrame, hmsma: pd.DataFrame, panel_dir: Path, data_dir: Path) -> Path:
    ig_uid = "hotspot_hepatocytes_f05c535ae5bbc0b9"
    bi_uid = "hotspot_hepatocytes_48f39dd4d817a10e"
    ig_gse, ig_vu = effect_row(effects, "GSE192741", ig_uid), effect_row(effects, "Vu_et_al_2025", ig_uid)
    bi_gse, bi_vu = effect_row(effects, "GSE192741", bi_uid), effect_row(effects, "Vu_et_al_2025", bi_uid)
    medians = hmsma.groupby("program_uid", observed=True)["residual_moran_i"].median().to_dict()
    examples = [
        {
            "example_id": "igfbp7_program",
            "display_label": "Stromal ECM\n(IGFBP7)",
            "primary_evidence_state": "supported",
            "source_qualification": "Vu remains source-dependent",
            "observed_result": f"GSE I={float(ig_gse['estimate']):.3f}, q={float(ig_gse['qvalue']):.4f}\nVu I={float(ig_vu['estimate']):.3f}, q={float(ig_vu['qvalue']):.4f} (S)",
            "inference_unit": "GSE: 4 donors\nVu: 10 arrays;\ndonor key unresolved",
            "interpretation": "physical organization\nreproduced; not disease\ndirection",
            "next_experiment": "donor-resolved,\nhistology-registered\nspatial RNA/protein replication",
            "gse_estimate": float(ig_gse["estimate"]),
            "gse_qvalue": float(ig_gse["qvalue"]),
            "vu_estimate": float(ig_vu["estimate"]),
            "vu_qvalue": float(ig_vu["qvalue"]),
            "hmsma_median_moran_i": np.nan,
        },
        {
            "example_id": "bicc1_program",
            "display_label": "Ductular injury\n(BICC1)",
            "primary_evidence_state": "indeterminate",
            "source_qualification": "complete two-program family retained",
            "observed_result": f"GSE I={float(bi_gse['estimate']):.3f}, q={float(bi_gse['qvalue']):.4f}\nVu I={float(bi_vu['estimate']):.3f}, q={float(bi_vu['qvalue']):.4f}",
            "inference_unit": "GSE: 4 donors\nVu: 10 arrays;\ndonor key unresolved",
            "interpretation": "no adequate organization\ncall; indeterminate\nis not negative",
            "next_experiment": "higher-coverage\ntargeted assay in\nductular regions",
            "gse_estimate": float(bi_gse["estimate"]),
            "gse_qvalue": float(bi_gse["qvalue"]),
            "vu_estimate": float(bi_vu["estimate"]),
            "vu_qvalue": float(bi_vu["qvalue"]),
            "hmsma_median_moran_i": np.nan,
        },
        {
            "example_id": "hmsma_clinical_gate",
            "display_label": "HMSMA\ndisease inference",
            "primary_evidence_state": "untestable",
            "source_qualification": "metadata pending",
            "observed_result": f"35-array descriptive medians:\nIGFBP7 {medians[ig_uid]:.3f}; BICC1 {medians[bi_uid]:.3f}",
            "inference_unit": "35 arrays; donor key absent\n130,097 spots =\ntechnical coverage",
            "interpretation": "label-blind organization\nonly; no clinical or\npopulation claim",
            "next_experiment": "authoritative array-to-donor\nand clinical key; registered H&E\nfor localization",
            "gse_estimate": np.nan,
            "gse_qvalue": np.nan,
            "vu_estimate": np.nan,
            "vu_qvalue": np.nan,
            "hmsma_median_moran_i": medians[ig_uid],
        },
    ]
    table = pd.DataFrame(examples)
    table.insert(0, "spatial_release_id", RESOURCE_RELEASE_ID)
    table.insert(0, "integration_release_id", INTEGRATION_RELEASE_ID)
    table.to_csv(data_dir / "fig5_spatial_passport_examples.tsv", sep="\t", index=False)

    headers = ["example", "spatial state", "inference unit", "interpretation", "next discriminating experiment"]
    widths = [1.10, 1.55, 1.48, 1.42, 1.85]
    edges = np.cumsum([0] + widths)
    state_colors = {"supported": CYAN, "indeterminate": GRAY_LIGHT, "untestable": GRAY}
    fig, ax = plt.subplots(figsize=(7.5, 2.75))
    ax.set_xlim(0, edges[-1])
    ax.set_ylim(0, 4.1)
    ax.invert_yaxis()
    ax.axis("off")
    for index, header in enumerate(headers):
        ax.text(edges[index] + 0.06, 0.30, header, ha="left", va="center", color=INK)
    for row_index, row in table.iterrows():
        y = 0.62 + row_index * 1.08
        values = [row["display_label"], row["observed_result"], row["inference_unit"], row["interpretation"], row["next_experiment"]]
        for column_index, value in enumerate(values):
            face = WHITE
            text_color = INK
            if column_index == 1:
                face = state_colors[str(row["primary_evidence_state"])]
                text_color = WHITE if row["primary_evidence_state"] in {"supported", "untestable"} else INK
            ax.add_patch(Rectangle((edges[column_index], y), widths[column_index], 0.96, facecolor=face, edgecolor=WHITE if column_index == 1 else GRAY_LIGHT, linewidth=0.6))
            y_text = y + 0.58 if column_index == 1 else y + 0.48
            ax.text(edges[column_index] + 0.06, y_text, str(value), ha="left", va="center", color=text_color)
        ax.text(edges[1] + widths[1] - 0.06, y + 0.09, str(row["primary_evidence_state"]).replace("_", " "), ha="right", va="top", color=WHITE if row["primary_evidence_state"] in {"supported", "untestable"} else INK)
    ax.text(0, 3.98, "I = abundance-adjusted Moran's I; S = source-dependent. Spatial organization is not disease direction.", ha="left", va="bottom", color=DARK_GRAY)
    path = panel_dir / "fig5_spatial_passport_examples.pdf"
    save_pdf(fig, path)
    return path


def panel_coverage(coverage: pd.DataFrame, panel_dir: Path, data_dir: Path) -> Path:
    source_columns = [
        "release_id", "program_release_id", "dataset_id", "assay_id", "program_uid", "program_label",
        "cell_type", "n_program_genes", "n_genes_detected", "retained_l1_weight", "coverage_status",
        "testability_reason", "dataset_gate", "biological_unit_resolution", "source_dependence",
    ]
    source = coverage[coverage["dataset_id"].isin(COVERAGE_DATASETS)][source_columns].copy()
    source.insert(0, "integration_release_id", INTEGRATION_RELEASE_ID)
    source.to_csv(data_dir / "figS_spatial_program_coverage.tsv", sep="\t", index=False)
    reference = coverage[coverage["dataset_id"] == COVERAGE_DATASETS[0]].copy()
    program_order = reference["program_uid"].tolist()
    if len(program_order) != 117 or len(set(program_order)) != 117:
        raise SpatialResourceError("coverage heatmap requires the 117-program frozen registry order")
    labels = reference.set_index("program_uid")["program_label"].to_dict()
    cell_types = reference.set_index("program_uid")["cell_type"].to_dict()
    code = {"untestable": 0, "partial": 1, "observable": 2}
    matrix = np.empty((117, len(COVERAGE_DATASETS)), dtype=int)
    for column_index, dataset in enumerate(COVERAGE_DATASETS):
        part = coverage[coverage["dataset_id"] == dataset].set_index("program_uid")
        if set(part.index) != set(program_order):
            raise SpatialResourceError(f"coverage heatmap family mismatch for {dataset}")
        matrix[:, column_index] = [code[str(part.loc[uid, "coverage_status"])] for uid in program_order]
    cmap = ListedColormap([GRAY, CYAN_LIGHT, CYAN])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    fig, ax = plt.subplots(figsize=(4.75, 5.5))
    ax.imshow(matrix, aspect="auto", interpolation="nearest", cmap=cmap, norm=norm)
    column_labels = []
    for dataset in COVERAGE_DATASETS:
        n_observable = int((coverage[(coverage["dataset_id"] == dataset)]["coverage_status"] == "observable").sum())
        column_labels.append(f"{DATASET_LABELS[dataset]}\n{n_observable}/117 obs")
    ax.set_xticks(range(len(COVERAGE_DATASETS)), column_labels, rotation=45, ha="left", rotation_mode="anchor")
    ax.xaxis.tick_top()
    ax.tick_params(length=0, pad=2)
    positions_by_type: dict[str, list[int]] = {}
    for index, uid in enumerate(program_order):
        positions_by_type.setdefault(str(cell_types[uid]), []).append(index)
    ticks, tick_labels = [], []
    for cell_type, positions in positions_by_type.items():
        ticks.append(float(np.mean(positions)))
        tick_labels.append(cell_type.replace("_", " "))
        if positions[0] > 0:
            ax.axhline(positions[0] - 0.5, color=WHITE, linewidth=0.7)
    ax.set_yticks(ticks, tick_labels)
    for uid, short_label in CONFIRMATORY_UIDS.items():
        row_index = program_order.index(uid)
        ax.add_patch(Rectangle((-0.5, row_index - 0.5), len(COVERAGE_DATASETS), 1, fill=False, edgecolor=INK, linewidth=0.55))
        ax.text(len(COVERAGE_DATASETS) - 0.35, row_index, short_label, ha="left", va="center", color=INK, clip_on=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xlim(-0.5, len(COVERAGE_DATASETS) + 0.75)
    legend = [
        Patch(facecolor=CYAN, edgecolor=GRAY, linewidth=0.35, label="observable"),
        Patch(facecolor=CYAN_LIGHT, edgecolor=GRAY, linewidth=0.35, label="partial"),
        Patch(facecolor=GRAY, edgecolor=GRAY, linewidth=0.35, label="untestable"),
    ]
    ax.legend(handles=legend, frameon=False, ncol=3, loc="upper left", bbox_to_anchor=(0, -0.035), borderaxespad=0, columnspacing=0.8, handlelength=1.0)
    ax.text(-0.5, 119.8, "Frozen registry order; no spatial outcome used for ordering. Outlines mark the two predeclared programs.", ha="left", va="top", color=DARK_GRAY)
    path = panel_dir / "figS_spatial_program_coverage.pdf"
    save_pdf(fig, path)
    return path


def natural_array_order(value: str) -> tuple[str, int]:
    prefix, _, suffix = value.rpartition("_")
    return prefix, int(suffix) if suffix.isdigit() else 10**9


def panel_hmsma(hmsma: pd.DataFrame, panel_dir: Path, data_dir: Path) -> Path:
    source = hmsma.sort_values(["array_id", "program_uid"], key=lambda column: column.map(natural_array_order) if column.name == "array_id" else column).copy()
    authorization = source["inferential_pvalue_authorized"].astype(str).str.lower()
    if not authorization.isin({"false", "0"}).all():
        raise SpatialResourceError("HMSMA label-blind source unexpectedly authorizes inferential P values")
    source["inferential_pvalue_authorized"] = "FALSE"
    source.insert(0, "integration_release_id", INTEGRATION_RELEASE_ID)
    source.to_csv(data_dir / "figS_hmsma_label_blind_organization.tsv", sep="\t", index=False)
    arrays = sorted(hmsma["array_id"].unique(), key=natural_array_order)
    if len(arrays) != 35 or len(hmsma) != 70:
        raise SpatialResourceError("HMSMA supplement requires 35 arrays by two programs")
    program_order = list(CONFIRMATORY_UIDS)
    colors = [CYAN, CYAN_MID]
    fig, ax = plt.subplots(figsize=(3.15, 2.80))
    by_array = hmsma.pivot(index="array_id", columns="program_uid", values="residual_moran_i").loc[arrays]
    offsets = np.linspace(-0.13, 0.13, len(arrays))
    for array_index, array_id in enumerate(arrays):
        ax.plot([0 + offsets[array_index], 1 + offsets[array_index]], by_array.loc[array_id, program_order], color=GRAY_LIGHT, linewidth=0.35, zorder=1)
    for index, (uid, color) in enumerate(zip(program_order, colors)):
        values = by_array[uid].to_numpy(dtype=float)
        ax.scatter(index + offsets, values, s=8, color=color, edgecolor=WHITE, linewidth=0.25, zorder=2)
        median = float(np.median(values))
        ax.plot([index - 0.20, index + 0.20], [median, median], color=INK, linewidth=0.8, zorder=3)
        ax.text(index, max(values) + 0.025, f"median {median:.3f}", ha="center", va="bottom", color=INK)
    ax.axhline(0, color=GRAY, linewidth=0.5, linestyle="--", zorder=0)
    ax.set_xticks([0, 1], ["Stromal ECM\n(IGFBP7)", "Ductular injury\n(BICC1)"])
    ax.set_ylabel("abundance-adjusted Moran's I")
    ax.set_xlim(-0.38, 1.38)
    lower = min(-0.04, float(hmsma["residual_moran_i"].min()) - 0.02)
    upper = float(hmsma["residual_moran_i"].max()) + 0.08
    ax.set_ylim(lower, upper)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(0.35)
    ax.spines["bottom"].set_linewidth(0.35)
    ax.tick_params(width=0.35, length=2)
    ax.text(-0.36, lower - (upper - lower) * 0.14, "35 physical arrays; donor identity unresolved; descriptive only; no population P value.", ha="left", va="top", color=DARK_GRAY)
    path = panel_dir / "figS_hmsma_label_blind_organization.pdf"
    save_pdf(fig, path)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--upstream-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    project = args.project_root.resolve()
    upstream = args.upstream_root.resolve()
    output = args.output_root.resolve()
    if output.exists():
        raise SpatialResourceError(f"immutable integration candidate already exists: {output}")
    require_ready(upstream)
    registry_payload = json.loads((upstream / "registry/spatial_dataset_registry.json").read_text(encoding="utf-8"))
    registry = registry_payload["datasets"]
    coverage = pd.read_parquet(upstream / "coverage/spatial_program_coverage.parquet")
    effects = pd.read_parquet(upstream / "effects/spatial_program_effects.parquet")
    hmsma = pd.read_csv(upstream / "hmsma_label_blind/per_array_program_organization.tsv", sep="\t")
    temp = Path(tempfile.mkdtemp(prefix=f".{output.name}.", dir=output.parent))
    try:
        panel_dir = temp / "panels"
        data_dir = temp / "data"
        panel_dir.mkdir(parents=True)
        data_dir.mkdir(parents=True)
        panels = [
            panel_fig1(registry, coverage, panel_dir, data_dir),
            panel_fig5(effects, hmsma, panel_dir, data_dir),
            panel_coverage(coverage, panel_dir, data_dir),
            panel_hmsma(hmsma, panel_dir, data_dir),
        ]
        write_tsv(
            temp / "BUILD_COMPLETE",
            ("integration_release_id", "spatial_release_id", "status", "n_panels", "upstream_validation_sha256", "canonical_figure_written"),
            [{
                "integration_release_id": INTEGRATION_RELEASE_ID,
                "spatial_release_id": RESOURCE_RELEASE_ID,
                "status": "rendered_candidate_awaiting_independent_validation",
                "n_panels": len(panels),
                "upstream_validation_sha256": sha256_file(upstream / "validation/READY"),
                "canonical_figure_written": "FALSE",
            }],
        )
        os.replace(temp, output)
    except Exception:
        shutil.rmtree(temp, ignore_errors=True)
        raise
    print(f"WROTE {len(panels)} candidate-only downstream spatial panels to {output}")


if __name__ == "__main__":
    main()
