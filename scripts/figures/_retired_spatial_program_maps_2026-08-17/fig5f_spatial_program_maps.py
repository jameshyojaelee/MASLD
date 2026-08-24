#!/usr/bin/env python3
"""Figure 5F: physical maps of two prespecified frozen programs.

KEY MESSAGE: Two frozen disease programs occupy different abundance-adjusted
spatial patterns, while matched-null support belongs to the program score and
not to any single member gene.

The map values and outcome-blind section selection come from the sealed v2
spatial-map source. One section or array per dataset was selected solely by
eligible spot count: closest to the dataset median, with lexical ID tie-break.
The same physical unit is shown for both programs. IGFBP7 and BICC1 are retained
only as registry provenance and are not used as program labels or inference
units.
"""

from __future__ import annotations

import os
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
        "axes.titleweight": "normal",
        "axes.labelweight": "normal",
    }
)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
import pandas as pd


BASE = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
CANDIDATE_ROOT = os.environ.get("FIGURE_CANDIDATE_ROOT", "")
PANELS = (
    Path(CANDIDATE_ROOT) / "figure5/panels"
    if CANDIDATE_ROOT
    else BASE / "figures/main/fig5_molecular_context/panels"
)
MAP_ROOT = (
    BASE
    / "Analysis/Multimodal_Program_Projection/candidates/"
    "program-context-v2-candidate-2026-08-07/spatial_context/map_source_v2"
)
MAP_SOURCE = MAP_ROOT / "per_spot_program_map.tsv.gz"
SELECTION_SOURCE = MAP_ROOT / "map_selection.tsv"
REGISTRY = BASE / "Analysis/Multimodal_Program_Projection/results/frozen_programs.tsv"
MEMBERSHIP = (
    BASE
    / "Analysis/Multimodal_Program_Projection/results/frozen_program_membership.tsv"
)

PROGRAMS = {
    "hepatocytes::8": "Stromal ECM program",
    "hepatocytes::20": "Ductular injury program",
}
PROGRAM_TITLE = {
    "hepatocytes::8": "Stromal ECM",
    "hepatocytes::20": "Ductular injury",
}
NAMED_GENE = {"hepatocytes::8": "IGFBP7", "hepatocytes::20": "BICC1"}
DATASETS = ["GSE192741", "Vu_et_al_2025"]
DATASET_LABEL = {
    "GSE192741": "GSE192741",
    "Vu_et_al_2025": "Vu et al.",
}


def load_and_validate_sources():
    maps = pd.read_csv(MAP_SOURCE, sep="\t", compression="gzip")
    selection = pd.read_csv(SELECTION_SOURCE, sep="\t")
    maps = maps[
        maps["legacy_program_id"].isin(PROGRAMS)
        & maps["dataset"].isin(DATASETS)
        & maps["graph_eligible"].astype(bool)
    ].copy()
    if set(maps["legacy_program_id"]) != set(PROGRAMS):
        raise RuntimeError("Sealed map source does not contain both frozen programs")
    if set(maps["dataset"]) != set(DATASETS):
        raise RuntimeError("Sealed map source does not contain both spatial datasets")
    if not (
        selection["program_outcomes_used_for_selection"]
        .astype(str)
        .str.upper()
        .eq("FALSE")
        .all()
    ):
        raise RuntimeError("Map section selection must be outcome-blind")
    if not (
        selection["selection_rule"]
        == "closest_to_dataset_median_source_spot_count_then_lexical_reporting_unit_id"
    ).all():
        raise RuntimeError("Unexpected map section-selection rule")
    selected_ids = dict(zip(selection["dataset"], selection["reporting_unit_id"]))
    for dataset in DATASETS:
        observed = set(maps.loc[maps["dataset"] == dataset, "reporting_unit_id"])
        if observed != {selected_ids[dataset]}:
            raise RuntimeError(f"Map source does not match sealed selection for {dataset}")
    return maps, selection


def build_audit(maps, selection):
    registry = pd.read_csv(REGISTRY, sep="\t").set_index("program_id")
    membership = pd.read_csv(MEMBERSHIP, sep="\t")
    rows = []
    for dataset in DATASETS:
        selected = selection.loc[selection["dataset"] == dataset].iloc[0]
        for program_id, display_name in PROGRAMS.items():
            part = maps[
                (maps["dataset"] == dataset)
                & (maps["legacy_program_id"] == program_id)
            ]
            named_gene = NAMED_GENE[program_id]
            gene_row = membership[
                (membership["program_id"] == program_id)
                & (membership["gene_symbol"] == named_gene)
            ]
            if len(gene_row) != 1:
                raise RuntimeError(f"Missing unique {named_gene} membership for {program_id}")
            rows.append(
                {
                    "dataset": dataset,
                    "sample_id": str(selected["reporting_unit_id"]),
                    "program_id": program_id,
                    "display_name": display_name,
                    "source_program_name": registry.loc[program_id, "program_name"],
                    "named_gene_provenance": named_gene,
                    "named_gene_l1_weight": float(gene_row.iloc[0]["original_l1_weight"]),
                    "program_n_genes": int(registry.loc[program_id, "n_symbol_mapped"]),
                    "n_spots": int(len(part)),
                    "n_measured": int(part["n_genes_measured"].iloc[0]),
                    "retained_l1_weight": float(part["retained_l1_weight"].iloc[0]),
                    "selection_rule": selected["selection_rule"],
                    "program_outcomes_used_for_selection": False,
                    "same_unit_for_all_programs": True,
                    "biological_unit_resolution": selected[
                        "biological_unit_resolution"
                    ],
                    "inference_unit": "weighted frozen program score; not named-gene expression",
                }
            )
    return pd.DataFrame(rows)


def main():
    maps, selection = load_and_validate_sources()
    audit = build_audit(maps, selection)
    PANELS.mkdir(parents=True, exist_ok=True)
    data_dir = PANELS / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    audit.to_csv(
        data_dir / "fig5f_spatial_program_maps_selection.tsv",
        sep="\t",
        index=False,
    )

    cmap = LinearSegmentedColormap.from_list(
        "program_score", ["#2166AC", "#F7F7F7", "#B2182B"]
    )
    norm = Normalize(vmin=-2, vmax=2, clip=True)
    fig, axes = plt.subplots(2, 2, figsize=(4.05, 3.05))
    for row_index, dataset in enumerate(DATASETS):
        for column_index, program_id in enumerate(PROGRAMS):
            ax = axes[row_index, column_index]
            part = maps[
                (maps["dataset"] == dataset)
                & (maps["legacy_program_id"] == program_id)
            ]
            ax.scatter(
                part["x"],
                -part["y"],
                c=part["residual_program_score_z"],
                cmap=cmap,
                norm=norm,
                s=2.0,
                linewidths=0,
                rasterized=True,
            )
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            for spine in ax.spines.values():
                spine.set_visible(False)
            if row_index == 0:
                ax.set_title(PROGRAM_TITLE[program_id], pad=2.5)
            if column_index == 0:
                ax.text(
                    -0.035,
                    0.5,
                    DATASET_LABEL[dataset],
                    transform=ax.transAxes,
                    ha="right",
                    va="center",
                    rotation=90,
                )
            sample = str(part["reporting_unit_id"].iloc[0])
            ax.text(
                0.01,
                0.01,
                sample,
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                color="#555555",
            )

    fig.subplots_adjust(
        left=0.12,
        right=0.99,
        top=0.92,
        bottom=0.14,
        wspace=0.04,
        hspace=0.08,
    )
    cax = fig.add_axes([0.34, 0.055, 0.34, 0.018])
    colorbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=cax,
        orientation="horizontal",
    )
    colorbar.set_ticks([-2, 0, 2])
    colorbar.ax.tick_params(width=0.35, length=1.5, pad=1)
    colorbar.outline.set_linewidth(0.35)
    colorbar.set_label("abundance-adjusted program score (z)", labelpad=1)
    output = PANELS / "fig5f_spatial_program_maps.pdf"
    fig.savefig(output, bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"[fig5f maps] saved: {output}")
    print(
        "CAPTION (Fig. 5F): Abundance-adjusted spatial scores for the prespecified "
        "Stromal ECM and Ductular injury programs. One outcome-blind section or "
        "array per source is shown for both programs. Scores are weighted frozen-program "
        "scores residualized for hepatocyte abundance, library size, and detected genes, "
        "then standardized within source for display. Matched-null inference is summarized "
        "in Panel 5E and applies to the program, not to IGFBP7 or BICC1 individually."
    )


if __name__ == "__main__":
    main()
