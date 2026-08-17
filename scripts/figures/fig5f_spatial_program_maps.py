#!/usr/bin/env python3
"""Figure 5F — physical maps of two representative replicated spatial programs.

KEY MESSAGE: Frozen disease-associated programs retain abundance-adjusted spatial organization in two independent tissue cohorts.

The two displayed programs are frozen Figure 3 modules selected as biologically
distinct examples after they passed the matched-null spatial robustness gate in
both GSE192741 and Vu. They are illustrative rather than an independent test.
For each cohort, one section is chosen without maximizing either program: the
section whose mean within-program Moran rank is closest to the cohort median.
Scores reproduce the primary spatial projection and are residualized for
matching-lineage abundance, library size, and detected-gene count, then
standardized within cohort for display.
"""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import anndata as ad
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
import numpy as np
import pandas as pd


BASE = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
ROOT = BASE / "Analysis/Multimodal_Program_Projection"
PANELS = BASE / "figures/main/fig5_molecular_context/panels"
PROGRAMS = {
    "hepatocytes::14": "Hep-14  Amino-acid metabolism",
    "fibroblasts::6": "Fib-6  Biliary-like",
}
DATASETS = {
    "GSE192741": BASE
    / "Analysis/Spatial/results/cell2location/spatial_model/spatial_deconvolved.h5ad",
    "Vu_et_al_2025": BASE
    / "Analysis/Spatial/results/cell2location/spatial_model_vu/spatial_deconvolved_vu.h5ad",
}
DATASET_LABEL = {"GSE192741": "GSE192741", "Vu_et_al_2025": "Vu et al."}
CT_ABUNDANCE = {
    "hepatocytes": "Hepatocytes",
    "fibroblasts": "Fibroblasts",
}


def load_projection_helpers():
    path = ROOT / "scripts/04_spatial_projection.py"
    spec = importlib.util.spec_from_file_location("fig5_spatial_projection", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import spatial projection helpers from {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def choose_median_sections(section_results: pd.DataFrame) -> pd.DataFrame:
    d = section_results[section_results["program_id"].isin(PROGRAMS)].copy()
    if d.empty:
        raise RuntimeError("No section-level spatial results for the requested programs")
    d["within_program_rank"] = d.groupby(["dataset", "program_id"])[
        "residual_moran_i"
    ].rank(method="average", pct=True)
    ranked = (
        d.groupby(["dataset", "sample_id"], as_index=False)
        .agg(
            mean_program_rank=("within_program_rank", "mean"),
            mean_residual_moran=("residual_moran_i", "mean"),
            individual=("individual", "first"),
            condition=("condition", "first"),
            n_programs=("program_id", "nunique"),
        )
    )
    ranked = ranked[ranked["n_programs"] == len(PROGRAMS)].copy()
    ranked["distance_to_median"] = (ranked["mean_program_rank"] - 0.5).abs()
    chosen = (
        ranked.sort_values(["dataset", "distance_to_median", "sample_id"])
        .groupby("dataset", as_index=False)
        .first()
    )
    if set(chosen["dataset"]) != set(DATASETS):
        raise RuntimeError("Failed to choose one representative section per cohort")
    chosen["selection_rule"] = (
        "same section for both programs; mean within-program residual-Moran rank closest to cohort median"
    )
    return chosen


def score_dataset(dataset: str, path: Path, selected_sample: str, helpers, registry, membership):
    adata = ad.read_h5ad(path)
    adata.var_names = adata.var_names.astype(str)
    obs = adata.obs.copy()
    obs["sample_id"] = obs["sample_id"].astype(str)
    counts = helpers.dense_counts(adata)
    library = pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(float)
    norm = helpers.log_normalize(counts, library)
    detection = np.asarray((counts > 0).mean(axis=0)).ravel()
    gene_to_index = {g: i for i, g in enumerate(adata.var_names)}

    factor_names = [
        str(x).replace("means_per_cluster_mu_fg_", "")
        for x in adata.uns["mod"]["factor_names"]
    ]
    abundance = np.asarray(adata.obsm["q05_cell_abundance_w_sf"])
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    keep_spots = obs["sample_id"].to_numpy() == selected_sample
    if keep_spots.sum() < 50:
        raise RuntimeError(f"Selected section {dataset}/{selected_sample} has too few spots")

    output = {}
    audit = []
    for program_id in PROGRAMS:
        row = registry.loc[registry["program_id"] == program_id].iloc[0]
        genes = membership.loc[membership["program_id"] == program_id].copy()
        genes = genes.groupby("gene_symbol", as_index=False)["original_l1_weight"].sum()
        genes["present"] = genes["gene_symbol"].isin(gene_to_index)
        genes["detected"] = genes["gene_symbol"].map(
            lambda g: detection[gene_to_index[g]] >= 0.01 if g in gene_to_index else False
        )
        measured = genes[genes["present"] & genes["detected"]].copy()
        retained = float(measured["original_l1_weight"].sum())
        if len(measured) < 8 or retained < 0.20:
            raise RuntimeError(f"Program {program_id} is not spatially testable in {dataset}")
        weights = measured["original_l1_weight"].to_numpy(float) / retained
        indices = np.array([gene_to_index[g] for g in measured["gene_symbol"]], dtype=int)
        z = helpers.extract_z(norm, indices)
        score = z @ weights
        lineage = str(row["cell_type"])
        lineage_abundance = abundance[:, factor_names.index(CT_ABUNDANCE[lineage])]
        design = helpers.design_matrix(obs, lineage_abundance)
        residual = helpers.residualize(score, design)
        residual = (residual - np.nanmean(residual)) / np.nanstd(residual, ddof=1)
        output[program_id] = pd.DataFrame(
            {
                "x": coords[keep_spots, 0],
                "y": coords[keep_spots, 1],
                "score_z": residual[keep_spots],
            }
        )
        audit.append(
            {
                "dataset": dataset,
                "sample_id": selected_sample,
                "program_id": program_id,
                "n_spots": int(keep_spots.sum()),
                "n_measured": int(len(measured)),
                "retained_l1_weight": retained,
                "score_definition": (
                    "weighted gene-z score residualized for lineage abundance, log library size, and detected genes; z within cohort"
                ),
            }
        )
    del adata, counts, norm
    return output, audit


def main():
    helpers = load_projection_helpers()
    registry = pd.read_csv(ROOT / "results/frozen_programs.tsv", sep="\t")
    membership = pd.read_csv(ROOT / "results/frozen_program_membership.tsv", sep="\t")
    membership = membership[
        membership["mapped_symbol"].fillna(False)
        & membership["gene_symbol"].notna()
    ].copy()
    spatial = pd.read_csv(ROOT / "results/spatial/spatial_program_results.tsv", sep="\t")
    robust = spatial[
        spatial["program_id"].isin(PROGRAMS)
        & spatial["dataset"].isin(DATASETS)
    ]
    if len(robust) != len(PROGRAMS) * len(DATASETS) or not robust["robust"].all():
        raise RuntimeError("Representative map programs must be robust in both spatial cohorts")

    sections = pd.read_csv(ROOT / "results/spatial/spatial_section_results.tsv", sep="\t")
    chosen = choose_median_sections(sections)
    maps = {}
    audit_rows = []
    for dataset, path in DATASETS.items():
        sample = str(chosen.loc[chosen["dataset"] == dataset, "sample_id"].iloc[0])
        maps[dataset], audit = score_dataset(
            dataset, path, sample, helpers, registry, membership
        )
        audit_rows.extend(audit)

    PANELS.mkdir(parents=True, exist_ok=True)
    data_dir = PANELS / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    audit = pd.DataFrame(audit_rows).merge(
        chosen[
            [
                "dataset",
                "sample_id",
                "individual",
                "condition",
                "mean_program_rank",
                "mean_residual_moran",
                "selection_rule",
            ]
        ],
        on=["dataset", "sample_id"],
        how="left",
    )
    audit.to_csv(data_dir / "fig5f_spatial_program_maps_selection.tsv", sep="\t", index=False)

    cmap = LinearSegmentedColormap.from_list(
        "program_score", ["#2166AC", "#F7F7F7", "#B2182B"]
    )
    norm = Normalize(vmin=-2, vmax=2, clip=True)
    fig, axes = plt.subplots(2, 2, figsize=(4.05, 3.05))
    for r, dataset in enumerate(DATASETS):
        sample = str(chosen.loc[chosen["dataset"] == dataset, "sample_id"].iloc[0])
        for c, (program_id, title) in enumerate(PROGRAMS.items()):
            ax = axes[r, c]
            d = maps[dataset][program_id]
            ax.scatter(
                d["x"],
                -d["y"],
                c=d["score_z"],
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
            if r == 0:
                ax.set_title(title, fontsize=6, fontfamily="Helvetica", pad=2.5)
            if c == 0:
                ax.text(
                    -0.035,
                    0.5,
                    DATASET_LABEL[dataset],
                    transform=ax.transAxes,
                    ha="right",
                    va="center",
                    rotation=90,
                    fontsize=6,
                    fontfamily="Helvetica",
                )
            ax.text(
                0.01,
                0.01,
                sample,
                transform=ax.transAxes,
                ha="left",
                va="bottom",
                    fontsize=6,
                color="#555555",
                fontfamily="Helvetica",
            )

    fig.subplots_adjust(left=0.07, right=0.99, top=0.92, bottom=0.14, wspace=0.04, hspace=0.08)
    cax = fig.add_axes([0.34, 0.055, 0.34, 0.018])
    cb = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, orientation="horizontal")
    cb.set_ticks([-2, 0, 2])
    cb.ax.tick_params(labelsize=6, width=0.35, length=1.5, pad=1)
    cb.outline.set_linewidth(0.35)
    cb.set_label("abundance-adjusted program score (z)", fontsize=6, labelpad=1)
    out = PANELS / "fig5f_spatial_program_maps.pdf"
    fig.savefig(out, bbox_inches="tight", dpi=400)
    plt.close(fig)
    print(f"[fig5f maps] saved: {out}")
    print(
        "CAPTION (Fig. 5F): Abundance-adjusted spatial scores for two representative frozen programs "
        "selected for visualization after passing the matched-gene spatial robustness gate in both cohorts. "
        "These maps are illustrative rather than an independent test. The same section is shown "
        "for both programs within each cohort and is selected as the section whose mean within-program "
        "Moran rank is closest to the cohort median, avoiding selection of an extreme map. Scores are "
        "residualized for matching-lineage abundance, library size, and detected genes and standardized "
        "within cohort for display; inferential matched-null results are summarized in Panel 4E."
    )


if __name__ == "__main__":
    main()
