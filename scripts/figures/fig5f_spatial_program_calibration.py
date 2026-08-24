#!/usr/bin/env python3
"""Figure 5F: matched-gene calibration and registered tissue context.

KEY MESSAGE: After adjustment for all inferred cell abundances, the Stromal
ECM program exceeds matched-gene spatial background in both sources, whereas
the Ductular injury program does not.

Artwork uses author citations. Dataset accessions remain only in the source
sidecar. Moran excess is the all-cell-adjusted observed residual Moran's I
minus the matched-gene null mean. BH q-values are inherited from the sealed
9,999-draw matched-null analysis and are not recomputed during rendering.

The calibration is paired with one outcome-blind representative Guilliams et
al. section. Its registered H&E overlay shows the same all-16-cell-type-
adjusted Stromal ECM score used for the calibration analysis. The map is
illustrative; it does not assign the signal to a histologic structure.
"""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
import sys

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
from PIL import Image


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
SOURCE = (
    BASE
    / "Analysis/Multimodal_Program_Projection/candidates/"
    "spatial-full-composition-sensitivity-candidate-2026-08-11-r2/"
    "full_composition_program_results.tsv"
)
FULL_COMPOSITION_ROOT = SOURCE.parent
FULL_COMPOSITION_DESIGN_AUDIT = FULL_COMPOSITION_ROOT / "full_composition_design_audit.tsv"
FULL_COMPOSITION_SOURCE_MANIFEST = FULL_COMPOSITION_ROOT / "source_manifest.tsv"
MAP_SELECTION_SOURCE = (
    BASE
    / "Analysis/Multimodal_Program_Projection/candidates/"
    "program-context-v2-candidate-2026-08-07/spatial_context/map_source_v2/"
    "map_selection.tsv"
)
SPATIAL_V2_SCRIPTS = (
    BASE / "Analysis/Multimodal_Program_Projection/scripts/spatial_context_v2"
)
SPACERANGER_ROOT = BASE / "Analysis/Spatial/results/spaceranger/GSE192741"
OVERLAY_PROGRAM_UID = "hotspot_hepatocytes_f05c535ae5bbc0b9"
EXPECTED_FACTOR_COUNT = 16

PROGRAMS = {
    "hotspot_hepatocytes_f05c535ae5bbc0b9": {
        "program_id": "hepatocytes::8",
        "display_name": "Stromal ECM",
    },
    "hotspot_hepatocytes_48f39dd4d817a10e": {
        "program_id": "hepatocytes::20",
        "display_name": "Ductular injury",
    },
}
SOURCES = {
    "GSE192741": {
        "display_source": "Guilliams et al.",
        "citation": "Guilliams et al. 2022",
        "order": 0,
    },
    "Vu_et_al_2025": {
        "display_source": "Vu et al.",
        "citation": "Vu et al. 2025",
        "order": 1,
    },
}
SUPPORTED = "#007C91"
INDETERMINATE = "#8C8C8C"
ZERO = "#BDBDBD"
INK = "#222222"
PROGRAM_CMAP = LinearSegmentedColormap.from_list(
    "adjusted_program", ["#1565C0", "#F7F7F7", "#C9265E"]
)
PROGRAM_NORM = Normalize(vmin=-2.0, vmax=2.0, clip=True)


def q_label(value: float) -> str:
    if value < 0.001:
        return "q<.001"
    return f"q={value:.3f}".replace("0.", ".")


def load_source() -> pd.DataFrame:
    data = pd.read_csv(SOURCE, sep="\t")
    data = data[
        data["program_id"].isin(PROGRAMS)
        & data["dataset"].isin(SOURCES)
    ].copy()
    if len(data) != 4:
        raise RuntimeError(f"Expected four program-by-source rows, found {len(data)}")
    if data.duplicated(["program_id", "dataset"]).any():
        raise RuntimeError("Duplicate program-by-source rows in spatial source")
    if not data["testable"].astype(bool).all() or not (data["n_null"] == 9999).all():
        raise RuntimeError("Spatial source is not the sealed testable 9,999-null family")

    data["program_uid"] = data["program_id"]
    data["legacy_program_id"] = data["program_uid"].map(
        lambda value: PROGRAMS[value]["program_id"]
    )
    data["display_name"] = data["program_uid"].map(
        lambda value: PROGRAMS[value]["display_name"]
    )
    data["source_dataset_id"] = data["dataset"]
    data["display_source"] = data["dataset"].map(
        lambda value: SOURCES[value]["display_source"]
    )
    data["source_citation"] = data["dataset"].map(
        lambda value: SOURCES[value]["citation"]
    )
    data["source_order"] = data["dataset"].map(
        lambda value: SOURCES[value]["order"]
    )
    data["moran_excess"] = data["residual_moran_i"] - data["residual_null_mean"]
    data["bh_q"] = data["residual_padj"]
    data["supported_within_source"] = data["within_source_support"].astype(bool)
    data["display_state"] = np.where(
        data["supported_within_source"], "supported", "indeterminate"
    )
    data["inference_note"] = (
        "all-cell-adjusted weighted program Moran excess over matched-gene null; "
        "not named-gene expression"
    )
    return data.sort_values(["program_uid", "source_order"])


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(8 * 1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def standardize_design(columns: list[np.ndarray]) -> np.ndarray:
    """Reproduce the sealed full-composition nuisance-design convention."""
    design = np.column_stack(columns).astype(float, copy=False)
    for column in range(1, design.shape[1]):
        mean = np.nanmean(design[:, column])
        sd = np.nanstd(design[:, column], ddof=1)
        design[:, column] = (
            (design[:, column] - mean) / sd
            if np.isfinite(sd) and sd > 0
            else 0.0
        )
    design[~np.isfinite(design)] = 0.0
    return design


def graph_indices(graphs: dict, sample: str | None = None) -> np.ndarray:
    selected = graphs if sample is None else {sample: graphs[sample]}
    return np.asarray(
        sorted(
            {
                int(index)
                for components in selected.values()
                for component in components
                for index in component[0]
            }
        ),
        dtype=int,
    )


def load_histology_overlay() -> tuple[pd.DataFrame, Image.Image, dict[str, object]]:
    """Recompute the exact all-cell-adjusted display score for one sealed section."""
    sys.path.insert(0, str(SPATIAL_V2_SCRIPTS))
    try:
        from visium_rerun_lib import (
            build_paths,
            load_legacy_engine,
            prepare_engine_inputs,
            verify_hotspot_ready,
        )
    finally:
        sys.path.pop(0)

    selection = pd.read_csv(MAP_SELECTION_SOURCE, sep="\t")
    selection = selection[selection["dataset"] == "GSE192741"].copy()
    if len(selection) != 1:
        raise RuntimeError("Expected one sealed Guilliams et al. map selection")
    selected = selection.iloc[0]
    expected_rule = (
        "closest_to_dataset_median_source_spot_count_then_lexical_reporting_unit_id"
    )
    if (
        selected["reporting_unit_id"] != "JBO019"
        or selected["selection_rule"] != expected_rule
        or str(selected["program_outcomes_used_for_selection"]).upper() != "FALSE"
        or str(selected["same_unit_for_all_programs"]).upper() != "TRUE"
    ):
        raise RuntimeError("Guilliams et al. representative-section contract drifted")
    sample_id = str(selected["reporting_unit_id"])

    paths = build_paths(BASE)
    verify_hotspot_ready(paths)
    registry, membership, _ = prepare_engine_inputs(paths, "v2")
    engine = load_legacy_engine(paths)
    program = registry[registry["program_id"] == OVERLAY_PROGRAM_UID]
    if len(program) != 1:
        raise RuntimeError("Stromal ECM program is not unique in the frozen display registry")

    design_audit = pd.read_csv(FULL_COMPOSITION_DESIGN_AUDIT, sep="\t")
    design_audit = design_audit[design_audit["dataset"] == "GSE192741"]
    if len(design_audit) != 1:
        raise RuntimeError("Missing unique Guilliams et al. full-composition design audit")
    factor_order = str(design_audit.iloc[0]["factor_order"]).split("|")
    if len(factor_order) != EXPECTED_FACTOR_COUNT or len(set(factor_order)) != EXPECTED_FACTOR_COUNT:
        raise RuntimeError("Full-composition factor order is not the sealed 16-factor design")

    source_manifest = pd.read_csv(FULL_COMPOSITION_SOURCE_MANIFEST, sep="\t")
    input_manifest = source_manifest[
        source_manifest["source_role"] == "GSE192741_input_h5ad"
    ]
    if len(input_manifest) != 1:
        raise RuntimeError("Missing unique Guilliams et al. h5ad source manifest row")
    input_path = BASE / str(input_manifest.iloc[0]["relative_path"])
    if input_path.resolve() != paths.datasets["GSE192741"].resolve():
        raise RuntimeError("Full-composition and figure input h5ad paths disagree")

    import anndata as ad

    adata = ad.read_h5ad(input_path)
    try:
        adata.var_names = adata.var_names.astype(str)
        obs = adata.obs.copy()
        obs["sample_id"] = obs["sample_id"].astype(str)
        obs["individual"] = obs["individual"].astype(str)
        if set(obs.loc[obs["sample_id"] == sample_id, "individual"]) != {"H37"}:
            raise RuntimeError("Selected section no longer maps uniquely to donor H37")
        if set(obs.loc[obs["sample_id"] == sample_id, "condition"].astype(str)) != {
            "Steatotic"
        }:
            raise RuntimeError("Selected section no longer has the sealed Steatotic label")

        raw_factor_names = [
            str(value).replace("means_per_cluster_mu_fg_", "")
            for value in adata.uns["mod"]["factor_names"]
        ]
        abundance = np.asarray(adata.obsm["q05_cell_abundance_w_sf"], dtype=float)
        if set(raw_factor_names) != set(factor_order):
            raise RuntimeError("Cell2location factors no longer match the sealed design")
        abundance = abundance[:, [raw_factor_names.index(name) for name in factor_order]]
        total = pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(float)
        detected_genes = pd.to_numeric(
            obs["n_genes_by_counts"], errors="coerce"
        ).to_numpy(float)
        design_columns = [np.ones(len(obs), dtype=float)]
        design_columns.extend(abundance[:, column] for column in range(abundance.shape[1]))
        design_columns.extend([np.log1p(total), np.log1p(detected_genes)])
        design = standardize_design(design_columns)
        if design.shape[1] != int(design_audit.iloc[0]["n_design_columns"]):
            raise RuntimeError("Reconstructed full-composition design has the wrong width")

        counts = engine.dense_counts(adata)
        norm = engine.log_normalize(counts, total)
        detection = np.asarray((counts > 0).mean(axis=0)).ravel()
        gene_to_index = {gene: index for index, gene in enumerate(adata.var_names)}
        measured = membership[membership["program_id"] == OVERLAY_PROGRAM_UID].copy()
        measured = measured.groupby("gene_symbol", as_index=False)[
            "original_l1_weight"
        ].sum()
        measured["present"] = measured["gene_symbol"].isin(gene_to_index)
        measured["detected"] = measured["gene_symbol"].map(
            lambda gene: (
                detection[gene_to_index[gene]] >= 0.01
                if gene in gene_to_index
                else False
            )
        )
        measured = measured[measured["present"] & measured["detected"]].copy()
        retained = float(measured["original_l1_weight"].sum())
        if len(measured) < 8 or retained < 0.20:
            raise RuntimeError("Stromal ECM no longer passes the native Visium testability gate")
        weights = measured["original_l1_weight"].to_numpy(float) / retained
        indices = np.asarray(
            [gene_to_index[gene] for gene in measured["gene_symbol"]], dtype=int
        )
        score = engine.extract_z(norm, indices) @ weights
        residual = np.asarray(engine.residualize(score, design), dtype=float)

        coords = np.asarray(adata.obsm["spatial"], dtype=float)
        graphs, _ = engine.section_graphs(obs, coords, k=6)
        if sample_id not in graphs:
            raise RuntimeError("Selected section lacks a valid sealed tissue-island graph")
        observed = float(
            engine.donor_collapse(
                engine.moran_columns(residual, graphs), obs, "GSE192741"
            )[0]
        )
        sealed_result = pd.read_csv(SOURCE, sep="\t")
        sealed_result = sealed_result[
            (sealed_result["dataset"] == "GSE192741")
            & (sealed_result["program_id"] == OVERLAY_PROGRAM_UID)
        ]
        if len(sealed_result) != 1 or not np.isclose(
            observed,
            float(sealed_result.iloc[0]["residual_moran_i"]),
            rtol=0,
            atol=1e-10,
        ):
            raise RuntimeError("Overlay score does not reproduce the sealed residual Moran's I")

        eligible = graph_indices(graphs)
        selected_indices = graph_indices(graphs, sample_id)
        center = float(np.mean(residual[eligible]))
        spread = float(np.std(residual[eligible], ddof=1))
        if not np.isfinite(spread) or spread <= 0:
            raise RuntimeError("All-cell-adjusted display score is degenerate")
        residual_z = (residual - center) / spread

        spatial_dir = SPACERANGER_ROOT / sample_id / "outs/spatial"
        image_path = spatial_dir / "tissue_hires_image.png"
        scale_path = spatial_dir / "scalefactors_json.json"
        with scale_path.open("r", encoding="utf-8") as handle:
            scale_factors = json.load(handle)
        hires_scale = float(scale_factors["tissue_hires_scalef"])
        image = Image.open(image_path).convert("RGB")
        hires_x = coords[selected_indices, 0] * hires_scale
        hires_y = coords[selected_indices, 1] * hires_scale
        if (
            np.min(hires_x) < 0
            or np.max(hires_x) >= image.width
            or np.min(hires_y) < 0
            or np.max(hires_y) >= image.height
        ):
            raise RuntimeError("Registered spot coordinates fall outside the H&E image")

        overlay = pd.DataFrame(
            {
                "spot_id": obs.index[selected_indices].astype(str),
                "fullres_x": coords[selected_indices, 0],
                "fullres_y": coords[selected_indices, 1],
                "hires_x": hires_x,
                "hires_y": hires_y,
                "all_cell_adjusted_program_score_z": residual_z[selected_indices],
            }
        )
        overlay.insert(0, "legacy_program_id", "hepatocytes::8")
        overlay.insert(1, "program_uid", OVERLAY_PROGRAM_UID)
        overlay.insert(2, "display_name", "Stromal ECM")
        overlay.insert(3, "source_dataset_id", "GSE192741")
        overlay.insert(4, "display_source", "Guilliams et al.")
        overlay.insert(5, "reporting_unit_id", sample_id)
        overlay.insert(6, "source_individual_label", "H37")
        overlay.insert(7, "condition", "Steatotic")
        overlay["display_value_clipped_to"] = "[-2,2]"

        metadata = {
            "legacy_program_id": "hepatocytes::8",
            "program_uid": OVERLAY_PROGRAM_UID,
            "display_name": "Stromal ECM",
            "source_dataset_id": "GSE192741",
            "display_source": "Guilliams et al.",
            "reporting_unit_id": sample_id,
            "source_individual_label": "H37",
            "condition": "Steatotic",
            "n_graph_eligible_spots": len(selected_indices),
            "n_genes_measured": len(measured),
            "retained_l1_weight": retained,
            "selection_rule": expected_rule,
            "program_outcomes_used_for_selection": False,
            "residualization": (
                "all_16_cell2location_q05_plus_log1p_total_counts_and_"
                "log1p_detected_genes"
            ),
            "display_standardization": (
                "z_over_all_graph_eligible_spots_within_GSE192741"
            ),
            "rederived_residual_moran_i": observed,
            "sealed_residual_moran_i": float(
                sealed_result.iloc[0]["residual_moran_i"]
            ),
            "histology_image_relative_path": image_path.relative_to(BASE).as_posix(),
            "histology_image_sha256": sha256_file(image_path),
            "histology_image_width": image.width,
            "histology_image_height": image.height,
            "scalefactors_relative_path": scale_path.relative_to(BASE).as_posix(),
            "scalefactors_sha256": sha256_file(scale_path),
            "tissue_hires_scalef": hires_scale,
            "spot_diameter_fullres": float(scale_factors["spot_diameter_fullres"]),
            "input_h5ad_relative_path": input_path.relative_to(BASE).as_posix(),
            "input_h5ad_expected_sha256": str(input_manifest.iloc[0]["sha256"]),
            "map_role": "illustrative_registered_histology_overlay_not_inference",
            "histology_annotation_status": "no_pathologist_region_masks_available",
            "inference_note": (
                "weighted_program_score_not_named_gene_expression; matched-null_"
                "inference_is_in_the_adjacent_calibration_plot"
            ),
        }
        return overlay, image, metadata
    finally:
        del adata


def write_sidecar(data: pd.DataFrame) -> None:
    columns = [
        "legacy_program_id",
        "program_uid",
        "display_name",
        "source_dataset_id",
        "display_source",
        "source_citation",
        "source_order",
        "biological_unit",
        "evidence_state",
        "residualization",
        "residual_moran_i",
        "residual_null_mean",
        "residual_null_sd",
        "moran_excess",
        "bh_q",
        "n_null",
        "supported_within_source",
        "display_state",
        "matched_set_sha256",
        "inference_note",
    ]
    output = PANELS / "data/fig5f_spatial_program_calibration.tsv"
    output.parent.mkdir(parents=True, exist_ok=True)
    data[columns].to_csv(output, sep="\t", index=False)


def write_overlay_sidecars(
    overlay: pd.DataFrame, metadata: dict[str, object]
) -> None:
    data_dir = PANELS / "data"
    data_dir.mkdir(parents=True, exist_ok=True)
    overlay.to_csv(
        data_dir / "fig5f_spatial_histology_overlay.tsv.gz",
        sep="\t",
        index=False,
        compression={"method": "gzip", "mtime": 0},
    )
    pd.DataFrame([metadata]).to_csv(
        data_dir / "fig5f_spatial_histology_overlay_metadata.tsv",
        sep="\t",
        index=False,
    )


def render(
    data: pd.DataFrame,
    overlay: pd.DataFrame,
    histology: Image.Image,
) -> Path:
    fig = plt.figure(figsize=(4.85, 1.62))
    map_ax = fig.add_axes([0.015, 0.28, 0.32, 0.54])
    axes = [
        fig.add_axes([0.455, 0.255, 0.225, 0.55]),
        fig.add_axes([0.755, 0.255, 0.225, 0.55]),
    ]
    y_positions = {"GSE192741": 1.0, "Vu_et_al_2025": 0.0}
    x_max = 0.075

    for index, (program_uid, metadata) in enumerate(PROGRAMS.items()):
        ax = axes[index]
        part = data[data["program_uid"] == program_uid]
        ax.axvline(0, color=ZERO, linewidth=0.55, linestyle=(0, (2, 2)), zorder=0)
        for row in part.itertuples(index=False):
            y = y_positions[row.source_dataset_id]
            color = SUPPORTED if row.supported_within_source else INDETERMINATE
            # Horizontal bars, not stem+ball: lollipops are barred by house style.
            ax.barh(
                y,
                row.moran_excess,
                height=0.42,
                left=0,
                color=color,
                edgecolor="none",
                zorder=1,
            )
            ax.text(
                min(row.moran_excess + 0.0022, x_max - 0.020),
                y,
                q_label(row.bh_q),
                ha="left",
                va="center",
                color=INK,
            )

        ax.set_title(metadata["display_name"], pad=5, color=INK)
        ax.set_xlim(-0.003, x_max)
        ax.set_ylim(-0.45, 1.45)
        ax.set_xticks([0, 0.02, 0.04, 0.06])
        ax.set_xticklabels(["0", ".02", ".04", ".06"])
        ax.tick_params(axis="x", length=2.0, width=0.35, pad=1.5, colors=INK)
        ax.tick_params(axis="y", length=0, pad=2.5, colors=INK)
        for spine in ("top", "right", "left"):
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_linewidth(0.35)
        ax.spines["bottom"].set_color(INK)

    axes[0].set_yticks([1.0, 0.0])
    axes[0].set_yticklabels(["Guilliams et al.", "Vu et al."])
    axes[1].tick_params(labelleft=False)
    fig.text(
        0.72,
        0.055,
        "Moran excess over matched genes",
        ha="center",
        va="bottom",
        color=INK,
    )

    histology_rgb = np.asarray(histology, dtype=float) / 255.0
    luminance = (
        0.2126 * histology_rgb[:, :, 0]
        + 0.7152 * histology_rgb[:, :, 1]
        + 0.0722 * histology_rgb[:, :, 2]
    )
    histology_muted = 0.68 * histology_rgb + 0.32 * luminance[:, :, None]
    margin = 45
    x0 = max(0, int(np.floor(overlay["hires_x"].min())) - margin)
    x1 = min(histology.width, int(np.ceil(overlay["hires_x"].max())) + margin)
    y0 = max(0, int(np.floor(overlay["hires_y"].min())) - margin)
    y1 = min(histology.height, int(np.ceil(overlay["hires_y"].max())) + margin)
    crop = histology_muted[y0:y1, x0:x1]
    map_ax.imshow(crop, origin="upper", interpolation="bilinear", zorder=0)
    map_ax.scatter(
        overlay["hires_x"] - x0,
        overlay["hires_y"] - y0,
        c=overlay["all_cell_adjusted_program_score_z"],
        cmap=PROGRAM_CMAP,
        norm=PROGRAM_NORM,
        s=2.4,
        alpha=0.68,
        edgecolor="none",
        linewidth=0,
        zorder=1,
    )
    map_ax.set_xlim(0, x1 - x0)
    map_ax.set_ylim(y1 - y0, 0)
    map_ax.set_aspect("equal")
    map_ax.set_xticks([])
    map_ax.set_yticks([])
    for spine in map_ax.spines.values():
        spine.set_visible(False)
    map_ax.set_title("Stromal ECM", pad=5, color=INK)
    colorbar_ax = fig.add_axes([0.045, 0.125, 0.25, 0.025])
    colorbar = fig.colorbar(
        plt.cm.ScalarMappable(norm=PROGRAM_NORM, cmap=PROGRAM_CMAP),
        cax=colorbar_ax,
        orientation="horizontal",
    )
    colorbar.set_ticks([-2, 0, 2])
    colorbar.ax.tick_params(width=0.35, length=1.5, pad=1.0, colors=INK)
    colorbar.outline.set_linewidth(0.35)
    colorbar.set_label("Adjusted program score (z)", labelpad=1.0, color=INK)

    output = PANELS / "fig5f_spatial_program_calibration.pdf"
    output.parent.mkdir(parents=True, exist_ok=True)
    # NO bbox_inches="tight": it retrims a hand-placed add_axes layout, so the
    # emitted page box stops matching the declared figsize and the panel places
    # at the wrong size in Illustrator. Without it the page is exactly figsize.
    fig.savefig(output, dpi=400)
    plt.close(fig)
    return output


def main() -> None:
    data = load_source()
    overlay, histology, metadata = load_histology_overlay()
    write_sidecar(data)
    write_overlay_sidecars(overlay, metadata)
    output = render(data, overlay, histology)
    print(f"[fig5f calibration] saved: {output}")
    print(
        "CAPTION (Fig. 5F): All-cell-adjusted Moran excess over matched genes for "
        "the Stromal ECM and Ductular injury programs in Guilliams et al. and Vu et al. "
        "Bar length is observed residual Moran's I minus the matched-null mean; labels are "
        "BH-adjusted empirical matched-null q-values from 9,999 draws. Stromal ECM is "
        "supported in both sources, with Vu et al. source-dependent because its arrays "
        "lack a public donor key. Ductular injury is indeterminate in both sources. "
        "The registered H&E overlay shows the all-cell-adjusted Stromal ECM score in "
        "one outcome-blind representative Guilliams et al. section; it is illustrative "
        "and is not a pathologist-annotated histologic localization."
    )


if __name__ == "__main__":
    main()
