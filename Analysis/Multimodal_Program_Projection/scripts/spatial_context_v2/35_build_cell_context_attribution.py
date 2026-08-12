#!/usr/bin/env python3
# KEY MESSAGE: the two frozen programs occupy reproducible but composition-linked
# spatial contexts; these descriptive associations do not establish lineage origin.
"""Build a candidate-only, descriptive Visium cell-context association lane.

The producer deliberately reuses the byte-pinned spatial engine and the frozen
two-program registry.  It computes the existing hepatocyte-abundance/QC-
residualized score, then reports within-section/array Spearman correlations with
all 16 shared cell2location q05 abundance factors.  No spot-level P values are
calculated.  GSE192741 repeated sections are collapsed donor-first; Vu remains
source-dependent at the physical-array level.
"""

from __future__ import annotations

import argparse
import json
import os
import platform
import shutil
import sys
from pathlib import Path

import anndata as ad
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
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
import numpy as np
import pandas as pd
from scipy.stats import rankdata

from spatial_resource_lib import sha256_file, write_tsv
from visium_rerun_lib import (
    RELEASE_ID as PROGRAM_RELEASE_ID,
    build_paths,
    load_legacy_engine,
    prepare_engine_inputs,
    read_tsv,
    utc_now,
)


CANDIDATE_RELEASE_ID = "spatial-cell-context-attribution-candidate-2026-08-11-r2"
RESOURCE_RELEASE_ID = "spatial-resource-candidate-2026-08-11"
NATIVE_V2_READY_SHA256 = "fcff53888cd8a29adc817752be4f3798e7fb3bcd7903f02622a90af729089898"
EXPECTED_SOURCE_HASHES = {
    "GSE192741": "35f9f39848942f3b1c4c4d7670bd35ec90fccd45ece3a919b52d44ef02c66794",
    "Vu_et_al_2025": "36d6e69d8653c60103d09331497c8eab12083c3c864fc96ede641ef575ea3905",
}
PROGRAM_SHORT = {
    "hotspot_hepatocytes_f05c535ae5bbc0b9": "IGFBP7 program",
    "hotspot_hepatocytes_48f39dd4d817a10e": "BICC1 program",
}
FACTOR_ORDER = [
    "Hepatocytes",
    "Cholangiocytes",
    "Fibroblasts",
    "Endothelial cells",
    "Macrophages",
    "Mono+mono derived cells",
    "Neutrophils",
    "cDC1s",
    "cDC2s",
    "pDCs",
    "B cells",
    "Plasma cells",
    "T cells",
    "Circulating NK/NKT",
    "Resident NK",
    "Basophils",
]
FACTOR_COLUMNS = [f"q05__{index:02d}" for index in range(len(FACTOR_ORDER))]
CYAN_DARK = "#003A4F"
MAGENTA_DARK = "#7C256F"
GRAY = "#9E9E9E"
INK = "#222222"


def graph_indices(graphs: dict[str, list[tuple[np.ndarray, np.ndarray, np.ndarray]]]) -> np.ndarray:
    parts = [component[0] for arrays in graphs.values() for component in arrays]
    if not parts:
        raise RuntimeError("no graph-eligible tissue spots")
    return np.unique(np.concatenate(parts).astype(int))


def sample_graph_indices(
    graphs: dict[str, list[tuple[np.ndarray, np.ndarray, np.ndarray]]], sample_id: str
) -> np.ndarray:
    parts = [component[0] for component in graphs.get(sample_id, [])]
    if not parts:
        raise RuntimeError(f"no graph-eligible spots for {sample_id}")
    return np.unique(np.concatenate(parts).astype(int))


def spearman_rho(left: np.ndarray, right: np.ndarray) -> float:
    """Tie-aware Spearman rho without calculating an inferential P value."""
    x = np.asarray(left, dtype=float)
    y = np.asarray(right, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if int(keep.sum()) < 3:
        return float("nan")
    rx = rankdata(x[keep], method="average")
    ry = rankdata(y[keep], method="average")
    if np.std(rx, ddof=1) <= 0 or np.std(ry, ddof=1) <= 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def sign_label(value: float) -> str:
    if not np.isfinite(value):
        return "not_estimable"
    if value > 0:
        return "positive"
    if value < 0:
        return "negative"
    return "zero"


def require_upstream(paths, resource_root: Path) -> None:
    native_ready = paths.native_root / "v2_candidate/READY"
    if sha256_file(native_ready) != NATIVE_V2_READY_SHA256:
        raise RuntimeError("native two-program READY drift")
    for relative in ("validation/READY", "figures/READY", "effects/READY"):
        if not (resource_root / relative).is_file():
            raise RuntimeError(f"spatial Resource candidate is not sealed: {relative}")
    ready = pd.read_csv(resource_root / "validation/READY", sep="\t", dtype=str)
    if len(ready) != 1 or ready.iloc[0]["canonical_promotion_authorized"] != "FALSE":
        raise RuntimeError("cell-context lane requires the validated, unpromoted Resource candidate")
    if ready.iloc[0]["release_id"] != RESOURCE_RELEASE_ID:
        raise RuntimeError("spatial Resource release ID drift")
    for dataset, expected in EXPECTED_SOURCE_HASHES.items():
        observed = sha256_file(paths.datasets[dataset])
        if observed != expected:
            raise RuntimeError(f"{dataset} source h5ad drift: {observed}")


def program_identity(paths) -> dict[str, dict[str, str]]:
    rows = read_tsv(paths.hotspot_root / "program_registry_v2.tsv")
    selected = {row["program_uid"]: row for row in rows if row["program_uid"] in PROGRAM_SHORT}
    if set(selected) != set(PROGRAM_SHORT):
        raise RuntimeError("frozen two-program identities are incomplete")
    return selected


def score_dataset(paths, engine, registry, membership, identities, dataset: str):
    adata = ad.read_h5ad(paths.datasets[dataset])
    try:
        adata.var_names = adata.var_names.astype(str)
        obs = adata.obs.copy()
        for column in ("sample_id", "individual", "condition"):
            obs[column] = obs[column].astype(str)
        counts = engine.dense_counts(adata)
        library = pd.to_numeric(obs["total_counts"], errors="coerce").to_numpy(float)
        norm = engine.log_normalize(counts, library)
        detection = np.asarray((counts > 0).mean(axis=0)).ravel()
        gene_to_index = {gene: index for index, gene in enumerate(adata.var_names)}
        factor_names = [
            str(value).replace("means_per_cluster_mu_fg_", "")
            for value in adata.uns["mod"]["factor_names"]
        ]
        if len(factor_names) != 16 or set(factor_names) != set(FACTOR_ORDER):
            raise RuntimeError(f"{dataset} does not have the sealed 16-factor cell2location axis")
        abundance = np.asarray(adata.obsm["q05_cell_abundance_w_sf"], dtype=float)
        if abundance.shape != (adata.n_obs, 16):
            raise RuntimeError(f"{dataset} q05 abundance shape drift: {abundance.shape}")
        abundance = abundance[:, [factor_names.index(name) for name in FACTOR_ORDER]]
        if not np.isfinite(abundance).all() or (abundance < 0).any():
            raise RuntimeError(f"{dataset} q05 abundance contains invalid values")
        graphs, graph_audit = engine.section_graphs(obs, np.asarray(adata.obsm["spatial"]), k=6)
        eligible_all = graph_indices(graphs)
        technical_ids = sorted(graphs)
        if dataset == "GSE192741" and len(technical_ids) != 5:
            raise RuntimeError("GSE192741 must contain five eligible sections")
        if dataset == "Vu_et_al_2025" and len(technical_ids) != 10:
            raise RuntimeError("Vu must contain ten eligible physical arrays")

        abundance_rows = []
        for index in eligible_all:
            row = {
                "candidate_release_id": CANDIDATE_RELEASE_ID,
                "program_release_id": PROGRAM_RELEASE_ID,
                "dataset": dataset,
                "technical_id": str(obs.iloc[index]["sample_id"]),
                "source_individual_label": str(obs.iloc[index]["individual"]),
                "spot_id": str(obs.index[index]),
            }
            row.update({column: f"{abundance[index, j]:.16g}" for j, column in enumerate(FACTOR_COLUMNS)})
            abundance_rows.append(row)

        score_rows = []
        correlation_rows = []
        scoring_audit = []
        for program in registry.itertuples(index=False):
            members = membership.loc[membership["program_id"] == program.program_id].copy()
            members = members.groupby("gene_symbol", as_index=False)["original_l1_weight"].sum()
            members["present"] = members["gene_symbol"].isin(gene_to_index)
            members["detected"] = members["gene_symbol"].map(
                lambda gene: detection[gene_to_index[gene]] >= 0.01 if gene in gene_to_index else False
            )
            measured = members.loc[members["present"] & members["detected"]].copy()
            retained = float(measured["original_l1_weight"].sum())
            if len(measured) < 8 or retained < 0.20:
                raise RuntimeError(f"{dataset}/{program.program_id} is not testable")
            weights = measured["original_l1_weight"].to_numpy(float) / retained
            indices = np.asarray([gene_to_index[gene] for gene in measured["gene_symbol"]], dtype=int)
            raw_score = engine.extract_z(norm, indices) @ weights
            design = engine.design_matrix(obs, abundance[:, FACTOR_ORDER.index("Hepatocytes")])
            residual = np.asarray(engine.residualize(raw_score, design), dtype=float)
            center = float(np.mean(residual[eligible_all]))
            spread = float(np.std(residual[eligible_all], ddof=1))
            if not np.isfinite(spread) or spread <= 0:
                raise RuntimeError(f"{dataset}/{program.program_id} residual score is degenerate")
            residual_z = (residual - center) / spread
            identity = identities[program.program_id]
            for index in eligible_all:
                score_rows.append({
                    "candidate_release_id": CANDIDATE_RELEASE_ID,
                    "program_release_id": PROGRAM_RELEASE_ID,
                    "dataset": dataset,
                    "technical_id": str(obs.iloc[index]["sample_id"]),
                    "source_individual_label": str(obs.iloc[index]["individual"]),
                    "spot_id": str(obs.index[index]),
                    "program_uid": program.program_id,
                    "program_label": identity["module_name"],
                    "membership_sha256": identity["membership_sha256"],
                    "residual_program_score_z": f"{residual_z[index]:.16g}",
                })
            for technical_id in technical_ids:
                idx = sample_graph_indices(graphs, technical_id)
                donor_values = obs.iloc[idx]["individual"].astype(str).unique().tolist()
                if len(donor_values) != 1:
                    raise RuntimeError(f"{dataset}/{technical_id} maps to multiple source labels")
                for factor_index, factor in enumerate(FACTOR_ORDER):
                    rho = spearman_rho(residual_z[idx], abundance[idx, factor_index])
                    if not np.isfinite(rho):
                        raise RuntimeError(f"non-estimable correlation: {dataset}/{technical_id}/{factor}")
                    correlation_rows.append({
                        "candidate_release_id": CANDIDATE_RELEASE_ID,
                        "program_release_id": PROGRAM_RELEASE_ID,
                        "dataset": dataset,
                        "source_dependence": "independent" if dataset == "GSE192741" else "source_dependent",
                        "biological_unit_resolution": "resolved_human_donor" if dataset == "GSE192741" else "unresolved_physical_array",
                        "technical_unit": "Visium_section" if dataset == "GSE192741" else "Visium_array",
                        "technical_id": technical_id,
                        "source_individual_label": donor_values[0],
                        "program_uid": program.program_id,
                        "program_label": identity["module_name"],
                        "factor_order": factor_index + 1,
                        "cell2location_factor": factor,
                        "association_method": "within_physical_unit_spearman_rho",
                        "spearman_rho": f"{rho:.16g}",
                        "direction": sign_label(rho),
                        "n_graph_eligible_spots": len(idx),
                        "spot_level_inferential_pvalue_authorized": "FALSE",
                        "interpretation": "descriptive_context_covariation_from_shared_visium_rna_matrix_not_lineage_origin_or_mechanism",
                    })
            scoring_audit.append({
                "candidate_release_id": CANDIDATE_RELEASE_ID,
                "dataset": dataset,
                "program_uid": program.program_id,
                "program_label": identity["module_name"],
                "membership_sha256": identity["membership_sha256"],
                "n_graph_eligible_spots": len(eligible_all),
                "n_genes_measured": len(measured),
                "retained_l1_weight": f"{retained:.16g}",
                "score_definition": "frozen_positive_weight_gene_z_residualized_for_hepatocyte_q05_log_umi_log_detected_genes_then_dataset_z",
                "association_method": "within_physical_unit_spearman_rho_no_pvalue",
            })
        return abundance_rows, score_rows, correlation_rows, scoring_audit, graph_audit
    finally:
        del adata


def build_unit_summary(correlations: pd.DataFrame) -> pd.DataFrame:
    columns = ["dataset", "program_uid", "program_label", "factor_order", "cell2location_factor"]
    rows = []
    for dataset in ("GSE192741", "Vu_et_al_2025"):
        part = correlations[correlations["dataset"] == dataset].copy()
        if dataset == "GSE192741":
            group_columns = columns + ["source_individual_label"]
            for key, group in part.groupby(group_columns, sort=False, observed=True):
                value = float(group["spearman_rho"].mean())
                rows.append({
                    "candidate_release_id": CANDIDATE_RELEASE_ID,
                    "program_release_id": PROGRAM_RELEASE_ID,
                    "dataset": key[0],
                    "source_dependence": "independent",
                    "summary_unit": "human_donor",
                    "summary_unit_id": key[-1],
                    "program_uid": key[1],
                    "program_label": key[2],
                    "factor_order": key[3],
                    "cell2location_factor": key[4],
                    "spearman_rho": value,
                    "direction": sign_label(value),
                    "n_technical_units_collapsed": group["technical_id"].nunique(),
                    "collapse_rule": "equal_weight_mean_of_within_section_spearman_rho",
                    "inferential_pvalue_authorized": "FALSE",
                })
        else:
            for row in part.itertuples(index=False):
                rows.append({
                    "candidate_release_id": CANDIDATE_RELEASE_ID,
                    "program_release_id": PROGRAM_RELEASE_ID,
                    "dataset": dataset,
                    "source_dependence": "source_dependent",
                    "summary_unit": "unresolved_physical_array",
                    "summary_unit_id": row.technical_id,
                    "program_uid": row.program_uid,
                    "program_label": row.program_label,
                    "factor_order": row.factor_order,
                    "cell2location_factor": row.cell2location_factor,
                    "spearman_rho": float(row.spearman_rho),
                    "direction": row.direction,
                    "n_technical_units_collapsed": 1,
                    "collapse_rule": "none_array_is_the_unresolved_reporting_unit",
                    "inferential_pvalue_authorized": "FALSE",
                })
    return pd.DataFrame(rows)


def build_dataset_summary(unit_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    keys = ["dataset", "program_uid", "program_label", "factor_order", "cell2location_factor"]
    for key, group in unit_summary.groupby(keys, sort=False, observed=True):
        values = group["spearman_rho"].to_numpy(float)
        source_dependence = "independent" if key[0] == "GSE192741" else "source_dependent"
        rows.append({
            "candidate_release_id": CANDIDATE_RELEASE_ID,
            "program_release_id": PROGRAM_RELEASE_ID,
            "dataset": key[0],
            "source_dependence": source_dependence,
            "summary_unit": "human_donor" if key[0] == "GSE192741" else "unresolved_physical_array",
            "program_uid": key[1],
            "program_label": key[2],
            "factor_order": key[3],
            "cell2location_factor": key[4],
            "median_spearman_rho": float(np.median(values)),
            "q25_spearman_rho": float(np.quantile(values, 0.25)),
            "q75_spearman_rho": float(np.quantile(values, 0.75)),
            "minimum_spearman_rho": float(np.min(values)),
            "maximum_spearman_rho": float(np.max(values)),
            "n_reporting_units": len(values),
            "n_positive": int(np.sum(values > 0)),
            "n_negative": int(np.sum(values < 0)),
            "n_zero": int(np.sum(values == 0)),
            "median_direction": sign_label(float(np.median(values))),
            "inferential_pvalue_authorized": "FALSE",
            "interpretation": "descriptive_distribution_across_reporting_units_with_shared_transcriptome_dependency",
        })
    return pd.DataFrame(rows)


def build_concordance(dataset_summary: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for uid in PROGRAM_SHORT:
        for order, factor in enumerate(FACTOR_ORDER, start=1):
            selected = dataset_summary[
                (dataset_summary["program_uid"] == uid)
                & (dataset_summary["cell2location_factor"] == factor)
            ]
            if set(selected["dataset"]) != {"GSE192741", "Vu_et_al_2025"}:
                raise RuntimeError(f"cross-source summary incomplete for {uid}/{factor}")
            gse = selected[selected["dataset"] == "GSE192741"].iloc[0]
            vu = selected[selected["dataset"] == "Vu_et_al_2025"].iloc[0]
            gse_value = float(gse["median_spearman_rho"])
            vu_value = float(vu["median_spearman_rho"])
            rows.append({
                "candidate_release_id": CANDIDATE_RELEASE_ID,
                "program_release_id": PROGRAM_RELEASE_ID,
                "program_uid": uid,
                "program_label": gse["program_label"],
                "factor_order": order,
                "cell2location_factor": factor,
                "gse_donor_median_spearman_rho": gse_value,
                "gse_direction": sign_label(gse_value),
                "vu_array_median_spearman_rho": vu_value,
                "vu_direction": sign_label(vu_value),
                "direction_concordant": str(sign_label(gse_value) == sign_label(vu_value)).upper(),
                "vu_source_dependent": "TRUE",
                "inferential_pvalue_authorized": "FALSE",
                "interpretation": "descriptive_cross_source_direction_only_not_replication_test",
            })
    return pd.DataFrame(rows)


def external_source_audit(project_root: Path, resource_root: Path) -> pd.DataFrame:
    context_root = project_root / "Analysis/Multimodal_Program_Projection/candidates" / PROGRAM_RELEASE_ID / "spatial_context"
    cosmx_path = context_root / "real_adapters/govaere2026_cosmx_il32_context/source_native_context.tsv"
    cosmx = pd.read_csv(cosmx_path, sep="\t", dtype=str)
    cosmx["is_mash_bool"] = cosmx["is_mash"].str.upper().eq("TRUE")
    mash = cosmx[cosmx["is_mash_bool"]]
    effects = pd.read_csv(resource_root / "effects/spatial_program_effects.tsv", sep="\t", dtype=str)
    cosmx_effects = effects[effects["dataset_id"] == "Govaere2026_CosMx"]
    yak_effects = effects[effects["dataset_id"] == "Yakubovsky2026"]
    yak_root = project_root / "Analysis/Spatial/candidates" / PROGRAM_RELEASE_ID / "yakubovsky2026"
    yak_gate = pd.read_csv(yak_root / "gate_status.tsv", sep="\t", dtype=str).iloc[0]
    existing_cosmx_panel = project_root / "scripts/figures/figS4c_il32_macrophage_cosmx.R"
    rows = [
        {
            "candidate_release_id": CANDIDATE_RELEASE_ID,
            "dataset": "Govaere2026_CosMx",
            "source_rows_audited": len(cosmx),
            "eligible_biological_or_physical_units": len(mash),
            "audit_metric_1": "MASH_arrays_with_IL32_high_hepatocytes_closer_to_macrophages",
            "audit_value_1": int((pd.to_numeric(mash["estimate"]) < 0).sum()),
            "audit_metric_2": "MASH_arrays_with_positive_IL32_vs_neighbour_CD74_rho",
            "audit_value_2": int((pd.to_numeric(mash["secondary_estimate"]) > 0).sum()),
            "resource_program_states": ";".join(sorted(cosmx_effects["evidence_state"].unique())),
            "existing_panel_or_output": str(existing_cosmx_panel.relative_to(project_root)),
            "additional_analysis_decision": "no_new_analysis_or_panel",
            "reason": "four exact array-level source-native rows and an existing supplementary panel already expose the narrow IL32 proximity vignette; whole-program scoring is not applicable and cell-level inference would pseudoreplicate",
            "permitted_claim": "source_dependent_IL32_macrophage_proximity_context",
            "prohibited_claim": "whole_program_validation_or_cell_level_inferential_replication",
        },
        {
            "candidate_release_id": CANDIDATE_RELEASE_ID,
            "dataset": "Yakubovsky2026",
            "source_rows_audited": len(yak_effects),
            "eligible_biological_or_physical_units": int(yak_gate["n_source_gate_donors"]),
            "audit_metric_1": "confirmatory_programs_indeterminate_in_resource_release",
            "audit_value_1": int((yak_effects["evidence_state"] == "indeterminate").sum()),
            "audit_metric_2": "programs_supported_in_source_defined_lipid_context",
            "audit_value_2": int((yak_effects["evidence_state"] == "supported").sum()),
            "resource_program_states": ";".join(sorted(yak_effects["evidence_state"].unique())),
            "existing_panel_or_output": str((yak_root / "zonation_reference.tsv").relative_to(project_root)),
            "additional_analysis_decision": "no_new_analysis_or_panel",
            "reason": "the three-donor source-defined binary lipid analysis is already represented as indeterminate normal-liver context; another panel would invite MASLD-validation or lipid-free-zone overinterpretation",
            "permitted_claim": "normal_liver_zonation_and_source_defined_lipid_context_boundary",
            "prohibited_claim": "MASLD_validation_dose_response_or_non_lipid_zone_is_lipid_free",
        },
    ]
    return pd.DataFrame(rows)


def render_heatmap(summary: pd.DataFrame, output: Path) -> pd.DataFrame:
    column_order = [
        ("GSE192741", "hotspot_hepatocytes_f05c535ae5bbc0b9"),
        ("GSE192741", "hotspot_hepatocytes_48f39dd4d817a10e"),
        ("Vu_et_al_2025", "hotspot_hepatocytes_f05c535ae5bbc0b9"),
        ("Vu_et_al_2025", "hotspot_hepatocytes_48f39dd4d817a10e"),
    ]
    records = []
    matrix = np.zeros((len(FACTOR_ORDER), len(column_order)), dtype=float)
    for row_index, factor in enumerate(FACTOR_ORDER):
        for column_index, (dataset, uid) in enumerate(column_order):
            row = summary[
                (summary["dataset"] == dataset)
                & (summary["program_uid"] == uid)
                & (summary["cell2location_factor"] == factor)
            ]
            if len(row) != 1:
                raise RuntimeError(f"heatmap family incomplete for {dataset}/{uid}/{factor}")
            record = row.iloc[0].to_dict()
            record["heatmap_column_order"] = column_index + 1
            record["heatmap_column_label"] = PROGRAM_SHORT[uid]
            records.append(record)
            matrix[row_index, column_index] = float(record["median_spearman_rho"])
    source = pd.DataFrame(records)
    source.to_csv(output / "data/figS_cell_context_association.tsv", sep="\t", index=False)

    limit = max(0.1, float(np.max(np.abs(matrix))))
    cmap = LinearSegmentedColormap.from_list("spatial_context_diverging", [MAGENTA_DARK, "#FFFFFF", CYAN_DARK])
    norm = TwoSlopeNorm(vmin=-limit, vcenter=0.0, vmax=limit)
    fig, ax = plt.subplots(figsize=(4.75, 4.35))
    image = ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(4), ["IGFBP7", "BICC1", "IGFBP7", "BICC1"])
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(len(FACTOR_ORDER)), FACTOR_ORDER)
    ax.tick_params(length=0, pad=2)
    for row_index in range(matrix.shape[0]):
        for column_index in range(matrix.shape[1]):
            value = matrix[row_index, column_index]
            color = "#FFFFFF" if abs(value) > limit * 0.58 else INK
            ax.text(column_index, row_index, f"{value:+.2f}", ha="center", va="center", color=color)
    ax.axvline(1.5, color="#FFFFFF", linewidth=1.2)
    ax.text(0.25, 1.105, "GSE192741 · 4 donors", transform=ax.transAxes, ha="center", va="bottom", color=INK)
    ax.text(0.75, 1.105, "Vu et al. · 10 arrays\nsource-dependent", transform=ax.transAxes, ha="center", va="bottom", color=INK)
    for spine in ax.spines.values():
        spine.set_visible(False)
    colorbar = fig.colorbar(image, ax=ax, fraction=0.04, pad=0.04)
    colorbar.set_label("median within-unit Spearman rho", rotation=90)
    colorbar.ax.tick_params(length=2)
    ax.text(
        0,
        -0.14,
        "Residual program expression vs q05 abundance from the same Visium RNA matrix.\nDescriptive only; no spot-level P values or lineage-origin inference.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        color=GRAY,
    )
    fig.subplots_adjust(left=0.34, right=0.90, top=0.82, bottom=0.15)
    fig.savefig(output / "panels/figS_cell_context_association.pdf", bbox_inches="tight", dpi=400)
    plt.close(fig)
    return source


def build(project_root: Path, output_root: Path) -> Path:
    project_root = project_root.resolve()
    output_root = output_root.resolve()
    if output_root.exists():
        raise RuntimeError(f"refusing to overwrite cell-context candidate: {output_root}")
    staging = output_root.with_name(f".{output_root.name}.incomplete.{os.environ.get('SLURM_JOB_ID', os.getpid())}")
    if staging.exists():
        raise RuntimeError(f"staging directory exists: {staging}")
    staging.mkdir(parents=True)
    (staging / "data").mkdir()
    (staging / "panels").mkdir()
    try:
        paths = build_paths(project_root)
        resource_root = project_root / "Analysis/Multimodal_Program_Projection/candidates" / RESOURCE_RELEASE_ID
        require_upstream(paths, resource_root)
        engine = load_legacy_engine(paths)
        registry, membership, universe = prepare_engine_inputs(paths, "v2")
        if len(registry) != 2 or len(universe) != 2 or set(registry["program_id"]) != set(PROGRAM_SHORT):
            raise RuntimeError("cell-context lane requires the complete sealed two-program family")
        identities = program_identity(paths)

        abundance_rows, score_rows, correlation_rows, audit_rows, graph_rows = [], [], [], [], []
        for dataset in ("GSE192741", "Vu_et_al_2025"):
            abundance, scores, correlations, scoring, graph = score_dataset(
                paths, engine, registry, membership, identities, dataset
            )
            abundance_rows.extend(abundance)
            score_rows.extend(scores)
            correlation_rows.extend(correlations)
            audit_rows.extend(scoring)
            graph = graph.copy()
            graph.insert(0, "dataset", dataset)
            graph_rows.append(graph)

        abundance_table = pd.DataFrame(abundance_rows)
        score_table = pd.DataFrame(score_rows)
        correlation_table = pd.DataFrame(correlation_rows)
        correlation_table["spearman_rho"] = pd.to_numeric(correlation_table["spearman_rho"])
        abundance_table.to_parquet(staging / "data/per_spot_cell2location_q05.parquet", index=False)
        score_table.to_parquet(staging / "data/per_spot_residual_program_scores.parquet", index=False)
        correlation_table.to_csv(staging / "data/per_physical_unit_cell_context.tsv", sep="\t", index=False)
        pd.DataFrame(audit_rows).to_csv(staging / "data/scoring_audit.tsv", sep="\t", index=False)
        pd.concat(graph_rows, ignore_index=True).to_csv(staging / "data/graph_audit.tsv", sep="\t", index=False)

        unit_summary = build_unit_summary(correlation_table)
        unit_summary.to_csv(staging / "data/per_reporting_unit_cell_context.tsv", sep="\t", index=False)
        dataset_summary = build_dataset_summary(unit_summary)
        dataset_summary.to_csv(staging / "data/dataset_cell_context_summary.tsv", sep="\t", index=False)
        concordance = build_concordance(dataset_summary)
        concordance.to_csv(staging / "data/cross_source_direction_concordance.tsv", sep="\t", index=False)
        audit = external_source_audit(project_root, resource_root)
        audit.to_csv(staging / "data/external_source_audit.tsv", sep="\t", index=False)
        render_heatmap(dataset_summary, staging)

        source_paths = {
            "producer": Path(__file__).resolve(),
            "pinned_scientific_engine": paths.v1_engine,
            "candidate_library": Path(__file__).with_name("visium_rerun_lib.py"),
            "v2_registry": paths.hotspot_root / "program_registry_v2.tsv",
            "v2_membership": paths.hotspot_root / "program_membership_v2.tsv",
            "native_v2_READY": paths.native_root / "v2_candidate/READY",
            "spatial_resource_READY": resource_root / "validation/READY",
            "GSE192741_source_h5ad": paths.datasets["GSE192741"],
            "Vu_source_h5ad": paths.datasets["Vu_et_al_2025"],
        }
        source_rows = []
        for role, path in source_paths.items():
            source_rows.append({
                "candidate_release_id": CANDIDATE_RELEASE_ID,
                "source_role": role,
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": sha256_file(path),
            })
        write_tsv(staging / "source_manifest.tsv", tuple(source_rows[0]), source_rows)
        write_tsv(
            staging / "execution_manifest.tsv",
            (
                "candidate_release_id", "program_release_id", "resource_release_id", "completed_utc",
                "python", "platform", "producer", "producer_sha256", "slurm_job_id", "n_datasets",
                "n_programs", "n_cell2location_factors", "n_physical_unit_rows", "n_reporting_unit_rows",
                "n_dataset_summary_rows", "n_cross_source_rows", "n_spot_score_rows", "n_spot_abundance_rows",
                "association_method", "pvalues_computed", "vu_source_dependent", "canonical_output_written",
            ),
            [{
                "candidate_release_id": CANDIDATE_RELEASE_ID,
                "program_release_id": PROGRAM_RELEASE_ID,
                "resource_release_id": RESOURCE_RELEASE_ID,
                "completed_utc": utc_now(),
                "python": sys.version.replace("\n", " "),
                "platform": platform.platform(),
                "producer": str(Path(__file__).resolve()),
                "producer_sha256": sha256_file(Path(__file__).resolve()),
                "slurm_job_id": os.environ.get("SLURM_JOB_ID", "not_slurm"),
                "n_datasets": 2,
                "n_programs": 2,
                "n_cell2location_factors": 16,
                "n_physical_unit_rows": len(correlation_table),
                "n_reporting_unit_rows": len(unit_summary),
                "n_dataset_summary_rows": len(dataset_summary),
                "n_cross_source_rows": len(concordance),
                "n_spot_score_rows": len(score_table),
                "n_spot_abundance_rows": len(abundance_table),
                "association_method": "within_physical_unit_spearman_rho_then_donor_first_or_array_summary",
                "pvalues_computed": "FALSE",
                "vu_source_dependent": "TRUE",
                "canonical_output_written": "FALSE",
            }],
        )
        (staging / "BUILD_COMPLETE").write_text(
            json.dumps({
                "candidate_release_id": CANDIDATE_RELEASE_ID,
                "status": "built_unvalidated",
                "canonical_output_written": False,
            }, sort_keys=True) + "\n",
            encoding="utf-8",
        )
        os.replace(staging, output_root)
        return output_root
    except Exception:
        if staging.exists():
            shutil.rmtree(staging)
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    args = parser.parse_args()
    try:
        result = build(args.project_root, args.output_root)
        print(f"BUILT {result}")
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
