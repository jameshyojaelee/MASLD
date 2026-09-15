#!/usr/bin/env python3
"""Build read-only source contracts for the aging-inspired Figure 4–6 rebuild."""

from __future__ import annotations

import hashlib
import json
import math
import os
import shutil
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
AMBIENT = (
    ROOT
    / "Analysis/SingleCell/candidates"
    / "ambient-program-recalibration-all117-candidate-2026-08-15-v2"
)
CROSS = (
    ROOT
    / "Analysis/SingleCell/candidates"
    / "cross-lineage-specificity-complete-atlas-candidate-2026-08-15-v2"
)
OUT = Path(os.environ["FIG_CAND_ROOT"])
SOURCE = OUT / "source_tables"

PROGRAM_RELEASE = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / "program-context-v2-stage-corrected-candidate-2026-08-13-v3"
)
REGISTRY = PROGRAM_RELEASE / "hotspot/program_registry_v2.tsv"
BULK = (
    ROOT
    / "figures/candidates/fig4-stage-terminology-corrected-2026-08-13-v5"
    / "source_tables/fig4e_bulk_projection.tsv"
)
SPATIAL_ROOT = (
    ROOT
    / "Analysis/Multimodal_Program_Projection/candidates"
    / "spatial-impact-figures-candidate-2026-08-11"
)
SPATIAL = SPATIAL_ROOT / "data/fig4f_matched_null_discrimination.tsv"
SPATIAL_MAP = SPATIAL_ROOT / "panels/fig4e_spatial_program_maps.pdf"
STAGE_ROOT = ROOT / "figures/candidates/pi-figure-redesign-2026-08-13-v8"
STAGE = STAGE_ROOT / "analysis/stage_extensions/stage_all_gene_results.tsv"
COHORT = STAGE_ROOT / "analysis/stage_extensions/cohort_disease_all_gene_results.tsv"
BLOCKERS = STAGE_ROOT / "BLOCKERS.tsv"
PROTEIN_ROOT = ROOT / "figures/main/fig5_molecular_context/data"
PROTEIN = PROTEIN_ROOT / "composite_mrna_protein_corrected.csv"
PROTEIN_HIST = PROTEIN_ROOT / "composite_protein_histology_partial_corrected.csv"

SUPPLEMENT_ASSETS = {
    "supplementary/figureS4/panels/figs4a_program_landscape_detail.pdf":
        STAGE_ROOT / "figure4/panels/fig4c_program_landscape.pdf",
    "supplementary/figureS4/panels/figs4b_fixed_communication.pdf":
        STAGE_ROOT / "figure4/panels/fig4d_communication.pdf",
    "supplementary/figureS4/panels/figs4c_tf_activity.pdf":
        STAGE_ROOT / "figure4/panels/fig4f_tf_activity.pdf",
    "supplementary/figureS4/source_tables/fig4c_all_117_programs.tsv":
        STAGE_ROOT / "source_tables/fig4c_all_117_programs.tsv",
    "supplementary/figureS4/source_tables/fig4c_selected_program_profiles.tsv":
        STAGE_ROOT / "source_tables/fig4c_selected_program_profiles.tsv",
    "supplementary/figureS4/source_tables/fig4d_fixed_13_pair_communication.tsv":
        STAGE_ROOT / "source_tables/fig4d_fixed_13_pair_communication.tsv",
    "supplementary/figureS4/source_tables/fig4f_all_tf_activity.tsv":
        STAGE_ROOT / "source_tables/fig4f_all_tf_activity.tsv",
    "supplementary/figureS5/panels/figs5a_multimodal_program_source_matrix.pdf":
        SPATIAL_ROOT / "panels/figS_multimodal_program_source_matrix.pdf",
    "supplementary/figureS5/panels/figs5b_full_composition_sensitivity.pdf":
        SPATIAL_ROOT / "panels/figS_full_composition_sensitivity.pdf",
    "supplementary/figureS5/source_tables/figS_multimodal_program_source_matrix.tsv":
        SPATIAL_ROOT / "data/figS_multimodal_program_source_matrix.tsv",
    "supplementary/figureS5/source_tables/figS_full_composition_sensitivity.tsv":
        SPATIAL_ROOT / "data/figS_full_composition_sensitivity.tsv",
}

HEROES = {
    "hotspot_hepatocytes_f05c535ae5bbc0b9": "ECM/IGFBP7",
    "hotspot_hepatocytes_48f39dd4d817a10e": "Ductular-injury/BICC1",
}
TRIAD = ("GNMT", "MAT1A", "CYP2C19")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def refuse(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)


def write_tsv(frame: pd.DataFrame, name: str) -> Path:
    path = SOURCE / name
    refuse(path)
    frame.to_csv(path, sep="\t", index=False)
    return path


def provenance(path: Path) -> dict[str, object]:
    return {
        "source_path": str(path.resolve()),
        "source_sha256": sha256(path),
    }


def clean_q(value: object) -> float:
    number = float(value)
    return max(number, np.finfo(float).tiny)


def compact_number(value: float, digits: int = 2) -> str:
    if not np.isfinite(value):
        return "NA"
    if abs(value) < 0.001 and value != 0:
        return f"{value:.1e}"
    return f"{value:.{digits}f}"


def evidence_row(
    example: str,
    node_order: int,
    node: str,
    modality: str,
    molecular_object: str,
    metric: str,
    estimate: float | None,
    unit: str,
    state: str,
    state_reason: str,
    source_path: Path,
    unresolved_alternative: str,
    next_experiment: str,
    short_display: str,
    auxiliary: str = "",
) -> dict[str, object]:
    return {
        "example": example,
        "node_order": node_order,
        "node": node,
        "modality": modality,
        "molecular_object": molecular_object,
        "metric": metric,
        "estimate": np.nan if estimate is None else estimate,
        "unit": unit,
        "state": state,
        "state_reason": state_reason,
        "unresolved_alternative": unresolved_alternative,
        "next_experiment": next_experiment,
        "short_display": short_display,
        "auxiliary": auxiliary,
        **provenance(source_path),
    }


def main() -> None:
    require((AMBIENT / "VALIDATED").exists(), "ambient candidate is not validated")
    require((CROSS / "VALIDATED").exists(), "complete cross-lineage candidate is not validated")
    inputs = [
        REGISTRY, BULK, SPATIAL, SPATIAL_MAP, STAGE, COHORT, BLOCKERS, PROTEIN,
        PROTEIN_HIST,
    ]
    inputs += list(SUPPLEMENT_ASSETS.values())
    inputs += [
        AMBIENT / "results/ambient_program_effects.tsv",
        AMBIENT / "results/ambient_dataset_lineage_qc.tsv",
        CROSS / "results/hero_lineage_contrasts.tsv",
        CROSS / "results/hero_lineage_promotion.tsv",
        AMBIENT / "results/validation_report.tsv",
        CROSS / "results/validation_report.tsv",
    ]
    for path in inputs:
        require(path.exists(), f"missing input: {path}")

    supplement_rows = []
    for relative_destination, source_path in SUPPLEMENT_ASSETS.items():
        destination = OUT / relative_destination
        refuse(destination)
        shutil.copy2(source_path, destination)
        require(sha256(destination) == sha256(source_path), "supplement copy hash drift")
        supplement_rows.append(
            {
                "relative_destination": relative_destination,
                "source_path": str(source_path.resolve()),
                "source_sha256": sha256(source_path),
                "destination_sha256": sha256(destination),
                "action": "byte_identical_move_to_supplement",
            }
        )

    registry = pd.read_csv(REGISTRY, sep="\t")
    effects_all = pd.read_csv(AMBIENT / "results/ambient_program_effects.tsv", sep="\t")
    effects = effects_all[
        effects_all["analysis_universe"] == "complete_case_common_universe"
    ].copy()
    require(len(effects) == 117, "primary ambient family is not 117 programs")
    require(set(effects["program_uid"]) == set(registry["program_uid"]), "program key drift")

    # 4A. Vector schematic contract. This is deliberately editable in Illustrator.
    nodes = pd.DataFrame(
        [
            (1, "Frozen biological object", "117 programs; memberships and weights fixed", 0.0, 1.0),
            (2, "Donor collapse", "Repeated runs pooled; donor is the replicate", 1.5, 1.0),
            (3, "Disease association", "score ~ cross-sectional stage + dataset", 3.0, 1.0),
            (4, "Ambient recalibration", "decontX by dataset; GSE189600 passthrough", 4.5, 1.0),
            (5, "Identity stress test", "paired corrected-minus-raw stage effect", 6.0, 1.0),
            (6, "Same-atlas localization", "fixed-weight lineage-minus-hepatocyte contrasts", 4.5, 0.0),
            (7, "Tissue transport", "fixed genes into bulk stage and spatial context", 6.0, 0.0),
        ],
        columns=["node_id", "title", "detail", "x", "y"],
    )
    edges = pd.DataFrame(
        [(1, 2), (2, 3), (3, 4), (4, 5), (3, 6), (5, 7), (6, 7)],
        columns=["from", "to"],
    )
    write_tsv(nodes, "fig4a_analysis_logic_nodes.tsv")
    write_tsv(edges, "fig4a_analysis_logic_edges.tsv")
    write_tsv(pd.DataFrame(supplement_rows), "supplement_asset_manifest.tsv")

    # 4C. Complete all-program disease skyline.
    direct_labels = {
        "positive_control_leukocyte_DOCK2": "DOCK2 control",
        "positive_control_endothelial_STAB2": "STAB2 control",
        "intrinsic_control_secretory_ALB": "ALB control",
        "intrinsic_control_xenobiotic_CYP": "CYP control",
        "hero_ECM_IGFBP7": "ECM/IGFBP7",
        "hero_ductular_BICC1": "Ductular-injury/BICC1",
    }
    skyline = effects.copy()
    skyline["disease_beta"] = skyline["raw_beta"]
    skyline["disease_qvalue"] = skyline["raw_qvalue"]
    skyline["family_evidence"] = -np.log10(skyline["disease_qvalue"].map(clean_q))
    skyline["direct_label"] = skyline["display_role"].map(direct_labels).fillna("")
    skyline["label_requested"] = (
        skyline["display_role"].ne("") | skyline["robust_display_frozen"].astype(bool)
    )
    skyline["hc3_badge"] = skyline["raw_hc3_qvalue"].astype(float) < 0.05
    keep4c = [
        "program_uid", "cell_type", "module", "module_name", "disease_beta",
        "raw_se", "raw_ci_low", "raw_ci_high", "disease_qvalue",
        "raw_hc3_qvalue", "family_evidence", "evidence_state", "display_role",
        "label_requested", "direct_label", "hc3_badge", "n_model_donors",
        "failed_decontx_datasets", "ambient_burden_median",
    ]
    write_tsv(skyline[keep4c], "fig4c_all117_disease_skyline.tsv")

    # 4D. Ambient effect transport and explicit dataset status strip.
    transport = skyline[
        [
            "program_uid", "cell_type", "module", "module_name", "raw_beta",
            "raw_ci_low", "raw_ci_high", "raw_qvalue", "corrected_beta",
            "corrected_ci_low", "corrected_ci_high", "corrected_qvalue",
            "delta_beta", "delta_ci_low", "delta_ci_high", "delta_qvalue",
            "raw_hc3_qvalue", "corrected_hc3_qvalue", "delta_hc3_qvalue",
            "evidence_state", "display_role", "direct_label", "label_requested",
            "ambient_burden_median", "ambient_burden_min", "ambient_burden_max",
            "n_model_donors", "failed_decontx_datasets",
        ]
    ].copy()
    write_tsv(transport, "fig4d_ambient_effect_transport.tsv")
    qc = pd.read_csv(AMBIENT / "results/ambient_dataset_lineage_qc.tsv", sep="\t")
    status_rows = []
    for dataset, group in qc.groupby("dataset", sort=True):
        require(group["correction_status"].nunique() == 1, f"status varies by lineage: {dataset}")
        reason_column = "failure_detail" if "failure_detail" in group else "failure_reason"
        reasons = sorted(set(group[reason_column].dropna().astype(str)) - {""})
        status_rows.append(
            {
                "dataset": dataset,
                "correction_status": group["correction_status"].iloc[0],
                "failure_reason": "; ".join(reasons),
                "n_program_universe_cells": int(group["n_cells"].sum()),
                "n_biological_donors_max_across_lineages": int(group["n_biological_donors"].max()),
                "explicit_passthrough_label": (
                    "GSE189600: uncorrected passthrough"
                    if dataset == "GSE189600"
                    else ""
                ),
            }
        )
    status = pd.DataFrame(status_rows)
    require(len(status) == 7 and "GSE189600" in set(status["dataset"]), "dataset QC drift")
    write_tsv(status, "fig4d_dataset_status_strip.tsv")

    # 4E. All six lineages are visible; hepatocytes are the zero reference.
    contrasts = pd.read_csv(CROSS / "results/hero_lineage_contrasts.tsv", sep="\t")
    primary_contrasts = contrasts[
        (contrasts["annotation_filter"] == "all_annotated_cells")
        & (contrasts["minimum_cells"] == 50)
    ].copy()
    require(len(primary_contrasts) == 10, "primary lineage family is not ten")
    promotion = pd.read_csv(CROSS / "results/hero_lineage_promotion.tsv", sep="\t")
    require(len(promotion) == 2, "promotion table must contain two programs")
    family_destination = promotion["destination"].unique()
    require(len(family_destination) == 1, "lineage destination is inconsistent")
    baseline = []
    for uid, label in HEROES.items():
        baseline.append(
            {
                "program_uid": uid,
                "program_name": label,
                "comparison": "Hepatocytes reference",
                "comparison_lineage": "Hepatocytes",
                "annotation_filter": "all_annotated_cells",
                "minimum_cells": 50,
                "beta": 0.0,
                "se": 0.0,
                "pvalue": np.nan,
                "hc3_se": 0.0,
                "hc3_pvalue": np.nan,
                "ci_low": 0.0,
                "ci_high": 0.0,
                "hc3_ci_low": 0.0,
                "hc3_ci_high": 0.0,
                "n_paired_donors": np.nan,
                "n_datasets": np.nan,
                "eligible_datasets": "",
                "failure_reason": "reference lineage",
                "qvalue": np.nan,
                "hc3_qvalue": np.nan,
            }
        )
    lineage = pd.concat([primary_contrasts, pd.DataFrame(baseline)], ignore_index=True)
    lineage = lineage.merge(
        promotion[
            [
                "program_uid", "target_lineage", "program_pass",
                "main_figure_family_eligible", "destination",
            ]
        ],
        on="program_uid",
        how="left",
        validate="many_to_one",
    )
    require(len(lineage) == 12, "lineage plotting contract must have 12 rows")
    lineage["program_name"] = lineage["program_uid"].map(HEROES)
    lineage["claim_boundary"] = "same-atlas descriptive localization"
    write_tsv(lineage, "fig4e_same_atlas_lineage_specificity.tsv")

    # 4F. Assay-native units stay on separate aligned axes.
    sc_rows = []
    for row in effects[effects["program_uid"].isin(HEROES)].itertuples(index=False):
        for estimate_type in ("raw", "corrected"):
            sc_rows.append(
                {
                    "program_uid": row.program_uid,
                    "program_name": HEROES[row.program_uid],
                    "estimate_type": estimate_type,
                    "beta": getattr(row, f"{estimate_type}_beta"),
                    "ci_low": getattr(row, f"{estimate_type}_ci_low"),
                    "ci_high": getattr(row, f"{estimate_type}_ci_high"),
                    "qvalue": getattr(row, f"{estimate_type}_qvalue"),
                    "n_biological_donors": row.n_model_donors,
                    "unit": "Hotspot score per cross-sectional stage ordinal",
                    "analysis_universe": row.analysis_universe,
                }
            )
    write_tsv(pd.DataFrame(sc_rows), "fig4f_singlecell_estimates.tsv")
    bulk = pd.read_csv(BULK, sep="\t")
    require(len(bulk) == 8 and set(bulk["program_uid"]) == set(HEROES), "bulk projection drift")
    bulk["program_name"] = bulk["program_uid"].map(HEROES)
    bulk["unit"] = "frozen L1-weighted member-gene log2FC, fibrosis stage versus F0"
    bulk["claim_boundary"] = "tissue-state transport; not a program-level significance test"
    write_tsv(bulk, "fig4f_bulk_stage_transport.tsv")

    # 5F. Existing maps remain upstream; only the presentation contract changes.
    spatial = pd.read_csv(SPATIAL, sep="\t")
    require(len(spatial) == 4 and set(spatial["program_id"]) == set(HEROES), "spatial family drift")
    spatial["map_panel_path"] = str(SPATIAL_MAP.resolve())
    spatial["map_panel_sha256"] = sha256(Path(spatial["map_panel_path"].iloc[0]))
    spatial["claim_boundary_rebuild"] = (
        "matched-null tissue organization; not lineage origin and not a spatial disease contrast"
    )
    write_tsv(spatial, "fig5f_spatial_map_and_matched_null.tsv")

    # 6C. Detailed fingerprints, followed by one display node per assay-native claim.
    fingerprints: list[dict[str, object]] = []
    triad_example = "GNMT–MAT1A–CYP2C19"
    triad_next = "Allele-aware hepatocyte perturbation at the promoted GNMT shared signal"
    fingerprints.append(
        evidence_row(
            triad_example, 1, "Inherited shared signal", "genetics", "GNMT",
            "corrected COLOC posterior", None, "PP.H4 for an exact eQTL–GWAS signal pair",
            "blocked_pending_corrected_coloc",
            "The corrected 50-study COLOC release is active and not promoted.",
            BLOCKERS,
            "The prior candidate signal may change after corrected signal-pair adjudication.",
            triad_next, "Corrected COLOC pending",
        )
    )
    cohort = pd.read_csv(COHORT, sep="\t")
    cohort = cohort[cohort["gene_name"].isin(TRIAD)].copy()
    require(len(cohort) == 15, "triad five-cohort rows must be 3 x 5")
    for gene, group in cohort.groupby("gene_name", sort=False):
        n_down = int((group["logFC"] < 0).sum())
        n_q = int((group["FDR"] < 0.05).sum())
        for row in group.itertuples(index=False):
            fingerprints.append(
                evidence_row(
                    triad_example, 2, "Five-cohort RNA remodeling", "bulk RNA-seq", gene,
                    "cohort disease-versus-control logFC", float(row.logFC), "cohort-specific log2FC",
                    "supported_multicohort_remodeling" if n_down >= 4 else "heterogeneous",
                    f"{n_down}/5 cohort effects are negative; {n_q}/5 have within-cohort gene-family BH q<0.05.",
                    COHORT, "Cell composition and disease severity can contribute to bulk RNA effects.",
                    "Cell-type-restricted perturbation with RNA and metabolic-flux readouts",
                    f"{gene}: {n_down}/5 down; {n_q}/5 q<0.05", f"dataset={row.dataset}; q={row.FDR:.3g}",
                )
            )
    stage = pd.read_csv(STAGE, sep="\t")
    stage = stage[stage["gene_name"].isin(TRIAD)].copy()
    require(len(stage) == 12, "triad adjacent-stage rows must be 3 x 4")
    for gene, group in stage.groupby("gene_name", sort=False):
        n_down = int((group["logFC"] < 0).sum())
        n_q = int((group["FDR"] < 0.05).sum())
        for row in group.itertuples(index=False):
            fingerprints.append(
                evidence_row(
                    triad_example, 3, "Cross-sectional histologic stage", "bulk RNA-seq", gene,
                    "adjacent fibrosis-stage logFC", float(row.logFC), "log2FC for adjacent fibrosis stages",
                    "supported_stage_association" if n_q else "directional_only",
                    f"{n_down}/4 adjacent-stage effects are negative; {n_q}/4 pass BH over 23,370 genes.",
                    STAGE, "Adjacent-stage contrasts are cross-sectional, not longitudinal progression.",
                    "Prospective or longitudinal sampling with repeated histology",
                    f"{gene}: {n_down}/4 down; {n_q}/4 q<0.05", f"contrast={row.contrast}; q={row.FDR:.3g}",
                )
            )
    protein = pd.read_csv(PROTEIN)
    protein = protein[protein["gene"].isin(TRIAD)].copy()
    require(len(protein) == 3, "triad protein rows missing")
    for row in protein.itertuples(index=False):
        fingerprints.append(
            evidence_row(
                triad_example, 4, "Liver protein decrease", "liver DIA-MS", row.gene,
                "adjusted MASLD-versus-control protein effect", float(row.protein_logFC), "protein log2FC",
                "supported_selection_conditioned"
                if str(row.protein_sig).strip().lower() == "true"
                else "not_supported",
                f"Adjusted protein q={row.protein_padj:.3g}; same-cohort selection-conditioned display.",
                PROTEIN, "Selection in the same protein cohort precludes independent validation.",
                "Independent liver proteomics with a frozen triad test",
                f"{row.gene}: β={row.protein_logFC:.2f}; q={row.protein_padj:.3g}",
            )
        )
    hist = pd.read_csv(PROTEIN_HIST)
    hist = hist[hist["gene"].isin(TRIAD)].copy()
    require(len(hist) == 3, "triad protein-histology rows missing")
    for row in hist.itertuples(index=False):
        for feature in ("Steatosis", "Ballooning", "Inflammation", "Fibrosis", "NAS"):
            rho = float(getattr(row, feature))
            fingerprints.append(
                evidence_row(
                    triad_example, 5, "Adjusted protein covariation", "liver DIA-MS", row.gene,
                    f"partial Spearman association with {feature}", rho,
                    "partial Spearman rho after batch, age, BMI, and sex residualization",
                    "descriptive_selection_conditioned",
                    "The displayed table contains effect estimates; multiplicity-adjusted test results are not exported here.",
                    PROTEIN_HIST, "Covariation does not establish direction or molecular mediation.",
                    "Perturb the triad and measure methionine-cycle flux plus histology-linked phenotypes",
                    f"{row.gene}: ρ {min(row[1:]):.2f} to {max(row[1:]):.2f}", f"feature={feature}",
                )
            )

    bulk_by_uid = bulk.groupby("program_uid")
    spatial_by_uid = spatial.groupby("program_id")
    for uid, example in HEROES.items():
        hero_effect = effects[effects["program_uid"] == uid]
        require(len(hero_effect) == 1, f"hero effect missing: {uid}")
        row = hero_effect.iloc[0]
        next_experiment = (
            "Fibroblast-restricted IGFBP7 perturbation with matrix and hepatocyte-state readouts"
            if "IGFBP7" in example
            else "Cholangiocyte-restricted BICC1 perturbation in a ductular-injury tissue model"
        )
        fingerprints.append(
            evidence_row(
                example, 1, "Donor-level disease association", "single-cell RNA", uid,
                "raw frozen-program stage coefficient", float(row.raw_beta),
                "Hotspot score per cross-sectional stage ordinal", "supported" if row.raw_qvalue < 0.05 else "not_supported",
                f"BH q={row.raw_qvalue:.3g} across 117 programs; n={int(row.n_model_donors)} biological donors.",
                AMBIENT / "results/ambient_program_effects.tsv",
                "Association can reflect transcript redistribution between cells.", next_experiment,
                f"raw β={row.raw_beta:.2f}; q={row.raw_qvalue:.3g}",
            )
        )
        fingerprints.append(
            evidence_row(
                example, 2, "Ambient recalibration", "single-cell RNA", uid,
                "corrected frozen-program stage coefficient", float(row.corrected_beta),
                "Hotspot score per cross-sectional stage ordinal", str(row.evidence_state),
                f"Corrected BH q={row.corrected_qvalue:.3g}; paired change β={row.delta_beta:.3g}, q={row.delta_qvalue:.3g}.",
                AMBIENT / "results/ambient_program_effects.tsv",
                "Ambient correction separates sensitivity from origin but cannot establish cell autonomy.", next_experiment,
                f"corrected β={row.corrected_beta:.2f}; {row.evidence_state}",
            )
        )
        promo = promotion[promotion["program_uid"] == uid].iloc[0]
        target = promo.target_lineage
        contrast = primary_contrasts[
            (primary_contrasts["program_uid"] == uid)
            & (primary_contrasts["comparison_lineage"] == target)
        ]
        require(len(contrast) == 1, f"target contrast missing: {uid}")
        contrast = contrast.iloc[0]
        lineage_state = (
            "descriptive_localization_supported"
            if str(promo.main_figure_family_eligible).strip().lower() == "true"
            else "not_promoted_family_gate"
        )
        fingerprints.append(
            evidence_row(
                example, 3, "Same-atlas localization", "single-cell RNA", uid,
                f"{target}-minus-hepatocyte fixed-weight contrast", float(contrast.beta),
                "within-dataset standardized donor-paired program score difference", lineage_state,
                f"BH q={contrast.qvalue:.3g} in the ten-comparison family; destination={promo.destination}.",
                CROSS / "results/hero_lineage_contrasts.tsv",
                "Same-atlas localization is descriptive and does not prove transcript origin.", next_experiment,
                f"{target}−hep β={contrast.beta:.2f}; q={contrast.qvalue:.3g}",
            )
        )
        hero_bulk = bulk_by_uid.get_group(uid)
        for bulk_row in hero_bulk.itertuples(index=False):
            fingerprints.append(
                evidence_row(
                    example, 4, "Bulk tissue-state transport", "bulk RNA-seq", uid,
                    f"fixed member-gene projection {bulk_row.stage}-versus-F0", float(bulk_row.effect),
                    "frozen L1-weighted member-gene log2FC", "directional_transport",
                    f"Direction agreement={bulk_row.direction_agreement:.2f}; this is not a program-level significance test.",
                    BULK, "Bulk transport mixes cell abundance and within-cell state.", next_experiment,
                    f"F1–F4 effect {hero_bulk.effect.min():.2f} to {hero_bulk.effect.max():.2f}",
                )
            )
        hero_spatial = spatial_by_uid.get_group(uid)
        for spatial_row in hero_spatial.itertuples(index=False):
            fingerprints.append(
                evidence_row(
                    example, 5, "Matched-null tissue organization", "spatial transcriptomics", uid,
                    "residual Moran I", float(spatial_row.residual_moran_i), "residual Moran I",
                    str(spatial_row.primary_evidence_state),
                    f"{spatial_row.dataset_label}: matched-null q={spatial_row.primary_qvalue_spatial:.3g}; {spatial_row.unit_sign_label}.",
                    SPATIAL,
                    "Spatial organization is not lineage origin or a spatial disease-versus-healthy contrast.",
                    next_experiment,
                    f"{spatial_row.dataset_label}: I={spatial_row.residual_moran_i:.2f}; {spatial_row.display_state}",
                )
            )

    fingerprints_frame = pd.DataFrame(fingerprints)
    write_tsv(fingerprints_frame, "fig6c_evidence_fingerprints.tsv")

    # Display nodes summarize only within an assay node; they never count votes or rank examples.
    node_rows = []
    for (example, order, node), group in fingerprints_frame.groupby(
        ["example", "node_order", "node"], sort=False
    ):
        states = list(dict.fromkeys(group["state"].astype(str)))
        displays = list(dict.fromkeys(group["short_display"].astype(str)))
        if len(displays) > 3:
            display = "; ".join(displays[:3])
        else:
            display = "; ".join(displays)
        state = states[0] if len(states) == 1 else "assay_native_mixed_states"
        node_rows.append(
            {
                "example": example,
                "node_order": order,
                "node": node,
                "state": state,
                "display": display,
                "native_units": "; ".join(dict.fromkeys(group["unit"].astype(str))),
                "unresolved_alternative": group["unresolved_alternative"].iloc[0],
                "next_experiment": group["next_experiment"].iloc[0],
                "n_native_rows": len(group),
            }
        )
    nodes6 = pd.DataFrame(node_rows)
    require(len(nodes6) == 15, f"expected 15 evidence nodes, found {len(nodes6)}")
    write_tsv(nodes6, "fig6c_evidence_nodes.tsv")
    experiment_router = fingerprints_frame[
        ["example", "node", "unresolved_alternative", "next_experiment"]
    ].drop_duplicates().sort_values(["example", "node"])
    require(
        experiment_router["unresolved_alternative"].str.len().gt(0).all()
        and experiment_router["next_experiment"].str.len().gt(0).all(),
        "every unresolved alternative must map directly to an experiment",
    )
    write_tsv(experiment_router, "fig6d_experiment_router.tsv")

    claims = pd.DataFrame(
        [
            (
                "Figure 4",
                "Donor-associated liver programs differ in how much their disease effect survives ambient correction, separating disease association from cellular origin.",
                "candidate_pending_synchronized_promotion",
            ),
            (
                "Figure 5",
                "Disease association and tissue organization are different claims, as shown by supported matched-null organization for ECM/IGFBP7 and an indeterminate state for BICC1.",
                "candidate_pending_synchronized_promotion",
            ),
            (
                "Figure 6",
                "Recurring examples connect assay-native evidence states to discriminating experiments without collapsing modalities into a score or rank.",
                "candidate_pending_corrected_coloc_and_catalog_rebuild",
            ),
        ],
        columns=["figure", "sentence_level_discovery", "release_state"],
    )
    write_tsv(claims, "figure_level_claims.tsv")

    # Assembly routing makes retained and supplementary panels explicit.
    routing = pd.DataFrame(
        [
            ("4A", "new", "panels/fig4a_analysis_logic.pdf", "editable vector schematic"),
            ("4B", "retain", str((STAGE_ROOT / "figure4/panels/fig4b_scrna_umap_embeddable.pdf").resolve()), "directly labelled UMAP and donor density"),
            ("4C", "new", "panels/fig4c_all117_disease_skyline.pdf", "all 117 programs"),
            ("4D", "new", "panels/fig4d_ambient_effect_transport.pdf", "GSE189600 visible"),
            ("4E", "conditional", f"panels/{'fig4e' if family_destination[0] == 'Figure_4E' else 'figS4e'}_same_atlas_lineage_specificity.pdf", family_destination[0]),
            ("4F", "new", "panels/fig4f_tissue_state_transport.pdf", "assay-native aligned axes"),
            ("5F", "rebuild+retain", "panels/fig5f_spatial_maps_and_matched_null.pdf", "representative tissue maps paired with observed residual Moran I versus matched null"),
            ("6C", "new", "panels/fig6c_evidence_ribbons.pdf", "no combined score or rank"),
            ("6D", "source-rebuild", "source_tables/fig6d_experiment_router.tsv", "each unresolved alternative maps directly to a discriminating experiment"),
            ("S4-program-detail", "move", "supplementary/figureS4/panels/figs4a_program_landscape_detail.pdf", "current detailed heatmaps"),
            ("S4-communication", "move", "supplementary/figureS4/panels/figs4b_fixed_communication.pdf", "13-pair roster"),
            ("S4-TF", "move", "supplementary/figureS4/panels/figs4c_tf_activity.pdf", "20-TF display"),
            ("S5-source-state", "move", "supplementary/figureS5/panels/figs5a_multimodal_program_source_matrix.pdf", "complete source-state matrix"),
            ("S5-composition", "move", "supplementary/figureS5/panels/figs5b_full_composition_sensitivity.pdf", "complete composition sensitivity"),
        ],
        columns=["panel", "action", "artifact", "reason"],
    )
    write_tsv(routing, "panel_routing.tsv")

    caption_dir = OUT / "captions"
    caption_dir.mkdir(parents=True, exist_ok=True)
    caption4 = caption_dir / "figure4.md"
    caption5 = caption_dir / "figure5.md"
    caption6 = caption_dir / "figure6.md"
    for path in (caption4, caption5, caption6):
        refuse(path)
    caption4.write_text(
        "Figure 4. Ambient correction separates disease association from cellular origin. "
        "A, fixed-program analysis logic. B, retained directly labelled atlas UMAP and donor-level disease density. "
        "C, disease effects for all 117 frozen programs; direct labels are restricted to the two recurring examples and four prespecified controls. "
        "D, complete-case raw versus ambient-corrected effects; GSE189600 is shown explicitly as an uncorrected passthrough only in the full-universe sensitivity. "
        "The 11 frozen T-cell programs are displayed as untestable because the complete-atlas integer transport does not reproduce their native non-integer scoring substrate. "
        "E, same-atlas descriptive localization for ECM/IGFBP7 and ductular-injury/BICC1; donors are paired between each lineage and hepatocytes, and this panel does not establish transcript origin or cell autonomy. "
        "F, raw and corrected single-cell estimates aligned with fixed member-gene bulk fibrosis-stage transport on separate assay-native axes. "
        "Ordinary and HC3 uncertainty are retained in the source tables; BH correction is applied across the stated 117-program or ten-comparison family.\n",
        encoding="utf-8",
    )
    caption5.write_text(
        "Figure 5F. Disease association and tissue organization are different claims. Representative tissue maps are paired with observed residual Moran I and the corresponding matched-gene null distribution. ECM/IGFBP7 exceeds matched spatial background, whereas BICC1 remains indeterminate. These statistics describe tissue organization, not lineage origin and not a spatial MASLD-versus-healthy contrast.\n",
        encoding="utf-8",
    )
    caption6.write_text(
        "Figure 6C. Assay-native evidence ribbons for GNMT–MAT1A–CYP2C19, ECM/IGFBP7, and ductular-injury/BICC1. Every node retains its native effect unit, provenance, evidence state, unresolved alternative, and discriminating experiment. No combined score, vote count, or rank is calculated. The GNMT inherited shared-signal node is withheld until the corrected 50-study COLOC release is audited and adopted.\n",
        encoding="utf-8",
    )

    manifest_rows = []
    for path in inputs:
        manifest_rows.append({"role": "input", **provenance(path), "bytes": path.stat().st_size})
    input_manifest = OUT / "input_manifest.tsv"
    refuse(input_manifest)
    pd.DataFrame(manifest_rows).to_csv(input_manifest, sep="\t", index=False)
    build_meta = OUT / "build_metadata.json"
    refuse(build_meta)
    build_meta.write_text(
        json.dumps(
            {
                "candidate_id": OUT.name,
                "ambient_candidate": str(AMBIENT.resolve()),
                "ambient_validated_sha256": sha256(AMBIENT / "VALIDATED"),
                "cross_lineage_candidate": str(CROSS.resolve()),
                "cross_lineage_validated_sha256": sha256(CROSS / "VALIDATED"),
                "lineage_panel_destination": family_destination[0],
                "genetics_state": "blocked_pending_corrected_coloc",
                "combined_score_calculated": False,
                "seed": 20260815,
            },
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    print(f"Built {len(list(SOURCE.glob('*.tsv')))} source tables in {SOURCE}")


if __name__ == "__main__":
    main()
