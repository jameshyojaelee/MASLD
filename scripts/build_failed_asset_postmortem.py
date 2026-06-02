"""Build failed_asset_postmortem.csv joining historical Phase 2/3 MASLD
trial data with atlas evidence (bulk LFC, per-cell-type LFC,
cascade-stage class, zonation, COLOC PP4 PolyFun-canonical).

RP3a Top-2 deliverable from 2026-05-12 deliberation. Anchored on the
audit finding that NONE of the failed targets (cenicriviroc CCR2/5,
simtuzumab LOXL2, selonsertib MAP3K5, aramchol SCD, OCA NR1H4,
elafibranor PPARA/D) clear PolyFun PP4>0.5 — vs THRB anchor 1.00 —
making "the target was never genetically real" the atlas's
strongest in-licensing deliverable.

Output: RNA-seq/results/drug_repurposing/failed_asset_postmortem.csv

Run under rnaseq env:
    micromamba activate rnaseq
    python scripts/build_failed_asset_postmortem.py
"""

from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd


PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
ATLAS_PATH = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
CELLTYPE_PATH = PROJECT_ROOT / "RNA-seq/results/celltype_attribution/celltype_attribution_matrix.csv"
PROGRESSION_PATH = PROJECT_ROOT / "RNA-seq/results/stratified_causal/onset_vs_progression_genes.csv"
ZONE_PATH = PROJECT_ROOT / "RNA-seq/results/stratified_causal/zone_causal_scores.csv"
OUT_PATH = PROJECT_ROOT / "RNA-seq/results/drug_repurposing/failed_asset_postmortem.csv"


# Curated Phase 2/3 MASLD drug roster spanning failures + active assets +
# the FDA-approved positive control (Resmetirom). Each row = one
# (drug, target_gene) pair so multi-target compounds have one row per
# target.
#
# Failure mode tags:
#   - failed_p3_efficacy : Phase 3 missed primary endpoint
#   - failed_p3_terminated : Phase 3 sponsor-terminated for futility
#   - failed_p2_efficacy : Phase 2 missed primary endpoint
#   - fda_rejected      : FDA CRL after Phase 3
#   - approved_p3        : FDA-approved (positive control)
#   - active_p2 / active_p3 : ongoing
#
# Extrahepatic flag = TRUE when the target's known therapeutic
# mechanism is non-hepatic (gut/CNS/adipose/systemic-metabolic).
DRUG_ROSTER = [
    # ---- failed Phase 2/3 (target post-mortems) ----
    dict(drug="Cenicriviroc", target_gene="CCR2",
         nct_id="NCT03028740 AURORA / NCT02217475 CENTAUR",
         sponsor="Allergan / AbbVie",
         stage="Phase 3 failed (2020)",
         moa="CCR2/CCR5 dual antagonist",
         failure_mode="failed_p3_efficacy",
         extrahepatic_mechanism=False,
         notes="AURORA Phase 3 missed fibrosis-improvement endpoint at week 12 (2020). CENTAUR Phase 2 was promising but did not replicate."),
    dict(drug="Cenicriviroc", target_gene="CCR5",
         nct_id="NCT03028740 AURORA / NCT02217475 CENTAUR",
         sponsor="Allergan / AbbVie",
         stage="Phase 3 failed (2020)",
         moa="CCR2/CCR5 dual antagonist",
         failure_mode="failed_p3_efficacy",
         extrahepatic_mechanism=False,
         notes="Second target of cenicriviroc."),
    dict(drug="Simtuzumab", target_gene="LOXL2",
         nct_id="NCT01672879 / NCT01672866",
         sponsor="Gilead",
         stage="Phase 2 failed (2016-2017)",
         moa="Humanized anti-LOXL2 mAb",
         failure_mode="failed_p2_efficacy",
         extrahepatic_mechanism=False,
         notes="Two Phase 2 trials (NASH-fibrosis + PSC) terminated 2016-2017 for futility on histologic endpoints."),
    dict(drug="Selonsertib", target_gene="MAP3K5",
         nct_id="NCT03053050 STELLAR-3 / NCT03053063 STELLAR-4",
         sponsor="Gilead",
         stage="Phase 3 failed (2019)",
         moa="ASK1 inhibitor (MAP3K5)",
         failure_mode="failed_p3_efficacy",
         extrahepatic_mechanism=False,
         notes="STELLAR-3 (F3) and STELLAR-4 (F4) both missed primary fibrosis-improvement endpoint (2019)."),
    dict(drug="Aramchol", target_gene="SCD",
         nct_id="NCT04104321 ARMOR",
         sponsor="Galmed",
         stage="Phase 3 terminated",
         moa="SCD1 inhibitor (fatty-acid-bile-acid conjugate)",
         failure_mode="failed_p3_terminated",
         extrahepatic_mechanism=False,
         notes="ARMOR Phase 3 paused 2022; conditional resumption stalled. Phase 2 NCT02684591 missed primary."),
    dict(drug="Obeticholic_acid", target_gene="NR1H4",
         nct_id="NCT03439254 REGENERATE",
         sponsor="Intercept",
         stage="Phase 3 FDA rejected (2023)",
         moa="FXR agonist",
         failure_mode="fda_rejected",
         extrahepatic_mechanism=False,
         notes="REGENERATE Phase 3 N=919 met histology at week 18 but FDA issued CRL in 2023 over benefit-risk. Pruritus dose-limiting."),
    dict(drug="Elafibranor", target_gene="PPARA",
         nct_id="NCT02704403 RESOLVE-IT",
         sponsor="Genfit / Ipsen",
         stage="Phase 3 terminated (2020)",
         moa="PPAR alpha/delta agonist",
         failure_mode="failed_p3_terminated",
         extrahepatic_mechanism=False,
         notes="RESOLVE-IT Phase 3 NASH terminated 2020 for futility. Approved 2024 for PBC."),
    dict(drug="Elafibranor", target_gene="PPARD",
         nct_id="NCT02704403 RESOLVE-IT",
         sponsor="Genfit / Ipsen",
         stage="Phase 3 terminated (2020)",
         moa="PPAR alpha/delta agonist",
         failure_mode="failed_p3_terminated",
         extrahepatic_mechanism=False,
         notes="Second target of elafibranor."),
    # ---- positive control: FDA approved ----
    dict(drug="Resmetirom", target_gene="THRB",
         nct_id="NCT03900429 MAESTRO-NASH / NCT05500222 MAESTRO-NASH-OUTCOMES",
         sponsor="Madrigal",
         stage="FDA approved (Mar 2024)",
         moa="THR-beta selective agonist",
         failure_mode="approved_p3",
         extrahepatic_mechanism=False,
         notes="MAESTRO-NASH met both NAS resolution and fibrosis improvement primary endpoints. First FDA-approved MASH drug."),
    # ---- active Phase 3 (ongoing) ----
    dict(drug="Lanifibranor", target_gene="PPARA",
         nct_id="NCT04849728 NATiV3",
         sponsor="Inventiva",
         stage="Phase 3 ongoing",
         moa="Pan-PPAR alpha/delta/gamma agonist",
         failure_mode="active_p3",
         extrahepatic_mechanism=False,
         notes="NATiV3 Phase 3 active. Phase 2 NATIVE met primary 2020."),
    dict(drug="Lanifibranor", target_gene="PPARD",
         nct_id="NCT04849728 NATiV3",
         sponsor="Inventiva",
         stage="Phase 3 ongoing",
         moa="Pan-PPAR alpha/delta/gamma agonist",
         failure_mode="active_p3",
         extrahepatic_mechanism=False,
         notes=""),
    dict(drug="Lanifibranor", target_gene="PPARG",
         nct_id="NCT04849728 NATiV3",
         sponsor="Inventiva",
         stage="Phase 3 ongoing",
         moa="Pan-PPAR alpha/delta/gamma agonist",
         failure_mode="active_p3",
         extrahepatic_mechanism=False,
         notes=""),
    dict(drug="Efruxifermin", target_gene="FGFR1",
         nct_id="NCT06215430 SYMMETRY",
         sponsor="Akero",
         stage="Phase 3 ongoing",
         moa="FGF21 analog (fusion Fc)",
         failure_mode="active_p3",
         extrahepatic_mechanism=True,
         notes="EXTRAHEPATIC mechanism: FGFR1+KLB co-receptor effect is largely adipose / systemic-metabolic. Atlas hepatic-intrinsic absence is mechanistically consistent with positive Phase 2 efficacy (HARMONY) via non-hepatocyte axis."),
    dict(drug="Efruxifermin", target_gene="KLB",
         nct_id="NCT06215430 SYMMETRY",
         sponsor="Akero",
         stage="Phase 3 ongoing",
         moa="FGF21 analog",
         failure_mode="active_p3",
         extrahepatic_mechanism=True,
         notes="EXTRAHEPATIC co-receptor."),
    dict(drug="Pegozafermin", target_gene="FGFR1",
         nct_id="NCT06318169 ENLIGHTEN-Fibrosis",
         sponsor="89bio",
         stage="Phase 3 ongoing",
         moa="FGF21 analog (pegylated)",
         failure_mode="active_p3",
         extrahepatic_mechanism=True,
         notes="ENLIVEN Phase 2b met endpoints (2023). Same caveat as efruxifermin."),
    dict(drug="Survodutide", target_gene="GLP1R",
         nct_id="NCT06251544",
         sponsor="Boehringer",
         stage="Phase 3 ongoing",
         moa="Dual GLP-1 / glucagon receptor agonist",
         failure_mode="active_p3",
         extrahepatic_mechanism=True,
         notes="EXTRAHEPATIC: GLP1R is gut/CNS/pancreas. Effect on MASH likely weight-mediated."),
    dict(drug="Survodutide", target_gene="GCGR",
         nct_id="NCT06251544",
         sponsor="Boehringer",
         stage="Phase 3 ongoing",
         moa="Dual GLP-1 / glucagon receptor agonist",
         failure_mode="active_p3",
         extrahepatic_mechanism=True,
         notes="Glucagon receptor — hepatic but lipogenic-regulation context, not direct anti-fibrotic mechanism."),
    dict(drug="Semaglutide", target_gene="GLP1R",
         nct_id="NCT04822181 ESSENCE",
         sponsor="Novo Nordisk",
         stage="FDA approved Aug 2025 (MASH expansion)",
         moa="GLP-1 receptor agonist",
         failure_mode="approved_p3",
         extrahepatic_mechanism=True,
         notes="EXTRAHEPATIC. Approved 2025 for MASH via weight-loss-mediated mechanism. Atlas hepatic-intrinsic absence is consistent."),
    dict(drug="Tirzepatide", target_gene="GLP1R",
         nct_id="NCT05205577 SYNERGY-NASH",
         sponsor="Eli Lilly",
         stage="Phase 3 ongoing",
         moa="GLP-1 + GIP dual agonist",
         failure_mode="active_p3",
         extrahepatic_mechanism=True,
         notes=""),
    dict(drug="Tirzepatide", target_gene="GIPR",
         nct_id="NCT05205577 SYNERGY-NASH",
         sponsor="Eli Lilly",
         stage="Phase 3 ongoing",
         moa="GLP-1 + GIP dual agonist",
         failure_mode="active_p3",
         extrahepatic_mechanism=True,
         notes="GIPR mainly adipose."),
    # ---- Phase 2 active ----
    dict(drug="ION224", target_gene="DGAT2",
         nct_id="NCT04483947",
         sponsor="Ionis",
         stage="Phase 2b ongoing",
         moa="GalNAc ASO (DGAT2 knockdown)",
         failure_mode="active_p2",
         extrahepatic_mechanism=False,
         notes=""),
    dict(drug="Denifanstat", target_gene="FASN",
         nct_id="NCT06129084 FASCINATE-3",
         sponsor="Sagimet",
         stage="Phase 3 planned",
         moa="FASN inhibitor",
         failure_mode="active_p2",
         extrahepatic_mechanism=False,
         notes=""),
    dict(drug="AZD2693", target_gene="PNPLA3",
         nct_id="NCT05809934",
         sponsor="AstraZeneca",
         stage="Phase 2b completed (data 2025)",
         moa="GalNAc ASO (PNPLA3 I148M carriers)",
         failure_mode="active_p2",
         extrahepatic_mechanism=False,
         notes="Genotype-stratified (I148M G/G)."),
    dict(drug="Rapirosiran", target_gene="HSD17B13",
         nct_id="NCT06122337",
         sponsor="Arrowhead / Pfizer",
         stage="Phase 3 planned",
         moa="GalNAc-siRNA (HSD17B13 silencing)",
         failure_mode="active_p2",
         extrahepatic_mechanism=False,
         notes=""),
]


# ---- atlas columns to extract ----
ATLAS_COLS = [
    "human_symbol", "ensembl_id", "gene_biotype",
    "dream_logFC", "dream_padj", "dream_tstat",
    "mouse_meta_logFC", "mouse_meta_padj", "n_diets_sig",
    "primary_category", "is_conserved",
    "coloc_best_pp4_polyfun", "coloc_best_susie_pp4_polyfun",
    "coloc_best_gwas_polyfun", "coloc_best_susie_gwas_polyfun",
    "coloc_susie_best_pp4", "coloc_susie_best_pp4_v1_sghatan",
    "coloc_abf_best_pp4", "coloc_abf_best_gwas",
    "broadaway_coloc_pp4",
    "ggt_coloc_pp4", "ast_coloc_pp4", "pdff_coloc_pp4", "ukbb_alt_coloc_pp4",
    "bbj_ggt_coloc_pp4",
    "finngen_nafld_coloc_pp4",
    "n_coloc_sources", "n_ancestry_gwas",
    "twas_z", "twas_pval",
    "essentiality_chronos", "is_essential",
    "attribution_class",
    "dgidb_druggable", "lincs_reversal",
    "f2_inflection_logFC", "f2_inflection_padj",
    "fibrosis_ordinal_coef", "inflammation_ordinal_coef",
    "sex_class", "sex_interaction_padj",
    "ferroptosis_class", "zonation_class",
    "is_onset_specific", "is_progression_specific",
]

# Per-CT LFC + padj columns we keep from celltype_attribution_matrix.
CT_NAMES = [
    "Hepatocytes", "Cholangiocytes", "Macrophages",
    "Mono+mono_derived_cells", "Endothelial_cells", "Fibroblasts",
    "B_cells", "T_cells", "Plasma_cells",
    "Resident_NK", "Circulating_NK_NKT",
]


def _to_float(v, default=0.0) -> float:
    """Coerce CSV/NaN/blank to float."""
    try:
        f = float(v)
        return default if np.isnan(f) else f
    except (TypeError, ValueError):
        return default


def verdict_from_evidence(row: pd.Series) -> tuple[str, str]:
    """Compute (target_real_verdict, rationale) per (drug, target_gene)."""
    pp4 = max(
        _to_float(row.get("coloc_best_susie_pp4_polyfun")),
        _to_float(row.get("coloc_best_pp4_polyfun")),
        _to_float(row.get("coloc_susie_best_pp4")),
        _to_float(row.get("coloc_abf_best_pp4")),
    )
    bulk_lfc = _to_float(row.get("dream_logFC"))
    bulk_padj = _to_float(row.get("dream_padj"), default=1.0)
    is_deg = (bulk_padj < 0.05) and (abs(bulk_lfc) > 0.5)
    extra = bool(row.get("extrahepatic_mechanism", False))

    if extra:
        if pp4 < 0.5 and not is_deg:
            return ("EXTRAHEPATIC_CONSISTENT",
                    f"Extrahepatic mechanism (GLP1R/FGFR1/GCGR/KLB/GIPR-class). "
                    f"Atlas hepatic-intrinsic absence (PP4={pp4:.3f}, DEG={is_deg}) "
                    f"is mechanistically consistent with non-hepatocyte axis.")
        return ("EXTRAHEPATIC_PARTIAL",
                f"Extrahepatic-class target with detectable hepatic signal "
                f"(PP4={pp4:.3f}, DEG={is_deg}).")

    if pp4 >= 0.9:
        return ("TARGET_REAL_STRONG",
                f"Strong COLOC support (PP4={pp4:.3f}). Trial enrolment / "
                f"dose / endpoint are likely failure modes; mechanism is "
                f"genetically validated.")
    if pp4 >= 0.5:
        return ("TARGET_REAL_MODERATE",
                f"Moderate COLOC support (PP4={pp4:.3f}). Likely "
                f"enrolment / stage-window or dose issue.")
    if pp4 >= 0.1 and is_deg:
        return ("TARGET_PROXIMITY",
                f"Bulk DEG without COLOC support (PP4={pp4:.3f}, "
                f"logFC={bulk_lfc:.2f}, padj={bulk_padj:.1e}). "
                f"Likely co-expression signature, not causal gene.")
    if is_deg:
        return ("TARGET_NEVER_REAL_REACTIVE",
                f"Bulk DEG (logFC={bulk_lfc:.2f}, padj={bulk_padj:.1e}) "
                f"but PP4={pp4:.3f}<<0.5. Target was likely a reactive "
                f"marker, NOT a causal driver. Failure mode = target dead.")
    return ("TARGET_NEVER_REAL_SILENT",
            f"No bulk DEG and PP4={pp4:.3f}. Target was silent at both "
            f"layers. Failure mode = target dead.")


def main() -> None:
    print(f"Reading atlas {ATLAS_PATH}")
    atlas = pd.read_csv(ATLAS_PATH, low_memory=False, usecols=ATLAS_COLS)
    atlas_index = atlas.set_index("human_symbol")

    print(f"Reading celltype attribution {CELLTYPE_PATH}")
    ct = pd.read_csv(CELLTYPE_PATH, low_memory=False)
    ct_cols = ["symbol", "bulk_lfc", "bulk_padj"]
    for cn in CT_NAMES:
        ct_cols += [f"lfc_{cn}", f"padj_{cn}"]
    ct = ct[[c for c in ct_cols if c in ct.columns]]
    ct = ct.rename(columns={"symbol": "human_symbol"}).set_index("human_symbol")

    print(f"Reading progression class {PROGRESSION_PATH}")
    prog = pd.read_csv(PROGRESSION_PATH, low_memory=False,
                       usecols=["symbol", "stages_sig", "n_stages_sig",
                                "earliest_stage", "latest_stage",
                                "progression_class", "progression_causal_score",
                                "best_stage", "best_logfc"])
    prog = prog.rename(columns={"symbol": "human_symbol",
                                "stages_sig": "progression_stages_sig",
                                "n_stages_sig": "progression_n_stages_sig",
                                "best_stage": "progression_best_stage",
                                "best_logfc": "progression_best_logfc"})
    prog_index = prog.set_index("human_symbol")

    print(f"Reading zone scores {ZONE_PATH}")
    zone = pd.read_csv(ZONE_PATH, low_memory=False)
    zone = zone[["human_symbol", "zonation_class", "spatial_zonation_class",
                 "zone_direction", "zone_causal_score", "zone_coloc_class"]]
    zone = zone.rename(columns={
        "zonation_class": "zone_classification_atlas",
        "spatial_zonation_class": "zone_classification_spatial",
    }).set_index("human_symbol")

    out_rows = []
    for d in DRUG_ROSTER:
        sym = d["target_gene"]
        row = {**d}
        # Atlas joinflag = atlas_index.index == sym
        if sym in atlas_index.index:
            a = atlas_index.loc[sym]
            if isinstance(a, pd.DataFrame):  # symbol collision; take first
                a = a.iloc[0]
            for c in ATLAS_COLS:
                if c == "human_symbol":
                    continue
                row[f"atlas_{c}"] = a.get(c, np.nan)
        else:
            for c in ATLAS_COLS:
                if c == "human_symbol":
                    continue
                row[f"atlas_{c}"] = np.nan
        # Per-CT join
        if sym in ct.index:
            cc = ct.loc[sym]
            if isinstance(cc, pd.DataFrame):
                cc = cc.iloc[0]
            for c in ct.columns:
                row[f"ct_{c}"] = cc.get(c, np.nan)
        else:
            for c in ct.columns:
                row[f"ct_{c}"] = np.nan
        # Progression join
        if sym in prog_index.index:
            pp = prog_index.loc[sym]
            if isinstance(pp, pd.DataFrame):
                pp = pp.iloc[0]
            for c in prog_index.columns:
                row[f"prog_{c}"] = pp.get(c, np.nan)
        else:
            for c in prog_index.columns:
                row[f"prog_{c}"] = np.nan
        # Zone join
        if sym in zone.index:
            zz = zone.loc[sym]
            if isinstance(zz, pd.DataFrame):
                zz = zz.iloc[0]
            for c in zone.columns:
                row[f"zone_{c}"] = zz.get(c, np.nan)
        else:
            for c in zone.columns:
                row[f"zone_{c}"] = np.nan

        # Verdict
        verdict, rationale = verdict_from_evidence({
            **row,
            "coloc_best_susie_pp4_polyfun": row.get("atlas_coloc_best_susie_pp4_polyfun", np.nan),
            "coloc_best_pp4_polyfun": row.get("atlas_coloc_best_pp4_polyfun", np.nan),
            "coloc_susie_best_pp4": row.get("atlas_coloc_susie_best_pp4", np.nan),
            "coloc_abf_best_pp4": row.get("atlas_coloc_abf_best_pp4", np.nan),
            "dream_logFC": row.get("atlas_dream_logFC", np.nan),
            "dream_padj": row.get("atlas_dream_padj", 1.0),
        })
        row["atlas_verdict"] = verdict
        row["atlas_verdict_rationale"] = rationale

        best_pp4 = max(
            _to_float(row.get("atlas_coloc_best_susie_pp4_polyfun")),
            _to_float(row.get("atlas_coloc_best_pp4_polyfun")),
            _to_float(row.get("atlas_coloc_susie_best_pp4")),
            _to_float(row.get("atlas_coloc_abf_best_pp4")),
        )
        row["max_PP4_across_methods"] = round(best_pp4, 4)
        out_rows.append(row)

    out_df = pd.DataFrame(out_rows)

    # Reorder columns: drug-level metadata, verdict, max_PP4, then atlas, ct, prog, zone
    meta_cols = ["drug", "target_gene", "nct_id", "sponsor", "stage", "moa",
                 "failure_mode", "extrahepatic_mechanism", "notes",
                 "atlas_verdict", "atlas_verdict_rationale",
                 "max_PP4_across_methods"]
    atlas_cols_out = [c for c in out_df.columns if c.startswith("atlas_") and c not in meta_cols]
    ct_cols_out = [c for c in out_df.columns if c.startswith("ct_")]
    prog_cols_out = [c for c in out_df.columns if c.startswith("prog_")]
    zone_cols_out = [c for c in out_df.columns if c.startswith("zone_")]
    other = [c for c in out_df.columns
             if c not in (meta_cols + atlas_cols_out + ct_cols_out
                          + prog_cols_out + zone_cols_out)]
    final_cols = meta_cols + atlas_cols_out + ct_cols_out + prog_cols_out + zone_cols_out + other
    out_df = out_df[final_cols]

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(OUT_PATH, index=False)
    print(f"Wrote {OUT_PATH}")
    print(f"Rows: {len(out_df)}; Cols: {len(out_df.columns)}")
    print()
    print("Verdict breakdown:")
    print(out_df["atlas_verdict"].value_counts().to_string())
    print()
    print("Failed-asset target-was-never-real verdicts:")
    failed_mask = out_df["failure_mode"].str.startswith("failed_") | (out_df["failure_mode"] == "fda_rejected")
    print(out_df.loc[failed_mask, ["drug", "target_gene", "failure_mode",
                                   "max_PP4_across_methods",
                                   "atlas_verdict"]].to_string(index=False))


if __name__ == "__main__":
    main()
