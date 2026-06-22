#!/usr/bin/env python3
"""
generate_gene_profiles.py
=========================
Generates per-gene JSON profile files for the MASLD Atlas v2 web app.
Each gene produces a compact JSON at public/data/genes/{SYMBOL}.json
containing all available evidence data.

Usage:
  micromamba run -n spatial python scripts/generate_gene_profiles.py --output-dir public/data
"""

import argparse
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)

ATLAS_CSV = PROJECT_ROOT / "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv"
PER_STUDY_DIR = (
    PROJECT_ROOT
    / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study"
)
STAGING_DIR = (
    PROJECT_ROOT
    / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/staging_classifier"
)
FIBROSIS_STAGE_CSV = STAGING_DIR / "one_vs_rest_fibrosis_dream.csv"
FIBROSIS_VS_F0_CSV = STAGING_DIR / "fibrosis_stage_trajectory_vs_F0.csv"
FIBROSIS_VS_REST_CSV = STAGING_DIR / "fibrosis_stage_trajectory_vs_rest.csv"
# F_i vs Healthy (diagnosis_harmonized == Control). Shim produced by
# scripts/portal/build_stage_trajectories.py over fibrosis_vs_healthy_dream.csv
# (from 05d_fibrosis_vs_healthy_dream.R).
FIBROSIS_VS_HEALTHY_CSV = STAGING_DIR / "fibrosis_stage_trajectory_vs_healthy.csv"
FIBROSIS_VS_HEALTHY_DREAM_CSV = STAGING_DIR / "fibrosis_vs_healthy_dream.csv"
NAS_TRAJECTORY_CSV = STAGING_DIR / "nas_trajectory.csv"
NAS_ONE_VS_REST_CSV = STAGING_DIR / "one_vs_rest_nas_dream.csv"
# NAS-vs-fixed-baseline trajectories (Phase 7D). Shim files produced by
# scripts/portal/build_stage_trajectories.py from dream output (05c_nas_vs_baseline_dream.R).
NAS_VS_NAS0_CSV = STAGING_DIR / "nas_trajectory_vs_nas0.csv"
NAS_VS_NAS0_DREAM_CSV = STAGING_DIR / "nas_vs_nas0_dream.csv"
NAS_VS_HEALTHY_CSV = STAGING_DIR / "nas_trajectory_vs_healthy.csv"
NAS_VS_HEALTHY_DREAM_CSV = STAGING_DIR / "nas_vs_healthy_dream.csv"
PER_STUDY_CONTRASTS_DIR = (
    PROJECT_ROOT
    / "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/per_study_contrasts"
)
PER_COHORT_CONTRAST_NAMES = (
    "disease_vs_control",
    "MASH_vs_control",
    "MASL_vs_control",
    "MASH_vs_MASL",
)

# Cohorts to drop from per-cohort display ONLY (retained in integrated analyses).
# PRJNA512027 has heavy control-arm batch effect from library-prep confounding
# (L0 = all normals+steatosis; S0 = all inflammation+fibrosis). Excluded from
# per-cohort display; still included in dream mega-analysis / multi-evidence atlas.
EXCLUDED_FROM_PER_COHORT = {"PRJNA512027"}
COLOC_GENE_CSV = (
    PROJECT_ROOT / "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"
)
COLOC_ALL_GWAS_CSV = (
    PROJECT_ROOT / "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"
)
DECONV_CSV = (
    PROJECT_ROOT / "RNA-seq/results/causal_inference/deconv_attribution_scores.csv"
)
SEX_CAUSAL_CSV = (
    PROJECT_ROOT / "RNA-seq/results/stratified_causal/sex_causal_scores.csv"
)
DRUG_PROGRESSION_CSV = (
    PROJECT_ROOT
    / "RNA-seq/results/stratified_causal/drug_target_progression_classification.csv"
)
CLINICAL_DRUGS_CSV = (
    PROJECT_ROOT / "RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv"
)
LINCS_RANKED_CSV = (
    PROJECT_ROOT / "RNA-seq/results/drug_repurposing/lincs_final_ranked.csv"
)
SUBTYPE_MARKERS_CSV = PROJECT_ROOT / "RNA-seq/results/subtypes/subtype_markers.csv"
PROTEOMICS_DE_CSV = PROJECT_ROOT / "Analysis/Proteomics/results/protein_differential_results_v3.csv"
PSEUDOBULK_DIR = (
    PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/pseudobulk_de"
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def safe(v):
    """Convert NaN/Inf to None; numpy scalars to Python primitives."""
    if v is None:
        return None
    if isinstance(v, (np.floating, np.integer)):
        v = v.item()
    if isinstance(v, np.bool_):
        return bool(v)
    if isinstance(v, float):
        if math.isnan(v) or math.isinf(v):
            return None
    return v


def safe_round(v, digits=4):
    """Round a value if numeric, else return None."""
    v = safe(v)
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return round(v, digits)
    return v


def strip_version(ensembl_id: str) -> str:
    """ENSG00000121410.12 -> ENSG00000121410"""
    if isinstance(ensembl_id, str) and "." in ensembl_id:
        return ensembl_id.split(".")[0]
    return ensembl_id


def coerce_r_bool(val) -> bool:
    """Convert R-style TRUE/FALSE to Python bool."""
    return str(val).strip().upper() in ("TRUE", "1", "1.0")


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def load_atlas():
    """Load atlas CSV, build symbol->row and ensembl->symbol maps."""
    print(f"Loading atlas: {ATLAS_CSV}")
    df = pd.read_csv(ATLAS_CSV, low_memory=False)
    print(f"  {df.shape[0]} genes x {df.shape[1]} columns")
    assert {"bulk_padj", "bulk_logFC"} <= set(df.columns), (
        "C2: atlas missing bulk_* — rebuild 27a"
    )
    return df


def load_per_study_de():
    """Load per-study DE files, return dict: versioned_ensembl -> list of records."""
    print(f"Loading per-study DE from {PER_STUDY_DIR}")
    frames = []
    for fp in sorted(PER_STUDY_DIR.glob("*_de_results.csv")):
        try:
            chunk = pd.read_csv(fp, usecols=["gene", "dataset", "logFC", "P.Value"])
            frames.append(chunk)
        except Exception as e:
            print(f"  WARNING: skipping {fp.name}: {e}")
    if not frames:
        print("  WARNING: no per-study DE files found")
        return {}

    df = pd.concat(frames, ignore_index=True)
    # Drop cohorts excluded from per-cohort display (still included in integrated analyses)
    pre_n = len(df)
    df = df[~df["dataset"].isin(EXCLUDED_FROM_PER_COHORT)].copy()
    if pre_n > len(df):
        print(f"  excluded {pre_n - len(df):,} rows from {EXCLUDED_FROM_PER_COHORT}")
    df["ensembl_base"] = df["gene"].apply(strip_version)
    print(f"  {len(df)} total per-study DE records, {df['dataset'].nunique()} datasets")

    # Group by base ensembl, keep list of {dataset, logfc, pval}
    result = {}
    for ens, grp in df.groupby("ensembl_base"):
        records = []
        for _, r in grp.iterrows():
            lfc = safe_round(r["logFC"], 3)
            pv = safe_round(r["P.Value"], 4)
            if lfc is not None and pv is not None:
                records.append({
                    "dataset": str(r["dataset"]),
                    "logfc": lfc,
                    "pval": pv,
                })
        if records:
            result[ens] = records
    return result


def load_fibrosis_stage():
    """Load one-vs-rest fibrosis dream results. Return dict: ensembl_base -> list of stage records."""
    print(f"Loading fibrosis stage trajectory: {FIBROSIS_STAGE_CSV}")
    df = pd.read_csv(
        FIBROSIS_STAGE_CSV,
        usecols=["gene", "stage_label", "logFC", "padj"],
    )
    df["ensembl_base"] = df["gene"].apply(strip_version)
    print(f"  {len(df)} rows, {df['stage_label'].nunique()} stages")

    result = {}
    for ens, grp in df.groupby("ensembl_base"):
        records = []
        for _, r in grp.iterrows():
            lfc = safe_round(r["logFC"], 3)
            padj = safe_round(r["padj"], 4)
            if lfc is not None:
                records.append({
                    "stage": str(r["stage_label"]),
                    "logfc": lfc,
                    "padj": padj,
                })
        if records:
            # Sort by stage label for consistent order
            records.sort(key=lambda x: x["stage"])
            result[ens] = records
    return result


def load_stage_vs_control():
    """Load fibrosis stage vs F0 control trajectory.

    Prefer aggregated per-cohort meta (fibrosis_stage_trajectory_vs_F0.csv).
    Return dict: ensembl_base -> list of {stage, logfc, padj, n_cohorts}.
    Returns {} if pipeline output not yet present.
    """
    if not FIBROSIS_VS_F0_CSV.exists():
        print(f"[stage_vs_control] file not present yet ({FIBROSIS_VS_F0_CSV.name}); skipping")
        return {}
    print(f"Loading fibrosis stage vs F0: {FIBROSIS_VS_F0_CSV}")
    df = pd.read_csv(FIBROSIS_VS_F0_CSV)
    df["ensembl_base"] = df["gene"].apply(strip_version)
    print(f"  {len(df)} rows, {df['stage'].nunique()} stages")
    result = {}
    for ens, grp in df.groupby("ensembl_base"):
        records = []
        for _, r in grp.iterrows():
            lfc = safe_round(r["logFC"], 3)
            padj = safe_round(r["padj"], 4)
            if lfc is not None:
                rec = {"stage": str(r["stage"]), "logfc": lfc, "padj": padj}
                nc = safe(r.get("n_cohorts"))
                if nc is not None:
                    rec["n_cohorts"] = nc
                records.append(rec)
        if records:
            records.sort(key=lambda x: x["stage"])
            result[ens] = records
    return result


def load_stage_vs_healthy():
    """Load fibrosis stage vs Healthy control trajectory (F_i vs Healthy).

    Prefers the shim (fibrosis_stage_trajectory_vs_healthy.csv), falls back to
    the raw dream CSV (fibrosis_vs_healthy_dream.csv). Returns
    dict: ensembl_base -> list of {stage, logfc, padj, [n_cohorts]}.
    Returns {} if neither file is present (05d job still running).
    """
    if FIBROSIS_VS_HEALTHY_CSV.exists():
        print(f"[fib_vs_healthy] loading shim: {FIBROSIS_VS_HEALTHY_CSV}")
        df = pd.read_csv(FIBROSIS_VS_HEALTHY_CSV)
        col = "stage" if "stage" in df.columns else "stage_label"
    elif FIBROSIS_VS_HEALTHY_DREAM_CSV.exists():
        print(f"[fib_vs_healthy] loading dream: {FIBROSIS_VS_HEALTHY_DREAM_CSV}")
        df = pd.read_csv(
            FIBROSIS_VS_HEALTHY_DREAM_CSV,
            usecols=["gene", "stage_label", "logFC", "padj"],
        )
        df = df.rename(columns={"stage_label": "stage"})
        col = "stage"
    else:
        print("[fib_vs_healthy] no source present yet — skipping")
        return {}
    df["ensembl_base"] = df["gene"].apply(strip_version)
    print(f"  {len(df)} rows, {df[col].nunique()} stages")
    result = {}
    for ens, grp in df.groupby("ensembl_base"):
        records = []
        for _, r in grp.iterrows():
            lfc = safe_round(r["logFC"], 3)
            padj = safe_round(r["padj"], 4)
            if lfc is not None:
                rec = {"stage": str(r[col]), "logfc": lfc, "padj": padj}
                nc = safe(r.get("n_cohorts")) if "n_cohorts" in df.columns else None
                if nc is not None:
                    rec["n_cohorts"] = nc
                records.append(rec)
        if records:
            records.sort(key=lambda x: x["stage"])
            result[ens] = records
    return result


def load_nas_stage():
    """Load NAS trajectory. Prefer nas_trajectory.csv, else one_vs_rest_nas_dream.csv.

    Return dict: ensembl_base -> list of {nas, logfc, padj}.
    """
    if NAS_TRAJECTORY_CSV.exists():
        print(f"Loading NAS trajectory: {NAS_TRAJECTORY_CSV}")
        df = pd.read_csv(NAS_TRAJECTORY_CSV)
        col = "nas_score" if "nas_score" in df.columns else "stage"
    elif NAS_ONE_VS_REST_CSV.exists():
        print(f"Loading NAS fallback: {NAS_ONE_VS_REST_CSV}")
        df = pd.read_csv(NAS_ONE_VS_REST_CSV, usecols=["gene", "stage_label", "logFC", "padj"])
        df = df.rename(columns={"stage_label": "nas_score"})
        col = "nas_score"
    else:
        print("[nas] no NAS file present; skipping")
        return {}
    df["ensembl_base"] = df["gene"].apply(strip_version)
    print(f"  {len(df)} rows, {df[col].nunique()} NAS levels")
    result = {}
    for ens, grp in df.groupby("ensembl_base"):
        records = []
        for _, r in grp.iterrows():
            lfc = safe_round(r["logFC"], 3)
            padj = safe_round(r["padj"], 4)
            if lfc is not None:
                records.append({
                    "nas": str(r[col]),
                    "logfc": lfc,
                    "padj": padj,
                })
        if records:
            records.sort(key=lambda x: x["nas"])
            result[ens] = records
    return result


def _load_nas_baseline_csv(shim_path: Path, dream_path: Path, tag: str):
    """Generic loader for NAS-vs-fixed-baseline trajectory files.

    Prefers the shim (nas_trajectory_vs_*.csv), falls back to the raw dream
    CSV (nas_vs_*_dream.csv). Returns dict: ensembl_base -> list of
    {nas, logfc, padj}. Returns {} if neither file is present (job still
    running).
    """
    if shim_path.exists():
        print(f"[{tag}] loading shim: {shim_path}")
        df = pd.read_csv(shim_path)
        col = "nas_score" if "nas_score" in df.columns else "stage_label"
    elif dream_path.exists():
        print(f"[{tag}] loading dream: {dream_path}")
        df = pd.read_csv(dream_path, usecols=["gene", "stage_label", "logFC", "padj"])
        df = df.rename(columns={"stage_label": "nas_score"})
        col = "nas_score"
    else:
        print(f"[{tag}] no source present yet; skipping")
        return {}
    df["ensembl_base"] = df["gene"].apply(strip_version)
    print(f"  {len(df)} rows, {df[col].nunique()} NAS levels")
    result = {}
    for ens, grp in df.groupby("ensembl_base"):
        records = []
        for _, r in grp.iterrows():
            lfc = safe_round(r["logFC"], 3)
            padj = safe_round(r["padj"], 4)
            if lfc is not None:
                records.append({
                    "nas": str(r[col]),
                    "logfc": lfc,
                    "padj": padj,
                })
        if records:
            records.sort(key=lambda x: x["nas"])
            result[ens] = records
    return result


def load_nas_vs_nas0():
    """NAS_i vs NAS_0 baseline trajectory (05c_nas_vs_baseline_dream.R)."""
    return _load_nas_baseline_csv(
        NAS_VS_NAS0_CSV, NAS_VS_NAS0_DREAM_CSV, "nas_vs_nas0"
    )


def load_nas_vs_healthy():
    """NAS_i vs Healthy (diagnosis_harmonized == Control) trajectory."""
    return _load_nas_baseline_csv(
        NAS_VS_HEALTHY_CSV, NAS_VS_HEALTHY_DREAM_CSV, "nas_vs_healthy"
    )


def load_per_cohort_contrasts():
    """Load per-cohort contrast DE files.

    Returns dict: contrast_name -> {ensembl_base -> list of {dataset, logfc, pval, padj}}.
    Only includes cohorts with a present *_de.csv under per_study_contrasts/{contrast}/.
    """
    result = {c: {} for c in PER_COHORT_CONTRAST_NAMES}
    if not PER_STUDY_CONTRASTS_DIR.is_dir():
        print(f"[per_cohort_contrasts] dir not present: {PER_STUDY_CONTRASTS_DIR}")
        return result
    for contrast in PER_COHORT_CONTRAST_NAMES:
        cdir = PER_STUDY_CONTRASTS_DIR / contrast
        if not cdir.is_dir():
            continue
        frames = []
        for fp in sorted(cdir.glob("*_de.csv")):
            cohort = fp.name.replace("_de.csv", "")
            if cohort in EXCLUDED_FROM_PER_COHORT:
                continue
            try:
                d = pd.read_csv(fp, usecols=["gene", "logFC", "P.Value", "padj"])
                d["dataset"] = cohort
                frames.append(d)
            except Exception as e:
                print(f"  WARNING: skip {fp}: {e}")
        if not frames:
            continue
        df = pd.concat(frames, ignore_index=True)
        df["ensembl_base"] = df["gene"].apply(strip_version)
        by_gene = {}
        for ens, grp in df.groupby("ensembl_base"):
            records = []
            for _, r in grp.iterrows():
                lfc = safe_round(r["logFC"], 3)
                pv = safe_round(r["P.Value"], 4)
                padj = safe_round(r["padj"], 4)
                if lfc is not None and pv is not None:
                    records.append({
                        "dataset": str(r["dataset"]),
                        "logfc": lfc,
                        "pval": pv,
                        "padj": padj,
                    })
            if records:
                by_gene[ens] = records
        result[contrast] = by_gene
        print(f"[per_cohort_contrasts] {contrast}: {len(by_gene):,} genes across "
              f"{df['dataset'].nunique()} cohorts")
    return result


def load_coloc_gene_level():
    """Load gene-level COLOC summary. Return dict: symbol -> record."""
    print(f"Loading COLOC gene-level: {COLOC_GENE_CSV}")
    df = pd.read_csv(COLOC_GENE_CSV)
    print(f"  {len(df)} genes")

    result = {}
    for _, r in df.iterrows():
        sym = str(r.get("gene", "")).strip()
        if not sym or sym == "nan" or sym == "":
            continue
        result[sym] = {
            "coloc_pp4_max": safe_round(r.get("coloc_best_pp4"), 4),
            "coloc_best_gwas": safe(r.get("coloc_best_gwas")),
            "coloc_n_gwas_05": safe(r.get("coloc_n_gwas_h4_05")),
            "coloc_n_gwas_08": safe(r.get("coloc_n_gwas_h4_08")),
        }
    return result


def load_coloc_per_gwas():
    """Load per-GWAS PP.H4 from all_gwas file. Return dict: symbol -> {gwas: pp4}."""
    print(f"Loading per-GWAS COLOC: {COLOC_ALL_GWAS_CSV}")
    df = pd.read_csv(
        COLOC_ALL_GWAS_CSV,
        usecols=["gene", "gwas_name", "PP.H4.abf"],
    )
    print(f"  {len(df)} gene-GWAS entries, {df['gwas_name'].nunique()} GWAS")

    # Only keep entries with PP.H4 > 0.01 to reduce noise
    df = df[df["PP.H4.abf"] > 0.01]
    print(f"  {len(df)} entries with PP.H4 > 0.01")

    result = {}
    for sym, grp in df.groupby("gene"):
        sym = str(sym).strip()
        if not sym or sym == "nan":
            continue
        by_gwas = {}
        for _, r in grp.iterrows():
            gwas = str(r["gwas_name"])
            pp4 = safe_round(r["PP.H4.abf"], 4)
            if pp4 is not None and pp4 > 0.01:
                by_gwas[gwas] = pp4
        if by_gwas:
            result[sym] = by_gwas
    return result


def load_deconv_attribution():
    """Load deconvolution attribution. Return dict: ensembl_base -> category."""
    print(f"Loading deconvolution attribution: {DECONV_CSV}")
    df = pd.read_csv(DECONV_CSV, usecols=["gene", "category"])
    df["ensembl_base"] = df["gene"].apply(strip_version)
    print(f"  {len(df)} genes")

    result = {}
    for _, r in df.iterrows():
        cat = str(r["category"]).strip()
        if cat and cat != "nan" and cat != "Not_significant":
            result[r["ensembl_base"]] = cat
    return result


def load_sex_causal():
    """Load sex causal scores. Return dict: symbol -> record."""
    print(f"Loading sex causal scores: {SEX_CAUSAL_CSV}")
    df = pd.read_csv(SEX_CAUSAL_CSV)
    print(f"  {len(df)} genes")

    result = {}
    for _, r in df.iterrows():
        sym = str(r.get("human_symbol", "")).strip()
        if not sym or sym == "nan":
            continue
        sex_class = safe(r.get("sex_class"))
        if sex_class and str(sex_class).strip().lower() in ("", "nan", "not_significant"):
            sex_class = None
        result[sym] = {
            "sex_class": sex_class,
            "logfc_female": safe_round(r.get("logFC_F"), 3),
            "logfc_male": safe_round(r.get("logFC_M"), 3),
            "padj_female": safe_round(r.get("padj_F"), 4),
            "padj_male": safe_round(r.get("padj_M"), 4),
        }
    return result


def load_drug_progression():
    """Load drug target progression classification. Return dict: symbol -> record."""
    print(f"Loading drug target progression: {DRUG_PROGRESSION_CSV}")
    df = pd.read_csv(DRUG_PROGRESSION_CSV)
    print(f"  {len(df)} drug targets")

    result = {}
    for _, r in df.iterrows():
        sym = str(r.get("target_gene", "")).strip()
        if not sym or sym == "nan":
            continue
        result[sym] = {
            "drug": safe(r.get("drug")),
            "stage": safe(r.get("stage")),
            "moa": safe(r.get("moa")),
            "target_class": safe(r.get("target_class")),
        }
    return result


def load_clinical_drugs():
    """Load clinical drug validation. Return dict: target_gene -> list of drug records."""
    print(f"Loading clinical drugs: {CLINICAL_DRUGS_CSV}")
    df = pd.read_csv(CLINICAL_DRUGS_CSV, low_memory=False)
    print(f"  {len(df)} drug entries")

    result = {}
    for _, r in df.iterrows():
        sym = str(r.get("target_gene", "")).strip()
        if not sym or sym == "nan":
            continue
        entry = {
            "drug": safe(r.get("drug")),
            "stage": safe(r.get("stage")),
            "moa": safe(r.get("moa")),
            "support": safe(r.get("atlas_support")),
        }
        result.setdefault(sym, []).append(entry)
    return result


def load_lincs_ranked():
    """Load LINCS ranked compounds. Return dict: target_gene -> list of compound records."""
    print(f"Loading LINCS ranked: {LINCS_RANKED_CSV}")
    df = pd.read_csv(LINCS_RANKED_CSV, low_memory=False)
    print(f"  {len(df)} compounds")

    result = {}
    for _, r in df.iterrows():
        # Targets can come from target.x or dgidb_targets columns
        targets = set()
        for col in ("target.x", "dgidb_targets"):
            val = r.get(col)
            if pd.notna(val):
                val_str = str(val).strip()
                if val_str and val_str != "nan":
                    # Targets may be semicolon or comma separated
                    for t in val_str.replace(";", ",").split(","):
                        t = t.strip()
                        if t and t != "nan":
                            targets.add(t)

        name = safe(r.get("display_name")) or safe(r.get("cmap_name")) or safe(r.get("pert_iname"))
        score = safe_round(r.get("composite_score"), 3)
        moa = safe(r.get("moa.x")) or safe(r.get("MOAss"))

        entry = {"name": name, "score": score}
        if moa and str(moa).strip().lower() not in ("", "nan", '""'):
            entry["moa"] = str(moa).strip().strip('"')

        for tg in targets:
            result.setdefault(tg, []).append(entry)

    return result


def load_subtype_markers():
    """Load NMF subtype markers. Return dict: ensembl_base -> {subtype, direction, logfc}."""
    print(f"Loading subtype markers: {SUBTYPE_MARKERS_CSV}")
    df = pd.read_csv(SUBTYPE_MARKERS_CSV)
    df["ensembl_base"] = df["gene"].apply(strip_version)
    print(f"  {len(df)} markers")

    result = {}
    for _, r in df.iterrows():
        ens = r["ensembl_base"]
        result[ens] = {
            "subtype": safe(r.get("subtype")),
            "direction": safe(r.get("direction")),
            "logfc": safe_round(r.get("logFC"), 3),
        }
    return result


def load_pseudobulk_de():
    """Load pseudobulk DE results.

    Returns two dicts:
      - ensembl_map: ensembl_base -> list of cell-type records (for ENSG* gene IDs)
      - symbol_map: gene_symbol -> list of cell-type records (for symbol gene IDs)

    The gene column in pseudobulk DE files is mixed: some rows have Ensembl IDs,
    some have gene symbols. We split into two lookup dicts accordingly.
    """
    print(f"Loading pseudobulk DE from {PSEUDOBULK_DIR}")
    frames = []
    for fp in sorted(PSEUDOBULK_DIR.glob("*_de.csv")):
        if fp.name.startswith("archive"):
            continue
        try:
            chunk = pd.read_csv(
                fp, usecols=["gene", "cell_type", "logFC", "padj"]
            )
            frames.append(chunk)
        except Exception as e:
            print(f"  WARNING: skipping {fp.name}: {e}")

    if not frames:
        print("  WARNING: no pseudobulk DE files found")
        return {}, {}

    df = pd.concat(frames, ignore_index=True)
    # Filter to significant or notable entries to reduce noise
    df = df[df["padj"].notna() & (df["padj"] < 0.1)]
    print(f"  {len(df)} significant (padj<0.1) pseudobulk DE records")

    # Split by ID type: ENSG* = ensembl, everything else = symbol
    is_ensembl = df["gene"].str.startswith("ENSG", na=False)
    df_ens = df[is_ensembl].copy()
    df_sym = df[~is_ensembl].copy()
    df_ens["key"] = df_ens["gene"].apply(strip_version)
    df_sym["key"] = df_sym["gene"].str.strip()

    print(f"    {len(df_ens)} with Ensembl IDs, {len(df_sym)} with gene symbols")

    def _group_to_map(sub_df):
        result = {}
        for key, grp in sub_df.groupby("key"):
            records = []
            for _, r in grp.iterrows():
                records.append({
                    "cell_type": str(r["cell_type"]),
                    "logfc": safe_round(r["logFC"], 3),
                    "padj": safe_round(r["padj"], 4),
                })
            if records:
                records.sort(key=lambda x: x["cell_type"])
                result[key] = records
        return result

    ensembl_map = _group_to_map(df_ens)
    symbol_map = _group_to_map(df_sym)
    return ensembl_map, symbol_map


def load_proteomics_de():
    """Load proteomics DE results. Returns dict: gene_symbol -> best record.
    Header: logFC, AveExpr, t, pvalue, padj, B, gene, dataset, ID
    For genes with entries in multiple datasets, pick the one with lowest padj.
    """
    if not PROTEOMICS_DE_CSV.exists():
        print(f"  WARNING: proteomics DE file not found: {PROTEOMICS_DE_CSV}")
        return {}
    print(f"Loading proteomics DE: {PROTEOMICS_DE_CSV}")
    df = pd.read_csv(PROTEOMICS_DE_CSV)
    df = df[df["gene"].notna() & df["padj"].notna()].copy()
    # For each gene, keep row with smallest padj
    df = df.sort_values("padj").drop_duplicates("gene", keep="first")
    print(f"  {len(df)} unique proteomics genes")
    result = {}
    for _, r in df.iterrows():
        sym = str(r["gene"]).strip()
        if not sym or sym.lower() == "nan":
            continue
        rec = {}
        lfc = safe_round(r.get("logFC"), 3)
        padj = safe_round(r.get("padj"), 6)
        tstat = safe_round(r.get("t"), 3)
        if lfc is not None:
            rec["logfc"] = lfc
        if padj is not None:
            rec["padj"] = padj
        if tstat is not None:
            rec["tstat"] = tstat
        ds = safe(r.get("dataset"))
        if ds:
            rec["dataset"] = ds
        if rec:
            result[sym] = rec
    return result


# ---------------------------------------------------------------------------
# Profile builder
# ---------------------------------------------------------------------------


def build_profile(
    row,
    ensembl_base,
    per_study_map,
    fibrosis_map,
    coloc_gene_map,
    coloc_gwas_map,
    deconv_map,
    sex_map,
    drug_prog_map,
    clinical_drug_map,
    lincs_map,
    subtype_map,
    pseudobulk_ens_map,
    pseudobulk_sym_map,
    stage_vs_control_map=None,
    stage_vs_healthy_map=None,
    nas_map=None,
    nas_vs_nas0_map=None,
    nas_vs_healthy_map=None,
    per_cohort_contrasts_map=None,
    proteomics_map=None,
):
    """Build a single gene profile dict. Omits empty sections."""
    symbol = str(row["human_symbol"]).strip()
    ensembl_id = safe(row.get("ensembl_id"))
    biotype = safe(row.get("gene_biotype"))

    profile = {
        "symbol": symbol,
        "ensembl_id": ensembl_id,
    }
    if biotype:
        profile["biotype"] = biotype

    # ------- EXPRESSION -------
    expression = {}
    bulk_lfc = safe_round(row.get("bulk_logFC"), 3)
    bulk_padj = safe_round(row.get("bulk_padj"), 4)
    bulk_tstat = safe_round(row.get("bulk_tstat"), 3)
    if bulk_lfc is not None:
        expression["bulk_logfc"] = bulk_lfc
    if bulk_padj is not None:
        expression["bulk_padj"] = bulk_padj
    if bulk_tstat is not None:
        expression["bulk_tstat"] = bulk_tstat

    # Legacy: single dream contrast per cohort
    per_cohort = per_study_map.get(ensembl_base)
    if per_cohort:
        expression["per_cohort_dream"] = per_cohort

    # New: per-cohort per-contrast (disease/MASH/MASL/MASH_vs_MASL)
    if per_cohort_contrasts_map:
        pc_out = {}
        for cname, gene_map in per_cohort_contrasts_map.items():
            rec = gene_map.get(ensembl_base)
            if rec:
                pc_out[cname] = rec
        if pc_out:
            expression["per_cohort"] = pc_out

    # Fibrosis stage trajectory — primary = vs F0 control, secondary = vs rest
    if stage_vs_control_map:
        vs_ctrl = stage_vs_control_map.get(ensembl_base)
        if vs_ctrl:
            expression["stage_trajectory_vs_control"] = vs_ctrl

    stage_rest = fibrosis_map.get(ensembl_base)
    if stage_rest:
        expression["stage_trajectory_vs_rest"] = stage_rest

    # Fibrosis stage vs Healthy controls (05d_fibrosis_vs_healthy_dream.R).
    if stage_vs_healthy_map:
        vs_healthy = stage_vs_healthy_map.get(ensembl_base)
        if vs_healthy:
            expression["stage_trajectory_vs_healthy"] = vs_healthy

    # NAS trajectory (one-vs-rest; primary)
    if nas_map:
        nas_rec = nas_map.get(ensembl_base)
        if nas_rec:
            expression["nas_trajectory"] = nas_rec
    # NAS trajectory vs NAS_0 baseline
    if nas_vs_nas0_map:
        rec = nas_vs_nas0_map.get(ensembl_base)
        if rec:
            expression["nas_trajectory_vs_nas0"] = rec
    # NAS trajectory vs Healthy baseline
    if nas_vs_healthy_map:
        rec = nas_vs_healthy_map.get(ensembl_base)
        if rec:
            expression["nas_trajectory_vs_healthy"] = rec

    if expression:
        profile["expression"] = expression

    # ------- CAUSAL -------
    causal = {}
    coloc_gene = coloc_gene_map.get(symbol, {})
    if coloc_gene:
        for k, v in coloc_gene.items():
            if v is not None:
                causal[k] = v

    coloc_by_gwas = coloc_gwas_map.get(symbol)
    if coloc_by_gwas:
        causal["coloc_by_gwas"] = coloc_by_gwas

    twas_z = safe_round(row.get("twas_z"), 3)
    twas_pval = safe_round(row.get("twas_pval"), 4)
    if twas_z is not None:
        causal["twas_z"] = twas_z
    if twas_pval is not None:
        causal["twas_pval"] = twas_pval

    if causal:
        profile["causal"] = causal

    # ------- CELL TYPE -------
    celltype = {}
    attr_class = deconv_map.get(ensembl_base)
    if attr_class:
        celltype["attribution_class"] = attr_class

    # Pseudobulk DE: check both ensembl-keyed and symbol-keyed maps
    pseudobulk = None
    if ensembl_base:
        pseudobulk = pseudobulk_ens_map.get(ensembl_base)
    if not pseudobulk:
        pseudobulk = pseudobulk_sym_map.get(symbol)
    if pseudobulk:
        celltype["pseudobulk_de"] = pseudobulk

    if celltype:
        profile["celltype"] = celltype

    # ------- SEX & SUBTYPE -------
    sex_subtype = {}
    sex_data = sex_map.get(symbol)
    if sex_data:
        for k, v in sex_data.items():
            if v is not None:
                sex_subtype[k] = v

    # Also include atlas sex columns if sex_map entry is missing
    if not sex_data:
        atlas_sex_class = safe(row.get("sex_class"))
        if atlas_sex_class and str(atlas_sex_class).lower() not in ("", "nan", "not_significant"):
            sex_subtype["sex_class"] = atlas_sex_class
        atlas_sex_padj = safe_round(row.get("sex_interaction_padj"), 4)
        if atlas_sex_padj is not None:
            sex_subtype["sex_interaction_padj"] = atlas_sex_padj
        atlas_lfc_f = safe_round(row.get("bulk_logFC_F"), 3)
        atlas_lfc_m = safe_round(row.get("bulk_logFC_M"), 3)
        if atlas_lfc_f is not None:
            sex_subtype["logfc_female"] = atlas_lfc_f
        if atlas_lfc_m is not None:
            sex_subtype["logfc_male"] = atlas_lfc_m

    subtype = subtype_map.get(ensembl_base)
    if subtype:
        for k, v in subtype.items():
            if v is not None:
                sex_subtype[f"nmf_{k}"] = v

    if sex_subtype:
        profile["sex_subtype"] = sex_subtype

    # ------- THERAPEUTIC -------
    therapeutic = {}
    dgidb = coerce_r_bool(row.get("dgidb_druggable"))
    ot = coerce_r_bool(row.get("opentargets_drug"))
    if dgidb:
        therapeutic["dgidb_druggable"] = True
    if ot:
        therapeutic["opentargets_drug"] = True

    clinical_drugs = clinical_drug_map.get(symbol)
    if clinical_drugs:
        therapeutic["drugs"] = clinical_drugs

    lincs_compounds = lincs_map.get(symbol)
    if lincs_compounds:
        therapeutic["lincs_compounds"] = lincs_compounds

    drug_prog = drug_prog_map.get(symbol)
    if drug_prog:
        prog_class = drug_prog.get("target_class")
        if prog_class:
            therapeutic["progression_class"] = prog_class
        prog_drug = drug_prog.get("drug")
        if prog_drug:
            therapeutic["progression_drug"] = prog_drug

    if therapeutic:
        profile["therapeutic"] = therapeutic

    # ------- PROTEOMICS -------
    if proteomics_map:
        prot_rec = proteomics_map.get(symbol)
        if prot_rec:
            profile["proteomics"] = prot_rec

    # ------- EXTERNAL LINKS -------
    external = {
        "genecards": f"https://www.genecards.org/cgi-bin/carddisp.pl?gene={symbol}",
        "pubmed": f"https://pubmed.ncbi.nlm.nih.gov/?term={symbol}+MASLD",
    }
    if ensembl_id:
        external["opentargets"] = f"https://platform.opentargets.org/target/{ensembl_id}"
    external["gtex"] = f"https://gtexportal.org/home/gene/{symbol}"
    profile["external"] = external

    return profile


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Generate per-gene JSON profiles for MASLD Atlas v2."
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=os.path.join(os.path.dirname(__file__), "..", "public", "data"),
        help="Output directory (genes/ subdir will be created inside)",
    )
    parser.add_argument(
        "--symbols",
        type=str,
        default=None,
        help="Comma-separated list of gene symbols to regenerate only. "
             "If omitted, all atlas genes are processed.",
    )
    args = parser.parse_args()

    output_dir = Path(args.output_dir).resolve()
    genes_dir = output_dir / "genes"
    genes_dir.mkdir(parents=True, exist_ok=True)
    print(f"Output directory: {genes_dir}")

    t0 = time.time()

    # ------ Load all data ------
    atlas = load_atlas()
    per_study_map = load_per_study_de()
    fibrosis_map = load_fibrosis_stage()
    stage_vs_control_map = load_stage_vs_control()
    stage_vs_healthy_map = load_stage_vs_healthy()
    nas_map = load_nas_stage()
    nas_vs_nas0_map = load_nas_vs_nas0()
    nas_vs_healthy_map = load_nas_vs_healthy()
    per_cohort_contrasts_map = load_per_cohort_contrasts()
    coloc_gene_map = load_coloc_gene_level()
    coloc_gwas_map = load_coloc_per_gwas()
    deconv_map = load_deconv_attribution()
    sex_map = load_sex_causal()
    drug_prog_map = load_drug_progression()
    clinical_drug_map = load_clinical_drugs()
    lincs_map = load_lincs_ranked()
    subtype_map = load_subtype_markers()
    pseudobulk_ens_map, pseudobulk_sym_map = load_pseudobulk_de()
    proteomics_map = load_proteomics_de()

    t_load = time.time()
    print(f"\nData loading: {t_load - t0:.1f}s")

    # ------ Generate profiles ------
    n_generated = 0
    n_skipped = 0
    total_bytes = 0

    symbols_filter = None
    if args.symbols:
        symbols_filter = {s.strip() for s in args.symbols.split(",") if s.strip()}
        print(f"Restricting to {len(symbols_filter)} requested symbols")

    for idx, row in atlas.iterrows():
        symbol = row.get("human_symbol")
        if pd.isna(symbol) or str(symbol).strip() == "":
            n_skipped += 1
            continue
        if symbols_filter is not None and str(symbol).strip() not in symbols_filter:
            continue

        symbol = str(symbol).strip()
        ensembl_id = str(row.get("ensembl_id", "")).strip()
        ensembl_base = ensembl_id if ensembl_id and ensembl_id != "nan" else None

        profile = build_profile(
            row=row,
            ensembl_base=ensembl_base,
            per_study_map=per_study_map,
            fibrosis_map=fibrosis_map,
            coloc_gene_map=coloc_gene_map,
            coloc_gwas_map=coloc_gwas_map,
            deconv_map=deconv_map,
            sex_map=sex_map,
            drug_prog_map=drug_prog_map,
            clinical_drug_map=clinical_drug_map,
            lincs_map=lincs_map,
            subtype_map=subtype_map,
            pseudobulk_ens_map=pseudobulk_ens_map,
            pseudobulk_sym_map=pseudobulk_sym_map,
            stage_vs_control_map=stage_vs_control_map,
            stage_vs_healthy_map=stage_vs_healthy_map,
            nas_map=nas_map,
            nas_vs_nas0_map=nas_vs_nas0_map,
            nas_vs_healthy_map=nas_vs_healthy_map,
            per_cohort_contrasts_map=per_cohort_contrasts_map,
            proteomics_map=proteomics_map,
        )

        # Sanitize filename: replace slashes and other unsafe chars
        safe_symbol = symbol.replace("/", "_").replace("\\", "_").replace(" ", "_")
        out_path = genes_dir / f"{safe_symbol}.json"

        with open(out_path, "w") as f:
            json.dump(profile, f, separators=(",", ":"), default=safe)

        total_bytes += out_path.stat().st_size
        n_generated += 1

        if n_generated % 5000 == 0:
            elapsed = time.time() - t_load
            print(
                f"  [{n_generated:,}/{atlas.shape[0]:,}] "
                f"generated ({elapsed:.1f}s, {total_bytes / (1024*1024):.1f} MB)"
            )

    t_end = time.time()
    total_mb = total_bytes / (1024 * 1024)

    print(f"\n{'=' * 60}")
    print(f"Done in {t_end - t0:.1f}s total ({t_end - t_load:.1f}s generation)")
    print(f"  Generated: {n_generated:,} gene profiles")
    print(f"  Skipped:   {n_skipped:,} (missing symbol)")
    print(f"  Total size: {total_mb:.1f} MB")
    print(f"  Avg size:   {total_bytes / max(n_generated, 1):.0f} bytes")
    print(f"  Output dir: {genes_dir}")


if __name__ == "__main__":
    main()
