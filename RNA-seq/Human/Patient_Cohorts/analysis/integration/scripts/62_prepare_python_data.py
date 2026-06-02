#!/usr/bin/env python3
"""
62_prepare_python_data.py
Convert R-generated feature matrices and metadata to HDF5 for Python modeling.

Outputs:
  - prepared_data.h5: HDF5 file containing:
    - rank_expression: (genes x samples) rank-transformed matrix
    - zscore_expression: (genes x samples) z-score-transformed matrix
    - raw_logcpm: (genes x samples) raw logCPM matrix
    - gene_names: gene identifiers
    - sample_ids: sample identifiers
    - metadata: per-sample metadata (dataset, NAS, fibrosis, sex, diagnosis)
    - loco_folds_fibrosis: per-sample fold assignment for fibrosis LOCO
    - loco_folds_nas: per-sample fold assignment for NAS LOCO
    - feature_candidates: top 3000 gene names and stats

Usage: python 62_prepare_python_data.py
SLURM: cpu, 4 CPUs, 64GB RAM, 48h
"""

import os
import sys
import numpy as np
import pandas as pd
import h5py
from pathlib import Path

# --- Paths ---
BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR = os.path.join(INT, "results")
OUTDIR = os.path.join(RDIR, "staging_classifier")
os.makedirs(OUTDIR, exist_ok=True)

print("=== 62: Prepare Python Data ===")
print(f"Output dir: {OUTDIR}")

# --- Load R-saved RDS files via CSV export ---
# Strategy: Script 61 already saved the key matrices as RDS.
# We use a small R subprocess to convert RDS -> CSV for Python consumption.

import subprocess
import tempfile

RSCRIPT = "/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript"

def rds_to_numpy(rds_path, tmpdir):
    """Convert an RDS matrix to numpy array + row/col names via R subprocess."""
    csv_path = os.path.join(tmpdir, "matrix.csv")
    r_code = f"""
    m <- readRDS("{rds_path}")
    write.csv(m, "{csv_path}", quote=FALSE)
    """
    result = subprocess.run(
        [RSCRIPT, "-e", r_code],
        capture_output=True, text=True, timeout=600
    )
    if result.returncode != 0:
        print(f"R error: {result.stderr}")
        raise RuntimeError(f"Failed to convert {rds_path}")

    df = pd.read_csv(csv_path, index_col=0)
    return df

# --- Load matrices ---
print("\nLoading expression matrices from R...")

with tempfile.TemporaryDirectory() as tmpdir:
    # Rank expression (top 3000 genes)
    rank_path = os.path.join(OUTDIR, "rank_expression_matrix.rds")
    if os.path.exists(rank_path):
        print(f"  Loading rank matrix: {rank_path}")
        rank_df = rds_to_numpy(rank_path, tmpdir)
        print(f"  Rank matrix: {rank_df.shape}")
    else:
        sys.exit(f"ERROR: {rank_path} not found. Run Script 61 first.")

    # Z-score expression
    zscore_path = os.path.join(OUTDIR, "zscore_expression_matrix.rds")
    if os.path.exists(zscore_path):
        print(f"  Loading z-score matrix: {zscore_path}")
        zscore_df = rds_to_numpy(zscore_path, tmpdir)
        print(f"  Z-score matrix: {zscore_df.shape}")
    else:
        print("  Z-score matrix not found, will skip")
        zscore_df = None

    # Raw logCPM
    raw_path = os.path.join(OUTDIR, "raw_logcpm_matrix.rds")
    if os.path.exists(raw_path):
        print(f"  Loading raw logCPM: {raw_path}")
        raw_df = rds_to_numpy(raw_path, tmpdir)
        print(f"  Raw logCPM: {raw_df.shape}")
    else:
        print("  Raw logCPM not found, will skip")
        raw_df = None

# --- Load metadata ---
print("\nLoading metadata...")
meta = pd.read_csv(os.path.join(INT, "metadata/unified_metadata.csv"))
print(f"  Unified metadata: {meta.shape}")

# Load QC report for pass_technical
qc = pd.read_csv(os.path.join(INT, "qc/sample_qc_report.csv"))
pass_samples = set(qc[qc["pass_technical"] == True]["sample_id"])
print(f"  QC-passing: {len(pass_samples)}")

# Merge metadata for samples in the expression matrix
sample_ids = rank_df.columns.tolist()
meta_sub = meta[meta["sample_id"].isin(sample_ids)].copy()
meta_sub = meta_sub.set_index("sample_id").loc[sample_ids].reset_index()
print(f"  Samples with metadata: {len(meta_sub)} of {len(sample_ids)}")

# --- Define LOCO folds ---
print("\nDefining LOCO folds...")

# Fibrosis LOCO: 6 cohorts with individual fibrosis staging
FIB_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478",
                 "GSE193066", "GSE240729"]

# NAS LOCO: 5 cohorts with per-sample NAS scores
NAS_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066"]

# Assign fold = dataset name for LOCO (fold = held-out cohort)
meta_sub["loco_fold_fibrosis"] = meta_sub.apply(
    lambda r: r["dataset"] if (r["dataset"] in FIB_DATASETS and
                                pd.notna(r.get("fibrosis_stage"))) else "excluded",
    axis=1
)
meta_sub["loco_fold_nas"] = meta_sub.apply(
    lambda r: r["dataset"] if (r["dataset"] in NAS_DATASETS and
                                pd.notna(r.get("nas_score"))) else "excluded",
    axis=1
)

fib_fold_counts = meta_sub[meta_sub["loco_fold_fibrosis"] != "excluded"]["loco_fold_fibrosis"].value_counts()
nas_fold_counts = meta_sub[meta_sub["loco_fold_nas"] != "excluded"]["loco_fold_nas"].value_counts()
print(f"\nFibrosis LOCO folds ({len(fib_fold_counts)} cohorts):\n{fib_fold_counts.to_string()}")
print(f"\nNAS LOCO folds ({len(nas_fold_counts)} cohorts):\n{nas_fold_counts.to_string()}")

# CRITICAL ASSERTIONS: Ensure sufficient cohort coverage for LOCO
n_fib_folds = len(fib_fold_counts)
n_nas_folds = len(nas_fold_counts)
n_fib_samples = fib_fold_counts.sum()
n_nas_samples = nas_fold_counts.sum()
print(f"\nFibrosis: {n_fib_folds} folds, {n_fib_samples} samples")
print(f"NAS: {n_nas_folds} folds, {n_nas_samples} samples")

if n_fib_folds < 3:
    print(f"CRITICAL WARNING: Only {n_fib_folds} fibrosis LOCO folds (expected >=5).")
    print("Check that all fibrosis-annotated cohorts are present in the expression matrix.")
    print("Fibrosis datasets expected:", FIB_DATASETS)
    for ds in FIB_DATASETS:
        n_ds = (meta_sub["dataset"] == ds).sum()
        n_fib = ((meta_sub["dataset"] == ds) & (meta_sub["loco_fold_fibrosis"] != "excluded")).sum()
        print(f"  {ds}: {n_ds} total samples, {n_fib} with fibrosis annotation in LOCO")
if n_nas_folds < 3:
    print(f"CRITICAL WARNING: Only {n_nas_folds} NAS LOCO folds (expected >=4).")
    for ds in NAS_DATASETS:
        n_ds = (meta_sub["dataset"] == ds).sum()
        n_nas_ds = ((meta_sub["dataset"] == ds) & (meta_sub["loco_fold_nas"] != "excluded")).sum()
        print(f"  {ds}: {n_ds} total samples, {n_nas_ds} with NAS annotation in LOCO")

# --- Create NAS group and fibrosis stage columns for modeling ---
meta_sub["nas_group"] = meta_sub["nas_score"].apply(
    lambda x: int(min(x, 7)) if pd.notna(x) else -1
).astype(int)

# Collapse NAS into 4 groups for ordinal modeling
def nas_to_group4(nas):
    if pd.isna(nas) or nas < 0:
        return -1
    elif nas <= 2:
        return 0  # Low
    elif nas <= 4:
        return 1  # Moderate
    elif nas <= 6:
        return 2  # High
    else:
        return 3  # Very High

meta_sub["nas_group4"] = meta_sub["nas_score"].apply(nas_to_group4).astype(int)

meta_sub["fib_stage"] = meta_sub["fibrosis_stage"].apply(
    lambda x: int(x) if pd.notna(x) and x in [0,1,2,3,4] else -1
).astype(int)

# Binary clinical thresholds
meta_sub["fib_ge3"] = (meta_sub["fib_stage"] >= 3).astype(int)
meta_sub.loc[meta_sub["fib_stage"] < 0, "fib_ge3"] = -1

meta_sub["nas_ge5"] = -1
valid_nas = pd.notna(meta_sub["nas_score"]) & (meta_sub["nas_score"] >= 0)
meta_sub.loc[valid_nas, "nas_ge5"] = (meta_sub.loc[valid_nas, "nas_score"] >= 5).astype(int)

# Disease vs Control
meta_sub["is_disease"] = (meta_sub["group_binary"] == "Disease").astype(int)

# Severity 4-class (for Plan 2 Tier 2)
def assign_severity(row):
    if row.get("group_binary") == "Control":
        return 0  # Healthy
    diag = row.get("diagnosis_harmonized", "")
    fib = row.get("fib_stage", -1)
    if diag == "NAFL":
        return 1  # MASL
    elif diag in ["NASH", "Borderline"]:
        if fib >= 0 and fib <= 1:
            return 2  # MASH-noFibrosis
        elif fib >= 2:
            return 3  # MASH-Fibrosis
        else:
            return -1  # Unknown fibrosis
    return -1  # Unknown

meta_sub["severity4"] = meta_sub.apply(assign_severity, axis=1).astype(int)
sev_counts = meta_sub[meta_sub["severity4"] >= 0]["severity4"].value_counts().sort_index()
print(f"\nSeverity 4-class distribution:\n  0 (Healthy): {sev_counts.get(0, 0)}")
print(f"  1 (MASL): {sev_counts.get(1, 0)}")
print(f"  2 (MASH-noFib): {sev_counts.get(2, 0)}")
print(f"  3 (MASH-Fib): {sev_counts.get(3, 0)}")

# --- Load feature candidates ---
print("\nLoading feature candidates...")
fc_path = os.path.join(OUTDIR, "feature_candidates_3000.csv")
if os.path.exists(fc_path):
    feature_candidates = pd.read_csv(fc_path)
    print(f"  Feature candidates: {len(feature_candidates)} genes")
else:
    print("  Feature candidates not found, using all genes in rank matrix")
    feature_candidates = pd.DataFrame({"gene": rank_df.index.tolist()})

# --- Load multi-evidence atlas for Plan 3 fusion ---
atlas_path = os.path.join(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
atlas_df = None
if os.path.exists(atlas_path):
    atlas_df = pd.read_csv(atlas_path)
    print(f"  Multi-evidence atlas: {atlas_df.shape}")

# --- Save to HDF5 ---
print("\nSaving to HDF5...")
h5_path = os.path.join(OUTDIR, "prepared_data.h5")

with h5py.File(h5_path, "w") as h5:
    # Expression matrices (genes x samples, float32)
    h5.create_dataset("rank_expression", data=rank_df.values.astype(np.float32),
                       compression="gzip", compression_opts=4)
    h5.create_dataset("gene_names", data=np.array(rank_df.index.tolist(), dtype="S"))
    h5.create_dataset("sample_ids", data=np.array(rank_df.columns.tolist(), dtype="S"))

    if zscore_df is not None:
        h5.create_dataset("zscore_expression", data=zscore_df.values.astype(np.float32),
                           compression="gzip", compression_opts=4)

    if raw_df is not None:
        h5.create_dataset("raw_logcpm", data=raw_df.values.astype(np.float32),
                           compression="gzip", compression_opts=4)

    # Metadata (as separate arrays for each column)
    meta_grp = h5.create_group("metadata")
    for col in ["dataset", "condition", "group_binary", "sex",
                 "diagnosis_harmonized", "loco_fold_fibrosis", "loco_fold_nas"]:
        if col in meta_sub.columns:
            vals = meta_sub[col].fillna("NA").values.astype(str)
            meta_grp.create_dataset(col, data=np.array(vals, dtype="S"))

    for col in ["age", "fibrosis_stage", "nas_score", "nas_group", "nas_group4",
                 "fib_stage", "fib_ge3", "nas_ge5", "is_disease", "severity4"]:
        if col in meta_sub.columns:
            vals = meta_sub[col].fillna(-1).values.astype(np.float32)
            meta_grp.create_dataset(col, data=vals)

    # Feature candidates
    fc_grp = h5.create_group("feature_candidates")
    fc_grp.create_dataset("gene", data=np.array(feature_candidates["gene"].tolist(), dtype="S"))
    if "max_t_abs" in feature_candidates.columns:
        fc_grp.create_dataset("max_t_abs", data=feature_candidates["max_t_abs"].values.astype(np.float32))
    if "n_sources" in feature_candidates.columns:
        fc_grp.create_dataset("n_sources", data=feature_candidates["n_sources"].values.astype(np.int32))

    # Multi-evidence atlas features (for Plan 3 fusion)
    if atlas_df is not None:
        atlas_grp = h5.create_group("atlas_features")
        # Match atlas genes to our gene list
        # Atlas uses "human_symbol" as gene identifier
        gene_col = "human_symbol" if "human_symbol" in atlas_df.columns else "gene"
        atlas_matched = atlas_df[atlas_df[gene_col].isin(rank_df.index)]
        atlas_grp.create_dataset("gene", data=np.array(atlas_matched[gene_col].tolist(), dtype="S"))
        # Select numeric columns only
        numeric_cols = atlas_matched.select_dtypes(include=[np.number]).columns.tolist()
        atlas_numeric = atlas_matched[numeric_cols].fillna(0).values.astype(np.float32)
        atlas_grp.create_dataset("features", data=atlas_numeric, compression="gzip")
        atlas_grp.create_dataset("feature_names", data=np.array(numeric_cols, dtype="S"))
        print(f"  Atlas features: {atlas_numeric.shape}")

h5_size_mb = os.path.getsize(h5_path) / 1e6
print(f"\nHDF5 saved: {h5_path} ({h5_size_mb:.1f} MB)")

# --- Save metadata as CSV too (easier for R scripts to read) ---
meta_out_path = os.path.join(OUTDIR, "modeling_metadata.csv")
meta_sub.to_csv(meta_out_path, index=False)
print(f"Metadata CSV saved: {meta_out_path}")

# --- Summary ---
print("\n=== Summary ===")
print(f"Genes: {rank_df.shape[0]}")
print(f"Samples: {rank_df.shape[1]}")
print(f"Fibrosis-annotated: {(meta_sub['fib_stage'] >= 0).sum()}")
print(f"NAS-annotated: {(meta_sub['nas_group'] >= 0).sum()}")
print(f"Both NAS+fibrosis: {((meta_sub['fib_stage'] >= 0) & (meta_sub['nas_group'] >= 0)).sum()}")
print(f"Unannotated (no staging): {((meta_sub['fib_stage'] < 0) & (meta_sub['nas_group'] < 0)).sum()}")
print(f"LOCO fibrosis folds: {len(fib_fold_counts)}")
print(f"LOCO NAS folds: {len(nas_fold_counts)}")

print(f"\n=== 62_prepare_python_data.py completed ===")
