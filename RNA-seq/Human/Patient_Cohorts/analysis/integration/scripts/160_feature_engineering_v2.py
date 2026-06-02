#!/usr/bin/env python3
"""
160_feature_engineering_v2.py — Rigorous Feature Engineering v2
---------------------------------------------------------------
Patched version of the feature matrix for the v2 prognosis pipeline:
  1. Copies the existing feature matrix from results/prognosis/
  2. Removes any velocity columns (if present)
  3. Creates a clean label file with fibrosis>=3 as primary, LOCO-NMF S2,
     and metadata covariates -- NO CPS
  4. Saves feature_matrix_v2.csv and labels_v2.csv to results/prognosis_v2/

Input:
  - results/prognosis/prognosis_feature_matrix.csv         (253 features x 1,444 samples)
  - results/staging_classifier/modeling_metadata.csv        (fibrosis, folds, covariates)
  - results/prognosis_v2/loco_nmf_labels/loco_nmf_combined.csv  (LOCO-NMF + global S1/S2)

Output (to results/prognosis_v2/):
  - feature_matrix_v2.csv   (feature matrix, velocity columns dropped if any)
  - labels_v2.csv           (sample_id + labels + covariates, no CPS)

SLURM: cpu, 4 CPUs, 8G RAM, 48h
Env:   micromamba activate spatial
"""

import logging
import os
import sys
import time

import numpy as np
import pandas as pd

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ── Paths ────────────────────────────────────────────────────────────────────
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INTEG = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")

FEATURE_MATRIX = os.path.join(INTEG, "results/prognosis/prognosis_feature_matrix.csv")
MODELING_META = os.path.join(INTEG, "results/staging_classifier/modeling_metadata.csv")
LOCO_NMF = os.path.join(
    INTEG, "results/prognosis_v2/loco_nmf_labels/loco_nmf_combined.csv"
)

OUT_DIR = os.path.join(INTEG, "results/prognosis_v2")


def main():
    t0 = time.time()
    log.info("=== Script 160: Rigorous Feature Engineering v2 ===")

    # ── Create output directory ──────────────────────────────────────────────
    os.makedirs(OUT_DIR, exist_ok=True)

    # ── Step 1: Load and copy feature matrix ─────────────────────────────────
    log.info("Loading feature matrix from %s", FEATURE_MATRIX)
    if not os.path.isfile(FEATURE_MATRIX):
        log.error("Feature matrix not found: %s", FEATURE_MATRIX)
        sys.exit(1)

    feat = pd.read_csv(FEATURE_MATRIX)
    log.info("  Loaded: %d samples x %d columns", feat.shape[0], feat.shape[1])

    # ── Step 2: Remove velocity columns (if any) ────────────────────────────
    velocity_cols = [c for c in feat.columns if "velocity" in c.lower()]
    if velocity_cols:
        log.info("  Dropping %d velocity column(s): %s", len(velocity_cols), velocity_cols)
        feat = feat.drop(columns=velocity_cols)
    else:
        log.info("  No velocity columns found -- nothing to drop")

    n_features = feat.shape[1] - 1  # exclude sample_id
    log.info("  Feature matrix v2: %d samples x %d features", feat.shape[0], n_features)

    feat_out = os.path.join(OUT_DIR, "feature_matrix_v2.csv")
    feat.to_csv(feat_out, index=False)
    log.info("  Saved: %s", feat_out)

    # ── Step 3: Build labels_v2.csv ──────────────────────────────────────────
    log.info("Loading modeling metadata from %s", MODELING_META)
    meta = pd.read_csv(MODELING_META)
    log.info("  Loaded: %d samples x %d columns", meta.shape[0], meta.shape[1])

    # Replace -1 sentinel with NaN for fibrosis fields
    meta["fib_stage"] = meta["fib_stage"].replace(-1, np.nan)
    meta["fib_ge3"] = meta["fib_ge3"].replace(-1, np.nan)

    # Load LOCO-NMF labels
    log.info("Loading LOCO-NMF labels from %s", LOCO_NMF)
    if not os.path.isfile(LOCO_NMF):
        log.warning("LOCO-NMF labels file not found: %s", LOCO_NMF)
        log.warning("  loco_s2 and global_s2 will be set to NaN for all samples")
        nmf = pd.DataFrame({"sample_id": meta["sample_id"]})
        nmf["loco_nmf_subtype"] = np.nan
        nmf["global_nmf_subtype"] = np.nan
    else:
        nmf = pd.read_csv(LOCO_NMF)
        log.info("  Loaded: %d samples with LOCO-NMF labels", nmf.shape[0])

    # Merge NMF onto metadata
    labels = meta[["sample_id"]].copy()

    # fib_ge3: 1 if fibrosis_stage >= 3, 0 if < 3, NaN if missing
    labels["fib_ge3"] = meta["fib_ge3"].copy()

    # fib_ordinal: 0-4 fibrosis stage, NaN if missing
    labels["fib_ordinal"] = meta["fib_stage"].copy()

    # LOCO-NMF S2: S1=0, S2=1
    nmf_merge = nmf[["sample_id", "loco_nmf_subtype", "global_nmf_subtype"]].copy()
    labels = labels.merge(nmf_merge, on="sample_id", how="left")

    labels["loco_s2"] = labels["loco_nmf_subtype"].map({"S1": 0, "S2": 1})
    labels["global_s2"] = labels["global_nmf_subtype"].map({"S1": 0, "S2": 1})
    labels = labels.drop(columns=["loco_nmf_subtype", "global_nmf_subtype"])

    # Covariates and fold info
    labels["loco_fold_fibrosis"] = meta["loco_fold_fibrosis"]
    labels["dataset"] = meta["dataset"]
    labels["sex"] = meta["sex"]
    labels["age"] = meta["age"]

    log.info("  Labels v2 assembled: %d samples x %d columns", labels.shape[0], labels.shape[1])

    # ── Summary statistics ───────────────────────────────────────────────────
    log.info("  Label summary:")
    log.info("    fib_ge3:  0=%d, 1=%d, NaN=%d",
             (labels["fib_ge3"] == 0).sum(),
             (labels["fib_ge3"] == 1).sum(),
             labels["fib_ge3"].isna().sum())
    log.info("    fib_ordinal distribution:")
    for stage in sorted(labels["fib_ordinal"].dropna().unique()):
        log.info("      F%d: %d", int(stage), (labels["fib_ordinal"] == stage).sum())
    log.info("      NaN: %d", labels["fib_ordinal"].isna().sum())
    log.info("    loco_s2:  0 (S1)=%d, 1 (S2)=%d, NaN=%d",
             (labels["loco_s2"] == 0).sum(),
             (labels["loco_s2"] == 1).sum(),
             labels["loco_s2"].isna().sum())
    log.info("    global_s2: 0 (S1)=%d, 1 (S2)=%d, NaN=%d",
             (labels["global_s2"] == 0).sum(),
             (labels["global_s2"] == 1).sum(),
             labels["global_s2"].isna().sum())

    labels_out = os.path.join(OUT_DIR, "labels_v2.csv")
    labels.to_csv(labels_out, index=False)
    log.info("  Saved: %s", labels_out)

    # ── Done ─────────────────────────────────────────────────────────────────
    elapsed = time.time() - t0
    log.info("=== Done in %.1f seconds ===", elapsed)
    log.info("Outputs:")
    log.info("  %s (%d x %d)", feat_out, feat.shape[0], feat.shape[1])
    log.info("  %s (%d x %d)", labels_out, labels.shape[0], labels.shape[1])


if __name__ == "__main__":
    main()
