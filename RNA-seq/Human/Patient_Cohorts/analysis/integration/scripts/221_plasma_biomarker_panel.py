#!/usr/bin/env python3
"""
221_plasma_biomarker_panel.py
Tissue-Informed Plasma Biomarker Panel with Random Baselines

Tests whether tissue DEG-informed protein selection outperforms random
protein selection for plasma-based fibrosis prediction (F3/F4 vs F0-2).

Panels:
  A: All 1,461 Olink proteins (full discovery)
  B: Tissue-informed — proteins whose gene is a tissue DEG (812 genes)
  C: Top-20 tissue-informed (ranked by |dream_logFC| among tissue DEGs in plasma)
  D: Random-20 baseline (100 draws from all 1,461 proteins)

CV strategy:
  - 5-fold stratified CV (primary)
  - Leave-one-etiology-out (LOEO) CV when etiology labels are available

Inputs:
  - olink.qc.finished.mendeley.data.txt  (1,461 proteins × 436 subjects, NPX)
  - gse276114_disease_metadata.csv        (177 subjects with fibrosis staging)
  - tissue_plasma_bridge.csv              (33,943 genes, is_deg + in_plasma flags)

Outputs (in results/multiprogram/):
  - plasma_panel_results.csv        — per-panel AUROC (5-fold + LOEO)
  - plasma_random_baselines.csv     — 100 random draws × AUROC
  - plasma_panel_summary.csv        — tissue-informed AUROC, random mean, p-value
  - plasma_selected_proteins.csv    — proteins in each panel + tissue evidence

Usage:
  micromamba run -n spatial python 221_plasma_biomarker_panel.py

SLURM: cpu, 8 CPUs, 64GB, 48h
"""

import os
import sys
import time
import logging
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from scipy.stats import rankdata

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================
# Logging
# ============================================================
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ============================================================
# Paths
# ============================================================
BASE = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)
INT = BASE / "RNA-seq/Human/Patient_Cohorts/analysis/integration"
OUTDIR = INT / "results/multiprogram"
OLINK_FILE = BASE / "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt"
META_FILE = BASE / "Analysis/Proteomics/results/gse276114_disease_metadata.csv"
BRIDGE_FILE = INT / "results/staging_classifier/tissue_plasma_bridge.csv"

OUTDIR.mkdir(parents=True, exist_ok=True)

# ============================================================
# Elastic net model — sklearn equivalent of glmnet alpha=0.5
# ============================================================
ENET_PARAMS = dict(
    penalty="elasticnet",
    solver="saga",
    l1_ratio=0.5,
    max_iter=2000,
    C=1.0,
    class_weight="balanced",
    random_state=42,
    n_jobs=1,
)

N_RANDOM_DRAWS = 100
PANEL_SIZE = 20
N_FOLDS = 5
RANDOM_SEED = 42


# ============================================================
# Helper: run 5-fold stratified CV, return per-fold + mean AUROC
# ============================================================
def cv_auroc(X: np.ndarray, y: np.ndarray, seed: int = RANDOM_SEED) -> dict:
    """5-fold stratified CV; returns mean AUROC and per-fold list."""
    skf = StratifiedKFold(n_splits=N_FOLDS, shuffle=True, random_state=seed)
    fold_aucs = []
    for train_idx, test_idx in skf.split(X, y):
        pipe = Pipeline(
            [
                ("scaler", StandardScaler()),
                ("clf", LogisticRegression(**ENET_PARAMS)),
            ]
        )
        pipe.fit(X[train_idx], y[train_idx])
        proba = pipe.predict_proba(X[test_idx])[:, 1]
        if len(np.unique(y[test_idx])) < 2:
            continue  # skip degenerate fold
        fold_aucs.append(roc_auc_score(y[test_idx], proba))
    if len(fold_aucs) == 0:
        return {"auroc_cv_mean": np.nan, "auroc_cv_sd": np.nan, "n_folds_cv": 0}
    return {
        "auroc_cv_mean": float(np.mean(fold_aucs)),
        "auroc_cv_sd": float(np.std(fold_aucs)),
        "n_folds_cv": len(fold_aucs),
    }


# ============================================================
# Helper: LOEO CV — hold out one etiology at a time
# ============================================================
def loeo_auroc(
    X: np.ndarray, y: np.ndarray, etiology: np.ndarray
) -> dict:
    """Leave-one-etiology-out CV; returns mean AUROC across held-out groups."""
    groups = np.unique(etiology)
    all_proba = np.full(len(y), np.nan)
    held_out_aucs = {}

    for grp in groups:
        test_mask = etiology == grp
        train_mask = ~test_mask

        if train_mask.sum() < 10:
            log.warning("  LOEO: group %s has too few training samples, skipping", grp)
            continue
        if len(np.unique(y[train_mask])) < 2:
            log.warning("  LOEO: group %s train set is single-class, skipping", grp)
            continue

        pipe = Pipeline(
            [
                ("scaler", StandardScaler()),
                ("clf", LogisticRegression(**ENET_PARAMS)),
            ]
        )
        pipe.fit(X[train_mask], y[train_mask])
        proba = pipe.predict_proba(X[test_mask])[:, 1]
        all_proba[test_mask] = proba

        if len(np.unique(y[test_mask])) >= 2:
            held_out_aucs[grp] = float(roc_auc_score(y[test_mask], proba))
        log.info(
            "  LOEO held-out %s (n=%d): AUROC=%.3f",
            grp,
            test_mask.sum(),
            held_out_aucs.get(grp, np.nan),
        )

    # Pooled AUROC over all held-out subjects
    valid_mask = ~np.isnan(all_proba)
    if valid_mask.sum() < 10 or len(np.unique(y[valid_mask])) < 2:
        return {"auroc_loeo_pooled": np.nan, "n_loeo_subjects": int(valid_mask.sum())}

    pooled = float(roc_auc_score(y[valid_mask], all_proba[valid_mask]))
    return {
        "auroc_loeo_pooled": pooled,
        "n_loeo_subjects": int(valid_mask.sum()),
        "loeo_per_group": held_out_aucs,
    }


# ============================================================
# STEP 1: Load Olink data
# ============================================================
log.info("=== STEP 1: Loading Olink data ===")
olink_raw = pd.read_csv(OLINK_FILE, sep="\t")
log.info("  Raw Olink: %d proteins × %d subjects", olink_raw.shape[0], olink_raw.shape[1] - 1)

# Transpose to subjects × proteins
protein_names = olink_raw["Assay"].values
subject_cols = [c for c in olink_raw.columns if c != "Assay"]
npx = olink_raw[subject_cols].values.T.astype(float)  # (n_subjects, n_proteins)
subject_labels = subject_cols  # "Subject 1", ..., "Subject 436"
log.info("  Transposed NPX: %d subjects × %d proteins", npx.shape[0], npx.shape[1])

# Extract integer subject number
subject_nums = np.array(
    [int(s.replace("Subject", "").strip()) for s in subject_labels]
)

# Impute NAs with column median
na_count = np.sum(np.isnan(npx))
if na_count > 0:
    col_medians = np.nanmedian(npx, axis=0)
    for j in range(npx.shape[1]):
        mask = np.isnan(npx[:, j])
        if mask.any():
            npx[mask, j] = col_medians[j]
    log.info("  Imputed %d NAs with column medians", na_count)
else:
    log.info("  No NAs found")

# ============================================================
# STEP 2: Load metadata and match subjects
# ============================================================
log.info("=== STEP 2: Loading metadata ===")
meta = pd.read_csv(META_FILE)
log.info("  Metadata: %d rows", len(meta))
log.info("  Etiology: %s", dict(meta["disease"].value_counts()))
log.info("  Fibrosis stage: %s", dict(meta["disease_group"].value_counts()))

# Map Olink subject number → fibrosis stage
meta_map = meta[["sample_number", "disease_group", "disease"]].copy()
meta_map = meta_map.dropna(subset=["disease_group"])

subj_df = pd.DataFrame({"subject_num": subject_nums, "array_idx": np.arange(len(subject_nums))})
merged = subj_df.merge(meta_map, left_on="subject_num", right_on="sample_number", how="inner")
log.info("  Matched %d subjects with fibrosis staging", len(merged))

# Build arrays for matched subjects
matched_idx = merged["array_idx"].values
npx_matched = npx[matched_idx]  # (n_matched, n_proteins)
# Binary target: advanced (F3/F4) = 1, early (F0-2) = 0
y = (merged["disease_group"].isin(["F3", "F4"])).astype(int).values
etiology = merged["disease"].values

log.info(
    "  Advanced (F3/F4): %d  |  Early (F0-2): %d",
    int(y.sum()),
    int((y == 0).sum()),
)
log.info("  Etiology breakdown:")
for eti in np.unique(etiology):
    mask = etiology == eti
    log.info(
        "    %s: n=%d, advanced=%d, early=%d",
        eti,
        int(mask.sum()),
        int(y[mask].sum()),
        int((y[mask] == 0).sum()),
    )

# ============================================================
# STEP 3: Load tissue-plasma bridge and define panels
# ============================================================
log.info("=== STEP 3: Defining protein panels ===")
bridge = pd.read_csv(BRIDGE_FILE)
log.info("  Bridge: %d genes", len(bridge))

# Tissue-informed: is_deg & in_plasma
tissue_degs_in_plasma = bridge[(bridge["is_deg"] == True) & (bridge["in_plasma"] == True)].copy()
tissue_degs_in_plasma = tissue_degs_in_plasma.sort_values("dream_logFC", key=abs, ascending=False)
log.info("  Tissue DEGs also in plasma (Panel B candidates): %d", len(tissue_degs_in_plasma))

# Map gene symbols to Olink column indices
protein_name_to_idx = {p: i for i, p in enumerate(protein_names)}

# Panel A: all Olink proteins
panel_a_idx = np.arange(len(protein_names))

# Panel B: tissue-informed DEGs in plasma
panel_b_symbols = [
    sym for sym in tissue_degs_in_plasma["human_symbol"].values
    if sym in protein_name_to_idx
]
panel_b_idx = np.array([protein_name_to_idx[s] for s in panel_b_symbols])
log.info("  Panel A (all Olink): %d proteins", len(panel_a_idx))
log.info("  Panel B (tissue-informed DEGs): %d proteins", len(panel_b_idx))

# Panel C: top-20 tissue-informed by |dream_logFC|
panel_c_symbols = panel_b_symbols[:PANEL_SIZE]
panel_c_idx = np.array([protein_name_to_idx[s] for s in panel_c_symbols])
log.info("  Panel C (top-%d tissue-informed): %s", PANEL_SIZE, panel_c_symbols)

# Panel D: 100 random draws of 20 proteins
rng = np.random.default_rng(RANDOM_SEED)
all_protein_idx = np.arange(len(protein_names))
random_panels = [
    rng.choice(all_protein_idx, size=PANEL_SIZE, replace=False)
    for _ in range(N_RANDOM_DRAWS)
]

# ============================================================
# STEP 4: Fit panels and evaluate
# ============================================================
log.info("=== STEP 4: Fitting panels ===")
results_rows = []
panel_configs = [
    ("A_full_olink", panel_a_idx, "All 1461 Olink proteins"),
    ("B_tissue_informed_all", panel_b_idx, f"Tissue DEGs in plasma (n={len(panel_b_idx)})"),
    ("C_tissue_informed_top20", panel_c_idx, f"Top-{PANEL_SIZE} tissue-informed by |logFC|"),
]

for panel_name, pidx, description in panel_configs:
    if len(pidx) == 0:
        log.warning("  %s: no proteins — skipping", panel_name)
        continue

    X_panel = npx_matched[:, pidx]
    log.info("  Panel %s (%d proteins)...", panel_name, len(pidx))

    t0 = time.time()
    cv_res = cv_auroc(X_panel, y)
    loeo_res = loeo_auroc(X_panel, y, etiology)
    elapsed = time.time() - t0

    row = {
        "panel": panel_name,
        "description": description,
        "n_proteins": len(pidx),
        **cv_res,
        "auroc_loeo_pooled": loeo_res.get("auroc_loeo_pooled", np.nan),
        "n_loeo_subjects": loeo_res.get("n_loeo_subjects", np.nan),
        "elapsed_s": round(elapsed, 1),
    }
    results_rows.append(row)
    log.info(
        "    CV AUROC=%.3f ± %.3f  |  LOEO AUROC=%.3f",
        cv_res["auroc_cv_mean"],
        cv_res["auroc_cv_sd"],
        loeo_res.get("auroc_loeo_pooled", np.nan),
    )

# ============================================================
# STEP 5: Random baselines for Panel C (top-20)
# ============================================================
log.info("=== STEP 5: Random baselines (n=%d draws, size=%d) ===", N_RANDOM_DRAWS, PANEL_SIZE)
random_rows = []
random_cv_aurocs = []
random_loeo_aurocs = []

for i, pidx in enumerate(random_panels):
    X_rand = npx_matched[:, pidx]
    cv_res = cv_auroc(X_rand, y, seed=RANDOM_SEED + i)
    loeo_res = loeo_auroc(X_rand, y, etiology)

    cv_auc = cv_res["auroc_cv_mean"]
    loeo_auc = loeo_res.get("auroc_loeo_pooled", np.nan)
    random_cv_aurocs.append(cv_auc)
    if not np.isnan(loeo_auc):
        random_loeo_aurocs.append(loeo_auc)

    random_rows.append(
        {
            "draw": i + 1,
            "proteins": ",".join(protein_names[pidx]),
            "auroc_cv_mean": cv_auc,
            "auroc_cv_sd": cv_res["auroc_cv_sd"],
            "auroc_loeo_pooled": loeo_auc,
        }
    )
    if (i + 1) % 20 == 0:
        log.info(
            "  Draw %d/%d complete  (cv mean so far=%.3f)",
            i + 1,
            N_RANDOM_DRAWS,
            float(np.nanmean(random_cv_aurocs)),
        )

# Tissue-informed Panel C AUROC (already computed above)
panel_c_row = next((r for r in results_rows if r["panel"] == "C_tissue_informed_top20"), None)

# Empirical p-values
def empirical_pval(observed: float, null: list) -> float:
    """Fraction of null draws >= observed (one-tailed, tissue beats random)."""
    null_arr = np.array([x for x in null if not np.isnan(x)])
    if len(null_arr) == 0 or np.isnan(observed):
        return np.nan
    return float(np.mean(null_arr >= observed))

summary_rows = []
if panel_c_row is not None:
    obs_cv = panel_c_row["auroc_cv_mean"]
    obs_loeo = panel_c_row["auroc_loeo_pooled"]
    p_cv = empirical_pval(obs_cv, random_cv_aurocs)
    p_loeo = empirical_pval(obs_loeo, random_loeo_aurocs)
    summary_rows.append(
        {
            "metric": "5fold_CV",
            "tissue_informed_auroc": obs_cv,
            "random_mean_auroc": float(np.nanmean(random_cv_aurocs)),
            "random_sd_auroc": float(np.nanstd(random_cv_aurocs)),
            "empirical_pval": p_cv,
            "n_random_draws": N_RANDOM_DRAWS,
        }
    )
    summary_rows.append(
        {
            "metric": "LOEO_pooled",
            "tissue_informed_auroc": obs_loeo,
            "random_mean_auroc": float(np.nanmean(random_loeo_aurocs)) if random_loeo_aurocs else np.nan,
            "random_sd_auroc": float(np.nanstd(random_loeo_aurocs)) if random_loeo_aurocs else np.nan,
            "empirical_pval": p_loeo,
            "n_random_draws": len(random_loeo_aurocs),
        }
    )
    log.info(
        "  Panel C vs random — CV: obs=%.3f, null_mean=%.3f, p=%.3f",
        obs_cv,
        float(np.nanmean(random_cv_aurocs)),
        p_cv,
    )
    log.info(
        "  Panel C vs random — LOEO: obs=%.3f, null_mean=%.3f, p=%.3f",
        obs_loeo,
        float(np.nanmean(random_loeo_aurocs)) if random_loeo_aurocs else np.nan,
        p_loeo,
    )

# ============================================================
# STEP 6: Build plasma_selected_proteins.csv
# ============================================================
log.info("=== STEP 6: Building protein selection table ===")
protein_table_rows = []

for i, pname in enumerate(protein_names):
    in_panel_a = True
    in_panel_b = pname in set(panel_b_symbols)
    in_panel_c = pname in set(panel_c_symbols)

    bridge_row = bridge[bridge["human_symbol"] == pname]
    dream_logfc = float(bridge_row["dream_logFC"].iloc[0]) if len(bridge_row) > 0 else np.nan
    dream_padj = float(bridge_row["dream_padj"].iloc[0]) if len(bridge_row) > 0 else np.nan
    is_deg = bool(bridge_row["is_deg"].iloc[0]) if len(bridge_row) > 0 else False
    classification = bridge_row["classification"].iloc[0] if len(bridge_row) > 0 else "unknown"

    protein_table_rows.append(
        {
            "protein": pname,
            "in_panel_a": in_panel_a,
            "in_panel_b": in_panel_b,
            "in_panel_c": in_panel_c,
            "is_tissue_deg": is_deg,
            "dream_logFC": dream_logfc,
            "dream_padj": dream_padj,
            "tissue_classification": classification,
        }
    )

protein_df = pd.DataFrame(protein_table_rows)

# ============================================================
# STEP 7: Save outputs
# ============================================================
log.info("=== STEP 7: Saving outputs ===")

results_df = pd.DataFrame(results_rows)
results_df.to_csv(OUTDIR / "plasma_panel_results.csv", index=False)
log.info("  Saved plasma_panel_results.csv (%d rows)", len(results_df))

random_df = pd.DataFrame(random_rows)
random_df.to_csv(OUTDIR / "plasma_random_baselines.csv", index=False)
log.info("  Saved plasma_random_baselines.csv (%d rows)", len(random_df))

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(OUTDIR / "plasma_panel_summary.csv", index=False)
log.info("  Saved plasma_panel_summary.csv (%d rows)", len(summary_df))

protein_df.to_csv(OUTDIR / "plasma_selected_proteins.csv", index=False)
log.info("  Saved plasma_selected_proteins.csv (%d rows)", len(protein_df))

# ============================================================
# Final summary
# ============================================================
log.info("")
log.info("=== FINAL SUMMARY ===")
log.info("Panel results:")
for _, row in results_df.iterrows():
    log.info(
        "  %-35s  n=%4d  CV=%.3f±%.3f  LOEO=%.3f",
        row["panel"],
        int(row["n_proteins"]),
        row["auroc_cv_mean"],
        row["auroc_cv_sd"],
        row["auroc_loeo_pooled"] if not pd.isna(row["auroc_loeo_pooled"]) else float("nan"),
    )

if summary_rows:
    log.info("Random baseline comparison (Panel C top-%d):", PANEL_SIZE)
    for _, row in summary_df.iterrows():
        log.info(
            "  %s: tissue=%.3f  random_mean=%.3f  p=%.3f",
            row["metric"],
            row["tissue_informed_auroc"],
            row["random_mean_auroc"],
            row["empirical_pval"],
        )

log.info("Done. All outputs in: %s", OUTDIR)
