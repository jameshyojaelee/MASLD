#!/usr/bin/env python3
"""192_genetics_only_prediction.py — COLOC gene panels for multi-task ordinal prediction.

Tests whether genetically validated genes (COLOC) predict each clinical axis.
Three panels: cross-ancestry (66 genes), EUR-only (309), all COLOC (375).
Compares each to matched-size random gene baselines (100 draws).

SLURM:
  sbatch --job-name=192_genetics --partition=cpu --cpus-per-task=8 --mem=64G --time=48:00:00 \\
         --output=logs/192_%j.out --error=logs/192_%j.err \\
         --wrap="bash -c 'eval \"\\$(micromamba shell hook --shell bash)\" && \\
                 micromamba activate spatial && \\
                 cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts && \\
                 python 192_genetics_only_prediction.py'"
"""

import os
import sys
import numpy as np
import pandas as pd
import h5py
import logging
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import cohen_kappa_score, mean_absolute_error, accuracy_score

SEED = 42
np.random.seed(SEED)

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
STAGING = os.path.join(INT, "results/staging_classifier")
OUTDIR = os.path.join(INT, "results/multiprogram")
os.makedirs(OUTDIR, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

MISSING = -1
NAS_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066"]
FIB_DATASETS = NAS_DATASETS + ["GSE240729"]

TASK_CFG = {
    "fibrosis": (5, True, "fib_stage", "loco_fold_fibrosis", FIB_DATASETS),
    "nas_composite": (9, True, "nas_score", "loco_fold_nas", NAS_DATASETS),
    "severity": (4, True, "severity4", "loco_fold_fibrosis", FIB_DATASETS),
}

log.info("=== 192: Genetics-Only Prediction ===")

# Load data
meta = pd.read_csv(os.path.join(STAGING, "modeling_metadata.csv"))
sample_ids = meta["sample_id"].values
datasets = meta["dataset"].values

with h5py.File(os.path.join(STAGING, "prepared_data.h5"), "r") as h5:
    expr_samples = [s.decode() if isinstance(s, bytes) else s for s in h5["sample_ids"][:]]
    gene_names = [g.decode() if isinstance(g, bytes) else g for g in h5["gene_names"][:]]
    expr_mat = h5["zscore_expression"][:].T  # (genes x samples) -> (samples x genes)

expr_order = {s: i for i, s in enumerate(expr_samples)}
expr_idx = [expr_order[s] for s in sample_ids if s in expr_order]
expr_mat = expr_mat[expr_idx]

# Load COLOC results
coloc_df = pd.read_csv(os.path.join(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv"))

# Define panels — match COLOC Ensembl IDs to expression matrix (versioned Ensembl)
gene_base = {g.split(".")[0]: i for i, g in enumerate(gene_names)}
all_coloc_ensembl = coloc_df[coloc_df["coloc_best_pp4"] > 0.5]["ensembl"].dropna().unique()
high_coloc_ensembl = coloc_df[coloc_df["coloc_best_pp4"] > 0.8]["ensembl"].dropna().unique()

panels = {
    "coloc_high_confidence": ([gene_base[e] for e in high_coloc_ensembl if e in gene_base], 0.8),
    "coloc_all": ([gene_base[e] for e in all_coloc_ensembl if e in gene_base], 0.5),
}
# panels now map to expression matrix column indices, not gene names

for name, (idx, thresh) in panels.items():
    log.info(f"Panel '{name}': {len(idx)} genes (PP4>{thresh})")

# Labels
labels = {}
for task, (n_cls, is_ord, col, fold_col, ds_list) in TASK_CFG.items():
    labels[task] = meta[col].fillna(MISSING).astype(int).values


def run_loco(X, task_name, panel_name):
    """Run LOCO-CV elastic net for one panel × one target."""
    n_cls, is_ord, col, fold_col, ds_list = TASK_CFG[task_name]
    y = labels[task_name]
    results = []

    for fold_ds in ds_list:
        ds_mask = np.isin(datasets, ds_list)
        test_mask = (datasets == fold_ds) & (y != MISSING) & ds_mask
        train_mask = (datasets != fold_ds) & (y != MISSING) & ds_mask

        if test_mask.sum() == 0 or train_mask.sum() == 0:
            continue

        X_tr, X_te = X[train_mask], X[test_mask]
        y_tr, y_te = y[train_mask], y[test_mask]

        scaler = StandardScaler()
        X_tr = scaler.fit_transform(X_tr)
        X_te = scaler.transform(X_te)

        if is_ord and n_cls > 2:
            preds = np.zeros(len(y_te))
            for threshold in range(1, n_cls):
                y_bin = (y_tr >= threshold).astype(int)
                if len(np.unique(y_bin)) < 2:
                    continue
                lr = LogisticRegression(
                    penalty="elasticnet", solver="saga", l1_ratio=0.5,
                    max_iter=2000, C=1.0, random_state=SEED
                )
                lr.fit(X_tr, y_bin)
                preds += lr.predict_proba(X_te)[:, 1]
            preds = np.clip(np.round(preds), 0, n_cls - 1).astype(int)
        else:
            lr = LogisticRegression(
                penalty="elasticnet", solver="saga", l1_ratio=0.5,
                max_iter=2000, C=1.0, random_state=SEED, multi_class="multinomial"
            )
            lr.fit(X_tr, y_tr)
            preds = lr.predict(X_te)

        qwk = cohen_kappa_score(y_te, preds, weights="quadratic")
        mae = mean_absolute_error(y_te, preds)
        results.append({"task": task_name, "panel": panel_name, "fold": fold_ds,
                        "n_genes": X_tr.shape[1], "qwk": qwk, "mae": mae})
    return results


# Run COLOC panels
all_results = []
for panel_name, (gene_idx, _) in panels.items():
    X_panel = expr_mat[:, gene_idx]
    for task_name in TASK_CFG:
        results = run_loco(X_panel, task_name, panel_name)
        all_results.extend(results)
        mean_qwk = np.mean([r["qwk"] for r in results]) if results else np.nan
        log.info(f"  {panel_name} × {task_name}: QWK={mean_qwk:.3f}")

# Random gene baselines (matched size per panel)
log.info("Running random gene baselines...")
rng = np.random.RandomState(SEED)
baseline_results = []

for panel_name, (gene_idx, _) in panels.items():
    n_genes = len(gene_idx)
    for draw in range(100):
        rand_idx = rng.choice(expr_mat.shape[1], size=n_genes, replace=False)
        X_rand = expr_mat[:, rand_idx]
        for task_name in TASK_CFG:
            results = run_loco(X_rand, task_name, f"random_{panel_name}")
            mean_qwk = np.mean([r["qwk"] for r in results]) if results else np.nan
            baseline_results.append({"task": task_name, "panel": panel_name,
                                     "draw": draw, "mean_qwk": mean_qwk})
        if (draw + 1) % 25 == 0:
            log.info(f"  Random baseline {panel_name}: {draw + 1}/100")

# Save results
results_df = pd.DataFrame(all_results)
results_df.to_csv(os.path.join(OUTDIR, "genetics_only_results.csv"), index=False)

baselines_df = pd.DataFrame(baseline_results)
baselines_df.to_csv(os.path.join(OUTDIR, "genetics_random_baselines.csv"), index=False)

# Summary with p-values
summary_rows = []
for panel_name, (gene_idx, _) in panels.items():
    for task_name in TASK_CFG:
        panel_qwk = results_df[(results_df["panel"] == panel_name) & (results_df["task"] == task_name)]["qwk"].mean()
        rand_qwks = baselines_df[(baselines_df["panel"] == panel_name) & (baselines_df["task"] == task_name)]["mean_qwk"].values
        p_val = (np.sum(rand_qwks >= panel_qwk) + 1) / (len(rand_qwks) + 1) if len(rand_qwks) > 0 else np.nan
        summary_rows.append({"panel": panel_name, "task": task_name,
                             "n_genes": len(gene_idx), "mean_qwk": panel_qwk,
                             "random_mean_qwk": rand_qwks.mean() if len(rand_qwks) > 0 else np.nan,
                             "p_value": p_val})
        log.info(f"  {panel_name} × {task_name}: QWK={panel_qwk:.3f}, random={rand_qwks.mean():.3f}, p={p_val:.3f}")

pd.DataFrame(summary_rows).to_csv(os.path.join(OUTDIR, "genetics_only_summary.csv"), index=False)
log.info("Done.")
