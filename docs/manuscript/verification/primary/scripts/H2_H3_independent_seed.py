#!/usr/bin/env python3
"""H2/H3 Independent Seed Verification
===========================================
Re-runs the multiprogram random-gene baseline with an independent random seed
(SEED=2026 instead of 42) and an expanded draw count (N=1000) to address the
pre-registered concern that the 100-draw p=0.010 is right at the boundary of
what a 100-draw test can resolve.

OUTPUT: multiprogram_random_baselines_seed2026.csv
        h2_h3_independent_seed_summary.csv
"""
import os
import sys
import time
import numpy as np
import pandas as pd
import h5py
from pathlib import Path

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import cohen_kappa_score, mean_absolute_error, accuracy_score

import warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# Use DIFFERENT seed
SEED = 2026
MODEL_SEED = 42  # keep model seed for fair comparison, only shuffle gene draws
np.random.seed(SEED)

BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR = os.path.join(INT, "results")
STAGING = os.path.join(RDIR, "staging_classifier")
OUTDIR = os.path.join(RDIR, "multiprogram")
VER_OUT = os.path.join(BASE, "docs/manuscript/verification/primary")

H5_PATH = os.path.join(STAGING, "prepared_data.h5")
META_PATH = os.path.join(STAGING, "modeling_metadata.csv")

MISSING = -1
N_DRAWS = int(os.environ.get("N_DRAWS", "1000"))
PANEL_SIZE = 500

print("=" * 70)
print(f"H2/H3 Independent seed test: SEED={SEED}, N_DRAWS={N_DRAWS}")
print("=" * 70)
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")

# Load metadata
meta = pd.read_csv(META_PATH)
sample_ids = meta["sample_id"].values
datasets = meta["dataset"].values
n_samples = len(sample_ids)
print(f"  Metadata: {n_samples} samples")

# Load expression
with h5py.File(H5_PATH, "r") as h5:
    expr_samples = [s.decode() if isinstance(s, bytes) else s for s in h5["sample_ids"][:]]
    gene_names = [g.decode() if isinstance(g, bytes) else g for g in h5["gene_names"][:]]
    expr_mat = h5["zscore_expression"][:].T

expr_order = {s: i for i, s in enumerate(expr_samples)}
expr_idx = [expr_order[s] for s in sample_ids if s in expr_order]
expr_mat = expr_mat[expr_idx]
print(f"  Expression: {expr_mat.shape}")

# Task config
NAS_DATASETS = ["GSE130970", "GSE135251", "GSE162694", "GSE174478", "GSE193066"]
FIB_DATASETS = NAS_DATASETS + ["GSE240729"]
TASK_CFG = {
    "fibrosis": (5, True, "fib_stage", "loco_fold_fibrosis", FIB_DATASETS),
    "nas_composite": (9, True, "nas_score", "loco_fold_nas", NAS_DATASETS),
    "severity": (4, True, "severity4", "loco_fold_fibrosis", FIB_DATASETS),
}

labels = {}
for task, (n_cls, is_ord, col, fold_col, ds_list) in TASK_CFG.items():
    y = meta[col].fillna(MISSING).astype(int).values
    labels[task] = y


def run_elastic_net_loco(X, task_name):
    n_cls, is_ord, col, fold_col, ds_list = TASK_CFG[task_name]
    y = labels[task_name]
    results = []

    for fold_ds in ds_list:
        test_mask = (datasets == fold_ds) & (y != MISSING)
        train_mask = (datasets != fold_ds) & (y != MISSING)
        ds_mask = np.isin(datasets, ds_list)
        test_mask = test_mask & ds_mask
        train_mask = train_mask & ds_mask
        if test_mask.sum() == 0 or train_mask.sum() == 0:
            continue

        X_train, X_test = X[train_mask], X[test_mask]
        y_train, y_test = y[train_mask], y[test_mask]

        # Within-fold feature selection: top 500 by variance
        if X_train.shape[1] > 500:
            var = X_train.var(axis=0)
            top_idx = np.argsort(var)[-500:]
            X_train = X_train[:, top_idx]
            X_test = X_test[:, top_idx]

        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        # Ordinal (cumulative logit)
        preds = np.zeros(len(y_test))
        for threshold in range(1, n_cls):
            y_bin_train = (y_train >= threshold).astype(int)
            if len(np.unique(y_bin_train)) < 2:
                continue
            lr = LogisticRegression(
                penalty="elasticnet", solver="saga", l1_ratio=0.5,
                max_iter=2000, C=1.0, random_state=MODEL_SEED,
            )
            lr.fit(X_train, y_bin_train)
            preds += lr.predict_proba(X_test)[:, 1]
        preds = np.round(preds).astype(int)
        preds = np.clip(preds, 0, n_cls - 1)

        qwk = cohen_kappa_score(y_test, preds, weights="quadratic")
        results.append({"task": task_name, "fold": fold_ds, "qwk": qwk})
    return results


def run_random_gene_baseline(n_draws, seed):
    rng = np.random.RandomState(seed)
    out = []
    for draw in range(n_draws):
        gene_idx = rng.choice(expr_mat.shape[1], size=PANEL_SIZE, replace=False)
        X_random = expr_mat[:, gene_idx]
        for task_name in TASK_CFG:
            res = run_elastic_net_loco(X_random, task_name)
            if res:
                mean_qwk = np.mean([r["qwk"] for r in res])
            else:
                mean_qwk = np.nan
            out.append({"task": task_name, "draw": draw, "seed": seed, "mean_qwk": mean_qwk})
        if (draw + 1) % 50 == 0:
            print(f"  [{time.strftime('%H:%M:%S')}] Draw {draw+1}/{n_draws}")
    return pd.DataFrame(out)


print(f"\nRunning {N_DRAWS} random draws with seed={SEED} (model seed fixed at {MODEL_SEED})...")
df = run_random_gene_baseline(n_draws=N_DRAWS, seed=SEED)

os.makedirs(VER_OUT, exist_ok=True)
out_path = os.path.join(VER_OUT, f"multiprogram_random_baselines_seed{SEED}_n{N_DRAWS}.csv")
df.to_csv(out_path, index=False)
print(f"\nSaved raw draws: {out_path}")

# Compute summary p-values
ablation = pd.read_csv(os.path.join(OUTDIR, "multiprogram_ablation.csv"))
rows = []
for task_name in TASK_CFG:
    curated = ablation[
        (ablation["task"] == task_name) & (ablation["modality"] == "expression")
    ]["mean_qwk"].values[0]
    rand_q = df[df["task"] == task_name]["mean_qwk"].values
    rand_q = rand_q[~np.isnan(rand_q)]
    p_val = (np.sum(rand_q >= curated) + 1) / (len(rand_q) + 1)
    rows.append({
        "task": task_name,
        "seed": SEED,
        "n_draws": len(rand_q),
        "curated_qwk": curated,
        "random_mean": rand_q.mean(),
        "random_sd": rand_q.std(ddof=1),
        "random_min": rand_q.min(),
        "random_max": rand_q.max(),
        "p_value_random_ge_curated": p_val,
        "effect_size_sd": (curated - rand_q.mean())/rand_q.std(ddof=1),
        "n_random_beating_curated": int(np.sum(rand_q >= curated)),
    })
summary = pd.DataFrame(rows)
sum_path = os.path.join(VER_OUT, f"h2_h3_independent_seed_summary_seed{SEED}_n{N_DRAWS}.csv")
summary.to_csv(sum_path, index=False)
print(f"\nSaved summary: {sum_path}")
print("\n" + summary.to_string())

print(f"\nCompleted: {time.strftime('%Y-%m-%d %H:%M:%S')}")
