#!/usr/bin/env python3
"""
Script 177: Integrated Patient Fingerprint
============================================
Creates a compact multi-modal fingerprint per patient and predicts multiple
endpoints simultaneously, then clusters patients by their fingerprint.

Approach:
  1. Build multi-modal feature matrix from:
       - TF activity (PCA top 50)
       - Cell-type proportions (13)
       - Divergence genes (top 100 from feature_matrix_v2)
       Total: ~163 features
  2. Train 4 separate elastic nets per LOCO fold for:
       - fib_ge3 (binary fibrosis >= F3)
       - loco_s2 (binary LOCO-NMF subtype)
       - stellate_high (median split on stellate proportion)
       - severe (fib_ge3 OR stellate_high -- composite endpoint)
  3. The "fingerprint" = predicted probability vector per patient
  4. Cluster patients by fingerprint (k=2-5)
  5. Compare multi-modal model vs Script 161 single-task baselines

Key question: Do patients cluster into meaningful subtypes beyond S1/S2
when viewed through a multi-task lens?

Inputs:
  - results/staging_classifier/tf_activity_features.csv
  - results/staging_classifier/modeling_metadata.csv
  - results/progression/cibersortx_celltype_expression/bayesprism_proportions.csv
  - results/prognosis_v2/labels_v2.csv
  - results/prognosis_v2/feature_matrix_v2.csv   (div_* genes)

Outputs to results/novel_ml/patient_fingerprint/:
  - fingerprint_multitask_results.csv  -- AUROC per target per fold
  - fingerprint_vectors.csv            -- N x K predicted probability fingerprints
  - fingerprint_clusters.csv           -- K-means cluster assignments (k=2-5)
  - fingerprint_vs_singletask.csv      -- comparison with Script 161 results
  - fingerprint_summary.csv            -- summary statistics
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
from pathlib import Path
from sklearn.linear_model import LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.metrics import (
    roc_auc_score,
    average_precision_score,
    brier_score_loss,
    silhouette_score,
    adjusted_rand_score,
)

warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
))
INTEGRATION = BASE / "RNA-seq/Human/Patient_Cohorts/analysis/integration"
STAGING = INTEGRATION / "results/staging_classifier"
PROGRESSION = INTEGRATION / "results/progression/cibersortx_celltype_expression"
PROGNOSIS = INTEGRATION / "results/prognosis_v2"
OUTDIR = INTEGRATION / "results/novel_ml/patient_fingerprint"
OUTDIR.mkdir(parents=True, exist_ok=True)

# Input files
TF_FILE = STAGING / "tf_activity_features.csv"
META_FILE = STAGING / "modeling_metadata.csv"
PROP_FILE = PROGRESSION / "bayesprism_proportions.csv"
LABEL_FILE = PROGNOSIS / "labels_v2.csv"
FEATURE_FILE = PROGNOSIS / "feature_matrix_v2.csv"

# Script 161 results for comparison
BASELINE_FILE = PROGNOSIS / "nested_loco_summary.csv"

# Parameters
RANDOM_STATE = 42
MAX_ITER = 5000
N_TF_PCA = 50
L1_RATIOS = [0.1, 0.5, 0.9]
N_INNER_FOLDS = 5
K_RANGE = range(2, 6)  # k=2,3,4,5 for clustering


def log(msg):
    print(f"[177] {msg}", flush=True)


def compute_metrics(y_true, y_prob):
    """Compute AUROC, AUPRC, Brier."""
    metrics = {}
    try:
        metrics["auroc"] = roc_auc_score(y_true, y_prob)
    except ValueError:
        metrics["auroc"] = np.nan
    try:
        metrics["auprc"] = average_precision_score(y_true, y_prob)
    except ValueError:
        metrics["auprc"] = np.nan
    try:
        metrics["brier"] = brier_score_loss(y_true, y_prob)
    except ValueError:
        metrics["brier"] = np.nan
    return metrics


def fit_elastic_net_cv(X_train, y_train, X_test, n_inner=5):
    """Fit elastic net with inner CV, return test predictions."""
    scaler = StandardScaler().fit(X_train)
    X_tr = scaler.transform(X_train)
    X_te = scaler.transform(X_test)

    # Ensure inner folds do not exceed minority class count
    min_class = min(y_train.sum(), (1 - y_train).sum())
    inner_cv = min(n_inner, int(min_class))
    if inner_cv < 2:
        inner_cv = 2

    model = LogisticRegressionCV(
        penalty="elasticnet",
        solver="saga",
        l1_ratios=L1_RATIOS,
        Cs=10,
        cv=inner_cv,
        scoring="roc_auc",
        max_iter=MAX_ITER,
        random_state=RANDOM_STATE,
        n_jobs=-1,
        class_weight="balanced",
    )
    model.fit(X_tr, y_train)
    y_prob = model.predict_proba(X_te)[:, 1]
    return y_prob, model, scaler


def main():
    t0 = time.time()
    log("Loading data...")

    # ---- 1. Load all data sources ----
    # TF activity (last column is sample_id)
    tf_raw = pd.read_csv(TF_FILE)
    if "sample_id" in tf_raw.columns:
        tf_raw = tf_raw.set_index("sample_id")
    else:
        # Fallback: align with modeling_metadata row order
        meta = pd.read_csv(META_FILE)
        tf_raw.index = meta["sample_id"].values

    sample_ids = tf_raw.index.values
    n_total = len(sample_ids)
    log(f"  Total samples: {n_total}")

    tf_df = tf_raw.copy()
    tf_cols = list(tf_df.columns)

    # Cell-type proportions
    prop_df = pd.read_csv(PROP_FILE)
    prop_df = prop_df.set_index("sample_id").reindex(sample_ids).fillna(0.0)
    ct_cols = list(prop_df.columns)

    # Divergence genes from feature_matrix_v2
    feat_df = pd.read_csv(FEATURE_FILE)
    feat_df = feat_df.set_index("sample_id").reindex(sample_ids).fillna(0.0)
    div_cols = [c for c in feat_df.columns if c.startswith("div_")]
    div_df = feat_df[div_cols]

    # Labels + LOCO folds
    labels = pd.read_csv(LABEL_FILE)
    labels = labels.set_index("sample_id").reindex(sample_ids)

    log(f"  TF activity: {tf_df.shape[1]} TFs")
    log(f"  Cell-type proportions: {len(ct_cols)} cell types")
    log(f"  Divergence genes: {len(div_cols)} genes")

    # ---- 2. Derive additional targets ----
    # stellate_high: median split on stellate proportion
    stellate_vals = prop_df["Stellate"].values
    stellate_median = np.median(stellate_vals[stellate_vals > 0])
    # Use median of non-zero values; if all non-zero, use overall median
    if stellate_median == 0:
        stellate_median = np.median(stellate_vals)
    stellate_high = (stellate_vals > stellate_median).astype(float)
    labels["stellate_high"] = stellate_high

    # severe: fib_ge3 OR stellate_high (composite)
    fib_vals = labels["fib_ge3"].values.copy()
    severe = np.where(
        np.isnan(fib_vals),
        np.nan,
        np.maximum(fib_vals, stellate_high).astype(float),
    )
    labels["severe"] = severe

    # ---- 3. Build combined feature matrix ----
    # TF PCA
    log("Building feature matrix...")
    scaler_tf_all = StandardScaler().fit(tf_df.values)
    tf_scaled = scaler_tf_all.transform(tf_df.values)
    n_components = min(N_TF_PCA, tf_scaled.shape[0], tf_scaled.shape[1])
    pca_tf_all = PCA(n_components=n_components, random_state=RANDOM_STATE)
    tf_pca = pca_tf_all.fit_transform(tf_scaled)
    tf_pca_cols = [f"tf_pc{i}" for i in range(n_components)]
    tf_pca_df = pd.DataFrame(tf_pca, index=sample_ids, columns=tf_pca_cols)

    log(f"  TF PCA: {n_components} components "
        f"(explained variance: {pca_tf_all.explained_variance_ratio_.sum():.3f})")

    # Combined: cell-type proportions + TF-PCA + divergence genes
    combined_df = pd.concat([prop_df, tf_pca_df, div_df], axis=1)
    combined_cols = list(combined_df.columns)
    log(f"  Combined features: {len(combined_cols)} total")

    # ---- 4. Filter to valid LOCO samples ----
    valid_mask = (labels["loco_fold_fibrosis"] != "excluded") & labels["fib_ge3"].notna()
    valid_ids = sample_ids[valid_mask.values]
    log(f"  Valid LOCO samples: {len(valid_ids)}")

    labels_valid = labels.loc[valid_ids]
    combined_valid = combined_df.loc[valid_ids]
    folds = labels_valid["loco_fold_fibrosis"].values
    unique_folds = sorted(set(folds))

    # ---- 5. Multi-task LOCO-CV ----
    target_list = ["fib_ge3", "loco_s2", "stellate_high", "severe"]
    all_results = []

    # Storage for fingerprint vectors (out-of-fold predictions)
    fingerprint_probs = {t: np.full(len(valid_ids), np.nan) for t in target_list}

    log(f"\n=== Multi-task LOCO-CV ({len(unique_folds)} folds, {len(target_list)} targets) ===")

    for fold_name in unique_folds:
        test_mask = folds == fold_name
        train_mask = ~test_mask

        test_idx = np.where(test_mask)[0]
        train_idx = np.where(train_mask)[0]

        X_train_all = combined_valid.values[train_idx]
        X_test_all = combined_valid.values[test_idx]

        log(f"\nFold {fold_name}: train={train_mask.sum()}, test={test_mask.sum()}")

        for target_name in target_list:
            y_col = labels_valid[target_name].values

            # Filter to samples with valid target
            train_valid_t = train_mask & ~np.isnan(y_col)
            test_valid_t = test_mask & ~np.isnan(y_col)

            if test_valid_t.sum() < 5 or train_valid_t.sum() < 20:
                log(f"  {target_name}: too few samples, skipping")
                continue

            y_train = y_col[train_valid_t].astype(int)
            y_test = y_col[test_valid_t].astype(int)

            if len(np.unique(y_test)) < 2:
                log(f"  {target_name}: single class in test, skipping")
                continue

            if min(y_train.sum(), (1 - y_train).sum()) < 2:
                log(f"  {target_name}: <2 minority in train, skipping")
                continue

            X_tr = combined_valid.values[train_valid_t]
            X_te = combined_valid.values[test_valid_t]

            try:
                y_prob, model, scaler = fit_elastic_net_cv(X_tr, y_train, X_te)
                metrics = compute_metrics(y_test, y_prob)
                log(f"  {target_name:16s}: AUROC={metrics['auroc']:.3f}  "
                    f"AUPRC={metrics['auprc']:.3f}")

                # Store out-of-fold predictions for fingerprint
                test_global_idx = np.where(test_valid_t)[0]
                for i, idx in enumerate(test_global_idx):
                    fingerprint_probs[target_name][idx] = y_prob[i]

                all_results.append({
                    "target": target_name,
                    "fold": fold_name,
                    "auroc": metrics["auroc"],
                    "auprc": metrics["auprc"],
                    "brier": metrics["brier"],
                    "n_train": len(y_train),
                    "n_test": len(y_test),
                    "train_pos_rate": float(y_train.mean()),
                    "test_pos_rate": float(y_test.mean()),
                })
            except Exception as e:
                log(f"  {target_name}: FAILED - {e}")

    # ---- 6. Save multi-task results ----
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(OUTDIR / "fingerprint_multitask_results.csv", index=False)
    log(f"\nSaved multi-task results: {len(results_df)} rows")

    # Summary
    summary_rows = []
    for target_name in target_list:
        subset = results_df[results_df["target"] == target_name]
        if len(subset) == 0:
            continue
        summary_rows.append({
            "target": target_name,
            "n_folds": len(subset),
            "mean_auroc": subset["auroc"].mean(),
            "std_auroc": subset["auroc"].std(),
            "mean_auprc": subset["auprc"].mean(),
            "std_auprc": subset["auprc"].std(),
            "mean_brier": subset["brier"].mean(),
        })
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(OUTDIR / "fingerprint_summary.csv", index=False)

    log("\n=== Multi-Task Summary ===")
    for _, row in summary_df.iterrows():
        log(f"  {row['target']:16s}: AUROC={row['mean_auroc']:.3f} +/- {row['std_auroc']:.3f} | "
            f"AUPRC={row['mean_auprc']:.3f}")

    # ---- 7. Build fingerprint vectors ----
    fp_data = {}
    for t in target_list:
        fp_data[f"p_{t}"] = fingerprint_probs[t]

    fp_df = pd.DataFrame(fp_data, index=valid_ids)
    fp_df.index.name = "sample_id"

    # Drop samples where all fingerprints are NaN
    fp_complete = fp_df.dropna()
    log(f"\nFingerprint vectors: {len(fp_complete)} complete (of {len(fp_df)} total)")
    fp_df.to_csv(OUTDIR / "fingerprint_vectors.csv")

    # ---- 8. Cluster patients by fingerprint ----
    if len(fp_complete) < 20:
        log("WARNING: Too few complete fingerprints for clustering")
    else:
        log("\n=== Fingerprint Clustering ===")
        fp_scaled = StandardScaler().fit_transform(fp_complete.values)

        cluster_results = []
        best_k = 2
        best_silhouette = -1

        for k in K_RANGE:
            if k >= len(fp_complete):
                continue
            km = KMeans(n_clusters=k, random_state=RANDOM_STATE, n_init=10)
            cluster_labels = km.fit_predict(fp_scaled)

            sil = silhouette_score(fp_scaled, cluster_labels)
            cluster_results.append({"k": k, "silhouette": sil, "inertia": km.inertia_})
            log(f"  k={k}: silhouette={sil:.3f}, inertia={km.inertia_:.1f}")

            if sil > best_silhouette:
                best_silhouette = sil
                best_k = k

        log(f"  Best k={best_k} (silhouette={best_silhouette:.3f})")

        # Save cluster assignments for all k values
        cluster_df = fp_complete.copy()
        for k in K_RANGE:
            if k >= len(fp_complete):
                continue
            km = KMeans(n_clusters=k, random_state=RANDOM_STATE, n_init=10)
            cluster_df[f"cluster_k{k}"] = km.fit_predict(fp_scaled)

        # Add metadata for interpretation
        for col in ["fib_ge3", "loco_s2", "loco_fold_fibrosis", "fib_ordinal"]:
            if col in labels_valid.columns:
                cluster_df[col] = labels_valid.loc[fp_complete.index, col].values

        cluster_df.to_csv(OUTDIR / "fingerprint_clusters.csv")

        # Characterize best-k clusters
        best_col = f"cluster_k{best_k}"
        log(f"\n=== Cluster Characterization (k={best_k}) ===")
        for c in range(best_k):
            mask = cluster_df[best_col] == c
            n = mask.sum()
            mean_probs = cluster_df.loc[mask, [f"p_{t}" for t in target_list]].mean()
            fib_rate = cluster_df.loc[mask, "fib_ge3"].mean() if "fib_ge3" in cluster_df else np.nan
            s2_rate = cluster_df.loc[mask, "loco_s2"].mean() if "loco_s2" in cluster_df else np.nan
            log(f"  Cluster {c} (n={n}):")
            log(f"    Fibrosis rate: {fib_rate:.3f}, S2 rate: {s2_rate:.3f}")
            for t in target_list:
                log(f"    P({t}): {mean_probs[f'p_{t}']:.3f}")

        # ARI vs S1/S2 subtypes
        if "loco_s2" in cluster_df.columns:
            valid_s2 = cluster_df["loco_s2"].dropna()
            if len(valid_s2) > 10:
                s2_labels = valid_s2.values.astype(int)
                km_labels = cluster_df.loc[valid_s2.index, f"cluster_k2"].values
                ari = adjusted_rand_score(s2_labels, km_labels)
                log(f"\n  ARI(fingerprint k=2 vs S1/S2): {ari:.3f}")
                if ari > 0.5:
                    log("    -> Fingerprint clusters largely recapitulate S1/S2")
                elif ari > 0.1:
                    log("    -> Partial overlap with S1/S2 -- fingerprint captures additional structure")
                else:
                    log("    -> Fingerprint clusters are ORTHOGONAL to S1/S2")

    # ---- 9. Compare with Script 161 baselines ----
    log("\n=== Comparison with Script 161 Baselines ===")
    comparison_rows = []

    if BASELINE_FILE.exists():
        baseline = pd.read_csv(BASELINE_FILE)

        for target_name in ["fib_ge3", "loco_s2"]:
            # Our result
            ours = summary_df[summary_df["target"] == target_name]
            if len(ours) == 0:
                continue
            our_auroc = ours["mean_auroc"].values[0]

            # Script 161 results
            for _, brow in baseline[baseline["target"] == target_name].iterrows():
                comparison_rows.append({
                    "target": target_name,
                    "source": f"161_{brow['config']}",
                    "mean_auroc": brow["mean_auroc"],
                    "std_auroc": brow["std_auroc"],
                })

            comparison_rows.append({
                "target": target_name,
                "source": "177_multimodal",
                "mean_auroc": float(our_auroc),
                "std_auroc": float(ours["std_auroc"].values[0]),
            })

        if comparison_rows:
            comp_df = pd.DataFrame(comparison_rows)
            comp_df.to_csv(OUTDIR / "fingerprint_vs_singletask.csv", index=False)

            for target_name in ["fib_ge3", "loco_s2"]:
                log(f"\n  {target_name}:")
                t_rows = comp_df[comp_df["target"] == target_name].sort_values(
                    "mean_auroc", ascending=False
                )
                for _, row in t_rows.iterrows():
                    log(f"    {row['source']:25s}: AUROC={row['mean_auroc']:.3f} "
                        f"+/- {row['std_auroc']:.3f}")
    else:
        log("  Script 161 baseline file not found -- skipping comparison")

    elapsed = time.time() - t0
    log(f"\nDone in {elapsed:.0f}s ({elapsed / 60:.1f} min)")


if __name__ == "__main__":
    main()
