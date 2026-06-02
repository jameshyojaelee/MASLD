#!/usr/bin/env python3
"""
74_label_imputation.py
Impute staging labels for unannotated cohorts with ensemble prediction
and conformal prediction for calibrated uncertainty.

Phase 5 (Shared Outputs) -- label imputation for samples without fibrosis/NAS.

Inputs:
  - results/staging_classifier/prepared_data.h5
  - results/staging_classifier/ordinal_loco_fibrosis.csv   (Plan 1)
  - results/staging_classifier/tier3b_fibrosis_results.csv  (Plan 2)
  - results/staging_classifier/embeddings_all_samples.csv   (Plan 3)
  - results/staging_classifier/modeling_metadata.csv

Outputs:
  - imputed_labels_all_plans.csv
  - conformal_prediction_sets.csv
  - imputation_confidence.csv
  - imputation_validation.csv
  - imputation_summary.csv

Usage: python 74_label_imputation.py
SLURM: cpu, 8 CPUs, 64GB RAM, 48h
Env:   micromamba activate rapids_singlecell
"""

import os
import sys
import time
import warnings
import numpy as np
import pandas as pd
import h5py
from collections import defaultdict

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score, balanced_accuracy_score, cohen_kappa_score,
    roc_auc_score, f1_score
)
from sklearn.neighbors import NearestNeighbors

RANDOM_STATE = 42
np.random.seed(RANDOM_STATE)

# =============================================================================
# PATHS
# =============================================================================
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR = os.path.join(INT, "results")
OUTDIR = os.path.join(RDIR, "staging_classifier")
os.makedirs(OUTDIR, exist_ok=True)

H5_PATH = os.path.join(OUTDIR, "prepared_data.h5")

print("=" * 60)
print("74: Label Imputation with Conformal Prediction")
print("=" * 60)
print(f"Output: {OUTDIR}")
print(f"Started: {time.strftime('%Y-%m-%d %H:%M:%S')}")
print()

# =============================================================================
# LOAD DATA
# =============================================================================
print("Loading data...")

if not os.path.exists(H5_PATH):
    sys.exit(f"ERROR: {H5_PATH} not found. Run Script 62 first.")

with h5py.File(H5_PATH, "r") as h5:
    rank_expr = h5["rank_expression"][:].T.astype(np.float32)  # samples x genes
    gene_names = np.array([g.decode() if isinstance(g, bytes) else g
                           for g in h5["gene_names"][:]])
    sample_ids = np.array([s.decode() if isinstance(s, bytes) else s
                           for s in h5["sample_ids"][:]])

    meta = {}
    for key in h5["metadata"].keys():
        vals = h5["metadata"][key][:]
        if vals.dtype.kind == "S" or vals.dtype.kind == "O":
            meta[key] = np.array([v.decode() if isinstance(v, bytes) else v
                                  for v in vals])
        else:
            meta[key] = vals.astype(np.float32)

N_SAMPLES, N_GENES = rank_expr.shape
print(f"  Expression: {N_SAMPLES} samples x {N_GENES} genes")

meta_df = pd.DataFrame(meta)
meta_df["sample_id"] = sample_ids

# Load modeling metadata for additional columns
meta_path = os.path.join(OUTDIR, "modeling_metadata.csv")
if os.path.exists(meta_path):
    meta_full = pd.read_csv(meta_path)
    print(f"  Metadata (full): {meta_full.shape}")
else:
    meta_full = meta_df.copy()
    print("  WARNING: modeling_metadata.csv not found; using HDF5 metadata")

# Extract labels
fib_stage = meta["fib_stage"].astype(np.int64)
nas_group4 = meta["nas_group4"].astype(np.int64)
datasets = meta["dataset"] if "dataset" in meta else np.array(["unknown"] * N_SAMPLES)

labeled_fib = fib_stage >= 0
unlabeled_fib = fib_stage < 0
labeled_nas = nas_group4 >= 0
unlabeled_nas = nas_group4 < 0

print(f"  Fibrosis labeled: {labeled_fib.sum()}, unlabeled: {unlabeled_fib.sum()}")
print(f"  NAS labeled: {labeled_nas.sum()}, unlabeled: {unlabeled_nas.sum()}")
print()

# =============================================================================
# LOAD PLAN-SPECIFIC OUTPUTS (graceful degradation if missing)
# =============================================================================
print("Loading plan-specific outputs...")

plan_outputs = {}

# Plan 1: ordinal_loco_fibrosis.csv
p1_path = os.path.join(OUTDIR, "ordinal_loco_fibrosis.csv")
if os.path.exists(p1_path):
    plan_outputs["plan1"] = pd.read_csv(p1_path)
    print(f"  Plan 1 (ordinal): {plan_outputs['plan1'].shape}")
else:
    print("  Plan 1 ordinal results not found; will use fresh model")

# Plan 2: tier3b_fibrosis_results.csv
p2_path = os.path.join(OUTDIR, "tier3b_fibrosis_results.csv")
if os.path.exists(p2_path):
    plan_outputs["plan2"] = pd.read_csv(p2_path)
    print(f"  Plan 2 (tier3b): {plan_outputs['plan2'].shape}")
else:
    print("  Plan 2 tier3b results not found; will use fresh model")

# Plan 3: embeddings_all_samples.csv
p3_path = os.path.join(OUTDIR, "embeddings_all_samples.csv")
has_embeddings = False
if os.path.exists(p3_path):
    emb_df = pd.read_csv(p3_path, index_col="sample_id")
    emb_cols = [c for c in emb_df.columns if c.startswith("z")]
    if len(emb_cols) > 0:
        has_embeddings = True
        print(f"  Plan 3 embeddings: {emb_df.shape[0]} samples x {len(emb_cols)} dims")
    else:
        print("  Plan 3 embeddings found but no z columns")
else:
    print("  Plan 3 embeddings not found; will use fresh model")

n_plans = len(plan_outputs) + (1 if has_embeddings else 0)
print(f"\n  Available plans for ensemble: {n_plans}")
if n_plans == 0:
    print("  WARNING: No plan outputs found. Will train fresh models only.\n")
print()


# =============================================================================
# HELPER: Train a classifier and predict probabilities
# =============================================================================
def train_and_predict(X_train, y_train, X_test, model_type="rf"):
    """Train a classifier and return probability matrix for test set."""
    classes = np.sort(np.unique(y_train))
    n_classes = len(classes)

    scaler = StandardScaler()
    X_tr = scaler.fit_transform(X_train)
    X_te = scaler.transform(X_test)

    if model_type == "rf":
        model = RandomForestClassifier(
            n_estimators=500, max_features="sqrt",
            class_weight="balanced", random_state=RANDOM_STATE, n_jobs=-1
        )
    elif model_type == "lr":
        model = LogisticRegression(
            penalty="elasticnet", l1_ratio=0.5, solver="saga",
            max_iter=5000, class_weight="balanced",
            random_state=RANDOM_STATE
        )
    else:
        model = RandomForestClassifier(
            n_estimators=500, class_weight="balanced",
            random_state=RANDOM_STATE, n_jobs=-1
        )

    model.fit(X_tr, y_train)
    probs = model.predict_proba(X_te)  # (n_test, n_classes)
    return probs, model.classes_, model


def prob_to_full_matrix(probs, model_classes, all_classes):
    """Expand probability matrix to cover all classes (fill missing with 0)."""
    n_test = probs.shape[0]
    full = np.zeros((n_test, len(all_classes)))
    for i, cls in enumerate(model_classes):
        idx = np.where(all_classes == cls)[0]
        if len(idx) > 0:
            full[:, idx[0]] = probs[:, i]
    # Renormalize rows
    row_sums = full.sum(axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1
    full = full / row_sums
    return full


# =============================================================================
# SECTION 1: ENSEMBLE PREDICTION FOR FIBROSIS
# =============================================================================
print("=" * 60)
print("SECTION 1: Ensemble Fibrosis Prediction")
print("=" * 60)

all_fib_classes = np.array([0, 1, 2, 3, 4])

# Train models on ALL labeled data, predict on unlabeled
X_labeled = rank_expr[labeled_fib]
y_labeled = fib_stage[labeled_fib]
X_unlabeled = rank_expr[unlabeled_fib]
unlabeled_ids = sample_ids[unlabeled_fib]

print(f"Training set: {X_labeled.shape[0]} labeled samples")
print(f"Prediction set: {X_unlabeled.shape[0]} unlabeled samples")

ensemble_probs_list = []
plan_names = []

# Fresh RF model (always available)
print("\n  Training fresh Random Forest...")
rf_probs, rf_classes, _ = train_and_predict(X_labeled, y_labeled, X_unlabeled, "rf")
rf_full = prob_to_full_matrix(rf_probs, rf_classes, all_fib_classes)
ensemble_probs_list.append(rf_full)
plan_names.append("rf_fresh")

# Fresh Logistic Regression
print("  Training fresh Logistic Regression...")
try:
    lr_probs, lr_classes, _ = train_and_predict(X_labeled, y_labeled, X_unlabeled, "lr")
    lr_full = prob_to_full_matrix(lr_probs, lr_classes, all_fib_classes)
    ensemble_probs_list.append(lr_full)
    plan_names.append("lr_fresh")
except Exception as e:
    print(f"    WARNING: LR failed: {e}")

# Embedding-based model (if available)
if has_embeddings:
    print("  Training RF on VAE embeddings...")
    try:
        # Align embeddings with labeled/unlabeled
        emb_sample_ids = emb_df.index.values
        labeled_ids_set = set(sample_ids[labeled_fib])
        unlabeled_ids_set = set(sample_ids[unlabeled_fib])

        emb_labeled_mask = np.array([s in labeled_ids_set for s in emb_sample_ids])
        emb_unlabeled_mask = np.array([s in unlabeled_ids_set for s in emb_sample_ids])

        if emb_labeled_mask.sum() > 10 and emb_unlabeled_mask.sum() > 0:
            X_emb_train = emb_df.loc[emb_labeled_mask, emb_cols].values
            # Need to align y with embedding order
            emb_labeled_ids = emb_sample_ids[emb_labeled_mask]
            y_emb_train = np.array([
                fib_stage[np.where(sample_ids == sid)[0][0]]
                for sid in emb_labeled_ids
            ])
            X_emb_test = emb_df.loc[emb_unlabeled_mask, emb_cols].values

            emb_probs, emb_classes, _ = train_and_predict(
                X_emb_train, y_emb_train, X_emb_test, "rf"
            )
            emb_full = prob_to_full_matrix(emb_probs, emb_classes, all_fib_classes)

            # Reorder to match unlabeled_ids
            emb_unlabeled_ids = emb_sample_ids[emb_unlabeled_mask]
            id_map = {sid: i for i, sid in enumerate(emb_unlabeled_ids)}
            reorder_idx = [id_map.get(sid, -1) for sid in unlabeled_ids]
            valid_reorder = [i for i in reorder_idx if i >= 0]

            if len(valid_reorder) == len(unlabeled_ids):
                emb_full = emb_full[reorder_idx]
                ensemble_probs_list.append(emb_full)
                plan_names.append("emb_rf")
            else:
                print(f"    WARNING: embedding ID mismatch ({len(valid_reorder)}"
                      f"/{len(unlabeled_ids)}); skipping")
    except Exception as e:
        print(f"    WARNING: embedding model failed: {e}")

print(f"\n  Ensemble components: {plan_names}")

# Compute ensemble average
if len(ensemble_probs_list) > 0 and X_unlabeled.shape[0] > 0:
    ensemble_probs = np.mean(ensemble_probs_list, axis=0)
    ensemble_preds = all_fib_classes[np.argmax(ensemble_probs, axis=1)]
    max_probs = np.max(ensemble_probs, axis=1)

    # Per-plan predictions
    per_plan_preds = {}
    for pname, probs in zip(plan_names, ensemble_probs_list):
        per_plan_preds[pname] = all_fib_classes[np.argmax(probs, axis=1)]

    # Build output table
    imputed_df = pd.DataFrame({
        "sample_id": unlabeled_ids,
        "ensemble_pred": ensemble_preds,
        "ensemble_max_prob": max_probs,
    })

    for pname in plan_names:
        imputed_df[f"pred_{pname}"] = per_plan_preds[pname]

    # Add probability columns
    for i, cls in enumerate(all_fib_classes):
        imputed_df[f"prob_F{cls}"] = ensemble_probs[:, i]

    # Compute agreement across plans
    if len(plan_names) > 1:
        pred_matrix = np.column_stack([per_plan_preds[p] for p in plan_names])
        n_agree = np.array([
            np.max(np.bincount(row.astype(int), minlength=5))
            for row in pred_matrix
        ])
        imputed_df["n_plans_agree"] = n_agree
    else:
        imputed_df["n_plans_agree"] = 1

    # Add dataset info
    unlabeled_datasets = datasets[unlabeled_fib]
    imputed_df["dataset"] = unlabeled_datasets

    print(f"\n  Imputed {len(imputed_df)} samples")
    print(f"  Ensemble prediction distribution:")
    for cls in all_fib_classes:
        n = (ensemble_preds == cls).sum()
        print(f"    F{cls}: {n} ({100*n/len(ensemble_preds):.1f}%)")
else:
    imputed_df = pd.DataFrame()
    ensemble_probs = np.array([])
    print("  No unlabeled samples or no models; skipping imputation")

# Save
imputed_path = os.path.join(OUTDIR, "imputed_labels_all_plans.csv")
if len(imputed_df) > 0:
    imputed_df.to_csv(imputed_path, index=False)
    print(f"  Saved imputed_labels_all_plans.csv")
print()


# =============================================================================
# SECTION 2: CONFORMAL PREDICTION
# =============================================================================
print("=" * 60)
print("SECTION 2: Conformal Prediction (Split-Conformal)")
print("=" * 60)

N_CONFORMAL_REPS = 10
ALPHA = 0.10  # 90% coverage guarantee
conformal_results = []

if labeled_fib.sum() >= 50 and unlabeled_fib.sum() > 0:
    for rep in range(N_CONFORMAL_REPS):
        rng = np.random.RandomState(RANDOM_STATE + rep)

        # Stratified split: 80% train, 20% calibration
        cal_indices = []
        train_indices = []

        for cls in np.unique(y_labeled):
            cls_idx = np.where(y_labeled == cls)[0]
            n_cal = max(1, int(0.2 * len(cls_idx)))
            perm = rng.permutation(len(cls_idx))
            cal_indices.extend(cls_idx[perm[:n_cal]])
            train_indices.extend(cls_idx[perm[n_cal:]])

        cal_indices = np.array(cal_indices)
        train_indices = np.array(train_indices)

        X_tr = X_labeled[train_indices]
        y_tr = y_labeled[train_indices]
        X_cal = X_labeled[cal_indices]
        y_cal = y_labeled[cal_indices]

        # Train model on training portion
        cal_probs, cal_classes, model = train_and_predict(X_tr, y_tr, X_cal, "rf")
        cal_full = prob_to_full_matrix(cal_probs, cal_classes, all_fib_classes)

        # Compute nonconformity scores on calibration set: 1 - P(true class)
        nonconf_scores = np.zeros(len(y_cal))
        for i in range(len(y_cal)):
            cls_idx = np.where(all_fib_classes == y_cal[i])[0]
            if len(cls_idx) > 0:
                nonconf_scores[i] = 1.0 - cal_full[i, cls_idx[0]]
            else:
                nonconf_scores[i] = 1.0

        # Quantile threshold (for 1-alpha coverage)
        n_cal = len(nonconf_scores)
        q_level = np.ceil((1 - ALPHA) * (n_cal + 1)) / n_cal
        q_level = min(q_level, 1.0)
        threshold = np.quantile(nonconf_scores, q_level)

        # Predict on unlabeled samples
        scaler = StandardScaler()
        X_tr_s = scaler.fit_transform(X_tr)
        X_unlabeled_s = scaler.transform(X_unlabeled)

        rf = RandomForestClassifier(
            n_estimators=500, max_features="sqrt",
            class_weight="balanced", random_state=RANDOM_STATE, n_jobs=-1
        )
        rf.fit(X_tr_s, y_tr)
        test_probs = rf.predict_proba(X_unlabeled_s)
        test_full = prob_to_full_matrix(test_probs, rf.classes_, all_fib_classes)

        # For each test sample, prediction set = {k : P(k) >= 1 - threshold}
        for i in range(len(unlabeled_ids)):
            pred_set = all_fib_classes[test_full[i, :] >= (1.0 - threshold)]
            if len(pred_set) == 0:
                # Fallback: include the most probable class
                pred_set = all_fib_classes[[np.argmax(test_full[i, :])]]

            conformal_results.append({
                "sample_id": unlabeled_ids[i],
                "rep": rep,
                "prediction_set": ";".join(map(str, pred_set)),
                "set_size": len(pred_set),
                "max_prob": np.max(test_full[i, :]),
                "threshold": threshold,
            })

    conformal_df = pd.DataFrame(conformal_results)
    print(f"  Conformal prediction: {len(conformal_df)} entries "
          f"({N_CONFORMAL_REPS} reps x {unlabeled_fib.sum()} samples)")

    # Aggregate across repetitions
    conformal_agg = conformal_df.groupby("sample_id").agg(
        mean_set_size=("set_size", "mean"),
        min_set_size=("set_size", "min"),
        max_set_size=("set_size", "max"),
        mean_max_prob=("max_prob", "mean"),
    ).reset_index()

    # Most common prediction set per sample
    def most_common_set(group):
        sets = group["prediction_set"].values
        from collections import Counter
        return Counter(sets).most_common(1)[0][0]

    common_sets = conformal_df.groupby("sample_id").apply(most_common_set).reset_index()
    common_sets.columns = ["sample_id", "consensus_prediction_set"]
    conformal_agg = conformal_agg.merge(common_sets, on="sample_id")

    # Confidence level
    def assign_confidence(row):
        if row["mean_max_prob"] > 0.7 and row["mean_set_size"] <= 1.5:
            return "High"
        elif row["mean_max_prob"] >= 0.5 or row["mean_set_size"] <= 2.5:
            return "Medium"
        else:
            return "Low"

    conformal_agg["confidence_level"] = conformal_agg.apply(assign_confidence, axis=1)

    conformal_agg.to_csv(os.path.join(OUTDIR, "conformal_prediction_sets.csv"), index=False)
    print(f"  Saved conformal_prediction_sets.csv ({len(conformal_agg)} samples)")

    # Confidence summary
    conf_counts = conformal_agg["confidence_level"].value_counts()
    for level in ["High", "Medium", "Low"]:
        n = conf_counts.get(level, 0)
        print(f"    {level}: {n} samples ({100*n/len(conformal_agg):.1f}%)")

else:
    conformal_agg = pd.DataFrame()
    print("  Insufficient labeled data for conformal prediction")
print()


# =============================================================================
# SECTION 3: IMPUTATION CONFIDENCE TABLE
# =============================================================================
print("=" * 60)
print("SECTION 3: Confidence Table")
print("=" * 60)

if len(imputed_df) > 0 and len(conformal_agg) > 0:
    confidence_df = imputed_df[["sample_id", "ensemble_pred", "ensemble_max_prob",
                                 "n_plans_agree", "dataset"]].copy()

    # Merge conformal info
    confidence_df = confidence_df.merge(
        conformal_agg[["sample_id", "mean_set_size", "confidence_level",
                       "consensus_prediction_set"]],
        on="sample_id", how="left"
    )

    confidence_df.to_csv(os.path.join(OUTDIR, "imputation_confidence.csv"), index=False)
    print(f"  Saved imputation_confidence.csv ({len(confidence_df)} samples)")
elif len(imputed_df) > 0:
    # No conformal results; use ensemble confidence only
    confidence_df = imputed_df[["sample_id", "ensemble_pred", "ensemble_max_prob",
                                 "n_plans_agree", "dataset"]].copy()

    def assign_conf_no_conformal(row):
        if row["ensemble_max_prob"] > 0.7:
            return "High"
        elif row["ensemble_max_prob"] >= 0.5:
            return "Medium"
        else:
            return "Low"

    confidence_df["confidence_level"] = confidence_df.apply(
        assign_conf_no_conformal, axis=1)
    confidence_df.to_csv(os.path.join(OUTDIR, "imputation_confidence.csv"), index=False)
    print(f"  Saved imputation_confidence.csv ({len(confidence_df)} samples, no conformal)")
else:
    confidence_df = pd.DataFrame()
    print("  No imputed samples; skipping confidence table")
print()


# =============================================================================
# SECTION 4: EMBEDDING-BASED IMPUTATION (kNN in VAE space)
# =============================================================================
print("=" * 60)
print("SECTION 4: Embedding-Based (kNN) Imputation")
print("=" * 60)

if has_embeddings and unlabeled_fib.sum() > 0:
    try:
        emb_sample_ids = emb_df.index.values
        emb_values = emb_df[emb_cols].values

        # Align labeled and unlabeled in embedding space
        labeled_ids_list = sample_ids[labeled_fib]
        unlabeled_ids_list = sample_ids[unlabeled_fib]

        emb_labeled_idx = [i for i, s in enumerate(emb_sample_ids) if s in set(labeled_ids_list)]
        emb_unlabeled_idx = [i for i, s in enumerate(emb_sample_ids) if s in set(unlabeled_ids_list)]

        if len(emb_labeled_idx) > 10 and len(emb_unlabeled_idx) > 0:
            X_emb_lab = emb_values[emb_labeled_idx]
            X_emb_unl = emb_values[emb_unlabeled_idx]
            y_emb_lab = np.array([
                fib_stage[np.where(sample_ids == emb_sample_ids[i])[0][0]]
                for i in emb_labeled_idx
            ])
            unl_ids_emb = emb_sample_ids[emb_unlabeled_idx]

            # k-NN in embedding space
            K = 10
            nn = NearestNeighbors(n_neighbors=K, metric="euclidean")
            nn.fit(X_emb_lab)
            distances, indices = nn.kneighbors(X_emb_unl)

            knn_results = []
            for i in range(len(unl_ids_emb)):
                neighbor_stages = y_emb_lab[indices[i]]
                neighbor_dists = distances[i]

                # Inverse distance weighting (add small epsilon to avoid div-by-zero)
                weights = 1.0 / (neighbor_dists + 1e-8)
                weights = weights / weights.sum()

                # Weighted vote
                stage_probs = np.zeros(5)
                for j in range(K):
                    stage_probs[int(neighbor_stages[j])] += weights[j]

                knn_pred = int(np.argmax(stage_probs))
                knn_results.append({
                    "sample_id": unl_ids_emb[i],
                    "knn_pred": knn_pred,
                    "knn_max_prob": stage_probs[knn_pred],
                    "knn_mean_dist": float(np.mean(neighbor_dists)),
                    **{f"knn_prob_F{s}": stage_probs[s] for s in range(5)}
                })

            knn_df = pd.DataFrame(knn_results)

            # Merge with main imputation table
            if len(imputed_df) > 0:
                imputed_df = imputed_df.merge(
                    knn_df[["sample_id", "knn_pred", "knn_max_prob"]],
                    on="sample_id", how="left"
                )
                # Update saved file
                imputed_df.to_csv(imputed_path, index=False)

            print(f"  kNN imputation: {len(knn_df)} samples")
            print(f"  Mean kNN max probability: {knn_df['knn_max_prob'].mean():.3f}")

            # Compare kNN vs model-based
            if len(imputed_df) > 0 and "knn_pred" in imputed_df.columns:
                agree = (imputed_df["ensemble_pred"].values ==
                         imputed_df["knn_pred"].values)
                print(f"  Agreement ensemble vs kNN: {agree.sum()}/{len(agree)} "
                      f"({100*agree.mean():.1f}%)")
        else:
            print("  Insufficient embeddings for kNN imputation")
    except Exception as e:
        print(f"  WARNING: kNN imputation failed: {e}")
else:
    print("  Embeddings not available or no unlabeled samples; skipping kNN")
print()


# =============================================================================
# SECTION 5: VALIDATION AGAINST PARTIAL METADATA
# =============================================================================
print("=" * 60)
print("SECTION 5: Validation Against Partial Metadata")
print("=" * 60)

validation_results = []

if len(imputed_df) == 0:
    print("  No imputed samples; skipping validation")
else:
    # Define partially-annotated datasets and their expected severity mapping
    partial_datasets = {
        "GSE126848": {
            "description": "57 samples, NAFL/NASH labels",
            "severity_mapping": {
                # NAFL -> likely F0-F1, NASH -> likely F2-F4
                "expected": "NASH samples should have higher predicted fibrosis than NAFL"
            }
        },
        "GSE167523": {
            "description": "98 samples, NAFL/NASH labels",
            "severity_mapping": {
                "expected": "NASH samples should have higher predicted fibrosis than NAFL"
            }
        },
    }

    for ds_name, ds_info in partial_datasets.items():
        ds_imputed = imputed_df[imputed_df["dataset"] == ds_name]
        if len(ds_imputed) == 0:
            print(f"\n  {ds_name}: no imputed samples (may be fully annotated)")
            validation_results.append({
                "dataset": ds_name,
                "metric": "n_imputed",
                "value": 0,
                "note": "No imputed samples for this dataset"
            })
            continue

        print(f"\n  {ds_name}: {len(ds_imputed)} imputed samples")

        # Try to get diagnosis labels from metadata
        ds_meta = meta_full[meta_full["sample_id"].isin(ds_imputed["sample_id"])]
        if len(ds_meta) == 0:
            ds_meta = meta_df[meta_df["sample_id"].isin(ds_imputed["sample_id"])]

        validation_results.append({
            "dataset": ds_name,
            "metric": "n_imputed",
            "value": len(ds_imputed),
            "note": ds_info["description"]
        })

        # Check if we have diagnosis or disease_stage information
        for diag_col in ["diagnosis", "diagnosis_harmonized", "condition", "disease_stage"]:
            if diag_col not in ds_meta.columns:
                continue
            diag_vals = ds_meta[diag_col].dropna().unique()
            if len(diag_vals) <= 1:
                continue

            print(f"    Using '{diag_col}' column: {diag_vals}")

            # Merge predictions with diagnosis
            merged = ds_imputed.merge(ds_meta[["sample_id", diag_col]], on="sample_id")
            if len(merged) < 5:
                continue

            # Check if NAFL/NASH-type labels
            nafl_labels = [v for v in diag_vals
                          if any(x in str(v).upper() for x in ["NAFL", "MASL", "STEATOSIS"])]
            nash_labels = [v for v in diag_vals
                          if any(x in str(v).upper() for x in ["NASH", "MASH", "FIBROSIS"])]

            if len(nafl_labels) > 0 and len(nash_labels) > 0:
                nafl_preds = merged[merged[diag_col].isin(nafl_labels)]["ensemble_pred"]
                nash_preds = merged[merged[diag_col].isin(nash_labels)]["ensemble_pred"]

                if len(nafl_preds) > 0 and len(nash_preds) > 0:
                    mean_nafl = nafl_preds.mean()
                    mean_nash = nash_preds.mean()
                    diff = mean_nash - mean_nafl

                    validation_results.append({
                        "dataset": ds_name,
                        "metric": "mean_pred_NAFL",
                        "value": mean_nafl,
                        "note": f"n={len(nafl_preds)}"
                    })
                    validation_results.append({
                        "dataset": ds_name,
                        "metric": "mean_pred_NASH",
                        "value": mean_nash,
                        "note": f"n={len(nash_preds)}"
                    })
                    validation_results.append({
                        "dataset": ds_name,
                        "metric": "NASH_minus_NAFL",
                        "value": diff,
                        "note": f"Expected positive (NASH > NAFL)"
                    })

                    # Mann-Whitney U test
                    from scipy.stats import mannwhitneyu
                    try:
                        stat, pval = mannwhitneyu(
                            nash_preds.values, nafl_preds.values, alternative="greater"
                        )
                        validation_results.append({
                            "dataset": ds_name,
                            "metric": "mannwhitney_pval",
                            "value": pval,
                            "note": f"U={stat:.1f}, NASH > NAFL"
                        })
                        print(f"    NASH vs NAFL: mean {mean_nash:.2f} vs {mean_nafl:.2f}, "
                              f"p={pval:.4f}")
                    except Exception as e:
                        print(f"    Mann-Whitney failed: {e}")

            # Ordinal correlation (if disease_stage is numeric-like)
            if diag_col == "disease_stage":
                try:
                    stage_numeric = pd.to_numeric(merged[diag_col], errors="coerce")
                    valid = ~stage_numeric.isna()
                    if valid.sum() > 5:
                        from scipy.stats import spearmanr
                        rho, pval = spearmanr(
                            stage_numeric[valid],
                            merged.loc[valid, "ensemble_pred"]
                        )
                        validation_results.append({
                            "dataset": ds_name,
                            "metric": "spearman_rho",
                            "value": rho,
                            "note": f"p={pval:.4f}, n={valid.sum()}"
                        })
                        print(f"    Spearman(disease_stage, pred_fib): rho={rho:.3f}, p={pval:.4f}")
                except Exception as e:
                    print(f"    Spearman failed: {e}")

            break  # use first available diagnosis column

    val_df = pd.DataFrame(validation_results)
    val_df.to_csv(os.path.join(OUTDIR, "imputation_validation.csv"), index=False)
    print(f"\n  Saved imputation_validation.csv ({len(val_df)} entries)")
print()


# =============================================================================
# SECTION 6: IMPUTATION SUMMARY
# =============================================================================
print("=" * 60)
print("SECTION 6: Summary")
print("=" * 60)

summary_rows = []
if len(imputed_df) > 0:
    for ds in imputed_df["dataset"].unique():
        ds_sub = imputed_df[imputed_df["dataset"] == ds]
        conf_sub = (confidence_df[confidence_df["dataset"] == ds]
                    if len(confidence_df) > 0 and "dataset" in confidence_df.columns
                    else pd.DataFrame())

        n_high = (conf_sub["confidence_level"] == "High").sum() if len(conf_sub) > 0 else 0
        n_medium = (conf_sub["confidence_level"] == "Medium").sum() if len(conf_sub) > 0 else 0
        n_low = (conf_sub["confidence_level"] == "Low").sum() if len(conf_sub) > 0 else 0

        summary_rows.append({
            "dataset": ds,
            "n_imputed": len(ds_sub),
            "n_high_confidence": n_high,
            "n_medium_confidence": n_medium,
            "n_low_confidence": n_low,
            "mean_ensemble_max_prob": ds_sub["ensemble_max_prob"].mean(),
            "pred_F0": (ds_sub["ensemble_pred"] == 0).sum(),
            "pred_F1": (ds_sub["ensemble_pred"] == 1).sum(),
            "pred_F2": (ds_sub["ensemble_pred"] == 2).sum(),
            "pred_F3": (ds_sub["ensemble_pred"] == 3).sum(),
            "pred_F4": (ds_sub["ensemble_pred"] == 4).sum(),
        })

summary_df = pd.DataFrame(summary_rows)
if len(summary_df) > 0:
    summary_df.to_csv(os.path.join(OUTDIR, "imputation_summary.csv"), index=False)
    print(f"  Saved imputation_summary.csv ({len(summary_df)} datasets)")
    print()
    print(summary_df.to_string(index=False))
else:
    print("  No imputed samples to summarize")

print(f"\nFinished: {time.strftime('%Y-%m-%d %H:%M:%S')}")
