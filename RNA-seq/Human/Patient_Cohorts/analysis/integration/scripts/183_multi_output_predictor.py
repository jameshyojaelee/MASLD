#!/usr/bin/env python3
"""
183_multi_output_predictor.py
Unified Multi-Output Prediction: Molecular Fate + Therapeutic Pathway + Active Transition

The flagship prediction model for the Cell Metabolism paper. Predicts 3 clinically
meaningful outcomes simultaneously using shared multi-task feature selection:
  Output 1: S2 progressor fate (binary) — "where is this patient headed?"
  Output 2: Therapeutic pathway (3-class) — "which drug class?"
    CAVEAT: Therapeutic pathway labels are derived from Script 149's
    classify_patient_modality(), which uses the SAME feature space (cell-type
    proportions, TF activity, CCC scores) as this predictor. The F1 score for
    output2 reflects definitional consistency, NOT genuine predictive power.
    Do NOT report this F1 as a standalone predictive result.
  Output 3: Active transition (3-class) — "which disease program is running?"

SLURM: cpu, 16 CPUs, 32G, 48h
Env: micromamba activate spatial
"""
import os, sys, logging, warnings
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV, LogisticRegression
from sklearn.metrics import roc_auc_score, average_precision_score, f1_score, accuracy_score
from sklearn.model_selection import StratifiedKFold
from scipy.stats import spearmanr

warnings.filterwarnings("ignore")
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INTEG = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
OUTDIR = os.path.join(INTEG, "results/prognosis_v2/multi_output")
os.makedirs(OUTDIR, exist_ok=True)

SEED = 42
np.random.seed(SEED)
N_PERM = 50
K_GRID = [25, 50, 75, 100]

# =========================================================================
# PART 1: Feature Assembly
# =========================================================================
log.info("=== PART 1: Feature Assembly ===")

# 1a. Cell-type proportions (17)
bp = pd.read_csv(os.path.join(INTEG, "results/progression/cibersortx_celltype_expression/bayesprism_proportions.csv"))
bp = bp.set_index("sample_id")
bp.columns = ["ct_" + c for c in bp.columns]
bp["ct_stellate_gt4pct"] = (bp.get("ct_Stellate", bp.iloc[:, 0]) > 0.04).astype(float)
stel_col = [c for c in bp.columns if "stellate" in c.lower() and "gt4" not in c]
hep_col = [c for c in bp.columns if "hepatocyte" in c.lower() or "Hepatocyte" in c]
if stel_col and hep_col:
    bp["ct_stellate_hep_ratio"] = bp[stel_col[0]] / (bp[hep_col[0]] + 1e-6)
immune_cols = [c for c in bp.columns if any(x in c for x in ["Macrophage", "Monocyte", "NK", "T_cell"])]
if immune_cols:
    bp["ct_immune_infiltration"] = bp[immune_cols].sum(axis=1)
chol_col = [c for c in bp.columns if "Cholangiocyte" in c]
if stel_col and chol_col:
    bp["ct_fibrogenic_axis"] = bp[stel_col[0]] + bp[chol_col[0]]
log.info(f"  Cell-type features: {bp.shape[1]} columns, {bp.shape[0]} samples")

# 1b. TF activity (top 50)
tf = pd.read_csv(os.path.join(INTEG, "results/staging_classifier/tf_activity_features.csv"))
if "sample_id" in tf.columns:
    tf = tf.set_index("sample_id")
else:
    meta_tmp = pd.read_csv(os.path.join(INTEG, "results/staging_classifier/modeling_metadata.csv"))
    tf.index = meta_tmp["sample_id"].values
tf_var = tf.var().sort_values(ascending=False)
top50_tf = tf_var.head(50).index.tolist()
tf = tf[top50_tf]
tf.columns = ["tf_" + c for c in tf.columns]
log.info(f"  TF activity: {tf.shape[1]} columns (top 50 by variance)")

# 1c. Divergence genes (top 50)
feat_v2 = pd.read_csv(os.path.join(INTEG, "results/prognosis_v2/feature_matrix_v2.csv"))
feat_v2 = feat_v2.set_index("sample_id")
div_cols = [c for c in feat_v2.columns if c.startswith("div_")]
div_var = feat_v2[div_cols].var().sort_values(ascending=False)
top50_div = div_var.head(50).index.tolist()
div_feats = feat_v2[top50_div]
log.info(f"  Divergence genes: {div_feats.shape[1]} columns (top 50 by variance)")

# 1d. CCC axis scores (8)
log.info("  Loading CCC scores (chunked)...")
ccc_agg = {}
for chunk in pd.read_csv(os.path.join(INTEG, "results/progression/ccc_lr_scores.csv"), chunksize=500000):
    grp = chunk.groupby(["sample_id", "axis"])["score"].mean()
    for (sid, axis), val in grp.items():
        if sid not in ccc_agg:
            ccc_agg[sid] = {}
        if axis not in ccc_agg[sid]:
            ccc_agg[sid][axis] = []
        ccc_agg[sid][axis].append(val)
# Average across chunks
ccc_rows = []
for sid, axes in ccc_agg.items():
    row = {"sample_id": sid}
    for axis, vals in axes.items():
        row[f"ccc_{axis}"] = np.mean(vals)
    ccc_rows.append(row)
ccc = pd.DataFrame(ccc_rows).set_index("sample_id")
log.info(f"  CCC axes: {ccc.shape[1]} columns, {ccc.shape[0]} samples")

# 1e. Clinical (4)
meta = pd.read_csv(os.path.join(INTEG, "results/staging_classifier/modeling_metadata.csv"))
meta = meta.set_index("sample_id")
clin = pd.DataFrame(index=meta.index)
clin["clin_sex"] = meta["sex"].map({"F": 1, "M": 0}).astype(float)
clin["clin_age"] = pd.to_numeric(meta["age"], errors="coerce")
clin["clin_fib_stage"] = pd.to_numeric(meta["fibrosis_stage"], errors="coerce").replace(-1, np.nan)
clin["clin_nas_score"] = pd.to_numeric(meta["nas_score"], errors="coerce").replace(-1, np.nan)
log.info(f"  Clinical: {clin.shape[1]} columns")

# Merge all features
common = sorted(set(bp.index) & set(tf.index) & set(div_feats.index) & set(ccc.index) & set(clin.index))
log.info(f"  Common samples: {len(common)}")
X = pd.concat([bp.loc[common], tf.loc[common], div_feats.loc[common], ccc.loc[common], clin.loc[common]], axis=1)
X = X.apply(pd.to_numeric, errors="coerce")
log.info(f"  Final feature matrix: {X.shape}")

# =========================================================================
# PART 2: Label Assembly
# =========================================================================
log.info("\n=== PART 2: Label Assembly ===")

labels_v2 = pd.read_csv(os.path.join(INTEG, "results/prognosis_v2/labels_v2.csv")).set_index("sample_id")
pathway = pd.read_csv(os.path.join(INTEG, "results/progression/patient_pathway_class.csv")).set_index("sample_id")
assign = pd.read_csv(os.path.join(INTEG, "results/progression/transition_subtype_assignments.csv")).set_index("sample_id")

# Output 1: S2 fate
y1 = labels_v2["loco_s2"].reindex(common)
log.info(f"  Output 1 (S2 fate): {y1.notna().sum()} labeled ({(y1==1).sum()} S2, {(y1==0).sum()} S1)")

# Output 2: Therapeutic pathway
pmap = {"anti_fibrotic": 0, "vascular": 1}
y2_raw = pathway["dominant_modality"].reindex(common)
y2 = y2_raw.map(lambda x: pmap.get(x, 2) if pd.notna(x) else np.nan)
for cls, name in [(0, "anti_fibrotic"), (1, "vascular"), (2, "other")]:
    log.info(f"    Output 2 class {cls} ({name}): {(y2==cls).sum()}")

# Output 3: Active transition (3-class)
tmap = {"quiescent": 0, "F0_to_F1": 1, "F1_to_F2": 1, "F2_to_F3": 2, "F3_to_F4": 2}
y3_raw = assign["fib_dominant"].reindex(common)
y3 = y3_raw.map(lambda x: tmap.get(x, np.nan) if pd.notna(x) else np.nan)
for cls, name in [(0, "quiescent"), (1, "early_active"), (2, "late_active")]:
    log.info(f"    Output 3 class {cls} ({name}): {(y3==cls).sum()}")

# LOCO folds
folds = labels_v2["loco_fold_fibrosis"].reindex(common)
valid_folds = sorted(folds[folds.notna() & (folds != "excluded")].unique())
log.info(f"  LOCO folds: {valid_folds}")

# =========================================================================
# PART 3: Utility Functions
# =========================================================================

def median_impute(X_train, X_test):
    """Impute NaN with training column medians."""
    medians = X_train.median()
    return X_train.fillna(medians), X_test.fillna(medians)

def univariate_auroc(X, y, binary=True):
    """Compute univariate AUROC per feature."""
    scores = {}
    valid = y.notna()
    Xv, yv = X[valid], y[valid]
    for col in Xv.columns:
        vals = Xv[col].values
        nan_mask = ~np.isnan(vals)
        if nan_mask.sum() < 20 or len(np.unique(yv[nan_mask])) < 2:
            scores[col] = 0.5
            continue
        try:
            if binary:
                scores[col] = roc_auc_score(yv[nan_mask], vals[nan_mask])
            else:
                # For multi-class, use mean of OvR AUROCs
                from sklearn.preprocessing import label_binarize
                yb = label_binarize(yv[nan_mask], classes=sorted(yv[nan_mask].unique()))
                if yb.shape[1] == 1:
                    scores[col] = roc_auc_score(yb, vals[nan_mask])
                else:
                    aucs = []
                    for k in range(yb.shape[1]):
                        try:
                            aucs.append(roc_auc_score(yb[:, k], vals[nan_mask]))
                        except:
                            aucs.append(0.5)
                    scores[col] = np.mean(aucs)
        except:
            scores[col] = 0.5
    return pd.Series(scores).sort_values(ascending=False)

def fit_predict(X_tr, y_tr, X_te, binary=True):
    """Fit elastic net and predict."""
    if binary:
        model = LogisticRegressionCV(
            penalty="elasticnet", l1_ratios=[0.5], solver="saga",
            max_iter=5000, class_weight="balanced", cv=3, scoring="roc_auc",
            random_state=SEED, n_jobs=-1
        )
        model.fit(X_tr, y_tr)
        prob = model.predict_proba(X_te)[:, 1]
    else:
        model = LogisticRegressionCV(
            penalty="elasticnet", l1_ratios=[0.5], solver="saga",
            max_iter=5000, class_weight="balanced", cv=3,
            scoring="f1_macro", random_state=SEED, n_jobs=-1
        )
        model.fit(X_tr, y_tr)
        prob = model.predict_proba(X_te)
    return model, prob

def compute_metrics(y_true, y_prob, binary=True):
    """Compute metrics for one output."""
    metrics = {}
    if binary:
        metrics["auroc"] = roc_auc_score(y_true, y_prob)
        metrics["auprc"] = average_precision_score(y_true, y_prob)
    else:
        pred = np.argmax(y_prob, axis=1) if y_prob.ndim > 1 else y_prob
        metrics["macro_f1"] = f1_score(y_true, pred, average="macro", zero_division=0)
        metrics["accuracy"] = accuracy_score(y_true, pred)
        # Adjacent accuracy (within ±1 class for ordinal)
        metrics["adj_accuracy"] = np.mean(np.abs(y_true - pred) <= 1)
    return metrics

# =========================================================================
# PART 4: Multi-Task LOCO-CV
# =========================================================================
log.info("\n=== PART 4: Multi-Task LOCO-CV ===")

all_results = []
all_predictions = []
all_importance = []

outputs = [
    ("output1_s2_fate", y1, True),
    ("output2_pathway", y2, False),
    ("output3_transition", y3, False),
]

for fold in valid_folds:
    log.info(f"\n--- Fold: {fold} ---")
    train_mask = (folds != fold) & folds.notna() & (folds != "excluded")
    test_mask = folds == fold

    X_train_raw = X[train_mask]
    X_test_raw = X[test_mask]

    # Impute + scale
    X_tr_imp, X_te_imp = median_impute(X_train_raw, X_test_raw)
    scaler = StandardScaler()
    X_tr_scaled = pd.DataFrame(scaler.fit_transform(X_tr_imp), index=X_tr_imp.index, columns=X_tr_imp.columns)
    X_te_scaled = pd.DataFrame(scaler.transform(X_te_imp), index=X_te_imp.index, columns=X_te_imp.columns)

    # Multi-task feature ranking: max AUROC across tasks
    log.info("  Multi-task feature ranking...")
    task_scores = []
    for oname, y_out, is_binary in outputs:
        y_tr = y_out.reindex(X_tr_scaled.index)
        valid_tr = y_tr.notna()
        if valid_tr.sum() > 20:
            sc = univariate_auroc(X_tr_scaled[valid_tr], y_tr[valid_tr], binary=is_binary)
            task_scores.append(sc)

    if task_scores:
        combined = pd.concat(task_scores, axis=1).max(axis=1).sort_values(ascending=False)
    else:
        combined = pd.Series(0.5, index=X_tr_scaled.columns)

    # Try K values, pick best via quick inner CV
    best_k = 50
    best_score = 0
    for k in K_GRID:
        top_feats = combined.head(k).index.tolist()
        # Quick evaluation on first output (S2 fate)
        y_tr_o1 = y1.reindex(X_tr_scaled.index)
        valid_o1 = y_tr_o1.notna()
        if valid_o1.sum() > 30:
            try:
                inner_cv = StratifiedKFold(3, shuffle=True, random_state=SEED)
                inner_scores = []
                for itr, ite in inner_cv.split(X_tr_scaled[valid_o1][top_feats], y_tr_o1[valid_o1]):
                    Xitr = X_tr_scaled[valid_o1][top_feats].iloc[itr]
                    Xite = X_tr_scaled[valid_o1][top_feats].iloc[ite]
                    yitr = y_tr_o1[valid_o1].iloc[itr]
                    yite = y_tr_o1[valid_o1].iloc[ite]
                    m = LogisticRegression(penalty="l2", C=1, solver="saga", max_iter=3000, class_weight="balanced", random_state=SEED)
                    m.fit(Xitr, yitr)
                    inner_scores.append(roc_auc_score(yite, m.predict_proba(Xite)[:, 1]))
                mean_sc = np.mean(inner_scores)
                if mean_sc > best_score:
                    best_score = mean_sc
                    best_k = k
            except:
                pass

    selected_features = combined.head(best_k).index.tolist()
    log.info(f"  Selected {best_k} features (inner AUROC={best_score:.3f})")

    # Fit and predict for each output
    for oname, y_out, is_binary in outputs:
        y_tr = y_out.reindex(X_tr_scaled.index)
        y_te = y_out.reindex(X_te_scaled.index)

        valid_tr = y_tr.notna()
        valid_te = y_te.notna()

        n_tr = valid_tr.sum()
        n_te = valid_te.sum()

        if n_tr < 10 or n_te < 5:
            log.warning(f"    {oname}: skipped (n_train={n_tr}, n_test={n_te})")
            continue

        # Check class balance
        classes_tr = y_tr[valid_tr].unique()
        if len(classes_tr) < 2:
            log.warning(f"    {oname}: skipped (single class in training)")
            continue

        try:
            model, prob = fit_predict(
                X_tr_scaled[valid_tr][selected_features].values,
                y_tr[valid_tr].values.astype(int),
                X_te_scaled[valid_te][selected_features].values,
                binary=is_binary
            )

            metrics = compute_metrics(y_te[valid_te].values.astype(int), prob, binary=is_binary)
            metrics["fold"] = fold
            metrics["output"] = oname
            metrics["config"] = "multi_task"
            metrics["n_train"] = n_tr
            metrics["n_test"] = n_te
            metrics["n_features"] = best_k
            all_results.append(metrics)

            # Store predictions
            if is_binary:
                for sid, p in zip(X_te_scaled[valid_te].index, prob):
                    all_predictions.append({"sample_id": sid, "output": oname, "config": "multi_task", "fold": fold, "prob_class1": p, "true_label": int(y_te[sid])})
            else:
                for sid, p in zip(X_te_scaled[valid_te].index, prob):
                    pred_dict = {"sample_id": sid, "output": oname, "config": "multi_task", "fold": fold, "true_label": int(y_te[sid])}
                    for ci in range(p.shape[0] if p.ndim > 0 else 1):
                        pred_dict[f"prob_class{ci}"] = p[ci] if p.ndim > 0 else p
                    all_predictions.append(pred_dict)

            # Feature importance
            if hasattr(model, "coef_"):
                coefs = model.coef_.flatten() if is_binary else model.coef_.mean(axis=0)
                for fi, fname in enumerate(selected_features):
                    if fi < len(coefs):
                        all_importance.append({"feature": fname, "output": oname, "fold": fold, "coefficient": coefs[fi]})

            primary = metrics.get("auroc", metrics.get("macro_f1", 0))
            log.info(f"    {oname}: {'AUROC' if is_binary else 'F1'}={primary:.3f} (n_train={n_tr}, n_test={n_te})")

        except Exception as e:
            log.error(f"    {oname}: FAILED — {e}")

# =========================================================================
# PART 5: Single-Task Baselines
# =========================================================================
log.info("\n=== PART 5: Single-Task Baselines ===")

for fold in valid_folds:
    train_mask = (folds != fold) & folds.notna() & (folds != "excluded")
    test_mask = folds == fold

    X_tr_imp, X_te_imp = median_impute(X[train_mask], X[test_mask])
    scaler = StandardScaler()
    X_tr_s = pd.DataFrame(scaler.fit_transform(X_tr_imp), index=X_tr_imp.index, columns=X_tr_imp.columns)
    X_te_s = pd.DataFrame(scaler.transform(X_te_imp), index=X_te_imp.index, columns=X_te_imp.columns)

    for oname, y_out, is_binary in outputs:
        y_tr = y_out.reindex(X_tr_s.index)
        y_te = y_out.reindex(X_te_s.index)
        valid_tr = y_tr.notna()
        valid_te = y_te.notna()

        if valid_tr.sum() < 10 or valid_te.sum() < 5:
            continue
        if len(y_tr[valid_tr].unique()) < 2:
            continue

        # Task-specific feature selection
        sc = univariate_auroc(X_tr_s[valid_tr], y_tr[valid_tr], binary=is_binary)
        st_feats = sc.head(50).index.tolist()

        try:
            model, prob = fit_predict(
                X_tr_s[valid_tr][st_feats].values,
                y_tr[valid_tr].values.astype(int),
                X_te_s[valid_te][st_feats].values,
                binary=is_binary
            )
            metrics = compute_metrics(y_te[valid_te].values.astype(int), prob, binary=is_binary)
            metrics["fold"] = fold
            metrics["output"] = oname
            metrics["config"] = "single_task"
            metrics["n_train"] = valid_tr.sum()
            metrics["n_test"] = valid_te.sum()
            metrics["n_features"] = 50
            all_results.append(metrics)
        except:
            pass

    # Clinical-only baseline
    clin_cols = [c for c in X.columns if c.startswith("clin_")]
    for oname, y_out, is_binary in outputs:
        y_tr = y_out.reindex(X_tr_s.index)
        y_te = y_out.reindex(X_te_s.index)
        valid_tr = y_tr.notna()
        valid_te = y_te.notna()

        if valid_tr.sum() < 10 or valid_te.sum() < 5:
            continue
        if len(y_tr[valid_tr].unique()) < 2:
            continue

        try:
            model, prob = fit_predict(
                X_tr_s[valid_tr][clin_cols].values,
                y_tr[valid_tr].values.astype(int),
                X_te_s[valid_te][clin_cols].values,
                binary=is_binary
            )
            metrics = compute_metrics(y_te[valid_te].values.astype(int), prob, binary=is_binary)
            metrics["fold"] = fold
            metrics["output"] = oname
            metrics["config"] = "clinical_only"
            metrics["n_train"] = valid_tr.sum()
            metrics["n_test"] = valid_te.sum()
            metrics["n_features"] = len(clin_cols)
            all_results.append(metrics)
        except:
            pass

# =========================================================================
# PART 6: Permutation Null
# =========================================================================
log.info(f"\n=== PART 6: Permutation Null ({N_PERM} permutations) ===")

# NULL FIX (mega-review A7): the permutation null must RE-SELECT features
# INSIDE each permutation. Previously it used X_tr_s[:, :50] (the first 50
# columns of the full matrix) and never paid the feature-selection cost, so
# the null underestimated chance performance and inflated the apparent
# significance of the real AUROC. We now rank features with the SAME
# univariate_auroc screen the real model uses, but computed on the SHUFFLED
# training labels, and take the top PERM_K features per fold/permutation.
PERM_K = 50
for oname, y_out, is_binary in outputs:
    perm_metrics = []
    for pi in range(N_PERM):
        fold_metrics = []
        for fold in valid_folds:
            train_mask = (folds != fold) & folds.notna() & (folds != "excluded")
            test_mask = folds == fold

            X_tr_imp, X_te_imp = median_impute(X[train_mask], X[test_mask])
            scaler = StandardScaler()
            X_tr_s = pd.DataFrame(scaler.fit_transform(X_tr_imp), index=X_tr_imp.index, columns=X_tr_imp.columns)
            X_te_s = pd.DataFrame(scaler.transform(X_te_imp), index=X_te_imp.index, columns=X_te_imp.columns)

            y_tr = y_out.reindex(X_tr_s.index)
            y_te = y_out.reindex(X_te_s.index)
            valid_tr = y_tr.notna()
            valid_te = y_te.notna()

            if valid_tr.sum() < 10 or valid_te.sum() < 5:
                continue

            # Shuffle training labels (positions aligned to valid training rows)
            y_tr_perm = pd.Series(y_tr[valid_tr].values.copy(), index=X_tr_s[valid_tr].index)
            shuffled = y_tr_perm.values.copy()
            np.random.shuffle(shuffled)
            y_tr_perm[:] = shuffled

            if len(np.unique(y_tr_perm.values)) < 2:
                continue

            # Re-select features INSIDE the permutation, on the shuffled labels.
            try:
                perm_rank = univariate_auroc(X_tr_s[valid_tr], y_tr_perm, binary=is_binary)
                perm_feats = perm_rank.head(PERM_K).index.tolist()
            except Exception:
                perm_feats = list(X_tr_s.columns[:PERM_K])

            try:
                _, prob = fit_predict(
                    X_tr_s[valid_tr][perm_feats].values,
                    y_tr_perm.values.astype(int),
                    X_te_s[valid_te][perm_feats].values,
                    binary=is_binary,
                )
                m = compute_metrics(y_te[valid_te].values.astype(int), prob, binary=is_binary)
                fold_metrics.append(m.get("auroc", m.get("macro_f1", 0)))
            except:
                pass

        if fold_metrics:
            perm_metrics.append(np.mean(fold_metrics))

    if perm_metrics:
        all_results.append({
            "output": oname, "config": "permutation_null",
            "auroc" if is_binary else "macro_f1": np.mean(perm_metrics),
            "fold": "all", "n_train": 0, "n_test": 0, "n_features": PERM_K,
            "perm_std": np.std(perm_metrics), "n_perms": len(perm_metrics)
        })
        log.info(f"  {oname} permutation null: {np.mean(perm_metrics):.3f} ± {np.std(perm_metrics):.3f}")

# =========================================================================
# PART 7: Summary and Save
# =========================================================================
log.info("\n=== PART 7: Summary and Save ===")

results_df = pd.DataFrame(all_results)
results_df.to_csv(os.path.join(OUTDIR, "multi_output_results.csv"), index=False)
log.info(f"  Saved multi_output_results.csv ({len(results_df)} rows)")

pred_df = pd.DataFrame(all_predictions)
pred_df.to_csv(os.path.join(OUTDIR, "multi_output_predictions.csv"), index=False)
log.info(f"  Saved multi_output_predictions.csv ({len(pred_df)} rows)")

imp_df = pd.DataFrame(all_importance)
imp_df.to_csv(os.path.join(OUTDIR, "multi_output_feature_importance.csv"), index=False)
log.info(f"  Saved multi_output_feature_importance.csv ({len(imp_df)} rows)")

# Summary
summary_rows = []
for oname in ["output1_s2_fate", "output2_pathway", "output3_transition"]:
    for config in ["multi_task", "single_task", "clinical_only"]:
        sub = results_df[(results_df["output"] == oname) & (results_df["config"] == config)]
        if len(sub) == 0:
            continue
        row = {"output": oname, "config": config, "n_folds": len(sub)}
        for metric in ["auroc", "auprc", "macro_f1", "accuracy", "adj_accuracy"]:
            if metric in sub.columns and sub[metric].notna().any():
                row[f"mean_{metric}"] = sub[metric].mean()
                row[f"std_{metric}"] = sub[metric].std()
        summary_rows.append(row)

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(os.path.join(OUTDIR, "multi_output_summary.csv"), index=False)
log.info(f"  Saved multi_output_summary.csv")

# Print summary
log.info("\n=== SUMMARY ===")
for _, row in summary_df.iterrows():
    primary = row.get("mean_auroc", row.get("mean_macro_f1", "N/A"))
    std = row.get("std_auroc", row.get("std_macro_f1", "N/A"))
    log.info(f"  {row['output']:25s} | {row['config']:15s} | metric={primary:.3f} ± {std:.3f} | n_folds={row['n_folds']}")

log.info("\n=== 183: COMPLETE ===")
