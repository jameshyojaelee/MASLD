"""
226a_plasma_shap_panels.py
Feature importance (SHAP / permutation) + panel-size curves for best plasma models.

For 3 targets (T1 binary fibrosis, T4 etiology 3-class, T5 etiology binary):
  1. Re-run best model in 10×5-fold nested CV (seed=42)
  2. Compute feature importance (SHAP or permutation) per fold
  3. Evaluate panel sizes {2, 5, 10, 20, 50, 100, 200, all} via importance ranking
  4. Evaluate NFASC+GDF15 2-protein baseline
  5. Random protein baselines: 100 draws × {20, 50, 100} per fold
  6. Track per-fold top-10 protein stability

SLURM: sbatch --job-name=226a --partition=cpu --cpus-per-task=16 --mem=64G --time=48:00:00
Env: micromamba activate spatial
"""

import importlib.util
import logging
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.inspection import permutation_importance
from sklearn.linear_model import LogisticRegression, RidgeClassifier
from sklearn.metrics import f1_score, roc_auc_score
from sklearn.model_selection import (
    GridSearchCV,
    RepeatedStratifiedKFold,
    StratifiedKFold,
)
from sklearn.multiclass import OneVsRestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Import shared module (filename starts with digit)
# ---------------------------------------------------------------------------

_common_path = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), "225_plasma_common.py"
)
_spec = importlib.util.spec_from_file_location("plasma_common_225", _common_path)
_plasma_common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_plasma_common)
load_olink_data = _plasma_common.load_olink_data
impute_train_fold = _plasma_common.impute_train_fold

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SEED = 42
np.random.seed(SEED)

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
OUTDIR = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram",
)
os.makedirs(OUTDIR, exist_ok=True)

PANEL_SIZES = [2, 5, 10, 20, 50, 100, 200]  # "all" appended dynamically
N_RANDOM_DRAWS = 100
RANDOM_PANEL_SIZES = [2, 5, 10, 20, 50, 100, 200]  # match all PANEL_SIZES

# ---------------------------------------------------------------------------
# Optional imports
# ---------------------------------------------------------------------------

try:
    import shap

    HAS_SHAP = True
    log.info("shap available")
except ImportError:
    HAS_SHAP = False
    log.info("shap not available — will fall back to permutation importance")

try:
    import xgboost as xgb

    HAS_XGB = True
    log.info("xgboost available")
except ImportError:
    HAS_XGB = False
    log.error("xgboost NOT available — T1 XGBoost will fail")

try:
    from imblearn.ensemble import BalancedBaggingClassifier

    HAS_IMBLEARN = True
    log.info("imbalanced-learn available")
except ImportError:
    HAS_IMBLEARN = False
    log.error("imbalanced-learn NOT available — T5 BalancedBagging will fail")


# ---------------------------------------------------------------------------
# Helper: evaluate a subset of features with LogisticRegression
# ---------------------------------------------------------------------------


def eval_panel(X_train, y_train, X_test, y_test, feature_idx, metric_name):
    """Train LogisticRegression on a feature subset and return the metric.

    Parameters
    ----------
    X_train, X_test : ndarray — already scaled/imputed
    y_train, y_test : ndarray
    feature_idx : array-like of int — column indices to use
    metric_name : str — "auroc" or "macro_f1"

    Returns
    -------
    float — metric value
    """
    if len(feature_idx) == 0:
        return np.nan

    Xtr = X_train[:, feature_idx]
    Xte = X_test[:, feature_idx]

    # Scale the subset
    scaler = StandardScaler()
    Xtr = scaler.fit_transform(Xtr)
    Xte = scaler.transform(Xte)

    n_classes = len(np.unique(y_train))

    if n_classes <= 2:
        lr = LogisticRegression(
            class_weight="balanced", max_iter=2000, random_state=SEED, solver="lbfgs"
        )
        lr.fit(Xtr, y_train)
        if metric_name == "auroc":
            y_prob = lr.predict_proba(Xte)[:, 1]
            return roc_auc_score(y_test, y_prob)
        else:
            y_pred = lr.predict(Xte)
            return f1_score(y_test, y_pred, average="macro")
    else:
        # Multiclass: OVR logistic regression
        lr = LogisticRegression(
            class_weight="balanced",
            max_iter=2000,
            random_state=SEED,
            solver="lbfgs",
        )
        lr.fit(Xtr, y_train)
        if metric_name == "macro_f1":
            y_pred = lr.predict(Xte)
            return f1_score(y_test, y_pred, average="macro")
        else:
            # AUROC OVR
            y_prob = lr.predict_proba(Xte)
            try:
                return roc_auc_score(
                    y_test, y_prob, multi_class="ovr", average="macro"
                )
            except ValueError:
                return np.nan


# ---------------------------------------------------------------------------
# Target configurations
# ---------------------------------------------------------------------------


def get_target_configs(y_targets, n_proteins):
    """Return dict of target name → config dict."""
    # Compute scale_pos_weight for T1
    y_bin = y_targets["binary"]
    n_pos = int(y_bin.sum())
    n_neg = int(len(y_bin) - n_pos)
    spw = n_neg / max(n_pos, 1)

    configs = {}

    # T1: binary fibrosis — XGBoost
    if HAS_XGB:
        configs["T1_binary"] = {
            "y_key": "binary",
            "base_model": xgb.XGBClassifier(
                scale_pos_weight=spw,
                random_state=SEED,
                eval_metric="logloss",
                use_label_encoder=False,
                n_estimators=300,
                max_depth=5,
                learning_rate=0.05,
                subsample=0.8,
                colsample_bytree=0.8,
                reg_lambda=10,
                n_jobs=1,
            ),
            "param_grid": {},  # Best params already fixed
            "importance_method": "shap_tree",
            "metric_name": "auroc",
            "scoring": "roc_auc",
        }
    else:
        log.warning("Skipping T1_binary — xgboost not available")

    # T4: etiology 3-class — Ridge OVR with GridSearch
    configs["T4_etiology_3"] = {
        "y_key": "etiology_3",
        "base_model": OneVsRestClassifier(RidgeClassifier(class_weight="balanced")),
        "param_grid": {
            "model__estimator__alpha": [0.001, 0.01, 0.1, 1, 10, 100, 1000]
        },
        "importance_method": "shap_linear",
        "metric_name": "macro_f1",
        "scoring": "f1_macro",
    }

    # T5: etiology binary — BalancedBagging
    if HAS_IMBLEARN:
        configs["T5_etiology_binary"] = {
            "y_key": "etiology_binary",
            "base_model": BalancedBaggingClassifier(
                estimator=LogisticRegression(max_iter=2000, random_state=SEED),
                n_estimators=50,
                random_state=SEED,
            ),
            "param_grid": {},  # Best params already fixed
            "importance_method": "permutation",
            "metric_name": "auroc",
            "scoring": "roc_auc",
        }
    else:
        log.warning("Skipping T5_etiology_binary — imbalanced-learn not available")

    return configs


# ---------------------------------------------------------------------------
# Importance computation
# ---------------------------------------------------------------------------


def compute_importance_shap_tree(pipeline, X_test_scaled, X_train_scaled):
    """SHAP TreeExplainer for tree-based models. Returns (n_test, n_features)."""
    model = pipeline.named_steps["model"]
    explainer = shap.TreeExplainer(model)
    shap_vals = explainer.shap_values(X_test_scaled)
    # For binary classifiers, shap_values may return a list of 2 arrays
    if isinstance(shap_vals, list):
        shap_vals = shap_vals[1]  # Positive class
    return shap_vals  # (n_test, n_features)


def compute_importance_shap_linear(pipeline, X_test_scaled, X_train_scaled):
    """SHAP LinearExplainer for OVR Ridge. Returns (n_test, n_features, n_classes)."""
    ovr_model = pipeline.named_steps["model"]
    n_classes = len(ovr_model.estimators_)
    n_test, n_feat = X_test_scaled.shape
    all_shap = np.zeros((n_test, n_feat, n_classes), dtype=np.float32)

    for c_idx, estimator in enumerate(ovr_model.estimators_):
        explainer = shap.LinearExplainer(estimator, X_train_scaled)
        sv = explainer.shap_values(X_test_scaled)
        if isinstance(sv, list):
            sv = sv[0]
        all_shap[:, :, c_idx] = sv

    return all_shap  # (n_test, n_features, n_classes)


def compute_importance_permutation(pipeline, X_test, y_test, scoring):
    """Permutation importance. Returns (n_features,) mean importance."""
    result = permutation_importance(
        pipeline,
        X_test,
        y_test,
        n_repeats=10,
        scoring=scoring,
        random_state=SEED,
        n_jobs=4,
    )
    return result.importances_mean  # (n_features,)


def get_importance_ranking(importance_values, method):
    """Convert importance values to a ranking (descending |importance|).

    Returns
    -------
    ranked_indices : ndarray of int — feature indices sorted by importance (most important first)
    importance_flat : ndarray of float — per-feature importance score used for ranking
    """
    if method == "shap_tree":
        # (n_test, n_features) → mean |SHAP| per feature
        importance_flat = np.abs(importance_values).mean(axis=0)
    elif method == "shap_linear":
        # (n_test, n_features, n_classes) → mean |SHAP| across test samples and classes
        importance_flat = np.abs(importance_values).mean(axis=(0, 2))
    elif method == "permutation":
        # Already (n_features,) — may have negative values
        importance_flat = np.abs(importance_values)
    else:
        raise ValueError(f"Unknown method: {method}")

    ranked_indices = np.argsort(importance_flat)[::-1]
    return ranked_indices, importance_flat


# ---------------------------------------------------------------------------
# Main analysis
# ---------------------------------------------------------------------------


def main():
    log.info("=== 226a: Plasma SHAP & Panel Size Analysis ===")
    t_start = time.time()

    # ------------------------------------------------------------------
    # Load data
    # ------------------------------------------------------------------
    X_raw, protein_names_list, meta, y_targets = load_olink_data()
    protein_names = np.array(protein_names_list)
    n_subjects, n_proteins = X_raw.shape
    log.info(f"Data loaded: {n_subjects} subjects × {n_proteins} proteins")

    # Locate NFASC and GDF15
    nfasc_idx = np.where(protein_names == "NFASC")[0]
    gdf15_idx = np.where(protein_names == "GDF15")[0]
    benchmark_idx = np.concatenate([nfasc_idx, gdf15_idx])
    log.info(
        f"Benchmark proteins: NFASC idx={nfasc_idx}, GDF15 idx={gdf15_idx}, "
        f"combined {len(benchmark_idx)} proteins"
    )
    if len(benchmark_idx) != 2:
        log.warning("Expected exactly 2 benchmark proteins (NFASC, GDF15)")

    # ------------------------------------------------------------------
    # Configure targets and models
    # ------------------------------------------------------------------
    target_configs = get_target_configs(y_targets, n_proteins)
    log.info(f"Targets to process: {list(target_configs.keys())}")

    # CV setup
    outer_cv = RepeatedStratifiedKFold(
        n_splits=5, n_repeats=10, random_state=SEED
    )
    inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=SEED)

    # Extended panel sizes including "all"
    panel_sizes_full = PANEL_SIZES + [n_proteins]

    # ------------------------------------------------------------------
    # Accumulators
    # ------------------------------------------------------------------
    all_panel_rows = []
    all_fold_top20_rows = []
    all_importance_rows = []  # per-protein aggregate

    # Per-target SHAP accumulators (for .npz output)
    shap_accumulators = {}

    # ------------------------------------------------------------------
    # Process each target
    # ------------------------------------------------------------------
    for target_name, cfg in target_configs.items():
        log.info(f"\n{'='*60}")
        log.info(f"Target: {target_name}")
        log.info(f"{'='*60}")

        y = y_targets[cfg["y_key"]]
        base_model = cfg["base_model"]
        param_grid = cfg["param_grid"]
        importance_method = cfg["importance_method"]
        metric_name = cfg["metric_name"]
        scoring = cfg["scoring"]

        # Class distribution
        unique, counts = np.unique(y, return_counts=True)
        log.info(
            f"Class distribution: "
            + ", ".join([f"{c}={n}" for c, n in zip(unique, counts)])
        )

        # SHAP accumulation: track per-subject across folds
        # Each subject appears in exactly 10 folds (10 repeats × 1/5 test fraction)
        # Always initialise both accumulators — SHAP may fall back to permutation
        # on individual folds, and we need somewhere to store those values.
        perm_importance_sum = np.zeros(n_proteins, dtype=np.float64)
        n_perm_folds = 0  # track how many folds actually used permutation
        if importance_method in ("shap_tree", "shap_linear"):
            if importance_method == "shap_tree":
                shap_sum = np.zeros((n_subjects, n_proteins), dtype=np.float64)
            else:
                n_classes = len(np.unique(y))
                shap_sum = np.zeros(
                    (n_subjects, n_proteins, n_classes), dtype=np.float64
                )
            shap_count = np.zeros(n_subjects, dtype=int)
        else:
            shap_sum = None
            shap_count = None

        # Fold-level importance for stability tracking
        fold_top10_tracker = np.zeros(
            (50, n_proteins), dtype=bool
        )  # 50 folds × n_proteins

        fold_primary_metrics = []

        t_target_start = time.time()

        for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X_raw, y)):
            t_fold_start = time.time()

            X_train_raw, X_test_raw = X_raw[train_idx], X_raw[test_idx]
            y_train, y_test = y[train_idx], y[test_idx]

            # Per-fold imputation
            X_train, X_test = impute_train_fold(X_train_raw, X_test_raw)

            # ----- 1. Train best model -----
            pipe = Pipeline(
                [("scaler", StandardScaler()), ("model", clone(base_model))]
            )

            if param_grid:
                grid = GridSearchCV(
                    pipe,
                    param_grid,
                    cv=inner_cv,
                    scoring=scoring,
                    n_jobs=4,
                    refit=True,
                )
                grid.fit(X_train, y_train)
                best_pipe = grid.best_estimator_
                if fold_idx < 3:
                    log.info(
                        f"  Fold {fold_idx}: best params = {grid.best_params_}"
                    )
            else:
                pipe.fit(X_train, y_train)
                best_pipe = pipe

            # Compute primary metric on full feature set
            if metric_name == "auroc":
                if hasattr(best_pipe, "predict_proba"):
                    y_prob = best_pipe.predict_proba(X_test)
                    if y_prob.ndim == 2 and y_prob.shape[1] == 2:
                        primary_metric = roc_auc_score(y_test, y_prob[:, 1])
                    else:
                        try:
                            primary_metric = roc_auc_score(
                                y_test, y_prob, multi_class="ovr", average="macro"
                            )
                        except ValueError:
                            primary_metric = np.nan
                elif hasattr(best_pipe, "decision_function"):
                    y_dec = best_pipe.decision_function(X_test)
                    primary_metric = roc_auc_score(y_test, y_dec)
                else:
                    primary_metric = np.nan
            else:
                y_pred = best_pipe.predict(X_test)
                primary_metric = f1_score(y_test, y_pred, average="macro")

            fold_primary_metrics.append(primary_metric)

            # ----- 2. Compute feature importance -----
            X_train_scaled = best_pipe.named_steps["scaler"].transform(X_train)
            X_test_scaled = best_pipe.named_steps["scaler"].transform(X_test)

            importance_values = None
            actual_method = importance_method

            if importance_method == "shap_tree" and HAS_SHAP:
                try:
                    importance_values = compute_importance_shap_tree(
                        best_pipe, X_test_scaled, X_train_scaled
                    )
                except Exception as e:
                    log.warning(
                        f"  Fold {fold_idx}: SHAP TreeExplainer failed ({e}), "
                        f"falling back to permutation"
                    )
                    actual_method = "permutation"

            elif importance_method == "shap_linear" and HAS_SHAP:
                try:
                    importance_values = compute_importance_shap_linear(
                        best_pipe, X_test_scaled, X_train_scaled
                    )
                except Exception as e:
                    log.warning(
                        f"  Fold {fold_idx}: SHAP LinearExplainer failed ({e}), "
                        f"falling back to permutation"
                    )
                    actual_method = "permutation"

            elif importance_method == "permutation":
                pass  # Will compute below
            else:
                # SHAP not available
                log.warning(
                    f"  Fold {fold_idx}: SHAP not available, using permutation"
                )
                actual_method = "permutation"

            # Fallback to permutation if SHAP failed or was unavailable
            if importance_values is None and actual_method == "permutation":
                try:
                    importance_values = compute_importance_permutation(
                        best_pipe, X_test, y_test, scoring
                    )
                    actual_method = "permutation"
                except Exception as e:
                    log.error(
                        f"  Fold {fold_idx}: Permutation importance also failed ({e})"
                    )
                    continue

            if importance_values is None:
                # Last resort: try permutation
                try:
                    importance_values = compute_importance_permutation(
                        best_pipe, X_test, y_test, scoring
                    )
                    actual_method = "permutation"
                except Exception as e:
                    log.error(
                        f"  Fold {fold_idx}: All importance methods failed ({e})"
                    )
                    continue

            # Accumulate SHAP / permutation values
            if actual_method in ("shap_tree", "shap_linear") and shap_sum is not None:
                for i, subj_idx in enumerate(test_idx):
                    shap_sum[subj_idx] += importance_values[i]
                    shap_count[subj_idx] += 1
            else:
                # Permutation importance (primary or fallback)
                perm_importance_sum += importance_values
                n_perm_folds += 1

            # ----- Rank features for panel evaluation -----
            ranked_indices, importance_flat = get_importance_ranking(
                importance_values, actual_method
            )

            # Track top-10 stability
            top10_idx = ranked_indices[:10]
            if fold_idx < 50:
                fold_top10_tracker[fold_idx, top10_idx] = True

            # Record top-20 for this fold
            for rank_pos in range(min(20, len(ranked_indices))):
                fidx = ranked_indices[rank_pos]
                all_fold_top20_rows.append(
                    {
                        "target": target_name,
                        "fold": fold_idx,
                        "rank": rank_pos + 1,
                        "protein": protein_names[fidx],
                        "importance": float(importance_flat[fidx]),
                    }
                )

            # ----- 3. Panel size evaluation -----
            for k in panel_sizes_full:
                top_k_idx = ranked_indices[:k]
                metric_val = eval_panel(
                    X_train, y_train, X_test, y_test, top_k_idx, metric_name
                )
                all_panel_rows.append(
                    {
                        "target": target_name,
                        "panel_size": k,
                        "fold": fold_idx,
                        "metric_value": metric_val,
                        "metric_name": metric_name,
                        "is_random": False,
                        "draw_id": 0,
                    }
                )

            # ----- 4. NFASC + GDF15 baseline -----
            if len(benchmark_idx) > 0:
                bm_metric = eval_panel(
                    X_train, y_train, X_test, y_test, benchmark_idx, metric_name
                )
                all_panel_rows.append(
                    {
                        "target": target_name,
                        "panel_size": len(benchmark_idx),
                        "fold": fold_idx,
                        "metric_value": bm_metric,
                        "metric_name": metric_name,
                        "is_random": False,
                        "draw_id": -1,  # sentinel for NFASC+GDF15
                    }
                )

            # ----- 5. Random protein baselines -----
            for rk in RANDOM_PANEL_SIZES:
                for draw_id in range(N_RANDOM_DRAWS):
                    rng = np.random.RandomState(fold_idx * 1000 + draw_id)
                    random_idx = rng.choice(n_proteins, size=rk, replace=False)
                    rand_metric = eval_panel(
                        X_train,
                        y_train,
                        X_test,
                        y_test,
                        random_idx,
                        metric_name,
                    )
                    all_panel_rows.append(
                        {
                            "target": target_name,
                            "panel_size": rk,
                            "fold": fold_idx,
                            "metric_value": rand_metric,
                            "metric_name": metric_name,
                            "is_random": True,
                            "draw_id": draw_id,
                        }
                    )

            elapsed_fold = time.time() - t_fold_start
            if fold_idx % 5 == 0 or fold_idx == 49:
                log.info(
                    f"  Fold {fold_idx:2d}/50  "
                    f"{metric_name}={primary_metric:.3f}  "
                    f"method={actual_method}  "
                    f"({elapsed_fold:.1f}s)"
                )

        # ------------------------------------------------------------------
        # Post-fold aggregation for this target
        # ------------------------------------------------------------------
        elapsed_target = time.time() - t_target_start
        mean_metric = np.nanmean(fold_primary_metrics)
        sd_metric = np.nanstd(fold_primary_metrics)
        log.info(
            f"  {target_name} complete: mean {metric_name}="
            f"{mean_metric:.3f} ± {sd_metric:.3f}  ({elapsed_target:.0f}s)"
        )

        # Compute fold stability: fraction of 50 folds each protein was in top-10
        fold_stability = fold_top10_tracker.sum(axis=0) / 50.0  # (n_proteins,)

        # Compute aggregate importance per protein
        # Decide which accumulator to use: prefer SHAP if any folds used it
        n_shap_subjects = int(shap_count.sum()) if shap_count is not None else 0
        # Each fold has ~35 test subjects, so n_shap_folds ≈ n_shap_subjects / 35
        n_shap_folds_est = n_shap_subjects // max(1, n_subjects // 5) if n_shap_subjects > 0 else 0
        used_shap = (
            importance_method in ("shap_tree", "shap_linear")
            and shap_count is not None
            and n_shap_subjects > 0
        )
        log.info(
            f"  Importance folds: ~{n_shap_folds_est} SHAP ({n_shap_subjects} subject-appearances), "
            f"{n_perm_folds} permutation"
        )

        if used_shap:
            # Average SHAP across appearances (each subject seen ~10 times)
            safe_count = np.maximum(shap_count, 1)
            if importance_method == "shap_tree":
                shap_mean = shap_sum / safe_count[:, None]
                mean_abs_importance = np.abs(shap_mean).mean(axis=0)
            else:
                shap_mean = shap_sum / safe_count[:, None, None]
                mean_abs_importance = np.abs(shap_mean).mean(axis=(0, 2))
        else:
            # Permutation: average across folds that used it
            denom = max(n_perm_folds, 1)
            mean_abs_importance = perm_importance_sum / denom

        # Rank proteins
        rank_order = np.argsort(mean_abs_importance)[::-1]
        protein_ranks = np.empty(n_proteins, dtype=int)
        protein_ranks[rank_order] = np.arange(1, n_proteins + 1)

        for p_idx in range(n_proteins):
            all_importance_rows.append(
                {
                    "protein": protein_names[p_idx],
                    "target": target_name,
                    "mean_abs_shap": float(mean_abs_importance[p_idx]),
                    "rank": int(protein_ranks[p_idx]),
                    "fold_stability": float(fold_stability[p_idx]),
                }
            )

        # ------------------------------------------------------------------
        # Save per-target SHAP / importance .npz
        # ------------------------------------------------------------------
        npz_path = os.path.join(OUTDIR, f"226a_shap_{target_name}.npz")

        if used_shap and importance_method == "shap_tree":
            safe_count = np.maximum(shap_count, 1)
            shap_avg = shap_sum / safe_count[:, None]
            np.savez_compressed(
                npz_path,
                shap_values=shap_avg.astype(np.float32),
                sample_indices=np.arange(n_subjects),
                protein_names=protein_names,
            )
            log.info(f"  Saved SHAP values: {npz_path}  shape={shap_avg.shape}")

        elif used_shap and importance_method == "shap_linear":
            safe_count = np.maximum(shap_count, 1)
            shap_avg = shap_sum / safe_count[:, None, None]
            np.savez_compressed(
                npz_path,
                shap_values=shap_avg.astype(np.float32),
                sample_indices=np.arange(n_subjects),
                protein_names=protein_names,
            )
            log.info(f"  Saved SHAP values: {npz_path}  shape={shap_avg.shape}")

        else:
            # Permutation importance (primary or all-fallback)
            denom = max(n_perm_folds, 1)
            avg_perm = perm_importance_sum / denom
            np.savez_compressed(
                npz_path,
                importance_values=avg_perm.astype(np.float32),
                protein_names=protein_names,
            )
            log.info(
                f"  Saved permutation importance: {npz_path}  "
                f"shape={avg_perm.shape}"
            )

    # ------------------------------------------------------------------
    # Save aggregate CSVs
    # ------------------------------------------------------------------

    # 226a_protein_importance.csv
    df_importance = pd.DataFrame(all_importance_rows)
    imp_path = os.path.join(OUTDIR, "226a_protein_importance.csv")
    df_importance.to_csv(imp_path, index=False)
    log.info(f"Saved protein importance: {imp_path}  ({len(df_importance)} rows)")

    # 226a_panel_curves.csv
    df_panels = pd.DataFrame(all_panel_rows)
    panels_path = os.path.join(OUTDIR, "226a_panel_curves.csv")
    df_panels.to_csv(panels_path, index=False)
    log.info(f"Saved panel curves: {panels_path}  ({len(df_panels)} rows)")

    # 226a_fold_top20.csv
    df_top20 = pd.DataFrame(all_fold_top20_rows)
    top20_path = os.path.join(OUTDIR, "226a_fold_top20.csv")
    df_top20.to_csv(top20_path, index=False)
    log.info(f"Saved fold top-20: {top20_path}  ({len(df_top20)} rows)")

    # ------------------------------------------------------------------
    # Summary statistics
    # ------------------------------------------------------------------
    log.info("\n=== Summary ===")
    for target_name in target_configs:
        sub = df_importance[df_importance["target"] == target_name].sort_values(
            "rank"
        )
        top5 = sub.head(5)
        log.info(f"\n{target_name} — top 5 proteins:")
        for _, row in top5.iterrows():
            log.info(
                f"  #{int(row['rank']):3d}  {row['protein']:15s}  "
                f"|SHAP|={row['mean_abs_shap']:.4f}  "
                f"stability={row['fold_stability']:.2f}"
            )

        # Panel curve summary (non-random)
        sub_panels = df_panels[
            (df_panels["target"] == target_name)
            & (~df_panels["is_random"])
            & (df_panels["draw_id"] == 0)
        ]
        if not sub_panels.empty:
            panel_summary = (
                sub_panels.groupby("panel_size")["metric_value"]
                .agg(["mean", "std"])
                .reset_index()
            )
            log.info(f"  Panel sizes ({sub_panels['metric_name'].iloc[0]}):")
            for _, ps_row in panel_summary.iterrows():
                log.info(
                    f"    k={int(ps_row['panel_size']):5d}: "
                    f"{ps_row['mean']:.3f} ± {ps_row['std']:.3f}"
                )

        # NFASC+GDF15 baseline
        bm_panels = df_panels[
            (df_panels["target"] == target_name) & (df_panels["draw_id"] == -1)
        ]
        if not bm_panels.empty:
            bm_mean = bm_panels["metric_value"].mean()
            bm_std = bm_panels["metric_value"].std()
            log.info(f"  NFASC+GDF15 baseline: {bm_mean:.3f} ± {bm_std:.3f}")

        # Random baselines
        for rk in RANDOM_PANEL_SIZES:
            rand_sub = df_panels[
                (df_panels["target"] == target_name)
                & (df_panels["is_random"])
                & (df_panels["panel_size"] == rk)
            ]
            if not rand_sub.empty:
                # Average metric across draws within each fold, then across folds
                fold_means = rand_sub.groupby("fold")["metric_value"].mean()
                log.info(
                    f"  Random k={rk}: {fold_means.mean():.3f} ± {fold_means.std():.3f}"
                )

    elapsed_total = time.time() - t_start
    log.info(f"\nTotal elapsed: {elapsed_total:.0f}s ({elapsed_total/3600:.1f}h)")
    log.info("=== 226a complete ===")


if __name__ == "__main__":
    main()
