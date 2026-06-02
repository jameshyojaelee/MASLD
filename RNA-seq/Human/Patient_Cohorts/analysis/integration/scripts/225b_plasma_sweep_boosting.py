"""
225b_plasma_sweep_boosting.py
Plasma biomarker model sweep — Boosting & SVM models (10 models)

Models: GradientBoosting, HistGradientBoosting, XGBoost, LightGBM,
        LinearSVC (SVC kernel=linear), RBF-SVM, NuSVC, AdaBoost,
        EasyEnsemble, BalancedBagging

Nested CV: 10-repeat x 5-fold outer, 3-fold inner GridSearchCV.
Output: plasma_sweep_225b_boosting.csv in results/multiprogram/
"""

import os, sys, time, warnings, logging, importlib.util
import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold, GridSearchCV, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, f1_score, balanced_accuracy_score
from sklearn.pipeline import Pipeline
warnings.filterwarnings("ignore")

# Import per-fold imputation helper from common module
_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("plasma_common", os.path.join(_here, "225_plasma_common.py"))
_pc = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_pc)
impute_train_fold = _pc.impute_train_fold

SEED = 42
np.random.seed(SEED)

BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram")
os.makedirs(OUTDIR, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Load Olink
# ---------------------------------------------------------------------------
olink = pd.read_csv(os.path.join(BASE, "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt"), sep="\t")
protein_names = olink["Assay"].values
npx = olink.drop("Assay", axis=1).T.values.astype(np.float32)  # NaN preserved for per-fold imputation

meta = pd.read_csv(os.path.join(BASE, "Analysis/Proteomics/results/gse276114_disease_metadata.csv"))
meta["subject_idx"] = meta["sample_number"].astype(int) - 1
valid = meta[meta["subject_idx"] < npx.shape[0]].copy()
valid["y"] = (valid["disease_group"].isin(["F3", "F4"])).astype(int)
valid["etiology"] = valid["disease"].values

X = npx[valid["subject_idx"].values]
y = valid["y"].values
etiologies = valid["etiology"].values
log.info(f"Data: {X.shape[0]} subjects x {X.shape[1]} proteins, {y.sum()} advanced / {(1-y).sum()} early")

# ---------------------------------------------------------------------------
# Model definitions
# ---------------------------------------------------------------------------
from sklearn.ensemble import GradientBoostingClassifier, HistGradientBoostingClassifier, AdaBoostClassifier
from sklearn.svm import SVC, LinearSVC, NuSVC
import xgboost as xgb
import lightgbm as lgb
try:
    from imblearn.ensemble import EasyEnsembleClassifier, BalancedBaggingClassifier
    HAS_IMBLEARN = True
    log.info("imbalanced-learn available — EasyEnsemble and BalancedBagging enabled")
except ImportError:
    HAS_IMBLEARN = False
    log.info("imbalanced-learn not available — EasyEnsemble and BalancedBagging skipped")

scale_pos = float(np.sum(y == 0)) / np.sum(y == 1)  # ~0.283

MODELS = {
    "gradient_boosting": {
        "model": GradientBoostingClassifier(random_state=SEED),
        "params": {
            "model__n_estimators": [100, 300],
            "model__learning_rate": [0.01, 0.1],
            "model__max_depth": [3, 5],
            "model__subsample": [0.8],
            "model__min_samples_leaf": [10],
        }
    },
    "hist_gradient_boosting": {
        "model": HistGradientBoostingClassifier(
            class_weight="balanced", random_state=SEED,
            early_stopping=True, validation_fraction=0.15),
        "params": {
            "model__max_iter": [100, 300],
            "model__max_depth": [3, 5, 7],
            "model__learning_rate": [0.01, 0.1],
            "model__l2_regularization": [0, 1, 10],
        }
    },
    "xgboost": {
        "model": xgb.XGBClassifier(
            scale_pos_weight=scale_pos, random_state=SEED,
            eval_metric="logloss", use_label_encoder=False),
        "params": {
            "model__learning_rate": [0.01, 0.05, 0.1],
            "model__max_depth": [3, 5, 7],
            "model__n_estimators": [100, 300],
            "model__subsample": [0.6, 0.8],
            "model__colsample_bytree": [0.5, 0.8],
            "model__reg_lambda": [1, 10],
        }
    },
    "lightgbm": {
        "model": lgb.LGBMClassifier(is_unbalance=True, random_state=SEED, verbose=-1, n_jobs=1),
        "params": {
            "model__num_leaves": [15, 31],
            "model__learning_rate": [0.01, 0.1],
            "model__n_estimators": [100, 300],
            "model__reg_lambda": [1, 10],
            "model__feature_fraction": [0.5, 0.8],
        }
    },
    "linear_svc": {
        "model": SVC(kernel="linear", class_weight="balanced", probability=True, random_state=SEED),
        "params": {
            "model__C": [0.01, 0.1, 1, 10, 100],
        }
    },
    "rbf_svm": {
        "model": SVC(kernel="rbf", class_weight="balanced", probability=True, random_state=SEED),
        "params": {
            "model__C": [0.1, 1, 10, 100],
            "model__gamma": ["scale", "auto", 0.001, 0.01],
        }
    },
    "nu_svc": {
        "model": NuSVC(kernel="rbf", class_weight="balanced", probability=True, random_state=SEED),
        "params": {
            "model__nu": [0.1, 0.3, 0.5],
            "model__gamma": ["scale", "auto"],
        }
    },
    "adaboost": {
        "model": AdaBoostClassifier(random_state=SEED),
        "params": {
            "model__n_estimators": [50, 100, 200],
            "model__learning_rate": [0.1, 0.5, 1.0],
        }
    },
}

if HAS_IMBLEARN:
    MODELS["easy_ensemble"] = {
        "model": EasyEnsembleClassifier(random_state=SEED),
        "params": {
            "model__n_estimators": [10, 20, 50],
        }
    }
    from sklearn.linear_model import LogisticRegression
    MODELS["balanced_bagging"] = {
        "model": BalancedBaggingClassifier(
            estimator=LogisticRegression(max_iter=2000, random_state=SEED),
            random_state=SEED),
        "params": {
            "model__n_estimators": [20, 50, 100],
        }
    }

# ---------------------------------------------------------------------------
# Nested CV
# ---------------------------------------------------------------------------
def run_nested_cv(X, y, model_name, model, param_grid, n_repeats=10, n_splits=5):
    outer_cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=SEED)
    inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=SEED)

    aurocs, f1s, bal_accs = [], [], []

    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        X_train_raw, X_test_raw = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Per-fold imputation: statistics from training data only
        X_train, X_test = impute_train_fold(X_train_raw, X_test_raw)

        pipe = Pipeline([("scaler", StandardScaler()), ("model", model)])

        if param_grid:
            # LightGBM deadlocks with parallel joblib workers — use n_jobs=1 for lgb
            _gn = 1 if "lightgbm" in str(type(model)).lower() or "lgbm" in str(type(model)).lower() else 4
            grid = GridSearchCV(pipe, param_grid, cv=inner_cv, scoring="roc_auc", n_jobs=_gn, refit=True)
            grid.fit(X_train, y_train)
            best_pipe = grid.best_estimator_
        else:
            pipe.fit(X_train, y_train)
            best_pipe = pipe

        if hasattr(best_pipe, "predict_proba"):
            y_prob = best_pipe.predict_proba(X_test)[:, 1]
            auroc = roc_auc_score(y_test, y_prob)
        elif hasattr(best_pipe, "decision_function"):
            y_dec = best_pipe.decision_function(X_test)
            auroc = roc_auc_score(y_test, y_dec)
        else:
            auroc = np.nan

        y_pred = best_pipe.predict(X_test)
        f1s.append(f1_score(y_test, y_pred))
        bal_accs.append(balanced_accuracy_score(y_test, y_pred))
        aurocs.append(auroc)

    return {
        "model": model_name,
        "mean_auroc": np.nanmean(aurocs), "sd_auroc": np.nanstd(aurocs),
        "mean_f1": np.mean(f1s), "mean_bal_acc": np.mean(bal_accs),
        "n_folds": len(aurocs),
    }

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # Support MODEL_SUBSET env var for parallel splitting (comma-separated model names)
    _subset_str = os.environ.get("MODEL_SUBSET", "")
    if _subset_str:
        _subset = set(_subset_str.split(","))
        MODELS = {k: v for k, v in MODELS.items() if k in _subset}
    # Support SKIP_MODELS env var to skip already-completed models
    _skip_str = os.environ.get("SKIP_MODELS", "")
    if _skip_str:
        _skip = set(_skip_str.split(","))
        MODELS = {k: v for k, v in MODELS.items() if k not in _skip}

    # Support custom output suffix for parallel runs
    _out_suffix = os.environ.get("OUT_SUFFIX", "boosting")

    log.info("=== 225b: Plasma Model Sweep — Boosting & SVM ===")
    log.info(f"Models to run: {list(MODELS.keys())}")
    results = []
    t0_total = time.time()

    for name, cfg in MODELS.items():
        log.info(f"Running {name}...")
        t0 = time.time()
        try:
            from sklearn.base import clone
            model = clone(cfg["model"])
            res = run_nested_cv(X, y, name, model, cfg["params"])
            results.append(res)
            elapsed = time.time() - t0
            log.info(f"  {name}: AUROC={res['mean_auroc']:.3f} ± {res['sd_auroc']:.3f}  ({elapsed:.1f}s)")
        except Exception as e:
            log.error(f"  {name} FAILED: {e}")
            results.append({
                "model": name, "mean_auroc": np.nan, "sd_auroc": np.nan,
                "mean_f1": np.nan, "mean_bal_acc": np.nan, "n_folds": 0,
            })

    df = pd.DataFrame(results)
    out_path = os.path.join(OUTDIR, f"plasma_sweep_225b_{_out_suffix}.csv")
    df.to_csv(out_path, index=False)
    total_elapsed = time.time() - t0_total
    log.info(f"\nResults:\n{df.to_string()}")
    log.info(f"Saved: {out_path}")
    log.info(f"Total elapsed: {total_elapsed:.1f}s")
    log.info("Done.")
