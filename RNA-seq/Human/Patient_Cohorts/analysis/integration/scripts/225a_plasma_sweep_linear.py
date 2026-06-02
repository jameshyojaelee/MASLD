"""
225a_plasma_sweep_linear.py
Plasma model sweep — Linear & Classical Models (8 models)

Elastic net, Ridge, Lasso, SGD, Bayesian Ridge, RF, Balanced RF, ExtraTrees.
10x repeated 5-fold nested CV with inner GridSearchCV.

Part of a 3-script parallel sweep (225a: linear/classical, 225b: boosting, 225c: deep/other).
All scripts use identical data loading and output to plasma_sweep_225{a,b,c}_*.csv.
Results are combined downstream for final model selection.
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
# Data loading (identical across all 3 sweep scripts)
# ---------------------------------------------------------------------------

# Load Olink (NaN preserved — imputation happens per CV fold)
olink = pd.read_csv(os.path.join(BASE, "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt"), sep="\t")
protein_names = olink["Assay"].values
npx = olink.drop("Assay", axis=1).T.values.astype(np.float32)  # 218 subjects x 1461 proteins

# Load metadata (subjects 1-177 have fibrosis staging)
meta = pd.read_csv(os.path.join(BASE, "Analysis/Proteomics/results/gse276114_disease_metadata.csv"))
meta["subject_idx"] = meta["sample_number"].astype(int) - 1  # 0-indexed
valid = meta[meta["subject_idx"] < npx.shape[0]].copy()
valid["y"] = (valid["disease_group"].isin(["F3", "F4"])).astype(int)  # 1=advanced, 0=early
valid["etiology"] = valid["disease"].values

X = npx[valid["subject_idx"].values]  # 177 x 1461
y = valid["y"].values
etiologies = valid["etiology"].values
log.info(f"Data: {X.shape[0]} subjects x {X.shape[1]} proteins, {y.sum()} advanced / {(1-y).sum()} early")

# ---------------------------------------------------------------------------
# Model definitions
# ---------------------------------------------------------------------------

from sklearn.linear_model import LogisticRegression, RidgeClassifier, SGDClassifier, BayesianRidge
from sklearn.ensemble import RandomForestClassifier, ExtraTreesClassifier
try:
    from imblearn.ensemble import BalancedRandomForestClassifier
    HAS_IMBLEARN = True
    log.info("imblearn available — BalancedRandomForestClassifier enabled")
except ImportError:
    HAS_IMBLEARN = False
    log.info("imblearn not available — BalancedRandomForestClassifier skipped")

MODELS = {
    "elastic_net": {
        "model": LogisticRegression(solver="saga", max_iter=5000, class_weight="balanced", random_state=SEED),
        "params": {"model__C": [0.001, 0.01, 0.1, 1, 10, 100], "model__l1_ratio": [0.1, 0.3, 0.5, 0.7, 0.9], "model__penalty": ["elasticnet"]}
    },
    "ridge": {
        "model": RidgeClassifier(class_weight="balanced"),
        "params": {"model__alpha": [0.001, 0.01, 0.1, 1, 10, 100, 1000]}
    },
    "lasso": {
        "model": LogisticRegression(penalty="l1", solver="liblinear", max_iter=5000, class_weight="balanced", random_state=SEED),
        "params": {"model__C": [0.001, 0.01, 0.1, 1, 10, 100]}
    },
    "sgd_elasticnet": {
        "model": SGDClassifier(loss="log_loss", penalty="elasticnet", class_weight="balanced", random_state=SEED, max_iter=2000),
        "params": {"model__alpha": [1e-5, 1e-4, 1e-3, 1e-2, 0.1], "model__l1_ratio": [0.1, 0.3, 0.5, 0.7, 0.9]}
    },
    "random_forest": {
        "model": RandomForestClassifier(class_weight="balanced", random_state=SEED),
        "params": {"model__n_estimators": [100, 500, 1000], "model__max_depth": [5, 10, None], "model__min_samples_leaf": [5, 10, 20]}
    },
    "extra_trees": {
        "model": ExtraTreesClassifier(class_weight="balanced", random_state=SEED),
        "params": {"model__n_estimators": [100, 500], "model__max_depth": [5, 10, None], "model__min_samples_leaf": [5, 10]}
    },
}

# Bayesian Ridge wrapper for classification
from sklearn.base import BaseEstimator, ClassifierMixin

class BayesianRidgeClassifier(BaseEstimator, ClassifierMixin):
    def __init__(self):
        self.model = BayesianRidge()

    def fit(self, X, y):
        self.model.fit(X, y)
        self.classes_ = np.array([0, 1])
        return self

    def predict(self, X):
        return (self.model.predict(X) >= 0.5).astype(int)

    def predict_proba(self, X):
        p = np.clip(self.model.predict(X), 0, 1)
        return np.column_stack([1 - p, p])

MODELS["bayesian_ridge"] = {
    "model": BayesianRidgeClassifier(),
    "params": {}  # Auto-tuned internally
}

if HAS_IMBLEARN:
    MODELS["balanced_rf"] = {
        "model": BalancedRandomForestClassifier(random_state=SEED),
        "params": {"model__n_estimators": [100, 500], "model__max_depth": [5, 10, None], "model__min_samples_leaf": [5, 10]}
    }

# ---------------------------------------------------------------------------
# Nested CV
# ---------------------------------------------------------------------------

def run_nested_cv(X, y, model_name, model, param_grid, n_repeats=10, n_splits=5):
    """10x repeated 5-fold nested CV with inner GridSearchCV."""
    outer_cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=SEED)
    inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=SEED)

    aurocs, f1s, bal_accs = [], [], []
    best_params_list = []

    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        X_train_raw, X_test_raw = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Per-fold imputation: statistics from training data only
        X_train, X_test = impute_train_fold(X_train_raw, X_test_raw)

        pipe = Pipeline([("scaler", StandardScaler()), ("model", model)])

        if param_grid:
            grid = GridSearchCV(pipe, param_grid, cv=inner_cv, scoring="roc_auc", n_jobs=4, refit=True)
            grid.fit(X_train, y_train)
            best_pipe = grid.best_estimator_
            best_params_list.append(grid.best_params_)
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
        f1 = f1_score(y_test, y_pred)
        bal_acc = balanced_accuracy_score(y_test, y_pred)

        aurocs.append(auroc)
        f1s.append(f1)
        bal_accs.append(bal_acc)

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
    log.info("=== 225a: Plasma Model Sweep — Linear & Classical ===")
    log.info(f"Models to run: {list(MODELS.keys())}")

    results = []
    for name, cfg in MODELS.items():
        log.info(f"Running {name}...")
        t0 = time.time()
        try:
            from sklearn.base import clone
            model = clone(cfg["model"])
            res = run_nested_cv(X, y, name, model, cfg["params"])
            elapsed = time.time() - t0
            res["elapsed_sec"] = round(elapsed, 1)
            results.append(res)
            log.info(f"  {name}: AUROC={res['mean_auroc']:.3f} ± {res['sd_auroc']:.3f}  ({elapsed:.0f}s)")
        except Exception as e:
            log.error(f"  {name} FAILED: {e}")
            results.append({"model": name, "mean_auroc": np.nan, "sd_auroc": np.nan,
                            "mean_f1": np.nan, "mean_bal_acc": np.nan, "n_folds": 0,
                            "elapsed_sec": np.nan})

    df = pd.DataFrame(results)
    out_path = os.path.join(OUTDIR, "plasma_sweep_225a_linear.csv")
    df.to_csv(out_path, index=False)
    log.info(f"\nResults:\n{df.to_string()}")
    log.info(f"Saved to {out_path}")
    log.info("Done.")
