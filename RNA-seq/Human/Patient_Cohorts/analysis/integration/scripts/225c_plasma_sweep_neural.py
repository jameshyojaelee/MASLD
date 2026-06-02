"""
225c_plasma_sweep_neural.py
Plasma biomarker model sweep — Neural, AutoML & Ensemble models (7 models)

Models: MLP, TabNet, TabPFN, FLAML (2 time budgets), VotingClassifier, StackingClassifier.
Nested CV: 10-repeat x 5-fold outer, 3-fold inner GridSearchCV for Pipeline-compatible models.
Non-Pipeline models (TabNet, TabPFN, FLAML) run directly in outer CV loop.
TabPFN uses top-100 variance features (input-size constraint).

Part of a 3-script parallel sweep (225a: linear/classical, 225b: boosting, 225c: neural/ensemble).
All scripts use identical data loading and output to plasma_sweep_225{a,b,c}_*.csv.
Results are combined downstream for final model selection.

Output: plasma_sweep_225c_neural.csv in results/multiprogram/
"""

import os, sys, time, warnings, logging, importlib.util
import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold, GridSearchCV, StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, f1_score, balanced_accuracy_score
from sklearn.pipeline import Pipeline
from sklearn.base import clone
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

olink = pd.read_csv(
    os.path.join(BASE, "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt"),
    sep="\t"
)
protein_names = olink["Assay"].values
npx = olink.drop("Assay", axis=1).T.values.astype(np.float32)  # 218 subjects x 1461 proteins
# NaN preserved — imputation happens per CV fold

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
# Optional dependency imports
# ---------------------------------------------------------------------------

# TabNet
try:
    from pytorch_tabnet.tab_model import TabNetClassifier as _TabNet
    import torch
    HAS_TABNET = True
    log.info("pytorch-tabnet available — TabNet enabled")
except ImportError:
    HAS_TABNET = False
    log.info("pytorch-tabnet not available — TabNet skipped")

# TabPFN
try:
    from tabpfn import TabPFNClassifier
    HAS_TABPFN = True
    log.info("tabpfn available — TabPFN enabled")
except ImportError:
    HAS_TABPFN = False
    log.info("tabpfn not available — TabPFN skipped")

# FLAML
try:
    from flaml import AutoML
    HAS_FLAML = True
    log.info("flaml available — FLAML enabled")
except ImportError:
    HAS_FLAML = False
    log.info("flaml not available — FLAML skipped")

# ---------------------------------------------------------------------------
# TabNet sklearn-compatible wrapper
# ---------------------------------------------------------------------------

if HAS_TABNET:
    from sklearn.base import BaseEstimator, ClassifierMixin

    class TabNetWrapper(BaseEstimator, ClassifierMixin):
        """Sklearn-compatible wrapper around pytorch-tabnet's TabNetClassifier.

        Does NOT fit into sklearn Pipeline (requires numpy array, not pipeline steps).
        Handles its own train/validation split for early stopping.
        """

        def __init__(self, n_steps=3, n_a=8, n_d=8, random_state=42):
            self.n_steps = n_steps
            self.n_a = n_a
            self.n_d = n_d
            self.random_state = random_state
            self.model = None
            self.classes_ = np.array([0, 1])

        def fit(self, X, y):
            self.model = _TabNet(
                n_steps=self.n_steps,
                n_a=self.n_a,
                n_d=self.n_d,
                seed=self.random_state,
                verbose=0,
                device_name="cpu",
            )
            # Internal validation split for early stopping
            n_val = max(int(0.15 * len(y)), 5)
            rng = np.random.RandomState(self.random_state)
            idx = np.arange(len(y))
            rng.shuffle(idx)
            val_idx, train_idx = idx[:n_val], idx[n_val:]
            self.model.fit(
                X[train_idx].astype(np.float32), y[train_idx],
                eval_set=[(X[val_idx].astype(np.float32), y[val_idx])],
                max_epochs=100,
                patience=15,
                batch_size=64,
            )
            return self

        def predict(self, X):
            return self.model.predict(X.astype(np.float32))

        def predict_proba(self, X):
            return self.model.predict_proba(X.astype(np.float32))

        def get_params(self, deep=True):
            return {
                "n_steps": self.n_steps,
                "n_a": self.n_a,
                "n_d": self.n_d,
                "random_state": self.random_state,
            }

        def set_params(self, **params):
            for k, v in params.items():
                setattr(self, k, v)
            return self

# ---------------------------------------------------------------------------
# FLAML sklearn-compatible wrapper
# ---------------------------------------------------------------------------

if HAS_FLAML:
    from sklearn.base import BaseEstimator, ClassifierMixin

    class FLAMLWrapper(BaseEstimator, ClassifierMixin):
        """Sklearn-compatible wrapper around FLAML AutoML.

        Does NOT fit into sklearn Pipeline. Runs AutoML with fixed time budget.
        """

        def __init__(self, time_budget=30, random_state=42):
            self.time_budget = time_budget
            self.random_state = random_state
            self.model = None
            self.classes_ = np.array([0, 1])

        def fit(self, X, y):
            self.model = AutoML()
            self.model.fit(
                X, y,
                task="classification",
                metric="roc_auc",
                time_budget=self.time_budget,
                verbose=0,
                seed=self.random_state,
            )
            return self

        def predict(self, X):
            return self.model.predict(X)

        def predict_proba(self, X):
            return self.model.predict_proba(X)

        def get_params(self, deep=True):
            return {"time_budget": self.time_budget, "random_state": self.random_state}

        def set_params(self, **params):
            for k, v in params.items():
                setattr(self, k, v)
            return self

# ---------------------------------------------------------------------------
# Model definitions
# ---------------------------------------------------------------------------

from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import (
    RandomForestClassifier,
    HistGradientBoostingClassifier,
    VotingClassifier,
    StackingClassifier,
)
from sklearn.linear_model import LogisticRegression

# Model 16: MLP — Pipeline-compatible, inner GridSearchCV
MODELS = {
    "mlp": {
        "model": MLPClassifier(
            random_state=SEED,
            max_iter=500,
            early_stopping=True,
            validation_fraction=0.15,
        ),
        "params": {
            "model__hidden_layer_sizes": [(64,), (128, 64), (256, 128)],
            "model__alpha": [0.001, 0.01, 0.1],
        },
        "is_pipeline": True,
    },
}

# Model 17: TabNet — non-Pipeline, fixed config
if HAS_TABNET:
    MODELS["tabnet"] = {
        "model": TabNetWrapper(n_steps=3, n_a=8, n_d=8, random_state=SEED),
        "params": {},
        "is_pipeline": False,
    }

# Model 18: TabPFN — handled separately in main (feature-limited)
# (added to results list directly, not to MODELS dict)

# Models 22-23: FLAML with two time budgets — non-Pipeline
if HAS_FLAML:
    MODELS["flaml_30s"] = {
        "model": FLAMLWrapper(time_budget=30, random_state=SEED),
        "params": {},
        "is_pipeline": False,
    }
    MODELS["flaml_60s"] = {
        "model": FLAMLWrapper(time_budget=60, random_state=SEED),
        "params": {},
        "is_pipeline": False,
    }

# Model 24: VotingClassifier (soft) — Pipeline-compatible, no inner tuning
_lr_base = LogisticRegression(
    penalty="l2", C=1.0, class_weight="balanced", max_iter=5000, random_state=SEED
)
_rf_base = RandomForestClassifier(
    n_estimators=500, class_weight="balanced", random_state=SEED
)
_hgb_base = HistGradientBoostingClassifier(
    class_weight="balanced", random_state=SEED
)

MODELS["voting_lr_rf_hgb"] = {
    "model": VotingClassifier(
        estimators=[("lr", _lr_base), ("rf", _rf_base), ("hgb", _hgb_base)],
        voting="soft",
    ),
    "params": {},
    "is_pipeline": True,
}

# Model 25: StackingClassifier — Pipeline-compatible, no inner tuning
MODELS["stacking_lr_rf_hgb"] = {
    "model": StackingClassifier(
        estimators=[
            ("lr", LogisticRegression(
                penalty="l2", C=1.0, class_weight="balanced", max_iter=5000, random_state=SEED
            )),
            ("rf", RandomForestClassifier(
                n_estimators=500, class_weight="balanced", random_state=SEED
            )),
            ("hgb", HistGradientBoostingClassifier(
                class_weight="balanced", random_state=SEED
            )),
        ],
        final_estimator=LogisticRegression(
            class_weight="balanced", max_iter=5000, random_state=SEED
        ),
        cv=3,
    ),
    "params": {},
    "is_pipeline": True,
}

# ---------------------------------------------------------------------------
# Nested CV — handles both Pipeline and non-Pipeline models
# ---------------------------------------------------------------------------

def run_nested_cv(X, y, model_name, model, param_grid, is_pipeline=True,
                  n_repeats=10, n_splits=5):
    """10x repeated 5-fold nested CV.

    Pipeline-compatible models: inner GridSearchCV with StandardScaler.
    Non-Pipeline models (TabNet, FLAML): scale externally, fit/predict directly.
    """
    outer_cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=SEED)
    inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=SEED)

    aurocs, f1s, bal_accs = [], [], []

    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        X_train_raw, X_test_raw = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Per-fold imputation: statistics from training data only
        X_train, X_test = impute_train_fold(X_train_raw, X_test_raw)

        if is_pipeline:
            # Build pipeline; StandardScaler inside pipeline so it's tuned correctly
            pipe = Pipeline([("scaler", StandardScaler()), ("model", clone(model))])

            if param_grid:
                grid = GridSearchCV(
                    pipe, param_grid, cv=inner_cv,
                    scoring="roc_auc", n_jobs=4, refit=True
                )
                grid.fit(X_train, y_train)
                best = grid.best_estimator_
            else:
                pipe.fit(X_train, y_train)
                best = pipe

            if hasattr(best, "predict_proba"):
                y_prob = best.predict_proba(X_test)[:, 1]
            else:
                y_prob = best.decision_function(X_test)
            y_pred = best.predict(X_test)

        else:
            # Non-Pipeline: scale manually, then fit the raw model
            scaler = StandardScaler()
            X_train_s = scaler.fit_transform(X_train)
            X_test_s = scaler.transform(X_test)

            try:
                m = clone(model) if hasattr(model, "get_params") else model
                m.fit(X_train_s, y_train)
                y_prob = m.predict_proba(X_test_s)[:, 1]
                y_pred = m.predict(X_test_s)
            except Exception as e:
                log.warning(f"  Fold {fold_idx} failed for {model_name}: {e}")
                aurocs.append(np.nan)
                f1s.append(np.nan)
                bal_accs.append(np.nan)
                continue

        aurocs.append(roc_auc_score(y_test, y_prob))
        f1s.append(f1_score(y_test, y_pred))
        bal_accs.append(balanced_accuracy_score(y_test, y_pred))

    return {
        "model": model_name,
        "mean_auroc": np.nanmean(aurocs),
        "sd_auroc": np.nanstd(aurocs),
        "mean_f1": np.nanmean(f1s),
        "mean_bal_acc": np.nanmean(bal_accs),
        "n_folds": int(np.sum(~np.isnan(aurocs))),
    }

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    log.info("=== 225c: Plasma Model Sweep — Neural, AutoML & Ensemble ===")
    log.info(f"Models in MODELS dict: {list(MODELS.keys())}")
    if HAS_TABPFN:
        log.info("TabPFN will run separately (top-100 variance features)")

    results = []
    t0_total = time.time()

    # --- Standard models (Pipeline + non-Pipeline) ---
    for name, cfg in MODELS.items():
        log.info(f"Running {name}...")
        t0 = time.time()
        try:
            res = run_nested_cv(
                X, y,
                model_name=name,
                model=cfg["model"],
                param_grid=cfg["params"],
                is_pipeline=cfg["is_pipeline"],
            )
            elapsed = time.time() - t0
            res["elapsed_sec"] = round(elapsed, 1)
            results.append(res)
            log.info(f"  {name}: AUROC={res['mean_auroc']:.3f} ± {res['sd_auroc']:.3f}  ({elapsed:.1f}s)")
        except Exception as e:
            log.error(f"  {name} FAILED: {e}")
            results.append({
                "model": name,
                "mean_auroc": np.nan,
                "sd_auroc": np.nan,
                "mean_f1": np.nan,
                "mean_bal_acc": np.nan,
                "n_folds": 0,
                "elapsed_sec": np.nan,
            })

    # --- TabPFN: top-100 variance features, no inner tuning (pretrained) ---
    # Feature selection by variance must use training data only per fold
    # to avoid data leakage.
    if HAS_TABPFN:
        log.info("Running tabpfn_top100 (top-100 features by variance, per-fold selection)...")
        t0 = time.time()
        try:
            tabpfn_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=SEED)
            tabpfn_aurocs, tabpfn_f1s, tabpfn_bal_accs = [], [], []

            for fold_idx, (train_idx, test_idx) in enumerate(tabpfn_cv.split(X, y)):
                X_train_raw, X_test_raw = X[train_idx], X[test_idx]
                y_train, y_test = y[train_idx], y[test_idx]

                # Per-fold imputation
                X_train_imp, X_test_imp = impute_train_fold(X_train_raw, X_test_raw)

                # Feature selection from training fold only
                var = X_train_imp.var(axis=0)
                top100 = np.argsort(var)[-100:]
                X_tr_sub = X_train_imp[:, top100]
                X_te_sub = X_test_imp[:, top100]

                scaler = StandardScaler()
                X_tr_sc = scaler.fit_transform(X_tr_sub)
                X_te_sc = scaler.transform(X_te_sub)

                m = clone(TabPFNClassifier(device="cpu", N_ensemble_configurations=32))
                m.fit(X_tr_sc, y_train)
                y_prob = m.predict_proba(X_te_sc)[:, 1]
                y_pred = m.predict(X_te_sc)

                tabpfn_aurocs.append(roc_auc_score(y_test, y_prob))
                tabpfn_f1s.append(f1_score(y_test, y_pred))
                tabpfn_bal_accs.append(balanced_accuracy_score(y_test, y_pred))

            elapsed = time.time() - t0
            res = {
                "model": "tabpfn_top100",
                "mean_auroc": np.nanmean(tabpfn_aurocs),
                "sd_auroc": np.nanstd(tabpfn_aurocs),
                "mean_f1": np.nanmean(tabpfn_f1s),
                "mean_bal_acc": np.nanmean(tabpfn_bal_accs),
                "n_folds": len(tabpfn_aurocs),
                "elapsed_sec": round(elapsed, 1),
            }
            results.append(res)
            log.info(f"  tabpfn_top100: AUROC={res['mean_auroc']:.3f} ± {res['sd_auroc']:.3f}  ({elapsed:.1f}s)")
        except Exception as e:
            log.error(f"  tabpfn_top100 FAILED: {e}")
            results.append({
                "model": "tabpfn_top100",
                "mean_auroc": np.nan,
                "sd_auroc": np.nan,
                "mean_f1": np.nan,
                "mean_bal_acc": np.nan,
                "n_folds": 0,
                "elapsed_sec": np.nan,
            })
    else:
        log.info("TabPFN not available — skipping tabpfn_top100")

    # --- Save results ---
    df = pd.DataFrame(results)
    out_path = os.path.join(OUTDIR, "plasma_sweep_225c_neural.csv")
    df.to_csv(out_path, index=False)

    total_elapsed = time.time() - t0_total
    log.info(f"\nResults:\n{df.to_string()}")
    log.info(f"Saved: {out_path}")
    log.info(f"Total elapsed: {total_elapsed:.1f}s")
    log.info("Done.")
