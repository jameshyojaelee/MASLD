"""
225_plasma_common.py
Shared plasma proteomics infrastructure for scripts 225d, 225e, and 225f.

Provides:
  - load_olink_data()       : Olink NPX matrix + multi-target labels
  - load_ms_data()          : DIA-MS PXD052937 matrix + metadata
  - OrdinalClassifierWrapper: Cumulative logit ordinal regression
  - get_model_registry()    : 35 models × 5 target types
  - run_nested_cv()         : Generalised nested CV (all target types)
  - compute_metrics()       : Target-type-appropriate metric dict
"""

import os
import sys
import logging
import warnings

# Set OMP_NUM_THREADS=1 BEFORE importing torch/tabnet to prevent OpenMP
# deadlock between PyTorch and LightGBM (both use OpenMP thread pools).
# Must be set before any library import that initialises OpenMP.
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, ClassifierMixin, RegressorMixin, clone
from sklearn.calibration import CalibratedClassifierCV
from sklearn.ensemble import (
    AdaBoostClassifier,
    ExtraTreesClassifier,
    GradientBoostingClassifier,
    HistGradientBoostingClassifier,
    RandomForestClassifier,
    VotingClassifier,
    StackingClassifier,
)
from sklearn.gaussian_process import GaussianProcessClassifier
from sklearn.gaussian_process.kernels import RBF
from sklearn.linear_model import (
    BayesianRidge,
    ElasticNet,
    Lasso,
    LogisticRegression,
    Ridge,
    RidgeClassifier,
    SGDClassifier,
)
from sklearn.metrics import (
    balanced_accuracy_score,
    cohen_kappa_score,
    f1_score,
    make_scorer,
    mean_absolute_error,
    roc_auc_score,
)
from sklearn.model_selection import (
    GridSearchCV,
    RepeatedKFold,
    RepeatedStratifiedKFold,
    StratifiedKFold,
    KFold,
)
from sklearn.multiclass import OneVsRestClassifier
from sklearn.multioutput import MultiOutputClassifier
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC, LinearSVC, NuSVC, SVR
from sklearn.kernel_ridge import KernelRidge
from scipy.stats import spearmanr

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
# Constants
# ---------------------------------------------------------------------------

SEED = 42
np.random.seed(SEED)

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)

OLINK_PATH = os.path.join(
    BASE,
    "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt",
)
OLINK_META_PATH = os.path.join(
    BASE, "Analysis/Proteomics/results/gse276114_disease_metadata.csv"
)
MS_MATRIX_PATH = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram/ms_plasma_matrix.csv",
)
# G1-002 fix (2026-06-01): the old pxd052937_disease_metadata.csv is a SEPARATE
# Spectronaut export (T-series file names) with ZERO overlapping sample IDs vs the
# matrix (P-series IDs); load_ms_data() formerly joined them positionally. The
# companion metadata below lives in the SAME export directory as the matrix and its
# sample_id matches the matrix index 100% (verified A/B/C/D counts: 7/17/38/10).
MS_META_PATH = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram/ms_plasma_metadata.csv",
)

# ---------------------------------------------------------------------------
# Optional dependency imports (graceful)
# ---------------------------------------------------------------------------

try:
    import xgboost as xgb
    HAS_XGB = True
    log.info("xgboost available")
except ImportError:
    HAS_XGB = False
    log.info("xgboost not available — XGBoost skipped")

try:
    import lightgbm as lgb
    HAS_LGB = True
    log.info("lightgbm available")
except ImportError:
    HAS_LGB = False
    log.info("lightgbm not available — LightGBM skipped")

try:
    from catboost import CatBoostClassifier, CatBoostRegressor
    HAS_CATBOOST = True
    log.info("catboost available")
except ImportError:
    HAS_CATBOOST = False
    log.info("catboost not available — CatBoost skipped")

try:
    import mord
    HAS_MORD = True
    log.info("mord available — OrdinalForest/CORAL enabled")
except ImportError:
    HAS_MORD = False
    log.info("mord not available — mord ordinal models skipped")

try:
    from tabpfn import TabPFNClassifier
    HAS_TABPFN = True
    log.info("tabpfn available")
except ImportError:
    HAS_TABPFN = False
    log.info("tabpfn not available — TabPFN skipped")

try:
    from flaml import AutoML as _FLAMLAutoML
    HAS_FLAML = True
    log.info("flaml available")
except ImportError:
    HAS_FLAML = False
    log.info("flaml not available — FLAML skipped")

try:
    from pytorch_tabnet.tab_model import TabNetClassifier as _TabNet
    import torch
    HAS_TABNET = True
    log.info("pytorch-tabnet available")
except ImportError:
    HAS_TABNET = False
    log.info("pytorch-tabnet not available — TabNet skipped")

try:
    from imblearn.ensemble import (
        BalancedRandomForestClassifier,
        EasyEnsembleClassifier,
        BalancedBaggingClassifier,
    )
    from imblearn.over_sampling import SMOTE
    HAS_IMBLEARN = True
    log.info("imbalanced-learn available")
except ImportError:
    HAS_IMBLEARN = False
    log.info("imbalanced-learn not available — imblearn models skipped")

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------


def impute_train_fold(X_train, X_test):
    """Compute column medians from X_train, apply to both splits.

    Prevents data leakage: imputation statistics are derived from training
    data only, then used to fill NaN in both training and test arrays.

    Parameters
    ----------
    X_train : ndarray (n_train, p) — may contain NaN
    X_test  : ndarray (n_test, p)  — may contain NaN

    Returns
    -------
    X_train_imp : ndarray (n_train, p) — NaN filled with train medians
    X_test_imp  : ndarray (n_test, p)  — NaN filled with train medians
    """
    X_train_imp = X_train.copy()
    X_test_imp = X_test.copy()
    for j in range(X_train.shape[1]):
        med = float(np.nanmedian(X_train[:, j]))
        train_mask = np.isnan(X_train_imp[:, j])
        test_mask = np.isnan(X_test_imp[:, j])
        if train_mask.any():
            X_train_imp[train_mask, j] = med
        if test_mask.any():
            X_test_imp[test_mask, j] = med
    return X_train_imp, X_test_imp


def load_olink_data():
    """Load Olink plasma data and metadata.

    Returns
    -------
    X : ndarray (177, 1461) float32 — raw values, may contain NaN
        (use impute_train_fold() inside CV loops to avoid leakage)
    protein_names : list of 1461 strings
    meta : DataFrame with geo_accession, disease_group, disease, sample_number
    y_targets : dict of arrays:
        "binary"          — (177,) int, F0-2=0, F3/F4=1
        "ordinal"         — (177,) int, F0-2=0, F3=1, F4=2
        "continuous"      — (177,) float, same encoding as ordinal
        "etiology_3"      — (177,) int, MASLD=0, CVH=1, ARLD=2
        "etiology_binary" — (177,) int, MASLD=1 else 0
        "etiologies"      — (177,) str, raw etiology label (for LOEO)
    """
    # ------------------------------------------------------------------
    # QUARANTINED 2026-06-01 — DO NOT USE for any reported result.
    # The label assignment below joins the Olink PLASMA matrix (anonymized
    # "Subject 1..218") to GSE276114 LIVER-tissue diagnosis/fibrosis-stage by a
    # POSITIONAL index (sample_number-1) with NO shared subject key. There is no
    # obtainable per-subject Olink crosswalk (Mendeley deposit = NPX only;
    # Cell Rep Med supplementary has no per-subject metadata table), and the two
    # cohorts differ in composition. Every supervised Olink result is therefore
    # WITHDRAWN. See docs/reviews/2026-06-01-proteomics-extreme-review.md
    # (findings G1-001/G3-*/G4-001/G5-001/G7-002) and
    # memory/proteomics_extreme_review_2026_06_01.md.
    # Set OLINK_LABELS_OK=1 ONLY after a real subject crosswalk
    # (Analysis/Proteomics/data/olink_plasma/subject_metadata.csv) is obtained
    # and this function is rewritten to key-join on it.
    if os.environ.get("OLINK_LABELS_OK") != "1":
        raise RuntimeError(
            "load_olink_data() is QUARANTINED (2026-06-01): Olink plasma labels "
            "were assigned by an unrecoverable positional join to GSE276114 liver "
            "metadata. No real subject crosswalk exists, so all supervised Olink "
            "results are withdrawn. See docs/reviews/2026-06-01-proteomics-extreme-review.md."
        )
    log.info(f"Loading Olink data from {OLINK_PATH}")
    olink = pd.read_csv(OLINK_PATH, sep="\t")
    protein_names = list(olink["Assay"].values)

    # Subjects × proteins matrix (rows = "Subject N" columns in file)
    npx = olink.drop("Assay", axis=1).T.values.astype(np.float32)

    # NaN values are preserved — imputation must happen per CV fold
    # via impute_train_fold() to avoid data leakage.
    n_nan = int(np.isnan(npx).sum())
    log.info(f"Olink matrix: {npx.shape} ({n_nan} NaN values retained for per-fold imputation)")

    log.info(f"Loading Olink metadata from {OLINK_META_PATH}")
    meta_raw = pd.read_csv(OLINK_META_PATH)
    meta_raw["subject_idx"] = meta_raw["sample_number"].astype(int) - 1  # 0-indexed

    # Keep only subjects within matrix bounds
    meta = meta_raw[meta_raw["subject_idx"] < npx.shape[0]].copy().reset_index(drop=True)

    X = npx[meta["subject_idx"].values]  # (177, 1461)

    # Ensure required columns exist with sensible defaults
    for col in ("geo_accession", "disease_group", "disease", "sample_number"):
        if col not in meta.columns:
            meta[col] = np.nan

    # Build target dict
    dg = meta["disease_group"].values
    disease = meta["disease"].values

    binary = (np.isin(dg, ["F3", "F4"])).astype(int)

    ordinal = np.zeros(len(dg), dtype=int)
    ordinal[dg == "F3"] = 1
    ordinal[dg == "F4"] = 2

    continuous = ordinal.astype(float)

    # Etiology 3-class: MASLD=0, CVH=1, ARLD=2 (fallback to 0 for unknowns)
    etiology_map = {"MASLD": 0, "CVH": 1, "ARLD": 2}
    etiology_3 = np.array(
        [etiology_map.get(str(d).upper(), 0) for d in disease], dtype=int
    )
    etiology_binary = (
        np.array([str(d).upper() for d in disease]) == "MASLD"
    ).astype(int)
    etiologies = np.array([str(d) for d in disease])

    y_targets = {
        "binary": binary,
        "ordinal": ordinal,
        "continuous": continuous,
        "etiology_3": etiology_3,
        "etiology_binary": etiology_binary,
        "etiologies": etiologies,
    }

    n_adv = int(binary.sum())
    log.info(
        f"Olink: {X.shape[0]} subjects × {X.shape[1]} proteins | "
        f"{n_adv} advanced (F3/F4), {X.shape[0] - n_adv} early"
    )
    return X, protein_names, meta, y_targets


def load_ms_data():
    """Load DIA-MS PXD052937 proteomics data.

    Returns
    -------
    X_ms : ndarray (n_samples, n_proteins) float32 — raw values, may contain NaN
        (use impute_train_fold() inside CV loops to avoid leakage)
    protein_names_ms : list of protein name strings
    ms_meta : DataFrame aligned row-for-row to ms_df.index (key-joined on sample_id),
              with columns: sample_id, condition, run_label, condition_letter,
              condition_label, condition_detail, is_masld. The A/B/C/D -> diagnosis
              label mapping is TENTATIVE (carried over from the PXD052937 design).
    """
    log.info(f"Loading DIA-MS matrix from {MS_MATRIX_PATH}")
    ms_df = pd.read_csv(MS_MATRIX_PATH, index_col=0)
    protein_names_ms = list(ms_df.columns)
    X_ms = ms_df.values.astype(np.float32)

    # NaN values are preserved — imputation must happen per CV fold
    # via impute_train_fold() to avoid data leakage.
    n_nan = int(np.isnan(X_ms).sum())
    log.info(f"DIA-MS matrix: {X_ms.shape} ({n_nan} NaN values retained for per-fold imputation)")

    log.info(f"Loading DIA-MS metadata from {MS_META_PATH}")
    meta_raw = pd.read_csv(MS_META_PATH)

    # G1-002 fix (2026-06-01): KEY-JOIN metadata to the matrix index on sample_id.
    # The companion metadata shares the matrix's P-series sample_id, so we index it
    # by sample_id and reindex to ms_df.index — NO positional join. A hard assertion
    # guarantees every matrix sample carried a metadata row (no silent NaN labels).
    assert "sample_id" in meta_raw.columns, (
        "ms_plasma_metadata.csv is missing the 'sample_id' key column required for "
        "the matrix join (G1-002)."
    )
    meta_indexed = meta_raw.set_index("sample_id")
    missing = ms_df.index.difference(meta_indexed.index)
    assert len(missing) == 0, (
        f"G1-002: {len(missing)} DIA-MS matrix samples have no metadata row "
        f"(key-join on sample_id failed). First unmatched: {list(missing[:5])}. "
        "MS_META_PATH must point at the companion export sharing the matrix sample_id."
    )
    ms_meta = meta_indexed.reindex(ms_df.index).reset_index()

    # A/B/C/D -> diagnosis label. TENTATIVE: the letter->diagnosis correspondence is
    # carried over from the original PXD052937 Spectronaut design (A=Normal, B=MASL,
    # C=MASH, D=Cirrhosis); per-letter counts (7/17/38/10) match the prior export, but
    # the underlying clinical key for these P-series runs is not independently confirmed.
    _letter_to_label = {"A": "Normal", "B": "MASL", "C": "MASH", "D": "Cirrhosis"}
    _letter_to_detail = {
        "A": "Healthy control",
        "B": "Simple steatosis (MASL/NAFL)",
        "C": "Steatohepatitis (MASH/NASH)",
        "D": "Advanced fibrosis / Cirrhosis",
    }
    _letter = ms_meta["condition"].astype(str).str.strip()
    ms_meta["condition_letter"] = _letter
    ms_meta["condition_label"] = _letter.map(_letter_to_label)
    ms_meta["condition_detail"] = _letter.map(_letter_to_detail)
    ms_meta["is_masld"] = _letter.isin({"B", "C", "D"})

    n_unmapped = int(ms_meta["condition_label"].isna().sum())
    if n_unmapped:
        log.warning(
            f"DIA-MS: {n_unmapped} samples have a condition letter outside A/B/C/D "
            "(condition_label left as NaN)."
        )

    n_gene_symbols = sum(1 for p in protein_names_ms if not p.startswith("sp|") and "_" not in p[:6])
    log.info(
        f"DIA-MS: {X_ms.shape[0]} samples × {X_ms.shape[1]} proteins "
        f"(~{n_gene_symbols} gene symbols, ~{X_ms.shape[1] - n_gene_symbols} UniProt IDs). "
        f"Low Olink overlap (~48 proteins) reflects DIA-MS enrichment for "
        f"immunoglobulin/complement proteins not in the Olink panel."
    )
    return X_ms, protein_names_ms, ms_meta


# ---------------------------------------------------------------------------
# OrdinalClassifierWrapper
# ---------------------------------------------------------------------------


class OrdinalClassifierWrapper(BaseEstimator, ClassifierMixin):
    """Cumulative logit ordinal regression via K-1 binary classifiers.

    For an ordinal outcome with classes 0, 1, …, K-1, fits K-1 binary
    classifiers: classifier k predicts P(Y >= k) for k in {1, …, K-1}.
    Prediction is the class c that maximises the estimated probability, where
    class probabilities are derived from adjacent cumulative probabilities.

    Parameters
    ----------
    base_estimator : sklearn estimator
        Binary classifier to clone for each threshold.
    n_classes : int, default 3
        Number of ordinal classes (K). Fits K-1 binary classifiers.
    """

    def __init__(self, base_estimator=None, n_classes=3):
        self.base_estimator = base_estimator
        self.n_classes = n_classes

    def fit(self, X, y):
        y = np.asarray(y, dtype=int)
        self.classes_ = np.arange(self.n_classes)
        self._clfs = []
        for k in range(1, self.n_classes):
            y_bin = (y >= k).astype(int)
            clf = clone(self.base_estimator)
            clf.fit(X, y_bin)
            self._clfs.append(clf)
        return self

    def _cumulative_probs(self, X):
        """Return array (n, K-1) of P(Y >= k) for k = 1, …, K-1."""
        cum = np.column_stack(
            [
                clf.predict_proba(X)[:, 1] if hasattr(clf, "predict_proba")
                else _decision_to_prob(clf.decision_function(X))
                for clf in self._clfs
            ]
        )
        return cum

    def predict_proba(self, X):
        """Return class probability matrix (n, K).

        P(Y = 0) = 1 - P(Y >= 1)
        P(Y = k) = P(Y >= k) - P(Y >= k+1)  for 0 < k < K-1
        P(Y = K-1) = P(Y >= K-1)
        """
        cum = self._cumulative_probs(X)
        n = X.shape[0]
        K = self.n_classes
        proba = np.zeros((n, K), dtype=float)
        proba[:, 0] = 1.0 - cum[:, 0]
        for k in range(1, K - 1):
            proba[:, k] = cum[:, k - 1] - cum[:, k]
        proba[:, K - 1] = cum[:, K - 2]
        # Clip to [0, 1] to absorb floating-point drift
        proba = np.clip(proba, 0.0, 1.0)
        # Renormalise rows
        row_sums = proba.sum(axis=1, keepdims=True)
        proba /= np.where(row_sums == 0, 1.0, row_sums)
        return proba

    def predict(self, X):
        return np.argmax(self.predict_proba(X), axis=1)


def _decision_to_prob(scores):
    """Convert raw decision scores to [0, 1] via sigmoid."""
    return 1.0 / (1.0 + np.exp(-scores))


# ---------------------------------------------------------------------------
# Helper wrappers for non-standard estimators
# ---------------------------------------------------------------------------


class _BayesianRidgeClassifier(BaseEstimator, ClassifierMixin):
    """BayesianRidge wrapped as binary classifier."""

    def __init__(self):
        self.model_ = BayesianRidge()
        self.classes_ = np.array([0, 1])

    def fit(self, X, y):
        self.model_.fit(X, y.astype(float))
        return self

    def predict(self, X):
        return (self.model_.predict(X) >= 0.5).astype(int)

    def predict_proba(self, X):
        p = np.clip(self.model_.predict(X), 0.0, 1.0)
        return np.column_stack([1 - p, p])

    def get_params(self, deep=True):
        return {}

    def set_params(self, **params):
        return self


class _TabNetWrapper(BaseEstimator, ClassifierMixin):
    """Sklearn-compatible wrapper around pytorch-tabnet TabNetClassifier.

    NOT Pipeline-compatible. Handles its own train/validation split.
    """

    def __init__(self, n_steps=3, n_a=8, n_d=8, random_state=SEED):
        self.n_steps = n_steps
        self.n_a = n_a
        self.n_d = n_d
        self.random_state = random_state

    def fit(self, X, y):
        self._model = _TabNet(
            n_steps=self.n_steps,
            n_a=self.n_a,
            n_d=self.n_d,
            seed=self.random_state,
            verbose=0,
            device_name="cpu",
        )
        n_val = max(int(0.15 * len(y)), 5)
        rng = np.random.RandomState(self.random_state)
        idx = np.arange(len(y))
        rng.shuffle(idx)
        val_idx, train_idx = idx[:n_val], idx[n_val:]
        self._model.fit(
            X[train_idx].astype(np.float32),
            y[train_idx],
            eval_set=[(X[val_idx].astype(np.float32), y[val_idx])],
            max_epochs=100,
            patience=15,
            batch_size=64,
        )
        self.classes_ = np.array([0, 1])
        return self

    def predict(self, X):
        return self._model.predict(X.astype(np.float32))

    def predict_proba(self, X):
        return self._model.predict_proba(X.astype(np.float32))

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


class _FLAMLWrapper(BaseEstimator, ClassifierMixin):
    """Sklearn-compatible wrapper around FLAML AutoML.

    NOT Pipeline-compatible.
    """

    def __init__(self, time_budget=30, random_state=SEED):
        self.time_budget = time_budget
        self.random_state = random_state

    def fit(self, X, y):
        self._automl = _FLAMLAutoML()
        self._automl.fit(
            X,
            y,
            task="classification",
            metric="roc_auc",
            time_budget=self.time_budget,
            verbose=0,
            seed=self.random_state,
        )
        self.classes_ = np.array([0, 1])
        return self

    def predict(self, X):
        return self._automl.predict(X)

    def predict_proba(self, X):
        return self._automl.predict_proba(X)

    def get_params(self, deep=True):
        return {"time_budget": self.time_budget, "random_state": self.random_state}

    def set_params(self, **params):
        for k, v in params.items():
            setattr(self, k, v)
        return self


# ---------------------------------------------------------------------------
# Scoring helpers
# ---------------------------------------------------------------------------


def _qwk_scorer(y_true, y_pred):
    return cohen_kappa_score(y_true, y_pred, weights="quadratic")


_qwk_make = make_scorer(_qwk_scorer)


def _scoring_for_target(target_type):
    return {
        "binary": "roc_auc",
        "ordinal": _qwk_make,
        "regression": "neg_mean_absolute_error",
        "multiclass": "f1_macro",
        "multitask": "f1_macro",  # fallback composite
    }[target_type]


# ---------------------------------------------------------------------------
# Model registry
# ---------------------------------------------------------------------------


def get_model_registry(target_type: str) -> dict:
    """Return dict of model_name → {"model": estimator, "params": grid, "is_pipeline": bool}.

    Parameters
    ----------
    target_type : str
        One of "binary", "ordinal", "regression", "multiclass", "multitask".
    """
    assert target_type in (
        "binary", "ordinal", "regression", "multiclass", "multitask"
    ), f"Unknown target_type: {target_type!r}"

    is_reg = target_type == "regression"
    is_ordinal = target_type == "ordinal"
    is_multi = target_type in ("multiclass", "multitask")

    # ------------------------------------------------------------------
    # Helper to wrap a classifier for ordinal/multiclass/multitask
    # ------------------------------------------------------------------
    def _remap_params(params, prefix):
        """Remap 'model__param' to 'model__prefix__param' for wrapped estimators.

        GridSearchCV cannot resolve params like ``model__max_depth`` when the
        Pipeline step ``model`` is itself a wrapper (OrdinalClassifierWrapper,
        OneVsRestClassifier, MultiOutputClassifier).  The inner estimator is
        exposed via ``base_estimator`` or ``estimator``, so params must be
        prefixed accordingly.
        """
        remapped = {}
        for key, val in params.items():
            if key.startswith("model__"):
                suffix = key[len("model__"):]
                remapped[f"model__{prefix}__{suffix}"] = val
            else:
                remapped[key] = val
        return remapped

    def _wrap_clf(clf):
        """Wrap classifier for ordinal/multiclass/multitask. Returns (wrapped, True)
        or (clf, False) so callers can track whether remapping is needed."""
        if is_ordinal:
            return OrdinalClassifierWrapper(base_estimator=clf, n_classes=3)
        if target_type == "multitask":
            return MultiOutputClassifier(clf)
        if target_type == "multiclass":
            # Most sklearn classifiers handle multiclass natively, but OVR is explicit
            return OneVsRestClassifier(clf)
        return clf

    def _is_wrapped_model(model, is_ordinal, target_type):
        """Check if a model object is wrapped by OrdinalClassifierWrapper/OneVsRestClassifier/MultiOutputClassifier."""
        if is_ordinal:
            return isinstance(model, OrdinalClassifierWrapper)
        if target_type == "multiclass":
            return isinstance(model, OneVsRestClassifier)
        if target_type == "multitask":
            return isinstance(model, MultiOutputClassifier)
        return False

    def _pfx(name, is_pipeline=True):
        """Return param key prefix (empty for non-pipeline, 'model__' for pipeline)."""
        return "model__" if is_pipeline else ""

    # Prefix helper for pipeline param grids
    P = "model__"

    # Class-balance weight string (not applicable to regression)
    CW = None if is_reg else "balanced"

    # ------------------------------------------------------------------
    # Scale pos weight for XGBoost (approximate; recomputed at fit time
    # in run_nested_cv from training fold data)
    # ------------------------------------------------------------------
    _spw = 3.5  # approximate F0-2 : F3/F4 ratio placeholder

    registry = {}

    # ---- 1. ElasticNet / LogisticRegression elastic ----
    if is_reg:
        registry["elastic_net"] = {
            "model": Pipeline([("scaler", StandardScaler()), ("model", ElasticNet(max_iter=5000, random_state=SEED))]),
            "params": {f"{P}alpha": [0.001, 0.01, 0.1, 1, 10], f"{P}l1_ratio": [0.1, 0.5, 0.9]},
            "is_pipeline": True,
        }
    else:
        _en_clf = LogisticRegression(
            solver="saga", max_iter=5000, class_weight=CW, random_state=SEED, penalty="elasticnet", l1_ratio=0.5
        )
        registry["elastic_net"] = {
            "model": _wrap_clf(_en_clf),
            "params": {f"{P}C": [0.001, 0.01, 0.1, 1, 10, 100], f"{P}l1_ratio": [0.1, 0.3, 0.5, 0.7, 0.9]},
            "is_pipeline": True,
        }

    # ---- 2. Ridge ----
    if is_reg:
        registry["ridge"] = {
            "model": Pipeline([("scaler", StandardScaler()), ("model", Ridge())]),
            "params": {f"{P}alpha": [0.001, 0.01, 0.1, 1, 10, 100, 1000]},
            "is_pipeline": True,
        }
    else:
        registry["ridge"] = {
            "model": _wrap_clf(RidgeClassifier(class_weight=CW)),
            "params": {f"{P}alpha": [0.001, 0.01, 0.1, 1, 10, 100, 1000]},
            "is_pipeline": True,
        }

    # ---- 3. Lasso ----
    if is_reg:
        registry["lasso"] = {
            "model": Pipeline([("scaler", StandardScaler()), ("model", Lasso(max_iter=5000, random_state=SEED))]),
            "params": {f"{P}alpha": [0.001, 0.01, 0.1, 1, 10]},
            "is_pipeline": True,
        }
    else:
        _lasso_clf = LogisticRegression(penalty="l1", solver="liblinear", max_iter=5000, class_weight=CW, random_state=SEED)
        registry["lasso"] = {
            "model": _wrap_clf(_lasso_clf),
            "params": {f"{P}C": [0.001, 0.01, 0.1, 1, 10, 100]},
            "is_pipeline": True,
        }

    # ---- 4. SGD ----
    if is_reg:
        from sklearn.linear_model import SGDRegressor
        registry["sgd_elasticnet"] = {
            "model": Pipeline([("scaler", StandardScaler()), ("model", SGDRegressor(penalty="elasticnet", max_iter=2000, random_state=SEED))]),
            "params": {f"{P}alpha": [1e-5, 1e-4, 1e-3, 1e-2, 0.1], f"{P}l1_ratio": [0.1, 0.5, 0.9]},
            "is_pipeline": True,
        }
    else:
        registry["sgd_elasticnet"] = {
            "model": _wrap_clf(SGDClassifier(loss="log_loss", penalty="elasticnet", class_weight=CW, random_state=SEED, max_iter=2000)),
            "params": {f"{P}alpha": [1e-5, 1e-4, 1e-3, 1e-2, 0.1], f"{P}l1_ratio": [0.1, 0.3, 0.5, 0.7, 0.9]},
            "is_pipeline": True,
        }

    # ---- 5. BayesianRidge ----
    if is_reg:
        registry["bayesian_ridge"] = {
            "model": Pipeline([("scaler", StandardScaler()), ("model", BayesianRidge())]),
            "params": {},
            "is_pipeline": True,
        }
    else:
        registry["bayesian_ridge"] = {
            "model": _wrap_clf(_BayesianRidgeClassifier()),
            "params": {},
            "is_pipeline": True,
        }

    # ---- 6. RandomForest ----
    _rf = RandomForestClassifier(class_weight=CW, random_state=SEED) if not is_reg else RandomForestClassifier(random_state=SEED)
    if is_reg:
        from sklearn.ensemble import RandomForestRegressor
        registry["random_forest"] = {
            "model": Pipeline([("scaler", StandardScaler()), ("model", RandomForestRegressor(random_state=SEED))]),
            "params": {f"{P}n_estimators": [100, 500], f"{P}max_depth": [5, 10, None], f"{P}min_samples_leaf": [5, 10]},
            "is_pipeline": True,
        }
    else:
        registry["random_forest"] = {
            "model": _wrap_clf(_rf),
            "params": {f"{P}n_estimators": [100, 500, 1000], f"{P}max_depth": [5, 10, None], f"{P}min_samples_leaf": [5, 10, 20]},
            "is_pipeline": True,
        }

    # ---- 7. BalancedRF (imblearn) ----
    if HAS_IMBLEARN and not is_reg:
        registry["balanced_rf"] = {
            "model": _wrap_clf(BalancedRandomForestClassifier(random_state=SEED)),
            "params": {f"{P}n_estimators": [100, 500], f"{P}max_depth": [5, 10, None], f"{P}min_samples_leaf": [5, 10]},
            "is_pipeline": True,
        }

    # ---- 8. ExtraTrees ----
    if is_reg:
        from sklearn.ensemble import ExtraTreesRegressor
        registry["extra_trees"] = {
            "model": Pipeline([("scaler", StandardScaler()), ("model", ExtraTreesRegressor(random_state=SEED))]),
            "params": {f"{P}n_estimators": [100, 500], f"{P}max_depth": [5, 10, None], f"{P}min_samples_leaf": [5, 10]},
            "is_pipeline": True,
        }
    else:
        registry["extra_trees"] = {
            "model": _wrap_clf(ExtraTreesClassifier(class_weight=CW, random_state=SEED)),
            "params": {f"{P}n_estimators": [100, 500], f"{P}max_depth": [5, 10, None], f"{P}min_samples_leaf": [5, 10]},
            "is_pipeline": True,
        }

    # ---- 9. GradientBoosting ----
    registry["gradient_boosting"] = {
        "model": (
            Pipeline([("scaler", StandardScaler()), ("model", GradientBoostingClassifier(random_state=SEED))])
            if not is_reg
            else Pipeline([("scaler", StandardScaler()), ("model", GradientBoostingClassifier(random_state=SEED))])
        ),
        "params": {
            f"{P}n_estimators": [100, 300],
            f"{P}learning_rate": [0.01, 0.1],
            f"{P}max_depth": [3, 5],
            f"{P}subsample": [0.8],
            f"{P}min_samples_leaf": [10],
        },
        "is_pipeline": True,
    }
    if is_reg:
        from sklearn.ensemble import GradientBoostingRegressor
        registry["gradient_boosting"]["model"] = Pipeline([("scaler", StandardScaler()), ("model", GradientBoostingRegressor(random_state=SEED))])

    # ---- 10. HistGBM ----
    if is_reg:
        from sklearn.ensemble import HistGradientBoostingRegressor
        registry["hist_gradient_boosting"] = {
            "model": Pipeline([("scaler", StandardScaler()), ("model", HistGradientBoostingRegressor(random_state=SEED))]),
            "params": {f"{P}max_iter": [100, 300], f"{P}max_depth": [3, 5, 7], f"{P}learning_rate": [0.01, 0.1]},
            "is_pipeline": True,
        }
    else:
        registry["hist_gradient_boosting"] = {
            "model": _wrap_clf(HistGradientBoostingClassifier(class_weight=CW, random_state=SEED, early_stopping=True, validation_fraction=0.15)),
            "params": {f"{P}max_iter": [100, 300], f"{P}max_depth": [3, 5, 7], f"{P}learning_rate": [0.01, 0.1], f"{P}l2_regularization": [0, 1, 10]},
            "is_pipeline": True,
        }

    # ---- 11. XGBoost ----
    if HAS_XGB:
        if is_reg:
            registry["xgboost"] = {
                "model": Pipeline([("scaler", StandardScaler()), ("model", xgb.XGBRegressor(random_state=SEED, eval_metric="rmse"))]),
                "params": {f"{P}learning_rate": [0.01, 0.05, 0.1], f"{P}max_depth": [3, 5], f"{P}n_estimators": [100, 300], f"{P}subsample": [0.8]},
                "is_pipeline": True,
            }
        else:
            registry["xgboost"] = {
                "model": _wrap_clf(xgb.XGBClassifier(scale_pos_weight=_spw, random_state=SEED, eval_metric="logloss", use_label_encoder=False)),
                "params": {
                    f"{P}learning_rate": [0.01, 0.05, 0.1],
                    f"{P}max_depth": [3, 5, 7],
                    f"{P}n_estimators": [100, 300],
                    f"{P}subsample": [0.6, 0.8],
                    f"{P}colsample_bytree": [0.5, 0.8],
                    f"{P}reg_lambda": [1, 10],
                },
                "is_pipeline": True,
            }

    # ---- 12. LightGBM ----
    if HAS_LGB:
        if is_reg:
            registry["lightgbm"] = {
                "model": Pipeline([("scaler", StandardScaler()), ("model", lgb.LGBMRegressor(random_state=SEED, verbose=-1, n_jobs=1))]),
                "params": {f"{P}num_leaves": [15, 31], f"{P}learning_rate": [0.01, 0.1], f"{P}n_estimators": [100, 300]},
                "is_pipeline": True,
                "serial_grid": True,
            }
        else:
            registry["lightgbm"] = {
                "model": _wrap_clf(lgb.LGBMClassifier(is_unbalance=True, random_state=SEED, verbose=-1, n_jobs=1)),
                "params": {
                    f"{P}num_leaves": [15, 31],
                    f"{P}learning_rate": [0.01, 0.1],
                    f"{P}n_estimators": [100, 300],
                    f"{P}reg_lambda": [1, 10],
                    f"{P}feature_fraction": [0.5, 0.8],
                },
                "is_pipeline": True,
                "serial_grid": True,  # LightGBM deadlocks with joblib parallel GridSearchCV
            }

    # ---- 13. LinearSVC ----
    if not is_reg:
        registry["linear_svc"] = {
            "model": _wrap_clf(SVC(kernel="linear", class_weight=CW, probability=True, random_state=SEED)),
            "params": {f"{P}C": [0.01, 0.1, 1, 10, 100]},
            "is_pipeline": True,
        }

    # ---- 14. RBF-SVM ----
    if is_reg:
        registry["rbf_svm"] = {
            "model": Pipeline([("scaler", StandardScaler()), ("model", SVR(kernel="rbf"))]),
            "params": {f"{P}C": [0.1, 1, 10, 100], f"{P}gamma": ["scale", "auto"]},
            "is_pipeline": True,
        }
    else:
        registry["rbf_svm"] = {
            "model": _wrap_clf(SVC(kernel="rbf", class_weight=CW, probability=True, random_state=SEED)),
            "params": {f"{P}C": [0.1, 1, 10, 100], f"{P}gamma": ["scale", "auto", 0.001, 0.01]},
            "is_pipeline": True,
        }

    # ---- 15. NuSVC ----
    if not is_reg:
        registry["nu_svc"] = {
            "model": _wrap_clf(NuSVC(kernel="rbf", class_weight=CW, probability=True, random_state=SEED)),
            "params": {f"{P}nu": [0.1, 0.3, 0.5], f"{P}gamma": ["scale", "auto"]},
            "is_pipeline": True,
        }

    # ---- 16. MLP ----
    if not is_reg:
        registry["mlp"] = {
            "model": _wrap_clf(MLPClassifier(random_state=SEED, max_iter=500, early_stopping=True, validation_fraction=0.15)),
            "params": {f"{P}hidden_layer_sizes": [(64,), (128, 64), (256, 128)], f"{P}alpha": [0.001, 0.01, 0.1]},
            "is_pipeline": True,
        }
    else:
        from sklearn.neural_network import MLPRegressor
        registry["mlp"] = {
            "model": Pipeline([("scaler", StandardScaler()), ("model", MLPRegressor(random_state=SEED, max_iter=500, early_stopping=True))]),
            "params": {f"{P}hidden_layer_sizes": [(64,), (128, 64)], f"{P}alpha": [0.001, 0.01]},
            "is_pipeline": True,
        }

    # ---- 17. TabNet ----
    if HAS_TABNET and not is_reg:
        registry["tabnet"] = {
            "model": _TabNetWrapper(n_steps=3, n_a=8, n_d=8, random_state=SEED),
            "params": {},
            "is_pipeline": False,
        }

    # ---- 18. TabPFN ----
    if HAS_TABPFN and not is_reg:
        registry["tabpfn"] = {
            "model": TabPFNClassifier(device="cpu", n_estimators=32),
            "params": {},
            "is_pipeline": False,  # needs top-100 features, handled in run_nested_cv
            "top100": True,
        }

    # ---- 19. EasyEnsemble ----
    if HAS_IMBLEARN and not is_reg:
        registry["easy_ensemble"] = {
            "model": _wrap_clf(EasyEnsembleClassifier(random_state=SEED)) if not is_ordinal else EasyEnsembleClassifier(random_state=SEED),
            "params": {f"{P}n_estimators": [10, 20, 50]},
            "is_pipeline": True,
        }

    # ---- 20. BalancedBagging ----
    if HAS_IMBLEARN and not is_reg:
        _bb_base = LogisticRegression(max_iter=2000, random_state=SEED)
        registry["balanced_bagging"] = {
            "model": _wrap_clf(BalancedBaggingClassifier(estimator=_bb_base, random_state=SEED)),
            "params": {f"{P}n_estimators": [20, 50, 100]},
            "is_pipeline": True,
        }

    # ---- 21. AdaBoost ----
    if not is_reg:
        registry["adaboost"] = {
            "model": _wrap_clf(AdaBoostClassifier(random_state=SEED)),
            "params": {f"{P}n_estimators": [50, 100, 200], f"{P}learning_rate": [0.1, 0.5, 1.0]},
            "is_pipeline": True,
        }
    else:
        from sklearn.ensemble import AdaBoostRegressor
        registry["adaboost"] = {
            "model": Pipeline([("scaler", StandardScaler()), ("model", AdaBoostRegressor(random_state=SEED))]),
            "params": {f"{P}n_estimators": [50, 100, 200], f"{P}learning_rate": [0.1, 0.5, 1.0]},
            "is_pipeline": True,
        }

    # ---- 22-23. FLAML ----
    if HAS_FLAML and not is_reg:
        registry["flaml_30s"] = {"model": _FLAMLWrapper(time_budget=30, random_state=SEED), "params": {}, "is_pipeline": False}
        registry["flaml_60s"] = {"model": _FLAMLWrapper(time_budget=60, random_state=SEED), "params": {}, "is_pipeline": False}

    # ---- 24. VotingClassifier ----
    if not is_reg:
        _v_lr = LogisticRegression(penalty="l2", C=1.0, class_weight=CW, max_iter=5000, random_state=SEED)
        _v_rf = RandomForestClassifier(n_estimators=500, class_weight=CW, random_state=SEED)
        _v_hgb = HistGradientBoostingClassifier(class_weight=CW, random_state=SEED)
        registry["voting_lr_rf_hgb"] = {
            "model": VotingClassifier(estimators=[("lr", _v_lr), ("rf", _v_rf), ("hgb", _v_hgb)], voting="soft"),
            "params": {},
            "is_pipeline": True,
        }

    # ---- 25. StackingClassifier ----
    if not is_reg:
        registry["stacking_lr_rf_hgb"] = {
            "model": StackingClassifier(
                estimators=[
                    ("lr", LogisticRegression(penalty="l2", C=1.0, class_weight=CW, max_iter=5000, random_state=SEED)),
                    ("rf", RandomForestClassifier(n_estimators=500, class_weight=CW, random_state=SEED)),
                    ("hgb", HistGradientBoostingClassifier(class_weight=CW, random_state=SEED)),
                ],
                final_estimator=LogisticRegression(class_weight=CW, max_iter=5000, random_state=SEED),
                cv=3,
            ),
            "params": {},
            "is_pipeline": True,
        }

    # ---- 26. CatBoost ----
    if HAS_CATBOOST:
        if is_reg:
            registry["catboost"] = {
                "model": Pipeline([("scaler", StandardScaler()), ("model", CatBoostRegressor(random_seed=SEED, verbose=0, allow_writing_files=False))]),
                "params": {f"{P}iterations": [100, 300], f"{P}learning_rate": [0.01, 0.1], f"{P}depth": [4, 6]},
                "is_pipeline": True,
            }
        else:
            registry["catboost"] = {
                "model": _wrap_clf(CatBoostClassifier(random_seed=SEED, verbose=0, auto_class_weights="Balanced", allow_writing_files=False)),
                "params": {f"{P}iterations": [100, 300], f"{P}learning_rate": [0.01, 0.1], f"{P}depth": [4, 6]},
                "is_pipeline": True,
            }

    # ---- 27. OrdinalForest (mord) ----
    if HAS_MORD and is_ordinal:
        registry["ordinal_forest"] = {
            "model": mord.OrdinalRidge(),
            "params": {f"{P}alpha": [0.001, 0.01, 0.1, 1, 10]},
            "is_pipeline": True,  # mord models are sklearn-compatible; wrapped in Pipeline by run_nested_cv
        }

    # ---- 28. CORAL (mord) ----
    if HAS_MORD and is_ordinal:
        registry["coral"] = {
            "model": mord.LogisticAT(),
            "params": {f"{P}alpha": [0.001, 0.01, 0.1, 1, 10]},
            "is_pipeline": True,
        }

    # ---- 29. Platt-SVM (calibrated) ----
    if not is_reg:
        registry["platt_svm"] = {
            "model": _wrap_clf(CalibratedClassifierCV(LinearSVC(class_weight=CW, max_iter=5000, random_state=SEED), method="sigmoid", cv=3)),
            # Params already have nested prefix — mark skip_remap to avoid double-prefixing
            "params": {f"{P}estimator__base_estimator__C": [0.01, 0.1, 1, 10]},
            "is_pipeline": True,
            "skip_remap": True,
        }

    # ---- 30. GaussianProcess ----
    if not is_reg:
        # Subset to top-100 variance features (handled in run_nested_cv via "top100" flag)
        registry["gaussian_process"] = {
            "model": GaussianProcessClassifier(kernel=1.0 * RBF(1.0), random_state=SEED, n_jobs=1),
            "params": {},
            "is_pipeline": False,
            "top100": True,
        }

    # ---- 31. KernelRidge ----
    registry["kernel_ridge"] = {
        "model": Pipeline([("scaler", StandardScaler()), ("model", KernelRidge(kernel="rbf"))]),
        "params": {f"{P}alpha": [0.001, 0.01, 0.1, 1], f"{P}gamma": [0.001, 0.01, 0.1]},
        "is_pipeline": True,
    }

    # ---- 32. ElasticNetReg (explicit regression variant even for classification) ----
    registry["elasticnet_reg"] = {
        "model": Pipeline([("scaler", StandardScaler()), ("model", ElasticNet(max_iter=5000, random_state=SEED))]),
        "params": {f"{P}alpha": [0.001, 0.01, 0.1, 1, 10], f"{P}l1_ratio": [0.1, 0.5, 0.9]},
        "is_pipeline": True,
    }

    # ---- 33. MultiOutput ----
    if target_type == "multitask":
        registry["multioutput_rf"] = {
            "model": MultiOutputClassifier(RandomForestClassifier(class_weight=CW, random_state=SEED)),
            "params": {f"estimator__n_estimators": [100, 500], f"estimator__max_depth": [5, 10, None]},
            "is_pipeline": True,
        }

    # ---- 34. CalibratedEnsemble ----
    if not is_reg:
        _cal_base = RandomForestClassifier(n_estimators=500, class_weight=CW, random_state=SEED)
        registry["calibrated_ensemble"] = {
            "model": _wrap_clf(CalibratedClassifierCV(_cal_base, method="isotonic", cv=3)),
            "params": {},
            "is_pipeline": True,
        }

    # ---- 35. SMOTE + best (imblearn) ----
    if HAS_IMBLEARN and not is_reg:
        from imblearn.pipeline import Pipeline as ImbPipeline
        registry["smote_rf"] = {
            "model": ImbPipeline([
                ("scaler", StandardScaler()),
                ("smote", SMOTE(random_state=SEED)),
                ("model", RandomForestClassifier(random_state=SEED)),
            ]),
            "params": {f"model__n_estimators": [100, 500], f"model__max_depth": [5, 10, None]},
            "is_pipeline": True,
        }

    # ------------------------------------------------------------------
    # Post-process: remap param grid keys for wrapped estimators so that
    # GridSearchCV can route params through the wrapper's estimator attribute.
    #   OrdinalClassifierWrapper  → base_estimator__<param>
    #   OneVsRestClassifier       → estimator__<param>
    #   MultiOutputClassifier     → estimator__<param>
    # This must happen after the full registry is built because _wrap_clf is
    # called inline during entry construction.
    # ------------------------------------------------------------------
    if is_ordinal:
        _wrapper_prefix = "base_estimator"
    elif target_type in ("multiclass", "multitask"):
        _wrapper_prefix = "estimator"
    else:
        _wrapper_prefix = None

    if _wrapper_prefix is not None:
        # Only remap models that were actually wrapped via _wrap_clf().
        # Models that build their own Pipeline (e.g. gradient_boosting, kernel_ridge,
        # voting, stacking, mord models, ImbPipeline models) already have correct
        # param prefixes and must NOT be remapped.
        _wrapped_models = {
            name for name, cfg in registry.items()
            if _is_wrapped_model(cfg["model"], is_ordinal, target_type)
        }
        for name, cfg in registry.items():
            if name in _wrapped_models and cfg.get("params") and not cfg.get("skip_remap"):
                cfg["params"] = _remap_params(cfg["params"], _wrapper_prefix)

    log.info(f"Model registry for target_type={target_type!r}: {len(registry)} models")
    return registry


# ---------------------------------------------------------------------------
# Nested CV
# ---------------------------------------------------------------------------


def run_nested_cv(
    X: np.ndarray,
    y: np.ndarray,
    model_name: str,
    model_cfg: dict,
    target_type: str,
    n_repeats: int = 10,
    n_splits: int = 5,
) -> dict:
    """Generalised nested CV for all target types.

    Parameters
    ----------
    X : ndarray (n, p)
    y : ndarray (n,) or (n, t) for multitask
    model_name : str
    model_cfg : dict — {"model": estimator, "params": grid, "is_pipeline": bool, ...}
    target_type : one of "binary", "ordinal", "regression", "multiclass", "multitask"
    n_repeats : int
    n_splits : int

    Returns
    -------
    dict with model, target_type, mean/sd per metric, n_folds
    """
    scoring = _scoring_for_target(target_type)
    is_reg = target_type == "regression"
    needs_top100 = model_cfg.get("top100", False)
    is_pipeline_model = model_cfg.get("is_pipeline", True)

    # Choose outer CV
    if is_reg:
        outer_cv = RepeatedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=SEED)
    else:
        # Use first column of y for stratification if multitask
        y_strat = y[:, 0] if y.ndim > 1 else y
        outer_cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=SEED)

    inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=SEED) if not is_reg else KFold(n_splits=3, shuffle=True, random_state=SEED)

    # Metric collectors
    metric_lists = {}  # metric_name → list of fold values

    fold_iter = outer_cv.split(X, y_strat if not is_reg else y)

    for fold_idx, (train_idx, test_idx) in enumerate(fold_iter):
        X_train_raw, X_test_raw = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Per-fold imputation: statistics from training data only
        X_train, X_test = impute_train_fold(X_train_raw, X_test_raw)

        # Feature subset for memory-intensive models
        if needs_top100:
            var = np.var(X_train, axis=0)
            top_idx = np.argsort(var)[-100:]
            X_tr_use = X_train[:, top_idx]
            X_te_use = X_test[:, top_idx]
        else:
            X_tr_use = X_train
            X_te_use = X_test

        estimator = clone(model_cfg["model"])
        param_grid = model_cfg.get("params", {})

        if is_pipeline_model:
            # Wrap in scaler Pipeline if model is a raw estimator (not already a Pipeline)
            from sklearn.pipeline import Pipeline as _Pipeline
            if not isinstance(estimator, _Pipeline):
                estimator = _Pipeline([("scaler", StandardScaler()), ("model", estimator)])
                # Prepend "model__" prefix is already in param_grid keys from registry

            if param_grid:
                # LightGBM deadlocks with joblib parallel workers; use n_jobs=1 for it
                _grid_njobs = 1 if model_cfg.get("serial_grid", False) else 4
                grid = GridSearchCV(estimator, param_grid, cv=inner_cv, scoring=scoring, n_jobs=_grid_njobs, refit=True, error_score="raise")
                try:
                    grid.fit(X_tr_use, y_train)
                    best_est = grid.best_estimator_
                except Exception as e:
                    log.warning(f"  [{model_name}] fold {fold_idx} GridSearchCV failed: {e}")
                    continue
            else:
                try:
                    estimator.fit(X_tr_use, y_train)
                    best_est = estimator
                except Exception as e:
                    log.warning(f"  [{model_name}] fold {fold_idx} fit failed: {e}")
                    continue
        else:
            # Non-Pipeline: scale externally
            scaler = StandardScaler()
            X_tr_sc = scaler.fit_transform(X_tr_use)
            X_te_sc = scaler.transform(X_te_use)
            try:
                estimator.fit(X_tr_sc, y_train)
                best_est = estimator
                X_te_use = X_te_sc
            except Exception as e:
                log.warning(f"  [{model_name}] fold {fold_idx} fit failed: {e}")
                continue

        # Predict
        try:
            y_pred = best_est.predict(X_te_use)
            y_prob = None
            if hasattr(best_est, "predict_proba"):
                y_prob = best_est.predict_proba(X_te_use)
            elif hasattr(best_est, "decision_function"):
                df = best_est.decision_function(X_te_use)
                y_prob = _decision_to_prob(df) if df.ndim == 1 else df
        except Exception as e:
            log.warning(f"  [{model_name}] fold {fold_idx} predict failed: {e}")
            continue

        fold_metrics = compute_metrics(y_test, y_pred, y_prob, target_type)
        for k, v in fold_metrics.items():
            metric_lists.setdefault(k, []).append(v)

    if not metric_lists:
        log.warning(f"  [{model_name}] all folds failed — returning NaN results")
        return {"model": model_name, "target_type": target_type, "n_folds": 0}

    result = {"model": model_name, "target_type": target_type, "n_folds": len(next(iter(metric_lists.values())))}
    for metric, vals in metric_lists.items():
        arr = np.array(vals, dtype=float)
        result[f"mean_{metric}"] = float(np.nanmean(arr))
        result[f"sd_{metric}"] = float(np.nanstd(arr))
    return result


# ---------------------------------------------------------------------------
# Compute metrics
# ---------------------------------------------------------------------------


def compute_metrics(
    y_true: np.ndarray,
    y_pred: np.ndarray,
    y_prob,
    target_type: str,
) -> dict:
    """Return metric dict appropriate to target_type.

    Parameters
    ----------
    y_true : ground-truth labels/values
    y_pred : predicted labels/values
    y_prob : probability array (n, K) or None
    target_type : str

    Returns
    -------
    dict of metric_name → float
    """
    metrics = {}

    if target_type == "binary":
        metrics["bal_acc"] = balanced_accuracy_score(y_true, y_pred)
        metrics["f1"] = f1_score(y_true, y_pred, zero_division=0)
        if y_prob is not None:
            prob_pos = y_prob[:, 1] if y_prob.ndim == 2 else y_prob
            try:
                metrics["auroc"] = roc_auc_score(y_true, prob_pos)
            except Exception:
                metrics["auroc"] = float("nan")

    elif target_type == "ordinal":
        try:
            metrics["qwk"] = cohen_kappa_score(y_true, y_pred, weights="quadratic")
        except Exception:
            metrics["qwk"] = float("nan")
        metrics["mae"] = float(mean_absolute_error(y_true, y_pred))
        metrics["bal_acc"] = balanced_accuracy_score(y_true, y_pred)

    elif target_type == "regression":
        metrics["mae"] = float(mean_absolute_error(y_true, y_pred))
        try:
            rho, _ = spearmanr(y_true, y_pred)
            metrics["spearman"] = float(rho)
        except Exception:
            metrics["spearman"] = float("nan")

    elif target_type == "multiclass":
        metrics["f1_macro"] = f1_score(y_true, y_pred, average="macro", zero_division=0)
        metrics["f1_weighted"] = f1_score(y_true, y_pred, average="weighted", zero_division=0)
        metrics["bal_acc"] = balanced_accuracy_score(y_true, y_pred)
        if y_prob is not None and y_prob.ndim == 2:
            try:
                metrics["auroc_ovr"] = roc_auc_score(y_true, y_prob, multi_class="ovr", average="macro")
            except Exception:
                metrics["auroc_ovr"] = float("nan")

    elif target_type == "multitask":
        # y_true and y_pred may be 2-D (n, t)
        y_true = np.asarray(y_true)
        y_pred = np.asarray(y_pred)
        if y_true.ndim == 1:
            y_true = y_true[:, None]
            y_pred = y_pred[:, None]
        per_task_f1 = []
        for t in range(y_true.shape[1]):
            try:
                per_task_f1.append(f1_score(y_true[:, t], y_pred[:, t], average="macro", zero_division=0))
            except Exception:
                per_task_f1.append(float("nan"))
        metrics["f1_macro_mean"] = float(np.nanmean(per_task_f1))
        for i, v in enumerate(per_task_f1):
            metrics[f"f1_task{i}"] = v

    return metrics
