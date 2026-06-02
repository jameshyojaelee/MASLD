"""
225h_diams_within_cohort.py
DIA-MS plasma proteomics — within-cohort 4-class classification (N=72)

Independent platform validation: Normal / MASL / MASH / Cirrhosis from
PXD052937 plasma mass spectrometry data.

8 models × 3 targets × 10x5-fold RepeatedStratifiedKFold (no nested CV —
N=72 is too small for GridSearchCV; sensible fixed defaults).

Targets:
  1. four_class  : Normal(7) / MASL(17) / MASH(38) / Cirrhosis(10) — macro F1
  2. binary_dis  : Normal(7) vs Disease(65)                         — AUROC
  3. binary_sev  : Normal+MASL(24) vs MASH+Cirrhosis(48)           — AUROC

Output: results/multiprogram/plasma_diams_within_cohort.csv
"""

import os, sys, time, warnings, logging
import numpy as np
import pandas as pd
from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.impute import SimpleImputer
from sklearn.pipeline import Pipeline
from sklearn.metrics import roc_auc_score, f1_score
warnings.filterwarnings("ignore")

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
OUT_CSV = os.path.join(OUTDIR, "plasma_diams_within_cohort.csv")

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

MATRIX_PATH = os.path.join(OUTDIR, "ms_plasma_matrix.csv")
# Use Script 220's metadata (sample IDs match the protein matrix)
META_220_PATH = os.path.join(OUTDIR, "ms_plasma_metadata.csv")
# Detailed condition labels from original metadata
META_DETAIL_PATH = os.path.join(
    BASE, "Analysis/Proteomics/results/pxd052937_disease_metadata.csv"
)

log.info("Loading DIA-MS protein matrix …")
prot = pd.read_csv(MATRIX_PATH, index_col=0)  # samples × proteins
log.info(f"  Protein matrix: {prot.shape[0]} samples × {prot.shape[1]} proteins")

log.info("Loading metadata …")
meta_220 = pd.read_csv(META_220_PATH)  # sample_id, condition (letter), run_label
meta_detail = pd.read_csv(META_DETAIL_PATH)  # sample_id, condition_letter, condition_label, ...

# Map condition letters to detailed labels via the detail metadata
letter_to_label = dict(zip(meta_detail["condition_letter"], meta_detail["condition_label"]))
meta_220["condition_label"] = meta_220["condition"].map(letter_to_label)
unmapped = meta_220.loc[meta_220["condition_label"].isna(), "condition"].unique()
assert len(unmapped) == 0, f"Unmapped condition letters: {unmapped}"
log.info(f"  Metadata: {meta_220.shape[0]} rows")
log.info(f"  Condition counts:\n{meta_220['condition_label'].value_counts().to_string()}")

# Align on sample_id (both use the same run-label IDs)
meta = meta_220.set_index("sample_id")
common = prot.index.intersection(meta.index)
log.info(f"  Common samples: {len(common)}")
prot = prot.loc[common]
meta = meta.loc[common]
assert len(prot) == len(meta), "Sample mismatch after alignment"
log.info(f"  Aligned: {len(prot)} samples")

X_raw = prot.values.astype(np.float32)
labels = meta["condition_label"].values  # Normal / MASL / MASH / Cirrhosis

# ---------------------------------------------------------------------------
# Target construction
# ---------------------------------------------------------------------------

label_order = ["Normal", "MASL", "MASH", "Cirrhosis"]
label_to_int = {l: i for i, l in enumerate(label_order)}

y_four = np.array([label_to_int[l] for l in labels])  # 0-3
y_dis = (labels != "Normal").astype(int)               # 0=Normal, 1=Disease
y_sev = np.where(
    np.isin(labels, ["MASH", "Cirrhosis"]), 1,
    np.where(np.isin(labels, ["Normal", "MASL"]), 0, -1),
).astype(int)                                          # 0=early, 1=advanced

log.info(
    f"  four_class counts: {dict(zip(*np.unique(y_four, return_counts=True)))}"
)
log.info(
    f"  binary_dis counts: Normal={( y_dis==0).sum()}, Disease={(y_dis==1).sum()}"
)
log.info(
    f"  binary_sev counts: early={(y_sev==0).sum()}, advanced={(y_sev==1).sum()}"
)

TARGETS = {
    "four_class": (y_four, "macro_f1"),
    "binary_dis": (y_dis,  "auroc"),
    "binary_sev": (y_sev,  "auroc"),
}

# ---------------------------------------------------------------------------
# Model definitions (fixed defaults — no GridSearchCV)
# ---------------------------------------------------------------------------

from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from sklearn.svm import SVC
from xgboost import XGBClassifier
from lightgbm import LGBMClassifier

try:
    from imblearn.ensemble import BalancedRandomForestClassifier
    HAS_IMBLEARN = True
    log.info("imblearn available — BalancedRandomForestClassifier enabled")
except ImportError:
    HAS_IMBLEARN = False
    log.info("imblearn not available — BalancedRandomForestClassifier skipped")

N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))

def make_models(n_classes: int):
    """Return dict of model name → sklearn estimator (not in pipeline)."""
    models = {}

    # 1. Elastic net logistic
    models["elastic_net"] = LogisticRegression(
        penalty="elasticnet", solver="saga", l1_ratio=0.5,
        C=0.1, max_iter=5000, class_weight="balanced", random_state=SEED,
    )

    # 2. Random Forest
    models["random_forest"] = RandomForestClassifier(
        n_estimators=300, max_features="sqrt",
        class_weight="balanced", random_state=SEED, n_jobs=N_JOBS,
    )

    # 3. XGBoost — scale_pos_weight only applies for binary; for multi we use sample_weight
    xgb_kwargs = dict(
        n_estimators=200, max_depth=3, learning_rate=0.05,
        subsample=0.8, colsample_bytree=0.8,
        use_label_encoder=False, eval_metric="logloss",
        random_state=SEED, n_jobs=N_JOBS,
    )
    if n_classes == 2:
        xgb_kwargs["scale_pos_weight"] = 1  # will be set per fold
    else:
        xgb_kwargs["objective"] = "multi:softprob"
        xgb_kwargs["num_class"] = n_classes
    models["xgboost"] = XGBClassifier(**xgb_kwargs)

    # 4. LightGBM
    lgbm_kwargs = dict(
        n_estimators=200, max_depth=4, learning_rate=0.05,
        num_leaves=15, is_unbalance=True,
        random_state=SEED, n_jobs=N_JOBS, verbose=-1,
    )
    if n_classes > 2:
        lgbm_kwargs["objective"] = "multiclass"
        lgbm_kwargs["num_class"] = n_classes
    models["lightgbm"] = LGBMClassifier(**lgbm_kwargs)

    # 5. SVM-RBF
    models["svm_rbf"] = SVC(
        kernel="rbf", C=1.0, gamma="scale",
        class_weight="balanced", probability=True, random_state=SEED,
    )

    # 6. HistGradientBoosting
    models["hist_gradient_boosting"] = HistGradientBoostingClassifier(
        max_iter=200, max_depth=4, learning_rate=0.05,
        class_weight="balanced", random_state=SEED,
    )

    # 7. Balanced Random Forest (optional)
    if HAS_IMBLEARN:
        models["balanced_random_forest"] = BalancedRandomForestClassifier(
            n_estimators=300, max_features="sqrt",
            random_state=SEED, n_jobs=N_JOBS,
        )

    # 8. Lasso (L1 logistic)
    models["lasso"] = LogisticRegression(
        penalty="l1", solver="saga", C=0.1, max_iter=5000,
        class_weight="balanced", random_state=SEED,
    )

    return models


# ---------------------------------------------------------------------------
# Evaluation helpers
# ---------------------------------------------------------------------------

def score_fold(model, X_tr, y_tr, X_te, y_te, metric: str, n_classes: int):
    """Fit model, predict, return scalar score."""
    # Compute sample weights for XGBoost binary (scale_pos_weight equivalent)
    sw = None
    if isinstance(model, XGBClassifier) and n_classes == 2 and y_tr.sum() > 0:
        neg = (y_tr == 0).sum()
        pos = (y_tr == 1).sum()
        sw = np.where(y_tr == 1, neg / pos, 1.0)

    fit_kwargs = {}
    if sw is not None:
        fit_kwargs["sample_weight"] = sw

    model.fit(X_tr, y_tr, **fit_kwargs)

    if metric == "auroc":
        proba = model.predict_proba(X_te)
        if n_classes == 2:
            score = roc_auc_score(y_te, proba[:, 1])
        else:
            score = roc_auc_score(y_te, proba, multi_class="ovr", average="macro")
    else:  # macro_f1
        pred = model.predict(X_te)
        score = f1_score(y_te, pred, average="macro", zero_division=0)

    return score


# ---------------------------------------------------------------------------
# Main CV loop
# ---------------------------------------------------------------------------

results = []
N_REPEATS = 10
N_SPLITS = 5

for target_name, (y, metric) in TARGETS.items():
    n_classes = len(np.unique(y))
    log.info(f"\n{'='*60}")
    log.info(f"Target: {target_name}  |  metric: {metric}  |  n_classes: {n_classes}")

    models = make_models(n_classes)

    for model_name, model_proto in models.items():
        t0 = time.time()
        fold_scores = []

        rskf = RepeatedStratifiedKFold(
            n_splits=N_SPLITS, n_repeats=N_REPEATS, random_state=SEED
        )

        for fold_i, (tr_idx, te_idx) in enumerate(rskf.split(X_raw, y)):
            X_tr_raw, X_te_raw = X_raw[tr_idx], X_raw[te_idx]
            y_tr, y_te = y[tr_idx], y[te_idx]

            # Within-fold imputation then scaling
            imputer = SimpleImputer(strategy="median")
            X_tr_imp = imputer.fit_transform(X_tr_raw)
            X_te_imp = imputer.transform(X_te_raw)

            scaler = StandardScaler()
            X_tr_sc = scaler.fit_transform(X_tr_imp)
            X_te_sc = scaler.transform(X_te_imp)

            import copy
            m = copy.deepcopy(model_proto)

            try:
                s = score_fold(m, X_tr_sc, y_tr, X_te_sc, y_te, metric, n_classes)
                fold_scores.append(s)
            except Exception as exc:
                log.warning(f"  [{target_name}/{model_name}] fold {fold_i} FAILED: {exc}")
                fold_scores.append(np.nan)

        mean_score = np.nanmean(fold_scores)
        sd_score   = np.nanstd(fold_scores)
        elapsed    = time.time() - t0
        log.info(
            f"  {model_name:30s}  {metric}={mean_score:.4f} ± {sd_score:.4f}  "
            f"({N_REPEATS*N_SPLITS} folds, {elapsed:.1f}s)"
        )

        results.append(
            {
                "target": target_name,
                "metric": metric,
                "model": model_name,
                f"mean_{metric}": mean_score,
                f"sd_{metric}": sd_score,
                "n_folds": N_REPEATS * N_SPLITS,
                "n_samples": len(y),
                "n_classes": n_classes,
                "n_proteins": X_raw.shape[1],
            }
        )

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

df_out = pd.DataFrame(results)
df_out.to_csv(OUT_CSV, index=False)
log.info(f"\nResults saved to: {OUT_CSV}")
log.info(df_out.to_string(index=False))
