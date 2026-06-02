"""
225_lgb_minimal.py
LightGBM sweep for all targets — minimal imports (NO torch, tabpfn, flaml, catboost).
Avoids PyTorch/OpenMP deadlock by not loading the full 225_plasma_common.py module.
"""
import os, sys, time, warnings
import numpy as np
import pandas as pd
import lightgbm as lgb
from sklearn.model_selection import RepeatedStratifiedKFold, GridSearchCV, StratifiedKFold, RepeatedKFold, KFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, f1_score, cohen_kappa_score
from scipy.stats import spearmanr
warnings.filterwarnings("ignore")

SEED = 42
np.random.seed(SEED)

BASE = os.environ.get("MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
OUTDIR = os.path.join(BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram")
os.makedirs(OUTDIR, exist_ok=True)

import logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)])
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Load Olink data directly (no common module)
# ---------------------------------------------------------------------------
olink = pd.read_csv(os.path.join(BASE,
    "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt"), sep="\t")
protein_names = np.array(olink["Assay"].values)
npx = olink.drop("Assay", axis=1).T.values.astype(np.float32)
for j in range(npx.shape[1]):
    col = npx[:, j]; mask = np.isnan(col)
    if mask.any(): npx[mask, j] = np.nanmedian(col)

meta = pd.read_csv(os.path.join(BASE,
    "Analysis/Proteomics/results/gse276114_disease_metadata.csv"))
meta["subject_idx"] = meta["sample_number"].astype(int) - 1
valid = meta[meta["subject_idx"] < npx.shape[0]].copy()

X = npx[valid["subject_idx"].values]
y_binary    = (valid["disease_group"].isin(["F3","F4"])).astype(int).values
y_ordinal   = valid["disease_group"].map({"F0-2":0,"F3":1,"F4":2}).fillna(0).astype(int).values
y_eti3      = valid["disease"].map({"MASLD":0,"CVH":1,"ARLD":2}).values
y_etib      = (valid["disease"] == "MASLD").astype(int).values
log.info(f"Data: {X.shape[0]} subjects x {X.shape[1]} proteins")

# ---------------------------------------------------------------------------
# LightGBM param grid
# ---------------------------------------------------------------------------
PARAMS = {
    "model__num_leaves": [15, 31],
    "model__learning_rate": [0.01, 0.1],
    "model__n_estimators": [100, 300],
    "model__reg_lambda": [1, 10],
    "model__feature_fraction": [0.5, 0.8],
}

def run_cv(X, y, target_label, target_type):
    is_reg = (target_type == "regression")
    is_multi = (target_type == "multiclass")
    is_ord = (target_type == "ordinal")

    if is_reg:
        outer_cv = RepeatedKFold(n_splits=5, n_repeats=10, random_state=SEED)
        inner_cv = KFold(n_splits=3, shuffle=True, random_state=SEED)
        scoring = "neg_mean_absolute_error"
    else:
        outer_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=SEED)
        inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=SEED)
        scoring = "roc_auc" if not (is_multi or is_ord) else "f1_macro"

    metrics = []
    y_strat = y if not is_reg else None

    for fold_idx, (train_idx, test_idx) in enumerate(
            outer_cv.split(X, y_strat) if not is_reg else outer_cv.split(X)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # Per-fold median imputation
        for j in range(X_train.shape[1]):
            col = X_train[:, j]; mask = np.isnan(col)
            if mask.any():
                med = np.nanmedian(col)
                X_train[mask, j] = med
                X_test[np.isnan(X_test[:, j]), j] = med

        if is_reg:
            model = lgb.LGBMRegressor(random_state=SEED, verbose=-1, n_jobs=1)
            params = {k.replace("model__","model__"): v for k,v in PARAMS.items()
                      if "feature_fraction" not in k}  # regression subset
            params = {"model__num_leaves":[15,31], "model__learning_rate":[0.01,0.1],
                      "model__n_estimators":[100,300]}
        elif is_multi or is_ord:
            model = lgb.LGBMClassifier(is_unbalance=True, random_state=SEED, verbose=-1,
                                        n_jobs=1, num_class=3, objective="multiclass")
        else:
            scale_pos = float((y_train==0).sum()) / max((y_train==1).sum(), 1)
            model = lgb.LGBMClassifier(scale_pos_weight=scale_pos, random_state=SEED,
                                        verbose=-1, n_jobs=1)

        pipe = Pipeline([("scaler", StandardScaler()), ("model", model)])
        # n_jobs=1 in GridSearchCV to avoid LightGBM/joblib deadlock
        grid = GridSearchCV(pipe, PARAMS, cv=inner_cv, scoring=scoring, n_jobs=1, refit=True)
        try:
            grid.fit(X_train, y_train)
            best = grid.best_estimator_
        except Exception as e:
            log.warning(f"  Fold {fold_idx} failed: {e}")
            continue

        if is_reg:
            y_pred = best.predict(X_test)
            rho, _ = spearmanr(y_test, y_pred)
            metrics.append({"spearman": rho})
        elif is_multi or is_ord:
            y_pred = best.predict(X_test)
            f1 = f1_score(y_test, y_pred, average="macro", zero_division=0)
            qwk = cohen_kappa_score(y_test, y_pred, weights="quadratic") if is_ord else np.nan
            metrics.append({"f1_macro": f1, "qwk": qwk})
        else:
            y_prob = best.predict_proba(X_test)[:,1]
            auroc = roc_auc_score(y_test, y_prob)
            metrics.append({"auroc": auroc})

        if fold_idx % 10 == 0:
            m = metrics[-1]
            log.info(f"  Fold {fold_idx}/50: {m}")

    df_m = pd.DataFrame(metrics)
    result = {"model":"lightgbm","target":target_label}
    for col in df_m.columns:
        result[f"mean_{col}"] = df_m[col].mean()
        result[f"sd_{col}"] = df_m[col].std()
    result["n_folds"] = len(metrics)
    log.info(f"  {target_label}: {df_m.mean().to_dict()}  (n={len(metrics)})")
    return result

if __name__ == "__main__":
    log.info("=== LightGBM minimal sweep (all 5 targets) ===")
    results = []
    for y, lbl, ttype in [
        (y_binary,  "T1_binary",         "binary"),
        (y_ordinal, "T2_ordinal",         "ordinal"),
        (y_ordinal.astype(float), "T3_continuous", "regression"),
        (y_eti3,    "T4_etiology_3class", "multiclass"),
        (y_etib,    "T5_masld_binary",    "binary"),
    ]:
        valid_mask = ~np.isnan(y) if y.dtype == float else np.ones(len(y), dtype=bool)
        if not valid_mask.all():
            log.warning(f"  {lbl}: {(~valid_mask).sum()} NaN subjects excluded")
        Xv, yv = X[valid_mask], y[valid_mask].astype(int if ttype!="regression" else float)
        t0 = time.time()
        res = run_cv(Xv, yv, lbl, ttype)
        res["elapsed_sec"] = round(time.time()-t0, 1)
        results.append(res)

    df = pd.DataFrame(results)
    out = os.path.join(OUTDIR, "plasma_sweep_lgb_all_targets.csv")
    df.to_csv(out, index=False)
    log.info(f"\nResults:\n{df.to_string()}")
    log.info(f"Saved: {out}")
    log.info("Done.")
