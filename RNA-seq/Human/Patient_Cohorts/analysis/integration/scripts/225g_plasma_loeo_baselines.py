"""
225g_plasma_loeo_baselines.py
Leave-one-etiology-out (LOEO) evaluation, random protein baselines,
and NFASC+GDF15 2-protein benchmark on the Olink plasma data.

Parts:
  1. LOEO evaluation — top-5 models (RF, Lasso, BalancedRF, XGBoost, HistGBM)
     held out on MASLD and CVH independently; ARLD excluded (n=14, all F4).
  2. Random protein baselines — 100 draws × {20, 50, 100} proteins via
     XGBoost (best model); compared to tissue-informed panels from bridge CSV.
  3. NFASC+GDF15 2-protein benchmark — 10x repeated 5-fold CV with
     LogisticRegression; reproduce Script 82 LOEO AUROC 0.743 benchmark.

Outputs (in results/multiprogram/):
  plasma_loeo_baselines.csv
  plasma_random_baselines_v2.csv
  plasma_benchmark_comparison.csv
"""

import importlib.util, os, sys, time, warnings, logging
import numpy as np
import pandas as pd
from sklearn.model_selection import (
    RepeatedStratifiedKFold,
    StratifiedKFold,
    cross_val_score,
)
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression, Lasso
from sklearn.ensemble import RandomForestClassifier, HistGradientBoostingClassifier
from scipy.stats import mannwhitneyu

# Import per-fold imputation helper (used in LOEO and CV loops)
_here_g = os.path.dirname(os.path.abspath(__file__))
_spec_g = importlib.util.spec_from_file_location(
    "plasma_common_g", os.path.join(_here_g, "225_plasma_common.py")
)
_pc_g = importlib.util.module_from_spec(_spec_g)
_spec_g.loader.exec_module(_pc_g)
impute_train_fold = _pc_g.impute_train_fold

warnings.filterwarnings("ignore")

SEED = 42
N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))
np.random.seed(SEED)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Load shared infrastructure via importlib
# ---------------------------------------------------------------------------
_here = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "plasma_common", os.path.join(_here, "225_plasma_common.py")
)
plasma_common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plasma_common)
load_olink_data = plasma_common.load_olink_data

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
OUTDIR = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram",
)
os.makedirs(OUTDIR, exist_ok=True)

TISSUE_BRIDGE_PATH = os.path.join(
    OUTDIR, "tissue_plasma_bridge.csv"
)

# ---------------------------------------------------------------------------
# Load Olink data (replicate 225b inline logic for robustness)
# ---------------------------------------------------------------------------
log.info("Loading Olink data ...")
OLINK_PATH = os.path.join(
    BASE,
    "Analysis/Proteomics/data/olink_plasma/olink.qc.finished.mendeley.data.txt",
)
OLINK_META_PATH = os.path.join(
    BASE, "Analysis/Proteomics/results/gse276114_disease_metadata.csv"
)

olink_raw = pd.read_csv(OLINK_PATH, sep="\t")
protein_names = olink_raw["Assay"].values
npx = olink_raw.drop("Assay", axis=1).T.values.astype(np.float32)
# NaN preserved — imputation happens per train/test split

meta = pd.read_csv(OLINK_META_PATH)
meta["subject_idx"] = meta["sample_number"].astype(int) - 1
valid = meta[meta["subject_idx"] < npx.shape[0]].copy()
valid["y"] = (valid["disease_group"].isin(["F3", "F4"])).astype(int)
valid["etiology"] = valid["disease"].values

X = npx[valid["subject_idx"].values]
y = valid["y"].values
etiologies = valid["etiology"].values

log.info(
    f"Data: {X.shape[0]} subjects x {X.shape[1]} proteins, "
    f"{y.sum()} advanced / {(1 - y).sum()} early"
)
for etiol in np.unique(etiologies):
    mask = etiologies == etiol
    log.info(
        f"  {etiol}: n={mask.sum()}, F>=3={y[mask].sum()}, F<3={(1-y[mask]).sum()}"
    )

# ---------------------------------------------------------------------------
# Model definitions (top-5 from 225a/225b results)
# ---------------------------------------------------------------------------
scale_pos = float(np.sum(y == 0)) / np.sum(y == 1)

try:
    import xgboost as xgb
    HAS_XGB = True
    log.info("xgboost available")
except ImportError:
    HAS_XGB = False
    log.warning("xgboost not available — XGBoost replaced with HistGBM")

try:
    from imblearn.ensemble import BalancedRandomForestClassifier
    HAS_IMBLEARN = True
    log.info("imbalanced-learn available — BalancedRF enabled")
except ImportError:
    HAS_IMBLEARN = False
    log.warning("imbalanced-learn not available — BalancedRF replaced with weighted RF")


def make_top_models():
    models = {}

    # 1. Random Forest
    models["RF"] = Pipeline([
        ("scaler", StandardScaler()),
        ("model", RandomForestClassifier(
            n_estimators=500, class_weight="balanced",
            random_state=SEED, n_jobs=N_JOBS,
        )),
    ])

    # 2. Lasso logistic
    models["Lasso"] = Pipeline([
        ("scaler", StandardScaler()),
        ("model", LogisticRegression(
            penalty="l1", solver="saga", C=0.01,
            class_weight="balanced", max_iter=2000, random_state=SEED,
        )),
    ])

    # 3. BalancedRF (or weighted RF fallback)
    if HAS_IMBLEARN:
        models["BalancedRF"] = Pipeline([
            ("scaler", StandardScaler()),
            ("model", BalancedRandomForestClassifier(
                n_estimators=500, random_state=SEED, n_jobs=N_JOBS,
            )),
        ])
    else:
        models["BalancedRF"] = Pipeline([
            ("scaler", StandardScaler()),
            ("model", RandomForestClassifier(
                n_estimators=500, class_weight="balanced_subsample",
                random_state=SEED, n_jobs=N_JOBS,
            )),
        ])

    # 4. XGBoost (best binary model, AUROC 0.790)
    if HAS_XGB:
        models["XGBoost"] = Pipeline([
            ("scaler", StandardScaler()),
            ("model", xgb.XGBClassifier(
                n_estimators=300, max_depth=4, learning_rate=0.05,
                scale_pos_weight=scale_pos,
                use_label_encoder=False, eval_metric="logloss",
                random_state=SEED, n_jobs=N_JOBS,
                verbosity=0,
            )),
        ])
    else:
        # HistGBM fallback
        models["XGBoost"] = Pipeline([
            ("scaler", StandardScaler()),
            ("model", HistGradientBoostingClassifier(
                max_iter=300, max_depth=4, learning_rate=0.05,
                class_weight="balanced", random_state=SEED,
            )),
        ])

    # 5. HistGradientBoosting
    models["HistGBM"] = Pipeline([
        ("scaler", StandardScaler()),
        ("model", HistGradientBoostingClassifier(
            max_iter=300, max_depth=4, learning_rate=0.05,
            class_weight="balanced", random_state=SEED,
        )),
    ])

    return models


# ---------------------------------------------------------------------------
# Part 1: LOEO Evaluation
# ---------------------------------------------------------------------------
log.info("=" * 60)
log.info("Part 1: Leave-One-Etiology-Out (LOEO) evaluation")
log.info("=" * 60)

LOEO_ETIOLOGIES = ["MASLD", "CVH"]  # ARLD excluded (n=14, all F4)
top_models = make_top_models()

loeo_rows = []
for model_name, pipeline in top_models.items():
    log.info(f"  Model: {model_name}")
    holdout_preds = []
    holdout_true = []

    for holdout in LOEO_ETIOLOGIES:
        train_mask = etiologies != holdout
        test_mask = etiologies == holdout

        X_train_raw, y_train = X[train_mask], y[train_mask]
        X_test_raw, y_test = X[test_mask], y[test_mask]

        # Per-split imputation: statistics from training data only
        X_train, X_test = impute_train_fold(X_train_raw, X_test_raw)

        # Skip if test set has only one class
        if len(np.unique(y_test)) < 2:
            log.warning(
                f"    {holdout}: single-class test set (skipping AUROC)"
            )
            loeo_rows.append({
                "model": model_name,
                "holdout_etiology": holdout,
                "n_train": int(train_mask.sum()),
                "n_test": int(test_mask.sum()),
                "n_test_positive": int(y_test.sum()),
                "auroc": np.nan,
                "note": "single_class_test",
            })
            continue

        from sklearn.base import clone
        pipe = clone(pipeline)
        pipe.fit(X_train, y_train)

        if hasattr(pipe, "predict_proba"):
            prob = pipe.predict_proba(X_test)[:, 1]
        else:
            prob = pipe.decision_function(X_test)
            prob = (prob - prob.min()) / (prob.max() - prob.min() + 1e-9)

        auroc = roc_auc_score(y_test, prob)
        log.info(
            f"    {holdout}: n_train={train_mask.sum()}, "
            f"n_test={test_mask.sum()}, AUROC={auroc:.3f}"
        )

        loeo_rows.append({
            "model": model_name,
            "holdout_etiology": holdout,
            "n_train": int(train_mask.sum()),
            "n_test": int(test_mask.sum()),
            "n_test_positive": int(y_test.sum()),
            "auroc": float(auroc),
            "note": "",
        })

        holdout_preds.append(prob)
        holdout_true.append(y_test)

    # Pooled MASLD + CVH predictions
    if len(holdout_preds) >= 2:
        pooled_prob = np.concatenate(holdout_preds)
        pooled_true = np.concatenate(holdout_true)
        if len(np.unique(pooled_true)) == 2:
            pooled_auroc = roc_auc_score(pooled_true, pooled_prob)
            log.info(f"    POOLED (MASLD+CVH): AUROC={pooled_auroc:.3f}")
            loeo_rows.append({
                "model": model_name,
                "holdout_etiology": "POOLED_MASLD_CVH",
                "n_train": np.nan,
                "n_test": int(len(pooled_true)),
                "n_test_positive": int(pooled_true.sum()),
                "auroc": float(pooled_auroc),
                "note": "pooled_loeo",
            })

loeo_df = pd.DataFrame(loeo_rows)
loeo_path = os.path.join(OUTDIR, "plasma_loeo_baselines.csv")
loeo_df.to_csv(loeo_path, index=False)
log.info(f"LOEO results saved: {loeo_path}")

# ---------------------------------------------------------------------------
# Part 2: Random Protein Baselines
# ---------------------------------------------------------------------------
log.info("=" * 60)
log.info("Part 2: Random protein baselines (XGBoost, 100 draws)")
log.info("=" * 60)

N_DRAWS = 100
PANEL_SIZES = [20, 50, 100]
CV_FOLDS = 5  # single repeat (not 10×) for speed across 300 draws
rng = np.random.default_rng(SEED)

# Best model: XGBoost pipeline (fresh clone each call)
best_model_template = top_models["XGBoost"]


def run_cv_auroc(X_sub, y_arr, template, n_folds=5):
    """Single stratified k-fold CV, return mean AUROC."""
    from sklearn.base import clone
    cv = StratifiedKFold(n_splits=n_folds, shuffle=True, random_state=SEED)
    aurocs = []
    for train_idx, test_idx in cv.split(X_sub, y_arr):
        X_tr_raw, X_te_raw = X_sub[train_idx], X_sub[test_idx]
        # Per-fold imputation: statistics from training data only
        X_tr, X_te = impute_train_fold(X_tr_raw, X_te_raw)
        pipe = clone(template)
        pipe.fit(X_tr, y_arr[train_idx])
        if hasattr(pipe, "predict_proba"):
            prob = pipe.predict_proba(X_te)[:, 1]
        else:
            prob = pipe.decision_function(X_te)
            prob = (prob - prob.min()) / (prob.max() - prob.min() + 1e-9)
        if len(np.unique(y_arr[test_idx])) < 2:
            continue
        aurocs.append(roc_auc_score(y_arr[test_idx], prob))
    return float(np.mean(aurocs)) if aurocs else np.nan


# Load tissue-plasma bridge for tissue-informed panels
tissue_informed_aurocs = {}
if os.path.exists(TISSUE_BRIDGE_PATH):
    bridge = pd.read_csv(TISSUE_BRIDGE_PATH)
    log.info(f"Tissue bridge loaded: {bridge.shape[0]} rows")
    # Filter to DEGs present in plasma
    if "is_deg" in bridge.columns and "in_plasma" in bridge.columns:
        tissue_informed = bridge[
            (bridge["is_deg"] == True) & (bridge["in_plasma"] == True)
        ].copy()
    elif "dream_logFC" in bridge.columns:
        tissue_informed = bridge.dropna(subset=["dream_logFC"]).copy()
    else:
        tissue_informed = bridge.copy()

    sort_col = "dream_logFC" if "dream_logFC" in tissue_informed.columns else tissue_informed.columns[1]
    tissue_ranked = tissue_informed.reindex(
        tissue_informed[sort_col].abs().sort_values(ascending=False).index
    )

    prot_set = set(protein_names)
    symbol_col = "human_symbol" if "human_symbol" in tissue_ranked.columns else tissue_ranked.columns[0]
    tissue_proteins_all = [p for p in tissue_ranked[symbol_col] if p in prot_set]
    log.info(f"Tissue-informed proteins matching Olink: {len(tissue_proteins_all)}")
else:
    log.warning(
        f"Tissue bridge not found at {TISSUE_BRIDGE_PATH} — "
        "tissue-informed comparison will be skipped"
    )
    tissue_proteins_all = []
    tissue_ranked = pd.DataFrame()

random_rows = []
for panel_size in PANEL_SIZES:
    log.info(f"  Panel size = {panel_size}")

    # Tissue-informed panel AUROC
    if tissue_proteins_all:
        tissue_proteins = tissue_proteins_all[:panel_size]
        tissue_idx = [
            i for i, p in enumerate(protein_names) if p in set(tissue_proteins)
        ]
        if len(tissue_idx) >= 2:
            tissue_auroc = run_cv_auroc(X[:, tissue_idx], y, best_model_template)
            log.info(f"    Tissue-informed ({len(tissue_idx)} proteins): AUROC={tissue_auroc:.3f}")
        else:
            tissue_auroc = np.nan
            log.warning(f"    Too few tissue proteins matched ({len(tissue_idx)})")
    else:
        tissue_auroc = np.nan

    # Random draws
    random_aurocs = []
    for draw_idx in range(N_DRAWS):
        rand_idx = rng.choice(X.shape[1], size=panel_size, replace=False)
        rand_auroc = run_cv_auroc(X[:, rand_idx], y, best_model_template)
        random_aurocs.append(rand_auroc)
        if (draw_idx + 1) % 20 == 0:
            log.info(
                f"    Draw {draw_idx + 1}/{N_DRAWS}: "
                f"running mean AUROC={np.nanmean(random_aurocs):.3f}"
            )

    # One-sided p-value: P(random >= tissue)
    valid_aurocs = [a for a in random_aurocs if not np.isnan(a)]
    if not np.isnan(tissue_auroc) and valid_aurocs:
        p_value = (sum(r >= tissue_auroc for r in valid_aurocs) + 1) / (
            len(valid_aurocs) + 1
        )
        stat, mwu_p = mannwhitneyu(
            valid_aurocs,
            [tissue_auroc],
            alternative="greater",
        )
    else:
        p_value = np.nan
        mwu_p = np.nan

    log.info(
        f"    Random mean={np.nanmean(random_aurocs):.3f} ± {np.nanstd(random_aurocs):.3f}, "
        f"tissue={tissue_auroc:.3f}, p(random>=tissue)={p_value:.3f}"
    )

    random_rows.append({
        "panel_size": panel_size,
        "tissue_informed_auroc": tissue_auroc,
        "n_tissue_proteins_matched": len(tissue_idx) if tissue_proteins_all else np.nan,
        "random_mean_auroc": float(np.nanmean(random_aurocs)),
        "random_std_auroc": float(np.nanstd(random_aurocs)),
        "random_median_auroc": float(np.nanmedian(random_aurocs)),
        "random_q25_auroc": float(np.nanpercentile(random_aurocs, 25)),
        "random_q75_auroc": float(np.nanpercentile(random_aurocs, 75)),
        "p_random_ge_tissue": float(p_value) if not np.isnan(p_value) else np.nan,
        "n_draws": N_DRAWS,
        "random_aurocs_json": str(random_aurocs),
    })

random_df = pd.DataFrame(random_rows)
random_path = os.path.join(OUTDIR, "plasma_random_baselines_v2.csv")
random_df.to_csv(random_path, index=False)
log.info(f"Random baselines saved: {random_path}")

# ---------------------------------------------------------------------------
# Part 3: NFASC+GDF15 2-protein benchmark
# ---------------------------------------------------------------------------
log.info("=" * 60)
log.info("Part 3: NFASC+GDF15 2-protein benchmark")
log.info("=" * 60)

prot_list = list(protein_names)
nfasc_idx = prot_list.index("NFASC") if "NFASC" in prot_list else None
gdf15_idx = prot_list.index("GDF15") if "GDF15" in prot_list else None

log.info(f"NFASC index: {nfasc_idx}")
log.info(f"GDF15 index: {gdf15_idx}")

benchmark_rows = []

if nfasc_idx is not None and gdf15_idx is not None:
    X_2p = X[:, [nfasc_idx, gdf15_idx]]

    lr_2p = LogisticRegression(
        class_weight="balanced", max_iter=1000, random_state=SEED
    )

    # 10x repeated 5-fold CV with per-fold scaling (no global fit_transform)
    rskf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=SEED)
    aurocs_2p = []
    for fold_i, (tr_idx, te_idx) in enumerate(rskf.split(X_2p, y)):
        X_2p_tr, X_2p_te = impute_train_fold(X_2p[tr_idx], X_2p[te_idx])
        scaler_fold = StandardScaler()
        X_2p_tr_s = scaler_fold.fit_transform(X_2p_tr)
        X_2p_te_s = scaler_fold.transform(X_2p_te)
        lr_fold = LogisticRegression(
            class_weight="balanced", max_iter=1000, random_state=SEED
        )
        lr_fold.fit(X_2p_tr_s, y[tr_idx])
        prob = lr_fold.predict_proba(X_2p_te_s)[:, 1]
        if len(np.unique(y[te_idx])) == 2:
            aurocs_2p.append(roc_auc_score(y[te_idx], prob))

    mean_auroc_2p = float(np.mean(aurocs_2p))
    std_auroc_2p = float(np.std(aurocs_2p))
    log.info(
        f"NFASC+GDF15 10x5-fold CV AUROC: {mean_auroc_2p:.3f} ± {std_auroc_2p:.3f}"
    )

    # LOEO for NFASC+GDF15 (reproduce Script 82)
    loeo_2p_rows = []
    holdout_preds_2p = []
    holdout_true_2p = []
    for holdout in LOEO_ETIOLOGIES:
        train_mask = etiologies != holdout
        test_mask = etiologies == holdout
        X_tr_raw = X_2p[train_mask]
        X_te_raw = X_2p[test_mask]
        y_tr = y[train_mask]
        y_te = y[test_mask]

        if len(np.unique(y_te)) < 2:
            loeo_2p_rows.append({
                "protein_pair": "NFASC+GDF15",
                "holdout_etiology": holdout,
                "auroc": np.nan,
                "note": "single_class_test",
            })
            continue

        # Per-split imputation then scaling
        X_tr, X_te = impute_train_fold(X_tr_raw, X_te_raw)
        sc = StandardScaler()
        X_tr_s = sc.fit_transform(X_tr)
        X_te_s = sc.transform(X_te)
        lr = LogisticRegression(
            class_weight="balanced", max_iter=1000, random_state=SEED
        )
        lr.fit(X_tr_s, y_tr)
        prob = lr.predict_proba(X_te_s)[:, 1]
        auroc = roc_auc_score(y_te, prob)
        log.info(f"  NFASC+GDF15 LOEO {holdout}: AUROC={auroc:.3f}")
        loeo_2p_rows.append({
            "protein_pair": "NFASC+GDF15",
            "holdout_etiology": holdout,
            "auroc": float(auroc),
            "note": "",
        })
        holdout_preds_2p.append(prob)
        holdout_true_2p.append(y_te)

    if len(holdout_preds_2p) >= 2:
        pooled_prob_2p = np.concatenate(holdout_preds_2p)
        pooled_true_2p = np.concatenate(holdout_true_2p)
        if len(np.unique(pooled_true_2p)) == 2:
            pooled_auroc_2p = roc_auc_score(pooled_true_2p, pooled_prob_2p)
            log.info(f"  NFASC+GDF15 POOLED LOEO AUROC={pooled_auroc_2p:.3f} (Script 82 reference: 0.743)")
            loeo_2p_rows.append({
                "protein_pair": "NFASC+GDF15",
                "holdout_etiology": "POOLED_MASLD_CVH",
                "auroc": float(pooled_auroc_2p),
                "note": "pooled_loeo; script82_reference=0.743",
            })

    # Best model 10x5-fold CV AUROC for comparison (per-fold imputation)
    rskf_best = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=SEED)
    aurocs_best = []
    from sklearn.base import clone
    for tr_idx, te_idx in rskf_best.split(X, y):
        X_tr_best, X_te_best = impute_train_fold(X[tr_idx], X[te_idx])
        pipe = clone(best_model_template)
        pipe.fit(X_tr_best, y[tr_idx])
        if hasattr(pipe, "predict_proba"):
            prob = pipe.predict_proba(X_te_best)[:, 1]
        else:
            prob = pipe.decision_function(X_te_best)
        if len(np.unique(y[te_idx])) == 2:
            aurocs_best.append(roc_auc_score(y[te_idx], prob))
    mean_auroc_best = float(np.mean(aurocs_best))
    log.info(f"Best model (XGBoost) 10x5-fold AUROC: {mean_auroc_best:.3f}")

    benchmark_rows.extend(loeo_2p_rows)
    benchmark_rows.append({
        "protein_pair": "XGBoost_full_1461",
        "holdout_etiology": "10x5fold_CV",
        "auroc": mean_auroc_best,
        "note": "full_panel_best_model; literature_reference=0.790",
    })
    benchmark_rows.append({
        "protein_pair": "NFASC+GDF15",
        "holdout_etiology": "10x5fold_CV",
        "auroc": mean_auroc_2p,
        "note": f"2-protein; std={std_auroc_2p:.3f}",
    })
else:
    log.warning(
        "NFASC or GDF15 not found in Olink panel — "
        f"NFASC: {nfasc_idx}, GDF15: {gdf15_idx}. Skipping benchmark."
    )
    benchmark_rows.append({
        "protein_pair": "NFASC+GDF15",
        "holdout_etiology": "N/A",
        "auroc": np.nan,
        "note": "proteins_not_in_panel",
    })

bench_df = pd.DataFrame(benchmark_rows)
bench_path = os.path.join(OUTDIR, "plasma_benchmark_comparison.csv")
bench_df.to_csv(bench_path, index=False)
log.info(f"Benchmark comparison saved: {bench_path}")

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
log.info("=" * 60)
log.info("Summary")
log.info("=" * 60)
log.info(f"LOEO results:      {loeo_path}")
log.info(f"Random baselines:  {random_path}")
log.info(f"Benchmark:         {bench_path}")

if not loeo_df.empty:
    pooled = loeo_df[loeo_df["holdout_etiology"] == "POOLED_MASLD_CVH"]
    if not pooled.empty:
        log.info("LOEO pooled AUROCs (top models):")
        for _, row in pooled.iterrows():
            log.info(f"  {row['model']}: {row['auroc']:.3f}")

log.info("Done.")
