#!/usr/bin/env python3
"""
225f_plasma_sweep_multitask.py
Plasma model sweep — Multi-Task Prediction (T6) & Cross-Platform Validation

Part A — Multi-task (T6):
    Joint prediction of fibrosis ordinal (3-class) + etiology 3-class.
    Custom 10×5 nested CV (MultiOutputClassifier does not support GridSearchCV
    param routing via the pipeline's "model__" prefix, so we run without inner
    hyperparameter search for multitask models).

Part B — Cross-Platform Validation (PXD052937 DIA-MS):
    Train top-3 single-task models (from T2 ordinal results if available,
    else RF + XGBoost + Lasso) on all 177 Olink subjects, then predict
    the 72 DIA-MS subjects after projecting to the shared protein space.

Output:
    results/multiprogram/plasma_sweep_225f_multitask.csv

#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --job-name=plasma_225f
#SBATCH --output=logs/plasma_225f_%j.out
"""

import importlib.util
import os
import sys
import time
import warnings
import logging

import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.metrics import cohen_kappa_score, f1_score, accuracy_score
from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore")

SEED = 42
np.random.seed(SEED)

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
# Paths
# ---------------------------------------------------------------------------

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
SCRIPTS_DIR = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram",
)
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Load shared plasma module
# ---------------------------------------------------------------------------

spec = importlib.util.spec_from_file_location(
    "plasma_common",
    os.path.join(SCRIPTS_DIR, "225_plasma_common.py"),
)
plasma_common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(plasma_common)

load_olink_data = plasma_common.load_olink_data
load_ms_data = plasma_common.load_ms_data
get_model_registry = plasma_common.get_model_registry

# ---------------------------------------------------------------------------
# Part A — Multi-task (T6)
# ---------------------------------------------------------------------------


def run_multitask_cv(X, y_multi, y_combined, model_name, model_obj, n_repeats=10, n_splits=5):
    """10×5 repeated stratified CV for multi-output models.

    Stratification uses the combined label (y_combined = y_ordinal * 3 + y_etiology_3).
    Inner hyperparameter search is skipped because MultiOutputClassifier wraps do
    not propagate GridSearchCV param routing through the pipeline "model__" prefix.
    A StandardScaler is applied in the outer fold only.

    Returns a dict with mean ± sd of:
        - QWK for fibrosis ordinal (column 0)
        - Macro-F1 for etiology 3-class (column 1)
        - Combined (mean of the two)
    """
    outer_cv = RepeatedStratifiedKFold(
        n_splits=n_splits, n_repeats=n_repeats, random_state=SEED
    )

    qwks, f1s, combined_scores = [], [], []

    for fold_idx, (train_idx, test_idx) in enumerate(outer_cv.split(X, y_combined)):
        X_train, X_test = X[train_idx], X[test_idx]
        y_train, y_test = y_multi[train_idx], y_multi[test_idx]

        scaler = StandardScaler()
        X_train_s = scaler.fit_transform(X_train)
        X_test_s = scaler.transform(X_test)

        try:
            m = clone(model_obj)
            m.fit(X_train_s, y_train)
            preds = m.predict(X_test_s)

            # preds shape: (n_test, 2) — col 0 = fibrosis, col 1 = etiology
            qwk = cohen_kappa_score(
                y_test[:, 0], preds[:, 0], weights="quadratic"
            )
            f1 = f1_score(y_test[:, 1], preds[:, 1], average="macro", zero_division=0)
            comb = (qwk + f1) / 2.0

            qwks.append(qwk)
            f1s.append(f1)
            combined_scores.append(comb)
        except Exception as exc:
            log.warning(f"  fold {fold_idx} failed for {model_name}: {exc}")
            qwks.append(np.nan)
            f1s.append(np.nan)
            combined_scores.append(np.nan)

    n_ok = int(np.sum(~np.isnan(combined_scores)))
    return {
        "model": model_name,
        "part": "A_multitask",
        "metric_primary": "combined_qwk_f1",
        "mean_combined": float(np.nanmean(combined_scores)),
        "sd_combined": float(np.nanstd(combined_scores)),
        "mean_qwk_fibrosis": float(np.nanmean(qwks)),
        "sd_qwk_fibrosis": float(np.nanstd(qwks)),
        "mean_f1_etiology": float(np.nanmean(f1s)),
        "sd_f1_etiology": float(np.nanstd(f1s)),
        "n_folds_ok": n_ok,
        "n_folds_total": n_repeats * n_splits,
    }


def run_part_a(X, y_targets):
    log.info("=" * 60)
    log.info("Part A: Multi-task prediction (T6)")
    log.info("=" * 60)

    y_ordinal = y_targets["ordinal"]       # (177,) — 0/1/2 fibrosis
    y_etiology_3 = y_targets["etiology_3"] # (177,) — 0/1/2 etiology

    y_multi = np.column_stack([y_ordinal, y_etiology_3])  # (177, 2)
    # Combined label for stratified splits (6 classes max: 0*3+0 … 2*3+2)
    y_combined = y_ordinal * 3 + y_etiology_3

    log.info(
        f"Multi-task targets: fibrosis ordinal classes={np.unique(y_ordinal)}, "
        f"etiology classes={np.unique(y_etiology_3)}, "
        f"combined classes={np.unique(y_combined)}"
    )

    registry = get_model_registry("multitask")
    log.info(f"Models in multitask registry: {list(registry.keys())}")

    results = []
    for name, cfg in registry.items():
        log.info(f"Running {name} ...")
        t0 = time.time()
        try:
            res = run_multitask_cv(X, y_multi, y_combined, name, cfg["model"])
            elapsed = time.time() - t0
            res["elapsed_sec"] = round(elapsed, 1)
            results.append(res)
            log.info(
                f"  {name}: combined={res['mean_combined']:.3f}±{res['sd_combined']:.3f} "
                f"(QWK_fib={res['mean_qwk_fibrosis']:.3f}, "
                f"F1_etio={res['mean_f1_etiology']:.3f})  [{elapsed:.0f}s]"
            )
        except Exception as exc:
            elapsed = time.time() - t0
            log.error(f"  {name} FAILED: {exc}")
            results.append({
                "model": name,
                "part": "A_multitask",
                "metric_primary": "combined_qwk_f1",
                "mean_combined": np.nan,
                "sd_combined": np.nan,
                "mean_qwk_fibrosis": np.nan,
                "sd_qwk_fibrosis": np.nan,
                "mean_f1_etiology": np.nan,
                "sd_f1_etiology": np.nan,
                "n_folds_ok": 0,
                "n_folds_total": 50,
                "elapsed_sec": round(elapsed, 1),
            })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Part B — Cross-Platform Validation (DIA-MS PXD052937)
# ---------------------------------------------------------------------------


def _map_ms_labels(condition_labels):
    """Map DIA-MS condition_label → fibrosis ordinal (0/1/2).

    Normal, MASL → 0 (early / no fibrosis)
    MASH         → 1 (intermediate)
    Cirrhosis    → 2 (advanced)
    """
    mapping = {
        "normal": 0,
        "masl": 0,
        "mash": 1,
        "cirrhosis": 2,
    }
    return np.array(
        [mapping.get(str(c).strip().lower(), np.nan) for c in condition_labels],
        dtype=float,
    )


def _get_top3_ordinal_models():
    """Return top-3 model names to use for cross-platform validation.

    Prefer models from plasma_sweep_225*ordinal* results if available;
    fall back to RF + XGBoost + Lasso.
    """
    fallback = ["random_forest", "xgboost", "lasso"]

    ordinal_results_candidates = [
        os.path.join(
            BASE,
            "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram",
            f,
        )
        for f in [
            "plasma_sweep_225a_linear.csv",
            "plasma_sweep_225b_boosting.csv",
            "plasma_sweep_225c_neural.csv",
        ]
    ]

    frames = []
    for p in ordinal_results_candidates:
        if os.path.exists(p):
            try:
                df = pd.read_csv(p)
                frames.append(df)
            except Exception:
                pass

    if not frames:
        log.info("No prior ordinal sweep results found — using fallback top-3 models.")
        return fallback

    combined = pd.concat(frames, ignore_index=True)

    # These scripts evaluate binary AUROC; use mean_auroc as proxy for ordinal ranking
    if "mean_auroc" in combined.columns and combined["mean_auroc"].notna().any():
        top = (
            combined.sort_values("mean_auroc", ascending=False)
            .dropna(subset=["mean_auroc"])
            .head(3)["model"]
            .tolist()
        )
        if len(top) >= 1:
            log.info(f"Top-3 models from prior results: {top}")
            return top

    log.info("Could not extract top-3 from prior results — using fallback.")
    return fallback


def run_part_b(X_olink, protein_names_olink, y_targets):
    """Train on all 177 Olink subjects; predict 72 DIA-MS subjects."""
    log.info("=" * 60)
    log.info("Part B: Cross-platform validation (Olink → DIA-MS PXD052937)")
    log.info("=" * 60)

    # ---- Load DIA-MS data ----
    try:
        X_ms, protein_names_ms, ms_meta = load_ms_data()
    except Exception as exc:
        log.warning(
            f"load_ms_data() failed ({exc}). "
            "ms_plasma_matrix.csv may not exist yet (Script 220 required). "
            "Skipping Part B."
        )
        return pd.DataFrame()

    # ---- Find overlapping proteins ----
    olink_set = set(protein_names_olink)
    ms_set = set(protein_names_ms)
    overlap = sorted(olink_set & ms_set)

    if len(overlap) == 0:
        log.warning(
            "No protein name overlap between Olink and DIA-MS. "
            "DIA-MS may use UniProt IDs instead of gene symbols. "
            "Skipping Part B."
        )
        return pd.DataFrame()

    log.info(
        f"Protein overlap: {len(overlap)} shared "
        f"(Olink={len(olink_set)}, DIA-MS={len(ms_set)})"
    )

    olink_idx = [protein_names_olink.index(p) for p in overlap]
    ms_idx = [protein_names_ms.index(p) for p in overlap]

    X_olink_shared = X_olink[:, olink_idx]
    X_ms_shared = X_ms[:, ms_idx]

    # ---- Map DIA-MS labels ----
    if "condition_label" not in ms_meta.columns:
        log.warning("ms_meta missing 'condition_label' column. Skipping Part B.")
        return pd.DataFrame()

    y_ms_float = _map_ms_labels(ms_meta["condition_label"].values)
    valid_mask = ~np.isnan(y_ms_float)
    if valid_mask.sum() == 0:
        log.warning("No DIA-MS samples with mappable condition labels. Skipping Part B.")
        return pd.DataFrame()

    y_ms = y_ms_float[valid_mask].astype(int)
    X_ms_valid = X_ms_shared[valid_mask]
    log.info(
        f"DIA-MS after label mapping: {int(valid_mask.sum())} samples "
        f"(classes: {np.unique(y_ms, return_counts=True)})"
    )

    # ---- Scale: fit on Olink, transform DIA-MS ----
    scaler = StandardScaler()
    X_olink_scaled = scaler.fit_transform(X_olink_shared)
    X_ms_scaled = scaler.transform(X_ms_valid)

    # ---- Train on Olink, predict DIA-MS ----
    y_olink = y_targets["ordinal"]

    top3_names = _get_top3_ordinal_models()
    # Build a minimal registry of single-task ordinal models for the top-3
    ordinal_registry = get_model_registry("ordinal")

    results = []
    for model_name in top3_names:
        if model_name not in ordinal_registry:
            log.warning(f"  {model_name} not in ordinal registry — skipping.")
            continue
        cfg = ordinal_registry[model_name]
        log.info(f"Cross-platform: training {model_name} on Olink ({X_olink_scaled.shape[0]} samples) ...")
        t0 = time.time()
        try:
            m = clone(cfg["model"])
            # Models in registry may already be Pipeline with scaler; since we
            # pre-scaled, pass X already scaled. Models without internal scaler
            # will still work; models with internal scaler will double-scale,
            # which is acceptable for rank-based methods. For correctness, wrap
            # only if the model is NOT already a Pipeline.
            if isinstance(m, Pipeline):
                # Re-create pipeline without the scaler step (data already scaled)
                steps_no_scaler = [(n, s) for n, s in m.steps if n != "scaler"]
                if steps_no_scaler:
                    m = Pipeline(steps_no_scaler)
                # else fall through — fit as-is

            m.fit(X_olink_scaled, y_olink)
            preds = m.predict(X_ms_scaled)

            qwk = cohen_kappa_score(y_ms, preds, weights="quadratic")
            acc = accuracy_score(y_ms, preds)
            f1 = f1_score(y_ms, preds, average="macro", zero_division=0)
            elapsed = time.time() - t0

            log.info(
                f"  {model_name}: QWK={qwk:.3f}, Acc={acc:.3f}, F1_macro={f1:.3f}  [{elapsed:.0f}s]"
            )
            results.append({
                "model": model_name,
                "part": "B_crossplatform",
                "metric_primary": "qwk",
                "n_olink_train": X_olink_scaled.shape[0],
                "n_ms_test": int(valid_mask.sum()),
                "n_shared_proteins": len(overlap),
                "qwk": qwk,
                "accuracy": acc,
                "f1_macro": f1,
                "elapsed_sec": round(elapsed, 1),
            })
        except Exception as exc:
            elapsed = time.time() - t0
            log.error(f"  {model_name} FAILED: {exc}")
            results.append({
                "model": model_name,
                "part": "B_crossplatform",
                "metric_primary": "qwk",
                "n_olink_train": X_olink_scaled.shape[0],
                "n_ms_test": int(valid_mask.sum()),
                "n_shared_proteins": len(overlap),
                "qwk": np.nan,
                "accuracy": np.nan,
                "f1_macro": np.nan,
                "elapsed_sec": round(elapsed, 1),
            })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    log.info("=== 225f: Plasma Model Sweep — Multi-Task & Cross-Platform Validation ===")

    # ---- Load Olink data (shared across both parts) ----
    X, protein_names, meta, y_targets = load_olink_data()

    # ---- Part A ----
    df_a = run_part_a(X, y_targets)

    # ---- Part B ----
    df_b = run_part_b(X, protein_names, y_targets)

    # ---- Combine and save ----
    frames = [df for df in [df_a, df_b] if not df.empty]
    if frames:
        df_out = pd.concat(frames, ignore_index=True)
    else:
        df_out = pd.DataFrame()

    out_path = os.path.join(OUTDIR, "plasma_sweep_225f_multitask.csv")
    df_out.to_csv(out_path, index=False)
    log.info(f"\nResults saved to {out_path}")

    if not df_a.empty:
        log.info("\n--- Part A Summary (Multi-task T6) ---")
        cols_a = [c for c in ["model", "mean_combined", "sd_combined",
                               "mean_qwk_fibrosis", "mean_f1_etiology", "elapsed_sec"]
                  if c in df_a.columns]
        log.info(df_a[cols_a].sort_values("mean_combined", ascending=False).to_string(index=False))

    if not df_b.empty:
        log.info("\n--- Part B Summary (Cross-platform DIA-MS) ---")
        cols_b = [c for c in ["model", "qwk", "accuracy", "f1_macro",
                               "n_shared_proteins", "elapsed_sec"]
                  if c in df_b.columns]
        log.info(df_b[cols_b].sort_values("qwk", ascending=False).to_string(index=False))

    log.info("Done.")
