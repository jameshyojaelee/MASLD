#!/usr/bin/env python
"""V6-a Test 4: F5 etiology AUROC via SIMPLER models.

Our claim: T5 MASLD-vs-other AUROC=0.840 ± 0.058 (BalancedBagging,
50-fold nested CV).
Competitors:
  - Logistic Regression with L2 penalty (default C=1.0, balanced class weights)
  - Random Forest with default parameters (class_weight=balanced)

Test: Does a simpler model achieve comparable AUROC?
- If random_forest (default) ~ 0.840 → "BalancedBagging critical" framing is moot
- If LR-L2 ~ 0.840 → even simpler

Note: the existing 226f_sweep_leaderboard.csv already has random_forest
at AUROC=0.808 (default sklearn parameters) under nested CV. We add
explicit logistic_l2 here.
"""
from __future__ import annotations
import logging
import os
import sys
import time
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
INT = f"{BASE}/RNA-seq/Human/Patient_Cohorts/analysis/integration"
sys.path.insert(0, f"{INT}/scripts")

import importlib
plasma_common = importlib.import_module("225_plasma_common")
load_olink_data = plasma_common.load_olink_data
run_nested_cv = plasma_common.run_nested_cv

OUT_DIR = f"{BASE}/docs/manuscript/verification/method_comparison"
os.makedirs(OUT_DIR, exist_ok=True)

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger(__name__)

SEED = 42
np.random.seed(SEED)
P = "model__"

# -- Load data --
log.info("Loading Olink data...")
X, protein_names, meta, y_targets = load_olink_data()
y_etiology_binary = y_targets["etiology_binary"]  # MASLD=1 vs other=0
log.info(f"X: {X.shape}, y: {y_etiology_binary.shape}, "
         f"MASLD={int((y_etiology_binary==1).sum())}, "
         f"other={int((y_etiology_binary==0).sum())}")

# -- Define competitor models (simple LR-L2 + default RF) --
models_to_run = {
    "logistic_l2_simple": {
        "model": Pipeline([
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(
                penalty="l2", C=1.0, class_weight="balanced",
                max_iter=5000, random_state=SEED, solver="lbfgs",
            )),
        ]),
        "params": {f"{P}C": [0.01, 0.1, 1.0, 10.0]},
        "is_pipeline": True,
    },
    "logistic_l2_no_grid": {  # truly simple — no hyperparameter tuning
        "model": Pipeline([
            ("scaler", StandardScaler()),
            ("model", LogisticRegression(
                penalty="l2", C=1.0, class_weight="balanced",
                max_iter=5000, random_state=SEED, solver="lbfgs",
            )),
        ]),
        "params": {},   # No grid → uses C=1.0 default
        "is_pipeline": True,
    },
    "random_forest_default": {  # default sklearn RF, no tuning
        "model": Pipeline([
            ("scaler", StandardScaler()),
            ("model", RandomForestClassifier(
                class_weight="balanced", random_state=SEED, n_jobs=1,
            )),
        ]),
        "params": {},
        "is_pipeline": True,
    },
}

# -- Run --
results = []
for model_name, model_cfg in models_to_run.items():
    log.info(f"\n=== Running {model_name} (T5 MASLD-binary, nested CV) ===")
    t0 = time.time()
    try:
        metrics = run_nested_cv(
            X, y_etiology_binary, model_name, model_cfg, target_type="binary"
        )
        elapsed = time.time() - t0
        metrics["elapsed_sec"] = round(elapsed, 1)
        results.append(metrics)
        log.info(
            f"  {model_name}: AUROC={metrics.get('mean_auroc', float('nan')):.3f} "
            f"(SD={metrics.get('sd_auroc', float('nan')):.3f}) "
            f"F1={metrics.get('mean_f1', float('nan')):.3f} ({elapsed:.0f}s)"
        )
    except Exception as exc:
        elapsed = time.time() - t0
        log.error(f"  {model_name} FAILED: {exc}")
        import traceback
        traceback.print_exc()
        results.append({
            "model": model_name, "n_folds": 0,
            "elapsed_sec": round(elapsed, 1), "error": str(exc),
        })

# -- Summary --
df = pd.DataFrame(results)
out_csv = f"{OUT_DIR}/test4_F5_simple_models.csv"
df.to_csv(out_csv, index=False)
log.info(f"\nSaved: {out_csv}")
log.info("\n=== Summary ===")
print(df.to_string(index=False))

log.info("\n=== Reference (from 226f_sweep_leaderboard.csv) ===")
log.info("balanced_bagging   AUROC=0.840 (claim)")
log.info("nu_svc             AUROC=0.834")
log.info("random_forest (already in leaderboard, default but with grid) AUROC=0.808")

# Decision rule
ours = 0.840
ref_random_forest = 0.808
ref_nu_svc = 0.834
log.info(f"\nLogistic L2 (simple, no tuning): AUROC = {df.iloc[1]['mean_auroc']:.3f}" if len(df) > 1 else "")
log.info(f"Default RF (no tuning):          AUROC = {df.iloc[2]['mean_auroc']:.3f}" if len(df) > 2 else "")
log.info(f"\nDelta vs claim ({ours:.3f}):")
for _, r in df.iterrows():
    if "mean_auroc" in r and not pd.isna(r["mean_auroc"]):
        d = r["mean_auroc"] - ours
        log.info(f"  {r['model']}: {r['mean_auroc']:.3f}  (Δ={d:+.3f})")

log.info("=== Done ===")
