"""
F5 ARLD-removal sensitivity (V5 verification).

Re-runs plasma T5 etiology prediction (MASLD vs CVH+ARLD) after excluding
14 ARLD samples, so the classification becomes MASLD (81) vs CVH (82).

Uses the same BalancedBagging model and the same nested CV as Script 225e.
Also runs top 5 T5 models for robustness (if present in registry).

Pass: AUROC >= 0.80 for balanced_bagging on MASLD-vs-CVH-only.

Output: verification/controls/f5_arld_removal.csv
"""
import importlib.util
import logging
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
COMMON = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/225_plasma_common.py",
)
OUT_CSV = os.path.join(
    BASE, "docs/manuscript/verification/controls/f5_arld_removal.csv"
)
os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)

# Import shared module (filename starts with digit)
_spec = importlib.util.spec_from_file_location("plasma_common", COMMON)
plasma_common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(plasma_common)
load_olink_data = plasma_common.load_olink_data
get_model_registry = plasma_common.get_model_registry
run_nested_cv = plasma_common.run_nested_cv

SEED = 42
np.random.seed(SEED)

# Candidate models to run (prioritise ones cited in manuscript and robust
# to class imbalance). balanced_bagging is the primary claim.
MODELS_TO_RUN = [
    "balanced_bagging",
    "nu_svc",
    "logistic_l2",
    "random_forest",
    "gradient_boosting",
]


def main():
    log.info("=== F5 ARLD-removal sensitivity ===")

    X, protein_names, meta, y_targets = load_olink_data()
    y_etiology_3 = y_targets["etiology_3"]       # MASLD=0, CVH=1, ARLD=2

    log.info(f"Full data: X={X.shape}, meta={meta.shape}")
    log.info(
        f"  etiology_3 counts: MASLD={np.sum(y_etiology_3==0)}, "
        f"CVH={np.sum(y_etiology_3==1)}, ARLD={np.sum(y_etiology_3==2)}"
    )

    # Remove ARLD (etiology_3 == 2)
    keep_mask = y_etiology_3 != 2
    X_sub = X[keep_mask]
    y_full = y_etiology_3[keep_mask]
    # Binary: MASLD=1, CVH=0
    y_bin = (y_full == 0).astype(int)
    log.info(f"After ARLD removal: X={X_sub.shape}, n_MASLD={y_bin.sum()}, "
             f"n_CVH={len(y_bin) - y_bin.sum()}")

    registry = get_model_registry("binary")
    available = [m for m in MODELS_TO_RUN if m in registry]
    missing = [m for m in MODELS_TO_RUN if m not in registry]
    if missing:
        log.warning(f"Missing from registry: {missing}")
    log.info(f"Running models: {available}")

    results = []
    for model_name in available:
        log.info(f"\n--- {model_name} (MASLD vs CVH only) ---")
        t0 = time.time()
        try:
            metrics = run_nested_cv(
                X_sub, y_bin, model_name, registry[model_name],
                target_type="binary",
            )
            elapsed = time.time() - t0
            metrics["elapsed_sec"] = round(elapsed, 1)
            metrics["task"] = "T5_masld_vs_cvh_ARLD_removed"
            metrics["n_MASLD"] = int(y_bin.sum())
            metrics["n_CVH"] = int(len(y_bin) - y_bin.sum())
            metrics["n_ARLD_excluded"] = int(np.sum(~keep_mask))
            results.append(metrics)
            log.info(
                f"  {model_name}: AUROC={metrics.get('mean_auroc', float('nan')):.3f} "
                f"(SD={metrics.get('sd_auroc', float('nan')):.3f}) "
                f"F1={metrics.get('mean_f1', float('nan')):.3f}  ({elapsed:.0f}s)"
            )
        except Exception as exc:
            elapsed = time.time() - t0
            log.error(f"  {model_name} FAILED: {exc}")
            results.append({
                "model": model_name,
                "task": "T5_masld_vs_cvh_ARLD_removed",
                "error": str(exc),
                "elapsed_sec": round(elapsed, 1),
            })

    # Save
    df = pd.DataFrame(results)
    log.info(f"\n=== SUMMARY ===")
    show_cols = [c for c in ("model", "mean_auroc", "sd_auroc", "mean_f1",
                              "mean_balanced_accuracy", "n_folds")
                 if c in df.columns]
    if show_cols:
        log.info(df[show_cols].to_string(index=False))
    df.to_csv(OUT_CSV, index=False)
    log.info(f"\nSaved: {OUT_CSV}")


if __name__ == "__main__":
    main()
