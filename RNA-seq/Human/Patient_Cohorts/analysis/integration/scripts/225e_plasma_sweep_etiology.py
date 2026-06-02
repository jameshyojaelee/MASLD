"""
225e_plasma_sweep_etiology.py
Plasma model sweep — Etiology Classification (T4, T5)

Phase 1 (T4): 3-class etiology prediction — MASLD vs CVH vs ARLD.
              Primary metric: macro F1. Note: ARLD n=14 (tiny class).
Phase 2 (T5): Binary MASLD-detection — MASLD=1 vs other (CVH+ARLD)=0.
              Primary metric: AUROC.

~35 models × 2 targets via shared get_model_registry() + run_nested_cv().
10× repeated 5-fold nested CV with inner GridSearchCV. class_weight='balanced'
enforced throughout.

Output: results/multiprogram/plasma_sweep_225e_etiology.csv
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

# ---------------------------------------------------------------------------
# Bootstrap logging before common module import
# ---------------------------------------------------------------------------

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Import shared module (filename starts with digit — can't use normal import)
# ---------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "plasma_common",
    os.path.join(_HERE, "225_plasma_common.py"),
)
plasma_common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(plasma_common)

load_olink_data = plasma_common.load_olink_data
get_model_registry = plasma_common.get_model_registry
run_nested_cv = plasma_common.run_nested_cv

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

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
OUT_PATH = os.path.join(OUTDIR, "plasma_sweep_225e_etiology.csv")

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main():
    log.info("=== 225e: Plasma Model Sweep — Etiology Classification (T4, T5) ===")

    # Support MODEL_SUBSET env var for parallel per-model SLURM jobs
    _subset_str = os.environ.get("MODEL_SUBSET", "")
    _model_subset = set(_subset_str.split(",")) if _subset_str else None
    _skip_str = os.environ.get("SKIP_MODELS", "")
    _skip_set = set(_skip_str.split(",")) if _skip_str else set()
    _out_suffix = os.environ.get("OUT_SUFFIX", "")

    # Load data
    X, protein_names, meta, y_targets = load_olink_data()
    y_etiology_3 = y_targets["etiology_3"]       # MASLD=0, CVH=1, ARLD=2
    y_etiology_binary = y_targets["etiology_binary"]  # MASLD=1, other=0

    # Class counts for reference
    unique_3, counts_3 = np.unique(y_etiology_3, return_counts=True)
    unique_b, counts_b = np.unique(y_etiology_binary, return_counts=True)
    log.info(
        f"T4 (3-class) class distribution: "
        + ", ".join(
            [f"{['MASLD','CVH','ARLD'][c]}={n}" for c, n in zip(unique_3, counts_3)]
        )
    )
    log.info(
        f"T5 (binary) class distribution: "
        + ", ".join([f"{c}={n}" for c, n in zip(unique_b, counts_b)])
    )
    log.info("NOTE: ARLD n≈14 — many models will struggle; bal_acc and macro F1 penalise class failures.")

    results = []

    # -----------------------------------------------------------------------
    # Phase 1: T4 — 3-class etiology (multiclass, primary metric: macro F1)
    # -----------------------------------------------------------------------
    log.info("\n--- Phase 1: T4 — 3-class etiology (MASLD / CVH / ARLD) ---")
    registry_mc = get_model_registry("multiclass")
    if _model_subset:
        registry_mc = {k: v for k, v in registry_mc.items() if k in _model_subset}
    if _skip_set:
        registry_mc = {k: v for k, v in registry_mc.items() if k not in _skip_set}
    log.info(f"Models registered for multiclass: {list(registry_mc.keys())}")

    for model_name, model_cfg in registry_mc.items():
        log.info(f"  [T4] Running {model_name} ...")
        t0 = time.time()
        try:
            metrics = run_nested_cv(
                X, y_etiology_3, model_name, model_cfg, target_type="multiclass"
            )
            elapsed = time.time() - t0
            metrics["phase"] = "T4_etiology_3class"
            metrics["elapsed_sec"] = round(elapsed, 1)
            results.append(metrics)
            f1_str = f"F1_macro={metrics.get('mean_f1_macro', float('nan')):.3f}"
            auroc_str = f"AUROC_OvR={metrics.get('mean_auroc_ovr', float('nan')):.3f}"
            log.info(f"    {model_name}: {f1_str}  {auroc_str}  ({elapsed:.0f}s)")
        except Exception as exc:
            elapsed = time.time() - t0
            log.error(f"    {model_name} FAILED: {exc}")
            results.append(
                {
                    "model": model_name,
                    "target_type": "multiclass",
                    "phase": "T4_etiology_3class",
                    "n_folds": 0,
                    "elapsed_sec": round(elapsed, 1),
                }
            )

    # -----------------------------------------------------------------------
    # Phase 2: T5 — binary MASLD vs other (binary, primary metric: AUROC)
    # -----------------------------------------------------------------------
    log.info("\n--- Phase 2: T5 — binary MASLD detection (MASLD=1 vs CVH+ARLD=0) ---")
    registry_bin = get_model_registry("binary")
    if _model_subset:
        registry_bin = {k: v for k, v in registry_bin.items() if k in _model_subset}
    if _skip_set:
        registry_bin = {k: v for k, v in registry_bin.items() if k not in _skip_set}
    log.info(f"Models registered for binary: {list(registry_bin.keys())}")

    for model_name, model_cfg in registry_bin.items():
        log.info(f"  [T5] Running {model_name} ...")
        t0 = time.time()
        try:
            metrics = run_nested_cv(
                X, y_etiology_binary, model_name, model_cfg, target_type="binary"
            )
            elapsed = time.time() - t0
            metrics["phase"] = "T5_masld_binary"
            metrics["elapsed_sec"] = round(elapsed, 1)
            results.append(metrics)
            auroc_str = f"AUROC={metrics.get('mean_auroc', float('nan')):.3f}"
            f1_str = f"F1={metrics.get('mean_f1', float('nan')):.3f}"
            log.info(f"    {model_name}: {auroc_str}  {f1_str}  ({elapsed:.0f}s)")
        except Exception as exc:
            elapsed = time.time() - t0
            log.error(f"    {model_name} FAILED: {exc}")
            results.append(
                {
                    "model": model_name,
                    "target_type": "binary",
                    "phase": "T5_masld_binary",
                    "n_folds": 0,
                    "elapsed_sec": round(elapsed, 1),
                }
            )

    # -----------------------------------------------------------------------
    # Save
    # -----------------------------------------------------------------------
    df = pd.DataFrame(results)

    # Reorder columns: identifiers first, then metrics
    id_cols = ["phase", "model", "target_type", "n_folds", "elapsed_sec"]
    metric_cols = [c for c in df.columns if c not in id_cols]
    df = df[id_cols + metric_cols]

    _sfx = f"_225e_etiology{'_' + _out_suffix if _out_suffix else ''}"
    _out_path = os.path.join(OUTDIR, f"plasma_sweep{_sfx}.csv")
    df.to_csv(_out_path, index=False)

    # Summary tables per phase
    for phase_label in ["T4_etiology_3class", "T5_masld_binary"]:
        sub = df[df["phase"] == phase_label].copy()
        if sub.empty:
            continue
        if phase_label == "T4_etiology_3class":
            sort_col = "mean_f1_macro" if "mean_f1_macro" in sub.columns else None
            log.info(f"\n=== T4 — 3-class etiology results (sorted by macro F1) ===")
        else:
            sort_col = "mean_auroc" if "mean_auroc" in sub.columns else None
            log.info(f"\n=== T5 — binary MASLD detection results (sorted by AUROC) ===")

        if sort_col and sort_col in sub.columns:
            sub = sub.sort_values(sort_col, ascending=False)

        display_cols = ["model"] + [
            c for c in sub.columns
            if c.startswith("mean_") and c not in ("mean_elapsed_sec",)
        ]
        display_cols = [c for c in display_cols if c in sub.columns]
        log.info(sub[display_cols].to_string(index=False))

    log.info(f"\nSaved {len(df)} rows to {_out_path}")
    log.info("Done.")


if __name__ == "__main__":
    main()
