"""
225d_plasma_sweep_ordinal.py
Plasma biomarker model sweep — Ordinal 3-class fibrosis (T2) + Continuous fibrosis (T3)

T2: ordinal 3-class (F0-2=0, F3=1, F4=2) — evaluated by QWK
T3: continuous fibrosis score (same encoding as float) — evaluated by Spearman rho

All models in 225_plasma_common.py model registry are swept for each target type.
Results written to plasma_sweep_225d_ordinal.csv in results/multiprogram/.

Part of the 4-script parallel sweep (225a: linear/classical, 225b: boosting,
225c: neural/ensemble, 225d: ordinal & continuous). Results are combined downstream
for final model selection.

SLURM: sbatch --job-name=225d --partition=cpu --cpus-per-task=16 --mem=64G --time=48:00:00
Env: micromamba activate spatial
"""

import os
import sys
import time
import logging
import warnings

import numpy as np
import pandas as pd

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

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Import common module — filename starts with a digit so use importlib
# ---------------------------------------------------------------------------

import importlib.util

_common_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "225_plasma_common.py")
_spec = importlib.util.spec_from_file_location("plasma_common_225", _common_path)
_plasma_common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_plasma_common)

load_olink_data   = _plasma_common.load_olink_data
get_model_registry = _plasma_common.get_model_registry
run_nested_cv     = _plasma_common.run_nested_cv

# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------

log.info("=== 225d: Plasma Model Sweep — Ordinal 3-class (T2) + Continuous Fibrosis (T3) ===")

# Support MODEL_SUBSET env var for parallel per-model SLURM jobs
_subset_str = os.environ.get("MODEL_SUBSET", "")
_model_subset = set(_subset_str.split(",")) if _subset_str else None
# Support SKIP_MODELS to avoid re-running completed models
_skip_str = os.environ.get("SKIP_MODELS", "")
_skip_set = set(_skip_str.split(",")) if _skip_str else set()
# Output suffix for parallel runs
_out_suffix = os.environ.get("OUT_SUFFIX", "")

X, protein_names, meta, y_targets = load_olink_data()
y_ordinal    = y_targets["ordinal"]    # (177,) int  F0-2=0, F3=1, F4=2
y_continuous = y_targets["continuous"] # (177,) float same encoding

log.info(f"Ordinal classes: {dict(zip(*np.unique(y_ordinal, return_counts=True)))}")
log.info(f"Continuous range: [{y_continuous.min():.1f}, {y_continuous.max():.1f}]")

# ---------------------------------------------------------------------------
# Phase 1: T2 — ordinal 3-class fibrosis (QWK primary metric)
# ---------------------------------------------------------------------------

log.info("\n--- Phase 1: T2 ordinal 3-class fibrosis ---")
models_ord = get_model_registry("ordinal")
if _model_subset:
    models_ord = {k: v for k, v in models_ord.items() if k in _model_subset}
if _skip_set:
    models_ord = {k: v for k, v in models_ord.items() if k not in _skip_set}
log.info(f"Ordinal model registry: {len(models_ord)} models — {list(models_ord.keys())}")

results = []
t0_total = time.time()

for name, cfg in models_ord.items():
    log.info(f"T2 ordinal: {name} ...")
    t0 = time.time()
    try:
        res = run_nested_cv(X, y_ordinal, name, cfg, "ordinal")
        elapsed = time.time() - t0
        res["target"] = "T2_ordinal"
        res["elapsed_sec"] = round(elapsed, 1)
        results.append(res)
        qwk = res.get("mean_qwk", np.nan)
        log.info(f"  {name}: QWK={qwk:.3f}  ({elapsed:.1f}s)")
    except Exception as e:
        elapsed = time.time() - t0
        log.error(f"  {name} FAILED: {e}")
        results.append({
            "model": name,
            "target": "T2_ordinal",
            "target_type": "ordinal",
            "mean_qwk": np.nan,
            "n_folds": 0,
            "elapsed_sec": round(elapsed, 1),
        })

# ---------------------------------------------------------------------------
# Phase 2: T3 — continuous fibrosis score (Spearman primary metric)
# ---------------------------------------------------------------------------

log.info("\n--- Phase 2: T3 continuous fibrosis score ---")
models_reg = get_model_registry("regression")
if _model_subset:
    models_reg = {k: v for k, v in models_reg.items() if k in _model_subset}
if _skip_set:
    models_reg = {k: v for k, v in models_reg.items() if k not in _skip_set}
log.info(f"Regression model registry: {len(models_reg)} models — {list(models_reg.keys())}")

for name, cfg in models_reg.items():
    log.info(f"T3 continuous: {name} ...")
    t0 = time.time()
    try:
        res = run_nested_cv(X, y_continuous, name, cfg, "regression")
        elapsed = time.time() - t0
        res["target"] = "T3_continuous"
        res["elapsed_sec"] = round(elapsed, 1)
        results.append(res)
        spearman = res.get("mean_spearman", np.nan)
        log.info(f"  {name}: Spearman={spearman:.3f}  ({elapsed:.1f}s)")
    except Exception as e:
        elapsed = time.time() - t0
        log.error(f"  {name} FAILED: {e}")
        results.append({
            "model": name,
            "target": "T3_continuous",
            "target_type": "regression",
            "mean_spearman": np.nan,
            "n_folds": 0,
            "elapsed_sec": round(elapsed, 1),
        })

# ---------------------------------------------------------------------------
# Save
# ---------------------------------------------------------------------------

df = pd.DataFrame(results)
_suffix = f"_225d_ordinal{'_' + _out_suffix if _out_suffix else ''}"
out_path = os.path.join(OUTDIR, f"plasma_sweep{_suffix}.csv")
df.to_csv(out_path, index=False)

total_elapsed = time.time() - t0_total
log.info(f"\nResults:\n{df.to_string()}")
log.info(f"Saved: {out_path}")
log.info(f"Total elapsed: {total_elapsed:.1f}s")
log.info("Done.")
