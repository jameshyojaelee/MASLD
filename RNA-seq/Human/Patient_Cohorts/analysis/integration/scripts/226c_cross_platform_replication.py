"""
226c_cross_platform_replication.py
Cross-platform transfer: Olink-trained severity models tested on DIA-MS (PXD052937)

DIA-MS has no CVH/ARLD, so only SEVERITY (binary: early vs advanced) can be
validated cross-platform.

Four experiments:
  1. Full Olink baseline — XGBoost on all 1,461 Olink proteins, 10x5-fold nested CV
  2. Shared-only Olink  — XGBoost restricted to shared proteins, 10x5-fold nested CV
  3. Cross-platform transfer — Train on ALL Olink (shared, rank-normalised),
     predict on ALL DIA-MS (shared, rank-normalised). XGBoost + LogisticRegression.
  4. Protein-level effect-size concordance — Cohen's d correlation (model-free)

Outputs (results/multiprogram/):
  226c_cross_platform.csv      — per-experiment AUROC summary
  226c_effect_concordance.csv  — per-protein Cohen's d on both platforms
  226c_shared_proteins.csv     — shared protein inventory

SLURM: sbatch --job-name=226c --partition=cpu --cpus-per-task=8 --mem=32G --time=04:00:00
Env: micromamba activate spatial
"""

import os
import sys
import time
import warnings
import logging
import importlib.util

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, ttest_ind

from sklearn.model_selection import RepeatedStratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline

import xgboost as xgb

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

N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", 8))

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
OUTDIR = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram",
)
os.makedirs(OUTDIR, exist_ok=True)

# ---------------------------------------------------------------------------
# Import per-fold imputation helper from common module
# ---------------------------------------------------------------------------

_here = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "plasma_common_225", os.path.join(_here, "225_plasma_common.py")
)
_plasma_common = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_plasma_common)
load_olink_data = _plasma_common.load_olink_data
impute_train_fold = _plasma_common.impute_train_fold

# ---------------------------------------------------------------------------
# Data loading — Olink
# ---------------------------------------------------------------------------

log.info("=" * 60)
log.info("Loading Olink data ...")
X_olink, olink_protein_names, olink_meta, y_targets = load_olink_data()
y_olink = y_targets["binary"]  # F0-2 = 0 (early, n=39), F3+F4 = 1 (advanced, n=138)
olink_protein_names = list(olink_protein_names)

log.info(
    f"Olink: {X_olink.shape[0]} subjects x {X_olink.shape[1]} proteins | "
    f"early={int((y_olink == 0).sum())} advanced={int((y_olink == 1).sum())}"
)

# ---------------------------------------------------------------------------
# Data loading — DIA-MS
# ---------------------------------------------------------------------------

log.info("Loading DIA-MS data ...")

MS_MATRIX_PATH = os.path.join(OUTDIR, "ms_plasma_matrix.csv")
MS_META_220_PATH = os.path.join(OUTDIR, "ms_plasma_metadata.csv")
MS_META_DETAIL_PATH = os.path.join(
    BASE, "Analysis/Proteomics/results/pxd052937_disease_metadata.csv"
)

prot_df = pd.read_csv(MS_MATRIX_PATH, index_col=0)
diams_protein_names = list(prot_df.columns)
log.info(f"DIA-MS protein matrix: {prot_df.shape[0]} samples x {prot_df.shape[1]} proteins")

meta_220 = pd.read_csv(MS_META_220_PATH)
meta_detail = pd.read_csv(MS_META_DETAIL_PATH)

# Map condition letters to detailed labels
letter_to_label = dict(zip(meta_detail["condition_letter"], meta_detail["condition_label"]))
meta_220["condition_label"] = meta_220["condition"].map(letter_to_label)
unmapped = meta_220.loc[meta_220["condition_label"].isna(), "condition"].unique()
assert len(unmapped) == 0, f"Unmapped condition letters: {unmapped}"

# Align samples
meta_ms = meta_220.set_index("sample_id")
common_samples = prot_df.index.intersection(meta_ms.index)
prot_df = prot_df.loc[common_samples]
meta_ms = meta_ms.loc[common_samples]
assert len(prot_df) == len(meta_ms), "Sample mismatch after alignment"

X_diams = prot_df.values.astype(np.float32)
labels_ms = meta_ms["condition_label"].values

# DIA-MS severity target: Normal+MASL -> 0 (early), MASH+Cirrhosis -> 1 (advanced)
y_diams = np.where(
    np.isin(labels_ms, ["MASH", "Cirrhosis"]), 1,
    np.where(np.isin(labels_ms, ["Normal", "MASL"]), 0, -1),
).astype(int)
assert (y_diams >= 0).all(), "Unmapped DIA-MS condition labels"

log.info(
    f"DIA-MS: {X_diams.shape[0]} samples x {X_diams.shape[1]} proteins | "
    f"early={int((y_diams == 0).sum())} advanced={int((y_diams == 1).sum())}"
)
log.info(f"  Condition counts: {dict(zip(*np.unique(labels_ms, return_counts=True)))}")

# ---------------------------------------------------------------------------
# Shared proteins
# ---------------------------------------------------------------------------

shared = sorted(set(olink_protein_names) & set(diams_protein_names))
log.info(f"Shared proteins: {len(shared)} / Olink {len(olink_protein_names)} / DIA-MS {len(diams_protein_names)}")
log.info(f"  First 20: {shared[:20]}")

# Build index arrays for subsetting each platform to shared proteins
olink_shared_idx = [olink_protein_names.index(p) for p in shared]
diams_shared_idx = [diams_protein_names.index(p) for p in shared]

# Save shared protein inventory
shared_df = pd.DataFrame({
    "protein": sorted(set(olink_protein_names) | set(diams_protein_names)),
})
shared_df["in_olink"] = shared_df["protein"].isin(olink_protein_names)
shared_df["in_diams"] = shared_df["protein"].isin(diams_protein_names)

# Compute per-protein means for informational purposes
olink_means = {}
for i, p in enumerate(olink_protein_names):
    olink_means[p] = float(np.nanmean(X_olink[:, i]))
diams_means = {}
for i, p in enumerate(diams_protein_names):
    diams_means[p] = float(np.nanmean(X_diams[:, i]))

shared_df["olink_mean_npx"] = shared_df["protein"].map(olink_means)
shared_df["diams_mean_log2"] = shared_df["protein"].map(diams_means)

shared_df.to_csv(os.path.join(OUTDIR, "226c_shared_proteins.csv"), index=False)
log.info(f"Saved shared protein inventory: {len(shared_df)} total proteins")

# ---------------------------------------------------------------------------
# Helper: rank normalisation per column
# ---------------------------------------------------------------------------


def rank_normalize(X):
    """Convert each column to percentile ranks in [0, 1].

    NaN values are preserved: rankdata assigns them the mean rank, then they
    are replaced by the column median rank.
    """
    X_out = np.empty_like(X, dtype=np.float64)
    for j in range(X.shape[1]):
        col = X[:, j].copy()
        mask = np.isnan(col)
        if mask.all():
            X_out[:, j] = 0.5
            continue
        # Impute NaN temporarily with column median for ranking
        med = float(np.nanmedian(col))
        col[mask] = med
        X_out[:, j] = rankdata(col) / len(col)
    return X_out.astype(np.float32)


# ---------------------------------------------------------------------------
# Helper: XGBoost with class imbalance
# ---------------------------------------------------------------------------

scale_pos_olink = float(np.sum(y_olink == 0)) / max(np.sum(y_olink == 1), 1)


def make_xgb(scale_pos_weight=None):
    """Return an XGBClassifier with sensible defaults."""
    return xgb.XGBClassifier(
        n_estimators=300,
        max_depth=5,
        learning_rate=0.05,
        subsample=0.8,
        colsample_bytree=0.8,
        scale_pos_weight=scale_pos_weight or scale_pos_olink,
        eval_metric="logloss",
        use_label_encoder=False,
        random_state=SEED,
        n_jobs=N_JOBS,
    )


# ---------------------------------------------------------------------------
# Collect results
# ---------------------------------------------------------------------------

results = []

# ===================================================================
# Experiment 1: Full Olink baseline (XGBoost, all 1,461 proteins)
# ===================================================================

log.info("\n" + "=" * 60)
log.info("EXPERIMENT 1: Full Olink baseline (XGBoost, all 1,461 proteins)")
log.info("=" * 60)

t0 = time.time()
outer_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=SEED)
aurocs_exp1 = []

for fold_idx, (tr_idx, te_idx) in enumerate(outer_cv.split(X_olink, y_olink)):
    X_tr_raw, X_te_raw = X_olink[tr_idx], X_olink[te_idx]
    y_tr, y_te = y_olink[tr_idx], y_olink[te_idx]

    # Per-fold imputation
    X_tr, X_te = impute_train_fold(X_tr_raw, X_te_raw)

    # Scale
    scaler = StandardScaler()
    X_tr_sc = scaler.fit_transform(X_tr)
    X_te_sc = scaler.transform(X_te)

    model = make_xgb(scale_pos_weight=scale_pos_olink)
    model.fit(X_tr_sc, y_tr)
    y_prob = model.predict_proba(X_te_sc)[:, 1]
    aurocs_exp1.append(roc_auc_score(y_te, y_prob))

mean_auroc_1 = np.mean(aurocs_exp1)
sd_auroc_1 = np.std(aurocs_exp1)
elapsed_1 = time.time() - t0

log.info(f"  Full Olink XGBoost: AUROC = {mean_auroc_1:.4f} +/- {sd_auroc_1:.4f} ({elapsed_1:.1f}s)")

results.append({
    "experiment": "full_olink",
    "auroc": round(mean_auroc_1, 4),
    "sd_auroc": round(sd_auroc_1, 4),
    "n_train": int(X_olink.shape[0]),
    "n_test": int(X_olink.shape[0]),
    "n_proteins": int(X_olink.shape[1]),
    "notes": "10x5-fold nested CV, all 1461 Olink proteins",
})

# ===================================================================
# Experiment 2: Shared-only Olink (XGBoost, shared proteins)
# ===================================================================

log.info("\n" + "=" * 60)
log.info(f"EXPERIMENT 2: Shared-only Olink (XGBoost, {len(shared)} shared proteins)")
log.info("=" * 60)

X_olink_shared = X_olink[:, olink_shared_idx]

t0 = time.time()
outer_cv = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=SEED)
aurocs_exp2 = []

for fold_idx, (tr_idx, te_idx) in enumerate(outer_cv.split(X_olink_shared, y_olink)):
    X_tr_raw, X_te_raw = X_olink_shared[tr_idx], X_olink_shared[te_idx]
    y_tr, y_te = y_olink[tr_idx], y_olink[te_idx]

    X_tr, X_te = impute_train_fold(X_tr_raw, X_te_raw)

    scaler = StandardScaler()
    X_tr_sc = scaler.fit_transform(X_tr)
    X_te_sc = scaler.transform(X_te)

    model = make_xgb(scale_pos_weight=scale_pos_olink)
    model.fit(X_tr_sc, y_tr)
    y_prob = model.predict_proba(X_te_sc)[:, 1]
    aurocs_exp2.append(roc_auc_score(y_te, y_prob))

mean_auroc_2 = np.mean(aurocs_exp2)
sd_auroc_2 = np.std(aurocs_exp2)
elapsed_2 = time.time() - t0

log.info(f"  Shared-only Olink XGBoost: AUROC = {mean_auroc_2:.4f} +/- {sd_auroc_2:.4f} ({elapsed_2:.1f}s)")
log.info(f"  Degradation from full: {mean_auroc_1 - mean_auroc_2:+.4f}")

results.append({
    "experiment": "shared_olink",
    "auroc": round(mean_auroc_2, 4),
    "sd_auroc": round(sd_auroc_2, 4),
    "n_train": int(X_olink.shape[0]),
    "n_test": int(X_olink.shape[0]),
    "n_proteins": len(shared),
    "notes": f"10x5-fold nested CV, {len(shared)} shared proteins only",
})

# ===================================================================
# Experiment 3: Cross-platform transfer (rank-normalised)
# ===================================================================

log.info("\n" + "=" * 60)
log.info("EXPERIMENT 3: Cross-platform transfer (rank-normalised)")
log.info("=" * 60)

# Subset to shared proteins
X_olink_sh = X_olink[:, olink_shared_idx]
X_diams_sh = X_diams[:, diams_shared_idx]

# Impute NaN before rank normalisation (column median)
def impute_median(X):
    """Fill NaN with column median."""
    X_out = X.copy()
    for j in range(X.shape[1]):
        mask = np.isnan(X_out[:, j])
        if mask.any():
            med = float(np.nanmedian(X_out[:, j]))
            X_out[mask, j] = med
    return X_out


X_olink_sh_imp = impute_median(X_olink_sh)
X_diams_sh_imp = impute_median(X_diams_sh)

# Rank-normalise each platform independently
X_olink_rank = rank_normalize(X_olink_sh_imp)
X_diams_rank = rank_normalize(X_diams_sh_imp)

log.info(f"  Rank-normalised Olink: {X_olink_rank.shape}")
log.info(f"  Rank-normalised DIA-MS: {X_diams_rank.shape}")

# --- 3a: XGBoost transfer ---
log.info("  3a: XGBoost transfer ...")
t0 = time.time()

scale_pos_diams = float(np.sum(y_diams == 0)) / max(np.sum(y_diams == 1), 1)
xgb_transfer = make_xgb(scale_pos_weight=scale_pos_olink)
xgb_transfer.fit(X_olink_rank, y_olink)

y_prob_xgb = xgb_transfer.predict_proba(X_diams_rank)[:, 1]
auroc_xgb_transfer = roc_auc_score(y_diams, y_prob_xgb)
elapsed_3a = time.time() - t0

log.info(f"  XGBoost transfer AUROC = {auroc_xgb_transfer:.4f} ({elapsed_3a:.1f}s)")

results.append({
    "experiment": "transfer_xgb",
    "auroc": round(auroc_xgb_transfer, 4),
    "sd_auroc": np.nan,
    "n_train": int(X_olink.shape[0]),
    "n_test": int(X_diams.shape[0]),
    "n_proteins": len(shared),
    "notes": "Train all Olink, test all DIA-MS, rank-normalised, shared proteins",
})

# --- 3b: Logistic Regression transfer ---
log.info("  3b: LogisticRegression transfer ...")
t0 = time.time()

lr_transfer = Pipeline([
    ("scaler", StandardScaler()),
    ("model", LogisticRegression(
        penalty="l2", C=0.1, max_iter=5000,
        class_weight="balanced", solver="lbfgs", random_state=SEED,
    )),
])
lr_transfer.fit(X_olink_rank, y_olink)

y_prob_lr = lr_transfer.predict_proba(X_diams_rank)[:, 1]
auroc_lr_transfer = roc_auc_score(y_diams, y_prob_lr)
elapsed_3b = time.time() - t0

log.info(f"  LogisticRegression transfer AUROC = {auroc_lr_transfer:.4f} ({elapsed_3b:.1f}s)")

results.append({
    "experiment": "transfer_logreg",
    "auroc": round(auroc_lr_transfer, 4),
    "sd_auroc": np.nan,
    "n_train": int(X_olink.shape[0]),
    "n_test": int(X_diams.shape[0]),
    "n_proteins": len(shared),
    "notes": "Train all Olink, test all DIA-MS, rank-normalised, shared proteins, L2 LR",
})

# ===================================================================
# Experiment 4: Protein-level effect size concordance
# ===================================================================

log.info("\n" + "=" * 60)
log.info("EXPERIMENT 4: Protein-level effect-size concordance (Cohen's d)")
log.info("=" * 60)

def cohens_d(x0, x1):
    """Cohen's d: (mean_advanced - mean_early) / pooled_sd."""
    n0, n1 = len(x0), len(x1)
    if n0 < 2 or n1 < 2:
        return np.nan
    pooled_sd = np.sqrt((np.std(x0, ddof=1) ** 2 + np.std(x1, ddof=1) ** 2) / 2)
    if pooled_sd == 0:
        return 0.0
    return (np.mean(x1) - np.mean(x0)) / pooled_sd


concordance_rows = []

for i, prot_name in enumerate(shared):
    # Olink — early (0) vs advanced (1)
    ok_early = X_olink_sh_imp[y_olink == 0, i]
    ok_adv = X_olink_sh_imp[y_olink == 1, i]

    # DIA-MS — early (0) vs advanced (1)
    ms_early = X_diams_sh_imp[y_diams == 0, i]
    ms_adv = X_diams_sh_imp[y_diams == 1, i]

    olink_d = cohens_d(ok_early, ok_adv)
    diams_d = cohens_d(ms_early, ms_adv)

    # t-test p-values (two-sided, unequal variance)
    if len(ok_early) >= 2 and len(ok_adv) >= 2:
        _, olink_p = ttest_ind(ok_adv, ok_early, equal_var=False)
    else:
        olink_p = np.nan

    if len(ms_early) >= 2 and len(ms_adv) >= 2:
        _, diams_p = ttest_ind(ms_adv, ms_early, equal_var=False)
    else:
        diams_p = np.nan

    concordance_rows.append({
        "protein": prot_name,
        "olink_d": round(olink_d, 4) if not np.isnan(olink_d) else np.nan,
        "diams_d": round(diams_d, 4) if not np.isnan(diams_d) else np.nan,
        "olink_pvalue": olink_p,
        "diams_pvalue": diams_p,
    })

concordance_df = pd.DataFrame(concordance_rows)
concordance_df.to_csv(os.path.join(OUTDIR, "226c_effect_concordance.csv"), index=False)

# Compute Spearman correlation of Cohen's d across shared proteins
valid_mask = concordance_df["olink_d"].notna() & concordance_df["diams_d"].notna()
if valid_mask.sum() >= 3:
    rho, pval = spearmanr(
        concordance_df.loc[valid_mask, "olink_d"],
        concordance_df.loc[valid_mask, "diams_d"],
    )
    log.info(f"  Effect-size concordance: Spearman rho = {rho:.4f}, p = {pval:.2e} (n={valid_mask.sum()} proteins)")

    # Direction concordance
    same_dir = (
        (concordance_df.loc[valid_mask, "olink_d"] > 0)
        == (concordance_df.loc[valid_mask, "diams_d"] > 0)
    ).mean()
    log.info(f"  Direction concordance: {same_dir:.1%} of shared proteins")
else:
    rho, pval = np.nan, np.nan
    log.warning("  Too few valid proteins for correlation")

# Add concordance summary to results
results.append({
    "experiment": "effect_concordance",
    "auroc": round(rho, 4) if not np.isnan(rho) else np.nan,
    "sd_auroc": pval,
    "n_train": int(X_olink.shape[0]),
    "n_test": int(X_diams.shape[0]),
    "n_proteins": int(valid_mask.sum()),
    "notes": f"Spearman rho of Cohen's d (auroc col = rho, sd_auroc col = p-value)",
})

# ===================================================================
# Save summary
# ===================================================================

results_df = pd.DataFrame(results)
results_df.to_csv(os.path.join(OUTDIR, "226c_cross_platform.csv"), index=False)

log.info("\n" + "=" * 60)
log.info("SUMMARY")
log.info("=" * 60)
log.info(results_df.to_string(index=False))
log.info(f"\nOutputs saved to {OUTDIR}:")
log.info(f"  226c_cross_platform.csv     — {len(results_df)} rows")
log.info(f"  226c_effect_concordance.csv  — {len(concordance_df)} rows")
log.info(f"  226c_shared_proteins.csv     — {len(shared_df)} rows")
log.info("Done.")
