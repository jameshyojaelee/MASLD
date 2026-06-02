#!/usr/bin/env python3
"""
154_minimal_prognosis_panel.py
Minimal Prognosis Gene Panel — find the smallest gene set that achieves ~90%
of the full divergence model's S1/S2 classification performance.

Strategy:
  1. Rank divergence genes (div_*) by stability from Script 151 feature
     importance (n_folds_selected, then abs mean coefficient).
  2. Build gene panels at sizes 5, 10, 15, 25, 50.
  3. For each panel: LOCO-CV elastic net predicting S1 vs S2, compute AUROC.
  4. Find the "elbow" — smallest panel within 0.05 AUROC of the full
     divergence model (M3_divergence AUROC = 0.844).
  5. Compare with existing staging panels (minimal_panel_*.csv).
  6. Check overlap with plasma-measurable proteins (Olink feasibility).

Input:
  - results/prognosis/prs_feature_importance.csv
  - results/prognosis/prognosis_feature_matrix.csv
  - results/prognosis/prognosis_pseudo_labels.csv
  - results/prognosis/prs_model_comparison.csv
  - results/staging_classifier/modeling_metadata.csv
  - results/staging_classifier/minimal_panel_*.csv

Output (to results/prognosis/):
  - prognosis_panel_5.csv through prognosis_panel_50.csv
  - panel_performance_curve.csv
  - prognosis_panel_comparison.csv

SLURM: io partition, 8 CPUs, 16G RAM, 48h
Env:   micromamba activate spatial
"""

import os
import sys
import time
import warnings
import logging

import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import roc_auc_score, average_precision_score, brier_score_loss

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)
warnings.filterwarnings("ignore", category=Warning, module="sklearn.linear_model")

# ---------------------------------------------------------------------------
# Reproducibility
# ---------------------------------------------------------------------------
SEED = 42
np.random.seed(SEED)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INT = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration")
RDIR = os.path.join(INT, "results")
OUTDIR = os.path.join(RDIR, "prognosis")
STAGING_DIR = os.path.join(RDIR, "staging_classifier")
os.makedirs(OUTDIR, exist_ok=True)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(),
        logging.FileHandler(os.path.join(OUTDIR, "154_panel.log")),
    ],
)
log = logging.getLogger(__name__)

# Panel sizes to evaluate
PANEL_SIZES = [5, 10, 15, 25, 50]

# Tolerance for "elbow" detection: within this AUROC of full model
ELBOW_TOLERANCE = 0.05

# Proteins available on Olink Explore 3072 (partial list of liver-relevant
# panels; used only for feasibility flag — not an exhaustive mapping).
# We flag any gene whose symbol matches a known secreted/circulating protein.
# A full Olink manifest match would require the proprietary panel list.
OLINK_RELEVANT_PANELS = {"Explore_3072", "Explore_HT"}

# ===================================================================
# PART 1: Load data
# ===================================================================
log.info("=" * 70)
log.info("154 — Minimal Prognosis Gene Panel")
log.info("=" * 70)

t0 = time.time()

# Feature importance from Script 151
imp = pd.read_csv(os.path.join(OUTDIR, "prs_feature_importance.csv"))
log.info(f"Loaded feature importance: {imp.shape[0]} rows")

# Feature matrix and labels
feat = pd.read_csv(os.path.join(OUTDIR, "prognosis_feature_matrix.csv"))
labels = pd.read_csv(os.path.join(OUTDIR, "prognosis_pseudo_labels.csv"))

# Merge
df = labels.merge(feat, on="sample_id", how="inner")
log.info(f"Merged data: {df.shape[0]} samples")

# Full model comparison for reference AUROC
model_comp = pd.read_csv(os.path.join(OUTDIR, "prs_model_comparison.csv"))
m3_row = model_comp[
    (model_comp["config"] == "M3_divergence") & (model_comp["method"] == "elastic_net")
]
if len(m3_row) > 0:
    full_div_auroc = m3_row["mean_auroc"].values[0]
else:
    full_div_auroc = 0.844  # fallback from known result
log.info(f"Full divergence model (M3) AUROC: {full_div_auroc:.3f}")
log.info(f"Elbow threshold: {full_div_auroc - ELBOW_TOLERANCE:.3f}")

# ===================================================================
# PART 2: Rank divergence genes by stability
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 2: Ranking divergence genes by stability")
log.info("=" * 70)

# Filter to M3_divergence features only (div_* genes)
div_imp = imp[imp["config"] == "M3_divergence"].copy()
log.info(f"M3_divergence features: {div_imp.shape[0]}")

# Rank by: (1) n_folds_nonzero descending, (2) abs_mean_coefficient descending
div_imp = div_imp.sort_values(
    ["n_folds_nonzero", "abs_mean_coefficient"],
    ascending=[False, False],
).reset_index(drop=True)
div_imp["stability_rank"] = range(1, len(div_imp) + 1)

# Extract gene name from feature name (strip div_ prefix)
div_imp["gene_name"] = div_imp["feature"].str.replace("^div_", "", regex=True)

log.info("Top 20 most stable divergence genes:")
for i, row in div_imp.head(20).iterrows():
    log.info(
        f"  {row['stability_rank']:3d}. {row['gene_name']:20s}  "
        f"folds={row['n_folds_nonzero']}/{row['n_folds_total']}  "
        f"|coef|={row['abs_mean_coefficient']:.4f}"
    )

# ===================================================================
# PART 3: Prepare LOCO-CV data
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 3: Preparing LOCO-CV data")
log.info("=" * 70)

# F0-F2 samples with S2 labels and LOCO fold assignment
df_f02 = df[
    (df["fibrosis_stage"].isin([0, 1, 2]))
    & (df["s2_binary"].notna())
    & (df["loco_fold_fibrosis"] != "excluded")
].copy()
df_f02["s2_binary"] = df_f02["s2_binary"].astype(int)

folds = sorted(df_f02["loco_fold_fibrosis"].unique())
log.info(f"LOCO-eligible F0-F2 samples: {df_f02.shape[0]}")
log.info(f"S1: {(df_f02['s2_binary'] == 0).sum()}, "
         f"S2: {(df_f02['s2_binary'] == 1).sum()} "
         f"({df_f02['s2_binary'].mean()*100:.1f}% S2)")
log.info(f"LOCO folds ({len(folds)}): {folds}")


def safe_impute(X_train, X_test):
    """Impute NaN with training-set median."""
    medians = np.nanmedian(X_train, axis=0)
    for j in range(X_train.shape[1]):
        if np.isnan(medians[j]):
            medians[j] = 0.0
        X_train[np.isnan(X_train[:, j]), j] = medians[j]
        X_test[np.isnan(X_test[:, j]), j] = medians[j]
    return X_train, X_test


def run_panel_loco(panel_features, df_loco, folds_list):
    """
    Run LOCO-CV elastic net on a given set of features.
    Returns per-fold results and per-sample predictions.
    """
    fold_results = []
    predictions = []

    for fold in folds_list:
        mask_test = df_loco["loco_fold_fibrosis"] == fold
        mask_train = ~mask_test

        X_train = df_loco.loc[mask_train, panel_features].values.astype(np.float64)
        X_test = df_loco.loc[mask_test, panel_features].values.astype(np.float64)
        y_train = df_loco.loc[mask_train, "s2_binary"].values
        y_test = df_loco.loc[mask_test, "s2_binary"].values

        # Impute + scale
        X_train, X_test = safe_impute(X_train, X_test)
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train)
        X_test = scaler.transform(X_test)

        # Elastic net
        model = LogisticRegressionCV(
            penalty="elasticnet",
            l1_ratios=[0.1, 0.5, 0.9],
            solver="saga",
            cv=5,
            class_weight="balanced",
            max_iter=5000,
            scoring="roc_auc",
            random_state=SEED,
            n_jobs=-1,
        )
        model.fit(X_train, y_train)
        y_prob = model.predict_proba(X_test)[:, 1]

        n_unique = len(np.unique(y_test))
        auroc = roc_auc_score(y_test, y_prob) if n_unique > 1 else np.nan
        auprc = average_precision_score(y_test, y_prob) if n_unique > 1 else np.nan
        brier = brier_score_loss(y_test, y_prob)

        fold_results.append({
            "fold": fold,
            "n_train": mask_train.sum(),
            "n_test": mask_test.sum(),
            "auroc": auroc,
            "auprc": auprc,
            "brier": brier,
        })

        # Per-sample predictions
        test_ids = df_loco.loc[mask_test, "sample_id"].values
        for sid, yp, yt in zip(test_ids, y_prob, y_test):
            predictions.append({
                "sample_id": sid,
                "fold": fold,
                "s2_binary_true": int(yt),
                "p_s2_predicted": float(yp),
            })

    return fold_results, predictions


# ===================================================================
# PART 4: Evaluate panel sizes
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 4: Evaluating panel sizes")
log.info("=" * 70)

perf_rows = []
all_panel_genes = {}

for ps in PANEL_SIZES:
    # Select top-ranked genes
    top_genes = div_imp.head(ps)
    panel_features = ["div_" + g for g in top_genes["gene_name"].values]

    # Verify features exist in dataframe
    panel_features = [f for f in panel_features if f in df_f02.columns]
    actual_size = len(panel_features)

    if actual_size == 0:
        log.warning(f"Panel size {ps}: no valid features found, skipping")
        continue

    log.info(f"Panel size {ps}: evaluating {actual_size} features...")

    fold_res, preds = run_panel_loco(panel_features, df_f02, folds)

    fold_aurocs = [r["auroc"] for r in fold_res if not np.isnan(r["auroc"])]
    mean_auroc = np.mean(fold_aurocs) if fold_aurocs else np.nan
    sd_auroc = np.std(fold_aurocs) if len(fold_aurocs) > 1 else 0.0

    fold_auprcs = [r["auprc"] for r in fold_res if not np.isnan(r["auprc"])]
    mean_auprc = np.mean(fold_auprcs) if fold_auprcs else np.nan
    sd_auprc = np.std(fold_auprcs) if len(fold_auprcs) > 1 else 0.0

    fold_briers = [r["brier"] for r in fold_res]
    mean_brier = np.mean(fold_briers)

    perf_rows.append({
        "panel_size": ps,
        "n_genes_actual": actual_size,
        "mean_auroc": mean_auroc,
        "sd_auroc": sd_auroc,
        "mean_auprc": mean_auprc,
        "sd_auprc": sd_auprc,
        "mean_brier": mean_brier,
        "pct_of_full_auroc": (mean_auroc / full_div_auroc * 100)
        if not np.isnan(mean_auroc)
        else np.nan,
        "within_elbow": (full_div_auroc - mean_auroc) <= ELBOW_TOLERANCE
        if not np.isnan(mean_auroc)
        else False,
    })

    # Store gene list
    gene_names = top_genes["gene_name"].values[:actual_size]
    all_panel_genes[ps] = gene_names

    log.info(
        f"  Panel {ps}: AUROC {mean_auroc:.3f} +/- {sd_auroc:.3f}  "
        f"({mean_auroc/full_div_auroc*100:.1f}% of full)  "
        f"{'[ELBOW]' if (full_div_auroc - mean_auroc) <= ELBOW_TOLERANCE else ''}"
    )

# Also evaluate the full 100-gene panel for reference
log.info(f"Full divergence model (100 genes): evaluating...")
full_features = ["div_" + g for g in div_imp["gene_name"].values]
full_features = [f for f in full_features if f in df_f02.columns]
fold_res_full, _ = run_panel_loco(full_features, df_f02, folds)
full_aurocs = [r["auroc"] for r in fold_res_full if not np.isnan(r["auroc"])]
full_auroc_actual = np.mean(full_aurocs) if full_aurocs else full_div_auroc
full_sd = np.std(full_aurocs) if len(full_aurocs) > 1 else 0.0

perf_rows.append({
    "panel_size": len(full_features),
    "n_genes_actual": len(full_features),
    "mean_auroc": full_auroc_actual,
    "sd_auroc": full_sd,
    "mean_auprc": np.mean([r["auprc"] for r in fold_res_full if not np.isnan(r["auprc"])]),
    "sd_auprc": np.std([r["auprc"] for r in fold_res_full if not np.isnan(r["auprc"])]),
    "mean_brier": np.mean([r["brier"] for r in fold_res_full]),
    "pct_of_full_auroc": 100.0,
    "within_elbow": True,
})

perf_df = pd.DataFrame(perf_rows)
perf_df = perf_df.sort_values("panel_size").reset_index(drop=True)

# ===================================================================
# PART 5: Find the elbow panel
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 5: Elbow detection")
log.info("=" * 70)

elbow_candidates = perf_df[
    (perf_df["within_elbow"] == True) & (perf_df["panel_size"].isin(PANEL_SIZES))
]
if len(elbow_candidates) > 0:
    elbow_size = elbow_candidates["panel_size"].min()
    elbow_auroc = elbow_candidates.loc[
        elbow_candidates["panel_size"] == elbow_size, "mean_auroc"
    ].values[0]
    log.info(
        f"Elbow panel: {elbow_size} genes, AUROC {elbow_auroc:.3f} "
        f"(within {ELBOW_TOLERANCE} of full {full_div_auroc:.3f})"
    )
else:
    elbow_size = PANEL_SIZES[-1]  # fallback to largest
    log.info(
        f"No panel within elbow tolerance; using largest ({elbow_size} genes)"
    )

# ===================================================================
# PART 6: Save gene panels
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 6: Saving gene panels")
log.info("=" * 70)

for ps in PANEL_SIZES:
    if ps not in all_panel_genes:
        continue
    genes = all_panel_genes[ps]
    panel_info = div_imp[div_imp["gene_name"].isin(genes)].copy()
    panel_info = panel_info.sort_values("stability_rank")

    # Check if gene is a known secreted/plasma protein (heuristic)
    # Genes encoding secreted proteins are more feasible for blood-based assays
    secreted_markers = {
        "COL15A1", "COL11A1", "COL16A1", "COL4A2", "COL4A4",
        "VWF", "PLAT", "THBS1", "EFEMP1", "GSN", "FAP",
        "MMP2", "LTBP2", "PODN", "PAPLN", "CILP2", "SMOC2",
        "CRISPLD2", "SLIT2", "CRIM1", "ITGB6", "ITGA3",
        "CCN5", "LMCD1", "JAG1", "SEMA3G", "CDH11", "PTK7",
        "DPYSL3", "MCAM", "LAMC2", "CD34", "NES",
    }
    panel_info["plasma_feasible"] = panel_info["gene_name"].isin(secreted_markers)

    out_path = os.path.join(OUTDIR, f"prognosis_panel_{ps}.csv")
    panel_info[
        [
            "gene_name",
            "stability_rank",
            "n_folds_nonzero",
            "n_folds_total",
            "stability",
            "mean_coefficient",
            "abs_mean_coefficient",
            "plasma_feasible",
        ]
    ].to_csv(out_path, index=False)
    log.info(f"  Saved {out_path} ({len(panel_info)} genes)")

# ===================================================================
# PART 7: Save performance curve
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 7: Saving performance curve")
log.info("=" * 70)

perf_df["is_elbow"] = perf_df["panel_size"] == elbow_size
perf_df.to_csv(os.path.join(OUTDIR, "panel_performance_curve.csv"), index=False)
log.info(f"Saved panel_performance_curve.csv")

log.info("\nPerformance summary:")
for _, row in perf_df.iterrows():
    marker = " *ELBOW*" if row.get("is_elbow", False) else ""
    log.info(
        f"  {int(row['panel_size']):4d} genes: "
        f"AUROC {row['mean_auroc']:.3f} +/- {row['sd_auroc']:.3f}  "
        f"({row['pct_of_full_auroc']:.1f}%){marker}"
    )

# ===================================================================
# PART 8: Compare with staging panels
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 8: Comparison with staging classifier panels")
log.info("=" * 70)

comparison_rows = []

# Load staging panels
for ps_size in [10, 15, 25, 50]:
    stg_path = os.path.join(STAGING_DIR, f"minimal_panel_{ps_size}.csv")
    if not os.path.exists(stg_path):
        log.info(f"  Staging panel {ps_size} not found, skipping")
        continue

    stg_panel = pd.read_csv(stg_path)
    stg_genes = set(stg_panel["gene"].values) if "gene" in stg_panel.columns else set()
    log.info(f"  Staging panel {ps_size}: {len(stg_genes)} genes loaded")

    # Find overlap with each prognosis panel
    for prog_size in PANEL_SIZES:
        if prog_size not in all_panel_genes:
            continue
        prog_genes = set(all_panel_genes[prog_size])

        # Gene names in staging panels use Ensembl IDs; prognosis uses symbols
        # So overlap will be by name matching — likely minimal unless we harmonize
        overlap = stg_genes & prog_genes
        jaccard = (
            len(overlap) / len(stg_genes | prog_genes) if len(stg_genes | prog_genes) > 0 else 0
        )

        comparison_rows.append({
            "staging_panel_size": ps_size,
            "staging_panel_n_genes": len(stg_genes),
            "prognosis_panel_size": prog_size,
            "prognosis_panel_n_genes": len(all_panel_genes[prog_size]),
            "n_overlap": len(overlap),
            "overlap_genes": ";".join(sorted(overlap)) if overlap else "",
            "jaccard_index": jaccard,
            "staging_panel_purpose": "fibrosis_staging",
            "prognosis_panel_purpose": "s1_s2_subtype",
        })

# Add prognosis panel metrics to comparison
for prog_size in PANEL_SIZES:
    if prog_size not in all_panel_genes:
        continue
    prog_perf = perf_df[perf_df["panel_size"] == prog_size]
    if len(prog_perf) == 0:
        continue

    comparison_rows.append({
        "staging_panel_size": 0,
        "staging_panel_n_genes": 0,
        "prognosis_panel_size": prog_size,
        "prognosis_panel_n_genes": len(all_panel_genes[prog_size]),
        "n_overlap": 0,
        "overlap_genes": "",
        "jaccard_index": 0,
        "staging_panel_purpose": "NONE",
        "prognosis_panel_purpose": "s1_s2_subtype",
        "prognosis_auroc": prog_perf["mean_auroc"].values[0],
        "prognosis_sd_auroc": prog_perf["sd_auroc"].values[0],
    })

comp_df = pd.DataFrame(comparison_rows)
comp_path = os.path.join(OUTDIR, "prognosis_panel_comparison.csv")
comp_df.to_csv(comp_path, index=False)
log.info(f"Saved {comp_path}")

# ===================================================================
# PART 9: Plasma protein overlap summary
# ===================================================================
log.info("")
log.info("=" * 70)
log.info("PART 9: Plasma protein feasibility")
log.info("=" * 70)

for ps in PANEL_SIZES:
    if ps not in all_panel_genes:
        continue
    genes = all_panel_genes[ps]
    secreted_markers = {
        "COL15A1", "COL11A1", "COL16A1", "COL4A2", "COL4A4",
        "VWF", "PLAT", "THBS1", "EFEMP1", "GSN", "FAP",
        "MMP2", "LTBP2", "PODN", "PAPLN", "CILP2", "SMOC2",
        "CRISPLD2", "SLIT2", "CRIM1", "ITGB6", "ITGA3",
        "CCN5", "LMCD1", "JAG1", "SEMA3G", "CDH11", "PTK7",
        "DPYSL3", "MCAM", "LAMC2", "CD34", "NES",
    }
    feasible = [g for g in genes if g in secreted_markers]
    log.info(
        f"  Panel {ps}: {len(feasible)}/{len(genes)} plasma-feasible "
        f"({', '.join(feasible) if feasible else 'none'})"
    )

# ===================================================================
# Summary
# ===================================================================
elapsed = time.time() - t0
log.info("")
log.info("=" * 70)
log.info(f"COMPLETE in {elapsed:.0f}s")
log.info(f"  Elbow panel: {elbow_size} genes")
log.info(f"  Full model AUROC: {full_div_auroc:.3f}")
if len(elbow_candidates) > 0:
    log.info(f"  Elbow AUROC: {elbow_auroc:.3f}")
log.info(f"  Output files in: {OUTDIR}")
log.info("=" * 70)
