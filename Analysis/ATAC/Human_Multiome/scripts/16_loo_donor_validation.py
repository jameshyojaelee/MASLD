#!/usr/bin/env python3
"""
16_loo_donor_validation.py — Internal LOO + cell-subset jackknife validation
of GSE244832 scATAC regulon and chromVAR calls (B4-Fork-1)

Since GSE296875 was confirmed unsuitable for external scATAC MASLD replication
(no disease staging), this script provides the defensive replacement: internal
leave-one-donor-out (LOO) and cell-subset jackknife stability quantification
for the 18-donor GSE244832 multiome cohort.

DONOR-LEVEL TESTING (pseudoreplication fix, F079, 2026-05-30)
  All three regimes now test at the DONOR level, not the cell level. Cells are
  collapsed to donor means (utils_pseudobulk.aggregate_cells_to_donors) and the
  disease-vs-healthy contrast is a donor-level Mann-Whitney U (n <= 11 / group)
  via utils_pseudobulk.donor_groupwise_test. Previously each LOO/jackknife
  iteration ran a CELL-level Mann-Whitney over the ~15k-cell regulon activity
  matrix, so the reported "retention" was the stability of a cell-pooled p-value,
  not donor-level reproducibility. This mirrors 09_chromvar_propagation.py.

Three validation regimes:

  (1) LOO donor — Hepatocyte SCENIC+ disease-regulon test
      For each of 18 donors, drop the donor's hepatocytes from the precomputed
      regulon_activity_scores.csv, aggregate the REMAINING cells to per-donor
      mean activity, and re-run the donor-level MASLD vs Normal Mann-Whitney test
      on each of the 67 candidate regulons. Report:
        - retention rate of the 24 disease regulons (padj < 0.05, same sign)
        - per-donor concordance of the full disease-regulon set
        - the donor whose removal degrades the call set the most (outlier flag)

  (2) Cell-subset jackknife — Hepatocyte SCENIC+ regulons
      Drop a random 10%, 25%, 50% of hepatocyte cells (5 replicates per rate),
      aggregate the retained cells to per-donor means, and recompute the
      donor-level Mann-Whitney padj table. Report retention rate per drop level.

  (3) LOO donor — chromVAR top-20 TFs per cell type
      For Hepatocyte, Macrophage, Fibroblast (Stellate), Endothelial,
      Cholangiocyte, hold out one donor at a time, aggregate the remaining cells
      to per-donor mean deviations, re-rank TFs by absolute donor-mean
      difference (MASLD - Normal), and report the Spearman correlation of the
      held-out ranking vs the full-cohort top-20.

Inputs (all pre-existing, no GRN refit needed):
  - scenic_plus/regulon_activity_scores.csv          (14,998 hep cells × 67 regulons)
  - scenic_plus/disease_regulons.csv                  (24 baseline disease regulons)
  - results/label_transfer/snapatac2_label_transferred.h5ad
       (88,814 cells, obs: donor_id + condition + cell_type, backed mode)
  - results/chromvar_v2/chromvar_deviations.h5ad     (88,814 cells × 1019 motifs)

Outputs (Analysis/ATAC/Human_Multiome/results/loo_validation/):
  - loo_regulon_per_donor.csv         (18 donors × 24 regulons table)
  - loo_regulon_summary.csv           (per-donor n_retained, jaccard, etc.)
  - jackknife_regulon.csv             (3 drop rates × 5 reps × 24 regulons)
  - jackknife_regulon_summary.csv     (per drop rate aggregate stats)
  - loo_chromvar_per_celltype.csv     (cell-type × donor × top20-rho)
  - loo_chromvar_summary.csv          (per cell type aggregate stats)
  - loo_outlier_flag.csv              (donors with anomalous influence)
  - loo_validation_report.txt          (concise headline report)

Environment:
  micromamba activate snapatac2
  export PYTHONNOUSERSITE=1
  export HDF5_USE_FILE_LOCKING=FALSE   # see B2 lesson learned

Usage:
  python scripts/16_loo_donor_validation.py
  (or via run_16_loo_validation.sbatch)

Author: MASLD-Atlas pipeline (B4-Fork-1 of the ATAC improvement plan)
"""
from __future__ import annotations

import logging
import os
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from utils_pseudobulk import aggregate_cells_to_donors, donor_groupwise_test


# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
ROOT = Path(
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/"
    "Analysis/ATAC/Human_Multiome"
)
REG_ACT_PATH = ROOT / "scenic_plus" / "regulon_activity_scores.csv"
DISEASE_REG_PATH = ROOT / "scenic_plus" / "disease_regulons.csv"
LABEL_H5AD = (
    ROOT / "results" / "label_transfer" / "snapatac2_label_transferred.h5ad"
)
CHROMVAR_H5AD = (
    ROOT / "results" / "chromvar_v2" / "chromvar_deviations.h5ad"
)
OUT_DIR = ROOT / "results" / "loo_validation"

MASLD_LABELS = ("MASLD", "MASH", "MASL", "NASH", "NAFLD", "Steatosis")
NORMAL_LABELS = ("NORMAL", "Normal", "Healthy", "Control")

# chromVAR cell-type labels (note: chromvar h5ad uses plural / fibroblast)
CHROMVAR_CELLTYPES = (
    "Hepatocytes",
    "Macrophages",
    "Fibroblasts",
    "Endothelial_cells",
    "Cholangiocytes",
)

JACKKNIFE_RATES = (0.10, 0.25, 0.50)
JACKKNIFE_REPS = 5
JACKKNIFE_SEED = 42
TOP_K_TF = 20
ALPHA = 0.05

log = logging.getLogger("loo_validation")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _classify_condition(s: pd.Series) -> pd.Series:
    s2 = s.astype(str)
    out = pd.Series("Other", index=s2.index)
    out[s2.isin(MASLD_LABELS)] = "MASLD"
    out[s2.isin(NORMAL_LABELS)] = "Normal"
    return out


def _donor_group_map(donor_ids: np.ndarray, cond: np.ndarray) -> dict:
    """Build {donor -> 'disease'|'healthy'} from per-cell condition labels.

    `cond` is the per-cell _classify_condition output ('MASLD'|'Normal'|'Other').
    Condition is constant within a donor, so the per-donor label is the modal
    (here: first non-'Other') value. 'MASLD' -> 'disease', 'Normal' -> 'healthy';
    donors with only 'Other' cells are left unmapped (excluded from the test).
    """
    donor_ids = np.asarray(donor_ids)
    cond = np.asarray(cond)
    mapping: dict = {}
    for d in pd.unique(donor_ids):
        vals = cond[donor_ids == d]
        labels = pd.Series(vals[vals != "Other"])
        if labels.empty:
            continue
        top = labels.mode().iloc[0]
        if top == "MASLD":
            mapping[d] = "disease"
        elif top == "Normal":
            mapping[d] = "healthy"
    return mapping


def _test_all_regulons(
    activity: pd.DataFrame,
    donor_ids: np.ndarray,
    cond: np.ndarray,
) -> pd.DataFrame:
    """DONOR-level Mann-Whitney per regulon; BH-FDR over the 67 candidate regulons.

    Cells are collapsed to per-donor mean activity (the DONOR is the unit), then
    disease vs healthy is compared with utils_pseudobulk.donor_groupwise_test.
    Column schema preserved for downstream consumers (tf_name, regulon_activity_diff
    signed disease-minus-healthy, activity_pval, activity_padj) plus n_donors.
    """
    cols = list(activity.columns)
    donor_to_group = _donor_group_map(donor_ids, cond)
    donors, M, _ = aggregate_cells_to_donors(activity.values, donor_ids, agg="mean")
    res = donor_groupwise_test(
        M, donors, donor_to_group, "disease", "healthy", feature_names=cols
    )
    n_disease = int(res["n_donors_a"].iloc[0]) if len(res) else 0
    n_healthy = int(res["n_donors_b"].iloc[0]) if len(res) else 0
    return pd.DataFrame({
        "tf_name": res["feature"].astype(str).values,
        "mean_activity_masld": res["mean_disease"].values,
        "mean_activity_normal": res["mean_healthy"].values,
        # mean_diff = mean_disease - mean_healthy (positive = up in MASLD), matches legacy sign
        "regulon_activity_diff": res["mean_diff"].values,
        "activity_pval": res["pvalue"].values,
        "activity_padj": res["padj"].values,
        "n_donors_disease": n_disease,
        "n_donors_normal": n_healthy,
    })


# ---------------------------------------------------------------------------
# Step 1: Load + attach donor_id to regulon activity matrix
# ---------------------------------------------------------------------------
def load_regulon_activity_with_donor() -> tuple[pd.DataFrame, pd.Series, pd.Series, list[str]]:
    """Returns (activity_matrix, donor_id_per_cell, condition_per_cell, baseline_tfs)."""
    log.info("Loading hepatocyte regulon activity matrix...")
    activity = pd.read_csv(REG_ACT_PATH, index_col=0)
    log.info("  Activity matrix: %d cells x %d regulons", *activity.shape)

    log.info("Attaching donor_id + condition from label-transferred h5ad...")
    a = ad.read_h5ad(LABEL_H5AD, backed="r")
    obs = a.obs[["donor_id", "condition", "cell_type"]].copy()
    obs = obs.loc[obs.index.isin(activity.index)]
    log.info("  Matched %d / %d cells to h5ad", len(obs), len(activity))

    # Reindex activity to obs order
    obs = obs.reindex(activity.index.intersection(obs.index))
    activity = activity.loc[obs.index]

    donor_id = obs["donor_id"].astype(str)
    cond = _classify_condition(obs["condition"])
    log.info("  Donors: %d unique; condition counts: %s",
             donor_id.nunique(), dict(cond.value_counts()))

    log.info("Loading baseline disease regulons (24 expected)...")
    dis = pd.read_csv(DISEASE_REG_PATH)
    baseline_tfs = list(dis["tf_name"])
    log.info("  Baseline disease regulons: %d", len(baseline_tfs))

    try:
        a.file.close()
    except Exception:
        pass
    return activity, donor_id, cond, baseline_tfs


# ---------------------------------------------------------------------------
# Step 2: LOO donor — regulon stability
# ---------------------------------------------------------------------------
def run_loo_regulon(
    activity: pd.DataFrame,
    donor_id: pd.Series,
    cond: pd.Series,
    baseline_tfs: list[str],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    log.info("=" * 60)
    log.info("Regime 1: LOO donor regulon stability (n=%d donors)", donor_id.nunique())
    log.info("=" * 60)

    # Full-cohort reference test (for sign of effect) — DONOR-level
    donor_arr = donor_id.values
    cond_arr = cond.values
    full_df = _test_all_regulons(activity, donor_arr, cond_arr)
    full_df = full_df.set_index("tf_name")
    full_sig = set(full_df.index[(full_df["activity_padj"] < ALPHA)])
    log.info("  Full-cohort significant regulons: %d (baseline file: %d)",
             len(full_sig), len(baseline_tfs))

    donors = sorted(donor_id.unique())

    per_donor_records = []  # rows: donor × regulon
    summary_records = []     # rows: donor (aggregate)

    for d in donors:
        keep = (donor_id != d).values
        sub_act = activity.iloc[keep]
        sub_cond = cond.iloc[keep]
        sub_donor = donor_arr[keep]
        n_drop = (~keep).sum()
        # cell-level condition counts retained as diagnostics only
        n_m = int((sub_cond == "MASLD").sum())
        n_n = int((sub_cond == "Normal").sum())

        # DONOR-level test on the remaining cells (aggregate cells -> donors inside)
        df = _test_all_regulons(sub_act, sub_donor, sub_cond.values).set_index("tf_name")
        n_donor_m = int(df["n_donors_disease"].iloc[0]) if len(df) else 0
        n_donor_n = int(df["n_donors_normal"].iloc[0]) if len(df) else 0

        # Compare to baseline disease regulons
        for tf in baseline_tfs:
            if tf in df.index:
                row = df.loc[tf]
                base = full_df.loc[tf] if tf in full_df.index else None
                same_sign = (
                    base is not None
                    and np.sign(row["regulon_activity_diff"]) == np.sign(base["regulon_activity_diff"])
                )
                retained = (row["activity_padj"] < ALPHA) and same_sign
                per_donor_records.append({
                    "donor": d, "tf_name": tf,
                    "activity_padj_loo": row["activity_padj"],
                    "activity_pval_loo": row["activity_pval"],
                    "diff_loo": row["regulon_activity_diff"],
                    "same_sign_as_full": same_sign,
                    "retained_padj_lt_alpha": retained,
                })

        # Per-donor summary
        loo_sig = set(df.index[df["activity_padj"] < ALPHA])
        # Among the baseline 24
        retained_baseline = sum(
            1 for tf in baseline_tfs
            if tf in df.index
            and df.loc[tf, "activity_padj"] < ALPHA
            and tf in full_df.index
            and np.sign(df.loc[tf, "regulon_activity_diff"]) == np.sign(full_df.loc[tf, "regulon_activity_diff"])
        )
        # Jaccard vs full significant set
        union = full_sig | loo_sig
        jacc = len(full_sig & loo_sig) / len(union) if union else np.nan
        summary_records.append({
            "donor": d,
            "n_cells_dropped": int(n_drop),
            "n_cells_used": int(keep.sum()),
            "n_masld_cells": int(n_m),
            "n_normal_cells": int(n_n),
            "n_donors": int(n_donor_m + n_donor_n),
            "n_donors_disease": int(n_donor_m),
            "n_donors_normal": int(n_donor_n),
            "n_baseline_retained": int(retained_baseline),
            "frac_baseline_retained": retained_baseline / len(baseline_tfs) if baseline_tfs else np.nan,
            "n_loo_sig": len(loo_sig),
            "jaccard_full_vs_loo": jacc,
        })
        log.info(
            "  Drop %s: %d cells → %d donors (%d disease, %d healthy) → retained %d/%d (%.1f%%); jaccard=%.3f",
            d, n_drop, n_donor_m + n_donor_n, n_donor_m, n_donor_n,
            retained_baseline, len(baseline_tfs),
            100 * retained_baseline / len(baseline_tfs) if baseline_tfs else 0.0,
            jacc,
        )

    per_df = pd.DataFrame(per_donor_records)
    sum_df = pd.DataFrame(summary_records)
    return per_df, sum_df


# ---------------------------------------------------------------------------
# Step 3: Cell-subset jackknife
# ---------------------------------------------------------------------------
def run_jackknife_regulon(
    activity: pd.DataFrame,
    donor_id: pd.Series,
    cond: pd.Series,
    baseline_tfs: list[str],
    full_df: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    log.info("=" * 60)
    log.info("Regime 2: Cell-subset jackknife (rates=%s, reps=%d)",
             JACKKNIFE_RATES, JACKKNIFE_REPS)
    log.info("=" * 60)

    rng = np.random.RandomState(JACKKNIFE_SEED)
    n_cells = len(activity)
    donor_arr = donor_id.values
    full_sig = set(full_df.index[(full_df["activity_padj"] < ALPHA)])

    detail_records = []
    summary_records = []

    for rate in JACKKNIFE_RATES:
        retentions = []
        jaccards = []
        n_donors_seen = []
        for rep in range(JACKKNIFE_REPS):
            keep_n = int(round(n_cells * (1 - rate)))
            keep_idx = rng.choice(n_cells, keep_n, replace=False)
            sub_act = activity.iloc[keep_idx]
            sub_cond = cond.iloc[keep_idx]
            sub_donor = donor_arr[keep_idx]
            # DONOR-level test: aggregate the retained cells to per-donor means
            df = _test_all_regulons(sub_act, sub_donor, sub_cond.values).set_index("tf_name")
            n_donor = (
                int(df["n_donors_disease"].iloc[0] + df["n_donors_normal"].iloc[0])
                if len(df) else 0
            )
            n_donors_seen.append(n_donor)

            for tf in baseline_tfs:
                if tf in df.index:
                    row = df.loc[tf]
                    base = full_df.loc[tf] if tf in full_df.index else None
                    same_sign = (
                        base is not None
                        and np.sign(row["regulon_activity_diff"]) == np.sign(base["regulon_activity_diff"])
                    )
                    retained = (row["activity_padj"] < ALPHA) and same_sign
                    detail_records.append({
                        "drop_rate": rate, "rep": rep, "tf_name": tf,
                        "activity_padj_jk": row["activity_padj"],
                        "diff_jk": row["regulon_activity_diff"],
                        "retained_padj_lt_alpha": retained,
                        "n_donors": n_donor,
                    })
            n_ret = sum(
                1 for tf in baseline_tfs
                if tf in df.index
                and df.loc[tf, "activity_padj"] < ALPHA
                and tf in full_df.index
                and np.sign(df.loc[tf, "regulon_activity_diff"]) == np.sign(full_df.loc[tf, "regulon_activity_diff"])
            )
            retentions.append(n_ret)
            loo_sig = set(df.index[df["activity_padj"] < ALPHA])
            union = full_sig | loo_sig
            jacc = len(full_sig & loo_sig) / len(union) if union else np.nan
            jaccards.append(jacc)
            log.info("  Drop %.0f%% rep %d: kept %d/%d (n_donors=%d), jaccard=%.3f",
                     100 * rate, rep, n_ret, len(baseline_tfs), n_donor, jacc)

        summary_records.append({
            "drop_rate": rate,
            "n_reps": JACKKNIFE_REPS,
            "n_donors": int(np.max(n_donors_seen)) if n_donors_seen else 0,
            "mean_n_retained": float(np.mean(retentions)),
            "sd_n_retained": float(np.std(retentions, ddof=1)) if len(retentions) > 1 else 0.0,
            "mean_frac_retained": float(np.mean(retentions) / len(baseline_tfs)) if baseline_tfs else np.nan,
            "mean_jaccard_full": float(np.mean(jaccards)),
            "sd_jaccard_full": float(np.std(jaccards, ddof=1)) if len(jaccards) > 1 else 0.0,
        })

    return pd.DataFrame(detail_records), pd.DataFrame(summary_records)


# ---------------------------------------------------------------------------
# Step 4: LOO donor — chromVAR top-20 TFs per cell type
# ---------------------------------------------------------------------------
def run_loo_chromvar() -> tuple[pd.DataFrame, pd.DataFrame]:
    log.info("=" * 60)
    log.info("Regime 3: LOO chromVAR top-20 TFs per cell type (n=%d CTs)",
             len(CHROMVAR_CELLTYPES))
    log.info("=" * 60)

    log.info("Loading chromvar_deviations.h5ad (memory-backed)...")
    a = ad.read_h5ad(CHROMVAR_H5AD)
    log.info("  Loaded chromVAR: %d cells x %d motifs", a.n_obs, a.n_vars)
    if not a.var_names.is_unique:
        a.var_names_make_unique()
        log.info("  De-duplicated motif names")

    donors = sorted(a.obs["donor_id"].astype(str).unique())
    motifs = list(a.var_names)
    X = a.X.astype(np.float32) if hasattr(a.X, "astype") else np.asarray(a.X, dtype=np.float32)

    per_records = []
    summary_records = []

    for ct in CHROMVAR_CELLTYPES:
        ct_mask = (a.obs["cell_type"].astype(str) == ct).values
        if ct_mask.sum() < 200:
            log.warning("  Skipping %s (only %d cells)", ct, ct_mask.sum())
            continue
        log.info("Cell type %s: %d cells", ct, ct_mask.sum())

        cond_ct = _classify_condition(a.obs.loc[ct_mask, "condition"]).values
        donor_ct = a.obs.loc[ct_mask, "donor_id"].astype(str).values
        X_ct = X[ct_mask]

        # Donor -> disease/healthy map for this cell type's cells
        donor_to_group = _donor_group_map(donor_ct, cond_ct)

        def _donor_mean_diff(cell_mask: np.ndarray):
            """Aggregate the masked cells to per-donor means, then return
            (mean over disease donors - mean over healthy donors) per motif,
            plus the donor counts used. The DONOR is the unit (pseudorep fix)."""
            ds, M, _ = aggregate_cells_to_donors(
                X_ct[cell_mask], donor_ct[cell_mask], agg="mean"
            )
            grp = np.array([donor_to_group.get(x, None) for x in ds], dtype=object)
            a_idx = np.where(grp == "disease")[0]
            b_idx = np.where(grp == "healthy")[0]
            if a_idx.size < 2 or b_idx.size < 2:
                return None, a_idx.size, b_idx.size
            diff = M[a_idx, :].mean(axis=0) - M[b_idx, :].mean(axis=0)
            return np.asarray(diff).ravel(), a_idx.size, b_idx.size

        # Full ranking (donor-level mean difference)
        full_diff, n_dis_full, n_heal_full = _donor_mean_diff(
            np.ones(X_ct.shape[0], dtype=bool)
        )
        if full_diff is None:
            log.warning("  Skipping %s (insufficient donors: disease=%d / healthy=%d)",
                        ct, n_dis_full, n_heal_full)
            continue
        log.info("  %s donor counts: %d disease, %d healthy", ct, n_dis_full, n_heal_full)
        full_abs = np.abs(full_diff)
        full_rank = pd.Series(full_abs, index=motifs).sort_values(ascending=False)
        top20 = list(full_rank.head(TOP_K_TF).index)
        log.info("  Full top-20 TFs (by |Δdev|, donor-level): %s", ", ".join(top20[:5]) + ", ...")

        # LOO per donor
        rhos = []
        kendalls = []  # overlap@20
        for d in donors:
            keep = (donor_ct != d)
            sub_diff, n_dis, n_heal = _donor_mean_diff(keep)
            if sub_diff is None:
                continue
            sub_abs = np.abs(sub_diff)
            sub_rank = pd.Series(sub_abs, index=motifs).sort_values(ascending=False)
            sub_top20 = list(sub_rank.head(TOP_K_TF).index)

            # Spearman of rank vectors on the full top-20
            full_ranks_on_top = full_rank.rank(ascending=False)[top20].values
            sub_ranks_on_top = sub_rank.rank(ascending=False)[top20].values
            if np.std(full_ranks_on_top) == 0 or np.std(sub_ranks_on_top) == 0:
                rho = np.nan
            else:
                rho, _ = spearmanr(full_ranks_on_top, sub_ranks_on_top)
            overlap20 = len(set(sub_top20) & set(top20))

            rhos.append(rho)
            kendalls.append(overlap20)
            per_records.append({
                "cell_type": ct, "donor": d,
                "n_cells_dropped": int((~keep).sum()),
                "n_donors": int(n_dis + n_heal),
                "n_donors_disease": int(n_dis),
                "n_donors_normal": int(n_heal),
                "spearman_rho_top20": rho,
                "overlap_top20": overlap20,
                "top20_held_out": ";".join(sub_top20[:10]),
            })

        rhos_arr = np.array([r for r in rhos if not np.isnan(r)])
        ks_arr = np.array(kendalls)
        summary_records.append({
            "cell_type": ct,
            "n_donors_evaluated": len(rhos),
            "n_donors": int(n_dis_full + n_heal_full),
            "mean_spearman_top20": float(np.mean(rhos_arr)) if len(rhos_arr) else np.nan,
            "sd_spearman_top20": float(np.std(rhos_arr, ddof=1)) if len(rhos_arr) > 1 else 0.0,
            "min_spearman_top20": float(np.min(rhos_arr)) if len(rhos_arr) else np.nan,
            "mean_overlap_top20": float(np.mean(ks_arr)) if len(ks_arr) else np.nan,
            "min_overlap_top20": float(np.min(ks_arr)) if len(ks_arr) else 0,
            "full_top20_TFs": ";".join(top20),
        })
        log.info("  %s: mean ρ=%.3f, mean overlap@20=%.1f/20",
                 ct, np.mean(rhos_arr) if len(rhos_arr) else np.nan,
                 np.mean(ks_arr) if len(ks_arr) else np.nan)

    return pd.DataFrame(per_records), pd.DataFrame(summary_records)


# ---------------------------------------------------------------------------
# Step 5: Outlier flag
# ---------------------------------------------------------------------------
def derive_outlier_flag(
    regulon_summary: pd.DataFrame,
    chromvar_per: pd.DataFrame,
) -> pd.DataFrame:
    """Flag donors whose removal degrades the call set significantly.

    Heuristic: a donor is an outlier if
      - its frac_baseline_retained is below mean - 1.5×SD across donors, OR
      - its mean Spearman ρ across cell types is below mean - 1.5×SD.
    """
    log.info("=" * 60)
    log.info("Outlier flagging")
    log.info("=" * 60)

    flags = regulon_summary[["donor", "frac_baseline_retained", "jaccard_full_vs_loo"]].copy()
    mu_r = flags["frac_baseline_retained"].mean()
    sd_r = flags["frac_baseline_retained"].std(ddof=1)
    flags["regulon_outlier"] = (
        flags["frac_baseline_retained"] < (mu_r - 1.5 * sd_r)
    )

    if not chromvar_per.empty:
        cv = (
            chromvar_per.groupby("donor")["spearman_rho_top20"]
            .mean()
            .reset_index()
            .rename(columns={"spearman_rho_top20": "mean_chromvar_rho"})
        )
        flags = flags.merge(cv, on="donor", how="left")
        mu_c = flags["mean_chromvar_rho"].mean()
        sd_c = flags["mean_chromvar_rho"].std(ddof=1)
        flags["chromvar_outlier"] = (
            flags["mean_chromvar_rho"] < (mu_c - 1.5 * sd_c)
        )
    else:
        flags["mean_chromvar_rho"] = np.nan
        flags["chromvar_outlier"] = False

    flags["any_outlier"] = flags["regulon_outlier"] | flags["chromvar_outlier"]
    n_out = flags["any_outlier"].sum()
    log.info("  %d / %d donors flagged as potential outliers", n_out, len(flags))
    if n_out:
        log.info("    %s", ", ".join(flags.loc[flags["any_outlier"], "donor"].tolist()))
    return flags


# ---------------------------------------------------------------------------
# Step 6: Report
# ---------------------------------------------------------------------------
def write_report(
    out_dir: Path,
    baseline_n: int,
    regulon_summary: pd.DataFrame,
    jackknife_summary: pd.DataFrame,
    chromvar_summary: pd.DataFrame,
    outlier_flag: pd.DataFrame,
):
    lines = []
    lines.append("=" * 70)
    lines.append("B4-Fork-1 INTERNAL VALIDATION REPORT — GSE244832 scATAC (18 donors)")
    lines.append("=" * 70)
    lines.append("")
    lines.append("REGIME 1 — LOO donor regulon stability (DONOR-level Mann-Whitney)")
    lines.append("-" * 70)
    lines.append(f"  Baseline disease regulons (full cohort): {baseline_n}")
    if "n_donors" in regulon_summary.columns and len(regulon_summary):
        lines.append(f"  Donors per LOO test (after holdout): "
                     f"{int(regulon_summary['n_donors'].min())}-"
                     f"{int(regulon_summary['n_donors'].max())} "
                     f"({int(regulon_summary['n_donors_disease'].max())} disease / "
                     f"{int(regulon_summary['n_donors_normal'].max())} healthy at full)")
    mu = regulon_summary["frac_baseline_retained"].mean()
    sd = regulon_summary["frac_baseline_retained"].std(ddof=1)
    lines.append(f"  Mean retention across 18 LOO iterations: "
                 f"{mu*100:.1f}% ± {sd*100:.1f}% "
                 f"(min {regulon_summary['frac_baseline_retained'].min()*100:.1f}%, "
                 f"max {regulon_summary['frac_baseline_retained'].max()*100:.1f}%)")
    lines.append(f"  Mean Jaccard(LOO ∩ Full): "
                 f"{regulon_summary['jaccard_full_vs_loo'].mean():.3f} "
                 f"± {regulon_summary['jaccard_full_vs_loo'].std(ddof=1):.3f}")
    lines.append("")

    lines.append("REGIME 2 — Cell-subset jackknife")
    lines.append("-" * 70)
    for _, r in jackknife_summary.iterrows():
        lines.append(
            f"  Drop {r['drop_rate']*100:.0f}% × {int(r['n_reps'])} reps: "
            f"retention {r['mean_frac_retained']*100:.1f}% ± {r['sd_n_retained']/baseline_n*100:.1f}%; "
            f"jaccard {r['mean_jaccard_full']:.3f} ± {r['sd_jaccard_full']:.3f}"
        )
    lines.append("")

    lines.append("REGIME 3 — LOO chromVAR top-20 TF rank stability")
    lines.append("-" * 70)
    if not chromvar_summary.empty:
        for _, r in chromvar_summary.iterrows():
            lines.append(
                f"  {r['cell_type']:>20s}: ρ = {r['mean_spearman_top20']:.3f} "
                f"± {r['sd_spearman_top20']:.3f} (min {r['min_spearman_top20']:.3f}); "
                f"overlap@20 = {r['mean_overlap_top20']:.1f}/20 "
                f"(min {int(r['min_overlap_top20'])}/20) "
                f"across {int(r['n_donors_evaluated'])} donors"
            )
    else:
        lines.append("  (no chromVAR results)")
    lines.append("")

    lines.append("OUTLIER FLAGS")
    lines.append("-" * 70)
    out_donors = outlier_flag.loc[outlier_flag["any_outlier"], "donor"].tolist()
    if out_donors:
        lines.append(f"  Donors flagged (≥1.5 SD below mean on regulon or chromVAR): "
                     f"{', '.join(out_donors)}")
    else:
        lines.append("  No donors flagged (no anomalous influence detected).")
    lines.append("")
    lines.append("=" * 70)

    path = out_dir / "loo_validation_report.txt"
    path.write_text("\n".join(lines))
    log.info("Report written: %s", path)
    print("\n".join(lines))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    # --- Step 1: load activity + metadata ---
    activity, donor_id, cond, baseline_tfs = load_regulon_activity_with_donor()

    # --- Step 2: LOO regulon ---
    per_donor_df, regulon_summary = run_loo_regulon(activity, donor_id, cond, baseline_tfs)
    per_donor_df.to_csv(OUT_DIR / "loo_regulon_per_donor.csv", index=False)
    regulon_summary.to_csv(OUT_DIR / "loo_regulon_summary.csv", index=False)
    log.info("Saved: loo_regulon_per_donor.csv (%d rows)", len(per_donor_df))
    log.info("Saved: loo_regulon_summary.csv (%d rows)", len(regulon_summary))

    # Cache full reference regulon test for jackknife (avoid recompute) — DONOR-level
    full_df = _test_all_regulons(activity, donor_id.values, cond.values).set_index("tf_name")

    # --- Step 3: jackknife ---
    jk_detail, jk_summary = run_jackknife_regulon(activity, donor_id, cond, baseline_tfs, full_df)
    jk_detail.to_csv(OUT_DIR / "jackknife_regulon.csv", index=False)
    jk_summary.to_csv(OUT_DIR / "jackknife_regulon_summary.csv", index=False)
    log.info("Saved: jackknife_regulon.csv (%d rows)", len(jk_detail))
    log.info("Saved: jackknife_regulon_summary.csv (%d rows)", len(jk_summary))

    # --- Step 4: LOO chromVAR ---
    chromvar_per, chromvar_summary = run_loo_chromvar()
    chromvar_per.to_csv(OUT_DIR / "loo_chromvar_per_celltype.csv", index=False)
    chromvar_summary.to_csv(OUT_DIR / "loo_chromvar_summary.csv", index=False)
    log.info("Saved: loo_chromvar_per_celltype.csv (%d rows)", len(chromvar_per))
    log.info("Saved: loo_chromvar_summary.csv (%d rows)", len(chromvar_summary))

    # --- Step 5: outlier flag ---
    outlier_flag = derive_outlier_flag(regulon_summary, chromvar_per)
    outlier_flag.to_csv(OUT_DIR / "loo_outlier_flag.csv", index=False)
    log.info("Saved: loo_outlier_flag.csv (%d rows)", len(outlier_flag))

    # --- Step 6: report ---
    write_report(
        OUT_DIR,
        baseline_n=len(baseline_tfs),
        regulon_summary=regulon_summary,
        jackknife_summary=jk_summary,
        chromvar_summary=chromvar_summary,
        outlier_flag=outlier_flag,
    )

    log.info("=" * 60)
    log.info("B4-Fork-1 LOO/jackknife validation COMPLETE in %.1f minutes",
             (time.time() - t0) / 60)
    log.info("Outputs: %s", OUT_DIR)
    log.info("=" * 60)


if __name__ == "__main__":
    main()
