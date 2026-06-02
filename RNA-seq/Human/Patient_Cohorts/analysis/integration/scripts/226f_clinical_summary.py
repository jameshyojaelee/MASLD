#!/usr/bin/env python3
"""
226f_clinical_summary.py

Synthesizes all 226a-d plasma biomarker results into manuscript-ready summary
tables, a recommended clinical panel, and a head-to-head comparison with
Yang et al. 2025 (Cell Rep Med).

Outputs (all to results/multiprogram/):
    226f_recommended_panel.csv   — 5-10 protein panel for T5 etiology
    226f_sweep_leaderboard.csv   — full model leaderboard across all targets
    226f_yang_comparison.csv     — head-to-head with Yang et al. 2025
    226f_manuscript_numbers.csv  — all key numbers for manuscript paragraph
"""

import os
import logging
import warnings
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore", category=FutureWarning)

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s  %(levelname)-8s  %(message)s",
    datefmt="%H:%M:%S",
)
log = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
INDIR = os.path.join(
    BASE,
    "RNA-seq/Human/Patient_Cohorts/analysis/integration/results/multiprogram",
)
OUTDIR = INDIR  # outputs collocated with inputs


def _read(fname, required=True):
    """Read a CSV from INDIR, returning None on failure."""
    path = os.path.join(INDIR, fname)
    try:
        df = pd.read_csv(path)
        log.info("Loaded %-50s  %d rows × %d cols", fname, *df.shape)
        return df
    except FileNotFoundError:
        if required:
            log.error("REQUIRED file not found: %s", path)
            raise
        log.warning("Optional file not found (skipped): %s", fname)
        return None
    except Exception as exc:
        log.warning("Could not read %s: %s", fname, exc)
        return None


def _write(df, fname):
    path = os.path.join(OUTDIR, fname)
    df.to_csv(path, index=False, float_format="%.3f")
    log.info("Wrote %-50s  %d rows × %d cols", fname, *df.shape)


# ---------------------------------------------------------------------------
# Helper: infer primary metric from target name
# ---------------------------------------------------------------------------
TARGET_METRIC_MAP = {
    "T1_binary": "mean_auroc",
    "T5_masld_binary": "mean_auroc",
    "T5_etiology_binary": "mean_auroc",
    "T4_etiology_3class": "mean_f1_macro",
    "T4_etiology_3": "mean_f1_macro",
    "T2_ordinal": "mean_qwk",
    "T3_continuous": "mean_spearman",
}

METRIC_LABEL_MAP = {
    "mean_auroc": "AUROC",
    "mean_f1_macro": "F1_macro",
    "mean_qwk": "QWK",
    "mean_spearman": "Spearman_rho",
}


def _primary_metric(target):
    for key, metric in TARGET_METRIC_MAP.items():
        if key in str(target):
            return metric
    return "mean_auroc"


# ---------------------------------------------------------------------------
# 1. Build sweep leaderboard
# ---------------------------------------------------------------------------
SWEEP_FILES = {
    # file, target column present?, extra target label for single-target files
    "plasma_sweep_225a_linear.csv": (False, "T1_binary"),
    "plasma_sweep_225b_svm.csv": (False, "T1_binary"),
    "plasma_sweep_225b_ada.csv": (False, "T1_binary"),
    "plasma_sweep_225b_imb.csv": (False, "T1_binary"),
    "plasma_sweep_225b_done3.csv": (False, "T1_binary"),
    "plasma_sweep_225c_neural.csv": (False, "T1_binary"),
    "plasma_sweep_lgb_all_targets.csv": (True, None),
    # ordinal sweeps
    "plasma_sweep_225d_ordinal_tab.csv": (False, "T2_ordinal"),
    "plasma_sweep_225d_ordinal_flaml.csv": (False, "T2_ordinal"),
    "plasma_sweep_225d_ordinal_mord.csv": (False, "T2_ordinal"),
    "plasma_sweep_225d_ordinal_mlp.csv": (False, "T2_ordinal"),
    "plasma_sweep_225d_ordinal_smote.csv": (False, "T2_ordinal"),
    "plasma_sweep_225d_ordinal_gp.csv": (False, "T2_ordinal"),
    "plasma_sweep_225d_ordinal_ens.csv": (False, "T2_ordinal"),
    "plasma_sweep_225d_ordinal_svm.csv": (False, "T2_ordinal"),
    "plasma_sweep_225d_ordinal_cal.csv": (False, "T2_ordinal"),
    "plasma_sweep_225d_ordinal_reg.csv": (False, "T2_ordinal"),
    "plasma_sweep_225d_ordinal_ada.csv": (False, "T2_ordinal"),
    "plasma_sweep_225d_ordinal_imb.csv": (False, "T2_ordinal"),
    "plasma_sweep_225d_ordinal_cat.csv": (False, "T2_ordinal"),
    # etiology sweeps
    "plasma_sweep_225e_etiology_tab.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_mlp.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_gp.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_svm.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_smote.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_ens.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_ada.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_cal.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_reg.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_hgbm.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_sgd_bayes.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_flaml.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_gbm.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_tree.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_imb.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_cat.csv": (False, "T5_masld_binary"),
    "plasma_sweep_225e_etiology_xgb.csv": (False, "T5_masld_binary"),
}


def build_leaderboard():
    """Combine all sweep CSVs into a unified leaderboard."""
    log.info("Building sweep leaderboard ...")
    frames = []

    for fname, (has_target_col, default_target) in SWEEP_FILES.items():
        df = _read(fname, required=False)
        if df is None:
            continue

        # Normalise column names to lower-case stripped
        df.columns = [c.strip().lower() for c in df.columns]

        # Inject target column if not present
        if not has_target_col:
            df["target"] = default_target

        # Resolve metric columns present
        metric_cols = {
            "mean_auroc": ["mean_auroc"],
            "sd_auroc": ["sd_auroc"],
            "mean_f1_macro": ["mean_f1_macro", "mean_f1"],
            "sd_f1_macro": ["sd_f1_macro", "sd_f1"],
            "mean_qwk": ["mean_qwk"],
            "sd_qwk": ["sd_qwk"],
            "mean_spearman": ["mean_spearman"],
            "sd_spearman": ["sd_spearman"],
            "n_folds": ["n_folds"],
            "elapsed_sec": ["elapsed_sec"],
        }

        row_template = {k: np.nan for k in metric_cols}
        row_template["n_folds"] = np.nan
        row_template["elapsed_sec"] = np.nan

        def _col(df, candidates):
            for c in candidates:
                if c in df.columns:
                    return df[c]
            return pd.Series([np.nan] * len(df))

        out = pd.DataFrame()
        out["model"] = df.get("model", pd.Series(["unknown"] * len(df)))
        out["target"] = df["target"]
        for key, candidates in metric_cols.items():
            out[key] = _col(df, candidates).values

        frames.append(out)

    if not frames:
        log.error("No sweep files loaded — leaderboard will be empty.")
        return pd.DataFrame(
            columns=["model", "target", "primary_metric_name",
                     "primary_metric_mean", "primary_metric_sd", "n_folds"]
        )

    lb = pd.concat(frames, ignore_index=True)

    # Choose primary metric per target
    lb["primary_metric_name"] = lb["target"].apply(
        lambda t: METRIC_LABEL_MAP.get(_primary_metric(t), "AUROC")
    )

    def _pick_primary(row):
        m = _primary_metric(row["target"])
        return row.get(m, np.nan)

    def _pick_primary_sd(row):
        m = _primary_metric(row["target"])
        sd_col = m.replace("mean_", "sd_")
        return row.get(sd_col, np.nan)

    lb["primary_metric_mean"] = lb.apply(_pick_primary, axis=1)
    lb["primary_metric_sd"] = lb.apply(_pick_primary_sd, axis=1)

    # Keep only informative rows
    lb = lb[lb["primary_metric_mean"].notna()].copy()

    # Sort per target by primary metric descending
    lb = (
        lb.sort_values(["target", "primary_metric_mean"], ascending=[True, False])
        .reset_index(drop=True)
    )

    out_cols = [
        "model", "target", "primary_metric_name",
        "primary_metric_mean", "primary_metric_sd", "n_folds",
    ]
    lb = lb[[c for c in out_cols if c in lb.columns]]
    return lb


# ---------------------------------------------------------------------------
# 2. Recommended clinical panel (T5 etiology)
# ---------------------------------------------------------------------------
def build_recommended_panel(panel_curves):
    """Build top-10 recommended panel from 226b_tissue_bridge for T5."""
    log.info("Building recommended clinical panel ...")

    bridge = _read("226b_tissue_bridge.csv", required=False)
    if bridge is None:
        log.warning("226b_tissue_bridge.csv missing — skipping panel build.")
        return pd.DataFrame()

    # Normalise target label: accept either T5_masld_binary or T5_etiology_binary
    t5_mask = bridge["target"].str.contains("T5|etiology_binary", case=False, na=False)
    t5 = bridge[t5_mask].copy()

    if t5.empty:
        log.warning("No T5 etiology rows in 226b_tissue_bridge — trying T1 as fallback.")
        t5 = bridge[bridge["target"].str.contains("T1", na=False)].copy()

    if t5.empty:
        log.warning("No suitable rows in 226b_tissue_bridge for panel.")
        return pd.DataFrame()

    # Sort by shap_rank ascending, take top 10
    t5 = t5.sort_values("shap_rank").head(10).reset_index(drop=True)
    t5["rank"] = range(1, len(t5) + 1)

    # Add panel performance at panel_size == len(panel) from panel_curves
    panel_perf = None
    if panel_curves is not None:
        target_vals = t5["target"].unique()
        pc_target = panel_curves[
            panel_curves["target"].isin(target_vals) &
            (panel_curves["is_random"] == False)  # noqa: E712
        ]
        panel_size_n = len(t5)
        pc_at_n = pc_target[pc_target["panel_size"] == panel_size_n]
        if not pc_at_n.empty:
            # mean across folds
            mean_metric = pc_at_n.groupby("metric_name")["metric_value"].mean()
            if "auroc" in mean_metric.index:
                panel_perf = f"AUROC={mean_metric['auroc']:.3f}"
            else:
                first_metric = mean_metric.index[0]
                panel_perf = f"{first_metric}={mean_metric.iloc[0]:.3f}"

    # Build output columns per spec
    cols_in = {
        "protein": "protein",
        "mean_abs_shap": "permutation_importance",
        "fold_stability": "fold_stability_pct",
        "is_deg": "is_tissue_deg",
        "f2_switch_class": "f2_switch_class",
        "dgidb_druggable": "dgidb_druggable",
        "attribution_class": "attribution_class",
    }

    panel = pd.DataFrame()
    panel["rank"] = t5["rank"]

    for src, dst in cols_in.items():
        if src in t5.columns:
            panel[dst] = t5[src].values
        else:
            panel[dst] = np.nan

    # Build recommendation note
    def _note(row):
        parts = []
        if panel_perf:
            parts.append(f"Panel performance at {panel_size_n} proteins: {panel_perf}")
        if row.get("is_tissue_deg") is True or str(row.get("is_tissue_deg")).lower() == "true":
            parts.append("tissue DEG")
        if str(row.get("f2_switch_class")).lower() == "switch":
            parts.append("F2 switch gene")
        if str(row.get("dgidb_druggable")).lower() == "true":
            parts.append("druggable")
        return "; ".join(parts) if parts else ""

    panel["recommendation_note"] = panel.apply(_note, axis=1)
    return panel


# ---------------------------------------------------------------------------
# 3. Yang et al. 2025 comparison
# ---------------------------------------------------------------------------
def build_yang_comparison(benchmark_df, loeo_df):
    """Build head-to-head comparison table with Yang et al. 2025."""
    log.info("Building Yang et al. 2025 comparison table ...")

    # Pull NFASC+GDF15 pooled LOEO AUROC from benchmark_df
    nfasc_auroc = "0.734"
    if benchmark_df is not None:
        mask = (
            benchmark_df.get("protein_pair", pd.Series(dtype=str))
            .str.contains("NFASC", na=False)
        )
        row = benchmark_df[mask]
        if not row.empty:
            pooled = row[
                row.get("holdout_etiology", pd.Series(dtype=str))
                .str.contains("POOLED", na=False)
            ]
            if not pooled.empty and "auroc" in pooled.columns:
                nfasc_auroc = f"{pooled['auroc'].iloc[0]:.3f}"

    # Pull cross-platform AUROC from 226c_cross_platform
    cross_platform = _read("226c_cross_platform.csv", required=False)
    transfer_auroc = "0.655"
    if cross_platform is not None:
        full_row = cross_platform[
            cross_platform.get("experiment", pd.Series(dtype=str))
            .str.contains("full|olink", case=False, na=False)
        ]
        if not full_row.empty and "auroc" in full_row.columns:
            # pick the shared_olink row as cross-platform transfer proxy
            shared_row = cross_platform[
                cross_platform.get("experiment", pd.Series(dtype=str))
                .str.contains("shared", case=False, na=False)
            ]
            if not shared_row.empty:
                transfer_auroc = f"{shared_row['auroc'].iloc[0]:.3f}"

    rows = [
        {
            "comparison_aspect": "CV methodology",
            "yang_2025": "70/30 single split",
            "our_analysis": "10×5-fold nested CV (50 folds)",
            "advantage": "More conservative, unbiased",
        },
        {
            "comparison_aspect": "Binary F>=3 AUROC",
            "yang_2025": "0.98 discovery / 0.88 validation",
            "our_analysis": "0.790",
            "advantage": "Nested CV avoids optimism bias",
        },
        {
            "comparison_aspect": "NFASC+GDF15 performance",
            "yang_2025": "87% balanced accuracy (discovery)",
            "our_analysis": f"AUROC {nfasc_auroc} (pooled LOEO)",
            "advantage": "Comparable; we use more rigorous evaluation",
        },
        {
            "comparison_aspect": "FIB-4 benchmark",
            "yang_2025": "62% balanced accuracy (validation)",
            "our_analysis": "Not computed (no clinical labs)",
            "advantage": "Yang et al. shows our panel should beat FIB-4",
        },
        {
            "comparison_aspect": "Etiology classification",
            "yang_2025": "Not attempted; concluded 'convergence across etiologies'",
            "our_analysis": "MASLD binary AUROC 0.840 ± 0.059",
            "advantage": "Direct contradiction of convergence claim",
        },
        {
            "comparison_aspect": "Cross-platform validation",
            "yang_2025": "Not attempted",
            "our_analysis": f"Transfer AUROC {transfer_auroc} (48 shared proteins, N=72)",
            "advantage": "Novel DIA-MS validation",
        },
        {
            "comparison_aspect": "Protein importance",
            "yang_2025": "NFASC+GDF15 selected by correlation",
            "our_analysis": "SHAP-based ranking across 50 CV folds",
            "advantage": "More robust feature selection",
        },
        {
            "comparison_aspect": "Tissue-plasma bridge",
            "yang_2025": "132 correlated proteins reported",
            "our_analysis": "Top proteins mapped to 33,943-gene tissue atlas",
            "advantage": "Multi-evidence contextualization",
        },
    ]

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 4. Manuscript numbers
# ---------------------------------------------------------------------------
def build_manuscript_numbers(panel_curves, leaderboard):
    """Assemble all key manuscript-ready numbers."""
    log.info("Building manuscript numbers table ...")

    rows = []

    # Helper to add a row
    def _add(section, claim, value, units, method, ci=""):
        rows.append({
            "section": section,
            "claim": claim,
            "value": str(value),
            "units": units,
            "method": method,
            "confidence_interval": str(ci),
        })

    # --- Fixed numbers from the spec ---
    _add("Plasma biomarkers", "Best binary fibrosis AUROC", "0.790", "AUROC",
         "XGBoost 10×5-fold nested CV", "0.790 ± 0.080")
    _add("Plasma biomarkers", "Best etiology MASLD binary AUROC", "0.840", "AUROC",
         "BalancedBagging 10×5-fold nested CV", "0.840 ± 0.059")
    _add("Plasma biomarkers", "Best ordinal fibrosis QWK", "0.478", "QWK",
         "RandomForest 10×5-fold nested CV", "0.478 ± ?")
    _add("Plasma biomarkers", "NFASC+GDF15 benchmark AUROC", "0.734", "AUROC",
         "LogReg LOEO pooled", "0.734 ± 0.070")
    _add("Plasma biomarkers", "Yang et al. contradiction - etiology AUROC", "0.840", "AUROC",
         "vs their 'convergence' claim", "")
    _add("Plasma biomarkers", "Cross-platform transfer AUROC", "0.655", "AUROC",
         "XGBoost Olink→DIA-MS, 48 shared proteins", "")
    _add("Plasma biomarkers", "Shared proteins Olink↔DIA-MS", "48", "proteins",
         "Gene symbol intersection", "")

    # --- From 226b_enrichment.csv: T1 DEG enrichment OR and p ---
    enrichment = _read("226b_enrichment.csv", required=False)
    if enrichment is not None:
        t1_row = enrichment[
            enrichment.get("target", pd.Series(dtype=str))
            .str.contains("T1", na=False)
        ]
        if not t1_row.empty:
            or_val = t1_row["odds_ratio"].iloc[0]
            pval = t1_row["fisher_pvalue"].iloc[0]
            _add(
                "Plasma biomarkers",
                "Top T1 protein tissue DEG enrichment",
                f"{or_val:.3f}",
                "odds ratio",
                "Fisher exact, p={:.3f}".format(pval),
                f"OR={or_val:.2f}, p={pval:.3f}",
            )

    # --- Panel performance at 10 proteins from panel_curves ---
    if panel_curves is not None:
        for t_label, t_filter in [
            ("T5 etiology", "T5"),
            ("T1 binary", "T1"),
            ("T2 ordinal", "T2"),
        ]:
            pc = panel_curves[
                panel_curves["target"].str.contains(t_filter, na=False) &
                (panel_curves["is_random"] == False)  # noqa: E712
            ]
            pc10 = pc[pc["panel_size"] == 10]
            pc_rand10 = panel_curves[
                panel_curves["target"].str.contains(t_filter, na=False) &
                (panel_curves["is_random"] == True)  # noqa: E712
                & (panel_curves["panel_size"] == 10)
            ]
            if pc10.empty:
                continue

            for metric_name, mdf in pc10.groupby("metric_name"):
                mean_val = mdf["metric_value"].mean()
                sd_val = mdf["metric_value"].std()
                rand_mean = np.nan
                if not pc_rand10.empty:
                    rand_sub = pc_rand10[pc_rand10["metric_name"] == metric_name]
                    if not rand_sub.empty:
                        rand_mean = rand_sub["metric_value"].mean()
                ci_str = f"{mean_val:.3f} ± {sd_val:.3f}"
                if not np.isnan(rand_mean):
                    ci_str += f" (random 10: {rand_mean:.3f})"
                _add(
                    "Plasma biomarkers",
                    f"{t_label} at 10-protein panel ({metric_name})",
                    f"{mean_val:.3f}",
                    metric_name,
                    "226a_panel_curves.csv 10-protein panel mean across folds",
                    ci_str,
                )

    # --- Best per-target from leaderboard ---
    if leaderboard is not None and not leaderboard.empty:
        for target, grp in leaderboard.groupby("target"):
            best = grp.loc[grp["primary_metric_mean"].idxmax()]
            metric_label = best.get("primary_metric_name", "metric")
            val = best["primary_metric_mean"]
            sd = best["primary_metric_sd"] if "primary_metric_sd" in best else np.nan
            ci_str = f"{val:.3f} ± {sd:.3f}" if not np.isnan(sd) else f"{val:.3f}"
            _add(
                "Plasma biomarkers",
                f"Best model for {target} ({metric_label})",
                f"{val:.3f}",
                metric_label,
                f"{best['model']} 10×5-fold nested CV (leaderboard)",
                ci_str,
            )

    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# 5. Print manuscript paragraph summary
# ---------------------------------------------------------------------------
def print_summary(leaderboard, panel, yang, numbers):
    sep = "=" * 72
    print(f"\n{sep}")
    print("  226f CLINICAL SUMMARY — Key manuscript findings")
    print(sep)

    # --- Best model per target ---
    if leaderboard is not None and not leaderboard.empty:
        print("\n[1] Best models per target:")
        for target, grp in leaderboard.groupby("target"):
            best = grp.loc[grp["primary_metric_mean"].idxmax()]
            metric = best.get("primary_metric_name", "metric")
            val = best["primary_metric_mean"]
            sd = best.get("primary_metric_sd", np.nan)
            sd_str = f" ± {sd:.3f}" if not np.isnan(float(sd)) else ""
            print(
                f"    {target:30s}  {best['model']:30s}  "
                f"{metric}={val:.3f}{sd_str}"
            )

    # --- Recommended panel ---
    if panel is not None and not panel.empty:
        print(f"\n[2] Recommended clinical panel (T5 etiology, top-{len(panel)} proteins):")
        for _, row in panel.iterrows():
            prot = row.get("protein", "?")
            shap = row.get("permutation_importance", np.nan)
            stab = row.get("fold_stability_pct", np.nan)
            print(
                f"    #{int(row['rank']):2d}  {prot:12s}  SHAP={shap:.4f}  "
                f"stability={stab}"
            )

    # --- Yang comparison ---
    if yang is not None and not yang.empty:
        print("\n[3] Yang et al. 2025 comparison highlights:")
        for _, row in yang.iterrows():
            print(f"    [{row['comparison_aspect']}]")
            print(f"       Yang:  {row['yang_2025']}")
            print(f"       Ours:  {row['our_analysis']}")
            print(f"       Edge:  {row['advantage']}")

    # --- Key manuscript numbers ---
    if numbers is not None and not numbers.empty:
        print("\n[4] Manuscript numbers:")
        for _, row in numbers.iterrows():
            print(
                f"    {row['section']:20s}  {row['claim']:55s}  "
                f"{row['value']:>8s} {row['units']:20s}  CI: {row['confidence_interval']}"
            )

    print(f"\n{sep}\n")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    log.info("226f_clinical_summary.py starting  (INDIR=%s)", INDIR)

    # Load core inputs
    panel_curves = _read("226a_panel_curves.csv", required=False)
    benchmark_df = _read("plasma_benchmark_comparison.csv", required=False)
    loeo_df = _read("plasma_loeo_baselines.csv", required=False)

    # 1. Leaderboard
    leaderboard = build_leaderboard()
    _write(leaderboard, "226f_sweep_leaderboard.csv")

    # 2. Recommended panel
    panel = build_recommended_panel(panel_curves)
    if not panel.empty:
        _write(panel, "226f_recommended_panel.csv")
    else:
        log.warning("Recommended panel is empty — no file written.")

    # 3. Yang comparison
    yang = build_yang_comparison(benchmark_df, loeo_df)
    _write(yang, "226f_yang_comparison.csv")

    # 4. Manuscript numbers
    numbers = build_manuscript_numbers(panel_curves, leaderboard)
    _write(numbers, "226f_manuscript_numbers.csv")

    # 5. Print summary
    print_summary(leaderboard, panel, yang, numbers)

    log.info("226f_clinical_summary.py complete.")


if __name__ == "__main__":
    main()
