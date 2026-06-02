#!/usr/bin/env python
"""
343d_fstage_external_validation.py

External validation of the scVI-inferred F-stage predictions (Script 343b) across
the MASLD scRNA donor roster (269 donors; 260 with predictions: 58 documented
training + 202 scvi_predicted).

The classifier was trained on a single cohort (Andrews/GSE202379, n=58).
This script does NOT re-train; it tests transferability against external
metadata (disease_stage_coarse / disease_stage_numeric) and biological priors
(pseudotime, progressor_frac, cNMF programs, cell-type fractions).

Tasks
-----
1. Cross-tab predicted F-stage vs disease_stage_coarse (all 260).
2. Per-dataset Spearman rho(F_stage_inferred, disease_stage_numeric).
3. Continuous biology axes: pseudotime, progressor_frac, cnmf_global_k16_P11.
4. Biological-prior tests with bootstrap CIs (frac_Hepatocytes, frac_Fibroblasts,
   frac_Macrophages).
5. Method-agreement (5 methods, 162 donors with all calls).
6. Confidence-stratified analysis (max posterior).
7. Compile single multi-panel PDF + summary text + per-donor tsv.

Constraints: PDF only (no PNG); MASLD palette only; do not re-train; do not
modify 343/343b/343c/344.
"""
from __future__ import annotations
import os
import sys
import json
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap

from scipy import stats

PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

STAGE_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
META_TSV  = STAGE_DIR / "donor_metadata_extended.tsv"
SCVI_TSV  = STAGE_DIR / "donor_fstage_scvi_predicted.tsv"
ALT_TSV   = STAGE_DIR / "donor_fstage_alt_predictions.tsv"

OUT_DONOR_TSV = STAGE_DIR / "fstage_scvi_validation.tsv"
OUT_SUMMARY   = STAGE_DIR / "fstage_scvi_validation_summary.txt"
OUT_PDF       = PROJECT_ROOT / "figures/supplementary/stage_ccc/figS_fstage_validation.pdf"

# MASLD palette (Magenta/Pink/Blue family). Stay aligned with publication_theme.R.
# Fibrosis ramp: F0 light blue -> F4 deep blue.
FSTAGE_COLORS = {
    0: "#E3F2FD", 1: "#90CAF9", 2: "#42A5F5", 3: "#1565C0", 4: "#0D47A1",
}
COARSE_COLORS = {
    "Healthy":          "#9E9E9E",   # neutral gray for control
    "Steatosis":        "#F57F17",   # amber
    "Steatohepatitis":  "#C2185B",   # magenta (disease)
    "Cirrhosis":        "#7B1FA2",   # deep purple (end-stage)
}
COARSE_ORDER = ["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"]
DATASET_COLOR = "#C2185B"            # primary MASLD magenta
DATASET_LO    = "#1565C0"            # blue for negative

# Custom 2-color heatmap: white -> deep magenta
HEAT_CMAP = LinearSegmentedColormap.from_list("masld_heat", ["#FFFFFF", "#C2185B"])


# --- helpers --------------------------------------------------------------

def spearman_safe(x, y):
    """Return (rho, p, n) or (nan, nan, n) if insufficient variance."""
    x = pd.Series(x); y = pd.Series(y)
    mask = x.notna() & y.notna()
    if mask.sum() < 3:
        return np.nan, np.nan, int(mask.sum())
    xv = x[mask].values; yv = y[mask].values
    if np.unique(xv).size < 2 or np.unique(yv).size < 2:
        return np.nan, np.nan, int(mask.sum())
    rho, p = stats.spearmanr(xv, yv)
    return float(rho), float(p), int(mask.sum())


def bootstrap_spearman_ci(x, y, n_boot=2000, seed=42, alpha=0.05):
    rng = np.random.default_rng(seed)
    x = np.asarray(x); y = np.asarray(y)
    mask = ~(np.isnan(x) | np.isnan(y))
    x = x[mask]; y = y[mask]
    n = x.size
    if n < 3:
        return np.nan, np.nan, np.nan, n
    point, _ = stats.spearmanr(x, y)
    boots = np.empty(n_boot)
    for i in range(n_boot):
        idx = rng.integers(0, n, n)
        xb = x[idx]; yb = y[idx]
        if np.unique(xb).size < 2 or np.unique(yb).size < 2:
            boots[i] = np.nan
        else:
            boots[i] = stats.spearmanr(xb, yb)[0]
    lo = np.nanpercentile(boots, 100 * alpha / 2)
    hi = np.nanpercentile(boots, 100 * (1 - alpha / 2))
    return float(point), float(lo), float(hi), n


def format_contingency(idx_rows, col_order, df, row_col, col_col):
    tab = pd.crosstab(df[row_col], df[col_col]).reindex(index=idx_rows, columns=col_order, fill_value=0)
    return tab


# --- load -----------------------------------------------------------------

def main():
    print(f"[load] meta:  {META_TSV}")
    meta = pd.read_csv(META_TSV, sep="\t")
    print(f"[load] scvi:  {SCVI_TSV}")
    scvi = pd.read_csv(SCVI_TSV, sep="\t")
    print(f"[load] alt:   {ALT_TSV}")
    alt = pd.read_csv(ALT_TSV, sep="\t")

    df = meta.copy()
    # Bring posteriors in if missing
    if "P_F0" not in df.columns:
        df = df.merge(scvi[["sample", "P_F0", "P_F1", "P_F2", "P_F3", "P_F4"]], on="sample", how="left")
    # max posterior confidence
    p_cols = ["P_F0", "P_F1", "P_F2", "P_F3", "P_F4"]
    df["posterior_max"] = df[p_cols].max(axis=1)

    # Subset donors with a non-null inferred F-stage
    pred_mask = df["F_stage_inferred"].notna()
    dpred = df.loc[pred_mask].copy()
    dpred["F_stage_inferred"] = dpred["F_stage_inferred"].astype(int)
    print(f"[load] n donors total = {len(df)}; n with F_stage_inferred = {len(dpred)}")

    # ---------- Task 1: cross-tab + overall correlation ------------------
    df1 = dpred.dropna(subset=["disease_stage_coarse", "disease_stage_numeric"])
    tab_overall = format_contingency([0, 1, 2, 3, 4], COARSE_ORDER, df1, "F_stage_inferred", "disease_stage_coarse")
    rho_overall, p_overall, n_overall = spearman_safe(df1["F_stage_inferred"], df1["disease_stage_numeric"])

    # Biologically implausible cells: F0 -> Cirrhosis; F4 -> Healthy
    implaus_cells = {
        "F0_Cirrhosis": int(tab_overall.loc[0, "Cirrhosis"]),
        "F4_Healthy":   int(tab_overall.loc[4, "Healthy"]),
        "F0_Steatohepatitis": int(tab_overall.loc[0, "Steatohepatitis"]),
        "F4_Steatosis": int(tab_overall.loc[4, "Steatosis"]),
    }

    # ---------- Task 2: per-dataset ----------
    perds_rows = []
    perds_confusion = {}
    for ds, g in df1.groupby("dataset"):
        rho, p, n = spearman_safe(g["F_stage_inferred"], g["disease_stage_numeric"])
        n_stages = g["disease_stage_coarse"].nunique()
        perds_rows.append(dict(dataset=ds, n_donors=int(n), rho=rho, p_value=p, n_stages=int(n_stages)))
        perds_confusion[ds] = format_contingency([0, 1, 2, 3, 4], COARSE_ORDER, g, "F_stage_inferred", "disease_stage_coarse")
    perds = pd.DataFrame(perds_rows).sort_values("n_donors", ascending=False).reset_index(drop=True)

    # ---------- Task 3: continuous biology ----------
    cont_axes = [
        ("macrophage_pseudotime_mean", "Macrophage pseudotime"),
        ("hepatocyte_pseudotime_mean", "Hepatocyte pseudotime"),
        ("progressor_frac",            "Hep Progressor fraction"),
        ("cnmf_global_k16_P11",        "cNMF program P11 (inflammatory)"),
    ]
    cont_rows = []
    for col, lbl in cont_axes:
        if col in dpred.columns:
            rho, p, n = spearman_safe(dpred["F_stage_inferred"], pd.to_numeric(dpred[col], errors="coerce"))
        else:
            rho, p, n = np.nan, np.nan, 0
        cont_rows.append(dict(axis=col, label=lbl, rho=rho, p_value=p, n=n))
    cont_df = pd.DataFrame(cont_rows)

    # ---------- Task 4: biological priors w/ bootstrap CI ----------
    priors = [
        ("frac_Hepatocytes", "frac_Hepatocytes (expect NEG)", -1),
        ("frac_Fibroblasts", "frac_Fibroblasts (expect POS)", +1),
        ("frac_Macrophages", "frac_Macrophages (expect POS)", +1),
    ]
    prior_rows = []
    for col, lbl, expected_sign in priors:
        if col not in dpred.columns:
            prior_rows.append(dict(axis=col, label=lbl, expected_sign=expected_sign,
                                    rho=np.nan, ci_lo=np.nan, ci_hi=np.nan, n=0,
                                    direction_matches=False))
            continue
        x = pd.to_numeric(dpred[col], errors="coerce").values
        y = dpred["F_stage_inferred"].astype(float).values
        rho, lo, hi, n = bootstrap_spearman_ci(x, y, n_boot=2000, seed=42)
        match = (rho > 0 and expected_sign > 0) or (rho < 0 and expected_sign < 0)
        prior_rows.append(dict(axis=col, label=lbl, expected_sign=expected_sign,
                                rho=rho, ci_lo=lo, ci_hi=hi, n=n,
                                direction_matches=bool(match)))
    prior_df = pd.DataFrame(prior_rows)

    # ---------- Task 5: method agreement ----------
    methods = ["F_stage_scvi", "F_stage_knn_scvi", "F_stage_cnmf", "F_stage_pseudotime", "F_stage_celltype"]
    alt_sub = alt[["sample"] + methods].copy()
    # 5-method agreement requires all non-NaN
    all5 = alt_sub.dropna(subset=methods).copy()
    n_all5 = len(all5)
    pair_mat = pd.DataFrame(np.nan, index=methods, columns=methods)
    for i, a in enumerate(methods):
        for j, b in enumerate(methods):
            if j < i:
                continue
            if a == b:
                pair_mat.loc[a, b] = 1.0
                continue
            rho, _, _ = spearman_safe(all5[a], all5[b])
            pair_mat.loc[a, b] = rho
            pair_mat.loc[b, a] = rho

    # Donors where scVI logistic vs scVI kNN diverge by >=2 F-stages
    div = alt_sub.dropna(subset=["F_stage_scvi", "F_stage_knn_scvi"]).copy()
    div["scvi_vs_knn_abs_diff"] = (div["F_stage_scvi"] - div["F_stage_knn_scvi"]).abs()
    flagged = div[div["scvi_vs_knn_abs_diff"] >= 2][["sample", "F_stage_scvi", "F_stage_knn_scvi", "scvi_vs_knn_abs_diff"]]
    flagged = flagged.merge(meta[["sample", "dataset", "disease_stage_coarse"]], on="sample", how="left")
    flagged_sorted = flagged.sort_values("scvi_vs_knn_abs_diff", ascending=False)

    # ---------- Task 6: confidence-stratified ----------
    df_conf = df1.copy()
    df_conf["confidence_bin"] = pd.cut(df_conf["posterior_max"],
                                       bins=[-0.01, 0.4, 0.7, 1.01],
                                       labels=["low", "medium", "high"])
    conf_rows = []
    conf_tabs = {}
    for b in ["low", "medium", "high"]:
        g = df_conf[df_conf["confidence_bin"] == b]
        rho, p, n = spearman_safe(g["F_stage_inferred"], g["disease_stage_numeric"])
        # match: predicted F maps to coarse stage in expected mapping
        # F0 -> Healthy/Steatosis ; F4 -> Cirrhosis ; F1-F3 -> SH
        def expected_match(row):
            f = int(row["F_stage_inferred"]); c = row["disease_stage_coarse"]
            if f == 0 and c in ("Healthy", "Steatosis"): return 1
            if f in (1, 2, 3) and c == "Steatohepatitis": return 1
            if f == 4 and c == "Cirrhosis": return 1
            return 0
        if len(g) > 0:
            match_rate = g.apply(expected_match, axis=1).mean()
        else:
            match_rate = np.nan
        conf_rows.append(dict(bin=b, n=int(len(g)), rho=rho, p_value=p, expected_match_rate=match_rate))
        conf_tabs[b] = format_contingency([0, 1, 2, 3, 4], COARSE_ORDER, g, "F_stage_inferred", "disease_stage_coarse")
    conf_df = pd.DataFrame(conf_rows)

    # ---------- Build per-donor validation tsv ----------
    out_cols = ["sample", "dataset", "disease_stage_coarse", "disease_stage_numeric",
                "F_stage_documented", "F_stage_inferred", "F_stage_source",
                "P_F0", "P_F1", "P_F2", "P_F3", "P_F4", "posterior_max"]
    out_cols = [c for c in out_cols if c in df.columns]
    donor_out = df[out_cols].copy()
    donor_out["confidence_bin"] = pd.cut(donor_out.get("posterior_max", pd.Series(np.nan, index=donor_out.index)),
                                          bins=[-0.01, 0.4, 0.7, 1.01],
                                          labels=["low", "medium", "high"])
    # add expected_match flag for donors with both inferred + coarse
    def _em(row):
        if pd.isna(row["F_stage_inferred"]) or pd.isna(row["disease_stage_coarse"]): return np.nan
        f = int(row["F_stage_inferred"]); c = row["disease_stage_coarse"]
        if f == 0 and c in ("Healthy", "Steatosis"): return 1
        if f in (1, 2, 3) and c == "Steatohepatitis": return 1
        if f == 4 and c == "Cirrhosis": return 1
        return 0
    donor_out["expected_coarse_match"] = donor_out.apply(_em, axis=1)
    # scvi vs knn flag
    knn_map = alt.set_index("sample")["F_stage_knn_scvi"].to_dict()
    scv_map = alt.set_index("sample")["F_stage_scvi"].to_dict()
    donor_out["F_stage_scvi_logistic"] = donor_out["sample"].map(scv_map)
    donor_out["F_stage_scvi_knn"]      = donor_out["sample"].map(knn_map)
    donor_out["scvi_method_abs_diff"] = (donor_out["F_stage_scvi_logistic"] - donor_out["F_stage_scvi_knn"]).abs()
    donor_out["flagged_method_disagreement"] = (donor_out["scvi_method_abs_diff"] >= 2).astype("Int64")

    donor_out.to_csv(OUT_DONOR_TSV, sep="\t", index=False)
    print(f"[write] {OUT_DONOR_TSV} ({len(donor_out)} rows)")

    # ---------- Write summary text ----------
    lines = []
    L = lines.append
    L("# scVI F-stage external validation summary")
    L(f"# inputs:")
    L(f"#   meta = {META_TSV.relative_to(PROJECT_ROOT)}")
    L(f"#   scvi = {SCVI_TSV.relative_to(PROJECT_ROOT)}")
    L(f"#   alt  = {ALT_TSV.relative_to(PROJECT_ROOT)}")
    L(f"# donors loaded: {len(df)}  |  predicted: {len(dpred)}  |  documented: {int((df['F_stage_source']=='documented').sum())}")
    L("")

    L("## 1. Overall cross-tab F_stage_inferred x disease_stage_coarse")
    L(tab_overall.to_string())
    L(f"\nOverall Spearman rho(F_stage_inferred, disease_stage_numeric) = {rho_overall:.3f}  (p={p_overall:.2e}, n={n_overall})")
    L(f"Biologically implausible cells (counts):")
    for k, v in implaus_cells.items():
        L(f"   {k}: {v}")
    L("")

    L("## 2. Per-dataset Spearman rho(F_stage_inferred, disease_stage_numeric)")
    L(perds.to_string(index=False))
    L("")
    L("Per-dataset confusion (F_stage_inferred rows x disease_stage_coarse cols):")
    for ds, t in perds_confusion.items():
        L(f"\n### {ds}  (n={int(t.values.sum())})")
        L(t.to_string())
    L("")

    L("## 3. Continuous biology axes")
    L(cont_df.to_string(index=False))
    L("")

    L("## 4. Biological-prior tests (Spearman rho with 95% bootstrap CI; 2000 boot)")
    L(prior_df.to_string(index=False))
    L("")

    L(f"## 5. Method agreement (Spearman rho, n_all5={n_all5})")
    L(pair_mat.round(3).to_string())
    L(f"\nDonors with scVI-logistic vs scVI-kNN diff >= 2 F-stages: n={len(flagged_sorted)}")
    if len(flagged_sorted) > 0:
        L(flagged_sorted.to_string(index=False))
    L("")

    L("## 6. Confidence-stratified analysis (max posterior bins)")
    L(conf_df.to_string(index=False))
    for b in ["low", "medium", "high"]:
        L(f"\n### bin = {b}  (n={int(conf_tabs[b].values.sum())})")
        L(conf_tabs[b].to_string())
    L("")

    L("## Recommendation drafting context")
    L(f"  - overall rho = {rho_overall:.3f}")
    L(f"  - n datasets with rho > 0.5 (transferable): "
       f"{int((perds['rho']>0.5).sum())} / {len(perds)}")
    L(f"  - n datasets with rho <= 0.0 (anti-correlated): "
       f"{int((perds['rho']<=0).sum())} / {len(perds)}")
    matching = prior_df["direction_matches"].sum()
    L(f"  - biological priors matching expected direction: {matching} / 3")

    OUT_SUMMARY.write_text("\n".join(lines) + "\n")
    print(f"[write] {OUT_SUMMARY}")

    # ---------- PDF (multi-panel) ----------
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(OUT_PDF) as pdf:

        # ---- Page 1: overall cross-tab heatmap + per-dataset rho bar ----
        fig = plt.figure(figsize=(11, 5.5))
        gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.2], wspace=0.4)

        ax1 = fig.add_subplot(gs[0, 0])
        mat = tab_overall.values.astype(float)
        # row-normalize for color
        row_sums = mat.sum(axis=1, keepdims=True); row_sums[row_sums==0] = 1
        mat_norm = mat / row_sums
        im = ax1.imshow(mat_norm, cmap=HEAT_CMAP, vmin=0, vmax=1, aspect="auto")
        ax1.set_xticks(range(len(COARSE_ORDER))); ax1.set_xticklabels(COARSE_ORDER, rotation=30, ha="right")
        ax1.set_yticks(range(5)); ax1.set_yticklabels([f"F{i}" for i in range(5)])
        for i in range(5):
            for j in range(len(COARSE_ORDER)):
                ax1.text(j, i, int(mat[i, j]), ha="center", va="center",
                         color="white" if mat_norm[i, j] > 0.5 else "black",
                         fontsize=8)
        ax1.set_title(f"All 260 donors: F_stage_inferred vs disease_stage_coarse\nrow-normalised; rho_overall = {rho_overall:.2f}",
                      fontsize=9)
        cbar = fig.colorbar(im, ax=ax1, fraction=0.045, pad=0.04); cbar.ax.tick_params(labelsize=7)
        cbar.set_label("row fraction", fontsize=8)

        ax2 = fig.add_subplot(gs[0, 1])
        bars = perds.copy()
        colors = [DATASET_COLOR if (not pd.isna(v) and v >= 0) else DATASET_LO for v in bars["rho"]]
        ax2.barh(bars["dataset"], bars["rho"].fillna(0), color=colors, edgecolor="black", linewidth=0.4)
        for i, (rho_v, n_v) in enumerate(zip(bars["rho"], bars["n_donors"])):
            txt = f"n={n_v}, rho={rho_v:.2f}" if not pd.isna(rho_v) else f"n={n_v}, NA"
            ax2.text(0.02 if pd.isna(rho_v) or rho_v < 0 else rho_v + 0.02, i, txt,
                     va="center", fontsize=7)
        ax2.axvline(0, color="black", linewidth=0.4)
        ax2.axvline(0.5, color="#9E9E9E", linewidth=0.4, linestyle="--")
        ax2.set_xlabel("Spearman rho(F_stage_inferred, disease_stage_numeric)", fontsize=8)
        ax2.set_xlim(-1, 1)
        ax2.set_title("Per-dataset transferability (rho > 0.5 = solid)", fontsize=9)
        ax2.tick_params(axis="both", labelsize=8)

        fig.suptitle("scVI F-stage external validation", fontsize=11, y=1.02)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # ---- Page 2: scatter of F-stage vs 4 continuous axes ----
        fig, axes = plt.subplots(2, 2, figsize=(9, 8))
        for (col, lbl), ax in zip(cont_axes, axes.flat):
            if col not in dpred.columns:
                ax.set_visible(False); continue
            v = pd.to_numeric(dpred[col], errors="coerce")
            mask = v.notna() & dpred["F_stage_inferred"].notna()
            xs = dpred.loc[mask, "F_stage_inferred"].astype(int).values
            ys = v[mask].values
            # jitter F-stage horizontally
            rng = np.random.default_rng(7)
            xj = xs + rng.uniform(-0.18, 0.18, size=xs.size)
            cols = [FSTAGE_COLORS[int(f)] for f in xs]
            ax.scatter(xj, ys, c=cols, edgecolor="black", linewidth=0.3, s=22, alpha=0.85)
            # boxplot underneath
            for f in range(5):
                yy = ys[xs == f]
                if yy.size > 1:
                    ax.boxplot([yy], positions=[f], widths=0.45, showfliers=False,
                               boxprops=dict(color="black", linewidth=0.5),
                               medianprops=dict(color="black", linewidth=0.8),
                               whiskerprops=dict(color="black", linewidth=0.5),
                               capprops=dict(color="black", linewidth=0.5))
            rho, p, n = spearman_safe(xs, ys)
            ax.set_title(f"{lbl}\nrho = {rho:.2f}  (n={n}, p={p:.1e})", fontsize=8)
            ax.set_xlabel("F_stage_inferred", fontsize=8)
            ax.set_ylabel(col, fontsize=8)
            ax.set_xticks([0,1,2,3,4])
            ax.set_xticklabels([f"F{i}" for i in range(5)])
            ax.tick_params(axis="both", labelsize=7)
            for s in ("top","right"): ax.spines[s].set_visible(False)
        fig.suptitle("F-stage vs continuous biology axes", fontsize=11)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # ---- Page 3: biological priors (forest) + cell-type fractions scatter ----
        fig = plt.figure(figsize=(11, 5.5))
        gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.4], wspace=0.4)

        ax1 = fig.add_subplot(gs[0, 0])
        for i, row in prior_df.iterrows():
            ec = "black"
            fc = DATASET_COLOR if row["direction_matches"] else DATASET_LO
            ax1.errorbar(row["rho"], i, xerr=[[row["rho"]-row["ci_lo"]], [row["ci_hi"]-row["rho"]]],
                         fmt="o", ms=8, mfc=fc, mec=ec, ecolor="black", elinewidth=1.0, capsize=3)
        ax1.axvline(0, color="black", linewidth=0.5)
        ax1.set_yticks(range(len(prior_df)))
        ax1.set_yticklabels(prior_df["label"], fontsize=8)
        ax1.set_xlabel("Spearman rho (95% bootstrap CI)", fontsize=8)
        ax1.set_xlim(-1, 1)
        ax1.set_title("Biological-prior tests", fontsize=9)
        for s in ("top","right"): ax1.spines[s].set_visible(False)
        ax1.tick_params(axis="x", labelsize=7)

        ax2 = fig.add_subplot(gs[0, 1])
        # confidence-stratified accuracy
        bin_names = ["low", "medium", "high"]
        match_rates = []
        bin_ns = []
        for b in bin_names:
            sub = conf_df[conf_df["bin"] == b]
            mr = float(sub["expected_match_rate"].iloc[0]) if len(sub) and not pd.isna(sub["expected_match_rate"].iloc[0]) else 0.0
            n  = int(sub["n"].iloc[0]) if len(sub) else 0
            match_rates.append(mr); bin_ns.append(n)
        bar_colors = ["#90CAF9", "#42A5F5", "#1565C0"]
        bars = ax2.bar(bin_names, match_rates, color=bar_colors, edgecolor="black", linewidth=0.4)
        for b, mr, n in zip(bars, match_rates, bin_ns):
            ax2.text(b.get_x()+b.get_width()/2, mr+0.02, f"{mr:.2f}\nn={n}",
                     ha="center", va="bottom", fontsize=7)
        ax2.set_ylim(0, 1.1)
        ax2.set_ylabel("expected_coarse_match rate", fontsize=8)
        ax2.set_xlabel("max posterior bin", fontsize=8)
        ax2.set_title("Confidence-stratified concordance with disease_stage_coarse", fontsize=9)
        for s in ("top","right"): ax2.spines[s].set_visible(False)
        ax2.tick_params(axis="both", labelsize=7)

        fig.suptitle("Biological-prior tests + confidence stratification", fontsize=11, y=1.02)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

        # ---- Page 4: pairwise method-agreement heatmap + flagged-donor table ----
        fig = plt.figure(figsize=(11, 5.5))
        gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.2], wspace=0.45)

        ax1 = fig.add_subplot(gs[0, 0])
        mat = pair_mat.values.astype(float)
        im = ax1.imshow(mat, cmap=HEAT_CMAP, vmin=0, vmax=1, aspect="auto")
        labels_short = ["scVI-logit", "scVI-kNN", "cNMF", "pseudotime", "celltype"]
        ax1.set_xticks(range(5)); ax1.set_xticklabels(labels_short, rotation=30, ha="right", fontsize=8)
        ax1.set_yticks(range(5)); ax1.set_yticklabels(labels_short, fontsize=8)
        for i in range(5):
            for j in range(5):
                ax1.text(j, i, f"{mat[i,j]:.2f}", ha="center", va="center",
                         color="white" if mat[i,j] > 0.5 else "black", fontsize=7)
        ax1.set_title(f"Method-agreement (Spearman rho, n_all5={n_all5})", fontsize=9)
        cbar = fig.colorbar(im, ax=ax1, fraction=0.045, pad=0.04); cbar.ax.tick_params(labelsize=7)

        ax2 = fig.add_subplot(gs[0, 1])
        ax2.axis("off")
        ax2.set_title(f"scVI-logistic vs scVI-kNN disagreement >=2 (n={len(flagged_sorted)})", fontsize=9, loc="left")
        if len(flagged_sorted) > 0:
            tbl = flagged_sorted.head(18).copy()
            cell_text = tbl.astype(str).values.tolist()
            col_labels = list(tbl.columns)
            t = ax2.table(cellText=cell_text, colLabels=col_labels, loc="upper left", cellLoc="center")
            t.auto_set_font_size(False); t.set_fontsize(7); t.scale(1.0, 1.1)
        else:
            ax2.text(0.05, 0.5, "no donors flagged", fontsize=10, transform=ax2.transAxes)

        fig.suptitle("Cross-method agreement and flagged donors", fontsize=11, y=1.02)
        fig.tight_layout()
        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    print(f"[write] {OUT_PDF}")
    print("[done] 343d external validation")


if __name__ == "__main__":
    main()
