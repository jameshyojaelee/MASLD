#!/usr/bin/env python
"""
343s_dubious_threshold_sensitivity.py

Question
--------
How sensitive is the augmented scVI ordinal F-stage classifier's external
Spearman rho to the dubious-healthy filter threshold used in Script 343j?

Sweep grid
----------
  F4_frac cutoff:        0.20, 0.30, 0.40 (current), 0.50, 0.60
  dominant_fstage cutoff: >=2, >=3 (current), >=4

For each (F4_cut, dom_cut) cell:
  - dubious-healthy := disease_stage_coarse=="Healthy" AND
                       (dominant_fstage >= dom_cut OR F4_frac > F4_cut)
  - clean-healthy   := Healthy AND not dubious  (anchor -> F0)
  - cirrhosis       := disease_stage_coarse=="Cirrhosis"  (anchor -> F4)
  - documented      := F_stage_documented in {0..4}       (kept as-is)

Train the same augmented scVI ordinal logistic from 343m on
(documented U clean-healthy@F0 U cirrhosis@F4) using donor_hep_scvi_mean.tsv
latents and report:
  - n_train, n_dubious, n_clean
  - external Spearman rho vs disease_stage_numeric
  - clean-healthy Healthy->F4 LOOCV count
  - SH F1-F4 predicted fractions

Outputs
-------
  - results/.../dubious_threshold_sensitivity.tsv  (15 rows)
  - figures/supplementary/stage_ccc/figS_dubious_threshold_sensitivity.pdf
"""

from __future__ import annotations
import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from sklearn.metrics import cohen_kappa_score
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    import mord
    HAVE_MORD = True
except ImportError:
    HAVE_MORD = False

# ------------------------------------------------------------------ paths
PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
ST_DIR    = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"

DONOR_FRAC = ST_DIR / "donor_fstage_cellfrac.tsv"
DONOR_META = ST_DIR / "donor_metadata_extended.tsv"
HEP_PB     = ST_DIR / "donor_hep_scvi_mean.tsv"
DOC_F      = ST_DIR / "donor_fstage_documented.tsv"

OUT_TSV    = ST_DIR / "dubious_threshold_sensitivity.tsv"
FIG_DIR    = PROJECT_ROOT / "figures/supplementary/stage_ccc"
OUT_PDF    = FIG_DIR / "figS_dubious_threshold_sensitivity.pdf"
FIG_DIR.mkdir(parents=True, exist_ok=True)

LATENT_COLS = [f"z{i}" for i in range(20)]

# Sweep grid (rows -> F4_frac cut, cols -> dominant_fstage cut)
F4_CUTS  = [0.20, 0.30, 0.40, 0.50, 0.60]
DOM_CUTS = [2, 3, 4]


# ------------------------------------------------------------------ I/O
def load_inputs() -> pd.DataFrame:
    """Build the canonical donor table with latents + documented + cellfrac."""
    frac = pd.read_csv(DONOR_FRAC, sep="\t")
    meta = pd.read_csv(DONOR_META, sep="\t", low_memory=False)
    doc  = pd.read_csv(DOC_F, sep="\t")
    pb   = pd.read_csv(HEP_PB, sep="\t")

    doc["F_stage_documented"] = pd.to_numeric(
        doc["F_stage_documented"], errors="coerce")

    # Canonical donor backbone: meta order, 269 rows
    df = meta[["sample", "dataset", "disease_stage_coarse",
               "disease_stage_numeric"]].copy()

    keep_cols = ["sample", "dataset"] + LATENT_COLS
    if "n_hepatocytes_for_latent" in pb.columns:
        keep_cols.append("n_hepatocytes_for_latent")
    df = df.merge(pb[keep_cols], on=["sample", "dataset"], how="left")

    df = df.merge(doc[["sample", "dataset", "F_stage_documented"]],
                  on=["sample", "dataset"], how="left")

    # Bring in dominant_fstage and F4_frac fresh from cellfrac (definitive)
    f = frac[["sample", "dataset", "dominant_fstage", "F4_frac"]].copy()
    f["dominant_fstage"] = pd.to_numeric(f["dominant_fstage"], errors="coerce")
    f["F4_frac"]         = pd.to_numeric(f["F4_frac"], errors="coerce")
    df = df.merge(f, on=["sample", "dataset"], how="left")

    return df


# ----------------------------------------------------------- training set
def assemble_training(df: pd.DataFrame,
                      f4_cut: float,
                      dom_cut: int) -> tuple[pd.DataFrame, dict]:
    """
    For the (f4_cut, dom_cut) cell:
      dubious-healthy: Healthy AND (dominant_fstage >= dom_cut OR F4_frac > f4_cut)
      clean-healthy:   Healthy AND not dubious             -> F0 anchor
      cirrhosis:                                            -> F4 anchor
      documented:      kept as-is from donor_fstage_documented.tsv

    Returns (train_set, info_dict).
    """
    t = df.copy()
    t["y_train"] = np.nan
    t["origin"]  = "none"

    is_healthy   = t["disease_stage_coarse"].astype(str) == "Healthy"
    is_cirrhosis = t["disease_stage_coarse"].astype(str) == "Cirrhosis"

    dom = pd.to_numeric(t["dominant_fstage"], errors="coerce")
    f4  = pd.to_numeric(t["F4_frac"], errors="coerce")
    dubious = is_healthy & ((dom >= dom_cut) | (f4 > f4_cut))
    clean_healthy = is_healthy & (~dubious)

    # 1. Documented Andrews (takes precedence)
    mask_doc = t["F_stage_documented"].notna()
    t.loc[mask_doc, "y_train"] = t.loc[mask_doc, "F_stage_documented"]
    t.loc[mask_doc, "origin"]  = "documented"

    # 2. Clean-healthy => F0
    mask_chealthy = clean_healthy & (~mask_doc)
    t.loc[mask_chealthy, "y_train"] = 0
    t.loc[mask_chealthy, "origin"]  = "clean_healthy_anchor"

    # 3. Cirrhosis => F4
    mask_cirr = is_cirrhosis & (~mask_doc)
    t.loc[mask_cirr, "y_train"] = 4
    t.loc[mask_cirr, "origin"]  = "cirrhosis_anchor"

    has_latent = t[LATENT_COLS].notna().all(axis=1)
    train = t[t["y_train"].notna() & has_latent].copy()
    train["y_train"] = train["y_train"].astype(int)

    info = {
        "n_dubious":       int(dubious.sum()),
        "n_clean_healthy": int(clean_healthy.sum()),
        "n_clean_in_train":int(((train["origin"] == "clean_healthy_anchor")).sum()),
        "n_cirr_in_train": int((train["origin"] == "cirrhosis_anchor").sum()),
        "n_doc_in_train":  int((train["origin"] == "documented").sum()),
        "n_train":         int(len(train)),
    }
    return train, info


# ----------------------------------------------------------- classifier
def make_ordinal():
    """Ordinal logistic. mord.LogisticIT when available; LR fallback."""
    if HAVE_MORD:
        return Pipeline([("scale", StandardScaler()),
                         ("clf",   mord.LogisticIT(alpha=1.0))])
    return Pipeline([("scale", StandardScaler()),
                     ("clf",   LogisticRegression(solver="lbfgs", C=1.0,
                                                   max_iter=2000,
                                                   class_weight="balanced"))])


def predict_argmax(clf, X):
    if HAVE_MORD and isinstance(clf.named_steps["clf"], mord.LogisticIT):
        return clf.predict(X).astype(int)
    proba = clf.predict_proba(X)
    return clf.named_steps["clf"].classes_[proba.argmax(axis=1)].astype(int)


# ----------------------------------------------------------- single cell
def evaluate_cell(df: pd.DataFrame, f4_cut: float, dom_cut: int) -> dict:
    train, info = assemble_training(df, f4_cut, dom_cut)

    X = train[LATENT_COLS].to_numpy()
    y = train["y_train"].to_numpy().astype(int)
    origin = train["origin"].to_numpy()

    # Full LOOCV for QWK + clean-healthy Healthy->F4 count
    n = len(y)
    ypred = np.full(n, -1, dtype=int)
    for i in range(n):
        mask = np.ones(n, dtype=bool); mask[i] = False
        clf = make_ordinal()
        clf.fit(X[mask], y[mask])
        ypred[i] = predict_argmax(clf, X[i:i+1])[0]

    qwk_full = float(cohen_kappa_score(
        y, ypred, weights="quadratic", labels=list(range(5))))

    andrews_idx = np.where(origin == "documented")[0]
    qwk_andrews = (float(cohen_kappa_score(
        y[andrews_idx], ypred[andrews_idx],
        weights="quadratic", labels=list(range(5))))
        if len(andrews_idx) > 0 else np.nan)

    chealthy_idx = np.where(origin == "clean_healthy_anchor")[0]
    n_hF4 = int((ypred[chealthy_idx] == 4).sum()) if len(chealthy_idx) else 0
    n_clean_loo = int(len(chealthy_idx))

    # Full-train, predict everyone -> external rho + SH breakdown
    clf = make_ordinal()
    clf.fit(X, y)
    has_latent = df[LATENT_COLS].notna().all(axis=1)
    X_all = df.loc[has_latent, LATENT_COLS].to_numpy()
    pred_all = predict_argmax(clf, X_all)
    df_pred = df.loc[has_latent].copy()
    df_pred["F_pred"] = pred_all

    valid = df_pred.dropna(subset=["disease_stage_numeric"])
    if len(valid) >= 3:
        rho, p_rho = stats.spearmanr(
            valid["F_pred"], valid["disease_stage_numeric"])
    else:
        rho, p_rho = np.nan, np.nan

    sh = df_pred[df_pred["disease_stage_coarse"] == "Steatohepatitis"]
    sh_pred = sh["F_pred"].astype(int)
    sh_counts = sh_pred.value_counts().reindex(range(5), fill_value=0)
    sh_frac = (sh_counts / max(len(sh_pred), 1)).round(4)

    return {
        "f4_cut":   f4_cut,
        "dom_cut":  dom_cut,
        "n_dubious":       info["n_dubious"],
        "n_clean_healthy": info["n_clean_healthy"],
        "n_doc_train":     info["n_doc_in_train"],
        "n_clean_train":   info["n_clean_in_train"],
        "n_cirr_train":    info["n_cirr_in_train"],
        "n_train":         info["n_train"],
        "qwk_full_loocv":  round(qwk_full, 4),
        "qwk_andrews_loocv": round(qwk_andrews, 4) if not np.isnan(qwk_andrews) else np.nan,
        "external_rho":    round(float(rho), 4) if not np.isnan(rho) else np.nan,
        "external_rho_p":  "{:.3e}".format(p_rho) if not np.isnan(p_rho) else "",
        "n_external":      int(len(valid)),
        "n_healthy_to_F4": n_hF4,
        "n_chealthy_loo":  n_clean_loo,
        "sh_F0_frac":      float(sh_frac.iloc[0]),
        "sh_F1_frac":      float(sh_frac.iloc[1]),
        "sh_F2_frac":      float(sh_frac.iloc[2]),
        "sh_F3_frac":      float(sh_frac.iloc[3]),
        "sh_F4_frac":      float(sh_frac.iloc[4]),
        "n_sh":            int(len(sh_pred)),
    }


# ----------------------------------------------------------- figure
def make_figure(grid: pd.DataFrame):
    """Heatmap of external rho + heatmap of Healthy->F4 + line plot."""
    with PdfPages(OUT_PDF) as pdf:
        # 2 x 2 layout
        fig, axes = plt.subplots(2, 2, figsize=(13, 11))

        # ---- (a) heatmap external rho ------------------------------------
        ax = axes[0, 0]
        rho_mat = grid.pivot(index="f4_cut", columns="dom_cut",
                             values="external_rho").reindex(
            index=F4_CUTS, columns=DOM_CUTS)
        im = ax.imshow(rho_mat.values, cmap="viridis", aspect="auto",
                       vmin=np.nanmin(rho_mat.values),
                       vmax=np.nanmax(rho_mat.values))
        ax.set_xticks(range(len(DOM_CUTS)))
        ax.set_xticklabels([f">={c}" for c in DOM_CUTS])
        ax.set_yticks(range(len(F4_CUTS)))
        ax.set_yticklabels([f"{c:.2f}" for c in F4_CUTS])
        ax.set_xlabel("dominant_fstage cutoff")
        ax.set_ylabel("F4_frac cutoff")
        ax.set_title(r"(a) External Spearman $\rho$ vs disease_stage_numeric")
        for i, f4 in enumerate(F4_CUTS):
            for j, dc in enumerate(DOM_CUTS):
                v = rho_mat.iloc[i, j]
                txt = f"{v:.3f}" if not np.isnan(v) else "NA"
                mark = " *" if (f4 == 0.40 and dc == 3) else ""
                ax.text(j, i, txt + mark, ha="center", va="center",
                        color="white" if v < np.nanmean(rho_mat.values) else "black",
                        fontsize=9, fontweight="bold" if mark else "normal")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        # ---- (b) heatmap Healthy->F4 count -------------------------------
        ax = axes[0, 1]
        hf4_mat = grid.pivot(index="f4_cut", columns="dom_cut",
                             values="n_healthy_to_F4").reindex(
            index=F4_CUTS, columns=DOM_CUTS)
        im = ax.imshow(hf4_mat.values, cmap="Reds", aspect="auto",
                       vmin=0, vmax=max(np.nanmax(hf4_mat.values), 1))
        ax.set_xticks(range(len(DOM_CUTS)))
        ax.set_xticklabels([f">={c}" for c in DOM_CUTS])
        ax.set_yticks(range(len(F4_CUTS)))
        ax.set_yticklabels([f"{c:.2f}" for c in F4_CUTS])
        ax.set_xlabel("dominant_fstage cutoff")
        ax.set_ylabel("F4_frac cutoff")
        ax.set_title(r"(b) Clean-healthy Healthy$\to$F4 hallucinations (LOOCV)")
        for i, f4 in enumerate(F4_CUTS):
            for j, dc in enumerate(DOM_CUTS):
                v = hf4_mat.iloc[i, j]
                mark = " *" if (f4 == 0.40 and dc == 3) else ""
                ax.text(j, i, f"{int(v)}{mark}", ha="center", va="center",
                        color="black", fontsize=9,
                        fontweight="bold" if mark else "normal")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

        # ---- (c) line plot rho vs F4 cut, one line per dom cut ------------
        ax = axes[1, 0]
        colors = {2: "#7A0E0E", 3: "#2E8B57", 4: "#5B7FB1"}
        for dc in DOM_CUTS:
            sub = grid[grid["dom_cut"] == dc].sort_values("f4_cut")
            ax.plot(sub["f4_cut"], sub["external_rho"],
                    "o-", color=colors[dc], lw=2,
                    label=f"dominant_fstage >= {dc}"
                          + ("  (current)" if dc == 3 else ""))
        ax.axvline(0.40, ls="--", color="gray", lw=0.8,
                   alpha=0.6, label="current F4 cut=0.40")
        ax.set_xlabel("F4_frac cutoff")
        ax.set_ylabel(r"External Spearman $\rho$")
        ax.set_title(r"(c) $\rho$ across F4_frac cutoff, by dominant_fstage cutoff")
        ax.legend(fontsize=8, loc="best")
        ax.grid(alpha=0.25)

        # ---- (d) Healthy->F4 vs F4 cut ----------------------------------
        ax = axes[1, 1]
        for dc in DOM_CUTS:
            sub = grid[grid["dom_cut"] == dc].sort_values("f4_cut")
            ax.plot(sub["f4_cut"], sub["n_healthy_to_F4"],
                    "o-", color=colors[dc], lw=2,
                    label=f"dominant_fstage >= {dc}"
                          + ("  (current)" if dc == 3 else ""))
        ax.axvline(0.40, ls="--", color="gray", lw=0.8, alpha=0.6,
                   label="current F4 cut=0.40")
        ax.set_xlabel("F4_frac cutoff")
        ax.set_ylabel(r"clean-healthy Healthy$\to$F4 (LOOCV count)")
        ax.set_title(r"(d) Hallucinations across thresholds")
        ax.legend(fontsize=8, loc="best")
        ax.grid(alpha=0.25)

        fig.suptitle(
            "Dubious-healthy threshold sensitivity (Script 343s)\n"
            "augmented scVI ordinal logistic; * marks current default (F4>0.40 | dom>=3)",
            fontsize=12, fontweight="bold")
        fig.tight_layout(rect=[0, 0, 1, 0.95])
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"[fig] wrote {OUT_PDF}")


# ----------------------------------------------------------- main
def main():
    print("=" * 70)
    print("343s: dubious-healthy threshold sensitivity sweep")
    print(f"HAVE_MORD={HAVE_MORD}")
    print(f"grid: F4_CUTS={F4_CUTS}  DOM_CUTS={DOM_CUTS}")
    print("=" * 70)

    df = load_inputs()
    print(f"[load] {len(df)} donors total")
    print(f"[load] {int(df[LATENT_COLS].notna().all(axis=1).sum())} "
          f"donors with scVI latents")
    print(f"[load] disease_stage_coarse:")
    print(df["disease_stage_coarse"].value_counts().to_string())

    rows = []
    for f4 in F4_CUTS:
        for dc in DOM_CUTS:
            print(f"\n>>> cell: F4_cut={f4:.2f}  dom_cut>={dc}")
            r = evaluate_cell(df, f4, dc)
            print(f"    n_dubious={r['n_dubious']}  "
                  f"n_clean_train={r['n_clean_train']}  "
                  f"n_train={r['n_train']}  "
                  f"rho={r['external_rho']}  "
                  f"H->F4={r['n_healthy_to_F4']}/{r['n_chealthy_loo']}")
            rows.append(r)

    grid = pd.DataFrame(rows)
    grid["is_current"] = ((grid["f4_cut"] == 0.40) & (grid["dom_cut"] == 3))
    grid.to_csv(OUT_TSV, sep="\t", index=False)
    print(f"\n[output] {OUT_TSV}  ({len(grid)} rows)")

    print("\n" + "=" * 70)
    print("FULL GRID")
    print("=" * 70)
    print(grid[["f4_cut", "dom_cut", "n_dubious", "n_clean_train",
                "n_train", "external_rho", "n_healthy_to_F4",
                "qwk_andrews_loocv", "qwk_full_loocv",
                "sh_F3_frac", "sh_F4_frac",
                "is_current"]].to_string(index=False))

    # ---- verdict --------------------------------------------------------
    cur = grid[grid["is_current"]].iloc[0]
    cur_rho = float(cur["external_rho"])
    best = grid.loc[grid["external_rho"].idxmax()]
    best_rho = float(best["external_rho"])
    delta = best_rho - cur_rho

    print("\n" + "=" * 70)
    print("VERDICT")
    print("=" * 70)
    print(f"current default:  F4_cut=0.40, dom_cut>=3  "
          f"-> rho={cur_rho:.4f}, H->F4={int(cur['n_healthy_to_F4'])}/"
          f"{int(cur['n_chealthy_loo'])}")
    print(f"best in grid:     F4_cut={best['f4_cut']:.2f}, "
          f"dom_cut>={int(best['dom_cut'])}  -> rho={best_rho:.4f}, "
          f"H->F4={int(best['n_healthy_to_F4'])}/"
          f"{int(best['n_chealthy_loo'])}")
    print(f"delta rho:        {delta:+.4f}")
    if delta > 0.05:
        print(f">> RECOMMEND SWITCHING to F4_cut={best['f4_cut']:.2f}, "
              f"dom_cut>={int(best['dom_cut'])} "
              f"(delta rho={delta:+.4f} > 0.05)")
    else:
        print(f">> DEFENSIBLE: current default within 0.05 rho of grid max "
              f"(delta={delta:+.4f}); stop tuning.")
    print("=" * 70)

    make_figure(grid)


if __name__ == "__main__":
    main()
