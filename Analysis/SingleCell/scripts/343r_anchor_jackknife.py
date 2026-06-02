#!/usr/bin/env python
"""
343r_anchor_jackknife.py

Anchor-donor leave-one-out jackknife stability for the MASLD augmented
F-stage classifier (companion to 343m_augmented_fstage.py).

Question
--------
Is the augmented training set's effect on SH predictions driven by a small
number of influential anchor donors, or is it broadly supported?

Procedure
---------
For each of the 127 augmented training donors
  (58 Andrews documented + 60 clean-healthy + 9 cirrhosis)
do LOO:
  1. Train augmented scVI ordinal logistic on the remaining 126 anchors
  2. Predict F-stage probabilities (and argmax) for all 101 SH donors
  3. Record per-SH-donor delta vs full-training prediction

Metrics
-------
  - Per-anchor influence = mean |delta argmax F| across 101 SH donors
  - Per-anchor probability influence = mean L1 distance between P(F) vectors
  - Top-5 most influential anchors, and per-class concentration
  - Per-class total + per-donor mean influence (class-size-normalised)
  - Decision: top-5 share of total influence
       <30%  -> broadly supported (anchor-robust)
       >60%  -> over-reliant on a few anchors (caveat flag)

Outputs
-------
  - anchor_jackknife_influence.tsv          (127 rows x metrics)
  - anchor_jackknife_top_influencers.tsv    (top-10 rows)
  - figS_anchor_jackknife.pdf               (ranked + per-class breakdown)
"""

from __future__ import annotations
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

from joblib import Parallel, delayed

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

try:
    import mord
    HAVE_MORD = True
except ImportError:
    HAVE_MORD = False


# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

STAGE_DIR  = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
FSTAGE_DOC = STAGE_DIR / "donor_fstage_documented.tsv"
DUBIOUS    = STAGE_DIR / "dubious_healthy_donors.tsv"
DONOR_META = STAGE_DIR / "donor_metadata_extended.tsv"
HEP_PB     = STAGE_DIR / "donor_hep_scvi_mean.tsv"

OUT_INFL   = STAGE_DIR / "anchor_jackknife_influence.tsv"
OUT_TOP    = STAGE_DIR / "anchor_jackknife_top_influencers.tsv"
OUT_PDF    = PROJECT_ROOT / "figures/supplementary/stage_ccc/figS_anchor_jackknife.pdf"

LATENT_COLS = [f"z{i}" for i in range(20)]

N_JOBS = int(os.environ.get("SLURM_CPUS_PER_TASK", 16))

PAL = {
    "documented":           "#2E8B57",
    "clean_healthy_anchor": "#5B7FB1",
    "cirrhosis_anchor":     "#7A0E0E",
}


# ----------------------------------------------------------------------
# Data assembly  (mirrors 343m)
# ----------------------------------------------------------------------
def load_donors() -> pd.DataFrame:
    meta = pd.read_csv(DONOR_META, sep="\t")
    doc  = pd.read_csv(FSTAGE_DOC,  sep="\t")
    doc["F_stage_documented"] = pd.to_numeric(doc["F_stage_documented"], errors="coerce")
    pb   = pd.read_csv(HEP_PB,     sep="\t")
    dub  = pd.read_csv(DUBIOUS,    sep="\t")

    df = meta[["sample", "dataset", "disease_stage_coarse",
               "disease_stage_numeric"]].copy()
    df = df.merge(pb[["sample", "dataset"] + LATENT_COLS], on=["sample", "dataset"], how="left")
    df = df.merge(doc[["sample", "dataset", "F_stage_documented"]],
                  on=["sample", "dataset"], how="left")
    df["dubious"] = df["sample"].isin(dub["sample"].values)
    return df


def assemble_training(df: pd.DataFrame) -> pd.DataFrame:
    """Same logic as 343m.assemble_training -> 127 anchors expected."""
    train = df.copy()
    train["y_train"] = np.nan
    train["origin"]  = "none"

    mask_doc = train["F_stage_documented"].notna()
    train.loc[mask_doc, "y_train"] = train.loc[mask_doc, "F_stage_documented"]
    train.loc[mask_doc, "origin"]  = "documented"

    mask_chealthy = ((train["disease_stage_coarse"] == "Healthy") &
                     (~train["dubious"]) & (~mask_doc))
    train.loc[mask_chealthy, "y_train"] = 0
    train.loc[mask_chealthy, "origin"]  = "clean_healthy_anchor"

    mask_cirr = ((train["disease_stage_coarse"] == "Cirrhosis") & (~mask_doc))
    train.loc[mask_cirr, "y_train"] = 4
    train.loc[mask_cirr, "origin"]  = "cirrhosis_anchor"

    train_set = train[train["y_train"].notna() &
                      train[LATENT_COLS].notna().all(axis=1)].copy()
    train_set["y_train"] = train_set["y_train"].astype(int)
    return train_set


# ----------------------------------------------------------------------
# Classifier  (matches 343m.make_ordinal: scVI -> StandardScaler -> mord.LogisticIT)
# ----------------------------------------------------------------------
def make_ordinal():
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


def predict_proba_full(clf, X):
    """Return (n, 5) probability matrix on F0..F4 ordering, padding absent classes with zeros.

    mord.LogisticIT exposes predict_proba even when not all 5 classes are present
    in y; class ordering follows clf.classes_.
    """
    proba = clf.predict_proba(X)
    classes = clf.named_steps["clf"].classes_
    full = np.zeros((proba.shape[0], 5), dtype=float)
    for j, c in enumerate(classes):
        full[:, int(c)] = proba[:, j]
    rowsum = full.sum(axis=1, keepdims=True)
    rowsum[rowsum == 0] = 1.0
    return full / rowsum


# ----------------------------------------------------------------------
# Jackknife core
# ----------------------------------------------------------------------
def fit_and_predict_sh(X_train: np.ndarray, y_train: np.ndarray,
                       X_sh: np.ndarray):
    clf = make_ordinal()
    clf.fit(X_train, y_train)
    yhat = predict_argmax(clf, X_sh).astype(int)
    proba = predict_proba_full(clf, X_sh)
    return yhat, proba


def loo_one(i: int, X_anchor: np.ndarray, y_anchor: np.ndarray,
            X_sh: np.ndarray):
    mask = np.ones(len(y_anchor), dtype=bool)
    mask[i] = False
    yhat, proba = fit_and_predict_sh(X_anchor[mask], y_anchor[mask], X_sh)
    return i, yhat, proba


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    print(f"[init] HAVE_MORD={HAVE_MORD}  N_JOBS={N_JOBS}")
    # NOTE: In production 343m the rapids_singlecell env does NOT ship mord, so
    # `make_ordinal` falls back to class-balanced LogisticRegression. The
    # jackknife mirrors that exact behaviour by reusing `make_ordinal` here,
    # so HAVE_MORD=False is the production-equivalent path.

    t0 = time.time()

    df_all = load_donors()
    print(f"[load] donor roster: {len(df_all)} donors")

    train_set = assemble_training(df_all)
    print(f"[load] augmented training set: {len(train_set)} anchors")
    print(train_set.groupby(["origin", "y_train"]).size().unstack(fill_value=0))

    # SH targets for prediction
    sh_mask = (df_all["disease_stage_coarse"] == "Steatohepatitis") & \
              (df_all[LATENT_COLS].notna().all(axis=1))
    df_sh = df_all.loc[sh_mask, ["sample", "dataset"] + LATENT_COLS].reset_index(drop=True)
    print(f"[load] SH donors with latents: {len(df_sh)}")

    # Matrices
    X_anchor = train_set[LATENT_COLS].to_numpy(dtype=float)
    y_anchor = train_set["y_train"].to_numpy(dtype=int)
    X_sh     = df_sh[LATENT_COLS].to_numpy(dtype=float)
    n_anchor = len(y_anchor)
    n_sh     = len(df_sh)

    anchor_meta = train_set[["sample", "dataset", "origin", "y_train",
                              "disease_stage_coarse"]].reset_index(drop=True)
    anchor_meta = anchor_meta.rename(columns={"y_train": "anchor_F"})

    # ------------------------------------------------------------------
    # Full-training reference predictions on SH
    # ------------------------------------------------------------------
    print("[ref] fitting on full 127 anchors...")
    y_full, proba_full = fit_and_predict_sh(X_anchor, y_anchor, X_sh)
    print(f"[ref] SH F-distribution: {np.bincount(y_full, minlength=5).tolist()}")

    # ------------------------------------------------------------------
    # Parallel LOO fits
    # ------------------------------------------------------------------
    print(f"[loo] running 127 LOO fits in parallel (n_jobs={N_JOBS})...")
    results = Parallel(n_jobs=N_JOBS, backend="loky", verbose=5)(
        delayed(loo_one)(i, X_anchor, y_anchor, X_sh) for i in range(n_anchor)
    )

    # Aggregate
    yhat_loo  = np.zeros((n_anchor, n_sh), dtype=int)
    proba_loo = np.zeros((n_anchor, n_sh, 5), dtype=float)
    for i, yhat, proba in results:
        yhat_loo[i]  = yhat
        proba_loo[i] = proba

    # ------------------------------------------------------------------
    # Per-anchor influence
    # ------------------------------------------------------------------
    delta_arg = yhat_loo - y_full[None, :]                            # (127, 101)
    abs_delta_arg = np.abs(delta_arg).astype(float)
    influence_arg = abs_delta_arg.mean(axis=1)                        # (127,)

    delta_proba = proba_loo - proba_full[None, :, :]                  # (127, 101, 5)
    l1_proba = np.abs(delta_proba).sum(axis=2)                        # (127, 101)
    influence_prob_l1 = l1_proba.mean(axis=1)                         # (127,)

    n_flipped = (delta_arg != 0).sum(axis=1).astype(int)
    n_flipped_up = (delta_arg > 0).sum(axis=1).astype(int)
    n_flipped_down = (delta_arg < 0).sum(axis=1).astype(int)
    max_abs_delta = np.abs(delta_arg).max(axis=1).astype(int)

    infl_df = anchor_meta.copy()
    infl_df["influence_argmax_mean_abs"] = influence_arg
    infl_df["influence_proba_l1_mean"]   = influence_prob_l1
    infl_df["n_sh_flipped"]              = n_flipped
    infl_df["n_sh_flipped_higher"]       = n_flipped_up
    infl_df["n_sh_flipped_lower"]        = n_flipped_down
    infl_df["max_abs_delta_F"]           = max_abs_delta

    # Rank by argmax influence
    infl_df = infl_df.sort_values("influence_argmax_mean_abs",
                                  ascending=False).reset_index(drop=True)
    infl_df.insert(0, "rank", np.arange(1, len(infl_df) + 1))

    infl_df.to_csv(OUT_INFL, sep="\t", index=False)
    print(f"[out] {OUT_INFL}")

    # ------------------------------------------------------------------
    # Top-N + decision rule
    # ------------------------------------------------------------------
    top10 = infl_df.head(10).copy()
    top10.to_csv(OUT_TOP, sep="\t", index=False)
    print(f"[out] {OUT_TOP}")

    total_infl = infl_df["influence_argmax_mean_abs"].sum()
    if total_infl > 0:
        top5_share = (infl_df.head(5)["influence_argmax_mean_abs"].sum()
                      / total_infl)
    else:
        top5_share = 0.0

    if top5_share < 0.30:
        verdict = "ANCHOR-ROBUST"
        verdict_desc = ("Top-5 anchors account for <30% of total influence "
                        "-> broadly supported.")
    elif top5_share > 0.60:
        verdict = "ANCHOR-DEPENDENT"
        verdict_desc = ("Top-5 anchors account for >60% of total influence "
                        "-> flag for sensitivity caveat.")
    else:
        verdict = "INTERMEDIATE"
        verdict_desc = (f"Top-5 anchors account for {top5_share*100:.1f}% "
                        "of total influence -> intermediate (30-60%).")

    # Per-class breakdown
    class_summary = infl_df.groupby("origin").agg(
        n=("sample", "size"),
        total_influence=("influence_argmax_mean_abs", "sum"),
        mean_influence=("influence_argmax_mean_abs", "mean"),
        median_influence=("influence_argmax_mean_abs", "median"),
        max_influence=("influence_argmax_mean_abs", "max"),
    ).reset_index()
    class_summary["share_of_total"] = class_summary["total_influence"] / total_infl
    class_summary = class_summary.sort_values("mean_influence", ascending=False)

    print("\n" + "=" * 72)
    print("ANCHOR JACKKNIFE SUMMARY")
    print("=" * 72)
    print(f"n_anchors        : {n_anchor}")
    print(f"n_sh_targets     : {n_sh}")
    print(f"mean influence   : {influence_arg.mean():.4f} F-units / SH donor")
    print(f"max influence    : {influence_arg.max():.4f} "
          f"(anchor: {infl_df.iloc[0]['sample']} / "
          f"{infl_df.iloc[0]['dataset']} / {infl_df.iloc[0]['origin']})")
    print(f"total influence  : {total_infl:.4f}")
    print(f"top-5 share      : {top5_share*100:.1f}%")
    print(f"VERDICT          : {verdict}")
    print(f"                   {verdict_desc}")
    print("\nPer-class summary:")
    print(class_summary.to_string(index=False))
    print("\nTop-5 most influential anchors:")
    print(infl_df.head(5)[["rank", "sample", "dataset", "origin",
                          "anchor_F", "influence_argmax_mean_abs",
                          "n_sh_flipped"]].to_string(index=False))
    print("=" * 72)

    # ------------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------------
    make_figure(infl_df, class_summary, total_infl, top5_share,
                verdict, verdict_desc, n_anchor, n_sh)

    print(f"\n[done] elapsed {time.time()-t0:.1f}s")


# ----------------------------------------------------------------------
def make_figure(infl_df: pd.DataFrame, class_summary: pd.DataFrame,
                total_infl: float, top5_share: float,
                verdict: str, verdict_desc: str,
                n_anchor: int, n_sh: int):
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)

    with PdfPages(OUT_PDF) as pdf:
        fig = plt.figure(figsize=(14, 11))
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 1],
                              hspace=0.42, wspace=0.32)

        # ----- Panel A: ranked per-anchor influence -----
        axA = fig.add_subplot(gs[0, :])
        x = np.arange(len(infl_df))
        colors = [PAL.get(o, "#888888") for o in infl_df["origin"]]
        axA.bar(x, infl_df["influence_argmax_mean_abs"], color=colors,
                edgecolor="black", linewidth=0.2)
        axA.set_xlabel("Anchor rank (sorted by influence)")
        axA.set_ylabel("Mean |Δ F| over 101 SH donors")
        axA.set_title("(a) Per-anchor LOO influence on SH F-stage predictions  "
                      f"(n={n_anchor} anchors -> n={n_sh} SH targets)")
        axA.axhline(0, color="black", lw=0.4)
        axA.grid(axis="y", alpha=0.25)

        # Annotate top-5
        for k in range(min(5, len(infl_df))):
            row = infl_df.iloc[k]
            axA.annotate(
                f"{row['sample']}\n{row['origin'].split('_')[0]} F={int(row['anchor_F'])}",
                xy=(k, row["influence_argmax_mean_abs"]),
                xytext=(k, row["influence_argmax_mean_abs"] + 0.03),
                ha="left", fontsize=7, rotation=90,
                arrowprops=dict(arrowstyle="-", color="black", lw=0.4),
            )

        # Legend
        from matplotlib.patches import Patch
        legend_elements = [Patch(facecolor=PAL[k], edgecolor="black", label=k)
                           for k in ["documented", "clean_healthy_anchor",
                                     "cirrhosis_anchor"] if k in PAL]
        axA.legend(handles=legend_elements, loc="upper right", fontsize=9)

        # ----- Panel B: cumulative share + top-5 mark -----
        axB = fig.add_subplot(gs[1, 0])
        cum = infl_df["influence_argmax_mean_abs"].cumsum() / total_infl
        axB.plot(np.arange(1, len(cum) + 1), cum * 100, color="#333333", lw=1.5)
        axB.axvline(5, color="red", ls="--", lw=1.0, alpha=0.7,
                    label=f"top-5 = {top5_share*100:.1f}%")
        axB.axhline(30, color="green", ls=":", lw=0.7, alpha=0.6,
                    label="robust (<30%)")
        axB.axhline(60, color="darkred", ls=":", lw=0.7, alpha=0.6,
                    label="caveat (>60%)")
        axB.set_xlabel("Number of top-ranked anchors")
        axB.set_ylabel("Cumulative share of total influence (%)")
        axB.set_title("(b) Concentration curve")
        axB.set_xlim(0, len(cum))
        axB.set_ylim(0, 100)
        axB.legend(fontsize=8, loc="lower right")
        axB.grid(alpha=0.3)

        # Verdict box
        axB.text(0.02, 0.97,
                 f"VERDICT: {verdict}\n{verdict_desc}",
                 transform=axB.transAxes,
                 va="top", ha="left", fontsize=9,
                 bbox=dict(boxstyle="round,pad=0.4",
                           fc="#F8F8F8", ec="black", lw=0.6))

        # ----- Panel C: per-class mean influence (normalised) -----
        axC = fig.add_subplot(gs[1, 1])
        cs = class_summary.copy()
        x = np.arange(len(cs))
        bars = axC.bar(x, cs["mean_influence"],
                       color=[PAL.get(o, "#888888") for o in cs["origin"]],
                       edgecolor="black", lw=0.4)
        axC.set_xticks(x)
        axC.set_xticklabels([f"{o}\n(n={int(n)})"
                             for o, n in zip(cs["origin"], cs["n"])],
                            fontsize=8)
        axC.set_ylabel("Mean per-anchor influence")
        axC.set_title("(c) Per-class mean influence (size-normalised)")
        axC.grid(axis="y", alpha=0.3)

        # Annotate share-of-total
        for b, share in zip(bars, cs["share_of_total"]):
            axC.text(b.get_x() + b.get_width() / 2.0,
                     b.get_height() + 0.001,
                     f"{share*100:.1f}%\nof total",
                     ha="center", va="bottom", fontsize=8)

        fig.suptitle(
            "Anchor-donor leave-one-out jackknife of augmented F-stage classifier\n"
            f"(127 anchors x 101 SH donors; top-5 share {top5_share*100:.1f}%)",
            fontsize=12, fontweight="bold", y=0.995)

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

    print(f"[fig] wrote {OUT_PDF}")


if __name__ == "__main__":
    main()
