#!/usr/bin/env python
"""
343q_bootstrap_stability.py

Bootstrap-stability gate for the MASLD F-stage augmented classifier (Script 343m).

Decision rule (pre-registered)
------------------------------
Adversarial synthesizer (brainstorm_adversarial_synthesis.md) concluded the
augmented scVI ordinal method (343m) is manuscript-ready ONLY if:
    (a) >= 80% of Steatohepatitis donors are stable within +/- 1 F-stage under
        B=500 stratified bootstrap resampling of the training set, AND
    (b) the SH F-stage distribution mode (F3) is preserved in >= 95% of
        bootstraps.

If either fails -> FALL_BACK to scVI ordinal vanilla (rho = 0.391, 343b).

Methods compared (B=500 each, same RNG seed across methods)
-----------------------------------------------------------
1. scVI ordinal AUGMENTED  (LogisticRegression, balanced;
                            train = 58 Andrews + clean-healthy F0 + cirrhosis F4)
2. scVI ordinal VANILLA    (same classifier; train = 58 Andrews only)
3. scANVI vanilla          (donor-level posterior bootstrap; multinomial sample
                            of the 5-class posterior per donor; no retraining)
4. scANVI augmented        (donor-level posterior bootstrap, same approach)

Outputs
-------
  - bootstrap_stability_per_donor.tsv     (101 SH x 4 methods x stability metrics)
  - bootstrap_sh_distribution_stability.tsv (5 stages x B for each method)
  - bootstrap_decision.txt                (SHIP or FALL_BACK with numbers)
  - figures/.../figS_bootstrap_stability.pdf  (3-panel)
  - BOOTSTRAP_FALLBACK_REQUIRED.flag       (only if FALL_BACK)

Constraints
-----------
- DO NOT modify any existing scripts.
- Same bootstrap RNG seed across all 4 methods => comparable resamples.
- mord is not installed in rapids_singlecell; we use LogisticRegression
  (class_weight='balanced'), matching the existing 343m fallback path.
"""

from __future__ import annotations
import os
import sys
import json
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

# ----------------------------------------------------------------------
# Paths
# ----------------------------------------------------------------------
PROJECT_ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))
STAGE_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"

FSTAGE_DOC          = STAGE_DIR / "donor_fstage_documented.tsv"
DUBIOUS             = STAGE_DIR / "dubious_healthy_donors.tsv"
DONOR_META          = STAGE_DIR / "donor_metadata_extended.tsv"
HEP_PB              = STAGE_DIR / "donor_hep_scvi_mean.tsv"
SH_PRED             = STAGE_DIR / "donor_fstage_augmented_predicted.tsv"
SCANVI_VAN          = STAGE_DIR / "donor_fstage_scanvi_predicted.tsv"
SCANVI_AUG          = STAGE_DIR / "donor_fstage_scanvi_augmented_predicted.tsv"

OUT_PER_DONOR       = STAGE_DIR / "bootstrap_stability_per_donor.tsv"
OUT_SH_DIST         = STAGE_DIR / "bootstrap_sh_distribution_stability.tsv"
OUT_DECISION        = STAGE_DIR / "bootstrap_decision.txt"
OUT_FLAG            = STAGE_DIR / "BOOTSTRAP_FALLBACK_REQUIRED.flag"
OUT_PDF             = PROJECT_ROOT / "figures/supplementary/stage_ccc/figS_bootstrap_stability.pdf"

LATENT_COLS = [f"z{i}" for i in range(20)]

B = 500
RNG_SEED = 42

# Method labels
M_AUG_ORD = "scvi_ord_augmented"
M_VAN_ORD = "scvi_ord_vanilla"
M_SCAN_V  = "scanvi_vanilla"
M_SCAN_A  = "scanvi_augmented"
METHODS = [M_AUG_ORD, M_VAN_ORD, M_SCAN_V, M_SCAN_A]


# ----------------------------------------------------------------------
# Data assembly (mirrors 343m logic)
# ----------------------------------------------------------------------
def load_data():
    meta = pd.read_csv(DONOR_META, sep="\t")
    doc  = pd.read_csv(FSTAGE_DOC, sep="\t")
    doc["F_stage_documented"] = pd.to_numeric(doc["F_stage_documented"],
                                              errors="coerce")
    pb   = pd.read_csv(HEP_PB, sep="\t")
    dub  = pd.read_csv(DUBIOUS, sep="\t")
    sh_tbl = pd.read_csv(SH_PRED, sep="\t")

    df = meta[["sample", "dataset", "disease_stage_coarse",
               "disease_stage_numeric"]].copy()
    df = df.merge(pb[["sample", "dataset"] + LATENT_COLS],
                  on=["sample", "dataset"], how="left")
    df = df.merge(doc[["sample", "dataset", "F_stage_documented"]],
                  on=["sample", "dataset"], how="left")
    df["dubious"] = df["sample"].isin(dub["sample"].values)
    return df, sh_tbl


def assemble_training(df: pd.DataFrame, augmented: bool) -> pd.DataFrame:
    train = df.copy()
    train["y_train"] = np.nan
    train["origin"]  = "none"

    # 1. Documented Andrews (used in both vanilla + augmented)
    mask_doc = train["F_stage_documented"].notna()
    train.loc[mask_doc, "y_train"] = train.loc[mask_doc, "F_stage_documented"]
    train.loc[mask_doc, "origin"]  = "documented"

    if augmented:
        # 2. Clean-healthy => F0 (where not already documented)
        mask_chealthy = ((train["disease_stage_coarse"] == "Healthy") &
                         (~train["dubious"]) &
                         (~mask_doc))
        train.loc[mask_chealthy, "y_train"] = 0
        train.loc[mask_chealthy, "origin"]  = "clean_healthy_anchor"

        # 3. Cirrhosis => F4
        mask_cirr = ((train["disease_stage_coarse"] == "Cirrhosis") &
                     (~mask_doc))
        train.loc[mask_cirr, "y_train"] = 4
        train.loc[mask_cirr, "origin"]  = "cirrhosis_anchor"

    train_set = train[train["y_train"].notna() &
                      train[LATENT_COLS].notna().all(axis=1)].copy()
    train_set["y_train"] = train_set["y_train"].astype(int)
    return train_set


def make_ordinal():
    """Match 343m fallback: scaled multinomial LR with balanced class weight."""
    return Pipeline([("scale", StandardScaler()),
                     ("clf",   LogisticRegression(solver="lbfgs",
                                                   C=1.0,
                                                   max_iter=2000,
                                                   class_weight="balanced"))])


def predict_argmax_lr(clf, X):
    proba = clf.predict_proba(X)
    return clf.named_steps["clf"].classes_[proba.argmax(axis=1)].astype(int)


# ----------------------------------------------------------------------
# Bootstrap routines
# ----------------------------------------------------------------------
def stratified_bootstrap_indices(y: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Stratified sample-with-replacement WITHIN each class. Preserves class sizes."""
    idx_out = []
    for c in np.unique(y):
        cls_idx = np.where(y == c)[0]
        n_c = len(cls_idx)
        sampled = rng.choice(cls_idx, size=n_c, replace=True)
        idx_out.append(sampled)
    return np.concatenate(idx_out)


def boot_one_scvi(b: int, X_tr: np.ndarray, y_tr: np.ndarray,
                  X_sh: np.ndarray, seed_base: int) -> np.ndarray:
    """One bootstrap iter for scVI ordinal: refit on resampled training set,
    predict SH donors. Returns 1D int array of length n_sh."""
    rng = np.random.default_rng(seed_base + b)
    boot_idx = stratified_bootstrap_indices(y_tr, rng)
    Xb = X_tr[boot_idx]
    yb = y_tr[boot_idx]
    # If a bootstrap resample drops a class (rare given stratified by class),
    # fall back to fewer classes; predict_argmax still picks the best of those.
    clf = make_ordinal()
    clf.fit(Xb, yb)
    return predict_argmax_lr(clf, X_sh)


def boot_one_scanvi(b: int, posterior_matrix: np.ndarray,
                    seed_base: int) -> np.ndarray:
    """One bootstrap iter for scANVI: per-donor multinomial sample from posterior.
    posterior_matrix: (n_donors, 5). Returns 1D int array of length n_donors."""
    rng = np.random.default_rng(seed_base + b)
    n = posterior_matrix.shape[0]
    out = np.empty(n, dtype=int)
    for i in range(n):
        # safety: if row sums to 0 (no cells / NaN), force F0 fallback
        p = posterior_matrix[i]
        if not np.all(np.isfinite(p)) or p.sum() <= 0:
            out[i] = 0
            continue
        p = p / p.sum()
        out[i] = rng.choice(5, p=p)
    return out


# ----------------------------------------------------------------------
# Stability metrics
# ----------------------------------------------------------------------
def per_donor_stability(boot_preds: np.ndarray,
                         original_pred: np.ndarray) -> pd.DataFrame:
    """
    boot_preds: (B, n_donors) int matrix of predicted F-stages per iter.
    original_pred: (n_donors,) the original (non-bootstrap) F-stage prediction.
    """
    B_, n = boot_preds.shape
    rows = []
    for i in range(n):
        preds_i = boot_preds[:, i]
        counts = np.bincount(preds_i, minlength=5)
        mode_f = int(counts.argmax())
        within_1_pct = float(np.mean(np.abs(preds_i - mode_f) <= 1)) * 100.0
        q25, q75 = np.percentile(preds_i, [25, 75])
        majority_match = float(np.mean(preds_i == original_pred[i])) * 100.0
        rows.append(dict(
            mode_fstage=mode_f,
            within_1_pct=within_1_pct,
            q25=float(q25),
            q75=float(q75),
            q25_q75_range=float(q75 - q25),
            bootstrap_majority_pct=majority_match,
            std_pred=float(np.std(preds_i)),
        ))
    return pd.DataFrame(rows)


def sh_distribution_per_boot(boot_preds: np.ndarray) -> np.ndarray:
    """Return (5, B) matrix: F-stage fractions across SH donors per iter."""
    B_, n = boot_preds.shape
    out = np.zeros((5, B_), dtype=float)
    for b in range(B_):
        c = np.bincount(boot_preds[b], minlength=5)
        out[:, b] = c / c.sum()
    return out


# ----------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------
MASLD_PAL = {
    M_AUG_ORD: "#2E8B57",
    M_VAN_ORD: "#5B7FB1",
    M_SCAN_V:  "#9B59B6",
    M_SCAN_A:  "#E66F51",
}
FSTAGE_PAL = ["#9E9E9E", "#FFD16B", "#FFB454", "#E66F51", "#7A0E0E"]


def make_figure(stability: dict[str, pd.DataFrame],
                sh_dist: dict[str, np.ndarray],
                decision: dict):
    OUT_PDF.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(OUT_PDF) as pdf:
        fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))

        # Panel A: per-donor within_1_pct distribution per method
        ax = axes[0]
        for m in METHODS:
            vals = stability[m]["within_1_pct"].values
            ax.hist(vals, bins=np.arange(0, 105, 5),
                    histtype="step", linewidth=2.0,
                    color=MASLD_PAL[m], label=m)
        ax.axvline(80, ls="--", lw=1.0, color="black", alpha=0.6,
                   label="80% threshold")
        ax.set_xlabel(r"Per-donor within $\pm$1 stage (%)")
        ax.set_ylabel("Number of SH donors")
        ax.set_title("(a) Per-donor bootstrap stability")
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(axis="y", alpha=0.25)

        # Panel B: SH F-stage distribution 95% CI ribbons for the gating method
        ax = axes[1]
        x = np.arange(5)
        w = 0.20
        for j, m in enumerate(METHODS):
            d = sh_dist[m]  # (5, B)
            lo = np.percentile(d, 2.5,  axis=1)
            hi = np.percentile(d, 97.5, axis=1)
            mid = np.percentile(d, 50.0, axis=1)
            xj = x + (j - 1.5) * w
            ax.bar(xj, mid, w, color=MASLD_PAL[m], alpha=0.85, label=m,
                   yerr=[mid - lo, hi - mid], capsize=2, ecolor="black",
                   error_kw={"lw": 0.6})
        ax.set_xticks(x)
        ax.set_xticklabels(["F0", "F1", "F2", "F3", "F4"])
        ax.set_ylabel("Fraction of SH donors (median, 95% CI)")
        ax.set_title("(b) SH F-stage distribution stability (B=500)")
        ax.legend(fontsize=8, loc="upper left")
        ax.grid(axis="y", alpha=0.25)

        # Panel C: method comparison bar chart of headline metrics
        ax = axes[2]
        metrics = {
            r"median within$\pm$1 (%)": [decision[m]["median_within_1"] for m in METHODS],
            r"% donors $\geq$80% stable":  [decision[m]["pct_donors_ge_80"] for m in METHODS],
            "% boots F3-mode preserved": [decision[m]["pct_boots_F3_mode"] for m in METHODS],
        }
        labels = list(metrics.keys())
        x = np.arange(len(labels))
        w = 0.20
        for j, m in enumerate(METHODS):
            vals = [metrics[k][j] for k in labels]
            ax.bar(x + (j - 1.5) * w, vals, w, color=MASLD_PAL[m], label=m)
        ax.axhline(80, ls="--", lw=1.0, color="black", alpha=0.6)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=8, rotation=12, ha="right")
        ax.set_ylabel("Percent")
        ax.set_ylim(0, 105)
        ax.set_title("(c) Method comparison (bootstrap headline stats)")
        ax.legend(fontsize=8, loc="lower right")
        ax.grid(axis="y", alpha=0.25)

        # Decision banner
        verdict = decision["verdict"]
        banner_color = "#2E8B57" if verdict == "SHIP" else "#7A0E0E"
        fig.suptitle(
            f"Bootstrap-stability gate  |  B={B}  |  "
            f"VERDICT: {verdict}  |  "
            f"AUG: {decision[M_AUG_ORD]['pct_donors_ge_80']:.1f}% donors ge 80%, "
            f"{decision[M_AUG_ORD]['pct_boots_F3_mode']:.1f}% boots F3-mode",
            fontsize=12, color=banner_color, y=1.00, fontweight="bold",
        )
        plt.tight_layout()
        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)
    print(f"[fig] wrote {OUT_PDF}")


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():
    print(f"[start] 343q bootstrap stability  |  B={B}  |  seed={RNG_SEED}")

    df, sh_tbl = load_data()

    # SH targets: 101 donors with disease_stage_coarse == Steatohepatitis
    sh_df = sh_tbl[sh_tbl["disease_stage_coarse"] == "Steatohepatitis"].copy()
    # require latents
    sh_df = sh_df.merge(df[["sample", "dataset"] + LATENT_COLS],
                        on=["sample", "dataset"], how="left")
    sh_df = sh_df[sh_df[LATENT_COLS].notna().all(axis=1)].reset_index(drop=True)
    sh_keys = sh_df[["sample", "dataset"]].copy()
    n_sh = len(sh_df)
    print(f"[sh] {n_sh} SH donors with full latents")
    assert n_sh == 101, f"Expected 101 SH donors, got {n_sh}"
    X_sh = sh_df[LATENT_COLS].to_numpy()

    # =====================================================================
    # Methods 1 + 2: scVI ordinal augmented + vanilla
    # =====================================================================
    boot_preds = {}
    orig_preds = {}

    n_jobs = int(os.environ.get("SLURM_CPUS_PER_TASK", "32"))
    print(f"[joblib] n_jobs = {n_jobs}")

    for augmented, name in [(True, M_AUG_ORD), (False, M_VAN_ORD)]:
        train_set = assemble_training(df, augmented=augmented)
        X_tr = train_set[LATENT_COLS].to_numpy()
        y_tr = train_set["y_train"].to_numpy().astype(int)
        print(f"[{name}] training n={len(y_tr)}, class counts="
              f"{np.bincount(y_tr, minlength=5).tolist()}")

        # Original (non-bootstrap) fit -> reference predictions
        clf = make_ordinal()
        clf.fit(X_tr, y_tr)
        orig = predict_argmax_lr(clf, X_sh)
        orig_preds[name] = orig
        print(f"[{name}] original SH counts F0..F4 = "
              f"{np.bincount(orig, minlength=5).tolist()}")

        # B=500 bootstrap fits in parallel
        results = Parallel(n_jobs=n_jobs, verbose=5)(
            delayed(boot_one_scvi)(b, X_tr, y_tr, X_sh, RNG_SEED)
            for b in range(B)
        )
        boot_preds[name] = np.vstack(results)
        print(f"[{name}] bootstrap matrix shape: {boot_preds[name].shape}")

    # =====================================================================
    # Methods 3 + 4: scANVI vanilla + augmented
    # =====================================================================
    for path, name, prefix, argmax_col in [
        (SCANVI_VAN, M_SCAN_V, "posterior_F", "F_stage_scanvi_argmax"),
        (SCANVI_AUG, M_SCAN_A, "posterior_F", "F_stage_scanvi_augmented_argmax"),
    ]:
        scan = pd.read_csv(path, sep="\t")
        # align to SH order
        scan_sh = sh_keys.merge(scan, on=["sample", "dataset"], how="left")
        post_cols = [f"{prefix}{k}_mean" for k in range(5)]
        post_mat = scan_sh[post_cols].to_numpy(dtype=float)
        # original argmax of posterior
        orig = scan_sh[argmax_col].to_numpy()
        # rows with missing posterior -> fallback to 0
        orig = np.where(np.isnan(orig), 0, orig).astype(int)
        orig_preds[name] = orig
        print(f"[{name}] original SH counts F0..F4 = "
              f"{np.bincount(orig, minlength=5).tolist()}")

        results = Parallel(n_jobs=n_jobs, verbose=5)(
            delayed(boot_one_scanvi)(b, post_mat, RNG_SEED)
            for b in range(B)
        )
        boot_preds[name] = np.vstack(results)
        print(f"[{name}] bootstrap matrix shape: {boot_preds[name].shape}")

    # =====================================================================
    # Stability metrics
    # =====================================================================
    stability = {}
    sh_dist_by_method = {}
    decision = {}
    for m in METHODS:
        stab = per_donor_stability(boot_preds[m], orig_preds[m])
        stab.insert(0, "method", m)
        stab.insert(0, "dataset", sh_keys["dataset"].values)
        stab.insert(0, "sample",  sh_keys["sample"].values)
        stab["original_pred"] = orig_preds[m]
        stability[m] = stab

        sh_dist = sh_distribution_per_boot(boot_preds[m])  # (5, B)
        sh_dist_by_method[m] = sh_dist

        # Headline numbers
        median_within_1 = float(np.median(stab["within_1_pct"]))
        pct_donors_ge_80 = float(np.mean(stab["within_1_pct"] >= 80.0) * 100)
        # F3 mode preservation across bootstraps
        argmax_per_boot = sh_dist.argmax(axis=0)  # 1D length B
        pct_boots_F3_mode = float(np.mean(argmax_per_boot == 3) * 100)

        # 95% CIs per stage
        ci = {}
        for k in range(5):
            ci[f"F{k}_lo"]  = float(np.percentile(sh_dist[k], 2.5))
            ci[f"F{k}_med"] = float(np.percentile(sh_dist[k], 50.0))
            ci[f"F{k}_hi"]  = float(np.percentile(sh_dist[k], 97.5))

        decision[m] = dict(
            median_within_1=median_within_1,
            pct_donors_ge_80=pct_donors_ge_80,
            pct_boots_F3_mode=pct_boots_F3_mode,
            **ci,
        )
        print(f"[stability] {m}: median within+/-1 = {median_within_1:.1f}%, "
              f"% donors >=80% stable = {pct_donors_ge_80:.1f}%, "
              f"% boots F3-mode = {pct_boots_F3_mode:.1f}%")

    # =====================================================================
    # Decision logic
    # =====================================================================
    aug = decision[M_AUG_ORD]
    pass_a = aug["pct_donors_ge_80"] >= 80.0
    pass_b = aug["pct_boots_F3_mode"] >= 95.0
    verdict = "SHIP" if (pass_a and pass_b) else "FALL_BACK"
    decision["verdict"] = verdict
    decision["pass_a_pct_donors_ge_80"] = pass_a
    decision["pass_b_pct_boots_F3_mode"] = pass_b

    # =====================================================================
    # Write outputs
    # =====================================================================
    per_donor_full = pd.concat([stability[m] for m in METHODS], ignore_index=True)
    per_donor_full.to_csv(OUT_PER_DONOR, sep="\t", index=False)
    print(f"[out] wrote {OUT_PER_DONOR}  ({len(per_donor_full)} rows)")

    # SH distribution stability: per method, per stage, per bootstrap
    rows = []
    for m in METHODS:
        d = sh_dist_by_method[m]
        for b in range(B):
            for k in range(5):
                rows.append(dict(method=m, boot=b, fstage=k, fraction=d[k, b]))
    pd.DataFrame(rows).to_csv(OUT_SH_DIST, sep="\t", index=False)
    print(f"[out] wrote {OUT_SH_DIST}")

    # Decision text
    with open(OUT_DECISION, "w") as fh:
        fh.write(f"Bootstrap-stability gate verdict\n")
        fh.write(f"=================================\n")
        fh.write(f"B = {B}, seed = {RNG_SEED}\n")
        fh.write(f"SH donors = {n_sh}\n\n")
        fh.write(f"Pre-registered thresholds:\n")
        fh.write(f"  (a) >=80% of SH donors with within_1_pct >= 80\n")
        fh.write(f"  (b) >=95% of bootstraps with SH F3 mode preserved\n\n")
        fh.write(f"Augmented scVI ordinal (gating method):\n")
        fh.write(f"  median within +/-1 = {aug['median_within_1']:.2f}%\n")
        fh.write(f"  %% donors >=80%% stable = {aug['pct_donors_ge_80']:.2f}%\n")
        fh.write(f"  %% boots F3-mode preserved = {aug['pct_boots_F3_mode']:.2f}%\n")
        fh.write(f"  pass (a) >=80%% donors: {pass_a}\n")
        fh.write(f"  pass (b) >=95%% boots F3-mode: {pass_b}\n\n")
        fh.write(f"VERDICT: {verdict}\n\n")
        fh.write("All methods (headline):\n")
        for m in METHODS:
            d = decision[m]
            fh.write(f"  {m}:\n")
            fh.write(f"    median within +/-1 = {d['median_within_1']:.2f}%\n")
            fh.write(f"    %% donors >=80%% stable = {d['pct_donors_ge_80']:.2f}%\n")
            fh.write(f"    %% boots F3-mode preserved = {d['pct_boots_F3_mode']:.2f}%\n")
            for k in range(5):
                fh.write(f"    F{k} median frac = {d[f'F{k}_med']:.3f} "
                         f"(95%% CI {d[f'F{k}_lo']:.3f} - {d[f'F{k}_hi']:.3f})\n")
        fh.write("\nIf VERDICT = FALL_BACK, BOOTSTRAP_FALLBACK_REQUIRED.flag "
                "is written and downstream should use scVI ordinal vanilla "
                "(rho=0.391).\n")
    print(f"[out] wrote {OUT_DECISION}")

    # Flag
    if verdict == "FALL_BACK":
        with open(OUT_FLAG, "w") as fh:
            fh.write(f"FALL_BACK required (verdict={verdict}). "
                     f"Use scVI ordinal vanilla (343b).\n")
            fh.write(f"AUG pct_donors_ge_80 = {aug['pct_donors_ge_80']:.2f}% "
                     f"(needed >=80%)\n")
            fh.write(f"AUG pct_boots_F3_mode = {aug['pct_boots_F3_mode']:.2f}% "
                     f"(needed >=95%)\n")
        print(f"[out] wrote FALLBACK FLAG: {OUT_FLAG}")
    else:
        # Remove stale flag if present
        if OUT_FLAG.exists():
            OUT_FLAG.unlink()

    # Figure
    make_figure(stability, sh_dist_by_method, decision)

    print(f"[done] verdict={verdict}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
