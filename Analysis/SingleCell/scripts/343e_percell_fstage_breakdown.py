#!/usr/bin/env python
"""
343e_percell_fstage_breakdown.py

Per-cell-resolution F-stage classifier. Train the SAME LogisticRegression
(multinomial, lbfgs) used in 343b on donor-level scVI pseudobulks for the 58
Andrews donors, but APPLY it to each individual hepatocyte's 20-d scVI
latent. Then aggregate the per-cell F-stage predictions per donor and per
disease_stage_coarse stratum to evaluate within-donor heterogeneity and
cross-dataset calibration.

This is the per-cell generalization of the donor-level argmax done by
Script 343b. We do NOT modify 343 / 343b / 343c / 343d.

Inputs
------
- Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad
    obsm['X_scVI'] 657,804 x 20; obs has sample, dataset
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv
    (has disease_stage_coarse)

Outputs
-------
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/percell_fstage_predictions.tsv.gz
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_cellfrac.tsv
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/cellfrac_by_stage.tsv
- Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/cellfrac_by_documented_fstage.tsv
- figures/supplementary/stage_ccc/figS_fstage_percell_breakdown.pdf
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import anndata as ad

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


PROJECT_ROOT = Path(
    os.environ.get(
        "MASLD_PROJECT_ROOT",
        "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
    )
)

HEP_H5AD = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad"
FSTAGE_DOC = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv"
DONOR_META = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv"
OUT_DIR = PROJECT_ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory"
FIG_DIR = PROJECT_ROOT / "figures/supplementary/stage_ccc"

PERCELL_TSV = OUT_DIR / "percell_fstage_predictions.tsv.gz"
DONOR_FRAC_TSV = OUT_DIR / "donor_fstage_cellfrac.tsv"
STAGE_TSV = OUT_DIR / "cellfrac_by_stage.tsv"
DOC_TSV = OUT_DIR / "cellfrac_by_documented_fstage.tsv"
FIG_PDF = FIG_DIR / "figS_fstage_percell_breakdown.pdf"

# MASLD palette (mirrors scripts/figures/publication_theme.R fibrosis_stage_colors)
FSTAGE_COLORS = {
    0: "#E3F2FD",  # F0 — palest blue
    1: "#90CAF9",  # F1
    2: "#42A5F5",  # F2
    3: "#1565C0",  # F3
    4: "#0D47A1",  # F4 — deepest blue
}
STAGE_ORDER = ["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"]
# Control = gray (CLAUDE.md control_gray rule); disease = magenta system
STAGE_COLORS = {
    "Healthy": "#9E9E9E",
    "Steatosis": "#F48FB1",
    "Steatohepatitis": "#C2185B",
    "Cirrhosis": "#880E4F",
}


def load_scvi_and_meta() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return (Z, samples, datasets) for 657K hepatocytes."""
    print(f"[load] {HEP_H5AD.name} (backed)")
    a = ad.read_h5ad(HEP_H5AD, backed="r")
    if "X_scVI" not in a.obsm:
        sys.exit(f"FATAL: X_scVI not in obsm. keys={list(a.obsm.keys())}")
    Z = np.asarray(a.obsm["X_scVI"], dtype=np.float32)
    samples = a.obs["sample"].astype(str).to_numpy()
    datasets = a.obs["dataset"].astype(str).to_numpy()
    print(f"[load] n_cells={Z.shape[0]} latent_dim={Z.shape[1]}")
    a.file.close()
    return Z, samples, datasets


def donor_pseudobulk(Z: np.ndarray, samples: np.ndarray, datasets: np.ndarray) -> pd.DataFrame:
    """Per-(sample,dataset) mean of the latent + cell count."""
    cols = [f"z{i}" for i in range(Z.shape[1])]
    df = pd.DataFrame(Z, columns=cols)
    df["sample"] = samples
    df["dataset"] = datasets
    pb = df.groupby(["sample", "dataset"], observed=True)[cols].mean().reset_index()
    n = df.groupby(["sample", "dataset"], observed=True).size().rename("n_hepatocytes").reset_index()
    return pb.merge(n, on=["sample", "dataset"]), cols


def fit_classifier(pb: pd.DataFrame, latent_cols: list[str]):
    """Train LogisticRegression on the 58 Andrews donors (Script 343b L127 spec).

    Returns (scaler, classifier, train_classes, X_train_scaled, y_train, train_idx).
    """
    fs = pd.read_csv(FSTAGE_DOC, sep="\t")
    fs["F_stage_documented"] = pd.to_numeric(fs["F_stage_documented"], errors="coerce")
    pb_lab = pb.merge(fs[["sample", "dataset", "F_stage_documented"]], on=["sample", "dataset"], how="left")
    train = pb_lab["F_stage_documented"].notna()
    print(f"[train] {int(train.sum())} Andrews-labeled donors")
    X_train_raw = pb_lab.loc[train, latent_cols].to_numpy(dtype=np.float64)
    y_train = pb_lab.loc[train, "F_stage_documented"].astype(int).to_numpy()
    print(f"[train] class counts: {dict(zip(*np.unique(y_train, return_counts=True)))}")

    scaler = StandardScaler().fit(X_train_raw)
    X_train = scaler.transform(X_train_raw)
    clf = LogisticRegression(
        solver="lbfgs", C=1.0, max_iter=5000, class_weight="balanced"
    )
    clf.fit(X_train, y_train)
    print(f"[train] classes_ = {clf.classes_.tolist()}")
    return scaler, clf, pb_lab, train


def pad_proba(proba: np.ndarray, classes: np.ndarray) -> np.ndarray:
    """Pad LR proba to full F0..F4 columns; renormalize."""
    full = np.zeros((proba.shape[0], 5), dtype=np.float32)
    for j, c in enumerate(classes):
        full[:, int(c)] = proba[:, j]
    rs = full.sum(axis=1, keepdims=True)
    rs[rs == 0] = 1.0
    return full / rs


def predict_percell(
    Z: np.ndarray, scaler: StandardScaler, clf: LogisticRegression
) -> tuple[np.ndarray, np.ndarray]:
    """Apply scaler+classifier to every cell. Returns (argmax, proba_5col)."""
    n = Z.shape[0]
    chunk = 200_000
    proba_full = np.empty((n, 5), dtype=np.float32)
    classes = clf.classes_
    for i in range(0, n, chunk):
        j = min(i + chunk, n)
        Xc = scaler.transform(Z[i:j].astype(np.float64))
        p = clf.predict_proba(Xc)
        proba_full[i:j] = pad_proba(p, classes)
        print(f"[predict] cells {i}-{j}/{n}")
    argmax = proba_full.argmax(axis=1).astype(np.int8)
    return argmax, proba_full


def entropy_of(vec: np.ndarray) -> float:
    """Shannon entropy in bits of a length-K probability vector."""
    v = np.clip(vec, 1e-12, 1.0)
    v = v / v.sum()
    return float(-np.sum(v * np.log2(v)))


def aggregate(
    samples: np.ndarray,
    datasets: np.ndarray,
    argmax: np.ndarray,
    donor_meta: pd.DataFrame,
    pb_lab: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Build donor cell-frac, by-stage, by-documented-F-stage tables."""
    df = pd.DataFrame({"sample": samples, "dataset": datasets, "F_argmax": argmax})
    # per-donor counts per F-stage
    counts = (
        df.groupby(["sample", "dataset", "F_argmax"], observed=True)
        .size()
        .unstack(fill_value=0)
    )
    # ensure all 5 columns
    for k in range(5):
        if k not in counts.columns:
            counts[k] = 0
    counts = counts[[0, 1, 2, 3, 4]]
    counts.columns = [f"F{k}_count" for k in range(5)]
    counts["n_cells"] = counts.sum(axis=1)
    fracs = counts[[f"F{k}_count" for k in range(5)]].div(counts["n_cells"], axis=0)
    fracs.columns = [f"F{k}_frac" for k in range(5)]
    donor = pd.concat([counts, fracs], axis=1).reset_index()
    # dominant + entropy
    donor["dominant_fstage"] = fracs.values.argmax(axis=1)
    donor["entropy"] = [entropy_of(v) for v in fracs.values]

    meta_cols = ["sample", "dataset", "disease_stage_coarse"]
    donor = donor.merge(donor_meta[meta_cols], on=["sample", "dataset"], how="left")
    # also attach documented F-stage
    donor = donor.merge(
        pb_lab[["sample", "dataset", "F_stage_documented"]], on=["sample", "dataset"], how="left"
    )

    # Aggregate by disease_stage_coarse: average per-donor F-stage cell-fractions
    rows = []
    for stage in STAGE_ORDER:
        sub = donor[donor["disease_stage_coarse"] == stage]
        if len(sub) == 0:
            continue
        means = sub[[f"F{k}_frac" for k in range(5)]].mean().values
        medians = sub[[f"F{k}_frac" for k in range(5)]].median().values
        q1 = sub[[f"F{k}_frac" for k in range(5)]].quantile(0.25).values
        q3 = sub[[f"F{k}_frac" for k in range(5)]].quantile(0.75).values
        rows.append({
            "disease_stage_coarse": stage,
            "n_donors": len(sub),
            "mean_entropy": float(sub["entropy"].mean()),
            "median_entropy": float(sub["entropy"].median()),
            **{f"F{k}_mean_frac": float(means[k]) for k in range(5)},
            **{f"F{k}_median_frac": float(medians[k]) for k in range(5)},
            **{f"F{k}_Q1": float(q1[k]) for k in range(5)},
            **{f"F{k}_Q3": float(q3[k]) for k in range(5)},
        })
    by_stage = pd.DataFrame(rows)

    # Aggregate by documented F-stage (Andrews only)
    rows_doc = []
    andrews = donor[donor["F_stage_documented"].notna()].copy()
    andrews["F_stage_documented"] = andrews["F_stage_documented"].astype(int)
    for k in range(5):
        sub = andrews[andrews["F_stage_documented"] == k]
        if len(sub) == 0:
            continue
        means = sub[[f"F{j}_frac" for j in range(5)]].mean().values
        rows_doc.append({
            "documented_fstage": k,
            "n_donors": len(sub),
            **{f"pred_F{j}_mean_frac": float(means[j]) for j in range(5)},
        })
    by_doc = pd.DataFrame(rows_doc)

    return donor, by_stage, by_doc


def make_figure(donor: pd.DataFrame, by_stage: pd.DataFrame, by_doc: pd.DataFrame) -> None:
    FIG_DIR.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({
        "pdf.fonttype": 42, "ps.fonttype": 42,
        "font.family": "Helvetica", "font.size": 8,
        "axes.labelsize": 9, "axes.titlesize": 10,
        "axes.spines.right": False, "axes.spines.top": False,
    })

    with PdfPages(FIG_PDF) as pdf:
        # Page 1: Panel A (stacked bars per stratum) + Panel B (Andrews heatmap) + Panel C (entropy box)
        fig = plt.figure(figsize=(13, 8.5))
        gs = fig.add_gridspec(2, 3, height_ratios=[1, 1], hspace=0.45, wspace=0.45)

        # Panel A: one stacked bar per stratum
        axA = fig.add_subplot(gs[0, 0:2])
        x_pos = np.arange(len(STAGE_ORDER))
        bottom = np.zeros(len(STAGE_ORDER))
        for k in range(5):
            heights = []
            for s in STAGE_ORDER:
                row = by_stage[by_stage["disease_stage_coarse"] == s]
                heights.append(float(row[f"F{k}_mean_frac"].iloc[0]) if len(row) else 0.0)
            heights = np.asarray(heights)
            axA.bar(
                x_pos, heights, bottom=bottom,
                color=FSTAGE_COLORS[k], edgecolor="black", linewidth=0.4,
                label=f"F{k}",
            )
            bottom += heights
        axA.set_xticks(x_pos)
        labels = []
        for s in STAGE_ORDER:
            row = by_stage[by_stage["disease_stage_coarse"] == s]
            n = int(row["n_donors"].iloc[0]) if len(row) else 0
            labels.append(f"{s}\n(n={n})")
        axA.set_xticklabels(labels)
        axA.set_ylim(0, 1.0)
        axA.set_ylabel("Mean fraction of hepatocytes\nassigned to each F-stage")
        axA.set_title("A. Per-cell F-stage breakdown by disease_stage_coarse",
                      loc="left", fontweight="bold")
        axA.legend(title="Predicted F-stage", bbox_to_anchor=(1.02, 1.0), loc="upper left", frameon=False)

        # Panel B: Andrews diagonal-check heatmap
        axB = fig.add_subplot(gs[0, 2])
        if len(by_doc) > 0:
            M = np.zeros((5, 5))
            row_idx = {int(r): i for i, r in enumerate(by_doc["documented_fstage"].values)}
            for _, r in by_doc.iterrows():
                i = int(r["documented_fstage"])
                for j in range(5):
                    M[i, j] = r[f"pred_F{j}_mean_frac"]
            im = axB.imshow(M, cmap="magma", vmin=0, vmax=1, aspect="auto")
            for i in range(5):
                for j in range(5):
                    txt = f"{M[i, j]:.2f}"
                    color = "white" if M[i, j] < 0.5 else "black"
                    axB.text(j, i, txt, ha="center", va="center", color=color, fontsize=7)
            axB.set_xticks(range(5)); axB.set_yticks(range(5))
            axB.set_xticklabels([f"pred F{k}" for k in range(5)], rotation=45, ha="right")
            n_per = by_doc.set_index("documented_fstage")["n_donors"].to_dict()
            axB.set_yticklabels([f"doc F{k}\n(n={int(n_per.get(k, 0))})" for k in range(5)])
            axB.set_title("B. Andrews diagonal-check\n(rows=documented, cols=predicted)",
                          loc="left", fontweight="bold", fontsize=9)
            plt.colorbar(im, ax=axB, fraction=0.045, pad=0.04, label="mean cell-fraction")

        # Panel C: entropy boxplot per stratum
        axC = fig.add_subplot(gs[1, 0:2])
        data = []
        positions = []
        colors = []
        for i, s in enumerate(STAGE_ORDER):
            sub = donor[donor["disease_stage_coarse"] == s]["entropy"].dropna().values
            if len(sub) == 0:
                continue
            data.append(sub); positions.append(i); colors.append(STAGE_COLORS[s])
        bp = axC.boxplot(data, positions=positions, widths=0.55, patch_artist=True,
                         showfliers=False, medianprops=dict(color="black", lw=1.2))
        for patch, c in zip(bp["boxes"], colors):
            patch.set_facecolor(c); patch.set_alpha(0.7); patch.set_edgecolor("black")
        for i, d in zip(positions, data):
            jitter = np.random.RandomState(42 + i).uniform(-0.12, 0.12, size=len(d))
            axC.scatter(np.full(len(d), i) + jitter, d, s=8, color="black", alpha=0.45, lw=0)
        axC.set_xticks(range(len(STAGE_ORDER)))
        axC.set_xticklabels(STAGE_ORDER)
        axC.set_ylabel("Per-donor entropy of\nper-cell F-stage distribution (bits)")
        axC.set_title("C. Within-donor heterogeneity (entropy across F0-F4)",
                      loc="left", fontweight="bold")
        axC.axhline(np.log2(5), color="grey", ls=":", lw=0.8)
        axC.text(len(STAGE_ORDER) - 0.5, np.log2(5) - 0.05, "uniform max", color="grey",
                 fontsize=7, ha="right", va="top")

        # Panel D blank used for legend already, skip
        axD = fig.add_subplot(gs[1, 2])
        axD.axis("off")
        axD.text(0.0, 0.95,
                 "Training: LogisticRegression(multinomial,\nlbfgs, class_weight=balanced)\non 58 Andrews donor scVI pseudobulks.",
                 fontsize=7, va="top")
        axD.text(0.0, 0.55,
                 "Inference: per-cell scVI z-scored\nwith same scaler, predict_proba\non each of the 657K hepatocytes.",
                 fontsize=7, va="top")
        axD.text(0.0, 0.18,
                 "Aggregation: per-donor cell-fraction\nargmax + mean per disease_stage_coarse.",
                 fontsize=7, va="top")

        pdf.savefig(fig, bbox_inches="tight"); plt.close(fig)

    print(f"[fig] {FIG_PDF}")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    Z, samples, datasets = load_scvi_and_meta()
    pb, latent_cols = donor_pseudobulk(Z, samples, datasets)
    print(f"[pseudobulk] {len(pb)} donor x dataset rows")

    scaler, clf, pb_lab, train_mask = fit_classifier(pb, latent_cols)

    print("[predict] applying classifier to all 657K cells in chunks...")
    argmax, proba = predict_percell(Z, scaler, clf)

    # Per-cell table
    percell = pd.DataFrame({
        "sample": samples,
        "dataset": datasets,
        "F_stage_argmax": argmax.astype(int),
        "P_F0": proba[:, 0].round(5),
        "P_F1": proba[:, 1].round(5),
        "P_F2": proba[:, 2].round(5),
        "P_F3": proba[:, 3].round(5),
        "P_F4": proba[:, 4].round(5),
    })
    percell.to_csv(PERCELL_TSV, sep="\t", index=False, compression="gzip")
    print(f"[output] {PERCELL_TSV} ({len(percell)} rows)")

    donor_meta = pd.read_csv(DONOR_META, sep="\t")
    keep_cols = ["sample", "dataset", "disease_stage_coarse"]
    donor_meta = donor_meta[keep_cols].drop_duplicates()

    donor_frac, by_stage, by_doc = aggregate(
        samples, datasets, argmax, donor_meta, pb_lab
    )
    donor_frac.to_csv(DONOR_FRAC_TSV, sep="\t", index=False)
    by_stage.to_csv(STAGE_TSV, sep="\t", index=False)
    by_doc.to_csv(DOC_TSV, sep="\t", index=False)
    print(f"[output] {DONOR_FRAC_TSV} ({len(donor_frac)} donors)")
    print(f"[output] {STAGE_TSV}")
    print(f"[output] {DOC_TSV}")

    # Mode-of-predictions vs donor-pseudobulk argmax comparison
    pb_pred_X = scaler.transform(pb[latent_cols].to_numpy(dtype=np.float64))
    pb_pred = pad_proba(clf.predict_proba(pb_pred_X), clf.classes_)
    pb_argmax = pb_pred.argmax(axis=1)
    pb_donor = pb[["sample", "dataset"]].copy()
    pb_donor["donor_pseudobulk_argmax"] = pb_argmax
    cmp = donor_frac.merge(pb_donor, on=["sample", "dataset"], how="left")
    n_agree = (cmp["dominant_fstage"] == cmp["donor_pseudobulk_argmax"]).sum()
    print(f"\n[compare] donor pseudobulk argmax vs cell-mode argmax agreement: "
          f"{n_agree}/{len(cmp)} = {100*n_agree/len(cmp):.1f}%")
    cmp[["sample", "dataset", "disease_stage_coarse", "dominant_fstage",
         "donor_pseudobulk_argmax", "entropy", "n_cells"]].to_csv(
        OUT_DIR / "donor_argmax_vs_cellmode.tsv", sep="\t", index=False
    )

    # Quick console summary
    print("\n[summary] mean cell-fraction per stratum (4x5):")
    print(by_stage[["disease_stage_coarse", "n_donors",
                    "F0_mean_frac", "F1_mean_frac", "F2_mean_frac",
                    "F3_mean_frac", "F4_mean_frac", "mean_entropy"]].to_string(index=False))
    print("\n[summary] Andrews diagonal-check (5x5):")
    print(by_doc.to_string(index=False))

    make_figure(donor_frac, by_stage, by_doc)


if __name__ == "__main__":
    main()
