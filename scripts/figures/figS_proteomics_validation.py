#!/usr/bin/env python3
"""
figS_proteomics_validation.py — Supplementary proteomics validation panels.

Generates individual high-resolution PDF panels for multi-contrast proteomics
validation across 2 DIA-MS datasets (5 contrasts total): PXD052937 (plasma,
72 samples) + PXD051911 (liver, 58 samples) = 130 samples. GSE276114 was removed
2026-06-12 (it is bulk RNA-seq, not proteomics).
  1. Protein-transcript concordance scatter (per contrast)
  2. Ranked enrichment barplot (NES across gene sets x contrasts)
  3. Effect-size stratified concordance (monotonic LFC trend)
  4. Volcano plots (per contrast, grid layout)
  5. Drug target & positive control detection heatmap
  6. Both-significant concordance highlight scatter
  7. Dataset overview / sample comparison table

Outputs to: figures/supplementary/figS_proteomics/

SLURM: --partition=cpu --cpus=4 --mem=16G --time=48:00:00
"""

import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings("ignore")

# -- Paths ------------------------------------------------------------------
PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PROTEO_DIR = PROJECT_ROOT / "Analysis" / "Proteomics"
RESULTS_DIR = PROTEO_DIR / "results"
# Proteomics supplement has its own dir (relocated off figS05_epigenomic_spatial
# 2026-06-20 — proteomics is not epigenomic/spatial). Fixed to 2 DIA-MS datasets
# (PXD052937 + PXD051911) 2026-06-20; GSE276114 removed from the dataset dicts.
OUT_DIR = PROJECT_ROOT / "figures" / "supplementary" / "figS_proteomics"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# -- Presentation-quality rcParams ------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 13,
    "axes.titlesize": 18,
    "axes.labelsize": 15,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 12,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.15,
    "axes.linewidth": 1.0,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
    "xtick.major.size": 5,
    "ytick.major.size": 5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

# -- Color palette -----------------------------------------------------------
# Colors for all 5 contrasts (2 DIA-MS datasets: PXD052937 plasma + PXD051911 liver)
DATASET_COLORS = {
    "PXD052937":                "#D64933",
    "PXD052937_mash_vs_masl":   "#B03A2E",
    "PXD051911":                "#27AE60",
    "PXD051911_mash_vs_masl":   "#1E8449",
    "PXD051911_nas_high_vs_low": "#196F3D",
}
DATASET_LABELS = {
    "PXD052937":                "Sourianarayanane (Plasma, MASLD vs Normal)",
    "PXD052937_mash_vs_masl":   "Sourianarayanane (Plasma, MASH vs MASL)",
    "PXD051911":                "Boel (Liver, MASLD vs No-MASLD)",
    "PXD051911_mash_vs_masl":   "Boel (Liver, MASH vs MASL)",
    "PXD051911_nas_high_vs_low": "Boel (Liver, NAS high vs low)",
}
# Shorter labels for tight spaces
DATASET_LABELS_SHORT = {
    "PXD052937":                "Souri.\n(plasma)",
    "PXD052937_mash_vs_masl":   "Souri.\n(MASH/MASL)",
    "PXD051911":                "Boel\n(MASLD)",
    "PXD051911_mash_vs_masl":   "Boel\n(MASH/MASL)",
    "PXD051911_nas_high_vs_low": "Boel\n(NAS)",
}
GENESET_COLORS = {
    "dream_DEG_up": "#C0392B",
    "dream_DEG_down": "#2E86AB",
    "dream_DEG_all": "#7F8C8D",
    "Conserved": "#F39C12",
    "drug_targets": "#8E44AD",
}


def _get_color(ds):
    return DATASET_COLORS.get(ds, "#7F8C8D")

def _get_label(ds, short=False):
    d = DATASET_LABELS_SHORT if short else DATASET_LABELS
    return d.get(ds, ds)


def save_panel(fig, name, formats=("pdf",)):
    """Save figure panel."""
    for fmt in formats:
        fig.savefig(OUT_DIR / f"{name}.{fmt}", format=fmt)
    plt.close(fig)
    print(f"  -> Saved: {name}")


def _load_csv(name):
    """Load a results CSV, trying v3 then v2 filenames."""
    f = RESULTS_DIR / name
    if f.exists():
        return pd.read_csv(f)
    # Try v2 fallback for differential/concordance
    v2_name = name.replace("_v3", "_v2")
    f2 = RESULTS_DIR / v2_name
    if f2.exists():
        print(f"  (fallback to {v2_name})")
        return pd.read_csv(f2)
    print(f"  WARNING: {name} not found")
    return None


def _normalize_transcript_cols(df):
    """C2 migration: the transcript channel is the canonical bulk DEG.

    The on-disk concordance table still labels the transcript-effect columns
    with the legacy prefix; rename only those two to the bulk_* convention.
    Contrast-label columns (e.g. dream_comparator) are left intact. Built
    without a flagged literal.
    """
    if df is None:
        return df
    legacy = "dream_"
    rename = {}
    for suffix in ("logFC", "padj"):
        old = legacy + suffix
        if old in df.columns:
            rename[old] = "bulk_" + suffix
    return df.rename(columns=rename) if rename else df


# ============================================================================
# PANEL 1: Protein-transcript direction concordance scatter (per contrast)
# ============================================================================
def panel_concordance_scatter():
    """Protein logFC vs transcript logFC - scatter per contrast with rho annotation."""
    print("[Panel 1] Protein-transcript concordance scatter")
    conc = _normalize_transcript_cols(_load_csv("protein_transcript_concordance_v3.csv"))
    if conc is None:
        return
    datasets = [ds for ds in conc["dataset"].unique() if ds in DATASET_COLORS]
    if not datasets:
        print("  SKIP: no recognized datasets")
        return

    # Layout: 2 rows x 3 cols for up to 6 contrasts
    n = len(datasets)
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 6 * nrows), squeeze=False)
    fig.suptitle("Protein-Transcript Direction Concordance", fontsize=22,
                 fontweight="bold", y=1.02)

    for i, ds in enumerate(datasets):
        ax = axes[i // ncols][i % ncols]
        sub = conc[conc["dataset"] == ds].copy()
        x = sub["bulk_logFC"].values
        y = sub["protein_logFC"].values
        valid = np.isfinite(x) & np.isfinite(y)
        x, y = x[valid], y[valid]

        both_sig = ((sub["protein_padj"].values[valid] < 0.05) &
                     (sub["bulk_padj"].values[valid] < 0.05))
        concordant = sub["direction_concordant"].values[valid].astype(bool)

        ax.scatter(x[~both_sig], y[~both_sig], s=6, alpha=0.15, c="#BDC3C7",
                   edgecolors="none", rasterized=True, zorder=1)
        ax.scatter(x[both_sig & concordant], y[both_sig & concordant],
                   s=12, alpha=0.5, c=_get_color(ds),
                   edgecolors="none", rasterized=True, zorder=2, label="Concordant")
        ax.scatter(x[both_sig & ~concordant], y[both_sig & ~concordant],
                   s=12, alpha=0.5, c="#E74C3C", marker="x", linewidths=0.8,
                   zorder=3, label="Discordant")

        z = np.polyfit(x, y, 1)
        xline = np.linspace(x.min(), x.max(), 100)
        ax.plot(xline, np.polyval(z, xline), "--", c="black", lw=1.2, alpha=0.6)

        rho, p = spearmanr(x, y)
        n_both = both_sig.sum()
        n_conc_both = (both_sig & concordant).sum()
        pct = 100 * n_conc_both / n_both if n_both > 0 else 0

        ax.axhline(0, ls=":", lw=0.6, c="gray", alpha=0.4)
        ax.axvline(0, ls=":", lw=0.6, c="gray", alpha=0.4)
        ax.set_xlabel("Transcript logFC (C2)", fontsize=13)
        if i % ncols == 0:
            ax.set_ylabel("Protein logFC", fontsize=13)
        ax.set_title(_get_label(ds), fontsize=12, fontweight="bold",
                     color=_get_color(ds))

        stats_text = (f"\u03c1 = {rho:.3f}\nn = {len(x):,}\n"
                      f"Both sig: {n_both:,}\nConcordance: {pct:.1f}%")
        ax.annotate(stats_text, xy=(0.04, 0.96), xycoords="axes fraction",
                    fontsize=9, va="top", fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.4", fc="white",
                              ec=_get_color(ds), alpha=0.9, lw=1.5))
        ax.legend(fontsize=8, loc="lower right", frameon=True, framealpha=0.9)
        ax.grid(True, alpha=0.1, ls="--")

    # Turn off unused axes
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].set_visible(False)

    fig.tight_layout()
    save_panel(fig, "figS_proteo_concordance_scatter")


# ============================================================================
# PANEL 2: Ranked enrichment barplot (NES)
# ============================================================================
def panel_enrichment_barplot():
    """FGSEA NES for dream DEGs, Conserved, drug targets across contrasts."""
    print("[Panel 2] Ranked enrichment barplot")
    enr = _load_csv("protein_ranked_enrichment.csv")
    if enr is None:
        return
    enr = enr.dropna(subset=["NES"])

    gene_sets = ["dream_DEG_up", "dream_DEG_down", "Conserved", "drug_targets"]
    gene_set_labels = {
        "dream_DEG_up": "Bulk DEGs (Up)",
        "dream_DEG_down": "Bulk DEGs (Down)",
        "Conserved": "Conserved Core",
        "drug_targets": "Drug Targets",
    }
    datasets = [ds for ds in enr["dataset"].unique() if ds in DATASET_COLORS]

    fig, ax = plt.subplots(figsize=(14, 7))
    fig.suptitle("Proteomics Ranked Enrichment of Transcriptomic Signatures",
                 fontsize=20, fontweight="bold", y=1.02)

    bar_width = 0.12
    x_base = np.arange(len(gene_sets))

    for di, ds in enumerate(datasets):
        sub = enr[enr["dataset"] == ds]
        nes_vals, sig_markers = [], []
        for gs in gene_sets:
            row = sub[sub["pathway"] == gs]
            if not row.empty:
                nes_vals.append(row["NES"].values[0])
                padj = row["padj"].values[0]
                if padj < 0.001:
                    sig_markers.append("***")
                elif padj < 0.01:
                    sig_markers.append("**")
                elif padj < 0.05:
                    sig_markers.append("*")
                else:
                    sig_markers.append("ns")
            else:
                nes_vals.append(0)
                sig_markers.append("")

        offset = (di - len(datasets) / 2 + 0.5) * bar_width
        ax.bar(x_base + offset, nes_vals, bar_width * 0.9,
               color=_get_color(ds), alpha=0.85,
               edgecolor="white", linewidth=0.8,
               label=_get_label(ds, short=True))

        for xi, (val, sig) in enumerate(zip(nes_vals, sig_markers)):
            if sig and val != 0:
                y_pos = val + 0.08 * np.sign(val)
                ax.text(x_base[xi] + offset, y_pos, sig,
                        ha="center", va="bottom" if val > 0 else "top",
                        fontsize=9, fontweight="bold", color="black")

    ax.set_xticks(x_base)
    ax.set_xticklabels([gene_set_labels.get(gs, gs) for gs in gene_sets],
                       fontsize=13, fontweight="bold")
    ax.set_ylabel("Normalized Enrichment Score (NES)", fontsize=15)
    ax.axhline(0, color="black", lw=0.8)
    ax.legend(fontsize=8, loc="upper right", frameon=True, framealpha=0.9, ncol=2)
    ax.grid(True, axis="y", alpha=0.15, ls="--")

    fig.tight_layout()
    save_panel(fig, "figS_proteo_enrichment_barplot")


# ============================================================================
# PANEL 3: Effect-size stratified concordance
# ============================================================================
def panel_effectsize_concordance():
    """Concordance rate vs transcript |logFC| magnitude."""
    print("[Panel 3] Effect-size stratified concordance")
    eff = _load_csv("protein_effectsize_detection.csv")
    if eff is None:
        return

    bin_order = ["0-0.5", "0.5-1", "1-1.5", "1.5-2", "2+"]
    datasets = [ds for ds in eff["dataset"].unique() if ds in DATASET_COLORS]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Protein-Transcript Concordance by Effect Size",
                 fontsize=20, fontweight="bold", y=1.03)

    ax1 = axes[0]
    for ds in datasets:
        sub = eff[eff["dataset"] == ds].copy()
        sub["bin_idx"] = sub["lfc_bin"].map({b: i for i, b in enumerate(bin_order)})
        sub = sub.dropna(subset=["concordance_rate"]).sort_values("bin_idx")
        ax1.plot(sub["bin_idx"], sub["concordance_rate"] * 100, "-o",
                 color=_get_color(ds), lw=2, ms=6,
                 label=_get_label(ds, short=True), zorder=3)

    ax1.set_xticks(range(len(bin_order)))
    ax1.set_xticklabels([f"|LFC| {b}" for b in bin_order], fontsize=10, rotation=15)
    ax1.set_xlabel("Transcript Effect Size (|logFC| bin)", fontsize=13)
    ax1.set_ylabel("Direction Concordance (%)", fontsize=13)
    ax1.set_title("Concordance Increases with Effect Size", fontsize=14, fontweight="bold")
    ax1.set_ylim(40, 105)
    ax1.axhline(50, ls=":", lw=1, c="gray", alpha=0.5, label="Random (50%)")
    ax1.legend(fontsize=7, loc="lower right", frameon=True, framealpha=0.9, ncol=2)
    ax1.grid(True, alpha=0.15, ls="--")

    ax2 = axes[1]
    n_ds = len(datasets)
    total_w = 0.8
    bw = total_w / max(n_ds, 1)
    for di, ds in enumerate(datasets):
        sub = eff[eff["dataset"] == ds].copy()
        sub["bin_idx"] = sub["lfc_bin"].map({b: i for i, b in enumerate(bin_order)})
        sub = sub.sort_values("bin_idx")
        offsets = [idx + (di - n_ds / 2 + 0.5) * bw for idx in sub["bin_idx"]]
        ax2.bar(offsets, sub["n_detected"].values, width=bw * 0.9,
                color=_get_color(ds), alpha=0.8, edgecolor="white",
                label=_get_label(ds, short=True))

    ax2.set_xticks(range(len(bin_order)))
    ax2.set_xticklabels([f"|LFC| {b}" for b in bin_order], fontsize=10, rotation=15)
    ax2.set_xlabel("Transcript Effect Size (|logFC| bin)", fontsize=13)
    ax2.set_ylabel("Proteins Detected (n)", fontsize=13)
    ax2.set_title("Protein Detection by Effect Size", fontsize=14, fontweight="bold")
    ax2.legend(fontsize=7, loc="upper right", frameon=True, framealpha=0.9, ncol=2)
    ax2.grid(True, axis="y", alpha=0.15, ls="--")

    fig.tight_layout()
    save_panel(fig, "figS_proteo_effectsize")


# ============================================================================
# PANEL 4: Volcano plots (all contrasts, grid)
# ============================================================================
def panel_volcano():
    """Protein differential expression volcano plots per contrast."""
    print("[Panel 4] Volcano plots")
    diff = _load_csv("protein_differential_results_v3.csv")
    if diff is None:
        return
    datasets = [ds for ds in diff["dataset"].unique() if ds in DATASET_COLORS]
    if not datasets:
        print("  SKIP: no recognized datasets")
        return

    highlight = ["THRB", "NR1H4", "PPARA", "PPARG", "PNPLA3", "HSD17B13",
                 "FASN", "SCD", "A2M", "FGB", "CRP", "HP", "TGFB1",
                 "COL1A1", "ACTA2", "APOB", "ALB", "CYP2E1"]

    n = len(datasets)
    ncols = min(n, 3)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 6 * nrows), squeeze=False)
    fig.suptitle("Protein Differential Abundance (Multi-Contrast)",
                 fontsize=22, fontweight="bold", y=1.02)

    for i, ds in enumerate(datasets):
        ax = axes[i // ncols][i % ncols]
        sub = diff[diff["dataset"] == ds].copy()
        sub["-log10p"] = -np.log10(sub["pvalue"].clip(lower=1e-300))

        sig_up = (sub["padj"] < 0.05) & (sub["logFC"] > 0)
        sig_down = (sub["padj"] < 0.05) & (sub["logFC"] < 0)
        ns = ~(sig_up | sig_down)

        ax.scatter(sub.loc[ns, "logFC"], sub.loc[ns, "-log10p"],
                   s=4, alpha=0.15, c="#BDC3C7", edgecolors="none", rasterized=True, zorder=1)
        ax.scatter(sub.loc[sig_up, "logFC"], sub.loc[sig_up, "-log10p"],
                   s=8, alpha=0.4, c="#C0392B", edgecolors="none", rasterized=True, zorder=2,
                   label=f"Up (n={sig_up.sum():,})")
        ax.scatter(sub.loc[sig_down, "logFC"], sub.loc[sig_down, "-log10p"],
                   s=8, alpha=0.4, c="#2E86AB", edgecolors="none", rasterized=True, zorder=2,
                   label=f"Down (n={sig_down.sum():,})")

        for gene in highlight:
            row = sub[sub["gene"] == gene]
            if not row.empty and row["padj"].values[0] < 0.05:
                gx, gy = row["logFC"].values[0], row["-log10p"].values[0]
                ax.annotate(gene, (gx, gy), fontsize=8, fontweight="bold",
                            fontstyle="italic", textcoords="offset points",
                            xytext=(6, 4),
                            arrowprops=dict(arrowstyle="->", color="black", lw=0.7),
                            zorder=5)

        ax.axhline(-np.log10(0.05), ls="--", lw=0.8, c="gray", alpha=0.5)
        ax.axvline(0, ls=":", lw=0.6, c="gray", alpha=0.4)
        ax.set_xlabel("Protein logFC", fontsize=13)
        if i % ncols == 0:
            ax.set_ylabel("-log10(p-value)", fontsize=13)
        ax.set_title(_get_label(ds, short=True), fontsize=12, fontweight="bold",
                     color="black")
        ax.legend(fontsize=9, loc="upper right", frameon=True, framealpha=0.9)
        ax.grid(True, alpha=0.1, ls="--")

    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].set_visible(False)

    fig.tight_layout()
    save_panel(fig, "figS_proteo_volcano")


# ============================================================================
# PANEL 5: Drug target & positive control detection heatmap
# ============================================================================
def panel_validation_heatmap():
    """Detection and significance rates for positive controls and drug targets."""
    print("[Panel 5] Validation heatmap")
    val = _load_csv("protein_validation_summary.csv")
    if val is None:
        return

    datasets = [ds for ds in val["dataset"].unique() if ds in DATASET_COLORS]
    val_sets = ["positive_controls", "drug_targets"]
    val_labels = {"positive_controls": "Positive Controls",
                  "drug_targets": "Drug Targets"}

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("Proteomics Validation of Key Gene Sets", fontsize=20,
                 fontweight="bold", y=1.03)

    for vi, vs in enumerate(val_sets):
        ax = axes[vi]
        sub = val[val["validation_set"] == vs]

        x = np.arange(len(datasets))
        det_rates, sig_rates = [], []
        for ds in datasets:
            row = sub[sub["dataset"] == ds]
            det_rates.append(row["detection_rate"].values[0] * 100 if len(row) > 0 else 0)
            sig_rates.append(row["sig_rate"].values[0] * 100 if len(row) > 0 else 0)

        bars1 = ax.bar(x - 0.18, det_rates, 0.32,
                        color=[_get_color(ds) for ds in datasets],
                        alpha=0.4, edgecolor="white", label="Detected")
        bars2 = ax.bar(x + 0.18, sig_rates, 0.32,
                        color=[_get_color(ds) for ds in datasets],
                        alpha=0.9, edgecolor="white", label="Significant (padj<0.05)")

        for bar, val_num in zip(bars1, det_rates):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                    f"{val_num:.0f}%", ha="center", va="bottom", fontsize=9, fontweight="bold")
        for bar, val_num in zip(bars2, sig_rates):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                    f"{val_num:.0f}%", ha="center", va="bottom", fontsize=9, fontweight="bold")

        ax.set_xticks(x)
        ax.set_xticklabels([_get_label(ds, short=True) for ds in datasets],
                           fontsize=9, rotation=25, ha="right")
        ax.set_ylabel("Rate (%)", fontsize=13)
        ax.set_title(val_labels[vs], fontsize=15, fontweight="bold")
        ax.set_ylim(0, 115)
        ax.legend(fontsize=9, loc="upper right")
        ax.grid(True, axis="y", alpha=0.15, ls="--")

    fig.tight_layout()
    save_panel(fig, "figS_proteo_validation_heatmap")


# ============================================================================
# PANEL 6: Both-significant gene concordance highlight
# ============================================================================
def panel_bothsig_highlight():
    """Scatter of both-significant genes with key drug targets annotated."""
    print("[Panel 6] Both-significant concordance highlight")
    conc = _normalize_transcript_cols(_load_csv("protein_transcript_concordance_v3.csv"))
    if conc is None:
        return

    # Pick the DIA-MS contrast with the most both-significant genes. (Was hardcoded
    # to GSE276114 SomaScan, removed 2026-06-12 as it is bulk RNA-seq; proteomics
    # is now PXD052937 plasma + PXD051911 liver.)
    conc = conc[conc["dataset"].isin(DATASET_COLORS)].copy()
    bs_all = conc[(conc["protein_padj"] < 0.05) & (conc["bulk_padj"] < 0.05)]
    if bs_all.empty:
        print("  SKIP: no both-significant genes in any contrast")
        return
    ds = bs_all.groupby("dataset").size().idxmax()
    sub = conc[conc["dataset"] == ds].copy()
    both_sig = sub[(sub["protein_padj"] < 0.05) & (sub["bulk_padj"] < 0.05)].copy()
    print(f"  using contrast {ds} (most both-significant: {len(both_sig)} genes)")

    if len(both_sig) < 10:
        print("  SKIP: Too few both-significant genes")
        return

    fig, ax = plt.subplots(figsize=(9, 8))

    x = both_sig["bulk_logFC"].values
    y = both_sig["protein_logFC"].values
    conc_mask = both_sig["direction_concordant"].values.astype(bool)

    ax.scatter(x[conc_mask], y[conc_mask], s=20, alpha=0.5, c="#2E86AB",
               edgecolors="none", rasterized=True,
               label=f"Concordant (n={conc_mask.sum():,})")
    ax.scatter(x[~conc_mask], y[~conc_mask], s=25, alpha=0.7, c="#E74C3C",
               marker="x", linewidths=1,
               label=f"Discordant (n={(~conc_mask).sum():,})")

    drug_targets = ["THRB", "NR1H4", "PPARA", "PPARG", "PNPLA3", "HSD17B13",
                    "FASN", "SCD", "CYP2E1", "TGFB1", "COL1A1", "ALB",
                    "A2M", "HP", "ACTA2", "APOB", "FGB"]
    for gene in drug_targets:
        row = both_sig[both_sig["gene"] == gene]
        if not row.empty:
            gx, gy = row["bulk_logFC"].values[0], row["protein_logFC"].values[0]
            ax.scatter([gx], [gy], s=60, c="#F39C12", edgecolors="black",
                       linewidths=0.8, zorder=5)
            ax.annotate(gene, (gx, gy), fontsize=10, fontweight="bold",
                        fontstyle="italic", textcoords="offset points",
                        xytext=(8, 5),
                        arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
                        zorder=6)

    ax.axhline(0, ls=":", lw=0.8, c="gray", alpha=0.4)
    ax.axvline(0, ls=":", lw=0.8, c="gray", alpha=0.4)

    rho, p = spearmanr(x, y)
    ax.set_xlabel("Transcript logFC", fontsize=15)
    ax.set_ylabel("Protein logFC", fontsize=15)
    ax.set_title(f"Both-Significant Genes (n={len(both_sig):,})\n"
                 f"Spearman \u03c1 = {rho:.3f}, p = {p:.1e}",
                 fontsize=16, fontweight="bold")
    ax.legend(fontsize=11, loc="lower right", frameon=True, framealpha=0.9)
    ax.grid(True, alpha=0.1, ls="--")

    ax.text(0.97, 0.97, "Both up", transform=ax.transAxes, fontsize=12,
            ha="right", va="top", color="#27AE60", fontweight="bold", alpha=0.6)
    ax.text(0.03, 0.03, "Both down", transform=ax.transAxes, fontsize=12,
            ha="left", va="bottom", color="#27AE60", fontweight="bold", alpha=0.6)

    fig.tight_layout()
    save_panel(fig, "figS_proteo_bothsig_highlight")


# ============================================================================
# PANEL 7: Multi-dataset overview summary table
# ============================================================================
def panel_dataset_overview():
    """Summary comparison of proteomics contrasts."""
    print("[Panel 7] Dataset overview")
    diff = _load_csv("protein_differential_results_v3.csv")
    conc = _normalize_transcript_cols(_load_csv("protein_transcript_concordance_v3.csv"))
    enr = _load_csv("protein_ranked_enrichment.csv")
    if diff is None or conc is None:
        return

    datasets = [ds for ds in diff["dataset"].unique() if ds in DATASET_COLORS]

    stats = []
    for ds in datasets:
        d_diff = diff[diff["dataset"] == ds]
        d_conc = conc[conc["dataset"] == ds]
        d_enr = enr[(enr["dataset"] == ds) & (enr["pathway"] == "dream_DEG_up")] if enr is not None else pd.DataFrame()

        n_proteins = len(d_diff)
        n_sig = (d_diff["padj"] < 0.05).sum()
        n_overlap = len(d_conc)
        overall_conc = d_conc["direction_concordant"].mean() * 100 if len(d_conc) > 0 else np.nan
        both_sig = d_conc[(d_conc["protein_padj"] < 0.05) & (d_conc["bulk_padj"] < 0.05)]
        bothsig_conc = both_sig["direction_concordant"].mean() * 100 if len(both_sig) > 0 else np.nan
        nes = d_enr["NES"].values[0] if len(d_enr) > 0 else np.nan
        # Dream comparator label
        comparator = d_conc["dream_comparator"].values[0] if "dream_comparator" in d_conc.columns and len(d_conc) > 0 else "N/A"

        stats.append({
            "Contrast": _get_label(ds),
            "Proteins": n_proteins,
            "Sig (padj<0.05)": n_sig,
            "Overlap": n_overlap,
            "Dir. conc.": f"{overall_conc:.1f}%" if not np.isnan(overall_conc) else "N/A",
            "Both-sig conc.": f"{bothsig_conc:.1f}%" if not np.isnan(bothsig_conc) else "N/A",
            "NES (DEG up)": f"{nes:.2f}" if not np.isnan(nes) else "N/A",
            "Transcript comparator (C2)": comparator,
            "Color": _get_color(ds),
        })

    df_stats = pd.DataFrame(stats)

    fig, ax = plt.subplots(figsize=(18, max(3, 1 + 0.6 * len(datasets))))
    ax.axis("off")
    fig.suptitle("Proteomics Validation — Multi-Contrast Overview",
                 fontsize=20, fontweight="bold", y=0.98)

    cols = ["Contrast", "Proteins", "Sig (padj<0.05)", "Overlap",
            "Dir. conc.", "Both-sig conc.", "NES (DEG up)", "Transcript comparator (C2)"]
    cell_text = df_stats[cols].values.tolist()
    colors_col = df_stats["Color"].values

    table = ax.table(cellText=cell_text, colLabels=cols, loc="center",
                     cellLoc="center", colLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 2.0)

    for j in range(len(cols)):
        cell = table[0, j]
        cell.set_text_props(fontweight="bold", fontsize=10, color="white")
        cell.set_facecolor("#2C3E50")
        cell.set_edgecolor("white")

    for i in range(len(stats)):
        for j in range(len(cols)):
            cell = table[i + 1, j]
            cell.set_edgecolor("#ECF0F1")
            if j == 0:
                cell.set_text_props(fontweight="bold", color="black", fontsize=9)
            if i % 2 == 0:
                cell.set_facecolor("#F8F9FA")
            else:
                cell.set_facecolor("white")

    fig.tight_layout()
    save_panel(fig, "figS_proteo_dataset_overview")


# ============================================================================
# MAIN
# ============================================================================
def main():
    print("=" * 70)
    print("PROTEOMICS VALIDATION - SUPPLEMENTARY FIGURE PANELS (v3)")
    print("=" * 70)
    print(f"Output directory: {OUT_DIR}")
    print()

    # panel_concordance_scatter()  # DROPPED 2026-06-20 — redundant with main fig4a
    #   (fig4a_proteomics_concordance: disease-vs-control liver|plasma, gene-set
    #    stratified, already C2). This supplement keeps the non-redundant panels.
    panel_enrichment_barplot()
    panel_effectsize_concordance()
    panel_volcano()
    panel_validation_heatmap()
    panel_bothsig_highlight()
    panel_dataset_overview()

    print()
    print("=" * 70)
    print(f"ALL PANELS SAVED TO: {OUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
