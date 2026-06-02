#!/usr/bin/env python3
"""
proteomics_presentation_panels.py — Presentation-quality proteomics validation panels.

Generates individual high-resolution PDF panels AND a combined composite:
  1. Protein–transcript concordance scatter (3 datasets)
  2. Ranked enrichment barplot (NES across gene sets × datasets)
  3. Effect-size stratified concordance (monotonic LFC trend)
  4. Volcano plots (3 datasets, side-by-side)
  5. Drug target & positive control detection heatmap
  6. Both-significant concordance highlight scatter
  7. Dataset overview / sample comparison

Outputs to: Analysis/Proteomics/results/presentation_panels/

SLURM: --partition=cpu --cpus=4 --mem=16G --time=48:00:00
"""

import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from scipy.stats import spearmanr
from adjustText import adjust_text
import warnings
warnings.filterwarnings("ignore")

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
PROTEO_DIR = PROJECT_ROOT / "Analysis" / "Proteomics"
RESULTS_DIR = PROTEO_DIR / "results"
OUT_DIR = RESULTS_DIR / "presentation_panels"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# ── Presentation-quality rcParams ──────────────────────────────────────────────
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

# ── Color palette ──────────────────────────────────────────────────────────────
DATASET_COLORS = {
    "GSE276114_fibrosis": "#2E86AB",   # Zeybel (liver SomaScan)
    "PXD052937": "#D64933",            # Sourianarayanane (plasma DIA-MS)
    "PXD051911": "#27AE60",            # Boel (liver DIA-MS)
}
DATASET_LABELS = {
    "GSE276114_fibrosis": "Zeybel (Liver SomaScan)",
    "PXD052937": "Sourianarayanane (Plasma DIA-MS)",
    "PXD051911": "Boel (Liver DIA-MS)",
}
GENESET_COLORS = {
    "dream_DEG_up": "#C0392B",
    "dream_DEG_down": "#2E86AB",
    "dream_DEG_all": "#7F8C8D",
    "Conserved": "#F39C12",
    "drug_targets": "#8E44AD",
}


def save_panel(fig, name, formats=("pdf",)):
    """Save figure panel."""
    for fmt in formats:
        fig.savefig(OUT_DIR / f"{name}.{fmt}", format=fmt)
    plt.close(fig)
    print(f"  ✓ Saved: {name}")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 1: Protein–transcript direction concordance scatter
# ═══════════════════════════════════════════════════════════════════════════════
def panel_01_concordance_scatter():
    """Protein logFC vs transcript logFC — scatter per dataset with rho annotation."""
    print("[Panel 1] Protein–transcript concordance scatter")
    conc = pd.read_csv(RESULTS_DIR / "protein_transcript_concordance_v2.csv")
    datasets = [ds for ds in conc["dataset"].unique() if ds in DATASET_COLORS]

    fig, axes = plt.subplots(1, len(datasets), figsize=(6 * len(datasets), 6))
    if len(datasets) == 1:
        axes = [axes]
    fig.suptitle("Protein–Transcript Direction Concordance", fontsize=22, fontweight="bold", y=1.03)

    for i, ds in enumerate(datasets):
        ax = axes[i]
        sub = conc[conc["dataset"] == ds].copy()
        x = sub["dream_logFC"].values
        y = sub["protein_logFC"].values
        valid = np.isfinite(x) & np.isfinite(y)
        x, y = x[valid], y[valid]

        # Both significant
        both_sig = (sub["protein_padj"].values[valid] < 0.05) & (sub["dream_padj"].values[valid] < 0.05)
        concordant = sub["direction_concordant"].values[valid]

        # Background: all genes
        ax.scatter(x[~both_sig], y[~both_sig], s=6, alpha=0.15, c="#BDC3C7",
                   edgecolors="none", rasterized=True, zorder=1)
        # Both significant: colored by concordance
        ax.scatter(x[both_sig & concordant], y[both_sig & concordant],
                   s=12, alpha=0.5, c=DATASET_COLORS[ds],
                   edgecolors="none", rasterized=True, zorder=2, label="Concordant")
        ax.scatter(x[both_sig & ~concordant], y[both_sig & ~concordant],
                   s=12, alpha=0.5, c="#E74C3C", marker="x", linewidths=0.8,
                   zorder=3, label="Discordant")

        # Regression line
        z = np.polyfit(x, y, 1)
        xline = np.linspace(x.min(), x.max(), 100)
        ax.plot(xline, np.polyval(z, xline), "--", c="black", lw=1.2, alpha=0.6)

        rho, p = spearmanr(x, y)
        n_both = both_sig.sum()
        n_conc_both = (both_sig & concordant).sum()
        pct = 100 * n_conc_both / n_both if n_both > 0 else 0

        ax.axhline(0, ls=":", lw=0.6, c="gray", alpha=0.4)
        ax.axvline(0, ls=":", lw=0.6, c="gray", alpha=0.4)

        ax.set_xlabel("Transcript logFC (dream)", fontsize=15)
        if i == 0:
            ax.set_ylabel("Protein logFC", fontsize=15)
        ax.set_title(DATASET_LABELS[ds], fontsize=14, fontweight="bold",
                     color=DATASET_COLORS[ds])

        stats_text = (f"ρ = {rho:.3f}\n"
                      f"n = {len(x):,}\n"
                      f"Both sig: {n_both:,}\n"
                      f"Concordance: {pct:.1f}%")
        ax.annotate(stats_text, xy=(0.04, 0.96), xycoords="axes fraction",
                    fontsize=11, va="top", fontweight="bold",
                    bbox=dict(boxstyle="round,pad=0.4", fc="white",
                              ec=DATASET_COLORS[ds], alpha=0.9, lw=1.5))
        ax.legend(fontsize=9, loc="lower right", frameon=True, framealpha=0.9)
        ax.grid(True, alpha=0.1, ls="--")

    fig.tight_layout()
    save_panel(fig, "panel_01_concordance_scatter")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 2: Ranked enrichment barplot (NES)
# ═══════════════════════════════════════════════════════════════════════════════
def panel_02_enrichment_barplot():
    """FGSEA NES for dream DEGs, Conserved, drug targets across datasets."""
    print("[Panel 2] Ranked enrichment barplot")
    enr = pd.read_csv(RESULTS_DIR / "protein_ranked_enrichment.csv")
    enr = enr.dropna(subset=["NES"])

    # Prepare display
    gene_sets = ["dream_DEG_up", "dream_DEG_down", "Conserved", "drug_targets"]
    gene_set_labels = {
        "dream_DEG_up": "Dream DEGs (Up)",
        "dream_DEG_down": "Dream DEGs (Down)",
        "Conserved": "Conserved Core",
        "drug_targets": "Drug Targets",
    }
    datasets = [ds for ds in ["GSE276114_fibrosis", "PXD051911", "PXD052937"]
                if ds in enr["dataset"].unique()]

    fig, ax = plt.subplots(figsize=(12, 7))
    fig.suptitle("Proteomics Ranked Enrichment of Transcriptomic Signatures",
                 fontsize=20, fontweight="bold", y=1.02)

    bar_width = 0.22
    x_base = np.arange(len(gene_sets))

    for di, ds in enumerate(datasets):
        sub = enr[enr["dataset"] == ds]
        nes_vals = []
        sig_markers = []
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
        bars = ax.bar(x_base + offset, nes_vals, bar_width * 0.9,
                      color=DATASET_COLORS[ds], alpha=0.85,
                      edgecolor="white", linewidth=0.8,
                      label=DATASET_LABELS[ds])

        # Significance stars above bars
        for xi, (val, sig) in enumerate(zip(nes_vals, sig_markers)):
            if sig and val != 0:
                y_pos = val + 0.08 * np.sign(val)
                ax.text(x_base[xi] + offset, y_pos, sig,
                        ha="center", va="bottom" if val > 0 else "top",
                        fontsize=12, fontweight="bold", color=DATASET_COLORS[ds])

    ax.set_xticks(x_base)
    ax.set_xticklabels([gene_set_labels.get(gs, gs) for gs in gene_sets],
                       fontsize=13, fontweight="bold")
    ax.set_ylabel("Normalized Enrichment Score (NES)", fontsize=15)
    ax.axhline(0, color="black", lw=0.8)
    ax.legend(fontsize=11, loc="upper right", frameon=True, framealpha=0.9)
    ax.grid(True, axis="y", alpha=0.15, ls="--")

    # Annotate NES direction
    ax.text(0.02, 0.98, "Enriched in\nprotein ↑", transform=ax.transAxes,
            fontsize=10, va="top", color="#C0392B", fontweight="bold")
    ax.text(0.02, 0.02, "Enriched in\nprotein ↓", transform=ax.transAxes,
            fontsize=10, va="bottom", color="#2E86AB", fontweight="bold")

    fig.tight_layout()
    save_panel(fig, "panel_02_enrichment_barplot")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 3: Effect-size stratified concordance
# ═══════════════════════════════════════════════════════════════════════════════
def panel_03_effectsize_concordance():
    """Concordance rate vs transcript |logFC| magnitude — monotonic trend."""
    print("[Panel 3] Effect-size stratified concordance")
    eff = pd.read_csv(RESULTS_DIR / "protein_effectsize_detection.csv")

    bin_order = ["0-0.5", "0.5-1", "1-1.5", "1.5-2", "2+"]
    datasets = [ds for ds in eff["dataset"].unique() if ds in DATASET_COLORS]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    fig.suptitle("Protein–Transcript Concordance by Effect Size",
                 fontsize=20, fontweight="bold", y=1.03)

    # (a) Concordance rate
    ax1 = axes[0]
    for ds in datasets:
        sub = eff[eff["dataset"] == ds].copy()
        sub["bin_idx"] = sub["lfc_bin"].map({b: i for i, b in enumerate(bin_order)})
        sub = sub.dropna(subset=["concordance_rate"]).sort_values("bin_idx")
        ax1.plot(sub["bin_idx"], sub["concordance_rate"] * 100, "-o",
                 color=DATASET_COLORS[ds], lw=2.5, ms=8,
                 label=DATASET_LABELS[ds], zorder=3)

    ax1.set_xticks(range(len(bin_order)))
    ax1.set_xticklabels([f"|LFC| {b}" for b in bin_order], fontsize=11, rotation=15)
    ax1.set_xlabel("Transcript Effect Size (|logFC| bin)", fontsize=14)
    ax1.set_ylabel("Direction Concordance (%)", fontsize=14)
    ax1.set_title("Concordance Increases with Effect Size", fontsize=15, fontweight="bold")
    ax1.set_ylim(40, 105)
    ax1.axhline(50, ls=":", lw=1, c="gray", alpha=0.5, label="Random (50%)")
    ax1.legend(fontsize=9, loc="lower right", frameon=True, framealpha=0.9)
    ax1.grid(True, alpha=0.15, ls="--")

    # (b) Detection rate
    ax2 = axes[1]
    for ds in datasets:
        sub = eff[eff["dataset"] == ds].copy()
        sub["bin_idx"] = sub["lfc_bin"].map({b: i for i, b in enumerate(bin_order)})
        sub = sub.sort_values("bin_idx")
        ax2.bar([i + list(datasets).index(ds) * 0.25 - 0.25
                 for i in sub["bin_idx"]],
                sub["n_detected"].values,
                width=0.22, color=DATASET_COLORS[ds], alpha=0.8,
                edgecolor="white", label=DATASET_LABELS[ds])

    ax2.set_xticks(range(len(bin_order)))
    ax2.set_xticklabels([f"|LFC| {b}" for b in bin_order], fontsize=11, rotation=15)
    ax2.set_xlabel("Transcript Effect Size (|logFC| bin)", fontsize=14)
    ax2.set_ylabel("Proteins Detected (n)", fontsize=14)
    ax2.set_title("Protein Detection by Effect Size", fontsize=15, fontweight="bold")
    ax2.legend(fontsize=9, loc="upper right", frameon=True, framealpha=0.9)
    ax2.grid(True, axis="y", alpha=0.15, ls="--")

    fig.tight_layout()
    save_panel(fig, "panel_03_effectsize_concordance")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 4: Volcano plots (3 datasets)
# ═══════════════════════════════════════════════════════════════════════════════
def panel_04_volcano():
    """Protein differential expression volcano plots per dataset."""
    print("[Panel 4] Volcano plots")
    diff = pd.read_csv(RESULTS_DIR / "protein_differential_results_v2.csv")
    datasets = [ds for ds in diff["dataset"].unique() if ds in DATASET_COLORS]

    fig, axes = plt.subplots(1, len(datasets), figsize=(6 * len(datasets), 6))
    if len(datasets) == 1:
        axes = [axes]
    fig.suptitle("Protein Differential Abundance", fontsize=22, fontweight="bold", y=1.03)

    # Known MASLD genes to highlight
    highlight = ["THRB", "NR1H4", "PPARA", "PPARG", "PNPLA3", "HSD17B13",
                 "FASN", "SCD", "A2M", "FGB", "CRP", "HP", "TGFB1",
                 "COL1A1", "ACTA2", "APOB", "ALB", "CYP2E1"]

    for i, ds in enumerate(datasets):
        ax = axes[i]
        sub = diff[diff["dataset"] == ds].copy()
        sub["-log10p"] = -np.log10(sub["pvalue"].clip(lower=1e-300))

        # Categories
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

        # Highlight known genes
        for gene in highlight:
            row = sub[sub["gene"] == gene]
            if not row.empty and row["padj"].values[0] < 0.05:
                gx, gy = row["logFC"].values[0], row["-log10p"].values[0]
                ax.annotate(gene, (gx, gy), fontsize=9, fontweight="bold",
                            fontstyle="italic", textcoords="offset points",
                            xytext=(6, 4),
                            arrowprops=dict(arrowstyle="->", color="black", lw=0.7),
                            zorder=5)

        ax.axhline(-np.log10(0.05), ls="--", lw=0.8, c="gray", alpha=0.5)
        ax.axvline(0, ls=":", lw=0.6, c="gray", alpha=0.4)
        ax.set_xlabel("Protein logFC", fontsize=14)
        if i == 0:
            ax.set_ylabel("-log₁₀(p-value)", fontsize=14)
        ax.set_title(DATASET_LABELS[ds], fontsize=14, fontweight="bold",
                     color=DATASET_COLORS[ds])
        ax.legend(fontsize=10, loc="upper right", frameon=True, framealpha=0.9)
        ax.grid(True, alpha=0.1, ls="--")

    fig.tight_layout()
    save_panel(fig, "panel_04_volcano")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 5: Drug target & positive control detection heatmap
# ═══════════════════════════════════════════════════════════════════════════════
def panel_05_validation_heatmap():
    """Detection and significance rates for positive controls and drug targets."""
    print("[Panel 5] Validation heatmap")
    val = pd.read_csv(RESULTS_DIR / "protein_validation_summary.csv")

    datasets = [ds for ds in val["dataset"].unique() if ds in DATASET_COLORS]
    val_sets = ["positive_controls", "drug_targets"]
    val_labels = {"positive_controls": "Positive Controls\n(n=67)",
                  "drug_targets": "Drug Targets\n(n=102)"}

    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    fig.suptitle("Proteomics Validation of Key Gene Sets", fontsize=20,
                 fontweight="bold", y=1.03)

    for vi, vs in enumerate(val_sets):
        ax = axes[vi]
        sub = val[val["validation_set"] == vs]

        x = np.arange(len(datasets))
        det_rates = []
        sig_rates = []
        for ds in datasets:
            row = sub[sub["dataset"] == ds]
            det_rates.append(row["detection_rate"].values[0] * 100 if len(row) > 0 else 0)
            sig_rates.append(row["sig_rate"].values[0] * 100 if len(row) > 0 else 0)

        bars1 = ax.bar(x - 0.18, det_rates, 0.32, color=[DATASET_COLORS[ds] for ds in datasets],
                       alpha=0.4, edgecolor="white", label="Detected")
        bars2 = ax.bar(x + 0.18, sig_rates, 0.32, color=[DATASET_COLORS[ds] for ds in datasets],
                       alpha=0.9, edgecolor="white", label="Significant (padj<0.05)")

        # Value labels on bars
        for bar, val_num in zip(bars1, det_rates):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                    f"{val_num:.0f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")
        for bar, val_num in zip(bars2, sig_rates):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 1,
                    f"{val_num:.0f}%", ha="center", va="bottom", fontsize=10, fontweight="bold")

        ax.set_xticks(x)
        ax.set_xticklabels([DATASET_LABELS[ds].split("(")[0].strip() for ds in datasets],
                           fontsize=11, rotation=15, ha="right")
        ax.set_ylabel("Rate (%)", fontsize=14)
        ax.set_title(val_labels[vs], fontsize=15, fontweight="bold")
        ax.set_ylim(0, 110)
        ax.legend(fontsize=10, loc="upper right")
        ax.grid(True, axis="y", alpha=0.15, ls="--")

    fig.tight_layout()
    save_panel(fig, "panel_05_validation_heatmap")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 6: Both-significant gene concordance highlight
# ═══════════════════════════════════════════════════════════════════════════════
def panel_06_bothsig_highlight():
    """Scatter of both-significant genes with key drug targets annotated."""
    print("[Panel 6] Both-significant concordance highlight")
    conc = pd.read_csv(RESULTS_DIR / "protein_transcript_concordance_v2.csv")

    # Focus on GSE276114 (largest overlap)
    ds = "GSE276114_fibrosis"
    sub = conc[conc["dataset"] == ds].copy()
    both_sig = sub[(sub["protein_padj"] < 0.05) & (sub["dream_padj"] < 0.05)].copy()

    if len(both_sig) < 10:
        print("  SKIP: Too few both-significant genes")
        return

    fig, ax = plt.subplots(figsize=(9, 8))

    x = both_sig["dream_logFC"].values
    y = both_sig["protein_logFC"].values
    conc_mask = both_sig["direction_concordant"].values

    ax.scatter(x[conc_mask], y[conc_mask], s=20, alpha=0.5, c="#2E86AB",
               edgecolors="none", rasterized=True, label=f"Concordant (n={conc_mask.sum():,})")
    ax.scatter(x[~conc_mask], y[~conc_mask], s=25, alpha=0.7, c="#E74C3C",
               marker="x", linewidths=1, label=f"Discordant (n={(~conc_mask).sum():,})")

    # Highlight key genes with manual offsets to avoid overlap
    drug_targets = ["THRB", "NR1H4", "PPARA", "PPARG", "PNPLA3", "HSD17B13",
                    "FASN", "SCD", "CYP2E1", "TGFB1", "COL1A1", "ALB",
                    "A2M", "HP", "ACTA2", "APOB", "FGB"]
    # Manual offsets (dx, dy in points) — fan out the tight bottom-left cluster
    # Cluster: THRB(-0.27,-0.52) NR1H4(-0.11,-0.55) PPARA(-0.23,-0.64) PNPLA3(-0.28,-0.66) CYP2E1(-0.35,-0.99)
    label_offsets = {
        "THRB":    (-70, 25),     # upper-left
        "NR1H4":   (60, 25),      # upper-right
        "PPARA":   (70, -10),     # right
        "PNPLA3":  (-75, -10),    # left
        "CYP2E1":  (-55, -35),    # lower-left
    }
    for gene in drug_targets:
        row = both_sig[both_sig["gene"] == gene]
        if not row.empty:
            gx, gy = row["dream_logFC"].values[0], row["protein_logFC"].values[0]
            ax.scatter([gx], [gy], s=60, c="#F39C12", edgecolors="black",
                       linewidths=0.8, zorder=5)
            offset = label_offsets.get(gene, (8, 5))
            ax.annotate(gene, (gx, gy), fontsize=10, fontweight="bold",
                        fontstyle="italic", textcoords="offset points",
                        xytext=offset,
                        arrowprops=dict(arrowstyle="->", color="black", lw=0.8),
                        zorder=6)

    # Quadrant shading
    ax.axhspan(0, ax.get_ylim()[1], xmin=0.5, alpha=0.03, color="#2ECC71", zorder=0)
    ax.axhspan(ax.get_ylim()[0], 0, xmax=0.5, alpha=0.03, color="#2ECC71", zorder=0)
    ax.axhline(0, ls=":", lw=0.8, c="gray", alpha=0.4)
    ax.axvline(0, ls=":", lw=0.8, c="gray", alpha=0.4)

    rho, p = spearmanr(x, y)
    ax.set_xlabel("Transcript logFC", fontsize=15)
    ax.set_ylabel("Protein logFC", fontsize=15)
    ax.set_title(f"Both-Significant Genes (n={len(both_sig):,})\n"
                 f"Spearman ρ = {rho:.3f}, p = {p:.1e}",
                 fontsize=16, fontweight="bold")
    ax.legend(fontsize=11, loc="lower right", frameon=True, framealpha=0.9)
    ax.grid(True, alpha=0.1, ls="--")

    # Quadrant labels
    ax.text(0.97, 0.97, "Both ↑", transform=ax.transAxes, fontsize=12,
            ha="right", va="top", color="#27AE60", fontweight="bold", alpha=0.6)
    ax.text(0.03, 0.03, "Both ↓", transform=ax.transAxes, fontsize=12,
            ha="left", va="bottom", color="#27AE60", fontweight="bold", alpha=0.6)
    ax.text(0.97, 0.03, "Discordant", transform=ax.transAxes, fontsize=12,
            ha="right", va="bottom", color="#E74C3C", fontweight="bold", alpha=0.4)
    ax.text(0.03, 0.97, "Discordant", transform=ax.transAxes, fontsize=12,
            ha="left", va="top", color="#E74C3C", fontweight="bold", alpha=0.4)

    fig.tight_layout()
    save_panel(fig, "panel_06_bothsig_highlight")


# ═══════════════════════════════════════════════════════════════════════════════
# PANEL 7: Multi-dataset overview summary
# ═══════════════════════════════════════════════════════════════════════════════
def panel_07_dataset_overview():
    """Summary comparison of 3 proteomics datasets — platform, tissue, key stats."""
    print("[Panel 7] Dataset overview")
    diff = pd.read_csv(RESULTS_DIR / "protein_differential_results_v2.csv")
    conc = pd.read_csv(RESULTS_DIR / "protein_transcript_concordance_v2.csv")
    enr = pd.read_csv(RESULTS_DIR / "protein_ranked_enrichment.csv")

    datasets = ["GSE276114_fibrosis", "PXD051911", "PXD052937"]

    # Gather stats per dataset
    stats = []
    for ds in datasets:
        d_diff = diff[diff["dataset"] == ds]
        d_conc = conc[conc["dataset"] == ds]
        d_enr = enr[(enr["dataset"] == ds) & (enr["pathway"] == "dream_DEG_up")]

        n_proteins = len(d_diff)
        n_sig = (d_diff["padj"] < 0.05).sum()
        n_overlap = len(d_conc)
        overall_conc = d_conc["direction_concordant"].mean() * 100
        both_sig = d_conc[(d_conc["protein_padj"] < 0.05) & (d_conc["dream_padj"] < 0.05)]
        bothsig_conc = both_sig["direction_concordant"].mean() * 100 if len(both_sig) > 0 else np.nan
        nes = d_enr["NES"].values[0] if len(d_enr) > 0 else np.nan

        stats.append({
            "Dataset": DATASET_LABELS[ds],
            "Proteins": n_proteins,
            "Sig (padj<0.05)": n_sig,
            "Overlap w/ dream": n_overlap,
            "Direction conc.": f"{overall_conc:.1f}%",
            "Both-sig conc.": f"{bothsig_conc:.1f}%" if not np.isnan(bothsig_conc) else "N/A",
            "NES (DEG up)": f"{nes:.2f}" if not np.isnan(nes) else "N/A",
            "Color": DATASET_COLORS[ds],
        })

    df_stats = pd.DataFrame(stats)

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.axis("off")
    fig.suptitle("Proteomics Dataset Comparison", fontsize=22, fontweight="bold", y=0.98)

    cols = ["Dataset", "Proteins", "Sig (padj<0.05)", "Overlap w/ dream",
            "Direction conc.", "Both-sig conc.", "NES (DEG up)"]
    cell_text = df_stats[cols].values.tolist()
    colors_col = df_stats["Color"].values

    table = ax.table(cellText=cell_text, colLabels=cols, loc="center",
                     cellLoc="center", colLoc="center")
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1, 2.2)

    # Style header
    for j in range(len(cols)):
        cell = table[0, j]
        cell.set_text_props(fontweight="bold", fontsize=12, color="white")
        cell.set_facecolor("#2C3E50")
        cell.set_edgecolor("white")

    # Style data rows
    for i in range(len(stats)):
        for j in range(len(cols)):
            cell = table[i + 1, j]
            cell.set_edgecolor("#ECF0F1")
            if j == 0:
                cell.set_text_props(fontweight="bold", color=colors_col[i])
            if i % 2 == 0:
                cell.set_facecolor("#F8F9FA")
            else:
                cell.set_facecolor("white")

    fig.tight_layout()
    save_panel(fig, "panel_07_dataset_overview")


# ═══════════════════════════════════════════════════════════════════════════════
# COMPOSITE
# ═══════════════════════════════════════════════════════════════════════════════
def make_composite():
    """Assemble panels into a composite overview."""
    print("\n[Composite] Assembling overview...")
    try:
        from PIL import Image
    except ImportError:
        print("  SKIP: PIL not available")
        return

    panel_files = sorted(OUT_DIR.glob("panel_*.png"))
    if len(panel_files) < 3:
        # Fall back to re-rendering as PNG for composite
        print("  No PNGs available for composite (PDF-only mode)")
        return

    imgs = []
    for pf in panel_files:
        try:
            imgs.append((pf.stem, Image.open(pf)))
        except Exception:
            pass
    if not imgs:
        return

    n_cols = 2
    n_rows = (len(imgs) + 1) // 2
    cell_w, cell_h, margin = 2400, 1600, 60
    canvas_w = n_cols * cell_w + (n_cols + 1) * margin
    canvas_h = n_rows * cell_h + (n_rows + 1) * margin + 120
    canvas = Image.new("RGB", (canvas_w, canvas_h), "white")

    for idx, (name, img) in enumerate(imgs):
        row, col = idx // n_cols, idx % n_cols
        scale = min(cell_w / img.width, cell_h / img.height)
        new_w, new_h = int(img.width * scale), int(img.height * scale)
        img_resized = img.resize((new_w, new_h), Image.LANCZOS)
        x = col * cell_w + (col + 1) * margin + (cell_w - new_w) // 2
        y = row * cell_h + (row + 1) * margin + 120 + (cell_h - new_h) // 2
        canvas.paste(img_resized, (x, y))

    canvas.save(OUT_DIR / "proteomics_composite.png", dpi=(300, 300))
    print(f"  ✓ Saved: proteomics_composite.png")


# ═══════════════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 70)
    print("PROTEOMICS VALIDATION — PRESENTATION PANELS")
    print("=" * 70)
    print(f"Output directory: {OUT_DIR}")
    print()

    panel_01_concordance_scatter()
    panel_02_enrichment_barplot()
    panel_03_effectsize_concordance()
    panel_04_volcano()
    panel_05_validation_heatmap()
    panel_06_bothsig_highlight()
    panel_07_dataset_overview()

    make_composite()

    print()
    print("=" * 70)
    print(f"ALL PANELS SAVED TO: {OUT_DIR}")
    print("=" * 70)


if __name__ == "__main__":
    main()
