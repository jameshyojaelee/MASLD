#!/usr/bin/env python3
"""
fig4_proteomics_summary.py — Compact multi-panel proteomics validation summary.

6-panel figure:
  (a) Spearman rho per contrast (horizontal bar)
  (b) NES heatmap (contrasts x gene sets)
  (c) Direction concordance: raw vs filtered (|LFC|>0.5)
  (d) Effect-size stratified concordance (line plot)
  (e) Significant proteins per contrast (stacked: up/down)
  (f) Drug target enrichment + validation (combined)

Output: figures/fig4/fig4_proteomics_summary.pdf
"""

import pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import FancyBboxPatch
from scipy.stats import spearmanr

PROJECT_ROOT = pathlib.Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
RESULTS_DIR = PROJECT_ROOT / "Analysis" / "Proteomics" / "results"
OUT_DIR = PROJECT_ROOT / "figures" / "fig4"
OUT_DIR.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 6.5,
    "figure.dpi": 300,
    "savefig.dpi": 300,
    "savefig.bbox": "tight",
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.5,
    "ytick.major.width": 0.5,
    "xtick.major.size": 3,
    "ytick.major.size": 3,
})

# -- Contrast metadata --------------------------------------------------------
CONTRAST_ORDER = [
    "GSE276114_fibrosis",
    "PXD051911",
    "PXD051911_mash_vs_masl",
    "PXD051911_nas_high_vs_low",
    "PXD052937",
    "PXD052937_mash_vs_masl",
]
CONTRAST_LABELS = {
    "GSE276114_fibrosis":        "Zeybel\nFibrosis (Adv/Early)",
    "PXD051911":                 "Boel\nMASLD vs No-MASLD",
    "PXD051911_mash_vs_masl":    "Boel\nMASH vs MASL",
    "PXD051911_nas_high_vs_low": "Boel\nNAS high vs low",
    "PXD052937":                 "Souri.\nMASLD vs Normal",
    "PXD052937_mash_vs_masl":    "Souri.\nMASH vs MASL",
}
CONTRAST_LABELS_SHORT = {
    "GSE276114_fibrosis":        "Zeybel\nFibrosis",
    "PXD051911":                 "Boel\nMASLD",
    "PXD051911_mash_vs_masl":    "Boel\nMASH/MASL",
    "PXD051911_nas_high_vs_low": "Boel\nNAS",
    "PXD052937":                 "Souri.\nMASLD",
    "PXD052937_mash_vs_masl":    "Souri.\nMASH/MASL",
}
CONTRAST_COLORS = {
    "GSE276114_fibrosis":        "#2E86AB",
    "PXD051911":                 "#27AE60",
    "PXD051911_mash_vs_masl":    "#1E8449",
    "PXD051911_nas_high_vs_low": "#196F3D",
    "PXD052937":                 "#D64933",
    "PXD052937_mash_vs_masl":    "#B03A2E",
}
TISSUE_MAP = {
    "GSE276114_fibrosis": "Liver (SomaScan)",
    "PXD051911": "Liver (DIA-MS)",
    "PXD051911_mash_vs_masl": "Liver (DIA-MS)",
    "PXD051911_nas_high_vs_low": "Liver (DIA-MS)",
    "PXD052937": "Plasma (DIA-MS)",
    "PXD052937_mash_vs_masl": "Plasma (DIA-MS)",
}


def load():
    conc = pd.read_csv(RESULTS_DIR / "protein_transcript_concordance_v3.csv")
    diff = pd.read_csv(RESULTS_DIR / "protein_differential_results_v3.csv")
    enr = pd.read_csv(RESULTS_DIR / "protein_ranked_enrichment.csv")
    eff = pd.read_csv(RESULTS_DIR / "protein_effectsize_detection.csv")
    val = pd.read_csv(RESULTS_DIR / "protein_validation_summary.csv")
    return conc, diff, enr, eff, val


def panel_a(ax, conc):
    """Spearman rho per contrast — horizontal bars."""
    contrasts = [c for c in CONTRAST_ORDER if c in conc["dataset"].unique()]
    rhos, labels, colors = [], [], []
    for c in contrasts:
        sub = conc[conc["dataset"] == c]
        x = sub["dream_logFC"].values
        y = sub["protein_logFC"].values
        valid = np.isfinite(x) & np.isfinite(y)
        rho, _ = spearmanr(x[valid], y[valid])
        rhos.append(rho)
        labels.append(CONTRAST_LABELS_SHORT[c])
        colors.append(CONTRAST_COLORS[c])

    y_pos = np.arange(len(contrasts))
    bars = ax.barh(y_pos, rhos, height=0.65, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels)
    ax.set_xlabel("Spearman \u03c1")
    ax.set_title("(a) Protein-transcript concordance", fontweight="bold", loc="left")
    ax.set_xlim(0, max(rhos) * 1.2)
    ax.axvline(0, color="gray", lw=0.3)
    ax.invert_yaxis()

    for i, (bar, rho) in enumerate(zip(bars, rhos)):
        ax.text(rho + 0.008, bar.get_y() + bar.get_height() / 2,
                f"{rho:.3f}", va="center", fontsize=6.5, fontweight="bold")

    # Tissue annotations
    for i, c in enumerate(contrasts):
        tissue = TISSUE_MAP[c]
        ax.text(-0.02, y_pos[i], tissue, ha="right", va="center",
                fontsize=5, color="gray", transform=ax.get_yaxis_transform())


def panel_b(ax, enr):
    """NES heatmap — contrasts (rows) x gene sets (cols)."""
    gene_sets = ["dream_DEG_up", "dream_DEG_down", "Conserved", "drug_targets"]
    gs_labels = ["DEG up", "DEG down", "Conserved\nCore", "Drug\ntargets"]
    contrasts = [c for c in CONTRAST_ORDER if c in enr["dataset"].unique()]

    mat = np.full((len(contrasts), len(gene_sets)), np.nan)
    sig_mask = np.full_like(mat, False, dtype=bool)

    for i, c in enumerate(contrasts):
        for j, gs in enumerate(gene_sets):
            row = enr[(enr["dataset"] == c) & (enr["pathway"] == gs)]
            if not row.empty and pd.notna(row["NES"].values[0]):
                mat[i, j] = row["NES"].values[0]
                if pd.notna(row["padj"].values[0]) and row["padj"].values[0] < 0.05:
                    sig_mask[i, j] = True

    cmap = plt.cm.RdBu_r
    vmax = np.nanmax(np.abs(mat[np.isfinite(mat)])) if np.any(np.isfinite(mat)) else 3
    im = ax.imshow(mat, cmap=cmap, vmin=-vmax, vmax=vmax, aspect="auto")

    # Annotations
    for i in range(len(contrasts)):
        for j in range(len(gene_sets)):
            val = mat[i, j]
            if np.isnan(val):
                ax.text(j, i, "NA", ha="center", va="center", fontsize=5.5, color="gray")
            else:
                color = "white" if abs(val) > vmax * 0.6 else "black"
                weight = "bold" if sig_mask[i, j] else "normal"
                star = "*" if sig_mask[i, j] else ""
                ax.text(j, i, f"{val:.2f}{star}", ha="center", va="center",
                        fontsize=5.5, color=color, fontweight=weight)

    ax.set_xticks(range(len(gene_sets)))
    ax.set_xticklabels(gs_labels, fontsize=6.5)
    ax.set_yticks(range(len(contrasts)))
    ax.set_yticklabels([CONTRAST_LABELS_SHORT[c] for c in contrasts], fontsize=6)
    ax.set_title("(b) Ranked enrichment (NES)", fontweight="bold", loc="left")

    cbar = plt.colorbar(im, ax=ax, shrink=0.7, pad=0.04, aspect=15)
    cbar.ax.tick_params(labelsize=5.5)
    cbar.set_label("NES", fontsize=6.5)


def panel_c(ax, conc):
    """Direction concordance: raw (all genes) vs filtered (|LFC|>0.5 both)."""
    contrasts = [c for c in CONTRAST_ORDER if c in conc["dataset"].unique()]
    raw_pcts, filt_pcts, labels = [], [], []

    for c in contrasts:
        sub = conc[conc["dataset"] == c]
        n_total = len(sub)
        raw_conc = sub["direction_concordant"].sum() / n_total * 100 if n_total > 0 else 0
        filt = sub[sub["direction_concordant_filtered"] == True]
        filt_eligible = sub[(abs(sub["protein_logFC"]) > 0.5) & (abs(sub["dream_logFC"]) > 0.5)]
        filt_conc = len(filt) / len(filt_eligible) * 100 if len(filt_eligible) > 0 else 0
        raw_pcts.append(raw_conc)
        filt_pcts.append(filt_conc)
        labels.append(CONTRAST_LABELS_SHORT[c])

    x = np.arange(len(contrasts))
    w = 0.35
    bars1 = ax.bar(x - w/2, raw_pcts, w, color="#BDC3C7", edgecolor="white",
                    linewidth=0.5, label="All genes")
    bars2 = ax.bar(x + w/2, filt_pcts, w, color="#2E86AB", edgecolor="white",
                    linewidth=0.5, label="|LFC|>0.5 both")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("Direction concordance (%)")
    ax.set_ylim(0, 109)
    ax.axhline(50, ls=":", lw=0.5, c="gray", alpha=0.5)
    ax.set_title("(c) Concordance: all vs effect-filtered", fontweight="bold", loc="left")
    ax.legend(fontsize=6, loc="upper right", frameon=True, framealpha=0.9)

    for bar, val in zip(bars1, raw_pcts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f"{val:.0f}%", ha="center", va="bottom", fontsize=5, color="gray")
    for bar, val in zip(bars2, filt_pcts):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1,
                f"{val:.0f}%", ha="center", va="bottom", fontsize=5, fontweight="bold")


def panel_d(ax, eff):
    """Effect-size stratified concordance — line plot."""
    bin_order = ["0-0.5", "0.5-1", "1-1.5", "1.5-2", "2+"]
    primary_ds = ["GSE276114_fibrosis", "PXD051911", "PXD052937"]

    for ds in primary_ds:
        sub = eff[eff["dataset"] == ds].copy()
        if sub.empty:
            continue
        sub["bin_idx"] = sub["lfc_bin"].map({b: i for i, b in enumerate(bin_order)})
        sub = sub.dropna(subset=["concordance_rate"]).sort_values("bin_idx")
        ax.plot(sub["bin_idx"], sub["concordance_rate"] * 100, "-o",
                color=CONTRAST_COLORS[ds], lw=1.5, ms=4, markeredgecolor="white",
                markeredgewidth=0.3, label=CONTRAST_LABELS_SHORT[ds].replace("\n", " "))

    ax.set_xticks(range(len(bin_order)))
    ax.set_xticklabels([f"|LFC|\n{b}" for b in bin_order], fontsize=6)
    ax.set_xlabel("Dream |logFC| bin")
    ax.set_ylabel("Direction concordance (%)")
    ax.set_ylim(40, 105)
    ax.axhline(50, ls=":", lw=0.5, c="gray", alpha=0.5)
    ax.set_title("(d) Concordance by effect size", fontweight="bold", loc="left")
    ax.legend(fontsize=5.5, loc="lower right", frameon=True, framealpha=0.9)
    ax.grid(True, alpha=0.1, ls="--")


def panel_e(ax, diff):
    """Significant proteins per contrast — stacked up/down bars."""
    contrasts = [c for c in CONTRAST_ORDER if c in diff["dataset"].unique()]
    n_up, n_down, labels, colors = [], [], [], []

    for c in contrasts:
        sub = diff[diff["dataset"] == c]
        sig = sub[sub["padj"] < 0.05]
        n_up.append((sig["logFC"] > 0).sum())
        n_down.append((sig["logFC"] < 0).sum())
        labels.append(CONTRAST_LABELS_SHORT[c])
        colors.append(CONTRAST_COLORS[c])

    x = np.arange(len(contrasts))
    ax.bar(x, n_up, 0.6, color=[mcolors.to_rgba(c, 0.9) for c in colors],
           edgecolor="white", linewidth=0.5, label="Up")
    ax.bar(x, [-v for v in n_down], 0.6,
           color=[mcolors.to_rgba(c, 0.4) for c in colors],
           edgecolor="white", linewidth=0.5, label="Down")

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=6)
    ax.set_ylabel("Significant proteins (padj<0.05)")
    ax.axhline(0, color="black", lw=0.5)
    ax.set_title("(e) Differential protein abundance", fontweight="bold", loc="left")

    # Total annotations
    for i, (u, d) in enumerate(zip(n_up, n_down)):
        total = u + d
        ax.text(i, u + max(n_up) * 0.03, str(u), ha="center", va="bottom",
                fontsize=5.5, fontweight="bold")
        ax.text(i, -d - max(n_down) * 0.03, str(d), ha="center", va="top",
                fontsize=5.5, color="gray")

    ax.text(0.98, 0.95, "Up", transform=ax.transAxes, ha="right", va="top",
            fontsize=6, color="#C0392B", fontweight="bold")
    ax.text(0.98, 0.05, "Down", transform=ax.transAxes, ha="right", va="bottom",
            fontsize=6, color="#2E86AB", fontweight="bold")


def panel_f(ax, enr, val):
    """Drug target NES + validation recovery — combined dot + bar."""
    contrasts = [c for c in CONTRAST_ORDER if c in enr["dataset"].unique()]

    # Drug target NES
    nes_vals, nes_sig = [], []
    for c in contrasts:
        row = enr[(enr["dataset"] == c) & (enr["pathway"] == "drug_targets")]
        if not row.empty and pd.notna(row["NES"].values[0]):
            nes_vals.append(row["NES"].values[0])
            nes_sig.append(row["padj"].values[0] < 0.05 if pd.notna(row["padj"].values[0]) else False)
        else:
            nes_vals.append(0)
            nes_sig.append(False)

    # Drug target sig rate from validation
    sig_rates = []
    for c in contrasts:
        row = val[(val["dataset"] == c) & (val["validation_set"] == "drug_targets")]
        if not row.empty:
            sig_rates.append(row["sig_rate"].values[0] * 100)
        else:
            sig_rates.append(0)

    x = np.arange(len(contrasts))
    colors = [CONTRAST_COLORS[c] for c in contrasts]

    # NES bars
    bars = ax.bar(x, nes_vals, 0.55, color=colors, alpha=0.7, edgecolor="white", linewidth=0.5)
    for i, (v, sig) in enumerate(zip(nes_vals, nes_sig)):
        marker = "*" if sig else ""
        ax.text(i, v + 0.05, f"{v:.1f}{marker}", ha="center", va="bottom",
                fontsize=5.5, fontweight="bold" if sig else "normal")

    ax.set_ylabel("Drug target NES", color="#2C3E50")
    ax.set_ylim(0, max(nes_vals) * 1.3 if nes_vals else 3)

    # Overlay: sig rate as diamonds on secondary axis
    ax2 = ax.twinx()
    ax2.plot(x, sig_rates, "D-", color="#E67E22", ms=5, lw=1.2,
             markeredgecolor="white", markeredgewidth=0.5, zorder=5)
    ax2.set_ylabel("Drug targets sig (%)", color="#E67E22", fontsize=7)
    ax2.set_ylim(0, max(sig_rates) * 1.3 if sig_rates else 60)
    ax2.tick_params(axis="y", colors="#E67E22")
    ax2.spines["right"].set_visible(True)
    ax2.spines["right"].set_color("#E67E22")
    ax2.spines["right"].set_linewidth(0.6)

    ax.set_xticks(x)
    ax.set_xticklabels([CONTRAST_LABELS_SHORT[c] for c in contrasts], fontsize=6)
    ax.set_title("(f) Drug target enrichment & significance", fontweight="bold", loc="left")
    ax.axhline(0, color="gray", lw=0.3)

    # Legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], marker="s", color="w", markerfacecolor="#7F8C8D",
               markersize=6, label="NES (bars)"),
        Line2D([0], [0], marker="D", color="#E67E22", lw=1,
               markersize=4, label="Sig rate (%)"),
    ]
    ax.legend(handles=legend_elements, fontsize=5.5, loc="upper left", framealpha=0.9)


def save_individual_panel(panel_func, panel_args, name, panel_dir, figsize=(5, 3.5)):
    """Save a single panel as standalone PDF."""
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    panel_func(ax, *panel_args)
    fig.tight_layout()
    out = panel_dir / f"{name}.pdf"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"  -> {name}")


def main():
    print("Loading data...")
    conc, diff, enr, eff, val = load()

    # -- Combined figure ---------------------------------------------------
    fig = plt.figure(figsize=(10, 7.5), constrained_layout=False)
    gs = fig.add_gridspec(3, 2, hspace=0.55, wspace=0.45,
                          left=0.08, right=0.95, top=0.94, bottom=0.06)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    ax_e = fig.add_subplot(gs[2, 0])
    ax_f = fig.add_subplot(gs[2, 1])

    print("Drawing panels...")
    panel_a(ax_a, conc)
    panel_b(ax_b, enr)
    panel_c(ax_c, conc)
    panel_d(ax_d, eff)
    panel_e(ax_e, diff)
    panel_f(ax_f, enr, val)

    fig.suptitle("Multi-Contrast Proteomics Validation of Transcriptomic Atlas",
                 fontsize=11, fontweight="bold", y=0.98)

    out = OUT_DIR / "fig4_proteomics_summary.pdf"
    fig.savefig(out, dpi=300)
    plt.close(fig)
    print(f"Saved: {out}")

    # -- Individual panels -------------------------------------------------
    panel_dir = OUT_DIR / "panels"
    panel_dir.mkdir(exist_ok=True)
    print(f"\nSaving individual panels to {panel_dir}/")

    save_individual_panel(panel_a, (conc,), "panel_a_concordance_rho", panel_dir, figsize=(5, 3.5))
    save_individual_panel(panel_b, (enr,), "panel_b_nes_heatmap", panel_dir, figsize=(5.5, 3.5))
    save_individual_panel(panel_c, (conc,), "panel_c_direction_concordance", panel_dir, figsize=(5.5, 3.5))
    save_individual_panel(panel_d, (eff,), "panel_d_effectsize_concordance", panel_dir, figsize=(5, 3.5))
    save_individual_panel(panel_e, (diff,), "panel_e_significant_proteins", panel_dir, figsize=(5.5, 3.5))
    save_individual_panel(panel_f, (enr, val), "panel_f_drug_target_enrichment", panel_dir, figsize=(5.5, 3.5))

    print(f"\nAll panels saved to: {panel_dir}")


if __name__ == "__main__":
    main()
