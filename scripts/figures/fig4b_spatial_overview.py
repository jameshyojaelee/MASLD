#!/usr/bin/env python3
"""
Fig 4B — Spatial overview (two panels):
  Left:   Representative tissue sections (Healthy | Steatohepatitis) with
          per-spot periportal→pericentral zone annotation
  Right:  Disease-emergent SVG score distribution across ALL spots from
          both cohorts (GSE192741 + Vu et al.), split by condition

Output: figures/main/fig4_validation/panels/fig4b_spatial_overview.pdf
Environment: spatial (conda)
"""
import os, sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import ListedColormap
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import scanpy as sc

BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
)
sc.settings.verbosity = 0

# Donor-aware helpers (ensure_lognorm for raw-count score_genes).
sys.path.insert(0, os.path.join(BASE, "Analysis/Spatial/scripts"))
from spatial_stats import ensure_lognorm

OUT_DIR      = os.path.join(BASE, "figures/main/fig4_validation/panels")
DATA_DIR     = os.path.join(BASE, "Analysis/Spatial/data/gsmap_input")
SVG_CSV      = os.path.join(BASE, "Analysis/Spatial/results/integration/spatial_consensus.csv")
# Real disease-emergent SVG call (F250/F196): genes that gained spatial structure
# in disease (category == disease_emergent_SVG), NOT a generic Moran's-I threshold.
DSVG_CSV     = os.path.join(BASE, "Analysis/Spatial/results/svg/differential_svgs.csv")
ZONE_CSV     = os.path.join(BASE, "Analysis/Spatial/results/zonation/zonation_scores.csv")

# ── Publication color constants (matching R theme) ────────────────────────────
ZONE_COLORS = {
    "PP1": "#1565C0",   # periportal — deep blue
    "PP2": "#5E92C8",
    "Mid": "#9E9E9E",   # mid-lobular — gray
    "PC2": "#E07B39",
    "PC1": "#E65100",   # pericentral — deep orange
}
COND_COLORS = {
    "Healthy":         "#9E9E9E",   # control gray
    "Steatohepatitis": "#C9265E",   # disease magenta
    "Vu (MASLD)":      "#E91E63",   # Tier2 bright magenta
}

# ── Disease-emergent SVG gene list ────────────────────────────────────────────
# F250/F196: use the actual disease-emergent SVG call (genes that gained spatial
# autocorrelation in MASLD vs healthy; category == disease_emergent_SVG, 184
# genes), ranked by delta_I so the top_n scored genes are the strongest. The
# previous definition (generic Moran's-I > 0.03 in >=1 cohort, alphabetical) was
# a generic spatial-variability signature with no disease direction.
dsvg_df = pd.read_csv(DSVG_CSV)
gene_col = "gene" if "gene" in dsvg_df.columns else dsvg_df.columns[0]
disease_emergent = (
    dsvg_df[dsvg_df["category"] == "disease_emergent_SVG"]
    .sort_values("delta_I", ascending=False)[gene_col]
    .tolist()
)

# ── Load zonation scores (GSE192741 only) ─────────────────────────────────────
zone_df = pd.read_csv(ZONE_CSV, index_col=0)
zone_bin_order = ["PP1", "PP2", "Mid", "PC2", "PC1"]

# ── Helper: compute SVG score for one adata ───────────────────────────────────
def compute_svg_score(adata, gene_list, top_n=50):
    genes = [g for g in gene_list[:top_n] if g in adata.var_names]
    if genes:
        # gsmap_input .X is RAW integer counts; score on a log1p-CPM layer so the
        # per-spot library size does not confound the score (F193 / ensure_lognorm).
        ln = ensure_lognorm(adata)
        scored = adata.copy()
        scored.X = scored.layers[ln]
        sc.tl.score_genes(scored, gene_list=genes,
                          score_name="svg_score", use_raw=False)
        return scored.obs["svg_score"].values
    return np.zeros(adata.n_obs)

# ── Helper: auto spot size from coordinate range and figure panel width ───────
def auto_spot_size(coords, panel_width_in):
    x = coords[:, 0]
    x_range = x.max() - x.min()
    scale = (panel_width_in * 72) / x_range
    spacing = x_range / np.sqrt(len(x))
    radius_pts = spacing * scale * 0.52
    return np.pi * radius_pts ** 2

# ═══════════════════════════════════════════════════════════════════════════════
# Panel A/B: tissue sections with zone overlay
# ═══════════════════════════════════════════════════════════════════════════════
tissue_samples = {
    "Healthy\n(GSE192741 JBO018)":         "gse192741_JBO018.h5ad",
    "Steatohepatitis\n(GSE192741 JBO019)": "gse192741_JBO019.h5ad",
}

tissue_adatas = {}
for label, fname in tissue_samples.items():
    fpath = os.path.join(DATA_DIR, fname)
    if os.path.exists(fpath):
        adata = sc.read_h5ad(fpath)
        sample_id = fname.split("_")[1].replace(".h5ad", "")
        # Join zonation bin
        zs = zone_df[zone_df["sample_id"] == sample_id][["zonation_bin"]]
        adata.obs = adata.obs.join(zs, how="left")
        tissue_adatas[label] = adata
    else:
        print(f"WARNING: {fname} not found")

# ═══════════════════════════════════════════════════════════════════════════════
# Panel C: SVG score violin — all spots, all samples
# ═══════════════════════════════════════════════════════════════════════════════
violin_data = []

# GSE192741
gse_map = {
    "gse192741_JBO018.h5ad": "Healthy",
    "gse192741_JBO022.h5ad": "Healthy",
    "gse192741_JBO014.h5ad": "Steatohepatitis",
    "gse192741_JBO015.h5ad": "Steatohepatitis",
    "gse192741_JBO019.h5ad": "Steatohepatitis",
}
for fname, cond in gse_map.items():
    fpath = os.path.join(DATA_DIR, fname)
    if os.path.exists(fpath):
        adata = sc.read_h5ad(fpath)
        scores = compute_svg_score(adata, disease_emergent)
        violin_data.append(pd.DataFrame({
            "svg_score": scores, "condition": cond,
            "sample": fname.split(".")[0]
        }))

# Vu et al.
vu_files = [f for f in os.listdir(DATA_DIR) if f.startswith("vu_") and f.endswith(".h5ad")]
for fname in sorted(vu_files):
    fpath = os.path.join(DATA_DIR, fname)
    adata = sc.read_h5ad(fpath)
    scores = compute_svg_score(adata, disease_emergent)
    violin_data.append(pd.DataFrame({
        "svg_score": scores, "condition": "Vu (MASLD)",
        "sample": fname.split(".")[0]
    }))

vdf = pd.concat(violin_data, ignore_index=True)

# Condition order and stats.
# F193: spots are pseudoreplicates of ~5 GSE192741 donors (Healthy n=2:
# JBO018/JBO022; Steatohepatitis n=3: JBO014/JBO015/JBO019) + Vu slices. Label
# the x-axis with the true DONOR/slice n, not the spot n, and overlay the
# per-donor median scores. The violins remain spot-level DESCRIPTIVE
# distributions (no pooled-spot significance is claimed).
cond_order = ["Healthy", "Steatohepatitis", "Vu (MASLD)"]
# Per-donor (per-sample) median SVG score — the donor-level unit of inference.
donor_med = (vdf.groupby(["condition", "sample"])["svg_score"]
             .median().reset_index())
cond_ndonor = donor_med.groupby("condition")["sample"].nunique()
cond_labels = [f"{c}\n({cond_ndonor.get(c,0)} donors)" for c in cond_order]

# ═══════════════════════════════════════════════════════════════════════════════
# Figure layout: [tissue H | tissue SH | violin] 1:1:1.1
# ═══════════════════════════════════════════════════════════════════════════════
PANEL_W = 2.2   # tissue panel width (inches)
VIO_W   = 2.4
FIG_H   = 2.8
fig = plt.figure(figsize=(PANEL_W * 2 + VIO_W + 0.5, FIG_H))
gs  = gridspec.GridSpec(1, 3, width_ratios=[PANEL_W, PANEL_W, VIO_W],
                         wspace=0.08, left=0.01, right=0.98,
                         top=0.88, bottom=0.08)

# ── Tissue panels ─────────────────────────────────────────────────────────────
for idx, (label, adata) in enumerate(tissue_adatas.items()):
    ax = fig.add_subplot(gs[0, idx])
    coords = adata.obsm["spatial"]
    x, y   = coords[:, 0], coords[:, 1]

    # Tight crop
    xpad = (x.max() - x.min()) * 0.03
    ypad = (y.max() - y.min()) * 0.03
    xlim = (x.min() - xpad, x.max() + xpad)
    ylim = (y.min() - ypad, y.max() + ypad)

    s = auto_spot_size(coords, PANEL_W)

    if "zonation_bin" in adata.obs.columns:
        c_vals = adata.obs["zonation_bin"].map(ZONE_COLORS).fillna("#CCCCCC").values
    else:
        c_vals = "#CCCCCC"

    ax.scatter(x, -y, c=c_vals, s=s, linewidths=0, rasterized=True)
    ax.set_xlim(xlim)
    ax.set_ylim(-ylim[1], -ylim[0])
    ax.set_aspect("equal")
    ax.axis("off")
    title_short = label.split("\n")[0]
    ax.set_title(title_short, fontsize=7, fontweight="bold", pad=3)

# Zone legend on first tissue panel
ax0 = fig.axes[0]
legend_elements = [Patch(facecolor=ZONE_COLORS[z], edgecolor="none",
                         label=z if z in ("PP1", "PC1") else z)
                   for z in zone_bin_order]
ax0.legend(handles=legend_elements, fontsize=5, title="Zone",
           title_fontsize=5.5, loc="lower left",
           frameon=True, framealpha=0.85, edgecolor="none",
           handlelength=1, handleheight=0.8)

# ── Violin panel ──────────────────────────────────────────────────────────────
ax_v = fig.add_subplot(gs[0, 2])

positions = list(range(len(cond_order)))
for pos, cond in zip(positions, cond_order):
    data = vdf.loc[vdf["condition"] == cond, "svg_score"].values
    if len(data) == 0:
        continue
    parts = ax_v.violinplot(data, positions=[pos], widths=0.65,
                             showmedians=False, showextrema=False)
    color = COND_COLORS[cond]
    for pc in parts["bodies"]:
        pc.set_facecolor(color)
        pc.set_alpha(0.75)
        pc.set_edgecolor("none")
    # Median line
    med = np.median(data)
    ax_v.hlines(med, pos - 0.2, pos + 0.2,
                color="white", linewidth=1.2, zorder=5)
    # Overlay per-donor median scores (the donor-level unit of inference, F193).
    dvals = donor_med.loc[donor_med["condition"] == cond, "svg_score"].values
    if len(dvals):
        jitter = (np.random.RandomState(42).rand(len(dvals)) - 0.5) * 0.18
        ax_v.scatter(np.full(len(dvals), pos) + jitter, dvals,
                     s=14, color="#222222", edgecolor="white",
                     linewidth=0.4, zorder=6)

ax_v.set_xticks(positions)
ax_v.set_xticklabels(cond_labels, fontsize=5.5)
ax_v.set_ylabel("Disease-emergent SVG score", fontsize=6)
ax_v.set_title("Disease-emergent SVG score per spot\n(points = per-donor medians)",
               fontsize=7, fontweight="bold")
ax_v.tick_params(axis="y", labelsize=5.5)
ax_v.spines["top"].set_visible(False)
ax_v.spines["right"].set_visible(False)
ax_v.spines["bottom"].set_linewidth(0.5)
ax_v.spines["left"].set_linewidth(0.5)
ax_v.tick_params(width=0.5)

# Figure caption — report the DONOR n, not the pooled spot n (F193). The violin
# shows the spot-level score distribution (descriptive); black points are the
# per-donor medians (the unit of inference: Healthy n=2, Steatohepatitis n=3
# GSE192741 donors; Vu n=10 slices).
n_total = len(vdf)
fig.text(0.5, 0.01,
         "Left: spots colored by periportal (blue) → pericentral (orange) zonation zone  |  "
         "Right: disease-emergent SVG score per spot (violins, descriptive); "
         "black points = per-donor medians (Healthy n=2, Steatohepatitis n=3 donors; Vu n=10 slices)",
         ha="center", fontsize=5.5, color="#444444")

out = os.path.join(OUT_DIR, "fig4b_spatial_overview.pdf")
plt.savefig(out, bbox_inches="tight", dpi=300)
print(f"Saved: {out}")
print(f"  GSE192741: {sum(1 for c in gse_map.values() if c=='Healthy')} healthy, "
      f"{sum(1 for c in gse_map.values() if c=='Steatohepatitis')} steatohepatitis samples")
print(f"  Vu et al.: {len(vu_files)} samples")
print(f"  Total spots in violin: {n_total:,}")
