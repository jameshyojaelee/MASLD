#!/usr/bin/env python
"""
351_scvi_fstage_umap_figure.py

Publication-quality 4-panel figure showing that the scVI latent space of
hepatocytes organizes by F-stage biology and that an ordinal classifier
trained on the documented Andrews 58 generalizes to all 260 donors.

Panel A  Per-cell hepatocyte UMAP coloured by donor's *documented* F-stage
         (Andrews-only; rest gray).
Panel B  Per-cell hepatocyte UMAP coloured by *inferred* F-stage argmax
         (all 260 donors with prediction; rest gray).
Panel C  Per-cell hepatocyte UMAP coloured by inferred F-stage posterior
         confidence (max-probability per donor).
Panel D  Donor-level UMAP of donor-mean hep-scVI (269x20) coloured by
         inferred F-stage; shape = documented (filled) vs scvi_predicted (open).

Output : figures/supplementary/stage_ccc/figS_scvi_fstage_umap.pdf
Conv.  : PDF only, sans-serif fonts (Arial), 8pt labels, 10pt titles,
         rasterized scatter for the per-cell panels (300k+ points).
"""
from __future__ import annotations

import os
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D

import anndata as ad
import scanpy as sc

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

# ------------------------------------------------------------------ paths
ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
))

H5AD = ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad"
PRED = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_scvi_predicted.tsv"
DOC  = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_fstage_documented.tsv"
DMEAN = ROOT / "Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_hep_scvi_mean.tsv"

OUT_DIR = ROOT / "figures/supplementary/stage_ccc"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_PDF = OUT_DIR / "figS_scvi_fstage_umap.pdf"

# ------------------------------------------------------------------ matplotlib defaults
# Mirrors theme_masld(): minimal, no top/right spines, sans-serif, 7-8pt.
mpl.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 8,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "legend.frameon": False,
    "axes.linewidth": 0.4,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": False,
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",
    "savefig.dpi": 300,
})

# ------------------------------------------------------------------ palettes
# MASLD project palette (publication_theme.R `masld_colors`):
#   control / healthy = warm gray (project rule: control is always gray, NEVER red)
#   disease gradient  = white -> rose -> magenta -> dark magenta (F1 -> F4)
# Same fibrosis_stage_colors used in fig2 family but recolored on the MASH
# magenta axis instead of the blue axis to keep "disease = red/magenta".
STAGE_COLORS = {
    0: "#BDBDBD",  # F0 / Healthy   — masld_colors$ns (warm gray)
    1: "#F8BBD0",  # F1             — pale pink
    2: "#F48FB1",  # F2             — masld_colors$masl
    3: "#C2185B",  # F3             — masld_colors$mash
    4: "#880E4F",  # F4             — masld_colors$fibrosis (dark magenta)
}
GRAY = "#E0E0E0"   # unlabeled cells — slightly lighter than F0 so the two read distinctly

stage_cmap = ListedColormap([STAGE_COLORS[i] for i in range(5)])
stage_norm = BoundaryNorm(np.arange(-0.5, 5.5, 1.0), stage_cmap.N)

# ------------------------------------------------------------------ load
print("[load] reading donor tables", flush=True)
doc_df  = pd.read_csv(DOC,   sep="\t")
pred_df = pd.read_csv(PRED,  sep="\t")
dmean   = pd.read_csv(DMEAN, sep="\t")

# Map sample -> documented F-stage (numeric)
sample_to_doc = doc_df.set_index("sample")["F_stage_documented"].to_dict()
# Map sample -> inferred F-stage argmax
sample_to_pred = pred_df.set_index("sample")["F_stage_predicted_argmax"].to_dict()
# Map sample -> max posterior
prob_cols = ["P_F0","P_F1","P_F2","P_F3","P_F4"]
pred_df["max_p"] = pred_df[prob_cols].max(axis=1)
sample_to_conf = pred_df.set_index("sample")["max_p"].to_dict()
sample_to_dataset = pred_df.set_index("sample")["dataset"].to_dict()

n_doc = doc_df["F_stage_documented"].notna().sum()
n_pred = pred_df["F_stage_predicted_argmax"].notna().sum()
print(f"[load] documented F-stage on {n_doc} donors; predicted on {n_pred} donors", flush=True)

# Load h5ad in backed mode and pull only what we need
print("[load] opening h5ad backed='r'", flush=True)
t0 = time.time()
adata = ad.read_h5ad(H5AD, backed="r")
print(f"[load] backed h5ad: {adata.shape} in {time.time()-t0:.1f}s", flush=True)

# Build minimal per-cell dataframe (no expression matrix loaded)
obs_use = adata.obs[["sample", "dataset"]].copy()
print(f"[load] obs slice: {obs_use.shape}", flush=True)

if "X_umap" in adata.obsm:
    print("[load] reusing existing X_umap from h5ad", flush=True)
    umap = np.asarray(adata.obsm["X_umap"])
else:
    # fallback: compute fresh UMAP from X_scVI; requires loading scVI matrix
    print("[load] X_umap missing; computing fresh UMAP from X_scVI", flush=True)
    scvi = np.asarray(adata.obsm["X_scVI"])
    tmp = ad.AnnData(X=np.zeros((scvi.shape[0], 1), dtype=np.float32))
    tmp.obsm["X_scVI"] = scvi
    sc.pp.neighbors(tmp, use_rep="X_scVI", n_neighbors=30, random_state=42)
    sc.tl.umap(tmp, min_dist=0.3, random_state=42)
    umap = np.asarray(tmp.obsm["X_umap"])
    del tmp, scvi

assert umap.shape[0] == obs_use.shape[0], (umap.shape, obs_use.shape)
obs_use = obs_use.reset_index(drop=True)
obs_use["umap1"] = umap[:, 0]
obs_use["umap2"] = umap[:, 1]

# Per-cell colour vectors
obs_use["doc_stage"]  = obs_use["sample"].map(sample_to_doc)
obs_use["pred_stage"] = obs_use["sample"].map(sample_to_pred)
obs_use["pred_conf"]  = obs_use["sample"].map(sample_to_conf)

print("[load] per-cell coverage:",
      f"documented={obs_use['doc_stage'].notna().sum():,}/{len(obs_use):,}",
      f"predicted={obs_use['pred_stage'].notna().sum():,}/{len(obs_use):,}",
      flush=True)

# ------------------------------------------------------------------ donor-level UMAP
print("[donor] computing UMAP on donor-mean hep-scVI", flush=True)
zcols = [c for c in dmean.columns if c.startswith("z")]
Z = dmean[zcols].values.astype(np.float32)
tmp = ad.AnnData(X=np.zeros((Z.shape[0], 1), dtype=np.float32))
tmp.obsm["X_scVI"] = Z
n_neigh = min(15, Z.shape[0] - 1)
sc.pp.neighbors(tmp, use_rep="X_scVI", n_neighbors=n_neigh, random_state=42)
sc.tl.umap(tmp, min_dist=0.3, random_state=42)
donor_umap = np.asarray(tmp.obsm["X_umap"])
del tmp

donor_df = dmean[["sample", "dataset"]].copy().reset_index(drop=True)
donor_df["umap1"] = donor_umap[:, 0]
donor_df["umap2"] = donor_umap[:, 1]
donor_df["doc_stage"]  = donor_df["sample"].map(sample_to_doc)
donor_df["pred_stage"] = donor_df["sample"].map(sample_to_pred)
donor_df["pred_conf"]  = donor_df["sample"].map(sample_to_conf)
# source = documented if doc not NA, else scvi_predicted
donor_df["source"] = np.where(donor_df["doc_stage"].notna(), "documented", "scvi_predicted")
donor_df.loc[donor_df["pred_stage"].isna() & donor_df["doc_stage"].isna(), "source"] = "none"
print(f"[donor] {donor_df.shape[0]} donors; doc={donor_df['doc_stage'].notna().sum()} "
      f"pred={donor_df['pred_stage'].notna().sum()}", flush=True)

# ------------------------------------------------------------------ figure
print("[plot] building 2x2 figure", flush=True)
fig = plt.figure(figsize=(7.5, 7.5))
gs = fig.add_gridspec(2, 2, hspace=0.22, wspace=0.18)

# Helper -----------------------------------------------------------
def _scatter_cells(ax, df, color_col, vmin=None, vmax=None, cmap=None, norm=None,
                   gray_mask=None, point_size=0.6, alpha=0.7, title="",
                   continuous=False):
    # Plot gray cells underneath
    if gray_mask is not None and gray_mask.any():
        ax.scatter(df.loc[gray_mask, "umap1"], df.loc[gray_mask, "umap2"],
                   s=point_size, c=GRAY, alpha=0.35, linewidths=0,
                   rasterized=True)
    colored_mask = ~gray_mask if gray_mask is not None else np.ones(len(df), dtype=bool)
    # Shuffle so no single stage sits on top
    idx = np.where(colored_mask)[0]
    rng = np.random.default_rng(42)
    rng.shuffle(idx)
    sub = df.iloc[idx]
    if continuous:
        sc_obj = ax.scatter(sub["umap1"], sub["umap2"],
                            s=point_size, c=sub[color_col], cmap=cmap,
                            vmin=vmin, vmax=vmax, alpha=alpha, linewidths=0,
                            rasterized=True)
    else:
        sc_obj = ax.scatter(sub["umap1"], sub["umap2"],
                            s=point_size, c=sub[color_col].astype(float),
                            cmap=cmap, norm=norm, alpha=alpha, linewidths=0,
                            rasterized=True)
    ax.set_title(title, loc="left", fontweight="bold", fontsize=8.5, pad=2)
    ax.set_xlabel("UMAP 1", labelpad=1)
    ax.set_ylabel("UMAP 2", labelpad=1)
    ax.set_xticks([]); ax.set_yticks([])
    for spine in ("left", "bottom"):
        ax.spines[spine].set_linewidth(0.4)
    return sc_obj

def _stage_legend(ax, loc="lower right"):
    handles = [Line2D([0], [0], marker="o", color="w",
                      markerfacecolor=STAGE_COLORS[k], markersize=4.5,
                      markeredgecolor="none",
                      label=f"F{k}") for k in range(5)]
    handles.append(Line2D([0], [0], marker="o", color="w",
                          markerfacecolor=GRAY, markersize=4.5,
                          markeredgecolor="none",
                          label="unlabelled"))
    ax.legend(handles=handles, loc=loc, frameon=False, ncol=1,
              handletextpad=0.25, borderpad=0.1, labelspacing=0.18,
              fontsize=6)

# Panel A -----------------------------------------------------------
ax = fig.add_subplot(gs[0, 0])
gray_a = obs_use["doc_stage"].isna().values
_scatter_cells(ax, obs_use, "doc_stage", cmap=stage_cmap, norm=stage_norm,
               gray_mask=gray_a,
               title="A  Documented F-stage (Andrews subset)")
_stage_legend(ax)

# Panel B -----------------------------------------------------------
ax = fig.add_subplot(gs[0, 1])
gray_b = obs_use["pred_stage"].isna().values
_scatter_cells(ax, obs_use, "pred_stage", cmap=stage_cmap, norm=stage_norm,
               gray_mask=gray_b,
               title="B  Inferred F-stage (scVI classifier, all donors)")
_stage_legend(ax)

# Panel C -----------------------------------------------------------
# Posterior confidence is a magnitude (0.2-1.0). Use a perceptually-uniform
# sequential colormap (matplotlib built-in viridis) instead of RdBu_r which
# would imply diverging-from-zero semantics. Project rule: do not use red
# for control / magnitude scales unless centered at zero. (seaborn's 'mako'
# is not registered in this env, so we stick to the matplotlib built-in.)
ax = fig.add_subplot(gs[1, 0])
gray_c = obs_use["pred_conf"].isna().values
sc_c = _scatter_cells(ax, obs_use, "pred_conf", cmap="viridis",
                      vmin=0.2, vmax=1.0, continuous=True,
                      gray_mask=gray_c,
                      title="C  Inferred posterior confidence (max P)")
cbar = fig.colorbar(sc_c, ax=ax, fraction=0.035, pad=0.015, shrink=0.55)
cbar.set_label("max P(F-stage)", fontsize=6)
cbar.ax.tick_params(labelsize=5.5, length=2, width=0.4)
cbar.outline.set_linewidth(0.3)

# Panel D -----------------------------------------------------------
ax = fig.add_subplot(gs[1, 1])
# choose colour by inferred argmax (fall back to documented)
donor_df["stage_for_color"] = donor_df["pred_stage"].fillna(donor_df["doc_stage"])
# split by source for marker shape
for src, marker, fc_mode in [
    ("documented",     "o", "filled"),
    ("scvi_predicted", "o", "open"),
]:
    sub = donor_df[donor_df["source"] == src]
    if sub.empty:
        continue
    colors = [STAGE_COLORS.get(int(s), GRAY) if pd.notna(s) else GRAY
              for s in sub["stage_for_color"]]
    if fc_mode == "filled":
        ax.scatter(sub["umap1"], sub["umap2"], s=28, c=colors,
                   marker=marker, edgecolors="black", linewidths=0.3,
                   alpha=0.95, rasterized=True)
    else:
        ax.scatter(sub["umap1"], sub["umap2"], s=28, facecolors="none",
                   edgecolors=colors, marker=marker, linewidths=0.8,
                   alpha=0.95, rasterized=True)
# orphans
orph = donor_df[donor_df["source"] == "none"]
if not orph.empty:
    ax.scatter(orph["umap1"], orph["umap2"], s=20, c=GRAY,
               marker="x", linewidths=0.6, alpha=0.7)
ax.set_title("D  Donor-mean hep-scVI UMAP", loc="left", fontweight="bold",
             fontsize=8.5, pad=2)
ax.set_xlabel("UMAP 1", labelpad=1)
ax.set_ylabel("UMAP 2", labelpad=1)
ax.set_xticks([]); ax.set_yticks([])
for spine in ("left", "bottom"):
    ax.spines[spine].set_linewidth(0.4)

# composite legend for panel D
stage_handles = [Line2D([0], [0], marker="o", color="w",
                        markerfacecolor=STAGE_COLORS[k], markersize=4.5,
                        markeredgecolor="black", markeredgewidth=0.25,
                        label=f"F{k}") for k in range(5)]
src_handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor="dimgray",
           markeredgecolor="black", markeredgewidth=0.25, markersize=4.5,
           label="documented"),
    Line2D([0], [0], marker="o", color="w", markerfacecolor="none",
           markeredgecolor="dimgray", markersize=4.5, markeredgewidth=0.8,
           label="scVI-inferred"),
]
leg1 = ax.legend(handles=stage_handles, loc="upper left", frameon=False,
                 ncol=1, handletextpad=0.25, borderpad=0.1, labelspacing=0.18,
                 fontsize=6, title="F-stage", title_fontsize=6)
leg1.get_title().set_fontweight("bold")
ax.add_artist(leg1)
leg2 = ax.legend(handles=src_handles, loc="lower right", frameon=False,
                 handletextpad=0.3, borderpad=0.1, labelspacing=0.18,
                 fontsize=6, title="source", title_fontsize=6)
leg2.get_title().set_fontweight("bold")

# ------------------------------------------------------------------ save
fig.suptitle("scVI latent space organizes by F-stage biology",
             fontsize=9.5, fontweight="bold", y=0.995, x=0.06, ha="left")
print(f"[save] writing {OUT_PDF}", flush=True)
fig.savefig(OUT_PDF, bbox_inches="tight", dpi=400)
plt.close(fig)

# ------------------------------------------------------------------ summary stats for caller report
print("[stats] computing per-stage composition for report", flush=True)
comp_doc = (obs_use[obs_use["doc_stage"].notna()]
            .groupby(["doc_stage", "dataset"]).size()
            .unstack(fill_value=0))
comp_pred = (obs_use[obs_use["pred_stage"].notna()]
             .groupby(["pred_stage", "dataset"]).size()
             .unstack(fill_value=0))
conf_by_stage = (obs_use[obs_use["pred_stage"].notna()]
                 .groupby("pred_stage")["pred_conf"]
                 .agg(["mean", "median", "min", "max"]))
print("[stats] documented composition by dataset (cells)")
print(comp_doc.to_string())
print("[stats] inferred composition by dataset (cells)")
print(comp_pred.to_string())
print("[stats] confidence by inferred stage")
print(conf_by_stage.to_string())
print(f"[stats] mean per-cell max-P across all predicted cells: "
      f"{obs_use['pred_conf'].dropna().mean():.3f}")
fsize = OUT_PDF.stat().st_size / 1e6
print(f"[done] {OUT_PDF}  ({fsize:.2f} MB)", flush=True)
