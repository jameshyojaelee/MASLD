#!/usr/bin/env python
# ============================================================================
# hep_hotspot_composite.py
#
# "Parallel Hotspot modules" left-panel analog, built for HEPATOCYTES.
# The reference figure is a CAR-T perturbation screen; our condition axis is
# DISEASE STAGE (Healthy -> Steatosis -> Steatohepatitis -> Cirrhosis), and the
# up/down-in-disease modules are chosen by sign of disease_stage_beta:
#   UP   = hep__20 Ductular injury (BICC1) (beta = +0.67, q = 6e-20) — CANONICAL
#          up-hero, consistent with fig3h/fig3i (changed from hep__17 2026-07-01
#          so all Fig3 hep-module panels feature the same 19/20 pair)
#   DOWN = hep__19 Fatty-acid/peroxisomal (FAO) (beta = -0.63, q = 1e-8)
#
# Individual PDFs (assembled in Illustrator) ->
#   figures/main/fig3_RNAseq/panels/figs3_hep_hotspot_composite/
#     A_local_corr_triangle.pdf   gene x gene local-autocorr, lower triangle
#     B_umap_stage.pdf            hep UMAP colored by disease stage
#     C_umap_up_hep20.pdf         hep UMAP colored by hep__20 score (up-hero)
#     D_umap_down_hep19.pdf       hep UMAP colored by hep__19 score
#     F_umap_up_hep26.pdf         hep UMAP colored by hep__26 (AP-1 stress) score
#     G_umap_up_hep24.pdf         hep UMAP colored by hep__24 (NRF2 antioxidant) score
#     H_umap_down_hep27.pdf       hep UMAP colored by hep__27 (Complement/CFH) score
#     E_module_stage_heatmap.pdf  stage-significant modules x stage (dScore)
#
# All data pre-exists (no new analysis): hotspot pkl local_corr, per-cell
# module scores (cell_scores.parquet), hep-specific X_umap + stage labels
# (hepatocyte_atlas_annotated.h5ad), and all_modules.tsv annotations.
# ============================================================================
import os
import pickle
from pathlib import Path

import numpy as np
import pandas as pd
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
mpl.rcParams.update({
    'font.size': 6, 'font.family': 'Helvetica', 'font.weight': 'normal',
    'axes.titlesize': 6, 'axes.titleweight': 'normal', 'axes.labelsize': 6,
    'axes.labelweight': 'normal', 'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'legend.fontsize': 6, 'legend.title_fontsize': 6, 'figure.titlesize': 6,
    'figure.titleweight': 'normal', 'pdf.fonttype': 42, 'ps.fonttype': 42,
    'savefig.dpi': 600, 'figure.dpi': 600,   # high-res rasterized scatters/heatmaps (2026-07-01)
})
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy.cluster.hierarchy import linkage, leaves_list

ROOT = Path(os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"))
HS = ROOT / "Analysis/SingleCell/results_gpu_v2/hotspot_modules"
HEPA = ROOT / "Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_atlas_annotated.h5ad"
OUT = ROOT / "figures/main/fig3_RNAseq/panels/figs3_hep_hotspot_composite"
DATA = OUT / "data"
for d in (OUT, DATA):
    d.mkdir(parents=True, exist_ok=True)

CT = "hepatocytes"
UP_MOD = int(os.environ.get("HEP_UP_MOD", 20))    # hep-20 Ductular injury (BICC1) — canonical up-hero (was 17; 2026-07-01, matches fig3h/fig3i)
DOWN_MOD = int(os.environ.get("HEP_DOWN_MOD", 19))

# 5 manuscript hepatocyte modules featured across Fig 3 (fig3h/fig3i) + this
# composite: direction (up/down in disease) + labels. Panel A brackets ALL 5 on
# the correlation structure; the UMAP loop (C/D/F/G/H) shows each per-cell.
# Fallback labels for the 5 manuscript hepatocyte modules; these are OVERRIDDEN
# below from the canonical module_names.tsv (`name`) once it is loaded, so they
# apply only if that file is unavailable. Kept canonical-consistent regardless.
MOD_NAME   = {19: "Fatty-acid / peroxisomal metab", 20: "Ductular injury (BICC1)",
              24: "NRF2 antioxidant (TXNRD1)", 26: "AP-1 stress (ATF3/GDF15)",
              27: "Complement (CFH)"}
MOD_SHORT  = {19: "Fatty-acid / peroxisomal metab", 20: "Ductular injury",
              24: "NRF2 antioxidant", 26: "AP-1 stress", 27: "Complement"}
HERO_ORDER = [UP_MOD, 26, 24, DOWN_MOD, 27]        # ups then downs (panel-A stack order)
HERO_DIR   = {UP_MOD: "up", 26: "up", 24: "up", DOWN_MOD: "down", 27: "down"}

GRAY = "#9E9E9E"
RED = "#C9265E"
BLUE = "#1565C0"
STAGE_ORDER = ["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis"]
STAGE_COL = {"Healthy": GRAY, "Steatosis": "#F4A582",
             "Steatohepatitis": "#C9265E", "Cirrhosis": "#6E0B36"}
DIV = mcolors.LinearSegmentedColormap.from_list("bwr_liang", [BLUE, "#FFFFFF", RED])

# ---------------------------------------------------------------------------
# Load shared inputs
# ---------------------------------------------------------------------------
print("[load] hotspot pkl")
d = pickle.load(open(HS / CT / "hotspot_obj.pkl", "rb"))
modules = d["modules"]
lc = d["local_corr"]

print("[load] all_modules.tsv (hep annotations)")
am = pd.read_csv(HS / "all_modules.tsv", sep="\t")
am = am[am["cell_type"] == CT].copy()
am["module"] = am["module"].astype(int)
beta = dict(zip(am["module"], am["disease_stage_beta"]))
qval = dict(zip(am["module"], am["disease_stage_q"]))
# functional names from 509_module_pathway_names (persistent authority, survives
# a 508 rebuild); fall back to best_match_program if the annotation isn't present.
prog = dict(zip(am["module"], am["best_match_program"].astype(str)))
try:
    mn = pd.read_csv(HS / "module_names.tsv", sep="\t")
    mn = mn[mn["cell_type"] == CT]
    name = dict(zip(mn["module"].astype(int), mn["module_name"].astype(str)))
except Exception:
    name = {}

# Align the 5 manuscript-module labels to the canonical Hotspot module_name
# values so Panel A brackets + the UMAP titles match fig3h/fig3i (and don't
# contradict the other Fig 3 panels). Full canonical name -> UMAP titles; the
# same name with its trailing gene parenthetical stripped -> Panel-A brackets.
for _m in list(MOD_NAME):
    _canon = name.get(_m)
    if _canon:
        MOD_NAME[_m]  = _canon
        MOD_SHORT[_m] = _canon.split(" (")[0].strip()
# hep__19's canonical name has no gene parenthetical to strip, so its derived
# Panel-A bracket label stays long and crowds the colorbar; use a compact,
# canonical-consistent abbreviation for that bracket only. The C/D UMAP title
# keeps the full canonical "Fatty-acid / peroxisomal metab" (it has room).
MOD_SHORT[19] = "Fatty-acid oxidation"

print("[load] hep UMAP + stage from hep-subtype atlas")
with h5py.File(HEPA, "r") as f:
    ik = f["obs"].attrs.get("_index", "_index")
    cell_ids = np.array([x.decode() if isinstance(x, bytes) else x
                         for x in f["obs"][ik][:]])
    umap = f["obsm"]["X_umap"][:, :2].astype(np.float32)
    sc = f["obs"]["disease_stage_coarse"]
    stage_cats = [x.decode() if isinstance(x, bytes) else x
                  for x in sc["categories"][:]]
    stage_codes = sc["codes"][:]
stage = np.array([stage_cats[c] if c >= 0 else "NA" for c in stage_codes])
umap_df = pd.DataFrame({"cell_id": cell_ids, "u1": umap[:, 0], "u2": umap[:, 1],
                        "stage": stage})
print(f"       {len(umap_df):,} hep cells; stages={sorted(set(stage))}")

print("[load] per-cell module scores")
cs = pd.read_parquet(HS / CT / "cell_scores.parquet")
cs["module"] = cs["module"].astype(int)

# ---------------------------------------------------------------------------
# Panel A — lower-triangle gene x gene local-autocorrelation heatmap
# ---------------------------------------------------------------------------
print("[A] triangle heatmap")
keep = modules > 0
mods_k = modules[keep]
lc_k = lc.loc[keep, keep]
order = []
for m in sorted(mods_k.unique()):
    gs = mods_k[mods_k == m].index.tolist()
    if len(gs) > 1:
        sub = lc_k.loc[gs, gs].values
        sim = sub.copy()
        np.fill_diagonal(sim, sim.max())
        dist = sim.max() - sub
        np.fill_diagonal(dist, 0)
        dist = (dist + dist.T) / 2
        dist[dist < 0] = 0
        Z = linkage(dist[np.triu_indices_from(dist, 1)], method="average")
        order.extend([gs[i] for i in leaves_list(Z)])
    else:
        order.extend(gs)
lc_ord = lc_k.loc[order, order].values
mods_ord = mods_k.loc[order].values
tri = lc_ord.astype(float).copy()
tri[np.triu_indices_from(tri, k=1)] = np.nan          # keep lower triangle
vmax = max(float(np.nanpercentile(np.abs(lc_ord[lc_ord != 0]), 98)), 4.0)

fig, ax = plt.subplots(figsize=(3.7, 3.7))
im = ax.imshow(tri, cmap=DIV,
               norm=mcolors.TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax),
               aspect="equal", interpolation="nearest", rasterized=True)
uniq = sorted(mods_k.unique())
cum = 0
bounds = {}
for m in uniq:
    n = int((mods_ord == m).sum())
    bounds[m] = (cum, cum + n)
    cum += n
for m in uniq[:-1]:
    s = bounds[m][1]
    ax.axhline(s - .5, color="0.35", lw=.12)
    ax.axvline(s - .5, color="0.35", lw=.12)
# mark the up/down module blocks: bracket the block corner + label anchored in
# the roomy empty upper-right triangle (stays clear of the colorbar).
N = len(mods_ord)
# bracket ALL 5 manuscript module blocks; labels stacked in the empty upper-right
# triangle (ups red, downs blue) — consistent with the 5 per-cell UMAP panels.
for i, mod_id in enumerate(HERO_ORDER):
    if mod_id not in bounds:
        continue
    col   = RED if HERO_DIR[mod_id] == "up" else BLUE
    arrow = "↑" if HERO_DIR[mod_id] == "up" else "↓"
    a, b  = bounds[mod_id]
    cmid  = (a + b) / 2.0
    ly    = N * (0.05 + 0.092 * i)
    ax.plot([a - .5, b - .5, b - .5], [a - .5, a - .5, b - .5],
            color=col, lw=1.0, clip_on=False)
    ax.annotate(f"{arrow} Hep-{mod_id} · {MOD_SHORT.get(mod_id, '')}", xy=(cmid, cmid),
                xytext=(N * 0.58, ly),
                ha="left", va="center", fontsize=6, fontweight="normal",
                color=col, annotation_clip=False,
                arrowprops=dict(arrowstyle="-", color=col, lw=0.45))
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_title("Hepatocytes", fontsize=6, fontweight="normal", pad=6)
cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, shrink=0.6)
cb.set_label("Local correlation Z", fontsize=6)
cb.ax.tick_params(labelsize=6, length=2, width=.3)
cb.outline.set_linewidth(.3)
fig.savefig(OUT / "A_local_corr_triangle.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# UMAP helpers
# ---------------------------------------------------------------------------
def score_on_umap(mod):
    s = cs[cs["module"] == mod][["cell_id", "score"]]
    m = umap_df.merge(s, on="cell_id", how="left")
    return m

def umap_scatter(ax, x, y, c, cmap=None, norm=None, categorical=None,
                 order_by=None):
    if categorical is not None:
        for lev in STAGE_ORDER:            # Healthy first (background)
            sel = categorical == lev
            if sel.sum() == 0:
                continue
            ax.scatter(x[sel], y[sel], s=1.2, c=STAGE_COL[lev], linewidths=0,
                       rasterized=True, label=lev)
    else:
        idx = np.argsort(np.nan_to_num(order_by, nan=-1e9))
        sc = ax.scatter(x[idx], y[idx], s=1.2, c=c[idx], cmap=cmap, norm=norm,
                        linewidths=0, rasterized=True)
        return sc
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

# ---------------------------------------------------------------------------
# Panel B — UMAP by disease stage
# ---------------------------------------------------------------------------
print("[B] UMAP by stage")
fig, ax = plt.subplots(figsize=(2.7, 2.7))
umap_scatter(ax, umap_df["u1"].values, umap_df["u2"].values, None,
             categorical=umap_df["stage"].values)
ax.set_xticks([]); ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)
ax.set_title("Disease stage", fontsize=6, fontweight="normal", pad=4)
leg = ax.legend(loc="upper right", frameon=False, fontsize=6,
                markerscale=3, handletextpad=0.2, borderpad=0.1)
fig.savefig(OUT / "B_umap_stage.pdf", bbox_inches="tight")
plt.close(fig)

# ---------------------------------------------------------------------------
# Panels C / D / F / G / H — per-cell UMAP by module score for the 5 manuscript
# hepatocyte modules featured in the main Fig 3 (fig3h/fig3i). C/D = up/down
# heroes; F/G/H = the other three fig3i cascade modules.
# ---------------------------------------------------------------------------
UMAP_PANELS = [("C", "up",   UP_MOD),   ("D", "down", DOWN_MOD),
               ("F", "up",   26),       ("G", "up",   24),
               ("H", "down", 27)]
for panel, tag, mod in UMAP_PANELS:
    arrow = "↑" if tag == "up" else "↓"
    fn    = f"{panel}_umap_{tag}_hep{mod}.pdf"
    title = f"{arrow} Hep-{mod} · {MOD_NAME.get(mod, '')}  (β={beta.get(mod):+.2f})"
    print(f"[{panel}/{tag}] UMAP by hep__{mod} score")
    m = score_on_umap(mod)
    sv = m["score"].values.astype(float)
    lim = float(np.nanpercentile(np.abs(sv[~np.isnan(sv)]), 98))
    fig, ax = plt.subplots(figsize=(2.7, 2.7))
    sc = umap_scatter(ax, m["u1"].values, m["u2"].values, sv, cmap=DIV,
                      norm=mcolors.TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim),
                      order_by=np.abs(sv))
    ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(title, fontsize=6, fontweight="normal", pad=4)
    cb = fig.colorbar(sc, ax=ax, fraction=0.04, pad=0.02, shrink=0.55)
    cb.set_label("Module score", fontsize=6)
    cb.ax.tick_params(labelsize=6, length=2, width=.3)
    cb.outline.set_linewidth(.3)
    fig.savefig(OUT / fn, bbox_inches="tight")
    plt.close(fig)

# ---------------------------------------------------------------------------
# Panel E — stage-significant modules x stage (mean score, Healthy-centered)
# ---------------------------------------------------------------------------
print("[E] module x stage heatmap")
cid2stage = dict(zip(umap_df["cell_id"], umap_df["stage"]))
cs["stage"] = cs["cell_id"].map(cid2stage)
agg = (cs[cs["stage"].isin(STAGE_ORDER)]
       .groupby(["module", "stage"])["score"].mean().unstack("stage"))
agg = agg.reindex(columns=STAGE_ORDER)
agg_c = agg.sub(agg["Healthy"], axis=0)             # Healthy-centered dScore
# stage-significant hep modules, ordered by disease beta (up on top).
# Cirrhosis is excluded from the quantitative panel: the 19 cirrhotic donors
# are GSE202379 snRNA survivor-selected hepatocytes (fibrotic replacement
# depletes hepatocytes), matching the 3-stage convention used elsewhere
# (e.g. figS_hotspot_soft_novel). The descriptive UMAPs (B–D) still show all cells.
STAGE_VIZ = [c for c in STAGE_ORDER if c not in ("Healthy", "Cirrhosis")]
sig = sorted([m for m in am["module"] if qval.get(m, 1) < 0.05],
             key=lambda m: -beta.get(m, 0))
agg_c = agg_c.loc[[m for m in sig if m in agg_c.index], STAGE_VIZ]
def _short(s, n=30):
    return s if len(s) <= n else s[:n - 1] + "…"
labels = [f"Hep-{m} · {_short(name.get(m) or prog.get(m, 'Unresolved'))}"
          for m in agg_c.index]
vmax_e = float(np.nanpercentile(np.abs(agg_c.values), 98)) or 0.2

fig, ax = plt.subplots(figsize=(2.5, max(2.2, 0.135 * len(agg_c))))
im = ax.imshow(agg_c.values, cmap=DIV,
               norm=mcolors.TwoSlopeNorm(vmin=-vmax_e, vcenter=0, vmax=vmax_e),
               aspect="auto")
ax.set_xticks(range(agg_c.shape[1]))
ax.set_xticklabels(agg_c.columns, rotation=30, ha="right", fontsize=6)
ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels, fontsize=6)
for m_i, m in enumerate(agg_c.index):       # flag up/down heroes (keep normal weight)
    if m in (UP_MOD, DOWN_MOD):
        ax.get_yticklabels()[m_i].set_fontweight("normal")
ax.tick_params(length=2, width=.3)
for sp in ax.spines.values():
    sp.set_linewidth(.3)
ax.set_title("Module score vs Healthy", fontsize=6, fontweight="normal", pad=4)
cb = fig.colorbar(im, ax=ax, fraction=0.05, pad=0.02, shrink=0.5)
cb.set_label("Δ mean module score", fontsize=6)
cb.ax.tick_params(labelsize=6, length=2, width=.3)
cb.outline.set_linewidth(.3)
fig.savefig(OUT / "E_module_stage_heatmap.pdf", bbox_inches="tight")
plt.close(fig)

# data exports
agg_c.assign(module_label=labels).to_csv(DATA / "module_stage_dscore.csv")
am[["module", "disease_stage_beta", "disease_stage_q", "best_match_program",
    "is_novel", "progression_module"]].to_csv(DATA / "hep_module_annot.csv", index=False)
print("\n[CAPTION] Hepatocyte Hotspot module landscape. "
      "(A) Lower-triangle gene x gene local-autocorrelation Z for the "
      f"{len(umap_df):,}-cell hepatocyte compartment; diagonal blocks = modules; "
      "the 5 manuscript hepatocyte modules bracketed (ups red / downs blue). "
      "(B) hepatocyte-specific scVI UMAP colored by disease stage. "
      f"(C,D,F,G,H) same UMAP colored by per-cell score for the 5 manuscript "
      f"hepatocyte modules of fig3i: Hep-{UP_MOD} Ductular injury (up) / "
      f"Hep-{DOWN_MOD} Fatty-acid / peroxisomal (down) / Hep-26 AP-1 stress / "
      "Hep-24 NRF2 antioxidant / Hep-27 Complement. "
      "(E) stage-significant modules (disease q<0.05) x stage, "
      "mean module score minus Healthy; Cirrhosis excluded (GSE202379 snRNA "
      "survivor-selected). Modules labeled by best-match reference program.")
print("done — 8 panels written to", OUT)
