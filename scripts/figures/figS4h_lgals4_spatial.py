#!/usr/bin/env python
"""
figS4h — LGALS4 spatial expression over registered H&E (GSE192741 human Visium).

Replaces the retired CYP3A4 spatial panel. Rationale: CYP3A4 is DOWN in disease
(bulk mRNA -0.42, protein -0.69) yet its per-spot log1p-CP10K *rises* in steatotic
sections purely from the compositional effect of lipid-engorged hepatocytes — the
image contradicted the quantification. LGALS4 is the opposite: a certified disease-UP
DEG whose spatial image honestly AGREES with its direction.

LGALS4 (human, no mouse used):
  bulk mRNA   logFC +0.97  padj 7.6e-12   (passes effect-size interval-null gate, treat_fdr 1.5e-6)
  protein     logFC +1.64  padj 0.050
  scRNA hep   logFC +0.81  (direction concordant)
  genetic     COLOC PP.H4 = 0.92
  Visium      53% spots detected; steatotic mean 0.73 vs healthy 0.27 (+0.46), all 5 sections agree.

The disease-UP claim rests on bulk+protein+COLOC; the spatial map is an illustrative
confirmation (single cohort), not the primary evidence.

Convention: 6 pt text, black, no bold, no colored/stat-laden title, gene symbol italic,
caption emitted to stdout (PI directive). PDF only.

Output: figures/main/fig5_molecular_context/panels/figS4h.pdf
"""
import os, json
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.environ.get("MASLD_PROJECT_ROOT",
                      "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SR   = os.path.join(BASE, "Analysis/Spatial/results/spaceranger/GSE192741")
OUT  = os.path.join(BASE, "figures/main/fig5_molecular_context/panels/figS4h.pdf")
GENE = "LGALS4"
# All 5 condition-labeled GSE192741 sections (2 Healthy, 3 MASLD), ordered low->high
# expression within each row. Showing every section (nothing dropped) makes the
# monotonic increase honest: every MASLD section exceeds every Healthy section
# (per-section mean log1p: 0.11, 0.43 healthy | 0.55, 0.62, 1.02 MASLD).
HEALTHY = ["JBO018", "JBO022"]
MASLD   = ["JBO014", "JBO015", "JBO019"]

plt.rcParams.update({
    "font.size": 6, "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.5, "savefig.dpi": 300,
})

def load(sample):
    outs = os.path.join(SR, sample, "outs"); spdir = os.path.join(outs, "spatial")
    ad = sc.read_10x_h5(os.path.join(outs, "filtered_feature_bc_matrix.h5"))
    ad.var_names_make_unique()
    pos_path = next(os.path.join(spdir, fn) for fn in
                    ("tissue_positions_list.csv", "tissue_positions.csv")
                    if os.path.exists(os.path.join(spdir, fn)))
    hdr = 0 if pos_path.endswith("tissue_positions.csv") else None
    pos = pd.read_csv(pos_path, header=hdr)
    if hdr is None:
        pos.columns = ["barcode", "in_tissue", "array_row", "array_col",
                       "pxl_row_in_fullres", "pxl_col_in_fullres"]
    pos = pos.set_index("barcode").reindex(ad.obs_names)
    sf = json.load(open(os.path.join(spdir, "scalefactors_json.json")))
    img = plt.imread(os.path.join(spdir, "tissue_hires_image.png"))
    scalef = sf["tissue_hires_scalef"]
    sc.pp.normalize_total(ad, target_sum=1e4); sc.pp.log1p(ad)
    x = pos["pxl_col_in_fullres"].values * scalef
    y = pos["pxl_row_in_fullres"].values * scalef
    v = ad[:, GENE].X
    v = np.asarray(v.todense()).ravel() if hasattr(v, "todense") else np.asarray(v).ravel()
    return dict(sample=sample, img=img, x=x, y=y, v=v)

healthy = [load(s) for s in HEALTHY]
mash    = [load(s) for s in MASLD]
allv = np.concatenate([d["v"] for d in healthy + mash])
vmax = float(np.percentile(allv, 99))

def draw(ax, d):
    ax.imshow(d["img"], origin="upper")
    order = np.argsort(d["v"])                          # draw high-expr spots on top
    sca = ax.scatter(d["x"][order], d["y"][order], c=d["v"][order],
                     cmap="magma", vmin=0, vmax=vmax, s=2.6, linewidths=0,
                     alpha=0.9, rasterized=True)
    x0, x1 = np.nanmin(d["x"]), np.nanmax(d["x"])
    y0, y1 = np.nanmin(d["y"]), np.nanmax(d["y"])
    dx, dy = (x1 - x0) * 0.04, (y1 - y0) * 0.04
    ax.set_xlim(x0 - dx, x1 + dx); ax.set_ylim(y1 + dy, y0 - dy)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.set_title(d["sample"], fontsize=6, color="black", pad=1.5)
    return sca

fig = plt.figure(figsize=(4.9, 3.4))
gs = fig.add_gridspec(2, 3)
sca = None
for j, d in enumerate(healthy):                          # row 0: 2 Healthy (cols 0,1)
    sca = draw(fig.add_subplot(gs[0, j]), d)
    if j == 0:
        fig.axes[-1].set_ylabel("Healthy", fontsize=6, color="black", rotation=90, labelpad=2)
for j, d in enumerate(mash):                              # row 1: 3 MASLD
    ax = fig.add_subplot(gs[1, j]); sca = draw(ax, d)
    if j == 0:
        ax.set_ylabel("MASLD", fontsize=6, color="black", rotation=90, labelpad=2)

fig.subplots_adjust(left=0.05, right=0.88, top=0.93, bottom=0.02, wspace=0.06, hspace=0.14)
cax = fig.add_axes([0.905, 0.20, 0.02, 0.58])
cb = fig.colorbar(sca, cax=cax)
cb.ax.tick_params(labelsize=6, width=0.4, length=2)
cb.outline.set_linewidth(0.3)
cb.set_label(f"$\\it{{{GENE}}}$  log1p", fontsize=6, color="black")

fig.savefig(OUT, bbox_inches="tight")
print(f"saved {OUT}")

# ---- caption / stats to stdout (figure legend; PI directive: not on-panel) ----
hmean = np.mean([d["v"].mean() for d in healthy])
smean = np.mean([d["v"].mean() for d in mash])
print(f"[caption] {GENE} expression (log1p normalized counts) over registered H&E, "
      f"GSE192741 human Visium. Top: two Healthy sections; bottom: three MASLD (steatotic) "
      f"sections; shared color scale. {GENE} is UP in MASLD across all sections "
      f"(mean log1p {hmean:.2f} healthy vs {smean:.2f} MASLD). The spatial map illustrates a "
      f"disease-up gene whose direction is established by bulk mRNA (logFC +0.97, padj 7.6e-12), "
      f"proteomics (logFC +1.64, padj 0.050) and genetic colocalization (COLOC PP.H4 = 0.92); "
      f"single-cohort spatial data is confirmatory, not the primary evidence.")
