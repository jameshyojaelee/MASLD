#!/usr/bin/env python
"""
figS4h — LGALS4 combined spatial + cross-assay panel (merges former figS4h images + figS4i bars).

One panel, all human:
  left  2x2 LGALS4 spatial expression over registered H&E, GSE192741 Visium — ONE section per
        patient (4 patients: 2 Healthy, 2 MASLD). GSE192741 has only these 4 human patients;
        H35 contributed two sections (JBO014+JBO015) — JBO015 dropped as a redundant same-patient
        replicate, JBO014 kept (larger, 1248 spots).
  right cross-assay concordance: LGALS4 log2FC (MASLD vs control) across bulk RNA-seq, proteomics,
        scRNA hepatocyte, GeoMx spatial (pooled/canonical). dark = padj<0.05 (*).

Direction is established by bulk RNA-seq (+0.97, padj 7.6e-12, 846 samples) + proteomics; the
spatial images illustrate it. No COLOC (atlas 0.92 is an off-trait platelet GWAS). No Vu 2025
(MASLD-spectrum only, no healthy controls / no per-biopsy labels).

Convention: 6 pt text, black, gene symbol italic, caption to stdout. control=#9E9E9E, up/sig=#E65100.

Output: figures/main/fig5_molecular_context/panels/figS4h.pdf
"""
import os, json
import numpy as np, pandas as pd, scanpy as sc
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE  = os.environ.get("MASLD_PROJECT_ROOT",
                       "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
SR    = os.path.join(BASE, "Analysis/Spatial/results/spaceranger/GSE192741")
ATLAS = os.path.join(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
OUT   = os.path.join(BASE, "figures/main/fig5_molecular_context/panels/figS4h.pdf")
GENE  = "LGALS4"
# one section per patient: H36, H38 (Healthy); H35->JBO014, H37 (MASLD)
IMAGES = [[("JBO018", "Healthy"), ("JBO022", "Healthy")],
          [("JBO014", "MASLD"),   ("JBO019", "MASLD")]]
GRAY, DARK, LIGHT = "#9E9E9E", "#E65100", "#FFB74D"

plt.rcParams.update({
    "font.size": 6, "font.family": "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.5, "savefig.dpi": 300,
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

data = [[load(s) for s, _ in row] for row in IMAGES]
vmax = float(np.percentile(np.concatenate([d["v"] for row in data for d in row]), 99))

# atlas cross-assay values
a = pd.read_csv(ATLAS, low_memory=False)
r = a[a.human_symbol == GENE].iloc[0]
def num(c): return float(r[c])
assays = [
    ("Bulk\nRNA-seq",    num("bulk_logFC"),                              num("bulk_padj")),
    ("Proteo-\nmics",    num("best_protein_logFC"),                      num("best_protein_padj")),
    ("scRNA\nhep.",      num("sc_hepatocyte_logFC"),                     num("sc_hepatocyte_padj")),
    ("Spatial\n(GeoMx)", num("spatial_govaere2026_geomx_sh_vs_pt_logfc"),num("spatial_govaere2026_geomx_sh_vs_pt_padj")),
]

# ---------------------------------------------------------------- figure
fig = plt.figure(figsize=(4.9, 2.35))
# col 3 is an empty spacer so the colorbar's right-side label doesn't collide
# with the bar panel's y-axis label.
gs = fig.add_gridspec(2, 5, width_ratios=[1, 1, 0.08, 0.48, 1.05], wspace=0.05, hspace=0.10)

def draw_img(ax, d, ylabel=None):
    ax.imshow(d["img"], origin="upper")
    order = np.argsort(d["v"])
    sca = ax.scatter(d["x"][order], d["y"][order], c=d["v"][order], cmap="magma",
                     vmin=0, vmax=vmax, s=2.4, linewidths=0, alpha=0.9, rasterized=True)
    x0, x1 = np.nanmin(d["x"]), np.nanmax(d["x"]); y0, y1 = np.nanmin(d["y"]), np.nanmax(d["y"])
    dx, dy = (x1 - x0) * 0.04, (y1 - y0) * 0.04
    ax.set_xlim(x0 - dx, x1 + dx); ax.set_ylim(y1 + dy, y0 - dy)
    ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
    for sp in ax.spines.values(): sp.set_visible(False)
    ax.set_title(d["sample"], fontsize=6, color="black", pad=1.5)
    if ylabel: ax.set_ylabel(ylabel, fontsize=6, color="black", rotation=90, labelpad=2)
    return sca

sca = None
for i, row in enumerate(data):
    for j, d in enumerate(row):
        sca = draw_img(fig.add_subplot(gs[i, j]), d,
                       ylabel=(("Healthy", "MASLD")[i] if j == 0 else None))

# shared colorbar — label as a short title ON TOP (rotated side-label would collide
# with the bar panel's y-axis label in the compact layout)
cax = fig.add_subplot(gs[:, 2])
cb = fig.colorbar(sca, cax=cax); cb.ax.tick_params(labelsize=5.5, width=0.4, length=2)
cb.outline.set_linewidth(0.3)
cax.set_title("log1p", fontsize=6, pad=3)

# cross-assay bars (spans both rows; col 3 is an empty spacer)
axb = fig.add_subplot(gs[:, 4])
labs = [x[0] for x in assays]; lfcs = [x[1] for x in assays]; padjs = [x[2] for x in assays]
xs = np.arange(len(assays))
axb.bar(xs, lfcs, width=0.66, color=[DARK if p < 0.05 else LIGHT for p in padjs])
axb.axhline(0, color="gray", lw=0.4)
for x, lf, p in zip(xs, lfcs, padjs):
    if p < 0.05:
        axb.text(x, lf + max(lfcs) * 0.03, "*", ha="center", va="bottom", fontsize=7)
axb.set_xticks(xs); axb.set_xticklabels(labs, fontsize=5.2)
axb.set_ylabel("log$_2$FC (MASLD vs control)", fontsize=6)
axb.set_ylim(0, max(lfcs) * 1.22)
axb.tick_params(axis="y", labelsize=6, length=2, width=0.4); axb.tick_params(axis="x", length=0)
for sp in ("top", "right"): axb.spines[sp].set_visible(False)

fig.suptitle(GENE, fontstyle="italic", fontsize=8, x=0.27, y=1.01)  # gene label over the images
fig.subplots_adjust(left=0.055, right=0.985, top=0.90, bottom=0.12)
fig.savefig(OUT, bbox_inches="tight")
print(f"saved {OUT}")

hmean = np.mean([d["v"].mean() for d in data[0]]); smean = np.mean([d["v"].mean() for d in data[1]])
print(f"[caption] {GENE} spatial + cross-assay. Left: LGALS4 expression (log1p normalized counts) over "
      f"registered H&E, GSE192741 human Visium — one section per patient (2 Healthy top, 2 MASLD bottom; "
      f"mean log1p {hmean:.2f} vs {smean:.2f}); GSE192741 has only these 4 human patients. Right: LGALS4 "
      f"log2FC (MASLD vs control) across assays — bulk RNA-seq {assays[0][1]:+.2f} (padj {assays[0][2]:.1e}, "
      f"pooled 846 samples/5 cohorts), proteomics {assays[1][1]:+.2f} (padj {assays[1][2]:.3f}), scRNA hep "
      f"{assays[2][1]:+.2f} (n.s.), GeoMx spatial {assays[3][1]:+.2f} (n.s.); * = padj<0.05. Disease direction "
      f"is established by bulk+proteomics; images illustrate it. No COLOC (atlas 0.92 off-trait platelet GWAS).")
