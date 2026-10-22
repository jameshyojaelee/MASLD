"""Fig 4B geometry/label check: legend colour -> label; per-group drawn point counts vs atlas proportions (50,000-cell
uniform subsample expectation, binomial SD); affine fit of atlas per-group median UMAP to drawn per-group medians
(residuals); density overlay extent and colour-bar range. Usage: python g4b_geometry.py <pdf> <tag>"""
import sys
import numpy as np
from pdfmarks import *
pdfp, tag = sys.argv[1], sys.argv[2]
out = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
P = pd.read_csv(f"{out}/r4b_celltype_props.tsv", sep="\t", index_col=0)
umb = {"Hepatocytes": "Hepatocytes", "T cells": "T cells", "Endothelial cells": "Endothelial cells", "Mono+mono derived cells": "Monocytes",
       "Fibroblasts": "Fibroblasts", "Macrophages": "Macrophages", "Cholangiocytes": "Cholangiocytes", "Circulating NK/NKT": "NK cells",
       "Resident NK": "NK cells", "Neutrophils": "Granulocytes", "Basophils": "Granulocytes", "B cells": "B cells", "cDC2s": "Dendritic cells",
       "pDCs": "Dendritic cells", "cDC1s": "Dendritic cells", "Plasma cells": "Plasma cells"}
P["group"] = P.index.map(umb)
P["wu1"], P["wu2"] = P.u1 * P.n_cells, P.u2 * P.n_cells
Gp = P.groupby("group").agg(n=("n_cells", "sum"), wu1=("wu1", "sum"), wu2=("wu2", "sum"))
Gp["frac"] = Gp.n / Gp.n.sum(); Gp["u1"], Gp["u2"] = Gp.wu1 / Gp.n, Gp.wu2 / Gp.n
pdf, pg, ch, w, r, l, c = load(pdfp)
keys = c[c.w.between(2, 2.3)].sort_values("yc")
# legend labels: first word of each legend row
rows = w[w.x0 > keys.x1.max()].copy(); rows["row"] = rows.top.round(0)
lab = rows.groupby("row").apply(lambda g: " ".join(g.sort_values("x0").text)).reset_index(name="label")
km = {}
for _, k in keys.iterrows():
    rr = lab.iloc[(lab.row - k.top).abs().argmin()]
    km[k.fill] = rr.label
pts = c[c.w.between(0.9, 1.1)].copy(); pts["group"] = pts.fill.map(km)
print("legend:", km)
R = pts.groupby("group").agg(n_drawn=("xc", "size"), x_med=("xc", "median"), y_med=("yc", "median")).join(Gp[["frac", "u1", "u2"]])
R["n_exp"] = R.frac * len(pts); R["z"] = (R.n_drawn - R.n_exp) / np.sqrt(R.n_exp * (1 - R.frac))
# affine fit atlas (u1,u2) -> page (x,y) weighted equally over groups
A = np.column_stack([R.u1, R.u2, np.ones(len(R))])
cx, *_ = np.linalg.lstsq(A, R.x_med, rcond=None); cy, *_ = np.linalg.lstsq(A, R.y_med, rcond=None)
R["resid_pt"] = np.hypot(A @ cx - R.x_med, A @ cy - R.y_med)
print(R.round(3).to_string())
print("points", len(pts), "| unmapped colours", int(pts.group.isna().sum()), "| max |z| count", round(R.z.abs().max(), 2), "| max centroid resid (pt)", round(R.resid_pt.max(), 2))
tiles = r[(r.w.between(0.2, 0.22)) & (r.h < 1)]
cbar = r[(r.w.between(0.2, 0.22)) & (r.h > 4)]
print("density tiles", len(tiles), "| non-white", int((tiles.fill != "#FFFFFF").sum()), "| x", round(tiles.x0.min(), 2), round(tiles.x1.max(), 2), "| y", round(tiles.top.min(), 2), round(tiles.bottom.max(), 2))
cbt = l[l.stroke == "#FFFFFF"].drop_duplicates("x0"); num = w[w.text.str.fullmatch(r"-?\d+")].copy(); num["v"] = num.text.astype(float)
mc = [(t.x0, num.loc[(num.xc - t.x0).abs().idxmin(), "v"]) for _, t in cbt.iterrows()]
a, b, _ = axis_fit(*zip(*mc))
print("colour bar ticks", mc, "| range", round(a * cbar.x0.min() + b, 3), round(a * cbar.x1.max() + b, 3), "| segments", len(cbar))
# fraction of non-white tiles more than 3 pt from any drawn point
from scipy.spatial import cKDTree
t = cKDTree(pts[["xc", "yc"]].values)
nw = tiles[tiles.fill != "#FFFFFF"]
d, _ = t.query(np.column_stack([(nw.x0 + nw.x1) / 2, (nw.top + nw.bottom) / 2]))
print("coloured density tiles farther than 3 pt from any drawn cell:", int((d > 3).sum()), "of", len(nw), f"({(d > 3).mean():.1%})")
txt = " ".join(w.text)
print("axis titles:", "UMAP 1" in txt, "".join(ch[ch.x1 < 11].sort_values("top", ascending=False).text), "| n donors printed:", any(k in txt for k in ["14", "36", "donor"]))
R.to_csv(f"{out}/g4b_{tag}.tsv", sep="\t")
