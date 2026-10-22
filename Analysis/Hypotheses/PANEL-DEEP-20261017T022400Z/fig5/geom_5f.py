"""Figure 5F geometry check.
Right: invert dot x positions (ticks 0,.02,.04,.06 per facet) to Moran excess; compare to recompute/5f_recomputed_program_level.tsv;
       q labels vs recomputed BH q; colours per mark.
Left: parse each spot marker (fill colour + cumulative translation) from the content stream, invert colour through the
      colorbar image (I2, spans -2..2 by its tick marks), match spots to the overlay table by a fitted affine map of hires
      coordinates, compare to the table z clipped to [-2, 2]. Then Moran's I of the table z on a Visium hex graph for the section.
Usage: python geom_5f.py <pdf> <tag>
"""
import sys, re
import numpy as np
import pandas as pd
import pdfplumber, pypdf
from scipy.spatial import cKDTree

B = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5"
P = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
pdf, tag = sys.argv[1], sys.argv[2]
rec = pd.read_csv(f"{B}/recompute/5f_recomputed_program_level.tsv", sep="\t")
rec = rec[rec.design == "all16"]

def hexc(c):
    if c is None or c == "": return ""
    if isinstance(c, (int, float)): c = (c,)
    if len(c) == 1: c = (c[0],) * 3
    return "#%02X%02X%02X" % tuple(int(round(v * 255)) for v in c)

pl = pdfplumber.open(pdf).pages[0]
H = pl.height
words = pl.extract_words(extra_attrs=["size", "non_stroking_color"])
# facets from the dashed zero lines
zero = sorted([l for l in pl.lines if abs(l["x1"] - l["x0"]) < 0.01 and l["bottom"] - l["top"] > 30], key=lambda l: l["x0"])
ticks = [l for l in pl.lines if abs(l["x1"] - l["x0"]) < 0.01 and 1.5 < l["bottom"] - l["top"] < 2.5 and l["x0"] > 120]
labs = [w for w in words if re.fullmatch(r"\.?\d+|0", w["text"]) and w["x0"] > 120 and w["top"] > 80]
heads = sorted([w for w in words if w["top"] < 16 and w["x0"] > 120], key=lambda w: w["x0"])
facets = []
for z in zero:
    tk = sorted([t for t in ticks if z["x0"] - 1 <= t["x0"] <= z["x0"] + 80], key=lambda t: t["x0"])
    pr = np.array([(t["x0"], float(min(labs, key=lambda w: abs((w["x0"] + w["x1"]) / 2 - t["x0"]))["text"])) for t in tk])
    fit = np.polyfit(pr[:, 0], pr[:, 1], 1)
    head = " ".join(w["text"] for w in heads if z["x0"] - 30 < w["x0"] < z["x0"] + 80)
    facets.append({"x0": z["x0"], "fit": fit, "res": float(np.abs(np.polyval(fit, pr[:, 0]) - pr[:, 1]).max()), "head": head, "ticks": pr[:, 1].tolist()})
rowlab = {"GSE192741": None, "Vu_et_al_2025": None}
for w in words:
    if w["text"] == "Guilliams": rowlab["GSE192741"] = (w["top"] + w["bottom"]) / 2
    if w["text"] == "Vu": rowlab["Vu_et_al_2025"] = (w["top"] + w["bottom"]) / 2
dots = [c for c in pl.curves if 3 < (c["x1"] - c["x0"]) < 7]
out = []
for c in dots:
    cx, cy = (c["x0"] + c["x1"]) / 2, (c["top"] + c["bottom"]) / 2
    f = max([f for f in facets if f["x0"] - 2 <= cx], key=lambda f: f["x0"])
    ds = min(rowlab, key=lambda k: abs(rowlab[k] - cy))
    prog = "Stromal ECM (IGFBP7)" if f["head"].startswith("Stromal") else "Ductular injury (BICC1)"
    r = rec[(rec.dataset == ds) & (rec.program == prog)].iloc[0]
    val = np.polyval(f["fit"], cx)
    rng = f["ticks"][-1] - f["ticks"][0]
    stem = [l for l in pl.lines if abs(l["top"] - cy) < 0.5 and abs(l["bottom"] - cy) < 0.5 and l["x1"] - l["x0"] > 1]
    ql = [w for w in words if w["text"].startswith("q") and abs((w["top"] + w["bottom"]) / 2 - cy) < 3 and abs(w["x0"] - c["x1"]) < 4]
    qtxt = ql[0]["text"] if ql else ""
    qexp = r.q
    qprint = ("q<.001" if qexp < 0.001 else "q=" + ("%.3f" % qexp).lstrip("0"))
    out.append({"dataset": ds, "program": prog, "excess_drawn": val, "excess_rec": r.excess, "diff": val - r.excess, "tol": 0.005 * rng,
                "off": abs(val - r.excess) > 0.005 * rng, "stem_from_zero": bool(stem) and abs(np.polyval(f["fit"], min(l["x0"] for l in stem))) < 1e-4,
                "fill": hexc(c["non_stroking_color"]), "q_printed": qtxt, "q_expected": qprint, "q_ok": qtxt == qprint,
                "q_text_colour": hexc(ql[0]["non_stroking_color"]) if ql else ""})
D = pd.DataFrame(out)
print(tag, "facets", [(f["head"], f["ticks"], round(f["res"], 6)) for f in facets])
print(D.to_string())
D.to_csv(f"{B}/geom_checks/5f_dots_{tag}.tsv", sep="\t", index=False)

# ---------------- overlay spots
rd = pypdf.PdfReader(pdf)
pg = rd.pages[0]
cs = pg.get_contents().get_data().decode("latin1")
toks = re.findall(r"([-\d.]+) ([-\d.]+) ([-\d.]+) rg|1 0 0 1 ([-\d.]+) ([-\d.]+) cm /P0 Do", cs)
x = y = 0.0; col = None; spots = []
for t in toks:
    if t[0]:
        col = (float(t[0]), float(t[1]), float(t[2]))
    else:
        x += float(t[3]); y += float(t[4])
        spots.append((x, H - y) + col)
S = pd.DataFrame(spots, columns=["px", "py", "r", "g", "b"])
print(tag, "spot markers drawn", len(S))
imgs = pl.images
cb = [i for i in imgs if i["srcsize"][1] < 30]
cbimg = None
for im in pg.images:
    if im.image.size[1] < 30:
        cbimg = np.asarray(im.image.convert("RGB")).astype(float) / 255
row = cbimg[cbimg.shape[0] // 2]
# colorbar placement and tick mapping
cbx0, cbx1 = cb[0]["x0"], cb[0]["x1"]
ctk = sorted([l for l in pl.lines if abs(l["x1"] - l["x0"]) < 0.01 and l["top"] > 95], key=lambda l: l["x0"])
ctl = [w for w in words if w["top"] > 98 and re.fullmatch(r"[−-]?\d", w["text"])]
cpr = np.array([(t["x0"], float(min(ctl, key=lambda w: abs((w["x0"] + w["x1"]) / 2 - t["x0"]))["text"].replace("−", "-"))) for t in ctk])
cfit = np.polyfit(cpr[:, 0], cpr[:, 1], 1)
px_x = cbx0 + (np.arange(row.shape[0]) + 0.5) / row.shape[0] * (cbx1 - cbx0)
zgrid = np.polyval(cfit, px_x)
tree = cKDTree(row)
dist, idx = tree.query(S[["r", "g", "b"]].values)
S["z_drawn"] = zgrid[idx]; S["colour_dist"] = dist
ov = pd.read_csv(f"{P}/figures/main/fig5_molecular_context/panels/data/fig5f_spatial_histology_overlay.tsv.gz", sep="\t")
A = np.c_[ov.hires_x, ov.hires_y, np.ones(len(ov))]
# affine fit by iterative nearest-neighbour (start from bounding box map)
def bbox_map(src, dst):
    s = (dst.max(0) - dst.min(0)) / (src.max(0) - src.min(0)); return lambda q: (q - src.min(0)) * s + dst.min(0)
src = ov[["hires_x", "hires_y"]].values; dst = S[["px", "py"]].values
m = bbox_map(src, dst); pred = m(src)
for _ in range(5):
    tr = cKDTree(dst); dd, ii = tr.query(pred)
    coef, *_ = np.linalg.lstsq(A, dst[ii], rcond=None); pred = A @ coef
dd, ii = cKDTree(dst).query(pred)
print(tag, "affine residual max (pt)", round(float(dd.max()), 4), "unique matches", len(set(ii)), "of", len(ov))
ov["z_drawn"] = S.z_drawn.values[ii]
ov["z_exp"] = ov.all_cell_adjusted_program_score_z.clip(-2, 2)
ov["zdiff"] = ov.z_drawn - ov.z_exp
step = 4 / row.shape[0]
off = ov[ov.zdiff.abs() > max(0.02, step * 1.5)]
print(tag, "overlay spots checked", len(ov), "off by > max(0.02, 1.5 colorbar px)", len(off), "max |dz|", round(float(ov.zdiff.abs().max()), 4), "median", round(float(ov.zdiff.abs().median()), 4))
print(tag, "z clipped spots: >2", int((ov.all_cell_adjusted_program_score_z > 2).sum()), "<-2", int((ov.all_cell_adjusted_program_score_z < -2).sum()),
      "z range", ov.all_cell_adjusted_program_score_z.min().round(2), ov.all_cell_adjusted_program_score_z.max().round(2),
      "section mean z", ov.all_cell_adjusted_program_score_z.mean().round(3))
ov.to_csv(f"{B}/geom_checks/5f_overlay_spots_{tag}.tsv", sep="\t", index=False)
print(tag, "colorbar endpoints", row[0].round(3), row[len(row) // 2].round(3), row[-1].round(3))
