"""Figure 5B geometry check.

Axis transform from vector tick marks + tick labels; OLS line endpoints inverted;
stacked bar heights; raster scatter points (embedded image) checked against the
recomputed gene table by colour-class presence tests at expected positions and by
blob detection for stray marks.
Usage: python geom_5b.py <pdf> <tag>
"""
import sys
import io
import numpy as np
import pandas as pd
import pdfplumber
from pypdf import PdfReader
from PIL import Image
from scipy import ndimage

BASE = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5"
pdf, tag = sys.argv[1], sys.argv[2]
rec = pd.read_csv(f"{BASE}/recompute/5b_recomputed_genes.tsv", sep="\t")
out = []

with pdfplumber.open(pdf) as doc:
    p = doc.pages[0]
    words = p.extract_words(extra_attrs=["size"])
    lines = p.lines
    rects = p.rects
    img = p.images[0]

# ---- axis ticks: short black lines at x 24.66-26.156 (y ticks) and y 112.844-114.336 (x ticks)
yt = sorted([l for l in lines if abs(l["x1"] - l["x0"]) < 2 and abs(l["top"] - l["bottom"]) < 0.01 and l["x0"] < 30],
            key=lambda l: l["top"])
xt = sorted([l for l in lines if abs(l["x1"] - l["x0"]) < 0.01 and 0.5 < abs(l["bottom"] - l["top"]) < 3],
            key=lambda l: l["x0"])
# tick labels: words left of y axis, below x axis
import re
ylab = [w for w in words if w["x1"] < 25 and w["top"] < 116 and re.fullmatch(r"[-\u2212]?\d+(\.\d+)?", w["text"])]
xlab = [w for w in words if 114 < w["top"] < 123 and re.fullmatch(r"[-\u2212]?\d+(\.\d+)?", w["text"])]

def match(ticks, labels, key_t, key_l):
    pairs = []
    for t in ticks:
        tv = key_t(t)
        best = min(labels, key=lambda w: abs(key_l(w) - tv))
        pairs.append((tv, float(best["text"].replace("−", "-"))))
    return np.array(pairs)

ypairs = match(yt, ylab, lambda t: t["top"], lambda w: (w["top"] + w["bottom"]) / 2)
xpairs = match(xt, xlab, lambda t: t["x0"], lambda w: (w["x0"] + w["x1"]) / 2)
ay = np.polyfit(ypairs[:, 0], ypairs[:, 1], 1)  # data = a*top + b
ax = np.polyfit(xpairs[:, 0], xpairs[:, 1], 1)
yres = np.abs(np.polyval(ay, ypairs[:, 0]) - ypairs[:, 1]).max()
xres = np.abs(np.polyval(ax, xpairs[:, 0]) - xpairs[:, 1]).max()
print(tag, "y tick pairs", ypairs.tolist(), "max resid", yres)
print(tag, "x tick pairs", xpairs.tolist(), "max resid", xres)
xr = (np.polyval(ax, 169.07) - np.polyval(ax, 26.156))
yr = (np.polyval(ay, 7.473) - np.polyval(ay, 112.844))
print(tag, "axis ranges x", xr, "y", yr)

# ---- OLS line: the longest dark line
ols = max(lines, key=lambda l: np.hypot(l["x1"] - l["x0"], l["bottom"] - l["top"]))
(x1, t1), (x2, t2) = ols["pts"][0], ols["pts"][-1]
X1, Y1, X2, Y2 = np.polyval(ax, x1), np.polyval(ay, t1), np.polyval(ax, x2), np.polyval(ay, t2)
slope = (Y2 - Y1) / (X2 - X1)
icpt = Y1 - slope * X1
summ = pd.read_csv(f"{BASE}/recompute/5b_recomputed_summary.tsv", sep="\t", header=None, index_col=0)[1]
out.append(("OLS slope", slope, float(summ["ols_slope"]), abs(slope - float(summ["ols_slope"])) <= 0.005 * 1.0))
out.append(("OLS intercept", icpt, float(summ["ols_intercept"]), abs(icpt - float(summ["ols_intercept"])) <= 0.005 * yr))

# ---- stacked bar
bars = [r for r in rects if r.get("fill") and r["width"] < 12 and r["height"] > 5 and r["x0"] > 180]
bars = sorted(bars, key=lambda r: r["top"])
tot = sum(r["height"] for r in bars)
nsig = float(summ["n_protein_sig"])
exp = {"top": float(summ["n_discordant"]) / nsig, "bottom": float(summ["n_concordant"]) / nsig}
out.append(("bar discordant fraction (top)", bars[0]["height"] / tot, exp["top"], abs(bars[0]["height"] / tot - exp["top"]) < 0.005))
out.append(("bar concordant fraction (bottom)", bars[1]["height"] / tot, exp["bottom"], abs(bars[1]["height"] / tot - exp["bottom"]) < 0.005))
bar_cols = [tuple(np.round(np.array(r["non_stroking_color"]) * 255).astype(int)) for r in bars]
print(tag, "bar colours", bar_cols)

# ---- raster points
reader = PdfReader(pdf)
page = reader.pages[0]
ims = page.images
pil = ims[0].image.convert("RGB")
arr = np.asarray(pil).astype(int)
H, W = arr.shape[:2]
sx = (img["x1"] - img["x0"]) / W
sy = (img["bottom"] - img["top"]) / H

def px_to_data(col, row):
    x_pt = img["x0"] + (col + 0.5) * sx
    top_pt = img["top"] + (row + 0.5) * sy
    return np.polyval(ax, x_pt), np.polyval(ay, top_pt)

def data_to_px(x, y):
    # invert linear maps
    x_pt = (x - ax[1]) / ax[0]
    top_pt = (y - ay[1]) / ay[0]
    return (x_pt - img["x0"]) / sx - 0.5, (top_pt - img["top"]) / sy - 0.5

blue = np.array([0x3B, 0x6E, 0xA5]); red = np.array([0xB0, 0x5A, 0x5A])
def near(c, tol=40):
    return np.abs(arr - c).sum(axis=2) < tol
mb, mr = near(blue), near(red)
nonwhite = (arr.sum(axis=2) < 3 * 245)
res = []
rad = 3  # px (~0.54 pt)
for _, g in rec.iterrows():
    cx, cy = data_to_px(g.rna_logFC, g.logFC)
    c0, r0 = int(round(cx)), int(round(cy))
    if not (0 <= c0 < W and 0 <= r0 < H):
        res.append((g.gene, g["class"], "outside_image"))
        continue
    sl = (slice(max(r0 - rad, 0), r0 + rad + 1), slice(max(c0 - rad, 0), c0 + rad + 1))
    if g["class"] == "concordant":
        ok = mb[sl].any()
    elif g["class"] == "discordant":
        ok = mr[sl].any()
    else:
        ok = nonwhite[sl].any()
    res.append((g.gene, g["class"], "ok" if ok else "missing"))
res = pd.DataFrame(res, columns=["gene", "class", "status"])
print(tag, res.groupby(["class", "status"]).size().to_dict())

# blobs of sig colours not explained by any expected point of that class
def stray(mask, cls):
    lab, n = ndimage.label(mask)
    cents = ndimage.center_of_mass(mask, lab, range(1, n + 1))
    exp = rec[rec["class"] == cls]
    ex = np.array([data_to_px(a, b) for a, b in zip(exp.rna_logFC, exp.logFC)])
    bad = []
    for (r, c) in cents:
        d = np.sqrt(((ex[:, 0] - c) ** 2 + (ex[:, 1] - r) ** 2)).min()
        if d > 12:  # px; blob centre > ~2 pt from any expected point of the class
            bad.append((r, c, d, *px_to_data(c, r)))
    return n, bad
nb, badb = stray(mb, "concordant")
nr, badr = stray(mr, "discordant")
print(tag, "blue blobs", nb, "stray", badb[:10])
print(tag, "red blobs", nr, "stray", badr[:10])

# axis limits vs data range
xlim = (np.polyval(ax, 26.156), np.polyval(ax, 169.07))
ylim = (np.polyval(ay, 112.844), np.polyval(ay, 7.473))
clipped = rec[(rec.rna_logFC < xlim[0]) | (rec.rna_logFC > xlim[1]) | (rec.logFC < ylim[0]) | (rec.logFC > ylim[1])]
print(tag, "xlim", xlim, "ylim", ylim, "points outside plotted range", len(clipped), clipped[["gene", "rna_logFC", "logFC", "class"]].to_dict("records"))
for o in out:
    print(tag, "CHECK", o)
res.to_csv(f"{BASE}/geom_checks/5b_points_{tag}.tsv", sep="\t", index=False)
pd.DataFrame(out, columns=["mark", "drawn", "expected", "ok"]).to_csv(f"{BASE}/geom_checks/5b_marks_{tag}.tsv", sep="\t", index=False)
