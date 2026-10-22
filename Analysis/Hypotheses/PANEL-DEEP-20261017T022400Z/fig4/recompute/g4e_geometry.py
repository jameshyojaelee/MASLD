"""Fig 4E geometry/label check: 224 points (112 programs x fibrosis/NAS) inverted through each facet's ticks and matched
(optimal assignment) to my recomputed REML+Knapp-Hartung meta (r4e_meta_recomputed.tsv, job 9842): x = beta, y = -log10 BH q.
Usage: python g4e_geometry.py <pdf> <tag>"""
import sys, subprocess, os
import numpy as np
from PIL import Image
from scipy.optimize import linear_sum_assignment
from pdfmarks import *
pdfp, tag = sys.argv[1], sys.argv[2]
out = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
M = pd.read_csv(f"{out}/r4e_meta_recomputed.tsv", sep="\t")
M["mlq"] = -np.log10(M.q_kh)
pdf, pg, ch, w, r, l, c = load(pdfp)
num = w[w.text.str.fullmatch(r"-?\d+(\.\d+)?")].copy(); num["v"] = num.text.astype(float)
# facet split: gap between the two facet frames (left frame ends at max x of left spine-bottom line)
frames = l[(l.top == l.bottom) & ((l.x1 - l.x0) > 50) & (l.stroke.astype(str).isin(["0.0", "0", "#000000"]))].sort_values("x0")
mid = (frames.x1.iloc[0] + frames.x0.iloc[-1]) / 2
png = f"{out}/../render/_{tag}_e.png"
subprocess.run(["pdftoppm", "-r", "300", "-png", "-singlefile", pdfp, png[:-4]], check=True)
im = np.asarray(Image.open(png).convert("RGB")).astype(int); s = 300 / 72
def samp(xc, yc):
    X, Y = int(round(xc * s)), int(round(yc * s)); return "#%02X%02X%02X" % tuple(np.median(im[Y - 1:Y + 2, X - 1:X + 2].reshape(-1, 3), 0).astype(int))
allres = []
for facet, axname, left in [("Fibrosis", "fibrosis", True), ("NAS", "nas", False)]:
    sel = (lambda d, col: d[col] < mid) if left else (lambda d, col: d[col] >= mid)
    tx = l[sel(l, "x0") & (l.x0 == l.x1) & ((l.bottom - l.top) < 2.5) & (l.top > 70)]
    ty = l[sel(l, "x0") & (l.top == l.bottom) & ((l.x1 - l.x0) < 2.5)]
    xl = num[sel(num, "xc") & (num.top > 80) & (num.top < 88)]
    yl = num[(num.xc < 33) & (num.top < 85)]
    mx = [(t.x0, xl.loc[(xl.xc - t.x0).abs().idxmin(), "v"]) for _, t in tx.iterrows()]
    my = [(t.top, yl.loc[(yl.yc - t.top).abs().idxmin(), "v"]) for _, t in ty.iterrows()]
    ax, bx, rx = axis_fit(*zip(*mx)); ay, by, ry = axis_fit(*zip(*my))
    pts = c[sel(c, "xc") & c.w.between(1.8, 2.2)].copy()
    pts["x"], pts["y"] = ax * pts.xc + bx, ay * pts.yc + by
    pts["colour"] = [samp(a, b) for a, b in zip(pts.xc, pts.yc)]
    E = M[M.axis == axname].reset_index(drop=True)
    D = np.hypot(pts.x.values[:, None] - E.beta.values[None], pts.y.values[:, None] - E.mlq.values[None])
    ri, ci = linear_sum_assignment(D)
    xr, yr = 0.4 - (-0.3), 2.4
    tolx, toly = 0.005 * xr, 0.005 * yr
    for i, j in zip(ri, ci):
        p, e = pts.iloc[i], E.iloc[j]
        allres.append(dict(tag=tag, facet=facet, feature_id=e.feature_id, x_drawn=p.x, x_exp=e.beta, y_drawn=p.y, y_exp=e.mlq,
                           dx=p.x - e.beta, dy=p.y - e.mlq, off=(abs(p.x - e.beta) > tolx) or (abs(p.y - e.mlq) > toly), colour=p.colour))
    thr = l[sel(l, "x0") & (l.top == l.bottom) & ((l.x1 - l.x0) > 30) & (l.top > 20) & (l.top < 75)]
    print(facet, "x ticks", mx, "y ticks", my, "| points", len(pts), "expected", len(E), "| threshold line y =", [round(ay * t + by, 4) for t in thr.top],
          "| zero line x =", [round(ax * v + bx, 4) for v in l[sel(l, "x0") & (l.x0 == l.x1) & ((l.bottom - l.top) > 30) & (l.stroke != "0.0")].x0])
    print("  drawn x range", round(pts.x.min(), 3), round(pts.x.max(), 3), "; recomputed", round(E.beta.min(), 3), round(E.beta.max(), 3),
          "| y max drawn", round(pts.y.max(), 3), "recomputed", round(E.mlq.max(), 3), "| calls q<0.05:", int((E.q_kh < 0.05).sum()))
R = pd.DataFrame(allres)
print("points checked", len(R), "off", int(R.off.sum()), "max|dx|", R.dx.abs().max(), "max|dy|", R.dy.abs().max())
print(R.colour.value_counts().head())
print(R[R.off].head(20).to_string())
txt = " ".join(w.text)
print("axis titles: x fibrosis", "stage-code" in txt, "x NAS", "NAS point" in txt, "; y title chars:", "".join(ch[ch.x1 < 28].sort_values("top", ascending=False).text),
      "; n printed:", any(k in txt for k in ["364", "n ="]))
R.to_csv(f"{out}/g4e_{tag}.tsv", sep="\t", index=False)
os.remove(png)
