"""Fig 4H geometry/label check. Points (curves with ~3.9 pt bbox) inverted through the tick transform and matched to my
recomputed feature-first system coordinates (r4h_featurefirst_panel_estimand.tsv). Point colours sampled from a 300-dpi render.
Usage: python g4h_geometry.py <pdf> <tag>"""
import sys, subprocess, os, re
import numpy as np
from PIL import Image
from scipy.optimize import linear_sum_assignment
from scipy import stats
from pdfmarks import *
pdfp, tag = sys.argv[1], sys.argv[2]
out = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
E = pd.read_csv(f"{out}/r4h_featurefirst_panel_estimand.tsv", sep="\t")
pdf, pg, ch, w, r, l, c = load(pdfp)
num = w[w.text.str.fullmatch(r"-?\d+\.\d+;?")].copy(); num["v"] = num.text.str.rstrip(";").astype(float)
xt = num[(num.top > 140) & (num.yc < 156)]; yt = num[num.x1 < 33]
ticks_x = l[(l.top > 140) & (l.x0 == l.x1) & ((l.bottom - l.top) < 3)]
ticks_y = l[(l.x1 - l.x0 < 3) & (l.top == l.bottom)]
mx = [(t.x0, xt.loc[(xt.xc - t.x0).abs().idxmin(), "v"]) for _, t in ticks_x.iterrows()]
my = [(t.top, yt.loc[(yt.yc - t.top).abs().idxmin(), "v"]) for _, t in ticks_y.iterrows()]
ax, bx, rx = axis_fit([p for p, _ in mx], [v for _, v in mx]); ay, by, ry = axis_fit([p for p, _ in my], [v for _, v in my])
print("x ticks", mx, "resid", rx, "\ny ticks", my, "resid", ry)
pts = c[(c.w.between(3.5, 4.3)) & (c.h.between(3.5, 4.3)) & (c.yc > 30) | ((c.w.between(3.5, 4.3)) & (c.h.between(3.5, 4.3)) & (c.xc > 60))]
pts = pts[~((pts.w.round(1) == 3.0))]
key = c[(c.w.round(1) == 3.0) & (c.h.round(1) == 3.0)]
print("points", len(pts), "legend keys", len(key), key.fill.tolist())
pts = pts.assign(x=ax * pts.xc + bx, y=ay * pts.yc + by)
# render to sample colours
png = f"{out}/../render/_{tag}_h.png"
subprocess.run(["pdftoppm", "-r", "300", "-png", "-singlefile", pdfp, png[:-4]], check=True)
im = np.asarray(Image.open(png).convert("RGB")).astype(int); s = 300 / 72
def samp(xc, yc):
    X, Y = int(round(xc * s)), int(round(yc * s)); patch = im[Y - 2:Y + 3, X - 2:X + 3].reshape(-1, 3)
    return "#%02X%02X%02X" % tuple(np.median(patch, 0).astype(int))
pts["colour"] = [samp(a, b) for a, b in zip(pts.xc, pts.yc)]
# assignment: drawn point -> recomputed system
D = np.hypot(pts.x.values[:, None] - E.x.values[None], pts.y.values[:, None] - E.y.values[None])
ri, ci = linear_sum_assignment(D)
rng_x, rng_y = 0.33, 0.85
tol_x, tol_y = max(0.005 * rng_x, 0.002), max(0.005 * rng_y, 0.002)
res = []
for i, j in zip(ri, ci):
    p, e = pts.iloc[i], E.iloc[j]
    cls = "pos" if (e.x > 0 and e.y > 0) else ("neg" if (e.x < 0 and e.y < 0) else "opp")
    res.append(dict(tag=tag, community_id=e.community_id, label=e.figure_label, x_drawn=p.x, x_exp=e.x, y_drawn=p.y, y_exp=e.y,
                    dx=p.x - e.x, dy=p.y - e.y, off=(abs(p.x - e.x) > tol_x) or (abs(p.y - e.y) > tol_y), colour=p.colour, expected_class=cls))
R = pd.DataFrame(res)
print("tolerance x", tol_x, "y", tol_y, "; points checked", len(R), "; off", int(R.off.sum()), "; max|dx|", R.dx.abs().max(), "max|dy|", R.dy.abs().max())
print(R.groupby(["expected_class", "colour"]).size())
print(R[R.off].to_string())
# labels: each S-label word's nearest point should be that system
lab = w[w.text.str.fullmatch(r"S\d\d(S\d\d)?")]
lr = []
for _, t in lab.iterrows():
    for k, sl in enumerate(re.findall(r"S\d\d", t.text)):
        e = E[E.figure_label == sl]
        if len(e) == 0:
            lr.append(dict(label=sl, ok=False, note="no such system")); continue
        e = e.iloc[0]; px, py = (e.x - bx) / ax, (e.y - by) / ay
        dist = np.hypot(t.xc - px, t.yc - py)
        lr.append(dict(label=sl, label_xc=t.xc, label_yc=t.yc, point_xc=px, point_yc=py, dist_pt=dist))
L = pd.DataFrame(lr); print(L.round(2).to_string())
lead = l[(l.stroke == "#777777") | (l.stroke.astype(str).str.startswith("0.466"))]
txt = " ".join(w.text)
m = re.search(r"rho\s*=\s*(-?[\d.]+);?\s*(\d+)/(\d+)", txt)
print("printed rho/same/n:", m.groups() if m else None, "; recomputed rho", round(stats.spearmanr(E.x, E.y).correlation, 4),
      "same", int((np.sign(E.x) == np.sign(E.y)).sum()), "n", len(E))
ytitle = "".join(ch[ch.x1 < 18].sort_values(["x0", "top"], ascending=[True, False]).text)
print("x title:", " ".join(w[w.top > 155].sort_values(["top", "x0"]).text), "| y title:", ytitle, "| n units printed:", any(k in txt for k in ["367", "109", "258"]))
R.to_csv(f"{out}/g4h_{tag}.tsv", sep="\t", index=False); L.to_csv(f"{out}/g4h_{tag}_labels.tsv", sep="\t", index=False)
os.remove(png)
