"""Fig 4F geometry/label check: 10 trajectories (2 programs x 5 cohorts), 9 vertices each, inverted through each facet's
ticks and compared to my recomputed window means (r4f_windows.tsv, last window closed). Cohort identity from the legend
key colour + dash. Usage: python g4f_geometry.py <pdf> <tag>"""
import sys
import numpy as np
import pdfplumber
from pdfmarks import *
pdfp, tag = sys.argv[1], sys.argv[2]
out = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
W = pd.read_csv(f"{out}/r4f_windows.tsv", sep="\t"); W = W[W.rule == "last_closed"]
prog = {"Ductular": "hotspot_hepatocytes_48f39dd4d817a10e", "Stromal": "hotspot_hepatocytes_f05c535ae5bbc0b9"}
pdf, pg, ch, w, r, l, c = load(pdfp)
raw = pdfplumber.open(pdfp).pages[0]
num = w[w.text.str.fullmatch(r"-?\d+(\.\d+)?%?")].copy(); num["v"] = num.text.str.rstrip("%").astype(float)
num.loc[num.text.str.endswith("%"), "v"] /= 100
# legend: short coloured key lines and the cohort word to their right
keys = l[(l.top == l.bottom) & ((l.x1 - l.x0) < 8) & (l.top > 70)]
leg = {}
for _, k in keys.iterrows():
    lab = w[(w.text.str.startswith("GSE")) & ((w.yc - k.top).abs() < 3) & (w.x0 > k.x1) & (w.x0 - k.x1 < 6)]
    leg[k.stroke] = lab.text.iloc[0] if len(lab) else None
dash = {}
for cv in raw.curves:
    if (cv["x1"] - cv["x0"]) > 80:
        pass
print("legend", leg)
facets = {"Ductular": w[w.text == "Ductular"].xc.iloc[0], "Stromal": w[w.text == "Stromal"].xc.iloc[0]}
traj = [cv for cv in raw.curves if (cv["x1"] - cv["x0"]) > 80]
res = []
for cv in traj:
    stroke = col(cv["stroking_color"]); cohort = leg.get(stroke)
    is_dashed = bool(cv.get("dash") and cv["dash"][0])
    xs = np.array([p[0] for p in cv["pts"]]); ys = np.array([p[1] for p in cv["pts"]])
    fac = min(facets, key=lambda f: abs(facets[f] - xs.mean()))
    # facet ticks
    lo, hi = xs.min() - 5, xs.max() + 5
    tx = l[(l.x0 == l.x1) & (l.bottom - l.top < 2) & (l.x0 > lo) & (l.x0 < hi)]
    xl = num[num.text.str.endswith("%") & (num.xc > lo) & (num.xc < hi)]
    mx = [(t.x0, xl.loc[(xl.xc - t.x0).abs().idxmin(), "v"]) for _, t in tx.iterrows()]
    yaxis_x = l[(l.x0 == l.x1) & (l.bottom - l.top > 30) & (l.x0 < xs.min() + 2) & (l.x0 > xs.min() - 5)].x0.iloc[0]
    ty = l[(l.top == l.bottom) & (l.x1 - l.x0 < 2) & ((l.x1 - yaxis_x).abs() < 0.1)]
    yl = num[(~num.text.str.endswith("%")) & (num.x1 < yaxis_x) & (num.x1 > yaxis_x - 15)]
    my = [(t.top, yl.loc[(yl.yc - t.top).abs().idxmin(), "v"]) for _, t in ty.iterrows()]
    ax, bx, rx = axis_fit(*zip(*mx)); ay, by, ry = axis_fit(*zip(*my))
    xv, yv = ax * xs + bx, ay * ys + by
    E = W[(W.program_uid == prog[fac]) & (W.dataset == cohort)].sort_values("centre")
    yr = (ay * 50.86 + by) - (ay * 8.99 + by)
    for k in range(len(xv)):
        e = E.iloc[k]
        res.append(dict(tag=tag, facet=fac, cohort=cohort, colour=stroke, dashed=is_dashed, centre=e.centre, x_drawn=xv[k], y_drawn=yv[k], y_exp=e.mean_outcome_z,
                        n_window=e.n, se=e.se, dx=xv[k] - e.centre, dy=yv[k] - e.mean_outcome_z,
                        off=(abs(xv[k] - e.centre) > 0.005 * 0.8) or (abs(yv[k] - e.mean_outcome_z) > max(0.005 * abs(yr), 0.005))))
R = pd.DataFrame(res)
print("trajectories", len(traj), "| vertices checked", len(R), "| off", int(R.off.sum()), "| max|dx|", R.dx.abs().max(), "max|dy|", R.dy.abs().max())
print(R.groupby(["cohort", "colour", "dashed"]).size())
print(R[R.off].to_string())
txt = " ".join(w.text)
print("x title:", "Continuum percentile" in txt, "; y title:", "".join(ch[ch.x1 < 16].sort_values(["x0", "top"], ascending=[True, False]).text), "; n printed:", "n =" in txt or "n=" in txt)
print("printed y range (Ductular) ticks", sorted(num[(num.x1 < 29)].v.tolist()), "| Stromal", sorted(num[(num.x1 > 130) & (num.x1 < 146)].v.tolist()))
R.to_csv(f"{out}/g4f_{tag}.tsv", sep="\t", index=False)
