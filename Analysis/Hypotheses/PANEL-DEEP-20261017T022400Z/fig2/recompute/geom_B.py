"""Fig 2B geometry: recover (PIP, PP.H4) of every data circle from tick positions; match to recomputed genes.
usage: geom_B.py <dump prefix> <tag>"""
import sys, json, numpy as np, pandas as pd
from scipy.optimize import linear_sum_assignment
pre, tag = sys.argv[1], sys.argv[2]
cv = pd.read_csv(f"{pre}.curves.tsv", sep="\t"); ln = pd.read_csv(f"{pre}.lines.tsv", sep="\t"); wd = pd.read_csv(f"{pre}.words.tsv", sep="\t")
# axis ticks: short black lines; x ticks vertical below the axis, y ticks horizontal left of the axis
tk = ln[(ln.stroke == "#000000") & (((ln.w < 0.01) & ln.h.between(1, 2)) | ((ln.h < 0.01) & ln.w.between(1, 2)))]
xt = sorted(tk[tk.w < 0.01].x0.tolist()); yt = sorted(tk[tk.h < 0.01].top.tolist(), reverse=True)
xlab = wd[(wd.top > 155) & (wd.top < 166) & (wd.x0 > 20) & wd.text.astype(str).str.match(r"^[0-9]+\.[0-9]+$")].sort_values("x0"); ylab = wd[(wd.x0 > 10) & (wd.x0 < 21) & wd.text.astype(str).str.match(r"^[0-9]+\.[0-9]+$")].sort_values("top", ascending=False)
xv = [float(t) for t in xlab.text]; yv = [float(t) for t in ylab.text]
ax = np.polyfit(xt, xv, 1); ay = np.polyfit(yt, yv, 1)
dat = cv[cv.top > 36].copy()  # legend circles sit above the plot frame
dat["cx"] = (dat.x0 + dat.x1) / 2; dat["cy"] = dat.top
dat["pip"] = np.polyval(ax, dat.cx); dat["pp4"] = np.polyval(ay, dat.cy)
rec = pd.read_csv("out_ABCEF/B_shown_both.tsv", sep="\t")
# cost in axis units
D = np.hypot(dat.pip.values[:, None] - rec.pip_max.values[None, :], dat.pp4.values[:, None] - rec.s4.values[None, :])
r, c = linear_sum_assignment(D)
d = D[r, c]
tolx = 0.005; toly = 0.005
dx = np.abs(dat.pip.values[r] - rec.pip_max.values[c]); dy = np.abs(dat.pp4.values[r] - rec.s4.values[c])
off = (dx > tolx) | (dy > toly)
m = pd.DataFrame({"gene": rec.gene.values[c], "rec_pip": rec.pip_max.values[c], "rec_pp4": rec.s4.values[c],
                  "drawn_pip": dat.pip.values[r], "drawn_pp4": dat.pp4.values[r], "fill": dat.fill.values[r], "w": dat.w.values[r], "off": off})
m.to_csv(f"out_geom/B_{tag}_matched.tsv", sep="\t", index=False)
out = {"tag": tag, "x_ticks": list(zip(xt, xv)), "y_ticks": list(zip(yt, yv)), "n_data_marks": len(dat), "n_recomputed": len(rec),
       "n_matched": len(r), "n_off": int(off.sum()), "max_dx": float(dx.max()), "max_dy": float(dy.max()),
       "fill_counts": dat.fill.value_counts().to_dict(), "header_words": wd[wd.top < 35].text.tolist(),
       "axis_titles": " ".join(wd[(wd.top > 166)].text.astype(str)),
       "off_genes": m[m.off].gene.tolist()[:20]}
json.dump(out, open(f"out_geom/B_{tag}.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
