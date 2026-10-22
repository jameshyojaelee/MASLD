"""Fig 4C geometry check: recover every data mark value from the PDF and compare to r4c_rows.tsv.
Usage: python g4c_geometry.py <pdf> <tag>"""
import sys, re
from pdfmarks import *
pdfp, tag = sys.argv[1], sys.argv[2]
out = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute"
rec = pd.read_csv(f"{out}/r4c_rows.tsv", sep="\t")
rec["lin"] = rec.cell_type.map({"hepatocytes": "Hep", "fibroblasts": "Fib", "macrophages": "Mac", "cholangiocytes": "Chol"})
rec["short"] = rec.module_name.str.replace(r"\s*\(.*\)$", "", regex=True)
pdf, pg, ch, w, r, l, c = load(pdfp)
res = []

# --- row labels: words left of the bracket, grouped by line
lab = w[w.x1 < 106].copy()
lab["yr"] = lab.yc.round(0)
labels = lab.sort_values("x0").groupby("yr").text.apply(lambda s: " ".join(s)).reset_index()
labels = labels[labels.text.str.contains("·")]
# strip lineage tag appended on the same line
labels["text"] = labels.text.str.replace(r"\s+(Hep|Fib|Mac|Chol)$", "", regex=True)
labels["module"] = labels.text.str.extract(r"^(\d+)\s*·").astype(float)
labels["name"] = labels.text.str.replace(r"^\d+\s*·\s*", "", regex=True)
labels = labels.sort_values("yr").reset_index(drop=True)
# lineage by bracket position: brackets are vertical #404040 lines
br = l[(l.stroke == "#404040") & (l.w < 0.1)].sort_values("top")
br_lin = []
lw_ = w[w.text.isin(["Hep", "Fib", "Mac", "Chol"]) & (w.x0 > 60) & (w.x1 < 110)]
for _, b in br.iterrows():
    t = lw_[(lw_.yc > b.top) & (lw_.yc < b.bottom)].text.tolist()
    br_lin.append((b.top, b.bottom, t[0] if t else None))
def lin_of(y):
    for t0, t1, name in br_lin:
        if t0 - 1 <= y <= t1 + 1:
            return name
labels["lin"] = labels.yr.apply(lin_of)
print(labels.to_string())

# --- x axis of forest: ticks are short vertical black lines below axis, labels -1 0 1
ticks = l[(l.stroke == "#000000") & (l.w < 0.1) & (l.h < 4)].sort_values("x0")
tl = w[w.text.isin(["-1", "0", "1"]) & (w.top > ticks.top.min())].sort_values("x0")
a, b, rr = axis_fit(ticks.x0.values, tl.text.astype(float).values)
print("forest x-axis: ticks", ticks.x0.round(2).tolist(), "labels", tl.text.tolist(), "resid", rr)
X = lambda px: a * px + b
pt_per_unit = 1 / a

# --- points (curves), black CI, gray interval
pts = c.copy(); pts["val"] = X(pts.xc)
blk = l[(l.stroke == "#262626") & (l.h < 0.1) & (l.x0 < 230)].copy()
gry = l[(l.stroke == "#D1D1D1") & (l.h < 0.1) & (l.x0 < 230)].copy()
# heat tiles
tiles = r[(r.w > 12) & (r.w < 13.5) & (r.h > 5) & (r.h < 6)].copy()
colx = sorted(tiles.x0.round(1).unique())
stage_of = {colx[0]: "Healthy", colx[1]: "Steatosis", colx[2]: "Steatohepatitis"}
# colour bar
cb = r[(r.w > 6.9) & (r.w < 7.2) & (r.x0 > 230)].sort_values("top")
cbt = l[(l.stroke == "#FFFFFF") & (l.x0 > 230) & (l.x0 < 236)].sort_values("top")
cbl = w[w.text.isin(["1", "0.5", "0", "-0.5", "-1"]) & (w.x0 > 240) & (w.top < 110)].sort_values("top")
ca, cbb, cr = axis_fit(cbt.top.values, cbl.text.astype(float).values)
cb["v"] = ca * cb.yc + cbb
cb["rgb"] = cb.fill.apply(hex2rgb)
sat_hi = cb[cb.fill == cb.fill.iloc[0]].v.min(); sat_lo = cb[cb.fill == cb.fill.iloc[-1]].v.max()
print(f"colour bar: span {cb.v.max():.3f} to {cb.v.min():.3f}; colour saturates at >= {sat_hi:.3f} and <= {sat_lo:.3f}; tick resid {cr:.2g}")

def tile_val(hexc):
    rgb = hex2rgb(hexc)
    d = np.array([np.abs(rgb - x).sum() for x in cb.rgb])
    i = d.argmin()
    # interpolate between neighbours on the gradient for precision
    return cb.v.iloc[i], d[i]

for _, lr in labels.iterrows():
    y = lr.yr
    m = rec[(rec.lin == lr.lin) & (rec.module == lr.module)]
    if len(m) != 1:
        res.append(dict(tag=tag, mark="row_label", row=lr.text, expected="no unique recomputed row", drawn=lr.text, off=True)); continue
    m = m.iloc[0]
    name_ok = m.short.lower().startswith(lr["name"].lower()[:12])
    res.append(dict(tag=tag, mark="row_label", row=lr.text, expected=f"{m.lin} {int(m.module)} {m.short}", drawn=lr.text, off=not name_ok))
    p = pts[(pts.yc - y).abs() < 1.5]
    if len(p) == 1:
        v = p.val.iloc[0]
        res.append(dict(tag=tag, mark="point_beta", row=lr.text, expected=round(m.gls_beta, 4), drawn=round(v, 4),
                        off=abs(v - m.gls_beta) > 0.016, colour=p.fill.iloc[0]))
    else:
        res.append(dict(tag=tag, mark="point_beta", row=lr.text, expected=round(m.gls_beta, 4), drawn=f"{len(p)} points", off=True))
    bl = blk[(blk.yc - y).abs() < 1.5]
    if len(bl) == 1:
        lo, hi = X(bl.x0.iloc[0]), X(bl.x1.iloc[0])
        res.append(dict(tag=tag, mark="ci95_low", row=lr.text, expected=round(m.gls_ci_low, 4), drawn=round(lo, 4), off=abs(lo - m.gls_ci_low) > 0.016))
        res.append(dict(tag=tag, mark="ci95_high", row=lr.text, expected=round(m.gls_ci_high, 4), drawn=round(hi, 4), off=abs(hi - m.gls_ci_high) > 0.016))
    else:
        res.append(dict(tag=tag, mark="ci95", row=lr.text, expected="1 black bar", drawn=f"{len(bl)}", off=True))
    gr = gry[(gry.yc - y).abs() < 1.5]
    if len(gr) == 1:
        lo, hi = X(gr.x0.iloc[0]), X(gr.x1.iloc[0])
        res.append(dict(tag=tag, mark="gray_low", row=lr.text, expected=round(m.gray_lo, 4), drawn=round(lo, 4), off=abs(lo - m.gray_lo) > 0.016))
        res.append(dict(tag=tag, mark="gray_high", row=lr.text, expected=round(m.gray_hi, 4), drawn=round(hi, 4), off=abs(hi - m.gray_hi) > 0.016))
    else:
        res.append(dict(tag=tag, mark="gray", row=lr.text, expected="1 gray bar", drawn=f"{len(gr)}", off=True))
    tt = tiles[(tiles.top < y) & (tiles.bottom > y)]
    for _, t in tt.iterrows():
        st = stage_of[round(t.x0, 1)]
        exp = m[f"med_c_mean_{st}"]
        v, dist = tile_val(t.fill)
        exp_disp = min(max(exp, sat_lo), sat_hi)
        v_disp = min(max(v, sat_lo), sat_hi)
        res.append(dict(tag=tag, mark=f"tile_{st}", row=lr.text, expected=round(exp, 4), drawn=round(v_disp, 4),
                        off=abs(v_disp - exp_disp) > 0.05, colour=t.fill, absdiff=abs(v_disp - exp_disp),
                        note=("saturated colour" if (exp > sat_hi or exp < sat_lo) else "")))
R = pd.DataFrame(res)
R.to_csv(f"{out}/g4c_{tag}.tsv", sep="\t", index=False)
print(R.groupby("mark").off.agg(["count", "sum"]))
print(R[R.off == True].to_string())
if "note" in R: print("tiles with recomputed |value| beyond colour saturation:\n", R[R.note == "saturated colour"][["row", "mark", "expected", "drawn"]].to_string())
if "colour" in R: print(R.groupby(["mark","colour"]).size().to_string())
