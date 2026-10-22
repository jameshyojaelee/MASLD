"""Figure 5C geometry check against recomputed values.

Marks: printed NES / q labels (5 processes), right-column dots (25 mRNA + 25 protein; x position and
fill state), pie wedges (125; angular fraction = |rho|, colour sign), heatmap raster cells (25 x 58;
colour inverted through the vector colour bar), top annotation (fibrosis class per column).
Usage: python geom_5c.py <pdf> <tag>
"""
import sys, re
import numpy as np
import pandas as pd
import pdfplumber
from pypdf import PdfReader
from scipy.optimize import linear_sum_assignment

B = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5"
pdf, tag = sys.argv[1], sys.argv[2]
rc = pd.read_csv(f"{B}/recompute/5c_recomputed_right_column.tsv", sep="\t")
ps = pd.read_csv(f"{B}/recompute/5c_recomputed_partial_spearman_raw.tsv", sep="\t")
pe = pd.read_csv(f"{B}/recompute/5c_process_enrichment_copy.tsv", sep="\t")
zt = pd.read_csv(f"{B}/recompute/5c_recomputed_heatmap_z_from_raw.tsv", sep="\t").set_index("gene")
cols = pd.read_csv(f"{B}/recompute/5c_recomputed_heatmap_columns.tsv", sep="\t")
genes = list(rc.gene)
feats = ["Steatosis", "Ballooning", "Inflammation", "Fibrosis", "NAS"]
findings = []

doc = pdfplumber.open(pdf)
p = doc.pages[0]
words = p.extract_words(extra_attrs=["size"])
# gene label rows -> y centre
ylab = {w["text"]: (w["top"] + w["bottom"]) / 2 for w in words if w["text"] in genes and w["x1"] < 116}
missing_labels = [g for g in genes if g not in ylab]
order_ok = [g for g, _ in sorted(ylab.items(), key=lambda kv: kv[1])] == genes
findings.append(("gene labels present and in fixed order", len(missing_labels) == 0 and order_ok, f"missing={missing_labels}"))

# ---- 1. NES / q labels
nes_words = [w for w in words if re.fullmatch(r"[+-]\d\.\d", w["text"]) and w["x0"] < 90]
q_words = [w for w in words if w["text"].startswith("q=") and w["x0"] < 90]
nes_words.sort(key=lambda w: w["top"]); q_words.sort(key=lambda w: w["top"])
rep = pe[pe.representable].sort_values("display_order")
lab_rows = []
for (_, r), wn, wq in zip(rep.iterrows(), nes_words, q_words):
    nes_p = float(wn["text"]); q_p = float(wq["text"][2:])
    ok_n = abs(nes_p - r.NES) <= 0.05 + 1e-9
    # one significant figure
    q_exp = float(f"{r.padj:.0e}")
    ok_q = abs(q_p - q_exp) / q_exp < 1e-6
    lab_rows.append((r.program, wn["text"], round(r.NES, 3), ok_n, wq["text"], f"{r.padj:.3g}", ok_q))
lab = pd.DataFrame(lab_rows, columns=["process", "NES_printed", "NES_recomputed", "NES_ok", "q_printed", "q_recomputed", "q_ok"])
print(lab.to_string())

# ---- 2. right column dots
lines = p.lines
xt = sorted([l for l in lines if abs(l["x1"] - l["x0"]) < 0.01 and 1 < l["bottom"] - l["top"] < 2 and l["x0"] > 370], key=lambda l: l["x0"])
xl = [w for w in words if w["x0"] > 370 and w["top"] > 180 and re.fullmatch(r"-?\d", w["text"])]
pairs = np.array([(t["x0"], float(min(xl, key=lambda w: abs((w["x0"] + w["x1"]) / 2 - t["x0"]))["text"])) for t in xt])
ax = np.polyfit(pairs[:, 0], pairs[:, 1], 1)
xrng = np.polyval(ax, 467.535) - np.polyval(ax, 374.805)
def hexc(c):
    return "#%02X%02X%02X" % tuple(int(round(v * 255)) for v in c) if c else ""
dots = [c for c in p.curves if c["x0"] > 370 and 30 < c["top"] < 175]
dot_rows = []
for c in dots:
    cx = (c["x0"] + c["x1"]) / 2; cy = c["top"]
    stroke = hexc(c.get("stroking_color")); fill = hexc(c.get("non_stroking_color"))
    g = min(ylab, key=lambda k: abs(ylab[k] - cy))
    kind = "mrna" if stroke == "#8EC7E2" else ("protein" if stroke == "#E8853A" else "other")
    dot_rows.append((g, kind, np.polyval(ax, cx), fill != "#FFFFFF", stroke))
dd = pd.DataFrame(dot_rows, columns=["gene", "kind", "drawn", "filled", "stroke"])
print(tag, "non-data (legend) dots excluded:", int((dd.kind == "other").sum()))
dd = dd[dd.kind != "other"]
dd = dd.merge(rc[["gene", "mrna_log2fc", "mrna_filled", "protein_log2fc", "protein_filled"]], on="gene")
dd["expected"] = np.where(dd.kind == "mrna", dd.mrna_log2fc, dd.protein_log2fc)
dd["exp_filled"] = np.where(dd.kind == "mrna", dd.mrna_filled, dd.protein_filled)
dd["pos_ok"] = (dd.drawn - dd.expected).abs() <= 0.005 * xrng + 0.01
dd["fill_ok"] = dd.filled == dd.exp_filled
print(tag, "right-column dots", len(dd), "pos off", (~dd.pos_ok).sum(), "fill off", (~dd.fill_ok).sum(), "max |diff|", (dd.drawn - dd.expected).abs().max())
xlim = (np.polyval(ax, 374.805), np.polyval(ax, 467.535))
print(tag, "right column x range", xlim)

# ---- 3. pies
circles = [c for c in p.curves if 300 < c["x0"] < 370 and c["top"] > 30 and hexc(c.get("non_stroking_color")) == "#FFFFFF"]
wedges = [c for c in p.curves if 300 < c["x0"] < 370 and c["top"] > 30 and hexc(c.get("non_stroking_color")) != "#FFFFFF"]
fx = sorted(set(round((c["x0"] + c["x1"]) / 2, 1) for c in circles))
pie_rows = []
for w in wedges:
    pts = np.array(w["pts"])
    centre = pts[0]
    arc = pts[1:-1]
    ang = np.unwrap(np.arctan2(arc[:, 0] - centre[0], -(arc[:, 1] - centre[1])))  # clockwise from 12 o'clock (screen)
    frac = abs(ang[-1] - ang[0]) / (2 * np.pi)
    g = min(ylab, key=lambda k: abs(ylab[k] - centre[1]))
    f = feats[int(np.argmin([abs(centre[0] - x) for x in fx]))]
    col = w["non_stroking_color"]
    sign = 1 if col[0] > col[2] else -1
    pie_rows.append((g, f, frac, sign, hexc(col)))
pp = pd.DataFrame(pie_rows, columns=["gene", "feature", "frac", "sign", "fill"])
pp = ps.merge(pp, on=["gene", "feature"], how="left")
pp["abs_rho"] = pp.rho.abs()
pp["frac_ok"] = (pp.frac - pp.abs_rho).abs() <= 0.02
pp["sign_ok"] = pp.sign == np.sign(pp.rho)
print(tag, "pies expected", len(pp), "wedges drawn", len(wedges), "circles", len(circles), "missing wedge", pp.frac.isna().sum(),
      "frac off", (~pp.frac_ok & pp.frac.notna()).sum(), "sign off", (~pp.sign_ok & pp.frac.notna()).sum(),
      "max |frac-|rho||", (pp.frac - pp.abs_rho).abs().max())

# ---- 4. heatmap raster
img_meta = sorted(p.images, key=lambda im: im["name"])
reader = PdfReader(pdf)
ims = {im.name.strip("/").split(".")[0]: im.image for im in reader.pages[0].images}
# heatmap = layer with the larger opaque area
arrs = {k: np.asarray(v.convert("RGBA")).astype(int) for k, v in ims.items()}
hk = max(arrs, key=lambda k: (arrs[k][:, :, 3] > 0).sum())
ak = [k for k in arrs if k != hk][0]
im0 = [im for im in p.images][0]
W = arrs[hk].shape[1]; H = arrs[hk].shape[0]
sx = (im0["x1"] - im0["x0"]) / W; sy = (im0["bottom"] - im0["top"]) / H
A = arrs[hk]
ys, xs = np.where(A[:, :, 3] > 0)
x0p, x1p, y0p, y1p = xs.min(), xs.max() + 1, ys.min(), ys.max() + 1
ncol, nrow = 58, 25
cw = (x1p - x0p) / ncol; rh = (y1p - y0p) / nrow
# colour bar (vector rects, width 0.8) left one = z, right one = rho
cb = [r for r in p.rects if round(r["width"], 1) == 0.8]
cbz = sorted([r for r in cb if r["x0"] < 250], key=lambda r: r["x0"])
zt_l = [w for w in words if 160 < w["x0"] < 230 and 180 < w["top"] < 188 and re.fullmatch(r"-?\d", w["text"])]
zpairs = np.array([((w["x0"] + w["x1"]) / 2, float(w["text"])) for w in zt_l])
azb = np.polyfit(zpairs[:, 0], zpairs[:, 1], 1)
cb_val = np.array([np.polyval(azb, (r["x0"] + r["x1"]) / 2) for r in cbz])
cb_rgb = np.array([np.array(r["non_stroking_color"]) * 255 for r in cbz])
drawn = np.full((nrow, ncol), np.nan); isgray = np.zeros((nrow, ncol), bool)
for i in range(nrow):
    for j in range(ncol):
        px = A[int(y0p + (i + 0.5) * rh), int(x0p + (j + 0.5) * cw), :3]
        if abs(px[0] - px[1]) < 6 and abs(px[1] - px[2]) < 6 and px[0] < 235:
            isgray[i, j] = True
            continue
        d = np.abs(cb_rgb - px).sum(1)
        k = np.argmin(d)
        drawn[i, j] = cb_val[k]
exp = zt.loc[genes].values.astype(float)
exp_c = np.clip(exp, -2, 2)
# column matching (ties within fibrosis/NAS may be ordered differently)
cost = np.zeros((ncol, ncol))
for a in range(ncol):
    for b in range(ncol):
        m = ~np.isnan(drawn[:, a]) & ~np.isnan(exp_c[:, b])
        cost[a, b] = np.abs(drawn[m, a] - exp_c[m, b]).mean() + 5 * (np.isnan(drawn[:, a]) != np.isnan(exp_c[:, b])).sum()
ri, ci = linear_sum_assignment(cost)
same_pos = (ci == np.arange(ncol)).sum()
diff = np.abs(drawn[:, ri] - exp_c[:, ci])
na_mismatch = (np.isnan(drawn[:, ri]) != np.isnan(exp_c[:, ci])).sum()
cell_off = int(np.nansum(diff > 0.15))
print(tag, "heatmap cells", nrow * ncol, "gray(NA) drawn", int(isgray.sum()), "NA expected", int(np.isnan(exp).sum()),
      "NA mismatch", int(na_mismatch), "cells off >0.15 z", cell_off, "median |diff|", float(np.nanmedian(diff)),
      "columns matched to same position", same_pos, "cells beyond |2| (saturated)", int((np.abs(exp) > 2).sum()))
# groups of matched columns: do matched participants share fibrosis/NAS with column position?
fib_pos = cols.Fibrosis.values; nas_pos = cols.NAS.values
fib_match = (fib_pos[ri] == fib_pos[ci]).mean(); nas_match = (nas_pos[ri] == nas_pos[ci]).mean()
print(tag, "matched column has same fibrosis as position", fib_match, "same NAS", nas_match)

# ---- 5. top annotation: fibrosis classes
Aa = arrs[ak]
ys2, xs2 = np.where(Aa[:, :, 3] > 0)
yb0, yb1 = ys2.min(), ys2.max() + 1
xb0, xb1 = xs2.min(), xs2.max() + 1
rowh = (yb1 - yb0) / 2
cwa = (xb1 - xb0) / ncol
fib_cols = [tuple(Aa[int(yb0 + 0.5 * rowh), int(xb0 + (j + 0.5) * cwa), :3]) for j in range(ncol)]
nas_cols = [tuple(Aa[int(yb0 + 1.5 * rowh), int(xb0 + (j + 0.5) * cwa), :3]) for j in range(ncol)]
# class runs
def runs(seq):
    out = []; prev = None
    for s in seq:
        if s != prev:
            out.append([s, 1]); prev = s
        else:
            out[-1][1] += 1
    return out
fr = runs(fib_cols)
print(tag, "fibrosis annotation runs (drawn)", [n for _, n in fr], "expected", [int((fib_pos == k).sum()) for k in sorted(set(fib_pos))])
nr = runs(nas_cols)
print(tag, "NAS annotation runs (drawn)", [n for _, n in nr])
exp_nas_runs = [n for _, n in runs(list(nas_pos))]
print(tag, "NAS runs expected", exp_nas_runs)

out = {"labels": lab, "dots": dd, "pies": pp}
lab.to_csv(f"{B}/geom_checks/5c_labels_{tag}.tsv", sep="\t", index=False)
dd.to_csv(f"{B}/geom_checks/5c_dots_{tag}.tsv", sep="\t", index=False)
pp.to_csv(f"{B}/geom_checks/5c_pies_{tag}.tsv", sep="\t", index=False)
pd.DataFrame(drawn, index=genes).to_csv(f"{B}/geom_checks/5c_heatmap_drawn_{tag}.tsv", sep="\t")
print(tag, "SUMMARY", findings)
