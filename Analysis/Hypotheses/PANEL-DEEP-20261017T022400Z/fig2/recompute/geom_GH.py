"""Fig 2G/2H geometry: recover GWAS -log10 p, SuSiE PIP and eQTL -log10 p of every drawn point and match to
recomputed upstream values (GWAS sumstats, fine-mapping tables, Broadaway eQTL marginals).
usage: geom_GH.py <dump prefix> <G|H> <tag>"""
import sys, json, numpy as np, pandas as pd
from scipy.spatial import cKDTree
pre, panel, tag = sys.argv[1], sys.argv[2], sys.argv[3]
cv = pd.read_csv(f"{pre}.curves.tsv", sep="\t"); ln = pd.read_csv(f"{pre}.lines.tsv", sep="\t"); wd = pd.read_csv(f"{pre}.words.tsv", sep="\t")
U = "upstream"
if panel == "G":
    gw = pd.read_csv(f"{U}/gwas_window_MVP_NAFLD_EUR_chr15.tsv", sep="\t"); eq = pd.read_csv(f"{U}/eqtl_RORA.tsv", sep="\t")
    coloc_snp, gwas_hit1 = 60878030, 60899031
else:
    gw = pd.read_csv(f"{U}/gwas_window_UKBB_ALT_chr2.tsv", sep="\t"); eq = pd.read_csv(f"{U}/eqtl_FABP1.tsv", sep="\t")
    coloc_snp, gwas_hit1 = 88424066, 88423510
gw["lp"] = -np.log10(gw.pval); eq["lp"] = -np.log10(eq.PVAL)
gwm = gw.groupby("position").lp.max().reset_index()
eqm = eq.groupby("POS").lp.max().reset_index().rename(columns={"POS": "position"})
# y transforms per track from tick marks (x0 = 45, w = 3.602) and their numeric labels (x1 < 44)
tk = ln[(ln.x0.round(1) == 45.0) & (ln.w.round(2) == 3.6)].top.values
lab = wd[(wd.x1 < 44) & (wd.x0 > 30) & wd.text.astype(str).str.match(r"^[0-9.]+$")].copy(); lab["cy"] = (lab.top + lab.bottom) / 2
lab["tick"] = lab.cy.apply(lambda y: tk[np.abs(tk - y).argmin()]); lab["v"] = lab.text.astype(float)
def track(lo, hi):
    l = lab[(lab.tick >= lo) & (lab.tick <= hi)]; return np.polyfit(l.tick, l.v, 1), l
(ag, lg), (ap, lpp), (ae, le) = track(0, 70), track(70, 95), track(95, 150)
yv = lambda a, y: np.polyval(a, y)
cv["cx"] = (cv.x0 + cv.x1) / 2
dia = cv[cv.npts == 5].copy(); dia["cy"] = (dia.top + dia.bottom) / 2
lead_d = dia[dia.fill == "#C2185B"].iloc[0]; gold = dia[dia.fill == "#FFD600"].sort_values("cy")
pts = cv[cv.npts == 3].copy(); pts["cy"] = pts.top
gwas_pts = pts[(pts.cy < 63) & (pts.fill != "#FFFFFF")].copy(); gwas_pts["v"] = yv(ag, gwas_pts.cy)
sus_pts = pts[(pts.cy > 70) & (pts.cy < 90)].copy(); sus_pts["v"] = yv(ap, sus_pts.cy)
eq_pts = pts[(pts.fill == "#E07B39") & (pts.cy > 95)].copy(); eq_pts["v"] = yv(ae, eq_pts.cy)
# x transform: fit scale/offset by matching GWAS points; anchor first on the gold diamond = coloc top SNP
best = None
for scale in np.linspace(0.25e-3, 0.5e-3, 2501):
    off = gold.cx.iloc[0] - scale * coloc_snp
    x = off + scale * gwm.position.values
    t = cKDTree(np.c_[x, gwm.lp.values * 5])
    d, _ = t.query(np.c_[gwas_pts.cx.values, gwas_pts.v.values * 5])
    n = int((d < 0.3).sum())
    if best is None or n > best[0]: best = (n, scale, off)
n, scale, off = best
pos_of = lambda x: (x - off) / scale
def match(drawn, ref, tolx_bp, toly):
    t = cKDTree(np.c_[ref.position.values / tolx_bp, ref.lp.values / toly])
    d, i = t.query(np.c_[pos_of(drawn.cx.values) / tolx_bp, drawn.v.values / toly])
    return d, i
dg, ig = match(gwas_pts, gwm, 0.3 / scale, 0.005 * (lg.v.max() - lg.v.min()))
de, ie = match(eq_pts, eqm, 0.3 / scale, 0.005 * (le.v.max() - le.v.min()))
out = {"panel": panel, "tag": tag, "x_scale_pt_per_bp": scale, "x_window_bp": [pos_of(48.602), pos_of(227.16)],
       "gwas_axis": lg[["v", "tick"]].values.tolist(), "pip_axis": lpp[["v", "tick"]].values.tolist(), "eqtl_axis": le[["v", "tick"]].values.tolist(),
       "gws_line_value": None, "n_gwas_marks": len(gwas_pts), "gwas_matched_within_tol": int((dg <= 1).sum()),
       "n_eqtl_marks": len(eq_pts), "eqtl_matched_within_tol": int((de <= 1).sum()),
       "gwas_variants_in_window_recomputed": int(((gwm.position >= pos_of(48.602)) & (gwm.position <= pos_of(227.16))).sum()),
       "eqtl_variants_in_window_recomputed": int(((eqm.position >= pos_of(48.602)) & (eqm.position <= pos_of(227.16))).sum()),
       "lead_diamond": {"pos": pos_of(lead_d.cx), "lp_drawn": float(yv(ag, lead_d.cy)),
                        "gwas_hit1": gwas_hit1, "lp_hit1": float(gwm.set_index("position").lp.get(gwas_hit1, np.nan)),
                        "lp_coloc_snp": float(gwm.set_index("position").lp.get(coloc_snp, np.nan)),
                        "window_min_p_pos": int(gwm.loc[gwm.lp.idxmax(), "position"]), "window_max_lp": float(gwm.lp.max())},
       "gold_diamonds": [{"pos": pos_of(r.cx), "y_value_on_track_axis": float(yv(ap if r.cy < 95 else ae, r.cy)), "track": "PIP" if r.cy < 95 else "eQTL"} for _, r in gold.iterrows()],
       "eqtl_at_coloc_snp": float(eqm.set_index("position").lp.get(coloc_snp, np.nan)),
       "susie_open": int((sus_pts.fill == "#FFFFFF").sum()), "susie_filled": int((sus_pts.fill == "#1565C0").sum()),
       "susie_filled_values": sorted([round(float(v), 4) for v in sus_pts[sus_pts.fill == "#1565C0"].v], reverse=True),
       "susie_max_open": float(sus_pts[sus_pts.fill == "#FFFFFF"].v.max()),
       "black_points": [(round(float(pos_of(r.cx))), round(float(yv(ag, r.cy)), 3)) for _, r in pts[(pts.fill == "#000000")].iterrows()]}
gl = ln[(ln.stroke == "#000000") & (ln.w > 170) & (ln.top < 62) & (ln.top > 8)]
if len(gl): out["gws_line_value"] = float(yv(ag, gl.top.iloc[0]))
ref = ln[(ln.stroke == "#7F7F7F")]
if len(ref): out["pip_reference_line"] = float(yv(ap, ref.top.iloc[0]))
# SuSiE track vs available fine-mapping for this window
fm = pd.concat([pd.read_csv(f"{U}/combined_finemapping_GH_windows.tsv.gz", sep="\t"), pd.read_csv(f"{U}/combined_finemapping_coloc_targeted_GH_windows.tsv.gz", sep="\t")])
chrom = 15 if panel == "G" else 2
fm = fm[fm.chromosome == chrom]
sp = sus_pts.assign(pos=pos_of(sus_pts.cx))
res = {}
for (s, l), f in fm.groupby(["study", "locus"]):
    f = f.groupby("position").susie_pip.max().reset_index().rename(columns={"susie_pip": "lp"})
    t = cKDTree(np.c_[f.position.values * scale / 0.3, f.lp.values / 0.005])
    d, i = t.query(np.c_[sp.pos.values * scale / 0.3, sp.v.values / 0.005])
    res[f"{s}|{l}"] = {"n_fm": len(f), "drawn_matched": int((d <= 1).sum()), "of": len(sp), "fm_max_pip": float(f.lp.max())}
out["susie_track_vs_finemapping"] = res
json.dump(out, open(f"out_geom/{panel}_{tag}.json", "w"), indent=1, default=str)
print(json.dumps(out, indent=1, default=str))
