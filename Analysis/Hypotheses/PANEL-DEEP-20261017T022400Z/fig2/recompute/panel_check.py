"""Label and geometry checks for Fig 2 A, C, D, D-ancestry, D-subclass, E, G/H text, I against the recomputed values.
usage: panel_check.py <dump dir> <tag>   (dump dir holds <panel>.{words,rects,lines,curves}.tsv)"""
import sys, json, re, numpy as np, pandas as pd
D, tag = sys.argv[1], sys.argv[2]
rd = lambda p, k: pd.read_csv(f"{D}/{p}.{k}.tsv", sep="\t")
A_loci = json.load(open("out_ABCEF/A_loci.json"))
dm = pd.read_csv("out_D/D_main_fractions.tsv", sep="\t").set_index("class")
reg = pd.read_csv("upstream/registry_t12.tsv", sep="\t"); reg["a"] = reg.ancestry.replace({"CSA": "SAS"})
res = {}
def words(p): return rd(p, "words").text.astype(str).tolist()
def phrase(p): return " ".join(words(p))
def need(p, toks):
    w = phrase(p); return [t for t in toks if t not in w]
ANC = {"#4C72B0": "EUR", "#55A868": "AFR", "#C44E52": "EAS", "#DD8452": "AMR", "#8172B3": "SAS"}
# ---------- A ----------
p = "fig2A_gwas_cascade"; r = rd(p, "rects"); out = {}
r = r[(r.h < 100) & (r.top < 155)]
cols = sorted(r.x0.round(1).unique())
c1 = r[r.x0.round(1) == cols[0]]; c2 = r[r.x0.round(1) == cols[1]]; c3 = r[r.x0.round(1) == cols[2]]; c4 = r[r.x0.round(1) == cols[3]]
u1 = c1.h.sum() / 35
out["col1_counts"] = {ANC.get(f, f): round(h / u1, 2) for f, h in zip(c1.fill, c1.h)}
exp1 = reg.a.value_counts().to_dict()
out["col1_off"] = [k for k, v in out["col1_counts"].items() if abs(v - exp1.get(k, 0)) > 0.05]
u2 = c2.h.min()  # one GWAS
out["col2_total_gwas"] = round(c2.h.sum() / u2, 2)
# trait blocks in column 3 are separated by gaps; assign by order top->bottom using expected trait order from words
u3 = c3.h.sum() / 360; out["col3_unit_pt_per_locus"] = u3
c3 = c3.sort_values("top").copy(); c3["n"] = (c3.h / u3).round(2); c3["anc"] = c3.fill.map(ANC)
gap = c3.top.values[1:] - c3.bottom.values[:-1]
blk = np.r_[0, np.cumsum(gap > 1.0)]; c3["blk"] = blk
order = ["NAFLD", "NASH", "PDFF", "ALT", "AST", "GGT"]
c3["trait"] = [order[b] if b < len(order) else "?" for b in c3.blk]
got3 = {f"{a}|{t}": n for a, t, n in zip(c3.anc, c3.trait, c3.n)}
exp3 = json.load(open("out_ABCEF/A_loci_tuple_sort.json"))  # drawn allocation: max PIP, ties broken by tuple sort (a_loci_ties.py)
out["col3_off"] = {k: (got3.get(k, 0), exp3.get(k, 0)) for k in set(got3) | set(exp3) if abs(got3.get(k, 0) - exp3.get(k, 0)) > 0.5}
out["col3_trait_totals"] = c3.groupby("trait").n.sum().round(1).to_dict()
u4 = c4.h.sum() / 996; c4 = c4.sort_values("top")
out["col4_counts"] = (c4.h / u4).round(1).tolist()
out["col4_off"] = [(g, e) for g, e in zip(out["col4_counts"], [244, 89, 103, 67, 89, 404]) if abs(g - e) > 0.5]
out["missing_labels"] = need(p, ["EUR (17)", "AFR (6)", "EAS (6)", "AMR (3)", "SAS (3)", "NAFLD (11)", "NASH (1)", "PDFF (3)", "ALT (8)", "AST (8)", "GGT (4)",
                                 "122", "115", "93", "union 996", "244", "89", "103", "67", "404", "436", "560", "5 ancestries", "35 GWAS", "Fine-mapped", "COLOC"])
ch = rd(p, "chars"); W = 254.16
out["chars_beyond_page"] = ch[ch.x1 > W + 0.05].text.tolist()
out["font_sizes"] = sorted(ch["size"].unique().tolist())
res["A"] = out
# ---------- C ----------
p = "Fig2C_coding_noncoding_split"; r = rd(p, "rects"); r = r[r.w < 30].sort_values("top")
tot = r.h.sum(); cnt = (r.h / tot * 996).round(1).tolist()
exp = [14, 21, 70, 104, 450, 337]
res["C"] = {"segment_counts_if_total_996": cnt, "fills": r.fill.tolist(),
            "off": [(g, e) for g, e in zip(cnt, exp) if abs(g - e) > 0.5],
            "missing_labels": need(p, ["996", "(14)", "(21)", "(70)", "(104)", "(450)", "(337)", "coding, this gene", "coding, other gene", "UTR", "promoter", "intron", "intergenic"]),
            "chars_beyond_page": rd(p, "chars").query("x1 > 93.05").text.tolist()}
# ---------- D main ----------
p = "Fig2D_pip_architecture_by_trait_directness"; r = rd(p, "rects"); ln = rd(p, "lines"); wd = rd(p, "words")
yt = ln[(ln.h < 0.01) & ln.w.between(2, 4) & (ln.stroke == "#000000")].top.sort_values(ascending=False).values
yl = wd[wd.text.astype(str).isin(["0", "50", "100"])].copy(); yl["cy"] = (yl.top + yl.bottom) / 2; yl = yl.sort_values("cy", ascending=False)
ay = np.polyfit(yt[:len(yl)], yl.text.astype(float).values, 1)
bars = r[(r.w > 10) & (r.w < 20)].copy(); bars["v"] = bars.h * abs(ay[0])
FILL = {"#C9265E": "protein_altering_pip_mass", "#8E3B76": "canonical_splice_pip_mass", "#4DB6AC": "synonymous_or_utr_pip_mass", "#B0BEC5": "other_noncoding_pip_mass"}
xs = sorted(bars.x0.round(1).unique()); cls = {xs[0]: "direct", xs[1]: "enzyme"}
offs = []; got = {}
for _, b in bars.iterrows():
    c = cls[round(b.x0, 1)]; k = FILL[b.fill]; e = dm.loc[c, k]; got[f"{c}|{k}"] = round(b.v, 3)
    if abs(b.v - e) > 0.5: offs.append((c, k, round(b.v, 3), round(e, 3)))
for c in ["direct", "enzyme"]:
    for k in FILL.values():
        if f"{c}|{k}" not in got and dm.loc[c, k] > 0.5: offs.append((c, k, "absent", round(dm.loc[c, k], 3)))
res["D"] = {"segments": got, "off": offs, "y_axis_title_chars": "".join(rd(p, "chars").query("x0 < 8").sort_values("top", ascending=False).text),
            "missing_labels": need(p, ["23.1", "3.6", "73.5", "3.9", "92.5", "Direct", "Liver", "enzymes", "Protein-alt.", "Splice", "Noncoding"]),
            "n_printed": bool(re.search(r"\b(207|6,?208|233|19|228)\b", phrase(p)))}
# ---------- D ancestry / subclass ----------
def forest(p, rows_exp, ticks):
    ln = rd(p, "lines"); cv = rd(p, "curves"); wd = rd(p, "words")
    xt = ln[(ln.w < 0.01) & ln.h.between(1, 2)].x0.sort_values().values
    lab = wd[wd.text.astype(str).isin([str(t) for t in ticks]) & (wd.top > ln[(ln.w < 0.01) & ln.h.between(1, 2)].top.min())].copy()
    lab["cx"] = (lab.x0 + lab.x1) / 2; lab = lab.sort_values("cx")
    ax = np.polyfit(xt[:len(lab)], lab.text.astype(float).values, 1); rng = ticks[-1] - ticks[0]
    pts = cv.copy(); pts["v"] = np.polyval(ax, (pts.x0 + pts.x1) / 2); pts = pts.sort_values("top")
    thick = ln[(ln.lw > 1) & (ln.h < 0.01)].sort_values("top"); thin = ln[(ln.stroke == "#999999") & (ln.h < 0.01)].sort_values("top")
    out = {"rows": []}; off = []
    for i, (name, pe, lo, hi, lo2, hi2) in enumerate(rows_exp):
        if i >= len(pts): off.append((name, "missing row")); continue
        y = pts.top.iloc[i]; row = {"row": name, "point": round(pts.v.iloc[i], 2)}
        th = thick[(thick.top - y).abs() < 0.5]; tn = thin[(thin.top - y).abs() < 0.5]
        if len(th): row["ci95"] = (round(np.polyval(ax, th.x0.iloc[0]), 2), round(np.polyval(ax, th.x1.iloc[0]), 2))
        if len(tn): row["ci98"] = (round(np.polyval(ax, tn.x0.iloc[0]), 2), round(np.polyval(ax, tn.x1.iloc[0]), 2))
        tol = 0.005 * rng
        if abs(row["point"] - pe) > tol: off.append((name, "point", row["point"], pe))
        if "ci95" in row and (abs(row["ci95"][0] - lo) > tol or abs(row["ci95"][1] - hi) > tol): off.append((name, "ci95", row["ci95"], (lo, hi)))
        if lo2 is not None and "ci98" not in row: off.append((name, "ci98 missing"))
        if lo2 is not None and "ci98" in row and (abs(row["ci98"][0] - lo2) > tol or abs(row["ci98"][1] - hi2) > tol): off.append((name, "ci98", row["ci98"], (lo2, hi2)))
        out["rows"].append(row)
    out["n_points"] = len(pts); out["off"] = off; return out
rows = [("pooled", 19.4695, 6.94, 34.01, None, None), ("EUR+AMR", 14.9727, 3.879, 27.180, 1.334, 30.272), ("EUR", 22.961, 5.498, 42.365, 1.789, 47.034),
        ("AMR", 2.1915, -4.891, 10.655, -6.329, 13.510), ("nonEUR", 3.3098, -3.830, 12.745, None, None)]
p = "Fig2D_coding_by_ancestry"; o = forest(p, rows, [-10, 0, 10, 20, 30, 40, 50])
o["missing_labels"] = need(p, ["23.3% vs 3.9%; 19.5 (6.9 to 34.0)", "21.0% vs 6.0%; 15.0 (3.9 to 27.2)", "29.9% vs 6.9%; 23.0 (5.5 to 42.4)",
                               "6.7% vs 4.5%; 2.2 (-4.9 to 10.7)", "7.0% vs 3.7%; 3.3 (-3.8 to 12.7)", "percentage points"])
res["D_ancestry"] = o
rows = [("all direct", 19.4695, 6.95, 33.87, None, None), ("PDFF", 36.2646, 18.35, 54.55, None, None), ("NAFLD+NASH", 8.5667, -1.64, 20.96, None, None)]
p = "Fig2D_direct_subclass_split"; o = forest(p, rows, [0, 10, 20, 30, 40, 50])
o["missing_labels"] = need(p, ["23.3% vs 3.9%; 207 sets / 19 clusters", "40.1% vs 3.9%; 37 sets / 12 clusters", "12.4% vs 3.9%; 170 sets / 13 clusters", "percentage points"])
res["D_subclass"] = o
# ---------- E ----------
p = "Fig2E_ancestry_unique_coloc_GWS"; r = rd(p, "rects"); r = r[r.w < 140].copy(); r["area"] = r.w * r.h
tot = r.area.sum(); EF = {"#D6DBDE": ("EUR-only", 340), "#A7B6BE": ("shared", 64), "#E9ECEE": ("sub-threshold non-EUR", 23), "#F79268": ("GWS non-EUR", 9)}
e = {EF.get(f, (f, None))[0]: round(a / tot * 436, 2) for f, a in zip(r.fill, r.area)}
res["E"] = {"areas_as_genes_if_total_436": e, "off": [(k, v, ex) for f, (k, ex) in EF.items() for kk, v in e.items() if kk == k and abs(v - ex) > 0.5] + [(k, "absent") for f, (k, ex) in EF.items() if k not in e],
            "missing_labels": need(p, ["340", "(78%)", "64", "(15%)", "23", "(5%)", "9", "(2%)", "8 EAS, 1 AFR"])}
# ---------- G/H text ----------
for p, pp, s4 in [("Fig2G_RORA_NAFLD", "0.38", "0.98"), ("Fig2H_FABP1_ALT", "0.56", "1.00")]:
    res[p[:5]] = {"missing_labels": need(p, [f"SNP.PP.H4 = {pp}", f"PP4 = {s4}", "GWAS", "SuSiE", "eQTL", "Genes"]),
                  "gene_text_fills": rd(p, "chars").query("font.str.contains('Oblique')", engine="python").groupby("fill").text.agg("".join).to_dict()}
# ---------- I ----------
p = "Fig2I_locus_legend"; wd = rd(p, "words"); r = rd(p, "rects"); cv = rd(p, "curves")
ns = wd[wd.text == "n.s."]; sw = r[r.fill == "#9E9E9E"]
res["I"] = {"missing_labels": need(p, ["r² to lead", "Lead SNP", "EUR (in CS)", "not in CS", "coloc shared variant", "log2FC", "+0.3", "n.s."]),
            "ns_swatch_overlap_pt": float(max(0, sw.x1.max() - ns.x0.min())) if len(ns) and len(sw) else None,
            "marker_fills": cv.fill.tolist(), "eqtl_marker_keyed": bool((cv.fill == "#E07B39").any())}
json.dump(res, open(f"out_geom/panel_check_{tag}.json", "w"), indent=1, default=str)
print(json.dumps(res, indent=1, default=str))
