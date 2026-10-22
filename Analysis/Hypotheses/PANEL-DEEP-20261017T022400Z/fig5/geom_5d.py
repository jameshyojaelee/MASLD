"""Figure 5D geometry check: stacked 100% bars -> gene counts; printed labels -> segment counts.
Usage: python geom_5d.py <pdf> <tag>
"""
import sys, re
import numpy as np
import pandas as pd
import pdfplumber

B = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5"
pdf, tag = sys.argv[1], sys.argv[2]
cnt = pd.read_csv(f"{B}/recompute/5d_recomputed_gene_level_counts.tsv", sep="\t").set_index(["cls", "lineage"])
COL = {"#D9D9D9": "under5", "#8EC7E2": "partial", "#7B1FA2": "one", "#00695C": "both"}
LIN = {"Hepatocyte": "hepatocyte", "Stellate": "stellate", "Macrophage": "macrophage", "Cholangiocyte": "cholangiocyte", "T/NK": "t_nk"}

def hexc(c):
    return "#%02X%02X%02X" % tuple(int(round(v * 255)) for v in c) if c else ""

p = pdfplumber.open(pdf).pages[0]
words = p.extract_words(extra_attrs=["size"])
# facets: x-axis lines (horizontal black lines at bottom)
axes = sorted([l for l in p.lines if abs(l["top"] - l["bottom"]) < 0.01 and hexc(l["stroking_color"]) == "#000000" and l["x1"] - l["x0"] > 50], key=lambda l: l["x0"])
ticks = [l for l in p.lines if abs(l["x1"] - l["x0"]) < 0.01 and hexc(l["stroking_color"]) == "#000000"]
tick_labels = [w for w in words if re.fullmatch(r"\d+", w["text"]) and w["top"] > 95]
facet = []
for a in axes:
    tk = sorted([t for t in ticks if a["x0"] - 0.5 <= t["x0"] <= a["x1"] + 0.5], key=lambda t: t["x0"])
    pr = np.array([(t["x0"], float(min(tick_labels, key=lambda w: abs((w["x0"] + w["x1"]) / 2 - t["x0"]))["text"])) for t in tk])
    facet.append({"x0": a["x0"], "x1": a["x1"], "fit": np.polyfit(pr[:, 0], pr[:, 1], 1)})
# facet n from header "(n = 22 genes)"
ns = [int(w["text"]) for w in words if re.fullmatch(r"\d+", w["text"]) and w["top"] < 15]
headers = ["direct", "enzyme"]
rowlab = {LIN[w["text"]]: (w["top"] + w["bottom"]) / 2 for w in words if w["text"] in LIN}
segs = [r for r in p.rects if hexc(r.get("non_stroking_color")) in COL and r["top"] < 95 and r.get("fill")]
# drop duplicate rect drawn twice (fill + stroke)
rows = []
for r in segs:
    fi = 0 if r["x0"] < 120 else 1
    yc = (r["top"] + r["bottom"]) / 2
    lin = min(rowlab, key=lambda k: abs(rowlab[k] - yc))
    f = facet[fi]
    v0, v1 = np.polyval(f["fit"], r["x0"]), np.polyval(f["fit"], r["x1"])
    rows.append((headers[fi], lin, COL[hexc(r["non_stroking_color"])], v0, v1, (v1 - v0)))
sd = pd.DataFrame(rows, columns=["cls", "lineage", "state", "pct0", "pct1", "pct"]).drop_duplicates(["cls", "lineage", "state"])
nmap = dict(zip(headers, ns))
sd["n_drawn"] = sd.apply(lambda r: r.pct / 100 * nmap[r.cls], axis=1)
sd["n_expected"] = sd.apply(lambda r: cnt.loc[(r.cls, r.lineage), r.state], axis=1)
sd["n_total_expected"] = sd.apply(lambda r: cnt.loc[(r.cls, r.lineage), "n_genes"], axis=1)
sd["off"] = (sd.n_drawn - sd.n_expected).abs() > max(0.5, 0)  # tolerance: 0.5 gene (0.5% axis = 0.5 pct)
# zero-count states should be absent
absent = []
for (c, l), r in cnt.iterrows():
    for s in ["under5", "partial", "one", "both"]:
        drawn = ((sd.cls == c) & (sd.lineage == l) & (sd.state == s)).any()
        if r[s] > 0 and not drawn:
            absent.append((c, l, s, r[s]))
        if r[s] == 0 and drawn:
            absent.append((c, l, s, "drawn but expected 0"))
print(tag, "header n", nmap, "expected", cnt.groupby(level=0).n_genes.first().to_dict())
print(tag, "segments", len(sd), "off", int(sd.off.sum()), "max |n diff|", float((sd.n_drawn - sd.n_expected).abs().max()), "absent/extra", absent)
# printed labels inside segments
labs = [w for w in words if re.fullmatch(r"\d+", w["text"]) and 15 < w["top"] < 92 and w["x0"] > 37]
lab_rows = []
for w in labs:
    xc = (w["x0"] + w["x1"]) / 2; yc = (w["top"] + w["bottom"]) / 2
    fi = 0 if xc < 120 else 1
    lin = min(rowlab, key=lambda k: abs(rowlab[k] - yc))
    pct = np.polyval(facet[fi]["fit"], xc)
    hit = sd[(sd.cls == headers[fi]) & (sd.lineage == lin) & (sd.pct0 <= pct) & (sd.pct1 >= pct)]
    st = hit.state.iloc[0] if len(hit) else "none"
    exp = cnt.loc[(headers[fi], lin), st] if st != "none" else None
    lab_rows.append((headers[fi], lin, st, int(w["text"]), exp, exp == int(w["text"])))
ld = pd.DataFrame(lab_rows, columns=["cls", "lineage", "state", "printed", "expected", "ok"])
print(tag, "printed labels", len(ld), "wrong", int((~ld.ok).sum()))
print(ld[~ld.ok].to_string() if (~ld.ok).any() else "")
unl = sd.merge(ld[["cls", "lineage", "state"]], how="left", indicator=True)
unl = unl[unl._merge == "left_only"]
print(tag, "unlabelled segments:", unl[["cls", "lineage", "state", "n_expected", "pct"]].round(1).to_dict("records"))
sd.to_csv(f"{B}/geom_checks/5d_segments_{tag}.tsv", sep="\t", index=False)
ld.to_csv(f"{B}/geom_checks/5d_labels_{tag}.tsv", sep="\t", index=False)
