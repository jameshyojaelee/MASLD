"""Figure 5E geometry check: row set, bar lengths/shades, interval endpoints, point positions, colours.
Compared with recompute/5e_recomputed_117_strict.tsv (panel membership rule) and
recompute/5e_recomputed_117.tsv (ALT-copy-inclusive membership).
Usage: python geom_5e.py <pdf> <tag>
"""
import sys, re
import numpy as np
import pandas as pd
import pdfplumber

B = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5"
pdf, tag = sys.argv[1], sys.argv[2]
rs = pd.read_csv(f"{B}/recompute/5e_recomputed_117_strict.tsv", sep="\t")
ra = pd.read_csv(f"{B}/recompute/5e_recomputed_117.tsv", sep="\t")

def hexc(c):
    return "#%02X%02X%02X" % tuple(int(round(v * 255)) for v in c) if c else ""

p = pdfplumber.open(pdf).pages[0]
words = p.extract_words(extra_attrs=["size"])
lines = p.lines
# group separators: long horizontal grey lines in the label column
seps = sorted(set(round(l["top"], 2) for l in lines if hexc(l["stroking_color"]) == "#9A9A9A" and l["x0"] < 10))
groups = ["hepatocytes", "cholangiocytes", "fibroblasts", "macrophages", "tcells"]
def group_of(y):
    return groups[sum(y > s for s in seps)]
# axes
def axis_fit(xmin, xmax, labels_re):
    tk = sorted([l for l in lines if abs(l["x1"] - l["x0"]) < 0.01 and hexc(l["stroking_color"]) == "#222222" and xmin <= l["x0"] <= xmax and l["top"] > 370], key=lambda l: l["x0"])
    lab = [w for w in words if w["top"] > 378 and xmin - 6 <= (w["x0"] + w["x1"]) / 2 <= xmax + 6 and re.fullmatch(labels_re, w["text"])]
    pr = np.array([(t["x0"], float(min(lab, key=lambda w: abs((w["x0"] + w["x1"]) / 2 - t["x0"]))["text"])) for t in tk])
    return np.polyfit(pr[:, 0], pr[:, 1], 1), pr
fb, prb = axis_fit(120, 179, r"\d")
fr, prr = axis_fit(181, 260, r"-?\d(\.\d)?")
print(tag, "bar axis ticks", prb.tolist(), "rho axis ticks", prr.tolist())
rho_rng = np.polyval(fr, 259.523) - np.polyval(fr, 180.555)
# rows from background tracks
tracks = sorted([r for r in p.rects if hexc(r.get("non_stroking_color")) == "#F1F1F1"], key=lambda r: r["top"])
rows = []
for t in tracks:
    yc = (t["top"] + t["bottom"]) / 2
    lw = sorted([w for w in words if 25 < w["x0"] < 120 and abs((w["top"] + w["bottom"]) / 2 - yc) < 2.5], key=lambda w: w["x0"])
    label = " ".join(w["text"] for w in lw)
    mod = int(lw[0]["text"]) if lw and lw[0]["text"].isdigit() else None
    bars = [r for r in p.rects if hexc(r.get("non_stroking_color")) in ("#C79A3E", "#7A5A1E") and abs((r["top"] + r["bottom"]) / 2 - yc) < 2]
    nbar = {"#7A5A1E": 0.0, "#C79A3E": 0.0}
    for b in bars:
        nbar[hexc(b["non_stroking_color"])] += np.polyval(fb, b["x1"]) - np.polyval(fb, b["x0"])
    ci = [l for l in lines if abs(l["top"] - yc) < 1.5 and abs(l["top"] - l["bottom"]) < 0.01 and l["x0"] > 180 and abs(l["linewidth"] - 1.024299) < 0.01]
    pt = [l for l in lines if abs((l["top"] + l["bottom"]) / 2 - yc) < 1.5 and abs(l["x1"] - l["x0"]) < 0.01 and l["x0"] > 180 and abs(l["linewidth"] - 1.387072) < 0.01]
    nots = [w for w in words if w["text"] == "not" and abs((w["top"] + w["bottom"]) / 2 - yc) < 2.5]
    rows.append({"y": yc, "group": group_of(yc), "module": mod, "label": label,
                 "bar_direct": nbar["#7A5A1E"], "bar_enzyme": nbar["#C79A3E"],
                 "ci_low_d": np.polyval(fr, ci[0]["x0"]) if ci else np.nan, "ci_high_d": np.polyval(fr, ci[0]["x1"]) if ci else np.nan,
                 "rho_d": np.polyval(fr, pt[0]["x0"]) if pt else np.nan,
                 "colour": hexc(ci[0]["stroking_color"]) if ci else "", "not_scored_label": bool(nots)})
d = pd.DataFrame(rows)
print(tag, "rows drawn", len(d), d.group.value_counts().to_dict())
def compare(rec, name):
    disp = rec[rec.display]
    m = d.merge(rec, left_on=["group", "module"], right_on=["cell_type", "module"], how="left", indicator=True)
    missing_from_panel = disp.merge(d, left_on=["cell_type", "module"], right_on=["group", "module"], how="left", indicator=True)
    missing_from_panel = missing_from_panel[missing_from_panel._merge == "left_only"]
    extra = m[~m.display.fillna(False).astype(bool)]
    m["bar_total_d"] = m.bar_direct + m.bar_enzyme
    m["bar_ok"] = (m.bar_total_d - m.n_open).abs() < 0.05
    m["shade_ok"] = ((m.bar_direct > 0.5) == (m.n_open_direct > 0))
    sc = m.testable.astype(str).str.upper() == "TRUE"
    tol = 0.005 * rho_rng + 0.005
    m["rho_ok"] = ~sc | ((m.rho_d - m.rho).abs() <= tol)
    m["ci_ok"] = ~sc | (((m.ci_low_d - m.ci_low).abs() <= tol) & ((m.ci_high_d - m.ci_high).abs() <= tol))
    m["colour_exp"] = np.where(~sc, "", np.where(m.robust_re & (m.rho > 0), "#B2182B", np.where(m.robust_re & (m.rho < 0), "#2166AC", "#C4C4C4")))
    m["colour_ok"] = m.colour == m.colour_exp
    m["notscored_ok"] = m.not_scored_label == ~sc
    off = m[~(m.bar_ok & m.shade_ok & m.rho_ok & m.ci_ok & m.colour_ok & m.notscored_ok)]
    print(tag, name, "displayed expected", len(disp), "drawn", len(d), "not drawn but expected", len(missing_from_panel),
          [(r.cell_type, r.module, r["name"]) for _, r in missing_from_panel.iterrows()], "drawn but not expected", len(extra))
    print(tag, name, "rows with any mark off", len(off))
    if len(off):
        print(off[["group", "module", "label", "bar_total_d", "n_open", "bar_direct", "n_open_direct", "rho_d", "rho", "ci_low_d", "ci_low", "ci_high_d", "ci_high", "colour", "colour_exp"]].round(3).to_string())
    print(tag, name, "max |rho diff|", float((m.rho_d - m.rho).abs().max()), "max |ci diff|", float(np.nanmax([(m.ci_low_d - m.ci_low).abs().max(), (m.ci_high_d - m.ci_high).abs().max()])))
    return m
ms = compare(rs, "strict")
ma = compare(ra, "alt_inclusive")
# label text vs registry name (truncation check)
ms["label_name"] = ms.label.str.replace(r"^\d+\s+", "", regex=True)
ms["name_short"] = ms["name"].str.replace(r"\s*\(.*\)$", "", regex=True)
trunc = ms[ms.label_name != ms.name_short][["group", "module", "label_name", "name"]]
print(tag, "labels differing from registry name (minus parenthetical):", trunc.to_dict("records"))
ms.to_csv(f"{B}/geom_checks/5e_rows_{tag}.tsv", sep="\t", index=False)
