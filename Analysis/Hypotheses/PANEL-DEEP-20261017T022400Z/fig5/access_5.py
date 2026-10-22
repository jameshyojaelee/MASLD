"""Accessibility checks for Figure 5 panels.
1. Text: font sizes, bold fonts, text colours, WCAG 2.x contrast of each character against the filled rectangle under it
   (last drawn rect containing the glyph centre; page white otherwise).
2. Colour-vision simulation, Machado, Oliveira & Fernandes 2009 (severity 1.0) matrices applied in linear sRGB, for the
   mark-class colours of each panel; reports CIE76 delta E (CIELAB, D65) between each pair under normal, deuteranopia, protanopia.
Usage: python access_5.py <pdf_dir> <tag>
"""
import sys, itertools, collections
import numpy as np
import pdfplumber

d, tag = sys.argv[1], sys.argv[2]
B = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5"
PANELS = ["fig5a_input_firewall", "fig5b_protein_triage", "fig5c_mrna_protein_composite", "fig5d_snatac_accessibility",
          "fig5e_multimodal_program_summary", "fig5f_spatial_program_calibration"]
MARKS = {
    "fig5a_input_firewall": {"gene-context box": "#6A51A3", "program box": "#C9265E", "protein box": "#00838F"},
    "fig5b_protein_triage": {"concordant": "#3B6EA5", "discordant": "#B05A5A", "not significant": "#9E9E9E"},
    "fig5c_mrna_protein_composite": {"mRNA dot": "#8EC7E2", "protein dot": "#E8853A", "heatmap low": "#2166AC", "heatmap high": "#B2182B", "NA cell": "#BDBDBD"},
    "fig5d_snatac_accessibility": {"under 5%": "#D9D9D9", "5-<50%": "#8EC7E2", "one source": "#7B1FA2", "both sources": "#00695C"},
    "fig5e_multimodal_program_summary": {"robust up": "#B2182B", "robust down": "#2166AC", "not robust": "#C4C4C4", "enzyme": "#C79A3E", "direct": "#7A5A1E"},
    "fig5f_spatial_program_calibration": {"supported": "#007C91", "not supported": "#8C8C8C", "overlay low": "#1465C0", "overlay high": "#C9265E"},
}
PROT = np.array([[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]])
DEUT = np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]])

def h2rgb(h):
    return np.array([int(h[i:i + 2], 16) for i in (1, 3, 5)]) / 255.0
def lin(c):
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
def delin(c):
    c = np.clip(c, 0, 1); return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)
def lab(rgb):
    l = lin(rgb)
    M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
    xyz = M @ l / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > (6 / 29) ** 3, np.cbrt(xyz), xyz / (3 * (6 / 29) ** 2) + 4 / 29)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])
def sim(rgb, M):
    return delin(M @ lin(rgb))
def lum(rgb):
    l = lin(rgb); return 0.2126 * l[0] + 0.7152 * l[1] + 0.0722 * l[2]
def contrast(a, b):
    la, lb = sorted([lum(a), lum(b)], reverse=True); return (la + 0.05) / (lb + 0.05)
def tohex(c):
    if c is None or c == "": return "#000000"
    if isinstance(c, (int, float)): c = (c,)
    if len(c) == 1: c = (c[0],) * 3
    return "#%02X%02X%02X" % tuple(int(round(v * 255)) for v in c[:3])

out = []
for name in PANELS:
    p = pdfplumber.open(f"{d}/{name}.pdf").pages[0]
    fills = [r for r in p.rects if r.get("fill") and r["width"] < p.width * 0.95]
    rows = collections.defaultdict(lambda: [0, 99.0, ""])
    sizes, bold = set(), set()
    for ch in p.chars:
        if not ch["text"].strip(): continue
        cx, cy = (ch["x0"] + ch["x1"]) / 2, (ch["top"] + ch["bottom"]) / 2
        under = [r for r in fills if r["x0"] <= cx <= r["x1"] and r["top"] <= cy <= r["bottom"]]
        bg = tohex(under[-1]["non_stroking_color"]) if under else "#FFFFFF"
        fg = tohex(ch["non_stroking_color"])
        cr = contrast(h2rgb(fg), h2rgb(bg))
        k = (fg, bg); rows[k][0] += 1
        if cr < rows[k][1]: rows[k][1] = cr; rows[k][2] = ch["text"]
        sizes.add(round(ch["size"], 2) if ch["upright"] else round(ch["height"], 2))
        if "bold" in ch["fontname"].lower(): bold.add(ch["fontname"])
    print(f"== {tag} {name}: glyph sizes {sorted(sizes)} bold fonts {sorted(bold)}")
    for (fg, bg), (n, cr, ex) in sorted(rows.items(), key=lambda kv: kv[1][1]):
        flag = "LOW" if cr < 4.5 else ""
        print(f"   text {fg} on {bg}: {n} glyphs, contrast {cr:.2f} {flag}")
        out.append(("text", name, fg, bg, n, round(cr, 2), flag))
    cols = MARKS[name]
    for (a, ha), (b, hb) in itertools.combinations(cols.items(), 2):
        ra, rb = h2rgb(ha), h2rgb(hb)
        dn = np.linalg.norm(lab(ra) - lab(rb))
        dd = np.linalg.norm(lab(sim(ra, DEUT)) - lab(sim(rb, DEUT)))
        dp = np.linalg.norm(lab(sim(ra, PROT)) - lab(sim(rb, PROT)))
        flag = "LOW" if min(dd, dp) < 15 else ""
        print(f"   CVD {a} {ha} vs {b} {hb}: dE normal {dn:.1f} deut {dd:.1f} prot {dp:.1f} {flag}")
        out.append(("cvd", name, f"{a} {ha}", f"{b} {hb}", round(dn, 1), round(min(dd, dp), 1), flag))
import pandas as pd
pd.DataFrame(out, columns=["kind", "panel", "a", "b", "n_or_dE_normal", "contrast_or_min_dE_cvd", "flag"]).to_csv(f"{B}/geom_checks/access_{tag}.tsv", sep="\t", index=False)
