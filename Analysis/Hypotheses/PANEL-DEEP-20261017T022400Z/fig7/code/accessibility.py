"""Colour-vision and contrast checks for the Figure 7 panels.

1. Machado, Oliveira & Fernandes 2009 (IEEE TVCG 15:1291) severity-1.0 matrices for protanopia,
   deuteranopia and tritanopia, applied in linear sRGB. Pairwise CIEDE2000 distance between the colours
   each panel uses to separate categories, before and after simulation.
2. WCAG 2.x relative-luminance contrast: black text on every filled rectangle that has text drawn on it
   (7A bar segments), and every mark colour against the white background (WCAG 1.4.11 asks 3:1 for
   graphical objects).
Positive controls: three CIEDE2000 pairs from Sharma, Wu & Dalal 2005 (expected 2.0425, 2.3669, 27.1492);
each Machado matrix maps neutral gray to itself (rows sum to 1).

Usage: python accessibility.py <out_dir>
"""
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

OUT = Path(sys.argv[1])
GEOM = OUT

MACHADO = {
    "protanopia": np.array([[0.152286, 1.052583, -0.204868],
                            [0.114503, 0.786281, 0.099216],
                            [-0.003882, -0.048116, 1.051998]]),
    "deuteranopia": np.array([[0.367322, 0.860646, -0.227968],
                              [0.280085, 0.672501, 0.047413],
                              [-0.011820, 0.042940, 0.968881]]),
    "tritanopia": np.array([[1.255528, -0.076749, -0.178779],
                            [-0.078411, 0.930809, 0.147602],
                            [0.004733, 0.691367, 0.303900]]),
}


def hex2rgb(h):
    h = h.lstrip("#")
    return np.array([int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)])


def to_lin(c):
    c = np.asarray(c, float)
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def to_srgb(c):
    c = np.clip(c, 0, 1)
    return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)


def simulate(h, kind):
    return to_srgb(MACHADO[kind] @ to_lin(hex2rgb(h)))


def lab(rgb):
    lin = to_lin(rgb)
    M = np.array([[0.4124564, 0.3575761, 0.1804375],
                  [0.2126729, 0.7151522, 0.0721750],
                  [0.0193339, 0.1191920, 0.9503041]])
    xyz = M @ lin / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > (6 / 29) ** 3, np.cbrt(xyz), xyz / (3 * (6 / 29) ** 2) + 4 / 29)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])


def de2000(l1, l2):
    L1, a1, b1 = l1
    L2, a2, b2 = l2
    C1, C2 = np.hypot(a1, b1), np.hypot(a2, b2)
    Cb = (C1 + C2) / 2
    G = 0.5 * (1 - np.sqrt(Cb ** 7 / (Cb ** 7 + 25 ** 7)))
    a1p, a2p = (1 + G) * a1, (1 + G) * a2
    C1p, C2p = np.hypot(a1p, b1), np.hypot(a2p, b2)
    h1p = np.degrees(np.arctan2(b1, a1p)) % 360 if C1p else 0.0
    h2p = np.degrees(np.arctan2(b2, a2p)) % 360 if C2p else 0.0
    dLp, dCp = L2 - L1, C2p - C1p
    if C1p * C2p == 0:
        dhp = 0.0
    else:
        dhp = h2p - h1p
        if dhp > 180:
            dhp -= 360
        elif dhp < -180:
            dhp += 360
    dHp = 2 * np.sqrt(C1p * C2p) * np.sin(np.radians(dhp / 2))
    Lbp, Cbp = (L1 + L2) / 2, (C1p + C2p) / 2
    if C1p * C2p == 0:
        hbp = h1p + h2p
    elif abs(h1p - h2p) <= 180:
        hbp = (h1p + h2p) / 2
    elif h1p + h2p < 360:
        hbp = (h1p + h2p + 360) / 2
    else:
        hbp = (h1p + h2p - 360) / 2
    T = (1 - 0.17 * np.cos(np.radians(hbp - 30)) + 0.24 * np.cos(np.radians(2 * hbp))
         + 0.32 * np.cos(np.radians(3 * hbp + 6)) - 0.20 * np.cos(np.radians(4 * hbp - 63)))
    dtheta = 30 * np.exp(-(((hbp - 275) / 25) ** 2))
    Rc = 2 * np.sqrt(Cbp ** 7 / (Cbp ** 7 + 25 ** 7))
    Sl = 1 + 0.015 * (Lbp - 50) ** 2 / np.sqrt(20 + (Lbp - 50) ** 2)
    Sc = 1 + 0.045 * Cbp
    Sh = 1 + 0.015 * Cbp * T
    Rt = -np.sin(np.radians(2 * dtheta)) * Rc
    return float(np.sqrt((dLp / Sl) ** 2 + (dCp / Sc) ** 2 + (dHp / Sh) ** 2 + Rt * (dCp / Sc) * (dHp / Sh)))


def rel_lum(h):
    lin = to_lin(hex2rgb(h))
    return float(0.2126 * lin[0] + 0.7152 * lin[1] + 0.0722 * lin[2])


def contrast(h1, h2):
    a, b = sorted([rel_lum(h1), rel_lum(h2)], reverse=True)
    return (a + 0.05) / (b + 0.05)


# positive controls
controls = []
for l1, l2, exp in [((50, 2.6772, -79.7751), (50, 0, -82.7485), 2.0425),
                    ((50, 0, 0), (50, -1, 2), 2.3669),
                    ((50, 2.5, 0), (73, 25, -18), 27.1492)]:
    got = de2000(np.array(l1, float), np.array(l2, float))
    controls.append({"check": f"CIEDE2000 {l1} vs {l2}", "expected": exp, "got": round(got, 4), "pass": abs(got - exp) < 1e-3})
for k, M in MACHADO.items():
    g = simulate("#808080", k)
    controls.append({"check": f"{k} keeps #808080 gray", "expected": 0.50196, "got": [round(float(v), 5) for v in g],
                     "pass": bool(np.allclose(g, hex2rgb("#808080"), atol=2e-5))})
controls.append({"check": "WCAG black on white", "expected": 21.0, "got": round(contrast("#000000", "#FFFFFF"), 3),
                 "pass": abs(contrast("#000000", "#FFFFFF") - 21) < 1e-9})

# category colour sets per panel, from the colour census of the rendered PDFs
SETS = {
    "7B point/interval: above / below / crosses zero": ["#C9265E", "#1565C0", "#9E9E9E"],
    "7C fill: exposure resolved / includes training-seen": ["#C9265E", "#F4A674"],
    "7D: achieved (magenta) / bound or chance (gray)": ["#C9265E", "#9E9E9E"],
    "7E: supported / untestable / tested negative": ["#C9265E", "#9E9E9E", "#1565C0"],
    "7F heatmap: confounded / clean / unresolved; forest: above / crosses zero": ["#C9265E", "#1565C0", "#9E9E9E"],
    "7A bar segments": ["#1565C0", "#C9265E", "#BDBDBD", "#9E9E9E", "#6E6E6E"],
}
pair_rows = []
for name, cols in SETS.items():
    for c1, c2 in combinations(cols, 2):
        row = {"set": name, "c1": c1, "c2": c2, "normal": de2000(lab(hex2rgb(c1)), lab(hex2rgb(c2)))}
        for k in MACHADO:
            row[k] = de2000(lab(simulate(c1, k)), lab(simulate(c2, k)))
        row["min_cvd"] = min(row[k] for k in MACHADO)
        pair_rows.append(row)

sim_rows = []
for h in sorted({c for v in SETS.values() for c in v}):
    sim_rows.append({"colour": h, **{k: "#%02X%02X%02X" % tuple(int(round(x * 255)) for x in simulate(h, k)) for k in MACHADO},
                     "contrast_vs_white": round(contrast(h, "#FFFFFF"), 3), "black_text_contrast": round(contrast(h, "#000000"), 3)})

# text-on-fill: which characters sit inside a filled rectangle, per panel (main and calib)
text_rows = []
for f in sorted(GEOM.glob("geom_*_fig7*.json")):
    g = json.load(open(f))
    rects = [r for r in g["rects"] if r.get("filled") and r["fill"] not in ("#FFFFFF", None)
             and (r["x1"] - r["x0"]) > 8 and (r["bottom"] - r["top"]) > 4]
    for r in rects:
        inside = [c for c in g["chars"] if c["x0"] >= r["x0"] - 0.5 and c["x1"] <= r["x1"] + 0.5
                  and c["top"] >= r["top"] - 0.5 and c["bottom"] <= r["bottom"] + 0.5 and c["text"].strip()]
        if inside:
            txt = "".join(c["text"] for c in sorted(inside, key=lambda c: (round(c["top"]), c["x0"])))
            tf = inside[0]["fill"]
            text_rows.append({"file": f.name, "rect_fill": r["fill"], "text_fill": tf, "text": txt[:60],
                              "wcag_contrast": round(contrast(r["fill"], tf), 3),
                              "passes_4.5_normal_text": contrast(r["fill"], tf) >= 4.5})

import pandas as pd  # noqa: E402

pd.DataFrame(pair_rows).to_csv(OUT / "accessibility_cvd_pairs.tsv", sep="\t", index=False, float_format="%.2f")
pd.DataFrame(sim_rows).to_csv(OUT / "accessibility_colours.tsv", sep="\t", index=False)
pd.DataFrame(text_rows).to_csv(OUT / "accessibility_text_on_fill.tsv", sep="\t", index=False)
json.dump(controls, open(OUT / "accessibility_controls.json", "w"), indent=1, default=str)
print(pd.DataFrame(controls).to_string())
print(pd.DataFrame(pair_rows).to_string(float_format=lambda v: "%.1f" % v))
print(pd.DataFrame(sim_rows).to_string())
print(pd.DataFrame(text_rows).to_string())
