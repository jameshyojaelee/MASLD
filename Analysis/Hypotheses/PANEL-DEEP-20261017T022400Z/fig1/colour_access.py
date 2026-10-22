"""Colour-blind simulation (Machado, Oliveira & Fernandes 2009, severity 1.0, applied in linear sRGB)
and WCAG 2.x contrast ratios for the mark and text colours read from the Figure 1 PDFs.
Matrices are the published severity-1.0 values as commonly reproduced (recalled, not re-derived here).
"""
import itertools
import json
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent
M = {
    "protanopia": np.array([[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]]),
    "deuteranopia": np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]]),
    "tritanopia": np.array([[1.255528, -0.076749, -0.178779], [-0.078411, 0.930809, 0.147602], [0.004733, 0.691367, 0.303900]]),
}


def hex2rgb(h):
    h = h.lstrip("#"); return np.array([int(h[i:i + 2], 16) / 255 for i in (0, 2, 4)])


def lin(c):
    return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)


def delin(c):
    c = np.clip(c, 0, 1); return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)


def lab(rgb):
    l = lin(rgb)
    X, Y, Z = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]]) @ l
    xn, yn, zn = 0.95047, 1.0, 1.08883
    f = lambda t: np.where(t > (6 / 29) ** 3, np.cbrt(t), t / (3 * (6 / 29) ** 2) + 4 / 29)
    fx, fy, fz = f(X / xn), f(Y / yn), f(Z / zn)
    return np.array([116 * fy - 16, 500 * (fx - fy), 200 * (fy - fz)])


def lum(rgb):
    return float(np.dot([0.2126, 0.7152, 0.0722], lin(rgb)))


def contrast(a, b):
    la, lb = lum(hex2rgb(a)), lum(hex2rgb(b)); hi, lo = max(la, lb), min(la, lb)
    return (hi + 0.05) / (lo + 0.05)


def sim(h, kind):
    return delin(M[kind] @ lin(hex2rgb(h)))


panels = {
    "1B": {"Both (coloc+DEG)": "#00695C", "Colocalized only": "#1565C0", "Disease-state only": "#C9265E"},
    "1C": {"Liver enzyme only": "#F4A674", "Both": "#00695C", "Direct disease/PDFF only": "#1565C0"},
}
res = {}
for pn, cols in panels.items():
    res[pn] = {}
    for kind in ["normal"] + list(M):
        simc = {k: (hex2rgb(v) if kind == "normal" else sim(v, kind)) for k, v in cols.items()}
        pairs = {}
        for (a, ca), (b, cb) in itertools.combinations(simc.items(), 2):
            pairs[f"{a} vs {b}"] = round(float(np.linalg.norm(lab(ca) - lab(cb))), 1)
        res[pn][kind] = dict(sim_hex={k: "#%02X%02X%02X" % tuple(int(round(255 * x)) for x in v) for k, v in simc.items()},
                             deltaE76=pairs)
texts = {
    "1B black '26 observed' on #1565C0": ("#000000", "#1565C0"),
    "1B black '26 observed' on #C9265E": ("#000000", "#C9265E"),
    "1B white '400/93.9%' on #1565C0": ("#FFFFFF", "#1565C0"),
    "1B white '194/88.2%' on #C9265E": ("#FFFFFF", "#C9265E"),
    "1B black dashed line on #00695C": ("#000000", "#00695C"),
    "black text on white": ("#000000", "#FFFFFF"),
}
res["wcag_contrast"] = {k: round(contrast(*v), 2) for k, v in texts.items()}
json.dump(res, open(OUT / "colour_access.json", "w"), indent=1)
print(json.dumps(res, indent=1))
