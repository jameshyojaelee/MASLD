"""Colour-vision simulation (Machado, Oliveira & Fernandes 2009, severity 1.0, applied in linear sRGB) for the mark colours
used in Figure 4, with CIEDE2000 between pairs that must be told apart; WCAG 2.1 contrast of text on its background."""
import itertools
import numpy as np
MACHADO = {"deuteranopia": np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]]),
           "protanopia": np.array([[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]])}
def h2r(h): h = h.lstrip("#"); return np.array([int(h[i:i + 2], 16) for i in (0, 2, 4)]) / 255
def lin(c): return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
def delin(c): c = np.clip(c, 0, 1); return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)
def sim(h, kind): return delin(MACHADO[kind] @ lin(h2r(h))) if kind != "normal" else h2r(h)
def lab(c):
    M = np.array([[0.4124564, 0.3575761, 0.1804375], [0.2126729, 0.7151522, 0.0721750], [0.0193339, 0.1191920, 0.9503041]])
    xyz = M @ lin(np.asarray(c, float)) / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > (6 / 29) ** 3, np.cbrt(xyz), xyz / (3 * (6 / 29) ** 2) + 4 / 29)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])
def dE(a, b):
    """CIEDE2000 (Sharma, Wu & Dalal 2005)."""
    L1, a1, b1 = lab(a); L2, a2, b2 = lab(b)
    C1, C2 = np.hypot(a1, b1), np.hypot(a2, b2); Cb = (C1 + C2) / 2
    G = 0.5 * (1 - np.sqrt(Cb ** 7 / (Cb ** 7 + 25 ** 7)))
    a1p, a2p = (1 + G) * a1, (1 + G) * a2
    C1p, C2p = np.hypot(a1p, b1), np.hypot(a2p, b2)
    h1p, h2p = np.degrees(np.arctan2(b1, a1p)) % 360, np.degrees(np.arctan2(b2, a2p)) % 360
    dLp, dCp = L2 - L1, C2p - C1p
    dhp = 0 if C1p * C2p == 0 else (h2p - h1p if abs(h2p - h1p) <= 180 else (h2p - h1p - 360 if h2p - h1p > 180 else h2p - h1p + 360))
    dHp = 2 * np.sqrt(C1p * C2p) * np.sin(np.radians(dhp / 2))
    Lbp, Cbp = (L1 + L2) / 2, (C1p + C2p) / 2
    if C1p * C2p == 0: hbp = h1p + h2p
    elif abs(h1p - h2p) <= 180: hbp = (h1p + h2p) / 2
    else: hbp = (h1p + h2p + 360) / 2 if h1p + h2p < 360 else (h1p + h2p - 360) / 2
    T = 1 - 0.17 * np.cos(np.radians(hbp - 30)) + 0.24 * np.cos(np.radians(2 * hbp)) + 0.32 * np.cos(np.radians(3 * hbp + 6)) - 0.20 * np.cos(np.radians(4 * hbp - 63))
    dth = 30 * np.exp(-((hbp - 275) / 25) ** 2); Rc = 2 * np.sqrt(Cbp ** 7 / (Cbp ** 7 + 25 ** 7))
    SL = 1 + 0.015 * (Lbp - 50) ** 2 / np.sqrt(20 + (Lbp - 50) ** 2); SC = 1 + 0.045 * Cbp; SH = 1 + 0.015 * Cbp * T
    RT = -np.sin(np.radians(2 * dth)) * Rc
    return float(np.sqrt((dLp / SL) ** 2 + (dCp / SC) ** 2 + (dHp / SH) ** 2 + RT * (dCp / SC) * (dHp / SH)))
def lum(h): c = lin(h2r(h)); return 0.2126 * c[0] + 0.7152 * c[1] + 0.0722 * c[2]
def contrast(a, b): la, lb = sorted([lum(a), lum(b)], reverse=True); return (la + 0.05) / (lb + 0.05)
sets = {
    "4C/4H direction (pos, neg, opposite)": {"pos": "#C9265E", "neg4C": "#1565C0", "neg4H": "#2166AC", "opposite": "#9E9E9E"},
    "4C lineage bars": {"Hep": "#0D47A1", "Fib": "#F57F17", "Mac": "#C2185B", "Chol": "#1565C0"},
    "4D cell fills": {"diag": "#E6E6E6", "offdiag": "#F2B6CD", "zero": "#FFFFFF"},
    "4F cohorts": {"GSE126848": "#BDBDBD", "GSE130970": "#7A7A7A", "GSE135251": "#4D4D4D", "GSE162694": "#C9265E", "GSE213621": "#1565C0"},
    "4B cell groups": {"Hepatocytes": "#1B5E20", "Cholangiocytes": "#00BCD4", "Endothelial": "#7CB342", "Plasma": "#8D6E63", "NK": "#00897B",
                        "Macrophages": "#212121", "T": "#1A237E", "Monocytes": "#EC407A", "Fibroblasts": "#FF6F00", "Granulocytes": "#F9A825",
                        "B": "#5D4037", "Dendritic": "#6A1B9A"},
    "4G/4B diverging ends": {"low": "#1565C0", "high": "#C9265E", "mid": "#F7F7F7"},
}
for name, d in sets.items():
    print(f"== {name}")
    for kind in ["normal", "deuteranopia", "protanopia"]:
        pairs = [(a, b, dE(sim(d[a], kind), sim(d[b], kind))) for a, b in itertools.combinations(d, 2)]
        low = sorted(pairs, key=lambda t: t[2])[:4]
        print(f"  {kind:13s} min dE2000 pairs:", ", ".join(f"{a}/{b} {v:.1f}" for a, b, v in low))
print("== WCAG contrast (text on background; AA needs 4.5 for normal text)")
for t, bg, what in [("#000000", "#E6E6E6", "4D black on diagonal gray"), ("#000000", "#F2B6CD", "4D black on off-diagonal pink"),
                    ("#9E9E9E", "#FFFFFF", "4D gray zero-count text on white"), ("#4D4D4D", "#FFFFFF", "4A dark-gray text on white"),
                    ("#C9265E", "#FFFFFF", "4A magenta text on white"), ("#000000", "#FFFFFF", "4G labels on white boxes")]:
    print(f"  {what}: {contrast(t, bg):.2f}")
