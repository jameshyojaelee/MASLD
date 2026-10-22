"""Machado 2009 (severity 1.0) deuteranopia/protanopia simulation of Figure 6 mark colours,
CIE76 Delta E between classes before/after simulation, and WCAG 2.1 contrast of black text on fills."""
import numpy as np, itertools, json
MAT = {"protanopia": np.array([[0.152286, 1.052583, -0.204868], [0.114503, 0.786281, 0.099216], [-0.003882, -0.048116, 1.051998]]),
       "deuteranopia": np.array([[0.367322, 0.860646, -0.227968], [0.280085, 0.672501, 0.047413], [-0.011820, 0.042940, 0.968881]])}
def h2rgb(h): h = h.lstrip("#"); return np.array([int(h[i:i+2], 16) / 255 for i in (0, 2, 4)])
def lin(c): return np.where(c <= 0.04045, c / 12.92, ((c + 0.055) / 1.055) ** 2.4)
def delin(c): c = np.clip(c, 0, 1); return np.where(c <= 0.0031308, 12.92 * c, 1.055 * c ** (1 / 2.4) - 0.055)
def sim(h, kind): return delin(MAT[kind] @ lin(h2rgb(h)))
def lab(rgb):
    l = lin(rgb); M = np.array([[0.4124, 0.3576, 0.1805], [0.2126, 0.7152, 0.0722], [0.0193, 0.1192, 0.9505]])
    xyz = M @ l / np.array([0.95047, 1.0, 1.08883])
    f = np.where(xyz > (6/29) ** 3, np.cbrt(xyz), xyz / (3 * (6/29) ** 2) + 4/29)
    return np.array([116 * f[1] - 16, 500 * (f[0] - f[1]), 200 * (f[1] - f[2])])
def lum(rgb): l = lin(rgb); return 0.2126 * l[0] + 0.7152 * l[1] + 0.0722 * l[2]
def contrast(a, b): la, lb = sorted([lum(a), lum(b)], reverse=True); return (la + 0.05) / (lb + 0.05)
sets = {
 "6B cell fills (tested vs untested)": {"tested classes #E8EEF5": "#E8EEF5", "untested #EEEEEE": "#EEEEEE"},
 "6D bar fills (untested-no-RNA/other classes vs untested-no-genetic)": {"light #C9D4E0": "#C9D4E0", "dark #7F95AD": "#7F95AD"},
 "6A/6E box outlines": {"genetics #1565C0": "#1565C0", "RNA #C9265E": "#C9265E", "object #00695C": "#00695C", "neutral #9E9E9E": "#9E9E9E"},
}
res = {}
for name, cols in sets.items():
    print("==", name)
    for (n1, h1), (n2, h2) in itertools.combinations(cols.items(), 2):
        d0 = np.linalg.norm(lab(h2rgb(h1)) - lab(h2rgb(h2)))
        dd = np.linalg.norm(lab(sim(h1, "deuteranopia")) - lab(sim(h2, "deuteranopia")))
        dp = np.linalg.norm(lab(sim(h1, "protanopia")) - lab(sim(h2, "protanopia")))
        print(f"  {n1} vs {n2}: dE normal {d0:.1f}, deutan {dd:.1f}, protan {dp:.1f}")
        res[f"{name}: {n1} vs {n2}"] = dict(normal=round(d0, 1), deutan=round(dd, 1), protan=round(dp, 1))
print("== WCAG contrast, black text #000000 on fill")
for h in ["#E8EEF5", "#EEEEEE", "#C9D4E0", "#7F95AD", "#FFFFFF"]:
    cr = contrast(h2rgb("#000000"), h2rgb(h)); print(f"  {h}: {cr:.1f}:1"); res[f"black on {h}"] = round(cr, 1)
print("== outline/arrow contrast against white background (WCAG 1.4.11 non-text needs 3:1)")
for h in ["#9E9E9E", "#1565C0", "#C9265E", "#00695C"]:
    cr = contrast(h2rgb(h), h2rgb("#FFFFFF")); print(f"  {h} on white: {cr:.2f}:1"); res[f"{h} on white"] = round(cr, 2)
print("== fill vs white background (does the cell/bar boundary read without its outline?)")
for h in ["#E8EEF5", "#EEEEEE", "#C9D4E0"]:
    cr = contrast(h2rgb(h), h2rgb("#FFFFFF")); print(f"  {h} vs white: {cr:.2f}:1")
json.dump(res, open("cvd_contrast.json", "w"), indent=1)
