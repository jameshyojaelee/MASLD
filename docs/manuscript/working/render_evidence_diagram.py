#!/usr/bin/env python3
"""atlas_evidence_diagram.pdf — compact, larger font."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

OUT = "docs/manuscript/working/atlas_evidence_diagram.pdf"

# ── compact node layout ────────────────────────────────────────────────────
#   3 rows:   RNA y=7.5  |  gene-classes y=5.2  |  modalities y=2.5
#   5 x-cols: 1.5 · 4.0 · 6.8 · 9.5 · 12.3

NODES = {
    "RNA":  (6.8,  7.5, 2.8, 1.1,
             "Bulk RNA-seq DEGs\nmixed signal",         "#EBF5FB"),
    "CAU":  (1.5,  5.2, 2.4, 1.0,  "Causal genes",    "#FDEBD0"),
    "DRUG": (4.0,  5.2, 2.2, 1.0,  "Drug targets",    "#E8F8F5"),
    "REA":  (6.8,  5.2, 2.4, 1.0,  "Reactive genes",  "#EBF5FB"),
    "BIO":  (9.5,  5.2, 2.2, 1.0,  "Biomarkers",      "#FEF9E7"),
    "COM":  (12.3, 5.2, 2.4, 1.1,  "Composition-\ndriven", "#F2F3F4"),
    "GEN":  (1.5,  2.5, 3.4, 1.65,
             "Genetic causal layer\nGWAS · eQTL COLOC\nFinemapping · scATAC-seq",
             "#EBF5EB"),
    "VAL":  (6.8,  2.5, 3.4, 1.1,
             "Cross-modal validation\nProteomics · Spatial",
             "#F5EEF8"),
    "CELL": (12.3, 2.5, 3.4, 1.1,
             "Cell-resolution layer\nscRNA-seq · Deconvolution",
             "#FEF9E7"),
}

BORDER  = "#37474F"
ARROW_C = "#37474F"
FS_MOD  = 13.0
FS_GENE = 16.0
FS_LBL  = 12.0

fig, ax = plt.subplots(figsize=(14, 8))
ax.set_xlim(-0.3, 14.5)
ax.set_ylim(1.2, 8.8)
ax.axis("off")
fig.patch.set_facecolor("white")

def draw_node(name):
    cx, cy, w, h, text, color = NODES[name]
    ax.add_patch(FancyBboxPatch(
        (cx - w/2, cy - h/2), w, h,
        boxstyle="round,pad=0.09",
        facecolor=color, edgecolor=BORDER, linewidth=1.6, zorder=3))
    fs = FS_MOD if name in ("GEN", "VAL", "CELL") else FS_GENE
    ax.text(cx, cy, text, ha="center", va="center",
            fontsize=fs, multialignment="center", zorder=4)

for n in NODES:
    draw_node(n)

def pt(name, side, dx=0.0, dy=0.0):
    cx, cy, w, h = NODES[name][:4]
    if side == "t": return (cx + dx, cy + h/2 + dy)
    if side == "b": return (cx + dx, cy - h/2 + dy)
    if side == "l": return (cx - w/2 + dx, cy + dy)
    if side == "r": return (cx + w/2 + dx, cy + dy)

def arr(p1, p2, rad=0.0, dashed=False, lw=1.5):
    ax.annotate("", xy=p2, xytext=p1,
                arrowprops=dict(
                    arrowstyle="->", color=ARROW_C, lw=lw,
                    linestyle="--" if dashed else "-",
                    connectionstyle=f"arc3,rad={rad}",
                    shrinkA=0, shrinkB=3),
                zorder=2)

def lbl(text, x, y):
    ax.text(x, y, text, fontsize=FS_LBL, ha="center", va="center",
            color=BORDER, style="italic",
            bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                      edgecolor="none", alpha=0.93),
            zorder=5)

# RNA fans to the three gene classes
arr(pt("RNA","b"), pt("CAU","t"))
arr(pt("RNA","b"), pt("REA","t"))
arr(pt("RNA","b"), pt("COM","t"))

# Gene-class → output (short horizontal)
arr(pt("CAU","r"), pt("DRUG","l"))
arr(pt("REA","r"), pt("BIO","l"))

# GEN → CAU  (straight up, left column)
arr(pt("GEN","t"), pt("CAU","b"))
lbl("identifies + mechanism", 0.35, 3.9)

# CELL → COM  (straight up, right column)
arr(pt("CELL","t"), pt("COM","b"))
lbl("removes artifact", 13.8, 3.9)

# CELL → REA  (diagonal, x ∈ [6.8, 12.3])
arr(pt("CELL","t", dx=-0.35), pt("REA","b", dx=+0.45), rad=0.10)
lbl("cell-intrinsic refinement", 10.4, 4.25)

# VAL → REA  (straight up, centre column, dashed)
arr(pt("VAL","t"), pt("REA","b", dx=-0.25), dashed=True)

# VAL → CAU  (diagonal, x ∈ [1.5, 6.8], bows down to clear Drug targets, dashed)
arr(pt("VAL","t", dx=-0.45), pt("CAU","b", dx=+0.25), rad=+0.18, dashed=True)

# Single label for both VAL validation arrows
lbl("protein + tissue context", 4.5, 3.25)

plt.tight_layout(pad=0.2)
fig.savefig(OUT, format="pdf", bbox_inches="tight", dpi=150)
print(f"Saved → {OUT}")
