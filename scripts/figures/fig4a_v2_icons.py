#!/usr/bin/env python3
"""
Fig 4A v2 -- CANDIDATE 6: ICON + ASSAYED-vs-VALIDATED BARS (PI's suggestion).

PI's brief (paraphrased): "icons for RNA-seq and genetics; bar graphs; assayed vs
validated." Read as: give each data modality a recognizable ICON (not a text label),
and per validation assay show how many prioritized targets it ASSAYED (measured) vs
how many VALIDATED (were significant). RNA-seq + genetics are the two prioritization
INPUTS (icons + nominated counts) that define the 9,882-gene target list; the three
orthogonal assays (proteomics / spatial / snATAC) carry the assayed-vs-validated bars.

STRICT definition (default in fig4a_v2_common): validated = per-assay disease
SIGNIFICANCE TEST only, so the assayed-vs-validated yield is genuine for all three:
proteomics 3,144->714 (~23%), spatial (Visium-hep DE) 164->112 (~68%), snATAC
(hepatocyte DA) 501->233 (~46%). The curated GeoMx zonation signature and static
promoter/TF-activity flags are excluded (under the inclusive definition they inflate
spatial to 447 and snATAC to 560, saturating the yield). All numbers read live from d.

House style (fig4a_v2_common sets rcParams): PDF, Helvetica 6, all text INK, no bold,
NO gene names, colour on marks only from the modality palette, control/absent GRAY,
ASCII-only on-panel. Icons are hand-drawn matplotlib patches (no external art),
ASPECT-CORRECTED: each icon takes separate x/y half-sizes (sx in gene-count units,
sy in row units) derived from the axes aspect, so lines/dots are not stretched.
Every bar/label encodes a live count. Full claim in the printed caption.
"""
from fig4a_v2_common import (get_data, save_pdf, INK, FS, GRAY, C_NS, C_CONV,
    C_BULK, C_GWAS, C_PROTEO, C_SPATIAL, C_ATAC, LENS_COL, LENS_NAME,
    tint, fold_ci, _ribbon, plt, np, Rectangle, PathPatch, FancyBboxPatch,
    Circle, Wedge, Ellipse, FancyArrowPatch, Polygon, Path)


# ---- hand-drawn, aspect-corrected modality icons ---------------------------
# Each takes (ax, cx, cy, sx, sy, c): sx = x half-size (gene-count units),
# sy = y half-size (row units). Built from lines + aspect-correct ellipses only
# (no RegularPolygon, which cannot be aspect-corrected in unequal data scales).
def _dot(ax, cx, cy, sx, sy, c, fill=True):
    ax.add_patch(Ellipse((cx, cy), 2 * sx, 2 * sy,
                         facecolor=(c if fill else "none"), edgecolor=c,
                         lw=0.8, zorder=7))


def ic_rnaseq(ax, cx, cy, sx, sy, c):
    """mRNA transcript: a single wavy strand."""
    xs = np.linspace(cx - sx, cx + sx, 40)
    ys = cy + 0.72 * sy * np.sin((xs - cx) / sx * 2.7)
    ax.plot(xs, ys, color=c, lw=1.1, solid_capstyle="round", zorder=7)


def ic_genetics(ax, cx, cy, sx, sy, c):
    """DNA double helix: two offset sine strands + rungs."""
    xs = np.linspace(cx - sx, cx + sx, 40)
    y1 = cy + 0.85 * sy * np.sin((xs - cx) / sx * 3.0)
    y2 = cy + 0.85 * sy * np.sin((xs - cx) / sx * 3.0 + np.pi)
    ax.plot(xs, y1, color=c, lw=1.0, zorder=7)
    ax.plot(xs, y2, color=c, lw=1.0, zorder=7)
    for xr in np.linspace(cx - 0.7 * sx, cx + 0.7 * sx, 4):
        yr1 = cy + 0.85 * sy * np.sin((xr - cx) / sx * 3.0)
        yr2 = cy + 0.85 * sy * np.sin((xr - cx) / sx * 3.0 + np.pi)
        ax.plot([xr, xr], [yr1, yr2], color=c, lw=0.5, alpha=0.8, zorder=6)


def ic_proteomics(ax, cx, cy, sx, sy, c):
    """DIA-MS spectrum: baseline + peaks of varying height."""
    base = cy - 0.85 * sy
    ax.plot([cx - sx, cx + sx], [base, base], color=c, lw=0.8, zorder=6)
    hs = [0.7, 1.5, 0.5, 1.2, 0.9]
    xs = np.linspace(cx - 0.8 * sx, cx + 0.8 * sx, len(hs))
    for x, h in zip(xs, hs):
        ax.plot([x, x], [base, base + h * sy], color=c, lw=1.1, zorder=7)


def ic_spatial(ax, cx, cy, sx, sy, c):
    """Tissue capture: a small grid of Visium-like spots (aspect-correct dots)."""
    for gx in (-0.6, 0.0, 0.6):
        for gy in (-0.45, 0.45):
            _dot(ax, cx + gx * sx, cy + gy * sy, 0.16 * sx, 0.16 * sy, c, fill=False)


def ic_atac(ax, cx, cy, sx, sy, c):
    """Open chromatin: DNA line threaded through 3 nucleosome beads."""
    ax.plot([cx - sx, cx + sx], [cy, cy], color=c, lw=0.7, zorder=6)
    for bx in (-0.62, 0.0, 0.62):
        _dot(ax, cx + bx * sx, cy, 0.26 * sx, 0.26 * sy, c, fill=False)


ICON = {"rnaseq": ic_rnaseq, "genetics": ic_genetics, "proteomics": ic_proteomics,
        "spatial": ic_spatial, "snATAC": ic_atac}


def build(d, stats):
    uni = int(d["universe"])
    sub = int(d["substrate"])          # RNA-seq nominated
    col = int(d["coloc"])              # genetics nominated

    # validation assays: (name, colour, icon key, assayed=measured, validated)
    ROWS = [
        ("proteomics", C_PROTEO, "proteomics", int(d["P_meas"]), int(d["P_val"])),
        ("spatial",    C_SPATIAL, "spatial",   int(d["S_meas"]), int(d["S_val"])),
        ("snATAC",     C_ATAC,    "snATAC",    int(d["A_meas"]), int(d["A_val"])),
    ]
    xmax = max(m for _, _, _, m, _ in ROWS)             # shared count axis (=3,144)

    # ---- geometry / canvas ----
    ROWDY = 1.10
    ys = [0.15 - i * ROWDY for i in range(len(ROWS))]   # proteomics/spatial/snATAC
    TH = 0.52
    GUT = xmax * 0.48                                   # left gutter (icons + names)
    XLO, XHI = -GUT - xmax * 0.03, xmax * 1.50
    YLO, YHI = ys[-1] - 0.95, 2.35

    figw, figh = 4.3, 3.0
    fig, ax = plt.subplots(figsize=(figw, figh))
    ax.set_xlim(XLO, XHI)
    ax.set_ylim(YLO, YHI)
    ax.axis("off")

    # aspect-correct icon half-sizes: choose sx in x-units, derive sy so the icon is
    # visually square (sx/x_per_in == sy/y_per_in).
    sx = xmax * 0.028
    sy = sx * ((YHI - YLO) * figw) / ((XHI - XLO) * figh)
    icx = -GUT + xmax * 0.055          # icon centre x (in the gutter)
    namex = -GUT + xmax * 0.135        # name/label left edge (clear of icon)

    # ---- header: RNA-seq + genetics inputs (icons + nominated counts) ----
    ic_rnaseq(ax, icx, 2.02, sx, sy, C_BULK)
    ax.text(namex, 2.02, f"RNA-seq {sub:,}", ha="left", va="center", color=INK, zorder=8)
    ic_genetics(ax, icx, 1.50, sx, sy, C_GWAS)
    ax.text(namex, 1.50, f"genetics {col:,}", ha="left", va="center", color=INK, zorder=8)

    # ---- section header (folds in the 9,882 target list; no separate arrow) ----
    ax.text(-GUT, 0.92,
            f"orthogonal validation of {uni:,} prioritized targets: assayed vs validated",
            ha="left", va="center", color=INK, zorder=8)

    # ---- assayed-vs-validated bars, one per assay, with icons ----
    for (name, c, ic, meas, val), yr in zip(ROWS, ys):
        yld = (100.0 * val / meas) if meas else 0.0
        # assayed (measured) bar = light body + thin outline
        ax.add_patch(Rectangle((0, yr - TH / 2), meas, TH, facecolor=C_NS,
                               edgecolor=GRAY, lw=0.5, zorder=3))
        # validated fill from the left = modality colour
        ax.add_patch(Rectangle((0, yr - TH / 2), val, TH, facecolor=c,
                               edgecolor="none", zorder=4))
        ICON[ic](ax, icx, yr, sx, sy, c)
        ax.text(namex, yr, name, ha="left", va="center", color=INK, zorder=8)
        ax.text(meas + xmax * 0.02, yr, f"{val:,} / {meas:,}  ({yld:.0f}%)",
                ha="left", va="center", color=INK, zorder=8)

    # ---- x baseline + count ticks (bar length reads as a real count) ----
    yb = ys[-1] - 0.52
    ax.plot([0, xmax], [yb, yb], color=INK, lw=0.5, zorder=3)
    for t in (0, 1000, 2000, 3000):
        if t <= xmax * 1.02:
            ax.plot([t, t], [yb, yb - 0.08], color=INK, lw=0.5, zorder=3)
            ax.text(t, yb - 0.16, f"{t:,}", ha="center", va="top", color=INK, zorder=8)
    ax.text(xmax / 2.0, yb - 0.40, "genes (of the prioritized universe)",
            ha="center", va="top", color=INK, zorder=8)

    # ---- legend (top-right): assayed vs validated ----
    lx = xmax * 1.06
    for ly, fc, ec, txt in ((2.02, C_NS, GRAY, "assayed (measured)"),
                            (1.50, "#6E6E6E", "none", "validated (significant)")):
        ax.add_patch(Rectangle((lx, ly - 0.13), xmax * 0.085, 0.26, facecolor=fc,
                               edgecolor=ec, lw=0.5, zorder=4))
        ax.text(lx + xmax * 0.11, ly, txt, ha="left", va="center", color=INK, zorder=8)

    save_pdf(fig, "fig4a_icons.pdf")

    print("[caption:icons] Assay throughput per modality (STRICT definition: validated = "
          "disease SIGNIFICANCE TEST only). RNA-seq ({sub:,}) and genetics ({col:,}) "
          "nominate the {uni:,}-gene prioritized target list (icons, top). Each orthogonal "
          "assay reports ASSAYED (prioritized targets it measured/tested) vs VALIDATED "
          "(significant): proteomics {pm:,} -> {pv:,} ({py:.0f}%; DIA-MS DE -- liver "
          "{plm:,}->{plv:,}, plasma {ppm:,}->{ppv:,}), spatial {sm:,} -> {sv:,} ({sy:.0f}%; "
          "Visium-hep Wilcoxon DE), snATAC {am:,} -> {av:,} ({ay:.0f}%; hepatocyte "
          "differential accessibility). All three are genuine measured->significant yields; "
          "the curated GeoMx MASH zonation signature and static promoter/TF-activity flags "
          "are EXCLUDED (they would inflate spatial to 447 and snATAC to 560, saturating the "
          "yield). Every bar/label is a live count; no gene names.".format(
              sub=sub, col=col, uni=uni,
              pm=int(d["P_meas"]), pv=int(d["P_val"]),
              py=100.0 * d["P_val"] / max(d["P_meas"], 1),
              sm=int(d["S_meas"]), sv=int(d["S_val"]),
              sy=100.0 * d["S_val"] / max(d["S_meas"], 1),
              am=int(d["A_meas"]), av=int(d["A_val"]),
              ay=100.0 * d["A_val"] / max(d["A_meas"], 1),
              plm=int(d["P_liver_meas"]), plv=int(d["P_liver_val"]),
              ppm=int(d["P_plasma_meas"]), ppv=int(d["P_plasma_val"])))


if __name__ == "__main__":
    d, stats = get_data()
    build(d, stats)
