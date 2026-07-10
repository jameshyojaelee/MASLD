#!/usr/bin/env python3
"""
Fig 4A v2 -- CANDIDATE 4: PHYSICAL-LAYER SCHEMATIC (the therapeutic-engagement map).

The claim this geometry welds on (vs the retired display-only Sankey):
  a MASLD target prioritized from DNA and mRNA is CORROBORATED at the protein,
  chromatin and tissue layers a therapeutic actually engages -- and of those, the
  protein-and-tissue convergence is the ONE that exceeds chance.

The device: a rigorous DATA-schematic (NOT a cartoon) of one stylized hepatocyte
cross-section in LIGHT-GRAY anatomy, with each evidence modality pinned to the
PHYSICAL LAYER it measures, read top-left -> bottom-right along biology:

  DNA (genetics, C_GWAS)           -- a short double-helix inside the nucleus;
                                      badge d["coloc"] prioritized causal genes.
  OPEN CHROMATIN gate (snATAC)     -- nucleosome beads unspooling at the nuclear
                                      edge; badge d["A_val"] ("gate" = access, not a
                                      central-dogma rung).
  mRNA (transcriptomics, C_BULK)   -- a wavy transcript through a nuclear pore into
                                      cytoplasm; badge d["substrate"]. A gray bracket
                                      over DNA+mRNA prints the prioritized universe
                                      d["universe"] (the two origins that define it).
  PROTEIN (proteomics, C_PROTEO)   -- a folded ribbon at a ribosome (badge liver
                                      d["P_liver_val"]) with one ribbon crossing the
                                      membrane into the vessel lumen = secreted
                                      (badge plasma d["P_plasma_val"]).
  TISSUE frame (spatial, C_SPATIAL)-- a hex-tiled tissue strip under the cell; badge
                                      d["S_val"] ("frame" = spatial context, not a rung).

  CONVERGENT CORE -- a filled STAR where the secreted-protein ribbon meets the tissue
  strip (protein and tissue physically coincide): d["_part"]["PS"] genes, tagged live
  with fold + p from stats["conv"]["PS"] -- the one convergence beating chance. A
  compact tally box carries the honest corroboration ladder (>=1 / >=2 / all-3), and a
  faint GRAY honesty badge by the DNA glyph records that genetics-alone (genetic-only
  d["_part"]["genetic_only"]) enriches in NOTHING (n.s.).

House style (fig4a_v2_common sets rcParams -- do NOT re-set): PDF, Helvetica 6, all
text INK, no bold, NO gene names, saturated colour ONLY on the count badges + the one
star (anatomy is light-gray line/faint-fill so the data ink dominates), control/absent
GRAY, CB-safe, ASCII-only on-panel ("open chromatin (gate)", "tissue (frame)", ">=2",
"all 3", "1.36x", "p=2e-3", "n.s."). Every badge prints a LIVE count; fold/p guarded so
the panel always renders. Full claim in the printed caption, not on-panel.
"""
from fig4a_v2_common import (get_data, save_pdf, INK, FS, GRAY, C_NS, C_CONV,
    C_BULK, C_GWAS, C_PROTEO, C_SPATIAL, C_ATAC, LENS_COL, LENS_NAME,
    tint, fold_ci, _ribbon, plt, np, Rectangle, PathPatch, FancyBboxPatch,
    Circle, Wedge, Ellipse, FancyArrowPatch, Polygon, Path)

# Light-gray anatomy palette (line + faint fill) so schematic outlines never compete
# with the coloured data badges -- house rule "anatomy = light gray line/shade only".
ANAT_LINE = "#BDBDBD"
ANAT_FILL = "#F6F6F6"
# white-backed label box so on-cell text lifts off the faint anatomy fill.
_WBOX = dict(boxstyle="round,pad=0.16", fc="white", ec="#E2E2E2", lw=0.3, alpha=0.94)


# -- small helpers ------------------------------------------------------------
def _pfmt(p):
    """Compact ASCII p-string (Helvetica lacks fancy glyphs): '2e-3', '0.30', 'n/a'."""
    if p is None or not np.isfinite(p):
        return "n/a"
    if p >= 0.01:
        return f"{p:.2f}"
    m, e = f"{p:.0e}".split("e")
    return f"{m}e{int(e)}"


def _iv(d, key, default=0):
    """Live int read with a guard so a missing/NaN source never breaks the panel."""
    try:
        v = d.get(key, default)
        return int(v) if v is not None and np.isfinite(float(v)) else int(default)
    except (TypeError, ValueError):
        return int(default)


def _on(cx, cy, r, deg):
    """Point on a circle at `deg` degrees (for nucleus-edge anchors)."""
    a = np.deg2rad(deg)
    return cx + r * np.cos(a), cy + r * np.sin(a)


def _wave(x0, y0, x1, y1, amp, cyc, n=64):
    """Sine-perturbed polyline from (x0,y0)->(x1,y1) (mRNA transcript / secreted
    ribbon). Perturbation is perpendicular to the segment; returns (X, Y) arrays."""
    t = np.linspace(0.0, 1.0, n)
    dx, dy = x1 - x0, y1 - y0
    L = float(np.hypot(dx, dy)) or 1.0
    nx, ny = -dy / L, dx / L
    off = amp * np.sin(2.0 * np.pi * cyc * t)
    return x0 + dx * t + nx * off, y0 + dy * t + ny * off


def _star_xy(cx, cy, ro, ri, n=5, rot=np.pi / 2.0):
    """Vertices of an n-point star (point up by default) for the convergent core."""
    ang = rot + np.linspace(0.0, 2.0 * np.pi, 2 * n, endpoint=False)
    r = np.empty(2 * n)
    r[0::2] = ro
    r[1::2] = ri
    return np.column_stack([cx + r * np.cos(ang), cy + r * np.sin(ang)])


def _hex(cx, cy, r):
    """Pointy-top hexagon vertices (tissue strip tiles)."""
    a = np.pi / 2.0 + np.linspace(0.0, 2.0 * np.pi, 6, endpoint=False)
    return np.column_stack([cx + r * np.cos(a), cy + r * np.sin(a)])


def _blob(cx, cy, rx, ry, wob=0.05, k=5, phase=0.0, n=160):
    """Organic closed outline (the hepatocyte membrane) -- an ellipse with a small
    sinusoidal wobble so it reads as a cell, not a box. Returns Nx2 vertices."""
    th = np.linspace(0.0, 2.0 * np.pi, n)
    rr = 1.0 + wob * np.sin(k * th + phase) + 0.6 * wob * np.cos((k + 3) * th)
    return np.column_stack([cx + rx * rr * np.cos(th), cy + ry * rr * np.sin(th)])


def _badge(ax, x, y, color, l1, l2, side="right", chip=0.52):
    """A colour-coded count badge = a saturated modality chip (the data mark) + a
    two-line INK label ('physical layer / assay' then the LIVE count) on a white box.
    Colour lives ONLY on the chip; the text stays INK and readable."""
    ax.add_patch(FancyBboxPatch(
        (x, y), chip, chip, boxstyle="round,pad=0.0,rounding_size=0.12",
        facecolor=color, edgecolor="white", lw=0.6, zorder=8, clip_on=False))
    if side == "right":
        tx, ha = x + chip + 0.16, "left"
    else:
        tx, ha = x - 0.16, "right"
    ax.text(tx, y + chip / 2.0, f"{l1}\n{l2}", ha=ha, va="center",
            color=INK, fontsize=FS, zorder=9, clip_on=False, linespacing=1.15,
            bbox=_WBOX)


def _leader(ax, x0, y0, x1, y1):
    """Hairline INK leader from a floated badge to its glyph (used sparingly)."""
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle="-", mutation_scale=1, color=INK,
        lw=0.4, connectionstyle="arc3,rad=0.06", zorder=7, clip_on=False))


# -- main builder -------------------------------------------------------------
def build(d, stats):
    part = d.get("_part", {}) or {}
    conv = (stats or {}).get("conv", {}) or {}

    # ---- LIVE counts (never hard-coded; guarded ints) ----
    universe = _iv(d, "universe")
    coloc    = _iv(d, "coloc")          # DNA layer (genetics)
    substrate = _iv(d, "substrate")     # mRNA layer (transcriptomics)
    A_val    = _iv(d, "A_val")          # open-chromatin gate (snATAC)
    P_liver  = _iv(d, "P_liver_val")    # protein, liver
    P_plasma = _iv(d, "P_plasma_val")   # protein, secreted -> plasma
    S_val    = _iv(d, "S_val")          # tissue frame (spatial)
    ge1  = _iv(part, "ge1")
    ge2  = _iv(part, "ge2")
    all3 = _iv(part, "all3")
    ps_n = _iv(part, "PS")              # convergent core: protein & tissue
    gen_only = _iv(part, "genetic_only")

    # PS enrichment vs the 10k-permutation null (guarded)
    ps = conv.get("PS", {}) or {}
    ps_fold = ps.get("fold", float("nan"))
    ps_p = ps.get("p", float("nan"))
    ps_exp = ps.get("exp", float("nan"))
    fstr = f"{ps_fold:.2f}x" if np.isfinite(ps_fold) else "n/a"
    ps_tag = f"{fstr}  p={_pfmt(ps_p)}"

    # ============================ canvas ============================
    fig, ax = plt.subplots(figsize=(3.9, 3.1))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-0.4, 13.7)
    ax.set_ylim(-0.5, 10.7)

    # ---- faint reading-path guide: top-left (DNA) -> bottom-right (tissue) ----
    ax.add_patch(FancyArrowPatch(
        (1.7, 8.4), (10.0, 1.7), arrowstyle="-|>", mutation_scale=6,
        color=ANAT_LINE, lw=0.8, alpha=0.40,
        connectionstyle="arc3,rad=-0.22", zorder=0.6))

    # ---- ANATOMY (light gray only) --------------------------------------
    ncx, ncy, nr = 3.1, 7.05, 1.35           # nucleus
    # hepatocyte membrane (organic blob)
    ax.add_patch(Polygon(_blob(5.95, 5.55, 5.30, 3.85, wob=0.045, phase=0.7),
                         closed=True, facecolor=ANAT_FILL, edgecolor=ANAT_LINE,
                         lw=1.0, zorder=1.0))
    # sinusoidal vessel lumen, lower-right (plasma compartment)
    ax.add_patch(Ellipse((10.05, 2.85), 2.75, 2.65, facecolor="white",
                         edgecolor=ANAT_LINE, lw=0.9, zorder=1.2))
    ax.text(10.35, 2.35, "vessel\nlumen", ha="center", va="center", color=INK,
            fontsize=FS, zorder=1.5, alpha=0.8)
    # nucleus
    ax.add_patch(Circle((ncx, ncy), nr, facecolor="white", edgecolor=ANAT_LINE,
                        lw=0.9, zorder=1.4))

    # ---- hex-tiled TISSUE strip along the bottom ----
    hr, hy = 0.50, 0.58
    step = np.sqrt(3.0) * hr
    hx = 0.55
    while hx < 13.5:
        ax.add_patch(Polygon(_hex(hx, hy, hr), closed=True, facecolor=ANAT_FILL,
                             edgecolor=ANAT_LINE, lw=0.7, zorder=1.1))
        hx += step

    # ============================ LAYER GLYPHS ============================
    # -- DNA double helix inside the nucleus --
    ylo, yhi = ncy - 1.0, ncy + 1.0
    yv = np.linspace(ylo, yhi, 90)
    amp, wv = 0.46, 2.0 * np.pi / 1.0        # ~2 turns over the helix height
    xa = ncx + amp * np.sin(wv * (yv - ylo))
    xb = ncx - amp * np.sin(wv * (yv - ylo))
    ax.plot(xa, yv, color=ANAT_LINE, lw=1.0, zorder=2.2)
    ax.plot(xb, yv, color=ANAT_LINE, lw=1.0, zorder=2.2)
    for yy in np.linspace(ylo + 0.14, yhi - 0.14, 8):
        xr1 = ncx + amp * np.sin(wv * (yy - ylo))
        xr2 = ncx - amp * np.sin(wv * (yy - ylo))
        ax.plot([xr1, xr2], [yy, yy], color=ANAT_LINE, lw=0.6, zorder=2.1)

    # -- OPEN CHROMATIN gate: nucleosome beads unspooling at the nuclear edge --
    bx0, by0 = _on(ncx, ncy, nr, -50)        # lower-right rim (toward cytoplasm)
    bx1, by1 = 5.35, 5.30
    tb = np.linspace(0.0, 1.0, 4)            # 4 beads, spacing widens -> "opening"
    tb = tb ** 1.35
    fibx, fiby = _wave(bx0, by0, bx1, by1, amp=0.10, cyc=1.5, n=40)
    ax.plot(fibx, fiby, color=ANAT_LINE, lw=0.8, zorder=2.0)
    for t in tb:
        cxb = bx0 + (bx1 - bx0) * t
        cyb = by0 + (by1 - by0) * t
        ax.add_patch(Circle((cxb, cyb), 0.17, facecolor=ANAT_FILL,
                            edgecolor=ANAT_LINE, lw=0.8, zorder=2.3))

    # -- mRNA: wavy transcript through a nuclear pore into the cytoplasm --
    px, py = _on(ncx, ncy, nr, 33)           # pore on the upper-right rim
    ax.add_patch(Circle((px, py), 0.10, facecolor="white", edgecolor=ANAT_LINE,
                        lw=0.9, zorder=1.6))  # pore
    mx, my = _wave(px, py, 6.30, 7.85, amp=0.16, cyc=2.4, n=70)
    ax.plot(mx, my, color=ANAT_LINE, lw=1.0, zorder=2.2)

    # -- PROTEIN: ribosome + folded ribbon in cytoplasm --
    rib_x, rib_y = 6.75, 6.20
    ax.add_patch(Ellipse((rib_x, rib_y - 0.10), 0.62, 0.40, facecolor=ANAT_FILL,
                        edgecolor=ANAT_LINE, lw=0.8, zorder=2.2))   # large subunit
    ax.add_patch(Ellipse((rib_x, rib_y + 0.22), 0.46, 0.26, facecolor=ANAT_FILL,
                        edgecolor=ANAT_LINE, lw=0.8, zorder=2.2))   # small subunit
    # folded ribbon (globular protein): faint blob + a coiled ribbon inside
    pbx, pby = 7.55, 5.55
    ax.add_patch(Ellipse((pbx, pby), 0.80, 0.66, facecolor=ANAT_FILL,
                        edgecolor=ANAT_LINE, lw=0.8, zorder=2.2))
    fx, fy = _wave(pbx - 0.30, pby - 0.12, pbx + 0.30, pby + 0.12, amp=0.18, cyc=2.5, n=44)
    ax.plot(fx, fy, color=ANAT_LINE, lw=0.9, zorder=2.4)

    # -- SECRETED ribbon: crosses the membrane into the vessel lumen --
    sx, sy = _wave(pbx + 0.10, pby - 0.30, 9.55, 3.35, amp=0.16, cyc=2.2, n=70)
    ax.plot(sx, sy, color=ANAT_LINE, lw=1.0, zorder=2.2)
    ax.add_patch(Ellipse((9.55, 3.35), 0.60, 0.50, facecolor=ANAT_FILL,
                        edgecolor=ANAT_LINE, lw=0.8, zorder=2.3))   # secreted protein
    scx, scy = _wave(9.30, 3.28, 9.80, 3.42, amp=0.13, cyc=2.2, n=32)
    ax.plot(scx, scy, color=ANAT_LINE, lw=0.8, zorder=2.4)

    # ---- "prioritized universe" gray bracket over DNA + mRNA (the two origins) ----
    bxa, bxb, byt = 1.85, 7.35, 9.45
    ax.plot([bxa, bxa, bxb, bxb], [byt - 0.20, byt, byt, byt - 0.20],
            color=GRAY, lw=0.8, zorder=3.4, clip_on=False)
    ax.text((bxa + bxb) / 2.0, byt + 0.10, f"prioritized  {universe:,}",
            ha="center", va="bottom", color=INK, fontsize=FS, zorder=9)

    # ============================ CONVERGENT CORE ============================
    # filled STAR where the secreted-protein ribbon meets the tissue strip
    stx, sty = 9.35, 1.38
    ax.add_patch(Polygon(_star_xy(stx, sty, 0.55, 0.22), closed=True,
                         facecolor=C_CONV, edgecolor="white", lw=0.6, zorder=6.0))
    _leader(ax, stx + 0.5, sty + 0.05, 10.02, 1.60)
    ax.text(10.15, 1.55, f"protein + tissue\n{ps_n}   {ps_tag}", ha="left",
            va="center", color=INK, fontsize=FS, zorder=9, linespacing=1.15,
            clip_on=False, bbox=_WBOX)

    # ============================ COUNT BADGES (colour = data) ============================
    # DNA (genetics) -- floated top-left, leader into the helix
    _badge(ax, 0.05, 8.45, C_GWAS, "DNA  genetics", f"{coloc:,}")
    _leader(ax, 0.05 + 0.9, 8.62, 2.70, 7.95)
    # honesty: genetics-alone enriches in nothing (faint GRAY, near the DNA glyph)
    ax.text(0.05, 6.35, f"genetic-only  {gen_only:,}   n.s.", ha="left", va="center",
            color=INK, fontsize=FS, zorder=9, clip_on=False,
            bbox=dict(boxstyle="round,pad=0.18", fc="white", ec=GRAY, lw=0.6, alpha=0.95))

    # open chromatin gate (snATAC) -- below the beads, short leader up
    _badge(ax, 3.55, 4.30, C_ATAC, "open chromatin (gate)", f"snATAC  {A_val:,}")
    _leader(ax, 4.55, 4.85, 4.85, 5.40)

    # mRNA (transcriptomics) -- at the transcript tip
    _badge(ax, 6.45, 7.55, C_BULK, "mRNA  transcriptome", f"{substrate:,}")

    # protein, liver (proteomics) -- by the ribosome/ribbon
    _badge(ax, 7.95, 5.85, C_PROTEO, "protein  proteome", f"liver  {P_liver:,}")
    _leader(ax, 7.95, 6.00, 7.75, 5.75)

    # secreted protein, plasma (proteomics) -- above the vessel
    _badge(ax, 9.75, 4.35, C_PROTEO, "secreted protein", f"plasma  {P_plasma:,}")
    _leader(ax, 9.90, 4.35, 9.65, 3.70)

    # tissue frame (spatial) -- on the left of the hex strip
    _badge(ax, 0.30, 0.30, C_SPATIAL, "tissue (frame)", f"spatial  {S_val:,}")

    # ---- corroboration tally box (honest ladder; INK on white) ----
    # top-right corner, lifted clear of the mRNA badge label that runs beneath it
    # (mRNA label reaches ~x10.4 at y~7.5-8.2; this box sits at y>=9.0 with a margin)
    tlx, tly, tlw, tlh = 9.55, 9.02, 3.95, 1.42
    ax.add_patch(FancyBboxPatch(
        (tlx, tly), tlw, tlh, boxstyle="round,pad=0.06,rounding_size=0.10",
        facecolor="white", edgecolor="#DBDBDB", lw=0.5, zorder=8.5, clip_on=False))
    ax.text(tlx + 0.18, tly + tlh - 0.14,
            f"corroborated\n>=1   {ge1:,}\n>=2    {ge2}\nall 3    {all3}",
            ha="left", va="top", color=INK, fontsize=FS, zorder=9, linespacing=1.25)

    save_pdf(fig, "fig4a_physlayers.pdf")

    # ---- printed caption (full claim + definitions; NOT an on-panel subtitle) ----
    print("[caption:physlayers] A MASLD target prioritized from DNA and mRNA is "
          "corroborated at the protein, chromatin and tissue layers a therapeutic "
          "engages; the protein-and-tissue convergence is the one exceeding chance. "
          "One stylized hepatocyte cross-section (light-gray anatomy) maps each "
          "evidence modality to the PHYSICAL LAYER it measures, read top-left to "
          "bottom-right along biology: DNA (genetics; a double-helix in the nucleus) "
          "prioritizes {coloc:,} causal genes; open chromatin (snATAC; nucleosomes "
          "unspooling at the nuclear edge -- an accessibility GATE, not a "
          "central-dogma rung) validates {A:,}; mRNA (transcriptomics; a transcript "
          "through a nuclear pore) is the {sub:,}-gene substrate -- the gray bracket "
          "marks the DNA-and-mRNA prioritized universe (n={U:,}); protein "
          "(proteomics; a folded ribbon at a ribosome) is recovered in liver "
          "({Pl:,}) with a secreted ribbon crossing into the vessel lumen "
          "detected in plasma ({Pp:,}); and tissue (spatial; the hex strip -- a "
          "spatial FRAME, not a rung) validates {S:,}. The filled star where the "
          "secreted-protein ribbon meets the tissue strip is the convergent core: "
          "{ps} genes validated in BOTH protein and tissue -- the only convergence "
          "above a 10,000-draw label-permutation null (seed 42; {tag}). The tally "
          "box gives the honest corroboration ladder (>=1 assay {ge1:,}; >=2 {ge2}; "
          "all-3 {a3}), and the faint gray badge by the DNA glyph records that "
          "genetics-alone ({go:,} genetic-only) enriches in NOTHING (n.s.). "
          "Saturated colour marks only the count badges and the star; anatomy is "
          "light-gray; every badge prints a live count; snATAC = single-nucleus "
          "ATAC; no gene names.".format(
              coloc=coloc, A=A_val, sub=substrate, U=universe, Pl=P_liver,
              Pp=P_plasma, S=S_val, ps=ps_n, tag=ps_tag, ge1=ge1, ge2=ge2,
              a3=all3, go=gen_only))


if __name__ == "__main__":
    d, stats = get_data()
    build(d, stats)
