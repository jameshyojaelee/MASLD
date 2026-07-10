#!/usr/bin/env python3
"""
Fig 4A v2 -- CANDIDATE 3: CORE-IN-A-FIELD (area-faithful waffle + per-gene zoom).

The claim this geometry welds on (vs the retired display-only Sankey):
  of the prioritized universe, corroboration is a DEFINED ~16% MINORITY that
  distills to a handful of >=2-assay genes and a single tiny above-chance
  proteome-and-spatial core -- the gray:colour AREA RATIO (~84:16) IS the claim.

Two coupled devices (shared data layer = fig4a_v2_common.get_data()):

  MAIN FIELD -- the whole universe as a POPULATION, not a value heatmap:
    * an area-faithful WAFFLE of discrete unit marks with visible white gutters,
      one mark ~= k genes (k chosen so marks fill a 20-wide grid exactly; the
      exact k is printed in the mandatory unit legend "1 square = ~k genes").
    * the vast UNVALIDATED remainder (d["_part"]["unvalidated"]) is drawn in the
      light neutral C_NS so the corroborated minority pops.
    * the CORROBORATED marks (ge1) are clustered as a quarter-disc WEDGE hugging
      the bottom-left corner, filled nearest-corner-first so convergence sits at
      the apex: multi-assay marks (>=2 lenses; the EXCLUSIVE PS/PA/SA/all3) at the
      very corner with a dark INK ring (identity dropped -- convergence encoded by
      ring + position), then single-assay bands proteomics (C_PROTEO) / spatial
      (C_SPATIAL) / snATAC (C_ATAC). Spatial (pink) is banded BETWEEN proteomics
      and snATAC so the two adjacent warm hues never touch.
    * every mark count is apportioned live by largest-remainder from d["_part"], so
      areas are faithful and the marks exactly fill the grid (no ragged remainder).

  ZOOM INSET (per-gene, 1 dot = 1 gene) -- magnify the convergent apex:
    * the ge2 (>=2-assay) genes drawn as individual small dots: the proteome-and-
      spatial core as a blend block, the other pairwise convergences (PA+SA) as
      muted GRAY dots (they sit at chance), and all-3 as larger INK-ringed dots.
    * the PS core carries a dashed GRAY GHOST OUTLINE sized to the permutation-
      EXPECTED count (stats["conv"]["PS"]["exp"]); the dots protruding past the
      ghost are the above-chance EXCESS, coloured C_CONV and annotated live with
      "<fold>x  p=<p>". So chance is a shape on the panel, not a sentence.
    * a labelled leader carries the reader from the main wedge apex to the inset.

House style (fig4a_v2_common sets rcParams -- do NOT re-set): PDF, Helvetica 6, all
text INK, no bold, NO gene names, colour on marks only from the modality palette,
control/absent/unvalidated GRAY/C_NS, ASCII-only on-panel (">=2", "1 square = ~25
genes", "1.36x", "p=2e-3"). Counts read live; expected/fold guarded so the panel
always renders. Full claim in the printed caption, not on-panel.
"""
from fig4a_v2_common import (get_data, save_pdf, INK, FS, GRAY, C_NS, C_CONV,
    C_BULK, C_GWAS, C_PROTEO, C_SPATIAL, C_ATAC, LENS_COL, LENS_NAME,
    tint, fold_ci, _ribbon, plt, np, Rectangle, PathPatch, FancyBboxPatch,
    Circle, Wedge, Ellipse, FancyArrowPatch, Polygon, Path)


# -- small helpers ------------------------------------------------------------
def _pfmt(p):
    """Compact ASCII p-string (Helvetica lacks fancy glyphs): '2e-3', '0.30', 'n/a'."""
    if p is None or not np.isfinite(p):
        return "n/a"
    if p >= 0.01:
        return f"{p:.2f}"
    m, e = f"{p:.0e}".split("e")
    return f"{m}e{int(e)}"


def _mix(a, b, t=0.5):
    """Blend two hex colours (the proteo/spatial blend for the convergent core)."""
    import matplotlib.colors as mc
    ra, ga, ba = mc.to_rgb(a)
    rb, gb, bb = mc.to_rgb(b)
    return (ra + (rb - ra) * t, ga + (gb - ga) * t, ba + (bb - ba) * t)


def _largest_remainder(counts, total):
    """Apportion `total` discrete marks across classes proportional to `counts`
    (a dict name->int gene count), summing EXACTLY to `total` (Hamilton / largest-
    remainder). Guarantees area-faithful marks that fill the grid with no drift."""
    tot_genes = float(sum(max(v, 0) for v in counts.values())) or 1.0
    quota = {k: max(v, 0) * total / tot_genes for k, v in counts.items()}
    floor = {k: int(np.floor(q)) for k, q in quota.items()}
    rem = total - sum(floor.values())
    # hand the leftover marks to the largest fractional parts
    order = sorted(counts, key=lambda k: quota[k] - floor[k], reverse=True)
    for k in order[:max(rem, 0)]:
        floor[k] += 1
    return floor


# -- main builder -------------------------------------------------------------
def build(d, stats):
    part = d["_part"]
    universe = int(d["universe"])
    ge1 = int(part["ge1"])
    ge2 = int(part["ge2"])
    unval = int(part["unvalidated"])

    # ---- exclusive single-assay classes (colour) + multi-assay ring class ----
    # P/S/A_only are disjoint; the multi-assay class = ge2 (PS/PA/SA/all3 pooled),
    # so P_only+S_only+A_only+multi = ge1 and +unval = universe (mutually consistent).
    P_only = int(part["P_only"])
    S_only = int(part["S_only"])
    A_only = int(part["A_only"])
    multi = ge2

    # ---- waffle sizing: k ~ 25 genes/mark, marks fill a 20-wide grid exactly ----
    NCOL = 20
    approx = universe / 25.0
    NROW = max(1, int(round(approx / NCOL)))
    TOTAL = NCOL * NROW                       # marks fill the grid with no blanks
    k_eff = universe / TOTAL                  # genes per mark (~24.7 -> "~25")
    k_disp = int(round(k_eff))

    # apportion every mark live (area-faithful, sums to TOTAL by construction)
    marks_n = _largest_remainder(
        {"M": multi, "P": P_only, "S": S_only, "A": A_only, "U": unval}, TOTAL)

    # ---- build the ordered mark list: convergent apex -> single-assay -> field --
    # Order the coloured classes M, P, S, A first (apex outward); spatial (pink)
    # sits between proteomics and snATAC so the two warm hues never abut. Then the
    # gray unvalidated field. (colour, ring?) per mark.
    CLASS_STYLE = {
        "M": (C_CONV,    True),               # >=2 assays: ringed, identity dropped
        "P": (C_PROTEO,  False),
        "S": (C_SPATIAL, False),
        "A": (C_ATAC,    False),
        "U": (C_NS,      False),              # unvalidated field
    }
    seq = []
    for cls in ("M", "P", "S", "A", "U"):
        seq.extend([CLASS_STYLE[cls]] * marks_n[cls])
    # guard: length must equal the grid; pad/trim defensively (never expected)
    seq = (seq + [CLASS_STYLE["U"]] * TOTAL)[:TOTAL]

    # cells sorted nearest-corner-first -> the coloured prefix forms a corner wedge
    cells = [(r, c) for r in range(NROW) for c in range(NCOL)]
    cells.sort(key=lambda rc: (float(np.hypot(rc[1] + 0.5, rc[0] + 0.5)),
                               rc[0] + rc[1], rc[0], rc[1]))

    # ============================ canvas ============================
    fig, ax = plt.subplots(figsize=(7.0, 3.7))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-8.5, 49.5)
    ax.set_ylim(-6.5, 23.5)

    # ---- MAIN FIELD: draw every unit mark with a visible white gutter ----
    S_MARK, GUT = 0.80, 0.20                  # square size + gutter (pitch = 1.0)
    for (r, c), (col, ring) in zip(cells, seq):
        ax.add_patch(Rectangle(
            (c + GUT / 2.0, r + GUT / 2.0), S_MARK, S_MARK,
            facecolor=col, edgecolor=(INK if ring else "none"),
            lw=(0.6 if ring else 0.0), zorder=4))

    # field frame (hairline) + short title (a label, not a subtitle)
    ax.add_patch(Rectangle((0, 0), NCOL, NROW, facecolor="none",
                           edgecolor=C_NS, lw=0.5, zorder=2))
    ax.text(0, NROW + 0.9, f"prioritized targets   n={universe:,}",
            ha="left", va="bottom", color=INK, fontsize=FS, zorder=9)

    # ---- field callouts: the 84:16 area ratio IS the claim ----
    pct_corr = int(round(100.0 * ge1 / max(universe, 1)))
    pct_unval = 100 - pct_corr
    _wb = dict(boxstyle="round,pad=0.16", fc="white", ec="#E6E6E6", lw=0.3, alpha=0.9)
    # unvalidated: labelled over the gray field (upper-right, away from the wedge)
    ax.text(NCOL * 0.62, NROW * 0.80, f"unvalidated\n{unval:,}  ({pct_unval}%)",
            ha="center", va="center", color=INK, fontsize=FS, zorder=9, bbox=_wb)
    # corroborated: labelled below the coloured corner wedge
    ax.text(NCOL * 0.30, -1.35,
            f"corroborated {ge1:,}  ({pct_corr}%, >=1 assay)",
            ha="center", va="top", color=INK, fontsize=FS, zorder=9)

    # ---- legends (below the field): mandatory unit legend + class key ----
    ly = -3.4
    # unit legend -- one sample square
    ax.add_patch(Rectangle((0.0, ly - S_MARK / 2.0), S_MARK, S_MARK,
                           facecolor=C_NS, edgecolor="none", zorder=4))
    ax.text(1.3, ly, f"1 square = ~{k_disp} genes", ha="left", va="center",
            color=INK, fontsize=FS, zorder=9)
    # class key -- swatches on the row below
    key = [("proteomics", C_PROTEO, False), ("spatial", C_SPATIAL, False),
           ("snATAC", C_ATAC, False), (">=2 assays", C_CONV, True),
           ("unvalidated", C_NS, False)]
    kx, ky = 0.0, ly - 1.55
    for name, col, ring in key:
        ax.add_patch(Rectangle((kx, ky - S_MARK / 2.0), S_MARK, S_MARK,
                               facecolor=col, edgecolor=(INK if ring else "none"),
                               lw=(0.6 if ring else 0.0), zorder=4))
        ax.text(kx + 1.15, ky, name, ha="left", va="center",
                color=INK, fontsize=FS, zorder=9)
        kx += 1.15 + 0.62 * (len(name) + 1)   # advance past the label width

    # ============================ ZOOM INSET ============================
    IX0, IY0, IW, IH = 25.0, 1.0, 24.0, 18.0
    ax.add_patch(FancyBboxPatch(
        (IX0, IY0), IW, IH, boxstyle="round,pad=0.0,rounding_size=0.6",
        facecolor="white", edgecolor=GRAY, lw=0.6, zorder=3))
    ax.text(IX0 + 0.6, IY0 + IH - 0.55, "convergent core   1 dot = 1 gene",
            ha="left", va="top", color=INK, fontsize=FS, zorder=9)

    # live PS stats (guarded) -- exclusive PS count from the partition; exp/fold/p
    # from the permutation null (both refer to the same exclusive proteo&spatial set)
    conv = (stats or {}).get("conv", {}) or {}
    ps = conv.get("PS", {}) or {}
    ps_obs = int(part["PS"])
    ps_exp = ps.get("exp", float("nan"))
    ps_fold = ps.get("fold", float("nan"))
    ps_p = ps.get("p", float("nan"))
    all3 = int(part["all3"])
    other2 = max(ge2 - ps_obs - all3, 0)      # PA+SA -- the at-chance pairs
    exp_r = int(round(ps_exp)) if np.isfinite(ps_exp) and ps_exp > 0 else 0
    excess = max(ps_obs - exp_r, 0) if exp_r else 0
    blend = _mix(C_PROTEO, C_SPATIAL)

    PD = 0.72                                 # dot pitch
    R_DOT = 0.26                              # small-dot radius (per-gene)

    # -- PS core block: choose ncol so the expected count fills whole rows if it can
    #    (clean ghost boundary); else fall back to a fractional expected line.
    pref = [13, 12, 11, 10, 14, 9, 15, 8]
    ncol_ps = next((n for n in pref if exp_r and exp_r % n == 0), 13)
    bx, by = IX0 + 2.0, IY0 + 2.3
    nrow_ps = int(np.ceil(ps_obs / ncol_ps)) if ps_obs else 0
    for i in range(ps_obs):
        r, c = divmod(i, ncol_ps)
        col = C_CONV if (exp_r and i >= exp_r) else blend
        ax.add_patch(Circle((bx + c * PD, by + r * PD), R_DOT,
                            facecolor=col, edgecolor="none", zorder=6))
    # dashed GRAY ghost outline at the permutation-EXPECTED count
    if exp_r:
        clean = (exp_r % ncol_ps == 0)
        g_rows = (exp_r // ncol_ps) if clean else (exp_r / ncol_ps)
        left = bx - 0.5 * PD - 0.10
        right = bx + (ncol_ps - 0.5) * PD + 0.10
        bot = by - 0.5 * PD - 0.10
        top = by + (g_rows - 0.5) * PD + 0.10
        ax.add_patch(Rectangle((left, bot), right - left, top - bot,
                               facecolor="none", edgecolor=GRAY, lw=0.8,
                               ls=(0, (2.4, 1.6)), zorder=7))
        ax.text(right + 0.15, bot + (top - bot) / 2.0,
                f"expected\n{ps_exp:.1f}", ha="left", va="center",
                color=INK, fontsize=FS, zorder=9)
        # excess callout: a short tag to the RIGHT of the protruding teal (above-
        # chance) dots, pointing back at them -- placed low so it stays clear of the
        # label block above and the "other >=2" block to its right.
        if excess > 0 and nrow_ps > 0:
            tgt_r = min(nrow_ps - 1, exp_r // ncol_ps + 1)   # a clearly-above-ghost row
            ax.annotate(f"{excess} above chance",
                        xy=(bx + (ncol_ps - 1) * PD, by + tgt_r * PD),
                        xytext=(right + 0.5, by + (nrow_ps - 1) * PD - 0.3),
                        ha="left", va="center", color=INK, fontsize=FS, zorder=9,
                        arrowprops=dict(arrowstyle="-", color=INK, lw=0.5))

    # PS block label + live fold/p ABOVE the dots, stacked with clear vertical gaps
    # (label on top, stats below) so neither line touches the other or the top dots.
    ps_top = by + (max(nrow_ps, 1)) * PD
    ax.text(bx - 0.5, ps_top + 1.55, "proteome and spatial", ha="left", va="bottom",
            color=INK, fontsize=FS, zorder=9)
    fstr = f"{ps_fold:.2f}x" if np.isfinite(ps_fold) else "n/a"
    ax.text(bx - 0.5, ps_top + 0.55,
            f"{ps_obs} vs {ps_exp:.1f}   {fstr}  p={_pfmt(ps_p)}"
            if np.isfinite(ps_exp) else f"{ps_obs} observed",
            ha="left", va="bottom", color=INK, fontsize=FS, zorder=9)

    # -- other pairwise convergences (PA+SA): muted gray dots, at chance --
    ncol_o = 6
    ox = bx + ncol_ps * PD + 6.5       # pushed right so the excess tag has clear space
    oy = by
    for i in range(other2):
        r, c = divmod(i, ncol_o)
        ax.add_patch(Circle((ox + c * PD, oy + r * PD), R_DOT,
                            facecolor=GRAY, edgecolor="none", zorder=6))
    if other2:
        o_rows = int(np.ceil(other2 / ncol_o))
        ax.text(ox + (ncol_o - 1) * PD / 2.0, oy + o_rows * PD + 0.25,
                f"other >=2\n(at chance) {other2}", ha="center", va="bottom",
                color=INK, fontsize=FS, zorder=9)

    # -- all-3 convergence: larger INK-ringed dots at the inset apex --
    a3y = IY0 + IH - 2.4
    a3x0 = IX0 + 2.2
    for i in range(all3):
        ax.add_patch(Circle((a3x0 + i * 1.5, a3y), 0.46, facecolor=C_CONV,
                            edgecolor=INK, lw=0.8, zorder=6))
    if all3:
        ax.text(a3x0 - 0.9, a3y, "all 3", ha="right", va="center",
                color=INK, fontsize=FS, zorder=9)

    # ---- labelled leader: main wedge apex -> inset (the per-gene zoom) ----
    ax.add_patch(FancyArrowPatch(
        (1.6, 1.6), (IX0 - 0.3, IY0 + IH * 0.45), arrowstyle="-|>",
        mutation_scale=6, color=INK, lw=0.6,
        connectionstyle="arc3,rad=-0.16", zorder=8))
    ax.text(NCOL + (IX0 - NCOL) / 2.0, IY0 + IH * 0.45 + 2.4, f"{k_disp}x per-gene",
            ha="center", va="bottom", color=INK, fontsize=FS, zorder=9,
            bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="#E6E6E6",
                      lw=0.3, alpha=0.9))

    save_pdf(fig, "fig4a_corefield.pdf")

    # ---- printed caption (full claim + definitions; NOT an on-panel subtitle) ----
    print("[caption:corefield] Corroboration is a defined minority of the prioritized "
          "universe. MAIN: an area-faithful waffle of the n={U:,} prioritized targets, "
          "one square ~= {k} genes (marks apportioned live by largest-remainder so areas "
          "are faithful and fill the 20x{nr} grid). The light-gray field is the "
          "unvalidated remainder ({uv:,}, {pu}%); the coloured quarter-disc wedge in the "
          "corner is the corroborated minority ({g1:,}, {pc}%, validated in >=1 orthogonal "
          "assay), filled nearest-corner-first so the >=2-assay marks (dark ring; identity "
          "dropped) sit at the apex, then single-assay bands proteomics ({Po}), spatial "
          "({So}) and snATAC ({Ao}) genes -- spatial banded between the two warm hues so "
          "they do not abut. INSET (1 dot = 1 gene): the {g2} >=2-assay genes -- the "
          "proteome-and-spatial core as a blend block against a gray dashed GHOST at its "
          "10,000-draw permutation-EXPECTED count (seed 42), the other pairwise "
          "convergences (proteome-and-snATAC + spatial-and-snATAC, {o2}) as gray dots at "
          "chance, and all-3 ({a3}) as ringed dots. Of {U:,} prioritized targets, "
          "corroboration is a defined ~{pc}% minority ({g1:,} in >=1 assay) distilling to "
          "{g2} in >=2 and a tiny above-chance proteome-and-spatial core "
          "(~{po} vs ~{pe:.1f}, {pf}x, p={pp}) -- the ~{pu}:{pc} gray:colour area ratio is "
          "the claim. snATAC = single-nucleus ATAC; every mark encodes a live count; no "
          "gene names.".format(
              U=universe, k=k_disp, nr=NROW, uv=unval, pu=pct_unval, g1=ge1, pc=pct_corr,
              Po=P_only, So=S_only, Ao=A_only, g2=ge2, o2=other2, a3=all3,
              po=ps_obs, pe=(ps_exp if np.isfinite(ps_exp) else float("nan")),
              pf=(f"{ps_fold:.2f}" if np.isfinite(ps_fold) else "n/a"),
              pp=_pfmt(ps_p)))


if __name__ == "__main__":
    d, stats = get_data()
    build(d, stats)
