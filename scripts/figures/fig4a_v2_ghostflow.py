#!/usr/bin/env python3
"""
Fig 4A v2 -- CANDIDATE 2: NULL-GHOST ALLUVIAL (the redeemed cascade).

The claim this geometry welds on (vs the retired display-only Sankey):
  orthogonal validation exceeds chance ONLY where two independent tissue assays
  agree -- proteome and spatial. The other pairwise convergences and the
  three-way convergence sit AT OR BELOW chance.

The device: keep the Fig-2A cubic-Bezier alluvial idiom the PI approved, then make
"beats chance" a VISIBLE EVENT by DRAWING CHANCE. Three stages, left -> right, all
built with the shared _ribbon():

  STAGE 0  origins  : two stacked nodes -- transcriptomics (C_BULK, height propto
                      d["substrate"]) + genetics (C_GWAS, height propto d["coloc"]);
                      an inset neutral seam encodes the shared overlap (coloc_shared).
  STAGE 1  lenses   : three nodes -- proteomics / spatial / snATAC (heights propto
                      P_val / S_val / A_val); ribbons from both origins fan into them.
  STAGE 2  outcomes : the GHOST DEVICE. Each pairwise convergence (proteome and
                      spatial, proteome and snATAC, spatial and snATAC) plus a small
                      all-three is a horizontal bar (length = observed genes) fed by
                      ribbons from its parent lanes, drawn against a GRAY (#9E9E9E)
                      DASHED GHOST rectangle at the UPPER 95% BOUND of its permutation
                      null -- the chance ENVELOPE, not the mean (seed 42, 10k-draw null).
                      An outcome beats chance only when its bar OVERFLOWS that 95% ghost
                      (equivalently p<0.05). proteome-and-spatial is the one overflow
                      (obs > 95% bound): its core is a proteo/spatial blend and the excess
                      beyond the ghost line is a saturated C_CONV "overflow cap",
                      annotated "1.36x  p=2e-3" (live). The others nest INSIDE their 95%
                      ghosts -- so chance is a shape on the panel, not a caption sentence.

House style (fig4a_v2_common sets rcParams -- do NOT re-set): PDF, Helvetica 6, all
text INK, no bold, NO gene names, colour on marks only from the modality palette,
control/absent GRAY, ASCII-only on-panel ("proteome and spatial", "1.36x", "p=2e-3").
Expected sizes are guarded against zero/NaN so the panel always renders. Full claim in
the printed caption, not on-panel.
"""
from fig4a_v2_common import (get_data, save_pdf, INK, FS, GRAY, C_NS, C_CONV,
    C_BULK, C_GWAS, C_PROTEO, C_SPATIAL, C_ATAC, LENS_COL, LENS_NAME,
    tint, fold_ci, _ribbon, plt, np, Rectangle, PathPatch, FancyBboxPatch,
    Circle, Wedge, Ellipse, FancyArrowPatch, Polygon, Path)

# Contract-sanctioned neutral gray for the shared-origin seam + at/below-chance bars.
NEUTRAL = "#B9BEC4"
# white-backed label box so on-flow text lifts off the ribbons (Fig-2A idiom).
_WBOX = dict(boxstyle="round,pad=0.12", fc="white", ec="#D9D9D9", lw=0.3, alpha=0.9)


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
    """Blend two hex colours (proteo/spatial blend for the overflowing core)."""
    import matplotlib.colors as mc
    ra, ga, ba = mc.to_rgb(a)
    rb, gb, bb = mc.to_rgb(b)
    return (ra + (rb - ra) * t, ga + (gb - ga) * t, ba + (bb - ba) * t)


def _stack_top(keys, wt, height, gap, minh):
    """Proportional slots filling `height` from the TOP down (first key at top),
    each >= minh. Returns {key: (y_low, y_high)} with y increasing upward."""
    n = len(keys)
    avail = height - gap * (n - 1)
    w = {k: max(float(wt[k]), 1e-9) for k in keys}
    h = {k: avail * w[k] / sum(w.values()) for k in keys}
    for _ in range(6):
        lo = [k for k in keys if h[k] < minh]
        if not lo:
            break
        free = [k for k in keys if k not in lo]
        fw = sum(w[k] for k in free)
        rem = avail - minh * len(lo)
        for k in lo:
            h[k] = minh
        for k in free:
            h[k] = (rem * w[k] / fw) if fw > 0 else rem / max(len(free), 1)
    slots = {}
    y = height
    for k in keys:
        slots[k] = (y - h[k], y)
        y = y - h[k] - gap
    return slots


# -- main builder -------------------------------------------------------------
def build(d, stats):
    part = d["_part"]
    conv = (stats or {}).get("conv", {}) or {}

    # ---- live scalars (never hard-coded) ----
    sub = int(d["substrate"])          # transcriptomics origin
    col = int(d["coloc"])              # genetics origin
    ovl = int(d.get("coloc_shared", part.get("overlap", 0)))
    Pv, Sv, Av = int(d["P_val"]), int(d["S_val"]), int(d["A_val"])

    def gv(key, fb):
        """observed / expected / fold / p for a convergence outcome, guarded."""
        c = conv.get(key, {})
        o = int(c.get("obs", fb))
        try:
            e = float(c.get("exp", float("nan")))
        except (TypeError, ValueError):
            e = float("nan")
        f = c.get("fold", float("nan"))
        p = c.get("p", float("nan"))
        return o, e, f, p

    # outcomes: (conv key, parent lanes top->bottom on the bar face, display name)
    OUT = [
        ("PS",   ["P", "S"],      "proteome\nand spatial"),
        ("PA",   ["P", "A"],      "proteome\nand snATAC"),
        ("SA",   ["S", "A"],      "spatial\nand snATAC"),
        ("all3", ["P", "S", "A"], "all three"),
    ]
    obs = {k: gv(k, int(part.get(k, 0)))[0] for k, _, _ in OUT}
    exp = {k: gv(k, int(part.get(k, 0)))[1] for k, _, _ in OUT}
    fld = {k: gv(k, int(part.get(k, 0)))[2] for k, _, _ in OUT}
    pvl = {k: gv(k, int(part.get(k, 0)))[3] for k, _, _ in OUT}
    # upper 95% bound of the permutation null (the CHANCE ENVELOPE the ghost draws)
    hiv = {}
    for k, _, _ in OUT:
        c = conv.get(k, {})
        try:
            hiv[k] = float(c.get("hi", float("nan")))
        except (TypeError, ValueError):
            hiv[k] = float("nan")

    # ============================ geometry ============================
    H = 10.0
    BW = 0.14
    X0, X1, X2 = 0.0, 1.90, 4.70      # origins / lenses / outcome baseline
    LANES = ["P", "S", "A"]
    LVAL = {"P": Pv, "S": Sv, "A": Av}
    LTOT = float(sum(LVAL.values())) or 1.0

    # --- STAGE 0 origins (heights strictly proportional; no floor) ---
    o_gap = 0.30
    scale0 = (H - o_gap) / max(sub + col, 1)
    h_tx = sub * scale0
    h_ge = col * scale0
    tx_lo, tx_hi = H - h_tx, H            # transcriptomics on top
    ge_hi = tx_lo - o_gap
    ge_lo = ge_hi - h_ge                  # genetics below (~0)
    tx_c = (tx_lo + tx_hi) / 2.0
    ge_c = (ge_lo + ge_hi) / 2.0
    seam_c = (tx_lo + ge_hi) / 2.0        # boundary between the two origins
    seam_h = ovl * scale0                 # overlap encoded as seam height

    # --- STAGE 1 lens slots (heights ~ proportional; tiny floor) ---
    slotL = _stack_top(LANES, LVAL, H, gap=0.80, minh=0.80)
    laneH = {l: slotL[l][1] - slotL[l][0] for l in LANES}
    f_tx = sub / max(sub + col, 1)        # per-lane origin split (both feed all)
    f_ge = col / max(sub + col, 1)

    # --- STAGE 2 outcome bars (length = observed count; ghost = expected) ---
    # A generous ROW_PITCH gives each outcome a clearly separated horizontal band so
    # its 2-line label + bar + count never touch a neighbour (coordinate check below).
    # ghost reaches the 95% bound (hi), so the x-scale must fit max(obs, hi) per row.
    def _rowmax(k):
        ref = (hiv[k] if np.isfinite(hiv[k])
               else (exp[k] if np.isfinite(exp[k]) else obs[k]))
        return max(obs[k], ref)
    draw_max = max([_rowmax(k) for k, _, _ in OUT] + [1])
    BAR_LMAX = 2.30
    scaleL = BAR_LMAX / draw_max
    MINL = 0.15
    BAR_TH = 0.58
    ROW_PITCH = 2.10                        # vertical distance between outcome rows
    n_out = len(OUT)
    block = (n_out - 1) * ROW_PITCH
    top_c = H / 2.0 + block / 2.0
    center = {OUT[i][0]: top_c - i * ROW_PITCH for i in range(n_out)}

    # "beats chance" = observed exceeds the UPPER 95% NULL BOUND (hi) AND the permutation
    # p clears 0.05. ONLY such a row (proteome-and-spatial) CROSSES its ghost and is
    # coloured + capped + annotated; every other row NESTS inside its 95% ghost and is
    # muted/gray with its observed COUNT ONLY -- a marginal obs>mean that is not
    # significant (e.g. proteome-and-snATAC, fold ~1.1, p ~0.15) still sits within its
    # envelope. Lengths are never distorted: bar = true observed, ghost = true 95% bound.
    geo = {}
    for k, _, _ in OUT:
        o = obs[k]
        e = exp[k]
        p = pvl[k]
        # ghost width = upper 95% bound of the null (chance envelope); fall back to the
        # mean, then observed, if hi is missing/NaN.
        ceil = (hiv[k] if np.isfinite(hiv[k]) and hiv[k] > 0
                else (e if np.isfinite(e) and e > 0 else o))
        beats = bool(np.isfinite(ceil) and o > ceil and np.isfinite(p) and p < 0.05)
        L_obs = max(o * scaleL, MINL)
        L_gh = max(ceil * scaleL, MINL)
        geo[k] = dict(L_obs=L_obs, L_gh=L_gh, beats=beats, o=o, e=e, p=p, ceil=ceil)

    # ============================ canvas ============================
    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    ax.set_xlim(-2.7, 8.3)
    ax.set_ylim(-0.7, 11.0)
    ax.axis("off")

    def node(x, y0, y1, c, z=5):
        ax.add_patch(Rectangle((x - BW / 2, y0), BW, y1 - y0, facecolor=c,
                               edgecolor="white", lw=0.2, zorder=z))

    # ---- STAGE 0 -> STAGE 1 ribbons: both origins fan into every lens ----
    # transcriptomics origin -> each lens's TOP portion (share ~ lens size)
    y = tx_hi
    for l in LANES:
        w = h_tx * LVAL[l] / LTOT
        ltop = slotL[l][1]
        lbot = ltop - laneH[l] * f_tx
        _ribbon(ax, X0 + BW / 2, y - w, y, X1 - BW / 2, lbot, ltop,
                tint(C_BULK, 0.18), alpha=0.22)
        y -= w
    # genetics origin -> each lens's BOTTOM portion
    y = ge_hi
    for l in LANES:
        w = h_ge * LVAL[l] / LTOT
        lbot = slotL[l][0]
        ltop = lbot + laneH[l] * f_ge
        _ribbon(ax, X0 + BW / 2, y - w, y, X1 - BW / 2, lbot, ltop,
                tint(C_GWAS, 0.18), alpha=0.20)
        y -= w

    # ---- STAGE 1 -> STAGE 2 ribbons: each lens sends its convergent slivers ----
    # A lens feeds only the outcomes it belongs to; slices ~ obs on the lens scale
    # (so the small convergent fraction reads as small); the rest falls away.
    lane_feeds = {"P": ["PS", "PA", "all3"],
                  "S": ["PS", "SA", "all3"],
                  "A": ["PA", "SA", "all3"]}
    parents_of = {k: par for k, par, _ in OUT}
    for l in LANES:
        s1 = laneH[l] / max(LVAL[l], 1)     # lens-local gene -> unit scale
        y = slotL[l][1]                     # anchor convergent slivers at lens top
        for k in lane_feeds[l]:
            w = max(obs[k] * s1, 0.012)
            par = parents_of[k]
            i = par.index(l)
            sub_h = BAR_TH / len(par)
            fy_hi = center[k] + BAR_TH / 2.0
            t_hi = fy_hi - i * sub_h
            t_lo = fy_hi - (i + 1) * sub_h
            _ribbon(ax, X1 + BW / 2, y - w, y, X2, t_lo, t_hi,
                    LENS_COL[l], alpha=0.26)
            y -= w

    # ---- STAGE 0 origin nodes + seam + labels ----
    node(X0, tx_lo, tx_hi, C_BULK)
    node(X0, ge_lo, ge_hi, C_GWAS)
    # shared-overlap seam: inset neutral band straddling the two origins
    ax.add_patch(Rectangle((X0 - BW * 0.28, seam_c - seam_h / 2.0), BW * 0.56, seam_h,
                           facecolor=NEUTRAL, edgecolor="white", lw=0.3, zorder=5.5))
    lx = X0 - BW / 2 - 0.16
    ax.text(lx, tx_c, f"transcriptomics\n{sub:,}", ha="right", va="center",
            color=INK, zorder=9)
    ax.text(lx, ge_c, f"genetics\n{col:,}", ha="right", va="center",
            color=INK, zorder=9)
    ax.text(lx, seam_c, f"shared {ovl:,}", ha="right", va="center",
            color=INK, zorder=9, bbox=_WBOX)

    # ---- STAGE 1 lens nodes + labels (centred above each node) ----
    for l in LANES:
        y0, y1 = slotL[l]
        node(X1, y0, y1, LENS_COL[l])
        ax.text(X1, y1 + 0.08, f"{LENS_NAME[l]} {LVAL[l]:,}", ha="center",
                va="bottom", color=INK, zorder=9, bbox=_WBOX)

    # ---- STAGE 2 ghost device: observed bar vs dashed expected ghost ----
    # Per-row layout (no two text boxes overlap -- coordinates reasoned in the report):
    #   LEFT  : 2-line outcome label, RIGHT-anchored at LAB_R, vertically CENTRED.
    #   MIDDLE: faint expected body + observed bar + dashed gray ghost (= expected).
    #   RIGHT : observed count immediately right of the bar, vertically centred.
    #   ABOVE : fold/p on its OWN line, ONLY for the beats-chance row, clear of the bar.
    LAB_R = X2 - 0.20
    for k, par, name in OUT:
        g = geo[k]
        yc = center[k]
        y_lo, y_hi = yc - BAR_TH / 2.0, yc + BAR_TH / 2.0
        L_obs, L_gh, beats = g["L_obs"], g["L_gh"], g["beats"]
        o = g["o"]
        # faint "expected volume" body behind the ghost outline
        ax.add_patch(Rectangle((X2, y_lo), L_gh, BAR_TH, facecolor=C_NS,
                               edgecolor="none", alpha=0.55, zorder=2.4))
        # observed bar
        if beats:
            core = _mix(C_PROTEO, C_SPATIAL)
            ax.add_patch(Rectangle((X2, y_lo), L_gh, BAR_TH, facecolor=core,
                                   edgecolor="none", alpha=0.95, zorder=5))
            # saturated cap = the part BEYOND the 95% envelope (observed - hi)
            cap_w = max(L_obs - L_gh, 0.0)
            if cap_w > 1e-6:
                ax.add_patch(Rectangle((X2 + L_gh, y_lo), cap_w, BAR_TH,
                                       facecolor=C_CONV, edgecolor="white", lw=0.4,
                                       zorder=7))
        else:
            ax.add_patch(Rectangle((X2, y_lo), L_obs, BAR_TH, facecolor=NEUTRAL,
                                   edgecolor="none", alpha=0.90, zorder=5))
        # GRAY dashed ghost outline (drawn on top so the chance line stays visible)
        ax.add_patch(Rectangle((X2, y_lo), L_gh, BAR_TH, facecolor="none",
                               edgecolor=GRAY, lw=0.8, ls=(0, (2.2, 1.6)), zorder=6))
        # LEFT: 2-line outcome label, right-anchored, vertically centred on the row
        ax.text(LAB_R, yc, name, ha="right", va="center", color=INK,
                zorder=9, bbox=_WBOX)
        # RIGHT: observed count immediately right of the mark (bar or, when the bar
        # nests, its wider 95% ghost) so the number never sits inside the ghost box.
        ax.text(X2 + max(L_obs, L_gh) + 0.10, yc, f"{o:,}", ha="left", va="center",
                color=INK, zorder=9)
        # ABOVE (beats row only): fold/p on its own line, clear of the bar + count
        if beats:
            cap_mid = X2 + L_gh + (L_obs - L_gh) / 2.0
            f = fld[k]
            fstr = f"{f:.2f}x" if np.isfinite(f) else "n/a"
            ax.text(cap_mid, y_hi + 0.34, f"{fstr}  p={_pfmt(pvl[k])}",
                    ha="center", va="bottom", color=INK, zorder=9)

    # ---- compact legend teaching the ghost device (top of the outcome column,
    #      above the top row's fold/p line -> no collision) ----
    lgx = X2
    for lgy, kind, txt in ((10.50, "ghost", "chance (95% null)"),
                           (9.85, "solid", "observed")):
        if kind == "ghost":
            ax.add_patch(Rectangle((lgx, lgy - 0.13), 0.46, 0.26, facecolor="none",
                                   edgecolor=GRAY, lw=0.8, ls=(0, (2.2, 1.6)),
                                   zorder=6))
        else:
            ax.add_patch(Rectangle((lgx, lgy - 0.13), 0.46, 0.26, facecolor=C_CONV,
                                   edgecolor="none", zorder=6))
        ax.text(lgx + 0.56, lgy, txt, ha="left", va="center", color=INK, zorder=9)

    # ---- subtle stage tags along the bottom ----
    for xtag, txt in ((X0, "origins"), (X1, "validation lenses"),
                      (X2 + BAR_LMAX / 2.0, "convergence vs chance")):
        ax.text(xtag, -0.42, txt, ha="center", va="top", color=INK, zorder=9)

    save_pdf(fig, "fig4a_ghostflow.pdf")

    # ---- printed caption (full claim + definitions; NOT an on-panel subtitle) ----
    ps_o, ps_e, ps_f, ps_p = gv("PS", part.get("PS", 0))
    ps_hi = hiv["PS"]

    def _fx(v):
        return f"{v:.2f}x" if np.isfinite(v) else "n/a"

    print("[caption:ghostflow] Orthogonal validation exceeds chance only where two "
          "independent tissue assays agree. Left: the prioritized universe splits into "
          "transcriptomic {sub:,} and genetic {col:,} candidates ({ovl:,} shared, drawn "
          "as the inset seam). Middle: three orthogonal lenses validate proteomics "
          "{Pv:,}, spatial {Sv:,} and snATAC {Av:,} genes (heights on each stage's own "
          "scale; ribbons carry only the convergent slivers -- non-convergent genes fall "
          "away). Right: each pairwise convergence outcome is an OBSERVED bar (length = "
          "genes) against a gray dashed GHOST drawn at the UPPER 95% BOUND of its "
          "permutation null -- the chance ENVELOPE, not the mean (seed 42, 10,000-draw "
          "label-permutation null over the co-measured universe). A bar that OVERFLOWS "
          "its 95% ghost beats chance (equivalently p<0.05) and is the only row coloured, "
          "capped and annotated. proteome-and-spatial is that row ({ps_o} observed vs "
          "{ps_e:.1f} expected, 95% null upper {ps_hi:.1f}, {ps_f:.2f}x, p={ps_p}); the "
          "saturated cap is the excess beyond the 95% envelope (observed minus upper "
          "bound). The other rows NEST inside their 95% ghosts -- none is significant, so "
          "each shows its observed count only: proteome-and-snATAC {pa_o} (fold {pa_f}, "
          "p={pa_p}), spatial-and-snATAC {sa_o} (fold {sa_f}, p={sa_p}), all-three {a3_o} "
          "(fold {a3_f}, p={a3_p}). A marginal observed-above-mean that stays within the "
          "95% envelope is left uncalled (no cap, no fold/p). Direction is annotated, not "
          "filtered; snATAC = single-nucleus ATAC; every mark encodes a live count; no "
          "gene names.".format(
              sub=sub, col=col, ovl=ovl, Pv=Pv, Sv=Sv, Av=Av,
              ps_o=int(ps_o), ps_e=float(ps_e) if np.isfinite(ps_e) else float("nan"),
              ps_hi=float(ps_hi) if np.isfinite(ps_hi) else float("nan"),
              ps_f=float(ps_f) if np.isfinite(ps_f) else float("nan"), ps_p=_pfmt(ps_p),
              pa_o=int(obs["PA"]), pa_f=_fx(fld["PA"]), pa_p=_pfmt(pvl["PA"]),
              sa_o=int(obs["SA"]), sa_f=_fx(fld["SA"]), sa_p=_pfmt(pvl["SA"]),
              a3_o=int(obs["all3"]), a3_f=_fx(fld["all3"]), a3_p=_pfmt(pvl["all3"])))


if __name__ == "__main__":
    d, stats = get_data()
    build(d, stats)
