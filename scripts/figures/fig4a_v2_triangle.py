#!/usr/bin/env python3
"""
Fig 4A v2 -- CANDIDATE 5: COUNT-EXACT VALIDATION TRIANGLE (node-link).

The claim this geometry welds on (vs the retired display-only Sankey):
  the three orthogonal validation assays share a convergent core, but only ONE
  pairwise bridge -- proteome-and-spatial -- exceeds chance; the other two pairwise
  cores and the three-way core sit at or below chance.

The device: a fixed, semantic node-link diagram (NO spring layout) where every node
AREA and every edge WIDTH encodes a live count, so the picture is quantitatively
auditable rather than decorative (an earlier area-decorative "wheel" was rejected).

  LEFT -- prioritization funnel:
    two input nodes -- transcriptomics (C_BULK, area propto d["substrate"]) and
    genetics (C_GWAS, area propto d["coloc"]) -- funnel via tapered ribbons into ONE
    neutral UNION node (GRAY, area propto d["universe"], labelled "prioritized N").
    Node area propto count => radius propto sqrt(count); a size legend gives two
    reference counts. The two streams overlap (d["coloc_shared"]) -- a thin neutral
    seam at the merge encodes that shared count.

  RIGHT -- validation triangle:
    three hubs on a FIXED equilateral triangle -- proteomics (C_PROTEO, area propto
    d["P_val"]), spatial (C_SPATIAL, d["S_val"]), snATAC (C_ATAC, d["A_val"]) -- with
    snATAC at the apex nearest the union so the proteome-and-spatial edge is the far
    (right) side, uncrossed by the incoming funnels. Funnel EDGES from the union to
    each hub have width propto that hub's validated count (714/447/560).

  THE CLAIM -- three internal triangle edges = pairwise shared-gene cores, width
    propto the EXCLUSIVE pairwise counts d["_part"] (PS / PA / SA). The proteome-and-
    spatial edge is drawn THICK + SOLID (light "expected" core + a saturated excess
    margin on each side) against a thin GRAY (#9E9E9E) DASHED GHOST outline at the
    permutation-EXPECTED width (stats["conv"]["PS"]["exp"]) -- so the observed band
    visibly OVERSHOOTS its ghost on both sides. The other two internal edges are
    faded/hairline and NEST INSIDE their dashed expected envelopes (at/below chance).
    A small central cluster = the all-3 core at the centroid (1 mark = 1 gene). Only
    the proteome-and-spatial edge is annotated ("<fold>x  p=<p>", live). "Above
    chance" is thus triple-encoded -- width + solidity + ghost-overshoot -- so it
    survives grayscale and colour-blind reproduction.

House style (fig4a_v2_common sets rcParams -- do NOT re-set): PDF, Helvetica 6, all
text INK, no bold, NO gene names, colour on marks only from the modality palette,
control/absent GRAY, ASCII-only on-panel ("proteomics", "spatial", "snATAC",
"proteome and spatial", "1.36x", "p=2e-3", "prioritized 9,882"). Counts read live;
expected/fold guarded against zero/NaN so the panel always renders. Full claim in the
printed caption, not on-panel.
"""
from fig4a_v2_common import (get_data, save_pdf, INK, FS, GRAY, C_NS, C_CONV,
    C_BULK, C_GWAS, C_PROTEO, C_SPATIAL, C_ATAC, LENS_COL, LENS_NAME,
    tint, fold_ci, _ribbon, plt, np, Rectangle, PathPatch, FancyBboxPatch,
    Circle, Wedge, Ellipse, FancyArrowPatch, Polygon, Path)

# Contract-sanctioned neutral gray for the shared-origin seam + at/below-chance bands.
NEUTRAL = "#B9BEC4"
# white-backed label box so on-flow text lifts off the ribbons (Fig-2A idiom).
_WBOX = dict(boxstyle="round,pad=0.10", fc="white", ec="#D9D9D9", lw=0.3, alpha=0.9)


# -- small helpers ------------------------------------------------------------
def _pfmt(p):
    """Compact ASCII p-string (Helvetica lacks fancy glyphs): '2e-3', '0.30', 'n/a'."""
    if p is None or not np.isfinite(p):
        return "n/a"
    if p >= 0.01:
        return f"{p:.2f}"
    m, e = f"{p:.0e}".split("e")
    return f"{m}e{int(e)}"


def _fin(v, fb=float("nan")):
    """Coerce to float; return fallback on None/non-numeric (never raises)."""
    try:
        f = float(v)
        return f if np.isfinite(f) else fb
    except (TypeError, ValueError):
        return fb


def _seg(A, B, rA, rB):
    """Trim segment A->B by the endpoint node radii and return the trimmed
    endpoints A2,B2 plus the unit PERPENDICULAR (vx,vy). Bands are built by
    offsetting A2,B2 along +/- perp, so an edge connects node RIM to node RIM."""
    ax_, ay_ = A
    bx_, by_ = B
    dx, dy = bx_ - ax_, by_ - ay_
    L = float(np.hypot(dx, dy)) or 1e-9
    ux, uy = dx / L, dy / L
    A2 = (ax_ + ux * rA, ay_ + uy * rA)
    B2 = (bx_ - ux * rB, by_ - uy * rB)
    return A2, B2, (-uy, ux)


def _band_poly(A2, B2, perp, hw):
    """Rectangle corners for a band of half-width hw centred on segment A2-B2."""
    vx, vy = perp
    return [(A2[0] + vx * hw, A2[1] + vy * hw), (B2[0] + vx * hw, B2[1] + vy * hw),
            (B2[0] - vx * hw, B2[1] - vy * hw), (A2[0] - vx * hw, A2[1] - vy * hw)]


def _strip_poly(A2, B2, perp, h1, h2):
    """Quad between signed offsets h1 and h2 on one side (the above-chance excess)."""
    vx, vy = perp
    return [(A2[0] + vx * h1, A2[1] + vy * h1), (B2[0] + vx * h1, B2[1] + vy * h1),
            (B2[0] + vx * h2, B2[1] + vy * h2), (A2[0] + vx * h2, A2[1] + vy * h2)]


# -- main builder -------------------------------------------------------------
def build(d, stats):
    part = d["_part"]
    conv = (stats or {}).get("conv", {}) or {}

    # ---- live scalars (never hard-coded) ----
    uni = int(d["universe"])
    sub = int(d["substrate"])                          # transcriptomics origin
    col = int(d["coloc"])                              # genetics origin
    shared = int(d.get("coloc_shared", part.get("overlap", 0)))
    Pv, Sv, Av = int(d["P_val"]), int(d["S_val"]), int(d["A_val"])

    def gv(key):
        """observed(exclusive) / expected / fold / p for a convergence edge, guarded.
        obs falls back to the exclusive partition; exp/fold/p from the seeded null."""
        c = conv.get(key, {}) or {}
        o = int(c.get("obs", part.get(key, 0)))
        return o, _fin(c.get("exp")), _fin(c.get("fold")), c.get("p", float("nan"))

    PS_o, PS_e, PS_f, PS_p = gv("PS")                  # the one above-chance bridge
    PA_o, PA_e, _PA_f, _PA_p = gv("PA")                # at/below chance
    SA_o, SA_e, _SA_f, _SA_p = gv("SA")                # at/below chance
    all3 = int(part.get("all3", d.get("all3", 0)))

    # ============================ scales =============================
    # NODE AREA channel: r = R_SCALE * sqrt(count) (one scale for all six nodes;
    # taught by the size legend). Union is largest, hubs smallest -- faithful.
    TARGET_UNION_R = 1.50
    R_SCALE = TARGET_UNION_R / float(np.sqrt(max(uni, 1)))

    def rad(n):
        return max(R_SCALE * float(np.sqrt(max(n, 0))), 0.04)

    r_uni = rad(uni)
    r_tx, r_ge = rad(sub), rad(col)
    r_P, r_S, r_A = rad(Pv), rad(Sv), rad(Av)

    # EDGE WIDTH channel is LINEAR in count. The dynamic range (hundreds for
    # validation flows vs tens for shared cores) is >10x, so a single linear map
    # would render the crucial ~70-vs-51.5 overshoot invisible OR blow the funnels
    # off-canvas. We therefore use TWO internally-proportional linear scales -- one
    # for the three validation funnels, one for the three shared-core edges -- and
    # print the exact count on/beside every edge so width is never the sole carrier
    # (the size legend covers node area; the caption states both scales).
    w_fun = 0.90 / max(Pv, Sv, Av, 1)                  # validation-funnel width/gene
    w_con = 0.85 / max(PS_o, PA_o, SA_o, 1)            # shared-core-edge width/gene

    # ============================ geometry ===========================
    # Inputs stacked left; union mid-left; equilateral triangle (apex = snATAC,
    # nearest the union) to the right so the proteome-and-spatial edge is the far
    # vertical side, clear of the incoming funnels.
    tx_c = (1.40, 6.30)
    ge_c = (1.40, 2.95)
    uni_c = (5.50, 5.00)

    apex = (9.60, 5.00)                                # snATAC (nearest the union)
    s_side = 3.80                                      # triangle side length
    dx_v = s_side * np.sqrt(3.0) / 2.0                 # apex -> right-vertex x-offset
    P_c = (apex[0] + dx_v, apex[1] + s_side / 2.0)     # proteomics (top-right)
    S_c = (apex[0] + dx_v, apex[1] - s_side / 2.0)     # spatial (bottom-right)
    centroid = ((apex[0] + P_c[0] + S_c[0]) / 3.0,
                (apex[1] + P_c[1] + S_c[1]) / 3.0)
    HUB = {"P": (P_c, r_P), "S": (S_c, r_S), "A": (apex, r_A)}

    # ============================ canvas =============================
    fig, ax = plt.subplots(figsize=(5.0, 3.2))
    ax.set_aspect("equal")
    ax.axis("off")

    def node(center, r, c, z=8, ec="white"):
        ax.add_patch(Circle(center, r, facecolor=c, edgecolor=ec, lw=0.5, zorder=z))

    # ---- (z2) input funnels: transcriptomics + genetics -> union left face ----
    # Widths follow the node-area (sqrt) channel: each ribbon leaves its node at the
    # node's full face and tapers onto the union rim, allocated by substrate:coloc.
    uf_lo, uf_hi = uni_c[1] - r_uni, uni_c[1] + r_uni  # union left-face y-extent
    frac_tx = sub / float(max(sub + col, 1))
    bnd = uf_hi - (uf_hi - uf_lo) * frac_tx            # tx/genetics boundary on rim
    xL_uni = uni_c[0] - r_uni
    _ribbon(ax, tx_c[0], tx_c[1] - r_tx, tx_c[1] + r_tx,
            xL_uni, bnd, uf_hi, tint(C_BULK, 0.32), alpha=0.50)
    _ribbon(ax, ge_c[0], ge_c[1] - r_ge, ge_c[1] + r_ge,
            xL_uni, uf_lo, bnd, tint(C_GWAS, 0.32), alpha=0.48)

    # ---- (z2) validation funnels: union -> each hub, width propto validated count -
    # Mouths stack behind the union node (hidden by it, z8) and emerge at its rim;
    # each ribbon tapers to fit its hub face (min(width, ~1.5*hub_r)).
    fw = {"P": w_fun * Pv, "A": w_fun * Av, "S": w_fun * Sv}
    order = ["P", "A", "S"]                            # top -> bottom by target
    gap = 0.12
    tot = sum(fw.values()) + gap * (len(order) - 1)
    x0 = uni_c[0] + 0.30                               # inside union -> emerges at rim
    y = uni_c[1] + tot / 2.0
    for l in order:
        (hx, hy), hr = HUB[l]
        m_lo, m_hi = y - fw[l], y                      # mouth band (behind union)
        aw = min(fw[l], 1.5 * hr)                      # arrival width (taper to hub)
        _ribbon(ax, x0, m_lo, m_hi, hx, hy - aw / 2.0, hy + aw / 2.0,
                tint(LENS_COL[l], 0.32), alpha=0.50)
        y -= fw[l] + gap

    # ---- (z5-6) THE CLAIM: three shared-core edges with the ghost device ----
    def conv_edge(A, rA, B, rB, obs, exp):
        """Draw one shared-core edge. Above chance -> light expected core + saturated
        excess margins + solid material poking past a dashed gray ghost. At/below
        chance -> faded neutral band nested inside its dashed expected envelope."""
        A2, B2, perp = _seg(A, B, rA, rB)
        ow = w_con * max(obs, 0)                                    # observed width
        ew = w_con * exp if (np.isfinite(exp) and exp > 0) else ow  # expected width
        above = bool(np.isfinite(exp) and exp > 0 and obs > exp)
        if above:
            # light "expected volume" core
            ax.add_patch(Polygon(_band_poly(A2, B2, perp, ew / 2.0), closed=True,
                                 facecolor=tint(C_CONV, 0.52), edgecolor="none",
                                 zorder=5, alpha=0.95))
            # saturated above-chance excess on BOTH sides (ew/2 -> ow/2)
            for s in (1.0, -1.0):
                ax.add_patch(Polygon(_strip_poly(A2, B2, perp, s * ew / 2.0,
                                                 s * ow / 2.0), closed=True,
                                     facecolor=C_CONV, edgecolor="none",
                                     zorder=5.4, alpha=0.97))
            # dashed gray ghost outline at the expected width = the chance boundary
            ax.add_patch(Polygon(_band_poly(A2, B2, perp, ew / 2.0), closed=True,
                                 facecolor="none", edgecolor=GRAY, lw=0.8,
                                 ls=(0, (2.2, 1.6)), zorder=6))
        else:
            ax.add_patch(Polygon(_band_poly(A2, B2, perp, ow / 2.0), closed=True,
                                 facecolor=NEUTRAL, edgecolor="none",
                                 zorder=5, alpha=0.65))
            ax.add_patch(Polygon(_band_poly(A2, B2, perp, ew / 2.0), closed=True,
                                 facecolor="none", edgecolor=GRAY, lw=0.7,
                                 ls=(0, (2.0, 1.6)), zorder=6, alpha=0.9))
        mid = ((A2[0] + B2[0]) / 2.0, (A2[1] + B2[1]) / 2.0)
        return mid, perp, ow, above

    # proteome-and-spatial = the far vertical (right) edge -- the one claim edge
    ps_mid, ps_perp, ps_ow, ps_above = conv_edge(P_c, r_P, S_c, r_S, PS_o, PS_e)
    # proteome-and-snATAC + spatial-and-snATAC -- faded, nested (at/below chance)
    pa_mid, pa_perp, pa_ow, _ = conv_edge(P_c, r_P, apex, r_A, PA_o, PA_e)
    sa_mid, sa_perp, sa_ow, _ = conv_edge(S_c, r_S, apex, r_A, SA_o, SA_e)

    # ---- (z9) all-3 core at the centroid: 1 mark = 1 gene (tiny count -> unit dots)
    def unit_cluster(cx, cy, n, r=0.11):
        if n <= 0:
            return
        if n == 1:
            pts = [(cx, cy)]
        else:
            rr = 0.15
            pts = [(cx + rr * np.cos(t), cy + rr * np.sin(t))
                   for t in np.linspace(np.pi / 2, np.pi / 2 + 2 * np.pi, n,
                                        endpoint=False)]
        for px, py in pts:
            ax.add_patch(Circle((px, py), r, facecolor=C_CONV, edgecolor=INK,
                                lw=0.5, zorder=9))
    unit_cluster(*centroid, all3)
    ax.text(centroid[0], centroid[1] + 0.42, f"all 3: {all3}", ha="center",
            va="bottom", color=INK, fontsize=FS, zorder=10, bbox=_WBOX)

    # ---- (z8) nodes on top (cap the ribbon/edge junctions cleanly) ----
    node(tx_c, r_tx, C_BULK)
    node(ge_c, r_ge, C_GWAS)
    node(uni_c, r_uni, GRAY)
    node(P_c, r_P, C_PROTEO)
    node(S_c, r_S, C_SPATIAL)
    node(apex, r_A, C_ATAC)

    # shared-overlap seam at the input merge (neutral, height propto shared count)
    seam_h = shared / float(max(uni, 1)) * (2.0 * r_uni)
    ax.add_patch(Rectangle((xL_uni - 0.07, bnd - seam_h / 2.0), 0.14, seam_h,
                           facecolor=NEUTRAL, edgecolor="white", lw=0.3, zorder=8.4))
    ax.text(xL_uni - 0.16, bnd, f"shared\n{shared:,}", ha="right", va="center",
            color=INK, fontsize=FS, zorder=10)

    # ---- (z10) node labels: modality + live count (no gene names) ----
    ax.text(tx_c[0], tx_c[1] + r_tx + 0.14, f"transcriptomics\n{sub:,}", ha="center",
            va="bottom", color=INK, fontsize=FS, zorder=10, bbox=_WBOX)
    ax.text(ge_c[0], ge_c[1] - r_ge - 0.14, f"genetics\n{col:,}", ha="center",
            va="top", color=INK, fontsize=FS, zorder=10, bbox=_WBOX)
    ax.text(uni_c[0], uni_c[1], f"prioritized\n{uni:,}", ha="center", va="center",
            color=INK, fontsize=FS, zorder=10)
    ax.text(P_c[0], P_c[1] + r_P + 0.12, f"proteomics {Pv:,}", ha="center",
            va="bottom", color=INK, fontsize=FS, zorder=10, bbox=_WBOX)
    ax.text(S_c[0], S_c[1] - r_S - 0.12, f"spatial {Sv:,}", ha="center",
            va="top", color=INK, fontsize=FS, zorder=10, bbox=_WBOX)
    ax.text(apex[0], apex[1] + r_A + 0.12, f"snATAC {Av:,}", ha="center",
            va="bottom", color=INK, fontsize=FS, zorder=10, bbox=_WBOX)

    # ---- (z10) annotate ONLY the proteome-and-spatial edge (live fold + p) ----
    ann_x = ps_mid[0] + ps_perp[0] * (ps_ow / 2.0 + 0.18)
    ann_y = ps_mid[1] + ps_perp[1] * (ps_ow / 2.0 + 0.18)
    fstr = f"{PS_f:.2f}x" if np.isfinite(PS_f) else "n/a"
    estr = f"{PS_e:.1f}" if np.isfinite(PS_e) else "n/a"
    ax.text(ann_x, ann_y,
            f"proteome and spatial\nobs {PS_o} vs exp {estr}\n{fstr}  p={_pfmt(PS_p)}",
            ha="left", va="center", color=INK, fontsize=FS, zorder=10, bbox=_WBOX)

    # small count tags on the at/below-chance edges (backing their widths)
    for mid, perp, ow, cnt in ((pa_mid, pa_perp, pa_ow, PA_o),
                               (sa_mid, sa_perp, sa_ow, SA_o)):
        # offset the tag to the OUTER side (away from the centroid)
        s = 1.0 if ((mid[0] - centroid[0]) * perp[0]
                    + (mid[1] - centroid[1]) * perp[1]) >= 0 else -1.0
        ax.text(mid[0] + s * perp[0] * (ow / 2.0 + 0.14),
                mid[1] + s * perp[1] * (ow / 2.0 + 0.14),
                f"{cnt}", ha="center", va="center", color=INK, fontsize=FS,
                zorder=10, bbox=_WBOX)

    # ---- (z10) size legend: node area = gene count (two reference counts) ----
    refs = [5000, 500]
    base_y = -1.30
    rL = rad(refs[0])
    rS = rad(refs[1])
    xb = 4.00                                          # bottom-center, clear of nodes
    xs = xb + rL + rS + 0.45
    ax.text(xb - rL, base_y + 2 * rL + 0.30, "node area = gene count",
            ha="left", va="bottom", color=INK, fontsize=FS, zorder=10)
    for cx, rr, n in ((xb, rL, refs[0]), (xs, rS, refs[1])):
        ax.add_patch(Circle((cx, base_y + rr), rr, facecolor="none", edgecolor=INK,
                            lw=0.6, zorder=10))
        ax.text(cx, base_y - 0.14, f"{n:,}", ha="center", va="top", color=INK,
                fontsize=FS, zorder=10)

    # ---- (z10) ghost-device key (bottom-right) ----
    gx, gy = 8.30, -0.42
    ax.add_patch(Rectangle((gx, gy - 0.13), 0.52, 0.26, facecolor=C_CONV,
                           edgecolor="none", zorder=10))
    ax.text(gx + 0.64, gy, "observed (above chance)", ha="left", va="center",
            color=INK, fontsize=FS, zorder=10)
    ax.add_patch(Rectangle((gx, gy - 0.60 - 0.13), 0.52, 0.26, facecolor="none",
                           edgecolor=GRAY, lw=0.8, ls=(0, (2.2, 1.6)), zorder=10))
    ax.text(gx + 0.64, gy - 0.60, "expected (chance)", ha="left", va="center",
            color=INK, fontsize=FS, zorder=10)

    # ---- (z10) stage tags ----
    # header over the funnel corridor (x~4), clear of the transcriptomics node label
    ax.text(4.00, 8.05, "prioritization", ha="center", va="bottom", color=INK,
            fontsize=FS, zorder=10)
    ax.text(centroid[0], 8.05, "orthogonal validation", ha="center", va="bottom",
            color=INK, fontsize=FS, zorder=10)

    # ---- fixed limits (aspect equal; bbox_inches='tight' crops on save) ----
    ax.set_xlim(-0.6, 15.4)
    ax.set_ylim(-1.95, 8.35)

    save_pdf(fig, "fig4a_triangle.pdf")

    # ---- printed caption (full claim + definitions; NOT an on-panel subtitle) ----
    print("[caption:triangle] The three orthogonal validation assays share a "
          "convergent core, but only the proteome-and-spatial link exceeds chance. "
          "LEFT (prioritization): two evidence streams -- transcriptomics {sub:,} and "
          "genetics {col:,} ({sh:,} shared, drawn as the neutral merge seam) -- funnel "
          "into one prioritized universe of {uni:,} genes; node area is proportional to "
          "gene count (radius propto sqrt(count); see the size legend). RIGHT "
          "(validation triangle): the universe feeds three orthogonal lenses -- "
          "proteomics {Pv:,}, spatial {Sv:,} and snATAC {Av:,} validated genes -- with "
          "snATAC at the apex nearest the union so the proteome-and-spatial edge is the "
          "far side; funnel-edge width is proportional to each hub's validated count. "
          "The three internal edges are the EXCLUSIVE pairwise shared-gene cores (edge "
          "width proportional to the shared count): proteome-and-spatial {PSo} is drawn "
          "solid with a light expected core and a saturated excess margin, overshooting "
          "its gray dashed permutation-EXPECTED width ({PSe} genes); proteome-and-snATAC "
          "{PAo} and spatial-and-snATAC {SAo} are faded and nest INSIDE their expected "
          "envelopes (at or below chance). The central marks are the {a3} genes shared by "
          "all three assays (one mark = one gene). Only proteome-and-spatial beats chance "
          "({PSo} observed vs {PSe} expected, {PSf}x, p={PSp}; seed-42 10,000-draw "
          "label-permutation null over the co-measured universe) -- above chance is "
          "encoded by edge width + solidity + the ghost-overshoot so it survives "
          "grayscale. Funnel-flow widths and shared-core-edge widths are drawn on two "
          "internally-proportional linear scales (their count magnitudes differ >10x) and "
          "every edge is labelled with its exact count. snATAC = single-nucleus ATAC; "
          "every node area and edge width encodes a live count; no gene names.".format(
              sub=sub, col=col, sh=shared, uni=uni, Pv=Pv, Sv=Sv, Av=Av,
              PSo=PS_o, PAo=PA_o, SAo=SA_o, a3=all3,
              PSe=(f"{PS_e:.1f}" if np.isfinite(PS_e) else "n/a"),
              PSf=(f"{PS_f:.2f}" if np.isfinite(PS_f) else "n/a"),
              PSp=_pfmt(PS_p)))


if __name__ == "__main__":
    d, stats = get_data()
    build(d, stats)
