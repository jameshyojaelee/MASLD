#!/usr/bin/env python3
"""
Fig 4A v2 -- CANDIDATE 1: CHANCE-ANCHORED ENRICHMENT FOREST (+ count strip).

The claim this geometry welds on (vs the retired display-only Sankey):
  prioritized targets reproduce independent disease signal ABOVE CHANCE;
  proteome-and-spatial is the one robust convergent core, while overall >=2 and
  genetics-alone sit at/below chance.

Two stacked parts (shared data layer = fig4a_v2_common.get_data()):

  TOP STRIP (opener) -- introduces the layers + modalities + counts:
    * one slim 100%-stacked bar of the prioritized universe, split by ORIGIN into
      transcriptomic-only (C_BULK) + shared (neutral gray) + genetic-only (C_GWAS),
      widths proportional to the EXCLUSIVE partition d["_part"], with ASCII counts.
    * directly below, three short modality-validated segments -- proteomics P_val
      (C_PROTEO) / spatial S_val (C_SPATIAL) / snATAC A_val (C_ATAC) -- drawn on the
      SAME universe scale (so they read as the small subsets they are), with counts.

  BOTTOM FOREST (the claim) -- fold-enrichment on a LOG2 axis, chance (fold=1) as a
  grey dashed reference. Points + a 95% NULL BAND whisker (NOT lollipops; no stem to
  the baseline). Three groups, hairline-separated, short left-side group labels:
    (a) each lens vs a genome-wide hypergeometric background (stats["enrich"]
        ["universe"]); point coloured by lens; whisker = normal-approx 95% null band
        from that entry's own N/K/k (fold_ci on conv[] would be the WRONG null here).
    (b) EXCLUSIVE intersections vs the 10k-permutation null (stats["conv"]): proteome
        and spatial (the one core, sits right of chance) + any >=2 of 3 (straddles it);
        whisker = fold_ci(...) 2.5-97.5% band.
    (c) genetics-alone contrast (stats["enrich"]["coloc_only"]): each lens on the
        genetic-only set -- points at/left of chance (enriches in nothing), kept ON the
        panel so the one real signal is credible.
  Filled marker = beats chance (p<0.05, fold>1); open marker = at/below chance.

House style (fig4a_v2_common sets rcParams -- do NOT re-set): PDF, Helvetica 6, all
text INK, no bold, NO gene names, colour on marks only, ASCII-only on-panel, folds
guarded against NaN/inf so the panel always renders. Full claim in the printed caption.
"""
from fig4a_v2_common import (get_data, save_pdf, INK, FS, GRAY, C_NS, C_CONV,
    C_BULK, C_GWAS, C_PROTEO, C_SPATIAL, C_ATAC, LENS_COL, LENS_NAME,
    tint, fold_ci, _ribbon, plt, np, Rectangle, PathPatch, FancyBboxPatch,
    Circle, Wedge, Ellipse, FancyArrowPatch, Polygon, Path)

# Contract-sanctioned neutral gray for the shared-origin segment (task spec).
NEUTRAL = "#B9BEC4"


# -- small helpers ------------------------------------------------------------
def _lg(v):
    """log2, guarding NaN / inf / non-positive so the panel never raises."""
    try:
        return float(np.log2(v)) if (v is not None and np.isfinite(v) and v > 0) else float("nan")
    except (TypeError, ValueError):
        return float("nan")


def _pfmt(p):
    """Compact ASCII p-string (Helvetica lacks fancy glyphs): '2e-3', '0.30', 'n/a'."""
    if p is None or not np.isfinite(p):
        return "n/a"
    if p >= 0.01:
        return f"{p:.2f}"
    m, e = f"{p:.0e}".split("e")
    return f"{m}e{int(e)}"


def _perlens_row(label, e, color):
    """Enrichment row from a stats['enrich'][...][lens] entry {obs,exp,fold,p,N,K,k}.
    The whisker is a NORMAL approximation to the hypergeometric 95% null band built
    from THIS entry's own N/K/k (its native null), NOT fold_ci on a conv[] entry --
    the per-lens hypergeometric and the permutation intersection null are different
    nulls, so borrowing the latter's band would misstate this row."""
    e = e or {}
    fold = e.get("fold", float("nan"))
    exp = e.get("exp", float("nan"))
    N, K, k = e.get("N", 0), e.get("K", 0), e.get("k", 0)
    lo_fold = hi_fold = float("nan")
    if np.isfinite(exp) and exp > 0 and N > 1 and k > 0 and 0 < K < N:
        var = k * (K / N) * (1.0 - K / N) * ((N - k) / (N - 1))
        sd = float(np.sqrt(max(var, 0.0)))
        lo_fold = max(exp - 1.96 * sd, 0.0) / exp
        hi_fold = (exp + 1.96 * sd) / exp
    p = e.get("p", float("nan"))
    sig = bool(np.isfinite(fold) and np.isfinite(p) and p < 0.05 and fold > 1.0)
    return dict(label=label, fold=fold, lo=lo_fold, hi=hi_fold,
                obs=e.get("obs", 0), exp=exp, p=p, color=color, sig=sig)


def _conv_row(label, st, color):
    """Intersection row from a stats['conv'][key] entry {obs,exp,fold,p,lo,hi}.
    Whisker = fold_ci(...) -> the 2.5-97.5% permutation null band expressed in fold."""
    st = st or {}
    fold, lo_fold, hi_fold = fold_ci(st)
    p = st.get("p", float("nan"))
    sig = bool(np.isfinite(fold) and np.isfinite(p) and p < 0.05 and fold > 1.0)
    return dict(label=label, fold=fold, lo=lo_fold, hi=hi_fold,
                obs=st.get("obs", 0), exp=st.get("exp", float("nan")), p=p,
                color=color, sig=sig)


# -- top strip ----------------------------------------------------------------
def _draw_strip(ax, d):
    """100%-stacked universe origin bar + three modality-validated segments."""
    ax.axis("off")
    U = max(int(d["universe"]), 1)
    part = d["_part"]
    # full-width strip so the two narrow origin segments (shared, genetic) have
    # room for centred name+count labels without colliding.
    ax.set_xlim(-0.28 * U, 1.26 * U)
    ax.set_ylim(-0.7, 4.8)

    # --- prioritized universe, split by ORIGIN (exclusive partition) ---
    segs = [("transcriptomic", int(part["tx_only"]),      C_BULK),
            ("shared",         int(part["overlap"]),      NEUTRAL),
            ("genetic",        int(part["genetic_only"]), C_GWAS)]
    yb, hb, x = 3.35, 0.82, 0.0
    for name, w, col in segs:
        ax.add_patch(Rectangle((x, yb - hb / 2), w, hb, facecolor=col,
                               edgecolor="white", lw=0.5, zorder=3))
        xm = x + w / 2.0
        ax.text(xm, yb + hb / 2 + 0.16, name, ha="center", va="bottom",
                color=INK, fontsize=FS, zorder=5)
        ax.text(xm, yb - hb / 2 - 0.16, f"{w:,}", ha="center", va="top",
                color=INK, fontsize=FS, zorder=5)
        x += w
    ax.text(-0.03 * U, yb, "prioritized", ha="right", va="center",
            color=INK, fontsize=FS, zorder=5)
    ax.text(U * 1.01, yb, f"n={U:,}", ha="left", va="center",
            color=INK, fontsize=FS, zorder=5)

    # --- three modality-validated segments on the SAME universe scale ---
    lens = [("proteomics", int(d["P_val"]), C_PROTEO),
            ("spatial",    int(d["S_val"]), C_SPATIAL),
            ("snATAC",      int(d["A_val"]), C_ATAC)]
    ys, hl = (1.55, 0.80, 0.05), 0.52
    for (name, cnt, col), yy in zip(lens, ys):
        ax.add_patch(Rectangle((0, yy - hl / 2), max(cnt, 0), hl, facecolor=col,
                               edgecolor="none", zorder=3))
        ax.text(cnt + 0.012 * U, yy, f"{name} {cnt:,}", ha="left", va="center",
                color=INK, fontsize=FS, zorder=5)
    ax.text(-0.03 * U, ys[1], "validated", ha="right", va="center",
            color=INK, fontsize=FS, zorder=5)


# -- main builder -------------------------------------------------------------
def build(d, stats):
    if not stats or "conv" not in stats or "enrich" not in stats:
        print("[build:forest] stats missing conv/enrich -- skipping render"); return
    conv = stats["conv"]
    enrU = stats["enrich"].get("universe", {}) or {}
    enrC = stats["enrich"].get("coloc_only", {}) or {}

    # ---- rows, grouped ----
    groups = [
        ("lens vs genome", [
            _perlens_row("proteomics", enrU.get("proteomics"), C_PROTEO),
            _perlens_row("spatial",    enrU.get("spatial"),    C_SPATIAL),
            _perlens_row("snATAC",     enrU.get("scATAC"),     C_ATAC)]),
        ("convergence", [
            _conv_row("proteome and spatial", conv.get("PS"),  C_CONV),
            _conv_row(">=2 of 3 lenses",      conv.get("ge2"), C_CONV)]),
        ("genetics alone", [
            _perlens_row("proteomics", enrC.get("proteomics"), C_GWAS),
            _perlens_row("spatial",    enrC.get("spatial"),    C_GWAS),
            _perlens_row("snATAC",     enrC.get("scATAC"),     C_GWAS)]),
    ]

    # ---- x extent from finite fold / band ends (guarded) ----
    xs = []
    for _, rows in groups:
        for r in rows:
            for v in (r["fold"], r["lo"], r["hi"]):
                lv = _lg(v)
                if np.isfinite(lv):
                    xs.append(lv)
    if not xs:
        xs = [0.0]
    xmin, xmax = min(xs), max(xs)
    b_lo = min(-1.0, xmin - 0.05)      # spine/tick span (fold 0.5 .. 2 at least)
    b_hi = max(1.0, xmax + 0.05)
    X_LAB = b_lo - 0.55                 # left label anchor (ha=right, overflows out)
    XANN = b_hi + 0.55                  # right annotation anchor (ha=left, overflows out)

    # ---- y layout (top -> bottom): header, rows, gap+separator ----
    header_gap, row_gap, group_gap = 0.75, 1.0, 0.8
    rows_xy, headers_xy, seps = [], [], []
    y, first = 0.0, True
    for htext, rows in groups:
        if not first:
            seps.append(y - group_gap * 0.5)
            y -= group_gap
        first = False
        headers_xy.append((y, htext))
        y -= header_gap
        for r in rows:
            rows_xy.append((y, r))
            y -= row_gap
    y_bottom = y + row_gap             # last actually-used row

    fig = plt.figure(figsize=(3.7, 3.0))
    # strip = full width (segments fit); forest = indented so its left row labels /
    # group labels and right annotations sit in the reserved margins.
    ax_top = fig.add_axes([0.05, 0.72, 0.90, 0.22])
    ax = fig.add_axes([0.30, 0.055, 0.52, 0.60])
    _draw_strip(ax_top, d)

    # ---- forest marks ----
    for yy, r in rows_xy:
        lo, hi, fold = r["lo"], r["hi"], r["fold"]
        # 95% null band whisker (grey), guarded
        xl, xh = _lg(lo), _lg(hi)
        if np.isfinite(xl) and np.isfinite(xh):
            ax.plot([xl, xh], [yy, yy], color=GRAY, lw=0.8, alpha=0.75,
                    solid_capstyle="round", zorder=3)
            for xe in (xl, xh):
                ax.plot([xe, xe], [yy - 0.12, yy + 0.12], color=GRAY, lw=0.8,
                        alpha=0.75, zorder=3)
        # observed-fold point: filled = beats chance, open = at/below chance
        xp = _lg(fold)
        if np.isfinite(xp):
            if r["sig"]:
                ax.scatter([xp], [yy], s=24, facecolor=r["color"],
                           edgecolor="white", linewidths=0.4, zorder=6)
            else:
                ax.scatter([xp], [yy], s=20, facecolor="white",
                           edgecolor=r["color"], linewidths=0.9, zorder=6)
        # left row label
        ax.text(X_LAB, yy, r["label"], ha="right", va="center",
                color=INK, fontsize=FS, zorder=7)
        # right annotation: obs / exp  p=...
        if np.isfinite(r["exp"]):
            ann = f"{int(r['obs']):,} / {r['exp']:.1f}  p={_pfmt(r['p'])}"
        else:
            ann = f"{int(r['obs']):,}  null n/a"
        ax.text(XANN, yy, ann, ha="left", va="center",
                color=INK, fontsize=FS, zorder=7)

    # group labels (left side) + hairline separators
    for hy, htext in headers_xy:
        ax.text(X_LAB, hy, htext, ha="right", va="center",
                color=INK, fontsize=FS, zorder=7)
    for sy in seps:
        ax.plot([b_lo, b_hi], [sy, sy], color=C_NS, lw=0.5, zorder=1)

    # chance reference (fold = 1)
    ax.axvline(0.0, color=GRAY, lw=0.8, ls=(0, (3, 2)), zorder=2)
    ax.text(0.0, 0.55, "chance", ha="center", va="bottom",
            color=INK, fontsize=FS, zorder=7)

    # ---- axis cosmetics: bottom spine only, log2 ticks at fold 0.5 / 1 / 2 ----
    ax.set_xlim(b_lo - 0.45, b_hi + 0.45)
    ax.set_ylim(y_bottom - 0.55, 1.15)
    for s in ("top", "left", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(INK)
    ax.spines["bottom"].set_linewidth(0.5)
    ax.spines["bottom"].set_bounds(b_lo, b_hi)
    ax.set_yticks([])
    ax.set_xticks([_lg(0.5), 0.0, _lg(2.0)])
    ax.set_xticklabels(["0.5", "1", "2"])
    ax.tick_params(axis="x", length=2.5, width=0.5, labelsize=FS, colors=INK)
    ax.set_xlabel("fold enrichment (observed / expected)", fontsize=FS, color=INK)

    save_pdf(fig, "fig4a_forest.pdf")

    # ---- printed caption (full claim + definitions; NOT an on-panel subtitle) ----
    ps, ge2 = conv.get("PS", {}), conv.get("ge2", {})
    part = d["_part"]
    print("[caption:forest] Prioritized targets reproduce independent disease signal "
          "above chance. TOP: the prioritized universe (n={U:,}) splits by origin into "
          "transcriptomic-only {tx:,}, shared {sh:,} and genetic-only {go:,}; three "
          "orthogonal validation lenses recover proteomics {P:,}, spatial {S:,} and "
          "snATAC {A:,} genes (widths on the universe scale). FOREST: fold-enrichment "
          "(observed/expected) on a log2 axis; chance = 1 (grey dashed); filled marker "
          "beats chance (p<0.05, fold>1), open does not. 'lens vs genome' tests each "
          "lens's disease signal in the universe against a genome-wide hypergeometric "
          "background (grey whisker = normal-approx 95% null band from that lens's "
          "N/K/k; obs/exp and empirical p at right). 'convergence' tests the EXCLUSIVE "
          "proteome-and-spatial and any->=2-of-3 intersections against a 10,000-draw "
          "label-permutation null (seed 42; whisker = 2.5-97.5% null band). "
          "proteome-and-spatial is the one robust convergent core "
          "({ps_o} vs {ps_e:.1f}, {ps_f:.2f}x, p={ps_p}), sitting right of chance; "
          "overall >=2 ({g_o} vs {g_e:.1f}, {g_f:.2f}x, p={g_p}) straddles chance. "
          "'genetics alone' tests the genetic-only set against each lens: all points sit "
          "at/left of chance (enriches in nothing), so convergence is not an artefact of "
          "the genetic axis. snATAC = single-nucleus ATAC; every mark encodes a live "
          "count; no gene names.".format(
              U=int(d["universe"]), tx=int(part["tx_only"]), sh=int(part["overlap"]),
              go=int(part["genetic_only"]), P=int(d["P_val"]), S=int(d["S_val"]),
              A=int(d["A_val"]),
              ps_o=int(ps.get("obs", 0)), ps_e=float(ps.get("exp", float("nan"))),
              ps_f=float(ps.get("fold", float("nan"))), ps_p=_pfmt(ps.get("p")),
              g_o=int(ge2.get("obs", 0)), g_e=float(ge2.get("exp", float("nan"))),
              g_f=float(ge2.get("fold", float("nan"))), g_p=_pfmt(ge2.get("p"))))


if __name__ == "__main__":
    d, stats = get_data()
    build(d, stats)
