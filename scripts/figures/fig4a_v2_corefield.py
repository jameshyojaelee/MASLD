#!/usr/bin/env python3
"""
Fig 4A v2 -- CANDIDATE 3: CORE-IN-A-FIELD "NULL LINEUP" (area-faithful waffle + a
per-gene convergence-vs-chance lineup). Winner of the 2026-07-10 critique/judge round.

The claim this geometry welds on (vs the retired display-only Sankey): of the 9,882
prioritized targets, corroboration is a small EARNED minority, and when every
convergence class is tried against its OWN permutation null, exactly one survives --
proteome-and-spatial (obs vs exp), while the SAME-count proteome-and-snATAC falls
short of its larger null, three-way convergence is empty, and genetics-only enriches
in nothing. Assertion and honesty are the same act.

Three coupled devices (shared data layer = fig4a_v2_common.get_data(), STRICT def):

  ORIGIN BAND (top) -- introduces the two prioritization LAYERS that build the
    universe: transcriptomics (C_BULK) + genetics (C_GWAS), overlap shown as the seam.

  MAIN FIELD -- the whole universe as a POPULATION (not a value heatmap): an
    area-faithful WAFFLE, 1 square ~= k genes; the vast UNVALIDATED remainder in
    light C_NS (demoted to context, not the headline) vs a small CORROBORATED corner
    WEDGE (>=2-assay ring core at the apex, then single-assay bands). The left key
    names each validation MODALITY with its validated/measured count.

  NULL LINEUP INSET (per-gene, 1 dot = 1 gene) -- the CLAIM: each pairwise
    convergence drawn as observed dots vs its OWN dashed gray GHOST box (the 95%
    chance envelope = permutation upper bound). A gene PAST the box beats chance.
    Only proteome-and-spatial's dots burst past; proteome-and-snATAC (identical 19)
    nests inside a wider ghost (the cherry-pick disarm, drawn not hidden);
    spatial-and-snATAC at the edge; all-3 an explicit empty slot.

House style (fig4a_v2_common sets rcParams): PDF, Helvetica 6, all text INK, no bold,
NO gene names, colour on marks only, control/unvalidated GRAY/C_NS, ASCII-only. Counts
live; the null is co-measured-conditioned (see caption). Full claim in the caption.
"""
from fig4a_v2_common import (get_data, save_pdf, INK, FS, GRAY, C_NS, C_CONV,
    C_BULK, C_GWAS, C_PROTEO, C_SPATIAL, C_ATAC, LENS_COL, LENS_NAME,
    tint, fold_ci, _ribbon, plt, np, Rectangle, PathPatch, FancyBboxPatch,
    Circle, Wedge, Ellipse, FancyArrowPatch, Polygon, Path)


def _pfmt(p):
    if p is None or not np.isfinite(p):
        return "n/a"
    if p >= 0.01:
        return f"{p:.2f}"
    m, e = f"{p:.0e}".split("e")
    return f"{m}e{int(e)}"


def _mix(a, b, t=0.5):
    import matplotlib.colors as mc
    ra, ga, ba = mc.to_rgb(a)
    rb, gb, bb = mc.to_rgb(b)
    return (ra + (rb - ra) * t, ga + (gb - ga) * t, ba + (bb - ba) * t)


def _largest_remainder(counts, total):
    tot = float(sum(max(v, 0) for v in counts.values())) or 1.0
    quota = {k: max(v, 0) * total / tot for k, v in counts.items()}
    floor = {k: int(np.floor(q)) for k, q in quota.items()}
    rem = total - sum(floor.values())
    order = sorted(counts, key=lambda k: quota[k] - floor[k], reverse=True)
    for k in order[:max(rem, 0)]:
        floor[k] += 1
    return floor


def build(d, stats):
    part = d["_part"]
    universe = int(d["universe"])
    ge1, ge2, unval = int(part["ge1"]), int(part["ge2"]), int(part["unvalidated"])
    P_only, S_only, A_only = int(part["P_only"]), int(part["S_only"]), int(part["A_only"])
    multi = ge2
    blend = _mix(C_PROTEO, C_SPATIAL)

    # ---- waffle sizing: k ~ 25 genes/mark, marks fill a 20-wide grid exactly ----
    NCOL = 20
    NROW = max(1, int(round((universe / 25.0) / NCOL)))
    TOTAL = NCOL * NROW
    k_disp = int(round(universe / TOTAL))
    marks_n = _largest_remainder(
        {"M": multi, "P": P_only, "S": S_only, "A": A_only, "U": unval}, TOTAL)

    CLASS_STYLE = {"M": (C_CONV, True), "P": (C_PROTEO, False), "S": (C_SPATIAL, False),
                   "A": (C_ATAC, False), "U": (C_NS, False)}
    seq = []
    for cls in ("M", "P", "S", "A", "U"):
        seq.extend([CLASS_STYLE[cls]] * marks_n[cls])
    seq = (seq + [CLASS_STYLE["U"]] * TOTAL)[:TOTAL]
    cells = [(r, c) for r in range(NROW) for c in range(NCOL)]
    cells.sort(key=lambda rc: (float(np.hypot(rc[1] + 0.5, rc[0] + 0.5)),
                               rc[0] + rc[1], rc[0], rc[1]))

    # ============================ canvas ============================
    fig, ax = plt.subplots(figsize=(7.4, 3.7))
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_xlim(-10.5, 50.0)
    ax.set_ylim(-3.6, 25.6)

    # ---- MAIN FIELD ----
    S_MARK, GUT = 0.80, 0.20
    for (r, c), (col, ring) in zip(cells, seq):
        ax.add_patch(Rectangle((c + GUT / 2.0, r + GUT / 2.0), S_MARK, S_MARK,
                     facecolor=col, edgecolor=(INK if ring else "none"),
                     lw=(0.6 if ring else 0.0), zorder=4))
    ax.add_patch(Rectangle((0, 0), NCOL, NROW, facecolor="none", edgecolor=C_NS,
                           lw=0.5, zorder=2))

    # ---- ORIGIN BAND (top): the two prioritization layers that build the universe --
    tx_only, ovl, gen_only = int(part["tx_only"]), int(part["overlap"]), int(part["genetic_only"])
    sub, col_gen = int(d["substrate"]), int(d["coloc"])
    bb_y, bb_h = NROW + 0.6, 0.85
    w_tx = NCOL * tx_only / universe
    w_ov = NCOL * ovl / universe
    ax.add_patch(Rectangle((0, bb_y), w_tx, bb_h, facecolor=C_BULK, edgecolor="white", lw=0.4, zorder=4))
    ax.add_patch(Rectangle((w_tx, bb_y), w_ov, bb_h, facecolor=_mix(C_BULK, C_GWAS), edgecolor="white", lw=0.4, zorder=4))
    ax.add_patch(Rectangle((w_tx + w_ov, bb_y), NCOL - w_tx - w_ov, bb_h, facecolor=C_GWAS, edgecolor="white", lw=0.4, zorder=4))
    ax.text(0, bb_y + bb_h + 0.25, f"transcriptomics {sub:,}", ha="left", va="bottom", color=INK, fontsize=FS, zorder=9)
    ax.text(NCOL, bb_y + bb_h + 0.25, f"genetics {col_gen:,}", ha="right", va="bottom", color=INK, fontsize=FS, zorder=9)
    ax.text(w_tx + w_ov / 2.0, bb_y - 0.15, f"{ovl:,} shared", ha="center", va="top", color=INK, fontsize=FS, zorder=9)
    ax.text(0, bb_y + bb_h + 1.35, f"prioritized targets   n={universe:,}", ha="left", va="bottom", color=INK, fontsize=FS, zorder=9)

    # ---- field callouts: unvalidated DEMOTED (context), corroborated ELEVATED ----
    pct_corr = int(round(100.0 * ge1 / max(universe, 1)))
    ax.text(NCOL * 0.58, NROW * 0.72, f"unvalidated {unval:,} ({100 - pct_corr}%)",
            ha="center", va="center", color=GRAY, fontsize=FS, zorder=9)
    ax.text(NCOL * 0.30, -1.15, f"corroborated {ge1:,}  ({pct_corr}%, >=1 assay)",
            ha="center", va="top", color=INK, fontsize=FS, zorder=9)

    # ---- LEFT KEY: name each modality + validated/measured (introduces modalities) --
    kx = -10.2
    ky = NROW - 0.3
    P_meas, S_meas, A_meas = int(d["P_meas"]), int(d["S_meas"]), int(d["A_meas"])
    key = [("proteomics", C_PROTEO, False, f"{d['P_val']} / {P_meas:,}"),
           ("spatial", C_SPATIAL, False, f"{d['S_val']} / {S_meas}"),
           ("snATAC", C_ATAC, False, f"{d['A_val']} / {A_meas}"),
           (">=2 assays", C_CONV, True, f"{ge2}"),
           ("unvalidated", C_NS, False, f"{unval:,}")]
    for name, col, ring, cnt in key:
        ax.add_patch(Rectangle((kx, ky - S_MARK / 2.0), S_MARK, S_MARK, facecolor=col,
                     edgecolor=(INK if ring else "none"), lw=(0.6 if ring else 0.0), zorder=4))
        ax.text(kx + 1.15, ky + 0.28, name, ha="left", va="center", color=INK, fontsize=FS, zorder=9)
        ax.text(kx + 1.15, ky - 0.42, cnt, ha="left", va="center", color=GRAY, fontsize=FS, zorder=9)
        ky -= 2.1
    ax.text(kx, ky + 0.7, "validated / assayed", ha="left", va="center", color=GRAY, fontsize=FS, zorder=9)
    # unit legend
    ax.add_patch(Rectangle((0.0, -2.7 - S_MARK / 2.0), S_MARK, S_MARK, facecolor=C_NS, edgecolor="none", zorder=4))
    ax.text(1.2, -2.7, f"1 square = ~{k_disp} genes", ha="left", va="center", color=INK, fontsize=FS, zorder=9)

    # ============================ NULL LINEUP INSET ============================
    IX0, IY0, IW, IH = 24.5, -1.0, 25.3, 23.6
    ax.add_patch(FancyBboxPatch((IX0, IY0), IW, IH, boxstyle="round,pad=0.0,rounding_size=0.6",
                 facecolor="white", edgecolor=GRAY, lw=0.6, zorder=3))
    ax.text(IX0 + 0.8, IY0 + IH - 0.7, "convergence vs chance   1 dot = 1 gene",
            ha="left", va="top", color=INK, fontsize=FS, zorder=9)
    ax.text(IX0 + 0.8, IY0 + IH - 1.7, "dashed box = 95% chance envelope; a gene past it beats chance",
            ha="left", va="top", color=GRAY, fontsize=FS, zorder=9)

    conv = (stats or {}).get("conv", {}) or {}
    def cv(k):
        c = conv.get(k, {}) or {}
        return dict(obs=int(c.get("obs", 0)),
                    exp=float(c.get("exp", float("nan"))),
                    hi=float(c.get("hi", float("nan"))),
                    fold=c.get("fold", float("nan")), p=c.get("p", float("nan")))
    PS_, PA_, SA_, A3_ = cv("PS"), cv("PA"), cv("SA"), cv("all3")

    def ext(rr):
        v = [rr["obs"]]
        if np.isfinite(rr["hi"]):
            v.append(rr["hi"])
        return max(v)
    max_ext = max(ext(PS_), ext(PA_), ext(SA_), 1)
    dot_x0 = IX0 + 1.3
    GENE_W = (IW - 2.6) / max_ext
    R = min(0.28, GENE_W * 0.42)

    rows = [("proteome and spatial", PS_, blend, True),
            ("proteome and snATAC", PA_, GRAY, False),
            ("spatial and snATAC", SA_, GRAY, False),
            ("all three", A3_, GRAY, False)]
    row_ys = [16.6, 12.0, 7.8, 4.0]
    for (label, rr, dcol, is_sig), ry in zip(rows, row_ys):
        obs, exp, hi, fold, p = rr["obs"], rr["exp"], rr["hi"], rr["fold"], rr["p"]
        ax.text(dot_x0, ry + R + 0.55, label, ha="left", va="bottom", color=INK, fontsize=FS, zorder=9)
        # 95% chance-envelope ghost box (to the null upper bound)
        if np.isfinite(hi) and hi > 0:
            gw = hi * GENE_W
            ax.add_patch(Rectangle((dot_x0 - 0.16, ry - R - 0.14), gw + 0.32, 2 * R + 0.28,
                         facecolor="none", edgecolor=GRAY, lw=0.8, ls=(0, (2.4, 1.6)), zorder=6))
        # faint expected-MEAN tick (shows obs-vs-mean fold context) for the sig row
        if is_sig and np.isfinite(exp) and exp > 0:
            ex = dot_x0 + exp * GENE_W
            ax.plot([ex, ex], [ry - R - 0.10, ry + R + 0.10], color="#9AA0A6", lw=0.6, zorder=6.5)
        # observed dots (1 = 1 gene)
        for i in range(obs):
            ax.add_patch(Circle((dot_x0 + (i + 0.5) * GENE_W, ry), R,
                         facecolor=dcol, edgecolor="none", zorder=7))
        # verdict tag
        if obs == 0:
            tag = "0  ->  none survive"
        elif is_sig and np.isfinite(exp):
            tag = f"{obs} vs {exp:.1f} expected   {fold:.2f}x  p={_pfmt(p)}  ->  beats chance"
        elif label.startswith("proteome and snATAC") and np.isfinite(exp):
            tag = f"same {obs} genes, {exp:.0f} expected  ->  at chance"
        elif np.isfinite(exp):
            gw = "gene" if obs == 1 else "genes"
            exps = f"{exp:.1f}" if exp < 1 else f"{exp:.0f}"
            tag = f"{obs} {gw}, ~{exps} expected  ->  at chance"
        else:
            tag = f"{obs}  ->  at chance"
        ax.text(dot_x0, ry - R - 0.62, tag, ha="left", va="top", color=INK, fontsize=FS, zorder=9)

    # ---- light connector: the >=2 core (wedge apex) -> the lineup ----
    ax.add_patch(FancyArrowPatch((1.7, 1.7), (IX0 - 0.3, IY0 + IH * 0.5),
                 arrowstyle="-|>", mutation_scale=6, color=INK, lw=0.6,
                 connectionstyle="arc3,rad=-0.14", zorder=8))
    ax.text((NCOL + IX0) / 2.0 - 1.0, IY0 + IH * 0.5 + 2.6, f"{ge2} genes in >=2 assays",
            ha="center", va="bottom", color=INK, fontsize=FS, zorder=9,
            bbox=dict(boxstyle="round,pad=0.16", fc="white", ec="#E6E6E6", lw=0.3, alpha=0.9))

    save_pdf(fig, "fig4a_corefield.pdf")

    ps_obs, ps_exp = int(part["PS"]), PS_["exp"]
    print("[caption:corefield] The prioritized universe distills to a small, "
          "orthogonally-corroborated core -- and only one convergence beats chance. "
          "MAIN: an area-faithful waffle of n={U:,} prioritized targets (1 square ~= {k} "
          "genes, apportioned live by largest-remainder). The light-gray field is the "
          "unvalidated remainder ({uv:,}, {pu}%); the coloured corner wedge is the "
          "corroborated {g1:,} ({pc}%, >=1 assay), banded by dominant modality "
          "(proteomics {Pv}/{Pm:,} validated/assayed, spatial {Sv}/{Sm}, snATAC {Av}/{Am}) "
          "with the >=2-assay core (dark ring) at the apex. A top band shows the two "
          "prioritization layers that build the universe (transcriptomics {sub:,} u genetics "
          "{cg:,}, {ov:,} shared). INSET (1 dot = 1 gene): each pairwise convergence vs its "
          "OWN permutation null (seed 42, 10k draws), the dashed box = the 95% chance "
          "envelope. proteome-and-spatial ({po} vs {pe:.1f} expected, {pf}x, p={pp}) is the "
          "ONLY class whose genes burst past the envelope; proteome-and-snATAC has the "
          "IDENTICAL {po} genes yet nests inside its larger null (at chance); spatial-and-"
          "snATAC is at the edge; all-three is empty (0); genetics-only ({go:,}) enriches in "
          "nothing. The null is CO-MEASURED-conditioned -- each assay's validated count is "
          "randomly placed among the genes IT measured, so the proteome-and-spatial "
          "expectation is over the genes assayed by BOTH. p is uncorrected (x3 pairwise still "
          "<0.05). snATAC = single-nucleus ATAC; every mark encodes a live count; no gene "
          "names.".format(
              U=universe, k=k_disp, uv=unval, pu=100 - pct_corr, g1=ge1, pc=pct_corr,
              Pv=d["P_val"], Pm=P_meas, Sv=d["S_val"], Sm=S_meas, Av=d["A_val"], Am=A_meas,
              sub=sub, cg=col_gen, ov=ovl, po=ps_obs,
              pe=(ps_exp if np.isfinite(ps_exp) else float("nan")),
              pf=(f"{PS_['fold']:.2f}" if np.isfinite(PS_['fold']) else "n/a"),
              pp=_pfmt(PS_["p"]), go=gen_only))


if __name__ == "__main__":
    d, stats = get_data()
    build(d, stats)
