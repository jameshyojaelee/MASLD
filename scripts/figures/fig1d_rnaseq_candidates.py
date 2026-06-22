#!/usr/bin/env python3
"""
Fig 1d candidates — RNA-seq cohort-summary panel (pairs with the GWAS alluvial 1c).
Each candidate renders as its own PDF in figures/misc/ (user assembles in Illustrator).

Verified data (9 paper cohorts, 1,284 samples). Honesty: 392 unstaged-disease as a hatched
node/segment; staged subset 638/422; F4 rare (46); sex inferred for Govaere/Chen/Verschuren;
control = gray. scRNA geometry from the 269-donor / 895,542-cell / 5-lineage table (1.23M/12
= headline text only). Spearman(fibrosis, NAS) = 0.40 (n=638), verified live.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
from matplotlib.gridspec import GridSpec
import csv, os, math
from collections import defaultdict, Counter

plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size": 9, "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.spines.top": False, "axes.spines.right": False,
    "figure.dpi": 150, "savefig.dpi": 300, "figure.facecolor": "white",
    "axes.facecolor": "white", "axes.grid": False, "legend.frameon": False,
})

# ---- palette (shared with 1c) ----
CONTROL = "#9E9E9E"; HATCH_FACE = "#ECECEC"; COHORT_NODE = "#D4D4D4"
DX_ORDER = ["Control", "NAFL", "Borderline", "NASH", "NotAssessed"]
# NotAssessed = group_binary Disease but diagnosis_harmonized NaN (Chen + Verschuren stage by
# FIBROSIS not NAS; + ~26 Bril disease samples missing NAS) — MASLD patients lacking a NAS
# subtype, NOT missing/control.
# Colored disease-peach (hatched) so they group with disease; they still carry fibrosis stage.
DX_COLOR = {"Control": CONTROL, "NAFL": "#FEE0D2", "Borderline": "#FCBBA1",
            "NASH": "#FB8050", "NotAssessed": "#F6C3AE"}
DX_LABEL = {"Control": "Control", "NAFL": "MASL", "Borderline": "Borderline",
            "NASH": "MASH", "NotAssessed": "MASLD\n(NAS n.a.)"}
FIB_ORDER = ["F0", "F1", "F2", "F3", "F4", "NotStaged"]
FIB_COLOR = {"F0": "#EDF8E9", "F1": "#BAE4B3", "F2": "#A1D99B", "F3": "#74C476",
             "F4": "#238B45", "NotStaged": HATCH_FACE}
SEX_COLOR = {"F": "#AD1457", "M": "#1A237E"}
# scRNA coarse stages -> shared severity palette
SC_ORDER = ["Healthy", "Steatosis", "Steatohepatitis", "Cirrhosis", ""]
SC_COLOR = {"Healthy": CONTROL, "Steatosis": "#FEE0D2", "Steatohepatitis": "#FB8050",
            "Cirrhosis": "#B30000", "": HATCH_FACE}
SC_LABEL = {"Healthy": "Healthy", "Steatosis": "MASL", "Steatohepatitis": "MASH",
            "Cirrhosis": "Cirrhosis", "": "stage n.a."}
COHORT_PI = {"GSE126848": "GSE126848", "GSE130970": "GSE130970", "GSE135251": "GSE135251",
             "GSE162694": "GSE162694", "GSE167523": "GSE167523", "GSE174478": "GSE174478",
             "GSE193066": "GSE193066", "GSE213621": "GSE213621", "GSE240729": "GSE240729"}
PAPER = list(COHORT_PI.keys())


def _root(): return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
def _out(fn):
    d = os.path.join(_root(), "figures", "misc"); os.makedirs(d, exist_ok=True); return os.path.join(d, fn)
def _v(x): return x if x not in ("", "NA", "NaN", "nan", None) else None


def _spearman(xy):
    xs = [a for a, _ in xy]; ys = [b for _, b in xy]
    def rk(v):
        s = sorted(range(len(v)), key=lambda i: v[i]); r = [0.0] * len(v); i = 0
        while i < len(v):
            j = i
            while j + 1 < len(v) and v[s[j + 1]] == v[s[i]]: j += 1
            for k in range(i, j + 1): r[s[k]] = (i + j) / 2 + 1
            i = j + 1
        return r
    rx, ry = rk(xs), rk(ys); n = len(xs); mx = sum(rx) / n; my = sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    den = (sum((rx[i] - mx) ** 2 for i in range(n)) * sum((ry[i] - my) ** 2 for i in range(n))) ** 0.5
    return num / den if den else 0.0


def load_bulk():
    f = os.path.join(_root(), "RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv")
    rows = [r for r in csv.DictReader(open(f)) if r["dataset"] in PAPER]
    # Restrict to the analysis cohort (pass_technical), matching 03_integrate_counts.R — the
    # raw metadata is pre-QC (1,284); pass_technical keeps the QC-passing set (~1,259).
    qcf = os.path.join(_root(), "RNA-seq/Human/Patient_Cohorts/analysis/integration/qc/sample_qc_report.csv")
    passq = {q["sample_id"]: str(q.get("pass_technical", "")).strip().upper() in ("TRUE", "T", "1", "YES")
             for q in csv.DictReader(open(qcf))}
    rows = [r for r in rows if passq.get(r["sample_id"], False)]
    for r in rows:
        r["_dx"] = _v(r["diagnosis_harmonized"]) or "NotAssessed"
        fb = _v(r["fibrosis_stage"]); r["_fib"] = ("F" + str(int(float(fb)))) if fb is not None else "NotStaged"
        r["_fibn"] = int(float(fb)) if fb is not None else None
        nas = _v(r["nas_score"]); r["_nas"] = int(float(nas)) if nas is not None else None
        age = _v(r["age"]); r["_age"] = float(age) if age is not None else None
    return rows


def load_sex():
    f = _out("_sex_source.csv")
    by_sample, src = {}, defaultdict(set)
    if os.path.exists(f):
        for r in csv.DictReader(open(f)):
            s = _v(r.get("inferred_sex")) or _v(r.get("sex"))
            if s: by_sample[r["sample_id"]] = s[0].upper()
            src[r["dataset"]].add(r.get("sex_source", ""))
    return by_sample, src


def load_scrna():
    f = os.path.join(_root(), "Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv")
    return list(csv.DictReader(open(f), delimiter="\t"))


def _ribbon(ax, x0, y0b, y0t, x1, y1b, y1t, color, alpha=0.45, z=2):
    cx = (x0 + x1) / 2
    verts = [(x0, y0b), (cx, y0b), (cx, y1b), (x1, y1b), (x1, y1t), (cx, y1t), (cx, y0t), (x0, y0t), (x0, y0b)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(mpatches.PathPatch(Path(verts, codes), facecolor=color, edgecolor="none", alpha=alpha, zorder=z))


def _stack(items, top, gap):
    pos, y = {}, top
    for k, h in items:
        pos[k] = (y - h, y); y -= h + gap
    return pos


# ─────────────────────────────────────────────────────────────────────────────
# CONCEPT 1 — Clinical Cascade (cohort -> diagnosis -> fibrosis)
# ─────────────────────────────────────────────────────────────────────────────
def _cascade_core(ax, rows, y_top, height_frac=1.0, label=True, alpha=0.45):
    coh_n = Counter(r["dataset"] for r in rows)
    cohorts = sorted(PAPER, key=lambda c: -coh_n[c])
    coh_dx = defaultdict(Counter); dx_fib = defaultdict(Counter); dx_tot = Counter(); fib_tot = Counter()
    for r in rows:
        coh_dx[r["dataset"]][r["_dx"]] += 1; dx_fib[r["_dx"]][r["_fib"]] += 1
        dx_tot[r["_dx"]] += 1; fib_tot[r["_fib"]] += 1
    N = sum(coh_n.values()); GAP = N * 0.012
    H = N + GAP * (len(cohorts) - 1)
    coh_pos = _stack([(c, coh_n[c]) for c in cohorts], y_top, GAP)
    def span(items):
        g = (H - sum(h for _, h in items)) / max(len(items) - 1, 1); return _stack(items, y_top, g)
    dx_pos = span([(d, dx_tot[d]) for d in DX_ORDER if dx_tot[d] > 0])
    fib_pos = span([(f, fib_tot[f]) for f in FIB_ORDER if fib_tot[f] > 0])
    X0, X1, X2, NW = 0.0, 1.0, 2.0, 0.055

    def node(x, lo, hi, color, hatch=None):
        ax.add_patch(mpatches.Rectangle((x, lo), NW, hi - lo, facecolor=color,
                     edgecolor=("#BDBDBD" if hatch else "white"), linewidth=0.8, hatch=hatch, zorder=4))
    for c in cohorts:
        lo, hi = coh_pos[c]; node(X0, lo, hi, COHORT_NODE)
        if label: ax.text(X0 - 0.04, (lo + hi) / 2, f"{COHORT_PI[c]}  {coh_n[c]}", ha="right", va="center",
                          fontsize=7.4, color="#333", fontweight="bold")
    for d in dx_pos:
        lo, hi = dx_pos[d]; node(X1, lo, hi, DX_COLOR[d], hatch=("///" if d == "NotAssessed" else None))
        if label: ax.text(X1 + NW + 0.03, (lo + hi) / 2, f"{DX_LABEL[d]}  {dx_tot[d]}", ha="left", va="center",
                          fontsize=7.2, color=("#9C5A3C" if d == "NotAssessed" else "#333"), fontweight="bold")
    for f in fib_pos:
        lo, hi = fib_pos[f]; node(X2, lo, hi, FIB_COLOR[f], hatch=("///" if f == "NotStaged" else None))
        if label: ax.text(X2 + NW + 0.03, (lo + hi) / 2, f"{'not staged' if f=='NotStaged' else f}  {fib_tot[f]}",
                          ha="left", va="center", fontsize=7.4, color=("#999" if f == "NotStaged" else "#333"),
                          fontweight=("normal" if f == "NotStaged" else "bold"))
    cc = {c: coh_pos[c][1] for c in cohorts}; di = {d: dx_pos[d][1] for d in dx_pos}
    for c in cohorts:
        for d in DX_ORDER:
            h = coh_dx[c][d]
            if h and d in dx_pos:
                a = cc[c]; cc[c] -= h; b = di[d]; di[d] -= h
                _ribbon(ax, X0 + NW, a - h, a, X1, b - h, b, DX_COLOR[d], alpha)
    do = {d: dx_pos[d][1] for d in dx_pos}; fi = {f: fib_pos[f][1] for f in fib_pos}
    for d in DX_ORDER:
        if d not in dx_pos: continue
        for f in FIB_ORDER:
            h = dx_fib[d][f]
            if h and f in fib_pos:
                a = do[d]; do[d] -= h; b = fi[f]; fi[f] -= h
                _ribbon(ax, X1 + NW, a - h, a, X2, b - h, b, DX_COLOR[d], alpha)
    return H, GAP, (X0, X1, X2, NW)


def build_cascade(path):
    rows = load_bulk()
    fig, ax = plt.subplots(figsize=(7.4, 5.6)); fig.patch.set_facecolor("white")
    H, GAP, (X0, X1, X2, NW) = _cascade_core(ax, rows, sum(Counter(r["dataset"] for r in rows).values()) + 0)
    H = sum(1 for _ in rows)  # not used for limits below; recompute via core return
    # column captions
    N = len(rows); GAP = N * 0.012; Htot = N + GAP * 8
    for x, t in [(0.027, "9 cohorts (n)"), (1.027, "diagnosis"), (2.027, "fibrosis stage")]:
        ax.text(x, -GAP * 3.5, t, ha="center", va="top", fontsize=7.5, color="#555", style="italic")
    ax.set_xlim(-0.95, 2.9); ax.set_ylim(-GAP * 6, Htot + GAP * 2); ax.axis("off")
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white", pad_inches=0.03); plt.close(fig)
    print("Saved", os.path.basename(path))


def build_cascade_scrna(path):
    """Cascade + a thin desaturated scRNA mini-flow band on top."""
    rows = load_bulk(); N = len(rows); GAP = N * 0.012; Htot = N + GAP * 8
    fig, ax = plt.subplots(figsize=(7.4, 6.6)); fig.patch.set_facecolor("white")
    _cascade_core(ax, rows, Htot)
    # scRNA band above: donors -> coarse stage, ribbon mass proportional to CELLS (scaled to bulk width)
    sc = load_scrna()
    cells = Counter();
    for r in sc: cells[r["disease_stage_coarse"]] += int(float(r["n_cells"]))
    tot_cells = sum(cells.values()); band_h = N * 0.5  # scale cells to ~half the bulk height
    band_y = Htot + GAP * 9
    src_lo = band_y; src_hi = band_y + band_h
    ax.add_patch(mpatches.Rectangle((0.0, src_lo), 0.055, band_h, facecolor="#BBBBBB", edgecolor="white", zorder=4))
    ax.text(-0.04, (src_lo + src_hi) / 2, "scRNA\n269 donors\n0.9M cells", ha="right", va="center",
            fontsize=6.8, color="#666")
    # stage nodes at X1
    items = [(s, cells[s]) for s in SC_ORDER if cells.get(s, 0) > 0]
    g = (band_h - sum(h for _, h in items) / tot_cells * band_h) / max(len(items) - 1, 1)
    y = src_hi; spos = {}
    for s, c in items:
        h = c / tot_cells * band_h; spos[s] = (y - h, y); y -= h + g
    cur_src = src_hi
    for s, c in items:
        h = c / tot_cells * band_h; lo, hi = spos[s]
        ax.add_patch(mpatches.Rectangle((1.0, lo), 0.055, hi - lo, facecolor=SC_COLOR[s],
                     edgecolor=("#BDBDBD" if s == "" else "white"), hatch=("///" if s == "" else None), zorder=4))
        _ribbon(ax, 0.055, cur_src - h, cur_src, 1.0, lo, hi, SC_COLOR[s] if s != "" else "#CFCFCF", 0.4)
        cur_src -= h
        ax.text(1.0 + 0.055 + 0.03, (lo + hi) / 2, SC_LABEL[s], ha="left", va="center", fontsize=6.6, color="#777")
    for x, t in [(0.027, "9 cohorts (n)"), (1.027, "diagnosis"), (2.027, "fibrosis stage")]:
        ax.text(x, -GAP * 3.5, t, ha="center", va="top", fontsize=7.5, color="#555", style="italic")
    ax.set_xlim(-0.95, 2.9); ax.set_ylim(-GAP * 6, src_hi + GAP * 3); ax.axis("off")
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white", pad_inches=0.03); plt.close(fig)
    print("Saved", os.path.basename(path))


# ─────────────────────────────────────────────────────────────────────────────
# CONCEPT 2 — Severity Ridge (fibrosis | NAS ridgelines per cohort)
# ─────────────────────────────────────────────────────────────────────────────
def build_ridge(path):
    rows = load_bulk()
    by = defaultdict(list)
    for r in rows: by[r["dataset"]].append(r)
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 5.2)); fig.patch.set_facecolor("white")
    for ax, (field, vmax, title) in zip(axes, [("_fibn", 4, "Fibrosis stage"), ("_nas", 8, "NAS score")]):
        present = [(c, [r[field] for r in by[c] if r[field] is not None]) for c in PAPER]
        present = [(c, v) for c, v in present if v]
        absent = [c for c in PAPER if c not in [p[0] for p in present]]
        present.sort(key=lambda cv: sum(cv[1]) / len(cv[1]))     # by mean severity, low at bottom
        rowh = 1.0
        for i, (c, vals) in enumerate(present):
            base = i * rowh
            hist = [0] * (vmax + 1)
            for v in vals: hist[v] += 1
            mx = max(hist) or 1
            xs = list(range(vmax + 1)); ys = [base + h / mx * rowh * 1.7 for h in hist]
            med = sum(vals) / len(vals)
            col = plt.get_cmap("OrRd")(0.25 + 0.6 * med / vmax)
            ax.fill_between(xs, base, ys, color=col, alpha=0.7, zorder=10 - i * 0.0, lw=0)
            ax.plot(xs, ys, color="white", lw=0.8, zorder=11)
            ax.text(-0.35, base + 0.15, f"{COHORT_PI[c]} ({len(vals)})", ha="right", va="bottom", fontsize=7, color="#333")
        # absent cohorts as a hatched 'not assessed' lane at top
        if absent:
            n_abs = sum(len(by[c]) for c in absent)
            ytop = len(present) * rowh
            ax.add_patch(mpatches.Rectangle((-0.0, ytop + 0.25), vmax, 0.4, facecolor=HATCH_FACE,
                         edgecolor="#BDBDBD", hatch="///", zorder=5))
            ax.text(vmax / 2, ytop + 0.78,
                    f"not assessed (n={n_abs}): " + ", ".join(COHORT_PI[c] for c in absent),
                    ha="center", va="bottom", fontsize=6.2, color="#666")
        ax.set_xlim(-0.6, vmax + 0.4); ax.set_ylim(-0.3, len(present) * rowh + 1.4)
        ax.set_yticks([]); ax.set_xticks(range(vmax + 1)); ax.tick_params(labelsize=7.5, bottom=False)
        ax.set_xlabel(title, fontsize=8.5)
        for sp in ["top", "right", "left"]: ax.spines[sp].set_visible(False)
        ax.spines["bottom"].set_linewidth(0.4)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white", pad_inches=0.03); plt.close(fig)
    print("Saved", os.path.basename(path))


# ─────────────────────────────────────────────────────────────────────────────
# CONCEPT 3 — Clinical Phase-Space (fibrosis x NAS bubble grid + marginals)
# ─────────────────────────────────────────────────────────────────────────────
def build_phasespace(path):
    rows = load_bulk(); sex_by, _ = load_sex()
    paired = [r for r in rows if r["_fibn"] is not None and r["_nas"] is not None]
    cell = defaultdict(Counter)                              # (fib,nas) -> dx counts
    for r in paired: cell[(r["_fibn"], r["_nas"])][r["_dx"]] += 1
    fig = plt.figure(figsize=(6.6, 6.2)); fig.patch.set_facecolor("white")
    gs = GridSpec(2, 2, width_ratios=[4, 1.1], height_ratios=[1.1, 4], hspace=0.04, wspace=0.04,
                  left=0.1, right=0.97, top=0.97, bottom=0.09)
    axc = fig.add_subplot(gs[1, 0]); axt = fig.add_subplot(gs[0, 0], sharex=axc); axr = fig.add_subplot(gs[1, 1], sharey=axc)
    mx = max(sum(c.values()) for c in cell.values())
    for (fb, na), c in cell.items():
        n = sum(c.values()); modal = c.most_common(1)[0][0]
        axc.scatter([fb], [na], s=20 + 360 * n / mx, color=DX_COLOR[modal], alpha=0.82,
                    edgecolor="white", lw=0.6, zorder=3)
    # ghost row: fibrosis-only (no NAS) below the grid
    fibonly = Counter(r["_fibn"] for r in rows if r["_fibn"] is not None and r["_nas"] is None)
    gmx = max(fibonly.values()) if fibonly else 1
    for fb, n in fibonly.items():
        axc.scatter([fb], [-1.1], s=20 + 360 * n / mx, color="#D9D9D9", alpha=0.7, edgecolor="#BDBDBD", lw=0.6, zorder=2)
    axc.axhline(-0.5, color="#E0E0E0", lw=0.7, zorder=1)
    axc.text(-0.7, -1.1, f"fibrosis-staged,\nno NAS  (n={sum(fibonly.values())})",
             ha="right", va="center", fontsize=5.8, color="#999", style="italic")
    axc.set_xlim(-0.8, 4.5); axc.set_ylim(-1.8, 8.5)
    axc.set_xticks(range(5)); axc.set_yticks(range(9)); axc.tick_params(labelsize=7.5)
    axc.set_xlabel("Fibrosis stage", fontsize=8.5); axc.set_ylabel("NAS score", fontsize=8.5)
    for sp in ["top", "right"]: axc.spines[sp].set_visible(False)
    rho = _spearman([(r["_fibn"], r["_nas"]) for r in paired])
    axc.text(0.03, 0.97, f"ρ = {rho:.2f}   n = {len(paired)}", transform=axc.transAxes, fontsize=7.5,
             color="#555", va="top")
    # top marginal: fibrosis x diagnosis stacked
    fibdx = defaultdict(Counter)
    for r in paired: fibdx[r["_fibn"]][r["_dx"]] += 1
    for fb in range(5):
        base = 0
        for d in DX_ORDER:
            h = fibdx[fb][d]
            if h: axt.bar(fb, h, bottom=base, color=DX_COLOR[d], width=0.8, edgecolor="white", lw=0.4); base += h
    axt.set_yticks([]); axt.tick_params(bottom=False, labelbottom=False)
    for sp in ["top", "right", "left"]: axt.spines[sp].set_visible(False)
    # right marginal: NAS x sex stacked
    nassex = defaultdict(Counter)
    for r in paired: nassex[r["_nas"]][sex_by.get(r["sample_id"], "?")] += 1
    for na in range(9):
        base = 0
        for s in ["F", "M"]:
            h = nassex[na][s]
            if h: axr.barh(na, h, left=base, color=SEX_COLOR[s], height=0.8, edgecolor="white", lw=0.4, alpha=0.85); base += h
    axr.set_xticks([]); axr.tick_params(left=False, labelleft=False)
    for sp in ["top", "right", "bottom"]: axr.spines[sp].set_visible(False)
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white", pad_inches=0.03); plt.close(fig)
    print("Saved", os.path.basename(path))


# ─────────────────────────────────────────────────────────────────────────────
# CONCEPT 4 — Two-Lane Cascade (bulk + scRNA on one severity axis)
# ─────────────────────────────────────────────────────────────────────────────
def build_twolane(path):
    rows = load_bulk(); sc = load_scrna()
    # shared 5-stop severity axis
    stops = ["Control", "MASL", "Borderline", "MASH", "Cirrhosis"]
    bulk = Counter()
    bmap = {"Control": "Control", "NAFL": "MASL", "Borderline": "Borderline", "NASH": "MASH"}
    for r in rows:
        if r["_dx"] in bmap: bulk[bmap[r["_dx"]]] += 1
    bulk_na = sum(1 for r in rows if r["_dx"] == "NotAssessed")
    sc_cells = Counter(); smap = {"Healthy": "Control", "Steatosis": "MASL", "Steatohepatitis": "MASH", "Cirrhosis": "Cirrhosis"}
    for r in sc:
        st = r["disease_stage_coarse"]
        if st in smap: sc_cells[smap[st]] += int(float(r["n_cells"]))
    sc_na_cells = sum(int(float(r["n_cells"])) for r in sc if r["disease_stage_coarse"] not in smap)

    fig, ax = plt.subplots(figsize=(8.0, 3.6)); fig.patch.set_facecolor("white")
    xpos = {s: i for i, s in enumerate(stops)}
    SEVC = {"Control": CONTROL, "MASL": "#FEE0D2", "Borderline": "#FCBBA1", "MASH": "#FB8050", "Cirrhosis": "#B30000"}
    # top lane: bulk barh (counts)
    bmax = max(bulk.values())
    for s in stops:
        if bulk[s]: ax.add_patch(mpatches.FancyBboxPatch((xpos[s] - 0.38, 1.15), 0.76, 0.7 * bulk[s] / bmax + 0.05,
                    boxstyle="round,pad=0.005", facecolor=SEVC[s], edgecolor="white", lw=0.8, alpha=0.85))
        if bulk[s]: ax.text(xpos[s], 1.15 + 0.7 * bulk[s] / bmax + 0.10, f"{bulk[s]}", ha="center", va="bottom", fontsize=7, color="#555")
    # bottom lane: scRNA dots, area ∝ cells
    cmax = max(sc_cells.values())
    for s in stops:
        if sc_cells.get(s, 0):
            ax.scatter([xpos[s]], [0.5], s=120 + 1500 * sc_cells[s] / cmax, color=SEVC[s], alpha=0.8, edgecolor="white", lw=0.8)
            ax.text(xpos[s], 0.18, f"{sc_cells[s]//1000}k", ha="center", va="top", fontsize=6.5, color="#777")
    # NA markers at right
    ax.add_patch(mpatches.Rectangle((len(stops) - 0.4, 1.15), 0.76, 0.7 * bulk_na / bmax + 0.05,
                 facecolor=DX_COLOR["NotAssessed"], edgecolor="#BDBDBD", hatch="///"))
    ax.text(len(stops), 1.15 + 0.7 * bulk_na / bmax + 0.10, f"MASLD\nNAS n.a.\n{bulk_na}", ha="center", va="bottom", fontsize=6, color="#9C5A3C")
    ax.scatter([len(stops)], [0.5], s=120 + 1500 * sc_na_cells / cmax, color=HATCH_FACE, edgecolor="#BDBDBD", lw=0.8, hatch="///")
    ax.text(-0.7, 1.5, f"bulk\n{len(rows):,} pts\n(QC-pass)", ha="right", va="center", fontsize=7.5, color="#333", fontweight="bold")
    ax.text(-0.7, 0.5, "scRNA\n0.9M cells", ha="right", va="center", fontsize=7.5, color="#333", fontweight="bold")
    for s in stops: ax.text(xpos[s], -0.15, s, ha="center", va="top", fontsize=7.8, fontweight="bold", color="#333")
    ax.text(len(stops), -0.15, "n.a.", ha="center", va="top", fontsize=7.5, color="#999")
    ax.set_xlim(-1.4, len(stops) + 0.8); ax.set_ylim(-0.5, 2.1); ax.axis("off")
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white", pad_inches=0.03); plt.close(fig)
    print("Saved", os.path.basename(path))


# ─────────────────────────────────────────────────────────────────────────────
# CONCEPT 5 — Marimekko (width ∝ cohort n, height = diagnosis proportion)
# ─────────────────────────────────────────────────────────────────────────────
def build_marimekko(path):
    rows = load_bulk()
    coh_n = Counter(r["dataset"] for r in rows); coh_dx = defaultdict(Counter)
    for r in rows: coh_dx[r["dataset"]][r["_dx"]] += 1
    cohorts = sorted(PAPER, key=lambda c: -coh_n[c]); N = sum(coh_n.values())
    fig, ax = plt.subplots(figsize=(7.6, 4.2)); fig.patch.set_facecolor("white")
    x = 0.0; gap = N * 0.004
    for c in cohorts:
        w = coh_n[c]; base = 0.0
        for d in DX_ORDER:
            h = coh_dx[c][d] / coh_n[c]
            if h <= 0: continue
            ax.add_patch(mpatches.Rectangle((x, base), w, h, facecolor=DX_COLOR[d],
                         edgecolor="white", lw=0.8, hatch=("///" if d == "NotAssessed" else None), zorder=2))
            base += h
        ax.text(x + w / 2, -0.03, f"{COHORT_PI[c]}\n{coh_n[c]}", ha="center", va="top", fontsize=7, color="#333")
        x += w + gap
    # legend
    leg = [mpatches.Patch(facecolor=DX_COLOR[d], label=DX_LABEL[d],
           hatch=("///" if d == "NotAssessed" else None)) for d in DX_ORDER]
    ax.legend(handles=leg, fontsize=7.5, loc="lower center", bbox_to_anchor=(0.5, -0.22), ncol=5, frameon=False)
    ax.set_xlim(-N * 0.01, x); ax.set_ylim(-0.12, 1.02); ax.axis("off")
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white", pad_inches=0.03); plt.close(fig)
    print("Saved", os.path.basename(path))


# ─────────────────────────────────────────────────────────────────────────────
# Add-on — Covariate-availability gutter
# ─────────────────────────────────────────────────────────────────────────────
def build_gutter(path):
    rows = load_bulk(); sex_by, src = load_sex()
    by = defaultdict(list)
    for r in rows: by[r["dataset"]].append(r)
    cohorts = sorted(PAPER, key=lambda c: -len(by[c]))
    cols = ["Sex", "Age", "Fibrosis", "NAS"]
    fig, ax = plt.subplots(figsize=(4.4, 4.4)); fig.patch.set_facecolor("white")
    n = len(cohorts)
    for i, c in enumerate(cohorts):
        y = n - 1 - i; rs = by[c]
        avail = {
            "Sex": ("inferred" if ("inferred_kmeans" in src.get(c, set())) else
                    ("yes" if any(r["sample_id"] in sex_by for r in rs) else "no")),
            "Age": "yes" if any(r["_age"] is not None for r in rs) else "no",
            "Fibrosis": "yes" if any(r["_fibn"] is not None for r in rs) else "no",
            "NAS": "yes" if any(r["_nas"] is not None for r in rs) else "no",
        }
        ax.text(-0.3, y, f"{COHORT_PI[c]} ({len(rs)})", ha="right", va="center", fontsize=8, color="#333", fontweight="bold")
        for j, col in enumerate(cols):
            st = avail[col]
            if st == "yes":
                ax.add_patch(mpatches.Circle((j, y), 0.27, facecolor="#4F7DB3", edgecolor="white", lw=0.8, zorder=3))
            elif st == "inferred":
                ax.add_patch(mpatches.Circle((j, y), 0.27, facecolor="#4F7DB3", edgecolor="white", lw=0.8, hatch="///", zorder=3))
            else:
                ax.add_patch(mpatches.Circle((j, y), 0.16, facecolor="none", edgecolor="#C0C0C0", lw=1.2, zorder=3))
    for j, col in enumerate(cols):
        ax.text(j, n - 0.4, col, ha="center", va="bottom", fontsize=8, color="#333", fontweight="bold")
    # staged-subset counts computed live (annotated sex only, since inferred cohorts have no per-sample annotation)
    staged = [r for r in rows if r["_fibn"] is not None and r["_nas"] is not None]
    ann_cohorts = {d for d, ss in src.items() if "inferred_kmeans" not in ss}
    staged_annsex = [r for r in staged if r["dataset"] in ann_cohorts and r["sample_id"] in sex_by]
    # legend + staged bracket (separate lines, no special glyphs)
    ax.text(-1.85, -0.9, "filled = measured   ·   hatched = inferred   ·   open = absent",
            ha="left", va="center", fontsize=6.6, color="#666")
    ax.text(-1.85, -1.35, f"staged subset: {len(staged)} (fib + NAS)  ·  {len(staged_annsex)} (+ annotated sex)",
            ha="left", va="center", fontsize=6.6, color="#666", style="italic")
    ax.set_xlim(-2.0, len(cols) - 0.3); ax.set_ylim(-1.7, n + 0.2); ax.set_aspect("equal"); ax.axis("off")
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white", pad_inches=0.03); plt.close(fig)
    print("Saved", os.path.basename(path))


if __name__ == "__main__":
    build_cascade(_out("fig1d_cascade.pdf"))
    build_cascade_scrna(_out("fig1d_cascade_scrna.pdf"))
    build_ridge(_out("fig1d_ridge.pdf"))
    build_phasespace(_out("fig1d_phasespace.pdf"))
    build_twolane(_out("fig1d_twolane.pdf"))
    build_marimekko(_out("fig1d_marimekko.pdf"))
    build_gutter(_out("fig1d_gutter.pdf"))
    print("all done")
