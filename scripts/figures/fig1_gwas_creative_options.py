#!/usr/bin/env python3
"""
Fig 2A GWAS Portfolio — alluvial, evidence cascade, and per-study bars (breadth, no N encoding).
KEY MESSAGE: 50 colocalising GWAS spanning 5 ancestries (EUR/AFR/AMR/EAS/SAS, MVP included)
x liver-disease traits — show coverage breadth. All counts are DATA-DRIVEN at runtime from
the GWAS registry, per-study fine-mapping summary, and per-gene colocalisation tables (no
hardcoded manifest).

Canonical Fig 2A panel outputs (figures/main/fig2_genetics/panels/):
  fig2A_gwas_alluvial.pdf     (build_B)
  fig2A_gwas_cascade.pdf      (build_cascade)
  FigS2A_gwas_perstudy.pdf    (build_study_bars)
Also writes exploratory options A–F to figures/misc/.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.path import Path
import numpy as np
import os
from collections import defaultdict

# ── rcParams ──────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.sans-serif":   ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":         9,
    "pdf.fonttype":      42,
    "ps.fonttype":       42,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "figure.dpi":        150,
    "savefig.dpi":       300,
    "figure.facecolor":  "white",
    "axes.facecolor":    "white",
    "axes.grid":         False,
    "legend.frameon":    False,
})

# ── Data (DATA-DRIVEN: read GWAS registry + fine-mapping + coloc at runtime) ────
# The portfolio (50 colocalising GWAS across 5 ancestries, incl. the MVP / Million
# Veteran Program multi-ancestry strata) is NOT hardcoded here — every count is
# recomputed at runtime from the three canonical GWAS files below.
import csv as _csv
import math as _math


def _project_root():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _out_path(fname):
    d = os.path.join(_project_root(), "figures", "misc")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, fname)


_ROOT          = os.environ.get("MASLD_PROJECT_ROOT", _project_root())
_REGISTRY_TSV  = os.path.join(_ROOT, "GWAS/finemapping/config/gwas_registry.tsv")
_STUDY_SUMMARY = os.path.join(_ROOT, "GWAS/finemapping/results/study_summary.csv")
_GENE_COLOC    = os.path.join(_ROOT, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")

# Trait derived from a study/gwas name — longer/specific tokens first so e.g.
# "ChronLiver"/"Albumin" match before the enzyme tokens.
_TRAIT_TOKENS = [
    ("ChronLiver", "Chronic liver disease"), ("Cirrhosis", "Cirrhosis"),
    ("Albumin", "Albumin"), ("Platelet", "Platelet"), ("PDFF", "PDFF"),
    ("NAFLD", "NAFLD"), ("NASH", "NASH"), ("HCC", "HCC"),
    ("ALT", "ALT"), ("AST", "AST"), ("GGT", "GGT"),
]


def _trait_of(name):
    low = name.lower()
    for tok, lab in _TRAIT_TOKENS:
        if tok.lower() in low:
            return lab
    return None


def _fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def _load_manifest():
    """Read the 50-GWAS registry, per-study fine-mapping yield, and per-gene coloc.
    Returns data-driven aggregates used by every panel — no frozen literals."""
    # 1) registry -> per-(ancestry, trait) GWAS-dataset counts
    reg = list(_csv.DictReader(open(_REGISTRY_TSV), delimiter="\t"))
    anc_of_study = {r["study_name"]: r["ancestry"] for r in reg}
    at_gwas, anc_gwas = defaultdict(int), defaultdict(int)
    studies = set()
    for r in reg:
        a, t = r["ancestry"], _trait_of(r["study_name"])
        at_gwas[(a, t)] += 1
        anc_gwas[a] += 1
        # distinct-study count collapses the Sveinbjornsson 2023 (PMID 36280732)
        # deCODE/Intermountain/UKBB cohort-arms of one study into a single study.
        base = "2023_36280732_NAFLD" if r["study_name"].startswith("2023_36280732") else r["study_name"]
        studies.add(base)

    # 2) study_summary -> per-study & per-(ancestry, trait) fine-mapped loci
    at_loci, study_loci = defaultdict(int), {}
    for r in _csv.DictReader(open(_STUDY_SUMMARY)):
        s = r["study"]
        a = anc_of_study.get(s)
        if a is None:
            continue
        nl = int(r["n_loci"])
        at_loci[(a, _trait_of(s))] += nl
        study_loci[s] = nl

    # 3) gene_level_coloc -> union (SuSiE OR ABF, PP.H4>0.5) colocalising genes.
    # Each gene belongs to the ancestry/trait of its single best colocalisation.
    at_genes, anc_genes, study_genes = defaultdict(int), defaultdict(int), defaultdict(int)
    union_total = susie_total = 0
    for r in _csv.DictReader(open(_GENE_COLOC)):
        abf, su = _fnum(r["coloc_best_pp4"]), _fnum(r["coloc_best_susie_pp4"])
        abf_ok, su_ok = abf > 0.5, su > 0.5
        if su_ok:
            susie_total += 1
        if not (abf_ok or su_ok):
            continue
        union_total += 1
        if su_ok and (_math.isnan(abf) or su >= abf):
            a, gw = r["coloc_best_susie_ancestry"], r["coloc_best_susie_gwas"]
        else:
            a, gw = r["coloc_best_ancestry"], r["coloc_best_gwas"]
        at_genes[(a, _trait_of(gw))] += 1
        anc_genes[a] += 1
        study_genes[gw] += 1

    return dict(reg=reg, anc_of_study=anc_of_study, at_gwas=at_gwas, anc_gwas=anc_gwas,
                n_gwas=len(reg), n_studies=len(studies), at_loci=at_loci, study_loci=study_loci,
                at_genes=at_genes, anc_genes=anc_genes, study_genes=study_genes,
                union_total=union_total, susie_total=susie_total)


_M          = _load_manifest()
N_GWAS      = _M["n_gwas"]        # 50 colocalising GWAS datasets
N_STUDIES   = _M["n_studies"]     # 48 distinct studies (Sveinbjornsson arms collapsed)
UNION_GENES = _M["union_total"]   # 1527 union SuSiE-or-ABF PP.H4>0.5 coloc genes
SUSIE_GENES = _M["susie_total"]   # 736 SuSiE PP.H4>0.5 coloc genes

# Ancestry order = descending GWAS count (EUR, AFR, EAS, AMR, SAS).
ANCESTRIES     = sorted(_M["anc_gwas"], key=lambda a: (-_M["anc_gwas"][a], a))
ANCESTRY_LABEL = {"EUR": "European", "AFR": "African", "AMR": "Admixed American",
                  "EAS": "East Asian", "SAS": "South Asian"}
ANCESTRY_N     = {a: _M["anc_gwas"][a] for a in ANCESTRIES}
# Canonical ancestry palette shared across the Fig 2 genetics panels (colorblind-safe).
ANCESTRY_COLOR = {"EUR": "#4C72B0", "AFR": "#55A868", "AMR": "#DD8452",
                  "EAS": "#C44E52", "SAS": "#8172B3"}

# Traits present in the portfolio, grouped disease -> imaging -> enzymes -> biomarkers.
TRAITS      = ["NAFLD", "NASH", "Cirrhosis", "Chronic liver disease", "PDFF",
               "ALT", "AST", "GGT", "Albumin", "Platelet"]
TRAIT_SHORT = {"Chronic liver disease": "Chr. liver"}
# Trait group: (label, col_start, col_end_excl, band_color) — greyscale (colour = ancestry).
TRAIT_GROUPS = [
    ("Disease diagnoses", 0, 4,  "#444444"),
    ("Imaging",           4, 5,  "#7F7F7F"),
    ("Liver enzymes",     5, 8,  "#C0C0C0"),
    ("Other biomarkers",  8, 10, "#9C9C9C"),
]
TRAIT_GROUP_OF = {}
for _lbl, _c0, _c1, _gc in TRAIT_GROUPS:
    for _t in TRAITS[_c0:_c1]:
        TRAIT_GROUP_OF[_t] = _gc

ABSENT_COLOR = "#E8E8E8"


def _cell_data():
    """Per-(ancestry, trait) GWAS-dataset counts, data-driven from the registry."""
    mat = defaultdict(lambda: defaultdict(int))
    for (a, t), n in _M["at_gwas"].items():
        mat[a][t] = n
    return mat


def _group_bands(ax, alpha=0.07):
    for _, c0, c1, color in TRAIT_GROUPS:
        ax.axvspan(c0 - 0.5, c1 - 0.5, color=color, alpha=alpha, zorder=0, linewidth=0)


def _group_top_labels(ax):
    ax2 = ax.twiny()
    ax2.set_xlim(ax.get_xlim())
    ax2.set_xticks([(c0 + c1 - 1) / 2 for _, c0, c1, _ in TRAIT_GROUPS])
    ax2.set_xticklabels([lbl for lbl, *_ in TRAIT_GROUPS],
                        fontsize=7, color="#666666", style="italic")
    ax2.tick_params(top=False, labeltop=True, left=False, right=False, bottom=False)
    for sp in ax2.spines.values():
        sp.set_visible(False)


def _anc_strips(ax, xoff=-0.55, w=0.07):
    for i, anc in enumerate(ANCESTRIES):
        ax.add_patch(mpatches.Rectangle(
            (xoff, i - 0.45), w, 0.9,
            transform=ax.transData, color=ANCESTRY_COLOR[anc],
            clip_on=False, zorder=6,
        ))


def _grid_axis(ax):
    n_trait = len(TRAITS)
    ax.set_xlim(-0.5, n_trait - 0.5)
    ax.set_ylim(-0.5, len(ANCESTRIES) - 0.5)
    ax.set_xticks(range(n_trait))
    ax.set_xticklabels([TRAIT_SHORT.get(t, t) for t in TRAITS], fontsize=6, rotation=45, ha="right")
    ax.set_yticks(range(len(ANCESTRIES)))
    ax.set_yticklabels([ANCESTRY_LABEL[a] for a in ANCESTRIES], fontsize=8.5)
    ax.tick_params(left=False, bottom=False)
    for sp in ax.spines.values():
        sp.set_visible(False)


# ─────────────────────────────────────────────────────────────────────────────
# A: Binary presence grid
# ─────────────────────────────────────────────────────────────────────────────
def build_A(path):
    """Binary grid: cell = ancestry color if covered, gray if absent. Count annotated."""
    cd = _cell_data()
    fig, ax = plt.subplots(figsize=(8.5, 3.0))
    fig.patch.set_facecolor("white")
    _group_bands(ax, alpha=0.05)

    for i, anc in enumerate(ANCESTRIES):
        for j, trait in enumerate(TRAITS):
            n = cd[anc][trait]
            color = ANCESTRY_COLOR[anc] if n > 0 else ABSENT_COLOR
            ax.add_patch(mpatches.FancyBboxPatch(
                (j - 0.43, i - 0.43), 0.86, 0.86,
                boxstyle="round,pad=0.04",
                facecolor=color, edgecolor="white", linewidth=1.0,
                alpha=0.80 if n > 0 else 1.0, zorder=2,
            ))
            if n > 1:
                ax.text(j, i, str(n), ha="center", va="center",
                        fontsize=6, color="white", zorder=3)

    for x in np.arange(-0.5, len(TRAITS), 1):
        ax.axvline(x, color="white", linewidth=1.5, zorder=1)
    for y in np.arange(-0.5, len(ANCESTRIES), 1):
        ax.axhline(y, color="white", linewidth=1.5, zorder=1)

    _grid_axis(ax)
    _anc_strips(ax)
    _group_top_labels(ax)

    leg = [mpatches.Patch(facecolor=ANCESTRY_COLOR[a], alpha=0.80,
                          label=f"{ANCESTRY_LABEL[a]} ({ANCESTRY_N[a]} studies)")
           for a in ANCESTRIES]
    leg.append(mpatches.Patch(facecolor=ABSENT_COLOR, label="Not covered"))
    ax.legend(handles=leg, fontsize=7, loc="lower right",
              bbox_to_anchor=(1.0, -0.26), ncol=3, frameon=False)

    ax.set_title("Option A — Binary coverage grid  |  number = study count  |  color = ancestry",
                 fontsize=6, pad=16, loc="left", color="black")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Saved: {path}")


# ─────────────────────────────────────────────────────────────────────────────
# B: Alluvial flow  (1 unit = 1 GWAS study)
# ─────────────────────────────────────────────────────────────────────────────
def _ribbon(ax, x0, y0b, y0t, x1, y1b, y1t, color, alpha=0.40):
    cx = (x0 + x1) / 2
    verts = [
        (x0, y0b), (cx, y0b), (cx, y1b), (x1, y1b),
        (x1, y1t), (cx, y1t), (cx, y0t), (x0, y0t),
        (x0, y0b),
    ]
    codes = [Path.MOVETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO,
             Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.CLOSEPOLY]
    ax.add_patch(mpatches.PathPatch(
        Path(verts, codes),
        facecolor=color, edgecolor="none", alpha=alpha, zorder=2,
    ))


def build_B(path, extra_paths=None):
    """Alluvial: unit width = 1 GWAS study. Left = ancestries, Right = traits."""
    cd = _cell_data()
    GAP = 0.40   # reduced from 0.55 to make layout more compact
    BW  = 0.055   # block width

    # Trait order: disease -> imaging -> enzymes -> biomarkers (TRAITS is already in this order)
    trait_total = {t: sum(cd[a][t] for a in ANCESTRIES) for t in TRAITS}
    right_order = list(TRAITS)

    # Ancestry blocks (left, top-to-bottom, y increasing downward after invert_yaxis)
    left_start = {}
    y = 0
    for anc in ANCESTRIES:
        left_start[anc] = y
        y += ANCESTRY_N[anc] + GAP
    total_left = y - GAP

    # Trait blocks (right, same order as right_order)
    right_start = {}
    y = 0
    for t in right_order:
        right_start[t] = y
        y += trait_total[t] + GAP
    total_right = y - GAP

    total_h = max(total_left, total_right) + 0.4
    TAG_MARGIN = 0.1

    fig, ax = plt.subplots(figsize=(7.0, 5.0))  # compact: 6pt text vs canvas
    fig.patch.set_facecolor("white")
    ax.set_xlim(-0.28, 1.28)
    ax.set_ylim(0, total_h + TAG_MARGIN)
    ax.invert_yaxis()
    ax.axis("off")

    # Left ancestry blocks + labels
    for anc in ANCESTRIES:
        y0 = left_start[anc]
        h  = ANCESTRY_N[anc]
        ax.add_patch(mpatches.Rectangle(
            (-BW / 2, y0), BW, h,
            facecolor=ANCESTRY_COLOR[anc], edgecolor="none", alpha=1.0, zorder=4,
        ))
        ax.text(-BW / 2 - 0.02, y0 + h / 2,
                f"{ANCESTRY_LABEL[anc]}\n({ANCESTRY_N[anc]})",
                ha="right", va="center", fontsize=6,
                color="black")

    # Right trait blocks + labels — each block coloured by its trait group (greyscale;
    # colour is reserved for ancestry).
    trait_colors = dict(TRAIT_GROUP_OF)

    for t in right_order:
        y0 = right_start[t]
        h  = trait_total[t]
        ax.add_patch(mpatches.Rectangle(
            (1.0 - BW / 2, y0), BW, h,
            facecolor=trait_colors[t], edgecolor="none", alpha=1.0, zorder=4,
        ))
        ax.text(1.0 + BW / 2 + 0.02, y0 + h / 2,
                f"{t}  ({h})", ha="left", va="center",
                fontsize=6, color="black")

    # Ribbons: iterate traits in right_order within each ancestry block
    l_off = {anc: left_start[anc] for anc in ANCESTRIES}
    r_off = {t: right_start[t] for t in right_order}

    for anc in ANCESTRIES:
        for t in right_order:
            n = cd[anc][t]
            if n == 0:
                continue
            _ribbon(ax, 0.0, l_off[anc], l_off[anc] + n,
                        1.0, r_off[t],   r_off[t] + n,
                    ANCESTRY_COLOR[anc], alpha=0.38)
            l_off[anc] += n
            r_off[t]   += n

    # Category legend (right side top)
    cat_patches = [mpatches.Patch(facecolor=gc, alpha=1.0, label=lbl)
                   for lbl, _, _, gc in TRAIT_GROUPS]
    ax.legend(handles=cat_patches, fontsize=6, loc="upper right",
              bbox_to_anchor=(1.26, 1.0), frameon=False)

    # Annotation tag: just below the bottom of the existing plot nodes
    # After invert_yaxis(), total_h is visually at the bottom; we place tag slightly below it.
    # Use DejaVu Sans for this text element to ensure the arrow glyph renders correctly.


    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    print(f"Saved: {path}")
    if extra_paths:
        for ep in extra_paths:
            fig.savefig(ep, bbox_inches="tight", dpi=300, facecolor="white")
            print(f"Saved: {ep}")
    print(f"CAPTION (Fig2A alluvial): GWAS portfolio — {N_GWAS} colocalising GWAS datasets "
          f"({N_STUDIES} distinct studies) spanning {len(ANCESTRIES)} ancestries (European, "
          "African, Admixed American, East Asian, South Asian) and "
          f"{len(TRAITS)} liver-disease traits, with the MVP (Million Veteran Program) "
          "multi-ancestry strata (AFR/AMR/EAS/EUR) included. Left = ancestry (n GWAS), "
          "right = trait; unit = 1 GWAS dataset.")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# C: Radial coverage portrait
# ─────────────────────────────────────────────────────────────────────────────
def build_C(path):
    """Radial: outer ring = each trait arc split into 4 ancestry segments (filled/absent).
       Inner ring = trait category colors. Center = study count."""
    cd = _cell_data()

    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    ax.set_aspect("equal")
    ax.set_xlim(-1.55, 1.55)
    ax.set_ylim(-1.55, 1.55)
    ax.axis("off")
    fig.patch.set_facecolor("white")

    R_OUT     = 1.05   # outer ancestry ring outer radius
    R_IN      = 0.62   # outer ancestry ring inner radius
    R_CAT_OUT = 0.55   # inner category ring outer radius
    R_CAT_IN  = 0.30   # inner category ring inner radius
    R_CENTER  = R_CAT_IN - 0.01

    N_TRAITS      = len(TRAITS)
    DEG_PER_TRAIT = 360.0 / N_TRAITS
    TRAIT_GAP     = 2.0   # degrees between traits
    ANC_GAP       = 0.5   # degrees between ancestry segments within a trait

    # Inner category ring
    cat_angle = {"Liver enzymes": (0, 135), "Disease endpoints": (135, 315), "Imaging": (315, 360)}
    for lbl, c0, c1, gc in TRAIT_GROUPS:
        t1, t2 = c0 * DEG_PER_TRAIT, c1 * DEG_PER_TRAIT
        ax.add_patch(mpatches.Wedge((0, 0), R_CAT_OUT, t1, t2,
                                    width=R_CAT_OUT - R_CAT_IN,
                                    facecolor=gc, edgecolor="white",
                                    linewidth=0.8, alpha=0.80, zorder=1))
        # Category label
        mid_rad = np.radians((t1 + t2) / 2)
        r_mid = (R_CAT_OUT + R_CAT_IN) / 2
        lx, ly = r_mid * np.cos(mid_rad), r_mid * np.sin(mid_rad)
        short = lbl.replace("Disease endpoints", "Disease\nendpoints").replace("Liver enzymes", "Liver\nenzymes")
        ax.text(lx, ly, short, ha="center", va="center",
                fontsize=5.5, color="#444444", style="italic", zorder=5,
                multialignment="center")

    # Outer ancestry coverage ring
    for i, trait in enumerate(TRAITS):
        t_start = i * DEG_PER_TRAIT + TRAIT_GAP / 2
        t_end   = (i + 1) * DEG_PER_TRAIT - TRAIT_GAP / 2
        eff     = t_end - t_start
        n_anc   = len(ANCESTRIES)
        anc_span = (eff - (n_anc - 1) * ANC_GAP) / n_anc

        for k, anc in enumerate(ANCESTRIES):
            a1 = t_start + k * (anc_span + ANC_GAP)
            a2 = a1 + anc_span
            covered = cd[anc][trait] > 0
            color = ANCESTRY_COLOR[anc] if covered else ABSENT_COLOR
            alpha = 0.85 if covered else 0.90
            ax.add_patch(mpatches.Wedge((0, 0), R_OUT, a1, a2,
                                        width=R_OUT - R_IN,
                                        facecolor=color, edgecolor="white",
                                        linewidth=0.7, alpha=alpha, zorder=2))

        # Trait label outside ring
        mid_deg = i * DEG_PER_TRAIT + DEG_PER_TRAIT / 2
        mid_rad = np.radians(mid_deg)
        lx = (R_OUT + 0.17) * np.cos(mid_rad)
        ly = (R_OUT + 0.17) * np.sin(mid_rad)
        # Text alignment based on position
        ha = "center"
        if np.cos(mid_rad) > 0.4:
            ha = "left"
        elif np.cos(mid_rad) < -0.4:
            ha = "right"
        va = "center"
        if np.sin(mid_rad) > 0.4:
            va = "bottom"
        elif np.sin(mid_rad) < -0.4:
            va = "top"
        ax.text(lx, ly, trait, ha=ha, va=va,
                fontsize=8.5, color="#222222", fontweight="bold", zorder=5)

    # White center circle
    ax.add_patch(mpatches.Circle((0, 0), R_CENTER,
                                  facecolor="white", edgecolor="#CCCCCC",
                                  linewidth=0.6, zorder=4))
    ax.text(0, 0.02, str(N_GWAS), ha="center", va="center",
            fontsize=14, color="black", zorder=6)
    ax.text(0, -0.13, "GWAS", ha="center", va="center",
            fontsize=7, color="black", zorder=6)

    # Ancestry legend
    anc_leg = [mpatches.Patch(facecolor=ANCESTRY_COLOR[a], alpha=0.85,
                               label=f"{ANCESTRY_LABEL[a]} ({ANCESTRY_N[a]})")
               for a in ANCESTRIES]
    absent_leg = mpatches.Patch(facecolor=ABSENT_COLOR, label="Not covered")
    ax.legend(handles=anc_leg + [absent_leg], fontsize=7.5,
              loc="lower center", bbox_to_anchor=(0.5, -0.07),
              ncol=3, frameon=False, handlelength=1.2)

    ax.set_title("Option C — Radial coverage  |  outer ring: ancestry × trait  |  filled = covered",
                 fontsize=8, pad=8, loc="center", color="#555555")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Saved: {path}")


# ─────────────────────────────────────────────────────────────────────────────
# D: Unit chart  (stacked squares)
# ─────────────────────────────────────────────────────────────────────────────
def build_D(path):
    """Unit chart: each small square = 1 GWAS study, stacked bottom-up per cell."""
    cd = _cell_data()
    n_anc, n_trait = len(ANCESTRIES), len(TRAITS)

    SQ  = 0.115   # square side in data units
    GAP = 0.018   # gap between stacked squares

    fig, ax = plt.subplots(figsize=(8.5, 3.8))
    fig.patch.set_facecolor("white")
    _group_bands(ax, alpha=0.05)

    for i, anc in enumerate(ANCESTRIES):
        for j, trait in enumerate(TRAITS):
            n = cd[anc][trait]
            if n == 0:
                # Absent placeholder
                ax.add_patch(mpatches.FancyBboxPatch(
                    (j - 0.43, i - 0.43), 0.86, 0.86,
                    boxstyle="round,pad=0.04",
                    facecolor=ABSENT_COLOR, edgecolor="white",
                    linewidth=0.8, alpha=0.6, zorder=1,
                ))
                continue
            cell_bot = i - 0.44
            for k in range(n):
                sq_x = j - SQ / 2
                sq_y = cell_bot + k * (SQ + GAP)
                ax.add_patch(mpatches.FancyBboxPatch(
                    (sq_x, sq_y), SQ, SQ,
                    boxstyle="round,pad=0.012",
                    facecolor=ANCESTRY_COLOR[anc],
                    edgecolor="white", linewidth=0.6,
                    alpha=0.82, zorder=2,
                ))

    for x in np.arange(-0.5, n_trait, 1):
        ax.axvline(x, color="white", linewidth=1.8, zorder=0)
    for y in np.arange(-0.5, n_anc, 1):
        ax.axhline(y, color="white", linewidth=1.8, zorder=0)

    _grid_axis(ax)
    _anc_strips(ax)
    _group_top_labels(ax)

    # Legend
    leg = [mpatches.Patch(facecolor=ANCESTRY_COLOR[a], alpha=0.82,
                          label=f"{ANCESTRY_LABEL[a]} ({ANCESTRY_N[a]})")
           for a in ANCESTRIES]
    leg.append(mpatches.Patch(facecolor=ABSENT_COLOR, label="Not covered"))
    ax.legend(handles=leg, fontsize=7, loc="lower right",
              bbox_to_anchor=(1.0, -0.28), ncol=3, frameon=False)

    ax.set_title("Option D — Unit chart  |  each square = 1 GWAS study  |  color = ancestry",
                 fontsize=8, pad=16, loc="left", color="#555555")
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    plt.close(fig)
    print(f"Saved: {path}")


# ─────────────────────────────────────────────────────────────────────────────
# Extended cascade: Ancestry -> Trait -> Fine-mapped loci -> Coloc genes
# ─────────────────────────────────────────────────────────────────────────────
def build_cascade(path, extra_paths=None):
    """4-stage cascade: Ancestry -> Trait (GWAS datasets) -> fine-mapped loci -> coloc genes.
    Trait = spine for the 3 right stages; colour = ancestry; each stage normalised to equal
    height (totals differ) with the total labelled. FULLY DATA-DRIVEN: GWAS from the registry,
    per-(ancestry,trait) loci from study_summary n_loci, genes = union SuSiE-or-ABF PP.H4>0.5
    with each gene assigned to the ancestry/trait of its single best colocalisation."""
    r_order = list(TRAITS)
    at_g, at_l, at_ge = _M["at_gwas"], _M["at_loci"], _M["at_genes"]
    DAT = {}
    for a in ANCESTRIES:
        for t in r_order:
            g, l, ge = at_g.get((a, t), 0), at_l.get((a, t), 0), at_ge.get((a, t), 0)
            if g or l or ge:
                DAT[(a, t)] = (g, l, ge)
    MIDX = {"trait": 0, "loci": 1, "genes": 2}
    H, GAP, BW = 13.0, 0.42, 0.06
    # equal horizontal width for every section (ancestry@0 -> trait -> loci -> genes)
    SEC = 0.82
    XCOL = {"trait": SEC, "loci": 2*SEC, "genes": 3*SEC}

    def layout(metric):
        tot = sum(v[metric] for v in DAT.values())
        scale = (H - GAP*(len(r_order)-1)) / tot if tot else 0.0
        seg, tspan = {}, {}; y = 0.0
        for t in r_order:
            y0t = y
            for a in ANCESTRIES:
                v = DAT.get((a, t), (0, 0, 0))[metric]
                if v > 0:
                    seg[(a, t)] = (y, y + v*scale); y += v*scale
            tspan[t] = (y0t, y); y += GAP
        return seg, tspan, scale
    L = {s: layout(MIDX[s]) for s in XCOL}
    tot = {s: sum(v[MIDX[s]] for v in DAT.values()) for s in XCOL}

    fig, ax = plt.subplots(figsize=(5.6, 4.3))
    ax.set_xlim(-0.62, 3*SEC + 0.5); ax.set_ylim(-1.35, H + 1.1); ax.invert_yaxis(); ax.axis("off")

    # Ancestry column (x=0)
    anc_tot = {a: sum(DAT[k][0] for k in DAT if k[0] == a) for a in ANCESTRIES}
    s0 = (H - GAP*(len(ANCESTRIES)-1)) / sum(anc_tot.values())
    anc_y = {}; y = 0.0
    for a in ANCESTRIES:
        h = anc_tot[a]*s0; anc_y[a] = (y, y+h)
        ax.add_patch(mpatches.Rectangle((-BW/2, y), BW, h, facecolor=ANCESTRY_COLOR[a], edgecolor="none", zorder=4))
        ax.text(-BW/2-0.03, y+h/2, f"{ANCESTRY_LABEL[a]} ({anc_tot[a]})", ha="right", va="center", fontsize=6, color="black")
        y += h + GAP
    ax.text(0, -0.7, f"{len(ANCESTRIES)} ancestries", ha="center", va="bottom", fontsize=6, color="black")

    # stage blocks + headers
    for s in XCOL:
        seg = L[s][0]; x = XCOL[s]
        for (a, t), (y0, y1) in seg.items():
            ax.add_patch(mpatches.Rectangle((x-BW/2, y0), BW, y1-y0, facecolor=ANCESTRY_COLOR[a], edgecolor="white", lw=0.2, zorder=4))
        hdr = {"trait": f"{tot['trait']} GWAS datasets",
               "loci": f"{tot['loci']} fine-mapped loci",
               "genes": f"{tot['genes']} coloc genes"}[s]
        ax.text(x, -0.7, hdr, ha="center", va="bottom", fontsize=6, color="black")
    # trait name + dataset count ABOVE each trait block, e.g. "NAFLD (11)"
    for t in r_order:
        y0, y1 = L["trait"][1][t]
        n_gw = sum(DAT[(a, t)][MIDX["trait"]] for a in ANCESTRIES if (a, t) in DAT)
        if n_gw == 0:
            continue
        ax.text(XCOL["trait"], y0 - 0.06, f"{TRAIT_SHORT.get(t, t)} ({n_gw})",
                ha="center", va="bottom", fontsize=6, color="black", zorder=8)
    # per-trait count ABOVE the loci / genes blocks (no name to pair with there)
    for s in ("loci", "genes"):
        tspan = L[s][1]; x = XCOL[s]
        for t in r_order:
            y0, y1 = tspan[t]
            cnt = sum(DAT[(a, t)][MIDX[s]] for a in ANCESTRIES if (a, t) in DAT)
            if cnt == 0:
                continue
            ax.text(x, y0 - 0.06, str(cnt), ha="center", va="bottom", fontsize=6, color="black", zorder=8)

    # ribbons: Anc -> Trait
    seg_tr = L["trait"][0]; off = {a: anc_y[a][0] for a in ANCESTRIES}
    for a in ANCESTRIES:
        for t in r_order:
            if (a, t) in seg_tr:
                n = DAT[(a, t)][0]*s0; hy0 = off[a]; off[a] += n
                ry0, ry1 = seg_tr[(a, t)]
                _ribbon(ax, BW/2, hy0, hy0+n, XCOL["trait"]-BW/2, ry0, ry1, ANCESTRY_COLOR[a], alpha=0.30)
    # ribbons: Trait->Loci, Loci->Genes
    for sf, sg in [("trait", "loci"), ("loci", "genes")]:
        segf, segg = L[sf][0], L[sg][0]
        for a in ANCESTRIES:
            for t in r_order:
                if (a, t) in segf and (a, t) in segg:
                    f0, f1 = segf[(a, t)]; g0, g1 = segg[(a, t)]
                    _ribbon(ax, XCOL[sf]+BW/2, f0, f1, XCOL[sg]-BW/2, g0, g1, ANCESTRY_COLOR[a], alpha=0.30)
    # fallback: coloc genes whose GWAS has no fine-mapped-loci entry (the 7 MVP-EUR strata,
    # fine-mapped under a separate PolyFun-LD pipeline) -> faint ribbon straight from the
    # trait block across the empty loci column to the gene block, so no gene block floats.
    seg_tr2, seg_lo, seg_ge = L["trait"][0], L["loci"][0], L["genes"][0]
    for a in ANCESTRIES:
        for t in r_order:
            if (a, t) in seg_ge and (a, t) not in seg_lo and (a, t) in seg_tr2:
                g0, g1 = seg_ge[(a, t)]; t0, t1 = seg_tr2[(a, t)]
                _ribbon(ax, XCOL["trait"]+BW/2, t0, t1, XCOL["genes"]-BW/2, g0, g1, ANCESTRY_COLOR[a], alpha=0.12)

    fig.tight_layout(); fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    print(f"Saved: {path}")
    if extra_paths:
        for ep in extra_paths:
            fig.savefig(ep, bbox_inches="tight", dpi=300, facecolor="white"); print(f"Saved: {ep}")
    print(f"CAPTION (Fig2A cascade): GWAS evidence cascade — {N_GWAS} colocalising GWAS datasets "
          f"({N_STUDIES} distinct studies; {len(ANCESTRIES)} ancestries x {len(TRAITS)} liver traits, "
          f"MVP included) -> {tot['loci']} fine-mapped loci -> {tot['genes']} colocalising genes "
          "(union SuSiE-or-ABF PP.H4>0.5), coloured by ancestry. Non-European panels (AFR/AMR/EAS) yield "
          "the most fine-mapped loci but far fewer colocalising genes, whereas European GWAS yield few "
          "loci yet most coloc genes — the European-eQTL ancestry bottleneck. Fine-mapped-loci counts are "
          "per-study locus windows (LD-panel dependent, so non-European panels report more windows); EUR "
          "loci exclude the 7 MVP-EUR strata fine-mapped under a separate PolyFun-LD pipeline (their "
          "colocalising genes are retained, shown as faint trait->gene ribbons).")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
# Companion to the cascade: per-study fine-mapped loci & coloc genes
# ─────────────────────────────────────────────────────────────────────────────
def build_study_bars(path, extra_paths=None):
    """Two aligned horizontal bar panels (loci | genes) sharing one row per GWAS study,
    grouped by ancestry. DATA-DRIVEN: rows = the studies carrying a fine-mapping summary
    (study_summary.csv); loci = n_loci; genes = union SuSiE-or-ABF PP.H4>0.5 assigned to
    that study; the Sveinbjornsson 2023 cohort-arms are collapsed to one study. The
    per-study loci-vs-genes contrast shows the European-eQTL bottleneck directly."""
    # Pretty display names for the PMID-keyed legacy EUR studies (curated); every other
    # study gets an auto-generated "<source> <trait>" label.
    NAME_MAP = {
        "2019_31311600": "Namjou 2019 NAFLD",   "2020_32298765": "Anstee 2020 NAFLD",
        "2021_34128465": "Liu 2021 PDFF",         "2021_34841290": "Ghodsian 2021 NAFLD",
        "2021_34957434": "Haas 2021 PDFF",         "2022_36402844": "van der Meer 2022 PDFF",
        "2023_36280732": "Sveinbjornsson 2022 NAFLD",
    }

    def _label(study):
        for pref, lab in NAME_MAP.items():
            if study.startswith(pref):
                return lab
        t = _trait_of(study)
        for pref, src in (("MVP_", "MVP"), ("BBJ_", "BBJ"), ("UKBB_", "UKBB"),
                          ("FinnGen_", "FinnGen"), ("PanUKBB_", "PanUKBB")):
            if study.startswith(pref):
                return f"{src} {t}"
        return study.replace("_", " ")

    # Aggregate per study (collapsing Sveinbjornsson arms): loci from study_summary,
    # genes from the union coloc assignment.
    agg = {}
    for s, loci in _M["study_loci"].items():
        a = _M["anc_of_study"][s]
        key = "2023_36280732_NAFLD" if s.startswith("2023_36280732") else s
        genes = _M["study_genes"].get(s, 0)
        if key not in agg:
            agg[key] = [a, 0, 0, _label(s)]
        agg[key][1] += loci
        agg[key][2] += genes
    # order rows by ancestry (EUR, AFR, EAS, AMR, SAS) then genes descending
    ROWS = []
    for a in ANCESTRIES:
        rws = sorted([(v[3], a, v[1], v[2]) for v in agg.values() if v[0] == a],
                     key=lambda z: -z[3])
        ROWS.extend(rws)

    n = len(ROWS); ys = list(range(n-1, -1, -1))  # first row at top
    loci_mx  = max(r[2] for r in ROWS) * 1.35
    genes_mx = max(r[3] for r in ROWS) * 1.12
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(7.8, 9.6), sharey=True,
                                   gridspec_kw={"wspace": 0.5})
    # Fine-mapped loci span ~3 orders of magnitude (MVP AFR/AMR panels are LD-inflated),
    # so the loci axis is log-x; coloc genes stay linear.
    for ax, idx, ttl, mx, logx in [(axL, 2, "Fine-mapped loci (log)", loci_mx, True),
                                   (axR, 3, "Coloc genes (union PP.H4 > 0.5)", genes_mx, False)]:
        for i, row in enumerate(ROWS):
            v = row[idx]
            ax.barh(ys[i], v, color=ANCESTRY_COLOR[row[1]], height=0.66, zorder=3)
            if v > 0:
                xoff = v * 1.08 if logx else v + mx*0.013
                ax.text(xoff, ys[i], str(v), va="center", ha="left", fontsize=6, color="black")
        if logx:
            ax.set_xscale("log"); ax.set_xlim(0.8, mx)
        else:
            ax.set_xlim(0, mx)
        ax.set_ylim(-0.7, n-0.3)
        ax.set_title(ttl, fontsize=7, color="black")
        for sp in ("top", "right", "left"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(axis="x", labelsize=6, colors="black")
        ax.set_yticks([])
    # study labels + ancestry-group separators on the left panel
    axL.set_yticks(ys); axL.set_yticklabels([r[0] for r in ROWS], fontsize=6, color="black")
    bounds = []; prev = ROWS[0][1]
    for i, row in enumerate(ROWS):
        if row[1] != prev:
            bounds.append(i - 0.5); prev = row[1]
    for b in bounds:
        for ax in (axL, axR):
            ax.axhline(n-1-b, color="#CCCCCC", lw=0.5, zorder=1)

    fig.savefig(path, bbox_inches="tight", dpi=300, facecolor="white")
    print(f"Saved: {path}")
    if extra_paths:
        for ep in extra_paths:
            fig.savefig(ep, bbox_inches="tight", dpi=300, facecolor="white"); print(f"Saved: {ep}")
    print(f"CAPTION (Fig2A per-study): Fine-mapped loci (log scale) and colocalising genes "
          f"(union SuSiE-or-ABF PP.H4>0.5) per GWAS study, for the {n} studies carrying a "
          f"fine-mapping summary, grouped by {len(ANCESTRIES)} ancestries (incl. MVP AFR/AMR/EAS). "
          "The loci-vs-genes contrast exposes the European-eQTL ancestry bottleneck at the study "
          "level: non-European studies (MVP/BBJ/PanUKBB AFR/AMR/EAS) yield many fine-mapped loci "
          "but few coloc genes, whereas European UKBB studies yield few loci but many genes. The 7 "
          "MVP-EUR strata (fine-mapped under a separate PolyFun-LD pipeline) are not shown here but "
          f"contribute to the {UNION_GENES} coloc-gene total.")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    build_A(_out_path("gwas_creative_A_binary_grid.pdf"))
    PANELS = os.path.join(_project_root(), "figures", "main", "fig2_genetics", "panels")
    build_B(_out_path("gwas_creative_B_alluvial.pdf"),
            extra_paths=[
                _out_path("fig1c_gwas_alluvial.pdf"),
                os.path.join(PANELS, "fig2A_gwas_alluvial.pdf"),
            ])
    build_cascade(_out_path("gwas_creative_E_cascade.pdf"),
                  extra_paths=[os.path.join(PANELS, "fig2A_gwas_cascade.pdf")])
    build_study_bars(_out_path("gwas_creative_F_perstudy.pdf"),
                     extra_paths=[os.path.join(PANELS, "FigS2A_gwas_perstudy.pdf")])
    print("All done.")
