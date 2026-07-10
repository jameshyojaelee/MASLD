#!/usr/bin/env python3
"""
Fig 4A validation-overview redesign — SHARED data + style layer for the 5 v2
candidates (see docs/plan: figures-main-fig4-validation-panels-fig4-concurrent).

Why this module exists
----------------------
The PI rejected the current Sankey (`fig4a_overview_cascade.pdf`) because it
DISPLAYS the whole evidence base without ASSERTING a claim. The 5 v2 candidates
each weld an honest claim onto the geometry. They MUST all speak the same data +
style so they are interchangeable at the data layer and obey house style.

This module:
  * reuses the LIVE data plumbing of the canonical generator — imports and runs
    `compute_counts()` + `permutation_null()` from `fig4a_overview_candidates.py`
    (no number is hard-coded; a broken source RAISES rather than mislabels), and
  * adds the exclusive convergence partition derived FROM the live gene sets
    (P-only/S-only/A-only, exactly-two PS/PA/SA, all-3) so every area/edge/mark is
    mutually consistent by construction (avoids the 141-vs-135 inclusion-exclusion
    trap that bites if you sum reported pairwise overlaps), and
  * exposes the frozen house style (PDF, Helvetica 6, near-black ink #231f20, no
    bold, NO gene names, colour on marks only from the canonical modality palette,
    control/empty grays, ASCII-only on-panel) + a couple of shared helpers.

Honest statistics (the claim every candidate asserts, from the seed-42 / 10k null):
  the ONLY convergence beating chance is proteo∩spatial (obs vs exp in
  stats["conv"]["PS"]); ≥2 overall is at chance; genetic-only enriches in nothing.
Do NOT cite the retired "3.76x on a 13,390 universe" figure — it does not survive
on the current 9,882 universe.

Env: rnaseq python (pandas/numpy). matplotlib import can hang on the login node ->
render via SLURM (run_fig4a_v2.sh). No scipy (rnaseq ABI risk) — the null is pure numpy.
"""
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import (Rectangle, PathPatch, FancyBboxPatch, Circle,
                                 Wedge, Ellipse, FancyArrowPatch, Polygon)
from matplotlib.path import Path

# -- Paths --------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# Reuse the CANONICAL live data plumbing (single source of truth for every count)
# and the modality palette. compute_counts/permutation_null read the atlas, DEG,
# proteomics and universe files live and RAISE on any unreadable source.
from fig4a_overview_candidates import (compute_counts, permutation_null,  # noqa: E402
                                       _ribbon)
from fig1_palette import MODALITY_COLORS, CONTROL_GRAY  # noqa: E402

# -- Canonical palette (marks only; text is always INK) -----------------------
C_BULK    = MODALITY_COLORS["bulk"]      # #0072B2  transcriptomics origin
C_GWAS    = MODALITY_COLORS["gwas"]      # #C2185B  genetics origin
C_PROTEO  = MODALITY_COLORS["proteo"]    # #D55E00  proteomics
C_SPATIAL = MODALITY_COLORS["spatial"]   # #CC79A7  spatial
C_ATAC    = MODALITY_COLORS["atac"]      # #E69F00  snATAC
GRAY      = CONTROL_GRAY                  # #9E9E9E  control/absent (locked invariant)
C_NS      = "#E4E4E4"                     # empty / unvalidated field (lighter neutral)
C_CONV    = "#2E6E8E"                     # convergent-core accent (deep teal, from parent)
INK       = "#231f20"                     # near-black ink (2026-07-09 house re-skin)
FS        = 6                             # SINGLE font size for every text element

# modality -> colour, for builders that iterate lenses
LENS_COL   = {"P": C_PROTEO, "S": C_SPATIAL, "A": C_ATAC}
LENS_NAME  = {"P": "proteomics", "S": "spatial", "A": "snATAC"}

plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":       FS,
    "pdf.fonttype":    42,   # editable text in Illustrator (Type42, NOT Type3)
    "ps.fonttype":     42,
    "figure.dpi":      150,
    "axes.linewidth":  0.5,
})

# v2 candidates render to a dedicated subdir so they never touch the live panel.
OUT_DIR = os.path.join(BASE, "figures/main/fig4_validation", "panels", "candidates")


# -- Data ---------------------------------------------------------------------
def get_data():
    """Run the canonical live compute + null once, and attach the EXCLUSIVE
    convergence partition derived from the live gene sets (so every mark a builder
    draws is mutually consistent). Returns (d, stats).

      d               : the full compute_counts() dict (scalar counts + res["_sets"])
      d["_part"]      : exclusive partition (see below), all disjoint, summing to ge1
      stats["conv"]   : per-EXCLUSIVE-intersection null — keys ge2/all3/PS/PA/SA/
                        P_only/S_only/A_only, each {obs,exp,fold,p,lo,hi}
      stats["enrich"] : {"universe"|"coloc_only": {lens: {obs,exp,fold,p,N,K,k}}}
                        (per-lens hypergeometric enrichment vs genome-wide background)
    """
    d = compute_counts()
    stats = permutation_null(d)

    S = d["_sets"]
    def _clean(s):
        return {x for x in s if isinstance(x, str) and x and x.lower() != "nan"}
    P, Sp, A = _clean(S["Pval"]), _clean(S["Sval"]), _clean(S["Aval"])

    part = dict(
        P_only=len(P - Sp - A),
        S_only=len(Sp - P - A),
        A_only=len(A - P - Sp),
        PS=len((P & Sp) - A),          # exactly-two, proteo & spatial (not snATAC)
        PA=len((P & A) - Sp),
        SA=len((Sp & A) - P),
        all3=len(P & Sp & A),
    )
    part["ge1"] = len(P | Sp | A)                                   # any-assay validated
    part["ge2"] = part["PS"] + part["PA"] + part["SA"] + part["all3"]  # exactly->=2 (clean)
    part["unvalidated"] = d["universe"] - part["ge1"]
    # origin split of the prioritized universe (for opener strips / marimekko)
    rna, coloc = _clean(S["rna"]), _clean(S["coloc"])
    part["tx_only"]      = len(rna - coloc)
    part["overlap"]      = len(rna & coloc)
    part["genetic_only"] = len(coloc - rna)
    d["_part"] = part

    banner(d, part)
    return d, stats


def banner(d, part):
    """Print the live counts prominently so a data drift is loud in the SLURM log
    (we deliberately do NOT hard-assert frozen numbers — the universe can change
    legitimately, and the panels read every count live so labels track the data)."""
    print("\n" + "=" * 72)
    print("[fig4a_v2] LIVE COUNTS (labels/marks read these — no frozen locks):")
    print(f"  universe {d['universe']:,}  = transcriptomics {d['substrate']:,} "
          f"(tx-only {part['tx_only']:,}) + genetics {d['coloc']:,} "
          f"(genetic-only {part['genetic_only']:,}); overlap {part['overlap']:,}")
    print(f"  validated: proteomics {d['P_val']} (liver {d['P_liver_val']} + "
          f"plasma {d['P_plasma_val']}) | spatial {d['S_val']} "
          f"(Visium {d['S_hep_val']} + GeoMx {d['S_zon_val']}) | snATAC {d['A_val']}")
    print(f"  convergence (exclusive): >=1 {part['ge1']:,} | >=2 {part['ge2']} | "
          f"all-3 {part['all3']} | proteo&spatial-only {part['PS']} | "
          f"proteo&snATAC-only {part['PA']} | spatial&snATAC-only {part['SA']}")
    print(f"  unvalidated {part['unvalidated']:,}")
    print("=" * 72 + "\n")


# -- Shared drawing helpers ---------------------------------------------------
def save_pdf(fig, name):
    """Save one candidate PDF into panels/candidates/ (never the live panel dir)."""
    os.makedirs(OUT_DIR, exist_ok=True)
    if not name.endswith(".pdf"):
        name += ".pdf"
    p = os.path.join(OUT_DIR, name)
    fig.savefig(p, bbox_inches="tight", pad_inches=0.02, facecolor="white")
    plt.close(fig)
    print("Saved:", p)
    return p


def tint(hexc, f):
    """Lighten a hex colour toward white by fraction f in [0,1] (sub-node tints)."""
    import matplotlib.colors as mcolors
    r, g, b = mcolors.to_rgb(hexc)
    return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)


def fold_ci(st):
    """(fold, lo_fold, hi_fold) from a conv[...] entry's obs/exp and 95% band,
    guarding exp==0. Used by the forest/dumbbell CIs."""
    exp = st.get("exp", float("nan"))
    if not np.isfinite(exp) or exp <= 0:
        return float("nan"), float("nan"), float("nan")
    obs = st["obs"]
    lo, hi = st.get("lo", exp), st.get("hi", exp)
    # fold of observed vs expected; CI expressed as the null band relative to exp
    return obs / exp, (lo / exp if exp else float("nan")), (hi / exp if exp else float("nan"))


if __name__ == "__main__":
    # smoke test: compute + print, write nothing
    _d, _s = get_data()
    print("[fig4a_v2_common] OK — conv keys:", sorted(_s["conv"]) if _s else None)
