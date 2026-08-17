#!/usr/bin/env python3
"""
Fig 4A validation-overview redesign — SHARED data + style layer for the 5 v2
candidates (see docs/plan: figures-main-fig4-validation-panels-fig4-concurrent).

Why this module exists
----------------------
The PI rejected the current Sankey (`fig5a_overview_cascade_candidate.pdf`) because it
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
OUT_DIR = os.path.join(BASE, "figures/main/fig5_molecular_context", "panels", "candidates")


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
    # DEFAULT = STRICT (significance-test-only) definition. The inclusive definition
    # (curated GeoMx zonation signature + static promoter-accessible / chromVAR-SCENIC
    # TF flags counted as "validated") inflates spatial 112->447 and snATAC 233->560 and
    # makes measured ~= validated for those assays. Strict keeps only per-assay disease
    # SIGNIFICANCE TESTS; under it proteome-and-spatial is stronger (2.07x, p=9e-4 vs
    # 1.36x, p=2.6e-3). Set FIG4A_INCLUSIVE=1 for the original inclusive definition.
    if os.environ.get("FIG4A_INCLUSIVE"):
        return _get_data_inclusive()
    return _compute_strict()


def _partition(d):
    """Attach the exclusive convergence partition (all disjoint, summing to ge1) +
    the origin split, derived from the live gene sets. Shared by both definitions."""
    S = d["_sets"]
    def _clean(s):
        return {x for x in s if isinstance(x, str) and x and x.lower() != "nan"}
    P, Sp, A = _clean(S["Pval"]), _clean(S["Sval"]), _clean(S["Aval"])
    part = dict(
        P_only=len(P - Sp - A), S_only=len(Sp - P - A), A_only=len(A - P - Sp),
        PS=len((P & Sp) - A),          # exactly-two, proteo & spatial (not snATAC)
        PA=len((P & A) - Sp), SA=len((Sp & A) - P), all3=len(P & Sp & A),
    )
    part["ge1"] = len(P | Sp | A)
    part["ge2"] = part["PS"] + part["PA"] + part["SA"] + part["all3"]
    part["unvalidated"] = d["universe"] - part["ge1"]
    rna, coloc = _clean(S["rna"]), _clean(S["coloc"])
    part["tx_only"], part["overlap"], part["genetic_only"] = \
        len(rna - coloc), len(rna & coloc), len(coloc - rna)
    d["_part"] = part
    return part


def _get_data_inclusive():
    """INCLUSIVE definition (original): validated = significant in ANY of a modality's
    layers, incl. curated signatures + static flags. Kept for comparison (FIG4A_INCLUSIVE=1)."""
    d = compute_counts()
    stats = permutation_null(d)
    part = _partition(d)
    banner(d, part)
    return d, stats


# -- STRICT definition (significance-test-only) -------------------------------
_UDIR  = os.path.join(BASE, "Analysis/Spatial/results/universe_validation")
_ATLAS = os.path.join(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
_PROT  = os.path.join(BASE, "Analysis/Proteomics/results/protein_transcript_concordance_v3.csv")
_DEG   = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/"
                            "results/integration/canonical_deg_results.csv")


def _sclean(s):
    return {x for x in s if isinstance(x, str) and x and x.lower() != "nan"}


def _compute_strict():
    """Validated = per-assay disease SIGNIFICANCE TEST only (no signatures / static flags):
    proteomics = protein DE padj<0.05 (measured = DIA-MS detected); spatial = Visium-hep
    Wilcoxon DE (measured = tested); snATAC = hepatocyte DA (measured = tested). Same
    schema as _get_data_inclusive so every builder reads strict numbers unchanged."""
    import pandas as pd

    def _rd(f):
        return _sclean(set(open(os.path.join(_UDIR, f)).read().split()))
    rna, coloc = _rd("universe_transcriptomic.txt"), _rd("universe_genetic.txt")
    target = _rd("prioritized_universe_FINAL.txt")
    coloc_only = coloc - rna

    deg = pd.read_csv(_DEG)
    deg["sym"] = deg["symbol"].fillna(deg["gene"])
    treat = _sclean(set(deg.loc[deg["treat_fdr"] < 0.05, "sym"]))

    c = pd.read_csv(_PROT)
    P_meas_all = _sclean(set(c["gene"].dropna()))
    P_val_all = _sclean(set(c.loc[c["protein_padj"] < 0.05, "gene"].dropna()))

    a = pd.read_csv(_ATLAS, low_memory=False).rename(columns={"human_symbol": "sym"})
    def _nn(col):   return _sclean(set(a.loc[a[col].notna(), "sym"])) if col in a.columns else set()
    def _padj(col): return _sclean(set(a.loc[a[col] < 0.05, "sym"])) if col in a.columns else set()
    S_meas_all, S_val_all = _nn("spatial_hep_wilcoxon_padj_bh"), _padj("spatial_hep_wilcoxon_padj_bh")
    A_meas_all, A_val_all = _nn("hepatocyte_da_padj"), _padj("hepatocyte_da_padj")
    P_val_all &= P_meas_all; S_val_all &= S_meas_all; A_val_all &= A_meas_all

    Pmeas, Pval = P_meas_all & target, P_val_all & target
    Smeas, Sval = S_meas_all & target, S_val_all & target
    Ameas, Aval = A_meas_all & target, A_val_all & target

    def _pds(prefix):
        x = c[c["dataset"].str.startswith(prefix)]
        m = _sclean(set(x["gene"].dropna())) & target
        v = _sclean(set(x.loc[x["protein_padj"] < 0.05, "gene"].dropna())) & target
        return len(m), len(v)
    P_liver_meas, P_liver_val = _pds("PXD051911")
    P_plasma_meas, P_plasma_val = _pds("PXD052937")

    P, Sp, A = Pval, Sval, Aval
    d = dict(
        substrate=len(rna), coloc=len(coloc), coloc_shared=len(rna & coloc),
        coloc_only=len(coloc_only), universe=len(target), confident=len(rna & treat),
        P_meas=len(Pmeas), P_val=len(Pval), S_meas=len(Smeas), S_val=len(Sval),
        A_meas=len(Ameas), A_val=len(Aval),
        P_liver_meas=P_liver_meas, P_liver_val=P_liver_val,
        P_plasma_meas=P_plasma_meas, P_plasma_val=P_plasma_val,
        S_hep_meas=len(Smeas), S_hep_val=len(Sval), S_zon_meas=0, S_zon_val=0,
        S_cos_meas=0, S_cos_val=0,
        A_acc_meas=len(Ameas), A_acc_val=len(Aval), A_tf_meas=0, A_tf_val=0,
        PS=len(P & Sp), PA=len(P & A), SA=len(Sp & A),
        P_only=len(P - Sp - A), S_only=len(Sp - P - A), A_only=len(A - P - Sp),
        all3=len(P & Sp & A), ge1=len(P | Sp | A),
    )
    d["ge2"] = len(P & Sp) + len(P & A) + len(Sp & A) - 2 * d["all3"]  # genes in >=2
    d["none"] = d["universe"] - d["ge1"]
    d.update(proteo=d["P_val"], spatial=d["S_val"], scatac=d["A_val"], PSA=d["all3"])
    d["_sets"] = dict(
        universe=target, rna=rna, coloc=coloc, coloc_only=coloc_only, treat=treat,
        Pmeas=Pmeas, Pval=Pval, Smeas=Smeas, Sval=Sval, Ameas=Ameas, Aval=Aval,
        P_meas_all=P_meas_all, P_val_all=P_val_all, S_meas_all=S_meas_all,
        S_val_all=S_val_all, A_meas_all=A_meas_all, A_val_all=A_val_all)
    part = _partition(d)

    stats = _strict_null(d)
    banner(d, part)
    print("[fig4a_v2] STRICT definition: validated = significance test only "
          "(spatial=Visium-hep DE, snATAC=hepatocyte DA; signatures/flags dropped). "
          "FIG4A_INCLUSIVE=1 reverts.")
    return d, stats


def _strict_null(d, n_perm=10000, seed=42):
    """Same random-placement convergence null + hypergeometric per-lens enrichment as
    the canonical permutation_null, on the STRICT sets. Fresh rng for the convergence
    loop (independent of the enrichment draws) so it reproduces the strict recompute."""
    S = d["_sets"]
    rng_e = np.random.default_rng(seed)

    def hyper(pset, meas_all, val_all):
        Nn, K = len(meas_all), len(val_all)
        k, obs = len(pset & meas_all), len(pset & val_all)
        if k == 0 or K == 0 or K >= Nn:
            return dict(obs=obs, exp=float("nan"), fold=float("nan"), p=1.0, N=Nn, K=K, k=k)
        draws = rng_e.hypergeometric(K, Nn - K, k, size=n_perm)
        exp = float(draws.mean())
        return dict(obs=obs, exp=exp, fold=(obs / exp if exp > 0 else float("inf")),
                    p=(int(np.sum(draws >= obs)) + 1) / (n_perm + 1), N=Nn, K=K, k=k)

    lenses = (("proteomics", S["P_meas_all"], S["P_val_all"]),
              ("spatial",    S["S_meas_all"], S["S_val_all"]),
              ("scATAC",     S["A_meas_all"], S["A_val_all"]))
    enrich = {nm: {ln: hyper(P, m, v) for ln, m, v in lenses}
              for nm, P in (("universe", S["universe"]), ("coloc_only", S["coloc_only"]))}

    rng_c = np.random.default_rng(seed)
    uni = sorted(g for g in S["universe"] if isinstance(g, str))
    idx = {g: i for i, g in enumerate(uni)}
    U = len(uni)
    meas_idx = {L: np.fromiter((idx[g] for g in ms if g in idx), dtype=np.int64)
                for L, ms in (("P", S["Pmeas"]), ("S", S["Smeas"]), ("A", S["Ameas"]))}
    val_n = {"P": len(S["Pval"]), "S": len(S["Sval"]), "A": len(S["Aval"])}
    keys = ("ge2", "all3", "P_only", "S_only", "A_only", "PS", "PA", "SA")
    null = {k: np.empty(n_perm, dtype=np.int32) for k in keys}
    for it in range(n_perm):
        mP = np.zeros(U, bool); mS = np.zeros(U, bool); mA = np.zeros(U, bool)
        for L, m in (("P", mP), ("S", mS), ("A", mA)):
            if val_n[L] and len(meas_idx[L]) >= val_n[L]:
                m[rng_c.choice(meas_idx[L], size=val_n[L], replace=False)] = True
        both = (mP & mS) | (mP & mA) | (mS & mA)
        null["ge2"][it] = int(both.sum())
        null["all3"][it] = int((mP & mS & mA).sum())
        null["P_only"][it] = int((mP & ~mS & ~mA).sum())
        null["S_only"][it] = int((mS & ~mP & ~mA).sum())
        null["A_only"][it] = int((mA & ~mP & ~mS).sum())
        null["PS"][it] = int((mP & mS & ~mA).sum())
        null["PA"][it] = int((mP & mA & ~mS).sum())
        null["SA"][it] = int((mS & mA & ~mP).sum())

    def summ(obs, nul):
        exp = float(nul.mean())
        lo, hi = np.percentile(nul, [2.5, 97.5])
        return dict(obs=int(obs), exp=exp, fold=(obs / exp if exp > 0 else float("inf")),
                    p=(int(np.sum(nul >= obs)) + 1) / (n_perm + 1), lo=float(lo), hi=float(hi))

    obs_excl = dict(ge2=d["ge2"], all3=d["all3"], P_only=d["P_only"],
                    S_only=d["S_only"], A_only=d["A_only"],
                    PS=d["PS"] - d["all3"], PA=d["PA"] - d["all3"], SA=d["SA"] - d["all3"])
    conv = {k: summ(obs_excl[k], null[k]) for k in keys}
    return dict(enrich=enrich, conv=conv, lenses=[l[0] for l in lenses])


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
    # smoke test: compute + print the key convergence verdicts, write nothing
    _d, _s = get_data()
    for _k in ("PS", "PA", "SA", "ge2", "all3"):
        _c = _s["conv"][_k]
        print(f"  conv {_k:>4}: obs={_c['obs']:>4}  exp={_c['exp']:>7.2f}  "
              f"[{_c['lo']:.1f},{_c['hi']:.1f}]  {_c['fold']:.2f}x  p={_c['p']:.2e}")
    print("[fig4a_v2_common] OK")
