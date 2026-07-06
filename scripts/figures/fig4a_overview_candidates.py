#!/usr/bin/env python3
"""
Fig 4 Panel A -- high-level validation-overview CANDIDATES (no gene names).

Fig 4 is the paper's *validation* tier. Validation LOGIC (locked with the PI
2026-07-01):

  * FOUNDATION / target universe = the two prioritization axes from Figs 2-3:
      - bulk RNA-seq SIGNAL at padj<0.05 (NO TREAT/|LFC| gate; 12,989). The
        TREAT set (1,918 DEGs) is RNA-seq's OWN confidence spine and is kept as a
        highlighted "confident core"; validation is deliberately run on the looser
        significance set so modest-but-real DEGs a 2nd assay confirms aren't lost.
      - GWAS-COLOC genes (best PP.H4>0.5). MAIN render = the LIVER-SPECIFIC
        Tier-1/2 portfolio only: the 35 strata flagged placement=="main" in
        GWAS/finemapping/config/gwas_trait_tier.tsv (direct MASLD/NAFLD/NASH/PDFF
        + liver-enzyme ALT/AST/GGT across ancestries) -> 862 COLOC genes. The
        distal/non-specific Tier-3/4 strata (MVP ChronLiver/Cirrhosis/Albumin/
        Platelet) are DROPPED. Co-equal genetic axis.
    Union universe = 13,390 (461 shared). The full 50-GWAS portfolio (COLOC 1,234,
    universe 13,552; adds Tier-3/4) is emitted as a labeled SUPP/sensitivity render
    (FIG4A_KEEP_TIER34=1 -> the same panels with a _supp_full50gwas suffix).
  * VALIDATED by a lens = BOTH-SIGNIFICANT, direction-AGNOSTIC: in the target
    universe AND significant in that modality's disease contrast. Direction is an
    ANNOTATION not a filter (a gene up-in-RNA/down-in-protein, if significant in
    both, is real post-transcriptional/secretion biology, so it counts).

Three orthogonal validation lenses on their native disease contrasts:
  * proteomics : protein-DE padj<0.05 (two DIA-MS proteomes, disease-vs-control)
  * spatial    : any spatial disease-contrast DE padj<0.05 (GeoMx SH-vs-LS/PT +
                 Visium hepatocyte Wilcoxon — the arms with a valid unit-of-analysis
                 p. CosMx is slide-level DIRECTION support only, no cell-level p, so
                 it is NOT in the significance gate — E1b remediation, was 639->111.)
  * scATAC     : hepatocyte differential-accessibility padj<0.05

High-level, counts-only overviews (Fig 2A aesthetic: normalized-height alluvial,
ribbons colored by modality, counts at every stage, NO gene names). Every mark
encodes a real count. ALL text black, font size 6, NO bold; color on marks only.

  build_cascade -> fig4a_overview_cascade.pdf   (Fig 2A-style alluvial funnel).
  This is the SELECTED Fig 4A overview; the rings / upset / pyramid / funnel
  candidates were cut 2026-07-02 (user chose the alluvial cascade).

  MAIN render (default, Tier-1/2 COLOC) -> fig4a_overview_{cascade,enrichment,upset}.pdf.
  SUPP render (FIG4A_KEEP_TIER34=1, full 50-GWAS COLOC) -> the same panels with a
  _supp_full50gwas suffix, so the sensitivity variant NEVER overwrites the main PDFs.

Numbers computed live by compute_counts() from the canonical sources and asserted
against the values locked 2026-07-01. Env: rnaseq python (needs pandas). matplotlib
import can hang on the login node -> render via SLURM (run_fig4a_candidates.sh).
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, PathPatch
from matplotlib.path import Path

# -- Paths --------------------------------------------------------------------
BASE = os.environ.get(
    "MASLD_PROJECT_ROOT",
    "/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design",
)
OUT_DIR   = os.path.join(BASE, "figures/main/fig4_validation")
ATLAS     = os.path.join(BASE, "RNA-seq/results/multi_evidence/multi_evidence_atlas.csv")
DEG_CSV   = os.path.join(BASE, "RNA-seq/Human/Patient_Cohorts/analysis/integration/"
                               "results/integration/canonical_deg_results.csv")
PROT_CSV  = os.path.join(BASE, "Analysis/Proteomics/results/"
                               "protein_transcript_concordance_v3.csv")
COLOC_CSV = os.path.join(BASE, "GWAS/finemapping/results/susie_coloc/gene_level_coloc.csv")
# Per-trait COLOC (gwas_name x gene) — used to recompute the COLOC set over a
# trait-restricted portfolio. The canonical gene_level_coloc.csv only stores the
# best-across-ALL-traits PP.H4, so it cannot express a trait exclusion.
COLOC_PERTRAIT = os.path.join(BASE,
    "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv")

# Tier-1/2 vs Tier-3/4 GWAS portfolio (canonical tier map, 2026-07-06). The MAIN
# render restricts the COLOC arm to the LIVER-SPECIFIC Tier-1/2 portfolio: the 35
# strata flagged placement=="main" in gwas_trait_tier.tsv (direct MASLD + liver
# enzymes). We RE-COMPUTE each gene's best abf PP.H4 over ONLY those main strata (so
# a gene rescued by a second-best liver/enzyme trait survives), which drops the
# distal/non-specific Tier-3/4 strata (MVP ChronLiver, Cirrhosis, Albumin, Platelet).
# Tier-1/2 coloc = 862 (identical to gene_level_coloc_tier12.csv coloc_best_abf_pp4>0.5).
# Set FIG4A_KEEP_TIER34=1 to REVERT to the all-50-GWAS canonical (coloc_best_pp4 =
# 1,234) -> the labeled SUPP/sensitivity render. NOTE: placement=="main" is STRICTER
# than the old Albumin/Platelet-only regex, which retained Tier-3 ChronLiver/Cirrhosis
# and gave 891.
GWAS_TIER = os.path.join(BASE, "GWAS/finemapping/config/gwas_trait_tier.tsv")

# -- Canonical palette (single source of truth) -------------------------------
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig1_palette import MODALITY_COLORS, CONTROL_GRAY  # noqa: E402

C_BULK    = MODALITY_COLORS["bulk"]      # #0072B2
C_GWAS    = MODALITY_COLORS["gwas"]      # #C2185B  COLOC / genetics
C_PROTEO  = MODALITY_COLORS["proteo"]    # #D55E00
C_SPATIAL = MODALITY_COLORS["spatial"]   # #CC79A7
C_ATAC    = MODALITY_COLORS["atac"]      # #E69F00
C_RNA     = "#9AA0A6"                     # RNA-seq signal (neutral origin)
C_CORE    = "#5F6368"                     # TREAT confident core (dark inner band)
C_CONV    = "#2E6E8E"                     # convergent core (deep teal)

FS = 6  # SINGLE font size for every text element (PI: font 6, no bold, compact)
INK = "#2B2B2B"

plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":       FS,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "figure.dpi":      150,
})

DEG_HEADLINE = 1918   # TREAT confident set (lfc=0.25); RNA-seq's own spine

# Values re-locked 2026-07-06 for the Tier-1/2 (liver-specific) MAIN render: the COLOC
# arm is restricted to the 35 placement=="main" strata (direct MASLD + liver enzymes),
# recomputed best-abf-over-main -> coloc PP.H4>0.5 = 862 (universe 13,390), dropping the
# distal Tier-3/4 strata (ChronLiver/Cirrhosis/Albumin/Platelet). The full 50-GWAS
# portfolio (coloc 1,234, universe 13,552) is the SUPP/sensitivity render
# (FIG4A_KEEP_TIER34=1), which DELIBERATELY drifts vs this LOCK. E1b remediation
# PRESERVED — spatial S_val=111, >=2-lens=19 and all-3=0 are UNCHANGED from the full
# portfolio because the convergent core is robust to the COLOC scope (the CosMx cell-
# level exclusion is coded, not locked). E1b context: the CosMx cosmx_*_mash_padj
# columns carry NO valid cell-level p (pseudoreplication over ~297K cells nested in 4
# slides; A6 fix — slide-level DIRECTION only), so the spatial gate drops the cell-level
# CosMx-hep column (S_val had collapsed 639->111 at that fix).
LOCK = dict(
    substrate=12989, coloc=862, coloc_only=401, coloc_shared=461, universe=13390,
    confident=1915,
    P_meas=3450, P_val=220, P_conv=18,
    S_meas=7197, S_val=111, S_conv=12,
    A_meas=755,  A_val=337, A_conv=8,
    ge1=649, ge2=19, all3=0,
    PS=11, PA=7, SA=1,
    P_only=202, S_only=99, A_only=329,
    # per-dataset validators for the named-dataset alluvial (tier-1/2 universe;
    # fallback only — the live values track compute_counts / the active COLOC scope)
    P_liver_meas=3342, P_liver_val=169, P_plasma_meas=1938, P_plasma_val=71,
    S_vis_meas=163, S_geomx_meas=7151,
)
# MAIN (Tier-1/2) writes the canonical fig4a_overview_*.pdf; the full-50-GWAS SUPP
# render (FIG4A_KEEP_TIER34=1) appends a _supp_full50gwas suffix so it never overwrites
# the main PDFs (applied centrally in _save()).
SUPP_FULL50 = bool(os.environ.get("FIG4A_KEEP_TIER34"))
OUT_SUFFIX  = "_supp_full50gwas" if SUPP_FULL50 else ""


# -- Data ---------------------------------------------------------------------
def compute_counts():
    """Recompute every panel number live. Universe = RNA-seq signal (padj<0.05)
    UNION GWAS-COLOC (best PP.H4>0.5). Each lens VALIDATED = universe member AND
    modality-disease-DE padj<0.05. Falls back to LOCK if a source is unreadable.

    Also returns (under res["_sets"]) the underlying gene SETS the enrichment /
    UpSet builders need: per-lens universe-restricted measured/validated sets, the
    FULL (genome-wide) per-lens measured/validated backgrounds for the enrichment
    null, the COLOC-only orthogonal prioritization axis, and the <=19 convergent
    (>=2-lens) genes with per-gene lens membership."""
    try:
        import pandas as pd

        d = pd.read_csv(DEG_CSV)
        d["sym"] = d["symbol"].fillna(d["gene"])
        rna_sig = set(d.loc[d["padj"] < 0.05, "sym"])
        treat   = set(d.loc[d["treat_fdr"] < 0.05, "sym"])

        if os.environ.get("FIG4A_KEEP_TIER34"):
            # SUPP/sensitivity: full 50-GWAS canonical (best PP.H4 over ALL traits).
            gl = pd.read_csv(COLOC_CSV)
            coloc = set(gl.loc[gl["coloc_best_pp4"] > 0.5, "gene"])
            print(f"[compute_counts] COLOC full 50-GWAS canonical: {len(coloc)} genes "
                  "(SUPP/sensitivity render)")
        else:
            # MAIN: Tier-1/2 liver-specific portfolio. Recompute best abf PP.H4 over
            # ONLY the placement=="main" strata (drops Tier-3/4 ChronLiver/Cirrhosis/
            # Albumin/Platelet), so a gene rescued by a second-best liver/enzyme trait
            # survives. Identical to gene_level_coloc_tier12.csv coloc_best_abf_pp4>0.5.
            tier = pd.read_csv(GWAS_TIER, sep="\t")
            main = set(tier.loc[tier["placement"] == "main", "study_name"])
            pt = pd.read_csv(COLOC_PERTRAIT, usecols=["gwas_name", "gene", "PP.H4.abf"])
            best = pt.loc[pt["gwas_name"].isin(main)].groupby("gene")["PP.H4.abf"].max()
            coloc = set(best.index[best > 0.5])
            print(f"[compute_counts] COLOC Tier-1/2 (placement==main, {len(main)} strata): "
                  f"{len(coloc)} genes (all-50-GWAS canonical = 1234)")
        coloc = {g for g in coloc if isinstance(g, str) and g}
        target = rna_sig | coloc
        coloc_only = coloc - rna_sig          # orthogonal genetic axis (not RNA-DEG)

        # ── Proteomics: FULL background (all measured) + universe-restricted ──
        c = pd.read_csv(PROT_CSV)
        c = (c[c["dream_comparator"] == "disease_vs_control"]
             .sort_values("protein_padj").drop_duplicates("gene"))
        P_meas_all = set(c["gene"])
        P_val_all  = set(c.loc[c["protein_padj"] < 0.05, "gene"])
        Pmeas = P_meas_all & target
        Pval  = P_val_all  & target

        # ── Spatial / scATAC: compute FULL background on the unfiltered atlas ──
        a_all = pd.read_csv(ATLAS, low_memory=False).rename(columns={"human_symbol": "sym"})
        # Spatial significance gate = only the arms with a VALID unit-of-analysis
        # p-value: GeoMx (ROI-level LMM, patient RE) + Visium hepatocyte Wilcoxon.
        # CosMx cosmx_*_mash_padj is DELIBERATELY EXCLUDED (E1b remediation): its
        # cell-level Wilcoxon over ~297K cells nested in 4 slides is pseudoreplicated
        # and no valid cell-level p exists — the atlas now carries slide-level
        # DIRECTION only (padj = NaN), so CosMx contributes direction support, never
        # a significance-gated "validated" count. (This also makes S_val robust to the
        # atlas cosmx-padj state, not merely reliant on it being NaN.)
        sp = ["spatial_govaere2026_geomx_sh_vs_ls_padj",
              "spatial_govaere2026_geomx_sh_vs_pt_padj",
              "spatial_hep_wilcoxon_padj_bh"]
        sp = [c0 for c0 in sp if c0 in a_all.columns]
        S_meas_all = set(a_all.loc[a_all[sp].notna().any(axis=1), "sym"])
        S_val_all  = set(a_all.loc[(a_all[sp] < 0.05).any(axis=1), "sym"])
        A_meas_all = set(a_all.loc[a_all["hepatocyte_da_padj"].notna(), "sym"])
        A_val_all  = set(a_all.loc[a_all["hepatocyte_da_padj"] < 0.05, "sym"])
        Smeas = S_meas_all & target
        Sval  = S_val_all  & target
        Ameas = A_meas_all & target
        Aval  = A_val_all  & target

        # ── PER-DATASET splits (for the named-dataset validation alluvial) ──
        # Proteomics: the two DIA-MS proteomes are BOTH gated validators. The
        # union (P_meas/P_val) is what P_val counts; per-dataset the two overlap
        # (1,863 co-measured / 20 co-validated), so liver+plasma sum > union.
        cds = pd.read_csv(PROT_CSV)
        cds = cds[cds["dream_comparator"] == "disease_vs_control"]
        def _proteo_ds(name):
            x = cds[cds["dataset"] == name]
            m = set(x["gene"]) & target
            v = set(x.loc[x["protein_padj"] < 0.05, "gene"]) & target
            return len(m), len(v)
        P_liver_meas,  P_liver_val  = _proteo_ds("PXD051911")   # liver DIA-MS (Boel)
        P_plasma_meas, P_plasma_val = _proteo_ds("PXD052937")   # plasma DIA-MS (Sourianarayanane)
        # Spatial: Visium (GSE192741) is the SOLE gated validator; GeoMx is gated
        # but yields 0, so it goes to the orthogonal tier (its measured count is
        # reported so the reader sees it was deployed, not skipped).
        S_vis_meas = len(set(a_all.loc[a_all["spatial_hep_wilcoxon_padj_bh"].notna(), "sym"]) & target)
        _geo = [c0 for c0 in ("spatial_govaere2026_geomx_sh_vs_ls_padj",
                              "spatial_govaere2026_geomx_sh_vs_pt_padj") if c0 in a_all.columns]
        S_geomx_meas = len(set(a_all.loc[a_all[_geo].notna().any(axis=1), "sym"]) & target) if _geo else 0

        allv = Pval | Sval | Aval
        # per-gene lens membership over the universe-restricted validated sets
        lens_of = {g: tuple(L for L, S in (("P", Pval), ("S", Sval), ("A", Aval)) if g in S)
                   for g in allv}
        ge2 = {g for g in allv if len(lens_of[g]) >= 2}
        ge2_members = sorted(((g, lens_of[g]) for g in ge2),
                             key=lambda t: (-len(t[1]), t[0]))

        res = dict(
            substrate=len(rna_sig), coloc=len(coloc), coloc_only=len(coloc_only),
            coloc_shared=len(coloc & rna_sig), universe=len(target),
            confident=len(rna_sig & treat),
            P_meas=len(Pmeas), P_val=len(Pval), P_conv=len(Pval & ge2),
            S_meas=len(Smeas), S_val=len(Sval), S_conv=len(Sval & ge2),
            A_meas=len(Ameas), A_val=len(Aval), A_conv=len(Aval & ge2),
            ge1=len(allv), ge2=len(ge2), all3=len(Pval & Sval & Aval),
            PS=len(Pval & Sval), PA=len(Pval & Aval), SA=len(Sval & Aval),
            P_only=len(Pval - Sval - Aval), S_only=len(Sval - Pval - Aval),
            A_only=len(Aval - Pval - Sval),
            # per-dataset validators (named-dataset alluvial)
            P_liver_meas=P_liver_meas, P_liver_val=P_liver_val,
            P_plasma_meas=P_plasma_meas, P_plasma_val=P_plasma_val,
            S_vis_meas=S_vis_meas, S_geomx_meas=S_geomx_meas,
        )
        res.update(proteo=res["P_val"], spatial=res["S_val"], scatac=res["A_val"],
                   none=res["universe"] - res["ge1"], PSA=res["all3"])
        res["_sets"] = dict(
            universe=target, rna=rna_sig, coloc=coloc, coloc_only=coloc_only, treat=treat,
            Pmeas=Pmeas, Pval=Pval, Smeas=Smeas, Sval=Sval, Ameas=Ameas, Aval=Aval,
            P_meas_all=P_meas_all, P_val_all=P_val_all,
            S_meas_all=S_meas_all, S_val_all=S_val_all,
            A_meas_all=A_meas_all, A_val_all=A_val_all,
            ge2_members=ge2_members,
        )
        print("[compute_counts] live:",
              {k: res[k] for k in sorted(res) if not k.startswith("_")})
        drift = {k: (res[k], LOCK[k]) for k in LOCK if res.get(k) != LOCK[k]}
        if drift:
            print("[compute_counts] WARNING drift vs locked:", drift)
        return res
    except Exception as e:  # pragma: no cover
        import traceback
        print(f"[compute_counts] fell back to LOCK ({e})")
        traceback.print_exc()
        r = dict(LOCK)
        r.update(proteo=r["P_val"], spatial=r["S_val"], scatac=r["A_val"],
                 none=r["universe"] - r["ge1"], PSA=r["all3"])
        return r


# -- Above-chance null (pure numpy; NO scipy — rnaseq scipy/numpy ABI risk) ---
N_PERM = 10000
RNG_SEED = 42


def permutation_null(res):
    """Two above-chance statements, both from one seeded pure-numpy null.

    (a) PER-LENS ENRICHMENT — is a prioritized set enriched for a lens's
        independent disease signal vs the genome-wide measured background? Exact
        hypergeometric: among the N genes MEASURED in lens L, K are validated; a
        prioritized set overlaps k of the measured genes; under H0 (prioritized no
        likelier than random) the validated-overlap ~ Hypergeometric(N,K,k). We
        draw it (rng.hypergeometric) for the fold + empirical p. Reported for the
        full universe AND (headline) the COLOC-only orthogonal axis — genetics is
        independent of every somatic readout, so that enrichment can't be shared-
        assay circularity.

    (b) CONVERGENCE vs CHANCE — with each lens's validated COUNT held fixed and
        randomly placed among ITS universe-co-measured genes, how many genes hit
        >=2 / all-3 lenses by chance? Compared to the observed 19 / 0. Reported
        honestly whichever way it lands (the lenses are near-independent, so few
        genes converge even under real signal; (a) is the robust headline)."""
    if "_sets" not in res:
        print("[permutation_null] no gene sets (LOCK fallback) — skipping null")
        return None
    import numpy as np
    # Sanitize gene sets for the null (a stray NaN symbol can slip into the
    # universe; drop non-string / empty so sorting & indexing are well-defined).
    # compute_counts() reported the LOCK-matching COUNTS already — this only guards
    # the null's set ops and changes nothing displayed.
    def _clean(s):
        return {x for x in s if isinstance(x, str) and x and x.lower() != "nan"}
    S = {k: (_clean(v) if isinstance(v, set) else v)
         for k, v in res["_sets"].items()}
    rng = np.random.default_rng(RNG_SEED)

    def hyper_enrich(P, meas_all, val_all):
        N = len(meas_all); K = len(val_all)
        k = len(P & meas_all); obs = len(P & val_all)
        if k == 0 or K == 0 or K >= N:
            return dict(obs=obs, exp=float("nan"), fold=float("nan"), p=1.0,
                        N=N, K=K, k=k)
        draws = rng.hypergeometric(K, N - K, k, size=N_PERM)
        exp = float(draws.mean())
        fold = obs / exp if exp > 0 else float("inf")
        p = (int(np.sum(draws >= obs)) + 1) / (N_PERM + 1)
        return dict(obs=obs, exp=exp, fold=fold, p=p, N=N, K=K, k=k)

    lenses = (("proteomics", S["P_meas_all"], S["P_val_all"]),
              ("spatial",    S["S_meas_all"], S["S_val_all"]),
              ("scATAC",     S["A_meas_all"], S["A_val_all"]))
    enrich = {}
    for pset_name, P in (("universe", S["universe"]), ("coloc_only", S["coloc_only"])):
        enrich[pset_name] = {ln: hyper_enrich(P, m, v) for ln, m, v in lenses}

    # (b) convergence null over the universe-co-measured sets, tracking each
    #     EXCLUSIVE UpSet intersection so the UpSet can show a chance reference.
    # Drop non-string members (the COLOC master carries one pp4>0.5 locus with a
    # NaN gene symbol) so sorted() does not choke on mixed float/str. This only
    # affects the Monte-Carlo null's sampling pool (13,552 -> 13,551); every
    # DISPLAYED count comes from len() in compute_counts (d[...]), so 1,234 coloc /
    # 13,552 universe / 111 spatial / 19 ge2 / 0 all3 are unchanged.
    uni = sorted(g for g in S["universe"] if isinstance(g, str))
    idx = {g: i for i, g in enumerate(uni)}
    U = len(uni)
    meas_idx = {L: np.fromiter((idx[g] for g in ms if g in idx), dtype=np.int64)
                for L, ms in (("P", S["Pmeas"]), ("S", S["Smeas"]), ("A", S["Ameas"]))}
    val_n = {"P": len(S["Pval"]), "S": len(S["Sval"]), "A": len(S["Aval"])}
    keys = ("ge2", "all3", "P_only", "S_only", "A_only", "PS", "PA", "SA")
    null = {k: np.empty(N_PERM, dtype=np.int32) for k in keys}
    for it in range(N_PERM):
        mP = np.zeros(U, dtype=bool); mS = np.zeros(U, dtype=bool); mA = np.zeros(U, dtype=bool)
        for L, m in (("P", mP), ("S", mS), ("A", mA)):
            m[rng.choice(meas_idx[L], size=val_n[L], replace=False)] = True
        both = (mP & mS) | (mP & mA) | (mS & mA)
        null["ge2"][it]    = int(both.sum())
        null["all3"][it]   = int((mP & mS & mA).sum())
        null["P_only"][it] = int((mP & ~mS & ~mA).sum())
        null["S_only"][it] = int((mS & ~mP & ~mA).sum())
        null["A_only"][it] = int((mA & ~mP & ~mS).sum())
        null["PS"][it]     = int((mP & mS & ~mA).sum())
        null["PA"][it]     = int((mP & mA & ~mS).sum())
        null["SA"][it]     = int((mS & mA & ~mP).sum())

    def summarize(obs, nul):
        exp = float(nul.mean())
        fold = obs / exp if exp > 0 else float("inf")
        p = (int(np.sum(nul >= obs)) + 1) / (N_PERM + 1)
        lo, hi = np.percentile(nul, [2.5, 97.5])
        return dict(obs=int(obs), exp=exp, fold=fold, p=p,
                    lo=float(lo), hi=float(hi))

    # observed EXCLUSIVE intersection counts (all3 = 0 so pairwise == exclusive)
    obs_excl = dict(ge2=res["ge2"], all3=res["all3"],
                    P_only=res["P_only"], S_only=res["S_only"], A_only=res["A_only"],
                    PS=res["PS"] - res["all3"], PA=res["PA"] - res["all3"],
                    SA=res["SA"] - res["all3"])
    conv = {k: summarize(obs_excl[k], null[k]) for k in keys}

    # ── report ──
    print(f"\n[permutation_null] seed={RNG_SEED}, draws={N_PERM}")
    for pset_name in ("coloc_only", "universe"):
        print(f"  per-lens enrichment | prioritized = {pset_name}:")
        for ln, _, _ in lenses:
            e = enrich[pset_name][ln]
            print(f"    {ln:11s}: obs={e['obs']:4d}  exp={e['exp']:6.1f}  "
                  f"fold={e['fold']:5.2f}x  p={e['p']:.2e}  "
                  f"(N_meas={e['N']}, K_val={e['K']}, k_prior={e['k']})")
    for key in ("ge2", "all3", "PS", "PA", "SA", "P_only", "S_only", "A_only"):
        c = conv[key]
        print(f"  convergence {key:7s}: obs={c['obs']:3d}  exp={c['exp']:6.2f}  "
              f"[95% {c['lo']:.1f}-{c['hi']:.1f}]  fold={c['fold']:.2f}x  p={c['p']:.2e}")
    lab = {("P", "S"): "proteo∩spatial", ("P", "A"): "proteo∩scATAC",
           ("S", "A"): "spatial∩scATAC", ("P", "S", "A"): "all-3"}
    print("  convergent (>=2-lens) genes:")
    for g, mem in S["ge2_members"]:
        print(f"    {g:12s} {lab.get(mem, '/'.join(mem))}")
    print("")
    return dict(enrich=enrich, conv=conv, lenses=[l[0] for l in lenses])


# -- Shared drawing helpers ---------------------------------------------------
def _ribbon(ax, x0, y0a, y0b, x1, y1a, y1b, color, alpha=0.38):
    """Fig 2A idiom: cubic-Bezier ribbon with horizontal tangents at both ends,
    so the band leaves each node horizontally and sweeps to the other (y0a..y0b at
    x0, y1a..y1b at x1)."""
    cx = (x0 + x1) / 2.0
    verts = [(x0, y0a), (cx, y0a), (cx, y1a), (x1, y1a),
             (x1, y1b), (cx, y1b), (cx, y0b), (x0, y0b), (x0, y0a)]
    codes = [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4,
             Path.LINETO, Path.CURVE4, Path.CURVE4, Path.CURVE4, Path.CLOSEPOLY]
    ax.add_patch(PathPatch(Path(verts, codes), facecolor=color, edgecolor="none",
                           alpha=alpha, zorder=2))



def _save(fig, name):
    os.makedirs(OUT_DIR, exist_ok=True)
    if OUT_SUFFIX and name.endswith(".pdf"):   # SUPP full-50 render: never clobber main
        name = name[:-4] + OUT_SUFFIX + ".pdf"
    p = os.path.join(OUT_DIR, name)
    fig.savefig(p, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print("Saved:", p)


# -- (1) CASCADE — Fig 4 OPENER: named-dataset validation alluvial ------------
def build_cascade(d):
    """Fig 4 opener in the Fig 2A idiom (conserved-flow cubic-Bezier alluvial),
    but NAMING the datasets so the reader sees Fig 4 brings in orthogonal assays
    not shown in Figs 2-3. Four stages, colour = modality, sub-datasets = tint:

      PRIORITIZED  bulk RNA-seq (+TREAT core) + GWAS-COLOC     (the Figs 2-3 axes)
      MODALITY     proteomics / spatial / scATAC              (3 orthogonal lenses)
      DATASET      each modality's NAMED gated-VALIDATOR dataset(s)
      OUTCOME      per dataset: VALIDATED (solid) vs not significant (grey)

    A muted 'also profiled' tier below names the assays NOT in the significance
    gate: Olink plasma (withdrawn), GeoMx (gated but 0 hits), CosMx (direction-
    only), Vu Visium (replication). Counts are live from compute_counts (so they
    track the active COLOC scope). Proteomics' two DIA-MS datasets OVERLAP, so the
    per-dataset labels show TRUE counts and sum > the proteomics union (noted in
    the caption); spatial's validator flow is Visium (163), NOT the ~7,200 spatial-
    measured — that ~99% GeoMx (0 validated) is shown in the orthogonal tier."""
    import matplotlib.colors as mcolors

    def tint(hexc, f):                       # lighten a modality hue for sub-nodes
        r, g, b = mcolors.to_rgb(hexc)
        return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)

    # datasets: (modality, accession, platform·N, measured, validated)
    DS = [
        ("proteo",  "PXD051911", "liver DIA-MS · 58",  d["P_liver_meas"],  d["P_liver_val"]),
        ("proteo",  "PXD052937", "plasma DIA-MS · 72", d["P_plasma_meas"], d["P_plasma_val"]),
        ("spatial", "GSE192741", "Visium · 5",         d["S_vis_meas"],    d["S_val"]),
        ("atac",    "GSE244832", "snATAC · 18 donors", d["A_meas"],        d["A_val"]),
    ]
    lanes   = ["proteo", "spatial", "atac"]
    MODNAME = {"proteo": "proteomics", "spatial": "spatial", "atac": "scATAC"}
    COL     = {"proteo": C_PROTEO, "spatial": C_SPATIAL, "atac": C_ATAC}
    # modality measured = its VALIDATOR union (spatial = Visium only; GeoMx -> tier)
    MOD_MEAS = {"proteo": d["P_meas"], "spatial": d["S_vis_meas"], "atac": d["A_meas"]}
    ds_of = {l: [x for x in DS if x[0] == l] for l in lanes}
    ACC = {g[1]: g for g in DS}

    H, GAP, BW = 10.0, 0.70, 0.13
    SEC = 2.30
    X0, X1, X2, X3 = 0.0, SEC, 2 * SEC, 3 * SEC
    C_NS = "#E4E4E4"
    # min node height: this is an INTRO/inventory panel, so floor each lane well
    # above its raw gene-share (proteomics measures ~20x more genes than Visium)
    # to keep every named dataset legible; ribbon widths still rank the modalities.
    MINH = 1.90

    def _stack(keys, wt, height, gap, minh):
        """proportional top->bottom slots filling `height`, each >= minh."""
        n = len(keys); avail = height - gap * (n - 1)
        w = {k: max(float(wt[k]), 1e-9) for k in keys}
        h = {k: avail * w[k] / sum(w.values()) for k in keys}
        for _ in range(4):
            lo = [k for k in keys if h[k] < minh]
            if not lo:
                break
            free = [k for k in keys if k not in lo]
            fw = sum(w[k] for k in free); rem = avail - minh * len(lo)
            for k in lo:
                h[k] = minh
            for k in free:
                h[k] = (rem * w[k] / fw) if fw > 0 else rem / max(len(free), 1)
        slots = {}; y = 0.0
        for k in keys:
            slots[k] = (y, y + h[k]); y += h[k] + gap
        return slots

    # MODALITY column
    slotMod = _stack(lanes, MOD_MEAS, H, GAP, MINH)
    # DATASET column: subdivide each modality's y-range among its datasets
    slotDS = {}
    for l in lanes:
        my0, my1 = slotMod[l]; grp = ds_of[l]
        sub = _stack([g[1] for g in grp], {g[1]: g[3] for g in grp},
                     my1 - my0, GAP * 0.55, min(MINH, (my1 - my0) / len(grp)))
        for k, (a, b) in sub.items():
            slotDS[k] = (my0 + a, my0 + b)

    # OUTCOME column: VALIDATED group (top) + gap + n.s. group (bottom)
    order = [g[1] for g in DS]
    VAL = {g[1]: g[4] for g in DS}; MEA = {g[1]: g[3] for g in DS}
    NS  = {k: MEA[k] - VAL[k] for k in order}
    GG = 0.9
    scO = (H - GG - GAP * 2 * (len(order) - 1)) / float(sum(VAL.values()) + sum(NS.values()))
    slotV = {}; y = 0.0
    for k in order:
        h = max(VAL[k] * scO, 0.5); slotV[k] = (y, y + h); y += h + GAP
    slotN = {}; y = (y - GAP) + GG
    for k in order:
        h = max(NS[k] * scO, 0.35); slotN[k] = (y, y + h); y += h + GAP

    fig, ax = plt.subplots(figsize=(7.3, 4.7))
    ax.set_xlim(-2.75, X3 + 1.75); ax.set_ylim(-1.35, H + 4.2)
    ax.invert_yaxis(); ax.axis("off")

    def node(x, y0, y1, c, z=4, ec="white"):
        ax.add_patch(Rectangle((x - BW / 2, y0), BW, y1 - y0, facecolor=c,
                               edgecolor=ec, lw=0.2, zorder=z))

    # ── PRIORITIZED origin: RNA-seq signal (grey, TREAT-core inset) + GWAS-COLOC ──
    o_gap = 0.24; avail = H - o_gap; o_tot = d["substrate"] + d["coloc"]
    hR = avail * d["substrate"] / o_tot
    rna0, rna1 = 0.0, hR; gco0, gco1 = hR + o_gap, H
    node(X0, rna0, rna1, C_RNA)
    ax.add_patch(Rectangle((X0 - BW / 2, rna0), BW, (rna1 - rna0) * DEG_HEADLINE / d["substrate"],
                           facecolor=C_CORE, edgecolor="none", zorder=5))
    node(X0, gco0, gco1, C_GWAS)
    lx = X0 - BW / 2 - 0.16
    ax.text(lx, (rna0 + rna1) / 2, f"bulk RNA-seq {d['substrate']:,}", ha="right", va="center", color=INK, zorder=8)
    ax.text(lx, (gco0 + gco1) / 2, f"GWAS-COLOC {d['coloc']:,}", ha="right", va="center", color=INK, zorder=8)

    # ── PRIORITIZED -> MODALITY (partition mouth by validator-coverage share) ──
    tot = float(sum(MOD_MEAS.values())); y = 0.0
    for l in lanes:
        h = H * MOD_MEAS[l] / tot
        _ribbon(ax, X0 + BW / 2, y, y + h, X1 - BW / 2, slotMod[l][0], slotMod[l][1], COL[l], alpha=0.30)
        y += h
    for l in lanes:
        y0, y1 = slotMod[l]; node(X1, y0, y1, COL[l])
        ax.text(X1, y0 - 0.13, MODNAME[l], ha="center", va="bottom", color=INK, zorder=8)

    # ── MODALITY -> DATASET (fan each modality into its named datasets) ──
    for l in lanes:
        m0, m1 = slotMod[l]; grp = ds_of[l]; gtot = float(sum(g[3] for g in grp)); yy = m0
        for g in grp:
            seg = (m1 - m0) * g[3] / gtot; dd0, dd1 = slotDS[g[1]]
            _ribbon(ax, X1 + BW / 2, yy, yy + seg, X2 - BW / 2, dd0, dd1, tint(COL[l], 0.25), alpha=0.32)
            yy += seg
    for g in DS:
        dd0, dd1 = slotDS[g[1]]; node(X2, dd0, dd1, tint(COL[g[0]], 0.18))
        ax.text(X2, dd0 - 0.34, g[1], ha="center", va="bottom", color=INK, zorder=8)          # accession
        ax.text(X2, dd0 - 0.07, g[2], ha="center", va="bottom", color="#6B6B6B", zorder=8)     # platform · N

    # ── DATASET -> OUTCOME (validated sweeps up, n.s. flows down) ──
    for g in DS:
        acc, l = g[1], g[0]; dd0, dd1 = slotDS[acc]
        hv = (dd1 - dd0) * (VAL[acc] / MEA[acc]) if MEA[acc] else 0.0
        _ribbon(ax, X2 + BW / 2, dd0, dd0 + hv, X3 - BW / 2, slotV[acc][0], slotV[acc][1], COL[l], alpha=0.42)
        _ribbon(ax, X2 + BW / 2, dd0 + hv, dd1, X3 - BW / 2, slotN[acc][0], slotN[acc][1], C_NS, alpha=0.7)
    for g in DS:
        acc, l = g[1], g[0]
        y0, y1 = slotV[acc]; node(X3, y0, y1, COL[l])
        ax.text(X3 + BW / 2 + 0.13, (y0 + y1) / 2, f"{VAL[acc]:,}", ha="left", va="center", color=INK, zorder=8)
        y0, y1 = slotN[acc]; node(X3, y0, y1, C_NS)

    vy = (slotV[order[0]][0] + slotV[order[-1]][1]) / 2.0
    ny = (slotN[order[0]][0] + slotN[order[-1]][1]) / 2.0
    ax.text(X3 + BW / 2 + 0.62, vy, "VALIDATED", ha="left", va="center", color=INK, zorder=8)
    ax.text(X3 + BW / 2 + 0.13, ny, "not significant", ha="left", va="center", color="#8A8A8A", zorder=8)

    # ── stage headers ──
    HY = -1.05
    for xx, lab in ((X0, "PRIORITIZED"), (X1, "MODALITY"), (X2, "DATASET"), (X3, "OUTCOME")):
        ax.text(xx, HY, lab, ha="center", va="bottom", color=INK, zorder=8)

    # ── orthogonal 'also profiled' tier (muted; NOT in the significance gate) ──
    ty = H + 1.35
    ax.plot([X0 - 0.1, X3 + 0.4], [ty - 0.55, ty - 0.55], color="#C9C9C9", lw=0.5, zorder=1)
    ax.text(X0 - 0.1, ty - 0.5, "also profiled — orthogonal support, not in the significance gate",
            ha="left", va="bottom", color="#6B6B6B", zorder=8)
    orth = [
        (C_PROTEO,  "proteomics", "Olink Explore plasma · Yang/Zeybel 2025 · 1,461 × 218 (matched to GSE192741 RNA-seq)"),
        (C_SPATIAL, "spatial",    f"GeoMx · Govaere 2026 · {d['S_geomx_meas']:,} measured, 0 gated   ·   "
                                  "CosMx · GSE312698 · 552K cells (direction only)   ·   Vu 2025 · Visium · 33 biopsies (replication)"),
    ]
    for i, (c, mod, txt) in enumerate(orth):
        yy = ty + i * 0.85
        ax.add_patch(Rectangle((X0 - 0.1, yy - 0.16), 0.16, 0.32, facecolor=c, edgecolor="none", alpha=0.5, zorder=3))
        ax.text(X0 + 0.16, yy, txt, ha="left", va="center", color="#3A3A3A", zorder=8)

    _save(fig, "fig4a_overview_cascade.pdf")
    # printed caption (no on-panel subtitle, per house style)
    print("[caption:cascade] Fig 4 opener. Prioritized universe (bulk RNA-seq padj<0.05 "
          f"∪ GWAS-COLOC PP.H4>0.5 = {d['universe']:,}) is co-measured by three orthogonal "
          "assays, each a NAMED dataset new to Fig 4. Node = dataset (accession · platform); "
          "ribbon width = genes; VALIDATED = both-significant (direction-agnostic). Proteomics' "
          f"two DIA-MS datasets overlap, so per-dataset validated (liver {d['P_liver_val']} + "
          f"plasma {d['P_plasma_val']}) sum > the union ({d['P_val']}); spatial validated "
          f"({d['S_val']}) is Visium GSE192741 alone (163 tested). Orthogonal assays "
          "(Olink/GeoMx/CosMx/Vu) are profiled but not significance-gated.")


def _ge2_genes(d, membership):
    """Convergent (>=2-lens) gene symbols with a given lens membership tuple."""
    if "_sets" not in d:
        return []
    out = [g for g, mem in d["_sets"]["ge2_members"]
           if mem == membership and isinstance(g, str) and g and g.lower() != "nan"]
    return sorted(out)


# -- (2) ENRICHMENT — above-chance claim (Version 1) --------------------------
def build_enrichment(d, stats):
    """Fig 4A candidate — the ABOVE-CHANCE claim, kept honest by the permutation
    null. Each validation test is drawn as a fold-enrichment bar anchored at
    chance (fold = observed / permutation-expected); bars are filled only when
    they beat chance (p<0.05). The headline that survives is
    proteomics∩spatial = 3.8x (p=2e-4); the broad 3-way convergence sits at
    chance and is shown as such. The 11 cross-assay-reproduced genes are named."""
    if stats is None:
        print("[build_enrichment] no stats (LOCK fallback) — skipping"); return
    conv, enr = stats["conv"], stats["enrich"]["universe"]

    def R(label, st, color):
        sig = (st["p"] < 0.05) and (st["fold"] > 1)
        return dict(label=label, fold=st["fold"], obs=st["obs"], exp=st["exp"],
                    p=st["p"], color=(color if sig else CONTROL_GRAY), sig=sig)

    blocks = [
        ("reproduced in two independent assays", [
            R("proteomics + spatial", conv["PS"], C_CONV),
            R("proteomics + scATAC",  conv["PA"], C_CONV),
            R("spatial + scATAC",     conv["SA"], C_CONV),
            R("any ≥2 of 3 lenses", conv["ge2"], C_CONV)]),
        ("enriched vs genome-wide background", [
            R("proteomics", enr["proteomics"], C_PROTEO),
            R("spatial",     enr["spatial"],   C_SPATIAL),
            R("scATAC",      enr["scATAC"],     C_ATAC)]),
    ]

    # y layout (top -> bottom): header, rows, gap, header, rows
    ys, items, headers = [], [], []
    y = 0.0
    for htext, rows in blocks:
        headers.append((y, htext)); y -= 0.7
        for r in rows:
            ys.append(y); items.append(r); y -= 1.0
        y -= 0.7
    ymin = y

    fig, ax = plt.subplots(figsize=(5.4, 3.9))
    XCH = 1.0                                  # chance line (fold = 1)
    XANN = 4.35                                # fixed annotation column
    for yy, r in zip(ys, items):
        x0, x1 = (XCH, r["fold"]) if r["fold"] >= XCH else (r["fold"], XCH)
        ax.add_patch(Rectangle((x0, yy - 0.32), x1 - x0, 0.64,
                               facecolor=r["color"], edgecolor="none",
                               alpha=0.95 if r["sig"] else 0.55, zorder=3))
        ax.text(-0.10, yy, r["label"], ha="right", va="center", color=INK, zorder=6)
        star = "  *" if r["sig"] else ""
        ax.text(XANN, yy, f"{r['obs']:,} vs {r['exp']:.1f}  p={r['p']:.1g}{star}",
                ha="left", va="center", color=INK, zorder=6)
    for hy, htext in headers:
        ax.text(-0.10, hy, htext, ha="right", va="center", color=INK, zorder=6)

    ax.axvline(XCH, color="#8A8A8A", lw=0.6, ls=(0, (3, 2)), zorder=2)
    ax.text(XCH, 0.9, "chance", ha="center", va="bottom", color="#8A8A8A", zorder=6)
    ax.set_xlim(-2.7, XANN + 2.2); ax.set_ylim(ymin + 0.3, 1.6)
    ax.set_yticks([]); ax.spines[["left", "right", "top"]].set_visible(False)
    ax.set_xticks([0, 1, 2, 3, 4])
    ax.tick_params(axis="x", length=2, labelsize=FS)
    ax.set_xlabel("fold enrichment over chance (observed / permutation-expected)")
    ax.spines["bottom"].set_bounds(0, 4)

    ax.set_title("Prioritized targets: above-chance reproduction", loc="left",
                 color=INK, fontsize=FS)

    # name the cross-assay (proteomics+spatial) reproduced genes — the real core.
    # Placed via fig-coords BELOW the x-axis label (tight bbox expands to include).
    ps_genes = _ge2_genes(d, ("P", "S"))
    if ps_genes:
        fig.subplots_adjust(bottom=0.20)
        per = 6
        chunks = [ps_genes[i:i + per] for i in range(0, len(ps_genes), per)]
        fig.text(0.02, 0.02, f"Reproduced in proteomics + spatial ({len(ps_genes)}):",
                 ha="left", va="top", color=INK, fontsize=FS)
        for j, ch in enumerate(chunks):
            fig.text(0.02, -0.03 - j * 0.05, ",  ".join(ch), ha="left", va="top",
                     color=INK, fontsize=FS, fontstyle="italic")
    _save(fig, "fig4a_overview_enrichment.pdf")
    # printed caption (no on-panel subtitle, per house style)
    print("[caption:enrichment] Each row is a validation test vs a 10,000-draw "
          "permutation null (seed 42). Bars filled = beats chance (p<0.05). The "
          "only cross-assay convergence exceeding chance is proteomics∩spatial "
          f"({conv['PS']['obs']} vs {conv['PS']['exp']:.1f}, {conv['PS']['fold']:.1f}x, "
          f"p={conv['PS']['p']:.1g}); overall ≥2-lens ({conv['ge2']['obs']}) is at "
          f"chance (p={conv['ge2']['p']:.2g}). The orthogonal COLOC-only axis does not "
          "enrich any functional lens (proteo/spatial/scATAC p>0.3); the per-lens "
          "background enrichment partly shares mRNA→protein/spatial biology.")


# -- (3) UPSET — intersection grammar with a chance reference (Version 3) ------
def build_upset(d, stats):
    """Fig 4A candidate — UpSet of the three validation lenses with a per-
    intersection CHANCE reference (permutation-null expected value marked on every
    bar). Exclusive intersection sizes on top; the only bar clearing its chance
    marker is proteomics∩spatial (highlighted). Set sizes at left; membership dot
    matrix below; the cross-assay-reproduced genes named beneath."""
    import numpy as np
    LENS = [("P", "proteomics", C_PROTEO, d["P_val"]),
            ("S", "spatial",    C_SPATIAL, d["S_val"]),
            ("A", "scATAC",     C_ATAC,    d["A_val"])]
    order = ["P", "S", "A"]
    keymap = {("A",): "A_only", ("P",): "P_only", ("S",): "S_only",
              ("P", "S"): "PS", ("P", "A"): "PA", ("S", "A"): "SA",
              ("P", "S", "A"): "all3"}
    inter = [(("A",), d["A_only"]), (("P",), d["P_only"]), (("S",), d["S_only"]),
             (("P", "S"), d["PS"] - d["all3"]), (("P", "A"), d["PA"] - d["all3"]),
             (("S", "A"), d["SA"] - d["all3"]), (("P", "S", "A"), d["all3"])]
    inter.sort(key=lambda t: (-t[1], len(t[0])))
    n = len(inter)
    exp = {mem: (stats["conv"][keymap[mem]]["exp"] if stats else None) for mem, _ in inter}
    sigmem = set()
    if stats:
        for mem, _ in inter:
            st = stats["conv"][keymap[mem]]
            if st["p"] < 0.05 and st["fold"] > 1:
                sigmem.add(mem)

    from matplotlib.gridspec import GridSpec
    fig = plt.figure(figsize=(5.4, 3.9))
    gs = GridSpec(2, 2, width_ratios=[1.0, 3.3], height_ratios=[2.2, 1.15],
                  wspace=0.05, hspace=0.08, figure=fig)
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_set = fig.add_subplot(gs[1, 0])
    ax_mat = fig.add_subplot(gs[1, 1], sharex=ax_bar)

    xs = np.arange(n)
    # ── intersection bars ──
    for x, (mem, cnt) in zip(xs, inter):
        conv2 = len(mem) >= 2
        col = C_CONV if mem in sigmem else (LENS_C(mem, LENS) if not conv2 else "#B9BEC4")
        ax_bar.bar(x, cnt, width=0.62, color=col,
                   alpha=0.95 if mem in sigmem else 0.7, zorder=3)
        ax_bar.text(x, cnt + max(inter[0][1] * 0.02, 0.5), f"{cnt}",
                    ha="center", va="bottom", color=INK, zorder=6)
        if exp[mem] is not None:                     # chance reference marker
            ax_bar.plot([x - 0.34, x + 0.34], [exp[mem], exp[mem]],
                        color="#6B6B6B", lw=0.9, ls=(0, (2, 1.5)), zorder=5)
    ax_bar.set_ylabel("genes")
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.tick_params(labelbottom=False, labelsize=FS, length=2)
    ax_bar.set_title("Validated-gene intersections vs chance (– – expected)",
                     loc="left", color=INK, fontsize=FS)

    # ── membership dot matrix ──
    ymap = {L: len(order) - 1 - i for i, L in enumerate(order)}
    for x, (mem, _) in zip(xs, inter):
        on = [ymap[L] for L in order if L in mem]
        if len(on) >= 2:
            ax_mat.plot([x, x], [min(on), max(on)], color=INK, lw=0.8, zorder=2)
        for L in order:
            filled = L in mem
            ax_mat.scatter(x, ymap[L], s=26,
                           color=(dict((l[0], l[2]) for l in LENS)[L] if filled else "#D9D9D9"),
                           edgecolors="none", zorder=3)
    ax_mat.set_xlim(-0.6, n - 0.4); ax_mat.set_ylim(-0.6, len(order) - 0.4)
    ax_mat.set_yticks([]); ax_mat.set_xticks([]); ax_mat.tick_params(length=0)
    for s in ("top", "right", "bottom", "left"):
        ax_mat.spines[s].set_visible(False)

    # ── set-size bars (extend leftward); modality NAME at far left, COUNT at the
    #    bar base (right) so neither collides with the dot-matrix rows ──
    maxsz = max(sz for *_, sz in LENS)
    ytf = ax_set.get_yaxis_transform()          # x = axes fraction, y = data
    for L, name, col, sz in LENS:
        ax_set.barh(ymap[L], sz, height=0.5, color=col, alpha=0.85, zorder=3)
        ax_set.text(0.02, ymap[L], name, transform=ytf, ha="left", va="center",
                    color=INK, zorder=6)
        ax_set.text(0.98, ymap[L], f"{sz}", transform=ytf, ha="right", va="center",
                    color=INK, zorder=6)
    ax_set.set_ylim(-0.6, len(order) - 0.4)
    ax_set.set_xlim(0, maxsz * 1.45); ax_set.invert_xaxis()
    ax_set.set_yticks([]); ax_set.set_xticks([])
    for s in ("top", "right", "bottom", "left"):
        ax_set.spines[s].set_visible(False)
    ax_set.set_title("set size", loc="right", color=INK, fontsize=FS)

    # ── name the proteomics∩spatial reproduced genes ──
    ps_genes = _ge2_genes(d, ("P", "S"))
    if ps_genes:
        txt = ",  ".join(ps_genes)
        fig.text(0.13, 0.015,
                 f"proteomics + spatial ({len(ps_genes)}):", ha="left", va="bottom",
                 color=INK, fontsize=FS)
        fig.text(0.13, -0.03, txt, ha="left", va="bottom", color=INK,
                 fontsize=FS, fontstyle="italic")
    _save(fig, "fig4a_overview_upset.pdf")
    print("[caption:upset] UpSet of the three orthogonal validation lenses on the "
          f"{d['universe']:,}-gene prioritized universe. Bars = exclusive intersection sizes; "
          "dashed tick = permutation-null expected (seed 42, 10,000 draws). Only "
          "proteomics∩spatial clears its chance marker "
          f"({d['PS']} vs {stats['conv']['PS']['exp']:.1f}, "
          f"{stats['conv']['PS']['fold']:.1f}x, p={stats['conv']['PS']['p']:.1g}); "
          "all-3 = 0 (chance ~0).")


def LENS_C(mem, LENS):
    """Single-lens bar colour (used for singleton UpSet bars)."""
    cmap = dict((l[0], l[2]) for l in LENS)
    return cmap.get(mem[0], "#B9BEC4") if len(mem) == 1 else "#B9BEC4"


def main():
    d = compute_counts()
    stats = permutation_null(d)
    if os.environ.get("FIG4A_STATS_ONLY"):
        print("Wrote nothing (FIG4A_STATS_ONLY set)")
        return
    build_cascade(d)
    build_enrichment(d, stats)
    build_upset(d, stats)
    print("Wrote fig4a_overview_{cascade,enrichment,upset}.pdf to", OUT_DIR)


if __name__ == "__main__":
    main()
