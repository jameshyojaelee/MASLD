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
OUT_DIR   = os.path.join(BASE, "figures/main/fig4_validation", "panels")
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
INK = "#231f20"  # 2026-07-09 house re-skin: near-black ink (was #2B2B2B)

plt.rcParams.update({
    "font.family":     "sans-serif",
    "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
    "font.size":       FS,
    "pdf.fonttype":    42,
    "ps.fonttype":     42,
    "figure.dpi":      150,
})

DEG_HEADLINE = 1918   # TREAT confident set (lfc=0.25); RNA-seq's own spine

# NO hard-coded "lock" of expected counts. The atlas / universe / DE sources change
# legitimately over time (new datasets, re-processed contrasts, broadened universe),
# so a frozen reference only goes stale and can silently mislabel a panel. compute_counts()
# is the SINGLE source of truth — every number is recomputed live from the source files,
# and if a source is unreadable it RAISES (fails the render) rather than falling back to
# stale numbers with current labels.
# MAIN (Tier-1/2) writes the canonical fig4a_overview_*.pdf; the full-50-GWAS SUPP
# render (FIG4A_KEEP_TIER34=1) appends a _supp_full50gwas suffix so it never overwrites
# the main PDFs (applied centrally in _save()).
SUPP_FULL50 = bool(os.environ.get("FIG4A_KEEP_TIER34"))
OUT_SUFFIX  = "_supp_full50gwas" if SUPP_FULL50 else ""


# -- Data ---------------------------------------------------------------------
def compute_counts():
    """Recompute every panel number live. Universe = RNA-seq signal (padj<0.05)
    UNION GWAS-COLOC (best PP.H4>0.5). Each lens VALIDATED = universe member AND
    modality-disease-DE padj<0.05. RAISES if a source is unreadable (no stale fallback).

    Also returns (under res["_sets"]) the underlying gene SETS the enrichment /
    UpSet builders need: per-lens universe-restricted measured/validated sets, the
    FULL (genome-wide) per-lens measured/validated backgrounds for the enrichment
    null, the COLOC-only orthogonal prioritization axis, and the <=19 convergent
    (>=2-lens) genes with per-gene lens membership."""
    try:
        import pandas as pd

        d = pd.read_csv(DEG_CSV)
        d["sym"] = d["symbol"].fillna(d["gene"])
        treat   = set(d.loc[d["treat_fdr"] < 0.05, "sym"])   # bulk TREAT confident core (1,915)

        # ── BROADENED, TREAT-consistent PRIORITIZED universe (locked 2026-07-06) ──
        # The first column is NOT just bulk padj<0.05. It is the union of ALL
        # transcriptomic candidates + ALL genetic candidates:
        #   transcriptomic (rna_sig, 8,088) = every DE contrast (bulk disease-vs-control,
        #     stage steatosis/SH/cirrhosis, fibrosis-gradient F0->F4, MASH-vs-MASL, sc
        #     pseudobulk per cell type) gated by the SAME interval-null TREAT @ lfc=0.25
        #     as the core Fig-3 DEG  ∪  hotspot / LIANA / SVG / conserved-core membership.
        #   genetic (coloc, 3,038) = all coloc (ABF + SuSiE, every GWAS) + fine-map
        #     credible sets (SuSiEx/MESuSiE) + regulatory-gwas-driven + TWAS + burden + ClinVar.
        # Precomputed gene lists (provenance in Analysis/Spatial/results/universe_validation/).
        UDIR = os.path.join(BASE, "Analysis/Spatial/results/universe_validation")
        rna_sig = set(open(os.path.join(UDIR, "universe_transcriptomic.txt")).read().split())
        coloc   = set(open(os.path.join(UDIR, "universe_genetic.txt")).read().split())
        target  = set(open(os.path.join(UDIR, "prioritized_universe_FINAL.txt")).read().split())
        coloc_only = coloc - rna_sig          # genetic-only axis (not a transcriptomic candidate)
        print(f"[compute_counts] BROADENED universe: transcriptomic {len(rna_sig)} "
              f"∪ genetic {len(coloc)} = {len(target)}")

        # ── VALIDATION = INCLUSIVE across ALL biological layers per modality ──
        # A gene is "validated" by a modality if it is significant in ANY of that
        # modality's disease/progression/co-localization layers — NOT just a single
        # disease-vs-control gate. This mirrors the multi-contrast prioritized universe.
        # Proteomics: every protein contrast across BOTH DIA-MS proteomes (liver
        # PXD051911 + PLASMA PXD052937): disease-vs-control + MASH-vs-MASL + NAS-high-
        # vs-low. (Olink is withdrawn 2026-06-01 — labels unrecoverable — so excluded.)
        c = pd.read_csv(PROT_CSV)
        P_meas_all = set(c["gene"].dropna())
        P_val_all  = set(c.loc[c["protein_padj"] < 0.05, "gene"].dropna())
        Pmeas = P_meas_all & target
        Pval  = P_val_all  & target

        a_all = pd.read_csv(ATLAS, low_memory=False).rename(columns={"human_symbol": "sym"})
        def _u(m):    return set(a_all.loc[m.fillna(False), "sym"])
        def _padj(col): return _u(a_all[col] < 0.05) if col in a_all.columns else set()
        def _flag(col): return _u((a_all[col] == True) | (a_all[col].astype(str).str.lower() == "true")) if col in a_all.columns else set()
        def _nn(col):   return _u(a_all[col].notna() & (a_all[col].astype(str).str.strip() != "") & (a_all[col].astype(str).str.lower() != "nan")) if col in a_all.columns else set()
        # Spatial: hepatocyte spatial localization + Govaere 2026 MASH zonation
        # signatures (centrilobular + periportal). SVG/Moran's-I is a TRANSCRIPTOMIC
        # prioritization candidate (excluded here — circular); GeoMx niche + zonation
        # yield 0 in-universe; CosMx = slide-level DIRECTION only (no valid cell-level p).
        S_val_all = (_padj("spatial_hep_wilcoxon_padj_bh")
                     | _padj("signature_govaere2026_epithelia_mash_centrilobular_padj")
                     | _padj("signature_govaere2026_epithelia_mash_periportal_padj"))
        # measured base MUST cover EVERY column the gate uses (hep + centrilobular +
        # periportal) — a periportal-only-validated gene missing from measured makes
        # S_val ⊄ S_meas, and the permutation null's rng.choice(meas, size=val_n)
        # raises "sample larger than population". `| S_val_all` guarantees the
        # superset invariant regardless of per-column NA/typing quirks.
        S_meas_all = set(a_all.loc[a_all["spatial_hep_wilcoxon_padj_bh"].notna()
                                   | a_all.get("signature_govaere2026_epithelia_mash_centrilobular_padj", pd.Series(index=a_all.index)).notna()
                                   | a_all.get("signature_govaere2026_epithelia_mash_periportal_padj", pd.Series(index=a_all.index)).notna(), "sym"]) | S_val_all
        # scATAC / epigenomic: hepatocyte differential accessibility + chromVAR TF
        # motif + SCENIC+ regulon target + promoter accessibility. (SCENIC+ enhancer-
        # gene link excluded: ~every gene has one -> trivially true, not validation.)
        A_val_all = (_padj("hepatocyte_da_padj") | _nn("chromvar_top_tf")
                     | _nn("scenic_regulon_tf") | _flag("human_promoter_accessible"))
        # same superset invariant as spatial: cover DA + chromVAR + SCENIC measured
        # bases, then `| A_val_all` folds in the promoter-accessible flagged genes
        # (bool column: notna == genome-wide, so we take the flagged set via A_val_all
        # rather than notna to keep the measured base — hence the cascade n.s.
        # remainder — sane instead of exploding to the whole atlas).
        A_meas_all = set(a_all.loc[a_all["hepatocyte_da_padj"].notna()
                                   | a_all.get("chromvar_top_tf", pd.Series(index=a_all.index)).notna()
                                   | a_all.get("scenic_regulon_tf", pd.Series(index=a_all.index)).notna(), "sym"]) | A_val_all
        Smeas = S_meas_all & target; Sval = S_val_all & target
        Ameas = A_meas_all & target; Aval = A_val_all & target

        # ── PER-DATASET proteomics splits (inclusive: union across a dataset's contrasts) ──
        def _proteo_ds(prefix):
            x = c[c["dataset"].str.startswith(prefix)]
            m = set(x["gene"].dropna()) & target
            v = set(x.loc[x["protein_padj"] < 0.05, "gene"].dropna()) & target
            return len(m), len(v)
        P_liver_meas,  P_liver_val  = _proteo_ds("PXD051911")   # liver DIA-MS (Boel)
        P_plasma_meas, P_plasma_val = _proteo_ds("PXD052937")   # PLASMA DIA-MS (Sourianarayanane)
        S_vis_meas = len(Smeas)
        _geo = [c0 for c0 in ("spatial_govaere2026_geomx_sh_vs_ls_padj",
                              "spatial_govaere2026_geomx_sh_vs_pt_padj") if c0 in a_all.columns]
        S_geomx_meas = len(set(a_all.loc[a_all[_geo].notna().any(axis=1), "sym"]) & target) if _geo else 0

        # ── PER-LAYER bundles for the VALIDATION-LAYER column (2 per modality) ──
        # Spatial splits into hepatocyte spatial-localization (Visium GSE192741) and
        # MASH zonation signatures (Govaere GeoMx CL+PP) — each with its OWN measured
        # base, so no node's significant count exceeds its measured base.
        S_hep_meas = len(set(a_all.loc[a_all["spatial_hep_wilcoxon_padj_bh"].notna(), "sym"]) & target)
        S_hep_val  = len(_padj("spatial_hep_wilcoxon_padj_bh") & target)
        _zc = [c0 for c0 in ("signature_govaere2026_epithelia_mash_centrilobular_padj",
                             "signature_govaere2026_epithelia_mash_periportal_padj") if c0 in a_all.columns]
        S_zon_meas = len(set(a_all.loc[a_all[_zc].notna().any(axis=1), "sym"]) & target) if _zc else 0
        S_zon_val  = len((_padj("signature_govaere2026_epithelia_mash_centrilobular_padj")
                          | _padj("signature_govaere2026_epithelia_mash_periportal_padj")) & target)
        # Govaere CosMx: NO valid cell-level p (E1b pseudoreplication) -> DIRECTION only.
        # A universe gene is CosMx-concordant if |MASH logFC| > 0.25 (same effect-size
        # floor as the DEG TREAT test) in ANY profiled cell type. Shown as a spatial
        # validator (direction, daggered), never a significance-gated count.
        _cx = [c0 for c0 in a_all.columns
               if c0.startswith("spatial_govaere2026_cosmx_") and c0.endswith("_logfc")]
        S_cos_meas = len(set(a_all.loc[a_all[_cx].notna().any(axis=1), "sym"]) & target) if _cx else 0
        S_cos_val  = len(set(a_all.loc[a_all[_cx].abs().gt(0.25).any(axis=1), "sym"]) & target) if _cx else 0
        # scATAC splits into disease chromatin ACCESSIBILITY (hepatocyte DA + promoter)
        # and regulatory TF ACTIVITY (chromVAR motif + SCENIC+ regulon).
        A_acc_meas = len(set(a_all.loc[a_all["hepatocyte_da_padj"].notna(), "sym"]) & target)
        A_acc_val  = len((_padj("hepatocyte_da_padj") | _flag("human_promoter_accessible")) & target)
        A_tf_val   = len((_nn("chromvar_top_tf") | _nn("scenic_regulon_tf")) & target)
        A_tf_meas  = max(A_tf_val, len(set(a_all.loc[a_all["chromvar_top_tf"].notna(), "sym"]) & target))

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
            # per-dataset / per-layer validators (validation-layer alluvial)
            P_liver_meas=P_liver_meas, P_liver_val=P_liver_val,
            P_plasma_meas=P_plasma_meas, P_plasma_val=P_plasma_val,
            S_vis_meas=S_vis_meas, S_geomx_meas=S_geomx_meas,
            S_hep_meas=S_hep_meas, S_hep_val=S_hep_val,
            S_zon_meas=S_zon_meas, S_zon_val=S_zon_val,
            S_cos_meas=S_cos_meas, S_cos_val=S_cos_val,
            A_acc_meas=A_acc_meas, A_acc_val=A_acc_val,
            A_tf_meas=A_tf_meas, A_tf_val=A_tf_val,
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
        return res
    except Exception as e:
        # NO stale fallback: a broken source must FAIL the render, never silently
        # produce a mislabeled panel with frozen numbers. Re-raise with context.
        import traceback
        traceback.print_exc()
        raise RuntimeError(f"[compute_counts] source unreadable — refusing to render "
                           f"stale/mislabeled numbers: {e}") from e


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
        print("[permutation_null] no gene sets — skipping null")
        return None
    import numpy as np
    # Sanitize gene sets for the null (a stray NaN symbol can slip into the
    # universe; drop non-string / empty so sorting & indexing are well-defined).
    # compute_counts() already reported the live COUNTS — this only guards
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
    # Drop non-string members (a COLOC locus may carry a pp4>0.5 NaN gene symbol)
    # so sorted() does not choke on mixed float/str. This only affects the
    # Monte-Carlo null's sampling pool; every DISPLAYED count comes from len() in
    # compute_counts (d[...]), so the 9,882 universe / 447 spatial / 141 ge2 /
    # 3 all3 values are unchanged.
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
    fig.savefig(p, bbox_inches="tight", pad_inches=0.02, facecolor="white")
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
      OUTCOME      TWO groups only: VALIDATED (>=1 assay, deep teal) vs UNVALIDATED
                   (grey; no >=1-assay support — includes genes never measured by any
                   assay, so "unvalidated" not "not significant") — every dataset's
                   validated flow merges into ONE Validated node, labelled with the
                   UNIQUE union count (assays overlap, so inflowing ribbons sum higher).

    Counts are live from compute_counts (so they track the active COLOC scope).
    Because the two DIA-MS proteomes (and the three spatial arms) share genes, the
    per-assay validated ribbons entering the Validated node sum to MORE than the
    d['ge1'] unique validated genes — that overlap is exactly why the outcome is
    collapsed to a single deduplicated Validated node rather than 6 disjoint grooves."""
    import matplotlib.colors as mcolors

    def tint(hexc, f):                       # lighten a modality hue for sub-nodes
        r, g, b = mcolors.to_rgb(hexc)
        return (r + (1 - r) * f, g + (1 - g) * f, b + (1 - b) * f)

    BLACK = "#231f20"  # 2026-07-09 house re-skin: near-black ink (was pure #000000)
    # DATASET nodes: (modality, dataset key [unique], platform/analysis, measured, validated).
    # EVERY dataset actually used to validate is a node here — no separate 'also profiled'
    # tier. Each has its OWN measured base so no node's validated count exceeds it. CosMx
    # is direction-only (no valid cell-level p, E1b); explained in the caption, not on-panel.
    # Olink (withdrawn 2026-06-01) and Vu 2025 (0 universe genes validated) are genuinely
    # not validators, so they simply do not appear.
    # DATASET column labels the ASSAY/MODALITY only — no GSE/PXD accession (dropped
    # from the panel; accessions live in the printed caption). 6th tuple field = the
    # on-panel display label. The atac modality is single-NUCLEUS ATAC (snATAC), so the
    # MODALITY lane is named "snATAC" too — NOT "scATAC" — to stay consistent with the
    # dataset (sc and sn are distinct assays; GSE244832 is snATAC).
    DS = [
        ("proteo",  "PXD051911", "liver DIA-MS",  d["P_liver_meas"],  d["P_liver_val"],  "liver DIA-MS"),
        ("proteo",  "PXD052937", "plasma DIA-MS", d["P_plasma_meas"], d["P_plasma_val"], "plasma DIA-MS"),
        ("spatial", "GSE192741", "Visium",        d["S_hep_meas"],    d["S_hep_val"],    "Visium"),
        ("spatial", "GeoMx",     "Govaere",       d["S_zon_meas"],    d["S_zon_val"],    "GeoMx"),
        ("spatial", "CosMx",     "Govaere",       d["S_cos_meas"],    d["S_cos_val"],    "CosMx"),
        ("atac",    "GSE244832", "snATAC",        d["A_meas"],        d["A_val"],        "snATAC"),
    ]
    lanes   = ["proteo", "spatial", "atac"]
    MODNAME = {"proteo": "proteomics", "spatial": "spatial", "atac": "snATAC"}
    COL     = {"proteo": C_PROTEO, "spatial": C_SPATIAL, "atac": C_ATAC}
    # modality measured = sum of its 2 layer-bundles' measured (fan conserves)
    MOD_MEAS = {l: float(len([x for x in DS if x[0] == l])) for l in lanes}  # inventory: room ∝ DATASET COUNT
    ds_of = {l: [x for x in DS if x[0] == l] for l in lanes}

    H, GAP, BW = 11.5, 0.85, 0.13
    # SEC widened (2.30->2.65, 2026-07-08) to use the room freed by wrapping the
    # left-most "transcriptomics"/"genetics" labels onto 2 lines (below) -- keeps
    # the same overall tight-cropped panel footprint but stretches the alluvial
    # itself into the reclaimed space instead of leaving it as empty left margin.
    SEC = 2.65
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
        sub = _stack([g[1] for g in grp], {g[1]: 1.0 for g in grp},
                     my1 - my0, GAP * 1.0, min(0.9, (my1 - my0) / len(grp) / 1.4))
        for k, (a, b) in sub.items():
            slotDS[k] = (my0 + a, my0 + b)

    # OUTCOME column: TWO groups only — Validated (any assay) vs Unvalidated.
    # The former 6 per-dataset validated + 6 n.s. endpoints collapse into ONE
    # Validated node (top) and ONE Not-significant node (bottom). Each dataset's
    # validated portion sweeps UP into the shared Validated node, its non-validated
    # portion DOWN into the shared Not-significant node; inside each group the ribbons
    # keep a per-dataset sub-slot (stacked in DS order) so they land without crossing.
    # The Validated node's LABEL is the UNIQUE union d['ge1'] (genes validated by >=1
    # assay); the inflowing per-assay validated ribbons sum higher than that because
    # assays share genes — the very overlap the 6-groove version double-counted.
    order = [g[1] for g in DS]
    VAL = {g[1]: g[4] for g in DS}; MEA = {g[1]: max(g[3], g[4]) for g in DS}
    dsh   = {k: (slotDS[k][1] - slotDS[k][0]) for k in order}     # dataset node heights
    vfrac = {k: (VAL[k] / MEA[k] if MEA[k] else 0.0) for k in order}
    val_w = {k: dsh[k] * vfrac[k] for k in order}                 # validated ribbon width
    ns_w  = {k: dsh[k] * (1.0 - vfrac[k]) for k in order}         # n.s. ribbon width
    GG = 0.9                                                      # gap between the 2 groups
    vtot = sum(val_w.values()); ntot = sum(ns_w.values())
    scO  = (H - GG) / (vtot + ntot) if (vtot + ntot) else 1.0
    subV = {}; y = 0.0
    for k in order:
        h = val_w[k] * scO; subV[k] = (y, y + h); y += h
    vy0, vy1 = 0.0, y
    ny0 = vy1 + GG; subN = {}; y = ny0
    for k in order:
        h = ns_w[k] * scO; subN[k] = (y, y + h); y += h
    ny1 = y

    # Compact canvas (Fig 2A idiom): the data content is drawn on a small figure so the
    # fixed 6 pt text reads large relative to the panel, instead of tiny on a 7x5.4 sheet.
    # Sized + tightly cropped to the Fig 4 layout slot (panel A, 3.03x2.13in; 2026-07-08).
    fig, ax = plt.subplots(figsize=(3.50, 2.71))
    ax.set_xlim(-1.3, X3 + 1.55); ax.set_ylim(-0.55, H + 0.15)
    ax.invert_yaxis(); ax.axis("off")

    def node(x, y0, y1, c, z=4, ec="white"):
        ax.add_patch(Rectangle((x - BW / 2, y0), BW, y1 - y0, facecolor=c,
                               edgecolor=ec, lw=0.2, zorder=z))

    # ── PRIORITIZED origin: RNA-seq signal (grey, TREAT-core inset) + GWAS-COLOC ──
    o_gap = 0.24; avail = H - o_gap; o_tot = d["substrate"] + d["coloc"]
    hR = avail * d["substrate"] / o_tot
    rna0, rna1 = 0.0, hR; gco0, gco1 = hR + o_gap, H
    node(X0, rna0, rna1, C_RNA)
    node(X0, gco0, gco1, C_GWAS)
    lx = X0 - BW / 2 - 0.16
    ax.text(lx, (rna0 + rna1) / 2, f"transcriptomics\n{d['substrate']:,}", ha="right", va="center", color=BLACK, zorder=8)
    ax.text(lx, (gco0 + gco1) / 2, f"genetics\n{d['coloc']:,}", ha="right", va="center", color=BLACK, zorder=8)

    # ── PRIORITIZED -> MODALITY (partition mouth by validator-coverage share) ──
    tot = float(sum(MOD_MEAS.values())); y = 0.0
    for l in lanes:
        h = H * MOD_MEAS[l] / tot
        _ribbon(ax, X0 + BW / 2, y, y + h, X1 - BW / 2, slotMod[l][0], slotMod[l][1], COL[l], alpha=0.30)
        y += h
    for l in lanes:
        y0, y1 = slotMod[l]; node(X1, y0, y1, COL[l])
        ax.text(X1 - BW / 2 - 0.10, (y0 + y1) / 2, MODNAME[l], ha="right", va="center", color=BLACK, zorder=8)

    # ── MODALITY -> DATASET (fan each modality into its named datasets) ──
    for l in lanes:
        m0, m1 = slotMod[l]; grp = ds_of[l]; yy = m0
        for g in grp:
            dd0, dd1 = slotDS[g[1]]; seg = dd1 - dd0
            _ribbon(ax, X1 + BW / 2, yy, yy + seg, X2 - BW / 2, dd0, dd1, tint(COL[l], 0.25), alpha=0.32)
            yy += seg
    _wbox = dict(boxstyle="round,pad=0.14", fc="white", ec="#D9D9D9", lw=0.3, alpha=0.9)
    for g in DS:
        dd0, dd1 = slotDS[g[1]]; node(X2, dd0, dd1, tint(COL[g[0]], 0.18))
        # single white-backed assay-name label ABOVE the node — lifts off the ribbons
        ax.text(X2, dd0 - 0.10, g[5], ha="center", va="bottom", color=BLACK,
                zorder=10, bbox=_wbox)

    # ── DATASET -> OUTCOME (each dataset's validated sweeps UP into the single
    #    Validated node, non-validated DOWN into the single Not-significant node) ──
    for g in DS:
        acc, l = g[1], g[0]; dd0, dd1 = slotDS[acc]
        hv = (dd1 - dd0) * vfrac[acc]
        a0, a1 = subV[acc]
        _ribbon(ax, X2 + BW / 2, dd0, dd0 + hv, X3 - BW / 2, a0, a1, COL[l], alpha=0.52)
        b0, b1 = subN[acc]
        _ribbon(ax, X2 + BW / 2, dd0 + hv, dd1, X3 - BW / 2, b0, b1, C_NS, alpha=0.45)
    # single Validated node (deep teal) + single Not-significant node (control grey)
    node(X3, vy0, vy1, C_CONV)
    node(X3, ny0, ny1, C_NS)
    ax.text(X3 + BW / 2 + 0.13, (vy0 + vy1) / 2, f"validated\n{d['ge1']:,}",
            ha="left", va="center", color=BLACK, zorder=8)
    # Unvalidated node left unlabelled (per PI 2026-07-07): the grey node speaks for itself.

    # Stage headers (Prioritized / Modality / Dataset / Outcome) removed per PI
    # (2026-07-07) — the stages read from the node labels + caption.

    _save(fig, "fig4a_overview_cascade.pdf")
    # printed caption (no on-panel subtitle, per house style)
    print(f"[caption:cascade] Fig 4 opener. PRIORITIZED universe = {d['universe']:,} genes = "
          f"transcriptomic candidates {d['substrate']:,} (all bulk + single-cell DE contrasts at the "
          f"same TREAT interval-null test, lfc=0.25, as the core DEG — disease-vs-control, MASH-vs-MASL, "
          f"stage, fibrosis-gradient — ∪ hotspot ∪ LIANA ∪ conserved-core) ∪ genetic candidates "
          f"{d['coloc']:,} (coloc [ABF + SuSiE, all GWAS] ∪ fine-map credible sets ∪ regulatory-GWAS ∪ "
          f"TWAS ∪ burden ∪ ClinVar); {d.get('coloc_shared',0):,} shared. Each prioritized gene is tested "
          "in EVERY dataset we hold, INCLUSIVELY across all biological layers (not just disease-vs-control): "
          "proteomics = disease-vs-control + MASH-vs-MASL + NAS-high-vs-low across liver (PXD051911) + "
          "plasma (PXD052937) DIA-MS; spatial = hepatocyte spatial-localization (Visium GSE192741) + MASH "
          "zonation signatures (Govaere GeoMx CL+PP) + CosMx MASH direction-concordance "
          "(|logFC|>0.25, no valid cell-level p — slide-level pseudoreplication, E1b); scATAC = hepatocyte "
          "differential accessibility + promoter accessibility + chromVAR TF-motif + SCENIC+ regulon "
          "(GSE244832). 'VALIDATED' = significant/positive/direction-concordant in that dataset. Per-modality "
          f"union: proteomics {d['P_val']}, spatial {d['S_val']} (Visium+GeoMx gated) + CosMx {d['S_cos_val']} "
          f"direction, scATAC {d['A_val']}; any-assay {d['ge1']} of {d['universe']:,}. Olink (Yang/Zeybel) "
          "withdrawn 2026-06-01 (unrecoverable labels) and Vu 2025 Visium (0 universe genes validated) are "
          "not validators and are omitted. OUTCOME collapses to two deduplicated groups — "
          f"VALIDATED = {d['ge1']:,} unique genes hit in >=1 assay (of {d['universe']:,}), UNVALIDATED = "
          f"{d['none']:,} (no >=1-assay support — includes genes not measured by any of the three assays); "
          f"because assays share genes the per-dataset validated ribbons entering the Validated "
          f"node sum to {sum(g[4] for g in DS):,} (> the {d['ge1']:,} unique union), so the node is labelled "
          "with the unique count.")


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
        print("[build_enrichment] no stats — skipping"); return
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
        undef = not np.isfinite(r["exp"])   # null undefined (background saturated: meas≈val)
        if not undef:
            x0, x1 = (XCH, r["fold"]) if r["fold"] >= XCH else (r["fold"], XCH)
            ax.add_patch(Rectangle((x0, yy - 0.32), x1 - x0, 0.64,
                                   facecolor=r["color"], edgecolor="none",
                                   alpha=0.95 if r["sig"] else 0.55, zorder=3))
        ax.text(-0.10, yy, r["label"], ha="right", va="center", color=INK, zorder=6)
        if undef:
            # meas == val -> hypergeometric expectation undefined; report count only
            ax.text(XANN, yy, f"{r['obs']:,}  null n/a", ha="left", va="center",
                    color=INK, zorder=6)
        else:
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

    _save(fig, "figS4f.pdf")   # demoted from main Fig 4a to Supp Fig S4F (2026-07-07)
    # caption fully live (no hard-coded counts): the exclusive P∩S bar = PS - all3,
    # matching both the panel bar and the permutation stat (which is computed on the
    # exclusive intersection); all-3 is reported live with its own chance expectation.
    ps_excl = d["PS"] - d["all3"]
    print("[caption:upset] UpSet of the three orthogonal validation lenses on the "
          f"{d['universe']:,}-gene prioritized universe. Bars = exclusive intersection sizes; "
          "dashed tick = permutation-null expected (seed 42, 10,000 draws). Only "
          "proteomics∩spatial clears its chance marker "
          f"({ps_excl} vs {stats['conv']['PS']['exp']:.1f}, "
          f"{stats['conv']['PS']['fold']:.1f}x, p={stats['conv']['PS']['p']:.1g}); "
          f"all-3 = {d['all3']} (chance {stats['conv']['all3']['exp']:.1f}, "
          f"p={stats['conv']['all3']['p']:.2g}).")


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
    # build_enrichment intentionally NOT rendered (dropped 2026-07-07 per PI):
    # fig4a_overview_enrichment.pdf is no longer generated.
    build_upset(d, stats)
    print("Wrote fig4a_overview_{cascade}.pdf (+ figS4f upset) to", OUT_DIR)


if __name__ == "__main__":
    main()
