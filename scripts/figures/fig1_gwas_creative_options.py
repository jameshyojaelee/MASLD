#!/usr/bin/env python3
"""
Fig 2A GWAS Portfolio — alluvial, evidence cascade, and per-study bars (breadth, no N encoding).
KEY MESSAGE: 35 Tier-1/2 (liver-specific) GWAS spanning 5 ancestries (EUR/AFR/AMR/EAS/SAS,
MVP NAFLD/ALT/AST strata included) x direct liver-disease/enzyme traits — coverage breadth.
Scoped 2026-07-06 to the Tier-1/2 MAIN portfolio ONLY (placement=="main" in
GWAS/finemapping/config/gwas_trait_tier.tsv): the Tier-3/4 supp strata (MVP Cirrhosis /
Chronic-liver-disease / Albumin / Platelet) are DROPPED here and moved to a supplementary
figure. All counts are DATA-DRIVEN at runtime from the GWAS registry, per-study fine-mapping
summary, and the per-gene-per-GWAS colocalisation table (no hardcoded manifest).

Canonical Fig 2A panel outputs (figures/main/fig2_genetics/panels/):
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
# The portfolio (35 Tier-1/2 liver-specific GWAS across 5 ancestries, incl. the MVP /
# Million Veteran Program NAFLD/ALT/AST strata) is NOT hardcoded here — every count is
# recomputed at runtime from the canonical GWAS files below, RESTRICTED to the MAIN
# (Tier-1/2) strata via the placement=="main" allowlist in gwas_trait_tier.tsv.
import csv as _csv
import json as _json
import math as _math


def _project_root():
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def _out_path(fname):
    d = os.path.join(_project_root(), "figures", "misc")
    os.makedirs(d, exist_ok=True)
    return os.path.join(d, fname)


_ROOT          = os.environ.get("MASLD_PROJECT_ROOT", _project_root())
_REGISTRY_TSV  = os.path.join(_ROOT, "GWAS/finemapping/config/gwas_registry.tsv")
_TIER_TSV      = os.path.join(_ROOT, "GWAS/finemapping/config/gwas_trait_tier.tsv")
_PHENOTYPE_REGISTRY = os.path.join(
    _ROOT, "Analysis/Multimodal_Program_Projection/candidates",
    "program-context-v2-candidate-2026-08-07/genetics_context/phenotype_registry.tsv",
)
_STUDY_SUMMARY = os.path.join(_ROOT, "GWAS/finemapping/results/study_summary.csv")
# Per-gene-per-GWAS long COLOC table. Override only with an audited, promoted release.
# MAIN strata directly from here: the pre-aggregated wide gene_level_coloc.csv and even the
# gene_level_coloc_tier12.csv `driving_gwas` are FULL-portfolio drivers (a MAIN-union gene
# can still be driven by a Tier-3/4 stratum), which would reintroduce the supp traits.
_GENE_COLOC    = os.environ.get(
    "FIG2_COLOC_INPUT",
    os.path.join(_ROOT, "GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv"),
)
# Per-variant credible sets: source of the UNIQUE physical fine-mapped loci (Change 1). We merge
# per-(study,locus) credible-set windows across ALL MAIN strata so the same physical locus tagged
# by several correlated enzymes/ancestries is counted ONCE (the old cascade summed per-study
# n_loci, ~2,824, which double-counts).
_CREDIBLE_SETS = os.path.join(_ROOT, "GWAS/finemapping/results/credible_sets.csv")
# r2_summary.json of the COLOC release (11_assemble_coloc_r2.py). Its count of
# single-signal-only genes whose SuSiE test was untestable goes into the caption.
_COLOC_SUMMARY = os.environ.get("FIG2_COLOC_SUMMARY", "")
_UNTESTABLE    = "susie_untestable_insufficient_shared_posterior"
_LOCI_PAD      = 250_000   # +/- window padding (bp) before merging overlapping loci
_FIG2_SIZE_TSV = os.path.join(_ROOT, "figures/layout_specs/figure2_panel_sizes.tsv")
_CASCADE_REL   = "main/fig2_genetics/panels/fig2A_gwas_cascade.pdf"
_CASCADE_EXPECTED_SIZE = (3.53, 2.35)


def _contract_size(pdf_rel):
    """Return the exact Illustrator placement size from the Figure 2 contract."""
    with open(_FIG2_SIZE_TSV) as fh:
        rows = {r["pdf"].strip(): r for r in _csv.DictReader(fh, delimiter="\t")}
    if pdf_rel not in rows:
        raise ValueError(f"No size-contract row for {pdf_rel}: {_FIG2_SIZE_TSV}")
    size = float(rows[pdf_rel]["width_in"]), float(rows[pdf_rel]["height_in"])
    if pdf_rel == _CASCADE_REL and any(
            abs(got - expected) > 1e-12
            for got, expected in zip(size, _CASCADE_EXPECTED_SIZE)):
        raise ValueError("Figure 2A size contract must be 3.53 x 2.35 in before rendering")
    return size

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


def _merge_padded(intervals, pad):
    """Merge a list of (min_pos, max_pos, pip, ancestry) windows padded by +/-`pad`.
    Returns one (best_pip, ancestry) per merged interval — the strongest-PIP contributor."""
    lst = sorted((mn - pad, mx + pad, pip, a) for mn, mx, pip, a in intervals)
    merged = []   # each: [running_padded_end, best_pip, ancestry]
    for ps, pe, pip, a in lst:
        if merged and ps <= merged[-1][0]:
            if pe > merged[-1][0]:
                merged[-1][0] = pe            # extend the running padded end
            if pip > merged[-1][1]:           # keep the strongest-PIP contributor's ancestry
                merged[-1][1], merged[-1][2] = pip, a
        else:
            merged.append([pe, pip, a])
    return [(m[1], m[2]) for m in merged]


def _unique_loci_by_at(MAIN, anc_of_study, pad=_LOCI_PAD):
    """Fine-mapped loci from the per-variant credible sets (Change 1 + within-trait refinement).

    For each (study, locus) build the genomic window [min(position), max(position)] spanned by its
    credible-set variants (MAIN strata only), then merge overlapping windows padded by +/-`pad`.
    Returns (at_loci_within, n_unique_global):
      * at_loci_within[(ancestry, trait)] = WITHIN-TRAIT-unique loci — windows are merged within
        each (trait, chromosome) across that trait's studies/ancestries, and each merged locus is
        assigned to its strongest-PIP ancestry. A physical locus tagged by several enzymes therefore
        RECURS across those traits, so the per-(ancestry, trait) counts sum to ~360 (NAFLD 13 /
        NASH 4 / PDFF 13 / ALT 122 / AST 115 / GGT 93). This is the illustrative per-trait funnel.
      * n_unique_global = distinct physical loci after merging across ALL traits per chromosome
        (global de-dup) = ~265 — the summary count (some of the 360 per-trait loci are the same
        physical locus shared across traits)."""
    win = {}   # (study, locus) -> [chrom, min_pos, max_pos, best_pip, ancestry, trait]
    with open(_CREDIBLE_SETS) as fh:
        for r in _csv.DictReader(fh):
            s = r["study"]
            if s not in MAIN:
                continue
            pos = int(float(r["position"]))
            pip = _fnum(r.get("recommended_pip"))
            if pip != pip:
                pip = _fnum(r.get("max_pip"))
            if pip != pip:
                pip = 0.0
            key = (s, r["locus"])
            w = win.get(key)
            if w is None:
                win[key] = [r["chromosome"], pos, pos, pip, anc_of_study.get(s), _trait_of(s)]
            else:
                if pos < w[1]:
                    w[1] = pos
                if pos > w[2]:
                    w[2] = pos
                if pip > w[3]:
                    w[3] = pip
    # (1) WITHIN-TRAIT-unique per (ancestry, trait) — merge windows within each (trait, chrom).
    by_trait_chr = defaultdict(list)
    for chrom, mn, mx, pip, a, t in win.values():
        by_trait_chr[(t, chrom)].append((mn, mx, pip, a))
    at_within = defaultdict(int)
    for (t, chrom), ivs in by_trait_chr.items():
        for _pip, a in _merge_padded(ivs, pad):
            at_within[(a, t)] += 1
    # (2) GLOBAL-unique across ALL traits per chromosome (summary count only).
    by_chr = defaultdict(list)
    for chrom, mn, mx, pip, a, t in win.values():
        by_chr[chrom].append((mn, mx, pip, a))
    n_unique_global = sum(len(_merge_padded(ivs, pad)) for ivs in by_chr.values())
    return at_within, n_unique_global


def _load_manifest():
    """Read the Tier-1/2 (MAIN) GWAS registry, per-study fine-mapping yield, and per-gene
    coloc. Returns data-driven aggregates used by every panel — no frozen literals.
    RESTRICTED to the MAIN (Tier-1/2, liver-specific) strata: placement=="main" in the tier
    map. The Tier-3/4 supp strata (MVP Cirrhosis/ChronLiver/Albumin/Platelet) are excluded."""
    # 0) MAIN (Tier-1/2) allowlist from the tier map — the single filter that scopes the panel.
    MAIN = {r["study_name"] for r in _csv.DictReader(open(_TIER_TSV), delimiter="\t")
            if r["placement"] == "main"}
    phenotype_rows = [
        r for r in _csv.DictReader(open(_PHENOTYPE_REGISTRY), delimiter="\t")
        if r["placement"] == "main"
    ]
    if {r["study_name"] for r in phenotype_rows} != MAIN:
        raise ValueError("Phenotype registry MAIN studies disagree with gwas_trait_tier.tsv")
    phenotype_of_study = {r["study_name"]: r["phenotype_stratum"] for r in phenotype_rows}
    eqtl_status_of_study = {
        r["study_name"]: r["regulatory_ancestry_status"] for r in phenotype_rows
    }
    trait_classes = defaultdict(set)
    for r in phenotype_rows:
        trait_classes[r["trait"]].add(r["phenotype_stratum"])
    if any(len(classes) != 1 for classes in trait_classes.values()):
        raise ValueError("A displayed GWAS trait maps to multiple phenotype classes")
    trait_phenotype_class = {trait: next(iter(classes))
                             for trait, classes in trait_classes.items()}

    # 1) registry -> per-(ancestry, trait) GWAS-dataset counts (MAIN strata only)
    reg = [r for r in _csv.DictReader(open(_REGISTRY_TSV), delimiter="\t")
           if r["study_name"] in MAIN]
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

    # 2) study_summary -> per-study & per-(ancestry, trait) fine-mapped loci.
    # anc_of_study is now MAIN-only, so non-MAIN studies are skipped automatically.
    at_loci, study_loci = defaultdict(int), {}
    for r in _csv.DictReader(open(_STUDY_SUMMARY)):
        s = r["study"]
        a = anc_of_study.get(s)
        if a is None:
            continue
        nl = int(r["n_loci"])
        at_loci[(a, _trait_of(s))] += nl
        study_loci[s] = nl

    # 3) per-gene MAIN-restricted colocalisation, recomputed from the long per-gene-per-GWAS
    # table. Multi-signal genes are assigned to their strongest SuSiE result;
    # single-signal-only genes are assigned to their strongest ABF result.
    g_bs, g_ba = defaultdict(float), defaultdict(float)
    g_arg_susie, g_arg_abf = {}, {}
    g_untestable = set()   # a MAIN study where every signal pair failed the shared-posterior check
    for r in _csv.DictReader(open(_GENE_COLOC)):
        gw = r["gwas_name"]
        if gw not in MAIN:
            continue
        gene = r["gene"]
        if r.get("method") == _UNTESTABLE:
            g_untestable.add(gene)
        su, abf = _fnum(r["PP.H4.susie"]), _fnum(r["PP.H4.abf"])
        if su == su and su > g_bs[gene]:
            g_bs[gene] = su
        if su == su:
            prev = g_arg_susie.get(gene)
            if prev is None or su > prev[0] or (su == prev[0] and gw < prev[2]):
                g_arg_susie[gene] = (su, anc_of_study.get(gw), gw)
        if abf == abf and abf > g_ba[gene]:
            g_ba[gene] = abf
        if abf == abf:
            prev = g_arg_abf.get(gene)
            if prev is None or abf > prev[0] or (abf == prev[0] and gw < prev[2]):
                g_arg_abf[gene] = (abf, anc_of_study.get(gw), gw)

    # at_genes = per-(ancestry,trait) union count; at_genes_susie / at_genes_abf split the SAME
    # genes into the multi-signal COLOC block (PP.H4>0.5) versus the single-signal
    # COLOC-only block (ABF PP.H4>0.5 without multi-signal support).
    at_genes, anc_genes, study_genes = defaultdict(int), defaultdict(int), defaultdict(int)
    at_genes_susie, at_genes_abf = defaultdict(int), defaultdict(int)
    at_genes_susie_band, at_genes_abf_band = defaultdict(int), defaultdict(int)
    union_total = susie_total = abf_only_total = 0
    union_gene_state = {}   # gene -> (signal model drawn, SuSiE state)
    susie_b9 = susie_b89 = susie_b58 = 0   # SuSiE PP.H4 confidence bands: >0.9 / 0.8-0.9 / 0.5-0.8
    abf_b9 = abf_b89 = abf_b58 = 0         # ABF-only PP.H4 confidence bands (same cut-points)
    for gene in set(g_bs) | set(g_ba):
        # Rows whose gene symbol is blank collapse into a single empty-string key
        # and are not a named gene. Exclude them before deriving every release-specific
        # endpoint count; no historical union literal belongs in this data-driven path.
        if not gene.strip():
            continue
        su, abf = g_bs.get(gene, 0.0), g_ba.get(gene, 0.0)
        su_ok, abf_ok = su > 0.5, abf > 0.5
        if su_ok:
            susie_total += 1
            if su > 0.9:
                susie_b9 += 1
            elif su >= 0.8:
                susie_b89 += 1
            else:
                susie_b58 += 1
        if not (su_ok or abf_ok):
            continue
        union_total += 1
        # An untestable SuSiE test is not a multi-signal negative; the gene is still
        # drawn as single-signal only when its ABF evidence stands.
        union_gene_state[gene] = (
            "multi_signal" if su_ok else "single_signal_only",
            "positive" if su_ok else "untestable" if gene in g_untestable
            else "tested_le_0.5" if gene in g_arg_susie else "no_susie_result")
        _, a, gw = (g_arg_susie[gene] if su_ok else g_arg_abf[gene])
        t = _trait_of(gw)
        at_genes[(a, t)] += 1
        anc_genes[a] += 1
        study_genes[gw] += 1
        if su_ok:
            at_genes_susie[(a, t)] += 1
            band = ">0.9" if su > 0.9 else ("0.8-0.9" if su >= 0.8 else "0.5-0.8")
            at_genes_susie_band[(a, t, band)] += 1
        else:
            abf_only_total += 1
            at_genes_abf[(a, t)] += 1
            band = ">0.9" if abf > 0.9 else ("0.8-0.9" if abf >= 0.8 else "0.5-0.8")
            at_genes_abf_band[(a, t, band)] += 1
            if abf > 0.9:              # ABF-only tail is bottom-heavy (most in the weakest bin)
                abf_b9 += 1
            elif abf >= 0.8:
                abf_b89 += 1
            else:
                abf_b58 += 1

    # Fine-mapped loci (Change 1) — within-trait-unique per (ancestry,trait) [~360, the per-trait
    # funnel] plus the global de-duplicated unique count [~265, the summary label].
    at_loci_within, n_unique_loci = _unique_loci_by_at(MAIN, anc_of_study)

    return dict(reg=reg, anc_of_study=anc_of_study,
                phenotype_of_study=phenotype_of_study,
                eqtl_status_of_study=eqtl_status_of_study,
                trait_phenotype_class=trait_phenotype_class,
                at_gwas=at_gwas, anc_gwas=anc_gwas,
                n_gwas=len(reg), n_studies=len(studies), at_loci=at_loci, study_loci=study_loci,
                at_loci_within=at_loci_within, n_unique_loci=n_unique_loci,
                at_genes=at_genes, anc_genes=anc_genes, study_genes=study_genes,
                at_genes_susie=at_genes_susie, at_genes_abf=at_genes_abf,
                at_genes_susie_band=at_genes_susie_band,
                at_genes_abf_band=at_genes_abf_band,
                union_total=union_total, susie_total=susie_total, abf_only_total=abf_only_total,
                union_gene_state=union_gene_state,
                susie_bands=(susie_b9, susie_b89, susie_b58),
                abf_bands=(abf_b9, abf_b89, abf_b58))


_M          = _load_manifest()
N_GWAS      = _M["n_gwas"]        # 35 Tier-1/2 (liver-specific) GWAS datasets
N_STUDIES   = _M["n_studies"]     # distinct studies (Sveinbjornsson arms collapsed)
UNION_GENES = _M["union_total"]      # promoted Tier-1/2 union at PP.H4 > 0.5
SUSIE_GENES = _M["susie_total"]      # multi-signal COLOC genes at PP.H4 > 0.5
ABF_ONLY_GENES = _M["abf_only_total"]  # single-signal COLOC-only genes at PP.H4 > 0.5
N_UNIQUE_LOCI  = _M["n_unique_loci"]   # ~265 globally-unique physical loci (per-trait blocks sum to ~360)

# Multi-signal endpoint = sequential blue; single-signal-only endpoint = hatched grey.
# Single shared sequential PP.H4 palette (dark->light = >0.9 / 0.8-0.9 / 0.5-0.8), used for BOTH the
# SuSiE-COLOC and ABF-only coloc-gene blocks — the two are separated by label + a light hatch on ABF,
# NOT by hue. Side effect: SuSiE reads dark-dominated (confident), ABF light-dominated (weak).
PP4_TINTS    = ["#08519C", "#4292C6", "#9ECAE1"]   # PP.H4 >0.9 / 0.8-0.9 / 0.5-0.8
SUSIE_RIBBON = "#6BAED6"   # loci -> SuSiE ribbons (single hue, matches the block)
ABF_COLOR    = "#98A2B3"   # loci -> single-signal-only ribbons

# Ancestry order = descending GWAS count (EUR, AFR, EAS, AMR, SAS).
ANCESTRIES     = sorted(_M["anc_gwas"], key=lambda a: (-_M["anc_gwas"][a], a))
ANCESTRY_LABEL = {"EUR": "European", "AFR": "African", "AMR": "Admixed American",
                  "EAS": "East Asian", "SAS": "South Asian"}
ANCESTRY_N     = {a: _M["anc_gwas"][a] for a in ANCESTRIES}
# Canonical ancestry palette shared across the Fig 2 genetics panels (colorblind-safe).
ANCESTRY_COLOR = {"EUR": "#4C72B0", "AFR": "#55A868", "AMR": "#DD8452",
                  "EAS": "#C44E52", "SAS": "#8172B3"}

# Traits present in the Tier-1/2 (MAIN) portfolio, grouped disease -> imaging -> enzymes.
# The Tier-3/4 supp traits (Cirrhosis / Chronic liver disease / Albumin / Platelet) are
# NOT shown here — they move to the supplementary full-portfolio figure.
TRAITS      = ["NAFLD", "NASH", "PDFF", "ALT", "AST", "GGT"]
TRAIT_SHORT = {}
PHENOTYPE_CLASSES = [
    ("direct_masld_mash_diagnosis", "Direct diagnosis"),
    ("mri_pdff_or_histologic_steatosis", "MRI-PDFF /\nsteatosis"),
    ("alt_ast_or_ggt", "Liver enzymes"),
]
PHENOTYPE_CLASS_OF = _M["trait_phenotype_class"]
# Trait group: (label, col_start, col_end_excl, band_color) — greyscale (colour = ancestry).
TRAIT_GROUPS = [
    ("Disease diagnoses", 0, 2, "#444444"),
    ("Imaging",           2, 3, "#7F7F7F"),
    ("Liver enzymes",     3, 6, "#C0C0C0"),
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


def _write_cascade_source(path):
    """Write every count displayed in the Figure 2A cascade in long form."""
    fields = [
        "record_type", "study", "ancestry", "eqtl_ancestry_status", "phenotype_class",
        "trait", "signal_model", "pp_h4_band", "metric", "value", "source",
        "gene", "susie_state",
    ]
    rows = []
    class_eqtl_counts = defaultdict(int)
    for registry_row in sorted(_M["reg"], key=lambda row: row["study_name"]):
        study = registry_row["study_name"]
        phenotype_class = _M["phenotype_of_study"][study]
        eqtl_status = _M["eqtl_status_of_study"][study]
        class_eqtl_counts[(phenotype_class, eqtl_status)] += 1
        rows.append({
            "record_type": "study_registry", "study": study,
            "ancestry": registry_row["ancestry"],
            "eqtl_ancestry_status": eqtl_status,
            "phenotype_class": phenotype_class,
            "trait": _trait_of(study), "signal_model": "", "pp_h4_band": "",
            "metric": "gwas_stratum", "value": 1, "source": _PHENOTYPE_REGISTRY,
        })

    for phenotype_class, _short in PHENOTYPE_CLASSES:
        for eqtl_status in ("ancestry_matched_eur", "cross_ancestry_eqtl_limited"):
            rows.append({
                "record_type": "phenotype_class_eqtl_summary", "study": "",
                "ancestry": "", "eqtl_ancestry_status": eqtl_status,
                "phenotype_class": phenotype_class, "trait": "",
                "signal_model": "", "pp_h4_band": "", "metric": "gwas_strata",
                "value": class_eqtl_counts[(phenotype_class, eqtl_status)],
                "source": _PHENOTYPE_REGISTRY,
            })

    for ancestry in ANCESTRIES:
        eqtl_status = ("ancestry_matched_eur" if ancestry == "EUR"
                       else "cross_ancestry_eqtl_limited")
        for trait in TRAITS:
            phenotype_class = PHENOTYPE_CLASS_OF[trait]
            metrics = (
                ("gwas_strata", _M["at_gwas"].get((ancestry, trait), 0), _REGISTRY_TSV),
                ("fine_mapped_loci_within_trait",
                 _M["at_loci_within"].get((ancestry, trait), 0), _CREDIBLE_SETS),
                ("multi_signal_coloc_genes",
                 _M["at_genes_susie"].get((ancestry, trait), 0), _GENE_COLOC),
                ("single_signal_coloc_only_genes",
                 _M["at_genes_abf"].get((ancestry, trait), 0), _GENE_COLOC),
            )
            for metric, value, source in metrics:
                rows.append({
                    "record_type": "ancestry_trait",
                    "study": "",
                    "ancestry": ancestry,
                    "eqtl_ancestry_status": eqtl_status,
                    "phenotype_class": phenotype_class,
                    "trait": trait,
                    "signal_model": "",
                    "pp_h4_band": "",
                    "metric": metric,
                    "value": value,
                    "source": source,
                })
            for model, counts in (("multi_signal", _M["at_genes_susie_band"]),
                                  ("single_signal_only", _M["at_genes_abf_band"])):
                for band in (">0.9", "0.8-0.9", "0.5-0.8"):
                    rows.append({
                        "record_type": "ancestry_trait_posterior_band", "study": "",
                        "ancestry": ancestry, "eqtl_ancestry_status": eqtl_status,
                        "phenotype_class": phenotype_class, "trait": trait,
                        "signal_model": model, "pp_h4_band": band,
                        "metric": "coloc_genes", "value": counts.get((ancestry, trait, band), 0),
                        "source": _GENE_COLOC,
                    })

    for model, counts in (("multi_signal", _M["susie_bands"]),
                          ("single_signal_only", _M["abf_bands"])):
        for band, value in zip((">0.9", "0.8-0.9", "0.5-0.8"), counts):
            rows.append({
                "record_type": "posterior_band", "study": "", "ancestry": "",
                "eqtl_ancestry_status": "", "phenotype_class": "", "trait": "",
                "signal_model": model, "pp_h4_band": band,
                "metric": "coloc_genes", "value": value, "source": _GENE_COLOC,
            })

    for metric, value, source in (
        ("gwas_strata", N_GWAS, _REGISTRY_TSV),
        ("globally_unique_fine_mapped_loci", N_UNIQUE_LOCI, _CREDIBLE_SETS),
        ("multi_or_single_signal_coloc_gene_union", UNION_GENES, _GENE_COLOC),
    ):
        rows.append({
            "record_type": "global_summary", "study": "", "ancestry": "",
            "eqtl_ancestry_status": "", "phenotype_class": "", "trait": "",
            "signal_model": "", "pp_h4_band": "", "metric": metric,
            "value": value, "source": source,
        })

    # One row per union gene with its SuSiE state, so a single-signal-only gene
    # whose multi-signal test was untestable is visible as such.
    state_counts = defaultdict(int)
    for gene, (model, state) in sorted(_M["union_gene_state"].items()):
        state_counts[(model, state)] += 1
        rows.append({
            "record_type": "union_gene", "study": "", "ancestry": "",
            "eqtl_ancestry_status": "", "phenotype_class": "", "trait": "",
            "signal_model": model, "pp_h4_band": "", "metric": "coloc_gene",
            "value": 1, "source": _GENE_COLOC, "gene": gene, "susie_state": state,
        })
    for (model, state), value in sorted(state_counts.items()):
        rows.append({
            "record_type": "susie_state_summary", "study": "", "ancestry": "",
            "eqtl_ancestry_status": "", "phenotype_class": "", "trait": "",
            "signal_model": model, "pp_h4_band": "", "metric": "coloc_genes",
            "value": value, "source": _GENE_COLOC, "gene": "", "susie_state": state,
        })

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as handle:
        writer = _csv.DictWriter(handle, fieldnames=fields, delimiter="\t", restval="")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Saved: {path}")


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
# Extended cascade: ancestry -> GWAS trait -> fine-mapped loci -> COLOC genes
# ─────────────────────────────────────────────────────────────────────────────
def build_cascade(path, extra_paths=None):
    """Four-stage cascade: ancestry -> GWAS trait -> loci -> COLOC genes.

    Stages 1-3 are ancestry-coloured; the loci column shows the
    WITHIN-TRAIT-unique fine-mapped loci (merged credible-set windows within each trait; Change 1),
    NOT the old ~2,824 per-study n_loci sum. A physical locus tagged by several correlated enzymes
    recurs across those traits, so the per-trait blocks sum to ~360 (illustrative funnel); the
    globally de-duplicated unique count is ~265 (annotated). The genes ENDPOINT is restructured
    (Change 2) into the multi-signal COLOC block, shaded by PP.H4 confidence,
    plus a hatched single-signal COLOC-only block with the same PP.H4 bands.
    Loci -> gene ribbons route into the multi-signal (blue) versus single-signal-only
    (grey, faint) groups. FULLY DATA-DRIVEN from the GWAS registry, per-variant credible
    sets, and the per-gene-per-GWAS coloc table (each gene assigned to the ancestry/trait of its
    single strongest colocalisation)."""
    r_order = list(TRAITS)
    at_g, at_l, at_ge = _M["at_gwas"], _M["at_loci_within"], _M["at_genes"]
    at_gs, at_ga = _M["at_genes_susie"], _M["at_genes_abf"]
    at_gsb, at_gab = _M["at_genes_susie_band"], _M["at_genes_abf_band"]
    DAT = {}
    for a in ANCESTRIES:
        for t in r_order:
            g, l, ge = at_g.get((a, t), 0), at_l.get((a, t), 0), at_ge.get((a, t), 0)
            if g or l or ge:
                DAT[(a, t)] = (g, l, ge)
    MIDX = {"trait": 0, "loci": 1}   # genes endpoint uses its own method-split layout below
    H, GAP, BW = 13.0, 0.6, 0.06   # compact canvas; GAP widened so per-segment count labels clear each other
    # Equal horizontal width for every section (ancestry@0 -> trait -> loci -> genes).
    # The 1.0-unit spacing makes the four-stage data region occupy at least 2.40 in
    # in the 3.53-in contracted panel after the compact left/right furniture.
    SEC = 1.00
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
    L = {s: layout(MIDX[s]) for s in ("trait", "loci")}
    tot = {s: sum(v[MIDX[s]] for v in DAT.values()) for s in ("trait", "loci")}

    # Genes endpoint: two signal-model blocks normalised so the union spans the same H.
    xg = XCOL["genes"]
    METHOD_GAP = GAP
    gscale = (H - METHOD_GAP) / UNION_GENES
    susie_y0, susie_y1 = 0.0, SUSIE_GENES * gscale
    abf_y0,   abf_y1   = susie_y1 + METHOD_GAP, susie_y1 + METHOD_GAP + ABF_ONLY_GENES * gscale

    candidate_mode = bool(os.environ.get("FIG2_CANDIDATE_DIR"))
    figure_size = _contract_size(_CASCADE_REL)
    fig, ax = plt.subplots(figsize=figure_size)
    # Reserve 10 pt top / 8 pt left for the Illustrator panel letter and use
    # the requested compact 3 pt right/bottom margins. Everything else is
    # inside the fixed page box; never recover whitespace with bbox_inches.
    fig.subplots_adjust(
        left=(8 / 72) / figure_size[0],
        right=1 - (3 / 72) / figure_size[0],
        bottom=(3 / 72) / figure_size[1],
        top=1 - (10 / 72) / figure_size[1],
    )
    x_left = -0.43
    x_right = 3 * SEC + 0.75
    ax.set_xlim(x_left, x_right); ax.set_ylim(-0.6, H + 1.5); ax.invert_yaxis(); ax.axis("off")

    # Ancestry column (x=0)
    anc_tot = {a: sum(DAT[k][0] for k in DAT if k[0] == a) for a in ANCESTRIES}
    s0 = (H - GAP*(len(ANCESTRIES)-1)) / sum(anc_tot.values())
    anc_y = {}; y = 0.0
    for a in ANCESTRIES:
        h = anc_tot[a]*s0; anc_y[a] = (y, y+h)
        ax.add_patch(mpatches.Rectangle((-BW/2, y), BW, h, facecolor=ANCESTRY_COLOR[a], edgecolor="none", zorder=4))
        ax.text(-BW/2-0.03, y+h/2, f"{a} ({anc_tot[a]})", ha="right", va="center", fontsize=6, color="black")
        y += h + GAP
    ax.text(0, H + 0.12, f"{len(ANCESTRIES)} ancestries", ha="center", va="top",
            fontsize=6, color="black")

    # stage blocks + headers for the ancestry-coloured trait & loci columns (angled titles below)
    for s in ("trait", "loci"):
        seg = L[s][0]; x = XCOL[s]
        for (a, t), (y0, y1) in seg.items():
            ax.add_patch(mpatches.Rectangle((x-BW/2, y0), BW, y1-y0, facecolor=ANCESTRY_COLOR[a], edgecolor="white", lw=0.2, zorder=4))
        hdr = {"trait": f"{tot['trait']} GWAS", "loci": "Fine-mapped\nloci"}[s]
        ax.text(x, H + 0.12, hdr, ha="center", va="top",
                fontsize=6, color="black", linespacing=0.9)
    ax.text(xg, H + 0.12, "COLOC\ngenes", ha="center", va="top",
            fontsize=6, color="black", linespacing=0.9)
    # trait name + dataset count ABOVE each trait block, e.g. "NAFLD (11)"
    for t in r_order:
        y0, y1 = L["trait"][1][t]
        n_gw = sum(DAT[(a, t)][MIDX["trait"]] for a in ANCESTRIES if (a, t) in DAT)
        if n_gw == 0:
            continue
        ax.text(XCOL["trait"], y0 - 0.06, f"{TRAIT_SHORT.get(t, t)} ({n_gw})",
                ha="center", va="bottom", fontsize=6, color="black", zorder=8)
    # per-trait WITHIN-TRAIT-unique loci count ABOVE each loci block (blocks sum to ~360, since a
    # physical locus tagged by several enzymes recurs across those traits)
    for t in r_order:
        y0, y1 = L["loci"][1][t]
        cnt = sum(DAT[(a, t)][MIDX["loci"]] for a in ANCESTRIES if (a, t) in DAT)
        if cnt == 0:
            continue
        ax.text(XCOL["loci"], y0 - 0.06, str(cnt), ha="center", va="bottom", fontsize=6, color="black", zorder=8)

    # ── Genes endpoint blocks (Change 2) ────────────────────────────────────────
    # SuSiE-COLOC primary block: 3 PP.H4 tints (dark->light) = one-hue confidence gradient. Every
    # band carries its explicit PP.H4 range + count (no bare/unlabelled band).
    _PP4_RANGES = [">0.9", "0.8-0.9", "0.5-0.8"]
    b9, b89, b58 = _M["susie_bands"]
    bands = [(b9, PP4_TINTS[0], str(b9), ">0.9"),
             (b89, PP4_TINTS[1], str(b89), "0.8-0.9"),
             (b58, PP4_TINTS[2], str(b58), "0.5-0.8")]
    susie_targets = {}
    yy = susie_y0
    for cnt, col, lab, band in bands:
        hh = cnt*gscale
        susie_targets[band] = (yy, yy+hh)
        ax.add_patch(mpatches.Rectangle((xg-BW/2, yy), BW, hh, facecolor=col,
                                        edgecolor="white", lw=0.25, zorder=5))
        ax.text(xg+BW/2+0.015, yy+hh/2, lab, ha="left", va="center", fontsize=6, color="black", zorder=8)
        yy += hh
    # Single-signal-only block: the same PP.H4 bands in a hatched grey treatment
    # distinguish signal-model availability without assigning a causal hierarchy.
    ab9, ab89, ab58 = _M["abf_bands"]
    abf_bands = [(ab9, PP4_TINTS[0], str(ab9), ">0.9"),
                 (ab89, PP4_TINTS[1], str(ab89), "0.8-0.9"),
                 (ab58, PP4_TINTS[2], str(ab58), "0.5-0.8")]
    abf_targets = {}
    yy = abf_y0
    for cnt, col, lab, band in abf_bands:
        hh = cnt*gscale
        abf_targets[band] = (yy, yy+hh)
        ax.add_patch(mpatches.Rectangle((xg-BW/2, yy), BW, hh, facecolor=col,
                                        edgecolor="#6B7280", lw=0.3, hatch="////", alpha=0.35, zorder=5))
        ax.text(xg+BW/2+0.015, yy+hh/2, lab, ha="left", va="center", fontsize=6, color="black", zorder=8)
        yy += hh

    # Brackets explicitly bind each aggregate label to its three PP.H4 bands.
    def endpoint_bracket(y0, y1, label):
        bx, tick = xg + BW/2 + 0.19, 0.045
        ax.plot([bx, bx], [y0, y1], color="black", lw=0.6, clip_on=False, zorder=9)
        ax.plot([bx-tick, bx], [y0, y0], color="black", lw=0.6, clip_on=False, zorder=9)
        ax.plot([bx-tick, bx], [y1, y1], color="black", lw=0.6, clip_on=False, zorder=9)
        ax.text(bx+0.005, (y0+y1)/2, label, ha="left", va="center",
                fontsize=6, fontstretch="condensed", color="black", zorder=9)

    endpoint_bracket(susie_y0, susie_y1, f"Multi-signal\nCOLOC {SUSIE_GENES}")
    endpoint_bracket(abf_y0, abf_y1,
                     f"Single-signal\nonly {ABF_ONLY_GENES}")
    # Small union annotation above the column top, subordinate to the multi-signal total.
    ax.text(xg-BW/2, susie_y0-0.30, f"union {UNION_GENES:,}", ha="left", va="bottom",
            fontsize=6, color="black", zorder=8)
    # Shared PP.H4 legend in one compact horizontal row at the bottom-left.
    _pp4_leg = [mpatches.Patch(facecolor=PP4_TINTS[i], edgecolor="none", label=lab)
                for i, lab in enumerate([">0.9", "0.8–0.9", "0.5–0.8"])]
    ax.text(x_left + 0.02, H + 1.10, "PP.H4", ha="left", va="center",
            fontsize=6, color="black", zorder=10)
    _lg = ax.legend(handles=_pp4_leg, fontsize=6, loc="center left", frameon=False,
                    ncol=3, handlelength=0.60, handleheight=0.60,
                    handletextpad=0.25, columnspacing=0.45,
                    borderpad=0.0, borderaxespad=0.0,
                    bbox_to_anchor=(x_left + 0.31, H + 1.10),
                    bbox_transform=ax.transData)
    ax.add_artist(_lg)

    # ribbons: ancestry -> individual GWAS trait
    seg_tr = L["trait"][0]
    off = {a: anc_y[a][0] for a in ANCESTRIES}
    for a in ANCESTRIES:
        for t in r_order:
            if (a, t) in seg_tr:
                n = DAT[(a, t)][0]*s0
                hy0 = off[a]
                off[a] += n
                ry0, ry1 = seg_tr[(a, t)]
                _ribbon(ax, BW/2, hy0, hy0+n, XCOL["trait"]-BW/2,
                        ry0, ry1, ANCESTRY_COLOR[a], alpha=0.30)
    # ribbons: Trait -> Loci (ancestry-coloured)
    segf, segg = L["trait"][0], L["loci"][0]
    for a in ANCESTRIES:
        for t in r_order:
            if (a, t) in segf and (a, t) in segg:
                f0, f1 = segf[(a, t)]; g0, g1 = segg[(a, t)]
                _ribbon(ax, XCOL["trait"]+BW/2, f0, f1, XCOL["loci"]-BW/2, g0, g1, ANCESTRY_COLOR[a], alpha=0.30)
    # ribbons: Loci -> Genes, split & recoloured by method (SuSiE blue vs ABF-only grey/faint).
    # Each (ancestry,trait) loci block feeds a SuSiE sub-ribbon (into the shaded block) and an
    # ABF-only sub-ribbon (into the faded block); the source loci block is split proportionally.
    # (ancestry,trait) whose coloc genes come from an MVP-EUR stratum have no unique-loci block ->
    # the ribbon sources from the trait block across the empty loci column, drawn fainter.
    seg_lo = L["loci"][0]
    band_order = (">0.9", "0.8-0.9", "0.5-0.8")
    target_offset = {
        **{("susie", band): susie_targets[band][0] for band in band_order},
        **{("abf", band): abf_targets[band][0] for band in band_order},
    }
    for t in r_order:
        for a in ANCESTRIES:
            key = (a, t)
            ns, na = at_gs.get(key, 0), at_ga.get(key, 0)
            if ns == 0 and na == 0:
                continue
            ntot = ns + na
            if key in seg_lo:
                src0, src1, sx, amult = seg_lo[key][0], seg_lo[key][1], XCOL["loci"]+BW/2, 1.0
            else:
                src0, src1, sx, amult = seg_tr[key][0], seg_tr[key][1], XCOL["trait"]+BW/2, 0.5
            span = src1 - src0
            source_y = src0
            parts = ([('susie', band, at_gsb.get((a, t, band), 0)) for band in band_order] +
                     [('abf', band, at_gab.get((a, t, band), 0)) for band in band_order])
            for model, band, count in parts:
                if not count:
                    continue
                next_source_y = source_y + span*(count/ntot)
                target_y = target_offset[(model, band)]
                color = SUSIE_RIBBON if model == "susie" else ABF_COLOR
                alpha = (0.32 if model == "susie" else 0.16) * amult
                _ribbon(ax, sx, source_y, next_source_y, xg-BW/2,
                        target_y, target_y+count*gscale, color, alpha=alpha)
                target_offset[(model, band)] += count*gscale
                source_y = next_source_y

    # Never use bbox_inches="tight": it changes the PDF page box and was the
    # source of the repeated off-contract Figure 2A drift.
    save_args = {"dpi": 300, "facecolor": "white"}
    fig.savefig(path, **save_args)
    print(f"Saved: {path}")
    _write_cascade_source(os.path.splitext(path)[0] + "_source.tsv")
    if extra_paths:
        for ep in extra_paths:
            fig.savefig(ep, **save_args); print(f"Saved: {ep}")
            _write_cascade_source(os.path.splitext(ep)[0] + "_source.tsv")
    print(f"CAPTION (Fig2A cascade): GWAS evidence cascade — {N_GWAS} prespecified Tier-1/2 "
          f"GWAS strata across {len(ANCESTRIES)} ancestries and six liver traits -> "
          f"{tot['loci']} per-trait fine-mapped "
          f"loci ({N_UNIQUE_LOCI} globally unique; merged credible-set windows within each trait, +/-250 kb; "
          "a physical locus tagged by several enzymes recurs across those traits) -> colocalising genes. The "
          f"genes endpoint is split by signal model: the multi-signal COLOC block = {SUSIE_GENES} "
          "genes analyzed with SuSiE, shaded by PP.H4 (dark->light for >0.9 / 0.8-0.9 / 0.5-0.8 = "
          f"{b9} / {b89} / {b58}); the faded/hatched block = {ABF_ONLY_GENES} single-signal COLOC-only "
          "genes (ABF PP.H4>0.5 without multi-signal support), shaded the same way in grey (>0.9 / 0.8-0.9 / "
          f"0.5-0.8 = {ab9} / {ab89} / {ab58}); union = {UNION_GENES:,}. "
          "Stages 1-3 are coloured by GWAS ancestry; loci->gene ribbons are recoloured into the SuSiE "
          "(blue) versus single-signal-only (grey, faint) groups and terminate in their matching PP.H4 band. "
          "The molecular-QTL panel is European, so European GWAS are ancestry-matched and non-European GWAS "
          "are cross-ancestry eQTL-limited. The MVP-EUR strata (separate PolyFun-LD pipeline) carry coloc "
          "genes but no credible-set loci here, shown as faint trait->gene ribbons."
          + _untestable_caption())
    plt.close(fig)


def _untestable_caption():
    """Caption sentence for single-signal-only genes whose SuSiE test was untestable.

    n is read from the release's r2_summary.json and must equal the count in the
    rendered source table; no sentence when the release has no such gene."""
    n_drawn = sum(1 for model, state in _M["union_gene_state"].values()
                  if model == "single_signal_only" and state == "untestable")
    if not _COLOC_SUMMARY:
        if n_drawn:
            raise ValueError(f"{n_drawn} single-signal-only genes are SuSiE-untestable; "
                             "set FIG2_COLOC_SUMMARY to the release's r2_summary.json")
        return ""
    with open(_COLOC_SUMMARY) as handle:
        n = _json.load(handle)["fig2_susie_untestable"]["drawn_as_single_signal_only"]
    if n != n_drawn:
        raise ValueError(f"r2_summary.json gives {n} untestable single-signal-only genes; "
                         f"the cascade counts {n_drawn}")
    if not n:
        return ""
    return (f" {n} genes shown as single-signal only could not be tested by multi-signal "
            "colocalization (shared posterior below 0.5); they are not multi-signal negatives.")


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

    def _write_source(pdf_path):
        source_path = os.path.splitext(pdf_path)[0] + "_source.tsv"
        with open(source_path, "w", newline="") as handle:
            writer = _csv.writer(handle, delimiter="\t")
            writer.writerow([
                "study_label", "ancestry", "fine_mapped_loci", "coloc_genes",
                "coloc_definition", "main_gwas_strata", "displayed_studies",
                "coloc_input",
            ])
            for label, ancestry, loci, genes in ROWS:
                writer.writerow([
                    label, ancestry, loci, genes,
                    "named gene; multi-signal or single-signal-only PP.H4 > 0.5",
                    N_GWAS, len(ROWS), _GENE_COLOC,
                ])

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
    _write_source(path)
    print(f"Saved: {path}")
    if extra_paths:
        for ep in extra_paths:
            fig.savefig(ep, bbox_inches="tight", dpi=300, facecolor="white"); print(f"Saved: {ep}")
            _write_source(ep)
    print(f"CAPTION (Fig2A per-study): Fine-mapped loci (log scale) and colocalising genes "
          f"(multi-signal or single-signal-only COLOC, PP.H4>0.5) per Tier-1/2 "
          f"(liver-specific) GWAS study, for the {n} "
          f"MAIN studies carrying a fine-mapping summary, grouped by {len(ANCESTRIES)} ancestries "
          "(incl. MVP NAFLD/ALT/AST AFR/AMR/EAS). The loci-vs-genes contrast exposes the "
          "European-eQTL ancestry bottleneck at the study level: non-European studies (MVP/BBJ/"
          "PanUKBB AFR/AMR/EAS) yield many fine-mapped loci but few coloc genes, whereas European "
          "UKBB studies yield few loci but many genes. The MVP-EUR strata (fine-mapped under a "
          f"separate PolyFun-LD pipeline) are not shown here but contribute to the {UNION_GENES} "
          "MAIN coloc-gene total.")
    plt.close(fig)


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    candidate_dir = os.environ.get("FIG2_CANDIDATE_DIR")
    panel_mode = os.environ.get("FIG2_PANEL", "cascade")
    if candidate_dir:
        os.makedirs(candidate_dir, exist_ok=True)
        if panel_mode not in ("cascade", "perstudy", "all"):
            raise SystemExit("FIG2_PANEL must be 'cascade', 'perstudy', or 'all'")
        if panel_mode in ("cascade", "all"):
            build_cascade(os.path.join(candidate_dir, "fig2A_gwas_cascade.pdf"))
        if panel_mode in ("perstudy", "all"):
            build_study_bars(os.path.join(candidate_dir, "FigS2A_gwas_perstudy.pdf"))
        raise SystemExit(0)

    build_A(_out_path("gwas_creative_A_binary_grid.pdf"))
    PANELS = os.path.join(_project_root(), "figures", "main", "fig2_genetics", "panels")
    build_cascade(_out_path("gwas_creative_E_cascade.pdf"),
                  extra_paths=[os.path.join(PANELS, "fig2A_gwas_cascade.pdf")])
    build_study_bars(_out_path("gwas_creative_F_perstudy.pdf"),
                     extra_paths=[os.path.join(PANELS, "FigS2A_gwas_perstudy.pdf")])
    print("All done.")
