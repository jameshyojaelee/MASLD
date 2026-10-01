#!/usr/bin/env python3
"""B-COLOC: external signed direction for the COLOC genes of two releases.

Inputs (defaults below; all hg19 unless stated)
  adopted release 2026-08-17   GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv
    its signal pairs (replay)  GWAS/finemapping/results/coloc_eligibility_offline_20260923T211525Z/eqtl_overlap_check.tsv
  candidate corrected release  GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv
    its signal pairs           .../susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz
  study tiers, ancestry        GWAS/finemapping/config/gwas_trait_tier.tsv, gwas_registry.tsv, gwas_af_sources.tsv
  GWAS summary statistics      the file COLOC read (sumstats_af/ copy when listed)
  Broadaway liver eQTL meta    data/broadaway_eqtl/chr*_marginal_summary_results.tsv; Beta is the EA effect
  eQTL SuSiE fits COLOC read   GWAS/finemapping/results/eqtl_susie_polyfun (read by eqtl_conditional_sign.R)
  1kGP high coverage, phased   data/external/allelic_refs/kgp_highcov (GRCh38); founders (no parent in the ped file)
  hg19 -> GRCh38 chain         data/broadaway_eqtl/hg19ToHg38.over.chain (UCSC liftOver)
  GTEx v8 liver signif pairs   check column only
  Model A tag lists            TAG_LISTS below: T3 tags_final, T3' kept, T1 v2, T3' v2 kept (positions only)

Outputs (--out-dir, must not exist)
  direction_by_trait.tsv       one row per gene x release x trait (the deliverable)
  direction_by_study.tsv       one row per gene x release x study
  direction_signal_pairs.tsv   every colocalized signal pair (PP.H4 > 0.5) with its call
  gwas_orientation_audit.tsv   per study: which allele the file's beta belongs to, and the evidence
  tag_ld_to_coloc_lead.tsv     EUR r between each Model A tag and each colocalized Broadaway lead
  tag_ld_gene_summary.tsv      per gene x release and tag list: any tag with r2 >= 0.8 to a colocalized
                               lead (trait-level flags are columns of direction_by_trait.tsv)
  positive_controls_gwas.tsv   sign of known trait-raising alleles in each scanned GWAS file
  summary.json, params.json, input_manifest.tsv

Rules
  * Genes: tier-1/2 studies, SuSiE PP.H4 > 0.5, Ensembl ids (adopted 462, candidate 436). Every
    signal pair with PP.H4 > 0.5 is evaluated (adopted: every pair, as that run had no
    eligibility check; candidate: eligible pairs only).
  * hit1 = GWAS credible-set lead, hit2 = eQTL credible-set lead. At hit1 the GWAS row is the
    lowest-p row at that position whose alleles match a Broadaway row for the gene (COLOC's merge).
  * GWAS effect allele. The files' `beta` is assumed to belong to `allele1` only after a
    per-study audit: sign agreement with a reference study of the same trait at independent
    (one per Mb) non-palindromic SNVs with p < 5e-8 in either file and p < 0.01 in the other,
    alleles matched by letters.
    >= 0.8 of >= 10 loci -> allele1; <= 0.2 -> allele2 (the beta belongs to allele2);
    otherwise, or when a known raising allele at p < 5e-8 (PNPLA3 rs738409 G, TM6SF2
    rs58542926 T; ALT/AST/NAFLD/NASH/PDFF only) contradicts the call, the study is unresolved
    and none of its rows is directional. The beta is never negated; the audit picks the letters.
  * Orientation is applied once: sign = sign(GWAS beta of its effect allele at hit1)
    x sign(haplotype r between that allele and the eQTL effect allele at hit2)
    x sign(eQTL Beta of its effect allele at hit2). Alleles are matched by letters in both
    orders; no REF/ALT column of any file is assumed to be reference-oriented.
  * r comes from phased 1kGP haplotypes of the study's superpopulation (PanUKBB CSA -> SAS),
    founders (no recorded parent) minus the KING-related samples T1 v2 drops (the same EUR set
    that selects the T1 v2 tags). Same lead position: r = +1 or -1 from the letters. Different
    leads: r2 >= 0.5, and for a non-EUR study also r2 >= 0.5 with the same sign in EUR (the
    eQTL cohort is European).
  * Technical blocks (the lead cannot be oriented):
    - either lead (GRCh38) inside excluded_spans_mhc_ig_tr.bed (MHC, IG, TR), as in the tag layer;
    - a palindromic SNV lead (A/T, C/G) needs a 1kGP record, MAF <= 0.4 in every available
      frequency, and strand evidence: the file frequency fits the 1kGP strand (distance <= 0.15,
      and not closer to the opposite strand by >= 0.10); or, when the file has no frequency, the
      file's allele1 follows the study's measured allele-order convention (UKBB ALT first; MVP,
      BBJ REF first; p0-orientation-check study_allele_conventions.tsv). Otherwise blocked;
    - different leads only (the 1kGP record is used for r): a non-palindromic lead whose source
      frequency differs from the 1kGP record by > 0.2 is blocked (the record may not be the
      file's variant). Source frequency: the GWAS file at hit1 (Broadaway EAF for an EUR study
      without a file frequency), Broadaway EAF at hit2.
  * Evidence blocks (the marginal sign is not evidence of the colocalized signal's sign;
    prespecified 2026-09-29 before the rerun):
    - GWAS marginal p >= 1e-3 at hit1, or Broadaway marginal p >= 1e-3 at hit2;
    - the colocalized eQTL SuSiE component's posterior effect at hit2 (mu[l, hit2]) and the
      fit's marginal there (XtXr) have opposite signs. The GWAS SuSiE fits were not saved, so
      the GWAS side has the p floor only.
    Each pair also carries sign_before_evidence_rules (the sign with only technical blocks).
  * A gene x study takes the agreed sign of its directional signal pairs, and a gene x trait
    the agreed sign of its directional studies; disagreement -> not_directional. Blocked pairs
    do not vote, with one exception (review round 1, 2026-09-30): a pair blocked only by the
    eQTL component rule (both marginals pass the floor, no technical block) whose sign before
    that rule opposes the agreed sign makes the gene x study, and the gene x trait, not
    directional. The reported leads come from the highest-PP.H4 directional pair (study), else
    the highest-PP.H4 one.
  * Labels: enzyme traits read 'liver-enzyme-raising haplotype'.
  * Provenance of every row: 'external eQTL, not MASLD tissue'.
"""
import argparse
import hashlib
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path("/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design")
FM = REPO / "GWAS/finemapping"
KGP = REPO / "data/external/allelic_refs/kgp_highcov"
HERE = Path(__file__).resolve().parent
RELEASES = {
    "adopted_2026-08-17": dict(
        master=FM / "results/susie_coloc/susie_coloc_all_gwas.csv",
        pairs=FM / "results/coloc_eligibility_offline_20260923T211525Z/eqtl_overlap_check.tsv",
        expected_genes=462, eligible_only=False),
    "candidate_r2_20260924T053244Z": dict(
        master=FM / "results/susie_coloc_r2_20260924T053244Z/susie_coloc_all_gwas.csv",
        pairs=FM / "results/susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz",
        expected_genes=436, eligible_only=True),
}
PAIR_SOURCE = {
    "adopted_2026-08-17": "hit1/hit2 from the 2026-09-23 replay (coloc_eligibility_offline_20260923T211525Z/"
                          "eqtl_overlap_check.tsv): the adopted run saved no pair records; pairs ineligible "
                          "under the corrected rule are included by design",
    "candidate_r2_20260924T053244Z": "hit1/hit2 from susie_coloc_r2_20260924T053244Z/susie_coloc_signal_pairs.tsv.gz; "
                                     "eligible pairs only",
}
PP4_MIN, MIN_R2, MAX_PAL_MAF, PAL_MARGIN, PAL_MAX_DIST, TAG_R2, MIN_HAPS = 0.5, 0.5, 0.4, 0.10, 0.15, 0.8, 50
# evidence floor at both leads (prespecified 2026-09-29 before rerun v4); 1e-5 is reported as sensitivity only
MARGINAL_P_MAX, MARGINAL_P_SENSITIVITY = 1e-3, 1e-5
FREQ_MAX_DIFF = 0.2          # source vs 1kGP frequency of a non-palindromic lead used for LD
TAG_R2_BAND = (0.75, 0.85)   # reported: genes whose best tag r2 sits near TAG_R2
EXEC = REPO / "Analysis/MASLD_Model_Benchmark/executions"
EXCLUDED_SPANS = EXEC / "in-biopsy-allelic-tags-20260925T194248Z/excluded_spans_mhc_ig_tr.bed"
ALLELE_CONVENTIONS = FM / "results/alphagenome_atlas/p0-orientation-check-20260914T230807Z/tables/study_allele_conventions.tsv"
KING_EXCLUDED = EXEC / "in-biopsy-allelic-tags-v2-20260929T134419Z/founders/exclude_related.txt"
TAG_LISTS = {   # name: (tag file, description); t3prime_kept is derived from t3_tags_final
    "t3_tags_final": (EXEC / "in-biopsy-allelic-tags-20260925T194248Z/t3/tags_final.tsv",
                      "T1 v1 + T2 + T3 tags_final (2026-09-25)"),
    "t3prime_kept": (EXEC / "in-biopsy-allelic-obs-20260926T155651Z/t3prime/tag_exclusions.tsv",
                     "t3_tags_final rows with kept == True after T3' (706-contig mappability, 2026-09-29)"),
    "t1v2": (EXEC / "in-biopsy-allelic-tags-v2-20260929T134419Z/tags",
             "T1 v2 tags (corrected founders + KING exclusions, 2026-09-29), before T2/T3 filters"),
    "t3prime_v2_kept": (EXEC / "in-biopsy-allelic-obs-20260926T155651Z/t3prime/v2-20260929T200244Z/tag_exclusions_v2.tsv",
                        "t1v2 rows with kept_steps_1_4 == True after T3' v2 (spans, gene/site filters, simulated "
                        "read mappability; 2026-09-29)"),
}
SUPERPOPS = ("EUR", "EAS", "AFR", "AMR", "SAS")
PROVENANCE = "external eQTL, not MASLD tissue"
COMPONENT_REASON = "eQTL marginal sign at hit2 opposes the colocalized SuSiE component's effect"
OPPOSING_COMPONENT_REASON = "a pair blocked only by the eQTL component rule has the opposite sign"
TRAIT_PHRASE = {"ALT": "liver-enzyme", "AST": "liver-enzyme", "GGT": "liver-enzyme",
                "NAFLD": "NAFLD-risk", "NASH": "NASH-risk", "PDFF": "liver-fat (PDFF)"}
# (label, hg19 chromosome, hg19 position, trait-raising allele, other allele)
POSITIVE_CONTROLS = [("TM6SF2 rs58542926 E167K", "19", 19379549, "T", "C"),
                     ("PNPLA3 rs738409 I148M", "22", 44324727, "G", "C")]
POSITIVE_CONTROL_TRAITS = {"ALT", "AST", "NAFLD", "NASH", "PDFF"}
# GWAS orientation audit: each study against a reference of its trait; each reference against a
# second study, so every reference is itself audited.
AUDIT_REFERENCE_BY_TRAIT = {"ALT": "UKBB_ALT", "AST": "UKBB_AST", "GGT": "UKBB_GGT", "NAFLD": "MVP_NAFLD_EUR",
                            "NASH": "MVP_NAFLD_EUR", "PDFF": "2021_34128465_PDFF_EUR"}
AUDIT_REFERENCE_OF_REFERENCE = {"UKBB_ALT": "MVP_ALT_EUR", "UKBB_AST": "MVP_AST_EUR", "UKBB_GGT": "UKBB_ALT",
                                "MVP_NAFLD_EUR": "FinnGen_NAFLD", "2021_34128465_PDFF_EUR": "MVP_NAFLD_EUR"}
AUDIT_REF_P, AUDIT_STUDY_P, AUDIT_MIN_LOCI, AUDIT_AGREE, AUDIT_WINDOW = 5e-8, 1e-2, 10, 0.8, 1_000_000
GWAS_COLS = ["chromosome", "position", "allele1", "allele2", "beta", "se", "pval", "af"]
COMPLEMENT = str.maketrans("ACGT", "TGCA")
EQTL_COLS = ["ENSG", "CHR", "POS", "NEA", "EA", "EAF", "Beta", "SE", "PVAL"]


# ----------------------------------------------------------------------------- rules

def is_palindromic(a, b):
    return {a, b} in ({"A", "T"}, {"C", "G"})


def finite(x):
    return x is not None and np.isfinite(x)


def palindrome_block(a, b, file_freq, panel_freq, lead="lead", convention_ok=None):
    """Reason a palindromic SNV lead cannot be oriented, or None.

    file_freq and panel_freq are the frequency of allele `a` in the source file and in the 1kGP
    superpopulation (either may be NaN). Strand evidence is needed besides the 1kGP MAF: the
    file frequency against 1kGP, or, without a file frequency, convention_ok (True when the
    file's allele1 follows the study's measured allele-order convention at this 1kGP record,
    False when it does not, None when the study has no convention)."""
    if not is_palindromic(a, b):
        return None
    if not finite(panel_freq):
        return f"palindromic {lead} without a 1kGP frequency (strand unknown)"
    freqs = [f for f in (file_freq, panel_freq) if finite(f)]
    if max(min(f, 1 - f) for f in freqs) > MAX_PAL_MAF:
        return f"palindromic {lead} with MAF > {MAX_PAL_MAF}"
    if finite(file_freq):
        d_same, d_opp = abs(file_freq - panel_freq), abs(file_freq - (1 - panel_freq))
        if min(d_same, d_opp) > PAL_MAX_DIST:
            return f"palindromic {lead}: file and 1kGP frequencies disagree"
        if d_same - d_opp >= PAL_MARGIN:
            return f"palindromic {lead}: file frequency fits the opposite strand"
        return None
    if convention_ok is True:
        return None
    if convention_ok is False:
        return f"palindromic {lead}: allele1 breaks the file's allele-order convention (strand unknown)"
    return f"palindromic {lead} with no file frequency and no allele-order convention (strand unknown)"


def frequency_block(a, b, source_freq, panel_freq, lead="lead"):
    """Reason the 1kGP record matched by letters may not be the source file's variant, or None.

    Applied to non-palindromic leads whose 1kGP record is used for r (different leads)."""
    if is_palindromic(a, b) or not (finite(source_freq) and finite(panel_freq)):
        return None
    if abs(source_freq - panel_freq) > FREQ_MAX_DIFF:
        return f"{lead}: source and 1kGP frequencies differ by > {FREQ_MAX_DIFF} (1kGP record may not be the file variant)"
    return None


def marginal_block(p, lead, source):
    """Evidence floor: the marginal at the lead must have p < MARGINAL_P_MAX."""
    if not finite(p) or p >= MARGINAL_P_MAX:
        return f"{source} marginal p >= {MARGINAL_P_MAX:g} at {lead}"
    return None


def marginal_floor_label(gwas_p, eqtl_p):
    fail = [n for n, p in (("GWAS", gwas_p), ("eQTL", eqtl_p)) if not finite(p) or p >= MARGINAL_P_MAX]
    return "pass" if not fail else "fail: " + " and ".join(fail)


def convention_followed(convention, effect_column, panel_match):
    """True/False: does the file's allele1 sit where the study's convention puts it; None if no convention.

    panel_match is 'a=ALT'/'a=REF' for the GWAS effect allele at the 1kGP record."""
    if convention not in ("allele1_is_alt", "allele1_is_ref") or panel_match not in ("a=ALT", "a=REF") \
            or effect_column not in ("allele1", "allele2"):
        return None
    allele1_is_alt = (panel_match == "a=ALT") == (effect_column == "allele1")
    return allele1_is_alt if convention == "allele1_is_alt" else not allele1_is_alt


def excluded_span(spans, chrom, pos):
    """Name of the excluded span (BED, GRCh38) holding chrom:pos, or ''."""
    for c, start, end, name in spans:
        if c == chrom and start < pos <= end:
            return name
    return ""


def linkage(same_lead, a1, a2, ea, nea, r_study, r_eur, study_pop, ld_note=""):
    """(r between GWAS effect allele a1 at hit1 and the eQTL effect allele at hit2, blocking reason or None)."""
    if same_lead:
        if (a1, a2) == (ea, nea):
            return 1.0, None
        if (a1, a2) == (nea, ea):
            return -1.0, None
        return np.nan, "alleles differ between GWAS and eQTL at the shared lead"
    if not np.isfinite(r_study):
        return np.nan, f"no LD estimate between leads in {study_pop} ({ld_note or 'unknown'})"
    if r_study ** 2 < MIN_R2:
        return r_study, f"weak LD between leads in {study_pop} (r2 < {MIN_R2})"
    if study_pop != "EUR":
        if not np.isfinite(r_eur) or r_eur ** 2 < MIN_R2:
            return r_study, f"weak or missing LD between leads in EUR, the eQTL population (r2 < {MIN_R2})"
        if np.sign(r_eur) != np.sign(r_study):
            return r_study, "LD sign between leads differs between the study superpopulation and EUR"
    return r_study, None


def direction_sign(beta_effect, r_effect_ea, beta_ea, blocks):
    """(+1 / -1 / 0, reason). The only place the three signs are combined.

    +1 means the haplotype carrying the trait-raising allele at the GWAS lead carries the
    expression-raising allele at the eQTL lead."""
    for b in blocks:
        if b:
            return 0, b
    for name, v in (("GWAS beta", beta_effect), ("LD r", r_effect_ea), ("eQTL beta", beta_ea)):
        if v is None or not np.isfinite(v) or v == 0:
            return 0, f"{name} missing or zero"
    return int(np.sign(beta_effect) * np.sign(r_effect_ea) * np.sign(beta_ea)), ""


def direction_label(trait, sign):
    if sign == 0:
        return "not_directional"
    return f"{TRAIT_PHRASE.get(trait, trait)}-raising haplotype {'increases' if sign > 0 else 'decreases'} expression"


def gwas_effect_alleles(allele1, allele2, af_allele1, column):
    """(effect allele, other allele, effect-allele frequency) for the audited effect-allele column.

    The af column of the sumstats_af copies is the frequency of allele1."""
    if column == "allele2":
        return allele2, allele1, 1 - af_allele1
    return allele1, allele2, af_allele1


def conditional_sign_block(conditional_equals_marginal):
    """Blocking reason when the eQTL SuSiE component's effect at hit2 opposes the marginal there."""
    if conditional_equals_marginal is False:
        return COMPONENT_REASON
    return None


def orientation_concordance(ref, study):
    """(n loci, share of loci where `study` agrees in sign with `ref`).

    ref and study hold non-palindromic SNV rows (chromosome, position, allele1, allele2, beta_num,
    pval_num). Loci: positions with the same two letters in both files, p < AUDIT_REF_P in either
    file and p < AUDIT_STUDY_P in the other, thinned to the lowest-p locus per chromosome x
    AUDIT_WINDOW. Counting hits of either file keeps a small reference (few genome-wide hits)
    from leaving a large study unresolved; sign agreement is symmetric in the two files."""
    m = ref.merge(study, on=["chromosome", "position"], suffixes=("_r", "_s"))
    hit = ((m.pval_num_r < AUDIT_REF_P) & (m.pval_num_s < AUDIT_STUDY_P)) | \
          ((m.pval_num_s < AUDIT_REF_P) & (m.pval_num_r < AUDIT_STUDY_P))
    same = (m.allele1_r == m.allele1_s) & (m.allele2_r == m.allele2_s)
    swap = (m.allele1_r == m.allele2_s) & (m.allele2_r == m.allele1_s)
    m = m.assign(beta_s_aligned=np.where(same, m.beta_num_s, -m.beta_num_s))[hit & (same | swap)]
    m = m.assign(window=m.position // AUDIT_WINDOW, p_min=np.minimum(m.pval_num_r, m.pval_num_s)) \
         .sort_values(["p_min", "chromosome", "position"], kind="stable")
    m = m.drop_duplicates(["chromosome", "window"])
    if m.empty:
        return 0, np.nan
    return int(len(m)), float((np.sign(m.beta_num_r) == np.sign(m.beta_s_aligned)).mean())


def effect_allele_column(n_loci, agree, control_signs):
    """'allele1', 'allele2' or 'unresolved' for one study.

    control_signs: sign of the allele1-convention beta of a known raising allele at p < 5e-8.
    With control_signs = [] this is the call from sign concordance alone."""
    if n_loci >= AUDIT_MIN_LOCI and agree >= AUDIT_AGREE:
        column, expected = "allele1", 1
    elif n_loci >= AUDIT_MIN_LOCI and agree <= 1 - AUDIT_AGREE:
        column, expected = "allele2", -1
    else:
        return "unresolved"
    if any(s != expected for s in control_signs):
        return "unresolved"
    return column


def combine_rows(d, keys, rank, tiebreak, label):
    """One row per `keys` with the agreed sign of its directional rows.

    The reported row is the highest-`rank` directional row, else the highest-`rank` row. Rows
    whose directional members disagree become not_directional with the reason 'label disagree'."""
    d = d.assign(_directional=d.direction_sign != 0) \
         .sort_values(["_directional", rank] + tiebreak, ascending=[False, False] + [True] * len(tiebreak),
                      kind="stable")
    agg = d.groupby(keys).direction_sign.agg(
        n_total="size", n_up=lambda s: int((s > 0).sum()), n_down=lambda s: int((s < 0).sum())).reset_index()
    top = d.drop_duplicates(keys).drop(columns="_directional").merge(agg, on=keys)
    split = (top.n_up > 0) & (top.n_down > 0)
    top.loc[split, "direction_sign"] = 0
    top.loc[split, "not_directional_reason"] = f"{label} disagree"
    if "eqtl_beta_of_hit2_allele_on_trait_raising_haplotype" in top:
        top.loc[split, "eqtl_beta_of_hit2_allele_on_trait_raising_haplotype"] = np.nan
    return top


def same_lead_check(d):
    """At a shared lead the 1kGP r must be +1/-1 with the sign the letters give (panel matching check)."""
    s = d[d.same_lead & d.ld_r_effect_ea_study.notna() & d.ld_r_used.notna()]
    return {"n": int(len(s)), "agree": int((np.sign(s.ld_r_effect_ea_study) == np.sign(s.ld_r_used)).sum()),
            "min_abs_r": float(s.ld_r_effect_ea_study.abs().min()) if len(s) else None}


def release_agreement(trait_tab):
    """Gene x trait rows present in both releases: how many are directional in both, and agree."""
    rel = sorted(trait_tab.release.unique())
    if len(rel) != 2:
        return {}
    x = trait_tab[trait_tab.release == rel[0]].merge(trait_tab[trait_tab.release == rel[1]],
                                                     on=["ensembl", "trait"], suffixes=("_a", "_b"))
    both = x[(x.direction_sign_a != 0) & (x.direction_sign_b != 0)]
    same_inputs = (both.reported_study_a == both.reported_study_b) & (both.hit1_a == both.hit1_b) & \
        (both.hit2_a == both.hit2_b)
    diff = both[~same_inputs]
    return {"shared_rows": int(len(x)), "directional_in_both": int(len(both)),
            "same_sign": int((both.direction_sign_a == both.direction_sign_b).sum()),
            "directional_in_both_same_study_and_leads (agreement by construction)": int(same_inputs.sum()),
            "directional_in_both_different_study_or_leads": int(len(diff)),
            "same_sign_where_study_or_leads_differ (informative)": int((diff.direction_sign_a == diff.direction_sign_b).sum())}


# ----------------------------------------------------------------------------- IO helpers

def read_str(path, **kw):
    """Read everything as text; 'NA', 'None' and '' stay literal until converted on purpose."""
    return pd.read_csv(path, sep=kw.pop("sep", "\t"), dtype=str, keep_default_na=False, **kw)


def to_num(s):
    return pd.to_numeric(s, errors="coerce")


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 22), b""):
            h.update(block)
    return h.hexdigest()


def awk_filter(src, keys, dst, key_expr):
    """Header plus rows of `src` whose key_expr (an awk expression) is in `keys`."""
    keyfile = Path(str(dst) + ".keys")
    keyfile.write_text("\n".join(sorted(keys)) + "\n")
    prog = f'BEGIN{{FS=OFS="\\t"}} NR==FNR{{k[$1]=1; next}} FNR==1{{print; next}} (({key_expr}) in k)'
    with open(dst, "w") as fh:
        subprocess.run(["awk", prog, str(keyfile), str(src)], stdout=fh, check=True)
    keyfile.unlink()
    return dst


def awk_p_below(src, dst, pmax):
    """Header plus rows of a GWAS_COLS file whose pval (column 7) is numeric and below pmax."""
    with open(src) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if header != GWAS_COLS:
        raise SystemExit(f"unexpected GWAS header in {src}: {header}")
    prog = f'BEGIN{{FS=OFS="\\t"}} FNR==1{{print; next}} ($7 ~ /^[0-9.]+([eE][-+]?[0-9]+)?$/) && ($7 + 0 < {pmax})'
    with open(dst, "w") as fh:
        subprocess.run(["awk", prog, str(src)], stdout=fh, check=True)
    return dst


def gwas_path(study, registry, af_sources):
    """The file 06_susie_coloc.R read: the sumstats_af copy when listed and present."""
    hit = af_sources[af_sources.study_name == study]
    if len(hit) == 1 and Path(hit.af_sumstats_path.iloc[0]).exists():
        return Path(hit.af_sumstats_path.iloc[0])
    return FM / registry.loc[registry.study_name == study, "sumstats_path"].iloc[0]


def founders(ped, superpop):
    """Founders of one superpopulation (no father or mother in the ped file), trio parents kept,
    as in tags/t1_tags.py v2. KING-related samples are removed through Panel's masks."""
    p = pd.read_csv(ped, sep=r"\s+", dtype=str)
    f = p[(p["Superpopulation"] == superpop) & (p["FatherID"] == "0") & (p["MotherID"] == "0")]
    return list(f["SampleID"])


# ----------------------------------------------------------------------------- inputs

def load_release(name, spec, tier):
    """Positive gene x study rows and their signal pairs for one release."""
    m = read_str(spec["master"], sep=",", usecols=["gwas_name", "gene", "ensembl", "chr", "PP.H4.susie"])
    m["pp4"] = to_num(m["PP.H4.susie"])
    m = m[m.gwas_name.map(tier).isin([1, 2]) & m.ensembl.str.match(r"^ENSG\d{11}$") & (m.pp4 > PP4_MIN)]
    m = m.rename(columns={"gwas_name": "study"}).drop(columns="PP.H4.susie")
    if m.duplicated(["study", "ensembl"]).any():
        raise SystemExit(f"{name}: duplicate study x gene rows in {spec['master']}")
    n_genes = m.ensembl.nunique()
    if n_genes != spec["expected_genes"]:
        raise SystemExit(f"{name}: {n_genes} genes, expected {spec['expected_genes']}")
    p = read_str(spec["pairs"])
    p = p.rename(columns={"gwas_name": "study", "gwas": "study", "replay_PP.H4": "pair_pp4", "PP.H4": "pair_pp4"})
    p["pair_pp4"] = to_num(p["pair_pp4"])
    p["eligible_corrected_rule"] = p["eligible"]
    if spec["eligible_only"]:
        p = p[p.eligible == "TRUE"]
    p = p[["study", "ensembl", "idx1", "idx2", "hit1", "hit2", "pair_pp4", "eligible_corrected_rule"]]
    p = p.merge(m[["study", "ensembl"]], on=["study", "ensembl"])
    best = p.groupby(["study", "ensembl"]).pair_pp4.max().rename("max_pair_pp4").reset_index()
    m = m.merge(best, on=["study", "ensembl"], how="left")
    m["pp4_reproduced_by_pairs"] = (m.max_pair_pp4 - m.pp4).abs() < 1e-3
    p = p.merge(m[["study", "ensembl", "gene", "chr", "pp4", "max_pair_pp4"]], on=["study", "ensembl"])
    p = p[p.pair_pp4 > PP4_MIN].copy()
    p["is_highest_pp4_pair"] = False
    top = p.sort_values(["pair_pp4", "idx1", "idx2"], ascending=[False, True, True]) \
           .drop_duplicates(["study", "ensembl"]).index
    p.loc[top, "is_highest_pp4_pair"] = True
    p.insert(0, "release", name)
    m.insert(0, "release", name)
    return m, p


def liftover(positions, chain, liftover_bin, work):
    """{(chr, hg19 pos): (GRCh38 chrom, pos, strand)} for positions that stay on their chromosome."""
    bed_in, bed_out, unm = work / "hg19.bed", work / "hg38.bed", work / "unmapped.bed"
    with open(bed_in, "w") as fh:
        for c, pos in sorted(positions, key=lambda x: (int(x[0]), x[1])):
            fh.write(f"chr{c}\t{pos - 1}\t{pos}\t{c}:{pos}\t0\t+\n")
    subprocess.run([liftover_bin, str(bed_in), str(chain), str(bed_out), str(unm)], check=True,
                   capture_output=True)
    out = {}
    for line in open(bed_out):
        f = line.rstrip("\n").split("\t")
        c, pos = f[3].split(":")
        if f[0] == f"chr{c}":
            out[(c, int(pos))] = (f[0], int(f[2]), f[5])
    return out


def panel_haplotypes(chrom, positions, samples, work):
    """{pos: [(REF, ALT, haplotypes)]} with haplotypes interleaved per sample, ALT = 1, -1 missing."""
    vcf = KGP / f"1kGP_high_coverage_Illumina.{chrom}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    reg = work / f"regions_{chrom}.tsv"
    reg.write_text("".join(f"{chrom}\t{p}\n" for p in sorted(set(positions))))
    cmd = ["bcftools", "query", "-R", str(reg), "-s", ",".join(samples), "-i", "N_ALT=1",
           "-f", "%POS\t%REF\t%ALT[\t%GT]\n", str(vcf)]
    code = {"0": 0, "1": 1}
    out = {}
    wanted = set(positions)
    for line in subprocess.run(cmd, capture_output=True, text=True, check=True).stdout.splitlines():
        f = line.split("\t")
        pos = int(f[0])
        if pos not in wanted:
            continue
        h = np.array([[code.get(g[0], -1), code.get(g[2], -1)] if len(g) == 3 and g[1] == "|" else [-1, -1]
                      for g in f[3:]], dtype=np.int8).reshape(-1)
        out.setdefault(pos, []).append((f[1], f[2], h))
    return out


class Panel:
    """Allele-indicator haplotypes from the 1kGP panel, matched by letters in both orders."""

    def __init__(self, haps, samples, sample_pop, excluded=frozenset()):
        """Masks: each superpopulation without `excluded` samples; 'EUR_all_founders' keeps them."""
        self.haps = haps
        pop = np.repeat(np.array([sample_pop[s] for s in samples]), 2)
        kept = np.repeat(np.array([s not in excluded for s in samples]), 2)
        self.mask = {sp: (pop == sp) & kept for sp in SUPERPOPS}
        self.mask["EUR_all_founders"] = pop == "EUR"
        self.n_samples = {k: int(v.sum() // 2) for k, v in self.mask.items()}

    def indicator(self, chrom, pos, strand, a, b):
        """(indicator of allele `a` or None, 'a=ALT'/'a=REF' or reason)."""
        if strand == "-":
            a, b = a.translate(COMPLEMENT)[::-1], b.translate(COMPLEMENT)[::-1]
        for ref, alt, h in self.haps.get(chrom, {}).get(pos, []):
            if (a, b) == (alt, ref):
                return h.copy(), "a=ALT"
            if (a, b) == (ref, alt):
                x = h.copy()
                x[h >= 0] = 1 - h[h >= 0]
                return x, "a=REF"
        return None, "allele pair not in 1kGP panel"

    def freq(self, x, pop):
        v = x[self.mask[pop] & (x >= 0)]
        return float(v.mean()) if len(v) >= MIN_HAPS else np.nan

    def r(self, x, y, pop):
        ok = self.mask[pop] & (x >= 0) & (y >= 0)
        if ok.sum() < MIN_HAPS:
            return np.nan
        xv, yv = x[ok].astype(float), y[ok].astype(float)
        if xv.std() == 0 or yv.std() == 0:
            return np.nan
        return float(np.corrcoef(xv, yv)[0, 1])


# ----------------------------------------------------------------------------- GWAS orientation audit

def audit_studies(pair_studies, trait):
    """{study: reference study}, closed so that every reference is audited too."""
    ref_of = {}
    todo = list(pair_studies)
    while todo:
        s = todo.pop()
        if s in ref_of:
            continue
        ref_of[s] = AUDIT_REFERENCE_OF_REFERENCE.get(s) or AUDIT_REFERENCE_BY_TRAIT[trait[s]]
        todo.append(ref_of[s])
    return ref_of


def load_snv_extract(path):
    d = read_str(path)
    d = d[d.allele1.str.fullmatch("[ACGT]") & d.allele2.str.fullmatch("[ACGT]")]
    d = d[~np.array([is_palindromic(a, b) for a, b in zip(d.allele1, d.allele2)], dtype=bool)]
    return pd.DataFrame(dict(chromosome=d.chromosome, position=d.position.astype(int), allele1=d.allele1,
                             allele2=d.allele2, beta_num=to_num(d.beta), pval_num=to_num(d.pval)))


def run_orientation_audit(ref_of, trait, controls, registry, af_sources, work, threads, log):
    """Per-study effect-allele column with its evidence (DataFrame, one row per study)."""
    def scan(st):
        dst = work / f"audit_p_lt_{AUDIT_STUDY_P:g}_{st}.tsv"
        awk_p_below(gwas_path(st, registry, af_sources), dst, AUDIT_STUDY_P)
        return st, load_snv_extract(dst)

    with ThreadPoolExecutor(threads) as ex:
        extracts = dict(ex.map(scan, sorted(ref_of)))
    log(f"orientation audit: {len(extracts)} studies scanned at p < {AUDIT_STUDY_P:g}")
    rows = []
    for st, ref in sorted(ref_of.items()):
        n, agree = orientation_concordance(extracts[ref], extracts[st])
        c = controls[(controls.study == st) & (controls.p < AUDIT_REF_P) & controls.trait.isin(POSITIVE_CONTROL_TRAITS)]
        signs = [int(np.sign(b)) for b in c.beta_known_raising_allele_if_allele1]
        rows.append(dict(study=st, trait=trait[st], reference_study=ref, n_loci=n, share_same_sign=agree,
                         n_positive_controls_p_lt_5e8=len(c),
                         positive_controls_allele1_convention_sign=",".join(f"{s:+d}" for s in signs),
                         effect_allele_column_concordance_only=effect_allele_column(n, agree, []),
                         effect_allele_column=effect_allele_column(n, agree, signs)))
    return pd.DataFrame(rows)


def positive_control_rows(gw, studies, trait):
    """Known raising alleles in each scanned file, read under the allele1 convention."""
    pc = []
    for label, c, p, up, other in POSITIVE_CONTROLS:
        for st in sorted(studies):
            rows = gw.get((st, c, p))
            hit = None if rows is None else rows[((rows.allele1 == up) & (rows.allele2 == other)) |
                                                 ((rows.allele1 == other) & (rows.allele2 == up))]
            if hit is None or hit.empty:
                pc.append(dict(control=label, study=st, trait=trait.get(st), found=False))
                continue
            h = hit.iloc[0]
            pc.append(dict(control=label, study=st, trait=trait.get(st), found=True, allele1=h.allele1,
                           allele2=h.allele2, beta=h.beta_num, p=h.pval_num,
                           beta_known_raising_allele_if_allele1=h.beta_num * (1 if h.allele1 == up else -1)))
    return pd.DataFrame(pc, columns=["control", "study", "trait", "found", "allele1", "allele2", "beta", "p",
                                     "beta_known_raising_allele_if_allele1"])


# ----------------------------------------------------------------------------- per-pair evaluation

def resolve_hit1(g_rows, e_rows):
    """COLOC merge at the GWAS lead: lowest-p GWAS row whose alleles match a Broadaway row."""
    if g_rows is None or g_rows.empty:
        return None, None, "GWAS lead not in summary statistics", 0
    if e_rows is None or e_rows.empty:
        return None, None, "no Broadaway row for the gene at the GWAS lead", 0
    g_rows = g_rows.sort_values("pval_num", kind="stable")
    for rank, (_, g) in enumerate(g_rows.iterrows(), start=1):
        hit = e_rows[((e_rows.EA == g.allele1) & (e_rows.NEA == g.allele2)) |
                     ((e_rows.EA == g.allele2) & (e_rows.NEA == g.allele1))]
        if len(hit):
            return g, hit.iloc[0], "", rank
    return None, None, "GWAS and Broadaway alleles differ at the GWAS lead", 0


def evaluate_pair(row, gw, eq, lifted, panel, pop, gtex, effect_column="allele1", conditional=None,
                  convention="", spans=()):
    """All columns of one colocalized signal pair.

    effect_column: audited effect-allele column of the study's file ('allele1', 'allele2',
    'unresolved'). conditional: True/False/None, whether the eQTL SuSiE component's effect at
    hit2 has the sign of the fit's marginal there. convention: the study's measured allele-order
    convention ('allele1_is_alt', 'allele1_is_ref', or anything else for none). spans: excluded
    spans as (chrom, start, end, name), GRCh38."""
    c, p1 = row.hit1.split(":")
    c2, p2 = row.hit2.split(":")
    p1, p2 = int(p1), int(p2)
    if c2 != c:
        raise SystemExit(f"hit1 and hit2 on different chromosomes: {row}")
    same_lead = p1 == p2
    out = dict(study_superpop=pop, same_lead=same_lead, gwas_effect_allele_column=effect_column,
               gwas_file_allele_order_convention=convention,
               eqtl_conditional_sign_equals_marginal=conditional)
    g, e1, why1, rank = resolve_hit1(gw.get((row.study, c, p1)), eq.get((row.ensembl, c, p1)))
    out["hit1_gwas_row_rank_by_p"] = rank
    tech = [why1]          # technical blocks: the lead cannot be oriented
    if effect_column == "unresolved":
        tech.append("GWAS file effect-allele convention unresolved by the study audit")
    e2_rows = eq.get((row.ensembl, c, p2))
    e2 = None
    if same_lead and e1 is not None:
        e2 = e1
    elif e2_rows is None or e2_rows.empty:
        tech.append("eQTL lead not in Broadaway marginal file")
    elif e2_rows[["EA", "NEA"]].drop_duplicates().shape[0] > 1:
        tech.append("eQTL lead position is multi-allelic in Broadaway")
    else:
        e2 = e2_rows.iloc[0]
    ga = gb = None
    gaf = np.nan
    if g is not None:
        ga, gb, gaf = gwas_effect_alleles(g.allele1, g.allele2, g.af_num, effect_column)
        raising = ga if g.beta_num > 0 else gb if g.beta_num < 0 else ""     # NaN gives ""
        if effect_column == "unresolved":
            raising = ""
        out.update(gwas_allele1=g.allele1, gwas_allele2=g.allele2, gwas_beta=g.beta_num, gwas_p=g.pval_num,
                   gwas_af_allele1=g.af_num, gwas_effect_allele=ga, gwas_other_allele=gb,
                   gwas_af_effect_allele=gaf, trait_raising_allele=raising)
    if e2 is not None:
        out.update(eqtl_ea_hit2=e2.EA, eqtl_nea_hit2=e2.NEA, eqtl_beta_ea_hit2=e2.beta_num,
                   eqtl_p_hit2=e2.p_num, eqtl_eaf_hit2=e2.eaf_num)
    l1, l2 = lifted.get((c, p1)), lifted.get((c, p2))
    out["hit1_grch38"] = f"{l1[0]}:{l1[1]}:{l1[2]}" if l1 else ""
    out["hit2_grch38"] = f"{l2[0]}:{l2[1]}:{l2[2]}" if l2 else ""
    # excluded spans (MHC, IG, TR), as in the tag layer
    out["hit1_excluded_span"] = excluded_span(spans, l1[0], l1[1]) if l1 else ""
    out["hit2_excluded_span"] = excluded_span(spans, l2[0], l2[1]) if l2 else ""
    for lead in ("hit1", "hit2"):
        if out[f"{lead}_excluded_span"]:
            tech.append(f"{lead} in excluded span ({out[f'{lead}_excluded_span']})")
    x1 = y2 = None
    o1 = ""
    ld_note = ""
    if g is not None:
        if l1:
            x1, o1 = panel.indicator(l1[0], l1[1], l1[2], ga, gb)
            out["hit1_panel_match"] = o1.replace("a=", "gwas_effect=")
        else:
            out["hit1_panel_match"] = "not lifted to GRCh38"
    if e2 is not None:
        if l2:
            y2, o2 = panel.indicator(l2[0], l2[1], l2[2], e2.EA, e2.NEA)
            out["hit2_panel_match"] = o2.replace("a=", "EA=")
        else:
            out["hit2_panel_match"] = "not lifted to GRCh38"
    # palindromes (frequency of the GWAS effect allele in the study superpopulation; of EA in EUR)
    f1 = f2 = np.nan
    if g is not None:
        f1 = panel.freq(x1, pop) if x1 is not None else np.nan
        out["panel_freq_gwas_effect_allele_hit1"] = f1
        conv_ok = convention_followed(convention, effect_column, o1) if x1 is not None else None
        out["gwas_allele1_follows_file_convention"] = conv_ok
        tech.append(palindrome_block(ga, gb, gaf, f1, "GWAS lead", conv_ok))
    if e2 is not None:
        f2 = panel.freq(y2, "EUR") if y2 is not None else np.nan
        out["panel_eur_freq_ea_hit2"] = f2
        tech.append(palindrome_block(e2.EA, e2.NEA, e2.eaf_num, f2, "eQTL lead"))
    # different leads: the 1kGP records give r, so they must be the files' variants
    if not same_lead:
        if g is not None:
            src, src_from = gaf, "GWAS file"
            if not finite(gaf) and pop == "EUR" and e1 is not None:
                src = e1.eaf_num if ga == e1.EA else 1 - e1.eaf_num
                src_from = "Broadaway EAF (no GWAS file frequency)"
            out["hit1_source_freq_effect_allele"] = src
            out["hit1_source_freq_from"] = src_from if finite(src) else ""
            out["hit1_freq_diff_source_vs_1kgp"] = abs(src - f1) if finite(src) and finite(f1) else np.nan
            tech.append(frequency_block(ga, gb, src, f1, "GWAS lead"))
        if e2 is not None:
            e2f = e2.eaf_num
            out["hit2_freq_diff_source_vs_1kgp"] = abs(e2f - f2) if finite(e2f) and finite(f2) else np.nan
            tech.append(frequency_block(e2.EA, e2.NEA, e2f, f2, "eQTL lead"))
    r_study = r_eur = np.nan
    if x1 is not None and y2 is not None:
        r_study, r_eur = panel.r(x1, y2, pop), panel.r(x1, y2, "EUR")
        if not np.isfinite(r_study):
            ld_note = "monomorphic or too few haplotypes"
    elif not same_lead:
        ld_note = "; ".join(v for k, v in out.items() if k.endswith("_panel_match") and "=" not in v) or \
                  "lead row missing"
    out.update(ld_r_effect_ea_study=r_study, ld_r_effect_ea_eur=r_eur)
    r_used = np.nan
    if g is not None and e2 is not None:
        r_used, lblock = linkage(same_lead, ga, gb, e2.EA, e2.NEA, r_study, r_eur, pop, ld_note)
        tech.append(lblock)
    out["ld_r_used"] = r_used
    # evidence blocks: the marginal sign must be evidence of the colocalized signal's sign
    gwas_p = g.pval_num if g is not None else np.nan
    eqtl_p = e2.p_num if e2 is not None else np.nan
    evid = []
    if g is not None:
        evid.append(marginal_block(gwas_p, "hit1", "GWAS"))
    if e2 is not None:
        evid.append(marginal_block(eqtl_p, "hit2", "Broadaway"))
    evid.append(conditional_sign_block(conditional))
    beta_g = g.beta_num if g is not None else np.nan
    beta_e = e2.beta_num if e2 is not None else np.nan
    sign_tech, _ = direction_sign(beta_g, r_used, beta_e, tech)
    sign, reason = direction_sign(beta_g, r_used, beta_e, tech + evid)
    out.update(direction_sign=sign, not_directional_reason=reason, sign_before_evidence_rules=sign_tech,
               evidence_blocks="; ".join(b for b in evid if b),
               all_block_reasons="; ".join(b for b in tech + evid if b),
               marginal_floor=marginal_floor_label(gwas_p, eqtl_p),
               **{"marginals_both_p_lt_1e-5": bool(finite(gwas_p) and finite(eqtl_p)
                                                  and gwas_p < MARGINAL_P_SENSITIVITY
                                                  and eqtl_p < MARGINAL_P_SENSITIVITY)})
    if sign != 0:
        out["eqtl_beta_of_hit2_allele_on_trait_raising_haplotype"] = abs(e2.beta_num) * sign
    # display column: r between the trait-raising allele and the eQTL EA (the call above is the only orientation)
    if g is not None and effect_column != "unresolved" and np.isfinite(r_used) and np.isfinite(g.beta_num) \
            and g.beta_num != 0:
        out["ld_r_trait_raising_allele_vs_ea"] = r_used * np.sign(g.beta_num)
    # check 1 (LD step only; the GWAS beta sign cancels): Broadaway marginal of the trait-raising
    # allele at the GWAS lead itself
    if g is not None and e1 is not None and np.isfinite(g.beta_num) and g.beta_num != 0:
        same = 1 if (ga == e1.EA) else -1
        out["check_eqtl_at_hit1_sign"] = int(np.sign(g.beta_num) * same * np.sign(e1.beta_num))
        out["check_eqtl_at_hit1_p"] = e1.p_num
    # check 2: GTEx v8 liver slope (ALT, GRCh38) at the eQTL lead, same allele by letters
    if e2 is not None and l2:
        ea, nea = (e2.EA, e2.NEA) if l2[2] == "+" else (e2.EA.translate(COMPLEMENT)[::-1],
                                                         e2.NEA.translate(COMPLEMENT)[::-1])
        for ref, alt, slope in gtex.get((row.ensembl, l2[0], l2[1]), []):
            if {ref, alt} == {ea, nea}:
                gsign = np.sign(slope) * (1 if ea == alt else -1)
                out["check_gtex_liver_same_sign_hit2"] = int(gsign == np.sign(e2.beta_num))
                out["check_gtex_liver_slope_ea_hit2"] = slope * (1 if ea == alt else -1)
    return out, x1, y2


def eqtl_conditional_signs(pairs, fit_dir, rscript, work):
    """{(ensembl, idx2, hit2): True/False/None} from eqtl_conditional_sign.R, plus its table."""
    inp, outp = work / "eqtl_conditional_input.tsv", work / "eqtl_conditional_output.tsv"
    pairs[["ensembl", "chr", "idx2", "hit2"]].drop_duplicates().to_csv(inp, sep="\t", index=False)
    subprocess.run([rscript, str(HERE / "eqtl_conditional_sign.R"), str(inp), str(fit_dir), str(outp)],
                   check=True)
    t = read_str(outp)
    val = t.conditional_sign_equals_fitted_marginal.map({"TRUE": True, "FALSE": False})
    return {k: (None if pd.isna(v) else bool(v)) for k, v in zip(zip(t.ensembl, t.idx2, t.hit2), val)}, t


# ----------------------------------------------------------------------------- main

def load_tag_lists(genes):
    """{list name: tags (gene_key, chrom, pos_i, ref, alt, tag_is_lead)} for the COLOC genes, plus log lines."""
    t3 = read_str(TAG_LISTS["t3_tags_final"][0])
    ex = read_str(TAG_LISTS["t3prime_kept"][0])
    key = ["gene_id", "chrom", "pos", "ref", "alt"]
    if len(ex) != len(t3) or ex[key].drop_duplicates().shape[0] != t3[key].drop_duplicates().shape[0]:
        raise SystemExit("T3' tag_exclusions.tsv does not cover tags_final row for row")
    t3p = t3.merge(ex.loc[ex.kept == "True", key].drop_duplicates(), on=key)
    v2 = pd.concat([read_str(f) for f in sorted(Path(TAG_LISTS["t1v2"][0]).glob("chr*.tags.tsv"))], ignore_index=True)
    v2["gene_key"] = v2.gene_id.str.split(".").str[0]
    ex2 = read_str(TAG_LISTS["t3prime_v2_kept"][0])
    if len(ex2) != len(v2) or ex2[key].drop_duplicates().shape[0] != v2[key].drop_duplicates().shape[0] \
            or len(v2[key].drop_duplicates().merge(ex2[key].drop_duplicates(), on=key)) != v2[key].drop_duplicates().shape[0]:
        raise SystemExit("T3' v2 tag_exclusions_v2.tsv does not cover the T1 v2 tags row for row")
    v2p = v2.merge(ex2.loc[ex2.kept_steps_1_4 == "True", key].drop_duplicates(), on=key)
    lists = (("t3_tags_final", t3), ("t3prime_kept", t3p), ("t1v2", v2), ("t3prime_v2_kept", v2p))
    notes = [f"tag list {n}: {len(d)} rows, {d.gene_key.nunique()} genes (before restricting to COLOC genes)"
             for n, d in lists]
    out = {}
    for name, d in lists:
        d = d[d.gene_key.isin(genes)].assign(pos_i=lambda x: x.pos.astype(int))
        out[name] = d[["gene_key", "chrom", "pos_i", "ref", "alt", "tag_is_lead"]].drop_duplicates(
            ["gene_key", "chrom", "pos_i", "ref", "alt"])
    return out, notes


def tag_flags(tag_ld, units, keys, genes_with_tags, list_name, r2col, prefix):
    """Per unit: tags of `list_name` in r2 >= TAG_R2 (EUR) with any colocalized lead of the unit.

    A lead is not estimable when no tag of the list has a finite r with it (lead or tags
    absent from the 1kGP panel)."""
    t = tag_ld[tag_ld[f"in_{list_name}"]]
    per_lead = t.groupby(keys + ["hit2"])[r2col].max().reset_index()
    agg = per_lead.groupby(keys)[r2col].agg(
        n_leads="size", n_not_est=lambda s: int(s.isna().sum()), max_r2="max").reset_index()
    ntag = t.groupby(keys).tag_pos.nunique().rename("n_tags").reset_index()
    m = units[keys].drop_duplicates().merge(agg, on=keys, how="left").merge(ntag, on=keys, how="left")

    def flag(r):
        if r.ensembl not in genes_with_tags:
            return "no Model A tag"
        if r.max_r2 >= TAG_R2:
            return "yes"
        if pd.isna(r.n_leads) or r.n_not_est == r.n_leads:
            return "not estimable"
        if r.n_not_est > 0:
            return f"no among estimable ({int(r.n_not_est)} of {int(r.n_leads)} leads not estimable)"
        return "no"

    m[f"{prefix}tag_r2_ge_{TAG_R2}__{list_name}"] = [flag(r) for r in m.itertuples()]
    return m.rename(columns={"n_leads": f"{prefix}n_leads__{list_name}", "n_not_est": f"{prefix}n_leads_not_estimable__{list_name}",
                             "max_r2": f"{prefix}max_r2_eur__{list_name}", "n_tags": f"{prefix}n_tags__{list_name}"})


def blocked_pair_counts(pair_tab, calls, keys, sign_col):
    """Per unit: pairs blocked, and evidence-blocked pairs whose technical-only sign opposes / agrees with the call.

    n_pairs_component_blocked_opposing_call counts the opposing pairs blocked by the eQTL
    component rule alone (both marginals pass the floor); opposing_component_pair_block uses it."""
    d = pair_tab.merge(calls[keys + [sign_col]].rename(columns={sign_col: "_call"}), on=keys)
    ev = (d.direction_sign == 0) & (d.sign_before_evidence_rules != 0)
    opp = ev & (d._call != 0) & (d.sign_before_evidence_rules != d._call)
    d = d.assign(_blocked=d.direction_sign == 0, _ev_opp=opp,
                 _ev_agree=ev & (d._call != 0) & (d.sign_before_evidence_rules == d._call),
                 _comp_opp=opp & (d.evidence_blocks.fillna("") == COMPONENT_REASON))
    return d.groupby(keys).agg(n_pairs_blocked=("_blocked", "sum"),
                               n_pairs_evidence_blocked_opposing_call=("_ev_opp", "sum"),
                               n_pairs_evidence_blocked_agreeing_call=("_ev_agree", "sum"),
                               n_pairs_component_blocked_opposing_call=("_comp_opp", "sum")).reset_index()


def opposing_component_pair_block(tab):
    """Directional units with a component-only-blocked pair of the opposite sign become not_directional.

    The component rule says the pair's marginal sign may not be the colocalized component's sign;
    with both marginals past the floor that pair is unresolved evidence against the call, so the
    call is not kept (R-GEN review round 1)."""
    hit = (tab.direction_sign != 0) & (tab.n_pairs_component_blocked_opposing_call > 0)
    tab.loc[hit, "direction_sign"] = 0
    tab.loc[hit, "not_directional_reason"] = OPPOSING_COMPONENT_REASON
    if "eqtl_beta_of_hit2_allele_on_trait_raising_haplotype" in tab:
        tab.loc[hit, "eqtl_beta_of_hit2_allele_on_trait_raising_haplotype"] = np.nan
    return tab


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--work-dir", required=True, type=Path)
    ap.add_argument("--liftover-bin", required=True)
    ap.add_argument("--rscript", default="/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript")
    ap.add_argument("--eqtl-fit-dir", type=Path, default=FM / "results/eqtl_susie_polyfun")
    ap.add_argument("--studies", default="", help="comma list; smoke-test subset")
    ap.add_argument("--chroms", default="", help="comma list of hg19 chromosomes; smoke-test subset")
    ap.add_argument("--threads", type=int, default=8)
    a = ap.parse_args()
    if a.out_dir.exists():
        raise SystemExit(f"refusing to overwrite {a.out_dir}")
    a.out_dir.mkdir(parents=True)
    a.work_dir.mkdir(parents=True, exist_ok=True)
    t0 = datetime.now(timezone.utc)
    log = lambda *m: print(f"[{datetime.now(timezone.utc):%H:%M:%S}]", *m, flush=True)

    tiers = read_str(FM / "config/gwas_trait_tier.tsv")
    tier = dict(zip(tiers.study_name, tiers.tier.astype(int)))
    trait = dict(zip(tiers.study_name, tiers.trait))
    tier_label = dict(zip(tiers.study_name, tiers.tier_label))
    registry = read_str(FM / "config/gwas_registry.tsv")
    ancestry = dict(zip(registry.study_name, registry.ancestry))
    af_sources = read_str(FM / "config/gwas_af_sources.tsv")
    conv_tab = read_str(ALLELE_CONVENTIONS)
    convention = dict(zip(conv_tab.study, conv_tab.convention))
    spans = [(f[0], int(f[1]), int(f[2]), f[3]) for f in
             (line.rstrip("\n").split("\t") for line in open(EXCLUDED_SPANS)) if len(f) >= 4]
    excluded_samples = frozenset(Path(KING_EXCLUDED).read_text().split())

    studies_rows, pairs = [], []
    for name, spec in RELEASES.items():
        m, p = load_release(name, spec, tier)
        log(name, f"{len(m)} gene x study rows, {m.ensembl.nunique()} genes, {len(p)} pairs with PP.H4 > 0.5,"
                  f" pp4 reproduced by pairs {int(m.pp4_reproduced_by_pairs.sum())}/{len(m)}")
        studies_rows.append(m)
        pairs.append(p)
    study_rows, pairs = pd.concat(studies_rows, ignore_index=True), pd.concat(pairs, ignore_index=True)
    full_gene_counts = study_rows.groupby("release").ensembl.nunique().to_dict()
    full_gene_union = int(study_rows.ensembl.nunique())
    if a.studies:
        keep = a.studies.split(",")
        study_rows, pairs = study_rows[study_rows.study.isin(keep)], pairs[pairs.study.isin(keep)]
    if a.chroms:
        keep = a.chroms.split(",")
        study_rows, pairs = study_rows[study_rows.chr.isin(keep)], pairs[pairs.chr.isin(keep)]
    pairs = pairs.reset_index(drop=True)
    pairs["pair_source"] = pairs.release.map(PAIR_SOURCE)
    unmatched = study_rows.merge(pairs[["release", "study", "ensembl"]].drop_duplicates(), how="left",
                                 on=["release", "study", "ensembl"], indicator=True)
    no_pair = unmatched[unmatched._merge == "left_only"]
    log(f"evaluating {len(pairs)} pairs over {len(study_rows)} gene x study rows; {len(no_pair)} rows without a pair")

    pos = lambda s: (s.split(":")[0], int(s.split(":")[1]))
    # --- GWAS rows at hit1 and at the positive-control positions, for every audited study -------
    ref_of = audit_studies(set(pairs.study), trait)
    need_gwas = {st: {f"{c}:{p}" for _, c, p, _, _ in POSITIVE_CONTROLS} for st in ref_of}
    for r in pairs.itertuples():
        need_gwas[r.study].add(r.hit1)

    def scan_gwas(st):
        dst = a.work_dir / f"gwas_{st}.tsv"
        awk_filter(gwas_path(st, registry, af_sources), need_gwas[st], dst, '$1":"$2')
        d = read_str(dst)
        d.insert(0, "study", st)
        return d

    with ThreadPoolExecutor(a.threads) as ex:
        gw_all = pd.concat(list(ex.map(scan_gwas, sorted(need_gwas))), ignore_index=True)
    for col in ("allele1", "allele2"):
        gw_all[col] = gw_all[col].str.upper()
    gw_all["beta_num"], gw_all["pval_num"], gw_all["af_num"] = \
        to_num(gw_all.beta), to_num(gw_all.pval), to_num(gw_all.af)
    gw = {k: v for k, v in gw_all.groupby(["study", "chromosome", gw_all.position.astype(int)])}
    log(f"GWAS rows at requested positions: {len(gw_all)}")

    # --- GWAS orientation audit ------------------------------------------------------------------
    pc = positive_control_rows(gw, ref_of, trait)
    audit = run_orientation_audit(ref_of, trait, pc, registry, af_sources, a.work_dir, a.threads, log)
    audit["file_allele_order_convention"] = audit.study.map(convention).fillna("not measured")
    effect_col = dict(zip(audit.study, audit.effect_allele_column))
    log("GWAS effect-allele column: " + ", ".join(f"{s} {c}" for s, c in sorted(effect_col.items())))
    audit.to_csv(a.out_dir / "gwas_orientation_audit.tsv", sep="\t", index=False)
    pc["effect_allele_column"] = pc.study.map(effect_col)
    pc["beta_known_raising_allele_audited"] = np.where(pc.effect_allele_column == "allele2",
                                                       -pc.beta_known_raising_allele_if_allele1,
                                                       pc.beta_known_raising_allele_if_allele1)
    pc.loc[pc.effect_allele_column == "unresolved", "beta_known_raising_allele_audited"] = np.nan
    pc.to_csv(a.out_dir / "positive_controls_gwas.tsv", sep="\t", index=False)

    # --- Broadaway rows at hit1 and hit2 --------------------------------------------------------
    need_eq = {}
    for r in pairs.itertuples():
        for h in (r.hit1, r.hit2):
            c, p = pos(h)
            need_eq.setdefault(c, set()).add(f"{r.ensembl}:{p}")

    def scan_eqtl(c):
        dst = a.work_dir / f"eqtl_chr{c}.tsv"
        awk_filter(REPO / f"data/broadaway_eqtl/chr{c}_marginal_summary_results.tsv", need_eq[c], dst, '$14":"$4')
        return read_str(dst)[EQTL_COLS]

    with ThreadPoolExecutor(a.threads) as ex:
        eq_all = pd.concat(list(ex.map(scan_eqtl, sorted(need_eq, key=int))), ignore_index=True)
    for col in ("EA", "NEA"):
        eq_all[col] = eq_all[col].str.upper()
    eq_all["beta_num"], eq_all["p_num"], eq_all["eaf_num"] = to_num(eq_all.Beta), to_num(eq_all.PVAL), to_num(eq_all.EAF)
    eq = {k: v for k, v in eq_all.groupby(["ENSG", "CHR", eq_all.POS.astype(int)])}
    log(f"Broadaway rows at requested gene x positions: {len(eq_all)}")

    # --- eQTL SuSiE component sign at hit2 ------------------------------------------------------
    cond, cond_tab = eqtl_conditional_signs(pairs, a.eqtl_fit_dir, a.rscript, a.work_dir)
    log(f"eQTL SuSiE component signs read for {len(cond_tab)} (gene, component, lead) rows")

    # --- liftover, tags, 1kGP haplotypes -------------------------------------------------------
    hg19 = {pos(h) for h in pd.concat([pairs.hit1, pairs.hit2])}
    lifted = liftover(hg19, REPO / "data/broadaway_eqtl/hg19ToHg38.over.chain", a.liftover_bin, a.work_dir)
    log(f"lifted {len(lifted)} of {len(hg19)} lead positions to the same GRCh38 chromosome")
    tag_lists, tag_notes = load_tag_lists(set(pairs.ensembl) | set(study_rows.ensembl))
    for n in tag_notes:
        log(n)
    tkey = ["gene_key", "chrom", "pos_i", "ref", "alt"]
    tags = pd.concat(tag_lists.values(), ignore_index=True).drop_duplicates(tkey).reset_index(drop=True)
    for name, d in tag_lists.items():
        tags[f"in_{name}"] = tags[tkey].merge(d[tkey].assign(_x=True), on=tkey, how="left")._x.notna().values
    by_chrom = {}
    for chrom, p, _ in lifted.values():
        by_chrom.setdefault(chrom, set()).add(p)
    for t in tags.itertuples():
        by_chrom.setdefault(t.chrom, set()).add(t.pos_i)
    samples, sample_pop = [], {}
    for sp in SUPERPOPS:
        for s in founders(KGP / "20130606_g1k_3202_samples_ped_population.txt", sp):
            samples.append(s)
            sample_pop[s] = sp
    with ThreadPoolExecutor(a.threads) as ex:
        chroms = sorted(by_chrom)
        haps = dict(zip(chroms, ex.map(lambda ch: panel_haplotypes(ch, by_chrom[ch], samples, a.work_dir), chroms)))
    panel = Panel(haps, samples, sample_pop, excluded_samples)
    log("1kGP samples used (founders minus KING-related): " + ", ".join(f"{k} {v}" for k, v in panel.n_samples.items()))

    gtex = {}
    gt = pd.read_csv(REPO / "data/external/allelic_refs/gtex_v8/GTEx_Analysis_v8_eQTL/Liver.v8.signif_variant_gene_pairs.txt.gz",
                     sep="\t", usecols=["variant_id", "gene_id", "slope"], dtype={"variant_id": str, "gene_id": str})
    gt["ensg"] = gt.gene_id.str.split(".").str[0]
    gt = gt[gt.ensg.isin(set(pairs.ensembl))]
    for r in gt.itertuples():
        chrom, p, ref, alt, _ = r.variant_id.split("_")
        gtex.setdefault((r.ensg, chrom, int(p)), []).append((ref, alt, r.slope))

    # --- evaluate every colocalized signal pair ------------------------------------------------
    recs, lead_haps = [], {}
    for r in pairs.itertuples():
        pop = ancestry.get(r.study, "EUR")
        res, x1, y2 = evaluate_pair(r, gw, eq, lifted, panel, pop, gtex, effect_col[r.study],
                                    cond.get((r.ensembl, r.idx2, r.hit2)), convention.get(r.study, ""), spans)
        recs.append(res)
        if y2 is not None:
            lead_haps[(r.release, r.study, r.ensembl, r.idx1, r.idx2)] = y2
    pair_tab = pd.concat([pairs, pd.DataFrame(recs)], axis=1)
    pair_tab = pair_tab.merge(cond_tab.rename(columns=lambda c: c if c in ("ensembl", "idx2", "hit2") else f"eqtl_susie_{c}"),
                              on=["ensembl", "idx2", "hit2"], how="left")
    pair_tab["trait"] = pair_tab.study.map(trait)
    pair_tab["tier"] = pair_tab.study.map(tier)
    pair_tab["tier_label"] = pair_tab.study.map(tier_label)
    pair_tab["direction"] = [direction_label(t, s) for t, s in zip(pair_tab.trait, pair_tab.direction_sign)]
    pair_tab["provenance"] = PROVENANCE

    # --- gene x release x study: agreed sign of the colocalized signal pairs --------------------
    study_tab = combine_rows(pair_tab, ["release", "study", "ensembl"], "pair_pp4", ["idx1", "idx2"],
                             "colocalized signal pairs of this study") \
        .rename(columns={"n_total": "n_pairs_pp4_gt_05", "n_up": "n_pairs_up", "n_down": "n_pairs_down"})
    study_tab["reported_pair_is_highest_pp4"] = study_tab.is_highest_pp4_pair
    study_tab = study_tab.merge(blocked_pair_counts(pair_tab, study_tab, ["release", "study", "ensembl"], "direction_sign"),
                                on=["release", "study", "ensembl"], how="left")
    study_tab = opposing_component_pair_block(study_tab)
    if len(no_pair):
        extra = no_pair.drop(columns="_merge").assign(not_directional_reason="no signal-pair record", direction_sign=0,
                                                      pair_source=lambda x: x.release.map(PAIR_SOURCE))
        study_tab = pd.concat([study_tab, extra], ignore_index=True)
        study_tab["trait"] = study_tab.study.map(trait)
    study_tab["direction"] = [direction_label(t, s) for t, s in zip(study_tab.trait, study_tab.direction_sign)]
    study_tab["provenance"] = PROVENANCE

    # --- Model A tags in LD with the colocalized Broadaway leads (EUR) --------------------------
    tag_rows = []
    tag_haps = {t.Index: panel.indicator(t.chrom, t.pos_i, "+", t.alt, t.ref) for t in tags.itertuples()}
    in_cols = [f"in_{n}" for n in TAG_LISTS]
    for r in pair_tab.itertuples():
        y = lead_haps.get((r.release, r.study, r.ensembl, r.idx1, r.idx2))
        for t in tags[tags.gene_key == r.ensembl].itertuples():
            x, o = tag_haps[t.Index]
            ok = x is not None and y is not None
            rr = panel.r(x, y, "EUR") if ok else np.nan
            ra = panel.r(x, y, "EUR_all_founders") if ok else np.nan
            tag_rows.append(dict(release=r.release, ensembl=r.ensembl, gene=r.gene, study=r.study, trait=r.trait,
                                 idx1=r.idx1, idx2=r.idx2, pair_pp4=r.pair_pp4, hit2=r.hit2, hit2_grch38=r.hit2_grch38,
                                 tag_chrom=t.chrom, tag_pos=t.pos_i, tag_ref=t.ref, tag_alt=t.alt,
                                 tag_is_gtex_lead=t.tag_is_lead, **{c: getattr(t, c) for c in in_cols},
                                 tag_panel_match=o, r_eur_tag_alt_vs_hit2_ea=rr,
                                 r2_eur=rr ** 2 if np.isfinite(rr) else np.nan,
                                 r2_eur_all_founders=ra ** 2 if np.isfinite(ra) else np.nan))
    tag_ld = pd.DataFrame(tag_rows, columns=["release", "ensembl", "gene", "study", "trait", "idx1", "idx2", "pair_pp4",
                                             "hit2", "hit2_grch38", "tag_chrom", "tag_pos", "tag_ref", "tag_alt",
                                             "tag_is_gtex_lead"] + in_cols +
                                            ["tag_panel_match", "r_eur_tag_alt_vs_hit2_ea", "r2_eur", "r2_eur_all_founders"])
    genes_rel = study_tab[["release", "ensembl", "gene"]].drop_duplicates()
    genes_with_tags = {n: set(d.gene_key) for n, d in tag_lists.items()}
    trait_units = study_tab[["release", "ensembl", "trait"]].drop_duplicates()
    sens = {}
    for n in TAG_LISTS:
        genes_rel = genes_rel.merge(tag_flags(tag_ld, genes_rel, ["release", "ensembl"], genes_with_tags[n], n,
                                              "r2_eur", "gene_"), on=["release", "ensembl"], how="left")
        alt = tag_flags(tag_ld, genes_rel, ["release", "ensembl"], genes_with_tags[n], n, "r2_eur_all_founders", "gene_")
        col = f"gene_tag_r2_ge_{TAG_R2}__{n}"
        cmp_ = genes_rel[["release", "ensembl", col]].merge(alt[["release", "ensembl", col]], on=["release", "ensembl"],
                                                             suffixes=("", "_all"))
        mx = genes_rel[f"gene_max_r2_eur__{n}"]
        sens[n] = {"genes_flag_differs_with_all_525_eur_founders": int((cmp_[col] != cmp_[col + "_all"]).sum()),
                   f"genes_best_tag_r2_in_[{TAG_R2_BAND[0]}, {TAG_R2_BAND[1]})": {
                       rel: int(((g >= TAG_R2_BAND[0]) & (g < TAG_R2_BAND[1])).sum())
                       for rel, g in mx.groupby(genes_rel.release)}}
        trait_units = trait_units.merge(tag_flags(tag_ld, trait_units, ["release", "ensembl", "trait"], genes_with_tags[n],
                                                  n, "r2_eur", "trait_"), on=["release", "ensembl", "trait"], how="left")

    # --- gene x release x trait: agreed sign of the studies of the trait -----------------------
    trait_tab = combine_rows(study_tab.rename(columns={"pp4": "study_pp4_susie"}), ["release", "ensembl", "trait"],
                             "study_pp4_susie", ["study"], "studies of this trait") \
        .rename(columns={"study": "reported_study", "study_pp4_susie": "reported_study_pp4_susie",
                         "n_total": "n_studies", "n_up": "n_studies_up", "n_down": "n_studies_down"})
    top_pp4 = study_tab.sort_values(["pp4", "study"], ascending=[False, True]).drop_duplicates(["release", "ensembl", "trait"])
    trait_tab = trait_tab.merge(top_pp4[["release", "ensembl", "trait", "study", "not_directional_reason"]]
                                .rename(columns={"study": "highest_pp4_study",
                                                 "not_directional_reason": "highest_pp4_study_not_directional_reason"}),
                                on=["release", "ensembl", "trait"], how="left")
    trait_tab = trait_tab.drop(columns=["n_pairs_blocked", "n_pairs_evidence_blocked_opposing_call",
                                        "n_pairs_evidence_blocked_agreeing_call",
                                        "n_pairs_component_blocked_opposing_call"], errors="ignore")
    trait_tab = trait_tab.merge(blocked_pair_counts(pair_tab, trait_tab, ["release", "ensembl", "trait"], "direction_sign"),
                                on=["release", "ensembl", "trait"], how="left")
    trait_tab = opposing_component_pair_block(trait_tab)
    trait_tab["direction"] = [direction_label(t, s) for t, s in zip(trait_tab.trait, trait_tab.direction_sign)]
    trait_tab = trait_tab.merge(pair_tab.groupby(["release", "ensembl", "trait"]).size().rename("n_pairs_trait").reset_index(),
                                on=["release", "ensembl", "trait"], how="left")
    trait_tab = trait_tab.merge(trait_units, on=["release", "ensembl", "trait"], how="left")
    trait_tab["provenance"] = PROVENANCE

    flag_cols = [f"trait_tag_r2_ge_{TAG_R2}__{n}" for n in TAG_LISTS]
    first = ["release", "gene", "ensembl", "trait", "tier_label", "reported_study", "reported_study_pp4_susie",
             "study_superpop", "hit1", "gwas_allele1", "gwas_allele2", "gwas_beta", "gwas_p",
             "gwas_effect_allele_column", "gwas_effect_allele", "trait_raising_allele", "hit2", "eqtl_ea_hit2",
             "eqtl_nea_hit2", "eqtl_beta_ea_hit2", "eqtl_p_hit2", "hit1_panel_match", "hit2_panel_match", "same_lead",
             "ld_r_effect_ea_study", "ld_r_effect_ea_eur", "ld_r_used", "ld_r_trait_raising_allele_vs_ea",
             "marginal_floor", "marginals_both_p_lt_1e-5", "eqtl_conditional_sign_equals_marginal",
             "eqtl_susie_lbf_argmax_is_hit2", "eqtl_susie_alpha_hit2", "hit1_excluded_span", "hit2_excluded_span",
             "direction", "direction_sign", "not_directional_reason", "all_block_reasons", "sign_before_evidence_rules",
             "eqtl_beta_of_hit2_allele_on_trait_raising_haplotype", "n_studies", "n_studies_up", "n_studies_down",
             "n_pairs_trait", "n_pairs_blocked", "n_pairs_evidence_blocked_opposing_call",
             "n_pairs_evidence_blocked_agreeing_call", "n_pairs_component_blocked_opposing_call",
             "highest_pp4_study", "highest_pp4_study_not_directional_reason",
             "check_eqtl_at_hit1_sign", "check_gtex_liver_same_sign_hit2"] + flag_cols + ["pair_source", "provenance"]
    order = lambda d, lead: d[[c for c in lead if c in d.columns] + [c for c in d.columns if c not in lead]]
    trait_tab = order(trait_tab, first).sort_values(["release", "gene", "trait"])
    trait_tab.to_csv(a.out_dir / "direction_by_trait.tsv", sep="\t", index=False)
    order(study_tab.rename(columns={"pp4": "study_pp4_susie"}), ["release", "gene", "ensembl", "study", "trait"] + first[7:]) \
        .sort_values(["release", "gene", "study"]).to_csv(a.out_dir / "direction_by_study.tsv", sep="\t", index=False)
    order(pair_tab, ["release", "gene", "ensembl", "study", "trait", "idx1", "idx2", "pair_pp4", "is_highest_pp4_pair"]) \
        .sort_values(["release", "gene", "study", "idx1", "idx2"]).to_csv(a.out_dir / "direction_signal_pairs.tsv", sep="\t", index=False)
    tag_ld.to_csv(a.out_dir / "tag_ld_to_coloc_lead.tsv", sep="\t", index=False)
    genes_rel["ld_reference_eur"] = f"1kGP EUR founders minus KING-related, n = {panel.n_samples['EUR']} (as T1 v2)"
    genes_rel.to_csv(a.out_dir / "tag_ld_gene_summary.tsv", sep="\t", index=False)

    # --- summary --------------------------------------------------------------------------------
    norm = lambda s: s.str.replace(r" in (EUR|EAS|AFR|AMR|SAS) ", " in <pop> ", regex=True)

    def counts(d):
        return {rel: {"rows": int(len(g)), "genes": int(g.ensembl.nunique()),
                      "increases": int((g.direction_sign > 0).sum()), "decreases": int((g.direction_sign < 0).sum()),
                      "not_directional": int((g.direction_sign == 0).sum()),
                      "not_directional_reasons": norm(g.loc[g.direction_sign == 0, "not_directional_reason"])
                      .value_counts().to_dict()}
                for rel, g in d.groupby("release")}

    def by_rel(d, f):
        return {rel: f(g) for rel, g in d.groupby("release")}

    directional = trait_tab[trait_tab.direction_sign != 0]
    nonpal = pair_tab[pair_tab.gwas_effect_allele.notna() & pair_tab.gwas_af_effect_allele.notna()
                      & pair_tab.panel_freq_gwas_effect_allele_hit1.notna()]
    nonpal = nonpal[[not is_palindromic(x, y) for x, y in zip(nonpal.gwas_effect_allele, nonpal.gwas_other_allele)]]
    nonpal = nonpal[(nonpal.panel_freq_gwas_effect_allele_hit1 - 0.5).abs() > 0.1]     # informative frequencies only
    has_g = pair_tab.gwas_effect_allele.notna()
    is_pal1 = np.array([is_palindromic(x, y) if isinstance(x, str) and isinstance(y, str) else False
                        for x, y in zip(pair_tab.gwas_effect_allele, pair_tab.gwas_other_allele)], dtype=bool)
    pal1 = pair_tab[has_g.values & is_pal1]
    pal1_reason = pal1.all_block_reasons.str.extract(r"(palindromic GWAS lead[^;]*)")[0]
    conv_studies = pair_tab.gwas_file_allele_order_convention.isin(["allele1_is_alt", "allele1_is_ref"])
    np_conv = pair_tab[(conv_studies & has_g & pair_tab.gwas_allele1_follows_file_convention.notna()).values & ~is_pal1]
    first_reason_is = lambda pat: pair_tab.not_directional_reason.str.contains(pat, regex=False)
    any_reason_is = lambda pat: pair_tab.all_block_reasons.str.contains(pat, regex=False)
    # evidence-blocked pairs inside a directional gene x trait: does the technical-only sign agree with the call?
    ev = pair_tab.merge(trait_tab[["release", "ensembl", "trait", "direction_sign"]].rename(columns={"direction_sign": "call"}),
                        on=["release", "ensembl", "trait"])
    ev = ev[(ev.direction_sign == 0) & (ev.sign_before_evidence_rules != 0) & (ev.call != 0)]

    def agree_split(d):
        comp_only = d.evidence_blocks == COMPONENT_REASON
        out = {}
        for lab, s in (("component rule only", d[comp_only]), ("marginal floor (any)", d[d.evidence_blocks.str.contains("marginal p")])):
            strong = s.eqtl_p_hit2 < 1e-8
            out[lab] = {"n": int(len(s)), "agree_with_call": int((s.sign_before_evidence_rules == s.call).sum()),
                        "broadaway_p_lt_1e-8": f"{int((s[strong].sign_before_evidence_rules == s[strong].call).sum())}/{int(strong.sum())}",
                        "broadaway_p_ge_1e-8": f"{int((s[~strong].sign_before_evidence_rules == s[~strong].call).sum())}/{int((~strong).sum())}"}
        return out

    top_differs = trait_tab[trait_tab.reported_study != trait_tab.highest_pp4_study]
    diff_lead = pair_tab[~pair_tab.same_lead]
    ctrl_studies = audit[audit.n_positive_controls_p_lt_5e8 > 0]
    summary = {
        "created_utc": t0.strftime("%Y-%m-%dT%H:%M:%SZ"),
        "subset": {"studies": a.studies, "chroms": a.chroms},
        "genes_per_release_full": full_gene_counts,
        "genes_union_full": full_gene_union,
        "pair_source": {**PAIR_SOURCE,
                        "adopted_pairs_ineligible_under_corrected_rule_included_by_design":
                            by_rel(pair_tab, lambda g: f"{int((g.eligible_corrected_rule == 'FALSE').sum())}/{len(g)}")},
        "gene_study_rows_without_pair": int(len(no_pair)),
        "pp4_reproduced_by_pairs (max pair PP.H4 vs gene x study PP.H4, not lead identity)": {
            rel: f"{int(g.pp4_reproduced_by_pairs.sum())}/{len(g)}" for rel, g in study_rows.groupby("release")},
        "gwas_effect_allele_column": audit.set_index("study")[["reference_study", "n_loci", "share_same_sign",
                                                               "positive_controls_allele1_convention_sign",
                                                               "effect_allele_column_concordance_only",
                                                               "effect_allele_column", "file_allele_order_convention"]].to_dict("index"),
        "gwas_audit_controls": {
            "studies_resolved": int((audit.effect_allele_column != "unresolved").sum()),
            "studies_with_a_positive_control_at_p_lt_5e-8": int(len(ctrl_studies)),
            "studies_made_unresolved_by_a_control_conflict": int(((audit.effect_allele_column == "unresolved")
                                                                  & (audit.effect_allele_column_concordance_only != "unresolved")).sum()),
            "note": "controls read positive in every resolved study by construction (a conflict makes the study unresolved)"},
        "signal_pairs": counts(pair_tab), "by_study": counts(study_tab), "by_trait": counts(trait_tab),
        "rules_prespecified_2026-09-29": {
            "marginal_floor": f"GWAS p < {MARGINAL_P_MAX:g} at hit1 and Broadaway p < {MARGINAL_P_MAX:g} at hit2",
            "frequency_concordance": f"|source - 1kGP frequency| <= {FREQ_MAX_DIFF} for non-palindromic leads of different-lead pairs",
            "excluded_spans": str(EXCLUDED_SPANS),
            "palindrome_without_file_frequency": "allele1 must follow the study's measured allele-order convention",
            "opposing_component_pair (added 2026-09-30, before rerun v5)":
                "a pair blocked only by the eQTL component rule whose technical-only sign opposes the agreed sign "
                "makes the gene x study and gene x trait not_directional"},
        "pairs_blocked_by_rule": by_rel(pair_tab, lambda g: {
            lab: {"first_reason": int(g.not_directional_reason.str.contains(pat, regex=False).sum()),
                  "any_reason": int(g.all_block_reasons.str.contains(pat, regex=False).sum()),
                  "only_reason": int((g.all_block_reasons == g.not_directional_reason).mul(
                      g.not_directional_reason.str.contains(pat, regex=False)).sum())}
            for lab, pat in (("GWAS marginal floor", "GWAS marginal p"), ("Broadaway marginal floor", "Broadaway marginal p"),
                             ("eQTL SuSiE component", "SuSiE component"), ("excluded span", "excluded span"),
                             ("frequency concordance", "1kGP record may not be the file variant"),
                             ("palindrome", "palindromic"))}),
        "frequency_concordance_blocks_any_reason": {
            "GWAS lead": int(any_reason_is("GWAS lead: source and 1kGP").sum()),
            "eQTL lead": int(any_reason_is("eQTL lead: source and 1kGP").sum()),
            "of which indel leads": int((any_reason_is("1kGP record may not be the file variant") &
                                         ((pair_tab.gwas_allele1.str.len() > 1) | (pair_tab.gwas_allele2.str.len() > 1) |
                                          (pair_tab.eqtl_ea_hit2.str.len() > 1) | (pair_tab.eqtl_nea_hit2.str.len() > 1))).sum())},
        "palindromic_gwas_leads": {
            "pairs": int(len(pal1)),
            "passed_on_file_frequency": int((pal1_reason.isna() & pal1.gwas_af_effect_allele.notna()).sum()),
            "passed_on_allele_order_convention (no file frequency)": int((pal1_reason.isna() & pal1.gwas_af_effect_allele.isna()).sum()),
            "blocked": norm(pal1_reason.dropna()).value_counts().to_dict()},
        "allele_order_convention_check_nonpalindromic_hit1": {
            "pairs": int(len(np_conv)), "allele1_follows_convention": int((np_conv.gwas_allele1_follows_file_convention == True).sum()),
            **{f"{lab}_leads": f"{int((s.gwas_allele1_follows_file_convention == True).sum())}/{int(len(s))}"
               for lab, s in (("snv", np_conv[np_conv.gwas_allele1.str.len().eq(1) & np_conv.gwas_allele2.str.len().eq(1)]),
                              ("indel", np_conv[np_conv.gwas_allele1.str.len().gt(1) | np_conv.gwas_allele2.str.len().gt(1)]))}},
        "excluded_span_trait_rows": by_rel(trait_tab, lambda g: int(((g.hit1_excluded_span.fillna("") != "") |
                                                                    (g.hit2_excluded_span.fillna("") != "")).sum())),
        "evidence_blocked_pairs_in_directional_gene_trait_rows": agree_split(ev),
        "trait_rows_with_evidence_blocked_pair_opposing_call": [
            f"{r.release} {r.gene} {r.trait} call {int(r.direction_sign):+d}" for r in
            trait_tab[trait_tab.n_pairs_evidence_blocked_opposing_call > 0].itertuples()],
        "directional_trait_rows_marginal_strength": by_rel(directional, lambda g: {
            "rows": int(len(g)),
            f"both_marginals_p_lt_{MARGINAL_P_MAX:g} (all, by rule)": int((g.marginal_floor == "pass").sum()),
            "both_marginals_p_lt_1e-5 (sensitivity)": int(g["marginals_both_p_lt_1e-5"].sum()),
            "gwas_p_at_hit1_ge_5e-8": int((g.gwas_p >= 5e-8).sum()),
            "hit2_not_component_lbf_argmax": int((g.eqtl_susie_lbf_argmax_is_hit2 == "FALSE").sum())}),
        "combine_rule_audit": by_rel(trait_tab, lambda g: {
            "trait_rows": int(len(g)),
            "reported_study_is_not_highest_pp4_study": int((g.reported_study != g.highest_pp4_study).sum()),
            "of_which_directional": int(((g.reported_study != g.highest_pp4_study) & (g.direction_sign != 0)).sum()),
            "highest_pp4_study_blocked_because": norm(g.loc[(g.reported_study != g.highest_pp4_study) & (g.direction_sign != 0),
                                                              "highest_pp4_study_not_directional_reason"]).value_counts().to_dict(),
            "studies_disagree": int((g.not_directional_reason == "studies of this trait disagree").sum())}),
        "hit1_panel_match": pair_tab.hit1_panel_match.value_counts(dropna=False).to_dict(),
        "hit2_panel_match": pair_tab.hit2_panel_match.value_counts(dropna=False).to_dict(),
        "same_lead_pairs": int(pair_tab.same_lead.sum()),
        "self_consistency_same_lead_panel_r_sign_equals_letter_sign (cannot fail: same panel record)": same_lead_check(pair_tab),
        "hit1_row_not_lowest_p": int((pair_tab.hit1_gwas_row_rank_by_p > 1).sum()),
        "gwas_af_nearer_panel_freq_than_its_complement_nonpalindromic_hit1_panel_maf_lt_0.4": {
            "n": int(len(nonpal)),
            "yes": int(((nonpal.gwas_af_effect_allele - nonpal.panel_freq_gwas_effect_allele_hit1).abs()
                        < (nonpal.gwas_af_effect_allele - (1 - nonpal.panel_freq_gwas_effect_allele_hit1)).abs()).sum())},
        "eqtl_susie_component_check": {
            "pairs": int(len(pair_tab)),
            "fit_found": int((pair_tab.eqtl_susie_fit_found == "TRUE").sum()),
            "hit2_is_component_lbf_argmax": int((pair_tab.eqtl_susie_lbf_argmax_is_hit2 == "TRUE").sum()),
            "conditional_sign_equals_marginal": pair_tab.eqtl_conditional_sign_equals_marginal
            .map({True: "TRUE", False: "FALSE"}).fillna("not available").value_counts().to_dict()},
        "ld_step_check_broadaway_at_hit1_agrees_with_call_different_lead_pairs (GWAS sign cancels; tests LD step only)": {
            "n": int(((diff_lead.direction_sign != 0) & diff_lead.check_eqtl_at_hit1_sign.notna()).sum()),
            "agree": int(((diff_lead.direction_sign != 0) & (diff_lead.check_eqtl_at_hit1_sign == diff_lead.direction_sign)).sum())},
        "check_gtex_liver_same_sign_as_broadaway_at_hit2 (independent cohort)": {
            "n": int(pair_tab.get("check_gtex_liver_same_sign_hit2", pd.Series(dtype=float)).notna().sum()),
            "same_sign": int((pair_tab.get("check_gtex_liver_same_sign_hit2", pd.Series(dtype=float)) == 1).sum())},
        "tag_ld_reference": {"EUR samples": panel.n_samples["EUR"], "EUR_all_founders": panel.n_samples["EUR_all_founders"],
                             "definition": "1kGP founders (no recorded parent, trio parents kept) minus KING-related "
                                           f"samples in {KING_EXCLUDED}; the EUR set T1 v2 uses to select tags"},
        "tag_lists": {n: d for n, (_, d) in TAG_LISTS.items()},
        "tag_genes_gene_level": {n: genes_rel[f"gene_tag_r2_ge_{TAG_R2}__{n}"].groupby(genes_rel.release).value_counts()
                                 .unstack(fill_value=0).to_dict("index") for n in TAG_LISTS},
        "tag_flag_sensitivity": sens,
        "positive_controls_audited_raising_allele_positive_where_p_lt_5e-8 (by construction)": {
            lab: f"{int((g.beta_known_raising_allele_audited > 0).sum())}/{len(g)}"
            for lab, g in pc[(pd.to_numeric(pc.p, errors="coerce") < 5e-8) & pc.trait.isin(POSITIVE_CONTROL_TRAITS)
                             & (pc.effect_allele_column != "unresolved")].groupby("control")},
        "lifted_positions": f"{len(lifted)}/{len(hg19)}",
        "releases_agree_on_shared_gene_trait_rows": release_agreement(trait_tab),
    }
    (a.out_dir / "summary.json").write_text(json.dumps(summary, indent=2, default=str))
    params = dict(PP4_MIN=PP4_MIN, MIN_R2=MIN_R2, MAX_PAL_MAF=MAX_PAL_MAF, PAL_MARGIN=PAL_MARGIN,
                  PAL_MAX_DIST=PAL_MAX_DIST, TAG_R2=TAG_R2, MIN_HAPS=MIN_HAPS, SUPERPOPS=SUPERPOPS,
                  MARGINAL_P_MAX=MARGINAL_P_MAX, MARGINAL_P_SENSITIVITY=MARGINAL_P_SENSITIVITY,
                  FREQ_MAX_DIFF=FREQ_MAX_DIFF, TAG_R2_BAND=TAG_R2_BAND,
                  AUDIT_REF_P=AUDIT_REF_P, AUDIT_STUDY_P=AUDIT_STUDY_P, AUDIT_MIN_LOCI=AUDIT_MIN_LOCI,
                  AUDIT_AGREE=AUDIT_AGREE, AUDIT_WINDOW=AUDIT_WINDOW,
                  AUDIT_REFERENCE_BY_TRAIT=AUDIT_REFERENCE_BY_TRAIT,
                  AUDIT_REFERENCE_OF_REFERENCE=AUDIT_REFERENCE_OF_REFERENCE,
                  POSITIVE_CONTROLS=POSITIVE_CONTROLS, POSITIVE_CONTROL_TRAITS=sorted(POSITIVE_CONTROL_TRAITS),
                  EXCLUDED_SPANS=str(EXCLUDED_SPANS), ALLELE_CONVENTIONS=str(ALLELE_CONVENTIONS),
                  KING_EXCLUDED=str(KING_EXCLUDED), TAG_LISTS={k: [str(v[0]), v[1]] for k, v in TAG_LISTS.items()},
                  PAIR_SOURCE=PAIR_SOURCE, eqtl_fit_dir=str(a.eqtl_fit_dir),
                  releases={k: {kk: str(vv) for kk, vv in v.items()} for k, v in RELEASES.items()},
                  argv=sys.argv)
    (a.out_dir / "params.json").write_text(json.dumps(params, indent=2))
    manifest = [("script", Path(__file__).resolve()), ("eqtl_conditional_sign", HERE / "eqtl_conditional_sign.R")] + \
        [(f"{k}_{kk}", v[kk]) for k, v in RELEASES.items() for kk in ("master", "pairs")] + \
        [("tiers", FM / "config/gwas_trait_tier.tsv"), ("registry", FM / "config/gwas_registry.tsv"),
         ("af_sources", FM / "config/gwas_af_sources.tsv"), ("excluded_spans", EXCLUDED_SPANS),
         ("allele_conventions", ALLELE_CONVENTIONS), ("king_excluded", KING_EXCLUDED),
         ("tags_t3_tags_final", TAG_LISTS["t3_tags_final"][0]), ("tags_t3prime_exclusions", TAG_LISTS["t3prime_kept"][0]),
         ("tags_t3prime_v2_exclusions", TAG_LISTS["t3prime_v2_kept"][0])] + \
        [(f"tags_t1v2_{f.name}", f) for f in sorted(Path(TAG_LISTS["t1v2"][0]).glob("chr*.tags.tsv"))] + \
        [("chain", REPO / "data/broadaway_eqtl/hg19ToHg38.over.chain")]
    with open(a.out_dir / "input_manifest.tsv", "w") as fh:
        fh.write("role\tpath\tsha256\n")
        for role, p in manifest:
            fh.write(f"{role}\t{p}\t{sha256(p)}\n")
        for st in sorted(need_gwas):
            p = gwas_path(st, registry, af_sources)
            fh.write(f"gwas_{st}\t{p}\tsize={p.stat().st_size};mtime={int(p.stat().st_mtime)}\n")
    log("wrote", a.out_dir)
    print(json.dumps(summary, indent=2, default=str))


if __name__ == "__main__":
    main()
