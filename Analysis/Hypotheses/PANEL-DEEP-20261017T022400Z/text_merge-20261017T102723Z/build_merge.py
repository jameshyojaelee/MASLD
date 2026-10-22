"""APPLIED to the live files by team-lead on 2026-10-17 06:56; it no longer runs against live text.
Follow-up edits are in followup_check.py.

Merge the staged text fixes (edits.json, 64 edits) with the panel-deep legend
changes and the decision-dependent text outside the legends into one ordered
list of exact old/new edits.

Each new edit gives its `old` text with line breaks written as single spaces;
`locate` finds the one live substring that matches it with any run of
whitespace and returns that exact substring, so `old` in edits.json is always a
verbatim slice of the live file. Edits from the 64-edit list that a merged edit
replaces are dropped and their ids are recorded in `source_ids` of the merged
edit. No live file is written.
"""
import json
import os
import re
from collections import defaultdict

ROOT = "/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas"
OUT = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/text_merge-20261017T102723Z")
TF = os.path.join(ROOT, "Analysis/Hypotheses/PANEL-QUALITY-20261017T020832Z/text_fixes-20261017T025625Z")
CAND = "figures/candidates/panel-deep-fixes-20261017T035633Z"
os.chdir(ROOT)

LEG = "docs/manuscript/05_figure_legends.md"
RES = "docs/RESULTS.md"
STA = "docs/STATUS.md"
MET = "docs/manuscript/METHODS.md"
IDX = "figures/main/INDEX.md"
D = "docs/manuscript/draft/"

_cache = {}


def live(path):
    if path not in _cache:
        _cache[path] = open(path, encoding="utf-8").read()
    return _cache[path]


def locate(path, text):
    """Return the unique live substring of `path` equal to `text` up to whitespace."""
    parts = text.split()
    pat = r"\s+".join(re.escape(p) for p in parts)
    hits = [m.group(0) for m in re.finditer(pat, live(path))]
    if len(hits) != 1:
        raise SystemExit(f"locate: {len(hits)} matches in {path} for: {text[:90]!r}")
    return hits[0]


def lc(fig):
    return f"{CAND}/{fig}/legend_changes.md"


NEW = []


def add(id_, path, old, new, reason, decision, source_ids, producing_file):
    NEW.append(dict(id=id_, file=path, old=locate(path, old), new=new, reason=reason,
                    decision=decision, source_ids=source_ids, producing_file=producing_file))


# --------------------------------------------------------------------------
# Figure 1 legend
# --------------------------------------------------------------------------
add("M-1B-1", LEG,
    "**(B)** Interface between jointly testable liver-trait regulatory evidence and established disease-state "
    "transcription in the native 774-unit interface, restricted to the 3,769 genes whose multi-signal COLOC state "
    "is positive (SuSiE PP.H4 > 0.5) or tested at PP.H4 ≤ 0.5: 426 colocalized, 220 canonical disease-state and 26 shared.",
    "**(B)** Overlap between genes colocalized with liver eQTL and disease vs control DEGs. The comparison covers "
    "the 3,769 genes whose multi-signal COLOC state is positive (SuSiE PP.H4 > 0.5) or tested at PP.H4 ≤ 0.5. Of "
    "these, 426 are colocalized and 220 are canonical DEGs from the 774-cluster bulk fit (774 genetic-identity "
    "clusters: 138 control, 636 disease; not verified donors). The unit is the gene, and 26 genes are on both lists.",
    "1B legend Change 1, merged with L-01 (unit of the 774 bulk units).",
    "redraw", ["legend_changes:1B-1", "L-01"],
    f"{CAND}/fig1/source_tables/fig1b.tsv; Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1/REVIEW.md (138/636)")

add("M-1B-2", LEG,
    "Each bar is a fraction of its own gene list. Dashed lines mark the overlap expected within 10 AveExpr quantile "
    "bins (expression deciles), 21.21 genes (4.98% of 426 and 9.64% of 220); Mantel-Haenszel OR 1.295 (95% CI 0.835 "
    "to 2.007), exact two-sided stratified P=0.2755 (one prespecified test, no correction).",
    "Gray bars show the exact null distribution of the overlap within 10 AveExpr quantile bins, computed as the "
    "convolution of 10 hypergeometric distributions. The darker bars are the central 95%, 13 to 30 genes. The dashed "
    "line marks the expected overlap, 21.21 genes. The filled point is the observed overlap (26) and the open circle "
    "is the expected overlap without stratification (24.87). The lower row shows the Mantel-Haenszel OR, 1.295 (95% "
    "CI 0.835 to 2.007), on a log axis. Exact two-sided stratified P = 0.2755 (one prespecified test, no correction).",
    "1B legend Change 2.", "redraw", ["legend_changes:1B-2"], f"{CAND}/fig1/source_tables/fig1b.tsv")

add("M-1B-3", LEG,
    "Each of 1,000 prespecified lists of 426 non-colocalized genes, matched one-to-one to the colocalized genes on "
    "AveExpr, carried 13 DEGs against 21.3 expected within the same bins, so expression is labelled not fully "
    "controlled and no direction is read from this comparison. Without stratification 24.87 are expected (5.84% and 11.30%).",
    "The gray diamond is the overlap for 426 non-colocalized genes matched one-to-one to the colocalized genes on "
    "AveExpr: 13 DEGs (exact P = 0.052 under the same null). All 1,000 prespecified matched lists carried 13 DEGs. "
    "The lists are nearly identical: across 200 reruns of the matching, 411 of 426 genes appeared in every list. So "
    "13 is one control count, not a second expectation. Expression is labelled not fully controlled, and no direction "
    "is read from this comparison.",
    "1B legend Changes 3 and 4 (the unstratified expectation is now drawn and stated in Change 2).",
    "redraw", ["legend_changes:1B-3", "legend_changes:1B-4"],
    f"{CAND}/fig1/source_tables/fig1b.tsv; Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig1/nc1_list_diversity.json")

add("M-1B-5", LEG,
    "A further 10,928 genes with a bulk test and a tier-1/2 COLOC record, including 935 canonical DEGs, have only an "
    "ABF result and no multi-signal COLOC result; they are not counted as negatives.",
    "A further 10,928 genes with a bulk test and a tier-1/2 COLOC record, including 935 canonical DEGs, have only an "
    "ABF result and no multi-signal COLOC result. They are not counted as negatives. Another 699 canonical DEGs have "
    "no tier-1/2 COLOC record (503 lncRNA, 63 IG_V genes, 49 protein-coding). The comparison therefore covers 220 of "
    "the 1,854 canonical DEGs.",
    "1B legend Change 5.", "redraw", ["legend_changes:1B-5"],
    f"{CAND}/fig1/source_tables/fig1b_degs_outside_comparison.tsv")

add("M-1C", LEG,
    "**(C)** Phenotype support of the 436 primary tier-1/2 SuSiE genes in the adopted COLOC r2 release: 414 "
    "liver-enzyme-only, 12 direct-disease/liver-fat-only and 10 supported by both classes. Support uses any eligible "
    "SuSiE gene–trait pair, not only each gene's strongest driver.",
    "**(C)** Trait-class support of the 436 genes colocalized with liver eQTL (SuSiE PP.H4 > 0.5 in any of 35 main "
    "tier-1/2 GWAS) in the adopted COLOC r2 release. The unit is the gene. Each bar is one trait class: 424 genes are "
    "colocalized in ALT/AST/GGT GWAS (20 studies) and 22 in NAFLD/NASH/PDFF GWAS (15 studies). Ten genes are on both "
    "bars. Segments give each gene's SuSiE state in the other trait class: also supported (10 in each bar); tested at "
    "PP.H4 ≤ 0.5 (46 enzyme genes, 11 direct genes); or not SuSiE-tested in any study of that class (outlined; 368 "
    "enzyme genes, 1 direct gene). Not SuSiE-tested means the gene has no multi-signal COLOC result in that class. It "
    "is untested, not negative. Of the 414 genes supported only in enzyme GWAS, 46 were tested in a NAFLD/NASH/PDFF "
    "GWAS. Of the 12 supported only in direct GWAS, 11 were tested in an enzyme GWAS. The two classes had unequal "
    "testing: 3,888 genes were SuSiE-tested in at least one enzyme study and 392 in at least one direct study. Among "
    "the 362 genes tested in both, 56 were enzyme-supported and 21 direct-supported. Support uses any eligible SuSiE "
    "gene–trait pair, not only each gene's strongest driver. PNPLA3 and TM6SF2 have no SuSiE result. Their loci enter "
    "as SAMM50 (NAFLD, PP.H4 up to 0.62) and MAU2 (ALT, AST and NAFLD, PP.H4 1.00), whose top SNPs are at the PNPLA3 "
    "and TM6SF2 coding lead positions (22:44324727 and 19:19379549, hg19).",
    "1C legend Change 7: untested drawn as its own state.", "1", ["legend_changes:1C-7"],
    f"{CAND}/fig1/source_tables/fig1c.tsv; {CAND}/fig1/source_tables/fig1c_common_denominator.tsv")

add("M-1-6", LEG,
    "Percentages in (B) use each gene list as the denominator; technical units are never counted as donors.",
    "Counts in (B) and (C) are genes. Technical units are never counted as donors.",
    "1B legend Change 6 (the panel no longer draws per-list percentages).", "redraw", ["legend_changes:1B-6"],
    f"{CAND}/fig1/source_tables/fig1b.tsv")

# --------------------------------------------------------------------------
# Figure 2 legend
# --------------------------------------------------------------------------
add("M-2A", LEG,
    "Brackets bind the endpoint bands to their separate multi-signal COLOC (436 genes) and single-signal COLOC-only "
    "(560 genes) denominators.",
    "The three displays count different units and are not linked: GWAS studies by ancestry within each trait (35); "
    "fine-mapped loci per trait by the ancestry with the highest PIP (360), with the 54 loci whose top PIP is shared "
    "by more than one ancestry, because the cross-ancestry joint PIP is copied to every ancestry, drawn as their own "
    "segment; and the 996 colocalized genes by PP.H4 band, in the multi-signal (436) and single-signal-only (560) "
    "sets. Each gene is counted once, in the band of its largest PP.H4 (SuSiE for multi-signal genes, ABF for "
    "single-signal-only genes): above 0.9 (244 multi-signal and 67 single-signal-only genes), 0.8 to 0.9 inclusive "
    "(89 and 89) and above 0.5 but below 0.8 (103 and 404). The 560 single-signal-only genes split into 334 with no "
    "SuSiE result, 223 with SuSiE PP.H4 ≤ 0.5 and 3 that SuSiE could not test.",
    "2A legend change merged with L-02 (PP.H4 band rule and counts; same values in the redrawn 2A source table).",
    "redraw", ["legend_changes:2A", "L-02"],
    f"{CAND}/fig2/source_tables/2A.tsv; {CAND}/fig2/source_tables/prep_log.json")

add("M-2B", LEG,
    "Of 3,918 genes with a multi-signal posterior, 866 have a SuSiE PIP at that lead variant and 327 are shown (PP.H4 "
    "> 0.05 or PIP ≥ 0.99). The three MVP European strata were not fine-mapped, so none of the 1,348 genes whose "
    "largest posterior comes from them is shown. Eighteen genes have PIP ≥ 0.9 with PP.H4 < 0.5 and 18 have PIP ≥ 0.9 "
    "with PP.H4 > 0.5. Color gives the hg19 consequence of the plotted lead variant: noncoding (311), coding sequence "
    "or splice site of the plotted gene (5), or of another gene only (11); coding includes synonymous positions. "
    "Absence from this panel is not a negative colocalization call.",
    "Of 3,918 genes with a multi-signal posterior, 866 have a SuSiE PIP at that lead variant and all 866 are drawn: "
    "the 327 that pass the display filter (PP.H4 > 0.05 or PIP ≥ 0.99) in colour and the other 539 in light gray. The "
    "three MVP European strata were not fine-mapped, so the 1,348 genes whose largest posterior comes from them have "
    "no PIP and are not drawn; for example SAMM50, whose lead is the PNPLA3 variant 22:44324727 (PP.H4 0.619 in MVP "
    "NAFLD EUR). Among the 327 coloured genes, 18 have PIP ≥ 0.9 with PP.H4 < 0.5 and 18 have PIP ≥ 0.9 with PP.H4 > "
    "0.5; among all 866, 25 have PIP ≥ 0.9 with PP.H4 < 0.5. The dashed line marks PP.H4 = 0.5 and the dotted line "
    "PIP = 0.9. Colour and shape give the hg19 consequence of the plotted lead variant: noncoding (311, filled "
    "circle), coding sequence or splice site of the plotted gene (5, filled triangle), or of another gene only (11, "
    "open circle); coding includes synonymous positions. Labels are limited to a fixed five-gene roster (HSD17B13, "
    "TM6SF2, MBOAT7, GCKR and RORA); TM6SF2, MBOAT7 and GCKR have no lead SuSiE PIP, so only HSD17B13 and RORA are "
    "labeled. Absence from this panel is not a negative colocalization call.",
    "2B legend change merged with L-03 (label roster; the dashed-line sentence of L-03 is covered by the 2B change).",
    "redraw", ["legend_changes:2B", "L-03"],
    f"{CAND}/fig2/source_tables/2B.tsv; {CAND}/fig2/source_tables/prep_log.json; figures/main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc.pdf (labels)")

add("M-2C", LEG,
    "The 996-gene multi-signal-or-single-signal COLOC union partitioned by the strongest consequence annotation of "
    "one lead COLOC variant per gene: 35 coding, 70 UTR, 104 promoter, 450 intronic, and 337 intergenic.",
    "The 996-gene multi-signal-or-single-signal COLOC union partitioned by the strongest consequence annotation of "
    "one lead COLOC variant per gene (the variant with the highest PP.H4 from SuSiE or ABF; 933 distinct lead "
    "variants), with intronic and intergenic assigned by direct overlap with hg19 knownGene transcripts: 35 coding, "
    "70 UTR, 104 promoter, 524 intronic, and 263 intergenic.",
    "2C legend change.", "redraw", ["legend_changes:2C"], f"{CAND}/fig2/source_tables/2C.tsv")

add("M-2D", LEG,
    "Because 106 of 207 direct but 82 of 6,208 enzyme credible sets are European, this pooled value is not a "
    "trait-class contrast.",
    "Because 106 of 207 direct but 82 of 6,208 enzyme credible sets are European, this pooled value is not a "
    "trait-class contrast. Each bar is labelled with its credible sets and 1-Mb clusters (direct 207 sets, 19 "
    "clusters; enzyme 6,208 sets, 228 clusters), and the inset shows the direct-minus-enzyme difference in "
    "protein-altering plus canonical-splice probability, 19.5 percentage points (95% cluster bootstrap 6.9 to 34.0). "
    "The value depends on the weighting: averaging over credible sets instead of over clusters gives 15.8% versus 3.3%.",
    "2D legend change.", "redraw", ["legend_changes:2D"],
    f"{CAND}/fig2/source_tables/2D.tsv; Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig2/recompute/out_D/D_main_fractions.tsv")

add("M-2D-anc", LEG,
    "Thick bars are 95% and thin bars 98.333% intervals (Bonferroni over the three prespecified members EUR+AMR, EUR "
    "and AMR; 10,000 shared draws over 81 pooled clusters, seed 20260923); the pooled and all-non-EUR rows are "
    "descriptive. No member receives a trait-class call: the EUR enzyme arm was fine-mapped at 8 of 767 lead windows, "
    "the EUR member also failed the power check, and all three members failed their random-split calibration controls.",
    "Thick black bars are 95% and thin gray bars 98.3% cluster-bootstrap intervals (Bonferroni over the three "
    "prespecified members EUR+AMR, EUR and AMR; 10,000 shared draws over 81 pooled clusters, seed 20260923); the "
    "pooled and all-non-EUR rows are descriptive. The n column gives credible sets and, in parentheses, 1-Mb clusters "
    "per arm; clusters are re-derived within each stratum, and EUR+AMR is the sum of the EUR and AMR rows. Using "
    "pooled clusters instead would give EUR direct 31.8% rather than 29.9%. No member receives a trait-class call "
    "(\"no call\" on the panel): the EUR enzyme arm was fine-mapped at 8 of 767 lead windows, the EUR member also "
    "failed the power check (minimum detectable effect 0.304, larger than the pooled difference 0.195), and all three "
    "members failed their random-split calibration controls.",
    "2D ancestry-panel legend change.", "redraw", ["legend_changes:2D-ancestry"],
    f"{CAND}/fig2/source_tables/2D-ancestry.tsv")

add("M-2D-sub", LEG,
    "**(D, direct-subclass panel)** Three EUR PDFF GWAS (37 sets, 12 clusters) 40.14% and NAFLD plus NASH (170 sets, "
    "13 clusters) 12.44%, each against all enzyme traits (3.87%); 95% intervals from 10,000 physical-cluster draws, "
    "seed 20260923.",
    "**(D, direct-subclass panel)** Three EUR PDFF GWAS (37 sets, 12 clusters) 40.14% and NAFLD plus NASH (170 sets, "
    "13 clusters) 12.44%, each against all enzyme traits (3.87%): differences 36.3 (18.3 to 54.6) and 8.6 (−1.6 to "
    "21.0) percentage points, 95% intervals from 10,000 physical-cluster draws, seed 20260923. PDFF minus NAFLD plus "
    "NASH is 27.7 (7.2 to 47.7; subclasses resampled independently, 20,000 draws, seed 20261017). Because all PDFF "
    "sets are EUR, PDFF is also compared with EUR enzyme sets only (82 sets, 14 clusters, 6.92%): 33.2 (14.4 to 52.2; "
    "10,000 draws, seed 20261017). NAFLD plus NASH sets are 69 EUR, 95 AMR, 5 AFR and 1 EAS; enzyme sets are 82 EUR, "
    "715 AMR, 1,664 AFR, 3,694 EAS and 53 SAS.",
    "2D direct-subclass legend change.", "redraw", ["legend_changes:2D-subclass"],
    f"{CAND}/fig2/source_tables/2D-subclass.tsv; {CAND}/fig2/source_tables/prep_log.json")

add("M-2E", LEG,
    "**(E)** Global ancestry partition of the 436-gene multi-signal COLOC evidence set: 340 EUR-only, 64 shared with "
    "EUR, 23 non-European-only but sub-threshold at the lead variant, and 9 non-European-only genes passing "
    "lead-variant genome-wide significance (8 EAS and 1 AFR). The partition is per gene, not per locus: 19 of the 32 "
    "non-European-only genes have their lead variant within 1 Mb of the lead variant of a gene with European "
    "multi-signal PP.H4 > 0.5. The 9 genome-wide-significant genes fall in 8 loci (COQ6 and ENTPD5 share one); 5 of "
    "them, in 4 loci (COQ6/ENTPD5, ELAPOR2, LUZP1, RGS12), have no such European neighbor.",
    "**(E)** Ancestry partition of the 436 genes with multi-signal COLOC support (SuSiE PP.H4 > 0.5 in at least one "
    "study). Of the 340 genes supported only in EUR studies, 127 were tested in at least one non-EUR study (finite "
    "SuSiE PP.H4) without support, and 213 were not tested in any non-EUR study, so their EUR-only status is untested "
    "rather than unsupported. Of the 127 and the 213, 12 and 22 respectively have a non-EUR single-signal ABF PP.H4 > "
    "0.5 (black strip). 64 genes are supported in EUR and at least one non-EUR study. Of the 32 genes supported only "
    "in non-EUR studies, 9 have a lead p < 5e-8 in that non-EUR GWAS (8 EAS and 1 AFR) and 23 do not. The partition "
    "is per gene, not per locus: 19 of the 32 non-European-only genes have their lead variant within 1 Mb of the lead "
    "variant of a gene with European multi-signal PP.H4 > 0.5. Of the 9 with lead p < 5e-8, 4 (ACACB, GLTP, SUPT7L, "
    "ZFYVE21) have such a European neighbor, and 5, in 4 loci (COQ6/ENTPD5, ELAPOR2, LUZP1, RGS12), do not.",
    "2E legend change (decision 1); 'Of these, 12 and 22' reworded so each count names its group.",
    "1", ["legend_changes:2E"],
    f"{CAND}/fig2/source_tables/2E.tsv; {CAND}/fig2/source_tables/2E_gene_level.tsv")

add("M-2F", LEG,
    "Circle area denotes the displayed PP.H4.",
    "Circle area denotes the displayed PP.H4, and the size key uses the same scale as the data circles. A gene is "
    "shown if PP.H4 > 0.5, by SuSiE or ABF, in four of the five ancestry groups; by SuSiE alone only GGT1 and PANX1 "
    "reach four groups (asterisk). Each column is the maximum over a different number of studies: EUR 17, AFR 6, EAS "
    "6, AMR 3 and SAS 3.",
    "2F legend change.", "redraw", ["legend_changes:2F"], f"{CAND}/fig2/source_tables/2F.tsv")

add("M-2GH", LEG,
    "**(G-H)** Matched locus views showing GWAS association, fine-mapping, COLOC, expression, and local gene context. "
    "The canonical examples are RORA and FABP1. Both use EUR GWAS and do not provide evidence for "
    "population-dependent biology.",
    "**(G-H)** Locus views over 500-kb windows (hg19 positions in Mb) showing GWAS association coloured by r² to the "
    "lead in the EUR LD reference, SuSiE fine-mapping PIP where a fine-mapping record exists, liver eQTL association "
    "for the colocalized gene (Broadaway, N = 1,183) and local genes, whose body fill gives the bulk log2FC. **(G)** "
    "RORA in MVP NAFLD EUR (N = 432,709; 29,342 cases). The lead variant 15:60899031 has p = 7.6e-8, and no variant in "
    "the window reaches p < 5e-8. This GWAS has no stored fine-mapping record, so no PIP track is drawn. COLOC PP.H4 = "
    "0.98 (SuSiE); the top shared variant 15:60878030 has SNP.PP.H4 = 0.38. **(H)** FABP1 in UKBB ALT (N = 343,850). "
    "The PIP track is the stored COLOC-targeted SuSiE fine-mapping of this GWAS. The lead variant 2:88423510 "
    "(rs1545223, p = 1.8e-12) has PIP 0.57 and the COLOC top variant 2:88424066 has PIP 0.43; both are drawn as hollow "
    "diamonds over their points. COLOC PP.H4 = 1.00 (SuSiE) and SNP.PP.H4 = 0.56. FABP1 is not differentially "
    "expressed in bulk liver (log2FC −0.04, padj 0.51), and KRCC1 is lower in MASLD (log2FC −0.10, padj 3.1e-5). Both "
    "use EUR GWAS and do not provide evidence for population-dependent biology.",
    "2G/2H legend change merged with L-04. The redraw removes the RORA display-only SuSiE fit that L-04 described, so "
    "only the FABP1 stored-PIP statement of L-04 is kept (variant count not repeated: the redrawn window holds 1,289 "
    "PIP rows, not the 1,135 of the earlier window).",
    "redraw", ["legend_changes:2G-2H", "L-04"],
    f"{CAND}/fig2/source_tables/2G.tsv; {CAND}/fig2/source_tables/2H.tsv")

add("M-2I", LEG,
    "**(I)** Shared legend for the two locus panels.",
    "**(I)** Shared key for G and H. Gene-body fill is the bulk log2FC for MASLD versus control from the adopted fit "
    "over 774 genetic-identity clusters, capped at ±0.3; gray genes have padj ≥ 0.05, and white genes with a black "
    "outline are outside the bulk DE gene universe.",
    "2I legend change.", "redraw", ["legend_changes:2I"], f"{CAND}/fig2/source_tables/2I.tsv")

# --------------------------------------------------------------------------
# Figure 3 legend
# --------------------------------------------------------------------------
add("M-3AB", LEG,
    "**(A–B)** Five-cohort coverage (326, 207, 140, 71 and 30 units) and documented fibrosis-by-NAS support. "
    "GSE213621 records fibrosis only as source bins (F0–F1, F2, F3–F4) and enters no exact-stage contrast. The joint "
    "histology grid contains 364 units with both an exact fibrosis stage and NAS.",
    "**(A–B)** Five-cohort coverage in genetic-identity clusters (326, 207, 140, 71 and 30 units; disease/control "
    "258/68, 197/10, 109/31, 65/6 and 7/23) and the records each cohort supplies. Sex is inferred from expression in "
    "GSE213621 and GSE135251, which record no sex. GSE162694 records NAS for 117 and fibrosis for 109 of its 140 "
    "units; the 31 without fibrosis are its controls. GSE213621 records fibrosis only as source bins (F0–F1, F2, "
    "F3–F4) and enters no exact-stage contrast. The grid in B counts the 364 units with both an exact fibrosis stage "
    "and NAS, controls included; the top row gives the 387 units with an exact stage (103/102/88/69/25) used in E and "
    "F. The 23 staged GSE162694 units without NAS are in the top row only.",
    "3A-B legend change.", "redraw", ["legend_changes:3A-B"],
    f"{CAND}/fig3/source_tables/fig3a.tsv; {CAND}/fig3/source_tables/fig3b.tsv")

add("M-3C", LEG,
    "Ribbon widths represent 20,054 cohort memberships.",
    "Bars give genes, not cohort memberships. The top bars give calls per cohort split by direction (GSE135251 "
    "7,381; GSE126848 6,147; GSE162694 3,199; GSE213621 2,045; GSE130970 1,282). The middle bars give genes per "
    "overlap class, darker where at least one call comes from GSE126848. The bottom bars give the mixed-sign share of "
    "multi-cohort genes: 1,865 of 4,696 (39.7%) with all five cohorts and 346 of 2,565 (13.5%) when GSE126848 is left "
    "out of the overlap count.",
    "3C legend change (decision 8).", "8", ["legend_changes:3C"], f"{CAND}/fig3/source_tables/fig3c.tsv")

add("M-3D", LEG,
    "Both density and scatter use the same 392 units. Controls are open gray circles; advanced disease follows the "
    "pale-red-to-magenta fibrosis ramp.",
    "Both density and scatter use the same 392 units (113 strict controls, 279 advanced). Controls are open gray "
    "circles. Units with an exact fibrosis stage are filled circles, one colour per stage (F0 7, F1 42, F2 59, F3 69, "
    "F4 25); the 108 F0–F2 units qualify by NAS ≥5 alone. The 77 GSE213621 units recorded only as the F3–F4 source "
    "bin are open triangles. Density curves split advanced units by entry rule: fibrosis F3 or higher (171) and NAS "
    "≥5 only (108). PC1 separates advanced units from controls within each cohort (AUC 0.83 to 0.90). Among advanced "
    "units, PC1 rank correlates with fibrosis stage (Spearman 0.48, n = 279; 0.48 among the 202 with an exact stage).",
    "3D legend change.", "redraw", ["legend_changes:3D"], f"{CAND}/fig3/source_tables/fig3d.tsv")

add("M-3E", LEG,
    "Five available k6 axes show continuous F0-centered activity; k6-P6 is unavailable.",
    "The lighter inner bars count the genes that also have |log2FC| >0.5: 182 / 798 / 1,756 / 4,019 for F1–F4 (up "
    "130 / 615 / 1,265 / 2,030; down 52 / 183 / 491 / 1,989). Stage n in units is printed under each stage. Five "
    "available k6 axes show F0-centered activity with 95% intervals (estimate ± 1.96 SE; points offset slightly along "
    "x so intervals do not overlap); k6-P6 is unavailable. Intervals include 0 for P3 at F1 and F2, P4 at F1 and F2, "
    "and P5 at F1, F3 and F4.",
    "3E legend change.", "redraw", ["legend_changes:3E"], f"{CAND}/fig3/source_tables/fig3e.tsv")

add("M-3F", LEG,
    "Colors show higher-minus-lower stage log2FC, clipped at ±0.688 for display; black dots mark BH significance.",
    "Colours show higher-minus-lower stage log2FC, clipped at ±0.688 for display (the 98th percentile of |log2FC|); a "
    "triangle marks each of the five clipped cells. Filled dots mark BH q <0.0125 (the 5% budget split over four "
    "contrasts; 16 marks in 15 genes) and open dots mark 0.0125 ≤ q <0.05. The top block holds the 25 genes with a "
    "mark at F1→F2 or later; the bottom block holds the 29 selected only through F0→F1. The trait column names the "
    "GWAS of the strongest tier-1/2 driver: all 54 are liver-enzyme traits (ALT, AST or GGT). The PP.H4 lane runs "
    "from 0.5 to 1.",
    "3F legend change.", "redraw", ["legend_changes:3F"], f"{CAND}/fig3/source_tables/fig3f.tsv")

add("M-3G", LEG,
    "**(G)** Raw disease-state RNA log2FC versus unsigned SuSiE COLOC PP.H4 for 3,776 finite native RNA–genetics "
    "pairs. Because PP.H4 is unsigned, its association is tested against |log2FC| (Spearman rho=0.001, nominal "
    "P=0.93). Points are plotted without a fitted line. Colors distinguish Neither (3,156, including seven "
    "SuSiE-untestable genes; 3,149 among the 3,769 genes of Figure 1B), Transcriptomics sig. (194), Genetics sig. "
    "(400) and Both (26). Transcriptomics sig. uses the",
    "**(G)** Disease-state RNA |log2FC| versus unsigned SuSiE COLOC PP.H4 for the 3,769 genes with a SuSiE test. "
    "Dashed lines mark the two class rules (|log2FC| = 0.5 and PP.H4 = 0.5). The lower strip expands PP.H4 below "
    "0.05; there, the 2,939 Neither genes are drawn as counts in 0.1 × 0.005 bins (darker = more genes) and the other "
    "genes as points. Colours distinguish Neither (3,149), Transcriptomics only (194), Genetics only (400) and Both "
    "(26). The seven SuSiE-untestable genes are open circles at their recorded values and are not counted as Neither "
    "or tested. Both holds 26 genes against 24.9 expected under independence (Fisher OR 1.06, P = 0.83). Spearman rho "
    "between |log2FC| and PP.H4 is 0.001 (P = 0.95). Including the seven untestable genes (3,776 genes) gives OR "
    "1.06, P = 0.74 and rho 0.0014, P = 0.93. Points are plotted without a fitted line. Transcriptomics only uses the",
    "3G legend change merged with T-20 (rho 0.0014 for the 3,776-gene universe, the precision used in RESULTS and "
    "fig3.md; 0.001016 for the 3,769-gene universe prints as 0.001).",
    "redraw", ["legend_changes:3G", "T-20"], f"{CAND}/fig3/source_tables/fig3g.tsv (universe rows)")

add("M-3L", LEG,
    "The lncRNA inset uses the 774-unit RNA fit: 520 one-to-one covariate-matched lncRNA/protein-coding pairs, "
    "selected at an outcome-blind 0.15-SD caliper, give a direction odds ratio of 0.132 (95% CI 0.090–0.193).",
    "The lncRNA inset uses the 774-unit RNA fit. Left: odds of a call for lncRNA versus protein-coding genes over "
    "21,501 genes, unadjusted (open, OR 1.12, 95% CI 1.00–1.25) and adjusted for expression and detection (filled, "
    "OR 0.43, 0.38–0.49), on a log scale. Top right: direction of calls, 83.0% up among 1,157 protein-coding calls "
    "and 47.9% up among 549 lncRNA calls. Bottom right: 520 one-to-one covariate-matched lncRNA/protein-coding pairs, "
    "selected at an outcome-blind 0.15-SD caliper, give a conditional direction odds ratio of 0.132 (30 discordant "
    "pairs with the lncRNA member up against 227 with the protein-coding member up; Wald 95% CI 0.090–0.193).",
    "lncRNA inset legend change merged with L-05 (discordant-pair counts and interval type).",
    "redraw", ["legend_changes:lncRNA", "L-05"],
    f"{CAND}/fig3/source_tables/fig3_lncrna.tsv; figures/main/fig3_bulk_transcriptomics/source_tables/fig3_lncrna_role/summary.tsv (30, 227)")

# --------------------------------------------------------------------------
# Figure 4 legend
# --------------------------------------------------------------------------
add("M-4ABC", LEG,
    "**(A)** Seven single-cell datasets, 104 descriptive biological donors, 1.23 million cells, one annotated atlas, "
    "102 analyzed donors, and 117 prespecified programs across five tested lineages. **(B)** Integrated UMAP with a "
    "cell-type color legend and donor-averaged disease density; sequencing runs and sort fractions are not treated as "
    "donors. **(C)** The three positive and three negative Hotspot effects with the smallest 117-program permutation "
    "BH q in each of the four non-T lineages. Heatmap columns are median dataset-centered donor scores in Healthy, "
    "Steatosis and Steatohepatitis, using only GSE202379 and GSE244832, the two datasets with more than one stage. "
    "Points are ordinal stage β from generalized least squares with a separate residual variance per dataset, on the "
    "60 donors whose datasets could be ambient-corrected; fitted donors are hepatocytes 55, cholangiocytes 38, "
    "fibroblasts 43 and macrophages 25. Black bars are nominal 95% confidence intervals. Gray bars are β ± the "
    "program's within-dataset permutation critical |t| at 0.05/117 times its standard error. No program passes the "
    "117-program permutation BH/stability rule (smallest q=0.323), and no gray bar excludes zero. T cell programs are "
    "omitted from this compact main panel; the complete 117-program registry view, including T cell programs, is "
    "shown in Figure S4.",
    "The single-cell atlas integrates 275 libraries and 1,232,318 cells from 7 datasets. 273 libraries (1,232,285 "
    "cells) from 102 donors are analyzed. 117 fixed Hotspot programs are defined in 5 lineages. **(B)** Integrated "
    "UMAP of 50,000 cells drawn uniformly from the 1,232,285 analyzed cells of 102 donors, coloured by cell type. The "
    "key gives the number of analyzed cells per type. Dendritic cells include plasmacytoid and conventional dendritic "
    "cells. Sequencing runs and sort fractions are not treated as donors. **(C)** Ordinal stage β for all 117 Hotspot "
    "programs, ranked within each of the five lineages, from generalized least squares with a separate residual "
    "variance per dataset, on the 60 donors whose datasets could be ambient-corrected (GSE202379 and GSE244832). "
    "Fitted donors are hepatocytes 55, fibroblasts 43, cholangiocytes 38, macrophages 25 and T cells 53. Thick gray "
    "bars are β ± the program's within-dataset permutation critical |t| at 0.05/117 times its standard error. Thin "
    "black bars are nominal 95% confidence intervals. No program passes the 117-program permutation BH/stability rule "
    "(smallest q = 0.323), and no gray bar excludes zero. Median dataset-centered donor scores by stage are shown in "
    "Figure S4.",
    "4A (caption only), 4B and 4C legend changes; 4A text becomes the unlettered opening of the Figure 4 legend.",
    "redraw", ["legend_changes:4A", "legend_changes:4B", "legend_changes:4C"],
    f"{CAND}/fig4/source_tables/4B.tsv; {CAND}/fig4/source_tables/4C.tsv; Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute/r4a_atlas_library_counts.tsv")

add("M-4D", LEG,
    "In 1,000 donor bootstrap draws (seed 20260923), 76 matches and 34 mismatches keep the same call in at least 95% "
    "of draws and 7 are unstable;",
    "In 1,000 donor bootstrap draws (seed 20260923), 76 matches and 34 mismatches keep the same call in at least 95% "
    "of draws and 7 are unstable. Cells that hold unstable calls print their number (hepatocyte/hepatocyte 1, "
    "hepatocyte/cholangiocyte 1, fibroblast/fibroblast 1, fibroblast/T cell 2, cholangiocyte/cholangiocyte 2). "
    "Matching cells are shaded;",
    "4D legend change.", "redraw", ["legend_changes:4D"], f"{CAND}/fig4/source_tables/4D.tsv")

add("M-4E", LEG,
    "Each point is one of 112 eligible programs. The x axis shows the effect in score SD per stage-code unit or NAS "
    "point from a three-cohort random-effects meta-analysis (REML, Knapp-Hartung); the y axis shows -log10 BH q over "
    "224 tests (112 programs × two axes), and the dotted line marks q = 0.05.",
    "The panel shows the fibrosis axis. NAS effects are in Figure S4. Each column is one of 112 eligible programs, "
    "ranked by its fibrosis effect (cohort-standardized program score SD per recorded stage-code unit). Within each "
    "cohort the program score is fitted by OLS on sex, fibrosis stage and NAS with HC3 standard errors, and the three "
    "cohort estimates are combined by random-effects meta-analysis (REML, Knapp-Hartung). Filled points are coloured "
    "by discovery lineage. Gray bars are Knapp-Hartung 95% confidence intervals (t with 2 df). Open marks are the "
    "per-cohort estimates (circle GSE130970, n = 71; triangle GSE135251, n = 207; square GSE162694, n = 86). The five "
    "largest fibrosis effects are labelled: Activated stellate (PDGFRA) 0.417, Versican / GAG synthesis 0.387, "
    "Ductular injury (BICC1) 0.376, Stromal ECM (IGFBP7) 0.320 and ECM organization 0.289. Each has the same sign in "
    "all three cohorts. BH correction is over 224 tests (112 programs × fibrosis and NAS). With three cohorts, the "
    "result depends on the Knapp-Hartung adjustment: z-based REML inference gives 37 fibrosis and 14 NAS associations "
    "at q < 0.05.",
    "4E legend change merged with L-06 (within-cohort OLS covariates and HC3 standard errors).",
    "redraw", ["legend_changes:4E", "L-06"],
    f"{CAND}/fig4/source_tables/4E.tsv; Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute/job_9842_22276639.log")

add("M-4F", LEG,
    "scores across nine prespecified fixed-projection windows in all five Resource cohorts. GSE126848, GSE130970, and "
    "GSE135251 overlap continuum construction and are shown as dashed gray descriptive trajectories; GSE162694 and "
    "GSE213621 are the competitor-independent evaluation cohorts.",
    "scores across nine prespecified fixed-projection windows in all five Resource cohorts (774 genetic-identity "
    "clusters: 30 GSE126848, 71 GSE130970, 207 GSE135251, 140 GSE162694 and 326 GSE213621). GSE126848, GSE130970, and "
    "GSE135251 overlap continuum construction and are shown as dashed gray descriptive trajectories. GSE162694 and "
    "GSE213621 are the competitor-independent evaluation cohorts (solid lines). Lines join window means. Shaded bands "
    "are ± 1 standard error. The key gives the genetic-identity clusters per window (GSE126848 6, GSE130970 14-15, "
    "GSE135251 41-42, GSE162694 28, GSE213621 65-66). Windows are 0.2 wide on the within-cohort fixed-projection "
    "percentile, at steps of 0.1, so adjacent windows share units. A unit on a window boundary is assigned to the "
    "upper window, and the last window includes its upper end. The Stromal ECM score excludes 9 continuum-signature "
    "genes that carry 19.5% of its weight (C7, CLSTN2, COL4A1, COL4A2, DCN, EPHA3, IGFBP7, LAMA2, PLCXD3).",
    "4F legend change merged with L-07. L-07 gave GSE130970 13 to 15 per window; the redraw assigns boundary units to "
    "the upper window and its source table gives 14 to 15, so the redraw value is used.",
    "redraw", ["legend_changes:4F", "L-07"],
    f"{CAND}/fig4/source_tables/4F.tsv (n_units by dataset); hotspot_testability.tsv (9 excluded genes)")

add("M-4G", LEG,
    "**(G)** Stage-adjusted continuum coefficients overlaid on the fixed, signature-excluded pathway/Hotspot "
    "membership geometry and its 43 systems.",
    "**(G)** Stage-adjusted continuum coefficient for each of the 5,435 testable nodes (gray points; one point per "
    "pathway or Hotspot program) in the 43 fixed, signature-excluded systems. There is one row per system, in "
    "registry order by size, with the node count in brackets. Open circles are system medians with 95% "
    "donor-bootstrap intervals (disease units resampled within cohort, B = 1,000, seed 20261017). The median is above "
    "zero in 38 of 43 systems. Four nodes are untestable and are counted at the right of their rows (S02 1, S04 2, "
    "S05 1). The membership-geometry view is in Figure S4.",
    "4G legend change.", "redraw", ["legend_changes:4G"], f"{CAND}/fig4/source_tables/4G.tsv")

add("M-4H", LEG,
    "**(H)** Fibrosis-stage coefficient (score SD per cohort stage-code SD) versus stage-adjusted continuum "
    "coefficient (score SD per continuum SD) for all 43 systems. Spearman rho is 0.63 and 32/43 effects share "
    "direction; both summaries are descriptive.",
    "**(H)** Fibrosis-stage coefficient (system score SD per cohort stage-code SD) versus stage-adjusted continuum "
    "coefficient (system score SD per continuum SD) for all 43 systems. Both are estimated donor-first: each "
    "participant's collection-balanced system score is fitted within each cohort, and the two cohorts are combined by "
    "fixed-effect meta-analysis. These are the estimates in the system inference tables. Bars are 95% donor-bootstrap "
    "intervals (B = 1,000, seed 20261017). Colour marks systems whose two intervals both exclude zero (19 both "
    "positive, 2 both negative, 1 opposite signs); 21 systems are not resolved (gray). Spearman rho is 0.67 (bootstrap "
    "95% interval 0.38 to 0.84) and 34/43 point estimates share direction; both summaries are descriptive. The "
    "feature-first summary (median node coefficient within system and collection, averaged over collections) gives "
    "rho 0.63 and 32/43.",
    "4H legend change (decision 7: donor-first estimand).", "7", ["legend_changes:4H"],
    f"{CAND}/fig4/source_tables/4H.tsv")

# --------------------------------------------------------------------------
# Figure 5 legend
# --------------------------------------------------------------------------
add("M-5A", LEG,
    "**(A)** Input separation between gene-level context, the prespecified program display set from the complete "
    "117-program registry, and the fixed descriptive protein display.",
    "Each panel has its own input set. Panel B uses the 10,044-gene display universe locked on 2026-07-10 from the "
    "July COLOC release, of which 3,054 genes are measured as proteins. Panel D uses the 436 genes in the 767 primary "
    "gene-study pairs of the corrected COLOC r2 release; 48 of these genes are outside the panel B universe and 189 "
    "are among the 3,054 genes plotted in panel B. Panel E shows programs from the 117-program registry that have a "
    "lineage-open colocalized member gene or a robust protein-histology association. Panel F uses two programs, "
    "Stromal ECM and Ductular injury, that were prespecified separately and are not in the panel E display set. Panel "
    "C uses 25 proteins fixed before the bulk refit, with the nine mRNA calls from the 774-unit refit.",
    "5A legend change (caption only); panel letters are kept, relettering is an assembly step.",
    "redraw", ["legend_changes:5A"],
    f"{CAND}/fig5/source_tables/5B.tsv; {CAND}/fig5/source_tables/5C_right_column.tsv; Analysis/Spatial/results/universe_validation/prioritized_universe_FINAL.txt")

add("M-5B-1", LEG,
    "The scatter retains the nonsignificant background; the adjacent bar includes protein-significant genes only and "
    "shows concordant and discordant counts and within-significant percentages.",
    "The scatter shows all 3,054 measured prioritized genes; blue points are protein-significant (limma BH-adjusted P "
    "< 0.05, BH over all 6,071 quantified proteins) with the same sign as mRNA, magenta points are "
    "protein-significant with the opposite sign, and gray points are not protein-significant. The adjacent bars split "
    "the 211 protein-significant genes by the RNA call: 170 of 177 genes with RNA BH padj < 0.05 have the same sign "
    "in mRNA and protein (96.0%), against 18 of 34 genes with RNA padj 0.05 or higher (52.9%). The dashed line marks 50%.",
    "5B legend change merged with L-09 (protein significance rule and BH family).",
    "redraw", ["legend_changes:5B-1", "L-09"],
    f"{CAND}/fig5/source_tables/5B_bars.tsv; {CAND}/fig5/source_tables/5B_summary.tsv; Analysis/Multimodal_Program_Projection/results/proteomics/protein_de_adjusted.tsv")

add("M-5B-2", LEG,
    "The current 774-unit PDF contains 3,054 measured genes and prints 211 protein-significant genes: 188 concordant "
    "(89.1%) and 23 discordant (10.9%) with mRNA; 2,843 are not protein-significant.",
    "The panel prints n = 3,054 genes. 211 genes are protein-significant: 188 have the same sign as mRNA and 23 the "
    "opposite sign; 2,843 are not protein-significant.",
    "5B legend change (the pooled 89.1% is no longer drawn).", "redraw", ["legend_changes:5B-2"],
    f"{CAND}/fig5/source_tables/5B_summary.tsv")

add("M-5C", LEG,
    "The histology column shows partial Spearman rho over all participants with the protein quantified (58 for 22 of "
    "25 proteins; 56 for ABHD5, 43 for P4HA2 and 40 for HLA-DQA1), including controls with fibrosis 0-1 and NAS 0-1, "
    "so it combines disease status with severity.",
    "The histology columns show partial Spearman rho between each protein and each histology score within the 46 "
    "MASLD participants (rank residuals adjusted for acquisition batch, age, BMI and sex), with participant bootstrap "
    "95% intervals (2,000 resamples). Filled points have BH q < 0.05 over the 125 tests (25 proteins x 5 scores); 60 "
    "of 125 pass. All 46 MASLD participants have each protein quantified except ABHD5 (45), P4HA2 (37) and HLA-DQA1 "
    "(32); these n are printed beside the rows. Controls are excluded so that the association reflects severity "
    "within disease, not disease status. Over all 58 participants, 111 of 125 tests pass and 9 rho values change sign "
    "relative to the MASLD-only estimate (source table column rho_all58). Process NES and q come from preranked fgsea "
    "on the adjusted protein moderated t statistic (gene sets of 8 to 600 proteins), with BH over all 1,145 tested "
    "sets (1,091 Reactome, 50 Hallmark and 4 curated). Mito/OXPHOS has the largest process change (adjusted NES -2.6, "
    "q = 2e-14) and is listed below the fixed rows because none of the 25 fixed proteins belongs to it. The heatmap "
    "strip marks the 46 MASLD and 12 control participants; columns are ordered by fibrosis stage and then NAS.",
    "5C legend change (decision 6) merged with L-10 (fgsea family; the L-10 BH-over-125 sentence is already in the "
    "decision-6 text).",
    "6", ["legend_changes:5C", "L-10"],
    f"{CAND}/fig5/source_tables/5C.tsv; {CAND}/fig5/source_tables/5C_nes.tsv; figures/candidates/bulk-context-774-20261014T140904Z_22236311/proteomics/protein_pathway_enrichment.tsv")

add("M-5D-1", LEG,
    "Bars separate direct MASLD/liver-fat genes (n=22; 28 gene-study pairs) from liver-enzyme genes (n=424; 739 pairs)",
    "Facets separate direct MASLD/liver-fat genes (n = 22 genes; 28 gene-study pairs) from liver-enzyme genes (n = "
    "424 genes; 739 pairs)",
    "5D legend change (facets, not bars).", "redraw", ["legend_changes:5D-2"], f"{CAND}/fig5/source_tables/5D.tsv")

add("M-5D-2", LEG,
    "The display retains states with at least 50% of the posterior in peaks of both sources, at least 50% in one "
    "source, 5% to under 50% in the better source, and under 5% in both sources (including no overlap); the 0.05 "
    "floor keeps trace overlap out of the partial state. Segments under 6.5% of genes are not labeled; all counts are "
    "in the source table.",
    "For each lineage and trait class, points show the percent of genes whose best gene-study pair places at least "
    "50% of the primary-signal posterior in peaks of at least one source (circle) or of both sources (triangle), with "
    "gene-bootstrap 95% intervals (10,000 resamples within trait class). Open diamonds show the expectation if the "
    "posterior were spread uniformly over each pair's tested variants; under that reference no gene reaches 50% in "
    "any lineage. Genes with 5% to under 50% in peaks are counted as under 50%; their number is in the source table, "
    "because under the uniform reference the median expected mass in peaks is already 1% to 6% by lineage, close to "
    "the 5% boundary.",
    "5D legend change.", "redraw", ["legend_changes:5D-1"],
    f"{CAND}/fig5/source_tables/5D.tsv; Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig5/recompute/5d_sensitivity_uniform_reference.tsv")

add("M-5E-1", LEG,
    "registry of 117 fixed programs; rows are not selected by disease association. Rows show the 45 programs with at "
    "least one member gene in an adopted COLOC r2 primary pair (the 767 pairs of panel D) that places at least half "
    "of its colocalization posterior in open chromatin of the program's own lineage in at least one snATAC source (28 "
    "programs), or with a robust protein-histology association (26 programs); 9 programs meet both.",
    "registry of 117 programs. Rows show 45 programs in two blocks. The upper block holds the 31 programs with at "
    "least one member gene, counted once per primary-assembly gene, in a COLOC r2 primary pair (the 767 pairs of "
    "panel D) that places at least half of its colocalization posterior in open chromatin of the program's own "
    "lineage in at least one snATAC source; every protein interval is shown for these. The lower block holds the 14 "
    "programs displayed only because their protein-histology association is robust; these rows are selected on the "
    "statistic plotted in the right column. 22 programs are robust and 8 meet both criteria.",
    "5E legend change, touch-up 2026-10-17 version (decision 9).", "9", ["legend_changes:5E-1"],
    f"{CAND}/fig5/source_tables/5E.tsv (display, block, n_open_genes, robust)")

add("M-5E-2", LEG,
    "and light shading a gene that is lineage-open only through liver-enzyme pairs.",
    "and light shading a gene that is lineage-open only through liver-enzyme pairs. Intervals that exclude 0 but fail "
    "the robustness rule are drawn as gray outline boxes (5 programs); intervals that include 0 are light gray.",
    "5E legend change (interval styles).", "9", ["legend_changes:5E-2"],
    f"{CAND}/fig5/source_tables/5E.tsv (interval_style)")

add("M-5E-3", LEG,
    "Freedman-Lane permutation P (9,999) is BH-corrected over the 94 programs with at least eight measured member "
    "proteins and 20% retained weight; \"not scored\" marks the two displayed programs below these thresholds "
    "(Elastic-fiber ECM, 6 measured proteins; Memory/MAIT, 3). Robust means q<0.05 with the equal-weight score, the "
    "score without its top gene and both acquisition batches agreeing in sign.",
    "Program protein scores use the same members as the gene column, counted once per primary-assembly gene. "
    "Freedman-Lane permutation P (9,999) is BH-corrected over the 95 programs with at least eight measured member "
    "proteins and 20% retained weight; \"not scored\" marks the two displayed programs below these thresholds "
    "(Elastic-fiber ECM, 7 measured proteins; Memory/MAIT, 3). Robust means q<0.05 with the equal-weight score, the "
    "score without its top gene and both acquisition batches agreeing in sign. One robust program, cholangiocyte "
    "Mitochondrial / metabolic (ETFB), sits at the BH boundary (P = 0.0121, q = 0.04998).",
    "5E legend change (decision 9), plus the Chol 5 boundary sentence that the decision-9 recommendation asks for.",
    "9", ["legend_changes:5E-3", "REPORT.md decision 9"],
    f"{CAND}/fig5/source_tables/5E.tsv (testable, n_measured, pvalue, padj); {CAND}/fig5/work/protein_histology_117_primary.tsv (n_permutations)")

add("M-5F-1", LEG,
    "labels are BH-adjusted empirical q-values from 9,999 matched-gene draws, BH over the two programs within each "
    "source. The null resamples genes with donors and arrays held fixed, so the q-values are not donor-level tests; "
    "donor-level values are in Results. Stromal ECM exceeds its null in both sources (supported in Guilliams et al.; "
    "source-dependent in Vu et al. because its arrays lack a public donor key), whereas Ductular injury is "
    "indeterminate in both.",
    "Small gray points show each donor (Guilliams et al., 4 donors over 5 sections) or array (Vu et al., 10 arrays); "
    "the larger point is the pooled excess, and the gray band is the matched-null 95% range. Labels give the "
    "BH-adjusted empirical q-value from 9,999 matched-gene draws (BH over the two programs within each source) and, "
    "beside it, the donor-level result. The gene null holds donors and arrays fixed, so it is not a donor-level test. "
    "Stromal ECM exceeds its gene null in both sources (q = .001 and q < .001). At donor level, 3 of 4 Guilliams et "
    "al. donors have positive excess (exact one-sided sign-flip P = 0.19; the smallest possible P with 4 donors is "
    "0.0625). In Vu et al. 10 of 10 arrays are positive, but the arrays have no public donor key, so the result is "
    "source-dependent and drawn as an open point. Ductular injury is indeterminate in both sources (q = .079 and "
    ".316; 3 of 4 donors and 7 of 10 arrays positive). Row labels give genes measured over program genes (Stromal ECM "
    "26/49 in Guilliams et al. and 25/49 in Vu et al.; Ductular injury 9/22 in both).",
    "5F legend change (decision 5), with the row-label sentence.", "5", ["legend_changes:5F-1", "legend_changes:5F-3"],
    f"{CAND}/fig5/source_tables/5F.tsv (pooled_point rows)")

add("M-5F-2", LEG,
    "of the five sections, it has the largest Stromal ECM excess (Results).",
    "of the five sections, it has the largest Stromal ECM excess (Results). The map title prints this section's "
    "excess (0.102, 1,164 spots) and the range over the five sections (−0.022 to 0.102).",
    "5F legend change (map title).", "5", ["legend_changes:5F-2"], f"{CAND}/fig5/source_tables/5F.tsv (section_excess rows)")

# --------------------------------------------------------------------------
# Figure 6 legend
# --------------------------------------------------------------------------
add("M-6A", LEG,
    "**(A)** Catalog structure separates inherited-trait and disease-state RNA evidence before stable-identity, "
    "testability and provenance qualification, molecular-object interpretation and a proposed discriminating "
    "experiment. Chromatin and spatial context are attached at program level for two programs (Stromal ECM and "
    "Ductular injury), not as gene-level calls; liver protein (PXD051911) was untestable for both. No record is a "
    "regulatory DNA element object.",
    "The Catalog keeps the eQTL-colocalization call and the bulk RNA call separate, qualifies each by identity, "
    "testability and provenance, assigns one molecular object per gene identity and proposes one discriminating "
    "experiment. The 28,233 molecular-object records are 14,176 protein-coding genes, 13,576 untestable nominations, "
    "423 lncRNA transcripts and 58 ambiguous lncRNA loci; 0 are regulatory DNA elements. Liver protein was not "
    "assessed at gene level, and chromatin and spatial context were not joined at gene level, for all 28,233 genes. "
    "Chromatin and spatial context are attached at program level for two programs (Stromal ECM and Ductular injury); "
    "liver protein (PXD051911) was untestable for both.",
    "6A legend change (caption only).", "redraw", ["legend_changes:6A"],
    "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig6/recompute/fig6a_object_type_by_route.tsv; catalog_r9d_gene_table.parquet")

add("M-6B-1", LEG,
    "**(B)** Genetics-call by RNA-call matrix over 28,233 selected r8 gene identities.",
    "**(B)** Liver-trait eQTL-colocalization call by bulk RNA call over 28,233 gene identities with r8 classes, "
    "unchanged in Catalog r9d (adopted 2026-10-16). The genetic call is supported when SuSiE-COLOC PP.H4 > 0.5 in any "
    "of 35 Tier-1/2 GWAS (15 direct MASLD, 20 liver enzyme). The RNA call is supported when BH q < 0.05 and |log2FC| "
    "> 0.5 in one limma-voom fit of 774 genetic-identity clusters (138 control, 636 disease) with dataset and "
    "inferred sex as covariates; BH is over the 24,115 genes in the fit. Tested-in-both cells are blue and untested "
    "cells gray. \"Coding PIP ≥ 0.5\" counts genes with a fine-mapped coding variant at SuSiE PIP ≥ 0.5 (Catalog r9d "
    "field; 31 genes: 20 unresolved, 7 genetically anchored, 3 established-state associated, 1 untested). This field "
    "does not enter the genetic call, so PNPLA3, TM6SF2 and MTARC1 stay unresolved. The untestable-by-untestable cell "
    "is 0 by construction: every Catalog identity enters from at least one input.",
    "6B legend change (decision 4) merged with I-04 (r9d release label).", "4", ["legend_changes:6B-1", "I-04"],
    f"{CAND}/fig6/source_tables/6B.tsv")

add("M-6B-2", LEG,
    "Genetics `untestable` means no Tier-1/2 COLOC record (9,312 genes) or SuSiE-untestable for insufficient shared "
    "posterior in one Tier-1/2 study (7); RNA `untestable` means not detected (4,232).",
    "Genetics `untestable` means no Tier-1/2 COLOC record (9,312 genes) or SuSiE-untestable for insufficient shared "
    "posterior in one Tier-1/2 study (7). These 7 genes have finite SuSiE PP.H4 in 1 to 8 other Tier-1/2 studies "
    "(maximum 0.164); they are kept untestable under the r8 rule, and under the rule used for other tested genes they "
    "would be indeterminate (unresolved 13,134). RNA `untestable` means not detected (4,232).",
    "6B legend change (decision 4: keep the 7 genes untestable and footnote them).", "4", ["legend_changes:6B-2"],
    "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig6/recompute/diagnostics.json; sensitivity_matrix_and_classes.tsv")

add("M-6C", LEG,
    "*HKDC1* is supported in both domains: SuSiE PP.H4 0.988 (UKBB AST) and log2FC +1.07 (BH q 8.0e-7). *MTARC1* is "
    "unresolved: no Tier-1/2 SuSiE result (best ABF PP.H4 0.046, FinnGen NASH) and log2FC -0.15 (q 1.2e-4). *THRB* is "
    "genetically anchored: SuSiE PP.H4 > 0.999 (UKBB GGT) and log2FC -0.20 (q 4.2e-10).",
    "*HKDC1* is supported in both domains: SuSiE PP.H4 0.988 (UK Biobank AST; > 0.5 in 3 of 5 Tier-1/2 studies with a "
    "SuSiE result, all liver enzyme) and log2FC +1.07 (95% CI 0.68 to 1.47; BH q 8.0e-7). *MTARC1* is unresolved: no "
    "Tier-1/2 SuSiE result (best ABF PP.H4 0.046, FinnGen NASH), so the eQTL-colocalization call is indeterminate, "
    "although the gene carries a fine-mapped missense variant (PIP 0.76; PIP ≥ 0.5 in 3 studies, top study a PDFF "
    "GWAS) that does not enter the call; log2FC -0.15 (-0.22 to -0.08; q 1.2e-4). *THRB* is genetically anchored: "
    "SuSiE PP.H4 > 0.999 (UK Biobank GGT; > 0.5 in 1 of 4 studies) and log2FC -0.20 (-0.26 to -0.14; q 4.2e-10). ABF "
    "PP.H4 is shown for reference and does not make a call. Intervals are log2FC ± 1.96 SE. \"Both supported\" "
    "requires a supported call in each domain and does not compare direction. Routes use the gene-level strings of "
    "(D); the molecular-object route for protein-coding genes is gene perturbation with an orthogonal rescue.",
    "6C legend change (decision 4).", "4", ["legend_changes:6C"], f"{CAND}/fig6/source_tables/6C.tsv")

add("M-6D", LEG,
    "The printed experiments are proposals, and the counts are not a target-priority ranking.",
    "The printed experiments are proposals, and the counts are not a target-priority ranking. Bar lengths follow from "
    "the class rules in (B). Untested segments are gray and tested classes blue. Hatched: the 708 untested genes with "
    "one supported call (698 RNA in the no-genetic-call segment; 10 genetics in the no-RNA-call segment, narrower than "
    "a line). The unresolved bar includes 20 genes with a fine-mapped coding variant at PIP ≥ 0.5, among them PNPLA3, "
    "TM6SF2 and MTARC1; their printed route was assigned from the eQTL-colocalization class and does not use the "
    "coding variant.",
    "6D legend change (decision 4).", "4", ["legend_changes:6D"], f"{CAND}/fig6/source_tables/6D.tsv")

add("M-6E", LEG,
    "companion. **(E)** Actual local access to worked gene entries, molecular-object uncertainty, source references, "
    "gene calls and molecular-object downloads; no combined score or rank determines the interface.",
    "companion.",
    "6E legend change: the (E) panel moves to Data availability (edit D-6E in METHODS.md).", "redraw",
    ["legend_changes:6E"], "catalog_r9d_gene_table.parquet (shape 28,233 x 126)")

add("M-6-I06", LEG,
    "All five individual panels are installed on October 14;",
    "Panels B, C and D are the October 17 redraws, and the content of the October 14 panels A and E is now in this "
    "legend and in Data availability;",
    "Supersedes I-06: B, C and D were replaced by the panel-deep redraws on 2026-10-17; A and E became text.",
    "redraw", ["I-06"],
    "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/install-20261017T102551Z/install_log.json")

# --------------------------------------------------------------------------
# Figure 7 legend
# --------------------------------------------------------------------------
add("M-7A", LEG,
    "**(A)** Evaluation framework. Every model is scored on a task whose own baseline is defined in advance, and "
    "every score is reported beside a per-cohort record of what that model saw during its own training. Where the "
    "exposure record is unresolved, the score is shown as uninterpretable rather than as evidence of skill. The "
    "framework census is the historical September 1 snapshot: 106 model directories and 104 audit files. Its bar "
    "contains 89 checkpoint findings (58 unexposed, 20 unknown, 9 other and 2 seen) plus 37 audit-file placeholders "
    "without a checkpoint finding, totaling 126 entries; it is not a count of independent models. These checkpoint "
    "findings are separate from the updated study-specific development exposure in panels C and F.",
    "Every model is scored on a task whose own baseline is defined in advance, and every score is reported beside a "
    "per-study record of what that model saw during its own training. A score whose exposure record is unresolved is "
    "treated as uninterpretable, not as evidence of skill. At the September 1 census, 106 model directories held 104 "
    "exposure audit files. These record 89 checkpoint findings: 58 unexposed (65%), 20 unknown, 9 in other recorded "
    "states and 2 seen by the encoder in pretraining. A further 37 audit entries carry no checkpoint finding and are "
    "not counted as findings. Eight audit files changed after the census, and the counts above are not updated for "
    "them. These checkpoint findings are separate from the study-specific exposure drawn in panels C and F.",
    "7A legend change (caption only).", "redraw", ["legend_changes:7A"],
    f"{CAND}/fig7/source_tables/fig7A_caption_trace.tsv; panel_7a_framework.tsv")

add("M-7B", LEG,
    "**(B)** Gain over each task's own baseline. On a synthetic reporter assay (GSE281364, 1,033 elements, HepG2 "
    "control and palmitate-oleate contexts only), one of seven pretrained sequence-model rows (six models), HyenaDNA, "
    "exceeds a ridge on allele identity alone after Bonferroni correction across the seven rows (gain +0.082, 95% "
    "block-bootstrap interval +0.027 to +0.139; corrected lower bound +0.009 from the same bootstrap draws). "
    "Sequence-based reporter baselines (deltaSVM, gkm-SVM, a sequence CNN or transformer, and the source study's "
    "MPRALegNet) were not run, so this is a gain over a sequence-free control. On the endogenous variant-to-gene |z| "
    "task, adding either pretrained embedding to distance lowers the score (HyenaDNA −0.022, Nucleotide Transformer "
    "−0.029; intervals exclude zero). No pretrained sequence model was scored on tissue-level severity: the severity "
    "rows are a composition-plus-expression model and a routed multi-modal system.",
    "**(B)** Gain over each task's own baseline. Each task has its own x axis, and the axis names the baseline and its "
    "score. On a synthetic reporter assay (GSE281364, 1,033 elements in 239 blocks, HepG2 control and "
    "palmitate-oleate contexts only), one of seven pretrained sequence-model rows (six models), HyenaDNA, exceeds a "
    "ridge on allele identity alone (Spearman 0.134) after Bonferroni correction across the seven rows. Its gain is "
    "+0.082, with a Bonferroni-7 interval of +0.009 to +0.158 (99.29% paired block bootstrap, 10,000 draws). The "
    "intervals of the other six rows include zero; Nucleotide Transformer is +0.054 (−0.018 to +0.130). Sequence-based "
    "reporter baselines (deltaSVM, gkm-SVM, a sequence CNN or transformer, and the source study's MPRALegNet) were not "
    "run, so this is a gain over a sequence-free control. On the endogenous variant-to-gene |z| task (held-out LD fold "
    "4: 111,426 rows in 55 LD blocks), adding either pretrained embedding to distance lowers the held-out Spearman "
    "below a distance ridge (Spearman 0.234). HyenaDNA changes it by −0.022 and Nucleotide Transformer by −0.029 "
    "(paired 95% LD-block bootstrap intervals, 2,000 resamples, seed 20260831: −0.033 to −0.010 and −0.042 to −0.016). On "
    "the same rows the best embedding head predicts eQTL sign at AUC 0.502 (0.490 to 0.513). A ridge on allele "
    "features alone scores 0.503, and within-block permutation chance is 0.500. For graded fibrosis in 521 bulk RNA "
    "samples from five cohorts, samples are the unit because no donor key exists. Composition plus "
    "lineage-residualized expression scores 0.128 macro Spearman below plain expression with matched degrees of "
    "freedom (95% sample bootstrap −0.176 to −0.082; 1 of 5 cohorts positive). It still adds +0.173 (+0.078 to "
    "+0.268) over hepatocyte fraction alone. No pretrained sequence model was scored on tissue-level severity. Bulk "
    "RNA predicts H3K27ac in 99 GSE267145 participants with a log2 CPM squared-error skill of 0.122 (95% participant "
    "bootstrap 0.083 to 0.159). The same rank-4 map, transferred to GSE296875 ATAC in 39 donors, has a multinomial "
    "deviance skill of 0.130 (donor jackknife 0.098 to 0.163). The routed multi-modal system is in panels D and E, "
    "because its output is identical to the reference model.",
    "7B legend change (decision 2) merged with L-12 (held-out fold, n, metric, paired LD-block bootstrap). Intervals "
    "are the producer values (team-lead 2026-10-17): eQTL from dnalm_heads 03_intervals.json, seed 20260831; MPRA "
    "from round-4 mpra_bonferroni.tsv (0.082008, Bonferroni 0.008947 to 0.158465). The panel is being redrawn to these.",
    "2", ["legend_changes:7B", "L-12"],
    "Analysis/MASLD_Model_Benchmark/executions/eqtl_sign_reliability-20260831T195702Z/logs/03_strata_intervals.json "
    "(strata/all: best_auc, best_auc_ci, chance_measured.median; headline_provenance.allele_feature_ridge_sign_auc_recomputed); "
    "Analysis/MASLD_Model_Benchmark/executions/dnalm_heads-20260831T183647Z/logs/03_intervals.json "
    "(paired_deltas/abs_z m1_dist_emb MINUS b1_dist); "
    "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/repair/F7a/exec/mpra/mpra_bonferroni.tsv "
    "(F7_displayed_rows)")

add("M-7C", LEG,
    "**(C)** Score beside training exposure. Exposure is resolved per study from each model's own audit and "
    "development crosswalk, never transcribed. TranscriptFormer TF-Sapiens has five source-qualified clean studies "
    "and two training-seen studies; its pooled broad cell-class macro F1 of 0.9389 includes both and is marked "
    "accordingly, not interpreted as clean-only skill. Geneformer V2 has zero qualified clean studies among the seven "
    "displayed cohorts and therefore no clean-only delta.",
    "**(C)** Score beside training exposure. Exposure is resolved per study from each model's own audit and "
    "development crosswalk, never transcribed. Each encoder's donor-class-balanced macro F1 for broad cell classes is "
    "drawn twice. The hollow mark uses all seven held-out studies (102 donors, 50,000 cells). The filled mark uses "
    "only that encoder's clean studies. Intervals are 95% study-stratified donor-cluster bootstraps (donors resampled "
    "within study; 2,000 draws; seed 9873). TranscriptFormer TF-Sapiens, with five clean studies (80 donors), scores "
    "0.939 (0.927 to 0.949) with the two-layer head on all seven studies and 0.963 (0.957 to 0.969) on its clean "
    "studies. Geneformer V2 has no clean study among the seven, so its scores (0.930 linear, 0.944 two-layer) are not "
    "evidence of skill. HVG+PCA with the same head is the gray reference (0.880 linear, 0.926 two-layer). Scores on "
    "different study sets are not compared with each other; panel F compares each encoder with HVG+PCA on the same "
    "studies. Five task-native classifiers on HVG+PCA features are not drawn; they score 0.702 (nearest centroid) to "
    "0.908 (k-nearest neighbors).",
    "7C legend change; replaces T-27 (n now stated in the new text) and L-13 (Geneformer V2 is drawn in the redrawn 7C).",
    "redraw", ["legend_changes:7C", "T-27", "L-13"], f"{CAND}/fig7/source_tables/fig7C.tsv")

add("M-7D", LEG,
    "The severity-transfer block shows the released instrument's Spearman, with its 95% interval, beside each "
    "cohort's own tie-structure ceiling in the inferential units of its three transfer cohorts: 109 GSE268273 "
    "participant records, 263 GSE213621 genetic identities (one record each; 11 identities whose two libraries carry "
    "different stage labels removed) and GSE276114 participants drawn as two etiology rows (81 MASLD and 82 chronic "
    "viral hepatitis).",
    "The signed-eQTL block shows the best embedding head's held-out sign AUC in each |z| stratum, with the number of "
    "held-out rows and LD blocks. Its 95% LD-block bootstrap interval is not corrected across the five strata. Each "
    "AUC is drawn beside within-block permutation chance (gray square) and the sequence-attributable bound (open "
    "diamond), which assumes at most one causal variant per gene and LD block. The severity-transfer block shows the "
    "released expression instrument (released weights, no refit), with its Spearman correlation with stage and a 95% "
    "bootstrap interval (units resampled within stage, 10,000 draws, seed 20260923). Each value is drawn beside that "
    "cohort's own stage tie ceiling (open diamond). Units are those of the three transfer cohorts: 109 GSE268273 "
    "participant records, 263 GSE213621 genetic identities (one record each; 11 identities whose two libraries carry "
    "different stage labels removed) and GSE276114 participants in two etiology rows (81 MASLD and 82 chronic viral "
    "hepatitis). GSE268273 also chose the instrument's expression transform. In the LD-block join, a row's effect "
    "correlates with its gene's mean effect at r = 0.867 when the row's own LD block is left in (644,556 "
    "training-fold rows in 208 LD blocks, 7,576 genes). It is 0.847 on the 15,192 rows (198 genes) that can be held "
    "out, and 0.021 on those same rows with the own block held out. The fitted routing rule returns output identical "
    "to the reference model (gain 0). An oracle router could add 0.0071 macro Spearman, below the pre-declared floor "
    "of 0.093 (521 bulk RNA samples, 5 cohorts; no donor key exists, so samples are the unit).",
    "7D legend change; replaces L-14 (interval method) and L-15 (n of each block; its 208 LD blocks and sample unit "
    "are kept). Severity interval is the F7b repair producer (seed 20260923; GSE268273 0.5685 to 0.7483), per "
    "team-lead 2026-10-17; the panel is being redrawn to it.",
    "redraw", ["legend_changes:7D", "L-14", "L-15"],
    "figures/candidates/review-r3-fixes-20261015T112148Z/F7b/repair/fig7_transfer_primary_units.tsv (ci_lo, ci_hi, ci_B, ci_seed)")

add("M-7E", LEG,
    "Signed eQTL direction keeps its prespecified verdict, tested negative: in every |z| stratum the headroom between "
    "the label-reliability ceiling and measured chance exceeds twice the test's resolution. If each gene and LD block "
    "holds at most one causal variant, a sequence model could gain at most 0.006 AUC above chance on all rows, below "
    "the 0.012 resolution, so the test cannot separate a model failure from a target that is not sequence-determined "
    "at these variants.",
    "Every row is placed on one scale, (estimate − null)/(threshold − null), so 0 is the null and 1 is the threshold "
    "its verdict uses. That threshold is the minimum detectable effect for the severity arms and the zone task, and "
    "the S-LDSC usability bar of h2 z = 7 for the GWAS rows. For signed eQTL it is twice the test's resolution, and "
    "for routing the pre-declared floor. Verdicts are those of the producing analyses. The severity rows draw all "
    "seven arms of the cross-modal stratification analysis, each labelled precommitted or post hoc as its producing "
    "file records it. Activity ≥ 4 on all participants is recorded there as a post hoc sensitivity, not a "
    "pre-registered arm. Intervals are 95% stratified bootstraps of the arm difference (10,000 draws). The producer "
    "records them for precision only, and verdicts use the label permutation null. As precommitted, an underpowered "
    "arm is labelled untestable, never negative. Signed eQTL direction is drawn as untestable. If each gene and LD "
    "block holds at most one causal variant, a sequence model could gain at most 0.006 AUC above chance on all rows. "
    "That bound (open diamond) is 0.24 of the threshold, so the test cannot separate a model failure from a target "
    "that is not sequence-determined at these variants. Routing is tested negative because even an oracle router "
    "reaches only 0.077 of the floor. Program heritability for ALT and AST is a count and is not drawn on this scale: "
    "0 of 117 programs have conditional tau significant at BH 0.05 within each trait.",
    "7E legend change (decision 3: all seven arms labelled; signed eQTL moved to untestable). Settles needs_user item 1. "
    "Signed-eQTL bound from the producer: (0.505453 - 0.499789) / (2 x 0.011627) = 0.24; the earlier recompute gave 0.26.",
    "3", ["legend_changes:7E", "needs_user:1"],
    f"{CAND}/fig7/source_tables/fig7E.tsv; Analysis/MASLD_Model_Benchmark/executions/eqtl_sign_reliability-20260831T195702Z/logs/03_strata_intervals.json "
    "(strata/all: sequence_attributable_bound, chance_measured.median, resolution_half_width)")

add("M-7F", LEG,
    "Per-study deltas are not drawn in this panel; the gains drawn are all-study values. The source correction "
    "changes only TranscriptFormer exposure states to five clean and two training-seen studies, preserving all "
    "displayed all-study gains and confidence intervals; these are not clean-only estimates.",
    "The grid shows, for each encoder and held-out study, whether the study was clean, seen in the encoder's "
    "pretraining, or unresolved. Column labels give the donors in each study, and Liver_Atlas is GSE192742. Gains over "
    "HVG+PCA with the same head are drawn on all seven studies (hollow) and on each encoder's clean studies (filled). "
    "Both come from the same donor-cluster bootstrap draws (95% intervals, uncorrected across rows, 2,000 draws, seed "
    "9873). On their clean studies, the four encoders that have any clean study gain +0.046 to +0.054 over HVG+PCA "
    "with the linear head and +0.010 to +0.015 with the two-layer head. The lower block compares each encoder with "
    "the scVI encoder trained here, on the encoder's clean studies. All four gain with the linear head (+0.012 to "
    "+0.018). With the two-layer head, SCimilarity gains +0.010 (+0.003 to +0.017), TranscriptFormer is level at "
    "−0.001 (−0.004 to +0.002), and UCE-4L and UCE-33L score lower (−0.005, −0.009 to −0.000; −0.007, −0.011 to "
    "−0.002).",
    "7F legend change; replaces L-16 (the drawn gains are no longer all-study only).", "redraw",
    ["legend_changes:7F", "L-16"], f"{CAND}/fig7/source_tables/fig7F.tsv")

# --------------------------------------------------------------------------
# Outside the legends: decisions 1, 5, 7, 9
# --------------------------------------------------------------------------
add("D1-fig1", D + "fig1.md",
    "414 enzyme-only, 12 direct-disease/liver-fat-only and 10 supported by both, using any eligible SuSiE pair rather "
    "than the strongest trait driver.",
    "414 enzyme-only, 12 direct-disease/liver-fat-only and 10 supported by both, using any eligible SuSiE pair rather "
    "than the strongest trait driver. Of the 414, 46 were SuSiE-tested in a NAFLD/NASH/PDFF GWAS and 368 in none; of "
    "the 12, 11 were tested in an enzyme GWAS. Among the 362 genes SuSiE-tested in both trait classes, 56 are "
    "enzyme-supported and 21 direct-supported.",
    "Decision 1: text quoting 414 enzyme-only gains the tested counts.", "1", ["REPORT.md decision 1"],
    f"{CAND}/fig1/source_tables/fig1c.tsv; {CAND}/fig1/source_tables/fig1c_common_denominator.tsv")

add("D1-fig2", D + "fig2.md",
    "(12 direct only and 10 with both classes; Figure 1C)",
    "(12 direct only, 11 of them SuSiE-tested in an enzyme GWAS, and 10 with both classes; among the 362 genes "
    "SuSiE-tested in both trait classes, 56 are enzyme-supported and 21 direct-supported; Figure 1C)",
    "Decision 1: direct-only count gains its tested count and the common denominator.", "1", ["REPORT.md decision 1"],
    f"{CAND}/fig1/source_tables/fig1c.tsv; {CAND}/fig1/source_tables/fig1c_common_denominator.tsv")

add("D1-RES-2E", RES,
    "Panel E ancestry partition: 340 EUR-only, 64 shared,",
    "Panel E ancestry partition: 340 EUR-only (127 SuSiE-tested in at least one non-EUR study, 213 not tested in "
    "any), 64 shared,",
    "Decision 1: text quoting 340 EUR-only gains 'of which N were tested'.", "1", ["REPORT.md decision 1"],
    f"{CAND}/fig2/source_tables/2E.tsv")

add("R-RES-2C", RES,
    "Panel C consequence counts: coding 35, UTR 70, promoter 104, intron 450, intergenic 337 (996);",
    "Panel C consequence counts: coding 35, UTR 70, promoter 104, intron 450, intergenic 337 (996; the 2026-10-17 "
    "redraw assigns intronic and intergenic by direct overlap with hg19 knownGene transcripts and prints intron 524, "
    "intergenic 263);",
    "Installed 2C now prints 524/263; the round-3 values stay as the record of that candidate.", "redraw",
    ["legend_changes:2C"], f"{CAND}/fig2/source_tables/2C.tsv")

add("R-STA-2CE", STA,
    "Installed Figure 2C prints 14/21/70/104/450/337 (coding split into own gene 14 and other gene 21); E prints "
    "340/64/23/9.",
    "Installed Figure 2C prints 14/21/70/104/524/263 since 2026-10-17 (coding split into own gene 14 and other gene "
    "21; intronic and intergenic by direct transcript overlap, 450/337 before); E prints 340/64/23/9, with the 340 "
    "split into 127 SuSiE-tested in a non-EUR study and 213 not tested.",
    "2C redraw value change and decision 1 (2E untested state).", "1", ["REPORT.md decision 1", "legend_changes:2C"],
    f"{CAND}/fig2/source_tables/2C.tsv; {CAND}/fig2/source_tables/2E.tsv")

add("D5-RES", RES,
    "Stromal ECM exceeded its matched-gene null in both sources and Ductular injury did not;",
    "Stromal ECM exceeded its gene-level matched null in both sources and Ductular injury did not; neither source "
    "gives a donor-level test (4 donors in GSE192741, arrays without a donor key in Vu);",
    "Decision 5: the both-sources sentence is a gene-level result, not a donor-level one.", "5",
    ["REPORT.md decision 5"], f"{CAND}/fig5/source_tables/5F.tsv (n_units_positive, signflip_p)")

add("D7-fig4-a", D + "fig4.md",
    "Among all 43 systems, stage and continuum effects have Spearman rho=0.63 and 32 share direction.",
    "Among all 43 systems, the donor-first stage and continuum estimates of the system inference tables have Spearman "
    "rho=0.67 (bootstrap 95% interval 0.38 to 0.84) and 34 share direction (feature-first summaries: 0.63 and 32).",
    "Decision 7: 4H plots the donor-first estimand.", "7", ["REPORT.md decision 7"], f"{CAND}/fig4/source_tables/4H.tsv")

add("D7-fig4-b", D + "fig4.md",
    "of the recorded stage code with the stage/sex-adjusted continuum coefficient in score SD per continuum SD, both "
    "estimated in the 367 disease units.",
    "of the recorded stage code with the stage/sex-adjusted continuum coefficient in score SD per continuum SD, both "
    "estimated donor-first in the 367 disease units, with 95% donor-bootstrap intervals on both axes; 22 systems have "
    "both intervals excluding zero.",
    "Decision 7: 4H plots donor-first estimates with intervals (19 + 2 + 1 = 22 resolved systems).", "7",
    ["REPORT.md decision 7"], f"{CAND}/fig4/source_tables/4H.tsv (state)")

add("D7-IDX", IDX,
    "Figure 4H prints rho=0.63 and 32/43 same direction, with no outline.",
    "Figure 4H plots the donor-first system estimates, rho=0.67 and 34/43 same direction (feature-first 0.63 and "
    "32/43), with no outline.",
    "Decision 7.", "7", ["REPORT.md decision 7"], f"{CAND}/fig4/source_tables/4H.tsv")

add("D7-RES-320", RES,
    "The current main panel uses the 367 disease units (rho 0.63, 32/43 same direction) and has no outline;",
    "The current main panel uses the 367 disease units and plots the donor-first system estimates (rho 0.670, 34/43 "
    "same direction; feature-first 0.631 and 32/43) and has no outline;",
    "Decision 7.", "7", ["REPORT.md decision 7"], f"{CAND}/fig4/source_tables/4H.tsv")

add("D7-RES-1713", RES,
    "0.631 and 32/43 on the 367 disease units; see Main Figure 4H below);",
    "0.631 and 32/43 on the 367 disease units; Main Figure 4H plots the donor-first estimates, 0.670 and 34/43);",
    "Decision 7: the equal-collection values are feature-first; 4H now draws donor-first.", "7",
    ["REPORT.md decision 7"], f"{CAND}/fig4/source_tables/4H.tsv")

# T-29 rewritten for decision 7: the 0.631 / 32/43 values are the feature-first 4H of 2026-10-16.
add("D7-RES-T29", RES,
    "(printed as 0.73) and 37/43 in the same direction; with the 68 GSE213621 controls removed from both axes (367 "
    "disease units) the values are 0.631 and 32/43 (candidate "
    "`figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/`)",
    "(printed as 0.73 on the earlier panel) and 37/43 in the same direction; with the 68 GSE213621 controls removed "
    "from both axes (367 disease units) the values are 0.631 and 32/43 (feature-first 4H installed 2026-10-16; "
    "`figures/candidates/review-r3-fixes-20261015T112148Z/F4efgh/fig4h/`). Since 2026-10-17 Figure 4H plots the "
    "donor-first estimates on the same 367 units: rho 0.670 (bootstrap 95% interval 0.382 to 0.842) and 34/43 "
    "(`figures/main/fig4_singlecell_programs/panels/data/panel_deep_20261017/4H.tsv`)",
    "T-29 merged with decision 7.", "7", ["T-29", "REPORT.md decision 7"], f"{CAND}/fig4/source_tables/4H.tsv")

add("D7-STA", STA,
    "4H has no outline, rho 0.631 and 32/43 same direction (0.7265 and 37/43 with the controls).",
    "4H has no outline and since 2026-10-17 plots the donor-first system estimates, rho 0.670 and 34/43 same "
    "direction (feature-first 0.631 and 32/43; 0.7265 and 37/43 with the controls).",
    "Decision 7.", "7", ["REPORT.md decision 7"], f"{CAND}/fig4/source_tables/4H.tsv")

add("D9-fig5-a", D + "fig5.md",
    "Twenty-eight programs contain at least one gene from an adopted colocalized pair",
    "Thirty-one programs contain at least one gene, counted once per primary-assembly gene, from an adopted "
    "colocalized pair",
    "Decision 9 consumer fig5.md:57 (28 to 31 from the 5 rejoined genes).", "9", ["REPORT.md decision 9"],
    f"{CAND}/fig5/source_tables/5E.tsv (n_open_genes)")

add("D9-fig5-b", D + "fig5.md",
    "in the other 25, open genes are open only through liver-enzyme pairs.",
    "in the other 28, open genes are open only through liver-enzyme pairs.",
    "Decision 9 consumer fig5.md:63-64.", "9", ["REPORT.md decision 9"],
    f"{CAND}/fig5/source_tables/5E.tsv (n_open_direct_genes)")

add("D9-fig5-c", D + "fig5.md",
    "26 of 94 scorable programs are robustly associated with a histologic-burden axis across participants (BH over "
    "94): 16 negatively,",
    "22 of 95 scorable programs are robustly associated with a histologic-burden axis across participants (BH over "
    "95; members counted once per primary-assembly gene): 14 negatively,",
    "Decision 9 consumer fig5.md:66-67.", "9", ["REPORT.md decision 9"],
    f"{CAND}/fig5/source_tables/5E.tsv (testable, robust, rho)")

add("D9-fig5-d", D + "fig5.md",
    "and 10 positively, including fibroblast ECM proteoglycan",
    "and 8 positively, including fibroblast ECM proteoglycan",
    "Decision 9 consumer fig5.md:69; the named examples all stay robust.", "9", ["REPORT.md decision 9"],
    f"{CAND}/fig5/source_tables/5E.tsv (robust, rho)")

add("D9-MET-a", MET,
    "programs is scored per participant as the L1-weighted mean of its measured member proteins; a program is "
    "scorable with at least eight measured members and at least 20% of its L1 weight retained (94 of 117).",
    "programs is scored per participant as the L1-weighted mean of its measured member proteins, with members "
    "counted once per primary-assembly gene, the rule used for the gene column of Figure 5E; a program is scorable "
    "with at least eight measured members and at least 20% of its L1 weight retained (95 of 117).",
    "Decision 9 consumer METHODS:1604-1609 (member rule; 94 to 95).", "9", ["REPORT.md decision 9"],
    f"{CAND}/fig5/work/protein_histology_117_primary.tsv (testable)")

add("D9-MET-b", MET,
    "and Benjamini-Hochberg correction over the 94 scorable programs.",
    "and Benjamini-Hochberg correction over the 95 scorable programs.",
    "Decision 9 consumer METHODS:1616.", "9", ["REPORT.md decision 9"],
    f"{CAND}/fig5/work/protein_histology_117_primary.tsv (correction_family)")

add("D9-MET-c", MET,
    "plus the program's position, from 1 to 94, among the scorable programs in registry order.",
    "plus the program's position, from 1 to 94, among the 94 programs scorable under the earlier symbol-level "
    "membership, in registry order. Portal fibroblast (NR1H4), scorable only under the primary-assembly membership, "
    "uses 1000 plus its registry row, 1032.",
    "Decision 9 consumer METHODS:1619-1620 (seed rule).", "9", ["REPORT.md decision 9"],
    f"{CAND}/fig5/work/protein_histology_117_primary.tsv (perm_seed, boot_seed); {CAND}/fig5/code/fig5e_protein_histology_117_primary_assembly.R:353")

add("D9-STA", STA,
    "45 rows with a lineage-open COLOC r2 gene (28 programs, 3 with direct-trait RORA) or a robust protein-histology "
    "association (26 of 94 scorable, BH over 94), and no spatial column.",
    "45 rows with a lineage-open COLOC r2 gene (31 programs, 3 with direct-trait RORA) or a robust protein-histology "
    "association (22 of 95 scorable, BH over 95; since 2026-10-17 members are counted once per primary-assembly gene "
    "in both columns), and no spatial column.",
    "Decision 9 consumer STATUS:860.", "9", ["REPORT.md decision 9"], f"{CAND}/fig5/source_tables/5E.tsv")

add("R-STA-5AD", STA,
    "5D counts genes with a 5% posterior floor: 22 direct and 424 enzyme genes, from 767 r2 gene-study pairs (28 "
    "direct, 739 enzyme); a locus-matched background companion uses the same genes. 5A retains assay-native context.",
    "5D gives, per lineage and trait class, the percent of genes with at least 50% of the posterior in peaks of one "
    "or both sources: 22 direct and 424 enzyme genes, from 767 r2 gene-study pairs (28 direct, 739 enzyme); a "
    "locus-matched background companion uses the same genes. 5A is legend text since 2026-10-17.",
    "5A and 5D redraws: the 5% floor states and 5A artwork are no longer drawn.", "redraw",
    ["legend_changes:5A", "legend_changes:5D-1"], f"{CAND}/fig5/source_tables/5D.tsv")

add("D-6E", MET,
    "a data dictionary, and the source-dependency graph used by the Gene Catalog.",
    "a data dictionary, and the source-dependency graph used by the Gene Catalog. The Catalog gene table (28,233 "
    "rows, 126 columns, none of which is a combined score or rank), the 28,233 molecular-object records and the "
    "column dictionary are deposited at [accession or path to be supplied]. No external portal is published.",
    "6E legend change: panel E content moves to Data availability.", "redraw", ["legend_changes:6E"],
    "Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/catalog_r9d-20261016T105550Z/build_out/catalog_r9d_gene_table.parquet (28,233 x 126)")

# --------------------------------------------------------------------------
# Outside the legends: draft text that describes the redrawn panels
# --------------------------------------------------------------------------
add("R-fig3-3G", D + "fig3.md",
    "across 3,776 finite native RNA–SuSiE pairs, Spearman rho between PP.H4 and |log2FC| is 0.0014 (nominal P=0.93), "
    "shown without a fitted line. The 426 genes with PP.H4 >0.5 split 218 up and 208 down. Neither, Transcriptomics "
    "sig., Genetics sig. and Both contain 3,156 / 194 / 400 / 26 genes; Neither includes seven SuSiE-untestable "
    "genes, so among the 3,769 genes evaluated by multi-signal COLOC (Figure 1B) Neither is 3,149.",
    "across the 3,769 genes evaluated by multi-signal COLOC (Figure 1B), Spearman rho between PP.H4 and |log2FC| is "
    "0.001 (nominal P=0.95), shown without a fitted line; with the seven SuSiE-untestable genes (3,776 genes) it is "
    "0.0014 (P=0.93). The 426 genes with PP.H4 >0.5 split 218 up and 208 down. Neither, Transcriptomics only, "
    "Genetics only and Both contain 3,149 / 194 / 400 / 26 genes; Both holds 26 against 24.9 expected under "
    "independence (Fisher OR 1.06, P=0.83). The seven SuSiE-untestable genes are drawn separately and are not counted "
    "as Neither.",
    "Redrawn 3G uses the 3,769 SuSiE-tested genes and draws the untestable seven apart.", "redraw",
    ["legend_changes:3G"], f"{CAND}/fig3/source_tables/fig3g.tsv (universe rows)")

add("R-fig2-3G", D + "fig2.md",
    "Across the 3,776 genes with finite SuSiE PP.H4 in the installed Figure 3G interface, PP.H4 and the absolute "
    "disease-state log2 fold change have Spearman rho=0.0014 (nominal P=0.93);",
    "Across the 3,769 SuSiE-tested genes drawn in Figure 3G, PP.H4 and the absolute disease-state log2 fold change "
    "have Spearman rho=0.001 (nominal P=0.95; 0.0014 and P=0.93 when the seven SuSiE-untestable genes are included);",
    "Redrawn 3G universe.", "redraw", ["legend_changes:3G"], f"{CAND}/fig3/source_tables/fig3g.tsv (universe rows)")

add("R-fig4-note", D + "fig4.md",
    "> Figure 4A and Figure 4F were installed on 2026-10-14; Figure 4C, 4E, 4G and > 4H were installed on 2026-10-15 "
    "and match the text below.",
    "> Figures 4B to 4H were redrawn and installed on 2026-10-17 and match the text\n> below; the content of Figure "
    "4A is now in the legend.",
    "Install note for the redrawn Figure 4 panels.", "redraw", ["install_log"],
    "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/install-20261017T102551Z/install_log.json")

add("R-fig4-4A", D + "fig4.md",
    "Seven datasets comprising 1,232,318 cells from 104 descriptive biological donors were integrated into one "
    "annotated liver atlas (Figure 4A); 102 donors remain in the analyzed atlas after two prespecified low-cell "
    "library exclusions.",
    "Seven datasets comprising 275 libraries and 1,232,318 cells were integrated into one annotated liver atlas; 273 "
    "libraries (1,232,285 cells) from 102 donors remain in the analyzed atlas after two prespecified low-cell library "
    "exclusions.",
    "4A moved to the legend; the 104-donor count counted the two excluded GSE189600 libraries as donors.", "redraw",
    ["legend_changes:4A"], "Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/fig4/recompute/r4a_atlas_library_counts.tsv; "
    f"{CAND}/fig4/source_tables/4B.tsv")

add("R-fig4-4B", D + "fig4.md",
    "The integrated UMAP identifies cell types with a color legend (Figure 4B). Disease-density surfaces are "
    "estimated per biological donor after sequencing runs are mapped to donor identity, then averaged within each "
    "disease group. A separate density scale shows the disease comparison. The UMAP is descriptive",
    "The integrated UMAP shows 50,000 cells drawn uniformly from the analyzed cells, coloured by cell type (Figure "
    "4B). The donor-averaged disease density is not drawn, because its Healthy and Steatohepatitis groups come from "
    "different dataset mixes. The UMAP is descriptive",
    "4B redraw removed the disease-density overlay.", "redraw", ["legend_changes:4B"], f"{CAND}/fig4/README.md (4B)")

add("R-fig4-4C", D + "fig4.md",
    "Figure 4C shows the three positive and three negative program effects with the smallest 117-program permutation "
    "BH q per lineage, with all 117 programs retained in the supplementary view.",
    "Figure 4C shows the ordinal stage effect of all 117 programs, ranked within each of the five lineages.",
    "4C redraw draws all 117 programs.", "redraw", ["legend_changes:4C"], f"{CAND}/fig4/source_tables/4C.tsv")

add("R-fig4-4C-ecm", D + "fig4.md",
    " Stromal ECM is no longer among the rows shown in Figure 4C.",
    "",
    "4C redraw draws all 117 programs, so Stromal ECM is shown again.", "redraw", ["legend_changes:4C"],
    f"{CAND}/fig4/source_tables/4C.tsv")

add("R-fig4-4E", D + "fig4.md",
    "Figure 4E shows linear fibrosis-category and NAS associations for all 112 eligible weighted programs, with BH "
    "correction over 224 tests,",
    "Figure 4E shows linear fibrosis-category associations for all 112 eligible weighted programs (NAS associations "
    "are in Figure S4), with BH correction over 224 tests on both axes,",
    "4E redraw draws the fibrosis axis only.", "redraw", ["legend_changes:4E"], f"{CAND}/fig4/source_tables/4E.tsv")

add("R-fig5-5A", D + "fig5.md",
    "each in its own units (Figure 5A). Model specifications",
    "each in its own units. Model specifications",
    "5A moved to the legend.", "redraw", ["legend_changes:5A"], f"{CAND}/fig5/README.md (5A)")

add("R-disc-6A", D + "fig5_discussion.md",
    "and deterministic next-experiment rule (Figure 6A).",
    "and deterministic next-experiment rule (Figure 6 legend).",
    "6A moved to the legend.", "redraw", ["legend_changes:6A"], f"{CAND}/fig6/legend_changes.md (6A)")

add("R-disc-6E", D + "fig5_discussion.md",
    "Figure 6E closes the Resource loop by linking each displayed Catalog entry to its",
    "The deposited Catalog (Data availability) links each displayed Catalog entry to its",
    "6E moved to Data availability.", "redraw", ["legend_changes:6E"], f"{CAND}/fig6/legend_changes.md (6E)")

add("R-fig7-7CF", D + "fig7.md",
    "Its pooled broad cell-class mapping macro F1 (0.9389) includes both, and the marker records that mixed "
    "exposure; all displayed gains and intervals remain all-study estimates.",
    "Its pooled broad cell-class mapping macro F1 (0.9389) includes both. Figures 7C and 7F draw each encoder's score "
    "and gain twice, on all seven held-out studies and on that encoder's clean studies; on its five clean studies "
    "TranscriptFormer scores 0.963 with the same two-layer head.",
    "7C/7F redraws add clean-only marks.", "redraw", ["legend_changes:7C", "legend_changes:7F"],
    f"{CAND}/fig7/source_tables/fig7C.tsv")

add("R-IDX-7C", IDX,
    "Figure 7C is broad cell-class mapping; its all-study score includes two TranscriptFormer training-seen studies. "
    "The five-clean-study comparison against inductive scVI is reported separately.",
    "Figure 7C is broad cell-class mapping, drawn on all seven held-out studies and on each encoder's clean studies; "
    "the all-study TranscriptFormer score includes two training-seen studies. Figure 7F draws the clean-study "
    "comparison against the scVI encoder trained here.",
    "7C/7F redraws.", "redraw", ["legend_changes:7C", "legend_changes:7F"], f"{CAND}/fig7/source_tables/fig7F.tsv")


# --------------------------------------------------------------------------
# docs/PAPER.md (team-lead 2026-10-17): 4H donor-first (decision 7), 5F gene-level null (decision 5)
# --------------------------------------------------------------------------
PAP = "docs/PAPER.md"
add("P-4H", PAP,
    "Spearman rho is 0.63 and 32/43 effects share direction (0.7265 and 37/43 when the controls entered as stage 0);",
    "On the donor-first estimates of the system inference tables, Spearman rho is 0.67 (bootstrap 95% interval 0.38 "
    "to 0.84) and 34/43 effects share direction (feature-first summaries: 0.63 and 32/43; 0.7265 and 37/43 when the "
    "controls entered as stage 0);",
    "Decision 7: 4H plots the donor-first estimand; worded as D7-fig4-a.", "7", ["REPORT.md decision 7"],
    f"{CAND}/fig4/source_tables/4H.tsv")
add("P-5F-results", PAP,
    "Stromal ECM exceeded its matched-gene null in both sources (excess 0.040, q=0.001; 0.050, q=0.0006) and Ductular "
    "injury did not (0.015, q=0.079; 0.003, q=0.32);",
    "Stromal ECM exceeded its matched-gene null in both sources (excess 0.040, q=0.001; 0.050, q=0.0006) and Ductular "
    "injury did not (0.015, q=0.079; 0.003, q=0.32). That null resamples genes with donors and arrays held fixed, so "
    "it is not a donor-level test: 3 of 4 GSE192741 donors have positive Stromal ECM excess (exact one-sided "
    "sign-flip P=0.19), and the 10 positive Vu arrays have no public donor key;",
    "Decision 5: the both-sources result is gene-level; print the donor-level result.", "5", ["REPORT.md decision 5"],
    f"{CAND}/fig5/source_tables/5F.tsv (n_units_positive, signflip_p)")
add("P-5F-claims", PAP,
    "Stromal ECM exceeds its matched-gene null in both Visium sources and Ductular injury does not;",
    "Stromal ECM exceeds its gene-level matched null in both Visium sources and Ductular injury does not; at donor "
    "level, 3 of 4 GSE192741 donors are positive (sign-flip P=0.19) and Vu arrays have no donor key, so neither "
    "source gives a donor-level test;",
    "Decision 5.", "5", ["REPORT.md decision 5"], f"{CAND}/fig5/source_tables/5F.tsv (n_units_positive, signflip_p)")

# --------------------------------------------------------------------------
# Merge with the 64 staged edits
# --------------------------------------------------------------------------
old64 = json.load(open(os.path.join(TF, "edits.json")))
check = {}
with open(os.path.join(TF, "check.tsv")) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    for line in f:
        row = dict(zip(hdr, line.rstrip("\n").split("\t")))
        check[row["id"]] = row["source_path"]

superseded = {sid for e in NEW for sid in e["source_ids"] if re.fullmatch(r"[TLIN]-\d+[ab]?", sid)}
kept = []
for e in old64:
    if e["id"] in superseded:
        continue
    kept.append(dict(id=e["id"], file=e["file"], old=e["old"], new=e["new"], reason=e["reason"],
                     decision="text_fixes", source_ids=[e["id"]], producing_file=check.get(e["id"], e["source"])))

allx = kept + NEW
# Order each file's edits by position of `old` in the live file, so the list applies top to bottom.
for e in allx:
    e["_pos"] = live(e["file"]).find(e["old"])
    assert e["_pos"] >= 0, e["id"]
allx.sort(key=lambda e: (e["file"], e["_pos"]))

# Overlap check: no two edits of one file may touch overlapping spans of the live text.
by_file = defaultdict(list)
for e in allx:
    by_file[e["file"]].append((e["_pos"], e["_pos"] + len(e["old"]), e["id"]))
problems = []
for f, spans in by_file.items():
    spans.sort()
    for (a0, a1, ai), (b0, b1, bi) in zip(spans, spans[1:]):
        if b0 < a1:
            problems.append(f"overlap in {f}: {ai} [{a0},{a1}) and {bi} [{b0},{b1})")
for e in allx:
    if live(e["file"]).count(e["old"]) != 1:
        problems.append(f"{e['id']}: old occurs {live(e['file']).count(e['old'])} times in {e['file']}")
for e in allx:
    del e["_pos"]

json.dump(allx, open(os.path.join(OUT, "edits.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
print(f"{len(allx)} edits ({len(kept)} kept from text_fixes, {len(NEW)} new or merged); "
      f"superseded from text_fixes: {sorted(superseded)}")
print("problems:", problems if problems else "none")
