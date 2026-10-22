#!/usr/bin/env python
"""Round-2 fig2d edit list: every printed Figure 2D / S2S / S2V / Note 4 value at genome-wide leads (D2).
Checks each `old` occurs exactly once in RUN/staged/round1/<file> (or the repo source when no copy exists),
that edits in one file do not overlap, and writes RUN/staged/main_fig2d.json."""
import json
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
RUN = ROOT / "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z"
AN2 = "figures/candidates/fig2d-gws-leads-20261021T124551Z"
R2 = "figures/candidates/fig2d-gws-leads-r2-20261021T173128Z"
AN2R = "Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/analyses/AN2/RESULT.md"
V = f"{AN2}/values/gws_p"
W = f"{R2}/values/gws_p"
FIG2 = "docs/manuscript/draft/fig2.md"
LEG = "docs/manuscript/05_figure_legends.md"
SN = "docs/manuscript/draft/supplementary_notes.md"
MET = "docs/manuscript/METHODS.md"
E = []


def add(id_, file, old, new, items, sources, reason, decision="D2"):
    E.append(dict(id=id_, file=file, old=old, new=new, items=items, number_sources=sources,
                  decision=decision, reason=reason))


# ------------------------------------------------------------------ fig2.md
add("fig2d-01", FIG2,
    "Of 6,467 study-by-locus credible-set instances, 6,415 passed the prespecified PIP-sum criterion.",
    "Of 6,467 study-by-locus credible-set instances, 6,415 passed the prespecified PIP-sum criterion. Figure 2D uses "
    "the 4,107 of these whose lead variant has p < 5 × 10⁻⁸ in its own GWAS (164 direct-trait and 3,943 enzyme "
    "credible sets; Methods).",
    ["B-G4.1"], [f"{V}/summary.json:n_sets", f"{V}/composition.tsv:n_credible_sets"],
    "States the genome-wide lead filter and the analyzed set size.")
add("fig2d-02", FIG2,
    "Direct MASLD and liver-fat loci then allocated 23.15% of the probability to protein-altering variants, 0.20% to "
    "canonical splice variants, 3.16% to synonymous or UTR variants, and 73.50% to other noncoding variants. "
    "Liver-enzyme loci allocated 3.87%, 0.00%, 3.60%, and 92.52% (Figure 2D).",
    "Direct MASLD and liver-fat loci then allocated 29.42% of the probability to protein-altering variants, 0.23% to "
    "canonical splice variants, 3.92% to synonymous or UTR variants, and 66.43% to other noncoding variants. "
    "Liver-enzyme loci allocated 6.67%, 0.00%, 3.50%, and 89.83% (Figure 2D).",
    ["B-G4.1"], [f"{V}/composition.tsv:locus_weighted_fraction"], "Four-class composition at genome-wide leads.")
add("fig2d-03", FIG2,
    "protein-altering variants carried 15.77% of the posterior for direct traits and 3.26% for liver enzymes.",
    "protein-altering variants carried 17.36% of the posterior for direct traits and 3.43% for liver enzymes.",
    ["B-G4.1"], [f"{V}/composition.tsv:instance_weighted_fraction"], "Set-weighted fractions at genome-wide leads.")
add("fig2d-04", FIG2,
    "In the Admixed American (AMR) GWAS, 10 of 37 direct and 84 of 271 enzyme lead windows have a credible set.",
    "In the Admixed American (AMR) GWAS, all 4 direct and all 35 enzyme lead windows at p < 5 × 10⁻⁸ have a "
    "credible set.",
    ["B-G4.1"], [f"{V}/G7_coverage.tsv:AMR rows"], "AMR coverage counted over genome-wide lead windows.")
add("fig2d-05", FIG2,
    "carry 23.34% and 3.87% of their posterior probability, respectively (difference 0.1947, 95% physical-cluster "
    "bootstrap interval 0.0692 to 0.3399; 233 physical clusters: 19 direct, 228 enzyme, and 14 shared).",
    "carry 29.65% and 6.67% of their posterior probability, respectively (difference 0.2298, 95% physical-cluster "
    "bootstrap interval 0.0808 to 0.3958; 120 physical clusters: 16 direct, 115 enzyme, and 11 shared).",
    ["B-G4.1"], [f"{V}/pooled_coding.tsv"], "Pooled coding fraction and interval at genome-wide leads.")
add("fig2d-06", FIG2,
    "Of the direct credible sets, 106 of 207 are European, and 15 of the 19 direct clusters contain a European "
    "direct set, whereas 6,126 of 6,208 enzyme credible sets are non-European.",
    "Of the direct credible sets, 106 of 164 are European, and 15 of the 16 direct clusters contain a European "
    "direct set, whereas 3,861 of 3,943 enzyme credible sets are non-European.",
    ["B-G4.1"], [f"{V}/summary.json:n_EUR_direct,n_nonEUR_enzyme", f"{AN2R}:fig2.md 37 row (15 of 16)"],
    "Ancestry make-up of the filtered arms.")
add("fig2d-07", FIG2,
    "To separate trait class from ancestry, we stratified the comparison by ancestry (Figure 2D, ancestry panel; "
    "10,000 shared draws over 81 pooled clusters). European studies gave 29.88% against 6.92% (difference 0.2296, "
    "95% interval 0.0550 to 0.4236), and AMR studies gave 6.73% against 4.54% (0.0219, −0.0489 to 0.1065). The "
    "prespecified fixed-weight European-plus-AMR estimate (weights 16/26 and 10/26) was 0.1497 (0.0388 to 0.2718). "
    "All non-European studies pooled gave 6.99% against 3.68% (0.0331, −0.0383 to 0.1274), a descriptive value with "
    "no calibration control. None of these estimates receives a trait-class call. The European estimate failed "
    "three checks: coverage (enzyme credible sets at 8 of 767 lead windows by locus identifier and 23 by position, "
    "direct arm at 47 of 47), power (minimum detectable difference 0.304 against the pooled 0.195), and its "
    "random-split calibration control. In this control, we split East Asian enzyme clusters (East Asian and African "
    "clusters for the European-plus-AMR estimate) at random into two pseudo-classes with the observed cluster "
    "counts, 200 times. Both pseudo-classes come from one trait class, so a calibrated 95% interval should cover "
    "zero in at least 180 splits, and the 98.333% interval (Bonferroni over the three estimates) in at least 190. "
    "The European 98.333% interval covered zero in 188 splits. The AMR and European-plus-AMR estimates also failed "
    "this control (141 and 173 splits at 95%; 150 and 184 at 98.333%), and the European-plus-AMR estimate also "
    "carries the European coverage failure. In these failed controls the intervals excluded zero more often than a "
    "calibrated interval would, so the AMR and European-plus-AMR intervals are too narrow.",
    "To separate trait class from ancestry, we stratified the comparison by ancestry (Figure 2D, ancestry panel; "
    "10,000 shared draws over 44 pooled clusters). European studies gave 29.88% against 6.92% (difference 0.2296, "
    "95% interval 0.0563 to 0.4211), and AMR studies gave 8.11% against 6.91% (0.0120, −0.1246 to 0.1620). The "
    "prespecified fixed-weight European-plus-AMR estimate (weights 16/20 and 4/20, the direct-cluster counts) was "
    "0.1861 (0.0503 to 0.3337). All non-European studies pooled gave 8.71% against 6.48% (0.0223, −0.0831 to "
    "0.1686), a descriptive value with no calibration control. None of these estimates receives a trait-class call. "
    "The European estimate failed two checks: coverage (enzyme credible sets at 8 of 767 lead windows by locus "
    "identifier and 23 by position, direct arm at 47 of 47) and power (minimum detectable difference 0.303 against "
    "the pooled 0.230). It passed its random-split calibration control. In this control, we split East Asian enzyme "
    "clusters (East Asian and African clusters for the European-plus-AMR estimate) at random into two pseudo-classes "
    "with the observed cluster counts, 200 times. Both pseudo-classes come from one trait class, so a calibrated 95% "
    "interval should cover zero in at least 180 splits, and the 98.333% interval (Bonferroni over the three "
    "estimates) in at least 190. The European intervals covered zero in 183 and 191 splits. AMR coverage passed (4 "
    "of 4 direct and 35 of 35 enzyme lead windows), but the AMR and European-plus-AMR estimates failed the control "
    "(122 and 172 splits at 95%; 135 and 188 at 98.333%), so their intervals are too narrow. Both also narrowly "
    "failed the power check (minimum detectable differences 0.234 and 0.232), and in 156 of the 10,000 draws the "
    "AMR direct arm, 4 clusters, was empty, more than the prespecified limit of 100. The European-plus-AMR estimate "
    "also carries the European coverage failure.",
    ["B-G4.1"],
    [f"{V}/family_C.tsv:C1-C3", f"{V}/summary.json:family_C_design(units 44, weights 0.8/0.2, clusters)",
     f"{V}/descriptive_refits.tsv:nonEUR_pooled", f"{V}/N1_coverage.tsv", f"{V}/G7_coverage.tsv",
     f"{V}/family_C.tsv:n_finite 9844 (156 dropped)"],
    "Ancestry estimates, checks and call reasons at genome-wide leads; C3 also fails G2 (0.2337 > 0.2298), "
    "which AN2's call list omitted.")
add("fig2d-08", FIG2,
    "carry 40.14% against 3.87% for all enzyme traits (difference 0.3626, 95% interval 0.1835 to 0.5455).",
    "carry 40.14% against 6.67% for all enzyme traits (difference 0.3347, 95% interval 0.1492 to 0.5280).",
    ["B-G4.1"], [f"{V}/subclass_split.tsv:PDFF/ENZ"], "PDFF subclass row at genome-wide leads.")
add("fig2d-09", FIG2,
    "NAFLD plus NASH (170 sets, 13 clusters) carry 12.44% (0.0857, −0.0164 to 0.2096) (Figure S2S).",
    "NAFLD plus NASH (127 sets, 8 clusters) carry 15.92% (0.0925, −0.0596 to 0.2844) (Figure S2S).",
    ["B-G4.1"], [f"{V}/subclass_split.tsv:NAFLD_NASH/ENZ"], "NAFLD+NASH subclass row at genome-wide leads.")
add("fig2d-10", FIG2,
    "Five of the 19 direct clusters carry 76.0% of the summed direct cluster fractions, and the largest, one PDFF "
    "credible set at chr11:823,586, carries 22.5%. The 57 robustness comparisons (study and direct-cluster "
    "omissions, five grouping windows, and one LD-block comparison) gave ratios from 4.37 to 8.28 around the pooled "
    "reference ratio of 6.03, a descriptive value.",
    "Five of the 16 direct clusters carry 73.9% of the summed direct cluster fractions, and the largest, one PDFF "
    "credible set at chr11:823,586, carries 21.1%. The 50 robustness comparisons (study and direct-cluster "
    "omissions, five grouping windows, and one LD-block comparison) gave ratios from 3.44 to 5.18 around the pooled "
    "reference ratio of 4.45, a descriptive value.",
    ["B-G4.1"], [f"{V}/summary.json:concentration", f"{V}/summary.json:robustness", f"{V}/robustness.tsv"],
    "Concentration and robustness at genome-wide leads; 50 counts comparisons without the reference fit.")
add("fig2d-11", FIG2,
    "The two trait classes in the 23.34% against 3.87% comparison also differ in GWAS size and effect size. "
    "Direct-trait credible sets come from smaller GWAS than enzyme credible sets (median effective n 21,893 against "
    "160,000 over all 6,467 credible sets) and carry larger standardized lead-variant effects (median 0.0368 against "
    "0.0183). Supplementary Note 4 reports a sensitivity analysis matched on standardized lead-variant effect, over "
    "all 6,467 credible sets including those that failed the PIP-sum criterion.",
    "The two trait classes in the 29.65% against 6.67% comparison also differ in GWAS size and effect size. "
    "Direct-trait credible sets come from smaller GWAS than enzyme credible sets (median effective n 21,893 against "
    "160,000 over the 4,134 credible sets at genome-wide leads) and carry larger standardized lead-variant effects "
    "(median 0.0423 against 0.0190). Supplementary Note 4 reports a sensitivity analysis matched on standardized "
    "lead-variant effect, over the same 4,134 credible sets, including those that failed the PIP-sum criterion.",
    ["B-G4.1"], [f"{V}/ascertainment_medians.tsv", f"{W}/note4_medians.tsv", f"{W}/summary.json:note4_sets"],
    "Ascertainment medians at genome-wide leads.")
add("fig2d-12", FIG2,
    "Summed over all credible sets of each class, 13.84% of direct-trait and 14.67% of liver-enzyme noncoding "
    "posterior probability lies within 2 kb of a TSS, 9.26% and 8.44% in liver-lineage open chromatin, 0.43% and "
    "3.70% in a liver ABC enhancer, and 77.90% and 77.69% in sequence with none of these contexts. These pooled "
    "values mix ancestries. AMR is the one ancestry with comparable fine-mapping coverage in both classes. Within "
    "AMR, we detected no difference in the share of noncoding probability with no deposited context (direct minus "
    "enzyme −0.0957, 99% interval −0.344 to 0.126; Bonferroni over five contexts). Promoter, ABC, and "
    "accessible-chromatin differences received no call because their random-split calibration controls failed.",
    "Summed over all credible sets of each class at genome-wide leads, 14.93% of direct-trait and 14.01% of "
    "liver-enzyme noncoding posterior probability lies within 2 kb of a TSS, 10.73% and 7.84% in liver-lineage open "
    "chromatin, 0.30% and 3.85% in a liver ABC enhancer, and 75.39% and 78.28% in sequence with none of these "
    "contexts. These pooled values mix ancestries. AMR is the one ancestry with comparable fine-mapping coverage in "
    "both classes. Within AMR, no context difference received a call (99% intervals, Bonferroni over five "
    "contexts). All four random-split calibration controls failed, and in 167 of the 10,000 draws the AMR direct "
    "arm, 4 clusters, was empty, more than the prespecified limit of 100. For the share with no deposited context, "
    "direct minus enzyme was −0.209 (99% interval −0.599 to 0.166).",
    ["B-G4.1"], [f"{W}/context_cs_weighted.tsv", f"{W}/family_X.tsv", f"{W}/N1x_coverage.tsv",
                 f"{W}/summary.json:family_x.amr_context_clusters"],
    "Context fractions and family X rerun at genome-wide leads; X4 loses its call (G3 175/185, 167 dropped draws).")

# ------------------------------------------------------------------ legends
add("fig2d-13", LEG,
    "for direct MASLD/liver-fat traits (207 credible sets, 19 clusters) and liver-enzyme traits (6,208 sets, 228 "
    "clusters). Inset: direct minus enzyme protein-altering plus splice probability, 19.5 percentage points (95% "
    "cluster bootstrap 6.9 to 34.0).",
    "for direct MASLD/liver-fat traits (164 credible sets, 16 clusters) and liver-enzyme traits (3,943 sets, 115 "
    "clusters); only credible sets whose lead variant has p < 5 × 10⁻⁸ are included. Inset: direct minus enzyme "
    "protein-altering plus splice probability, 23.0 percentage points (95% cluster bootstrap 8.1 to 39.6).",
    ["B-G4.1"], [f"{R2}/panels/check/Fig2D_coding_mass_by_trait_class.txt", f"{V}/pooled_coding.tsv"],
    "2D legend matches the redrawn panel.")
add("fig2d-14", LEG,
    "**(S)** Coding plus canonical-splice share of fine-mapped posterior (Figure 2D), split within the direct class.",
    "**(S)** Coding plus canonical-splice share of fine-mapped posterior (Figure 2D; credible sets at leads with "
    "p < 5 × 10⁻⁸), split within the direct class.",
    ["B-G4.1"], [], "S2S uses the Figure 2D set.")
add("fig2d-15", LEG,
    "as a share of noncoding posterior summed over all credible sets of a class. Contexts overlap and are not "
    "summed. Pooled values mix ancestries. Within AMR, only the no-context difference passed its calibration "
    "control, and no difference was detected.",
    "as a share of noncoding posterior summed over all credible sets of a class at leads with p < 5 × 10⁻⁸. "
    "Contexts overlap and are not summed. Pooled values mix ancestries. Within AMR, no context difference received "
    "a call, because every calibration control failed.",
    ["B-G4.1"], [f"{W}/family_X.tsv:call"], "S2V legend: X4 no longer receives a call.")

# ------------------------------------------------------------------ Supplementary Note 4
add("fig2d-16", SN,
    "credible sets for direct MASLD and liver-fat traits carry 23.34% of their posterior probability in "
    "protein-altering or canonical-splice variants, against 3.87% for the liver enzymes ALT, AST and GGT. The "
    "difference is 0.1947 (95% physical-cluster bootstrap interval 0.0692 to 0.3399).",
    "credible sets for direct MASLD and liver-fat traits carry 29.65% of their posterior probability in "
    "protein-altering or canonical-splice variants, against 6.67% for the liver enzymes ALT, AST and GGT. The "
    "difference is 0.2298 (95% physical-cluster bootstrap interval 0.0808 to 0.3958). Only the 4,107 reliable "
    "credible sets whose lead variant has p < 5 × 10⁻⁸ enter these estimates.",
    ["B-G4.1"], [f"{V}/pooled_coding.tsv", f"{V}/summary.json:n_sets"], "Note 4 pooled values at genome-wide leads.")
add("fig2d-17", SN,
    "Within European studies, the fractions are 29.88% and 6.92%, and within admixed American (AMR) studies, 6.73% "
    "and 4.54%. No stratum supports a calibrated trait-class call. The European estimate failed the coverage check "
    "(UKBB enzyme credible sets lie at 8 of 767 lead windows by locus identifier and within 500 kb of 23) and the "
    "power check. The European, AMR and combined estimates all failed their random-split calibration controls. In "
    "these controls the intervals excluded zero more often than a calibrated interval would. The European interval "
    "failed only at the Bonferroni-corrected 98.333% level, but the AMR and combined intervals failed already at "
    "95%, so these two intervals are too narrow. The pooled direct-minus-enzyme difference comes mainly from the "
    "three PDFF GWAS and 5 of the 19 direct clusters. The ancestry-stratified estimates used 10,000 shared draws "
    "over 81 pooled clusters (seed 20260923). The 57 robustness comparisons gave ratios from 4.37 to 8.28.",
    "Within European studies, the fractions are 29.88% and 6.92%, and within admixed American (AMR) studies, 8.11% "
    "and 6.91%. No stratum supports a calibrated trait-class call. The European estimate failed the coverage check "
    "(UKBB enzyme credible sets lie at 8 of 767 lead windows by locus identifier and within 500 kb of 23) and the "
    "power check, and passed its random-split calibration control. The AMR and combined estimates failed their "
    "random-split calibration controls at 95% and at the Bonferroni-corrected 98.333% level, so their intervals "
    "are too narrow. Both also narrowly failed the power check, and the AMR direct arm, 4 clusters, was empty in "
    "156 of 10,000 draws, above the prespecified limit of 100. The pooled direct-minus-enzyme difference comes "
    "mainly from the three PDFF GWAS and 5 of the 16 direct clusters. The ancestry-stratified estimates used 10,000 "
    "shared draws over 44 pooled clusters (seed 20260923). The 50 robustness comparisons gave ratios from 3.44 to "
    "5.18.",
    ["B-G4.1"], [f"{V}/family_C.tsv", f"{V}/N1_coverage.tsv", f"{V}/summary.json:family_C_design,robustness"],
    "Note 4 stratified results and call reasons at genome-wide leads.")
add("fig2d-18", SN,
    "In a sensitivity analysis, we used all 6,467 credible sets, including those that failed the PIP-sum criterion "
    "(seed 20260817). Matching on standardized lead-variant effect gave a direct-to-enzyme ratio of 5.78 (95% "
    "interval 1.99 to 48.05; 27 matched enzyme loci). Direct credible sets come from smaller GWAS (median effective "
    "n 21,893 against 160,000, over the 208 direct and 6,259 enzyme credible sets in this analysis) and carry larger "
    "standardized effects (median 0.0368 against 0.0183).",
    "In a sensitivity analysis, we used all 4,134 credible sets at genome-wide leads, including those that failed "
    "the PIP-sum criterion (seed 20260817). Matching on standardized lead-variant effect gave a direct-to-enzyme "
    "ratio of 3.80 (95% interval 1.29 to 39.62; 19 matched enzyme loci). Direct credible sets come from smaller "
    "GWAS (median effective n 21,893 against 160,000, over the 165 direct and 3,969 enzyme credible sets in this "
    "analysis) and carry larger standardized effects (median 0.0423 against 0.0190).",
    ["B-G4.1"], [f"{W}/note4_effect_matched.tsv:matched_on_standardized_effect", f"{W}/note4_medians.tsv",
                 f"{W}/summary.json:note4_sets"],
    "Effect-matched sensitivity rerun at genome-wide leads (job 22343298).")

# ------------------------------------------------------------------ intro, discussion
add("fig2d-19", "docs/manuscript/draft/01_intro.md",
    "Pooled across ancestries, fine-mapped credible sets for direct MASLD and liver-fat phenotypes carry 23.34% of "
    "their posterior probability in protein-altering or canonical-splice variants. The value is 3.87% for the liver "
    "enzymes",
    "Pooled across ancestries, fine-mapped credible sets at genome-wide-significant leads for direct MASLD and "
    "liver-fat phenotypes carry 29.65% of their posterior probability in protein-altering or canonical-splice "
    "variants. The value is 6.67% for the liver enzymes",
    ["B-G4.1"], [f"{V}/pooled_coding.tsv"], "Intro pooled values at genome-wide leads.")
add("fig2d-20", "docs/manuscript/draft/fig5_discussion.md",
    "protein-altering plus canonical-splice variants carry 23.34% of fine-mapped\nposterior probability for direct "
    "MASLD and liver-fat traits and 3.87% for enzyme\ntraits.",
    "protein-altering plus canonical-splice variants carry 29.65% of fine-mapped\nposterior probability for direct "
    "MASLD and liver-fat traits and 6.67% for enzyme\ntraits, counting credible sets at genome-wide-significant "
    "leads.",
    ["B-G4.1"], [f"{V}/pooled_coding.tsv"], "Discussion pooled values at genome-wide leads.")

# ------------------------------------------------------------------ Methods
add("fig2d-21", MET,
    "Of the 6,415 reliable credible sets, 2,363 lie at suggestive non-European leads (43 of 207 direct-trait sets "
    "and 2,320 of 6,208 enzyme sets).",
    "Of the 6,415 reliable credible sets, 6,227 come from non-European strata, and 2,363 of these lie at leads "
    "labeled suggestive (43 of 207 direct-trait sets and 2,320 of 6,208 enzyme sets).",
    ["B-G4.1"], [f"{AN2R}:Non-European reliable set count: 6,227"], "Prints the corrected 6,227 non-European count.")
add("fig2d-22", MET,
    "Figure 2D and its ancestry rows used only credible sets whose lead variant was labeled genome-wide "
    "significant (p < 5 × 10⁻⁸; 3.24 × 10⁻⁸ in AFR strata).",
    "Figure 2D, its ancestry rows, Figure S2S, the noncoding-context comparison and the Supplementary Note 4 "
    "sensitivity used only credible sets whose lead variant has p < 5 × 10⁻⁸ in its own GWAS, one threshold in "
    "every stratum: 4,107 of the 6,415 reliable credible sets (164 direct-trait and 3,943 enzyme). This keeps 55 "
    "AFR enzyme sets whose leads have p between 3.24 × 10⁻⁸ and 5 × 10⁻⁸ and are labeled suggestive. For 48 "
    "European sets, the restored lead position has no exact match in the current lead table. They were kept "
    "because each has a sentinel p < 5 × 10⁻⁸ in its own summary statistics.",
    ["B-G4.1"], [f"{AN2R}:Answer (55 AFR sets); Handling of the 48 unmatched European sets", f"{V}/summary.json:n_sets"],
    "Methods states the D2 rule (p < 5e-8 in every stratum) instead of the tier label.")
add("fig2d-23", MET, "(weights 16/26 and 10/26)",
    "(weights 16/20 and 4/20, each stratum's share of direct clusters)",
    ["B-G4.1"], [f"{V}/summary.json:family_C_design.weights,clusters"], "Fixed weights at genome-wide leads.")
add("fig2d-24", MET, "draws over the 81 pooled clusters holding EUR or AMR credible sets",
    "draws over the 44 pooled clusters holding EUR or AMR credible sets",
    ["B-G4.1"], [f"{V}/summary.json:family_C_design.units"], "Resampled units at genome-wide leads.")
add("fig2d-25", MET, "was at most the pooled difference (0.195).",
    "was at most the pooled difference (0.230).",
    ["B-G4.1"], [f"{V}/pooled_coding.tsv:difference"], "Power threshold at genome-wide leads.")
add("fig2d-26", MET,
    "and the 98.333% interval in at least 190. A difference was called only if",
    "and the 98.333% interval in at least 190. A member with more than 100 of the 10,000 draws in which a required "
    "arm was empty received no call. A difference was called only if",
    ["B-G4.1"], [], "States the prespecified dropped-draw limit now cited in the Results (PRESPEC_H58.md:292).")

# ------------------------------------------------------------------ PAPER.md
add("fig2d-27", "docs/PAPER.md",
    "(H58,\n  2026-10-15: 106 of 207 direct but 82 of 6,208 enzyme credible sets are\n  European; EUR 29.88% vs 6.92%, "
    "AMR 6.73% vs 4.54%; no call in any stratum;\n  the pooled excess is carried mainly by three PDFF GWAS). Pooled "
    "values:\n  23.3428% against 3.8732%, a 6.026671-fold ratio after\n  restoring 715 truncated positions among 6,415 "
    "gated credible sets. The\n  pooled single-linkage clustering of lead positions within 1 Mb yields 233 clusters: "
    "19 direct, 228 enzyme and 14\n  shared.",
    "(H58,\n  2026-10-15, recomputed 2026-10-21 on the 4,107 gated credible sets whose lead has\n  p < 5e-8 in every "
    "stratum: 106 of 164 direct but 82 of 3,943 enzyme credible sets are\n  European; EUR 29.88% vs 6.92%, AMR 8.11% "
    "vs 6.91%; no call in any stratum;\n  the pooled excess is carried mainly by three PDFF GWAS). Pooled values:\n  "
    "29.65% against 6.67%, a 4.45-fold ratio (physical-cluster bootstrap\n  2.05–9.88), after restoring 715 truncated "
    "positions among 6,415 gated credible sets. The\n  pooled single-linkage clustering of lead positions within 1 Mb "
    "yields 120 clusters: 16 direct, 115 enzyme and 11\n  shared.",
    ["B-G4.1"], [f"{V}/summary.json:pooled", f"{V}/summary.json:n_sets,n_EUR_direct", f"{V}/family_C.tsv:C2,C3"],
    "PAPER.md claim record at genome-wide leads.")
add("fig2d-28", "docs/PAPER.md",
    "The adopted October 13 corrected\n  rerun gives a physical-cluster bootstrap interval of 2.57459–13.51016 for\n  the ratio. Its complete robustness menu has 57 comparisons: one reference,\n  13 direct-study omissions, 18 enzyme-study omissions, 19 direct-cluster\n  omissions, five clustering windows and one LD-block comparison. Restoring\n  19 direct clusters instead of the historical 24 explains the change from\n  62 comparisons; every comparison type was retained. Ratios range from\n  4.36813 to 8.27569, and every difference interval has a positive lower bound\n  (minimum 0.02237). None of the 57 comparisons",
    "At genome-wide leads, the\n  robustness menu has 51 rows: one reference, 12 direct-study omissions,\n  16 enzyme-study omissions, 16 direct-cluster omissions, five clustering\n  windows and one LD-block comparison. Over the 50 non-reference comparisons,\n  ratios range from 3.44 to 5.18, and every difference interval has a positive\n  lower bound (minimum 0.0211). None of these comparisons",
    ["B-G4.1"], [f"{V}/robustness.tsv (51 rows; min diff_boot_lo 0.021135)", f"{V}/summary.json:robustness.by_arm"],
    "PAPER.md robustness record at genome-wide leads; the October 13 ratio interval is replaced by the pooled 2.05 to 9.88 in fig2d-27.")
add("fig2d-29", SN,
    "It rests on 19 direct-trait regions and is pooled across ancestries.",
    "It rests on 16 direct-trait regions at genome-wide-significant leads and is pooled across ancestries.",
    ["B-G4.1"], [f"{V}/pooled_coding.tsv:n_clusters_direct"], "Direct cluster count at genome-wide leads.")


def src(file):
    p = RUN / "staged/round1" / file
    return (p if p.exists() else ROOT / file).read_text()


bad = []
spans = {}
for e in E:
    t = src(e["file"])
    n = t.count(e["old"])
    if n != 1:
        bad.append((e["id"], n))
        continue
    i = t.index(e["old"])
    spans.setdefault(e["file"], []).append((i, i + len(e["old"]), e["id"]))
    for w in ["—", "TREAT", "causal", "progression"]:
        assert w not in e["new"], (e["id"], w)
for f, s in spans.items():
    s.sort()
    for a, b in zip(s, s[1:]):
        assert a[1] <= b[0], ("overlap", f, a[2], b[2])
assert not bad, bad
panels = [
    {"candidate_pdf": f"{R2}/panels/Fig2D_coding_mass_by_trait_class.pdf",
     "install_target": "figures/main/fig2_genetics/panels/Fig2D_pip_architecture_by_trait_directness.pdf"},
    {"candidate_pdf": f"{R2}/panels/Fig2D_coding_by_ancestry.pdf",
     "install_target": "figures/main/fig2_genetics/panels/Fig2D_coding_by_ancestry.pdf"},
    {"candidate_pdf": f"{R2}/panels/FigS2S_direct_subclass_split.pdf",
     "install_target": "figures/main/fig2_genetics/panels/FigS2S_direct_subclass_split.pdf"},
    {"candidate_pdf": f"{R2}/panels/FigS2V_noncoding_context.pdf",
     "install_target": "figures/main/fig2_genetics/panels/Fig2_companion_noncoding_annotation_gap.pdf"},
]
out = dict(edits=E, panels=panels)
(RUN / "staged/main_fig2d.json").write_text(json.dumps(out, indent=1, ensure_ascii=False))
print(len(E), "edits OK; written", RUN / "staged/main_fig2d.json")
