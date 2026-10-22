# F2d repair: final text edits, written against the current files, then checked.
# Each old_text must occur exactly once in its file; two edits in one file must not overlap.
import json
import sys
from pathlib import Path

ROOT = Path("/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas")
OUT = Path(sys.argv[1])
MS = "docs/manuscript/draft/"
LG = "docs/manuscript/05_figure_legends.md"

E = {}  # fid -> list of edits


def add(fid, file, line, old, new, claim=False):
    E.setdefault(fid, []).append(dict(file=file, line=line, old_text=old, new_text=new,
                                      claim_level=claim))


# ------------------------------------------------------------------ R13-08 (H57)
add("R13-08", MS + "abstract.md", "8-10",
    "Inherited genetic variation and the molecular state of established disease\n"
    "capture different dimensions of metabolic dysfunction-associated steatotic\n"
    "liver disease (MASLD), but",
    "Inherited genetic variation and the molecular state of established disease\n"
    "have each been measured in metabolic dysfunction-associated steatotic\n"
    "liver disease (MASLD), but", True)
add("R13-08", MS + "fig2.md", "9",
    "## Trait-associated inherited variation defines a genetic map distinct from established disease-state transcription",
    "## Trait-associated inherited variation, qualified by phenotype, ancestry and fine-mapping coverage", True)
add("R13-08", MS + "fig2.md", "86-87",
    "The distinction between inherited genetic evidence and established tissue\n"
    "state has two implications. First,",
    "The colocalized and disease-state gene sets share 26 genes, against 21.2\n"
    "expected for genes of the same expression level among the 3,769 genes\n"
    "evaluated by multi-signal colocalization (odds ratio 1.29, 95% interval 0.84\n"
    "to 2.01; exact P = 0.28, one prespecified test). Expression level was not\n"
    "fully controlled, so this is a description, not evidence that the two layers\n"
    "are biologically separate or convergent. Two practical points follow. First,", True)
add("R13-08", MS + "fig2.md", "116-117",
    "colocalized liver-trait regulatory candidates. It complements, but does not\n"
    "substitute for, the established disease-state transcriptome.",
    "colocalized liver-trait regulatory candidates. It does not substitute for the\n"
    "established disease-state transcriptome; among genes evaluated by multi-signal\n"
    "colocalization, the two gene sets overlap at about the rate expected for genes\n"
    "of the same expression level.", True)

# ------------------------------------------------------------------ R02-03
add("R02-03", LG, "34-36",
    "ancestries and six liver traits, followed by 360 trait-specific fine-mapped\n"
    "loci (265 globally unique loci) and 996 COLOC genes. Brackets bind the",
    "ancestries and six liver traits, with 360 trait-specific fine-mapped loci\n"
    "(265 globally unique loci) and 996 COLOC genes counted separately. The COLOC\n"
    "gene column is not derived from the fine-mapped loci: the MVP European ALT,\n"
    "AST and NAFLD GWAS were not fine-mapped, and the UKBB ALT, AST and GGT GWAS\n"
    "were fine-mapped at 6, 9 and 15 loci, which match 8 of their 767 lead windows\n"
    "by locus identifier (23 within 500 kb by position). Of the 436 multi-signal\n"
    "genes, 100 have a colocalized lead within 250 kb of a fine-mapped window from\n"
    "the same study and 94 are supported only by the unfine-mapped MVP European\n"
    "strata (Figure 2A support companion, not assembled in canonical A-I). Brackets bind the")
add("R02-03", MS + "fig2.md", "101",
    "We therefore quantify DNA architecture directly from all reliable Tier-1/2 SuSiE credible-set instances rather than",
    "We therefore quantify DNA architecture directly from every reliable Tier-1/2 SuSiE credible set that fine-mapping produced rather than")
add("R02-03", MS + "fig2.md", "101",
    "The grouping limits repeated-study dominance and does not assert LD or statistical independence. "
    "Credible-set-weighted estimates remain a sensitivity analysis in the source table.",
    "The grouping limits repeated-study dominance and does not assert LD or statistical independence. "
    "Credible-set-weighted estimates remain a sensitivity analysis in the source table. "
    "Fine-mapping coverage differs by stratum. Every European direct-trait lead window has a credible set (47 of 47). "
    "The European enzyme arm is 82 credible sets from UKBB ALT, AST and GGT, fine-mapped at 6, 9 and 15 loci that match "
    "8 of their 767 lead windows by locus identifier (23 of 767 within 500 kb by position); it has no credible set at "
    "PNPLA3 or TM6SF2. The MVP European ALT, AST and NAFLD GWAS were not fine-mapped. In AMR, 10 of 37 direct and 84 of "
    "271 enzyme lead windows have a credible set.")

# ------------------------------------------------------------------ R02-01
add("R02-01", MS + "abstract.md", "24-32",
    "Credible sets for\n"
    "direct MASLD and liver-fat phenotypes carried\n"
    "23.34% of their posterior probability in protein-altering or canonical-splice\n"
    "variants, against 3.87% for the liver enzymes ALT, AST and GGT, a 6.026671-fold\n"
    "ratio (95% physical-cluster bootstrap interval 2.574593 to 13.510165;\n"
    "difference 0.194695, 0.069229 to 0.339940). Across 57 leave-one-study-out,\n"
    "leave-one-cluster-out, window-size and LD-block comparisons the ratio ranged\n"
    "from 4.368132 to 8.275691 and every\n"
    "difference interval stayed above zero.",
    "Pooled across\n"
    "ancestries, credible sets for direct MASLD and liver-fat phenotypes carried\n"
    "23.34% of their posterior probability in protein-altering or canonical-splice\n"
    "variants, against 3.87% for the liver enzymes ALT, AST and GGT. This pooled\n"
    "contrast mixes trait class with GWAS ancestry and fine-mapping coverage: 106 of\n"
    "207 direct-trait but 82 of 6,208 enzyme credible sets are European, the\n"
    "European enzyme sets cover 8 of the 767 UKBB ALT, AST and GGT lead windows,\n"
    "and the MVP European enzyme GWAS were not fine-mapped. Within European studies\n"
    "the fractions were 29.88% and 6.92% (difference 0.2296, 95% interval 0.0550 to\n"
    "0.4236); within AMR, the one ancestry with comparable fine-mapping coverage in\n"
    "both classes, they were 6.73% and 4.54% (difference 0.0219, −0.0489 to\n"
    "0.1065). No stratum passed all of its coverage, power and calibration checks,\n"
    "so no trait-class call is made, and the pooled direct-trait excess comes mainly\n"
    "from three PDFF GWAS. A trait-class comparison needs fine-mapping of all UKBB\n"
    "enzyme lead windows and of the MVP European ALT, AST and NAFLD GWAS.", True)
add("R02-01", MS + "01_intro.md", "77-84",
    "Fine-mapping the trait-associated regions shows that credible sets for direct\n"
    "MASLD and liver-fat phenotypes carry 23.34% of their posterior probability in\n"
    "protein-altering or canonical-splice variants, against 3.87% for the liver\n"
    "enzymes ALT, AST and GGT (ratio 6.026671, 95% physical-cluster bootstrap interval\n"
    "2.574593 to 13.510165). Separately, the strongest evidence for 423 of 436\n"
    "colocalized genes comes from a liver-enzyme trait and for only 13 from a direct\n"
    "disease or liver-fat trait, so the genetic map that expression-based\n"
    "colocalization returns is largely a map of liver-enzyme regulation.",
    "Pooled across ancestries, credible sets for direct MASLD and liver-fat\n"
    "phenotypes carry 23.34% of their posterior probability in protein-altering or\n"
    "canonical-splice variants, against 3.87% for the liver enzymes ALT, AST and GGT.\n"
    "The pooled contrast mixes trait class with GWAS ancestry and fine-mapping\n"
    "coverage. Within European studies the fractions are 29.88% and 6.92%, within\n"
    "AMR 6.73% and 4.54%, and no stratum supports a calibrated trait-class call.\n"
    "Separately, the strongest evidence for 423 of 436 colocalized genes comes from a\n"
    "liver-enzyme trait and for 13 from a direct disease or liver-fat trait, and 22\n"
    "genes have any direct-trait colocalization. Far more enzyme than direct-trait\n"
    "gene-by-GWAS pairs were tested by multi-signal colocalization (12,121 against\n"
    "771).", True)
add("R02-01", MS + "fig2.md", "103",
    "Taking protein-altering and canonical-splice variants together, the two trait classes carry 23.3428% against "
    "3.8732% of their posterior probability, a 6.026671-fold ratio with a 95% physical-cluster bootstrap interval of "
    "2.574593 to 13.510165 and a difference of 0.194695 (0.069229 to 0.339940). The interval resamples 233 physical "
    "1-Mb clusters: 19 direct, 228 enzyme and 14 shared. Native credible-set member coordinates uniquely restore 715 "
    "sentinels whose decimal encoding had lost trailing zeros; historical September intervals from the faulty "
    "coordinates are superseded. The complete corrected robustness menu contains 57 comparisons: one reference, 13 "
    "direct-study omissions, 18 enzyme-study omissions, 19 direct-cluster omissions, five grouping windows between "
    "0.25 and 5 Mb, and one ancestry-matched LD-block comparison (ratio 6.084745). The reduction from 62 comparisons "
    "follows from restoring 19 rather than 24 direct clusters. Ratios range from 4.368132 to 8.275691, and every "
    "difference interval has a positive lower bound (minimum 0.022373). These correlated sensitivity comparisons are "
    "not independent replications. No permutation p-value is reported, and bootstrap tail diagnostics are not "
    "calibrated trait-class hypothesis tests. COLOC r2 was separately adopted on October 14.",
    "Taking protein-altering and canonical-splice variants together, the pooled trait classes carry 23.3428% against "
    "3.8732% of their posterior probability (difference 0.194695, 95% physical-cluster bootstrap interval 0.069229 to "
    "0.339940; 233 physical 1-Mb clusters: 19 direct, 228 enzyme and 14 shared). This pooled value is not a "
    "trait-class contrast. Of the direct credible sets, 106 of 207 are European and 15 of the 19 direct clusters "
    "contain a European direct set, whereas 6,126 of 6,208 enzyme credible sets are non-European, so trait class is "
    "confounded with GWAS ancestry, LD reference and fine-mapping coverage. Stratified by ancestry (Figure 2D, "
    "ancestry panel; 10,000 shared draws over 81 pooled clusters, seed 20260923), European studies give 29.88% "
    "against 6.92% (difference 0.2296, 95% interval 0.0550 to 0.4236), AMR gives 6.73% against 4.54% (0.0219, "
    "−0.0489 to 0.1065), and the prespecified fixed-weight European-plus-AMR estimate (weights 16/26 and 10/26) is "
    "0.1497 (0.0388 to 0.2718). All non-European studies pooled give 6.99% against 3.68% (0.0331, −0.0383 to "
    "0.1274), a descriptive value with no calibration control. None of these receives a trait-class call. The "
    "European estimate failed the coverage check (enzyme arm fine-mapped at 8 of 767 lead windows, direct arm at 47 "
    "of 47), the power check (minimum detectable difference 0.304 against the pooled 0.195) and its random-split "
    "calibration control; the AMR and European-plus-AMR estimates failed their random-split calibration controls, "
    "and the European-plus-AMR estimate also carries the European coverage failure. Within the direct class, the "
    "three European PDFF GWAS (37 credible sets, 12 clusters) carry 40.14% against 3.87% for all enzymes "
    "(difference 0.3626, 0.1835 to 0.5455), and NAFLD plus NASH (170 sets, 13 clusters) carry 12.44% (0.0857, "
    "−0.0164 to 0.2096). The direct excess is also concentrated: five of 19 direct clusters carry 76.0% of the "
    "summed direct cluster fractions, and the largest, one PDFF credible set at chr11:823,586, carries 22.5%; the "
    "PNPLA3 cluster, with 124 direct sets, has a mean of 0.117. In an ungated sensitivity analysis (6,467 credible "
    "sets, seed 20260817), matching on standardized lead-variant effect gives ratio 5.78 (95% interval 1.99 to "
    "48.05; 27 matched enzyme loci); direct loci come from smaller GWAS (median effective n 21,893 against 160,000) "
    "and carry larger standardized effects (median 0.0368 against 0.0183). Native credible-set member coordinates "
    "uniquely restore 715 sentinels whose decimal encoding had lost trailing zeros; historical September intervals "
    "from the faulty coordinates are superseded. The 57 adopted robustness comparisons (study and direct-cluster "
    "omissions, five grouping windows and one LD-block comparison) give ratios 4.368132 to 8.275691 but include no "
    "ancestry or trait stratification and are correlated sensitivities, not replications. No permutation p-value "
    "is reported. Comparing the trait classes within one ancestry requires SuSiE fine-mapping of all 767 UKBB ALT, "
    "AST and GGT lead windows and of the MVP European ALT, AST and NAFLD GWAS. COLOC r2 was separately adopted on "
    "October 14.", True)
add("R02-01", LG, "60-66",
    "Coding plus splice mass is 23.34% versus 3.87%, ratio 6.026671, with a 95%\n"
    "physical-cluster bootstrap interval of 2.57459–13.51016 (10,000 draws, seed\n"
    "20260817). The 233 pooled clusters include 19 direct, 228 enzyme and 14 shared;\n"
    "a sampled shared cluster carries both classes. The complete 57-comparison\n"
    "robustness menu gives ratios 4.36813–8.27569, with positive lower bounds for\n"
    "every difference interval. These are correlated sensitivities rather than\n"
    "independent replications. No calibrated trait-class p-value is reported.",
    "Pooled across ancestries, the proportion of fine-mapped posterior probability\n"
    "carried by protein-altering or canonical-splice variants is 23.34% for direct\n"
    "traits and 3.87% for enzyme traits (difference 0.1947, 95% physical-cluster\n"
    "bootstrap interval 0.0692 to 0.3399; 10,000 draws, seed 20260817; 233 pooled\n"
    "clusters: 19 direct, 228 enzyme, 14 shared; a sampled shared cluster carries\n"
    "both classes). Because 106 of 207 direct but 82 of 6,208 enzyme credible sets\n"
    "are European, this pooled value is not a trait-class contrast. **(D, ancestry\n"
    "panel)** Each row prints the direct and enzyme fractions and the\n"
    "direct-minus-enzyme difference in percentage points with its 95% interval:\n"
    "pooled 23.3% vs 3.9%, 19.5 (6.9 to 34.0); EUR+AMR with fixed weights 16/26 and\n"
    "10/26, 21.0% vs 6.0%, 15.0 (3.9 to 27.2); EUR 29.9% vs 6.9%, 23.0 (5.5 to\n"
    "42.4); AMR 6.7% vs 4.5%, 2.2 (−4.9 to 10.7); all non-EUR 7.0% vs 3.7%, 3.3\n"
    "(−3.8 to 12.7). Thick bars are 95% and thin bars 98.333% intervals (Bonferroni\n"
    "over the three prespecified members EUR+AMR, EUR and AMR; 10,000 shared draws\n"
    "over 81 pooled clusters, seed 20260923); the pooled and all-non-EUR rows are\n"
    "descriptive. No member receives a trait-class call: the EUR enzyme arm was\n"
    "fine-mapped at 8 of 767 lead windows, the EUR member also failed the power\n"
    "check, and all three members failed their random-split calibration controls.\n"
    "**(D, direct-subclass panel)** Three EUR PDFF GWAS (37 sets, 12 clusters)\n"
    "40.14% and NAFLD plus NASH (170 sets, 13 clusters) 12.44%, each against all\n"
    "enzyme traits (3.87%); 95% intervals from 10,000 physical-cluster draws, seed\n"
    "20260923. Five of 19 direct clusters carry 76.0% of the summed direct cluster\n"
    "fractions. The 57-comparison robustness menu (ratios 4.36813 to 8.27569) has no\n"
    "ancestry or trait stratification; these are correlated sensitivities rather\n"
    "than independent replications. No calibrated trait-class p-value is reported.")
add("R02-01", MS + "fig5_discussion.md", "105-121",
    "The fine-mapped architecture differs descriptively between trait classes.\n"
    "Member-supported coordinate restoration gives protein-altering plus\n"
    "canonical-splice posterior fractions of 23.34% for direct MASLD/liver-fat\n"
    "traits and 3.87% for enzyme traits, ratio 6.026671 (95% physical-cluster\n"
    "bootstrap interval 2.574593–13.510165). The 233 physical clusters include\n"
    "19 direct, 228 enzyme and 14 shared; sampled shared clusters carry both\n"
    "classes. The 57 correlated robustness comparisons cover study and direct-cluster\n"
    "omission, clustering windows and LD blocks, giving ratios 4.36813–8.27569 with\n"
    "positive lower bounds for every difference interval. They are sensitivities,\n"
    "not independent replications; no calibrated trait-class p-value is reported.\n"
    "The earlier 6.69 ratio, its interval and 62-refit menu are historical and\n"
    "withdrawn after coordinate correction. Separately, the strongest native\n"
    "colocalization evidence for 423 of 436 genes comes from liver-enzyme traits\n"
    "and for 13 from direct disease/liver-fat traits. Some protein-altering variants\n"
    "may act through protein function without changing steady-state RNA, while\n"
    "others change RNA abundance through nonsense-mediated decay or splicing\n"
    "(Rivas et al., 2015).",
    "Pooled across ancestries, protein-altering plus canonical-splice variants\n"
    "carry 23.34% of fine-mapped posterior probability for direct MASLD/liver-fat\n"
    "traits and 3.87% for enzyme traits. This pooled contrast is confounded with\n"
    "GWAS ancestry and fine-mapping coverage. Within European studies the fractions\n"
    "are 29.88% and 6.92%; within AMR they are 6.73% and 4.54%. No stratum supports\n"
    "a calibrated trait-class call: the European estimate failed the coverage check\n"
    "(the UKBB enzyme GWAS were fine-mapped at 8 of 767 lead windows) and the power\n"
    "check, and the European, AMR and combined estimates all failed their\n"
    "random-split calibration controls. The pooled direct-trait excess comes mainly\n"
    "from the three PDFF GWAS and five of 19 direct clusters. Comparing the trait\n"
    "classes within one ancestry needs fine-mapping of all UKBB enzyme lead windows\n"
    "and of the MVP European ALT, AST and NAFLD GWAS. The earlier 6.69 ratio, its\n"
    "interval and 62-refit menu are historical and withdrawn after coordinate\n"
    "correction. Separately, the strongest native colocalization evidence for 423\n"
    "of 436 genes comes from liver-enzyme traits and for 13 from direct\n"
    "disease/liver-fat traits; 22 genes have any direct-trait colocalization, and\n"
    "far more enzyme than direct-trait gene-by-GWAS pairs were tested (12,121\n"
    "against 771). Some protein-altering variants may act through protein function\n"
    "without changing steady-state RNA, while others change RNA abundance through\n"
    "nonsense-mediated decay or splicing (Rivas et al., 2015).", True)
add("R02-01", "docs/manuscript/METHODS.md", "1006-1007",
    "liver lineage-accessible chromatin were nonexclusive context annotations and\n"
    "were not summed as probabilities.",
    "liver lineage-accessible chromatin were nonexclusive context annotations and\n"
    "were not summed as probabilities. Context was assigned by GRCh38 position after\n"
    "1:1 liftOver of the 16,324 sites to be annotated, because lookup by variant key\n"
    "had failed for sites absent from the annotation table, almost all in\n"
    "non-European studies. Direct-minus-enzyme coding plus canonical-splice\n"
    "differences were estimated within EUR, within AMR and as a fixed-weight EUR+AMR\n"
    "combination (weights 16/26 and 10/26), with 10,000 shared physical-cluster\n"
    "bootstrap draws over the 81 pooled clusters holding EUR or AMR credible sets\n"
    "(seed 20260923) and 98.333% percentile intervals (Bonferroni over these three\n"
    "members). A member could carry a trait-class call only if it passed three\n"
    "checks. Coverage: in each contributing ancestry, the share of GWAS lead windows\n"
    "with a credible set differed by at most twofold between classes. Power: the\n"
    "minimum detectable difference, 3.236 times the bootstrap standard deviation,\n"
    "was at most the pooled difference (0.195). Calibration: in 200 random splits of\n"
    "East Asian enzyme clusters (East Asian and African for the combined member)\n"
    "into pseudo-classes with the same cluster counts, the 95% interval covered zero\n"
    "in at least 180 splits and the 98.333% interval in at least 190. A difference\n"
    "was called only if, in addition, its 98.333% lower bound was above zero; a\n"
    "member that passed all checks with an interval containing zero was reported as\n"
    "no difference detected. Noncoding-context differences within AMR used the same\n"
    "calibration rule with 99% intervals (Bonferroni over five contexts).")
add("R02-01", "docs/PAPER.md", "545-547",
    "- Direct MASLD and liver-fat credible sets have a larger descriptive\n"
    "  protein-altering plus canonical-splice posterior fraction than liver-enzyme\n"
    "  proxy credible sets: 23.3428% against 3.8732%, a 6.026671-fold ratio after",
    "- Pooled across ancestries, direct MASLD and liver-fat credible sets have a\n"
    "  larger descriptive protein-altering plus canonical-splice posterior fraction\n"
    "  than liver-enzyme proxy credible sets. The pooled contrast is confounded with\n"
    "  GWAS ancestry and fine-mapping coverage and is not a trait-class claim (H58,\n"
    "  2026-10-15: 106 of 207 direct but 82 of 6,208 enzyme credible sets are\n"
    "  European; EUR 29.88% vs 6.92%, AMR 6.73% vs 4.54%; no call in any stratum;\n"
    "  the pooled excess is carried mainly by three PDFF GWAS). Pooled values:\n"
    "  23.3428% against 3.8732%, a 6.026671-fold ratio after", True)
add("R02-01", "docs/PAPER.md", "558-560",
    "  4.36813 to 8.27569, and every difference interval has a positive lower bound\n"
    "  (minimum 0.02237). Noncoding-context, trait-directness and ascertainment\n"
    "  reruns also completed and were adopted for Figure 2 only. The historical",
    "  4.36813 to 8.27569, and every difference interval has a positive lower bound\n"
    "  (minimum 0.02237). None of the 57 comparisons stratifies by ancestry or trait\n"
    "  or has matched fine-mapping coverage, so they are correlated sensitivities of\n"
    "  the pooled value, not replications of a trait-class difference.\n"
    "  Noncoding-context, trait-directness and ascertainment reruns also completed\n"
    "  and were adopted for Figure 2 only; H58 (2026-10-15) withdrew the\n"
    "  noncoding-context unresolved fractions, which reflected failed lookup\n"
    "  (repaired unresolved 0% in both classes). The historical", True)
add("R02-01", "docs/STATUS.md", "44",
    "The descriptive direction persists, conditional on the studied GWAS portfolio.",
    "The descriptive direction persists for the pooled value, conditional on the\n"
    "studied GWAS portfolio. H58 (2026-10-15) found the pooled contrast confounded\n"
    "with ancestry and fine-mapping coverage: 106 of 207 direct but 82 of 6,208\n"
    "enzyme credible sets are European; EUR 29.88% vs 6.92% (difference 0.2296,\n"
    "95% interval 0.0550 to 0.4236), AMR 6.73% vs 4.54% (0.0219, -0.0489 to\n"
    "0.1065); no trait-class call in any stratum. The 36.49%/77.85%\n"
    "(cluster-weighted) and 52.11%/60.43% (summed over credible sets)\n"
    "unresolved-context values reflected failed lookup and are withdrawn (repaired\n"
    "unresolved 0%).", True)
add("R02-01", "docs/RESULTS.md", "5",
    "fraction of 23.34% for direct traits and 3.87% for enzyme traits, ratio 6.03.",
    "fraction of 23.34% for direct traits and 3.87% for enzyme traits, ratio 6.03,\n"
    "pooled across ancestries. H58 (2026-10-15) shows this pooled value mixes trait\n"
    "class with ancestry and fine-mapping coverage: EUR 29.88% vs 6.92% (difference\n"
    "0.2296 [0.0550, 0.4236]), AMR 6.73% vs 4.54% (0.0219 [-0.0489, 0.1065]); no\n"
    "stratum receives a trait-class call.", True)
add("R02-01", "docs/RESULTS.md", "2466-2468",
    "Direct-trait credible sets carry a larger share of posterior probability in\n"
    "coding variants than liver-enzyme credible sets, and the strongest colocalization\n"
    "evidence for most genes comes from a liver-enzyme trait.",
    "Pooled across ancestries, direct-trait credible sets carry a larger share of\n"
    "posterior probability in coding variants than liver-enzyme credible sets; H58\n"
    "(2026-10-15) found this pooled contrast confounded with ancestry and\n"
    "fine-mapping coverage, with no trait-class call in any stratum. The strongest\n"
    "colocalization evidence for most genes comes from a liver-enzyme trait, and far\n"
    "more enzyme than direct-trait gene-by-GWAS pairs were tested (12,121 against\n"
    "771).", True)
add("R02-01", "docs/README.md", "23-25",
    "**Figure 2 coordinate correction adopted, 2026-10-13.** The coding/splice\n"
    "posterior ratio is 6.026671 (95% cluster-bootstrap interval 2.574593–13.510165),\n"
    "with 233 physical clusters and 57 robustness comparisons.",
    "**Figure 2 coordinate correction adopted, 2026-10-13.** The coding/splice\n"
    "posterior ratio is 6.026671 (95% cluster-bootstrap interval 2.574593–13.510165),\n"
    "with 233 physical clusters and 57 robustness comparisons. This value is pooled\n"
    "across ancestries and confounded with GWAS ancestry and fine-mapping coverage\n"
    "(H58, 2026-10-15): within EUR 29.88% vs 6.92%, within AMR 6.73% vs 4.54%, no\n"
    "trait-class call.", True)
add("R02-01", "figures/main/fig2_genetics/README.md", "5-12",
    "Figure 2D's coordinate-corrected results and regenerated PDF were adopted by\n"
    "the user on 2026-10-13. Its coding/splice posterior ratio is 6.026671 (95%\n"
    "physical-cluster bootstrap interval 2.574593–13.510165), with 19 direct and\n"
    "228 enzyme clusters. Corrected fractions, captions and interval tables come\n"
    "from `figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/` and\n"
    "the rendering job `22219710`. The regenerated noncoding-context companion is\n"
    "at `panels/Fig2_companion_noncoding_annotation_gap.pdf`, outside the A-I sequence;\n"
    "it displays credible-set-weighted fractions.",
    "Figure 2D's coordinate-corrected results and regenerated PDF were adopted by\n"
    "the user on 2026-10-13. Its coding/splice posterior ratio is 6.026671 (95%\n"
    "physical-cluster bootstrap interval 2.574593–13.510165), with 19 direct and\n"
    "228 enzyme clusters. This ratio is pooled across ancestries and confounded with\n"
    "GWAS ancestry and fine-mapping coverage (H58, 2026-10-15): within EUR 29.88% vs\n"
    "6.92%, within AMR 6.73% vs 4.54%, no trait-class call. Corrected fractions,\n"
    "captions and interval tables come\n"
    "from `figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/` and\n"
    "the rendering job `22219710`. The regenerated noncoding-context companion is\n"
    "at `panels/Fig2_companion_noncoding_annotation_gap.pdf`, outside the A-I sequence;\n"
    "it displays fractions summed over all credible sets of a class. Its unresolved\n"
    "fractions, 52.11% and 60.43%, reflected failed context lookup and are withdrawn\n"
    "(H58 repaired value 0% in both classes).", True)
# conflict 5: title.md 22-24 (supersedes Claims E197, same old_text)
add("R02-01", MS + "title.md", "22-24",
    "“Architecture” refers to the distribution of fine-mapped posterior probability\n"
    "across variant consequence classes, which is measured and adopted. It does not\n"
    "claim a resolved causal variant at any locus.",
    "“Architecture” refers to how fine-mapped posterior probability at liver-trait\n"
    "loci divides between protein-altering, splice and noncoding variants. The\n"
    "direct-trait versus liver-enzyme contrast is descriptive, rests on 19\n"
    "direct-trait regions and is pooled across ancestries; it is confounded with\n"
    "GWAS ancestry and fine-mapping coverage and receives no trait-class call in any\n"
    "ancestry stratum. After position re-annotation, most noncoding probability in\n"
    "both classes lies outside the deposited promoter, open-chromatin and ABC\n"
    "annotations. It does not claim a resolved causal variant at any locus.", True)

# ------------------------------------------------------------------ R14-02 / R14-06 (RESULTS note outside the adopted table)
add("R14-02", "docs/RESULTS.md", "2536-2537",
    "The legacy trait-directness table also carries bootstrap tail diagnostics;\n"
    "these must not be reported as calibrated trait-class p-values.",
    "The legacy trait-directness table also carries bootstrap tail diagnostics;\n"
    "these must not be reported as calibrated trait-class p-values.\n"
    "\n"
    "**Ancestry and subclass stratification (H58 and review round 3, 2026-10-15;\n"
    "not part of the adopted table above).** Ancestry composition: 106 of 207\n"
    "direct and 82 of 6,208 enzyme credible sets are European, and 15 of 19 direct\n"
    "clusters contain a European direct set. H58 family C (10,000 shared draws over\n"
    "81 pooled clusters, seed 20260923; 98.333% intervals, Bonferroni over three):\n"
    "EUR+AMR 0.1497 [0.0388, 0.2718], EUR 0.2296 [0.0550, 0.4236], AMR 0.0219\n"
    "[-0.0489, 0.1065] (95% intervals shown); no call for any member. EUR failed\n"
    "coverage (enzyme 8 of 767 lead windows, direct 47 of 47), power (MDE 0.304 >\n"
    "0.195) and calibration (N1a: 181 and 188 of 200 at 95% and 98.333%); AMR failed\n"
    "calibration (N1b: 141 and 150 of 200); EUR+AMR carries the EUR coverage failure\n"
    "and failed calibration (N1c: 173 and 184 of 200). Direct subclasses against\n"
    "all enzyme clusters (3.87%): PDFF 40.14% (37 sets, 12 clusters; difference\n"
    "0.3626 [0.1835, 0.5455]) and NAFLD+NASH 12.44% (170 sets, 13 clusters; 0.0857\n"
    "[-0.0164, 0.2096]); 10,000 physical-cluster draws, seed 20260923. Direct-cluster\n"
    "concentration: the top 5 of 19 direct clusters carry 76.0% of the summed direct\n"
    "cluster fractions and one PDFF credible set at chr11:823,586 carries 22.5%;\n"
    "ratios after dropping the direct rows of the top 1 to 5 clusters are 4.93,\n"
    "4.15, 3.45, 2.69 and 1.96 (a biased diagnostic, not a test). The ungated\n"
    "ascertainment medians below (effective n 21,893 vs 160,000; standardized effect\n"
    "0.0368 vs 0.0183) cover 208 and 6,259 credible sets, including the 52 that fail\n"
    "the PIP-sum gate.")

# ------------------------------------------------------------------ R02-02 / R01-10
add("R02-02", MS + "fig2.md", "107",
    "The adopted coordinate-corrected cluster-weighted analysis also makes missing context visible: 36.49% of "
    "direct-trait and 77.85% of liver-enzyme noncoding posterior probability lacks resolved deposited regulatory "
    "context. Direct-trait ABC overlap is only 0.126%, and its bootstrap ratio interval reaches zero; sparse boundary "
    "behavior does not establish biological depletion. Consequently, the DNA result supports a broad noncoding "
    "architecture and an assay-coverage map, not a claim that every noncoding variant has a known regulatory "
    "mechanism.",
    "Context lookup had failed for sites absent from the annotation table, almost all in non-European studies (0.7% "
    "of gated credible-set member rows unresolved in EUR, 98.2% in AMR); re-annotating by GRCh38 position, with all "
    "16,324 sites lifting 1:1, leaves no noncoding posterior probability without a context call. Summed over all "
    "credible sets of each class, as in the companion panel, 13.84% of direct-trait and 14.67% of liver-enzyme "
    "noncoding posterior probability lies within 2 kb of a TSS, 9.26% and 8.44% in liver-lineage open chromatin, "
    "0.43% and 3.70% in a liver ABC enhancer, and 77.90% and 77.69% in sequence with none of these contexts. These "
    "pooled values mix ancestries. Within AMR, the one ancestry with comparable fine-mapping coverage in both "
    "classes, no difference was detected in the share of noncoding probability with no deposited context (direct "
    "minus enzyme −0.0957, 99% interval −0.344 to 0.126; Bonferroni over five contexts); promoter, ABC and "
    "accessible differences received no call because their random-split calibration controls failed. Most noncoding "
    "fine-mapped probability therefore lies outside the deposited promoter, open-chromatin and ABC annotations, and "
    "the DNA result does not assign a regulatory mechanism to most noncoding signals.", True)
add("R02-02", LG, "92-100",
    "**(Noncoding-context companion, not assembled in canonical A-I)** Deposited\n"
    "annotations of noncoding posterior mass. Displayed credible-set-weighted\n"
    "fractions are 52.11% versus 60.43% unresolved, 7.87% versus 4.46% within 2 kb\n"
    "of a TSS, 7.32% versus 4.15% open in a liver lineage, and 0.23% versus 1.57%\n"
    "overlapping a liver ABC enhancer (direct versus enzyme). Context categories\n"
    "overlap and must not be summed. The adopted region-balanced context intervals\n"
    "are separate source-table estimates; the sparse ABC overlap does not establish\n"
    "biological depletion. Unresolved context does not justify assigning the nearest\n"
    "gene, and regulatory-DNA context is distinct from transcript biotype.",
    "**(Noncoding-context companion, not assembled in canonical A-I)** Deposited\n"
    "annotations of noncoding posterior probability after re-annotating, by GRCh38\n"
    "position, sites whose context lookup had failed (16,324 sites, all lifted 1:1);\n"
    "no noncoding posterior remains unresolved. Each value is the context posterior\n"
    "summed over all credible sets of a class divided by the noncoding posterior\n"
    "summed over the same sets (direct versus enzyme): 13.84% versus 14.67% within\n"
    "2 kb of a TSS, 9.26% versus 8.44% open in a liver lineage, 0.43% versus 3.70%\n"
    "overlapping a liver ABC enhancer, and 77.90% versus 77.69% resolved with none of\n"
    "these contexts. Context categories overlap and must not be summed. The earlier\n"
    "displayed 52.11% and 60.43% unresolved fractions reflected failed lookup, not\n"
    "missing annotation, and are withdrawn. Pooled values mix ancestries; within AMR\n"
    "only the no-context difference passed its calibration control (no difference\n"
    "detected). Regulatory-DNA context is distinct from transcript biotype.")
add("R02-02", "docs/RESULTS.md", "1837-1840",
    "| Noncoding PIP mass within 2 kb of a GENCODE TSS | 7.87% | 4.46% |\n"
    "| Noncoding PIP mass in a liver ABC enhancer | 0.23% | 1.57% |\n"
    "| Noncoding PIP mass open in at least one liver lineage | 7.32% | 4.15% |\n"
    "| Noncoding PIP mass with unresolved deposited context | 52.11% | 60.43% |",
    "| Noncoding PIP mass within 2 kb of a GENCODE TSS (pre-repair; H58 repaired 13.84% / 14.67%) | 7.87% | 4.46% |\n"
    "| Noncoding PIP mass in a liver ABC enhancer (pre-repair; H58 repaired 0.43% / 3.70%) | 0.23% | 1.57% |\n"
    "| Noncoding PIP mass open in at least one liver lineage (pre-repair; H58 repaired 9.26% / 8.44%) | 7.32% | 4.15% |\n"
    "| Noncoding PIP mass with unresolved deposited context (withdrawn 2026-10-15: failed lookup; H58 repaired 0% / 0%) | 52.11% | 60.43% |")
add("R02-02", "docs/RESULTS.md", "1842-1845",
    "Regulatory-context percentages overlap and must not be summed. Credible-set\n"
    "instances are descriptive study-by-locus units, not independent biological\n"
    "replicates. The large unresolved fractions prevent a claim that every\n"
    "noncoding signal has a localized regulatory mechanism.",
    "Regulatory-context percentages overlap and must not be summed. Credible-set\n"
    "instances are descriptive study-by-locus units, not independent biological\n"
    "replicates. The unresolved fractions came from failed lookup of sites absent\n"
    "from the context table, almost all in non-European studies (H58, 2026-10-15).\n"
    "After re-annotation by GRCh38 position (16,324 sites, all lifted 1:1), no\n"
    "noncoding posterior is unresolved, and 77.90% (direct) and 77.69% (enzyme) of\n"
    "noncoding posterior lies outside the deposited promoter, open-chromatin and ABC\n"
    "annotations (summed over all credible sets of each class). The DNA result\n"
    "therefore does not assign a localized regulatory mechanism to most noncoding\n"
    "signals.", True)
add("R02-02", "docs/RESULTS.md", "2539-2545",
    "The five context differences (direct minus enzyme; 95% physical-cluster\n"
    "intervals) are promoter 0.138659 [0.018337, 0.292354], ABC enhancer -0.005811\n"
    "[-0.011888, -0.001089], accessible 0.059563 [-0.004648, 0.135635], resolved with\n"
    "no context 0.226098 [0.076323, 0.379006], and unresolved context -0.413643\n"
    "[-0.616094, -0.199384]. Direct ABC fraction is only 0.126% and its ratio\n"
    "interval reaches zero; sparse/boundary behavior prevents a biological\n"
    "depletion claim.",
    "The five adopted context differences (promoter 0.138659, ABC -0.005811,\n"
    "accessible 0.059563, resolved-no-context 0.226098, unresolved -0.413643) are\n"
    "withdrawn: H58 (2026-10-15) showed the unresolved mass came entirely from failed\n"
    "lookup, almost all of it in non-European studies (summed noncoding PIP 3,604.4\n"
    "from sites absent from the lookup table and 47.8 from other allele order, of\n"
    "which 0.27 came from European rows). After position re-annotation (16,324\n"
    "sites, 1:1 liftOver) the pooled cluster-weighted differences are promoter\n"
    "0.0995 [-0.028, 0.248], ABC -0.0338 [-0.052, -0.017], accessible -0.0312\n"
    "[-0.102, 0.046], resolved-no-context -0.0896 [-0.229, 0.038], and unresolved 0.\n"
    "Pooled values mix ancestries; within AMR (99% intervals, Bonferroni over five)\n"
    "only resolved-no-context passed calibration (-0.0957 [-0.344, 0.126], no\n"
    "difference detected); promoter, ABC and accessible received no call.", True)
add("R02-02", "docs/RESULTS.md", "2566-2568",
    "remains outside A-I and displays credible-set-weighted fractions: unresolved\n"
    "52.11% / 60.43%, promoter 7.87% / 4.46%, accessible 7.32% / 4.15%, ABC\n"
    "0.23% / 1.57% (direct / enzyme).",
    "remains outside A-I and displays fractions summed over all credible sets of a\n"
    "class: unresolved 52.11% / 60.43%, promoter 7.87% / 4.46%, accessible\n"
    "7.32% / 4.15%, ABC 0.23% / 1.57% (direct / enzyme). H58 (2026-10-15) withdrew\n"
    "the unresolved values (failed lookup). Repaired values, same estimator\n"
    "(candidate panel\n"
    "`figures/candidates/review-r3-fixes-20261015T112148Z/F2d/repair/Fig2_companion_noncoding_context_repair.pdf`):\n"
    "unresolved 0% / 0%, promoter 13.84% / 14.67%, accessible 9.26% / 8.44%, ABC\n"
    "0.43% / 3.70%, resolved with none of these contexts 77.90% / 77.69%.")
add("R02-02", "docs/STATUS.md", "75-76",
    "The companion retains instance-weighted unresolved fractions\n"
    "52.11% and 60.43% and remains outside the assembled A–I panel set.",
    "The installed companion still prints unresolved fractions of 52.11% and\n"
    "60.43%, computed by summing context posterior over all credible sets of a\n"
    "class; H58 (2026-10-15) withdrew them (failed lookup; repaired value 0%). A\n"
    "corrected candidate companion is\n"
    "`figures/candidates/review-r3-fixes-20261015T112148Z/F2d/repair/Fig2_companion_noncoding_context_repair.pdf`.\n"
    "It remains outside the assembled A–I panel set.")

# ------------------------------------------------------------------ R14-03
add("R14-03", MS + "abstract.md", "32-35",
    "The strongest colocalization evidence for\n"
    "423 of 436 genes came from a liver-enzyme trait and for only 13 from a direct\n"
    "disease or liver-fat trait, so expression-based mapping of liver traits returns\n"
    "predominantly a liver-enzyme regulatory map.",
    "The strongest colocalization evidence for\n"
    "423 of 436 genes came from a liver-enzyme trait and for 13 from a direct\n"
    "disease or liver-fat trait; 22 genes have any direct-trait colocalization, and\n"
    "far more enzyme than direct-trait gene-by-GWAS pairs were tested (12,121\n"
    "against 771).", True)
add("R14-03", MS + "fig2.md", "105",
    "Separately, the strongest colocalization evidence for 423 of 436 genes comes from a liver-enzyme trait (ALT 137, "
    "AST 158, GGT 128) and for only 13 from a direct disease or liver-fat trait (NAFLD 8, PDFF 4, NASH 1). Because "
    "the table records one row per gene, this describes which trait class supplies each gene's best evidence rather "
    "than excluding weaker evidence from the other class. Expression-based colocalization applied to liver traits "
    "therefore returns predominantly liver-enzyme evidence. Whether the coding difference between trait classes "
    "contributes to this pattern, or to the small overlap with the disease-state gene set, is a hypothesis that this "
    "comparison does not test.",
    "Separately, the strongest colocalization evidence for 423 of 436 genes comes from a liver-enzyme trait (ALT 137, "
    "AST 158, GGT 128) and for 13 from a direct disease or liver-fat trait (NAFLD 8, PDFF 4, NASH 1). Because the "
    "table records one row per gene, this describes which trait class supplies each gene's best evidence; 22 genes "
    "have a direct-trait SuSiE PP.H4 above 0.5 (12 direct only and 10 with both classes; Figure 1C). Far more enzyme "
    "than direct-trait gene-by-GWAS pairs were tested by multi-signal colocalization: 12,121 against 771, of which "
    "739 (6.1%) and 28 (3.6%) colocalized. These per-pair rates are descriptive; pairs within a locus are not "
    "independent and no test is made. The shared-posterior rerun selected only first-pass positives (816 "
    "gene-by-GWAS pairs), so the 171 direct and 3,844 enzyme credible-set pairs eligible in that rerun do not "
    "measure per-signal yield, and 124 of the 171 direct pairs come from MVP European NAFLD, which has no "
    "fine-mapped credible sets. Expression-based colocalization applied to liver traits therefore returns "
    "predominantly liver-enzyme evidence. Whether the pooled coding difference contributes to this pattern, or to "
    "the overlap with the disease-state gene set, is a hypothesis that this comparison does not test.", True)
add("R14-03", "docs/PAPER.md", "567-574",
    "  Separately, the colocalized map is predominantly a\n"
    "  liver-enzyme map: the strongest evidence for 423 of 436 genes comes from a\n"
    "  liver-enzyme trait and for 13 from a direct disease or liver-fat trait. Some\n"
    "  protein-altering variants may act through protein function without changing\n"
    "  steady-state RNA, while others change RNA through nonsense-mediated decay or\n"
    "  splicing (Rivas et al. 2015, PMID 25954003). Any role of the coding difference\n"
    "  in the small genetic/disease-state overlap is a hypothesis this comparison did\n"
    "  not test.",
    "  Separately, the colocalized map is predominantly a\n"
    "  liver-enzyme map: the strongest evidence for 423 of 436 genes comes from a\n"
    "  liver-enzyme trait and for 13 from a direct disease or liver-fat trait, and\n"
    "  22 genes have any direct-trait colocalization. Far more enzyme than\n"
    "  direct-trait gene-by-GWAS pairs were tested by multi-signal colocalization\n"
    "  (12,121 against 771; 6.1% and 3.6% colocalized, descriptive, pairs not\n"
    "  independent). The 171 direct and 3,844 enzyme credible-set pairs of the\n"
    "  shared-posterior rerun were selected on first-pass colocalization, so their\n"
    "  hit rates do not measure per-signal yield. Some protein-altering variants may\n"
    "  act through protein function without changing steady-state RNA, while others\n"
    "  change RNA through nonsense-mediated decay or splicing (Rivas et al. 2015,\n"
    "  PMID 25954003). Any role of the pooled coding difference in the\n"
    "  genetic/disease-state overlap is a hypothesis this comparison did not test.", True)

# ------------------------------------------------------------------ check
bad = []
by_file = {}
for fid, eds in E.items():
    for ed in eds:
        txt = (ROOT / ed["file"]).read_text()
        n = txt.count(ed["old_text"])
        if n != 1:
            bad.append((fid, ed["file"], ed["line"], n))
        else:
            s = txt.index(ed["old_text"])
            by_file.setdefault(ed["file"], []).append((s, s + len(ed["old_text"]), fid, ed["line"]))
        for ch in ("—",):
            if ch in ed["new_text"]:
                bad.append((fid, ed["file"], ed["line"], "em dash"))
        for tag in ("[Inference]", "[Unverified]", "[Estimate]", "[Speculation]"):
            if tag in ed["new_text"]:
                bad.append((fid, ed["file"], ed["line"], tag))
for f, spans in by_file.items():
    spans.sort()
    for a, b in zip(spans, spans[1:]):
        if a[1] > b[0]:
            bad.append(("overlap", f, a, b))
print("edits:", sum(len(v) for v in E.values()), {k: len(v) for k, v in E.items()})
print("problems:", bad)
(OUT / "final_text_edits.json").write_text(json.dumps(E, indent=1, ensure_ascii=False))
