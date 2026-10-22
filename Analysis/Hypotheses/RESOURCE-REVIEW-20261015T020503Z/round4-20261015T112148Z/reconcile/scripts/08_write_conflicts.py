"""Write the reconciled conflict list (ranked by effect on a scientific claim) to conflicts.json and conflicts.tsv."""
import json, csv
C = []
def add(kind, where, detail, resolution):
    C.append(dict(kind=kind, where=where, detail=detail, resolution=resolution))

add("edit_conflict",
    "docs/manuscript/draft/title.md:10-11; docs/PAPER.md:64-65; docs/README.md:101-102; docs/STATUS.md:236-238; title.md:14-16 and :26",
    "Two different titles are drafted for the decided row (H52 untestable / H54 supports, coordination clause removed). Claims R13-01 (E187, also E191/E193/E194 in PAPER, README, STATUS) gives 'describes the coding and noncoding architecture of liver-trait genetic signals and maps cross-sectional multicellular remodeling'. F4abcd R13-04 (E327) gives 'resolves the architecture of inherited risk and maps cross-sectional multicellular remodeling' and edits title.md only, so installing it leaves PAPER, README and STATUS with the Claims title. Claims E189 (title.md:26) still contains the placeholders '[H52 outcome] / [H54 outcome]'. title.md:14-16 still says the 2026-09-07 title 'converts two topics into two findings', but the new second clause 'maps ...' is a topic again.",
    "Install the Claims set (E187, E189, E191, E193, E194) and drop all six F4abcd R13-04 title rows (E323-E328). 'Resolves the architecture' is the wording the review document lists as unsupported, and H58 found the coding/splice contrast confounded with ancestry (C3 untestable). Fill E189 with 'H52 untestable / H54 supports'. Rewrite title.md:14-16 to give the 2026-10-15 reason (clause removed under the joint rule; 'resolves' replaced because of H58). title.md:3-4 requires an explicit paper-level decision, so the user must approve before install.")

add("edit_conflict",
    "docs/manuscript/05_figure_legends.md:508-511; docs/RESULTS.md:3175; docs/STATUS.md:831; docs/manuscript/draft/fig7.md:37-40; docs/manuscript/METHODS.md:835-840",
    "Fig 7 severity transfer has two incompatible primary sets on the same legend sentence. F7a R08-12 (E147) prints 585 participant records and 0.68 to 0.79 of the tie ceiling (published record sets). F7b R08-02 (E237) prints 453 participants with MASLD (109 + 263 GSE213621 genetic identities + 81 GSE276114 MASLD) and 0.70 to 0.87. METHODS (E240) makes the identity-level and MASLD-only values primary, while RESULTS 3175 (E238) and STATUS 831 (E239) still lead with the record-level 0.677/0.774/0.794. RESULTS 3175, STATUS 831 and fig7.md print 0.681 for chronic viral hepatitis; CVH alone is 0.680 (0.680496, n=82); 0.681 is CVH plus ARLD.",
    "Use the F7b set as primary everywhere: genetic identities and MASLD participants are the inferential units, and CVH is not MASLD. Drop E147 (keep F7a E146 for panel D). Lead RESULTS 3175 and STATUS 831 with 0.709 and 0.866 and give the record-level values as the published-record sensitivity. Print 0.680 for CVH in RESULTS 3175, STATUS 831 and fig7.md. Panel 7D must show the same values as the legend; if the installed 7D shows record-level values, either re-render or have the legend say so. The new values come from figures/candidates/review-r3-fixes-20261015T112148Z/F7b, not a named release.")

add("edit_conflict",
    "Fig 4E arm: 05_figure_legends.md:248-258; draft/fig4.md:5-6, 113-121, 193-195; PAPER.md:444-450; draft/01_intro.md:97-99; RESULTS.md:1635-1636, 2146-2147, 1727-1729; STATUS.md:6, 586, 609-610, 888; docs/README.md:7; ROADMAP.md:299; figures/main/fig4_singlecell_programs/README.md:8-9, 37",
    "F4efgh R13-03 rewrites Figure 4E as the 364-unit three-cohort arm (0 fibrosis, 0 NAS calls; smallest q 0.079/0.159). F4abcd R10-01 (fig4.md:193-195, RESULTS 1635-1636 and 2146-2147) calls the 457-record arm 'the Figure 4E arm' (20 of 114 vs 26.9 expected, P=0.871). F4abcd R04-10 (E316, intro 97-99) states the ECM/IGFBP7 disease association as 'unadjusted bulk fibrosis (BH q=0.045) ... composition adjustment (q=0.230)', which are the 457-record values (RESULTS 1722); in the 364-unit arm they are q=0.214 and q=0.421. README, STATUS 6/586/888, ROADMAP 299 and the Fig 4 README still name the 457-record arm. RESULTS 1727-1729 records that the user chose the four-cohort arm on 2026-10-13.",
    "Pick one arm and make all listed places name it. Recommendation: the 364-unit arm (R13-03). It stays inside the identity-resolved 774-unit substrate, and all 20 four-cohort calls need 93 non-identity-resolved GSE174478 records plus the untruncated Knapp-Hartung variance. Then fig4.md:193-195 reads '1 against 1.24 (p=0.473)' and intro 97-99 reads 'unadjusted bulk fibrosis association is q=0.214 in the 364 units (q=0.045 only with 93 GSE174478 source records)'. This reverses a recorded user choice, so it needs the user's decision before install.")

add("abstract_title_intro",
    "docs/manuscript/draft/abstract.md:35-38; draft/01_intro.md:62-64 and 85-89; draft/fig5_discussion.md:124, 128, 165; docs/PAPER.md:120",
    "The H57 fix changes Fig 1B to 26 shared against 21.2 expected at matched expression (O/E 1.23, MH OR 1.29, P 0.28), and F1's own text says two sets of these sizes are mostly disjoint whether or not the biology differs. Untouched sentences still treat the overlap as smaller than expected and offer explanations for it: abstract 37-38 'whether this contributes to the small overlap is a hypothesis these data do not test'; intro 62 'Connecting the two maps exposes a limit...' followed by F1 E112 'share few genes'; intro 87-88 'the sparse overlap'; discussion 124, 128, 165 'the sparse overlap'; PAPER 120 'Their sparse overlap'. No edit covers them (F3cdg E022 starts after the abstract sentence).",
    "Abstract: delete 'Some protein-altering variants may act ... these data do not test.' Intro 62-64: replace the lead with 'The two maps can be compared gene by gene:' and drop 'share few genes'. Intro 87-89 and discussion 124, 128, 165: replace 'the sparse overlap' with 'which genes appear in both sets' and keep the 'not tested' clause. PAPER 120: 'Their overlap (26 genes, 21.2 expected at matched expression) does not establish...'.")

add("uncovered_stale_number",
    "docs/manuscript/draft/fig2.md:86-87 (F2d R13-08, E235)",
    "The new text says the colocalized and disease-state sets overlap 'at about the rate expected for sets of these sizes (26 shared against 33.46 expected)'. 33.46 is the unstratified expectation in the withdrawn universe that counted 10,928 unevaluated genes as negatives (H57 D5). F1 replaces it with 21.2 at matched expression among 3,769 evaluated genes everywhere else.",
    "Use '(26 shared against 21.2 expected for genes of the same expression level among the 3,769 genes evaluated by multi-signal colocalization)'. The 'about the rate expected' wording at fig2.md:116-117 (E236) is consistent with 21.2 and needs no change.")

add("abstract_title_intro",
    "docs/manuscript/draft/title.md:22-24 (Claims R13-02, E197)",
    "New text says 'Most noncoding probability (74.9% direct, 93.7% enzyme) has no deposited regulatory context.' These are pre-repair values (resolved-no-context plus unresolved: 0.3842+0.3649 and 0.1581+0.7785). H58 showed the unresolved mass was failed lookup; after repair unresolved is 0 and resolved-no-context is 68.5% direct and 77.4% enzyme (cluster-weighted; credible-set-weighted 77.9% and 77.7%). F2d R02-02 withdraws the unresolved values elsewhere.",
    "Print '(68.5% direct, 77.4% enzyme after the H58 context repair)', or drop the parenthesis.")

add("abstract_title_intro",
    "docs/manuscript/draft/abstract.md:32-35 (E231); draft/01_intro.md:77-84 (E214); docs/PAPER.md:567-569 (E232); draft/fig2.md:105 (E333); draft/fig5_discussion.md:105-121 (E217)",
    "F2d R14-03 adds 'The colocalized map is dominated by liver-enzyme evidence because far more enzyme signals were tested' and, in the intro, 'the enzyme majority reflects the 22-fold larger number of enzyme signals tested', based on 17.5% vs 21.1% per eligible pair. check:F2d found the eligible pairs come only from the 816 gene-trait rows already positive in the first pass. Per tested row the rates are 3.6% direct vs 6.1% enzyme, per tested credible-set pair 1.5% vs 1.9%, and enzyme signals are 15.7-fold (rows) or 20.8-fold (pairs) more numerous, not 22-fold.",
    "Keep the counts (423 of 436; 13; 22 genes with any direct-trait colocalization; 171 and 3,844 eligible pairs) and drop the causal clause. Abstract: 'Far more enzyme than direct-trait signals were tested.' Intro: delete 'so the enzyme majority reflects the 22-fold larger number of enzyme signals tested' and the 17.5%/21.1% rates, or report the per-tested-row rates.")

add("abstract_title_intro",
    "docs/manuscript/draft/abstract.md:19-23; draft/01_intro.md:68-74; draft/fig5_discussion.md:102, 135-145, 160-165; docs/PAPER.md:589-592",
    "The hepatocyte-specificity comparison (281 expression-matched pairs from 399 genetic-only and 1,129 state-only genes) is stated in the abstract and intro as a failed prespecified prediction with an opposite-direction post hoc difference. H57 D6/D7 show that 935 of the 1,129 state-only genes have no multi-signal COLOC result (mean AveExpr 0.232 vs 1.319 for the 194 tested). The 'state-only' pool is therefore mostly genes that could not be tested genetically, not genes without a genetic signal. F1 E125 adds this caveat in RESULTS only; the comparison was not rerun.",
    "Rerun the 281-pair comparison with the 194 tested state-only genes (not run). Until then add to abstract, intro and discussion: 'most state-only genes (935 of 1,129) could not be tested by multi-signal colocalization.'")

add("edit_conflict",
    "docs/manuscript/draft/abstract.md:24-32 and draft/01_intro.md:77-84",
    "Claims R13-02 (E195, E196) and F2d R02-01 (E213, E214) rewrite the same sentences. Claims keeps the 6.026671-fold ratio and adds credible-set-weighted 15.81% vs 3.26%; F2d drops the ratio and adds the H58 ancestry strata (EUR 29.88% vs 6.92%, AMR 6.73% vs 4.54%, no call). Both use 'pooled' with different meanings (across ancestries vs over credible sets).",
    "Install F2d's text and drop E195/E196. If the credible-set weighting is wanted, add one clause to F2d's sentence: '23.34% (cluster-weighted; 15.81% over credible sets) against 3.87% (3.26%)'. Do not print the ratio in the abstract.")

add("edit_conflict",
    "docs/manuscript/05_figure_legends.md:432; draft/abstract.md:60-63 (E208); draft/01_intro.md:111-115 (E209)",
    "Claims R07-01 (E211) and F6 R15-05 (E282) replace the same legend sentence differently. E211 says protein, chromatin and spatial context are attached at program level for Stromal ECM and Ductular injury; check:Claims found the only protein row (PXD051911) is untestable (insufficient gene coverage) for both programs. The abstract adds 'protein, chromatin and spatial context enter at program level only' and 'links genes to the 117 cell programs by membership', but only 3,422 of 28,233 identities are program members.",
    "One legend sentence: 'Chromatin and spatial context are attached for two programs (Stromal ECM and Ductular injury); liver protein was untestable for both. No gene-level protein, chromatin or spatial call is joined to the r8 Catalog, and no record is a regulatory DNA element object.' Abstract: 'links its 3,422 program-member genes to the 117 cell programs ... chromatin and spatial context are attached for two programs.' Intro 111-115 the same.")

add("edit_conflict",
    "docs/PAPER.md:312-318 and docs/manuscript/METHODS.md:1079-1085",
    "Methods R11-06 (E084, E085) rewrites the paragraph and keeps 'zero strict tested_negative calls'; F6 R07-02 (E259, E264) replaces that clause with 'defines no negative-call rule, so tested_negative is not assessed'. After R11-06, F6's old text no longer exists. R11-06 also adds 'the liver proteome meets the coverage rule for 26 of 117', taken from an unadopted module candidate (RESULTS 4476), not the spatial observability rule, and a reader can confuse it with Fig 5E's 26 of 94 protein programs. 'Which is not a property of the assay' (GeoMx) is [Unverified]: no gene axis was loaded.",
    "Merge into R11-06's text with F6's 'defines no negative-call rule, so tested_negative is not assessed'. Delete the proteome '26 of 117' clause. For GeoMx write 'no gene axis loaded; direction and normalization unresolved'.")

add("edit_conflict",
    "docs/manuscript/draft/fig5.md:41-50 (E066); docs/PAPER.md:333-346 (E069); docs/RESULTS.md:2994-3017 (E078) and 3003-3004 (E283)",
    "H90 is untestable (G3: 14 genes, weight 0.636 < 0.80), so only the 'H90 fail or untestable' rows E066/E069 apply; E067/E070 must not be installed. E066 says 'neither was selected by a stage-association rule', but both were primary-selected in the 2026-08-07 registry (q 0.0012, 0.0007) and Ductular injury still is in the 2026-10-14 registry (q 0.0376). E069 says the programs 'differ' and the nulls 'separated them at the program level'; E078 says the calibration 'separates the programs'. The abstract (F5ef G3-04, E072) says the difference between programs does not reach P<0.05, and no prespecified between-program test exists. E078 replaces 2994-3017 and removes the span that F6 R15-09 (E283) appends to, which records that the r8 Catalog program-context row for Vu Stromal ECM reads 'supported' and should read 'source_dependent'.",
    "Install E066 and E069 only. Replace 'neither was selected by a stage-association rule' with 'fixed in advance from the superseded August stage fit; neither passes the corrected 117-program rule'. In PAPER and RESULTS write 'Stromal ECM exceeded its matched-gene null in both sources and Ductular injury did not; no prespecified test compared the two programs.' Append E283's Catalog sentence to the Vu paragraph of E078.")

add("uncovered_stale_number",
    "figures/main/fig4_singlecell_programs/README.md:35",
    "The 4C row says 'Ductular injury is the only selected primary call'. Under H59 (adopted text everywhere else via F4abcd R04-01) 0 of 117 programs are selected; the 62-donor Ductular injury call (q 0.0376) is a sensitivity. No edit covers this line.",
    "Replace with 'H59 primary fit (60 donors, per-dataset variance, permutation BH over 117): no program selected; the 62-donor Ductular injury call is a sensitivity only.'")

add("uncovered_stale_number",
    "figures/main/fig4_singlecell_programs/README.md:10; docs/ROADMAP.md:282",
    "Fig 4H is moved by F4efgh R05-02 from rho 0.7265 and 37/43 (435 units including 68 GSE213621 controls at stage 0) to 0.631 and 32/43 on 367 disease units. The Fig 4 README still says 'Figure 4H prints rho=0.73 and 37/43'; ROADMAP 282 says 'stage-complete inference 435 units ... rho 0.7265 and 37/43 same direction'. All other occurrences are edited.",
    "README 10: 'Figure 4H prints rho=0.63 and 32/43 same direction on 367 disease units, with no outline.' ROADMAP 282: same values, and say the 435-unit maxT/competitive counts include 68 controls.")

add("uncovered_stale_number",
    "docs/manuscript/05_figure_legends.md:203; docs/RESULTS.md:1279-1284",
    "Supplementary S3 (I-L) legend says 'stage-complete anchor models retain 435 genetic-identity clusters' with no mention that 68 GSE213621 controls without histology are placed at stage 0; the main-figure legends now state this. RESULTS 1279-1284 (S3 stage-vs-continuum layer table) prints effect rho on that roster: Hotspot 0.823, Hallmark 0.745, KEGG 0.669, Reactome 0.687, GO BP 0.740, GO MF 0.683. F4efgh computed disease-only values 0.629, 0.662, 0.675, 0.612, 0.671, 0.523 but drafted no edit. Fig 4G (435 incl. controls, not re-rendered) and 4H (367) are now on different rosters within one figure.",
    "Legend 203: add '(including 68 GSE213621 controls with no recorded histology, placed at stage 0)'. RESULTS 1279-1284: add a disease-only rho column with the six values and note that the stage-supported counts were not recomputed. State the 4G/4H roster difference in the Fig 4 legend, or re-render 4G on 367 units (not run).")

add("uncovered_stale_number",
    "figures/main/fig1_atlas_overview/README.md:8-10",
    "Fig 1B is described as 'the adopted 774-unit stable-Ensembl interface (14,704 jointly testable genes, 426 genetic, 1,155 canonical disease-state and 26 shared)'. F1 R01-01 changes Fig 1B to the 3,769 genes evaluated by multi-signal COLOC: 426 genetic, 220 disease-state, 26 shared against 21.2 expected. No edit covers the README.",
    "Replace with '3,769 genes evaluated by multi-signal COLOC (426 genetic, 220 canonical disease-state, 26 shared against 21.2 expected at matched expression; H57)'.")

add("uncovered_stale_number",
    "docs/manuscript/draft/fig2.md:103 (F2d R02-01, E215)",
    "New text says 'Direct credible sets are mostly European (16 of 19 direct clusters contain a European set)'. check:F2d recomputed 15 of 19 (check/F2d/out/02_coding.log:40); 16 counts clusters formed within EUR only, two of which fall in one pooled cluster. At the credible-set level 106 of 207 direct sets (51%) are European, so the abstract's 'most direct credible sets are European' (E213) rests on a bare majority.",
    "Print '15 of 19 direct clusters contain a European set; 106 of 207 direct credible sets are European'. Keep 'most' in the abstract only with the 51% figure, or write 'about half'.")

add("uncovered_stale_number",
    "docs/RESULTS.md:1989-1990",
    "After F1 E118 replaces the overlap table above it with the H57 U_eval values, the next sentence reads 'SuSiE-or-ABF on the same basis: 959 genetic, 57 overlap, 94.1% non-overlap, OR 0.73, nominal p 0.022'. Those values are on the withdrawn 14,720 universe that counts unevaluated genes as negatives, so 'same basis' is now false and the OR reads as depletion. H57 D8 (ABF map, U_full, 10 expression bins): 49 shared against 43.6 expected, MH OR 1.15 (0.85 to 1.57), P 0.41.",
    "Either label the sentence as the withdrawn 14,720-universe value, or replace it with the H57 D8 values marked as a diagnostic.")

add("uncovered_stale_number",
    "docs/PAPER.md:552-558; docs/README.md:23-26; figures/main/fig2_genetics/README.md:5-12",
    "F2d R02-01 relabels the coding/splice contrast as a pooled, ancestry-confounded descriptive value (H58). PAPER 545-547 gets the qualifier, but 552-558 keep 'interval of 2.57459-13.51016 for the ratio' and 'every difference interval has a positive lower bound' with no statement that none of the 57 comparisons stratifies by ancestry or trait (the legend and fig2.md now say so). docs/README.md:24 and the Fig 2 README print 'coding/splice posterior ratio is 6.026671 (95% interval 2.574593-13.510165)' with no H58 qualifier. The Fig 2 README also describes the installed noncoding companion, which prints the withdrawn unresolved fractions 52.11%/60.43%.",
    "PAPER 558: add 'none of the 57 stratifies by ancestry or trait, so they are correlated sensitivities of the pooled value, not replications of a trait-class difference.' Both READMEs: add 'pooled across ancestries; confounded with GWAS ancestry and fine-mapping coverage (H58); within EUR 29.88% vs 6.92%, within AMR 6.73% vs 4.54%, no trait-class call.' Fig 2 README: say the installed companion's 52.11%/60.43% are withdrawn (repaired 0%/0%).")

add("uncovered_stale_number",
    "docs/RESULTS.md:557, 1109, 1120, 1129; docs/STATUS.md:539, 817",
    "F3ef G2-01/G2-04 add to the legend and fig3.md that the F1-vs-F0 count (1,338) rests on GSE130970 (14 of 1,338 remain without it) and that GNMT's F0-to-F1 decrease (q=0.00913) does not hold without GSE130970 (q=0.207) or with joint library-complexity covariates; fig3.md now says 'the earlier steps do not hold'. RESULTS 557 (GNMT hero-gene row) still says the F0-to-F1 step 'stays supported', and RESULTS 1109/1120/1129 and STATUS 539/817 print 1,338 and q=0.00913 without the caveat. F3ef flagged these and drafted no edit.",
    "RESULTS 557: replace 'which stays supported on disease-only F0 ... ' with 'supported on disease-only F0 (q=0.0426) but not without GSE130970 (q=0.207) or with joint library-complexity covariates'. Add 'rests mainly on GSE130970 (14 of 1,338 without it)' at RESULTS 1109/1120 and STATUS 539/817.")

add("edit_conflict",
    "docs/manuscript/METHODS.md:1121-1122; docs/manuscript/05_figure_legends.md:403-409; figures/main/fig5_molecular_context/README.md:27-28",
    "F5bcd R11-02 (E159) rewrites the start of METHODS 1121, which removes the anchor of F5ef R11-03 (E065). E065's appended text says 'The count uses the same 767 COLOC r2 pairs as Figure 5D', while F5bcd changes Fig 5D to gene units (22 direct, 424 enzyme genes from 767 pairs). The installed 5D PDF and the Fig 5 README (27-28) still print pair counts (28 direct, 739 enzyme), so the R06-08/R11-02 legend text describes a candidate panel, not the installed one.",
    "Apply E159, then append E065's sentences after 'not a negative functional test.' and change the last one to 'uses the same 767 COLOC r2 pairs that Figure 5D summarizes per gene (22 direct, 424 enzyme)'. Install the gene-level 5D legend only together with the candidate 5D panel, and update the Fig 5 README in the same step.")

add("edit_conflict",
    "docs/manuscript/draft/fig5.md:35-39",
    "F5bcd R06-08 (E153) rewrites 35-39; F5ef R01-03 (E061) inserts the Fig 5E paragraph after line 39 using the old last line as anchor. The anchor survives in E153's text but with a different line break, so a literal string match fails. E061 also names 'amino-acid metabolism programs' (the robust row is 'Amino-acid degradation (TAT)', not 'Amino-acid metabolism (GOT1)'), says RORA programs carry a direct MASLD trait (the RORA direct pair is partial, 0.37, in T/NK) and says protein programs 'track a histologic-burden axis' (cross-sectional association).",
    "Apply E153, then append E061's paragraph after its last sentence. In E061: 'amino-acid degradation (TAT)'; 'are associated with histologic burden across participants'; describe RORA as lineage-open through a direct-trait pair in T/NK at 0.37 posterior.")

add("edit_conflict",
    "docs/manuscript/05_figure_legends.md:147-152; docs/manuscript/draft/fig3.md:124-127",
    "F3cdg R16-02 (E019, E020) and F3ef G2-03 (E175, E176) overlap: G2-03's old text begins with the last words of R16-02's span ('COLOC r2 input.' / 'contrast.'). Applying R16-02 first removes G2-03's anchor. The numbers do not contradict (54 genes; 15 under a split BH budget; 26 of 59 marks in 22 genes hold under both checks; 29 genes selected only through F0-to-F1). R16-02's 'bounds the FDR of the union at 5%' holds over all discoveries in the four 24,115-gene families under BH assumptions, not for the 436-gene subset.",
    "Install one merged legend sentence and one merged fig3.md sentence (R16-02 text followed by G2-03 text). Write 'which bounds the FDR over all four 24,115-gene families at 5%'.")

add("edit_conflict",
    "docs/RESULTS.md:1767-1768",
    "F5ef R06-02 (E068) writes 'Both prespecified programs (Stromal ECM and Ductular injury) passed this'; F4abcd R04-01 (E293) writes 'Both focal programs passed this lncRNA-removal sensitivity'. Same old text, different wording; both remove 'primary-selected'.",
    "Install E068 (names the programs) and drop E293.")

add("edit_conflict",
    "docs/manuscript/draft/fig1.md:44-45",
    "F4efgh R13-11 (E059): 'Figure 4 resolves donor-level single-cell programs and their bulk-tissue context; transcription-factor activity is supplementary.' F4abcd R13-04 (E321): 'Figure 4 reports donor-level program scores in five cell lineages and regulatory activity.' Same old text.",
    "One sentence: 'Figure 4 reports donor-level program scores in five cell lineages and their bulk-tissue context; transcription-factor activity is supplementary.'")

add("edit_conflict",
    "docs/RESULTS.md:2051-2053, 2221-2222; docs/STATUS.md:998-1002",
    "Methods R16-07 (E101, E103, E104) and F1 R01-01 (E119, E120, E121) rewrite the same overlap sentences. F1's text carries the H57 values (26 vs 21.21, MH OR 1.295); R16-07 predates them.",
    "Install F1 E119-E121 and drop E101, E103, E104. Keep R16-07's other edits (for example E102).")

add("abstract_title_intro",
    "docs/manuscript/draft/abstract.md:16-18 (F1 E111)",
    "Abstract prints 'Mantel-Haenszel odds ratio 1.30'. H57 results.json MH_OR is 1.2947, which rounds to 1.29; 1.30 comes from rounding the printed 1.295 a second time. PAPER and RESULTS print 1.295.",
    "Print 1.29 in the abstract (and in any other two-decimal copy).")

add("abstract_title_intro",
    "docs/manuscript/draft/abstract.md:65-67",
    "The last sentence, 'This Resource provides a reusable map from inherited variation to multicellular MASLD biology', is untouched. After the fixes, the colocalized and disease-state sets share no more genes than expected (H57), no program passes the 117-program stage rule (H59), and the coordination clause is removed from the title. The sentence states a link the corrected results do not show.",
    "Replace with 'This Resource records inherited-variation and disease-state evidence side by side for each gene, with the experiment needed to connect them.'")

add("uncovered_stale_number",
    "docs/STATUS.md:579, 834; docs/RESULTS.md:2832",
    "F6 relabels the 26-gene Catalog class as 'both supported (direction not compared)'. E273's old text 'genetically anchored 400, concordant 26, established-state associated 1,129,' occurs at STATUS 463 and 834, and a string install changes only 463. STATUS 579 and RESULTS 2832 also still say 'concordant 26'.",
    "Apply E273 at both STATUS 463 and 834, and change STATUS 579 and RESULTS 2832 the same way.")

add("uncovered_stale_number",
    "docs/manuscript/draft/fig5_discussion.md:266-269",
    "'the severity instrument at 0.6631 against 0.6094 for the in-cohort ridge' is untouched. F7b adds elsewhere that GSE268273 also chose the shipped expression transform (highest of four eligible arms, 0.6506 to 0.6631) and that the paired difference (+0.054) has a 95% interval that includes zero.",
    "Append '(GSE268273 also chose the expression transform; the paired difference, +0.054, has a 95% interval that includes zero)'.")

add("edit_conflict",
    "F3cdg R03-02 new text in 05_figure_legends.md, draft/fig3.md, RESULTS.md, PAPER.md, METHODS.md (E007-E011)",
    "Each new text contains the literal tag '[Inference]' before 'Unequal eligibility limits representativeness'. Installed as drafted, the tag would print in the figure legend and Methods.",
    "Remove the tag and write 'Unequal eligibility may limit representativeness.'")

add("edit_conflict",
    "docs/manuscript/05_figure_legends.md:43-45 and :60",
    "F2abce R02-05 (E128) chains on R02-04's new 2B text (E126) and cannot be installed against the current file. F2d R13-07 (E334) is a placeholder note ('replaced by the R02-01 legend edit above') whose old text sits inside F2d R02-01's span (E216).",
    "Install E126 then E128. Drop E334.")

json.dump(C, open('conflicts.json', 'w'), indent=1, ensure_ascii=False)
with open('conflicts.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=['kind', 'where', 'detail', 'resolution'], delimiter='\t')
    w.writeheader(); w.writerows(C)
print(len(C))
import re
for c in C:
    for k, v in c.items():
        if '—' in v: print('EM DASH in', c['where'][:40], k)
