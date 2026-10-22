# Current paper and codebase status

**Release 2026-10-21-r9 (built 2026-10-21; supersedes 2026-10-20-r8, r7, r6, r5, r4, r3 and 2026-10-17-r2).** Its manifest is
`RNA-seq/results/manuscript_release/2026-10-21-r9/release_manifest.tsv`. r9 applies the triage of two external reviews of r8 (`Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/`; user decisions D1-D4 in DECISIONS.md): Figure 2D restricted to genome-wide leads, Figure 3C with the GSE126848 filtered refit, the shared-posterior check of carried negatives reported as a sensitivity, spatial specificity wording corrected against a 114-program family, identity-screen Methods, and about 270 text corrections; numbers in RESULTS.md, section 'Review triage 2026-10-21'. Fig 7 items are handed off in FIG7_HANDOFF.md and not applied. r8 replaces the figure legends with a concise rewrite (17,442 to 7,213 words; per-figure writers, then blind checks against the rendered panels; `Analysis/Hypotheses/LEGEND-CONCISE-20261020T234202Z/`). Results that only the old legends stated are listed in that folder's NOT_COVERED.md. Methods gains the Figure 2C variant-class rule and five bootstrap seeds. r7: r7 applies the triage of `docs/technical/MASLD_RESOURCE_GOAL_AND_EVIDENCE_REVIEW_2026-10-20.md` (checks in `Analysis/Hypotheses/REVIEW-TRIAGE-20261020T141232Z/CHECKS.md`): GNMT's NAFLD association is reported as sub-threshold (anchor p 2.6e-4; PP.H4 0.33/0.71/0.83 at p12 1e-6/5e-6/1e-5), the four Guilliams spatial donors are described as non-MASLD surgical samples sharing H37/H38 with the atlas, the aspect and variant-combination questions are stated as indeterminate in the Discussion, and Supplementary Note 14 records that no Figure 7B duplicate crosses cohort folds. The review's proposed goal rewrite (its section 7) was not adopted; PAPER.md is unchanged. r6 removes, by user decision, the evidence-weighted score, the UniCure Discussion example and four Methods analyses with no reported results (sex-stratified enzyme effects, nearest-gene comparison, plasma-liver protein linkage, assay-native class tests). r5 replaced r4's manuscript text with a sentence-level rewrite for concision and accuracy, with cross-section number conflicts resolved against the source and the author block filled (`Analysis/Hypotheses/MS-CONCISE-20261020T121448Z/`). Since r2 it adds the rewritten manuscript text (Summary, Supplementary Notes, references; `Analysis/Hypotheses/MS-WRAPUP-20261019T231208Z/`) and all supplementary Figure S2-S4 panels redrawn on COLOC r2 and the 774 units, each independently re-derived (`figures/candidates/supp-redraw-20261020T0300Z/INSTALL_PLAN.tsv`). Figure S2E (ACADS) and S2N are withdrawn. k4 P1 is displayed as "Inflammatory" (user approval 2026-10-20; membership unchanged). r2 was promoted 2026-10-17 with user approval and superseded 2026-10-17-r1. Panel-deep review changes
(`Analysis/Hypotheses/PANEL-DEEP-20261017T022400Z/REPORT.md`):

- Decision 1: Figures 1C and 2E draw genes never SuSiE-tested in the other trait class or in a non-EUR study as their own untested state (368 of the 414 enzyme-only genes; 213 of the 340 EUR-only genes).
- Decision 2: Figure 7B draws graded fibrosis against plain expression with matched degrees of freedom (−0.128 macro Spearman), and the Nucleotide Transformer reporter-assay gain is drawn as including zero.
- Decision 3: Figure 7E draws all seven severity arms, each labelled precommitted or post hoc, and draws signed eQTL direction as untestable.
- Decision 4: Figures 6B to 6D label the genetic axis eQTL colocalization and mark the 31 genes with a fine-mapped coding variant at PIP ≥ 0.5; r8 classes are unchanged, and the 7 one-study-untestable genes stay untestable with a footnote.
- Decision 5: Figure 5F prints the donor-level result beside the gene-level null (Stromal ECM, 3 of 4 donors positive, sign-flip P 0.19) and draws Vu et al. as source-dependent.
- Decision 6: Figure 5C computes protein-histology associations within the 46 MASLD participants (60 of 125 tests pass BH).
- Decision 7: Figure 4H plots the donor-first system estimates (rho 0.67, 34/43 same direction).
- Decision 8 (revised 2026-10-21 by user decision D3, `Analysis/Hypotheses/REVIEW2-TRIAGE-20261021T123627Z/DECISIONS.md`): Figure 3C uses five cohorts, with GSE126848 entering through its within-cohort filterByExpr refit (BH over its 15,592 genes; the other four cohorts use the fixed 24,115-gene family). It shows the mixed-sign share of multi-cohort genes with and without GSE126848 (31.5% and 13.5%).
- Decision 9: Figure 5E counts protein members once per primary-assembly gene, as its gene column does (95 testable programs, 22 robust).
- Figure 4D stability counts use the joint donor bootstrap, which keeps each donor's lineages together (76 stable matches, 34 stable mismatches, 7 unstable).
- The ancestry companion of Figure 2D is adopted in the main figure, and the direct-subclass split moves to Figure S2S.

**Current correction state, 2026-10-14.** The adopted 774-unit bulk fit is
now reflected in working Figure 3/4/S3/S4 and Figure 5B/C panels. Figure 3C
includes all five cohort contrasts; Figure 3F has 54 genes and 216 effects.
Figure 4A is complete; after review round 3, 4E uses the three-cohort 364-unit
arm (the 457-record arm is a sensitivity), and 4G and 4H use the 367 disease
units, with no outline in 4H. Verified COLOC r2 is selected for the
genetics panels (436 primary, 560 ABF-only, 996 union). The selected r8
Catalog has completed local review, a molecular-object companion, all
five individual Figure 6 panels and a native-size review assembly.
The Figure 5 vector assembly is rebuilt with current 5B/C/D. Eight missing
PDFF correlations and the remaining 27 recoverable native twelve-trait matrix
pairs are complete: all 66 pairs have terminal finite estimates and positive
native heritability. The matrix preserves successful historical fits, identifies
preprocessing differences and reports nominal p-values. Six unconstrained
correlations above one are qualified by weak heritability and their SEs.
N4's different-visit null is untestable and historical 87 random-null
and 69 stage-BH calls are withdrawn. GSE193066 remains identity-unqualified.
The class comparison's prespecified directional prediction fails; the
opposite-direction hepatocyte difference is post hoc and attenuated.
Original-source diagnosis checking additionally corrects the two GSE174748
disease explants to Cirrhosis. The unchanged three-stage rule excludes them.
H59 (2026-10-15) refits all 117 programs on the 60-donor ambient-correctable
universe with a per-dataset residual variance and a within-dataset
permutation null: zero programs selected (smallest q=0.323375). The 62-donor
equal-variance fit (BICC1 alone, q=0.037633; zero HC3-selected) becomes a
sensitivity. Figure 4C and S4C show the H59 effects; their source panels are
`figures/candidates/review-r3-fixes-20261015T112148Z/F4abcd/repair/` and
replace the 62-donor panels together with this text. Liver Atlas retains all 19 stage exclusions and
descriptively separates four Healthy, three Steatosis, one Cirrhosis and
eleven Unresolved donors. Old 64/62-donor effects are historical.
The complete communication refits are checked: six of the fixed thirteen
stage-context pairs are testable (all q=1), seven untestable. The fixed
26-test program screen has 22 estimable tests and two signed Ductular injury
associations, each on 58 donors. Old stage-adjusted q values are historical;
RNA covariation does not establish a signaling mechanism.
Completed outputs and checks are recorded in
[the October 14 outcomes](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#remaining-resource-corrections-2026-10-14).
Earlier dated paragraphs retain historical inspection states. They do not
override the current figures or imply that completed reruns are pending.

**Scientific review corrections, 2026-10-13.** Figure 2's historical 6.69
coding/splice posterior ratio used 715 truncated positions among 6,415 gated
credible sets. Member-supported coordinate restoration gives 23.34% versus
3.87%, ratio 6.03, with 233 clusters (19 direct, 228 enzyme, 14 shared).
The old bootstrap intervals, noncoding-context uncertainty and 62 robustness
refits are superseded by the adopted corrected-coordinate results below.
The descriptive direction persists for the pooled value, conditional on the
studied GWAS portfolio. H58 (2026-10-15) found that the available fine-mapping cannot separate trait
class from ancestry and fine-mapping coverage: 106 of 207 direct but 82 of 6,208
enzyme credible sets are European; EUR 29.88% vs 6.92% (difference 0.2296,
95% interval 0.0550 to 0.4236), AMR 6.73% vs 4.54% (0.0219, -0.0489 to
0.1065); no trait-class call in any stratum. The 36.49%/77.85%
(cluster-weighted) and 52.11%/60.43% (summed over credible sets)
unresolved-context values reflected failed lookup and are withdrawn (repaired
unresolved 0%).

**Figure 2 corrected-coordinate results adopted by the user, 2026-10-13.**
Job 22218816 completed all four analyses under
`figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/`:
the six coding/noncoding-context comparisons, coding robustness,
trait-directness descriptive estimates and the ungated ascertainment control.
The gated coding/splice fractions are 23.34% versus 3.87%, ratio 6.026671
(pooled and descriptive only; H58 makes no trait-class call);
the 95% physical-cluster bootstrap interval is 2.574593–13.510165 (seed
20260817, 10,000 draws). The corrected robustness menu has 57 comparisons:
one reference, 13 direct-study omissions, 18 enzyme-study omissions, 19
direct-cluster omissions, five window sizes and one LD-block comparison.
The historical 62 included 24 direct clusters from erroneous coordinates;
correcting their count to 19 changes the menu size without dropping a test
type. Ratios across the 57 comparisons span 4.368132–8.275691; every difference
interval has a positive lower endpoint (minimum 0.022373), with 10,000 finite
draws in each. These overlapping sensitivity comparisons are not independent
replications. The ungated ascertainment baseline and matched ratios are
6.178788 and 5.782896; they are separate from the gated 6.026671 estimate.
Ascertainment's native lead-position lookup was already correct; only its
locus clustering needed repair. Permutation diagnostics do not establish a
calibrated hypothesis test. All 21 focused checks passed, and all 57 input/snapshot
checksum records remained unchanged. The user adopted these corrected Figure 2
results on 2026-10-13. Their source directory retains its candidate path; no
named historical release is overwritten. The corrected
[Figure 2D](../figures/main/fig2_genetics/panels/Fig2D_pip_architecture_by_trait_directness.pdf)
and [noncoding annotation companion](../figures/main/fig2_genetics/panels/Fig2_companion_noncoding_annotation_gap.pdf)
are installed from `figures/candidates/fig2-coordinate-panels-20261013T220927Z_22219710/`.
Figure 2D prints 23.1% protein-altering and 73.5% noncoding for direct traits,
and 3.9% protein-altering, 3.6% synonymous/UTR and 92.5% noncoding for enzymes.
These displayed component percentages differ from the combined coding/splice
ratio above. The installed companion, summed over all credible sets of a class, prints
promoter 13.84% / 14.67%, liver-lineage open chromatin 9.26% / 8.44%, liver ABC
enhancer 0.43% / 3.70% and no deposited context 77.90% / 77.69% (direct / enzyme).
The earlier unresolved fractions of 52.11% and 60.43% reflected failed lookup
and were withdrawn by H58 (2026-10-15; repaired value 0%).
Both remain outside the assembled A–I panel set. Source-table,
PDF and visual checks passed; prior affected files were archived, and unrelated
PDFs were unchanged by that October 13 decision. COLOC r2 and other corrected
panels are selected separately on October 14; overall publication-bundle
approval remains separate.

**Figure 5B/5C and class comparison corrected, 2026-10-14.** Both Figure 5
panels are installed from the completed 774-unit rerun: 5B has 3,054 measured
genes, 211 protein-significant (188 concordant, 23 discordant), Spearman
rho=0.457835 and 66% rounded direction agreement. The fixed July display
universe is unchanged. Figure 5C has nine canonical RNA and 23 protein
markers; CD74 enters and HLA-DQA1 leaves the RNA call set relative to
September. Both HLA-DQA1 markers are open. Its primary RNA annotation has
log2FC=0.379801, padj=0.024969 and fails the canonical effect-size floor.
The class comparison now has 281 hepatocyte-specificity pairs (mean
genetic-minus-state difference −0.265515 log2, post hoc nominal sign-flip
p=0.01482985) and 73 residualized protein pairs (difference −0.015601,
nominal sign-flip p=0.487465). The prespecified direction fails. The retained
relabel/rematch null is mismatched and its p/q values are not calibrated
significance. Historical 844-record panel values and the 272-pair test are
superseded. Job 22236311 and independent native-count verification job
22236559 completed; exact memberships and 24,115-gene BH values reproduce.
See the [amended scientific review](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#review-of-this-review-2026-10-13).

**Historical Figure 3F built, 2026-10-05, superseded October 14:** the main PDF showed 61 genes and
244 adjacent-stage effects alongside corrected exact-pair PP.H4. The September
genetics candidate is complete; waiting for COLOC no longer explains the missing
artwork. The PIP column was removed; historical stage fits still
require counting and identity correction. The panel is conditional working
artwork, not a synchronized scientific release. The 61 and 244 are on the 844
records and are historical. On the 774 units the installed source
`figures/candidates/fig3f-no-pip-774-20261012T011128Z/` shows
54 genes and 216 effects, 59 at BH FDR < 0.05 (30 / 13 / 11 / 5 by contrast;
65 and 32 / 21 / 7 / 5 before), with 24,115-gene BH families and 205 / 190 /
157 / 94 units; ten genes leave and three enter (list in `docs/RESULTS.md`). Twenty-nine of the 54 genes are selected only through an F0→F1 mark; of the 30 F0→F1 marks, 18 keep BH q < 0.05 with disease-only F0 and none without GSE130970.

**Historical Figures 3–4 identity qualification, 2026-09-29:** the genotype ledger identifies
repeated bulk libraries and conflicting source labels; GSE193066 identity calls
remain unresolved. Historical bulk roster counts and newly reconstructed fits
are conditional, not verified independent-donor inference. The correction
candidate includes an explicit roster join and a separate provisional identity
sensitivity. That restriction applied to the September candidates; the user
adopted the corrected-count bulk refit on genetic-identity clusters on October
11. Neither that decision nor the historical network analysis verifies
independent donors or qualifies GSE193066; see
[the disposition](technical/FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md).

**Figures 3–4 correction candidate, 2026-09-29.** Participant-level continuum,
content-lineage and protein analyses and fixed-network resampling were run on
compute nodes. Candidate source tables and PDFs are under
`figures/candidates/fig34-scientific-corrections-20260929-v1`;
use the successful directories identified in the
[technical disposition](technical/FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md).
The original program scores/effects reproduce numerically. Corrected histology
labels, cohort effects, uncertainty and narrower claims are incorporated into
working drafts. New independent protein/lineage validation remains unresolved
after the bounded metadata screen. Twenty-five corrected PDFs were installed at
working manuscript paths on 2026-09-30; fixed memberships are unchanged and
the synchronized scientific release is not adopted.

**Historical Figure 4 artwork, superseded for 4E–H on 2026-10-14.**
Figure 4B–H individual PDFs were updated on 2026-10-05 at the user's request
from the reviewed compact assembly candidate. C retains the tracked 24 modules;
D now compares discovery and content-associated lineage across all 117 programs;
E shows linear fibrosis/NAS associations for all 113 eligible weighted programs.
H retains the 18 strict-maxT outlines. Embedded Helvetica 6 pt and fixed page
boxes are verified. E remains conditional on historical bulk identity and does
not summarize nonlinear NAS remodeling. The 113 eligible programs (27 fibrosis
and 4 NAS calls) belong to the 469-record working roster of four cohorts
(GSE130970, GSE135251, GSE162694, GSE174478); Figure 4E has never used the
five-cohort fit, and GSE174478 is outside the 774 units. Two 774-unit candidate
arms exist in `figures/candidates/fig4e-program-axis-774-20261012T010830Z/`
(112 eligible programs each; values in the "Open, in order" paragraph and in
`docs/RESULTS.md`); neither is adopted. Previous PDFs and checks are retained
in `figures/candidates/fig4-assembly-review-20261005-v1/`.

**Correction implementation, 2026-10-05.** The user selected exclusion of
unresolved identity/source-label records, complete ALT/scaffold annotation
correspondence before the primary bulk fit, and preparation of both Figure 7
candidates. The 1,725-DEG fit and 1,833-DEG name-fold sensitivity are historical
library-roster candidates, not corrected donor-level results. A new Figure 1C
candidate gives 414 enzyme-only, 12 direct-disease/liver-fat-only and 10 shared
genes from the completed September genetics (436 total). Final bulk-dependent
Catalog classes and Figure 6 counts wait for counting and identity resolution;
the unchanged September Catalog must not be finalized. The recovered identity
candidate selects 774 source-qualified genetic clusters for bulk analysis
(738 in the strict-depth sensitivity), not verified independent donors.
GSE193066 remains unqualified. Full annotation resolution is still incomplete:
855 correspondences are unresolved in `annotation_resolution_r7/`. After the
earlier r5 review, the user authorized exclusion of seven contamination-scaffold
and 56 complete false-duplication annotations. Sixteen further unresolved
records were classified using exact curated RNA transcript identities and
independent complete exon-coordinate matches: 13 equivalents and three copies.
Three earlier equivalences also became explicit same-haplotype copies.
TRPM3 stays included because none of its exons overlap the invalid region.
On 2026-10-07 the user decided that the 855 remaining annotations are reported
as not testable instead of being resolved one by one. `annotation_resolution_r8/`
is r7 with those 855 relabeled `held_untestable` (291 lncRNA, 260 protein-coding,
remainder pseudogenes and small RNAs; all on alternate or patch contigs); the
other 85,514 annotations and their counting targets are unchanged. APOB, COL3A1,
PKLR, TDO2 and GLUD2 are matched, not held. Counting rule: in the primary fit a
read shared between a tested gene and a held annotation is not counted; a
sensitivity fit credits it to the primary gene when the held annotation has
exactly one candidate (403 of 855); genes whose call differs between the two
fits are flagged. Held genes are counted and reported but never tested. The
776-library recount (`verified_counts_r8_20261007/`, 79,095 genes) and the
refit have run. The user adopted the 774-unit refit on 2026-10-11; corrected
Figure 3/4/S3/S4 panels are installed on 2026-10-14. Canonical calls (BH padj < 0.05,
|logFC| > 0.5), primary source-control labels: 1,854 in 774 units
(`bulk-corrected-controls-20261008T162308Z`) and 1,906 in 738 strict-depth
units, versus 1,725 in the 844-record corrected-controls fit. The 774 units are
genetic-identity clusters, not verified donors. Among the 23,121 genes tested in
both the 774-unit and 844-record fits, 1,553 calls are shared, 126 of the
844-record calls are not called and 189 are new; over all genes 172 of the 1,725
earlier calls are absent and 301 of the 1,854 are new, the difference being
genes tested in only one fit. One gene called in both, IGHA2, changes sign. The
738-unit fit is a four-cohort disease contrast, because GSE126848 keeps 5
controls and no disease unit there. The sensitivity
count arm (single-candidate credit) changes no call among commonly tested genes
(logFC Pearson 0.9999). Reproduction runs agree for all four fits, and an
independent re-derivation on 2026-10-11 reproduces the 774-unit, 738-unit and
leave-one-cohort-out numbers
(`Analysis/Hypotheses/refit-adoption-review-20261011T151011Z/out/`). Leave-one-
cohort-out has run on the 774-unit fit (2026-10-09 block under "What is
official, updated, or open"). Comparison tables:
`figures/candidates/resource-corrections-20261005-v1/refit_comparison_20261008T191739Z/`.
Both six-panel Figure 7 review sets exist;
their exposure and source-table checks have been repeated without new fitting.
All five cohort counting pilots reproduce the native unique gene counts exactly;
they are diagnostic samples, not the full corrected recount. Identity-only fits
on unchanged historical counts give 1,756 canonical annotation-row calls in
774 supported units and 1,807 in 738 strict-depth units. They retain duplicate
counting targets and unresolved identities, so neither is a corrected bulk
gene total.
Provisional Figure 6A/C/E and 28,017 molecular-object records are complete for
local review; their source classes still precede the corrected bulk refit.
The revised Catalog candidate retains 28,017 records, including 14 affected
reference-artifact records with non-applicable bulk RNA evidence; historical
classes remain unchanged. The authorized revised NMF protocol retains 4,517
unchanged counting objects for stage projections and 4,438 for the separate
signature-excluded continuum variant. Fixed k4/k6 weights remain unchanged.
Historical-count normalization checks support eight of ten stage axes and all
ten continuum axes under the fixed reconstruction criteria; stage k4-P2 and
k6-P6 remain unavailable in the revised protocol. Fixed-NMF scoring on the
774-unit counts and the continuum and molecular-system steps listed in the
2026-10-09 block under "What is official, updated, or open" have run; the other
downstream models have not. The revised definition is a sensitivity for review
and has not replaced the original axes.
Execution evidence and
the change log are in the
[existing outcomes document](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#resource-corrections-implementation-2026-10-05).

**Verified:** 2026-08-17 EDT  
**Framing updated:** 2026-08-14 EDT  
**Continuum integration verified:** 2026-08-21 EDT  
**Scope:** standalone MASLD Resource paper  
**Important:** live SLURM state must be queried directly; the snapshot below is dated.

## Executive assessment

The title, approved by the user on 2026-10-16, is **A multimodal human MASLD Resource records liver-trait genetic evidence and
disease-state liver expression side by side for each gene**. The coordination clause was removed under the joint H52/H54
rule (row: H52 untestable / H54 supports), and "resolves the architecture" was
dropped after H58. The 2026-09-08 source-roster analysis is historical: it reported the coordination result in the
already-used GSE193066 evaluation cohort with fixed
program definitions and cohort-fitted nuisance components (105 source participants,
whose independent identities remain unresolved;
89.7% edge sign concordance against a 70.9% baseline, a statistic whose excess also rose under H52's program-group permutation control (mean excess 0.0251, 95th percentile 0.0809); fixed-partition modularity
above both nulls; the ductular-portal community retaining elevated density). The
cohort cannot supply untouched confirmation; further analyses are retrospective,
and GSE193066 counts toward no replication. The two
primary maps are connected to donor-aware cell programs and assay-native
protein, chromatin, and spatial tissue context, then delivered through the
MASLD Gene Catalog. Observability, provenance, and evidence states remain the
methodological guardrails, not the headline. The immediate task is
**reconciliation, not expansion**: propagate COLOC r2, adopted 2026-10-14, through its downstream
consumers, propagate the chosen
five-cohort fragment-count bulk substrate through dependent analyses, regenerate
figures that have changed, and settle one synchronized release. The paper does not need a
new public-data mechanism or follow-up screen result to be complete. The recurring
biological examples are now fixed by question. A provisional direct-NAFLD GNMT
locus bridges inherited variation to the severity-linked
GNMT–MAT1A–CYP2C19 RNA–stage–protein state. The Stromal ECM versus
Ductular injury program contrast separates disease association, tissue organization, and
lineage attribution. Current official COLOC does not support a primary
SuSiE-colocalized mechanism. The aligned ambient-RNA and donor-paired
cross-lineage candidates are independently validated but not adopted. A
Kamzolas-adapted molecular continuum complements recorded fibrosis stage with
stage-adjusted donor-level resolution. Its additive Figure 4F–4H and complete
S3/S4 panel package were adopted into the working figure authorities on
2026-08-24; the synchronized Resource release has not been adopted.

## What is official, updated, or open

Historical selection: Validated update selected for the paper. The
selected bulk analysis is the corrected 774-unit fit; current panels and
October 14 checks supersede earlier candidate-installation notes.

**Public-data hypothesis tests, 2026-10-07 to 2026-10-09: run; results
provisional unless explicitly selected.** The user adopted the 774-unit
five-cohort refit on 2026-10-11; corrected r2 genetics are selected on
2026-10-14. The 42 rows below record the tests and checks
from the hypothesis review that ran on open public data. Every other result is a
candidate. The "Depends on" column names the rows that use the
five-cohort refit (`bulk-corrected-controls-20261008T162308Z`, adopted
2026-10-11) or the selected
COLOC r2 release (`susie_coloc_r2_20260924T053244Z`). Dated candidate headers
retain their earlier adoption state; choosing their substrate does not
independently adopt every hypothesis result. The refit's numbers entered
`docs/RESULTS.md` on 2026-10-11; at that update no other row had entered it,
apart from the PXD051911 wording listed under data corrections. V4's numerical
qualifications were added on October 13 as candidate-only. The 117
Hotspot programs and weights were not changed and the GSE289173 holdout was not
used by those analyses. H90 opened all 48 GSE289173 expression matrices on
2026-10-15; its test was untestable at the gene-coverage gate (RESULTS.md). Classes: Usable = ran, controls passed, can be used as a candidate or
validation. Null = ran, controls passed, no effect detected within the stated
power. Uninformative = a named control failed, so the result is not read either
way. Untestable = a prespecified stop rule fired before the outcome was
analysed. Note = sensitivity result only. The class is a judgment call for H2,
H27, H47, H49 and N3. Paths without a top-level directory are under
`Analysis/Hypotheses/`.

| ID | Class | Question, unit, n, correction, result | Depends on | Path |
|---|---|---|---|---|
| H1 refit | Usable | Does the disease signal hold with one record per selected genetic cluster? Pooled five-cohort limma-voom with quality weights; unit: genetic-identity cluster, n = 774 (138 control, 636 disease), 24,115 genes tested; dataset + inferred sex + disease, seed 20260923. BH padj < 0.05 over tested genes and absolute logFC > 0.5: 1,854 calls (1,333 up, 521 down); interval-null (lfc 0.25) 2,131. Arms: strict NAS 0 controls 1,907; sex-discordant persons dropped (771 units) 1,897; 738 strict-depth units 1,906. Compare: adopted 844-sample fit 1,347 calls (157 control, 687 disease; 23,370 genes; interval-null 1,616); 844-record corrected-controls fit 1,725. With the old labels the pipeline reproduces the adopted fit (26 of 26 checks). GSE213621 supplies 326 of the 774 units, 68 of the 138 controls and 45.1% of the summed quality weight. The 738-unit fit is a four-cohort disease contrast (GSE126848 keeps 5 controls and no disease unit there). Per-cohort disease membership (Figure 3C) now uses all five cohorts: GSE126848 was restored on 2026-10-14 with constant sex omitted (23 control, 7 disease; 6,147 calls in the fixed 24,115-gene BH family at FDR <0.05 and |log2FC| >0.5). Native fragments and an independent group-means refit reproduced all call members. Since 2026-10-21 (user decision D3, `REVIEW2-TRIAGE-20261021T123627Z/DECISIONS.md`), Figure 3C uses the GSE126848 within-cohort filterByExpr refit instead (3,440 calls; BH over its 15,592 genes; `REVIEW2-TRIAGE-20261021T123627Z/analyses/AN3/run-20261021T124631Z/`). The five-cohort union contains 12,180 genes: 8,341 called in one cohort and 3,839 in two or more; 1,211 of the 3,839 have opposite signs between cohorts (946 involve GSE126848), leaving 1,820 / 582 / 186 / 40 single-sign genes in two to five cohorts. With the 6,147 calls the union was 13,761 genes and 1,865 of 4,696 (39.7%) were mixed. The source has 25 male controls and 31 disease records, including ten women; the identity-site screen excludes two controls and 24 disease records (all ten women and 14 men). [Inference] Unequal eligibility may affect representativeness. Units are genetic-identity clusters, not verified donors. Against the 844-record corrected-controls fit: 1,553 calls shared, 126 lost and 189 new among the 23,121 genes tested in both; 172 and 301 over all genes; IGHA2 changes sign (-0.774 to +0.615); lncRNA calls go from 56.3% up to 47.9% up. Re-derived 2026-10-11 in `refit-adoption-review-20261011T151011Z/`. | Adopted by the user on 2026-10-11; seven panels installed 2026-10-12; Figure 3F, 3G, 4E and all five Figure 6 panels installed on 2026-10-14; continuum composition-adjusted coefficients rebuilt and current supplementary panels installed (`figures/candidates/continuum-composition-774-20261012T012137Z/`) | `RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20261008T162308Z/arm_deg_counts.tsv` |
| H1-LOCO | Usable | How many refit calls survive removing each cohort? Unit: cohort removed (5 refits). Calls kept with the same sign: GSE126848 96.8% (1,795 of 1,854), GSE130970 94.4%, GSE135251 94.1%, GSE162694 88.6%, GSE213621 66.7% (1,236). Removing GSE213621 loses 618 calls (510 no longer called, 108 untested); 467 of the 618 are lost only under that removal; 5 calls have the opposite sign at any p. 1,026 calls (55%) survive all five removals. No call called in both fits changes sign. Of the 510 lost calls still tested without GSE213621, 505 keep their sign and effects shrink (slope 0.87 of leave-out on full log2FC). These full refits (66.7% to 96.8% kept) are the values for text; Figure S3C plots the fixed 24,115-gene version, 70.7% to 96.6% kept, Spearman 0.870 to 0.995. | Refit adopted 2026-10-11 | `RNA-seq/results/manuscript_release/candidates/bulk-loco-summary-20261009T143029Z/loco_summary.tsv` |
| H79 | Usable | Do the 1,854 calls hold when genetic ancestry enters the model? Unit: genetic-identity cluster, n = 774 (138 control, 636 disease); adopted design plus PC1 to PC4 from RNA-called genotypes (T1) and the adopted design on 614 EUR-called units (T2); family F1 = {T1, T2}, Holm alpha 0.05; seed 20260923. T1 changed 233 calls (140 lost, 93 gained) against 200 within-cell PC permutations (median 120, 95th percentile 161), Holm p 0.0199. T2 changed 570 (101 lost, 469 gained) against 100 matched subsets (median 409.5, 95th percentile 501), Holm p 0.0199. log2FC over the calls Pearson 0.9972, slope 0.984 (descriptive). PC-0, PC-A, PC-B, PC-C, NC-1 and NC-2 passed. 64 genes meet the descriptive ancestry-sensitive rule (expected false flags median 6, 95th percentile 17). GSE135251 per-cohort membership cannot be ancestry matched (2 EUR controls). The adopted fit stays the release; a Catalog flag is a user decision. | Adopted 774-unit refit; ancestry `rna_genotypes_restricted/obs-20260926T155651Z/ancestry/full-20260929T142807Z` | `H79-20261015T125125Z/results/run-20261015T132031Z/RESULT_H79.md` |
| H81 | Usable | What does each bulk evidence state mean when known changes are added? Unit inside each fit: genetic-identity cluster, n = 774; 50 replicates of within-cohort label permutation plus binomial thinning, 3,200 spiked genes each; no test, bootstrap intervals over replicates (seed 20260923 + r; bootstrap 20261923). Pass. Recall of a uniform 1.0 log2FC change 0.507 to 1.000 per decile; detection limit at 80% recall 0.649 to 0.732 log2FC in deciles 2 to 10, above 1.5 in decile 1; false `supported` 0.0004. The L09-01 equivalence rule is usable only in deciles 7 to 9 (not adopted; Catalog `tested_negative` stays 0). One-cohort effects are called in the pooled fit mainly through GSE213621 (0.8745 at 1.0 log2FC, at most 0.0025 elsewhere); 99.8% of supported one-cohort effects lose a LOCO call; uniform 1.0 log2FC effects keep all five LOCO calls in 90.6%. No call relabelled; a per-decile detection-limit Catalog field is a user decision. | Adopted 774-unit refit; LOCO targets `refit-adoption-review-20261011T151011Z`; per-cohort fits `bulk-context-774-20261014T140904Z_22236311` | `H81-20261015T125233Z/run-20261015T141242Z/RESULT_H81.md` |
| H1-continuum | Usable | Do continuum and molecular-system results hold on the 774-unit counts? Unit: selected genetic-identity cluster (466 in the two primary evaluation cohorts, was 501); 43 systems, BH over 43 per aggregation method, axis and model; 10,000 permutations per axis and cohort. Unchanged from the 844 run: 33 collection-balanced continuum systems, 27 stage systems, 24 supported by both, 18 by maxT. Changed: feature-balanced continuum 29 to 31, permutation BH on all four tests 22 to 20, joint maxT 35 to 33, effective independent systems 15 to 14, stage-adjusted participants 470 to 435. Pooled-window Spearman with the fixed 844 windows 0.849 (not a 774 reconciliation). Run 2026-10-09: setup, core 00 to 05, molecular 10, 20, 31, systems 70, stress 73. Run 2026-10-11 (`logs/continuum_774/runC_20261011T152106Z/stage_status.tsv`): molecular 21, 30, 40, 50, 60 to 65, systems 67 to 72, the 74 renderer, compact 4G and 4H, core 06 to 09; outputs are candidates under roots stamped 20261011T152106Z and panels in `figures/candidates/continuum-774-panels-20261011T152106Z/`, `...165944Z/` and `...170854Z/`. Corrected Fig 4G/4H and S4Q to S4U panels were installed on 2026-10-14; after review round 3 (R05-02), main 4G and 4H use the 367 disease units and 4H has no outline, while the 18 strict maxT calls (435 units, including the 68 GSE213621 controls) remain a distinct supplementary criterion. | Refit counts adopted 2026-10-11; continuum panels installed 2026-10-14 | `logs/continuum_774/run_20261009T171017Z/comparison_844_vs_774.tsv` |
| H1-NMF | Usable | Do fixed NMF axes transfer to corrected counts? Unit: analysis unit, n = 774; fixed k4 and k6 weights, bases not refitted, outcomes not read, seed 20260923. Stage axes available 8 of 10 (k4-P2 0.9897 and k6-P6 0.9874 fall below the 0.99 historical Spearman rule); continuum axes 10 of 10. 4,517 counting objects (stage), 4,438 (continuum). The 855 held-untestable annotations are counted, never tested; 0 unresolved. | Candidate `fixed-nmf-corrected-scores-20261009T141312Z`; original axes not replaced | `RNA-seq/results/manuscript_release/candidates/fixed-nmf-corrected-scores-20261009T141312Z/axis_availability.tsv` |
| V1 | Usable | How much of the bulk signal is perioperative dexamethasone response? Unit: gene; liver donors n = 26 (14 with dexamethasone, 12 without) give 70 response genes (difference in differences). 13 of the 1,854 refit calls are among them: 0.70% (95% CI 0.37 to 1.20%). Overlap exceeds chance: odds ratio 3.04 (1.52 to 5.65), Holm 0.0021 over 3 tests; matched-null p 0.0002 (Holm 0.0006). Sign agreement 2 of 13 (Holm 0.077). 857 calls cannot be tested in the dexamethasone data. 6 of 116 programs flagged (BH q < 0.05). Post hoc sensitivity refit on the 774 units with a 70-gene glucocorticoid-response score as covariate, run 2026-10-12, outside the Holm family (the first prespecification said no refit): 1,900 calls, 1,732 of 1,854 kept (93.4%), 122 lost (119 by |log2FC| alone, none changes sign), 168 gained; log2FC Pearson 0.9936; all 26 genes in both the genetics and transcriptomics classes keep their call; 200 matched random 70-gene scores lose a median of 193 calls (29 to 753), so the change is within what an expression-derived covariate produces. 9 of 521 down calls are in the set (odds ratio 6.78, 2.94 to 13.8, nominal p 1.9e-05), all genes dexamethasone raises; 8 of the 13 set calls are lost under adjustment (FKBP5 -0.587 to -0.340). The disease-versus-control signal does not depend on the score; the FKBP5, TAT, KLF9 group of calls may reflect perioperative exposure differences. No cohort has an exposure variable. | 774-unit refit (adopted 2026-10-11) calls | `V1-20261008T191919Z/results/run-20261008T193139Z/summary.tsv` (no RESULT file); `V1-774-20261012T012943Z/results/run-20261012T013138Z/` |
| V4 | Usable | Does stromal ECM spatial clustering replicate in more diseased Visium tissue? Mean centered residual Moran's I, Holm over 2. HMSMA (HRA007511), unit donor, n = 32: 0.098 (bootstrap 95% CI 0.083 to 0.114), 32 of 32 positive, one-sided exact sign-flip nominal p 2.3e-10 (Holm 4.7e-10 over HMSMA and Vu). Vu 2025, unit array (no donor key), n = 10: 0.147 (0.119 to 0.176), 10 of 10 positive, exact p 0.00098. GSE192741, 4 donors: 0.139, exact p 0.0625, the smallest possible with 4 donors, outside the family. These values adjust for the cell2location hepatocyte abundance (GSE192741, Vu) or a 15-marker hepatocyte proxy (HMSMA), not the 16 abundances of Figure 5F (Stromal ECM excess 0.040 and 0.050), so V4 does not replicate the 5F statistic; Vu arrays are not donors. [Inference] The claimed HMSMA roster of 28 MASLD/MASH participants and four controls depends on inferred exclusion of three controls, without an authoritative person-to-clinical-label join. Stage is not testable (no array-to-clinical key). | Fixed programs; no refit or COLOC input | `V4-20261008T015029Z/summary_donor_tests.tsv` (no RESULT file) |
| H8 | Usable | Which transcripts track ballooning after adjusting for steatosis, fibrosis and inflammation? Unit: person, 3 cohorts (GSE130970 n = 71, GSE267145 n = 99, E-MTAB-9815 n = 58). DerSimonian-Laird meta-analysis, BH over 12,896 genes: 88 at q < 0.05; 62 replicate (same sign in all three, nominal p < 0.05 in at least two; 37 up, 25 down); stratified permutation p 0.006. No gene passes BH within any single cohort; 30 of 62 keep the call with categorical covariates. Positive controls (COL1A1, COL1A2, PLIN2) pass narrowly. | GSE130970 roster and counts from the unadopted r8 refit candidate; the other two cohorts use source counts | `H8-20261009T015529Z/RESULT_H8.md` |
| H26 | Usable | Do eQTL direction and rare-variant direction agree for anchored genes? Unit: gene, 383 anchored genes. Common-variant eQTL direction: 269 calls (138 inhibit, 131 activate), 103 untestable, 11 ambiguous. Exome burden for the same trait: 368 testable, 2 at BH q < 0.05 (CLEC4M, CTSD). CTSD is the only gene supported by both lines and they agree. Sign agreement between lines among 254 genes: 135 (0.53), p 0.35, descriptive. 4 of 5 calibration genes have the expected burden sign at nominal p < 0.05. Burden for ALT or AST in the 1,054 non-colocalized disease-state genes (BH over up to 2,108 tests): 1 gene (NEB, AST) at q < 0.05. | Unadopted COLOC r2 | `H26-20261007T210504Z/run_22092011/summary.json` (no RESULT file) |
| H27 | Usable | Do anchored liver loci colocalize with kidney or heart failure traits? Unit: independent locus, 274 tested per trait; call = PP.H4 > 0.5, SuSiE primary, ABF second arm. Union calls: eGFR 48, CKD 8, HF 10; 51 of 66 locus-trait calls (45 of 56 loci) are new relative to the Jamialahmadi 2024 and Vujkovic 2022 lists. ABF calls exceed all 200 seeded permutations per trait (empirical FDR 0.13, 0.11, 0.18). PC2 (TM6SF2 x eGFR) passed; PC1 (SHROOM3 x eGFR) and PC3 (ABO x HF) failed as distinct signals; no anchored HF positive control exists. Sample overlap is not modelled; 31 of 51 new calls use a UKB liver GWAS. L107 direction values: see data corrections. | Unadopted COLOC r2 | `H27-20261008T230909Z/RESULT_H27.md` |
| N1 | Usable | Is the GSE193066 "possibly mixed" call a technical artifact? Unit: library; 164 libraries, the existing rule flags 132. Mixture method accepted: sensitivity 0.88 at mixture fraction 0.10 (mid depth); false-call rate 0.0 at fraction 0; 0.01 and 0.0 for the deamination and symmetric-flip negative controls. Primary rule: 83 of the 132 not supported as mixed, 49 remain called mixed (strict rule: 99 and 33). | No refit or COLOC input | `N1-20261007T210855Z/summary.json` (no RESULT file) |
| N4 | Retrospective descriptive evidence; different-visit null untestable | Eligibility correction completed 2026-10-14, job 22236285. Unit: source-declared pair/genetic cluster, not verified person. Different-visit pairs, n = 22 (GSE130970 3, GSE193066 19): stage ICC 0.07 (95% CI -0.37 to 0.48); median program ICC 0.71 (0.38 to 0.84); best program activated stellate (PDGFRA) 0.91 (0.77 to 0.94); 102/108 program point estimates above stage. These remain descriptions only. Historical October 7 claims of 69 stage-BH calls are withdrawn because raw bootstrap tails were uncalibrated; 87 random-null calls are withdrawn because eligibility differed. Corrected null requires equal N1 non-mixture eligibility, distinct source/genotype units and exact cohort/stage matching: GSE193066 stage 2 requires 12 endpoints but has 8, stage 3 requires 19 but has 2, stage 4 requires 1 but has none. No valid different-visit null p/BH calls exist. All 164 GSE193066 libraries remain identity-unqualified; no primary replication, ICC or ceiling claim follows. Pairs of unestablished type, n = 28 (GSE135251 3, GSE213621 25): stage ICC 0.64 (0.32 to 0.81), median program ICC 0.11; zero calls against the eligible exact-stage null, BH over all 117 fixed programs. The 19 GSE193066 pairs remain the N1-eligible historical subset; the other 39 source-declared pairs were not analysed. | Historical fragment counts; no corrected primary refit or COLOC input; fixed 117 programs; no additional N4 job pending | [Corrected N4 results](RESULTS.md#n4-within-person-agreement-eligibility-correction-2026-10-14); `Analysis/Hypotheses/N4-eligibility-corrected-20261014T140730Z/summary_by_group.tsv` |
| H2 | Null | Does disease move anchored genes the way the risk allele does? Unit: gene; 16 anchored genes with a disease call (of 285 tested); one pooled test, no correction. 11 risk-concordant: 0.69 (95% Wilson CI 0.44 to 0.86), exact binomial p 0.21, permutation p 0.20 (B = 100,000, seed 20260923). The test cannot separate 0.5 from about 0.85. Genes sharing a locus are counted separately, so both p values are optimistic. No informative literature control exists; calibration is the permutation null. | 774-unit refit (adopted 2026-10-11) and COLOC r2 | `H2-20261008T224908Z/RESULT_H2.md` |
| H3 | Null | Do liver risk-variant genotypes predict program scores beyond fibrosis stage? Unit: person (up to 969, 7 cohorts, for PNPLA3). 5 variants x 2 programs = 10 tests, Holm, plus a permutation family-wise p. Smallest pooled p 0.083 (MTARC1 x lipid droplet), Holm 0.83, permutation family-wise p 0.49. Effects above about 0.1 to 0.2 SD per allele are excluded for PNPLA3, GCKR and TM6SF2, and above about 0.15 to 0.27 SD for MTARC1 and HSD17B13. Positive control passed (each PNPLA3 G allele raises stage by 0.20, one-sided p 0.009). HSD17B13 is callable in 49% of persons: lower confidence. | Four cohorts use 774-unit counts (refit adopted 2026-10-11); three use historical fragment counts | `H3-20261009T180331Z/RESULT_H3.md` |
| H10 | Null | Do enzyme-only loci have PNPLA3-sized cirrhosis effects per unit ALT? Unit: independent locus. 436 anchored genes in 289 loci: 270 enzyme-only (3 MHC, not classified) and 19 with a direct-trait gene. Cirrhosis log-OR per unit ALT effect, inverse-variance weighted, locus-block bootstrap (10,000, seed 20261007): enzyme-only 0.20 (95% CI -0.08 to 0.49); direct-involved (n = 9) 1.58 (0.84 to 3.76); 5 positive controls 3.35 (2.37 to 3.61). Of 267 non-MHC enzyme-only loci, 72 are evaluable (BH over the 72): 1 PNPLA3-class, 17 no cirrhosis effect, 53 uninformative, 1 opposite direction; 195 not evaluable. Ghouse replication outcome 0.21 (-0.09 to 0.53), 2 PNPLA3-class. GPT negative control passes. | Unadopted COLOC r2 | `H10-20261007T210859Z/summary.json` (no RESULT file) |
| H17 | No call (in vivo arm cannot detect canonical lipid-droplet genes) | Source-corrected refit completed 2026-10-14 (job 22237956): do iPSC and in vivo hepatocytes share a lipid-loading response? Unit: participant for iPSC-hepatocytes (193, paired); donor for in vivo hepatocytes (55: 9 control, 46 steatotic; GSE244832 and GSE202379); 14,975 genes. Same sign among 7,628 iPSC q < 0.05 genes: 0.447 (donor-label permutation 95th percentile 0.603, lower-tail p 0.21, 1,000 permutations). Spearman rho -0.093 (donor bootstrap 95% CI -0.177 to 0.036). Significant in both: 3 against 3.06 expected (p 0.67); core 1 gene (MCTP2). All five prespecified support criteria failed; iPSC positive control passed. The in vivo arm calls 7 genes at q < 0.05 (6 in the 14,975-gene shared universe) and misses CIDEC (q 0.71), PLIN2 (q 0.9985) and PLIN4 (q 0.63), so it cannot detect canonical lipid-droplet genes. No opposite-direction claim: the gene-shuffle p 1e-4 treats correlated genes as independent and the steatosis-only sensitivity reverses sign. | Bulk signature from the 774-unit refit (adopted 2026-10-11); in vivo stage labels from the 2026-10-14 source-diagnosis correction roster, the correction used for Figure 4C | `H17-source-stage-corrected-20261014_22237956/decision.tsv` |
| H22 | Null | Do liver-fat and cirrhosis heritability coefficients differ within liver cell types? Unit: annotation (cell-type peak set); difference of `tau*` between PDFF and cirrhosis, S-LDSC block jackknife. Four prespecified families at BH q < 0.05: 0 of 8, 0 of 7, 0 of 8, 0 of 8. The SE of a difference is 0.94 to 1.16 `tau*` units, so only differences of about 2.6 to 3.2 were detectable [Estimate]. Controls passed (ALT hepatocyte-specific `tau*` z 2.61; height z -0.65; null calibration 0.050 and 0.0275 at absolute z of 1.96 or more). | Unadopted candidates (H24 GSE296875 pseudobulk, endpoint lock of 2026-08-25); no refit or COLOC input | `H22-20261009T015913Z/RESULT_H22.md` |
| H23 | Null | Do non-EUR ancestries change the fine-mapped lead variant beyond sample size? Unit: locus window per trait and study set (20 to 26 paired loci per test); 12 tests, BH. The lead variant is replaced less often with multi-ancestry data than between two EUR studies: 0.346 vs 0.731 (MVP ALT), 0.150 vs 0.800 (MVP AST), 0.238 vs 0.524 and 0.385 vs 0.654 (4WAY ALT, AST); one-sided p 0.981 to 1.000, q = 1.000. Against size-matched EUR (exploratory, nominal): smallest BH q 0.36 over 24 tests, so a moderate excess is not excluded. Controls partial: TM6SF2 and MBOAT7 not evaluable, ALDH2 East Asian control not testable, gene endpoint untestable. | Unadopted COLOC r2 and fine-mapping releases | `H23-20261009T140416Z/RESULT_H23.md` |
| H31 | Null | Do SERPINA1 Pi*MZ carriers show a raised ER-stress program? Unit: person, 21 `Pi*MZ` vs 461 `Pi*MM` in 3 cohorts; UPR score ~ MZ + stage + sex, HC3 errors, pooled (I-squared 0). Beta -0.248 SD (95% CI -0.608 to +0.113), p 0.18, permutation p 0.19; one prespecified test. Low power: only shifts of about 0.5 SD or more were detectable at 80% power. The stop rule of 20 MZ persons is cleared by one. No control failed. | Four cohorts use 774-unit counts (refit adopted 2026-10-11) (as H3) | `H3-20261009T180331Z/RESULT_H3.md` |
| H34 | Null; candidate-only | Source-diagnosis stage refit 2026-10-14 uses unchanged Y-negative fractions, sex calls and 50 matched autosomal control sets. OLS on stage + dataset, HC3, BH over two immune cell types. Macrophages +1.487719 percentage points per stage (95% HC3 interval -1.413582 to 4.389019), nominal p=0.314886, q=0.482539; T cells +0.486050 (-0.870552 to 1.842651), p=q=0.482539. Actual fitted n=18 and 22 across three datasets; eligible n=19/23 includes one singleton-dataset donor dropped by the original HC3 rule in both old and new fits. Zero primary calls. Within-donor scores and cache hashes are unchanged. Earlier stage estimates and outside-family slopes are historical. | Qualified source diagnoses; no bulk/COLOC input | `H34-stage-source-corrected-20261014_22237596/stage_association.tsv`, `negative_control_stage_slopes.tsv`, `stage_estimate_comparison.tsv`, `CHECKED.txt` |
| H38 | Null | Does bulk RNA show an IgA-high subgroup of patients? Unit: person, 774 (5 cohorts); beta-binomial with cohort-specific dispersion. Subgroup likelihood ratio 3.14, bootstrap p 0.24 (B = 500, seed 20260923); odds ratio 1.5 against a prespecified minimum of 2.7: not supported. IgA fraction rises with stage (713 persons, 4 cohorts, one prespecified test): 0.083 per step (95% CI 0.007 to 0.159), p 0.033, permutation p 0.035, I-squared 76%, and depends on all four cohorts. Secondary, BH over 2: NAS q 0.0104, disease vs control q 0.0104. Controls passed. | Roster from the 774-unit refit (adopted 2026-10-11) | `H38-20261009T180100Z/RESULT_H38.md` |
| H47 | Null | Does signature reversal in LINCS recover PPAR and FXR agonists? Unit: compound (4,186 with liver-line signatures); one-sided Mann-Whitney, BH over 3 testable sets, matched null of 1,000 random queries. AUC (95% CI): PPAR agonists (19) 0.52 (0.37 to 0.67), q 0.38; Lefever set (9) 0.55, q 0.38; Khazaei-Poul set (12) 0.66 (0.51 to 0.82), q 0.062, matched-null p 0.045. Recovery needed q < 0.05 and matched-null p < 0.05; none met it. FXR agonists: 3 compounds, untestable. This calibrates the method; it is not a drug-discovery result. The Khazaei-Poul genes derive from two of our five cohorts. | Query from the 774-unit refit (adopted 2026-10-11) | `H47-20261008T214453Z/RESULT_H47.md` |
| N3 | No qualified corrected calls | Source-qualified refit completed October 14: 70 biological donor records across four datasets, Healthy/Steatosis/Steatohepatitis/Cirrhosis 14/11/34/11; only SRR14602057/58 change from Steatosis to source-documented Cirrhosis. Native limma-voom model retains dataset + inferred sex + ordinal stage, fixed 51 targets and controls, and 1,000 within-dataset permutations per raw/corrected fit (seed 20260923). BH families contain 45 raw and 40 corrected expression-eligible targets. Ordinary-BH calibration fails in both fits (mean nominal-positive fractions 0.078889/0.078125), so primary evidence uses empirical permutation BH: nine raw calls, zero corrected calls. Ordinary-BH counts 12/six are diagnostic only. Corrected LEPR beta=0.5440, ordinary q=0.01486, primary q=0.07421; the earlier 'LEPR robust' conclusion is superseded. Both fixed control composites pass. All 208 included libraries were corrected using native empty-droplet profiles; no included donor has passthrough. The fixed historical nine-target ambient diagnostic roster is preserved. Contamination rises with stage in the nominal primary diagnostic (beta=0.1773 logit per step, 95% CI 0.04938–0.3052, p=0.007348). Limits remain incomplete ambient removal, no doublet analysis and no second correction method. Zero detected corrected calls is not absent effects, a powered negative or clinical target evidence. | Source categories corrected; both full fits and independent terminal checks complete. October 9/13 numbers are historical. | `N3-source-diagnosis-20261014/verification.json`; `corrected/ambient/old_vs_new_all_targets.tsv`; historical `N3-20261009T140326Z/RESULT_N3.md` and `N3-ambient-20261013T060951Z/RESULT.md` |
| V3 | Null | Do stromal ECM and ductular programs replicate in new donors? Unit: donor, GSE256398 hepatocyte pseudobulk. C1 stage in non-cirrhotic donors (n = 13); C2 disease vs healthy (n = 17). 4 tests, Holm, HC3. Stromal ECM: C1 0.30 (95% CI -0.71 to 1.31), C2 0.18 (-0.84 to 1.19), Holm 1.0. Ductular injury: C1 0.69 (0.06 to 1.33), p 0.035, Holm 0.139; C2 1.08 (-0.09 to 2.25), p 0.068, Holm 0.203. No test passes Holm; power is low at these donor counts. No positive control was named. GSE289173 was not used. | Fixed programs; no refit or COLOC input | `V3-20261008T015239Z/results/v3_primary_family.tsv` (no RESULT file) |
| V0 | Usable (audit) | Does V3 hold under the H90 detection rule? Rerun of V3 labelling and pseudobulk reproduces V3 exactly only on Intel v4 nodes (labels depend on CPU type). Hepatocyte coverage at >= 1% detection: Ductular injury 20 genes / 0.910, Stromal ECM 41 / 0.841; both stay testable, V3 verdicts unchanged. | Fixed programs | `Analysis/Hypotheses/V0-20261021T122533Z/RESULT_V0.md` |
| V5 | Run 2026-10-21; not adopted | Is the Ductular injury program higher in MASLD in a new Flex cohort? GSE298719, unit donor, 6 normal v 7 F2-F4 MASLD, `score ~ disease + sex + log10(n_hep)`, HC3, Freedman-Lane 19,999. Beta 1.35 SD (95% CI 0.29 to 2.42), permutation p 0.006; robust to plate (within-plate beta 1.69). Not specific: NC3 at the 86th percentile of 1,000 matched random sets (rule: 95th). NC1 not available; Stromal ECM untestable (coverage 0.424). Figure rule not met. V5P (post hoc, V3 known, CIs only): fixed-effect summary with V3 on 18 shared genes 1.24 (0.71 to 1.78). | Fixed programs | `Analysis/Hypotheses/V5-20261021T122723Z/RESULT_V5.md` |
| V6 | Untestable | Stromal ECM Vinculin+ vs Albumin+ within ROI, GSE255312 GeoMx. Compartment orientation fails (Albumin 9/24, Vinculin 4/24, CD45 11/24 ROIs; rule 20/24); ECM coverage 0.154. | Fixed programs | `Analysis/Hypotheses/V6-20261021T234432Z/RESULT_V6.md` |
| V7 | Descriptive (phase 1) | CosMx observability and in situ lineage of on-panel members (Fukumoto, GSE292268, Govaere). All 117 programs not measured (top panel weight 0.45). FGF13 hepatocyte expression exceeds segmentation spillover in 8 of 9 evaluable sections, MIR4435-2HG in 7 of 9; Stromal ECM members mostly mixed or spillover. Sections are units; no p-values. SCP3359 Xenium (18 donors, 477-gene panel): gate fails, top program weight 0.164, Ductular 0 members (`V7-20261021T235717Z/D4_SCP3359_GATE.md`). | Fixed programs | `Analysis/Hypotheses/V7-20261021T235717Z/RESULT_V7_phase1.md` |
| V8 | Untestable | GSE344087 CITE-seq gate: no-fibrosis donors have 0, 57 and 7 hepatocytes (rule: 30 in all 3). | Fixed programs | `Analysis/Hypotheses/V8-20261021T235438Z/RESULT_V8.md` |
| V9 | Untestable | GSE309979 GeoMx: no count matrix, no patient key, condition codes not decoded. | none | `Analysis/Hypotheses/V9-20261021T235557Z/RESULT_V9.md` |
| H9 | Uninformative | Do published locus classes predict Hotspot program or cell type? Unit: classed locus (Jamialahmadi 2024: 10 discordant, 13 concordant). 610 tests in three arms, BH within arm, class and level: 0 supported. PC1 failed (436 anchored genes against cell-type unions: hepatocytes OR 1.30, second to cholangiocytes 1.34). PC2 failed (0 of 20 lipid-class genes in any of the 5 lipid hepatocyte programs). With 10 vs 13 loci a contrast needs at least 8 hit loci to pass BH over 117 programs; the most observed was 2. Weng 2026 arm not run (no open locus-class table). | Arm A: unadopted COLOC r2; programs fixed | `H9-20261009T140001Z/RESULT_H9.md` |
| H21 | Uninformative | Which Hotspot programs carry cirrhosis signal beyond liver fat? Unit: gene set (MAGMA competitive test); 116 testable programs x 4 traits = 464 tests, BH per trait: 0 at q < 0.05 in any trait (smallest q: ALT 0.142, PDFF 0.142, NAFLD 0.616, cirrhosis 0.533). Cirrhosis-minus-PDFF difference: 0 of 116. Positive control failed: the best of 5 lipid hepatocyte programs for PDFF has p 0.041 (required p < 0.01), so "not liver fat" cannot be read. | No refit or COLOC input | `H21-20261007T211640Z/results/SUMMARY.json` (no RESULT file) |
| H35 | Uninformative | Does the Alu editing index change with fibrosis stage? Unit: person, n = 552 (5 cohorts); beta-binomial per cohort, DerSimonian-Laird; one confirmatory test. Odds ratio per stage step 0.990 (95% CI 0.978 to 1.003), p 0.143, permutation p 0.133; same sign in 4 of 5 cohorts; I-squared 67%. Positive control failed: ADAR expression did not track the editing index (+0.011 per log2CPM, one-sided p 0.45, I-squared 94%). | Roster, refit and COLOC release listed as unadopted in the file header | `H35-20261009T180343Z/RESULT_H35.md` |
| H41 | Uninformative | Does any hepatocyte identity TF knockdown match the disease signature? Unit: knockdown gene (HepG2 only; hairpins are replicates); 12 of 17 candidate TFs present; BH over 12; 10,000 label permutations. Positive control failed: HNF4A weighted connectivity score 0 (percentile 0.497, rule 0.90 or more). Negative control passed. No TF called. CEBPB passes BH (q 0.028, plate-matched p 0.0054) and is reported, not called. | Disease signature from the 774-unit refit (adopted 2026-10-11) | `H41-20261009T015336Z/RESULT_H41.md` |
| H42 | Uninformative | Does PNPLA3 I148M dosage shift zonation toward portal? Unit: person, 758 callable of 774 (442 carriers); within-cohort OLS with sex, hepatocyte fraction, stage and group; DerSimonian-Laird; one primary test. One extra G allele: -0.003 within-cohort SD (SE 0.078), p 0.97; Wald 95% CI -0.156 to +0.151. Portal shifts above about 0.16 SD per allele are excluded. PC-B failed as written (stage coefficient on the score -0.017, p 0.66); PC-A passed (80.8% of zonation genes in the expected zone). | Expression and composition from the 774-unit refit (adopted 2026-10-11) | `H42-20261009T180511Z/RESULT_H42.md` |
| H48 | Uninformative | Which mouse hepatocyte knockouts move toward the human MASLD signature? Unit: knockout gene (Saunders 2025: 202 targets, 187 tested). Primary arm S1 (single-cell hepatocyte signature, 88 donors, 3 genes at BH < 0.05): 3 knockouts toward and 15 away at N1; 0 toward and 14 away at N2. N1 calibration 8.6% of null draws at p < 0.05 (5% expected). No named control passes N1: Pten and Cpt1a (expected toward) q 0.096 and 0.59; Srebf1 and Plin2 (expected away) point toward. The expected directions were recalled, not re-checked. S2 (bulk): 18 toward, 5 away. | S1: single-cell atlas. S2: historical 844-sample bulk fit, not the 774 refit | `H48-20261008T014825Z/step3/summary.json` (no RESULT file) |
| H49 | Uninformative | Are early-MASLD genes (Zhao 2026) induced by steatogenic drugs in hepatocytes? Unit: culture condition (2 treated vs 2 vehicle per drug; Open TG-GATEs primary human hepatocytes, high dose, 24 h). 10 drugs x 2 directions = 20 tests, BH; 10,000 random gene sets. 7 of 20 concordant at q < 0.05 (5 of 10 drugs). 6 non-steatogenic controls shift the same genes: steatogenic minus control 0.124 (up) and 0.169 (down), exact p 0.59 and 0.61. The vehicle-versus-vehicle control failed (45.5% and 25.6% with p < 0.05; rule 10% or less), so the prespecified fallback threshold was applied. Positive control passed through amiodarone only. | No adopted release; early-gene list is Zhao 2026 | `H49-20261009T015504Z/RESULT_H49.md` |
| H7 | Untestable | Is program coordination present at F0 and F1 in single cells? Stop rule fired: 14 donors with F0 or F1 (6 F0, 8 F1) among 46 donors with a documented stage; the rule needs 15. No coordination statistic was computed. | No refit or COLOC input | `H7-20261008T015239Z/stop_rule.tsv` (no RESULT file) |
| H24 | Untestable | Do non-hepatocyte regulatory layers explain direct loci that lack bulk eQTL colocalization? Unit: donor (GSE296875, 39 donors, pseudobulk); family m = 97. 0 of 8 loci (0 of 3 strict) have a two-layer call. Stop rule S1: no lineage passed the informativeness check, so layer 2 never ran. Stop rule S2: GTEx v10 liver recovered 0 of 10 testable control loci at PP.H4 > 0.5. | Unadopted COLOC r2 | `H24-20261008T211750Z/RESULT_H24.md` |
| H25 | Untestable | Is H3K27ac allelic imbalance higher in high-grade livers? Unit: participant (GSE267119 CUT&RUN; 35 high grade, 64 low grade). 0 of 35 and 0 of 64 participants have 10 or more callable heterozygous credible-set variants in disease-gained regions; the rule needs at least 20 per group. Median callable variants per participant 2. No imbalance statistic was computed. | Credible sets from the unadopted COLOC r2 | `H25-20261009T020111Z/RESULT_H25.md` |
| H30 | Untestable | Are FOXO1 hotspot mutations detectable in bulk biopsy RNA? Unit: person, 774. Rule: minimum detectable allele fraction of 0.10 or less (alpha 0.001, power 0.8) in at least 30 persons overall, 10 in F0 and 10 in F3 to F4. Observed 10, 5 and 1; at 0.20 it still fails (55, 19, 8). Median depth at the hotspot sites is 1. No alt-allele count was used and no stage test was run. | Roster from the 774-unit refit (adopted 2026-10-11) | `H30-20261009T175911Z/RESULT_H30.md` |
| H36 | Untestable | Do stage-associated splice clusters colocalize with liver splicing QTLs? Unit: person, 713 staged (4 cohorts); BH over 33,711 clusters. 73 clusters change usage with stage (q < 0.05 and shift of at least 0.10); 10 within-cohort permutations gave 4 calls in total. Splicing-QTL arm: stop rule fired, 15 colocalized introns map to a tested cluster and 20 are required; none of the 15 and none of 260 tested pairs lies in a stage-called cluster. The joint hypothesis is not supported, and these data could not have supported it. | 774-unit refit (adopted 2026-10-11) and COLOC r2 | `H36-20261009T180345Z/RESULT_H36.md` |
| H40 | Untestable | Does an alcohol-specific signature exist in MASLD-labeled biopsies? Stop rule S4 fired. Training cohorts GSE142530 (n = 27) and GSE143318 (n = 25); the first signature has 1,032 genes. Held-out residual AUC after removing HNF4A-target and collagen scores: 0.79 and 0.64; the rule needs 0.75 with lower bound above 0.5 in both folds. Alcohol status is close to collinear with advanced fibrosis in these cohorts [Inference]. The 774 refit units were not scored. | Refit units never scored | `H40-20261008T210215Z/RESULT_H40.md` |
| H46 | Untestable | Do ten adipose genes replicate in GTEx by steatosis, adjusted for BMI? Stopped before any adipose expression was read: BMI is absent from the open GTEx subject files (v8 and v10 hold SUBJID, SEX, AGE band, DTHHRDY only). Exposure counts only: v8 has 151 donors with a liver pathology note and visceral RNA-seq (72 steatosis, 79 not); v10 has 163 (75, 88). BMI is probably in dbGaP phs000424 [Unverified]; dbGaP is excluded by the open-access rule. | No refit or COLOC input | No result file |
| H50 | Untestable | Do ductular-program plasma proteins add to published fibrosis panels? Stop rule: at least 5 of 22 Ductular injury genes must be detected. Detected: 0 in PXD051911 plasma (578 protein groups), 3 in its liver proteome (GLS, ITGAV, MAP2), 2 in the Olink assay list (no stage labels). No stage statistic was computed. Positive control (published 5-protein panel; n = 136, 33 with F2 or higher): cross-validated AUC 0.80 (95% interval 0.69 to 0.88), gain over covariates 0.175 (0.082 to 0.275); not independent, the panel was derived in this cohort. | Fixed programs; no refit or COLOC input | `H50-20261009T135841Z/RESULT_H50.md` |
| N2a | Untestable | Does X-inactivation skew differ between two biopsies of one woman? Unit: woman. Primary set (different visits), 17 women: 9 have at least 5 informative X sites (the rule needs 15) and the median minimum detectable change is 0.212 (rule 0.20 or less); both fail. Secondary set, 15 women: 14 have at least 5 informative sites (needs 15); median 0.078. No X-linked ratio or test was computed. | Sex calls from N1 or the refit's inferred sex | `N2a-20261009T202200Z/RESULT_N2a.md` |
| N5 | Untestable | Does the protective-genotype signature rank compounds and knockouts? Stopped at the positive-control check; nothing was scored. No HSD17B13 or PNPLA3 perturbation exists in LINCS GSE92742 or GSE70138 (0 signatures, any line); 0 of 39,321 compounds name either gene; neither is among the 202 Saunders knockouts. H3 carrier signatures are weak for a 50-up/50-down landmark query (of 932 landmark genes, 23 for PNPLA3 and 26 for HSD17B13 reach FDR < 0.05). | H3 signatures on the 774 roster (refit adopted 2026-10-11) | `N5-20261009T201448Z/RESULT_N5.md` |
| N6 | Note | Does a statin-exposure score change the refit calls? Unit: patient (GSE130991, 910 arrays: 173 statin, 737 not). Out-of-fold AUC 0.783 (95% CI 0.751 to 0.814) for the 50-gene rank score and 0.923 for ridge. In the five RNA-seq cohorts the score rises with fibrosis stage (0.14 SD per stage, Holm p 9e-6), not with disease vs control (0.10 SD, p 0.31). As a covariate it changes 55 of 1,854 refit calls (20 gained, 35 lost, no sign flips) against 6 to 36 for 20 random covariates; all but 4 changed genes lie within 0.05 of the logFC 0.5 cutoff. 46 of 50 genes are up (a cholesterol-synthesis response), so the score does not show that a person took a statin. Array-to-RNA-seq transfer cannot be measured. | RNA-seq step uses the 774-unit refit (adopted 2026-10-11) | `N6-20261009T135702Z/RESULT_N6.md` |

H3 and H31 share one result file. Dropped before running: H4 as a discovery
(kept as V1), H5, H13, H14, H15, H16, H18 as a hypothesis (kept as V4), H19,
H28, H29, H32 as its own test (folded into H3), H33, H37, H43, H44, H45, N7, the
N6 smoking arm, the N2 chromosome-arm arm, and the N3 thyroid-signaling score as
a claim. Merged: H11 and H20 into H26, H12 into H24, N8 into H26, H39v into V2.
Run on 2026-10-11 as candidates (dated paragraph below): nine-cohort bulk
composition acceptance and the Fig S3 stage summaries on the 774 units, and the
continuum steps named in the H1-continuum row. The historical N3 ambient-RNA check was run on 2026-10-13
(`Analysis/Hypotheses/N3-ambient-20261013T060951Z/`; N3 row above). Its
source-diagnosis refit and terminal checks completed October 14; zero corrected
targets pass the required empirical-BH fallback. The earlier estimates are
historical, not current stage results. The V1
glucocorticoid-covariate refit was run on 2026-10-12 as a post hoc sensitivity
analysis (`Analysis/Hypotheses/V1-774-20261012T012943Z/`; V1 row above). The LDSC
genetic-correlation rerun with the corrected PDFF file was run on 2026-10-12
(`GWAS/ldsc/rg_pdff_corrected_20261012T012936Z/`; values in the Pazoki bullet
below). Run on 2026-10-12 as
candidates: Figure 3F on the 774-unit stage table, the Figure 4E program models
with composition adjustment, Figure 3G in the working layout, Figure 6 from r8,
the continuum composition-adjusted coefficients and the NAS0-reference NMF refit ("Open, in order" paragraph below).

**Data corrections, 2026-10-09 and 2026-10-10, made at the user's instruction.**

- *PXD051911 wording in `docs/RESULTS.md` (lines 1050 to 1057).* The 27
  initial-to-surgery-visit pairs differ in age by 0 years in 12 pairs and by 1
  year in 15. The 30 follow-up pairs without a surgery-visit grade differ in age
  by 2 to 4 years; 2 of the 30 lost more than 10 BMI points, and the other 28
  give 0.4818 stages (interval not computed). The 0.5092 and 0.4655 values are
  unchanged.
- *Pazoki 2022 PDFF (GCST90267352) allele labels were reversed on disk.* With no
  allele-frequency column, a label swap equals negating beta. Corrected in
  `GWAS/MR_Data/GCST90267352_PDFF_Pazoki2022.tsv.gz`,
  `GWAS/MR_Data/hg19/GCST90267352_PDFF_Pazoki2022_hg19.tsv.gz`,
  `GWAS/finemapping/data/sumstats/2022_36402844_PDFF_EUR_preprocessed.tsv` and
  `GWAS/ldsc/munged/pdff.sumstats.gz`; the originals are kept as
  `<name>.as_downloaded_20261010T024403Z` and the record is
  `GWAS/orientation_checks/pazoki_pdff-20261010T022251Z/`. The copies derived
  from them were corrected on 2026-10-10 (record
  `GWAS/orientation_checks/pazoki_pdff_derived_copies-20261010T145917Z/`; the
  previous file sits beside each as `<name>.before_..._20261010T145917Z`): labels
  exchanged in the fine-mapping input
  `GWAS/finemapping/data/sumstats_af/2022_36402844_PDFF_EUR_af.tsv`, in the 11-row
  lead-SNP table
  `GWAS/finemapping/data/lead_snps/2022_36402844_PDFF_EUR_preprocessed_leadSNPs_merged.tsv`
  and in `GWAS/finemapping/data/sumstats/Pazoki_PDFF_reformatted_hg19.tsv`. The
  source positions are GRCh37 and had been lifted a second time in both hg19
  copies. The position is restored in 8,895,197 of 8,975,721 rows of each copy and
  `hg19_chr` in 2,091 rows of the MR_Data copy. The share of positions that also
  occur in the other two PDFF studies (`2021_34128465`, `2021_34957434`) is 92.65%
  after the fix and 1.20% before. The sha256 on line 7 of
  `GWAS/finemapping/config/gwas_af_sources.tsv` now names the corrected `_af.tsv`
  (`aaaadfe6...`, was `3f0942c9...`). The uniform35_v4 run manifest still holds the
  old hash because it records what that run read; `01_prepare_study.R` asserts
  the hash, so a rerun of that run's prepare stage stops until a new manifest is
  written. The two hg19 copies keep 8,975,721 of the 9,057,071 source rows and the
  81,350 missing variants were not restored. In those 8,975,721 rows, 327,184 share
  chromosome, alleles, beta, se and p with another row at a different position
  (not investigated).
  Actual SuSiE and September COLOC r2 inputs used the coordinate-correct
  AF-augmented copy, not either double-lifted hg19 copy. Its correction changed
  allele labels only; unsigned PP.H4 and PIPs are invariant to the uniform sign
  reversal. C2orf16 and RABGEF1 are the two of 436 primary genes with Pazoki
  support, and each also passes the primary threshold in another study.
  Coordinate uncertainty is withdrawn for those inspected pipelines;
  direction-sensitive/model outputs beyond them remain untraced.
  Colocalization PP.H4, SuSiE PIPs, S-LDSC and MAGMA do not depend on the sign.
  The LDSC genetic-correlation logs for PDFF in `GWAS/ldsc/rg/` were computed from
  the reversed file and keep the old sign. Rerun with the corrected file on
  2026-10-12 (`GWAS/ldsc/rg_pdff_corrected_20261012T012936Z/`: `ma_pdff.log`,
  `ma_finngen_nafld.log`, `ma_ghodsian_nafld.log`, `h2_pdff.log`): PDFF with
  Ghouse cirrhosis rg +0.5713 (SE 0.07325, nominal p 6.243e-15), cross-trait
  intercept +0.007736 (SE 0.009169); PDFF with FinnGen NAFLD rg +0.7436 (SE
  0.112); PDFF with Ghodsian NAFLD rg +1.1925 (SE 0.2360), of which only the
  sign is usable because that trait's observed-scale h2 is 0.0016 (SE 0.0006).
  The old logs give the same sizes with the opposite sign. Eight other PDFF
  pairs failed in both the old and October 12 runs. The installed LDSC code
  filtered incompatible variants but reused the unfiltered allele vector
  during alignment. On October 14, native-rule prefiltering resolved all
  eight: MVP NAFLD rg 0.7389 (SE 0.0687), MVP cirrhosis 0.6155 (0.0891),
  FinnGen cirrhosis 0.5890 (0.1719), UKBB ALT 0.4979 (0.0466), MVP ALT
  0.3790 (0.0512), UKBB AST 0.1860 (0.0493), MVP AST 0.2219 (0.0522) and
  FinnGen obesity 0.2556 (0.0498). FinnGen cirrhosis remains weakly heritable
  (observed h2 0.0026, SE 0.0012), so its rg is imprecise. The candidate drops
  12,123 duplicated rsIDs from each input and 118,338–125,501 incompatible or
  ambiguous shared IDs, preserving native allele labels and Z scores.
  Job 22236252 completed all eight fits; job 22236390 verified all 256 allele
  combinations against installed LDSC and checked unchanged input hashes.
  Values and nominal p-values are in
  `GWAS/ldsc/rg_pdff_compatible_20261014T140601Z_22236252/pairwise_estimates.tsv`;
  these are not adjusted for the original correlation-matrix hypothesis family.
  The remaining 27 non-PDFF allele-compatibility failures were also repaired
  by job 22237000 and checked by 22237084: all 66 fixed-family pairs now have
  finite estimates and positive native observed-scale heritability. The terminal
  candidate preserves successful historical fits and labels every source;
  repaired and historical preprocessing are distinguished, and p-values are
  nominal. Six unconstrained estimates above one involve weakly heritable
  Ghodsian NAFLD or FinnGen cirrhosis; they retain their SE rather than being
  clipped. Source: `GWAS/ldsc/rg_missing_compatible_20261014T150911Z_22237000/rg_12trait_66pair_terminal.tsv`.
  H27 rows for L107 (RABGEF1; CKD, HF, eGFR) used the reversed labels, so their
  direction values flip on rerun (eGFR to concordant harm) [Inference, not
  rerun]; L107 is on a published list, so the direction table for new calls is
  unchanged. Not traced: the AlphaGenome atlas inputs and any output column that
  carries a Pazoki beta.
- *Refit candidate README.* Arm sizes now read 138 / 636 (774 units), 136 / 638
  and 771 participants (137 / 634), matching `arm_deg_counts.tsv`. The writer
  scripts `scripts/analysis/alt_contig_recount/r5_refit_corrected_controls_alt.R`
  (lines 21 to 25 and 720 to 723) and
  `scripts/manuscript/resource_f_five/refit_corrected_controls.R` (lines 14 to 18
  and 773 to 776) still print the 844-era arm sizes.

**Gene Catalog r8 candidate, 2026-10-09 (historical build record; local review
completed on October 14 below):
`RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-r8-20261009T142027Z`.**
Built on the 774-unit refit, adopted by the user on 2026-10-11
(`primary_source_controls/deg_results.csv`
of the refit reproduction, 24,115 rows), and on COLOC r2, adopted 2026-10-14 after r8 was built. The user approved the
PASS-00 selection in session (decision `2026-10-09-gene-catalog-r8-selection`;
the signature reads 2026-10-10T02:17:24Z because the UTC date had rolled over).
PASS-00 and PASS-02 to PASS-05 pass. 28,233 genes, 56,466 evidence rows, 117
programs with 6,633 memberships; no negative-call rule is defined, so tested-negative is not assessed. Primary classes:
genetically anchored 400, concordant (both supported; direction not compared) 26, established-state associated 1,129,
unresolved 13,127, untested 13,551 (the 2026-09-25 candidate had 383, 23, 1,054,
12,700 and 13,857 over 28,017 genes). Identity adjudication: 28,292 source rows,
28,233 genes, 59 quarantined rows, 37 irreducible multi-ID, 22 absent from
GENCODE v49, 0 identity conflicts. The program, Hotspot and MYOJIN inputs are
the 2026-08-07 candidate's and were not rerun; the GEN_CLOSURE and MYOJIN_CLASS
limitation texts are marked stale. On 2026-10-10 the user chose r8 as the
candidate that receives the PASS-06 attestation. The five images
(`portal_candidate/review/`: overview, MTARC1, HKDC1, GLP1R, THRB) were checked
on 2026-10-11 as recorded below. The portal review was pending at that date;
it did not prevent use of r8 in the manuscript under that user decision. The
refit r8 is built on was adopted on 2026-10-11; its COLOC r2 input was selected
on October 14 after independent completion checks.

**Gene Catalog review and sidecar completed, 2026-10-14.** The
[reviewed copy](../figures/candidates/catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8/PASS06_VALIDATED)
passed the completed local review and read-only verification in job `22236566`.
The actual reviewer is Codex remaining_catalog scientific review; this is not
a fabricated user signature. All 28,233 identities, 56,466 gene evidence rows,
117 programs and 6,633 memberships preserve their exact source values. Classes
remain 400 / 26 / 1,129 / 13,127 / 13,551. The historical October 9 bundle was
not changed, and no external portal was published.

Job `22236527` built the [current local interface](../figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/access_v2/index.html)
and [complete molecular-object companion](../figures/candidates/catalog-r8-completion-20261014T141645Z_22236527/molecular_objects/index.html).
All five rebuilt images were visually reviewed, all 20 source-native spot
checks pass, and 55 recorded inputs remain unchanged. The interface explicitly
uses 774 conservative genetic-identity clusters, not verified independent
donors, and records the 855 held annotations as untestable. Historical Myojin
class tests remain historical; no program or assay record supplies a gene-level
negative. The [Catalog results](RESULTS.md#catalog-review-and-molecular-object-companion-2026-10-14)
distinguish the 14,720 symbol-grouped, 14,704 native stable-ID and 14,682
Catalog joint universes. The earlier 24/8 universe split is not reproduced;
the explicit native source-symbol comparison gives 16/0 and the
Catalog-versus-native-ID comparison gives 6/28.

**Decided by the user, 2026-10-10:** r8 receives the PASS-06 attestation; the
remaining Pazoki copies and the `gwas_af_sources.tsv` hash are corrected (data
corrections above); N2a stays untestable; H46 stays dropped.

**Decided by the user, 2026-10-11:** the 774-unit five-cohort refit
(`bulk-corrected-controls-20261008T162308Z`) is adopted as the bulk analysis of
the paper.

**Work of 2026-10-11 (all outputs are candidates; no panel is installed in the
manuscript figure folders).**
*Refit re-derivation*
(`Analysis/Hypotheses/refit-adoption-review-20261011T151011Z/`): the 774-unit,
738-unit, sensitivity-count and leave-one-cohort-out numbers reproduce; values
are in the H1 refit and H1-LOCO rows.
*Continuum remainder* (molecular 21, 30, 40, 50, 60 to 65, systems 67 to 72, 74;
core 06 to 09): roots
`RNA-seq/results/histology_anchored_continuum/candidates/hac-continuum-corrected774-20261011T152106Z`
and `.../molecular_layers/hac-molecular-layers-corrected774-20261011T152106Z`;
panels `figures/candidates/continuum-774-panels-20261011T152106Z/`,
display-corrected panels in `continuum-774-panels-20261011T165944Z/` and the
compact 4H in `continuum-774-panels-20261011T170854Z/`. The same 18 systems pass
strict maxT; 0 of 43 exceed the connected-subgraph null and 1 of 43 (system_01)
the size-matched null; stage and continuum system effects have Spearman rho
0.73 with 37 of 43 the same sign (0.63 and 32 of 43 on the 367 disease units
once the 68 GSE213621 controls coded as stage 0 are removed; review round 3,
R05-02) (`fig4h_competitive_outline/provenance/audit.tsv`);
7 of the 8 NMF axes with a stage estimate are stage-supported (k4-P2 and k6-P6
have none); the paired step 40 uses 54 GSE193066 participants whose identity
calls are unresolved.
*Composition* (`figures/candidates/composition-774-20261011T151616Z/`): the
five pooled cohorts are rebuilt on the 774 units; the other four cohorts (413
records) are unchanged and not identity-resolved; 1,187 records in all. The
method criteria pass (no zero inflation, median hepatocyte fraction 0.887,
direction correct in 3 lineages, 1 lineage eligible as a covariate,
Macrophages). Mean hepatocyte fraction is 0.886, 0.885, 0.867, 0.851 and 0.825
for F0 to F4 (103, 102, 88, 69 and 25 units) and 0.903, 0.892, 0.880 and 0.859
for NAS 0, 1-2, 3-4 and 5-8; Spearman with ordered stage is -0.351 for fibrosis
(387 units, 3 cohorts) and -0.419 for NAS (395 units), BH p below 1e-11 for
both over the two tests. The composition-adjusted program models were rebuilt
as candidates on 2026-10-12 (Figure 4E, below), as were the continuum
composition-adjusted coefficients (item 5 below).
*Figure 3E and S3* (`figures/candidates/fig3e-s3-774-20261011T165944Z/`; S3C in
`figs3c-774-20261011T165944Z/`): calls at BH q < 0.05 versus F0 are 1,338,
4,266, 6,190 and 7,043 for F1 to F4 (387 units in GSE130970, GSE135251 and GSE162694; the F1 count is 14 without GSE130970), and versus NAS0 5,338, 9,337 and 9,841 for
NAS 1-2, 3-4 and 5-8 (NAS0 is 45 source controls and 1 disease unit, so these are mainly disease-versus-control counts). P6 is not shown (unavailable on the corrected counts).
The NAS NMF lines are relative to NAS 1-2; a refit with NAS0 as the reference
has not been run.
*Gene Catalog r8* (`Analysis/Hypotheses/r8-review-20261011T150915Z/`): the
five review images were checked against the tables; the stored spot check
matches 20 of 20 values (`spot_check_20.tsv`). 28 of the 41 Myojin class rows
test the earlier class sizes and are labelled stale
(`passport_assay_status.tsv`). On 2026-10-11 the user decided that r8 is the
Gene Catalog used by the paper on the strength of these checks; the PASS-06
attestation and seal are a portal build step, change no number or class, and
are not required for the manuscript. The NAS NMF activity model was refitted
with NAS0 as reference (46 NAS0 units, all with scores; candidate
`figures/candidates/fig3e-s3-774-nas0ref-20261012T003338Z/`).

**Historical installation, 2026-10-12 01:06 UTC, superseded October 14** (`logs/display_fixes_774/install_774_panels.sh`;
previous PDFs and record in `figures/candidates/pre-774-install-20261012T010624Z/`):
Figure 3D, 3E (fibrosis direct-label layout), 4G, 4H (18 strict maxT outlines),
S3C, S3G and the S3 endpoint concordance panel, all on the 774 units. Not
installed: Figure 3G (candidate is 176 x 176 pt, working panel 244 x 200 pt), the
S4Q to S4W and renderer-62 S3/S4 panels and the full S3 stage figure (no
same-named working file; candidates in `continuum-774-panels-20261011T152106Z/`,
`continuum-774-panels-20261011T165944Z/` and
`fig3e-s3-774-nas0ref-20261012T003338Z/`). Assembled Figure 3 and 4 PDFs were
not rebuilt.

**Historical build state, 2026-10-12; panels installed October 14** (values in
`docs/RESULTS.md`):
(1) Figure 3F on the 774 stage table,
`figures/candidates/fig3f-no-pip-774-20261012T011128Z/`: 54 genes, 216 effects,
59 at BH FDR < 0.05 (61, 244 and 65 on the 844 records); GNMT F1 to F2 is q
0.0562 on 774.
(2) Figure 3G in the working layout,
`figures/candidates/fig3g-working-layout-774-20261012T011006Z/`: 3,776 points
(Neither 3,156, Transcriptomics 194, Genetics 400, Both 26; 3,622 points
before); a legend-placement fix is being rendered separately.
(3) Figure 6 from Gene Catalog r8,
`figures/candidates/fig6-catalog-r8-20261012T011202Z/`: 6B and 6D over 28,233
genes (13,551 untested, 13,127 unresolved, 1,129 established-state associated,
400 genetically anchored, 26 concordant, meaning both supported with direction
not compared); the working 6B and 6D were on the
30,754-gene catalog of 2026-08-07; wording fixes on 6B, 6C and 6E are being
rendered separately.
(4) Figure 4E program models,
`figures/candidates/fig4e-program-axis-774-20261012T010830Z/`: working roster
469 records, 113 eligible, 27 fibrosis and 4 NAS calls, 0 and 0 after
composition adjustment; 774 units only (three cohorts, 364 units, 112 eligible;
Figure 4E after review round 3) 0 and 0 calls, also after composition
adjustment with either fraction set; 774 units plus GSE174478 (457 records, 112
eligible; sensitivity) 16 and 1 calls, 3 and 0 after composition adjustment
with the rebuilt fractions and 1 and 0 with the earlier fractions (H85), and 0
and 0 under the truncated Knapp-Hartung variance (R12-01).
(5) Continuum composition-adjusted coefficients,
`figures/candidates/continuum-composition-774-20261012T012137Z/`: 367 units
(109 GSE162694 + 258 GSE213621; 402 before). Stromal ECM 0.760 [0.654, 0.866],
composition-adjusted 0.294 [0.154, 0.434], BH q 1.4e-4 in arm A with the rebuilt
fractions (0.203 to 0.479 across arms and fraction sets; 0.726 and 0.488 on the
844 records); Ductular injury 0.707 and 0.601 (0.60 to 0.73 across arms; 0.677
and 0.639 on the 844 records). Adjusted Stromal ECM in GSE162694 alone is 0.172
[-0.032, 0.377], BH q 0.233, in arm A with the rebuilt fractions (0.097 to 0.217
across arms); with the earlier fractions it is 0.344 to 0.433. The change comes from the rebuilt cell-type fractions: with
the 774 scores and the earlier fractions adjusted Stromal ECM is 0.479, so this
estimate depends on the deconvolution input. Calls: programs 91 to 94
unadjusted, 67 to 67 adjusted (8 gained, 8 lost); systems 38 to 37 unadjusted,
11 to 24 adjusted (14 gained, 1 lost).

**Previously open, completed 2026-10-14:** the round-three installation,
corrected S3/S4 panels, selected COLOC r2, Figure 4A schematic, reviewed r8
Catalog and molecular-object companion, and eight PDFF correlations.
The missing correlations failed in installed LDSC's allele-alignment path;
explicit filtering using its native allele rule produced eight checked
estimates without changing original inputs. Figure 3/4 assemblies use native
panel sizes and have passed source, font and visual checks. The final
publication checklist remains separate from these completed corrections.
On 2026-10-13 the user told the main session to proceed with its
recommendations: Figure 4E is the four-cohort arm (774 units plus GSE174478,
named in the legend), and Figure 4H carries a single outline on the one system
exceeding the size-matched null, with no maxT count printed (`docs/PAPER.md`
now says this in both passages). The other 4E and 4H versions are kept in the
same directory. Review round 3 replaced both choices. Figure 4E uses the
three-cohort 364-unit arm (0 and 0 calls), because all 20 four-cohort calls
need GSE174478 and the untruncated Knapp-Hartung variance (R13-03, R12-01).
Figures 4G and 4H use the 367 disease units, with no outline in 4H, because
the GSE213621 stage code placed its 68 controls at stage 0 (R05-02). The Figure 3 and Figure 4
candidates are in `figures/candidates/fig3-fig4-layout-774-20261013T135430Z/`:
`Figure3.pdf` (468 x 648 pt) and
`Figure4_four_cohort_4E_4H_single_outline.pdf` (468 x 576 pt). Figure 3A there
shows the five cohorts (326, 207, 140, 71, 30; 774 units) with Sex, NAS and
Fibrosis availability tiles and a bar axis labelled "Units"; the Control tile
column is not drawn because it read "Partial" in every row. Figure 3B has a
grid total of 364 (units with both an exact fibrosis stage and a NAS:
GSE135251 207, GSE162694 86, GSE130970 71).
Historical assembly discrepancies, superseded by the 2026-10-14 installations:
Figure 3C at the working path prints "Unique to study (9,777)" and
five cohorts; the 774-unit panel has four cohorts and 8,065 genes unique to one
cohort, 1,958 in two, 502 in three and 105 in four
(`figures/candidates/bulk-release-corrected-controls-20261008T162308Z/figure3/panels/fig3c_cohort_alluvial.pdf`;
compact version in the layout directory). Figure 3A at the working path prints
"9 cohorts, 1,284 samples" and Figure 3B a grid total of 630; the 774-unit
versions show the five cohorts (774 units) and a grid total of 364 (same
bulk-release candidate; compact versions are in the layout directory). Figure 4F at
the working path predates the 774 refit and clips its "90%" tick; the 774
version at slot size is in the layout directory. In 4F the GSE126848 Stromal ECM
line departs from the other four cohorts (that cohort is 23 control and 7
disease units, all male). Figure 3 at native panel sizes was 468 x 1,339 pt, so
C, E, F and G were re-rendered smaller with 6 pt text; no panel is scaled.
Examined and closed, kept as a record: the rise from 11 to 24 composition-adjusted
systems on the 774 units was examined on 2026-10-12
(`figures/candidates/continuum-composition-774-diagnosis-20261012T013141Z/`):
it follows from the rebuilt cell-type fractions (11, 10 and 24 systems for the
844 records with earlier fractions, the 774 units with earlier fractions and the
774 units with rebuilt fractions; adjusted Stromal ECM 0.488, 0.479 and 0.294)
and no defect was found. Reporting rule agreed by the user on 2026-10-12:
composition-adjusted continuum results are reported as dependent on the
deconvolution input and no individual adjusted system or program is named as a
finding. The marker-expression check, the attribution to cell types and the
alternative composition codings were done on 2026-10-12
(`figures/candidates/continuum-fraction-markers-substitution-20261012T014841Z/`,
`figures/candidates/continuum-rare-type-arms-20261012T020132Z/`): neither
fraction set agrees better with marker expression overall, no single cell type
and no alternative coding removes the difference between the sets, Ductular
injury holds after adjustment in every estimable coding (0.60 to 0.73;
unadjusted 0.707), and adjusted Stromal ECM is called at BH q < 0.05 in every
estimable coding but runs from 0.203 to 0.479 (unadjusted 0.760). Installing
these candidates falls under item 3. The S4G panels with the GSE213621
Control marker and whisker in gray are done as candidates, not installed, in
`figures/candidates/continuum-fraction-markers-substitution-20261012T014841Z/basis_774_s4g_control_gray_whisker/continuum/`
(panel text and stage table identical to the previous render). The V1
glucocorticoid-covariate refit, historical N3 ambient-RNA check (2026-10-13)
and LDSC genetic-correlation rerun are done. The N3 source-diagnosis refit and
independent checks are also complete; no corrected primary calls are detected
(V1 and N3 rows, Pazoki bullet above).
Not resolved: the V2 hepatitis-screen row and notes were
removed on 2026-10-10 by an edit outside this work; the V2 tables under
`Analysis/Hypotheses/V2-20261007T210714Z/` are intact.

**Historical scientific review corrections, 2026-09-23.** The Figure 2
uncertainty claims below are superseded by the coordinate correction at the
top of this document. Figure 5B/5C use historical 844-record inputs, and the
September Catalog is superseded by r8 selected on 2026-10-11. The corrections to the
2026-09-22 scientific review are recorded in the
[outcomes document](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md).
Items 3 and 4 (Figure 2 coding-mass inference, reported as physical-cluster
bootstrap intervals with no p-value), item 6 (the coding explanation for sparse
overlap, removed from the text), items 7 and 8 (Figure 5C, rendered in the
Figure 5 v7 candidate) and items 12 to 14 (Figure 7 wording) are complete as
candidates. The GSE130970 control relabel and bulk refit (item 1) is complete:
844 participants, 140 controls and 704 disease, 1,725 canonical DEGs (Ensembl
gene rows) and 2,036 interval-null DEGs. The COLOC signal-eligibility and trait-driver
rerun (items 2 and 5) is complete: 436 tier-1/2 SuSiE genes, 996 in the union
and 708 in the full portfolio. Figure 2 A–I, Figure 3 bulk panels and 3G, the F0
arms, Figure 4E, the continuum transcript layer, Figure 5 v7 and the evidence
classes are rebuilt on these inputs as candidates; the documents now quote the
corrected values. The Figure 3 lncRNA inset is rebuilt as a candidate
(`fig3-lncrna-corrected-controls-20260925T195928Z`: 358 high-confidence lncRNAs;
matched direction OR 0.210 over 499 pairs, with an outcome-blind 0.10-SD
caliper). The Figure 4G–H 43-system values stand: that analysis reads no control
label or disease-versus-control result and uses the adopted 23,370-gene
expression-filtered universe, which the caption states. On 2026-09-25 the user
approved adoption of the corrected Gene Catalog release
(`RNA-seq/results/evidence_passports/candidates/gene-catalog-corrected-20260925T175215Z/`).
The build passes PASS-00 to PASS-05: 28,017 genes, 56,034 evidence rows, 117
programs unchanged, zero tested-negative rows, and primary classes genetically
anchored 383, concordant 23, established-state associated 1,054, unresolved
12,700 and untested 13,857. The PASS-06 manual review attestation and seal are
pending with the user at that date; this September Catalog is superseded by
r8 and its completed October 14 local review. Figures 6B and 6D were to
re-render after that review. No other
corrected figure or release is adopted. Items 9 to 11 (cross-assay module candidate) are
deferred by user decision, and that candidate stays out of the paper until it
is corrected.

**H&E-to-program external test, 2026-10-03: run once; primary endpoint not met; not adopted.**
The GTEx-fitted image model was applied without refitting to 50 TCGA-LIHC donors
with a non-tumor liver cryosection and RNA-seq from the same sample. Median partial
Spearman over 14 testable programs is 0.204 (within-site shuffle p = 0.103); one of
16 programs (Myeloid / immune, DOCK2) meets the rule in both the primary and the
negative-control-residualized versions, and the hypoxia negative control (0.313)
exceeds the family median. The claim stays at "readable from GTEx postmortem H&E
beyond recorded covariates". Numbers and limits are in
[RESULTS](RESULTS.md#he-to-molecular-programs-external-test-in-non-postmortem-liver-2026-10-03-not-adopted).

**Mouse diet models against human programs, 2026-10-03: candidate, not adopted.**
Seventeen public diet-versus-lean-control contrasts (15 GEO series) were scored on the
fixed programs through one-to-one orthologs and compared with the human fibrosis-stage
trend (four cohorts, leave-one-cohort-out Spearman 0.80 to 0.83). The prespecified
primary endpoint was not met: series of the same diet class are not more alike than
series of different classes (difference −0.012, permutation p=0.479). Median agreement
with the human trend is 0.726; 2 of 16 testable contrasts pass BH q<0.05. 27 of 62
stage-associated human programs are untestable in mouse. This uses mouse data and so
sits outside the human-only Resource scope unless the user adopts it. Numbers are in
RESULTS; output in `Analysis/MASLD_Model_Benchmark/executions/preclinical-fidelity-20261003T153626Z/`

**Two-model Figure 7 candidate, 2026-09-22: development only.** The
99-participant GSE267145 held-participant RNA-to-H3K27ac comparison now includes
a training-fold-tuned histology, sex and RNA-technical comparator. Its local RNA
model has 0.14133 skill, versus 0.01971 for that comparator, on the same target.
An exact-histology/sex comparison recovers differences in 22 disjoint pairs
(local-model difference skill 0.11070), limited to fibrosis 0–1. At the primary
20 additional measurements out of 99, beyond 45 already measured, RNA diversity and model disagreement find more
candidate-pool error than random, but neither has a Holm-significant advantage
in held-participant refit MSE. The single sequence-by-RNA residual experiment
missed its +0.01 complexity margin (interaction − additive original-target
skill −0.0000159),
so the simpler chromatin form remains selected. The allele comparison completed
on 2026-09-23: on 24,430 matched leads, native-1-Mb MSE is 0.300161 and
combined MSE is 0.290705. The 3.15% reduction misses the fixed 5% complexity
margin despite positive paired contrasts, so calibrated native 1 Mb remains
the selected simple model. Fresh sequence inference remains experimental
because of between-process numerical variation; the exact-table query is
reproducible. All six candidate panels are rendered and independently checked
at `figures/candidates/figure7_two_models_20260922_21849207/`; their numerical
reconstruction passed. They remain separate from the working figure. An inspected
17-donor paired public source is not yet eligible for deployed H3K27ac validation;
Wenz is not yet eligible for allele-effect transfer. Exact numbers and limits are
in [RESULTS](RESULTS.md#two-model-molecular-value-development-2026-09-22-not-adopted).
No model or figure is adopted from this candidate.

**Regulatory-model initial week, 2026-09-15: implementation in progress.**
The user authorized the bounded campaign in
[research program §19.4](technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md#194-authorized-initial-week-2026-09-15).
The user removed campaign compute-hour and concurrent-GPU caps on 2026-09-16;
the 500-GB file allowance, scheduler constraints and scientific holdouts remain.
Endpoint tables, scalar-effect training, training-calibrated donor comparisons,
P4 life-stage/completeness accounting and a local Catalog specimen are implemented
as candidates. Reviewed initial reporter and donor comparisons are recorded in
[RESULTS.md](RESULTS.md#initial-regulatory-context-comparisons-2026-09-15-not-adopted):
no interaction fit beats zero by point RMSE, and the donor pilot establishes
no held-region improvement. The checked six-entry Catalog extension contains
436,838 evidence records, adding measured Currin associations and all twelve
matched out-of-fold prediction recipes while preserving the initial specimen.
Twelve fixed-representation recipes and twelve adaptation
arms are complete; all fitted reporter interaction arms remain worse than zero
by point RMSE on 4,358 source constructs. The 6,834-variant adaptation pilot
supports six adapter gains over their matching fixed controls under paired
uncertainty, with only five validation chromosomes and one seed. It does not
establish superiority to native AlphaGenome. Full native attempt accounting
is complete: 30,611 scored, 1,711 excluded, zero unresolved among 32,322 eligible
rows. The repaired matched comparison and precision strata are complete; independent
numerical reconstruction remains pending. Its preliminary shared-validation table
puts calibrated native ATAC ahead of the best adapter by RMSE. The independently
reviewed zero-adapter output comparison passed with exact agreement on 32 development sequences;
the separate archived/new reporter-feature difference remains unresolved.
The separate measured ATAC check is complete: original and ten adapted
backbones were assessed on 512 held regions after calibration on 512 training
regions, using the 138-donor Currin count matrix. The strongest allele-effect
adapter's ATAC error difference remains uncertain; both partial-update pilots
increase ATAC error. This is regional same-source assessment, not broad
native-function preservation. The stronger 99-donor comparison and twenty
uncertainty-aware allele-effect fits are complete and independently checked.
RNA improves held-donor predictions at trained regions; improvement across
held regions remains unestablished. Among the loss-training arms, the lowest
observed allele-effect RMSE remains the prior squared-error adapter. These are development results with
fixed-prediction uncertainty, not external confirmation or adopted findings.
Removing 16 exact extra copies from the historical signed-splice table leaves
621 scored pairs in 247 blocks: concordance 0.5765, interval 0.5383–0.6135.
Independent reconstruction preserves the conclusion of weak direction below
the original ranking criteria; the exact-event inference remains unfinished.
Protected outcomes remain closed. GSE267145/GSE296875
derivative weights are releasable under the adopted owner decision, and
GSE267145 counted-interval sequence use carries the recorded 1-bp uncertainty.
The 117 programs and adopted genetic memberships remain unchanged.

**Continuum translation extension, 2026-09-06: implemented candidate.** The
complete Gene Catalog sidecar, integrated 117-program figure, signed genetic
records, composition sensitivity and matched liver/plasma protein analysis are
available. This extension used the 2026-08-17 release (462 SuSiE and 822 ABF
genes); the corrected COLOC candidate has 436 SuSiE and 822 ABF genes (see the
[2026-09-23 scientific review block](#what-is-official-updated-or-open)).
Only 83 of 925 supported main-trait signal pairs resolve direction under the
declared conservative marginal-effect rules; multiple signals and ambiguous
alleles remain unresolved. Of 117 programs, 93 meet the existing stage/sex
two-axis support rule, 68 meet that rule with composition adjustment, and 26
also pass BH separately in each evaluation cohort and axis after composition
adjustment. These counts are different criteria, not independent replications.
The public protein data link 41 participants; 49 of 345 shared proteins pass
the adjusted matched-association BH test. Genetic cross-study direction
conflicts remain explicit. Details and completed checks are in the
[technical review](technical/continuum_translation_20260906T135538Z/IMPLEMENTATION.md).
The candidate is not an adopted figure or a treatment-response claim.

| Area | Current state | What may be used now | Next dependency |
|---|---|---|---|
| Paper scope | **Settled** | Resource/framework framing and seven-figure story | Keep all drafts inside [PAPER.md](PAPER.md) |
| Bulk disease-vs-control | **774-unit refit adopted 2026-10-11; current Figure 3 panels installed 2026-10-14** | The bulk analysis is 774 genetic-identity clusters (not verified donors) across five cohorts: 138 controls and 636 disease, 24,115 genes tested, 1,854 canonical DEGs (1,333 up, 521 down), 2,131 interval-null. GSE213621 supplies 326 units, 68 controls and 45.1% of the quality weight; without it 1,236 of 1,854 calls are kept (66.7%), against 88.6% to 96.8% for the other four removals. Main Figure 3C shows all five per-cohort fits: GSE126848's identity-selected male-only roster (23 control, 7 disease) is modeled with constant sex omitted, independently verified against native counts. Figure 3C uses its within-cohort filterByExpr refit in place of its 6,147 fixed-family calls (3,440 calls, BH over 15,592 genes; user decision D3, 2026-10-21), and the refit joins a 12,180-gene per-cohort union (8,341 called in one cohort and 3,839 in two or more; 1,211 of the 3,839 have opposite signs between cohorts), distinct from pooled membership. Source and exclusion counts are recorded in the H1 row above. The four-cohort 392-unit endpoint analysis is sensitivity evidence only, recovering 1,755/1,854 canonical DEGs with all-gene log2FC rho 0.953. Historical: 844-record corrected-controls fit 140 / 704, 1,725 DEGs, endpoint 416, 1,630/1,725, rho 0.940; 2026-08-10 fit 157 / 687, 1,347; 418, 1,307/1,347, 0.938 | Current panels installed and visually checked; after review round 3, Figure 4E uses the three-cohort 364-unit arm of the 774-unit substrate; the four-cohort 457-record arm, which adds GSE174478 from outside the pooled 774-unit fit, is a sensitivity |
| Noncoding RNA source check | **Passed under a documented amended cohort-level rule** | All five cohorts retain every library; the stricter all-sample rule remains preserved as a failed sensitivity | Use only the passed, hashed source-check candidate for lncRNA outcome work |
| lncRNA identity and program content | **Validated candidates; not adopted** | Versioned GENCODE v49 identities, seven genomic classes, mapping exclusions, unchanged 117-program content, exact pooled biotype subdivisions, 271 high-confidence disease-associated lncRNAs, and 48/49 passing native leave-all-lncRNA-out sensitivities | Integrate these candidates with composition and the synchronized release |
| Stage analysis | **Five-cohort synchronized candidate validated** | Clean job `19966029` produced checksum-verified adjacent fibrosis, adjacent NAS, F-versus-F0, and per-cohort models with complete 23,370-gene BH families from the validated 844-participant DGE. The corrected-controls stage rerun uses 23,702 Ensembl gene rows; F0→F1 gives 1,426 genes with all F0, 929 with disease-only F0 and 975 with a control-status term, and the earlier control-confounding explanation is withdrawn. On the 774 units (`bulk-release-corrected-controls-20261008T162308Z`, `f0-reference-arms-r8-20261008T193309Z`) the families are 24,115 genes and F0→F1 gives 1,338 genes with all F0 (n=103), 783 with disease-only F0 (n=89) and 814 with a control-status term; GNMT F1→F2 has BH q 0.0562 and is no longer below 0.05. The all-F0 F0→F1 count rests mainly on GSE130970 (14 of 1,338 calls without it, 155 units), and GNMT and MAT1A F0→F1 have q 0.207 and 0.0892 without GSE130970. Candidate v1 remains failed execution provenance after job `19957749` was edited during execution | Corrected 774-unit stage panels and sources installed and visually checked 2026-10-14; historical v1/v2 records remain preserved |
| Histology-anchored molecular continuum | **Corrected 774-unit main and supplementary panels installed 2026-10-14** | Primary stage/sex-adjusted inference uses 435 selected units (109 GSE162694, 326 GSE213621); descriptive windows include all five cohorts. Fixed 43 systems: 33 collection-balanced BH/direction continuum calls, 31 feature-balanced calls, 20 all-four permutation-BH calls, 18 strict maxT calls; joint cross-cohort maxT 33. Competitive specificity is separate: 0 connected, 1 uniform. The 326 GSE213621 units include 68 controls coded as stage 0. After review round 3 (R05-02), main 4G and 4H use the 367 disease units; 4H has no outline and since 2026-10-17 plots the donor-first system estimates, rho 0.670 and 34/43 same direction (feature-first 0.631 and 32/43; 0.7265 and 37/43 with the controls). Complete-source and window displays remain descriptive; stage k4-P2/k6-P6 unavailable, EMT/RELA/SMAD3 context untestable. | Preserve fixed memberships, separate families and qualified units; no adjusted individual system nominated |
| Noncoding RNA figure package | **Nine-panel candidate independently validated and visually reviewed** | Revised Figure 3 pooled/stage panels, six lncRNA/program supplementary panels, and the Figure 6 molecular-object routing panel under `noncoding-resource-2026-08-12/plot-upgrade-v3`; final repository validation job `19779385` passed | Assemble only after synchronized adoption |
| Noncoding DNA figure package | **Corrected Figure 2D and unassembled noncoding annotation companion adopted and installed 2026-10-13; broader package remains candidate material** | The corrected locus-balanced consequence panel replaces the August 17 Figure 2D. The companion is installed outside the assembled A–I set; its unresolved fractions (52.11% and 60.43%) reflect failed lookup and were withdrawn by H58 (repaired value 0%). Both use `fig2-coordinate-panels-20261013T220927Z_22219710`. Other promoter/ABC/open-chromatin panels, six-lineage accessibility, the 30-variant heatmap and the full noncoding-mass distribution remain candidates under `noncoding-resource-2026-08-12/noncoding-dna-v5` | The adopted gene-level target companion is installed separately; retain the broader DNA package as supplementary candidate material until synchronized selection |
| Cell composition | **Synchronized candidate rebuilt; render validation active** | Donor-aware hepatocyte and non-parenchymal summaries use the accepted composition source and the same stage sample manifests as Figure S3 | Validate candidate v3 and keep canonical output unchanged |
| Hotspot programs | **Diagnosis-qualified refits and corrected 4C/S4C installed 2026-10-14** | All 117 program memberships and weights remain fixed. Original-source correction excludes two cirrhotic GSE174748 explants under the existing rule: H59 primary fit (60-donor universe, 55 donors fitted in GSE202379 and GSE244832; per-dataset variance, permutation BH over 117) selects zero programs (Ductular injury q=0.323375, ECM q=0.577711); the 62-donor equal-variance fit (Ductular injury/BICC1 q=0.037633, ECM q=0.085458, zero HC3-selected) is a sensitivity. The ambient common universe has maximum n=60; neither focal paired change passes the 117-program BH family. Five datasets have decontX output; GSE189600/Liver Atlas remain failed, 106 programs pass integer-count transport and 11 T-cell programs are untestable. Ten same-atlas donor-paired localization contrasts and the separate 26-donor content-lineage meta analysis are unchanged. | Preserve source exclusions and distinct testing families; do not describe localization as transcript origin or cell autonomy |
| Continual single-cell integration | **Ran to a terminal decision; primary adoption REJECTED on 2026-08-15; permitted only as supplementary sensitivity** | The corrected scvi-tools 1.3.3 replay-plus-EWC port ran end to end on real GPU fits (V100/L40S, seed 17, one all-lineage SCANVI plus five lineage models) against the read-only 275/104 descriptive and 273/102 analyzed census, the Harmony export built to the current requirements, and the fixed 1,232,318-cell x 4,000-HVG count object; the final regression suite passed 136 tests with 0 failures and 1 skip. `candidates/final-adoption-decision-v42-2026-08-15T215000Z/promotion_decision.json` is terminal: 13/14 criteria met pass and `independent_external_routing` FAILS on the untouched GSE212837 cohort (12 donors, 128,851 cells; donor-balanced macro-F1 0.2756 against the predefined 0.70 floor). The method identity is a **replay-plus-EWC-guided geometry-preserving adapter**, not the unmodified preprint method. `candidates/router-v2-final-decision-v46-2026-08-17T233700Z/decision.json` then rejects the pretrained-CellTypist Router V2: on the untouched GSE296875 cohort (39 donors, 68,398 cells) overall donor-balanced macro-F1 0.7071 clears the 0.70 floor and paired improvement over V40 is 0.2033 (95% CI 0.1847-0.2219), but B cells 0.2888 and NK-T 0.4229 fail the predefined 0.50 per-class floor. Harmony remains the primary representation, and the strict seven-donor reference is retained after both the local nine-donor arm and the common-universe 26-donor HLiCA Andrews_2022/Andrews_2024 expansion failed the reference-preservation criteria; wholesale HLiCA ingestion remains forbidden. The canonical H5AD, manuscript, paper authorities, and the byte-identical 117-program separation are unchanged. ⚠ the documented full-scale sbatch path never executed. B6K jobs `20208726` and `20215473` both sat pending from 2026-08-13 and were cancelled at 00:00:00 elapsed, so every fit ran instead inside `gpu_interactive` V100/L40S allocations and the persistent CPU allocation `19664223`. ⚠ `secondary_technical_stress` passed on technical stress alone with `independent_confirmation_pass: false` because GSE136103 is not independent and GSE174748 is underpowered. ⚠ reference fits use the prespecified 2,000-cell per donor-by-label cap, so 61,960 of 216,957 strict-reference cells entered training while embeddings still cover all 216,957 | Nothing on the critical path. Cite the candidate only as supplementary sensitivity and never as the Figure 4 representation. GSE212837 and GSE296875 are diagnostic-only and cannot tune any router or threshold. Any continuation is a new hierarchical immune-aware router, settled in advance, validated on a third untouched cohort, which is a new method study rather than a rescue of Router V2 |
| Recurring biological examples | **Genetic states checked on adopted r2; molecular mechanism remains provisional** | GNMT supplies an ABF-only direct-NAFLD locus-to-state bridge; MAT1A and CYP2C19 define the severity-linked metabolic state. Stromal ECM/Ductular injury separate disease association, ambient qualification, descriptive lineage localization and matched-null tissue organization; IGFBP7 and BICC1 are program member-gene provenance only | r2 retains GNMT direct-NAFLD ABF PP.H4 0.705832 and SuSiE 0.128219; none of the five has primary SuSiE support. MAT1A, CYP2C19, IGFBP7 and BICC1 have no finite SuSiE result and ABF maxima 0.070472, 0.037590, 0.082400 and 0.025806. Retain exact pair and unresolved-haplotype qualifications |
| Single-cell ATAC context | **ATAC v3 remains candidate; its r2 primary-pair context feeds adopted Figure 5D** | The isolated v3 candidate contains a 676,658-peak condition-blind atlas, exact donor fragment-record counts, matched DA, all 117 program projections in both cohorts and official chromVAR. Historical replay exported 5,362 signal pairs and 22,356,697 variant posterior rows over 816 gene-study pairs/462 genes. The installed October 14 r2 Figure 5D restricts this context to 767 eligible primary pairs over 436 genes: 28 direct MASLD/liver-fat and 739 enzyme pairs, across five shared lineages. Every retained pair meets the 0.95 mapped-posterior criterion | Peak overlap is regulatory-DNA context, not target assignment or functional validation. Per-variant posterior vectors retain the documented replay-tolerance limitation; adopting 5D does not adopt the entire ATAC candidate |
| Inherited genetic variation | **COLOC r2 independently checked and adopted 2026-10-14; Figure 1B/C, 2A/B/C/E and 5D updated** | Selected input `GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/`: primary tier-1/2 SuSiE 436, ABF 822, single-signal-only 560, union 996, full-portfolio SuSiE 708. All 50 studies and 948,503 unique study/chromosome/gene rows are present; all five recorded output hashes match. Shared-posterior correction reran 1,307 previously positive gene-study rows across 32 studies; other results and the other 18 studies were carried unchanged. Of these rows, 1,233 retained, 49 eligible at PP.H4 <=0.5 and 25 untestable; seven tier-1/2 genes remain SuSiE-untestable. Installed Figure 2C prints 14/21/70/104/524/263 since 2026-10-17 (coding split into own gene 14 and other gene 21; intronic and intergenic by direct transcript overlap, 450/337 before); E prints 340/64/23/9, with the 340 split into 127 SuSiE-tested in a non-EUR study and 213 not tested. F/G/H/I are render-identical and retain their artwork with r2 sources. Figure 2D retains the October 13 coordinate correction. Historical August membership was 462/1,013/745; it is superseded | The separate target-gene companion (424 protein-coding, 12 lncRNA) and native exact-example source rows are installed; separate supplementary selection and broader Resource adoption remain separate |
| Evidence classes | **Current r8 classes selected; exact Catalog copy reviewed 2026-10-14** | Native joint stable-ID comparison: 14,704 genes, 426 genetic, 1,155 RNA calls, 26 both; The repaired Figure 1B candidate (`figures/candidates/review-r3-fixes-20261015T112148Z/F1/repair/fig1b_complementary_maps_U_eval_repair.pdf`) uses the 3,769 genes evaluated by multi-signal COLOC (426 genetic, 220 RNA calls, 26 both; 21.21 expected within AveExpr deciles; expression not fully controlled, no direction stated; H57); the installed Figure 1B prints this universe (26 of 220 DEGs, 11.8%, against 9.6% expected within expression deciles; 935 further canonical DEGs had no multi-signal COLOC result); raw symbol-grouped comparison has 14,720 rows. Adjudicated Catalog jointly eligible identities 14,682; full classes 400 genetically anchored, 26 concordant (both supported; direction not compared), 1,129 state-associated, 13,127 unresolved and 13,551 untested. No jointly eligible quarantined identity has positive support. September 844-record 409/1,082/24 and 14,243-background values are historical. | Current selected bulk/r2 and reviewed r8 calls; distinguish missing testability from negative evidence |
| Spatial context | **Validated candidate awaiting synchronized release** | Complete 117-program observability, the two-program candidate family, label-blind HMSMA coverage, and categorical portal requirements | Keep HMSMA clinical inference withheld; rebuild genetics-dependent spatial classes from COLOC r2 (adopted 2026-10-14) |
| Figure 5 context | **Corrected working panels installed 2026-10-14** | 5B uses 774-unit RNA effects in the fixed July display universe: 3,054 measured, 211 protein-significant (188 concordant, 23 discordant), rho 0.457835, 66% rounded direction agreement. 5C retains 25 proteins with nine RNA and 23 protein calls; CD74 enters and HLA-DQA1 leaves the RNA set, both HLA-DQA1 markers open. 5D gives, per lineage and trait class, the percent of genes with at least 50% of the posterior in peaks of one or both sources: 22 direct and 424 enzyme genes, from 767 r2 gene-study pairs (28 direct, 739 enzyme); a locus-matched background companion uses the same genes. 5A is legend text since 2026-10-17. Installed 5E (2026-10-15) replaces the 13 rows chosen from 22 programs selected on the superseded July disease fit with the complete 117-program registry: 45 rows with a lineage-open COLOC r2 gene (31 programs, 3 with direct-trait RORA) or a robust protein-histology association (22 of 95 scorable, BH over 95; since 2026-10-17 members are counted once per primary-assembly gene in both columns), and no spatial column. 5F keeps the all-cell-adjusted calibration; its candidate source table moves Vu Ductular injury from `source_dependent` to `indeterminate`. September 844-record 5B/C and earlier v6 numbers are historical. | Preserve selection-conditioned protein and source-dependent spatial qualifications; never upgrade peak overlap to target assignment or the illustrative overlay to named histologic localization |
| Figure 3G genetics–state interface | **Corrected 774-unit working panel installed 2026-10-14** | Native RNA/SuSiE scatter contains 3,776 finite pairs: Neither 3,156, Transcriptomics sig. 194, Genetics sig. 400 and Both 26; Spearman rho of PP.H4 with \|log2FC\| 0.00141, nominal p 0.931 (signed rho 0.0217507 withdrawn as the wrong form for an unsigned predictor). The joint native stable-ID universe is 14,704 with 426 genetic and 1,155 RNA positives. Separate S3G interval-null sensitivity uses 426 expression-matched pairs (37 genetic-only, 28 control-only, 360 neither, one both; paired OR 1.3214, nominal exact McNemar p 0.3211). Historical 844-record and August values are superseded. | Unsigned PP.H4 versus signed RNA is descriptive; sparse overlap and failed significance do not establish biological opposition, complementarity or equivalence |
| Figure 7 model framework | **Validated candidate; working panels 7A–7F on disk, with 7C/F source exposure corrected October 14; not synchronized; no authority row existed before 2026-09-07** | Released severity instrument `release/masld-severity-v1.2/` (severity bitwise = v1.1; activity head added; joint-space refusal bound measured and NOT shipped): out-of-cohort Spearman 0.6631 (109 participants, ceiling 0.9523), two additional external cohorts, Spearman/ceiling 0.709 (263 GSE213621 identities), 0.866 (81 GSE276114 MASLD) and 0.680 (82 GSE276114 chronic viral hepatitis) as drawn in Figure 7D (record-level 0.677 to 0.794) (no demonstrated advantage over the recorded transductive native comparator; not equivalence), Historical source-reported ICC 0.7658 vs the Kleiner label's 0.1886 is identity-unqualified GSE193066 and out-of-contract scoring (only 3/58 pairs meet shipped bounds); it cannot support primary ICC or label-error calibration. Exposure: 109/109 GSE268273 samples are BulkFormer pretraining rows. RNA↔H3K27ac coupling skill 0.1221/0.1348 (rank 4 is a floor); late fusion +0.0030 [−0.0754, +0.0806]; DNA-LM signed direction at chance (0.5016 vs allele ridge 0.5033). Development lane 2026-09-07 (GSE267145, not adopted): a steatosis-selective chromatin component is stable across folds (mean \|cos\| 0.711) and RNA predicts it beyond the four grades, sex and library descriptors (partial 0.552 [0.351, 0.704], 0/500 null) — one cohort, age/BMI unmeasured, identity of the axis open | Owner review of the development lane; derivative-weight rights resolved by the 2026-09-07 owner decision; re-derive panels by someone who did not produce them; enter prose only after `docs/RESULTS.md` numbers are adopted in the synchronized release. **2026-09-22 development lane, not adopted** (`docs/RESULTS.md`, "Ceilings, matched baselines and the fold-0 adapter confirmation"): the closed fold 0 of the Currin caQTL set was spent on its single authorised confirmation and both prespecified contrasts cleared Holm (+0.03303 vs native ATAC at matched 2 kb, +0.02459 vs the matched fixed head); fold 0 is now spent. Five ceiling and matched-baseline measurements accompany it: the chromatin cis-specific increment is flat past 31 participants per fold against reduced-rank regression, which weakens the standing case that a second paired H3K27ac cohort would buy the claimed increment; AlphaGenome zero-shot already reaches 0.85 of the caQTL LD-tagging bound; the LX-2 interaction head was under-penalised and is null; severity sits at 0.82-0.95 of its label ceiling. A 7E exposure candidate moves the panel to 54 clean / 22 confounded / 8 unresolved. The remaining Liver_Atlas aliases were corrected to Guilliams GSE192742 on October 14; LangCell exposure is unknown, historical scBERT clean by chronology, mutable checkpoints unknown. Current model scores and gains are unchanged. The exact-checkpoint Figure 7C/F correction uses TranscriptFormer five clean/two training-seen studies; its pooled broad-class F1 includes all seven. Geneformer V2 alone has zero qualified clean studies among the displayed models. The five-clean-study companion uses 80 source donor identifiers and the original saved donor bootstrap, conditional on fixed models/studies. The historical truncated caQTL gain is absent from current framework artwork, not an uncorrected current-panel value |
| Figures 1/6 and spatial supplement | **Impact-layout candidate validated; not adopted** | Compact spatial assay-capability and Gene Catalog insets, the retained 117-program coverage heatmap, simplified HMSMA panel, and supplementary source-state matrix under `spatial-impact-figures-candidate-2026-08-11` | Do not add an adjusted coding-versus-lncRNA coverage panel until complete versioned assay eligibility exists; assemble the remaining panels after inputs are settled |
| Spatial robustness diagnostics | **Validated candidates; supplementary only** | Full source tables plus compact composition, unit/matched-null/weight, focused cell-context, and HMSMA mask panels; one bundle preserves the complete diagnostic provenance | Use the impact candidate for assembly and `spatial-diagnostic-bundle-candidate-2026-08-11` for complete provenance; do not change the two-program family or main Figure 5 calls |
| MASLD Gene Catalog | **r8 selected for the paper; local review and molecular-object companion completed 2026-10-14** | `gene-catalog-corrected-r8-20261009T142027Z` uses the adopted 774-unit bulk fit and COLOC r2: 28,233 genes, 56,466 evidence rows, 117 programs (6,633 memberships); classes genetically anchored 400, concordant (both supported; direction not compared) 26, established-state associated 1,129, unresolved 13,127, untested 13,551. The September 844-record Catalog is historical. Five rebuilt images and 20/20 independently rederived source values passed. The exact scientific tables are preserved in `catalog-r8-reviewed-20261014T142026Z_22236566/catalog_reviewed_r8`; the current local interface and 28,233 molecular-object records are under `catalog-r8-completion-20261014T141645Z_22236527` | Completed local review changes no gene call or class and publishes no external portal. Distinguish the 14,720 symbol-grouped, 14,704 native stable-ID and 14,682 Catalog joint universes. Synchronized Resource adoption remains a separate overall decision |
| Public perturbation plans | **Complete supplementary boundaries** | Valid results and nulls in supplement/provenance | No replacement-dataset search |
| Manuscript | **Working draft** | Seven-figure Resource story and bounded claims | Reconcile prose and legends to the synchronized release |
| Experimental follow-up screen | **Separate future paper** | Resource-to-screen prospective rationale only | No screen result may feed back into this paper |
| AlphaGenome–MASLD program | **Executed; candidate research layer, nothing adopted** | Nothing for the paper. Results are in [technical/ALPHAGENOME_PROGRAM_RESULTS_2026-09-14.md](technical/ALPHAGENOME_PROGRAM_RESULTS_2026-09-14.md) with a Gene Catalog attachment under `GWAS/finemapping/results/alphagenome_program/g-integration-20260915T072500Z/` in which every row is marked `candidate_not_adopted`. It changes no adopted number, gene class, evidence state, threshold, hero gene or figure | Reconcile remaining source-use and authority conflicts with current evidence. Local weights, local scoring and the reporter probe are now available; the original weight-access blocker is superseded. Broader training and independent evaluation remain proposed in the [integrated plan](technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md), alongside public-data expansion, six campaigns and the database specification |

The **2026-10-22 ATAC statistics correction** is available as candidate tables
in `Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-statistics-correction-20261022/`.
It corrects 195 program SE/p/q triples, 79,828 primary peak standard errors and
80,878 legacy sensitivity standard errors, including dependent numerical copies.
Effects, significance calls and evidence states are unchanged. The original
ATAC candidate, figures and official release remain unchanged; see the
[ATAC results](RESULTS.md#validated-single-cell-atac-context-v3-candidate-not-adopted).

## Active compute

**Historical disposition: Corrected rerun PROMOTED to canonical 2026-08-17.**
Its 462 primary genes and 1,013-gene union are superseded by the October 14
r2 adoption. The historical array finished and was adopted; its signal
eligibility and trait drivers are now under correction (see the
[2026-09-23 scientific review block](#what-is-official-updated-or-open)). Array `19648525` held
1,100 study-by-chromosome tasks and every one has reached a terminal state
(`sacct -j 19648525 -X`, read 2026-08-17):

`1,100/1,100  ██████████  1,057 COMPLETED · 42 FAILED · 1 TIMEOUT`

The 42 scheduler-failed tasks have structurally complete result files in the
read-only ATAC handoff audit (40 of them from a 35-minute window on 2026-08-11
when `06_susie_coloc.R` was transiently corrupted): 17 have no exceptional rows
and 23 contain only one anonymous, low-support `abf_only` row with no SuSiE
posterior, so none can enter the primary SuSiE replay. The read-only audit is
under `atac-context-v3-candidate-2026-08-11-r1/upstream_monitor/coloc-failure-audit-2026-08-12-r2/`.

Task `19648525_176` (2023_36280732_NAFLD_UKBB_EUR chr1) hit TIMEOUT after 3d18h
on 2026-08-14 and was re-executed as `20633389_176` (COMPLETED 2026-08-16, 3 min).
The canonical chr1 file is byte-identical to the rerun output (md5
aa913ed3b92c65d0c2724b7dd7dfc0dd), and COLOC r2 carries 2,035 chr1 rows for this
study (6 with a SuSiE posterior, none above 0.5).

Running as of 2026-08-17 10:20 EDT (`squeue -u $USER`):

| Job | Name | Resources | What it is |
|---|---|---|---|
| `20871418` | `cabA2bnull` | 8 CPU / 140 GB / cpu | Permutation null, ComBat arm; perm 5/20 at last read |
| `20871423` | `ambient` | 16 CPU / 200 GB / cpu | Stromal ECM program ambient sensitivity rescore (`dec_depth`); registry and membership hashes verified against READY, stored BH family reproduced exactly |

```bash
squeue -u "$USER" -o '%.10i %.9P %.24j %.2t %.10M %R'
sacct -j 19648525 -X --format=JobID%20,State,Elapsed,End%20
```

The adopted cache optimization changes only the `lambda_s_locus` diagnostic at
approximately floating-point precision; all tested PP values and reported
fields were bit-identical. The pre-cache products do not need to be rerun, but
the final release must record the two execution code versions. See the dated
technical note in the archived progress index.

## Figure readiness

| Figure | Status | Remaining work |
|---|---|---|
| **1. Multimodal MASLD Resource** | Resource architecture, two-map interface, seven-state vocabulary, spatial capability matrix, and experiment-routing grammar validated; coding-versus-lncRNA adjusted matching excluded by the source check | Lead with inherited variation, multicellular remodeling, and their cell/tissue context; retain coverage and evidence states as interpretive guardrails |
| **2. Inherited genetic variation** | Main A-I is canonical, with D superseded by the adopted and installed October 13 coordinate correction; its noncoding annotation companion is installed but remains outside the assembled A–I set. Early and late supplementary sets were independently rederived and validated as candidates; supplementary E, K, and N are explicitly withheld rather than rendered from unsupported adopted inputs | The separate gene-level target companion is installed; adjudicate supplementary adoption separately; the canonical RORA/FABP1 examples are EUR and support method-resolution evidence, not ancestry-specific biology |
| **3. Bulk remodeling and genetic overlap** | Corrected 774-unit panels and vector review assembly installed, including all five cohort contrasts and 54-gene/216-effect 3F; source values and biological-unit qualifications checked | No bulk/stage rerun remains pending; preserve selection and genetic-cluster limitations |
| **4. Single-cell programs and multicellular context** | Atlas schematic and corrected 4E–H/S4 panels installed; 4D compares all 117 discovery/content lineages, 4E uses 112 programs in the three-cohort 364-unit arm (the 457-record arm is a sensitivity), 4G and 4H use the 367 disease units with no outline in 4H (review round 3); vector review assembly complete | Unavailable NMF axes and untestable EMT/RELA/SMAD3 are explicit; no negative result inferred from missing observability |
| **5. Molecular/physical context** | Six working panels use corrected 774-unit 5B/C and r2 5D alongside unchanged assay-native 5A/E/F | Keep fixed July gene universe, selection-conditioned protein and source-dependent spatial states explicit |
| **6. MASLD Gene Catalog** | All five individual r8 panels and native-size vector review assembly complete; local review and molecular-object companion checked | Final publication layout remains; no external portal was published |

## What is still missing

The October 14 corrections listed above are complete. Remaining manuscript
work and unavailable evidence are distinct:

1. Assemble the complete publication set from the current individual panels,
   including the final layout of Figure 6's completed native-size review assembly, then complete the ROADMAP's single
   read-only output bundle and overall approval. Figure 3/4 review assemblies
   are already complete; older assemblies must not supply current numbers.
2. HMSMA donor, histology and MALDI-registration joins remain unavailable.
   Its coverage and label-blind spatial summaries stay descriptive; these
   missing joins do not require rerunning the completed analyses.
3. GSE193066 identities remain unqualified. Primary repeat-biopsy replication,
   ICC and measurement-error calibration remain excluded; the corrected N4
   null is untestable. Qualifying the cohort requires new identity evidence,
   not another resampling run.
4. Independent confirmation of the selection-conditioned protein example,
   additional paired RNA/H3K27ac donors and variant-combination experiments
   are future evidence, not unfinished fits needed to support this Resource.

Round-3 scientific review (40 findings verified) and prespecified analyses H51 to H60, run on 2026-10-15: nine completed (H52 untestable) and H56 stopped by its prespecified rule. Results are in RESULTS.md; the H52, H54, H57, H58 and H59 outcomes are summarized above: [RESOURCE_REVIEW_AND_HYPOTHESES_2026-10-15](technical/RESOURCE_REVIEW_AND_HYPOTHESES_2026-10-15.md).

Round-4 changes installed on 2026-10-16 with the user's approval of items A to N and P to S of the [approval package](../Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/wave4/APPROVAL_PACKAGE.md): the new title and central sentence, Fig 1B on the 3,769 COLOC-evaluated genes, Fig 2B to 2D, Fig 3 to 7 panel values, Catalog r9d (r9c plus the H79 and H81 fields) and rebuilt Fig 5 and Fig 6 assemblies. Item O was applied on 2026-10-16: the benchmark registry now records that H90 opened all 48 GSE289173 expression matrices and the disease and sex labels, so its cell-state task is no longer blind; the source cell-state labels and eQTL/ieQTL summaries stayed unopened, and no Figure 7 number used GSE289173. Overwritten files are archived under `Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/round4-20261015T112148Z/wave4/install-20261016T104914Z/archive/`.

### Historical missing-work inventory, superseded 2026-10-14

The following dated inventory retains earlier inspection states. It is not
the current work queue; current genetics, bulk, Catalog and panel values are
given above and in RESULTS.

1. **Corrected COLOC propagation.** The rerun is complete and adopted
   (2026-08-17), so this is no longer a compute dependency. Its signal eligibility
   and trait drivers are now under correction (see the
   [2026-09-23 scientific review block](#what-is-official-updated-or-open)). What is missing now is that
   some downstream consumers may still carry July values. **This list is
   measured, not assumed — several of its original entries proved already
   current.** Verified on 2026-08-20/21:
   - **Atlas: was already current.** 0 of 27,187 genes differed from the
     adopted release before any rebuild. The apparent "706 vs master 745" gap
     is a gene-universe difference (38 adopted SuSiE-positive symbols are
     absent from the atlas universe), not staleness. The chain
     `27a→45a→75→217→46d→46e→27b` was rerun anyway on 2026-08-20 (jobs
     `20951364`-`20951370`), refreshing 12 `multi_evidence/` products.
   - **Number ledger: already current.** Every 473/736 in
     [RESULTS.md](RESULTS.md) sits in an explicitly superseded row; 462/745 are
     the stated values for the 2026-08-17 release (corrected candidate: 436/708). There is no `docs/manuscript/NUMBERS.md` at the path the
     adoption note names.
   - **`relayreg`/`orientsub`: inputs verified equivalent.** The master
     aggregate and 8/8 sampled shards are byte-identical between
     `susie_coloc_rerun/` and canonical, and the relay reads shards only, never
     `gene_level_coloc.csv` (which does differ: `b1c890c4…` vs `fcb593df…`).
     Their outputs are stale only because they *ran* on 2026-08-15, before the
     08-16 eQTL merge. ⚠ Regeneration is **not scheduled**: steps 00-07 head a
     64-step pipeline whose steps 14 and 22-64 are experimental-study Stage A execution,
     guide design, and collaborator/clinical outreach — the separate-paper
     boundary. Needs an explicit scope decision before rerunning.
   - Still genuinely open: `fig4a_overview_candidates.py`, evidence classes,
     Gene Catalog genetics fields, and the worked-example posteriors. The latter
     two additionally wait on the bulk adoption, not on COLOC.
   - The "~30 consumers" figure understates *references*: 234 live code files
     mention `gene_level_coloc` or `results/susie_coloc` (GWAS 80, scripts 75,
     RNA-seq 48, Analysis 20, other 11), excluding archived snapshot copies. A
     reference is not an active consumer, and a token/mtime sweep cannot tell
     them apart — only a per-consumer value-level join against the adopted
     release can. Do not trust a blanket staleness list.
   Rebuilt already: the portal ancestry export, the Figure 1 GWAS panels,
   canonical Figure 2 A-I, and independently validated early/late supplementary
   Figure 2 candidates.
2. **One synchronized bulk adoption.** The updated five-cohort result is the
   intended paper substrate, but stage, composition, figures, and ledgers must
   be regenerated from the same read-only inputs. The 774-unit fit
   (`bulk-corrected-controls-20261008T162308Z`, 1,854 canonical DEGs) was
   adopted by the user on 2026-10-11 and replaces the 844-record
   corrected-controls fit (`bulk-corrected-controls-20260923T231939Z`, 1,725
   canonical DEGs). Figure 3, the stage and F0 arms, composition, the continuum
   layers and the evidence classes are rebuilt on it as candidates, as are
   Figure 3F, the Figure 4E program models and the continuum
   composition-adjusted coefficients (2026-10-12, not installed).
3. **HMSMA metadata.** Clinical, histologic, donor, and MALDI claims remain
   unavailable until an authoritative array-to-donor and registration table arrives;
   this does not hold up label-blind coverage or the rest of Figure 5.
4. **Final release assembly.** The July release is historical and incomplete for
   the present Resource requirements.
5. **Noncoding genetics synchronization.** The lncRNA strand rule and a complete
   DNA-only candidate are resolved. On COLOC r2 (adopted 2026-10-14), the
   target-gene companion (424 protein-coding, 12 lncRNA), Figure 2C and the
   exact-example source rows are installed; Gene Catalog r8 was built on r2
   before its adoption.
6. **Paired RNA and H3K27ac donors.** The Figure 7 cis head is still limited by donor
   count: out-of-fold skill rises from 0.032 at 20 to 0.143 at 78 training participants
   per fold without a plateau, and its gain over the no-cis recipe rises from 0.020 to
   0.095. A second paired liver cohort would add skill; more severity training samples
   would add little (0.6357 to 0.6631 over a fourfold pool). See [RESULTS.md](RESULTS.md).

## Canonical DEG threshold change (2026-08-12)

The canonical bulk DEG definition changed from the interval-null threshold
(`treat_fdr < 0.05` at `lfc = 0.25`) to the conventional `padj < 0.05` with
`|log2FC| > 0.50`, for comparability with the published RNA-seq literature.
On the current named release this moves the count from 1,918 to **1,853**
(1,438 up / 415 down). The genetic non-overlap that accompanies it depends on
which COLOC release is in the denominator, so quote it with its release:

| | July COLOC (447 / 14,931) | adopted 2026-08-17 (437 / 14,920) | corrected candidate, 2026-09-25 (409 / 14,243) |
|---|---|---|---|
| non-overlap, canonical threshold | 94.85% | **95.4%** | 94.1% |
| overlap odds ratio | 0.690, p = 0.094 | **0.616, nominal p = 0.037** | 0.75, nominal p = 0.22 |
| same test, fragment substrate | — | 0.661, p = 0.101 | 0.72, nominal p = 0.15 (Figure 3G interface universe, 406 / 14,183) |

The 2026-09-25 column is on the 844-record bulk fit. On the adopted 774-unit
refit with COLOC release `susie_coloc_r2_20260924T053244Z` (adopted 2026-10-14),
restricted to the 3,769 genes evaluated by multi-signal COLOC (H57): 426
genetic, 220 disease-state, 26 overlap against 21.21 expected within 10
AveExpr quantile bins, MH OR 1.295 (95% CI 0.835 to 2.007), exact P 0.2755;
expression not fully controlled, no direction stated. The 14,720-row values
(33.4 expected, OR 0.76) are withdrawn.

The earlier figure of 0.752 (p = 0.214) in this section predated the 2026-08-13
background correction and is retired. **Do not read the nominal p = 0.037 as a
finding**: it is one test in a family, it is not significant on the intended
fragment substrate, and the expression-matched ladder is flat, so colocalized
genes must still not be described as buffered or constrained. See the DEG-threshold
sensitivity section of [RESULTS.md](RESULTS.md).

- The single definition lives in `load_figure_data.R::is_canonical_deg()`.
  `is_dream_deg()` is now an alias for it, so scripts calling the shared helper
  follow automatically. `is_treat_deg()` exposes the retired threshold for the
  sensitivity arm. Do not re-derive a threshold inline.
- The `treat_*` columns stay in `canonical_deg_results.csv`. They back the
  calibrated comparator and the recorded F_five and lncRNA candidate
  disagreement assertions; deleting them breaks those checks.
- **Current candidate state:** the fragment-substrate primary-threshold count is
  validated. Figure redesign v2 stopped during supplementary composition
  rendering; v3 render job `20161509` failed closed after identifying that the
  UMAP renderer had not applied the two prespecified low-cell GSE189600 library
  exclusions. Candidate v4 validated that census and then failed closed on a
  program-renderer syntax error. Candidate v5 job `20161525` rendered every
  panel that was not withheld but failed closed in Ghostscript log decoding. Candidate v6
  completed PDF, checksum, cardinality, composite-proof, and scope validation
  in job `20161528`; manual proof review identified three label collisions.
  Candidate v7 resolved those collisions; candidate v8 narrows the final 4D
  legend to prevent edge clipping. Canonical PDFs remain
  unchanged until review and adoption.
- Deliberately left on other thresholds, by design, not oversight: the
  master-regulator and drug-target panels (`lfsr < 0.05`, so small-effect
  targets such as THRB survive) and the experimental follow-up library (its own analytic interval-null
  at `lfc = 0.2` for library sizing).
- RORA was convergent under interval-null (logFC −0.381) and fails the 0.50 floor. It
  must no longer be described as a convergent hit.

## Historical defect inventory and scientific dispositions

This inventory records the earlier defects and their then-current effects.
The alternate-contig recount, 774-unit stage refits, r2 posteriors, Catalog
and current panels are now corrected. GSE193066 remains unqualified. Legacy
July-release checks and historical hash records below do not certify the
October working figures; current sources and focused checks are recorded in
the [October 14 outcomes](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#remaining-resource-corrections-2026-10-14).

| Defect | Location | Effect |
|---|---|---|
| Bulk counts drop reads that also map to ALT/patch contigs | The bulk STAR BAMs were aligned to the 706-contig GRCh38.p14 and featureCounts ran without `-M`; 4,115 genes that share a name with an ALT/patch gene keep about a quarter of their reads | The paper-candidate five-cohort fit (corrected GSE130970 controls) has **1,833 canonical DEGs, not 1,725**, when those reads are folded back (featureCounts `-M --fraction`, same model; 128 gained, 20 lost; 24,650 genes tested instead of 23,702). The adopted 2026-08-10 fit would have 1,419, not 1,347. APOB average log-CPM is −1.04 in the current table and 13.6 after folding; COL3A1 −1.93 → 7.49; PKLR, TDO2 and GLUD2 are tested only after folding. Genes without an ALT/patch copy move by at most 0.019 in logFC; 18 of their 19 call changes sit within 0.01 of the |logFC| = 0.5 cut. Every figure, stage arm and Catalog field built on the bulk DEG table inherits this. Candidate rebuild: `Analysis/MASLD_Model_Benchmark/executions/alt-contig-recount-f-five-rerun-20260925T203044Z/` (`bulk-corrected-controls-alt-20260925T223034Z/`, gene list `paper_candidate_primary_deg_changes.tsv`); detail in `docs/technical/NEW_MODELS_PROGRAM_2026-09-23.md`. Adoption is a user decision. |
| Release-consistency checker requires two output files the current official release does not contain | `scripts/manuscript/check_release_consistency.py:78` requires `strict_drug_verdict.json` and `credible_set_coding_summary.tsv`; both exist only in the July 10 release, not the July 15 one | Retired for release 2026-10-17-r2. Its required-file list was written for the July 10 release and fails on later releases. Release 2026-10-17-r2 is checked by `scripts/manuscript/check_release_manifest.py` (checksums plus main-panel source tables) and the scope validator. |
| Two rendered stage panels carry superseded numbers | `figures/main/fig3_bulk_transcriptomics/_legacy/2026-08-11_pre_fragment_stage/` | Moved out of the main panel directory on 2026-08-11; see the stage section of RESULTS. |
| Source comment quotes a retired diagnostic margin | `GWAS/finemapping/src/06_susie_coloc.R:43-45` still justifies the LD cache with a 2.5e12× `lambda_s` margin; the true value is 1.66e9× | Overstates the safety margin by ~1,500×. Corrected in `docs/archive/progress/2026-08-11_coloc_ld_cache.md`; fix the comment when that file is next touched. |
| One ledger row has no recorded producer | `docs/RESULTS.md`, the interval-null-comparator overlap row (33 / 92.5%) | The 2026-08-17 interface recompute emits the canonical threshold only, so this number is not reproduced by any run on the adopted release. Recompute or drop it before it enters prose. |
| The finalized ATAC replay submitter assumes a node-local scratch parent that does not exist | `Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/slurm/08_genetics_batch.sbatch` calls `mktemp -d /scratch/$USER/...`, but `/scratch/jameslee` is absent on the compute nodes and `mktemp` does not create parents | Every replay batch dies on its first command. The file is hash-pinned in `genetics/implementation_validation/source_manifest.tsv`, so it was not edited; batches were run on 2026-08-19 through `candidates/atac-v3-genetics-upstream-manifest-20260819T144901Z/01_replay_batch_wrapper.sbatch`, which re-checks the recorded sha256, creates the directory, and `exec bash`es the recorded script unmodified. Fix the recorded script and record its new hash, or keep the wrapper as the documented submission path. |
| ATAC replay reproduction tolerance was amended after it was finalized | `Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/22_validate_replay_batch.py:130` now asserts `rel_tol=1e-4` (was `1e-7`); `GENETICS_IMPLEMENTATION_READY` and `genetics/implementation_validation/` still carry the pre-amendment hash `31837b0c…` rather than `d5ad5db2…` | Non-EUR strata fall back to rebuilt 1000G Gram LD panels that are not positive semidefinite, so `susie_rss` zeroes negative eigenvalues and the adopted run is not bit-reproducible: across the complete set of all 816 validated pairs the median difference is 5.4e-12, the max relative 6.6e-5 and max absolute 4.5e-5, with 42/816 over 1e-7 and 0 over 1e-4 (an earlier 116-pair mid-run sample reported max 1.9e-5 and is superseded — it understated the tail 3.5×). The tolerance is bounded by the selection margin — the adopted pair closest to the 0.5 cut is 0.503716, **82×** the largest absolute difference — so no gene-level membership call can move. ⚠ Two residual risks: the observed max is only 1.5× below the 1e-4 assertion, and the assertion covers pair-level `PP.H4` only — the per-variant `SNP.PP.H4` vectors that Figure 5D consumes have no adopted reference to check against at any tolerance. Recorded in `REPLAY_TOLERANCE_AMENDMENT.md`, `replay_tolerance_amendment.tsv`, and `replay_tolerance_amendment_hashes.tsv` in the candidate root, following the `PEAK_CALL_AMENDMENT.md` convention. `20_validate_genetics_implementation.py:140` refuses to overwrite an existing hash record, so re-recording requires deleting the 2026-08-12 fixture record and explicit approval. |
| ATAC liftover unpacking was vectorized after the same finalization | `Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/08_liftover_variants.R:126-149` replaced a per-element `lapply` over the `liftOver` `GRangesList` with `elementNROWS`/`unlist`/`cumsum`; sha256 `90ac8e71…` → `cc74fea8…`, and the recorded hash is still the pre-amendment one | The loop cost 5.62 ms per unique variant (962 s for 171,301 real variants) and retained 1,845,061 single-row `data.table`s before `rbindlist`. Over 22,356,697 posterior rows that killed job `20918763` OUT_OF_MEMORY at 30.7 GB of 32 GB after 15:03:40 with no output written. The replacement is `identical()`-verified against the recorded implementation on synthetic edge cases and on 171,301 real variants, and is 50,635× faster; its one uncovered path is multi-mapping variants, absent from the real-data test and covered synthetically only. Recorded in `LIFTOVER_VECTORIZATION_AMENDMENT.md` and `liftover_vectorization_amendment.tsv`. Rerun is job `20948247` at 128 GB on `--qos=interactive`. |
| Genetic figure page-size assertion contradicted its own renderer | `11_validate_genetics.py:275-276` asserted 388.8 x 201.6 pts within `abs_tol=0.2`, but `13_render_genetic.R:110-111` renders via `cairo_pdf`, which truncates the page box to whole points and wrote 388 x 201 | The check could never pass on a cairo host; no fixture rendered a real PDF, so it surfaced only on first real-data execution (job `20948356`). Amended to `abs_tol=1.0` (1/72 inch), recorded in `FIGURE_PAGESIZE_AMENDMENT.md` and `figure_pagesize_amendment.tsv`, sha256 `78376e1c…` → `503e98be…`. The renderer was deliberately left on `cairo_pdf`: base `pdf()` emits exact fractional points but risks the Type 3 fonts the next assertion rejects, and the recorded non-genetic validator asserts no numeric page dimensions at all (`10_validate_non_genetic.py:375`). |
| Retracted probe run wrote real readiness files | A probe copy of `11_validate_genetics.py`, differing only in the page-size tolerance, was run to enumerate remaining failures; it passed and wrote `GENETIC_READY`, `FULL_READY`, and `genetics/validation/` | Those readiness files were produced by an unmanifested ad-hoc file and did not record the code that produced them. They were moved to `candidates/atac-v3-genetics-upstream-manifest-20260819T144901Z/superseded_probe_run_20260820T192500Z/` and superseded by the recorded rerun `20948527`. Comparing the two runs afterwards showed every data file byte-identical, so the retraction cost nothing scientifically, but a validator that writes readiness files must not be used as a dry run. |
| "844 participants" and "704 disease" in the paper-candidate five-cohort bulk fit are at most 813 people and at most 673 disease people, and its stage contrasts place the same person on both sides (4 people in NAS 3–4→5–8; 1 in F0→F1 and in F1 vs F0): genotypes place 32 same-person library pairs and 1 first-degree pair inside the fit (added 2026-09-26; revised 2026-09-29) | Sample list `RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20260923T231939Z/primary_source_controls/model_design.tsv` (844 `sample_id`; the 2026-08-10 fit `candidates/resource-f-five-coloc-v3-candidate-2026-08-10/workstreams/BULK-POOLED-REPRO/model_design.tsv` lists the same 844) and stage roster `figures/candidates/bulk-release-corrected-controls-20260925T175215Z/analysis/stage_extensions/stage_extension_sample_manifest.tsv`, each joined to the genotype identity tables `Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/identity/kinship_related_pairs.tsv` and `sample_to_individual.tsv` (job 0989; fit pairs counted 2026-09-26, Figure 3D and stage rosters 2026-09-29). Participant counts (844, 704, 687, 416, 303) are printed in this file at lines 53, 167, 170 and 404; `docs/RESULTS.md:55-57`, `:84-85`, `:116`, `:581`, `:637`; `docs/manuscript/draft/fig3.md:9-10`, `:23-25`, `:32`, `:39`; `docs/manuscript/METHODS.md:45`, `:49`; `docs/manuscript/05_figure_legends.md:89`, `:96-99`; `docs/HIGH_IMPACT_PATTERNS.md:44`. DEG counts that treat each library as a participant (1,725; 1,347) are printed in this file at lines 53, 167 and 286; `docs/RESULTS.md:59`, `:88`, `:103`, `:353`, `:710`, `:1099`; `docs/manuscript/draft/fig3.md:16`, `:36`; `docs/manuscript/METHODS.md:50`, `:54`; `docs/manuscript/05_figure_legends.md:85`; `docs/HIGH_IMPACT_PATTERNS.md:44`. Stage-contrast participant counts: `docs/RESULTS.md:117-118` and the table at `:646-658`. No main-figure PDF prints 844, 704, 687, 416, 303, 1,725 or 1,347 (`pdftotext` over the 649 PDFs under `figures/main/`, `_legacy` excluded) | Genotypes place 32 same-person pairs inside the fit (GSE213621 24, GSE135251 5, GSE130970 3; 30 people, one of them with three libraries) and 1 first-degree pair (one GSE213621 library and one GSE130970 library, kinship 0.28); accessions in `Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/seal-20260926T155512Z/resource_identity_defects.tsv`. The 831 fit libraries with an identity call come from 800 people; 13 GSE126848 libraries have too few genotype sites to check. So the "844 participants" are at most 813 people. All 31 repeat libraries are in the disease group under both label sets, so "704 disease" is at most 673 (2026-08-10 labels: 687, at most 656). Five passages state that each participant enters once or is the inferential unit, and the repeat libraries contradict them: `docs/RESULTS.md:108-109`, `docs/manuscript/05_figure_legends.md:88`, `docs/manuscript/draft/fig3.md:19-20`, `docs/manuscript/draft/fig3.md:61-63` (Figure 3's stated contribution includes "biological participants as the inferential unit in every fit") and `docs/manuscript/METHODS.md:126-127` ("Each analysis retained one record per participant", said of the stage contrasts and contradicted by them below). Numbers it could change: the 1,725 canonical DEGs (`RNA-seq/results/manuscript_release/candidates/bulk-corrected-controls-20260923T231939Z/arm_deg_counts.tsv`, arm `primary_source_controls`), 1,833 after the ALT/patch fold-back (`Analysis/MASLD_Model_Benchmark/executions/alt-contig-recount-f-five-rerun-20260925T203044Z/bulk-corrected-controls-alt-20260925T223034Z/arm_deg_counts.tsv`, same arm) and 1,347 in the adopted fit (first file, arm `adopted_2026-08-10`), all of which count the 31 repeat libraries as independent people. Figure 3D roster (`figures/candidates/bulk-release-corrected-controls-20260925T175215Z/source_tables/fig3d_endpoint_pca_coordinates.tsv`): the 416 participants are 404 people by genotype and the 303 advanced-disease participants are 291 (12 repeat libraries: GSE213621 7, GSE130970 3, GSE135251 2); the 113 strict controls are 113 people, and no person is in both groups. Stage contrasts (manifest above; every manifest library has an identity call), libraries → people: adjacent fibrosis F0→F1, F1→F2, F2→F3, F3→F4 209, 195, 165, 100 → 206, 192, 161, 97 (`docs/RESULTS.md:117-118`, `:648-651`); fibrosis versus F0 for F1, F2, F3, F4 209, 194, 179, 129 → 206, 193, 176, 129 (`:652-655`; the 75 F3 libraries are 72 people); adjacent NAS 0→1–2, 1–2→3–4, 3–4→5–8 111, 186, 296 → 111, 186, 289 (`:656-658`; the 175 NAS 5–8 libraries are 172 people). The NAS-versus-NAS 0 contrasts in the same manifest, not printed in RESULTS, are 111, 167, 221 → 111, 167, 218. People with libraries on both sides of one contrast: in NAS 3–4→5–8, four people: GSE130970 SRR9036308 (NAS 3–4) with SRR9036312 (NAS 5–8), manifest lines 1744 and 1748; GSE135251 SRR9883105 (NAS 3–4) with SRR9882960 and SRR9883106 (NAS 5–8), lines 1919, 1806 and 1920; GSE135251 SRR9883022 (NAS 3–4) with SRR9883014 (NAS 5–8), lines 1853 and 1846; GSE135251 SRR9883159 (NAS 3–4) with SRR9883149 (NAS 5–8), lines 1962 and 1954. In F0→F1 and F1 vs F0, one person: SRR9882960 labelled F0, SRR9883105 and SRR9883106 labelled F1 (manifest lines 118, 184-185 and 787, 853-854). The same-person call for these three libraries is at `resource_identity_defects.tsv:7-9`; those rows mark the libraries stage-concordant because the `stage_coding()` used there puts F0 and F1 in one class (S0). A refit with one library per person has not been run, so the size of any change is not known. A refit with one library per genetic individual restores the donor unit; if these numbers stay on the public sample labels, `docs/PAPER.md:579-582` (a main-paper analysis must not weaken biological unit control) and `:586` (stop or move to the supplement when the public metadata do not identify donors) apply to the 844-participant numbers and the stage contrasts. Refitting and adopting a one-library-per-person fit is the user's decision. |
| The 361 GSE213621 participants in the stage-adjusted program and system models are at most 342 people: genotypes place 25 same-person library pairs in GSE213621, and 11 of the pairs carry two different fibrosis stage labels (added 2026-09-26; revised 2026-09-29) | `identity/kinship_related_pairs.tsv` and `identity/sample_to_individual.tsv` in the execution above; stage labels are the GSE213621 `condition` values in `RNA-seq/Human/Patient_Cohorts/analysis/integration/metadata/unified_metadata.csv`, coded by `stage_coding()` in `scripts/analysis/in_biopsy_allelic/tags/t5a_allelic_table.py:57-71`. The 361 is printed at `docs/RESULTS.md:27-30`, `:729-730` (continuum-window figures, GSE213621 n=361), `:802` and `docs/manuscript/METHODS.md:341-342` ("109 and 361 complete participants"); the 299 at `docs/RESULTS.md:2035` and `:2207` (label-ceiling table). `docs/manuscript/METHODS.md:226-227` states that "Every continuum analysis used biological participants as the inferential unit" | The source paper's Methods give 368 human subjects, 69 without chronic liver disease and 299 with biopsy-proven NAFLD, sampled at a diagnostic liver biopsy (PMID 36535507, PMC9957752, doi 10.1016/j.jcmgh.2022.12.008); it counts one sample per patient and describes no repeat biopsies. GEO GSE213621 lists 368 samples with 368 distinct titles (series matrix `https://ftp.ncbi.nlm.nih.gov/geo/series/GSE213nnn/GSE213621/matrix/GSE213621_series_matrix.txt.gz`, retrieved 2026-09-29, field `fibrotic stage`: Control 69, F0F1 97, F2 107, F3F4 95). `unified_metadata.csv` has 367 GSE213621 libraries (Control 68, F0F1 97, F2 107, F3F4 95), so one GEO Control sample is not in the local metadata; that accounts for 368 against 367 and is separate from the same-person finding. Genotypes place the 367 in 342 people (317 with one library, 25 with two). All 25 pairs are in the staged group, so the 299 staged libraries come from 274 people; the 68 Control libraries are 68 people. Of the 25 pairs, 14 carry the same stage label and 11 do not: 6 F2 with F3F4, 4 F0F1 with F2, 1 F0F1 with F3F4. The source design does not explain 25 genotype-identical pairs; these data cannot separate a duplicated sample from a sample swap or a label error. 24 of the 25 pairs, including all 11 discordant ones, are in the 844-participant fit. Numbers it could change: the 361 GSE213621 participants in the stage-adjusted program and system models (93/117 programs over 470 participants). Any 361 of the 367 libraries keep at least 19 of the 25 repeat libraries, so they are at most 342 people, and 337 if they are the 361 fit libraries [Inference: the model roster was not joined]; the 299- and 367-record GSE213621 external severity checks (`docs/RESULTS.md:2035`); the GSE213621 row of the label-ceiling table (299 participants; simulated ceiling 0.779 and 0.873, achieved / ceiling 0.820 and 0.731; `docs/RESULTS.md:2207`); the GSE213621 continuum-window figures (n=361, `docs/RESULTS.md:729-730`); Model B's 299 GSE213621 records, which are 274 people (candidate only; `docs/technical/NEW_MODELS_PROGRAM_2026-09-23.md`); and any stage contrast built on GSE213621 labels. A refit with one library per genetic individual restores the donor unit; if these numbers stay on the public sample labels, the `docs/PAPER.md:586` rule (stop or move to the supplement when the public metadata do not identify donors) applies to them. Model A-final excludes same-person pairs with discordant labels. Applying an exclusion to Resource analyses is the user's decision. |
| The 105-participant GSE193066 coordination evaluation cohort (edge sign concordance 0.8974), the 58-participant repeat-biopsy numbers (within-participant SD 0.6433, ICC 0.7658 vs 0.1886) and the label-ceiling scenario built on SD 0.6433 rest on a cohort in which 132 of 164 libraries meet the possibly-mixed rule; donor identity there cannot be checked until the contamination-aware identity rerun (added 2026-09-26; revised 2026-09-29) | `identity/sample_to_individual.tsv` in the execution above; the rule (expression sex "both", or other-allele fraction at homozygous sites > 0.02) as run is `rules.possibly_mixed` in `Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/seal-20260926T155512Z/a2_summary.json` (132 at `per_cohort.GSE193066.libraries_excluded_by_reason.possibly_mixed`, job 1000); the current rerun, job 1002, gives the same rule string and the same 132 (`Analysis/MASLD_Model_Benchmark/executions/in-biopsy-allelic-a1-full-20260924T025858Z/seal-20260929T042934Z/a2_summary.json`); marker thresholds at `scripts/analysis/in_biopsy_allelic/a1_identity/03_identity.py:61-63`. The affected numbers are printed at `docs/RESULTS.md:1488-1494`, `:926-928`, `:1522`; this file lines 15-16; `docs/manuscript/draft/title.md:30-35` (89.7%, 105 participants); `docs/manuscript/draft/fig7.md:59-60` and `docs/manuscript/05_figure_legends.md:324-326` (within-participant SD 0.4655 to 0.6433); `docs/manuscript/METHODS.md:680` (58 paired participants); `docs/RESULTS.md:2036` and this file line 184 (repeat-biopsy ICC 0.7658 vs 0.1886, 58 participants); `docs/RESULTS.md:2194-2210` (label-ceiling simulations at SD 0.6433: simulated ceilings 0.799, 0.779, 0.710; achieved / ceiling 0.830, 0.820, 0.949; within-v1.2 spread 0.130); `docs/RESULTS.md:950` and `:1533` (range 0.4655 to 0.6433); `docs/RESULTS.md:981-983` (58 paired participants, stage variance 0.74 times that of the full 106); `docs/manuscript/METHODS.md:224-225` (paired GSE193066 models) | 130 libraries have an other-allele fraction above 0.02 (cohort median 0.029, against 0.0044 across all 1,466 libraries in `identity/identity_summary.json`); 24 express both XIST (log2 CPM > 5) and Y-chromosome genes (log2 CPM > 3); 22 meet both rules. Only 141 of 164 have enough sites for an identity call. An elevated other-allele fraction, whether from mixed RNA or a cohort-wide technical cause not yet separated (per-library error estimate pending), makes identity calls in this cohort unreliable. Provisional identity calls, not genotype statements until the contamination-aware rerun: against the source crosswalk `RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-STAGE/PREFLIGHT-v1/audits/gse193066_participant_crosswalk.tsv`, 29 of the 58 recorded repeat-biopsy participants are called one person, 11 pairs with calls on both libraries are not, and 18 cannot be checked; 7 of the 36 same-person calls join libraries of different recorded participants. Numbers it could change: the sealed coordination replication (105 participants, edge sign concordance 0.8974; `docs/RESULTS.md:1488-1494`, this file lines 15-16), the repeat-biopsy stage SD of 0.6433 from 58 participants (`docs/RESULTS.md:926-928`), the 58-pair within-person test (Spearman 0.236; `docs/RESULTS.md:1522`), the repeat-biopsy ICC 0.7658 vs 0.1886 (`docs/RESULTS.md:2036`, this file line 184) and the SD-0.6433 label-ceiling scenario (`docs/RESULTS.md:2194-2210`). Only the 0.6433 end of the 0.4655 to 0.6433 range comes from GSE193066; 0.4655 comes from 30 PXD051911 pairs (`docs/RESULTS.md:2196-2197`). If the rerun confirms the crosswalk conflicts, the `docs/PAPER.md:586` rule (stop or move to the supplement when the public metadata do not identify donors) applies to the 105-participant evaluation result and the repeat-biopsy numbers. GSE193066 is not in the 844-participant fit. Model A-final excludes GSE193066 from allelic statistics. Acting on this for Resource numbers is the user's decision. |

## Work that can proceed now

Numerical and panel corrections listed in the historical inventory are
complete. Publication assembly, references and the overall ROADMAP checklist
remain. The bullets below retain the earlier work allocation and are not
instructions to repeat completed fits.

- Rewrite introduction, discussion, and non-numerical Results around the
  Resource contribution and claim boundaries.
- Finalize Methods sections whose inputs and models are already settled.
- Prepare figure layouts, legends, and data dictionaries with placeholders for
  values that are explicitly not yet available.
- Rebuild bulk-dependent stage and composition analyses in isolated candidate
  paths.
- Review the 117-program Figure 5E candidate and the corrected Figure 5F source
  table (`repair/fig5f_spatial_program_calibration_corrected_r2.tsv`) in
  `figures/candidates/review-r3-fixes-20261015T112148Z/F5ef/`; positive Moran
  direction alone is not the main spatial result.
- Draft supplementary legends for the validated composition, unit/null,
  cell-context, and HMSMA proxy-mask diagnostics without upgrading their claim
  boundaries.
- Audit portal/Gene Catalog schemas without publishing genetics-dependent entries;
  the spatial candidate products are ready for synchronized integration.
- Retain the finalized ATAC v3 replay and canonical Figure 5D against the adopted
  2026-08-17 COLOC manifest; do not convert peak overlap into target assignment.
- Use the validated GENCODE identity, pooled biotype subdivisions,
  fragment-native stage counts, and prespecified-program lncRNA annotations in
  candidate figures and Methods; do not adopt them before release
  synchronization.
- Reconcile references and competitor positioning.
- Integrate the validated histology-anchored continuum qualitatively as a
  complementary molecular-resolution layer; keep exact candidate values tied
  to the checksummed result ledger until synchronized adoption.

Exact genetic counts and overlap percentages are now authorized **from COLOC r2,
adopted 2026-10-14**, as ledgered in [RESULTS.md](RESULTS.md); the July and
2026-08-17 values are superseded. Gene examples and Gene Catalog classes use
their r2 rebuilds.

## Known interpretation boundaries

- The adopted bulk analysis contains **774** pooled disease-vs-control units
  (genetic-identity clusters, one library each; 844 records before identity
  resolution);
  the complete corrected fragment source has **1,257** QC-passing samples.
  These are different denominators for different purposes.
- The prior six-cohort fragment-native stage candidate passed its own checks but
  is provenance only. Figure 3 requires the active five-cohort synchronized
  rebuild and one shared release ID.
- The single-cell analysis must use biological donors, not sequencing runs or sort
  fractions. Cross-cohort inferred fibrosis stage is not authorized.
- NMF k4/k6 factors are continuous interpretive axes, not reproducible discrete
  patient subtypes.
- The three disease aspects of metabolism, inflammation, and fibrosis are **not
  resolved as separate outcome axes**. NAS enters as a graded sum and as an
  adjustment set; fibrosis stage is the severity anchor. Separability of the
  three NASH-CRN components is **`indeterminate` in GSE267145**. In 99
  participants the histology grades alone fail two of the three prespecified
  criteria (Spearman 0.68/0.73/0.82; label participation ratio 1.42 against a 2.0
  floor), so the rule could not return a positive; the RNA criterion also fails
  (association-vector correlations 0.83/0.87/0.87 against 0.80). In simulation
  on these grades with independent genes, one shared program and three programs
  that each follow one recorded grade give the same values; programs tied to the
  part of each grade not shared with the other two give values the data do not
  show (RNA criterion below 0.66; H3K27ac participation ratio 1.48 to 2.08
  against the observed 1.15). Aspect **attribution stays `indeterminate`**;
  nothing can be assigned to a single aspect on this evidence. The result
  concerns one cohort and does not address the SAF partition of
  steatosis, activity and fibrosis (SAF grades in GSE202379 and PXD051911 remain
  unused; GSE135251 carries only the NAS sum).
- **Variant combinations are untested** beyond a single characterised local
  haplotype at GNMT, which is reported as unresolved. The colocalization design
  is single-variant and does not test combinations or interaction.
- **Chromatin supplies context, not effect.** The large majority of jointly
  testable peaks resolve to `indeterminate` and a minority to `source_dependent`.
  No epigenetic effect or causal claim is available.
- Yakubovsky is a normal liver zonation/lipid-context reference, not a MASLD
  disease-validation cohort.
- Vu arrays remain unresolved physical units; directional replication is
  `source_dependent`, not donor-level inference.
- HMSMA contributes label-blind coverage and descriptive per-array organization
  only until the authoritative metadata and image-registration joins arrive.
- GSE287826 cannot support donor-level inference without a public donor key.
- No currently usable spatial cohort supplies every desired MASLD control,
  donor, histology, and metabolite join.
- Proteomics is not applicable to lncRNA molecules. Chromatin observes a locus
  or regulatory context, not the RNA molecule itself.
- ATAC v3 separates static observability from MASH-associated accessibility.
  Four stellate-lineage peaks reach q < 0.05 in both cohorts with the same
  sign. In 45 of 2,002 GSE244832 donor relabelings at least four peaks do so
  (p = 0.0225, a single test added in review; null mean 0.17).
  Stellate directional replication has
  donor-relabeling BH q = 0.074, so the lineage does not meet the dynamic-main
  eligibility rule. The four peaks are candidates, not chromatin findings, and
  not causal-gene or causal-cell-type proof; three of the four promoters are
  unnamed ENSG loci. The GSE281367 T/NK label remains a
  combined, source-dependent observability context.

## Source of truth

- Scientific scope and claims: [PAPER.md](PAPER.md)
- Current numeric ledger: [RESULTS.md](RESULTS.md)
- Execution sequence: [ROADMAP.md](ROADMAP.md)
- Detailed release history: [technical/RELEASE_PROVENANCE.md](technical/RELEASE_PROVENANCE.md)
- Full superseded handoffs and dated logs: [archive/INDEX.md](archive/INDEX.md)
