# Resource review and hypotheses, round 3

**Date:** 2026-10-15
**Purpose:** report what the third review round found in the human MASLD Resource manuscript, and set out the analyses it proposes. This document is a review. It adopts nothing. It changes no result, figure, Gene Catalog entry or authority document. None of the proposed analyses has been run.

## Sources and how to read this document

RUN below means `Analysis/Hypotheses/RESOURCE-REVIEW-20261015T020503Z/`.

Inputs:

1. The five authorities: [PAPER](../PAPER.md), [STATUS](../STATUS.md), [RESULTS](../RESULTS.md), [ROADMAP](../ROADMAP.md), [README](../README.md). Also the draft manuscript, the figure legends and the installed figure PDFs.
2. Prior reviews: [combined review 2026-10-07](COMBINED_SCIENTIFIC_REVIEW_2026-10-07.md), [corrections ledger](SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md), `RUN/context/known_issues.tsv` (K01 to K21), and the read-only archive `RUN/context/archive/Hypothesis_Review_2026-10/`.
3. The outputs of this round's finders, verifiers, gates, planners, plan audits, judges and portfolio critic, all under `RUN/scratch/` and `RUN/jobs/`.

Labels:

- **Confirmed:** a separate verification agent re-derived the finding from the producing files. Its scripts and outputs are in `RUN/scratch/verify-<ID>-reproduce/` and, for 12 findings, also in `RUN/scratch/verify-<ID>-coverage/`.
- **Unverified:** one finder reported it with its own commands, and no second agent re-derived it. Do not act on it before verification.
- **New / Known:** judged against K01 to K21, the 2026-10-07 combined review, the 2026-09-23 corrections ledger and the archive.
- Nominal p values are labelled as nominal and are never called significant.

Short names: LEG = `docs/manuscript/05_figure_legends.md`; ABS = `docs/manuscript/draft/abstract.md`; METHODS = `docs/manuscript/METHODS.md`.

## 1. Summary

This round checked the manuscript against its producing files. Sixteen review dimensions (R01 to R16) and three focused gap finders (G1 to G3) together raised 175 findings. 43 were merged as duplicates, which leaves 132. 40 of the 132 were verified by blind agents in two passes. The main pass covered the 20 ranked first, and all 20 were confirmed. A second pass covered the 20 remaining high-severity findings: 13 were confirmed, 7 partly confirmed, and none refuted (Section 2B). All five planted decoys were classified correctly (one true, four false). The other 92 findings (62 medium, 30 low) remain unverified. I re-derived three confirmed findings myself and got the same values (Section 2C).

The confirmed findings affect both title clauses, the abstract, and main panels 1B, 2A, 2B, 2D, 3C, 3E, 3F, 4C, 4H, 5E and 6B. They fall into five patterns:

- **Wrong denominators.** The Fig 1B universe counts 10,928 genes that multi-signal COLOC never evaluated as non-colocalized.
- **One dataset or ancestry produces a pooled effect.** In Fig 4C the stage slope is 0.107 in GSE202379 and 0.742 in GSE244832. In Fig 2D, 6,126 of 6,208 enzyme instances are non-EUR, and within one ancestry the ratio is 1.7 to 4.3, not 6.03.
- **Control labels enter as a stage.** In Fig 4H, 68 GSE213621 healthy controls are coded as fibrosis stage 0.
- **Count artifacts.** 437 of the 6,147 GSE126848 calls in Fig 3C are genes with zero reads in all 30 libraries. 1,865 of the 4,696 "shared" genes are called in opposite directions.
- **Claims with no analysis behind them.** Title clause 2 ("the coordination of multicellular remodeling") has no Results section, legend or panel. The Catalog has no regulatory DNA object, although the abstract says it links them.

The second pass adds three more:

- **Results that depend on one cohort or on one analysis choice.** Without GSE130970 the Fig 3E F1 bar drops from 1,338 calls to 14. 29 of the 54 Fig 3F genes depend on the same cohort. Fig 4E's 16 fibrosis calls are 0 on the 774 units alone, and 0 under truncated Knapp-Hartung.
- **A Catalog label that measures less than it says.** "Genetics" is liver-eQTL colocalization only, so MTARC1, PNPLA3, TM6SF2 and HSD17B13 read as unresolved.
- **Wording beyond the numbers.** This includes "resolves the architecture", a "direct MASLD" coding enrichment that mostly reflects the three PDFF GWAS (NAFLD and NASH alone: ratio 3.21 [0.53, 8.25]), and four uncited prior papers.

The round proposes ten analyses. By the plans' own estimates they need about 252 CPU hours of expected use, and H51 alone could reach about 400 if single-thread fits are slower:

- Three correct printed claims: H57 (Fig 1B), H58 (Fig 2D and 2A) and H59 (Fig 4C).
- Four test claims the title or novelty statement depends on: H52 and H54 (coordination), H60 (signed link between genetic and disease-state maps) and H51 (stage replication).
- Three add reader-facing Catalog fields: H53, H55 and H56.

All ten plans passed a plan audit after patches. Running them, and adopting any result, needs the author's decision.

## 2. Confirmed findings

The ranking puts title and abstract claims first, then main-panel numbers, then Results text. The three Figure 7 findings come last because Figure 7 is a validated candidate that has not been synchronized into the paper. Every value below was re-derived by a verifier unless it is labelled otherwise.

| # | ID | Where | Printed | What the evidence supports | Action (plan) |
|---|---|---|---|---|---|
| 1 | R13-01 | Title clause 2 | "resolves ... the coordination of multicellular remodeling" | No Results section, legend or panel reports a covariance result. The only analysis is historical, on 469 identity-unqualified records. | Reword (H52, H54) |
| 2 | R02-01 | Fig 2D legend and text, abstract | Ratio 6.026671 [2.57459, 13.51016] read as a trait-class effect | Trait class is aliased with ancestry. EUR only 4.32. Non-EUR only 1.90, with a difference interval that includes 0. | Reword (H58) |
| 3 | R01-01 | Fig 1B, abstract, Discussion | 14,704 "jointly testable"; expected 7.9% and 2.9%; DEG-only 1,129 (97.7%) | 10,928 genes were never evaluated. In the 3,776 evaluated genes: OR 1.057, expected 5.83% and 11.28%, DEG-only 194 (88.2%). | Recompute (H57) |
| 4 | R07-03 | Fig 1B, Discussion | 33.46 expected; OR 0.757 | Expression-matched expectation is 22.1 to 22.8. Spline-adjusted OR 1.201 (0.789 to 1.828). | Recompute (H57) |
| 5 | R04-02 | Fig 4C beta and q | Pooled beta 0.402 with homoskedastic SE | Within-dataset slopes are 0.107 (n=37) and 0.742 (n=18), Cochran Q p 0.016. IVW slope 0.144. | Recompute (H59) |
| 6 | R04-01 | Abstract, Fig 4C star | "Ductular injury alone met the primary 117-program rule" | Without the 2 uncorrected GSE189600 donors: q 0.0985 and 0 calls. 12 of 62 single-donor deletions lose the call. | Reword (H59) |
| 7 | R06-02 | fig5.md:41, 45-46; PAPER:336-337 | "Two equally robust disease programs" | Stromal ECM q 0.0855 is not selected. Neither program passes HC3. | Reword |
| 8 | R07-01 | Abstract, intro, Fig 6A | Catalog links regulatory DNA elements and tissue context | 0 of 28,233 records are regulatory DNA objects. No gene has a chromatin, spatial or protein call. | Reword or add layers (H53) |
| 9 | R05-02 | Fig 4H, S4J, METHODS:372 | rho 0.73; 37/43 same direction | Without the 68 GSE213621 controls coded stage 0: rho 0.643, 31/43. | Recompute |
| 10 | R03-02 | Fig 3C GSE126848 ribbon | 6,147 calls; 3,131 unique | 437 calls are all-zero genes. A per-cohort filter gives 3,446 calls. Depth is aliased with group. | Recompute |
| 11 | R03-01 | Fig 3C tiers; PAPER:392 | 9,065 / 3,462 / 922 / 261 / 51 read as recurrence | 1,865 of 4,696 multi-cohort genes are called in opposite directions. Same-direction tiers: 2,223 / 675 / 225 / 44. | Recompute |
| 12 | R03-03 | LEG:113; Fig 3A, 3D, 3E, 3F | "All panels use the corrected 774-unit substrate" | 3E and 3F use 387 units from 3 cohorts. GSE213621 contributes 0 units. | Doc fix (H51) |
| 13 | R02-02 | Noncoding-context sentence and companion panel | Unresolved 36.49% vs 77.85% | This is a failed lookup for non-EUR studies. After position repair it is 0% vs 0%. | Recompute (H58) |
| 14 | R02-04 | Fig 2B point positions | 695 points; 152 in the high-PIP, low-PP.H4 corner | With SuSiE x at y's own GWAS and lead: 327 points, 18 in the corner. | Recompute |
| 15 | R02-03 | Fig 2A legend; 2D denominator sentence | 360 fine-mapped loci feed 996 COLOC genes | 100 of 436 multi-signal genes have a positive lead in a fine-mapped window of the same study. | Reword (H58) |
| 16 | R01-03 | Fig 5E rows | 13 programs | Selected on the superseded August fit. 0 of 13 pass the corrected q<0.05. | Recompute |
| 17 | R07-02 | Fig 6B; PAPER:316; RESULTS:2699 | Zero tested_negative presented as a result | No negative rule exists, so the zero is fixed by construction. | Reword |
| 18 | R08-01 | Fig 7B legend (candidate) | Sequence models help on MPRA and not on severity | No sequence model was scored on severity. 1 of 7 MPRA rows survives Bonferroni. | Reword |
| 19 | R08-02 | Fig 7 legend (candidate) | 585 participants; 0.68 to 0.79 of tie ceiling | 96 of 177 GSE276114 records are not MASLD. MASLD ratio 0.866. | Recompute |
| 20 | R08-03 | Fig 7 legend (candidate) | 299 GSE213621 participants | 274 identity clusters. One record per identity gives ratio 0.709. | Recompute |

### 2.1 R13-01. Title clause 2 has no supporting analysis in the paper (high)

- **Where:** title (`docs/manuscript/draft/title.md:10-11`), README:102-104, STATUS:236-246.
- **Printed:** the Resource "resolves ... the coordination of multicellular remodeling". title.md:26 defines coordination as tissue covariance. The producing VERDICT.md:52-53 says "the title's coordination clause is earned by" the GSE193066 holdout.
- **What the evidence supports:** 637 main and supplementary figure PDFs were scanned. No draft Results section, legend or figure reports a program-pair covariance, edge, modularity or community result. The only covariance numbers in the draft are in title.md:27-35, labelled historical. The one covariance analysis has these limits:
  - It used 469 source records from GSE130970, GSE135251, GSE162694 and GSE174478, and covered 113 programs, not the 117-program registry.
  - METHODS:55-56 says these counts "must not be read as verified independent-donor n".
  - It was never refit on the 774 identity clusters, and it excludes GSE213621.
  - Its evaluation cohort, GSE193066 (n=105), is identity-unqualified and excluded from primary replication.

  The other evidence called "coordinated" does not test covariance between cell types. Single-cell: 1/117 programs pass the primary rule and 0/117 pass HC3. Systems: 0/43 pass the connected-region null and 1/43 the uniform null. H7 was untestable (14 F0/F1 donors against a rule of 15).
- **Re-derived:**
  - Discovery records: 76 + 214 + 86 + 93 = 469.
  - 6,328 possible pairs over 113 programs, 2,426 discovery-supported edges, 4 communities.
  - Cross-content C3 covariance after histology and PC1: discovery 0.2877 [0.2360, 0.3462]; GSE193066 0.1890 [0.1166, 0.2886]. These match RESULTS:298-301 and are not in the draft.
- **New or known:** known in part. K01 records the GSE193066 identity limit (see also OUTCOMES:1297). New here: no Results section, legend or panel presents the analysis.
- **Action:** reword now. The claim the evidence supports is "and maps cross-sectional multicellular remodeling". To keep a covariance claim, it would have to read: bulk program scores from different lineages co-vary across samples after histology adjustment, in retrospective, identity-unqualified rosters. H52 and H54 test whether an identity-qualified replication exists.

### 2.2 R02-01. The Fig 2D ratio of 6.03 is mainly a contrast between ancestries (high)

- **Where:** Fig 2D legend and the DNA-architecture paragraph in Results. The same sentence is in the abstract, intro, PAPER and STATUS (H58 lists the lines).
- **Printed:** "Coding plus splice mass is 23.34% versus 3.87%, ratio 6.026671 [2.57459-13.51016]", followed by "every difference interval has a positive lower bound".
- **What the evidence supports:** the two classes come from different GWAS ancestries and LD panels.
  - 6,126 of the 6,208 enzyme instances are non-EUR (BBJ EAS alone 3,689).
  - Direct-trait instances are EUR 106, AMR 95, AFR 5 and EAS 1 (n=207).
  - The LD panel is fully aliased with ancestry: every EUR instance uses ukbb_eur and every non-EUR instance uses 1kg_*.

  Within one ancestry the ratio is smaller, and outside EUR the difference interval includes 0. The pooled 6.03 is larger than either stratum. The 57-comparison robustness menu has no ancestry-stratified arm. Its LD-block arm (6.084745) changes only the grouping. The ascertainment models adjust only for standardized beta and N.
- **Re-derived:** the pooled value reproduces exactly: 0.233428 vs 0.0387324, ratio 6.02667 [2.57459, 13.5102]; clusters 19 direct, 228 enzyme, 233 physical, 14 shared. Stratified with the producer's estimator:

  | Stratum | Direct vs enzyme | Ratio [95%] | Difference [95%] |
  |---|---|---|---|
  | EUR only | 29.88% vs 6.92% | 4.319 [1.447, 24.219] | 0.2296 [0.0548, 0.4104] |
  | Non-EUR only | 6.99% vs 3.68% | 1.900 [0, 5.51] | 0.0331 [-0.0388, 0.1318] |
  | MVP only | | 1.672 | 0.0281 [-0.0405, 0.1201] |
  | AMR only | | 1.482 | 0.0219 [-0.0486, 0.1083] |
  | Direct EUR vs enzyme non-EUR | | 8.657 [3.96, 18.92] | |
  | Direct non-EUR vs enzyme EUR | | 1.0099 | 0.0007 [-0.107, 0.114] |

  Reweighting the enzyme fractions to the direct arm's ancestry mix (0.615 EUR) gives 3.716 [Estimate, no interval].
- **Added by the H58 plan audit:** the EUR stratum is not a clean contrast either. The EUR direct arm has credible sets in 47 of 47 registry lead windows. UKBB ALT, AST and GGT have them in 8 of 767 by locus key. UKBB_ALT and UKBB_AST list significant leads at PNPLA3 and TM6SF2, but neither has a UKBB credible set within 500 kb. After the audit, only AMR passes the coverage gate, and its estimate is the 0.0219 [-0.0486, 0.1083] above. The failed candidate L10-02 (section 5) shows that a genome-wide EUR fine-mapping would confound trait class with GWAS power. With the GWAS in hand, a clean trait-class coding contrast may not be estimable.
- **New or known:** known in part. RESULTS:1909-1910 says the August ascertainment control did not control ancestry or LD panel, but gives no stratified estimate. Not in K01 to K21, the 10-07 reviews or the 09-23 ledger.
- **Action:** stop printing 6.03 and "every difference interval has a positive lower bound" as a trait-class effect. Report the stratified values. This bears on title clause 1 ("the architecture of inherited risk"). H58 supplies the formal intervals. Related unverified items: R14-02 (the direct-class enrichment comes from the three PDFF GWAS) and R13-02 (title clause 1).

### 2.3 R01-01. The Fig 1B universe counts 10,928 genes never evaluated by multi-signal COLOC as genetic negatives (high)

- **Where:** Fig 1B; the abstract sentence "Among 14,704 jointly testable..."; the first Discussion paragraph; LEG:21.
- **Printed:** LEG:21 gives 14,704 jointly testable genes, with 426 colocalized, 1,155 disease-state and 26 shared, against 33.46 expected. The panel prints expected lines of 7.9% (colocalized row) and 2.9% (DEG row), observed 6.1% and 2.3%, and "Disease-state only 1,129 (97.7%)". The x axis reads "Fraction of native evidence interface".
- **What the evidence supports:**
  - 10,928 of the 14,704 genes (74.3%) have no SuSiE result (all NA PP.H4), and 935 of the 1,155 DEGs are among them.
  - The paper's own rule (fig2.md:27-29) treats missing SuSiE output as an observability state, not a negative. These genes cannot be eligible negatives in the expectation.
  - Of the 10,928, 10,884 had a SuSiE-COLOC attempt that fell back to ABF, and 44 had no eQTL SuSiE fit. 3,720 have an ABF eQTL signal, and 309 have ABF PP.H4 > 0.5.
  - The DEG rate is 8.56% in the unevaluated genes and 5.83% in the evaluated genes.
  - Fig 3G already restricts to the 3,776 genes with a finite PP.H4, so Figs 1B and 3G show the same 26-gene overlap against different denominators.
  - The x-axis label is wrong: each bar is a fraction of its own gene list, not of an interface.
- **Re-derived:**
  - Full universe reproduces: OR 0.7570, nominal Fisher p 0.2001.
  - Restricted to the 3,776 SuSiE-evaluated genes: genetic 426, DEG 220, both 26, expected 24.82, OR 1.057, nominal p 0.743. Expected fractions are 5.83% (colocalized row) and 11.28% (DEG row). DEG-only is 194 (88.2%), not 1,129 (97.7%).
  - Fig 3G prints Neither 3,156, Transcriptomics 194, Genetics 400 and Both 26, which equals this restriction.
  - In S3G, 317 of the 426 matched controls are no_susie_result genes.
- **New or known:** new in this form. A related asymmetric-denominator candidate (archive 04_CANDIDATE_POOL.md:3868, resource_design-04) was never executed. K10 covers only the 14,720 / 14,704 / 14,682 universe mismatch.
- **Action:** recompute (H57). The wording "near chance" survives. The name "jointly testable", the count 1,129 and the depletion reading do not.

### 2.4 R07-03. The Fig 1B expectation treats genes of every expression level as exchangeable (medium)

- **Where:** the Fig 1B expected markers; the Discussion sentence with OR 0.757 and 93.90% non-overlap against 92.14% expected (fig5_discussion.md:93-97).
- **What the evidence supports:** DEGs are mostly low-expression genes (mean AveExpr 0.463 vs 3.155 for non-DEGs). Colocalized genes are mostly high-expression genes (4.096 vs 2.909). RESULTS:1023-1028 itself names size and expression profile as the right reference. With expression matched, the observed 26 is at or above expectation. "Close to chance" holds under both nulls. The depletion implied by OR 0.757 and the dashed markers does not.
- **Re-derived (14,704 universe):**
  - AveExpr-stratified expectation: 22.66, 22.32, 22.13 and 22.28 for 5, 10, 20 and 50 bins.
  - Seeded permutation (B=10,000, seed 20261015): mean 22.17 to 22.59, two-sided p 0.418 to 0.498, Mantel-Haenszel OR 1.17 to 1.21.
  - SE-decile stratification: expected 27.39, p 0.824, MH OR 0.932. The finder's 27.6 did not reproduce.
  - Logistic model with an AveExpr spline: OR 1.201 (0.789 to 1.828), nominal Wald p 0.394.
  - The repo's existing expression-matched ladder: 26 genetic DEGs vs 23 matched controls, nominal McNemar p 0.775.
  - The three K10 universes give unstratified expectations of 33.426, 33.462 and 33.512, all with Fisher p 0.200.
- **New or known:** known in part (K10 covers the universes, not the null).
- **Action:** recompute. H57 corrects both the universe and the null in one test.

### 2.5 R04-02. One 18-donor dataset produces most of the pooled Fig 4C stage effect, and the homoskedastic SE is not valid (high)

- **Where:** Fig 4C beta bars and primary q; the Fig 4 paragraph in RESULTS.
- **Printed:** one pooled ordinal-stage beta per program from `score ~ stage + factor(dataset)`, with homoskedastic OLS SE as the prespecified primary uncertainty. Ductular injury beta is 0.402.
- **What the evidence supports:** the two datasets that carry stage information have different slopes and different residual variance, so one homoskedastic SE does not describe the uncertainty. The 37-donor GSE202379 slope is about 0.1, with p 0.0986.
- **Re-derived:** the 62-donor fit reproduces all 117 betas and q values (max |dq| 1.4e-14).

  | Program | GSE202379 slope (n=37) | GSE244832 slope (n=18) | Cochran Q p | IVW slope | DerSimonian-Laird |
  |---|---|---|---|---|---|
  | Ductular injury | 0.1074 (SE 0.0633, p 0.0986) | 0.7423 (SE 0.2568, p 0.0107) | 0.0164 | 0.1438 (SE 0.0615) | 0.376 (SE 0.314, p 0.23) |
  | Stromal ECM | 0.0373 (p 0.633) | 0.7556 | 0.0068 | 0.0985 (SE 0.0741, p 0.183) | 0.356 (SE 0.357, p 0.32) |

  Residual SD is 0.327 in GSE202379 and 0.950 in GSE244832 (Brown-Forsythe p 0.033). GSE244832 contributes 75.6% (Ductular injury) and 88.0% (Stromal ECM) of the pooled OLS slope numerator. Without GSE244832 (44 donors), 0 of 117 programs reach BH q<0.05, with minimum q 0.4251.
- **New or known:** known in general form. The combined review (:256, item 1.6; :259, item 2.3) said within-cohort support "needs careful interpretation" and gave no numbers. Not in known_issues.tsv.
- **Action:** recompute (H59).

### 2.6 R04-01. The only primary-selected program depends on two GSE189600 donors that ambient correction could not process (medium)

- **Where:** ABS:48-51 ("Donor-aware single-cell analysis and a prespecified 117-program registry resolved remodeling into multicellular programs ... Ductular injury alone met the primary 117-program rule"); the Fig 4C star; PAPER claim 2.
- **What the evidence supports:** Ductular injury is a stage association that passes BH over 117 programs only under OLS variance and with two uncorrectable donors included. No family-level disease result that survives these checks supports "resolved remodeling into multicellular programs".
- **Re-derived:**
  - Full 62-donor fit: beta 0.402086, OLS p 0.000322, BH117 q 0.037633, HC3 q 0.155826. One program has q<0.05 and none has HC3 q<0.05. The prespecified robust-display rule needs HC3 q<0.05 (robust_display_frozen is FALSE).
  - Without the 2 GSE189600 donors (60 donors): beta 0.384613, q 0.098523, 0 programs pass. This equals the complete-case raw q in `ambient_program_effects.tsv`. The corrected q there is 0.111792. In the full passthrough universe the corrected q is 0.059262.
  - Single-donor deletions: q ≥ 0.05 in 12 of 62 (the 2 GSE189600 donors, 8 GSE244832 donors and 2 GSE202379 donors).
  - The two GSE189600 donors are SRR22677587 (Healthy) and SRR22677595 (Steatohepatitis). Cell Ranger median genes per cell are 42 and 62. They contribute 43 and 99 hepatocytes, their decontX status is uncorrected_passthrough, and they carry the highest leverage (h = 0.531).
  - Seeded within-dataset permutation null:

    | Draws, seed | OLS-t BH117 q | OLS-t maxT FWER | HC3-t BH117 q | HC3-t maxT FWER |
    |---|---|---|---|---|
    | 19,999, seed 20261015 | 0.1053 | 0.02195 | 0.0351 | 0.1517 |
    | 50,000 | 0.101 | 0.0355 | 0.030 | 0.161 |

    The finder's values (OLS q 0.082, FWER 0.0365) did not reproduce exactly. The direction is the same.
- **New or known:** known for HC3 sensitivity only (combined review :109, older 64-donor state). The dependence on dataset and donor deletion is new.
- **Action:** reword the abstract and the Fig 4C legend; H59 refits. Related unverified items: G1-01 (the call is lost under the project's own 50-cell gate) and R13-04.

### 2.7 R06-02. "Two equally robust disease programs" contradicts the corrected Fig 4C (high)

- **Where:** fig5.md:41 and 45-46; PAPER:336-337.
- **Printed:** "the two prespecified robust programs"; "Two equally robust disease programs therefore differ"; "The two equally prespecified, robust disease programs".
- **What the evidence supports:** since the 2026-10-14 correction, LEG:235-237 says only Ductular injury passes, and fig5_discussion.md:44-46 says the same. The program with spatial support (Stromal ECM) lacks disease-association support. The program with disease-association support (Ductular injury) lacks spatial support.
- **Re-derived (independent 62-donor refit):**

  | Program | Beta | p | BH q | HC3 BH q | Registry state |
  |---|---|---|---|---|---|
  | Stromal ECM | 0.351525 | 0.002191 | 0.085458 | 0.165240 | tested_nonsignificant |
  | Ductular injury | 0.402086 | 0.000322 | 0.037633 | 0.155826 | selected_hc3_fragile |

  Both programs are multi_dataset_supported (2 cohort slopes, 5 leave-one-dataset-out fits with the same direction).
- **New or known:** new.
- **Action:** reword fig5.md and PAPER now. H59 recomputes both programs.

### 2.8 R07-01. The abstract says the Catalog links regulatory DNA elements and tissue context; the r8 Catalog has none at gene level (high)

- **Where:** abstract, introduction, Fig 6A and 6E, fig5_discussion:65-66.
- **Printed:** "The MASLD Gene Catalog connects protein-coding genes, lncRNA transcripts, regulatory DNA elements, and unresolved loci to their genetic, disease-state, cell, and tissue contexts." Fig 6A prints "Molecular object: Gene / RNA / DNA element".
- **What the evidence supports:** gene-level evidence covers only genetics and bulk RNA. A reader who opens any gene finds no chromatin, spatial, protein or credible-set link. That includes GNMT, whose liver protein change the paper reports.
- **Re-derived:**
  - The 28,233 records are protein_coding_gene 14,176, untestable_nomination 13,576, lncrna_transcript 423, ambiguous_lncrna_locus 58 and regulatory DNA element 0.
  - All 28,233 have noncoding_dna_context = unresolved, credible_set_to_gene_link_status = no_authoritative_link, an empty credible_set_id and a null coding PIP fraction.
  - Applicability is indeterminate for chromatin and spatial in all 28,233. Proteomics is indeterminate in 21,305 and not_applicable in 6,928.
  - GNMT's proteomics applicability is "indeterminate", although protein_de_adjusted.tsv gives logFC -1.13581 and padj 0.0064196 (6,071 proteins).
  - bulk_rna applicability is "indeterminate" for all 1,853 RNA-supported genes (HKDC1 is one), which contradicts rna_call = supported in the same records.
  - Evidence rows: 56,466 = 28,233 transcriptomics + 28,233 genetics.
  - Spatial, proteomics and chromatin are populated only at program grain. Program context rows cover only Stromal ECM and Ductular injury (14 rows).
- **New or known:** new.
- **Action:** reword the abstract and intro to what the Catalog holds, or add the layers. H53 would add a counted regulatory-region field. Section 6 lists per-gene results the paper already has that could be added. Fix the bulk_rna applicability contradiction in either case.

### 2.9 R05-02. Fig 4H "fibrosis-stage coefficient" codes 68 GSE213621 healthy controls as stage 0 (high)

- **Where:** Fig 4H printed rho and same-direction count; x-axis label; S4J; METHODS:372-373.
- **Printed:** the panel reads "Spearman rho = 0.73; 37/43 same direction", with x axis "Stage coefficient (score SD / cohort stage-code SD)". LEG:271-274 says "Fibrosis-stage coefficient" over "435 stage-complete units (109 GSE162694, 326 GSE213621)". METHODS:372-373 says "No step producing a Figure 4G or 4H value reads a control label or disease-versus-control result".
- **What the evidence supports:**
  - In the five-cohort manifest, GSE213621 stage 0 is exactly the 68 Control units. Disease units are coded 1, 2 and 3 (87, 94 and 77).
  - participant_histology gives these controls fibrosis_stage_exact = NA, so the code 0 comes from the control label.
  - GSE162694's 31 controls have NA stage and are excluded, so the two cohorts are not treated the same way.
  - Under the stated estimand, rho is 0.643 and 31/43 systems have the same direction. The METHODS:372 statement is false.
- **Re-derived:**
  - The published values reproduce: rho 0.726518, 37/43, median system stage coefficient 0.10984.
  - Using only the 258 GSE213621 disease units, with the published continuum coefficients kept:
    - 43 systems: rho 0.643008, 31/43 same direction, median stage coefficient 0.041539.
    - Re-standardized variant: rho 0.6571, 31/43.
    - 116 Hotspot programs (S4J): rho 0.8228 to 0.6598; same direction 93 to 85.
    - Donor-level system models: stage-supported systems 27 to 21; rho against the continuum 0.7531 to 0.6797; same direction 39 to 34.
- **Not recomputed:** the continuum coefficients themselves. They also include the GSE213621 controls, as their own factor level.
- **New or known:** new.
- **Action:** recompute 4H, S4J and the x-axis label; correct METHODS:372. The same stage coding enters Fig 3D (2.12) and possibly the GSE213621 rows of the stage-trend table that feeds the 11,130 stage-supported genes [Inference, flagged in the H51 plan].

### 2.10 R03-02. The GSE126848 per-cohort model behind 6,147 Fig 3C calls includes 437 all-zero genes, and its contrast is aliased with depth (medium)

- **Where:** Fig 3C GSE126848 ribbon (6,147 calls, 3,131 unique to the study).
- **Printed:** the legend presents GSE126848 (23 control, 7 disease clusters) as a restored per-cohort model. It notes only that unequal eligibility limits representativeness.
- **What the evidence supports:** the per-cohort calls are not a usable disease contrast as printed.
  - The identity screen removed 24 of 31 disease libraries and 2 of 25 controls (Fisher p 1.33e-7; sample OR 39.4, conditional MLE 35.85).
  - In the retained 30, disease libraries have median 19.48M fragments against 35.30M for controls (AUC 0.031, exact Wilcoxon p 1.87e-5).
  - Group and log2 library size correlate at r = -0.709.
  - Verification found that 437 of the 6,147 calls are genes with zero reads in all 30 libraries. Each has logFC 0.554 and FDR 0.016. The per-cohort fit runs on all 24,115 pooled genes with no per-cohort expression filter.
- **Re-derived:**
  - 1,873 calls have fewer than 30 reads across the 30 libraries.
  - The lowest-expression quartile holds 2,285 calls, 99.8% of them up.
  - A per-cohort filterByExpr removes 2,730 calls. On the remaining 15,592 genes, the refit gives 3,446 calls (1,806 up, 1,640 down) instead of 6,147 (4,385 up, 1,762 down).
  - Removing only the 437 all-zero calls changes the 3C union from 13,761 to 13,447, unique-to-study from 9,065 to 8,843, and memberships from 20,054 to 19,617. 314 of the 3,131 GSE126848-unique genes are all-zero.
  - A depth-adjusted model is estimable (VIF 2.01). It gives 2,678 calls and keeps 2,025 of the 6,147 (32.9%).
  - The pooled 1,854 calls are little affected: leaving GSE126848 out keeps 96.8%.
- **New or known:** K14 records the screen imbalance only. The depth aliasing and the all-zero calls are new.
- **Action:** refit GSE126848 with a per-cohort expression filter, then recompute Fig 3C and reword the legend. Sub-claims that did not reproduce are in section 3.2.

### 2.11 R03-01. Fig 3C overlap tiers count opposite-direction calls as shared (medium)

- **Where:** Fig 3C overlap tiers; PAPER:392 ("a broad established-disease signal recurs across cohorts").
- **Printed:** union 9,065 / 3,462 / 922 / 261 / 51 genes called in one to five cohorts.
- **What the evidence supports:** the printed counts reproduce, but the tiers ignore direction.
  - Of the 4,696 genes called in two or more cohorts, 1,865 (39.7%) are up in one cohort and down in another. By tier: 1,529 of 3,462, 286 of 922, 43 of 261 and 7 of 51.
  - 1,632 of the 1,865 involve a GSE126848 call. GSE126848 and GSE135251 share 1,936 calls, and only 690 of them (35.6%) have the same sign.
  - The all-gene logFC Spearman between GSE126848 and the other cohorts is 0.125, -0.119, 0.012 and 0.179. Among the other four cohorts it ranges from 0.244 to 0.602.
  - The 7 discordant all-five genes are ABCB11, CD83, RGS16, JUN, S100A4, LINC01717 and ENSG00000279218.
- **Re-derived:** same-direction tiers are 2,223 / 675 / 225 / 44, or 3,167 genes in two or more cohorts instead of 4,696. Without GSE126848, 346 of 2,565 multi-cohort genes are discordant (13.5%).
- **New or known:** known (docs/HIGH_IMPACT_PATTERNS.md:53). Not corrected on the panel.
- **Action:** recompute the tiers by direction, after the 2.10 refit.

### 2.12 R03-03. Fig 3E and 3F use 387 units from three cohorts, and the legend implies the 774-unit, five-cohort roster (medium)

- **Where:** LEG:113 ("All panels use the corrected 774-unit substrate"); Fig 3A GSE213621 fibrosis cell; Fig 3D; fig3.md:58 ("eligible cohorts").
- **What the evidence supports:**
  - The 3E and 3F fits are hard-coded to GSE130970, GSE135251 and GSE162694 (`build_pi_stage_extensions.R:226-229`). That gives 387 units, with F0/F1/F2/F3/F4 = 103/102/88/69/25. GSE135251 supplies 207 (53.5%).
  - GSE213621 contributes 0 units to every stage contrast, although it carries 45.12% of the quality weight and has 258 staged disease units. Its source stages are binned.
  - The 774 metadata store those bins as integers: Control as 0 (68), F0-F1 as 1 (87), F2 as 2 (94), F3-F4 as 3 (77). That is why 3A shows the GSE213621 fibrosis cell as "Available".
  - 3D colors the 77 GSE213621 F3-F4 units as stage 3, and puts the 68 controls in "strict control" with stage 0.
  - The F0 reference is 14 source controls plus 89 disease units. GSE162694's F0 is 35 disease units only, because its 31 controls have no recorded stage.
  - Neither the 3E nor the 3F PDF prints a cohort name.
- **Re-derived:** 3E counts reproduce as 1,338 / 4,266 / 6,190 / 7,043 of 24,115 genes. n_cohorts = 3 for all fibrosis contrasts.
- **New or known:** new as a legend problem. The archive noted GSE213621's absence from the stage substrate (04_CANDIDATE_POOL.md:1584), but no legend says so.
- **Action:** fix the legends now: name the three cohorts and n = 387, and mark the GSE213621 3A cell as binned. H51 tests replication in the GSE213621 bins. Related unverified items G2-01 to G2-03 (section 3.3) say the F1-vs-F0 count of 1,338 rests mostly on one cohort.

### 2.13 R02-02. "Unresolved regulatory context" mostly records failed coordinate lookup in non-EUR studies (medium)

- **Where:** the noncoding-context sentence in Results; the noncoding-context companion panel; RESULTS:2540-2542.
- **Printed:** "36.49% of direct-trait and 77.85% of liver-enzyme noncoding posterior probability lacks resolved deposited regulatory context". RESULTS gives the unresolved difference as -0.413643 and the promoter difference as 0.138659. The direct ABC fraction is 0.126%. The credible-set-weighted panel prints unresolved 52.11 / 60.43, TSS 7.87 / 4.46, open 7.32 / 4.15 and ABC 0.23 / 1.57.
- **What the evidence supports:** "unresolved" is a processing failure that depends on the study, not missing deposited context.
  - Rows labelled unresolved_liftover: EUR 0.74%, AMR 98.39%, AFR 85.68%, EAS 61.28%, SAS 58.67%.
  - The hg38 lookup table holds 6,912 variants from 22 studies, with no MVP study and no PanUKBB_CSA_ALT. 25,944 of the 26,442 unresolved gated rows sit at sites absent from that table. The other 498 have swapped allele order.
  - liftOver of all 18,401 member sites fails for none.
  - TM6SF2 rs58542926 (19:19379549) is resolved in four EUR studies and unresolved in MVP_ALT_AFR and MVP_NAFLD_AMR.
  - ABC overlap is recorded as False, not missing. AMR has 0 of 9,210 True.
  - Consequence (VEP) annotation is complete in every ancestry, so the 2D coding fractions are unaffected.
- **Re-derived:** re-annotation by position matches the recorded values on 100% of the 10,670 already-resolved rows. Cluster-weighted values, direct vs enzyme:

  | Context | Adopted | Adopted difference [95%] | Position-repaired | Repaired difference [95%] |
  |---|---|---|---|---|
  | Unresolved | 36.49% vs 77.85% | -0.4136 [-0.616, -0.199] | 0% vs 0% | n/a |
  | Promoter | 16.93% vs 3.06% | 0.1387 [0.018, 0.292] | 22.53% vs 12.58% | 0.0995 [-0.025, 0.249] |
  | Accessible | 10.48% vs 4.52% | +0.0596 | 11.92% vs 15.04% | -0.0312 [-0.103, 0.050] |
  | Resolved, no context | 38.42% vs 15.81% | +0.2261 [0.076, 0.379] | 68.47% vs 77.43% | -0.0896 [-0.234, 0.038] |
  | ABC | 0.126% vs 0.707% | | 0.536% vs 3.913% | -0.0338 [-0.052, -0.018] |

  Repaired panel values: unresolved 0 / 0, TSS 13.84 / 14.68, open 9.26 / 8.44, ABC 0.43 / 3.70. The repaired values are still pooled across ancestries, so 2.2 applies to them too.
- **New or known:** new.
- **Action:** recompute with position-based annotation (H58 family X) and reword the sentence. One sub-claim is unresolved (section 3.2).

### 2.14 R02-04. The Fig 2B x axis takes a maximum over four fine-mapping methods and a different GWAS and lead than y (medium)

- **Where:** Fig 2B point positions and corner populations.
- **Printed:** the producer header says x is "within-ancestry ... NOT the contested cross-ancestry SuSiEx". There are 695 plotted points, and RORA is labelled at x = 1.
- **What the evidence supports:**
  - x is `pmax(recommended_pip, susie_pip)`, where recommended_pip is already a maximum over SuSiE and CARMA, and is further raised to SuSiEx and MESuSiE.
  - 62 of the 695 x values come from a cross-ancestry joint PIP.
  - 204 of the 695 are SuSiE PIPs from non-EUR strata. The aggregator had replaced these with CARMA because of 1KG-LD ghost credible sets, and the producer's extra pmax brings them back.
  - x and y come from different GWAS for 465 of 695 points, and from different lead SNPs for 359.
  - The visual count of "well-localized variants that fail to colocalize" therefore depends mostly on pairing x and y from different sources.
- **Re-derived:** an exact reproduction (2,547 rows, 695 plotted, max |dx| = 0) gives 152 points in the corner (PIP ≥ 0.9, PP.H4 < 0.5). Alternatives:

  | x definition | Points that change | Points plotted | Corner |
  |---|---|---|---|
  | SuSiE only, gene-level max | 223 | 662 | 116 |
  | SuSiE at y's own GWAS and lead | 559 | 327 | 18 |
  | Producer's max-method at y's own GWAS and lead | | 343 | 28 |

  RORA's x is 0.9999 (UKBB_GGT, 15:60883281). Its y is PP.H4 0.9977 (BBJ_GGT, 15:60878030). The SuSiE PIP at the y lead is 0.4402.
- **New or known:** K05 and the combined review :195-208 record the maximum over methods and say the count was never recomputed. The x/y source mismatch and the re-use of non-EUR SuSiE are new.
- **Action:** recompute with one method, and with x and y from the same GWAS and lead.

### 2.15 R02-03. The Fig 2A cascade implies the 996 COLOC genes come from the 360 fine-mapped loci (medium)

- **Where:** Fig 2A legend and panel; the Fig 2D denominator sentence.
- **Printed:** "followed by 360 trait-specific fine-mapped loci (265 globally unique loci) and 996 COLOC genes". The panel draws a fine-mapped loci column feeding a COLOC genes column. The 2D text reads "all reliable Tier-1/2 SuSiE credible-set instances".
- **What the evidence supports:** most COLOC genes, and most EUR enzyme genes, have no fine-mapped credible set in this release.
  - MVP_ALT_EUR, MVP_AST_EUR and MVP_NAFLD_EUR have no credible sets.
  - UKBB_ALT, UKBB_AST and UKBB_GGT were fine-mapped at only 6, 9 and 15 GWAS Catalog leads, although the current lead files list 564, 832 and 1,113 leads [Unverified cause].
  - The panel's EUR ALT row has 4 loci but 124 multi-signal genes.
  - The 2D direct arm has 13 studies and lacks MVP_NAFLD_EUR and NAFLD_Intermountain_EUR. The enzyme arm has 18 and lacks MVP_ALT_EUR and MVP_AST_EUR. The EUR enzyme arm is 82 instances (UKBB_GGT 61, AST 15, ALT 6).
- **Re-derived:** 360 and 265 reproduce. Of the 436 multi-signal genes:
  - 100 have a positive top SNP within ±250 kb of a same-study window (89 without padding).
  - 190 have one within any window (147 without padding).
  - 94 are supported only by MVP EUR strata.
  - 330 are supported only by MVP EUR plus UKBB enzyme strata.
- **New or known:** known (archive 01_HYPOTHESIS_INVENTORY.md:3150). Not corrected.
- **Action:** reword the cascade legend and the 2D denominator sentence. H58 part 3 produces the support count.

### 2.16 R01-03. The Fig 5E program roster was selected on the superseded August single-cell fit (medium)

- **Where:** Fig 5E.
- **Printed:** 13 rows (Hep 17/6/14/8/24; Fib 6/23/4; Mac 7/14/15; Chol 20/15). The legend (LEG:408-409) reads only "(E) Assay-native program evidence with complete-family correction and testability states."
- **What the evidence supports:** row membership rests on disease calls that no longer exist, and the one corrected disease program is absent.
  - The producer reads a 22-program universe frozen on 2026-08-06. All 22 had August q < 0.05 (max 0.0432).
  - Under the 2026-10-14 fit, 1 of the 22 has q < 0.05: Ductular injury (hepatocytes::20, q 0.037633). It is not displayed, because it has 0 open-promoter genes, untestable histology and no robust spatial result.
  - 0 of the 13 displayed rows has q < 0.05. Two have q < 0.10 (Stromal ECM and NRF2, both 0.085458). All 13 keep their August sign.
  - Three rows appear only because of the histology test (BH family of 17): Stress padj 0.02975, NRF2 0.0153 and Mac Adhesion/Rho 0.0323.
  - The open-promoter n column comes from `GWAS/finemapping/results/credible_sets.csv` (5 July; in-CS or PIP > 0.1, 35 Tier-1/2 studies). That file lacks the three MVP EUR studies that Fig 5D uses (see unverified R11-03). Fig 5D prints the COLOC r2 values 28 and 739.
- **New or known:** new.
- **Action:** recompute the roster under the corrected fit, and state the selection rule in the legend.

### 2.17 R07-02. Zero tested_negative calls is fixed by the Catalog rules but presented as a result (medium)

- **Where:** Fig 6B tested-negative row and column; PAPER:316; RESULTS:2699; the Fig 6 README ("a central result"); fig5_discussion:24-26 ("a valid Resource result").
- **What the evidence supports:** no negative rule was ever defined, so the zero could not have been otherwise. Genes where COLOC favors distinct causal variants carry the same label as genes with no signal.
- **Re-derived:**
  - negative_call_rule_id, decision_boundary, margin and passed are NA in all 56,466 evidence rows.
  - Genetics calls: 18,478 indeterminate, 436 supported, 9,319 untestable. Transcriptomics: 22,148 indeterminate, 1,853 supported, 4,232 untestable. tested_negative is 0 in both.
  - Fig 6B prints 7 zero cells that cannot be nonzero.
  - Data for a rule exist. In the signal-pair file (745 genes), 21 genes have an eligible SuSiE pair with PP.H3 > 0.8 and max PP.H4 ≤ 0.5, and 14 have PP.H3 > 0.9 and PP.H4 < 0.1 (VPS29, SLC17A1, PDHB and JAZF1 among them). The Catalog labels 19 of the 21 indeterminate, 1 untestable (CKB), and 1 is absent (ARL17B).
  - Caveats: 6 of the 21 have tier-1/2 ABF PP4 > 0.5. At gene-study level in the master table, 4,081 genes meet PP.H3 > 0.8 with max PP.H4 ≤ 0.5.
- **New or known:** related to archive candidate resource_design-07 (04_CANDIDATE_POOL.md:88), which called the zero a defect but did not record that it is fixed by construction.
- **Action:** reword: "no negative rule was defined; tested_negative was not assessed". A negative rule would be a new analysis. The bulk version (L09-01) failed the gate (section 5).

### 2.18 R08-01. Fig 7B: no pretrained sequence model was scored on severity, and 1 of 7 MPRA rows survives correction (medium; candidate figure)

- **Where:** legend 7B; draft fig7.md:23-27.
- **Printed:** pretrained sequence models "improve on the task-native baseline for variant-to-regulatory-activity ... and do not improve on it for tissue-level severity"; "Neither half alone is the result."
- **What the evidence supports:** the severity half of the claim has no data. The sequence half rests on one of seven rows, from one synthetic reporter assay. On the endogenous sequence task, the two pretrained embeddings score below distance alone.
- **Re-derived:**
  - The severity rows are "composition + lineage-residualised expression", +0.1733 [0.0785, 0.2683], and "routed multi-modal system", 0.0000. Neither is a pretrained sequence model.
  - MPRA rows (7 rows, 6 models):

    | Model | Gain [95%] |
    |---|---|
    | caduceus | +0.0107 [-0.0408, 0.0624] |
    | dnabert2 | -0.0378 [-0.1183, 0.0413] |
    | HyenaDNA | +0.0820 [0.0274, 0.1389] |
    | Nucleotide Transformer | +0.0537 [0.0011, 0.1108] |
    | Enformer SAD | -0.0512 [-0.1200, 0.0161] |
    | Enformer SAR | -0.0492 [-0.1140, 0.0141] |
    | sei | +0.0596 [-0.0035, 0.1195] |

  - With Bonferroni across the 7 rows, only HyenaDNA stays above 0 (lower bound +0.0055 by normal approximation; percentile [0.0089, 0.1585]).
  - Endogenous variant-to-gene |z|: HyenaDNA -0.0218 [-0.0329, -0.0100], NT -0.0294 [-0.0421, -0.0162].
  - The only pretrained-model severity comparison found is BulkFormer on GSE268273: -0.0006 [-0.0662, +0.0649], below the MDE of 0.0667 (n=109). All 109 samples are BulkFormer pretraining rows.
- **New or known:** new. Related unverified R09-01 (the only baseline run was an allele-identity ridge).
- **Action:** reword to what was measured.

### 2.19 R08-02. 96 of the 177 GSE276114 records counted among the 585 participants do not have MASLD (medium; candidate figure)

- **Where:** closing paragraph of the Fig 7 legend; RESULTS:3175 and 3353; METHODS:836-840.
- **Printed:** "585 participants in three cohorts it never trained on"; the instrument reaches "0.68 to 0.79" of the tie ceiling.
- **What the evidence supports:** GSE276114's 177 records are 81 MASLD, 82 chronic viral hepatitis (CVH) and 14 alcohol-related liver disease, all 14 at F4. The 0.79 at the top of the range comes from a mix of etiologies.
- **Re-derived:**

  | Group | n | rho | Tie ceiling | Ratio |
  |---|---|---|---|---|
  | Pooled | 177 | 0.673865 | 0.848327 | 0.79435 |
  | MASLD | 81 | 0.785980 | 0.907474 | 0.86612 |
  | CVH | 82 | 0.544323 | 0.799892 | 0.68050 |
  | All non-MASLD | 96 | 0.515797 | 0.757760 | 0.68069 |

  The 96 non-MASLD records are 16.4% of the 585.
- **New or known:** new.
- **Action:** report the composition, or restrict the MASLD transfer claim to the 81 MASLD participants.

### 2.20 R08-03. GSE213621's 299 "participants" are 274 genetic identities, and 11 identities carry conflicting stages (medium; candidate figure)

- **Where:** Fig 7 legend (585 participants); RESULTS (GSE213621 299 and 367 records).
- **What the evidence supports:** the Resource's own identity screen covers all 367 GSE213621 libraries (342 components, 0 spanning datasets).
  - The 299 staged records collapse to 274 clusters, 25 of them pairs.
  - 11 pairs carry different fibrosis categories: F2/F3F4 6, F0F1/F2 4, F0F1/F3F4 1.
  - The three cohorts therefore contain at most 109 + 274 + 177 = 560 distinct people.
  - A bootstrap that resamples the 299 records as independent people is too narrow [Inference].
- **Re-derived:** the published values reproduce (rho 0.638198, ratio 0.677263). Keeping one record per identity and dropping the conflicting clusters gives n = 263, rho 0.667979, ratio 0.709. Other ways of picking one record per cluster give ratios of 0.6906 to 0.7204.
- **New or known:** known (STATUS:1057, which does not cover the Figure 7 transfer counts).
- **Action:** recompute at identity level and bootstrap identities. Not run: a genotype check of GSE268273 and GSE276114 against the training libraries.

## 2B. Second verification pass: the 20 remaining high-severity findings

The main pass verified only the 20 findings ranked first. The 20 other findings the finders rated high were checked in a second pass (run `wf_00cf4a6c-d4e`). Each got one blind reproduce verifier, who was given only the claim text and the locations. Two false decoys were hidden among them under ordinary-looking IDs (G4-01, G4-02), and both were refuted. Results: 13 confirmed, 7 partly confirmed, 0 refuted. Severity below is the verifier's own rating, not the finder's. Scripts and outputs are in `RUN/scratch/verify2-<ID>/`, and the 10 sbatch jobs are in `RUN/jobs/40x1` to `RUN/jobs/42x1`. Full verdicts are in `RUN/verify2_result.json`.

Order: claims in the title or abstract first, then main-panel numbers, then the Gene Catalog, then text, then the Figure 7 candidate.

| # | ID | Verdict | Severity | Where | What the evidence supports |
|---|---|---|---|---|---|
| 1 | R13-03 | Confirmed | high | Fig 4E, LEG, fig4.md | On the 774-unit substrate alone (364 units, 3 cohorts), Fig 4E has 0 fibrosis and 0 NAS calls (smallest q 0.079 and 0.159). The printed 16 fibrosis and 1 NAS calls, and 3 and 0 after composition adjustment, need the 93 identity-unresolved GSE174478 records. Those records add a fourth cohort, so the Knapp-Hartung df goes from 2 to 3. Both arms came from one run, and the four-cohort arm was chosen on 2026-10-13 after both had been computed. |
| 2 | R12-01 | Confirmed | high | Fig 4E | All 20 Fig 4E calls have tau2 = 0 and Q of 0.09 to 1.6, against an expected 3. Under these conditions untruncated Knapp-Hartung shrinks the SE by up to 5.8-fold. With the truncated variant (metafor `test="adhoc"`), the counts are 0/0/0/0. A plain REML z-test gives more calls (47/28/21/10), so untruncated HK is not the most permissive rule. The call set depends on which accepted rule is used when only 4 cohorts are available. |
| 3 | R13-02 | Confirmed | high | Title clause 1, abstract 24-27, intro 77-80 | The 23.34% coding plus splice share for direct traits averages over 19 one-Mb regions. Pooled over credible-set instances it is 15.81% vs 3.26% (ratio 4.85). The direction holds under both weightings and all 57 sensitivity rows (ratios 4.37 to 8.28). 75% (direct) and 94% (enzyme) of noncoding probability has no deposited regulatory context. "Resolves the architecture" overstates the result. What is resolved is the coding versus noncoding split, and the abstract has to state the weighting. |
| 4 | R14-02 | Confirmed | high | Fig 2D, abstract, fig2.md | Most of the direct-trait enrichment comes from the three EUR PDFF imaging GWAS. They hold 37 of 207 sets but 12 of 19 direct clusters, with a 40.1% coding fraction and a ratio of 10.36 [5.04, 22.24]. The NAFLD and NASH GWAS alone (10 studies, 170 sets) give 12.4% and a ratio of 3.21 [0.53, 8.25]. The difference interval includes 0 in all 6 seeds tried. "Direct MASLD and liver-fat" mainly describes liver fat. |
| 5 | G2-03 | Confirmed | high | Fig 3F | 29 of the 54 genes are on the panel only because of an F0-to-F1 mark. Without GSE130970, none of those 30 marks stays at BH q<0.05. Under joint complexity covariates only 1 does (FAM177A1). Reselecting under the joint adjustment gives 41 genes: 22 stay, 32 leave (including HKDC1, SPP1, MTOR) and 19 enter (including SLC2A2, ARNT). |
| 6 | G2-01 | Confirmed | medium | Fig 3E F1 bar | Of the 1,338 F1-vs-F0 calls, 14 remain without GSE130970 (155 units). Without GSE135251 there are 929, and without GSE162694 there are 651. Fitted alone, the cohorts give 252 (GSE130970, 50 units), 3 (GSE135251) and 0 (GSE162694) calls. 88 to 91% of the effects have the same sign in the other cohorts, so an early shift exists, but its size at BH q<0.05 comes from one cohort. GNMT and MAT1A lose F0-to-F1 support without GSE130970 (q 0.207 and 0.089). |
| 7 | G2-02 | Confirmed | medium | Fig 3E F1 bar | The 1,338 calls fall to 280 when adjusted for the TMM norm factor, 485 for Gini, and 104 to 188 for library size, Gini and norm factor jointly. Library size alone leaves 1,350. Surviving effects are about 23% smaller (median ratio 0.773). The F2 to F4 bars are less affected: reviewer job 1721 kept 50 to 78% of those calls, but this was not re-run. |
| 8 | R15-01 | Confirmed | high | Fig 6C, Gene Catalog | The Catalog genetics call measures only liver-eQTL colocalization (PP.H4 > 0.5). So MTARC1, PNPLA3, TM6SF2 and HSD17B13 read as "indeterminate / Unresolved". In the Resource's own fine-mapping (Fig 2D), the MTARC1 missense variant has PIP 0.757. All 28,233 records have an empty credible_set_id and a null coding PIP. The column is labelled only "Genetics". |
| 9 | R15-02 | Confirmed | medium | Catalog PP.H4 column | One PP.H4 column holds SuSiE values for the 436 supported genes and ABF values for 15,231 indeterminate rows, with no method field. 555 indeterminate rows show ABF > 0.5, and 152 of them show ≥ 0.8. GNMT shows 0.706 (ABF) while its SuSiE value is 0.128. No class count changes. |
| 10 | G1-01 | Partly | high | Abstract, Fig 4C star | Ductular injury (q 0.0376, the only primary call) is lost under a 50-cell lineage gate (q 0.0765), when GSE189600 is dropped (q 0.0985), in the 60-donor common universe, and after 12 of 62 single-donor deletions. Not supported: the finder said "the project's own 50-cell gate" applies. METHODS:21-27 makes that gate primary only for the two targeted lineage tests, so the call passes the prespecified rule. The fragility is also broader than GSE189600: dropping GSE244832 gives q 0.259 and dropping GSE202379 gives q 0.118. |
| 11 | G3-01 | Partly | high | Fig 5E vs 5F | Fig 5E prints Stromal ECM spatial excess of 0.1436 (Guilliams) and 0.1405 (Vu). Fig 5F prints 0.0405 and 0.0503. 5E adjusts only for the program's own lineage, 5F adjusts for all 16 cell2location factors, and neither panel says so. The adjustment causes most of the drop (0.199 to 0.059); the membership change causes little. 6 of the 7 programs with a 5E bar have never been fitted with all-cell adjustment. Not supported: the finder said the other programs were "hepatocyte-only" (each is adjusted for its own lineage) and gave "about 71%" for both cohorts (Vu is 64 to 66%). |
| 12 | R11-01 | Confirmed | high | PAPER:187-189, STATUS:1116-1118 | The NAS-component separability state "tested_negative" cannot be reached through the decision rule. C1 and C2 fail on histology label correlations before any expression is loaded (participation ratio 1.42 from the 3x3 grade matrix). In simulation, the molecular statistic gives the same range whether the components share one axis or have separate programs (1.09 to 1.43 vs 1.10 to 1.44). The state should be "not resolvable with this design". |
| 13 | R11-02 | Confirmed | medium | Fig 5D | "Partial mass" covers 3,409 of 3,835 lineage rows because any mass above 0 counts. In 1,943 of those rows the mass is below 0.01, and in 1,422 it is below 1e-6. 239 of 767 pairs, including 13 of 28 direct pairs, have under 0.01 in every lineage. No background or expected overlap is computed, yet LEG:407 says peak overlap "localizes regulatory DNA". |
| 14 | R13-04 | Partly | medium | Abstract 48-51, intro 55-56 | "Resolved remodeling into multicellular programs" overstates 1/117 primary calls (HC3-fragile), 0/117 HC3 calls and 37/117 lineage mismatches. The abstract's own numbers are correct. Not supported: the finder's "0/224 in bulk Fig 4E". The Fig 4E table gives 17/224, or 3/224 after composition adjustment. Finding 1 explains why those 17 depend on GSE174478. |
| 15 | R13-05 | Partly | medium | fig4.md, Discussion | The main text does not say that composition adjustment lowers Stromal ECM in GSE162694 from 0.888 to 0.172 (HC3 q 0.23), or that the cross-cohort estimate falls from 0.760 to 0.294. The S3-S4 legend already reports these values, so only the main text omits them. Within-stage signal survives in every estimable arm. Fix: one sentence plus the effective number of independent systems (14). |
| 16 | R14-03 | Partly | medium | Abstract "423 of 436" | The 423/13 split mainly reflects how many signals were tested: 6,208 vs 207 credible sets, and 15.7-fold more SuSiE-tested gene-studies. The positive rate is also lower for direct traits: 28/771 (3.6%) vs 739/12,121 (6.1%), OR 0.58, nominal p 0.0036. "Only 13" hides that 22 genes reach direct-trait PP.H4 > 0.5. Not supported: the finder's "not lower yield" rested on a table that includes only gene-studies already positive. |
| 17 | R14-04 | Partly | medium | Intro 37-42, Discussion 213-222 | Four relevant papers are not cited. Piras and DiStefano 2024 (10 cohorts, 1,058 biopsies, integrated with GWAS by MAGMA and TWAS) contradicts the Discussion's claim that earlier bulk meta-analyses lack an inherited layer. The other three are Gobeil 2026 JCEM, Saliba-Gustafsson 2024 and Jamialahmadi 2024. The competitor file also records a wrong DOI for Piras; the correct one is 10.26508/lsa.202302517. The core novelty sentence survives if narrowed, because no single paper combines all three elements. |
| 18 | R08-05 | Confirmed | high | Fig 7F (candidate) | Geneformer has 0 clean held-out studies and an "uninterpretable" label. It is still drawn in the magenta "above zero" style, and it is the largest two-layer gain on the panel (+0.0173 vs +0.0124 for the next). This contradicts Legend A's rule. |
| 19 | R09-01 | Confirmed | high | Fig 7B (candidate), PAPER idea 8 | The only MPRA baseline was a ridge on allele identity (0.134). Every run receipt records `mandatory_task_native_baselines_complete = False`, with deltaSVM, gkm-SVM, sequence CNN, sequence transformer and MPRALegNet missing. 2 of the 7 models have intervals above 0. |
| 20 | R08-04 | Partly | medium | Fig 7E (candidate) | The treatment-response dot at 0.592 divides the diet-arm difference by the pooled MDE (0.100). With the diet arm's own MDE (0.138) the ratio is 0.430, and oriented to the prespecified direction it is -0.43, which is the wrong direction. Not supported: the finder's "should be untestable". The prespecification calls a reversed direction negative, and the interval excludes a correct-direction effect the size of the MDE (0 of 9,979 draws). |

What the second pass adds to the five patterns in Section 1:

- **One cohort or one analysis choice produces the result.** Fig 3E's F1 bar and half of Fig 3F's genes depend on GSE130970 and on library complexity. Fig 4E's calls depend on adding GSE174478 and on the Knapp-Hartung variant.
- **The Catalog's genetic layer measures less than its label says.** "Genetics" means liver-eQTL colocalization, so the best-known coding risk genes read as unresolved.
- **Wording that the numbers do not support.** This covers "resolves the architecture", "direct MASLD" enrichment that mostly reflects liver fat, "resolved remodeling into multicellular programs", and the uncited prior work.

## 2C. Re-derivation by the main loop

I re-derived three confirmed findings with my own code from the producing files, without reusing any agent script. Code and outputs are in `RUN/scratch/main-rederive/`.

- **R01-01 (Fig 1B).** `r01_01.py`. Joined `gene_level_coloc_tier12.csv` (r2) to the 774-unit `deg_results.csv` on base Ensembl ID, with DEG = padj < 0.05 and |logFC| > 0.5. The panel reproduces: n 14,704, genetic 426, DEG 1,155, both 26, expected 33.462 (7.86% and 2.90%), OR 0.757, Fisher p 0.200. The SuSiE states are no_susie_result 10,928 (935 DEGs), tested_le_0.5 3,343 (194), positive 426 (26) and untestable 7 (0). Restricted to SuSiE-evaluated genes: n 3,769, DEG 220, both 26, expected 24.866 (5.84% and 11.30%), OR 1.055, p 0.826. The H57 plan quotes 3,776 genes, OR 1.057 and p 0.743 because it keeps the 7 untestable genes. Argument against the finding: a gene with no SuSiE result often has no GWAS signal nearby, and a referee could call that a real negative. The paper's own rule (fig2.md:27-29) says missing SuSiE output is an observability state, so the inconsistency stands. The fix is to report both denominators. "Near chance" holds under both.
- **R03-03 (Fig 3E/3F).** `stage_extension_all_gene_results.tsv` in `figures/candidates/fig3e-s3-774-20261011T165944Z/` gives F1, F2, F3 and F4 vs F0 counts of 1,338, 4,266, 6,190 and 7,043. Every contrast has n_cohorts 3, with F0/F1/F2/F3/F4 = 103/102/88/69/25 (387 units). `scripts/figures/build_pi_stage_extensions.R:226-227` restricts to GSE130970, GSE135251 and GSE162694. LEG:113 says "All panels use the corrected 774-unit substrate". In the five-cohort manifest GSE213621 is Control 68 at stage 0 and Disease 87/94/77 at stages 1/2/3.
- **R05-02 (Fig 4H, Hotspot part).** `r05_02_hotspot.py`. A fresh OLS of outcome_z on standardized stage plus sex in GSE162694 and GSE213621, combined by inverse-variance meta, reproduces the published 116 Hotspot stage betas (max difference 7.5e-16). Rho with the continuum beta is 0.8228, with 93/116 in the same direction. After the 68 GSE213621 controls are removed, rho is 0.6598 and 85/116 are in the same direction. These are the finder's numbers. My first attempt included all four staged cohorts and did not reproduce the published betas (max difference 0.088). The published model uses only the two cohorts with n_cohorts = 2. I did not re-derive the 43-system values (rho 0.7265 to 0.643, 37/43 to 31/43); the finder's job 1051 reproduced the published system values to 1e-15.

## 3. Refuted, unresolved and unverified findings

### 3.1 Refuted

- None of the 20 verified findings was refuted.
- Decoy X-02 (planted false): refuted, correct.
- Decoy X-03 (planted false): refuted, correct.
- Second-pass decoy G4-01 (planted false: the Fig 3E F3 bar is 6,910, not 6,190): refuted, correct. The source table gives 3,144 + 3,046 = 6,190.
- Second-pass decoy G4-02 (planted false: r8 has 426 genetically anchored genes): refuted, correct. r8 has 400 genetically anchored and 26 concordant.

### 3.2 Sub-claims the verifiers did not reproduce, or on which they disagreed

| Finding | Finder's sub-claim | What verification found | State |
|---|---|---|---|
| R03-02 | Depth cannot be adjusted out | The depth-adjusted model is estimable (VIF 2.01, 2,678 calls) | Refuted |
| R03-02 | Quality weights 0.72 vs 1.14 | These are pooled 774-fit weights. The per-cohort model gives 0.935 vs 0.970. | Corrected |
| R03-02 | Depth slope correlates +0.113 / +0.189 / +0.195 / +0.060 in the other cohorts | Verifier 1: +0.158 / +0.221 / +0.252 / +0.110. Verifier 2: -0.153 / +0.038 / -0.069 / +0.033. | Unresolved; not used |
| R02-02 | MVP studies run 0.93 to 1.00 unresolved | Verifier 1: MVP_ALT_EAS 0.262 and MVP_NAFLD_EAS 0.477. Verifier 2: 0.9306 to 1.000. | Unresolved; conclusion unaffected |
| R04-01 | Permutation OLS q 0.082, maxT FWER 0.0365 | q 0.101 to 0.1053, FWER 0.02195 to 0.0355 | Same direction, values corrected |
| R07-03 | SE-decile expectation 27.6 | 27.39 | Corrected |
| R02-01 | 16 of 19 direct clusters contain EUR | 15 in the pooled clustering. 16 comes from re-clustering the EUR subset. | Corrected |
| R02-04 | 39 MESuSiE, 23 SuSiEx | Verifier 2: 23 SuSiEx, 46 MESuSiE, 7 ties. Both give 62 in total. | Split unresolved, total agreed |
| R08-01 | Bonferroni lower bound +0.0055 | This is the normal approximation. The percentile interval gives +0.0089. | Both above 0 |
| R01-03 | Open-promoter n from a 2026-08-08 table | Traced to credible_sets.csv (5 July) | Corrected source |

### 3.3 Unverified findings

92 findings remain unverified: 62 medium and 30 low. By status, 73 are new, 16 are known and still open, and 3 are known and resolved. One finder reported each one with commands, but no second agent re-derived it. Do not act on them before verification. The 43 raw findings not counted here were merged into others as duplicates (`RUN/review.json`, `merged`).

The 20 high-severity findings that this section listed before the second pass are now in Section 2B.

A further number from the portfolio critic was not verified separately. `fig5b_protein_triage_summary.tsv` gives n_full_measured 5,836, rho_full 0.34582 and direction_full 0.61361. Fig 5B prints rho 0.46 and 66% agreement on a 3,054-gene subset chosen from the superseded July COLOC universe. LEG:388 discloses the subset choice, but the full-set values are not printed.

### 3.4 Full list of findings not verified in the main pass

Status key: new; known, open; known, resolved; in flight. For the 20 high-severity rows, the status cell also gives the second-pass verdict (Section 2B).

| ID | Severity | Status | Title |
|---|---|---|---|
| R08-04 | high | new; second pass: partly confirmed | Panel 7E treatment response: ratio should be 0.430, not 0.592, and the verdict should be untestable, not tested negative |
| R08-05 | high | new; second pass: confirmed | Panel 7F plots Geneformer, which has no clean study, as an above-zero gain and the best two-layer encoder |
| R09-01 | high | known, open; second pass: confirmed | Fig 7B legend, Fig 7 README and PAPER idea 8 say sequence models beat the task-native baseline; the only baseline run was an allele-identity ridge |
| R11-01 | high | new; second pass: confirmed | NAS-component separability is recorded as tested_negative, but its decision rule could not return GO, and its statistic cannot tell three axes from one |
| R11-02 | high | new; second pass: confirmed | Fig 5D 'partial' state counts pairs with essentially zero posterior mass in peaks, and the legend's 'localizes regulatory DNA' has no background |
| R12-01 | high | new; second pass: confirmed | Figure 4E calls exist only under untruncated Knapp-Hartung SEs; with the standard ad hoc truncation no program passes |
| R13-02 | high | new; second pass: confirmed | Title clause 1 is partly supported: one descriptive coding-share contrast with a wide interval, while most noncoding probability lacks regulatory context |
| R13-03 | high | new; second pass: confirmed | Fig 4E: all 16 fibrosis and 1 NAS program calls appear only after 93 identity-unresolved GSE174478 records are added; the 774-unit arm gives 0/0, and the arm was chosen after both had been computed |
| R13-04 | high | new; second pass: partly confirmed | Abstract and intro say the 117-program registry 'resolved remodeling into multicellular programs'; the program tests do not support this |
| R13-05 | high | new; second pass: partly confirmed | Continuum 'coordinated within-stage variation' claims omit that results depend on composition adjustment; Stromal ECM is unsupported in GSE162694 after adjustment and 0/43 systems pass the connected null |
| R14-02 | high | new; second pass: confirmed | Fig 2D 'direct MASLD' coding enrichment comes from the three PDFF GWAS; NAFLD/NASH alone has an interval that includes 1 |
| R14-03 | high | new; second pass: partly confirmed | '423 of 436 enzyme-driven genes' reflects 22-fold more tested enzyme signal pairs, not lower colocalization yield for direct traits |
| R14-04 | high | new; second pass: partly confirmed | Novelty statement ('what remains missing') contradicted by uncited prior work linking MASLD GWAS to liver expression |
| R15-01 | high | new; second pass: confirmed | Gene Catalog genetic layer omits the paper's own fine-mapped coding evidence; PNPLA3, TM6SF2, MTARC1 read as 'unresolved' and Fig 6C prints MTARC1 genetics 'indeterminate' |
| R15-02 | high | known, open; second pass: confirmed | Catalog PP.H4 column mixes ABF and SuSiE values with no method field; GNMT shows 0.706 (ABF) instead of its SuSiE 0.128 |
| G1-01 | high | new; second pass: partly confirmed | The only primary program call (Ductular injury, Fig 4C and abstract) is lost under the project's own 50-cell gate; it rests on two GSE189600 donors with 43 and 99 hepatocyte nuclei |
| G2-01 | high | new; second pass: confirmed | Fig 3E F1-vs-F0 count (1,338) comes from one cohort: 14 calls without GSE130970 |
| G2-02 | high | new; second pass: confirmed | F1-vs-F0 count falls to 104 to 485 under composition covariates; library size alone has no effect |
| G2-03 | high | new; second pass: confirmed | Fig 3F: 29 of 54 genes enter only through F0-to-F1 marks, which vanish without GSE130970; joint refit keeps 22 of 54 genes |
| G3-01 | high | in flight; second pass: partly confirmed | Fig 5E and Fig 5F print different Stromal ECM Moran excess because they are two different analyses; 5F is the right one for the legend claim |
| R01-10 | medium | new | Fig 2 text and the noncoding-context companion print unresolved and ABC fractions under different weightings |
| R02-05 | medium | new | Fig 2B colours HSD17B13 and MBOAT7 'Coding' by a hard-coded list, not by the plotted lead variant |
| R02-06 | medium | new | Fig 2E 'non-European-only' genes are often at loci where EUR GWAS already colocalize a neighbouring gene |
| R02-07 | medium | new | Fig 2C 'coding' class counts genes whose lead variant is protein-altering in a different gene |
| R03-04 | medium | new | F4-vs-F0 effects differ strongly by cohort; 57% of the 7,043 calls show significant cohort-by-stage interaction |
| R03-05 | medium | new | 12 of the 54 Fig 3F rows enter only through F0-to-F1 calls that depend on the 14 source controls in F0 |
| R03-06 | medium | new | Legend 3E says F0-reference sensitivities are shown in Figure S3; no S3 panel shows them |
| R03-07 | medium | new | GSE213621's binned stages give an unused held-out check of 3E direction |
| R03-09 | medium | new | Pooled 1,854 calls are mostly robust to library complexity, but 16% are lost under joint adjustment (sensitivity not reported) |
| R04-03 | medium | new | '62 donors across five datasets' overstates the donors and datasets that inform the stage slope |
| R04-04 | medium | new | Stromal ECM ambient attenuation of 64% is called unsupported only because it is tested in the 117-program family |
| R04-05 | medium | new | Fig 4C confidence intervals are unadjusted homoskedastic intervals: 11 of 24 shown rows exclude zero while 1 passes BH |
| R04-06 | medium | new | S4 native replication and omission panel still shows the 64-donor design, including a GSE174748 slope built on cirrhotic explants |
| R04-07 | medium | new | Independent-donor replication in GSE256398 (V3, null after Holm) is not reported beside the Fig 4 program claims |
| R04-08 | medium | new | Draft text describes Figure 4D as donor-paired localization contrasts; the installed 4D is the 80/37 discovery-versus-content matrix |
| R05-04 | medium | new | Installed S4 focal continuum forest prints the superseded 844-record roster and contradicts the legend on GSE162694 Stromal ECM |
| R05-05 | medium | new | S3-S4 legend breaks the composition reporting rule: single-arm adjusted values, and a GSE162694 statement that holds for one fraction set only |
| R05-06 | medium | new | 14 of the 18 strict maxT systems do not hold across composition adjustments |
| R05-07 | medium | known, open | Second co-primary continuum axis is counted as extra support but adds no information at system level |
| R06-06 | medium | new | 5B 'concordant' needs no RNA evidence; rho is driven by RNA-selected genes |
| R06-08 | medium | new | 5D counts gene-study pairs, so genes with many GWAS dominate |
| R06-09 | medium | new | MAT1A is 'unresolved' in the Catalog but routed as an established-state gene |
| R07-04 | medium | known, open | RNA 'indeterminate' covers 11,301 genes whose disease direction is resolved at BH 5%, including THRB and MTARC1 in Fig 6C |
| R07-05 | medium | new | Catalog omits 114 bulk-tested and 61 COLOC-tested Ensembl genes, so it holds 1,853 of the 1,854 Fig 3 canonical calls |
| R07-06 | medium | new | 'Concordant' class has no direction check and is undefined in METHODS; the paper's own direction test was null |
| R07-07 | medium | new | Readers cannot tell 'untested', 'unresolved' and 'not applicable' apart in the Catalog tables |
| R08-06 | medium | new | Panel 7E signed eQTL: verdict should be untestable, because even the sequence-attributable bound is below the detectable effect |
| R08-07 | medium | new | 'Matches rather than beats' in-cohort ridge asserts an equivalence the intervals do not support |
| R08-08 | medium | new | TranscriptFormer clean-only sign depends on one 2-donor study and reverses under donor weighting |
| R08-09 | medium | new | Legend 7F range '-0.065 to -0.002' describes GSE189600, not 'every other study', and study-level deltas are not on the panel |
| R08-10 | medium | new | GSE268273 was used to choose the shipped transform; the model card says the opposite |
| R10-01 | medium | new | Hotspot vocabulary-null numbers in the Fig 4 text come from the superseded 469-record roster; on the installed 457-record Fig 4E arm the Hotspot count is below the matched-random mean |
| R10-02 | medium | new | Discussion says the public perturbation stress tests are visible in the supplement, but no supplementary legend or Methods text describes them |
| R11-03 | medium | new | Fig 5E open-promoter counts use credible sets that lack the three MVP EUR studies included in Fig 5D |
| R11-04 | medium | known, open | ATAC directional-replication q values treat peaks as independent; donor-label permutation gives p near 0.03 |
| R11-05 | medium | known, open | 'source_dependent' records significance in one cohort, mostly with the same sign in the other, not disagreement between sources |
| R11-06 | medium | new | The 'scale of assay limitation' paragraph counts unloaded data and low power as assay limits and leaves out the positive coverage |
| R11-07 | medium | new | GNMT/MAT1A/CYP2C19 protein result is a bulk MASLD-versus-control contrast but is described as severity-linked loss of hepatocyte function |
| R12-03 | medium | new | Figure 5F 'supported' for Stromal ECM is a gene-set test with the 4 GSE192741 donors held fixed; no donor-level test reaches 0.05 |
| R12-04 | medium | new | NAS-versus-NAS0 counts are disease-versus-control contrasts: 45 of 46 NAS0 units are source controls |
| R13-07 | medium | new | Banned wording remains: 'only 13' for the trait-driver split (three places) and 'coding plus splice mass' in the Fig 2D legend |
| R13-08 | medium | new | Chance-level overlap is still written as separate or complementary layers in the abstract opening and Figure 2 |
| R14-06 | medium | new | Fig 2D direct-class coding mass is concentrated in five clusters; the ascertainment-matched control is not reported in the paper |
| R14-07 | medium | new | Fig 3G tests unsigned PP.H4 against signed log2FC; the magnitude test gives the same null answer and should replace it |
| R14-09 | medium | known, open | GNMT 'worked example' rests on one single-signal ABF result in a stratum with no eligible SuSiE pairs |
| R14-10 | medium | new | Abstract reports a nominal post hoc P after the prespecified test failed; a simple sign test on the same pairs gives P=0.095 |
| R15-03 | medium | new | Catalog rows state no denominator or uncertainty; gene pages show no effect size |
| R15-05 | medium | new | Molecular-object companion contradicts the gene calls: bulk_rna 'indeterminate' for all 1,853 RNA-supported genes; no gene carries a protein call |
| R15-07 | medium | new | Readers cannot identify the GWAS and several other input datasets: Methods names none of the 35 GWAS sources and the Catalog source graph cites no study |
| R15-08 | medium | new | No public route to any result table; Code availability lists superseded 'historical' scripts instead of the producers of the current numbers |
| R15-09 | medium | known, open | Catalog program layer is the historical 14-record August context; the current 117-program stage results and HMSMA replication are absent |
| R16-02 | medium | known, open | K16: Figure 3F selects genes significant in any of four separate BH families; only 15 of 54 survive a split FDR budget |
| R16-04 | medium | known, open | K12: H17 source-corrected result is a null; the 'opposite direction' call rests on gene-level nulls and fails the donor-label null |
| R16-05 | medium | known, open | K01: GSE193066 still identity-unqualified; RESULTS still calls Ductular injury 'within-person supported' from 54 unqualified participants |
| G1-02 | medium | new | GSE189600 enters the stage fit from 2 of 6 replicate samples, one sequencing run each, and about 0.5% of its nuclei |
| G1-03 | medium | new | Donor scores in the 117-program family include donors with 1 to 49 cells of the tested lineage, mostly T cells |
| G1-04 | medium | new | Methods do not state the per-cell QC of the atlas used for program scoring (no MT filter, 200-count floor, no doublet step) |
| G2-04 | medium | new | GNMT F0-to-F1 decrease (q 0.00913) depends on GSE130970 and composition; F2-F4 decreases are robust |
| G2-05 | medium | new | F3-to-F4 adjacent count (281) ranges 115 to 643 depending on which complexity covariate is added |
| G3-03 | medium | new | Main-text spatial numbers come from the hepatocyte-only model while Fig 5F plots the all-16 model; '4/4 donors positive' fails under the 5F model |
| G3-04 | medium | new | No between-program test backs 'matched spatial nulls distinguished' Stromal ECM from Ductular injury; in Guilliams the difference is untestable at n=4 |
| G3-05 | medium | new | Spatial scores use 26/49 Stromal ECM genes (55% of weight) and 9/22 Ductular injury genes (41%) because the input is the cell2location-filtered gene matrix |
| R02-08 | low | new | Palindromic SNPs kept without an allele-frequency strand check in the UKBB strata that supply half of SuSiE-positive rows |
| R02-09 | low | known, open | Decimal locus parsing remains in supplementary and pipeline scripts; main Fig 2 numbers unaffected |
| R03-11 | low | new | Fig 3D endpoint model source table in figures/main prints historical 418/115/303 units and 23,370 genes |
| R04-09 | low | new | The 80/37 lineage-match split has no uncertainty; donor bootstrap gives 75 to 80 matches with 8 unstable calls |
| R04-10 | low | new | Intro calls Stromal ECM a 'disease-associated program'; support is unadjusted bulk fibrosis only |
| R05-08 | low | new | 4H panel key says 'Exceeds size-matched null' without saying it is the permissive uniform null |
| R05-09 | low | new | RESULTS pooled-window estimand agreement 0.997 is the 844 value; the 774 run gives 0.849 |
| R05-10 | low | new | 4E composition-adjusted count (3/0) was never tested against the earlier cell-type fraction set |
| R06-10 | low | known, open | V4 HMSMA 32/32 uses a weaker adjustment and cannot replicate the 5F statistic |
| R06-11 | low | new | 5F source table gives Vu Ductular injury the same state as supported Stromal ECM |
| R06-12 | low | known, resolved | 5C rows were chosen for RNA-protein concordance, so their agreement is by construction |
| R07-11 | low | known, open | Six Catalog genes carry a genetics call computed on a different Ensembl ID |
| R08-11 | low | new | Nonlinear-vs-ridge range '-0.046 to +0.001' includes two comparisons with unmatched features and is not on panel 7D |
| R08-12 | low | new | Panel 7D severity ceilings come from other cohorts, not the three transfer cohorts |
| R09-03 | low | new | Fig 7B MPRA task uses 2 of the 4 measured contexts (HepG2 only); legend does not say so |
| R09-04 | low | new | METHODS describes an AlphaGenome allele analysis whose results appear nowhere in the manuscript |
| R10-03 | low | known, open | Myojin palmitate class-test rows in the r8 Gene Catalog still use v1 evidence classes; refit with r8 classes leaves the nonconfirmatory verdict unchanged |
| R10-04 | low | new | METHODS describes a sex-concordance filter for sex-stratified analyses the manuscript does not report |
| R11-08 | low | new | STATUS calls the four stellate peaks 'chromatin findings'; their overlap is about three times chance at nominal P 0.04 |
| R11-09 | low | new | Protein measurability by Catalog class is not reported; 85% of established-state genes have no protein measurement |
| R12-06 | low | new | Genetic/disease-state overlap is described as 'approximately the independent-sampling rate', but the interval allows up to about half the expected overlap |
| R12-07 | low | new | Spatial matched-null P=0.0001 is the resolution floor of 9,999 draws, printed as an exact value |
| R13-11 | low | new | Bulk units called 'donor-level' and 'participants' in continuum text; Figure 4 summary claims 'regulatory activity' that is supplementary |
| R16-07 | low | known, open | K13: authority lines still stale on H17/N3, COLOC r2 adoption and the August release |
| R16-08 | low | known, resolved | K17: the timed-out COLOC task completed on re-execution; only the STATUS note is stale |
| R16-09 | low | known, open | K06: C8 anchor guard still in the combine script but no current number depends on it |
| R16-11 | low | known, resolved | K02, K10, K11: resolved by disclosure or reconciliation |
| R16-12 | low | known, open | K03, K04, K07, K09, K18, K19, K20, K21: open in code or docs, no current reader-visible number |
| G3-07 | low | new | The 5F overlay section is the strongest of five sections; RESULTS cites the 4-donor value as if it were the overlay's |
| G3-08 | low | new | Methods do not name the spatial adjustment model, and the stated BH family (2 programs) does not match 5E's 21-program gate |

## 4. New analyses

### 4.1 Ranked plans

The rank is the run order I recommend:

1. Correctives first, because they change printed numbers whatever their result, and two of them are cheap.
2. Then the two tests that decide title clause 2.
3. Then the only direct test of the stated novelty link.
4. Then replication and Catalog additions.

Judge means (three judges) are given in 4.2. They do not set this order. None of these analyses has been run. Each audited prespecification is in `RUN/plans/PRESPEC_<ID>.md`. The audit working copies are in `RUN/scratch/audit-<ID>/`. CPU h is expected use, not the reservation: every plan requests the 90 h default wall time.

| Rank | ID | Title | Kind | Unit and n | Primary test | Destination | Claim it changes | CPU h | Audit |
|---|---|---|---|---|---|---|---|---|---|
| 1 | H57 | Fig 1B overlap within genes COLOC evaluated, against an expression-matched expectation | Corrective | Gene; 3,769 evaluated genes (426 colocalized, 220 DEG); DEG calls from 774 units | Exact stratified permutation over 10 AveExpr bins; O/E and MH OR with CI; family of 1 | Fig 1B, LEG:21, abstract, intro, Discussion; Fig 3G "Neither" becomes 3,149 | 14,704 universe, 33.46 expected, OR 0.757, DEG-only 1,129 | 0.5 | Passed after patches; added a one-gene-per-500-kb check and a control-only stratifier |
| 2 | H59 | 117-program stage test with dataset-specific variance on 60 correctable donors | Corrective | Donor; 60-donor universe; fitted n 25 to 55 per program (GSE202379, GSE244832) | GLS with varIdent(dataset); within-dataset permutation, 19,999 draws with a REML refit each; BH over 117 | Fig 4C, S4, ABS:48-51, LEG:234-237, PAPER:266-269 and 336-337, fig5 text | "Ductular injury alone met the primary 117-program rule"; 5F "two robust programs" | 8 | Passed after patches; stage now from metadata (score files label 2 cirrhotic explants as Steatosis) |
| 3 | H58 | Fig 2D coding contrast within ancestry; position-repaired context; Fig 2A support | Corrective | 1-Mb cluster; 207 direct and 6,208 enzyme instances; only AMR (10 direct, 63 enzyme clusters) is called after the audit | Ancestry-standardized difference, 10,000-draw cluster bootstrap; Bonferroni over 3 (coding) and 5 (context) | Fig 2D panel and legend, abstract, intro, PAPER, STATUS, RESULTS; 2A legend | Ratio 6.03 and "every difference interval positive"; 36.49% vs 77.85% unresolved | about 4 (360 reserved) | Passed after patches; a new coverage gate fails EUR (ratio 100), so C1 and C2 become no call |
| 4 | H52 | Program-pair covariance replicated in 258 GSE213621 disease units | Positive discovery (title clause 2) | Identity cluster; 371 discovery, 258 replication disease units | S1 sign concordance, S2a/S2b modularity, S3 community density; Holm over 4, plus unit-bootstrap conditions | New Fig 4 panel; supplementary table | Title clause 2; RESULTS:2313-2345 | 6 | Passed after patches; sign-matched baseline added |
| 5 | H54 | Frozen coordinated pairs in GSE281797 and GSE225740, and in F0-F1 | Positive discovery (title clause 2) | Participant; 91 and 87 (secondary 66 and 53) | S1 at k=1 and k=3 against a lineage-label permutation (1,000 draws); participant bootstrap; Holm over 2 | Supplementary panel beside Fig 4; RESULTS row | Title clause 2; K01 | 8 | Passed after patches; label-null calibration (NC3) added |
| 6 | H60 | Signed agreement between ALT TWAS z and disease t | Positive discovery (novelty link) | Gene, 14,703 under a circular-shift null; DE from 774 units | Decile-stratified Spearman; 10,000 within-chromosome shifts; block-jackknife CI; composition qualifier; family of 1 | Fig 1 supplement; main panel only if a direction reading passes | First signed all-gene test of PAPER:104-107 and 372-375 | 6 | Passed after patches; direction now also needs a jackknife CI that excludes 0 |
| 7 | H51 | Grouped stage effects replicated in 258 GSE213621 disease units | Positive discovery (replication) | Identity cluster; 371 discovery, 258 replication | OLS slope of GSE213621 log2FC on discovery log2FC; stratified permutation, 1,000 draws; Holm over 2; ratio to held-out-cohort slope | Supplementary S3 panel; 3E legend sentence; candidate Catalog field | Fig 3E and 3F gain, or are shown to lack, replication outside three cohorts | 116 (up to about 400) | Passed after patches; exclusion statement now uses the bootstrap bound |
| 8 | H53 | caQTL sharing at colocalized liver-trait signals | Reader-facing (Catalog), one test | Locus; 238 confirmatory liver loci; T1 compares 185 liver-only and 101 platelet-only loci | coloc.bf_bf PP.H4 > 0.8 and LD r2 ≥ 0.5; T1 logistic arm effect | Catalog annotation columns; supplementary table and panel | Abstract "regulatory DNA elements" | 40 | Passed after patches; arms made symmetric; caQTL overlap handling fixed |
| 9 | H56 | Per-gene composition sensitivity of the 1,854 calls | Reader-facing (Catalog) | Participant (774 units); 1,854 calls; 8 coding-by-input cells | Conditional test on D = b_a - b_u; intersection-union T_keep and T_lose; BH over 1,854 each | Catalog sidecar column; supplementary table and panel | Adds a per-gene statement; changes no current number | 48 | Passed after patches; selection-bias fix |
| 10 | H55 | Histologic fat in GTEx postmortem livers and the disease-state genes | Positive discovery plus Catalog column | Donor; 222 GTEx livers; 1,002 genes | Spearman of bulk logFC against the GTEx fat slope; Freedman-Lane permutation, 10,000 draws in 6 strata; family of 1 | Catalog annotation column; supplementary figure | None current; adds one statement | 16 | Passed after patches; exposure check added for 104 unmeasured slides |

Total expected use: about 252 CPU hours (116 + 6 + 40 + 8 + 16 + 48 + 0.5 + 4 + 8 + 6). An earlier total of 608.5 counted H58's 4 CPU x 90 h reservation as use. The H58 prespecification expects under 1 hour.

### 4.2 Why each plan matters

**Correctives (H57, H58, H59).** These change printed numbers in the abstract and on main panels. The paper needs them whatever new results are added.

- **H57** (corrective, 0.5 CPU h). Fig 1B is the paper's first quantitative statement of how the two maps relate, and it is in the abstract. As printed, 74% of its universe was never evaluated, and its expectation ignores expression. H57 fixes both in one test. The 3,769-gene universe drops the 7 untestable genes that the 3,776 in 2.3 kept. Direction may be stated only if the control-only stratifier (D3) and the one-gene-per-500-kb version (D10, 293 clusters) agree. Otherwise the paper states the evaluated-universe counts and "close to the expression-matched expectation". The 935 never-evaluated DEGs are reported as unobserved.
- **H59** (corrective, 8 CPU h). Fig 4C holds the paper's only family-level single-cell disease association, and the abstract's "resolved remodeling into multicellular programs" rests on it. H59 lets each dataset have its own residual variance, and refits that variance on every permutation draw. A pass gives a call that survives both checks found in 2.5 and 2.6. A null means the abstract, the Fig 4C legend and the 5F premise are reworded. The audit fixed a stage-source error that would have rebuilt the superseded 64-donor universe.
- **H58** (corrective, about 4 CPU h expected; 360 reserved). This bears on title clause 1. After the audit the EUR stratum fails the coverage gate, so the coding contrast will be called only in AMR. The AMR point estimate is already known: 0.0219 [-0.0486, 0.1083]. Most of H58's value is the position-repaired context family, the Fig 2A same-study support count and formal intervals for the replacement sentence. Recommendation: reword the 2D legend and abstract now from the review numbers, and run H58 to supply intervals. If CPU is limited, run family X and the Fig 2A count first.

**Tests that decide title clause 2 (H52, H54).**

- **H52** (positive discovery, 6 CPU h, judge mean 7.33). It tests whether the discovery program-pair covariance recurs in the 258 GSE213621 disease units. Those units were not used to define programs, edges or communities. GSE213621 was used in the gene-level cross-assay co-expression graph. A pass gives an identity-qualified replication that can replace GSE193066, worded as bulk program-score covariance beyond fibrosis bin and sex. A null removes the clause's last support, and under the 2026-09-08 rule the clause leaves the title (user decision). The audit found that the 371-unit rediscovery does not recover community 3 as a separate community, so a passing S3 shows density only.
- **H54** (positive discovery, 8 CPU h, judge mean 6.67). The same clause, tested in two cohorts never used for programs, with participants as units. Its F0-F1 sets rescue H7, which stopped at 14 single-cell donors. GSE225740 participant identity cannot be checked (no reads, no deposited sex, age or batch).
- H52 and H54 need one joint decision rule, written into both prespecifications before either runs (see 6.1).

**Test of the stated novelty (H60).** 6 CPU h; not scored by the judges (judge mean recorded as 0; it came from the portfolio critic's candidate, `RUN/scratch/critic-cand-1`). It is the only plan that tests the paper's central link with a signed statistic over all genes. A grep for OTTERS or TWAS returns 0 hits in all five authority documents. The bound reading allows a definite statement: genome-wide signed agreement within |rho| 0.03. A stated direction needs the composition qualifier and a jackknife CI that excludes 0.

**Replication (H51).** 116 CPU h, judge mean 7.67 (highest). Fig 3E and 3F rest on three cohorts, with GSE135251 supplying 53.5% of units, while GSE213621's 258 binned disease units go unused for stage. H51 gives the first stage replication outside those cohorts. It is not blind, because GSE213621 stage labels already feed the continuum held-out check and the 11,130 stage-supported genes. It also flagged that the GSE213621 rows of `transcript_stage_trend_cohort.tsv.gz` use 326 units. If the controls carry stage 0 there [Inference], that slope mixes disease against control into stage.

**Reader-facing Catalog additions (H53, H56, H55).**

- **H53** (40 CPU h, judge mean 7.0). It either gives the Catalog a counted regulatory-region field, so the abstract phrase can stay with a number, or shows that the Currin caQTL source cannot, so the phrase is dropped. The liver-versus-platelet test adds at most one supplementary sentence.
- **H56** (48 CPU h, judge mean 5.0). No per-gene composition statement exists today. H56 adds three counts and a sidecar column, and changes no number or class. The audit corrected a selection bias: near the call threshold, the original test rejected 9.7% to 13.7% of the time when the true ratio was exactly 0.5.
- **H55** (16 CPU h, judge mean 5.33). A disease-state test outside clinical cohorts, plus an annotation column. GTEx records no BMI, alcohol use or cause of steatosis, so the result is about histologic fat, not MASLD. 104 of the 222 slides still need morphometry.

## 5. Candidates that did not become plans

**Pipeline counts.** 57 candidates were generated:

- 50 from lenses L01 to L09.
- 4 from the referee lens L10, written after review findings were confirmed.
- 3 from the portfolio critic (C-01 to C-03).

The code check dropped none. Triage of the 50 lens candidates kept 35, merged 14 into others and dropped 1. The dropped one is L02-05: it asks the same question as an archive candidate and has the H40 failure mode, where the estimand is nearly collinear with the nuisance. The gate cap admitted 14 of the 35, so 21 were never gated (table below). Of the 18 gated lens candidates (14 from L01 to L09 and 4 from L10), 7 passed: L01-01, L02-01, L03-04, L07-01, L10-01, L10-03 and L10-04. Six of the seven became plans H51 to H56. Of the critic's three, C-01 passed and became H60.

**Passed the gate, killed by the judges:** L01-01 (GGT-only vs ALT/AST-only colocalized genes in biliary vs hepatocyte compartments). All three judges scored it 2 and killed it. The MDE is 0.43 SD, but the whole 289-locus COLOC set departs from background by only 0.15 SD. GGT1 itself has an index of 0.07, and the class labels partly record which biobank measured which enzyme. It was the only candidate that tested the inherited-to-remodeling link directly.

**Failed the gate or validity check:**

| ID | Title | Reason |
|---|---|---|
| L02-02 | Promoter H3K27ac co-change in 99 GSE267145 participants | The COL1A1 positive control cannot pass (+0.020 log2, p 0.80, against an MDE of 0.241). The permutation null SD of 0.27 gives an MDE near 0.75. NOR and NAFL participants are all women. |
| L06-02 | Within-MVP cross-ancestry transfer of enzyme lead effects | MVP effects are in ancestry-specific SD units, so the ratio cannot separate per-allele effect from scale. Equivalence is showable at only 7/8/1 ALT loci. |
| L04-01 | UKBB vs MVP replication of COLOC calls | No MVP rows exist for the positive control on r2. Design effect 2.23 and 1.96. The bin outcomes had already been seen. |
| L01-04 | Same-nucleus enhancer-gene links at colocalized loci | 19 testable loci, not 51 (MDE OR 4 to 5). The SORT1 control falls below the 5% detection filter. COLOC genes are 1.61x more expressed. |
| L03-01 | Genetic observability fields for every Catalog gene | Wrong outcome file (pre-PolyFun fits; 407 genes change label). The control pool is selected on the outcome. |
| L09-01 | Equivalence-tested negatives for the bulk map | Two of three conditions remove nothing. 391 of the proposed negatives have F4-vs-F0 calls in Fig 3E. |
| L03-03 | Within-lineage donor stage slopes from snRNA-seq | The COL1A1 control fails in fibroblasts (-0.012, p 0.91). MDE 0.67 to 0.90. Ambient load rises with stage. |
| L03-05 | Held-out fibrosis replication column for Catalog calls | 397 vs 127 calls; MDE 0.21 to 0.24 against a plausible 0.04 to 0.07. Confounded by fibrosis relatedness. The gene-independent Wald test is invalid (63 to 89x variance). |
| L01-05 | GGT-class genes and the ductular compartment | Expected effect far below the MDE. Class labels are confounded by MVP lacking GGT. Depends on the killed L01-01. |
| L08-03 | GTEx v10 sQTL, apaQTL and small-RNA QTL classes | The primary estimand is 0 in the Resource credible sets. The null cannot fail. Low sensitivity (13 of 436). |
| L10-02 | Genome-wide EUR fine-mapping of enzyme and MVP NAFLD GWAS | Trait class is collinear with GWAS power. The pilot showed proxy choice and LD panel drive the coding fraction (PNPLA3 PIP moves to an intronic variant). lambda_s is NA for all 14 EUR studies. |

**Critic candidates.** C-01 (TWAS-predicted expression against disease direction) passed and became H60. C-02 (are bulk-coordinated program pairs co-patterned within HMSMA livers?) failed because its test is invalid for its claim: under random relabeling of the 662 pairs, the donor sign-flip test rejects at p<0.05 in 40.8% of 4,000 labelings. C-03 (Catalog per-gene columns for risk-allele direction, stage effects and liver protein) failed its concordance test on power: the MDE is 0.65 at 83 to 89 non-MHC loci, against a disclosed agreement of 0.60. Its column-building part needs no test, and it overlaps gap 6.1(1).

**Never gated (21).** These passed triage but were beyond the gate cap of 14. Nothing here says whether they are feasible. A follow-up gating round would cost about 21 gate agents plus validity and prior-art checks for the survivors.

| ID | Title |
|---|---|
| L05-04 | Within-person reversal of Catalog disease-state genes after bariatric surgery (GSE83452) |
| L06-06 | Does the adult established-disease signature appear in adolescent and pediatric MASLD? |
| L08-04 | Healthy-liver zonation columns for every Catalog gene from live-donor spatial mRNA and single-hepatocyte protein |
| L07-03 | Liver H3K27ac allelic imbalance at fine-mapped credible-set variants, pooled across 99 participants (rescue of H25) |
| L05-03 | Plasma detectability and direction of Catalog proteins in PXD051911 |
| L04-03 | Follow-up assay observability of each Catalog class: liver proteome, plasma proteome, CosMx panel and single-cell lineage |
| L02-06 | Do disease-acquired liver H3K27ac regions carry extra liver-trait heritability? S-LDSC on GSE267119 gained and lost regions |
| L01-03 | Does GGT heritability concentrate in cholangiocyte chromatin more than ALT heritability does? |
| L08-02 | Does MASLD flatten the hepatocyte zonation gradient? Healthy human protein zonation applied to donor-level hepatocyte RNA |
| L09-03 | Re-scoring published MASLD liver gene claims with the evidence-state vocabulary |
| L04-06 | Competitor gene lists (Kumak 2026, Elison 2025) requalified by Catalog evidence state |
| L06-03 | Sex heterogeneity of genetic effects at independent liver-enzyme loci, including chrX |
| L05-01 | Genetic direction label for MASLD trial targets from liver eQTL COLOC and OTTERS TWAS |
| L04-05 | Nearest-gene contradiction at liver-trait loci and 3,335 genetic bystander negatives |
| L06-05 | Menopause-age shift of the female liver transcriptome toward the male pattern in 910 severely obese livers |
| L06-04 | Stromal ECM program in women versus men at matched fibrosis stage |
| L07-02 | X-inactivation skew in biopsy RNA as a readout of clonal hepatocyte expansion across fibrosis stage (rescue of N2a and H30) |
| L03-06 | Constraint, secretion and essentiality fields for the 1,129 genes routed to stress perturbation |
| L06-01 | Genetic-ancestry imbalance between controls and disease in the 774-unit bulk fit |
| L09-06 | Recall of the SuSiE-COLOC pipeline at real loci, for reading genetic non-calls |
| L09-02 | Spike-in recall and false-call rates for the bulk evidence states |


## 6. What the portfolio still lacks

### 6.1 Gaps named by the critic and the judges

1. **The Catalog does not carry per-gene results the paper already has.** Its evidence rows are 56,466 = 28,233 bulk canonical + 28,233 COLOC (critic, reading `passport_evidence_long.parquet`). Fig 3E stage contrasts, Fig 5B liver protein (3,054 measured genes) and donor-level lineage results are absent. Adding them would add more for readers than any new hypothesis test here. No plan does it.
2. **No participant-level test links inherited variation to remodeling.** H3 bounds per-allele program effects at about 0.1 to 0.2 SD in 969 persons. The total effect expected through stage is about 0.07 SD [Estimate], below what 969 persons can detect. H60 is the better-powered gene-level route. Even if every plan passes, no analysis joins the title's two halves at the participant level.
3. **No test of coordination in tissue space.** H52 and H54 test bulk covariance only.
4. **No replication of the Fig 4C program-stage result in new donors at adequate n.** V3 had 13 to 17 donors. The critic does not recommend GSE296875: 22 donors have no significant fibrosis, 3 are missing, and 6 of the 8 "yes" donors sit in wells 7 and 8. Leave-one-well-out leaves almost no fibrosis contrast. GSE289173 is a held-out cohort with disease and control labels only, and opening it is a user decision.
5. **H52 and H54 need one joint decision rule.** Recommendation: before either job starts, write the same rule into both prespecifications. The clause stays as bulk program-score covariance only if H52 replicates and at least one primary H54 cohort passes. Otherwise it becomes "maps cross-sectional multicellular remodeling".
6. **A clean trait-class coding contrast may not be estimable** with the GWAS in hand. The H58 coverage gate and the L10-02 gate failure both point this way. The paper should be prepared to state the 2D result as ancestry-stratified and descriptive.
7. **Reader access** (unverified R15-07 and R15-08). Methods name none of the 35 GWAS sources, and no public route to result tables exists. No plan addresses this.
8. **Competitors.** R14 found that competitors offer things this Resource does not:
   - Kamzolas 2026: a validated 57-gene plasma panel.
   - LiRNA: an open web app over more than 3,000 biopsies.
   - Gobeil 2026: in-sample liver eQTL in 504 genotyped people.
   - Saliba-Gustafsson 2024: CRISPRi validation.
   - Elison 2025: chromatin-level fine-mapping to cell types.

   R14-04 (unverified) says the novelty statement does not cite prior work that links MASLD GWAS to liver expression.

### 6.2 Data applications to file in parallel

These are controlled-access sources from the brief. None may be a completion requirement. File them in parallel only if the author wants the gap closed in a later version. The mapping to gaps is my [Inference].

| Source | Gap it would address |
|---|---|
| Gaulton 86-liver multiome | In-sample chromatin for a regulatory DNA layer (2.8, H53); new-donor single-cell stage replication (gap 4) |
| SteatoSITE | Staged bulk with histology, for Fig 3E replication beyond three cohorts and for covariance replication (gaps 3 and 5) |
| NASH CRN DB2 / NIMBLE | Staged participants with histology (same use as SteatoSITE) |
| UKB-PPP and imaging; All of Us CT | Participant-level genotype with liver phenotype, for the participant-level link (gap 2) |
| CARDIA/LITMUS SomaScan; Quebec pQTL | A protein layer for Catalog genes (gap 1) |
| JGA | Listed in the brief as controlled; not mapped to a gap here |

Public leads that need no application (confirm accessions before use): SCP3359 Xenium (18 livers), GSE210077 MERFISH plus snRNA and Zenodo 14795740 (Yakubovsky live-donor spatial) for coordination in tissue space; FinnGen R13 SuSiE credible sets for EUR direct-trait fine-mapping.

## 7. Method of this round

**Finders.** Sixteen review dimensions (R01 to R16) covered Figs 1 to 7, statistics, wording, literature, the Catalog and the K01 to K21 items. Three gap finders covered the single-cell Fig 4C design (G1), the bulk stage Figs 3E and 3F (G2) and spatial Figs 5E and 5F (G3). Together they raised 175 raw findings. Outputs are in `RUN/scratch/review-R01/` to `review-R16/` and `review-G1/` to `review-G3/`.

**Verification design.**

- 20 findings went to verification. Each was re-derived from producing files by a reproduce agent (`RUN/scratch/verify-<ID>-reproduce/`). 12 also had a separate coverage agent (`verify-<ID>-coverage/`). That is 32 verification runs.
- All 20 were confirmed. Section 3.2 lists the sub-claims the verifiers corrected or disagreed on.
- Three decoys were mixed in, each checked by both agent types: X-01 (true) was confirmed, and X-02 and X-03 (false) were refuted. That is 3 of 3 correct.
- Second pass: the 20 remaining high-severity findings each got one blind reproduce verifier. Each verifier read only its own claim (title, claim, locations) from `RUN/context/verify2_claims_blind.json`, not the finder's evidence or severity. Two false decoys (G4-01, G4-02) were given ordinary-looking IDs, and both were refuted. Over both passes, decoys were 5 of 5 correct. With five decoys this shows the verifiers can refute a false claim. It does not estimate an error rate, and single-verifier verdicts were not downgraded.
- Main-loop re-derivation of R01-01, R03-03 and R05-02 (Hotspot part) with new code: all three match (Section 2C).

**Hypothesis pipeline.**

- 10 generation lenses (L01 to L10), one triage pass, 18 feasibility gates on L candidates and 3 on critic candidates.
- 7 validity checks and 7 prior-work checks.
- Three judges: hostile referee, P(informative) and resource utility.
- One portfolio critic, 10 planners and 10 plan audits.
- Every audit patched its plan and passed it. The patched prespecifications record each patch and the commands behind it.

**Jobs.** sbatch jobs were capped at 8 minutes on `cpu` with at most 4 CPUs and 16 GB. Each had a 4-digit numeric name, and all are listed in `RUN/job_index.tsv` (42 jobs: 32 main pass, 10 second pass). `sacct` on 2026-10-15 gave 36 COMPLETED, 3 TIMEOUT and 3 FAILED, with none still running. The non-completed jobs:

- 1031: TIMEOUT after the pooled models. Its per-cohort GSE126848 part was rerun on the login node (n=30, 13 s).
- 3081: TIMEOUT. The L03-04 positive-control refit printed no fit. The gate passed on the Stromal ECM membership counts.
- 3831: TIMEOUT after the H51 grouped discovery refit. The permutation and SE steps did not run.
- 3051: FAILED after 59 s. The L01-04 detection and testability counts finished, but the pilot link model did not.
- 3101: FAILED after 3 s, when `fread` read the GSE225740 header as data. The same script, with header=TRUE, then ran on the login node (18.7 s, 745 MB).
- 3841: FAILED at the final summary print of the L10-04 caQTL pilot, after both output tables were written. The summary was recomputed on the login node.

Other notes from the main pass:

- Job 1161 completed but keyed its join on the wrong GWAS, and a login-node recomputation replaced it.
- Job 1032 skipped F1 leave-one-cohort-out, the F2/F3 interaction and the complexity covariates. Second-pass job 4191 supplied the F1 leave-one-cohort-out values (Section 2B, G2-01).
- Job 3111 (H55 morphometry pilot) finished 118 of 222 slides.

**Agents.** 134 in the main workflow (`wf_a8d2ba5f-3d1`) and 22 in the second pass (`wf_00cf4a6c-d4e`). Two more started in a second-pass launch that was stopped within a minute because a placeholder was passed instead of the claims; they submitted no job and wrote no file. Total: 158, under the planned maximum of 160.

**Edits outside RUN.** One workflow agent imported a repo script in place inside a figure-candidate snapshot, which left two `.pyc` files and an empty `__pycache__` directory under `figures/candidates/fig2-coordinate-rerun-20261013T213349Z_22218816/snapshot/scripts/figures/`. I checked that they were created after the run started, then removed them; the second pass prompts forbade in-place imports. After both passes, the sha256 manifest of `docs/`, `figures/main/` and `docs/manuscript/` (3,441 files) matched the start of the run, and `git status` was identical to the start.

`integrated_atlas.h5ad` was not loaded. The 104 integrated donors and 1.23M cells were checked only against RESULTS and the 4B source table.

**Seeds.** R12 opened every permutation or bootstrap null behind a main claim and found none unseeded. One null prints its resolution floor as an exact value (spatial P = 0.0001; R12-07, unverified). This round's own nulls used seed 20261015 (review) or 20260923 (plans).

**Literature search (R14).** PubMed worked; one rate-limit error was retried after a pause. Confirmed by ID:

- Piras and DiStefano 2024 (PMID 38565287)
- Gobeil 2026, JCEM (41635218)
- Jamialahmadi 2024, Nat Med (39653778)
- Saliba-Gustafsson 2024, Hepatology (39190705)
- LiRNA 2026 (42264302)
- Kamzolas 2026, Nat Metab (42448794)

Retrieved from searches: Li 2025 Nat Genet, a 2026 Nat Genet macrophage multi-omics study, Nat Immunol 2025, Cell 2025. Elison 2025 was retrieved from medRxiv.

Not retrieved: Zhu 2026 Nat Genet (no PMID returned) and Kumak 2026 bioRxiv. HLiCA, Du 2026 and Yao 2026 are known only from an archived repo plan. bioRxiv has no keyword search. PDFF GWAS sample overlap was not checked.

**Subsystems no earlier review touched.**

- AlphaGenome (R09): no manuscript number uses it, so the K21 items change no printed value. R09 instead reviewed the Fig 7B MPRA claim. Job 1091 reproduced the stored gains and added Holm and max-|t| inference.
- RNA-seq scripts 200 to 297 and `Analysis/Perturbation` (R10): none feeds the manuscript, apart from the Hotspot membership file behind the vocabulary null (R10-01).
- Chromatin, proteomics, GSE287826 and GeoMx (R11): CosMx can observe 12/117 programs; GeoMx and GSE287826 can observe 0/117. Same-participant H3K27ac and RNA coupling is already described as not independent validation.

**Limits.** These were not checked:

- Fig 4B, 4F and 4G values; Fig 5C marker fills.
- Supplementary S2 to S4 numbers beyond pdftotext.
- METHODS and intro text clause by clause.
- Fig 2G to 2I locus views; Fig 3D PCA coordinates; 3E NMF activity values.
- The Fig 5B protein model (not refit).
- Whether the COLOC r2 release changes the Fig 5E open-promoter counts.
- Genotype overlap of GSE268273 and GSE276114 with the training libraries.
- A per-cohort complexity adjustment for GSE135251.
- The continuum coefficients without the GSE213621 controls.

The synthesis agent wrote Sections 1 to 7 from the finder, verifier, gate, audit and critic outputs. The main loop then added Sections 2B and 2C, the decoy and agent accounting, the job table and the never-gated list. It also corrected the candidate counts in Section 5 and the CPU total in Sections 1 and 4. Values in 2B come from `RUN/verify2_result.json`. Values in 2C come from `RUN/scratch/main-rederive/`.
