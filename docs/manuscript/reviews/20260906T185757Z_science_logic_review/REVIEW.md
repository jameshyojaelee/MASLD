# Science and logic review of the MASLD Resource manuscript (2026-09-06)

Scope: docs/PAPER.md, STATUS.md, RESULTS.md, ROADMAP.md, README.md, HIGH_IMPACT_PATTERNS.md, literature/REVIEW.md, decisions/*, the drafts in docs/manuscript/draft/, figures/main/*/README.md and panel indexes, and the Figure 7 authority in Analysis/MASLD_Model_Benchmark. Numbers marked [verified] were recomputed from the producing file in this review. Findings marked [L1] come from the one review agent that completed before the session usage limit stopped the 12-lens workflow; each of those cites its files and lines. Everything else is my direct reading. Literature statements not checked online are marked [Unverified]. No repository file other than this directory was modified.

## Verdict in five sentences

The manuscript is rigorous about what it cannot claim and, as drafted, has no claim it can. Its organizing observation, that the genetic and disease-state maps are "largely non-overlapping", is the independence null restated: two sets of 437 and 1,066 genes in a 14,920-gene background are 92.9% disjoint by arithmetic and the paper observes 95.4%. The two worked examples are known biology carried by weak or absent genetics, and Figures 4F to 6 are descriptive or null after the paper's own stress tests. As drafted this is a Cell Genomics or Genome Biology Resource paper (2025 impact factor about 9). Three positive, replicated, quotable results already exist in the repository and none is in the abstract; put them first, split the model benchmark into its own paper, and the realistic ceiling rises to Journal of Hepatology (40.1) or Nature Metabolism (27.5).

## What holds

- Participants and donors are the inferential unit everywhere I looked. Cells, spots, arrays, runs never inflate n.
- The overlap arithmetic (20 observed, 31.2 expected, OR 0.616, p 0.037; fragment arm 18 / 26.4 / 0.661 / 0.10), the flat expression-matched ladder, and the threshold-free rho 0.014 all reproduce from `RNA-seq/results/genetics_state_interface/20260817T133951Z_promoted_coloc/overlap_by_arm.tsv` and `figures/main/fig3_bulk_transcriptomics/source_tables/fig3g_raw_treat_interface.tsv` [L1, verified].
- The rule "never call colocalized genes buffered or constrained" is correct and is obeyed in the draft.
- The coding-mass asymmetry (RESULTS.md:922-995) has a class-label permutation null, BH over its six-test family, a difference-of-fractions statistic suited to zero-inflated masses, and an effect-size-matched ascertainment control. Retiring the ABC-enhancer depletion was right.
- The benchmark retractions of 2026-09-01 (sigma_e is a range, ICC does not transport, transfer is "not worse than native", intervals are not calibrated) were done honestly and are recorded.
- The 2026-08-30 Figure 7 decision keeps the prohibitions on a universal target rank and on biomarker performance.

## Ranked findings

### Must fix before submission (severity 5 to 4)

**1. The headline is the null. Severity 5.** [verified; L1-1]
Where: 01_intro.md:57-64, abstract.md:31-32, fig5_discussion.md:77-84, PAPER.md idea 1, Figure 1B.
As written: "The organizing observation is that ... inherited genetic evidence and the established disease-state transcriptome provide complementary, largely non-overlapping maps."

| quantity | value |
|---|---:|
| jointly testable genes | 14,920 |
| colocalized genes in background | 437 |
| DEGs in background (adopted release) | 1,066 |
| expected overlap under independence | 31.2 |
| observed overlap | 20 |
| non-overlap expected by chance | 92.9% |
| non-overlap observed | 95.4% |
| fragment substrate, Fisher p | 0.10 |

Problem: any two sets this small are more than 90% disjoint. RESULTS.md already says so ("what independent sampling of the transcriptome produces for two sets of this size") and forbids upgrading the depletion. The intro then calls it the organizing observation and the Discussion calls the separation "biologically informative". A referee writes that the central observation is arithmetic.
Fix: one Results sentence and a supplementary panel: "the two maps intersect at the rate expected under independence". Delete "organizing observation" and "biologically informative". If a two-map framing stays, replace membership overlap with a test that can fail: do colocalized and disease-state genes differ in cell-type expression, protein response, or coding architecture.

**2. The Discussion and Figure 1B carry retired July numbers. Severity 5.** [verified; L1]
Where: fig5_discussion.md:102-106 and :127; `figures/main/fig1_atlas_overview/fig1b_complementary_maps.pdf`.
As written: "adjusted odds ratio of 2.1 against 9.1", "the corrected gsMap analysis", "the strict drug benchmark ... does not outperform modality count", "covers only 75 of 473 primary genes".
Problem: none of these is in RESULTS.md, STATUS.md, or PAPER.md. They trace only to the July ledger, where the gsMap class gate and the coding-coverage gate are marked FAIL and the gsMap and drug numbers are listed under "Retired from the active manuscript" (`RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v5-candidate-2026-08-10/inputs/BG001-DECISION/source_snapshot/docs/manuscript/NUMBERS.md:18-21`). Figure 1B prints 34 convergent of 447 SuSiE genes and a 1,261-gene "Canonical DEGs (interval-null FDR gate)" set; the adopted release gives 20 of 437 under the canonical threshold. The "473" is the superseded count.
Fix: delete the three Discussion paragraphs or rebuild them on the adopted release; rebuild Figure 1B or drop it.

**3. The abstract has no positive quantitative result, and three exist. Severity 4.** [L1-2; my reading agrees]
Where: abstract.md, PAPER.md "Claims we can make".
Problem: of the abstract's 13 sentences, none gives an effect size, a q value, or a replicated positive. The repository holds three:

| result | number | where |
|---|---|---|
| direct MASLD/liver-fat credible sets are coding-heavy | 22.4% vs 3.4% protein-altering PIP mass, ratio 6.55, BH q 6e-4, survives effect-size matching | RESULTS.md:929, 943-995 |
| a bulk-RNA severity instrument transfers | 3 external cohorts, 585 participants, 0.68-0.79 of each cohort's ceiling, equal to a within-cohort ridge | MODEL_RESULTS_SUMMARY.md:59-62 |
| multicellular coordination beyond histology | 2,426/6,195 program pairs, modularity 0.358 vs 0.138 / 0.089 nulls | HIGH_IMPACT_PATTERNS.md rows 2-3; candidate only, see 5 |

The first is disclaimed by the Discussion's own "75 of 473" sentence. The second is not mentioned anywhere outside Figure 7.
Fix: rewrite the abstract around one or two of these numbers. Add the coding-mass asymmetry to "Claims we can make".

**4. The F0 reference mixes source controls and disease by cohort. Severity 4.** [verified]
Where: adjacent fibrosis contrasts, RESULTS.md:490-499; fig3.md:46-53; the GNMT F0 to F1 effect at RESULTS.md:78.
Crosstab of `figures/candidates/nas-vs-nas0-extensions-2026-08-17/analysis/stage_extensions/stage_sample_manifest.tsv` joined to `five_cohort_sample_manifest.tsv`:

| cohort | F0 controls | F0 disease |
|---|---:|---:|
| GSE130970 | 23 | 0 |
| GSE135251 | 8 | 38 |
| GSE162694 | 0 | 35 |

Problem: F0 is 31 controls and 73 disease. In GSE130970 the F0 to F1 contrast is control versus disease; in GSE162694 it is F0 disease versus F1. The cohort term cannot separate that, because control status is nested in cohort. The 1,471 F0-to-F1 DEGs, the "expands early and then stabilizes" pattern (HIGH_IMPACT row 20), and the GNMT and MAT1A F0-to-F1 stage effects all partly measure disease versus control.
Fix: rerun the adjacent contrasts with disease-only F0 (73 participants) or with a control-status term, and report both. The NAS0 group has the same problem in reverse: 32 of 46 NAS0 participants are GSE162694 (see finding 11).

**5. The 117-program registry does not beat random gene sets on fibrosis, and the paper does not say so. Severity 4.** [L1-5; not recomputed]
Where: `RNA-seq/results/manuscript_release/candidates/resource-f-five-coloc-v6-candidate-2026-08-10/workstreams/BULK-PROGRAM-MAP-v9/HANDOFF.md:69-73`; HIGH_IMPACT_PATTERNS.md row 50.
Problem: against size- and expression-matched random gene sets, 26 of 114 Hotspot programs are fibrosis-supported versus 18.4 expected (p 0.134), whereas Hallmark (3.43x, p 0.010) and cNMF (3.04x, p 0.005) beat their nulls. The paper's coordinate system for Figures 4 and 5 is therefore a vocabulary that performs like random on the paper's own severity anchor. Add that the two featured programs were hepatocyte-labelled, turn out to be fibroblast and cholangiocyte by content, and the ECM one is largely ambient RNA (beta 0.401 to 0.130). A referee asks why Hotspot, and why the registry's two headline programs are mislabelled.
Fix: report the vocabulary null in the Figure 4 supplement. Reframe 4C to 4E as "prespecified coordinate system", not discovery. Consider replacing the two exemplars with programs whose lineage label matches their content.

**6. The title promises a model-learning result that the abstract, intro, Discussion, legends, STATUS.md and README.md do not contain. Severity 4.** [L1-3]
Where: PAPER.md:9-11 (three-clause title) versus title.md:10-11 (two clauses); STATUS.md:38 and :60 and README.md:63 ("six-figure"); 05_figure_legends.md:8 and fig5_discussion.md:70 ("Across Figures 1-6").
Problem: the seventh figure was added on 2026-08-30 and nothing downstream was updated. The positive half of the Figure 7 claim ("a small margin over a baseline that uses allele identity alone" on a synthetic reporter assay) has no number in PAPER.md, and MODEL_RESULTS_SUMMARY.md records signed-direction accuracy at chance for the DNA language-model heads (0.5016 vs 0.5033 for the allele ridge).
Fix: see the strategy section. Either the model result gets an abstract sentence and a Discussion paragraph, or the title returns to two clauses and Figure 7 leaves.

**7. The worked examples are known biology on weak genetics. Severity 4.** [my reading; literature Unverified]
Where: PAPER.md narrative spine; fig2.md:52-73; fig3.md:73-91; RESULTS.md:207-234.
Problem: GNMT has no primary SuSiE colocalization (PP.H4 0.14 for ALT); the only support is a single-signal ABF PP.H4 0.71 in one direct-NAFLD stratum, retained as a "sensitivity observation", and the variant-to-protein-to-reporter link was published by Ji et al. 2012, which the draft credits. Loss of GNMT and MAT1A expression in human MASLD, and steatohepatitis in Mat1a and Gnmt knockout mice, are established [Unverified, from memory: Mato and Lu laboratories, 2000s]. The second example, an ECM program scoring highest in fibroblasts and a ductular program in cholangiocytes, is a property of the gene lists. Neither example produces a result a reader did not expect.
Fix: drop the "two recurring examples" device, or choose examples with primary SuSiE support and a non-obvious outcome. The GNMT locus can stay as a one-paragraph Catalog vignette.

**8. The continuum layer is three main panels of descriptive output on a borrowed axis. Severity 4.** [my reading of RESULTS.md:516-802]
Where: Figure 4F to 4H, S4Q to S4W, PAPER.md's long Figure 4G/4H paragraph.
Problem: the axis is the Kamzolas 145-gene signature. The paper's own competitive nulls find no system continuum-specific under the connected null and one under the permissive null; 43 systems have an effective number of 15; the 2.79x continuum-versus-stage ratio falls to 1.45 to 1.57 when compared fairly; every window is descriptive. What survives is one clean negative: within-stage variation is a single diffuse axis, not a set of system-specific programs. That is a quotable result and a contrast with Kamzolas, but the draft buries it under 33-of-43 and 35-of-43 counts.
Fix: keep 4F (two focal programs, stage-adjusted, two cohorts) as one panel; state the diffuse-axis conclusion in one sentence with the competitive-null numbers; move 4G, 4H and the window heatmaps to the supplement. Move the Kleiner test-retest paragraphs (RESULTS.md:745-802) out of the continuum section; they belong with the severity instrument.

**9. Hedge density leaves no sentence quotable. Severity 4.** [L1-4]
Abstract: 7 hedges in 13 sentences. Discussion: 36 in 59, with three paragraphs above one hedge per sentence. The abstract also contains "withheld" and "await" placeholders. Fix: one claim per paragraph, stated once with its number; one Limitations paragraph; no placeholders in a submitted abstract.

**10. Figure 7 dilutes the paper and one panel rests on a two-donor study. Severity 4.** [previous session, verified then; this session read only]
Where: figures/main/fig7_model_framework/; panel 7F; release/encoder-benchmark-v1/MODEL_CARD.md:35-36.
Problem: the UCE-4L versus liver-encoder contrast in 7F is an unweighted mean over five studies in which GSE189600 (2 donors, 149 cells) contributes -0.455 against -0.001 to -0.005 elsewhere; without it the delta crosses zero. More broadly, a hepatology or genetics editor gains nothing from a benchmark of AlphaGenome, Borzoi, UCE and scVI heads in a Resource paper, while a methods editor would value the exposure audit and the "a linear model of bulk RNA is at the ceiling; a nonlinear model with no pretraining adds 0.000 ± 0.067" result on their own.
Fix: split Figure 7 into a separate methods paper (strategy section). Fix 7F by weighting studies by donors or dropping two-donor studies before any use.

### Fix (severity 3)

**11. NAS0 is 70% one cohort.** [verified] NAS0 = 4 / 10 / 32 across GSE130970 / GSE135251 / GSE162694; NAS1-2 = 13 / 31 / 21. The middle contrast (NAS1-2 vs 3-4) is cohort-balanced, so its three DEGs are probably a real null; the outer contrasts are not balanced and the first includes source controls. Report per-cohort effects for each contrast before calling the pattern "punctuated".

**12. Ancestry.** 18 of 35 GWAS strata are non-European against a European eQTL panel. The 13-gene "four of five ancestry panels" matrix (RESULTS.md:1121-1124) tests cross-population sharing under one LD reference, not ancestry replication. State the LD panel per stratum in Methods and label the matrix accordingly. [Unverified which LD panel each stratum used]

**13. Spatial.** Four donors plus ten arrays without a donor key; the ECM effect falls from 0.139 to 0.041 after composition adjustment; an ECM program organized in fibrotic liver is expected. The result is a methods point (matched nulls are needed to read Moran's I), fine as a supplement, not a "flagship".

**14. Proteomics.** One cohort, 58 patients, panel 5C selection-conditioned, RNA-protein concordance conditioned on 211 significant proteins. State the conditioning in the legend and do not call it validation.

**15. Figure 1 leads with what could not be tested.** 12/117 CosMx, 0/117 GeoMx, 38,066 of 39,914 ATAC peaks indeterminate, zero tested-negative calls, as the reader's first quantitative encounter. Move to Figure 6 or S1.

**16. Two DEG thresholds with the calibrated one demoted.** The canonical threshold changed on 2026-08-12 to padj plus fold change "for comparability", and RESULTS.md:474-478 concedes it has no calibrated error rate. fig3.md:77-80 then uses the TREAT arm to keep MAT1A in the worked example. Pick one primary and never let the sensitivity arm rescue a headline gene.

**17. Aspect resolution.** The three-component negative (participation ratio 1.42) is one cohort; SAF grades in GSE202379 and PXD051911 sit unused. Either test them or drop "MASLD is not one axis" from the eight linked ideas, because the paper's own evidence says it is one axis.

### Polish (severity 1 to 2)

- Kendall et al. 2023 Nature Medicine and LiRNA 2026 are in the literature list and not in the intro. Cite them and say what they did not do.
- The in-background DEG count on the adopted release is 1,066; 1,079 is the July value (RESULTS.md:392).
- "Across Figures 1-6" in legends and Discussion.
- Title length: 24 words, three clauses.

## Strategy for the highest realistic journal

2025 impact factors (JCR, released June 2026; from web search this session, official-looking sources; recheck before submission):

| journal | IF |
|---|---:|
| Nature Medicine | 52.5 |
| Journal of Hepatology | 40.1 |
| Cell Metabolism | 37.0 |
| Nature Methods | 28.3 |
| Nature Metabolism | 27.5 |
| Nature Genetics | 25.5 |
| Gut | 24.6 |
| Nature Communications | 18.1 |
| Hepatology | 18.0 |
| Genome Medicine | 10.8 |
| Genome Biology | 9.2 |
| Cell Genomics | 9.0 |

The current plan (Cell Genomics first) targets the bottom of this table. Journal of Hepatology now outranks Nature Genetics and Nature Metabolism, and it publishes multi-cohort human MASLD transcriptomic and genetic work. That is the target if the paper is rewritten as a hepatology paper with a patient-level message; as a "Resource and analytical framework" it will be desk-rejected there.

**Committed recommendation: two papers.**

Paper 1, the biology, to Journal of Hepatology, fallback Nature Metabolism, then Nature Communications.

One-sentence claim: inherited MASLD risk and established-disease remodeling are separable layers: risk concentrates in protein-altering variation (22% of posterior mass at direct-disease loci versus 3% at liver-enzyme loci), whereas tissue remodeling is one diffuse severity axis that a linear bulk-RNA instrument captures at the reproducibility ceiling of histology itself (0.68 to 0.79 of ceiling in 585 external participants).

Title candidates: "Inherited risk and tissue remodeling are separable layers of human MASLD"; "Protein-altering variation and a diffuse severity axis partition human MASLD biology".

Figure plan:
1. Cohorts, phenotypes, and the two layers (current Figure 1 without the observability counts).
2. Inherited variation is coding-heavy at direct-disease loci (Figure 2D as the lead; ancestry and consequence context around it).
3. Established disease is one diffuse axis: five-cohort DEGs, corrected adjacent stage contrasts (disease-only F0), the severity instrument and its external transfer, and the Kleiner test-retest ceiling (sigma_e 0.47 to 0.64 stages) as the bound.
4. The two layers intersect at chance and carry different information (the null, stated once, plus a test that can fail: cell-type expression or protein response by class).
5. Donor-level programs and, if it reproduces in the sealed GSE193066 paired biopsies, the fibroblast-cholangiocyte coordination community; the vocabulary null in the supplement.
6. Gene Catalog with rebuilt genetics fields, one worked entry.
Cut or demote: Figure 4G/4H and window heatmaps, most of Figure 5 (keep 5B and one spatial calibration panel), the seven-state accounting as a main-figure element, the GNMT bridge as a narrative device.

Paper 2, the method, to Nature Methods (Brief Communication or Analysis), fallback Genome Biology.
Claim: a model score is uninterpretable without a per-cohort record of what the model already saw; with that control, pretrained sequence and single-cell models add nothing to a linear task-native baseline on tissue-level MASLD severity, and a nonlinear model without biological pretraining adds 0.000 ± 0.067. Ship the severity instrument, the exposure audit, and the benchmark. This reverses the 2026-08-30 "one paper" decision; the reason is that the two audiences do not overlap and each half is stronger alone.

Required analyses on existing data, in order:
1. Rebuild the Discussion and Figure 1B on the adopted release; delete the July claims (days).
2. Rerun adjacent fibrosis with disease-only F0 and NAS with per-cohort effects (days).
3. Coding-mass robustness: leave-one-study-out, LD-block clustering, window-size sensitivity, as HIGH_IMPACT_PATTERNS row 7 already lists as owed (one to two weeks).
4. GSE193066 paired biopsies: prespecify and test whether within-person change in the severity instrument tracks within-person change in Kleiner stage across the 58 pairs. This is the only longitudinal substrate in the project and it is used today only for a test-retest statistic (one week).
5. Frozen-community external test for the coordination result in GSE193066 and the 41 paired liver-plasma donors, without refitting (one to two weeks). If it fails, delete rows 2 and 3 of HIGH_IMPACT_PATTERNS.md.
6. SAF aspect test in GSE202379 and PXD051911, prespecified (one week).

New data or experiments: none required for Journal of Hepatology or Nature Metabolism. Nature Medicine or Nature Genetics would need a new cohort with clinical outcomes or a resolved locus with functional evidence, which the two-paper firewall places in the Cas13 paper.

Risks: Journal of Hepatology may still send a public-data-only paper to a lower tier; Nature Metabolism published Kamzolas and may read Paper 1 as a follow-on unless the diffuse-axis contrast with their systems claim is explicit; the reframing gives up the "Resource" identity the documentation was built around; splitting adds three to six months.

Presubmission pitch (Paper 1, under 150 words): Human MASLD genetics and tissue transcriptomics have been analyzed separately and are usually assumed to point at the same genes. Across 35 GWAS strata, 1,257 liver RNA-seq samples from nine cohorts, and liver proteomics, we show they are separable layers. Inherited risk at direct-disease loci concentrates in protein-altering variants (22% of fine-mapped posterior mass versus 3% at liver-enzyme loci), while established disease is one diffuse transcriptional severity axis: a linear instrument fitted in five cohorts transfers to 585 external participants at 68 to 79% of the ceiling set by histology's own repeat-biopsy reproducibility, and no foundation model improves on it. The two gene sets intersect at the rate expected by chance. We release the instrument, the harmonized cohorts, and a gene catalog that records what each assay could and could not test.

## What would raise the ceiling

| action | tier unlocked | cost |
|---|---|---|
| paired-biopsy change analysis in GSE193066 (58 pairs) | secures J Hepatol / Nat Metab | one week, data in hand |
| coordination community replicated in a frozen holdout | Nat Metab / Cell Metab biology headline | one to two weeks |
| an outcome-linked cohort (events, FIB-4 trajectories, treatment response) | Nature Medicine | new collaboration, months |
| one direct-disease locus resolved to variant and function | Nature Genetics | belongs to the Cas13 paper by the firewall |

## Coverage: what this review did not check

The 12-lens workflow (wf_cb787e44-053) was stopped by the session usage limit after one lens (headline logic) completed; the adversarial verification, the three strategists, and the judge did not run. Not re-derived here: COLOC posteriors and LD panels per stratum; the continuum competitive nulls; the spatial matched nulls; the ambient audit; the coordination result's permutation code; Figure 7 numbers beyond the summary header; literature claims (Ji 2012 is cited in the draft; Mato/Lu, Mostafavi 2023, Kendall 2023 are from memory). The workflow can be resumed after the limit resets to produce the verified version.

## Method

Read the five authorities, the drafts, the figure authorities, and the benchmark summary; recomputed the overlap null, the F0 and NAS cohort-by-group crosstabs, and the absence of the Discussion numbers from the ledger; launched a 12-lens review with per-finding adversarial verification and a three-strategist panel, of which one lens completed; checked 2025 impact factors online. Nothing in the repository outside this directory was changed.
