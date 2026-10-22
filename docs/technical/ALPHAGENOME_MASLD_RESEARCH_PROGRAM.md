# An integrated MASLD regulatory model and database

## Research program: sequence, cellular state, and inherited molecular effects

**Requested:** 2026-09-14.  
**Expanded:** 2026-09-15; revised against the Atlas, regulatory program and model-benchmark results, including the completed signed-splice evaluation.  
**Status:** initial one-week execution authorized 2026-09-15; development candidates only, with protected outcomes closed. T01–T24 and the six campaigns remain the scientific specification.  
**Literature expansion:** 2026-09-15; 30 annotated references and 16 bounded comparisons, with prior MASLD applications checked.  
**Working model name:** MASLD regulatory response model.  
**Governing decisions:** strengthen the existing atlas paper and preserve its overall structure; expand through suitable public data; deliver a usable public resource with source-specific data and model terms. The [2026-09-07 owner decision](../decisions/2026-09-07-derivative-weight-release.md) permits GSE267145/GSE296875 derivative weights, while source data remain unredistributed. AlphaGenome derivatives retain the applicable noncommercial terms. New experiments are not manuscript dependencies.

This program extends the scientific ambition beyond an AlphaGenome annotation exercise or an MPRA leaderboard. It builds on the completed analyses in the [September 9 protocol](MASLD_ALPHAGENOME_RESOURCE_PLAN_2026-09-09.md) and [working results synthesis](alphagenome_atlas_RESULTS_2026-09-11.md). Those records and their prespecified hypotheses remain intact. Current manuscript claims and numbers remain governed by [PAPER](../PAPER.md), [STATUS](../STATUS.md), [RESULTS](../RESULTS.md), and [ROADMAP](../ROADMAP.md).

**Start here:** [evidence and dispositions](#32-integrated-evidence-and-next-action-dispositions), [public-data expansion](#51-public-data-expansion-and-protected-evaluation), [six scientific goals](#171-six-scientific-goals), [24 tasks](#172-task-map), [parameter search](#18-how-to-search-this-space), [six campaigns](#19-campaign-sequence-and-deliverables), [database and public models](#20-the-masld-regulatory-database-and-public-models), and [orchestrator handoff](#21-one-orchestrator-and-bounded-scientific-assignments), and [literature-informed experiments](#22-literature-informed-expansion-and-experiments-worth-testing). T01–T24 retain their identifiers and parameter menus. Each task states its goal, supporting results, unresolved question, data, evaluation unit, variables, comparator, success criterion, dependencies, outputs and consequence for the paper.

The new work is intended to strengthen this atlas paper. Its scientific scope is defined by the questions below, and its eventual inclusion by the evidence obtained. It need not wait for submission to begin. The Resource remains independent of the future pooled mouse screen; that screen's outcomes cannot select this paper's genes, programs, thresholds, or claims. The subsequent user instruction authorizes bounded implementation and public-source acquisition under section 19.4; controlled data, new experiments, final adoption and publication are outside that authorization.

## 1. The objective and the result worth pursuing

**Central question:** Which inherited regulatory effects can we predict from sequence, which depend on cellular context, and how do those effects relate to the multicellular remodeling observed in human MASLD?

The model and database connect **variant → molecular effect → candidate target → cellular context → disease-associated state**, with the evidence supporting each connection visible. This is a set of separately testable relationships; joining the records does not establish a causal chain. The scientific ambition is a reusable account of human liver regulatory response, including the relationships that fail to transfer.

The intended result is a human MASLD atlas that links four levels of evidence:

1. **Inherited variation:** signal-specific genetic uncertainty, coding and noncoding mechanisms, and actual phased haplotypes where available.
2. **Molecular response:** assay-specific allele effects on accessibility, expression, RNA processing, and other measured regulatory outputs.
3. **Cellular context:** where an effect is observable, whether it changes across molecular states, and whether that change generalizes to new donors and studies.
4. **Disease interpretation:** how the supported molecular relationships connect to metabolism, inflammation, fibrosis, and the Resource's fixed programs, protein, and spatial observations.

[Inference] The combination could make a substantial contribution to both MASLD biology and computational genomics if the model predicts held-out molecular differences that existing approaches miss. Neither use of AlphaGenome nor a larger parameter count establishes that contribution.

The flagship deliverable should let a researcher submit a variant or local haplotype with an explicitly defined liver context and obtain molecular-effect predictions, calibrated uncertainty where supported, relevant measured evidence, and a clear statement of unsupported contexts. It should also reveal where predictions fail. A precise model of one demonstrably transferable biological relationship is more useful than an unvalidated claim to model all of MASLD.

**Three decisive demonstrations to pursue:**

- A model predicts signed endogenous variant effects at new loci and in an independent human source, beyond the strongest eligible specialist and simple comparator.
- Molecular context explains reproducible variation in allelic effects, beyond average cell-type profiles and sequence-only predictions.
- For an independently evaluated subset, predictions distinguish local reporter activity from endogenous response, or identify a reproducible haplotype-by-state interaction.

These are research objectives, not required positive results. The paper should explain the measured boundary if one of them fails.

## 2. Where the novelty would come from

Novelty can come from a new predictive method, a rigorous adaptation of an existing method to MASLD, or a useful new biological result made possible by that adaptation. A new architecture is not required. Methods developed in related metabolic or fibrotic diseases, cancer, immune disease, or other biological systems are eligible when they answer a specific MASLD question with suitable data and independent evidence. Usefulness and priority are separate: an impactful application can merit inclusion even if a literature check finds an earlier MASLD application.

The targeted literature checks on 2026-09-14/15 establish that the basic ingredients already exist; they are not an exhaustive priority search.

| Existing work | What is already established | What this program must add |
|---|---|---|
| AlphaGenome and Atlas | Multimodal sequence predictions and genome-wide precomputed single-nucleotide effects | New measured human MASLD knowledge and independently demonstrated transfer, beyond querying the public model |
| scooby | Adaptation of a sequence model with a cell-state-specific decoder | Tests of disease-associated allelic response, reporter-to-endogenous transfer, and held-donor/locus/study generalization in human liver |
| Corgi | Sequence plus trans-regulator expression for context-dependent functional profiles | Evidence that adaptation learns the molecular differences of interest in these human datasets, including failures against strong RNA-only models |
| AlphaGenome MPRA adaptation | Short-sequence encoder probing and fine-tuning for reporter assays | A measured connection between reporter effects and endogenous allele effects, rather than another reporter-only application |

Sources, revisited during planning on 2026-09-15: [AlphaGenome model](https://huggingface.co/google/alphagenome-all-folds), [Atlas description](https://blog.google/innovation-and-ai/models-and-research/google-deepmind/alphagenome-atlas/), [scooby](https://www.nature.com/articles/s41592-025-02854-5), [Corgi](https://www.nature.com/articles/s41467-026-75527-2), and [MPRA adaptation by Murphy, Durán and Koo](https://genomicsxai.github.io/blogs/2026-002/). The MPRA adaptation study is separate work, not a study by the AlphaGenome model's authors.

One proposed methodological contribution is an **assay-aware model of allelic response**, trained to preserve the distinction between shared regulatory activity, allele effects, and allele-by-context interactions. Its short-construct and endogenous-genome tasks may share sequence representations, but retain their own outputs and measurement models. Shared learning must earn its place against separate specialists. A well-supported MASLD adaptation of another method is an equally eligible contribution; its value is judged by the resulting biological insight or predictive capability.

The proposed biological contribution is a **measured map of transferable and context-dependent inherited regulation in human liver**. Calling it the first, best, or comprehensive requires an updated literature comparison and the corresponding empirical evidence. No such claim is made here.


### 2.1 Adapt methods from other diseases and biological systems

The rows below are proposed adaptations, not claims of novelty or completed MASLD findings. Scout both nearby disease settings, such as fibrosis, obesity and diabetes, and distant settings, such as cancer and neuroimmune disease. Prioritize the information a method adds over the similarity of the disease name. The original papers establish the source methods; the MASLD designs are proposals. The expanded [reading map](#222-annotated-reading-map) and [detailed comparisons](#223-detailed-comparisons-to-add-to-t01t24) give implementation choices and stopping rules. SCENIC+ already has a reported MASLD application (R29); its proposed value here is the additional allele-sensitive or independently validated capability.

| Candidate and original setting | MASLD goal and required adaptation | Decisive comparison and completion criterion | Tasks / campaign |
|---|---|---|---|
| [SCENIC+](https://www.nature.com/articles/s41592-023-01938-4): enhancer-driven networks from single-cell RNA and accessibility, including melanoma applications | Connect allele-sensitive regulatory regions to lineage-specific TF–region–gene relationships. Fit networks using eligible adult liver training donors, retain native RNA/ATAC models, then attach sequence-effect evidence without treating motifs as measured TF binding. | Compare proximity, co-accessibility and the same network without sequence evidence. Complete with independently supported target discrimination or reproducible network interpretation beyond expression/activity-matched controls. | T12, T16, T22 / B, D, E |
| [sc-linker](https://www.nature.com/articles/s41588-022-01187-9): single-cell programs, epigenomic SNP-to-gene maps and GWAS across multiple diseases | Test whether sequence-informed links clarify which fixed liver programs carry inherited evidence. Preserve fixed117 memberships/weights; use suitable genome-wide GWAS and ancestry-matched LD, not a selected credible-set table as a replacement for genome-wide inputs. | Compare the original linking strategy, proximity and sequence-augmented links with the same GWAS/LD/program universe. Complete with a replicated program relationship or a precise limit, conditional on baseline genomic annotations and co-expressed-set controls. | T16, T22 / D, E |
| [MIRA](https://www.nature.com/articles/s41592-022-01595-z): joint RNA/accessibility representations and local versus broader-state regulation | Ask whether paired-modality context explains donor regulation beyond RNA alone. Use training-only auxiliary representations and donor-aware evaluation; retain modality-native count models and exclude the evaluated target and relevant genomic features from context inputs. | Compare RNA-only, ATAC-only, additive and paired inputs on held donors/regions. Complete only if the receiving molecular endpoint improves with leakage controls; this does not reopen the closed atlas integration or redefine fixed117. | T13–T15 / B |
| [CellOracle](https://www.nature.com/articles/s41586-022-05688-9): regulatory-network perturbation modeling in development and hematopoiesis | Test whether inferred regulator responses explain measured MASLD-associated cellular states. Adapt lineage networks and evaluate against eligible public perturbation observations where available; no new experiments are required. | Compare expression-only regulator ranking and network-only predictions, then test whether allele-sensitive links add information. Without compatible public perturbation labels, return observational hypotheses only. Predicted state shifts do not establish longitudinal progression or treatment benefit. | T16, T22 / D |
| Context-conditioned sequence models such as [Corgi](https://www.nature.com/articles/s41467-026-75527-2) and [scooby](https://www.nature.com/articles/s41592-025-02854-5) | Adapt existing conditioning to adult liver donor differences and allele-by-state effects, including objectives and context inputs absent from the completed local comparison. | Beat sequence-only, RNA-only and additive models on the specified held-donor/locus task. A successful MASLD adaptation qualifies even if the conditioning architecture is established. | T06, T12–T15 / B, F |

For each additional candidate, write a short admission record before implementation: original disease and assay; precise MASLD question; previous MASLD applications; available public labels; required change to the original model; strongest eligible comparator; replication unit; useful improvement or biological decision; compute and reuse conditions. Assign it within T01–T24 rather than creating an unbounded parallel program. A method with no compatible endpoint remains a literature lead.

Keep the source method's measurement assumptions explicit. For example, an RNA count model and an accessibility count model need their respective likelihoods and library-size treatment; a methylation-variance method cannot be applied to RNA counts merely because its API accepts a matrix. If the modality changes, specify and validate the new observation model and null before using the result. Independent biological support, rather than a more elaborate network plot, determines whether an adaptation advances the paper.

## 3. Starting evidence: retain the positives and the negatives

The original starting observations below are retained beside the newer executed results in section 3.1. They are development evidence, not independent confirmation of the proposed context-dependent model. Arithmetic for the two HyenaDNA gains and the Corgi gain was checked during initial planning; underlying training and bootstrap analyses were not rerun for this document.

| Observation | Consequence for the research design |
|---|---|
| HyenaDNA allele-difference head: MPRA macro Spearman 0.2162 versus allele-identity ridge 0.1342; gain 0.0820, block-bootstrap interval 0.0274–0.1389 | Sequence representations can add information on this development task. Stronger specialists and external evaluation remain necessary. |
| HyenaDNA using only the assayed 107-bp sequence: 0.2063; gain 0.0720, interval 0.0196–0.1271 | Useful signal survives within the reporter sequence. The 4,096-bp version's advantage over this arm is unresolved. This does not test endogenous long-range regulation. |
| Corgi FiLM-plus-head on 99 participants: held-tile participant-differential H3K27ac skill improvement 0.000158 | One context-conditioning setup failed. Its 99 distinct RNA context vectors were present; the model did not use them effectively. Optimize and evaluate molecular differences explicitly. |
| Accessibility carries measured signal; liver RNA direction is weak; unsigned Atlas splice scores could not test direction; the later signed model-API route recovers weak direction below its ranking threshold | Separate expression direction, splice-event detection, magnitude and direction. The completed signed-splice result is recorded in E10 below. |
| Existing within-locus prioritization and disease-state tests have limited independent loci or donors | Add independent biological information. Additional correlated rows are not a remedy for limited replication. |
| Existing TF-convergence and competitive module analyses do not establish the proposed broad specificity claims | Preserve these results. Do not change thresholds, programs, or controls until a desired story appears. A materially different hypothesis requires new evidence and a new prespecification. |

Evidence: [MPRA head comparison](../../Analysis/MASLD_Model_Benchmark/executions/model-cpu-train-605-21099008/evaluation/model_head_summary.tsv), [window ablation](../../Analysis/MASLD_Model_Benchmark/executions/gse281364-hyenadna-window-ablation-20260905T192600Z/evaluation/arm_metrics.tsv), [Corgi limitations](../../Analysis/MASLD_Model_Benchmark/executions/corgi-film-h3k27ac-finetune-20260907T191151Z/RESULTS.md), [Atlas working results](alphagenome_atlas_RESULTS_2026-09-11.md), and [TF-convergence analysis](alphagenome_atlas_execution_2026-09-09.md).

### 3.1 What the executed program changes, 2026-09-15

All values here remain candidate research results. The wider task menu does not amend a completed experiment's primary contrast or turn a secondary result into confirmation.

| Observation | Meaning for the expanded program | Direct source |
|---|---|---|
| Local AlphaGenome ATAC versus archived Atlas: Spearman 0.9965, 3,627 variants in 1,676 blocks; 218 of the original 3,845 windows excluded | Local inference is available. Agreement between implementations is not evidence of unseen-donor or unseen-genome generalization. | [Local scoring](../../GWAS/finemapping/results/alphagenome_program/ad1-gate-20260915T114615Z/RESULTS.md) |
| Local zero-shot ATAC fold-macro Spearman 0.7441 versus HyenaDNA 0.2344 on the same 3,627 variants | Adaptation must be compared with the strong existing local score, as well as HyenaDNA. Neither number describes the full 32,322-variant universe. | [Matched model scores](../../GWAS/finemapping/results/alphagenome_program/ad1-gate-20260915T114615Z/tables/zeroshot_arms.tsv) |
| Frozen AlphaGenome reporter primary 0.2198 versus HyenaDNA 0.2162; difference +0.0036, 95% CI −0.0587 to +0.0635 | No primary superiority or equivalence established. | [Reporter probe](../../GWAS/finemapping/results/alphagenome_program/ad-arm2-mpra-20260915T114407Z/RESULTS.md) |
| Reporter variant-bin secondary 0.3292; difference versus HyenaDNA +0.1130, 95% CI +0.0537 to +0.1714 | Pooling and representation deserve a substantial development search. Independent confirmation remains necessary. | [Reporter configurations](../../GWAS/finemapping/results/alphagenome_program/ad-arm2-mpra-20260915T114407Z/tables/cell_scores.tsv) |
| Currin training leads: median distance to their own peak centre 198 bp; 23,218/32,322 within 500 bp | The 11.3-kb median in adaptation amendment 02 belongs to the reporter-overlap set. It cannot justify long context for these peak leads. The 198-bp calculation is the median of `abs(pos_hg38 − (peak_start_hg38 + peak_stop_hg38)/2)` over the deposited rows. | [Currin training labels](../../GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/inputs/currin_lead_labels.tsv.gz) |
| HyenaDNA's final training-size increment is +0.0173 in mean seed-level macro Spearman; the 2,000-locus macro exists for only 1/5 seeds | Explore label size and representation together; do not fit a smooth complete learning curve through undefined statistics. | [Learning curve](../../GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/tables/learning_curve.tsv) |
| Projected treatment-dependent allele-effect reliability: HepG2 0.0427; LX-2 0.7446 | Start state-interaction development with LX-2. Improving HepG2 interaction measurement is a separate data task. | [Replicate reliability](../../GWAS/finemapping/results/alphagenome_program/a1-mpra-reliability-v2-20260914T204500Z/tables/decision_rule_section3.tsv) |
| Reporter prediction added to the endogenous head: −0.0015, 95% CI −0.0070 to +0.0039, on 1,244 variants in 168 blocks | This one-column transfer route did not help. Joint training, encoder adaptation and better-matched labels remain distinct experiments. | [Transfer experiment](../../GWAS/finemapping/results/alphagenome_program/e2-transfer-20260914T203752Z/RESULTS.md) |
| Existing lineage allelic coverage: hepatocyte 292 sites/50 blocks; endothelium 64/23; stellate 27/15; cholangiocyte 3/3, over 30 donors | Build richer adult profile heads where donor data permit; acquire or recover more allele labels for sparse lineages before claiming lineage-specific allele prediction. | [Lineage coverage](../../GWAS/finemapping/results/alphagenome_program/d1-lineage-allelic-20260914T194213Z/tables/lineage_feasibility.tsv) |

Two interpretations from the acquisition proposal should not constrain this expansion. The post hoc Pearson reliability ceiling of 0.984 is explicitly optimistic because the caQTL leads were significance-selected; it is not a bound on the benefit of additional donors. A ridge penalty at the top of its grid is a reason to explore regularization, not proof that widening the representation will improve prediction. Independent donor subsampling requires individual-level inputs; summary statistics alone do not provide an empirical donor learning curve. See the [label-precision caveat](../../GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/tables/label_precision_posthoc.json).

The acquired local weights remove the earlier access obstacle. They do not remove the distinction between all-folds pretraining and a genuinely held genomic region, or between local outputs and hosted outputs. Cross-process numerical differences were also observed on the same GPU architecture; small allele effects need tolerance-based reproducibility checks, not a blanket bit-identical claim.

### 3.2 Integrated evidence and next-action dispositions

This table covers the original Atlas A–H packages and their follow-up analyses, the executed A–G regulatory program, local runtime/adaptation experiments, and the relevant sequence, chromatin, severity and module benchmarks. It is an experiment-level synthesis, not a claim that every archived read or bootstrap has been independently rerun for this document. The supporting reports, small result tables, score metadata and relevant code were inspected; scientific reruns were not performed.

**Evidence snapshot:** 2026-09-15, including the completed step-36 signed-splice result. Result availability is stated below; live job status is not inferred from dated reports. An E-number identifies evidence in this document, not an independent dataset or a new hypothesis family.

Completed work receives one or more explicit dispositions: **build on**, **correct**, **retain as a limitation**, or **revisit with a materially different design**. Missing outputs receive **unfinished**. Changes to a completed experiment's recipe are new development analyses and retain the original result beside them.

#### Atlas coverage, molecular effects and biological interpretation

| ID / analysis | Result and scientific unit | Disposition and next action | Source |
|---|---|---|---|
| E01. Track 0; packages A–H | 417,155 scored variants across 989 signals; 12,857 channel-profile rows and a 989-row Catalog sidecar. Multiple channels and repeated signal–variant rows are not independent evidence. | **Build on.** Reuse the existing molecular landscape, chromatin, RNA-processing, direction, saturation, tissue/ancestry, context and predictor-interpretation products; attach new capabilities rather than repeat annotation. | [Atlas synthesis, coverage](alphagenome_atlas_RESULTS_2026-09-11.md#5-coverage-delivered); [original A–H protocol](MASLD_ALPHAGENOME_RESOURCE_PLAN_2026-09-09.md#6-analysis-packages) |
| E02. RNA direction; package D | Liver eQTL-direction auROC 0.5304, interval 0.449–0.607, on 201 leads in the reconstructed benchmark. Other direction tables use different populations and track panels. | **Retain as a limitation; revisit with a materially different design.** Test target/event-aware readouts and measured context; do not substitute accessibility accuracy for expression accuracy. | [Atlas synthesis, benchmark definitions](alphagenome_atlas_RESULTS_2026-09-11.md#1-what-the-atlas-can-and-cannot-do-in-liver) |
| E03. Accessibility caQTL benchmark | Atlas auROC 0.7128 on 10,966 common leads; archived API comparison 0.71454 on those same leads. The 138 donors are the source caQTL cohort, not 10,966 independent participants. | **Build on.** Preserve the matched zero-shot benchmark and its source-exposure limitations; separate auROC from the later signed-Spearman endpoints. | [Atlas caQTL comparison](alphagenome_atlas_RESULTS_2026-09-11.md#the-caqtl-arm-reproduces-the-archive-exactly) |
| E04. Same-variant allelic accessibility | GSE281367: concordance 0.6605 at 430 sites in 48 blocks; GSE244832: 0.5668 at 277 sites in 44 blocks. These are distinct from later expanded matched-model populations. | **Build on.** Preserve donor-aware allele measurements and source-specific results; use them to test whether caQTL learning transfers to an endogenous allelic endpoint. | [Same-variant evidence](alphagenome_atlas_RESULTS_2026-09-11.md#2-the-new-evidence-same-variant-allelic-concordance) |
| E05. Disease-state allelic tests | Existing global tests did not show a state interaction; the small donor design could not support site-level inference. | **Retain as a limitation.** Recover more eligible loci from existing reads and seek public donors with documented state variation; synthetic substitutions do not add donors. | [Atlas claim boundaries](alphagenome_atlas_RESULTS_2026-09-11.md#3-what-is-not-established); [state analysis code](../../scripts/analysis/alphagenome_atlas/74_state_dependence.py) |
| E06. Allelic calibration | A magnitude trend was supported in one cohort and unresolved in the other. A high genome-wide quantile does not establish a large native effect or exceptional within-locus rank. | **Retain as a limitation; build on.** Report native effect, calibration population and local rank separately; evaluate calibration in a new compatible source. | [Calibration and score interpretation](alphagenome_atlas_RESULTS_2026-09-11.md#4-two-constraints-on-presenting-these-scores) |
| E07. Within-locus prioritization | 57 sets collapse to 11 independent blocks; intervals include chance and posterior ranking is not distinguishable from Atlas ranking. | **Retain as a limitation.** Increase independent informative loci and compare within-locus discrimination against genetic posterior and measured-region baselines. | [Prioritization analysis](../../scripts/analysis/alphagenome_atlas/81_prioritisation.py); [execution record](alphagenome_atlas_execution_2026-09-09.md) |
| E08. Anchors and target directions | Known anchors help check channels and orientation. Only three drug-target directions were available; no concordance rate was supported. Extreme quantiles can accompany small raw changes. | **Build on** as implementation checks; **retain as a limitation** for accuracy or therapeutic claims. Known anchors cannot be the independent test of a model selected to recover them. | [Atlas synthesis](alphagenome_atlas_RESULTS_2026-09-11.md); [anchor code](../../scripts/analysis/alphagenome_atlas/21_anchor_panel.py); [target-direction code](../../scripts/analysis/alphagenome_atlas/60_target_direction.py) |
| E09. Unsigned splicing; step 35 | Atlas SPLICE_JUNCTIONS and SPLICE_SITE_USAGE are unsigned. The earlier direction statistic was invalid. Target-intron quantile rank is 0.4873 versus proximity 0.6053, on 657 pairs in 255 blocks; contrast ranking uses 583 pairs in 228 blocks. | **Correct** the direction interpretation; **retain as a limitation** for event detection. Keep event detection, magnitude and direction as separate endpoints. | [Corrected step-35 table](../../GWAS/finemapping/results/alphagenome_atlas/p3a-splice-detection-20260915T142531Z/tables/splice_detection.json); [scorer declarations](../../GWAS/finemapping/results/alphagenome_atlas/run-20260909T153939Z/tables/scorer_metadata.tsv) |
| E10. Signed splicing; step 36, now complete | Model-API excision-ratio direction is 0.5777 versus marginal 0.5004; interval 0.5405–0.6152; 637 pairs, 247 blocks, block p=0.000546. It misses both the prespecified P1 threshold (>0.65, excess >0.05, p<0.01) and the inconsistent deposited D1 wording (>0.60, p<0.05); preserve that inconsistency. Magnitude correlation is 0.0357. Paired event detection remains below proximity and shows no gain over Atlas on 617 pairs in 242 blocks. | **Build on** the weak directional signal; **retain as a limitation** for ranking and magnitude. **Correct** the deposited decision text saying “shows no direction.” A tighter, independently defined event denominator is a new development comparison; it does not erase D1. | [Completed signed-splice results](../../GWAS/finemapping/results/alphagenome_atlas/p3a-splice-direction-20260915T141415Z/tables/splice_direction.json); [step-36 code](../../scripts/analysis/alphagenome_atlas/36_splice_direction.py) |
| E11. Hosted-model MPRA benchmark | The reported HepG2 H3K27ac auROC 0.685 uses 3,498 oligos. The separate liver_h3k27ac channel has zero tracks because its filter searched the wrong metadata field. AVI sign rows are not valid allele-direction tests. | **Build on** the correctly defined channels; **correct** the empty liver channel through a declared amendment and archived-data re-read. This does not invalidate the different HepG2 headline channel. | [Execution corrections](alphagenome_atlas_execution_2026-09-09.md); [MPRA evaluator](../../scripts/analysis/alphagenome_atlas/33_benchmark_eval.py) |
| E12. TF convergence | Across 761 signals, 18,233 variants and 502 factors, the prespecified general-machinery control ranks above the hepatocyte-regulator set. Earlier null and per-factor interpretations were corrected. | **Retain as a limitation.** No MASLD-specific TF convergence is established. A new context-specific hypothesis needs an appropriate active-region control and independent evidence. | [TF-convergence execution](alphagenome_atlas_execution_2026-09-09.md#tf-convergence-step-80-a-clean-negative-delivered-by-the-prespecified-control); [producer](../../scripts/analysis/alphagenome_atlas/80_tf_convergence.py) |
| E13. Ancestry and missing mass; P6b | 623 signals covered, 353 partial_indel, 13 uncovered_indel. Indels account for the main missing posterior mass. Core channels serve the queried SNV mass across ancestry strata. | **Build on.** Carry both row- and signal-level coverage. Coverage parity is not ancestry-specific predictive accuracy; no performance claim follows without measured labels. | [Coverage analysis](alphagenome_atlas_RESULTS_2026-09-11.md#6-the-coverage-gap-is-a-variant-class-not-an-ancestry-p6b-2026-09-13) |
| E14. Indel identities and orientation; steps 67–69 | 114,389 source IDs represent 96,143 distinct deferred variants. The final dbSNP procedure resolves 112,863 IDs; independent source-order checks support its resolution rules. Ambiguous and assembly-conflicting cases remain. | **Build on; correct** consumers of earlier orientation columns. Preserve source IDs and normalized alleles; identity recovery does not make an indel Atlas-served. | [Final orientation evidence](alphagenome_atlas_RESULTS_2026-09-11.md#7-the-indels-are-recoverable-through-the-model-api-steps-63-64-64b-68-68b-2026-09-14) |
| E15. Indel rescue; steps 63–64b | 156 selected variants scored with 526 SNV controls; paired comparison: 142 variants in 97 blocks. Larger predicted effects are model outputs in native track units, not measured biological superiority of indels. | **Build on.** Expand coverage and event validation using local sequences. HSD17B13's predicted splice magnitude is not validated by the SNV direction result in E10. | [Final rescue summary](../../GWAS/finemapping/results/alphagenome_atlas/p6f-indel-rescue-20260914T231617Z/tables/indel_rescue_summary.json); [paired analysis](../../scripts/analysis/alphagenome_atlas/64_indel_vs_snv_paired.py) |
| E16. Direct-MASLD dossier; steps 65–66 | 1,422 rows over 202 direct signals, 951 readable under its coverage rule. Only 1.9% overlap measured peaks; the direction check has 19 blocks and does not establish a directional claim. | **Build on.** Preserve variants at incompletely covered signals and show missing mass explicitly. Peak non-overlap is not a negative test of regulatory function. | [Final dossier](../../GWAS/finemapping/results/alphagenome_atlas/p6d-dossier-20260915T000218Z/tables/clinical_variant_dossier_summary.json); [Atlas interpretation](alphagenome_atlas_RESULTS_2026-09-11.md) |
| E17. Full P4 saturation | Registered universe: 96,460 H3K27ac regions plus 6,930 unique snATAC regions, 103,390 total. Only feature pilots, explicitly not results, were present at this revision. U2c contributed no genetic-overlap peaks under the original construction. | **Unfinished.** Reconcile completion of retrieval and existing analyses before submitting dependent work. Preserve original F1–F5 questions; explicitly flag the empty genetic-overlap arm rather than manufacturing it from a new selection. | [Saturation universe and execution](alphagenome_atlas_execution_2026-09-09.md); [feature reduction](../../scripts/analysis/alphagenome_atlas/43_saturation_features.py); [tests](../../scripts/analysis/alphagenome_atlas/44_saturation_tests.py) |
| E18. Package H and chromatin-predictor interpretation | Region covariates connect the released RNA→chromatin predictor to sequence features. The shipped reliable-region mask and annotation conventions require their own definitions; region-level joins do not license base-resolution claims. | **Build on.** Use the actual shipped-form mask and matched covariates in P4; preserve signal-strength and technical explanations when interpreting predictability. | [Package H protocol](MASLD_ALPHAGENOME_RESOURCE_PLAN_2026-09-09.md#h-biological-interpretation-of-the-released-chromatin-predictor); [producer](../../scripts/analysis/alphagenome_atlas/18_package_H_region_covariates.py) |
| E19. Original P5 haplotypes | 1,572 pair rows produced no BH calls, but the 300-draw null had insufficient resolution for sparse discoveries in that family. Predicted residuals are scale- and window-dependent. | **Correct** any interpretation as evidence of additivity. Carry the repaired F2/F3 results in E30 and require measured combinations for a biological interaction claim. | [Original haplotype summary](../../GWAS/finemapping/results/alphagenome_atlas/p5-haplotypes-20260910T004716Z/tables/haplotype_summary.json); [later program corrections](ALPHAGENOME_PROGRAM_RESULTS_2026-09-14.md#64-f1-phase-substrate) |

#### Executed regulatory program and local AlphaGenome

| ID / analysis | Result and scientific unit | Disposition and next action | Source |
|---|---|---|---|
| E20. A1 measurement reliability | 4,359 elements in 239 long-range blocks; treatment-interaction projected reliability 0.0427 in HepG2 and 0.7446 in LX-2. HepG2 PAOA allele-effect reliability is 0.6938 (interval 0.6356–0.7383), which contains 0.70. Experimental aliquots are not donors. | **Build on** reliable allele effects and LX-2 interactions; **retain as a limitation** the current HepG2 interaction endpoint. Public datasets with better documented contrasts are the expansion route. | [A1 results](../../GWAS/finemapping/results/alphagenome_program/a1-mpra-reliability-v2-20260914T204500Z/RESULTS.md) |
| E21. A2 access and exposure | Initial 403/no-runtime statements were true on September 14 but superseded by E32–E34. Hosted outputs and local-weight derivatives have different training routes. Dataset-level restrictions remain distinct. | **Correct** current summaries; **build on** the existing exposure inventory. Reconcile local weights, input data and derivative distribution separately. | [Original A2 inventory](../../GWAS/finemapping/results/alphagenome_program/a2-rights-exposure-20260914T193626Z/RESULTS.md); [adaptation amendment](../../scripts/analysis/alphagenome_program/ADAPTATION_ARM_AMENDMENT_01.md) |
| E22. B1/B2 and dbSNP refresh | 989 signals, 10,511 signal–variant rows; 3,822 readable rows and 623 covered signals are different denominators. Measurement distinguishes no explanation class at 734 signals. Life-stage filtering changes some RNA signs. | **Build on** the response layer; **correct** only the superseded indel columns using the refresh. Keep adult, pediatric, embryonic and transformed contexts explicit. Missing target/event evidence motivates T16–T17. | [B1/B2 results](../../GWAS/finemapping/results/alphagenome_program/b1-response-layer-v2-20260914T204700Z/RESULTS.md); [orientation refresh](../../GWAS/finemapping/results/alphagenome_program/b1-dbsnp-refresh-20260915T084700Z/RESULTS.md) |
| E23. B3 bounded saturation | 3,000 regions in 1,337 blocks from a chromosome-truncated U1 snapshot. CHIP_HISTONE Q5−Q1 sensitivity +0.00166, q=0.037; continuous slope unresolved; AVI reverses. Promoter/TSS expectations fail in direction. | **Retain as a limitation; build on** the controls. Do not treat this as completed P4 or a monotone sequence-sensitivity law. Examine concentration and signal abundance separately in full P4. | [B3 results](../../GWAS/finemapping/results/alphagenome_program/b3-saturation-20260914T193641Z/RESULTS.md) |
| E24. C1 matched models, including Borzoi | ChromBPNet−Atlas signed-Spearman contrast −0.0910 at 3,845 caQTL leads; −0.1481 on the weighted 495-site allelic endpoint. Track averaging and pretraining favor Atlas in this comparison. | **Build on** the stronger comparator. **Retain as a limitation** the comparison's asymmetries and allele-identity signal on the allelic endpoint; do not claim unconditional model superiority. | [Completed four-family comparison](../../GWAS/finemapping/results/alphagenome_program/c1-endpoint2-v2-20260915T084559Z/RESULTS.md) |
| E25. C2/C3 open heads and learning curves | HyenaDNA macro Spearman 0.1737 at 32,322 leads in 2,643 blocks; gain +0.1301 over allele identity. Final training-size increment +0.0173; some small-size macro statistics are undefined. | **Build on.** Cross label breadth with representation and trainable depth. Keep the selected-label reliability ceiling qualified; it cannot establish that donor expansion has no value. | [C2 results](../../GWAS/finemapping/results/alphagenome_program/c2-endogenous-head-20260914T194000Z/RESULTS.md) |
| E26. D1 lineages and mapping bias | 30 donors; hepatocyte 292 sites/50 blocks, endothelial 64/23, stellate 27/15, cholangiocyte 3/3. Matched remapping controls narrow but do not remove genotype-ascertainment uncertainty. | **Build on** hepatocyte labels; **retain as a limitation** sparse-lineage allele claims. Expand eligible measured sites and public genotype support rather than treating more cell rows as more donors. | [D1 results](../../GWAS/finemapping/results/alphagenome_program/d1-lineage-allelic-20260914T194213Z/RESULTS.md) |
| E27. E1/E2 bridge and transfer | 3,894 reporter variants have nominal caQTL tests, but only 209 are within 500 bp of the tested peak center. Added reporter-OOF feature changes correlation by −0.0015, interval −0.0070 to +0.0039, on 1,244 variants/168 blocks. | **Retain as a limitation** of this route; **revisit with a materially different design** using matched labels and shared learning. The positive-control ceiling is a property of that simulation/recipe, not a universal transfer bound. | [Bridge census](../../GWAS/finemapping/results/alphagenome_program/e1-bridge-20260914T194145Z/RESULTS.md); [transfer results](../../GWAS/finemapping/results/alphagenome_program/e2-transfer-20260914T203752Z/RESULTS.md) |
| E28. E3 reporter sequence background | Hosted-model native flanks did not reliably improve on N-padding in the evaluated random and selected-effect strata. API acceptance of N-padding says nothing about biological adequacy. | **Retain as a limitation; revisit with a materially different design.** Compare documented constructs, encoder adaptation and endogenous backgrounds within their own assay tasks. | [E3 results](../../GWAS/finemapping/results/alphagenome_program/e3-reporter-context-20260914T193515Z/RESULTS.md) |
| E29. F1 phase | The reference-panel pairs commonly have high LD, but this is not donor-specific read-backed phase. The GNMT single-variant arm occurs on only 2/758 reference-panel EUR haplotypes. | **Build on** population phase for hypothesis construction; **retain as a limitation** individual sequence claims and poorly represented combination arms. | [F1 results](../../GWAS/finemapping/results/alphagenome_program/f1-phase-20260914T193608Z/RESULTS.md); [program haplotype interpretation](ALPHAGENOME_PROGRAM_RESULTS_2026-09-14.md) |
| E30. Corrected F2/F3 | GNMT stays inside the repaired matched null. Broader apparent nonadditivity concentrates in short, overlapping readouts and attenuates with arm matching; pair rows also repeat across signals. | **Retain as a limitation.** Require distinct pairs, identical readout geometry, adequate null resolution and measured joint effects. No GNMT-positive result is required. | [Corrected F2/F3](../../GWAS/finemapping/results/alphagenome_program/f2-haplotype-v2-20260915T063000Z/RESULTS.md) |
| E31. G integration and acquisition proposal | Candidate Catalog attachment and complete correction register exist; every row remains candidate. Acquisition suggestions include experiments outside the public-data-only route chosen for this revision. | **Build on** the attachment; **correct** stale weight blockers and overgeneralized data ceilings. Convert missing-label priorities into public-source searches; unavailable labels remain explicit limits. | [Integration](../../GWAS/finemapping/results/alphagenome_program/g-integration-20260915T072500Z/RESULTS.md); [original data proposal](../../GWAS/finemapping/results/alphagenome_program/g-integration-20260915T072500Z/DATA_DEVELOPMENT_PROPOSAL.md) |
| E32. Local runtime H1 | Local weights and JAX inference work. All five standard input lengths run on L40S; 1-Mb inference uses about 17.6 GiB. Backward-pass memory is unmeasured. Cross-GPU differences reach 0.0388 log2 units on tested outputs. | **Build on.** Reuse the runtime; measure custom-head/loss training before scaling. Keep hardware fixed within small-effect comparisons. Static ontology selection does not supply a donor-context input. | [Runtime results](../../Analysis/MASLD_Model_Benchmark/executions/alphagenome-runtime-probe-20260915T104500Z/RESULTS.md) |
| E33. AD.1 local scoring | Local/archive ATAC correlation 0.9965 on 3,627 leads/1,676 blocks; 218 windows excluded. Local zero-shot macro 0.7441 versus matched HyenaDNA 0.2344. | **Build on.** Establish full eligible-population baselines and measure length/coverage tradeoffs. Implementation agreement is not external validation. | [AD.1 results](../../GWAS/finemapping/results/alphagenome_program/ad1-gate-20260915T114615Z/RESULTS.md) |
| E34. AD reporter probe | Primary pool-region macro 0.2198 versus HyenaDNA 0.2162, unresolved difference. Secondary variant-bin macro 0.3292, paired gain +0.1130. | **Build on** the pooling hypothesis; **retain as a limitation** independent superiority. Keep primary and development comparisons identified and confirm selected recipes elsewhere. | [Reporter results](../../GWAS/finemapping/results/alphagenome_program/ad-arm2-mpra-20260915T114407Z/RESULTS.md) |

#### Other MASLD models and tissue-level evidence

| ID / analysis | Result and scientific unit | Disposition and next action | Source |
|---|---|---|---|
| E35. Earlier DNA-model MPRA and window ablation | HyenaDNA 0.2162 versus allele ridge 0.1342; the assayed 107-bp arm retains 0.2063. These are reporter effects, not clinical severity or endogenous expression effects. | **Build on.** Keep assay-native baselines and short-sequence controls. Reconcile block definitions before comparing historical DNA-model tables. | [Original head summary](../../Analysis/MASLD_Model_Benchmark/executions/model-cpu-train-605-21099008/evaluation/model_head_summary.tsv); [window ablation](../../Analysis/MASLD_Model_Benchmark/executions/gse281364-hyenadna-window-ablation-20260905T192600Z/evaluation/arm_metrics.tsv) |
| E36. Corgi adaptation | On 99 paired participants, the historical reconstruction reports differential skill about 0.00024, gain over zero-shot 0.000158, and RNA reconstruction 0.1186. Its scale was fitted to held outcomes; these are not training-calibrated prediction errors. The RNA comparator trained on every output region and is not equally region-held. | **Correct** the comparison's interpretation; preserve the failed configuration. New profile, residual and donor-difference fits learn preprocessing and calibration from permitted training data, and distinguish donor-held from donor × region-held evaluation. | [Corgi results and producer](../../Analysis/MASLD_Model_Benchmark/executions/corgi-film-h3k27ac-finetune-20260907T191151Z/RESULTS.md); [new comparison](../../scripts/analysis/alphagenome_campaign/review_donor.py) |
| E37. RNA–chromatin coupling and technical controls | RNA→chromatin predictability recurs in H3K27ac and same-nucleus ATAC. Disease dependence is not independently established across the two cohorts; adding chromatin residuals to severity can be explained by technical features. | **Build on** cross-assay predictability; **retain as a limitation** disease-specific and mediation interpretations. Match the exact target, rank and preprocessing of each comparator. | [Model results and corrections](../../Analysis/MASLD_Model_Benchmark/MODEL_RESULTS_SUMMARY.md); [current Figure 7 numbers](../RESULTS.md#model-framework-and-transfer-limits-figure-7-validated-candidate-working-figure) |
| E38. Severity models, exposure and learning | RNA severity transfer and its native ridge comparisons already exist. BulkFormer exposure includes the external RNA cohort; missing-gene sensitivity matters. RNA→H3K27ac learning curves have not plateaued, but sample count and component rank rise together in one recipe. | **Build on** the measured baselines; **correct** any claim of wholly unseen samples or a pure donor-size effect. Do not make clinical prediction the new regulatory model's success endpoint. | [Figure 7 ledger and learning curves](../RESULTS.md#model-framework-and-transfer-limits-figure-7-validated-candidate-working-figure) |
| E39. Cross-assay modules v2 | 160 modules from a 4,468-gene common universe; association replication and competitive specificity differ. Only one module-by-assay cell exceeds the co-expressed-set reference after correction. Within-donor coupling does not establish module-specific regulation. | **Build on** eligible scoring coordinates; **retain as a limitation** gene-list stability and specificity. Keep v2 corrections and fixed117 separate; no module member inherits its module's result. | [Cross-assay module ledger](../RESULTS.md#cross-assay-modules-v2-validated-candidate-2026-09-12-not-adopted) |
| E40. Fixed programs, continuum, physical context and integration limits | The 117 programs remain fixed. Continuum associations are cross-sectional; spatial/protein observability is assay-specific. Closed single-cell integration and external routing results do not support replacing the atlas representation. | **Build on** the current atlas, fixed projections and donor-aware tissue context; **retain as a limitation** causal, longitudinal and unsupported aspect-specific interpretations. No new regulatory task reopens closed integration or gene-set optimization. | [Scientific identity](../PAPER.md); [current evidence](../RESULTS.md); [closed integration and release boundaries](../ROADMAP.md) |

### 3.3 Corrections and dependencies that change the new design

- **Signed splicing is complete, not pending.** E10 establishes weak direction, below the prespecified ranking threshold, and poor event detection. The deposited D1 sentence overstates the negative; the numeric result and original decision threshold both remain. The sGene-lead secondary result (0.6854 on 89 pairs) cannot replace the full-set primary endpoint.
- **Exact duplicate correction preserves the splice conclusion.** E10's 637 scored rows contain 16 extra identical copies. The corrected 621 unique pairs retain 247 blocks: concordance 0.5765 versus marginal 0.5002, interval 0.5383–0.6135, original nominal block-majority P=0.000301, and magnitude Spearman 0.0284. Direction, uncertainty and magnitude were independently reconstructed in `week1-20260915/coordinator-splice-21773779/` under the campaign results root. The historical outputs and P1/D1 inconsistency remain; this changes repeated-row weighting, not the loose event definition or foundation-model exposure.
- **Tighter splice clusters are an untested explanation.** The current denominator contains a median 200 model-proposed junctions. [Inference] Restricting it using source-native clusters, annotations or a reference-only usage rule may improve estimand alignment. More predicted junctions alone do not prove dilution or guarantee better direction. Select rules in development and evaluate independently.
- **The empty liver-H3K27ac MPRA channel is a repair task.** Its existing HepG2 counterpart remains a separate valid channel. Record the metadata correction before a new analysis; use the archived intervals rather than a rescore that cannot reconstruct omitted tracks.
- **P4 is a dependency, not a completed biological result.** Distinguish U1, U2, the B3 partial sample and feature pilots. Confirm final region coverage and original test outputs before interpreting disease grammar or adding a new arm. Additional saturation predictions identify candidates; they do not create observed allele labels.
- **Use the final units.** The indel comparison has 142 paired variants/97 blocks in its final report, not the older 91-block summary. Distinguish 989 signals, 10,511 response rows, 3,822 readable response rows and 951 readable direct-dossier rows.
- **Evidence terms retain their meanings.** A caQTL beta is a marginal association and may tag another variant. Allelic imbalance is a donor read-ratio measurement affected by heterozygote ascertainment and mapping. A reporter allele effect is measured in its construct. An endogenous perturbation measures a different intervention. Models may share representations, but labels, claims and likelihoods retain these distinctions.
- **Historical constraints do not describe current implementation.** Local weights and inference exist; a custom training head, assay-specific loss and measured training memory are still needed. The signed scorer declarations and coordinate/readout conventions are inputs to every future comparison, not optional narrative notes.

## 4. Use the breadth of Atlas deliberately

Full use of Atlas means exploiting its different biological capabilities where they answer a question. It does not require downloading the entire petabyte-scale resource or treating every output as useful supervision.

| Capability | Research use | Evidence needed beyond the prediction |
|---|---|---|
| Single-variant effects across a signal | Compare competing molecular explanations while retaining variant posterior uncertainty | Compatible fine-mapping, measured allele effects, and locus-level evaluation |
| Interval saturation across possible substitutions | Map local sequence sensitivity; compare naturally observed variants with their local background; select diverse, discriminating candidates | Reporter or endogenous measurements; saturation alone does not measure redundancy or interaction |
| ATAC, DNase, histone and TF outputs | Identify candidate regulatory responses and compare them with adult liver measurements | Assay-matched data; accessibility does not establish TF occupancy or histone modification |
| RNA, CAGE and PRO-cap | Separate gene-level abundance and promoter-associated predictions | Compatible RNA coverage or promoter assays; a UMI matrix is not a base-resolution coverage track |
| Splice sites, usage and junctions | Identify exact event hypotheses missed by total-gene summaries | Event-level RNA measurements and genotype-linked or perturbation evidence for an allele-effect claim |
| Polyadenylation scores | Investigate transcript-end and isoform hypotheses | Compatible transcript-end evidence; unsigned maxima do not identify direction |
| Contact predictions | Add a candidate long-range explanation at ambiguous loci | Independent enhancer–gene or contact evidence; a contact score is not proof of a target gene |
| Coding predictions and AVI components | Keep protein-altering mechanisms visible alongside regulatory ones | Separate coding and regulatory interpretation; AVI feature contributions are not mechanism probabilities |
| Tissue comparisons | Compare liver and relevant extrahepatic candidate contexts | Matched scorer definitions, biological tissue provenance, and tissue-specific human evidence |
| Motif resources and sequence interpretation | Examine known and potentially uncharacterized regulatory sequence patterns | Verify retrievable artifacts; motif matches suggest families, not a unique bound TF or a mechanism |

Atlas publicly describes a motif resource, but the current project's point-query route does not establish access to motif calls. Verify the available download/interface before scheduling an Atlas-motif analysis. Existing independent motif resources can support a separately identified annotation. [Atlas resources](https://deepmind.google/blog/alphagenome-atlas-a-predictive-map-of-every-possible-dna-letter-change-in-the-human-genome/).

Use normalized variants and explicit genome builds. Preserve allele orientation, indel representation, transcript/event identity, tissue age and source, scorer definition, and raw effect units. Do not silently pool embryonic hepatocyte, transformed cell-line, and adult liver tracks. Static tissue outputs do not become individual MASLD states by relabeling them.

Reuse archived queries before requesting more. Track unserved variants and unannotated posterior mass. Use the model on explicitly constructed sequences for supported indels and haplotypes; archive those separately from Atlas scores because the readouts and scales can differ. Possible substitutions are hypothetical sequences, not new patients or newly observed disease alleles. [Scoring definitions](https://www.alphagenomedocs.com/variant_scoring.html), [haplotype analysis](https://www.alphagenomedocs.com/colabs/haplotype_workaround.html).

## 5. Human data and the additional information needed

The counts below describe registered sources, not an additive participant total. Reconcile donor overlap, source families, phenotype definitions, assay pairing, and usable endpoint-specific denominators before fitting. Existing accessible outcomes are development data unless a specific external evaluation has genuinely been preserved.

| Source | Role in the ambitious program | Boundary |
|---|---|---|
| GSE281364 MPRA | Local allele effects and treatment interactions in HepG2 and LX-2 contexts | Four experimental replicates per deposited context; the established 1,033-element model comparison uses two HepG2 contexts. Audit the other contexts before extending it. No patient-level inference. |
| GSE296875, 39 donors | Same-nucleus RNA–ATAC development for sequence/context learning; recorded steatosis and fibrosis | Histology is recorded for 38 and 37 donors respectively; includes minors and lacks MASLD/MASH adjudication and inflammation components. It cannot supply those missing diagnoses. |
| GSE244832, 18 donors | Human liver RNA/ATAC transport and allele-level accessibility where supported | Different aliquots, not cell pairs. Do not assume a resolved cross-assay participant join or cohort independence. |
| GSE281367, 12 donors | Adult liver accessibility context and allelic development | ATAC only; case-label semantics and detailed pathology are limited. Shares a publication family with the MPRA and is not its independent confirmation. |
| GSE267145 paired cohort, 99 participants | RNA–H3K27ac learning, differential-profile objectives, and separately registered histology analyses | Honor coordinate and source-specific assay restrictions. Existing Corgi results remain a serious comparator. |
| Resource single-cell and bulk cohorts | Cell-state vocabulary, fixed program projection, donor-aware disease remodeling and phenotype context | No invented genotype joins; counts, coverage, and cell composition remain distinct quantities. |
| Currin liver caQTL and compatible eQTL/sQTL sources | Native allele-effect evaluation or development according to exposure | Existing inspected sources cannot simply be renamed untouched tests; molecular direction must refer to the same allele and event. |
| GSE289173 registered eQTL/interaction-eQTL source | Potential protected external evaluation | Usable endpoint schema, exposure, and eligibility remain unresolved. Do not open protected outcomes for feature or model selection. |
| Protein and spatial datasets | Test observability, localization, and cross-assay associations of supported candidate genes/programs | They are not variant-effect labels without an authoritative genotype or perturbation relationship. |

Source details: [GSE281364](../../Analysis/MASLD_Model_Benchmark/config/datasets/gse281364.toml), [GSE296875](../../Analysis/MASLD_Model_Benchmark/config/datasets/gse296875.toml), [GSE244832/GSE281367 pairing rules](../../Analysis/MASLD_Model_Benchmark/OVERALL_PLAN.md), [GSE267145](../../Analysis/MASLD_Model_Benchmark/config/datasets/gse267145_znf469_human_liver.toml), and [protected GSE289173 registry](../../Analysis/MASLD_Model_Benchmark/config/datasets/gse289173.toml).

The most valuable expansion is **public genotype-linked, donor-resolved human liver molecular data with adequate context variation**. Seek existing phased genotype information and RNA/accessibility or exact event measurements where public access and reuse permit. For reporter-to-endogenous transfer, prioritize overlap in variants, molecular endpoints and relevant contexts. Newly genotyping the existing donors is outside this revision's execution path.

For aspect-specific biology, seek separately recorded steatosis/metabolic measures, inflammatory components, and fibrosis. Metabolism is broader than steatosis, and an in vitro exposure is not a fibrosis grade. A recorded NAS sum cannot substitute for its components. Better phenotype data expand the question; current indeterminacy is not a permanent exclusion.

Rare-variant association and donor-haplotype validation remain conditional on existing public participant-level genotypes, suitable phenotypes and sufficient independent observations. If those inputs cannot be obtained publicly, retain the data requirement and the prediction-only result. Neither a new cohort nor new experiments become dependencies of this paper.

### 5.1 Public-data expansion and protected evaluation

The default is public-source expansion, with metadata and existing deposits checked before a download proposal. Choose training and evaluation roles before examining new molecular-effect outcomes. A public cohort can be useful without being independent of every model's pretraining.

| Priority / source | Intended role and concrete information sought | Eligibility and completion decision |
|---|---|---|
| P1. GSE277774 / Wenz liver caQTL | Candidate external accessibility evaluation: 189 profiled liver samples, 175 used for caQTL mapping, 14,076 reported caQTLs. Use existing Currin data for initial development. | An evaluator checks study/donor overlap, genome build, allele orientation, effect/SE fields, tested-region universe and model exposure before admitting outcomes. This source is reserved from tuning by this plan; it is not declared globally untouched. |
| P2. Registered GSE289173 | Conditional external cell-type eQTL/interaction-eQTL evaluation; registry expects 48 donors and jointly held-back outcomes. | Preserve the existing evaluator-only outcome access, annotation crosswalk and checkpoint restrictions. Disease/control classification is already ineligible as a protected endpoint. Do not inspect its outcome tables or linked archive during training-data selection. |
| P3. Existing Currin and compatible public QTL releases | Broaden measured variant–region coverage beyond selected leads; acquire effect estimates with uncertainty and clearly defined targets. | Reconcile source overlap, LD, effect-allele conventions and ascertainment. Distinct eQTL/caQTL/sQTL labels stay separate. Summary statistics do not provide empirical donor-subsampling curves. |
| P4. GSE281367 / GSE244832 existing reads | Recover additional observed allelic sites across independent loci, using the complete saturation universe to locate relevant regions but not to supply labels. | Reuse existing alignments and barcode joins; prespecify genotype support, mapping sensitivity and depth criteria. Public genotype availability is an admission question, not an assumed input. Keep all cells/aliquots of each donor together. |
| P5. Public reporter, endogenous perturbation and splicing studies | Add replicated effects and matched contexts; prioritize overlap with endogenous targets and assays over library size alone. | Search accession/metadata and source methods first. Record construct, pairing, assay units and foundation-model exposure before choosing a development role. No independent status is inferred from a different accession alone. |
| P6. Additional public paired liver molecular cohorts | Increase independent donor and context information for T12–T15, and identify a distributable context-model training route. | Require an authoritative biological-unit join, measurable state overlap and clear reuse terms. If a second compatible assay is absent, report cross-assay recurrence separately from exact replication. |
| P7. Existing GSE105127 zonated RNA/RRBS substrate | Optional context extension: 19 participants with repeated zones and neighboring-section assays; RNA and CpG preprocessing already exist. | Inspect the existing registry and prepared outputs before considering a fit. Use methylation-native count/fraction models, preserve repeated-zone structure, and do not describe this sole zonation source as independent MASLD mechanism validation. |

Source records: [GSE277774](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE277774), [Wenz study and public-data statement](https://pmc.ncbi.nlm.nih.gov/articles/PMC13012615/), [GSE289173 restrictions](../../Analysis/MASLD_Model_Benchmark/config/datasets/gse289173.toml), and [GSE105127 prepared-substrate record](../../Analysis/MASLD_Model_Benchmark/config/datasets/gse105127_zonated_rna_rrbs.toml). GSE277774's public summary-statistics location is reported by its authors; endpoint-level files and donor independence were not validated in this documentation revision.

For each candidate source, T23 produces a compact admission table: endpoint, independent donor/locus counts, assay identity, pairing, overlap/exposure, public access, training use, distribution use, estimated bytes/compute and assigned role. Mark missing facts explicitly. Prefer a source that adds independent observations of the required molecular effect over one adding many correlated measurements. No controlled-access acquisition or new experimental work is scheduled by this plan.

**Public core plus research models:** train the distributable model only on sources whose terms permit the intended derivative distribution; publish eligible tables and inference recipes under their actual terms. GSE267145/GSE296875 derivative weights are releasable under the [adopted owner decision](../decisions/2026-09-07-derivative-weight-release.md); source matrices and per-participant values are not rehosted. GSE267145 sequence extraction is permitted under [counted-interval semantics with a 1-bp residual](../decisions/2026-09-15-gse267145-sequence-use.md); this does not establish called-peak boundaries or global dataset admission. AlphaGenome derivatives retain their noncommercial terms.

## 6. Define the quantities the model must learn

Let `f_m(s, c, b)` predict molecular readout `m` from sequence `s`, observed or training-defined context `c`, and assay background `b` (reporter construct or endogenous locus). Each assay has its own native output scale and fixed transformation `g_m`, specified before evaluation.

For a variant `v`, define its predicted allele effect on that scale:

`D_m(v, c, b) = g_m(f_m(ALT_v, c, b)) - g_m(f_m(REF_v, c, b))`.

The state-dependent effect is:

`J_m(v; c1, c0, b) = D_m(v, c1, b) - D_m(v, c0, b)`.

For two variants on a specified haplotype, let `F00`, `F10`, `F01`, and `F11` denote predictions on the same readout scale, window, context, and sequence background. Define:

`I_m(v1, v2, c) = F11 - F10 - F01 + F00`.

The haplotype-by-state contrast is `I_m(v1, v2, c1) - I_m(v1, v2, c0)`. These are scale-dependent molecular interaction quantities. They are not population-level disease epistasis, and a log-scale interaction is not interchangeable with one on the original count scale.

Keep five endpoints separate:

| Endpoint | Primary evaluation unit | What success establishes |
|---|---|---|
| Average regulatory profile | Held genomic regions, with donor/source stratification | Sequence-to-profile prediction |
| Between-donor molecular differences | Biological donors crossed with held regions | Useful context dependence within the measured assay |
| Signed allele effect | Independent loci; donor-by-locus units where individual measurements exist | Variant-effect prediction |
| Allele-by-state interaction | Replicate/donor-level condition contrasts grouped by locus | Prediction of measured context-dependent effects |
| Haplotype interaction | Measured combinations grouped by locus and biological replicate | Joint effects beyond the prespecified additive baseline |

A model can succeed on one endpoint and fail on another. Predicting broad chromatin profiles cannot stand in for predicting allele effects. Cross-sectional state association cannot establish an intervention effect or longitudinal progression.

Holding the measured RNA context fixed while changing sequence estimates a conditional molecular prediction. It does not simulate downstream feedback that would change the RNA context after a perturbation. Make that distinction explicit in every model-derived counterfactual.

## 7. Aim 1: construct the human regulatory response map

**Question:** Which molecular explanations of inherited MASLD signals are supported, in which measured contexts, and which remain distinguishable alternatives?

Use the adopted genetic sources with separate direct-MASLD/liver-fat and liver-enzyme universes. Keep expression-colocalized loci separate from the broader fine-mapped set to avoid selecting an expression mechanism by construction. Propagate existing signal-specific variant uncertainty; never multiply a genetic posterior by an Atlas score and call the result a probability of a true mechanism.

Integrate Atlas effects and saturation profiles with measured adult liver regions, allele effects, transcript events, and independent gene links. The scientific product is a locus-by-context-by-readout map with effect direction, units, uncertainty, and missingness. Include ambiguous genes, opposing allele directions, unserved indels, and untestable cell types.

At each locus ask whether evidence distinguishes a coding effect, local regulatory effect, RNA-processing effect, candidate long-range relationship, or an unresolved alternative. Test incremental within-locus prediction against genetic posterior rank and measured regulatory context. Report independent loci rather than counting the same LD signal repeatedly.

This map also defines representative training and validation strata: well-covered and poorly covered loci, weak and strong measured effects, concordant and discordant models, coding and noncoding candidates, and different assay contexts. Atlas confidence must not be the sole entry criterion; retain outcome-independent controls and a random eligible sample so evaluation is not restricted to favorable predictions.

**Deliverable:** a complete research layer attached to the Gene Catalog, with no automatic change to the existing evidence classes or gene membership. Static predictions, adapted predictions, measured associations, and variant-level experimental evidence remain distinguishable.

## 8. Aim 2: learn a MASLD sequence-to-function model across cellular contexts

**Question:** Can measured liver context improve molecular allele-effect prediction at unseen loci and in independent human donors?

AlphaGenome is the principal candidate backbone. Its value must be demonstrated against strong alternatives, including the existing HyenaDNA model. The intended contribution is a shared representation of allele effects that transfers across contexts and assays, with assay-specific predictions. Merely adding a liver head or naming a model after MASLD is insufficient.

### Proposed architecture

1. **Sequence representation at two scales.** Encode the actual reporter sequence and relevant construct information for the reporter task. Any required padding or shortened-encoder adaptation must respect the implementation's supported inputs and remain consistent across alleles. Encode the native genomic sequence, variant position, target coordinates, and supported long-range context for endogenous tasks. Compare shared and separate encoders. Do not attach a 1-Mb endogenous sequence to a reporter label and describe all of that sequence as experimentally assayed.
2. **Explicit context input.** Condition on measured lineage and molecular state, with assay, study, and treatment metadata handled explicitly. Start with transparent donor/lineage RNA summaries and compare a learned context encoder. Include missing-modality masks. Cohort and donor identifiers are nuisance information, not substitutes for transferable biology.
3. **Adult liver adaptation.** Learn new output heads or adapters for measured hepatocyte, stellate/fibroblast, cholangiocyte, endothelial, and immune contexts when each has sufficient independent donors and relevant labels. State which cell types are actually supported. A trained adult stellate head would be a new capability; an embryonic hepatocyte track renamed as adult liver would not.
4. **Shared allele representation with native assay heads.** Use matched REF/ALT representations to predict local activity, chromatin, RNA, or an exact transcript event. Fit each output in its own units. Test whether sharing information improves held-out effects relative to separately trained assay specialists. Represent absent assays as missing, never as zero effect.
5. **Separate profile and differential objectives.** Train on average profiles, between-donor residuals, allelic differences, and treatment differences only where each target exists. Balance objectives using development data so accurate shared profiles cannot conceal failure on donor or allele differences. This directly addresses the limitation exposed by the Corgi run.

A useful candidate decomposition on each assay's declared prediction scale is:

`molecular profile = sequence contribution + context contribution + sequence–context interaction`.

The sequence contribution preserves useful native prediction; the context-only contribution supplies a strong RNA/lineage baseline; the interaction is the additional sequence-dependent context response. Compare all three components against their additive-only and single-input ablations, with regularization toward the simpler model where labels are sparse. Center components using training data and declare the constraints used: without constraints the decomposition need not be unique. Its terms do not automatically identify separate biological mechanisms. REF/ALT subtraction at fixed context cancels a purely context-only term; a model must predict that contrast before its context sensitivity becomes evidence of allele-by-context learning.

### Match the learning objective to the assay

| Training target | Appropriate supervision | Important restriction |
|---|---|---|
| MPRA | Assay-native RNA/DNA activity estimation with replicate uncertainty; signed allele and treatment contrasts | Preserve the actual construct, normalization, pairing, and count-generating design. Do not treat estimated activities as noiseless labels. |
| ATAC | Native count/profile objectives with library-depth handling; measured allelic counts or caQTL effects for allele learning | Accessibility reconstruction and allele-effect prediction are different tasks. Correct allele-specific mapping bias before using imbalance labels. |
| Single-cell/pseudobulk RNA | Gene-count objectives with donor-aware aggregation or an explicit hierarchical count model | Gene counts do not supervise nucleotide-resolution RNA coverage or splice-junction outputs. |
| H3K27ac | The deposited assay's documented count or continuous-signal likelihood and normalization | Establish what the values represent before choosing a loss; a continuous deposited signal is not automatically a raw count. |
| Splicing/isoforms | Measured junction/event counts or uncertainty-aware event-effect estimates | Match the event, transcript, allele, and direction. A maximum absolute splice score cannot supply a signed event label. |
| Cross-assay links | Actual paired donors, nuclei, variants, or validated region–gene links | Neither shared gene symbols nor computational integration creates an observed pair. |

Loss scaling, measurement-error weighting, and any weight caps are selected inside development folds. Cell counts may improve precision within a donor, but cannot turn that donor into hundreds of independent training or evaluation subjects.

Context must be available for the declared prediction task. When predicting RNA from a measured RNA-derived state, exclude the target measurements from the context construction and fit preprocessing within training folds. Report prediction with measured molecular context separately from prediction using lineage/treatment metadata alone. Holding RNA fixed while altering genotype remains the conditional estimand in Section 6.

### Compare adaptation strategies without an artificial entry requirement

Run a prespecified comparison of a frozen AlphaGenome probe **and** parameter-efficient adaptation. A frozen-probe win is not required to try adapters: adaptation could be necessary precisely because the frozen representation is insufficient. Conversely, fine-tuning is not guaranteed to help.

Include zero-shot AlphaGenome where the endpoint is served, a calibrated native-score model, simple sequence/motif and context baselines, the existing allele-identity control, HyenaDNA, and the strongest eligible existing sequence models. Give competitive baselines adequate tuning and report training-data exposure and computational budgets. Compare selected upper-layer or full-backbone tuning if development learning curves justify the additional cost. Repeated architecture searches cannot use the external test.

Implement against the verified official research release or an explicitly validated compatible implementation. Verify REF/ALT outputs and output-coordinate alignment before interpreting a port's fine-tuning results. The [official research code](https://github.com/google-deepmind/alphagenome_research) and existing [fine-tuning implementation](https://alphagenome-pytorch.readthedocs.io/en/latest/finetuning/index.html) establish practical starting points, not evidence that our adaptation will work.

**Deliverable:** a reproducible model comparison and, if supported, a context-conditioned molecular effect model with new adult liver outputs. Retain the strongest model for each scientifically distinct endpoint; AlphaGenome need not win every task.

## 9. Aim 3: learn when reporter effects transfer to endogenous liver regulation

**Question:** Which sequence effects persist in their native genomic setting, and when do cellular state, promoter compatibility, or surrounding sequence change the prediction?

This is the most important bridge between the available MPRA and a human disease atlas. An MPRA measures activity in its particular reporter background. Episomal and integrated reporters have different limitations; neither automatically measures the variant's original chromosomal neighborhood. A local reporter effect is valuable evidence, but it is not itself an endogenous target-gene effect.

Build a development set of variants with reporter measurements and independently measured endogenous accessibility, expression, or exact RNA-processing effects. Record which combinations actually overlap in allele, locus, target, and context. Use unpaired data for separate objectives or explicitly unpaired representation learning; do not manufacture paired transfer labels.

Test three levels of transfer:

| Level | Test | Scientific gain |
|---|---|---|
| Local sequence effect | Does a reporter-trained representation improve a measured endogenous effect at held-out loci? | Establishes useful transfer beyond predicting another reporter experiment. |
| Native genomic context | Does surrounding sequence or a measured regulatory link add predictive value within the endogenous assay? | Identifies where local sequence is insufficient for prediction. |
| Cellular context | Does the measured context improve endogenous allele effects or interactions in held donors/studies? | Establishes context-dependent transfer into human liver measurements. |

Use task-specific calibration learned only from development data. Subtracting an MPRA log activity from an ATAC or RNA effect does not quantify a biological context contribution. Evaluate each prediction against its own measured endpoint, then compare whether information from the other assay improved that prediction.

Compare short local sequence, longer native sequence, and controlled flanking-sequence ablations within the same endogenous task. These are tests of model dependence. A model's sensitivity to synthetic flanks is not proof that a particular endogenous enhancer–promoter relationship exists. Contact predictions and candidate gene links generate discriminating hypotheses; measured links and variant-level endogenous evidence determine how strongly they can be stated.

Prioritize validation of disagreements: reporter-active variants predicted to be inactive endogenously, context-stable versus context-changing effects, and variants for which long-range context changes the ranking. Include concordant predictions, weak predictions, and outcome-independent controls. Avoid selecting only dramatic examples after inspecting validation outcomes.

**Deliverable:** a held-out map of reporter-to-endogenous transfer and its limits, with native assay performance, context dependence, and representative measured examples. [Speculation] Discovering reproducible rules for this transfer would provide a stronger scientific contribution than an isolated increase in MPRA correlation.

## 10. Aim 4: predict local haplotypes and genotype-by-state interactions

**Question:** When do combinations of variants produce molecular effects that differ from the sum of their individual effects, and does that difference depend on cellular state?

Atlas single-substitution scores provide a search space. Joint predictions require explicit multi-variant sequences scored by the underlying model; adding Atlas single-variant scores is the additive comparator, not a measurement of interaction. The documented [haplotype approach](https://www.alphagenomedocs.com/colabs/haplotype_workaround.html) supports constructing and comparing sequences within a shared model window.

Start with phased, observed local haplotypes where possible. Keep the same genomic window, background, readout, and context when evaluating the four sequences in Section 6. If phase is unresolved, report alternative compatible haplotypes separately. A diploid genotype requires the relevant two haplotypes and an explicit rule relating their predictions to the measured assay. Do not imply that one 1-Mb prediction models genome-wide trans interactions.

Select combinations using development-independent genetic evidence, observed co-occurrence, regulatory proximity, and prespecified model disagreement. Include combinations expected to be approximately additive. Control the candidate family before examining interaction outcomes; an unrestricted combinatorial screen followed by a few extreme examples cannot establish epistasis.

Use the following progression of evidence:

1. **Predicted interaction:** the model produces a reproducible, scale-defined non-additive contrast. This nominates a hypothesis.
2. **Measured reporter interaction:** observed joint and single-variant activities support the contrast in that reporter context.
3. **Measured endogenous interaction:** genotype-linked or variant-level endogenous measurements support the contrast in a relevant liver context.
4. **Measured haplotype-by-state interaction:** the non-additive effect itself changes across adequately replicated states.

The last two are conditional public-data aims, not claims supported by the current MPRA comparison. Seek existing variant-combination measurements with appropriate molecular readouts. When public labels are absent, this arm ends with explicitly predicted interactions and a description of the missing measurement; it does not schedule new experiments. The future pooled mouse study remains separate.

**Deliverable:** a local haplotype response map containing the additive expectation, joint prediction, native-scale interaction, phase uncertainty, and validation state. A convincing small set of independently measured interactions is more informative than thousands of unsupported interaction scores.

## 11. Biological questions opened by the combined atlas and model

These questions are research objectives. Existing negative or indeterminate analyses remain valid for their original data and estimands.

| Question | What the ambitious program adds | Evidence needed for the stronger claim |
|---|---|---|
| Which adult liver cells express an inherited regulatory effect? | New measured adult cell heads and variant-by-context predictions, beyond fixed Atlas track availability | Independent cell-resolved allele effects; distinguish plausible regulatory cell of action from the cell showing disease remodeling. |
| Which effects depend on metabolic stress, inflammation, or fibrosis? | Separate recorded phenotype components and treatment-by-allele targets | Component-specific donor data or replicated treatments with explicit interpretation. A composite severity score cannot identify the component. |
| Do RNA-processing mechanisms explain otherwise unresolved loci? | Event-specific adaptation and exact junction/isoform evaluation | Measured event effects, including supported indels; E09–E10 separate unsigned Atlas scores, weak signed direction and failed event detection. |
| Do different loci converge on a context-specific regulatory mechanism? | Saturation-derived sequence hypotheses, measured chromatin/TF context, and independent locus replication | Prespecified motif/TF/gene-set families, appropriate matched controls, and replication. A motif match is not TF occupancy or causal regulation. |
| Does inherited molecular direction align with or oppose disease-associated remodeling? | Joint display of measured genetic directions, calibrated molecular predictions, and donor-aware state associations | Same allele and target, appropriate cell context, and explicit alternatives such as feedback or cell-composition changes. Opposition does not establish a therapeutic direction. |
| Why do some variants transfer between tissues or ancestry backgrounds and others do not? | Context-specific predictions on observed sequence backgrounds and evaluation across adequately represented cohorts | Separate LD/fine-mapping differences, phenotype/assay differences, and actual predictive transport. An ancestry label alone is not a molecular context. |
| Can rare variants be grouped by a reproducible molecular consequence? | Functionally informed groups defined independently of association outcomes | Participant-level genotypes and phenotypes, allele-count/power assessment, and suitable burden or variance-component tests. Scores alone cannot establish rare-variant association. |
| Do predicted effects connect to coordinated human tissue states? | Projection onto the existing 117 programs and measured spatial/protein context | Donor-aware, externally tested associations with fixed programs. Neither program membership nor spatial proximity establishes a variant-to-program causal chain. |

Convergence is not ruled out because one Atlas interface lacks motif calls or one heritability analysis is underpowered. It requires a new, identifiable test with the necessary data. Conversely, weak existing TF and cross-assay module specificity cannot be converted to discovery by changing the null or repeatedly reselecting gene sets.

For disease association, distinguish direct MASLD endpoints from liver enzymes and liver-fat measurements throughout. Component-specific and cross-population extensions need their own prespecified estimands and covariates. Do not pool unlike phenotypes to obtain an apparently larger independent sample.

## 12. Reliability, generalization, and decisive comparisons

### Measurement quality guides the next experiment

For the four-replicate MPRA contexts, evaluate all three unique balanced two-versus-two replicate partitions, preserving the actual treatment/control pairing where present. Estimate agreement for activity, allele effects, and treatment differences separately, with uncertainty grouped by locus. The three partitions reuse the same observations and are not independent experiments.

Split-half reliability is a diagnostic, not an exact universal ceiling for every metric or for a model trained on all replicates. Any attenuation correction requires stated assumptions about independent measurement error. Estimate native count/replicate uncertainty, label coverage, and learning curves before interpreting a small benchmark difference. Retain weak and uncertain measured effects in an uncertainty-aware evaluation rather than testing only significance-selected positives.

If a treatment interaction is not reliably measured, the current dataset cannot validate that endpoint. Examine the uncertainty, number of independent loci and replicates, and likely effect sizes to determine what additional information is needed. This is a reason to improve the measurements or revise the present claim, not to abandon the scientific question.

### Split at the level of the scientific claim

| Generalization claim | Required separation |
|---|---|
| New locus | Keep alleles, nearby/LD-linked variants, repeated constructs, and overlapping receptive fields in the same split; use buffers sufficient for the context length. |
| New donor | Keep every cell, aliquot, molecular profile, and repeated measurement from a donor together. |
| New study | Hold out an independent source and publication family, with compatible endpoint definitions. GSE281367 is not an independent study confirmation of GSE281364. |
| New cellular state | Hold out a defined context where the design has sufficient overlap to separate state from study and assay. Confounded state/cohort labels cannot establish state transport. |
| New donor and locus | Cross donor and genomic holdouts; succeeding on either marginal split is insufficient. |
| New haplotype | Keep related combinations at a locus together for the strict test; report interpolation among seen combinations separately. |

Fit normalization, context encoders, feature selection, assay calibration, and hyperparameters inside the permitted training/development partition. Count pretraining exposure separately from adaptation and test exposure. GTEx use in pretraining does not mean every GTEx variant-effect label was trained on; absence from our fine-tuning set does not prove absence from foundation-model pretraining.

Existing inspected public benchmarks remain useful development or exposed evaluations. Preserve a genuinely independent test where possible. GSE289173 is a candidate only after its endpoint schema, access rules, and prior exposure are resolved without using protected effect estimates to choose the model. An independent evaluator should execute the final locked comparison when protected data permit it.

### Prespecify what a model must beat

| Endpoint | Primary comparison | Supporting measures |
|---|---|---|
| Signed endogenous allele effect | Paired improvement over the strongest eligible baseline on held-out signed-effect correlation, with the aggregation rule fixed | Native-scale prediction error, slope/calibration, sign accuracy conditional on measurement precision, coverage |
| Allele-by-state interaction | Native-scale interaction prediction error versus a context-independent model and an additive state baseline | Interaction correlation where estimable, predictive likelihood, calibrated direction and uncertainty |
| Haplotype interaction | Error for measured joint effects versus the sum of measured/predicted singles, specifying each comparator | Error for the interaction contrast, uncertainty coverage, additive-control behavior |
| Reporter-to-endogenous transfer | Endogenous endpoint gain from adding reporter supervision to an otherwise matched model | Matched-variant sensitivity, assay specialists, context ablations, locus/study transport |
| Between-donor profiles | Improvement in held-donor differential prediction over simple RNA/context baselines | Mean-profile accuracy reported separately, donor and cohort robustness |

For confirmatory model claims, name one primary contrast per endpoint, define the hypothesis family in advance, and use Holm correction at family-wise 0.05. Secondary variant, event, or mechanism screens use Benjamini–Hochberg correction within explicitly named families; report nominal estimates and full tested universes. Genetic posterior probabilities are not p-values and do not enter these corrections as if they were.

Use confidence intervals that respect the crossed structure of the data: donors and locus blocks for donor-linked effects, experimental replicate and locus structure for reporter data, and independent studies where available. Choose the specific resampling or hierarchical model before the final comparison. With very few studies, report study-specific effects and acknowledge that a study bootstrap is poorly supported.

Use the existing comparison seeds `1103, 2909, 4721, 6673, 8111` where applicable. Seeds measure optimization variability; they do not increase biological n. Keep folds paired across models. Report all planned seeds and failed fits, and use nested development selection rather than selecting the best seed.

Define a scientifically useful improvement before opening the final test, using measurement uncertainty, baseline performance, and attainable sample size. Do not declare an arbitrarily tiny positive difference sufficient or invent a guaranteed target correlation. Learning curves should identify whether the next investment should be additional donors, loci, replicated conditions, or model capacity.

### Report where predictions are unsupported

Fit uncertainty calibration on separate development data and check empirical coverage on the held-out population. Variation across seeds alone is not a calibrated prediction interval. Report the applicable cell context, assay, genomic coverage, measurement precision, and degree of distribution shift with every prediction.

Distinguish a calibrated near-zero effect from an uncertain effect, an unserved assay, and an unsupported context. For a new lineage, phase, or exposure outside the training support, return that limitation explicitly. Do not turn an out-of-distribution flag into an invented probability that the biological prediction is wrong.

## 13. What would make the contribution distinctive

The intended result is **a human liver regulatory response atlas that predicts and tests how inherited sequence effects depend on cellular and genomic context**. The following claims require different evidence and should be assessed independently.

| Proposed advance | Evidence that would substantiate it | Insufficient substitute |
|---|---|---|
| MASLD-relevant sequence-to-function adaptation | Reliable improvement on independent human liver molecular allele effects, with competitive baselines and declared context support | Lower training loss or accurate shared profiles |
| Transfer between experimental backgrounds | Reporter supervision improves independently measured endogenous effects, with useful rules for when transfer fails | MPRA performance alone or agreement between two prediction tools |
| Context-dependent regulatory genetics | Prediction of replicated allele-by-state effects in a context not used for model selection | Disease classification from RNA or different average expression between cases and controls |
| Local haplotype response prediction | Measured joint effects outperforming the additive expectation, with a valid interaction contrast | Nonlinearity in model outputs alone |
| Mechanistic resolution at disease loci | Independent evidence distinguishes competing molecular explanations at held-out loci | Many annotations attached to the same LD signal |

No claim of being the first such model is made here. Before manuscript novelty language is finalized, check prior MASLD applications as well as the source methods, including AlphaGenome adaptations, scooby, Corgi and the candidates in section 2.1. A novel MASLD adaptation, a useful new disease interpretation, or a demonstrated predictive improvement can each justify a contribution with the corresponding evidence. Report a methodological invention only when one has actually been made.

### Scientific artifacts for the atlas paper

- **Human regulatory response map:** variant/locus, allele orientation, observed haplotype where available, target, molecular assay, measured context, signed effect, uncertainty, and support status. Keep measurements and model predictions in separate fields.
- **Adapted model and evaluation package:** reproducible inference and training recipes, permitted weights/adapters, independent splits, native assay evaluations, learning curves, and exposure documentation. Include a competitive baseline even if AlphaGenome becomes the selected model.
- **Reporter-to-endogenous benchmark:** actual overlap tables, held-out effects, transfer results, and failures. This is reusable by the broader AI/ML and computational biology community.
- **Local haplotype map:** measured versus predicted joint effects, additive comparators, phase limitations, and supported state interactions.
- **Gene Catalog integration:** a research layer that answers which effect was measured, which was predicted, what context it applies to, and which competing explanation remains unresolved.

An inference request should specify a single variant or local haplotype, a reference assembly/sequence background, and the available molecular context. Its output should be a molecular prediction in defined units with uncertainty and applicability, not a single unsupported MASLD-risk score. Research users should be able to compare contexts and inspect the evidence behind a prediction.

### Result-led figure story

Use the existing seven-figure manuscript structure as the starting point. The following are proposed panel roles, not changes to adopted figures:

| Figure role | Proposed content |
|---|---|
| Inherited variation, currently Figure 2 | Full locus universe with molecular alternatives, independent within-locus evaluation, and representative supported mechanisms |
| Human cellular context, currently Figures 4–5 | Measured adult liver contexts and held-donor/state effect comparisons, with donor counts and assay support visible |
| Public resource, currently Figure 6 | Variant/haplotype-to-context navigation in the Gene Catalog, including unresolved and untestable cases |
| Predictive framework, currently Figure 7 | Competitive external allele-effect results, reporter-to-endogenous transfer, context ablations, and haplotype results where measured |

The key panels should show signed effects, independent units, uncertainty, and held-out predictions. A model schematic supports those results. Complete assay, lineage, and negative-result tables belong alongside the headline panels so the community can see the limits of transfer.

Select examples by a recorded rule that includes representativeness and assay support. Existing public variant-level validation may support this Resource; new experimental studies and the future pooled mouse screen are not dependencies. Keep known source-owned mechanisms identified as such.

## 14. Work sequence, resources, and completion criteria

Work proceeds according to information dependencies. Preserve the current paper structure and let completed evidence determine additions. The table below retains the original work-package vocabulary for interpreting executed A–G deposits; the six forward campaigns in section 19 govern new assignments.

| Work package | Concrete output | Completion criterion | Main dependency |
|---|---|---|---|
| A. Establish the training/evaluation substrate | Rights and exposure table, exact assay/participant joins, replicate uncertainty, native endpoint definitions | Every training target and test unit is identifiable; missing labels have an explicit acquisition route | Existing source records and focused metadata inspection |
| B. Build the Atlas research layer | Selected locus scores, saturation maps, all available relevant modalities, measured-context overlap | Effect orientation, units, coverage, and model/source version are recoverable; missing modalities are explicit | A and lawful access to the relevant artifacts |
| C. Establish adaptation results | Frozen probes and adapter pilots, simple/strong sequence baselines, profile-versus-differential evaluation | Paired development comparisons and learning curves support a model/data investment decision | A; released weights and suitable compute access |
| D. Train and test the context model | Native adult liver heads, shared-versus-specialist comparison, protected external test | The prespecified human molecular effect comparisons are reported, including failures | C and genotype/context-linked labels |
| E. Test reporter-to-endogenous transfer | Measured bridge dataset and independent endogenous predictions | Transfer evaluated at unseen loci/studies; incompatible assay effects never pooled as one target | C–D and measured overlapping variants |
| F. Extend to haplotypes and discriminating validation | Phased sequence comparisons, additive controls, existing public joint-effect measurements where available | Joint effects and context interactions evaluated in their own native units; prediction-only status retained when measurements are absent | Reliable single-variant tasks, public phase, and combination measurements |
| G. Integrate the paper and public research product | Result-led panels, complete tables, model artifacts and Gene Catalog research layer | Claims match independent evidence; results and manuscript artifacts receive the ordinary scientific review | Completed evidence from B–F, with unsupported aims reported honestly |

A–B and the expanded tasks below can use the existing local-weight research route. Frozen and adapter pilots can be compared within C without requiring one to win first. Haplotype sequence preparation can proceed before a validated interaction model exists, but its outputs remain hypotheses until F supplies measurements. The initial A–F analyses measured feasibility and specific comparisons; they did not complete all of the broader scientific aims in this table.

### A deliberate data-development program

After A, prepare the public-source proposal in section 5.1. Prioritize independent genotype-linked adult liver measurements, matched reporter/endogenous observations, replicated context contrasts and measured local combinations. Use pilot uncertainty and independent-locus counts to estimate precision. Do not invent a universal required sample size before observing the variance and design.

Evaluation in newly admitted public sources should include fixed eligible populations, weak/additive controls, independent biological units and the readout that distinguishes competing hypotheses. Preserve complete eligible denominators alongside model-disagreement subsets. Reporter-only measurements validate the reporter task. New experiments can be described as the next discriminating measurement in a database record, but are outside the execution dependencies.

### Compute should follow measured learning curves

Reuse the existing benchmark substrate and targeted Atlas caches. Start with metadata, labels, and small output subsets rather than attempting to download the complete Atlas. Measure memory, runtime, I/O, and accuracy on representative supported sequence lengths before sizing a larger run. The local inference pilot measured approximately 17.6 GiB and 1.60 seconds per warm variant score at 1,048,576 bp on L40S; this is not a measurement of backward-pass memory, adapter-training throughput, or a useful training batch size. Those remain [Unverified].

Cache frozen embeddings only for frozen encoders. Trainable-backbone outputs must be recomputed consistently during adaptation. Compare the gain from longer context against both the short-context baseline and its computational cost; different sequence lengths also require appropriate split buffers.

Training and substantive analytical workloads belong on compute nodes. A later execution proposal should state job count, CPUs/GPUs, memory, wall time, and aggregate allocation before submission. Follow the project's partition and QOS rules, including the per-user memory cap; no arrays, new dependencies, or jobs are authorized by this document-writing request. Record seeds and the actual software/model versions with each result.

## 15. Model and Atlas usage terms

Treat the hosted Atlas/API and the released model as distinct sources with distinct terms. The official [AlphaGenome terms](https://github.com/google-deepmind/alphagenome#terms) restrict training other machine-learning models on hosted outputs unless an explicit permission applies to the particular artifact. The [released model card and license](https://huggingface.co/google/alphagenome-all-folds) govern the local weights, fine-tuned derivatives, and output-based model derivatives under noncommercial conditions.

The default design is therefore:

| Source | Planned use |
|---|---|
| Hosted Atlas/API outputs | Interpretation, hypothesis generation, and evaluation within the applicable terms; exclude from training inputs, targets, or learned calibration by default |
| Released AlphaGenome weights | Local inference and adaptation under the accepted model license |
| Outputs generated from those released weights | Eligible local feature/teacher route subject to the derivative-model terms; record this origin distinctly from hosted outputs |
| Human experimental measurements | Training and evaluation only with the dataset's applicable access, consent, and reuse conditions |
| Public motif or other Atlas-associated artifacts | Check the exact artifact and license before selecting a training use; do not infer permission from public availability |

Verify the current terms at acquisition and distribution. Noncommercial access is not an unrestricted open-source license, and an independently trained student can still be a covered model derivative if it was trained using restricted model outputs. Distribute checkpoints, adapters, derived features, and datasets only as the applicable licenses permit. These distinctions determine the lawful technical route; they do not imply that building a MASLD adaptation is scientifically infeasible.

## 16. Decisions to preserve when implementation begins

The scientific target is an experimentally assessable regulatory response atlas for the MASLD paper. Current zero-shot failures, missing tracks, and limited labels identify adaptation and data-development needs. They neither demonstrate that the proposed model will succeed nor establish that those questions are impossible.

The first implementation specification should fix the allele and state estimands, training rights, independent data units, endpoint-native losses, protected comparisons, and resource estimates described above. It should carry forward existing negative results and include both strong sequence baselines and independently justified adaptation. It should not quietly revert to a frozen-embedding MPRA contest as the whole objective.

The most consequential outcome would be evidence that we can predict **which inherited molecular effects change with cellular context, which reporter findings persist in native human liver, and which local haplotypes violate the additive expectation**. Those are the claims to design the work around. Their inclusion and prominence in the manuscript should follow the resulting evidence, with timing handled as a separate decision.

## 17. Expanded experimental agenda

### 17.1 Six scientific goals

The next level is a model that adds experimentally assessable information about human liver regulation. The ambition includes new representations, new molecular outputs and better data, as well as better tuning of existing heads.

| Goal | Scientific advance to pursue | Decisive evidence | Tasks |
|---|---|---|---|
| G1. Improve endogenous allele prediction | Predict the direction and magnitude of molecular allele effects more accurately than the strongest eligible current model | Paired improvement on identical loci, followed by an independent compatible source; both correlation and native-scale error reported | T01–T09, T20–T21 |
| G2. Learn cellular-state dependence | Predict which allele effects change across treatments, adult liver lineages and measured donor states | Improvement on measured allele-by-context contrasts beyond sequence-only and additive-context models | T12–T15, T18, T20 |
| G3. Transfer between assays | Learn when local reporter effects inform endogenous accessibility, expression or RNA processing | Endogenous prediction improves when reporter or another assay's supervision is added to an otherwise matched model | T10–T11, T16–T17 |
| G4. Resolve molecular alternatives | Distinguish local accessibility, target-gene expression, splicing and haplotype explanations at inherited liver-trait loci | Correct target/event/allele predictions and independent evidence discriminating among explanations | T16–T19, T22 |
| G5. Generalize with useful uncertainty | Make calibrated predictions for new loci, donors and studies, with explicit limits in poorly supported contexts | Locked external comparisons, interval coverage and error-versus-coverage curves; failures retained by subgroup | T18, T20–T21, T23 |
| G6. Deliver a reusable scientific resource | Connect measured and predicted molecular responses to the MASLD Gene Catalog and fixed disease programs | Reproducible per-endpoint models, complete effect tables and evidence-linked examples that stand independently of the separate follow-up study | T22–T24 |

**Planning ambitions, not observed effects or automatic success thresholds:** Campaign A's declared stretch target is +0.02 absolute Spearman over the strongest matched baseline; retain +0.03 to +0.05 as a higher aspiration where the endpoint permits it. For context interactions, seek a 10% reduction in native prediction error; for transfer, 5–10%. Seek replication across at least two compatible independent source families for the broadest transport claim and measured discrimination of molecular alternatives at several independent loci. Feasibility is [Unverified]. These ambitions guide development and precision calculations; final minimum useful effects must be fixed in native units before protected evaluation. A smaller precise improvement and an unresolved large point estimate have different interpretations.

### 17.2 Task map

The table describes scientific dependencies, not a queue of authorized jobs. A task is complete when its comparison and output exist, including an informative negative result. A positive biological claim has the additional success requirements stated below.

| Task | Question | Starting substrate | Dependency |
|---|---|---|---|
| T01 | Are the labels and evaluation units sufficient for the claim? | Existing reporter, caQTL, allelic and donor-assay records | Start now with inspection |
| T02 | What is the strongest fair comparator for each endpoint? | Existing local scores and model predictions | T01 |
| T03 | How much sequence and which genomic background matter? | Local sequences and existing labels | T01–T02 |
| T04 | Which positions, resolutions and layers carry the effect? | Frozen local representations | T03 |
| T05 | Can better heads use information already present? | Local scores and frozen representations | T02–T04 |
| T06 | Which small trainable updates improve allele effects? | Local backbone and assay-native labels | T01–T05; may start alongside T05 |
| T07 | Does broader unfreezing add value? | Frozen and adapter comparisons | Stable training pilot; T06 results inform scale |
| T08 | Which loss preserves the biological quantity? | Assay-native measurements and uncertainties | T01; developed alongside T05–T07 |
| T09 | Which training schedule and regularization improve transport? | Selected architectures and native losses | T05–T08 |
| T10 | Which assays benefit from shared learning? | Compatible assay targets with missingness recorded | T01, T05–T09 |
| T11 | Which reporter features transfer to native liver regulation? | Measured reporter/endogenous overlaps | T10 and an adequate matched label set |
| T12 | Can we learn new adult liver cell-type outputs? | Donor-aware RNA/ATAC profiles and eligible allele labels | T01; richer allele coverage for sparse lineages |
| T13 | Which cellular-context inputs are useful and transferable? | Donor RNA, fixed programs, treatment and lineage metadata | T12; can also start on the paired bulk cohort |
| T14 | Can we predict treatment-dependent allele effects? | LX-2 first; other contexts after measurement development | T08, T13 and reliable interaction labels |
| T15 | Can sequence plus context predict donor differences? | Paired GSE267145 and eligible adult multiome donors | T12–T13 |
| T16 | Can the model identify the affected gene and expression effect? | Measured gene-level effects and independent region–gene links | T02, T10, T12 |
| T17 | Can adaptation recover event-specific RNA-processing effects? | Junction/event measurements, sQTL and supported indels | T01–T08; event identity and exposure resolved |
| T18 | Does actual donor sequence improve prediction? | Genotype-linked donor assays | Genotypes, phase and assay joins |
| T19 | Can measured haplotype interactions be predicted? | Phased combinations with joint-effect measurements | T18 and adequate combination labels |
| T20 | Which gains survive new loci, donors, states and studies? | Protected evaluation populations | Begins with T01; final scoring after model selection |
| T21 | Are uncertainty and applicability useful? | Calibration data and locked predictions | T20 split design; sufficient calibration units |
| T22 | Which molecular mechanisms are independently supported? | Complete predictions plus orthogonal measurements | Relevant completed molecular tasks |
| T23 | Which additional data would improve the model most? | Learning curves, label precision and coverage | Starts with T01; revisited after each campaign |
| T24 | What enters the reusable Resource? | Completed comparisons, complete tables and supported examples | Scientific review of the relevant tasks |

### 17.3 How to read the variables

All numerical ranges below are **proposed development search ranges**, not established best settings, package defaults or a promise that every option is implemented. A supported input length does not establish that a particular training loop fits that length in memory. New architectural hooks, low-precision variants and alternative optimizers need implementation checks against the pinned local code before they enter a comparison.

Separate three types of experiment:

- **Use existing information better:** pooling, output heads, local-score calibration, regularization and loss.
- **Learn new information from existing measurements:** adapters, unfreezing, joint assays and context conditioning.
- **Acquire information absent from the current measurements:** independent genotypes, rare-lineage allele effects, reliable additional states, native target/event labels and measured haplotypes.

All three belong in the program through the public-data route in section 5.1. Missing information may remain unavailable publicly; that produces a documented capability limit rather than a new experimental dependency. More model capacity cannot be assumed to replace biological replication.

### T01. Endpoint and label development

**Goal:** produce the exact training and evaluation populations needed for G1–G5, including weak effects, missing assays and genuine biological replication.

**Supporting results:** E04–E05, E09–E10, E20, E25–E27.

**Unresolved question:** Which target definitions and ascertainment rules support the intended claim?

**Data:** Reporter replicates, Currin effect/SE tables, donor allelic counts, paired profiles and eligible events.

**Evaluation unit:** Locus for genetic prediction; donor for human context; experimental replicate for reporter measurement.

**Success criterion:** Every target has recoverable identity, native units, uncertainty, independent groups and missingness; no unmeasured target becomes a zero label.

**Dependencies:** Existing records; T20 split design developed alongside this task.

**Consequence for the paper:** Define exactly what each Figure 7 comparison measures and which database fields are measured.

**Inputs:** reporter replicate measurements; Currin lead and nominal tests; donor-level allelic counts; RNA/ATAC/H3K27ac; event-level RNA measurements where available.

**Variables to try:**

| Axis | Candidate comparisons |
|---|---|
| Label universe | Existing significance-selected caQTL leads; outcome-independent eligible variant–peak tests; precision strata retained for evaluation rather than used to remove difficult examples |
| Label estimation | Existing point estimate; replicate-aware shrinkage estimate; point estimate plus sampling variance; raw-count model when the necessary counts and design exist |
| Variant–target identity | One predeclared target per variant versus all explicitly tested targets with clustering by locus; keep conflicting target directions visible |
| Replication | Existing replicate means versus models retaining experimental replicate or donor variation; paired versus unpaired condition models according to the actual design |
| Mapping and genotype | Existing read-based heterozygote calls versus independently supported genotypes and allele-bias-controlled reads; explicitly defined sensitivity populations |
| Selection | Strong and weak measured effects, random eligible controls, model-disagreement examples, and all eligible loci; keep a representative evaluation population |
| Estimand alignment | Separate association, isolated variant effect and allelic-count targets; conditional LD-aware supervision only with matched LD, scaling and covariance; compare against independent native measurements (22.3.1) |

**Decisive comparison:** determine how each label definition changes reliability, coverage, calibration and model ordering on a common population. Evaluate precision without treating correlated variants or repeated assays as independent donors.

**Output and completion:** an endpoint table with assembly, alleles, target/event, assay units, uncertainty, source, donor/replicate keys and split assignment; plus a reliability and missing-data report. Report biological n separately from loci, sites, cells and reads. A model task lacking identifiable labels receives a concrete data requirement.

### T02. Strong and comparable baselines

**Goal:** identify what each new model must improve upon and separate representation improvement from favorable comparator choices.

**Supporting results:** E03, E24–E25, E33–E38.

**Unresolved question:** Does an apparent gain survive comparison with the strongest eligible baseline?

**Data:** Existing native scores, specialist predictions, RNA/context baselines and identical endpoint populations.

**Evaluation unit:** Matched loci and donor groups; score comparisons remain paired.

**Success criterion:** Reproduce each selected baseline on the common population and report exposure/coverage asymmetries; retain constant or failed fits.

**Dependencies:** T01; reconcile local score provenance before calibration.

**Consequence for the paper:** Ensure model comparisons support an incremental capability claim.

**Variables to try:**

| Baseline family | Candidate arms and purpose |
|---|---|
| Simple sequence | Allele identity, local sequence composition, motif features and a small convolutional model trained from scratch; measure the value beyond these features |
| Existing DNA representations | HyenaDNA and Caduceus with comparably tuned heads, supported sequence lengths and identical eligible examples |
| Local AlphaGenome | Native ATAC/DNase or endpoint-appropriate score; per-track scores; fixed track average; training-only linear or nonlinear calibration; score-plus-context model |
| Assay specialists | Existing eligible ChromBPNet/Borzoi comparisons for accessibility and a suitable event specialist for RNA processing, with the exact artifact and exposure named |
| Context | Training-mean profile, lineage/treatment-only predictor, RNA-only ridge/reduced-rank model and the existing Corgi recipe on its matched target |
| Architecture control | Small task-trained model, randomly initialized head on frozen features, and selected capacity-matched alternatives |
| Adapted methods | Original source-method implementation, its proposed MASLD adaptation and a simple endpoint-native model; section 2.1 supplies candidates and modality checks |
| Updated target/context specialists | ENCODE-rE2G and scE2G for compatible links; GET for its observed-chromatin RNA task; CREsted-style compact lineage models; Pangolin/MTSplice/APARENT2 only on matching events (R03–R04, R12–R13, R22–R24) |

**Decisive comparison:** paired performance on identical variants, targets, windows where feasible, contexts and folds. Report both equal-search-budget comparisons and each baseline's strongest development-selected version. Explain genuine differences in pretraining, track averaging and receptive field instead of calling the entire comparison perfectly controlled.

**Output and completion:** a baseline table per endpoint with correlation, native-scale error, calibration, coverage, exposure and compute. Local AlphaGenome's 0.7441 and HyenaDNA's 0.1737 must never be compared as if measured on the same variant set.

### T03. Sequence length and genomic background

**Goal:** discover the sequence context needed for local effects, distal regulation and native-versus-reporter differences.

**Supporting results:** E28, E32–E35; the 198-bp training-lead distance in section 3.1.

**Unresolved question:** Does extra context add information, or mainly change eligibility and pooling?

**Data:** Matched reference/alternate windows and exact reporter or endogenous target coordinates.

**Evaluation unit:** Held genomic groups, buffered for the longest compared input; independent donors when relevant.

**Success criterion:** A useful length improves measured-effect accuracy or coverage at justified cost on a common population; report the full eligible population separately.

**Dependencies:** T01–T02; validated local sequence extraction.

**Consequence for the paper:** Show which biological questions require local or distal context.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Native input length | 2,048; 16,384; 131,072; 524,288; 1,048,576 bp, using the verified local input interface |
| Reporter input | Actual validated reporter sequence; complete reconstructed construct when documented; reporter plus documented vector/promoter context; native-flank comparator. Short-encoder processing is a separate implementation, not an unsupported short input to the full trunk |
| Centre | Variant-centred; measured peak-centred; target promoter-centred; paired variant and target windows with explicit coordinates |
| Distal information | Native long window; two local windows fused after encoding; local plus pooled distal summaries; target-aware sparse selection from a longer window |
| Background | GRCh38; verified alternate haplotype background; donor sequence where genotypes are known; controlled flank ablations labelled as model-dependence tests |
| Sequence uncertainty | Strict ACGT subset as the current reference; an explicitly designed ambiguous-base/mask-aware extension evaluated separately if supported, never silent replacement |

**Decisive comparison:** matched-population short versus long predictions with the same model/head, plus performance on each length's full eligible population. Record excluded regions and their effect distribution. Compare variant-centred and target-centred models on the same biological target. Agreement with Atlas is an implementation diagnostic; measured endogenous effects decide which length is useful.

**Output and completion:** accuracy, coverage, memory and runtime by length, target distance and assay; nominate one length per primary endpoint before final evaluation. The 198-bp caQTL median favors testing a strong local model, while gene-expression and splicing tasks may justify different contexts.

### T04. Representation, layer and pooling search

**Goal:** extract allele-sensitive information without discarding spatial information or confusing a bin boundary with biological localization.

**Supporting results:** E34–E35.

**Unresolved question:** Which representation preserves the focal allele effect and its spatial location?

**Data:** Local encoder/trunk/decoder features from fixed weights; target/variant coordinates.

**Evaluation unit:** Held loci; all projections and feature selection fitted inside training folds.

**Success criterion:** A pooling/layer gain survives coordinate shifts, common-example comparison and nested selection; independent evaluation confirms any broad claim.

**Dependencies:** T03 initial length comparison; T01–T02.

**Consequence for the paper:** Explain the representation gain and its assay-specific limits.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Position | Variant bin; bins immediately to its left/right; symmetric two-bin or three-bin neighborhood; measured peak; promoter; exon/event; combined variant and target |
| Pooling | Mean, sum, maximum, concatenated bins, learned attention, target-distance-weighted pooling and a short convolution over neighboring bins |
| Window | 128, 256, 512, 1,024, 1,536 and 2,048 bp where aligned to the actual representation; compare genomic coverage explicitly |
| Layer | Selected encoder stages; last transformer output; learned mixture of selected layers; decoder output; output-track features; verify each hook in the installed model |
| Resolution | 128-bp representation first; targeted 1-bp slices for event or boundary tasks; multiresolution fusion. Avoid storing whole-genome-sized 1-bp tensors merely to pool them |
| Allele encoding | ALT−REF; REF/ALT concatenation; signed and absolute differences; per-allele nonlinear head followed by subtraction; explicit target-coordinate embeddings |
| Compression | None, centering/scaling, PCA/SVD, whitening versus no whitening, learned bottleneck; dimensions 32, 64, 128, 256, 512, 1,024 subject to training rank |
| Orientation and offsets | Forward only versus reverse-complement averaging; offsets around the 128-bp boundary; coordinate-aligned shift averaging |

**Decisive comparison:** initially cross length with pooling, then combine the best development candidates with layer/resolution changes. Fit every projection inside the training partition. Test whether a pooling gain survives re-centering and reverse-complement coordinate alignment.

**Output and completion:** a representation ablation table that explains where the predictive information resides. Carry promising variant-bin results forward as development choices; confirm their benefit on an independent test or through selection fully nested within a new evaluation design.

### T05. Frozen heads, native-score calibration and residual prediction

**Goal:** establish how much can be gained without changing the pretrained backbone, and preserve useful native predictions while learning their errors.

**Supporting results:** E25, E33–E34.

**Unresolved question:** Is useful information already present in native scores or fixed representations?

**Data:** Eligible local scores and embeddings; measured effects and sampling uncertainty.

**Evaluation unit:** Paired held-locus comparisons, with donor grouping for context targets.

**Success criterion:** Improve ranking or native error beyond the strongest matched score/head; distinguish calibration-only improvement from ranking improvement.

**Dependencies:** T02–T04; T08 loss design.

**Consequence for the paper:** Provide a competitive inexpensive model even if backbone adaptation adds no gain.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Linear head | Ridge, elastic net and reduced-rank multitask regression; ridge penalty 1e−4 through 1e7 on consistently standardized features, including an intercept-only fallback |
| Nonlinear head | One to three hidden layers; width 64, 128, 256, 512; GELU or ReLU; dropout 0, 0.1, 0.3, 0.5; a small convolution/attention head for spatial bins |
| Calibration | Identity, linear, monotone calibration and small nonlinear score-only head. Monotone calibration targets magnitude/error, not an artificial expectation of higher Spearman correlation |
| Residual model | `measured effect = calibrated local score + learned residual`; compare frozen-score plus embedding, score plus context and score plus both |
| Allele symmetry | Unconstrained paired head versus a shared per-allele function or explicitly antisymmetric difference head; swapped alleles reverse the effect and identical alleles give zero |
| Output structure | One head per assay/context; shared head with context input; shared trunk with assay-specific intercept, scale and residual |
| Feature selection | None; training-only selection; regularized bottleneck; comparable parameter counts across representation choices |
| Mean versus contrast preservation | Separate native mean-profile and native REF–ALT retention terms; compare their value for the same allele-effect task, with identical-allele and swap checks (22.3.3) |

**Decisive comparison:** each head against both the native-score/calibration baseline and the existing representation head. Use correlation and native-scale loss together. When a constant predictor is selected, keep its valid error statistic; record correlation as undefined and apply a prespecified aggregate-metric rule rather than silently dropping that fold.

**Output and completion:** selected frozen recipes, complete failed/constant-fit counts and paired increments. This task can succeed even if a simple calibrated score beats every embedding head.

### T06. Parameter-efficient backbone adaptation

**Goal:** test whether changing a small part of AlphaGenome learns endogenous allele information absent from the frozen representation.

**Supporting results:** E32–E34, E36.

**Unresolved question:** Can small backbone updates recover the molecular differences missed by fixed features?

**Data:** Local weights and custom assay-native batches from T01.

**Evaluation unit:** Held loci and, for context learning, crossed held donors.

**Success criterion:** Intended gradients and parameter updates work; adaptation improves the selected endpoint beyond its matched fixed model without concealing native-function loss.

**Dependencies:** T01–T04 and T08; T05 may run alongside the adapter pilot.

**Consequence for the paper:** Support an adaptation claim only for the measured endpoint and comparison.

Low-rank adaptation learns small matrix updates while retaining the original weights. It was developed as a general model-adaptation method; its use here requires assay-native losses and verified insertion points, not a borrowed language-model training objective. [LoRA paper](https://arxiv.org/abs/2106.09685)

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| LoRA rank | 4, 8, 16, 32, 64 |
| Update scale | Alpha/rank ratios 0.5, 1 and 2, with the implementation's scaling convention stated |
| Adapter dropout | 0, 0.05, 0.1 |
| Placement | Attention query/value; query/key/value/output; feed-forward projections; selected encoder or decoder projections where their tensor structure permits |
| Depth | Final quarter of eligible blocks, final half, all eligible blocks; compare depth at matched trainable-parameter budgets |
| Adapter family | LoRA; bottleneck residual adapters of width 32–256; learned channel scaling; bias/normalization-only updates |
| Sharing | One adapter shared across assays; assay-specific adapters; lineage-specific adapters; shared adapter plus a small context-dependent residual |
| Initialization | Identity/zero-effect initialization; small random update with identical initial predictions checked where required |
| Precision | Standard bfloat16 compute with float32 reductions; selective higher precision for small allele differences. Quantized-backbone adapters are a later implementation study, not assumed available in this JAX runtime |

**Decisive comparison:** adapter versus its frozen counterpart with the same input, head, loss, examples and tuning opportunity. Compare native-score retention and held-out allele accuracy. A frozen-probe win is not a prerequisite for trying adapters.

**Output and completion:** rank-by-placement-by-depth results, trainable parameter counts, memory and throughput, with a selected adapter or an explicit lack of gain. Before scaling, demonstrate nonzero intended gradients, unchanged frozen parameters, allele-coordinate consistency and reproducible effect tolerances on a small compute-node example.

### T07. Partial unfreezing, progressive tuning and full fine-tuning

**Goal:** determine whether useful adaptation requires changing broader sequence features rather than only low-rank updates or output heads.

**Supporting results:** E25, E32, E36.

**Unresolved question:** Is the representation itself limiting learning after head/adaptor changes?

**Data:** The same labels and inputs as T05–T06; eligible native replay data where used.

**Evaluation unit:** Identical held groups and declared training exposure.

**Success criterion:** Broader tuning improves the selected endpoint enough to justify measured cost, or yields a precise scoped limit; report retained native performance.

**Dependencies:** Stable training pilot; T06 informs scale but need not win; T08.

**Consequence for the paper:** Determine whether the paper's model needs substantial backbone modification.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Trainable components | New head only; head plus decoder; final quarter/half of transformer blocks; selected convolutional stages; encoder only; whole backbone |
| Schedule | Immediate partial unfreezing; head warm-up then unfreezing; progressively deeper unfreezing; adapter warm-start followed by broader tuning |
| Learning rates | Backbone 3e−7 to 3e−5; head/adapter 3e−5 to 3e−3; test ratios rather than applying one rate to every component |
| Layer-wise rate decay | 0.7, 0.85, 1.0 across successively earlier trainable blocks |
| Retention | No retention term; regularization toward pretrained parameters; local teacher-output consistency on training regions; limited replay of eligible native-profile training examples |
| Task order | Native profiles then allele effects; allele effects directly; joint profile/allele training; context interactions introduced after or jointly with allele learning |
| Initialization comparator | Same architecture frozen/pretrained; selected smaller model trained from scratch; task-head initialization controlled across runs |
| Specificity curriculum | Consensus regions followed by training-defined lineage-specific refinement versus either alone; test 0%, 25%, 50% refinement exposure at fixed update budget (22.3.4) |

**Decisive comparison:** partial/full tuning against the best competitive adapter and frozen/native-score baselines. Report allele gains beside degradation on a fixed, separate native-function assessment. Match training exposure or report the additional compute required.

**Output and completion:** a depth-versus-benefit curve and a decision on whether full tuning is justified. Full tuning is a legitimate research arm; scale it after measuring backward-pass memory and development learning curves, without making an adapter win an absolute prerequisite.

### T08. Assay-native losses and uncertainty-aware supervision

**Goal:** optimize the molecular difference of interest while representing the uncertainty and sampling process of each assay.

**Supporting results:** E20, E25–E27, E36.

**Unresolved question:** Does the objective reward the molecular difference rather than shared abundance?

**Data:** Native counts, replicate effects, QTL SEs and valid event counts; no inferred count units.

**Evaluation unit:** Donor/replicate and locus structure appropriate to each assay.

**Success criterion:** Lower native differential error with credible calibration on matched targets; improvement cannot be explained solely by high-depth or strong-effect selection.

**Dependencies:** T01; developed with T05–T07.

**Consequence for the paper:** Make the methodological contribution explicit in the target and likelihood.

**Variables to try:**

| Measurement | Candidate objectives and required interpretation |
|---|---|
| caQTL/eQTL effect estimates | Squared error; Huber; Student-t residuals; Gaussian effect likelihood with variance `SE² + learned residual variance`; cap precision weights within training data |
| Reporter RNA/DNA | Existing replicate-level log-ratio effect; replicate-aware effect likelihood; joint count model with DNA/plasmid abundance as an exposure or latent quantity when raw design information supports it |
| ATAC counts/profiles | Poisson or negative-binomial count model as justified by overdispersion; separate profile-shape and total-count objectives when real profiles exist |
| Allelic read counts | Binomial versus beta-binomial likelihood conditional on total allele-informative reads, with donor structure, genotype support and mapping-bias handling |
| RNA counts | Donor/lineage pseudobulk negative-binomial models with library offsets; continuous-effect losses for already estimated effects. Do not fit count likelihoods to arbitrary transformed values |
| H3K27ac | Documented count likelihood if counts; Gaussian/robust continuous loss if the deposited quantity is normalized signal; donor-differential objective evaluated separately |
| Splicing | Event-effect likelihood for estimated effects; beta-binomial inclusion counts for a valid binary event; a multinomial/overdispersed counterpart for genuinely competing junction categories |
| Allele/state differences | Direct loss on `D` and `J` from section 6; shared baseline plus interaction residual; uncertainty from the actual replicate design |
| Ranking and direction | Secondary pairwise ranking or sign loss alongside magnitude loss; uncertainty-aware labels for weak effects; avoid forcing noisy near-zero estimates into confident sign classes |
| Compatible context sharing | Independent effects versus shared-effect and multivariate shrinkage models; sampling covariance for shared donors/controls; never pool distinct assay scales by default (22.3.2) |
| Donor-contrast weighting | Train-mean residual versus pairwise difference versus joint objectives; 2 or 4 partners per donor per epoch, with donor-balanced loss and donor/locus inference (22.3.5) |

For continuous losses, test Huber transition scales at 0.5, 1 and 2 training residual standard deviations, and Student-t degrees of freedom 3, 5 and 10. For combined objectives, compare normalized auxiliary weights 0, 0.1, 0.3, 1 and 3; normalization is fitted in training data and recorded. These are proposed search values, not borrowed assay defaults.

**Decisive comparison:** matched architecture and data, changing only the objective first. Evaluate native-scale error, calibration, ranking and weak-effect performance. A good profile likelihood does not substitute for a good allele-effect or state-interaction result.

**Output and completion:** an objective recommendation per assay, complete uncertainty assumptions and evidence that improved loss does not merely emphasize high-depth donors or strong selected effects.

### T09. Optimization, augmentation and training schedules

**Goal:** determine whether training choices constrain a promising model and make optimization comparisons distinguishable from changes in data exposure.

**Supporting results:** E32, E34, E36.

**Unresolved question:** Does the chosen training recipe learn reliably at comparable data exposure?

**Data:** Selected model/loss families, measured memory and small-effect numerical diagnostics.

**Evaluation unit:** Independent training groups; optimization seeds reported separately.

**Success criterion:** A stable recipe survives all declared seeds and improves accuracy/cost without selecting on protected outcomes; failures remain counted.

**Dependencies:** T05–T08; hardware-specific pilot.

**Consequence for the paper:** Report reproducible model training and realistic inference requirements.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Optimizer | Existing stable optimizer versus AdamW; momentum-based alternative only if the gradient/cost pilot supports it |
| Effective batch | 16, 32, 64, 128 independent training examples or groups, reached through supported microbatches and accumulation; record the number of distinct loci/donors per update |
| Weight decay | 0, 1e−6, 1e−4, 1e−2, excluding parameters for which decay is inappropriate under the chosen model |
| Learning-rate schedule | Constant, cosine decay, validation plateau; warm-up over 0%, 5% or 10% of the planned updates |
| Stability | Gradient-norm clipping at 0.5, 1 or 5; float32 reductions; rematerialization; loss scaling only where the numeric format requires it |
| Stochastic regularization | Head/adapter dropout; modest input-feature dropout; context-modality dropout; early stopping versus fixed exposure |
| Sampling | Uniform examples, equal locus weight, equal source weight, donor-balanced sampling and task-balanced sampling |
| Sequence augmentation | Reverse complement with correct strand/channel transformation; allele swap with effect-sign reversal; shifts of 1, 8, 32, 64 or 128 bp with target coordinates preserved |
| Context augmentation | Missing-modality masking; small noise calibrated to measured context uncertainty; no invented treatment labels or synthetic donors |
| Curriculum | Profile-to-effect; strong-to-all effects versus all effects from the start; short-to-long context; equal final data exposure in the comparison |
| Evaluation-time averaging | Forward/reverse complement; a small fixed shift set; seed ensemble; report the single-model result and additional inference cost |
| Retention balance | Native mean-profile and allele-contrast consistency weights 0, 0.1, 1 after training-only loss normalization; assess whether contrast preservation prevents the intended correction (22.3.3) |

Sequence shifts cannot move the measured peak/event out of the readout. Random nucleotide changes must not inherit an unchanged measured label. Feature-space interpolation, if explored, is a separate regularization experiment; it does not create measured intermediate genomes or disease states.

**Decisive comparison:** paired development performance at matched numbers of seen training groups, with failure rate, calibration and wall time. Use one or two fixed seeds to identify broken runs; use three for development finalists and the existing five seeds for final internal comparisons. Seeds measure optimization variability.

**Output and completion:** a stable training recipe and an accuracy-versus-cost table. Any failed seed, constant output or nonfinite gradient is retained in the accounting; replacing only failed seeds until the mean improves is excluded.

### T10. Shared learning across assays

**Goal:** learn transferable regulatory information from multiple measurements while retaining each assay's biological meaning and strongest specialist.

**Supporting results:** E27–E28, E36–E37.

**Unresolved question:** Which shared representation helps a receiving assay rather than only the dominant task?

**Data:** Compatible measured assay targets with real pairing and explicit missingness.

**Evaluation unit:** Held loci, donors and publication families as appropriate.

**Success criterion:** The receiving assay improves over an equally tuned separate specialist; all assay losses and negative transfer are reported.

**Dependencies:** T01, T05–T09.

**Consequence for the paper:** Establish transferable regulatory information across assays.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Shared component | Frozen shared encoder; shared adapter; shared early encoder with private later layers; separate encoders with a shared low-dimensional effect representation |
| Assay heads | Independent heads; shared head plus assay conditioning; hierarchical shared/private heads; separate activity, allele-effect and interaction outputs |
| Assay combination | Reporter+ATAC; ATAC+RNA; RNA+splicing; RNA+H3K27ac; broader combination only after the relevant two-assay comparisons |
| Data relationship | Matched variants/targets; paired donors; unpaired assay-specific batches with no invented correspondence |
| Training schedule | Sequential pretraining and adaptation; alternating batches; simultaneous losses; small adapters trained per assay then jointly refined |
| Task balance | Equal normalized loss, equal source/assay sampling, measured-uncertainty weighting, and gradient-conflict handling when supported by the implementation |
| Missing labels | Explicit loss masks; missing-context masks; shared representation learned from available tasks without replacing absent effects by zero |
| Retention | Local teacher consistency on eligible training examples; specialist residual heads; limited replay to prevent one large assay from dominating |
| Measurement hierarchy | Treat episomal reporter, randomly integrated reporter and native-locus assays as different measurements; compare shared effects with assay-specific residuals on a fixed receiving population (22.3.9) |

**Decisive comparison:** leave one supervision source out while holding the target assay, test set, model capacity and training budget as close as feasible. A multitask model must beat the target-assay specialist on at least one prespecified transfer question without concealing losses on the other assays.

**Output and completion:** a directed transfer matrix: which training assay helps which measured endpoint, in which context and with how much uncertainty. If sharing fails, retain separate specialists and report the boundary.

### T11. Reporter-to-endogenous transfer beyond a single extra column

**Goal:** establish whether reporter learning improves native human liver effects and identify the contexts in which that transfer holds.

**Supporting results:** E20, E27–E28, E34–E35.

**Unresolved question:** Can matched reporter information improve native liver prediction?

**Data:** Peak-proximal reporter–caQTL pairs first; additional eligible public functional assays.

**Evaluation unit:** Independent loci; source-family holdout for a transport claim.

**Success criterion:** Improve the endogenous receiving endpoint beyond its native-only counterpart; retain the old extra-column null and quantify overlap/label agreement.

**Dependencies:** T10; T01 bridge census; T23 public overlap expansion.

**Consequence for the paper:** Explain when reporter evidence contributes to endogenous interpretation.

**Inputs:** exact reporter allele effects and measured endogenous variant–peak, variant–gene or variant–event effects. Separate a proximity match from a validated molecular target match.

**Variables to try:**

| Transfer route | Candidate experiment |
|---|---|
| Frozen feature transfer | Reporter-trained local representation appended to endogenous features; compare with the earlier one-column reporter-score route |
| Warm-start | Train a head/adapter on reporter labels, then adapt to endogenous labels; compare with the same amount of endogenous-only training |
| Joint supervision | Shared encoder/adapter with reporter and endogenous heads; vary reporter sampling and loss weight within development data |
| Shared/private decomposition | Shared local allele component plus assay-specific residual; compare with completely separate models |
| Native-context addition | Reporter sequence alone; short native context; long native context; target-aware paired windows and independently supported links |
| Context matching | Matched cell-line conditions where available; adult bulk versus lineage-resolved endogenous labels; source mismatch explicitly tested |
| Distance and target | Fixed peak-proximity strata; actual assayed target; promoter compatibility; gene/event identity rather than nearest gene by default |
| Transfer calibration | Reporter-only model, endogenous-only model, combined model, and measured-label association as a diagnostic rather than a universal ceiling |
| Construct identifiability | Full construct/promoter factors only when recorded and actually varied; compare matched versus equally sized unmatched supervision, preserving endogenous training size (22.3.9) |

**Decisive comparison:** held-locus endogenous error with versus without reporter supervision, including training loci for every auxiliary reporter head. Remove evaluation loci and their linked/overlapping sequence from all contributing training tasks. Do not force the two assays' effect signs or scales to agree.

**Output and completion:** a transferable subset defined before its independent test, a complete transfer-effect table and failure explanations grounded in measured assay/context differences. The current 209 variants within 500 bp in 89 blocks are development evidence, not enough by themselves to promise a precise interaction or external-transfer claim; size the next matched dataset using T23.

### T12. Adult liver cell-type heads

**Goal:** add measured adult hepatocyte, stellate/fibroblast, endothelial, cholangiocyte and immune outputs where data support them, including capabilities absent from the pretrained catalogue.

**Supporting results:** E04, E22, E26, E36–E37.

**Unresolved question:** Which adult lineages support profile learning and which also support allele effects?

**Data:** Donor/lineage RNA and ATAC profiles; eligible lineage allelic counts.

**Evaluation unit:** Donor × lineage, with held genomic targets; cells improve precision within donors.

**Success criterion:** New heads improve held-out profiles where measured; allele-effect claims require separate allele-label success, not a profile surrogate.

**Dependencies:** T01–T02; public source and track-age reconciliation.

**Consequence for the paper:** Add adult cellular context without overstating sparse-lineage validation.

Absence of a pretrained adult track does not prohibit training a new head from appropriate measurements. It does mean that the existing static model cannot already validate that cell-type claim. Multimodal cell-specific decoding is an established direction; the new question here is adult human liver allele and state transport. [scooby](https://www.nature.com/articles/s41592-025-02854-5)

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Supervision level | Donor-by-lineage pseudobulk; donor-aware hierarchical cell model; native bulk profiles as an auxiliary task |
| Decoder | One head per lineage; shared head with lineage embedding; low-rank lineage factors; shared baseline plus lineage residual |
| Parameter sharing | Shared backbone and private heads; shared adapter and private heads; lineage-specific adapters; hierarchical partial pooling for sparse lineages |
| Measurement | Profile/count prediction first where measured; separate allele-effect head where genotypes and allele counts/effects exist |
| Region sampling | Outcome-independent eligible peaks/regions, expression/accessibility strata and random eligible controls; do not restrict to Atlas-served variants alone |
| Cohort balance | Equal donor/source contribution, library-depth handling and a cohort-aware training sensitivity |
| Context support | Observed lineages; related-lineage transfer as an explicitly new test; missing lineage returned as unsupported until assessed |
| Consensus-to-lineage learning | Compact specialist and adapted backbone under the same consensus/specific-region curriculum; test broad regions and lineage-specific regions separately (22.3.4) |

**Decisive comparison:** new-donor and new-region performance versus the lineage mean, a simple native-assay model and the strongest existing head. For allele claims, evaluate actual allele labels separately. Thousands of cells from a few donors cannot substitute for a held-donor test.

**Output and completion:** an adult lineage-by-assay capability table, with distinct profile and allele-effect evidence. Broad profile learning may proceed while rare-lineage allele labels are developed; completion of one does not imply completion of the other.

### T13. Cellular-context inputs and conditioning architecture

**Goal:** identify which measured context makes sequence predictions more informative and which apparent gains are attributable to study, depth or composition.

**Supporting results:** E05, E22, E36–E39.

**Unresolved question:** Does measured context add molecular information beyond nuisance structure?

**Data:** Training-only RNA/regulator representations, unchanged program projections, lineage/treatment metadata and admissible paired modalities.

**Evaluation unit:** Crossed held donors and targets; source/well-preserving controls.

**Success criterion:** Beat both sequence-only and context-only baselines; the interaction increment survives additive, shuffled-context, technical and composition controls.

**Dependencies:** T12 or admitted paired bulk profiles; T08 objective.

**Consequence for the paper:** Distinguish context-conditioned regulation from donor/study recognition.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Minimal context | Lineage and experimental treatment labels; recorded histology components only where available |
| Molecular context | Transcription-factor/regulator expression; training-selected RNA genes; donor/lineage RNA components; the frozen 117 program scores; assay-specific chromatin summaries |
| Combination | RNA alone; programs alone; RNA plus lineage; RNA plus measured chromatin; metadata plus molecular context; explicit missing-modality masks |
| Encoder | Linear projection; small MLP; attention over a fixed regulator set; training-only low-dimensional bottleneck |
| Cross-system representations | Eligible paired-modality/topic inputs motivated by MIRA; compare original representation, donor-aware MASLD adaptation and RNA-only inputs with target exclusion; no atlas reintegration |
| Context dimension | 4, 8, 16, 32, 64, subject to donor rank and regularization; gene inputs of 128, 512 or 2,048 selected only within training |
| Fusion | Late concatenation; feature-wise affine modulation (FiLM); low-rank sequence×context interaction; target-to-context cross-attention; shared plus context-specific residual |
| Placement | Head only; selected encoder layers; selected transformer layers; decoder; compare early versus late conditioning |
| Missing/shifted context | Context dropout 0, 0.1, 0.3; measured-context versus metadata-only inference; training-mean substitution as an explicit baseline |
| Chromatin input contract | RNA-only, chromatin-only and paired context; distinguish observed-cis-chromatin RNA imputation from strict target-excluded context for donor/allele prediction (22.3.6) |

Corgi combines sequence features and regulator expression through FiLM, providing a concrete conditioning comparator. Its existence does not establish that our donor-differential or allele-interaction endpoint is learnable with the present data. [Corgi](https://www.nature.com/articles/s41467-026-75527-2)

**Decisive comparison:** sequence-only, context-only, additive sequence+context and sequence×context models on the same endpoint. Compare real context with outcome-blind shuffled or mismatched context under a design that preserves relevant source structure. Include measured depth and composition sensitivities; report unavailable age/BMI rather than inventing adjustment.

**Output and completion:** context ablations and evidence that the model uses transferable molecular information. Exclude target RNA and target-containing derived features when those would reveal the prediction target; fit every context transformation in training data. Cohort ID is not a transferable biological input. Holding RNA context fixed while changing genotype remains a conditional prediction, not a feedback simulation.

### T14. Allele-by-treatment and allele-by-state effects

**Goal:** predict `J = allele effect in context 1 − allele effect in context 0`, extending the program from average effects to measured regulatory response.

**Supporting results:** E05, E20, E26.

**Unresolved question:** Can the model predict a reliable change in allele effect across contexts?

**Data:** LX-2 treatment contrasts first; additional public state-interaction labels after reliability review.

**Evaluation unit:** Experimental replicate and locus for MPRA; donor and locus for human states.

**Success criterion:** Lower native interaction error than additive/no-interaction and separate-condition models; a human state claim additionally needs human replication.

**Dependencies:** T08, T13; T01 reliable interaction labels.

**Consequence for the paper:** Demonstrate allele-by-context capability without equating treatment with fibrosis stage.

**Starting endpoint:** the LX-2 treatment contrast, whose current reliability supports development. HepG2 interaction learning remains a measurement-development task unless a new, documented design supplies a reliable label.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Prediction target | Four allele×condition measurements jointly; two condition-specific allele effects; direct interaction residual with a shared allele baseline |
| Context | Treatment identity; measured regulator/program state; both; conditional effects within recorded histology strata when individual molecular labels exist |
| Model | Separate condition specialists; shared conditioned model; low-rank interaction head; adapter conditioned on treatment/state |
| Loss | Replicate-aware likelihood; direct interaction error; joint activity/effect/interaction loss with prespecified weights |
| Shrinkage | No interaction regularization versus regularization toward a shared allele effect; precision-aware partial pooling |
| Generalization | New loci within the two observed conditions; new treatment/state only after additional independently measured contexts exist |
| Specificity | Additive-context predictor, no-context predictor, matched treatment permutation and uncertain-effect strata |
| Effect-sharing model | Independent, shared and partially pooled context effects; covariance ranks 1, 2, 4 only when context count and replicate structure justify them (22.3.2) |

**Decisive comparison:** native interaction error versus a model that cannot express an allele×context term, plus comparison with independently fitted condition models. Evaluate correlation only where the interaction signal is estimable. Two measured conditions can support new-locus interaction prediction, but cannot validate broad unseen-treatment generalization.

**Output and completion:** a measured-versus-predicted interaction table and an independently evaluated improvement or a clearly scoped negative result. A treatment contrast in LX-2 is not a donor fibrosis-stage effect; distinguish those questions in the figures and follow-up design.

### T15. Donor-differential RNA and chromatin

**Goal:** predict biologically meaningful differences among donors beyond the shared genomic profile and a strong RNA/context-only model.

**Supporting results:** E36–E39.

**Unresolved question:** Does sequence explain donor-differential regulation beyond the strong RNA-only predictor?

**Data:** Paired GSE267145 and eligible GSE296875/GSE244832 profiles, with exact assay and donor joins.

**Evaluation unit:** Held participants × held genomic targets, with study-level assessment separately.

**Success criterion:** Lower donor-differential error beyond matched RNA-only/reduced-rank baselines; survive technical, well and composition controls.

**Dependencies:** T01, T12–T13; rights-specific public/research model route.

**Consequence for the paper:** Connect model development to measured remodeling while keeping genotype effects distinct.

**Inputs:** the paired GSE267145 RNA–H3K27ac cohort and other authoritative donor-linked modalities where eligible. Reconcile source overlap before describing replication.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Target | Native profile; residual from the training-donor mean; deviation within recorded histology strata; donor-paired cross-assay effects |
| Model decomposition | Mean sequence profile plus donor residual; RNA-only baseline; sequence+RNA late fusion; FiLM; bilinear or cross-attention conditioning |
| Objective | Total profile loss; explicit donor-difference loss; combined objective with normalized weights; error for donor contrasts as an auxiliary term |
| Covariates | Recorded histology components, sex, library descriptors, available age/BMI and cell composition; retain missingness and study-specific availability |
| RNA representation | Small training-only components; selected cis genes; regulator expression; fixed programs; cis-excluded and target-excluded sensitivities |
| Evaluation | Held donor; held region; crossed donor×region holdout; independent cohort when comparable data are available |
| Pair sampling and transport | 2/4 partners per donor versus affordable all-pairs training; same-source and technically matched pairs; seen/new donor crossed with seen/new locus assessment (22.3.5) |

**Decisive comparison:** donor-differential error against the strongest matched RNA-only/reduced-rank baseline and the existing Corgi recipe. Calculate training means and residual transformations without test donors. Compare models on exactly the same target definition; the several existing H3K27ac scores are not interchangeable.

**Output and completion:** separate mean-profile and donor-differential results, composition/technical sensitivities and crossed-holdout performance. A successful donor-context model is a new capability; it still does not establish genotype-dependent effects when individual genotypes are absent.

### T16. Target-gene assignment and expression effects

**Goal:** move from accessibility prediction to the affected molecular target, and test whether sequence/context information resolves target ambiguity at inherited loci.

**Supporting results:** E02, E07, E16, E22.

**Unresolved question:** Can predictions distinguish the target and signed expression effect within a locus?

**Data:** Measured gene-level effects and independently supported region–gene links; competing candidate genes retained.

**Evaluation unit:** Held loci/genes, with repeated targets clustered.

**Success criterion:** Improve correct-target discrimination and native expression-effect prediction beyond proximity, genetic evidence and the strongest eligible model, on independent labels.

**Dependencies:** T02, T10, T12; T23 target-label admission.

**Consequence for the paper:** Strengthen a locus explanation without equating accessibility with target-gene mediation.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Output | Gene expression effect, promoter activity, transcript-specific effect and a separately evaluated region–gene relationship |
| Sequence | Variant-centred local window; promoter window; long native context; paired variant/promoter encoders; selected distal elements |
| Target information | Distance-only baseline; measured accessibility; independent enhancer–gene/contact evidence; sequence-derived contact prediction recorded as a predictor rather than truth |
| Head | Target-aware attention; low-rank variant×target interaction; graph aggregation over independently defined candidate links; separate gene/assay heads |
| Network and genetic links | SCENIC+-style TF–region–gene networks; sc-linker-style SNP-to-gene maps with appropriate genome-wide inputs; original versus sequence-augmented links and proximity; fit links without held-out target outcomes |
| Supervision | Native eQTL effects with alleles reconciled; genotype-linked donor RNA; measured regulatory links; avoid using the same colocalization outcome as both feature and validation target |
| Mechanism class | Coding, regulatory and RNA-processing candidates assessed separately; distinguish lncRNA locus activity from mature-RNA function |
| Model information | Sequence-only and annotation-assisted versions reported separately; genetic posterior used for locus interpretation, not silently converted into supervised effect truth |
| Linker comparison | Proximity, current annotations, ENCODE-rE2G, scE2G-ATAC and eligible scE2G-Multiome; add sequence effects as separate features before considering larger graph models (22.3.7) |
| Network stability | Unaugmented versus allele-augmented TF–region–gene networks; donor-resampling stability and degree/activity/motif-matched feature controls (22.3.8) |

**Decisive comparison:** gene-level effect prediction and within-locus target discrimination beyond proximity and the strongest eligible native predictor, at held loci. Compare the target-aware model with the same sequence encoder lacking target information. Audit GTEx and other pretraining exposure for the actual endpoint.

**Output and completion:** candidate variant–target effects with signed assay units, uncertainty and competing targets, plus independent target/event tests. Preserve coding explanations and unresolved noncoding signals; a successful accessibility model does not establish the gene or RNA molecule mediating disease risk.

### T17. Splicing, transcript events and indels

**Goal:** develop a separately calibrated RNA-processing model that predicts the affected event and its signed effect, including variants not adequately described by total-gene expression.

**Supporting results:** E09–E10, E14–E15.

**Unresolved question:** Can event-aware supervision improve detection, direction and magnitude separately?

**Data:** Exact event/junction measurements; development GTEx results with exposure declared; independently eligible public event effects.

**Evaluation unit:** Held genes/loci; indels evaluated separately from SNVs.

**Success criterion:** Report three endpoint results; any new ranking claim must pass a new fixed independent test, while the original step-36 D1 remains recorded.

**Dependencies:** T01–T08; completed E10; event/annotation alignment and source-native denominator checks.

**Consequence for the paper:** Resolve RNA-processing alternatives with measured event evidence rather than unsigned scores.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Endpoint | Event detection, junction identity, usage change, intron retention, exon inclusion, alternative promoter or transcript-end use; each requires its own measured label |
| Input | Local donor/acceptor/exon context; exon plus flanking introns; larger transcript context; supported long input for competing events |
| Representation | Targeted 1-bp decoder features; 128-bp context features; multiresolution fusion; native splice-site/junction scores plus a learned residual |
| Adaptation | New event head; decoder-only tuning; encoder/decoder adapters; broader tuning if measured labels and cost support it |
| Allele type | SNVs, short insertions/deletions and explicitly supported complex alleles; matched eligibility and correct coordinate remapping |
| Target conditioning | Transcript/event coordinates, strand and competing splice sites; outcome-independent annotation version; multiple events retained where appropriate |
| Event denominator | Source-native measured cluster where recoverable; fixed annotation-defined competing junctions; reference-usage-filtered model junctions; retain the loose step-36 definition as the exposed baseline |
| Loss and calibration | Event-native effect likelihood; count-based usage model where counts exist; detection and signed-direction performance reported separately |
| Generalization | Held gene/locus; indel-versus-SNV strata; independent RNA source; exposure-matched event-specialist comparison |
| Assay-native specialists | Pangolin liver site usage; MTSplice-compatible cassette exons; APARENT2 cleavage/APA labels; each evaluated on its own identifiable event universe (22.3.10–22.3.11) |
| RNA context | TF-only versus externally defined RNA-binding-protein context where paired data support it; shared baseline plus tissue/event residual (22.3.11) |

**Decisive comparison:** correct-event signed prediction beyond the best eligible event-specific baseline on held genes/loci, together with detection precision–recall where labels support it. Include proximity for event detection and the E10 signed route for direction. Fix annotation, usage rule and competing-event universe before each new development comparison; REF and ALT use the same identifiable event universe, with newly created/lost junctions handled explicitly. A source-native cluster is not guaranteed to equal the set of model junctions sharing a boundary. Predicting a maximum unsigned splice score is not signed event prediction. Do not interpret a standardized sQTL effect as a change in percent spliced in unless that is the actual source unit.

**Output and completion:** an event-resolved benchmark and model, with coverage by allele class, sign convention and transcript. Single-cell gene-count matrices cannot supply junction coverage that was never measured. Existing inspected examples can be development checks; they cannot be the independent test of a model selected to recover them.

### T18. Donor sequence, genotype and haplotype background

**Goal:** measure whether real donor sequence improves molecular predictions and separate cis sequence information from cellular context.

**Supporting results:** E14–E15, E26, E29–E30.

**Unresolved question:** Does actual cis background improve prediction beyond dosage and measured context?

**Data:** Existing public donor genotypes/phased haplotypes linked to assays; source availability must be established.

**Evaluation unit:** Held donors × held loci; related haplotypes grouped.

**Success criterion:** Improve matched molecular prediction with real background information; propagate phase uncertainty and retain reference-sequence controls.

**Dependencies:** T01; public genotype/phase/assay admission through T23.

**Consequence for the paper:** Support personalized-sequence capability only where real linked data exist.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Sequence background | Reference-only; focal allele on reference background; verified local common variants; full supported local phased haplotype |
| Input decomposition | Two haplotypes separately; explicitly defined diploid aggregation; focal-allele contrast on the same donor background |
| Context | No donor context; measured RNA/regulator state; genotype plus context; context-only comparator |
| Phase | Verified phase; phase uncertainty propagated over plausible configurations; unphased sites reported separately |
| Genotype certainty | High-confidence genotypes; dosage/likelihood-aware sensitivity where available; read-based calls retained only as a separately named comparison |
| Neighborhood | Local peak, gene/transcript and longer cis region; matched coverage across background choices |
| Learning | Frozen scores; calibrated scores; background-aware residual head; adapter training on genotype-linked molecular effects |
| Diploid readout scale | Aggregate haplotype outputs in the assay's physical scale before logs/ratios; compare reference, local variants and actual supported haplotypes (22.3.12) |

**Decisive comparison:** molecular effects in held donors and held loci with versus without the actual sequence background. Use identical measured context for the focal cis comparison. Evaluate incremental value beyond genotype dosage and measured-context baselines, while acknowledging that context may itself reflect genotype and disease state.

**Output and completion:** a personalized-sequence evaluation, genotype/phase uncertainty accounting and a documented boundary on what can be predicted. The current 30-donor read-based allelic set does not by itself supply all required genotypes or phase. Add those inputs through T23 before claiming this capability.

### T19. Measured haplotype and haplotype-by-state interactions

**Goal:** predict joint molecular effects beyond a prespecified additive baseline and test whether those interactions change with context.

**Supporting results:** E19, E29–E30.

**Unresolved question:** Do measured joint effects contain a predictable interaction beyond additivity?

**Data:** Existing public joint/single-variant measurements with phase and identical readout geometry.

**Evaluation unit:** Distinct variant pairs and independent loci; donor/replicate pairing when measured.

**Success criterion:** Improve measured joint-effect and interaction prediction over additive controls; prediction-only analysis remains labeled if no eligible public labels exist.

**Dependencies:** T18 and adequate public combination labels; otherwise retain a conditional data requirement.

**Consequence for the paper:** Report supported molecular combinations without claiming population disease epistasis.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Combination | Observed two-variant haplotypes first; larger observed combinations only where measurements and phase provide enough independent information |
| Model | Sum of single effects; direct joint-sequence prediction; additive baseline plus learned interaction residual; low-rank pair model with distance/context input |
| Context | One measured context; multiple replicated contexts; haplotype-by-state residual only when all required contrasts are observed |
| Readout | Same peak, gene or event across all four allele combinations; native and transformed scales analysed separately with a primary scale fixed |
| Mechanism features | Separation, target identity, motif arrangement and measured context; compare features alone with sequence adaptation |
| Null and controls | Match readout geometry, baseline activity, single-arm magnitude, separation and ascertainment; match within a selected active-arm stratum if that stratum is tested |
| Dependence | Unique variant-pair/target hypotheses clustered by locus; repeated COLOC rows and several readouts of one pair are not independent discoveries |
| Adaptation | Joint-effect head; small adapter; sequential single-effect then interaction training; joint loss with strong additive regularization |
| Background versus interaction | Separate ordinary background effects from measured joint-interaction contrasts; propagate phase uncertainty and retain the same readout geometry (22.3.12) |

**Decisive comparison:** error on independently measured joint effects and on `I` from section 6, beyond the additive predictor. Evaluate the joint effect as well as the residual so a noisy subtraction does not conceal useful or poor prediction. Include uncertainty from all measured arms and the actual experimental pairing.

**Output and completion:** a measured interaction benchmark with phase, native units and independent-locus evaluation. A local in-silico null can characterize model nonlinearity, but it cannot substitute for measured combinations. No significant residual is required at GNMT or any other predetermined example; failure to reject non-additivity is not proof of additivity.

### T20. Generalization, source shift and adaptation to new contexts

**Goal:** determine which improvements survive the specific deployment setting claimed for the Resource.

**Supporting results:** E03–E07, E24, E33, E37–E40.

**Unresolved question:** Which generalization claim survives source overlap and prior exposure?

**Data:** Protected GSE277774 candidate, restricted GSE289173 evaluator role, and other prospectively admitted sources.

**Evaluation unit:** Explicit locus, donor, study, state or haplotype separation for each claim.

**Success criterion:** The locked comparison succeeds on the claimed independent units with a complete exposure/attrition account; otherwise narrow the claim.

**Dependencies:** Split design begins with T01; final evaluation follows fixed recipes from relevant tasks.

**Consequence for the paper:** Make Figure 7's transfer claims independently interpretable.

**Variables to try:**

| Generalization axis | Candidate experiments |
|---|---|
| Locus | Whole-chromosome or locus-group holdout; LD/overlapping-window buffering; gene-family sensitivity where homologous sequence could explain transfer |
| Donor | All cells, aliquots and repeated measurements kept together; crossed donor×locus holdout for claims requiring both |
| Study | Leave-one-source-family-out development; independently protected final source; compatible target definition required |
| State | Leave-context-out where at least enough distinct contexts exist to train and test the claim; interpolation and extrapolation reported separately |
| Ancestry/background | Evaluation by documented ancestry and local sequence/LD support; matched precision and allele-frequency analyses; ancestry is not itself a molecular state |
| Domain training | Source-balanced sampling, shared/private adapters, hierarchical source effects and carefully assessed domain-invariance objectives |
| New-context adaptation | Zero-shot transfer; small support-set head calibration; adapter update on a support set; held query donors/loci untouched by that update |
| Model combination | Fixed seed average; training-only stacked specialist scores; two- or four-expert mixture with context-dependent weights; compare with the best single model |
| Pretraining | All-folds local model with exposure declared; official fold-specific model where exact held genomic intervals can be reconciled; same adaptation/test boundaries within each comparison |

Domain-invariance training cannot resolve a design in which treatment, lineage or disease state is perfectly confounded with source. Official AlphaGenome pretraining folds and the project's five chromosome groups are different partitions; map their actual intervals rather than treating their fold numbers as interchangeable.

**Decisive comparison:** predeclared contrasts on the genuinely protected units. Report each source and subgroup, not just a pooled score. For few-shot adaptation, report the number of support donors and loci and keep query data out of every tuning step.

**Output and completion:** a generalization matrix, identifying which claims have locus, donor, study, state and pretraining support. A new prespecification after inspecting outcomes does not make the same outcomes untouched. If independent compatible data remain unavailable, report internal development and the exact missing evaluation rather than manufacturing external confirmation.

### T21. Uncertainty, calibration and selective prediction

**Goal:** provide useful uncertainty and identify contexts in which a model should return an unsupported or low-confidence prediction.

**Supporting results:** E06, E20, E25–E26, E32, E38.

**Unresolved question:** When are predicted effect intervals and applicability useful?

**Data:** Separate calibration data and locked native-effect predictions, with uncertainty in measured labels.

**Evaluation unit:** Independent calibration donors/loci; source-specific evaluation.

**Success criterion:** Evaluate interval coverage/width and error versus abstention coverage with uncertainty; retain unsupported contexts rather than imputing confidence.

**Dependencies:** T20 split design; trained models and adequate calibration units.

**Consequence for the paper:** Let database users assess the precision and supported domain of a prediction.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Effect uncertainty | Sampling-variance-aware head; heteroscedastic residual model; seed ensemble; training-group bootstrap where computationally feasible |
| Calibration | Separate training/calibration/test units; group-aware residual calibration; conformal approaches only with their exchangeability assumptions and grouping stated |
| Calibration target | Latent molecular effect versus a future noisy measurement; distinguish these intervals explicitly |
| Direction | Calibrated probability of positive effect; three-way uncertain/positive/negative summaries with thresholds fixed before test |
| Applicability | Distance in training context/feature space; assay/lineage coverage; genotype/phase uncertainty; missing input and sequence eligibility |
| Selective prediction | Error versus retained population coverage; effect-size and lineage-stratified coverage; predefined rejection rules |
| Shift | Calibration within source; transported calibration; limited new-source support calibration evaluated on separate query units |
| Interval and panel stress test | 80%, 90%, 95% nominal coverage with widths; recorded gene-panel masks and missing modalities; calibration in disjoint donor/locus groups (22.3.16) |

**Decisive comparison:** interval coverage and width, native-scale error, direction calibration and the cost of abstention on protected data. A nominal 95% interval is evaluated with uncertainty on its empirical coverage; five optimization seeds are not five biological replicates and do not establish calibrated intervals.

**Output and completion:** calibrated prediction tables and applicability rules, with the calibration population stated. Distribution shift can invalidate coverage guarantees; measure that failure instead of presenting an out-of-distribution score as a known probability of biological error.

### T22. Regulatory interpretation and discriminating biological tests

**Goal:** explain supported predictive gains in terms of testable molecular alternatives and identify where the model adds information beyond existing annotations.

**Supporting results:** E07–E08, E12, E16–E19, E22–E23, E30, E39–E40.

**Unresolved question:** Which molecular explanation predicts an independent observation better than its alternative?

**Data:** Full P4 when complete, measured molecular effects, fixed programs, eligible protein/spatial evidence and public functional results.

**Evaluation unit:** Independent loci for regulatory tests; donors for tissue associations.

**Success criterion:** A prespecified relationship survives active-region/co-expressed-set and technical controls, with independent support for any mechanism claim.

**Dependencies:** Relevant molecular tasks; E17 completion before full-P4 claims.

**Consequence for the paper:** Connect model gains to inherited regulation and cross-sectional remodeling.

**Variables to try:**

| Axis | Proposed settings |
|---|---|
| Interpretation | Local saturation; gradient/integrated-gradient attribution; sequence occlusion; comparisons across model seeds, shifts and architectures |
| Readout | Peak-local response; target-gene response; exact transcript event; interaction residual; each interpreted separately |
| Motifs | Known motif-family comparison; reproducible candidate sequence patterns; positional/co-occurrence analysis; distinguish motif compatibility from demonstrated TF occupancy |
| Baselines/nulls | GC, baseline activity, accessibility, promoter/TSS distance, region width and locus/LD structure; include measured signal mean/variance for saturation-skill questions |
| Regulatory response | Sequence sensitivity, positional concentration and donor-context dependence; test buffering/amplification as predictions of measured contrasts, keeping partial B3 and full P4 separate |
| External context | Independent TF/chromatin measurements, measured region–gene links, donor-aware protein/spatial observations and the fixed 117 programs |
| Example selection | A fixed diverse development set; model-disagreement cases; random eligible controls; independently locked validation examples |
| Biological alternatives | Local regulatory versus target-link effect; total RNA versus event effect; coding versus regulatory route; lncRNA locus versus RNA-product hypothesis |
| Cross-disease adaptation | Original network/program method, donor-aware MASLD version and addition of sequence-effect information; section 2.1 candidates require compatible public validation and unchanged fixed117; CellOracle-like regulator responses remain hypotheses without measured support |
| Multicellular hypothesis | Fixed117 cross-lineage prediction beyond stage/composition and within-lineage baselines; rank 1/2/4 shared components only with adequate donors (22.3.13) |
| Spatial support | Expression-only versus spatially supported signaling hypotheses; 50/100/200 micrometre neighborhoods only with valid coordinates and appropriate resolution; donor/compartment-preserving nulls (22.3.14) |
| Distal context dependence | Reuse P4, then native context versus nominated-region and matched-control perturbations on a bounded locus set; same output geometry and independent support (22.3.15) |

**Decisive comparison:** interpretations must be stable enough to nominate a discriminating measurement, and the chosen molecular explanation must predict an independent observation better than its competitor. A motif map or attractive saturation plot alone is not mechanism validation.

**Output and completion:** a bounded set of locus explanations with explicit alternatives, measured support and suitable follow-up readouts. Preserve source-owned discoveries and fixed Resource examples; new exploratory examples do not silently replace adopted manuscript choices.

### T23. Learning curves, label acquisition and experimental design

**Goal:** determine whether the next investment should be donors, loci, read depth, replicated states, matched assays or model capacity, and make ambitious data expansion quantitatively reviewable.

**Supporting results:** E05–E07, E20, E25–E27, E38; public priorities in section 5.1.

**Unresolved question:** Which publicly obtainable information most reduces the current error or uncertainty?

**Data:** Training-size curves, label precision, locus/donor coverage and public metadata; no protected outcome inspection.

**Evaluation unit:** Donors and independent loci separately; experimental replicates distinguished.

**Success criterion:** A ranked public-source proposal states expected precision, overlap, reuse, costs and role; unavailable information remains an explicit limit.

**Dependencies:** Starts with T01; update after each campaign's informative comparison.

**Consequence for the paper:** Expand the paper's evidence without introducing new experimental dependencies.

**Variables to try:**

| Information axis | Candidate study |
|---|---|
| Training loci | Existing 800/2,000/5,000/12,500/full comparisons plus nested fractions appropriate to the current task; hold evaluation loci fixed and report undefined metrics |
| Biological donors | Donor subsets only where individual-level molecular/genotype data exist; refit label estimation and preprocessing within each subset as required |
| Experimental replicates | Replicate subsets and uncertainty models for activity, allele effect and interaction separately; preserve real pairing |
| Depth | Read subsampling and donor-site coverage curves for each lineage; this estimates losses below observed depth, not guaranteed gains from arbitrarily deeper sequencing |
| Matched assays | Increase the number of correctly paired variants/targets and independent loci; compare gains with adding the same number of unpaired examples |
| State diversity | More reliably replicated contexts versus more observations of existing contexts; assess identifiability of new-state prediction |
| Population breadth | Independent source, ancestry/LD backgrounds and disease-component coverage rather than repeated rows from the same discovery population |
| Model capacity | Cross a small number of data sizes with frozen, adapter and partially unfrozen models to distinguish label limitation from representation limitation |
| Acquisition policy | Random eligible sampling; uncertainty; model disagreement; biological/context diversity; mixtures with a prespecified random component, for example 20–30% |
| Public information value | Prioritize independent target/event labels, assay bridges and donor-linked backgrounds from the source routes in 22.4; source admission precedes download or performance inspection |

**Data-development priorities to compare:** publicly available genotype support for existing allele-assay donors; independent adult liver molecular-effect data; correctly matched reporter/endogenous variants; reliable replicated treatment contrasts; wider rare-lineage allele coverage; and phased combinations with measured effects. The retained replicate/depth/design comparisons describe information requirements and subsampling of existing data. New sequencing, genotyping, libraries and prospective experiments are not dependencies. An unavailable public input stays a recorded limitation.

**Decisive comparison:** use paired learning curves and native uncertainty to estimate the improvement and precision expected from each next increment. For new-data design, report ranges under alternative effect distributions, ascertainment and locus dependence. The significance-selected Pearson ceiling and a fitted slope beyond observed training size are insufficient to dismiss donor acquisition or guarantee a gain.

**Output and completion:** a concrete ranked public-data proposal with biological units, target labels, expected precision, source/access feasibility, transfer/storage size and compute cost. Keep development acquisitions separate from the independent final test; model-selected examples alone cannot establish population-wide performance. Experimental requirements may be described as future research needs, not scheduled work.

### T24. Research models, response tables and Resource integration

**Goal:** turn the successful capabilities and informative limits into a usable molecular-response resource with reproducible models and appropriate experimental follow-up.

**Supporting results:** E01, E13–E16, E21–E22, E31–E40.

**Unresolved question:** Can a researcher retrieve and reproduce the supported regulatory answer?

**Data:** Existing Catalog/Atlas layers, completed model evaluations and eligible public model artifacts.

**Evaluation unit:** Database record identities plus the biological units underlying each evidence link.

**Success criterion:** All six query entry points return traceable effects, coverage, support and limitations; an eligible model reproduces declared predictions locally.

**Dependencies:** Completed relevant tasks; sections 20–21; source-specific distribution review.

**Consequence for the paper:** Deliver the database and models as usable contributions to the existing Resource paper.

**Variables to try:**

| Axis | Candidate deliverables/comparisons |
|---|---|
| Model organization | Best specialist per assay; shared model plus specialists; context-conditioned mixture where T20 supports it |
| Inference cost | Full model; frozen-feature head; compact student trained from eligible local outputs; distillation assessed on effect error, calibration and context support, under applicable derivative terms |
| Response record | Variant/haplotype, target/event, allele orientation, native effect, uncertainty, context, eligibility, measured support and alternative explanation |
| Research interaction | Single-variant and phased-background queries; comparison across supported contexts; explicit response when an input context is unsupported |
| Gene Catalog | Candidate attachment preserving existing genetic evidence classes, direct-disease versus enzyme traits, and object-specific follow-up rules |
| Fixed programs | Projection onto unchanged 117 memberships/weights; donor-aware context association; no claim that membership establishes variant-to-program causation |
| Publication | Result-led molecular comparisons and complete negative/coverage tables; actual data and model reuse instructions appropriate to their terms |
| External evidence identity | Preserve source dataset/event/score-set identifiers and units; expose duplicate origins and source-specific permitted local reproduction (R27–R28; 22.3.16) |

**Decisive comparison:** someone using the research artifacts can reproduce the declared effect and its uncertainty from the documented inputs, distinguish prediction from measurement, and identify what observation would resolve the remaining molecular alternative.

**Output and completion:** permitted model artifacts, inference and training recipes, complete response tables, a benchmark with explicit exposure, and a candidate Gene Catalog attachment. Manuscript adoption remains a separate explicit decision; no new gene class, threshold, program or figure is adopted by writing or executing this plan.

The complete database behavior, minimum record structure, public/research model distinction and acceptance examples are specified in section 20. Database and inference design is a scientific deliverable of T24, not an optional interface added after training.

## 18. How to search this space

### 18.1 Broad exploration with protected evaluation

The tables intentionally cover many plausible variables. They should not be multiplied into one exhaustive Cartesian grid. Most settings are conditional on the assay, architecture or data available, and a large search on a repeatedly inspected test would measure selection advantage rather than generalization.

Use the following experimental sequence:

1. **Identify the claim and its unit.** Write the target `D`, `J`, donor difference or event effect; list the actual donors, loci and assays; select the strongest eligible comparator.
2. **Separate development, calibration and protected evaluation.** Existing inspected benchmarks remain development. Preserve an independent source if available; otherwise distinguish nested internal evaluation from external evidence.
3. **Test broad architectural families on development data.** Include native-score models, frozen probes, adapters and at least a bounded partial-unfreezing pilot. A poor frozen probe is not a reason to omit all adaptation.
4. **Search promising families more deeply.** Fit normalization, compression, features, losses, hyperparameters and stopping inside the allowed training/development partitions. Search the same important variables for competing baselines.
5. **Test interactions between important variables.** Use the explicit factorial comparisons below; do not assume a pooling choice selected for one input or assay will remain best after fine-tuning.
6. **Select final recipes and margins.** Fix model, input, endpoint, primary comparator, minimum useful improvement, uncertainty method and multiplicity family before opening protected outcomes.
7. **Run the locked comparison once.** Report all planned seeds, exclusions, native-scale errors, correlations, confidence intervals and failures. Additional model selection starts a new development cycle and needs a new protected assessment.

### 18.2 Interactions worth testing explicitly

| Interaction | Why the combined experiment matters | Minimum informative design |
|---|---|---|
| Input length × pooling | The reporter result suggests the representation can discard a local effect; a longer input may help only some pooling schemes | Two native lengths × variant, symmetric local and target pooling, on matched examples |
| Layer × assay | Accessibility, total RNA and exact splice events may use different representations | Selected encoder/transformer/decoder features across the relevant native endpoints |
| Pooling × adaptation | A frozen-feature winner may change after the backbone learns the scalar effect | Same two pooling choices with frozen and adapter models |
| Adapter rank × depth | A larger update in few layers can differ from a small update throughout the network | Two ranks × two depths at comparable parameter budgets |
| Context input × fusion | Rich RNA inputs may require different conditioning than a lineage label | Minimal versus molecular context crossed with late versus internal conditioning |
| Objective × training universe | A precision-weighted loss may behave differently on selected leads and an outcome-independent eligible population | Two losses × both populations; common evaluation and ascertainment reported |
| Profile loss × effect loss | A dominant profile objective can obscure poor allele or donor prediction | Profile-only, effect-only and joint training with normalized weights |
| Data size × trainable depth | Full tuning may become useful only after enough independent labels exist | Small/intermediate/full training subsets crossed with frozen/adapter/partial tuning |
| Assay sharing × context match | Apparent negative transfer may arise from incompatible biological contexts | Matched versus mismatched contexts crossed with separate versus shared training |
| Genomic background × state | A context gain may change when actual cis variation is supplied | Reference versus donor sequence crossed with absent versus measured context |
| Precision × weak-effect accuracy | Numerically close profiles can differ in small allele contrasts | Selected higher-precision computation and standard computation evaluated by effect magnitude |
| Source balance × lineage scarcity | Large hepatocyte-rich sources can dominate a pooled objective | Natural-frequency versus donor/source-balanced training, reporting each lineage |

### 18.3 Candidate search budgets and escalation

These are proposed configuration counts for planning, not job counts or authorized allocations. A configuration means a defined model/input/loss recipe; it can require several fits across folds and seeds. Complete an informative comparison before expanding it.

| Development family | Initial planning range | Expansion criterion |
|---|---:|---|
| Local-score calibration and simple heads | 12–24 configurations | Reproducible development gain or a clear unresolved calibration issue |
| Length, pooling and representation | 12–24 structured combinations, then selected layer/compression variants | Gains survive matched-example and coordinate checks |
| Adapter placement/rank/depth | 24–48 sampled or structured configurations | Stable training and evidence about which update family is useful |
| Partial unfreezing | 8–16 configurations | Pilot gradients and memory support training; final scale follows learning curves |
| Full tuning | 2–6 carefully selected recipes initially | Adequate data and measured compute; compare with simpler alternatives at the same target |
| Context/fusion architectures | 12–24 configurations per sufficiently supported endpoint | Gain beyond both context-only and sequence-only baselines |
| Multitask/transfer models | 8–16 configurations per assay pair | Incremental benefit or a specific negative-transfer explanation worth testing |
| Event/haplotype specialist models | A small balanced family first; sized to independent labels | Labels, events/phase and controls support the proposed complexity |

Use structured random or space-filling sampling for wide continuous/log-scale ranges. A Bayesian search or successive-halving method is an optional implementation choice, not a new mandatory dependency. Early pruning uses development loss and numerical failures, not protected data. Compare enough training exposure to avoid automatically discarding a slower-learning model.

After finalists are selected, use the established five seeds `1103, 2909, 4721, 6673, 8111` and bootstrap seed `20260914` where compatible with the existing comparison. Record any new random seed before sampling. Report the ensemble and individual-seed distributions separately.

### 18.4 What qualifies as progress

| Result | Interpretation and next action |
|---|---|
| Better native error with unchanged correlation | Potentially useful calibration; report magnitude improvement rather than claiming a ranking advance |
| Better ranking with poor calibration | Useful prioritization candidate; perform independent calibration before reporting effect magnitudes |
| Gain over HyenaDNA but loss to local AlphaGenome scores | A representation result, but no demonstrated improvement over the strongest current predictor |
| Adapter gain, with broad degradation of native outputs | Test whether the task-specific benefit justifies a specialist; do not claim a generally improved backbone |
| Gain confined to one source or strong-effect subset | Limited supported domain; test the competing source/ascertainment explanation |
| No gain with broad intervals | Insufficient precision for the chosen improvement; use T23 rather than declaring equivalence |
| No gain with a narrow interval excluding the planned useful effect | An informative limit for this approach and population; retain it and move to a materially different information source or hypothesis |
| Better mean profiles but no donor/allele improvement | Completion of a profile task only; revise the objective if donor/allele differences remain the goal |
| Successful state prediction only from source/technical variables | No transferable molecular-state result established; repair the design or obtain overlapping states/sources |
| New lineage head without allele labels | Profile capability established to the extent tested; allele-effect capability remains untested |

For confirmatory claims, retain the endpoint family fixed before evaluation and apply Holm at family-wise 0.05. Exploratory mechanism/variant families use BH within named complete families. Reporting BH across configurations does not correct all the bias introduced by choosing a model on a repeatedly inspected test; protected evaluation and nested selection remain necessary.

### 18.5 Assumptions to challenge rather than inherit

- Longer input is not automatically better; eligibility, target geometry and useful information all change with length.
- Native-score calibration may outperform a large embedding head. Preserve it in the comparison.
- A poor frozen representation does not establish that adapters or broader tuning cannot work.
- Missing pretrained adult outputs motivate new measured heads; missing measured labels constrain the claim, not the existence of a possible architecture.
- More loci and more donors answer different information questions. A significance-selected reliability estimate does not determine which investment is best.
- A condition-aware model must predict an allele/state difference, not simply identify treatment, donor or cohort.
- Assays can share useful predictive information without having identical effect signs, scales or molecular mechanisms.
- Human haplotypes, synthetic combinations and population disease interactions are different objects and require different measurements.

### 18.6 Which uncertainty each variable resolves

| Variable family | Scientific uncertainty | Required paired control | Campaign |
|---|---|---|---|
| Length × pooling × adaptation | Does the gain come from additional context, preserved local information or trainable sequence features? | Same eligible examples, target geometry and native baseline at each length | A, C |
| Encoder versus transformer/decoder | Does the task need local motif information or long-range/event structure? | Competitive encoder-only model and comparable head/optimization opportunity | A, C, E |
| Context representation × fusion × objective | Is a context gain biological and sequence-dependent, or reconstruction of an average profile? | Sequence-only, RNA/context-only, additive and shuffled-context models | B |
| Profile, donor, allele and interaction loss weights | Which supervised quantity is actually improving? | Separate native metrics for every target; loss normalization fitted in training | A–C |
| Training-data size × trainable depth | Is the current limit independent information or representational capacity? | Nested data subsets; fixed evaluation; controlled compression rank | A, B, F |
| Assay sharing × context/target match | Does transferred information concern the same molecular relationship? | Receiving-assay specialist; matched and mismatched transfer; real overlap census | C |
| Event definition × label/source exposure | Is an apparent splice gain a better estimand or selection on known outcomes? | Retained E10 baseline; fixed competing-event rules; independent event assessment | C, E |
| Genomic background × measured state | Does actual cis variation alter the context prediction? | Same focal variant and measured context; dosage and reference-background baselines | E |
| Numerical precision × hardware × weak-effect magnitude | Are small allele contrasts larger than implementation variation? | Same sequence/readout on a fixed device class; higher-precision subset | A–F |
| Source balance × lineage coverage | Does pooled training hide poor transport or scarce-lineage labels? | Donor/source-balanced comparisons and lineage-specific coverage/error | B, F |
| Saturation scale × abundance × region selection | Is apparent regulatory grammar explained by active-region properties? | Native-scale analysis, matched activity/variability/GC/TSS/width, full eligible population | D |
| Compression/distillation × uncertainty | Does a cheaper public model retain the actual scientific capability? | Full-model and native baselines; calibration and unsupported-context behavior | F, T24 |
| Original method × MASLD adaptation × sequence evidence | Does adaptation add useful disease knowledge, and which information explains it? | Original method and simple modality-native comparator; same eligible population, independent biological endpoint and fixed programs | B, D, E, F |
| Label definition × LD × uncertainty | Does a larger association set improve molecular prediction or teach correlated marginal signals? | Separate native heads; matched LD/scales where required; independent measured allele endpoint | A, C, F |
| Consensus-to-lineage curriculum × trainable depth | Does specificity training retain broad grammar and improve scarce-lineage outputs? | Same total exposure; training-defined specific regions; full-universe and allele evaluations | A, B |
| Link model × allele-effect features | Which stage improves target discrimination? | Strong linker alone, sequence alone and calibrated combination; power-qualified labels and matched coverage | C, E |

This matrix and the T01–T24 menus define the search space. They are not a Cartesian grid or authorization to run every combination. First compare identifiable model families, then explore interactions within development data. Retain precise negative results; change the information source or hypothesis when a useful gain is excluded.

## 19. Campaign sequence and deliverables

The campaigns organize T01–T24 around scientific outcomes. They share data and comparisons; they do not create six separate models by default. Campaign F's evaluation design and T24's database specification begin alongside data reconciliation. A task can support more than one campaign, with one owner and one record of its result.

| Campaign | Primary tasks | Shared dependencies | Scientific deliverable |
|---|---|---|---|
| A. Transferable liver allele effects | T01–T09 | T20–T21, T23 | Best supported signed accessibility predictor |
| B. Cellular context | T12–T15 | T01–T02, T08–T09, T20–T21, T23 | Separately validated profile, donor-difference and allele-interaction capabilities |
| C. Assay transfer | T10–T11, T16–T17 | Eligible A/B recipes, T01, T20–T21 | Directed map of transferable molecular information |
| D. Sequence sensitivity and remodeling | T22 | T13, T16, existing P4, T20–T21 | Supported connections and separations between inherited regulation and disease state |
| E. Indels, targets and haplotypes | T16–T19 | T01–T02, T17, T20–T21, T23 | Event- and haplotype-aware molecular alternatives |
| F. Generalization and public use | T20–T24 | Selected recipes from A–E | Independent benchmark, uncertainty, public database and permitted models |

### Campaign A. Establish the strongest transferable liver allele-effect predictor

**Goal:** improve signed endogenous accessibility prediction beyond the strongest eligible existing model, then determine which information provides the gain.

**Supporting evidence:** E03–E04, E20, E24–E25 and E32–E34. Local zero-shot macro Spearman 0.7441 on 3,627 matched variants is a strong baseline for that population, not an estimate for the full 32,322-variant training universe.

**Design and variables:** native scores, calibration, embedding heads and residual corrections; local and long inputs; variant-bin, peak, target and multiscale pooling; encoder-only tuning, adapters, partial and full tuning; profile versus effect supervision. Cross length × pooling × adaptation and independent training loci × trainable depth. Retain lead-selected and broader eligible variant–region populations separately, with ascertainment and label uncertainty stated.

**Decisive comparison:** Currin caQTL and the two allelic cohorts retain distinct estimands. Compare AlphaGenome, ChromBPNet, Borzoi, existing heads and simple native-assay models on matched examples, splits, coordinates and eligible exposures. Different-length adaptation needs a same-length native baseline as well as the historical 1-Mb reference. Reserve GSE277774 as a candidate external source subject to section 5.1's overlap and field checks.

**Success and output:** an independently evaluated gain over the strongest eligible comparator, with locus-respecting uncertainty, native-scale calibration and coverage. **+0.02 absolute Spearman is a stretch ambition**, not an expected result or a substitute for prespecifying an endpoint-specific useful effect. Deliver the best supported model and complete comparison even if no adaptation wins.

**Paper consequence:** model comparison in Figure 7; supported molecular interpretation in existing genetic panels. Source admission and T20 design precede final evaluation; full saturation is not a dependency.

### Campaign B. Learn cellular context without losing allele sensitivity

**Goal:** predict reproducible donor and cellular-state differences beyond both sequence-only and context-only models.

**Supporting evidence:** E05, E20, E26 and E36–E38. Begin treatment-interaction development in LX-2. HepG2 remains measurement-limited. The tested Corgi configuration learned almost no donor-differential H3K27ac signal despite an informative RNA-only comparator.

**Design and variables:** test section 8's sequence + context + interaction decomposition with lineage labels, regulator RNA, broader RNA, fixed program scores and admissible paired modalities. Compare additive models, bilinear interactions, FiLM and attention; shared versus assay-specific components; profile-only, donor-difference, allele-effect and combined objectives. Vary loss weights, donor/source balance, context masking and missing-modality training. MIRA-like paired inputs and existing context-model adaptations enter through section 2.1's modality and evidence requirements.

**Decisive comparison:** cross held donors with held regions against sequence-only, RNA/context-only and additive baselines. Exclude held targets and relevant genomic features from context construction; learn representations within training partitions. Shuffled context must preserve relevant source/batch strata. Technical-covariate, composition and source controls test whether context supplies biological information. Evaluate cohort transport separately.

**Success and output:** improvement in the explicitly named donor or interaction endpoint on held donors and regions, beyond both strong single-input baselines and the controls. Report profile accuracy and allele-interaction accuracy independently. An accurate mean profile alone completes neither the donor-difference nor the interaction goal.

**Paper consequence:** supported context comparisons in Figure 7 and the existing cellular-context figures. Sparse-lineage label expansion can proceed alongside profile development; no allele claim precedes measured allele evaluation.

### Campaign C. Explain what regulatory information transfers between assays

**Goal:** establish when reporter activity, endogenous accessibility, expression and RNA processing share predictable information.

**Supporting evidence:** E02, E09–E11, E20 and E27–E28. Preserve the failed reporter extra-column transfer, weak RNA direction and completed signed-splice limitation. Reporter variant-bin pooling is promising secondary development evidence requiring independent evaluation.

**Design and variables:** prioritize peak-proximal reporter–caQTL pairs, reliable LX-2 contrasts and exact target/event identity. Compare shared encoders with assay-native heads, separate specialists, joint versus sequential training and context-matched versus mismatched transfer. Expand public functional-assay sources only after documenting overlap and foundation-model exposure. Distinguish proximity from an experimentally supported molecular link.

**Decisive comparison:** adding source-assay supervision must improve the receiving assay on independent loci/examples over its otherwise matched specialist. Exclude receiving test loci from every auxiliary head. For splicing, test event identity, magnitude and direction separately. E10 is complete: retain its weak direction below the ranking bar; a fixed, better-matched event denominator is a new hypothesis, not a reinterpretation of the old primary result.

**Success and output:** a directed assay-transfer map with the receiving endpoint, native units, eligible context, precision and uncertainty. Improvements confined to the source assay do not qualify. Reliable negative transfer also defines where a shared model should give way to specialists.

**Paper consequence:** Figure 7 transfer comparisons and molecular interpretation within existing figures. T01's overlap census and assay-native labels precede training; no downstream signed-splice ranking claim assumes success.

### Campaign D. Connect regulatory sequence sensitivity to MASLD remodeling

**Goal:** identify reproducible relationships linking inherited regulatory variation to disease-associated molecular states, including where the two remain separate.

**Supporting evidence:** E05–E08, E12, E16–E18, E23 and E37–E40. Full P4 remains unfinished in the inspected outputs; bounded B3 and feature pilots cannot stand in for its registered universe or F1–F5 tests.

**Design and variables:** complete and reuse the existing saturation universe before dependent full-universe analysis. Test sequence sensitivity against donor regulation, local versus distal links and measured buffering/amplification across contexts. Keep direct-MASLD and liver-enzyme signals separate. Project fixed117 unchanged; retain cross-assay module v2 results as a separate analysis. Compare eligible SCENIC+, sc-linker or regulator-response adaptations with their original method and simple alternatives when public evidence supports the added question.

**Decisive comparison:** sequence sensitivity must add information beyond promoter proximity, signal abundance, variability, GC, region width and regulatory annotation. Use matched active-region and co-expressed-set controls. Protein and spatial corroboration require compatible biological units and assay coverage; count donors, not spots. Observational RNA–chromatin coupling is not evidence of regulatory mediation.

**Success and output:** a reproducible relationship survives the relevant competing explanations and clarifies a locus or disease program. Deliver effect sizes, full tested families, uncertainty, and unresolved alternatives. For sparse genetic overlap or narrow assay support, report the limit rather than expand the selected subset until a positive result appears.

**Paper consequence:** stronger biological interpretation in the relevant existing figures and Gene Catalog; model-dependent findings also require Campaign F. No new clustering, fixed117 optimization or reopening of closed integration is implied.

### Campaign E. Resolve indels, targets and haplotypes

**Goal:** expand SNV annotations into independently assessable molecular alternatives.

**Supporting evidence:** E09–E10, E13–E16, E19 and E29–E30. Reuse the 156 already-scored indels and corrected orientation, with the final paired analysis of 142 variants in 97 blocks. Coverage recovery and effect accuracy remain different outcomes.

**Design and variables:** compare SNVs and indels on compatible native accessibility, expression and exact event readouts; reference backgrounds versus actual public phased haplotypes; candidate targets from proximity, existing genetic evidence and eligible regulatory-network links. Compare additive and interaction models on identical readout geometry and matched controls. Preserve unresolvable phase and orientation; synthetic combinations are hypothetical sequences.

**Decisive comparison:** independently measured molecular effects or targets distinguish alternatives beyond proximity and existing genetic evidence. Haplotype claims need measured combinations or matched natural backgrounds with identifiable contrasts. Predicted nonadditivity alone does not establish biological epistasis.

**Success and output:** improved posterior coverage plus independently supported effect or target discrimination where labels permit. Every event/haplotype record retains alternate mechanisms and names the observation needed to distinguish them. If compatible public labels are absent, complete the coverage/identifiability result and stop short of an accuracy claim.

**Paper consequence:** locus interpretation and queryable Gene Catalog alternatives. T01 identity checks precede inference; source admission and T23 determine whether the richer validation can run without new experiments.

### Campaign F. Demonstrate methodological and scientific generalization

**Goal:** determine which capabilities extend beyond one cohort, assay or convenient subset and make the supported results reusable.

**Design and variables:** evaluate new regions/loci, donors, batches, studies, cell types, measured states, variant classes and effect magnitudes. Assess ancestry-specific accuracy only where actual validation labels permit it. Challenge missing genes/modalities and unfamiliar contexts; report support and abstention behavior. Preserve GSE289173 outcome restrictions and assess its eligibility through the existing evaluator without exposing protected outcomes to development.

**Decisive comparison:** strongest eligible AlphaGenome and specialist baselines, simple task-native models, and the original version of any cross-disease method being adapted. Separate foundation-model exposure, adaptation exposure and evaluation inspection. Region holdout within an exposed study is not independent-study validation.

**Success and output:** the stated improvement or biological relationship survives an independent test with appropriate exposure accounting. If none is available, return an exposed benchmark and explicitly unresolved generalization. Deliver capability-specific uncertainty, failure domains and matched performance; T24 makes these part of the public resource.

**Paper consequence:** Figure 7 evidence and Gene Catalog integration. Novel MASLD adaptation can qualify under section 2.1; established context conditioning by itself does not justify claiming invention.

### 19.1 First bounded assignments to make ready for execution

| Assignment | Clear goal | Concrete design and completion criterion |
|---|---|---|
| Evidence and current-work reconciliation | Avoid repeating completed or active science | Resolve E01–E40 against tables/code; check final P4 and local adaptation outputs and authorized live jobs; return completed, unfinished and conflicting records with dispositions |
| Public-source admission, T01/T20/T23 | Identify usable independent information | Audit GSE277774, protected GSE289173 eligibility and public assay/phase sources without development access to protected labels; deliver participant/locus overlap, schema, exposure and reuse decisions |
| Native-score and representation comparison, T03–T05 | Establish the strongest baseline and what pooling retains | Same caQTL folds; calibration and existing heads; two lengths × variant/local/target pooling; native error and correlation; measured storage |
| Adaptation feasibility and depth, T06–T09 | Test whether updating sequence features adds allele information | Verify a local training route and effect objective; measure backward-pass memory; cross selected pooling with adapter depths/ranks, then bounded partial/full tuning justified by data curves |
| Context and interaction design, T12–T15 | Separate useful context from average-profile prediction | LX-2 interaction and explicit donor-difference objectives; matched RNA-only, sequence-only and additive controls; crossed donor/region and context-exclusion specification |
| Joint-assay and event design, T10–T11/T16–T17 | Test a transfer route and event definition not already rejected | Exact overlap census; receiving-assay primary endpoint; retain E10 and extra-column negatives; specify independent assessment before new fitting |
| Full-saturation reconciliation, T22 | Reuse P4 and define its biological interpretation | Identify complete retrieval/features/tests or exact missing outputs; inspect original F1–F5 results before specifying any dependent amendment |
| Cross-disease adaptation review, T13/T16/T22 | Select a method that adds interpretable MASLD knowledge | Source-method and prior-MASLD search; native modality assumptions; public endpoint; original-method and simple comparator; admit a bounded comparison or record why it is not yet testable |
| Database specimen, T24 | Make a queryable scientific deliverable concrete | Define linked records and representative queries in section 20; document public-core versus restricted-model inference; completion requires traceable evidence and honest unsupported answers |
| Independent verification, T20–T21 | Confirm that a claimed gain answers the specified question | Reconcile biological n, estimand, coordinates, signedness, overlap and native units; verify the primary contrast and uncertainty from deposited outputs; return disagreements to the coordinator |

These are assignments for the subsequently authorized execution, not submissions made by this document. Design dependencies can proceed in parallel. A final protected evaluation depends on a fixed recipe and complete exposure accounting; source admission and evaluation design do not wait for a trained winner.

### 19.2 Resource planning for execution

Each campaign needs a measured pilot and a concrete allocation before submission. Estimate storage as well as compute. For a frozen tensor with `N` variants, `A` allele/orientation forms, `P` retained positions and `D` channels, the uncompressed float32 size is `4 × N × A × P × D` bytes. Retaining 128-bp bins from a full 1-Mb window for 32,322 variants, four allele/orientation forms and 3,072 channels would be about **13.0 TB decimal**, whereas one pooled vector per form is about **1.59 GB decimal**. These are [Estimate] raw-array sizes, excluding metadata, temporary arrays and compression; choose retained positions deliberately.

Frozen representations can be cached; trainable-backbone features must be regenerated during training. Estimate memory for parameters, gradients, optimizer state and saved activations before deciding the batch size. Checkpoints, accumulation and distributed execution are implementation options, not evidence that the scientific comparison will improve. Numerical precision and cross-process reproducibility should be assessed on the small allele differences being interpreted.

The execution proposal should list configurations, folds, seeds, expected fits, GPUs/CPUs, RAM, wall time, storage and aggregate allocation. Follow the project's compute-node, partition, QOS and submission rules. No installation, job submission, new sample acquisition or experimental execution is part of the present document-writing task.

### 19.3 Completion standard for the expanded program

The program reaches its scientific objective when the evidence supports a useful, independently assessable capability: better native allele effects, measured context dependence, cross-assay transfer, target/event resolution or real haplotype response, with competitive baselines and calibrated scope. It need not win every task or use one architecture for every assay.

Every completed task should return a concise scientific record: its goal, exact biological population, native estimand, inputs, variable settings, comparator, split/exposure design, uncertainty and multiplicity method, results including failures, compute cost, and the next decision. The broad task menu remains available for materially different experiments; it is not permission to keep changing the same test until a positive result appears.

### 19.4 Authorized initial week, 2026-09-15

The implementation instruction supersedes this document's earlier documentation-only
execution statements within this bounded campaign. It preserves T01–T24, all six
campaigns and previous negative results. Week one is an initial checkpoint;
it is not completion of every scientific aim.

| Limit | Authorized amount and accounting |
|---|---|
| GPUs | Initial four-GPU / 300-GPU-hour envelope superseded by the explicit 2026-09-16 user instruction: no campaign compute or concurrent-GPU cap |
| CPUs and RAM | Initial 1,500-core-hour envelope superseded on 2026-09-16; allocate as scientifically needed within scheduler limits |
| New files | 500 GB decimal, including temporary files and checkpoints |
| Existing retrieval | P4 retrieval and its existing continuation jobs are accounted separately |
| Submissions | Individual `sbatch` jobs, `nslab` account/QOS; no arrays or interactive-QOS substitution |
| Partitions | CPU, bigmem, I/O and GPU authorized as needed by the workload; obey each partition's scheduler limits |
| Acquisition | Suitable targeted public sources; no controlled-access data or unbounded sequencing downloads |
| Protected outcomes | GSE277774/Wenz and GSE289173 remain closed to development during week one |

The executable journal is
`GWAS/finemapping/results/alphagenome_campaign/week1-20260915/budget.json`.
[`campaign_budget.py`](../../scripts/analysis/alphagenome_campaign/campaign_budget.py)
reserves requested wall time before submission and replaces it with measured
allocation only when a job is terminal. It records total CPU/GPU hours and checks
requested storage and current user memory. The 2026-09-16 compute amendment
preserves the earlier envelope in the journal's authorization history; the
500-GB file limit remains. P4 reduction, source processing and independent
scientific checks remain required work alongside model fitting.

**Initial comparisons.** Native 1-Mb scoring covers the complete existing
32,322-variant, historically selected Currin population, with all exclusions
reported. A separate all-pair reporter-position population avoids new minimum-P
target selection. Twelve frozen-head recipes cross 2,048/16,384 bp,
variant/symmetric/measured-target pooling and ridge/shared scalar heads.
Eight adapter configurations cross ranks 4/16, final three/five attention blocks
and variant/symmetric pooling. Two query/value partial-unfreezing pilots use
backbone learning rates 3e-6 and 1e-5. Head/adapter learning rate is 3e-4,
AdamW weight decay 1e-4 and gradient clipping 1. A measured 100-step pilot
precedes larger fits; initial single-split, single-seed results are feasibility
evidence. Formal finalists require the existing five seeds and nested selection.

**Context endpoints.** LX-2 allele-by-treatment comparisons preserve experimental
pairing and the operative long-range block folds. Donor comparisons use the
99 paired GSE267145 participants, training-only calibration, continuous RNA
inputs and counted-interval semantics with the recorded 1-bp residual. Held
donors at trained regions and held donors at held regions are different tests.
Local-gene, ±1-Mb and chromosome exclusions have dimension-matched random,
technical and shuffled controls. Without authoritative target links, overlap-gene
exclusion is a sensitivity; exact target removal remains untestable.

**Interpretation and outputs.** P4 retains historical groups and adds adult,
developmental and unresolved-life-stage groups as an amendment informed by
inspected results. Full original F1–F5 analyses require actual U1/U2 completeness
and their original minimum-group rules. Hosted saturation outputs are not
training features. Reporter transfer first distinguishes physical overlap from
nearest-peak relationships. Exact splice-event and natural-haplotype validation
require source-native event/phase identities and molecular joins. The local
SQLite specimen exposes variant, gene, region, event, program and context records
with explicit unsupported answers. No individual MASLD-risk score is supplied.

Numerical candidate outputs remain under the campaign directory; candidate
Figure 7 panels are in `figures/candidates/regulatory-week1-20260915-v4/`. Compute checks,
failures, source suitability and task dispositions are recorded there. Current
SLURM state must be queried directly. Figures and interpretations can be proposed
from reviewed results; final adoption and publication remain separate user decisions.

#### Initial implementation checkpoint, updated 2026-09-16

The following is the initial research checkpoint, not a declaration that the
one-week scientific campaign or T01–T24 is complete. Exact completed context
results and their uncertainty are in
[RESULTS.md](../RESULTS.md#initial-regulatory-context-comparisons-2026-09-15-not-adopted).

**Overall progress:** the initial model comparisons and local database specimen
are substantially implemented. The program's central claims remain open:
adaptation has not established improvement over the strongest native predictor;
context has not established improvement across held donors and held regions;
and new target, event or haplotype mechanisms have not been independently
validated. Completed fits and numerical checks are intermediate outputs, not
completion of those scientific goals. Section 22's sixteen added comparisons
remain an agenda; the initial week does not execute all of them.

| Campaign | Implemented comparison or artifact | Remaining scientific requirement |
|---|---|---|
| A: allele effects | Twelve fixed-representation recipes on 28,600 matched variants; twelve 5,000-step adaptation/control arms on 6,834 validation variants; twenty uncertainty-aware fits complete and independently checked. Six adapters improve their matching fixed controls, but the newly completed matched native table puts calibrated native ATAC ahead of the best adapter on 6,441 shared validation variants (RMSE 0.536015 versus 0.545425). This latest numerical comparison awaits independent reconstruction. Native coverage is 30,611 scored / 32,322 attempted. Zero-adapter checks and measured regional ATAC assessment are complete; both partial-update pilots increase regional ATAC error | Independently reconstruct the complete matched table before selecting another comparison. No native-superiority claim for adaptation is established. Regional ATAC is a same-source assessment. Conditional single-seed uncertainty does not replace nested selection, five seeds or independent-source evaluation |
| B: context | Replicate-paired treatment heads on 4,358/4,359 constructs and the stronger 99-participant/128-region donor comparison are complete and independently checked. No treatment fit beats zero by point RMSE. RNA improves held-donor prediction at trained regions (RMSE 0.352445 versus training means 0.382406); pairwise FiLM is worse than the selected RNA baseline. Crossed held-donor/held-region context benefit remains unestablished. Context exclusions and matched controls are complete for the available annotations | Adult lineage heads, new-source transport and a useful held-region context gain remain unfinished. Validated target links and batch annotations are still missing. Keep unrestricted 42,163 RNA features separate from the matched 42,126 annotated-gene exclusion universe. Archived and new reporter representations differ slightly |
| C: transfer | Outcome-independent target-aware shared/private comparison on 2,765 pairs; 112 physically overlapping variants in 63 blocks. Original zero-selection result retained; learned-only development amendment establishes no improvement across the fixed 12-stratum family | No biological information-limit claim: intervals condition on source effects and fitted predictions. A different measurement source or hypothesis is needed before repeating this design |
| D: remodeling interpretation | Life-stage amendment, 100/100 archive pilot, full-reduction implementation and four focused checks; source-native F4 count covariates prepared | Full U1/U2 retrieval and reduction. Original F5 gene-set null needs the unscored eligible promoter background; do not restrict that background to observed program peaks |
| E: events and haplotypes | Exact source identities resolve 630 historical variant-event pairs; 621 have scores after removing 16 identical extra copies. Independently checked direction remains weak and below original ranking criteria. Actual local single-SNV and synthetic phased-SNV inference are complete. Independent reconstruction confirms the author Currin 138-column GSE26105 accession map and sparse read co-coverage across six selected donor-by-pair cases (three molecules in Z08, two in D17, zero in the other four) | Complete source-cluster denominators, raw junction union and new reference-allele checks. The author accession map and sparse read counts do not independently verify genotypes or retained phase. No final phased-call/phase-quality asset was listed in the scoped Currin public deposits; natural-background validation remains unresolved |
| F: independent evaluation | Wenz metadata maps 189 public assays; ATAC uses GRCh38 but preprint genotyping Methods names a GRCh37 configuration. RASQUAL algebra supports Pi > 0.5 as ALT-increasing; Pi is not a Currin additive beta. Wenz/Dudek share a source family; protected outcomes remain closed | Final association-VCF reference/conversion record, included 175 donor-to-assay/genotype crosswalk, checkpoint-bound model exposure, useful-effect margins and an isolated evaluator. Collection chronology supports nonoverlap with Currin only as an inference |

The checked Catalog extension exposes six entry points and 436,838 aggregate
evidence records: the original 61,316, 32,322 measured Currin associations, and
343,200 matched out-of-fold development predictions from all twelve recipes.
Individual prediction uncertainty is unavailable and measurement SEs are not
reused for predictions; fitted heads are not deployed. Unknown event identities, unavailable assays,
changed allele orientation, competing targets and unsupported joint effects
remain explicit. Program membership is not a measured program effect.
The initial figures illustrate these candidate results and limitations;
they do not change Resource memberships or evidence classes.

#### T01–T24 progress and next scientific decisions, 2026-09-16

“Pilot complete” below describes the bounded comparison, not the full task's
success criterion. Earlier negative results remain part of every task.
The completed native report is `model/final_21775091/` and its precision strata
are `review/precision-21775092/` under the campaign directory. Their job completion
does not substitute for independent numerical reconstruction. Reviewed context,
adaptation and loss comparisons are recorded in RESULTS and the campaign's
`initial_assignment_disposition.tsv`.

| Task | Progress against the scientific question | Next requirement |
|---|---|---|
| T01: labels and units | Endpoint tables, identities, exclusions and cross-fold sequence overlap checked. Selected Currin leads remain explicitly selected; separate all-pair reporter-position labels exist | Broader outcome-independent accessibility labels, masked allelic identities and exact event denominators |
| T02: strongest baselines | Native scoring, matched model tables and precision strata complete; latest numerical reconstruction pending | Confirm matched results, exposure and calibration; retain the strongest predictor even if adaptation loses |
| T03: sequence context | Matched 2,048/16,384-bp fixed representations and native 1-Mb scoring complete | Same-length native comparisons for trainable models; justified longer-context adaptation and construct/native-background comparisons |
| T04: representation | Twelve recipes cover two lengths, three pooling choices and two head types | Broader layer/base-resolution searches only for a specific unresolved endpoint |
| T05: fixed heads and residuals | Sixty fold fits independently reconstructed; training-only native calibration complete | Native-score residual learning and reusable fitted head weights unfinished |
| T06: adapters | Eight rank/depth/pooling pilots complete; six improve matching fixed controls. Gradient, allele-symmetry, zero-adapter and reload checks passed | Gain over native scores, then nested five-seed confirmation of an eligible finalist |
| T07: broader tuning | Two query/value partial-update pilots complete; neither improves the matched fixed control, and both impair measured regional ATAC | No broad/full-backbone tuning run; expansion needs a different justified design |
| T08: assay-native loss | Twenty source-SE/constant-SE Gaussian fits checked against ten squared-error fits. Lowest observed allele RMSE remains the squared-error adapter | Count/allelic/event likelihoods require their corresponding native measurements |
| T09: optimization | Stable scalar trainer, measured throughput and 5,000-step pilot schedule implemented | Five seeds, nested selection and justified schedule/learning-curve comparisons not run |
| T10: shared assay learning | Target-aware shared/private development comparison complete on 2,765 pairs | No receiving-assay gain established; broaden only with different compatible information |
| T11: reporter transfer | Physical-overlap/distance/target strata evaluated; fixed 12-stratum family retains unresolved and untestable comparisons | New assay/construct information or a different hypothesis; no general biological transfer-limit claim |
| T12: adult lineage heads | Existing donor-aware label coverage and life-stage distinctions retained | New lineage profile heads and expanded lineage allele measurements not run |
| T13: context inputs | RNA, technical, shuffled, genomic exclusions and dimension-matched controls evaluated; stronger nonlinear/kernel models checked | Additional paired chromatin/program inputs and usable target/batch annotations |
| T14: treatment interactions | Full-source LX-2/HepG2 comparison complete; none of three fitted routes beats zero by point RMSE | A materially different input or measurement-supported design; no patient-response claim |
| T15: donor differences | Stronger comparison and 4,920 individual neural fits complete. RNA helps at trained regions; context benefit at held regions unestablished | Improve or precisely bound crossed donor/region prediction; independent-cohort transport unfinished |
| T16: target genes | Existing target identities and conflicts preserved in tables and queries | Independent enhancer–gene benchmark admission and target-discrimination modeling largely unstarted |
| T17: events/indels | Historical event identities reconciled and duplicates removed; weak signed-splice result reconstructed; earlier indel coverage retained | Exact source clusters, gained/lost raw junction union and new event-aligned inference missing |
| T18: natural backgrounds | Currin donor/accession join and six read co-coverage cases checked; synthetic phased-SNV inference works | Verified donor genotypes, usable phase and molecular joins; no natural-background accuracy result |
| T19: haplotype interactions | Corrected GNMT null preserved; sparse read support documented | Measured joint effects with adequate phase/contrasts unavailable; new biological interaction test not run |
| T20: generalization | Development folds/exposure records and Wenz metadata inspected; protected outcomes closed | Final source eligibility, fixed recipes/margins and isolated external evaluation |
| T21: uncertainty | Training-only calibration, paired prediction uncertainty, complete exploratory families and chromosome sensitivities reported | Reconstruct latest native table; fitting/selection uncertainty, predictive interval coverage and useful-effect margins unfinished |
| T22: biological interpretation | P4 life-stage/completeness pilot, candidate figures and existing fixed117 links prepared | Full U1/U2 and F1–F5 accounting; eligible F5 promoter background. New network, multicellular, spatial and buffering comparisons largely unstarted |
| T23: useful new data | Source-suitability records and measured Currin profiles acquired and checked | Formal label learning curves and targeted acquisitions filling event/target/phase gaps |
| T24: reusable resource | Six-entry SQLite specimen, 436,838 evidence records and query examples checked; local native inference works | Selected fitted models, supported prediction uncertainty and full public-model evaluation/distribution unfinished |

**Next decisions in order:** independently reconstruct the latest matched native
comparison; retain native ATAC as the provisional allele baseline while deciding
whether a native-score residual or same-length comparison supplies a distinct
test; complete P4 and its explicit missing-background accounting; and admit the
specific target/event/phase measurements needed for the underdeveloped branches.
Five-seed expansion should follow a defensible comparison, not a better result
against a weaker fixed head alone. Protected evaluation remains closed until its
source and outcome-handling requirements are met.

## 20. The MASLD regulatory database and public models

### 20.1 Four linked scientific components

T24 makes the database a primary deliverable of this paper, integrated with the Gene Catalog. Reuse the existing Atlas substrate of 417,155 variants across 989 signals and the candidate program attachment; identify the new capability added by every extension. This is a research resource for molecular interpretation, not an individual clinical risk predictor.

| Component | Scientific objects | Required distinction |
|---|---|---|
| Measured evidence | Genetic signals/posteriors; marginal molecular QTLs; allelic imbalance; measured reporter or endogenous perturbation effects; cross-sectional remodeling; protein and spatial observations | Source, assay, population, biological n, measurement units and uncertainty remain attached to every observation; repeated imports of one measurement are not independent support |
| Predictions | Archived Atlas scores, local AlphaGenome outputs and adapted/specialist outputs | Separate model/version, sequence background, readout, signedness, raw versus quantile units, exposure and evaluated capability; predictions are never stored as observations |
| Regulatory relationships | Variant–region, variant–event, candidate target, haplotype, context and fixed-program links | State how each link was obtained and its alternatives; genetic colocalization and observational coupling do not establish mediation |
| Evaluation | Matched comparisons, uncertainty/calibration, tested populations, failures and unsupported domains | A model's profile, allele-effect, interaction and target capabilities have separate evidence; a global performance number cannot certify every query |

Start from normalized tables and a documented local query interface, then expose the same supported records through the existing Gene Catalog. Downloads should include compact tabular data, a data dictionary and reproducible example queries. Publication of the service and artifacts is a later authorized action. Do not create a competing universal gene rank or change Resource evidence classes.

### 20.2 Identity, evidence and uncertainty carried by each answer

| Record family | Minimum fields |
|---|---|
| Variant or haplotype | Genome build; original and normalized IDs/alleles; reference verification; normalization/orientation status; interval convention; phased alleles and phase source if known; sequence-background identity |
| Molecular object | Region coordinates and assay; gene/transcript version; exact junction/event and strand; target/readout window; aggregation rule and output units |
| Context | Tissue, adult/developmental status, lineage, measured condition and cohort; donor-level representation definition; available/missing modalities; observed versus synthetic context |
| Observation | Source accession/table and producer; biological population and n; association/imbalance/experimental effect type; effect allele; native effect and standard error or other supported uncertainty; covariates; exposure/reuse conditions |
| Prediction | Model/checkpoint and permitted derivative identity; input/representation recipe; reference/alternate values and defined contrast; signedness; calibration population and uncertainty method; unsupported status/reason |
| Genetic support | Signal/trait and fine-mapping version; variant posterior; served and missing posterior mass with denominator; alternate variants/targets; direct-MASLD versus liver-trait category |
| Relationship and evaluation | Link type/source; competing targets/events; evaluation endpoint and split; compatible tested context; matched baseline/performance/interval; known source/checkpoint exposure; limits |
| Reuse | Source-specific terms and access route; permitted table/code/weight distribution; whether local reproduction needs separately obtained weights or context; source version and retrieval date |

Keep unserved, untested, untestable, non-applicable and measured near-zero outcomes distinct. Do not turn posterior weight multiplied by a model score into a probability of mechanism. Do not sum heterogeneous QTL, reporter and chromatin units. Summaries of posterior coverage retain their exact signal-specific denominator and missing alternatives.

For combined evidence, identify shared donors, experiments, loci and source publications so correlated records are visible. Display conflicting observations rather than averaging away different molecular estimands. A changed fine-mapping release produces a separately identified genetic view until adoption is authorized.

### 20.3 Queries the resource must answer

| Entry point | Supported answer and acceptance example |
|---|---|
| Variant | Show measured/predicted assay effects, context and signal coverage. A normalized indel can be locally scored yet remain absent from the hosted Atlas; both facts must remain visible. |
| Gene | List competing variants, regions and events with target-assignment evidence. Nearest-gene assignment is labeled as proximity, and alternative targets remain accessible. |
| Region | Compare measured accessibility/chromatin regulation with sequence sensitivity in the stated coordinates. Partial saturation is never displayed as complete regional coverage. |
| Transcript/event | Show exact competing-event definitions, strand, detection and any supported magnitude/direction. An unsigned Atlas score cannot produce an up/down label. |
| Program | Display fixed117 projections and their source evidence; provide separately labeled cross-assay module results. Do not optimize membership to improve a queried disease association. |
| Context | Return only measured or evaluated context-specific capabilities. A model with profile validation but no allele-interaction validation says so; missing genes/modalities trigger the assessed fallback or an unsupported answer. |

Also test ambiguous build/alleles, unavailable target, unknown phase, unsupported ancestry/lineage, duplicate measurements and conflicting candidate targets. Each example must be reproducible from the distributed tables/code; use independently evaluated examples for demonstrations of accuracy and synthetic fixtures only for query behavior.

### 20.4 A reproducible public core and distinct research models

1. **Public measured/predicted tables:** release only permitted source-derived records with their source conditions. Use accession references or allowed summaries when redistribution is restricted. Do not expose donor-level genotypes or private context vectors merely because an aggregate result is public.
2. **Public inference components:** provide permissible code, preprocessing definitions, normalization parameters, calibration, heads/adapters or standalone specialists trained from eligible sources. A user obtains any separately licensed base weights through the authorized provider and runs a documented variant-to-output example locally. Record expected hardware, supported lengths/readouts and small-effect precision tolerance.
3. **Source-specific research models:** identify terms separately for source data, base weights and derivatives. The adopted GSE267145/GSE296875 decision allows derivative weights, including distilled students; it does not permit rehosting their raw or processed data. Any other source's actual restrictions still apply. AlphaGenome-derived weights remain subject to its noncommercial terms.
4. **Reproducibility demonstration:** on an eligible, non-protected example, a clean local run using the distributable components and documented base-weight route reproduces the supported prediction within a prespecified numerical tolerance. Compare full and compact models before substituting a cheaper model. Hosted-service availability must not be the only way to reproduce the public capability.

Source terms apply separately to data, code, base weights, adapted weights, outputs and deployment. T02 records the exact applicable terms before execution and T24 identifies which capabilities belong in the public core. If a source cannot support the public model, train a distinct eligible-source model and report its own performance; do not imply that a restricted model's accuracy transfers to it.

**T24 completion:** all four components join through stable scientific identities; the six entry points and unsupported cases work on a bounded specimen; every accuracy statement links to a compatible evaluation; the public inference example is reproducible; permitted artifacts and separate research-only capabilities are enumerated. Figures follow completed evidence: model comparisons in Figure 7, molecular interpretation in the relevant existing figures, and database integration with the Gene Catalog.

## 21. One orchestrator and bounded scientific assignments

### 21.1 Coordination responsibilities

The executing coordinator owns the central scientific question, T01–T24 dependencies, evidence reconciliation and integrated interpretation. The user authorized one coordinator and three nonoverlapping workers on 2026-09-15. Workers read the current authorities, reconcile active work, and return bounded analyses for independent review.

Three concurrent worker roles are useful when the work is independent: data/evidence reconciliation; model development or one specified biological analysis; and independent statistical/scientific verification. Roles may rotate across assignments. The coordinator retains manuscript direction, integration of results and the decision to propose an amendment. Workers return evidence, failures and unresolved issues; they do not alter paper claims, thresholds, fixed117 or protected outcomes.

Before new submissions, inspect current authorized jobs through SLURM and match their inputs/outputs to this plan. Dated notes are not live job status. Reuse completed hosted calls, local predictions and feature archives. In particular, recheck full P4, the completed signed-splice outputs and later local adaptation results before assigning duplicate work. No stale “weights unavailable” statement may block a now-supported local experiment; training feasibility still requires the relevant gradient, loss and memory checks.

### 21.2 Assignment contract

Every assignment contains the following fields. The coordinator supplies the scientific objective and evaluation rules; a worker should not have to invent them.

| Field | Required content |
|---|---|
| Identity and question | Campaign, T identifiers, E evidence IDs, one unresolved hypothesis, and how its answer changes the paper or resource |
| Inputs | Exact tables/code/model versions; native modality/units; source exposure and reuse status; biological population and n; allowed train/development/evaluation partitions |
| Boundaries | Allowed input paths and unique output directory/files; protected outcomes and forbidden joins; whether the work is documentation, analysis, training or evaluation; authorized compute limits |
| Comparison | Variable settings and interactions; primary comparator/endpoint; useful-effect ambition versus formal criterion; controls and exact evaluation unit |
| Statistics | Donor/locus/study separation; ascertainment and label uncertainty; precision calculation; seeds; resampling/model; named multiple-testing family |
| Dependencies | Required completed inputs, design decisions that can proceed now, and final assessments that must wait |
| Completion | Exact outputs, numerical/scientific checks, stop/expansion rules, negative-result handling and conditions requiring coordinator resolution |
| Return | Result table and source paths; numerator/denominator and exclusions; completed versus unfinished steps; compute used; interpretation supported by the data; unresolved alternatives |

Model workers see development data only. Preserve the existing GSE289173 evaluator restrictions: do not copy protected effect estimates into the coordinator or worker context, feature selection, prompts or shared development files. Subagents with the same filesystem are not a technical isolation mechanism. Use the existing evaluator's allowed data interface and outcome-handling policy; if it cannot preserve the required separation, retain the endpoint as unavailable for development rather than assuming a role label provides protection.

The coordinator controls edits to this plan and the human-facing authorities. Workers write to nonoverlapping assigned outputs and preserve named results. New public acquisition, dependencies, compute and publication follow the subsequently authorized scope and project rules. A research assignment does not authorize changing adopted results.

### 21.3 Search, verification and decisions

- **Design before scale:** source admission, native estimands, split/exposure design and memory measurements come first. T20/T21 design accompanies model development; their final evaluation follows recipe selection. Avoid circular dependencies that require a final model before its evaluation population can be defined.
- **Broad but interpretable search:** retain T01–T24 ranges and section 18 interactions. Reuse pooled features; compare promising families fairly before deeper searches. Extra configurations do not create independent evidence. Stop an approach when uncertainty excludes its prespecified useful gain, or return an explicit information deficit.
- **Cross-disease adaptation:** use section 2.1's admission record and assign the selected comparison within existing tasks. Test the original method, the MASLD adaptation and a simple alternative; report a useful disease application without inventing architectural novelty.
- **Independent scientific review:** verify allele orientation, coordinates, event identity, signedness, readout geometry, native units, source overlap, biological n and relevant null/uncertainty. Recompute high-impact summary claims from small deposited outputs when practical. Successful execution alone is not evidence that the estimand is valid.
- **Result synthesis:** return positive, negative, underpowered and untestable outcomes distinctly. A changed primary endpoint or selected subset is a new development design, never a silent replacement of a completed comparison. Update E dispositions with the later source and preserve historical experimental records.

### 21.4 Acceptance of this documentation revision

The revised document is ready for the future coordinator when:

- Every relevant prior analysis in the inspected Atlas, program, local-model and MASLD-benchmark records has a source, interpretation and next-action disposition; unresolved final outputs are named.
- T01–T24 retain their parameter menus and each specifies a goal, supporting evidence, unresolved question, data, variables, decisive comparator, evaluation unit, measurable success criterion, dependencies, outputs and consequence for the paper.
- Completed results, proposed experiments and unfinished dependencies are distinguishable; the signed-splice result and incomplete full P4 have the correct status.
- Public-data expansion, public reuse, the four database components and local reproduction have explicit implementation paths.
- Each campaign and cross-disease adaptation can be assigned with its comparison and completion rule without inventing a new scientific direction.
- Local paths/anchors, document structure and consistency with the human authorities pass inspection; applicable authority changes pass `python3 scripts/manuscript/validate_resource_scope.py`.

This revision is documentation work. Scientific reruns, downloads, training, job submissions, publication and changes to adopted results require their subsequent authorized execution. New experiments remain possible future follow-up rather than dependencies of this program.

## 22. Literature-informed expansion and experiments worth testing

**Research date:** 2026-09-15. This expansion follows a broad online search of primary research, official data resources and author-maintained code, followed by targeted reading of methods, limitations and data availability for the most relevant candidates. It covers sequence models, enhancer–gene mapping, personal genomes, cross-context statistics, RNA processing, multicellular regulation, perturbation benchmarks and public evidence databases. It is a curated research review, not an exhaustive systematic review or a priority claim.

**What changed after reading:** the strongest additions are more precise supervision, a consensus-to-lineage training curriculum, donor-contrast learning, stronger enhancer–gene comparators and explicit molecular-event models. Multicellular and spatial adaptations are valuable when they test a defined link to inherited evidence. More architectures, gene-network plots or outcome-selected subsets alone do not advance the goal.

### 22.1 Admission and priority

An added experiment is worth testing when it meets all four conditions: it addresses an unresolved result in E01–E40; it adds information or tests a materially different explanation; its required public measurements exist or have a bounded admission check; and a comparison can distinguish success from a plausible simpler explanation. The detailed designs below supplement T01–T24. They are not additional top-level tasks or authorization to run every option.

| Priority | Meaning | Initial candidates |
|---|---|---|
| First development comparisons | Address existing bottlenecks with eligible current data or modest feature reuse | Label/LD alignment; donor-contrast objectives; allele-focused adaptation and retention; consensus-to-lineage curriculum; baseline enhancer–gene comparison |
| Conditional extensions | Potentially important, but exact labels, overlap, modality or independent support must be admitted first | GET-like measured-chromatin context; allele-aware networks; construct-aware transfer; event/APA models; real diploid backgrounds |
| Biological synthesis after prerequisites | Depends on completed molecular comparisons or full P4 and suitable donor-level validation | Fixed-program multicellular coupling; spatially supported signaling hypotheses; distal regulatory buffering |
| Do not expand by default | No identified information gain or valid evaluation at present | An optimizer/adapter catalogue, unrestricted giant-model sweep, synthetic-donor augmentation, universal disease score, generic gene-network discovery, or new molecular-dynamics/protein-structure work without a specific supported MASLD question |

**Novelty correction:** SCENIC+ has already been applied to human MASLD in the Elison et al. preprint. The proposed advance must therefore be a stronger allele-aware or independently validated capability, not “first SCENIC+ in MASLD.” The preprint status and prior use are explicit in R29. Likewise, AlphaGenome already combines sequence predictions with enhancer-link models, and existing methods already condition predictions on cell context. The specific MASLD adaptation and its evidence determine the contribution.

A source read during this review is recorded as literature-exposed. Before using its data for evaluation, reconcile its publication family, donor overlap and model exposure. Search failure is not evidence that a MASLD application is novel. Preserve the protected GSE289173 evaluator and assess any possible cohort match without importing protected outcomes into development.

### 22.2 Annotated reading map

Each reference has a concrete use or a reason to limit an approach. Publication years refer to the linked record; preprints are identified. These papers motivate tests, not numerical performance expectations in MASLD.

#### Sequence, context and individual variation

| Ref | Primary reference / status | What to take from it and where it matters |
|---|---|---|
| R01 | [Avsec et al., Advancing regulatory variant effect prediction with AlphaGenome](https://www.nature.com/articles/s41586-025-10014-0), Nature, 2026 | Native multimodal and variant-effect baselines; resolution/context ablations; existing integration with enhancer-link prediction. Molecular prediction alone does not establish gene-to-disease effects. T02–T07, T16–T17. |
| R02 | [Patel et al., DART-Eval](https://arxiv.org/abs/2412.05430), 2024 preprint record; [author benchmark code](https://github.com/kundajelab/DART-Eval) | Distinguish regulatory classification, profile/activity prediction and variant counterfactuals. Strong task-trained baselines can challenge large DNA language models. T02, T20; avoid size-driven model selection. |
| R03 | [Kempynck et al., CREsted](https://www.nature.com/articles/s41592-026-03057-2), Nature Methods, 2026 | Consensus-peak training followed by cell-specific refinement and comparisons against larger backbones. Its cross-tissue and cancer-state applications motivate a liver-lineage curriculum; profile specificity still needs separate allele validation. T07, T09, T12. |
| R04 | [Fu et al., A foundation model of transcription across human cell types (GET)](https://www.nature.com/articles/s41586-024-08391-z), Nature, 2025 | Sequence and chromatin context provide a distinct route to expression and regulator relationships. Compare observed-chromatin input with RNA-only and sequence-only tasks without confusing imputation with allele response. T13–T16. |
| R05 | [Corgi](https://www.nature.com/articles/s41467-026-75527-2), Nature Communications, 2026 | Established sequence plus regulator-expression conditioning. The local E36 result motivates differential objectives and context-use diagnostics, not a claim that context conditioning itself is new. T13–T15. |
| R06 | [scooby](https://www.nature.com/articles/s41592-025-02854-5), Nature Methods, 2025 | A cell-state decoder is a relevant alternative to a single bulk output. Test adult-liver generalization and missing-context behavior with separate profile and allele endpoints. T12–T13, T20. |
| R07 | [Rastogi et al., Fine-tuning sequence-to-expression models on personal genome and transcriptome data](https://link.springer.com/article/10.1186/s13059-026-04091-1), Genome Biology, 2026 | Pairwise donor-difference objectives are directly relevant. Reported gains for new individuals at trained genes did not extend to held-out genes; cross donor and locus holdouts explicitly. T08, T15, T18, T20. |
| R08 | [Sasse et al., Benchmarking of deep neural networks for predicting personal gene expression from DNA sequence highlights shortcomings](https://www.nature.com/articles/s41588-023-01524-6), Nature Genetics, 2023 | Accurate cross-gene averages can coexist with poor personal-expression prediction. Evaluate donor differences and sign rather than using mean-profile correlation as the main success metric. T15, T18. |
| R09 | [AlphaGenome Enhances Personal Gene Expression Prediction but Retains Key Limitations](https://pubmed.ncbi.nlm.nih.gov/42039486/), preprint, DOI 10.1101/2025.08.05.668750 | A newer AlphaGenome-specific personal-expression assessment prevents extrapolating older Enformer failures unchanged. Use native AlphaGenome and genotype-based comparators; account for GTEx exposure. T02, T18, T20. |

#### Labels, targets and regulatory networks

| Ref | Primary reference / status | What to take from it and where it matters |
|---|---|---|
| R10 | [Urbut et al., Flexible statistical methods for estimating and testing effects in genomic studies with multiple conditions](https://www.nature.com/articles/s41588-018-0268-8), Nature Genetics, 2019 | Multivariate adaptive shrinkage learns effect-sharing patterns and uses effect uncertainty. It can regularize compatible context effects, with sampling covariance and training-only estimation. It does not make different assay units interchangeable. T01, T08, T14. |
| R11 | [van de Geijn et al., WASP](https://www.nature.com/articles/nmeth.3582), Nature Methods, 2015 | Allelic mapping bias and joint allele/total-count information matter for labels. Its original remapping approach excludes indel-overlapping reads, so SNV success cannot certify indel allelic labels. T01, T08, T18. |
| R12 | [Gschwind et al., An encyclopedia of human enhancer–gene regulatory interactions](https://www.nature.com/articles/s41586-026-10781-4), Nature, 2026 | ENCODE-rE2G supplies CRISPR-trained comparators, complementary evaluation types and a public interaction resource. Use the 2026 publication/version rather than silently mixing earlier preprint counts or models. T16, T20, T24. |
| R13 | [Sheth et al., Mapping enhancer–gene regulatory interactions from single-cell data (scE2G)](https://www.nature.com/articles/s41588-026-02695-8), Nature Genetics, 2026; [full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC13447104/) | Separate scATAC and paired-multiome models, with CRISPR-supervised links and within-cell-type information. These are strong task-specific comparators for adult liver target mapping; transport from training cell lines remains to be assessed. T12, T16. |
| R14 | [SCENIC+](https://www.nature.com/articles/s41592-023-01938-4), Nature Methods, 2023 | TF–region–gene networks can organize the interpretation of allele effects. Test the extra information from sequence evidence over the same unaugmented network. Prior MASLD use exists, R29. T16, T22. |
| R15 | [MIRA](https://www.nature.com/articles/s41592-022-01595-z), Nature Methods, 2022 | Paired RNA/ATAC representations distinguish local chromatin and broader cellular information. Use as an auxiliary representation comparison with native observation models, not a replacement for the Resource's integration. T13–T15. |
| R16 | [Jagadeesh et al., Identifying disease-critical cell types and cellular processes by integrating single-cell RNA-sequencing and human genetics (sc-linker)](https://www.nature.com/articles/s41588-022-01187-9), Nature Genetics, 2022 | Test whether improved SNP-to-gene maps clarify inherited support for fixed liver programs. Retain genome-wide GWAS/LD requirements and conditional annotation controls. T16, T22. |
| R30 | [Cell-type-specific and disease-associated expression quantitative trait loci in the human lung](https://www.nature.com/articles/s41588-024-01702-0), Nature Genetics, 2024 | A fibrosis-related example of separating cell-type sharing from disease-associated effects. Borrow the design and uncertainty treatment; lung effects are not liver training labels unless explicit cross-tissue transfer is being tested. T08, T14–T15, T22. |

#### Multicellular remodeling, functional assays and scientific limits

| Ref | Primary reference / status | What to take from it and where it matters |
|---|---|---|
| R17 | [Jerby-Arnon and Regev, DIALOGUE](https://www.nature.com/articles/s41587-022-01288-0), Nature Biotechnology, 2022 | Coordinated programs across cell types in disease motivate a donor-level test using the existing fixed117 scores. The proposed constrained analysis borrows this question; it does not rediscover or rename those programs. T22. |
| R18 | [Dimitrov et al., LIANA+](https://www.nature.com/articles/s41556-024-01469-w), Nature Cell Biology, 2024 | Integrates intercellular and intracellular evidence, including a myocardial-infarction application. Its weak correspondence between expression-only ligand–receptor inference and spatial colocalization motivates independent spatial controls. T22. |
| R19 | [CellOracle](https://www.nature.com/articles/s41586-022-05688-9), Nature, 2023 | Network-based perturbation hypotheses can be useful when compared with compatible measured public perturbations. In-silico regulatory shifts cannot establish treatment efficacy or patient trajectories. T16, T22. |
| R20 | [Interpreting cis-regulatory interactions from large-scale deep neural networks (CREME)](https://www.nature.com/articles/s41588-024-01923-3), Nature Genetics, 2024 | Context removal, regional dependence and interactions suggest tests of distal buffering beyond local saturation. Keep model perturbations distinct from observed regulation. T03, T19, T22. |
| R21 | [Agarwal et al., Massively parallel characterization of transcriptional regulatory elements](https://www.nature.com/articles/s41586-024-08430-9), Nature, 2025 | Large lentiMPRA data and MPRALegNet provide task-specific information and comparators. Random integration does not preserve a variant's endogenous locus; retain construct and assay distinctions. T10–T11, T23. |
| R22 | [Zeng and Li, Predicting RNA splicing from DNA sequence using Pangolin](https://link.springer.com/article/10.1186/s13059-022-02664-4), Genome Biology, 2022 | Liver splice-site usage is a relevant specialist readout. Site usage is not the same object as every intron/exon event, so compare only compatible definitions. T17. |
| R23 | [MTSplice predicts effects of genetic variants on tissue-specific splicing](https://link.springer.com/article/10.1186/s13059-021-02273-7), Genome Biology, 2021 | A shared constitutive component plus tissue adjustment is an alternative to a fully conditioned model, especially for compatible cassette-exon labels. T13, T17. |
| R24 | [Deciphering the impact of genetic variation on human polyadenylation using APARENT2](https://pmc.ncbi.nlm.nih.gov/articles/PMC9636789/), Genome Biology, 2022 | RNA-end usage and cleavage can explain effects missed by accessibility or total RNA. Its tissue-resolved residual approach motivates an APA-specific specialist comparison rather than relabeling every RNA effect as splicing. T17. |
| R25 | [Ahlmann-Eltze et al., Deep-learning-based gene perturbation effect prediction does not yet outperform simple linear baselines](https://www.nature.com/articles/s41592-025-02772-6), Nature Methods, 2025 | In the tested settings, complex models did not beat simple perturbation baselines. Require no-change, additive and task-native linear controls; the result is not a universal impossibility claim. T14–T15, T19–T20. |
| R26 | [Benchmarking algorithms for generalizable single-cell perturbation response prediction](https://www.nature.com/articles/s41592-025-02980-0), Nature Methods, 2025 | Distinguish unseen perturbations from unseen contexts and examine several complementary endpoints. This motivates a factorized evaluation matrix, not another unrestricted benchmark sweep. T20. |
| R29 | [Elison et al., Single cell multiomics reveals drivers of metabolic dysfunction-associated steatohepatitis](https://pubmed.ncbi.nlm.nih.gov/40385416/), medRxiv preprint, 2025; [methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC12083587/) | Existing MASLD use of SCENIC+ and multiomic regulatory interpretation. Treat its findings as source-owned and identify cohort overlap before using any data. This directly narrows novelty claims. T01, T16, T22. |

#### Public evidence and source admission

| Ref | Primary reference / status | What to take from it and where it matters |
|---|---|---|
| R27 | [Kerimov et al., A compendium of uniformly processed human gene expression and splicing quantitative trait loci](https://www.nature.com/articles/s41588-021-00924-w), Nature Genetics, 2021; [current official data access](https://www.ebi.ac.uk/eqtl/Data_access/) | Source-native molecular objects, uniform metadata and regional access support targeted public expansion. A database reprocessing of an existing cohort is not a new independent cohort. T01, T17, T23–T24. |
| R28 | [MaveDB 2024: a curated community database with over seven million variant effects from multiplexed functional assays](https://pmc.ncbi.nlm.nih.gov/articles/PMC11753097/), Genome Biology, 2025 | Stable assay/score-set identifiers and explicit measurement types help discover suitable variant labels and expose their provenance. Admission is assay-specific; protein-abundance measurements do not become regulatory labels. T17–T19, T23–T24. |

Read R03/R07/R10/R12/R13 first for immediate design decisions; R14–R20 for mechanism and remodeling; R21–R24/R27–R28 for public label expansion; R02/R08–R09/R25–R26 for controls and limits. Existing section 2.1 candidates remain eligible, with the narrower novelty interpretation above.


### 22.3 Detailed comparisons to add to T01–T24

These designs use the existing task definitions for source admission, biological units, multiplicity and protected evaluation. All numerical settings below are **proposed development values**, not reported paper defaults or predictions of success. Select a small informative subset first. Finish with the original task's complete result table even when the outcome is negative.

#### 22.3.1 Align supervision with the molecular quantity (T01, T08)

**Goal:** determine whether broader training fails because marginal QTL associations, allelic imbalance and experimental allele effects measure different quantities. **Motivation:** E03–E04, E25 and R11/R27. This is a first development comparison, not a new fine-mapping release.

**Inputs and variables:** compare selected-lead effect/SE labels, broader eligible associations, and measured allelic counts in separate heads. Where complete summary statistics and reliable matched LD exist, add a model that predicts the marginal association implied by its underlying variant contributions. Derive the required genotype/phenotype scaling and sampling covariance before implementing that branch; a marginal effect is not automatically an isolated sequence effect. Compare uniform and capped precision weights, robust versus Gaussian residuals, and direct versus uncertainty-aware labels. Joint total-count/allelic likelihoods require actual donor pairing and a dependence model.

**Decisive test:** hold loci and, where applicable, donors fixed across arms. Evaluate against independent assay-native measurements, not against labels shrunk using the same predictions. Count a broader eligible population and the matched overlap separately. The output must show whether label definition, population expansion or model capacity explains any gain.

**Stop or escalate:** without compatible LD/scales, omit the LD-aware branch and retain association prediction as its own endpoint. New COLOC or fine-mapping memberships are outside this task. Success is a measured improvement with the estimand made clearer, not a larger apparent training set.

#### 22.3.2 Learn shared and context-specific effects with uncertainty (T08, T14–T15)

**Goal:** improve precision for sparse contexts while preserving real context differences. **Motivation:** E20, E26 and R10/R30.

**Inputs and variables:** within one compatible molecular endpoint, compare independent context estimates, a shared-effect model and multivariate adaptive shrinkage. Test diagonal, shared and training-learned covariance families; a low-rank component of rank 1, 2 or 4 is considered only when the number of measured contexts supports it. Estimate sampling-error covariance when contexts share donors or controls. Cross no shrinkage versus shrinkage with the direct interaction objective. Keep RNA, accessibility and reporter effect scales separate.

**Decisive test:** training-only sharing patterns must improve held-locus or held-donor likelihood/error and preserve calibrated context contrasts. Evaluate small and large differences; a model that predicts the shared mean well but suppresses all interactions does not meet the goal. A fitted sharing coefficient does not establish biological universality.

**Output and completion:** context-specific effect/uncertainty tables and an explicit estimate of when pooling helps. Two treatments provide limited information about a general covariance structure; in that case use simple partial pooling. Stop increasing covariance complexity when held-out performance fails to improve or uncertainty becomes unstable.

#### 22.3.3 Preserve allele effects during adaptation (T05–T09)

**Goal:** determine whether adaptation improves the small REF–ALT contrast while retaining useful native information. **Motivation:** E33–E34, R01/R02 and the distinction between profiles and differences in R07–R08.

**Inputs and variables:** compare a shared per-allele head with an explicitly antisymmetric contrast head, using the best supported local pooling. Cross fixed features versus encoder/adapter updates with effect-only and profile-plus-effect objectives. Separate two optional retention terms: native mean-profile consistency and native allele-contrast consistency. Test normalized weights 0, 0.1 and 1 in development; a strong contrast-retention term could prevent the intended correction and must not be assumed beneficial. Begin with two adapter depths or trainable-component choices, not the full Cartesian menu.

**Decisive test:** identical training exposure and target geometry; local zero-shot and calibrated-score baselines; independent allele error/correlation and separate native-profile retention. Verify allele-swap antisymmetry and zero change for identical alleles. Inspect weak effects with higher-precision reductions before interpreting small improvements.

**Output and completion:** a performance-versus-retention comparison identifying a transferable improvement or a clearly labeled task specialist. Increase depth only if the gain survives paired development comparisons and learning curves justify the additional trainable parameters.

#### 22.3.4 Train shared liver grammar before lineage-specific refinement (T07, T09, T12)

**Goal:** learn adult lineage specificity without discarding broadly useful sequence information. **Motivation:** E26 and the consensus-to-specific curriculum in R03.

**Inputs and variables:** use eligible donor/lineage ATAC profiles. Compare consensus-only training, lineage-specific-only training and consensus training followed by lineage refinement. Define lineage-specific regions from training donors only. Allocate 0%, 25% or 50% of a fixed update budget to refinement; compare shared heads with partially pooled lineage heads. Keep the original full region universe in evaluation and report the lineage-specific subset separately. Include a competitive compact sequence model and eligible ChromBPNet/Borzoi/AlphaGenome arms.

**Decisive test:** held donors and chromosomes/regions, with equal training exposure and depth-aware evaluation. Report broad profile fidelity, lineage discrimination and measured allele effects separately. Extra cells in a rare lineage improve its measurement but do not replace independent donors.

**Output and completion:** a lineage-by-assay comparison and a recipe whose improvement is not confined to the selected specific peaks. Advance to allele adaptation only where genotype-linked labels support it; do not claim that this curriculum fixes the existing 3-site cholangiocyte allele-validation limit.

#### 22.3.5 Train donor contrasts explicitly (T08, T15, T18, T20)

**Goal:** predict which donors differ at the same molecular target, beyond a target's average profile. **Motivation:** E36 and R07–R09.

**Inputs and variables:** compare a train-donor-mean residual objective, pairwise donor-difference loss and their combination with the profile loss. Fit target scaling in training data. Sample 2 or 4 partners per training donor per epoch versus an affordable all-pairs development arm; donor-balanced weighting prevents prolific pairs from dominating. Cross random same-source pairs with technical-stratum-matched pairs, while avoiding matching away the biological contrast being studied. Compare sequence-only, RNA-only, additive and interaction models.

**Decisive test:** a four-cell evaluation matrix: seen/new donors × seen/new loci. Neither donor in a held-out pair can enter training for the strict donor test. Pairs are dependent; uncertainty is based on donors and locus blocks, not the quadratic pair count. A new-donor gain at trained loci remains valuable but does not support new-locus generalization.

**Output and completion:** paired-difference error, calibration and the four-cell matrix, with mean profiles reported separately. Actual genotype-linked inputs are required for personal-genome claims; RNA-conditioned donor prediction is a different capability. More complex conditioning is justified only by a gain over the strong RNA-only comparator.

#### 22.3.6 Compare what RNA and chromatin context contribute (T13–T16)

**Goal:** distinguish sequence information from information supplied by an observed molecular context. **Motivation:** E36–E37 and R04–R06/R15.

**Inputs and variables:** compare regulator RNA, broader RNA components, fixed117 scores and eligible paired ATAC inputs; then a small combined representation. Begin with 8, 16 and 32 context dimensions where donor rank permits. Compare linear/additive, bilinear and FiLM conditioning before attention. Separate lineage means from donor residual context and measured from imputed modalities. GET-like expression prediction from observed cis chromatin is explicitly an imputation task; it does not establish prediction of accessibility or an allele effect.

**Decisive test:** target/cis-feature exclusion applies to the context encoders used for strict donor/allele prediction. Where observed cis chromatin is deliberately supplied for RNA imputation, report that distinct input contract and never reuse the target assay as its own predictor. Compare shuffled context within source, technical-only context and missing-modality fallbacks. For an allele test, keep the measured context fixed across REF and ALT.

**Output and completion:** an information-source ablation identifying a gain beyond RNA-only and sequence-only models on the stated held donor/region endpoint. Stop if the gain depends on target reconstruction or source identity. New paired-modality integration of the whole atlas is not required.

#### 22.3.7 Separate variant effects from enhancer–gene links (T16, T20, T24)

**Goal:** improve candidate target discrimination by combining two independently assessable relationships: variant→region effect and region→gene link. **Motivation:** E07–E08, E16 and R01/R12/R13.

**Inputs and variables:** compare proximity, current genetic/annotation evidence, ENCODE-rE2G, scE2G-ATAC and, where truly paired data exist, scE2G-Multiome. Then add AlphaGenome/local accessibility and expression contrasts as separate features in a regularized link/effect model. Compare population-aggregated and lineage-specific inputs on matched coverage; distinguish a missing assay from a measured absence of activity. Begin with linear/logistic or additive heads before a graph neural network.

**Decisive test:** hold out loci and source families; audit overlap with the public CRISPR training/evaluation sets used by the source models. Evaluate element perturbation, variant effects and genetic target support separately. Negatives require adequate assay power; untested links are not negative examples. Compare AUPRC and precision/coverage on the same candidate universe, with expression direction as a separate endpoint.

**Output and completion:** target alternatives with evidence for each edge and a reproducible gain over the best eligible linker. Do not interpret a product of posterior, link score and allele effect as a mechanism probability. If only non-liver validation is available, label the liver assignments as predictions pending target validation.

#### 22.3.8 Add allele information to regulatory networks (T16, T22)

**Goal:** test whether sequence-sensitive edges improve an existing MASLD regulatory network's explanatory or predictive value. **Motivation:** E12, E39 and R14/R29, including prior SCENIC+ use in MASLD.

**Inputs and variables:** fit or reuse eligible TF–region–gene networks within training data; compare the ordinary network, a network augmented with sequence-effect features and a version with matched shuffled features. Consider motif-family compatibility, region accessibility and allele-effect direction as distinct annotations. Retain ambiguous TF families and alternative targets. Compare no filtering with training-defined stability filtering; quantify edge stability by resampling donors rather than only cells.

**Decisive test:** improvements must concern held-out target links, measured regulator-response relationships or independently supported locus interpretation. Match controls on node degree, expression, region activity, motif frequency and target distance. Do not validate against a gene set used to build the network or optimize fixed117 to improve agreement.

**Output and completion:** a small set of supported network differences plus the full tested set and failures. A new network image or generic TF enrichment does not meet the goal. Without independent link/response evidence, retain an interpretation aid and avoid claiming a new regulatory mechanism.

#### 22.3.9 Model reporter assay context explicitly (T03, T10–T11)

**Goal:** determine which transfer failures arise from the reporter construct, genomic placement or receiving assay. **Motivation:** E27–E28, E34 and R21.

**Inputs and variables:** compare assayed insert only, documented full construct and native flanking sequence as separately defined inputs. Where the dataset actually varies these factors, add promoter identity, reporter delivery method as recorded in the source metadata, and cell condition. Compare assay-native heads with a shared effect representation and a residual transfer model. Hold the amount of endogenous training fixed while adding reporter data; compare matched target/context additions with equally sized unmatched additions.

**Decisive test:** the receiving endogenous assay must improve on held loci. Keep lentiMPRA distinguished from the endogenous native locus. If a promoter or vector is constant in a dataset, it cannot support a learned promoter/vector-effect claim. Report changes in assay coverage and baseline activity rather than attributing all transfer differences to biological context.

**Output and completion:** a transfer map identifying where shared information helps and where construct/native-context mismatch limits it. No label-generating synthetic reporter experiments are required. Expand the model only after the metadata support the factor being estimated.

#### 22.3.10 Resolve splice-event competition before increasing model size (T17)

**Goal:** improve a correctly defined event's detection, magnitude or direction beyond E10's completed signed-splice result. **Motivation:** E09–E10 and R22–R23.

**Inputs and variables:** compare the retained loose model-junction denominator, a source-native cluster and a fixed annotation/reference-defined competitor set. Fit any reference-usage threshold in development with the score's actual units documented. Compare local site scores, paired donor/acceptor features and targeted base-resolution features with long-context summaries. Include Pangolin for compatible liver site-usage labels and MTSplice for compatible cassette-exon labels; use a common evaluation subset and report excluded events.

**Decisive test:** prespecify rules for novel/lost junctions, REF/ALT coordinate correspondence and the union of eligible event identities. Do not filter on the observed ALT outcome or on correct prediction. Evaluate event identification, effect magnitude and direction separately, retaining proximity and E10 as appropriate baselines. Report all-pair and prespecified lead-only results separately.

**Output and completion:** an event-definition comparison showing whether alignment improves the intended endpoint independently. A better secondary subset cannot replace E10's primary negative threshold decision. Only then test a larger event head or decoder adaptation.

#### 22.3.11 Add a distinct RNA-end mechanism where labels exist (T13, T17)

**Goal:** determine whether alternative cleavage/polyadenylation explains regulatory effects missed by accessibility, total expression or splice-junction outputs. **Motivation:** weak RNA direction in E02 and the event-specific approach in R24.

**Inputs and variables:** admit public cleavage/APA reporter or molecular-QTL labels with exact strand, transcript-end coordinates and competing sites. Compare native AlphaGenome-compatible readouts, APARENT2 and calibrated specialist residuals. Test local polyadenylation-signal context versus additional transcript context. Where enough paired context exists, compare a small RNA-binding-protein panel with the TF-only context used for accessibility and with a shared-plus-tissue residual. Do not assume a pretrained liver-specific APA output exists.

**Decisive test:** held genes/loci and compatible source families; distinguish cleavage location, relative site usage and total transcript abundance. GTEx-derived evaluation remains exposure-qualified. A source endpoint is usable only if its native quantity can be matched, including transformed QTL units.

**Output and completion:** a bounded APA capability or a clear absence of usable MASLD-relevant labels, plus event records in the database. Do not infer RNA stability or protein effects merely from predicted APA changes. This is conditional on admitted labels, not an obligation to build another model.

#### 22.3.12 Use real diploid backgrounds with the correct readout (T18–T19)

**Goal:** establish whether actual cis backgrounds explain molecular variation beyond focal-variant and dosage baselines. **Motivation:** E29–E30 and R07–R09/R11.

**Inputs and variables:** where public genotype-linked molecular data permit, compare focal allele on reference, verified local variants and complete supported local haplotypes. Score both haplotypes. Define diploid aggregation in the assay's physical scale before applying logs or ratios; averaging log outputs is not generally equivalent to aggregating expected abundance. Compare observed phase with uncertainty over admissible phase configurations, and background-independent with background-aware residual heads.

**Decisive test:** new donors and loci, matched readout geometry and genotype-based linear/additive controls. Report background effects separately from a joint-interaction contrast. Mapping-bias filtering that supported SNV allelic analyses does not automatically validate indel labels; check the actual alignment method and exclusions.

**Output and completion:** better independently measured molecular predictions or a precise limit, with phase/coverage uncertainty retained. A genotype-derived prediction and an RNA-context-derived prediction need separate ablations. Stop at predicted alternatives if public data cannot identify the required contrast.

#### 22.3.13 Test multicellular coupling with fixed programs (T22)

**Goal:** identify whether genetically supported regulatory processes connect to coordinated remodeling across liver lineages. **Motivation:** E05, E37–E40 and R16–R17/R30.

**Inputs and variables:** use fixed117 donor-by-lineage scores and admitted genetic/target evidence. Borrow DIALOGUE's multicellular question through a constrained covariance/prediction analysis; do not rediscover programs. Compare within-lineage scores, cross-lineage predictors and a low-rank shared component of rank 1, 2 or 4 only if donor support permits. Adjust for recorded stage, source and composition; compare residualized and covariate-inclusive models without using the held target's information.

**Decisive test:** held-donor and, where possible, held-cohort prediction of a prespecified receiver-lineage program, beyond stage/composition and within-lineage baselines. Match co-expressed gene sets and preserve fixed117 coverage. Keep direct-MASLD genetic support separate from enzyme-trait support and common fibrosis-associated patterns.

**Output and completion:** a reproducible cross-lineage relationship with an incremental effect size and unresolved direction/mediation. This is an observational relationship. Repeating the completed module screen with a different cutoff, or adding a program score to its own predictor, does not meet the goal.

#### 22.3.14 Require spatial support for signaling hypotheses (T22)

**Goal:** distinguish spatially supported intercellular hypotheses from correlations caused by changing cell abundance. **Motivation:** existing protein/spatial coverage in E40 and R18–R19.

**Inputs and variables:** start with a small prespecified set of inherited-evidence-linked targets and measured sender/receiver lineages. Compare expression-only ligand–receptor evidence, spatially supported evidence and compatibility with downstream TF/program responses. If coordinates are calibrated, compare 50, 100 and 200 micrometre neighborhoods as development choices appropriate to the assay resolution. Analyse proteins only when the measured analyte and cellular attribution support that connection.

**Decisive test:** preserve donor, section and tissue-compartment structure in spatial controls; simple unrestricted spot permutations destroy the relevant null. Compare abundance-matched regions and donors. Target/ligand/receptor coverage must be explicit for targeted spatial panels; missing panel genes are not absent expression.

**Output and completion:** replicated donor-level spatial support for a specified regulatory hypothesis, beyond composition and expression-only inference. Neither proximity nor a sign-consistent network establishes signal transmission, inherited mediation or treatment benefit. If the required targets are unmeasured, return the observability limit.

#### 22.3.15 Test distal buffering with geometry-matched model perturbations (T03, T19, T22)

**Goal:** ask whether local allele sensitivity changes with surrounding regulatory context in a way supported by independent measurements. **Motivation:** E17/E23/E30 and R20.

**Inputs and variables:** reconcile full P4 first and reuse its features. On a bounded locus set, compare native context, removal of a nominated distal region and matched sequence-control perturbations. Use two prespecified regional scales appropriate to the actual target geometry; compare single versus paired region perturbations. Include matched background activity, GC/dinucleotide composition, distance and window boundaries. Distinguish local saturation from full-window context dependence.

**Decisive test:** all arms use the same molecular object and readout scale; assess stability across control sequences, shifts and selected model families. Compare predictions with public enhancer perturbations or measured natural-background effects where compatible. An arbitrary shuffle is a model diagnostic, not a naturally occurring allele or a measured enhancer deletion.

**Output and completion:** supported context dependence or a prediction-only hypothesis with the discriminating observation specified. Do not call model nonlinearity epistasis or biological buffering without matched measurements. No new full saturation scan is justified until existing P4 completion and reuse are resolved.

#### 22.3.16 Calibrate the public capability, including sparse contexts (T20–T24)

**Goal:** make the public resource accurate about both its predictions and their limits. **Motivation:** E21/E31/E38 and R02/R25–R28.

**Inputs and variables:** compare score-only calibration, measurement-variance-aware residual calibration and grouped conformal calibration where exchangeability is defensible. Evaluate nominal 80%, 90% and 95% coverage with interval width. Compare within-source calibration, transported calibration and limited new-source recalibration on disjoint donors/loci; size the calibration set by available independent groups, not a fixed percentage alone. Evaluate missing regulator genes and missing modalities using recorded panel masks rather than arbitrary gene dropout alone.

**Decisive test:** error-versus-retained-coverage curves, calibration by native endpoint and adequately supported context, and full-versus-compact public model comparisons. A restricted research model's performance cannot be assigned to a public student without testing the student. Distinguish an interval for a future noisy measurement from uncertainty about the latent molecular effect.

**Output and completion:** a reproducible query response carrying applicable context, source versions, native effect, supported uncertainty and an honest unsupported answer. Use R27/R28-style molecular/score-set identifiers to connect external records without counting a duplicate source as replication. No distribution-free guarantee is claimed under arbitrary cohort shift.


### 22.4 Targeted public-data searches with a concrete use

This table adds source routes to section 5.1, not admitted training sets. Inspect existing local copies and cohort registries first. A “public” record does not prove usable labels, independent donors or distributable model derivatives.

| Source route | Information worth seeking | Admission decision before acquisition |
|---|---|---|
| R12/R13 enhancer-link benchmarks and released maps | Power-qualified element–gene perturbations; exact training/test partitions; adult-liver or relevant lineage maps; available DNase/ATAC/Hi-C inputs | Reconcile overlap with AlphaGenome/GET/linker benchmarks. Element inhibition is a link label, not a single-nucleotide effect. A predicted liver map is not liver validation. |
| R21 lentiMPRA and existing reporter deposits | Actual sequence/construct, promoter and cell context; element activity versus variant effects; replicate and barcode structure | Admit only measurements matching T10/T11's target. Resolve publication-family overlap with existing reporters; never treat random integration as the native endogenous context. |
| R22–R24 and R28 functional RNA assays | Exact splice or cleavage events, quantitative variant effects, combinations where measured, strand and coordinate mapping | Count independent loci/events rather than variant rows alone. Distinguish assay-specific activity, usage and clinical classifications; check which source examples trained each specialist. |
| R27 eQTL Catalogue | Full eligible regional association fields, effect/SE, event identifiers and independent context/source metadata | Use a recorded release and dataset identifiers. A new processing pipeline for GTEx or another existing source adds no new participants. Keep all protected-source restrictions. |
| R29 MASLD multiomics study | Potential source maps, public processed files, donor joins and assay availability for method comparisons | Prior MASLD method use is already established. Verify actual public availability and reconcile overlap before assigning a data role; the preprint's existence does not establish release of every underlying assay. |
| [GSE292268 public spatial record](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE292268) | CosMx MASLD spatial measurements and fibrosis-morphology context for T22, if the target panel, donor identifiers and coordinates permit the stated hypothesis | Metadata candidate only: donor count, independent cohort status, panel coverage, quantitative outputs and reuse terms must be verified. Spatial recurrence cannot validate a sequence allele effect. |
| Existing public liver and permitted phased reference sources | Genotype-linked molecular measurements, actual haplotypes and the overlap needed for T18/T19 | Public reference haplotypes supply sequence backgrounds only; they do not supply donor-linked RNA/accessibility labels or measured epistasis. Prioritize a small compatible overlap over unpaired sample volume. |

Rank admitted sources by the independent information they add: new loci, new donors, measured contexts, target/event identity or an assay bridge. Record the most important expected gain and the bottleneck it addresses, rather than ranking by total file size or number of cells. No raw-data download is needed to complete this source-admission proposal.

### 22.5 A practical order for the expanded search

1. **Reconcile labels and baselines.** Resolve the T01/T08 estimand comparison and the latest local baseline. Inspect the public enhancer-link benchmark partitions and the existing full-P4 state. These checks determine which later comparisons are interpretable.
2. **Use the existing data for three focused model questions.** Compare allele-focused adaptation/retention, donor-contrast objectives and consensus-to-lineage refinement. Keep their original task-native comparators and existing development budgets.
3. **Test the most consequential added information.** Admit enhancer–gene models and measured-chromatin context where the inputs exist. Compare the simple original method with the proposed MASLD adaptation before adding graph or attention complexity.
4. **Expand molecular endpoints selectively.** Resolve splice-event geometry before scaling a model; add APA or construct-aware transfer only when labels and context metadata meet the admission rules.
5. **Interpret biological relationships after prerequisites.** Reuse completed P4 and fixed117 for the defined multicellular, spatial and buffering hypotheses. Carry all negative and underpowered results forward.
6. **Evaluate and release the supported capability.** Fix the recipe and endpoint, use independent evaluation where eligible, calibrate the public model itself, then attach complete evidence and limits to the database.

For every comparison, the coordinator records the decision before execution: what result would justify continuing, what precise negative would stop this branch, and what missing information would leave it unresolved. A gain that appears only after source-specific tuning, a changed eligible population or repeated outcome inspection cannot support the original generalization claim. Important interactions from section 18 remain in scope; this expansion does not multiply every variable into one enormous grid.

### 22.6 What these adaptations could add to the paper

| Potential contribution | Minimum scientific evidence | Natural home in the existing paper |
|---|---|---|
| Better liver allele-effect prediction from precise supervision or adaptation | Matched gain over the strongest eligible baseline; independent assessment, calibration and coverage | Figure 7 and the public model |
| A donor-context model that predicts differences rather than average profiles | New-donor and new-locus results beyond RNA-only and sequence-only controls | Figure 7, linked to cellular-context results |
| A liver variant-to-target resource with explicit alternatives | Better target discrimination or clearly qualified predictions, with separate allele-effect and region–gene evidence | Existing genetic interpretation and Gene Catalog |
| A useful RNA-processing extension | Exact event or APA prediction with assay-native validation, not an unsigned score relabeled as direction | Molecular interpretation panels and event queries |
| A supported connection between inherited regulation and multicellular remodeling | Reproducible donor-level relationship beyond state, composition, expression and region-property controls | Existing remodeling/program/spatial figures |
| A novel, useful adaptation from another disease | Clear prior-method/MASLD comparison and a new tested biological insight or capability, even without new architecture | The relevant biological result plus a concise methods comparison |

The ambitious outcome is a coherent answer to the central question using the strongest method for each supported relationship. It may be a shared model with specialists and a common evidence database. Architectural uniformity is not itself a scientific objective. Completed comparisons determine which of these contributions belongs in the manuscript.

## 23. Two useful models and a compact Figure 7 candidate

**Scope fixed 2026-09-22; candidate work, not adopted.** The near-term pair is
(1) liver RNA to measured H3K27ac profiles, for differences among participants
with comparable histology and for choosing the next chromatin measurements, and
(2) signed liver caQTL effect prediction, for variant–regulatory-element
experiments. The released severity/activity model remains a reference. The
fixed task, comparisons and margins are in
`../../scripts/analysis/masld_two_models/design.json`; the figure region and
example choices are in `acquisition_figure_addendum.json` beside it. These
choices do not change the fixed 117 programs or Resource evidence classes.

For chromatin, use the 99 GSE267145 biological participants and existing five
participant folds. The target is H3K27ac log2 CPM after fold-training removal
of four H3K27ac signal-concentration descriptors. Compare released global RRR
and global-plus-local predictions with an inner-tuned ridge using four
histology grades, sex and RNA technical descriptors; add estimated hepatocyte
mRNA fraction only in a separate, cohort-adapted sensitivity. This tests improvement over the
specified linear covariate model; it does not establish general conditional
independence from histology or a regulatory mechanism. A separate exact-grade,
sex- and fold-matched analysis recovered chromatin differences in 22 disjoint
pairs, with global-plus-local skill 0.11070 versus global skill 0.08859.
It represents fibrosis 0–1 only, mainly women with fibrosis 0; no fibrosis
2–4 pair was available. The released batch-marginal
profile requires at least ten scored samples; its fixed raw route scores one
individual with weaker observed ATAC transport. The released wide-RRR rank
grid already reached 56, and its five selected ranks were 12, 11, 11, 14 and
14. None met the rank ceiling, so the initially planned rank-80 repeat on the
same cohort is unnecessary; see `rank_resolution_addendum.json` beside the
candidate design. Matched held-participant input perturbations now show that
gene-count coverage alone misses large compositional differences: omitting
the most abundant 10% removes about 90% of RNA mass, defined here as summed
input abundance among the 42,163 modeled genes rather than physical RNA
mass or a whole-transcriptome fraction. Raw predictions also
deteriorate under random 60% coverage. Batch-marginal mapping is more robust
to these perturbations, but depends on batch size/composition and offers no
clear full-input advantage. Library scaling is invariant; observed read-depth
accuracy is not tested from fractional estimates. GSE296875 ATAC is an inspected
different-assay transport source, not external H3K27ac confirmation. A
numerical check of the released scorer on source-exposed RNA found exact raw
library-scale invariance and an RMS profile shift of 0.0481 when the same ten
people were batch-marginally scored with ten extra people. At 90% gene
coverage, dropping high-abundance genes moved profiles far more than dropping
low-abundance genes. This is numerical sensitivity, not measured prediction
accuracy under missing inputs. A separate held-participant test then confirmed
missing-input effects. An adaptive training-reference RNA-mass correction
improved the fixed individual route under the tested missing-gene patterns
(90/75/60% random retention skills 0.14651/0.13600/0.12177) while leaving
complete-input predictions unchanged. It is an optional candidate sensitivity,
not a universally superior route; the assumption is that retained genes carry
RNA mass comparable to the training reference. Its evaluation uses the bounded
rank-12, 1,000-region recipe, not full released weights. Details are in
[RESULTS](../RESULTS.md#two-model-molecular-value-development-2026-09-22-not-adopted).

One additional sequence-by-RNA residual experiment was run at 1,024 fixed
regions. Frozen AlphaGenome sequence features and training-only RNA/sequence
PCA entered a regularized interaction of ranks 4, 8 or 16, selected within
each outer training set. On held participants, its residual skill was
**−0.00056** against zero correction, compared with **−0.00054** for the
additive sequence-plus-RNA control and **−0.02401** for RNA-only. Its increment
over the best simple control is **−0.000018** (paired-participant bootstrap
interval [−0.00170, 0.00157]). This uses residual SSE as denominator, whereas
the +0.01 complexity margin concerns original-target skill. Independent
reconstruction gives −0.0000159 [−0.001520, 0.001375] on the latter scale.
Its negative increment cannot qualify on either scale; retain
the simpler chromatin model. This experiment did not test unseen-region
generalization. Sequence features retain the permitted `counted_interval`
semantics and recorded 1-bp uncertainty, without establishing called-peak
boundaries. Inner choices used cross-fitted base-model residuals generated
once within each outer training fold; the outer participants remain untouched,
but the comparison is not a fully nested test of the complete base-model
selection. Exact errors and choices are in the candidate
`sequence_rna_residual_21849236/`.

The bounded acquisition experiment used 1,000 outcome-blind,
nearby-gene-annotated regions and five held-participant folds. At the primary
20% additional budget (20 of 99 participants, on top of 45 base measurements),
RNA diversity and global/local disagreement found more candidate
pool error than random selection, but neither delivered a Holm-significant
refit MSE advantage on separate evaluation participants. The refits share
training participants. Fold sign flips are dependence sensitivities, not exact
biological-unit inference. A 99-participant paired bootstrap also leaves the
refit advantages uncertain (RNA diversity +0.001486, disagreement +0.001323;
both intervals include zero, Holm P=0.4962 within three policy contrasts).
This uncertainty is conditional on the fitted models. The initial histology policy deterministically
balances exact-grade strata; proportional stratified random sampling is a
separate completed check, with participant-weighted primary refit reduction
0.00721 versus random 0.00629; this is descriptive. Global rank-12 RRR uses all 96,460 responses, with local
correction and evaluation restricted to 1,000 regions. Exact values and conditional-comparator results
are in [RESULTS](../RESULTS.md#two-model-molecular-value-development-2026-09-22-not-adopted).

For allele effects, keep Currin fold 0 closed after its recorded same-input
adapter confirmation. Development folds 1–4 compare native AlphaGenome 1 Mb,
the five-seed 2-kb adapter, and a fixed-penalty combination fitted from
fold-excluding meta predictions. Native 2 kb and the matched fixed head are
controls. All report signed Currin FastQTL ALT-dosage beta, coverage, paired
chromosome-group uncertainty and historical 1-Mb-bin sensitivity. A 5%
relative MSE reduction over the better simple model is the development
complexity margin. External Wenz/GSE277774 admission still requires final
association assembly/alleles, donor overlap, checkpoint exposure and a usable
effect field; RASQUAL Pi cannot be pooled as Currin beta. A public 14,076-row
FDR-selected lead file now supports GRCh38 coordinates for most rows by REF
comparison, but 388 reported REF alleles mismatch the local GRCh38 reference,
and it lacks the all-tested denominator. Its source record is
`source_metadata/wenz_admission.json` under the candidate results. If unadmitted, the
claim stays single-source. The exact-pair variant query rejects arbitrary
peak intervals because the adapter does not consume peak coordinates; target
gene, individual predictive interval and MASLD effect remain unsupported.
Currin beta is an association per additional ALT allele, not an isolated
causal editing effect. Fresh sequence inference remains experimental because
unchanged scoring differed across fresh GPU processes. The reproducible
exact-pair query returns saved predictions; their intervals exclude runtime
variation as well as refitting and recipe-selection uncertainty. Timing remains
an implementation diagnostic without an admitted cost/performance ratio.

The comparison completed on 2026-09-23 at 24,430 matched development leads
(94.88% of 25,748; 17 chromosomes). Native 1-Mb MSE is 0.300161, adapter
MSE 0.302252 and combined MSE 0.290705. The paired gains have conditional
chromosome intervals above zero and Holm P=0.000666 in the two-contrast family,
but the 3.15% relative gain misses the fixed 5% complexity margin. Retain
calibrated native 1 Mb as the selected simple allele model. A descriptive
20% prioritization check captures 26.43% of measured absolute beta with
native 1 Mb versus 20.00% under uniform selection; 98.51% direction agreement
among selected association leads is not an experimental-success estimate.

The completed candidate Figure 7 has six individual vector panels. A and C show measured
and predicted participant chromatin variation, the histology comparator and
the negative sequence-by-RNA result. B and D show allele effects and the
matched native-length comparison, keeping spent confirmation separate. E
shows error capture and held-participant refit value at 10%, 20% and 30%.
F shows one development-selected GNMT-vicinity source pair; proximity is not
a target-gene link. Assay-native axes and distinct region sets stay visible.
The panels, six source tables and caption are in
`figures/candidates/figure7_two_models_20260922_21849207/` from the project root.
All six passed rendered text and visual inspection; the independent numerical
check also rederived all 60 stacking-input splits, four stacking fits, paired
allele intervals and the principal chromatin and figure estimates. Final Figure 7
adoption and any named model release are separate user decisions.
