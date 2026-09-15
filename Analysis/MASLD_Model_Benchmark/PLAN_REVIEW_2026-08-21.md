# Review: MASLD Multimodal Model Benchmark plan

**2026-09-06 update:** See [Multimodal research-model proposal: 2026-09-06](#multimodal-research-model-proposal-2026-09-06)
for the current recommendation using the paper's existing modalities. The September 5
model review and August review below retain earlier reasoning; they are not current status.

**Date:** 2026-08-21 · **Reviewed against:** repo state, GEO/PRIDE records, live SLURM, local reference files.

## Verdict

The machinery is right and the arithmetic is not. Holding back, leakage control, mandatory baselines,
and the preregistered null are better than most benchmark papers ever get. But six of the seven
tracks are disqualified from producing a best model by the plan's own rules, and the one that
survives reuses pass/fail floors this project has already failed twice. Expected outcome after
~120 model adapters is zero best models.

## Ranked issues

**1. One eligible track, with a failing in-project prior.** Variant-to-regulation: "no best model"
without sumstats. RNA-conditioned ATAC: retrospective without FNIH. Graph: "no best model without a
held-back source." Perturbation: exploratory. Unified: "no unified best model is possible." Spatial:
deferred. Bulk: dead on power (#4). That leaves cell-state mapping, whose check (macro-F1 ≥0.70,
every class F1 ≥0.50) is numerically identical to the floors that terminated ROADMAP WP6a: V42
scored 0.2756 on untouched GSE212837; Router V2 hit 0.7071 overall on GSE296875 but failed
per-class at B cells 0.2888 and NK-T 0.4229. The plan never cites this.
→ Cite WP6a as the pilot with those numbers. The binding constraint is rare immune-class F1, not
the overall floor, so either argue why 70 models change that specific behavior or swap the
conjunctive pass/fail check for a preregistered per-class estimation endpoint that is interpretable
either way. Drop "best model" as the deliverable.

**2. GSE289173 cannot be split into two holds back, and its outcomes are not on GEO.** The GEO deposit is
one `GSE289173_RAW.tar` of MTX/TSV matrices, no eQTL/ieQTL summary statistics, no genotypes; the
outcomes live in a single 6.6 GB Zenodo archive (10.5281/zenodo.14586466). Pulling task-1 cell-state
labels puts task-2 sumstats in hand at the same moment, which contradicts "separately held-back" and
"one-time held-back inference."
→ Inventory the Zenodo archive before fixing, confirm the sumstats are cell-type stratified (if not,
"no major lineage worse by >0.03" is uncomputable), and state that GSE289173 supports one
unblinding event.

**3. The hold-back has no enforcement and no external timestamp.** `id` returns uid=91825 only; no second
account or group; `setfacl`/`getfacl` absent; GPFS `mmputacl` is owner-revocable; no sudo.
"Independent evaluator" means an independent code module in the same uid's directory. And the one
primitive that would make it tamper-evident does not land: `git check-ignore -v` confirms
`.gitignore:8:*` swallows `locks/selection_lock.json`, `locks/predictions.sha256`, and
`locks/manifest.tsv`, so no selection record ever reaches origin.
→ Add `!Analysis/MASLD_Model_Benchmark/locks/*` re-includes and make an external stamp (pushed
commit SHA, or OSF/Zenodo write of the record hashes) a hard precondition of wave 9. If no second
party holds the archive, call it a fixed in advance retrospective evaluation, not a hold back.

**4. The bulk check is unreachable by 5-9×.** At n=18 with 5 ordinal stages, simulated SD of the
paired Spearman difference is 0.08-0.14; 80%-power MDE is 0.22-0.44 rank-correlation units, and
power at the registered 0.05 gain is 0.05-0.13, i.e. at alpha. Under the exact null the 0.05
effect-size half of the check is crossed 27-37% of the time. The plan's own "activate only after
prospective 80% power" clause therefore blocks the track by construction.
→ Delete the 0.05. Register the detectable effect (~0.25-0.35) or demote GSE267031 to a
CI-reporting descriptive consistency check.

**5. GSE296875 has donor-linked steatosis and fibrosis histology, but not an adjudicated MASLD
label.** The GSM pages expose only tissue and well, but the paper's Data S1 Table S1 contains BMI,
categorical and numeric macrovesicular steatosis, and categorical and free-text fibrosis for the
de-identified donors. The steatosis field is nonmissing for 38 of 39 analyzed donors and the
fibrosis field for 37 of 39; these counts do not mean all 38 or 37 donors are phenotype-positive.
The public fields do not include NAS, ballooning, lobular inflammation, alcohol exposure, or an
etiology adjudication, so steatosis must not be relabeled as confirmed MASLD or MASH. Donors were
pooled five per 10x inlet, so donor-grouped outer folds alone do not test well-level batch
transport. STATUS.md also records prior project use as a diagnostic-only router cohort, which
makes it project-exposed rather than a new external holdout for this paper.
→ Use the histology as development-only cross-sectional phenotype evidence, retain donor as the
inferential unit and well as a technical batch, freeze the public barcode-to-donor assignment, and
prespecify leave-one-well-out sensitivity. Do not claim disease etiology from these fields.

**6. The reference requirements pins a file nothing can read.** Both SHA-256s verify exactly, and the
truncation call on `gencode_v49/` is correct (3.7 MB genome, 71 KB transcripts, both fail `gzip -t`).
But `refdata-gex-GRCh38-2024-A/genome/GRCh38.p14.genome.fa.gz` is plain gzip (`1f8b 0808`, no BGZF
`BC` field) with no `.fai`/`.gzi`, so no sequence model can window-extract from it; the repo's 17
existing sequence scripts all read an indexed `fasta/genome.fa` that the requirements do not name.
It is also the patch/ALT assembly (706 contigs; 503 of 528 GTF seqnames non-primary; 2,945 symbols
on more than one contig), so a locus duplicated onto a PATCH contig escapes the chromosome/LD-block
holdout and breaks the PredictionBundle row-ID guarantee. `contig_policy` in `reference.py:207` is
inert free text.
→ Pin the primary-only analysis set plus `gencode.v49.primary_assembly.annotation.gtf`, with
fasta + .fai + chrom.sizes + blacklist + per-model bias files each hashed, and implement an actual
contig-membership check. Add `native_annotation_release` to DatasetManifest: the project runs three
live annotations (v49 bulk, v47 in refdata/salmon, v44 for the multiome-ATAC path).

**7. Compute is budgeted; engineering is not.** Campaign approval requires GPU-hours, CPU-hours,
memory, and a hash, and says nothing about person-effort. The census is ~92 checkpoints plus ~28
baselines, each needing an isolated env and a probe/prepare/fit/predict/export adapter. You have 26
micromamba envs, 4 of which touch this census, already spanning six mutually exclusive framework
pins. GPU reality: 18 L40S + 22 B6K institute-wide, 27 pending vs 15 running right now, and your own
B6K sbatch jobs 20208726/20215473 pended from 2026-08-13 and were cancelled at 0:00 elapsed.
[Estimate] adapter/env engineering, not GPU-hours, is the 10× oversubscription.
→ Add a per-model engineering budget line and cap the census by that budget, not by scientific
completeness. State at the top what Resource work this defers, given STATUS.md's standing
"reconciliation, not expansion" and ~30 unpropagated COLOC consumers.

## Cheap fixes before census freeze

- Corgi's released path rank-maps the trans vector to a packaged 2,891-gene reference
  (`predict.py` L61-78), so raw UMI, CPM, and log-CPM give byte-identical inputs. Two of your four
  mapper arms are the same experiment. Keep length-adjusted TPM and the learned encoder. The
  "log-transformed, quantile-normalized, batch-corrected" description of Corgi's training input is
  [Unverified] from open sources; source it or cut it.
- 16 RNA-ATAC integration models are included with no task row, no endpoint, and no check, and the
  latent-factor members cannot emit the only chromatin endpoint on offer. Add a retrospective
  same-nucleus retrieval row on GSE296875 or move the island to "registered but deferred."
- "scPRINT-2 small/medium": no medium-v2 checkpoint exists (small-v1, medium-v1, medium-v1.5,
  large-v1, small-v2 only). Change to small-v2 with SHA at freeze.
- GSE313774 has no resmetirom arm: all 33 GSMs carry only Condition 1/2 ± Ins across 3 batches, no
  donor or sex field. The CXCL1/IL8 exclusion guards a leak that cannot occur, and
  `expected_biological_units = 33` violates your own pseudoreplication rule (n = 3 batches).
- Device eligibility is per-environment, not global: cu13.0 builds do not run on L40S under driver
  550.107.02, cu12.4 cannot target sm_120, TF 2.15.1 has no sm_120 kernels. A cu12.8 torch build
  spans both pools. Standardize on cu12.8 and add `device_eligibility` to ModelManifest.
- The prospective-power rule fixes the variance source but not the effect estimate. State that the
  numerator uses a shrunk or cross-fitted gain (one-SE-below-max, or bootstrap-debiased), not the
  development maximum.
- Being read-only has no retention tier. A 32-bp genome-wide track is ~97M bins (~19 GB per
  model-fold-seed at 100 tracks fp16; ~1.5 TB for a full 7,611-track Borzoi pass), /gpfs is 89%
  full, and `mmlsquota -j sanjana_lab` returns no fileset. Tier retention (metrics, predictions, and
  hashes forever; embeddings and tracks with a declared expiry plus a hashed regeneration recipe)
  and get the real quota before campaign 1.

## Do not let anyone talk you out of these

- Preregistering the null and publishing "no foundation model beat the baseline" as an accepted outcome.
- Mandatory strong baselines retained through final selection, with a baseline allowed to win.
- Explicit pairing topologies and missingness states, with masks instead of missing-as-zero.
- Checkpoint-exposure classification screening best-model eligibility, and the frozen-census no-latest-tags rule.
- Scientific language is already compliant: donors as the unit, "cross-sectional stage-associated
  remodeling," "colocalized"/"genetically anchored candidate," no MR, no "TREAT."

## If you change one thing

Cut the census to the one or two tracks that have a real external endpoint and a reachable check, and
budget adapter and environment engineering explicitly. A 120-model sweep whose own rules permit at
most one best model is the failure mode, not the ambition.

---

## Model-development review: 2026-09-05

**Status:** Recommendations for discussion, not an adopted campaign amendment or model result.

**Goal:** Deliver a prediction or foundation model that measurably changes how MASLD
research is done, which other groups can run on their own liver samples and build on.

**User priorities:** Develop models rather than another tool; use public data only;
explore bulk, single-cell, and joint bulk–single-cell approaches in parallel.

**Review basis:** Current working files, including uncommitted work authorized for review;
September 4–5 result records; model cards; relevant fitting and evaluation code; and
the primary sources linked below. Scores are from stored results, not independently
reproduced in this review. No training, scoring, or compute jobs were run for this review.

### Recommendation

Build a MASLD-adapted RNA model through three parallel tracks. Start with
**BulkFormer for bulk RNA**, **SCimilarity with donor-level learning for single-cell RNA**,
and **a shared model connecting bulk and single-cell representations**.

The central opportunity remains open: the benchmark has not established whether
disease-focused pretraining or adaptation improves MASLD prediction. Much of the
completed work tests classical severity models, frozen cell encoders on broad cell
identity, or specialized regulatory tasks. Those are useful comparisons, but they
do not settle the proposed model-development question.

### What the existing work establishes

| Finding | Evidence | Consequence |
|---|---|---|
| Bulk RNA carries transferable severity information | Released fibrosis model: Spearman **0.6631 in 109 external participants** | Retain it as an incumbent to improve |
| The supervised bulk pool is modest | Latest pool: **557 fibrosis-labeled samples; 464 with NAS** | Prefer adaptation over assuming sufficient disease supervision for a large model from scratch |
| Single-cell disease labels are scarcer than atlas cells | **102 atlas donors**, but **64 donors across five cohorts** in the source-diagnosis disease analysis | Cells do not replace independent disease-labeled donors |
| Existing encoder comparisons principally measure annotation | Five broad cell classes, frozen pretrained representations, fitted heads | They cannot establish MASLD representation quality |
| Bulk foundation-model experiments are unfinished | September 5 synthesis marks Geneformer, TranscriptFormer, UCE, and the neural grid **NOT RUN** | A conclusion that foundation models do not help bulk MASLD prediction is unsupported |
| Cross-scale transfer has a promising starting point | Stored pseudobulk fibrosis result: **rho = 0.5409 in 35 donors** | Joint bulk–single-cell learning merits a bounded experiment |

Evidence: [RNA model card](release/masld-severity-v1.1/MODEL_CARD.md),
[expanded bulk results](executions/transportable-substrate-20260904T214259Z/RESULTS.md),
[donor census](../../docs/RESULTS.md#single-cell-programs-and-composition),
[encoder model card](release/encoder-benchmark-v1/MODEL_CARD.md),
[September 5 representation synthesis](executions/liver-repbench-synthesis/REPORT.md),
and [pseudobulk results](executions/pseudobulk-transfer-20260903T002550Z/results/analysis.json).

The expanded bulk training unit remains **sample**, because donor identity is unresolved
in that pool. Do not describe its 557 samples as 557 verified independent participants.

### Scientific improvements before the next comparison

1. **Make MASLD representation learning the main development question.** The plan's
   106-family census spans many unrelated tasks. Keep completed results as supporting
   benchmarks, but direct new development toward prediction, transfer, and label
   efficiency for MASLD phenotypes. An MPRA result does not establish that a model
   understands a new liver transcriptome.

2. **Separate selection from performance estimation consistently.** The latest bulk
   baseline chooses the best of 984 configurations using the same leave-one-cohort-out
   scores it reports. These are selection scores, not an unbiased estimate of the
   selected procedure. Some competitors tune inside the training cohorts, making
   the comparison unequal. Put feature selection, normalization, adaptation, and
   tuning inside inner study splits; evaluate the selected procedure on an outer
   held study. Bootstrap intervals on selected predictions do not account for the
   entire selection process. See
   [the baseline selection](executions/transportable-substrate-20260904T214259Z/s1_baselines.py)
   and [nested evaluation guidance](https://scikit-learn.org/stable/auto_examples/model_selection/plot_nested_cross_validation_iris.html).

3. **Do not interpret minimum detectable effect as evidence of a null.**
   [The increment evaluator](executions/transportable-substrate-20260904T214259Z/transport_core.py)
   assigns `BOUNDED_NULL_at_MDE` when the absolute observed difference is below its
   MDE. That is not an equivalence test. The latest NAS head adds **+0.1689,
   CI [0.0139, 0.3304]** over clinical covariates, yet fails the overall criterion
   because its increment is below **MDE 0.2260**. Preserve that prespecified stopping
   decision; do not reinterpret it as no NAS signal. Use MDE for planning and
   confidence intervals against a prespecified practical margin for superiority,
   noninferiority, or equivalence. See the
   [NAS decision](executions/two-axis-gate-20260904T230306Z/results/gate.json).

4. **Distinguish unseen during fitting from unused during development.** Repeatedly
   examining an external cohort and changing the next model makes it development
   evidence even when its labels never enter gradient descent. Under the current
   [GSE289173 rules](config/datasets/gse289173.toml), disease/control AUROC is no
   longer reserved after metadata exposure; its remaining reserved tasks differ.
   Final MASLD confirmation needs another eligible public test, or an explicit
   limitation to reused external development cohorts. Audit pretraining exposure
   per evaluated cohort and checkpoint rather than assuming public data are unseen.

5. **Test composition and complexity as competing explanations.** Fibrogenesis can
   change cell abundance and expression distributions. Library-complexity descriptors
   may therefore carry biology as well as technical variation. Compare molecular-only,
   composition-only, technical-descriptor-only, and combined predictors. Claim
   cell-state information beyond composition only when that comparison supports it;
   do not automatically remove useful tissue-level signal.

6. **Correct two interpretations before using them to close development.** A confidence
   interval crossing zero does not establish the RNA model card's phrase “not worse
   than native”; noninferiority needs a margin. The implemented tie ceiling describes
   a perfectly ordered tie-free prediction against tied outcomes, not the universal
   Spearman maximum for predictors that reproduce outcome ties. Keep raw scores
   primary and qualify ceiling-normalized comparisons. See
   [the model card](release/masld-severity-v1.1/MODEL_CARD.md) and
   [the ceiling definition](executions/standardized-metric-panel-20260902T130923Z/panel_metrics.py).

### Track 1: bulk MASLD prediction

**First candidate: BulkFormer-37M.** BulkFormer is pretrained on human bulk RNA-seq,
which matches this assay better than assuming that a single-cell encoder transfers
unchanged. The official release provides multiple sizes and links its pretraining
data, making a source-overlap audit potentially feasible. It explicitly incorporates
expression and missingness summaries, so complexity controls remain necessary.
Sources: [official model](https://github.com/KangBoming/BulkFormer) and
[pretraining resources](https://raw.githubusercontent.com/KangBoming/BulkFormer/main/data/README.md).

Test the following adaptation sequence:

- Frozen embeddings with ridge and a small MLP.
- Limited trainable encoder adaptation with separate fibrosis and NAS heads.
- Continued liver-expression pretraining on eligible development data, followed
  by the same supervised heads.

Use missing-label masks so the 557 fibrosis samples and 464 NAS samples contribute
to their applicable objectives. Predict NAS as NAS, not as an inflammation-specific
output. Do not infer three independent disease aspects from multiple output heads.

Retain the full-expression ridge, PCA with identical heads, and a randomly initialized
architecture control. Add **TabICLv2 on training-fold-fitted 30- or 100-dimensional
representations** as a nonlinear comparator without biological pretraining. It supports
classification and regression; do not feed it 26,629 raw genes by default. Its fitted
prediction context must also be accounted for when assessing reproducibility by another
group. Source: [official TabICLv2 implementation](https://github.com/soda-inria/tabicl).

Continue to larger BulkFormer checkpoints only after held-study improvement or
substantially better label efficiency appears. Pin the exact checkpoint, native
preprocessing, source exposure, and applicable weight terms before fitting.

### Track 2: single-cell MASLD prediction

**First candidate: SCimilarity embeddings plus MixMIL**, compared with liver-trained
scVI features passed through the identical donor-level head. SCimilarity supplies
cell representations; MixMIL was designed to predict sample phenotypes from
collections of cell embeddings. Sources:
[SCimilarity](https://www.nature.com/articles/s41586-024-08411-y) and
[MixMIL](https://proceedings.mlr.press/v238/engelmann24a.html).

Construct each example from one donor's cells, retaining lineage identity and observed
cell counts. Do not assign a donor disease label to individual cells as independent
supervised observations. Compare:

- Donor pseudobulk and lineage pseudobulk predictors.
- Mean cell embeddings and MixMIL.
- Equal-lineage-weight state-only pooling and composition-plus-state prediction.

Do not build disease bags from the existing 50,000-cell annotation screen. It selects
10,000 cells per broad class, so its frequencies cannot represent donor tissue
composition. Build from the original eligible count rows. See the
[sampling implementation](scripts/build_resource_atlas_frozen_screen_50000.py).

Predict source-supported disease states first, with cohort-compatible definitions.
Evaluate fibrosis separately where authoritative grades exist. GSE296875 fibrosis or
steatosis is not a MASLD diagnosis. Report eligible donor counts for each endpoint,
and keep scRNA/snRNA protocol, biopsy/explant provenance, and missing lineages explicit.

After the frozen-feature comparison, test limited encoder adaptation using the donor
objective. The hypothesis is that disease-focused adaptation preserves differences
within cell types that broad annotation overlooks. Demonstrate inference on one held
donor without training on the query cohort, and measure sensitivity to cell depth,
subsampling, and missing lineages.

**Secondary candidate: MrVI.** It is a count-native, sample-aware model for variation
within cell states. Its sample-specific representation is not automatically a frozen
predictor for a new donor. Demonstrate unseen-donor inference before treating it as
the released predictive model. Source:
[MrVI](https://www.nature.com/articles/s41592-025-02808-x).

### Track 3: joint bulk–single-cell learning

**[Hypothesis] This offers the strongest opportunity for a distinct MASLD contribution:**
bulk phenotype supervision could improve donor-level single-cell predictions, while
single-cell structure could improve the bulk representation.

Start with two assay-specific encoders connected through a small shared latent space:

- BulkFormer for real bulk expression.
- SCimilarity with donor pooling for single-cell data.
- Separate phenotype heads and assay-specific reconstruction objectives.

Learn representation consistency between a donor's cell collection and pseudobulk
generated from those same cells. These are two views of one donor, not additional
biological replicates. Real bulk samples contribute phenotype supervision without
invented matched single-cell counterparts. Do not create cross-cohort positive pairs
merely because two samples share a diagnosis.

Preserve measurement differences: suitable single-cell UMI counts can use count
likelihoods; fractional bulk estimates need continuous-expression objectives. Use
explicit masks for missing genes, unavailable labels, and absent lineages. Synthetic
pseudobulk is a training bridge whose transfer to real bulk must be measured.

Compare against **DEGAS**, an existing bulk–single-cell transfer framework, and both
independently trained specialists. Source:
[DEGAS](https://pubmed.ncbi.nlm.nih.gov/35105355/).

Required ablations:

- Shared versus separate representation learning.
- Same-donor consistency versus shuffled donor pairing.
- Pretrained versus randomly initialized encoders.
- Bulk phenotype supervision present versus removed.
- Single-cell information present versus removed.

Add ATAC, protein, methylation, or DNA only when an auxiliary objective improves a
specified RNA prediction task. Failed late fusion does not prove auxiliary learning
must fail, but it argues against an unrestricted fusion campaign. Respect existing
source-specific training and weight-release restrictions; public accessibility alone
does not settle those permissions.

### Evaluation and success criteria

Develop the three tracks in parallel under the same nested study-held-out evaluation.
Use three fixed screening seeds and five finalist seeds. Seeds measure optimization
variability, not biological replication. Model size increases require a measured
benefit in the smaller-model experiments, not a larger-cell-count argument alone.

| Capability | Required evidence |
|---|---|
| Better MASLD prediction | Paired improvement over a strong assay-native baseline on held studies |
| Better label efficiency | Learning curves using nested subsets of training donors or samples |
| Cross-assay transfer | Improvement over direct bulk-to-pseudobulk scoring and independently trained specialists |
| Reusable representation | Frozen features or limited adaptation help more than one downstream task |

Report raw Spearman for graded endpoints, donor-level discrimination for categorical
endpoints, and calibration when probabilities are emitted. Keep cohort-specific
results beside aggregate scores. Use donor-clustered uncertainty where donor identity
exists, and retain the bulk pool's unresolved sample-level replication limitation.
Keep all libraries, repeat biopsies, and assay views from an identified person in
the same partition.

For biological specificity, include molecular-only, metadata-only, and combined
comparisons. State age, recorded sex, BMI, diabetes, and any other included covariates
explicitly, with training-only missingness handling. Never input the endpoint or a
quantity derived from it. Avoid indiscriminate study-invariance objectives when study
and disease are confounded.

Apply Holm correction to the final prespecified superiority comparisons and BH
correction to exploratory biological association families. Define the endpoints,
comparison family, practical margins, and treatment of insufficient precision before
evaluating finalists. Preserve historic stopping decisions rather than retroactively
changing their criteria.

Keep the 117 programs fixed and post hoc. Attention weights, feature importance, and
simulated gene deletion may nominate hypotheses; none establishes perturbation effects.
The Resource and the future Cas13 screen remain separate.

Public data can establish predictive value, transfer, and reproducibility. Actual
adoption requires use by other groups. The immediate deliverable should be a reusable
checkpoint with a convincing MASLD learning result. Preserve the completed ordinal,
anchor-regression, pairwise, and fusion results; do not repeat them without a materially
different hypothesis. These suggestions do not authorize new dependencies, job
submissions, or changes to adopted scientific results.

## Multimodal research-model proposal: 2026-09-06

This is a research recommendation based on the current working files and checked
primary literature. It is not an implemented model, a replacement for the active
campaign specification, or an adopted manuscript claim. No training or new scoring
was run for this review. The two existing RNA instruments remain useful comparators.

### The recommendation

**Build a transferable model of liver molecular state that learns RNA–chromatin
relationships and helps researchers choose samples, mechanisms, and measurements
for follow-up. Make RNA alone the practical entry point.**

[Hypothesis] The strongest opportunity is to predict molecular differences between
livers assigned the same histologic stage, identify which differences have measured
regulatory support, and connect them to the paper's inherited evidence and protein
readouts. A successful model would answer: *Which of my samples represents the
biology I want to study, what evidence supports that interpretation, and what should
I measure next?*

The scientific advance must be demonstrated in withheld molecular measurements.
A larger network, a new embedding, or another association with fibrosis would not
establish that capability. The intended output is a continuous molecular profile,
with uncertainty and assay-specific support, rather than newly named patient subtypes.

The model should learn from multiple assays during development while accepting the
assays a laboratory actually has at inference. With RNA alone it can estimate the
recoverable part of chromatin state. With measured RNA and chromatin it can also
report departures from that expectation. It must not present chromatin information
that is unpredictable from RNA as if it had been measured or recovered.

### What the current data can actually connect

The critical design is the pattern of real assay links. Sharing genes or a fibrosis
label does not make two people a paired training example.

| Existing substrate | Verified scope | Proposed role and limitation |
|---|---|---|
| Bulk RNA | The Resource's validated update has 844 participants across five cohorts. The benchmark's separate training pool lacks a donor key. | Learn transportable tissue states and retain the existing severity output. Keep the Resource and benchmark denominators separate; resolve repeated-person membership before any donor-level training claim. |
| Single-cell/nucleus RNA atlas | 102 analyzed biological donors; 64 donors in the corrected disease-state analysis. | Learn lineage expression distributions and aggregate to donor-level states. Cells can train an expression model but cannot multiply disease-label replication. |
| GSE296875 RNA + ATAC | 39 donors, 68,398 same-nucleus profiles; fibrosis observed for 37 and steatosis for 38. | Main paired anchor for lineage-specific RNA–accessibility learning. Cardiometabolic histopathology is available; adjudicated MASLD/MASH, NAS, ballooning, and lobular inflammation are absent. |
| GSE267145 RNA + H3K27ac | 99 exact participant pairs from different aliquots of the same sample. | Direct bulk RNA–chromatin development experiment with histology. Fractional RNA estimates require continuous-expression treatment. Internal research is activated for specified tasks; the local source-use record currently withholds release of materially trained derivative weights pending review. |
| GSE244832 RNA and ATAC | The source declares 18 matched livers, with different nuclei. | Useful assay-specific MASLD development material. The current plan says the public MM/JB participant crosswalk is unresolved, so it cannot yet provide paired cross-assay validation. |
| GSE281367 ATAC | 12 source-declared individuals. | Assay-specific transport/context; no invented RNA input. Source diagnosis semantics and model-input eligibility have task-specific limits. |
| PXD051911 liver protein | 58 participants with liver proteomes and histology components. | Learn a separate protein-state representation and test molecular correspondence at the gene/state level. These participants are not paired with bulk RNA. |
| Matched liver/plasma protein extension | 41 initial-visit participants; 345 shared proteins; 49 pass the existing adjusted BH association analysis. | A real protein-to-protein measurement bridge. It can inform candidate readouts; it does not create an RNA-to-plasma training edge. The existing associations are development evidence, not independent prediction validation. |
| Genetics and regulatory context | Adopted 462 main SuSiE genes; the new candidate resolves marginal direction for 83 of 925 supported signal pairs. | Signal-specific inherited hypotheses and possible soft priors. These are population-level associations, not donor genotypes or drug-response labels. |
| Spatial RNA and tissue context | Broad program coverage, but restricted donor-level inference; GSE192741 has four inferential donors. | Localize fixed model outputs and test spatial structure where supported. Do not infer a new biopsy's spatial arrangement from its bulk transcriptome. |
| GSE105127 RNA + RRBS | 19 participants, three zones, adjacent sections. | Optional zonation/methylation method stress test after its evaluator and source-use requirements. It cannot carry a replicated MASLD-severity claim. |

Sources: [current numbers](../../docs/RESULTS.md), [model registry](MODELS.md),
[paired multiome record](config/datasets/gse296875.toml),
[paired bulk chromatin record](config/datasets/gse267145_znf469_human_liver.toml),
[GSE244832 boundary](OVERALL_PLAN.md#gse244832-different-aliquot-activation-boundary),
[protein record](config/datasets/pxd051911.toml),
[matched-protein extension](../../docs/technical/continuum_translation_20260906T135538Z/IMPLEMENTATION.md),
and [RRBS record](config/datasets/gse105127_zonated_rna_rrbs.toml).

This supports a model with several connected components. It does **not** support a
fully observed RNA–ATAC–protein–genotype–spatial training cohort. The disconnected
RNA/protein blocks cannot identify person-specific RNA-to-protein covariance from
their marginal distributions. A gene-linked prior can constrain that relationship,
but its predictions would remain assumption-dependent until measured pairs exist.

### Why this is a different scientific question

The current RNA instruments already order fibrosis severity. A small increase in
that correlation would improve an existing function, but it would not establish the
research capability above. The panel's RNA-seq performance also does not establish
qPCR practicality: the current abundance audit reports seven of its twenty genes
at or below detection in liver. Retain it as a comparator, not the required input
panel for the new molecular-state model. See [the registry](MODELS.md).

There is evidence to motivate a cross-assay learning experiment. The saved
[rank-four RNA-to-H3K27ac result](executions/gse267145-crossmodal-ceiling-20260830/results_rna2h3_v2_20260830T170000/rna_conditioned_results.json)
reports aggregate profile-deviance skill 0.13345 against a training-mean profile
across 99 participants, positive in all five folds. The producing code uses a
ratio of summed deviances, so this is not an equal-weight mean of participant
skills or a comparison against every strong baseline. This is an existing
development result, not performance of the proposed
model or proof of new within-stage biology. It makes reduced-rank regression a
mandatory comparator.

The existing fusion experiments found no detectable benefit for the tested
RNA/H3K27ac severity stack. That does not answer whether molecular supervision can
improve prediction of other measured features. Conversely, it does not justify
assuming a new fusion method will work. The summary also records failure of the
RNA → imputed chromatin → severity route. An imputed assay cannot add new
sample-specific observations beyond its input, although auxiliary training could
improve estimation in finite data. See [the recorded experiments](MODEL_RESULTS_SUMMARY.md).

There is a second novelty constraint inside the repository. The September 6
continuum extension already connects within-stage programs, inherited direction,
assay context, and protein readouts. Merely placing those tables behind an RNA
upload would not establish a new predictive model. The additional result must be
**accurate inference on a new donor and a measured improvement in an experimental
selection decision**. The existing 117 programs stay fixed and post hoc; neither
their scores nor their association calls become training targets or selection metrics.

### A concrete model, with a small first implementation

Use a **gene-linked latent factor model with shared and assay-specific components**.
Begin with regularized linear encoders/decoders. The first question is whether the
data support a transferable molecular state, not whether a transformer can fit them.

1. **Learn RNA–chromatin factors on real pairs.** In GSE296875, use verified same-nucleus
   measurements with donor-grouped splits and donor-balanced contribution. In
   GSE267145, use participant pairs and a separate H3K27ac decoder. ATAC accessibility
   and H3K27ac enrichment are different assays, not interchangeable target values.
   Start with four shared factors and a small nested comparison to eight. These are
   proposed design choices, not claims that four biological mechanisms exist.

2. **Keep composition and expression state explicit.** Aggregate single-cell
   representations to donor-by-lineage summaries. The bulk bridge models expected
   RNA abundance as a mixture of lineage expression distributions on the abundance
   scale, before logarithms. Its weights initially represent RNA contributions,
   not literal cell fractions. Variable RNA content and nucleus-versus-bulk sampling
   require sensitivity analyses. An unsampled lineage is unavailable, not inactive.

3. **Preserve assay-specific variation.** Include separate residual components for
   RNA and each chromatin assay, alongside technical covariates. Do not force all
   views into one severity direction. With RNA-only input, marginalize the
   unobserved chromatin-specific component rather than assigning it a confident
   value. With measured chromatin, quantify the departure from the predicted shared
   state. Large departures may be biological or technical; they are not automatically
   evidence of regulatory rewiring.

4. **Use molecular links as qualified constraints.** Begin with gene identity and
   source-supported promoter/peak relationships. A nearby peak is a candidate link,
   not a confirmed enhancer target. An optional genetic prior can shrink or weight
   specific links using the original signal uncertainty, but cannot force the
   disease-expression direction to equal the inherited-risk direction. Compare it
   against an otherwise identical model without genetics and with matched shuffled
   links. Retain genetic information as an interpretation layer if it adds no
   predictive value. Do not validate it by recovering the same COLOC labels supplied
   as inputs.

5. **Extend to protein without inventing paired people.** First fit the protein
   representation independently and compare gene-linked loadings and recorded
   phenotype relationships. A later joint model may use weak, explicitly tested
   constraints between RNA and protein loadings, with assay-private factors. It
   cannot claim a validated RNA-to-protein decoder. The real liver/plasma pairs
   support a separate, strongly regularized readout model with person-held-out
   evaluation. Feature selection must repeat within training folds; the 49 existing
   associations cannot be selected on all 41 participants and then cross-validated
   as if unselected.

6. **Make an unseen-sample transform part of the model.** Freeze the feature mapping,
   encoders, and decoder parameters. A new RNA sample must not require training on
   the user's cohort or observing its withheld assay. Expose a small continuous
   state vector, predicted assay summaries, coverage, and uncertainty. Compare
   molecular neighbors only within supported populations and report reference
   differences. Keep any existing severity output separate and clearly labeled.

Use assay-native observation models: UMI RNA can support count likelihoods; the
fractional bulk estimates cannot. Binary ATAC detection and fragment-count profiles
need different likelihoods and depth handling. H3K27ac can use a count model or a
validated transformed continuous model. LC-MS protein needs log-intensity and
missing-detection treatment, not a UMI count model; RRBS needs methylated and covered
read counts. Balance losses by donors and feature blocks so the largest assay does
not determine the representation merely by its number of cells or regions.

The synthetic pseudobulk bridge is a useful first experiment, but success on sums of
single-nucleus counts is not validation on real bulk biopsies. Lineage-specific
states inferred from bulk remain model-dependent until independently checked.

### What each modality contributes to a research decision

The deliverable should make three decisions inspectable:

| Researcher's decision | Model/evidence contribution | Necessary distinction |
|---|---|---|
| Which samples should enter a mechanistic experiment? | Continuous molecular profiles and uncertainty, compared within recorded histology where available. | A sample selected for a state is not a predicted treatment responder. |
| Which molecular object and cell context should be tested? | Learned RNA–chromatin relationships connected to signal-specific genetics and observed lineage/spatial context. | Regulatory DNA, gene expression, and the mature RNA molecule are different intervention objects; association does not establish the effect of perturbing them. |
| Which measurement could read out that state? | Tissue-protein observability, measured liver/plasma relationships, and assay coverage. | A correlated readout need not be the intervention target; shared program membership does not establish regulation. |

[Illustrative use case, not a result] A laboratory has twenty F2 biopsies. The
model identifies a continuous contrast between samples with stronger stromal
remodeling and samples with stronger hepatocyte metabolic changes. The lab selects
samples spanning that contrast, then measures predefined chromatin or protein
features. The model is useful if that selection captures more of the measured
molecular contrast than selection using fibrosis, the current RNA severity score,
or a simple expression signature. It fails this claim if the contrast only
reconstructs its input genes, sequencing quality, or broad cell composition.

The existing ECM/ductular and GNMT-related examples can explain the eventual report,
but cannot be chosen as the sole success tests after inspecting model results.
The complete frozen program family remains an interpretation surface. Proteomic
and spatial projections of that surface remain separate measured evidence; predicted
ATAC from RNA is not a second independent confirmation of the RNA signal.

### The first experiment that could invalidate the proposal

**Test donor-specific chromatin prediction within histology strata before scaling
the model or adding a new backbone.** Use the existing GSE267145 paired development
lane and an assay-native GSE296875 experiment in parallel as scientific comparisons,
not as a claim that two assays are identical external replications.

For each held participant, hide the entire target assay and predict it from RNA.
All preprocessing, factor fitting, gene/region selection, and tuning occur inside
the training partition. Evaluate both total prediction and the remaining molecular
variation after a training-fitted histology/covariate baseline. Use fibrosis and NAS
only where actually available; NAS is not an inflammation-specific score. For
GSE296875, use its observed pathology fields and prespecified adult-only and
leave-well-out sensitivities. Do not input a histology endpoint to its own prediction.

Mandatory comparators are the existing reduced-rank regression, RNA PCA with the
same decoder capacity, lineage-average/shrunken chromatin, and models using only
histology, age, recorded sex, BMI, assay quality, and estimated composition where
available. Compare composition-only, state-only, and combined models. Composition
may be disease biology, so its contribution should be reported rather than silently
removed.

Evaluate target features across donors, centered against training-derived feature
means. A high correlation across peaks within each sample can largely recover the
average liver accessibility pattern and miss donor differences. Use appropriate
held-assay deviance or normalized prediction error, donor-level uncertainty, and
within-stratum molecular retrieval as a secondary result. Define sparse strata and
their exclusions before examining results; do not search for the grouping that wins.

Test pairing by shuffling training donor links within the permitted histology and
technical strata and rerunning the fitting procedure. On single-cell data, break
donor pairing while preserving lineage. Otherwise, the model can appear successful
by learning broad cell identity. Whole-donor splits keep every assay, nucleus, and
repeat specimen from a person together.

**Suggested practical criterion, to be fixed before execution:** at least 5% lower
held-assay prediction error than the strongest prespecified simple comparator on
the donor-varying target, with an uncertainty interval supporting improvement. This
is a proposed usefulness margin, not an established power calculation. Also report
whether the interval clears the full 5% margin. Assess attainable precision from
existing development data before treating any experiment as decisive. Failure to
reach a powered detection threshold is not evidence of equivalence.

If donor-specific improvement disappears under these controls, stop the large-model
extension. The result would support a reference annotation function, but not the
proposed new molecular measurement capability. The strongest alternative explanation
is that RNA predicts generic tissue identity and composition while the interesting
chromatin differences remain unobservable from RNA.

### How to establish novelty, impact, and practicality

| Dimension | Evidence needed |
|---|---|
| Novelty | Frozen RNA-only inference predicts measured within-stage molecular variation beyond existing RNA scores and matched simple models. Shared/private factorization alone is not novel. |
| Multimodal value | Removing the paired assay supervision reduces held-assay or downstream selection performance. Test intact versus shuffled molecular links separately. |
| Impact | At a fixed sample-selection budget, the model captures more of a predefined measured molecular contrast than severity-matched, random, and RNA-only selection. This retrospective task is measurable now; improved biological discovery still needs subsequent use. |
| Transfer | Results survive cohort-held evaluation where matching input/output assays exist. The current data do not yet supply an unrestricted, untouched second RNA–chromatin cohort for the full claim. |
| Practicality | Frozen inference works on an external matrix without refitting, with realistic missing genes and low coverage. CPU bulk inference is a design target to measure, not a demonstrated runtime. |
| Adoption | Other groups run the released model on their own samples and use its outputs in an analysis or experiment. Public retrospective performance alone cannot establish adoption. |

Predefine the final primary comparison family and use Holm correction across it.
Use BH within explicitly named exploratory biological feature families. Estimate
uncertainty at the donor level and report each study; a pooled interval over many
cells or gene–peak rows cannot certify generalization to another cohort. Seeds
measure optimization variability, not independent replication. Keep all previously
inspected cohorts labeled development for this new objective and preserve genuinely
reserved endpoints.

Prediction uncertainty and usefulness of an additional assay are worth studying
after the core result. With held paired data, test whether predicted uncertainty
identifies samples for which revealing chromatin improves accuracy more than a
random assay allocation. That could help laboratories spend an assay budget, but it
requires a measured gain; uncertainty alone does not prove an optimal experiment.

### Model choices and the novelty boundary in current literature

Use the small factor model as the primary implementation. Treat **MultiVI** as the
native RNA–ATAC comparator and **scGLUE** as the graph-guided integration comparator.
Their existence means that integrating partially paired assays or adding regulatory
links cannot itself be the novelty claim. The local model registry also records
query-inference, missingness, and decoder limitations, so neither is a ready-made
substitute for the proposed donor-level predictor. See
[MultiVI](https://www.nature.com/articles/s41592-023-01909-9),
[GLUE](https://www.nature.com/articles/s41587-022-01284-4), and the
[local model specifications](config/models/multimodal_integration.toml).

**MIDAS** already addresses mosaic single-cell integration and knowledge transfer.
**MOFA+** is a useful development factor-analysis comparator, but the local audit
does not establish an official frozen RNA-to-factor transform for unseen donors.
Do not compare a model trained on the test donor with a genuinely inductive model
without declaring that difference. Sources:
[MIDAS](https://www.nature.com/articles/s41587-023-02040-y) and
[MOFA+](https://pmc.ncbi.nlm.nih.gov/articles/PMC7212577/).

**GET** already learns expression from accessibility and sequence. It is relevant
only to an observed-ATAC regulatory extension after its local artifact/source-use
requirements; it is not the first backbone for an RNA-only research tool. Source:
[GET](https://www.nature.com/articles/s41586-024-08391-z).

MASLD-specific novelty also faces existing molecular-continuum and spatial
multiomics work. The distinguishing result must be a reusable predictive function
on new samples, with independently measured molecular targets and a demonstrated
selection benefit. Do not claim first or unique from this focused literature check.
Sources: [Kamzolas et al.](https://www.nature.com/articles/s42255-026-01543-7) and
[spatial MASLD multiomics](https://www.nature.com/articles/s41588-025-02407-8).

### Scope of the next development step

Concentrate the next effort on the paired RNA–chromatin falsification experiment,
then the RNA/single-cell-to-bulk bridge. Add protein coupling or genetic priors only
with their corresponding ablations; retain spatial and frozen-program context for
interpretation until a separate predictive task is supported. Keep the existing
severity models as useful outputs and comparators.

One material constraint precedes a shareable final checkpoint: the richest paired
bulk chromatin source currently has a model-weight release restriction in the local
source-use record. A successful internal model would therefore not automatically
satisfy the adoption goal. Establish a releasable training subset or resolve that
specific restriction before making this cohort indispensable to the released model.
Do not substitute pseudo-pairs to avoid the limitation.

The current data are sufficient to test the central learning hypothesis and build a
bounded research prototype. They do not yet establish a universal liver simulator,
person-specific missing proteomes, causal perturbation effects, treatment response,
or independent multi-cohort validation of every output. The proposed development
should earn those capabilities individually, while the Resource and future Cas13
paper remain separate.
