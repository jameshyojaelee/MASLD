# Scientific review: correction plan

**Prepared:** 2026-09-23  
**Scope:** the findings from the read-only scientific and code review completed on
2026-09-22. This document specifies corrective analyses and decisive checks. It
does not replace [PAPER](../PAPER.md), [STATUS](../STATUS.md),
[RESULTS](../RESULTS.md), or the [ROADMAP](../ROADMAP.md), and it does not adopt
new results. The current authorities describe the corrected 50-study COLOC
output as adopted and the five-cohort bulk analysis as selected but awaiting a
synchronized paper release. Findings below challenge particular calculations in
those outputs; their current numbers remain identifiable until corrected results
are independently checked and explicitly adopted.

## Order of work and claim boundaries

1. **Resolve source identity first.** Reconcile GSE130970 control diagnoses
   (item 1) and repair the COLOC signal-overlap test (item 2). These affect who
   enters the disease comparison and which genetic signals can enter a gene's
   evidence state. Recompute dependent quantities only after each input is
   settled.
2. **Recompute central numerical claims.** Correct the Figure 2 architecture
   input and physical-cluster inference (items 3–4), the ABF trait split (item
   5), and Figure 5C (items 7–8). Revise the coding/overlap interpretation
   immediately (item 6); its categorical mechanism claim is unsupported even
   before a rerun.
3. **Repair candidates within their own scope.** Correct the cross-assay
   single-nucleus diagnoses, proteome batch model, and RNA–ATAC sensitivity
   null (items 9–11) before treating the associated module counts as evidence.
   Keep the 117 frozen Hotspot programs and weights unchanged. Resolve the
   Figure 7 development and ceiling claims (items 12–14) without reopening
   Currin fold 0.
4. **Reconcile the paper only after corrected analysis.** Produce one
   old-versus-new table for every affected gene membership, adjusted p-value,
   effect, panel symbol, and claim. Rebuild Figures 1–6, the Gene Catalog, and
   portal entries only where their source values change. Figure 7 remains a
   candidate decision. Apply the final checklist in [ROADMAP](../ROADMAP.md)
   before any release decision.

Until the corresponding analysis is corrected, do not use the following as
settled evidence: GSE130970's 23 putative healthy controls; the 1,471-to-7
F0-to-F1 interpretation; SuSiE signal pairs with inadequate shared posterior
mass; the Figure 2 coding-mass permutation q-value; the ABF 758/65 trait split;
the Figure 5C HLA-DQA1 significance symbol and protein–histology p-values;
the candidate cross-assay single-nucleus/proteome sensitivity counts; or the
claim that either Figure 7 model is close to a general label-imposed ceiling.
The directions and sizes of changes not quantified below are **unknown**.

## Shared analysis rules

- Preserve the original source files and named releases. Write any corrected
  analysis to a distinct result location and compare it with the original at
  the same eligible biological units, genes, variants, and hypotheses.
- Record participant/donor or physical-cluster counts, the estimand, effect
  units, adjustment variables, missingness, and the family used for each
  multiple-testing correction. Runs, cells, peaks, and repeated sections are
  not independent donors. Label a nominal p-value as nominal.
- Freeze the corrected source join and any analysis choices before examining
  changed outcomes. Retain negative, indeterminate, untestable, and
  source-dependent states rather than converting an ineligible test into a
  negative result.
- Set and record random seeds for permutation/bootstrap analyses. Record the
  R `sessionInfo()` or corresponding environment inventory beside rerun
  results. Run analytical refits on compute nodes under the existing project
  environments; this plan authorizes no job submission or new dependency.
- For each corrected figure, compare the complete plotted source table with
  the rendered PDF, extract printed text with `pdftotext`, render to PNG and
  inspect the image. Check the complete family, not only changed examples.

## 1. GSE130970 F0 biopsies are mislabeled as controls

**Observed and affected claim.** `RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/00_harmonize_metadata.R:158-163`
sets `condition` and `group_binary` to `Control` whenever fibrosis stage is
zero. The selected five-cohort manifest contains 23 such GSE130970 controls,
although 16 have steatosis above zero and 19 have NAS above zero. In
`figures/candidates/bulk-release-2026-09-09-bulk-r1/analysis/stage_extensions/five_cohort_sample_manifest.tsv:600`,
SRR9036353 is `Control`, F0, NAS 5, and `diagnosis_harmonized=NASH`.
The source SRA row at
`RNA-seq/Human/Patient_Cohorts/archive/old_data_metadata/metadata/GSE130970_SraRunTable.csv:48`
records steatosis 3 and lobular inflammation 2. The pooled disease/control
fit uses `group_binary`; the F0 sensitivity also treats these samples as source
controls. `docs/RESULTS.md:524-558` currently assumes GSE130970 F0 is entirely
control and reads the 1,471-to-7 F0-to-F1 call change as control confounding.

**Correction.** Build a per-biopsy diagnosis table from the study's source
clinical classification and SRA/GEO identifiers. Resolve disagreements
explicitly; do not define a true healthy control by `fibrosis_stage == 0` or
by NAS alone. Retain separate fields for fibrosis grade, steatosis/activity,
source diagnosis, and healthy-control status. Join by stable source IDs, fail
if any of the 78 source biopsies or retained 23 F0 samples has an unresolved
control label, and assert that a steatotic biopsy cannot silently enter the
healthy-control arm. Amend harmonization and the disease-control design, then
refit the selected fragment-count five-cohort contrast on the same QC
substrate. Recompute the conventional canonical DEG family (`padj < 0.05` and
`|log2FC| > 0.50`), the separate TREAT sensitivity, cohort/leave-one-out
contrasts, F0 reference arms, stage summaries that use control status, and
bulk-dependent Gene Catalog fields. Check the actual design rank and
within-cohort overlap of controls and disease before interpreting an adjusted
contrast.

**Decisive checks.** Reproduce a source-ID-by-diagnosis table for every retained
GSE130970 sample; verify SRR9036353 and the other steatotic F0 examples
manually; report old/new control counts by cohort. Refit the complete gene
family with its existing correction, report changed coefficients and DEGs,
then recompute the F0-to-F1 comparison under all-F0, source-disease-only F0,
and a source-control-status term where estimable. Withdraw the existing
1,471-to-7 explanation unless the corrected arms support it. Do not call any
cross-sectional contrast longitudinal progression.

## 2. COLOC signal pairs lose the shared-posterior eligibility check

**Observed and affected claim.** `GWAS/finemapping/src/06_susie_coloc.R:695-709`
reduces both full SuSiE `lbf_variable` matrices to shared SNPs before calling
`coloc.susie`. This makes the package's default 50% posterior-overlap check
appear to pass by construction. In saved EZR–UKBB_ALT fits, two eQTL signals
retain 0.0001226 and 8.97e-28 of their original posterior support yet receive
finite PP.H4 values. Their eligibility is wrong. The number of affected
primary-positive genes has not been measured; EZR's inspected primary positive
signal retains 0.998987 support. The installed `coloc` version used for the
rerun is 5.2.3; its behavior must be verified against that installed source.

**Correction.** Preserve the complete original GWAS and eQTL signal matrices
through the package call, with stable SNP names and the appropriate null
column. For every signal, calculate and export the fraction of its full
posterior mass on truly shared, allele-harmonized SNPs before any restriction.
Require the intended overlap threshold for both signals, and mark failures
`untestable_insufficient_shared_posterior`, not PP4-negative. Keep the
preexisting allele, ancestry, LD-reference, and convergence rules intact.
Use the installed `coloc.susie` behavior or an explicitly equivalent
eligibility calculation; do not renormalize the full posterior on the
intersection before testing support.

**Decisive checks.** Include synthetic fits with high and near-zero shared
support; both must exercise the intended branches. Reproduce the EZR boundary
pair and a passing primary pair. Across all 50-study-by-chromosome tasks,
tabulate eligible/ineligible signal pairs and genes, PP4 changes, threshold
crossings, exact variant anchors, and the adopted 462-gene primary membership.
Rederive the gene-level and trait-level COLOC summaries, genetic/state overlap,
Figure 2 genetics, variant-consistent Figure 5D, worked loci, and Gene Catalog
genetic evidence states from the corrected output. Explicitly report any
change to GNMT, MAT1A, CYP2C19, IGFBP7, or BICC1. Do not infer causality
from COLOC alone.

## 3. Figure 2 architecture includes rejected credible sets

**Observed and affected claim.** `scripts/figures/fig2_multiplicity_correction.py:230-245`
and `fig2_coding_mass_robustness.py:110-112` omit the source
`pip_sum_gate_passed` rule. One direct-trait and 51 enzyme-trait credible sets
that fail the rule enter the coding-mass results; all-missing consequences can
be summed as zero. An independent grouping of the saved source table gives:

| Input | Direct coding mass | Enzyme coding mass | Ratio | Enzyme clusters |
|---|---:|---:|---:|---:|
| Current combined-mass calculation | 22.542908% | 3.442006% | 6.549352 | 266 |
| Passing credible sets only | 22.680461% | 3.389893% | 6.690613 | 259 |

The original Figure 2D descriptive protein-altering segments are a different
quantity and should not be described as having this 22.54% defect.

**Correction.** Apply the passing-credible-set filter before constructing both
consequence and noncoding-context numerators, using the same stable set IDs on
the credible-set and member tables. Require finite numerator and denominator
mass for an eligible set; keep missing or unannotated mass explicit instead of
turning it into a negative observation. Recompute all six contrasts and the 62
recorded robustness fits with the corrected denominator. Preserve the
distinction between protein-altering mass and protein-altering-plus-canonical-
splice mass in panel labels and prose.

**Decisive checks.** Assert exactly 1 direct and 51 enzyme failures in the
reviewed input; independently reproduce the two corrected point estimates
above before rerunning inference. Produce old/new counts, fractions, ratios,
intervals, and six-test family-adjusted q-values. Visually verify every
affected Figure 2 panel and caption against the new source table.

## 4. Figure 2 permutation and bootstrap break physical-cluster dependence

**Observed and affected claim.** `scripts/figures/fig2_multiplicity_correction.py:154-199`
resamples one row per trait class and operational cluster. Nineteen physical
clusters appear in both classes, so the 290 current rows represent 271
physical clusters; after the credible-set filter they represent 264. The
independent row resampling does not preserve shared region context. Thus the
reported family BH q=0.00059994 and the robustness p-values have not been
established under a dependence-preserving null. This does not show that the
descriptive asymmetry disappears.

**Correction.** Define one physical-cluster ID across trait classes before
inferential resampling. Resample a physical cluster with all of its class rows
together. Specify the trait-label permutation's exchangeability assumptions:
shared clusters must retain their paired observations, and study/trait
ascertainment must be handled or shown not to change the null. A paired
within-cluster exchange, a justified stratified randomization, or a clustered
bootstrap interval may be suitable after inspecting the actual incidence
structure; a simple shuffle of 290 rows is not. If no defensible permutation
exists, report the effect and cluster-aware interval without a significance
claim. Retain the six-test family and apply BH to that entire prespecified
family, with a recorded seed.

**Decisive checks.** Reproduce the 19 shared clusters and the eligible-cluster
census; test a synthetic shared-cluster example in which paired rows always
move together. Report corrected confidence intervals, p-values only for a
defensible null, and six-test q-values. Recompute the 62 sensitivity fits
using the same corrected inferential unit.

## 5. The ABF trait split uses a SuSiE driver

**Observed and affected claim.** `GWAS/finemapping/src/perf/rebuild_tier12.py:68-73`
stores the strongest SuSiE study as `driving_gwas` whenever any finite SuSiE
result exists. `docs/RESULTS.md:1340-1359` then calls the 823 = 758 enzyme +
65 direct ABF split each gene's strongest evidence. In a recoverable subset of
772 ABF-positive records, 62 trait classes are wrong: 52 change from enzyme
to direct, 10 from direct to enzyme. GNMT's stored UKBB_ALT SuSiE PP4 is
0.128219, while its best ABF PP4, 0.705832, comes from
`2021_34841290_NAFLD_EUR`. The 823 also includes one blank gene identifier.
The separate 448/462 primary SuSiE split is not implicated by this defect.

**Correction.** Track `best_susie_gwas` and `best_abf_gwas` independently with
their PP4 values and the frozen Tier 1/2 trait registry. Apply the same
predeclared tie rule to each method and exclude blank/unmapped gene IDs from
named-gene counts, retaining their rows separately for diagnosis. Only the
corresponding method-specific driver may define its trait-class breakdown.
Recompute ABF counts after item 2 because signal eligibility or gene membership
may change.

**Decisive checks.** The pre-COLOC-correction reconstruction must identify
the 62 class disagreements and GNMT example; the final reconstruction must
report complete denominators, both method-specific splits, missing IDs, and
any downstream text or Figure 2 changes. A gene with a high ABF PP4 and low
finite SuSiE PP4 is a required regression case.

## 6. Coding enrichment is presented as a proven explanation for gene-set overlap

**Observed and affected claim.** `docs/RESULTS.md:1362-1369` says
protein-altering variants do not generate expression QTLs and makes this the
reason genetic and disease-state gene sets overlap at chance. Coding and
splice variants can alter RNA abundance through nonsense-mediated decay or
splicing; see [Rivas et al. 2015](https://pubmed.ncbi.nlm.nih.gov/25954003/).
The aggregate consequence fractions do not test that mechanism or distinguish
it from expression observability, phenotype choice, and study power.

**Correction.** Remove the categorical mechanistic explanation from current
claims. State the descriptive result in its own units: the relative coding PIP
mass differs between direct and liver-enzyme trait sets under the corrected
analysis; the maps' overlap is a separately estimated gene-level quantity.
It is reasonable to say that some coding effects may act through protein
function without changing steady-state RNA, while others can change RNA.
Present any proposed role in the observed overlap as a hypothesis, not an
established result. Do not start new mechanism searches just to preserve the
current sentence; a formal conditional comparison would be a separately
specified analysis with common gene eligibility, phenotype, power, and
observability controls.

**Decisive checks.** Read the abstract, Results, Figure 1/2 captions, and
Gene Catalog explanatory text for variants of “cannot generate an eQTL” or
“explains overlap.” The paper should make no causal attribution that the
current comparison did not estimate.

## 7. Figure 5C marks HLA-DQA1 mRNA nonsignificant

**Observed and affected claim.** The rendered
`figures/main/fig5_molecular_context/panels/fig5c_mrna_protein_composite.pdf`
has a hollow HLA-DQA1 mRNA marker. For the selected ENSG00000206305.12 row,
`canonical_deg_results.csv` gives log2FC 0.651447, conventional BH p 0.016089,
and TREAT q 0.425763. The current canonical definition is conventional
`padj < 0.05` with `|log2FC| > 0.50` (`docs/STATUS.md:261-267`).
`Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R:242-262`
renames `treat_fdr` to `bulk_padj`, and
`scripts/figures/composite_mrna_protein.R:65-69` passes it into the canonical
predicate. In this fixed 25-gene display, one mRNA significance marker is
known to change.

**Correction.** Export conventional `padj` and `treat_fdr` as separate,
unambiguous fields for the fixed selected gene ID and symbol; use the shared
canonical predicate on the conventional field. Verify the selected Ensembl
row because HLA-DQA1 has multiple rows in the bulk source. Keep TREAT only as
the explicitly labeled sensitivity; do not reselection the 25 protein rows.

**Decisive checks.** Recalculate all 25 mRNA marker states from the source
gene IDs and compare against the exported panel table. Confirm HLA-DQA1 is
filled in the new rendered panel and that its printed effect remains 0.651.
Extract PDF text and inspect the rasterized panel before replacing any
reader-facing Figure 5C artifact.

## 8. Figure 5C partial protein–histology associations use mismatched samples

**Observed and affected claim.** `Analysis/Multimodal_Program_Projection/scripts/02_proteomics_projection.R:280-316`
residualizes histology ranks over all 58 participants but protein ranks over
each protein's observed participants, then correlates the two residual sets.
HLA-DQA1, P4HA2, and ABHD5 have 40, 43, and 56 measured participants,
respectively, so their 15 displayed wedges are not a same-sample partial
Spearman estimate. The ordinary `cor.test` on residuals also uses n−2 degrees
of freedom despite four fitted nuisance terms: acquisition batch, age, BMI,
and sex. With saved rho held fixed and n−6 degrees of freedom, BH-significant
associations across the 125 tests decrease from 113 to 111. ABHD5–ballooning
changes from q=0.04303 to 0.05199 and CD74–steatosis from q=0.04828 to
0.05749. That calculation is a degrees-of-freedom diagnostic, not the fully
corrected refit.

**Correction.** For each of the fixed 25 proteins and five histology
features, form one complete participant set for protein, histology, and
covariates. Rank each variable within that same set, residualize both against
the identical full-rank nuisance design, then compute the partial Spearman
effect and its chosen uncertainty method. Use appropriate residual degrees
of freedom for a parametric approximation or a covariate-aware permutation
with participant exchangeability justified. Record each pair's n and missing
reason; adjust all 125 p-values together by BH.

**Decisive checks.** Recompute all 125 rho/p/q rows; verify that both rank
fits and the test use the same participant IDs. Report old/new rho and q for
the 15 missing-protein cases and the two boundary significance calls above.
Regenerate the Figure 5C wedges and inspect the PDF and PNG.

## 9. Cross-assay single-nucleus candidate uses nine wrong diagnoses

**Observed and affected claim.** `scripts/analysis/cross_assay_modules/08_score_snrna.R:51-58`
gets donor diagnoses from an old run-level metadata table. Among the 70 modeled
donors, P6, P51, P55, P56, and P84 are source NAFLD but are scored as
Steatohepatitis; P30, P98, PHL1, and PHL2 are source Healthy controls but
are scored as Steatosis. The source IDs are in
`data/GSE202379/metadata/donor_pairing.csv`; corrected donor metadata already
exist at
`Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/metadata/donor_metadata_stage_corrected_donor.tsv`.
The candidate reports six all-cell replications and 27/1,280 lineage calls.
Which calls change has not been determined. This is a 160-module cross-assay
candidate; the separate 117 Hotspot programs are fixed and unaffected by this
specific scoring error.

**Correction.** Join each collapsed donor to the corrected source-diagnosis
table by stable dataset and donor ID; reject conflicts, ties, and unpaired
runs rather than majority-voting historical run diagnoses. Keep the coarse
Healthy/Steatosis/Steatohepatitis/Cirrhosis endpoint explicitly distinct from
fibrosis grade. Preserve the donor as the inferential unit and the frozen
160-module gene memberships. Refit all-cell and all lineage tests, then
recompute the full all-cell and 1,280 module-by-lineage correction families
and downstream evidence states.

**Decisive checks.** Assert the nine corrected donor states against the source
table and report 197 runs collapsed to 70 unique donors. Report old/new
effects, q-values, and counts; do not simply relabel unchanged fitted results.
Check that no run or cell is counted as an independent patient.

## 10. Cross-assay proteome candidate omits acquisition batch

**Observed and affected claim.** `scripts/analysis/cross_assay_modules/07_score_proteome.R:101-123`
fits age, BMI, and sex with or without intensity descriptors but omits the
known 2019/2020 acquisition batch. The native proteomics model includes this
factor. The batches contain 14/33 versus 3/25 MASH participants and mean NAS
3.24 versus 2.00; the existing protein PC1 batch diagnostic shows substantial
separation. The candidate currently reports 31 of 34 nominally surviving
proteome associations after complexity adjustment. [Inference] Some may
reflect batch-associated disease differences; the changed count is unknown.

**Correction.** Include acquisition batch in the primary participant-level
NAS/fibrosis models and in the complexity sensitivity. Verify that both
outcome categories occur within each batch and that the design is full rank.
Keep the prespecified covariates, coverage rule, 160-module membership,
measurement units, and assay-specific family correction. Refit the relevant
competitive controls with the same batch adjustment where they are used for
claims; do not treat intensity descriptors as a proven substitute for batch.

**Decisive checks.** Publish a batch-by-diagnosis and batch-by-stage table,
complete-case n, design rank, old/new coefficients and intervals, and the
complete 160-test complexity-adjusted BH family. Reevaluate proteome evidence
states and figures only after the refit. If batch and phenotype cannot be
separated for a contrast, mark it indeterminate rather than significant or
negative.

## 11. RNA–ATAC sensitivity null does not preserve added nuisance relationships

**Observed and affected claim.** `scripts/analysis/cross_assay_modules/11_paired_transfer.R:202-214`
adjusts RNA and ATAC scores for well plus fibrosis or hepatocyte fraction,
then permutes raw RNA scores within wells. This changes RNA's relationship
with the additional covariate. The resulting null does not generally test
residual coupling conditional on that covariate. The candidate statements that
33 associations persist after fibrosis adjustment and 19 after hepatocyte
fraction adjustment therefore have uncalibrated p-values. Observed partial
correlations and the separate well-only primary test are not changed by this
diagnosis.

**Correction.** Keep the same 39 paired donors and the same well blocks. For
each sensitivity, fit the reduced RNA nuisance model, permute its residuals
within justified exchangeability blocks, reconstruct the score, and reapply
the full partial-correlation statistic against ATAC. Use the same scheme for
all eligible modules with fixed seeds; verify the treatment of rank transforms
and discrete fibrosis ties. The script's H3K27ac arm provides an existing
nuisance-aware permutation pattern. Apply BH over the full prespecified
sensitivity hypothesis family, not only modules selected by the primary test.

**Decisive checks.** A synthetic example with RNA–fibrosis association but no
conditional RNA–ATAC association should have calibrated null p-values.
Recompute both sensitivity families and their old/new q-values, report the
number retained and uncertainty across donors, and keep the well-only primary
33/150 result separate.

## 12. Adapter development selection is not nested across folds 1–4

**Observed and affected claim.** `scripts/analysis/alphagenome_campaign/model_nested_selection.py:235-258`
chooses an outer fold's recipe from fits validated on the other development
folds. Those fits were trained on the outer fold: a saved fold-2 split lists
training folds 1, 3, and 4. Thus the outer fold can influence its own recipe
choice indirectly. The development result is not an unbiased nested
selection estimate. Fold 0 was held out of this process and has already been
spent on its recorded confirmation; this finding does not undo that distinct
comparison.

**Correction.** Choose one of two scientifically honest dispositions. The
low-compute option is to relabel the four-fold selected-recipe results as
exploratory development estimates and base any prespecified confirmation claim
on the completed fold-0 contrast. If an unbiased development estimate is
needed, retrain each inner candidate for every outer fold using only the
remaining development folds, excluding both the outer and the inner validation
fold and always excluding fold 0. Perform recipe selection entirely inside
the outer training subset, then score the untouched outer fold once. Preserve
the block-level uncertainty unit and the candidate family correction.

**Decisive checks.** For every predicted row, assert that its block/fold is
absent from all training fits that selected its recipe. Publish the complete
outer/inner split table and revised mean MSE contrast with intervals. Do not
reuse fold 0 for tuning, retesting, or new confirmation.

## 13. The caQTL “LD-tagging bound” is a restricted simulation

**Observed and affected claim.** `Analysis/MASLD_Model_Benchmark/executions/fig7-step0-ceilings-20260922T175014Z/0b_caqtl_ceiling/eb_mixture.py:330-350`
samples whether each lead is causal using its PIP, gives selected leads an
oracle effect, and predicts every unselected lead with a fold constant.
Noncausal lead variants can have predictable marginal association through LD
with nearby causal variants; PIP alone does not make their marginal effects
unpredictable. The simulated 1-kb score of 0.807 is therefore conditional on
that information restriction, not a demonstrated ceiling for every sequence
model. The 0.7441/0.85 comparison is on 3,627 matched leads and uses their
0.8774 simulated value, not the 32,322-lead 0.807 value.

**Correction.** Rename the quantity as a restricted causal-lead oracle
sensitivity and state its causal-mask and constant-prediction assumptions. Use
the same 3,627-lead denominator for model-to-simulation comparisons, and show
the full 32,322-lead result separately. Remove claims that AlphaGenome is near
a general label-imposed limit or that remaining capacity gains are bounded by
0.13 Spearman. A genuine task-wide bound would require a defensible joint
model of causal effects, LD, observed marginal betas, and accessible sequence
information, or independent repeated labels; neither follows from the present
PIP draws.

**Decisive checks.** Verify every model/denominator pair and uncertainty
interval in Figure 7 and its text. If only wording changes, retain the achieved
model scores and show the 0.807/0.8774 simulations explicitly as conditional
sensitivity analyses.

## 14. Severity “label ceiling” does not identify measurement error alone

**Observed and affected claim.** The severity simulation in
`Analysis/MASLD_Model_Benchmark/executions/fig7-step0-ceilings-20260922T175014Z/0a_severity_ceiling/`
uses repeat biopsies separated by time to calibrate an independent label-error
term. These paired differences can include biological change and biopsy
sampling as well as reread error. `RESULTS.md:1928-1939` claims the released
instrument achieves 0.82–0.95 of its label ceiling and that GSE276114 has
almost no headroom. Even under the simulation's own assumptions, those ratios
use the higher error source (0.6433 stages); its technical report identifies
0.4655 as the conservative choice for an “approaches ceiling” claim. Saved
achieved/ceiling ratios for the latter are 0.757, 0.731, and 0.885 across the
three v1.2 cohorts.

**Correction.** Describe both paired sources and their transport assumptions
as sensitivity scenarios, not direct estimates of single-occasion reread
noise. Present observed Spearman, recorded-stage tie ceiling, each simulated
ceiling, and the matched achieved/ceiling ratio with separate denominators.
Remove “almost no headroom” unless an independently justified measurement
model supports it. Contemporaneous independent pathology rereads would be
needed to estimate reader error directly; serial biopsies answer a different
question. Do not use RNA-sample biopsy variability as though it were measured
by the current label data.

**Decisive checks.** Reproduce the three ratios from
`0a_severity_ceiling/out/label_ceilings.tsv`, with uncertainty and the source
sample sizes (58 and 30 repeat pairs). Verify that Figure 7 and its caption
name the assumptions next to each number. No rerun of the released severity
instrument is implied by a correction to this ceiling interpretation.

## Completion evidence for this correction plan

The correction is complete only when the affected families have new source
tables and independent numerical checks; old and new gene/variant membership,
effect sizes, q-values, evidence states, and printed panel values are compared;
and the applicable manuscript claims are revised or explicitly withheld. A
new result may be null, weaker, or untestable. Report all three outcomes
without replacing an original threshold or searching for a new example.

The high-priority downstream order is: corrected GSE130970 source labels to
bulk and stage claims; corrected COLOC signal eligibility to genetic members;
then Figure 2 inference, exact locus examples, the Gene Catalog, and dependent
Figures 1/5/6. Figure 5C can be corrected independently. Candidate cross-assay
modules and Figure 7 can be reviewed in parallel because they do not redefine
the Resource's 117 programs, genetic evidence classes, or separate follow-up paper. After
paper, figure, release, Catalog, or portal authorities change, run
`python3 scripts/manuscript/validate_resource_scope.py` and the final ROADMAP
checks. Final adoption remains a separate explicit user decision.
