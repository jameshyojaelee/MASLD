# Resource paper roadmap

Use [STATUS](STATUS.md) for current completion and [RESULTS](RESULTS.md) for
exact values. Dated execution notes below are records of their stated scope;
this documentation cleanup authorizes no new analysis or submission.

| Planning question | Section |
|---|---|
| Latest genetic-identity and transportability follow-up | [November 7 priorities](#variant-identity-and-transportability-follow-up-2026-11-07) |
| Prespecified scientific comparisons | [November 6 analysis plan](#scientific-analysis-plan-2026-11-06) |
| Dependencies and figure completion | [Critical path](#critical-path) · [Figure checklist](#figure-completion-checklist) |
| Final manuscript release requirements | [Final checklist](#final-release-checklist) |
| Closed or future analyses | [Closed analyses](#analyses-that-are-closed-for-this-paper) |
| Find the existing files | [Analysis routes](technical/DATASETS_AND_PIPELINES.md) · [Hypothesis index](../Analysis/Hypotheses/README.md) |

## Manuscript rewrite with the novelty analyses, 2026-11-12 (user-approved)

The ENZFAT1, DIR6, MQTL2, PROT5, TRAIL4, CCL1, JOIN7/TRIAL1, AG2X Stage A/B and GSET8
results are used in the manuscript text by user approval, 2026-11-12; formal adoption
under the [final release checklist](#final-release-checklist) is pending. Scope is the
2026-11-12 block at the top of [PAPER](PAPER.md); changed files are listed in
[STATUS](STATUS.md). Remaining items:

1. Place the Figure 6G candidate panel
   (`figures/candidates/fig6g-trial-mechanisms-20261112T140124Z/`; highest trial phase,
   not approval, and four endpoint columns) and build any new Figure 2 panels as
   individual PDFs with numbers in CAPTION.md; then assign panel letters in the text.
2. Write the Figure 5 Complex I and PLIN2 sentence from DIG1 values after they are
   recorded in RESULTS.
3. Bring `docs/manuscript/draft/*.md` in line with the LaTeX text, or mark it historical.
4. Run the final checklist before any of these numbers counts as adopted.

## Review-round restructure, 2026-11-09 (user-approved plan)

The user approved the review plan and asked for it to be carried through.
Scope authority is the 2026-11-09 block at the top of [PAPER](PAPER.md).
Work in order; each item keeps the final checklist below.

1. Step 0 corrections: rebuild the 90 duplicate-ID eQTL fits
   (`Analysis/Hypotheses/EQTL-DEDUP-*`) and finish the unique-nucleus recount
   (`Analysis/Hypotheses/NUCLEUS-FIX-20261107/`).
2. Figure 2: add the enzyme-only versus disease-liability panel and its
   out-of-sample check; replace the RORA and FABP1 loci with TNFSF10; move the
   coding-fraction contrast to the supplement.
3. Figure 3: add the SteatoSITE replication scatter, the composition adjustment
   (slope 1.462 to 0.256) and the within-hepatocyte partition; move the overlap
   null to the supplement.
4. Figure 4: lead with SteatoSITE programs (98/112) and their
   composition-adjusted effects, then within-hepatocyte results; move the
   0/117 program test, continuum and 43 systems to the supplement.
5. Figure 6: add coding fine-mapping and disease-liability fields, a known-gene
   panel and the 26 shared genes; rename `concordant`.
6. Figure 7 to the supplement; fix the text and panel mismatches listed in the
   plan (Part 1D).

Candidate panels and tables go to new timestamped directories; nothing in
`figures/main` or a named release is overwritten.

## Variant identity and transportability follow-up, 2026-11-07

The all-eQTL-fit check and manuscript transportability revision requested by
the user are complete. Ninety primary genetic genes have duplicated eQTL
variant identities; MPRIP additionally has contradictory MVP allele metadata
and invalid original reference LD. Do not infer a validated 346-gene core by
subtraction or reclassify the Catalog from one variant-deletion sensitivity.

1. **Resolve genetic events before further recurrence claims.** Start with
   MPRIP and DIP2B, whose duplicated positions enter supported eQTL credible
   sets, then the remaining 88 affected primary genes and 47 additional
   full-portfolio genes. Normalize reference-anchored alleles on the stated
   assembly; require unique mappings for both assays. MPRIP requires clarification
   of the original MVP repeat-event coding, not a choice based on H4 or p value.
   Reject ambiguous rows explicitly. Rebuild both full assay fits, not columns
   of an existing posterior. Include all variant-coverage losses in the report.
2. **Use valid, justified LD.** The six original-reference MPRIP arms fail the
   unchanged matrix criterion. The local genotype-derived sensitivity uses
   one EUR founder per pedigree family, common-sample Gram LD and all three
   predeclared event interpretations. It cannot establish source-matched MVP
   inference. All three GWAS fits converge on genotype-derived LD, but both
   eQTL fits fail the fixed 500-iteration limit, including repeat exclusion; no
   H4 is qualified. Obtain study-matched or ancestry-justified LD before asserting a
   mechanism; compare fixed control loci and both-assay LD diagnostics. Preserve
   L=10, priors, seed, convergence, full-fit shared-posterior checks and explicit
   untestable states. No full 50-study refit has been submitted in this task.
3. **Re-evaluate matched signals and dependent claims together.** Match GWAS
   credible sets and physical variant events across the fixed independent
   same-trait comparisons after correcting identities. Reusing an eQTL effect
   index is insufficient. Recalculate genetic membership, overlap, accessibility
   links and affected Catalog entries from complete corrected outputs before
   requesting adoption of any new gene list. Retain original counts as history.
4. **Strengthen tissue transportability without outcome-driven selection.**
   The four/five-cohort REML sensitivity is now in the manuscript. Next use
   RNA-count-native cohort models with their own filtering and donor/cluster
   units to assess the moderated-SE approximation, retaining the 24,115-gene
   BH family and the fixed 1,854-gene reporting set. Prespecify stage, specimen,
   control definition and ancestry checks before examining new cohort outcomes.
   A new confirmatory disease/control analysis requires an exposure-qualified
   external cohort; SteatoSITE's stage comparison cannot substitute for it.

Exact results and limits are in [RESULTS](RESULTS.md#variant-identity-and-transportability-2026-11-07).
Detailed commands and checks are under
`Analysis/Hypotheses/COLOC-IDENTITY-20261107T150017Z/`.

**Required count correction authorized, 2026-11-07.** Finish joint counting of
the six remaining repeated RNA captures, exclude the four source-confirmed
antibody-tag runs from RNA processing, verify unique barcode/molecule
inputs, and refit donor-level gene contrasts and composition. Assess all 117
fixed program projections with unchanged memberships and weights; retain the
distinction from historical Hotspot neighborhood/PCA scores. Reconcile affected
working figures and text only against verified corrected results. The model
registry now selects an explicit-assumption histology-error calculator and the
remaining unsupported cohort-disagreement sentence is removed. Source records,
commands and results belong under `Analysis/Hypotheses/NUCLEUS-FIX-20261107/`.

**Catalog adoption and priority analyses authorized, 2026-11-06.** The user
adopted the 188 priority-rule genetic-state corrections and 180 class changes.
The selected Catalog is `gene-catalog-eligibility-20261106`: 400 genetically
anchored, 26 both-supported, 1,118 state-associated, 12,958 unresolved and
13,731 untested genes. Figures 1B/3G use the matching 3,589-gene native overlap
subset, distinct from the 14,502 jointly eligible Catalog identities.
Work and checks are recorded in
`Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/`.

Priorities 1 and 3 below now start with a five-cohort transportability analysis,
a saved-signal matching inventory and one LD diagnostic locus with two controls.
All new jobs use `cpu`, account/QOS `nslab`, and 90-hour limits; requested
resources are respectively 2 CPUs/16 GB, 1 CPU/4 GB and 4 CPUs/64 GB. The LD
pilot follows input checks. This starts the highest-priority research without
changing the fixed 1,854-gene signature, 436 genetic positives, 117 programs or
43 systems. New research estimates remain candidates until evaluated; the
existing Catalog correction is adopted separately. Priorities 2, 4 and 5 remain
planned and depend on the specimen, exposure and assay checks below.

**Review corrections, 2026-11-06.** The authorized focused reruns and independent
checks are complete; working text and affected artwork use the corrected values
in [RESULTS](RESULTS.md#review-corrections-2026-11-06). The
[scientific analysis plan](#scientific-analysis-plan-2026-11-06) below specifies
five prioritized questions, inputs, inferential units, comparisons, correction
families and stopping rules. The subsequent user request above authorizes the
highest-priority analyses and adopts the eligibility-corrected Catalog.

**Main Figure 4 completed at the user's request, 2026-11-03.** Individual A–H
PDFs and the compact 468 × 572 pt assembly use embedded 6 pt Helvetica with
minimal labels and no footnotes. C uses the October 26 corrected 117-program
family (0 selected, minimum q=0.18135). The Figure 4 legend and Results text
include the November 2 discovery-dependence and matched random-axis caveats.
No program, system or continuum definition changed. Incorporating the
remaining source-dependent supplementary artwork and other October 26 outputs
is still required before a synchronized release; named r9 is preserved.

**Single-cell source-label correction, 2026-10-26: completed candidate.**
Gene contrasts, composition, communication, content-lineage contrasts, target
models and source-dependent hypothesis analyses now use the corrected RNA
diagnoses. Program permutations, QC/ambient/T-cell sensitivities, numeric-SAF
models and paired donor bootstraps are complete. Independent checks reproduce
donor membership, fixed weights, FDR families and intervals; five candidate
PDF panels were visually checked. The candidate includes matching manuscript
text, Methods and legends under
`Analysis/SingleCell/candidates/source-diagnosis-correction-20261026/`.
Historical label-derived portal scores remain withdrawn. Remaining action:
obtain explicit approval for manuscript adoption using the complete matching
candidate and the final checklist below. r9 is preserved.

**Main-panel corrections completed, 2026-10-14.** Corrected bulk-dependent
main and supplementary panels, Figure 4A, the restored five-cohort Figure 3C,
774-unit Figure 5B/C and class comparisons, all 66 fixed twelve-trait genetic
correlations and the
N4 eligibility correction are complete. Verified COLOC r2 is selected for
the genetics panels. The selected r8 Catalog has completed local review,
a molecular-object companion and all five individual Figure 6 panels.
Figure 3/4/6 vector review assemblies and the current Figure 5 assembly are
installed. The source-diagnosis H17 and N3 refits completed on 2026-10-14
(H17: no shared lipid-loading response detected, and the in vivo arm could not
detect canonical lipid-droplet genes; N3: no qualified corrected calls).
Remaining publication work includes final publication layout and the overall
checklist below; unavailable donor/histology joins remain explicit scientific
limits.

**Additional source-diagnosis correction, 2026-10-14.** The two GSE174748
disease explants are source-documented Cirrhosis, so the unchanged three-stage
program model excludes them. The fixed-score refit has maximum n=62
(BICC1 alone selected, zero HC3-selected). H59 refits the 60-donor ambient
common universe with per-dataset variances and selects zero programs, so the
62-donor call becomes a sensitivity. Figure 4C and all 117 supplementary effects are corrected.
Liver Atlas's exact 19-donor biopsy join retains its stage exclusion and
explicitly represents unresolved fibrosis categories. Complete communication
families are refitted: fixed thirteen pairs include six testable and seven
untestable; the fixed 26-program-pair screen has 22 estimable and two supported
signed covariations. Older stage q values remain historical.

**Figure 2 coordinate correction adopted by the user, 2026-10-13.** The
corrected coding/splice estimate, its bootstrap intervals, all 57 robustness
comparisons and the corrected descriptive context/ascertainment results are
selected for Figure 2. This resolves the Figure 2 coordinate-dependent analysis
work. Figure 2D and its unassembled noncoding-context companion are regenerated
and installed with matching source tables and captions. This does not adopt
COLOC r2 at that date or close the synchronized Resource checklist.
COLOC r2 is selected separately on October 14. Exact estimates are in
[RESULTS](RESULTS.md#adopted-corrected-coordinate-results-2026-10-13).

**Bulk refit adopted by the user, 2026-10-11.** The five-cohort fit on 774
genetic-identity clusters (138 control, 636 disease; 24,115 genes; 1,854
canonical calls; `bulk-corrected-controls-20261008T162308Z`) is the bulk
analysis of the paper; the 855 held annotations are reported as not testable.
Numbers and what is not yet rebuilt on it are in [STATUS](STATUS.md) and
[RESULTS](RESULTS.md). Stage, composition, continuum and Figure 3/4/S3/S4
panels were installed on October 14 from their exact candidate sources under
`figures/candidates/`. Catalog review and its companion, Figure 3F,
composition-adjusted models and r2-dependent panel corrections are complete.
Paragraphs below dated before 2026-10-11 that
describe the refit as waiting are history.

**Historical Figure 3F artwork completed, 2026-10-05, superseded October 14:** the 61-gene matrix was rendered in
the existing panel sequence. Remaining scientific work is the already-required
counting and identity correction followed by donor-stage fitting. The PIP
column was removed; 3F displays expression effects and exact-pair PP.H4. Rebuilding artwork does
not adopt the synchronized scientific release.

**Historical Figures 3–4 identity qualification, 2026-09-29:** the genotype ledger identifies
repeated bulk libraries and conflicting source labels; GSE193066 identity calls
remain unresolved. Historical bulk roster counts and newly reconstructed fits
are conditional, not verified independent-donor inference. The correction
candidate includes an explicit roster join and a separate provisional identity
sensitivity. Do not adopt its bulk/network significance claims until identity
is resolved; see [the disposition](technical/FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md).

**Historical Figures 3–4 correction preparation (2026-09-29).** The candidate analyses,
source-category displays, cohort/omission effects and revised scientific claims
are documented in the [claim disposition](technical/FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md).
Review the successful candidate PDFs and source tables in their existing panel
sequence. Unavailable independent validation narrows the claim; it is not a
reason to search repeatedly for a positive cohort. The remaining adoption step
is the existing final checklist and explicit user decision. Do not overwrite
adopted panels, alter fixed memberships or substitute competitive outlines for
maxT while retaining a maxT caption.

**Historical implementation order selected 2026-10-05, completed dependencies updated above.** Resolve analysis-specific identity
and full ALT/scaffold annotation correspondence in parallel with the two
Figure 7 review sets. Exclude unresolved identity/source-category components
from primary inference. The bulk refit waits for complete annotation and
fragment-count evidence; stage, composition, continuum inference and final
Catalog classes follow that fit. Recompute maxT calls rather than retaining 18
as an expected result. Review the provisional Figure 6A/C/E and refresh 6B/D before
finalizing the changed Catalog. Preserve named releases and fixed definitions.
The resumed checks reproduce the identity candidates and both six-panel
Figure 7 review sets. The remaining first dependency is the 855 unresolved
annotation correspondences, followed by complete fragment counting. Do not
start the primary bulk refit or finalize bulk-dependent Catalog classes from
the completed identity screen alone. Genetic clusters remain distinct from
verified independent participants.
The completed 774/738-unit fits use historical unique-only counts and remain
conditional identity sensitivities. Five cohort pilots verify native count
reproduction but do not replace complete fragment counting. Reference exclusions
for seven contamination-scaffold records and 56 complete false-duplication
overlaps were authorized and applied. Their records remain non-applicable;
the other outstanding gene identities still require source evidence.
The revised NMF projection protocol is prepared with unchanged fixed weights,
an explicit unchanged-object mask and one historical normalization reference.
After the corrected recount, validate transport before rebuilding NMF-dependent
panels. Historical normalization checks support eight stage axes and all ten
continuum axes; preserve stage k4-P2 and k6-P6 as unavailable under the existing
criteria. Historical scores are not corrected projections. The revised
projection is a sensitivity definition for review, not a replacement of the
original axes. Do not duplicate merged counts across old feature rows.
The [existing correction outcomes](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md#resource-corrections-implementation-2026-10-05)
contain the implementation log; final scientific adoption remains a separate
user decision.

**Updated:** 2026-10-14 (remaining scientific corrections and current artwork reconciled; final publication checklist retained).  
**Release goal:** finish one internally consistent Resource release. The
research program below is a separate scientific planning decision; manuscript
timing does not set its scope.

**Figure 1–6 framework:** every main figure must expose the claim, applicable
assay and biological unit, coverage/join gate, machine-readable evidence state,
unresolved alternative, and next discriminating experiment. No panel may imply
that modality count is confidence.

## Scientific analysis plan, 2026-11-06

The user requested correction of the confirmed review findings, verification of
affected results, and a concrete plan for the strongest scientific improvements.
The correction runs are complete. The subsequent user request authorizes
starting priorities 1 and 3 as specified above; the remaining comparisons are
planned. Work stays within existing public human data. The 117
programs, their weights, the 43 system memberships and the separate follow-up
study boundary remain fixed. Named release r9 is preserved.

### Immediate correction and verification

| Work | Changed quantity | Required check before using the corrected result |
|---|---|---|
| Lineage-specific pseudobulk eligibility | Donors, tested genes and gene-DE calls | Count cells from eligible libraries actually contributing to each lineage; require at least 50 per donor-lineage; compare full and P30-excluded fits, stage fits and dependent comparisons with their historical results |
| Xenium disease comparison | Permutation p values and uncertainty after carrier selection | Repeat all label-dependent eligibility and carrier selection under each label split and bootstrap; retain the complete 14- and 12-gene Holm families; healthy-liver location remains descriptive |
| Exon/intron ratio | Adjusted ratio effects, BH calls and highlighted genes | Include batch and inferred sex; reproduce the saved adjusted sensitivity and recompute its permutation reference; never infer negligible RNA-stability effects from a nonsignificant ratio |
| Genetic and Catalog qualifications | LD sensitivity, testability and method-specific posterior fields | Reproduce the 14-gene LD diagnostic and separate ABF/SuSiE values; adopt the user-approved priority-rule correction (188 states, 180 classes), preserve historical tables and leave the any-eligible-study alternative separate |
| Manuscript and artwork | Evidence states and printed values | Check corrected source tables, rendered panels and legends; distinguish an underpowered comparison from a tested negative and an MDE from a confidence bound |

Commands, input hashes, environments, result comparisons and visual checks belong
in `Analysis/Hypotheses/REVIEW-FIX-20261107T020001Z/`. Current numerical outcomes
belong in RESULTS and completion/remaining limits in STATUS. New scientific
memberships require review of the complete matching evidence before final release.

### 1. Quantify reproducibility of tissue remodeling across cohorts

**Question and estimand.** Which fixed disease-state effects transfer across
cohorts, and which fibrosis-stage effects transfer to SteatoSITE? These are two
different contrasts. Start from the 1,854 disease-state calls, complete
24,115-gene bulk family, existing per-cohort effects and the 493 staged
SteatoSITE donors. Discovery units remain genetic-identity clusters; SteatoSITE
uses source-qualified donors after its identity screen. Existing SteatoSITE
results are already observed and cannot become a new blind confirmation.

**Analysis.** First assemble the already completed cohort-effect classifications,
heterogeneity, ancestry sensitivity and external comparisons into one gene-level
effect table. Report effect and SE by cohort, the pooled estimand, and which
cohorts contribute information. Fit a prespecified random-effects sensitivity
with restricted maximum likelihood and a small-study interval; report a
prediction interval only where the number of cohorts and estimated variance
support interpretation. Do not select the model that maximizes calls. Compare
the existing pooled estimates with this sensitivity without redefining the
1,854-gene list. Treat missing covariates and nonoverlapping control populations
as limits, not variables recoverable by imputation from the outcome.

**Validation and decisions.** Retain the original SteatoSITE grouped-stage family
and its Holm correction. For a new disease/control confirmation, choose one
source-qualified public cohort after an exposure and identity inventory but
before reading its relevant expression outcomes. Freeze normalization,
eligibility, contrast and a gene-set slope statistic; test it by donor-label
permutation within valid source strata. A positive global slope and its
confidence interval support transfer of the fixed signature; gene-level
replication requires directional tests with BH over all eligible fixed genes.
Do not infer equivalence from an interval spanning zero. If no unexposed eligible
cohort exists, report retrospective transportability and stop the confirmatory
arm. Deliver a cohort-effect table, a compact forest/replication figure and
Catalog transfer fields. **Priority:** first biological synthesis; it reuses
completed analyses before requesting new computation.

**Started, 2026-11-06.** The five source-model reconstructions reproduce their
effects and moderated SEs to 2.7e-14 in 774 unique genetic-identity clusters.
The long table has 120,575 gene-by-cohort slots: 112,052 eligible and 8,523
untestable in the native GSE126848 gene universe. Job `22811742` fits REML
with modified Hartung–Knapp intervals (variance factor floored at one), using
the reconstructed limma moderated SEs. This is an empirical-Bayes approximation
to sampling uncertainty, not independent raw-variance estimation. BH covers all
24,115 genes, assigning p=1 to non-estimable fits for correction while keeping
their state untestable. Five-cohort prediction intervals with positive estimated
heterogeneity are exploratory conditional intervals, not validated future-cohort
coverage. The fixed 1,854-gene signature is unchanged.

**First candidate result.** Job `22811901` completes all 24,115 models after
correcting an annotation join that included 4,232 Catalog rows without native
IDs. Independent REML and full-family BH checks pass. There are zero BH calls;
all fixed 1,854 pooled calls have unresolved random-effects support, although
1,850 retain the pooled sign and 1,208 agree in sign in every available cohort.
Minimum full-family q is 0.17143; median I-squared is 60.35%. This is a sensitivity with four/five studies,
moderated input SEs and conservative small-study intervals, not a test proving
absence of tissue effects. Review transportability claims with the full cohort
effects before any incorporation; do not replace the signature or search for a
more favorable meta-analysis. No eligible unexposed disease/control cohort was
established, so the blind confirmation arm remains unstarted. Full specification,
checks and results: `Analysis/Hypotheses/CATALOG-ADOPTION-20261107T022216Z/transportability/`.

### 2. Separate tissue composition from expression within a lineage

**Question and estimand.** For the fixed GNMT, MAT1A and CYP2C19 tissue example,
how much evidence supports a within-hepatocyte association? The tissue-stage
replication slope of 1.462 and its composition-adjusted value of 0.256 motivate
this question but do not estimate causal mediation.

**Analysis.** Synthesize and refit the already observed RV1 within-hepatocyte comparisons
using corrected donor-lineage pseudobulks with the 50-cell rule,
source-native diagnoses and independently recorded numeric histology. Keep
binary disease and numeric-stage models separate. Estimate the three gene
effects within each eligible cohort before combining cohorts. The RNA-count
model is edgeR normalization plus limma-voom; include source/batch and available
prespecified covariates only where identifiable. State donor n and reference n
for every fit. Report captured-cell fractions separately from tissue fractions.
Use the existing bulk and protein effects as tissue measurements, not independent
estimates of hepatocyte abundance per cell.

**Independent measurement and controls.** Before a new tissue decomposition,
identify public donor-linked histologic cell-area/count measurements or an
orthogonal abundance assay with a verified specimen join. Without that join,
the composition-adjusted result stays a sensitivity. Use fixed hepatocyte
identity markers as controls, and exclude the three targets from any expression
signature used for adjustment. Report positive controls even if they fail.

**Decisions.** Holm over the three genes within each declared primary contrast;
donor bootstrap for paired differences between models. A within-lineage claim
needs an estimable effect with corrected support and agreement in an independent
eligible cohort. A planning equivalence margin of 0.2 within-lineage SD may be
used for power calculations, explicitly not as a clinical threshold; freeze its
justification before opening validation outcomes. If adequate controls or power
are unavailable, retain the tissue-level interpretation. Deliver one tissue/
lineage effect panel and a table distinguishing changed composition, within-lineage
association and unresolved evidence. **Dependency:** corrected single-cell fits.

### 3. Establish a replicated genetic core with matched signals

**Question and estimand.** Which gene-trait associations retain support after
appropriate LD checks and recur at the same molecular signal? Keep direct
disease/liver-fat traits distinct from liver enzymes. Do not describe recurrence
through an enzyme GWAS as direct-disease replication.

**Analysis.** Start from COLOC r2, its 35 primary strata, native eQTL SuSiE fits,
the 14 genes selected for LD diagnostics and the existing independent-GWAS comparisons. Resolve
alleles, genome build, variant coverage, ancestry and LD-panel provenance before
refitting. Compare available ancestry/source-matched LD with the present reference.
Run the same LD-consistency diagnostics in predefined multi-signal positive
controls; deleting diagnostic outliers is a sensitivity, not an automatic truth
rule. Retain explicit failures and absent credible sets. Display PP.H4 at the
existing p12 values 1e-6, 5e-6 and 1e-5 without optimizing the prior.

**Replication and decisions.** Match the eQTL signal/credible set and the GWAS
signal across studies, using documented variant correspondence and LD, before
calling signal-level replication. A different lead variant alone neither proves
nor disproves recurrence. Record an unavailable ancestry-matched panel as a
limitation. The output is a locus-by-study evidence table with posterior,
eligibility, LD sensitivity and replication identity; posterior thresholds are
not FDR estimates. Any new frequentist screen needs its family fixed before
analysis. No automatic removal of the eleven sensitivity losses and no change
to named r8/r9 memberships follows from this plan. **Priority:** resolve the
existing diagnostic uncertainty before expanding the GWAS portfolio.

**Started, 2026-11-06.** Saved-signal inventory job `22811613` completed: 383
unique unordered gene-study comparisons, 94 with recurrent gene-level support
and 88 sharing a recorded eQTL signal index. The 94 represent 78 genes and are
all enzyme comparisons; no direct NAFLD recurrence is demonstrated. Shared
eQTL indices and shared-variant posterior compatibility do not establish matched
GWAS causal signals because GWAS-only credible sets and adequate LD remain
unresolved. Pilot job `22811698` refits MPRIP/MVP ALT with RORA/MVP NAFLD and
FABP1/UKBB ALT controls using the original reference. No ancestry-matched
alternative panel is currently documented, and no variants are deleted to force
agreement. Input alleles/build/order and original fit reproduction must pass
before interpreting the three-fit pilot or expanding it.

**First candidate diagnostic.** The three original-reference fits reproduce
exactly, but MPRIP has opposite-allele indel records collapsed to one position
ID with incompatible LD. Resolve source allele identity before deduplication;
arbitrary row selection is not justified. The immediate code correction refuses
SuSiE inference for ambiguous aligned variant IDs and retains a distinct
untestable reason. Guarded pilot `22812217` passes: MPRIP is refused as ambiguous,
while RORA and FABP1 reproduce exactly. Source variants are retained. Resolve
the underlying allele identities and LD before expanding the pilot. This does
not resolve source-matched LD or automatically change the
436 retained positive memberships. The saved 88 matching eQTL indices remain
insufficient to claim matching GWAS mechanisms.

### 4. Test whether the fixed protein programs add specific information

**Question and estimand.** Do the 22 previously selected protein programs retain
their histology associations in a second biological cohort, and do the fixed
programs outperform appropriate gene-set references? The existing 46-participant
protein cohort is discovery for this question and cannot supply confirmation.

**Inputs and analysis.** First inventory public liver proteomes with participant
IDs, histology, raw or documented normalized abundance and batch information.
Exclude overlapping participants and cohorts used to select the programs.
Freeze the 22-program directional family and the source-native histology endpoint
before reading outcome-linked protein values. Use the assay's documented log
abundance model, participant-level inference and measured batch/covariates;
report missing proteins and unchanged weight coverage. Retain the current
minimum of eight measured members and 20% of the original L1 weight; count
members once per primary-assembly gene. BH is over all 22 prespecified
replication hypotheses, assigning p=1 to untestable members for family
correction while retaining their untestable state.

**Specificity comparison.** In a separately declared family, compare all fixed
testable programs with gene sets matched on size, protein coverage, abundance,
lineage specificity and coexpression. Estimate matching features in an external
reference or training partition, never from held-out histology. Match positive
and negative weights and preserve correlation structure; a zero-association
null alone does not establish program specificity. Resample participants for
biological uncertainty and gene sets for vocabulary specificity, reporting the
two uncertainties separately. Require agreement of effect direction and
corrected evidence in the independent cohort; a null result leaves the current
single-cohort claim unchanged. **Stop:** no eligible independent proteome, failed
positive controls, or insufficient fixed-weight coverage. Deliver a replication
forest and specificity table. New specimen collection is outside this Resource.

### 5. Measure the usefulness of the Catalog for defined decisions

**Question.** Does the evidence record improve a stated decision beyond using
RNA, genetics or simple annotations alone? Keep two decisions separate:
predicting an independently measured disease-state direction, and identifying
which inherited mechanism is supported by a held-out molecular assay.

**Design.** Correct method/testability fields first. Freeze the eligible genes,
decision rules and source exclusions before withholding a complete cohort or
assay. Report any prior exposure; already examined outcomes permit only a
retrospective benchmark. Compare the full record with RNA-only, genetics-only
and annotation-only rules on identical eligible genes. Preserve missingness as
abstention. Include coding-only mechanisms so an expression-QTL absence cannot
count as an incorrect negative. No universal gene rank or modality-vote score
is introduced.

**Outcomes and decisions.** Report coverage, directional accuracy among explicit
predictions, the rate of unsupported definitive calls, and exact reasons for
abstention. Use donor resampling for measured disease effects and locus-cluster
resampling for inherited evidence; do not treat multiple genes at one locus as
independent validation. Holm-correct the six decision-by-baseline comparisons (two decisions and
three baselines), with the family fixed before evaluation.
An improvement requires a paired interval favoring the full record at matched
coverage; report the tradeoff if coverage differs. If the rules mostly request
missing evidence, describe that honestly as an information resource. Deliver
worked examples selected by the fixed criteria and a complete evaluation table.
**Dependency:** corrected Catalog fields and a documented independent endpoint.

### Execution order and limits

Finish and check the correction runs first, then synthesize completed
transportability results and address LD-sensitive loci. Pursue the lineage and
protein questions only after data eligibility and power are established.
Evaluate Catalog decisions last, with those inputs fixed. Broad cell annotation
and generic sequence-model comparisons remain supporting analyses; no new model
family is a dependency for these biological questions.

For each new run, record the exact hypothesis family, participant or locus unit,
estimand, covariates, missingness rule, exposure history, positive/negative
controls, seed (20261106 unless an existing fixed analysis requires its original
seed), environment and success/stop rule before reading new outcomes. Existing
exploratory results cannot be renamed confirmatory. Size power calculations to
the biological unit and verify that proposed effects are detectable before a
large run. Use CPU/nslab allocations within 16 CPUs and 210 GB; determine memory
from a bounded pilot, show the actual header before submission, and do not submit
arrays or new research workloads merely because they appear in this plan.

## Critical path

```text
Updated five-cohort bulk substrate ──┬─> stage + composition rebuild ──> Figure 3
                                    └─> bulk-dependent tables/Gene Catalog fields
                                    └─> passed lncRNA source check ──> lncRNA disease/state tables

Corrected COLOC rerun ────────────────> genetics audit ──> evidence classes
                                                       └─> noncoding DNA / target-biotype split
                                                       └─> Figures 1, 2, and 5

Accepted assay-native source tables ─> validated candidate Figure 5 ─> adjudication

Fixed external continuum loadings + selected bulk substrate
                                    └─> continuum molecular layers
                                        └─> working Figure 4F–4G + complete S3/S4 comparisons

All settled release products ────────> synchronized release audit ──> submission draft
```

The corrected COLOC rerun adopted on 2026-08-17 is superseded by COLOC r2, adopted
on 2026-10-14; its signal eligibility and trait drivers were corrected (see the
[2026-09-23 scientific review block](STATUS.md#what-is-official-updated-or-open)).
The remaining work is
release reconciliation: rebuild downstream consumers, close the outstanding
scheduler-disposition audit, and adopt one synchronized set of output files.

The 2026-09-22 scientific review identified source-diagnosis, COLOC
signal-overlap, Figure 2 inference, and Figure 5C corrections that must be
resolved before relying on the affected claims. Their evidence, bounded
corrections, and completion checks are in the
[technical correction plan](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23.md).
Results are recorded in the
[outcomes document](technical/SCIENTIFIC_REVIEW_CORRECTIONS_2026-09-23_OUTCOMES.md):
the Figure 2 inference, coding-overlap wording, Figure 5C code and Figure 7
wording corrections are complete as candidates; the 774-unit bulk refit was
adopted on 2026-10-11 and COLOC r2 on 2026-10-14.

## Open dimensions not yet on the critical path

**Two useful molecular models and compact Figure 7, candidate.** The bounded
[design](../scripts/analysis/masld_two_models/design.json) uses the existing
99-participant RNA-to-H3K27ac model, one additional sequence-by-RNA residual
experiment, and a four-development-fold comparison of native 1-Mb AlphaGenome,
the five-seed 2-kb adapter and their regularized combination. Currin fold 0
remains spent. The primary chromatin acquisition adds 20 measurements to 45
already available per fold (20/99 total participants; 54.1–60.6% of the selectable
pool); adding 10 or 30 is a sensitivity. The six-panel vector Figure 7 candidate
now has matched source tables, a passing independent numerical reconstruction
and completed rendered-panel inspection at
`figures/candidates/figure7_two_models_20260922_21849207/`. The public paired-H3K27ac and Wenz
sources require explicit eligibility resolution; absence of either narrows the
claims rather than prompting substitution with previously inspected data.
Both complexity decisions are now fixed by the completed development results:
retain global-plus-local RNA, since the sequence interaction misses +0.01
skill, and calibrated native 1 Mb, since the combination's 3.15% relative MSE
gain misses 5%. Fresh sequence inference remains experimental; no new-source
confirmation or empirical RNA read-depth accuracy has been established.
The integrated details are in [the research program](technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md#23-two-useful-models-and-a-compact-figure-7-candidate)
and current observations in [RESULTS](RESULTS.md#two-model-molecular-value-development-2026-09-22-not-adopted).
This candidate is not on the synchronized Resource release path; final adoption
remains a separate user decision.

**AlphaGenome–MASLD research program (requested 2026-09-14; integrated 2026-09-15).** The
[research plan](technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md)
asks which inherited regulatory effects are predictable from sequence, which
depend on context, and how they relate to human MASLD remodeling. Its evidence
table assigns prior results explicit next actions, including the completed
weak signed-splice result and unfinished full saturation.
T01–T24 retain their input, fine-tuning, objective and parameter menus across
six campaigns: accessibility, cellular context, assay transfer, remodeling,
indels/targets/haplotypes, and generalization. Each task specifies its evidence,
unresolved question, data, comparison, evaluation unit, measurable completion
criterion, dependencies, outputs and paper consequence.
Public-data expansion supports a primary database deliverable and reproducible
public model components; restricted research models remain separately
identified. Rigorous MASLD adaptations of methods from other diseases are
eligible contributions, with source-method comparisons and independent
biological evidence. The literature expansion adds 30 annotated references and
16 bounded comparisons, prioritizing supervision, donor differences, lineage
training, target links and molecular events. It records existing MASLD uses
before proposing novelty claims. One coordinator and three workers now implement
the authorized [initial week](technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md#194-authorized-initial-week-2026-09-15):
full existing-lead native scoring, scalar-effect adaptation feasibility,
donor/interaction development, P4 completion accounting and a local query specimen.
On 2026-09-16 the user removed campaign compute-hour and concurrent-GPU caps
and authorized CPU, bigmem, I/O and GPU partitions as needed. Scheduler limits,
the 500-GB new-file allowance and individual nslab submissions remain. Existing
P4 retrieval is accounted separately.
New experiments are not dependencies. Scientific scope and manuscript timing
are separate decisions. Final adoption and publication remain separate decisions;
existing findings and the separation from the experimental follow-up study remain in force.

The initial development checkpoint supports six adapter gains over matching
fixed controls. The stronger donor comparison and uncertainty-aware fits are
complete and independently checked. The repaired matched native table is also
complete, with independent numerical reconstruction pending: its point estimates
favor native ATAC over the best adapter. RNA benefits held donors at trained
regions; crossed held-donor/held-region benefit remains unestablished. Review
the complete native comparison before choosing a distinct residual or
same-length comparison and any nested five-seed finalist. The measured regional
ATAC assessment and treatment-interaction negatives argue against extending the
same partial-update or reporter design without a changed hypothesis or
information source. Finish P4 and prioritize the missing target/event/phase
measurements for the less-developed branches. New lineage heads, independently
supported mechanisms and protected generalization remain unfinished; the
[T01–T24 progress table](technical/ALPHAGENOME_MASLD_RESEARCH_PROGRAM.md#t01t24-progress-and-next-scientific-decisions-2026-09-16)
distinguishes those goals from completed pilots. Reviewed values are in
[RESULTS.md](RESULTS.md#initial-regulatory-context-comparisons-2026-09-15-not-adopted).

**AlphaGenome–MASLD program, executed 2026-09-14/15 (candidate research layer).**
Work packages A, B, C, D, E and F of that program have now run and are reported
in [ALPHAGENOME_PROGRAM_RESULTS_2026-09-14.md](technical/ALPHAGENOME_PROGRAM_RESULTS_2026-09-14.md),
with a Gene Catalog attachment under
`GWAS/finemapping/results/alphagenome_program/g-integration-20260915T072500Z/`
in which every row is marked `candidate_not_adopted`. **Nothing in it is
adopted and it changes no adopted number, gene class, evidence state,
threshold, hero gene or figure.** It is not on the critical path and the
synchronized release does not wait on it. Local AlphaGenome weights are now
available, and local scoring plus the reporter probe with fixed weights have run; their
candidate results and subsequent corrections are linked in sections 3.1–3.3 of the integrated plan. The initial
work packages did not complete the broader context-dependent modeling aims.
The next model comparisons include the strong local zero-shot baseline,
matched inputs, adapters and broader tuning, while the data tasks develop
independent molecular-effect and state labels. Historical authority conflicts
identified by the program remain documented in its technical defect register;
this prospective plan changes no dataset permissions or weight-release policy.

The user-approved continuum translation extension now has an implemented and
checked candidate figure, complete Gene Catalog sidecar and Methods section;
see the [technical review](technical/continuum_translation_20260906T135538Z/IMPLEMENTATION.md).
The remaining manuscript decision is whether to adopt that concrete candidate
within the synchronized figure/table release. Conditional genetic directions,
cell-specific mechanisms, target-to-readout regulation and treatment-response
prediction remain unresolved experiments, not required positive results or
reasons to change thresholds. Coherence stays excluded until N1 reconstruction
and N2 calibration are resolved. This extension does not authorize new patient
recruitment or move the future experimental study into the Resource.

Both dimensions below are stated at goal level in `docs/PAPER.md` and neither carries a claim.
They are recorded here so they are scheduled rather than rediscovered.

**Aspect decomposition — metabolism, inflammation, fibrosis.** The three aspects
have never been resolved as separate outcome axes; NAS is analysed as a graded
sum and as an adjustment set. Evidence state `indeterminate`. The substrate is
already kept separate and summarised, and the pooled n across GSE267145 (99),
GSE202379 (40), PXD051911 (58) and GSE135251 (180) is **377**, which is a
different regime from any single cohort attempted so far. The re-test needs a
recorded prespecification before any component is examined, because a question
that has been looked at once informally cannot be re-opened casually.

**Variant combinations and chromatin effect.** Scoped and unstarted. One local
haplotype at GNMT is characterised and unresolved; the colocalization design is
single-variant throughout. Chromatin is evaluated as assay-native context with
the large majority of jointly testable peaks `indeterminate`. Any move from
context to effect requires its own design, not a reinterpretation of the
existing peaks.

## Work packages

| Order | Work package | Can start? | Completion test | Main outputs |
|---:|---|---|---|---|
| 1 | Settle documentation and scope | **Complete** | Five human entry points; archived originals; scope validator passes | `docs/README.md` through `ROADMAP.md` |
| 2 | Finish corrected COLOC | **r2 adopted 2026-10-14 after independent source checks** | All 50 studies present, 948,503 unique gene-study/chromosome rows, all five recorded output hashes matched; 1,307 shared-posterior corrections complete across 32 studies, other results carried unchanged; primary 436, ABF-only 560, union 996, full-portfolio 708 | Selected `susie_coloc_r2_20260924T053244Z` and synchronized genetics panels |
| 3 | Rebuild cross-sectional stage models | **Complete on 774 units; panels installed 2026-10-14** | Corrected genetic-identity-cluster fits use complete 24,115-gene BH families; stage, F0, endpoint, five per-cohort and leave-cohort-out sources checked; endpoint roster 392, five-cohort union 12,180 (GSE126848 filtered refit, user decision D3, 2026-10-21) | Current Figure 3 and S3 panels, sources and vector review assembly |
| 3a | Add genetics–state interface panel | **Priority eligibility adopted 2026-11-06** | 3,589 eligible native RNA/SuSiE pairs, 426 genetic and 209 RNA calls, 26 shared; the historical 426-pair interval-null sensitivity remains separate | Figures 1B/3G and S3G, sources and qualifications |
| 3b | Add histology-anchored molecular continuum | **774-unit panels installed; Figure 4 review assembly complete** | Fixed 117 programs and 43 systems preserved; stage-complete inference 435 units, including 68 GSE213621 controls with no recorded histology coded as stage 0, two co-primary scores; 33 primary BH/direction calls and 18 strict maxT calls distinct from competitive specificity (0 connected, 1 uniform), all on those 435 units. After review round 3 (R05-02), main 4G and 4H use the 367 disease units; 4H has no outline, rho 0.63 and 32/43 same direction. Complete S3/S4 sources retain failed NMF stage axes and untestable signatures. | Current Figure 4F–H and complete S3/S4 comparisons |
| 4 | Rebuild cell composition | **Complete on corrected 774-unit substrate** | Composition panels and matched 367-unit disease-only continuum adjustments checked; Stromal ECM depends on deconvolution input, no adjusted system nominated | Current S3 composition and qualified continuum effects |
| 5 | Reconcile bulk-dependent Gene Catalog fields | **Eligibility correction adopted 2026-11-06** | 28,233 identities retain the adopted RNA effects/calls and 436 native genetic positives; exactly 188 genetic states and 180 classes change under the existing priority rule; 855 held annotations remain untestable | `gene-catalog-eligibility-20261106` and matching molecular-object companion |
| 5a | Resolve the lncRNA strand rule | **Complete** | Five cohorts pass the documented amended rule; no cohort or sample was removed; stricter failed sensitivity retained | Hashed authorization for disease-associated lncRNA analysis |
| 5b | Add noncoding layers | **RNA package and five-panel DNA-only package validated; gene-level target companion installed; transcript isoforms unresolved** | Stable IDs, all seven genomic classes, all-gene BH families, credible-set PIP architecture, target biotype, assay applicability, and object-specific Catalog routes share one release | Cross-cutting additions to Figures 1–5 |
| 6 | Regenerate assay-observability figures | **Complete for Figure 5; supplementary/Catalog integration remains candidate-only** | Spatial candidates plus ATAC v3 exact counts, 234 program rows, official chromVAR, explicit spatial states, the adopted COLOC replay, and variant-consistent Figure 5D passed validation; Figure 5A–F and its vector assembly are canonical working outputs | Physical-context, observability, and robustness panels |
| 6a | Benchmark continual single-cell integration | **Complete and CLOSED; primary adoption rejected 2026-08-15; Harmony stays primary** | The V42 decision met 13/14 predefined criteria, including the read-only census, the control-only selection settled in advance, internal five-seed/held-study/label/disease/six-order matrix, donor-bootstrap alignment improvement, full 1,232,318-cell projection, six-model full-scale expansion, bit-identical GPU reproducibility, and the byte-identical 117-program separation. It FAILED `independent_external_routing`: untouched GSE212837 donor-balanced macro-F1 0.2756 against the predefined 0.70 floor. The prospectively replaced pretrained-CellTypist Router V2 was then rejected by V46 because untouched GSE296875 cleared the overall 0.70 floor at 0.7071 but missed the predefined 0.50 per-class floor for B cells (0.2888) and NK-T (0.4229). The strict seven-donor reference is retained; the local nine-donor arm and the common-universe 26-donor HLiCA arm both failed reference preservation. Note that the queued B6K sbatch fits (`20208726`, `20215473`) never ran and every fit was executed in interactive GPU allocations instead | No candidate Figure 4 representation. The supplementary benchmark may be reported only as a negative sensitivity result; a further router is a new method study, settled in advance, on a third untouched cohort |
| 6b | Compare snATAC disease estimates and regulatory links | **Complete; qualified manuscript inclusion approved 2026-10-26** | All six count reconstructions conserve donor counts and four corrected primary fits reproduce. Robust QL does not improve fitted-NB calibration. The 30 outer folds and 48 original-grid contrasts retain 0/24 endpoint-baseline successes. Wider-grid hepatocyte gains remain separate; Cicero retains a 0.29–0.34% advantage over distance across every unshuffled stress test, and similar gains survive pairing shuffles. All 18 RNA libraries are reconstructed; only hepatocyte disease RNA meets donor eligibility. Independent MASLD validation remains unavailable; native-method expansion is outside this completion | Figure 5 text, Supplementary Figure S5P–AA, Note 16 and complete source tables; primary voom, fixed programs and Catalog classes unchanged; synchronized release separate |
| 7 | Rebuild genetic/state classes | **Eligibility correction adopted 2026-11-06** | Native interface 14,704, qualified overlap 3,589 and jointly eligible Catalog 14,502 distinguished; classes 400/26/1,118/12,958/13,731; ATAC r2 context remains 767 pairs | Figures 1B/3G/6 and selected Gene Catalog |
| 8 | Rebuild the MASLD Gene Catalog and portal export | **Local Catalog correction adopted 2026-11-06; external portal unpublished** | Full index, evidence states, follow-up routes, molecular records and local downloads use the same approved 180 class changes, with separate SuSiE/ABF values and corrected RNA applicability | Figure 6 and local Catalog; [current numbers](RESULTS.md#catalog-eligibility-adoption-2026-11-06) |
| 9 | Assemble synchronized manuscript release | **After 3–8** | Numbers, claims, figures, methods, data dictionary, and hashes agree | Submission candidate |
| 10 | Scientific adjudication and adoption | **User decision** | Scope, biology, integrity, and reproducibility reviews pass | New current official release |

## Figure completion checklist

The October 24 matched-gene wider-ridge follow-up is complete: Cicero adds
0.306% hepatocyte RNA-error reduction beyond distance (Holm p = 0.000400 over
four tests); adapted DIRECT-NET worsens there. Both beat distance in Mesenchymal
cells without establishing superiority to promoter-only prediction. This is a
separate post hoc candidate, with original results unchanged.

External validation remains contingent on usable donor-level peak/RNA inputs.

The user-selected October 25 completion focuses on existing MASLD data:
the 42 fixed-fold wider-grid robustness refits and RNA qualification are complete
and combined with matched disease sensitivities and descriptive all-donor link
maps. Three unchanged omission folds reuse exact baseline predictions. External cohort searches, author contact and
native-method expansion are outside this completion. Its conclusion retains
untestability when cross-assay source keys are unresolved. Approved manuscript
inclusion does not change Catalog classes or qualify individual regulatory links.
All 18 RNA libraries pass exact per-cell/per-gene molecule accounting. Only
hepatocyte RNA meets the fixed donor threshold (5 normal/9 MASH); stellate and
macrophage RNA remain untestable. The RNA arm uses deposited diagnoses; two
library diagnoses differ from the local cross-assay crosswalk, which does not
authorize changing ATAC labels or claiming verified pairing. Donor-matched
associations remain untestable. Independent predictive reconstruction and all
twelve individual PDF reviews pass. On October 26 the user approved qualified
working-manuscript inclusion, including supplementary links and RNA candidates;
the synchronized release still follows the final checklist below.

The October 24 alternative-cohort audit identifies specific missing inputs:
SenNet peak counts/fragments, GSE223843 common-peak reconstruction, and GSE281574
barcode-to-donor assignments. No external predictive result is available.

Historical Figure 4B–H artwork was updated after user review on 2026-10-05.
The October 14 panels now use the selected 774-unit bulk effects and qualified
source diagnoses. D retains the complete discovery/content lineage comparison;
E uses the three-cohort 364-unit arm of the 774-unit substrate after review
round 3; the 457-record arm is a sensitivity. The synchronized publication
release remains pending the checklist below.

| Figure | Must be settled against the release before final assembly | May be drafted now |
|---|---|---|
| **1** | Final cohort/assay coverage, two-map interface, joint-testability, seven-state vocabulary, and evidence-state counts | Spatial capability inset and noncoding-aware schema are validated; lead with inherited variation, multicellular remodeling, and cell/tissue context, and omit adjusted coding-versus-lncRNA matching until complete versioned eligibility tables exist |
| **2** | Corrected COLOC, trait/ancestry provenance, method agreement, synchronized credible-set DNA architecture, target biotype and exact examples | COLOC r2 is adopted on 2026-10-14. Figure 2A/B/C/E use the checked r2 PDFs; F/G/H/I are render-identical with synchronized sources. Figure 2D and its unassembled context companion retain the adopted October 13 coordinate correction. The A–I sequence is unchanged. The separate target-gene companion (424 protein-coding, 12 lncRNA) and exact-example table are installed; separate supplementary selection remains open |
| **3** | Updated bulk counts, five-cohort cohort/leave-one-out validation, prespecified endpoint sensitivity, five-cohort adjacent-stage and NAS results, corrected composition, continuous NMF axes, the GNMT–MAT1A–CYP2C19 RNA–stage–protein state, within-MASLD protein covariation, and exact signal-pair genetics lanes | Figure 3C/D endpoint-validation candidate v3 is complete; the protein-covariation candidate still needs an independent biological cohort, despite computational reconstruction of its 46 participants, and 3F artwork is complete (2026-10-05), with historical stage inference conditional and the PIP column removed |
| **4** | Seven-dataset single-cell atlas, the fixed 117-program registry, complete signed disease forest, ambient-effect transport audit, donor-paired descriptive localization for the Stromal ECM and Ductular injury programs, bulk tissue-state transport, and the histology-anchored continuum; detailed heatmaps, communication, and TF activity remain supplementary | User-approved working main 4D–4G and complete S3/S4 continuum comparisons are adopted without changing Figure 3. The ambient audit remains required source-table provenance but is not displayed; synchronized Resource-release adoption remains open |
| **5** | Validated assay-native context tables/PDFs and synchronized adoption approval | All six panels and the vector assembly are working canonical outputs. 5D uses the finished primary-pair context from the adopted COLOC release across five shared lineages; 5E/5F use program-level labels; 5F pairs matched-gene calibration with one illustrative outcome-blind Guilliams et al. H&E overlay. Full Resource-release synchronization remains open |
| **6** | Selected r8 gene classes, molecular-object fields, evidence states, unresolved alternatives, deterministic next-experiment rules, and worked examples | r8 has 28,233 identities and completed local review; the molecular-object companion and interface preserve all gene calls, all 117 programs and explicit 774-unit/855-held restrictions. Historical Myojin class tests remain historical; no causal or signal-target claim is inferred |
| **7** | Task-specific scores, source exposure and transfer limits from the current Figure 7 section of `RESULTS`; model-development alternatives remain separate | Current 7A–7F are complete, with corrected audit-entry units and 7C/F exposure labels. TranscriptFormer has five clean and two seen studies; Geneformer has zero clean studies. Broad cell-class mapping and the five-clean-study donor-bootstrap comparison against inductive scVI are separate estimands. The embedded six-panel layout is a review proposal, with final journal sizing open. Earlier 7B and molecular-model alternatives retain their own candidate state. No follow-up screen quantity may appear |

## Final release checklist

The former “Plan 60” is now this human-readable checklist. A release is accepted
only when all boxes can be answered from read-only files:

- [x] One release ID names every paper output file: release 2026-10-21-r9, with manifest
      `RNA-seq/results/manuscript_release/2026-10-21-r9/release_manifest.tsv`.
      r4 supersedes 2026-10-20-r3 and 2026-10-17-r2 (manuscript rewrite, supplementary redraws); r2 superseded r1.
- [x] Supplementary Figure S2-S4 panels redrawn on the current release.
- [x] Source files, cited producer files and the release manifest are read-only.
      Environments are recorded by `pip freeze` and conda export in the release
      directory rather than locked.
- [x] The chosen five-cohort fragment-count bulk substrate is used consistently;
      the July Figure 5 display universe remains fixed and explicitly labeled.
- [x] Stage and composition analyses use selected genetic-identity clusters,
      one library per cluster, with donor verification limits stated.
- [x] Cross-cohort single-cell states reproduce source diagnoses under the
      Healthy/Steatosis/Steatohepatitis/Cirrhosis vocabulary; generic MASLD is
      not treated as a stage. Original-source checks correct GSE174748 to
      Cirrhosis, retain eleven Liver Atlas categories as Unresolved and preserve
      all nineteen Liver Atlas exclusions; exact joins and citations are in the
      [source correction](../Analysis/SingleCell/candidates/source-diagnosis-primary-correction-20261014/source/README.md).
- [x] All 50 r2 study outputs and 948,503 gene–study rows are checked;
      the 1,307 reevaluated rows have terminal dispositions: 1,233 retained,
      49 eligible below threshold and 25 untestable. The older 1,100-task
      execution accounting is historical, not the r2 verification criterion.
- [x] Genetic/state classes are rederived from the selected r2 COLOC and adopted
      774-unit bulk sets; raw, native-ID and adjudicated universes are distinct.
- [x] The adopted COLOC release rechecks exact gene–trait posteriors and
      shared-variant anchors for GNMT, MAT1A, CYP2C19, IGFBP7, and BICC1 without
      requiring genetic support for their disease-remodeling roles.
- [x] If GNMT remains the worked locus, its risk/eQTL direction, local LD
      haplotype, prior human liver protein evidence, and exact corrected-COLOC
      signal pair agree in one read-only source table.
- [x] All 117 prespecified programs, weights, stability fields, and testability rules are unchanged.
- [x] The complete 117-program ambient-RNA sensitivity uses the current
      source-diagnosis metadata, is independently validated, resolves or
      explicitly retains the GSE189600 decontX failure, and corrects atlas QC
      metadata that claims decontX was applied to the original substrate.
- [x] The five-cohort strand rule is explicitly resolved before any lncRNA disease claim.
- [x] Noncoding DNA context and colocalized target gene biotype are displayed separately; transcript isoforms remain unresolved.
- [x] Same-strand-overlap and complex lncRNAs cannot become named examples:
      all 335 such r8 records are ineligible for main-text/high-confidence
      examples; current named-example inspection passed.
- [x] Adjusted cross-assay lncRNA observability is omitted because complete,
      versioned eligibility denominators are not available.
- [x] Main-figure source tables reproduce every plotted value exactly (46 panels,
      panel-deep review 2026-10-17); supplementary panels were not rechecked.
- [x] Figure redesign v6 passes panel cardinality, one-page PDF, Ghostscript,
      Type 3 font, checksum, composite-proof, and Resource-scope validation;
      the corrected r2 panels and all five individual Figure 6 panels are now
      installed; complete publication assembly remains.
- [x] Figure 5B/C now use 774-unit RNA and 5D uses verified r2; the earlier
      v6 and September panel disagreements are retained as historical.
- [x] ATAC v3 non-genetic outputs have exact fragment recounts, complete prespecified-program families, official chromVAR provenance, and candidate-only integration completion markers.
- [x] ATAC v3 genetics replay, posterior export, batch validation, liftover normalization, path guards, and upstream fail-closed behavior passed compute-node fixtures without producing biological rows.
- [x] ATAC v3 variant-posterior replay and variant-consistent Figure 5D passed after the corrected COLOC release was adopted.
- [x] Release-linked Figure 1/Figure 6 spatial insets and supplement panels are validated; canonical assembly remains pending.
- [x] Full-composition, unit/null/weight, cell-context, and HMSMA mask diagnostics are validated; only the visually accepted read-only candidates are eligible for assembly.
- [x] Spatial publication panels use one decision per panel; the impact family foregrounds matched-null calibration, moves the complete source-state matrix to the supplement, and preserves exact source-table linkage.
- [x] Gene Catalog entries distinguish support, discordance, negative tests, untestability,
      non-applicability, source dependence, and indeterminacy.
- [x] Every indexed working Figure 1–6 source table carries or joins to the applicable evidence
      state, state reason, unresolved alternative, and next-experiment rule.
      The [checked source-context companion](../figures/candidates/figure-claim-context-20261014T185701Z_22240302/checks.json)
      preserves native values and explicitly retains unmatched identities.
- [x] CosMx 12/117 observability, GeoMx 0/117 testability, GSE287826 0/117
      testability, ATAC state counts, and Catalog `tested_negative` calls (not assessed) are presented as evidence-state guardrails around the inherited-
      variation and multicellular-remodeling maps, not as the paper headline.
- [x] Current manuscript numbers and figure legends agree with
      [RESULTS.md](RESULTS.md).
- [x] The local portal and Figure 6 use the same unchanged r8 scientific calls;
      the reviewed copy and complete molecular-object companion are recorded in
      [RESULTS](RESULTS.md#catalog-review-and-molecular-object-companion-2026-10-14).
      This check publishes no external portal.
- [x] No follow-up screen result, screen roster, phenotype, hit, or guide performance entered
      the Resource selection or narrative.
- [x] The scope validator and the release manifest check
      (`scripts/manuscript/check_release_manifest.py`: checksums plus main-panel
      source tables) pass; the July-only `check_release_consistency.py` is retired
      for this release.
- [x] The user explicitly approves promotion: approved by the user 2026-10-17.

## Reopened public-data hypothesis testing (user decision, 2026-10-07)

Public-data mechanism hunting is reopened. New hypotheses may add claims or
panels to the Resource. Limits: open-access data only; the 117 Hotspot programs
and weights stay fixed; no follow-up-screen result enters the Resource; donors
or participants are the unit of inference; each test has a prespecification
written before outcomes are read; nothing reaches a figure without a second
dataset or a passed negative control; adoption needs explicit user approval.
None of these is a submission dependency except Step 0.

Novelty was checked by literature search on 2026-10-07 (PubMed, Europe PMC,
GEO); "not found" means not found by those searches.

| Order | Item | Question |
|---|---|---|
| 0 | Bulk refit | Does five-cohort remodeling hold with one library per person and alternate-contig reads credited; 855 unmatched alternate-contig annotations are reported as not testable |
| Check | V1 | Share of bulk calls explained by the perioperative glucocorticoid response (GSE320593) |
| Check | V2 | Hepatitis B, C and E transcripts in the public cohorts |
| Check | V3 | ECM and Ductular programs in donors never used to define them, HC3 errors |
| Check | V4 | Spatial organization on diseased tissue (HRA007511, Vu 2025) |
| Check | V0, V5-V9 | Fixed programs in seven new public datasets (2026-10-21); plan `~/.claude/plans/i-agree-with-gse298719-whimsical-rivest.md`, D0 in `Analysis/Hypotheses/V-GATES-20261021T121655Z/D0_DECISION.md`. V7d (SCP3359 Xenium) panel gate failed 2026-10-22 |
| 1 | H26 | Inhibit or activate call per anchored gene from eQTL direction and rare predicted gene-disrupting variant burden |
| 1 | H10 | Cirrhosis effect per unit ALT at each enzyme-only locus |
| 1 | H21 | Programs with polygenic signal for cirrhosis but not liver fat (MAGMA competitive test) |
| Corrected 2026-10-14 | N4 | Retrospective within-person agreement on historical counts; 22 different-visit pairs remain identity-unqualified and their exact-stage null is untestable. Historical 69 stage-BH and 87 random-null calls withdrawn. The 28 untyped pairs yield zero eligible-null BH calls. No further N4 job is pending; primary ICC/ceiling claims remain excluded without qualified identities. See [RESULTS](RESULTS.md#n4-within-person-agreement-eligibility-correction-2026-10-14). |
| 1 | N1 | Effect of undisclosed duplicate individuals and mixed libraries on published cohort results |
| 1 | H40 | Alcohol-specific expression signature in MASLD-labeled biopsies |
| 1 | H34 | Loss of Y-chromosome expression in liver immune cells by stage |
| 1 | H48 | Mouse hepatocyte knockouts ranked against the human disease signature |
| 1 | H3 | HSD17B13, PNPLA3, TM6SF2, MTARC1, GCKR dosage and program score beyond stage |
| 1 | H24 | Target gene and cell type at direct-disease loci without bulk eQTL colocalization |
| 1 | H30 | FOXO1 hotspot mutation detectability in bulk biopsy RNA |
| 2 | H25, H22, H2, H8, H35, H36, H47, H49, H17, H9, H31, H46, H50, H41, H27, H7, H38, H23, H42, N3, N5, N2a, N6-statin | Narrowed forms only; run after the refit |

Round-3 plans H51 to H60 (2026-10-15; prespecified and run; H56 stopped at its positive-control rule) and the confirmed review findings they correct: [RESOURCE_REVIEW_AND_HYPOTHESES_2026-10-15](technical/RESOURCE_REVIEW_AND_HYPOTHESES_2026-10-15.md).

Round-4 corrections and Catalog r9d were installed on 2026-10-16 after user approval; see the STATUS entry of that date. Item O (GSE289173 registry exposure record) was applied the same day.

## Analyses that are closed for this paper

| Historical work | Resource disposition |
|---|---|
| Plans 40–44 public functional/mechanistic tests | Complete supplementary boundaries; valid nulls and limitations retained |
| Yakubovsky lipid/zonation work | Reference context only; not MASLD validation |
| GSE287826 | Could not be tested without a public donor key |
| HMSMA clinical/histology/MALDI inference | Cannot be used without authoritative donor, clinical, image, and metabolite registration; label-blind coverage remains valid |
| Existing 22-program spatial family | Regression/provenance only; not a biological result family in the new release |
| GeoMx spatial inference | Dropped because unresolved direction/normalization prevents trustworthy gene or program evidence |
| SteatoSITE EHR, NAS, outcomes and processed data | Excluded: permission-only (PMS-IC/EPCC application). The public ENA reads and NASH-CRN fibrosis stage are used instead as a held-out replication cohort (SS series, user decision 2026-10-25; `Analysis/Hypotheses/SS-20261025T135237Z/PRESPEC_SS.md`). S1-S4 complete 2026-10-28 (RESULT_SS.md; candidate Fig 3E replication strip in `figures_fig3e_strip_20261028T035939Z/`; Catalog sidecar candidate `RNA-seq/results/evidence_passports/candidates/gene-catalog-r8-steatosite-stage-replication-20261028T035950Z/`); adopted 2026-10-30 by user approval (strip, scatter, Methods, Results, Catalog field). Planned for later: a supplementary pooled sensitivity check (five cohorts plus SteatoSITE, dataset term, Fig 3E adjacent steps); main Fig 3 stays discovery-only (user decision 2026-10-27). SSX extension: exon-intron split and NMF axes adopted 2026-11-03 as Figure S3 panels; eQTL direction one Note 17 sentence |
| Plans 45–46C | Future experimental/collaborator protocols, outside the Resource |
| Historical ncRNA scripts 53–57 | Excluded from the Resource pipeline; ceRNA, live-query, cross-species, and outcome-selected analyses are not reused |
| New biomarker-performance headline | Excluded without prospective outcome validation |

Closed analyses may be cited as stress tests or limitations, but they cannot
select a new hero gene, program, contrast, or main-figure claim.

## Writing sequence

1. Settle title, abstract logic, Introduction, and Discussion around the Resource
   contribution in [PAPER.md](PAPER.md).
2. Finalize stable Methods and explicitly mark release-dependent subsections.
3. Draft Results by figure using qualitative statements for values not yet available.
4. Populate exact values only from the synchronized release ledger.
5. Run claim, number, figure-source, scope, and link checks.
6. Prepare journal-specific title/abstract/cover letter only after the scientific
   release is settled.

## Journal decision

Default to the strongest journal that evaluates the work as a reusable Resource
and analytical framework. Current working order:

1. Cell Genomics Resource-style submission.
2. Nature Communications or Genome Biology, depending on final biological
   coherence and portal/reuse strength.
3. Cell Metabolism or Nature Metabolism presubmission inquiry only if their
   current article format clearly accommodates this Resource contribution.

Do not add causal language or reopen a mechanistic dependency to chase a journal
category. Article types and author instructions are rechecked immediately before
submission.

## Handoff rule

Every agent should begin with [README.md](README.md), [PAPER.md](PAPER.md), and
[STATUS.md](STATUS.md). Update this roadmap only for a directional decision or a
completed dependency. Put detailed commands, job diagnostics, and dated results
in technical records or the archive—not in the human-facing roadmap.
