# Combined scientific review of the MASLD Resource

**Date:** 2026-10-07  
**Purpose:** reconcile the earlier bulk/single-cell review with the supplied
hypothesis review, evaluate their claims against current code and available HPC
results, and identify what the evidence supports. This document is a review;
it does not authorize new analyses, change the paper's scope, or adopt results.

## Sources and how to read this document

The inputs are:

1. [Bulk and single-cell RNA-seq review](BULK_SINGLE_CELL_RNA_SCIENTIFIC_REVIEW_2026-10-07.md),
   including its independent numerical checks from this conversation.
2. [Other agent's findings](02_WHAT_IS_BROKEN.md).
3. [Hypothesis-review archive](Hypothesis_Review_2026-10_rev1.tar.gz), especially
   its README, findings, top-20 execution plan, hypothesis inventory, and
   parallel-track proposal. Archive members were read without executing or
   unpacking their code.
4. Current [scientific scope](../PAPER.md), [status](../STATUS.md),
   [results](../RESULTS.md), [roadmap](../ROADMAP.md), and the actual code and
   saved outputs cited below.

The standalone `02_WHAT_IS_BROKEN.md` is byte-identical to the archive member.
It is one review supplied in two forms, not two independent confirmations.
Its README explicitly says that its audit had tracked code/configuration and
an August manuscript, but could not access HPC result tables or private docs.
That explains several differences from the current scientific state.

Here, **confirmed** means the stated mechanism or result is supported by the
inspected evidence. **Historical** means it applies to an older analysis or
manuscript. **Qualified** means the concern is valid but its stated scope or
conclusion needs narrowing. **Unverified** means the relevant numerical result
or current use was not independently established. These labels do not imply
that every other result in the affected subsystem has been reviewed.

## Combined assessment

There is substantial legitimate work here: appropriate RNA-seq models, donor
collapse in the current single-cell program analysis, complete testing
families, fixed program definitions, and explicit reporting of inconclusive
and untestable results. The descriptive atlas and conditional evidence for
broad disease-associated expression changes are defensible.

The strongest current limitations are participant identity, unfinished
counting correction, fragile native single-cell disease associations, and
interpretation of correlated measurements as independent confirmation.
Several older claims about a universal convergence ranking, discrete NMF
patient types, and a regulatory cascade are outside the current paper's
allowed claims. They must not be revived from old scripts or draft text.

The other review is useful as a historical defect inventory and a source of
specific code questions. Its severity ranking is not a current assessment of
the whole Resource: it mixes historical and current analyses, sometimes treats
code comments as measured results, and makes several stronger statistical
conclusions than the evidence supports. Its top-20 plan should not be executed
as a package under the current Resource scope.

## Findings that matter most now

### 1. Exact bulk gene and stage totals remain conditional

The historical disease/control model contains **844 library records**, not a
verified total of independent people. Current identity candidates retain
**774** supported source/genetic units, with **738** in the strict-depth
sensitivity. Genetic clusters are useful identity evidence but are not
independently verified people. GSE193066 remains unqualified for primary
replication, repeat-biopsy ICC, or derived reliability claims.

On historical unique-only counts, the identity sensitivities have **1,756** and
**1,807** canonical annotation-row calls. The 774-unit fit retains 1,572 of
1,725 historical calls and has logFC Pearson correlation **0.988**. This
supports persistence of a broad expression pattern while leaving exact
membership conditional. Changes also reflect sample composition, filtering,
normalization, and quality weights; they do not isolate the effect of repeated
libraries.

The r7 annotation correspondence still has **855 unresolved records**; **63**
reference artifacts have been explicitly excluded. Duplicate targets and
unresolved annotations occur among historical calls. The current code requires
completed annotation and count evidence before primary use. Complete resolution
is the user's selected analysis requirement, not a universal claim that any
unknown scaffold invalidates every primary-chromosome RNA result. A narrower
estimand would require its own explicit definition and scientific decision.

Evidence: [current counts and qualifications](../RESULTS.md),
[qualified bulk inputs](../../scripts/manuscript/resource_f_five/qualified_bulk_inputs.R),
[r7 summary](../../figures/candidates/resource-corrections-20261005-v1/annotation_resolution_r7/resolution_summary.json),
[identity-only comparison](../../figures/candidates/resource-corrections-20261005-v1/identity_only_comparison_r2/comparison_metrics.tsv).

The erroneous GSE130970 control labels were corrected in a completed
historical-roster fit: six source healthy controls were distinguished from F0
disease biopsies. The resulting **1,725** calls are a corrected-control result,
not a completed identity/counting correction. Label-dependent analyses must
use that correction; analyses that read neither the control label nor the
disease/control fit do not automatically change. In particular, the existing
43-system estimates are not recomputed merely by relabeling those controls.

### 2. Current native single-cell disease evidence is weaker than either old headline

The current atlas has **104 descriptive donors**, **102 analyzed donors**, and
**64 donors** in the focal native disease models. The other review's 73-donor
stage count and “two robust plus two fragile” result belong to earlier
metadata/analysis states.

| Native program | Ordinary BH q | HC3 BH q | Current interpretation |
|---|---:|---:|---|
| Stromal ECM | 0.0496 | 0.1053 | Ordinary selection; sensitive to the variance model |
| Ductular injury | 0.0187 | 0.1053 | Ordinary selection; sensitive to the variance model |

The complete current family has **two ordinary BH calls and zero HC3 BH
calls**. The prior review independently recalculated both families. Failure
of HC3 correction does not prove no effect; it prevents calling these robust
native disease discoveries.

The current native model is a common equal-step Healthy/Steatosis/
Steatohepatitis trend with dataset adjustment. It is a cross-sectional diagnosis
association, with limited demographic adjustment. It does not measure
within-person progression. The current hepatocyte design includes **26
steatohepatitis donors in GSE202379 and nine in GSE244832**, plus one in
GSE189600. Thus the historical statement that 99.94% of steatohepatitis cells
came from one dataset cannot be substituted for the current donor-level design.
Those cell counts refer to the old Progressor/subcluster analysis, which had
its own confounding problem.

Evidence: [current design](../../Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/stage_dataset_design_audit.tsv),
[current effects](../../Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/work/primary_effects.tsv),
[current producer](../../Analysis/SingleCell/scripts/hotspot_modules/512_hotspot_v2_donor_refit.R),
[historical Progressor producer](../../Analysis/SingleCell/scripts/361_cross_cohort_progressor.py).

### 3. Program origin and disease association are separate questions

**37/117** fixed programs have a highest content-associated lineage different
from their discovery label. In the matched ambient sensitivity, Stromal ECM
changes from β=0.401 to 0.130, with corrected BH q=0.312; Ductular injury changes
from β=0.396 to 0.389, with corrected q=0.129. Neither correction-induced change
passes its full BH family. These observations raise a lineage/source concern;
they do not establish that contamination explains every mismatch or effect.

The targeted correction-available fibroblast ECM comparison is positive in
**26 donors**: β=0.329, 95% CI 0.098–0.560, two-program BH q=0.0106. The
cholangiocyte Ductular injury contrast remains unresolved (q=0.449). These are
retrospective projected-score comparisons with limited cohorts and possible
substate mixture; they do not establish cell autonomy or transcript origin.

Fixed Hotspot membership and cohort-removal stability are useful protections.
They do not make discovery and testing within the same atlas independent.
Evidence: [ambient table](../../Analysis/SingleCell/candidates/ambient-program-recalibration-all117-candidate-2026-08-15-v2/results/ambient_program_effects.tsv)
and [lineage follow-up](FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md#content-lineage-follow-up-and-replication).

### 4. Continuum and network results support conditional tissue covariance

Signature exclusion reduces direct arithmetic circularity, but continuum and
program scores still share RNA measurements. The two co-primary continuum
scores correlate approximately **0.999**. Agreement between them is a
consistency check, with little independent corroboration. Grouped histology
also leaves residual true-severity variation within categories. An exact
sample-ID overlap check cannot prove biological independence across studies.

Figure 4H's **18 maxT** outlines test residual association under restricted
sample permutation. The competitive-background tests answer whether the
particular systems are exceptional compared with matched alternatives:
**0/43** pass the connected null and **1/43** the uniform null. Neither count
establishes intercellular regulation. GSE193066's unresolved mapping prevents
its current use as independent-person network replication.

The continuous-axis benchmark's negative controls are scientifically relevant
warnings about co-expression and ordering claims. They do not logically refute
every fixed external-score association under a different estimand. Its own
report also records limitations in the benchmark decision rules. Current
wording should remain cross-sectional and conditional.

Evidence: [claim disposition](FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md),
[signature-excluded producer](../../scripts/analysis/histology_anchored_continuum/molecular_layers/31_hotspot_continuum.R),
[exact-ID check](../../scripts/analysis/histology_anchored_continuum/03_build_projection_rf_axes.R),
[benchmark interpretation](../../scripts/analysis/continuous_axis_benchmark/11_assemble_report.R).

### 5. The other review identifies real genetic evidence issues, with narrower consequences

The **C8 biological-anchor requirement is real and still used** by the
September combine path. It requires selected ABF anchors to remain above or
below 0.5 and also checks EXOC3L4's SuSiE result. Such biological expectations
must not determine which scientific outcomes are admissible: protein-altering
mechanisms and expression-QTL associations are not mutually exclusive.
However, the requirement does not alter PP.H4 values, and the code writes the
gene table **before** stopping. The other review's statement that a failing run
“never produces a combined table” is therefore incorrect. The September log
records that the requirement passed; no unguarded sensitivity was run in this
review.

Evidence: [September combine script](../../GWAS/finemapping/src/perf/07_combine_rerun.R),
[execution script](../../GWAS/finemapping/src/perf/11b_assemble_coloc_r2.sbatch),
[saved combine log](../../GWAS/finemapping/results/susie_coloc_r2_20260924T053244Z/qc/combine.log).

The fine-mapping aggregator also really takes maxima across method-specific
PIPs and unions membership fields. A maximum of posterior probabilities from
different models is not itself a calibrated posterior probability, and the
union is not an intersection of method support. This concern reaches the
Figure 2B producer, which reads `recommended_pip` and takes a further maximum
with `susie_pip`. The August panel was inspected as PDF text and as a rendered
image, and its generating script was traced. The number of plotted points
whose value changes under a method-specific definition was **not recomputed**.
The scientifically necessary correction is an explicit method-specific
probability/estimand, not an assertion that all fine-mapping results are false.

Evidence: [PIP aggregation](../../GWAS/finemapping/src/04_aggregate_results.R),
[Figure 2B producer](../../scripts/figures/fig2_resolution_plane.R),
[plotted source](../../figures/main/fig2_genetics/panels/Fig2B_PIP_vs_SuSiE-coloc_source.csv).

LD mismatch and ancestry-mismatched molecular-QTL coverage are valid limits on
fine-mapping resolution and interpretation of absence. The 77/163 historical
orientation errors cited by the other review refer to a specific older
`seqfunc_v2` release. Current code contains revised abstention and ridge
requirements; those old counts do not demonstrate the same current error rate
in every COLOC arm. Population sample size alone also cannot prove that every
cross-ancestry resolution gain is unsupported; the actual locus-level
information and calibration must be examined. [SuSiE's own reference-mismatch
example](https://stephenslab.github.io/susieR/reference/rss_mismatch_example.html)
illustrates why this requires explicit checking.

### 6. Historical score calibration problems are real, but are not the current Catalog definition

The old clinical tiering explicitly includes a trial-positive term. It cannot
support “using no clinical information.” The old Govaere signature comparison
also shares a source cohort with its bulk-derived score component. The LOMO
helper has a real `S2_intact`/`S2_coloc` mismatch and labels a top-PP4 proxy as a
nearest-gene baseline. A leave-GSE135251-out implementation exists, but changing
to its reduced S1 definition requires a matched comparator; merely replacing
one AUROC does not validate the original full score.

Current manuscript methods retain the old score as supplementary/navigation
and report its failed strict clinical benchmark. Current evidence classes are
assigned from specified genetic and RNA calls, not the trial-positive tier.
The historical leakage therefore cannot be generalized to every current
Catalog class. The Catalog's RNA-dependent classes still await the corrected
bulk analysis for the separate reasons above.

Evidence: [clinical term](../../RNA-seq/40_drug_target_validation.R),
[LOMO helper](../../RNA-seq/46d_baselines_and_lomo.R),
[leave-cohort-out implementation](../../RNA-seq/scripts/wave1_d_lomo/lomo_validation.R),
[class assignment](../../scripts/manuscript/build_evidence_class_release.R),
[current score limitations](../manuscript/METHODS.md#historical-multi-evidence-score).

## Disposition of the supplied “what is broken” findings

The identifiers below match `02_WHAT_IS_BROKEN.md`. This table preserves the
concerns without treating the August draft as the current manuscript.

| Item | Assessment | Combined disposition |
|---|---|---|
| 1.1 Govaere calibration leakage | Confirmed historical concern | Cohort reuse and LOMO defects matter for that calibration. The current paper no longer relies on the old convergence ranking as its central deliverable. Existing leave-cohort-out code is not automatic validation of the original score. |
| 1.2 Biological-anchor requirement and palindrome handling | Confirmed requirement; qualified conclusion | It can condition acceptance on expected biology, but does not edit values and runs after output writing. Current palindrome handling differs from the unconditional historical exclusion; the old PNPLA3 number is not a current exact-pair result. |
| 1.3 Clinical outcome in clinical-free tiering | Confirmed historical concern | The C4 term makes the old wording false. It does not enter the current genetic/RNA evidence-class rule. |
| 1.4 Library counts reported as donors | Confirmed historical problem; supplied update stale | Use 104 descriptive, 102 analyzed, and 64 focal native-model donors. Current calls are two ordinary and zero HC3, not two robust plus two fragile. A calendar date alone cannot establish that every earlier p-value was pseudoreplicated. |
| 1.5 Two programs may reflect contamination | Qualified | Ambient/source concerns are real. Neither paired correction-change test is supported, so contamination is not established as the unique explanation. |
| 1.6 Stage may be a dataset axis | Confirmed for the cited historical Progressor analysis | The 99.94% cell fraction is not the design of the current donor-level 117-program model. Current within-cohort support still needs careful interpretation. |
| 1.7 Wrong GSE130970 controls | Confirmed; corrected fit already exists | The historical 1,347-call fit was corrected to 1,725 calls. Identity and counting correction remain unfinished. Not every analysis reads the affected control label. |
| 2.1 Benchmark contradicts progression narrative | Valid warning; overgeneralized | Ordering and shared-expression concerns matter. One negative benchmark under a different estimand cannot refute all external-score associations. Longitudinal claims remain unsupported, and GSE193066 is unqualified. |
| 2.2 Multiple evidence channels share sources | Confirmed historical interpretation problem | Derived channels are not independent modalities. Low rank correlation cannot prove independent errors or mechanisms. A nominal minimum across TWAS tests is not a calibrated gene-level call; the exact false-positive count is unverified. |
| 2.3 Fixed-effect pooling under heterogeneity | Qualified; automatic invalidity rejected | The cited I² value is an inherited historical comment. A pooled model with cohort intercepts differs from a fixed-effect meta-analysis. Heterogeneity limits a common-effect interpretation, but does not alone prove anti-conservative p-values or invalidate all effect directions. |
| 2.4 NMF instability and labels | Confirmed historical assignment instability; qualified for current axes | The saved legacy refit has ARI 0.295; a separate v2 multi-seed comparison has different basis-stability metrics and favors k4. These are different caches and measures. Fixed continuous axes can be retained with limitations; neither supports stable patient subtypes. |
| 2.5 Non-independent cohorts/donors | Confirmed concern; each use needs qualification | Genetic/source identities remain unresolved in primary bulk correction. Expression fingerprint similarity is an inference about reuse, not independently verified participant identity. Current COLOC aggregation already counts overlap-aware GWAS groups separately from raw study counts. |
| 2.6 LD/PIP and ancestry claims | Partly confirmed; historical rates not globally transferable | Maximum-PIP aggregation is a real interpretation issue. LD and molecular-QTL ancestry mismatch are real limitations. The historical orientation/ridge failures cannot be applied without checking the actual current arm and outputs. |
| 2.7 Small-n/spatial/protein examples | Mixed historical concerns; several exact claims unverified here | Sample size, statistic identity, and selection must be preserved. Old plasma and small-n spatial statements cannot be reused. The full set of old claimed numbers was not reconstructed from current outputs in this review. |
| 2.8 Four-level regulatory cascade | Unsupported as a current Resource claim | Do not revive the cascade. The cited 46.3% concordance and zero four-way calls were not independently rederived here; 46.3% alone does not establish a result below an appropriate null. |
| Tier 3 Plans 40–44 | Retain as terminal supplementary boundaries | Failure to confirm and insufficient precision are different from refutation or a demonstrated negative. Current Catalog strict tested-negative count is zero. |
| Tier 4.1 Numbers differ between drafts | Qualified | Different populations, definitions, and versions can yield different correct numbers. The problem is using incompatible numbers in one scientific claim, not the existence of historical differences. |
| Tier 4.2 Different genetic denominators | Qualified | ABF, SuSiE, their union, trait subsets, and jointly testable universes legitimately differ. Current reporting must name the denominator; the supplied old numbers are not current values. |
| Tier 4.3 Authored verdict constants | Valid caution | Constants are not independent numerical derivation. Some are legitimate expectations or scientific decisions; their existence does not establish that upstream computed results are fabricated or wrong. |
| Tier 4.4 Easy overlap/circular endpoint | Valid conceptual criticism | A base-rate overlap calculation or class-defining endpoint cannot establish independent biology. Current scope explicitly avoids inferring opposition from sparse overlap. |
| Tier 4.5 Post hoc rescue | Historical concerns require tracing | Changing thresholds after inspection affects confirmatory interpretation. The listed chains were not all traced to current results; no blanket invalidation of current testing families follows. |
| Tier 4.6 Selected examples | Confirmed limitation | The protein trio and focal programs are selection-conditioned examples. Current text already narrows protein confirmation and separates discovery, content lineage, and disease testing. |
| Tier 4.7 Listed code defects | Mixed; do not accept the list without checking use | Several plausibly affect old outputs. Each needs a current producer/consumer trace. The independently confirmed RNA code defects are listed below; other regex, length-offset, and null defects are unverified here. |
| Tier 4.8 Almost no results exist | Incorrect for this workspace | The other environment lacked result tables. HPC contains saved gene tables, donor effects, correspondence tables, numerical comparisons, and figures; their presence is not itself scientific validation. |

For item 2.3, the scientific target determines whether a common-effect or
random-effects summary is appropriate; model choice should not follow an I²
threshold alone. See [Cochrane's discussion of heterogeneity and model
interpretation](https://www.cochrane.org/authors/handbooks-and-manuals/handbook/current/chapter-10).

## Additional code findings retained from the RNA review

| Code path | Scientific consequence | Qualification |
|---|---|---|
| [Portal single-cell exporter](../../scripts/portal/generate_singlecell_data.py) | Reads the historical all-cell DE table instead of the donor-collapsed path. | Which public build contains that output was not checked. |
| [All-cell pseudobulk producer](../../Analysis/SingleCell/scripts/314_allcell_pseudobulk_de.R) | Its `~ condition` contrast remains dataset-confounded; captured/sorted fractions do not reconstruct natural tissue proportions. | The producer itself warns against using it as a validated disease axis. |
| [Lineage pseudobulk helper](../../Analysis/SingleCell/scripts/pseudobulk_de_donorcollapse.R) | The ≥50-cell filter counts all donor cells rather than cells in the tested lineage. | Historical/helper use; not established as an input to current 117-program effects. |
| [QC doublet handler](../../Analysis/SingleCell/scripts/v2/00_qc_harmonize_per_dataset.py) | A Scrublet exception becomes zero doublet scores. | Code vulnerability confirmed; an actual failed run was not demonstrated. |
| [Legacy ALT folding](../../scripts/analysis/alt_contig_recount/fold_lib.R) | Shared gene names and rounded sums do not establish locus equivalence. | Current authorities already reject it as the primary correction. |

The bulk canonical rule also combines a zero-null test with an observed
|logFC| >0.5 filter. That does not test a minimum true effect of 0.5. TREAT at
0.25 is a distinct test. BH within each contrast likewise does not give a
joint 5% guarantee for an “any stage” union. See [limma's reference
manual](https://bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf).

## Assessment of the proposed top-20 execution plan

The archive usefully emphasizes source identity, explicit contrasts, null
comparisons, power, and stopping rules. Its numerical thresholds and estimates
are proposals. They are not newly measured facts or prospectively registered
tests merely because they are written in this archive. The README also says
the JSON files were not updated with the Markdown revisions; the two cannot be
used interchangeably as an execution specification.

The following changes are necessary before considering any of that plan:

1. **Use the current identity/counting decision.** Item 1 retains unknown
   identities in primary inference, contradicting the user's current exclusion
   decision. It also targets an old 1,347-call/23,370-row baseline and omits the
   unfinished ALT correspondence. Current corrected-control and identity-only
   baselines are different. A similarity criterion cannot justify keeping
   obsolete numbers after a correction.
2. **Remove unjustified replication arithmetic.** Item 1 estimates a roughly
   1.5% standard-error change from residual degrees of freedom alone. Correlated
   repeated rows also change covariance and weighting; their effect is not
   determined by that degrees-of-freedom calculation. Repeated observations
   can also change point estimates when participants are weighted unequally.
3. **Do not estimate reader error directly from serial biopsies.** Items 15
   and 18 use or propose reliability values from repeat-biopsy differences.
   Such differences contain true change and sampling variation as well as
   reading error; GSE193066 additionally lacks qualified identity/source
   mapping. Even the 29 same-cluster pairs do not solve those problems.
   Classical additive-error formulas are not automatically Spearman bounds for
   grouped ordinal stages. Use the actual score/label definition and an
   explicit conditional error model; report raw correlations separately.
   Current RESULTS already contains raw correlations, tie ceilings, and
   conditional simulated ceilings.
4. **Keep participant linkage independent of the endpoint.** Item 5 properly
   prohibits a positional plasma/RNA join. Biological concordance may check
   an independently documented mapping; it cannot manufacture or select the
   mapping and then serve as independent evidence of coupling. No plasma
   modeling should proceed without authoritative participant correspondence.
5. **Keep the Resource's scientific endpoint separate from screen stability.**
   Items 6, 8, 9, 19, and 20 include mouse/Cas13 design, targets, or screen-slot
   retention. Those are separate follow-up questions under current scope.
   A corrected Catalog remains useful if its old routed list changes or its
   concordant class is empty. Catalog validity is not determined by retaining
   20 of 30 old screen slots or meeting a Jaccard threshold.
6. **Avoid unnecessary new work.** Several proposed corrections already have
   candidates or narrowed claims. Plans 40–44 are terminal supplementary
   analyses. The archive's approximately 56 postdoc-week total is an estimate,
   not a required Resource program. Controlled-access applications, additional
   readers, and wet-lab studies are separate decisions. Dataset-lead
   availability and sample counts were not independently rechecked here.

The parallel-track document also says all models are at their label ceiling.
That does not follow from a tie ceiling or an assumed reliability scenario.
Its displayed 86-donor multiome breakdown has 23 MASH donors but a stated
19+16=35 MASH fibrosis split; the denominator needs resolution before use.
Claims that one resource is the only route to an outcome or label-error test
are unverified and should not determine the paper's dependencies.

## Proportionate next actions for the current paper

The existing approved correction sequence remains appropriate: resolve the
selected identity/counting inputs; fit the corrected bulk and dependent
stage/composition analyses; reassess the actual effects and uncertainty; then
reconcile figures and Catalog fields. Preserve failed tests and unavailable
states rather than trying to recover a preferred count.

The revised NMF protocol already keeps **k4-P2 and k6-P6 unavailable** for
stage use after failed historical reconstruction criteria. The remaining
historical checks do not establish corrected-count transport. Keep NMF as
continuous axes and retain the distinction between fixed-score reproduction
and discovering a stable biological decomposition.

The newly identified biological-anchor acceptance rule, PIP interpretation,
and RNA exporter/QC/helper problems should be considered individually before
implementation. This review does not approve changes to methods, fixed
memberships, scientific thresholds, or result adoption.

## Verification, boundaries, and source identity

This combined review checked archive membership, compared the two copies of
the supplied findings, inspected the current authority documents, traced the
clinical/LOMO and September COLOC code paths, inspected saved NMF comparisons,
read the current donor-by-diagnosis table, and inspected Figure 2B as text and
a rendered image. It carries forward the earlier RNA review's independent BH,
call-count, annotation-disposition, and Figure 4 checks. Those earlier checks
are not presented as new model fits in this turn.

No complete raw-data audit, corrected recount, differential model, permutation
analysis, genetics rerun, or external-dataset acquisition was performed. No
analysis code or scientific result values were changed. The full inventory's
362 reconstructed entries, 166 candidates, 20 plans, and external dataset
leads were inspected structurally and selectively; this is not independent
validation of every hypothesis or citation in those large files.

Input SHA-256 values identify the reviewed copies:

| Input | SHA-256 |
|---|---|
| `Hypothesis_Review_2026-10_rev1.tar.gz` | `0bf87f50bb4c4f60731bfb8e250eabd0853f96b4f6b3a81adae3c20f09d83953` |
| `02_WHAT_IS_BROKEN.md` | `8154eef6a871d9251937ca0c261939ff40995dd9f99497a2e5163337f2387695` |
| Earlier RNA review | `a10ca1bdc004912a71245bd047c8a3ac7e9fa9b2cf79c55bdf1d72e23dc026f3` |
