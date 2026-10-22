# Focused statistical rigor integration: execution assessment

Date: 2026-09-08. Companion to the [approved handoff](OAKES_PACHTER_RIGOR_INTEGRATION_HANDOFF_2026-09-08.md).
The five paper authorities retain control of scientific scope and result selection.
This assessment does not adopt a release or declare the manuscript ready.

## Findings and changes

The selected fragment QC and all 67 inspected complete BH families reproduce.
The frozen severity scorer and three additional transfer correlations reproduce.
An existing stage-table defect omitted every confidence interval in 11 contrasts:
257,070 gene–contrast rows. A bounded repair replays those exact fits and fills
intervals only after reproducing the original effects, standard errors, t
statistics, P values and FDR. Original result files remain unchanged. See the
separate repair output and checks below; missing intervals in the original
remain an expected failed check, not an assertion that the source was correct.

Working Methods now distinguishes selected and historical denominators,
fragment QC recomputation from subset reuse, reported from modeled sex,
within-family BH from adaptive selection, and gene-level from isoform inference.
Discussion and Figure 7 prose qualify contamination, identity, model exposure,
transductive comparisons and uncertainty. STATUS and RESULTS correct two
additional external cohorts being described as three and remove an unsupported
equivalence interpretation. Fixed gene thresholds, 117 program definitions,
COLOC memberships and the separate follow-up paper's scope are unchanged.

## What was learned from the source paper

[Oakes and Pachter, 2026](https://link.springer.com/article/10.1186/s12863-026-01440-8),
its [analysis repository](https://github.com/pachterlab/OP_2025), and its
[public peer review](https://link.springer.com/article/10.1186/s12863-026-01440-8/peer-review)
motivated explicit identity/QC checks, attention to transcript resolution,
task-specific model comparisons and separation of observed expression from
physiological explanations. Their demonstrated analyses are evidence about the
rat experiment, not validation of our human data. Generalization beyond a
particular predictive task requires additional evidence.
Author reputation is not a substitute for checking assumptions or outputs.

The adaptations below are our own assay-specific safeguards, not a claim that
the source paper performed every procedure in this assessment. Its UMI handling,
transcript P-value aggregation and isoform analysis were not transplanted to
human gene-count data. Log-expression predictors are valid in linear regression;
their transformation changes the modeled relationship. A particular scVI
comparison cannot establish universal linear-model superiority.

## Evidence assessment

“Verified” below means the specified output-level check, not independent
reproduction of an entire analysis or source study. Disposition describes the
finding before the accompanying text changes.

| Safeguard and scientific purpose | Producing evidence and available result | Claim affected and disposition |
|---|---|---|
| Sample-key alignment: avoid multiplying observations in joins | Selected metadata, QC and count columns have unique keys; pooled count order agrees with saved design; donor/program scores join the donor roster without inflation | **Verified but inadequately reported.** Sample-key consistency is not genetic identity |
| Genetic identity and participant overlap | No genotype-based confirmation or complete cross-study identity crosswalk established in inspected inputs; model card lacks a complete training donor key | **Required check missing or not reproducible** for a genetically identity-confirmed or paper-wide independent-person claim. Such claims are withheld |
| Counting units and QC provenance | F_legacy manifest: fragments, QC recompute, 1,257 records; F_five: fragments, reuse_F_legacy, 844 records; existing reverse-strand audit described in Methods | **Verified but inadequately reported** for mode/denominators. No raw-read recount or independent strand audit in this task |
| Replication | Bulk selected analysis declares one biopsy per participant; actual keys and group eligibility checked. Fixed-program scores have unique donor/program keys. Chromatin 99 participant keys map one-to-one to five folds | **Verified but inadequately reported**, conditional on deposited identity. Do not call libraries, genes, cells or repeated biopsies independent people |
| Confounding and missingness | Pooled model full rank 7/7; stage memberships and F0 arms reconciled; PCA/sex/library flags retained; no new relabeling or exclusion | **Verified but inadequately reported.** Healthy/disease F0 mixture remains material; technical versus biological cause of each PCA outlier is unresolved |
| Mean–variance model, effects and uncertainty | Selected bulk uses limma-voom quality weights; pooled t = effect/SE. Stage source intervals all missing; later F0-arm intervals finite | **Required check missing or not reproducible** in the original stage CI table; bounded interval repair supplied separately. Full pooled refit and all sensitivities not rerun |
| Multiplicity | 67 complete saved families agree with BH, including point-null and TREAT, stage/F0, ambient/HC3, lineage and historical cell-type pseudobulk | **Verified but inadequately reported.** No global FDR guarantee from choosing a best cell type or protein assay |
| Prediction and exposure | Four transfer correlations rederived; 23 model/identity checks pass. Additional native baseline uses full-cohort expression scaling. Native individual predictions absent from inspected producer outputs | **Verified but inadequately reported** for scores; **required check missing or not reproducible** for an independently reproduced paired baseline interval or strictly inductive superiority claim |
| Contamination | Aligned existing program sensitivity: five corrected datasets; two failed decontX datasets; 106 transport-testable programs and 11 untestable T-cell programs. Raw/corrected/delta and HC3 families checked | **Verified and adequately reported** in Methods; stale Discussion wording corrected. Not proof of cell autonomy or transcript origin |
| Transcript resolution | Selected input is gene counts; primary isoform claims explicitly excluded | **Not applicable or outside this integration** for new transcript discovery/usage tests. Gene-level terminology and limitations clarified |

## Supplementary denominators and identities

Primary table: [supplementary_analysis_denominators.tsv](oakes_pachter_rigor_20260908/tables-v4/supplementary_analysis_denominators.tsv).
One row per cohort, assay and analysis; includes zero-contribution cohorts so
eligibility is visible. This is not a sum-across-analyses participant registry.

| Selected nine-cohort source | Count-matched records | Technical-QC records | Pooled records |
|---|---:|---:|---:|
| GSE126848 | 56 | 53 | 53 |
| GSE130970 | 78 | 76 | 76 |
| GSE135251 | 216 | 214 | 214 |
| GSE162694 | 143 | 140 | 140 |
| GSE213621 | 367 | 361 | 361 |
| GSE167523 | 98 | 96 | 0 |
| GSE174478 | 93 | 93 | 0 |
| GSE193066 | 164 | 160 | 0 |
| GSE240729 | 66 | 64 | 0 |
| Total records, not independently verified unique people | 1,281 | 1,257 | 844 |

The source metadata has 1,469 records: 185 belong to the excluded PRJNA512027
cohort, leaving 1,284 within the nine-cohort scope and three unmatched records.
The excluded 185 are not failed technical QC. Pooled records comprise 157
controls and 687 disease records. The same cohorts may be reused across assays
or model tasks; assay totals must not be summed as unique people.

Table interpretation and reproducing identifiers:

- `source_records` counts metadata rows, not necessarily distinct donors.
  `count_matched_records` counts keys with selected counts; `technical_qc_records`
  counts the PCA/library union rule. `final_analysis_records` counts the actual
  selected membership. `unique_participants` for bulk follows the selected
  one-biopsy declaration, not a new genetic check; it is NA for the broad QC
  inventory where repeated biopsies may occur. Program values use the existing
  biological-donor roster. Unknown upstream cell-assay totals remain NA.
- [bulk_sample_qc_membership.tsv](oakes_pachter_rigor_20260908/tables-v4/bulk_sample_qc_membership.tsv)
  and [bulk_analysis_membership.tsv](oakes_pachter_rigor_20260908/tables-v4/bulk_analysis_membership.tsv)
  contain the identifiers behind every bulk denominator. Program membership is
  in [program_donor_score_membership.tsv.gz](oakes_pachter_rigor_20260908/tables-v4/program_donor_score_membership.tsv.gz);
  apply the analysis-universe filters in the replay script, not a sum across programs.
- [QC flags by clinical group](oakes_pachter_rigor_20260908/tables-v4/qc_flags_by_clinical_group.tsv)
  separates overlapping PCA, library and sex flags. Technical failures equal
  PCA failures + library failures − their intersection. `technical_pass_not_in_analysis`
  is a subsequent eligibility loss, not another QC failure. Clinical groups and
  stage-specific memberships remain analysis-specific.
- [Sex annotation summary](oakes_pachter_rigor_20260908/tables-v4/sex_annotation_summary.tsv):
  632 records have reported sex and 649 lack it; 24 disagreements are flagged,
  including six retained in the pooled analysis. All matched records have
  inferred sex. The selected model uses `inferred_sex`, irrespective of the
  `sex_source` or `sex_final` display field. No identities or exclusions changed.
- [Program testability and results](oakes_pachter_rigor_20260908/tables-v4/program_testability_and_results.tsv)
  retains all 117 programs in both universes, including the 11 untestable
  programs omitted from fitted-membership denominators. GSE189600 and Liver
  Atlas are absent from the common corrected universe and explicitly uncorrected
  in the full-universe sensitivity. Missing or failed is not negative.
- Model evaluation populations have their own
  [participant membership table](oakes_pachter_rigor_20260908/models-v2/evaluation_membership.tsv)
  and [99-participant paired-assay folds](oakes_pachter_rigor_20260908/models-v2/chromatin_participant_folds.tsv).
  Upstream model QC and training donor counts not independently established here
  must not be reconstructed from these evaluation-only tables.

## Statistical checks and bounded repair

The [67-family table](oakes_pachter_rigor_20260908/tables-v4/statistical_family_checks.tsv)
reports family size, finite tests, missingness agreement, BH discrepancy and
significant counts. All pass the absolute 1e-10 comparison. Families comprise
two pooled, 11 original stage, 23 F0-sensitivity, 12 ambient, eight lineage and
11 historical pseudobulk families. The historical `_3dataset` pseudobulk files
are not asserted to be the adopted atlas's primary effects merely because
their BH calculations reproduce. Their exact producing run still needs linkage.

The pooled design is `~ dataset + inferred_sex + group_binary`, with Control
as the comparison reference and the saved 23,370-gene family. Primary selection
is q < 0.05 and |log2FC| > 0.50 (1,347 records); TREAT tests the interval null
at 0.25 with its own complete-family BH adjustment (1,616 records). These are
different statistical guarantees, not two interchangeable ways to label the
same evidence. Stage fits use the same assay-native voom model with the stated
contrast and available cohorts. No age/BMI adjustment is implied where those
terms are absent. Complete-case stage eligibility and inferred-sex modeling are
visible in the membership tables; no clinical imputation was added.

The F0 reference is 31 controls plus 73 disease biopsies. The existing B arm
removes only controls at F0, not controls at F1/F2; the C arm adds control
status. The C F3→F4 branch cannot estimate that term. The existing corrected
GNMT/MAT1A interpretation in RESULTS remains binding. The
[example sensitivity rows](oakes_pachter_rigor_20260908/tables-v4/f0_example_sensitivity.tsv)
retain all arms, rather than selecting a preferred P value.

The older stage producer computed a vector of t critical values without gene
names and then indexed it by gene identifier, yielding NA intervals. The fix
adds those names; it does not modify the fit. The separate
[interval candidate](oakes_pachter_rigor_20260908/stage-interval-repair/stage_extension_all_gene_results.tsv.gz)
and [11-contrast reproduction checks](oakes_pachter_rigor_20260908/stage-interval-repair/reproduction_checks.tsv)
retain all 23,370 genes per contrast and require each original effect, SE, t,
P and FDR to agree within 1e-8 before writing. These are pointwise moderated-t
intervals, not simultaneous confidence bounds. Repairing them neither removes
F0 confounding nor changes any gene's significance. No panel was regenerated
and no original result overwritten.

Existing cohort-omission, composition, equal-weight and ambient/HC3 sensitivities
remain the relevant sensitivity analyses. This task recomputed the listed
ambient/lineage BH families and eligibility, not every sensitivity coefficient.
No inference that a technical outlier is biologically invalid was added.

## Model evaluation and remaining evidence

The [model-evaluation summary](oakes_pachter_rigor_20260908/models-v2/model_evaluation_summary.tsv)
contains 13 task rows. Four transfer metrics were independently rederived;
the other rows explicitly retain source-reported evidence rather than invented
verification. [Model checks](oakes_pachter_rigor_20260908/models-v2/model_checks.tsv)
report 23/23 passing identifier, scorer and metric checks.

| Evaluation | n, source-declared participant records | Reproduced raw Spearman | Ratio to tie ceiling |
|---|---:|---:|---:|
| GSE268273 reference, frozen v1.2 scorer | 109 | 0.6631190837 | Not newly calculated here |
| GSE213621 disease-only three-level outcome | 299 | 0.6381976612 | 0.6772633146 |
| GSE213621 control-inclusive sensitivity | 367 | 0.7472935825 | 0.7736379732 |
| GSE276114 three-level outcome | 177 | 0.6738645294 | 0.7943454354 |

The model card's 521 training records are samples, not a verified donor count.
GSE193066's 58 paired participants are not 116 independent evaluation people.
The model card changed while this task was paused. Its updated section 4.2
reports only three pairs with both biopsies inside the shipped input bounds;
the ICC is consequently an out-of-contract sensitivity, not accepted-input
deployment performance. This source-reported qualification was added to Methods,
Figure 7 and RESULTS; no ICC replay was performed here.
GSE268273 is external to the severity fit but its 109 records occur in
BulkFormer pretraining; those are different exposure statements. GSE213621 is
also a Resource cohort, so it is not paper-wide untouched data. The chromatin
comparison is a development-exposed cohort, with 99 unique participant keys
in five folds; lack of age and BMI limits confounding assessment.

In the additional severity comparisons, the recorded native producer's
`native_oof` passes the entire expression matrix as scaling reference. Outer
test labels are not used for fitting that fold, but test expression contributes
to scaling; outer-training feature selection/PCA are not refitted in each
inner tuning split. This is transductive evaluation. No strictly inductive
equivalence/superiority claim is retained. Paired source intervals use 10,000
stage-stratified draws, seed 213621. Native individual predictions are not
serialized by the inspected producer, so those intervals were not independently
rederived. No replacement 2,000-draw bootstrap was warranted: existing valid
conditional intervals are preserved, while unavailable paired arrays prevent
a new paired calculation. A future missing-interval calculation should use
2,000 paired participant draws, seed 20260908, grouping repeat biopsies, and
state that the intervals condition on fitted models.

Unfinished checks, affected claims and missing evidence:

| Check not run / unresolved | Consequence for interpretation | Evidence needed to finish |
|---|---|---|
| Genotype-based sample identity and complete cross-study donor overlap | No genetic identity-confirmed or all-assay unique-person claim | Authorized donor crosswalk or consent-compatible genotype concordance outputs; expression similarity cannot substitute |
| Strictly inductive native baseline with nested preprocessing | Existing paired deltas describe the recorded transductive comparator only | Fit and serialize fold-native preprocessing and participant predictions in a new, explicitly scoped comparison; not a silent replacement |
| Independent paired baseline intervals | No new validation of the source aggregate intervals | Recover native participant predictions, task outcomes and split keys; source execution is `Analysis/MASLD_Model_Benchmark/executions/external-validation-expansion-20260901T190000Z/external_validation.py` |
| Chromatin skill/interval replay from individual losses | Pairing/folds checked; numerical skill remains source-reported | Participant-level out-of-fold losses/predictions from `executions/chromatin-state-thesis-20260906T124724Z/run_a2.py`; source output retains aggregate results |
| ICC, activity-head, gene-absence, TabICLv2, late-fusion and DNA-direction numerical replay | Prior evidence retained; this review does not certify every Figure 7 row | Sources are indexed separately in model summary; retrieve task-specific paired outputs and repetition/split keys before recalculation |
| Every cohort-omission/composition coefficient and current atlas pseudobulk provenance | Existing sensitivity records are reused, not newly independently reproduced | Exact producing inputs and saved complete families for each current claim; historical pseudobulk BH alone is insufficient |
| Raw-read strand/count reprocessing and isoform analyses | No new claim at those resolutions | Outside this integration; not run |

These missing inputs do not have a safe executable replay command that can be
invented from aggregate metrics. Do not run the historical producers against
their original output directories to fill the gaps. The exact available
reproduction commands below cover all checks actually executed by this task.

## Reproduction, scope and verification

Executed in the existing allocation 21327089 on compute node ne1dc4-015,
using the existing rnaseq environment and one numerical thread. No new jobs,
dependencies, GPU training, raw reads, branches or commits. Scripts require a
new output directory and refuse overwriting existing outputs. Reproducing uses
new destination names; the following commands match the completed executions
with destinations renamed to avoid overwriting them:

```bash
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/manuscript/oakes_pachter_rigor/audit_tables.R "$PWD" docs/technical/oakes_pachter_rigor_20260908/tables-replay
PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python3 scripts/manuscript/oakes_pachter_rigor/audit_models.py --root "$PWD" --out docs/technical/oakes_pachter_rigor_20260908/models-replay
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/manuscript/oakes_pachter_rigor/repair_stage_intervals.R "$PWD" docs/technical/oakes_pachter_rigor_20260908/stage-interval-replay
python3 scripts/manuscript/validate_resource_scope.py
PYTHONNOUSERSITE=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python3 scripts/manuscript/oakes_pachter_rigor/verify_outputs.py --root "$PWD" --out docs/technical/oakes_pachter_rigor_20260908/verification-replay
```

Input SHA256 inventories and R sessionInfo/Python environment records accompany
each output set. Table audit uses seed 20260908; interval replay preserves seed
42 from the original stage workflow. No new stochastic model fitting or
bootstrap was performed. `tables-v4` and `models-v2` are the assessment tables;
`tables-v1` through `tables-v3` are intermediate task outputs, not manuscript
evidence. `models-v1` records the initial review; `models-v2` repeats it against
the updated model documentation. The first final verification detected changed
MODEL_CARD.md and MODELS.md hashes; these were reread, not silently accepted as
unchanged. Numerical input and stage-source hashes matched. Original stage
missing-CI checks are intentionally preserved.

The manuscript edits are in `docs/manuscript/METHODS.md`,
`docs/manuscript/draft/fig5_discussion.md`, `docs/manuscript/draft/fig7.md`,
`docs/STATUS.md` and `docs/RESULTS.md`. Analytical edits comprise the four
replay scripts and the single gene-name assignment in
`scripts/figures/build_pi_stage_extensions.R`. The initial dirty working tree
was recorded and unrelated changes preserved. No figure changed, so figure
text extraction and visual checks were not applicable to this task.

### Final verification outcomes

- [Independent output verification](oakes_pachter_rigor_20260908/verification-v3/verification_checks.tsv):
  **24/24 passed**, including unique denominator keys, sequential losses,
  membership recounts, model evaluation sizes, unchanged stage statistics,
  finite repaired intervals, current input hashes and local assessment links.
- The supplementary table has **1,844 unique cohort/assay/analysis rows**;
  all 360 bulk rows recount from the saved membership table. Program-level
  eligibility was separately checked against donor-level scores by the R audit.
- The R audit has **355 passing checks and one expected source defect**:
  missing intervals in the unmodified original stage table. All **67 BH
  families pass**; maximum absolute BH discrepancy is 4.44e-15.
- All **11 stage fits reproduce**; maximum absolute discrepancy across the
  replayed effects, SE, t, P and FDR is 4.80e-14. The **257,070 repaired rows**
  have finite intervals. Independent serialization comparison found differences
  below 1e-16 in saved effects/P/FDR and exactly zero in SE. No significance
  membership or threshold changed.
- Refreshed model audit: **23/23 checks pass**, including all four transfer
  correlations. The other Figure 7 quantities remain source-reported as
  specified above, not newly reproduced.
- [Resource scope check](oakes_pachter_rigor_20260908/verification-v3/resource_scope.txt):
  **PASS**. Its initial retired-word finding was corrected before the final run.
  The final source hash checks also pass after reviewing the two concurrently
  updated model documents and producing `models-v2`.

The requested handoff, assessment, manuscript revisions, supplementary tables
and bounded interval correction are reviewable. Adoption and the explicitly
listed missing evidence remain outside this completed integration.
