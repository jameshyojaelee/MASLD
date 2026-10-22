# Scientific review of bulk and single-cell RNA-seq

**Review date:** 2026-10-07  
**Scope:** scientific validity of the bulk and single-cell RNA-seq analyses,
interpretations, code paths, and reported results. This is a review record, not
a new analysis release or a change to result readiness.

## Overall assessment

The work uses several appropriate and carefully implemented methods. In
particular, bulk RNA-seq is analyzed with voom and limma, the current single-cell
program analysis collapses sequencing runs to donors, multiple-testing families
are retained, and failed or untestable analyses are often reported explicitly.
These strengths support a descriptive atlas and conditional evidence for broad
disease-associated expression changes.

The current files do not yet establish corrected independent-participant bulk
DEG or stage totals. The two focal native single-cell disease associations are
fragile to HC3 standard errors, and the results do not establish hepatocyte
transcript origin or independent replication of the tissue coordination claim.
Some of these limits are already stated in the current paper authorities; this
review records their scientific implications and identifies additional code
problems without changing those authorities.

## Main findings

### Bulk identity and annotation limit exact gene counts

The historical bulk analysis uses 844 sample records. Genotype-based review
finds repeated libraries and source-label conflicts, so those rows cannot all
be treated as verified independent people. The newer 774-unit and 738-unit
rosters use supported source or conservative genetic identity clusters; they
are not independently verified people. GSE193066 remains excluded from primary
replication and repeat-biopsy inference because its source-person crosswalk is
unresolved. See [`RESULTS.md`](../RESULTS.md) and the
[identity disposition](FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md).

Identity-only fits reproduce **1,756** canonical annotation-row calls in the
774-unit roster and **1,807** in the 738-unit sensitivity. The primary fit
retains 1,572 calls from the historical 844-row model; its log-fold-change
correlation with that model is 0.988. This supports persistence of a broad
expression pattern under the roster sensitivity, but it does not isolate the
effect of removing repeated identities because filtering, cohort composition,
normalization, and quality weights also change.

Those fits retain duplicate counting targets and unresolved annotation rows.
The current annotation review has **855 unresolved correspondences** among
86,369 records. Therefore, the current 1,756 and 1,807 counts are not corrected
physical-gene DEG totals. The unresolved records do not establish that every
primary gene is wrong; they do block the exact corrected totals. The current
counter appropriately refuses primary use until annotation correspondence and
count evidence are complete:
[`count_verified_fragments.py`](../../scripts/analysis/alt_contig_recount/count_verified_fragments.py).

The underlying bulk model is a reasonable RNA-seq model: `~ dataset +
inferred_sex + group_binary`, followed by `voomWithQualityWeights` and limma
moderation. The pooled disease coefficient assumes a common adjusted effect
across cohorts. It is an association, not a causal effect or proof that each
cohort has the same response. The canonical rule (BH <0.05 and observed
|log2FC| >0.5) does not test that the true effect exceeds 0.5. TREAT is a
separate interval-null test at 0.25 and answers a different question.
Likewise, the 61 genes selected as significant in any of four adjacent-stage
contrasts do not form one jointly BH-controlled family; each contrast has its
own correction. These are valid display definitions when stated precisely,
not a jointly controlled “any stage” discovery set.

The newer fragment counter is conservative about complete NH/HI alignment
sets, mate integrity, and one common gene target across all reported
alignments. Five diagnostic pilots reproduced native unique counts exactly.
They do not substitute for the full corrected recount. A legacy helper,
[`fold_lib.R`](../../scripts/analysis/alt_contig_recount/fold_lib.R), folds
primary and alternate-locus counts by shared gene name and rounds the sum. A
gene-name match alone does not prove identity. Current authorities reject
that helper as the primary correction.

### The two focal single-cell disease programs are fragile

The corrected source-diagnosis refit uses **64 biological donors** across five
datasets; the atlas-wide cell count is descriptive and is not the inferential
sample size. It tests a linear, equal-step Healthy-to-Steatosis-to-
Steatohepatitis trend with dataset adjustment. This is a cross-sectional
diagnosis association, not longitudinal fibrosis progression.

An independent BH recalculation over all 117 programs reproduced the saved
results to numerical precision: two ordinary model-based BH calls and zero
HC3 BH calls.

| Program | Coefficient | Ordinary BH q | HC3 BH q |
|---|---:|---:|---:|
| Stromal ECM | 0.399 | 0.0496 | 0.1053 |
| Ductular injury | 0.412 | 0.0187 | 0.1053 |

The 117-row family and separate HC3 results are correctly implemented in
[`512_hotspot_v2_donor_refit.R`](../../Analysis/SingleCell/scripts/hotspot_modules/512_hotspot_v2_donor_refit.R)
and saved under
`Analysis/Multimodal_Program_Projection/candidates/program-context-v2-stage-corrected-candidate-2026-08-13-v3/hotspot/work/`.
The ordinary results are suggestive; they are not robust to the recorded
heteroskedasticity sensitivity. This is not evidence that both effects are
absent.

The original Hotspot programs were discovered on cell-neighbor graphs and
then tested in donor summaries from the same atlas. Leave-one-dataset-out
membership persistence is useful stability evidence, but it is not independent
disease replication or donor-balanced program discovery. The 117 fixed
memberships protect downstream comparisons from redrawing modules around a
desired result; they do not remove selection dependence from the original
atlas test.

### Cell content and ambient RNA limit origin claims

The original scoring atlas was not ambient-corrected. A donor-level content
audit finds **80/117** programs whose highest content lineage matches the
discovery label and **37/117** that differ. Thus discovery lineage is
provenance, not evidence of transcript origin.

In the matched ambient sensitivity, Stromal ECM changes from β=0.401 to
β=0.130; its corrected BH q is 0.312 and the paired correction-change q is
0.271. Ductular injury changes from β=0.396 to β=0.389; its corrected BH q is
0.129 and correction-change q is 1.00. Neither corrected association nor
correction-induced change passes BH. Ambient correction therefore raises a
real source concern without proving that ambient RNA explains the original
associations.

A targeted fibroblast ECM score comparison is positive in a small
correction-available sample: **26 donors**, β=0.329, 95% CI 0.098–0.560,
two-program BH q=0.0106. It supports a narrow projected-score association.
The corresponding cholangiocyte Ductular injury contrast is unresolved
(β=0.0848, q=0.449). Nonrandom correction coverage and broad-lineage substate
mixture remain possible explanations. Neither score proves cell autonomy.
See the [content-lineage disposition](FIG34_SCIENTIFIC_CORRECTIONS_2026-09-29.md#content-lineage-follow-up-and-replication).

### Continuum and coordination are conditional transcriptomic evidence

The continuum analysis removes the external signature genes from fixed Hotspot
memberships, which reduces direct arithmetic circularity. The remaining scores
and outcomes still use the same expression substrate, so this is evidence of
transcriptomic organization, not an independent assay validation. The two
co-primary continuum methods correlate approximately 0.999 in each evaluation
cohort, so their agreement is a consistency check rather than independent
corroboration. Grouped fibrosis categories also leave room for true severity
differences within category.

The historic GSE193066 coordination evaluation uses a source-person crosswalk
that current identity review does not qualify. It therefore cannot presently
establish independent-person replication. Figure 4H’s 18 maxT outlines test
residual sample association under their permutation scheme. They do not show
that these systems are unusually specific relative to matched alternative gene
sets; the separate competitive-background results are 0/43 for the connected
null and 1/43 for the uniform size-matched null. The two tests must retain
separate interpretations.

The revised NMF transport protocol keeps historical weights fixed and leaves
two stage axes unavailable (k4-P2 and k6-P6) after they miss the prespecified
reconstruction thresholds. The diagnostics check transport fidelity to
historical scores; corrected-count scores and downstream models have not run.

## Code issues and their current scope

| Finding | Consequence and scope |
|---|---|
| [`generate_singlecell_data.py`](../../scripts/portal/generate_singlecell_data.py) reads `pseudobulk_de/allcell_pseudobulk_de.csv`, while the donor-collapsed output is under `pseudobulk_de/donor_collapsed/`. | If this exporter is used for a current portal build, it exports a historical table. The donor-collapsed producer itself warns that its `~ condition` fit is dataset-confounded and must not be used as a disease axis. Whether a particular public portal build contains this table was not checked. |
| [`314_allcell_pseudobulk_de.R`](../../Analysis/SingleCell/scripts/314_allcell_pseudobulk_de.R) sums captured cells across preparation fractions and fits `~ condition`. | Donor collapse fixes pseudoreplication but not dataset/condition confounding or sorted-fraction sampling bias. This is not a valid independent whole-liver tissue comparison. |
| [`pseudobulk_de_donorcollapse.R`](../../Analysis/SingleCell/scripts/pseudobulk_de_donorcollapse.R) applies the ≥50-cell gate using each donor’s total cells across lineages. | A donor can pass a lineage-specific analysis using cells from other lineages. Bounded review found this helper in historical audits, not in the current 117-program producer; do not attribute the defect to that producer. |
| [`00_qc_harmonize_per_dataset.py`](../../Analysis/SingleCell/scripts/v2/00_qc_harmonize_per_dataset.py) catches Scrublet exceptions and sets doublet scores to zero. | A failure could be recorded as no doublets. No actual failure was demonstrated. The available QC audit has full metrics only for Liver Atlas; six other datasets are marked `skipped_existing`, limiting independent reconstruction of QC rates. |

## What is defensible now

- Bulk RNA-seq methods are appropriate for count data. Voom models the
  mean–variance relation and supplies observation-level precision weights;
  the limma model is a recognized RNA-seq workflow ([Law et al., 2014](https://doi.org/10.1186/gb-2014-15-2-r29)).
- Current single-cell program disease tests use donors, not cells, as
  replicates, and preserve the complete 117-program BH family. Donor-level
  pseudobulk or replicate-aware analysis is the appropriate basis for
  condition comparisons ([Squair et al., 2021](https://doi.org/10.1038/s41467-021-25960-2)).
- The atlas and fixed program registry are useful descriptive resources.
  Broad disease-associated RNA remodeling is supported conditionally, while
  exact corrected bulk counts, hepatocyte-origin claims, and independent
  coordination replication remain open.

## Review checks

The review independently recalculated the 117-program ordinary and HC3 BH
families, checked the two focal program values against their saved table,
recalculated the identity-only canonical DEG counts and BH values, and tallied
the r7 annotation dispositions. It also inspected extracted text and rendered
Figure 4C and 4H. No full corrected bulk recount, new differential model,
continuum model, coordination model, or QC run was performed. No analysis code
or result values were changed for this review.
