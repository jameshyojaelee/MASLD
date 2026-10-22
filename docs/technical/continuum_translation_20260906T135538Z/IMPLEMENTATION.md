# Continuum translation extension: implementation and scientific review

Implemented 2026-09-06 as a candidate for the Resource paper. This report
records the completed analysis, its limits, and reproducible commands. It does
not adopt a new main figure or change inherited evidence memberships.

## Deliverables

The completed result root is
`RNA-seq/results/histology_anchored_continuum/translation/hac-translation-20260906T230500Z/`.
The final integration is **`integrated-20260906T231500Z/`**. Intermediate
integration directories are not the completed deliverable.

- [Integrated PDF figure](../../../RNA-seq/results/histology_anchored_continuum/translation/hac-translation-20260906T230500Z/integrated-20260906T231500Z/integrated_therapeutic_hypotheses.pdf), with individual PDF panels and an external caption in the same directory.
- [Complete Gene Catalog sidecar](../../../RNA-seq/results/histology_anchored_continuum/translation/hac-translation-20260906T230500Z/integrated-20260906T231500Z/gene_catalog_translation_complete.tsv.gz): 28,234 unique gene identifiers, retaining all 23,370 bulk-tested genes and genetics-only genes.
- [Signal-specific hypothesis records](../../../RNA-seq/results/histology_anchored_continuum/translation/hac-translation-20260906T230500Z/integrated-20260906T231500Z/therapeutic_hypothesis_records_complete.tsv.gz): molecular object, perturbation hypothesis, context links, disease profile, candidate readout and unresolved experiment.
- [All 117 programs](../../../RNA-seq/results/histology_anchored_continuum/translation/hac-translation-20260906T230500Z/integrated-20260906T231500Z/all117_integrated_program_evidence.tsv), with the frozen nine windows, separate cohort effects, composition sensitivity and native modality tables in the result root.
- [Methods](../../manuscript/METHODS.md#continuum-extension-intervention-hypotheses-and-measurable-disease-states), stating the inferential unit, covariates, multiple-testing family, effect units and scope of each connection.

## What the implementation establishes

The continuum is the participant ordering or continuous score. Sliding windows
are overlapping summaries along that score. Here the two frozen scores support
participant-level regressions, while the nine windows remain descriptive. A
window is neither a biological replicate nor a disease transition.

The extension separates inherited hypotheses from disease readouts. It does
not turn modality count into evidence strength or make the genetic association
change across windows. It retains the frozen 117 program memberships and
signature-excluded scores. No model benchmark, classifier, new patient subtype,
coherence claim, clinical cutoff or separate follow-up study outcome enters the extension.

| Analysis | Verified output | Interpretation |
|---|---|---|
| Genetic replay | 5,362 pairs; 22,356,697 variant records | All exported pairs remain, including low-support signals |
| Supported main-trait signals | 925 pairs, of which 83 have resolved marginal direction | 842 remain direction-unresolved; these are signal pairs, not unique genes |
| Inherited memberships | 462 main SuSiE genes; 822 main ABF genes | Existing memberships reproduced; not redefined by expression |
| Stage/sex program support | 93/117 | Existing two-axis BH/direction rule, reproduced without changing the axis |
| Composition sensitivity | 68/117 | Same two-axis rule after adding estimated lineage log-ratios |
| Stronger separate-cohort criterion | 26/117 | Both axes pass BH separately in each evaluation cohort, with consistent direction |
| Participant denominators | 109 + 361 = 470 | Complete evaluation participants per testable program; source-overlap cohorts descriptive |
| Matched protein participants | 41, including 11 controls and 30 with MASLD | Exact initial-sample metadata links; no inferred pairing by disease stage |
| Matched protein family | 345, including six untestable genes | Whole eligible family retained in BH correction |
| Adjusted liver/plasma association | 49/345 | Participant-linked association after fibrosis and other covariates; not source-organ attribution |

The 93, 68 and 26 program counts represent different inferential criteria.
They must not be described as three independent sets of confirming modalities
or as a formal test that the coefficient changed with composition adjustment.
HC3 sensitivity results are available separately from the primary ordinary
least-squares results.

## Genetic direction and changed experimental recommendations

The signed source files use GWAS allele1 as the effect allele and eQTL EA as
the effect allele. The code aligns the complete allele pair at the same GRCh37
position and negates the eQTL beta for a swap. It leaves palindromes,
non-SNVs, allele mismatches, duplicate or missing effects unresolved. The
original posterior mass is retained, including unresolved mass. Direction is
called only when at least 95% of that original mass agrees. This conservative
criterion does not change the adopted colocalization probabilities.

The source effects are marginal. Multi-signal conditional effects are not
available in the replay exports, so a locus with multiple exported signals
does not acquire a signal-specific perturbation direction. The complete
signed variant records and direction-reason fields make the distinction
inspectable. A single exported pair can still conceal undetected signals;
resolved marginal direction is not a causal expression or pharmacological
effect. The replay posterior vectors also retain the existing documented
replay-tolerance limitation: there is no adopted per-variant posterior vector
against which to verify them independently.

Among the direction-resolved ordinary remodeling records, direct MASLD has
two same-direction signal pairs and no opposing pair; liver fat has one of
each; liver enzymes have 28 matching and 13 opposing pairs. The remaining
resolved pairs include 29 with no supported bulk remodeling, seven outside
the bulk observable universe and two axis-component records requiring separate
leave-one-out interpretation. Gene counts overlap across signals and traits.

The liver-fat pair of opposing classifications belongs to **SPATA31H1 in two
different PDFF studies**. Their individual directions are retained and a common
directional recommendation is suspended. This is not the replicated opposing
direct-MASLD example sought in the plan. The implementation does not select a
threshold or gene to manufacture that result. For each same-gene, same-trait,
same-ancestry group, the table records cross-study sign consistency and warns
that participant overlap has not been adjudicated. Distinct signals could
represent distinct mechanisms and are not pooled into one causal effect.

The practical change from expression-only interpretation is therefore
concrete but bounded: keep inherited candidates with no bulk differential
remodeling, test opposing relationships bidirectionally, and resolve inconsistent
genetic directions before choosing a common directional perturbation.
Liver-enzyme associations remain trait-specific and cannot become direct MASLD
risk evidence merely because their genes appear in liver tissue.

## What each additional modality contributes

The complete program-to-gene table connects fixed memberships to assay-native
evidence without transferring the statistical support of a program to every
member gene. The figure displays all programs in fixed lineage/module order,
not order chosen by a new association P value. It keeps “not tested,”
“untestable,” and “tested without support” distinguishable.

Single-cell source-lineage labels provide discovery context. They do not prove
transcript origin after ambient RNA correction. In particular, the 11 native
non-integer T-cell programs retain the known ambient-transport limitation.
Bulk composition sensitivity does not repair that limitation. Signal-specific
ATAC context and program ATAC tests retain separate lineage and donor fields.
Spatial coverage is complete; only the originally tested programs inherit
formal matched-background results. The all-program spatial table makes absent
tests explicit rather than treating them as negative results.

Existing paired-biopsy results retain the person as the inferential unit and
their original complete-family correction. The extension does not reinterpret
between-biopsy biological changes as technical measurement noise. Legacy
program-protein tests are supplied in their original subset and identifiers;
they are not assumed to be all-117 tests. The complete program-member readout
table instead records every member's tissue-protein observability and matched
plasma evidence. A different member protein can be a proposed readout for an
intervention experiment, but its regulation by the target is explicitly
`not_tested_shared_program_membership_only`.

The additional matched protein analysis was possible because the source
metadata explicitly link 41 unique participants' initial liver and plasma
files. It uses native log2 protein intensities, assay-separate normalization,
no missing-value imputation, and a complete 345-gene family. The association
model adjusts for categorical fibrosis, age, BMI, sex and both acquisition
batches, with at least 20 complete participants and 10 residual degrees of
freedom per test. A separate sensitivity omits fibrosis. Sample-level
denominators and absent assay links are exported. This public cohort analysis
is a measurement connection, not a new patient cohort, independent clinical
validation, secretion assay or proof that blood protein comes from the liver.

As an example of why the readout should not be assumed, **IGFBP7 does not show
adjusted matched liver/plasma support** in this analysis: beta −0.00276 SD/SD,
95% CI −0.493 to 0.487, nominal P=0.991, BH q=1, n=41. Its tissue/program
evidence therefore cannot establish that its blood level tracks the same
participants' liver protein. This is an example read from the complete family,
not a selected confirmatory hypothesis. All matched protein samples remain
unpaired with the bulk RNA cohorts; no bulk percentile is assigned to them.

## Verification and commands

All inferential computation and tests ran on Slurm compute nodes. No new
production dependency was installed. The successful evidence job used 4 CPUs,
64 GB requested RAM and `--time=90:00:00`; it completed in 4 min 23 s with
6,530,800 KB maximum recorded RSS. The matched-protein and integration jobs
used 1 CPU and 8 GB requested RAM. The standalone matched-protein analysis
completed in 17 s. Each output records `sessionInfo()` and seed 20260906.

The test script checks direct/swapped allele invariance, palindromic and
missing-allele handling, unresolved posterior mass, sign conflicts, complete
BH families, trait classification and composition-regression coefficients.
The full-data verification independently recalculates signed posterior masses
from **all 22,356,697 exported variant rows**, checks all direct/swapped effect
identities, and tests that multi-signal records remain unresolved. It also
checks complete testing families, unique gene identifiers, donor denominators,
paired-protein sample bounds and unchanged inherited memberships.

The stage/sex fixed-effect meta-analysis coefficients reproduce the frozen
source with maximum absolute difference **5.55×10⁻¹⁶**. The complete main-trait
gene counts reproduce **462/822**. Exact final counts are in
`integrated-20260906T231500Z/verified_numbers.tsv`; the verification result is
`VERIFICATION.txt`. Slurm logs are adjacent to this report: evidence job
`21478023`, matched-protein job `21478029`, final integration/verification job
`21478179`. The final scope check is recorded in `scope-check-final-v2.log`.

The final PDF was checked with `pdftotext -layout` and visually inspected after
rendering a PNG preview. The printed direction counts match the source table,
all 117 row labels are visible, and missing-test categories remain distinct.
The manuscript scope check initially rejected the word "frozen" in the new
Methods prose; replacing it with "fixed" preserves the scientific meaning.

To reproduce into a new directory from the repository root:

```bash
RUN_STAMP=$(date -u +%Y%m%dT%H%M%SZ)
sbatch scripts/analysis/histology_anchored_continuum/translation/run.sbatch \
  "RNA-seq/results/histology_anchored_continuum/translation/hac-translation-${RUN_STAMP}"
```

The submitted script runs `test_translation.R`, `01_build.R`,
`02_paired_protein.R`, `03_integrate_render.R`, `04_verify.R`, and
`scripts/manuscript/validate_resource_scope.py` sequentially. The combined
submission wrapper was updated after those constituent analyses completed;
the constituents were run and verified, while the updated all-in-one wrapper
was not submitted again. Existing output directories cause the analysis to
stop rather than overwrite results. The original failed attempts and
intermediate figures remain separate from the completed directory.

## Modified files and purpose

All paths below are relative to
`/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/`.

| File | Purpose |
|---|---|
| `scripts/analysis/histology_anchored_continuum/translation/lib_translation.R` | Allele, posterior-direction, trait, BH, cohort regression and meta-analysis rules |
| `scripts/analysis/histology_anchored_continuum/translation/test_translation.R` | Focused statistical and harmonization checks |
| `scripts/analysis/histology_anchored_continuum/translation/01_build.R` | Complete inherited evidence, composition sensitivity, Gene Catalog and native assay tables |
| `scripts/analysis/histology_anchored_continuum/translation/02_paired_protein.R` | Exact participant linkage and native matched protein analysis |
| `scripts/analysis/histology_anchored_continuum/translation/03_integrate_render.R` | Complete joined evidence records, cross-study direction qualification and vector figures |
| `scripts/analysis/histology_anchored_continuum/translation/04_verify.R` | Independent all-variant and complete-family verification |
| `scripts/analysis/histology_anchored_continuum/translation/inspect_inputs.R` | Diagnostic resolving blank genetic identifiers and checking paired assay availability |
| `scripts/analysis/histology_anchored_continuum/translation/run.sbatch` | Compute-node reproduction of the full extension |
| `docs/manuscript/METHODS.md` | Assay-specific methods and interpretation limits |
| `docs/README.md` | Navigation to the completed extension |
| `docs/PAPER.md` | Approved translational question and claim boundary |
| `docs/STATUS.md` | Implemented candidate status |
| `docs/RESULTS.md` | Exact verified candidate numbers |
| `docs/ROADMAP.md` | Concrete remaining manuscript decision and unresolved experiments |
| `docs/technical/continuum_translation_20260906T135538Z/IMPLEMENTATION.md` | This implementation review and reproducible commands |

No existing results, adopted main figures, model analyses, or genetic
memberships were overwritten. Final manuscript adoption remains a separate
decision under the repository's existing release instructions.
