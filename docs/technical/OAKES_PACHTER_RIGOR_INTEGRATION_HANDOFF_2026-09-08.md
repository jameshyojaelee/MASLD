# Integrate lessons from Oakes and Pachter into the MASLD Resource

## Purpose and handoff

User-approved execution plan, saved 2026-09-08. The five existing paper
authorities continue to govern scope and results. Execution findings belong in
the accompanying assessment, not in this plan.

Prepare a focused integration of statistical safeguards, transcriptomic QC, and
model reporting motivated by [Oakes and Pachter, *Transcriptomic responses to
endurance exercise training in rats*](https://link.springer.com/article/10.1186/s12863-026-01440-8).
Use the [analysis repository](https://github.com/pachterlab/OP_2025) and
[public peer review](https://link.springer.com/article/10.1186/s12863-026-01440-8/peer-review)
to distinguish demonstrated methods from exploratory interpretations.

The user selected **focused rigor integration**: inspect existing evidence,
address bounded QC and statistical gaps, and improve Methods and supplementary
tables. New isoform discovery, model-development studies, and wholesale
reprocessing are excluded.

## Starting evidence and operating boundaries

Read the five project authorities and applicable AGENTS.md instructions at
execution time. Record the working-tree baseline and preserve existing changes.
The repository contains extensive unrelated modifications.

Inspection established that:

- Bulk QC code already implements expression-based sex checks, within-cohort
  PCA outlier detection, and library-size checks. Its fragment-count workflow
  can preserve earlier QC decisions while recording updated count totals.
  Determine which mode actually generated the selected data before describing
  QC as recomputed.
- Single-cell Methods already describe donor pseudobulk, collapse of runs and
  sort fractions, and ambient-RNA sensitivity.
- Model cards and evaluation records already exist in the benchmark directory.
  The relevant severity card is
  `Analysis/MASLD_Model_Benchmark/release/masld-severity-v1.2/MODEL_CARD.md`.
- A genetic sample-identity check was not established by this inspection. A
  script named "patient fingerprint" concerns predictive phenotypes and does
  not establish donor identity.
- The manuscript explicitly excludes primary isoform claims. Existing
  gene-level outputs sometimes use "transcript" terminology.

These are inspection findings, not an independent reproduction of the analyses.
The earlier conversational conclusion that no major method was missing must
not become a manuscript claim.

Preserve the selected fragment-count analysis, fixed 117 programs, accepted
thresholds, and separate experimental-paper scope. Any corrected scientific outputs
remain separately identified candidates pending the existing adoption process.

## Execution work packages

### 1. Build a compact evidence assessment

Create one technical assessment alongside this handoff. For each safeguard,
record its scientific purpose, producing analysis, available result, relevant
manuscript claim, and disposition:

- Verified and adequately reported.
- Verified but inadequately reported.
- Required check missing or not reproducible.
- Not applicable or outside this integration.

Cover sample identity, counting units, biological replication, confounding,
multiple testing, predictive validation, contamination, and transcript
resolution. Existence of a script alone is insufficient to mark a safeguard
verified.

### 2. Reconcile sample identity, QC, and analysis denominators

Use the selected sample manifests, source metadata, QC outputs, and existing
donor crosswalks.

- Check unique cohort/sample keys, count-to-metadata alignment, repeated
  libraries, repeated biopsies, and explicitly documented participant overlap.
- Distinguish reported sex, expression-inferred sex, disagreement, and missing
  annotation. Trace which field enters each model and how disagreements affect
  eligibility.
- Establish whether genotype-based identity evidence already exists. If
  unavailable, report that limitation. Expression similarity or a sex
  discrepancy can flag a record for review but cannot independently authorize
  relabeling, merging donors, or exclusion.
- Verify counting units and strand settings from the selected outputs.
  Separate fragments, reads, and UMIs, and distinguish inherited QC decisions
  from QC recalculated on corrected fragments.
- Inspect exclusions by cohort and disease group. Determine whether an
  apparent outlier reflects technical failure, unusual biology, or unresolved
  evidence; preserve the existing inclusion rule.

Produce a supplementary table with one row per cohort, assay, and analysis:
source records, count-matched records, technical-QC records, unique participants,
final analyzed participants, clinical groups, and exclusion reasons. Report
overlapping exclusion flags separately from sequential losses so totals remain
interpretable. Do not sum assay-specific participant counts into a unique-person
total without a verified crosswalk.

Acceptance: every reported denominator reproduces from included identifiers,
and unresolved identity questions remain explicit.

### 3. Verify statistical safeguards supporting current claims

Review actual models and saved results supporting primary bulk, stage,
single-cell, and fixed-program claims.

- Confirm biological replication, model terms, reference groups, design
  estimability, assay units, effect estimates, uncertainty, and
  missing-covariate handling.
- Check that current stage results use the documented handling of healthy
  versus disease-only F0. Do not revive superseded early-stage claims.
- Verify normalization and mean–variance modeling for the relevant assay.
  Preserve limma-voom for the selected bulk analysis.
- Recompute BH adjustment from complete saved test families where possible.
  Explicitly inspect selection of the "best" cell type, contrast, or assay:
  within-family adjusted P values do not automatically correct selection
  across families.
- Preserve the distinction between a point-null test followed by a fold-change
  filter and the effect-size-aware TREAT comparator. Describe each according
  to its actual statistical guarantee.
- Reuse existing cohort-omission, composition, ambient-RNA, and
  robust-standard-error sensitivities. Report failed and untestable branches
  alongside successful ones.

If a demonstrated defect changes a claim, identify the affected gene
membership, effect, uncertainty, or evidence state. Correct the smallest
affected analysis using its complete original testing family. Do not change
thresholds to retain a preferred result. Where inputs or compute are
unavailable, qualify the claim and report the exact unfinished check.

Acceptance: each primary statistical statement has a traceable biological unit,
estimand, adjustment set, and multiplicity family.

### 4. Integrate predictive-model validation into the manuscript

Restrict work to models and tasks currently supporting Figure 7. Use existing
model cards, split assignments, prediction files, and exposure records.

For each reported comparison document:

- Training, tuning, and evaluation participants, including repeat-biopsy
  grouping and participant overlap.
- Whether preprocessing, feature selection, latent representations, and
  hyperparameters were fitted using evaluation samples or outcomes.
- The distinction between an external cohort, a development-exposed cohort,
  and a cohort encountered during pretraining.
- The task-native baseline, identical evaluation participants, primary metric,
  and uncertainty for the paired performance difference.
- Input coverage, refusals, and any use of the evaluation cohort to set
  normalization or acceptance thresholds.

Reproduce reported metrics from stored participant-level predictions where
available. If uncertainty is missing and predictions permit it, use 2,000
paired participant bootstrap draws with seed 20260908, keeping repeated
measurements together. Label these intervals as conditional on fitted models;
they do not capture uncertainty from retraining. Preserve existing valid
uncertainty analyses rather than replacing them.

Keep metrics aligned with the estimand. Rank correlation supports ordering;
agreement and calibration require measurements on a meaningful common scale.
Repeated-biopsy ICC also requires qualification for possible biological change
and shared technical effects.

Correct the earlier normalization explanation: linear regression can use log
expression. Its validity does not require unlogged counts; the transformation
changes the relationship being modeled. A particular scVI result cannot
establish general superiority of linear models.

Acceptance: every retained model-performance claim has a documented evaluation
population, appropriate baseline, exposure qualification, and supported
interpretation.

### 5. Incorporate transcript-resolution and contamination lessons

Add a concise limitations paragraph explaining that gene totals can conceal
changes among isoforms. Distinguish differential gene expression, differential
transcript expression, and differential transcript usage. Call existing
gene-count analyses gene-level wherever "transcript" could imply isoform
resolution.

Describe short-read assignment uncertainty and annotation dependence without
implying isoform switching has been tested. Do not transfer the rat study's UMI
processing or transcript P-value aggregation into the human pipeline.

Use existing ambient-RNA and technical-QC findings to qualify lineage and
disease interpretations.

## Deliverables, verification, and completion

Update `docs/manuscript/METHODS.md` and the active limitations text identified
through manuscript navigation. Deliver the technical assessment, reproducible
supplementary denominator table, and model-evaluation summary. Modify the five
paper authorities only where a verified finding changes current status,
numbers, or interpretation.

Use existing environments. Run analytical checks on an already authorized
compute allocation; otherwise prepare exact commands and report them as
**not run**. This handoff does not authorize new dependencies, raw-read
processing, GPU training, job arrays, or publication.

Verification must establish:

- Identifier joins do not inflate biological sample size.
- QC losses and final analysis counts reconcile.
- Saved effects, metrics, and multiple-testing corrections reproduce within
  recorded precision.
- Model and baseline comparisons use matched evaluation participants.
- Updated text accurately distinguishes completed checks from proposed work.
- Any changed figure agrees with its source table, passes text extraction,
  and is visually inspected.
- `python3 scripts/manuscript/validate_resource_scope.py` passes after
  authority or manuscript changes.

Completion means the handoff and assessment exist, manuscript changes and tables
are reviewable, and every unresolved item identifies the affected claim and
missing evidence. A new positive biological finding is not required.
