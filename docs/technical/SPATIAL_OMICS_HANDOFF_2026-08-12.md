# Spatial omics handoff — 2026-08-12

Supersedes the 2026-08-11 handoff, which was renamed to this file rather than
forked. The preferred candidate moved to
`spatial-impact-figures-candidate-2026-08-12` (presentation and claim rules
only), and the HMSMA class arm was corrected. See
*What 2026-08-12 changed against 2026-08-11* and *HMSMA class arm*.

This is a concise technical handoff, not a paper authority or release decision.
Current scientific scope, status, exact results, and next dependencies remain in
[`PAPER.md`](../PAPER.md), [`STATUS.md`](../STATUS.md),
[`RESULTS.md`](../RESULTS.md), and [`ROADMAP.md`](../ROADMAP.md).

## Headline finding

Positive spatial autocorrelation is not sufficient evidence of
program-specific tissue architecture. Both predeclared robust disease programs
had positive residual Moran statistics in both Visium sources, with positive
unit-level direction in 4/4 and 4/4 GSE192741 donors and 10/10 and 9/10 Vu
arrays. Matched expression-, detection-, biotype-, mitochondrial/ribosomal-,
and lineage-correlation-aware gene nulls separated them:

- ECM/IGFBP7 exceeded expected spatial background in GSE192741 and within Vu;
  the Vu result remains source-dependent because arrays cannot be joined to
  donors.
- Ductular injury/BICC1 overlapped expected spatial background and remains
  indeterminate, not negative.
- Full adjustment for all 16 shared cell-abundance factors attenuated IGFBP7
  but preserved the asymmetric IGFBP7-supported/BICC1-indeterminate result.

The paper contribution is therefore matched-null calibration of spatial
transportability, not a new causal mechanism or hero-gene claim.

## What was completed

- Registered every spatial source with biological-unit, join, observability,
  source-dependence, and permitted-claim criteria.
- Exported outcome-free coverage for all 117 frozen programs across seven
  assays: 819 program-by-assay rows.
- Restricted confirmatory inference to the two programs frozen before spatial
  outcomes were inspected.
- Reused the previously fixed 9,999-draw matched-gene nulls and donor-first GSE192741
  inference; Vu remains array-level and source-dependent.
- Completed label-blind HMSMA QC and coverage for 35 arrays and 130,097 retained
  spots without treating either as biological sample size.
- Added categorical spatial products for the portal and Gene Catalog; legacy
  numeric spatial scores are compatibility-only and excluded from rankings.
- Completed bounded composition, physical-unit, null, weight, cell-context,
  and HMSMA mask diagnostics without expanding the confirmatory family.
- Simplified the spatial publication family and replaced Figure 4A/4F with a
  result-led calibration sequence. The complete source-state matrix moved to
  the supplement.

## Preferred candidate

The preferred, not-yet-adopted assembly is:

`Analysis/Multimodal_Program_Projection/candidates/spatial-impact-figures-candidate-2026-08-12/`

Key files:

- `panels/fig4a_spatial_calibration_logic.pdf` — compact inferential logic.
- `panels/fig4f_matched_null_discrimination.pdf` — observed programs versus
  matched-null 5th–95th percentiles, unit signs, and evidence states.
- `panels/figS_multimodal_program_source_matrix.pdf` — complete supplementary
  assay/source state matrix.
- `data/spatial_transportability_profile.tsv` — exact four-row source table.
- `impact_claim_contract.tsv` — allowed and prohibited interpretation.
- `READY`, `panel_manifest.tsv`, and `source_manifest.tsv` — validation and
  provenance records.

The candidate contains 14 one-page PDFs, retains 12 accepted r4 panels
byte-for-byte, adds no inferential family or score, and does not write canonical
outputs.

### What 2026-08-12 changed against 2026-08-11

Presentation and claim rules only. No statistic, estimand, null, program
selection, or evidence state was recomputed, and the four-row transportability
profile is numerically identical.

- Figure 4F labelled each row `N/N donors +`. A bare plus sign beside an
  evidence state reads as calibrated support, when it is only the raw sign of
  the uncalibrated statistic. Rows now read `N/N donors I > 0` under a
  `before calibration` / `after calibration` column header pair, so the
  identical raw sign and the divergent verdicts are legible as a sequence.
- The claim rules gain a fifth prohibition: *a positive per-unit Moran sign
  is evidence of program-specific organization*. `50_validate_...py` expects
  4 allowed / 5 prohibited accordingly.

The 2026-08-11 candidate remains unchanged, valid, and byte-intact; it is
superseded as *preferred*, not invalidated.

## Validation state

For the 2026-08-12 candidate:

- Impact build and validation SLURM job `19790869`: completed, exit 0.
- Validator checks: 8/8 pass (`validation_report.tsv`), including
  `IMPACT-UNITS` rederiving the 4/4, 10/10, 4/4, 9/10 unit signs independently.
- Panel checksums: 14/14 matched.
- Source checksums: 11/11 matched.
- Resource-scope validation when the spatial candidate was finalized: passed.
- Visual review: **not yet done** for the two relabelled main panels. This is
  the open item in *Next actions*.
- Portal spatial semantic tests were not re-run: they exercise portal sources,
  which this candidate does not touch. The 35/35 result from the 2026-08-11
  candidate still stands.

For the superseded 2026-08-11 candidate: build job `19751610` and verification
job `19751830`, both exit 0, 35/35 semantic tests, 14/14 and 11/11 checksums,
visual review passed.

A current global scope-validator rerun now reports two stale stage-analysis
markers (`model not run` and `no stage DEG number is`) because `STATUS.md` and
`RESULTS.md` were independently advanced to a validated stage candidate after
the spatial candidate was finalized. Those concurrent stage changes were not
modified here; the validator requirements must be reconciled by that workstream.

The general release-consistency checker still exits 1 because the existing
July release lacks `strict_drug_verdict.json` and
`credible_set_coding_summary.tsv`. This is a pre-existing defect in the
non-spatial release requirements, documented in `STATUS.md`; it is not a failure
of the spatial candidate.

## HMSMA class arm

A separate line of work from the two-program confirmatory family: it tests
frozen *evidence classes* (`frozen_evidence_classes.tsv`), not programs, and is
in no finalized candidate. Current result is
`Analysis/Spatial/results/hmsma_class_spatial_v3`. Both
`hmsma_class_spatial/` (2026-08-10) and `hmsma_class_spatial_v2/` are
superseded and must not be cited.

Two defects were corrected, one of them material.

- **Unit of inference.** The original ran one-sample t-tests and Wilcoxon tests
  across the 35 arrays, which assumes the donor independence this dataset
  cannot establish — there is no array-to-donor key. Replaced by a
  9,999-draw matched-gene-set null with the arrays held fixed, matching the
  estimand the confirmatory arm already uses. Per-array values are now
  descriptive only, and the p-values are conditional on these arrays rather
  than generalizable to new donors.
- **Result change.** Genetic-only genes no longer exceed matched background
  (`+0.0018, p = 1.2e-05` → `-0.0001, p = 0.59`). The old positive came from a
  null whose expectation was mis-weighted by control-pool size instead of class
  composition, compounded by class genes in control-free cells entering the
  observed mean but not the expectation. Disease-state-only
  (`+0.0112, p < 1e-04`) and convergent (`+0.0075, p < 1e-04`) survive.
  The reading is now a dissociation, and the retired `7.1x gradient` framing
  must not be reused; with the genetic-only delta on zero the ratio is
  unstable and meaningless.
- **Hex geometry.** `obsm['spatial']` was raw array indices, which are
  anisotropic, so a 6-nearest-neighbour graph tied the two true `col ± 2`
  neighbours against `row ± 2` spots that are not neighbours.
  `53_build_hmsma_anndata.py` now emits micrometres. This changed effect
  magnitudes by 5–6% and no conclusion. It affected **only** this arm: the
  program arm builds its lattice from explicit integer offsets in
  `22_build_hmsma_label_blind.py::moran_hex` and never reads `obsm`, so its
  already-fixed values are unaffected and need no rebuild.

## Claim boundaries

- Do not call BICC1 spatially negative.
- Do not call IGFBP7 cell-autonomous, causal, or a proven fibrosis mechanism.
- Do not read a positive per-unit Moran sign as calibrated support; it is the
  raw sign of the uncalibrated statistic and is positive for both programs.
- Do not cite the retired HMSMA class-arm numbers — the `7.1x` gradient, or
  genetic-only genes exceeding their matched null.
- Do not describe Vu as donor-level replication.
- Do not perform HMSMA disease, NAS, fibrosis, histology, or MALDI inference
  until an authoritative array-to-donor and registration table arrives.
- Do not treat spots, sections, arrays, AOIs, or cells as independent donors.
- Do not adopt GeoMx, GSE287826, the old 22-program family, trajectory,
  communication, or imputed-chromatin outputs into main-paper evidence.
- Do not mutate the validated candidate; any changed analysis or artwork must
  use a new versioned candidate.

## Next actions

1. Scientifically adjudicate the two new Figure 4 panels and wording.
2. Assemble them with the non-spatial panels only after corrected COLOC, bulk,
   stage, composition, and Gene Catalog inputs share one synchronized release.
3. Adopt the spatial products only as part of that user-approved Resource
   release; do not adopt this candidate independently.
4. If authoritative HMSMA metadata arrives, freeze the joined analysis
   specification before reading outcomes and build a new versioned candidate.
