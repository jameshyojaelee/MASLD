# Number history and supersession guide

Use [../RESULTS.md](../RESULTS.md) for current writing. This document explains
why recurring values differ and points to the complete historical ledger.

## High-frequency number changes

| Topic | Older value | Newer or correct interpretation | Why it changed |
|---|---:|---|---|
| Pooled disease-vs-control samples | 846 | 844 in the validated fragment update | Two samples fail corrected fragment-based QC |
| Pooled genes tested | 27,638 | 23,370 in the validated fragment update | Corrected counting and filtering universe |
| TREAT DEGs | 1,918 | 1,616 in the validated fragment update | Fragment rather than read counting; synchronized adoption pending |
| Complete human QC pass | 1,260 | 1,257 in corrected fragment source | Three additional fragment-QC failures |
| Pre-recount human counts | 1,277 / 1,259 | 1,281 count-matched; release-specific QC totals above | Count registry audit |
| Andrews fibrosis sample size | 58 “donors” | 58 runs representing 46 donors | Run-to-donor collapse |
| Disease-significant Hotspot modules | 54 | Historical donor-level reanalysis: 22; present frozen selection: four primary/stable, two robust display | Run-level pseudoreplication and later registry criteria |
| Stage ligand–receptor hits | 13 | 0 of 13 pass donor-level family FDR | Donor-level correction |
| Sex-dimorphic genes | 290 | 8 under fixed-cohort model | Retired random-slope specification |
| NMF k=6 reproducibility | 92.8–99.3% phrasing | Matched factor cosine 0.976; hard sample partition ARI about 0.465 | Subspace stability had been conflated with discrete clustering |
| Genetic/state non-overlap | 92.4% | Current official release only; corrected value pending | Active COLOC membership rerun |

## Threshold authorities

| Analysis | Current rule |
|---|---|
| Human Tier-1 bulk DE | limma-voom quality-weighted, fixed cohort term, TREAT `FDR < 0.05` at interval null `|logFC| ≤ 0.25` |
| Adjacent cross-sectional stages | Contrast-native filtering and normalization; final model not yet run on fragment counts |
| Primary COLOC | SuSiE PP.H4 threshold under the frozen trait/provenance registry, with ABF retained as a distinct layer |
| Hotspot program selection | Complete 117-program family with donor-level association and frozen stability rules |
| Program testability | At least 8 mapped genes and 20% retained original L1 weight for external projection |

## Denominator checklist

Before using a percentage, record:

1. named genes or Ensembl IDs;
2. the release and tested universe;
3. joint assay testability;
4. whether the denominator is genes, loci, donors, cohorts, samples, or programs;
5. whether proxy-trait and ancestry strata are included; and
6. whether missing values mean negative evidence or lack of observability.

The genetic/state overlap percentage is especially sensitive to these choices.
Never divide by the full atlas when the question requires jointly testable
genetic genes.

## Full historical ledger

The previous `docs/manuscript/NUMBERS.md` is preserved byte-for-byte at:

`docs/archive/documentation_consolidation_2026-08-11/originals/docs/manuscript/NUMBERS.md`

It contains dated corrections, retired values, detailed pipeline outputs, and
source paths. It is an audit trail, not the default source for new writing.
