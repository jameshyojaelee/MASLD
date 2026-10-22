# Manuscript workspace

This directory contains the current drafting surfaces for the standalone MASLD
Resource paper. The drafts are editable; they are not independent sources of
scientific scope or numerical truth.

The October 14 drafts and legends use the selected 774-unit bulk fit,
verified COLOC r2 and r8 Catalog. Historical 844-record estimates are labeled;
GSE193066 remains excluded from primary repeat-biopsy agreement claims.

Before editing, read:

1. [../PAPER.md](../PAPER.md) for the paper’s identity, novelty, figure logic,
   and claim boundaries.
2. [../STATUS.md](../STATUS.md) for open dependencies.
3. [../RESULTS.md](../RESULTS.md) for every number used in prose or legends.

## Draft files

| File | Role | Current constraint |
|---|---|---|
| [draft/title.md](draft/title.md) | Binding title | Records liver-trait genetic evidence and disease-state liver expression side by side for each gene; no coordination clause; no proven molecular mechanism (new wording awaits user approval) |
| [draft/abstract.md](draft/abstract.md) | Abstract | Current counts and corrected Figure 2 intervals; retain biological-unit and class-comparison limitations |
| [draft/01_intro.md](draft/01_intro.md) | Introduction | Lead with inherited variation, multicellular remodeling, and their cell/tissue context |
| [draft/fig1.md](draft/fig1.md) | Figure 1 Results | Preserve native stable-ID and Catalog-adjudicated universe distinctions |
| [draft/fig2.md](draft/fig2.md) | Figure 2 Results | Corrected coordinates and r2 values; exact anchors do not establish causal variants |
| [draft/fig3.md](draft/fig3.md) | Figure 3 Results | All five cohort contrasts and current 54-gene 3F; inference uses genetic-identity clusters, not verified donors |
| [draft/fig4.md](draft/fig4.md) | Figure 4 Results | H59 60-donor per-dataset-variance fit selects zero of 117 programs; the 62-donor Ductular injury call is a sensitivity; fixed weights and unresolved lineage origin |
| [draft/fig5.md](draft/fig5.md) | Figure 5 Results | Current 774-unit RNA effects and r2 loci; protein examples remain selection-conditioned |
| [draft/fig5_discussion.md](draft/fig5_discussion.md) | Figure 6 Results and Discussion | Selected r8 entries and failed prespecified class-direction prediction; final publication checklist remains open |
| [draft/fig7.md](draft/fig7.md) | Figure 7 Results | Broad cell-class mapping and five-clean/two-seen TranscriptFormer studies; conditional clean-study comparison; identity-unqualified ICC and ceilings remain historical diagnostics |
| [METHODS.md](METHODS.md) | Working Methods | Clearly separate prespecified methods from release-dependent values |
| [05_figure_legends.md](05_figure_legends.md) | Figure legend requirements | Populate exact values only from read-only figure source tables |

The former `manuscript/working/` tree has been retired. Its original contents
remain checksummed in the documentation archive; all current drafting happens
in `draft/`.

Historical manuscript and figure reviews are indexed under
`docs/archive/reviews/`. The remaining `manuscript/reviews/` tree contains only
the tracked checkpoint assets from the 2026-07-09 impact campaign and a short
provenance note; it is not a current writing surface.

## Drafting rules

- Use “cross-sectional stage-associated” or “stage-ordered” for unrelated
  biopsies; do not imply within-person temporal change.
- Use “colocalized” or “genetically anchored candidate”; COLOC is not proof of
  a causal gene.
- State the biological unit for every analysis.
- Preserve assay-native effect units instead of combining modalities into one
  universal score.
- Distinguish supported, discordant, tested-negative, untestable,
  not-applicable, source-dependent, and indeterminate evidence.
- Present k4/k6 NMF as continuous axes, not discrete patient subtypes.
- Cite Kamzolas et al. when introducing the molecular continuum. Keep recorded
  histology in the model, call sliding windows fixed descriptive summaries,
  and do not interpret the continuum as within-person temporal ordering.
- The future experimental follow-up paper may be motivated prospectively, but no screen roster,
  target count, guide result, phenotype, hit, or mechanism belongs here.

## Release state

Everything in `draft/` remains a working manuscript until the final release
checklist in [../ROADMAP.md](../ROADMAP.md) passes and the user explicitly
adopts the synchronized release. The material in `manuscript/release/`
describes historical July packages and is not submission readiness for the
current Resource.
