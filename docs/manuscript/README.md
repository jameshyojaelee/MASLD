# Manuscript workspace

**LaTeX and reference-paper formatting, 2026-11-08.** The editable main paper is
[latex/main.tex](latex/main.tex); supplementary notes, legends and artwork are
in [latex/supplement.tex](latex/supplement.tex). Both use XeLaTeX. The supplied
Modformer PDF guided the single-column layout, centered title and author block,
numbered section headings, 10-point Times-style Nimbus Roman text and 9-point
captions. Regular, bold, italic and bold italic fonts are bundled with their
license. Each main figure is embedded in its Results section as one complete,
indivisible plate on a single page, with its full caption following the artwork. The original
individual PDF panels remain editable and unchanged.

The self-contained [updated Overleaf ZIP](masld-overleaf-complete-figures-20261109.zip)
contains the current project; the original preparation ZIP is retained.
Compiled previews are [main paper](latex/build/main.pdf) and
[supplement](latex/build/supplement.pdf). The initial source selection is
recorded in `scripts/manuscript/build_latex_project.py`; subsequent typography
and placement are edited directly in `latex/`. Figure 1A/D remain proposals;
scientific qualifications and unresolved citation/deposit fields are preserved.

The Summary and Introduction were rewritten around the biological problem,
approach, interpretation and use, without result-count or statistical-metric
inventories. Recent Sanjana-authored publication/preprint abstracts indexed in
PubMed for 2024–2026 were reviewed, with available introductions used to guide
structure, including [the lncRNA paper](https://doi.org/10.1016/j.xgen.2026.101253),
[the MYC regulatory study](https://doi.org/10.1038/s41467-025-56568-5) and
[the isoform-targeting study](https://doi.org/10.1038/s41467-025-62066-5).
The Google Doc retains the preceding opening text because the user redirected
the work to LaTeX before that rewrite was applied online.

The [existing Overleaf project](https://www.overleaf.com/project/6af1056703ca6c90d76b4225)
is connected through its Git integration. The local Overleaf checkout is
`/gpfs/commons/home/jameslee/overleaf-masld-6af1056703ca6c90d76b4225`.
Its initial files matched the prepared project exactly; the reference-paper
formatting is synchronized through that separate manuscript repository.
Both local compilations, input-file existence, all panel/source SHA256
comparisons and the Resource scope check passed. No scientific rerun was
performed for this formatting update.

**Google manuscript reconciled, 2026-11-08.** The
[main manuscript and section outline](https://docs.google.com/document/d/1n1PFbzUBWkKrI8kXXLvFtHvRRNFrzQ1alja9DLi8SBI/edit?tab=t.0)
now use the current Summary, Introduction, Figure 1–7 Results and Discussion.
The author block and existing reference entries were preserved; references
from the current local bibliography were added when their DOI was absent.
Citation placeholders and unassigned panel labels remain unresolved. Methods,
legends and Supplementary Notes remain in their local drafting files; they were
not added as new Google Doc tabs. This is a working-text update, not adoption of
new genetic memberships, unique-cell effects or a submission-ready package.

This review added missing qualifications to [Introduction](draft/01_intro.md),
[overview](draft/fig1.md), [genetics](draft/fig2.md),
[Catalog Results and Discussion](draft/fig5_discussion.md),
[Methods](METHODS.md) and [legends](05_figure_legends.md), with the same boundary
recorded in [STATUS](../STATUS.md). Genetic-dependent numbers remain provisional;
historical atlas-derived estimates are not verified unique-cell results.
`python3 scripts/manuscript/validate_resource_scope.py` passed. Google Docs
readback matched the prepared main text and all 21 outline corrections exactly,
and the author block's native structure was unchanged. Scientific reruns were
not run for this text review. The fixed-program run failed and its downstream
verification was cancelled in the scheduler; no corrected effects were adopted.

**Working corrections and Catalog adoption, 2026-11-06.** Current drafts and Methods use the checked 71-donor hepatocyte gene fit (90 BH calls), selection-aware Xenium inference and adjusted exon/intron ratio results (2,711 calls, 114 highlighted genes). The adopted Catalog eligibility release, `gene-catalog-eligibility-20261106`, changes 188 genetic states and 180 classes, retaining all RNA calls/effects and positive genetic genes: 13,731 untested, 12,958 unresolved, 1,118 established-state associated, 400 genetically anchored and 26 supported on both axes. Figures 1B and 3G use the qualified 3,589-gene native comparison; Figure 6 uses 14,502 jointly eligible Catalog identities within 28,233 total. Named r8/r9 files remain historical. New LD and signal-matched comparisons remain separate research candidates. See [current values](../RESULTS.md) and the [analysis plan](../ROADMAP.md#scientific-analysis-plan-2026-11-06).

This directory contains the current drafting surfaces for the standalone MASLD
Resource paper. The drafts are editable; they are not independent sources of
scientific scope or numerical truth.

The October 14 drafts and legends use the selected 774-unit bulk fit,
verified COLOC r2 and historical r8 Catalog; the November 6 eligibility release supersedes the working Catalog classes. Historical 844-record estimates are labeled;
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
| [draft/fig4.md](draft/fig4.md) | Figure 4 Results | Source-diagnosis-corrected per-dataset-variance fit (at most 59 donors) selects zero of 117 programs (min q 0.18135); equal-variance fit is a sensitivity; fixed weights and unresolved lineage origin |
| [draft/fig5.md](draft/fig5.md) | Figure 5 Results | Current 774-unit RNA effects and r2 loci; qualified disease accessibility and S5P–AA/Note 16 added October 26; protein examples remain selection-conditioned |
| [draft/fig5_discussion.md](draft/fig5_discussion.md) | Figure 6 Results and Discussion | Adopted November 6 eligibility classes; historical r8 class comparisons remain labeled; final publication checklist remains open |
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
