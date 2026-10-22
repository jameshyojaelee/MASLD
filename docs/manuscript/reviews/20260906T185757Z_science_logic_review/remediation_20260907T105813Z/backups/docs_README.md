# MASLD Resource paper: documentation guide

This repository supports a standalone, human-centered MASLD Resource and
analytical-framework paper. Start here; do not reconstruct the project from old
plans or dated progress logs.

## Read these five documents

| Document | Answers | Read when |
|---|---|---|
| [PAPER.md](PAPER.md) | What is the paper, what is novel, and what may we claim? | Deciding analyses, figures, framing, or journal strategy |
| [STATUS.md](STATUS.md) | What is official, what has been validated, and what remains open? | Starting any work session |
| [RESULTS.md](RESULTS.md) | Which numbers may be used now, and from which release? | Writing text, legends, slides, or reviews |
| [ROADMAP.md](ROADMAP.md) | What should happen next, in what order, and what is out of scope? | Planning or coordinating work |
| [manuscript/README.md](manuscript/README.md) | Where are the current manuscript drafts? | Editing the paper |

The order of authority is: **PAPER → STATUS → RESULTS → ROADMAP → technical
records → archived history**. A dated progress note never overrides a current
summary or a read-only release.

## The paper in one paragraph

**Title:** *A multimodal human MASLD Resource maps inherited genetic variation
and multicellular disease remodeling.* The Resource connects trait-associated
germline variation, multicohort cross-sectional disease remodeling, and
donor-aware cell and tissue context. It distinguishes coding, regulatory, and
unresolved noncoding genetic signals; localizes remodeling across cell
programs, protein, chromatin, and spatial assays; and delivers the results in a
MASLD Gene Catalog with appropriate follow-up experiments. Assay observability,
biological-unit control, provenance, and machine-readable evidence states make
these maps interpretable, but they are supporting rigor rather than the paper's
headline. A provisional direct-NAFLD GNMT locus provides the bridge from
inherited variation to a severity-linked GNMT–MAT1A–CYP2C19 RNA–stage–protein
state, without treating the three genes as one inherited mechanism. The
ECM/IGFBP7 and ductular-injury/BICC1 programs then separate disease association,
matched-null tissue organization, and lineage attribution. The aligned
117-program ambient candidate attenuates ECM/IGFBP7 more strongly than BICC1,
although neither paired effect change passes family correction. Same-atlas
donor-paired contrasts descriptively localize ECM/IGFBP7 to fibroblasts and
BICC1 to cholangiocytes, while only ECM/IGFBP7 exceeds matched spatial
background. Following Kamzolas et al. (2026;
doi:10.1038/s42255-026-01543-7), a validated histology-anchored molecular
continuum candidate adds stage-adjusted donor-level resolution without
replacing recorded fibrosis stage. Its sliding windows are fixed descriptive
summaries, not inferential or temporal bins. Exact GNMT genetics remains on hold. The paper must be scientifically useful without Cas13 results; the
future Cas13 screen prospectively tests Resource nominations but cannot
redefine its classes, figures, or claims.

The compact aging-inspired Figure 4–6 rebuild is a validated, visually reviewed
candidate under `figures/candidates/aging-inspired-resource-rebuild-candidate-2026-08-15-v8/`.
It has not been adopted into `figures/main`.

## Folder map

The implemented continuum translation extension is documented in the
[technical review](technical/continuum_translation_20260906T135538Z/IMPLEMENTATION.md)
and [Methods](manuscript/METHODS.md#continuum-extension-intervention-hypotheses-and-measurable-disease-states).
It connects inherited direction, within-stage program variation, assay context
and readouts as a candidate analysis; it does not replace the adopted evidence
classes or establish treatment-response prediction.

| Location | Contents | Authority |
|---|---|---|
| `docs/manuscript/` | Draft text, legends, and release notes | Writing surface; numbers come from RESULTS |
| `docs/technical/` | Dataset, pipeline, release, number-history, and figure indexes | Technical reference |
| `docs/literature/` | Competitor and literature map | Scientific context, not result authority |
| `docs/reference/` | Source notes and imported references | Supporting material |
| `docs/competitor_analysis/` | Detailed paper-by-paper analyses | Supporting material |
| `docs/archive/` | Superseded plans, progress, audits, and verbatim originals | Provenance only |
| `figures/main/INDEX.md` | Current six-figure directories, panel requirements, candidate states, and what has to be true before a panel is adopted | Figure-routing authority |

## Plain-language status labels

| Label | Meaning |
|---|---|
| **Current official release** | Read-only release currently allowed to support exact manuscript claims |
| **Validated update awaiting release** | Checked analysis intended for the paper, but downstream numbers and output files have not yet been adopted together |
| **Continuum working figures adopted** | The Kamzolas-adapted co-primary scores and complete molecular families met their candidate criteria; hash-pinned Figure 4F–4G and S3/S4 panels are in the working figure authorities, while Figure 3 and the synchronized release remain unchanged |
| **Validated spatial candidate** | Spatial registry, coverage, effects, portal tables, release-linked Figure 1/Figure 6 panels, and bounded diagnostics met their criteria; Figure 5A–F and its vector assembly are adopted working outputs, while the synchronized Resource release remains pending |
| **In progress** | Analysis or compute is still running; memberships and exact counts can change |
| **Supplementary boundary** | Valid stress test that did not support a main-paper upgrade |
| **Could not be tested** | Public data lacked an authoritative biological-unit or assay join |
| **Historical** | Preserved rationale or provenance; not an active instruction |
| **Aspect-specific** | An observation attributed to metabolism, inflammation, or fibrosis individually. Not available in the current release |
| **Aspect indeterminate** | The design resolves a summary severity score only, so which aspect an observation belongs to is unresolved rather than negative |
| **Variant combination** | Two or more variants at a locus considered jointly. Distinct from a single colocalized variant, and untested by the colocalization design |

Internal labels such as `F_five`, `BG-001`, candidate release IDs, job IDs, and
hashes remain in technical records and the archive. Human-facing documents use
“updated five-cohort bulk analysis,” “bulk-counting decision,” and “validated
update awaiting release.”

**A NAS sum is not aspect resolution.** NAS is the sum of steatosis, ballooning,
and lobular inflammation. A result that tracks the sum cannot say which of the
three it belongs to. Say "severity-associated" and record the aspect as
indeterminate; do not say "inflammation-associated" without a component-level
analysis.

## Common tasks

| If you need to… | Use |
|---|---|
| cite a sample, gene, DEG, program, or overlap count | [RESULTS.md](RESULTS.md) |
| understand why two numbers differ | [technical/NUMBERS_HISTORY.md](technical/NUMBERS_HISTORY.md) |
| check a live dependency or decide what can run now | [STATUS.md](STATUS.md) and [ROADMAP.md](ROADMAP.md) |
| edit the story or a figure section | [PAPER.md](PAPER.md), then [manuscript/README.md](manuscript/README.md) |
| locate a current main panel or generator | [`figures/main/INDEX.md`](../figures/main/INDEX.md), then the figure-specific `CANDIDATE_PANEL_INDEX.tsv` |
| inspect pipeline or dataset ownership | [technical/DATASETS_AND_PIPELINES.md](technical/DATASETS_AND_PIPELINES.md) |
| inspect the noncoding upgrade | [STATUS.md](STATUS.md) for the current checks, [RESULTS.md](RESULTS.md) for validated candidate counts, `scripts/manuscript/noncoding_resource/` for the producers, and `Analysis/Multimodal_Program_Projection/candidates/gene-catalog-v2-contract-2026-08-11/` for the human-readable Catalog rules |
| inspect or hand off the validated spatial candidate | Start with the concise [spatial omics handoff](technical/SPATIAL_OMICS_HANDOFF_2026-08-12.md), then inspect `Analysis/Multimodal_Program_Projection/candidates/spatial-resource-candidate-2026-08-11/`, `spatial-diagnostic-bundle-candidate-2026-08-11/`, and the result-led `spatial-impact-figures-candidate-2026-08-12/`; canonical Figure 5F now pairs the accepted matched-gene calibration with one outcome-blind registered Guilliams et al. H&E overlay, and accepted results and layout decisions are summarized in [RESULTS.md](RESULTS.md) |
| inspect the validated single-cell ATAC candidate | Use `Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/`; `NON_GENETIC_READY`, `GENETIC_READY`, and `FULL_READY` mark the non-genetic and adopted-COLOC replay outputs as complete. Canonical Figure 5D uses the completed primary-pair context without adopting the entire ATAC candidate package |
| inspect the histology-anchored molecular continuum | Start with `scripts/analysis/histology_anchored_continuum/README.md`. `hac-molecular-layers-20260818T173348Z` retains transcript/pathway/program inference; `hac-system-inference-20260824T232135Z` retains the participant-level 43-system models and permutations. The corrected atlas and pooled candidates retain Figure 4G/S4Q–S4U values. The active display revision uses `histology-continuum-molecular-systems-caption-free-20260825T130145Z` and `histology-continuum-molecular-systems-s4u-caption-free-20260825T130145Z`: embedded explanatory prose is removed from Figure 4G/S4Q–S4U and retained in the external caption authority. S4U keeps all 43 systems, nine descriptive windows, and 18 strict maxT asterisks; the late-minus-early contrast remains the retained display-only row-order statistic but is not displayed, and the 33 primary BH calls remain in the source registry. Mirrored source tables, hashes, and audits are in the Figure 4 authority. |
| audit a release or recover an old plan | [technical/RELEASE_PROVENANCE.md](technical/RELEASE_PROVENANCE.md) and [archive/INDEX.md](archive/INDEX.md) |
| compare recent competitors | [literature/REVIEW.md](literature/REVIEW.md) |

## Safety rules

- Never mix numbers from the current official release and the validated update
  in one claim without labeling both.
- Describe stage contrasts only as cross-sectional stage-associated remodeling.
- Treat the molecular continuum as complementary to recorded histology. Fixed
  sliding windows are descriptive only and cannot select features, define
  significance, orient an axis, or imply within-person temporal ordering.
- Cross-cohort plots use Healthy, Steatosis, Steatohepatitis, and Cirrhosis;
  retain MASL/MASH/NASH only for source-specific provenance or contrasts, and
  never treat generic MASLD as a stage.
- COLOC identifies shared-signal candidates; it does not prove a causal gene.
- Cells, spots, arrays, sections, guides, or sequencing runs are not biological
  donors.
- In spatial results, zero or missing coverage is not negative evidence;
  distinguish `untestable`, `indeterminate`, and `source_dependent` states.
- Do not use future Cas13 outcomes to select Resource genes, classes, figures,
  thresholds, or claims.
- Do not reopen broad public-data mechanism hunting as a submission dependency.

## What happened to the older documents?

The former current-state, scope, narrative, outline, number ledger, plan pack,
and progress logs were retired from the active documentation tree on
2026-08-11. Their verbatim, checksummed originals remain recoverable through
[archive/INDEX.md](archive/INDEX.md) and the migration manifest. Historical
files are provenance, not current instructions.
