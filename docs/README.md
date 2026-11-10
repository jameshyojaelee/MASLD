# MASLD Resource paper: documentation guide

Use this page to find the manuscript, scientific decisions and analyses.
The five documents below remain the scientific authorities; directory indexes
are navigation aids. Latest explicit qualifications take precedence over older
dated entries and historical tables.

## Read these five documents

| Document | Answers |
|---|---|
| [PAPER.md](PAPER.md) | Scientific identity, novelty, figure logic, permitted claims and the two-paper boundary |
| [STATUS.md](STATUS.md) | Current manuscript, selected outputs, unresolved evidence and remaining work |
| [RESULTS.md](RESULTS.md) | Exact numbers, source analyses, release states and current qualifications |
| [ROADMAP.md](ROADMAP.md) | Scientific priorities, authorized work, closed analyses and final release checklist |
| This README | Navigation, terminology and document roles |

The order of authority is **PAPER → STATUS → RESULTS → ROADMAP → technical
records → archived history**. A folder name, report file, timestamp or dated
progress note does not select a scientific result.

## The paper in one paragraph

This standalone human Resource records liver-trait genetic evidence and
cross-sectional disease-state expression side by side for each gene, adds
donor-aware cell and tissue context, and organizes the evidence in the MASLD
Gene Catalog with appropriate follow-up experiments. The fixed programs and
assay-native measurements retain their limits. Experimental follow-up outcomes
cannot redefine Resource evidence classes, thresholds, figures or claims. The
binding title, full goals and claim boundaries are in [PAPER.md](PAPER.md).

## Folder map

| Material | Entry point |
|---|---|
| Working paper, supplement and editable sources | [Manuscript workspace](manuscript/README.md) |
| Main figure panels and source selections | [Figures 1–7](../figures/main/INDEX.md) |
| Supplementary artwork and figure alternatives | [Figure directory](../figures/README.md) · [All candidate folders](../figures/candidates/README.md) |
| Analyses by modality or question | [Analysis directory](../Analysis/README.md) · [Hypothesis folders](../Analysis/Hypotheses/README.md) |
| Dataset, script and environment routing | [Dataset and pipeline index](technical/DATASETS_AND_PIPELINES.md) |
| Local Gene Catalog | [Selected Catalog](../RNA-seq/results/evidence_passports/releases/gene-catalog-eligibility-20261106/index.html) · [Current interpretation](RESULTS.md#catalog-eligibility-adoption-2026-11-06) |
| Data storage and environment specifications | [Data registry](../data/REGISTRY.md) · [Environment catalog](../envs/README.md) |
| Literature and competitors | [Literature review](literature/REVIEW.md) |
| Historical decisions and previous numbers | [Number history](technical/NUMBERS_HISTORY.md) · [Release history](technical/RELEASE_PROVENANCE.md) · [Archive](archive/INDEX.md) |

## Plain-language status labels

| Label | Meaning |
|---|---|
| **Selected / adopted** | Explicitly chosen for the stated analysis, figure or Catalog; does not adopt every related output or remove later qualifications |
| **Validated update / candidate** | Checked within its recorded scope; manuscript use still follows the explicit decision in STATUS and RESULTS |
| **Working figure** | Editable manuscript artwork; distinguish it from a named historical publication bundle |
| **In progress** | A dated analysis state; check SLURM for live compute |
| **Supplementary boundary** | A retained stress test or limitation with no upgrade to a main-paper discovery |
| **Could not be tested** | Required biological-unit identity, assay coverage or source join is unavailable |
| **Historical** | A prior result or decision; latest qualifications and selected outputs govern current use |
| **Aspect indeterminate** | A summary severity result does not resolve metabolism, inflammation and fibrosis individually |
| **Variant combination** | Joint consideration of multiple variants, which single-variant or colocalization evidence alone does not test |

## Common tasks

| Task | Read |
|---|---|
| Cite a count or effect | [RESULTS](RESULTS.md), including its latest qualifications and the relevant release state |
| Understand discrepant numbers | [Number history](technical/NUMBERS_HISTORY.md) and the dated source record |
| Edit text or a legend | [PAPER](PAPER.md), then the [manuscript workspace](manuscript/README.md) |
| Locate a current panel | [Main figure index](../figures/main/INDEX.md), then the figure-specific README and source index |
| inspect the histology-anchored molecular continuum | [Analysis overview](../scripts/analysis/histology_anchored_continuum/README.md), [current results](RESULTS.md#histology-anchored-molecular-continuum-candidate) and [Figure 4](../figures/main/fig4_singlecell_programs/README.md) |
| Locate a hypothesis or past figure alternative | [Hypothesis index](../Analysis/Hypotheses/README.md) or [candidate index](../figures/candidates/README.md); follow linked authorities for scientific state |
| Plan or resume work | [STATUS](STATUS.md) and [ROADMAP](ROADMAP.md); a dated job note is not live scheduler state |
| Inspect detailed methods or earlier decisions | [Dataset and pipeline index](technical/DATASETS_AND_PIPELINES.md) and [archive index](archive/INDEX.md) |

## Safety rules

- Use biological donors or participants as inferential units; runs, cells,
  spots and sections do not increase biological replication.
- Preserve native assay units and distinguish unsupported evidence from
  untestability, non-applicability and source dependence.
- Describe unrelated biopsy-stage contrasts as cross-sectional stage-associated
  remodeling. Fixed continuum windows are descriptive, and NMF axes are continuous.
- Keep the 117 program memberships and weights fixed. A fixed-weight projection
  is a different measurement from historical neighborhood-dependent Hotspot scores.
- Apply current variant-identity and LD qualifications to all dependent genetic
  claims. Colocalization alone does not establish a causal gene.
- Keep named release states explicit and do not combine historical and updated
  values silently. Scientific adoption follows the [roadmap checklist](ROADMAP.md#final-release-checklist).
- Keep future experimental follow-up outcomes outside Resource selection and claims.

## What happened to the older documents?

The previous guide's dated updates and detailed routes are preserved in the
[November 9 snapshot](archive/folder_navigation_2026-11-09/DOCS_README.md).
The [archive index](archive/INDEX.md) also links the earlier documentation
consolidation. Archived content explains prior work and does not override the
current five authorities. Data, named releases and analysis paths remain in place.
