# Release and provenance guide

## Three release states

| State | Human meaning | Allowed use |
|---|---|---|
| **Current official release** | Frozen package currently recognized as official | Exact claims are allowed if the release is named |
| **Validated update awaiting release** | Analysis passed its own checks but downstream outputs are not synchronized | Candidate rebuilding and clearly labeled comparison only |
| **Working candidate** | Exploratory or partially assembled result | Internal review only; no manuscript claim |

“Validated” does not mean “adopted.” Adoption is a coordinated change across
inputs, derived analyses, figures, manuscript numbers, Gene Catalog, portal data,
and release manifests.

## Current layers

| Layer | Status | Notes |
|---|---|---|
| July analysis package (`2026-07-15-r2`) | Current official historical baseline | Contains the 846-sample/1,918-DEG read-count substrate; incomplete for the present Resource requirements |
| Updated five-cohort bulk analysis | Validated update selected for the paper | 844 samples and 1,616 TREAT DEGs; downstream stage/composition/figures still need synchronization |
| Histology-anchored molecular continuum | Validated candidate; not adopted | `hac-molecular-layers-20260818T173348Z` contains complete molecular families and the Figure 4F checks; `histology-vs-continuum-comparison-20260818T192630Z` contains matched stage-versus-continuum figures. Figure 3, canonical DEGs, NMF definitions, and Hotspot weights are unchanged |
| Fragment stage source | Preflight only | 661 biological samples; no inferential model or DEG counts yet |
| Corrected 50-study COLOC portfolio | In progress | Writes outside current official products; exact memberships await audit |
| MASLD Gene Catalog export | Working candidate | Structure is useful; genetics-dependent entries require rebuild |
| Final Resource release | Not assembled | Must satisfy the checklist in [../ROADMAP.md](../ROADMAP.md) |

## Adoption sequence

1. Freeze source inputs and code by checksum.
2. Reproduce source sample counts, biological units, and expected headline
   directions before projecting Resource objects.
3. Run analysis-native validation and complete statistical families.
4. Perform independent number and biological-unit audits.
5. Rebuild all dependent products from the selected upstream release.
6. Generate figures only from read-only source tables.
7. Assemble numbers, claims, methods, data dictionary, Gene Catalog, portal, and
   environment manifests under one candidate release ID.
8. Run scope, consistency, and link validators.
9. Obtain explicit user approval, then adopt it without overwriting historical
   releases.

## What a complete manifest must record

- release ID and parent release;
- file path, byte size, SHA256, source URL/accession, and retrieval date;
- producing script, code commit or snapshot hash, environment, and command;
- dataset, assay, biological unit, sample count, exclusions, and contrast;
- tested universe, multiple-testing family, effect unit, and robustness state;
- upstream registry/bulk/COLOC hashes;
- figure and manuscript consumers;
- source dependence and known limitations.

## Rules for keeping releases unchanged

- Never overwrite a named official release or the v1 22-program multimodal
  package.
- Never mix source tables from two releases in one figure.
- Preserve failed, skipped, null, and indeterminate branches; biological
  positivity is not a software acceptance criterion.
- A live mutable path is not paper provenance. Copy accepted output files into a
  read-only release root and checksum them.
- Exact live job status comes from SLURM; dated progress prose is not a monitor.

## Historical package warning

`docs/manuscript/release/` describes July material. Its existing `PASS` means
only that the historical package met its then-current checks. It does not prove
that the current five-figure Resource, updated bulk substrate, corrected COLOC,
Gene Catalog, portal, or separation from the experimental follow-up study are complete.

## Recovery and audit

Superseded authority documents, plans, and progress logs were copied verbatim
to `docs/archive/documentation_consolidation_2026-08-11/originals/`. The archive
contains source SHA256 values and a migration manifest; see
[../archive/INDEX.md](../archive/INDEX.md).
