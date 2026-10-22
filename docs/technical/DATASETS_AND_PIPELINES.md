# Dataset and pipeline index

This is a routing document, not a replacement for code, configuration, or
release manifests. Use it to find the relevant producer without reading the
entire historical current-state handoff.

## Source authority

| Question | Authority |
|---|---|
| Which human or mouse cohort is retained? | `config/human_datasets.yaml` or `config/mouse_datasets.yaml` |
| Which thresholds and reference paths apply? | `config/pipeline_params.yaml` |
| Which exact inputs produced a paper value? | Named release/candidate manifest |
| Which number is currently usable? | [../RESULTS.md](../RESULTS.md) |
| Which work is still open? | [../STATUS.md](../STATUS.md) |

## Primary analysis routes

| Domain | Canonical producer or root | Current Resource role | Environment |
|---|---|---|---|
| Human bulk integration | `RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/05h_limma_voom_qw_canonical.R` | Disease-state bulk substrate; fragment-count successor is being synchronized | `rnaseq` |
| Cross-sectional stage | `scripts/manuscript/noncoding_resource/build_stage_lncrna_candidate.R` and fragment candidate root | Validated Figure 3 cross-sectional stage candidate; synchronized adoption pending | `rnaseq` |
| Histology-anchored molecular continuum | `scripts/analysis/histology_anchored_continuum/` and `scripts/analysis/histology_anchored_continuum/molecular_layers/` | Kamzolas-adapted co-primary scores; candidate Figure 4F and complete S3/S4 transcript, pathway, TF, NMF, Hotspot, and paired comparisons; Figure 3 unchanged | `rnaseq` |
| Cell composition | Human integration deconvolution outputs | Figure 3 donor/cohort-aware composition; rebuild pending | `rnaseq` |
| Gene-level integration table | `RNA-seq/27a → 45a → 75 → 217` and guarded downstream chain | MASLD Gene Catalog substrate; no universal validity score | `rnaseq` |
| Regulatory genetics | `GWAS/finemapping/src/` and `RNA-seq/results/causal_inference/` | Figure 2 shared-signal candidates and context | `rnaseq`; `motifbreakr` for motif step |
| Single-cell map | `Analysis/SingleCell/` | Donor-aware lineage and state context | `rapids_singlecell` plus assay-specific environments |
| Single-cell ATAC context | `Analysis/Multimodal_Program_Projection/scripts/atac_context_v3/` with candidate root `Analysis/Multimodal_Program_Projection/candidates/atac-context-v3-candidate-2026-08-11-r1/` | Condition-blind cross-cohort peaks, exact donor fragment-record counts, matched DA, all-program promoter projection, official chromVAR, and a variant-observability replay run on fixed fixtures and conditioned on the adopted COLOC | `snapatac2`, `rnaseq`, dedicated pinned chromVAR environment |
| Frozen Hotspot programs | `Analysis/SingleCell/scripts/hotspot_modules/500–509` | 117-program coordinate system | `rapids_singlecell`, `rnaseq` |
| Multicellular programs | `Analysis/SingleCell/scripts/400_mcp/` | Supplementary and contextual program analyses | `rapids_singlecell`, `celltype_bio` |
| Spatial context | `Analysis/Spatial/` and multimodal projection adapters | Figure 4 physical context | `spatial` |
| Protein context | Tissue and secretome candidate adapters | Figure 4 assay-native transportability | assay-specific |
| Mouse context | `RNA-seq/Mouse/Unified_Integration/scripts/M00–M07` | Contextual cross-species analysis, not a human causal proof | `rnaseq` |
| lncRNA and noncoding interpretation | `scripts/manuscript/noncoding_resource/` and `scripts/portal/gene_catalog_v2.py` | Five-cohort lncRNA measurement, genomic class, program sensitivity, fine-mapped DNA adapter, and molecular-object routing | `rnaseq`, `alphagenome` |
| MASLD Gene Catalog | `Analysis/Multimodal_Program_Projection/` and portal exporters | Figure 5 and gene follow-up guidance | existing project environments |

## Important boundaries

- The human bulk model is limma-voom quality-weighted with fixed cohort effect;
  the updated fragment-count substrate must be propagated coherently.
- TREAT Tier 1 uses an interval null at `|logFC| ≤ 0.25`; do not apply an extra
  post-hoc fold-change filter.
- Single-cell inference is donor-level. Runs, sorted fractions, and cells do not
  increase biological sample size.
- Documented Andrews fibrosis stages cover 46 donors represented by 58 runs;
  cross-cohort inferred fibrosis stage is not authorized.
- COLOC/TWAS are the active regulatory-genetic framework. MR and INTACT are
  retired. BBJ analyses retain the limitation that liver eQTLs are European.
- The 117 Hotspot modules are frozen; public outcomes cannot trigger discovery
  or reweighting.
- NMF k4/k6 is interpretive and continuous. It is not a patient-classifier claim.
- The molecular continuum is stage-adjusted and cross-sectional. Its nine
  fixed sliding windows are descriptive only; they cannot determine
  significance, select features, orient scores, or imply temporal ordering.
- The historical convergence heuristic may remain a supplementary navigation
  field but cannot determine default portal ranking or target validity.

## Current high-risk producer notes

| Risk | Rule |
|---|---|
| Old human wrapper runs retired dream and omits the canonical producer | Do not use `run_full_5cohort_slurm.sh` for a paper release until repaired |
| Old mouse wrapper carries retired five-diet logic | Use `run_integration_4model.sh` |
| Old integration/portal chain can publish stale ranking semantics | Keep legacy portal generators fail-closed until the Gene Catalog release |
| Historical `RNA-seq/53–57_ncrna_*` workflows mix ncRNA classes, live queries, and outcome-selected analyses | They fail closed unless `ALLOW_HISTORICAL_NCRNA_PIPELINE=1`; never use them to produce Resource results |
| Figure scripts can regenerate legacy panels | Respect `ALLOW_LEGACY_*` guards and read-only figure manifests |
| Updated stage source is only a preflight | Run contrast-native filtering, normalization, and models before any stage DEG claim |

## Detailed history

The full 1,141-line current-state handoff, including every historical pipeline
correction and result path, is preserved at
`docs/archive/documentation_consolidation_2026-08-11/originals/docs/CODEBASE_CURRENT_STATE.md`.
Use it only when a technical question is not answered by this index, the code,
or a frozen manifest.

The end-to-end single-cell implementation history is retained separately in
[SINGLE_CELL_HISTORY.md](SINGLE_CELL_HISTORY.md). Its overview counts predate
some donor-collapse corrections; use [../RESULTS.md](../RESULTS.md) for current
paper numbers.
