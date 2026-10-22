# Single-Cell Analysis — End-to-End Documentation

**Last update:** 2026-05-24 (post-protocol-contamination remediation)
**Status:** Architecture and remediation history. Supersedes `scrna_pipeline_rigor_plan.md` and `ccc_review/*` (archived under `docs/archive/single_cell_pre_remediation_2026-05-22/`), but its overview counts and run-level significance are not current.

> **CURRENTNESS NOTICE (2026-07-13):** The 2026-07-12 donor-collapse is complete
> and supersedes this document's 269-donor / 211-stage-donor overview. The
> canonical replication unit is **125 biological donors from 273 analyzed
> libraries**. Current donor-level counts, withdrawn claims, and release impact
> are summarized in `docs/STATUS.md`, `docs/RESULTS.md`, and
> `docs/technical/NUMBERS_HISTORY.md`. In particular, do not
> cite the run-level LIANA 4→10→16 count, the 13-pair CCC FDR claim, or the
> cross-cohort inferred/augmented F-stage axis. Keep this file for pipeline DAG,
> dataset provenance, and the 2026-05-22 protocol-contamination history.

Paths and “CURRENT” labels quoted later in this historical record describe its
2026-05-24 state. Retired manuscript paths are recoverable through
`docs/archive/DOCUMENT_MIGRATION_MANIFEST.tsv`; they are not active documents.

---

## 1. Overview

The MASLD single-cell atlas integrates 895,542 cells from 269 donors across seven datasets, processed through a unified scVI-based latent embedding with cell-type annotation, hotspot autocorrelation module discovery, ligand-receptor (LIANA) communication inference, multi-cellular program (cNMF / DIALOGUE) decomposition, hepatocyte subclustering, and pseudotime trajectory inference. After the 2026-05-22 protocol-contamination remediation, **211 donors (~796,000 cells) from five whole-tissue or paired-modality datasets** contribute to stage-stratified analyses; the remaining 58 donors (GSE136103 + Liver_Atlas) are flagged `exclude_stage_analysis = TRUE` in donor metadata and excluded from cross-stage comparisons.

The atlas centerpiece in the manuscript is **Figure 2** (hepatocyte hotspot cascade, MASH-vs-Healthy UMAP density, LIANA stage-CCC heatmap, cell-type composition trajectory).

---

## 2. Datasets and protocol groups

| Group | Dataset | Paper | Donors | Hep % | Status |
|-------|---------|-------|--------|-------|--------|
| **Whole-tissue** ✅ | GSE244832 | Kim, Rosenthal et al., *J Hepatol* 2025 (PMID 39522884) | 117 runs / **18 donors**† | 82% (44-100%) | VALID |
| **Whole-tissue** ✅ | GSE185477 | (unattributed in atlas) | 21 | 81% (70-90%) | VALID |
| **Whole-tissue** ✅ | GSE174748 | (unattributed in atlas) | 4 | 78% (72-84%) | VALID |
| **NPC-enriched / FACS** ❌ | GSE136103 | Ramachandran et al., *Nature* 2019 (PMID 31597160) | 20 | 0% (NPC enrichment by gravity sedimentation) | EXCLUDED |
| **NPC-enriched / FACS** ❌ | Liver_Atlas / GSE192740 | Guilliams et al., *Cell* 2022 (PMID 35021063) | 38 | 20% bimodal (29 FACS-sorted NPC, 9 hepatocyte-enriched) | EXCLUDED |
| **Paired sc + sn** ⚠️ | GSE202379 | Gribben et al., *Nature* 2024 (PMID 38778114) | 67 | 73% (29-100%, SD=21%) | VALID (snRNA-only per primary publication; previously mislabeled as Andrews 2022 sc+sn) |
| **Hep-specific** | GSE189600 | (unattributed in atlas) | 2 | 100% (Percoll/EpCAM hepatocyte isolation, 142 cells total) | VALID but n=2 |

> **† GSE244832 run-vs-donor correction (2026-07-09):** the 117 GSE244832 rows in `donor_metadata.tsv` are 117 SRR sequencing **runs** that map to only **18 biological donors** (D01–D18; verified 117/117 against the paired multiome `Analysis/ATAC/Human_Multiome/metadata/donor_metadata_curated.tsv`, ~6.5 runs/donor). This is the snRNA arm of the Kim/Rosenthal paired snRNA + snATAC HSC multiome — NOT a "Wang et al." dataset (prior mis-attribution). Treating runs as donors **pseudo-replicates this cohort ~6.5×**; because GSE244832 dominates the disease arm (~88% of progressor cells), every donor-level statistic that counts 117 GSE244832 "donors" is inflated and should aggregate to the 18 real donors. Compounds the megareview single-cohort/pseudoreplication findings for this dataset. The atlas totals below (269 donors, 211 stage donors) still count runs; not yet re-derived to donor level.

**Exclusion rationale**: NPC-enriched and FACS-sorted preparations deliberately remove hepatocytes (Ramachandran 2019) or select for specific immune subsets (Guilliams 2022, 22/38 donors are Kupffer-cell-sorted). Pooling these donors with whole-tissue captures distorts disease-stage comparisons of cell-type composition and donor-level module scores. The exclusion is implemented at the donor-metadata layer (column `exclude_stage_analysis`) and inherited by every downstream script that joins donor metadata.

**Full per-dataset table**: see `docs/reference/scrna_dataset_protocol_comparison.md` (per-donor stats, paper citations, sample-prep notes).

---

## 3. Pipeline architecture

### Script number ranges

| Range | Phase | Purpose |
|-------|-------|---------|
| `300–307` | Pseudotime | Palantir + DPT + CytoTRACE v1 + PAGA per cell type |
| `308–312` | Hepatocyte subclustering | scVI on hep subset, 5 meta-subtypes (Progressor/Moderate/Stable/Neutral/Healthy), cross-modal validation |
| `319–339` | Cell-cell communication pre-stage | Hep-stromal circuits, metabolic CCC, chromatin coupling, sex-stratified CCC, cross-species CCC |
| `343*` | F-stage harmonization | Documented (Andrews), scVI ordinal (343b), augmented (343m, conditional on the bootstrap check) |
| `344` | **Donor metadata builder** | Single source of truth for stage-stratified analyses (`donor_metadata_extended.tsv`) |
| `345–353` | LIANA stage-trajectory CCC | Per-donor LR scoring, merge, LMM, LOO, bulk concordance, bootstrap, figures |
| `400_mcp/` | Multi-cellular programs | cNMF (Kotliar 2019) global + per-cell-type, DIALOGUE (Jerby-Arnon 2022), changepoint detection |
| `500–510` (`hotspot_modules/`) | Hotspot autocorrelation | Per-cell-type module discovery, scoring, LOO stability, phenotype correlation, atlas integration, CRN transitions |

### Per-cell-type Hotspot run order

```
501 (discover modules per CT) → 502 (aggregate donor scores)
                              → 503 (LOO stability, dataset hold-outs)
                              → 503b (consolidate LOO)
                              → 504 (novelty vs cNMF/Hallmark/SCENIC+)
                              → 505 (phenotype correlations, LMM)
                              → 506 (bulk replication)
                              → 507 (Olink + Visium validation)
                              → 508 (atlas integration → all_modules.tsv)
                              → 510 (CRN transition scores)
                              → fig2_panel_hotspot_cascade.R
```

### Per-LR stage-CCC run order

```
345 (per-donor LIANA) → 345b (merge into all_donor_lr_scores.tsv.gz)
                     → 346 (stage LMM with frac_Hep covariate for hep pairs)
                     → 347 (LOO replication, leave-one-dataset-out)
                     → 348 (bulk concordance)
                     → 348b/348_pseudobulk (per-cell-type F-stage Dream)
                     → 350 (bootstrap CI)
                     → 349/352/353 (figures + trajectory heatmaps)
```

### Critical orchestrators

- `344_build_donor_metadata_extended.R` — rebuilds `donor_metadata_extended.tsv` with all stage / pseudotime / proportion / program / demographic columns. **Re-run after any upstream metadata change.**
- `502_score_modules_in_atlas.py` — iterates RUN_ORDER over per-CT directories and aggregates `donor_scores_all.tsv` + `cell_scores_all.parquet`. Skips missing directories with `[WARN]`.

---

## 4. Hotspot module discovery (cascade — Fig 2 centerpiece)

### Method

Per-cell-type Hotspot (DeTomaso & Yosef 2021) on scVI / Harmony-PCA latent: DANB model, k=30 (k=50 for macrophages post-remediation), autocorr Z ≥ 7 / FDR < 0.05, max 2000 features, ≥ 20-gene modules.

### Cell-type-specific runs (post-2026-05-22 remediation)

| Cell type | h5ad cells (clean) | Modules | k | Status |
|-----------|--------------------|---------|---|--------|
| Hepatocytes | ~657,412 | 30 | 30 | OLD outputs retained (NPC datasets contribute 0 hep cells; no contamination) |
| Macrophages | 19,595 | 19 | **50** | Re-discovered 2026-05-22 (k raised from 30→50 after Liver_Atlas exclusion dropped cell count from 42,371) |
| Fibroblasts | 38,439 | 29 | 30 | Re-discovered 2026-05-22 |
| Cholangiocytes | 32,402 | 28 | 30 | Re-discovered 2026-05-22 |
| T cells | 42,295 | 11 | 30 | Re-discovered 2026-05-22 |
| Endothelial | (deferred) | — | — | Old contaminated outputs moved aside as `endothelial_cells_old_skip_502/` |
| Global (cross-CT) | (deferred) | — | — | Old contaminated outputs moved aside as `global_old_skip_502/` |

Endothelial and global re-discoveries were deferred because (a) `fig2_panel_hotspot_cascade.R` filters to `cell_type == "hepatocytes"` and is unaffected, and (b) the local-correlation step on 93K (endothelial) or 855K (global) cells was running >1h with no completion bound. Old outputs are retained for downstream scripts that need any entry for these cell types; their modules are flagged as pre-remediation in `all_modules.tsv` via `mtime`.

### Focal hepatocyte modules (manuscript headline)

All five focal modules' β coefficients are quantitatively unchanged (rounding-level shifts) after the protocol exclusion — these modules were discovered on the (already protocol-clean) hepatocyte subset.

| Module | Direction | F-stage β | F-stage q | Disease-stage q | Top genes |
|--------|-----------|-----------|-----------|-----------------|-----------|
| **Hep-20** | GAIN | +0.36 | 2.6 × 10⁻¹² | **6.4 × 10⁻²⁰** | BICC1, GLS, ITGAV (glutamine catabolism + TGF-β integrin) |
| **Hep-24** | GAIN | +0.36 | 3.4 × 10⁻¹¹ | 5.0 × 10⁻¹⁴ | SOD2, TXNRD1, SQSTM1, AKR1B10, **HKDC1** (NRF2 antioxidant + GWAS-colocalizing hexokinase) |
| **Hep-26** | GAIN | +0.57 | 2.4 × 10⁻⁹ | 3.6 × 10⁻⁷ | JUN, ATF3, KLF6, SERPINE1, GDF15 (AP-1 acute injury) |
| **Hep-19** | LOSS | −0.36 | 4.4 × 10⁻⁷ | 1.2 × 10⁻⁸ | HNF4A-AS1, HNF1A-AS1, EHHADH, BDH1, PECR, ACOX2 (HNF identity + peroxisomal FAO) |
| **Hep-27** | LOSS | −0.48 | 9.3 × 10⁻⁹ | 3.7 × 10⁻³ | CFH, CFHR4, C8A, C5, NR1H4/FXR (complement + bile-acid synthesis) |

### Hotspot module summary (clean atlas)

- 30 total hepatocyte modules (autocorrelation Z ≥ 7, FDR < 0.05, ≥ 20 genes)
- 17 of 30 sig at F_stage q < 0.05
- 13 of these tagged as progression modules
- 9 of 13 progression modules bulk-replicated (≥ 80% concordant-directional genes in bulk mega-analysis)

---

## 5. LIANA cell-cell communication (stage-CCC)

### Method

Per-donor LIANA (4 methods: SCA, connectome, NATMI, LogFC) scored on the clean atlas (211 donors × ~518 LR pairs after MIN_DONORS_PER_LR = 30 filter). Mixed-effects LMM `(magnitude_rank ~ disease_stage + (1|dataset))` for each axis (`disease_stage_coarse`, `disease_stage_numeric` continuous, `F_stage_inferred`, macrophage_pseudotime_mean). For Hep-involving LR pairs (source or target == Hepatocytes), `frac_Hepatocytes` is added as a fixed-effect covariate to deconfound 0%-vs-80%-hep composition variation.

### Bootstrap stability check

The pre-registered bootstrap stability check (`BOOTSTRAP_FALLBACK_REQUIRED.flag` at `Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/`) flags that `F_stage_augmented` failed the F3-mode preservation criterion (27.6% < 95% target). When the flag is present, `346_stage_ccc_mixedmodel.R` and `510_crn_transition_scores.py` fall back to `F_stage_inferred` (Script 343b scVI ordinal classifier). Override env var: `ALLOW_F_STAGE_AUGMENTED_DESPITE_FLAG=TRUE` (requires peer sign-off; document in methods).

### Headline paracrine pairs (clean LMM, q < 0.05)

| LR pair | Cell-type pairs | Min q across axes |
|---------|-----------------|-------------------|
| **CDH1–PTPRM** | 3 coarse / 1 continuous / 2 fstage | 2.2 × 10⁻⁶ |
| **COL18A1–KDR** | 4 / 1 / 3 | 4.4 × 10⁻⁴ |
| **ST6GAL1–EGFR** | 6 / 7 / 2 | 1.7 × 10⁻⁹ |

All three survive at padj_within_ct < 0.05 in the clean LMM.

### Total clean LMM hits

- 518 LR pairs at coarse-stage axis (q < 0.05)
- 541 at continuous-pseudotime axis
- 249 at F-stage-numeric axis

### CRN transition counts (per cell-type module / LR pair)

- F0 → F1: 14 of 117 sig (q < 0.05)
- F1 → F2: 27 of 117
- F2 → F3: 74 of 117
- F3 → F4: 61 of 117 — **secondary axis** since 19 cirrhosis donors all derive from GSE202379

---

## 6. Multi-cellular programs (cNMF + DIALOGUE)

### cNMF (Kotliar 2019)

Global k=16 + per-cell-type cNMF programs at `Analysis/SingleCell/results_gpu_v2/mcp/`. Programs scored per donor; donor scores joined into `donor_metadata_extended.tsv` as `cnmf_global_k16_P1`, `_P5`, `_P11`, `_P13`, `_P14`.

### Post-remediation changepoint analysis

The clean re-run of `461_changepoint.R` (and `460_phenotype_screen.py`) on 211 donors gave a **MAJOR sensitivity verdict**:

- 4 of 16 programs (P2, P3, P11, P15) lost their `seg1` (one-breakpoint) classification — these were Steatosis → SH breakpoints driven by the contaminated GSE136103 Cirrhosis endpoint.
- 13 of 16 programs now best-fit by GAM (monotonic trajectories) — `seg1: 4 → 0`, `gam: 10 → 13`, `linear: 2 → 3`.
- 17 of 32 phenotype rows flip significance (9 lost q<0.05, 8 gained). Lost: P9 disease (q 2e-15 → 0.86), P11 disease (1.2e-4 → 0.86), P15 (8e-5 → 0.82).
- **P6 and P7 (strongest disease_stage_numeric programs) retain significance and gain larger effect sizes** in the clean run — best candidates for future fig2 cNMF expansion panels.

Reference: `Analysis/SingleCell/results_gpu_v2/mcp/protocol_sensitivity_summary.md`.

### DIALOGUE

Multi-cellular co-variation programs identified across 9 phenotypes × 4 k values (`Analysis/SingleCell/scripts/400_mcp/430_dialogue_run.R`). DIALOGUE outputs are not currently in fig2 manuscript text; supplementary potential.

---

## 7. Hepatocyte subclustering (Scripts 308-312)

### 5 meta-subtypes

| Meta-subtype | Biology |
|--------------|---------|
| **Progressor** | F2-expanding population; signatures include disease-associated injury and TGF-β signaling |
| Moderate | Mid-stage, transitional |
| Stable | Early-stage, partially preserved metabolic identity |
| Neutral | Cohort-mean state |
| Healthy | Preserved hepatocyte-intrinsic metabolic identity |

### Progressor expansion across stages

| Stage | Progressor fraction |
|-------|---------------------|
| Healthy | 0.0% |
| Steatosis | 2.25% |
| Steatohepatitis | **56.54%** (50-fold expansion vs Steatosis) |

This trajectory is **byte-identical** between contaminated and clean sets, because Steatosis and SH donors are 100% whole-tissue (GSE244832 [Kim/Rosenthal] + Andrews + GSE174748) — neither GSE136103 nor Liver_Atlas contributes Steatosis or SH donors.

### Sensitivity check (Phase 6.1, 2026-05-22)

- Healthy ↔ Steatohepatitis centroid cosine shift in scVI latent: **0.79%** (threshold for action: 5%) → **ROBUST**
- Output: `Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/protocol_sensitivity_check.{tsv,md}`

---

## 8. Pseudotime (Scripts 300-307)

Four methods consensus per cell type: **Palantir + DPT + CytoTRACE v1 + PAGA** (Monocle 3 dropped as weakest method). Output: `consensus_pseudotime_all.csv` (per-cell pseudotime + per-method aligned values).

### Per-cell-type Δ (max − min mean pseudotime across stages)

Pre-remediation values:

| Cell type | Δ | F3→F4 ρ |
|-----------|---|---------|
| Macrophages | 0.73 | 0.39 |
| Hepatocytes | 0.62 | (lower) |
| Endothelial | 0.24 | (lower) |

### Sensitivity check (Phase 6.2, 2026-05-22) — UPDATE REQUIRED

Re-computing donor-level macrophage pseudotime × F_stage Spearman correlation **after excluding GSE136103 + Liver_Atlas donors**:

| Set | n donors | ρ | p | Δ |
|-----|----------|---|---|---|
| Full (contaminated) | 242 | +0.017 | 0.79 | 0.347 |
| **Clean** | 193 | **−0.224** | **1.7 × 10⁻³** | 0.153 |

**Sign flip** indicates the original positive trajectory was driven by Liver_Atlas FACS-sorted Kupffer-cell preparations (mean pseudotime 0.384 vs 0.90 in whole-tissue macrophages) inflating the F0/Healthy bin with low-PT donors.

Caveat: the previously cited "F3→F4 ρ = 0.39" was a **gene-signature × cell-pseudotime** metric from Script 307g, not the donor-mean × F-stage ρ recomputed here. The donor-mean number was effectively never published.

Output: `Analysis/SingleCell/results_gpu_v2/pseudotime/macrophage_pseudotime_protocol_sensitivity.tsv`

---

## 9. Protocol contamination remediation (2026-05-22)

The remediation was triggered by investigating an anomalous macrophage Healthy → SH composition decrease (6.6% → 2.2%) which turned out to be a FACS-sort artifact. Full investigation revealed three incompatible sample-prep protocols had been co-mingled without flagging.

### Eight-phase remediation (24 tasks, all completed)

1. **Phase 0 — F-stage check fix**: `510_crn_transition_scores.py` now respects `BOOTSTRAP_FALLBACK_REQUIRED.flag` (falls back to `F_stage_inferred` when present).
2. **Phase 1 — Donor metadata**: `344` extended with three columns (`protocol_group`, `exclude_stage_analysis`, `F_stage_augmented_clean`). 58 donors flagged (20 GSE136103 + 38 Liver_Atlas).
3. **Phase 2 — CCC parquet archival**: 48 contaminated per-donor LIANA parquets moved to `excluded_protocol_contaminated/`; `345b` re-merged `all_donor_lr_scores.tsv.gz` (864,121 → 576,613 rows, 33% reduction).
4. **Phase 3 — Hotspot re-discovery**: `501_run_hotspot.py` gained `--exclude-datasets` CLI flag; 4 of 6 non-hep cell types re-discovered cleanly (mac, fib, chol, tcells); 2 deferred (global, endothelial).
5. **Phase 4 — LIANA stage-CCC**: All scripts in the 346 → 350 chain re-run with donor-level exclusion filter. `frac_Hepatocytes` covariate added to 346 for Hep-involving LR pairs.
6. **Phase 5 — UMAP / panels**: `gen_fig2e_embeddable.R` KDE filter expanded to drop GSE136103 + Liver_Atlas explicitly; NMF violin + composition panels filter to whole-tissue datasets.
7. **Phase 6 — Sensitivity checks**: 3 independent sensitivity analyses (Phase 6.1 hep subcluster ROBUST, Phase 6.2 mac pseudotime UPDATE-NEEDED, Phase 6.3 cNMF changepoint MAJOR).
8. **Phase 7 — Figure regeneration**: All Phase 3 + Phase 4 + Phase 7 figures regenerated. Key panel: `fig2_panel_hotspot_cascade.pdf` re-rendered with 3-stage (Healthy/ST/SH) cascade.
9. **Phase 8 — Manuscript audit + methods text**: `02_results_fig2_audit_2026_05_22.md` (14 change items) and `02_methods_scrna_protocol_audit.md` (methods + reviewer-defense paragraph).

### What changed in the manuscript

Six text patches applied to `docs/manuscript/working/02_results_fig2.md` (backup at `02_results_fig2.md.bak_pre_protocol_fix_patches`):
1. L6 — added 2.3% Healthy macrophage fraction; FACS-exclusion caveat
2. L8 — atlas count 275 → 269 donors / 211 clean / ~796K cells; per-stage breakdown
3. L8 (CCC) — 3-stage Healthy → ST → SH cascade; F3→F4 secondary axis caveat
4. L13 — 16/5 → 17/13/9 (sig modules / progression-tagged / bulk-replicated)
5. L15 — q-values refreshed for Hep-20/24/26
6. L17 — q-values refreshed for Hep-19/27
7. **New paragraph** after L17 — sensitivity findings (robust hep modules, mac pseudotime sign-flip, cNMF GAM-monotonic).

### Andrews/Gribben dataset correction

Internal project memory previously labeled GSE202379 as "Andrews 2022 paired sc + sn". Phase 3 Agent 6 verified against the primary publication (Gribben et al., *Nature* 2024, PMID 38778114, DOI 10.1038/s41586-024-07465-2) that the dataset is **snRNA-only** (all 47 patients, 67 SRRs). The frac_Hepatocytes range (29-100%) reflects tissue-sampling variation, not modality variation. Mapping at `Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/andrews_modality_mapping.tsv`.

---

## 10. Sensitivity findings — three-tier summary

| Finding | Verdict | Number |
|---------|---------|--------|
| 5 focal hepatocyte modules (Hep-19, -20, -24, -26, -27) β coefficients | **ROBUST** | All within rounding |
| 3 headline LIANA paracrine pairs (CDH1-PTPRM, COL18A1-KDR, ST6GAL1-EGFR) | **ROBUST** | All survive q < 0.05 |
| Hepatocyte subclustering Progressor expansion | **ROBUST** | 2.25% → 56.54% byte-identical |
| Healthy ↔ SH scVI centroid distance | **ROBUST** | Cosine shift 0.79% < 5% threshold |
| Donor-level macrophage pseudotime × F-stage correlation | **SIGN FLIP** | ρ +0.017 → −0.224 (p = 1.7×10⁻³) |
| cNMF program changepoints (P2, P3, P11, P15) | **REFRAMED** | seg1 → GAM monotonic (4 of 16 lost breakpoints) |
| cNMF programs P6, P7 effect sizes | **STRONGER** | Best candidates for future expansion |
| Donor count for stage analyses | **UPDATED** | 275 → 211 |
| Cell count for stage analyses | **UPDATED** | 1,232,318 → ~796,510 |
| Macrophage Healthy fraction | **UPDATED** | 6.6% → 2.30% |

---

## 11. Key output files (canonical paths)

### Donor metadata
- `Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/donor_metadata_extended.tsv` — 269 donors × 72 columns. **Single source of truth** for stage analyses. Built by `344_build_donor_metadata_extended.R`. Key columns added 2026-05-22: `protocol_group`, `exclude_stage_analysis`, `F_stage_augmented_clean`.
- `Analysis/SingleCell/results_gpu_v2/mcp/inputs/donor_metadata.tsv` — minimal donor metadata used by some 400_mcp scripts (does NOT have `exclude_stage_analysis` — affected scripts must join the extended file or apply equivalent filter).
- `Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/andrews_modality_mapping.tsv` — GSE202379 modality clarification (all snRNA per Gribben 2024).

### Atlas h5ads
- `Analysis/SingleCell/integration/output/human/scalesc_human_annotated_celltypist.h5ad` — Master integrated atlas (1.2M cells). Used as ATLAS_GLOBAL in `hotspot_modules/hotspot_io.py`.
- `Analysis/SingleCell/results_gpu_v2/mcp/inputs/atlas_cnmf_{hepatocytes,macrophages,fibroblasts,endothelial_cells,cholangiocytes,tcells}.h5ad` — Per-cell-type subsets (used by Hotspot, cNMF, hep subclustering).

### Hotspot module results
- `Analysis/SingleCell/results_gpu_v2/hotspot_modules/{cell_type}/` — Per-cell-type module discovery outputs (`autocorr.tsv`, `module_genes.tsv`, `cell_scores.parquet`, `donor_scores.tsv`, `hotspot_obj.pkl`, `run_metadata.json`).
- `Analysis/SingleCell/results_gpu_v2/hotspot_modules/all_modules.tsv` — Master module catalog with disease_stage_beta, F_stage_beta, novelty, bulk_replicated flags. Refreshed 2026-05-22 by `508_atlas_integration.R`.
- `Analysis/SingleCell/results_gpu_v2/hotspot_modules/donor_scores_all.tsv` — Per-donor per-module scores (namespaced as `<celltype>__<module_int>`). Refreshed 2026-05-22 by `502_score_modules_in_atlas.py`.
- `Analysis/SingleCell/results_gpu_v2/hotspot_modules/crn_transition_scores.tsv` — F0→F1, F1→F2, F2→F3, F3→F4 per module. Refreshed 2026-05-22 by `510_crn_transition_scores.py`.
- `Analysis/SingleCell/results_gpu_v2/hotspot_modules/phenotype_correlations.tsv` — LMM β + q per (cell_type, module, axis). Refreshed 2026-05-22.
- `Analysis/SingleCell/results_gpu_v2/hotspot_modules/loo_stability.tsv` — LOO stability per dataset (note: only mac/fib/chol/tcells were re-run; hep + endo + global retain pre-remediation values).
- `Analysis/SingleCell/results_gpu_v2/hotspot_modules/{cell_type}_old_skip_502/` — Archived pre-remediation outputs for global + endothelial.

### LIANA stage-CCC
- `Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/per_donor_lr/*.parquet` — Per-donor LR scores (active set, 234 parquets post-remediation). Excluded 48 archived under `excluded_protocol_contaminated/`.
- `Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/all_donor_lr_scores.tsv.gz` — Merged LR-pair × donor table (576,613 rows post-remediation). Built by `345b`. Backup: `all_donor_lr_scores.tsv.gz.bak_pre_protocol_fix`.
- `Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/stage_lr_lmm_{coarse,continuous,fstage}_v3.tsv` — Stage-LMM outputs per axis.
- `Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/bootstrap_ci.tsv` — Bootstrap CI per LR pair.
- `Analysis/SingleCell/results_gpu_v2/ccc/stage_trajectory/BOOTSTRAP_FALLBACK_REQUIRED.flag` — Pre-registered check flag. If present, `F_stage_inferred` is used in place of `F_stage_augmented`.

### cNMF / DIALOGUE
- `Analysis/SingleCell/results_gpu_v2/mcp/integration/program_donor_scores_wide.tsv.gz` — Donor-level program scores.
- `Analysis/SingleCell/results_gpu_v2/mcp/integration/phenotype_screen_long_clean.tsv` — Clean phenotype screen.
- `Analysis/SingleCell/results_gpu_v2/mcp/integration/program_changepoint_clean.tsv` — Clean changepoints.
- `Analysis/SingleCell/results_gpu_v2/mcp/protocol_sensitivity_{comparison_changepoint,comparison_phenotype,summary}.{tsv,md}` — Phase 6.3 sensitivity outputs.

### Hep subclustering
- `Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/hepatocyte_subtype_metadata.csv` — Per-cell subtype assignments.
- `Analysis/SingleCell/results_gpu_v2/hepatocyte_subtypes/protocol_sensitivity_{check,donor_table,stage_compare}.tsv` + `.md` — Phase 6.1 sensitivity outputs.

### Pseudotime
- `Analysis/SingleCell/results_gpu_v2/pseudotime/consensus_pseudotime_all.csv` — Per-cell consensus pseudotime.
- `Analysis/SingleCell/results_gpu_v2/pseudotime/macrophage_pseudotime_protocol_sensitivity.tsv` — Phase 6.2 sensitivity output (sign flip).

### Figures
- `figures/main/fig2_progression_sex/panels/fig2_hotspot_cascade.pdf` — Manuscript Fig 2 cascade panel.
- `figures/main/fig2_progression_sex/panels/fig2e_embeddable.pdf` — MASH-vs-Healthy UMAP density.
- `figures/main/fig2_progression_sex/panels/fig2_panel_nmf_sc_celltype_composition.pdf` — Stacked-bar composition.
- `figures/supplementary/figS_hotspot_*` — Hotspot supp panels.
- `figures/supplementary/stage_ccc/*` — LIANA stage-CCC heatmaps + trajectories.

---

## 12. Environments

| Env | Used for |
|-----|----------|
| `rapids_singlecell` | GPU scRNA workflows (scVI, scANVI, Hotspot, cNMF). Init PyTorch CUDA before RMM pool; use 30% GPU mem; `LD_LIBRARY_PATH` fix. Shared utility: `Analysis/SingleCell/scripts/gpu_utils.py`. |
| `rnaseq` | Default for R/Python scripts (LMMs, LIANA aggregation, figure generation, donor metadata) |
| `spatial` | Python spatial / scanpy (rnaseq has scipy/numpy ABI conflict for Python) |
| `snapatac2` | scATAC integration (separate scipy/numpy stack) |
| `motifbreakr` | GWAS-ATAC motif disruption (Script 56) |
| `celltype_bio` | DIALOGUE multi-cellular program discovery |

---

## 13. SLURM partition guidance (for scRNA jobs)

- Hotspot discovery (501): `bigmem 256-512G 96h` for large cell types (global ≥ 800K cells, hepatocytes ~660K cells); `cpu 200G` or `io 96G` for smaller cell types (≤ 50K cells)
- 345b merge, 502/505/508/510, figure scripts: `io 32G` (small jobs)
- 346 stage LMM: `bigmem 500G` (over-allocated to bypass cpu queue; typical runtime < 2 min)
- 503 LOO stability: `gpu 80G` array (per dataset hold-out)
- Avoid `--qos=interactive` for routine parallel work (4-job cap)

Job naming convention: single word (e.g., `Hotspot`, `LMM`, `score`, `cascade`, `KDE`).

---

## 14. Outstanding rigor work (post-remediation)

These are tracked items not yet addressed by the 2026-05-22 remediation. Listed in priority order.

1. **Pseudobulk DE rewrite for cell-type DE**: Switch from single-cell Wilcoxon to muscat / Dream / edgeR pseudobulk for headline cell-type DE. Per Rich et al. *Cell Systems* 2026 (Pachter), pseudobulk is the modern standard.
2. **Version pinning**: Capture Scanpy / scVI / Hotspot / LIANA versions for every output. `Analysis/SingleCell/PROVENANCE.tsv` already drafted in archived `scrna_pipeline_rigor_plan.md`.
3. **Global + endothelial Hotspot re-discovery**: Deferred during remediation. Required if any supplementary figure shows cross-CT or endothelial-specific modules from the clean atlas.
4. **GSE244832 SRR pooling**: Sequencing-run replicates currently treated as independent donors. `donor_pairing.csv` documents 8 SRRs per biological donor for some D01-tier donors. Either pool SRRs in 345 or add `(1|biological_donor)` to the LMM.
5. **F-stage axis decoupling**: `F_stage_augmented` derives from anchors defined via `disease_stage_coarse` — partly tautological. The bootstrap check already routes around this for now; long-term, re-train on a scVI-independent criterion.
6. **LMM family-count correction**: `N_CT_PAIR_FAMILIES <- 72L` in `347b` does not match observed family counts on every axis (coarse=34, continuous=63, etc.). Replace with adaptive `n_families <- length(unique(out$ct_pair))`.
7. **Permutation null for 3-way Jaccard**: `353_three_way_concordance.R` should compute a permutation null to substantiate the "low / high concordance" framing.
8. **cNMF program expansion in fig2**: Pending user decision. Phase 6.3 verdict suggests P6 / P7 as best candidates for new monotonic-trajectory panels.

---

## 15. Index of related docs

| Doc | Status | Purpose |
|-----|--------|---------|
| `single_cell_analysis.md` (this file) | CURRENT | End-to-end consolidated reference |
| `reference/scrna_dataset_protocol_comparison.md` | CURRENT | Per-dataset protocol comparison table |
| `dataset_labeling_and_harmonization.md` | CURRENT | Bulk + scRNA dataset labels + diagnosis harmonization |
| `progression_atlas.md` | CURRENT | Bulk-focused progression atlas (Scripts 115–167 + 130–146) — has scRNA touchpoints |
| `sex_stratified_analysis_framework.md` | CURRENT | Mostly bulk; small scRNA section (per-cell-type sex × COLOC) |
| `figure_redesign/fig2_surrounding_context.md` | CURRENT | Fig 2 narrative context |
| `manuscript/working/02_results_fig2.md` | CURRENT (patched 2026-05-22) | Fig 2 manuscript text |
| `manuscript/working/02_results_fig2_audit_2026_05_22.md` | CURRENT | 14-item audit + Section D text patches |
| `manuscript/working/02_methods_scrna_protocol_audit.md` | CURRENT | Methods text + reviewer-defense paragraph for the remediation |
| `manuscript/methods_draft/agent3_atlas_drugs.md` | CURRENT | Long-form atlas / drugs methods draft |
| `manuscript/methods_draft/agent4_sc_epigenomic_spatial.md` | CURRENT | Long-form sc + epigenomic + spatial methods draft |
| `archive/single_cell_pre_remediation_2026-05-22/scrna_pipeline_rigor_plan.md` | ARCHIVED | Pre-remediation reviewer-defense plan (v2, 2026-05-13) |
| `archive/single_cell_pre_remediation_2026-05-22/ccc_review/*` | ARCHIVED | 10-file multi-team review of stage-CCC pipeline (2026-05-16/17) |
