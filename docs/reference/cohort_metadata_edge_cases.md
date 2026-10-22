# Human RNA-seq cohort metadata edge cases (reference)

**Created 2026-06-24.** Durable catalog of cohort-specific metadata quirks that can contaminate stage/NAS/sex-stratified analyses, and how the codebase guards against them. Read this before building any analysis or figure that bins per-sample `fibrosis_stage`, `nas_score`, or `sex`.

## The guardrail (single source of truth)
`scripts/figures/load_figure_data.R` defines:
- Constants `COARSE_STAGE_COHORTS = c("GSE213621")`, `DROPPED_COHORTS = c("PRJNA512027")`, `INFERRED_SEX_COHORTS = c("GSE135251","GSE213621","GSE240729")`, `STAGED_FIB_COHORTS` (the 6 true-Kleiner cohorts = the Script 14b DE allowlist).
- `staged_disease_meta(disease_only=TRUE)` → per-sample metadata for the clean true-Kleiner staged set (excludes coarse + dropped cohorts) with a ready `stage` factor (F0–F4). Use `disease_only=FALSE` to keep controls as F0.
- `clean_nas_meta()` → NAS-stratified samples with non-NA NAS (drops Bril's unscored patients + Chen's wholly-absent NAS) and dropped cohorts.

Any figure/analysis that reads raw `unified_metadata.csv` and bins `fibrosis_stage`/`nas_score` MUST route through these (or apply `!dataset %in% c(COARSE_STAGE_COHORTS, DROPPED_COHORTS)`). The DEG/GSEA per-stage result tables (`fibrosis_stage_dream.csv`, `fibrosis_stage_vs_ctrl_dream.csv`, `nas_score_dream.csv`; Scripts 14b/14c) are **already clean** — they use the 6-cohort allowlist.

## Edge-case catalog

### CRITICAL — GSE213621 (Chen): coarse grouped fibrosis
The study reports grouped bins (`Fibrosis_F0F1`, `F2`, `F3F4`), and `00_harmonize_metadata.R:161-166` maps them to integers `Control→0, F0F1→1, F2→2, F3F4→3`. So Chen's `fibrosis_stage` is **NOT true Kleiner F0–F4**: "stage 1" = an F0+F1 mix, "stage 3" = an F3+F4 mix; it has **no true F0 or F4**. Chen is the **largest cohort** — at the metadata level it is ~⅓ of all F1/F2/F3 disease samples — so binning it as if it had real stages materially distorts stage-stratified composition/trajectory readouts. **Guard: exclude from all per-sample stage binning** (the 14b DE already does). NAS is wholly absent for Chen.

### HIGH (NAS only) — GSE162694 (Bril): NAS partial among patients
All 112 disease patients are precisely Kleiner-staged (fibrosis = Available); only the 31 controls are NA. But **26 of the 112 disease patients lack a NAS score** (NAS = Partial). Fibrosis is fine to use; NAS-stratified analyses must drop the unscored patients (`clean_nas_meta()` does).

### Handled — PRJNA512027 (Gerhard): dropped batch confound
Dropped from paper presentation 2026-05-15 (L0/S0 library-prep batch perfectly confounded with diagnosis). `include_in_mega=false`; already absent from the mega-DE, the per-stage DE allowlists, and the NMF/deconvolution inputs. Guarded in `DROPPED_COHORTS` for any figure that still reads raw metadata.

### MED (labeling) — inferred sex: GSE135251, GSE213621, GSE240729
`sex_annotated: false` in `config/human_datasets.yaml`; sex is inferred by XIST/DDX3Y k-means (Script 01, seed 42). Used as a legitimate model covariate (`inferred_sex`); not a contamination, but figures should not present it as annotated. `fig1_umap.R` / `fig1_atlas_overview_v2.R` already label these honestly.

### MED (handled) — control-less cohorts: GSE167523, GSE174478, GSE193066, GSE240729
Disease-only (no healthy controls); `include_in_mega=false`. Valid for stage-stratified disease analyses (they carry true Kleiner stages), excluded from Disease-vs-Control contrasts by the mega allowlist.

## Clean per-stage sample sizes (after Chen + Gerhard exclusion)
`staged_disease_meta()` disease samples, metadata level: **F0 95 / F1 196 / F2 208 / F3 151 / F4 46** (6 cohorts). Figures intersect this with their NMF/deconv inputs (e.g. fig3k NMF: F0 95 / F1 194 / F2 205 / F3 147 / F4 45; celltype_proportion incl. controls = 719 donors / 6 cohorts). **F4 (~45) is the smallest bin** — an honest power limitation, since Chen contributed no true F4 (its F4 patients were hidden in the F3F4 group) and Gerhard's high-grade samples were batch-confounded.

## Sensitivity (fig3k NMF dominant-program %, with vs without Chen)
Modest but real; trajectory shape preserved. Fibrogenic-dominant % rises at mid-stages after excluding Chen's coarse samples: F1 5→7, F2 9→11, **F3 32→39**; F0 7→7, F4 60→60 unchanged. Progression-Inflammatory ≈stable (F1 7→6, others unchanged). Qualitative monotonic story is unchanged; mid-stage percentages shift ≤7 pp.

## Affected scripts fixed (2026-06-24)
- **Figures** (`scripts/figures/`): `nmf_dominant_program_stage_composition` (fig3k), `progression_cascade` (Tracks 2&3), `celltype_proportion`, `figS08_nmf_subtyping`, `nmf_programs`, `nmf_programs_lines`, `figS_pathway_programs`, `figS_pi_concordance_panels(_lfsr)`, `fig_composition_features`, `figS_composition_mediation`, `fig_disease_signatures_presentation`, `nas_fib_grid` (robustness), and 3 clinical PCAs (`figS_pca_clinical`, `figS_supervised_disease_pca_top200DEG_clinical`, `figS_covariate_pc_association` — Chen's coarse fibrosis nulled → shown as "Unknown", not mislabeled F1/F3).
- **Analysis**: `250c_nmf_program_by_stage.R` (Chen exclusion; output not figure-consumed → re-run deferred).
- **Already clean (no change)**: all DEG/GSEA per-stage panels (`fib_stage_degs`, `fib_stage_upset`, `fib_stage_vs_ctrl_*`, `fibrosis_stage_directional_asymmetry`, `fibrosis_stage_hallmark_escalation`, cascade Track 1); `fig1_umap`, `fig1_atlas_overview_v2` (honest coarse labels).
