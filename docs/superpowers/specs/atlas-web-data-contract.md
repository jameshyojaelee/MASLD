> **LEGACY JULY 2026 DATA CONTRACT — SUPERSEDED FOR PUBLICATION.** This rank-first
> schema is retained for implementation provenance only. It is not the current
> Resource portal authority and must not trigger a rebuild or deployment.
> Convergence ranks/tiers are compatibility annotations, not default ordering or
> target-validity fields. The replacement must implement the MASLD Gene Catalog,
> including testability, provenance, source dependence, and follow-up guidance,
> with adoption owned by Plan 60.

# MASLD Atlas Web — Historical Data Contract (Phase 0)

**Status:** historical implementation schema; superseded for publication on 2026-08-10. Written 2026-07-08 from then-current generator scripts and source-file headers.

**Project root** (`$ROOT`): `/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design`
**Web data dir** (`$WEB`): `masld-atlas-v2/public/data` (local) → mirrored to HF Dataset `jameshyojaelee/masld-atlas-data` at `resolve/main/`.

## §0. Global conventions (apply to every generator)

1. **Env:** all Python generators run under `micromamba run -n spatial` (rnaseq has a NumPy2/pandas ABI break). R steps (27a/27b) run under `rnaseq`.
2. **Project root** read from `os.environ.get("MASLD_PROJECT_ROOT", <default>)`. Every generator must honor `--output-dir`.
3. **JSON:** compact (`json.dump(obj, f, separators=(",",":"))`). NaN/±Inf → `None` (recursively). R-style `"TRUE"/"FALSE"` strings → Python bool. NEVER emit `NaN`/`Infinity` (invalid JSON, breaks the client parser).
4. **Parquet:** engine pyarrow, `index=False`, `±Inf → NaN` before write. Symbol-keyed long-format parquets MUST be sorted by `symbol` (enables DuckDB predicate pushdown / row-group skipping on `WHERE symbol=?`).
5. **Gene key:** `human_symbol` (aka `symbol`) is the canonical join key across all web files. Ensembl IDs are version-stripped (`ENSG…12` → `ENSG…`) via `strip_version()` when joining to per-study/pseudobulk/trajectory sources whose gene column is a versioned Ensembl ID.
6. **Doctrine (non-negotiable):**
   - **Never emit the string "TREAT"** anywhere (column names, JSON keys, labels, provenance). The columns `treat_lfc/treat_p/treat_fdr` are used numerically only; describe them as the **effect-size-aware interval-null FDR gate** (`fdr < 0.05` at `lfc = 0.25`, McCarthy & Smyth 2009). A CI grep-guard fails the build if the banned word appears in `masld-atlas-v2/src/`. Keep the raw column names as-is inside parquet (they are data, not UI), but **do not** invent a "TREAT" label anywhere in JSON/prose.
   - Atlas DE columns are `bulk_*` (the C2 canonical swap renamed `dream_* → bulk_*` on 2026-06-08). **Any web file that still carries `dream_*` is stale and MUST be renamed to `bulk_*` in its generator** (this bites `network_nodes.csv` and `drug_pipeline.json` — see §4, §6).
   - Human-only disease framing. Cross-**ancestry** (EUR/AFR/EAS/AMR/SAS) is first-class; cross-**species**/mouse content (`cross_species.json`, `/species`) is flagged keep/retire, not auto-removed (§8).
   - Control/healthy color `#9E9E9E`; significance encoded redundantly (never color alone). No lollipop, no 3D. (Frontend rule; noted here so chart-data shapes carry a `significant` bool where relevant.)

## §1. `atlas.parquet` (full atlas)

Generator: `masld-atlas-v2/scripts/preprocess_atlas_data.py` → `build_atlas_parquet()`.
**Rule:** pass-through of **every column** of `$ROOT/RNA-seq/results/multi_evidence/multi_evidence_atlas.csv` (currently **27,187 rows × 443 cols**, 2026-07-06), with `±Inf → NaN`. No renames, no drops. Row order preserved (atlas gene order).

- Preflight assertion already in code: `{"bulk_padj","bulk_logFC"} ⊆ columns` (fails if atlas still on retired `dream_*`). **Keep this assert.**
- The 443 columns are grouped as: identity (4), bulk DE (`bulk_*`, incl. `bulk_treat_fdr`), mouse (`mouse_*`, `translatability_score`, `is_conserved`), genetic/COLOC (`coloc_*`, `twas_*`, `susiex_*`, `mesusie_*`, `*_coloc_pp4` per-GWAS, `sceqtl_*`, `ieqtl_*`, `pleiotropy_*`), essentiality (`essentiality_chronos`, `is_essential`), epigenomic (`gwas_atac_*`, `gwas_variant_*`, `scenic_*`, `chromvar_top_tf`, `ccre_*`, `caqtl_*`, `abc_*`), spatial (`spatial_*`, `signature_govaere2026_*`), single-cell (`sc_*`, `liana_*`), sex (`sex_class`, `bulk_logFC_M/F`, `sex_interaction_padj`, `sex_progression_class`), programs (`prog_*`, `dominant_program_*`, `subtype_coloc_enriched`), zonation/ferroptosis (`zonation_class`, `spatial_zonation_class`, `ferroptosis_class`), drug (`dgidb_druggable`, `opentargets_drug`, `lincs_reversal`, `clintrial_*`, `drug_dev_status`, `max_phase_*`, `pharos_*`, `minikel_approval_or`), derived (`layers_active`, `sources_active`, `n_regulatory_layers`, `concordance`/`primary_category`).
- **Kept as a `/downloads` file only** (443 cols is too wide for client load). The client loads `atlas_core.parquet` (§2) instead.

Sibling JSON emitted by the same script (all in `$WEB`):

| File | Shape | Notes |
|---|---|---|
| `gene_index.json` | `[ {symbol, ensembl_id, biotype, bulk_logfc, bulk_padj, is_deg?, is_conserved?, sex_class?, zonation_class?, ferroptosis_class?, dgidb_druggable?, layers_active?, evidence:{s1_human..s8_proteomics}} ]` sorted by symbol | Compact; null/false optional keys omitted; `evidence` values rounded 2dp, zeros dropped. **Post-processed in place** by `generate_convergence_data.py` (renames s-keys per §4, injects `s5_epigenomic`+`s8_proteomics`). ⚠ Redesign retires this 10.6 MB file from the gene page (replaced by DuckDB on `atlas_core`); keep only a slim `gene_symbols.json` for fuse.js (see §8 item 4). |
| `atlas_summary.json` | `{total_genes, total_degs, total_cohorts, total_samples, conserved_count, coloc_genes, drug_targets, mouse_datasets, evidence_sources}` | ⚠ **Stale-number risk.** `compute_deg_count()` currently uses `padj<0.05 & |logFC|>0.3` (=legacy). Must switch to the canonical gate = `treat_fdr<0.05` (target **1,918**). `total_cohorts`/`total_samples` hardcoded (10/1444) — reconcile to headline. |
| `featured_genes.json` | `[ {symbol, ensembl_id, tagline, category, bulk_logfc, bulk_padj, is_deg, is_conserved, layers_active, evidence:{s1..s7}, drugs:[{drug,stage,moa,atlas_support}]} ]` | 9 curated symbols in `FEATURED_GENES` dict (THRB, PNPLA3, RORA, HSD17B13, PPARA, DGAT2, ACSL4, HKDC1, CFLAR). Taglines are editorial — refresh PP.H4 numbers. |

## §2. `atlas_core.parquet` (NEW — the client's primary table)

**New output of `preprocess_atlas_data.py`.** ~42-column subset of `atlas.parquet` **plus** three joined columns from `canonical_deg_results.csv` and four from `convergence_evidence.csv` (both keyed by `human_symbol`/`symbol`). This is what the frontend loads for the gene page header + all tables (DuckDB-WASM whole-file GET). Target < ~8–12 MB.

**Chosen columns (exact names):**

| Column | Source | Type | Meaning |
|---|---|---|---|
| `human_symbol` | atlas | str | gene symbol (join key) |
| `ensembl_id` | atlas | str | versioned Ensembl gene id |
| `gene_biotype` | atlas | str | GENCODE biotype |
| `mouse_ortholog` | atlas | str | mouse symbol (Cas13 lib only; not disease-framed) |
| `bulk_logFC` | atlas | float | pooled C2 limma-voom-qw logFC |
| `bulk_padj` | atlas | float | BH padj |
| `bulk_tstat` | atlas | float | moderated t |
| `bulk_shrunk_logFC` | atlas | float | ashr posterior logFC |
| `bulk_lfsr` | atlas | float | ashr local false sign rate |
| `bulk_sig` | atlas | bool | atlas-carried DE flag |
| `treat_lfc` | canonical_deg_results.csv | float | effect-size-aware interval-null logFC estimate (join on `symbol`) |
| `treat_p` | canonical_deg_results.csv | float | interval-null p-value |
| `treat_fdr` | canonical_deg_results.csv | float | interval-null FDR (the canonical DEG gate value) |
| `is_deg` | derived | bool | `treat_fdr < 0.05` (canonical; DEG count must equal **1,918**) |
| `bulk_logFC_M` | atlas | float | male-stratum logFC |
| `bulk_logFC_F` | atlas | float | female-stratum logFC |
| `sex_class` | atlas | str | {Concordant, Male_biased, Female_biased} |
| `sex_interaction_padj` | atlas | float | sex×group interaction padj |
| `coloc_best_susie_pp4` | atlas | float | best SuSiE-COLOC PP.H4 (canonical) |
| `coloc_best_susie_gwas` | atlas | str | GWAS giving the best SuSiE PP.H4 |
| `coloc_abf_best_pp4` | atlas | float | best ABF-COLOC PP.H4 (fallback) |
| `coloc_abf_best_gwas` | atlas | str | GWAS giving best ABF PP.H4 |
| `n_coloc_sources` | atlas | int | # colocalizing evidence sources |
| `n_ancestry_gwas` | atlas | int | # ancestries with a COLOC hit |
| `coloc_cross_ancestry_replicated` | atlas | bool | PP.H4 replicated across ≥2 ancestries |
| `coloc_susie_conf_tier` | atlas | str | High/Suggestive/Nominal confidence label |
| `twas_z` | atlas | float | OTTERS TWAS z |
| `twas_pval` | atlas | float | TWAS p |
| `convergence_rank` | convergence_evidence.csv | int | 1-based convergence rank (join on `human_symbol`) |
| `convergence_score` | convergence_evidence.csv | float | evidence-weighted convergence score (NOT a posterior) |
| `convergence_tier` | convergence_evidence.csv (`tier`) | str | {1_Genetic_validated, 2_Convergent, 3_Suggestive, 4_Weak}; Tier-1 target **677** |
| `concordance_state` | convergence_evidence.csv | str | {Concordant-up, Concordant-down, Protective-LOF, Conflicted, Insufficient-evidence} |
| `essentiality_chronos` | atlas | float | DepMap Chronos (more negative = more essential) |
| `is_essential` | atlas | bool | DepMap essential flag |
| `spatial_is_svg` | atlas | bool | spatially variable gene |
| `spatial_morans_i` | atlas | float | Moran's I |
| `spatial_consensus_direction` | atlas | str | consensus spatial direction |
| `zonation_class` | atlas | str | {Periportal, Pericentral, …} |
| `ferroptosis_class` | atlas | str | ferroptosis driver/suppressor class |
| `dgidb_druggable` | atlas | bool | DGIdb druggable |
| `opentargets_drug` | atlas | bool | OpenTargets has a drug |
| `max_phase_masld` | atlas | float | max clinical phase in a MASLD indication |
| `drug_dev_status` | atlas | str | e.g. Approved / Clinical / Preclinical / None |
| `pharos_tdl` | atlas | str | Tclin/Tchem/Tbio/Tdark |
| `dominant_program_for_gene` | atlas | str | dominant NMF program label (P1–P6) |
| `dominant_program_logFC` | atlas | float | logFC in the dominant program |
| `layers_active` | atlas | int | # active evidence layers |
| `is_conserved` | atlas | bool | conserved-core flag (mouse; retained for filtering, not disease-framed) |

**Join notes:** `treat_*` join = `atlas.human_symbol == canonical_deg_results.symbol` (left join; NaN where absent). `convergence_*` join = `atlas.human_symbol == convergence_evidence.human_symbol` (exclude `excluded_from_ranking` rows before joining). Preflight must assert `convergence_evidence.csv` mtime ≥ atlas mtime.

## §3. Gene detail parquets (NEW — replace the 33,943 per-gene JSONs)

Generator: **new** `scripts/portal/build_gene_detail_parquets.py` (replaces `generate_gene_profiles.py`; **reuse its loader functions** verbatim — they already resolve every source path and shim). Each file is **long-format, keyed by `symbol`, sorted by `symbol`**, DuckDB-filterable `WHERE symbol = ?`. All numeric columns rounded as the loaders already do.

### 3.1 `gene_per_cohort_de.parquet`
Source: `results/per_study/*_de_results.csv` (legacy single contrast) + `results/per_study_contrasts/{disease_vs_control,MASH_vs_control,MASL_vs_control,MASH_vs_MASL}/*_de.csv`. Excludes cohort `PRJNA512027` (per-cohort display only).

| col | type | meaning |
|---|---|---|
| `symbol` | str | gene |
| `contrast` | str | one of {`legacy_dream`, `disease_vs_control`, `MASH_vs_control`, `MASL_vs_control`, `MASH_vs_MASL`} |
| `dataset` | str | cohort GSE/PRJNA id |
| `logfc` | float | per-cohort logFC (3dp) |
| `pval` | float | P.Value (4dp) |
| `padj` | float | BH padj (4dp; null for `legacy_dream`) |

### 3.2 `gene_trajectories.parquet`
Source: staging_classifier trajectory CSVs + shims (`one_vs_rest_fibrosis_dream.csv`, `fibrosis_stage_trajectory_vs_F0.csv`, `…_vs_healthy.csv`, `nas_trajectory*.csv`). One row per (symbol, axis, reference, stage).

| col | type | meaning |
|---|---|---|
| `symbol` | str | gene |
| `axis` | str | {`fibrosis`, `nas`} |
| `reference` | str | {`vs_rest`, `vs_F0`, `vs_healthy`, `vs_nas0`} |
| `stage` | str | stage label (`F0`..`F4` or `NAS0`..) |
| `logfc` | float | 3dp |
| `padj` | float | 4dp |
| `n_cohorts` | int | # cohorts contributing (null where not a meta shim) |

### 3.3 `gene_coloc_by_gwas.parquet`
Source: `GWAS/finemapping/results/susie_coloc/susie_coloc_all_gwas.csv` (927,225 data rows — read only needed cols), filtered `PP.H4.abf > 0.01`.

| col | type | meaning |
|---|---|---|
| `symbol` | str | gene (`gene` col in source) |
| `gwas` | str | `gwas_name` |
| `ancestry` | str | `ancestry` (EUR/AFR/EAS/AMR/SAS) |
| `method` | str | `method` (susie/abf) |
| `pp4` | float | `PP.H4.abf` (4dp) |
| `pp4_susie` | float | `PP.H4.susie` (4dp, nullable) |

### 3.4 `gene_pseudobulk_de.parquet`
Source: `Analysis/SingleCell/results_gpu_v2/pseudobulk_de/*_de.csv` (skip `archive*`), filtered `padj < 0.1`. Gene id column is mixed Ensembl/symbol — resolve to `symbol` via atlas map (Ensembl-keyed rows joined through `ensembl_id`).

| col | type | meaning |
|---|---|---|
| `symbol` | str | gene |
| `cell_type` | str | pseudobulk cell type |
| `logfc` | float | 3dp |
| `padj` | float | 4dp |

### 3.5 `gene_lincs.parquet`
Source: `RNA-seq/results/drug_repurposing/lincs_final_ranked.csv`. One row per (symbol-target, compound); targets parsed from `target.x`/`dgidb_targets` (`;`/`,` split).

| col | type | meaning |
|---|---|---|
| `symbol` | str | target gene |
| `compound` | str | `display_name`\|`cmap_name`\|`pert_iname` |
| `score` | float | `composite_score` (3dp) |
| `moa` | str | `moa.x`\|`MOAss` (nullable) |

### 3.6 `gene_drugs.parquet`
Source: `RNA-seq/results/drug_repurposing/clinical_drug_validation_table.csv` (+ `stratified_causal/drug_target_progression_classification.csv` for `target_class`).

| col | type | meaning |
|---|---|---|
| `symbol` | str | `target_gene` |
| `drug` | str | drug name |
| `stage` | str | clinical stage / status |
| `moa` | str | mechanism of action |
| `support` | str | `atlas_support` |
| `target_class` | str | progression class (nullable; from drug_target_progression_classification) |

## §4. Existing web JSON shapes (top-level schema per file)

All in `$WEB`. Generators as noted. These are the CURRENT shapes (redesign keeps the shapes; refreshes the numbers).

**`convergence_ranking.json`** — `scripts/portal/generate_convergence_evidence_json.py` from `convergence_evidence.csv`. Canonical 46d ranking.
```
{ genes:[ { rank:int, symbol:str, convergence_score:float, log_convergence_odds:float,
            evidence_total:float, tier:str, concordance_state:str, dominant_stage_S1:str|null,
            druggability_tier:str, coloc_best_gwas:str|null, coloc_best_pp4:float|null,
            n_modalities_active:int, percentile:float, is_deg:bool, is_coloc:bool, is_druggable:bool,
            govaere2026_evidence?:bool, govaere2026_n_evidence?:int,
            convergence_score_with_govaere2026?:float, convergence_rank_with_govaere2026?:int } ],
  metadata:{ n_genes, source, score_definition, tier_definitions{}, concordance_states[] } }
```
(No `cross_species_concordant` — removed 2026-07-05.)

**`convergence_matrix.json`** — `generate_convergence_data.py`. Per-gene 8-modality vector.
```
[ { gene:str, modalities:[8×float 0-1] (s1_human,s2_mouse,s3_genetic,s4_essential,s5_epigenomic,
    s6_spatial,s7_singlecell,s8_proteomics), count:int(>0.1), is_deg:bool, is_coloc:bool, is_druggable:bool } ]
```
sorted by count desc. ⚠ `is_deg`/M1 currently use `padj<0.05 & |logFC|>0.3` (legacy) — reconcile to canonical gate.

**`bayesian_ranking.json`** — `generate_convergence_data.py` (legacy 46b archetype similarity; **NOT** the canonical ranking).
```
{ genes:[ { rank:int, symbol:str, score:float, percentile:float, is_deg:bool, is_coloc:bool, is_druggable:bool } ],
  metadata:{ n_genes, score_source:str, notes:str } }
```
⚠ 13 MB; keyed off atlas fallback formula unless `bayesian_posterior.csv` exists. Redesign may drop this file entirely (superseded by convergence_ranking).

**`proteomics_summary.json`** — `scripts/portal/generate_proteomics_data.py`.
```
{ n_proteins_tested:int, n_proteins_significant:int, n_atlas_tested:int, n_atlas_significant:int,
  overall_rho:float, conserved_rho:float, n_validated_targets:int,
  concordance_strata:[ {stratum:str, n:int, rho:float, direction_pct:float} ],
  top_upregulated:[ {gene,logfc,padj,tstat,dataset} ], top_downregulated:[ …same… ],
  datasets:[str], provenance:{pxd052937,pxd051911,gse276114,olink} }
```

**`progression_journey.json`** — `generate_phase3_data.py`.
```
{ stages:[ { stage:"F0".."F4", label:str, n_samples:int, n_degs_up:int, n_degs_down:int,
             top_genes_up:[{symbol,logfc,padj}], top_genes_down:[…], s1_fraction:float|null,
             s2_fraction:float|null, key_pathways:[str] } ],
  f2_switch:{ total_late_degs, early_degs, new_in_late, new_inflammatory_genes, key_switch_genes:[{symbol,logfc_early,logfc_late}] },
  classifier_metrics:{ f_ge_3_auroc, nas_ge_5_auroc, fibrosis_qwk } }
```
⚠ "F2 switch" narrative is RETIRED in the thesis — redesign should reframe as multi-step cascade / de-emphasize `f2_switch`. NMF fractions use legacy `nmf_subtype` S1/S2 (replace with P1–P6 axes, §5 programs).

**`drug_pipeline.json`** — `generate_phase3_data.py`.
```
{ funnel:[ {stage:str, count:int, description:str} ],
  top_compounds:[ {name,target,moa,composite_score,score_reversal,score_network,score_dgidb} ],
  validated_drugs:[ {drug,target,stage,moa,support,is_deg,bulk_logfc} ],
  sex_stats:{female_biased,male_biased,balanced,total_compounds} }
```
⚠ `validated_drugs[].bulk_logfc` is read from source col `dream_logFC` (non-atlas table). Redesign: retire the funnel (plan B-drugs → stacked bars by phase; **2 approved / 208 clinical / 1,504 preclinical**).

**`cross_species.json`** — `generate_phase3_data.py`. ⚑ **Human-only-doctrine review (§8 item 1).**
```
{ genes:[ {symbol,human_logfc,mouse_logfc,category,n_concordant,n_diets_sig,translatability_score,diets_concordant} ],
  summary:{total_orthologs,conserved,human_specific,mouse_specific,moderate_concordance,diet_selective,discordant,not_significant,diet_models[]},
  pathway_concordance:[ {pathway,pathway_id,human_nes,mouse_nes,concordant} ] }
```

**`gwas_atac_browser.json`** — `generate_phase4_data.py`.
```
{ variants:[ {id,chr,pos,pip,rec_pip,cell_types:[str],nearest_gene,distance,
              motifs_disrupted:[{tf,effect,in_regulon,alleleDiff?}], coloc_pp4?,coloc_gwas?} ] (top 200 by PIP),
  summary:{total_variants,in_peaks,with_motif_disruption,disease_regulon_tfs},
  cell_type_enrichment:[ {cell_type,count,enrichment,padj} ] }
```

**`knowledge_graph.json`** — `generate_phase4_data.py`.
```
{ nodes:[ {id,type:"gene"|"drug"|"tf"|"pathway",label,…type-specific…} ],
  edges:[ {source,target,type:"drug_target"|"regulon"|"pathway_member"} ],
  stats:{n_nodes,n_edges,n_genes,n_drugs,n_tfs,n_pathways} }
```
gene node: `{evidence_score(=layers_active),is_deg,is_coloc,logFC}`; drug: `{stage,moa,atlas_support}`|`{stage:"LINCS",composite_score}`; tf: `{n_targets,activity_diff}`; pathway: `{count}`.

**`pathway_genesets.json`** — `generate_phase4_data.py`.
```
{ collections:[ {name:"Hallmark", gene_sets:[ {id,name,genes:[str],size} ]} ],
  total_sets:int, total_genes_in_sets:int, universe_size:int }
```
⚠ `universe_size` hardcoded 33943 — reconcile to atlas row count (27,187).

## §5. New-layer outputs (5 summary JSON + parquets)

New generators in `scripts/portal/` (env `spatial`, `--output-dir`, compact JSON, NaN→None).

### 5.1 Single-cell — `generate_singlecell_data.py`
Sources (all under `Analysis/SingleCell/results_gpu_v2/`):
- `cell_type_proportions.csv` — WIDE: `sample,dataset,condition,preparation_method,condition_harmonized,<16 cell-type fraction columns>` (275 samples). Cell types: B cells, Basophils, Cholangiocytes, Circulating NK/NKT, Endothelial cells, Fibroblasts, Hepatocytes, Macrophages, Mono+mono derived cells, Neutrophils, Plasma cells, Resident NK, T cells, cDC1s, cDC2s, pDCs.
- `disease_signatures/celltype_disease_activity.csv` — gene × `ct_<CellType>_MASLD` activity matrix (first col unnamed = gene; 11 cell-type columns; 16 genes = signature summary rows).
- `disease_signatures/fgsea_celltype_F_transitions.csv` — `cell_type,contrast,NES,pval,padj,padj_per_contrast,size,leadingEdge_n,leadingEdge_genes` (120 rows).
- `pseudobulk_de/allcell_pseudobulk_de.csv` — `gene,lfc,ave,t_stat,pvalue,padj,B` (all-cell pseudobulk; 19,510 rows).
- `hepatocyte_subtypes/subtype_markers.csv` — `names,scores,logfoldchanges,pvals,pvals_adj,pct_nz_group,pct_nz_reference,subtype`.
- `hepatocyte_subtypes/subtype_axis_scores.csv` — `subtype,axis_periportal,axis_pericentral,axis_lipid_accumulation,axis_inflammatory,axis_ferroptosis`.
- `hepatocyte_subtypes/subtype_disease_enrichment.csv` — `subtype,n_cells,n_masld,n_healthy,odds_ratio,fisher_pval,fisher_padj,enrichment_class,stage_spearman_rho,stage_spearman_p,stage_kw_stat,stage_kw_pval`.

**Output `singlecell_summary.json`:**
```
{ n_cells:int, n_datasets:int, n_samples:int, cell_types:[str],
  composition:[ { sample:str, dataset:str, condition:str, condition_harmonized:str,
                  fractions:{<cell_type>:float} } ],
  celltype_disease_activity:[ { cell_type:str, activity:float } ],   # from celltype_disease_activity summary
  fgsea_transitions:[ { cell_type:str, contrast:str, nes:float, padj:float, size:int,
                        leading_edge_n:int, significant:bool } ],
  hep_subtypes:[ { subtype:str, n_cells:int, odds_ratio:float, fisher_padj:float,
                   enrichment_class:str, stage_spearman_rho:float,
                   axes:{periportal,pericentral,lipid_accumulation,inflammatory,ferroptosis} } ] }
```
**Output `sc_pseudobulk_de.parquet`** (long, sorted by symbol) — this is the ALL-CELL table for the `/single-cell` DE tab; per-cell-type detail already lives in `gene_pseudobulk_de.parquet` (§3.4).

| col | type | meaning |
|---|---|---|
| `symbol` | str | gene (`gene` col) |
| `logfc` | float | `lfc` |
| `ave_expr` | float | `ave` |
| `tstat` | float | `t_stat` |
| `pval` | float | `pvalue` |
| `padj` | float | `padj` |

Hepatocyte marker dotplot data (`subtype_markers.csv`) → optionally a small `sc_hep_markers.parquet` `{symbol, subtype, logfoldchanges, pvals_adj, pct_nz_group}` (top-N per subtype) if the dotplot needs it client-side.

### 5.2 Spatial — `generate_spatial_data.py`
Sources (under `Analysis/Spatial/results/`):
- `zonation/zonation_scores.csv` — per-sample: `,sample_id,dataset,condition,zonation_score,zonation_bin,periportal_score,pericentral_score` (6,546 rows).
- `zonation/deg_zonation_classification.csv` — `,gene,zonation_class,kruskal_stat,kruskal_pval,spearman_rho,spearman_pval,n_donors,mean_PP1,mean_PP2,mean_Mid,mean_PC2,mean_PC1` (6,136 rows).
- `zonation/zonation_disruption_scores.csv` — `,gene,sample,rho_healthy,rho_masld,disruption_score,n_healthy_donors,underpowered` (42 rows).
- `svg/differential_svgs.csv` — first col unnamed = gene; `morans_I_healthy,morans_I_masld,svg_healthy,svg_masld,delta_I,category` (3,000 rows).

**Output `spatial_summary.json`:**
```
{ n_samples:int, datasets:[str],
  zonation_dist:[ { zonation_bin:str, n:int } ],       # counts per bin
  zonation_by_condition:[ { condition:str, mean_zonation_score:float, n:int } ],
  top_disrupted:[ { gene:str, rho_healthy:float, rho_masld:float, disruption_score:float, underpowered:bool } ],
  svg_categories:[ { category:str, n:int } ],
  top_differential_svgs:[ { gene:str, morans_i_healthy:float, morans_i_masld:float, delta_i:float, category:str } ] }
```
**Output `spatial_zonation.parquet`** (long, sorted by symbol) — the zonation heatmap / gene-lookup table:

| col | type | meaning |
|---|---|---|
| `symbol` | str | gene |
| `zonation_class` | str | Periportal/Pericentral/Mid/… |
| `spearman_rho` | float | zonation-axis correlation |
| `kruskal_pval` | float | KW p across zones |
| `n_donors` | int | donors supporting |
| `mean_PP1`/`mean_PP2`/`mean_Mid`/`mean_PC2`/`mean_PC1` | float | mean expr per zone bin (5 cols) |
| `morans_i_healthy` | float | from differential_svgs (join by gene; nullable) |
| `morans_i_masld` | float | from differential_svgs (nullable) |
| `delta_i` | float | Δ Moran's I (nullable) |
| `svg_category` | str | differential SVG category (nullable) |

### 5.3 Ancestry-COLOC — `generate_ancestry_coloc_data.py`
Sources (under `GWAS/finemapping/results/susie_coloc/`):
- `gene_level_coloc.csv` — 18,975 rows; cols incl. `gene,ensembl,coloc_best_pp4,coloc_best_gwas,coloc_n_gwas_h4_05,…,coloc_best_pp4_EUR,coloc_best_gwas_EUR,coloc_n_anc_h4_05,coloc_best_susie_pp4,coloc_best_susie_gwas,coloc_n_gwas_susie_h4_05/08/09,coloc_best_susie_pp4_EUR,coloc_n_anc_susie_h4_05,coloc_best_ancestry,coloc_best_susie_ancestry,coloc_abf_conf_tier,coloc_susie_conf_tier,coloc_abf_ancestry_matched,coloc_susie_ancestry_matched,coloc_abf_headline_cross_anc,coloc_susie_headline_cross_anc,coloc_cross_ancestry_replicated,…` (full 44-col header captured).
- `susie_coloc_all_gwas.csv` — 927,225 rows: `gwas_name,gene,ensembl,chr,PP.H0-4.abf,PP.H3/4.susie,n_cs_pairs,n_snps,method,top_snp,top_snp_PP,ancestry,ld_panel,ld_panel_n,…,lambda_s,ld_reliability`. **Read header only; filter PP.H4.abf > 0.01 when materializing.**
- `broadaway_coloc_comparison.csv` — 751 rows: `ensembl,gene,coloc_best_susie_pp4,coloc_best_pp4,our_best_pp4,our_drivers,via_EUR_ALT_GGT,via_EUR_AST,via_nonEUR,via_MASLD_disease,in_broadaway_747,broadaway_traits,novel_to_us,novel_source`.

**Output `coloc_ancestry_summary.json`:**
```
{ n_effector_genes_susie:int (target 473), n_effector_genes_union:int (target 1031),
  n_gwas:int (target ~50-55), ancestries:[EUR,AFR,EAS,AMR,SAS],
  by_ancestry:[ { ancestry:str, n_gwas:int, n_genes_h4_05:int, n_genes_h4_08:int } ],
  conf_tier_counts:[ { tier:"High"|"Suggestive"|"Nominal", n:int } ],
  cross_ancestry_replicated:[ { gene:str, best_pp4:float, ancestries:[str] } ],
  broadaway_overlap:{ n_in_broadaway_747:int, n_novel_to_us:int } }
```
**Output `coloc_by_gwas.parquet`** (long, sorted by symbol) — from `susie_coloc_all_gwas.csv` filtered `PP.H4.abf > 0.01`:

| col | type | meaning |
|---|---|---|
| `symbol` | str | gene |
| `ensembl` | str | ensembl id |
| `gwas` | str | `gwas_name` |
| `ancestry` | str | EUR/AFR/EAS/AMR/SAS |
| `method` | str | susie/abf |
| `pp4` | float | `PP.H4.abf` |
| `pp4_susie` | float | `PP.H4.susie` (nullable) |
| `ld_panel` | str | LD reference (polyfun/1kg/…) |
| `ld_reliability` | str | reliability flag |

(This is the same source as §3.3 `gene_coloc_by_gwas.parquet`. **Recommendation:** materialize ONCE as `coloc_by_gwas.parquet` and have the gene page query it too, dropping §3.3 as a separate file — noted in §8 item 5.)

### 5.4 Programs — `generate_programs_data.py`
Sources:
- `RNA-seq/results/subtypes/program_labels.csv` — `program_code,biological_label,label_category,top_pathway,top_genes` (6 rows = P1–P6). Labels: P1 Inflammatory-EMT, P2 Hepatocyte-Metabolic, P3 Fibrotic-ECM, P4 Unresolved, P5 Noncoding, P6 Skeletal-muscle.
- `RNA-seq/results/subtypes/nmf_assignments.csv` — `sample_id,dominant_program_code,dominant_program,nmf_subtype,P1..P6,P_<name>_score…` (1,104 samples).
- `RNA-seq/results/subtypes/subtype_markers.csv` — `gene,program_code,program,logFC,AveExpr,t,pvalue,padj,direction,subtype,adj.P.Val` (165,828 rows).
- `Analysis/SingleCell/results_gpu_v2/hotspot_modules/module_names.tsv` — `cell_type,module,module_hallmark,module_top_pathway,module_name` (117 modules).
- `Analysis/SingleCell/results_gpu_v2/hotspot_modules/all_modules.tsv` — 117 rows, 39 cols incl. `cell_type,module,best_match_*,matched,is_novel,stability_score,loo_*,disease_stage_beta,F_stage_beta,NAS_beta,*_q,progression_module,activity_module,bulk_module_score,n_genes_overlapping,mean_bulk_logFC,pct_concordant_up,bulk_replicated,module_hallmark,module_top_pathway,module_name`.
- `Analysis/SingleCell/results_gpu_v2/mcp/crn_changepoints_clean.csv` — `program,best_model,aic_linear,aic_gam,aic_seg1,delta_aic,breakpoint_stage,breakpoint_ci_low,breakpoint_ci_high,n` (16 cnmf_global_k16 programs).

**Framing:** k=6 NMF programs = **continuous interpretive axes P1–P6**, NOT discrete subtypes (2026-07-04 audit). Do not present as patient stratification.

**Output `programs_summary.json`:**
```
{ nmf_programs:[ { code:"P1".."P6", label:str, label_category:str, top_pathway:str, top_genes:[str],
                   mean_activity_by_stage:[ {stage:"F0".."F4", mean_score:float} ] } ],   # from nmf_assignments × metadata fibrosis_stage
  hotspot_modules:[ { cell_type:str, module:str, module_name:str, module_hallmark:str,
                      is_novel:bool, stability_score:float, bulk_replicated:bool,
                      n_genes:int, mean_bulk_logfc:float, progression_module:bool } ],
  n_hotspot_modules:int, n_novel:int,
  changepoints:[ { program:str, best_model:str, breakpoint_stage:float|null, delta_aic:float } ] }
```
**Output `program_gene_membership.parquet`** (long, sorted by symbol):

| col | type | meaning |
|---|---|---|
| `symbol` | str | gene |
| `program_code` | str | P1–P6 |
| `program_label` | str | biological label |
| `logfc` | float | marker logFC in that program |
| `padj` | float | `adj.P.Val` |
| `direction` | str | up/down |

(Hotspot module→gene membership, if surfaced per-gene, is a second table `hotspot_module_membership.parquet` `{symbol, cell_type, module, module_name, is_novel}` from the hotspot module gene lists — optional; flag in §8 if the module gene-list file must be located.)

### 5.5 Sex — `generate_sex_data.py`
**Source (PINNED — resolves Open Item #2):**
`$ROOT/RNA-seq/Human/Patient_Cohorts/analysis/integration/results/integration/sex_v3/lvqw_staging/sex_deg_classification_v3.csv`, column **`sex_class`** ∈ {`Male_biased`, `Female_biased`} → **exactly 7 Male + 1 Female = the canonical 8** (verified). Effect columns in the same file: `logFC_M`, `logFC_F`, `interaction_padj`, `beta_M`, `beta_F`, `gene_symbol`, `ensembl_base`, `chr`, `gene_biotype`.

⚠ **Do NOT use** the top-level `sex_v3/sex_deg_classification_v3.csv` (Jun 9, overwritten with a v6-consensus that gives 11F/4M) nor the atlas `sex_class` (maps only 6 = 5M/1F after dedup) nor `male_biased_characterization.csv` (5 rows, characterization not the canonical set). The `lvqw_staging/` copy is the frozen LVQW-collapse output that yields the headline 8.

**Output `sex_summary.json`:**
```
{ n_dimorphic:8, n_female_biased:1, n_male_biased:7,
  genes:[ { symbol:str, ensembl:str, chr:str, biotype:str, sex_class:"Male_biased"|"Female_biased",
            logfc_female:float, logfc_male:float, interaction_padj:float,
            beta_female:float, beta_male:float } ],
  note:"LVQW-fixed (C2) sex layer; female-progressor narrative retired 2026-06-08. Sex-biased framing only." }
```
**Output `sex_gene_classification.parquet`** — same 8 rows, columns `{symbol, ensembl, chr, biotype, sex_class, logfc_female, logfc_male, interaction_padj, beta_female, beta_male}`. (Full 27,638-row classification is NOT web-shipped; only the 8 dimorphic.)

## §6. Network parquets

Generator: **new** `scripts/portal/network_to_parquet.py`. Sources of truth = `RNA-seq/results/network/edge_annotation_atlas.parquet` (**15,295,605 rows × 65 cols**, 1.5 GB) + `RNA-seq/results/network/network_nodes.csv`. **The 16,601 per-gene-graph JSONs are retired.**

⚠ **The full edge parquet is 15.3M edges — must be FILTERED for the web.** Recommended web filter: keep edges with `in_string_ge700 == True` OR `coloc_linked == True` OR `|delta_max| ≥ threshold` OR `is_contested`, AND both endpoints in a bounded node set (e.g. top-N by `layers_active`, or the DE/COLOC/druggable union). Document the exact filter in the generator; target < ~5–10 MB.

### 6.1 `network_edges.parquet` (filtered subset)
Columns to carry (subset of the 65; keep the interpretable + provenance ones):

| col | type | meaning |
|---|---|---|
| `gene_a`, `gene_b` | str | edge endpoints (symbols) |
| `type` | str | edge channel (coexpr/coloc/lr/regulon/cerna/xspecies) |
| `in_string_ge700` | bool | STRING support ≥700 |
| `is_contested` | bool | contested edge |
| `string_score` | float | STRING combined score |
| `r_F01`,`r_F2`,`r_F34` | float | stage-binned correlation |
| `delta_max` | float | max Δ correlation across stages |
| `emergence_stage` | str | stage the edge emerges |
| `fdr_f2` | float | differential-edge FDR |
| `loco_replication_fraction` | float | LOCO stability |
| `pp4_min` | float | min COLOC PP4 across endpoints |
| `gwas_list` | str | GWAS supporting the coloc link |
| `score_diff_lr` | float | ligand–receptor differential score |
| `cell_type_pairs` | str | CCC cell-type pair(s) |
| `reg_tf`,`reg_target` | str | regulon TF→target |
| `regulon_activity_diff` | float | regulon activity Δ |
| `coloc_linked` | bool | endpoints share a coloc locus |
| `druggable_pair` | bool | both endpoints druggable |
| `community_macro_a`/`_b` | int/str | macro community assignment (for layout/coloring) |

### 6.2 `network_nodes.parquet`
Source `network_nodes.csv` cols: `human_symbol,ensembl_id,gene_biotype,dream_logFC,dream_padj,dream_tstat,is_deg,coloc_susie_best_pp4,is_conserved_core,sex_class,zonation_class,ferroptosis_class,dgidb_druggable,n_coloc_sources,attribution_class,mouse_ortholog`.

⚠ **`dream_*` → `bulk_*` rename required** (this CSV predates C2). Output schema:

| col | type | meaning |
|---|---|---|
| `symbol` | str | `human_symbol` |
| `ensembl_id` | str | ensembl id |
| `gene_biotype` | str | biotype |
| `bulk_logFC` | float | ← `dream_logFC` (RENAME) |
| `bulk_padj` | float | ← `dream_padj` (RENAME) |
| `bulk_tstat` | float | ← `dream_tstat` (RENAME) |
| `is_deg` | bool | DE flag |
| `coloc_susie_best_pp4` | float | best SuSiE PP4 |
| `is_conserved_core` | bool | conserved-core |
| `sex_class` | str | sex class |
| `zonation_class` | str | zonation |
| `ferroptosis_class` | str | ferroptosis |
| `dgidb_druggable` | bool | druggable |
| `n_coloc_sources` | int | # coloc sources |
| `attribution_class` | str | deconvolution attribution |
| `mouse_ortholog` | str | mouse symbol |

⚠ Ideally regenerate `network_nodes.csv` from the current atlas so `is_deg` uses the canonical gate and DE values are truly `bulk_*`; short-term the rename above is acceptable but flag as staleness risk (§8).

## §7. Page → data map

| Route | Data files consumed |
|---|---|
| `/atlas` | `atlas_core.parquet` (DuckDB table + volcano; threshold `|logFC|` guide at 0.25), `atlas_summary.json`, `featured_genes.json`, `gene_symbols.json` (fuse.js) |
| `/gene` (`?symbol=`) | `atlas_core.parquet` (header row `WHERE human_symbol=?`), `gene_per_cohort_de.parquet`, `gene_trajectories.parquet`, `coloc_by_gwas.parquet` (§5.3, replaces §3.3), `gene_pseudobulk_de.parquet`, `gene_lincs.parquet`, `gene_drugs.parquet`, `spatial_zonation.parquet`, `sex_gene_classification.parquet`, `program_gene_membership.parquet` — all `WHERE symbol=?` |
| `/single-cell` (absorbs `/resolution`) | `singlecell_summary.json`, `sc_pseudobulk_de.parquet`, `gene_pseudobulk_de.parquet` (per-cell-type), `sc_hep_markers.parquet` (dotplot), UMAP embedding (⚠ §8 item 6), `spatial_summary.json`+`spatial_zonation.parquet` (spatial tab) |
| `/genetics` (absorbs `/causal` + GWAS-ATAC) | `coloc_ancestry_summary.json`, `coloc_by_gwas.parquet`, `atlas_core.parquet` (effector-gene table: `coloc_best_susie_pp4`), `gwas_atac_browser.json` (regulatory sub-panel) |
| `/programs` (NEW) | `programs_summary.json`, `program_gene_membership.parquet`, `sex_summary.json`+`sex_gene_classification.parquet` (sex-biased diverging bars) |
| `/drugs` | `drug_pipeline.json` (stacked-by-phase; retire funnel), `gene_drugs.parquet`, `gene_lincs.parquet`, `atlas_core.parquet` (druggability cols) |
| `/progression` | `progression_journey.json`, `gene_trajectories.parquet`, `atlas_core.parquet` |
| `/translation` | `atlas_core.parquet` (bulk↔spatial↔proteomics corroboration), `proteomics_summary.json`, `spatial_summary.json` |
| `/proteomics` | `proteomics_summary.json` |
| `/network` | `network_edges.parquet`, `network_nodes.parquet` |
| `/causal` | (folded into `/genetics`) — same files |
| `/species` (⚑ review) | `cross_species.json` — keep only if human-only doctrine decision is "relabel-and-keep" (§8 item 1) |
| `/downloads` | `atlas.parquet` (full 443-col), `convergence_ranking.json`, plus a manifest of every parquet/JSON |
| `/explore` | `atlas_core.parquet`, `convergence_ranking.json`, `convergence_matrix.json`, `pathway_genesets.json`, `gene_symbols.json` |

## §8. Open / uncertain items

1. **`/species` + `cross_species.json`** — human-only doctrine leans retire/relabel. Data shape is defined (§4) but keep/retire is a product decision; not blocking. Cross-**ancestry** stays.
2. **Sex-8 source — RESOLVED.** `sex_v3/lvqw_staging/sex_deg_classification_v3.csv`, `sex_class` col = 7 Male + 1 Female (verified). Use this file, NOT the top-level v3 (11F/4M) or the atlas (6). No longer blocking.
3. **`atlas_summary.json` numbers are stale** (`total_degs` via legacy `padj<0.05,|logFC|>0.3`; hardcoded cohorts/samples; `pathway_genesets.universe_size=33943`). Backend must re-point `compute_deg_count()` to `treat_fdr<0.05` (=1,918) and reconcile counts to the atlas (27,187 genes) and headline sample/cohort totals. Purely a generator patch.
4. **`gene_index.json` (10.6 MB) retirement** — redesign drops it from the gene page (DuckDB on `atlas_core` supplies header + evidence fingerprint derived from atlas columns). A slim `gene_symbols.json` (`[str]` or `[{symbol,biotype}]`) must be emitted for fuse.js autocomplete. Generator for `gene_symbols.json` is unspecified — assign to `preprocess_atlas_data.py`.
5. **COLOC parquet duplication** — §3.3 `gene_coloc_by_gwas.parquet` and §5.3 `coloc_by_gwas.parquet` have the same source+filter. **Recommend materializing once** as `coloc_by_gwas.parquet` (in `generate_ancestry_coloc_data.py`) and dropping the §3.3 build from `build_gene_detail_parquets.py`. Confirm with Team A before both write it.
6. **UMAP embedding for `/single-cell`** — none of the listed CSVs carry per-cell 2-D coordinates. The regl-scatterplot UMAP needs a downsampled embedding (~50–75k cells) exported from the integrated `.h5ad` (env `rapids_singlecell`/`spatial`). Proposed output `sc_umap_downsampled.parquet` `{cell_id, umap_x, umap_y, cell_type, dataset, condition_harmonized}`. Source h5ad path + downsample step is UNSPECIFIED here — Team A/C must locate the atlas embedding. Flagged, non-blocking for the other 4 new layers.
7. **Network web filter threshold** — 15.3M edges must be reduced (§6). The exact predicate (STRING≥700 ∪ coloc_linked ∪ |delta_max|≥τ, bounded node set) is a generator decision; pick τ so the output is < ~10 MB and record it in the generator + dataset README.
8. **`hotspot_module_membership.parquet`** (§5.4) — per-gene hotspot module membership needs the module→gene-list file (not `all_modules.tsv`, which is module-level). Locate the gene-membership TSV under `hotspot_modules/` before building; optional for v1 of `/programs`.
9. **`network_nodes.csv` provenance** — carries retired `dream_*`; §6.2 renames them, but ideally regenerate from the current atlas so `is_deg`/DE values are canonical. Staleness risk if only renamed.
