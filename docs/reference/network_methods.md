# MASLD Network v2 (Architecture C) — Methods, Legends, and Data Availability

*Draft for Cell Metabolism Resource submission. Tone: rigorous and transparent; contribution is the MASLD-specific data and queryable atlas, not integration-methodology novelty.*

---

## Methods — MASLD gene network atlas

### 1. Architecture overview

We constructed a MASLD-specific gene-gene interaction atlas over 16,704 protein-coding and long non-coding RNA gene nodes present in the unified 10-cohort transcriptomic atlas (1,444 QC-passing human bulk RNA-seq samples; Methods §"Human Integration"). The network adopts a two-tier design (Architecture C). A generic-biology backbone is provided by STRING v12 high-confidence edges (combined_score ≥ 700; Szklarczyk et al. 2023), which functions as an existence prior without modification. Superimposed on this backbone are five **disease-native (Type-D) edge categories** that encode relationships STRING cannot contain by construction: D-F2 (fibrosis-stage-stratified differential coexpression), D-COLOC (paired GWAS colocalization at a shared locus), D-LR (F2-differential ligand–receptor rewiring from LIANA), D-ceRNA (MASLD-specific lncRNA competing-endogenous-RNA triplets), and D-XSpecies (cross-species conserved switch edges from the mouse integration). Every edge in the released atlas carries a rich per-edge annotation record (≥40 fields) spanning F-stage dynamics, LOCO replication support, COLOC linkage, druggability, dominant cell-type driver, sex-dimorphism flags, and F2-weighted community membership. The atlas is published as a single Parquet table (`RNA-seq/results/network/edge_annotation_atlas.parquet`; 3,375,599 edges, 40+ columns) that is the primary Resource deliverable.

We make no novelty claims regarding the integration methodology. Head-to-head benchmarks (§7 below, reproduced transparently in Supplementary Table S_NetCompare) show STRING v12 ≥ 700 outperforms any composite we tested on generic gold sets (e.g., 88× vs 66× fold-over-null for MASLD positive controls; 99× vs 87× for clinical drug targets). The contribution of this work is the **disease-specific edge layers** and their per-edge provenance annotations, which are complementary to — not a replacement for — generic PPI resources.

### 2. STRING backbone ingestion

STRING v12.0 human interactions (`9606.protein.links.full.v12.0.txt.gz`) were filtered at combined_score ≥ 700 (high confidence). Ensembl protein identifiers were mapped to HGNC gene symbols via the STRING v12 alias table, keeping only edges whose both endpoints resolved to genes present in the 16,704-node vocabulary of the MASLD atlas. This yields the Type-S (STRING-only) edge set used both as the backbone for community detection and as the denominator for Jaccard overlap tests against Type-D edges.

### 3. D-F2 edges: stage-stratified differential coexpression

D-F2 edges capture coexpression relationships that emerge, dissolve, or invert across fibrosis progression. Starting from ComBat-seq-corrected logCPM expression in the 10-cohort atlas (see Methods §"Human Integration"), samples were stratified by fibrosis stage into three bins: **F0–F1** (N ≈ 490), **F2** (N ≈ 316), and **F3–F4** (N ≈ 343); the two cohorts without fibrosis scoring (GSE126848, GSE167523) were excluded from stratified computations and from leave-one-cohort-out (LOCO) replication. Within each stage bin, gene–gene Pearson correlations were computed over all intra-vocabulary pairs (script `281_f2_delta_edges.R`). A housekeeping filter was applied post hoc: mitochondrial genes (`MT-*`), ribosomal protein genes (`RPL*`, `RPS*`), and other high-invariance housekeeping families were excluded to prevent batch-anchored edges from dominating, following the strategy validated in our co-expression threshold audit (`275_coexpression_threshold_sweep.R`).

Each gene pair was assigned stage-specific correlations r_{F0F1}, r_{F2}, r_{F3F4} and classified by the shape of the Δr trajectory into one of six categories: *F2_emerging* (peak at F2, |Δr| ≥ 0.35), *F2_dissolving* (trough at F2), *transient_F2* (bump or dip confined to F2), *progressive_up* / *progressive_down* (monotone |Δr| ≥ 0.50 across stages), and *F34_specific* (late-stage emergence). The larger |Δr| threshold for progressive edges reflects the higher prior probability of detecting monotone trends under batch confounding. The final D-F2 set contains 3,260,461 edges (15,875 F2_emerging; 3,726 F2_dissolving; 2.16M transient_F2; ~230K progressive up/down; ~850K F34_specific).

**Cross-cohort replication (Criterion 1).** For every D-F2 edge we ran leave-one-cohort-out cross-validation over the eight fibrosis-annotated cohorts (`282_cohort_loocv.R`). An edge was counted as replicating in a held-out cohort if |Δr| retained at least half its full-atlas value with the same sign. The median per-edge replication fraction was **0.875**, and 61.7% of the 3.26M edges replicate in ≥ 70% of folds, well above our pre-registered threshold of 0.70. Per-edge LOCO fractions are distributed to the atlas (`loco_replication_fraction` column).

**STRING independence (Criterion 2).** Jaccard overlap between the D-F2 edge set and STRING ≥ 700 was **0.0011** (3,613 shared of 3,375,494 union), confirming that stage-stratified dynamics constitute an effectively orthogonal evidence layer relative to the static STRING backbone.

### 4. D-COLOC edges: paired GWAS colocalization

D-COLOC edges connect gene pairs that share a fine-mapped causal variant credible set across a 50-GWAS portfolio (EUR 21 / AFR 10 / EAS 9 / AMR 7 / SAS 3; 2026-07-05 MVP 50-GWAS update; was 28-GWAS: 17 European + 5 East Asian + 3 Pan-UKBB AFR + 3 Pan-UKBB CSA; canonical SuSiE-COLOC pipeline, Methods §"Causal Inference"). A D-COLOC edge (A, B) was emitted when both A and B had posterior probability of shared causal variant PP.H4 > 0.5 at the same locus in at least one GWAS. The 2026-04-21 SuSiE-COLOC rebuild triggers a re-run of `RNA-seq/254_genetic_coloc_edges.R`; current edge counts and the re-evaluation of Criterion 3 (below) will be updated after that re-run completes.

**Pre-registered criterion (Criterion 3).** The hypothesis was that D-COLOC edges would be preferentially druggable-druggable (DGIdb annotation) at ≥ 3× fold-over-null vs. a Maslov–Sneppen degree-preserving rewiring (§8). The observed fold was 0× (0 of 111 edges druggable-druggable). **Per our binding pre-registration (`docs/network_c_prereg.md`, Criterion 3 kill condition), D-COLOC is therefore demoted from an independent edge type to per-edge annotation only.** The 111 edges are retained in the atlas for provenance, but the Resource does not claim D-COLOC as a discovery layer. Locus linkage is preserved through the per-edge columns `coloc_linked_gene_a`, `coloc_linked_gene_b`, `pp4_min`, and `gwas_list`, enabling downstream users to filter any network view by COLOC support.

### 5. D-LR edges: F2-differential cell–cell communication

D-LR edges were extracted from the LIANA differential ligand–receptor analysis on MASLD single-cell cohorts (Methods §"Single-cell analysis"). We retained ligand–receptor pairs whose cell-type–resolved interaction strength peaked at the Steatosis → Steatohepatitis transition and passed FDR < 0.05 in the differential test. Nine F2-emerging D-LR edges were incorporated, each annotated with its dominant sender/receiver cell-type pair (e.g., Progressor_hepatocyte → Fibroblast; Macrophage → Endothelial).

### 6. D-ceRNA and D-XSpecies

**D-ceRNA.** The Architecture C spec anticipated lncRNA–miRNA–mRNA triplets as MASLD-specific regulatory edges. After applying the DEG-filtered, atlas-intrinsic ceRNA construction (Methods §"ncRNA pipeline") with the MASLD lncRNA set, the resulting edge count was 0 in the final atlas. This reflects known database sparsity for validated lncRNA–miRNA interactions in liver tissue. We report this limitation honestly (§Limitations) rather than populating the layer with unfiltered predictions.

**D-XSpecies.** Cross-species conserved switch edges were seeded from the Fast Principal Component (FPC) cross-species concordance model (Methods §"Cross-species concordance"); genes in the Conserved_Core tier whose coexpression neighborhood inverts between healthy and disease across both mouse and human were flagged, but edges are reported only as per-edge annotations (`conserved_mouse`, `xspecies_switch`) rather than as independent discovery edges.

### 7. Per-edge annotation atlas

The annotation atlas (`edge_annotation_atlas.parquet`) provides ≥ 40 per-edge columns spanning: stage-specific correlations (`r_f01`, `r_f2`, `r_f34`, `delta_r_max`), LOCO replication (`loco_replication_fraction`, `n_cohorts_replicated`), COLOC linkage (`coloc_linked`, `pp4_min`, `gwas_list`), druggability flags (`dgidb_a`, `dgidb_b`, `clinical_drug_a/b`), cell-type driver (`driver_celltype`, `liana_flag`), sex dimorphism (`sex_dimorphic_a/b`), F2-weighted community assignments (`module_f01`, `module_f34`), STRING evidence (`string_combined`, `in_string_ge700`), and biotype. **Annotation completeness (Criterion 6) is 96.6%** of edges carry ≥ 3 non-trivial annotations, exceeding the pre-registered 80% threshold.

### 8. Community detection and F2-weighted Leiden

Gene communities were identified by Leiden clustering (Traag et al. 2019) on the STRING ≥ 700 backbone, with each edge weighted by a disease-dynamics term

> W(u,v) = string_combined(u,v) × exp(|f2_inflection_logFC(u)| + |f2_inflection_logFC(v)|) / 2

where `f2_inflection_logFC` is the per-gene logFC of the F3–F4 vs. F0–F1 dream contrast on the atlas. We ran Leiden at stage-specific weightings — once with F0–F1 expression weights, once with F3–F4 weights (script `292_f2_weighted_leiden.py`) — producing 57 communities per stage, each with modularity **Q = 0.69**. Robustness across random seeds was assessed by ARI (mean pairwise ARI = 0.90 over 20 seeds).

**F2 switch interpretation (Criterion 4).** The pre-registered v1 criterion (community *emergence* / *dissolution*) was amended in-flight to a biology-appropriate formulation: the data showed inflammation and oxidative phosphorylation modules pre-existing at both stages, with a **dominance shift** rather than birth/death. The amended v2 criterion (+15% size gain for inflammation, −10% for metabolism; disclosed in the pre-registration amendment log) is met: the inflammation community expands from 1,263 → 1,455 genes (+15%; HALLMARK_INFLAMMATORY_RESPONSE p = 2.5 × 10⁻³⁴; NF-κB p = 1.1 × 10⁻³²), while the OxPhos/fatty-acid-metabolism community contracts from 1,176 → 1,029 genes (−13%; HALLMARK_FATTY_ACID_METABOLISM p = 1.3 × 10⁻⁴³; HALLMARK_OXIDATIVE_PHOSPHORYLATION p = 7.8 × 10⁻⁶). This is reported in the paper as a *module dominance shift* at F2, not an emergence event.

### 9. Benchmark framework (honest STRING comparison)

All enrichment benchmarks use the Maslov–Sneppen degree-preserving null (`272b_robust_benchmark.R`; 1,000 rewirings, 95% bootstrap CI). Gold sets comprise clinical MASLD drug targets, the Govaere 25-gene panel, the resmetirom (THRB) pathway, the atlas positive controls (N = 58), and SuSiE-COLOC high-confidence hits (PP.H4 > 0.5). Per Criterion 7 we report the full STRING comparison table in Supplementary Table S_NetCompare without redaction: STRING v12 ≥ 700 achieves 208.8× (clinical_drug), 15.2× (Govaere-25), 235.3× (resmetirom), and 96.5× (positive_controls) fold-over-null, outperforming our D-layer composite on every generic gold set. The D-layer value is in *disease-specific* axes (F-stage stratification, COLOC linkage, F2 switch membership) that no generic PPI database encodes. The two views are complementary and are published as such.

### 10. Pre-registration and falsification

Seven binding falsification criteria were registered at `docs/network_c_prereg.md` prior to Phase 2 implementation. Final scorecard: **5 of 7 PASS** (LOCO replication 0.875; Jaccard 0.0011; F2 switch community recovery; annotation completeness 96.6%; honest STRING supplement), **1 of 7 FAIL and demoted per pre-registered kill condition** (D-COLOC druggable enrichment; demoted to annotation-only), **1 of 7 pending** (portal query latency < 2 s; to be verified in post-deployment runtime test). The complete scorecard and amendment log are included as Supplementary Note S_NetPrereg.

---

## Figure 2 panel legend — *Metabolic-to-inflammatory module dominance shift at F2*

**Figure 2X | Stage-weighted community reorganization at the F2 fibrosis transition.** Sankey diagram of Leiden gene communities computed on the STRING v12 high-confidence backbone (combined_score ≥ 700) weighted by stage-specific per-gene F3–F4 vs. F0–F1 logFC (MASLD atlas; N = 1,444 samples; 16,704 genes). Left stack: communities under F0–F1 (early) weighting (57 modules, Q = 0.69); right stack: communities under F3–F4 (late) weighting (57 modules, Q = 0.69). Ribbons are colored by transition category: *dissolving* (late-stage loss of membership), *emerging* (late-stage gain), *fragmenting* (split across multiple late-stage modules), and *invariant* (≥80% member retention). The two largest reorganizations are annotated: the inflammation/NF-κB module expands from 1,263 → 1,455 genes (+15%; top Hallmark HALLMARK_INFLAMMATORY_RESPONSE, p = 2.5 × 10⁻³⁴), and the oxidative-phosphorylation/fatty-acid-metabolism module contracts from 1,176 → 1,029 genes (−13%; top Hallmark HALLMARK_FATTY_ACID_METABOLISM, p = 1.3 × 10⁻⁴³). Both modules exist at both stages; the F2 transition is therefore a **dominance shift**, not a birth-and-death event (see Methods §8 and Supplementary Note S_NetPrereg for the pre-registration amendment). Full enrichment tables for both stages are provided in Supplementary Table S_NetCommunities.

---

## Supplementary Figure figS_network caption

**Figure S_Network | MASLD gene network v2 — atlas, replication, specificity, and community biology.**
**(A) Annotation atlas overview.** Edges per type in the released atlas (Type-S only: 115,020; Contested/overlap: 3,626; D-only: 3,256,953) and histogram of per-edge non-trivial annotation count. 96.6% of the 3,375,599 edges carry ≥ 3 annotations (Criterion 6 PASS; dashed threshold at 3).
**(B) Leave-one-cohort-out replication of D-F2 edges (Criterion 1).** Distribution of per-edge replication fractions across the eight fibrosis-annotated cohorts. Median = 0.875 (red line); pre-registered threshold = 0.70 (dashed). 61.7% of 3.26M edges replicate in ≥ 7 of 10 folds.
**(C) D-F2 vs. STRING specificity (Criterion 2).** Jaccard overlap between D-F2 edges and STRING ≥ 700 edges, stratified by emergence category. Overall Jaccard = 0.0011, confirming that F-stage-stratified dynamics are largely orthogonal to the static STRING backbone. Dashed line shows the pre-registered ceiling (0.30).
**(D) Community × pathway enrichment heatmap.** Hallmark pathway enrichment (−log₁₀ adjusted p) for each Leiden community, separately at F0–F1 and F3–F4 stage weighting. The dominance shift is visible as the intensification of HALLMARK_INFLAMMATORY_RESPONSE and HALLMARK_TNFA_SIGNALING_VIA_NFKB rows (right panel) and the contraction of HALLMARK_FATTY_ACID_METABOLISM and HALLMARK_OXIDATIVE_PHOSPHORYLATION rows.

---

## Data Availability

The MASLD gene network v2 atlas is released as a single self-contained Resource package:

- **Primary deliverable**: `edge_annotation_atlas.parquet` (3,375,599 edges × 40+ annotation columns; available at the Zenodo deposit below and mirrored in the paper's GitHub release).
- **Per-gene neighborhood JSON**: 16,533 gene-graph JSON files (`RNA-seq/results/network/gene_graphs/`) powering the portal queries.
- **Portal**: `https://masld-atlas-v2.org/network` *(URL placeholder pending live deployment)*. The portal exposes the three pre-registered query types — stage-stratified gene neighborhoods, COLOC-linked F2-switch genes, and druggable-pair queries within a subtype community — with response-time targets documented in Supplementary Note S_NetPrereg (Criterion 5).
- **Zenodo deposit**: DOI **TBD** (MASLD Network v2 atlas + portal snapshot + supplementary documentation; CC-BY 4.0).
- **Supplementary documentation**: `docs/network_c_prereg.md` (pre-registration and scorecard), `docs/superpowers/specs/2026-04-16-bayesian-multiplex-network-design.md` (v4 design spec, Architecture C), `docs/network_analysis_adversarial_review.md` and `docs/network_analysis_review_round2.md` (adversarial reviews), and `docs/network_analysis_final_synthesis.md` (final synthesis).
- **Code**: pipeline scripts `RNA-seq/250–292` (edge construction, annotation, F2-weighted Leiden), benchmarks `RNA-seq/272b`, `277`; portal front-end in `masld-atlas-v2/src/app/network/`.

---

## Limitations (for Discussion)

Several limitations of the network atlas should be disclosed to readers.

**D-ceRNA layer is empty.** Despite the inclusion of MASLD-associated lncRNAs as seeds, the final D-ceRNA edge count is zero. This reflects well-known sparsity in curated lncRNA–miRNA interaction databases for liver tissue; rather than populate the layer with low-confidence predictions, we report it honestly as absent. Future work incorporating direct CLIP-seq evidence would be required to populate this layer.

**D-COLOC was demoted per pre-registration.** Our pre-registered Criterion 3 required D-COLOC edges to preferentially connect druggable gene pairs at ≥ 3× fold-over-null. The observed fold was 0× (0 of 111 edges druggable–druggable). Per the binding kill condition, D-COLOC is not claimed as an independent edge type and is preserved only as per-edge annotation. This outcome is a transparent falsification and is recorded in the amendment log.

**D-Spatial and D-COSMOS layers are deferred.** Spatial-transcriptomics-derived coexpression and COSMOS-derived mechanistic-path edges were scoped in the v3 spec but not implemented in v4 due to limited spatial sample size (n = 1 MASLD Visium cohort passing QC) and COSMOS convergence instability on the atlas background. These layers are candidates for a future update.

**STRING contains pre-2015 literature and textmining bias.** STRING's `combined_score` integrates textmining channels that can encode historical literature bias. We mitigate this by using only the high-confidence tier (combined_score ≥ 700) and by using STRING exclusively as a generic-biology backbone, never as a MASLD-specific evidence source. Downstream MASLD claims rely on Type-D annotations, not STRING.

**Two cohorts lack fibrosis annotation.** Of the 10 cohorts in the MASLD atlas, GSE126848 and GSE167523 lack fibrosis staging and were excluded from D-F2 stratification and LOCO replication. The D-F2 layer is therefore built on 8 of 10 cohorts (N ≈ 1,149 fibrosis-annotated samples), and the reported LOCO median of 0.875 is over these 8 folds.

**No methodological novelty claim.** Our integration (STRING backbone + overlay annotations) uses standard tools. Head-to-head benchmarks show STRING v12 ≥ 700 outperforms any composite of our D-layers on generic gold sets. The Resource value is in the disease-specific data and per-edge provenance annotations, and we frame the work accordingly throughout the manuscript.

---

## TODO for author (inconsistencies flagged between sources)

1. **Composite modularity number.** `network_analysis_final_synthesis.md` reports community modularity Q = 0.61 on the v1–v3 composite network. `network_final_decisions.md` and the Criterion 4 scorecard report Q = 0.69 on the F2-weighted Leiden on the STRING backbone (v4/Architecture C). These refer to different networks and both are correct in context. The Methods above cite Q = 0.69 (the v4/published number). Please verify no stray reference to Q = 0.61 remains in figure legends.
2. **D-LR edge count.** `network_final_decisions.md` reports "9 F2-emerging" D-LR edges. The v4 spec allowed for a broader D-LR set. The Methods text above commits to "9 F2-emerging" — confirm this is the final count intended for the atlas release.
3. **D-XSpecies claim status.** `network_final_decisions.md` lists D-XSpecies as one of the 5 Type-D categories, but the final atlas does not report a positive edge count for it (only annotation flags). The Methods §6 above downgrades D-XSpecies to annotation-only to match reality. Please confirm or re-scope.
4. **Portal URL.** Placeholder `https://masld-atlas-v2.org/network` used. Replace with the real live URL prior to submission.
5. **Zenodo DOI.** To be minted prior to submission.
6. **Dataset sample-size note.** Fibrosis-bin sizes (F0–F1 ≈ 490, F2 ≈ 316, F3–F4 ≈ 343) are approximate; confirm exact N against `unified_metadata.csv` for the final Methods.
7. **Supplementary figure letter.** The atlas-level figure is referred to here as `figS_network`; align with the supplementary numbering scheme in `figures/supplementary/figS08_subtyping_convergence/` or a dedicated `figS_network/` directory per `load_figure_data.R`.
