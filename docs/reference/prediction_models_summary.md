# Prediction Models Summary

Single-source lookup for every predictor trained in this project. Post-leakage-fix numbers only.
Last updated: 2026-04-16.

## Conventions
- **Tier 1** = main-figure candidate (novel + defensible + properly nested CV).
- **Tier 2** = supplementary (honest but not headline-grade).
- **Tier 3** = retired / retracted / negative result. Cite only as methodological contrast.
- All metrics are nested / LOCO cross-validated unless noted. "Random" = random gene/protein panel of matched size.

---

## Tier 1 — Main figure candidates

### 1. Biology-only Fibrogenic-subtype classifier
| | |
|---|---|
| Script | `RNA-seq/.../scripts/183_multi_output_predictor.py` + `183v2_ablation.py` |
| Method | Multi-task elastic net predicting 3 clinical outcomes jointly (Fibrogenic subtype, therapeutic pathway, active transition) from 185 features across 6 biology layers (expression, cell-type, TF activity, COLOC, pathway scores, clinical). 6-fold LOCO-CV. |
| **Predicts** | **Fibrogenic NMF cluster membership** (binary: Fibrogenic subtype vs. Metabolic) — Fibrogenic = female-enriched, high-fibrosis, COL15A1/CDH11/LTBP2 ECM-stellate signature |
| **Input features** | **Biology-only config**: bulk RNA-seq DEG expression + BayesPrism cell-type proportions + SCENIC+ TF activity scores + COLOC gene scores + ssGSEA pathway scores (185 features total). Clinical covariates (fibrosis_stage, NAS, age, sex) excluded to prove molecular independence. |
| Headline | **AUROC 0.853** using **biology-only features** (clinical covariates fibrosis_stage/NAS/age/sex removed) |
| Full ablation | All features 0.939 → no fibrosis_stage 0.912 → **biology-only 0.853** → cell-type+clinical 0.879 → clinical-only 0.835 |
| Results CSV | `RNA-seq/results/staging_classifier/multi_output_results.csv` |
| Figure | **Not currently plotted** — needs standalone main-figure panel |
| Tier rationale | Kills circularity worry: removing clinical grade drops AUROC only −0.086 (0.939 → 0.853). Biology features identify Fibrogenic patients independently of fibrosis stage — defensible molecular-discovery claim. |

### 2. Cell-type composition → F≥3 (M2, cleanest baseline)
| | |
|---|---|
| Script | `RNA-seq/.../scripts/161_nested_loco_classifier.py` |
| Method | Elastic net, 17 BayesPrism cell-type proportions, nested LOCO-CV |
| **Predicts** | **Advanced fibrosis (F≥3)** — binary: bridging fibrosis + cirrhosis vs. none/mild fibrosis |
| **Input features** | **17 BayesPrism cell-type deconvolution proportions** estimated from bulk RNA-seq (e.g. hepatocytes, macrophages, HSCs, endothelial, cholangiocytes, etc.). No direct gene expression — composition only. |
| Headline | **AUROC 0.799 (M2)**; M3 expression 0.853 but does NOT beat 500 random genes (p=0.73) |
| Full results | M1 clinical 0.660 / **M2 0.799** / M3 expression 0.853 / M4 combined 0.827 |
| Results CSV | `results/prognosis_v2/nested_loco_summary.csv`, `nested_loco_baselines.csv` |
| Figure | `figures/supplementary/figS10_prediction/panels/figS_staging_classifier.pdf`, `fig_composition_features.pdf` |
| Tier rationale | Zero pre-screening bias, zero label leakage, proper nested LOCO — defensible as "cleanest single predictor" |

### 3. Cross-ancestry COLOC genetic classifier
| | |
|---|---|
| Script | `RNA-seq/.../scripts/172v2_cross_ancestry_classifier.py` |
| Method | Elastic net on 66 COLOC genes (PP.H4>0.5 in both EUR + EAS GWAS), 6-fold LOCO, permutation + random-gene null |
| **Predicts** | **Advanced fibrosis (F≥3)** — binary |
| **Input features** | **Bulk RNA-seq expression of 66 genetically-validated genes** — specifically genes with COLOC PP.H4>0.5 in BOTH European (FinnGen/Broadaway) AND East Asian (BBJ) GWAS. Gene set is genetically prioritized, not data-driven. |
| Headline | **AUROC 0.862 ± 0.053** (permutation p=0.010, beats random 66 genes p=0.020) |
| Results CSV | `RNA-seq/results/novel_ml/cross_ancestry/cross_ancestry_classifier_results_v2.csv` |
| Figure | `figures/main/fig4_causal_architecture/panels/fig3f_multiancestry.pdf`; `figures/supplementary/figS09_multi_ancestry/figS09_multi_ancestry_coloc.pdf` |
| Tier rationale | Genetic prioritization generalizes cross-ancestry; one of the few models that beats random gene panels |

### 4. Plasma MASLD binary (T5)
| | |
|---|---|
| Script | `RNA-seq/.../scripts/225e_plasma_sweep.py` + `226a_interpretation.py` |
| Method | BalancedBagging, 1,461 Olink proteins, 10×5-fold nested CV, N=177 (MASLD=81, CVH=82, ARLD=14) |
| **Predicts** | **MASLD diagnosis** (binary: MASLD vs. non-MASLD [CVH + ARLD combined]) — non-invasive blood-based diagnosis |
| **Input features** | **1,461 plasma proteins** measured by Olink proximity extension assay. Blood-only — no tissue biopsy, no imaging, no clinical staging needed. |
| Headline | **AUROC 0.840 ± 0.058** |
| Key finding | Directly contradicts Yang et al. 2025 Cell Rep Med "fibrosis biology converges across etiologies"; etiology signal is diffuse (needs >200 proteins); CHI3L1/YKL-40 dominates fibrosis (SHAP=0.398) |
| Results CSV | `RNA-seq/.../results/multiprogram/plasma_sweep_225e_*.csv`, `226a_protein_importance.csv`, `226b_tissue_bridge.csv` |
| Figure | `figures/supplementary/figS10_prediction/figS16_plasma_proteomics.pdf` + panels A–D; Yang panel = `figS16_panelD_yang_comparison.pdf` |
| Tier rationale | Novel contradiction of published claim + tight SD + nested CV |

---

## Tier 2 — Supplementary worthy

### 5. Concept bottleneck (fibrosis + disease)
| | |
|---|---|
| Script | `RNA-seq/.../scripts/106_concept_bottleneck.py` (post-leakage-fix 2026-04-08) |
| Method | Concept bottleneck with macrophage/cell-type concepts; per-fold concept definitions |
| **Predicts** | **Fibrosis stage F0–F4** (ordinal QWK) AND **disease binary** (healthy vs. MASLD, AUROC) — jointly via concept bottleneck |
| **Input features** | **Bulk RNA-seq expression** → first mapped through interpretable "concept" nodes (macrophage activation, hepatocyte function, steatosis, etc. derived from cell-type deconvolution) → then concept activations predict outcome. Human-interpretable intermediate layer. |
| Headline | **Fibrosis QWK 0.633; Disease AUROC 0.828** |
| Attribution | Macrophage proportions dominate NAS (15.9) vs fibrosis (4.6) |
| Results CSV | `RNA-seq/.../results/multiprogram/106_concept_bottleneck_*.csv` |
| Figure | `figures/supplementary/figS10_prediction/panels/figS10_panel_d_concept_attribution.pdf` |

### 6. Multi-program NAS ordinal (expression elastic net)
| | |
|---|---|
| Script | `RNA-seq/.../scripts/190_multiprogram_ablation.py` |
| Method | Elastic net, per-fold gene selection + standardization, NAS 0-8 (N=660, 5 folds) |
| **Predicts** | **NAS score 0–8** (ordinal; components: steatosis 0-3 + lobular inflammation 0-3 + ballooning 0-2) — continuous disease activity, not fibrosis stage |
| **Input features** | Ablation across 4 modalities: (1) **bulk RNA-seq expression** (per-fold DEG selection), (2) **COLOC gene expression** (genetically prioritized), (3) **ssGSEA pathway scores** (hallmark gene sets), (4) **combined** (all modalities). Expression is the only modality beating random. |
| Headline | **NAS QWK 0.556 (p=0.010 vs random)** — only modality beating random |
| Full ablation | Expression 0.556 / combined 0.547 / ssGSEA 0.434 / COLOC 0.384 |
| Results CSV | `RNA-seq/.../results/multiprogram/multiprogram_ablation.csv` |
| Figure | `figures/supplementary/figS10_prediction/panels/figS10_panel_a_modality_heatmap.pdf`, `figS10_panel_b_random_ceiling.pdf` |

### 7. Multi-task CORN (fibrosis + NAS joint ordinal)
| | |
|---|---|
| Script | `RNA-seq/.../scripts/92_multitask_corn.py` (post-leakage-fix) |
| Method | CORN ordinal multi-task loss, per-fold gene selection |
| **Predicts** | **Fibrosis F0–F4** (ordinal QWK) AND **NAS≥5** (binary AUROC) — simultaneously in one model with shared representation |
| **Input features** | **Bulk RNA-seq expression** of per-fold selected DEGs. Same feature space as Model 6 but jointly optimized for both fibrosis and NAS outcomes using CORN ordinal ranking loss. |
| Headline | Fibrosis QWK 0.369; NAS≥5 AUROC 0.718 |
| Results CSV | `RNA-seq/.../results/multiprogram/92_multitask_corn_results.csv` |
| Figure | `figures/supplementary/figS10_prediction/panels/figS_staging_classifier.pdf` |

### 8. Plasma F≥3 (XGBoost T1)
| | |
|---|---|
| Script | `225e_plasma_sweep.py` (T1 target) |
| Method | XGBoost, 1,461 Olink proteins, 10×5-fold nested CV |
| **Predicts** | **Advanced fibrosis (F≥3)** — binary, using only blood-based measurements |
| **Input features** | **1,461 plasma Olink proteins** (same panel as Model 4). Unlike Model 4 which predicts MASLD diagnosis, this predicts fibrosis severity within MASLD patients. |
| Headline | AUROC 0.790 ± 0.080 |
| Key finding | CHI3L1 SHAP 0.398 (90% fold stability, 9.7% of total signal); focal biology |
| Results CSV | `plasma_sweep_225e_T1.csv`, `226a_shap_T1_binary.npz` |
| Figure | `figS16_panelA_chi3l1_dominance.pdf`, `figS16_panelB_panel_curves.pdf` |

### 9. Plasma etiology (T4)
| | |
|---|---|
| Script | `225e_plasma_sweep.py` (T4 target) |
| Method | Ridge / linear SVC, 3-class etiology (MASLD / CVH / ARLD) |
| **Predicts** | **Liver disease etiology** — 3-class: MASLD (metabolic) vs. CVH (cardiovascular/healthy) vs. ARLD (alcohol-related). Blood-based etiology discrimination. |
| **Input features** | **1,461 plasma Olink proteins** (same panel as Models 4 and 8). Key insight: etiology signal is diffuse — random 20 proteins beats SHAP-selected 20 (0.724 vs 0.686), suggesting no dominant protein signature. |
| Headline | macro F1 0.657–0.668 |
| Key finding | Diffuse signal — random 20 proteins beats SHAP-selected 20 (0.724 vs 0.686); metabolic signature (NPY, REN, GCG, PRKCQ) |
| Results CSV | `plasma_sweep_225e_T4.csv`, `226d_etiology_pathways.csv` |
| Figure | `figS16_panelC_etiology_proteins.pdf` |

### 10. Sex-specific F≥3 classifier
| | |
|---|---|
| Script | `RNA-seq/.../scripts/26_sex_stratified_analysis.R` + staging classifier sex variant |
| Method | Sex-stratified elastic net, separate panels per sex |
| **Predicts** | **Advanced fibrosis (F≥3)** — binary, separately within female cohort and within male cohort |
| **Input features** | **Bulk RNA-seq expression** — but gene panels selected independently per sex via sex-stratified LOCO-CV. Zero gene overlap between panels demonstrates sex-specific transcriptomic fibrosis programs. |
| Headline | Female AUROC 0.898 / Male AUROC 0.808; **0 gene overlap** between panels |
| Results CSV | `sex_auroc_comparison.csv`, `sex_specific_panels.csv` |
| Figure | `figures/main/fig5_sex_dimorphism/panels/` + `figures/supplementary/figS10_prediction/panels/figS_staging_classifier.pdf` |

### 11. Plasma NFASC+GDF15 benchmark
| | |
|---|---|
| Script | `226a_interpretation.py` (benchmark subroutine) |
| Method | 2-protein logistic, LOEO pooled |
| **Predicts** | **Advanced fibrosis (F≥3)** — binary, replicating Yang et al. 2025 reported panel |
| **Input features** | **2 plasma proteins only: NFASC + GDF15** — the exact panel reported by Yang et al. as achieving 87% balanced accuracy. Our nested CV LOEO gives AUROC 0.734 ± 0.070, substantially lower than the published claim. |
| Headline | AUROC 0.734 ± 0.070 |
| Context | Yang et al. reported 87% balanced accuracy; our nested CV is more honest |
| Results CSV | `plasma_nfasc_gdf15_benchmark.csv` |
| Figure | Shared panels in `figS16_plasma_proteomics.pdf` |

### 12. NMF Fibrogenic subtype (expression M3)
| | |
|---|---|
| Script | `161_nested_loco_classifier.py` (loco_s2 target) |
| Method | Elastic net, top-500 div_* genes per fold, nested LOCO |
| **Predicts** | **Fibrogenic NMF cluster membership** (binary: Fibrogenic vs. Metabolic) — same target as Model 1 but using a pre-screened expression-only feature set |
| **Input features** | **Top-500 "divergence" genes** (div_* columns) selected per fold from bulk RNA-seq — genes showing maximal Metabolic vs Fibrogenic subtype divergence. Note: div_* pre-screening performed on full data before CV splits, introducing optimistic bias. Use Model 1 (0.853) as primary cite. |
| Headline | AUROC 0.904 (random gene p=0.04, barely significant) |
| Caveat | div_* pre-screened on full data; use Fibrogenic molecular-only from Script 183v2 (0.853) as primary cite |
| Results CSV | `nested_loco_summary.csv` (loco_s2 rows) |
| Figure | `figures/supplementary/figS10_prediction/panels/panel_e_bifurcation.pdf` |

---

## Tier 3 — Retired / retracted / negative results

| Model | Script | Predicts | Input features | Metric | Status |
|-------|--------|----------|----------------|--------|--------|
| F≥3 binary (any signature) | 61 / 70 | F≥3 binary | Bulk RNA-seq expression (curated DEG panel) | AUROC 0.853 curated vs 0.860 random (p=0.733) | **RETIRED** — fibrosis is transcriptome-pervasive (publishable as negative result) |
| Fibrosis QWK expression enet | 190 | Fibrosis F0-F4 ordinal | Bulk RNA-seq expression | QWK 0.460 vs random 0.469 (p=0.614) | Retired; pair with F≥3 as negative result |
| COLOC panels for ordinal prediction | 192 | Fibrosis F0-F4 / NAS 0-8 ordinal | Expression of COLOC-prioritized genes (PP4>0.5 or PP4>0.8 panels) | QWK 0.384–0.477 (all p>0.5 vs random) | Negative result — COLOC ≠ predictive gene set |
| Plasma T2 ordinal (3-class fibrosis) | 225 | Fibrosis 3-class ordinal (F0-1 / F2 / F3+) | 1,461 Olink plasma proteins | QWK 0.469 ± 0.137 | **DROPPED** — F3 N=24 too small |
| Plasma cross-platform Olink→DIA-MS | 226c | F≥3 binary (cross-platform transfer) | 48 proteins shared between Olink (N=177) and DIA-MS (N=72) | AUROC 0.655 | **DROPPED** — 48 shared proteins only (3.3%); platform incompatibility |
| Therapeutic pathway multi-task | 183 | Therapeutic pathway assignment (binary) | Same 185 features as Model 1 | F1 0.777 | Framing only: clinical staging drives it (molecular-only 0.451) |
| Active transition multi-task | 183 | Active disease transition (binary) | Same 185 features as Model 1 | F1 0.439 | Weak signal regardless of features |
| Stacking ensemble | — | F≥3 binary | VAE latent embeddings + expression | QWK 0.851 | **RETRACTED** — VAE leakage; clean version 0.343 |
| BulkFormer transfer (Script 167) | 167 | F≥3 binary | BulkFormer 643-dim embeddings from bulk RNA-seq | AUROC ~0.969 | **RETRACTED** — was silent VAE fallback |
| Phase VII v1 prognosis | 151, 153, 154, 155 | F≥3 / HCC risk / progression | Expression + cell-type (no proper nesting) | Cell-type 0.865, 5-gene 0.863, C-index 0.638 | **DEPRECATED** — no nested CV, inflated metrics; superseded by Script 161 |
| Cirrhosis risk head | 152 | Cirrhosis risk (binary) | Pseudo-labels derived from same data | AUROC 0.431 | **DROPPED** — below chance, circular pseudo-labels |

Machine-readable status: `results/prognosis_v2/model_status.csv`.

---

## Figure directory index

| Figure | Location | Source script |
|--------|----------|---------------|
| SFig 10 assembled (multi-program) | `figures/supplementary/figS10_prediction/figS10_multiprogram.pdf` | `scripts/figures/figS10_multiprogram.R` |
| SFig 10 panels (14 PDFs) | `figures/supplementary/figS10_prediction/panels/` | — |
| SFig 16 plasma (4 panels) | `figures/supplementary/figS10_prediction/figS16_plasma_proteomics.pdf` | `scripts/figures/figS16_plasma_proteomics.R` |
| Fig 4 multi-ancestry | `figures/main/fig4_causal_architecture/panels/fig3f_multiancestry.pdf` | — |
| Fig 5 sex dimorphism | `figures/main/fig5_sex_dimorphism/panels/fig5_panel_{A-E}.pdf` | — |
| SFig 9 multi-ancestry COLOC | `figures/supplementary/figS09_multi_ancestry/figS09_multi_ancestry_coloc.pdf` | — |
| **MISSING: biology-only progressor-classifier main panel** | — | Script 183v2 outputs exist; panel TBD |

---

## Suggested main-figure integration

- **Fig 3 (Deconfounding)**: adopt M2 cell-type AUROC 0.799 (Model 2) from SFig 10 as the unbiased-baseline headline.
- **Fig 4 (Causal Architecture)**: keep cross-ancestry COLOC 0.862 (Model 3); already present.
- **Fig 5 (Sex)**: keep sex-specific 0.898/0.808 (Model 10); already present.
- **Fig 6 (Therapeutic Windows)**: add new panel for biology-only progressor classifier 0.853 (Model 1) — currently has no dedicated figure.
- **Fig 2 (F2 switch) or Fig 6**: CHI3L1 focal + etiology diffuse (Models 4, 8, 9) for "F2 switch visible in plasma" story.
- **SFig 10**: keep Models 5–7, 11–12 + Tier 3 negative results (random-gene ceiling) as the honest-reporting section.
