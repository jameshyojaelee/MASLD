# Gandal 2022 (ASD cortex) vs. MASLD atlas — DE methods, step by step

**Reference paper:** Gandal, Haney et al., *Broad transcriptomic dysregulation occurs across the cerebral cortex in ASD*, Nature 2022 (doi 10.1038/s41586-022-05377-7).
**Their code:** https://github.com/dhglab/Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD
**Purpose:** map every differential-expression method we have *actually run* against Gandal's actual code, with the math at each pipeline step, and explain how they differ.

> One-line summary: our canonical **`dream`** is the direct successor to Gandal's core engine (`duplicateCorrelation` + `lmFit`, both from the `variancePartition`/limma lineage). We additionally ran a full method bake-off (`mega_validation/`) that *reconstructs* Gandal's exact engine in two of its arms (`07_limma_trend` = their moderation; `02_voomLmFit_block` = their random-effect structure). The biggest genuine divergences are **quantification** (STAR→RSEM vs STAR→featureCounts), **GC/length-bias correction** (cqn vs none), and **technical-covariate selection** (EARTH/MARS vs random-effect absorption).

---

## 0. Method roster — every DE method we've run, and its Gandal analog

| Our method | Script (`.../integration/scripts/`) | Engine / test | Gandal analog |
|---|---|---|---|
| **dream (CANONICAL)** | `05_dream_mega_analysis.R` | voom weights + per-gene LMM `(1\|dataset)`, moderated-t w/ Satterthwaite df | Successor of their `duplicateCorrelation`+`lmFit` |
| ashr shrinkage | `05b_ashr_shrinkage.R` | adaptive shrinkage of dream LFC → lfsr | (none — they keep raw moderated-t) |
| per-study limma-voom | `02_per_study_de.R` | voom → `lmFit` → `eBayes` per cohort | Same family; their analysis is single-dataset |
| metafor RE (reference) | `06_meta_analysis.R` | 2-stage REML random-effects meta, I²/τ²/Q | (none for DE; they use sLDSC for genetics) |
| edgeR-QL (FE cohort) | `mega_validation/01_edgeRql_cohortFE.R` | NB + quasi-likelihood F-test | — |
| **voom + duplicateCorrelation** | `mega_validation/02_voomLmFit_block.R` | 2-round voom + `block=dataset` consensus ρ | **≈ their exact RE structure** (block=Subject) |
| metafor REML + HKSJ | `mega_validation/03_metafor_REML_HKSJ.R` | 2-stage REML + small-K t correction | — |
| mashr | `mega_validation/04_mashr_sharing.R` | multivariate adaptive shrinkage, K cohorts | — |
| DESeq2 Wald (FE cohort) | `mega_validation/06_DESeq2_cohortFE.R` | NB GLM, Wald test | — |
| **limma-trend** | `mega_validation/07_limma_trend.R` | log-CPM → `eBayes(trend=TRUE, robust=TRUE)` | **= their exact moderation** (trend+robust eBayes) |
| DESeq2 LRT | `mega_validation/08_DESeq2_LRT.R` | NB GLM, likelihood-ratio test | — |
| edgeR LRT | `mega_validation/09_edgeR_LRT.R` | NB GLM, likelihood-ratio test | — |
| voom + SVA | `mega_validation/10_voom_sva.R` | surrogate-variable batch removal + voom | Conceptually parallels their MARS covariate hunt |
| mashr (sex) | `sex_v3/05_mashr_fit.R` | multivariate shrinkage across sex×stage | — |
| mouse metafor | `Mouse/.../M03_mouse_meta_analysis.R` | 2-stage REML across diet models | — |
| Cas13 library spine | `Cas13_Library_Design/.../06b_meta_ashr_shrinkage.R` | limma-voom + metafor REML + ashr | — |

**Gandal's own engine:** `02_DEGenesIsoforms/02_01_A_DEGenes.R` — `duplicateCorrelation(block=Subject)` → `lmFit(block=Subject, correlation=consensus)` → `contrasts.fit` → `eBayes(trend=TRUE, robust=TRUE)`. snRNA arm: `06_snRNAseq/NEBULA_pipeline_011122.R` — NEBULA negative-binomial mixed model.

---

## 1. Quantification (reads → gene/isoform counts)

| | Gandal 2022 | MASLD (canonical) |
|---|---|---|
| Aligner | STAR (GRCh37 / GENCODE v25), `--quantMode TranscriptomeSAM` | STAR (GRCh38.p14 / GENCODE v49) |
| Counting | **RSEM** `rsem-calculate-expression` (EM over transcriptome) → gene **and isoform** | **featureCounts** `-s 2` (reverse-stranded), gene-level |
| Sensitivity arm | — | kallisto (pseudoalignment) as a parallel arm |
| QC metrics | Picard `CollectRnaSeqMetrics`/`CollectGcBiasMetrics`/`CollectInsertSizeMetrics` + STAR logs | MultiQC over STAR/featureCounts |

**Math difference.** RSEM resolves multi-mapping reads with an EM that maximizes the transcriptome likelihood `L(θ) = Π_reads Σ_transcripts P(read | transcript)·θ_transcript`, giving fractional (expected) transcript counts and a native **isoform** layer. featureCounts assigns each read to exactly one gene by union-exon overlap (no isoform deconvolution). Consequence: Gandal can run a whole second **isoform/DTE** analysis; our bulk layer is gene-level only (we recover transcript resolution elsewhere via kallisto + scRNA).

---

## 2. Gene filtering (which genes enter the model)

| | Gandal 2022 | MASLD canonical (`03_integrate_counts.R:67`) | Bake-off arms |
|---|---|---|---|
| Rule | RSEM expression cutoff + biotype filtering on the cqn'd matrix | `filterByExpr(dge, design = ~0 + group_binary)` | same `filterByExpr`, sized on `group_binary` |
| Universe | per their cortex dataset | 9-cohort `filterByExpr` universe, then mega on 5 cohorts | identical |

**Math.** `filterByExpr` (Chen 2016) keeps gene g if it has ≥`min.count` (default 10) CPM-equivalent reads in ≥ the smallest group size n_min samples — i.e. `rowSums(cpm > k) ≥ n_min` with `k = min.count / median(libsize) · 1e6`. Gandal's filter is a manual expression floor rather than `filterByExpr`. Minor and not outcome-driving.

---

## 3. Normalization (library size, composition, and GC/length bias)

| | Gandal 2022 | MASLD canonical | Bake-off |
|---|---|---|---|
| Composition/library | `calcNormFactors` (TMM) | **`calcNormFactors(method="RLE")`** (`03:79`) | TMM (arms 01/02/07) |
| **GC + length bias** | **cqn** (Hansen 2012) → per-gene, per-sample offset | **none** | none |
| Per-obs precision | — (handled at moderation, see §4) | voom weights | voom (01/02), none (07), NB (06) |
| Final matrix | `cpm(log=TRUE, offset = cqn.dat$offset, prior.count=0.5)` | log2-CPM via voom | `cpm(log=TRUE, prior.count=3)` (arm 07) |

**Math difference — this is a real one.** cqn fits, per sample, a smooth (natural-spline) regression of log2 count on **GC content** and **gene length**, `log2 count_gi = s_i^{GC}(GC_g) + s_i^{len}(len_g) + offset`, and exports the systematic part as a per-gene/per-sample offset folded into the GLM. This removes within-sample GC/length trends that TMM/RLE (which are single global scale factors per sample) cannot touch. **We do not correct GC/length bias at all** — STAR→featureCounts→voom assumes it is negligible or absorbed by `(1|dataset)`. RLE vs TMM is a second-order difference (both are global per-sample scale factors; RLE = DESeq2 median-of-ratios, chosen here for robustness under asymmetric group sizes).

---

## 4. Mean–variance modeling (count heteroscedasticity)

| Method | How variance enters | Where in our code |
|---|---|---|
| **voom (dream, arms 01/02, per-study)** | per-observation precision weight `w_gi = ŝ(μ̂_gi)^{-4}`, where `ŝ(·)` is a LOWESS of √(residual SD) on mean log-CPM | `voomWithDreamWeights` (`05:132`), `voom` (`02:76`) |
| **limma-trend (arm 07) = GANDAL** | no per-obs weights; instead the eBayes prior variance `s0²` is made a function of mean log-expression (`trend=TRUE`) | `mega_validation/07_limma_trend.R:60` |
| **NB GLM (DESeq2 / edgeR)** | variance `= μ + α·μ²` modeled explicitly via gene-wise dispersion α, EB-shrunk toward a fitted trend | arms 01/06/08/09 |

**Difference.** Gandal handles heteroscedasticity through **limma-trend** (mean-variance trend baked into the empirical-Bayes prior) on a cqn-offset log-CPM matrix. Our canonical uses **voom** (mean-variance pushed into per-observation weights). Our arm `07_limma_trend` is the *direct* test of "does voom's weighting beat Gandal-style trend-only moderation on our data?" — same eBayes call (`trend=TRUE, robust=TRUE`), differing only in the absence of voom weights and the use of FE-cohort instead of cqn offsets.

---

## 5. Technical-covariate selection & batch handling

| | Gandal 2022 | MASLD canonical | Bake-off |
|---|---|---|---|
| Strategy | **EARTH/MARS** forward selection over STAR+Picard QC metrics, cross-validated (10× / 90% subsamples), selected by median per-gene R² | **random-effect absorption**: `(1\|dataset)` soaks up inter-study technical variance; lean fixed set `group + sex` | FE cohort `~ dataset + sex + group` (arms 01/06/07); SVA (arm 10) |
| Collinearity prune | drop one of any pair with adj-R² ≥ 0.95; keep the one most associated with top-5 expression PCs (`01_02_A_CountsProcessing.R:323–395`) | n/a (few fixed covariates) | n/a |
| Covariate vetting | `variancePartition` variance decomposition before fixing the model | `04_variance_partition.R`, `dream_residual_pvca.R` | — |
| Covariates carried | subject(RE) + diagnosis×region + SeqBatch + Sex + Ancestry + Age + Age² + PMI + RIN + ~13 Picard/STAR metrics | group + sex (+ `(1\|dataset)`); deconvolution covariates *deliberately excluded* (`05:88`) | + dataset, + SVs (arm 10) |

**Difference — the second real one.** Gandal *explicitly models* ~15 technical covariates chosen by a cross-validated MARS search. We *absorb* technical variance into a single random intercept per dataset and otherwise keep only biological fixed effects. Both are defensible; ours is the standard mega-analysis move (per-study technical covariates are not harmonized across 5 independent GEO cohorts, so a study-level RE is the principled summary). Arm 10 (`voom_sva`) is our closest analog to their covariate hunt: SVA *learns* latent technical factors instead of selecting observed QC metrics.

A subtle scope difference: their RE (`block=Subject`) models **repeated brain regions within one donor** in a single dataset. Our RE (`(1|dataset)`) models **repeated samples within a cohort** across a mega-analysis. Same compound-symmetry math (§6), different biological grain.

---

## 6. The repeated-measures / multi-study model (the core)

This is the step where `dream` and Gandal are siblings.

| Approach | Model for gene g | Random-effect treatment | Code |
|---|---|---|---|
| **Gandal: duplicateCorrelation** | `y_g = Xβ_g + ε_g`, GLS with compound-symmetric blocks | **single genome-wide consensus** intra-block correlation ρ (one ρ for ALL genes), `block=Subject` | `02_01_A_DEGenes.R:26–28` |
| **MASLD: dream (CANONICAL)** | `y_g = Xβ_g + Z b_g + ε_g`, `b_g ~ N(0, σ²_{RE,g})` | **per-gene** random-intercept variance `(1\|dataset)` fit by `lmer` (weighted by voom) | `05_dream_mega_analysis.R:112,135` |
| **MASLD arm 02 = Gandal idiom** | same GLS as Gandal | single consensus ρ, **`block=dataset`** (2-round voom) | `02_voomLmFit_block.R:74–85` |
| **MASLD arms 01/06/07: FE cohort** | `y_g = Xβ_g + ε_g` with **dataset as fixed factor** | none (cohort = fixed dummies) | arms 01/06/07 |
| **MASLD metafor (06 / arm 03): 2-stage** | stage-2 `β_{gk} = μ_g + u_{gk} + e_{gk}`, `u ~ N(0, τ²_g)` | **per-gene τ²_g** across cohorts (REML) | `06_meta_analysis.R:87` |

**Math difference (duplicateCorrelation vs dream).** `duplicateCorrelation` estimates one consensus correlation `ρ̄ = tanh(mean(atanh(ρ_g)))` and plugs the compound-symmetric correlation matrix `V = (1−ρ̄)I + ρ̄ J` into a single GLS for every gene. `dream` instead fits a separate mixed model per gene, so each gene gets its *own* random-effect variance — strictly more flexible, and the documented win case is genes whose cohort-RE variance departs from the genome-wide average (Hoffman & Roussos 2021). The FE-cohort arms drop the RE entirely; with K=5 balanced-ish cohorts the FE-vs-RE distinction is statistically minor (Hodges 2014), which our diagnostic confirmed (`05_diagnostic_5cohort_re_vs_fe.R`, ρ≈0.999 between 8- and 5-cohort dream fits).

**Mega-analysis vs meta-analysis.** dream and the FE/RE arms are **one-stage** (all samples in one model). metafor (`06`, arm `03`) and mashr (arm `04`) are **two-stage** (per-cohort summary stats → combine). Gandal is one-stage. Our reference metafor exists precisely to show the one-stage dream is robust to the one-stage/two-stage choice (concordance ρ = 0.887, direction 100%, recorded in `dream_vs_meta_concordance.csv`).

---

## 7. Estimation & test statistic

| Method | Estimator | Test statistic | df / null |
|---|---|---|---|
| Gandal | GLS (consensus ρ) + limma EB | moderated t on contrast, `eBayes(trend,robust)` | EB-augmented residual df |
| dream | weighted `lmer` per gene | moderated t on fixed effect | **Satterthwaite** approx df (`05:136–138` — no second `eBayes`) |
| arm 02 | GLS (consensus ρ) + EB | moderated t | EB residual df |
| arm 07 | OLS on log-CPM + EB | moderated t | EB residual df (trend prior) |
| edgeR-QL (01) | NB QL | **quasi-likelihood F** | EB-squeezed QL dispersion, F-df |
| DESeq2 Wald (06) | NB GLM | **Wald** z = β̂/SE | NB asymptotic |
| DESeq2 / edgeR LRT (08/09) | NB GLM | **likelihood-ratio** χ² | df = #dropped params |
| metafor REML (06) | inverse-variance RE | z = μ̂/SE(μ̂) | normal |
| metafor REML+HKSJ (03) | inverse-variance RE | **t on K−1 df**, SE inflated by q | small-K t (`03:70 test="knha"`) |

**Note on dream.** Our script deliberately does **not** call `eBayes()` after `dream()` — dream already returns moderated statistics via a Satterthwaite-df approximation; a second `eBayes` would double-shrink (`05_dream_mega_analysis.R:136–138`). Gandal *does* call `eBayes(trend=TRUE, robust=TRUE)` because their `lmFit` path produces *unmoderated* fits that need the limma EB step. So both end at a moderated-t, by different routes.

---

## 8. Effect-size shrinkage

| | Gandal 2022 | MASLD |
|---|---|---|
| LFC shrinkage | **none** — raw moderated logFC | **ashr** (`05b`): unimodal half-uniform prior, posterior-mean shrunk LFC + **lfsr** |
| Multivariate | — | **mashr** (arm 04 / `sex_v3`): K-condition covariance-learned shrinkage + sharing classes |
| Significance axis | adj.P.Val (BH) on moderated t | **lfsr < 0.05 & \|shrunk_logFC\| > 0.5** (canonical 2026-06-02) |

**Math.** ashr (Stephens 2017): prior `g(β)=Σ_k π_k Unif/N(0,σ_k²)`, likelihood `β̂_g | β_g ~ N(β_g, s_g²)`, EB-estimate the mixture weights π, report posterior mean and **lfsr** = P(wrong sign | data) — a calibrated alternative to a p-value that also down-weights small effects. Gandal never shrinks LFC; they threshold on BH-adjusted p alone. This is a genuine philosophical difference: we make the **effect size** the inferential object (lfsr + shrunk LFC), they make the **p-value** the object.

---

## 9. Multiple-testing & thresholds

| | Gandal 2022 | MASLD canonical |
|---|---|---|
| Correction | Benjamini–Hochberg FDR (`topTable`) | BH on dream p (atlas); **lfsr** for primary DEG calls |
| Primary cut | FDR < 0.05 | lfsr < 0.05, \|shrunk_logFC\| > 0.5 |
| Effect floor | none | 0.5 (LOO-CV stability plateau) |

Both use BH. We add an effect-size floor justified by leave-one-out CV stability; they rely on FDR alone (with separate region-vs-whole-cortex effect-size comparison by total-least-squares regression).

---

## 10. Heterogeneity quantification

Only the meta-analytic arms produce this; Gandal's one-stage model has no explicit τ²/I².

| | Definition | Where |
|---|---|---|
| τ²_g | REML between-cohort variance of LFC | `06:92`, `03:tau2` |
| I²_g | `τ² / (τ² + s²)` — fraction of total variance from heterogeneity | `06_meta_analysis.R` (median I² ≈ 73%) |
| Cochran Q | `Σ w_k (β_k − μ̂)²`, χ² on K−1 | `06`, `03` |
| HKSJ | t-test on K−1 df with q-inflated SE (small-K coverage) | `03:70` |

Gandal's analog of "is this effect consistent across units?" is their **bootstrap of region-specific slopes** (`02_04_A_RegSpecificASD_Bootstrap.R`) and total-least-squares regression of region vs whole-cortex LFC — a within-dataset cross-region consistency check rather than a cross-study τ².

---

## 11. Downstream (shared in spirit, different tools)

| Step | Gandal 2022 | MASLD |
|---|---|---|
| Co-expression modules | rWGCNA (resampled consensus TOM) + isoform `blockwiseModules` | WGCNA in Cross-Species Concordance pipeline |
| Module ↔ trait | limma on module eigengenes (`duplicateCorrelation` again) | NMF k-program eigengene associations |
| Cell-type enrichment | EWCE (bootstrap specificity) | MuSiC/BayesPrism deconvolution + scRNA pseudobulk |
| Variance decomposition | `variancePartition` | `variancePartition` (same package) + PVCA |
| Genetic integration | **sLDSC** partitioned heritability | COLOC / TWAS / INTACT (gene-level) |
| snRNA DE | **NEBULA** NB-mixed (`~Diagnosis+Age+Sex+nGene+percent_mito`, id=subject) | pseudobulk + LIANA |
| Effect-size comparison | total-least-squares (orthogonal) regression + bootstrap CI | cross-species/cohort LFC concordance scatter |

---

## 12. Bottom line — how to talk about it

1. **`dream` is the right primary answer.** It is the modern generalization of Gandal's `duplicateCorrelation`+`lmFit` (same `variancePartition`/limma lineage, same compound-symmetry math, upgraded from one genome-wide ρ to a per-gene random-effect variance).
2. **We already ran Gandal's exact engine as validation.** Arm `07_limma_trend` reproduces their *moderation* (`eBayes(trend=TRUE, robust=TRUE)`); arm `02_voomLmFit_block` reproduces their *random-effect structure* (`block=`, consensus correlation). Neither beats dream on our data; both serve as reviewer armor.
3. **Three genuine divergences, not similarities:**
   - **Quantification:** they use STAR→**RSEM** (with a full isoform/DTE layer); we use STAR→**featureCounts** gene-level (+ kallisto sensitivity).
   - **GC/length bias:** they correct it with **cqn**; we do not correct it at all.
   - **Technical covariates:** they **select ~15 via cross-validated MARS**; we **absorb** them into `(1|dataset)` and keep biology-only fixed effects (with SVA as arm 10 analog).
4. **Inferential philosophy differs:** we threshold on **lfsr + shrunk effect size** (ashr/mashr); they threshold on **BH-adjusted p** of unshrunk moderated-t.

### Two things worth borrowing as sensitivity checks (not adopting wholesale)
- **cqn** GC/length-bias correction — the one normalization step they have and we lack; cheap to add as a sensitivity arm on the merged matrix.
- **MARS / variancePartition covariate vetting** — confirm that `(1|dataset)` alone is absorbing the technical variance a MARS search would have flagged, rather than assuming it.
