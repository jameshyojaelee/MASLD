> **Historical competitor/source dossier only (2026-08-10).** Source-paper summaries may remain useful, but all “our work,” winner, figure, journal, target, and next-analysis sections predate the standalone Resource scope and may be stale. Reverify source facts, and use [`../PAPER.md`](../PAPER.md), [`../STATUS.md`](../STATUS.md), and [`../RESULTS.md`](../RESULTS.md) for current project claims.

# MASLD Competitor Paper Detailed Extractions

> 14 papers. Structured for LLM ingestion. Fields: Citation, DOI/Accession, Samples, Methods, Key Findings, Datasets Used, Limitations, Overlap with Our Work.

---

## 1. Govaere et al. 2020 — Transcriptomic Profiling Across the NAFLD Spectrum

**Citation:** Govaere O et al. Sci Transl Med 2020;12(572):eaba4448
**DOI:** 10.1126/scitranslmed.aba4448

**Samples:** Discovery n=206 (RNA-seq, European Liver Fibrosis Study / Newcastle / Bern); Replication n=175 (NanoString nCounter); Biomarker n=305 (serum ELF score)

**Methods:** Bulk RNA-seq (Illumina TruSeq, STAR alignment, featureCounts, DESeq2). NanoString validation on independent cohort. WGCNA for coexpression modules. ROC/AUC for diagnostic performance. Logistic regression classifiers. IHC for protein validation (THY1, CD34, LYVE1). Pathway enrichment (Reactome, KEGG).

**Key Findings:**
- 25-gene progression signature discriminating disease stages
- AUROCs: NAS>=4 = 0.948, F>=2 = 0.916, F3-F4 = 0.930
- Disease progression driven by: ECM remodeling, complement activation, metabolic reprogramming
- THY1 (CD90) marks portal myofibroblasts; IHC confirms protein-level upregulation
- LYVE1 downregulation = sinusoidal endothelial dedifferentiation marker
- Serum ELF score correlates with transcriptomic disease stage (r=0.69)
- WGCNA: 25 modules; turquoise module (ECM/fibrosis) most correlated with NAS and fibrosis
- Identified CD34+ neovascularization increasing with fibrosis
- Zonation shift: pericentral gene loss, periportal gene gain with progression
- Ductular reaction genes (KRT7, KRT19, SOX9) upregulated in advanced fibrosis

**Datasets:** GSE135251 (discovery RNA-seq, n=206)

**Limitations:** Single-ethnicity European cohort. Bulk RNA-seq only (no cell-type resolution). Cross-sectional design. NanoString replication platform less granular than RNA-seq. No causal inference or genetic integration.

**Overlap with Our Work:** GSE135251 is one of our ~10 discovery cohorts. Their 25-gene signature is a direct comparator for our mega-analysis DEGs. Their WGCNA modules can be compared against our pathway results. CRITICAL: GSE135251 also appears in Govaere 2023 proteomics, creating partial circularity in our validation design.

---

## 2. Govaere et al. 2023 — Proteo-Transcriptomic Profiling of MASLD

**Citation:** Govaere O et al. Nat Metab 2023;5:572-578 (approximate)
**DOI:** (Nature Metabolism 2023)

**Samples:** n=306 patients; SomaScan 4,730 proteins + matched transcriptomics. Subset overlaps with GSE135251.

**Methods:** SomaScan v4 proteomics (aptamer-based, 4,730 analytes). Integrated proteo-transcriptomic analysis. Elastic net regression for biomarker panel. ROC analysis. NMF (non-negative matrix factorization) for molecular subtyping. Pathway enrichment. Cox proportional hazards for outcomes.

**Key Findings:**
- 31-marker proteo-transcriptomic signature for NASH diagnosis
- Diagnostic model (ADAMTSL2 + AKR1B10 + CFHR4 + TREM2 + BMI + T2DM): AUC = 0.878
- NMF molecular subtypes: 3 distinct MASLD subtypes with different progression trajectories
- Circulating protein-tissue mRNA concordance varies by protein (not uniform)
- TREM2 (soluble) elevated in advanced fibrosis — macrophage activation marker
- AKR1B10 — hepatocyte stress marker, strong RNA-protein correlation
- Fibrosis-associated proteins include collagens, ECM glycoproteins, complement components
- Subtype 1 = metabolic/lipotoxic; Subtype 2 = inflammatory/fibrotic; Subtype 3 = mixed
- Pathway-level: oxidative phosphorylation downregulated, ECM/complement upregulated in severe disease

**Datasets:** GSE276114 (SomaScan, n=306); overlaps with GSE135251

**Limitations:** SomaScan measures circulating (plasma) proteins, NOT tissue-level abundance — concordance with tissue RNA is indirect. Aptamer-based platform has known binding specificity issues for some analytes. Cohort partially overlaps GSE135251 discovery set.

**Overlap with Our Work:** We use GSE276114 SomaScan as a proteomics validation layer. The partial overlap with GSE135251 creates circularity risk — we need leave-one-out sensitivity analysis excluding GSE135251. Their NMF subtypes are a comparator for any subtyping we do. SomaScan validates circulating not tissue protein levels — we must frame this distinction clearly.

---

## 3. Li et al. 2025 — Single-Cell Multi-Omics Atlas of MASLD

**Citation:** Li T et al. Nat Genet 2025 (advance online)
**DOI:** (Nature Genetics 2025)

**Samples:** 61 human livers (healthy + MASLD spectrum). 540,216 scRNA-seq cells. 47,864 Visium spatial transcriptomics spots. 841,534 MALDI-MSI (mass spectrometry imaging) points.

**Methods:** scRNA-seq (10x Chromium). Visium spatial transcriptomics. MALDI-MSI for spatial metabolomics/lipidomics. CellChat for ligand-receptor interactions. SCENIC+ for GRN inference. Pseudotime trajectory analysis. Integration of spatial metabolite data with transcriptomics. Deconvolution of spatial spots.

**Key Findings:**
- MITF identified as master transcription factor of lipid-associated macrophages (LAMs)
- RSPO3-LGR6 axis: novel profibrotic signaling between mesenchymal cells and hepatocytes
- PLA2G7: key enzyme in phospholipid metabolism, spatially co-localizes with steatotic regions
- LAM trajectory: resident Kupffer cells → intermediate → TREM2+GPNMB+ LAMs
- Spatial metabolomics: ceramides and diacylglycerols enriched in steatotic zones
- Hepatocyte zonation disrupted: zone 1 (periportal) genes expand, zone 3 (pericentral) genes contract
- Stellate cell activation trajectory with COL1A1/COL3A1 upregulation
- Cholangiocyte/ductular reaction expansion in advanced fibrosis
- Novel cell-cell communication: HSC-hepatocyte, LAM-HSC, endothelial-HSC crosstalk networks
- Identified 12 cell types with distinct spatial niches in MASLD progression

**Datasets:** Proprietary (61 livers, not yet in GEO at time of publication)

**Limitations:** Single-center cohort. No bulk RNA-seq integration for broader validation. No causal inference (no GWAS/TWAS integration). No cross-species comparison. No drug repurposing analysis. MALDI-MSI spatial resolution lower than Visium for some applications. Cost-prohibitive to replicate.

**Overlap with Our Work:** Their LAM/TREM2 findings intersect our deconvolution-resolved macrophage signatures. Hepatocyte zonation disruption is consistent with our bulk-level findings. Their spatial data could contextually validate our DEGs but they don't share data publicly yet. No methodological overlap (they are single-cell atlas; we are bulk mega-analysis with multi-evidence convergence).

---

## 4. Feng et al. 2026 — MegaMASLD Multi-Omics Portal

**Citation:** Feng R et al. bioRxiv 2026
**DOI:** 10.64898/2026.02.19.706502

**Samples:** 29 datasets, 2,640 samples (mixed RNA-seq + microarray), 26 matched phenotypes. Human bulk transcriptomics only.

**Methods:** 4 modules: (A) DE meta-analysis (limma + random-effects meta-analysis across stages), (B) correlation meta-analysis (gene-phenotype correlations across cohorts), (C) SMR+HEIDI and TWAS for causal inference, (D) in vitro validation (primary hepatocytes + HepG2 + THP-1 with FFA/LPS/TGFb). 6-layer binary scoring system (DE + correlation + SMR + TWAS + liver-specific eQTL + essentiality). Web portal: masldportal.net.

**Key Findings:**
- 391 early-stage DEGs, 697 all-stage DEGs, 315 late-stage DEGs
- 27-gene core progression signature (present across all stages)
- MASLD GWAS yielded 0 significant SMR/TWAS genes (underpowered GWAS) — matches our finding
- 6-layer binary scoring: EFHD1 ranked #1 overall
- MLIP experimentally validated in vitro (knockdown reduces lipid accumulation in hepatocytes)
- Correlation meta-analysis: genes correlated with NAS, fibrosis, steatosis, ballooning, lobular inflammation
- Stage-specific gene programs: lipid metabolism (early) → inflammation (mid) → ECM remodeling (late)
- Portal features: gene-level query, cross-dataset visualization, downloadable results

**Datasets:** 29 datasets including GSE135251, GSE130970, GSE167523, GSE162694, GSE126848, GSE89632, and others (mixed platforms)

**Limitations:** Mixed RNA-seq + microarray without systematic platform-effect modeling (batch effects). Binary scoring loses quantitative information (a gene barely passing 4 layers scores same as one strongly passing). No cell-type deconvolution. No spatial or proteomic validation. No mouse/cross-species data. In vitro validation only (no in vivo). Preprint, not peer-reviewed.

**Overlap with Our Work:** HIGHEST direct competitor. Shares ~7 of our ~10 cohorts. Their 0 SMR/TWAS result mirrors ours. Key differentiators for us: (1) we do mega-analysis (joint mixed model via dream) not just meta-analysis, (2) we have cell-type deconvolution, (3) we have cross-species mouse benchmarking, (4) we have spatial/proteomic/epigenomic validation layers, (5) we have pharmacotranscriptomics. Their 6-layer binary scoring vs our multi-evidence convergence framework is a direct head-to-head comparison. Their portal (masldportal.net) is live; ours is pending.

---

## 5. Piras & DiStefano 2024 — GeneMeta Random-Effects Meta-Analysis

**Citation:** Piras IS, DiStefano JK. Life Sci Alliance 2024;7(11):e202402791
**DOI:** 10.26508/lsa.202402791

**Samples:** 10 datasets, 1,058 samples. Human bulk transcriptomics.

**Methods:** GeneMeta algorithm (random-effects model for combining effect sizes across studies). MAGMA gene-set analysis across 4 GWAS (Ghodsian, Sveinbjornsson, Anstee, Namjou). TWAS with FUSION + GTEx v8 liver tissue. MEGENA (Multiscale Embedded Gene Co-Expression Network Analysis) for coexpression modules. Bayesian causal network inference. ssGSEA for pathway activity.

**Key Findings:**
- DEGs: MASLD vs CTL = 685; MASH vs CTL = 1,870; MASH vs MASLD = 3,284; 173-gene overlap across all comparisons
- MAGMA: 141 significant gene-level signals across 4 GWAS, 127 unique genes
- TWAS: 84 genes significant in liver tissue
- 8 genes overlap between TWAS-meta and transcriptomic meta-analysis
- MEGENA: 103 coexpression modules; key hubs identified per module
- Bayesian network: directional causal relationships among top genes
- Top genes: HSD17B13 (genetic + transcriptomic evidence), PNPLA3, TM6SF2
- Stage-specific patterns: immune/inflammatory genes dominate MASH vs MASLD contrast
- ssGSEA: oxidative phosphorylation down, complement/coagulation up in MASH

**Datasets:** 10 datasets including GSE135251, GSE130970, GSE126848, GSE167523, GSE162694, others

**Limitations:** GeneMeta is less powerful than joint mixed-model mega-analysis (doesn't model sample-level covariates). No cell-type deconvolution. No experimental validation (purely computational). No spatial, proteomic, or epigenomic data. No cross-species analysis. No drug repurposing. Coexpression modules not validated externally.

**Overlap with Our Work:** Shares ~6-7 of our cohorts. Their GeneMeta approach is a weaker version of meta-analysis vs our dream mega-analysis. Their TWAS finding (84 genes) contrasts with our 0 FDR-significant result — methodological differences in TWAS pipeline worth investigating. Their 8-gene TWAS-meta overlap is a comparator for our causal inference module. We surpass them on every axis except Bayesian causal network (which we don't do).

---

## 6. Teufel/Raverdy et al. 2024 — SteatoSITE Endotyping

**Citation:** Teufel A, Raverdy V et al. Nat Med 2024
**DOI:** (Nature Medicine 2024)

**Samples:** ABOS cohort: n=1,389 discovery (clinical phenotyping). Liver transcriptomics n=831. Plasma metabolomics n=1,322. External validation in multiple cohorts.

**Methods:** PAM (Partitioning Around Medoids) clustering on 6 clinical variables (BMI, ALT, GGT, HbA1c, TG, HOMA-IR). Consensus clustering stability. Transcriptomic profiling of cluster-specific liver signatures. Plasma metabolomics (targeted). Cox proportional hazards for long-term outcomes. Kaplan-Meier for event-free survival. Web tool: https://ulr-metrics.univ-lille.fr/masldclusters/

**Key Findings:**
- 3 PAM clusters → collapsed into 2 clinically meaningful MASLD endotypes:
  - Cardiometabolic endotype: dysglycemia, high TG, CVD risk, HR 6.82 for T2DM, older
  - Liver-specific endotype: younger, PNPLA3-enriched, high ALT, more hepatocellular injury
- Transcriptomic differences between endotypes: cardiometabolic = insulin signaling/lipid metabolism genes; liver-specific = ECM/inflammatory genes
- Metabolomic validation: distinct metabolite profiles per endotype (BCAAs, ceramides)
- Long-term outcomes differ: cardiometabolic = cardiovascular events; liver-specific = liver-related events
- Endotype assignment reproducible across external cohorts
- Web tool allows clinician input of 6 variables → endotype classification

**Datasets:** ABOS cohort (Lille, France); external validation cohorts not fully specified

**Limitations:** Clustering on only 6 clinical variables (not transcriptomic-driven). PAM clustering sensitive to variable selection. ABOS is a bariatric surgery cohort (selection bias toward severe obesity). Not directly comparable to non-surgical MASLD populations. Endotypes are clinical phenotypes, not molecular subtypes. Transcriptomics is secondary validation, not primary discovery.

**Overlap with Our Work:** Their endotyping is complementary but orthogonal — they cluster on clinical variables, we characterize molecular programs. Their liver-specific endotype's transcriptomic signature could be compared against our DEG results. Their ABOS transcriptomic data (n=831) is not publicly available as a standard GEO dataset. Competitive threat is low; they occupy a clinical stratification niche, not a molecular atlas niche.

---

## 7. Vacca et al. 2024 — Mouse Model Benchmarking (LITMUS)

**Citation:** Vacca M et al. Nat Metab 2024;6(6):1178–1196
**DOI:** 10.1038/s42255-024-01043-6 (PMID 38867022; open-access PMC11199145; CC BY — supplementary tables downloaded to `data/external/vacca_2024/` and reused under attribution)

**Samples:** 598 animals across 39 murine MASLD models. 15 diet types. Multi-site (LITMUS consortium). Phenotyping: histology + biochemistry + transcriptomics.

**Methods:** MASLD Human Proximity Score (MHPS) = composite of **PHPS (Phenotype: body weight, lipids, liver enzymes)** + **HHPS (Histopathology: centralized steatosis/inflammation/ballooning/fibrosis scoring)** + **DHPS (Disease transcriptome: DSEA — a GSEA/gene2drug-style enrichment of the human-MASH DEG + KEGG signature in each model's ranked transcriptome)**. [NOTE: earlier internal notes mislabeled DHPS as "Dyslipidemia/metabolic" — corrected 2026-06-21; DHPS is the transcriptomic score.] Each sub-score is computed for two arms (metabolic-relevance and MASH-fibrosis). PLSR for study-design parameter analysis. Human benchmark = UCAM/VCU + EPoS biopsy cohorts (EPoS ≈ Govaere GSE135251).

**Key Findings:**
- Western diets (WD) closest to human MASLD transcriptomically and phenotypically
- Choline-deficient diets: fast fibrosis but variable metabolic translatability — **plain MCD (6J) has the LOWEST metabolic DHPS (0.23 single model; ~0.55 group mean)**, whereas CDAHFD actually scores HIGH metabolic DHPS (~0.77); both are non-physiological acute-injury models (lean, no obesity/insulin resistance)
- Top-ranked models: WD0.2% cholesterol, WD2% cholesterol, GAN-C2
- MHPS scoring reveals no single model recapitulates all human features
- Strain matters: C57BL/6J vs C57BL/6N show different susceptibilities
- Duration matters: 12-16 week WD insufficient for fibrosis; 24+ weeks needed
- PLSR: diet composition (fat%, cholesterol%, sugar type) and duration are dominant design parameters
- Fructose supplementation in water enhances metabolic phenotype but not fibrosis
- ob/ob and db/db genetic models: good metabolic phenotype, poor fibrosis
- Gene-level concordance: ~40-60% overlap between top mouse models and human DEGs depending on threshold

**Datasets:** Not deposited in GEO as unified dataset. Individual studies referenced.

**Limitations:** MHPS weights are somewhat arbitrary (equal weighting of PHPS/HHPS/DHPS). Transcriptomic concordance computed on bulk liver (no cell-type resolution). Not all 39 models had RNA-seq data. Inter-lab variability acknowledged but not fully resolved. No single model scores >0.7 on all three sub-scores simultaneously.

**Overlap with Our Work — REFRAMED 2026-06-21 (competitor → external benchmark we ADOPT + extend):** Rather than re-deriving a mouse-model ranking and "competing," we **cite LITMUS's DHPS as the authoritative model-translatability benchmark** and show concordance with our existing cross-species results (canonical scripts `Analysis/Cross_Species_Concordance/scripts/11_vacca_q1_fibrotic_validation.R` + `12_vacca_q2_contrast_matrix.R`; diagnostic ladder `08/09/10`; figure `figS06_vacca_benchmark`):
- **Gene-level — CONTRAST-MATCHED (Q2, script 12):** LITMUS's 951 signature is itself stage-partitioned (271 Early / 526 Progression / 154 All). Matching the contrast resolves the modest raw overlap: our **progression** human signature (NAFL→NASH, LVQW) matches their **progression** bin best (**Fisher OR 31.7**, p=4e-165), while disease-vs-control matches the **early** bin only weakly (**OR 1.7** — a near-null specificity control), and the progression signature is *depleted* in the early bin (OR 0.46). A **threshold-free, NON-CIRCULAR anchor** — genes concordant in ≥1 mouse diet (no human-side cutoff) — is **4.3× enriched in their 951-signature** (Fisher p=1e-72, n=407), so the 172-gene raw overlap was a *contrast mismatch + cutoff*, not a biological ceiling. Separately, our 1,355-gene Conserved_Core is enriched **OR 2.8×** (3.3× for the progression subset; biotype-matched-null p<0.001) with **99% same-direction agreement**, and **1,151 conserved-core genes lie beyond their signature (our per-gene extension)**.
- **Model-level — FIBROTIC-arm validation (Q1, script 11; DHPS is two-armed + PATHWAY-level):** our initial GENE-level comparison anti-correlated with DHPS (Spearman −1). The diagnostic ladder (08/09/10) showed why: fgsea+signature on their own per-model data reproduces DHPS at only Spearman 0.06 gene-level but **0.73 at the KEGG-PATHWAY level** — DHPS is a pathway-program metric (consistent with their DSEA/gep2pep methods). DHPS is also **two-armed** (metabolic + fibrotic). Our pathway-PEP reproduces the **FIBROTIC arm across LITMUS's own 33 mouse models (Spearman 0.62 full / 0.54 shared-KEGG, p=1e-4)**; the **metabolic arm does NOT reproduce (0.29, NS)** and is **deferred to LITMUS**. On that validated fibrotic axis our **Western diet (FPC) is the most human-proximal of our 4 diets (PEP 0.92)** — reproducing "Western closest to human"; MCD's mid-rank is genuine fibrotic proximity (LITMUS CDD f_DHPS~0.75), its metabolic demerit being the un-reproduced arm. (FPC's point estimate is leverage-sensitive — few-pathway-driven, bootstrap 95% CI [0.44, 0.98] — so we present **consistency with, and defer the precise ranking to, LITMUS**.) This corrected a false claim in `fig4g`/`crossspecies_conservation.R` (old gene-level GSEA-NES panel mislabeled "the DHPS-style metric of Vacca" — replaced with the validated fibrotic pathway-PEP).
- ⚠ **Shared-cohort caveat & circularity hygiene:** their human benchmark (EPoS) ≈ our Govaere GSE135251 cohort, so gene-level overlap is a *consistency* check; the OR 31.7 progression match **survives dropping Govaere** (refit LOO, OR 26.9→27.7) and 1,000× label shuffle (p<0.001), and the non-circular anchor is mouse-anchored. We **dropped all circular human-vs-human metrics** (RRHO, signed Spearman ~0.60, 72.6% direction concordance) — they are signed by our own human logFC against LITMUS's (largely shared) human logFC.
- We still extend beyond them with **cell-type-resolved** concordance (mouse-bulk deconvolution) and **per-gene** concordance scores they lack. Net: LITMUS strengthens (not threatens) our cross-species story.

---

## 8. Wang et al. 2021 — scRNA-seq of Liver Non-Parenchymal Cells in Fibrosis

**Citation:** Wang ZY et al. Sci Rep 2021;11:13381
**DOI:** 10.1038/s41598-021-92564-1

**Samples:** 6 healthy human liver donors, 17,810 NPCs (non-parenchymal cells) by scRNA-seq. Mouse in vivo (CCl4 fibrosis model) and in vitro (primary HSC activation) comparators.

**Methods:** scRNA-seq (10x Chromium, Seurat v3). SCENIC for regulon/GRN analysis. Cross-species comparison (human vs mouse). Cross-organ comparison (liver vs pancreas stellate cells). Ligand-receptor interaction analysis. GO/KEGG enrichment. Pseudotime/trajectory analysis for HSC activation.

**Key Findings:**
- RUNX1 and CREB3L1 identified as master regulators of activated HSC (aHSC) program
- 61 translationally conserved fibrosis genes (overlap of human + mouse in vivo + in vitro)
- HSC activation trajectory: quiescent → intermediate → fully activated, with bifurcation point
- Cross-organ: liver and pancreatic stellate cells share ~40% of activation program
- Key aHSC markers: LARP6 (collagen chaperone), AEBP1 (ECM crosslinker), PRRX1 (TF)
- Novel ligand-receptor interactions: ITGAV-LAMC1, TNFRSF11B-VWF, NOTCH2-DLL4
- Quiescent HSC markers: RGS5, IGFBP7, DCN; distinguish from other mesenchymal cells
- Kupffer cell and endothelial cell subtypes also characterized
- Liver-specific stellate cell genes vs shared fibroblast genes delineated

**Datasets:** Not deposited as unified GEO. Individual references to existing atlases.

**Limitations:** Healthy donors only (no MASLD/fibrosis human tissue). NPC-focused, no hepatocytes sequenced. Small sample (n=6 donors). Mouse model is CCl4 (toxic, not metabolic). Cross-species comparison at gene level only, not pathway or regulatory level. SCENIC regulon analysis sensitive to parameter choices.

**Overlap with Our Work:** Their HSC activation genes and regulons are comparators for our deconvolution-resolved fibrosis signatures. The 61 conserved fibrosis genes can be checked against our cross-species concordance analysis. Their RUNX1/CREB3L1 TF findings could enrich our regulatory inference if we include GRN analysis. Limited direct competition — they focus on NPCs in healthy liver, we focus on disease progression across the full spectrum.

---

## 9. Suppli et al. 2019 — Hepatic Transcriptome in NAFLD Spectrum

**Citation:** Suppli MP et al. Am J Physiol Gastrointest Liver Physiol 2019;317(2):G157-G167
**DOI:** 10.1152/ajpgi.00358.2018

**Samples:** n=57 total: 14 normal-weight (NW), 12 obese-no-NAFLD, 15 NAFL, 16 NASH. Paired liver biopsy RNA-seq + histomorphometry + IHC.

**Methods:** Bulk RNA-seq (Illumina HiSeq). DESeq2 for differential expression. IHC for Sonic Hedgehog (SHH) pathway activation. Histomorphometry (digital image analysis of steatosis, inflammation, fibrosis area). WGCNA. GSEA. Stratification of NASH patients by SHH+ hepatocyte status.

**Key Findings:**
- Obesity alone does NOT significantly alter liver transcriptome (only 55 DEGs obese vs NW)
- NASH vs NAFL: moderate transcriptomic shift (282 DEGs) driven by ECM and inflammation
- SHH-positive hepatocyte staining stratifies NASH: SHH+ NASH vs NAFL = 282 DEGs; SHH- NASH vs NAFL = only 2 DEGs
- Implication: Hedgehog pathway activation marks a transcriptomically distinct NASH subtype
- Top NASH vs NAFL genes: MMP2, FMNL3, OTOA, CPXM1, MMP24, GPC3
- SHH+ genes enrich for wound healing, ECM organization, angiogenesis
- Histomorphometry: steatosis area correlates poorly with NAS steatosis grade (Spearman r=0.52)
- Perisinusoidal fibrosis area correlates with overall fibrosis stage
- WGCNA: fibrosis module distinct from steatosis module

**Datasets:** GSE126848 (n=57)

**Limitations:** Small cohort (n=57). Male-only controls (NW and obese groups). Single-center (Copenhagen). Low disease severity (mostly F0-F2). No advanced fibrosis representation. No validation cohort. SHH IHC semi-quantitative. No causal inference or genetic data.

**Overlap with Our Work:** GSE126848 may be one of our discovery cohorts. Their finding that obesity alone barely changes liver transcriptome supports our disease-specific mega-analysis approach. SHH stratification is a potential subtyping comparator. Their SHH+ NASH gene set could be tested for enrichment in our progression-associated DEGs. Low competitive threat given small size and narrow scope.

---

## 10. Kuchenhoff et al. 2026 — Cross-Organ Fibrosis scRNA-seq Meta-Analysis

**Citation:** Kuchenhoff L et al. bioRxiv 2026
**DOI:** 10.64898/2026.03.09.709232

**Samples:** 20 studies, 4 organs (liver, heart, lung, kidney), 953 samples, >5 million cells total. scRNA-seq meta-analysis.

**Methods:** scRNA-seq integration across organs using Harmony/scVI. ScDist (distance-based metric) for cross-study transcriptomic shifts. Multicellular factor analysis. Cell-type annotation harmonization. Cross-organ consensus gene identification. AUROC for cross-study agreement. Pathway enrichment per organ and cell type. Web resource: https://organfibrosis.saezlab.org/

**Key Findings:**
- Liver has FEWEST organ-consensus fibrosis genes (avg 36.8 upregulated per cell type vs heart 353.8)
- Liver shows LOWEST cross-study agreement (median AUROC 0.60 vs heart 0.75)
- Disease fibroblasts/myofibroblasts converge across all 4 organs: shared activation program
- PPARG consistently downregulated in fibroblasts across all organs
- Cross-organ mesenchymal consensus genes: SMOC2, LOXL1, FGF14, MOXD1
- Macrophage response partially shared (TREM2, SPP1, GPNMB) but liver macrophages are more heterogeneous
- Endothelial cells show organ-specific fibrosis responses (minimal cross-organ consensus)
- Hepatocytes are the most organ-specific cell type (no meaningful cross-organ comparison)
- Multicellular factor analysis reveals coordinated cell-type programs in fibrosis

**Datasets:** 20 published scRNA-seq datasets across organs (references listed in paper)

**Limitations:** Preprint, not peer-reviewed. Cross-organ comparison forced cell-type harmonization that may lose organ-specific subtypes. Liver's low consensus could reflect biological heterogeneity OR methodological challenges (different protocols, lower cell quality). ScDist metric not yet widely validated. No bulk RNA-seq integration. No causal inference.

**Overlap with Our Work:** Their finding that liver has the lowest cross-study agreement contextualizes why mega-analysis (pooling raw data with harmonized covariates) may outperform meta-analysis for liver. Their cross-organ fibroblast genes (SMOC2, LOXL1) could appear in our fibrosis DEGs. Their TREM2/SPP1/GPNMB macrophage findings intersect our deconvolution results. Low direct competition — they are cross-organ focused, not MASLD-specific. Their liver heterogeneity finding could be cited to justify our mega-analysis approach.

---

## 11. Tzouanas et al. 2026 — Longitudinal Multi-Omic MASLD Atlas (Cell)

**Citation:** Tzouanas CN et al. Cell 2026;189(3):435-460
**DOI:** 10.1016/j.cell.2025.11.031

**Samples:** Mouse: longitudinal HFD model at 6/12/15 months, tandem snRNA-seq + snATAC-seq (exact n per timepoint in paper). Human: Xenium spatial transcriptomics, n=18 MASH cirrhosis livers, 723,466 cells. Validation across 6 human cohorts (Govaere, Hoang, Pantano, Filliol, Niu, TCGA-LIHC).

**Methods:** snRNA-seq + snATAC-seq (10x Multiome). MATCHA algorithm (novel): TF-enhancer-gene regulatory program triads. Xenium spatial transcriptomics (10x). Hepatocyte stress program classification (4 categories). HMGCS2 LiverKO mouse model (CRISPR, AAV-Cre). In vitro perturbation (SOX4, RELB overexpression in hepatocytes). Kaplan-Meier survival. Pseudotime trajectory. ChromVAR for TF motif activity. Validation in 6 human bulk RNA-seq cohorts.

**Key Findings:**
- 4 hepatocyte stress response programs defined by temporal dynamics:
  - Longitudinal Increase (LI): progressive upregulation (inflammation, ECM)
  - Longitudinal Decrease (LD): progressive downregulation (metabolic function)
  - Sustained Upregulation (SU): early onset, maintained
  - Sustained Downregulation (SD): early loss, maintained
- HMGCS2 is in LD program; HMGCS2 LiverKO + HFD = increased tumorigenesis (HCC)
- HMGCS2 expression decreases across human MASLD cohorts (validated in 6 datasets)
- MATCHA identifies SOX4 and RELB as disease TFs; validated by overexpression → stress response
- Xenium spatial: CD9+TREM2+ macrophages co-localize with stress-high hepatocytes in fibrotic niches
- snATAC-seq: progressive chromatin remodeling at metabolic gene loci
- Enhancer-gene links disrupted in disease: metabolic enhancers close, inflammatory enhancers open
- TF motif activity: HNF4A/CEBPA decrease, AP-1/NFkB increase with progression
- Cross-cohort validation: stress programs replicate in all 6 human datasets

**Datasets:** Mouse data (likely to be deposited, GEO TBD). Xenium data (proprietary/TBD). Validated in GSE135251, GSE130970, others.

**Limitations:** Mouse model is standard HFD (not the highest-fidelity MASLD model per Vacca et al.). Xenium panel limited to ~300 genes (not transcriptome-wide). Human validation is correlation-based (no prospective validation). MATCHA algorithm is novel and not yet independently replicated. No drug repurposing analysis. No proteomic validation. Mouse-to-human translation limited to gene-level, not regulatory-level.

**Overlap with Our Work:** HIGH competitive relevance. Published in Cell — sets the bar for multi-omic MASLD atlases. Their hepatocyte stress programs are comparators for our DEG trajectory results. Their HMGCS2 finding intersects our hepatocyte-intrinsic signal. Their use of Govaere and Hoang cohorts for validation overlaps our discovery cohorts. Key differentiation: (1) they are mouse-first with human validation; we are human-first with mouse benchmarking, (2) they focus on epigenomic regulation; we focus on convergence across evidence types, (3) they have no pharmacotranscriptomics, (4) they have no causal inference / GWAS integration. Their MATCHA algorithm and in vivo validation (HMGCS2 KO) are strengths we cannot match without wet lab work.

---

## 12. Chen et al. 2023 — Smoothened/Hedgehog in Hepatocyte Metabolism

**Citation:** Chen J et al. Cell Mol Gastroenterol Hepatol 2023;15(5):1025-1046
**DOI:** 10.1016/j.jcmgh.2022.12.008

**Samples:** Mouse: hepatocyte-specific Smo conditional KO (SmoΔHep) + controls, chow-fed. scRNA-seq: 9,207 liver cells. Human validation: Duke cohort n=368, GSE135251 n=206, GSE167523.

**Methods:** Conditional Cre-lox knockout (Alb-Cre x Smo-flox). scRNA-seq (10x Chromium). Insulin/glucose tolerance tests. Hepatocyte zonation analysis (pericentral-periportal axis). Overlap with aging gene signatures. Oil Red O / H&E histology. Western blot. qRT-PCR. ROC analysis for human biomarker potential. Pathway analysis (IPA).

**Key Findings:**
- Smo KO in hepatocytes → rapid fatty liver + hyperinsulinemia + insulin resistance on normal chow
- Mechanism: loss of Smo → mTORC1-SREBP1c axis activation → de novo lipogenesis
- Smo-active hepatocytes localize to central-to-midzonal layers (zonation L3-L7)
- Smo-/- hepatocytes show mitochondrial dysfunction, DNA damage, premature aging phenotype
- 105-gene overlap between Smo-KO and aged hepatocyte signatures (P < 2.1e-122)
- Human validation: SMO expression declines with fibrosis stage across 3 cohorts
- SMO AUC for F3/F4 vs F0/F1: 0.75 (moderate discriminative ability)
- 12 Smo+ signature genes declining in human NAFLD: MAT1A, ALDH6A1, GPT2, PPM1K, HAL, RORA, ARRDC3, ACADSB, KLF9, KLF15, RETREG1, ETNK2
- Smo-active program = hepatocyte metabolic homeostasis; its loss = MASLD initiation
- Hedgehog pathway drugs (vismodegib etc.) may have hepatotoxic implications

**Datasets:** Duke cohort (not public), GSE135251, GSE167523

**Limitations:** Mouse model: constitutive KO (not inducible, not disease-model), chow-fed (no dietary MASLD model). scRNA-seq sample size small (1 KO + 1 WT). Human validation is correlative (SMO decline could be consequence not cause). AUC 0.75 is moderate. Single pathway focus. No drug repurposing beyond noting Hh pathway drug risk.

**Overlap with Our Work:** Uses GSE135251 and GSE167523 (our discovery cohorts) for validation. Their 12 Smo+ signature genes are testable against our DEG results. SMO / Hedgehog pathway may appear in our pathway enrichment. Minimal competitive threat — they are a mechanistic single-pathway study, not an atlas. Their aging-MASLD connection could enrich our discussion if we find overlapping genes.

---

## 13. Kozumi et al. 2021 — THBS2/TSP-2 as NASH Biomarker

**Citation:** Kozumi K et al. Hepatology 2021;74(5):2452-2468
**DOI:** 10.1002/hep.31995

**Samples:** Transcriptome cohort n=98 (51 NAFL, 47 NASH). Serum biomarker cohort n=213. Follow-up cohort n=164. Japanese single-center.

**Methods:** Microarray transcriptomics (Agilent). Bioinformatic prediction of secreted proteins (SignalP, TMHMM). ELISA for serum TSP-2 measurement. ROC analysis. Kaplan-Meier survival for liver events. Multivariate logistic regression. Spearman correlations between THBS2 mRNA and collagen genes. IHC for tissue localization.

**Key Findings:**
- 137 upregulated genes in NASH vs NAFL; 29 predicted secreted proteins among them
- THBS2 (thrombospondin-2 / TSP-2) identified as top secreted protein biomarker
- THBS2 mRNA AUROC: 0.915 (NASH vs NAFL), 0.957 (F3-F4 vs F0-F2)
- Serum TSP-2 AUROC: 0.776 (NASH), 0.856 (F3-F4) — lower than mRNA, expected for circulating biomarker
- 8 fibrosis-associated secreted genes: THBS2, LUM, LAMC3, LAMA2, AKR1B10, COL4A4, A2M, C7
- THBS2 strongly correlated with COL1A1 (r=0.899), COL1A2 (r=0.912) — suggests HSC-derived
- IHC: THBS2 protein localizes to fibrotic septa and sinusoidal spaces
- Kaplan-Meier: high serum TSP-2 associated with liver cancer and hepatic decompensation
- TSP-2 outperforms FIB-4 and APRI for advanced fibrosis detection in their cohort
- Validated directionally in GSE135251

**Datasets:** Proprietary Japanese cohort (microarray). Validated in GSE135251.

**Limitations:** Single-center Japanese cohort (ethnic specificity unknown). Microarray platform (lower sensitivity than RNA-seq). Small transcriptome cohort (n=98). Serum AUCs lower than mRNA AUCs (0.776 vs 0.915). No multi-cohort meta-analysis. No mechanistic validation (correlation not causation). No cell-type resolution (THBS2 source cell type inferred, not proven).

**Overlap with Our Work:** THBS2 likely appears in our DEG results given its strong effect size. Uses GSE135251 for validation (our discovery cohort). Their 8 secreted protein list can be cross-referenced against our fibrosis DEGs. AKR1B10 appears in both this paper and Govaere 2023. Low competitive threat — narrow biomarker focus, not an atlas. Could be cited as supporting evidence if THBS2 appears in our convergence gene list.

---

## 14. Hoang et al. 2019 — Gene Expression Signatures Across NAFLD Spectrum

**Citation:** Hoang SA et al. Sci Rep 2019;9:12541
**DOI:** 10.1038/s41598-019-48746-5

**Samples:** VCU (Virginia Commonwealth University) single-center cohort. Full NAFLD spectrum: NAS 0-6, fibrosis F0-F4. Exact n not specified in abstract but full-spectrum coverage.

**Methods:** Bulk RNA-seq. Ordinal regression for continuous severity scores (gNAS = genomic NAS, gFib = genomic fibrosis score). DEGs by FDR 1% threshold. PPI network analysis (STRING + Louvain community detection). Lasso regression for gene-based disease prediction. Gene shaving for molecular subtyping (unsupervised). GSEA/pathway enrichment. Correlation with histological features.

**Key Findings:**
- DEGs: 2,970 for NAS, 1,656 for fibrosis (FDR 1%)
- Ordinal regression enables continuous (not binary) disease scoring — gNAS and gFib
- PPI network: NAS = 2,139 nodes / 18 communities; Fibrosis = 1,225 nodes / fewer communities
- Lasso regression: ~20 genes predict NAS (R² = 0.96) and fibrosis (R² = 0.94)
- Gene shaving molecular subtypes: distinct patient clusters with different leading-edge processes
  - Subtype A: inflammation/apoptosis dominant
  - Subtype B: cell stress/metabolism dominant
- Fibrosis network hubs: COL1A1, COL1A2, COL3A1, TIMP1, THBS2, LUM
- NAS network hubs: lipid metabolism, immune cell activation, oxidative stress genes
- gNAS and gFib scores are orthogonal — steatohepatitis and fibrosis have partially independent molecular programs
- Validation: internal cross-validation only

**Datasets:** GSE130970

**Limitations:** Single-center (VCU). No external validation cohort. No causal inference. No cell-type deconvolution. No proteomic or spatial data. No cross-species analysis. Ordinal regression assumes linear gene-severity relationship. Gene shaving molecular subtypes not replicated externally. PPI network analysis depends on STRING database completeness.

**Overlap with Our Work:** GSE130970 is one of our ~10 discovery cohorts. Their gNAS/gFib continuous scoring approach is a methodological comparator — our dream mega-analysis uses histological scores differently but addresses the same question. Their 2,970 NAS DEGs and 1,656 fibrosis DEGs are direct comparators for our DEG lists. Their finding that NAS and fibrosis are partially orthogonal molecular programs is important context for our analysis design. Hoang is used as a validation cohort in Tzouanas et al. 2026.

---

## Cross-Paper Summary Statistics

| Paper | Year | Journal | Samples | Primary Data Type | Has Deconvolution | Has Causal Inference | Has Cross-Species | Has Drug Repurposing | Has Experimental Validation |
|---|---|---|---|---|---|---|---|---|---|
| Govaere 2020 | 2020 | Sci Transl Med | 206+175+305 | Bulk RNA-seq | No | No | No | No | IHC only |
| Govaere 2023 | 2023 | Nat Metab | 306 | SomaScan proteomics | No | No | No | No | No |
| Li 2025 | 2025 | Nat Genet | 61 livers | scRNA-seq + spatial + MALDI | Spatial deconv | No | No | No | No |
| Feng 2026 | 2026 | bioRxiv | 2,640 | Bulk mixed (RNA-seq+array) | No | SMR/TWAS | No | No | In vitro |
| Piras 2024 | 2024 | Life Sci Alliance | 1,058 | Bulk RNA-seq | No | MAGMA/TWAS | No | No | No |
| SteatoSITE 2024 | 2024 | Nat Med | 1,389+831+1,322 | Clinical+transcriptomics+metabolomics | No | No | No | No | No |
| Vacca 2024 | 2024 | Nat Metab | 598 animals | Mouse bulk RNA-seq | No | No | Yes (mouse only) | No | No |
| Wang 2021 | 2021 | Sci Rep | 6 donors | scRNA-seq (NPCs) | No | No | Mouse+human | No | No |
| Suppli 2019 | 2019 | AJP-GI | 57 | Bulk RNA-seq | No | No | No | No | IHC (SHH) |
| Kuchenhoff 2026 | 2026 | bioRxiv | 953 samples | scRNA-seq (4 organs) | No | No | No | No | No |
| Tzouanas 2026 | 2026 | Cell | Mouse+18 human | snRNA+snATAC+Xenium | No | No | Mouse→human | No | In vivo (HMGCS2 KO) + in vitro |
| Chen 2023 | 2023 | CMGH | Mouse+368+206 human | scRNA-seq + bulk | No | No | Mouse→human | No | In vivo (Smo KO) |
| Kozumi 2021 | 2021 | Hepatology | 98+213+164 | Microarray + ELISA | No | No | No | No | Serum ELISA |
| Hoang 2019 | 2019 | Sci Rep | Full spectrum | Bulk RNA-seq | No | No | No | No | No |

## Competitive Gap Analysis (What NO Competitor Has Done)

1. **Mega-analysis with dream mixed model** across >10 cohorts (Feng does meta-analysis, not mega-analysis)
2. **Cell-type deconvolution of bulk RNA-seq** applied to mega-analysis DEGs (MuSiC/BayesPrism)
3. **Cell-type-resolved drug repurposing** via LINCS L1000 (no competitor does pharmacotranscriptomics)
4. **Multi-evidence convergence scoring** combining transcriptomic + genetic + epigenomic + spatial + proteomic + essentiality + cross-species (Feng has 6-layer binary; we have quantitative multi-evidence)
5. **Cross-species concordance at cell-type resolution** (Vacca does bulk-level only)
6. **Integrated causal inference + deconvolution + drug repurposing** in single framework

## Shared Dataset Overlap Matrix

| Dataset | Our Work | Govaere 2020 | Govaere 2023 | Feng 2026 | Piras 2024 | Hoang 2019 | Suppli 2019 | Chen 2023 | Kozumi 2021 | Tzouanas 2026 |
|---|---|---|---|---|---|---|---|---|---|---|
| GSE135251 | Discovery | Discovery | Overlap | Yes | Yes | - | - | Validation | Validation | Validation |
| GSE130970 | Discovery | - | - | Yes | Yes | Discovery | - | - | - | Validation |
| GSE167523 | Discovery | - | - | Yes | Yes | - | - | Validation | - | - |
| GSE126848 | Discovery? | - | - | Yes | Yes | - | Discovery | - | - | - |
| GSE162694 | Discovery? | - | - | Yes | Yes | - | - | - | - | - |

---

## Our Novel Contributions vs. All Competitors

This section maps each of our unique analyses to the competitors it surpasses, organized by analysis type. Use this to identify what to emphasize in manuscript framing, introduction, and discussion.

### 1. Scale & Methodology
| Our Capability | What We Have | Closest Competitor | Why We Win |
|---|---|---|---|
| **Pure RNA-seq mega-analysis** | 1,444 samples, 10 cohorts, dream mixed model with `(1\|dataset)` | Feng: 2,640 samples BUT mixed RNA-seq+microarray, metafor meta-analysis only | Dream pools raw data (sample-level covariates); meta-analysis only combines summary statistics. Kuchenhoff shows liver has lowest cross-study agreement (AUROC 0.60) — mega-analysis is the right approach |
| **LOO-CV validation** | 88.1% mean recovery, rho=0.959, Jaccard=0.822 | None | No competitor reports leave-one-cohort-out cross-validation of their DEG results |
| **Cross-study AUROC** | 0.895 (integrated) vs 0.735 (single-cohort) | SteatoSITE: single-cohort only | Demonstrates integration value quantitatively |

### 2. Cell-Type Resolution
| Our Capability | What We Have | Closest Competitor | Why We Win |
|---|---|---|---|
| **Deconvolution-adjusted bulk DE** | MuSiC + BayesPrism dual-method; 1,635 hepatocyte-intrinsic DEGs (~29.8% of primary DEGs; cross-method Jaccard=0.774) | None at population scale | Li (Nat Genet) has scRNA deconv but only 61 livers; no competitor does deconvolution of bulk mega-analysis DEGs |
| **Sex-stratified integrated analysis** | 3,080 sex-dimorphic DEGs via interaction model: 1,978 female-biased, 675 male-biased, 427 divergent (matched-dataset, 4 cohorts) | None | No competitor reports sex-stratified DE at multi-cohort scale in MASLD |
| **Disease progression (pseudotime)** | Consensus pseudotime (rho_fib=0.443), S1/S2 bifurcation (1,797 genes), CCC rewiring (4,800 pairs) | Tzouanas (Cell): longitudinal HFD mouse + Xenium | We have population-scale bulk (1,444 samples); they have 3 mouse timepoints + 18 human spatial samples. Complementary approaches |

### 3. Causal Inference
| Our Capability | What We Have | Closest Competitor | Why We Win |
|---|---|---|---|
| **Multi-source COLOC** | 538 SuSiE-COLOC genes at PP.H4.susie > 0.5 (322 > 0.8, 226 > 0.9; PolyFun UKBB EUR LD, Phase 10 swap 2026-05-06) across a 28-GWAS portfolio × Broadaway eQTL (N=1,183); ABF fallback recovers 502/241/168 across the full 18,975 eGenes | Feng: SMR/TWAS only (no COLOC); Piras: MAGMA/TWAS only | SMR ≠ colocalization; SuSiE-COLOC fine-maps both sides under ancestry-matched LD and resolves multi-signal loci |
| **Cell-type-resolved causal** | sc-eQTL COLOC across 4 cell types (hepatocyte, endothelial, cholangiocyte, stellate) from 312 donors | None | No competitor does cell-type-specific COLOC in MASLD |
| **Multi-ancestry validation** | 28-GWAS portfolio: 17 EUR (incl. FinnGen R12) + 5 BBJ EAS (ALT/AST/GGT + 2020 Cirrhosis/HCC) + 3 Pan-UKBB AFR + 3 Pan-UKBB CSA (ALT/AST/GGT), each fine-mapped against ancestry-matched LD | Seagle: multi-ancestry GWAS MR but no transcriptomics | We integrate a 4-ancestry GWAS portfolio with the transcriptomic atlas; EUR+EAS cross-ancestry validated targets (67 genes) are the highest-confidence tier |
| **GWAS-ATAC regulatory variants** | 7,766 credible set variants → 533 in scATAC peaks, 380 with motif disruption, 12 disease regulon TFs | None | No competitor maps fine-mapped variants to cell-type-specific chromatin accessibility |

### 4. Pharmacotranscriptomics & Translation
| Our Capability | What We Have | Closest Competitor | Why We Win |
|---|---|---|---|
| **LINCS L1000 drug repurposing** | 1,107 compounds annotated; MOA-level reversal scores | None | No competitor does pharmacotranscriptomic drug matching against MASLD signatures |
| **Network proximity** | 625 drugs with z-scored proximity to disease module (STRING PPI) | None | Complementary to reversal scoring |
| **Clinical drug validation** | 10 MASLD drugs × full evidence profile (4 Strong, 4 Moderate) | None | Systematic benchmark of atlas predictions against known MASLD therapeutics |
| **Proteomics validation** | 3 datasets (Govaere SomaScan, plasma, tissue-level) | Govaere 2023 has SomaScan but uses it for discovery, not validation of bulk atlas | We use proteomics as independent validation layer |

### 5. Multi-Evidence Convergence
| Our Capability | What We Have | Closest Competitor | Why We Win |
|---|---|---|---|
| **7-source evidence atlas** | 33,943 genes × 101 columns; 7 independent sources + convergence integration | Feng: 6-layer binary scoring | Quantitative vs binary; demonstrated source non-redundancy (3.6% DEG-COLOC overlap, OR = 1.67); Convergence score + network propagation |
| **Network convergence** | 123 genes with 7/7 source convergence | None | No competitor demonstrates multi-source convergence at this scale |
| **ncRNA integration** | 6,348 lncRNA DEGs, 4,674 synteny-conserved, 11 ceRNA triplets, atlas expanded to include ncRNAs | None | No competitor systematically analyzes ncRNAs in MASLD at population scale |

### 6. Prediction & Staging
| Our Capability | What We Have | Closest Competitor | Why We Win |
|---|---|---|---|
| **Disease stage classifier** | F>=3 AUROC 0.881 (LOCO-CV); concept bottleneck QWK 0.624; multi-architecture comparison | Govaere 2020: 25-gene signature (single cohort); Hoang: Lasso panels (single cohort) | LOCO-CV across 6+ cohorts vs single-cohort internal CV; concept bottleneck provides interpretability |
| **Cross-species concordance** | 10,932 genes, 1,108 Conserved_Core, 5 diet models, FPC model | Vacca: MHPS across 39 models but bulk-level only | We add cell-type-resolved concordance; per-gene concordance scores integrated into atlas |

### Key Competitor Findings to Validate in Our Atlas

These specific claims from competitors should be checked against our 1,444-sample atlas as validation or replication:

1. **HMGCS2 decline** with progression (Tzouanas, Cell 2026) — check in dream DEGs + progression pseudotime
2. **Hepatocyte stress programs LI/LD/SU/SD** (Tzouanas) — test program gene sets for enrichment in our transition-specific DEGs
3. **RUNX1/CREB3L1 as HSC regulators** (Wang 2021) — check in SCENIC+ regulon analysis
4. **SHH+ stratification** (Suppli 2019) — test SHH pathway genes in our NMF subtypes
5. **THBS2 as fibrosis biomarker** (Kozumi 2021) — check dream DEG status, progression trajectory, Conserved_Core membership
6. **Cross-organ consensus genes** SMOC2/LOXL1/FGF14/MOXD1 (Kuchenhoff 2026) — check in our fibrosis DEGs
7. **MITF as LAM master TF** (Li 2025) — check in our SCENIC+ regulons and scRNA analysis
8. **gNAS/gFib orthogonality** (Hoang 2019) — compare with our NAS vs fibrosis program separation
9. **EFHD1 as top target** (Feng 2026) — already checked: Bayesian rank 3088 (top 22%), mediocre
10. **MLIP validation** (Feng 2026) — already checked: COLOC for ALT/GGT but not NAFLD/NASH-specific
