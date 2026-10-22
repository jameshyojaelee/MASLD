> **Historical competitor/source dossier only (2026-08-10).** Source-paper
> facts may remain useful, but project comparisons, figures, targets, journals,
> and next-analysis suggestions are non-authoritative. Reverify source facts;
> use [`../../PAPER.md`](../../PAPER.md),
> [`../../STATUS.md`](../../STATUS.md), and
> [`../../RESULTS.md`](../../RESULTS.md) for current claims.

# Govaere et al. 2026 — Summary: Integrated multi-omics identifies distinct macrophage alterations during MASH progression

**Citation**: Boesch M, Anak S, El Abyad D, et al. (Govaere O, last author). *Nature Genetics*, published online 18 May 2026. https://doi.org/10.1038/s41588-026-02600-3

**Data access**: snRNA-seq: EGA accession EGAS50000000768. Bulk RNA-seq: GSE135251. CosMx + GeoMx: GSE312698 + Zenodo 10.5281/zenodo.18153344. Code: github.com/govaerelab/OG_MB-MASLD-LIVER-2025.

---

## Study Design

Multi-modal integration across four data layers:

| Layer | n | Purpose |
|-------|---|---------|
| snRNA-seq (10x Chromium) | 18 liver biopsies (lean n=4, obese n=5, MASL n=4, MASH n=5); 176,150 high-quality nuclei | Discovery |
| Spatial transcriptomics GeoMx (WTA) | 8 MASH F3 biopsies | Validation |
| Spatial proteomics CosMx (1000-plex) | 10 MASLD + 2 end-stage explants | Validation |
| Spatial proteomics COMET Hyperplex | 6 MASLD biopsies | Validation |
| Spatial proteomics MILAN | 8 MASLD biopsies | Validation |
| Bulk RNA-seq | 206 MASLD patients + 10 controls (European NAFLD Registry, GSE135251) | Clinical translation |
| Serum proteomics (SomaScan v4.0) | 247 biopsy-proven MASLD patients | Biomarker validation |
| PheWAS (OLINK Explore 3072) | 53,030 UKB participants | Population-level association |

**Cohort**: Biopsies from bariatric surgery / cholecystectomy (KU Leuven + UZ Antwerp). Fibrosis staged F0–F4 (NASH-CRN). At-risk MASH = NAS ≥ 4 + F ≥ 2.

---

## Main Findings

### 1. snRNA-seq: 8 Hepatic Myeloid Subtypes

176,150 nuclei → 8 major cell types → 12,834 myeloid cells → **8 transcriptionally distinct myeloid subpopulations**:

| Subtype | Markers | Key features |
|---------|---------|-------------|
| **KC (Kupffer cells)** | *MARCO*, *CD5L* | Resident; progressively depleted in MASH |
| **MetMac (Metabolically active macrophages)** | *GPNMB*, *HS3ST2*, *LPL*, *FABP5* | Central finding; expanded in MASH |
| **TransMac (Transitional macrophages)** | *STAB1*, *MGAT4A* | Transcriptional similarity to monocytes + mature macrophages |
| **preMac** | *PCNX2*, *ADAM28*, *RUNX2* | Precursors |
| **Monocyte** | *VCAN*, *FCN1* | Expanded in MASL/MASH |
| **cDC1** | *CLNK*, *CLEC9A* | Conventional dendritic cell |
| **cDC2** | *CD1C* | Conventional dendritic cell |
| **migDC** | *LAMP3*, *CCR7* | Migratory dendritic cell |

**Key population dynamics (MASL → MASH)**:
- **KC proportion decreases** significantly (ANOVA p < 0.001)
- **MetMac and monocyte proportions increase** significantly (p < 0.001)
- migDC and TransMac increase in obesity; monocytes increase in both obesity and MASH

---

### 2. MetMacs (GPNMB+ Macrophages): The Central Finding

**GPNMB+ cells increase from 5.7% (non-MASH) to 12% (MASH)** of total myeloid population.

Within MetMacs, two transcriptional clusters were defined:
- **LPL+ MetMacs** (4.1% of myeloid cells): High *LPL* + *FABP5* + *SPP1*; enriched for LAM-like features; lipid-clearing phenotype
- **HS3ST2+ MetMacs** (2.8% of myeloid cells): High *HS3ST2* (heparan sulfate metabolism) + *HLA-DRA*; pro-inflammatory antigen-presenting phenotype

**Key distinction from prior LAM/SAM classification**:
- LAM signature (*TREM2, LIPA, LPL, CTSB, FABP4, FABP5, LGALS1, LGALS3, CD9, CD36*) detected only in a **fraction** of MetMacs
- GPNMB+ MetMacs are broader, more heterogeneous than Trem2-LAMs
- *TREM2*+ *GPNMB*+ cells = only 0.4% of myeloid cells in MASH
- *GPNMB*+ *LPL*+ *SPP1*+ = 1% of myeloid cells

**GPNMB functional roles**:
- Phagocytosis marker (lysosomal stress → GPNMB induction)
- Cleaved fragment exerts immunomodulatory effects on myeloid/lymphoid cells
- Expression increases under lipotoxic conditions
- Previously identified as MetMac marker in MASLD adipose tissue

---

### 3. Spatial Transcriptomics (GeoMx, at-risk MASH F3): Regional Macrophage Heterogeneity

**Three histological regions segmented** via CD68 (macrophage), CD45 (immune), panCK (epithelial) fluorescence:
- Steatohepatitis (SH) regions
- Parenchyma with low steatosis (LS) regions
- Portal tract (PT) regions

**Steatohepatitis (SH) vs. Low Steatosis (LS) — top DEGs in CD68 segments**:
- Upregulated in SH: *GPNMB*, *LYZ*, *FABP5*, *HLA-DRA*, *LPL*, *HLA-DRB1*, *MSR1*, *CD36* (FDR < 0.05)
- Downregulated in SH: *TIMD4*, *MARCO*, *CCR2*, *CD14*, *TREM2*

**Portal Tract (PT) vs. Steatohepatitis**: 
- PT upregulated: *CCL19*, *TCF7*, *CLEC10A*, *LYRM2* (lymphocyte-regulatory, immuno-monocyte phenotype)
- SH upregulated: *MSR1*, *CD68*, *SPP1*, *CD36*

**Pathway analysis** (SH macrophages vs. KC from LS): enrichment in phagosome, antigen presentation, CD4/CD64 antigen, ITGB3/integrin-beta3.

**Hepatocytes in SH regions** show elevated *Ubiquitin D*, *SQSTM1/p62* (ballooning markers) and cell remodeling genes (*ANXA2*, *MSN*, *KRT8*). Portal tract immune cells elevate *CCL21* and *TRAF1*.

---

### 4. GPNMB+ Macrophages: Spatial Protein Validation (CosMx + COMET + MILAN)

**IHC (n = 16 MASLD + 4 controls)**: GPNMB+ cells increase progressively from control → MASL → MASH F0/1 → MASH F2–4 (ANOVA linear trend p = 0.002).

**CosMx (end-stage liver)**: 30,000 myeloid cells annotated. Of myeloid cells:
- GPNMB+ = 18.1%
- GPNMB+ MARCO+ LYZ+ = 0.8% (KC-like, reduced in MASH)
- GPNMB+ HLA-DR+ (antigen-presenting) = largest GPNMB subpopulation
- GPNMB+ LPL+ FABP5+ SPP1+ = 1%

**MILAN (n = 8 biopsies)**: Increase in GPNMB+ HLA-DR+ cells MASL → MASH; CD163+ cells trend downward; LYZ+ and SPP1+ macrophages unchanged.

**Spatial proximity**: 9% of GPNMB+ macrophages within 10 µm of cholangiocytes → portal infiltration.

---

### 5. IL32-Hepatocyte Axis Drives GPNMB+ Macrophage Phenotype

**Key discovery**: *IL32* is differentially expressed in **hepatocytes** from MASH samples (identified in both snRNA-seq and GeoMx).

**Spatial validation (COMET)**:
- *IL32* mRNA + GPNMB protein co-localize in **end-stage MASH**, with spatial association increasing from MASL → MASH F1 → MASH F3
- *IL32* expression elevated in AKR1B10+ ballooned hepatocytes (ER stress marker)
- Progressive increase in *IL32*+ hepatocytes with disease advancement, especially pericentral zones

**In vitro (PCLS — precision-cut liver slices, ex vivo)**:
- Lipid loading (oleic + palmitic acid) → GPNMB+ aggregates at congested hepatocytes, enlarged lipid-containing cells, ceroid macrophage formation

**THP-1 macrophage conditioning experiments** (IL4, IL32, retinol, LPS + lipid loading):
- IL32 + lipid loading → **enhanced phagocytic uptake** of free lipid particles (endo-lysosomal degradation increased dose-dependently; p < 0.001)
- IL32 + PA → **increased *IL1B* expression** (pro-inflammatory) vs. reduced *IL10* (anti-inflammatory)
- IL32 + PA → **increased GPNMB and MSR1 vesicle formation**
- *TNF* and *GDF1* upregulated by IL32 + PA combination

**Conclusion**: Hepatocyte-derived IL32 (induced by lipotoxic stress) enhances GPNMB+ macrophage phagocytic capacity and shifts pro-inflammatory balance, linking metabolic stress to immune reprogramming.

---

### 6. Clinical Translation: Markers in Bulk RNA-seq and Serum

**Bulk RNA-seq cohort** (n = 206 MASLD patients, European NAFLD Registry, GSE135251):

Genes tracking disease progression (Normal → MASL → MASH F0/1 → F2 → F3 → F4):
- **Upregulated with progression**: *GPNMB*, *LYZ*, *LPL*, *HLA-DRA*, *HLA-DRB1*, *FABP5*, *CCL19*, *DP2*, *CD9*
- **Downregulated with progression**: *CD163*, *ACADSB*, *MSR1* (early-to-mid), *TIMD4*, *MARCO*, *CD5L*, *PIPOX*, *GPX1*

**CIBERSORTx deconvolution** from snRNA-seq signatures applied to bulk data:
- MetMac proportion increases significantly with fibrosis (F0 vs F2–4, p < 0.01) and NAS (NAS 0 vs 4–8, p < 0.01)
- KC proportion decreases significantly
- Monocyte proportion increases

**NAS ≥ 4 prediction (logistic regression)**:
- Combined model: age + *ACP5* + *CD163* + *LPL* + *MSR1* 
- *LPL* was the strongest independent predictor (AUROC shown in Extended Data Fig. 7c)

**Serum proteomics (SomaScan, n = 247)**:
- Soluble **GPNMB** elevated in: at-risk MASH (p < 0.05), advanced fibrosis F3–4 (p < 0.05), high disease activity NAS ≥ 4 (p < 0.05)
- Soluble **MSR1** elevated in at-risk MASH and high disease activity
- **CD163** elevated in high disease activity (p < 0.05) but did not significantly differ for fibrosis stage alone
- **LPL** elevated in high disease activity

**PheWAS (UKB, n = 53,030, OLINK)**:
- Soluble GPNMB: strongest associations with **type 1 diabetes** (-log10p ~80), diabetes mellitus, type 2 diabetes, hypoglycemia, chronic kidney disease

---

### 7. LPL Functional Validation

**LPL (lipoprotein lipase)** — triglyceride hydrolase essential for lipoprotein processing — marked MetMacs enriched in steatohepatitis.

**CRISPR-Cas9 LPL knockout** in THP-1 cells + lipid loading (2000 µM oleic acid, 24h):
- LPL KO → **foamy macrophage phenotype** (Oil Red O staining), demonstrating LPL is required to prevent lipid accumulation
- Consistent with mouse atherosclerosis data showing LPL + MSR1 independently predict high disease activity

---

## Key Gene Markers and Their Roles

| Gene | Subtype | Direction in MASH | Function |
|------|---------|------------------|----------|
| *GPNMB* | MetMac | ↑↑ | Phagocytosis, immunomodulation, lysosomal stress sensor |
| *LPL* | MetMac (LPL+ subset) | ↑ (SH > LS) | Triglyceride hydrolysis, prevents foamy macrophage formation |
| *FABP5* | MetMac | ↑ | Fatty acid binding/transport |
| *HS3ST2* | MetMac (HS3ST2+ subset) | ↑ | Heparan sulfate biosynthesis; pro-inflammatory in inflamed niches |
| *SPP1* | MetMac/monocyte | ↑ | Osteopontin; lipid droplet uptake, lipogranuloma |
| *MSR1* | MetMac | ↑ (SH) | Macrophage scavenger receptor; mediates lipid-induced inflammation |
| *HLA-DRA/DRB1* | MetMac, SH regions | ↑ | Antigen presentation |
| *LYZ* | GPNMB+ cells | ↑ | Lysozyme, antimicrobial |
| *CD163* | KC→MetMac transition | ↓ in KC, ↑ in some MetMacs | Scavenger receptor, hemoglobin clearance |
| *MARCO* | KC | ↓↓ | Scavenger receptor; KC identity/self-renewal |
| *TIMD4* | KC | ↓↓ | KC residence marker |
| *CD5L* | KC | ↓↓ | Lipid metabolism, anti-inflammatory |
| *TREM2* | LAM subset only | Low fraction | Fibrosis resolution (distinct from MetMacs) |
| *ACP5* | MetMac | ↑ | Acid phosphatase; independent predictor of NAS ≥ 4 |
| *CCL19* | Portal macrophages | ↑ (portal tract) | Lymphocyte recruitment |
| *CLEC10A* | Portal macrophages | ↑ (portal tract) | C-type lectin, immuno-regulatory |
| *IL32* | **Hepatocytes** (not macrophages) | ↑ with MASH | Multifunctional cytokine; drives GPNMB+ phenotype paracrine |
| *AKR1B10* | Hepatocytes | ↑ with ballooning | Retinol metabolism enzyme; marks ballooned hepatocytes |

---

## Key Numbers

| Metric | Value |
|--------|-------|
| Total nuclei profiled (snRNA-seq) | 176,150 |
| Myeloid cells subclustered | 12,834 |
| GPNMB+ cells (non-MASH) | 5.7% of myeloid |
| GPNMB+ cells (MASH) | 12% of myeloid |
| GeoMx segments analyzed | 77 (from 80; 3 failed QC) |
| DEGs SH vs LS (CD68 segments) | 207 (FDR<0.05) |
| DEGs PT vs SH (CD68 segments) | 348 (FDR<0.05) |
| CosMx myeloid cells (end-stage) | 30,000 |
| GPNMB+ fraction of myeloid (end-stage CosMx) | 18.1% |
| Bulk RNA-seq clinical cohort | n=206 MASLD + 10 controls |
| Serum proteomics cohort | n=247 |
| PheWAS participants (UKB) | n=53,030 |

---

## Conceptual Model (Extended Data Fig. 9)

**Portal tract in MASH**:
- Loss of MARCO+ CD5L+ TIMD4+ KCs
- Replacement by diverse GPNMB+ subtypes: Kupffer cell-like, HS3ST2+ MetMacs, LPL+ MetMacs (with LAM subset), monocyte-like
- Elevated portal inflammation and phenotypic heterogeneity

**Parenchyma in MASH**:
- Stressed hepatocytes (ballooned, ER stress, AKR1B10+) secrete **IL32** ↑
- IL32 drives GPNMB+ macrophages toward enhanced phagocytosis of free lipid droplets and cell debris
- Two outcomes: (1) Metabolic adaptation → lipogranuloma formation (GPNMB+ HLA-DR– LPL+ FABP5+ SPP1+); (2) Increased lipid-induced inflammation (GPNMB+ CD163+ HLA-DR+) → ceroid macrophage formation
- KC depletion: loss of CD163+ MARCO+ TIMD4+ and CD163+ CD5L+ MARCO+ TIMD4+ populations

---

## Relevance to Our MASLD Project

1. **Spatial data now integrated into our atlas** (2026-05-21): GeoMx SH vs PT logFC/padj and CosMx Hep/KC mash logFC/padj columns added to multi-evidence atlas (`spatial_govaere2026_*` columns). See `CLAUDE.md` S6 spatial sources.

2. **Bulk RNA-seq data**: GSE135251 (n=206, this Govaere cohort) is our **canonical bulk RNA-seq dataset** already in our 5-cohort mega-analysis. This paper's own bulk analysis therefore is directly concordant with our DEG results.

3. **Key gene overlaps to check against our DEGs**: *GPNMB*, *LPL*, *HS3ST2*, *SPP1*, *MSR1*, *ACP5*, *CCL19*, *CLEC10A*, *LYZ*, *FABP5* — all should be in our Tier 1 or Tier 2 DEG list. *IL32* is a hepatocyte-expressed cytokine that likely appears in our hepatocyte-intrinsic or spatial analysis.

4. **Macrophage subtyping**: Their MetMac concept (GPNMB+ as the dominant disease-associated macrophage) aligns with our cell-type attribution results. Their TREM2/LAM distinction (only 0.4% of MASH myeloid cells vs. GPNMB+ at 12%) is a key finding that challenges simple LAM = disease macrophage narratives.

5. **Biomarker angle**: Serum GPNMB and LPL as non-invasive disease activity markers — relevant for our drug repurposing and clinical translation sections. GPNMB's PheWAS association with T1D/T2D/metabolic syndrome has systemic implications.

6. **Competition note**: This paper (Nat Genetics, May 2026) directly overlaps with our macrophage/CCC axis and spatial transcriptomics findings. Our atlas is broader (19 bulk cohorts, cross-species, GWAS causal inference) but they have more spatial depth on macrophage subtypes specifically. Key differentiator: our paper focuses on the multi-evidence convergence of targets for therapeutic intervention, not macrophage biology per se.
