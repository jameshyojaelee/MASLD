> **Historical competitor/source dossier only (2026-08-10).** Source-paper
> facts may remain useful, but project comparisons, figures, targets, journals,
> and next-analysis suggestions are non-authoritative. Reverify source facts;
> use [`../../PAPER.md`](../../PAPER.md),
> [`../../STATUS.md`](../../STATUS.md), and
> [`../../RESULTS.md`](../../RESULTS.md) for current claims.

# Govaere et al. 2026 — Nature Genetics — Key Takeaways

**Full title:** Integrated multi-omics identifies distinct macrophage alterations during progression of metabolic dysfunction-associated steatohepatitis  
**Journal:** Nature Genetics (accepted April 14, 2026; published online May 18, 2026)  
**DOI:** https://doi.org/10.1038/s41588-026-02600-3  
**Corresponding authors:** Quentin M. Anstee, Schalk van der Merwe, Olivier Govaere  
**Data:** snRNA-seq = EGA EGAS50000000768; bulk RNA-seq = GSE135251; CosMx + GeoMx = GSE312698  
**Code:** github.com/govaerelab/OG_MB-MASLD-LIVER-2025 + github.com/MarkusBoesch93/OG_MB-MASLD-LIVER-2025

---

## Study Design & Cohorts

| Modality | n | Disease stages |
|---|---|---|
| snRNA-seq (10x Chromium, liver nuclei) | 18 biopsies → 176,150 nuclei | lean=4, obese=5, MASL=4, MASH=5 |
| Spatial transcriptomics GeoMx (WTA) | 8 MASH F3 biopsies | at-risk MASH (fibrosis ≥2, NAS ≥4) |
| CosMx spatial molecular imaging (1000-plex) | 8 MASLD + 2 normal explants | MASL + MASH |
| COMET Hyperplex proteomics | 6 sequential samples | MASL + MASH |
| MILAN multiplex proteomics | 8 biopsies | MASL–MASH |
| Bulk RNA-seq (European NAFLD Registry) | 206 MASLD + 10 controls | Normal–MASL–MASH F0-F4 (GSE135251) |
| SomaScan serum proteomics | 247 biopsy-proven MASLD patients | at-risk MASH spectrum |
| PheWAS (OLINK proteomics, UKB) | 53,030 participants | — |

**At-risk MASH** = MASH + fibrosis stage ≥2 + NAS ≥4 (at least 1 point per NAS component).

---

## Central Finding

Disease progression from MASL to MASH is marked by:
1. **Progressive depletion of resident Kupffer cells (KCs)**
2. **Emergence of phenotypically distinct GPNMB+ macrophage subsets** with metabolic-adaptive and pro-inflammatory phenotypes
3. These changes are tightly regulated by **spatial context** (portal tract vs. steatohepatitis parenchyma) and **disease stage**
4. **IL32-producing stressed hepatocytes** drive the phagocytic and pro-inflammatory phenotype of GPNMB+ macrophages

---

## Macrophage Subpopulations (snRNA-seq)

8 cell types identified from 176,150 nuclei. Myeloid compartment yielded 10 transcriptionally distinct subclusters:

| Cluster | Key markers | Disease behavior |
|---|---|---|
| KC (Kupffer cell) | MARCO, CD5L, TIMD4, CD163 | Progressively depleted from MASL→MASH |
| MetMac (metabolically active Mac) | GPNMB, HS3ST2, LPL, FABP5 | Expands in MASH; two sub-clusters |
| TransMac (transitional Mac) | STAB1, MGAT4A | Increases modestly |
| preMac | PCNX2, ADAM28, RUNX2 | Transcriptional similarity to monocytes and mature macs |
| Monocyte | VCAN, FCN1 | Increases in MASH |
| cDC1 | CLNK, CLEC9A | — |
| cDC2 | CD1C | — |
| migDC | LAMP3, CCR7 | — |

**Key quantitative shifts (MiloR differential abundance):**
- KC proportion: significantly decreased in MASH
- MetMac and monocyte proportions: significantly increased in MASH
- MARCO+ CD5L+ TIMD4+ KCs: gradually declining with fragmented sub-populations emerging

---

## GPNMB+ Macrophages: The Key Cell Type

### Abundance
- GPNMB+ cells rise from **5.7% of total myeloid population in non-MASH** to **12% in MASH**
- Within GPNMB+ cells:
  - GPNMB+LPL+: 2.8% of myeloid cells in MASH
  - GPNMB+HS3ST2+: 4.1% of myeloid cells in MASH
  - GPNMB+MARCO-: decreased from 1.3% (non-MASH) → 0.8% (MASH) — KC-like cells reducing
  - TREM2+GPNMB+: only 0.4% in MASH (distinct from described LAMs/SAMs)
  - GPNMB+LPL+SPP1+: 1% (MASH)

### Transcriptional profile in MASH vs non-MASH
- **Upregulated:** LPL, FABP5, SPP1 (and concomitant downregulation of CD5L, TIMD4)
- LPL = triglyceride hydrolase essential for lipoprotein processing → marks macrophages enriched in steatohepatitis; required to prevent foamy macrophage formation (CRISPR-Cas9 LPL KO → foamy macrophage phenotype)
- HS3ST2 = heparan sulfate biosynthesis regulator; associated with pro-inflammatory phenotype in atherosclerosis

### Spatial distribution (CosMx 1000-plex)
- GPNMB+ cells cluster in **portal tracts** and **hepatic parenchyma**
- In lipogranulomas: high LPL and FABP5 expression
- In portal tract: high HLA-DRA expression
- 9% of GPNMB+ macrophages located within 10 µm of cholangiocytes (portal infiltration)

### LAM signature
- Trem2-associated LAM signature (TREM2, LIPA, LPL, CTSB, FABP4, FABP5, LGALS1, LGALS3, CD9, CD36) detected in a **subset** of MetMac cluster, NOT globally enriched in GPNMB+ cells
- This distinguishes MetMacs from previously described LAMs/SAMs

---

## Spatial Transcriptomics (GeoMx, n=8 MASH F3 biopsies)

Three histological regions profiled per biopsy: steatohepatitis (SH), parenchyma with low steatosis (LS), portal tract (PT).

### Steatohepatitis (SH) vs. Low Steatosis (LS) macrophages
207 DEGs (FDR < 0.05). SH-macrophages upregulate:
- GPNMB, LYZ, FABP5, LPL, HLA-DRA, CD63, MSR1, CD36
- Enriched for phagosome, antigen-presenting, and phagocytic pathway genes
- Show features of both monocytes (high LSP1 and LYZ) and mature macrophages (high MSR1)

### Portal Tract (PT) macrophages
348 DEGs (PT vs SH). PT-macrophages express:
- CCL19, CD48, IL7R, CCL5 — lymphocyte-regulatory genes
- Typical immuno/monocyte phenotype lacking scavenger receptors (CD163, CD36, MSR1)

### Hepatocytes in steatohepatitis regions
- Increased SQSTM1/p62, Ubiquitin D (ballooning markers)
- ANXA2, MSN, KRT8 (cell remodeling)
- Elevated IL32 expression (ER stress marker) — validated spatially and across disease stages

---

## IL32–Macrophage Crosstalk (Mechanistic Core)

**Discovery:** IL32 identified as differentially expressed in MASH hepatocytes in both snRNA-seq and GeoMx datasets.

**Spatial validation:**
- IL32+ hepatocytes co-localize with GPNMB+ macrophages in steatohepatitis regions
- IL32 expression increases progressively from MASL→MASH F1→MASH F3 (CosMx and COMET)
- Elevated in periportal zones and areas with ballooned hepatocytes (AKR1B10-marked)

**In vitro (THP-1 + lipids + IL32):**
- Enhanced phagocytic uptake of free lipid particles (dose-dependent, IL32 50/100 ng/mL)
- Increased GPNMB aggregate formation in cytoplasm; enlarged MSR1 and GPNMB vesicles
- Enhanced endo-lysosomal degradation activity
- Ceroid macrophage formation (lipogranuloma)
- IL32 upregulates **IL1B** while reducing anti-inflammatory **IL10**
- TNF and GDF1 expression also increased with IL32 + palmitic acid

**LPL CRISPR-Cas9 KO:** → foamy macrophage phenotype after lipid loading (blocked clearance)

**Conclusion:** Hepatocyte-derived IL32 links lipotoxic stress to immune reprogramming by enhancing GPNMB+ macrophage phagocytic capacity and promoting a pro-inflammatory microenvironment.

---

## Clinical Translation

### Bulk RNA-seq (n=206 MASLD, European NAFLD Registry)
28-gene heatmap (Fig. 6a) tracks expression from MASL → MASH F0/1 → F2 → F3 → F4:
- Progressively upregulated: GPNMB, LYZ, HS3ST2, LPL, HLA-DRA, LSP1, FABP5, CCL19, DP2
- Progressively downregulated: CD163, ACADSB, MSR1, MTX, GPX1, TIMD4, MARCO, CD5L, CDSL

**Logistic regression model (NAS ≥ 4 classification):**
- Independent predictors: age, ACP5, CD163, LPL, MSR1
- LPL expression contributed most strongly

### Serum proteomics (SomaScan, n=247)
Significant differences (p < 0.05) between MASH vs non-MASH:
- **Soluble GPNMB**: elevated in advanced MASH (strongest signal)
- **Soluble CD163**: elevated in some MASH stages
- **MSR1**: elevated in at-risk/advanced MASH

### PheWAS (OLINK, n=53,030 UKB participants)
- **Soluble GPNMB** → strongly associated with diabetes mellitus, hypoglycemia, Type 2 diabetes, chronic kidney disease

---

## Spatial Proteomics Validation (COMET + MILAN)

**COMET (n=6 sequential patient samples):**
- Protein-level confirmation of snRNA-seq findings
- In MASH: GPNMB+ cells comprise **18.1% of myeloid population** (CD68+, CD163+, or CD14+)
- GPNMB+MARCO+ LYZ+ cells: 0.5% (KC-like, declining)
- GPNMB+HLA-DR+ macrophages significantly increase from MASL→MASH

**MILAN (n=8 biopsies, MASL + MASH):**
- 30,000 myeloid cells segmented across 8 biopsies
- GPNMB+ HLA-DR+ cells: significantly increased in MASH vs MASL
- CD163+ cells: trend toward reduction in MASH
- LYZ+ and SPP1+ macrophage fractions: unchanged between MASL and MASH

**IHC (n=20 biopsies, 16 MASLD + 4 controls):**
- Progressive rise in GPNMB immunopositivity from normal → MASL → MASH F0-1 → MASH F2-4
- ANOVA linear trend p = 0.002

---

## Key Marker Summary for Cas13 Library Design

| Gene | Direction in MASH | Evidence | Notes |
|---|---|---|---|
| GPNMB | UP | snRNA-seq, GeoMx, CosMx, IHC, serum proteomics, PheWAS | Master macrophage marker; 5.7%→12% of myeloid |
| LPL | UP | snRNA-seq, GeoMx, serum, logistic regression | LPL KO → foamy mac; blocks lipid clearance |
| FABP5 | UP | snRNA-seq, GeoMx, CosMx | Fatty acid binding in MetMacs |
| HS3ST2 | UP | snRNA-seq, GeoMx | Heparan sulfate biosynthesis |
| SPP1 | UP | snRNA-seq, GeoMx | Shared KC↓ and Mac↑ context |
| LYZ | UP | snRNA-seq, GeoMx, serum | Shared across myeloid subtypes |
| HLA-DRA | UP (SH) | snRNA-seq, GeoMx, CosMx | Antigen presentation in GPNMB+ macs |
| MSR1 | UP | snRNA-seq, bulk RNA-seq, serum | Macrophage scavenger receptor; lipid droplet uptake |
| ACP5 | UP | GeoMx, logistic regression | Portal enriched |
| CCL19 | UP (PT) | GeoMx, bulk RNA-seq | Portal macrophage marker |
| IL32 | UP (hepatocytes) | snRNA-seq, GeoMx, CosMx, COMET | Hepatocyte→macrophage paracrine signal |
| MARCO | DOWN | snRNA-seq, CosMx, bulk RNA-seq | KC depletion marker |
| CD5L | DOWN | snRNA-seq, bulk RNA-seq | KC depletion marker |
| TIMD4 | DOWN | snRNA-seq, bulk RNA-seq | KC depletion marker |
| CD163 | DOWN (tissue) / UP (serum) | snRNA-seq, bulk RNA-seq, serum, PheWAS | Complex; tissue CD163+ KCs decline, soluble shed |
| TREM2 | LOW | snRNA-seq, CosMx | Only 0.4% TREM2+GPNMB+ in MASH |

---

## Relevance to Our MASLD Atlas

### Overlap with our pipeline
- **GSE135251** (bulk RNA-seq n=216) is used in both their study and ours — their Fig. 6 uses this cohort
- **GeoMx data (n=8)** now deposited at GSE312698 and Zenodo 10.5281/zenodo.18153344 — directly relevant to our Govaere2026 GeoMx integration already in our spatial evidence layer (`spatial_govaere2026_geomx_*` columns in atlas, integrated 2026-05-21)
- **CosMx data (n=8-10)** at GSE312698 — we have integrated `spatial_govaere2026_cosmx_*` columns for hep and KC mash logFC/padj

### What they have that we don't
- Single-nucleus resolution macrophage subclustering (n=176,150 nuclei) — our scRNA data does not focus on myeloid subclusters at this resolution
- GPNMB+ MetMac subcluster functional validation (PCLS, CRISPR)
- IL32 paracrine hepatocyte→macrophage mechanism (in vitro)
- Protein-level spatial (COMET, MILAN) confirmation
- PheWAS for serum GPNMB and CD163

### Differentiation of our work from Govaere 2026
- Our atlas covers 9 human bulk cohorts (1,259 samples) + 5 mouse diet models + 7 modalities → broader transcriptomic coverage
- We use genetic causal inference (TWAS + SuSiE-COLOC across 28 GWAS) — they do not
- We have drug repurposing (LINCS/CGP), NMF subtyping (k=6 programs), cross-species concordance, sex stratification
- Their macrophage focus is orthogonal to our multi-evidence convergence approach — cite as complementary single-cell resolution

### Genes to check in our convergence atlas
GPNMB, LPL, FABP5, HS3ST2, SPP1, MSR1, LYZ, IL32, ACP5, CCL19, MARCO, TIMD4 — all relevant cross-reference targets between their macrophage findings and our DEG/COLOC evidence

---

## Limitations Acknowledged by Authors
- snRNA-seq subject to technical dropouts (lower mRNA abundance in nuclei vs whole cells)
- High GPNMB+ macrophage heterogeneity makes precise cluster boundary definition difficult
- GeoMx enrichment-based fluorescence approach may include neighboring cells upon cell segmentation
- COMET focal nuclear background staining (DNA binding artifact)
- Control samples were from healthy-by-BMI patients (bariatric/cholecystectomy), not truly healthy individuals
- Analyses performed on relatively small cohorts — larger independent cohort validation warranted

---

*Summary compiled 2026-05-22 from full paper (39 pages). For spatial data columns already integrated into our atlas, see `docs/progress.md` and `RNA-seq/results/multi_evidence/multi_evidence_atlas.csv`.*
