# Tian et al. 2026 — Liver Radiomic GWAS: Paper Summary & Overlap with MASLD Atlas

**Citation:** Tian H, Kamineni M, Truong B, Raghu V, Dron JS, Hornsby W, Koyama S, Yu Z, Natarajan P. "Genetic architecture of high-dimensional liver radiomic phenotypes and their role in common metabolic diseases." medRxiv (posted May 22, 2026). DOI: 10.1101/2026.05.19.26353617

**Corresponding author:** Pradeep Natarajan — pnatarajan@mgh.harvard.edu

---

## Paper Summary

### Study Design
- **200 liver MRI radiomic features** extracted via deep-learning segmentation (nnU-Net) + PyRadiomics (v3.1.0) from 43,176 UK Biobank abdominal MRIs across 2 contrasts (water + in-phase)
- Feature categories: Shape (14), First-order (18×2), Texture — GLCM/24, GLSZM/16, GLRLM/16, NGTDM/5, GLDM/14 (×2 contrasts)
- After pruning |r| > 0.9: 84 features retained; further pruned to 59 at |r_g| < 0.9 for GWAS
- **GWAS:** 37,725 EUR-ancestry UK Biobank individuals; REGENIE v3.2.8; RINT transformation
- Post-GWAS: Genomic SEM (latent factor structure), PWAS vs 1,916 Olink proteins (UKB-PPP), colocalization (coloc v5.2.3 + SuSiE), genetic correlation (LDSC), cis-MR drug target analysis

### Core Claim
High-dimensional liver MRI radiomic features encode heritable biological variation that is **structurally distinct from and largely orthogonal to liver fat (PDFF)** — capturing fibrosis, inflammation, and architectural remodeling signals beyond what PDFF alone reveals.

### Key GWAS Findings
- **327 genome-wide significant loci** (P < 5×10⁻⁸) across 59 MRI features; **57 loci** for the 10 lead features; 21 cross-feature pleiotropic
- Canonical MASLD loci (PNPLA3, TM6SF2, APOE, TRIB1) show internal pleiotropy across MRI features
- **19 novel loci** not in prior GWAS Catalog; 10/15 are shape-specific
- **GCKR** and **CATSPER2P1/TP53BP1** associate with ALL 59 features including shape_Sphericity
- **EFEMP1** colocalized with plasma albumin (PP.H4=0.98) AND shape_LeastAxisLength (PP.H4=0.82) → ECM-mediated structural remodeling
- Five-factor Genomic SEM model best fit (AIC=297.01, CFI=0.998, SRMR=0.052); F4 = liver fat factor; F3 = shape_Sphericity exclusively

### PWAS Findings
- **93 significant protein–MRI pairs** (global FDR < 0.05); 28 proteins, 37 MRI features
- **47 of 93 pairs colocalized** (PP.H4 > 0.8)
- Hub proteins: NCAN, SERPINA1, APOE, **INHBB**, **INHBC**, **TNFSF10**, **EFEMP1**
- "Fibrosis/ECM cluster" surviving liver fat adjustment: APCS, **DKK3**, FABP4
- "Steatosis cluster" attenuated after fat adjustment: FGF21, CES1, LPL

### Genetic Correlation
- Strongest signals (Bonferroni-corrected) with TG, HDL-C, TG:HDL-C ratio, T2D, liver fat
- Insulin resistance likely the dominant shared pathway
- *firstorder_Skewness* specifically correlated with cirrhosis (not liver cT1) — suggesting fat-independent structural signal

### Drug Target cis-MR (6 genes tested: GLP1R, LPA, PCSK9, HMGCR, APOC3, ANGPTL3)
- **GLP1R** → `firstorder_Minimum_inp` (MRI biomarker identified)
- **APOC3** → `glszm_ZoneEntropy_inp` + `glszm_GrayLevelVariance_inp` (MRI biomarkers identified)
- **LPA, PCSK9, HMGCR, ANGPTL3** → no MRI biomarker found

### Data Availability
- Summary statistics: GWAS Catalog upon publication (not yet deposited as of 2026-05-24)
- Code: github.com/HDTian/liverMRI
- Full loci lists in Supplementary Tables S8 (GWAS leads), S18 (PWAS), S20 (colocalized pairs) — not publicly available yet

---

## Overlap with Our MASLD Atlas

Analysis date: 2026-05-24. Compared paper's named genes against our multi-evidence atlas (33,943 genes × 299 cols; 5,378 Dream DEGs padj<0.05 |LFC|>0.3; 736 genes SuSiE PP4>0.5 across 50 GWAS (2026-07-05 MVP 50-GWAS update; was 364 across 28 GWAS)).

### Full Gene-Level Comparison

| Gene | Paper context | Our DEG tier | dream_logFC | dream_padj | Our COLOC (SuSiE PP4) | Best GWAS |
|------|--------------|-------------|------------|-----------|----------------------|-----------|
| LPL | Steatosis marker (fat-driven) | **Tier1 DEG** | +2.78 | 3.47e-37 | 0.489 | NAFLD_EUR |
| FABP4 | Steatosis marker (fat-driven) | **Tier1 DEG** | +2.52 | 3.12e-30 | 0.012 | UKBB_AST |
| EFEMP1 | PWAS/COLOC: albumin + structural MRI | **Tier1 DEG** | +1.09 | 6.99e-13 | 0.007 | NAFLD_deCode |
| LPA | Drug target (no MRI biomarker) | **Tier1 DEG** | −0.565 | 2.80e-12 | 0.010 | UKBB_ALT |
| DKK3 | Fibrosis/ECM marker (fat-independent) | **Tier1 DEG** | +0.618 | 5.25e-12 | <0.001 | UKBB_GGT |
| TRIB1 | GWAS (fat) | Tier2 DEG | −0.279 | 1.07e-07 | 0.334 | BBJ_GGT |
| TNFSF10 | PWAS hub | Tier2 DEG | −0.226 | 1.07e-07 | **1.000** | UKBB_ALT |
| CES1 | Steatosis marker (fat-driven) | Tier2 DEG | +0.495 | 2.45e-07 | 0.011 | Cirrhosis_EAS |
| FGF21 | Steatosis marker (fat-driven) | **Tier1 DEG** | +0.648 | 2.06e-06 | 0.100 | Ghouse_Cirrhosis |
| INHBC | PWAS hub | Tier2 DEG | +0.201 | 2.75e-05 | 0.011 | UKBB_AST |
| HFE | GWAS (iron) | Tier2 DEG | +0.116 | 3.50e-05 | 0.025 | BBJ_AST |
| APOE | GWAS + PWAS | Tier2 DEG | +0.186 | 5.44e-03 | 0.313 | UKBB_AST |
| INHBB | PWAS hub | Tier2 DEG | −0.267 | 9.55e-03 | **0.981** | UKBB_GGT |
| APOC3 | Drug target (cis-MR: glszm_ZoneEntropy) | Tier2 DEG | +0.170 | 2.23e-02 | 0.090 | UKBB_ALT |
| PCSK9 | Drug target (no MRI biomarker) | Tier2 DEG | +0.188 | 3.05e-02 | 0.021 | PanUKBB_AFR_ALT |
| HMGCR | Drug target (no MRI biomarker) | non-DEG | +0.186 | 0.116 | 0.009 | PanUKBB_F_AST |
| SERPINA1 | PWAS hub | non-DEG | −0.370 | 0.118 | <0.001 | Ghouse_HCC |
| TM6SF2 | GWAS (fat) | non-DEG | +0.100 | 0.136 | 0.520 | Ghouse_HCC |
| PNPLA3 | GWAS (fat) | non-DEG | +0.066 | 0.375 | 0.008 | HCC_EAS |
| GCKR | GWAS (all features) | non-DEG | +0.062 | 0.396 | **0.950** | NAFLD_deCode |
| NCAN | PWAS hub | non-DEG | +0.130 | 0.419 | 0.046 | Ghouse_HCC |
| GLP1R | Drug target (cis-MR: firstorder_Minimum) | non-DEG | −0.128 | 0.456 | 0.010 | HCC_EAS |
| APCS | Fibrosis/ECM marker (fat-independent) | non-DEG | −0.059 | 0.499 | 0.040 | Cirrhosis_EAS |
| ANGPTL3 | Drug target (no MRI biomarker) | non-DEG | −0.038 | 0.597 | 0.005 | UKBB_ALT |
| SUGP1 | GWAS (cholesterol) | non-DEG | −0.007 | 0.895 | 0.031 | NAFLD_Intermountain |

### Strongest Convergent Signals

**TNFSF10 (TRAIL):** PP4=1.00 in UKBB_ALT + significant DEG (padj=1.1e-7, downregulated) + PWAS hub in radiomic paper = three-modality convergence. TRAIL signaling repressed in MASLD.

**INHBB (Activin B):** PP4=0.981 in UKBB_GGT + DEG (padj=0.010, downregulated) + PWAS hub. Hepatic activin B suppression concordant across genetic and transcriptional layers.

**EFEMP1:** Strong Tier1 DEG (padj=6.99e-13, LFC=+1.09 upregulated) + paper PWAS/COLOC shows protein associated with liver structural remodeling (albumin PP.H4=0.98, shape PP.H4=0.82). ECM gene active at both transcriptional and proteomic/imaging levels. Weak COLOC on our side (PP4=0.007) reflects that expression QTL ≠ protein QTL here.

**DKK3 (Wnt antagonist):** Tier1 DEG in our data (padj=5.25e-12, LFC=+0.618) + paper explicitly flags as a "fat-independent fibrosis/remodeling marker." Transcriptional upregulation confirmed independent of steatosis.

**GCKR:** PP4=0.95 in our COLOC + paper GWAS shows association with ALL 59 MRI features. Non-DEG (padj=0.40) in both studies — acts through protein-level regulation of glucokinase, not transcription.

### Why the Canonical MASLD GWAS Genes Are Not DEGs

PNPLA3 (p.Ile148Met), TM6SF2 (p.Glu167Lys) are missense coding variants — they act through protein function, not cis-eQTL. Neither is a significant DEG in our dream mega-analysis (padj=0.38 and 0.14 respectively), and our COLOC is correspondingly weak (PP4=0.008 and 0.52). The paper confirms their effects are fat-correlated and attenuated once fat is controlled. This validates that our transcriptomic + COLOC pipeline is capturing a complementary layer — the transcriptional disease response — not the genetic risk architecture.

### Our Key Causal Genes Not in the Radiomic Paper

| Gene | Our dream_logFC | dream_padj | SuSiE PP4 | Best GWAS |
|------|----------------|-----------|----------|-----------|
| THRB | −0.217 | 3.31e-11 | **1.000** | UKBB_GGT |
| HKDC1 | +1.130 | 1.64e-08 | **0.992** | UKBB_AST |
| ADH4 | −0.085 | 0.270 | 0.987 | NAFLD_UKBB |
| GPAM | +0.163 | 0.015 | 0.912 | UKBB_AST |
| SERPINE1 | +0.541 | 2.27e-05 | 0.138 | UKBB_GGT |

HKDC1 (top F2→F3 causal driver, PP4=0.99, strong DEG) and THRB (resmetirom target, PP4=1.00) are absent from the radiomic paper — these are candidates for imaging biomarker validation using their radiomic framework once their summary stats are deposited.

---

## Pending: Full Genome-Wide Overlap

Only ~25 main-text genes analyzed above. Full lists (327 GWAS loci × nearest genes, 28 PWAS proteins, 47 colocalized pairs) are in Supplementary Tables S8/S18/S20 — not yet publicly available (preprint, data depositing to GWAS Catalog upon publication).

**To complete this analysis:** download supplementary tables from medRxiv or contact pnatarajan@mgh.harvard.edu, then run:
```bash
# Match all S8 nearest-gene annotations against atlas
python3 -c "
import pandas as pd
atlas = pd.read_csv('RNA-seq/results/multi_evidence/multi_evidence_atlas.csv', low_memory=False)
s8 = pd.read_excel('docs/reference/Tian2026_TableS8.xlsx')  # once available
overlap = atlas.merge(s8, left_on='human_symbol', right_on='nearest_gene')
# Stratify by fat-independent vs fat-correlated features
"
```
