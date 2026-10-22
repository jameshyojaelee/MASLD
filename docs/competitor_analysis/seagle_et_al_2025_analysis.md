> **Historical competitor/source dossier only (2026-08-10).** Source-paper summaries may remain useful, but all “our work,” winner, figure, journal, target, and next-analysis sections predate the standalone Resource scope and may be stale. Reverify source facts, and use [`../PAPER.md`](../PAPER.md), [`../STATUS.md`](../STATUS.md), and [`../RESULTS.md`](../RESULTS.md) for current project claims.

# Seagle et al. 2025 (AJHG) -- Detailed Results & Comparison to MASLD Atlas

**Paper**: Genomics-informed drug-repurposing strategy identifies two therapeutic targets for preventing liver disease associated with metabolic dysfunction  
**Authors**: Seagle, Akerele, DeCorte, Hellwege, Breeyear, Kim, Levin, Khodursky, Bress, Lee, Meiler, Gill, Lee, Heberer, Miller, Reaven, Chang, Lynch, VA MVP, Khankari, Shuey, Edwards, Vujkovic  
**Journal**: The American Journal of Human Genetics 112, 1778--1791, August 7, 2025  
**DOI**: https://doi.org/10.1016/j.ajhg.2025.06.014  
**License**: CC BY-NC-ND 4.0  

---

## 1. Study Design and Pipeline Architecture

### 1.1 Source GWAS

The study builds on Vujkovic et al. 2022 (Nat Genet), a multi-ancestry GWAS in the Million Veteran Program (MVP):

- **Cases**: 90,408 MASLD subjects
- **Controls**: 128,187 individuals
- **Total N**: ~218,595
- **Phenotype definition**: ALT-based proxy (elevated ALT >40 U/L males, >30 U/L females at 2+ timepoints 6 months apart within a 2-year window) + exclusion of other liver diseases, chronic conditions, and alcohol use disorders
- **Discovery output**: 77 independent MASLD-associated variants, 81 putative causal genes
- **External replication**: 17 of 77 variants replicated in an independent histologic cohort (biopsy-defined MASLD) and a radiologic cohort (CT/MRI hepatic triglyceride content)

**Critical limitation acknowledged by authors**: The ALT-based phenotype could capture signals related to enzyme kinetics or systemic perturbations rather than MASLD pathology per se.

### 1.2 Pipeline Overview (Figure 1)

The analytical framework has three sequential stages:

1. **Gene-Set Analysis**: GWAS variants -> putative causal genes via 8 lines of evidence from the source GWAS + new GPGE (S-PrediXcan) + pQTL MR/colocalization -> 212 total genes
2. **Drug Candidate Identification**: 212 genes -> drug-target mapping (OpenTargets, DGIdb, ChEMBL, DrugBank) -> directional concordance filtering -> 57 drug-gene pairs
3. **Drug Candidate Evaluation**: MR proxying drug effects via GPGE + protein structural modeling (RosettaLigand, AlphaFold3, MD) -> 2 prioritized candidates (icosapent ethyl, fingolimod)

---

## 2. Detailed Results

### 2.1 Gene Discovery (212 Putative Causal Genes)

Starting from 81 genes from the source GWAS, additional genes were identified through two approaches:

| Source | Genes Added | Notes |
|--------|-------------|-------|
| Source GWAS (8 lines of evidence) | 81 | Coding variant linkage, eQTL coloc, splicing, histone/chromatin, liver enrichment, PPI network |
| pQTL MR + colocalization | 3 | Using deCODE GWAS of 4,907 circulating proteins in 35,559 Icelanders; cis-acting variants within +/-250 kb; IVW MR with 1000G EUR reference; FDR 5% across 1,649 tests; coloc PP.H4 > 0.70 |
| S-PrediXcan GPGE analysis | 160 (131 net new) | 48 GTEx v7 tissues; Bonferroni threshold p < 2.03 x 10^-7 (217,391 tests); 29 overlapped with GWAS gene set |

**Final tally**: 212 unique putative causal genes, of which **158 (74.5%) were previously unreported** for MASLD or liver enzyme levels in the GWAS Catalog or OpenTargets.

**Gene-tissue pair results**: 412 statistically significant gene-tissue pairs (p < 2.03 x 10^-7) with colocalization PP.H4 > 0.80. Of these, 23 had significant MR results.

**Notable finding**: S1PR2 expression was strongly associated with MASLD in liver tissue specifically (p = 4.63 x 10^-9, PP.H4 = 0.90). This liver-specificity is relevant because it was the only gene with strong signal in liver tissue among the final prioritized targets.

### 2.2 Pathway Analysis (IPA)

Ingenuity Pathway Analysis on all S-PrediXcan-significant genes yielded:

| Pathway/Category | p-value | Context |
|------------------|---------|---------|
| FXR/RXR activation | 6.81 x 10^-12 | Top canonical pathway; consistent with bile acid signaling role in MASLD |
| LXR/RXR activation | 9.93 x 10^-8 | Second canonical pathway; lipid/cholesterol metabolism |
| Gastrointestinal/hepatic system disease | Significant | Expected for MASLD |
| Lipid metabolism | Significant | Core MASLD biology |
| Increased potassium, ALP, ALT, albumin | Significant | Disease-relevant biomarkers |

The prominence of rexinoid (RXR) activity supports the observed PPARG signal and suggests a broader nuclear receptor disease mechanism.

### 2.3 Drug Candidate Identification (57 Drug-Gene Pairs)

From 212 genes, 57 drugs with MASLD repurposing potential were identified:

- **Gene targets**: 13 unique genes (6.67% of TWAS genes) + 3 proteins (7.3% of pQTL genes)
- **Of 57 drugs**: 11 were local anesthetics (mapped to SCN2A)
- **Hepatotoxicity screening**:
  - 33 of 57 listed in LiverTox or FDA DILI database
  - 5 flagged as Class A (well-known cause of clinically apparent liver injury)
  - 3 flagged as Class B (likely rare cause)
  - 10 flagged as Class C (probable cause)
  - 8 classified by FDA as "most concern" for DILI with confirmed causal evidence
- **Cancer drugs excluded** due to severe adverse event profiles

**Directional concordance logic**: Inhibitors/blockers selected for genes where increased expression associates with increased MASLD risk. Agonists/activators selected for the converse. This is the core filtering step.

### 2.4 MR Results (Drug Candidate Evaluation)

The key MR results (Table 1) testing proxied drug effects on MASLD risk:

| Primary Indication | Gene Target | Proxied Drug | Drug Action | No. GPGE Tissues | MR IVW OR (95% CI) | p-value | Egger p |
|--------------------|-------------|-------------|-------------|------------------|--------------------|---------|---------|
| Atrial fibrillation | SCN2A | riluzole, etc. | inhibitor | 1 | 0.86 (0.24--3.14) | 0.82 | N/A |
| Epilepsy | SCN2A | erlosomide, etc. | inhibitor | 2 | 0.86 (0.23--3.26) | 0.83 | N/A |
| Hypertriglyceridemia | FADS3* | icosapent ethyl | inhibitor | 8 | 0.31 (0.22--0.44) | 0.86 | 0.92 |
| **Hypertriglyceridemia** | **FADS1** | **icosapent ethyl** | **agonist** | **23** | **0.28 (0.26--0.31)** | **6.8 x 10^-160** | **0.29** |
| **Hypertriglyceridemia** | **FADS2** | **icosapent ethyl** | **inhibitor** | **26** | **0.25 (0.25--0.28)** | **3.4 x 10^-127** | **0.97** |
| **Multiple sclerosis** | **S1PR2** | **fingolimod, etc.** | **agonist** | **2** | **0.91 (0.86--0.97)** | **4.0 x 10^-4** | **N/A** |
| Type 2 diabetes | PPARG | pioglitazone HCl | agonist | 3 | 1.93 (1.62--2.31) | 1.2 x 10^-14 | 0.90 |

*FADS3 standard error inflated for significant heterogeneity; random-effects IVW used.

**Key quantitative findings**:

- **FADS1 (icosapent ethyl proxy)**: OR = 0.28 per SD increase in expression. This is a large predicted protective effect. 23 tissue-specific GPGEs used as instruments. No evidence of directional pleiotropy (Egger p = 0.29).
- **FADS2 (icosapent ethyl proxy)**: OR = 0.25 per SD decrease in expression. Consistent with FADS1 direction. 26 tissues. No pleiotropy (Egger p = 0.97).
- **S1PR2 (fingolimod proxy)**: OR = 0.91 per SD increase in expression. Modest protective effect. Only 2 tissue-specific GPGEs (limited power for sensitivity analyses). No Egger test possible with 2 instruments.
- **PPARG (pioglitazone proxy)**: OR = 1.93 per SD increase in expression. Strong risk-increasing effect. Seemingly paradoxical given pioglitazone's clinical benefit in MASH (discussed in detail below).

**Sensitivity analyses**: MR Egger conducted for directional pleiotropy assessment. Q-statistic for heterogeneity. FADS3 showed significant heterogeneity (random-effects IVW used). No other significant pleiotropy detected.

### 2.5 Protein Structural Modeling

#### PPARG-Pioglitazone

- Used experimental cocrystal structure (>50 structures available)
- RosettaLigand docking: pioglitazone had superior IDX (binding energetics) compared to nearly all other docked compounds
- Binding mode: 2,4-thiazolidinedione ring engages H323, Y473, H449 via hydrogen bonds; central benzene ring and pyridine tail contact hydrophobic pocket interior
- **Implication**: Pioglitazone is a strong PPARG binder, but the MR suggests increased PPARG expression increases MASLD risk -- interpreted as a contra-indication rather than a therapeutic target

#### FADS1-Icosapent Ethyl

- No experimental structures exist for any FADS protein
- Generated AlphaFold3 model of FADS1 bound to unsaturated phospholipid
- Binding pocket confirmed via DiffDock (blind docking to apo-state model)
- RosettaLigand docking: lipid tail engages phenylalanyl and tyrosyl residues in hydrophobic pocket; ester head group forms hydrogen bond to ASN
- **MD simulations**: Three independent 100-ns unbiased trajectories in ER-like membrane (60% POPC, 30% POPE, 5% POPS, 5% cholesterol). Icosapent ethyl remains bound throughout. Head group interactions stable; tail shows expected flexibility (RMSD 1.5--2.75 A from initial pose)
- **Novel contribution**: First published computational model of FADS1 bound to icosapent ethyl

#### FADS2-Icosapent Ethyl

- Modeled via identical procedure to FADS1
- Near-identical binding pocket (high sequence and structural homology)
- Average IDX: FADS1 = -14.199, FADS2 = -14.058 (nearly equivalent binding)
- MD simulations not performed (justified by binding pocket homology)

#### S1PR2-Fingolimod

- Used experimental structure of S1PR2 bound to S1P (PDB: 7T6B)
- Fingolimod has structural similarity to S1P
- RosettaLigand docking: fingolimod's binding mode highly homologous to S1P
- Both engage hydrophobic core with unsaturated alkyl chains; both form hydrogen bonds via hydrophilic head groups
- **Key difference**: S1P's larger hydrophilic head group allows additional hydrogen bonds to S1PR2 that fingolimod cannot access
- Fingolimod has known efficacy through S1PR1 binding; low affinity for S1PR2
- **Implication**: Structural modeling supports binding but highlights that fingolimod is a weaker S1PR2 ligand than the endogenous ligand

### 2.6 The PPARG Paradox

This is the most nuanced finding in the paper and requires detailed treatment:

**The contradiction**: Increased PPARG expression predicts increased MASLD risk (OR 1.93, p = 1.2 x 10^-14), yet PPARG agonists (pioglitazone) show clinical benefit in MASH resolution.

**Author explanations**:

1. **Cell-type-specific effects**: Hepatocyte-specific PPARG acts as a steatogenic factor (mice with hepatocyte-specific Pparg knockout show reduced steatosis). PPARG overexpression increases steatogenic gene expression. The clinical benefit of TZDs comes from peripheral tissue insulin sensitization, not hepatic PPARG activation.
2. **Diet-dependent effects**: Pioglitazone/rosiglitazone increases liver steatosis in mice on high-fat diets. In overweight mice, pioglitazone decreases hepatic fat accumulation. The direction of effect depends on metabolic context.
3. **Pharmacological vs. intrinsic activation**: Differential activation patterns and cell-specific effects may explain why pharmacological PPARG activation (pioglitazone) can be beneficial while genetic predisposition to higher PPARG expression is detrimental.
4. **GPGE limitation**: The GPGE approach does not capture context-dependent expression changes, which is particularly problematic for a gene like PPARG where cell-type-specific expression matters.

**Clinical context**: EASL/EASD/EASO 2024 guidelines assign pioglitazone a weak recommendation for MASLD and do not recommend it as MASH-targeted therapy, partly due to adverse effects (edema, fluid retention, heart failure risk).

---

## 3. Stated Limitations

The authors acknowledge the following:

1. **ALT-based phenotype**: Signals may reflect enzyme kinetics rather than MASLD pathology
2. **European ancestry expression models**: S-PrediXcan models derived from GTEx (European ancestry), applied to multi-ancestry GWAS -- therapeutic targets may not be relevant across all populations
3. **Gene expression as drug proxy**: Drugs act on proteins, not gene expression; using GPGE to proxy drug effects is an approximation
4. **Missing positive controls**: FGF21 and GLP1R absent from findings due to the analytical framework prioritizing only GWAS-derived effector genes
5. **MR instrument limitations**: Acknowledged standard MR weaknesses (Table S1)
6. **FADS3 heterogeneity**: Tissue-specific effects were heterogeneous, requiring random-effects model
7. **S1PR2 instrument count**: Only 2 tissue-specific GPGEs available, precluding robust sensitivity analyses
8. **Hepatotoxicity not fully resolved**: Directional concordance filtering did not eliminate all hepatotoxic drug candidates (5 Class A, 3 Class B, 10 Class C drugs remained among the 57)
9. **No experimental validation**: Entirely computational; no cell-based or animal model confirmation

---

## 4. Comparison with Our MASLD Atlas

### 4.1 Fundamental Architectural Differences

| Dimension | Seagle et al. 2025 | Our MASLD Atlas |
|-----------|--------------------|--------------------|
| **Primary question** | Which existing drugs can be repurposed for MASLD prevention? | What are the transcriptomic, genetic, and multi-modal signatures of MASLD, and what targets emerge from their integration? |
| **Starting point** | GWAS summary statistics (no diseased tissue data) | Diseased-tissue RNA-seq (1,444 human samples, ~464 mouse) |
| **Unit of analysis** | Drug-gene pair | Gene across 7 evidence dimensions |
| **Source GWAS** | Vujkovic et al. 2022 (MVP; N=218K; ALT-based) | Ghodsian et al. + Chen et al. + UKBB + FinnGen + BBJ |
| **TWAS/GPGE** | S-PrediXcan across 48 GTEx tissues; gene discovery tool | S-PrediXcan (liver tissue); sc-TWAS across 4 cell types; causal inference module |
| **Colocalization** | eQTL + pQTL coloc with PP.H4 > 0.80 | coloc.abf against Broadaway liver eQTLs (N=1,183) and 4 European GWAS; 348 genes with PP.H4 > 0.5 |
| **Drug repurposing approach** | GWAS -> gene -> drug database mapping -> directional concordance -> MR proxy | LINCS L1000 reversal (51,000 perturbations) + CGP pathway reversal + network proximity + DGIdb |
| **Protein modeling** | RosettaLigand, AlphaFold3, DiffDock, 100-ns MD simulations | None |
| **Hepatotoxicity screening** | LiverTox + FDA DILIRank | Not systematically implemented |
| **Cell-type resolution** | None (bulk tissue-level only) | MuSiC/BayesPrism deconvolution; sc-TWAS (4 cell types); cell-type-specific LINCS signatures |
| **Disease progression** | Binary (MASLD yes/no) | NAS score, fibrosis stage, NAFL-vs-NASH, stage-transition programs, molecular subtypes |
| **Cross-species** | None | 10,932-gene concordance atlas across 5 mouse diet models |
| **Spatial transcriptomics** | None | Visium (HRA007511, GSE192741) |
| **Proteomics** | pQTL MR only (circulating proteins) | Govaere tissue proteomics (PXD052937) + SomaScan (GSE276114) |
| **Epigenomics** | None | snATAC-seq (GSE244832); SCENIC+ GRN; ATAC-seq (GSE246213) |
| **Machine learning** | None | RF/SVM/XGBoost target prediction; SHAP importance; source ablation |
| **Multi-evidence integration** | Sequential filtering (GWAS -> TWAS -> drug mapping -> MR -> structural) | Bayesian integration framework (7 independent sources, continuous evidence strength) |
| **Validation** | Computational (structural modeling + MR sensitivity analyses) | Leave-one-cohort-out CV (88.1% recovery); proteomics; spatial; cross-species |
| **Web resource** | Data/code at Zenodo | Streamlit interactive explorer (planned) |

### 4.2 Shared Methodological Components

| Method | Seagle et al. | Our Atlas | Key Differences |
|--------|--------------|-----------|-----------------|
| **S-PrediXcan** | 48 tissues, Bonferroni p < 2.03 x 10^-7, 160 new genes identified | Liver tissue focused; zero FDR-significant genes reported (flagged as weakness) | Their larger source GWAS (N=218K vs. our GWAS Ns) likely explains their success; our TWAS failure may reflect GWAS power rather than methodology |
| **Colocalization** | eQTL coloc (PP.H4 > 0.80) + pQTL coloc (PP.H4 > 0.70) | coloc.abf against Broadaway eQTLs; 348 genes PP.H4 > 0.5 | Different eQTL reference panels; different PP.H4 thresholds; different source GWAS |
| **MR** | Drug-proxied MR: tissue-specific GPGE as IV for drug effect on MASLD risk | Standard two-sample MR: gene expression as exposure, MASLD as outcome | Fundamentally different MR designs answering different questions |
| **Drug database mapping** | OpenTargets, DGIdb, ChEMBL, DrugBank | DGIdb (for LINCS annotation and druggability) | They use 4 databases; we use 1 primarily |
| **IPA pathway analysis** | Applied to S-PrediXcan genes | Not used (we use GSEA, ORA, ssGSEA, GSVA, decoupleR, COSMOS) | Different pathway analysis ecosystems |

### 4.3 Gene-Level Overlap Questions

The following should be checked against actual atlas results:

1. **Do FADS1, FADS2, FADS3 appear in our multi-evidence atlas?** If yes, what is their convergence ranking, dream LFC, cross-species concordance category, and COLOC PP.H4?
2. **Does S1PR2 appear?** They found it was the only gene with strong liver-tissue GPGE signal (p = 4.63 x 10^-9). Does it appear in our dream DEGs, COLOC results, or LINCS reversal hits?
3. **Does PPARG appear?** Given our deconvolution data, can we resolve their PPARG paradox by showing cell-type-specific PPARG expression patterns?
4. **How many of their 212 genes appear in our atlas?** What is the distribution across our evidence layers? This is a concrete external validation exercise.
5. **Does icosapent ethyl appear in our LINCS L1000 reversal hits?** If yes, at what WTCS rank? If no, this illustrates the complementarity of GWAS-based vs. transcriptome-based drug discovery.
6. **Does fingolimod appear in our LINCS results?**

### 4.4 TWAS Discrepancy

A significant methodological point: they successfully identified 160 genes via S-PrediXcan (Bonferroni p < 2.03 x 10^-7), while our TWAS module reports zero FDR-significant genes.

Possible explanations:

1. **GWAS power**: Their MVP GWAS (N=218K) is substantially larger than the GWAS we used (Ghodsian, Chen, etc.). TWAS/S-PrediXcan power is directly determined by the underlying GWAS power.
2. **Tissue scope**: They run across 48 tissues (increasing discovery space and test count, but also increasing multiple testing burden). We focus on liver tissue.
3. **GTEx version**: They use GTEx v7. We may use a different version.
4. **Multiple testing correction**: They use Bonferroni across 217,391 tests. We use FDR.

**Implication**: Our planned cTWAS upgrade should use the largest available GWAS (consider incorporating Vujkovic/MVP if summary statistics are accessible) and multiple tissues. Their result demonstrates that MASLD TWAS is feasible with sufficient GWAS power.

---

## 5. Actionable Implications for Our Atlas

### 5.1 High Priority

#### A. Check FADS1/FADS2/S1PR2/PPARG in Multi-Evidence Atlas

**Action**: Query `multi_evidence_atlas.csv` for these genes. Report:
- Dream LFC and padj (S1)
- Cross-species concordance category (S2)
- COLOC PP.H4 across all tested GWAS (S3)
- DepMap essentiality (S4)
- ATAC-seq accessibility (S5)
- Spatial expression pattern (S6)
- scRNA-seq pseudobulk DE (S7)
- Convergence ranking
- LINCS WTCS score if applicable

**Purpose**: Concrete convergence/divergence analysis. If FADS1 ranks high in our atlas through an entirely independent analytical path (diseased-tissue expression + LINCS reversal + cross-species concordance), this is strong cross-study validation. If it does not appear, this illustrates the complementarity we claim -- GWAS captures causal initiators, expression captures disease-responsive effectors.

#### B. Add Hepatotoxicity Screening to Drug Repurposing Module

**Action**: Cross-reference all LINCS-nominated compounds and multi-layer drug targets against:
- NCBI LiverTox database
- FDA DILIRank dataset (available at https://www.fda.gov/science-research/liver-toxicity-knowledge-base-ltkb/drug-induced-liver-injury-rank-dilirank-dataset)

**Rationale**: Seagle et al. found that 33/57 candidates (58%) had some hepatotoxicity signal, including 5 Class A (well-known hepatotoxicants). Reviewers will expect this filter in a translational drug repurposing analysis. This is low-effort, high-credibility.

**Implementation sketch**: Download DILIRank CSV; merge on drug name with our LINCS hits; add LiverTox/DILI class column to `multi_layer_drug_targets.csv`. Flag and optionally exclude Class A compounds.

#### C. Refine Discussion Paragraph on Seagle et al.

**Current text** (from `03_discussion.md`):
> "Seagle et al. (2025, AJHG) performed GWAS-only drug repurposing identifying targets including FADS1 and S1PR2. Our overlap analysis (3.6% DEG-COLOC overlap, Spearman rho = -0.031) demonstrates that such GWAS-only approaches capture fundamentally incomplete biology."

**Problem**: "GWAS-only" is reductive. They add protein structural modeling (RosettaLigand + AlphaFold3 + 100-ns MD), drug-proxied MR, pQTL integration, and hepatotoxicity screening. Reviewers familiar with the paper will note this.

**Suggested revision**:
> Seagle et al. (2025, AJHG) developed a GWAS-to-drug pipeline integrating S-PrediXcan, drug-proxied Mendelian randomization, and protein structural modeling to nominate icosapent ethyl (via FADS1) and fingolimod (via S1PR2) as MASLD prevention candidates. Their approach provides mechanistic depth on specific drug-target interactions but operates entirely within the GWAS signal space, using genetically predicted gene expression in healthy tissue (GTEx) rather than diseased-tissue transcriptomics. The near-zero overlap between our transcriptomic DEGs and COLOC-nominated causal genes (3.6%, Spearman rho = -0.031) demonstrates that GWAS-anchored approaches -- however methodologically sophisticated -- capture a fundamentally different and complementary biology to disease-responsive expression programs. [If confirmed: FADS1 ranks Xth in our convergence ranking / does not appear among our LINCS reversal hits, illustrating this complementarity.] Our atlas would identify their GWAS-nominated targets through the genetic causal evidence layer while providing the transcriptomic, epigenomic, and cell-type-resolved context necessary for therapeutic prioritization that their pipeline cannot access.

### 5.2 Medium Priority

#### D. Use Their 212-Gene List as External Benchmark

**Action**: Download their Table S3 (212 putative causal genes). Cross-reference against our multi-evidence atlas. Compute:
- Fraction appearing in atlas
- Distribution across evidence layers
- Mean convergence rank vs. random gene set
- Enrichment for specific evidence sources (are GWAS-derived genes enriched in our S3/causal layer but depleted in S1/transcriptomic layer?)

**Purpose**: Quantitative validation of the complementarity claim. This could become a supplementary figure showing that GWAS-nominated genes cluster in the genetic evidence layer but are underrepresented in transcriptomic and cross-species layers.

#### E. Resolve the PPARG Paradox with Deconvolution

**Action**: If PPARG shows cell-type-specific expression patterns in our MuSiC/BayesPrism results, this directly addresses Seagle et al.'s acknowledged limitation. Specifically:
- Does PPARG expression differ between hepatocyte-intrinsic and compositional components?
- Does PPARG appear in our sc-TWAS at cell-type resolution?
- Can we distinguish hepatocyte PPARG from macrophage/stellate PPARG effects?

**Rationale**: Their paper devotes ~2 paragraphs to the PPARG paradox but cannot resolve it because they lack cell-type resolution. If our deconvolution shows that the PPARG risk association is driven by hepatocyte-specific expression (consistent with mouse knockout data they cite), this is a concrete demonstration of the atlas's added value over GWAS-only approaches.

#### F. Consider pQTL MR as Additional Evidence Layer

**Action**: Evaluate whether incorporating deCODE pQTL data (4,907 circulating proteins, 35,559 Icelanders) is feasible as an additional evidence source in the multi-evidence atlas.

**Assessment**: This is a separate evidence modality from our existing eQTL-based COLOC. Circulating protein levels capture post-transcriptional regulation and secreted protein biology that eQTL-COLOC misses. However, circulating proteins are not tissue-specific, and many MASLD-relevant proteins (intracellular enzymes, membrane receptors) are poorly represented in blood-based assays. Seagle et al. identified only 3 additional genes from this approach, suggesting limited marginal yield. **Verdict**: Low priority unless the cTWAS pipeline naturally extends to pQTL integration.

### 5.3 Lower Priority / Monitoring

#### G. Evaluate MVP GWAS for cTWAS Upgrade

Their successful S-PrediXcan results (160 genes at Bonferroni significance) using the MVP GWAS suggest that incorporating Vujkovic et al. summary statistics into the planned cTWAS/SharePro/INTACT pipeline could resolve our current zero-FDR-significant TWAS weakness. Check whether MVP MASLD GWAS summary statistics are publicly available (dbGaP: phs001672.v7.p1).

**Caveat**: The MVP GWAS uses ALT-based phenotype definition, which introduces different biases than the GWAS we currently use. This could be an advantage (complementary signal) or a disadvantage (enzyme-kinetic confounding).

#### H. Drug-Proxied MR Design

Their MR design (using tissue-specific GPGE at drug targets as IVs to proxy drug effects) is methodologically distinct from standard two-sample MR and could be adopted for our top LINCS candidates. This would require:
1. Identifying the primary indication for each LINCS-nominated compound
2. Running S-PrediXcan against the primary indication GWAS
3. Using those GPGE results as IVs in MR with MASLD as outcome

**Assessment**: High effort, moderate yield. Worth considering if icosapent ethyl or other Seagle candidates also appear in our LINCS results -- in that case, we could directly validate their MR results using our independent analytical path.

#### I. Structural Biology Integration

Their protein structural modeling (RosettaLigand, AlphaFold3, MD simulations) represents a layer of mechanistic validation absent from our pipeline. This is appropriate for a drug-repurposing-focused paper but outside the scope of an atlas paper. **No action needed**, but worth acknowledging in discussion that structural validation of atlas-nominated drug targets is a natural follow-up.

---

## 6. Specific Gene Overlap Analysis Template

To be populated after checking actual atlas results:

### FADS1

| Evidence Layer | Seagle et al. Result | Our Atlas Result | Concordance |
|----------------|---------------------|------------------|-------------|
| GWAS association | MVP GWAS locus | Check S3 layer | |
| S-PrediXcan/TWAS | 23 significant tissue GPGEs | Check sc-TWAS | |
| COLOC | PP.H4 > 0.80 (multiple tissues) | Check coloc_results.csv | |
| Dream DEG | N/A (no expression data) | Check dream results | |
| Cross-species | N/A | Check concordance atlas | |
| LINCS reversal | N/A (not queried) | Check LINCS WTCS | |
| MR | OR = 0.28, p = 6.8 x 10^-160 | Check MR results | |
| Structural modeling | AlphaFold3 + RosettaLigand + 100-ns MD | N/A | |

### S1PR2

| Evidence Layer | Seagle et al. Result | Our Atlas Result | Concordance |
|----------------|---------------------|------------------|-------------|
| GWAS association | MVP GWAS locus | Check S3 layer | |
| S-PrediXcan (liver) | p = 4.63 x 10^-9, PP.H4 = 0.90 | Check sc-TWAS | |
| Dream DEG | N/A | Check dream results | |
| Cross-species | N/A | Check concordance atlas | |
| LINCS reversal | N/A | Check LINCS WTCS | |
| MR | OR = 0.91, p = 4.0 x 10^-4 | Check MR results | |
| Structural modeling | S1PR2-fingolimod docking (PDB: 7T6B) | N/A | |

### PPARG

| Evidence Layer | Seagle et al. Result | Our Atlas Result | Concordance |
|----------------|---------------------|------------------|-------------|
| GWAS association | MVP GWAS locus | Check S3 layer | |
| MR | OR = 1.93 (risk-increasing), p = 1.2 x 10^-14 | Check MR results | |
| Dream DEG | N/A | Check dream results | |
| Deconvolution | N/A (identified as limitation) | Check MuSiC/BayesPrism | |
| sc-TWAS | N/A | Check cell-type-specific effects | |
| Structural modeling | Pioglitazone cocrystal; best IDX of all ligands | N/A | |

---

## 7. Citation and Positioning Notes

### How to Cite in Introduction

Seagle et al. fits in the introduction alongside other GWAS-to-function studies as an example of genomics-informed drug repurposing for MASLD. It demonstrates the translational potential of GWAS data but also (through the PPARG paradox and ALT-phenotype limitations) illustrates why GWAS signals alone are insufficient for target prioritization.

### How to Cite in Discussion

Position as complementary: their genetic-to-drug depth vs. our multi-modal breadth. The key differentiating claim is that our atlas captures disease-responsive biology (AKR1B10, LPL, TREM2, etc.) that is invisible to any GWAS-anchored approach, while simultaneously incorporating GWAS-derived causal evidence as one of seven orthogonal layers. Their PPARG paradox serves as a specific example where cell-type resolution (which we have, they do not) is needed to interpret GWAS signals correctly.

### What Not to Claim

- Do not claim their work is methodologically inferior -- it is not. Their structural biology component (AlphaFold3, MD simulations) is more sophisticated than anything in our pipeline.
- Do not overstate our drug repurposing as "better" -- it is different. Their MR-based drug evaluation provides causal evidence that our LINCS reversal analysis does not.
- Do not dismiss the ALT-based phenotype as a fatal flaw -- they replicated in both histologic and radiologic cohorts.

---

## 8. Key Numbers Reference

For quick manuscript citation:

- 212 putative causal MASLD genes identified (158 = 74.5% novel)
- 412 significant gene-tissue pairs (p < 2.03 x 10^-7, PP.H4 > 0.80)
- 57 drug-gene pairs after directional concordance filtering
- FADS1 MR: OR = 0.28 (95% CI 0.26--0.31), p = 6.8 x 10^-160
- S1PR2 MR: OR = 0.91 (95% CI 0.86--0.97), p = 4.0 x 10^-4
- PPARG MR: OR = 1.93 (95% CI 1.62--2.31), p = 1.2 x 10^-14
- MVP GWAS: 90,408 cases / 128,187 controls
- deCODE pQTL: 4,907 proteins in 35,559 Icelanders
- 100-ns MD simulations (3 independent trajectories for FADS1-IPE)
- FADS1 IDX = -14.199, FADS2 IDX = -14.058

---

*Document generated for internal project reference. Not for distribution.*
*Date: 2026-04-14*
