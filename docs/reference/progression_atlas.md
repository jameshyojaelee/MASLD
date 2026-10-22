# Disease Progression Atlas (Consolidated)

> Sources: `docs/progression_executive_summary.md`, `docs/progression_review_findings.md`, `docs/progression_contrast_audit.md`. Optimized for LLM parsing.

## Current-use correction (2026-08-08; appended, older atlas preserved)

This document is a historical progression-pipeline reference. For the new
Figure 2 adjacent-fibrosis panel, the only accepted source is the isolated,
not-yet-adopted candidate
`RNA-seq/Human/Patient_Cohorts/analysis/integration/results/candidates/program-context-v2-candidate-2026-08-07/fibrosis-adjacent-true-kleiner-lvqw-v1/`.
It uses exact Kleiner stages in six cohorts, one pass-QC first biopsy per
GSE193066 donor, fixed-effect LVQW
`~ dataset + inferred_sex + fib_group`, BH across each complete 27,638-gene
family, and no LFC cutoff.

| Transition | Current candidate DEGs | Up / down | Donors (low / high); cohorts |
|---|---:|---:|---:|
| F0→F1 | **798** | 442 / 356 | 313 (126 / 187); 6 |
| F1→F2 | **2,190** | 1,510 / 680 | 361 (187 / 174); 6 |
| F2→F3 | **3,490** | 2,131 / 1,359 | 307 (174 / 133); 6 |
| F3→F4 | **783** | 504 / 279 | 152 (108 / 44); 5 |

Unique exact-stage biological donors are F0/F1/F2/F3/F4 =
126/187/174/133/44 (664 total); 105 are eligible GSE193066 first biopsies.
Job `19642832` passed 13/13 checks and an independent manifest/count/census
replay. `READY` SHA256 is
`f6a6d274124749b5c8bbdbf204629c847a33b45bc9d2a570ba98352d895fbcc9`;
gene-table SHA256 is
`01a942d33fb321631db0a3cb0579b2b89b96373c399d61f87030a49ed9788f18`.

All Script-115 Dream adjacent counts and June Script-14 counts below remain
visible for provenance but are retired from Figure 2. The June analysis mixed
binned-stage GSE213621 samples with exact stages and counted repeat GSE193066
biopsies independently. Transition-derived pathway, AUCell, carrier/routing,
CCC, drug-window, and cascade claims below were not recomputed from the clean
signatures and are not current manuscript evidence. Cross-sectional stage
contrasts must not be described as observed longitudinal progression.

## Pipeline Overview

Pipeline: Scripts 115-167 (Phases I-VI validated, Phase VII v2 corrected) + Scripts 130-146 (15-script contrast pipeline)

| Script | Phase | What it does | Key output | Status |
|--------|-------|-------------|------------|--------|
| 115 | I | Transition programs (dream DE per stage pair) | `results/progression/transition_programs.csv` (180K rows) | COMPLETE |
| 115b | I | Pathway enrichment (fgsea Hallmark+KEGG) | `results/progression/transition_pathway_enrichment.csv` (1,591 rows, 511 sig) | COMPLETE |
| 116 | I | Consensus pseudotime (NAS-VAE + BulkFormer + PCA) | `results/progression/consensus_pseudotime.csv` (1,444 samples) | COMPLETE |
| 117 | I | Bifurcation divergence (S1/S2 subtype fates) | `results/progression/divergence_genes.csv` (3K genes) | COMPLETE |
| 118a | II | Prepare scRNA reference (50K cells, 13 types) | Reference object | COMPLETE |
| 118b | II | BayesPrism deconvolution (Gibbs, 1,444 samples) | `results/progression/celltype_composition_changes.csv` | COMPLETE |
| 119 | II | Cell-type resolved progression (per-type DE) | `results/progression/celltype_transition_programs.csv` (38.6M) | COMPLETE |
| 124 | III | CCC L-R scoring + transition DE | `results/progression/ccc_transition_de.csv` (3.2M) | COMPLETE |
| 125 | III | CCC pseudotime dynamics (trajectory clustering) | `results/progression/ccc_pseudotime_dynamics.csv` (1.2M rows) | COMPLETE |
| 126 | IV | Driver scoring (6-dimension composite) | `results/progression/driver_scores.csv` (238K rows) | COMPLETE |
| 127 | IV | Spatial validation (5-test suite) | Validation metrics | COMPLETE |
| fig3 | V | 8-panel progression figure | `figures/fig3/fig3_progression.pdf` | COMPLETE |
| figS_ccc | V | CCC dynamics supplementary figure | Supplementary PDF | COMPLETE |
| 130 | VI | Core contrasts (14 total incl C11/C12/C13/C17) | Per-contrast DEG tables | COMPLETE (v3) |
| 131 | VI | NAS component ordinal + C17 | C7a-c + C17 results | COMPLETE (v3) |
| 132 | VI | Sex-stratified NASH-vs-NAFL | Interaction model results | COMPLETE |
| 133 | VI | Pathway enrichment (15 contrasts) | GSEA tables | COMPLETE (v3) |
| 134 | VI | Gene progression classification | 8 gene classes + tau | COMPLETE (v3) |
| 135 | VI | TWAS overlay | Ghodsian + Chen | COMPLETE |
| 138 | VI | Drug repurposing (CGP reversal) | 565 hits, 242 unique CGP sets | COMPLETE (v2) |
| 139 | VI | Deconvolution attribution | % hepatocyte-intrinsic per contrast | COMPLETE |
| 140 | VI | Cross-species concordance | Mouse model ranking | COMPLETE |
| 141 | VI | ATAC integration | Peak overlaps | COMPLETE |
| 142 | VI | TF activity (decoupleR) | Top progression TFs | COMPLETE |
| 143 | VI | Regulon-Drug-GWAS convergence | Triple convergence stats | COMPLETE |
| 144 | VI | Subtype x progression | S2 fibrosis OR=3.20 | COMPLETE (v2) |
| 145 | VI | Atlas assembly | 135-column progression atlas | COMPLETE (v4) |
| 146 | VI | Synthesis figures | 9 panels + 3 new figure scripts | COMPLETE (v4) |

### Phase VII v1 Scripts [DEPRECATED -- do not cite]

| Script | Issue |
|--------|-------|
| 147 | Permutation bugs V4/V5/V6 |
| 148 | Endothelial anti-correlated; stellate rho=0.12 |
| 149 | Built on inflated v1 risk scores |
| 150 | BulkFormer = VAE fallback; superseded by 160 |
| 151 | No nested CV; 0.865 superseded by 0.799 |
| 152 | AUROC 0.431 (below chance); circular pseudo-labels |
| 153 | CPS not independently validatable |
| 154 | Needs nested CV revalidation |
| 155 | Based on v1 numbers |

### Phase VII v2 Scripts [CURRENT]

| Script | Key result | Status |
|--------|-----------|--------|
| 160 | Feature engineering v2 (cell-type + expression + clinical) | COMPLETE |
| 160b | HCC score gene mapping fix: 56/59 mapped, rho=0.126 (p=2.0e-5) | COMPLETE |
| 161 | Nested LOCO classifier: M2 fib_ge3 AUROC 0.799 | COMPLETE |
| 162 | LOCO-NMF subtypes: cross-fold cosine >0.99 (batch correction leakage documented) | COMPLETE |
| 164 | TAS permutation null: 3/6 tests pass at p<0.05 | COMPLETE |
| 167 | ~~BulkFormer transfer~~ | [RETRACTED] VAE fallback |

## Key Results

### v2 Honest Numbers (Nested LOCO-CV)

| Metric | Value | Source | Assessment |
|--------|-------|--------|------------|
| F>=3 AUROC (cell-type, 17 feat) | **0.799** | Script 161 (M2) | Most defensible: no pre-screening, proper nested CV |
| F>=3 AUROC (expression, 500 div_* genes) | 0.853 | Script 161 (M3) | Pre-screened bias; does NOT beat random genes (p=0.73) |
| F>=3 AUROC (clinical, age+sex) | 0.660 | Script 161 (M1) | 2/6 folds degenerate |
| F>=3 AUROC (combined) | 0.827 | Script 161 (M4) | No improvement over M3 |
| S2 subtype AUROC (expression) | 0.904 | Script 161 (M3) | Random gene p=0.04 (barely significant) |
| TAS pseudotime ordering | rho=0.679, p=0.001 | Script 164 (V2) | PASS |
| TAS within-stage heterogeneity | p=0.001 | Script 164 (V4) | PASS |
| TAS pathway coherence | p=0.016 | Script 164 (V6) | PASS |
| TAS cell-type concordance | p=1.0 | Script 164 (V1) | FAIL |
| TAS fate probability | p=0.128 | Script 164 (V3) | FAIL |
| TAS LOCO cross-dataset | p=1.0 | Script 164 (V5) | FAIL |
| HCC score vs fibrosis | rho=0.126, p=2.0e-5 | Script 160b | 56/59 genes mapped |
| HCC Kruskal-Wallis | p=7.3e-4 | Script 160b | Differs by fibrosis group |
| HCC S2-vs-S1 | p=0.091 | Script 160b | Non-significant after correction |
| LOCO-NMF cross-fold cosine | >0.99 (S1 and S2) | Script 162 | Global concordance 0.717 |
| M2 top feature: fibrogenic_axis | mean |coef|=0.592 | Script 161 | Strongest predictor |

### v1 Retracted / Dropped [DO NOT CITE]

| Claim | v1 | Corrected | Reason |
|-------|-----|-----------|--------|
| ~~Cell-type AUROC~~ | 0.865 | 0.799 | Nested CV |
| ~~5-gene panel AUROC~~ | 0.863 | Needs revalidation | Nested CV required |
| ~~BulkFormer transfer AUROC~~ | 0.969 | [RETRACTED] | VAE fallback, not foundation model |
| ~~Pseudo-survival C-index~~ | 0.638 | [DEPRECATED] | CPS not independently validatable |
| ~~Personalized drug recs~~ | -- | [DEPRECATED] | Built on inflated v1 scores |
| ~~Cirrhosis risk AUROC~~ | 0.431 | [DROPPED] | Below chance, circular pseudo-labels |
| ~~HCC score rho~~ | 0.267 (p=6.5e-20) | 0.126 (p=2.0e-5) | 56/59 genes mapped (was 15-gene subset) |

### Contrast DEG Counts

| Contrast | N samples | N datasets | DEGs (padj<0.1, exploratory) | DEGs (padj<0.05, Tier 2) | DEGs (padj<0.05+\|LFC\|>0.3) |
|----------|-----------|------------|------------------------------|--------------------------|-------------------------------|
| C3 (Adv Fib F3-F4 vs F0-F2) | 1,077 | 7 | 14,486 | 12,758 | -- |
| C9 (F2 Inflection F2-F4 vs F0-F1) | 1,077 | 7 | 13,803 | 12,171 | 2,131 |
| C8 (Cirrhosis F4 vs F0-F3) | 609 | 6 | 11,720 | 9,813 | -- |
| C17 (Fibrosis ordinal) | 1,128 | -- | 11,508 | -- | -- |
| C11 (NASH vs Control) | ~609 | -- | 10,135 | 8,271 | -- |
| C12 (Early vs Late NASH) | ~160-200 | -- | 9,043 | 7,357 | -- |
| C6 (Extreme endpoints) | 234 | 5 | 8,023 | 6,578 | -- |
| C5 (NAS>=5 vs NAS<5) | 660 | 5 | 6,792 | 5,386 | -- |
| C4 (NAFL vs Ctrl) | 160 | 4 | 7,335 | 5,318 | -- |
| C13 (NASH vs NAFL fib-adj) | -- | -- | 497 | 268 | -- |

- At |logFC|>0.5 + padj<0.05, C3 drops to 1,153 DEGs
- Lambda values 9.2-10.2 for C3/C9 (high but consistent with C1 lambda=7.6 and large N)

### Positive Control Validation (Script 130)

| Gene | C3 logFC | C3 padj | C8 logFC | C8 padj |
|------|---------|---------|---------|---------|
| COL1A1 | +0.924 | 9.3e-36 | +0.899 | 5.9e-08 |
| ACTA2 | +0.461 | 5.6e-15 | +0.768 | 1.5e-08 |
| TGFB1 | +0.260 | 7.2e-09 | +0.347 | 0.003 |
| THRB | -0.228 | 2.7e-16 | -0.365 | 2.6e-09 |
| PPARA | -0.164 | 3.8e-11 | -0.169 | 0.003 |
| NR1H4 | -0.097 | 1.9e-06 | -0.127 | 0.003 |

All 7 positive controls validated with correct direction. Direction consistency: 98.5-100% among fibrosis contrasts, 63% for C3 vs C4 (biologically correct -- different biology).

### Phase I-VI Key Numbers

> **Historical pipeline values; not current Figure 2 inputs.** In particular,
> the Script-115 adjacent-transition rows below are superseded by the clean
> 2026-08-08 candidate table at the top of this document.

| Metric | Value | Source |
|--------|-------|--------|
| Transition DEGs: F0-F1 | 4,777 | Script 115 |
| Transition DEGs: F1-F2 | 2,265 | Script 115 |
| Transition DEGs: F2-F3 | 3,288 | Script 115 |
| Transition DEGs: F3-F4 | 2,244 | Script 115 |
| OT cost: F0-F1 / F1-F2 / F2-F3 / F3-F4 | 0.559 / 0.665 / 0.759 / 0.930 | Script 115 |
| S2 patients | 345 (27.5%) | Script 117 |
| S2 female fraction | 73.7% (p=8.2e-9) | Script 117 |
| S2 vs S1 fibrosis | F3 vs F2 (p=5.4e-14) | Script 117 |
| Divergence genes | 1,797 (210 with \|Cohen's d\|>0.5) | Script 117 |
| S2 ECM signature | COL15A1, CDH11, LTBP2, COL16A1, COL11A1 | Script 117 |
| CCC rewired L-R pairs | 4,800 (padj<0.1, \|LFC\|>0.25) | Script 124 |
| CCC feasible pairs | 381 of 1,586 | Script 124 |
| Most rewired axis (F2-F3) | Stellate-to-Macrophage (374 sig pairs) | Script 124 |
| GAS6::MERTK at F3-F4 | LFC=-1.0, padj=4.9e-8 | Script 124 |
| THRB driver rank (F3-F4) | 136 / 34,047 (top 0.4%) | Script 126 |
| Drug reversal hits | 565 total, 242 unique CGP, 85 progression-specific | Script 138 |
| Macrophage-dependent genes (C2) | 293 | Script 139 |
| Progression-specific CC genes | 829 | Script 140 |

## Code Review Findings (from progression_review_findings)

### 11 Publishable Findings (Ranked)

1. **Vast majority of NASH-vs-NAFL DEGs are fibrosis-confounded** (C13, Script 130): C2 unadjusted vs C13 fibrosis-adjusted at padj<0.05: 268 DEGs survive adjustment. True NASH core: TREM2, LPL, FABP4, MMP9, IL1B.
2. **S2 = fibrogenesis subtype, orthogonal to NASH** (Script 144): NASH/NAFL chi2 p=1.0 (no association). F3-F4 vs F0-F2 chi2=60.09 (p=9.06e-15, Cramer's V=0.242). OR=3.20 per fibrosis stage (p=1.15e-17). S2 fraction: F0=17%, F1=19.3%, F2=23.9%, **F3=43.3%**, F4=50%.
3. **FPC overtakes MCD as best mouse progression model** (Script 140): FPC progression rho=0.350 vs MCD 0.329. MCD better for onset (0.314 vs 0.264). Concordance rates increase onset->progression (~67% -> ~72%). Peak: MCD for NAS>=5 at 80.7%.
4. **STAT3 fibrosis-specific, STAT1 inflammation-specific** (Script 142): STAT3 NOT sig in C2 (NASH-vs-NAFL, padj=0.18), 7th strongest in C3 (padj=2.9e-25).
5. **NAFL-to-NASH transition is sex-uniform** (Script 132): 1 interaction DEG (ANOS2P, pseudogene). pi0=0.9919, SD(t)=1.015. 823 samples, NOT underpowered. 0 sex-divergent genes vs 591 at onset. 100% direction concordance (1,657/1,657). Female:Male DEG ratio drops from 3.84x (onset) to 1.79x (progression).
6. **HNF4A/PPARG drug reversal is progression-specific** (Script 138): SUMI_HNF4A_TARGETS reversal in all 4 progression contrasts (NES -2.04 to -2.24), absent from onset. PPARG (pioglitazone) reversal progression-specific (C3/C5/C6, not C1).
7. **Fibrosis shifts from transcriptional to compositional** (Script 139): Among DEGs: C1 ~29.8% hepatocyte-intrinsic (1,635 genes), C2 83.8%, C3 60.7% (39.3% composition-driven). [Note: C1 81.5% was from earlier pipeline; current verified number is 1,635 / 5,484 = 29.8%. C2/C3 percentages use per-contrast denominators and may differ.]
8. **829 progression-specific Conserved_Core genes** (Script 140): Cross-species conserved progression targets NOT in original C1 set. C9 (F2 inflection) yields largest CC set (2,000 genes).
9. **NF-kB pathway dominates progression** (Script 142): 5 family members (NFKB1, NFKB2, RELA, RELB, REL) all activated in 3/3 contrasts. NFKB1 largest effect (+0.548).
10. **85 progression-specific drug reversal signatures** (Script 138): Mitochondrial bioenergetics (MOOTHA_VOXPHOS NES=-2.41), TGF-beta fibrosis (COULOUARN), HNF4A loss.
11. **953 true multi-contrast progression genes** (Script 134): Of 4,154 "progression_only" genes, 3,052 sig only in C3 -> actually fibrosis-driven. Only 953 are multi-contrast progression genes.

### Per-Module Validation Status

| Script | Verdict | Key caveat |
|--------|---------|-----------|
| 130 (Contrasts) | PASS | DEG counts inflated by power (52% at padj<0.1) -- not a bug |
| 131 (NAS ordinal) | PASS | C7a underpowered (N=76); ~24% ballooning DEGs have \|logFC\|<0.1 |
| 132 (Sex-stratified) | PASS | Clean null is the finding |
| 133 (GSEA) | PASS | C7a unidirectional enrichment may reflect single-cohort confounding |
| 134 (Classification) | PASS | "progression_only" class conflates fibrosis-driven with true multi-contrast |
| 138 (Drugs) | PASS | Strongly publishable |
| 139 (Deconv) | PASS | CYP2E1 NAFL-to-NASH signal is composition-driven (macrophage-modulated, interaction p=0.02) |
| 140 (Concordance) | PASS | Strongly publishable |
| 142 (TF activity) | PASS | THRB absent from DoRothEA regulons (known limitation for nuclear receptors) |
| 143 (Convergence) | PASS | 61 TFs with dual GWAS+Drug convergence; triple convergence 2.48x enrichment (p<1e-4) |
| 144 (Subtype) | PASS | Strongly publishable |

### NAS Components (Script 131)

Steatosis 1,812 DEGs << Inflammation 5,202 << Ballooning 5,377. Single-dataset GSE130970 (N=76). IL1B/CCL2/CXCL8 validated; FASN/SCD/DGAT2 correct direction but underpowered.

### Top 10 Progression TFs (Script 142, sig 3/3 contrasts)

| TF | Dir | Mean Diff | Min padj | Role |
|----|-----|-----------|----------|------|
| TCF7L1 | UP | +0.162 | 3.9e-40 | Wnt signaling |
| NFKB2 | UP | +0.347 | 1.0e-34 | Non-canonical NF-kB |
| TCF7L2 | UP | +0.329 | 9.5e-33 | Wnt/metabolic (T2D locus) |
| SMAD1 | UP | +0.272 | 1.8e-32 | BMP signaling |
| SNAI1 | DOWN | -0.084 | 3.4e-26 | EMT repressor |
| NR1H4 | DOWN | -0.115 | 4.3e-20 | FXR (bile acids) |
| SMAD3 | UP | +0.298 | 9.3e-20 | TGF-beta/fibrosis |
| RELA | UP | +0.480 | 8.3e-18 | Canonical NF-kB |
| FOXM1 | UP | +0.299 | 1.1e-15 | Proliferation |
| NFKB1 | UP | +0.548 | 8.3e-15 | Canonical NF-kB |

Key TF patterns: STAT3 fibrosis-specific (not sig C2, 7th in C3). PPARA progressively lost (-0.121 C2 -> -0.261 C3). SMAD3 escalates (0.173 -> 0.262 -> 0.459). Hepatocyte identity loss: NR1H4, HNF1A, RXRA, BCL6, FOXK1.

### Convergence Permutation (Script 143)

| Test | Observed | Expected | Fold | p |
|------|----------|----------|------|---|
| GWAS enrichment | 448 | 210.7 | 2.13x | <1e-4 |
| Drug enrichment | 1,228 | 505.0 | 2.43x | <1e-4 |
| Triple convergence | 100 | 40.4 | 2.48x | <1e-4 |

Top TFs: FOXA1 (20 triple targets, Drug OR=6.04), CEBPA (Drug OR=6.56, 12 triple targets), HNF4A (C3-only -- fibrosis-specific).

### Pathway Contrasts (Script 133, selected NES)

| Pathway | C1 (MASLD) | C4 (Onset) | C3 (Fib) | C8 (Cirr) | C7b (Inflam) |
|---------|-----------|-----------|---------|-----------|-------------|
| FATTY_ACID_METABOLISM | +1.28 | -1.51 | -2.13 | -2.28 | +1.42 |
| EMT | +2.85 | +0.84 ns | +3.28 | +3.33 | +3.08 |
| OXIDATIVE_PHOSPHORYLATION | +1.23 | -1.94 | -1.42 | -2.16 | +2.36 |
| INFLAMMATORY_RESPONSE | +2.12 | -1.47 | +2.23 | +2.39 | +2.86 |
| BILE_ACID_METABOLISM | -1.39 | -1.44 | -2.27 | -2.37 | -1.54 |

### Gene Classification (Script 134)

| Class | N genes | Protein-coding | lncRNA | Mean tau |
|-------|---------|---------------|--------|----------|
| ubiquitous | 8,908 | 6,229 | 2,161 | 0.547 |
| not_significant | 8,550 | 2,394 | 3,932 | 0.832 |
| onset_only | 6,902 | 2,814 | 2,855 | 0.777 |
| both_onset_and_progression | 4,161 | 2,055 | 1,539 | 0.655 |
| progression_only | 4,154 | 2,048 | 1,486 | 0.678 |
| fibrosis_specific | 1,571 | 579 | 670 | 0.759 |
| inflammation_specific | 192 | 119 | 56 | 0.656 |
| extreme_only | 109 | 74 | 22 | 0.648 |

## Contrast Audit

### Redundancy Matrix (Spearman rho of logFC, N=21,703 genes)

```
                 C1     C2     C3     C4     C5     C6     C8     C9
C1_MASLD         1.000  0.658  0.680  0.391  0.778  0.588  0.541  0.790
C2_NASH_NAFL     0.658  1.000  0.684 -0.055  0.789  0.813  0.578  0.751
C3_AdvFib        0.680  0.684  1.000  0.108  0.710  0.840  0.859  0.930
C4_NAFL_Ctrl     0.391 -0.055  0.108  1.000  0.164  0.044  0.046  0.170
C5_NAS5          0.778  0.789  0.710  0.164  1.000  0.764  0.567  0.789
C6_Extreme       0.588  0.813  0.840  0.044  0.764  1.000  0.703  0.831
C8_Cirrhosis     0.541  0.578  0.859  0.046  0.567  0.703  1.000  0.799
C9_F2Inflect     0.790  0.751  0.930  0.170  0.789  0.831  0.799  1.000
```

C7 ordinal (N=23,187): C7a-C7b rho=0.633, C7a-C7c rho=0.669, C7b-C7c rho=0.751

### New Contrasts (Added 2026-03-30)

| ID | Definition | DEGs (padj<0.05) | Relationship |
|----|-----------|-------------------|-------------|
| C11 | NASH vs Control | 8,271 | rho=0.799 with C1 |
| C12 | Early vs Late NASH | 7,357 | rho=0.889 with C3 |
| C13 | NASH vs NAFL (fib-adj) | 268 | Dramatic reduction from C2 |
| C17 | Fibrosis ordinal | 11,508 | rho=0.962 with C3 |

### Recommended Portfolio

**Primary (main results)**:
- C1 (MASLD vs Control): disease presence
- C2 (NASH vs NAFL): disease severity/inflammation
- C3 (F3-F4 vs F0-F2): fibrosis progression [PRIMARY]
- C4 (NAFL vs Control): disease initiation [ORTHOGONAL, rho=0.04-0.39]
- C5 (NAS>=5 vs NAS<5): clinical threshold
- C10 (Sex-stratified NAFL-NASH): sex biology
- C11 (NASH vs Control): pure NASH signature [NEW]
- C12 (Early vs Late NASH): within-NASH progression [NEW, novel]

**Secondary (supplementary/sensitivity)**:
- C6: Extreme endpoints (rho=0.840 with C3)
- C8: Cirrhosis binary (rho=0.859 with C3, lower N)
- C9: F2 inflection hypothesis (rho=0.930 with C3)
- C7b-c: NAS components (single-dataset, GSEA-level only)
- C17: Fibrosis ordinal (rho=0.962 with C3, dose-response)

**Removed**: C7a (6 DEGs, N=76, underpowered)

### C4 Hidden Star

C4 anti-correlated with C2 (rho=-0.055): genes UP at onset are slightly DOWN during progression. Strongest quantitative evidence for "onset != progression" thesis.

## V1 vs V2 Overhaul

| v1 Bug | v2 Fix |
|--------|--------|
| No nested CV -> inflated AUROCs | Proper nested LOCO-CV for all models |
| BulkFormer = VAE fallback | Retracted entirely |
| V4/V5/V6 permutation bugs | Fixed -> 3/6 pass (was 6/6 with bugs) |
| HCC score used 15-gene subset | Expanded to 56/59 genes |
| Cirrhosis risk circular pseudo-labels | Dropped (AUROC 0.431) |
| CPS regression not validatable | Deprecated |

Final: M2 cell-type AUROC **0.799** (was 0.865). HCC rho=**0.126** (was 0.267). M3 expression (0.853) does NOT beat 500 random genes (p=0.73).

## Biological Findings

### k-Program Bifurcation (FINAL k=6 on 2026-04-21 after comprehensive confounder strip)

**Prior k=7 on partial strip SUPERSEDED** — contained artifact programs (P4 Y-pseudogene leakage + P7 non-hepatic contamination). Legacy caches preserved for A/B at `nmf_results_cache_nosex.rds`. Current primary cache is `nmf_results_cache_clean.rds`.

Clean pool built with GENCODE v49 metadata: dropped chrY (full chromosome, 46 genes), chrM (34), rRNA biotype (68), RPS*/RPL* subunits (476), MRPS*/MRPL* (73), IG[HKL][VJCD] variable/constant (69), hemoglobin (2), X-inactivation escapees (11) → 776 genes from 34,453 → 33,677 → top-5000 IQR.

**Rubric 6/6 pass at k=6**. Labels are directional along the MASLD spectrum (Option A); only "Fibrogenic" uses gene-level biology.

| Code | Label | N (% dom) | stage_ρ | Top genes | F2 dominance trend |
|------|-------|-----------|---------|-----------|---|
| P1 | **Progression-Inflammatory** | 236 (18.8%) | **+1.0** | AKR1B10, CA12, CCL20, TREM2, HKDC1, S100A9, LPL, TNFRSF12A, OLR1 | 3.2→5.5→**20.0**→34.1→39.1 (3.6× F1→F2 jump) |
| P2 | **Quiescent-Parenchyma_1** | 332 (26.5%) | **−1.0** | TIMD4, CXCR1/2, FCGR3B, NCAM2, S100A8-9, CPA3 | 42.6→41.7→24.5→8.9→3.1 |
| P3 | **Quiescent-Parenchyma_2** | 117 (9.3%) | −0.9 | SCX, CALML6, CCDC188, CLDN9 | 24.5→10.3→11.3→5.9→3.1 |
| P4 | **Stable** | 178 (14.2%) | −0.4 | MUC17, LINC00951, GABRG2, RBFOX1 | 4.3→15.5→14.8→11.1→0 |
| P5 | **Progression-Other** | 261 (20.8%) | +0.7 | ENSG-dominated lincRNA | 16→22→21→21→23 (flat) |
| **P6** | **Fibrogenic** | 130 (10.4%) | +0.7 | **DES, CNN1, MYOCD, FLNC, FHL5, PDPN, PTPRQ, ACTG2, JPH2, CASQ2** | 9.6→4.5→8.4→18.9→31.2 (2.25× F2→F3) |

**Legacy S1/S2 alias**: `nmf_subtype` column in nmf_assignments.csv retained as S1 (non-Fibrogenic) vs S2 (Fibrogenic dominant). All legacy downstream scripts (144, 146, 155, 215, fig5 panel D) work with this 2-way view. New work uses `dominant_program` / `dominant_program_code`.

**F2 switch story (k=6)**: At F4, Progression-Inflammatory (39%) + Fibrogenic (31%) + Progression-Other (23%) dominate 93%. At F0, Quiescent-Parenchyma_1 + _2 (67%) dominate. The inflammatory axis makes a sharp F1→F2 jump (3.6×); the fibrogenic axis activates F2→F3 (2.25×). Multinomial LOCO-CV best: P1 AUROC=0.739. P2 most COLOC-enriched (OR=1.71 p=5.3e-4).

**Sanity**: 0 confounder hits (chrY/MT/RPS/RPL/IG/HB/XIST) across all 6 programs' top-50 loadings. Max pairwise |H correlation| = 0.39 (orthogonal). Max |ρ_sex_F| = 0.17 (sex leakage eliminated).

### Stellate Cell Expansion

| Stage | Stellate % | Change |
|-------|-----------|--------|
| F0 | 1.7% | -- |
| F1 | 2.1% | +24% |
| F2 | 2.9% | +38% |
| F3 | 4.2% | +45% |
| F4 | 8.5% | +102% |

All transitions padj<0.004. 5-fold total expansion. Stellate dominates DEGs: F1-F2 99.4% (3,533/3,553), F2-F3 86.4% (4,253/4,923). Hepatocyte surpasses stellate only at F3-F4.

### CCC Rewiring

| Transition | Most rewired axis | N sig pairs | Top L-R pair | LFC |
|-----------|-------------------|-------------|-------------|-----|
| F2-F3 | Stellate-to-Macrophage | 374 (all increasing) | CCN1::CAV1 | +1.67 |
| F2-F3 | Stellate-to-Macrophage | -- | CXCL8::SDC2 | +2.12 |
| F2-F3 | Stellate-to-Macrophage | -- | CCN2::IGF2R | +1.70 |
| F3-F4 | Hepatocyte-to-Macrophage | 312 (all decreasing) | GAS6::MERTK | -1.0 (padj=4.9e-8) |
| F3-F4 | Hepatocyte-to-Macrophage | -- | FCN2::LRP1 | -1.22 |
| F3-F4 | Hepatocyte-to-Macrophage | -- | LIPC::LRP1 | -1.07 |

14 canonical pathways: dominant decreasing trend across pseudotime. HGF/MET (rho=-0.256) and Interleukin (rho=-0.247) decline most.

### Driver Scores (Script 126)

| Gene | Best transition | Rank / 34,047 | Top % | Causal class |
|------|----------------|--------------|-------|-------------|
| THRB | F3-F4 | 136 | 0.4% | Robust (3+ methods) |
| THRB | F2-F3 | 148 | 0.4% | Robust |
| THRB | F1-F2 | 157 | 0.5% | Robust |
| NR1H4 (FXR) | -- | -- | 4.8% | -- |
| PPARA (elafibranor) | -- | -- | 6.6% | -- |

6 dimensions: transition specificity, genetic causality, OT cost, bifurcation divergence, cell-type resolution, cross-species conservation.

### Mouse Model Ranking (Script 140)

| Diet | Onset rho | Progression rho | Delta |
|------|-----------|-----------------|-------|
| FPC | 0.264 | **0.350** | +33% |
| MCD | **0.314** | 0.329 | +5% |
| CDAHFD | 0.244 | 0.318 | +30% |
| LIDPAD | 0.224 | 0.296 | +32% |
| HFD | 0.204 | 0.236 | +16% |

### Deconvolution Attribution (Script 139)

| Contrast | % Hepatocyte-intrinsic | % Composition-driven |
|----------|----------------------|---------------------|
| C1 (Disease vs Control) | ~29.8% (1,635 genes) | ~70.2% | [Note: was 81.5%; corrected per RESULTS_CATALOG verification 2026-04-14]
| C2 (NAFL vs NASH) | 83.8% | 16.2% |
| C3 (Advanced fibrosis) | 60.7% | 39.3% |

## Output Locations

- Phase I-V: `RNA-seq/Human/Patient_Cohorts/analysis/integration/results/progression/`
- Phase VII v2: `results/prognosis_v2/` (25 files)
- Phase VII v1: `results/prognosis/` [DEPRECATED], `results/prognosis_v2/bulkformer_transfer_results_RETRACTED.csv` [RETRACTED]
- Figures: `figures/fig3/fig3_progression.pdf`, `figures/progression/panel_{c13,nash,landscape}_*.pdf`
- Key output files listed in Pipeline Overview table above (Key output column)

### Caveats

- C7a/b/c unidirectional GSEA may reflect single-cohort confounding
- THRB absent from DoRothEA; PPARA/NR1H4 limited by small regulon sizes
- Apply |logFC|>=0.25 for ballooning DEG interpretation
- ACTA2: classified "progression_only" but functionally fibrosis-specific (sig C3/C6/C8/C9, not C1/C2/C4/C5)
