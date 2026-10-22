> **HISTORICAL CAMPAIGN LEDGER — SUPERSEDED 2026-08-10.** This file contains
> July figure and narrative snapshots, including values later corrected or
> challenged. It is not the current numerical authority and must not supply
> Resource prose, figures, portal defaults, or Cas13 statements. Use
> `docs/manuscript/NUMBERS.md` for named-release numbers and Plan 60 for the
> future immutable Resource ledger.

# Campaign Numbers Ledger — historical 2026-07-13 snapshot

**Purpose.** ONE locked numbers table for the manuscript-polish campaign. All figure-writing teams and the synthesis pass key off this file. Values are frozen against active release **`2026-07-10-r1`** (single-cell rates carry the 2026-07-12 donor-correction, pending formal r2).

**Authority precedence.** (1) `docs/manuscript/release/claim_ledger.tsv` (wording) → (2) `NUMBERS.md` active-release header `2026-07-10-r1` → (3) `docs/paper_outline.md` Frozen Evidence Contract → (4) source CSVs. Where these conflicted, the drift resolutions below were locked by the team lead.

**Spot-check legend.** ✅ = independently reproduced from the source file this pass; ⧗ = ledger value, not re-derived this pass; n/a = definitional/narrative.

---

## Cross-figure shared constants (appear in ≥2 figures)

| Quantity | Canonical value (verbatim) | Source file/section | Spot-check | Notes / drift resolution |
|---|---|---|---|---|
| Canonical bulk DEGs | **1,918 significant DEG records / 1,915 unique nonempty gene symbols** (1,419 up / 499 down) | `RNA-seq/Human/.../results/integration/canonical_deg_results.csv`, gate `treat_fdr < 0.05` at lfc = 0.25 | ✅ 1918 records; 1419 up / 499 down; 1915 unique symbols | Class/symbol-level analyses use **1,915**; record-level counts use 1,918. Method = limma-voom-qw C2, `~ dataset + inferred_sex + group_binary`, dataset FIXED. |
| Genes tested (bulk DE) | **27,638** | same file, all rows | ✅ 27,638 data rows | Unchanged under C2 swap. |
| Primary SuSiE genetic set | **473 genes** (PP.H4.susie > 0.5); **447 jointly testable**; **34 overlap TREAT DEGs**; **92.4% not differentially expressed** | `GWAS/finemapping/results/susie_coloc/gene_level_coloc_tier12.csv`; `RNA-seq/results/manuscript_release/2026-07-10-r1/orthogonality_audit.tsv` (row `all / susie / 0.5`) | ✅ 473 (tier12); 473/447/34, 92.394% (audit) | Tier-1/2 (35 liver GWAS) headline. |
| Named SuSiE/ABF union | **1,030 named genes**; **971 jointly testable**; **74 overlap TREAT DEGs**; **92.4% not differentially expressed** | tier12 file (1,031 ensembl IDs) → `orthogonality_audit.tsv` (row `all / susie_or_abf / 0.5`) | ✅ 1,031 ensembl-level; 1,030/971/74 named (audit) | **DRIFT LOCK: 1,030 (named).** Historical 1,031 = ensembl-level, included one blank-symbol row; **1,031 is retired as the headline**. |
| Overlap NOT enriched | Primary SuSiE Fisher **OR 0.89, p = 0.60**; union **OR 0.89, p = 0.37** | `orthogonality_audit.tsv` | ✅ OR 0.8895 p 0.604; OR 0.8877 p 0.371 | Frame: "largely non-overlapping," NOT "statistically/biologically independent." |
| Human sample counts | **1,259 QC-passing** across **9 cohorts**; **846 samples / 5 control-bearing cohorts** enter the canonical pooled disease-vs-control model (1,277 collected total) | `config/human_datasets.yaml`; `NUMBERS.md` Sample Counts | ⧗ | PRJNA512027/Gerhard dropped (L0/S0 batch confound). |
| Convergence Tier-1 | **677 genes** (Tier-1 Genetic_validated) | `RNA-seq/results/multi_evidence/convergence_evidence.csv`, col `tier` | ✅ 677 (also 204 Tier-3, 25,668 Tier-4, 638 Excluded; sum 27,187) | Tier-1/2-gated 6-channel human-only build. Supersedes full-portfolio 1,018. |
| Atlas dimensions | **27,187 genes × 443 columns** | `RNA-seq/results/multi_evidence/multi_evidence_atlas.csv` | ✅ 443 cols; 27,187 rows | **DRIFT LOCK: 443** (donor-corrected live atlas). NOT 411/412/197/156/84/218/444. |
| Gene count (atlas/convergence universe) | **27,187** | atlas + convergence (677+204+25,668+638) | ✅ both = 27,187 | |
| Jointly-testable background | **14,931 genes** (genetic × bulk overlap universe); TREAT DEGs jointly testable = **1,261** | `orthogonality_audit.tsv` (`n_joint_testable`, `n_treat_deg_joint`) | ✅ | Denominator for evidence-class Fisher tests. |

---

## Fig 1 / Abstract + Intro (Atlas construction, complementary-maps thesis)

| Quantity | Canonical value (verbatim) | Source file/section | Spot-check | Notes / drift resolution |
|---|---|---|---|---|
| QC-passing samples / cohorts | 1,259 / 9 (human bulk RNA-seq) | `config/human_datasets.yaml` | ⧗ | |
| Canonical model / gate | limma-voom-qw C2, `~ dataset + inferred_sex + group_binary` (dataset FIXED); gate `TREAT FDR < 0.05` at `lfc = 0.25` | Frozen Evidence Contract §Human map | n/a | Method-id `limma_voom_qw__C2`. Do NOT print the term "TREAT" in prose (see banned wording). |
| Canonical DEGs | 1,918 records / 1,915 unique symbols (1,419 up / 499 down) | `canonical_deg_results.csv` | ✅ | shared constant |
| Atlas dimensions | 27,187 × 443 | `multi_evidence_atlas.csv` | ✅ | shared constant |
| Independent evidence sources | 7 (S1–S7) | atlas schema | n/a | |
| Primary orthogonality (Fig 1B) | 473 SuSiE / 447 joint / 34 overlap / **92.4% non-DE**; union 1,030 / 971 / 74 / 92.4% | `orthogonality_audit.tsv` | ✅ | Central result of the paper. |
| DEG–COLOC overlap | **NOT enriched** vs joint-testable background (SuSiE OR 0.89 p 0.60; union OR 0.89 p 0.37) | `orthogonality_audit.tsv` | ✅ | The strict Tier-1 set is not COLOC-enriched. Do NOT revive the retracted Fisher-OR-1.67 / 3.6% framing. |
| LOO-CV recovery | 81.6% (padj<0.1) / 80.8% (padj<0.05); mean Spearman ρ 0.952; direction concordance 99.9%; robust genes (5/5) **7,842** | `RNA-seq/Human/.../loo_cv_C2/loo_cv_C2_summary.csv` | ⧗ | Supersedes 8-fold 88.1%/7,128; and the non-comparable "0.895/0.735" 8-cohort dream AUROC. |

---

## Fig 2 (Liver-trait regulatory genetic map)

| Quantity | Canonical value (verbatim) | Source file/section | Spot-check | Notes / drift resolution |
|---|---|---|---|---|
| Primary GWAS portfolio | **35 Tier-1/2 liver GWAS strata**, 5 ancestries (EUR 17 / AFR 6 / EAS 6 / AMR 3 / SAS 3; 33 studies, 6 traits) | `GWAS/finemapping/config/gwas_trait_tier.tsv`; NUMBERS §2 | ⧗ | 50-GWAS full portfolio + Tier-3/4 = supplement only. Fig 2A study view = **48 distinct studies**. |
| SuSiE-COLOC (main) | **473 / 369 / 267** at PP.H4.susie > 0.5 / 0.8 / 0.9 | `gene_level_coloc_tier12.csv` col `coloc_best_susie_pp4` | ✅ 473 / 369 / 267 | supp full-50: 736 / 569 / 416 |
| ABF-COLOC (main) | **863 / 375 / 234** at PP.H4 > 0.5 / 0.8 / 0.9 | same file col `coloc_best_abf_pp4` | ✅ 863 / 375 / 234 | supp full-50: 1,234 / 547 / 316 |
| Named union (SuSiE ∪ ABF > 0.5) | **1,030 named genes** (1,031 ensembl-level) | tier12 union → `orthogonality_audit.tsv` | ✅ 1,031 ensembl; 1,030 named | supp full-50: 1,527. **Headline = 1,030 named; 1,031 retired.** |
| Direct-disease / PDFF subset | SuSiE 0.5 = **22 genes**; union = **172 named** (160 joint / 11 overlap / 93.1% non-DE) | `orthogonality_audit.tsv` (rows `direct_disease`) | ✅ 22; 172/160/11, 93.125% | Tier-1 = direct NAFLD/NASH/PDFF; Tier-2 = ALT/AST/GGT (liver-enzyme evidence, not direct susceptibility). |
| Coding lead-variant fraction | 34 / 1,030 (~3.3%) coding | `coloc_variant_annotation.csv` | ⧗ | **Descriptive locus annotation ONLY.** Do NOT frame as genome-wide architecture (see Coding-architecture gate, Fig 5/Discussion). |
| Cross-ancestry (main Fig 2G, 473 SuSiE set) | EUR-only 356 / shared 81 / non-EUR-unique 36 (EAS 26 / AFR 5 / AMR 5); significance-gated non-EUR discovery = **14 genome-wide-sig (p<5e-8)** / 20 suggestive (p<1e-6) | `scripts/figures/fig2_ancestry_unique_coloc_gated.py`; NUMBERS §2 | ⧗ | Ungated base panel RETIRED. Union-based "168 non-EUR" is superseded/supplementary. Non-EUR GWAS use EUR liver-eQTL (`cross_ancestry` flag = LABEL not gate). |
| Fig 2E boundary genes | PNPLA3, MTARC1 (protein-altering / protein-level); MTARC1 recovered by Gobeil 2026 pQTL (PPH4 = 1.00 vs eQTL PPH4 = 0.003), sits in `neither` class of expression map | paper_outline §Fig2 | n/a | |
| Broadaway overlap | 269 / 751 shared (36%); apples-to-apples coloc.abf≥0.7: 191/378 (51%) | `broadaway_coloc_comparison.csv`; NUMBERS §2b | ⧗ | Positioning vs source eQTL paper. |

---

## Fig 3 (Established disease-state transcriptome; stage-ordered remodeling)

| Quantity | Canonical value (verbatim) | Source file/section | Spot-check | Notes / drift resolution |
|---|---|---|---|---|
| Pooled DEG result (3A) | 1,918 records / 1,915 symbols (1,419 up / 499 down; 27,638 tested) | `canonical_deg_results.csv` | ✅ | shared constant |
| Fibrosis-adjustment attenuation | **2,724 → 268 = 90.2% reduction** (matched 615 samples, padj<0.05); envelope **87.3–96.1%** (padj 0.10/0.05/0.01) | `claim_C1_mediation.md` §5; NUMBERS §3 | ⧗ | **Never cite 95.5%** (composite of sample-set + fibrosis effects). Never cite 10,135 (that is C11 NASH-vs-Control, different contrast). |
| True-NASH core | 437 genes (C2 padj<0.1 ∩ C13 padj<0.1); TREM2, LPL, FABP4, MMP9, IL1B | NUMBERS Fig 3 | ⧗ | |
| Adjacent-stage contrasts | Label as `F1 versus F0` (adjacent-stage), NOT transitions in individual patients | paper_outline §Fig3 | n/a | Cross-sectional. Do NOT say "multicellular cascade" / "therapeutic window." |
| Hepatocyte-intrinsic (deconvolution) | **975 hepatocyte-intrinsic DEGs = 65.1% of 1,497 classifiable** (of TREAT-1918: 975 cell-autonomous + 522 composition-driven; 415 not-sig + 6 unmasked excluded) | `RNA-seq/results/celltype_attribution/c2_treat/`; NUMBERS Deconvolution | ⧗ | Supersedes 1,135/1,853=61.3% and 1,207/1,885=64.0%. This is the global Disease-vs-Control contrast; stage-progression signal is separately multi-cellular. |
| Fig 3D carrier shift (documented F-stage, Andrews n=58) | hepatocyte carrier share 14.5% (F0→F1) → 25.8% (F1→F2) → 20.2% (F2→F3) → 54.3% (F3→F4); F2→F3 dominated by non-parenchymal (Endo 28.1% / Chol 27.6% / Fib 24.1% / Hep 20.2%) | `routing_summary_by_transition.csv` | ⧗ | Restrict F-stage to 58 documented Andrews donors. Supersedes forbidden `F_stage_augmented` figures. |
| Fig 3G Hotspot disease-sig modules (donor-level) | **22 modules** (hep 7 / fib 6 / mac 6 / chol 3); T-cell arm does not survive donor-level | Fig 3G donor-level test; NUMBERS WS3 | ⧗ | **Supersedes run-level "54."** Multi-cellular hand-off preserved (β Pearson r 0.939 run-vs-donor). |
| Fig 3H stage CCC L-R heatmap (donor-level) | **DESCRIPTIVE only** — per-pair stage-progressive significance NOT claimed (0/13 pass family-FDR at donor resolution; was 9/13 at run level) | Fig 3H donor-level LMM; NUMBERS WS3 | ⧗ | **RETIRE "13 pairs pass all four gates."** |
| Single-cell atlas (2E) | 1,232,318 cells; **125 biological donors** (273 donor-mapped sequencing samples after run/CD45-fraction collapse); 7 integrated datasets | NUMBERS Fig 2e | ⧗ | **DRIFT LOCK: 125 donors** (donor census; earlier "121/269" was a curated-subset undercount). |
| Multi-cellular routing verdict | multi_cellular on all 7 transitions (0 hepatocyte-dominant); 3/3 stage-axis agreement | `stagedeg_celltype_verdict.csv` | ⧗ | Verdict SURVIVES donor collapse; per-pair/per-gene stage counts do not (see retired list). |

---

## Fig 4 (Independent tissue profiles of the evidence classes)

| Quantity | Canonical value (verbatim) | Source file/section | Spot-check | Notes / drift resolution |
|---|---|---|---|---|
| Evidence classes | Genetic-only / Disease-state-only / Convergent / Neither / Not-jointly-testable (defined from canonical DE + primary SuSiE only) | paper_outline §Evidence Classes | n/a | |
| Proteomics positive rates | neither 4.0% / genetic-only 7.0% / disease-state-only 20.0% / convergent 28.6%; adjusted vs neither OR: disease-state-only 6.90, convergent 10.69 | `evidence_class_validation_*.tsv`; paper_outline | ⧗ | Convergent vs disease-state-only NOT significant after MT. |
| Single-cell pseudobulk (donor-level, 2026-07-12) | **neither 1.625% (213/13,105) / genetic-only 2.670% (11/412) / disease-state-only 5.181% (63/1,216) / convergent 8.824% (3/34)**; adjusted vs neither OR: genetic-only 1.533, disease-state-only 3.967, convergent 5.755 | paper_outline §Single-cell; `docs/CODEBASE_CURRENT_STATE.md` | ⧗ | **DRIFT LOCK: donor-corrected values SUPERSEDE frozen r1 (3.5% / 3.9% / 8.8% / 11.8%).** Convergent-vs-genetic-only n.s. (OR 3.512, p 0.0825, q 0.2475). Release r2 pending. |
| Spatial SVG positive rates | neither 12.3% / genetic-only 17.1% / disease-state-only 35.1% / convergent 36.4% | `evidence_class_validation_*.tsv` | ⧗ | Disease-state-only strongly enriched; convergent NOT superior to disease-state-only. |
| Interpretation | Independent modalities corroborate the **disease-state axis**; convergent set NOT detectably superior to disease-state-only | paper_outline §Interpretation | n/a | Rejects "convergence is truth." |
| Fig 4A validation cascade (if retained) | universe 9,882 (transcriptomic 8,088 ∪ genetic 3,038; ∩ 1,244); proteo 714 / spatial 447 / scATAC 560; ≥1-assay 1,577, ≥2 = 141, all-3 = 3; proteomics enriched 1.2× P=1e-4, spatial n.s.; genetic-only subset enriched in NO assay | `scripts/figures/fig4a_overview_candidates.py`; NUMBERS Fig 4 | ⧗ | Supersedes 13,390-universe cascade. **Do NOT cite "3.76× proteomics∩spatial above-chance"** (13,390-universe, dead). |
| Cell-type COLOC enrichment (donor-level) | 11/12 lineages fold>1 (1.05–1.34×); **hepatocytes the SOLE non-enriched lineage (fold 0.858, n.s.)**; top = Circulating_NK_NKT 1.34× | `celltype_heritability_results.csv` (donor-level) | ⧗ | Qualitative claim (hep sole non-enriched) robust; exact hep fold snapshot-sensitive (~0.73–0.86). "NK 1.46×" is a myth. |
| Marquee COLOC hits | THRB PP.H4 = 0.9999 (GGT); RORA = 0.9977 (BBJ_GGT); HKDC1 = 0.9915 (AST); MTTP = 0.8353 (ALT, new PolyFun gain) | `gene_level_coloc.csv` | ⧗ | GCKR = **borderline 0.5453 susie / 0.4260 abf** (NOT a 0.95 hit); TM6SF2 = 0.6701 abf (susie fails). |

---

## Fig 5 + Discussion (Translational interpretation, model boundaries)

| Quantity | Canonical value (verbatim) | Source file/section | Spot-check | Notes / drift resolution |
|---|---|---|---|---|
| Strict drug benchmark | **22 clinical/approved targets vs 110 druggable controls**; convergence AUROC **0.526**, raw modality-count **0.582**, genetics-only **0.587**; only **1 approved target** in ranking universe | `strict_drug_benchmark.tsv`; `strict_drug_verdict.json`; claim `C_DRUG_RETIRED` = **failed_acceptance_gate** | ⧗ | FAIL gate. No ordered drug-development gradient. |
| Fig 5D held-out benchmark (5-panel) | Govaere expression 0.93 vs 0.95 (tie, p 0.62); **Moylan fibrosis 0.94 vs 0.77 (integration WINS, p 7e-8)**; protein biomarkers 0.87 vs 0.85 (tie, p 0.60); genetic loci 0.68 vs 0.85 (loses to matched channel, p 0.004); NIDDK clinical 0.69 vs 0.78 (n.s. worse, p 0.07) | `figures/main/fig5_convergence/panels/fig5d_delong_pairwise.csv` | ⧗ | Headline: integration MATCHES or EXCEEDS best single channel on independent expression/protein truth sets; bettered only by matched channel on same-modality sets. **Do NOT claim universal convergence superiority.** |
| gsMap exclusivity | **FAIL** — no direct MASLD/PDFF trait shows a significant convergent-vs-genetic-only advantage replicated in both cohorts | `gsmap_class_validation.tsv`; `gsmap_trait_replication.tsv`; claim `C_GSMAP_RETIRED` | ⧗ | Native gsMap retained as supplementary spatial annotation only. |
| Coding architecture | **FAIL** — credible-set coding estimable for **75 / 473 primary SuSiE genes (15.9% coverage)**, 207 reliable credible sets | `credible_set_coding_architecture.tsv`; claim `C_CODING_SCOPE` = insufficient | ⧗ | Insufficient for genome-wide architecture claim. **Do NOT say "~97% noncoding."** |
| GLP1R vs GCGR (5E, mechanism-only) | GLP1R hepatocytes **0.039%** (replicates <0.31% across 7 datasets); GCGR hepatocytes **17.9%** (6.4× next type) | `RNA-seq/results/glp1ra/replication/receptor_by_dataset.csv`; NUMBERS GLP-1RA | ⧗ | MECHANISM-ONLY. All reversal claims CUT (reach-test null p=1; selectivity null p=0.93). |
| Sex dimorphism | **8 sex-dimorphic genes** (1 Female-biased / 7 Male-biased / 0 Divergent); atlas maps 6 | NUMBERS Fig 5; LVQW-fixed canonical | ⧗ | **Female-enriched-progressor narrative RETIRED.** Supersedes dream's 290 and all v5/v6 (161 Suggestive) candidates. |
| Marquee framing | THRB = genetics-only clinically-validated mechanism; HKDC1 = high-ranked convergent hypothesis; RORA is NOT a top-ranked hero under the frozen score | paper_outline §Fig5 / §Discussion | n/a | Cite tier + PP.H4, NOT integer convergence rank (plateau-unstable). |
| Central thesis | Liver-trait regulatory genetics and the established MASLD transcriptome nominate **largely non-overlapping** genes — complementary, not interchangeable, target lists | paper_outline §Central Thesis | n/a | |

---

## RETIRED — must NOT appear anywhere

| Retired value / claim | Replace with |
|---|---|
| Named union **1,031** as headline | 1,030 named (1,031 = ensembl-level, one blank-symbol row) |
| **~97% (or ~96.7%) of inherited MASLD risk is noncoding** / 31,036-variant architecture count | 75/473 = 15.9% credible-set coverage; insufficient for architecture claim |
| **34/1,031 = 3.3% coding-led** as an architecture claim | 3.3% as descriptive locus annotation only (moves to supplement) |
| **gsMap OR 1.6–1.8** as headline / "inherited risk tissue-real only at convergence" / genetic-only depletion | gsMap exclusivity gate FAILED |
| **Convergence recovers both approved targets** / full drug-development gradient / "where colocalization alone did not" | strict benchmark: convergence AUROC 0.526, 1 approved target, gate FAILED |
| **Convergence is required / universally superior** ("convergence is truth") | complementary evidence classes; no universal convergence advantage |
| **Prospectively validated RORA/HKDC1** | external corroboration labeled by timing (not prospective) |
| Single-cell r1 rates **3.5% / 3.9% / 8.8% / 11.8%** | donor-corrected 1.625% / 2.670% / 5.181% / 8.824% |
| Fig 3G Hotspot **"54 modules"** (run-level) | 22 donor-level modules |
| Fig 3H **"13 L-R pairs pass all four gates"** (per-pair FDR) | 0/13 at donor resolution → descriptive only |
| LIANA **"4 → 10 → 16 stage-DEG L-R pairs"** | withdrawn (pseudoreplicated) |
| Coarse-stage per-gene **"19,043 significant"** | 1 donor-level survivor (multi-cellular verdict does not depend on it) |
| Deconfounding **95.5%** collapse | 90.2% (matched 615; envelope 87.3–96.1%) |
| Fig 4A **"3.76× proteomics∩spatial above-chance"** (13,390-universe) | 9,882-universe cascade; proteomics 1.2× P=1e-4 |
| Atlas **411 / 412 / 197 / 156 / 84 / 218 / 444** columns | 443 |
| DEG-COLOC categorical **Fisher OR 1.67 / p 5.5e-6 / 196-5,484 = 3.6%** | overlap NOT enriched (SuSiE OR 0.89 p 0.60) |
| Sex **290 dimorphic / female-enriched-progressor** | 8 sex-dimorphic genes |
| Convergence Tier-1 **1,018 / 1,049 / 513** | 677 (full-portfolio 1,018 = supp/sensitivity only) |
| Cross-ancestry **67** validated / union-based **168** non-EUR | 62 (pre-MVP, superseded) → main Fig 2G non-EUR-unique 36, gated 14 GWS |
| "10,135 apparent NASH DEGs" as the deconfounding numerator | C11 is a different contrast; use C2→C13 90.2% |

---

## Banned wording (do not use in prose)

- **"TREAT"** — say "effect-size-aware interval-null FDR gate" or "TREAT-style FDR gate ≡ report as `adj P < 0.05` at lfc = 0.25"; call the set "1,918 canonical DEGs" (1,915 unique symbols).
- **"novel" / "first" / "unprecedented" / "prospective"** as priority claims — the last only if a public/immutable pre-registered prediction predates the validation.
- **"multicellular cascade" / "therapeutic window"** for cross-sectional stage data — use "stage-ordered remodeling," "adjacent-stage contrast."
- **"validated target"** for observational evidence; **"causal gene"** from COLOC alone — use "colocalized gene," "candidate effector."
- **"acquired risk" / "tissue-real"** — banned per claim ledger.
- **"full drug-development gradient" / "recovers both approved targets" / "convergence is required"** — retired failed-gate claims.
- Frame overlap as **"largely non-overlapping,"** NOT "statistically/biologically independent."

**Preferred terminology:** colocalized gene / candidate effector · liver-trait regulatory genetics · disease-state-associated · stage-ordered remodeling · independent corroboration · evidence-weighted heuristic (the convergence score is NOT a posterior/probability) · source-level evidence vector · primary and sensitivity definitions.

---

*Spot-checks reproduced this pass (✅): 1,918/1,915/1,419/499/27,638 DEGs; SuSiE 473/369/267; ABF 863/375/234; union 1,031 ensembl → 1,030 named; orthogonality 473/447/34/92.4% + union 1,030/971/74/92.4% + direct 172/160/11/93.1% + Fisher ORs 0.89; convergence Tier-1 677 (204/25,668/638, Σ27,187); atlas 27,187 × 443. All matched their ledger values exactly.*
