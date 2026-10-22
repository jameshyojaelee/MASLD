> **Historical competitor/source dossier only (2026-08-10).** Source-paper summaries may remain useful, but all “our work,” winner, figure, journal, target, and next-analysis sections predate the standalone Resource scope and may be stale. Reverify source facts, and use [`../PAPER.md`](../PAPER.md), [`../STATUS.md`](../STATUS.md), and [`../RESULTS.md`](../RESULTS.md) for current project claims.

# Competitor Briefing — Elison/Chang/Xie et al. (Ren + Gaulton + Kisseleva + Brenner)

**Paper**: *Single cell multiomics reveals drivers of metabolic dysfunction-associated steatohepatitis*
**Posted**: medRxiv 2025-05-11 · DOI 10.1101/2025.05.09.25327043
**Co-firsts**: Weston Elison, Lei Chang, Yang Xie · **Corresponding**: Allen Wang, Tatiana Kisseleva, David Brenner, **Bing Ren**, Kyle Gaulton
**NYGC angle**: Bing Ren is now NYGC affiliate (#14). This is his entry-flag paper at NYGC.
**Code/data**: Code at `github.com/Gaulton-Lab/FNIH.Liver`. Zenodo 15298484 (available). GEO embargoed until publication. Processed portal: `epigenome.wustl.edu/MASLD/`.

---

## 1. What they did, in one paragraph

5-modality single-cell multiomics on **86 deceased-donor livers** from a single OPO (Lifesharing San Diego), pooled 15–30 per pool with genotype-based demultiplexing:
- **10x Multiome** (snRNA + snATAC) — 408,178 nuclei across 86 donors → 611,576 cCREs
- **Droplet Paired-Tag H3K27ac** — 94,125 nuclei (Ren-lab protocol, ab4729)
- **Droplet Paired-Tag H3K27me3** — 86,554 nuclei (ab192985)
- **Droplet single-cell Hi-C** (Chang & Ren droplet Hi-C, NOT sciHi-C) — 89,891 nuclei, 85 donors, scHiCluster loops at 10 kb
- **Visium HD** (FFPE, whole-transcriptome probe v2) — 7 donors (2 N / 2 MASL / 2 MASH / 1 MetALD), Bin2cell + StarDist → **1.73M segmented cells**

Disease groups (NASH-CRN, double-blinded pathology): Normal (26) / MASL (25) / MASH (23) / MetALD (12); fibrosis collapsed into Fib- (F0–F2, n=19) vs Fib+ (F3–F4, n=16) — **no F-stage continuum**.

Integration is **Seurat anchor label-transfer, not joint embedding**. Paired-Tag and scHi-C cells are NOT the same physical nuclei as multiome — modalities are glued together by donor + cell-type label via prediction-score-0.8 anchor mapping. *There is no scGLUE/MultiVI joint VAE.*

---

## 2. Their headline claims (with strength scoring)

| # | Claim | Strength | Caveat |
|---|---|---|---|
| 1 | 5-modality SC atlas of MASLD (incl. scHi-C + H3K27me3) | **★★★ unprecedented** | Crown jewel — we cannot match |
| 2 | 43.5–47.9% of cCREs are novel vs prior catalogs | ★★ partly methodological | Comparator is ENCODE3 + Cherry/Wang only; not ENCODE4 rDHS or Catlas. ~40% are retrotransposon-overlapping (filtered out by Cherry/Wang). |
| 3 | 88% of cALT-MASLD loci annotated (68/77); 1/3 in hepatocyte cCREs | ★★ | Rests on a **single GWAS** (Vujkovic 2022 cALT — a damage proxy, not histological MASH) |
| 4 | Stepwise hepatocyte GRN cascade: MASL (HNF4A/RORA/CUX2) → Fib- MASH (FOXO1/NR4A1) → Fib+ MASH (TEAD1/KLF6/CREB5/NFKB/RELB/ATF3) | **★★★ clean, mechanistic** | Implicitly a "switch" story — functionally equivalent to what we *retired* on 2026-05-09 |
| 5 | "Fib-a" / "Fib-adj" hepatocyte subpopulation with H3K27me3 rewiring + spatial co-localization at fibrotic foci + HSC→Fib-adj COL1A1–CD44 signaling | **★★★ headline** | **n=1 spatial donor** (HL2201019) for the Squidpy neighborhood enrichment |
| 6 | 10,591 caQTLs / 1,644 eQTLs / 418 H3K27ac-QTLs / 391 H3K27me3-QTLs at FDR<0.10; 18 loci colocalize with cALT GWAS PP4>0.8 | ★ at scale, **★★★ in concept** | **N=86 is severely underpowered for QTL** (published caQTL studies need 100–200+ donors). No replication, no permutation null reported. |
| 7 | ChromBPNet sequence model nominates causal variants at 24/26 hepatocyte-caQTL loci | ★★ | Methodologically competent, not novel |
| 8 | Therapeutic implications | ★ | **No LINCS / CGP / drug repurposing.** Names TFs + CD44; no compound, no clinical-trial cross-walk, no plasma biomarker. |

**No dedicated Limitations subsection** — three Discussion paragraphs are summary only. Reviewer bait.

---

## 3. Head-to-head matrix

| Axis | Elison/Ren | Us | Winner |
|---|---|---|---|
| Donor count (single-cell) | 86 | ~260 across 7 GEO datasets (F-stage augmented to 127) | **Us** (donor diversity) |
| Donor count (bulk) | 0 | **847 in 5-cohort dream mega + 4 control-less** (1,277 collected) | **Us** by 10× |
| Single-cell modalities | RNA, ATAC, H3K27ac, H3K27me3, Hi-C, Visium HD = 6 | RNA + 1 multiome (Wang) + bulk ATAC + Visium = 4 (bulk + spot scale) | **Them** decisively |
| Hepatocyte N | ~250K (subset of 408K) | **657,804** (Pipeline 10) / 522K (v2 atlas) | **Us** |
| Hepatocyte subtypes | 5 (Z1/Z2/Z3/MASH-a/Fib-a) | **43 fine → 5 meta** (Progressor/Moderate/Stable/Neutral/Healthy) | Comparable; ours finer |
| F-stage modeling | Fib-/Fib+ binary (F0–2 vs F3–4) | **Continuous F_stage_augmented** scVI ordinal, LOOCV QWK 0.862 | **Us** |
| GWAS portfolio | 1 (Vujkovic cALT) | **28** (17 EUR + 5 EAS + 3 AFR + 3 SAS) | **Us** by 5× |
| Coloc PP.H4>0.8 | 18 loci | **280 genes** (SuSiE) + 255 (ABF) | **Us** by 15× |
| Cell-type caQTL/hQTL | **10,591 / 809** at N=86 | none | **Them** (but underpowered) |
| scHi-C 3D loops | 4 cell types | none | **Them** |
| ChromBPNet sequence model | yes (24/26 caQTL loci) | none | **Them** |
| TWAS (S-PrediXcan / OTTERS) | none | **18,975 Broadaway eGenes tested** | **Us** |
| INTACT / HyPrColoc / sc-TWAS | none | yes (full stack) | **Us** |
| Multi-ancestry replication | trans-ancestry GWAS but single | **67 EUR+EAS cross-validated genes** | **Us** |
| Cross-species (mouse) | none | **9 mouse datasets × 5 diets (463 samples)**; FPC 48/48 TFs | **Us** entirely |
| Plasma proteomics chain | none | **Olink 1,461 × 177**; **87 triple-concordant ligands** (TIMP1, VCAN, SERPINE1, CCL3/4/5/14, GRN, NAMPT) | **Us** entirely |
| Drug repurposing | none (named TFs only) | **LINCS 1,107 reversal / 625 network / 20 multi-layer / 4 Strong** (resmetirom-THRB, OCA-FXR) | **Us** entirely |
| Cell-cell communication | CellChat only, Fib+ donors only, 1 axis (COL1A1–CD44) | **LIANA + CellChat + NicheNet + spatial-CCC multimethod + COSMOS** | **Us** by depth |
| Spatial resolution | Visium HD 2 µm bin (7 donors) | Visium 55 µm (GSE192741 + Vu, 2 cohorts) | **Them** on resolution, **us** on replication |
| Pseudotime trajectories | Monocle3 only, Zone3→Fib-a only | **Palantir + DPT + CytoTRACE + PAGA** across 5 cell types; macrophage Δ=0.73 | **Us** by breadth |
| Sex stratification | not reported | **3,080 sex-dimorphic genes** (interaction model); Progressor female-enriched | **Us** entirely |
| Multi-cellular programs | none | **DIALOGUE + cNMF k=16 + Hotspot 7-cell-type** | **Us** entirely |

---

## 4. Where they OBJECTIVELY beat us (be honest)

1. **scH3K27ac + scH3K27me3 + scHi-C in disease tissue** — the only such liver dataset that exists.
2. **Per-cell-type cCRE × ChromHMM lattice** — 47.7% of hepatocyte genome annotated to chromatin states vs 5.9% from cCREs alone.
3. **Cell-type caQTL/hQTL in the same disease donors** — no GTEx tissue mismatch; ETS/FOS in HSC, SPIC/IRF/STAT/CEBP in myeloid.
4. **ChromBPNet sequence models** for allele-resolved cCRE prediction.
5. **Visium HD 2 µm bin resolution** with 1.7M cells (beats our Visium 55 µm).
6. **Single prospective single-site cohort** — no cross-cohort batch concerns reviewers will hit us with.
7. **4-corresponding-author firepower** (Ren + Gaulton + Kisseleva + Brenner) — institutional weight.
8. **Catchy nomenclature** ("Fib-a" / "Fib-adj") sells better than "Progressor/Moderate/Stable."

## 5. Where we OBJECTIVELY beat them

1. **Cohort breadth**: 9 prospectively meta-analyzed cohorts × 847 mega samples × 5 control-bearing studies; LOO-CV recovery 81.6%. Replication is built-in.
2. **Genetic coverage**: 28 GWAS × 4 ancestries × SuSiE-COLOC (364 PP4>0.5) + TWAS + INTACT + HyPrColoc + sc-TWAS + monogenic Okur 2026 rare-variant layer. They have 1 GWAS, ABF only, no TWAS.
3. **Cross-species (mouse)**: 9 datasets × 5 diets × FPC TF concordance. They are human-only.
4. **Translational arc**: LINCS 1,107 reversals + 4 Strong-validated drugs (THRB/resmetirom) + Olink plasma chain (32 Progressor secreted markers) + COSMOS mechanistic ligand→TF→target reconstruction. They have TF names + CD44.
5. **F-stage continuous modeling** with scVI ordinal augmentation (QWK 0.862). They collapse to Fib-/Fib+ binary.
6. **Multi-step cascade narrative** (4 CRN transitions, Bayesian changepoint rejected the switch). Their stepwise MASL → Fib- → Fib+ trichotomy is functionally a switch claim they did not validate against monotonic alternatives.
7. **Sex stratification**: 3,080 sex-dimorphic genes via interaction model; Progressor female-enriched. They do not stratify by sex anywhere.
8. **Multi-cellular programs**: DIALOGUE + cNMF + Hotspot autocorrelation — they have none.
9. **Cell-type GWAS** (328 INTACT genes) + GWAS-ATAC motif disruption (7,766 CS variants → 380 motif-disrupted → 12 disease-regulon TFs: HNF4A:12, RORA:11, **THRB:8**).

---

## 6. Direct biological overlap to verify in person

| Their cell/program | Our equivalent | Question to confirm |
|---|---|---|
| **Fib-a hepatocytes** (NFKB1/STAT4/HIF1A/TEAD1/SOD2) | **Progressor meta-subtype** (COL15A1/CDH11/LTBP2/SERPINE1) | "Are Fib-a and Progressor the same population? Your markers are TF regulators; ours are ECM effectors." |
| **Zone 3 → Fib-a pseudotime** (Monocle3, ρ=-0.47) | We have not run within-hepatocyte sub-cluster pseudotime | "What's the lineage evidence beyond inverse-abundance correlation?" |
| **MASL GRN: HNF4A/RORA/CUX2/NR1I2** | Our **P5 Hepatic-metabolic** NMF program | Likely 1:1 mapping. HNF4A/RORA already in our 12-disease-regulon panel. |
| **Fib+ GRN: TEAD1/KLF6/ETS1/CREB5/NFKB1/RELB** | Our **P1 Pro-inflammatory + P6 Stellate-myofibroblast** | TEAD1 = canonical YAP/TAZ. Already in our 24-TF GWAS-ATAC panel. |
| **THRB/ONECUT1 decreased across all stages** | **THRB Strong drug target** (GGT SuSiE PP4=1.000, resmetirom benchmark) | Direct concordance! Their finding validates ours. |
| **HSC → Fib-adj COL1A1–CD44** | Our 87 triple-concordant LIANA ligands incl. TIMP1, VCAN; spatial Progressor↔TIMP1 ρ=0.196 | Different ligand pairs; possibly complementary. |
| **PPP1R3B / KRT8 / EFHD1 / CEBPA distal-linked loci** | Not specifically featured in our top-20 multi-evidence | Why these specific loci over PNPLA3 / TM6SF2 / HSD17B13 / MBOAT7? |
| **chromHMM 5-state** | We have ATAC-only chromatin states | They can call repressive states; we cannot. |

---

## 7. The 10 questions to ask Weston / Lei / Yang

### Strategic / collaboration (open with these)

1. **NYGC adjacency.** "Bing just moved to NYGC. Is there appetite for a coordinated submission window or a companion-paper structure? Our paper occupies the cross-cohort + cross-species + multi-ancestry + translational niche — orthogonal to your atlas-of-chromatin framing."
2. **Data availability.** "When does the GEO accession come off embargo? Is dbGaP individual-level genotype access available before publication, and would Sanjana lab be on the eligible-collaborator list for replication?"

### Scientific (probe their weakest spots)

3. **caQTL power at N=86.** "10,591 caQTLs at FDR<0.10 from 86 donors with 15 expression PCs — what's the lambda inflation, and did you do permutation null per cell type? Published caQTL studies usually need 200+ donors. How are you planning replication of the 18 colocalized loci?"
4. **Single GWAS.** "Why only Vujkovic 2022 cALT? cALT is a damage proxy — did you test FinnGen NAFLD/NASH/HCC, Ghodsian, Anstee, or Pan-UKBB AST/PDFF? How robust is your 88%-of-loci annotation rate to GWAS choice? Our 28-GWAS portfolio shows substantial phenotype-specific signal divergence."
5. **F-stage continuum.** "You collapse F0–F2 vs F3–F4. Have you tested ordinal modeling, or specifically F1→F2 as a transition? Our analysis finds F1→F2 is the *smallest* DEG-count transition (2,600 DEGs) and Bayesian changepoint rejects discrete anchors — your Fib-/Fib+ binary cannot distinguish."
6. **n=1 for Fib-adj neighborhood.** "Your Squidpy nhood_enrichment for Fib-adj is one donor (HL2201019). Have you replicated in your other MASH/MetALD spatial donor, or is replication deferred?"
7. **Sex stratification.** "Of 86 donors, what's the F:M ratio per stage? We see Progressor hepatocyte expansion as female-enriched ~3× (interaction model on 153 control donors). Have you tested sex × stage interaction for Fib-a abundance?"

### Defense (where they will hit us in review)

8. **Cross-cohort batch effects.** "How would you assess whether a cross-cohort DEG from a 9-study meta-analysis is real? We have LOCO recovery 81.6%, ComBat-sensitivity, dream mixed-models — what would convince you?"
9. **Chromatin gap.** "We don't have single-cell chromatin. We do have GWAS-ATAC motifbreakR over 7,766 fine-mapped credible-set variants, disrupting HNF4A 12×, RORA 11×, THRB 8×. Would direct comparison to your hQTL hits at those loci be a useful triangulation?"

### Bonus

10. **PNPLA3 / TM6SF2 / HSD17B13 / MBOAT7.** "Which canonical lipid-genetics loci did you recover at the caQTL/hQTL level, and at which cell type? The PPP1R3B / KRT8 / CEBPA examples in Fig 4 are unusual featured choices."

---

## 8. Strategic positioning for our paper

### Reposition explicitly
**Not** a single-cell multiomics atlas. **Yes** a causal + translational + cross-cohort + cross-species atlas. Their Nature submission is orthogonal to our Cell Metabolism Resource target.

### Differentiation pillars (lead with these in cover letter)
1. **9 cohorts × 847 samples** — cross-cohort replication is built-in (LOCO ρ=0.67)
2. **9 mouse × 5 diets** — cross-species (they have zero)
3. **28 GWAS × 4 ancestries** with SuSiE-COLOC (they have 1 GWAS)
4. **Drug repurposing at scale**: LINCS + CGP + 4 Strong-validated (resmetirom/THRB benchmark)
5. **Olink plasma chain**: 32 Progressor secreted markers, mid-stage inflection detectable in blood
6. **F-stage continuous** scVI-augmented (QWK 0.862) vs their Fib-/Fib+ binary
7. **Multi-step cascade** with M6 Bayesian changepoint validation vs their trichotomy

### Defensive add (recommended)
**1-figure ENCODE/ChIP-Atlas pseudo-chromatin supplement** — overlay our DEGs / COLOC hits against ENCODE liver H3K27ac/H3K27me3 + ChIP-Atlas TF tracks. Frame as cross-validation, not novelty. Preempts the "no chromatin" reviewer ding. <2 weeks; does not require their data.

### Journal trajectory (UNKNOWN — inferred only)
- **No direct evidence.** Paper has no journal disclosure. Inference is from scope (5 modalities + 4 corresponding authors + FNIH consortium funding → flagship-tier ambition) and Ren's track record (Cell 2021 chromatin atlas).
- **medRxiv (not bioRxiv) posting** is a real signal — argues toward translational/clinical journals: Nature Medicine, Cell Metabolism, JCI, Hepatology. Pure-Nature basic-science atlases usually go to bioRxiv.
- **Plausible targets, recalibrated**: Nature Genetics (GWAS framing fits) > Nature Medicine (medRxiv signal) > Cell (Ren track record) ≥ Nature ≥ Cell Genomics / Nat Comms (fallback).
- **Implications for our positioning**:
  - If Nat Genet → competing on causal/GWAS axis. Our 28-GWAS multi-ancestry portfolio is the rival.
  - If Nat Med → competing on translational axis. Our Olink + LINCS + 4 Strong drugs is the rival. **This is the worst case for us** because lanes overlap.
  - If Cell or Nature → competing on prestige. Lanes still mostly orthogonal.
- **Ours**: Cell Metabolism Resource (primary) → Nat Metab Resource (secondary). Mostly orthogonal unless they go Nat Med.
- **Ask Weston directly.** Open with "How's the submission process going?" / "Where are you sending it?" — natural conversation, real intel.

### Collaboration calculus
- **Offer**: our bulk-RNA-seq mega-analysis as their external replication cohort (low-cost goodwill; they get reviewer-shield, we get acknowledgement + relationship)
- **Share**: cohort dimensions, GWAS portfolio, modality list, journal target (signals non-competing)
- **Withhold**: specific drug candidates, sex-stratified gene lists, multi-step cascade transition specifics, NMF program biology — until you check with Neville
- **Frame**: "Our papers are orthogonal — atlas-of-chromatin vs causal-cross-cohort-translation. Companion-paper structure if timing aligns."
- **Test**: do not commit to data sharing in the first meeting. Say "let me check with my PI."

---

## 9. Key facts for the meeting (cheat sheet)

- **They have**: 86 donors, 5 single-cell modalities including H3K27me3 + scHi-C, 10,591 caQTLs, 47.7% chromHMM coverage of hepatocyte genome, Visium HD on 7 donors, 5-cluster hepatocyte map with Fib-a / Fib-adj as MASH-enriched.
- **They lack**: replication cohort, mouse, plasma, LINCS, sex stratification, multi-ancestry GWAS, TWAS, INTACT, F-stage continuous, dedicated Limitations.
- **We have**: 847 bulk + 260 sc donors across 9 cohorts, 28 GWAS × 4 ancestries × SuSiE-COLOC (364 hits), 4 Strong drugs (resmetirom-THRB), Olink chain (87 ligands), 5 mouse diets, F_stage_augmented QWK 0.862, 5 hepatocyte meta-subtypes incl. Progressor (1.3→63.1%), DIALOGUE + cNMF + Hotspot.
- **We lack**: single-cell chromatin in disease, scHi-C, ChromBPNet, Visium HD.
- **Direct overlap**: THRB-decreasing across stages (they validate our drug target); HNF4A/RORA in their MASL GRNs match our GWAS-ATAC disease regulons; TEAD1 in their Fib+ GRN matches our P6.

---

# ADDENDUM — Round 2 Deep-Cut Findings

Five additional agents probed angles not in the first pass. Highest-value findings, in priority order:

## A. The canonical-loci hole (BIGGEST FINDING)

**`grep` for PNPLA3 | TM6SF2 | HSD17B13 | MBOAT7 | GCKR | APOE | PCSK7 | TRIB1 | MTARC1 | TMC4 in their full text = 0 hits.** Not in main text, not in Methods, not in any reference. Unprecedented for a MASLD genetics paper of this scope.

Their featured loci are **PPP1R3B (8p23), KRT8/18 (12q13), CEBPA (19q13), XBP1 (22q12), EFHD1 (2q37), NFATC4, SCAI, RORC** — all caQTL-recoverable. The reason is structural: their pipeline is caQTL-first, and **PNPLA3 I148M (protein aggregation), HSD17B13 (splice variant), and MBOAT7 (largely sQTL/protein-level) are NOT caQTL-mediated**. Their cCRE/caQTL framework is **structurally blind to the two most famous MASLD risk variants in human genetics**.

What our atlas recovers at the same loci:

| Locus | Their paper | Our atlas (SuSiE PP.H4 / TWAS / dream) |
|---|---|---|
| PPP1R3B | featured Fig 4G | **PP4=0.958** PolyFun; dream LFC +0.39 padj 1.8e-6 — **concordant** |
| KRT8 | featured Fig 4F | **PP4=0.869**; dream LFC +0.75 padj 1.5e-36 — **concordant** |
| EFHD1 | featured | **PP4=0.984 ABF / 0.841 SuSiE**; HyPrColoc 0.976; dream LFC -1.01 padj 2.6e-13 — **concordant** |
| **PNPLA3** | **absent** | PP4=0.008 (correctly null — post-translational); sceQTL PP4=0.034 |
| **TM6SF2** | **absent** | PP4=0.520 SuSiE / 0.670 ABF; HyPrColoc 4.4e-5 — borderline, recoverable |
| **HSD17B13** | **absent** | PP4=0.049 (correctly null — splice variant; needs sQTL framework) |
| **MBOAT7** | **absent** | **TWAS z=−3.69, p=2.2e-4** (strong); PP4=0.069 — recovered by TWAS, not COLOC |
| **GCKR** | **absent** | **PP4=0.950 PolyFun SuSiE** — strong direct hit, missed entirely |
| CEBPA | featured | PP4=0.022 (we do NOT recover) |
| XBP1 | featured | PP4=0.063 (we do NOT recover) |
| APOE | absent | PP4=0.313; dream LFC +0.19 padj 5e-3 |
| TRIB1 | absent | PP4=0.334; dream LFC -0.28 padj 1e-7 |

**Differentiation pillar (cover-letter language)**: *"Our 28-GWAS multi-ancestry SuSiE-COLOC + TWAS + INTACT framework recovers the canonical MASLD loci (GCKR PP4=0.950; TM6SF2 borderline; PNPLA3/HSD17B13 correctly null at the COLOC level and handled via TWAS / monogenic / rare-variant layers), in addition to the novel caQTL hits (PPP1R3B, KRT8, EFHD1) shared with the Elison/Ren multiomic atlas."*

**Three sharp questions about loci** (probe gently — they may be defensive):
1. "Why PPP1R3B/KRT8/CEBPA/XBP1/EFHD1 in Fig 4 over PNPLA3/TM6SF2/HSD17B13/MBOAT7? Did the canonical loci fail to colocalize, or did you choose to emphasize the novel ones?"
2. "Does your caQTL framework miss the canonical variants by construction? PNPLA3 is post-translational, HSD17B13 is splice — neither would produce a chromatin signal."
3. "Of the 18 PP4>0.8 colocalizations, how many are already in the GTEx Liver eQTL catalog vs genuinely novel?"

---

## B. Statistical rigor — three concrete weaknesses

(1) **Fibrosis-as-covariate paradox** (Methods L2050–51). They regress out `fibrosis_stage + disease_label` before QTL discovery, then colocalize the residual signal against a cALT GWAS that is itself fibrosis-correlated. They cannot disentangle whether their 1,644 eQTLs are healthy-tissue cis-effects or disease-modifying. **They never stratify QTL discovery by Fib-/Fib+** (N too small to admit; the limitation is unstated).

(2) **Pool × disease confounding unmodeled**. 86 donors across 3–6 pools (15–30 each, L1518). Pool is a covariate, not an interaction. Paired-Tag antibody-lot drift (ab4729, ab192985) across pools is the natural worst case. They do not test pool × disease, pool × fibrosis, or pool × cell type.

(3) **No global FDR**. BH at FDR<0.10 within each axis independently. With 8 cell types × 4 modalities × 7 contrasts ≈ 224 parallel FDR universes. Family-wise expected false discoveries are substantial. No hierarchical or stratified FDR scheme.

Other rigor gaps: no held-out donors for QTL/GRN; no external bulk replication (no Govaere, no Liver_Atlas); no bootstrap CIs anywhere; chromBPNet train/validate is **chromosome-not-donor** split (allele leakage across train/validate possible at PPP1R3B 8p23 / KRT8 12q13 — the two loci they feature); SCENIC+ GRN run on one 50K subsample with no donor-LOO Jaccard.

**Reviewer baits we should preempt** (for our paper):
- "Cross-cohort batch effects" → defense: dream mixed-model `(1|cohort)`, ComBat-seq sensitivity ρ>0.99, LOCO 81.6% (disclose Chen-held-out floor 61% prominently)
- "Why dream over DESeq2 + two-tier LFC>0.5" → defense: dream random-effect absorbs cohort dispersion; LFC>0.5 chosen by LOO-CV plateau (NUMBERS.md)
- "Sex-augmentation circularity" → defense: full-LOOCV QWK 0.86 from Andrews-only 0.74 lift is out-of-sample, not anchor-confined

---

## C. Citation snubs — who they didn't cite

| Paper | Their citation | Signal |
|---|---|---|
| **Govaere 2020 Sci Transl Med + 2023 Nat Metab** | **ZERO** | Most-cited human MASLD bulk paper. **Uncontested territory for us.** |
| **Anstee 2020 Nat Genet GWAS** | absent | Used Vujkovic cALT proxy instead of MASH-histology GWAS |
| **Ghodsian 2021 / FinnGen R12** | absent | 1-GWAS choice is deliberate pruning |
| **Guilliams 2022 Cell Liver_Atlas** | absent | Gold-standard healthy liver scRNA reference |
| **Wang 2024 GSE244832 / Su 2024 Hepatology** | absent | Direct multiomic competitors at scale |
| Andrews 2022 | shallow (cholangiocyte only) | not framed as scRNA comparator |
| Ramachandran 2019 (SAM macrophages) | medium | **Nomenclature ignored** — they invented "Fib-a/Fib-adj" instead of using SAM/SAH |
| **Hoshida / Friedman / Sanyal / Loomba** | **all absent** | Notable UCSD lane separation — Loomba is *at UCSD* but not coauthor |

**Implication for our paper**: lead the cover letter with **Govaere 2020/2023 as our primary bulk anchor** — entirely uncontested. Cite Anstee 2020 + Ghodsian 2021 + Schipper 2025 prominently in our Fig 4. Cite Guilliams + Wang + Su as our scRNA reference panel. Position our **Progressor** explicitly relative to **Ramachandran SAM nomenclature** — molecular ECM-effector state vs their TF-regulator state.

**FNIH consortium signaling**: Acknowledgements L1201 confirms FNIH funding to B.R., K.J.G., D.A.B., T.K. The 4-corresponding-author roster IS the consortium PI cluster. Likely AMP-MASH or NIDDK MASH initiative. **Pfizer interlocking**: K.J.G. honoraria, R.M. Evans (Salk, listed) is Pfizer employee+shareholder. Reviewer-pool likely includes Ramachandran/Henderson (Edinburgh), Watson/Mullen, Aerts/Bravo, Kundaje (ChromBPNet). Our reviewer pool likely includes Govaere/Daly/Anstee (EPoS) — non-overlapping.

---

## D. Figure mechanics — what to steal, what to avoid

| Steal | Source | Adapt for |
|---|---|---|
| Browser-track + Hi-C arc composite (4 stacked tracks: ATAC + H3K27ac + H3K27me3 + chromHMM) | Fig 2C (CRP/CD163), Fig 4F/G (KRT8, PPP1R3B) | Our **THRB locus** hero panel with ENCODE bulk-liver H3K27ac/me3 + GWAS-ATAC variant pin + PP4=1.000 |
| GRN cluster heatmap × stage columns | Fig 3G | Our **24-TF GWAS-ATAC panel × 4 F-transitions** |
| Pseudotime-on-tissue ridge | Fig 5J bottom | Our Progressor along Palantir trajectory |
| Squidpy `nhood_enrichment` z-score matrix | Fig 5L | Our spatial CCC niche (cleaner than chord diagrams) |
| Cell-type × disease × modality matrix | Fig 3A | Our cell-type × F-transition × modality (RNA/Olink/SVG/GWAS-ATAC) |

| Avoid | Why |
|---|---|
| Featuring obscure loci (PPP1R3B/KRT8/CEBPA) | We MUST feature PNPLA3/TM6SF2/HSD17B13/MBOAT7 — anchor reviewer biology intuition first |
| Stage trichotomy MASL/Fib-/Fib+ | Our multi-step cascade explicitly disagrees |
| Single-method pseudotime (Monocle3) | We have Palantir+DPT+CytoTRACE+PAGA — lead with multi-method |

Their **5 most vulnerable panels** (don't attack publicly, but know the soft underbelly):
1. Fig 5K Visium HD Fib-adj — **n=1 donor (HL2201019)**
2. Fig 4A,B caQTL Manhattan — 10,591 hits at N=86, no permutation null
3. Fig 2E + Fig 4D–E enrichment — entire genetics arc on single Vujkovic cALT
4. Fig 4F–G featured loci — PPP1R3B/KRT8 instead of canonical lipid quartet
5. Fig 5A,B,J Fib-a + pseudotime — Monocle3 only; ρ=-0.47 inverse-abundance does the lineage heavy-lifting

---

## E. Data sharing posture

| Resource | Status | What's there |
|---|---|---|
| **Zenodo 15298484** | **LIVE NOW, CC-BY-NC 4.0** | Suppl Tables 1–18: donor metadata, per-cell-type DEGs, cCRE coords, SCENIC+ GRN edges, caQTL/hQTL/eQTL summary stats, loops, niches |
| **WashU Epigenome Browser** (epigenome.wustl.edu/MASLD/) | **LIVE NOW** | bigWig tracks per cell type × modality + bedpe loops + cCRE BED. No DAC restriction. Downloadable. |
| **GitHub Gaulton-Lab/FNIH.Liver** | **LIVE NOW** | Per-figure code + SCENIC+/chromBPNet wrappers. Paths hard-coded — mineable but not rerunnable without their VCF. |
| **GEO raw** | Embargoed "upon publication" → late 2025 / mid 2026 | Not usable |
| **dbGaP genotype** | **NOT MENTIONED anywhere** | They have WGS/imputed VCF for demuxlet + QTL but disclose no accession. Real vulnerability for them; reviewer-1 ammunition. |

**License**: CC-BY-NC 4.0 — academic atlas reuse compliant; commercial drug-discovery use forbidden.

**FNIH DAC** is not disclosed. **This is the single most valuable question to ask Weston in person**: which FNIH DAC oversees this consortium, and whether NYGC can register as a replication site for early-access genotype + matched bulk.

---

## F. The 7 concrete collaboration vehicles (ranked by leverage)

| Vehicle | Risk to us | Leverage | What we offer | What we get |
|---|---|---|---|---|
| **GWAS portfolio at canonical loci** | LOW | HIGH | SuSiE-COLOC PP.H4 across 28 GWAS at PNPLA3/TM6SF2/HSD17B13/MBOAT7/GCKR | Their caQTL/hQTL beta+SE at matched loci → joint table |
| **Mouse cross-species offering** | LOW | HIGH | Project their Fib-a/Fib-adj signature on 9 mouse × 5 diets | Acknowledgement + relationship; cross-species validation they lack |
| **HCC overlay** | LOW | HIGH | FinnGen/deCODE/Ghouse HCC COLOC vs their Fib-adj cCRE catalog | Their Fib-adj→HCC framing becomes empirical, not rhetorical |
| **Bulk replication exchange** | LOW (public data) | MED | Our 5-cohort dream DEG (already public) | Their per-cell-type pseudobulk DEG tables |
| **TF GRN drug-reversal benchmark** | MED | HIGH | LINCS reversal on their HNF4A/RORA/CUX2/TEAD1/KLF6/NFKB | Bridges their TF-only story to translation; companion-paper framing |
| **Olink projection** | MED (must embargo our top-5) | HIGH | Project their Fib-a list onto Olink 1,461×177 | Plasma supplement for their paper |
| **NYGC seminar** | NONE | LOW–MED | Sanjana seminar invitation | Academic warming; activates Bing's NYGC affiliation |

**Avoid** (until journal target is known): companion-paper structure if they're heading to **Nat Medicine** (lane overlap with our Cell Metab Resource target).

## G. HOLD BACK — do NOT share in first meeting

1. Specific drug candidates beyond resmetirom/THRB + OCA/FXR canon (the 4-Strong list)
2. Sex-stratified gene lists (sex_v5 interaction rebuild: 29 F-only, 4 M-only, 24+94 divergent)
3. Multi-step cascade verdicts M1–M8 (our explicit F2-switch rejection — counter to their MASL→Fib-→Fib+ trichotomy)
4. Multi-output predictor numbers (S2 AUROC 0.939 all-features → 0.853 molecular-only)
5. v2 chain Phase 0.5 retrain (F-stage LOOCV QWK 0.862)
6. Okur 2026 monogenic 30/90 panel + 31 P/LP variants (Zenodo embargoed)
7. 87 triple-concordant Olink ligands; 32 Progressor secreted markers
8. GWAS-ATAC THRB:10 / HNF4A:12 / RORA:11 motif-disruption counts (clarification: these are variants disrupting these TFs' **binding motifs** across the genome — not variants *at* the THRB/HNF4A/RORA gene loci. Narrative: "these TFs are upstream of GWAS variant action," not "these loci are dense with motif disruption.")

If asked "what's your top finding?" → say **"multi-modal convergence atlas at 28 GWAS × 9 cohorts × cross-species"** — Resource framing, differentiated, not extractable.

---

## H. THE THREE ASKS — exact phrasing for the meeting

**1. (HIGHEST PRIORITY — timeline + journal + collab framing)**
> *"Congrats on the preprint — really impressive multimodal stack. What's the journal target and where are you in revision? I ask because we have a Resource manuscript with a 9-cohort bulk + 28-GWAS COLOC stack that's complementary to yours, and if our windows are close I'd love to coordinate so we end up companion-positioned rather than competing — especially since Bing's now at NYGC and Neville and I would rather build than collide."*

**2. (MID — concrete bidirectional data exchange)**
> *"Two things that would be mutually high-leverage with low overhead: (a) at the canonical MASLD loci — PNPLA3, TM6SF2, HSD17B13, MBOAT7, GCKR — we have multi-ancestry SuSiE-COLOC across 28 GWAS portfolios. Could we trade per-locus posteriors against your caQTL/hQTL effect sizes for a joint locus table? (b) We have 9 mouse RNA-seq datasets across 5 diet models — happy to project your Fib-a/Fib-adj signature against them and send back cross-species concordance scores, on your timeline, attribution your call."*

**3. (LOWEST RISK — relationship building)**
> *"Separately — would you or Bing want to give a Sanjana-lab seminar at NYGC this summer or fall? Neville's series gets the whole genomics floor; nice way to put faces to the work and meet the team."*

**Bonus probe (one to drop casually)**:
> *"Which FNIH DAC oversees the FNIH.Liver consortium data? We're at NYGC; curious whether there's a replication-tier registration path for groups with orthogonal cohorts."*


---

# CORRECTIONS — applied 2026-05-15 from figure-build agents

## C1. THRB / HNF4A / RORA "motif disruption counts" semantics

Earlier briefing repeatedly cited *"GWAS-ATAC motifbreakR disrupts HNF4A 12×, RORA 11×, THRB 8×"* as if these were variant counts at the gene loci. **They are not.** These are counts of variants across the genome whose ALT allele disrupts the **binding motif** of these TFs (TF acting as a trans-regulator).

**Correct framing**: "HNF4A/RORA/THRB are the most-frequently-disrupted TF binding motifs across MASLD credible-set variants — i.e. these TFs are upstream of GWAS variant action genome-wide."

**Not**: "the THRB locus has 8 motif-disrupting variants in its credible set."

This is still a strong differentiator vs Elison/Ren — we have TF-binding-motif-level mechanism inference at scale, which they cover via ChromBPNet/motif enrichment in their hQTL/caQTL pipeline but not via direct cross-GWAS motif-disruption counting.

## C2. THRB locus genetic anchor is GGT, not MASH

Per the figure-build (`figS_thrb_locus.pdf`):
- THRB COLOC PP.H4 = **1.000 only for UKBB-GGT** across 16 GWAS tested
- All NAFLD/NASH/HCC/cirrhosis: PP.H4 < 0.04
- ALT/AST: PP.H4 < 0.03
- Top SNP (chr3:24,520,283, ~25 kb 3' of THRB) lies in a strong ENCODE-liver H3K27ac peak (~65× FC over input), low H3K27me3, hepatocyte-accessible
- SuSiE did not converge at THRB → top SNP is **ABF-fallback**, not a SuSiE credible-set variant

**Implication for meeting / paper**:
- The line "THRB PP.H4=1.000 validates resmetirom as a drug target" is partially defensible but **must be qualified**: the genetic anchor is bile-flow biology (GGT), not direct fibrosis/MASH severity. Resmetirom's MASH efficacy is mechanism-driven (thyromimetic THR-β agonism, hepatocyte cholesterol/lipid handling), not common-variant-validated at the MASH endpoint.
- Soft-collab framing: still strong — we have multi-GWAS triangulation at the locus across 16 phenotypes, which Elison/Ren do not (single Vujkovic cALT GWAS).
- Reviewer defense: explicitly state in the resmetirom-defense paragraph that the genetic anchor is GGT and that the MASH translational case rests on mechanism-of-action + Madrigal MAESTRO-NASH trial data, not direct MASH GWAS colocalization.

## C3. Per-transition decoupleR activity is missing for 6 of 24 disease-regulon TFs

Per the figure-build (`figS05_gwas_atac_tf_transition_heatmap.pdf`):
- **THRB, RORA, KLF15, MLXIPL, NR1H4, RORC have no rows in `RNA-seq/results/stratified_causal/transition_tf_activity.csv`** — i.e. these TFs were not estimated per-transition by Script 211/213 decoupleR runs.
- They appear on the right-side motif-disruption bar but render as NA cells in the heatmap.
- **HNF4A is the visual hero instead** — monotonic activity loss F0→F1 (z=−2.1) → F1→F2 (−4.6) → F2→F3 (−3.4) → F3→F4 (−6.0), all FDR<0.10.

**To fix properly**: re-run upstream Script 211/213 with the 24-TF disease-regulon set as an explicit subset target. Until then, the heatmap leads with HNF4A and the THRB / RORA story is told on the side-annotation only.

## C4. Spatial niche resolution is broad, not Progressor-level

Per the figure-build (`figS_spatial_niche.pdf`):
- The c2l reference uses 16 broad cell types — **Progressor / Moderate / Stable / Neutral / Healthy hepatocytes are not separately deconvoluted at spot level**.
- HSC merged into "Fibroblasts"; LSEC merged into "Endothelial".
- Disease axis is **Healthy → Steatotic only** (GSE192741 has no F3+ donors; HRA007511 embargoed; Vu et al has no control split).
- Headline finding: **Hepatocyte ↔ Fibroblast Δz = +6.53** (flips from −4.41 in Healthy to +2.12 in Steatotic). Macrophage ↔ B cell Δz = +5.93 (TLS-like).

**To upgrade to Progressor-resolved**: rebuild c2l reference from the 657K-hepatocyte 5-meta-subtype atlas (Pipeline 10). Separate pilot, not a same-day swap.
