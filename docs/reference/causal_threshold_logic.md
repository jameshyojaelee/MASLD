# Causal-pipeline threshold logic across ancestries (GWAS → fine-mapping → COLOC)

**Status: current (2026-07-02).** Single source of truth for *which threshold applies where* across the
three layers of the causal arm — GWAS summary-stat significance, fine-mapping, and colocalization —
and how each behaves across ancestries (EUR / AFR / AMR / EAS / SAS). Manuscript-Methods-ready.

**Relationship to other docs**
- `GWAS/MR_Data/MVP/THRESHOLD_FRAMEWORK.md` — full rationale + citations for the **Layer-1** per-ancestry
  GWAS significance scheme (this doc summarizes it and does not duplicate the derivation).
- `docs/review/2026-06-03-gwas-coloc-solidified.md` — COLOC pipeline solidification history.
- Memory: `project-mvp-gwas-thresholds-2026-06-30`, `project-coloc-ancestry-tiering-2026-07-02`.

## Organizing principle

> **We vary thresholds by ancestry at the GWAS-significance layer, but HOLD them constant at the COLOC layer.**

- At the **GWAS layer**, differing power across ancestries is a discovery/error-rate question, so the
  significance bar is ancestry-specific (panel-matched Kanai; a lenient non-EUR suggestive tier).
- At the **COLOC layer**, the field convention (MVP/Vujkovic 2022, our domain analog) is to **hold the same
  PP.H4 bar across ancestries** and capture the ancestry difference as a **provenance *label*** — not a moved
  number. The reason non-EUR COLOC is exploratory is the **European eQTL panel** (no ancestry-matched liver
  eQTL exists), not a weaker GWAS threshold.

---

## Layer 1 — GWAS summary-stat significance (lead-SNP nomination → fine-mapping input)

Per-ancestry. Implemented in `scripts/mvp/02_lead_snps_twotier_all.sh` (array over all registry GWAS).

| Ancestry | Tier-1 "significant" = `min(5e-8, Kanai)` | Tier-2 "suggestive" | Clump LD panel |
|---|---|---|---|
| **EUR** | 5e-8 | — (EUR clumps at Tier-1 only; no suggestive expansion → main Fig 2 locus set preserved) | `data/1kg_eur` / PolyFun |
| **AFR** | **3.24e-8** (Kanai AFR — the only stratum tightened below 5e-8) | 1e-5 | `ld_ref/1kg_afr` |
| **EAS** | 5e-8 | 1e-5 | `ld_ref/1kg_eas` |
| **AMR** | 5e-8 | 1e-5 | `ld_ref/1kg_amr` |
| **SAS** | 5e-8 | 1e-5 | `ld_ref/1kg_csa` |

- **Tier-1 = `min(5e-8, Kanai_ancestry)`**: field-standard 5e-8 everywhere, tightened only where the
  panel-matched Kanai 2016 value is *stricter* (AFR 3.24e-8). Kanai's non-AFR values (EUR 9.26e-8, EAS
  1.61e-7, AMR 1.83e-7, SAS 9.46e-8) are *more lenient* than 5e-8 (small-1000G-sub-panel artifacts) and are
  floored at 5e-8.
- **Tier-2 suggestive = 1e-5 (Lander-Kruglyak), NON-EUR only** — a discovery-expansion for underpowered
  strata. EUR stays strict.
- **Cross-ancestry replication of EUR-discovered loci** = Bonferroni `0.05/k`, one-sided,
  direction-concordant + binomial sign test (Vujkovic 2022 precedent).
- Gotcha (implemented): the non-EUR `ld_ref/1kg_<anc>` beds carry rsIDs but the clump input is
  `chr:pos:a1:a2` → the `.bim` col-2 must be re-keyed to `chr:pos:a1:a2` before `plink --clump`, else non-EUR
  returns 0 loci.

**Scope:** feeds lead-SNP tables + fine-mapping + the Fig-2 genetics panels. Does **not** filter the COLOC
atlas columns.

## Layer 2 — Fine-mapping (SuSiE; `GWAS/finemapping/src/03 → 04`)

- Consumes Layer-1 lead SNPs; inherits Layer-1 significance (no new p-value threshold).
- **LD panel is ancestry-matched on the GWAS side:**

  | Ancestry | LD reference | ~N |
  |---|---|---|
  | EUR | PolyFun (UKBB, production) / `1kg_eur` (sensitivity) | ~337K / 503 |
  | AFR | `1kg_afr` | ~660 |
  | EAS | `1kg_eas` | 504 |
  | AMR | `1kg_amr` | 347 |
  | SAS | `1kg_csa` | ~489 |

- Caveat (documented, drives the COLOC provenance rule): the non-EUR 1000G sub-panels are small
  (≤~660, at/below the ~500 "already-degrading" floor; Zou 2022) and out-of-sample → SuSiE fine-mapping is
  less reliable for non-EUR than for EUR.

## Layer 3 — COLOC (`GWAS/finemapping/src/06` → `gene_level_coloc.csv` → multi-evidence atlas)

**The GWAS-p-value side is THRESHOLD-FREE.** `06` tests *all* Broadaway eGenes and reports a posterior
(PP.H4); Layer-1 significance thresholds do **not** filter it. Ancestry enters two ways:

### (a) LD, per side
- **GWAS side:** ancestry-matched panel (`06:275-282`, `get_ld_per_locus(ancestry = gwas_ancestry)`), with
  `residual_variance = 1` (out-of-sample fix) + LD regularization + dual SuSiE convergence guards.
- **eQTL side:** **always EUR** (Broadaway liver eQTL, hg19) — no ancestry-matched liver eQTL panel exists.
  This EUR-eQTL substrate is the reason non-EUR COLOC is labeled exploratory.

### (b) PP.H4 confidence tier — UNIFORM across all ancestries (field-grounded)
| Tier | PP.H4 | Precedent |
|---|---|---|
| **High** | ≥ **0.8** (+ purity subflag `PP4/(PP3+PP4) ≥ 0.9`) | eQTL Catalogue (PMID 34493866); FinnGen; MVP purity (PMID 35654975) |
| **Suggestive** | **0.5**–0.8 | GTEx v8 enloc RCP > 0.5 (PMID 32913098) — deliberately loose atlas default |
| **Nominal** | **0.1**–0.5 | fastENLOC permissive floor |

The atlas keeps the **inclusive cross-stratum `max(PP.H4)`** as the headline (loose/suggestive by design);
the tiers are *labels*, not cutoffs. Modal field cutoff for a confident call is 0.8 (range 0.75–0.9); the repo
canonical operational line remains PP.H4.susie > 0.5.

### (c) Provenance label — orthogonal to confidence, NOT a threshold change
- `ancestry_matched` — EUR GWAS × EUR eQTL (headline-eligible).
- `cross_ancestry` — non-EUR GWAS × EUR eQTL → **exploratory regardless of PP.H4 magnitude**.
- `cross_ancestry_replicated` — same 0.5 bar met in **≥ 2 distinct ancestries** (matches the prior
  "PP.H4 > 0.5 in EUR AND EAS" rule).
- **Method matters within cross-ancestry:** `coloc.abf` is LD-free → cross-ancestry errs *conservative*
  (a positive is trustworthy corroboration); `coloc.susie` needs an LD matrix → small non-EUR panels risk
  *false-H4 inflation*, so susie-only non-EUR is the more exploratory arm. (Both posteriors are preserved
  side-by-side: `coloc_best_pp4` = ABF, `coloc_best_susie_pp4` = SuSiE.)
- Underpowered strata (e.g. EAS-NAFLD, 505 cases) need no special-case exclusion — they surface as
  `cross_ancestry` + a low confidence tier automatically. A low/absent PP.H4 in an underpowered stratum is
  *uninformative, not refuting* (Rosen 2025; Wallace 2021).

### (d) Effect-scale / ABF inputs — what actually drives PP.H4 (verified 2026-07-02)
`coloc.abf` computes Wakefield's ABF via `approx.bf.estimates(z = β/√V, V = varbeta, sdY)`, with prior SD
`0.15·sdY` (quant) or `0.2` (cc). **When `beta` + `varbeta`(=SE²) + `sdY` are supplied (and no MAF, so there
is no p-value/MAF reconstruction path — the case in `06`), the sample size `N` and the case-fraction `s` do
NOT enter the per-SNP ABF.** Consequences audited on the MVP run:
- The hardcoded `EQTL_N = 1183` (`06:38`) and the per-stratum `s = N_cases/N_tot` are numerically **inert** for
  the ABF — a low-N eQTL variant is *not* miscalibrated by using the panel-max N. (A raised "8× miscalibration"
  concern was checked and is a false alarm.) The one real lever is `sdY` (via the prior) and `type` (cc/quant).
- **GWAS side:** enzymes are `Mean_INT` (inverse-normal) → `sdY = 1` valid; phecodes → `type = "cc"`. ✓
- **eQTL side (`sdY = 1`):** Broadaway is a cross-cohort meta-analysis of **normalized** liver expression
  (chr22 `mean|β| ≈ 0.05`, `max|β| ≈ 1.34`, `mean SE ≈ 0.05`) — effects on a common standardized scale, so
  `sdY = 1` is valid. ✓

### Implemented columns (Script 07 → `gene_level_coloc.csv`; ingested by 27a → atlas)
Headline `coloc_best_pp4` / `coloc_best_susie_pp4` are UNCHANGED (inclusive max). Added:
`coloc_best_{,susie_}ancestry`, `coloc_{abf,susie}_conf_tier`, `coloc_{abf,susie}_ancestry_matched`,
`coloc_{abf,susie}_headline_cross_anc` (requires ≥ 0.5), `coloc_best_{,susie_}pp4_EUR` (+ `_gwas_EUR`,
the EUR-only companion for definitive-claim consumers), `coloc_n_anc_{,susie_}h4_05`,
`coloc_cross_ancestry_replicated`. Ancestry map: `BBJ_*`=EAS, `PanUKBB_AFR_*`=AFR, `PanUKBB_CSA_*`=SAS,
`MVP_*_{AFR,AMR,EAS,SAS}`=non-EUR, else EUR.

---

## One-line summary per layer
- **GWAS significance:** ancestry-*specific* (5e-8; AFR 3.24e-8; + 1e-5 suggestive non-EUR only).
- **Fine-mapping:** inherits Layer 1; ancestry-matched LD (non-EUR panels small → less reliable).
- **COLOC:** ancestry-*uniform* PP.H4 tiers (0.8 / 0.5 / 0.1); atlas inclusive at 0.5; ancestry captured as a
  **provenance label** (matched vs cross-ancestry) + method (ABF vs SuSiE), **not** a moved threshold.

## Citations (all PubMed/DOI-verified)
- Kanai 2016 (PMID 27305981) — per-ancestry genome-wide significance.
- Lander & Kruglyak 1995 (PMID 7581446) — suggestive 1e-5.
- Vujkovic 2022 (PMID 35654975) — MVP multi-ancestry NAFLD; ancestry-specific GWAS × EUR-GTEx COLOC; holds
  PP4/(PP3+PP4) > 0.9 uniformly across ancestries.
- Giambartolomei 2014 (PMID 24830394) — coloc.abf (PP4 > 0.75; no LD matrix).
- Wallace 2020 (PMID 32310995), Wallace 2021 (PMID 34587156) — abf LD-free vs coloc.susie needs LD.
- Zou 2022 (PMID 35853082) — SuSiE-RSS; small/mismatched-panel false credible sets; panel size > regularization.
- Kanai 2022 (PMID 36643910) — meta-analysis fine-mapping miscalibration (SLALOM).
- GTEx Consortium 2020 (PMID 32913098) — enloc RCP > 0.5 "colocalized gene".
- eQTL Catalogue: Kerimov 2021 (PMID 34493866) — PP.H4 ≥ 0.8.
- Rosen 2025 (PMID 40799583) — don't over-interpret absence of colocalization in underpowered studies.
- Weissbrod 2020 (PMID 33199916) — PolyFun (in-sample large-N LD).
