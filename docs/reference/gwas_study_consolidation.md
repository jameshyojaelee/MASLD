# GWAS Study Consolidation Report

**Date**: 2026-04-09 | **Result**: 24 EUR GWAS → 17 EUR GWAS (7 archived: Anstee_NAFLD, Pazoki_PDFF, 5× Whitfield; 3 FinnGen R12 RESTORED)

---

## Archival Summary

### Round 1 (2026-04-08): 4 studies archived

| Archived Study | Reason | Duplicate Of |
|---------------|--------|-------------|
| `Anstee_NAFLD` | Same PMID 32298765, reformatted differently | `2020_32298765_NAFLD_EUR` |
| ~~`FinnGen_NAFLD`~~ | ~~Probable FinnGen duplicate~~ | **RESTORED** — provenance verified, independent FinnGen R12 data |
| ~~`FinnGen_NASH`~~ | ~~Probable FinnGen duplicate~~ | **RESTORED** — provenance verified, independent FinnGen R12 data |
| ~~`FinnGen_HCC`~~ | ~~Probable FinnGen duplicate~~ | **RESTORED** — provenance verified, independent FinnGen R12 data |

### Round 2 (2026-04-09): 5 studies archived

| Archived Study | Reason | Evidence |
|---------------|--------|----------|
| `2023_36653562_NAFLD_EUR` | Provenance unverified | PMID 36653562 = Kurki et al. FinnGen resource paper, not a liver GWAS |
| `2023_36653562_NASH_EUR` | Provenance unverified | Same PMID issue; "Whitfield JB" not found on PubMed |
| `2023_36653562_Cirrhosis_EUR` | Provenance unverified + case overlap | Identical 6,541 cases as Ghouse_Cirrhosis |
| `2023_36653562_HCC_EUR` | Provenance unverified + case overlap | Identical 1,367 cases as Ghouse_HCC |
| `2023_36653562_Obesity_EUR` | Provenance unverified | Same PMID issue; no verifiable liver GWAS publication |

**Key evidence for Round 2:**
- Beta correlation between Whitfield and Ghouse: r=0.17 (Cirrhosis), r=0.26 (HCC) — not simple copies, but identical case counts prove shared patient sets
- Whitfield SEs are 3.8x larger than Ghouse — inferior statistical precision
- Ghouse studies are fully verifiable (PMID 38632349, 40823170; GCST IDs registered)

---

## Remaining 17 EUR GWAS

| Study ID | Display Label | Phenotype | Type | N | Source |
|----------|--------------|-----------|------|---|--------|
| `UKBB_ALT` | UKBB ALT | ALT | quant | 343,850 | UK Biobank |
| `UKBB_AST` | UKBB AST | AST | quant | 343,850 | UK Biobank |
| `UKBB_GGT` | UKBB GGT | GGT | quant | 343,850 | UK Biobank |
| `2019_31311600_NAFLD_EUR` | eMERGE NAFLD | NAFLD | binary | 8,434 | eMERGE network (EHR) |
| `2020_32298765_NAFLD_EUR` | Biopsy NAFLD | NAFLD | binary | 9,491 | Histologically confirmed |
| `2021_34841290_NAFLD_EUR` | EHR NAFLD | NAFLD | binary | 778,614 | 4-cohort EHR meta-analysis |
| `2023_36280732_NAFLD_deCode_EUR` | deCODE NAFLD | NAFLD | binary | 370,000 | deCODE Genetics |
| `2023_36280732_NAFLD_Intermountain_EUR` | Intermountain NAFLD | NAFLD | binary | 111,000 | Intermountain Healthcare |
| `2023_36280732_NAFLD_UKBB_EUR` | UKBB NAFLD | NAFLD | binary | 400,000 | Sveinbjornsson UKBB cohort |
| `2021_34128465_PDFF_EUR` | Abdominal MRI PDFF | PDFF | quant | 36,116 | Abdominal MRI (Liu/Agrawal) |
| `2021_34957434_PDFF_EUR` | ML-derived PDFF | PDFF | quant | 32,858 | ML-derived MRI (Haas) |
| `2022_36402844_PDFF_EUR` | Whole-body MRI PDFF | PDFF | quant | 44,867 | Whole-body MRI (van der Meer) |
| `Ghouse_Cirrhosis` | Ghouse Cirrhosis | Cirrhosis | binary | 431,122 | deCODE+UKBB meta (Nat Genet 2024) |
| `Ghouse_HCC` | Ghouse HCC | HCC | binary | 310,000 | deCODE+UKBB meta (JHEP Rep 2025) |
| `FinnGen_NAFLD` | FinnGen NAFLD | NAFLD | binary | 438,857 | FinnGen R12 (RESTORED) |
| `FinnGen_NASH` | FinnGen NASH | NASH | binary | ~438,000 | FinnGen R12 (RESTORED) |
| `FinnGen_HCC` | FinnGen HCC | HCC | binary | ~438,000 | FinnGen R12 (RESTORED) |

---

## Naming Convention

**Principle**: Biobank/consortium/method first, trait second, author only for disambiguation.

| Category | Convention | Examples |
|----------|-----------|----------|
| Biobank studies | `{Biobank} {Trait}` | UKBB ALT, deCODE NAFLD, BBJ GGT |
| Phenotype method | `{Method} {Trait}` | Abdominal MRI PDFF, ML-derived PDFF |
| Study design | `{Design} {Trait}` | EHR NAFLD, Biopsy NAFLD, eMERGE NAFLD |
| Meta-analyses | `{First author} {Trait}` | Ghouse Cirrhosis, Ghouse HCC |
| Catalog entries | `GWAS Cat {Trait}` | GWAS Cat PDFF |

---

## COLOC Impact

| Metric | 24 GWAS (original) | 20 GWAS (Round 1) | 15 GWAS (Round 2) | **17 GWAS (FINAL)** |
|--------|--------------------|--------------------|-------------------|-----------------|
| PP.H4>0.5 | 599 | 592 | 504 | **535** |
| PP.H4>0.8 | 275 | 270 | 249 | **255** |
| PP.H4>0.9 | 180 | 176 | 171 | **173** |
| Replicated ≥2 (overlap-corrected) | 179 | 129 | 120 | **121** |
| Replicated ≥3 | 51 | 41 | 36 | **39** |
| Gene-GWAS pairs | 436,460 | 360,578 | 265,763 | **322,423** |

Pazoki_PDFF removal had zero impact on PP.H4 thresholds (its only PP.H4>0.5 hit was already covered by `2022_36402844_PDFF_EUR`). Group replication (118/35) accounts for sample overlap: 3 PDFF studies → 1 group (UKBB_PDFF), Namjou+Ghodsian+Sveinbjornsson_UKBB → 1 group (Ghodsian_meta). 6 genes have inflated raw counts; worst: C2orf16 (9 GWAS → 6 groups).

88 genes lost at PP.H4>0.5 between Round 1 and Round 2 — mostly from the Obesity GWAS (N=818K, 61 hits) which had the highest power among the removed studies.

### Round 3 (2026-04-09): 1 study archived

| Archived Study | Reason | Duplicate Of |
|---------------|--------|-------------|
| `Pazoki_PDFF` | Same GCST90267352 | `2022_36402844_PDFF_EUR` |

---

## Files Updated

| Domain | Files |
|--------|-------|
| GWAS registry | `GWAS/finemapping/config/gwas_registry.tsv` |
| COLOC results | `gene_level_coloc.csv`, `susie_coloc_all_gwas.csv` regenerated |
| Atlas (27a) | FinnGen + Whitfield loading blocks archived; columns retained as NA |
| R scripts | 210, 211 (progression COLOC), 51, 52, 75 |
| Figure scripts | figS04_coloc, figS_sensitivity, figS_finemapping_concordance, fig3_compact, fig6_* |
| Documentation | CLAUDE.md, progress.md, paper_outline.md, NUMBERS.md, methods_draft.md, 04_methods.md, threshold report, 6 design specs, memory files |
| GWAS docs | GWAS_DATASETS.md, README.md, RESULTS.md, FIGURE_GUIDE.md |

---

## Sample Overlap Matrix

Studies are not independent — several share underlying cohort samples. This affects cross-study replication claims.

### Exact Duplicates (archived)
All resolved in Rounds 1-3 above.

### Subsumed Studies
| Parent Study | Child Study | Relationship |
|-------------|------------|--------------|
| Ghodsian (34841290) | Namjou (31311600) | Ghodsian includes eMERGE/Namjou as one of 4 arms (1,106 of 8,434 cases) |

### Major UKBB Sample Overlap
| Study A | Study B | Shared Samples | Impact |
|---------|---------|---------------|--------|
| Ghodsian (34841290) | Sveinbjornsson UKBB (36280732) | ~395K UKBB controls, ~2,558 cases | Cross-replication inflated |
| Liu PDFF (34128465) | Haas PDFF (34957434) | ~33K UKBB MRI cohort | Pseudoreplication |
| Liu PDFF (34128465) | van der Meer PDFF (36402844) | ~33K UKBB MRI cohort | Pseudoreplication |
| Haas PDFF (34957434) | van der Meer PDFF (36402844) | ~33K UKBB MRI cohort | Pseudoreplication |

### Partial Overlap
| Study A | Study B | Overlap |
|---------|---------|---------|
| Ghodsian (34841290) | Archived FinnGen R12 (36653562_*) | FinnGen DF4 (651 cases) ⊂ R12 (~3,504 cases) |

### Ghodsian Cohort Composition
The Ghodsian 2021 NAFLD GWAS (PMID 34841290, N=778,614) is a 4-cohort meta-analysis:
- UK Biobank: 2,558 cases + 395,241 controls
- Estonian Biobank: 4,119 cases + 190,120 controls
- eMERGE (= Namjou 2019): 1,106 cases + 8,571 controls
- FinnGen data freeze 4: 651 cases + 176,248 controls

### Mitigation
Script 07 now computes `coloc_n_groups_h4_05` using gwas_group assignments that collapse overlapping studies:
- 3 PDFF studies → 1 group (UKBB_PDFF)
- Namjou + Ghodsian + Sveinbjornsson_UKBB → 1 group (Ghodsian_meta)
- Other studies remain individual groups
