# Figure 1 Adobe Illustrator Manual Update Checklist (`2026-11-18`)

> [!IMPORTANT]
> Per User Decision A3, [figure1_JL.pdf](file:///gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/manuscript/latex/figures/main/figure1/figure1_JL.pdf) was **NOT modified programmatically**. Apply the text and numeric updates below manually in Adobe Illustrator before final submission so Figure 1 matches the synchronized manuscript, Figure 2–6 panels, and [RESULTS.md](file:///gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/docs/RESULTS.md).

## Exact Printed Text Audit (`pdftotext` of current `figure1_JL.pdf`) vs Synchronized Resource Values

| Panel / Location | Currently Printed on `figure1_JL.pdf` | Replace With (Synchronized `2026-11-18` Value) |
|---|---|---|
| Embedded top title (if retained) | `Figure 1: Assembly of a harmonized multi-modal MASLD atlas reveals reproducible disease-associated gene programs` | Remove embedded title (LaTeX provides the caption) or align to: `A multimodal human MASLD Resource records liver-trait genetic evidence and disease-state liver expression side by side` |
| Species legend | `Species: Human / Mouse` | **Remove `Mouse`** — the Resource is **human-only** (`Human` liver cohorts only). |
| Ancestry legend | `Ancestry: EUR / AFR / EAS / SAS` | Add **`AMR`**: `EUR, AFR, EAS, AMR, SAS` (**5 GWAS ancestries**). |
| Bulk RNA-seq ring / count | `17 RNA-seq · 358 samples` | **`5` pooled discovery cohorts (`774` genetic-identity clusters: `138` control, `636` disease) + `493` held-out SteatoSITE donors** (`9` transcriptomic cohorts / `1,257` QC libraries total). |
| Single-cell RNA-seq ring / count | `7 scRNA-seq · 117 donors · 1.2M` | **`7` sc/snRNA-seq datasets · `104` source donors · `650,083` unique physical cells** (`145` captures; `226,457` unsorted unique cells across `63` donors; `117` frozen Hotspot gene programs). |
| Genetics ring / count | `23 GWAS · 3 eQTL · 779K subjects` | **`35` Tier-1/2 liver-trait GWAS strata** (`5` ancestries, `6` traits: `15` direct MASLD/NASH/PDFF + `20` ALT/AST/GGT) + **`1` European liver-eQTL panel** (`Broadaway et al., N = 1,183`) + **`1` bulk liver caQTL map** (`Currin et al., N = 138`). |
| Epigenetics / ATAC-seq ring / count | `2 ATAC-seq · 18 samples` | **`3` liver snATAC / multiome cohorts** (`GSE244832`, `GSE281367`, `GSE296875`) + **bulk liver caQTL** (`N = 138` donors). |
| Proteomics ring / count | `2 Proteomics · 72 samples` | **Liver DIA-MS proteome (`PXD051911`: `58` participants — `46` MASLD, `12` controls; `5,836` uniquely mapped proteins)** + **plasma proteome (`143` participants)** + **UKB-PPP plasma panel (`PROT5`)**. |
| Spatial transcriptomics ring / count | `3 Spatial · 12 sections` | **Visium (`Guilliams`, `Vu`, `Yakubovsky`, `HMSMA`) + Xenium (`18` donors: `9` control, `9` MASH-cirrhosis) + CosMx panel annotation (`GSE312698`)**; note GeoMx `GSE287826` is excluded. |
| Panel B bottom-right output block | `Convergence matrix` / `Validations` | Replace with **`MASLD Gene Catalog (28,233 genes)`** (`401` genetically anchored, `26` both supported, `1,118` established-state, `12,955` unresolved, `13,733` untested; no combined rank/score) and **`Next discriminating experiment`** (remove `"Validations"` to respect the two-paper firewall). |

## Key Headline Counts for Reference (If Displayed Anywhere in Panel A/B Callouts)

- **Primary multi-signal SuSiE colocalization (`EQTL-DEDUP-PSD-20261118`)**: **`430`** genes (`417` liver-enzyme driver, `13` direct MASLD/PDFF/NASH driver); **`560`** single-signal-only; **`990`** union.
- **Genetics vs disease-state overlap (`AN1` Rule a, `Fig S1A`)**: **`N = 3,587`** evaluated genes (`420` colocalized, `209` DEGs, **`26` shared** vs `20.9` expected, `MH OR = 1.32 [0.85–2.05]`, exact `p = 0.222`).
- **Bulk disease-state & stage remodeling (`Fig 3`)**: **`1,854`** canonical bulk calls (`774` clusters); **`7,214`** SteatoSITE-replicated F3–F4 stage genes (`1,276` inside hepatocytes uncorrected `C3`; **`188`** ambient-corrected within hepatocytes `C3b`, `0/9` controls, `perm p = 0.000999`).
- **Liver proteomics (`Fig 5`)**: **`3,054`** prioritized universe genes measured (`rho = 0.46`, `188/211 = 89.1%` same sign); **`21 / 96`** scorable Hotspot programs robust to histologic burden (`13` lower, `8` higher).
- **Unified MASLD Gene Catalog (`Fig 6`)**: **`401`** `genetically_anchored` (`394` regulatory + `7` Tier-1 direct-MASLD/PDFF coding: `ADH1B, APOE, MTARC1, PNPLA2, PNPLA3, SERPINA1, TM6SF2`; `GCKR, SAMM50` already regulatory), **`26`** `both_supported`, **`1,118`** `established_state_associated`, **`12,955`** `unresolved`, **`13,733`** `untested`.
