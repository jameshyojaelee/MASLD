# External Data Access (Consolidated)

Sources: `docs/african_ancestry_outreach_template.md`, `docs/dbgap_mvp_application_notes.md`, `docs/hmsma_access_request.md`, `docs/vu_et_al_jhep_2025_data.md`

---

## Resource Overview

| Resource | Type | Status | Why unavailable | Timeline | Contact |
|----------|------|--------|---------|----------|---------|
| MVP GWAS (phs001672.v3.p1) | Multi-ancestry GWAS (90K cases) | NOT APPLIED | dbGaP application needed | 2-4 months | Gawrieh et al. |
| HMSMA (Li et al. Nat Gen 2025) | 35 Visium spatial (7 CTRL, 10 MASLD, 18 MASH) | BLOCKED | Embargoed until June 2026 | Unknown | jin.chai@cldcsw.org |
| Vu et al. (JHEP Reports 2025) | 33 FFPE Visium + CODEX (F0-F4) | BLOCKED | UQeSpace 403 for programmatic access | Browser download or email | Elizabeth.Powell@qimrb.edu.au, Quan.Nguyen@qimrb.edu.au |
| African-ancestry RNA-seq | Bulk RNA-seq | DOES NOT EXIST | No public AFR MASLD RNA-seq dataset | Ongoing monitoring | Multiple cohort PIs |
| All of Us | Multi-ancestry biobank | MONITORING | Liver RNA-seq expected 2026-2027 | Ongoing | allofus.nih.gov |
| GTEx AFR eQTLs | AFR-specific liver eQTLs | AVAILABLE | Extraction needed | 1-2 weeks | -- |

---

## MVP Multi-Ancestry GWAS

- **Accession**: phs001672.v3.p1 (parent) / phs002690 (liver sub-study)
- **Publication**: Gawrieh et al. 2023

| Ancestry | Cases | Controls | Total |
|----------|-------|----------|-------|
| European | ~58K | ~85K | ~143K |
| African American | ~18K | ~24K | ~42K |
| Hispanic | ~12K | ~16K | ~28K |
| Asian | ~2K | ~3K | ~5K |

- **Need**: Summary-level GWAS only (beta, SE, p per variant per ancestry stratum)
- **Application**: PI-level applicant, IRB approval/exemption, MVP-specific DUA
- **Analysis plan**: COLOC per ancestry x Broadaway EUR eQTLs; extend Script 48 cross-ancestry
- **Ancestry mismatch**: EUR eQTLs matched for EUR; attenuated for AFR (most diverse LD) and ASN; partial for HIS
- **Future**: GTEx ~103 AFR liver samples could provide AFR-specific eQTLs
- **Alternative**: Contact authors directly for ancestry-stratified summary stats

---

## HMSMA Spatial Data

- **Contact**: jin.chai@cldcsw.org (Jin Chai, corresponding author)
- **Portal**: https://db.genomics.cn/stomics/hmsma/download -- popup: "under control, contact jin.chai"
- **Value**: 35 Visium samples = 7x our current spatial (GSE192741, n=5); full MASLD spectrum
- **Data requested**: 175 Stomics files, 35 H&E images, Supplementary_Table_1.docx, adata_29samples.h5ad
- **Embargo**: Some data until June 30, 2026
- **Fallback**: Request pre-embargo batch (7 CTRL + 10 MASLD); or H&E images alone to unblock SpaceRanger on existing 794GB FASTQs from NGDC

---

## Vu et al. Spatial Data

- **Publication**: JHEP Reports 8(2):101657 (2026). DOI: 10.1016/j.jhepr.2025.101657
- **Technology**: 10x Visium (FFPE CytAssist), 33 biopsies across 8 arrays
- **Staging**: F0 (5), F1 (9), F2 (3), F3 (9, 8 unique patients), F4 (7)
- **Mean spots/sample**: ~401 (low -- FFPE core biopsies)
- **Unique value**: F0-F4 fibrosis granularity + pathologist per-spot annotations (steatosis, ballooning, inflammation)
- **Data access**: UQeSpace DOI 10.48610/e95155f -- 403 for all programmatic access (CloudFront WAF)
- **Code repo**: https://github.com/BiomedicalMachineLearning/Liver (cloned to `data/vu_et_al_2025/Liver/`)
- **scRNA reference**: GSE136103 (healthy + cirrhotic, not from this study)
- **Format**: Standard SpaceRanger `{data.dir}/{ARRAY}/outs/filtered_feature_bc_matrix.h5`
- **Action**: Browser download at UQeSpace; if fails, email authors (Powell, Nguyen, Irvine)
- **Also includes**: CODEX spatial proteomics (not needed for pipeline but potentially valuable)

---

## African-Ancestry Outreach

### Priority Contacts for Ethnicity Data

| Dataset | PI / Institution | Probability |
|---------|-----------------|-------------|
| GSE130970 | Hoang, UCLA | Moderate (diverse catchment) |
| GSE162694 | Bril, U. Florida | Moderate (bariatric, diverse) |
| GSE174478 | Govaere, Newcastle UK | Low (UK demographics) |

> PRJNA512027 (Gerhard, Geisinger; Low — 99.5% White) was removed from this outreach list 2026-05-15 after the cohort was dropped from the paper for an L0/S0 library-prep batch perfectly confounded with diagnosis label.

### Complementary Approaches (no new RNA-seq needed)
- **GTEx AFR eQTLs**: ~103 AFR liver samples (healthy only); AFR-specific eQTLs for COLOC
- **MVP AFR GWAS**: ~18K AFR NAFLD cases; COLOC with EUR eQTLs (conservative, LD mismatch)
- **Proxy**: AFR GWAS loci cross-referenced against EUR expression atlas

### Manuscript Limitation Text
> "Our atlas lacks African-ancestry liver transcriptomic data. African Americans have paradoxically lower MASLD prevalence despite higher metabolic risk (potentially PNPLA3 allele frequencies). Cross-ancestry COLOC with MVP AFR GWAS and European eQTLs (pending dbGaP) will partially address this gap."

### Emerging Resources
- **All of Us**: ~50% non-European by design; liver RNA-seq expected 2026-2027
- **SteatoSITE** (PMID 37735578): 940-patient resource; check ancestry composition
- **UK Biobank**: ~5% non-White; generating bulk RNA-seq from surgical samples
