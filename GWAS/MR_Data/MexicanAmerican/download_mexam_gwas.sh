#!/usr/bin/env bash
# =============================================================================
# Mexican American / Hispanic Liver Enzyme GWAS — Data Availability Summary
# =============================================================================
#
# PAPER 1 (requested): Sabotta et al. 2022, Front Genet, PMID 36386790
#   DOI: 10.3389/fgene.2022.995488
#   Title: "Genetic variants associated with circulating liver injury markers
#           in Mexican Americans, a population at risk for NAFLD"
#   Cameron County Hispanic Cohort (CCHC), N=564
#
#   STATUS: NOT a GWAS. This is a targeted REPLICATION STUDY.
#   Design: Tested 339 pre-selected variants (previously found in UK/Japan GWAS)
#           in 564 Mexican Americans. Evaluated 2 traits: AST and ALT.
#   Data availability statement: "The original contributions presented in the
#     study are included in the article/Supplementary Materials, further
#     inquiries can be directed to the corresponding author."
#   Supplementary Table S1 (600KB XLSX): Contains genotyping results for
#     the 339 variants per subject + AST/ALT values.
#   GWAS Catalog: NOT registered.
#   Full genome-wide summary statistics: NOT available (only 339 variants).
#   CONCLUSION: CANNOT BE USED FOR COLOC. Only 339 pre-selected SNPs, not
#   genome-wide. No per-SNP beta/SE/p published for full GWAS.
#
# =============================================================================
#
# PAPER 2 (actual GWAS): Young et al. 2019, Obesity, PMID 31219225
#   DOI: 10.1002/oby.22527, PMC: PMC6656610
#   Title: "Genome-Wide Association Study Identifies Loci for Liver Enzyme
#           Concentrations in Mexican Americans: The GUARDIAN Consortium"
#   GUARDIAN Consortium (7 family + non-family Mexican American cohorts)
#   Traits: AST (N=3,644), ALT (N=3,595), GGT (N=1,577)
#   Genome build: hg19 (GRCh37)
#   Array: Illumina HumanOmniExpress BeadChip (~660K SNPs)
#   Imputation: Not mentioned (appears to be chip-level only, not TOPMed)
#
#   GWAS Catalog accessions:
#     GCST008173 — ALT (N=3,595)
#     GCST008174 — AST (N=3,644)
#     GCST008175 — GGT (N=1,577)
#   Full summary statistics in GWAS Catalog FTP: NOT available
#     (fullPvalueSet=False for all three accessions)
#
#   Data availability statement (2019 paper, pre-FAIR era): No formal
#   data availability statement. Supplementary material is a PDF of cohort
#   descriptions only. No public deposition of summary statistics.
#
#   CONCLUSION: CANNOT BE USED FOR COLOC either. Full genome-wide summary
#   statistics were not deposited publicly. Only top hits in paper tables.
#
# =============================================================================
#
# RECOMMENDATION
# =============================================================================
# Neither paper provides genome-wide summary statistics required for COLOC.
#
# Alternative Hispanic/Latino GWAS with publicly available liver enzyme
# summary statistics to consider:
#
#   1. PAGE Study (Pankow et al. or follow-up):
#      Multi-ethnic including Hispanic/Latino; check GWAS Catalog
#
#   2. HCHS/SOL (Hispanic Community Health Study/Study of Latinos):
#      Large N~12,000 Hispanic/Latino; some metabolic traits available
#      Check: https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000810
#
#   3. Million Veteran Program (MVP):
#      Includes ~100K Hispanic/Latino participants; liver enzymes available
#      Summary stats may be available via dbGaP or VA research portal
#
#   4. UK Biobank Latin American subset:
#      Very small N; not recommended
#
#   5. TOPMed / NHLBI BioData Catalyst:
#      Hispanic enriched panels; requires access approval
#
# For immediate use in COLOC: the BBJ study (East Asian, already in pipeline)
# and FinnGen (European Finnish) are available. A well-powered Hispanic/Latino
# GWAS with public summary stats for liver enzymes does not currently exist
# in the GWAS Catalog as of 2026-04-10.
#
# =============================================================================
# If you still want to attempt to obtain data from Young et al. 2019:
# Contact: Lynne E. Wagenknecht, lwgnkcht@wakehealth.edu
# or: The Below Lab (Jennifer Below), Vanderbilt Genetics Institute
# =============================================================================

echo "No downloadable summary statistics available for Mexican American liver enzyme GWAS."
echo "See comments in this script for details."
echo ""
echo "Sabotta 2022 (PMID 36386790): Targeted replication study, NOT a GWAS."
echo "  - 339 pre-selected variants, N=564, no genome-wide data."
echo ""
echo "Young 2019 GUARDIAN (PMID 31219225): Actual GWAS but no public summary stats."
echo "  - GWAS Catalog GCST008173/74/75: fullPvalueSet=False"
echo "  - No FTP files on GWAS Catalog."
