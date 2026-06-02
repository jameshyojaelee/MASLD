#!/usr/bin/env bash
# =============================================================================
# Download: Korean MAFLD GWAS summary statistics
# =============================================================================
#
# PAPER 1 (PRIMARY — target for this download):
#   Title  : Genome-wide association study of metabolic dysfunction-associated
#             fatty liver disease in a Korean population
#   Journal: Scientific Reports (2024), vol.14, article 9753
#   DOI    : 10.1038/s41598-024-60152-0
#   PMID   : 38679617
#   PMCID  : PMC11056367
#   Authors: Lee Y, Cho EJ, Choe EK, Kwak MS, Yang JI, Oh SW, Yim JY, Chung GE
#
#   Design  : Case-control GWAS
#   Cohort  : H-PEACE / Gene-Environment Interaction and Phenotype (GENIE) study
#             at Seoul National University Hospital Healthcare System Gangnam Center
#   Array   : Affymetrix Axiom Customized Biobank Genotyping Array (Korean Chip,
#             Center for Genome Science, Korea National Institute of Health)
#   Build   : GRCh37 (NCBI build 37)
#   Ancestry: Korean (EAS)
#   Samples :
#     Discovery  — 2,282 MAFLD cases + 4,669 controls  (N=6,951)
#     Validation — 639 MAFLD cases + 1,578 controls    (N=2,217)
#     Total      — 2,921 cases + 6,247 controls         (N=9,168)
#   Top loci: rs738409 & rs3810622 (PNPLA3), rs59148799 (GATAD2A)
#
#   GWAS Catalog: NOT DEPOSITED (as of 2026-04-10)
#   Data availability: "available from the corresponding author on reasonable
#     request" (Goh Eun Chung, corresponding author)
#     Cohort access: http://en-healthcare.snuh.org/HPEACEstudy
#
# PAPER 2 (SECONDARY — larger but 2025 publication):
#   Title  : Genetic variants associated with metabolic dysfunction-associated
#             fatty liver diseases in a Korean population
#   Journal: European Journal of Medical Research (2025)
#   DOI    : 10.1186/s40001-025-02576-6
#   PMID   : 40264238
#   PMCID  : PMC12016408
#   Samples: 4,061 cases + 9,396 controls (N=13,457 total)
#   Build  : GRCh37 (hg19)
#   Top loci: PNPLA3, SAMM50, PARVB; 22q13.3, 19p13.11, 2p23.3
#   Data availability: "Complete dataset will not be made publicly available
#     due to ethics committee restrictions" — must request from corresponding author
#
# =============================================================================
# SUMMARY: Neither study has publicly downloadable GWAS summary statistics.
#   Both use institutional Korean cohort data subject to ethics committee
#   restrictions or author-mediated access. No GWAS Catalog accessions exist.
#   Full summary statistics must be requested directly from the authors.
# =============================================================================

set -euo pipefail

OUTDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG="${OUTDIR}/download.log"

echo "[$(date)] Korean MAFLD GWAS download script" | tee "$LOG"
echo "" | tee -a "$LOG"

# =============================================================================
# Attempt 1: Check GWAS Catalog (GCST accession) — none known as of 2026-04-10
# =============================================================================
echo "=== Checking GWAS Catalog for PMID 38679617 ===" | tee -a "$LOG"
CATALOG_RESPONSE=$(curl -sf \
  "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByPublicationIdPubmedId?pubmedId=38679617" \
  2>/dev/null || echo "CURL_FAILED")

if echo "$CATALOG_RESPONSE" | python3 -c "
import json,sys
d=json.load(sys.stdin)
n=d.get('page',{}).get('totalElements',0)
print(f'GWAS Catalog hits: {n}')
sys.exit(0 if n>0 else 1)
" 2>/dev/null; then
    echo "GWAS Catalog accession found — updating script with accession and downloading..." | tee -a "$LOG"
    # If a GCST accession ever appears, add download here:
    # curl -L "https://www.ebi.ac.uk/gwas/api/search/downloads/summary_stats?accessionId=GCSTXXXXXXX" \
    #   -o "${OUTDIR}/Korean_MAFLD_GWAS_sumstats.tsv.gz"
else
    echo "No GWAS Catalog entry found for PMID 38679617 (expected — not yet deposited)." | tee -a "$LOG"
fi

# =============================================================================
# Attempt 2: Check PMID 40264238 (larger 2025 study)
# =============================================================================
echo "" | tee -a "$LOG"
echo "=== Checking GWAS Catalog for PMID 40264238 ===" | tee -a "$LOG"
CATALOG_RESPONSE2=$(curl -sf \
  "https://www.ebi.ac.uk/gwas/rest/api/studies/search/findByPublicationIdPubmedId?pubmedId=40264238" \
  2>/dev/null || echo "CURL_FAILED")

if echo "$CATALOG_RESPONSE2" | python3 -c "
import json,sys
d=json.load(sys.stdin)
n=d.get('page',{}).get('totalElements',0)
print(f'GWAS Catalog hits: {n}')
sys.exit(0 if n>0 else 1)
" 2>/dev/null; then
    echo "GWAS Catalog accession found for PMID 40264238!" | tee -a "$LOG"
else
    echo "No GWAS Catalog entry found for PMID 40264238 (expected)." | tee -a "$LOG"
fi

# =============================================================================
# STATUS SUMMARY
# =============================================================================
cat <<'EOF' | tee -a "$LOG"

=============================================================================
STATUS: Full GWAS summary statistics are NOT publicly available for either
Korean MAFLD GWAS study. Both papers restrict data access to ethics-approved
requests directed to the corresponding authors.

RECOMMENDED ACTIONS:
  1. Email corresponding author of Paper 1 (Sci Reports 2024):
       Goh Eun Chung — Seoul National University Hospital Gangnam Center
       Request: Full GWAS summary statistics (GRCh37, discovery cohort,
       N=6,951; MAFLD vs. control, all tested variants with beta/OR, SE, P)
       Cohort data portal: http://en-healthcare.snuh.org/HPEACEstudy

  2. Email corresponding author of Paper 2 (Eur J Med Res 2025):
       DOI: 10.1186/s40001-025-02576-6
       Note: Ethics committee restriction means data may not be shareable
       even on request.

  3. Alternative — request data through Korea Biobank (KNIH):
       The Korean Chip genotyping data may be accessible via:
       https://koreabiobank.re.kr/eng/

  4. Alternative — check KoGES / Korea BioBank Array (K-BBAC) portal:
       https://www.nih.go.kr/ko/main/contents.do?menuNo=300489

  NOTE for liftover: Both studies used GRCh37 (hg19). If summary stats are
  obtained, liftover to GRCh38 is required for consistency with the project
  pipeline. Use the existing BBJ liftover script as a template:
    GWAS/MR_Data/BBJ/format_bbj_for_coloc.R
=============================================================================
EOF

echo "" | tee -a "$LOG"
echo "[$(date)] Script complete." | tee -a "$LOG"
