#!/bin/bash
# Download Sveinbjornsson et al. 2022 (Nature Genetics) NAFLD GWAS summary statistics
# GWAS Catalog accession: GCST90294579
# NOTE: Confirm exact file URL at https://www.ebi.ac.uk/gwas/studies/GCST90294579
# before running. Update GWAS_URL below with the actual .tsv.gz download link.
#
# Usage: bash download_sveinbjornsson_gwas.sh

set -euo pipefail

OUTDIR="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data"
mkdir -p "$OUTDIR"

# Primary URL (NHGRI-EBI GWAS Catalog FTP)
GWAS_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90294579/GCST90294579-buildGRCh38.tsv.gz"
# Alternative if primary fails — update to confirmed URL:
# GWAS_URL="https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/GCST90294579/harmonised/GCST90294579.h.tsv.gz"

echo "Downloading Sveinbjornsson 2022 NAFLD GWAS..."
wget -c -P "$OUTDIR" "$GWAS_URL"

echo "Download complete. File saved to $OUTDIR"
ls -lh "$OUTDIR"/GCST90294579*
