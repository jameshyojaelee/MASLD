#!/bin/bash
#SBATCH --job-name=dl_allofus
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --partition=cpu
#SBATCH --output=logs/dl_allofus_%j.out
#SBATCH --error=logs/dl_allofus_%j.err

# Download summary statistics for:
#   Chen R et al. (2025) "Trans-ancestral rare variant association study with
#   machine learning-based phenotyping for MASLD."
#   Genome Biology 26:50. doi:10.1186/s13059-025-03518-5 (PMID: 40065360)
#
# Data availability statement (verbatim from paper):
#   "Summary statistics and analysis code supporting the conclusions of this
#   article are available under a Creative Commons Attribution 4.0 International
#   license in Zenodo (10.5281/zenodo.14804238)."
#
# Zenodo record:  https://zenodo.org/records/14804239  (DOI resolves here)
# License:        CC-BY-4.0 (open access, no application required)
#
# Cohort composition (N = 736,010 total across 5 cohorts):
#   UK Biobank PDFF:          38,695 participants
#   UK Biobank MASLD ICD:    423,675 participants
#   All of Us:               229,710 participants (8,504 MASLD cases, 3.7%)
#   BioMe Sample 1:           29,545 participants
#   BioMe Sample 2:           14,388 participants
#
# Ancestry breakdown (pooled):
#   EUR: 546,699 (74.3%) | AFR: 73,225 (9.9%) | AMR: 53,520 (7.3%)
#   SAS:  12,907 (1.8%)  | EAS:  9,404 (1.3%) | MID: 2,374 (0.3%)
#
# Genome build:
#   GWAS summary statistics: GRCh37 (hg19)
#   All other files:         GRCh38 (hg38)
#
# Variant scope:
#   Rare and ultra-rare CODING variants (exome-wide; MAC >= 10 in processed files)
#   Single-variant + gene-level (SKAT/SKATO/ACAT burden tests via regenie)
#   Does NOT include common variant GWAS (those are in GWAS Catalog / other sources)
#
# IMPORTANT NOTE FOR COLOC USE:
#   These are RARE VARIANT (exome) summary statistics, NOT common variant GWAS.
#   Standard COLOC (ABF/SuSiE) requires common variant GWAS with LD structure.
#   These files are best suited for:
#     - Rare variant burden test replication of COLOC hits
#     - Cross-ancestry rare variant validation (AFR/AMR-specific variants)
#     - Gene-level overlap with TWAS/eQTL targets
#   For common variant cross-ancestry COLOC, see BBJ and FinnGen files instead.

DEST=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/MR_Data/AllOfUs
ZENODO_RECORD=14804239
BASE_URL="https://zenodo.org/records/${ZENODO_RECORD}/files"

mkdir -p "${DEST}/logs"
cd "${DEST}"

echo "[$(date)] Starting download from Zenodo record ${ZENODO_RECORD}"

# --- Processed summary statistics (271 MB, MAC>=10 filtered, GRCh37 for GWAS portion) ---
# Contains ancestry-pooled and ancestry-stratified results for:
#   true PDFF, true MASLD, predicted PDFF, predicted MASLD
# Format: .parquet (with Ensembl VEP annotations included)
wget -c "${BASE_URL}/Processed%20summary%20statistics.zip?download=1" \
     -O allofus_processed_sumstats.zip \
     --progress=bar:force 2>&1
echo "[$(date)] Processed summary statistics: $?"

# --- Raw summary statistics (1.9 GB, unfiltered, .tsv/.TBL) ---
# Contains per-cohort raw regenie output before MAC filtering
# Useful if you need cohort-specific (All of Us only) results
wget -c "${BASE_URL}/Raw%20summary%20statistics.zip?download=1" \
     -O allofus_raw_sumstats.zip \
     --progress=bar:force 2>&1
echo "[$(date)] Raw summary statistics: $?"

# --- PheWAS results (277 kB) ---
# Phenome-wide association results for top variants/genes
wget -c "${BASE_URL}/Phewas%20results.zip?download=1" \
     -O allofus_phewas_results.zip \
     --progress=bar:force 2>&1
echo "[$(date)] PheWAS results: $?"

# --- Ensembl VEP annotations (159 MB) ---
# Variant effect predictions for all tested variants
wget -c "${BASE_URL}/Ensembl%20VEP%20annotations.parquet?download=1" \
     -O allofus_vep_annotations.parquet \
     --progress=bar:force 2>&1
echo "[$(date)] VEP annotations: $?"

# --- Analysis scripts ---
wget -c "${BASE_URL}/Plink%20and%20Regenie%20scripts.docx?download=1" \
     -O allofus_plink_regenie_scripts.docx \
     --progress=bar:force 2>&1
echo "[$(date)] Analysis scripts: $?"

# --- ML phenotyping notebook ---
wget -c "${BASE_URL}/Predicting%20MASLD%20and%20PDFF%20in%20the%20UK%20Biobank.ipynb?download=1" \
     -O allofus_ml_phenotyping.ipynb \
     --progress=bar:force 2>&1
echo "[$(date)] ML notebook: $?"

# Unzip summary statistics for downstream use
echo "[$(date)] Unzipping processed summary statistics..."
unzip -n allofus_processed_sumstats.zip -d processed_sumstats/
echo "[$(date)] Unzip exit code: $?"

echo "[$(date)] All downloads complete. Files in ${DEST}:"
ls -lh "${DEST}"/*.zip "${DEST}"/*.parquet "${DEST}"/*.ipynb "${DEST}"/*.docx 2>/dev/null
