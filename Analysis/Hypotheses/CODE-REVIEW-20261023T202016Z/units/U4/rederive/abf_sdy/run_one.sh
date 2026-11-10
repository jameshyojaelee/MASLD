#!/bin/bash
# One (study, chr) ABF-only pass of the sdY-patched copy of 06_susie_coloc.R.
# Usage: run_one.sh <study> <chr> <sdY_gwas_alt> <sdY_eqtl_alt> <gene filter tsv>
set -eo pipefail
A=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/abf_sdy
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
LDREF=${FM}/data/ld_ref
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AFR_LD_DIR="${LDREF}/1kg_afr_gram" AMR_LD_DIR="${LDREF}/1kg_amr_gram"
export EAS_LD_DIR="${LDREF}/1kg_eas_gram" SAS_LD_DIR="${LDREF}/1kg_sas_gram"
export EQTL_SUSIE_DIR="${FM}/results/eqtl_susie_polyfun"
export SKIP_SUSIE=1 COLOC_OUT_ROOT=${A}/out SDY_GWAS_ALT=$3 SDY_EQTL_ALT=$4 COLOC_GENE_FILTER=$5
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
[[ -f ${A}/out/$1/susie_coloc_chr$2.csv ]] && { echo "skip $1 chr$2"; exit 0; }
mkdir -p ${A}/out/$1/logs
cd ${FM}/src
GWAS_NAME="$1" CHR_FILTER="$2" /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla ${A}/06_susie_coloc_sdy.R \
  > ${A}/out/$1/logs/chr$2.log 2>&1 || echo "FAILED $1 chr$2"
echo "[$(date -u +%FT%TZ)] done $1 chr$2"
