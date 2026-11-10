#!/bin/bash
# Run task <k> of tasks.tsv with the unmodified 06_susie_coloc.R (09b environment).
set -eo pipefail
D=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/CODE-REVIEW-20261023T202016Z/units/U4/rederive/neg_recheck_20261024T111335Z
K=$1
IFS=$'\t' read -r GWAS CHR TASK FILTER < <(awk -F'\t' -v k="$K" 'NR>1 && $3==k' ${D}/tasks.tsv)
FM=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/GWAS/finemapping
LDREF=${FM}/data/ld_ref
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export AFR_LD_DIR="${LDREF}/1kg_afr_gram" AMR_LD_DIR="${LDREF}/1kg_amr_gram"
export EAS_LD_DIR="${LDREF}/1kg_eas_gram" SAS_LD_DIR="${LDREF}/1kg_sas_gram"
export EQTL_SUSIE_DIR="${FM}/results/eqtl_susie_polyfun"
export COLOC_OUT_ROOT="${D}/out" COLOC_GENE_FILTER="${FILTER}"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
mkdir -p ${D}/out/${GWAS}/logs
echo "[$(date -u +%FT%TZ)] task ${K}: ${GWAS} chr${CHR} $(($(wc -l < ${FILTER}) - 1)) genes"
cd ${FM}/src
GWAS_NAME="${GWAS}" CHR_FILTER="${CHR}" /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla ${FM}/src/06_susie_coloc.R > ${D}/out/${GWAS}/logs/chr${CHR}.log 2>&1
echo "[$(date -u +%FT%TZ)] task ${K} done"
