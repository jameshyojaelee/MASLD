#!/bin/bash
#SBATCH --job-name=eqtl_parity_B
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/logs/coloc_refit_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z/logs/coloc_refit_%j.err
# SuSiE-COLOC on refit eQTL fits for one study group (1-3); not an array.
# Usage: sbatch coloc_refit.sbatch <group 1|2|3>
set -eo pipefail
GROUP=$1; [[ ${GROUP} =~ ^[1-3]$ ]] || { echo "group 1-3" >&2; exit 2; }
eval "$(micromamba shell hook -s bash)"; micromamba activate rnaseq
export RUN_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-PARITY-20261110T041858Z
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
FM=${MASLD_PROJECT_ROOT}/GWAS/finemapping; LDREF=${FM}/data/ld_ref
export LD_PANEL=polyfun
export AFR_LD_DIR=${LDREF}/1kg_afr_gram AMR_LD_DIR=${LDREF}/1kg_amr_gram
export EAS_LD_DIR=${LDREF}/1kg_eas_gram SAS_LD_DIR=${LDREF}/1kg_sas_gram
export EQTL_SUSIE_DIR=${RUN_DIR}/eqtl_susie
export COLOC_OUT_ROOT=${RUN_DIR}/coloc
export COLOC_GENE_FILTER=${RUN_DIR}/inputs/gene_study_filter.tsv
export IDENTITY_LOSS_DIR=${RUN_DIR}/identity_losses
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
Rscript -e 'library(coloc); library(susieR); library(data.table); sessionInfo()' > ${RUN_DIR}/env/sessionInfo_coloc_g${GROUP}.txt 2>&1
cd ${FM}
tail -n +2 ${COLOC_GENE_FILTER} | cut -f1 | sort -u | awk -v g=${GROUP} '(NR-1)%3==g-1' | while read -r GWAS; do
  for CHR in $(awk -F'\t' -v s=${GWAS} '$1==s{print $2}' ${COLOC_GENE_FILTER} | sort -n | uniq); do
    if [[ -f ${COLOC_OUT_ROOT}/${GWAS}/susie_coloc_chr${CHR}.csv && -f ${COLOC_OUT_ROOT}/${GWAS}/susie_coloc_pairs_chr${CHR}.tsv ]]; then
      echo "skip ${GWAS} chr${CHR}"; continue; fi
    echo "[$(date -u +%FT%TZ)] ${GWAS} chr${CHR}"
    GWAS_NAME=${GWAS} CHR_FILTER=${CHR} Rscript ${RUN_DIR}/code/run_coloc_refit.R || echo "FAILED ${GWAS} chr${CHR} exit $?"
  done
done
echo "[$(date -u +%FT%TZ)] group ${GROUP} done"
