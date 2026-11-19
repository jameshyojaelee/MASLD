#!/bin/bash
#SBATCH --job-name=eqtl_psd_Bp
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=05:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/logs/coloc_refit_p_%A_%a.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118/logs/coloc_refit_p_%A_%a.err
# 1 study per array task for remaining EQTL-DEDUP-PSD-20261118 studies (race-free resume).
set -eo pipefail
IDX=${SLURM_ARRAY_TASK_ID}
eval "$(micromamba shell hook -s bash)"; micromamba activate rnaseq
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export RUN_DIR=${MASLD_PROJECT_ROOT}/Analysis/Hypotheses/EQTL-DEDUP-PSD-20261118
STUDY_LIST=${RUN_DIR}/inputs/remaining_studies_20261118.txt
GWAS=$(sed -n "${IDX}p" "${STUDY_LIST}")
if [[ -z "${GWAS}" ]]; then
  echo "No study at line ${IDX} of ${STUDY_LIST}" >&2; exit 2
fi
FM=${MASLD_PROJECT_ROOT}/GWAS/finemapping; LDREF=${FM}/data/ld_ref
export LD_PANEL=polyfun
export AFR_LD_DIR=${LDREF}/1kg_afr_gram AMR_LD_DIR=${LDREF}/1kg_amr_gram
export EAS_LD_DIR=${LDREF}/1kg_eas_gram SAS_LD_DIR=${LDREF}/1kg_sas_gram
export EQTL_SUSIE_DIR=${RUN_DIR}/eqtl_susie
export COLOC_OUT_ROOT=${RUN_DIR}/coloc
export COLOC_GENE_FILTER=${RUN_DIR}/inputs/gene_study_filter_with_controls.tsv
export IDENTITY_LOSS_DIR=${RUN_DIR}/identity_losses
export LD_DIAG_DIR=${RUN_DIR}/ld_diagnostics
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
mkdir -p ${RUN_DIR}/env ${COLOC_OUT_ROOT} ${IDENTITY_LOSS_DIR} ${LD_DIAG_DIR}
cd ${FM}
echo "[$(date -u +%FT%TZ)] START task=${IDX} GWAS=${GWAS}"
for CHR in $(awk -F'\t' -v s=${GWAS} 'NR>1 && $1==s{print $2}' ${COLOC_GENE_FILTER} | sort -n | uniq); do
  if [[ -f ${COLOC_OUT_ROOT}/${GWAS}/susie_coloc_chr${CHR}.csv && -f ${COLOC_OUT_ROOT}/${GWAS}/susie_coloc_pairs_chr${CHR}.tsv ]]; then
    echo "skip ${GWAS} chr${CHR}"
    continue
  fi
  # Clean any partial files from a previously interrupted chr before recomputing cleanly
  rm -f "${COLOC_OUT_ROOT}/${GWAS}/susie_coloc_chr${CHR}.csv" \
        "${COLOC_OUT_ROOT}/${GWAS}/susie_coloc_pairs_chr${CHR}.tsv" \
        "${LD_DIAG_DIR}/${GWAS}_chr${CHR}_ld_diag.tsv" \
        "${IDENTITY_LOSS_DIR}/${GWAS}_chr${CHR}_identity_losses.tsv"
  echo "[$(date -u +%FT%TZ)] ${GWAS} chr${CHR}"
  GWAS_NAME=${GWAS} CHR_FILTER=${CHR} Rscript ${RUN_DIR}/code/run_coloc_refit.R || { echo "FAILED ${GWAS} chr${CHR} exit $?" >&2; exit 1; }
done
echo "[$(date -u +%FT%TZ)] DONE task=${IDX} GWAS=${GWAS}"
