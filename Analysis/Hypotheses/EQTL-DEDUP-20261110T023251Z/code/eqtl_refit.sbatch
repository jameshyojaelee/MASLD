#!/bin/bash
#SBATCH --job-name=eqtl_dedup_A
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-20261110T023251Z/logs/eqtl_refit_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-20261110T023251Z/logs/eqtl_refit_%j.err
# eQTL SuSiE refits (89 genes), one-to-one allelic matching; see PRESPEC.md.
set -eo pipefail
eval "$(micromamba shell hook -s bash)"; micromamba activate rnaseq
export RUN_DIR=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/EQTL-DEDUP-20261110T023251Z
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export LD_PANEL=polyfun
export EQTL_SUSIE_DIR=${RUN_DIR}/eqtl_susie
export TARGET_GENES=${RUN_DIR}/inputs/affected_primary_genes.tsv
export IDENTITY_LOSS_DIR=${RUN_DIR}/identity_losses
export OPENBLAS_NUM_THREADS=4 OMP_NUM_THREADS=4 MKL_NUM_THREADS=4
if [[ -e ${EQTL_SUSIE_DIR} && -z ${RESUME:-} ]]; then echo "EQTL_SUSIE_DIR exists; set RESUME=1 to resume" >&2; exit 1; fi
mkdir -p ${RUN_DIR}/env
Rscript -e 'library(coloc); library(susieR); library(data.table); sessionInfo()' > ${RUN_DIR}/env/sessionInfo_eqtl.txt 2>&1
cd ${MASLD_PROJECT_ROOT}/GWAS/finemapping
for CHR in $(awk -F'\t' 'NR>1 && $7=="TRUE"{print $3}' ${TARGET_GENES} | sort -n | uniq); do
  echo "[$(date -u +%FT%TZ)] chr${CHR}"
  CHR=${CHR} Rscript ${RUN_DIR}/code/run_eqtl_refit.R || echo "FAILED chr${CHR} exit $?"
done
echo "[$(date -u +%FT%TZ)] eQTL refits done"
