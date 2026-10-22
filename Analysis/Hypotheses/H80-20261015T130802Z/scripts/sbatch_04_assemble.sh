#!/bin/bash
#SBATCH --job-name=6204
#SBATCH --partition=cpu
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H80-20261015T130802Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H80-20261015T130802Z/logs/%j.err
set -euo pipefail
umask 027
H80=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H80-20261015T130802Z
OUT=${H80}/run_assemble            # one leaf per job; must not exist
test ! -e "${OUT}"; mkdir "${OUT}"; mkdir "${OUT}/tmp"
export TMPDIR=${OUT}/tmp XDG_CACHE_HOME=${OUT}/tmp
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
FM=${MASLD_PROJECT_ROOT}/GWAS/finemapping
export AFR_LD_DIR=${FM}/data/ld_ref/1kg_afr_gram AMR_LD_DIR=${FM}/data/ld_ref/1kg_amr_gram
export EAS_LD_DIR=${FM}/data/ld_ref/1kg_eas_gram SAS_LD_DIR=${FM}/data/ld_ref/1kg_sas_gram
export EQTL_SUSIE_DIR=${FM}/results/eqtl_susie_polyfun
unset LD_PANEL
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
echo "host $(hostname) job ${SLURM_JOB_ID} start $(date -Is)"
${ENV}/bin/Rscript ${H80}/scripts/04_assemble.R --setup ${H80}/run_setup --arm_a ${H80}/run_armA --arm_b ${H80}/run_armB_pc1 --readout ${H80}/run_readout --out "${OUT}" --seed 20260923
${ENV}/bin/Rscript -e 'sessionInfo()' > "${OUT}/sessionInfo.txt"
sha256sum ${H80}/scripts/* > "${OUT}/script_sha256.tsv"
echo "end $(date -Is)"
