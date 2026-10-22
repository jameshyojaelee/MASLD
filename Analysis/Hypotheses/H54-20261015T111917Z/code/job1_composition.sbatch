#!/bin/bash
# Job 1: composition (InstaPrism on GSE281797 n94 and GSE225740 n93, symbol-collapsed, integer-rounded counts)
#SBATCH --job-name=5401
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/logs/%x_%j.err
# env: RNA-seq/.mamba/music_deconv (R_LIBS_USER empty), as in run_composition_rebuild.sbatch; set.seed(20260923) before each run
# Output dir: $H54_COMP_OUT (default <run>/composition). Refuses if it exists.
set -eo pipefail
umask 027
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z
OUT="${H54_COMP_OUT:-${RUN}/composition}"
test ! -e "${OUT}"
PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
RNASEQ_ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
REF_SCE="${PROJECT_ROOT}/Analysis/Deconvolution/reference/reference_human_scalesc_sce.rds"
export TMPDIR="${RUN}/tmp_${SLURM_JOB_ID}"; mkdir "${TMPDIR}"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8
echo "prespec sha256: $(sha256sum ${RUN}/PRESPEC_H54.md)"
echo "step 1: symbol-collapsed integer counts ($(date))"
"${RNASEQ_ENV}/bin/Rscript" --vanilla "${RUN}/code/job1_prepare_counts.R" "${OUT}"

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate "${PROJECT_ROOT}/RNA-seq/.mamba/music_deconv"
set -u
export R_LIBS_USER=""
R_PREFIX="${PROJECT_ROOT}/RNA-seq/.mamba/music_deconv"
export R_LIBS_SITE="${R_PREFIX}/lib/R/library"
export LD_LIBRARY_PATH="${R_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
R_BIN="${R_PREFIX}/bin/R"
for DATASET in GSE281797 GSE225740; do
  echo "step 2: InstaPrism ${DATASET} ($(date))"
  "${R_BIN}" --no-save --no-restore --slave --file="${RUN}/code/job1_instaprism_wrapper.R" \
    --args "${REF_SCE}" "${OUT}/bulk/${DATASET}/${DATASET}_counts.tsv" "${OUT}/results/${DATASET}" "${DATASET}" celltype
done
"${R_BIN}" --no-save --no-restore --slave -e 'writeLines(capture.output(sessionInfo()), commandArgs(TRUE)[1])' --args "${OUT}/sessionInfo_instaprism.txt"
sha256sum "${RUN}"/code/job1_* > "${OUT}/code_sha256.txt"
echo "H54_COMPOSITION_COMPLETE ${OUT} ($(date))"
