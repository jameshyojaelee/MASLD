#!/bin/bash
# Job 2: stage 1 (controls PC1-PC3, NC1-NC3, diagnostics, power gate G1, stored permutations and bootstrap indices; no new-cohort S1)
#SBATCH --job-name=5402
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/logs/%x_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z/logs/%x_%j.err
# env: /gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript; OMP/OPENBLAS threads 4
# dependency: afterok on job 1 (diagnostic 11 and the CLR components need composition)
# Optional: H54_S1_OUT, H54_COMP_DIR, H54_DEBUG_DIVISOR (debug runs only).
set -euo pipefail
umask 027
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H54-20261015T111917Z
OUT="${H54_S1_OUT:-${RUN}/stage1}"
COMP="${H54_COMP_DIR:-${RUN}/composition}"
test ! -e "${OUT}"
export TMPDIR="${RUN}/tmp_${SLURM_JOB_ID}"; mkdir "${TMPDIR}"
export OMP_NUM_THREADS=4 OPENBLAS_NUM_THREADS=4 FLEXIBLAS_NUM_THREADS=4 MKL_NUM_THREADS=4
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript --vanilla "${RUN}/code/stage1.R" "${OUT}" "${COMP}"
echo "H54_STAGE1_COMPLETE ${OUT} ($(date))"
