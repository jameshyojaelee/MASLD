#!/bin/bash
#SBATCH --job-name=ambient
#SBATCH --partition=bigmem
#SBATCH --qos=interactive
#SBATCH --mem=900G
#SBATCH --cpus-per-task=16
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/logs/02_decontx_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/logs/02_decontx_%j.log

set -euo pipefail
umask 027

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export CAND_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "[sbatch] host=$(hostname) job=${SLURM_JOB_ID} start=$(date -u)"
micromamba run -n decontx Rscript "${CAND_ROOT}/scripts/02_decontx.R"
echo "[sbatch] DONE $(date -u)"
