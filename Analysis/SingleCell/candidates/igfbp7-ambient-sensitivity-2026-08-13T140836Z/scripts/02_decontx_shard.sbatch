#!/bin/bash
# Per-dataset decontX shard. DECONTX_ONLY_DATASET selects the dataset; each
# shard writes its own results files, merged afterwards by 02c_merge.py.
#SBATCH --job-name=ambient
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/logs/02_decontx_%x_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/logs/02_decontx_%x_%j.log

set -euo pipefail
umask 027

: "${DECONTX_ONLY_DATASET:?Set DECONTX_ONLY_DATASET}"
export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export CAND_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}

echo "[sbatch] host=$(hostname) job=${SLURM_JOB_ID} dataset=${DECONTX_ONLY_DATASET} start=$(date -u)"
micromamba run -n decontx Rscript "${CAND_ROOT}/scripts/02_decontx.R"
echo "[sbatch] DONE $(date -u)"
