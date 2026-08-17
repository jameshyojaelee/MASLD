#!/bin/bash
# Re-score with the DANB depth recomputed from the decontaminated matrix.
#
# Runs 03_rescore.py ONLY. 04_diagnostics.py is deliberately not chained: it reads
# the unsuffixed 03_program_stage_effects.tsv and writes 04_*.tsv, so chaining it
# would overwrite the raw_depth results this run exists to be compared against.
# The ambient-burden numbers 04 produces are depth-independent anyway.
#
# 03_pip_freeze.txt is likewise not rewritten; the environment is unchanged from
# the raw_depth run and that file records it.
#SBATCH --job-name=ambient
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/logs/03b_rescore_dec_depth_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z/logs/03b_rescore_dec_depth_%j.log

set -euo pipefail
umask 027

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export CAND_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/Analysis/SingleCell/candidates/igfbp7-ambient-sensitivity-2026-08-13T140836Z
export RESCORE_DEPTH_MODE=dec_depth
export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export MKL_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export PYTHONNOUSERSITE=1
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python

echo "[sbatch] host=$(hostname) job=${SLURM_JOB_ID} mode=${RESCORE_DEPTH_MODE} start=$(date -u)"
"${PY}" "${CAND_ROOT}/scripts/03_rescore.py"
echo "[sbatch] DONE $(date -u)"
