#!/bin/bash
#SBATCH --job-name=agp_f1-phase
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f1-phase_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f1-phase_%j.err

set -euo pipefail

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=$PROJ/GWAS/finemapping/results/alphagenome_program/f1-phase-20260914T193608Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python

export PYTHONNOUSERSITE=1
export F1_NPROC=8
export F1_TMP=$OUT/_tmp
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd "$PROJ"
echo "host $(hostname)  job ${SLURM_JOB_ID:-none}  start $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "cpus ${SLURM_CPUS_PER_TASK:-?}  mem ${SLURM_MEM_PER_NODE:-?}MB"

"$PY" "$PROJ/scripts/analysis/alphagenome_program/f1_phase.py" "$OUT"

echo "end $(date -u +%Y-%m-%dT%H:%M:%SZ)"
