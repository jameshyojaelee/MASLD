#!/bin/bash
#SBATCH --job-name=agp_f1-phase-ci
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f1-phase-ci_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_f1-phase-ci_%j.err

set -euo pipefail
PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUT=$PROJ/GWAS/finemapping/results/alphagenome_program/f1-phase-20260914T193608Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "$PROJ"
echo "host $(hostname) job ${SLURM_JOB_ID:-none} start $(date -u +%Y-%m-%dT%H:%M:%SZ)"
"$PY" "$PROJ/scripts/analysis/alphagenome_program/f1_phase_intervals.py" "$OUT"
echo "end $(date -u +%Y-%m-%dT%H:%M:%SZ)"
