#!/bin/bash
#SBATCH --job-name=agp_c1-unmatched
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1_unmatched_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1_unmatched_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUTDIR="$1"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1
export PYTHONPATH="$ROOT/scripts/analysis/alphagenome_program:${PYTHONPATH:-}"
export OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8

echo "host $(hostname)  job ${SLURM_JOB_ID}  start $(date -u)"
"$PY" "$ROOT/scripts/analysis/alphagenome_program/c1_unmatched_coverage.py" --outdir "$OUTDIR"
cp -p "$ROOT/scripts/analysis/alphagenome_program/c1_unmatched_coverage.py" "$OUTDIR/"
echo "end $(date -u)"
