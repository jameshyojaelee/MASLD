#!/bin/bash
#SBATCH --job-name=agp_a1-mpra-reliability
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_a1_reliability_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_a1_reliability_%j.err

set -euo pipefail
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
OUT="$1"

echo "host $(hostname)  start $(date -u +%Y-%m-%dT%H:%M:%SZ)  jobid ${SLURM_JOB_ID}"
mkdir -p "$OUT"
cp "$PROJ/scripts/analysis/alphagenome_program/a1_mpra_reliability.py" "$OUT/a1_mpra_reliability.py.run"
sha256sum "$PROJ/scripts/analysis/alphagenome_program/a1_mpra_reliability.py" > "$OUT/script.sha256"
"$PY" "$PROJ/scripts/analysis/alphagenome_program/a1_mpra_reliability.py" "$OUT"
echo "done $(date -u +%Y-%m-%dT%H:%M:%SZ)"
