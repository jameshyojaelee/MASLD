#!/bin/bash
#SBATCH --job-name=agp_c1-matched-endogenous
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1_matched_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/logs/agp_c1_matched_%j.err

set -euo pipefail
ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
OUTDIR="$1"
PY=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python
export PYTHONNOUSERSITE=1
export OMP_NUM_THREADS=8
export OPENBLAS_NUM_THREADS=8
export MKL_NUM_THREADS=8

mkdir -p "$OUTDIR/tables" "$OUTDIR/logs"
echo "host $(hostname)  job ${SLURM_JOB_ID}  start $(date -u)"
"$PY" -V
"$PY" -m pip freeze > "$OUTDIR/pip_freeze.txt" 2>/dev/null || true

"$PY" "$ROOT/scripts/analysis/alphagenome_program/c1_matched_endogenous.py" --outdir "$OUTDIR"

cp -p "$ROOT/scripts/analysis/alphagenome_program/c1_matched_endogenous.py" "$OUTDIR/"
sha256sum "$ROOT/scripts/analysis/alphagenome_program/c1_matched_endogenous.py" \
  > "$OUTDIR/analysis_script.sha256"
echo "end $(date -u)"
