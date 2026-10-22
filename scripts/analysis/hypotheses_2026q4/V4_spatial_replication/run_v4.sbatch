#!/usr/bin/env bash
#SBATCH --job-name=V4spatial
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=90:00:00
# Usage: sbatch --output=<outdir>/logs/%x_%j.out run_v4.sbatch <outdir> [--smoke]
set -euo pipefail
OUTDIR="$1"; shift
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
SCRIPT=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/V4_spatial_replication/v4_spatial_replication.py
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 MKL_NUM_THREADS=8
mkdir -p "$OUTDIR/env"
"$PY" -m pip freeze > "$OUTDIR/env/pip_freeze.txt"
"$PY" --version > "$OUTDIR/env/python_version.txt"
sha256sum "$SCRIPT" > "$OUTDIR/env/script_sha256.txt"
cd "$OUTDIR"
"$PY" -B "$SCRIPT" --outdir "$OUTDIR" "$@"
echo EXIT=0
