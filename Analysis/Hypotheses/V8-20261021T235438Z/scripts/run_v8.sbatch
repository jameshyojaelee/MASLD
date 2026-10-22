#!/bin/bash
#SBATCH --job-name=V8gate
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V8-20261021T235438Z
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V8-20261021T235438Z/logs/%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V8-20261021T235438Z/logs/%j.err
# V8: GSE344087 gate and observability row. No program score.
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1 PYTHONNOUSERSITE=1 PYTHONHASHSEED=0 OMP_NUM_THREADS=8
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/V8-20261021T235438Z
PY=/gpfs/commons/home/jameslee/micromamba/envs/spatial/bin/python
sha256sum "$RUN/PRESPEC_V8.md" "$RUN"/scripts/*.py > "$RUN/script_sha256.txt"
(cd "$RUN/downloads_untrusted" && sha256sum -c SHA256SUMS) > "$RUN/logs/input_sha256_check.txt"
"$PY" -I "$RUN/scripts/01_gate.py"
"$PY" -I -m pip freeze > "$RUN/pip_freeze_spatial.txt"
