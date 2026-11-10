#!/bin/bash
#SBATCH --job-name=rv9
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV9-20261102T151012Z/logs/%x_%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV9-20261102T151012Z
OUT=$RUN/out-$(date -u +%Y%m%dT%H%M%SZ)
mkdir "$OUT"
cd "$RUN" && sha256sum -c PRESPEC_RV9.md.sha256
sha256sum scripts/rv9_replication_state.py > "$OUT/script_sha256.txt"
export PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -I scripts/rv9_replication_state.py "$OUT"
