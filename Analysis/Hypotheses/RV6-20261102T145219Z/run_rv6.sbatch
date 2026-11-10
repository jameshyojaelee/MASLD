#!/bin/bash
#SBATCH --job-name=rv6
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV6-20261102T145219Z/logs/%x_%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV6-20261102T145219Z
OUT=$RUN/out-$(date -u +%Y%m%dT%H%M%SZ)
mkdir "$OUT"
cd "$RUN" && sha256sum -c PRESPEC_RV6.md.sha256
sha256sum scripts/rv6_allele_direction.py > "$OUT/script_sha256.txt"
export PYTHONNOUSERSITE=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -m pip freeze > "$OUT/pip_freeze.txt" 2>/dev/null || true
export OMP_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -I scripts/rv6_allele_direction.py "$OUT"
