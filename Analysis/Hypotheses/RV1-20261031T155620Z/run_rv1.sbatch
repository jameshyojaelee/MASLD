#!/bin/bash
#SBATCH --job-name=rv1
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV1-20261031T155620Z/logs/%x_%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV1-20261031T155620Z
OUT=$RUN/out-$(date -u +%Y%m%dT%H%M%SZ)
mkdir "$OUT"
cd "$RUN" && sha256sum -c PRESPEC_RV1.md.sha256
sha256sum scripts/rv1_export.py scripts/rv1_analysis.R > "$OUT/script_sha256.txt"
export PYTHONNOUSERSITE=1 OMP_NUM_THREADS=4
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
$ENV/python -m pip freeze > "$OUT/pip_freeze.txt"
$ENV/python scripts/rv1_export.py "$OUT"
$ENV/Rscript scripts/rv1_analysis.R "$OUT"
