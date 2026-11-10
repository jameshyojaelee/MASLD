#!/bin/bash
#SBATCH --job-name=rv3
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=48G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV3-20261102T144338Z/logs/%x_%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV3-20261102T144338Z
OUT=$RUN/out-$(date -u +%Y%m%dT%H%M%SZ)
mkdir "$OUT"
cd "$RUN" && sha256sum -c PRESPEC_RV3.md.sha256
sha256sum scripts/rv3_random_axis.R > "$OUT/script_sha256.txt"
export OMP_NUM_THREADS=2
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/rv3_random_axis.R "$OUT"
