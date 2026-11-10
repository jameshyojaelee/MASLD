#!/bin/bash
#SBATCH --job-name=rv5
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV5-20261102T144853Z/logs/%x_%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV5-20261102T144853Z
OUT=$RUN/out-$(date -u +%Y%m%dT%H%M%SZ)
mkdir "$OUT"
cd "$RUN" && sha256sum -c PRESPEC_RV5.md.sha256
sha256sum scripts/rv5_heterogeneity.R > "$OUT/script_sha256.txt"
export OMP_NUM_THREADS=2
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/rv5_heterogeneity.R "$OUT"
