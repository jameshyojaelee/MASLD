#!/bin/bash
#SBATCH --job-name=rv2
#SBATCH --account=nslab
#SBATCH --qos=nslab
#SBATCH --partition=cpu
#SBATCH --prefer=v6
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV2-20261102T144015Z/logs/%x_%j.out
set -euo pipefail
RUN=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/RV2-20261102T144015Z
NPERM=${NPERM:-1000}
OUT=$RUN/out-$(date -u +%Y%m%dT%H%M%SZ)-n${NPERM}
mkdir "$OUT"
cd "$RUN" && sha256sum -c PRESPEC_RV2.md.sha256
sha256sum scripts/rv2_replication.R > "$OUT/script_sha256.txt"
export OMP_NUM_THREADS=4
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin
$ENV/Rscript scripts/rv2_replication.R "$OUT" "$NPERM"
