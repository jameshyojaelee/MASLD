#!/bin/bash
#SBATCH --job-name=H24_inventory_v2
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=90:00:00
# H24 job A2: inventory v2 (ATAC count rule, physical loci, reference families). No outcome.
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err run_A2_inventory.sbatch <OUT>
set -euo pipefail
OUT="${1:?out dir required}"
HERE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H24_nonhepatocyte_targets
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export OMP_NUM_THREADS=2
$ENV/bin/Rscript "$HERE/04_inventory.R" "$OUT"
echo "H24 job A2 done"
