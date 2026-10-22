#!/bin/bash
#SBATCH --job-name=H27_summary
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=90:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/Analysis/Hypotheses/H27-20261008T230909Z/logs/slurm-%j.out

# Summary stage only. Job 22120968 wrote every results/ table and then stopped
# in the summary code (a data.table subsetting error). This reads those tables
# and writes summary.json and sessionInfo.txt; nothing is recomputed.
set -euo pipefail
REPO=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
OUT=$REPO/Analysis/Hypotheses/H27-20261008T230909Z
SC=$REPO/scripts/analysis/hypotheses_2026q4/H27_kidney_hf_coloc
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
test -f "$OUT/results/permutation_null.tsv"
test ! -e "$OUT/summary.json"
sha256sum "$SC/h27_kidney_hf_coloc.R" "$SC/run_h27_summary.sbatch" \
  | tee "$OUT/logs/script_sha256_${SLURM_JOB_ID:-NA}.txt"
"$RS" --vanilla "$SC/h27_kidney_hf_coloc.R" summary "$REPO" "$OUT"
