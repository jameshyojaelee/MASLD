#!/bin/bash -l
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --job-name=sex_coloc_agg
#SBATCH --output=RNA-seq/results/stratified_causal/logs/agg_%j.out
#SBATCH --error=RNA-seq/results/stratified_causal/logs/agg_%j.err

# ---------------------------------------------------------------------------
# Team B1 final step: aggregate sex-stratified COLOC outputs (A1) and run
# the pooled-sex × COLOC threshold sweep (A2).
# Depends on the 6 SuSiE-COLOC array jobs finishing for A1; A2 only depends
# on the canonical susie_coloc_all_gwas.csv (already on disk).
# ---------------------------------------------------------------------------

set -eo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
mkdir -p RNA-seq/results/stratified_causal/logs

# rnaseq env's binutils activate hook references unbound ADDR2LINE; relax `-u` around activate.
set +u
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

echo "=== Team B1 aggregation step ==="
echo "Start: $(date)"

echo ""
echo "--- A1: aggregate sex-stratified per-stratum COLOC ---"
Rscript RNA-seq/216a_aggregate_sex_strat_coloc.R || \
  echo "WARN: 216a failed (likely some strata still running)"

echo ""
echo "--- A2: pooled-sex × COLOC threshold sweep ---"
Rscript RNA-seq/216_sex_coloc_threshold_sweep.R

echo ""
echo "End: $(date)"
