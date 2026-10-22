#!/bin/bash
#SBATCH --job-name=H35_analysis
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/H35-20261009T180343Z/logs/analysis-%j.out
#SBATCH --error=Analysis/Hypotheses/H35-20261009T180343Z/logs/analysis-%j.err
# H35 step 2: editing index versus fibrosis stage (PRESPEC_H35.md).
# Usage: sbatch run_02_analysis.sbatch PERSON_TABLE OUT_DIR N_PERM
set -euo pipefail
PT="${1:?person table}"; OUT="${2:?output dir (must not exist)}"; NP="${3:?number of permutations}"
RS=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
"$RS" scripts/analysis/hypotheses_2026q4/H35_editing/02_analysis.R "$PT" "$OUT" "$NP"
sha256sum scripts/analysis/hypotheses_2026q4/H35_editing/02_analysis.R scripts/analysis/hypotheses_2026q4/H35_editing/lib_bb.R > "$OUT/script_sha256.txt"
