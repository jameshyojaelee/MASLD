#!/bin/bash
#SBATCH --job-name=H24_coloc_calls
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=24G
#SBATCH --time=90:00:00
# H24 job B2: rerun of steps 06 (coloc layers) and 07 (calls) after job 22115389
# finished step 05 and stopped in 06 on a locus without variant peaks (L4 code bug).
# Usage: sbatch --output=<OUT>/logs/%x_%j.out --error=<OUT>/logs/%x_%j.err run_B2_coloc_calls.sbatch <OUT>
set -euo pipefail
OUT="${1:?out dir required}"
HERE=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas/scripts/analysis/hypotheses_2026q4/H24_nonhepatocyte_targets
ENV=/gpfs/commons/home/jameslee/micromamba/envs/rnaseq
export OMP_NUM_THREADS=2 PATH="$ENV/bin:$PATH"
test -f "$OUT/results/l1_l2_link_tests.tsv"
$ENV/bin/Rscript "$HERE/06_qtl_coloc.R" "$OUT"
$ENV/bin/Rscript "$HERE/07_calls.R" "$OUT"
echo "H24 job B2 done"
