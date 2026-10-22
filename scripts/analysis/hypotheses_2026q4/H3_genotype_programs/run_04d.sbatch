#!/bin/bash
#SBATCH --job-name=H3_sig_hsd_depth
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/H3-20261009T180331Z/logs/%x-%j.out
set -euo pipefail
ROOT=Analysis/Hypotheses/H3-20261009T180331Z
CALLS=$(cat $ROOT/results/LATEST_CALLS.txt)
SCORES=$(cat $ROOT/results/LATEST_SCORES.txt)
GENO=$(cat $ROOT/results/LATEST_GENO.txt)
BASE=$(cat $ROOT/results/LATEST_SIGNATURE_hsd.txt)
OUT=$ROOT/results/signature-hsd-depthadj-$(date -u +%Y%m%dT%H%M%SZ)
SCRIPT=scripts/analysis/hypotheses_2026q4/H3_genotype_programs/04d_hsd_depth_sensitivity.R
export OMP_NUM_THREADS=1
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript $SCRIPT "$CALLS" "$SCORES" "$GENO" "$BASE" "$OUT"
sha256sum $SCRIPT > "$OUT/script_sha256.txt"
echo "$OUT" > $ROOT/results/LATEST_SIGNATURE_hsd_depthadj.txt
