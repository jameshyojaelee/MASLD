#!/bin/bash
#SBATCH --job-name=H3_verify
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=90:00:00
#SBATCH --chdir=/gpfs/commons/groups/sanjana_lab/Cas13/masld_atlas
#SBATCH --output=Analysis/Hypotheses/H3-20261009T180331Z/logs/%x-%j.out
set -euo pipefail
ROOT=Analysis/Hypotheses/H3-20261009T180331Z
CALLS=$(cat $ROOT/results/LATEST_CALLS.txt)
SCORES=$(cat $ROOT/results/LATEST_SCORES.txt)
TESTS=$(cat $ROOT/results/LATEST_TESTS.txt)
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/python -I scripts/analysis/hypotheses_2026q4/H3_genotype_programs/05_verify.py "$CALLS" "$SCORES" "$TESTS"
/gpfs/commons/home/jameslee/micromamba/envs/rnaseq/bin/Rscript scripts/analysis/hypotheses_2026q4/H3_genotype_programs/05_verify.R "$TESTS"
