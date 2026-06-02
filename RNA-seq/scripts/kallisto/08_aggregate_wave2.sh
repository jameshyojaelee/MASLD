#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/aggregate_w2_%j.out
#SBATCH --error=logs/aggregate_w2_%j.err

# Wave 2 - Aggregate kallisto results: tximport per cohort + merge all 9 cohorts.
# Submit AFTER all quant array tasks finish:
#   sbatch --dependency=afterok:QUANT_JOB_ID 08_aggregate_wave2.sh

set -eo pipefail
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/kallisto
Rscript 08_aggregate_wave2.R
