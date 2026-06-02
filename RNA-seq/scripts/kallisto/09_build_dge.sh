#!/bin/bash
#SBATCH --job-name=kallisto
#SBATCH --partition=cpu
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --output=logs/build_dge_%j.out
#SBATCH --error=logs/build_dge_%j.err

# Build merged_dge.rds + merged_counts_raw.rds from 9-cohort kallisto counts.
# Submit AFTER aggregation completes:
#   sbatch --dependency=afterok:AGG_JOB_ID 09_build_dge.sh

set -eo pipefail
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/scripts/kallisto
Rscript 09_build_kallisto_dge_9cohort.R
