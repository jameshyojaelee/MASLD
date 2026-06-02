#!/bin/bash
#SBATCH --job-name=sex_sub_agg
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/sex_subsample_aggregate_%j.out
#SBATCH --error=logs/sex_subsample_aggregate_%j.err

mkdir -p logs

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

Rscript 14.4b_aggregate_sex_subsampling.R
