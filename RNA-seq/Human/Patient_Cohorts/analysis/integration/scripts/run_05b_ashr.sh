#!/bin/bash
#SBATCH --job-name=ashr
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/ashr_%j.out
#SBATCH --error=RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/ashr_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/05b_ashr_shrinkage.R
