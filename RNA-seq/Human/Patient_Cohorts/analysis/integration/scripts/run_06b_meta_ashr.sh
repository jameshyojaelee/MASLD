#!/bin/bash
#SBATCH --job-name=ashr
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/meta_ashr_%j.out
#SBATCH --error=RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/logs/meta_ashr_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/06b_meta_ashr_shrinkage.R
