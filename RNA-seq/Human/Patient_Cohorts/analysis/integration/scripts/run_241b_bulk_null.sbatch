#!/bin/bash
#SBATCH --job-name=241b_bulk_null_multi
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=scripts/logs/241b_%j.out
#SBATCH --error=scripts/logs/241b_%j.err

set -eo pipefail
source ~/.bashrc
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

Rscript scripts/241b_bulk_null_multi_method.R
