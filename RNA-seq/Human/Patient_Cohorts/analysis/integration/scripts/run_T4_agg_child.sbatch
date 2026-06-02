#!/bin/bash
#SBATCH --job-name=sex_T4_agg
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_bootstrap/logs/T4_agg_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/audit_sensitivity/sex_bootstrap/logs/T4_agg_%j.err

set -eo pipefail

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

PROJ=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
Rscript ${PROJ}/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/14.4e_aggregate_bootstrap_T4.R
