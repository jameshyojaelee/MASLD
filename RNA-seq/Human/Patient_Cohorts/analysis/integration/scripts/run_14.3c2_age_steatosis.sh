#!/bin/bash
#SBATCH --job-name=limma
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --mem=120G
#SBATCH --cpus-per-task=8
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/14.3c2_age_steatosis_%j.log
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/logs/14.3c2_age_steatosis_%j.err

set -eo pipefail
source ~/.bashrc
micromamba activate rnaseq
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design
export SLURM_CPUS_PER_TASK=${SLURM_CPUS_PER_TASK:-8}
Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/14.3c2_age_steatosis_sensitivity_lvqw.R
