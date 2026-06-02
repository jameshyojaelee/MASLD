#!/bin/bash
#SBATCH --job-name=225h_diams
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=logs/225h_%j.out
#SBATCH --error=logs/225h_%j.err
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts
eval "$(micromamba shell hook --shell bash)" && micromamba activate spatial
python 225h_diams_within_cohort.py
