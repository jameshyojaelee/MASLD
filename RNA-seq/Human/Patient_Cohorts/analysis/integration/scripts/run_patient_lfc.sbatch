#!/bin/bash
#SBATCH --job-name=patient_lfc
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/patient_lfc_%j.log

eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

Rscript RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts/patient_level_lfc_analysis.R
