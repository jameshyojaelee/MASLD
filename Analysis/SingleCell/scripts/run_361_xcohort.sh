#!/bin/bash
#SBATCH --job-name=s3_361_xcohort
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --qos=nslab
#SBATCH --output=Analysis/SingleCell/scripts/logs/resolution_sweep/361_xcohort_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/logs/resolution_sweep/361_xcohort_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
micromamba run -n rapids_singlecell python Analysis/SingleCell/scripts/361_cross_cohort_progressor.py
