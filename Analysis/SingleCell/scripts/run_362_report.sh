#!/bin/bash
#SBATCH --job-name=s3_362_report
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --qos=nslab
#SBATCH --output=Analysis/SingleCell/scripts/logs/resolution_sweep/362_report_%j.out
#SBATCH --error=Analysis/SingleCell/scripts/logs/resolution_sweep/362_report_%j.err

set -eo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
micromamba run -n rnaseq Rscript Analysis/SingleCell/scripts/362_aggregate_resolution_report.R
