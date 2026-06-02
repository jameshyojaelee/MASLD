#!/bin/bash
#SBATCH --job-name=S4_aggregate
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/S4_aggregate_%j.out
#SBATCH --error=logs/S4_aggregate_%j.err

set -euo pipefail
cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design_wt/rna-validation
mkdir -p logs
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design

Rscript Analysis/SingleCell/scripts/350_aggregate_compositional_report.R
