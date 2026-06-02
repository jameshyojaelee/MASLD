#!/bin/bash
#SBATCH --job-name=metafor
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/metafor_%j.out
#SBATCH --error=logs/metafor_%j.err

set -eo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

echo "=== Script 06: Meta-analysis with unmoderated SE ==="
echo "Started: $(date)"
echo "Node: $(hostname)"

Rscript 06_meta_analysis.R

echo "=== Script 06 complete: $(date) ==="
