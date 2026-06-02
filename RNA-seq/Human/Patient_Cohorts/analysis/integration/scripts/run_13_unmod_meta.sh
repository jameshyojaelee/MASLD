#!/bin/bash
#SBATCH --job-name=dream-nafl
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/dream_nafl_%j.out
#SBATCH --error=logs/dream_nafl_%j.err

set -eo pipefail
eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts

echo "=== Script 13: NAFL vs NASH DE + meta with unmoderated SE ==="
echo "Started: $(date)"
echo "Node: $(hostname)"

Rscript 13_nafl_vs_nash_de.R

echo "=== Script 13 complete: $(date) ==="
