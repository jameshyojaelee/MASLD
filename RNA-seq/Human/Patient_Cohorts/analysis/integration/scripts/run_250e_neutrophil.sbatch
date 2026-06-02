#!/bin/bash
#SBATCH --job-name=250e_neutrophil_reg
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=48:00:00
#SBATCH --output=logs/250e_neutrophil_%j.out
#SBATCH --error=logs/250e_neutrophil_%j.err

set -e

PROJECT_ROOT="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
SCRIPT_DIR="$PROJECT_ROOT/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts"

cd "$SCRIPT_DIR"
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

export MASLD_PROJECT_ROOT="$PROJECT_ROOT"

echo "== 250e_neutrophil_regularization start: $(date) =="
Rscript 250e_neutrophil_regularization.R
echo "== 250e_neutrophil_regularization done: $(date) =="
