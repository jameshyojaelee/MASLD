#!/bin/bash
#SBATCH --job-name=LOCO
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=48:00:00
#SBATCH --output=logs/LOCO_pillarB_%j.out
#SBATCH --error=logs/LOCO_pillarB_%j.err

# Pillar B — Cross-cohort transferability (depends on LOO-CV completion).

set -eo pipefail
mkdir -p logs

eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq
set -u

cd '/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/scripts'
echo "[$(date)] Pillar B transferability on $(hostname)"
Rscript pillar_B_transferability.R
echo "[$(date)] Done Pillar B"
