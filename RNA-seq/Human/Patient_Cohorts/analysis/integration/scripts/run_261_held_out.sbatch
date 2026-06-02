#!/bin/bash
#SBATCH --job-name=261_held_out_repl
#SBATCH --partition=cpu
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/261_held_out_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration/logs/261_held_out_%j.err

set -eo pipefail
export PYTHONUNBUFFERED=1

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Human/Patient_Cohorts/analysis/integration

mkdir -p logs

echo "=== Script 261: Held-out replication of cascade hierarchy (Team M7) ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node:   $(hostname)"
echo "Start:  $(date)"
echo ""

eval "$(micromamba shell hook --shell bash)"
micromamba activate rnaseq

set -u

Rscript scripts/261_held_out_replication.R

echo ""
echo "End: $(date)"
