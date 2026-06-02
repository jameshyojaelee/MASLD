#!/bin/bash
#SBATCH --job-name=liver_multi_evidence_v2
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=logs/multi_evidence_v2_%j.log

set -euo pipefail

echo "=== Multi-Evidence Score Redesign v2 ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

# Activate environment
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq

Rscript 27_multi_evidence_score.R

echo "End: $(date)"
