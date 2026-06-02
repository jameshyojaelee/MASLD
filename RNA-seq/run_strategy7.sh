#!/bin/bash
#SBATCH --job-name=liver_multi_evidence
#SBATCH --partition=cpu
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/multi_evidence/logs/strategy7_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/results/multi_evidence/logs/strategy7_%j.err
set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
eval "$(/gpfs/commons/home/jameslee/.local/bin/micromamba shell hook --shell bash)"
micromamba activate rnaseq

mkdir -p results/multi_evidence/logs
echo "Job ID: ${SLURM_JOB_ID:-local}"
echo "Start: $(date)"

Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/27_multi_evidence_score.R

echo "Done: $(date)"
