#!/bin/bash
#SBATCH --job-name=NMF-nostrip
#SBATCH --partition=bigmem
#SBATCH --qos=nslab
#SBATCH --cpus-per-task=16
#SBATCH --mem=500G
#SBATCH --time=48:00:00
#SBATCH --output=logs/95e_nmf_nostrip_%j.out
#SBATCH --error=logs/95e_nmf_nostrip_%j.err

set -eo pipefail

echo "=== NMF no-strip sensitivity (R1 §9.2) ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node:   $(hostname)"
echo "Start:  $(date)"

eval "$(micromamba shell hook --shell=bash)"
micromamba activate rnaseq

set -u

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq
mkdir -p logs

Rscript 95e_nmf_nostrip_sensitivity.R

echo "End: $(date)"
