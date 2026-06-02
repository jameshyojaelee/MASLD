#!/bin/bash
#SBATCH --job-name=salmon_main
#SBATCH --output=logs/salmon_main_%j.out
#SBATCH --error=logs/salmon_main_%j.err
#SBATCH --time=48:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2
#SBATCH --partition=cpu

# Salmon Pseudoalignment Pipeline Submitter
# This script runs Snakemake which submits individual Salmon jobs to SLURM

set -euo pipefail

cd /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/other_diet_mouse_RNAseq

# Activate environment with snakemake and salmon
MAMBA_ENV="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/.mamba/masld_rnaseq_dl"

echo "=========================================="
echo "Starting Salmon Pipeline: $(date)"
echo "=========================================="

# Create log directories
mkdir -p logs/slurm/salmon_quant

# Run Snakemake with SLURM cluster profile
micromamba run -p $MAMBA_ENV snakemake \
    --snakefile workflow/Snakefile \
    --configfile workflow/config.yaml \
    --profile workflow/slurm_profile \
    --jobs 50 \
    --latency-wait 120 \
    --keep-going \
    --rerun-incomplete \
    "$@"

echo "=========================================="
echo "Salmon Pipeline Complete: $(date)"
echo "=========================================="
