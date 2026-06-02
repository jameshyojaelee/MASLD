#!/bin/bash
#SBATCH --job-name=scdrs
#SBATCH --partition=bigmem
#SBATCH --cpus-per-task=8
#SBATCH --mem=200G
#SBATCH --time=48:00:00
#SBATCH --output=RNA-seq/logs/scdrs_%j.out
#SBATCH --error=RNA-seq/logs/scdrs_%j.err

# Phase 4A: scDRS disease scoring using TWAS z-scores + scRNA atlas
# Submit: sbatch RNA-seq/run_scdrs.sbatch

set -euo pipefail

BASE="/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design"
cd "${BASE}"

eval "$(micromamba shell hook -s bash)"
micromamba activate causal_inference

export PYTHONNOUSERSITE=1

echo "=== scDRS sbatch ==="
echo "Job ID: ${SLURM_JOB_ID}"
echo "Node: $(hostname)"
echo "CPUs: ${SLURM_CPUS_PER_TASK}"
echo "Memory: ${SLURM_MEM_PER_NODE}MB"
echo "Start: $(date)"
echo ""

# Run scDRS for all GWAS (downsample to 200K cells to stay within memory)
python RNA-seq/69_scdrs_scoring.py --downsample 200000

echo ""
echo "=== Done: $(date) ==="
