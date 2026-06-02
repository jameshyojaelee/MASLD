#!/bin/bash
#SBATCH --job-name=GSE305484_DE
#SBATCH --partition=io
#SBATCH --qos=interactive
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets/GSE305484/logs/deseq2_%j.out
#SBATCH --error=/gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets/GSE305484/logs/deseq2_%j.err

echo "=== GSE305484 DESeq2 (WT WD vs Chow) ==="
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $(hostname)"
echo "Start: $(date)"

# Activate environment (allow unbound vars during conda activation)
set +u
eval "$(micromamba shell hook -s bash)"
micromamba activate rnaseq
set -euo pipefail

# Run DESeq2
Rscript /gpfs/commons/groups/sanjana_lab/Cas13/MASLD_library_design/RNA-seq/Mouse/Western_Diet_Datasets/GSE305484/scripts/01_run_deseq2_wt_only.R

echo "End: $(date)"
echo "=== Done ==="
